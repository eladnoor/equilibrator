import json
import logging

import util.django_utils
import gzip

from gibbs import models
from gibbs import constants

# Cache compounds so we can look them up faster.
COMPOUNDS_CACHE = {}
CITATIONS_CACHE = {}


def GetOrCreateNames(names_list):
    """Find all the names in the database.
    
    Create them if they are not present.
    """
    return [models.CommonName.GetOrCreate(n)
            for n in names_list] 
    
    
def GetReactions(rids_list):
    """Find all the given reactions in the database.
    
    Skip those that are not present.
    """
    if not rids_list:
        return []
    
    rxns = []
    for rid in rids_list:
        try:
            rxns.append(models.StoredReaction.objects.get(kegg_id=rid))
        except Exception:
            logging.warning('Failed to retrieve reaction %s', rid)
            continue
    return rxns
        

def GetCompounds(cids_list):
    """Find all the given compounds in the database.
    
    Skip those that are not present.
    """
    if not cids_list:
        return []
    
    global COMPOUNDS_CACHE
    
    compounds = []
    for kegg_id in cids_list:
        try:
            compounds.append(models.Compound.objects.get(kegg_id=kegg_id))
        except Exception:
            logging.warning('Failed to retrieve compound %s', kegg_id)
            continue
    return compounds
        

COMPOUND_FILE = 'data/kegg_compounds.json.gz'
REACTION_FILE = 'data/kegg_reactions.json.gz'
ENZYME_FILE = 'data/kegg_enzymes.json.gz'
GC_NULLSPACE_FILENAME = 'data/kegg_gc_nullspace.json.gz'

def GetSource(source_string):
    if not source_string:
        return None
    
    lsource = source_string.strip().lower()
    if lsource in CITATIONS_CACHE:
        return CITATIONS_CACHE[lsource]
    
    try:
        source = models.ValueSource.objects.get(name__iexact=lsource)
        CITATIONS_CACHE[lsource] = source
        return source
    except Exception, e:
        logging.warning('Failed to find source "%s"', source_string)
        logging.error(e)
        logging.fatal('Bailing!')
        
    return None

def AddPmapToCompound(pmap, compound):
    source_string = pmap.get('source')
    source = GetSource(source_string)
    if not source:
        logging.error('Failed to get source %s', source_string)
        return 

    logging.debug('Writing data from source %s', source.name)
    
    if 'priority' not in pmap or 'species' not in pmap:
        logging.error('Malformed pmap field for %s', compound.kegg_id)
        return
    
    sg = models.SpeciesGroup(kegg_id=compound.kegg_id,
                             priority=pmap['priority'],
                             formation_energy_source=source)
    sg.save()

    for sdict in pmap['species']:
        formation_energy = sdict.get('dG0_f', None)
        if formation_energy is None:
            logging.error('A specie of %s is missing its formation energy' %
                           compound.kegg_id)

        number_of_hydrogens = sdict.get('nH', None)
        number_of_mgs = sdict.get('nMg', None)
        net_charge = sdict.get('z', None)
        phase = sdict.get('phase', constants.DEFAULT_PHASE)
    

        if phase == constants.AQUEOUS_PHASE_NAME:
            if None in [number_of_hydrogens, number_of_mgs, net_charge]:
                logging.error('An aqueous specie of %s is missing essential info' %
                              compound.kegg_id)
                raise ValueError
        else:
            if number_of_hydrogens is None:
                number_of_hydrogens = 0
            if number_of_mgs is None:
                number_of_mgs = 0
            if net_charge is None:
                net_charge = 0
        
        specie = models.Specie(kegg_id=compound.kegg_id,
                               number_of_hydrogens=number_of_hydrogens,
                               number_of_mgs=number_of_mgs,
                               net_charge=net_charge,
                               formation_energy=formation_energy,
                               phase=phase)
        specie.save()
        sg.species.add(specie)
    
    sg.save()
    compound.species_groups.add(sg)

def LoadKeggCompounds(kegg_json_filename=COMPOUND_FILE, draw_thumbnails=True):
    parsed_json = json.load(gzip.open(kegg_json_filename, 'r'))
    
    for cd in parsed_json:
        try:
            cid = cd['CID']
            logging.debug('Handling compound %s', cid)
            
            formula = cd.get('formula')
            mass = cd.get('mass')
            if mass is not None:
                mass = float(mass)
            inchi = cd.get('InChI')
            num_electrons = cd.get('num_electrons')
            group_vector = cd.get('group_vector')
            name = models.CommonName.GetOrCreate(cd['name'])
        
            """
            if formula is None:
                raise KeyError('Missing formula for CID %s' % cid)
            
            if mass is None:
                raise KeyError('Missing mass for CID %s' % cid)

            if inchi is None:
                raise KeyError('Missing inchi for CID %s' % cid)
            """
             
            c = models.Compound(kegg_id=cid,
                                formula=formula,
                                inchi=inchi,
                                mass=mass,
                                name=name,
                                group_vector=group_vector)
            
            if num_electrons is not None:
                c.num_electrons = int(num_electrons)
            if draw_thumbnails:
                c.WriteStructureThumbnail()
            
            c.save()

            # Add the thermodynamic data.
            pmaps = cd.get('pmaps')
            if not pmaps:
                error = cd.get('error')
                if error:
                    c.no_dg_explanation = error
            else:
                for pmap in pmaps:
                    AddPmapToCompound(pmap, c)
            
            # Add the common names.
            names = GetOrCreateNames(cd['names'])
            for n in names:
                c.common_names.add(n)
            c.save()
        except Exception, e:
            logging.error(e)
            continue
        

def LoadKeggReactions(reactions_json_filename=REACTION_FILE):
    parsed_json = json.load(gzip.open(reactions_json_filename))

    for rd in parsed_json:
        try:
            rid = rd['RID']
            rxn = rd['reaction']
            substrates = []
            products = []
            for coeff, cid in rxn:
                reactant = models.Reactant.GetOrCreate(cid, abs(coeff))
                if coeff < 0:
                    substrates.append(reactant)
                else:
                    products.append(reactant)
                
            # Need to save once.
            rxn = models.StoredReaction(kegg_id=rid)
            rxn.save()
            
            for reactant in substrates:
                rxn.substrates.add(reactant)
            for product in products:
                rxn.products.add(product)
            rxn.hash = rxn.GetHash()
            rxn.save()
        except Exception, e:
            logging.warning('Missing data for rid %s', rid)
            logging.warning(e)
            continue


def LoadKeggEnzymes(enzymes_json_filename=ENZYME_FILE):
    parsed_json = json.load(gzip.open(enzymes_json_filename))

    for ed in parsed_json:
        try:
            ec = ed['EC']
            names = ed['names']
            reactions = ed['reaction_ids']
            
            if not ec:
                raise KeyError('Encountered an enzyme without an EC number.')
            if not names:
                raise KeyError('Common names are required for enzymes (EC) %s.' % ec )
            
            names = GetOrCreateNames(names)
            reactions = GetReactions(reactions)
            if not reactions:
                logging.info('Ignoring EC %s since we found no reactions.' % ec)
                continue
            
            substrates = GetCompounds(ed.get('substrates'))
            products = GetCompounds(ed.get('products'))
            cofactors = GetCompounds(ed.get('cofactors'))
            
            # Save first so we can do many-to-many mappings.
            enz = models.Enzyme(ec=ec)
            enz.save()
            
            # Add names, reactions, and compound mappings.
            map(enz.common_names.add, names)
            map(enz.reactions.add, reactions)
            map(enz.substrates.add, substrates)
            map(enz.products.add, products)
            map(enz.cofactors.add, cofactors)
            enz.save()
            
        except Exception, e:
            logging.warning('Missing data for ec %s', ec)           
            logging.warning(e)
            continue


#def LoadKeggGCNullspace(gc_nullspace_filename=GC_NULLSPACE_FILENAME):
#    parsed_json = json.load(gzip.open(gc_nullspace_filename))
#    
#    for rd in parsed_json:
#        claw = models.ConservationLaw()
#        claw.msg = rd['msg']
#        claw.save()
#        for coeff, cid in rd['reaction']:
#            reactant = models.Reactant.GetOrCreate(cid, coeff)
#            claw.reactants.add(reactant)
#        claw.save()


def LoadKeggGCNullspace(gc_nullspace_filename=GC_NULLSPACE_FILENAME):
    parsed_json = json.load(gzip.open(gc_nullspace_filename))
    
    for rd in parsed_json:
        claw = models.ConservationLaw()
        claw.msg = rd['msg']
        claw.reactants = json.dumps(rd['reaction'])
        claw.save()


def CheckData(filenames=(COMPOUND_FILE,
                         REACTION_FILE,
                         ENZYME_FILE)):
    for json_fname in filenames:
        json.load(gzip.open(json_fname, 'r'))


def LoadAllKeggData(draw_thumbnails=True):
    LoadKeggGCNullspace()
    LoadKeggCompounds(draw_thumbnails=draw_thumbnails)
    LoadKeggReactions()
    LoadKeggEnzymes()


if __name__ == '__main__':
    LoadAllKeggData()
