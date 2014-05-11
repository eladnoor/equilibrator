import json, csv
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
    """
        Find all the given compounds in the database.
    
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
        

COMPOUND_NAME_FILE = 'data/kegg_compound_names.tsv'
REACTION_FILE = 'data/kegg_reactions.json.gz'
ENZYME_FILE = 'data/kegg_enzymes.json.gz'
GC_NULLSPACE_FILENAME = 'data/kegg_gc_nullspace.json.gz'
CC_FILENAME = 'data/cc_compounds.json.gz'

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

def AddPmapToCompound(pmap, compound, priority=1):
    source_string = pmap.get('source')
    source = GetSource(source_string)
    if not source:
        logging.error('Failed to get source %s', source_string)
        return 

    logging.debug('Writing data from source %s', source.name)
    
    if 'species' not in pmap:
        logging.error('Malformed pmap field for %s', compound.kegg_id)
        return
    
    sg = models.SpeciesGroup(kegg_id=compound.kegg_id,
                             priority=priority,
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

def LoadKeggGCNullspace(gc_nullspace_filename=GC_NULLSPACE_FILENAME):
    parsed_json = json.load(gzip.open(gc_nullspace_filename))
    
    for rd in parsed_json:
        claw = models.ConservationLaw()
        claw.msg = rd['msg']
        claw.reactants = json.dumps(rd['reaction'])
        claw.save()

def LoadKeggCompoundNames(kegg_names_filename=COMPOUND_NAME_FILE):
    for row in csv.reader(open(kegg_names_filename, 'r'), delimiter='\t'):
        compound_id = row[0]
        name = models.CommonName.GetOrCreate(row[1])
        
        # split the list of common names (delimiter is |) and also add the 
        # KEGG compound ID as one of them
        common_names = GetOrCreateNames(row[2].split('|') + [compound_id])
        
        c = models.Compound(kegg_id=compound_id, name=name)
        c.save()

        for common_name in common_names:
            c.common_names.add(common_name)

        c.save()
    
def LoadFormationEnergies(energy_json_filenane=CC_FILENAME, priority=1):
    parsed_json = json.load(gzip.open(energy_json_filenane, 'r'))

    for cd in parsed_json:
        try:
            compound_id = cd['CID']
            logging.debug('Handling compound %s', compound_id)
            c = models.Compound.objects.get(kegg_id=compound_id)
            
            c.formula = cd.get('formula')
            c.inchi = cd.get('InChI')
            c.group_vector = cd.get('group_vector')

            mass = cd.get('mass')
            if mass is not None:
                c.mass = float(mass)
        
            num_electrons = cd.get('num_electrons')
            if num_electrons is not None:
                c.num_electrons = int(num_electrons)

            # Add the thermodynamic data.
            pmap = cd.get('pmap')
            if not pmap:
                error = cd.get('error')
                if error:
                    c.no_dg_explanation = error
            else:
                AddPmapToCompound(pmap, c, priority=priority)
            
            c.save()

        except Exception, e:
            logging.error(e)
            continue
            
def DrawThumbnails():
    for c in models.Compound.objects.all():
        c.WriteStructureThumbnail()
        c.save()

def LoadKeggReactions(reactions_json_filename=REACTION_FILE):
    parsed_json = json.load(gzip.open(reactions_json_filename))

    for rd in parsed_json:
        try:
            rid = rd['RID']
            rxn = models.StoredReaction.FromJson(rd)
            rxn.GenerateHash()
            rxn.save()
        except Exception, e:
            logging.warning('Missing data for reaction %s', rid)           
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
            
            # Save first so we can do many-to-many mappings.
            enz = models.Enzyme(ec=ec)
            enz.save()
            
            # Add names, reactions, and compound mappings.
            map(enz.common_names.add, names)
            map(enz.reactions.add, reactions)
            enz.save()
            
        except Exception, e:
            logging.warning('Missing data for ec %s', ec)           
            logging.warning(e)
            continue

def CheckData(filenames=(GC_NULLSPACE_FILENAME,
                         CC_FILENAME,
                         REACTION_FILE,
                         ENZYME_FILE)):
    for json_fname in filenames:
        json.load(gzip.open(json_fname, 'r'))

def LoadAllKeggData():
    #LoadKeggGCNullspace()
    LoadKeggCompoundNames()
    LoadFormationEnergies()
    #DrawThumbnails()
    LoadKeggReactions()
    LoadKeggEnzymes()


if __name__ == '__main__':
    LoadAllKeggData()
