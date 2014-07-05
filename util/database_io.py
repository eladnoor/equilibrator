import util.django_utils
import logging
import json, gzip, os, csv
from gibbs import models, constants, conditions, reaction

DEFAULT_CITATION_DATA_FILENAME = 'data/citation_data.json'
COMPOUND_NAME_FILE = 'data/kegg_compound_names.tsv'
COMPOUND_RENAME_FILE = 'data/kegg_compound_renaming.tsv'
REACTION_FILE = 'data/kegg_reactions.json.gz'
ENZYME_FILE = 'data/kegg_enzymes.json.gz'
GC_NULLSPACE_FILENAME = 'data/kegg_gc_nullspace.json.gz'
CC_FILENAME = 'data/cc_compounds.json.gz'
ALBERTY_FILENAME = 'data/alberty.json.gz'
DEFAULT_ADDITIONAL_DATA_FILENAME = 'data/additional_compound_data.json'
DOWNLOADS_COMPOUND_PREFIX = 'media/downloads/kegg_compounds'
DOWNLOADS_REACTION_PREFIX = 'media/downloads/kegg_reactions'
DOWNLOADS_PSEUDOISOMER_PREFIX = 'media/downloads/kegg_pseudoisomers'
DOWNLOADS_JSON_FNAME = 'media/downloads/kegg_compounds.json.gz'

# Cache compounds so we can look them up faster.
COMPOUNDS_CACHE = {}
CITATIONS_CACHE = {}

def CheckData():
    json.load(open(DEFAULT_CITATION_DATA_FILENAME))
    json.load(gzip.open(GC_NULLSPACE_FILENAME, 'r'))
    json.load(gzip.open(CC_FILENAME, 'r'))
    json.load(gzip.open(REACTION_FILE, 'r'))
    json.load(gzip.open(ENZYME_FILE, 'r'))
    json.load(open(DEFAULT_ADDITIONAL_DATA_FILENAME, 'r'))
    assert(os.path.exists(COMPOUND_NAME_FILE))
    
def LoadCitationData(json_filename=DEFAULT_CITATION_DATA_FILENAME):
    models.ValueSource.objects.all().delete()
    
    parsed_json = json.load(open(json_filename))

    for cd in parsed_json:
        try:
            data = json.dumps(cd)
            name = cd['name']
            source = models.ValueSource(name=name, data=data)
            source.save()
        except Exception, e:
            logging.error('Error parsing reference %s', cd)
            logging.error(e)
            continue

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
        if compound.kegg_id == 'C00001':
            phase = constants.LIQUID_PHASE_NAME
        else:
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

def LoadKeggCompoundNames(kegg_names_filename=COMPOUND_NAME_FILE,
                          kegg_renaming_filename=COMPOUND_RENAME_FILE):

    cid2names = {} # the first name on the list should be the 'common name'
    for row in csv.DictReader(open(kegg_names_filename, 'r'), delimiter='\t'):
        compound_id = row['CID']
        name = row['common name']
        names = row['all names'].split('|')
        if name not in names:
            raise ValueError('The common name \'%s\' is not in the name list for %s'
                             % (name, compound_id))
        if names[0] != name:
            names.remove(name)
            names = [name] + names
        cid2names[compound_id] = names
    
    cid_replace = {}
    for row in csv.DictReader(open(kegg_renaming_filename, 'r'), delimiter='\t'):
        compound_id = row['CID']
        if compound_id not in cid2names:
            raise ValueError('%s appears in the renaming file, but not in the KEGG list'
                             % compound_id)
        
        command = row['command']
        name = row['name']
        if command.lower() == 'remove':
            # remove 'name' from the list of names
            try:
                cid2names[compound_id].remove(name)
            except ValueError:
                logging.warning('The name %s is not one of the options for %s, '
                                'so it cannot be removed' % (name, compound_id))
        elif command.lower() == 'add':
            # put 'name' in the end of the list (or move it there if it is
            # already in the list)
            if name in cid2names[compound_id]:
                cid2names[compound_id].remove(name)
            cid2names[compound_id] = cid2names[compound_id] + [name]
        elif command.lower() == 'delete':
            del cid2names[compound_id]
        elif command.lower() == 'replace':
            del cid2names[compound_id]
            cid_replace[compound_id] = name
        else:
            raise ValueError('Unknown command: %s' % command)

    for compound_id, names in sorted(cid2names.iteritems()):
        names = map(models.CommonName.GetOrCreate, names + [compound_id])
        
        compound = models.Compound(kegg_id=compound_id, name=names[0])
        compound.save()

        # add the compound_id itself also as an optional name
        map(compound.common_names.add, names)

        compound.save()
    
    return cid_replace
    
def LoadFormationEnergies(energy_json_filenane=CC_FILENAME, priority=1):
    parsed_json = json.load(gzip.open(energy_json_filenane, 'r'))

    compound_id = ''
    for cd in parsed_json:
        try:
            compound_id = cd['CID']
            logging.debug('Handling compound %s', compound_id)
            compound = models.Compound.objects.get(kegg_id=compound_id)
            if compound is None:
                continue
            
            compound.formula = cd.get('formula')
            compound.inchi = cd.get('InChI')
            compound.index = cd.get('compound_index')
            compound.group_vector = cd.get('group_vector')

            mass = cd.get('mass')
            if mass is not None:
                compound.mass = float(mass)
        
            num_electrons = cd.get('num_electrons')
            if num_electrons is not None:
                compound.num_electrons = int(num_electrons)

            # Add the thermodynamic data.
            pmap = cd.get('pmap')
            if not pmap:
                error = cd.get('error')
                if error:
                    compound.no_dg_explanation = error
            else:
                AddPmapToCompound(pmap, compound, priority=priority)
            
            compound.save()

        except Exception, e:
            logging.error(e)
            logging.info('Last compound ID was %s' % compound_id)
            continue

def LoadAlbertyEnergies(alberty_json_filenane=ALBERTY_FILENAME, priority=2):
    parsed_json = json.load(gzip.open(alberty_json_filenane, 'r'))

    for cd in parsed_json:
        try:
            compound_id = 'C%05d' % int(cd['cid'])
            logging.debug('Handling compound %s', compound_id)
            compound = models.Compound.objects.get(kegg_id=compound_id)
            if compound is None:
                continue
            
            # Add the thermodynamic data.
            AddPmapToCompound(cd, compound, priority=priority)
            compound.save()

        except Exception, e:
            logging.error(e)
            continue

def DrawThumbnails():
    for c in models.Compound.objects.all():
        c.WriteStructureThumbnail()
        c.save()

def LoadKeggReactions(cid_replace, reactions_json_filename=REACTION_FILE):
    parsed_json = json.load(gzip.open(reactions_json_filename))

    for rd in parsed_json:
        try:
            rid = rd['RID']

            for coeff_cid_pair in rd['reaction']:
                if coeff_cid_pair[1] in cid_replace:
                    logging.info('replacing %s with %s in reaction %s' % 
                                 (coeff_cid_pair[1],
                                  cid_replace[coeff_cid_pair[1]],
                                  rid))
                    coeff_cid_pair[1] = cid_replace[coeff_cid_pair[1]]

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
            names = ed['names'] # add the EC number also as an optional name
            reactions = ed['reaction_ids']
            
            if not ec:
                raise KeyError('Encountered an enzyme without an EC number.')
            if not names:
                raise KeyError('Common names are required for enzymes (EC) %s.' % ec )
            
            reactions = GetReactions(reactions)
            if not reactions:
                logging.info('Ignoring EC %s since we found no reactions.' % ec)
                continue
            
            # Save first so we can do many-to-many mappings.
            enz = models.Enzyme(ec=ec)
            enz.save()
            
            # Add names, reactions, and compound mappings.
            names = map(models.CommonName.GetOrCreate, names + [ec])
            map(enz.common_names.add, names)
            map(enz.reactions.add, reactions)
            enz.save()
            
        except Exception, e:
            logging.warning('Missing data for ec %s', ec)           
            logging.warning(e)
            continue

def GenerateCompoundThumbnails():
    
    for compound in models.Compound.objects.all():
        # thumbnail starts as None when the compound is created
        # so an empty string will mean there was an error and the thumbnail
        # cannot be created
        compound.WriteStructureThumbnail()
        compound.save()

def LoadAdditionalCompoundData(json_filename=DEFAULT_ADDITIONAL_DATA_FILENAME):
    parsed_json = json.load(open(json_filename, 'r'))

    for cd in parsed_json:
        try:
            cid = cd['CID']
            
            compound = models.Compound.objects.get(kegg_id=cid)
            
            note = cd.get('note')
            preferred_name = cd.get('preferred name')
            details_link = cd.get('details_link')
            pmaps = cd.get('pmaps')
            names = cd.get('names')

            if note:
                compound.note = note
            if preferred_name:
                compound.preferred_name = preferred_name
            if details_link:
                compound.details_link = details_link
            if names:
                names = map(models.CommonName.GetOrCreate, names)
                map(compound.common_names.add, names)
            if pmaps:
                # override the pseudoisomer map that appears in the
                # kegg_compound.json file
                compound.species_groups.clear()
                for pmap in pmaps:
                    priority = pmap['priority']
                    AddPmapToCompound(pmap, compound, priority=priority)
                    
            compound.save()
        except Exception, e:
            logging.error('Error parsing cid %s', cid)
            logging.error(e)
            continue
    
def export_database():
    
    export_compounds(priority=2, name='Alberty',
                     ionic_strength=constants.DEFAULT_IONIC_STRENGTH,
                     pMg=constants.DEFAULT_PMG,
                     pH_list=constants.PH_RANGE_VALUES)

    export_reactions(priority=1, name='CC',
                     ionic_strength=constants.DEFAULT_IONIC_STRENGTH,
                     pMg=constants.DEFAULT_PMG,
                     pH_list=constants.PH_RANGE_VALUES)

    export_json()
    
def export_json():
    logging.info("Writing compound data to JSON file: %s" % DOWNLOADS_JSON_FNAME)
    rowdicts = []
    for c in models.Compound.objects.all():
        d = {'name': str(c.name), 
             'KEGG_ID': c.kegg_id,
             'InChI': c.inchi,
             'mass': c.mass,
             'formula': c.formula,
             'num_electrons': c.num_electrons}
        rowdicts.append(d)
    json.dump(rowdicts, gzip.open(DOWNLOADS_JSON_FNAME, 'w'), sort_keys=True, indent=4)

def export_reactions(priority, name, ionic_strength, pMg, pH_list):
    csv_reaction_dict = {}
    for pH in pH_list:
        reaction_fname = DOWNLOADS_REACTION_PREFIX + '_%s_ph%.1f.csv.gz' % (name, pH)
        logging.info("Writing dG'0_r for %s at pH %g to: %s"
                     % (name, pH, reaction_fname))
        csv_reaction_dict[pH] = csv.writer(gzip.open(reaction_fname, 'w'))
        csv_reaction_dict[pH].writerow(["!MiriamID::urn:miriam:kegg.reaction",
                                        "!dG0_prime (kJ/mol)", 
                                        "!sigma[dG0] (kJ/mol)",
                                        "!pH",
                                        "!I (mM)", "!T (Kelvin)", "!Note"])
    
    for r in models.StoredReaction.objects.all():
        rxn = r.ToReaction()
        try:
            dG0_std = rxn.DeltaGUncertainty()
            if dG0_std is not None:
                dG0_std = round(dG0_std, 1)
        except Exception as e:
            logging.debug(str(e))
            dG0_std = None
        
        if dG0_std is None or dG0_std > 200:
            for pH in pH_list:
                row = (r.kegg_id, None, None, pH, ionic_strength,
                       constants.DEFAULT_TEMP, 'uncertainty is too high')
                csv_reaction_dict[pH].writerow(row)
            continue
        
        for pH in pH_list:
            rxn.aq_params = conditions.AqueousParams(pH=pH, pMg=pMg,
                                                     ionic_strength=ionic_strength,
                                                     max_priority=priority)
                                                     
            try:
                dG0_prime = rxn.DeltaG0Prime()
                dG0_prime = round(dG0_prime, 1)
                comment = None
            except Exception:
                logging.warning(str(e))
                dG0_prime = None
                comment = rxn.NoDeltaGExplanation()
            
            row = (r.kegg_id, dG0_prime, dG0_std, pH, ionic_strength,
                   constants.DEFAULT_TEMP, comment)
            csv_reaction_dict[pH].writerow(row)
                
def export_compounds(priority, name, ionic_strength, pMg, pH_list):
    pseudoisomer_fname = DOWNLOADS_PSEUDOISOMER_PREFIX + '_%s.csv.gz' % name
    csv_pseudoisomers = csv.writer(gzip.open(pseudoisomer_fname, 'w'))
    csv_pseudoisomers.writerow(["!MiriamID::urn:miriam:kegg.compound",
                                "!Name", "!dG0 (kJ/mol)",
                                "!nH", "!charge", "!nMg", "!Note"])

    csv_compound_dict = {}
    for pH in pH_list:
        compound_fname = DOWNLOADS_COMPOUND_PREFIX + '_%s_ph%.1f.csv.gz' % (name, pH)
        csv_compound_dict[pH] = csv.writer(gzip.open(compound_fname, 'w'))
        csv_compound_dict[pH].writerow(["!MiriamID::urn:miriam:kegg.compound",
                                        "!Name", "!dG0_prime (kJ/mol)",
                                        "!pH", "!I (mM)", "!T (Kelvin)",
                                        "!Note"])
    
    logging.info("Writing chemical and biochemical formation energies for %s to: %s" %
                 (name, pseudoisomer_fname))
    for compound in models.Compound.objects.all():
        phase = compound.GetDefaultPhaseName()
        rows = compound.ToCSVdG0(priority, phase=phase)        
        csv_pseudoisomers.writerows(rows)
        for pH in pH_list:
            aq_params = conditions.AqueousParams(pH=pH, pMg=pMg,
                                                 ionic_strength=ionic_strength)
            rows = compound.ToCSVdG0Prime(priority, aq_params=aq_params, 
                                          phase=phase)
            csv_compound_dict[pH].writerows(rows)
