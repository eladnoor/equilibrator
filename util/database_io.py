import logging
import json
import gzip
import os
import csv
from django.apps import apps
from django.db.utils import DataError
from gibbs import conditions
from util import constants
import pandas as pd

DEFAULT_CITATION_DATA_PATH = 'data/citation_data.json'
COMPOUND_NAME_PATH = 'data/kegg_compound_names.tsv'
COMPOUND_RENAME_PATH = 'data/kegg_compound_renaming.tsv'
REACTION_PATH = 'data/kegg_reactions.json.gz'
ENZYME_PATH = 'data/kegg_enzymes.json.gz'
CC_PATH = 'data/cc_compounds.json.gz'
ALBERTY_PATH = 'data/alberty.json'
DEFAULT_ADDITIONAL_DATA_PATH = 'data/additional_compound_data.json'

DOWNLOADS_COMPOUND_TEMPLATE = 'static/downloads/kegg_compounds_%s_ph%.1f.csv.gz'
DOWNLOADS_REACTION_TEMPLATE = 'static/downloads/kegg_reactions_%s_ph%.1f.csv.gz'
DOWNLOADS_PSEUDOISOMER_TEMPLATE = 'static/downloads/kegg_pseudoisomers_%s.csv.gz'
DOWNLOADS_JSON_FNAME = 'static/downloads/kegg_compounds.json.gz'

# Cache compounds so we can look them up faster.
CITATIONS_CACHE = {}

def CheckData():
    json.load(open(DEFAULT_CITATION_DATA_PATH))
    json.load(open(ALBERTY_PATH, 'rt'))
    json.load(open(DEFAULT_ADDITIONAL_DATA_PATH, 'rt'))
    json.load(gzip.open(CC_PATH, 'rt'))
    json.load(gzip.open(REACTION_PATH, 'rt'))
    json.load(gzip.open(ENZYME_PATH, 'rt'))
    assert(os.path.exists(COMPOUND_NAME_PATH))
    assert(os.path.exists(COMPOUND_RENAME_PATH))


def LoadCitationData(citation_data_path=DEFAULT_CITATION_DATA_PATH):
    values_source_model = apps.get_model('gibbs.ValueSource')
    values_source_model.objects.all().delete()

    parsed_json = json.load(open(citation_data_path, 'rt'))

    for cd in parsed_json:
        data = json.dumps(cd)
        name = cd['name']
        source, created = values_source_model.objects.get_or_create(
                name=name, data=data)
        source.save()


def GetReactions(rids_list):
    """Find all the given reactions in the database.

    Skip those that are not present.
    """
    if not rids_list:
        return []

    reaction_model = apps.get_model('gibbs.StoredReaction')
    rxns = []
    for rid in rids_list:
        rxn = reaction_model.objects.filter(kegg_id=rid)
        rxns += list(rxn)
    return rxns


def GetCompounds(cids_list):
    """
        Find all the given compounds in the database.

        Skip those that are not present.
    """
    if not cids_list:
        return []

    compound_model = apps.get_model('gibbs.Compound')
    compounds = []
    for kegg_id in cids_list:
        try:
            compounds.append(compound_model.objects.get(kegg_id=kegg_id))
        except compound_model.DoesNotExist:
            logging.warning('Compound %s does not exist in the database', kegg_id)
            continue
    return compounds

def GetSource(source_string):
    if not source_string:
        return None

    source_model = apps.get_model('gibbs.ValueSource')
    lsource = source_string.strip().lower()
    if lsource in CITATIONS_CACHE:
        return CITATIONS_CACHE[lsource]

    try:
        source = source_model.objects.get(name__iexact=lsource)
        CITATIONS_CACHE[lsource] = source
        return source
    except source_model.DoesNotExist:
        logging.warning('Source "%s" does not exist in the database', source_string)
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

    sg, created = apps.get_model('gibbs.SpeciesGroup').objects.get_or_create(
            kegg_id=compound.kegg_id,
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

        specie, created = apps.get_model('gibbs.Specie').objects.get_or_create(
                kegg_id=compound.kegg_id,
                number_of_hydrogens=number_of_hydrogens,
                number_of_mgs=number_of_mgs,
                net_charge=net_charge,
                formation_energy=formation_energy,
                formation_energy_source=source,
                phase=phase)
        specie.save()
        sg.species.add(specie)

    sg.save()
    compound.species_groups.add(sg)


def LoadKeggCompoundNames(compound_name_path=COMPOUND_NAME_PATH,
                          compound_rename_path=COMPOUND_RENAME_PATH):

    kegg_names_df = pd.read_csv(compound_name_path, index_col=0, header=0,
                                delimiter='\t')
    cid2names = kegg_names_df['all names'].str.split('|').to_dict()
    
    # the first name on the list should be the 'common name'
    for cid, row in kegg_names_df.iterrows():
        if row['common name'] not in cid2names[cid]:
            raise ValueError('The common name \'%s\' is not in the name list for %s'
                             % (row['common name'], cid))
        # move the common name to be first on the list of "all names"
        if cid2names[cid][0] != row['common name']:
            cid2names[cid].remove(row['common name'])
            cid2names[cid] = [row['common name']] + cid2names[cid]

    kegg_rename_df = pd.read_csv(compound_rename_path, index_col=None,
                                 header=0, delimiter='\t')
    
    for _, row in kegg_rename_df.iterrows():
        if row.CID not in cid2names:
            raise ValueError('%s appears in the renaming file, but not in the KEGG list'
                             % row.CID)

        if row.command.lower() == 'remove':
            # remove 'name' from the list of names
            try:
                cid2names[row.CID].remove(row['name'])
            except ValueError:
                logging.warning('The name %s is not one of the options for %s, '
                                'so it cannot be removed' % (row['name'], row.CID))
        elif row.command.lower() == 'add':
            # put 'name' in the end of the list (or move it there if it is
            # already in the list)
            if row['name'] in cid2names[row.CID]:
                cid2names[row.CID].remove(row['name'])
            cid2names[row.CID] = cid2names[row.CID] + [row['name']]
        elif row.command.lower() == 'delete':
            del cid2names[row.CID]
        elif row.command.lower() == 'replace':
            del cid2names[row.CID]
        else:
            raise ValueError('Unknown command: %s' % row.command)

    # add the KEGG ID itself as one of the names, this enables us to
    # write down reactions using KEGG IDs in the search bar
    for cid in cid2names.keys():
        cid2names[cid].append(cid)

    cname_model = apps.get_model('gibbs.CommonName')
    compound_model = apps.get_model('gibbs.Compound')
    for compound_id, names in sorted(cid2names.items()):
        cnames = [cname_model.objects.get_or_create(name=n, enabled=True)[0]
                  for n in names]
        
        for cname in cnames:
            cname.save()

        compound, created = compound_model.objects.get_or_create(
                kegg_id=compound_id, name=cnames[0])
        compound.save()
        
        for name in cnames:
            compound.common_names.add(name)

        compound.save()

def LoadComponentContributionEnergies(cc_path=CC_PATH, priority=1):
    parsed_json = json.load(gzip.open(cc_path, 'rt'))

    compound_model = apps.get_model('gibbs.Compound')
    compound_id = ''
    for cd in parsed_json:
        compound_id = cd['CID']
        logging.debug('Adding formation energy to compound %s', compound_id)
        try:
            compound = compound_model.objects.get(kegg_id=compound_id)
        except compound_model.DoesNotExist:
            logging.warning('Compound with KEGG ID %s does not exist in database'
                            % compound_id)
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

def LoadAlbertyEnergies(alberty_path=ALBERTY_PATH, priority=2):
    parsed_json = json.load(open(alberty_path, 'rt'))
    compound_model = apps.get_model('gibbs.Compound')
    for cd in parsed_json:
        logging.debug('Adding Alberty energy to compound %s', cd['cid'])
        try:
            compound = compound_model.objects.get(kegg_id=cd['cid'])
        except compound_model.DoesNotExist:
            logging.warning('Compound with KEGG ID %s does not exist in database'
                            % cd['cid'])
            continue

        # Add the thermodynamic data.
        AddPmapToCompound(cd, compound, priority=priority)
        compound.save()

def DrawThumbnails():
    for c in apps.get_model('gibbs.Compound').objects.all():
        c.WriteStructureThumbnail()
        c.save()

def LoadKeggReactions(reactions_path=REACTION_PATH,
                      compound_rename_path=COMPOUND_RENAME_PATH):
    kegg_rename_df = pd.read_csv(compound_rename_path, index_col=None,
                                 header=0, delimiter='\t')
    cid_replace = kegg_rename_df[kegg_rename_df['command'].str.lower() == 'replace']
    cid_replace = cid_replace.set_index('CID')['name'].to_dict()

    parsed_json = json.load(gzip.open(reactions_path, 'rt'))

    reaction_model = apps.get_model('gibbs.StoredReaction')
    for rd in parsed_json:
        rid = rd['RID']

        # if one of the stoichiometric coefficients is 0, that usually
        # suggests that the reaction is not an explicit one (i.e. contains
        # template compounds, such as starch in R01790). We skip these 
        # reactions since they will not be chemically balanced in our
        # system, and we cannot estimate their dG'0 anyway.
        if 0 in map(lambda x: x[0], rd['reaction']):
            logging.warning('Skipping reaction %s: one of the reactants '
                            'has a coefficient of 0' % rid)
            continue

        for coeff_cid_pair in rd['reaction']:
            if coeff_cid_pair[1] in cid_replace:
                logging.debug('replacing %s with %s in reaction %s' %
                              (coeff_cid_pair[1],
                               cid_replace[coeff_cid_pair[1]],
                               rid))
                coeff_cid_pair[1] = cid_replace[coeff_cid_pair[1]]

        rxn = reaction_model.FromJson(rd)
        
        try:
            rxn.GenerateAttributes()
            rxn.save()
        except KeyError as e:
            logging.warning('Skipping reaction %s: %s' % (rid, str(e)))

def LoadKeggEnzymes(enzymes_path=ENZYME_PATH):
    parsed_json = json.load(gzip.open(enzymes_path, 'rt'))
    enzyme_model = apps.get_model('gibbs.Enzyme')
    cname_model = apps.get_model('gibbs.CommonName')

    for ed in parsed_json:
        ec = ed['EC']
        names = ed['names'] # add the EC number also as an optional name
        reactions = ed['reaction_ids']

        if not ec:
            logging.warning('Encountered an enzyme without an EC number.')
            continue
        if not names:
            logging.warning('Ignoring EC %s: has no common names' % ec )
            continue
        names.append(ec)
        
        reactions = GetReactions(reactions)
        if not reactions:
            logging.warning('Ignoring EC %s: no reactions found' % ec)
            continue

        # Save first so we can do many-to-many mappings.
        try:
            enz, _ = enzyme_model.objects.get_or_create(ec=ec)
        except DataError as e:
            logging.warning('Ignoring EC %s: %s' % (ec, str(e)))

        # Add names, reactions, and compound mappings.
        cnames = [cname_model.objects.get_or_create(name=n, enabled=True)[0]
                  for n in names]
        
        for n in cnames:
            enz.common_names.add(n)
        
        for r in reactions:
            enz.reactions.add(r)
        enz.save()

def GenerateCompoundThumbnails():
    compound_model = apps.get_model('gibbs.Compound')
    for compound in compound_model.objects.all():
        # thumbnail starts as None when the compound is created
        # so an empty string will mean there was an error and the thumbnail
        # cannot be created
        compound.WriteStructureThumbnail()
        compound.save()

def LoadAdditionalCompoundData(
        additional_data_path=DEFAULT_ADDITIONAL_DATA_PATH):
    parsed_json = json.load(open(additional_data_path, 'rt'))
    compound_model = apps.get_model('gibbs.Compound')
    cname_model = apps.get_model('gibbs.CommonName')

    for cd in parsed_json:
        cid = cd['CID']

        compound = compound_model.objects.get(kegg_id=cid)

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
            for n in names:
                cname, _ = cname_model.objects.get_or_create(name=n,
                                                             enabled=True)
                compound.common_names.add(cname)
        if pmaps:
            # override the pseudoisomer map that appears in the
            # kegg_compound.json file
            compound.species_groups.clear()
            for pmap in pmaps:
                priority = pmap['priority']
                AddPmapToCompound(pmap, compound, priority=priority)

        compound.save()


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
    for c in apps.get_model('gibbs.Compound').objects.all():
        d = {'name': str(c.name),
             'KEGG_ID': c.kegg_id,
             'InChI': c.inchi,
             'mass': c.mass,
             'formula': c.formula,
             'num_electrons': c.num_electrons}
        rowdicts.append(d)
    json.dump(rowdicts, gzip.open(DOWNLOADS_JSON_FNAME, 'wt'),
              sort_keys=True, indent=4)


def export_reactions(priority, name, ionic_strength, pMg, pH_list):
    csv_reaction_dict = {}
    for pH in pH_list:
        reaction_fname = DOWNLOADS_REACTION_TEMPLATE % (name, pH)
        logging.info("Writing dG'0_r for %s at pH %g to: %s"
                     % (name, pH, reaction_fname))
        csv_reaction_dict[pH] = csv.writer(gzip.open(reaction_fname, 'wt'))
        csv_reaction_dict[pH].writerow(["!!SBtab TableType='Reaction' TableName='Reaction Energies pH=5.0"])
        csv_reaction_dict[pH].writerow(["!Identifiers:kegg.reaction",
                                        "transformed standard Gibbs free energy of reaction [kJ/mol]:Mean",
                                        "transformed standard Gibbs free energy of reaction [kJ/mol]:Std",
                                        "pH",
                                        "ionic strength [M]", "temperature [K]", "!Comment"])

    for r in apps.get_model('gibbs.StoredReaction').objects.all():
        rxn = r.ToReaction()
        try:
            dG0_std = rxn.DeltaGUncertainty()
            if dG0_std is not None:
                dG0_std = round(dG0_std, 1)
        except Exception as e:
            logging.warning(str(e))
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
            except Exception as e:
                logging.warning(str(e))
                dG0_prime = None
                comment = rxn.NoDeltaGExplanation()

            row = (r.kegg_id, dG0_prime, dG0_std, pH, ionic_strength,
                   constants.DEFAULT_TEMP, comment)
            csv_reaction_dict[pH].writerow(row)


def export_compounds(priority, name, ionic_strength, pMg, pH_list):
    pseudoisomer_fname = DOWNLOADS_PSEUDOISOMER_TEMPLATE % name
    csv_pseudoisomers = csv.writer(gzip.open(pseudoisomer_fname, 'wt'))
    csv_pseudoisomers.writerow(["!MiriamID::urn:miriam:kegg.compound",
                                "!Name", "!dG0 (kJ/mol)",
                                "!nH", "!charge", "!nMg", "!Note"])

    csv_compound_dict = {}
    for pH in pH_list:
        compound_fname = DOWNLOADS_COMPOUND_TEMPLATE % (name, pH)
        csv_compound_dict[pH] = csv.writer(gzip.open(compound_fname, 'wt'))
        csv_compound_dict[pH].writerow(["!!SBtab TableType='Compound' TableName='Compound Formation Energies pH=5.0'"])
        csv_compound_dict[pH].writerow(["!Identifiers:kegg.compound",
                                        "!Compound", "transformed standard Gibbs free energy of formation [kJ/mol]",
                                        "pH", "ionic strength [M]", "temperature [K]",
                                        "!Comment"])

    logging.info("Writing chemical and biochemical formation energies for %s to: %s" %
                 (name, pseudoisomer_fname))
    for compound in apps.get_model('gibbs.Compound').objects.all():
        phase = compound.GetDefaultPhaseName()
        rows = compound.ToCSVdG0(priority, phase=phase)
        csv_pseudoisomers.writerows(rows)
        for pH in pH_list:
            aq_params = conditions.AqueousParams(pH=pH, pMg=pMg,
                                                 ionic_strength=ionic_strength)
            rows = compound.ToCSVdG0Prime(priority, aq_params=aq_params,
                                          phase=phase)
            csv_compound_dict[pH].writerows(rows)
