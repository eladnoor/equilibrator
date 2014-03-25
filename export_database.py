from util import django_utils
import csv
import gzip
import json
import logging

django_utils.SetupDjango()

from gibbs import models
from gibbs import constants

COMPOUND_PREFIX = 'media/downloads/kegg_compounds'
REACTION_PREFIX = 'media/downloads/kegg_reactions'
PSEUDOISOMER_PREFIX = 'media/downloads/kegg_pseudoisomers'
JSON_FNAME = 'media/downloads/kegg_compounds.json.gz'

def export_database():
    
    export_json()
    
    for priority, name in [(1, 'UGC'), (2, 'alberty')]:

        export_compounds(priority=priority, name=name,
                         ionic_strength=constants.DEFAULT_IONIC_STRENGTH,
                         pMg=constants.DEFAULT_PMG,
                         pH_list=constants.PH_RANGE_VALUES)
        
        export_reactions(priority=priority, name=name,
                         ionic_strength=constants.DEFAULT_IONIC_STRENGTH,
                         pMg=constants.DEFAULT_PMG,
                         pH_list=constants.PH_RANGE_VALUES)

def export_json():
    logging.info("Writing compound data to JSON file: %s" % JSON_FNAME)
    rowdicts = []
    for c in models.Compound.objects.all():
        d = {'name': str(c.name), 
             'KEGG_ID': c.kegg_id,
             'InChI': c.inchi,
             'mass': c.mass,
             'formula': c.formula,
             'num_electrons': c.num_electrons}
        rowdicts.append(d)
    json.dump(rowdicts, gzip.open(JSON_FNAME, 'w'), sort_keys=True, indent=4)

def export_reactions(priority, name, ionic_strength, pMg, pH_list):
    csv_reaction_dict = {}
    for pH in pH_list:
        reaction_fname = REACTION_PREFIX + '_%s_ph%.1f.csv.gz' % (name, pH)
        logging.info("Writing transformed reaction energies for %s at pH %g to: %s"
                     % (name, pH, reaction_fname))
        csv_reaction_dict[pH] = csv.writer(gzip.open(reaction_fname, 'w'))
        csv_reaction_dict[pH].writerow(["!MiriamID::urn:miriam:kegg.reaction",
                                        "!dG0_prime (kJ/mol)", "!pH",
                                        "!I (mM)", "!T (Kelvin)", "!Note"])
    
    for r in models.StoredReaction.objects.all():
        for pH in pH_list:
            rows = r.ToCSVdG0Prime(priority, pH=pH, pMg=pMg,
                                   ionic_strength=ionic_strength)
            for row in rows:
                csv_reaction_dict[pH].writerow(row)
                
def export_compounds(priority, name, ionic_strength, pMg, pH_list):
    pseudoisomer_fname = PSEUDOISOMER_PREFIX + '_%s.csv.gz' % name
    logging.info("Writing chemical formation energies for %s to: %s" %
                 (name, pseudoisomer_fname))
    csv_pseudoisomers = csv.writer(gzip.open(pseudoisomer_fname, 'w'))
    csv_pseudoisomers.writerow(["!MiriamID::urn:miriam:kegg.compound",
                                "!Name", "!dG0 (kJ/mol)",
                                "!nH", "!charge", "!nMg", "!Note"])

    csv_compound_dict = {}
    for pH in pH_list:
        compound_fname = COMPOUND_PREFIX + '_%s_ph%.1f.csv.gz' % (name, pH)
        logging.info("Writing transformed formation energies for %s at pH %g to: %s" %
                     (name, pH, compound_fname))
        csv_compound_dict[pH] = csv.writer(gzip.open(compound_fname, 'w'))
        csv_compound_dict[pH].writerow(["!MiriamID::urn:miriam:kegg.compound",
                                        "!Name", "!dG0_prime (kJ/mol)",
                                        "!pH", "!I (mM)", "!T (Kelvin)",
                                        "!Note"])
    
    for c in models.Compound.objects.all():
        for row in c.ToCSVdG0(priority):
            csv_pseudoisomers.writerow(row)
        for pH in pH_list:
            for row in c.ToCSVdG0Prime(priority, pH=pH, pMg=pMg,
                                       ionic_strength=ionic_strength):
                csv_compound_dict[pH].writerow(row)

if __name__ == '__main__':
    export_database()