#!/usr/bin/python

import csv
import sys
import json
import logging
import django
import os
os.environ['DJANGO_SETTINGS_MODULE'] = 'equilibrator.settings'
django.setup()
from django.apps import apps
from optparse import OptionParser
from gibbs.conditions import AqueousParams
from util import constants

# Column names
KEGG_ID = '!Identifiers:kegg.compound'
NAME = '!Compound'
INCHI = '!InChI'
SOURCE = '!Source'
FORMATION = '!FormationEnergy'
ROW_ORDER = [NAME, KEGG_ID, INCHI, FORMATION, SOURCE]

def GenFormationEnergyData(aq_params):
    """Returns a list of dictionaries of compound formation energies.
    
    TODO(flamholz): return data from multiple sources per compound when possible.
    
    Args:
        pH: the pH.
        ionic_strength: the ionic strength.
    """
    dicts = []
    for compound in apps.get_model('gibbs.Compound').objects.all():
        dG = compound.DeltaG0Prime(aq_params=aq_params)
        if dG:
            dG = round(dG, 3)
        
        if dG is not None and not compound.dg_source:
            logging.error('%s has a dg but no source' % compound)
            
        d = {KEGG_ID: compound.kegg_id,
             NAME: compound.FirstName(),
             FORMATION: dG,
             INCHI: compound.inchi,
             SOURCE: None}
        if dG is not None:
            d[SOURCE] = compound.dg_source.name
        
        dicts.append(d)
    return dicts

def MakeOpts():
    """Returns an OptionParser object with all the default options."""
    opt_parser = OptionParser()
    opt_parser.add_option("-p", "--ph", dest="pH", type="float",
                          default=constants.DEFAULT_PH, help="The pH")
    opt_parser.add_option("-i", "--ionic_strength", dest="ionic_strength",
                          type="float", default=constants.DEFAULT_IONIC_STRENGTH,
                          help="The ionic strength, M")
    opt_parser.add_option("-o", "--output_name", dest="output_name",
                          type="string", default="compound_formation_energies",
                          help="The name of the file to write csv output to.")
    return opt_parser

def WriteHeader(dict_writer, row_order=ROW_ORDER):
    """writeheader() is new in Python 2.7"""
    if hasattr(dict_writer, 'writeheader'):
        dict_writer.writeheader()
    else:
        dict_writer.writer.writerow(ROW_ORDER)

def main():
    options, _ = MakeOpts().parse_args(sys.argv)
    print('Using pH = %.2f and ionic strength = %.3f' % (options.pH,
                                                         options.ionic_strength))
    output_name = '%s_pH%.2f_is%.3f' % (options.output_name, options.pH,
                                        options.ionic_strength)
    output_tsv_name = output_name + '.tsv'
    output_json_name = output_name + '.json'
    print('Will write tsv output to %s' % output_tsv_name)
    print('Will write json output to %s' % output_json_name)

    aq_params = AqueousParams(pH=options.pH,
                              ionic_strength=options.ionic_strength)
    dicts = GenFormationEnergyData(aq_params)
    sorted_data = sorted(dicts, key=lambda x: (x[KEGG_ID], x[SOURCE]))
    csv_file = open(output_tsv_name, 'w')
    writer = csv.DictWriter(csv_file, ROW_ORDER, dialect=csv.excel_tab)
    WriteHeader(writer, ROW_ORDER)
    writer.writerows(sorted_data)
    csv_file.close()
    
    json_file = open(output_json_name, 'w')
    json.dump(sorted_data, json_file, sort_keys=True, indent=4)
    json_file.close()
    print('Done.')
            
if __name__ == '__main__':
    main()


