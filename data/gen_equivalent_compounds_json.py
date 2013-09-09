#!/usr/bin/python

import csv
import json
import logging


def main():
    equivalence_classes = {}
    names_by_id = {}
    reader = csv.reader(open('data/equivalent_compounds.csv'))
    
    logging.info('Reading equivalence data.')
    for row in reader:
        compound_a = row[1]
        compound_b = row[3]
        name_a = row[0]
        name_b = row[2]
         
        equivalence_classes.setdefault(compound_a, set()).add(compound_b)
        equivalence_classes.setdefault(compound_b, set()).add(compound_a)

        names_by_id.setdefault(compound_a, set()).add(name_a)
        names_by_id.setdefault(compound_b, set()).add(name_b)
    
    logging.info('Reducing equivalence data.')
    seen_compounds = set()
    reduced_classes = []
    
    for key, eq_class in equivalence_classes.iteritems():
        if key in seen_compounds:
            continue
        
        cur_class = set(eq_class)
        cur_class.add(key)
        for eq_item in eq_class:
            cur_class.update(equivalence_classes.get(eq_item, set()))
        
        reduced_classes.append(cur_class)
        seen_compounds.add(key)
        seen_compounds.update(eq_class)
    
    logging.info('Formatting for JSON output.')
    json_data = []
    for eq_class in reduced_classes:
        names = []
        for i in eq_class:
            names.extend(names_by_id[i])
            
        d = {'KEGG IDS': list(eq_class),
             'Names': names}
        json_data.append(d)
    
    json_filename = 'data/equivalent_compounds.json'
    logging.info('Writing to JSON file %s', json_filename)
    json.dump(json_data, open(json_filename, 'w'),
              sort_keys=True, indent=4)
    
if __name__ == '__main__':
    main()
        