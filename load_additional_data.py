import json
import logging

import util.django_utils
from load_kegg_json import GetOrCreateNames, AddPmapToCompound

from gibbs import models

DEFAULT_ADDITIONAL_DATA_FILENAME = 'data/additional_compound_data.json'

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
                for n in GetOrCreateNames(names):
                    compound.common_names.add(n)
            if pmaps:
                # override the pseudoisomer map that appears in the
                # kegg_compound.json file
                compound.species_groups.clear()
                for pmap in pmaps:
                    AddPmapToCompound(pmap, compound)
                    
            compound.save()
        except Exception, e:
            logging.error('Error parsing cid %s', cid)
            logging.error(e)
            continue


def CheckData(filenames=(DEFAULT_ADDITIONAL_DATA_FILENAME,)):
    for json_fname in filenames:
        json.load(open(json_fname, 'r'))


def LoadAllAdditionalData():
    LoadAdditionalCompoundData()

    
if __name__ == '__main__':
    LoadAdditionalCompoundData()