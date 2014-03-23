import json
import logging

from util import django_utils
import gzip

from gibbs import models

DEFAULT_COMPOUND_DATA_FILENAME = 'data/additional_compound_data.json.gz'


def GetOrCreateNames(names_list):
    """
        Find all the names in the database.
        Create them if they are not present.
    """
    return [models.CommonName.GetOrCreate(n)
            for n in names_list]
    

def LoadAdditionalCompoundData(json_filename=DEFAULT_COMPOUND_DATA_FILENAME):
    parsed_json = json.load(gzip.open(json_filename))

    for cd in parsed_json:
        try:
            cid = cd['CID']
            
            compound = models.Compound.objects.get(kegg_id=cid)
            
            note = cd.get('note')
            preferred_name = cd.get('preferred name')
            details_link = cd.get('details_link')
            if note:
                compound.note = note
            if preferred_name:
                compound.preferred_name = preferred_name
            if details_link:
                compound.details_link = details_link
            
            names = cd.get('names')
            if names:
                for n in GetOrCreateNames(names):
                    compound.common_names.add(n)
                    
            compound.save()
        except Exception, e:
            logging.error('Error parsing cid %s', cid)
            logging.error(e)
            continue


def CheckData(filenames=(DEFAULT_COMPOUND_DATA_FILENAME,)):
    for json_fname in filenames:
        json.load(gzip.open(json_fname))


def LoadAllAdditionalData():
    LoadAdditionalCompoundData()

    
if __name__ == '__main__':
    LoadAdditionalCompoundData()