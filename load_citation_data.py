import json
import logging

import util.django_utils

from gibbs import models

DEFAULT_CITATION_DATA_FILENAME = 'data/citation_data.json'


def CheckData(filenames=(DEFAULT_CITATION_DATA_FILENAME,)):
    for json_fname in filenames:
        json.load(open(json_fname))

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

if __name__ == '__main__':
    LoadCitationData()