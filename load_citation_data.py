import json
import logging

from util import django_utils

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
            name = cd['name']
            year = int(cd['year'])
            ref = cd['ref']
            link = cd.get('url', None)
            doi = cd.get('doi', None)
            pubmed_id = cd.get('pmid', None)
            
            logging.info('Loading source "%s"', name)
            source = models.ValueSource(name=name, citation=ref, year=year,
                                        pubmed_id=pubmed_id, doi=doi, link=link)
            source.save()
        except Exception, e:
            logging.error('Error parsing reference %s', cd)
            logging.error(e)
            continue

if __name__ == '__main__':
    LoadCitationData()