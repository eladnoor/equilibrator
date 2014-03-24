#!/usr/local/bin/python

from util import django_utils
import logging

import load_additional_data
import load_citation_data
import load_kegg_json
import load_compound_mappings

def main():
    load_citation_data.CheckData()
    load_kegg_json.CheckData()
    load_additional_data.CheckData()

    logging.info('Loading citation data')
    load_citation_data.LoadCitationData()
        
    logging.info('Loading KEGG data')
    load_kegg_json.LoadAllKeggData(draw_thumbnails=False)
    
    logging.info('Loading corrections/additions to KEGG')
    load_additional_data.LoadAdditionalCompoundData()
    
    logging.info('Loading compound mappings')
    load_compound_mappings.LoadEquivalentCompounds()
    
if __name__ == '__main__':
    main()
