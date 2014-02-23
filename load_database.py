#!/usr/bin/python

from util import django_utils
import export_database
import logging

django_utils.SetupDjango()

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
    draw_thumbnails = True
    ### uncomment this to make things faster (don't draw compounds):
    draw_thumbnails = False
    load_kegg_json.LoadAllKeggData(draw_thumbnails=draw_thumbnails)
    
    logging.info('Loading corrections/additions to KEGG')
    load_additional_data.LoadAdditionalCompoundData()    
    
    logging.info('Loading compound mappings')
    load_compound_mappings.LoadEquivalentCompounds()
    
    logging.info('Exporting database to JSON and CSV files')
    export_database.export_database()
    
if __name__ == '__main__':
    main()
