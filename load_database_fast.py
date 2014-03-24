import util.django_utils
import logging
import timeit

import load_additional_data
import load_citation_data
import load_kegg_json
import load_compound_mappings

from django.db import transaction

def main():
    transaction.set_autocommit(False)
    start = timeit.timeit()
    
    load_citation_data.CheckData()
    load_kegg_json.CheckData()
    load_additional_data.CheckData()
    transaction.commit()

    logging.info('Loading citation data')
    load_citation_data.LoadCitationData()
    transaction.commit()
        
    logging.info('Loading KEGG data')
    load_kegg_json.LoadAllKeggData(draw_thumbnails=False)
    transaction.commit()
    
    logging.info('Loading corrections/additions to KEGG')
    load_additional_data.LoadAdditionalCompoundData()
    transaction.commit()
    
    logging.info('Loading compound mappings')
    load_compound_mappings.LoadEquivalentCompounds()
    transaction.commit()
    
    end = timeit.timeit()
    logging.info('Elapsed loading time = %s' % str(end - start))
    
if __name__ == '__main__':
    main()
