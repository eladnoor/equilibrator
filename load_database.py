import util.django_utils
import export_database
import logging
import numpy as np
import time, datetime

import load_additional_data
import load_citation_data
import load_kegg_json
import load_compound_mappings

from django.db import transaction

def main():
    transaction.set_autocommit(False)

    load_citation_data.CheckData()
    load_kegg_json.CheckData()
    load_additional_data.CheckData()

    logging.info('Loading citation data')
    load_citation_data.LoadCitationData()
    transaction.commit()
        
    logging.info('Loading KEGG data')
    load_kegg_json.LoadAllKeggData(draw_thumbnails=True)
    transaction.commit()
    
    logging.info('Loading corrections/additions to KEGG')
    load_additional_data.LoadAdditionalCompoundData()
    transaction.commit()
    
    logging.info('Loading compound mappings')
    load_compound_mappings.LoadEquivalentCompounds()
    transaction.commit()
    
    logging.info('Exporting database to JSON and CSV files')
    export_database.export_database()
    
if __name__ == '__main__':
    start = time.time()
    main()
    end = time.time()
    elapsed = datetime.timedelta(seconds=np.floor(end - start))
    logging.info('Elapsed loading time = %s' % str(elapsed))
