import logging, time, datetime, sys
import numpy as np
from util import database_io as db
from django.db import transaction
from distutils.util import strtobool

def main(draw_thumb=False, export_csv=False):
    transaction.set_autocommit(False)

    db.CheckData()

    logging.info('Loading citation data')
    db.LoadCitationData()
    transaction.commit()
        
    #LoadKeggGCNullspace()

    logging.info('Loading KEGG compound names')
    db.LoadKeggCompoundNames()
    transaction.commit()

    logging.info('Loading KEGG compound thermodynamic data')
    db.LoadFormationEnergies()
    transaction.commit()

    logging.info('Loading KEGG reaction data')
    db.LoadKeggReactions()
    transaction.commit()

    logging.info('Loading KEGG enzyme data')
    db.LoadKeggEnzymes()
    transaction.commit()

    if draw_thumb:
        logging.info('Drawing thumbnails for all KEGG compounds')
        db.GenerateCompoundThumbnails()
        transaction.commit()
    
    logging.info('Loading corrections/additions to KEGG')
    db.LoadAdditionalCompoundData()
    transaction.commit()
    
    if export_csv:
        logging.info('Exporting database to JSON and CSV files')
        db.export_database()
    
def user_yes_no_query(question, default=False):
    if default:
        sys.stdout.write('%s? [(yes)/no] ' % question)
    else:
        sys.stdout.write('%s? [yes/(no)] ' % question)
    
    while True:
        try:
            ri = raw_input().lower()
            if not ri:
                return default
            return strtobool(ri)
        except ValueError:
            sys.stdout.write('Please respond with \'y\' or \'n\'.\n')
            
if __name__ == '__main__':
    print 'Welcome to the load_database script.'
    print 'Drawing thumbnails for all the compounds takes about 30 minutes'
    draw_thumb = user_yes_no_query('Draw thumbnails')

    print 'Exporting the raw data as CSV files takes about 1 hour'
    export_csv = user_yes_no_query('Export raw data')
    
    start = time.time()
    main(draw_thumb, export_csv)
    end = time.time()
    elapsed = datetime.timedelta(seconds=np.floor(end - start))
    logging.info('Elapsed loading time = %s' % str(elapsed))
