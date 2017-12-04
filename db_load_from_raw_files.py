#!/usr/bin/python3
import logging
import time
import datetime
import sys
import numpy
from util import database_io
from django.db import transaction
from distutils.util import strtobool
from django.core.management import execute_from_command_line

def main(draw_thumb=False, export_csv=False):
    logging.info('> Flushing SQL data')
    execute_from_command_line(['', 'flush', '--noinput'])

    transaction.set_autocommit(False)

    database_io.CheckData()

    logging.info('> Loading citation data')
    database_io.LoadCitationData()
    transaction.commit()

    logging.info('> Loading KEGG compound names')
    cid_replace = database_io.LoadKeggCompoundNames()
    transaction.commit()

    logging.info('> Loading KEGG compound thermodynamic data')
    database_io.LoadFormationEnergies()
    database_io.LoadAlbertyEnergies()
    transaction.commit()

    logging.info('> Loading KEGG reaction data')
    database_io.LoadKeggReactions(cid_replace)
    transaction.commit()

    logging.info('> Loading KEGG enzyme data')
    database_io.LoadKeggEnzymes()
    transaction.commit()

    if draw_thumb:
        logging.info('> Drawing thumbnails for all KEGG compounds')
        database_io.GenerateCompoundThumbnails()
        transaction.commit()

    logging.info('> Loading corrections/additions to KEGG')
    database_io.LoadAdditionalCompoundData()
    transaction.commit()

    if export_csv:
        logging.info('> Exporting database to JSON and CSV files')
        database_io.export_database()

    logging.info('> Clearing Solr index\n')
    execute_from_command_line(['', 'clear_index', '--noinput'])

    logging.info('> Building Solr index\n')
    execute_from_command_line(['', 'update_index'])

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
    logging.info('Welcome to the load_database script')
    logging.info('Drawing thumbnails for all the compounds takes about 30 minutes')
    draw_thumb = user_yes_no_query('Draw thumbnails', default=True)

    logging.info('Exporting the raw data as CSV files takes about 1 hour')
    export_csv = user_yes_no_query('Export raw data', default=False)

    start = time.time()
    main(draw_thumb, export_csv)
    end = time.time()
    elapsed = datetime.timedelta(seconds=numpy.floor(end - start))
    logging.info('Elapsed loading time = %s' % str(elapsed))
