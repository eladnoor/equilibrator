#!/usr/bin/python3
from equilibrator.settings import DATABASES, HAYSTACK_BACKEND
import os
import logging
import time
import datetime
import numpy
import django
from django.core.management import execute_from_command_line
from django.db import transaction
import argparse

def MakeParser():
    parser = argparse.ArgumentParser(
        description=('Initialize eQuilibrator MySQL database'))
    parser.add_argument('--raw', action='store_true',
                        help='Load from raw data files')
    parser.add_argument('--export_csv', action='store_true',
                        help='export final database to CSV file')
    parser.add_argument('--draw_thumb', action='store_true',
                        help='draw chemical structure thumbnails')
    return parser

def main():
    parser = MakeParser()
    args = parser.parse_args()

    os.environ.setdefault("DJANGO_SETTINGS_MODULE", "equilibrator.settings")
    django.setup()
    logging.getLogger().setLevel(logging.INFO)

    db_user, db_name, db_pass = map(DATABASES['default'].get, ['USER', 'NAME', 'PASSWORD'])
    os.environ['MYSQL_PWD'] = db_pass

    logging.info('> Drop MySQL database')
    cmd = "mysql -u %s -e 'DROP DATABASE IF EXISTS %s;'" % (db_user, db_name)
    os.system(cmd)

    logging.info('> Creating MySQL database')
    cmd = "mysql -u %s -e 'CREATE DATABASE %s;'" % (db_user, db_name)
    os.system(cmd)
        
    logging.info('> Creating tables if they do not exist')
    execute_from_command_line(['', 'migrate', '--noinput', '--run-syncdb', '-v', 0])

    logging.info('> Flushing MySQL data')
    execute_from_command_line(['', 'flush', '--noinput', '-v', 0])

    if args.raw:
        load_from_raw_files(args.draw_thumb, args.export_csv)
    else:
        load_from_sqldump(db_user, db_name)

    if HAYSTACK_BACKEND == 'solr':
        logging.info('> Clearing Solr index\n')
        execute_from_command_line(['', 'clear_index', '--noinput'])
    
        logging.info('> Building Solr index\n')
        execute_from_command_line(['', 'update_index'])

def load_from_sqldump(db_user, db_name):
    logging.info('> Loading data from sqldump into MySQL')
    cmd = "gunzip -c data/sqldump.txt.gz | mysql -u %s %s" % (db_user, db_name)
    os.system(cmd)

def load_from_raw_files(draw_thumb, export_csv):
    from util import database_io
    transaction.set_autocommit(False)

    database_io.CheckData()

    logging.info('> Loading citation data')
    database_io.LoadCitationData()
    transaction.commit()

    logging.info('> Loading KEGG compound names')
    database_io.LoadKeggCompoundNames()
    transaction.commit()

    logging.info('> Loading Component Contribution data')
    database_io.LoadComponentContributionEnergies()
    transaction.commit()

    logging.info('> Loading Alberty data')
    database_io.LoadAlbertyEnergies()
    transaction.commit()

    logging.info('> Loading KEGG reaction data')
    database_io.LoadKeggReactions()
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


if __name__ == '__main__':
    logging.getLogger().setLevel(logging.INFO)
    logging.info('Welcome to the init_db script')

    start = time.time()
    main()
    end = time.time()
    elapsed = datetime.timedelta(seconds=numpy.floor(end - start))
    logging.info('Elapsed loading time = %s' % str(elapsed))
