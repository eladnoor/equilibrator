#!/usr/bin/python3
import os
import logging
import time
import datetime
import numpy
from equilibrator.settings import DATABASES, USE_SOLR
import django
from django.core.management import execute_from_command_line

def main():
    os.environ.setdefault("DJANGO_SETTINGS_MODULE", "equilibrator.settings")
    django.setup()

    logging.getLogger().setLevel(logging.INFO)

    logging.info('> Flushing MySQL data')
    execute_from_command_line(['', 'flush', '--noinput'])

    logging.info('> Loading data from sqldump into MySQL')
    db_user, db_name, db_pass = map(DATABASES['default'].get, ['USER', 'NAME', 'PASSWORD'])
    os.environ['MYSQL_PWD'] = db_pass
    cmd = "gunzip -c data/sqldump.txt.gz | mysql -u %s %s" % (db_user, db_name)
    os.system(cmd)

    if USE_SOLR:
        logging.info('> Clearing Solr index\n')
        execute_from_command_line(['', 'clear_index', '--noinput'])
    
        logging.info('> Building Solr index\n')
        execute_from_command_line(['', 'update_index'])

if __name__ == '__main__':
    logging.info('Welcome to the sql_load script')

    start = time.time()
    main()
    end = time.time()
    elapsed = datetime.timedelta(seconds=numpy.floor(end - start))
    logging.info('Elapsed loading time = %s' % str(elapsed))
