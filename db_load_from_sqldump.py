#!/usr/bin/python
import os
import logging
import settings
import time
import datetime
import numpy
import django
import util.django_utils # keep this here, it initializes the DJANGO_SETTINGS_MODULE
from django.core.management import execute_from_command_line


def main():
    django.setup()
    logging.info('> Flushing SQL data')
    execute_from_command_line(['', 'flush', '--noinput'])

    logging.info('> Loading data from sqldump into MySQL')
    db_user, db_name, db_pass = map(settings.DATABASES['default'].get, ['USER', 'NAME', 'PASSWORD'])
    os.environ['MYSQL_PWD'] = db_pass
    cmd = "gunzip -c data/sqldump.txt.gz | mysql -u %s %s" % (db_user, db_name)
    os.system(cmd)

    logging.info('> Rebuilding Solr index\n')
    execute_from_command_line(['', 'rebuild_index', '--noinput'])

if __name__ == '__main__':
    logging.info('Welcome to the sql_load script')

    start = time.time()
    main()
    end = time.time()
    elapsed = datetime.timedelta(seconds=numpy.floor(end - start))
    logging.info('Elapsed loading time = %s' % str(elapsed))
