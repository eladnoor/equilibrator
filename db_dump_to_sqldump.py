#!/usr/bin/python
import os
import logging
import settings
import time
import datetime
import numpy
import django
import util.django_utils # keep this here, it initializes the DJANGO_SETTINGS_MODULE


def main():
    django.setup()

    logging.info('> Dumping MySQL database into sqldump file')
    db_user, db_name, db_pass = map(settings.DATABASES['default'].get, ['USER', 'NAME', 'PASSWORD'])
    os.environ['MYSQL_PWD'] = db_pass
    cmd = "mysqldump -u %s %s | gzip -c > data/sqldump.txt.gz" % (db_user, db_name)
    os.system(cmd)

if __name__ == '__main__':
    logging.info('Welcome to the sql_load script')

    start = time.time()
    main()
    end = time.time()
    elapsed = datetime.timedelta(seconds=numpy.floor(end - start))
    logging.info('Elapsed loading time = %s' % str(elapsed))
