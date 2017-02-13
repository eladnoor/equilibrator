#!/usr/bin/python
import os
from settings import DATABASES

db_user = DATABASES['default']['USER']
db_name = DATABASES['default']['NAME']

cmd = "gunzip -c data/sqldump.txt.gz | mysql -u %s -D %s -p" % (db_user, db_name)

print cmd
os.system(cmd)

print "Don't forget to run 'python manage.py rebuild_index' now"
