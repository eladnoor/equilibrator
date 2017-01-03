#!/usr/bin/python
from settings import DATABASES
import os

db_user = DATABASES['default']['USER']
db_name = DATABASES['default']['NAME']

cmd = "mysqldump -u %s -p %s | gzip -c > data/sqldump.txt.gz" % (db_user, db_name)

print cmd
os.system(cmd)
