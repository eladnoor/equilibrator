#!/usr/local/bin/python

import sys
import StringIO
from django.core.management import call_command

from util import django_utils

# Now we can import Django stuff!
from django.db import connection

print "Deleting Tables for gibbs... ",
sys.stdout = buff = StringIO.StringIO() # temporarily redirect the stdout to a stored buffer
call_command('sqlclear', 'gibbs')
sys.stdout = sys.__stdout__

queries = buff.getvalue().split(';')[1:-2]
cursor = connection.cursor()
for query in queries:
    print query,
    cursor.execute(query.strip())
print "Done"

print "Sync-ing Database... ",
call_command('syncdb')
print "Done"

