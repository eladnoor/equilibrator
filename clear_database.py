#!/usr/bin/python

import sys
import StringIO
from django.core.management import setup_environ, call_command

print "Setting Up Environment... ",
try:
    import settings # Assumed to be in the same directory.
except ImportError:
    sys.stderr.write("Error: Can't find the file 'settings.py' in the directory containing %r.\n"
                     "(If the file settings.py does indeed exist, it's causing an ImportError somehow.)\n" % __file__)
    sys.exit(1)

setup_environ(settings)
print "Done"

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

