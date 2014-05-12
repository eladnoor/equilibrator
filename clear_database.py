import sys
import StringIO
from django.core.management import call_command

import util.django_utils

# Now we can import Django stuff!
from django.db import connection
from django.db import transaction
import logging
logger = logging.getLogger(__name__)

@transaction.atomic
def delete_all_tables():
    logger.debug("Deleting Tables for gibbs... ")
    sys.stdout = buff = StringIO.StringIO() # temporarily redirect the stdout to a stored buffer
    call_command('sqlclear', 'gibbs')
    sys.stdout = sys.__stdout__
    
    queries = buff.getvalue().split(';')[1:-2]
    cursor = connection.cursor()
    for query in queries:
        print query,
        cursor.execute(query.strip())

    logger.debug("Sync-ing Database... ")
    call_command('syncdb')


if __name__ == '__main__':
    delete_all_tables()