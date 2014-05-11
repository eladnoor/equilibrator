import util.django_utils
import logging
import numpy as np
import time, datetime

from django.db import transaction
from gibbs import models
from load_kegg_json import DrawThumbnails

def GenerateCompoundThumbnails():
    
    for compound in models.Compound.objects.all():
        # thumbnail starts as None when the compound is created
        # so an empty string will mean there was an error and the thumbnail
        # cannot be created
        logging.info('Drawing %s' % compound.kegg_id)
        compound.WriteStructureThumbnail()
        compound.save()
    
if __name__ == '__main__':
    start = time.time()
    transaction.set_autocommit(True)
    DrawThumbnails()
    end = time.time()
    elapsed = datetime.timedelta(seconds=np.floor(end - start))
    logging.info('Elapsed loading time = %s' % str(elapsed))
