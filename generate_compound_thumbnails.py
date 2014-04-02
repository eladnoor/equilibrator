import util.django_utils
import logging
import numpy as np
import time, datetime

from django.db import transaction
from gibbs import models

def GenerateCompoundThumbnails():
    transaction.set_autocommit(True)
    
    for compound in models.Compound.objects.all():
        if compound.kegg_id in ['C09078', 'C09093', 'C09145', 'C09246',
                                'C10282', 'C10286', 'C10356', 'C10359',
                                'C10396', 'C16818', 'C16839', 'C16857']:
            # these compound cause a segmentation fault in Indigo.
            # I have no idea why. The only way is to skip them.
            logging.warning('Skipping compound: %s' % compound.kegg_id)            
            continue
        
        # thumbnail starts as None when the compound is created
        # so an empty string will mean there was an error and the thumbnail
        # cannot be created
        if compound.thumbnail == '':
            logging.info('Drawing %s' % compound.kegg_id)
            compound.WriteStructureThumbnail()
            compound.save()
    
if __name__ == '__main__':
    start = time.time()
    GenerateCompoundThumbnails()
    end = time.time()
    elapsed = datetime.timedelta(seconds=np.floor(end - start))
    logging.info('Elapsed loading time = %s' % str(elapsed))
