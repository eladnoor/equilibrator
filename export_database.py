import logging, time, datetime
import numpy as np
from util import database_io as db

def main():
    logging.info('Exporting database to JSON and CSV files')
    db.export_database()
                
if __name__ == '__main__':
    start = time.time()
    main()
    end = time.time()
    elapsed = datetime.timedelta(seconds=np.floor(end - start))
    logging.info('Elapsed loading time = %s' % str(elapsed))
