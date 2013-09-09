import logging
import sys

formatter = logging.Formatter('%(asctime)s %(filename)s %(lineno)s %(levelname)s  %(message)s')

stdout_handler = logging.StreamHandler(sys.stdout)
stdout_handler.setFormatter(formatter)

logger = logging.getLogger('')
logger.addHandler(stdout_handler)
logger.setLevel(logging.INFO)