import logging
import os
import numpy as np

logging.basicConfig(level=logging.DEBUG, filename=f'HE2 pid {os.getpid()}.log', format='%(asctime)s %(levelname)s %(funcName)s(): %(message)s')
np.seterr(invalid='raise')

def getLogger(name):
    return logging.getLogger(name)

def check_for_nan(**kwargs):
    pass
