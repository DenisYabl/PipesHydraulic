import logging
import os
import numpy as np
from Tools.HE2_ABC import SOLVER_VERSION

# logging.basicConfig(level=logging.INFO, filename=f'HE2 pid {os.getpid()}.log', format='%(asctime)s %(levelname)s %(funcName)s(): %(message)s')
logging.basicConfig(level=logging.DEBUG, filename=f'HE2 pid {os.getpid()}.log', format='%(asctime)s %(levelname)s %(funcName)s(): %(message)s')
np.seterr(invalid='raise')

def getLogger(name):
    return logging.getLogger(name)

def check_for_nan(**kwargs):
    return
    arr = np.array([v for v in kwargs.values()], dtype=float)
    if np.isnan(arr).any():
        msg = ', '.join([f'{k} = {v}' for k, v in kwargs.items()])
        raise ValueError(msg)

logger = getLogger('')
logger.error(f'Solver version is {SOLVER_VERSION}')
logger.error(f'Numpy version is {np.version.version}')
