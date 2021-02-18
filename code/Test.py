import time
from DFOperations.calculate_DF import calculate_DF
import pandas as pd
import logging

logger = logging.getLogger('Python debug')
formatter = logging.Formatter('%(asctime)s,%(msecs)d %(name)s %(levelname)s %(message)s')
filehandler = logging.FileHandler(filename='run.log', mode='w')
filehandler.setFormatter(formatter)
logger.addHandler(filehandler)
logger.setLevel(logging.DEBUG)

logging.warning('This is a debug logging')

dataset = pd.read_csv("../CommonData/showcaseDF.csv")
logger.debug("Df is read")
start_time = time.time()
mdf = calculate_DF(dataset, logger)
logger.debug("Df is calculated")
print("--- %s seconds ---" % (time.time() - start_time))
mdf.to_csv("../CommonData/calculatedDF.csv")
pass