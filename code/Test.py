import time
from DFOperations.calculate_DF import calculate_DF
import pandas as pd
import logging
import sys


dataset = pd.read_csv("../CommonData/showcaseDF.csv")
start_time = time.time()
mdf = calculate_DF(dataset)
print("--- %s seconds ---" % (time.time() - start_time))
mdf.to_csv("../CommonData/calculatedDF.csv")
pass