import streamlit as st
import time
import pandas as pd
import logging
import sys


from DFOperations.calculate_DF import calculate_DF

dataset = pd.read_csv("../CommonData/1well.csv")
print(dataset)
start_time = time.time()
mdf = calculate_DF(dataset)
print("--- %s seconds ---" % (time.time() - start_time))
mdf.to_csv("../CommonData/mdf.csv")


