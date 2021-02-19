from abc import ABC
import numpy as np
import os
import pathlib
import pandas as pd
import logging
import sys
from DFOperations.calculate_DF import calculate_DF
import logging

ROOT_DIR = os.path.dirname(os.path.abspath(__file__))

# create logger


class Plugin(ABC):
    def __init__(self):
        self.schema = None
        self.dataset = None

    def load_data(self, schema, dataset):
        self.schema = schema
        self.dataset = dataset

    def calc(self):
        pass

class Worker(Plugin):
    def calc(self):
        logger = logging.getLogger('Python debug')
        formatter = logging.Formatter('%(asctime)s,%(msecs)d %(name)s %(levelname)s %(message)s')
        filehandler = logging.FileHandler(filename='run.log', mode='w')
        filehandler.setFormatter(formatter)
        streamhandler = logging.StreamHandler(sys.stderr)
        streamhandler.setFormatter(formatter)

        logger.addHandler(filehandler)
        logger.addHandler(streamhandler)
        logger.setLevel(logging.DEBUG.encode())

        dataset = self.dataset
        schema = self.schema

        #transform input dataset to pandas
        logger.debug("Data is acquired")
        columns = [x.split(" ")[0][1:-1] for x in schema.split(", ")]
        df = pd.DataFrame.from_records(dataset).transpose()
        df.columns = columns
        logger.debug("Pandas DF is created".encode())

        df_result = calculate_DF(df, logger=logger)
        logger.debug("DF is calculated".encode())
        return (
            schema,
            [   df_result[column].values for column in df.columns
            ],
        )