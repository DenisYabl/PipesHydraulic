from abc import ABC
import numpy as np
import os
import pathlib
import pandas as pd
import logging
import sys
from DFOperations.calculate_DF import calculate_DF

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

        dataset = self.dataset
        schema = self.schema

        #transform input dataset to pandas

        columns = [x.split(" ")[0][1:-1] for x in schema.split(", ")]
        df = pd.DataFrame.from_records(dataset).transpose()
        df.columns = columns

        df_result = calculate_DF(df)

        return (
            schema,
            [   df_result[column].values for column in df.columns
            ],
        )