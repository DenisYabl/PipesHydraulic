import os

import pandas as pd
import sys
import logging
#from envyaml import EnvYAML

import time
from DFOperations.calculate_DF import calculate_DF
import pandas as pd
import logging
import sys

from Tools.HE2_schema_maker import make_oilpipe_schema_from_OT_dataset
from data.data_load import load_data_from_dataset

sys.path.append(os.path.abspath("."))
# create logger
from plugin import Plugin


class Worker(Plugin):
    def calc(self):
        #config = EnvYAML("params.yaml")

        dataset = self.input_dataset()
        schema = self.input_schema()

        logger = logging.getLogger('osr_hid')
        formatter = logging.Formatter('%(asctime)s,%(msecs)d %(name)s %(levelname)s %(message)s')
        streamhandler = logging.StreamHandler(sys.stderr)
        streamhandler.setFormatter(formatter)

        logger.addHandler(streamhandler)
        #logger.level = config["log"]["log_level"]

        #HKT_path = os.path.abspath(config["data"]["HKT"])
        #pump_curves_path = os.path.abspath(config["data"]["pump_curves"])
        #inclination_path = os.path.abspath(config["data"]["inclination"])

        #PumpChart = pd.read_parquet(pump_curves_path)
        #inclination = pd.read_parquet(inclination_path)
        #HKT = pd.read_parquet(HKT_path)

        df = load_data_from_dataset(dataset, schema)

        logger.debug(df)

        mdf = calculate_DF(df)
        mdf = mdf[['juncType', 'pipeline_purpose_id', 'simple_part_id', 'part_id', 'rs_schema_id', 'schema_id',
       'pipeline_id', 'node_id_end', 'node_id_start', 'L', 'simple_part_creation_date', 'node_name_start', 'altitude_start',
       'node_type_start', 'node_name_end', 'altitude_end', 'node_type_end', 'D', 'S', 'thread_number', 'uphillM', 'startIsSource', 'VolumeWater',
       'startKind', 'startValue', 'endIsOutlet', 'endKind', 'endValue', 'startP', 'startT', 'endP', 'endT', 'effectiveD', 'intD', 'roughness',
       'productivity', 'model', 'frequency', 'perforation', 'pumpDepth', 'wellNum', 'padNum']]
        mdf = mdf.fillna(0)
        mdf['startIsSource'] = mdf['startIsSource'].astype(int)
        mdf['endIsOutlet'] = mdf['endIsOutlet'].astype(int)
        output_rows = []
        for _, rows in mdf.iterrows():
            output_rows.append(
                rows.values.tolist()
            )
        self.output_dataset(output_rows)
        logger.debug(output_rows)
        return_schema = "'juncType' STRING, 'pipeline_purpose_id' STRING, 'simple_part_id' BIGINT, 'part_id' BIGINT, 'rs_schema_id' BIGINT, 'schema_id' BIGINT, 'pipeline_id' BIGINT, 'node_id_end' BIGINT, 'node_id_start' BIGINT, 'L' DOUBLE, 'simple_part_creation_date' BIGINT, 'node_name_start' STRING, 'altitude_start' DOUBLE, 'node_type_start' INT, 'node_name_end' STRING, 'altitude_end' DOUBLE, 'node_type_end' INT, 'D' DOUBLE, 'S' DOUBLE, 'thread_number' INT, 'uphillM' DOUBLE, 'startIsSource' INT, 'VolumeWater' DOUBLE, 'startKind' STRING, 'startValue' DOUBLE, 'endIsOutlet' INT, 'endKind' STRING, 'endValue' DOUBLE, 'startP' DOUBLE, 'startT' DOUBLE, 'endP' DOUBLE, 'endT' DOUBLE, 'effectiveD' DOUBLE, 'intD' DOUBLE, 'roughness' DOUBLE,'productivity' DOUBLE, 'model' STRING, 'frequency' DOUBLE, 'perforation' DOUBLE, 'pumpDepth' DOUBLE, 'wellNum' STRING, 'padNum' INT"
        self.output_schema(return_schema)
        logger.debug(return_schema)
        return 0