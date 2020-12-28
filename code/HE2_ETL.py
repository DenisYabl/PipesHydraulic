import psycopg2
import pandas as pd
from sqlalchemy import create_engine
from ot_simple_connector.connector import Connector
from datetime import datetime
import numpy as np

class HE2_ETL():
    """
    This class perform all steps to get all necessary data from sources and form properly dataframe for HE2 solve run
    These steps are:
    1. Get data from heterogeneous sources and put it to local postgres tables:
    1.1. EVA pipeline/pipeline
    1.2. Pad telemetry xls file
    1.3. KNS telemetry xls file
    1.4. Pipelines and pads schema xls file, includes Visio drawning
    :return:
    """

    def __init__(self, postgres_credentials, eva_credentials):
        self.pg_eng = create_engine(postgres_credentials)
        self.eva_eng = Connector(**eva_credentials)
        self.pipelines_tablename = 'pipeline/pipeline'
        self.PQ_xls_filename = '..\\data\\zak.xlsx'
        self.KNS_xls_filename = '..\\data\\TAYL-K1_K2_Out_PI.xlsx'
        self.Visio_schema_xls_filename = '..\\data\\Visio schema.xlsx'
        self.force_reload = False
        self.stage_name = ''

    def table_exists(self, table_name):
        rez = True
        try:
            probe_df = pd.read_sql(f'Select * from "HE2".{table_name} limit 10', con=self.pg_eng)
        except Exception as e:
            rez = False
        return rez

    def need_reload(self, table_name):
        rez = True
        if not self.force_reload:
            rez = not self.table_exists(table_name)
        return rez

    def move_data_from_xls_to_postgres(self, filename, tablename):
        if not self.need_reload(tablename):
            return
        self.stage_name = f'grab xls file, filename = {filename}'
        df = pd.read_excel(filename, header=0, engine='openpyxl', index_col=None, parse_dates=True,usecols='A:G,I', verbose=True, na_values=['', ' ', '\n'],
                dtype=dict(OIS_ID=np.int32, IDSKV=np.int32, Mestor=object, kust=object, skv=object, ValueDATE=datetime, DVALUE=np.float, priz=object))
        self.stage_name = f'upload {tablename} table to postgres'
        df.to_sql(name=tablename, con=self.pg_eng, schema='HE2', if_exists='replace', index=False)


    def grab_pipelines_from_eva(self):
        query = f'| __read__ path={self.pipelines_tablename}'
        self.stage_name = 'grab_pipelines_from_eva, query='+ query
        job = self.eva_eng.jobs.create(query_text=query, cache_ttl=60, tws=0, twf=0)
        print(job.status)
        res = job.dataset.load()
        df = pd.DataFrame(res)
        return df

    def fill_pipelines_table(self):
        if not self.need_reload('pipelines'):
            return

        df = self.grab_pipelines_from_eva()
        self.stage_name = 'upload pipelines table to postgres'
        df.to_sql(name='pipelines', con=self.pg_eng, schema='HE2', if_exists='replace', index=False)

    def main(self):
        try:
            self.fill_pipelines_table()
            self.move_data_from_xls_to_postgres(self.PQ_xls_filename, 'pads_telemetry')
            self.move_data_from_xls_to_postgres(self.KNS_xls_filename, 'kns_telemetry')
            self.move_data_from_xls_to_postgres(self.Visio_schema_xls_filename, 'visio_schema')
        except Exception as e:
            print(f'Something wrong on stage {self.stage_name}')
            raise e

if __name__ == '__main__':
    # engine = create_engine('postgresql://db_user:210106@localhost:5432/postgres')
    # test_df = pd.read_csv('..\\data\\input_df.csv')
    # test_df.to_sql(name='test',con=engine, schema='HE2', index=False)
    # conn =  psycopg2.connect(dbname='postgres', user='db_user', password='210106', host='localhost')

    pipeline = HE2_ETL('postgresql://db_user:210106@localhost:5432/postgres', dict(host="192.168.4.65", port='80', user="am", password="12345678"))
    pipeline.main()
