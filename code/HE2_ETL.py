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
        self.force_reload = False

        process_list = []
        process_list += [dict(kind='EVA', src='pipeline/pipeline', dst='pipelines')]
        process_list += [dict(kind='EVA', src='omds_well_wellop', dst='wellop_better_than_nothing',
                query='| __read__ path=omds_well_wellop | table _time, IDPad, padNum, WELL_ID, wellNum, nagWellopTechnAverageResponse')]
        process_list += [dict(kind='EVA', src='oms__ids', dst='oms_ids')]
        process_list += [dict(kind='xls', src='..\\data\\zak.xlsx', dst='pads_telemetry')]
        process_list += [dict(kind='xls', src='..\\data\\TAYL-K1_K2_Out_PI.xlsx', dst='kns_telemetry')]
        process_list += [dict(kind='xls', src='..\\data\\Visio schema.xlsx', dst='visio_schema')]
        self.process_list = process_list


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

    def move_data_from_eva_to_postgres(self, src, dst, **kwargs):
        query = kwargs.get('query', f'| __read__ path={src}')
        job = self.eva_eng.jobs.create(query_text=query, cache_ttl=60, tws=0, twf=0)
        print(job.status)
        res = job.dataset.load()
        df = pd.DataFrame(res)
        df.to_sql(name=dst, con=self.pg_eng, schema='HE2', if_exists='replace', index=False)


    def do_process_item(self, kind, src, dst, **kwargs):
        if not self.need_reload(dst):
            return
        if kind == 'xls':
            self.move_data_from_xls_to_postgres(src, dst)
        elif kind == 'EVA':
            self.move_data_from_eva_to_postgres(src, dst, **kwargs)



    def fill_pipelines_table(self):
        if not self.need_reload('pipelines'):
            return
        query = f'| __read__ path={self.pipelines_tablename}'
        self.stage_name = 'grab_pipelines_from_eva, query='+ query
        job = self.eva_eng.jobs.create(query_text=query, cache_ttl=60, tws=0, twf=0)
        print(job.status)
        res = job.dataset.load()
        df = pd.DataFrame(res)

        self.stage_name = 'upload pipelines table to postgres'
        df.to_sql(name='pipelines', con=self.pg_eng, schema='HE2', if_exists='replace', index=False)


    def fill_better_than_nothing_wellop_table(self):
        if not self.need_reload('wellop_better_than_nothing'):
            return

        # query = f'''| __read__ path=omds_well_wellop
        # | table _time, IDPad, padNum, WELL_ID, wellNum, nagWellopTechnAverageResponse'''

        self.stage_name = 'grab wellop from csv-file'
        df = pd.read_csv(self.wellop_csv_filename)

        self.stage_name = 'upload wellop table to postgres'
        df.to_sql(name='wellop_better_than_nothing', con=self.pg_eng, schema='HE2', if_exists='replace', index=False)

    def fill_wells_pads_ids_table(self):
        # В справочнике oms__ids есть поле padNum и поле pipePadNum. Второе - это номер куста в OIS Pipe.
        # P.S. Аналогично - wellNum и pipeWellNum

        if not self.need_reload('pipelines'):
            return
        query = f'| readFile format=parquet path=oms__ids'
        self.stage_name = 'grab_pipelines_from_eva, query='+ query
        job = self.eva_eng.jobs.create(query_text=query, cache_ttl=60, tws=0, twf=0)
        print(job.status)
        res = job.dataset.load()
        df = pd.DataFrame(res)

        self.stage_name = 'upload pipelines table to postgres'
        df.to_sql(name='pipelines', con=self.pg_eng, schema='HE2', if_exists='replace', index=False)

        pass

    def fill_all_neccesary_tables(self):
        for p_item in self.process_list:
            self.do_process_item(**p_item)
        # self.fill_pipelines_table()
        # self.move_data_from_xls_to_postgres(self.PQ_xls_filename, 'pads_telemetry')
        # self.move_data_from_xls_to_postgres(self.KNS_xls_filename, 'kns_telemetry')
        # self.move_data_from_xls_to_postgres(self.Visio_schema_xls_filename, 'visio_schema')
        # self.fill_better_than_nothing_wellop_table()
        # self.fill_wells_pads_ids_table()

    def main(self):
        try:
            self.fill_all_neccesary_tables()
# Два параметра на вход - какую КНС считаем 1 или 2, и дата-время
# По первому получаем rs_id, делаем выборку из pipelines. Это уже граф, но неполный
    # Выпиливаем трубы соединяющие куст и скважину
    # Выпиливаем двойные трубы по критерию даты создания или еще какому-то флагу (статусы трубопроводов)
# По телеметрии выбираем на дату-время закачку и давления по скважинам. Выкидываем аутлайеров
# Агрегируем скважины в кусты.
# Выделяем множество кустов, по которым есть телеметрия, но нет труб
    # Проверяем что эти кусты подключены на схеме.
    # Генерируем синтетические трубы
# Выделяем множество кустов, которые присоединены к сети ППД, но по ним нет телеметрии
    # Подтаскиваем по ним какие-то заглушки из режимной закачки
# Формируем и выгружаем датасет для расчета

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
    print('Finished')
