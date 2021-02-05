import pandas as pd
from sqlalchemy import create_engine
from ot_simple_connector.connector import Connector
from datetime import datetime
import numpy as np
from matplotlib import pyplot

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
        self.pg_eng = create_engine(postgres_credentials, connect_args={'options': '-c search_path=dbo,public,"HE2"'})
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
            probe_df = pd.read_sql(f'Select * from {table_name} limit 10', con=self.pg_eng)
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


    def fill_all_neccesary_tables(self):
        for p_item in self.process_list:
            self.do_process_item(**p_item)

    def get_wells_need_preprocess(self):
        query = f'''
            With
            q0 as (select kust, skv, "ValueDATE" dt, "DVALUE" pitch from pads_telemetry where priz='Q'),
            q1 as (select * from q0 where pitch <> FLOOR(pitch)),
            q2 as (select min(pitch) min_q, max(pitch) max_q, count(*) cnt_q, skv from q1 group by skv),
            q4 as (select distinct skv, kust from q1),
            q5 as (select q4.*, min_q, max_q, cnt_q from q4 left join q2 on q4.skv=q2.skv),
            good2 as (select * from q5 where max_q/min_q <= 2 and cnt_q >= 5)
            select skv from q5 except select skv from good2        
        '''
        df = pd.read_sql(query, self.pg_eng)
        return df

    def get_well_Q(self, well):
        query = f'''
            select "ValueDATE" dt, "DVALUE" pitch from pads_telemetry where priz='Q' and skv='{well}'             
        '''
        df = pd.read_sql(query, self.pg_eng)
        return df

    def estimate_outliers_by_window(self, pitches, window_width=5):
        n = len(pitches)
        ww = window_width
        hww = (window_width - 1) // 2
        mx = np.zeros([ww, n + 2 * hww])
        mx[:, :ww] = pitches[0]
        mx[:, -ww:] = pitches[-1]
        for i in range(ww):
            mx[i, i:i+n] = pitches
        mx = mx[:, hww:-hww]
        m = np.mean(mx,axis=0)
        sigma = np.std(mx, axis=0)
        u = m + 2*sigma
        d = m - 2*sigma
        mask1 = pitches > u
        mask2 = pitches < d
        mask = mask1 | mask2
        return sum(mask)


    def preprocess_wells_Q(self):
        df = self.get_wells_need_preprocess()
        wells_list = list(df.skv.values)
        for well in wells_list:
            df = self.get_well_Q(well)
            n = len(df)
            if n <= 5:
                print('I dont know what to do')
                continue
            out_cnt = self.estimate_outliers_by_window(df.pitch.values)
            print(well, out_cnt)



    def execute_query(self, query):
        conn = self.pg_eng.connect()
        conn.execute(query)
        conn.close()

    # По телеметрии выбираем на дату-время закачку и давления по скважинам. Выкидываем аутлайеров
    # Агрегируем скважины в кусты.
    # Выделяем множество кустов, по которым есть телеметрия, но нет труб
    # Проверяем что эти кусты подключены на схеме.
    # Генерируем синтетические трубы
    # Выделяем множество кустов, которые присоединены к сети ППД, но по ним нет телеметрии
    # Подтаскиваем по ним какие-то заглушки из режимной закачки
    # Формируем и выгружаем датасет для расчета

    def get_pipelines_and_drop_outliers(self, rs_id, tmp_table_name):
        query = f'''
            with
            q1 as (select * from pipelines p where rs_schema_id = {rs_id}),
            q2 as (select * from q1 where pipe_status='Действующий'),
            q3 as (select * from q2 where node_type_start in (1, 3, 8) and node_type_end in (1, 3, 8))
            select * into temporary {tmp_table_name} from q3
        '''
        self.execute_query(query)


    def make_dataframe_for_calculation(self, timespan, rs_id):
        # Два параметра на вход - какую КНС считаем 1 или 2, и дата-время
        # По первому получаем rs_id, делаем выборку из pipelines. Это уже граф, но неполный
        # Выпиливаем трубы соединяющие куст и скважину
        # Выпиливаем двойные трубы по критерию даты создания или еще какому-то флагу (статусы трубопроводов)
        self.get_pipelines_and_drop_outliers(rs_id, 't_1st_step')




if __name__ == '__main__':
    # conn =  psycopg2.connect(dbname='postgres', user='db_user', password='210106', host='localhost')
    rs_id = [1750000976, 1750000710][0]
    timespan = (datetime.fromisoformat('2019-02-01 01:00'), datetime.fromisoformat('2019-02-01 01:15'))

    pipeline = HE2_ETL('postgresql://db_user:210106@localhost:5432/postgres', dict(host="192.168.4.65", port='80', user="am", password="12345678"))
    pipeline.fill_all_neccesary_tables()
    pipeline.preprocess_wells_Q()
    pipeline.make_dataframe_for_calculation(timespan, rs_id)
    print('Finished')
