import time
from DFOperations.calculate_DF import calculate_DF
import pandas as pd
import logging
import sys
from Tools.HE2_schema_maker import make_oilpipe_schema_from_OT_dataset

"""
За полный пересчет датафрейма DF => DF отвечает функция DFOperations.calculate_DF.calculate_DF
Входной датафрейм должен соответствовать testDf по следующим полям:

node_id_start, node_id_end - для создания узлов графа и дуг между ними 

juncType - тип соединения узлов, pipe - нефтяная труба, wellpump - насос скважины, plast - призабойная часть пласта. 
Все поля, необязательные для данного juncType могут быть незаполнены

L - длина, обязательна для труб

altitude_start, altitude_end, uphillM - высотный перепад, обязателен для труб. Непосредственно расчет использует поле uphillM = altitude_end - altitude_start

D, S, intD - внутренний диаметр, обязателен для труб. Непосредственно расчет использует поле intD = (D - 2 * S) / 1000

effectiveD - эффективный диаметр трубы, учитывающий отложения, задается от 0 до 1, обязателен для труб

roughness - шероховатость, обязательна для труб

productivity - коэффициент продуктивности, обязателен для пластов

model, frequency - модель и частота работы насоса, обязательна для насосов. Модель должна входить в актуальнын НРХ, расположенные в "/CommonData/PumpChart.csv" 

startIsSource, startKind, startValue - признак источника (True или любое другое значение), тип источника ('P' или 'Q') и значение. 
Дополнительно по node_id проверяется уникальность узла, все узлы с уникальным node_id_start должны быть источниками, и у всех узлов с startIsSource == True 
должны быть заполнены поля  startKind и startValue, а также startT

endIsOutlet, endKind, endValue - выходы из системы, по аналогии с источниками

startP, startT, endP, endT - расчетные поля, заполняемые по мере расчета

Функция дополнительно отфильтровывает дуги, не связанные с основным графом (дуги с уникальным началом и концом, не соединенные с другими узлами)
"""

dataset = pd.read_csv("../CommonData/DNS2_with_wells.csv")
start_time = time.time()
mdf = calculate_DF(dataset)
mdf = mdf[['juncType', 'pipeline_purpose_id', 'simple_part_id', 'part_id', 'rs_schema_id', 'schema_id',
           'pipeline_id', 'node_id_end', 'node_id_start', 'L', 'simple_part_creation_date', 'node_name_start',
           'altitude_start',
           'node_type_start', 'node_name_end', 'altitude_end', 'node_type_end', 'D', 'S', 'thread_number', 'uphillM',
           'startIsSource', 'VolumeWater',
           'startKind', 'startValue', 'endIsOutlet', 'endKind', 'endValue', 'startP', 'startT', 'endP', 'endT',
           'effectiveD', 'intD', 'roughness',
           'productivity', 'model', 'frequency', 'perforation', 'pumpDepth', 'wellNum', 'padNum']]
mdf = mdf.fillna(0)
print("--- %s seconds ---" % (time.time() - start_time))
mdf.to_csv("../CommonData/mdf.csv")


"""
Также возможно создание расчетного графа из датафрейма без его расчета с помощью Tools.HE2_schema_maker.make_oilpipe_schema_from_OT_dataset
Требования к датафрейму такие же, как для calculate_DF
"""

# G, df, df_to_graph_edges_mapping = make_oilpipe_schema_from_OT_dataset(dataset)
pass