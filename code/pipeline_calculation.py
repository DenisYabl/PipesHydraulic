import streamlit as st
import time
import pandas as pd
import logging
import sys
import plotly.graph_objects as go

from DFOperations.calculate_DF import calculate_DF

start_time = time.time()


def main():
    st.title("Расчет трубопроводной сети")
    option = st.selectbox('Выберите участок для расчета',
    ('ДНС 1', 'ДНС 2', 'ДНС 3'))
    selection_dict = {"ДНС 1":  "DNS1_no_wells.csv", "ДНС 2":  "DNS2_no_wells.csv", "ДНС 3":  "DNS3_no_wells.csv"}
    dataset = pd.read_csv("../CommonData/" + selection_dict[option])
    EffD = st.sidebar.slider("Средний эффективный диаметр сети", min_value=0.01, max_value=1.0, value = 0.85, step = 0.01)
    dataset["effectiveD"] = EffD
    baseproperties = [('sat_P_bar', 66.7, "Давление насыщения"), ('plastT_C', 84.0, "Пластовая температура"), ('gasFactor', 39.0, "Газовый фактор"),
                     ('oildensity_kg_m3', 826.0, "Плотность сепарированной нефти"), ('waterdensity_kg_m3', 1015.0, "Плотность пластовой воды"),
                     ('gasdensity_kg_m3', 1.0, "Плотность газа"), ('oilviscosity_Pa_s', 35e-3, "Вязкость нефти"),
                     ('volumeoilcoeff', 1.015, "Объемный коэффициент")]
    property_dict = {}
    for property in baseproperties:
        property_dict.update({property[0]:st.sidebar.text_input(f"{property[2]}", property[1])})
    for key in list(property_dict.keys()):
        dataset[key] = float(property_dict[key])

    inlets_df = dataset[dataset["startIsSource"]]
    inlets_dict = {}
    for index, row in inlets_df.iterrows():
        inlets_dict.update({index:(
        st.sidebar.text_input(f"Давление на {row['node_name_start']}", row["startValue"]),
        st.sidebar.slider(f"Обводненность на {row['node_name_start']}", min_value=0.01, max_value=100.0, value = row["VolumeWater"] if pd.notna(row["VolumeWater"]) else 50.0, step=0.1))})

    for key in list(inlets_dict.keys()):
        dataset.loc[key, "startValue"] = float(inlets_dict[key][0])
        dataset.loc[key, "VolumeWater"] = inlets_dict[key][1]

    with st.spinner("Расчет сети"):
        mdf = calculate_DF(dataset)
    try:

        mdf = mdf[["simple_part_id", 'pipeline_id', 'L', 'node_name_start',
                   'node_name_end','uphillM', 'startP',
                   'endP', 'res_watercut_percent', 'res_liquid_density_kg_m3', 'X_kg_sec', 'velocity_m_sec']]
        mdf = mdf.rename(columns = {"simple_part_id": "ID простого участка", 'pipeline_id': "ID трубопровода", 'L': "Длина участка", 'node_name_start': "Название точки начала",
                   'node_name_end': "Название точки конца",'uphillM': "Перепад высот", 'startP' : "Давление в начале",
                   'endP': "Давление в конце", 'intD': "Внутренний диаметр", 'res_watercut_percent':"Расчетная обводненность",
                                    'res_liquid_density_kg_m3':"Расчетная плотность", 'X_kg_sec':"Массовый расход", 'velocity_m_sec':"Скорость потока" })
        st.write(mdf)
        st.success("Решение найдено")
    except:
        st.error("Решение для данной конфигурации сети не найдено")


if __name__ == "__main__":
    main()