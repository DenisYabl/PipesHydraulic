import streamlit as st
import streamlit.components.v1 as components
import time
import pandas as pd
import logging
import sys
import plotly.graph_objects as go

from DFOperations.calculate_DF import calculate_DF
from Visualisation.pyvis_v0 import draw_result_graph

start_time = time.time()


def main():
    st.set_page_config(layout="wide")
    st.title("Расчет трубопроводной сети")
    option = st.selectbox('Выберите участок для расчета',
    ('ДНС 1', 'ДНС 2', 'ДНС 3'))
    selection_dict = {"ДНС 1":  "DNS1_real_values.csv", "ДНС 2":  "DNS2_real_values.csv", "ДНС 3":  "DNS3_real_values.csv"}
    dataset = pd.read_csv("../CommonData/" + selection_dict[option])
    use_coord = st.sidebar.checkbox(label='Отрисовка сети по координатам')
    coordinate_scaling = float(st.sidebar.text_input(f"Масштаб координат", '3'))
    EffD = st.sidebar.slider("Средний эффективный диаметр сети", min_value=0.01, max_value=1.0, value = 0.85, step = 0.01)
    dataset["effectiveD"] = EffD
    with st.form("Свойства жидкости"):
        cols = st.beta_columns(8)
        baseproperties = [('sat_P_bar', 66.7, "Давление насыщения"), ('plastT_C', 84.0, "Пластовая температура"), ('gasFactor', 39.0, "Газовый фактор"),
                         ('oildensity_kg_m3', 826.0, "Плотность сепарированной нефти"), ('waterdensity_kg_m3', 1015.0, "Плотность пластовой воды"),
                         ('gasdensity_kg_m3', 1.0, "Плотность газа"), ('oilviscosity_Pa_s', 35e-3, "Вязкость нефти"),
                         ('volumeoilcoeff', 1.015, "Объемный коэффициент")]
        property_dict = {}
        for i, col in enumerate(cols):
            property = baseproperties[i]
            property_dict.update({property[0]:col.text_input(f"{property[2]}", property[1])})
        submitted = st.form_submit_button("Submit")
        if submitted:
            for key in list(property_dict.keys()):
                dataset[key] = float(property_dict[key])

    inlets_df = dataset[dataset["startIsSource"]]
    inlets_dict = {}
    for index, row in inlets_df.sort_values("node_name_start").iterrows():
        inlets_dict.update({index:(
        st.sidebar.text_input(f"Давление на {row['node_name_start']}", row["startValue"]),
        st.sidebar.slider(f"Обводненность на {row['node_name_start']}", min_value=0.01, max_value=100.0, value = row["VolumeWater"] if pd.notna(row["VolumeWater"]) else 50.0, step=0.1))})

    for key in list(inlets_dict.keys()):
        dataset.loc[key, "startValue"] = float(inlets_dict[key][0])
        dataset.loc[key, "VolumeWater"] = inlets_dict[key][1]

    outlets_df = dataset[dataset["endIsOutlet"]].copy().drop_duplicates(subset='node_id_end')
    outlets_dict = {}
    for index, row in outlets_df.sort_values("node_name_start").iterrows():
        outlets_dict.update({index:
        st.sidebar.text_input(f"P в узле {row['node_name_end']}", row["endValue"])})

    for key in list(outlets_dict.keys()):
        dataset.loc[key, "endValue"] = float(outlets_dict[key])

    with st.spinner("Расчет сети"):
        mdf, G = calculate_DF(dataset, return_graph=True)
        pyvis_graph = draw_result_graph(mdf, G, use_coordinates=use_coord, coordinate_scaling=coordinate_scaling)
    try:
        mdf = mdf[[ 'node_name_start',
                   'node_name_end','L','uphillM', 'startP',
                   'endP', 'res_watercut_percent', 'res_liquid_density_kg_m3', 'X_kg_sec', 'velocity_m_sec']]

        mdf[['L', 'uphillM', 'startP','endP', 'res_watercut_percent', 'res_liquid_density_kg_m3', 'X_kg_sec','velocity_m_sec']] = \
            mdf[['L', 'uphillM', 'startP','endP', 'res_watercut_percent', 'res_liquid_density_kg_m3', 'X_kg_sec','velocity_m_sec']].astype(float).round(1)



        mdf = mdf.rename(columns = {'L': "Длина участка", 'node_name_start': "Название точки начала",
                   'node_name_end': "Название точки конца",'uphillM': "Перепад высот", 'startP' : "Давление в начале",
                   'endP': "Давление в конце", 'intD': "Внутренний диаметр", 'res_watercut_percent':"Расчетная обводненность",
                                    'res_liquid_density_kg_m3':"Расчетная плотность", 'X_kg_sec':"Массовый расход", 'velocity_m_sec':"Скорость потока" })

        st.write(mdf)
        st.success("Решение найдено")
    except:
        st.error("Решение для данной конфигурации сети не найдено")

    pyvis_graph.save_graph('temp.html')
    HtmlFile = open("temp.html", 'r', encoding='utf-8')
    source_code = HtmlFile.read()
    components.html(source_code, height=900, width=1200)

if __name__ == "__main__":
    main()