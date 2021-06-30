import streamlit as st
import time
import pandas as pd
import logging
import sys
import plotly.graph_objects as go



from DFOperations.calculate_DF import calculate_DF

start_time = time.time()

def main():
    st.title("Расчет давлений")
    dataset = pd.read_csv("../CommonData/1well.csv")
    pump_curves  = pd.read_csv("../CommonData/PumpChart.csv").set_index("pumpModel")

    
    st.sidebar.write("Параметры пласта")
    P_plast = st.sidebar.slider("Пластовое давление", min_value=0, max_value=350)
    dataset["startValue"] = P_plast
    
    P_ust = st.sidebar.slider("Устьевое давление", min_value=0, max_value=50)
    dataset["endValue"] = P_ust
    
    productivity = st.sidebar.slider("Коэффициент продуктивности", min_value=0.1, max_value=3.0)
    dataset["productivity"] = productivity
    
    perforation = st.sidebar.slider("Глубина перфорации", min_value=2500, max_value=4000, value=3759)
    dataset["perforation"] = perforation
    
    st.sidebar.write("Параметры насоса")
    pumpDepth = st.sidebar.slider("Глубина спуска насоса", min_value=1500, max_value=3500, value=2100)
    dataset["pumpDepth"] = pumpDepth
    
    K_pump = st.sidebar.slider("Коэффициент заполнения насоса", min_value=0.0, max_value=1.0, value=0.7)
    dataset["K_pump"] = K_pump
    
    frequency = st.sidebar.slider("Частота вращения насоса", min_value=0, max_value=70, value=50)
    dataset["frequency"] = frequency
    
    st.sidebar.write("Параметры жидкости")
    sat_P_bar = st.sidebar.slider("Давление насыщения", min_value=0, max_value=100, value=63)
    dataset["sat_P_bar"] = sat_P_bar

    plastT_C = st.sidebar.slider("Пластовая температура", min_value=0, max_value=100, value=84)
    dataset["plastT_C"] = plastT_C
    
    gasFactor = st.sidebar.slider("Газовый фактор", min_value=0, max_value=100, value=49)
    dataset["gasFactor"] = gasFactor
    
    oildensity_kg_m3 = st.sidebar.slider("Плотность сепарированной нефти", min_value=700, max_value=1200, value=826)
    dataset["oildensity_kg_m3"] = oildensity_kg_m3
    
    waterdensity_kg_m3 = st.sidebar.slider("Плотность пластовой воды", min_value=900, max_value=1100, value=1015)
    dataset["waterdensity_kg_m3"] = waterdensity_kg_m3
    
    gasdensity_kg_m3 = st.sidebar.slider("Плотность попутного газа", min_value=0.1, max_value=2.1, value=1.0)
    dataset["gasdensity_kg_m3"] = gasdensity_kg_m3
    
    oilviscosity_Pa_s = st.sidebar.slider("Вязкость сепарированной нефти", min_value=0.0, max_value=100e-3, value=35e-3)
    dataset["oilviscosity_Pa_s"] = oilviscosity_Pa_s
    
    volumeoilcoeff = st.sidebar.slider("Объемный коэффициент нефти", min_value=0.0, max_value=2.0, value=1.015)
    dataset["volumeoilcoeff"] = volumeoilcoeff
  

    st.write(dataset)
    pump_curve = pump_curves.loc[dataset["model"]].sort_values(by="debit")
    st.write(pump_curve)
    mdf = calculate_DF(dataset)
    Q = mdf.iloc[-1]
    st.write(Q)
    Q =  Q["X_kg_sec"] * 86400/Q["res_liquid_density_kg_m3"]
    st.write(f"Текущий дебит - {Q} м^3/сутки")
    pressure=mdf["startP"][1:]
    st.write(pressure)

 
    fig = go.Figure(
	    data=[go.Scatter(x=[perforation, pumpDepth, pumpDepth-100, 0], y=pressure)],
	    layout=go.Layout(
		title=go.layout.Title(text="Градиент давления")
	    ))
    st.write(fig)
    
    fig = go.Figure()
    fig.add_trace(go.Scatter(x=pump_curve["debit"], y=pump_curve["pressure"]))
    fig.add_trace(go.Scatter(x=pump_curve["debit"], y=pump_curve["eff"]))
    fig.add_vline(x=Q, line_width=3, line_dash="dash", line_color="green")
    st.write(fig)
    
    st.write("Моделирование периодического режима")
    well_973 = pd.read_csv("~/well_973.csv", parse_dates=["index"])
    
    def calc_q(dataset):
        mdf = calculate_DF(dataset)
        Q = mdf.iloc[-1]
        Q =  Q["X_kg_sec"] * 86400/Q["res_liquid_density_kg_m3"]
        return Q
    Q_list = []
    if st.button("Запустить периодический режим"):
        for _, state in well_973.iterrows():
            temp_dataset = dataset.copy()
        
            if state["status"] == 1.0:
              q = calc_q(dataset)
            else:
                temp_dataset["frequency"] = 1.0
                q = calc_q(temp_dataset)
                if q < 0: 
                   q = 0
            Q_list.append(q)
        well_973["Q"] = Q_list
        st.write(well_973)
        fig = go.Figure()
        fig.add_trace(go.Scatter(x=well_973["index"], y=well_973["status"], line_shape="hv"))
        fig.add_trace(go.Scatter(x=well_973["index"], y=well_973["Q"], line_shape="hv"))
        st.write(fig)


if __name__=="__main__":
    main()
