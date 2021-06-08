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
    P_plast = st.sidebar.slider("Пластовое давление", min_value=0, max_value=350)
    dataset["startValue"] = P_plast
    P_ust = st.sidebar.slider("Устьевое давление", min_value=0, max_value=50)
    dataset["endValue"] = P_ust
    st.write(dataset)
    mdf = calculate_DF(dataset)
    st.write(mdf)
    pressure=mdf["startP"][1:]
    fig = go.Figure(
	    data=[go.Scatter(x=[2803, 2103, 1903, 0], y=pressure)],
	    layout=go.Layout(
		title=go.layout.Title(text="Градиент давления")
	    ))
    st.write(fig)
    
    Q = mdf.iloc[-1]
    st.write(Q)
    Q =  Q["X_kg_sec"] * 86400/Q["res_liquid_density_kg_m3"]
    st.write(Q)


if __name__=="__main__":
    main()
