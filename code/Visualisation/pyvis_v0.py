import time
from DFOperations.calculate_DF import calculate_DF
import pandas as pd
import logging
import sys

from GraphNodes.HE2_Vertices import HE2_Boundary_Vertex, HE2_Source_Vertex
from Tools.HE2_schema_maker import make_oilpipe_schema_from_OT_dataset
import numpy as np
from pyvis.network import Network

def draw_result_graph(dataset, G, use_coordinates = False, coordinate_scaling = 3):
    id_dict = pd.DataFrame()
    id_dict["ids"] = pd.concat((dataset["node_id_start"], dataset['node_id_end']))
    id_dict["names"] = pd.concat((dataset["node_name_start"], dataset['node_name_end']))

    if use_coordinates:
        id_dict["Xs"] = pd.concat((dataset["X_start"], dataset['X_end']))
        id_dict["Ys"] = pd.concat((dataset["Y_start"], dataset['Y_end']))
        id_dict["Xs"] = id_dict["Xs"] / coordinate_scaling
        id_dict["Ys"] = id_dict["Ys"] / coordinate_scaling

    id_dict = id_dict.drop_duplicates(subset="ids")
    dataset.to_csv("../CommonData/dataset.csv")

    nt = Network('768px', '1024px', directed=True)
    nt.toggle_physics(False)
    base_size = 15
    for n in G.nodes:
        propetries = id_dict.loc[id_dict["ids"] == n]
        obj = G.nodes[n]["obj"]
        if isinstance(obj, HE2_Source_Vertex):
            color = 'green'
            size = base_size
        elif isinstance(obj, HE2_Boundary_Vertex):
            color = 'red'
            size = 2 * base_size
        else:
            color = 'blue'
            size = base_size * 0.75
        if use_coordinates:
            nt.add_node(n, label=propetries['names'].iloc[0] + f'\n P в узле {np.round(obj.result["P_bar"], decimals=2)}',
                        size=size, title=propetries['names'].iloc[0] + f'<br> P в узле {np.round(G.nodes[n]["obj"].result["P_bar"], decimals=2)}',
                        x=propetries['Xs'].iloc[0], y=propetries['Ys'].iloc[0], color = color)
        else:
            nt.add_node(n, label=propetries['names'].iloc[0] + f'\n P в узле {np.round(obj.result["P_bar"], decimals=2)}',
                        size=size, title=propetries['names'].iloc[0] + f'<br> P в узле {np.round(G.nodes[n]["obj"].result["P_bar"], decimals=2)}',
                        color=color)
    for i, row in dataset.iterrows():
        start = row["node_id_start"] if row["X_kg_sec"] > 0 else row["node_id_end"]
        end = row["node_id_end"] if row["X_kg_sec"] > 0 else row["node_id_start"]
        volumewater = np.round(row["res_watercut_percent"], decimals=2)
        title = f'Расход {np.round(row["X_kg_sec"] * 86400 / row["res_liquid_density_kg_m3"], decimals=2)} м3/сутки ' \
                f'<br> Обв. {volumewater} %'
        nt.add_edge(start, to=end, title=title,
                    width=base_size / 2)  # value = abs(np.round(row["X_kg_sec"] *86400 / row["res_liquid_density_kg_m3"], decimals=2)))
    return nt