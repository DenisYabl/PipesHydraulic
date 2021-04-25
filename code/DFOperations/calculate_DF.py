from Solver.HE2_Solver import HE2_Solver
from Tools.HE2_schema_maker import make_oilpipe_schema_from_OT_dataset
import logging
import sys

def calculate_DF(dataframe):


    G = make_oilpipe_schema_from_OT_dataset(dataframe)
    solver = HE2_Solver(G)
    solver.solve()
    modified_dataframe = dataframe.copy()
    if solver.op_result.success == True:
        for n in G.nodes:
            modified_dataframe.loc[modified_dataframe["node_id_start"] == n, "startP"] = G.nodes[n]["obj"].result["P_bar"]
            modified_dataframe.loc[modified_dataframe["node_id_start"] == n, "startT"] = G.nodes[n]["obj"].result["T_C"]
            modified_dataframe.loc[modified_dataframe["node_id_end"] == n, "endP"] = G.nodes[n]["obj"].result["P_bar"]
            modified_dataframe.loc[modified_dataframe["node_id_end"] == n, "endT"] = G.nodes[n]["obj"].result["T_C"]

    return modified_dataframe