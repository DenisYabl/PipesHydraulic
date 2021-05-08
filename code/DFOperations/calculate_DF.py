from Solver.HE2_Solver import HE2_Solver
from Tools.HE2_schema_maker import make_oilpipe_schema_from_OT_dataset
from Tools.HE2_tools import check_solution
import logging
import sys

def calculate_DF(dataframe):


    G, calc_df = make_oilpipe_schema_from_OT_dataset(dataframe)
    solver = HE2_Solver(G)
    solver.solve(threshold = 0.2)
    modified_dataframe = dataframe.copy()
    if solver.op_result.success == True:
        validity = check_solution(G)
        print(validity)

        for n in G.nodes:
            calc_df.loc[calc_df["node_id_start"] == n, "startP"] = G.nodes[n]["obj"].result["P_bar"]
            calc_df.loc[calc_df["node_id_start"] == n, "startT"] = G.nodes[n]["obj"].result["T_C"]
            calc_df.loc[calc_df["node_id_end"] == n, "endP"] = G.nodes[n]["obj"].result["P_bar"]
            calc_df.loc[calc_df["node_id_end"] == n, "endT"] = G.nodes[n]["obj"].result["T_C"]

    return calc_df