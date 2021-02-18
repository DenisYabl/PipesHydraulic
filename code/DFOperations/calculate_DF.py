from Solver.HE2_Solver import HE2_Solver
from Tools.HE2_schema_maker import make_oilpipe_schema_from_OT_dataset

def calculate_DF(dataframe, logger):
    G = make_oilpipe_schema_from_OT_dataset(dataframe)
    logger.debug("Graph schema is created")
    solver = HE2_Solver(G)
    solver.solve()
    logger.debug("Graph schema is solved")
    modified_dataframe = dataframe.copy()
    for n in G.nodes:
        modified_dataframe.loc[modified_dataframe["node_id_start"] == n, "startP"] = G.nodes[n]["obj"].result["P_bar"]
        modified_dataframe.loc[modified_dataframe["node_id_end"] == n, "endP"] = G.nodes[n]["obj"].result["P_bar"]
    logger.debug("Dataframe is filled with calc results")

    return modified_dataframe