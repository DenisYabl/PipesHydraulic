from Solver.HE2_Solver import HE2_Solver
from Tools.HE2_schema_maker import make_oilpipe_schema_from_OT_dataset
import logging
import sys

def calculate_DF(dataframe, logger = None):
    if logger is None:
        logger = logging.getLogger('Python debug')
        formatter = logging.Formatter('%(asctime)s,%(msecs)d %(name)s %(levelname)s %(message)s')
        filehandler = logging.FileHandler(filename='run.log', mode='w')
        filehandler.setFormatter(formatter)
        streamhandler = logging.StreamHandler(sys.stderr)
        streamhandler.setFormatter(formatter)

        logger.addHandler(filehandler)
        logger.addHandler(streamhandler)

    G = make_oilpipe_schema_from_OT_dataset(dataframe)
    logger.debug("Graph schema is created".encode())
    solver = HE2_Solver(G)
    solver.solve()
    logger.debug("Graph schema is solved".encode())
    modified_dataframe = dataframe.copy()
    if solver.op_result.success == True:
        for n in G.nodes:
            modified_dataframe.loc[modified_dataframe["node_id_start"] == n, "startP"] = G.nodes[n]["obj"].result["P_bar"]
            modified_dataframe.loc[modified_dataframe["node_id_start"] == n, "startT"] = G.nodes[n]["obj"].result["T_C"]
            modified_dataframe.loc[modified_dataframe["node_id_end"] == n, "endP"] = G.nodes[n]["obj"].result["P_bar"]
            modified_dataframe.loc[modified_dataframe["node_id_end"] == n, "endT"] = G.nodes[n]["obj"].result["T_C"]
        logger.debug("Dataframe is filled with calc results".encode())
    else:
        logger.debug("Dataframe is not filled, solver failed to find solution".encode())
    return modified_dataframe