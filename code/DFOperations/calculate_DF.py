from Solver.HE2_Solver import HE2_Solver
from Tools.HE2_schema_maker import make_oilpipe_schema_from_OT_dataset
from Tools.HE2_tools import check_solution, print_solution
from Tools.HE2_Logger import getLogger
logger = getLogger(__name__)

def calculate_DF(dataframe):


    G, calc_df, df_to_graph_edges_mapping = make_oilpipe_schema_from_OT_dataset(dataframe)
    solver = HE2_Solver(G)
    solver.solve(threshold = 1.2)
    modified_dataframe = dataframe.copy()
    if solver.op_result.success == True:
        validity = check_solution(G)
#        print(validity)
#        print_solution(G)

        for n in G.nodes:
            calc_df.loc[calc_df["node_id_start"] == n, "startP"] = G.nodes[n]["obj"].result["P_bar"]
            calc_df.loc[calc_df["node_id_start"] == n, "startT"] = G.nodes[n]["obj"].result["T_C"]
            calc_df.loc[calc_df["node_id_end"] == n, "endP"] = G.nodes[n]["obj"].result["P_bar"]
            calc_df.loc[calc_df["node_id_end"] == n, "endT"] = G.nodes[n]["obj"].result["T_C"]

        calc_df['res_X_kg_sec'] = None
        calc_df['res_watercut_percent'] = None
        calc_df['res_liquid_density_kg_m3'] = None

        for index, row in calc_df.iterrows():
            start = row['node_id_start']
            end = row['node_id_end']
            u, v, k = df_to_graph_edges_mapping[index]
            if start != u or end != v:
                logger.error('broken df_to_graph_edges_mapping, or calc_df.index')
                logger.error(f'{start}, {end}, {u}, {v}, {k}')
                raise IndexError

            obj = G[u][v][k]['obj']
            calc_df.loc[index, 'X_kg_sec'] = obj.result['x']
            calc_df.loc[index, 'res_watercut_percent'] = obj.result['WC']
            calc_df.loc[index, 'res_liquid_density_kg_m3'] = obj.result['liquid_density']
            Q = obj.result['x'] / obj.result['liquid_density']
            Area = 3.1415926 * row['intD'] ** 2 / 4
            V = Q / Area
            calc_df.loc[index, 'velocity_m_sec'] = V

    return calc_df