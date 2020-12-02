import pandas as pd
import HE2_schema_maker as sm
import HE2_tools as tools
from HE2_Solver import HE2_Solver
import numpy as np

def process_error():
    pass

def split_input_df_to_pipes_and_boundaries(df):
    bnd_cols = ['kind', 'Q', 'is_source', 'P']
    start_cols = ['node_id_start'] + ['start_' + col for col in bnd_cols]
    end_cols = ['node_id_end'] + ['end_' + col for col in bnd_cols]
    df_bnd_start = input_df[start_cols]
    df_bnd_start.columns = ['id'] + bnd_cols
    df_bnd_end = input_df[end_cols]
    df_bnd_end.columns = ['id'] + bnd_cols
    df_bnds = pd.concat([df_bnd_start, df_bnd_end])
    df_bnds = df_bnds[~df_bnds.kind.isna()]
    to_drop = ['start_' + col for col in bnd_cols] + ['end_' + col for col in bnd_cols]
    df_pipes = input_df.drop(columns=to_drop)
    return df_pipes, df_bnds


def do_something_right_now_good_name_will_come_later(input_df):
    df_pipes, df_bnds = split_input_df_to_pipes_and_boundaries(input_df)

    G = sm.make_multigraph_schema_from_OISPipe_dataframes(df_pipes, df_bnds)

    solver = HE2_Solver(G)
    solver.solve()
    op_result = solver.op_result
    if op_result.fun > 1e-3:
        process_error()
        return None

    resd1, resd2 = tools.check_solution(G)
    if resd1 + resd2 > 1e-3:
        process_error()
        return None

    result_cols = ['mass_flow', 'flow_direction', 'result_start_P', 'result_start_T', 'result_end_P', 'result_end_T']
    result_df = input_df.copy()
    for col in result_cols:
        result_df[col] = 0

    for u, v, k in G.edges:
        idx = G[u][v][k]["idx_for_result"]
        obj = G.nodes[u]['obj']
        P1 = obj.result['P_bar']
        T1 = obj.result['T_C']
        obj = G.nodes[v]['obj']
        P2 = obj.result['P_bar']
        T2 = obj.result['T_C']
        obj = G[u][v][k]['obj']
        x = obj.result['x']
        dir = np.sign(x)
        result_df.loc[idx, result_cols] = [abs(x), dir, P1, T1, P2, T2]
    return result_df

input_df = pd.read_csv('..\\data\\ends_202011301340.csv')
result_df = do_something_right_now_good_name_will_come_later(input_df)
result_df.to_csv('..\\data\\rez1.csv')