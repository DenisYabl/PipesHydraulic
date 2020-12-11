import pandas as pd
import HE2_schema_maker as sm
import HE2_tools as tools
from HE2_Solver import HE2_Solver
import numpy as np
from matplotlib import pyplot as plt


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


def solve_and_put_results_to_dataframe(input_df):
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

    # tools.draw_solution(G, shifts=None, p_nodes=[], sources=[], sinks=[], juncs=[])

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

def do_upcase_columns_adhoc(columns):
    rez = []
    for col in columns:
        c = col
        if len(col) == 1:
            c = col.upper()
        if col[-1] in ('p', 'q'):
            c = col[:-1] + col[-1].upper()
        rez += [c]
    return rez

def split_result_df(result_df):
    edges_cols = ['rs_schema_id', 'creation_date', 'L', 'D', 'S', 'mass_flow', 'flow_direction']
    edges_df = result_df[edges_cols]
    nodes_cols = ['node_id', 'node_name', 'altitude', 'node_type', 'x', 'y', 'kind', 'Q', 'P', 'result_P', 'result_T']

    start_cols = ['node_id_start', 'node_name_start', 'altitude_start', 'node_type_start', 'x_start', 'y_start', 'start_kind', 'start_Q', 'start_P', 'result_start_P', 'result_start_T']
    end_cols = ['node_id_end', 'node_name_end', 'altitude_end', 'node_type_end', 'x_end', 'y_end', 'end_kind', 'end_Q', 'end_P', 'result_end_P', 'result_end_T']

    df1 = result_df[start_cols]
    df1.columns = nodes_cols
    df2 = result_df[end_cols]
    df2.columns = nodes_cols
    nodes_df = pd.concat([df1, df2]).drop_duplicates()
    return edges_df, nodes_df

def transform_measure_units(df):
    df['D'] = df.D / 1000
    df['S'] = df.S / 1000
    df['start_Q'] = df.start_Q * 1000 / 86400
    df['end_Q'] = df.end_Q * 1000 / 86400
    return df

def do_predict(input_df):
    df = input_df
    df.columns = do_upcase_columns_adhoc(df.columns)
    df = transform_measure_units(df)
    result_df = solve_and_put_results_to_dataframe(df)
    return result_df


if __name__ == '__main__':
    input_df = pd.read_csv('..\\data\\q4_202012041333.csv')
    result_df = do_predict(input_df)
    result_df.to_csv('..\\data\\rez1.csv')
