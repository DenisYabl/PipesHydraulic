import networkx as nx
import pandas as pd
from GraphEdges.HE2_Pipe import HE2_WaterPipe, HE2_OilPipe
from GraphEdges.HE2_Plast import HE2_Plast
from GraphEdges.HE2_WellPump import HE2_WellPump, create_HE2_WellPump_instance_from_dataframe
from GraphNodes import HE2_Vertices as vrtxs
from Fluids.HE2_Fluid import HE2_DummyOil

def make_oilpipe_schema_from_OT_dataset(dataset):
    pump_curves = pd.read_csv("../CommonData/PumpChart.csv")
    outlets = {}
    inlets = {}
    juncs = {}
    calc_df = dataset
    ids_count =  pd.concat((calc_df['node_id_start'], calc_df['node_id_end'])).value_counts()
    ids_count.rename('ids_count')
    calc_df = calc_df.join(ids_count.to_frame(), on='node_id_start', how='left')
    calc_df = calc_df.rename(columns = {0:'start_id_count'})
    calc_df = calc_df.join(ids_count.to_frame(), on='node_id_end', how='left')
    calc_df = calc_df.rename(columns = {0:'end_id_count'})

    calc_df['sourceByCount'] = calc_df['start_id_count'] == 1
    calc_df['outletByCount'] = calc_df['end_id_count'] == 1

    calc_df['sourceMistakes'] = calc_df['sourceByCount'] == calc_df['startIsSource']
    calc_df['outletMistakes'] = calc_df['outletByCount'] == calc_df['endIsOutlet']

    calc_df['sourceValueIsFilled'] = pd.notna(calc_df['startValue'])
    calc_df['outletValueIsFilled'] = pd.notna(calc_df['endValue'])

    calc_df['sourceKindIsFilled'] = pd.notna(calc_df['startKind'])
    calc_df['outletKindIsFilled'] = pd.notna(calc_df['endKind'])

    calc_df['inletBoundaryMistakes'] = True
    calc_df['outletBoundaryMistakes'] = True

    calc_df.loc[calc_df["startIsSource"], 'inletBoundaryMistakes'] = calc_df[calc_df["startIsSource"]]['sourceValueIsFilled'] & calc_df[calc_df["startIsSource"]]['sourceKindIsFilled']

    calc_df.loc[calc_df["endIsOutlet"], 'outletBoundaryMistakes'] = calc_df[calc_df["endIsOutlet"]]['outletValueIsFilled'] & calc_df[calc_df["endIsOutlet"]]['outletKindIsFilled']



    calc_df['sumOfCounts'] = calc_df['start_id_count'] + calc_df['end_id_count']
    calc_df = calc_df[calc_df['sumOfCounts'] > 2]

    mistakes_df = calc_df[(~calc_df['sourceMistakes']) | (~calc_df['outletMistakes']) | (~calc_df['inletBoundaryMistakes']) | (~calc_df['outletBoundaryMistakes'])]

    if not mistakes_df.empty :
        print(f"Following nodes: {mistakes_df[~mistakes_df['sourceMistakes']]['node_id_start'].values} should be sources")
        print(f"Following nodes: {mistakes_df[~mistakes_df['outletMistakes']]['node_id_start'].values} should be outlets")
        print(f"Start kInd and value for following nodes: {mistakes_df[~mistakes_df['inletBoundaryMistakes']]['node_id_start'].values} should be filled")
        print(f"End kInd and value for following nodes: {mistakes_df[~mistakes_df['outletBoundaryMistakes']]['node_id_end'].values} should be filled")
        assert False



    inlets_df = calc_df[dataset["startIsSource"]]
    outletsdf = calc_df[dataset["endIsOutlet"]]
    juncs_df = pd.concat((dataset["node_id_start"], dataset["node_id_end"])).unique()
    for index, row in inlets_df.iterrows():
        volumewater = row["VolumeWater"] if pd.notna(row["VolumeWater"]) else 50
        inlets.update({row["node_id_start"]:vrtxs.HE2_Source_Vertex(row["startKind"], row["startValue"], HE2_DummyOil(volumewater), row["startT"])})
    for index, row in outletsdf.iterrows():
        outlets.update({row["node_id_end"]: vrtxs.HE2_Boundary_Vertex(row["endKind"], row["endValue"])})
    for id in juncs_df:
        if (id not in list(inlets.keys()) + list(outlets.keys())):
            juncs.update({id:vrtxs.HE2_ABC_GraphVertex()})

    G = nx.DiGraph()  # Di = directed

    for k, v in {**inlets, **outlets, **juncs}.items():
        G.add_node(k, obj=v)

    for index, row in calc_df.iterrows():
        start = row["node_id_start"]
        end = row["node_id_end"]
        junctype = row["juncType"]
        volumewater = row["VolumeWater"] if pd.notna(row["VolumeWater"]) else 50
        if junctype == "pipe":
            L = row["L"]
            uphill = row["uphillM"]
            diam_coef = row["effectiveD"]
            D = row["intD"]
            roughness = row["roughness" ]
            G.add_edge(start, end, obj=HE2_OilPipe([L], [uphill], [D * diam_coef], [roughness], [HE2_DummyOil(volumewater)]))
        elif junctype == "plast":
            productivity = row["productivity" ]
            G.add_edge(start, end, obj=HE2_Plast(productivity=productivity, fluid=HE2_DummyOil(volumewater)))
        elif junctype == "wellpump":
            model = row["model"]
            frequency = row["frequency"]
            G.add_edge(start, end, obj = create_HE2_WellPump_instance_from_dataframe(full_HPX=pump_curves, model=model, fluid=HE2_DummyOil(volumewater),
                       frequency=frequency))
    return G


def make_schema_from_OISPipe_dataframes(df_pipes, df_boundaries):
    df = df_pipes[["node_id_start", "node_id_end"]]
    df.columns = ["source", "target"]
    G = nx.from_pandas_edgelist(df, create_using=nx.DiGraph)
    edge_list = list(G.edges())
    edge_set = set(edge_list)
    rev_set = set([(v, u) for u, v in edge_set])
    fin_edge_set = edge_set - rev_set
    G = nx.DiGraph(fin_edge_set)

    cmpnts = nx.algorithms.components.number_connected_components(nx.Graph(fin_edge_set))
    if cmpnts != 1:
        print('Not single component graph!')
        assert False

    pipes = dict()
    for u, v in G.edges():
        df = df_pipes
        df = df[(df.node_id_start == u) & (df.node_id_end == v)]
        d = df.iloc[0].to_dict()

        Ls = [d['L']]
        Hs = [d['altitude_end'] - d['altitude_start']]
        Ds = [d['D'] - 2 * d['S']]
        Rs = [1e-5]

        pipe = HE2_WaterPipe(Ls, Hs, Ds, Rs)
        pipes[(u, v)] = pipe
    nx.set_edge_attributes(G, name='obj', values=pipes)


    df = df_boundaries[['Q', 'P']].fillna(-1e9)
    df_boundaries['value'] = df.max(axis=1)

    nodes = dict()
    id_list = list(df_boundaries.id.values)
    for n in G.nodes():
        df = df_boundaries
        if n in id_list:
            df = df[df.id == n]
            d = df.iloc[0].to_dict()
        else:
            obj = vrtxs.HE2_ABC_GraphVertex()
            nodes[n] = obj
            continue
        # if d['kind']=='P':
        #     print(d)

        if d['is_source']:
            obj = vrtxs.HE2_Source_Vertex(d['kind'], d['value'], 'water', 20)
        elif d['kind']=='Q' and ((d['Q'] is None) or d['Q']==0):
            obj = vrtxs.HE2_ABC_GraphVertex()
        else:
            obj = vrtxs.HE2_Boundary_Vertex(d['kind'], d['value'])
        nodes[n] = obj

    nx.set_node_attributes(G, name='obj', values=nodes)


    for n in G.nodes():
        o = G.nodes[n]['obj']
        assert o is not None

    for u, v in G.edges():
        o = G[u][v]['obj']
        assert o is not None

    return G


def make_multigraph_schema_from_OISPipe_dataframes(df_pipes, df_boundaries):
    df = df_pipes[["node_id_start", "node_id_end"]]
    df["idx_for_result"] = df.index
    df.columns = ["source", "target", "idx_for_result"]
    G = nx.from_pandas_edgelist(df, create_using=nx.MultiDiGraph, edge_attr="idx_for_result")

    cmpnts = nx.algorithms.components.number_connected_components(nx.Graph(G))
    if cmpnts != 1:
        print('Not single component graph!')
        assert False

    pipes = dict()
    for u, v, k in G.edges:
        df = df_pipes
        df = df[(df.node_id_start == u) & (df.node_id_end == v)]
        d = df.iloc[k].to_dict()

        Ls = [d['L']]
        Hs = [d['altitude_end'] - d['altitude_start']]
        Ds = [d['D'] - 2 * d['S']]
        Rs = [1e-5]

        pipe = HE2_WaterPipe(Ls, Hs, Ds, Rs)
        pipes[(u, v, k)] = pipe
    nx.set_edge_attributes(G, name='obj', values=pipes)

    mask = df_boundaries.kind == 'Q'
    df_boundaries['value'] = -1e9
    df_boundaries.loc[mask, 'value'] = df_boundaries.loc[mask, 'Q']
    df_boundaries.loc[~mask, 'value'] = df_boundaries.loc[~mask, 'P']

    nodes = dict()
    id_list = list(df_boundaries.id.values)
    for n in G.nodes():
        df = df_boundaries
        if n in id_list:
            df = df[df.id == n]
            d = df.iloc[0].to_dict()
        else:
            obj = vrtxs.HE2_ABC_GraphVertex()
            nodes[n] = obj
            continue
        # if d['kind']=='P':
        #     print(d)

        if d['is_source']:
            obj = vrtxs.HE2_Source_Vertex(d['kind'], d['value'], 'water', 20)
        elif d['kind']=='Q' and ((d['Q'] is None) or d['Q']==0):
            obj = vrtxs.HE2_ABC_GraphVertex()
        else:
            obj = vrtxs.HE2_Boundary_Vertex(d['kind'], d['value'])
        nodes[n] = obj

    nx.set_node_attributes(G, name='obj', values=nodes)


    for n in G.nodes():
        o = G.nodes[n]['obj']
        assert o is not None

    for u, v, k in G.edges:
        o = G[u][v][k]['obj']
        assert o is not None

    return G
