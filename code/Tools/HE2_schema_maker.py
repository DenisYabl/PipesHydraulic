import networkx as nx
import pandas as pd
from GraphEdges.HE2_Pipe import HE2_WaterPipe, HE2_OilPipe
from GraphEdges.HE2_Plast import HE2_Plast
from GraphEdges.HE2_WellPump import HE2_WellPump, create_HE2_WellPump_instance_from_dataframe
from GraphNodes import HE2_Vertices as vrtxs
from Fluids.HE2_Fluid import gimme_dummy_BlackOil
from Tools.HE2_Logger import getLogger

logger = getLogger(__name__)

pump_curves = None
inclination = None
HKT = None

def make_oilpipe_schema_from_OT_dataset(dataset, folder="../CommonData/", calc_df = None, ignore_Watercut=False):
    global pump_curves
    if pump_curves is None:
        pump_curves = pd.read_csv(folder + "PumpChart.csv")

    dataset_cols = set(dataset.columns)
    well_cols = ['perforation', 'frequency', 'wellNum', 'padNum', 'pumpDepth', 'productivity', 'model']
    for col in well_cols:
        if not col in dataset_cols:
            dataset[col] = None

    if calc_df is None:
        calc_df = make_calc_df(dataset, folder)

    outlets = {}
    inlets = {}
    juncs = {}

    inlets_df = calc_df[calc_df["startIsSource"]]
    outletsdf = calc_df[calc_df["endIsOutlet"]]
    juncs_df = pd.concat((calc_df["node_id_start"], calc_df["node_id_end"])).unique()
    for index, row in inlets_df.iterrows():
        if ignore_Watercut:
            volumewater = 50
        elif pd.notna(row["VolumeWater"]):
            volumewater = row["VolumeWater"]
        else:
            logger.warning(f'Watercut should be known for source nodes: {row["node_id_start"]}')
            volumewater = 50

        inlets.update({row["node_id_start"]:vrtxs.HE2_Source_Vertex(row["startKind"], row["startValue"],
                                                                    gimme_dummy_BlackOil(VolumeWater=volumewater), row["startT"])})
    for index, row in outletsdf.iterrows():
        outlets.update({row["node_id_end"]: vrtxs.HE2_Boundary_Vertex(row["endKind"], row["endValue"])})
    for id in juncs_df:
        if (id not in list(inlets.keys()) + list(outlets.keys())):
            juncs.update({id:vrtxs.HE2_ABC_GraphVertex()})

    G = nx.MultiDiGraph()  # Di = directed
    df_to_graph_edges_mapping = dict()

    for k, v in {**inlets, **outlets, **juncs}.items():
        G.add_node(k, obj=v)

    for index, row in calc_df.iterrows():
        start = row["node_id_start"]
        end = row["node_id_end"]
        junctype = row["juncType"]
        VolumeWater = row["VolumeWater"] if pd.notna(row["VolumeWater"]) else 50
        fluid = gimme_dummy_BlackOil(VolumeWater=VolumeWater)
        if junctype == "pipe":
            L = row["L"]
            uphill = row["uphillM"]
            diam_coef = row["effectiveD"]
            D = row["intD"]
            roughness = row["roughness"]
            obj = HE2_OilPipe([L], [uphill], [D * diam_coef], [roughness], fluid)
        elif junctype == "plast":
            productivity = row["productivity" ]
            obj = HE2_Plast(productivity=productivity, fluid=fluid)
        elif junctype == "wellpump":
            model = row["model"]
            frequency = row["frequency"]
            pump = create_HE2_WellPump_instance_from_dataframe(full_HPX=pump_curves, model=model, fluid=fluid, frequency=frequency)
            if 'K_pump' in calc_df.columns:
                K_pump = row['K_pump']
                if K_pump > 0 and K_pump < 100500:
                    pump.change_stages_ratio(K_pump)
            obj = pump
        else:
            logger.warning(f'unknown type of graph edge in dataset. start, end id is {start} {end}. Cannot add this edge to graph, skip')
            continue

        k = G.add_edge(start, end, obj=obj) # edge index in MultiDiGraph. k in G[u][v][k] index
        df_to_graph_edges_mapping[index] = (start, end, k)

    cmpnts = nx.algorithms.components.number_weakly_connected_components(G)
    if cmpnts != 1:
        logger.error(f'Not single componented graph!')
        raise ValueError

    return G, calc_df, df_to_graph_edges_mapping


def make_calc_df(dataset, folder):
    global inclination, HKT
    if inclination is None:
        inclination = pd.read_parquet(folder + "inclination")
    if HKT is None:
        HKT = pd.read_parquet(folder + "HKT")
    wells_df = dataset[dataset["juncType"] == 'oilwell']
    dataset = dataset[dataset["juncType"] != 'oilwell']
    dataset[['node_id_start', 'node_id_end']] = dataset[['node_id_start', 'node_id_end']].astype(int).astype(str)
    dict_list = []
    tempdf = populate_wells_df(HKT, dict_list, inclination, wells_df)
    dataset = dataset.append(tempdf)
    calc_df = dataset
    calc_df["startIsSource"] = calc_df["startIsSource"].fillna(False)
    calc_df["endIsOutlet"] = calc_df["endIsOutlet"].fillna(False)
    calc_df[['node_id_start', 'node_id_end']] = calc_df[['node_id_start', 'node_id_end']].astype(str)
    ids_count = pd.concat((calc_df['node_id_start'], calc_df['node_id_end'])).value_counts()
    ids_count.rename('ids_count')
    calc_df = calc_df.join(ids_count.to_frame(), on='node_id_start', how='left')
    calc_df = calc_df.rename(columns={0: 'start_id_count'})
    ids_count = calc_df['node_id_start'].value_counts()
    ids_count.rename('ids_count')
    calc_df = calc_df.join(ids_count.to_frame().rename(columns={"node_id_start": 0}), on='node_id_end', how='left')
    calc_df = calc_df.rename(columns={0: 'end_id_count'})
    calc_df['end_id_count'] = calc_df["end_id_count"].fillna(0)
    calc_df['sourceByCount'] = calc_df['start_id_count'] == 1
    calc_df['outletByCount'] = calc_df['end_id_count'] == 0
    calc_df['sourceMistakes'] = calc_df['sourceByCount'] == calc_df['startIsSource']
    calc_df['outletMistakes'] = calc_df['outletByCount'] == calc_df['endIsOutlet']
    calc_df['sourceValueIsFilled'] = pd.notna(calc_df['startValue'])
    calc_df['outletValueIsFilled'] = pd.notna(calc_df['endValue'])
    calc_df['sourceKindIsFilled'] = pd.notna(calc_df['startKind'])
    calc_df['outletKindIsFilled'] = pd.notna(calc_df['endKind'])
    calc_df['inletBoundaryMistakes'] = True
    calc_df['outletBoundaryMistakes'] = True
    calc_df.loc[calc_df["startIsSource"], 'inletBoundaryMistakes'] = calc_df[calc_df["startIsSource"]][
                                                                         'sourceValueIsFilled'] & \
                                                                     calc_df[calc_df["startIsSource"]][
                                                                         'sourceKindIsFilled']
    calc_df.loc[calc_df["endIsOutlet"], 'outletBoundaryMistakes'] = calc_df[calc_df["endIsOutlet"]][
                                                                        'outletValueIsFilled'] & \
                                                                    calc_df[calc_df["endIsOutlet"]][
                                                                        'outletKindIsFilled']
    calc_df['sumOfCounts'] = calc_df['start_id_count'] + calc_df['end_id_count']
    calc_df = calc_df[calc_df['sumOfCounts'] >= 2]
    mistakes_df = calc_df[
        (~calc_df['sourceMistakes']) | (~calc_df['outletMistakes']) | (~calc_df['inletBoundaryMistakes']) | (
            ~calc_df['outletBoundaryMistakes'])]
    if not mistakes_df.empty:
        print(
            f"Following nodes: {mistakes_df[~mistakes_df['sourceMistakes']]['node_id_start'].values} should be sources")
        print(
            f"Following nodes: {mistakes_df[~mistakes_df['outletMistakes']]['node_id_start'].values} should be outlets")
        print(
            f"Start kind and value for following nodes: {mistakes_df[~mistakes_df['inletBoundaryMistakes']]['node_id_start'].values} should be filled")
        print(
            f"End kind and value for following nodes: {mistakes_df[~mistakes_df['outletBoundaryMistakes']]['node_id_end'].values} should be filled")
        assert False

    calc_df = calc_df.reset_index() # Damn you, pandas black emperor!
    return calc_df


def populate_wells_df(HKT, dict_list, inclination, wells_df):
    for i, row in wells_df.iterrows():
        dict_list = process_well_row(HKT, dict_list, inclination, row)
    tempdf = pd.DataFrame.from_dict(dict_list)
    return tempdf


def process_well_row(HKT, dict_list, inclination, row):
    try:
        wellNum = str(int(row["wellNum"]))
    except:
        wellNum = row["wellNum"]
    # padNum = str(int(row["padNum"]))
    padNum = row["padNum"]
    pumpdepth = row["pumpDepth"]
    perforation = row["perforation"]
    tubing = inclination[inclination["wellNum"] == wellNum]
    tubing = tubing.sort_values(by='depth')
    tubing["Roughness"] = 3e-5
    tubing["IntDiameter"] = 0.57
    tubing["NKTlength"] = 10
    local_HKT = HKT[HKT['wellNum'] == wellNum]
    fulldepth = 0
    stages = local_HKT['stageNum'].unique()
    for stageNum in stages:
        stage_HKT = local_HKT[local_HKT["stageNum"] == stageNum]
        stage_HKT = stage_HKT[stage_HKT["_time"] == stage_HKT["_time"].max()]
        fulldepth += stage_HKT["stageLength"].iloc[0]
        tubing.loc[tubing["depth"] <= fulldepth, "IntDiameter"] = (stage_HKT["stageDiameter"].iloc[0] - 16) / 1000
    pump_place = tubing[abs(tubing["depth"] - pumpdepth) == min(abs(tubing["depth"] - pumpdepth))].iloc[0]
    perforation_place = tubing[abs(tubing["depth"] - perforation) == min(abs(tubing["depth"] - perforation))].iloc[0]

    dict_list = make_well_parts_rows(dict_list, padNum, perforation_place, pump_place, row, wellNum)
    return dict_list


def make_well_parts_rows(dict_list, padNum, perforation_place, pump_place, row, wellNum):
    row_d = row.to_dict()
    # plast-zaboi
    dct = row_d.copy()
    dct["juncType"] = "plast"
    dct["effectiveD"] = 1
    dct["roughness"] = 3e-5
    dct["node_id_start"] = f"PAD_{padNum}_WELL_{wellNum}"
    dct["node_id_end"] = f"PAD_{padNum}_WELL_{wellNum}_zaboi"
    dct["endIsOutlet"] = False
    dict_list += [dct]
    # zaboi-intake
    dct = dct.copy()
    dct["startIsSource"] = False
    dct["juncType"] = "pipe"
    dct["node_id_start"] = f"PAD_{padNum}_WELL_{wellNum}_zaboi"
    dct["node_id_end"] = f"PAD_{padNum}_WELL_{wellNum}_pump_intake"
    absdiff = perforation_place["absMark"] - pump_place["absMark"]
    Ldiff = perforation_place["prolongation"] - pump_place["prolongation"]
    dct["L"] = Ldiff
    dct['uphillM'] = absdiff
    dct["intD"] = 0.127
    dict_list += [dct]
    # wellpump
    dct = dct.copy()
    dct["juncType"] = "wellpump"
    dct["node_id_start"] = f"PAD_{padNum}_WELL_{wellNum}_pump_intake"
    dct["node_id_end"] = f"PAD_{padNum}_WELL_{wellNum}_pump_outlet"
    dict_list += [dct]
    # outlet-wellhead
    dct = dct.copy()
    dct["juncType"] = "pipe"
    dct["node_id_start"] = f"PAD_{padNum}_WELL_{wellNum}_pump_outlet"
    dct["node_id_end"] = f"PAD_{padNum}_WELL_{wellNum}_wellhead"
    absdiff = pump_place["absMark"]
    Ldiff = pump_place["prolongation"]
    dct["L"] = Ldiff
    dct['uphillM'] = absdiff
    dct["intD"] = pump_place["IntDiameter"]
    dict_list += [dct]
    # wellhead-pad
    dct = dct.copy()
    dct["juncType"] = "pipe"
    dct["node_id_start"] = f"PAD_{padNum}_WELL_{wellNum}_wellhead"
    dct["node_id_end"] = row["node_id_end"]
    dct["endIsOutlet"] = row["endIsOutlet"]
    dct["L"] = 50
    dct['uphillM'] = 0
    dct["intD"] = 0.1
    dict_list += [dct]
    return dict_list


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