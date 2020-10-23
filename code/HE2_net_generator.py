from HE2_Pipe import HE2_WaterPipe
import HE2_Vertices as vrtxs
import networkx as nx
import numpy as np

def generate_random_net_v0(N=15, E=20, SRC=3, SNK=3, Q=20, P=200, D=0.5, H=50, L=1000, RGH=1e-4, SEGS=10,
                           randseed=None):
    '''
    :param N: total nodes count
    :param E: total edges count, cannot be less than N-1
    :param SRC: sources count
    :param SNK: sinks count
    :param Q: maximum boundary flow on source (not sink!)  node
    :param P: maximum boundary pressure on one of source/sink node
    :param D: maximum pipe segment inner diameter on graph edge
    :param H: maximum pipe segment slope
    :param L: maximum pipe segment length
    :param SEGS: maximum pipe segments count
    :return: DiGraph (not Multi)
    This method produce water pipelines network with one pressure constrained node, with constant temperature
    and with some sources/sinks. Network graph can contains cycles
    '''
    np.random.seed(randseed)
    E = max(E, N - 1)
    J = N - SRC - SNK
    juncs = {f'junc_{i}': vrtxs.HE2_ABC_GraphVertex() for i in range(J)}
    coin = np.random.randint(0, 2)
    pressure = np.random.randint(-P, P)
    if coin == 1:
        p_nodes = {'p_node_0': vrtxs.HE2_Source_Vertex('P', pressure, 'water', 20)}
    else:
        p_nodes = {'p_node_0': vrtxs.HE2_Boundary_Vertex('P', pressure)}

    src_q = np.random.uniform(0, Q, SRC)
    total_q = sum(src_q)
    sink_q = np.random.uniform(0, 1, SNK)
    sink_q = total_q * sink_q / sum(sink_q)
    np.testing.assert_almost_equal(total_q, sum(sink_q))

    SRC = SRC - coin
    SNK = SNK - (1 - coin)
    sources = {f'src_{i}': vrtxs.HE2_Source_Vertex('Q', src_q[i], 'water', 20) for i in range(SRC)}
    sinks = {f'sink_{i}': vrtxs.HE2_Boundary_Vertex('Q', -sink_q[i]) for i in range(SNK)}

    nodes = {**p_nodes, **sources, **sinks, **juncs}
    mapping = dict(zip(range(len(nodes)), nodes.keys()))
    UDRG = nx.generators.random_graphs.gnm_random_graph(N, E, directed=False, seed=randseed)
    while nx.algorithms.components.number_connected_components(nx.Graph(UDRG)) != 1:
        UDRG = nx.generators.random_graphs.gnm_random_graph(N, E, directed=False)
    edgelist = []
    for u, v in UDRG.edges:
        if np.random.randint(0, 2):
            edgelist += [(u, v)]
        else:
            edgelist += [(v, u)]
    DRG = nx.DiGraph(edgelist)
    G = nx.relabel.relabel_nodes(DRG, mapping)
    nx.set_node_attributes(G, name='obj', values=nodes)

    pipes = dict()
    for u, v in G.edges():
        segs = np.random.randint(SEGS) + 1
        Ls = np.random.uniform(1e-5, L, segs)
        Hs = np.random.uniform(-H, H, segs)
        Ds = np.random.uniform(1e-5, D, segs)
        Rs = np.random.uniform(0, RGH, segs)
        pipe = HE2_WaterPipe(Ls, Hs, Ds, Rs)
        pipes[(u, v)] = pipe
    nx.set_edge_attributes(G, name='obj', values=pipes)

    return G
