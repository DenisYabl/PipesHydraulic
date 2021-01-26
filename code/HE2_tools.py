from HE2_Pipe import HE2_WaterPipe
import HE2_Vertices as vrtxs
import networkx as nx
import numpy as np
import matplotlib.pyplot as plt

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

    # for n in G.nodes():
    #     print(n)

    pipes = dict()
    for u, v in G.edges():
        segs = np.random.randint(SEGS) + 1
        Ls = np.random.uniform(1e-5, L, segs)
        Hs = np.random.uniform(-H, H, segs)
        Ds = np.random.uniform(1e-5, D, segs)
        Rs = np.random.uniform(0, RGH, segs)
        # print(u, v, Ls, Hs, Ds, Rs)
        pipe = HE2_WaterPipe(Ls, Hs, Ds, Rs)
        pipes[(u, v)] = pipe
    nx.set_edge_attributes(G, name='obj', values=pipes)

    return G, dict(p_nodes=p_nodes, juncs=juncs, sources=sources, sinks=sinks)

def HE2_draw_node_labels(G, g_nodes, nodelist, keys, **kwargs):
    lbls = dict()
    for n in list(set(nodelist) & g_nodes):
        sss = [n]
        obj = G.nodes[n]['obj']
        for k in keys:
            try:
                if 'result' in obj.__dict__ and k in obj.result:
                    sss += [f'{obj.result[k]:.2f}']
                elif k in obj.__dict__:
                    sss += [f'{obj.__dict__[k]:.2f}']
            except:
                pass
        lbls.update({n: '\n'.join(sss)})
    nx.draw_networkx_labels(G, labels=lbls, **kwargs, font_size=7)

def draw_solution(G, shifts, p_nodes, sources, sinks, juncs):
    #TODO Не однообразно рисую узлы и дуги, просит рефакторинга
    #TODO Не однообразно формирую лейблы для узлов и для дуг, тоже просит рефакторинга
    fig = plt.figure(constrained_layout=True, figsize=(12, 8))
    ax = fig.add_subplot(1, 1, 1)
    pos = nx.drawing.layout.planar_layout(G)
    g_nodes = set(G.nodes)
    # params = zip([p_nodes, sources, sinks, juncs], [50, 50, 50, 10], ['red', 'blue','blue','black'], [[], ['Q'], ['Q'], []])
    params = zip([sources, sinks], [50, 10], ['blue','black'], [['Q'], ['Q']])
    label_pos = {k:(pos[k][0] + shifts[k][0], pos[k][1] + shifts[k][1]) for k in pos} if shifts is not None else pos
    for nodelist, node_size, node_color, ks in params:
        nx.draw_networkx_nodes(G, nodelist=list(set(nodelist) & g_nodes), node_size=node_size, node_color=node_color, ax=ax, pos=pos)
        HE2_draw_node_labels(G, g_nodes, list(set(nodelist) & g_nodes), keys=['P_bar']+ks, ax=ax, pos=label_pos)

    # edge_labels = {(u,v): str(G[u][v]['obj'])+f"\n{G[u][v]['obj'].result['x']:.2f}" for u, v in G.edges()}
    edge_labels = {(u,v): f"\n{G[u][v]['obj'].result['x']:.2f}" for u, v in G.edges()}
    nx.draw_networkx_edge_labels(G, pos=pos, edge_labels=edge_labels, font_size=7)
    nx.draw_networkx_edges(G, pos=pos, width=2, ax=ax, edge_color='black')
    plt.show()


def evaluate_1stCL_residual(G):
    Q_dict = {}
    X_sum_dict = dict(zip(G.nodes, [0]*len(G.nodes)))
    p_nodes = []
    for n in G.nodes:
        obj = G.nodes[n]['obj']
        if isinstance(obj, vrtxs.HE2_Boundary_Vertex) and obj.kind == 'P':
            p_nodes += [n]
            continue

        Q_dict[n] = 0
        if isinstance(obj, vrtxs.HE2_Boundary_Vertex) and obj.kind == 'Q':
            Q_dict[n] = obj.value if obj.is_source else -obj.value

    for u, v in G.edges:
        x = G[u][v]['obj'].result['x']
        X_sum_dict[u] -= x
        X_sum_dict[v] += x

    residual = 0
    for n in Q_dict:
        residual += abs(Q_dict[n] + X_sum_dict[n])

    Q_net_balance = sum(Q_dict.values())
    p_x_sum = 0
    for n in p_nodes:
        p_x_sum += X_sum_dict[n]
    residual += abs(p_x_sum - Q_net_balance)

    return residual

def evaluate_2ndCL_residual(G):
    residual = 0
    for (u, v) in G.edges():
        u_obj = G.nodes[u]['obj']
        v_obj = G.nodes[v]['obj']
        edge_obj = G[u][v]['obj']
        x = edge_obj.result['x']
        p_u = u_obj.result['P_bar']
        t_u = u_obj.result['T_C']
        p_v = v_obj.result['P_bar']
        t_v = v_obj.result['T_C']
        p, t = edge_obj.perform_calc_forward(p_u, t_u, x)
        residual += abs(p - p_v)
    return residual

def check_solution(G):
    res1 = evaluate_1stCL_residual(G)
    res2 = evaluate_2ndCL_residual(G)
    return res1 + res2