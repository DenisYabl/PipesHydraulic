import colorama
from colorama import Back, Style

from GraphEdges.HE2_Pipe import HE2_WaterPipe
from GraphEdges.HE2_WellPump import HE2_WellPump
from GraphEdges.HE2_Plast import HE2_Plast
from GraphNodes import HE2_Vertices as vrtxs
import networkx as nx
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
from collections import namedtuple
from Fluids.HE2_Fluid import HE2_BlackOil
from Tools.HE2_Logger import check_for_nan, getLogger

logger = getLogger(__name__)

CheckSolutionResults = namedtuple('CheckSolutionResults', ['first_CL_resd', 'second_CL_resd', 'negative_P', 'misdirected_flow', 'bad_directions', 'first_CL_OWG_resd'])


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

def generate_random_net_v1(N=15, E=20, SRC=3, SNK=3, P_CNT=1, Q=20, P=200, D=0.5, H=50, L=1000, RGH=1e-4, SEGS=1,
                           randseed=None):
    '''
    :param N: total nodes count
    :param E: total edges count, cannot be less than N-1
    :param SRC: sources count
    :param SNK: sinks count
    :param P_CNT: count of nodes with fixed pressure
    :param Q: maximum boundary flow on source (not sink!)  node
    :param P: maximum boundary pressure on one of source/sink node
    :param D: maximum pipe segment inner diameter on graph edge
    :param H: maximum pipe segment slope
    :param L: maximum pipe segment length
    :param SEGS: maximum pipe segments count
    :return: MultiDiGraph
    This method produce water pipelines network with some pressure constrained node, with constant temperature
    and with some sources/sinks. Result will be a non-tree, linked multigraph with objects attached to nodes and edges
    '''
    np.random.seed(randseed)
    E = max(E, N - 1)
    B = SRC + SNK
    J = N - B
    juncs = {f'junc_{i}': vrtxs.HE2_ABC_GraphVertex() for i in range(J)}

    kinds = ['P']*P_CNT + ['Q']*(B-P_CNT)
    srcs = [True] * SRC + [False] * SNK
    np.random.shuffle(kinds)
    total_q = np.random.uniform(Q * B)
    src_q = np.random.uniform(0, 1, SRC)
    src_q = total_q * src_q / sum(src_q)
    snk_q = np.random.uniform(0, 1, SNK)
    snk_q = total_q * snk_q / sum(snk_q)
    qs = list(src_q) + list(snk_q)

    p_nodes, sources, sinks = dict(), dict(), dict()
    for kind, src, q in zip(kinds, srcs, qs):
        if src and kind == 'P':
            node = vrtxs.HE2_Source_Vertex(kind, np.random.randint(-P, P), 'water', 20)
            name = f'p_node_{len(p_nodes)}'
            p_nodes[name] = node
        elif src and kind == 'Q':
            node = vrtxs.HE2_Source_Vertex(kind, q, 'water', 20)
            name = f'src_{len(sources)}'
            sources[name] = node
        elif not src and kind == 'P':
            node = vrtxs.HE2_Boundary_Vertex(kind, np.random.randint(-P, P))
            name = f'p_node_{len(p_nodes)}'
            p_nodes[name] = node
        elif not src and kind == 'Q':
            node = vrtxs.HE2_Boundary_Vertex(kind, q)
            name = f'snk_{len(sinks)}'
            sinks[name] = node

    nodes = {**p_nodes, **sources, **sinks, **juncs}
    assert len(nodes) == N
    mapping = dict(zip(range(N), nodes.keys()))
    RT = nx.generators.trees.random_tree(N, seed=randseed)
    edgelist = [tuple(np.random.choice([u, v], 2, replace=False)) for u, v in RT.edges]
    edgelist += [tuple(np.random.choice(range(N), 2, replace=False)) for i in range(E-(N-1))]

    MDRG = nx.MultiDiGraph(edgelist)
    G = nx.relabel.relabel_nodes(MDRG, mapping)
    nx.set_node_attributes(G, name='obj', values=nodes)

    pipes = dict()
    for u, v, k in G.edges:
        segs = np.random.randint(SEGS) + 1
        Ls = np.random.uniform(1e-5, L, segs)
        Hs = np.random.uniform(-H, H, segs)
        Ds = np.random.uniform(1e-5, D, segs)
        Rs = np.random.uniform(0, RGH, segs)
        # print(u, v, Ls, Hs, Ds, Rs)
        pipe = HE2_WaterPipe(Ls, Hs, Ds, Rs)
        pipes[(u, v, k)] = pipe
    nx.set_edge_attributes(G, name='obj', values=pipes)

    assert nx.algorithms.components.number_connected_components(nx.Graph(G)) == 1

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
    # pos = nx.drawing.circular_layout(G)
    g_nodes = set(G.nodes)
    params = zip([p_nodes, sources, sinks, juncs], [50, 50, 50, 10], ['red', 'blue','blue','black'], [[], ['Q'], ['Q'], []])
    # params = zip([sources, sinks], [50, 10], ['blue','black'], [['Q'], ['Q']])
    label_pos = {k:(pos[k][0] + shifts[k][0], pos[k][1] + shifts[k][1]) for k in pos} if shifts is not None else pos
    for nodelist, node_size, node_color, ks in params:
        if nodelist:
            nx.draw_networkx_nodes(G, nodelist=list(set(nodelist) & g_nodes), node_size=node_size, node_color=node_color, ax=ax, pos=pos, alpha=0.3)
            HE2_draw_node_labels(G, g_nodes, list(set(nodelist) & g_nodes), keys=['P_bar']+ks, ax=ax, pos=label_pos)


    # edge_labels = {(u,v): str(G[u][v]['obj'])+f"\n{G[u][v]['obj'].result['x']:.2f}" for u, v in G.edges()}
    edgs = [(u, v) for (u, v) in G.edges() if 'result' in G[u][v]['obj'].__dict__]
    edge_labels = {(u,v): f"\n{G[u][v]['obj'].result['x']:.2f}" for u, v in edgs}
    nx.draw_networkx_edge_labels(G, pos=pos, edge_labels=edge_labels, font_size=7)

    # if type(G) == nx.MultiDiGraph:
    #     edge_labels = {(u,v): str(G[u][v][k]['obj'])+f"\n{G[u][v][k]['obj'].result['x']:.2f}" for u, v, k in G.edges}
    # else:
    #     edge_labels = {(u, v): str(G[u][v]['obj']) + f"\n{G[u][v]['obj'].result['x']:.2f}" for u, v in G.edges}

    if type(G) == nx.MultiDiGraph:
        edge_labels = {(u,v): f"{G[u][v][k]['obj'].result['x']:.2f}" for u, v, k in G.edges}
    else:
        edge_labels = {(u, v): f"{G[u][v]['obj'].result['x']:.2f}" for u, v in edgs}

    nx.draw_networkx_edge_labels(G, pos=pos, edge_labels=edge_labels, font_size=9, alpha=0.3 )

    nx.draw_networkx_edges(G, pos=pos, width=2, ax=ax, edge_color='black')
    plt.show()


def evaluate_1stCL_residual(graph):
    G = nx.MultiDiGraph(graph)
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

    for u, v, k in G.edges:
        x = G[u][v][k]['obj'].result['x']
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

def evaluate_2ndCL_residual(graph):
    G = nx.MultiDiGraph(graph)
    residual = 0
    for (u, v, k) in G.edges:
        u_obj = G.nodes[u]['obj']
        v_obj = G.nodes[v]['obj']
        edge_obj = G[u][v][k]['obj']
        x = edge_obj.result['x']
        p_u = u_obj.result['P_bar']
        t_u = u_obj.result['T_C']
        p_v = v_obj.result['P_bar']
        t_v = v_obj.result['T_C']
        p, t = edge_obj.perform_calc_forward(p_u, t_u, x)
        edge_residul = abs(p - p_v)
        residual += edge_residul
    return residual

def evalute_pressures_below_zero(graph, P_threshold_bar = 0):
    G = nx.MultiDiGraph(graph)
    rez = 0
    for n in G.nodes:
        obj = G.nodes[n]['obj']
        P = obj.result['P_bar']
        if P < P_threshold_bar:
            rez += P_threshold_bar - P
            logger.info(f'node {n}, P={P}')
    return rez

def check_directions(graph):
    violations_count = 0
    violations_flow_sum = 0
    G = nx.MultiDiGraph(graph)
    for n in G.nodes:
        obj = G.nodes[n]['obj']
        if isinstance(obj, vrtxs.HE2_Boundary_Vertex) and obj.is_source:
            Q = obj.result['Q']
            if Q < 0:
                logger.info(f'source node {n}, Q={Q}')
                violations_flow_sum += abs(Q)
                violations_count += 1

    for (u, v, k) in G.edges:
        edge_obj = G[u][v][k]['obj']
        if isinstance(edge_obj, HE2_Plast) or isinstance(edge_obj, HE2_WellPump):
            x = edge_obj.result['x']
            if x < 0:
                logger.info(f'well edge {u} - {v}, x={x}')
                violations_flow_sum += abs(x)
                violations_count += 1

    return violations_flow_sum, violations_count


def split_mass_flow_to_OWG(x, fluid: HE2_BlackOil):
    oil_params = fluid.oil_params
    oil_ro, wat_ro, gas_ro = oil_params.oildensity_kg_m3, oil_params.waterdensity_kg_m3, oil_params.gasdensity_kg_m3
    wc = oil_params.volumewater_percent / 100
    gf = oil_params.gasFactor
    owg_mix_pseudo_density = oil_ro * (1 - wc) + wat_ro * wc + gas_ro * (1 - wc) * gf
    Q_owg = x / owg_mix_pseudo_density
    Qo = (1 - wc) * Q_owg
    Qw = wc * Q_owg
    Qg = (1 - wc) * gf* Q_owg
    Xo = Qo * oil_ro
    Xw = Qw * wat_ro
    Xg = Qg * gas_ro
    if abs(x - Xo - Xw - Xg) > 1e-4:
        logger.warning('OWG mass flows doesnt match total mass flow!')
        raise ValueError
    return Xo, Xw, Xg


def check_1stKL_by_OWG_separately(graph):
    G = nx.MultiDiGraph(graph)
    Qo_dict = {}
    Qw_dict = {}
    Qg_dict = {}
    Xo_sum_dict = dict(zip(G.nodes, [0]*len(G.nodes)))
    Xw_sum_dict = dict(zip(G.nodes, [0]*len(G.nodes)))
    Xg_sum_dict = dict(zip(G.nodes, [0]*len(G.nodes)))
    p_nodes = []
    for n in G.nodes:
        obj = G.nodes[n]['obj']
        if isinstance(obj, vrtxs.HE2_Boundary_Vertex) and obj.kind == 'P':
            p_nodes += [n]
            continue

        Qo_dict[n] = 0
        Qw_dict[n] = 0
        Qg_dict[n] = 0
        if isinstance(obj, vrtxs.HE2_Boundary_Vertex) and obj.kind == 'Q':
            Q = obj.value if obj.is_source else -obj.value
            Qo, Qw, Qg = split_mass_flow_to_OWG(Q, obj.fluid)
            Qo_dict[n] = Qo
            Qw_dict[n] = Qw
            Qg_dict[n] = Qg

    for u, v, k in G.edges:
        obj = G[u][v][k]['obj']
        x = obj.result['x']
        xo, xw, xg = split_mass_flow_to_OWG(x, obj.fluid)

        Xo_sum_dict[u] -= xo
        Xw_sum_dict[u] -= xw
        Xg_sum_dict[u] -= xg

        Xo_sum_dict[v] += xo
        Xw_sum_dict[v] += xw
        Xg_sum_dict[v] += xg

    res_o, res_w, res_g = 0, 0, 0
    for n in Qo_dict:
        res_o += abs(Qo_dict[n] + Xo_sum_dict[n])
        res_w += abs(Qw_dict[n] + Xw_sum_dict[n])
        res_g += abs(Qg_dict[n] + Xg_sum_dict[n])
        pass

    Qo_net_balance = sum(Qo_dict.values())
    Qw_net_balance = sum(Qw_dict.values())
    Qg_net_balance = sum(Qg_dict.values())
    p_xo_sum, p_xw_sum, p_xg_sum = 0, 0, 0
    for n in p_nodes:
        p_xo_sum += Xo_sum_dict[n]
        p_xw_sum += Xw_sum_dict[n]
        p_xg_sum += Xg_sum_dict[n]

    res_o += abs(p_xo_sum - Qo_net_balance)
    res_w += abs(p_xw_sum - Qw_net_balance)
    res_g += abs(p_xg_sum - Qg_net_balance)

    residual = res_o + res_w + res_g

    return residual


def check_solution(G):
    '''
    :param G: HE2 graph with results
    :return: CheckSolutionResults named tuple
    This method checks solution to violate physics and technical constraints, like as negative pressures and misdirected flows
    CheckSolutionResults contains next fields:
        first_CL_resd: this is residual of 1st Kirchhoff law. Expected value - very small positive number, 1e-6. 0.1 is very strong violation
        second_CL_resd: this is residual of 2nd Kirchhoff law. You can ignore this value, and use solver.op_result.succes and solver.op_result.fun instead
        negative_P: Solver looks for solution in negative pressure area, cause get solution with negative pressures can be much useful than no solution without any additional info.
            So, negative_P field is sum of abs() of all negative pressures at graph nodes.
        bad_directions: Some edges and nodes on graph, like plast or pump, have obvious flow directions, that really cannot be violated. So this field just counts nodes and edges with violated direction.
        misdirected_flow: This field sum all flows on all nodes and edges with violated directions
    '''
    res1 = evaluate_1stCL_residual(G)
    res2 = evaluate_2ndCL_residual(G)
    res3 = evalute_pressures_below_zero(G)
    res4, res5 = check_directions(G)
    res6 = check_1stKL_by_OWG_separately(G)
    rez = CheckSolutionResults(res1, res2, res3, res4, res5, res6)
    return rez

def check_fluid_mixation(G, x_dict, cocktails, sources):
    S = len(sources)
    for n in G.nodes:
        x_in, x_out = 0, 0
        for e in G.in_edges(n):
            x_in += x_dict.get(e, 0)
        for e in G.out_edges(n):
            x_out += x_dict.get(e, 0)

        out_cktls = []
        if n in cocktails:
            out_cktls += [cocktails[n]]
        for e in G.out_edges(n):
            if e in cocktails:
                out_cktls += [cocktails[e]]

        if len(out_cktls) == 0:
            assert x_in == 0
            assert x_out == 0
            continue

        for c1, c2 in zip(out_cktls, out_cktls[1:]):
            if np.linalg.norm(c1-c2) > 1e-5:
                return False

        out_cktl = out_cktls[0]
        in_cktl = np.zeros(S)
        if n in sources:
            assert x_out >= x_in
            idx = sources.index(n)
            in_cktl[idx] = x_out - x_in

        for e in G.in_edges(n):
            if e in cocktails:
                x = x_dict.get(e, 0)
                in_cktl += x * cocktails[e]

        in_cktl /= sum(in_cktl)
        if np.linalg.norm(in_cktl - out_cktl) > 1e-5:
            return False

    return True


def build_dual_schema_from_solved(schema, p_nodes, sources, sinks, juncs):
    G = nx.MultiDiGraph(schema, data=True)
    nodes_P, nodes_Q = {}, {}

    X_sum_dict = dict(zip(G.nodes, [0]*len(G.nodes)))
    for u, v, k in G.edges:
        x = G[u][v][k]['obj'].result['x']
        X_sum_dict[u] -= x
        X_sum_dict[v] += x

    for n in p_nodes:
        nodes_P[n] = G.nodes[n]['obj'].P
        nodes_Q[n] = -X_sum_dict[n]

    for n in juncs:
        nodes_Q[n] = 0
    for n in sources:
        nodes_Q[n] = G.nodes[n]['obj'].Q
    for n in sinks:
        nodes_Q[n] = -G.nodes[n]['obj'].Q

    for n in {**sources, **sinks, **juncs}:
        obj = G.nodes[n]['obj']
        nodes_P[n] = obj.result['P_bar']

    newschema = nx.MultiDiGraph()
    p_n = np.random.choice(G.nodes)
    for n in G.nodes:
        kind = np.random.choice(['P', 'Q'], p=[0.2, 0.8])
        if n == p_n:
            kind = 'P'
        P = nodes_P[n]
        Q = nodes_Q[n]
        value = abs(Q) if kind=='Q' else P
        obj = None

        if (Q<0) or (Q==0 and kind == 'P'):
            obj = vrtxs.HE2_Boundary_Vertex(kind, value)
        elif Q==0 and kind == 'Q':
            obj = vrtxs.HE2_ABC_GraphVertex()
        elif Q>0:
            obj = vrtxs.HE2_Source_Vertex(kind, value, 'water', 20)

        newschema.add_node(n, obj=obj)

    for u, v, k in G.edges:
        obj = G[u][v][k]['obj']
        dxs, dys, diams, rghs = [], [], [], []
        for seg in obj.segments:
            dxs += [seg.dx_m]
            dys += [seg.uphill_m]
            diams += [seg.inner_diam_m]
            rghs += [seg.roughness_m]

        newobj = HE2_WaterPipe(dxs, dys, diams, rghs)
        newschema.add_edge(u, v, obj=newobj)

    return newschema

def generate_superpositioned_colored_flows_graph(N=10, E=13, SRC=3, SNK=3, randseed=None):
    return
    # Bullshit. Main idea doesnt work
    np.random.seed(randseed)
    RT = nx.generators.trees.random_tree(N, seed=randseed)
    edge_list = [tuple(np.random.choice([u, v], 2, replace=False)) for u, v in RT.edges]
    edge_set = set(edge_list)
    while len(edge_set) < E:
        u, v = tuple(np.random.choice(range(N), 2, replace=False))
        if ((v, u) in edge_set) or ((v, u) in edge_set):
            continue
        edge_set |= {(u, v)}
    edge_list = list(edge_set)

    DRG = nx.DiGraph()
    DRG.add_nodes_from(RT.nodes)
    DRG.add_edges_from(edge_list)
    base = DRG
    can_be_source, can_be_sink = set(), set()
    for n in base.nodes:
        if len(base.in_edges(n)) > 0:
            can_be_sink |= {n}
        if len(base.out_edges(n)) > 0:
            can_be_source |= {n}
    A = can_be_source - can_be_sink
    B = can_be_source - A
    srcs = list(A) + list(B)
    sources = {srcs[i] for i in range(SRC)}

    can_be_sink -= sources
    A = can_be_sink - can_be_source
    B = can_be_sink - A
    snks = list(A) + list(B)
    sinks = {snks[i] for i in range(SNK)}

    scaffold = nx.DiGraph(base)
    scaffold_nodes = ['SUPERSOURCE', 'SUPERSINK']
    scaffold.add_nodes_from(scaffold_nodes)
    scaffold_edges = []

    for n in sources:
        scaffold_edges += [('SUPERSOURCE', n)]
    src_edges = scaffold_edges[:]
    sink_edges = []
    for n in sinks:
        sink_edges += [(n, 'SUPERSINK')]
    scaffold_edges += sink_edges
    scaffold.add_edges_from(scaffold_edges)
    zero_capacity = {e:0 for e in src_edges}
    all_flows, Q = {}, {}
    for n in sources:
        nx.set_edge_attributes(scaffold, name='capacity', values=zero_capacity)
        scaffold['SUPERSOURCE'][n]['capacity'] = 1000
        capacities = {e:np.random.randint(20, 100) for e in base.edges}
        # capacities.update({e:np.random.randint(0, 1000) for e in sink_edges})
        nx.set_edge_attributes(scaffold, name='capacity', values=capacities)
        flow_value, flow_dict = nx.algorithms.flow.maximum_flow(scaffold, 'SUPERSOURCE', 'SUPERSINK')
        all_flows[n] = flow_dict
        Q[n] = flow_value

    rez = {}
    for n, flow_dict in all_flows.items():
        for u, v_x_dict in flow_dict.items():
            for v, x in v_x_dict.items():
                xs = rez.get((u, v), [])
                xs += [x]
                rez[(u, v)] = xs

    result, x_dict = {}, {}
    for (u, v), xs in rez.items():
        if len({u, v} & {'SUPERSOURCE', 'SUPERSINK'}) > 0:
            continue
        assert len(xs) == len(sources)
        arr = np.array(xs)
        sss = sum(arr)
        x_dict[(u, v)] = sss
        if sss == 0:
            continue
        arr = arr / sss
        result[(u, v)] = arr
        print(u, v, arr)



    return result, sources, base, x_dict


def split_input_df_to_pipes_and_boundaries(df):
    bnd_cols1 = ['kind', 'Q', 'is_source', 'P']
    bnd_cols2 = ['node_id', 'node_name']
    start_cols = [col + '_start' for col in bnd_cols2] + ['start_' + col for col in bnd_cols1]
    end_cols = [col + '_end' for col in bnd_cols2] + ['end_' + col for col in bnd_cols1]
    df_bnd_start = df[start_cols]
    df_bnd_start.columns = bnd_cols2 + bnd_cols1
    df_bnd_end = df[end_cols]
    df_bnd_end.columns = bnd_cols2 + bnd_cols1
    df_bnds = pd.concat([df_bnd_start, df_bnd_end])
    df_bnds = df_bnds[~df_bnds.kind.isna()]
    to_drop = start_cols + end_cols
    df_pipes = df.drop(columns=to_drop)
    return df_pipes, df_bnds

def split_result_df_to_pipes_and_nodes(df):
    bnds_cols = ['kind', 'Q', 'is_source', 'P']
    node_cols = ['node_id', 'node_name', 'node_type', 'x', 'y']
    start_cols = [col + '_start' for col in node_cols] + ['result_start_P']
    end_cols = [col + '_end' for col in node_cols] + ['result_end_P']
    df_node_start = df[start_cols]
    df_node_start.columns = node_cols + ['result_P']
    df_node_end = df[end_cols]
    df_node_end.columns = node_cols + ['result_P']
    df_nodes = pd.concat([df_node_start, df_node_end]).drop_duplicates()
    node_id_set = set(df.node_id_start.values) | set(df.node_id_end.values)
    print(len(node_id_set))
    from collections import Counter
    cntr = Counter(df_nodes.node_id)
    print(cntr.most_common(10))



    to_drop = start_cols + end_cols
    df_pipes = df.drop(columns=to_drop)
    # assert len(df_nodes) == len(df_nodes.node_id.unique())
    return df_pipes, df_nodes

def make_oilupstream_graph_layout(G):
    nodes = list(G.nodes)
    long_wells = {}
    for n in nodes:
        if len(G.in_edges(n)) == 0 and len(G.out_edges(n)) == 1:
            long_wells[n]=[]

    for l in long_wells:
        n = l
        while True:
            if not (len(G.in_edges(n)) <= 1 and len(G.out_edges(n)) == 1):
                break
            u, n = list(G.out_edges(n))[0]
            long_wells[l] += [n]

    pads = dict()
    for l in long_wells:
        p = long_wells[l][-1]
        wells = pads.get(p, [])
        wells += [l]
        pads[p] = wells

    R, r = 1, 0.5
    rad_pos = dict()
    if len(long_wells) > 0:
        well_step = 2*np.pi / len(long_wells)
        well_a = 0
        for p in pads:
            for w in pads[p]:
                well_points = [w] + long_wells[w][:-1]
                nn = len(well_points)
                if nn > 1:
                    step = (R-r)/(nn-1)
                    for i, n in enumerate(well_points):
                        rad_pos[n] = (well_a, R-i*step)
                else:
                    rad_pos[w] = (well_a, R)
                well_a += well_step
            w0 = pads[p][0]
            w1 = pads[p][-1]
            pad_a = (rad_pos[w0][0] + rad_pos[w1][0])/2
            rad_pos[p] = (pad_a, r)

    wells_and_pads = list(rad_pos.keys())
    other_nodes = set(G.nodes) - set(wells_and_pads)
    other_G = nx.DiGraph()
    other_G.add_nodes_from(other_nodes)
    for u, v in G.edges():
        if u in other_nodes and v in other_nodes:
            other_G.add_edge(u, v)
    other_pos = nx.kamada_kawai_layout(other_G)

    pos = dict()
    for n in rad_pos:
        a, rr = rad_pos[n]
        x = np.sin(a) * rr
        y = np.cos(a) * rr
        pos[n] = np.array([x, y])

    for n in other_pos:
        pos[n] = other_pos[n] * r

    return pos


def print_wells_pressures(G, wells):
    colorama.init()
    table_header = '                                            bottom   intake   outlet  wellhead'
    print(table_header)
    for pad_well in wells:
        l = pad_well.split('_')
        if len(l) < 4:
            continue
        pad, well = l[1], l[3]
        well_subgraph, well_nodes = cut_single_well_subgraph(G, pad, well)
        row_header = 'pad ' + f' {pad}'[-2:] + ' well ' + f'  {well}'[-4:] + ', from plast:   '
        print(row_header, end=' ')
        for n in well_nodes:
            P = G.nodes[n]['obj'].result['P_bar']
            prefix = Back.RED if P <= 1 else Style.RESET_ALL
            print(prefix + f'{P:8.3f}', end= ' ')
        print(Style.RESET_ALL + '   up to pad collector')


def print_solution(G):
    colorama.init()
    table_header = f' {"start":>30} {"P_bar":>7} {"Q kg/s":>8}   {"end":>30} {"P_bar":>7} {"Q kg/s":>8}    {"X kg/s":>7}'
    print(table_header)
    for e in G.edges:
        if len(e) == 2:
            u, v = e
            obj = G[u][v]['obj']
        elif len(e) == 3:
            u, v, k = e
            obj = G[u][v][k]['obj']
        else:
            return

        u_obj = G.nodes[u]['obj']
        v_obj = G.nodes[v]['obj']
        x = obj.result['x']
        p_u = u_obj.result['P_bar']
        pu_str = f'{p_u:8.3f}'
        if p_u <= 1:
            pu_str = Back.RED + pu_str + Style.RESET_ALL

        p_v = v_obj.result['P_bar']
        pv_str = f'{p_v:8.3f}'
        if p_v <= 1:
            pv_str = Back.RED + pv_str + Style.RESET_ALL

        q_u = u_obj.result['Q']
        q_u_str = ''
        if abs(q_u) > 1e-5:
            q_u_str = f"{q_u:8.3f}"

        q_v = v_obj.result['Q']
        q_v_str = ''
        if abs(q_v) > 1e-5:
            q_v_str = f"{q_v:8.3f}"

        row = f' {u:>30} {pu_str:>8} {q_u_str:>8}   {v:>30} {pv_str:>8} {q_v_str:>8}    {x:7.3f}{Style.RESET_ALL}'
        print(row)


def cut_single_well_subgraph(G, pad_name, well, nodes = None):
    rez = nx.DiGraph()  # Di = directed
    if nodes == None:
        nms = dict()
        nms.update(plast = f'PAD_{pad_name}_well_{well}')
        nms.update(zaboi = f'Zaboi_{well}')
        nms.update(intake = f'Pump_intake_{well}')
        nms.update(outlet = f'Pump_outlet_{well}')
        nms.update(wellhead = f'Wellhead_{well}')
        nms.update(pad = f'PAD_{pad_name}')
        nodes = [nms['plast'], nms['zaboi'], nms['intake'], nms['outlet'], nms['wellhead'], nms['pad']]

    edgelist = []
    for i in range(len(nodes)-1):
        edgelist += [(nodes[i], nodes[i+1])]
    rez.add_nodes_from(nodes)
    rez.add_edges_from(edgelist)
    node_objs = {n:G.nodes[n]['obj'] for n in nodes[:-1]}
    node_objs[nodes[-1]] = vrtxs.HE2_Boundary_Vertex('P', 5)
    edge_objs = {}
    if isinstance(G, nx.MultiDiGraph):
        for u, v in edgelist:
            obj = G[u][v][0]['obj']
            edge_objs[(u, v)] = obj
    elif isinstance(G, nx.DiGraph):
        for u, v in edgelist:
            obj = G[u][v]['obj']
            edge_objs[(u, v)] = obj

    nx.set_node_attributes(rez, name='obj', values=node_objs)
    nx.set_edge_attributes(rez, name='obj', values=edge_objs)
    return rez, nodes