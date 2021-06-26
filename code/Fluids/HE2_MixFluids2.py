import networkx as nx
import numpy as np
from Tools.HE2_ABC import Root
from Tools.HE2_Logger import getLogger
from Tools.HE2_SolverInternalViewer import plot_neighbours_subgraph as plot_nghbs
from GraphNodes.HE2_Vertices import is_junction
from typing import Tuple, Dict

logger = getLogger(__name__)

reduced_graph = None
mappings = None

def evalute_network_fluids_with_root(solver_G, x_dict, have_to_reduce = False):
    # This is like a nested doll
    # Outer layer - solver call with solver graph and xs (flows)
    # Second layer - we remove root node and zero flows, also reverse negative flows
    # Third layer - we reduce graph removing 2-degree nodes and parallel edges, cause it does not affect to solultion
    # Internal layer - on reduced graph we build mass conservation equation matrix and solve it
    # Than, from deep to surface, we rebuild solution on each layer
    global reduced_graph, mappings
    if reduced_graph is None or have_to_reduce:
        reduced_graph, mappings = make_small_graph_for_mixer(solver_G, removeRoot=True)

    G2 = reduced_graph
    x_dict2 = build_x_dict_for_reduced_graph(G2, x_dict, **mappings)
    edges1, edges2, x_dict3 = revert_edges_to_make_all_flows_positive(G2, x_dict2)

    G3 = nx.DiGraph(edges2)
    cocktails3, srcs = evalute_network_fluids_wo_root(G3, x_dict3)

    cocktails2 = restore_reverted_edges_cocktails(cocktails3, edges1, edges2)

    cocktails = restore_cocktail_from_reduced_graph(cocktails2, srcs, x_dict2, **mappings)

    return cocktails, srcs


def revert_edges_to_make_all_flows_positive(G2, x_dict2):
    edges1 = [(u, v) for u, v in G2.edges]
    edges2 = []
    x_dict3 = {}
    for (u, v) in edges1:
        e = (u, v)
        if x_dict2[(u, v)] < 0:
            e = (v, u)
        x_dict3[e] = abs(x_dict2[(u, v)])
        edges2 += [e]
    return edges1, edges2, x_dict3


def restore_reverted_edges_cocktails(cocktails3, edges1, edges2):
    cocktails2 = {}
    for key, cktl in cocktails3.items():
        if (key in edges2) and not (key in edges1):
            u, v = key
            cocktails2[(v, u)] = cktl
        else:
            cocktails2[key] = cktl
    return cocktails2


def restore_cocktail_from_reduced_graph(cocktails3, srcs, x_dict, fluid_em, x_pos_em, x_neg_em):
    cocktails = {}

    # fluid_em = dict(), dict() # Key - edge in reduced graph. Value - list of edges of original graph
    for key, lst in fluid_em.items():
        if not key in cocktails3:
            continue
        cktl = np.around(cocktails3[key], 6)
        for edge in lst:
            cocktails[edge] = cktl

    return cocktails




def build_x_dict_for_reduced_graph(G2, x_dict, fluid_em, x_pos_em, x_neg_em):
    rez_x_dict = dict()
    # pos_em, neg_em = dict(), dict() # Key - edge in reduced graph. Value - list of edges of original graph
    for key, _lst in x_pos_em.items():
        lst = x_pos_em[key]
        s = 0
        for edge in lst:
            s += x_dict[edge]
        rez_x_dict[key] = s

    for key, _lst in x_neg_em.items():
        lst = x_neg_em[key]
        s = 0
        for edge in lst:
            s -= x_dict[edge]
        sss = rez_x_dict.get(key, 0)
        rez_x_dict[key] = sss + s

    return rez_x_dict

def gimme_any_multiple_edges_bunch(G: nx.MultiDiGraph):
    for u, v, k in G.edges:
        if k >= 1:
            return u, v
        if u in G[v]:
            return u, v
    return None

def gimme_junc_nodes_2deg(src_G):
    nodes_2deg = []
    for n in src_G.nodes:
        l1 = list(src_G.in_edges(n))
        l2 = list(src_G.out_edges(n))
        if len(l1) + len(l2) != 2:
            continue
        if not is_junction(src_G, n):
            continue
        nodes_2deg += [n]
    return nodes_2deg

def make_small_graph_for_mixer(solver_G: nx.DiGraph, removeRoot: bool = True) -> Tuple[nx.DiGraph, Dict]:
    new_G = nx.MultiDiGraph(solver_G)
    if removeRoot:
        if Root in new_G.nodes:
            new_G.remove_node(Root)
    fluid_em = dict()
    x_pos_em = dict()
    x_neg_em = dict()
    for u, v, k in new_G.edges:
        fluid_em[(u, v, k)] = [(u, v)] # key - edge in reduced graph, value - list of original graph edges
        x_pos_em[(u, v, k)] = [(u, v)]
        x_neg_em[(u, v, k)] = []
        assert k == 0

    nodes_2deg = gimme_junc_nodes_2deg(new_G)
    i, j = 0, 0
    while True:
        e = gimme_any_multiple_edges_bunch(new_G)
        if e is not None:
            i += 1
            # collapse bunch to one edge
            u, v = e
            k1 = len(new_G[u][v])
            k2 = 0
            if u in new_G[v]:
                k2 = len(new_G[v][u])
            bunch_uv = [(u, v, _k) for _k in range(k1)]
            bunch_vu = [(v, u, _k) for _k in range(k2)]
            new_G.remove_edges_from(bunch_uv + bunch_vu)
            fl = []
            pos, neg = [], []
            for e in bunch_uv:
                fl += fluid_em.pop(e)
                pos += x_pos_em.pop(e)
                neg += x_neg_em.pop(e)

            for e in bunch_vu:
                fl += fluid_em.pop(e)
                pos += x_neg_em.pop(e)
                neg += x_pos_em.pop(e)

            k = new_G.add_edge(u, v)
            fluid_em[(u, v, k)] = fl
            x_pos_em[(u, v, k)] = pos
            x_neg_em[(u, v, k)] = neg
            continue

        if len(nodes_2deg) == 0:
            break

        j += 1
        n = nodes_2deg.pop(0)
        in_edges = list(new_G.in_edges(n))
        out_edges = list(new_G.out_edges(n))
        edges = in_edges + out_edges
        if len(edges) !=2:
            continue
        e1 = edges[0]
        e2 = edges[1]
        u1, v1 = e1
        u2, v2 = e2

        wo_n = [u1, v1, u2, v2]
        wo_n.remove(n)
        wo_n.remove(n)
        if len(set(wo_n)) != 2: # two parallel edges
            continue
        u, v = tuple(wo_n)
        new_G.remove_edge(*e1)
        new_G.remove_edge(*e2)
        new_G.remove_node(n)
        k = new_G.add_edge(u, v)
        l1 = fluid_em.pop((u1, v1, 0))
        l2 = fluid_em.pop((u2, v2, 0))
        fluid_em[(u, v, k)] = l1+l2

        x_pos_em.pop((u2, v2, 0))
        x_neg_em.pop((u2, v2, 0))
        pos = x_pos_em.pop((u1, v1, 0))
        neg = x_neg_em.pop((u1, v1, 0))
        if n == u1:
            pos, neg = neg, pos
        x_pos_em[(u, v, k)] = pos
        x_neg_em[(u, v, k)] = neg
        continue


    rez_G = nx.DiGraph()
    rez_fluid_em = dict()
    rez_pos_em, rez_neg_em = dict(), dict()
    for u, v, k in new_G.edges:
        assert k == 0
        rez_G.add_edge(u, v)
        rez_fluid_em[(u, v)] = fluid_em[(u, v, 0)]
        rez_pos_em[(u, v)] = x_pos_em[(u, v, 0)]
        rez_neg_em[(u, v)] = x_neg_em[(u, v, 0)]


    return rez_G, dict(fluid_em=rez_fluid_em, x_pos_em=rez_pos_em, x_neg_em=rez_neg_em)

def evalute_network_fluids_wo_root(_G, x_dict):
    assert not (Root in _G.nodes)

    for u, v in _G.edges:
        assert x_dict[(u, v)] >= 0
    G = nx.DiGraph(_G)
    for (u, v) in x_dict:
        if abs(x_dict[(u, v)]) < 1e-7:
            G.remove_edge(u, v)
    for node in _G.nodes:
        if (len(G.in_edges(node)) == 0) and (len(G.out_edges(node)) == 0):
            G.remove_node(node)

    nodes = list(G.nodes)
    edges = list(G.edges)
    EN = len(edges)
    N = len(nodes)

    x = np.array([x_dict[e] for e in edges])
    A = -1 * nx.incidence_matrix(G, nodelist=nodes, edgelist=edges, oriented=True).toarray()
    Q = np.matmul(A, x)
    var_idx = {e:i for i, e in enumerate(edges)}
    sinks = [n for i, n in enumerate(nodes) if Q[i] < -1e-7]
    var_idx.update({n:i+EN for i, n in enumerate(sinks)})
    M = len(var_idx)
    mx1 = np.zeros((N, M)) # the first partition of matrix is for 1stCL
    mx2 = np.zeros((M-N, M)) # the second partition is to dictate condition: all fluids leaving one node have to be the same
    i2 = 0
    for i, node in enumerate(nodes):
        these_vars_are_the_same = []
        if node in sinks: # n is sink, so there is an unknown (equation variable) for outbound flow
            j = var_idx[node] # get the unknown index in matrix
            mx1[i, j] = Q[i]
            these_vars_are_the_same += [j] # remember the unknown for a while

        for u, v in G.in_edges(node):  # inlet flows
            j = var_idx[(u, v)]
            mx1[i, j] = x_dict[(u, v)]

        for u, v in G.out_edges(node): # all the outlet flows are the same
            j = var_idx[(u, v)]
            mx1[i, j] = -x_dict[(u, v)]
            these_vars_are_the_same += [j] # so we have to remember this unknown too

        if Q[i] > 0:
            mx1[i] /= -Q[i] # cause we need only ones and zeros on the right side

        if len(these_vars_are_the_same) > 1: # if there are more than one remembered outlet flow
            for fl1, fl2 in zip(these_vars_are_the_same[1:], these_vars_are_the_same[:-1]):
                mx2[i2, fl1] = 1  # we have to right an equation to matrix, to describe that these outlet fluids are equal
                mx2[i2, fl2] = -1
                i2 += 1
    mx = np.append(mx1, mx2, axis=0)
    assert mx.shape == (M, M)

    mx_inv = np.linalg.inv(mx)
    rez_mx = np.zeros((M, N))
    srcs = []
    S = 0
    for i, n in enumerate(nodes):
        if Q[i] > 1e-7:
            rez_mx[:, S] = mx_inv[:,i]
            srcs += [n]
            S += 1

    cocktails = {}
    for k, i in var_idx.items():
        cocktails[k] = rez_mx[i, :S]
        sck = cocktails[k].sum()
        if abs(sck - 1) > 1e-7:
            logger.error('Cocktail matrix is invalid')
            raise ValueError
    return cocktails, srcs
