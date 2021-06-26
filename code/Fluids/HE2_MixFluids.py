import networkx as nx
import numpy as np
from Tools.HE2_ABC import Root
from Tools.HE2_Logger import getLogger
from Tools.HE2_SolverInternalViewer import plot_neighbours_subgraph as plot_nghbs
from GraphNodes.HE2_Vertices import is_junction
from typing import Tuple, Dict

logger = getLogger(__name__)


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
        sck = np.sum(cocktails[k])
        if abs(sck - 1) > 1e-7:
            logger.error('Cocktail matrix is invalid')
            raise ValueError
    return cocktails, srcs


def evalute_network_fluids_with_root(G, x_dict):
# def evalute_network_fluids_with_root_old(G, x_dict):
    edges1 = [(u, v) for u, v in G.edges if (u != Root) and (v != Root)]
    edges2 = []
    x_dict2 = {}
    for (u, v) in edges1:
        e = (u, v)
        if x_dict[(u, v)] < 0:
            e = (v, u)
        x_dict2[e] = abs(x_dict[(u, v)])
        edges2 += [e]

    G2 = nx.DiGraph(edges2)
    cocktails, srcs = evalute_network_fluids_wo_root(G2, x_dict2)
    cocktails2 = {}
    for key, cktl in cocktails.items():
        cktl2 = np.around(cktl, 6)

        if (key in edges2) and not (key in edges1):
            u, v = key
            cocktails2[(v, u)] = cktl2
        else:
            cocktails2[key] = cktl2
    return cocktails2, srcs

def evalute_network_fluids_with_root_new_version_with_reduce(G, x_dict):
# def evalute_network_fluids_with_root(G, x_dict):
    # This is like a nested doll
    # Outer layer - solver call with solver graph and xs (flows)
    # Second layer - we remove root node and zero flows, also reverse negative flows
    # Third layer - we reduce graph removing 2-degree nodes and parallel edges, cause it does not affect to solultion
    # Internal layer - on reduced graph we build mass conservation equation matrix and solve it
    # Than, from deep to surface, we rebuild solution on each layer

    edges1 = [(u, v) for u, v in G.edges if (u != Root) and (v != Root)]
    edges2 = []
    x_dict2 = {}
    for (u, v) in edges1:
        e = (u, v)
        if x_dict[(u, v)] < 0:
            e = (v, u)
        x_dict2[e] = abs(x_dict[(u, v)])
        edges2 += [e]

    G2 = nx.DiGraph(edges2)
    G3, x_dict3, reducing_stack = reduce_graph2(G2, x_dict2)

    cocktails3, srcs = evalute_network_fluids_wo_root(G3, x_dict3)
    cocktails2, srcs = restore_cocktail_from_reduced_graph(cocktails3, srcs, reducing_stack)

    cocktails = {}
    for key, cktl in cocktails2.items():
        cktl2 = np.around(cktl, 6)

        if (key in edges2) and not (key in edges1):
            u, v = key
            cocktails[(v, u)] = cktl2
        else:
            cocktails[key] = cktl2
    return cocktails, srcs

# def evalute_network_fluids_with_root(G, x_dict):
def evalute_network_fluids_with_root_third_version_with_reduce(G, x_dict):
    # This is like a nested doll
    # Outer layer - solver call with solver graph and xs (flows)
    # Second layer - we remove root node and zero flows, also reverse negative flows
    # Third layer - we reduce graph removing 2-degree nodes and parallel edges, cause it does not affect to solultion
    # Internal layer - on reduced graph we build mass conservation equation matrix and solve it
    # Than, from deep to surface, we rebuild solution on each layer

    edges1 = [(u, v) for u, v in G.edges if (u != Root) and (v != Root)]
    edges2 = []
    x_dict2 = {}
    for (u, v) in edges1:
        e = (u, v)
        if x_dict[(u, v)] < 0:
            e = (v, u)
        x_dict2[e] = abs(x_dict[(u, v)])
        edges2 += [e]

    G2 = nx.DiGraph(edges2)
    G3, x_dict3, reducing_stack = reduce_graph2(G2, x_dict2)

    cocktails3, srcs = evalute_network_fluids_wo_root(G3, x_dict3)
    cocktails2, srcs = restore_cocktail_from_reduced_graph(cocktails3, srcs, reducing_stack)

    cocktails = {}
    for key, cktl in cocktails2.items():
        cktl2 = np.around(cktl, 6)

        if (key in edges2) and not (key in edges1):
            u, v = key
            cocktails[(v, u)] = cktl2
        else:
            cocktails[key] = cktl2
    return cocktails, srcs


NodeMapping = Dict[str, str]
EdgeMapping = Dict[Tuple[str, str], Tuple[str, str]]

def make_small_graph_for_mixer(src_G: nx.DiGraph) -> Tuple[nx.DiGraph, NodeMapping, EdgeMapping]:
    new_G = nx.MultiDiGraph(src_G)
    nodes_2deg = gimme_junc_nodes_2deg(src_G)
    node_mapping = dict()
    edge_mapping = dict()
    pos_em, neg_em = dict(), dict() # Key - edge in reduced graph. Value - list of edges of original graph
    for u, v in src_G.edges:
        pos_em[(u, v, 0)] = [(u, v)]
        neg_em[(u, v, 0)] = []

    i, j = 0, 0
    while True:
        if len(nodes_2deg) > 0:
            i += 1
            n = nodes_2deg.pop(0)
            in_edges = list(new_G.in_edges(n))
            assert len(in_edges) == 1
            e1 = in_edges[0]
            out_edges = list(new_G.out_edges(n))
            assert len(out_edges) == 1
            e2 = out_edges[0]
            u1, v1 = e1
            u2, v2 = e2

            assert u2 == n and v1 == n
            u, v = u1, v2
            new_G.remove_edge(*e1)
            new_G.remove_edge(*e2)
            new_G.remove_node(n)
            k = new_G.add_edge(u, v)
            pel = pos_em.pop((u1, v1, 0))
            pel += pos_em.pop((u2, v2, 0))
            nel = neg_em.pop((u1, v1, 0))
            nel += neg_em.pop((u2, v2, 0))
            pos_em[(u, v, k)] = pel
            neg_em[(u, v, k)] = nel


        e = gimme_any_multiple_edges_bunch(new_G)
        if e is None:
            break

        j += 1
        u, v = e
        k1 = len(new_G[u][v])
        k2 = 0
        if v in new_G[v]:
            k2 = len(new_G[v][u])
        bunch_uv = [(u, v, _k) for _k in range(k1)]
        bunch_vu = [(v, u, _k) for _k in range(k2)]
        new_G.remove_edges_from(bunch_uv + bunch_vu)
        pel, nel = [], []
        for key in bunch_uv:
            pel += pos_em.pop(key)
            nel += neg_em.pop(key)

        for key in bunch_vu:
            pel += neg_em.pop(key) #sic!
            nel += pos_em.pop(key) #sic!

        new_G.add_edge(u, v)
        edge_to_add = u, v, 0
        pos_em[edge_to_add] = pel
        neg_em[edge_to_add] = nel

    rez_G = nx.DiGraph()
    for u, v, k in new_G.edges:
        assert k == 0
        rez_G.add_edge(u, v)

    return rez_G, node_mapping, edge_mapping



def gimme_any_multiple_edges_bunch(G: nx.MultiDiGraph):
    for u in G:
        for v in G[u]:
            if len(G[u][v]) > 1:
                return u, v
    return None

def reduce_graph2(src_G, x_dict):
    new_G = nx.MultiDiGraph(src_G)
    new_x_dict = {(u, v, 0): x_dict[(u, v)] for (u, v) in x_dict}
    reducing_stack = []

    nodes_2deg = gimme_junc_nodes_2deg(src_G)

    i, j = 0, 0
    while True:
        if len(nodes_2deg) > 0:
            i += 1
            n = nodes_2deg.pop(0)
            e1 = list(new_G.in_edges(n))[0]
            e2 = list(new_G.out_edges(n))[0]
            u1, v1 = e1
            u2, v2 = e2

            x1 = new_x_dict.pop((u1, v1, 0))
            x2 = new_x_dict.pop((u2, v2, 0))
            if abs(x1-x2) > 1e-6:
                continue

            assert u2 == n and v1 == n
            u, v = u1, v2
            new_G.remove_edge(*e1)
            new_G.remove_edge(*e2)
            new_G.remove_node(n)
            k = new_G.add_edge(u, v)
            new_x_dict[(u, v, k)] = x1
            reducing_stack += [dict(op='2node', rem_node=n, rem_e1=(u1, v1, 0), rem_e2=(u2, v2, 0), add_e=(u,v,k))]
            continue

        e = gimme_any_multiple_edges_bunch(new_G)
        if e is None:
            break

        j += 1
        u, v = e
        k1 = len(new_G[u][v])
        k2 = 0
        if v in new_G[v]:
            k2 = len(new_G[v][u])
        bunch_uv = [(u, v, _k) for _k in range(k1)]
        bunch_vu = [(v, u, _k) for _k in range(k2)]
        new_G.remove_edges_from(bunch_uv + bunch_vu)

        x = 0
        for e in bunch_uv:
            x += new_x_dict.pop(e)
        for e in bunch_vu:
            x -= new_x_dict.pop(e)

        edge_to_add = None
        if x > 1e-6:
            new_G.add_edge(u, v)
            edge_to_add = u, v, 0
            new_x_dict[edge_to_add] = x
        elif x < -1e-6:
            new_G.add_edge(v, u)
            edge_to_add = v, u, 0
            new_x_dict[edge_to_add] = x

        reducing_stack += [dict(op='edges', u=u, v=v, uv_k=k1, vu_k=k2, add_e=edge_to_add)]

    rez_G = nx.DiGraph()
    rez_x_dict = dict()
    for u, v, k in new_G.edges:
        assert k == 0
        rez_G.add_edge(u, v)
        rez_x_dict[(u, v)] = new_x_dict[(u, v, 0)]

    return rez_G, rez_x_dict, reducing_stack


def gimme_junc_nodes_2deg(src_G):
    nodes_2deg = []
    for n in src_G.nodes:
        l1 = list(src_G.in_edges(n))
        l2 = list(src_G.out_edges(n))
        if len(l1) != 1 or len(l2) != 1:
            continue
        # if not is_junction(src_G, n):
        #     continue
        nodes_2deg += [n]
    return nodes_2deg


def reduce_graph(G, x_dict):
    new_G = nx.MultiDiGraph()
    edges = list(G.edges)
    new_G.add_edges_from(edges)
    new_x_dict = x_dict.copy()
    edges_mapping = {(u, v):(u, v, 0) for (u, v) in G.edges}
    nodes_mapping = {n:n for n in G.nodes}
    while True:
        found = False
        for n, d in new_G.degree():
            if d!=2:
                continue
            if not (len(new_G.in_edges(n))==1 and len(new_G.out_edges(n))==1):
                continue
            e1 = list(new_G.in_edges(n))[0]
            e2 = list(new_G.out_edges(n))[0]
            found = True
            u1, v1 = e1
            u2, v2 = e2
            assert u2 == n and v1 == n
            u, v = u1, v2

            new_G.remove_edge(*e1)
            new_G.remove_edge(*e2)
            k = G.add_edge(u, v)
            for _k, _v in edges_mapping.items():
                if _v == e1 or _v == e2:
                    edges_mapping[_k] = u, v, 0
            for _k, _v in nodes_mapping.items():
                if _v == n:
                    nodes_mapping[_k] = u
            new_x_dict[(u, v, k)] = x_dict[e1]
            new_x_dict.pop(e1)
            new_x_dict.pop(e2)
            break

        if found:
            continue

        for u, v, k in new_G.edges:
            if k == 0:
                continue
            found = True
            new_G.remove_edge(u, v, k)
            for _k, _v in edges_mapping.items():
                if _v == (u, v, k) or _v == (u, v, k):
                    edges_mapping[_k] = u, v, 0
            new_x_dict[(u, v, 0)] += new_x_dict[(u, v, k)]

        if not found:
            break

    rez_G = nx.DiGraph(new_G)
    rez_x_dict = {(u, v):new_x_dict[(u, v, k)] for u, v, k in new_x_dict}
    return rez_G, rez_x_dict, edges_mapping, nodes_mapping

def restore_cocktail_from_reduced_graph(cocktails_reduced, srcs, reducing_stack):
    cocktails = {} # cocktails_reduced.copy()
    for key in cocktails_reduced:
        if len(key) == 2:
            u, v = key
            cocktails[u, v, 0] = cocktails_reduced[key]
    while reducing_stack:
        item = reducing_stack.pop(-1)

        # reducing_stack += [dict(op='2node', rem_node=n, rem_e1=e1, rem_e2=e2, add_e=(u,v,k))]
        if item['op'] == '2node':
            try:
                cktl = cocktails.pop(item['add_e'])
            except:
                pass
            cocktails[item['rem_node']] = cktl
            cocktails[item['rem_e1']] = cktl
            cocktails[item['rem_e2']] = cktl
            continue

        assert item['op'] == 'edges'
        # reducing_stack += [dict(op='edges', u=u, v=v, uv_k=k1, vu_k=k2, added_edge=edge_to_add)]
        add_e = item['add_e']
        cktl = cocktails.pop(add_e)
        u = item['u']
        v = item['v']
        k1 = item['uv_k']
        k2 = item['vu_k']
        for k in range(k1):
            cocktails[(u, v, k)] = cktl
        for k in range(k2):
            cocktails[(v, u, k)] = cktl

    rez_cocktails = {}
    for key in cocktails:
        if len(key) <= 2:
            rez_cocktails[key] = cocktails[key]
        if len(key) == 3:
            u, v, k = key
            assert k == 0
            rez_cocktails[(u, v)] = cocktails[key]

    return rez_cocktails, srcs
