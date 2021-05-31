import networkx as nx
import numpy as np
from Tools.HE2_ABC import Root
from Tools.HE2_Logger import getLogger
from Tools.HE2_SolverInternalViewer import plot_neighbours_subgraph as plot_nghbs

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
        sck = sum(cocktails[k])
        if abs(sck - 1) > 1e-7:
            logger.error('Cocktail matrix is invalid')
            raise ValueError
    return cocktails, srcs


def evalute_network_fluids_with_root(G, x_dict):
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
    G3, x_dict3, edges_mapping, nodes_mapping = reduce_graph(G2, x_dict2)

    cocktails3, srcs = evalute_network_fluids_wo_root(G3, x_dict3)
    cocktails2, srcs = restore_cocktail_from_reduced_graph(cocktails3, srcs, edges_mapping, nodes_mapping, G2)

    cocktails = {}
    for key, cktl in cocktails2.items():
        cktl2 = np.around(cktl, 6)

        if (key in edges2) and not (key in edges1):
            u, v = key
            cocktails[(v, u)] = cktl2
        else:
            cocktails[key] = cktl2
    return cocktails, srcs

def reduce_graph2(G, x_dict):
    e_out = dict()
    e_in = dict()
    for u, v in G.edges():
        e_out[u] = [v]
        e_in[v] = [u]

    nodes = list(set(e_out.keys) | set(e_in.keys))
    while True:
        for n in nodes:
            if len(e_out[n]) == len(e_in[n]) == 1:
                found = True
                u = e_out[n][0]
                v = e_in[n]
                e_out.pop(n)
                e_in.pop(n)
                e_out[u] += [v]
                e_in[v] += [u]
    pass


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

def restore_cocktail_from_reduced_graph(cocktails_reduced, srcs, edges_mapping, nodes_mapping, original_G):
    cocktails = {}
    for u, v in original_G.edges:
        u, v, k = edges_mapping[(u, v)]
        if not (u, v, k) in cocktails_reduced:
            continue
        cocktails[(u, v)] = cocktails_reduced[(u, v, k)]

    for n in original_G.nodes:
        n2 = nodes_mapping[n]
        if not n in cocktails_reduced:
            continue
        cocktails[n] = cocktails_reduced[n2]

    return cocktails, srcs
