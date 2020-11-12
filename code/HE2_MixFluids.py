import networkx as nx
import numpy as np
from HE2_ABC import Root
import HE2_tools as tools
from HE2_Solver import HE2_Solver


def evalute_network_fluids_wo_root(_G, x_dict):
    for u, v in _G.edges:
        assert x_dict[(u, v)] >= 0
    G = nx.DiGraph(_G)
    for (u, v) in x_dict:
        if x_dict[(u, v)] == 0:
            G.remove_edge(u, v)
    for n in _G.nodes:
        if len(G[n])==0:
            G.remove_node(n)

    assert not (Root in G.nodes)

    nodes = list(G.nodes)
    edges = list(G.edges)
    x = np.array([x_dict[e] for e in edges])
    A = -1 * nx.incidence_matrix(G, nodelist=nodes, edgelist=edges, oriented=True).toarray()
    Q = np.matmul(A, x)
    variables = [i for i, e in enumerate(edges) if x_dict[e] != 0]
    EN = len(variables)
    variables += [i for i, n in enumerate(nodes) if Q[i] < 0]
    N = len(variables)
    QN = N - EN
    mx = np.zeros((N, N))
    for i, n in enumerate(nodes):
        if Q[i] < 0: # n is sink
            j = get_variable_idx(n)
            mx[i, j] = Q[i]
        for u, v in G.in_edges(n):
            j = get_variable_idx(u, v)
            mx[i, j] = x_dict[(u, v)]
        for u, v in G.out_edges(n):
            j = get_variable_idx(u, v)
            mx[i, j] = -x_dict[(u, v)]


if __name__ == '__main__':
    G0, nodes = tools.generate_random_net_v0(N=7, E=9, SNK=2, randseed=424242)
    solver = HE2_Solver(G0)
    solver.solve()
    x_dict = dict()
    G = nx.DiGraph()
    G.add_nodes_from(G0)
    for u, v in G0.edges:
        x = G0[u][v]['obj'].result['x']
        if x < 0:
            x_dict[(v,u)] = -x
            G.add_edge(v, u)
        else:
            x_dict[(u, v)] = x
            G.add_edge(u, v)
    evalute_network_fluids_wo_root(G, x_dict)

