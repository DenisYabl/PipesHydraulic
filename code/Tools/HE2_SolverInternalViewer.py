import numpy as np
# from Solver.HE2_Solver import HE2_Solver
import matplotlib.pyplot as plt
import networkx as nx

def plot_y_toward_gradient_from_actual_x(solver, start=-0.1, stop=1, points=111):
    x0 = solver.actual_x.copy()
    grad = solver.actual_dx
    steps = np.linspace(start, stop, points)
    ys = []

    for step in steps:
        x = x0 + step * grad
        y = solver.target(x)
        ys += [y]

    fig = plt.figure(constrained_layout=True, figsize=(12, 8))
    ax = fig.add_subplot(1, 1, 1)
    plt.scatter(steps, ys)
    plt.show()

def make_node_labels(G, solver, pos, keys_to_plot):
# def make_node_labels(G, solver : HE2_Solver, pos, keys_to_plot):
    nodes = set(G.nodes) & set(pos.keys())
    labels = {n:[] for n in nodes}
    if 'name' in keys_to_plot:
        for n in nodes:
            labels[n] += [n]

    if 'pos' in keys_to_plot:
        for n in nodes:
            x, y = pos[n][0], pos[n][1]
            sss = f'xy:{x:.3f} {y:.3f}'
            labels[n] += [sss]

    if 'P' in keys_to_plot:
        chord_ends = {v: (u, v) for (u, v) in solver.chordes}
        print(chord_ends)
        for n in nodes:
            pt = solver.pt_on_tree[n]
            sss = f'P:{pt[0]:.3f}'
            if n in chord_ends:
                chord = chord_ends[n]
                pt2 = solver.pt_on_chords_ends[chord]
                sss += f' P_с:{pt2[0]:.3f}'
            labels[n] += [sss]

    if 'Q' in keys_to_plot:
        for n in nodes:
            if n in solver.known_Q:
                Q = solver.known_Q[n]
                sss = f'Q:{Q:.3f}'
                labels[n] += [sss]

    rez_lbls = dict()
    for n in nodes:
        rez_lbls[n] = '\n'.join(labels[n])

    labels_pos = dict()
    for n in nodes:
        node_lbls = labels[n]
        ls = [len(s) for s in node_lbls]
        H = len(node_lbls)
        W = max(ls)
        xy_arr = pos[n]
        dxy_arr = np.zeros(xy_arr.shape)
        dxy_arr[0] = W * 0.01 + 0.02
        dxy_arr[1] = H * 0.025
        # labels_pos[n] = xy_arr + dxy_arr
        labels_pos[n] = xy_arr * 0.84

    return rez_lbls, labels_pos


def make_edge_labels(G, solver, keys_to_plot):
# def make_edge_labels(G, solver: HE2_Solver, keys_to_plot):
    edges = list(G.edges)
    labels = {e:[] for e in edges}
    if 'x' in keys_to_plot:
        for e in edges:
            x = solver.edges_x[e]
            sss = f'x:{x:.3f}'
            labels[e] += [sss]

    if 'dp/dx' in keys_to_plot:
        for e in edges:
            dpdx = solver.derivatives[e]
            sss = f'dp/dx:{dpdx:.6f}'
            labels[e] += [sss]

    if 'dP' in keys_to_plot:
        for e in edges:
            u, v = e
            P_u = solver.pt_on_tree[u][0]
            P_v = solver.pt_on_tree[v][0]
            dp = P_v - P_u
            sss = f'dP:{dp:.3f}'
            labels[e] += [sss]


    rez_lbls = dict()
    for e in edges:
        rez_lbls[e] = '\n'.join(labels[e])
    return rez_lbls


def plot_chord_cycle(solver, u, v, keys_to_plot=('name', 'P', 'Q', 'x', 'dp/dx', 'dP', 'WC', 'GOR', 'pos')):
    solver_tree = solver.span_tree
    solver_chordes = solver.chordes
    G = nx.Graph()
    G.add_edges_from(solver_tree)
    # G.add_edges_from(solver_chordes)
    path_nodelist = nx.shortest_path(G, u, v)
    path_edgelist = list(zip(path_nodelist[:-1], path_nodelist[1:]))
    path = path_edgelist
    path += [(u, v)]
    G = nx.DiGraph()
    G_tree = []
    G_chordes = []
    for (_u, _v) in path:
        if (_u, _v) in solver_tree:
            G.add_edge(_u, _v)
            G_tree += [(_u, _v)]
        if  (_u, _v) in solver_chordes:
            G.add_edge(_u, _v)
            G_chordes += [(_u, _v)]
        if (_v, _u) in solver_tree:
            G.add_edge(_v, _u)
            G_tree += [(_v, _u)]
        if (_v, _u) in solver_chordes:
            G.add_edge(_v, _u)
            G_chordes += [(_v, _u)]

    fig = plt.figure(constrained_layout=True, figsize=(8, 8))
    ax = fig.add_subplot(1, 1, 1)
    pos = nx.drawing.layout.circular_layout(G)
    nx.draw_networkx_nodes(G, node_size=15, node_color='black', ax=ax,pos=pos, alpha=0.9)
    nx.draw_networkx_edges(G, edgelist=G_tree, pos=pos, width=1, ax=ax, edge_color='black')
    nx.draw_networkx_edges(G, edgelist=G_chordes, pos=pos, width=1, ax=ax, edge_color='black', style='dashed')

    n_lbls, n_lbls_pos = make_node_labels(G, solver, pos, keys_to_plot)
    nx.draw_networkx_labels(G, labels=n_lbls, ax=ax, pos=n_lbls_pos, font_size=7)

    e_lbls = make_edge_labels(G, solver, keys_to_plot)
    nx.draw_networkx_edge_labels(G, pos=pos, edge_labels=e_lbls, font_size=7)

    plt.show()

