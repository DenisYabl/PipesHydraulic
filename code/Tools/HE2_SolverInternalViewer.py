import numpy as np
# from Solver.HE2_Solver import HE2_Solver
import matplotlib.pyplot as plt
import networkx as nx
import pandas as pd
from Tools.HE2_tools import make_oilupstream_graph_layout
from Tools.HE2_ABC import Root

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


# def plot_residuals_toward_gradient(solver : HE2_Solver, start=-0.1, stop=1, points=111):
def plot_residuals_toward_gradient(solver, start=-0.1, stop=1, points=111, filter=0.1):
    x0 = solver.actual_x.copy()
    grad = solver.actual_dx
    steps = np.linspace(start, stop, points)
    ys = []
    edgelist = solver.edge_list
    n = len(solver.edge_list)
    traces = np.zeros((points, n))
    xs = np.zeros((points, n))

    stash_here_solver_itermediate_results_flag = solver.save_intermediate_results
    solver.save_intermediate_results = True

    for i, step in enumerate(steps):
        x = x0 + step * grad
        y = solver.target(x)
        ys += [y]
        for j, e in enumerate(edgelist):
            edge_func_result = solver.edge_func_last_results[e]
            xx, unknown, p_kn, p_unk = edge_func_result
            traces[i, j] = p_unk - p_kn
            xs[i, j] = xx

    fig = plt.figure(constrained_layout=True, figsize=(12, 8))
    ax = fig.add_subplot(1, 1, 1)
    hndls = []
    for i, e in enumerate(edgelist):
        trace = traces[:,i]
        x = xs[:,i]
        if max(trace) - min(trace) < filter:
            continue
        hndl, = ax.plot(x, trace, linewidth=1, label=str(e), alpha=0.2)
        hndls += [hndl]
        ax.scatter(x, trace)
    ax.legend(handles=hndls)
    plt.show()

    solver.save_intermediate_results = stash_here_solver_itermediate_results_flag


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
        for n in nodes:
            if not solver.pt_on_tree:
                break
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
            if not solver.edges_x:
                break
            x = solver.edges_x[e]
            sss = f'x:{x:.3f}'
            labels[e] += [sss]

    if 'dp/dx' in keys_to_plot:
        for e in edges:
            if not solver.derivatives:
                break
            dpdx = solver.derivatives[e]
            sss = f'dp/dx:{dpdx:.6f}'
            labels[e] += [sss]

    if 'dP' in keys_to_plot:
        for e in edges:
            if not solver.pt_on_tree:
                break
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


def make_chorde_cycle_graph(solver, u, v):
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
    return G, G_tree, G_chordes

def make_node_neighbours_graph(solver, nodes, deep = 1):
    solver_tree = solver.span_tree
    solver_chordes = solver.chordes
    solverG = solver.graph
    neighb_edges = []
    neighb_nodes = set(nodes)
    for i in range(deep):
        new_neigbs = []
        for u, v in solverG.edges:
            if (u, v) in neighb_edges:
                continue
            if u == Root or v == Root:
                continue
            if u in neighb_nodes:
                new_neigbs += [v]
                neighb_edges += [(u, v)]
            elif v in neighb_nodes:
                new_neigbs += [u]
                neighb_edges += [(u, v)]
        neighb_nodes |= set(new_neigbs)
        # neighb_nodes -= set([Root])

    neighb_edges = list(set(neighb_edges))
    G = nx.DiGraph()
    G_tree = []
    G_chordes = []
    for (_u, _v) in neighb_edges:
        if (_u, _v) in solver_tree:
            G.add_edge(_u, _v)
            G_tree += [(_u, _v)]
        if  (_u, _v) in solver_chordes:
            G.add_edge(_u, _v)
            G_chordes += [(_u, _v)]
    return G, G_tree, G_chordes

def plot_some_subgraph(G, G_tree, G_chordes, solver, keys_to_plot):
    fig = plt.figure(constrained_layout=True, figsize=(8, 8))
    ax = fig.add_subplot(1, 1, 1)
    pos = nx.drawing.planar_layout(G)
    # pos = nx.drawing.layout.circular_layout(G)
    # pos = make_oilupstream_graph_layout(G)
    nx.draw_networkx_nodes(G, node_size=15, node_color='black', ax=ax,pos=pos, alpha=0.9)
    nx.draw_networkx_edges(G, edgelist=G_tree, pos=pos, width=1, ax=ax, edge_color='black')
    nx.draw_networkx_edges(G, edgelist=G_chordes, pos=pos, width=1, ax=ax, edge_color='black', style='dashed')

    n_lbls, n_lbls_pos = make_node_labels(G, solver, pos, keys_to_plot)
    nx.draw_networkx_labels(G, labels=n_lbls, ax=ax, pos=n_lbls_pos, font_size=7)

    e_lbls = make_edge_labels(G, solver, keys_to_plot)
    nx.draw_networkx_edge_labels(G, pos=pos, edge_labels=e_lbls, font_size=7)

    plt.show()


def plot_chord_cycle(solver, u, v, keys_to_plot=('name', 'P', 'Q', 'x', 'dp/dx', 'dP', 'WC', 'GOR', 'pos')):
    G, G_tree, G_chordes = make_chorde_cycle_graph(solver, u, v)
    plot_some_subgraph(G, G_tree, G_chordes, solver, keys_to_plot)

def plot_neighbours_subgraph(solver, nodes, deep = 1, keys_to_plot=('name', 'P', 'Q', 'x', 'dp/dx', 'dP', 'WC', 'GOR', 'pos')):
    G, G_tree, G_chordes = make_node_neighbours_graph(solver, nodes, deep)
    plot_some_subgraph(G, G_tree, G_chordes, solver, keys_to_plot)

def plot_all(solver, keys_to_plot=('name', 'P', 'Q', 'x', 'dp/dx', 'dP', 'WC', 'GOR', 'pos')):
    G = solver.graph
    G_tree = solver.span_tree
    G_chordes = solver.chordes
    plot_some_subgraph(G, G_tree, G_chordes, solver, keys_to_plot)

def plot_all_wo_root(solver, keys_to_plot=('name', 'P', 'Q', 'x', 'dp/dx', 'dP', 'WC', 'GOR', 'pos')):
    nodelist = solver.node_list[::]
    nodelist.remove(Root)
    G = nx.DiGraph()
    G.add_nodes_from(nodelist)
    edgelist = []
    for u, v in solver.edge_list:
        if u == Root or v == Root:
            continue
        edgelist += [(u, v)]
    G.add_edges_from(edgelist)
    G_tree = list(set(solver.span_tree) & set(edgelist))
    G_chordes = list(set(solver.chordes) & set(edgelist))
    plot_some_subgraph(G, G_tree, G_chordes, solver, keys_to_plot)

