import numpy as np
import networkx as nx
import scipy.optimize as scop

from GraphEdges.HE2_SpecialEdges import HE2_MockEdge
from GraphNodes import HE2_Vertices as vrtxs
from Tools import HE2_ABC as abc
from Tools.HE2_ABC import Root
import pandas as pd
from Tools.HE2_Logger import check_for_nan, getLogger
import Fluids.HE2_MixFluids as mixer
import Fluids.HE2_Fluid as fl
from Tools.HE2_SolverInternalViewer import plot_y_toward_gradient_from_actual_x as plot_y, plot_chord_cycle as plot_chord
from Tools.HE2_SolverInternalViewer import plot_neighbours_subgraph as plot_nghbs


logger = getLogger(__name__)

class HE2_Solver():
    def __init__(self, schema):
        logger.debug('New solver instance created')
        self.schema = schema
        self.graph = None
        self.op_result = None
        self.span_tree = None
        self.chordes = None
        self.edge_list = []
        self.node_list = []
        self.A_tree = None
        self.A_chordes = None
        self.A_inv = None
        self.Q_static = None
        self.edges_x = None
        self.pt_on_tree = None
        self.tree_travers = None
        self.mock_nodes = []
        self.mock_edges = []
        self.result_edges_mapping = dict()
        self.pt_on_chords_ends = None
        self.pt_residual_vec = None

        self.imd_rez_df = None
        self.save_intermediate_results = False
        self.initial_edges_x = None
        self.ready_for_solve = False

        self.derivatives = None
        self.forward_edge_functions = dict()
        self.backward_edge_functions = dict()

        self.fluids_move_rate = 0.5
        self.known_Q = dict()
        self.actual_x = None
        self.actual_dx = None

    def set_known_Q(self, Q_dict):
        self.known_Q = Q_dict


    def fill_known_Q(self):
        G = self.graph
        Q_dict = self.known_Q
        nodelist = [n for n in G.nodes()]
        known_src, unknown_src = dict(), dict()
        known_snk, unknown_snk = dict(), dict()
        for n in nodelist:
            obj = G.nodes[n]['obj']
            if isinstance(obj, vrtxs.HE2_Source_Vertex):
                d_kn, d_unk = known_src, unknown_src
            elif isinstance(obj, vrtxs.HE2_Boundary_Vertex):
                d_kn, d_unk = known_snk, unknown_snk
            else:
                continue

            if obj.Q:
                d_kn[n] = abs(obj.Q)
            elif n in Q_dict:
                d_kn[n] = abs(Q_dict[n])
            else:
                d_unk[n] = None

        if len(known_snk) == 0 and len(known_src) == 0:
            for n in unknown_src:
                known_src[n] = 1
            unknown_src = dict()

        known_src_sum = sum(known_src.values())
        known_snk_sum = sum(known_snk.values())
        if len(unknown_src) > 0 and len(unknown_snk) > 0: # Известны не все сорцы, и не все синки
            if known_src_sum > known_snk_sum: # В этом случае делаем все сорцы известными, дозаполняя средними значениями
                avg_src_Q = known_src_sum / len(known_src)
                for n in unknown_src:
                    known_src[n] = avg_src_Q
                unknown_src = dict()
            else: # тут делаем все синки известными, дозаполняя средними значениями
                avg_snk_Q = known_snk_sum / len(known_snk)
                for n in unknown_snk:
                    known_snk[n] = avg_snk_Q
                unknown_snk = dict()

        if len(unknown_src) == 0 and len(unknown_snk) > 0: # Известны все сорцы, надо дополнить все синки так чтобы сбилась сумма
            known_delta_Q = sum(known_src.values()) - sum(known_snk.values())
            avg_Q = known_delta_Q / len(unknown_snk)
            for n in unknown_snk:
                known_snk[n] = avg_Q
            unknown_snk = dict()
        elif len(unknown_snk) == 0 and len(unknown_src) > 0:
            known_delta_Q = sum(known_snk.values()) - sum(known_src.values())
            avg_Q = known_delta_Q / len(unknown_src)
            for n in unknown_src:
                known_src[n] = avg_Q
            unknown_src = dict()
        else:
            pass

        assert len(unknown_src) == 0
        assert len(unknown_snk) == 0

        known_src_sum = sum(known_src.values())
        known_snk_sum = sum(known_snk.values())
        if abs(known_src_sum - known_snk_sum) > 1e-3: # Такое может случится если снаружи пришли кривые данные. Выправим, чтобы суммы совпадали
            if known_src_sum > known_snk_sum:
                k = known_src_sum / known_snk_sum
                for n in known_snk:
                    known_snk[n] = known_snk[n] * k
            elif known_src_sum < known_snk_sum:
                k = known_snk_sum / known_src_sum
                for n in known_src:
                    known_src[n] = known_src[n] * k

        for n in known_snk:
            known_snk[n] = -1 * known_snk[n]

        self.known_Q = dict()
        self.known_Q.update(known_src)
        self.known_Q.update(known_snk)


    def make_initial_approximation(self):
        G = self.graph
        self.fill_known_Q()
        Q_dict = self.known_Q
        nodelist = [n for n in G.nodes()]
        assert not (Root in nodelist), 'better call this method before Root add'
        edgelist = [(u, v) for (u, v) in G.edges()]
        A_full = nx.incidence_matrix(G, nodelist=nodelist, edgelist=edgelist, oriented=True)
        A_full = -1 * A_full.toarray()
        q_vec = np.zeros((len(nodelist), 1))
        for i, node in enumerate(nodelist):
            if node in Q_dict:
                q_vec[i] = Q_dict[node]

        logger.debug(f'q_vec = {q_vec.flatten()}')
        A_inv = np.linalg.pinv(A_full)
        xs = np.matmul(A_inv, q_vec)
        self.initial_edges_x = dict(zip(edgelist, xs.flatten()))

        # cocktails, srcs = mixer.evalute_network_fluids_with_root(G, self.initial_edges_x)
        # src_fluids = [G.nodes[n]['obj'].fluid for n in srcs]
        # for key, cktl in cocktails.items():
        #     initial_fluid = fl.dot_product(list(zip(cktl, src_fluids)))
        #     if key in edgelist:
        #         u, v = key
        #         obj = G[u][v]['obj']
        #     else:
        #         obj = G.nodes[key]['obj']
        #     obj.fluid = initial_fluid

        logger.debug(f'initial_edges_x = {self.initial_edges_x}')

    def get_initial_approximation(self):
        x0 = np.zeros((len(self.chordes), 1))
        if self.initial_edges_x is None:
            logger.info('is finished. Graph is a tree')
            return x0

        for i, c in enumerate(self.chordes):
            x0[i] = self.initial_edges_x[c]

        logger.debug(f'x0 = {x0}')
        return x0

    def prepare_for_solve(self):
        logger.debug('is started')
        self.graph = self.transform_multi_di_graph_to_equal_di_graph(self.schema)
        self.make_initial_approximation()

        self.graph = self.add_root_to_graph(self.graph)
        self.span_tree, self.chordes = self.split_graph(self.graph)
        self.edge_list = self.span_tree + self.chordes
        self.node_list = list(self.graph.nodes())
        if self.node_list[-1] != Root:
            logger.error(f'Something wrong with graph restructure, Root shoold be last node in node_list')
            assert False
        self.tree_travers = self.build_tree_travers(self.span_tree, Root)
        self.A_tree, self.A_chordes = self.build_incidence_matrices()
        if self.A_tree.shape != (len(self.node_list)-1, len(self.node_list)-1):
            logger.error(f'Invalid spanning tree, inc.matrix shape is {self.A_tree.shape}, check graph structure.')
            assert False
        self.A_inv = np.linalg.inv(self.A_tree)
        self.B = self.build_circuit_matrix()
        self.Bt = np.transpose(self.B)
        self.Q_static = self.build_static_Q_vec(self.graph)
        for (u, v) in self.edge_list:
            obj = self.graph[u][v]['obj']
            if not isinstance(obj, abc.HE2_ABC_GraphEdge):
                assert False
            self.forward_edge_functions[(u, v)] = obj.perform_calc_forward
            self.backward_edge_functions[(u, v)] = obj.perform_calc_backward

        self.ready_for_solve = True

    def target(self, x_chordes):
        check_for_nan(x_chordes=x_chordes)

        Q = self.Q_static
        check_for_nan(Q_static=Q)

        x = x_chordes.reshape((len(x_chordes), 1))
        Q_dynamic = np.matmul(self.A_chordes, x)
        Q = Q - Q_dynamic
        check_for_nan(Q_dynamic=Q)

        x_tree = np.matmul(self.A_inv, Q)
        self.edges_x = dict(zip(self.span_tree, x_tree.flatten()))
        self.edges_x.update(dict(zip(self.chordes, x_chordes.flatten())))

        self.pt_on_tree = self.evalute_pressures_by_tree()
        self.pt_residual_vec, self.pt_on_chords_ends = self.evalute_chordes_pressure_residual()
        check_for_nan(chordes_pt_residual_vec=self.pt_residual_vec)

        rez = np.linalg.norm(self.pt_residual_vec)
        if self.save_intermediate_results:
            self.do_save_intermediate_results()

        return rez

    def solve(self, save_intermediate_results=False, threshold=0.05, it_limit = 100, step = 1):
        logger.info('is started')
        y_best, x_best, it_num, rnd_seed = 100500100500, None, 0, 42
        try:
            if not self.ready_for_solve:
                self.prepare_for_solve()

            x_chordes = self.get_initial_approximation()
            dx = np.zeros(x_chordes.shape)

            while True:
                it_num += 1
                self.actual_x = x_chordes
                self.actual_dx = dx
                # Best place to call plot_y(self) in debugger console
                # or maybe plot_chord(self, node1, node2)
                # plot_chord(self, 'Root', 'PAD_33')

                x_chordes = x_chordes + step * dx
                y = self.target(x_chordes)
                logger.debug(f'X = {x_chordes.flatten()}')
                logger.info(f'Y = {y}')

                # if y < y_best:
                #     self.evaluate_and_set_new_fluids()

                logger.info(f'it_num = {it_num}, y = {y}, step = {step}')
                if y < y_best:
                    logger.info(f'y {y} is better than y_best {y_best}')
                    y_best = y
                    x_best = x_chordes
                elif step > 0.1:
                    step = step/2
                else:
                    step = self.gimme_random_step(rnd_seed)

                if y_best < threshold:
                    logger.info(f'Solution is found, cause threshold {threshold} is touched')
                    break
                if it_num > it_limit:
                    logger.error(f'Solution is NOT found, iterations limit {it_limit} is exceed. y_best = {y_best} threshold = {threshold}')
                    break

                self.derivatives, der_vec = self.evaluate_derivatives_on_edges()
                check_for_nan(der_vec=der_vec)

                F_ = np.diag(der_vec)
                B_F_Bt = np.dot(np.dot(self.B, F_), self.Bt)
                det_B_F_Bt = np.linalg.det(B_F_Bt)
                rnd_seed = int(det_B_F_Bt) % 65536
                p_residuals = self.pt_residual_vec[:,0]
                logger.debug(f'det B = {det_B_F_Bt}')

                inv_B_F_Bt = np.linalg.inv(B_F_Bt)
                check_for_nan(inv_B_F_Bt=inv_B_F_Bt)
                check_for_nan(p_residuals=p_residuals)

                dx = -1 * np.matmul(inv_B_F_Bt, p_residuals).reshape((len(dx), 1))
                check_for_nan(dx=dx)

            self.attach_results_to_schema()
        except Exception as e:
            logger.error(e, exc_info=True)

        self.op_result = scop.OptimizeResult(success=y_best < threshold, fun=y_best, x=x_best, nfev=it_num)
        logger.info(f'Gradient descent result is {self.op_result}')


    def evaluate_derivatives_on_edges(self):
        rez = dict()
        rez_vec = np.zeros(len(self.edge_list))
        for i, (u, v) in enumerate(self.edge_list):
            p, t = self.pt_on_tree[u]
            x = self.edges_x[(u, v)]
            dx = 1e-3
            # if (u, v) in self.span_tree:
            #     p_, t_ = self.pt_on_tree[v]
            # elif (u, v) in self.chordes:
            #     p_, t_ = self.pt_on_chords_ends[(u, v)]
            # else:
            #     logger.error('Something wrong with graph, there is an edge neither in edges nor in chordes. It should not be')
            #     assert False

            edge_func = self.forward_edge_functions[(u, v)]
            p_, t_ =  edge_func(p, t, x)
            p__, t__ =  edge_func(p, t, x + dx)
            dpdx =  (p__ - p_) / dx
            rez[(u, v)] = dpdx
            rez_vec[i] = dpdx
        return rez, rez_vec

    def do_save_intermediate_results(self):
        if self.imd_rez_df is None:
            cols = list(self.pt_on_tree.keys())
            self.imd_rez_df = pd.DataFrame(columns=cols)
        row = {k:v[0] for k, v in self.pt_on_tree.items()}
        self.imd_rez_df = self.imd_rez_df.append(row, ignore_index=True)


    def perform_self_test_for_1stCL(self):
        resd_1stCL = self.evaluate_1stCL_residual()
        x_sum = sum(map(abs, self.edges_x.values()))
        if np.isnan(x_sum).any():
            logger.error(f'We cant check 1stCL cause edge flows vec contains NaN! edges_x = {self.edges_x}')

        if abs(resd_1stCL) > 1e-7 * x_sum:
            logger.error(f"Solution violates Kirchhoff''s first law, residual = {resd_1stCL}")
            assert False

    def build_circuit_matrix(self):
        # just for remember self.edge_list = self.span_tree + self.chordes
        c = len(self.chordes)
        if c==0:
            logger.debug('is finished, graph is a tree')
            return None
        m = len(self.edge_list)
        B = np.zeros((c, m))
        B[:,-c:] = np.identity(c)
        A_tree_inv = np.linalg.inv(self.A_tree)
        B_tree_transp = -1 * np.dot(A_tree_inv, self.A_chordes)
        B_tree = np.transpose(B_tree_transp)
        B[:,:m-c] = B_tree
        if np.isnan(B).any():
            logger.error(f'Something is wrong! Circuit matrix contains NaN! B = {B}')
            assert False

        return B

    def build_tree_travers(self, di_tree, root):
        di_edges = set(di_tree)
        undirected_tree = nx.Graph(di_tree)
        tree_travers = []
        for u, v in nx.algorithms.traversal.edgebfs.edge_bfs(undirected_tree, root):
            if (u, v) in di_edges:
                tree_travers += [(u, v, 1)]
            else:
                if not (v, u) in di_edges:
                    logger.error(f'Cannot find edge ({u}, {v}) in the tree')
                    assert False
                tree_travers += [(v, u, -1)]
        return tree_travers

    def split_graph(self, graph):
        G = nx.Graph(graph)

        t_ = nx.minimum_spanning_tree(G)
        te_ = set(t_.edges())

        tl, cl = [], []
        for e in self.graph.edges():
            e_ = (e[1], e[0])
            if e in te_ or e_ in te_:
                tl += [e]
            else:
                cl += [e]

        if (len(tl) != len(G.nodes)-1) or (len(tl) + len(cl) != len(G.edges)):
            logger.error(f'Cannot split graph! Tree edges = {len(tl)}, chordes = {cl}, G.nodes = {len(G.nodes)}, G.edges = {len(G.edges)}')
            assert False
        return tl, cl

    def transform_multi_di_graph_to_equal_di_graph(self, zzzz):
        MDG = nx.MultiDiGraph(zzzz, data=True)
        if type(zzzz) == nx.DiGraph:
            for u, v in zzzz.edges:
                assert zzzz[u][v]['obj'] is MDG[u][v][0]['obj']
        elif type(zzzz) == nx.MultiDiGraph:
            for u, v, k in zzzz.edges:
                assert zzzz[u][v][k]['obj'] is MDG[u][v][k]['obj']

        MUDG = nx.MultiGraph()
        MUDG.add_nodes_from(MDG)
        # obj_mdg = {id(MDG[u][v][k]['obj']) :(u, v, k) for (u, v, k) in MDG.edges}
        nodes_order = dict(zip(MDG.nodes, range(len(MDG.nodes))))
        edge_mapping = {}
        for (u, v, k) in MDG.edges:
            u_, v_ = u, v
            if nodes_order[u] > nodes_order[v]:
                u_, v_ = v, u
            k_ = MUDG.add_edge(u_, v_)
            edge_mapping[u_, v_, k_] = (u, v, k)
        assert len(MDG.edges) == len(MUDG.edges)

        rez = nx.DiGraph()
        rez.add_nodes_from(zzzz.nodes(data=True))
        for _u, _v, _k in MUDG.edges:
            u, v, k = edge_mapping[(_u, _v, _k)]
            e = MDG[u][v][k]
            if _k==0:
                # rez.add_edge(u, v, k=k, **e)
                rez.add_edge(u, v, **e)
                self.result_edges_mapping[(u, v, k)] = (u, v)
            else:
                mn = f'mock_node{len(self.mock_nodes)}'
                self.mock_nodes += [mn]
                rez.add_node(mn, obj=vrtxs.HE2_ABC_GraphVertex())
                rez.add_edge(u, mn, **e)
                rez.add_edge(mn, v, obj=HE2_MockEdge())
                self.mock_edges += [(mn, v)]
                self.result_edges_mapping[(u, v, k)] = (u, mn)
        return rez

    def add_root_to_graph(self, graph):
        p_node_found = False
        self.mock_nodes += [Root]
        G = nx.DiGraph(graph)
        G.add_node(Root, obj=None)
        for n in G.nodes:
            obj = G.nodes[n]['obj']
            if isinstance(obj, vrtxs.HE2_Boundary_Vertex) and obj.kind == 'P':
                new_obj = vrtxs.HE2_ABC_GraphVertex()
                G.nodes[n]['obj'] = new_obj
                G.add_edge(Root, n, obj=HE2_MockEdge(obj.value))
                self.initial_edges_x[(Root, n)] = self.known_Q[n]
                self.mock_edges += [(Root, n)]
                p_node_found = True
        if not p_node_found:
            logger.error('There must be a a node with constrained pressure! Solve cannot be performed')
            assert False

        logger.debug('is finished')
        return G

    def build_static_Q_vec(self, G):
        q_vec = np.zeros((len(self.node_list)-1, 1))
        for i, node in enumerate(G.nodes):
            obj = G.nodes[node]['obj']
            if isinstance(obj, vrtxs.HE2_Boundary_Vertex):
                assert obj.kind == 'Q'
                q_vec[i] = obj.value if obj.is_source else -obj.value
        logger.debug('is finished')
        return q_vec

    def build_incidence_matrices(self):
        nodelist = self.node_list
        tree_edgelist = self.span_tree
        chordes_edgelist = self.chordes

        A_full = nx.incidence_matrix(self.span_tree, nodelist=nodelist, edgelist=tree_edgelist, oriented=True)
        A_full = -1 * A_full.toarray()
        A_truncated = A_full[:-1]

        A_chordes_full = nx.incidence_matrix(self.chordes, nodelist=nodelist, edgelist=chordes_edgelist, oriented=True)
        A_chordes_full = -1 * A_chordes_full.toarray()
        A_chordes_truncated = A_chordes_full[:-1]
        logger.debug('is finished')
        return A_truncated, A_chordes_truncated

    def evalute_pressures_by_tree(self):
        pt = dict()
        pt[Root] = (0, 20)  # TODO: get initial T from some source

        for u, v, direction in self.tree_travers:
            obj = self.graph[u][v]['obj']
            if not isinstance(obj, abc.HE2_ABC_GraphEdge):
                logger.error(f'({u}, {v}) graph edge cannot evaluate its pressure drop, so we cannot evaluate pressures on the tree')
                assert False
            known, unknown = u, v
            if v in pt:
                known, unknown = v, u

            if (unknown in pt):
                logger.error(f'We dont know pressure in {unknown} node, so we cannot propagate pressure down by the tree')

            p_kn, t_kn = pt[known]
            if np.isnan(p_kn):
                logger.warning(f'P_known is NaN! Edge is ({u}, {v}), known is {known}')

            x = self.edges_x[(u, v)]
            if u == known:
                edge_func = self.forward_edge_functions[(u, v)]
            else:
                edge_func = self.backward_edge_functions[(u, v)]
            p_unk, t_unk = edge_func(p_kn, t_kn, x)

            # row = dict(known=known, unknown=unknown, x=x, p_known=p_kn, p_unknown=p_unk)
            # self.df_edge_func = self.df_edge_func.append(row, ignore_index=True)
            # if known == 'PAD_39' and unknown == 'intake_pad_59':
            #     print(x, p_kn, p_unk)

            if np.isnan(p_unk):
                logger.warning(f'edge_func returns NaN! Edge is ({u}, {v}), known is {known}')

            pt[unknown] = (p_unk, t_unk)
        return pt

    def evalute_chordes_pressure_residual(self):
        # if self.C == 0:
        #     return 0
        d = dict()
        pt_v1 = np.zeros((len(self.chordes), 2))
        pt_v2 = np.zeros((len(self.chordes), 2))
        for i, (u, v) in enumerate(self.chordes):
            x = self.edges_x[(u, v)]
            obj = self.graph[u][v]['obj']
            if not isinstance(obj, abc.HE2_ABC_GraphEdge):
                logger.error(f'{obj} on ({u}, {v}) graph edge cannot evaluate its pressure drop, so we cannot evaluate pressures on the tree')
                assert False
            p_u, t_u = self.pt_on_tree[u]
            p_v, t_v = obj.perform_calc_forward(p_u, t_u, x)

            # row = dict(known=u, unknown=v, x=x, p_known=p_u, p_unknown=p_v)
            # self.df_edge_func = self.df_edge_func.append(row, ignore_index=True)

            pt_v1[i,0] = p_v
            pt_v1[i,1] = t_v
            pt_v2[i,0] = self.pt_on_tree[v][0]
            pt_v2[i,1] = self.pt_on_tree[v][1]
            d[(u, v)] = p_v, t_v
        pt_residual_vec = pt_v1 - pt_v2
        return pt_residual_vec, d

    def attach_results_to_schema(self):
        if self.pt_on_tree is None:
            logger.warning('Cannot get results after solve')
            return

        for u, pt in self.pt_on_tree.items():
            if u in self.mock_nodes:
                continue
            obj = self.schema.nodes[u]['obj']
            obj.result = dict(P_bar=pt[0], T_C=pt[1])

            Q = 0
            if isinstance(obj, vrtxs.HE2_Boundary_Vertex) and obj.kind == 'P':
                Q = self.edges_x[(Root, u)]
            if isinstance(obj, vrtxs.HE2_Boundary_Vertex) and obj.kind == 'Q':
                Q = obj.value if obj.is_source else -obj.value
            obj.result['Q'] = Q

        if type(self.schema) == nx.DiGraph:
            for u, v in self.schema.edges:
                obj = self.schema[u][v]['obj']
                x = self.edges_x[(u, v)]
                obj.result = dict(x=x)
        elif isinstance(self.schema, nx.MultiDiGraph):
            for u, v, k in self.schema.edges:
                obj = self.schema[u][v][k]['obj']
                _u, _v = self.result_edges_mapping[(u, v, k)]
                x = self.edges_x[(_u, _v)]
                obj.result = dict(x=x)
        logger.debug('is finished')


    def evaluate_1stCL_residual(self):
        residual = 0
        G = self.graph
        nodes = set(G.nodes())
        nodes -= {Root}
        Q_net_balance = 0
        for n in list(nodes) + [Root]:
            if n != Root:
                Q = 0
                obj = G.nodes[n]['obj']
                if isinstance(obj, vrtxs.HE2_Boundary_Vertex) and obj.kind == 'P':
                    continue

                if isinstance(obj, vrtxs.HE2_Boundary_Vertex) and obj.kind == 'Q':
                    Q = obj.value if obj.is_source else -obj.value

                Q_net_balance += Q
            else:
                Q = -Q_net_balance

            X_sum = 0
            for u, v in G.in_edges(n):
                X_sum -= self.edges_x[(u, v)]
            for u, v in G.out_edges(n):
                X_sum += self.edges_x[(u, v)]
            residual += abs(Q - X_sum)

        return residual

    def evaluate_2ndCL_residual(self):
        residual = 0
        G = self.graph
        for (u, v) in G.edges():
            edge_obj = G[u][v]['obj']
            x = self.edges_x[(u, v)]
            p_u, t_u = self.pt_on_tree[u]
            p_v, t_v = self.pt_on_tree[v]
            p, t = edge_obj.perform_calc_forward(p_u, t_u, x)
            residual += abs(p - p_v)
        return residual

    def gimme_random_step(self, rand_seed):
        logger.info(f'randseed is {rand_seed}')
        np.random.seed(rand_seed)
        return np.random.uniform(0.1, 0.5)


    def evaluate_and_set_new_fluids(self):
        G = self.graph
        mr = self.fluids_move_rate
        cocktails, srcs = mixer.evalute_network_fluids_with_root(G, self.edges_x)
        src_fluids = [G.nodes[n]['obj'].fluid for n in srcs]
        for (u, v), cktl in cocktails.items():
            obj = G.nodes[u][v]['obj']
            fluidA = fl.dot_product(list(zip(cktl, src_fluids)))
            fluidB = obj.fluid
            fluidC = fl.dot_product([(1-mr, fluidA), (mr, fluidB)])
            obj.fluid = fluidC