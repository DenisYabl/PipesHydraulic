import numpy as np
import networkx as nx
import scipy.optimize as scop
from HE2_SpecialEdges import HE2_MockEdge
import HE2_Vertices as vrtxs
import networkx.algorithms.operators.binary as bool_ops
import scipy.sparse
import HE2_ABC as abc


Root = 'Root'

class HE2_Solver():
    def __init__(self, schema):
        #TODO have to implement take MultiDiGraph and convert it to equal DiGraph with some added mock edges
        self.schema = schema
        self.graph = None

    def solve(self):
        op_result = None
        self.graph = self.add_root_to_graph()
        self.N = len(self.graph)-1
        self.span_tree, self.chordes = self.split_graph(self.graph)
        self.C = len(self.chordes.edges())
        # span_tree and chordes are DiGraphs, which edges match with self.graph edges
        self.tree_travers = self.build_tree_travers(self.span_tree, Root)
        self.A_tree, self.A_chordes = self.buildIncMatrices()
        self.A_inv = np.linalg.inv(self.A_tree)
        self.Q_static = self.build_static_Q_vec(self.graph)

        def target(x_chordes):
            Q = self.Q_static
            if self.C:
                Q_dynamic = np.matmul(self.A_chordes, x_chordes).flatten()
                Q = Q - Q_dynamic
            x = np.matmul(self.A_inv, Q).flatten()
            self.tree_x = dict(zip(self.a_tree_edgelist, x))
            self.pt_on_tree =  self.evalute_pressures_by_tree()
            pt_residual_vec = self.evalute_chordes_pressure_residual(self.chordes, x_chordes, self.pt_on_tree)
            rez = np.linalg.norm(pt_residual_vec)
            return rez

        if self.C:
            x0 = np.zeros(self.C)
            # Newton-CG, dogleg, trust-ncg, trust-krylov, trust-exact не хочут, Jacobian is required
            # 'Nelder-Mead', 'CG', L-BFGS-B, TNC, COBYLA, SLSQP - плохо сходятся, fun > 2
            # 'BFGS' при низком tol плохо сходится, при тол < 1e-5 - норм, но долго, nfev 472
            # Powell сходится, 350 nfev при tol 1e-6, 550 при tol 1e-3
            # 'trust-constr' сходится, default tol, 540 nfev

            op_result = scop.minimize(target, x0, method='Powell')
            print(op_result)
            target(op_result.x)
            self.chord_x = dict(zip(self.chordes.edges(), op_result.x))
            # TODO Вот здесь надо забирать давления по ключу op_result.x из промежуточных результатов, когду они будут сохраняться
            # А пока может быть так что давления от одной итерации, а потоки от другой
        else:
            target(None)

        self.attach_results_to_schema()
        return op_result

    def build_tiers(self, tree):
        pass

    def build_tree_travers(self, di_tree, root):
        di_edges = set(di_tree.edges())
        undirected_tree = nx.Graph(di_tree)
        tree_travers = []
        for u, v in nx.algorithms.traversal.edgebfs.edge_bfs(undirected_tree, root):
            if (u, v) in di_edges:
                tree_travers += [(u, v, 1)]
            else:
                assert (v, u) in di_edges
                tree_travers += [(v, u, -1)]
        return tree_travers


    def split_graph(self, graph):
        G = nx.Graph(graph)

        t_ = nx.minimum_spanning_tree(G)
        te_ = set(t_.edges())

        tl, cl = [], []
        for u, v in self.graph.edges():
            if (u, v) in te_ or (v, u) in te_:
                tl += [(u, v)]
            else:
                cl += [(u, v)]

        T = nx.DiGraph(tl)
        C = nx.DiGraph(cl)
        assert len(T.edges()) + len(C.edges()) == len(self.graph.edges())

        return T, C

    def add_root_to_graph(self):
        self.mock_nodes = [Root]
        self.mock_edges = []
        G = nx.DiGraph(self.schema)
        G.add_node(Root, obj=None)
        for n in G.nodes:
            obj = G.nodes[n]['obj']
            if isinstance(obj, vrtxs.HE2_Boundary_Vertex) and obj.kind == 'P':
                new_obj = vrtxs.HE2_ABC_GraphVertex()
                G.nodes[n]['obj'] = new_obj
                G.add_edge(Root, n, obj=HE2_MockEdge(obj.value))
                self.mock_edges += [(Root, n)]
        return G


    def build_static_Q_vec(self, G):
        n = len(G.nodes)
        q_vec = np.zeros(n-1)
        for i, node in enumerate(G.nodes):
            obj = G.nodes[node]['obj']
            if isinstance(obj, vrtxs.HE2_Boundary_Vertex):
                assert obj.kind == 'Q'
                q_vec[i] = obj.value if obj.is_source else -obj.value
        rez = q_vec.reshape((n-1, 1))
        return rez

    def buildIncMatrices(self):
        nodelist = list(self.graph.nodes)
        assert nodelist[-1] == Root
        te_ = set(self.span_tree.edges())
        tree_edgelist, chordes_edgelist = [], []
        for e in self.graph.edges():
            er = tuple(reversed(e))
            if e in te_ or er in te_:
                tree_edgelist += [e]
            else:
                chordes_edgelist += [e]

        A_full = -1 * nx.incidence_matrix(self.span_tree, nodelist=nodelist, edgelist=tree_edgelist, oriented=True).toarray()
        A_truncated = A_full[:-1]
        self.a_tree_edgelist = tree_edgelist

        A_chordes_full = -1 * nx.incidence_matrix(self.chordes, nodelist=nodelist, edgelist=chordes_edgelist, oriented=True).toarray()
        A_chordes_truncated = A_chordes_full[:-1]
        return A_truncated, A_chordes_truncated

    def evalute_pressures_by_tree(self):
        pt = dict()
        pt[Root] = (0, 20) #TODO: get initial T from some source

        for u, v, direction in self.tree_travers:
            obj = self.graph[u][v]['obj']
            if not isinstance(obj, abc.HE2_ABC_GraphEdge):
                assert False
            known, unknown = u, v
            if v in pt:
                known, unknown = v, u

            assert not (unknown in pt)
            p_kn, t_kn = pt[known]
            x = self.tree_x[(u,v)]
            if u == known:
                p_unk, t_unk = obj.perform_calc_forward(p_kn, t_kn, x)
            else:
                p_unk, t_unk = obj.perform_calc_backward(p_kn, t_kn, x)
            pt[unknown] = (p_unk, t_unk)
        return pt

    def evalute_chordes_pressure_residual(self, chordes, chordes_x, tree_pressure):
        if (chordes is None) or (len(chordes)==0):
            return 0
        pt_u, pt_v = [], []
        for x, (u, v) in zip(chordes_x, chordes.edges()):
            obj = self.graph[u][v]['obj']
            p_u, t_u = tree_pressure[u]
            if not isinstance(obj, abc.HE2_ABC_GraphEdge):
                assert False
            p_v, t_v = obj.perform_calc_forward(p_u, t_u, x)
            pt_v += [(p_v, t_v)]
            pt_u += [(p_u, t_u)]
        pt_v_vec = np.array(pt_v)
        pt_u_vec = np.array(pt_u)
        pt_residual_vec = pt_v_vec - pt_u_vec
        return pt_residual_vec

    def attach_results_to_schema(self):
        for u, pt in self.pt_on_tree.items():
            if u in self.mock_nodes:
                continue
            obj = self.schema.nodes[u]['obj']
            obj.result = dict(P_bar=pt[0], T_C=pt[1])
        for u,v in self.schema.edges:
            obj = self.schema[u][v]['obj']
            x = None
            if (u, v) in self.tree_x:
                x = self.tree_x[(u, v)]
            elif (u, v) in self.chord_x:
                x = self.chord_x[(u, v)]

            obj.result = dict(x=x)

    def evalute_tier_pressure_drop(self, edges, p, x, fluids):
        pass

    def check_solution(self):
        G = self.schema
        for (u, v) in G.edges():
            edge_obj = G[u][v]['obj']
            u_obj = G.nodes[u]['obj']
            v_obj = G.nodes[v]['obj']
            x = edge_obj.result['x']
            p_u = u_obj.result['P_bar']
            t_u = u_obj.result['T_C']
            p_v = v_obj.result['P_bar']
            t_v = v_obj.result['T_C']
            print(u, v, f'{x:.2f}', f'{p_u:.2f}', f'{p_v:.2f}', edge_obj)
            p, t = edge_obj.perform_calc_forward(p_u, t_u, x)
            np.testing.assert_almost_equal(p, p_v)
            np.testing.assert_almost_equal(t, t_v)
