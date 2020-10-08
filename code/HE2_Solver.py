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
        self.graph = self.add_root_to_graph()
        self.span_tree, self.chordes = self.split_graph(self.graph)
        self.A_tree = self.buildTreeIncMatrix()
        self.A_inv = np.linalg.inv(self.A_tree)
        self.Q = self.build_Q_vec(self.graph)
        x = np.matmul(self.A_inv, self.Q)
        self.tree_x = dict(zip(self.a_tree_edgelist, x))
        self.evalute_pressures_by_tree()

        # self.tiers = self.build_tiers(self.span_tree)
        self.A_chordes = nx.incidence_matrix(self.chordes)
        return

        def target(x_chordes):
            assert(len(x_chordes) == len(self.chordes))
            Q_iter = self.Q - self.A_chordes * x_chordes
            x_tree = self.A_inv * Q_iter
            x = np.concatenate(x_tree, x_chordes)
            fluids = self.restore_fluids(Q_iter, x)
            p = np.zeros(self.graph.nodes_count)
            p[0] = self.boundaries.root_P
            for tier in self.tiers:
                p = self.evalute_tier_pressure_drop(tier, p, x, fluids)
            p_ = self.evalute_tier_pressure_drop(self.chordes, p, x, fluids)
            residual = np.sum(np.square(p_ - p))
            return residual

        x0 = np.zeros(len(self.chordes))
        op_result = scop.minimize(target, x0, method='Nelder-Mead')
        return op_result

    def build_tiers(self, tree):
        pass

    def split_graph(self, graph):
        G = nx.Graph(graph)

        t_ = nx.minimum_spanning_tree(G)
        self.tree_travers = list(nx.algorithms.traversal.edgebfs.edge_bfs(t_, Root))
        te_ = set(t_.edges())

        c_ = bool_ops.difference(G, t_)
        ce_ = set(c_.edges())

        assert len(te_ & ce_) == 0
        T = nx.DiGraph([e for e in G.edges() if e in te_ or reversed(e) in te_])
        C = nx.DiGraph([e for e in G.edges() if e in ce_ or reversed(e) in ce_])

        return T, C

    def add_root_to_graph(self):
        G = nx.DiGraph(self.schema)
        G.add_node(Root, obj=None)
        for n in G.nodes:
            obj = G.nodes[n]['obj']
            if isinstance(obj, vrtxs.HE2_Boundary_Vertex) and obj.kind == 'P':
                new_obj = vrtxs.HE2_ABC_GraphVertex()
                G.nodes[n]['obj'] = new_obj
                G.add_edge(Root, n, obj=HE2_MockEdge(obj.value))
        return G


    def build_Q_vec(self, G):
        n = len(G.nodes)
        q_vec = np.zeros(n-1)
        for i, node in enumerate(G.nodes):
            obj = G.nodes[node]['obj']
            if isinstance(obj, vrtxs.HE2_Boundary_Vertex):
                assert obj.kind == 'Q'
                q_vec[i] = obj.value if obj.is_source else -obj.value
        rez = q_vec.reshape((n-1, 1))
        return rez

    def buildTreeIncMatrix(self):
        # G = self.span_tree
        # for n in G.nodes:
        #     print(n)
        # for idx in G.edges:
        #     print(idx)
        #
        # print('-'*80)
        # G = self.graph
        # for n in G.nodes:
        #     print(n)
        # for idx in G.edges:
        #     print(idx)

        nodelist = list(self.graph.nodes)
        assert nodelist[-1] == Root
        te_ = set(self.span_tree.edges())
        edgelist =[]
        for e in self.graph.edges():
            er = tuple(reversed(e))
            if e in te_ or er in te_:
                edgelist += [e]

        A_full = -1 * nx.incidence_matrix(self.span_tree, nodelist=nodelist, edgelist=edgelist, oriented=True).toarray()
        A_truncated = A_full[:-1]
        self.a_tree_edgelist = edgelist
        return A_truncated

    def evalute_pressures_by_tree(self):
        pt = dict()
        pt[Root] = (0, 20) #TODO: get initial T from some source
        for u,v in self.tree_travers:
            obj = self.graph[u][v]['obj']
            if not isinstance(obj, abc.HE2_ABC_GraphEdge):
                assert False
            p_u, t_u = pt[u]
            x = self.tree_x[(u,v)]
            p_v, t_v = obj.perform_calc(p_u, t_u, x)
            pt[v] = (p_v, t_v)
        self.pt = pt


    def evalute_tier_pressure_drop(self, edges, p, x, fluids):
        pass

    def restore_fluids(self, Q, x):
        pass
