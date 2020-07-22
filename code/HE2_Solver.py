import numpy as np
import networkx as nx
import scipy.optimize as scop

class HE2_Solver():

    def __init__(self):
        pass

    def solve(self, task):
        self.schema = task.schema
        self.fluids = task.fluids
        self.graph = task.schema.graph
        self.boundaries = task.boundaries
        self.attach_fluids_to_inlets()
        self.attach_boundaries_to_nodes()
        self.add_root_to_graph()

        self.Q = self.build_Q_vec()
        self.span_tree, self.chordes = self.split_graph()
        self.tiers = self.build_tiers(self.span_tree)
        A_tree = nx.incidence_matrix(self.span_tree)
        self.A_inv = np.linalg.inv(A_tree)
        self.A_chordes = nx.incidence_matrix(self.chordes)

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

    def split_graph(self):
        return None, None

    def add_root_to_graph(self):
        pass

    def attach_fluids_to_inlets(self):
        pass

    def attach_boundaries_to_nodes(self):
        pass

    def build_Q_vec(selfself):
        pass

    def evalute_tier_pressure_drop(self, edges, p, x, fluids):
        pass

    def restore_fluids(self, Q, x):
        pass
