import unittest
from GraphEdges.HE2_Pipe import HE2_WaterPipe
import uniflocpy.uTools.uconst as uc
import Hydraulics.Methodics.Mukherjee_Brill as mb
import Hydraulics.Properties.Mishenko as msch
from Solver.HE2_Solver import HE2_Solver
from GraphNodes import HE2_Vertices as vrtxs
from itertools import product
import networkx as nx
import numpy as np
import pandas as pd
from Fluids import HE2_MixFluids as mixer
from Fluids.HE2_Fluid import HE2_BlackOil, HE2_DummyWater
from Solver import HE2_Fit
from Tools import HE2_Visualize as vis, HE2_tools as tools
from Tools.cachespline import create_lazy_spline_cache_f_wrapper
from Tools.HE2_tools import draw_solution, print_solution
import random


class TestWaterNet(unittest.TestCase):
    def setUp(self):
        pass

    def generate_water_graph(self):
        fluid = HE2_DummyWater()
        roughness = 1e-5
        real_diam_coefficient = 0.85

        G = nx.DiGraph()  # Di = directed
        bound_kind, bound_value = 'P', 4.8
        outlets = dict(DNS_2=vrtxs.HE2_Boundary_Vertex(bound_kind, float(bound_value)))
        inlets = dict()
        juncs_to_add = ['intake_pad_5', 'intake_pad_9_36', 'intake_pad_35', 'UDR_1', 'intake_pad_33', 'UDR_2',
                        'intake_pad_59', 'intake_pad_59_46']
        juncs_to_add += ['intake_pad_47', 'intake_pad_53', 'intake_pad_41', 'ZKL_98', 'intake_pad_49', 'intake_node_33',
                         'intake_pad_46', 'node_15']
        juncs_to_add += ['intake_node_49', 'intake_pad_57', 'intake_pad_49', 'intake_pad_134']
        juncs = {k: vrtxs.HE2_ABC_GraphVertex() for k in set(juncs_to_add)}
        for k, v in {**inlets, **outlets, **juncs}.items():
            G.add_node(k, obj=v)

        PAD_33 = 'P=3.9'
        PAD_34 = 'Q=12.7'
        PAD_5 = 'Q=2.6'
        PAD_39 = 'Q=6.5'
        PAD_49 = 'Q=16.9'
        PAD_57 = 'Q=12.5'
        pad_dict = dict(PAD_5=PAD_5, PAD_33=PAD_33, PAD_34=PAD_34, PAD_39=PAD_39, PAD_49=PAD_49, PAD_57=PAD_57)

        for pad_name, describing in pad_dict.items():
            bound_kind, bound_value = tuple(describing.split('='))
            obj = vrtxs.HE2_Source_Vertex(bound_kind, float(bound_value), fluid, 20)
            G.add_node(pad_name, obj=obj)
            inlets[pad_name] = obj

        # Трубопроводная система
        G.add_edge('PAD_5', 'intake_pad_5',
                   obj=HE2_WaterPipe([10], [-0.2], [0.143 * real_diam_coefficient], [roughness]))
        G.add_edge('intake_pad_5', 'intake_pad_9_36',
                   obj=HE2_WaterPipe([785], [-0.6], [0.199 * real_diam_coefficient], [roughness]))
        G.add_edge('intake_pad_9_36', 'intake_pad_35',
                   obj=HE2_WaterPipe([2326], [4.6], [0.199 * real_diam_coefficient], [roughness]))
        G.add_edge('intake_pad_35', 'UDR_1',
                   obj=HE2_WaterPipe([2648], [13.6], [0.199 * real_diam_coefficient], [roughness]))
        G.add_edge('UDR_1', 'DNS_2', obj=HE2_WaterPipe([200], [0.2], [0.305 * real_diam_coefficient], [roughness]))
        G.add_edge('PAD_33', 'intake_pad_33',
                   obj=HE2_WaterPipe([440], [-1.9], [0.139 * real_diam_coefficient], [roughness]))
        G.add_edge('intake_pad_33', 'UDR_2',
                   obj=HE2_WaterPipe([4394], [-4.6], [0.199 * real_diam_coefficient], [roughness]))
        G.add_edge('UDR_2', 'UDR_1', obj=HE2_WaterPipe([1087], [3.4], [0.253 * real_diam_coefficient], [roughness]))
        G.add_edge('PAD_34', 'intake_pad_134',
                   obj=HE2_WaterPipe([818], [0.6], [0.143 * real_diam_coefficient], [roughness]))
        G.add_edge('intake_pad_134', 'UDR_2',
                   obj=HE2_WaterPipe([3344], [10.2], [0.203 * real_diam_coefficient], [roughness]))
        G.add_edge('PAD_39', 'intake_pad_59',
                   obj=HE2_WaterPipe([7568], [1.6], [0.139 * real_diam_coefficient], [roughness]))
        G.add_edge('intake_pad_59', 'intake_pad_59_46',
                   obj=HE2_WaterPipe([4601], [4.4], [0.143 * real_diam_coefficient], [roughness]))
        G.add_edge('intake_pad_59_46', 'intake_pad_47',
                   obj=HE2_WaterPipe([2625], [18.4], [0.139 * real_diam_coefficient], [roughness]))
        G.add_edge('intake_pad_47', 'intake_pad_53',
                   obj=HE2_WaterPipe([2250], [-1.3], [0.203 * real_diam_coefficient], [roughness]))
        G.add_edge('intake_pad_53', 'intake_pad_41',
                   obj=HE2_WaterPipe([1391], [9], [0.203 * real_diam_coefficient], [roughness]))
        G.add_edge('intake_pad_41', 'ZKL_98',
                   obj=HE2_WaterPipe([9290], [-13.5], [0.257 * real_diam_coefficient], [roughness]))
        G.add_edge('ZKL_98', 'intake_pad_49',
                   obj=HE2_WaterPipe([600], [3], [0.257 * real_diam_coefficient], [roughness]))
        G.add_edge('intake_pad_49', 'intake_node_33',
                   obj=HE2_WaterPipe([4001], [0], [0.410 * real_diam_coefficient], [roughness]))
        G.add_edge('intake_node_33', 'intake_pad_46',
                   obj=HE2_WaterPipe([2920], [2.1], [0.410 * real_diam_coefficient], [roughness]))
        G.add_edge('intake_pad_46', 'node_15',
                   obj=HE2_WaterPipe([1181], [-5.9], [0.410 * real_diam_coefficient], [roughness]))
        G.add_edge('node_15', 'UDR_2', obj=HE2_WaterPipe([92], [0.2], [0.199 * real_diam_coefficient], [roughness]))
        G.add_edge('PAD_49', 'intake_node_49',
                   obj=HE2_WaterPipe([492], [-0.8], [0.139 * real_diam_coefficient], [roughness]))
        G.add_edge('intake_node_49', 'intake_pad_57',
                   obj=HE2_WaterPipe([3737], [-5.7], [0.143 * real_diam_coefficient], [roughness]))
        G.add_edge('intake_pad_57', 'intake_pad_49',
                   obj=HE2_WaterPipe([3852], [-4.3], [0.309 * real_diam_coefficient], [roughness]))
        G.add_edge('PAD_57', 'intake_pad_57',
                   obj=HE2_WaterPipe([370], [-1], [0.203 * real_diam_coefficient], [roughness]))
        return G

    def test_1(self):
        G = self.generate_water_graph()
        solver = HE2_Solver(G)
        solver.solve(threshold=0.05, it_limit=20)
        print_solution(G)

        G = random_reverse_edge_in_graph(G)
        solver = HE2_Solver(G)
        solver.solve(threshold=0.05, it_limit=20)
        print_solution(G)

    def test_2(self):
        G = self.generate_water_graph()
        G = add_random_edge(G)
        G = add_random_edge(G)

        solver = HE2_Solver(G)
        solver.solve(threshold=0.05, it_limit=20)
        print_solution(G)

        G = random_reverse_edge_in_graph(G)
        G = random_reverse_edge_in_graph(G)
        solver = HE2_Solver(G)
        solver.solve(threshold=0.05, it_limit=20)
        print_solution(G)

    def test_3(self):
        f = open('randoms.txt')

        diffs = []
        not_solved = []
        for i, sss in enumerate(f):
            print(f'-----------------------------------------{i}----------------------------------')
            params = sss[:-1].split(';')
            G = self.generate_water_graph()

            u, v = params[0:2]
            dx, dy, diam, rgh = tuple(map(float, params[2:6]))
            G.add_edge(u, v, obj=HE2_WaterPipe([dx],[dy],[diam],[rgh]))

            u, v = params[6:8]
            dx, dy, diam, rgh = tuple(map(float, params[8:12]))
            G.add_edge(u, v, obj=HE2_WaterPipe([dx],[dy],[diam],[rgh]))

            u1, v1, u2, v2 = params[12:16]

            solver = HE2_Solver(G)
            solver.solve(threshold=0.25, it_limit=100)
            # print_solution(G)
            x1 = solver.op_result.x
            ch1 = list(solver.chordes)
            if not solver.op_result.success:
                not_solved += [i]

            G = reverse_edge_in_graph(G, u1, v1)
            G = reverse_edge_in_graph(G, u2, v2)

            solver = HE2_Solver(G)
            solver.solve(threshold=0.25, it_limit=100)
            # print_solution(G)

            l1 = sorted(list(abs(x1.flatten())))
            l2 = []
            for u, v in ch1:
                if (u, v) in solver.edges_x:
                    l2 += [abs(solver.edges_x[(u, v)])]
                elif (v, u) in solver.edges_x:
                    l2 += [abs(solver.edges_x[(v, u)])]
                else:
                    assert False
            l2 = sorted(l2)

            diff = np.linalg.norm(np.array(l1) - np.array(l2))
            if diff > 5e-3:
                diffs += [i]

            print(diffs)
            print(not_solved)


def random_reverse_edge_in_graph(G):
    edgelist = list(G.edges)
    e = random.choice(edgelist)
    u, v = e
    G = reverse_edge_in_graph(G, u, v)
    return G

def reverse_edge_in_graph(G: nx.DiGraph, u, v):
    obj = G[u][v]['obj']
    if isinstance(obj, HE2_WaterPipe):
        dxs, dys, diams, rghs = [], [], [], []
        for seg in obj.segments:
            dxs += [seg.dx_m]
            dys += [-1 * seg.uphill_m]
            diams += [seg.inner_diam_m]
            rghs += [seg.roughness_m]
        new_obj = HE2_WaterPipe(dxs[::-1], dys[::-1], diams[::-1], rghs[::-1])
        G.remove_edge(u, v)
        G.add_edge(v, u, obj=new_obj)
    return G

def add_random_edge(G: nx.DiGraph):
    nodes = list(G.nodes)
    u = random.choice(nodes)
    nodes.remove(u)
    v = random.choice(list(G.nodes))
    dx = random.uniform(10, 1000)
    dy = random.uniform(-50, 50)
    diam = random.uniform(0.05, 0.5)
    rgh = 1e-5
    obj = HE2_WaterPipe([dx], [dy], [diam], [rgh])
    G.add_edge(u, v, obj=obj)
    return G

def generate_random_edges_and_swaps(G):
    edgelist = list(G.edges)
    nodelist = list(G.nodes)
    f = open('randoms.txt', 'w')
    for i in range(1000):
        u1 = random.choice(nodelist)
        v1 = random.choice(list(G.nodes))
        if u1==v1 or (u1, v1) in edgelist or (v1, u1) in edgelist:
            continue

        u2 = random.choice(nodelist)
        v2 = random.choice(list(G.nodes))
        if u2==v2 or (u2, v2) in edgelist or (v2, u2) in edgelist:
            continue

        if (u1, v1) == (u2, v2) or (v1, u1) == (u2, v2):
            continue

        u3, v3 = random.choice(edgelist)
        u4, v4 = random.choice(edgelist)

        if (u3, v3) == (u4, v4) or (v3, u3) == (u4, v4):
            continue

        dx = random.uniform(10, 1000)
        dy = random.uniform(-50, 50)
        diam = random.uniform(0.05, 0.5)
        rgh = 1e-5
        tpl1 = (u1, v1, dx, dy, diam, rgh)

        dx = random.uniform(10, 1000)
        dy = random.uniform(-50, 50)
        diam = random.uniform(0.05, 0.5)
        rgh = 1e-5
        tpl2 = (u2, v2, dx, dy, diam, rgh)

        tpl3 = (u3, v3, u4, v4)
        lst = list(tpl1) + list(tpl2) + list(tpl3)
        print(';'.join(map(str, lst)), file=f)

if __name__ == "__main__":

    test = TestWaterNet()
    # G = test.generate_water_graph()
    # generate_random_edges_and_swaps(G)
    test.test_3()