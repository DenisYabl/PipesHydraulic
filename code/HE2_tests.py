import unittest
from HE2_Pipe import HE2_WaterPipeSegment
from HE2_Pipe import HE2_WaterPipe
import uniflocpy.uWell.uPipe as uPipe
import uniflocpy.uTools.uconst as uc
import Hydraulics.Methodics.Mukherjee_Brill as mb
import Hydraulics.Properties.Mishenko as msch
from HE2_Solver import HE2_Solver
from HE2_Fluid import HE2_DummyWater
import HE2_Vertices as vrtxs
from itertools import product
from functools import reduce
import networkx as nx
import numpy as np

class TestWaterPipe(unittest.TestCase):
    def setUp(self):
        pass

    def test_5(self):
        xs = [1, 0.01, 0.1, 10, 100]
        ds = [0.1, 0.05, 0.25, 0.5]
        rs = [1e-5, 0, 5e-6, 5e-5, 1e-4]
        conditions = product(xs, ds, rs)
        max_bias = 0
        for x, d, r in conditions:
            p0_bar = 50
            t0_C = 20
            x_kgs = x
            pipe = HE2_WaterPipeSegment()
            pipe.inner_diam_m = d
            pipe.roughness_m = r
            p_grad = pipe.calc_P_friction_gradient_Pam(p0_bar, t0_C, x_kgs)

            Q_m3day = x_kgs * 86400 / pipe.fluid.rho_wat_kgm3
            oil_params = dict(CurrentP=p0_bar, PlastT=t0_C, CurrentT=t0_C, adkuLiquidDebit=Q_m3day)
            oil_params.update(wellopVolumeWater=100, GasFactor=1e-8, PlastWaterWeight=pipe.fluid.rho_wat_kgm3/1000, wellopLiquidDensityTm=pipe.fluid.rho_wat_kgm3/1000, adkuOilDebit=0)
            oil_params.update(VolumeOilCoeff=1.2, PlastOilDynamicViscosity=1,SepOilDynamicViscosity=1, OilSaturationP=100, PlastOilWeight=0.85, SepOilWeight=0.85, GasDensity=0.8)
            tubing = dict(IntDiameter=pipe.inner_diam_m, angle=90, Roughness=pipe.roughness_m)
            mishenko = msch.Mishenko.from_oil_params(oil_params, tubing)
            misch_grad = mb.calculate(mishenko, tubing)
            p_grad = abs(p_grad)
            misch_grad = abs(misch_grad)
            bias = 1 - min(p_grad, misch_grad)/max(p_grad, misch_grad)
            if bias > max_bias:
                max_bias = bias
        self.assertLess(max_bias, 0.1)

    def test_6(self):
        p0_bar = 50
        t0_C = 20
        pipe = HE2_WaterPipeSegment()
        pipe.inner_diam_m = 0.05
        pipe.roughness_m = 1e-5
        pipe.L_m = 100

        pipe.uphill_m = 0
        x_kgs = 10
        p_fwd_0, t = pipe.calc_segment_pressure_drop(p0_bar, t0_C, x_kgs, 1)

        pipe.uphill_m = -10
        x_kgs = 10
        p_fwd_downhill, t = pipe.calc_segment_pressure_drop(p0_bar, t0_C, x_kgs, 1)

        pipe.uphill_m = 10
        x_kgs = 10
        p_fwd_uphill, t = pipe.calc_segment_pressure_drop(p0_bar, t0_C, x_kgs, 1)

        pipe.uphill_m = 0
        x_kgs = -10
        p_back_0, t = pipe.calc_segment_pressure_drop(p0_bar, t0_C, x_kgs, 1)

        pipe.uphill_m = -10
        x_kgs = -10
        p_back_downhill, t = pipe.calc_segment_pressure_drop(p0_bar, t0_C, x_kgs, 1)

        pipe.uphill_m = 10
        x_kgs = -10
        p_back_uphill, t = pipe.calc_segment_pressure_drop(p0_bar, t0_C, x_kgs, 1)

        p10m_bar = uc. Pa2bar(9.81 * 10 * 1000)
        self.assertAlmostEqual(p_fwd_downhill, p_fwd_0 + p10m_bar, 3)
        self.assertAlmostEqual(p_fwd_uphill, p_fwd_0 - p10m_bar, 3)

        self.assertAlmostEqual(p_back_downhill, p_back_0 + p10m_bar, 3)
        self.assertAlmostEqual(p_back_uphill, p_back_0 - p10m_bar, 3)

        self.assertAlmostEqual(p_back_0 + p_fwd_0, 2*p0_bar, 3)
        self.assertAlmostEqual(p_back_downhill + p_fwd_uphill, 2*p0_bar, 3)
        self.assertAlmostEqual(p_back_uphill + p_fwd_downhill, 2*p0_bar, 3)

    def test_7(self):
        data = dict(dxs=[], dys = [], diams = [], rghs = [])
        pipeline = HE2_WaterPipe(**data)
        p0_bar, t0_C, x_kgs  = 50, 20, 10
        p, t = pipeline.perform_calc_forward(p0_bar, t0_C, x_kgs)
        self.assertEqual(p0_bar, p)
        self.assertEqual(t0_C, t)

    def test_8(self):
        p0_bar, t0_C, x_kgs  = 50, 20, 10
        pipe = HE2_WaterPipeSegment()
        pipe.inner_diam_m = 0.05
        pipe.roughness_m = 1e-5
        pipe.set_pipe_geometry(dx=100, dy=10)
        p1, t1 = pipe.calc_segment_pressure_drop(p0_bar, t0_C, x_kgs, 1)

        data = dict(dxs=[pipe.dx_m], dys = [pipe.uphill_m], diams = [pipe.inner_diam_m], rghs = [pipe.roughness_m])
        pipeline = HE2_WaterPipe(**data)
        p2, t2 = pipeline.perform_calc_forward(p0_bar, t0_C, x_kgs)
        self.assertEqual(p1, p2)
        self.assertEqual(t1, t2)


    def test_9(self):
        p0_bar, t0_C, x_kgs  = 50, 20, 10
        pipe = HE2_WaterPipeSegment()
        pipe.inner_diam_m = 0.05
        pipe.roughness_m = 1e-5
        pipe.set_pipe_geometry(dx=100, dy=10)
        p1, t1 = pipe.calc_segment_pressure_drop(p0_bar, t0_C, x_kgs, 1)

        data = dict(dxs=[pipe.dx_m/2]*2, dys = [pipe.uphill_m/2]*2, diams = [pipe.inner_diam_m]*2, rghs = [pipe.roughness_m]*2)
        pipeline = HE2_WaterPipe(**data)
        p2, t2 = pipeline.perform_calc_forward(p0_bar, t0_C, x_kgs)
        self.assertEqual(p1, p2)
        self.assertEqual(t1, t2)

    def test_10(self):
        for x_kgs in [-1, 0, 1]:
            for dy in [-10, 0, 10]:
                p0_bar, t0_C  = 50, 20
                pipe = HE2_WaterPipeSegment()
                pipe.inner_diam_m = 0.05
                pipe.roughness_m = 1e-5
                pipe.set_pipe_geometry(dx=100, dy=dy)
                p1, t1 = pipe.calc_segment_pressure_drop(p0_bar, t0_C, x_kgs, 1)
                p2, t2 = pipe.calc_segment_pressure_drop(p1, t1, x_kgs, -1)
                self.assertAlmostEqual(p0_bar, p2)
                self.assertAlmostEqual(t0_C, t2)

                data = dict(dxs=[pipe.dx_m/2]*2, dys = [pipe.uphill_m/2]*2, diams = [pipe.inner_diam_m]*2, rghs = [pipe.roughness_m]*2)
                pipeline = HE2_WaterPipe(**data)
                p3, t3 = pipeline.perform_calc_forward(p0_bar, t0_C, x_kgs)
                p4, t4 = pipeline.perform_calc_backward(p3, t3, x_kgs)

                self.assertAlmostEqual(p3, p1)
                self.assertAlmostEqual(t3, t1)
                self.assertAlmostEqual(p4, p2)
                self.assertAlmostEqual(t4, t2)


    def test_11(self):
        for x_kgs in [-1, 0, 1]:
            for dy in [-10, 0, 10]:
                p0_bar, t0_C  = 50, 20
                pipe = HE2_WaterPipeSegment()
                pipe.inner_diam_m = 0.05
                pipe.roughness_m = 1e-5
                pipe.set_pipe_geometry(dx=100, dy=dy)
                p1, t1 = pipe.calc_segment_pressure_drop(p0_bar, t0_C, x_kgs, 1)
                p2, t2 = pipe.calc_segment_pressure_drop(p1, t1, x_kgs, -1)
                self.assertAlmostEqual(p0_bar, p2)
                self.assertAlmostEqual(t0_C, t2)


    def test_12(self):
        gr = 98066 * 10**-5
        fr = 6648 * 10**-5
        p0_bar, t0_C  = 50, 20
        pipe_up = HE2_WaterPipeSegment()
        pipe_up.inner_diam_m = 0.05
        pipe_up.roughness_m = 1e-5
        pipe_up.set_pipe_geometry(dx=100, dy=10)

        p, t = pipe_up.calc_segment_pressure_drop(p0_bar, t0_C, 1, 1)
        self.assertAlmostEqual(p, p0_bar - gr - fr, 4)
        p, t = pipe_up.calc_segment_pressure_drop(p0_bar, t0_C, -1, 1)
        self.assertAlmostEqual(p, p0_bar - gr + fr, 4)
        p, t = pipe_up.calc_segment_pressure_drop(p0_bar, t0_C, 1, -1)
        self.assertAlmostEqual(p, p0_bar + gr + fr, 4)
        p, t = pipe_up.calc_segment_pressure_drop(p0_bar, t0_C, -1, -1)
        self.assertAlmostEqual(p, p0_bar + gr - fr, 4)

        pipe_down = HE2_WaterPipeSegment()
        pipe_down.inner_diam_m = 0.05
        pipe_down.roughness_m = 1e-5
        pipe_down.set_pipe_geometry(dx=100, dy=-10)

        p, t = pipe_down.calc_segment_pressure_drop(p0_bar, t0_C, 1, 1)
        self.assertAlmostEqual(p, p0_bar + gr - fr, 4)
        p, t = pipe_down.calc_segment_pressure_drop(p0_bar, t0_C, -1, 1)
        self.assertAlmostEqual(p, p0_bar + gr + fr, 4)
        p, t = pipe_down.calc_segment_pressure_drop(p0_bar, t0_C, 1, -1)
        self.assertAlmostEqual(p, p0_bar - gr + fr, 4)
        p, t = pipe_down.calc_segment_pressure_drop(p0_bar, t0_C, -1, -1)
        self.assertAlmostEqual(p, p0_bar - gr - fr, 4)


class TestWaterNet(unittest.TestCase):
    def setUp(self):
        pass

    def generate_random_net_v0(self, N=15, E=20, SRC=3, SNK=3, Q=20, P=200, D=0.5, H=50, L=1000, RGH=1e-4, SEGS=10):
        '''
        :param N: total nodes count
        :param E: total edges count, cannot be less than N-1
        :param SRC: sources count
        :param SNK: sinks count
        :param Q: maximum boundary flow on source (not sink!)  node
        :param P: maximum boundary pressure on one of source/sink node
        :param D: maximum pipe segment inner diameter on graph edge
        :param H: maximum pipe segment slope
        :param L: maximum pipe segment length
        :param SEGS: maximum pipe segments count
        :return: DiGraph (not Multi)
        This method produce water pipelines network with one pressure constrained node, with constant temperature
        and with some sources/sinks. Network graph can contains cycles
        '''
        E = max(E, N-1)
        J = N-SRC-SNK
        juncs = {f'junc_{i}':vrtxs.HE2_ABC_GraphVertex() for i in range(J)}
        coin = np.random.randint(0, 2)
        pressure = np.random.randint(-P, P)
        if coin == 1:
            p_nodes = {'p_node_0':vrtxs.HE2_Source_Vertex('P', pressure, 'water', 20)}
        else:
            p_nodes = {'p_node_0':vrtxs.HE2_Boundary_Vertex('P', pressure)}

        src_q = np.random.uniform(0, Q, SRC)
        total_q = sum(src_q)
        sink_q = np.random.uniform(0, 1, SNK)
        sink_q = total_q * sink_q / sum(sink_q)
        self.assertAlmostEqual(total_q, sum(sink_q))

        SRC = SRC - coin
        SNK = SNK - (1-coin)
        sources = {f'src_{i}':vrtxs.HE2_Source_Vertex('Q', src_q[i], 'water', 20) for i in range(SRC)}
        sinks = {f'sink_{i}':vrtxs.HE2_Boundary_Vertex('Q', -sink_q[i]) for i in range(SNK)}

        nodes = {**p_nodes, **sources, **sinks, **juncs}
        mapping = dict(zip(range(len(nodes)),nodes.keys()))
        while True:
            RG = nx.generators.random_graphs.gnm_random_graph(N, E, directed=True)
            UDG = nx.Graph(RG)
            if nx.algorithms.components.number_connected_components(UDG) == 1:
                break
        G = nx.relabel.relabel_nodes(RG, mapping)
        nx.set_node_attributes(G, name='obj', values=nodes)

        pipes = dict()
        for u, v in G.edges():
            segs = np.random.randint(SEGS) + 1
            Ls = np.random.uniform(1e-5, L, segs)
            Hs = np.random.uniform(-H, H, segs)
            Ds = np.random.uniform(1e-5, D, segs)
            Rs = np.random.uniform(0, RGH, segs)
            pipe = HE2_WaterPipe(Ls, Hs, Ds, Rs)
            pipes[(u,v)] = pipe
        nx.set_edge_attributes(G, name='obj', values=pipes)

        return G

    def test_10(self):
        inlets = dict(KNS_0=vrtxs.HE2_Source_Vertex('P', 200, 'water', 20))
        outlets = dict(well_0=vrtxs.HE2_Boundary_Vertex('Q', 10))
        outlets.update(well_1=vrtxs.HE2_Boundary_Vertex('Q', 10))
        juncs = dict(junc_0=vrtxs.HE2_ABC_GraphVertex())

        G = nx.DiGraph() # Di = directed
        for k, v in {**inlets, **outlets, **juncs}.items():
            G.add_node(k, obj=v)

        G.add_edge('KNS_0', 'junc_0', obj=HE2_WaterPipe([300], [10], [0.1], [1e-5]))
        G.add_edge('junc_0', 'well_0', obj=HE2_WaterPipe([200], [10], [0.1], [1e-5]))
        G.add_edge('junc_0', 'well_1', obj=HE2_WaterPipe([200], [10], [0.1], [1e-5]))

        solver = HE2_Solver(G)
        solver.solve()

        p0_bar, t0_C = 200, 20
        pipe = HE2_WaterPipeSegment()
        pipe.inner_diam_m = 0.1
        pipe.roughness_m = 1e-5
        pipe.L_m = (300*300 + 10*10)**0.5
        pipe.uphill_m = 10
        x_kgs = 20
        p_junc, t_junc = pipe.calc_segment_pressure_drop(p0_bar, t0_C, x_kgs, 1)
        pipe = HE2_WaterPipeSegment()
        pipe.inner_diam_m = 0.1
        pipe.roughness_m = 1e-5
        pipe.L_m = (200*200 + 10*10)**0.5
        pipe.uphill_m = 10
        x_kgs = 10
        p_well, t_well = pipe.calc_segment_pressure_drop(p_junc, t_junc, x_kgs, 1)

        self.assertAlmostEqual(G.nodes['KNS_0']['obj'].result['P_bar'], p0_bar)
        self.assertAlmostEqual(G.nodes['KNS_0']['obj'].result['T_C'], t0_C)

        self.assertAlmostEqual(G.nodes['junc_0']['obj'].result['P_bar'], p_junc)
        self.assertAlmostEqual(G.nodes['junc_0']['obj'].result['T_C'], t_junc)

        self.assertAlmostEqual(G.nodes['well_0']['obj'].result['P_bar'], p_well)
        self.assertAlmostEqual(G.nodes['well_0']['obj'].result['T_C'], t_well)

        self.assertAlmostEqual(G.nodes['well_1']['obj'].result['P_bar'], p_well)
        self.assertAlmostEqual(G.nodes['well_1']['obj'].result['T_C'], t_well)

        # for n in G.nodes:
        #     print(n, G.nodes[n]['obj'].result)
        # for u,v in G.edges:
        #     print(f'{u}->{v}', G[u][v]['obj'].result)


    def test_11(self):
        inlets = dict(KNS_0=vrtxs.HE2_Source_Vertex('P', 200, 'water', 20))
        outlets = dict(well_0=vrtxs.HE2_Boundary_Vertex('Q', 10))
        outlets.update(well_1=vrtxs.HE2_Boundary_Vertex('Q', 10))
        juncs = dict(junc_0=vrtxs.HE2_ABC_GraphVertex())

        G = nx.DiGraph() # Di = directed
        for k, v in {**inlets, **outlets, **juncs}.items():
            G.add_node(k, obj=v)

        G.add_edge('junc_0', 'KNS_0', obj=HE2_WaterPipe([300], [-10], [0.1], [1e-5]))
        G.add_edge('well_0', 'junc_0', obj=HE2_WaterPipe([200], [-10], [0.1], [1e-5]))
        G.add_edge('well_1', 'junc_0', obj=HE2_WaterPipe([200], [-10], [0.1], [1e-5]))

        solver = HE2_Solver(G)
        solver.solve()

        inlets = dict(KNS_0=vrtxs.HE2_Source_Vertex('P', 200, 'water', 20))
        outlets = dict(well_0=vrtxs.HE2_Boundary_Vertex('Q', 10))
        outlets.update(well_1=vrtxs.HE2_Boundary_Vertex('Q', 10))
        juncs = dict(junc_0=vrtxs.HE2_ABC_GraphVertex())

        G2 = nx.DiGraph() # Di = directed
        for k, v in {**inlets, **outlets, **juncs}.items():
            G2.add_node(k, obj=v)

        G2.add_edge('KNS_0', 'junc_0', obj=HE2_WaterPipe([300], [10], [0.1], [1e-5]))
        G2.add_edge('junc_0', 'well_0', obj=HE2_WaterPipe([200], [10], [0.1], [1e-5]))
        G2.add_edge('junc_0', 'well_1', obj=HE2_WaterPipe([200], [10], [0.1], [1e-5]))

        solver = HE2_Solver(G2)
        solver.solve()

        for n in G.nodes:
            res1 = G.nodes[n]['obj'].result
            res2 = G2.nodes[n]['obj'].result
            self.assertAlmostEqual(res1['P_bar'], res2['P_bar'])
            self.assertAlmostEqual(res1['T_C'], res2['T_C'])

        for u,v in G.edges:
            res1 = G[u][v]['obj'].result
            res2 = G2[v][u]['obj'].result
            self.assertAlmostEqual(res1['x'], -res2['x'])

    def test_12(self):
        inlets = dict(KNS_0=vrtxs.HE2_Source_Vertex('P', 200, 'water', 20))
        outlets = dict(well_0=vrtxs.HE2_Boundary_Vertex('Q', 10))
        outlets.update(well_1=vrtxs.HE2_Boundary_Vertex('Q', 10))
        juncs = dict(junc_0=vrtxs.HE2_ABC_GraphVertex())
        juncs.update(junc_1=vrtxs.HE2_ABC_GraphVertex())

        G = nx.DiGraph() # Di = directed
        for k, v in {**inlets, **outlets, **juncs}.items():
            G.add_node(k, obj=v)

        G.add_edge('junc_0', 'KNS_0', obj=HE2_WaterPipe([300], [-10], [0.1], [1e-5]))
        G.add_edge('junc_1', 'KNS_0', obj=HE2_WaterPipe([300], [-10], [0.1], [1e-5]))
        G.add_edge('junc_0', 'junc_1', obj=HE2_WaterPipe([300], [0], [0.1], [1e-5]))
        G.add_edge('well_0', 'junc_0', obj=HE2_WaterPipe([200], [-10], [0.1], [1e-5]))
        G.add_edge('well_1', 'junc_1', obj=HE2_WaterPipe([200], [-10], [0.1], [1e-5]))

        solver = HE2_Solver(G)
        solver.solve()

    def test_13(self):
        G = self.generate_random_net_v0()
        solver = HE2_Solver(G)
        solver.solve()




if __name__ == "__main__":
    # pipe_test = TestWaterNet()
    # pipe_test.test_12()

    unittest.main()
