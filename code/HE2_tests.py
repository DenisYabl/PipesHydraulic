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
        p_fwd_0, t = pipe.perform_calc(p0_bar, t0_C, x_kgs)

        pipe.uphill_m = -10
        x_kgs = 10
        p_fwd_downhill, t = pipe.perform_calc(p0_bar, t0_C, x_kgs)

        pipe.uphill_m = 10
        x_kgs = 10
        p_fwd_uphill, t = pipe.perform_calc(p0_bar, t0_C, x_kgs)

        pipe.uphill_m = 0
        x_kgs = -10
        p_back_0, t = pipe.perform_calc(p0_bar, t0_C, x_kgs)

        pipe.uphill_m = -10
        x_kgs = -10
        p_back_downhill, t = pipe.perform_calc(p0_bar, t0_C, x_kgs)

        pipe.uphill_m = 10
        x_kgs = -10
        p_back_uphill, t = pipe.perform_calc(p0_bar, t0_C, x_kgs)

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
        p, t = pipeline.perform_calc(p0_bar, t0_C, x_kgs)
        self.assertEqual(p0_bar, p)
        self.assertEqual(t0_C, t)

    def test_8(self):
        p0_bar, t0_C, x_kgs  = 50, 20, 10
        pipe = HE2_WaterPipeSegment()
        pipe.inner_diam_m = 0.05
        pipe.roughness_m = 1e-5
        pipe.set_pipe_geometry(dx=100, dy=10)
        p1, t1 = pipe.perform_calc(p0_bar, t0_C, x_kgs)

        data = dict(dxs=[pipe.dx_m], dys = [pipe.uphill_m], diams = [pipe.inner_diam_m], rghs = [pipe.roughness_m])
        pipeline = HE2_WaterPipe(**data)
        p2, t2 = pipeline.perform_calc(p0_bar, t0_C, x_kgs)
        self.assertEqual(p1, p2)
        self.assertEqual(t1, t2)


    def test_9(self):
        p0_bar, t0_C, x_kgs  = 50, 20, 10
        pipe = HE2_WaterPipeSegment()
        pipe.inner_diam_m = 0.05
        pipe.roughness_m = 1e-5
        pipe.set_pipe_geometry(dx=100, dy=10)
        p1, t1 = pipe.perform_calc(p0_bar, t0_C, x_kgs)

        data = dict(dxs=[pipe.dx_m/2]*2, dys = [pipe.uphill_m/2]*2, diams = [pipe.inner_diam_m]*2, rghs = [pipe.roughness_m]*2)
        pipeline = HE2_WaterPipe(**data)
        p2, t2 = pipeline.perform_calc(p0_bar, t0_C, x_kgs)
        self.assertEqual(p1, p2)
        self.assertEqual(t1, t2)

class TestWaterNet(unittest.TestCase):
    def setUp(self):
        pass

    def test_1(self):
        inlets = dict(KNS_0=vrtxs.HE2_Source_Vertex('P', 200, 'water', 20))
        outlets = dict(well_0=vrtxs.HE2_Boundary_Vertex('Q', 1000))
        outlets.update(well_1=vrtxs.HE2_Boundary_Vertex('Q', 1000))
        juncs = dict(junc_0=vrtxs.HE2_ABC_GraphVertex())

        G = nx.MultiDiGraph() # Di = directed
        for k, v in {**inlets, **outlets, **juncs}.items():
            G.add_node(k, obj=v)

        # for n in G.nodes:
        #     print(G.nodes[n])

        G.add_edge('KNS_0', 'junc_0', obj=HE2_WaterPipe([100], [10], [0.1], [1e-5]))
        G.add_edge('junc_0', 'well_0', obj=HE2_WaterPipe([200], [10], [0.1], [1e-5]))
        G.add_edge('junc_0', 'well_1', obj=HE2_WaterPipe([300], [10], [0.1], [1e-5]))

        # for u,v,k in G.edges:
        #     print(G[u][v][k]['obj'])

        solver = HE2_Solver(G)
        solver.solve()
        pass


if __name__ == "__main__":
    unittest.main()
