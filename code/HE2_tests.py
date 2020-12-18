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
import HE2_tools as tools
import pandas as pd
import HE2_schema_maker as maker
import HE2_MixFluids as mixer
import HE2_Fit


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
        G, n_dict = tools.generate_random_net_v0(randseed=42, N=6, E=8, SNK=1, SRC=2, SEGS=2)
        solver = HE2_Solver(G)
        solver.solve()
        shifts = dict(zip(G.nodes, [(0, -0.05), (0, -0.05), (0.125, 0), (0.15, 0.025), (0, 0.05), (0, 0.05)]))
        # tools.draw_solution(G, shifts, **n_dict)
        tools.check_solution(G)

    def handle_error(self, G, rs, msg, n_dict, op_result=None):
        print('-' * 80)
        print(f'{msg} Randseed=={rs}')
        if op_result:
            print(op_result)
        # tools.draw_solution(G, None, **n_dict)

    def test_14(self):
        errs = []
        cant_solve = []
        for rs in range(50):
            G, n_dict = tools.generate_random_net_v0(randseed=rs, N=6, E=8, SNK=1, SRC=2, SEGS=2)
            solver = HE2_Solver(G)
            solver.solve()
            op_result = solver.op_result
            if op_result.fun > 1e-3:
                self.handle_error(G, rs, 'Cant solve.', n_dict, op_result)
                cant_solve += [rs]
                continue

            resd1, resd2 = tools.check_solution(G)
            if resd1 + resd2 > 1e-3:
                self.handle_error(G, rs, f'1stCL residual is {resd1: .2f}, 2ndCL residual is {resd2: .2f}', n_dict)
                errs += [rs]
                continue

        print('-' * 80)
        print('errs', errs)
        print('cant_solve', cant_solve)

    # def test_15(self):
    #     methods = ['SLSQP','BFGS','L-BFGS-B','Powell','CG','trust-constr','Nelder-Mead','TNC','COBYLA']
    #     for m1, m2 in product(methods, methods):
    #         errs = []
    #         for rs in [23, 34]:
    #             G, n_dict = tools.generate_random_net_v0(randseed=rs, N=6, E=8, SNK=1, SRC=2, SEGS=2)
    #             solver = HE2_Solver(G)
    #             solver.solve_2methods(m1, m2)
    #             op_result = solver.op_result
    #             resd1 = solver.check_solution()
    #             if resd > 1e-3:
    #                 errs += [rs]
    #         print('-' * 80)
    #         print(m1, m2)
    #         print(len(errs), errs)

    def test_16(self):
        df_pipes = pd.read_csv('..\\data\\rs_1750000976.csv')
        df_bnds  = pd.read_csv('..\\data\\boundaries.csv')
        G = maker.make_schema_from_OISPipe_dataframes(df_pipes, df_bnds)

        solver = HE2_Solver(G)
        solver.solve()
        op_result = solver.op_result
        if op_result.fun > 1e-3:
            self.handle_error(G, None, 'Cant solve.', None, op_result)

        resd1, resd2 = tools.check_solution(G)
        assert resd1 + resd2 < 1e-3, f'1stCL residual is {resd1: .2f}, 2ndCL residual is {resd2: .2f}'


    def test_17(self):
        df_pipes = pd.read_csv('..\\data\\rs_1750000976.csv')
        df_bnds  = pd.read_csv('..\\data\\boundaries.csv')
        G = maker.make_multigraph_schema_from_OISPipe_dataframes(df_pipes, df_bnds)

        solver = HE2_Solver(G)
        solver.solve()
        op_result = solver.op_result
        assert op_result.fun < 1e-3, 'Cant solve'

        resd1, resd2 = tools.check_solution(G)
        assert resd1 + resd2 < 1e-3, f'1stCL residual is {resd1: .2f}, 2ndCL residual is {resd2: .2f}'

    def test_18(self):
        for rs in range(50):
            G, n_dict = tools.generate_random_net_v1(randseed=rs, P_CNT=2)

            solver = HE2_Solver(G)
            solver.solve()
            op_result = solver.op_result
            resd1, resd2 = tools.check_solution(G)
            print(f'Randseed is {rs}, fun is {op_result.fun: .2f}, 1stCL residual is {resd1: .2f}, 2ndCL residual is {resd2: .2f}', n_dict['p_nodes'])
        # tools.draw_solution(G, None, **n_dict)
        # assert op_result.fun < 1e-3, 'Cant solve'
        #
        # resd1, resd2 = tools.check_solution(G)
        # assert resd1 + resd2 < 1e-3, f'1stCL residual is {resd1: .2f}, 2ndCL residual is {resd2: .2f}'

    def test_19(self):
        for rs in list(range(20))[5:]:
            G, n_dict = tools.generate_random_net_v1(randseed=rs, P_CNT=2, N=15, E=17)

            solver = HE2_Solver(G)
            solver.solve()
            op_result = solver.op_result

            if op_result.fun > 1e-3:
                continue
            print(rs)
            print(op_result)

            G2 = tools.build_dual_schema_from_solved(G, **n_dict)
            solver = HE2_Solver(G2)
            solver.solve()
            print(op_result)

            for n in G.nodes:
                P1 = G.nodes[n]['obj'].result['P_bar']
                P2 = G2.nodes[n]['obj'].result['P_bar']
                if abs(P1 - P2) > 1e-2:
                    print(f'{n}, {P1: .3f}, {P2: .3f}, ')

            for u, v, k in G.edges:
                x1 = G[u][v][k]['obj'].result['x']
                x2 = G[u][v][k]['obj'].result['x']
                if abs(x1 - x2) > 1e-2:
                    print(f'{u} -> {v}, {x1: .3f}, {x2: .3f}, ')
            print('-'*80)


class TestFluidMixer(unittest.TestCase):
    def setUp(self):
        pass

    def transform_solution_for_mixer(self, solution_graph):
        x_dict = dict()
        G0 = solution_graph
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
        return G, x_dict

    def test_20(self):
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

        G, x_dict = self.transform_solution_for_mixer(G)

        cocktails, srcs = mixer.evalute_network_fluids_wo_root(G, x_dict)
        for k in cocktails:
            self.assertEqual(len(cocktails[k]), 1)
            self.assertAlmostEqual(cocktails[k][0], 1)

    def test_21(self):
        inlets = dict(KNS_0=vrtxs.HE2_Source_Vertex('P', 200, 'water', T=20))
        outlets = dict(well_0=vrtxs.HE2_Boundary_Vertex('Q', 30))
        outlets.update(well_1=vrtxs.HE2_Boundary_Vertex('Q', 10))
        juncs = dict(junc_0=vrtxs.HE2_ABC_GraphVertex())

        G = nx.DiGraph()  # Di = directed
        for k, v in {**inlets, **outlets, **juncs}.items():
            G.add_node(k, obj=v)

        G.add_edge('KNS_0', 'junc_0', obj=HE2_WaterPipe([300], [10], [0.1], [1e-5]))
        G.add_edge('junc_0', 'well_0', obj=HE2_WaterPipe([200], [10], [0.1], [1e-5]))
        G.add_edge('junc_0', 'well_1', obj=HE2_WaterPipe([200], [10], [0.1], [1e-5]))

        solver = HE2_Solver(G)
        solver.solve()

        G, x_dict = self.transform_solution_for_mixer(G)

        cocktails, srcs = mixer.evalute_network_fluids_wo_root(G, x_dict)
        for k in cocktails:
            self.assertEqual(len(cocktails[k]), 1)
            self.assertAlmostEqual(cocktails[k][0], 1)

    def test_22(self):
        inlets = dict(KNS_0=vrtxs.HE2_Source_Vertex('Q', 30, 'water', T=20))
        inlets.update(KNS_1=vrtxs.HE2_Source_Vertex('Q', 10, 'water', T=20))
        outlets = dict(well_0=vrtxs.HE2_Boundary_Vertex('P', 200))
        juncs = dict(junc_0=vrtxs.HE2_ABC_GraphVertex())

        G = nx.DiGraph()  # Di = directed
        for k, v in {**inlets, **outlets, **juncs}.items():
            G.add_node(k, obj=v)

        G.add_edge('KNS_0', 'junc_0', obj=HE2_WaterPipe([300], [10], [0.1], [1e-5]))
        G.add_edge('KNS_1', 'junc_0', obj=HE2_WaterPipe([300], [10], [0.1], [1e-5]))
        G.add_edge('junc_0', 'well_0', obj=HE2_WaterPipe([200], [10], [0.1], [1e-5]))

        solver = HE2_Solver(G)
        solver.solve()

        G, x_dict = self.transform_solution_for_mixer(G)

        cocktails, srcs = mixer.evalute_network_fluids_wo_root(G, x_dict)
        self.assertEqual(len(srcs), len(inlets))
        for k in cocktails:
            self.assertEqual(len(cocktails[k]), len(srcs))

        ethalon = {('KNS_0', 'junc_0'): np.array([1., 0.])}
        ethalon.update({('KNS_1', 'junc_0'):np.array([0., 1.])})
        ethalon.update({('junc_0', 'well_0'): np.array([0.75, 0.25])})
        ethalon.update({'well_0': np.array([0.75, 0.25])})

        self.assertEqual(len(cocktails), len(ethalon))
        for k in cocktails:
            self.assertAlmostEqual(np.linalg.norm(cocktails[k] - ethalon[k]), 0)

    def test_24(self):
        for rs in range(10):
            G, n_dict = tools.generate_random_net_v0(randseed=rs)
            solver = HE2_Solver(G)
            solver.solve()
            op_result = solver.op_result
            if op_result.fun > 1e-3:
                continue
            G, x_dict = self.transform_solution_for_mixer(G)
            cocktails, srcs = mixer.evalute_network_fluids_wo_root(G, x_dict)
            rez = tools.check_fluid_mixation(G, x_dict, cocktails, srcs)
            self.assertTrue(rez)

class TestFit(unittest.TestCase):
    def setUp(self):
        pass

    def test_25(self):
        input_df = pd.read_csv('..\\data\\input_df.csv')
        fitter = HE2_Fit.HE2_PMNetwork_Model(input_df)
        rez = fitter.fit()
        print(rez)

    def test_26(self):
        # Newton-CG, dogleg, trust-ncg, trust-krylov, trust-exact не хочут, Jacobian is required
        methods = ['SLSQP', 'BFGS', 'L-BFGS-B', 'Powell', 'CG', 'trust-constr', 'Nelder-Mead', 'TNC', 'COBYLA']
        input_df = pd.read_csv('..\\data\\input_df.csv')
        best_ys = []
        for meth in methods:
            print(meth)
            fitter = HE2_Fit.HE2_PMNetwork_Model(input_df, method=meth)
            rez = fitter.fit()
            print(meth, rez)
            best_ys += [rez.fun]
        for m, y in zip(methods, best_ys):
            print(m ,y)

    def test_27(self):
        # Bounds on variables for L - BFGS - B, TNC, SLSQP, Powell, and trust-constr methods.There are two ways to specify the bounds:
        # 2. Sequence of(min, max) pairs for each element in x. None is used to specify no bound.
        methods = ['SLSQP', 'L-BFGS-B', 'Powell', 'trust-constr', 'TNC']
        input_df = pd.read_csv('..\\data\\input_df.csv')
        best_ys = []
        for meth in methods:
            print(meth)
            fitter = HE2_Fit.HE2_PMNetwork_Model(input_df, method=meth, use_bounds=True)
            rez = fitter.fit()
            print(meth, rez)
            best_ys += [rez.fun]
        for m, y in zip(methods, best_ys):
            print(m ,y)

    def test_28(self):
        input_df = pd.read_csv('..\\data\\input_df.csv')
        boundless_methods = ['COBYLA', 'L-BFGS-B', 'trust-constr', ]
        bounded_methods = ['SLSQP', 'trust-constr']
        methods = []
        methods += list(zip(bounded_methods, [True]*100))
        methods += list(zip(boundless_methods, [False]*100))
        methods = methods[:1]
        best_ys = []
        for meth, b in methods:
            print(meth)
            fitter = HE2_Fit.HE2_PMNetwork_Model(input_df, method=meth, use_bounds=b)
            fitter.max_it = 2000
            rez = fitter.fit()
            print(meth, rez)
            best_ys += [rez.fun]
            rez_df = fitter.best_rez_df
            rez_df.to_csv(f'..\\data\\rez_df{meth}{b}.csv')

        for (m, b), y in zip(methods, best_ys):
            print(m, b, y)

    def test_29(self):
        input_df = pd.read_csv('..\\data\\input_df.csv')
        fitter = HE2_Fit.HE2_PMNetwork_Model(input_df, fit_version=1)
        fitter.max_it = 2000
        rez = fitter.fit()
        print(rez)

    def test_30(self):
        input_df = pd.read_csv('..\\data\\input_df.csv')
        fitter = HE2_Fit.HE2_PMNetwork_Model(input_df, fit_version=2)
        fitter.max_it = 2000
        rez = fitter.fit()
        print(rez)

    def test_31(self):
        input_df = pd.read_csv('..\\data\\input_df.csv')
        fitter = HE2_Fit.HE2_PMNetwork_Model(input_df, fit_version=3)
        fitter.max_it = 2000
        rez = fitter.fit()
        print(rez)


    def test_32(self):
        input_df = pd.read_csv('..\\data\\input_df.csv')
        fitter = HE2_Fit.HE2_PMNetwork_Model(input_df, fit_version=4)
        fitter.max_it = 20000
        rez = fitter.fit()
        print(rez)

    def test_33(self):
        input_df = pd.read_csv('..\\data\\input_df.csv')
        fitter = HE2_Fit.HE2_PMNetwork_Model(input_df, fit_version=5)
        fitter.max_it = 20000
        rez = fitter.fit()
        print(rez)

    def test_34(self):
        input_df = pd.read_csv('..\\data\\input_df.csv')
        fitter = HE2_Fit.HE2_PMNetwork_Model(input_df, method='SLSQP', use_bounds=True)
        fitter.max_it = 5
        rez = fitter.fit()
        rez_df = fitter.best_rez_df
        rez_df.to_csv('..\\data\\rez_df.csv')
        print(rez)


if __name__ == "__main__":
    pipe_test = TestFit()
    pipe_test.test_28()

    # unittest.main()
