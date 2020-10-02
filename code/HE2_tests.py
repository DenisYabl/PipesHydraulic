import unittest
from HE2_Pipe import HE2_WaterPipeSegment
import uniflocpy.uWell.uPipe as uPipe
import uniflocpy.uTools.uconst as uc
import Hydraulics.Methodics.Mukherjee_Brill as mb
import Hydraulics.Properties.Mishenko as msch
from itertools import product

class TestWaterPipeSegment(unittest.TestCase):
    def setUp(self):
        pass

    # def test_1(self):
    #     pipe = HE2_WaterPipeSegment()
    #     pipe.inner_diam_m = 0.1
    #     pipe.roughness_m = 1e-5
    #     pipe.L_m = 100
    #     pipe.slope_m = 0
    #
    #     p, t = pipe.perform_calc(50, 20, 100)
    #     self.assertAlmostEqual(p, 49.99978631383106)
    #     self.assertAlmostEqual(t, 20)
    #
    # def test_2(self):
    #     pipe = HE2_WaterPipeSegment()
    #     pipe.inner_diam_m = 0.1
    #     pipe.roughness_m = 1e-5
    #     pipe.L_m = 100
    #     pipe.slope_m = 0
    #
    #     p, t = pipe.perform_calc(50, 20, -1)
    #     self.assertAlmostEqual(p, 50.0002, 4)
    #     self.assertAlmostEqual(t, 20)
    #
    # def test_3(self):
    #     pipe = HE2_WaterPipeSegment()
    #     pipe.inner_diam_m = 0.1
    #     pipe.roughness_m = 1e-5
    #     pipe.L_m = 100
    #     pipe.slope_m = 10
    #
    #     p, t = pipe.perform_calc(50, 20, 0)
    #     self.assertAlmostEqual(p, 50.98, 2)
    #     self.assertAlmostEqual(t, 20)
    #
    # def test_4(self):
    #     pipe = HE2_WaterPipeSegment()
    #     pipe.inner_diam_m = 0.1
    #     pipe.roughness_m = 1e-5
    #     pipe.L_m = 100
    #     pipe.slope_m = 10
    #
    #     p, t = pipe.perform_calc(-50, 20, 0)
    #     self.assertAlmostEqual(p, -49.02, 2)
    #     self.assertAlmostEqual(t, 20)
    #
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
            # pipe.L_m = 100
            # pipe.slope_m = 0
            p_grad = pipe.calc_P_friction_gradient_Pam(p0_bar, t0_C, x_kgs)

            Q_m3day = x_kgs * 86400 / pipe.fluid.rho_wat_kgm3
            oil_params = dict(CurrentP=p0_bar, PlastT=t0_C, CurrentT=t0_C, adkuLiquidDebit=Q_m3day)
            oil_params.update(wellopVolumeWater=100, GasFactor=1e-8, PlastWaterWeight=pipe.fluid.rho_wat_kgm3/1000, wellopLiquidDensityTm=pipe.fluid.rho_wat_kgm3/1000, adkuOilDebit=0)
            oil_params.update(VolumeOilCoeff=1.2, PlastOilDynamicViscosity=1,SepOilDynamicViscosity=1, OilSaturationP=100, PlastOilWeight=0.85, SepOilWeight=0.85, GasDensity=0.8)
            tubing = dict(IntDiameter=pipe.inner_diam_m, angle=90, Roughness=pipe.roughness_m)
            mishenko = msch.Mishenko.from_oil_params(oil_params, tubing)
            misch_grad = mb.calculate(mishenko, tubing)
            # self.assertAlmostEqual(abs(misch_grad), p_grad, 1)
            # print(x, d, r, p_grad, misch_grad)
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

        pipe.downhill_m = 0
        x_kgs = 10
        p_fwd_0, t = pipe.perform_calc(p0_bar, t0_C, x_kgs)

        pipe.downhill_m = 10
        x_kgs = 10
        p_fwd_downhill, t = pipe.perform_calc(p0_bar, t0_C, x_kgs)

        pipe.downhill_m = -10
        x_kgs = 10
        p_fwd_uphill, t = pipe.perform_calc(p0_bar, t0_C, x_kgs)

        pipe.downhill_m = 0
        x_kgs = -10
        p_back_0, t = pipe.perform_calc(p0_bar, t0_C, x_kgs)

        pipe.downhill_m = 10
        x_kgs = -10
        p_back_downhill, t = pipe.perform_calc(p0_bar, t0_C, x_kgs)

        pipe.downhill_m = -10
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


if __name__ == "__main__":
    unittest.main()
