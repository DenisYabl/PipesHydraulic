import unittest
from HE2_Pipe import HE2_WaterPipeSegment
import uniflocpy.uWell.uPipe as uPipe
import uniflocpy.uTools.uconst as uc
import Hydraulics.Methodics.Mukherjee_Brill as mb
import Hydraulics.Properties.Mishenko as msch

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
        p0_bar = 50
        t0_C = 20
        x_kgs = 1
        pipe = HE2_WaterPipeSegment()
        pipe.inner_diam_m = 0.1
        pipe.roughness_m = 1e-5
        pipe.L_m = 100
        pipe.slope_m = 0
        p, t = pipe.perform_calc(p0_bar, t0_C, x_kgs)
        print(uc.bar2Pa(p-p0_bar), t)

        Q_m3day = x_kgs * 86400 / pipe.fluid.rho_wat_kgm3

        oil_params = dict(CurrentP=p0_bar, PlastT=t0_C, CurrentT=t0_C, adkuLiquidDebit=Q_m3day)
        oil_params.update(wellopVolumeWater=100, GasFactor=1e-8, PlastWaterWeight=pipe.fluid.rho_wat_kgm3/1000, wellopLiquidDensityTm=pipe.fluid.rho_wat_kgm3/1000, adkuOilDebit=0)
        oil_params.update(VolumeOilCoeff=1.2, PlastOilDynamicViscosity=1,SepOilDynamicViscosity=1, OilSaturationP=100, PlastOilWeight=0.85, SepOilWeight=0.85, GasDensity=0.8)
        tubing = dict(IntDiameter=pipe.inner_diam_m, angle=90, Roughness=pipe.roughness_m)
        mishenko = msch.Mishenko.from_oil_params(oil_params, tubing)
        misch_grad = mb.calculate(mishenko, tubing)
        misch_dp = misch_grad * pipe.L_m
        print(misch_dp)
        pass



if __name__ == "__main__":
    unittest.main()