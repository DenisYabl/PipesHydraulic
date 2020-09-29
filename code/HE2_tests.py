import unittest
from HE2_Pipe import HE2_WaterPipeSegment

class TestWaterPipeSegment(unittest.TestCase):
    def setUp(self):
        pipe = HE2_WaterPipeSegment()
        pipe.inner_diam_m = 0.1
        pipe.roughness_m = 1e-5
        pipe.L_m = 100
        pipe.delta_H_m = 0
        self.pipe = pipe

    def test_1(self):
        p, t = self.pipe.perform_calc(50, 20, 1)
        self.assertAlmostEqual(p, 49.99978631383106)
        self.assertAlmostEqual(t, 20)

    def test_2(self):
        p, t = self.pipe.perform_calc(50, 20, -1)
        self.assertAlmostEqual(p, 50.0003)
        self.assertAlmostEqual(t, 20)


if __name__ == "__main__":
    unittest.main()