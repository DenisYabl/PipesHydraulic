from Hydraulics.Properties.Mishenko import Mishenko, from_oil_params
from Tools.HE2_ABC import HE2_ABC_Fluid
from Fluids.oil_params import oil_params, dummy_oil_params

class HE2_DummyWater(HE2_ABC_Fluid):
    def __init__(self):
        self.rho_wat_kgm3 = 1000
        self.mu_wat_cp = 1

    def calc(self, P_bar, T_C):
        pass


class HE2_OilWater(HE2_ABC_Fluid):
    def __init__(self, oil_params):
        # check_for_nan(**oil_params)
        self.oil_params = oil_params
        self.CurrentLiquidDensity_kg_m3 = (self.oil_params.oildensity_kg_m3 * (1 - self.oil_params.volumewater_percent / 100) +
                                                      self.oil_params.waterdensity_kg_m3 * self.oil_params.volumewater_percent / 100)


    def calc(self, P_bar, T_C, X_kgsec, IntDiameter=None):
        P_for_PVT = max(abs(P_bar), 0.75)
        calc_params = self.oil_params

        tubing = {"IntDiameter": IntDiameter} if IntDiameter else None
        temp_mishenko = from_oil_params(P_for_PVT, T_C, None, calc_params=calc_params, tubing=tubing)
        #Side effects
        self.CurrentLiquidDensity_kg_m3 = temp_mishenko.CurrentLiquidDensity_kg_m3
        self.CurrentOilViscosity_Pa_s = temp_mishenko.CurrentOilViscosity_Pa_s
        #Return for pressure gradient calculation
        return temp_mishenko

class HE2_DummyOil(HE2_ABC_Fluid):
    def __init__(self, daily_Q=100, VolumeWater=50):
        self.oil_params = dummy_oil_params(Q_m3_day=daily_Q, volumeWater=VolumeWater)
        self.CurrentLiquidDensity_kg_m3 = (self.oil_params.oildensity_kg_m3 * (1 - self.oil_params.volumewater_percent / 100) +
                                                      self.oil_params.waterdensity_kg_m3 * self.oil_params.volumewater_percent / 100)


    def calc(self, P_bar, T_C, X_kgsec, IntDiameter=None):
        P_for_PVT = max(abs(P_bar), 0.75)
        calc_params = self.oil_params

        tubing = {"IntDiameter": IntDiameter} if IntDiameter else None
        temp_mishenko = from_oil_params(P_for_PVT, T_C, X_kgsec, calc_params=calc_params, tubing=tubing)
        #Side effects
        self.CurrentLiquidDensity_kg_m3 = temp_mishenko.CurrentLiquidDensity_kg_m3
        self.CurrentOilViscosity_Pa_s = temp_mishenko.CurrentOilViscosity_Pa_s
        #Return for pressure gradient calculation
        return temp_mishenko


if __name__ == '__main__':
    fl = HE2_DummyWater()
    fl.calc(100, 100)
    print(fl.rho_wat_kgm3)