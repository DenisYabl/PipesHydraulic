from Hydraulics.Properties.Mishenko import Mishenko
from Tools.HE2_ABC import HE2_ABC_Fluid

class HE2_DummyWater(HE2_ABC_Fluid):
    def __init__(self):
        self.rho_wat_kgm3 = 1000
        self.mu_wat_cp = 1

    def calc(self, P_bar, T_C):
        pass


class HE2_OilWater(HE2_ABC_Fluid):
    def __init__(self, oil_params):
        #Давление насыщения
        self.SaturationPressure = oil_params["OilSaturationP"]
        # Пластовая температура
        self.PlastT = oil_params["PlastT"]
        # Газовый фактор нефти
        self.GasFactor = oil_params["GasFactor"]
        # Доля углеводородных газов
        self.CarbGasAmount = 0.91  # Нет в данных
        # Доля неуглеводородных газов
        self.NonCarbGasAmount = 0.09  # Нет в данных
        # Плотность дегазированной нефти
        self.SepOilDensity = oil_params["SepOilWeight"]
        # Плотность газа
        self.GasDensity = oil_params["GasDensity"]
        # Динамическая вязкость сепарированной нефти
        self.SepOilDynamicViscosity = oil_params["SepOilDynamicViscosity"]
        # Обводненность нефти
        self.VolumeWater = oil_params["wellopVolumeWater"]
        # Объемный фактор нефти
        self.OilVolumeCoeff = oil_params["VolumeOilCoeff"]
        # Плотность пластовой воды
        self.PlastWaterDensity = oil_params["PlastWaterWeight"]
        self.CurrentLiquidDensity = 1000 * (self.SepOilDensity * (1 - self.VolumeWater / 100) + self.PlastWaterDensity * self.VolumeWater / 100)



    def calc(self, P_bar, T_C, Q_Liquid, IntDiameter):
        P_for_PVT = max(abs(P_bar), 0.75)
        oil_params = {
            "OilSaturationP": self.SaturationPressure,  # Давление насыщения нефти при стандартных условиях, исходные данные
            "PlastT": self.PlastT,  # Пластовая температура, исходные данные
            "GasFactor": self.GasFactor,  # Газовый фактор нефти, исходные данные
            "SepOilWeight": self.SepOilDensity,  # Плотность нефти, исходные данные
            "GasDensity": self.GasDensity,  # Плотность попутного газа, исходные данные
            "SepOilDynamicViscosity": self.SepOilDynamicViscosity,  # Динамическая вязкость нефти, исходные данные
            "wellopVolumeWater": self.VolumeWater,  # Обводненность нефти, исходные данные
            "VolumeOilCoeff": self.OilVolumeCoeff,  # Объемный коэффициент нефти, исходные данные
            "PlastWaterWeight": self.PlastWaterDensity,  # Плотность попутной воды, исходные данные
            "adkuLiquidDebit": Q_Liquid,  # Дебит скважины, исходные данные
            "CurrentP": P_for_PVT,
            "CurrentT": T_C
        }
        tubing = {"IntDiameter": IntDiameter}
        temp_mishenko = Mishenko.from_oil_params(oil_params=oil_params, tubing=tubing)
        #Side effects
        self. CurrentLiquidDensity = temp_mishenko.CurrentLiquidDensity
        self.CurrentOilViscosity = temp_mishenko.CurrentOilViscosity
        self.Q_m3sec = temp_mishenko.Q
        #Return for pressure gradient calculation
        return temp_mishenko

class HE2_DummyOil(HE2_ABC_Fluid):
    def __init__(self, VolumeWater):
        #Давление насыщения
        self.SaturationPressure = 67
        # Пластовая температура
        self.PlastT = 84
        # Газовый фактор нефти
        self.GasFactor = 39
        # Доля углеводородных газов
        self.CarbGasAmount = 0.91  # Нет в данных
        # Доля неуглеводородных газов
        self.NonCarbGasAmount = 0.09  # Нет в данных
        # Плотность дегазированной нефти
        self.SepOilDensity = 0.826
        # Плотность газа
        self.GasDensity = 1
        # Динамическая вязкость сепарированной нефти
        self.SepOilDynamicViscosity = 35
        # Обводненность нефти
        self.VolumeWater = VolumeWater
        # Объемный фактор нефти
        self.OilVolumeCoeff = 1.15
        # Плотность пластовой воды
        self.PlastWaterDensity = 1.015
        self.CurrentLiquidDensity = 1000 * (self.SepOilDensity * (1 - self.VolumeWater / 100) + self.PlastWaterDensity * self.VolumeWater / 100)

    def calc(self, P_bar, T_C, X_kgsec, IntDiameter):
        P_for_PVT = abs(P_bar)
        oil_params = {
            "OilSaturationP": self.SaturationPressure,  # Давление насыщения нефти при стандартных условиях, исходные данные
            "PlastT": self.PlastT,  # Пластовая температура, исходные данные
            "GasFactor": self.GasFactor,  # Газовый фактор нефти, исходные данные
            "SepOilWeight": self.SepOilDensity,  # Плотность нефти, исходные данные
            "GasDensity": self.GasDensity,  # Плотность попутного газа, исходные данные
            "SepOilDynamicViscosity": self.SepOilDynamicViscosity,  # Динамическая вязкость нефти, исходные данные
            "wellopVolumeWater": self.VolumeWater,  # Обводненность нефти, исходные данные
            "VolumeOilCoeff": self.OilVolumeCoeff,  # Объемный коэффициент нефти, исходные данные
            "PlastWaterWeight": self.PlastWaterDensity,  # Плотность попутной воды, исходные данные
            "adkuLiquidDebit": X_kgsec,  # Дебит скважины, исходные данные
            "CurrentP": P_for_PVT,
            "CurrentT": T_C
        }
        tubing = {"IntDiameter": IntDiameter}
        temp_mishenko = Mishenko.from_oil_params(oil_params=oil_params, tubing=tubing)
        #Side effects
        self. CurrentLiquidDensity = temp_mishenko.CurrentLiquidDensity
        self.CurrentOilViscosity = temp_mishenko.CurrentOilViscosity
        self.Q_m3sec = temp_mishenko.Q
        #Return for pressure gradient calculation
        return temp_mishenko


if __name__ == '__main__':
    fl = HE2_DummyWater()
    fl.calc(100, 100)
    print(fl.rho_wat_kgm3)