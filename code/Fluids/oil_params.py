oil_params = {
    "OilSaturationP": 66.7, #Давление насыщения нефти при стандартных условиях, исходные данные
    "PlastT": 84.0, #Пластовая температура, исходные данные
    "GasFactor": 34.8, #Газовый фактор нефти, исходные данные
    "SepOilWeight": 0.885, #Плотность нефти, исходные данные
    "GasDensity": 1.003, #Плотность попутного газа, исходные данные
    "SepOilDynamicViscosity": 43.03, #Динамическая вязкость нефти, исходные данные
    "wellopVolumeWater": 59.7, #Обводненность нефти, исходные данные
    "VolumeOilCoeff": 1.097, #Объемный коэффициент нефти, исходные данные
    "PlastWaterWeight": 1.015, #Плотность попутной воды, исходные данные
}


class oil_params():
    def __init__(self, dailyQ=None, saturationPressure=None, plastT=None, gasFactor=None, oilDensity=None,
                 waterDensity=None, gasDensity=None, oilViscosity=None, volumeWater=None, volumeoilcoeff=None):
        self.Q_m3_day = dailyQ
        self.Q_m3_sec = self.Q_m3_day / 86400
        self.sat_P_bar = saturationPressure
        self.plastT_C = plastT
        self.gasFactor = gasFactor
        self.oildensity_kg_m3 = oilDensity
        self.waterdensity_kg_m3 = waterDensity
        self.gasdensity_kg_m3 = gasDensity
        self.oilviscosity_Pa_s = oilViscosity
        self.volumewater_percent = volumeWater
        self.volumeoilcoeff = volumeoilcoeff
        self.currentP_bar = None
        self.currentT_C = None

class dummy_oil_params:
    def __init__(self, dailyQ=None, volumeWater=50):
        self.Q_m3_day = dailyQ
        self.Q_m3_sec = self.Q_m3_day / 86400
        self.sat_P_bar = 66.7
        self.plastT_C = 84
        self.gasFactor = 39
        self.oildensity_kg_m3 = 826
        self.waterdensity_kg_m3 = 1015
        self.gasdensity_kg_m3 = 1
        self.oilviscosity_Pa_s = 35e-3
        self.volumewater_percent = volumeWater
        self.volumeoilcoeff = 1.015
        self.currentP_bar = None
        self.currentT_C = None