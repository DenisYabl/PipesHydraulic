# oil_params = {
#     "OilSaturationP": 66.7, #Давление насыщения нефти при стандартных условиях, исходные данные
#     "PlastT": 84.0, #Пластовая температура, исходные данные
#     "GasFactor": 34.8, #Газовый фактор нефти, исходные данные
#     "SepOilWeight": 0.885, #Плотность нефти, исходные данные
#     "GasDensity": 1.003, #Плотность попутного газа, исходные данные
#     "SepOilDynamicViscosity": 43.03, #Динамическая вязкость нефти, исходные данные
#     "wellopVolumeWater": 59.7, #Обводненность нефти, исходные данные
#     "VolumeOilCoeff": 1.097, #Объемный коэффициент нефти, исходные данные
#     "PlastWaterWeight": 1.015, #Плотность попутной воды, исходные данные
# }

# TODO make oil_params checkable for nans

fieldlist = ['Q_m3_sec', 'Q_m3_day', 'sat_P_bar', 'plastT_C', 'gasFactor', 'oildensity_kg_m3', 'waterdensity_kg_m3', 'gasdensity_kg_m3',
             'oilviscosity_Pa_s', 'volumewater_percent', 'volumeoilcoeff', 'currentP_bar', 'currentT_C', 'CurrentLiquidDensity_kg_m3']

class oil_params():
    def __init__(self, Q_m3_day=None, sat_P_bar=None, plastT_C=None, gasFactor=None, oildensity_kg_m3=None,
                 waterdensity_kg_m3=None, gasdensity_kg_m3=None, oilviscosity_Pa_s=None, volumewater_percent=None, volumeoilcoeff=None):
        self.Q_m3_day = Q_m3_day
        self.Q_m3_sec = self.Q_m3_day / 86400
        self.sat_P_bar = sat_P_bar
        self.plastT_C = plastT_C
        self.gasFactor = gasFactor
        self.oildensity_kg_m3 = oildensity_kg_m3
        self.waterdensity_kg_m3 = waterdensity_kg_m3
        self.gasdensity_kg_m3 = gasdensity_kg_m3
        self.oilviscosity_Pa_s = oilviscosity_Pa_s
        self.volumewater_percent = volumewater_percent
        self.volumeoilcoeff = volumeoilcoeff
        self.currentP_bar = None
        self.currentT_C = None
        self.CurrentLiquidDensity_kg_m3 = None

class dummy_oil_params:
    def __init__(self, Q_m3_day=None, volumeWater=50):
        self.Q_m3_day = Q_m3_day
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
        self.CurrentLiquidDensity_kg_m3 = None
