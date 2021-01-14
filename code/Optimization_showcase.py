from HE2_Fluid import HE2_OilWater
from Optimization_test import model_DNS_2
oil_params = {
    "OilSaturationP": 65.7, #Давление насыщения нефти при стандартных условиях, исходные данные
    "PlastT": 84.0, #Пластовая температура, исходные данные
    "GasFactor": 34.8, #Газовый фактор нефти, исходные данные
    "SepOilWeight": 0.8815, #Плотность нефти, исходные данные
    "GasDensity": 1.003, #Плотность попутного газа, исходные данные
    "SepOilDynamicViscosity": 43.03, #Динамическая вязкость нефти, исходные данные
    "wellopVolumeWater": 56, #Обводненность нефти, исходные данные
    "VolumeOilCoeff": 1.097, #Объемный коэффициент нефти, исходные данные
    "PlastWaterWeight": 1.015, #Плотность попутной воды, исходные данные
    "adkuLiquidDebit": 240, #Дебит скважины, исходные данные
    "CurrentP": 90, #Текущее давление, меняем по необходимости
    "CurrentT": 84.0, #Текущая температура, меняем по необходисости
}

fluid = HE2_OilWater(oil_params)
pressures = {"PAD_5": 20, "PAD_33": 19, "PAD_34": 18, "PAD_39": 17, "PAD_49": 16, "PAD_57": 15}
G = model_DNS_2(pressures=pressures, daily_debit=8000, fluid=fluid)
print(G.nodes["DNS_2"]["obj"].result)