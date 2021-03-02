from Fluids.HE2_Fluid import HE2_OilWater
from Optimization_test import model_DNS_2, build_DNS2_graph, gimme_DNS2_inlets_outlets_Q
from Optimization_test import model_DNS_2_by_parts
import pandas as pd

from Solver.HE2_Solver import HE2_Solver

"""
oil_params - описание ФХС водонефтяной смеси, но данном этапе - усредненные ФХС Тайлаковского месторождения
pressures - перечень пластовых давлений всех используемых в модели скважин

plasts - описание свойств пласта используемых в модели скважин. В настоящий момент задается коэффициентом продуктивности.
В датайрейме для оптимизации он в поле productCoeff.
Коэффициент продуктивности можно вычислить по историческим данным (глубина данных 1-2 месяца должна быть оптимальной) 
как efficiency = Debit / (P_plast - P_zaboi),
где Debit - суточный дебит скважины, P_plast - пластовое давление, бар, P_zaboi - забойное давление, бар

pumps - описание насосов скважин, задается моделью насоса и частотой работы
"""
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

pump_curves = pd.read_csv("../CommonData/PumpChart.csv")

fluid = HE2_OilWater(oil_params)
pressures = {"PAD_5": {"WELL_1523" : 270.5, "WELL_146" : 270.5, "WELL_142" : 270.5, "WELL_1562" : 268.5},
             "PAD_33":  {"WELL_1385" : 268.5, "WELL_736" : 268.5, "WELL_739" : 270.5, "WELL_1383" : 270.5, "WELL_738" : 268.5,"WELL_725" : 268.5},
             "PAD_34": {"WELL_731" : 270.5, "WELL_196" : 268.5, "WELL_734" : 270.5, "WELL_198" : 270.5, "WELL_199" : 270.5,"WELL_197" : 268.5,
                        "WELL_195" : 268.5, "WELL_191" : 268.5, "WELL_729" : 268.5, "WELL_730" : 268.5, "WELL_192" : 268.5,"WELL_148" : 270.5},
             "PAD_39": {"WELL_3552" : 270.5, "WELL_617" : 270.5, "WELL_567" : 270.5, "WELL_614" : 270.5, "WELL_619" : 270.5,"WELL_609" : 268.5},
             "PAD_49":  {"WELL_1816" : 268.5, "WELL_2630" : 268.5, "WELL_1815" : 268.5, "WELL_676" : 268.5, "WELL_3270" : 268.5,
                         "WELL_3266" : 268.5, "WELL_1814" : 268.5, "WELL_1817" : 268.5, "WELL_4532" : 268.5, "WELL_2631" : 268.5,
                         "WELL_677" : 268.5},
             "PAD_57": {"WELL_3113" : 268.5, "WELL_3118" : 268.5, "WELL_3112" : 268.5, "WELL_4235" : 268.5, "WELL_3117" : 268.5,
                         "WELL_1493" : 268.5, "WELL_1574" : 268.5, "WELL_1579" : 268.5, "WELL_3116" : 268.5}}


plasts = {"PAD_5": {"WELL_1523":0.237, "WELL_146": 0.416, "WELL_142": 0.158,"WELL_1562": 0.276},
"PAD_33": {"WELL_1385" : 0.45, "WELL_736" : 0.3, "WELL_739" : 0.378, "WELL_1383" : 0.5, "WELL_738" : 0.33,"WELL_725" : 0.5},
"PAD_34": {"WELL_731" : 0.347, "WELL_196" : 0.6, "WELL_734" : 0.6, "WELL_198" : 0.4, "WELL_199" : 0.457,"WELL_197" : 0.236,
    "WELL_195" : 0.27, "WELL_191" : 0.227, "WELL_729" : 1.4, "WELL_730" : 0.189, "WELL_192" : 0.445,"WELL_148" : 0.232},
"PAD_39": {"WELL_3552": 0.98, "WELL_617": 0.348, "WELL_567": 0.034, "WELL_614": 0.98, "WELL_619": 0.31,"WELL_609": 0.235},
"PAD_49":  {"WELL_1816" : 1.66, "WELL_2630" : 0.05, "WELL_1815" : 0.516, "WELL_676" : 0.823, "WELL_3270" : 0.288,
    "WELL_3266" : 0.492, "WELL_1814" : 1.259, "WELL_1817" : 0.695, "WELL_4532" : 0.293, "WELL_2631" : 0.68, "WELL_677" : 0.37},
"PAD_57": {"WELL_3113" : 0.515, "WELL_3118" : 0.052, "WELL_3112" : 0.5, "WELL_4235" : 0.17, "WELL_3117" : 1.29,
    "WELL_1493" : 0.792, "WELL_1574" : 0.53, "WELL_1579" : 0.197, "WELL_3116" : 0.98}
}


pumps = {"PAD_5": {"WELL_1523":["ЭЦН5-80-2500", 50], "WELL_146": ["ЭЦН5-60-2450", 50], "WELL_142": ["ЭЦН5-45-2500", 50],"WELL_1562": ["ЭЦН5-50-2500", 50]},
"PAD_33": {"WELL_1385" : ["ЭЦН5-30-2600", 50], "WELL_736" : ["ЭЦН5-30-2600", 50], "WELL_739" : ["ЭЦН5-80-2200", 50],
    "WELL_1383" : ["ЭЦН5-80-2550", 50], "WELL_738" : ["ЭЦН5-80-2500", 50],"WELL_725" : ["ЭЦН5-80-2500", 50]},

"PAD_34": {"WELL_731" : ["ЭЦН5-80-2550", 50], "WELL_196" : ["ЭЦН5-125-2500", 50], "WELL_734" : ["ЭЦН5-80-2200", 50],
    "WELL_198" :["ЭЦН5-80-2550", 50], "WELL_199" : ["ЭЦН5-80-2600", 50],"WELL_197" : ["ЭЦН5-80-2550", 50],
    "WELL_195" : ["ЭЦН5-60-2500", 50], "WELL_191" : ["ЭЦН5-60-2450", 50], "WELL_729" :["ЭЦН5А-250-2400", 50],
    "WELL_730" : ["ЭЦН5-45-2400", 50], "WELL_192" : ["ЭЦН5-80-2650", 50],"WELL_148" : ["ЭЦН5-60-2450", 50]},

"PAD_39": {"WELL_3552": ["ЭЦН5А-160-2400", 50], "WELL_617": ["ЭЦН5-80-2500", 50], "WELL_567": ["ЭЦН5-50-2450", 50],
    "WELL_614": ["ЭЦН5А-160-2500", 50], "WELL_619": ["ЭЦН5-80-2550", 50],"WELL_609": ["ЭЦН5-80-2400", 50]},

"PAD_49":  {"WELL_1816" : ["ЭЦН5А-320-2400", 50], "WELL_2630" : ["ЭЦН5-80-2400", 50], "WELL_1815" : ["ЭЦН5А-160-2500", 50],
            "WELL_676" : ["ЭЦН5А-160-2600", 50], "WELL_3270" : ["ЭЦН5-45-2500", 50], "WELL_3266" : ["ЭЦН5-80-2700", 50],
            "WELL_1814" : ["ЭЦН5А-250-2400", 50], "WELL_1817" : ["ЭЦН5А-160-2600", 50], "WELL_4532" : ["ЭЦН5-60-2400", 50],
            "WELL_2631" : ["ЭЦН5А-160-2400", 50], "WELL_677" : ["ЭЦН5-80-2600", 50]},

"PAD_57": {"WELL_3113" : ["ЭЦН5-125-2500", 50], "WELL_3118" : ["ЭЦН5-80-2500", 50], "WELL_3112" : ["ЭЦН5-80-2400", 50],
           "WELL_4235" : ["ЭЦН5-80-2500", 50], "WELL_3117" : ["ЭЦН5А-250-2500", 50], "WELL_1493" : ["ЭЦН5А-160-2500", 50],
           "WELL_1574" : ["ЭЦН5-125-2450", 50], "WELL_1579" : ["ЭЦН5-80-2500", 50], "WELL_3116" : ["ЭЦН5А-250-2400", 50]}
}


def full_test():
    Full_system_daily_debit = 5200
    G, inlets, juncs, outlets = model_DNS_2(pressures=pressures, pumps=pumps, plasts=plasts,
                                            daily_debit=Full_system_daily_debit, pump_curves=pump_curves, fluid=fluid)
    for n in inlets:
        print(n, G.nodes[n]["obj"].result)

    for n in outlets:
        print(n, G.nodes[n]["obj"].result)


def part_test():
    Full_system_daily_debit = 4750
    # well_list = ['PAD_49_well_1816']
    model_DNS_2_by_parts(pressures=pressures, pumps = pumps, plasts= plasts, daily_debit=Full_system_daily_debit, pump_curves=pump_curves, fluid=fluid)

def test_with_change_graph_ont_the_fly():
    # Выходное условие заменено на выходное давление ДНС-2
    outputpressure = 4.8
    roughness = 1e-5
    real_diam_coefficient = 0.85

    #Получаем расчетный граф с заданными параметрами
    G, inlets, juncs, outlets = build_DNS2_graph(pressures = pressures, plasts = plasts, pumps = pumps, pump_curves = pump_curves, fluid = fluid, roughness = roughness, real_diam_coefficient = real_diam_coefficient,
                                                 DNS_pressure = outputpressure)
    # Создаем солвер и решаем полученный расчетный граф
    solver = HE2_Solver(G)
    solver.solve()

    for n in inlets:
        print(n, G.nodes[n]["obj"].result)

    for n in outlets:
        print(n, G.nodes[n]["obj"].result)

    #Меняем давление на пласте скважины 1523 куста 5

    G.nodes["PAD_5_well_1523"]['obj'].value = 250


    solver.solve()
    print("\nРешение с измененным давлением")
    for n in inlets:
        print(n, G.nodes[n]["obj"].result)

    for n in outlets:
        print(n, G.nodes[n]["obj"].result)

    #Меняем частоту насоса скважины 1562 куста 5

    G.edges._adjdict["Pump_intake_1562"]['Pump_outlet_1562']['obj'].changeFrequency(55)

    solver.solve()
    print("\nРешение с измененной частотой насоса")
    for n in inlets:
        print(n, G.nodes[n]["obj"].result)

    for n in outlets:
        print(n, G.nodes[n]["obj"].result)


if __name__ == '__main__':
    test_with_change_graph_ont_the_fly()
    #full_test
    #part_test()
