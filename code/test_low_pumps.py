from Fluids.HE2_Fluid import HE2_OilWater
import pandas as pd
from Optimization_test import model_DNS_2, build_DNS2_graph, gimme_DNS2_inlets_outlets_Q
from Fluids.oil_params import oil_params
from Solver.HE2_Solver import HE2_Solver
import shame_on_me
import json

pumps_json_str = r'{"PAD_33": {"WELL_1385": ["\u042d\u0426\u041d5\u0410-320-2400", 42.1], "WELL_736": ["\u042d\u0426\u041d5\u0410-250-2600", 43.7], "WELL_739": ["\u042d\u0426\u041d5\u0410-250-2450", 41.2], "WELL_1383": ["\u042d\u0426\u041d5\u0410-250-2300", 43.0], "WELL_738": ["\u042d\u0426\u041d5\u0410-250-2350", 40.8], "WELL_725": ["\u042d\u0426\u041d5\u0410-800-2200", 42.2]}, "PAD_34": {"WELL_731": ["\u042d\u0426\u041d5-125-2750", 41.1], "WELL_196": ["\u042d\u0426\u041d5\u0410-250-2600", 40.8], "WELL_734": ["\u042d\u0426\u041d5\u0410-320-2550", 43.4], "WELL_198": ["\u042d\u0426\u041d5\u0410-800-2100", 41.0], "WELL_199": ["\u042d\u0426\u041d5-125-2550", 44.6], "WELL_197": ["\u042d\u0426\u041d5-125-2750", 43.3], "WELL_195": ["\u042d\u0426\u041d5\u0410-250-2400", 45.1], "WELL_191": ["\u042d\u0426\u041d5\u0410-800-2100", 43.2], "WELL_729": ["\u042d\u0426\u041d5-125-2550", 41.7], "WELL_730": ["\u042d\u0426\u041d5\u0410-400-2400", 46.2], "WELL_192": ["\u042d\u0426\u041d5-200-2500", 41.9], "WELL_148": ["\u042d\u0426\u041d5\u0410-800-2200", 43.6]}}'
new_pumps = json.loads(pumps_json_str)

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
oil_params = oil_params(Q_m3_day=500, sat_P_bar=67, plastT_C=84, gasFactor=36, oildensity_kg_m3=826,
                        waterdensity_kg_m3=1015, gasdensity_kg_m3=1, oilviscosity_Pa_s=35e-3, volumewater_percent=50, volumeoilcoeff=1.017)

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


def test1():
    outputpressure = 4.8
    roughness = 1e-5
    real_diam_coefficient = 0.85

    graph_params = dict(
        pressures=pressures,
        plasts=plasts,
        pumps=new_pumps,
        pump_curves=pump_curves,
        fluid=fluid,
        roughness=roughness,
        real_diam_coefficient=real_diam_coefficient,
        DNS_pressure=outputpressure
    )

    G, inlets, juncs, outlets = shame_on_me.build_DNS2_graph_pads33_34(**graph_params)
    solver = HE2_Solver(G)
    inlets_Q = gimme_DNS2_inlets_outlets_Q()
    solver.prepare_initial_approximation(G, inlets_Q)
    solver.solve()

    for n in inlets:
        print(n, G.nodes[n]["obj"].result)

    for n in outlets:
        print(n, G.nodes[n]["obj"].result)

import json
def test2():
    outputpressure = 4.8
    roughness = 1e-5
    real_diam_coefficient = 0.85

    fn = r'pumps.txt'
    f = open(fn, 'r')
    succ, total = 0, 0
    l = list(f)
    for i, s in enumerate(l):
    # for i in [519, 665, 689]:
    #     s = l[i]

        if i%10 == 0:
            print(f'{i}/{len(l)}')
        pumps = json.loads(s)

        graph_params = dict(
            pressures=pressures,
            plasts=plasts,
            pumps=pumps,
            pump_curves=pump_curves,
            fluid=fluid,
            roughness=roughness,
            real_diam_coefficient=real_diam_coefficient,
            DNS_pressure=outputpressure
        )

        G, inlets, juncs, outlets = shame_on_me.build_DNS2_graph_pads33_34(**graph_params)
        solver = HE2_Solver(G)
        inlets_Q = gimme_DNS2_inlets_outlets_Q()
        solver.prepare_initial_approximation(G, inlets_Q)
        solver.solve(threshold=0.05)
        if not solver.op_result.success:
            print(i, solver.op_result.fun)
        succ += solver.op_result.success
        total += 1

    print(f'succ/total: {succ}/{total}')


if __name__ == '__main__':
    test2()