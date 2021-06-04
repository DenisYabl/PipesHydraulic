from Fluids.HE2_Fluid import HE2_BlackOil
import pandas as pd
from Tests.Optimization_test import gimme_DNS2_inlets_outlets_Q
from Tools.HE2_ABC import oil_params
from Solver.HE2_Solver import HE2_Solver
import shame_on_me
import json
from Tajlaki_DNS2_graph_example import build_DNS2_graph
import networkx as nx
import matplotlib.pyplot as plt
from Tools.HE2_tools import draw_solution, print_wells_pressures

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
oil_params = oil_params(sat_P_bar=67, plastT_C=84, gasFactor=36, oildensity_kg_m3=826,
                        waterdensity_kg_m3=1015, gasdensity_kg_m3=1, oilviscosity_Pa_s=35e-3, volumewater_percent=50, volumeoilcoeff=1.017)

pump_curves = pd.read_csv(r"../../CommonData/PumpChart.csv")

fluid = HE2_BlackOil(oil_params)
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

inclination = {
    '1523': (32.6, 2561, 0.125, 267.3, 202.7, 0.143),
    '146': (47.569, 2575.43, 0.125, 161.39, 210.56, 0.143),
    '142': (324.89, 2511.1, 0.125, 242.24, 249.7599, 0.143),
    '1562':(549.92, 2392.08, 0.125, 36.4499, 325.55, 0.143),
    '1385': (395.8, 2139.2, 0.125, 80, 584, 0.143),
    '736': (216, 2532, 0.125, 448, 155, 0.143),
    '739': (323.45, 2356.55, 0.125, 396.35, 290.6, 0.143),
    '1383': (526, 2353.99, 0.125, 30.16, 377.84, 0.143),
    '738': (243.83, 2128.17, 0.125, 402.45, 526.85, 0.143),
    '725': (495.76, 2329.24, 0.125, 28.7, 384.3, 0.143),
    '731': (285.7, 2414.27, 0.125, 13.5, 330, 0.143),
    '196': (172.6, 2504.43, 0.125, 7.38, 245.12, 0.143),
    '734': (32.21, 2362.8, 0.125, 282, 286.88, 0.143),
    '198': (244, 2586, 0.125, 4.82, 169.18, 0.143),
    '199': (228.11, 2520.89, 0.125, 198.34, 241.56, 0.143),
    '197': (668.11, 2345.89, 0.125, 52.91, 403.09, 0.143),
    '195': (610.82, 2372.96, 0.125, 50.06, 357.93, 0.143),
    '191': (254.37, 2585.63, 0.125, 251.19, 166.221, 0.143),
    '729': (92.78, 2345.22, 0.125, 247.95, 399.048, 0.143),
    '730': (75.76, 2218.24, 0.125, 99.72, 532.78, 0.143),
    '192': (337.39, 2177.61, 0.125, 255.33, 531.87, 0.143),
    '148': (337.39, 2177.61, 0.125, 255.33, 531.87, 0.143),
    '3552': (130.37, 2169.63, 0.125, 155.83, 687.06, 0.143),
    '617': (461.4, 2409.6, 0.125, 8, 333, 0.143),
    '567': (674.06, 2165.94, 0.125, 112.83, 554.17, 0.143),
    '614': (128.89, 2529.07, 0.125, 156.93, 215.51, 0.143),
    '619': (269.27, 2343.73, 0.125, 110.95, 384.75, 0.143),
    '609': (117.42, 2291.58, 0.125, 202.87, 535.23, 0.143),
    '1816': (61.69, 2508.31, 0.125, 811.92, 67.67, 0.143),
    '2630': (419.86, 2370.14, 0.125, 570.43, 310.38, 0.143),
    '1815': (113.77, 2516.23, 0.125, 482.82, 333.82, 0.143),
    '676': (248.29, 2433.71, 0.125, 537.87, 336.76, 0.143),
    '3270': (310.83, 2080.17, 0.125, 315.89, 637.11, 0.143),
    '3266': (119.76, 2500.24, 0.125, 280.05, 489.3, 0.143),
    '1814': (205.63, 2455.37, 0.125, 186, 600.24, 0.143),
    '1817': (324.58, 2443.42, 0.125, 690.26, 244.53, 0.143),
    '4532': (363.75, 2482.25, 0.125, 403.62, 509.5, 0.143),
    '2631': (67.5, 2394.5, 0.125, 626.2, 278.96, 0.143),
    '677': (248.2, 2351.8, 0.125, 728.15, 217.82, 0.143),
    '3113': (382.13, 2367.87, 0.125, 652.87, 230.96, 0.143),
    '3118': (152.43, 2458.57, 0.125, 670.7, 109, 0.143),
    '3112': (55.74, 2481.26, 0.125, 753.03, -3, 0.143),
    '4235': (351.06, 2514.94, 0.125, 511.21, 253.89, 0.143),
    '3117': (117.8, 2477.2, 0.125, 527.58, 255.96, 0.143),
    '1493': (403.96, 2442.04, 0.125, 409.48, 399.746, 0.143),
    '1574': (130.42, 2514.58, 0.125, 679.53, -18.46, 0.143),
    '1579': (425.39, 2494.6, 0.125, 345.25, 281.16, 0.143),
    '3116': (167.57, 2453.43, 0.125, 715.74, 81.184, 0.143)
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
    solver.set_known_Q(inlets_Q)
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
        solver.set_known_Q(inlets_Q)
        solver.solve(threshold=0.05)
        if not solver.op_result.success:
            print(i, solver.op_result.fun)
        succ += solver.op_result.success
        total += 1

    print(f'succ/total: {succ}/{total}')

def test3():
    outputpressure = 4.8
    roughness = 1e-5
    real_diam_coefficient = 0.85

    pumps = {'PAD_33': {'WELL_1385': ['ЭЦН5-125-2500', 41.7], 'WELL_736': ['ЭЦН5-125-2500', 41.7], 'WELL_739': ['ЭЦН5-125-2500', 41.7], 'WELL_1383': ['ЭЦН5-125-2500', 41.7], 'WELL_738': ['ЭЦН5-125-2500', 41.7], 'WELL_725': ['ЭЦН5-125-2500', 41.7]}, 'PAD_34': {'WELL_731': ['ЭЦН5-125-2500', 41.7], 'WELL_196': ['ЭЦН5-125-2500', 41.7], 'WELL_734': ['ЭЦН5-125-2500', 41.7], 'WELL_198': ['ЭЦН5-125-2500', 41.7], 'WELL_199': ['ЭЦН5-125-2500', 41.7], 'WELL_197': ['ЭЦН5-125-2500', 41.7], 'WELL_195': ['ЭЦН5-125-2500', 41.7], 'WELL_191': ['ЭЦН5-125-2500', 41.7], 'WELL_729': ['ЭЦН5-125-2500', 41.7], 'WELL_730': ['ЭЦН5-125-2500', 41.7], 'WELL_192': ['ЭЦН5-125-2500', 41.7], 'WELL_148': ['ЭЦН5-125-2500', 41.7]}}

    # pumps = json.loads(s)

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
    solver.set_known_Q(inlets_Q)
    solver.solve(threshold=0.05)
    print(solver.op_result)
    print_wells_pressures(G, inlets)

def test4():
    outputpressure = 4.8
    roughness = 1e-5
    real_diam_coefficient = 0.85

    pumps = {'PAD_33': {'WELL_1385': ['ЭЦН5-125-2500', 41.7], 'WELL_736': ['ЭЦН5-125-2500', 41.7], 'WELL_739': ['ЭЦН5-125-2500', 41.7], 'WELL_1383': ['ЭЦН5-125-2500', 41.7], 'WELL_738': ['ЭЦН5-125-2500', 41.7], 'WELL_725': ['ЭЦН5-125-2500', 41.7]}, 'PAD_34': {'WELL_731': ['ЭЦН5-125-2500', 41.7], 'WELL_196': ['ЭЦН5-125-2500', 41.7], 'WELL_734': ['ЭЦН5-125-2500', 41.7], 'WELL_198': ['ЭЦН5-125-2500', 41.7], 'WELL_199': ['ЭЦН5-125-2500', 41.7], 'WELL_197': ['ЭЦН5-125-2500', 41.7], 'WELL_195': ['ЭЦН5-125-2500', 41.7], 'WELL_191': ['ЭЦН5-125-2500', 41.7], 'WELL_729': ['ЭЦН5-125-2500', 41.7], 'WELL_730': ['ЭЦН5-125-2500', 41.7], 'WELL_192': ['ЭЦН5-125-2500', 41.7], 'WELL_148': ['ЭЦН5-125-2500', 41.7]}}

    # pumps = json.loads(s)

    graph_params = dict(
        pressures=pressures,
        plasts=plasts,
        pumps=pumps,
        pump_curves=pump_curves,
        fluid=fluid,
        roughness=roughness,
        real_diam_coefficient=real_diam_coefficient,
        inclination=inclination
    )

    G, inlets, juncs, outlets = build_DNS2_graph(PAD_33='P=5.9', PAD_34='Q=12.7', **graph_params)
    solver = HE2_Solver(G)
    # inlets_Q = gimme_DNS2_inlets_outlets_Q()
    # solver.prepare_initial_approximation(G, inlets_Q)
    solver.solve(threshold=0.05, it_limit=10)
    print(solver.op_result)
    print_wells_pressures(G, inlets)


def test5():
    outputpressure = 4.8
    roughness = 1e-5
    real_diam_coefficient = 0.85

    pumps = {'PAD_33': {'WELL_1385': ['ЭЦН5-125-2500', 41.7], 'WELL_736': ['ЭЦН5-125-2500', 41.7], 'WELL_739': ['ЭЦН5-125-2500', 41.7], 'WELL_1383': ['ЭЦН5-125-2500', 41.7], 'WELL_738': ['ЭЦН5-125-2500', 41.7], 'WELL_725': ['ЭЦН5-125-2500', 41.7]}, 'PAD_34': {'WELL_731': ['ЭЦН5-125-2500', 41.7], 'WELL_196': ['ЭЦН5-125-2500', 41.7], 'WELL_734': ['ЭЦН5-125-2500', 41.7], 'WELL_198': ['ЭЦН5-125-2500', 41.7], 'WELL_199': ['ЭЦН5-125-2500', 41.7], 'WELL_197': ['ЭЦН5-125-2500', 41.7], 'WELL_195': ['ЭЦН5-125-2500', 41.7], 'WELL_191': ['ЭЦН5-125-2500', 41.7], 'WELL_729': ['ЭЦН5-125-2500', 41.7], 'WELL_730': ['ЭЦН5-125-2500', 41.7], 'WELL_192': ['ЭЦН5-125-2500', 41.7], 'WELL_148': ['ЭЦН5-125-2500', 41.7]}}

    # pumps = json.loads(s)

    graph_params = dict(
        pressures=pressures,
        plasts=plasts,
        pumps=pumps,
        pump_curves=pump_curves,
        fluid=fluid,
        roughness=roughness,
        real_diam_coefficient=real_diam_coefficient,
        inclination=inclination
    )

    G, inlets, juncs, outlets = build_DNS2_graph(PAD_33='P=5.9', PAD_34='Q=12.7', **graph_params)
    solver = HE2_Solver(G)
    solver.solve(threshold=0.05, it_limit=100)
    draw_solution(G, None, None, inlets, outlets, juncs)
    print_wells_pressures(G, inlets)


if __name__ == '__main__':
    test5()