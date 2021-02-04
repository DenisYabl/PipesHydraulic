import networkx as nx

from HE2_Fluid import HE2_OilWater
import HE2_Vertices as vrtxs
import pandas as pd

from HE2_Pipe import HE2_OilPipe
from HE2_Solver import HE2_Solver

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
"""
inlets = dict(PAD_36=vrtxs.HE2_Source_Vertex('Q', 116 * 826 / 86400, fluid, 20),
            PAD_5=vrtxs.HE2_Source_Vertex('Q', 301.2 * 826 / 86400, fluid, 20),
            PAD_35=vrtxs.HE2_Source_Vertex('Q', 160 * 826 / 86400, fluid, 20),
            PAD_155=vrtxs.HE2_Source_Vertex('P', 19.5, fluid, 20))

#PAD_5=vrtxs.HE2_Source_Vertex('P', 19.5, fluid, 20),


juncs = dict(junc_2=vrtxs.HE2_ABC_GraphVertex(),
            junc_4=vrtxs.HE2_ABC_GraphVertex(),
             junc_7=vrtxs.HE2_ABC_GraphVertex(),
             junc_8=vrtxs.HE2_ABC_GraphVertex(),
             junc_9=vrtxs.HE2_ABC_GraphVertex())

outlets = dict(DNS_2 = vrtxs.HE2_Boundary_Vertex('Q', 3 * 767.2 * 826 / 86400))

G = nx.DiGraph()  # Di = directed
for k, v in {**inlets, **outlets, **juncs}.items():
    G.add_node(k, obj=v)

roughness = 3.5
real_diam_coefficient = 0.85

G.add_edge('PAD_36', 'junc_2', obj=HE2_OilPipe([2308], [7], [0.139 * real_diam_coefficient], [roughness]))  #+
#G.add_edge('junc_2', 'junc_3', obj=HE2_OilPipe([250], [0], [0.199 * real_diam_coefficient], [roughness]))
G.add_edge('junc_2', 'junc_4', obj=HE2_OilPipe([1478], [-16.5], [0.199 * real_diam_coefficient], [roughness])) #+
G.add_edge('PAD_5', 'junc_4', obj=HE2_OilPipe([600], [0.5], [0.09 * real_diam_coefficient], [roughness])) # +
G.add_edge('junc_4', 'junc_7', obj=HE2_OilPipe([1823], [3.7], [0.133 * real_diam_coefficient], [roughness])) #+
G.add_edge('PAD_155', 'junc_4', obj=HE2_OilPipe([703], [0.5], [0.141 * real_diam_coefficient], [roughness])) # +
#G.add_edge('PAD_87', 'junc_8', obj=HE2_OilPipe([2208], [1], [0.143], [1e-5])) # +
G.add_edge('junc_7', 'junc_8', obj=HE2_OilPipe([503], [0], [0.133  * real_diam_coefficient], [roughness])) #+
G.add_edge('PAD_35', 'junc_8', obj=HE2_OilPipe([499], [-3], [0.094 * real_diam_coefficient], [roughness])) # +
G.add_edge('junc_8', 'junc_9', obj=HE2_OilPipe([2648], [14.95], [0.199 * real_diam_coefficient], [roughness])) # +
G.add_edge('junc_9', 'DNS_2', obj=HE2_OilPipe([200], [2.7], [0.305 * real_diam_coefficient], [roughness]))  #+

solver = HE2_Solver(G)
solver.solve()
print(fDNS {G.nodes['DNS_2']['obj'].Q}
)b
"""

a = HE2_OilPipe([668], [-211.7], [0.073], [1e-5])
b = a.perform_calc_forward(22, 22, -5)
print(b)
print(a.segments[0].angle_dgr)

a = HE2_OilPipe([668], [211.7], [0.073], [1e-5])
b = a.perform_calc_backward(22, 22, 5)
print(b)