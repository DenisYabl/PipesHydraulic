import networkx as nx

from HE2_Fluid import HE2_OilWater
import HE2_Vertices as vrtxs
from HE2_Pipe import HE2_OilPipe
from HE2_Plast import HE2_Plast
from HE2_WellPump import HE2_WellPump
from HE2_Solver import HE2_Solver
import pandas as pd

oil_params = {
    "OilSaturationP": 66.7, #Давление насыщения нефти при стандартных условиях, исходные данные
    "PlastT": 84.0, #Пластовая температура, исходные данные
    "GasFactor": 34.8, #Газовый фактор нефти, исходные данные
    "SepOilWeight": 0.885, #Плотность нефти, исходные данные
    "GasDensity": 1.003, #Плотность попутного газа, исходные данные
    "SepOilDynamicViscosity": 43.03, #Динамическая вязкость нефти, исходные данные
    "wellopVolumeWater": 92.8, #Обводненность нефти, исходные данные
    "VolumeOilCoeff": 1.097, #Объемный коэффициент нефти, исходные данные
    "PlastWaterWeight": 1.025, #Плотность попутной воды, исходные данные
}
fluid = HE2_OilWater(oil_params)

inlets = dict(PLAST_1816=vrtxs.HE2_Source_Vertex('P', 270.5, fluid, 20))


pump_curves = pd.read_csv("../CommonData/PumpChart.csv")

juncs = dict(Pump_intake=vrtxs.HE2_ABC_GraphVertex(),
            Pump_outlet=vrtxs.HE2_ABC_GraphVertex(),
             ZABOI_ZONE = vrtxs.HE2_ABC_GraphVertex())

#outlets = dict(wellhead = vrtxs.HE2_Boundary_Vertex('P', 25))
outlets = dict(wellhead = vrtxs.HE2_Boundary_Vertex('Q',  380 * 890 / 86400))

G = nx.DiGraph()  # Di = directed
for k, v in {**inlets, **outlets, **juncs}.items():
    G.add_node(k, obj=v)

pump_3270 = ["ЭЦН5А-320-2400", 47.5]

roughness = 1e-5
real_diam_coefficient = 1
#Пласт
G.add_edge('PLAST_1816', 'ZABOI_ZONE', obj=HE2_Plast(productivity = 1.66, fluid=fluid))  #+
#Инклинометрия до устья

G.add_edge('Pump_outlet', 'wellhead',
           obj=HE2_OilPipe([61.69], [2508.31], [0.073 * real_diam_coefficient], [roughness]))
G.add_edge('ZABOI_ZONE', 'Pump_intake',
           obj=HE2_OilPipe([668], [211.7], [0.143 * real_diam_coefficient], [5 * roughness]))

#Насос
G.add_edge('Pump_intake', 'Pump_outlet',
           obj=HE2_WellPump(full_HPX=pump_curves, model=pump_3270[0], fluid=fluid, IntDiameter=0.12, frequency=pump_3270[1]))



solver = HE2_Solver(G)
solver.solve()
print(f"""Plast {G.nodes['PLAST_1816']['obj'].result}
ZABOI_ZONE {G.nodes['ZABOI_ZONE']['obj'].result}
Pump_intake {G.nodes['Pump_intake']['obj'].result}
Pump_outlet {G.nodes['Pump_outlet']['obj'].result}
wellhead {G.nodes['wellhead']['obj'].result}
""")


