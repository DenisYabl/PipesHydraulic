import networkx as nx

from Fluids.HE2_Fluid import HE2_OilWater
from Fluids.oil_params import oil_params
from GraphNodes import HE2_Vertices as vrtxs
from GraphEdges.HE2_Pipe import HE2_OilPipe
from GraphEdges.HE2_Plast import HE2_Plast
from GraphEdges.HE2_WellPump import HE2_WellPump
from Solver.HE2_Solver import HE2_Solver
import pandas as pd

calc_params = oil_params(Q_m3_day=500, sat_P_bar=67, plastT_C=84, gasFactor=36, oildensity_kg_m3=826,
                        waterdensity_kg_m3=1015, gasdensity_kg_m3=1, oilviscosity_Pa_s=35e-3, volumewater_percent=50, volumeoilcoeff=1.017)
fluid = HE2_OilWater(calc_params)


inlets = dict(PLAST_1816=vrtxs.HE2_Source_Vertex('P', 270.5, fluid, 20))

outlets = dict(wellhead = vrtxs.HE2_Boundary_Vertex('Q', 3.5))
#outlets = dict(wellhead = vrtxs.HE2_Boundary_Vertex('P', 11))


pump_curves = pd.read_csv("/home/denis/PycharmProjects/Pipeline/PipesHydraulic/CommonData/PumpChart.csv")

juncs = dict(Pump_intake=vrtxs.HE2_ABC_GraphVertex(),
            Pump_outlet=vrtxs.HE2_ABC_GraphVertex(),
             ZABOI_ZONE = vrtxs.HE2_ABC_GraphVertex())

#outlets = dict(wellhead = vrtxs.HE2_Boundary_Vertex('P', 8))

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
           obj=HE2_OilPipe([61.69], [2508.31], [0.073 * real_diam_coefficient], [roughness], [fluid]))
G.add_edge('ZABOI_ZONE', 'Pump_intake',
           obj=HE2_OilPipe([668], [211.7], [0.143 * real_diam_coefficient], [5 * roughness], [fluid]))

#Насос
G.add_edge('Pump_intake', 'Pump_outlet',
           obj=HE2_WellPump(full_HPX=pump_curves, model=pump_3270[0], fluid=fluid, frequency=pump_3270[1]))



solver = HE2_Solver(G)
solver.solve()
print(f"""Plast {G.nodes['PLAST_1816']['obj'].result}
ZABOI_ZONE {G.nodes['ZABOI_ZONE']['obj'].result}
Pump_intake {G.nodes['Pump_intake']['obj'].result}
Pump_outlet {G.nodes['Pump_outlet']['obj'].result}
wellhead {G.nodes['wellhead']['obj'].result}
""")


