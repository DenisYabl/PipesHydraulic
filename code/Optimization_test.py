import networkx as nx

from HE2_Fluid import HE2_OilWater
import HE2_Vertices as vrtxs
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


def model_DNS_3(daily_debit_55 = 0, pressure_88 = 0, daily_debit = 0, fluid = fluid, roughness = 3.5, real_diam_coefficient = 0.85 ):
    inlets = dict(PAD_88=vrtxs.HE2_Source_Vertex('P', pressure_88, fluid, 20),
                  PAD_55=vrtxs.HE2_Source_Vertex('Q', daily_debit_55 * fluid.SepOilDensity * 1000 / 86400, fluid = fluid, T = 20))
                  #PAD_55=vrtxs.HE2_Source_Vertex('P', pressure_55, fluid, 20))

    juncs = dict(intake_pad_40=vrtxs.HE2_ABC_GraphVertex(),
                 intake_pad_45=vrtxs.HE2_ABC_GraphVertex(),
                 intake_pad_43=vrtxs.HE2_ABC_GraphVertex(),
                 intake_pad_44=vrtxs.HE2_ABC_GraphVertex())

    outlets = dict(DNS_3=vrtxs.HE2_Boundary_Vertex('Q', daily_debit * fluid.SepOilDensity * 1000 / 86400))

    G = nx.DiGraph()  # Di = directed
    for k, v in {**inlets, **outlets, **juncs}.items():
        G.add_node(k, obj=v)


    G.add_edge('PAD_55', 'intake_pad_40', obj=HE2_OilPipe([4871], [2.4], [0.143 * real_diam_coefficient], [roughness]))
    G.add_edge('intake_pad_40', 'intake_pad_45', obj=HE2_OilPipe([3562], [5], [0.203 * real_diam_coefficient], [roughness]))
    G.add_edge('PAD_88', 'intake_pad_45', obj=HE2_OilPipe([5234], [0], [0.143 * real_diam_coefficient], [roughness]))
    G.add_edge('intake_pad_45', 'intake_pad_43',obj=HE2_OilPipe([3042], [-6.5], [0.203 * real_diam_coefficient], [roughness]))
    G.add_edge('intake_pad_43', 'intake_pad_44',obj=HE2_OilPipe([2652], [2.1], [0.203 * real_diam_coefficient], [roughness]))
    G.add_edge('intake_pad_44', 'DNS_3',obj=HE2_OilPipe([15], [-1.5], [0.203 * real_diam_coefficient], [roughness]))
    solver = HE2_Solver(G)
    solver.solve()
    return G

def model_DNS_2(pressures:dict = {}, daily_debit = 0, fluid = None, roughness = 3.5, real_diam_coefficient = 0.85 ):
    # Давления в источниках
    pressure_5 = pressures["PAD_5"]
    pressure_33 = pressures["PAD_33"]
    pressure_34 = pressures["PAD_34"]
    pressure_39 = pressures["PAD_39"]
    pressure_49 = pressures["PAD_49"]
    pressure_57 = pressures["PAD_57"]

    #Источники системы
    inlets = dict(PAD_5=vrtxs.HE2_Source_Vertex('P', pressure_5, fluid, 20), #+
                  PAD_33=vrtxs.HE2_Source_Vertex('P', pressure_33, fluid, 20), #+
                  PAD_34=vrtxs.HE2_Source_Vertex('P', pressure_34, fluid, 20), #+
                  PAD_39=vrtxs.HE2_Source_Vertex('P', pressure_39, fluid, 20), #+
                  PAD_49=vrtxs.HE2_Source_Vertex('P', pressure_49, fluid, 20), #+
                  PAD_57=vrtxs.HE2_Source_Vertex('P', pressure_57, fluid, 20), #+
                  )
                  #PAD_55=vrtxs.HE2_Source_Vertex('P', pressure_5, fluid, 20))
    # Узловые точки
    juncs = dict(intake_pad_5=vrtxs.HE2_ABC_GraphVertex(),
                 intake_pad_9_36=vrtxs.HE2_ABC_GraphVertex(),
                 intake_pad_35=vrtxs.HE2_ABC_GraphVertex(),
                 intake_pad_33=vrtxs.HE2_ABC_GraphVertex(),
                 intake_node_33=vrtxs.HE2_ABC_GraphVertex(),
                 intake_node_49=vrtxs.HE2_ABC_GraphVertex(),
                 intake_pad_59=vrtxs.HE2_ABC_GraphVertex(),
                 node_15=vrtxs.HE2_ABC_GraphVertex(),
                 intake_pad_46=vrtxs.HE2_ABC_GraphVertex(),
                 intake_pad_47=vrtxs.HE2_ABC_GraphVertex(),
                 intake_pad_57=vrtxs.HE2_ABC_GraphVertex(),
                 intake_pad_41=vrtxs.HE2_ABC_GraphVertex(),
                 intake_pad_49=vrtxs.HE2_ABC_GraphVertex(),
                 intake_pad_53=vrtxs.HE2_ABC_GraphVertex(),
                 intake_pad_134=vrtxs.HE2_ABC_GraphVertex(),
                 intake_pad_59_46=vrtxs.HE2_ABC_GraphVertex(),
                 UDR_1=vrtxs.HE2_ABC_GraphVertex(),
                 UDR_2=vrtxs.HE2_ABC_GraphVertex(),
                 ZKL_98=vrtxs.HE2_ABC_GraphVertex())

    outlets = dict(DNS_2=vrtxs.HE2_Boundary_Vertex('Q', daily_debit * fluid.SepOilDensity * 1000 / 86400))

    G = nx.DiGraph()  # Di = directed
    for k, v in {**inlets, **outlets, **juncs}.items():
        G.add_node(k, obj=v)

    #Описание трубопроводных соединений между узлами
    G.add_edge('PAD_5', 'intake_pad_5', obj=HE2_OilPipe([10], [-0.2], [0.143 * real_diam_coefficient], [roughness]))
    G.add_edge('intake_pad_5', 'intake_pad_9_36', obj=HE2_OilPipe([785], [-0.6], [0.199 * real_diam_coefficient], [roughness]))
    G.add_edge('intake_pad_9_36', 'intake_pad_35', obj=HE2_OilPipe([2326], [4.6], [0.199 * real_diam_coefficient], [roughness]))
    G.add_edge('intake_pad_35', 'UDR_1',obj=HE2_OilPipe([2648], [13.6], [0.199 * real_diam_coefficient], [roughness]))
    G.add_edge('UDR_1', 'DNS_2',obj=HE2_OilPipe([200], [0.2], [0.305 * real_diam_coefficient], [roughness]))
    G.add_edge('PAD_33', 'intake_pad_33',obj=HE2_OilPipe([440], [-1.9], [0.139 * real_diam_coefficient], [roughness]))
    G.add_edge('intake_pad_33', 'UDR_2', obj=HE2_OilPipe([4394], [-4.6], [0.199 * real_diam_coefficient], [roughness]))
    G.add_edge('UDR_2', 'UDR_1', obj=HE2_OilPipe([1087], [3.4], [0.253 * real_diam_coefficient], [roughness]))
    G.add_edge('PAD_34', 'intake_pad_134', obj=HE2_OilPipe([818], [0.6], [0.143 * real_diam_coefficient], [roughness]))
    G.add_edge('intake_pad_134', 'UDR_2', obj=HE2_OilPipe([3344], [10.2], [0.203 * real_diam_coefficient], [roughness]))
    G.add_edge('PAD_39', 'intake_pad_59', obj=HE2_OilPipe([7568], [1.6], [0.139 * real_diam_coefficient], [roughness]))
    G.add_edge('intake_pad_59', 'intake_pad_59_46', obj=HE2_OilPipe([4601], [4.4], [0.143 * real_diam_coefficient], [roughness]))
    G.add_edge('intake_pad_59_46', 'intake_pad_47',obj=HE2_OilPipe([2625], [18.4], [0.139 * real_diam_coefficient], [roughness]))
    G.add_edge('intake_pad_47', 'intake_pad_53',obj=HE2_OilPipe([2250], [-1.3], [0.203 * real_diam_coefficient], [roughness]))
    G.add_edge('intake_pad_53', 'intake_pad_41',obj=HE2_OilPipe([1391], [9], [0.203 * real_diam_coefficient], [roughness]))
    G.add_edge('intake_pad_41', 'ZKL_98',obj=HE2_OilPipe([9290], [-13.5], [0.257 * real_diam_coefficient], [roughness]))
    G.add_edge('ZKL_98', 'intake_pad_49',obj=HE2_OilPipe([600], [3], [0.257 * real_diam_coefficient], [roughness]))
    G.add_edge('intake_pad_49', 'intake_node_33', obj=HE2_OilPipe([4001], [0], [0.410 * real_diam_coefficient], [roughness]))
    G.add_edge('intake_node_33', 'intake_pad_46', obj=HE2_OilPipe([2920], [2.1], [0.410 * real_diam_coefficient], [roughness]))
    G.add_edge('intake_pad_46', 'node_15',obj=HE2_OilPipe([1181], [-5.9], [0.410 * real_diam_coefficient], [roughness]))
    G.add_edge('node_15', 'UDR_2',obj=HE2_OilPipe([92], [0.2], [0.199 * real_diam_coefficient], [roughness]))
    G.add_edge('PAD_49', 'intake_node_49', obj=HE2_OilPipe([492], [-0.8], [0.139 * real_diam_coefficient], [roughness]))
    G.add_edge('intake_node_49', 'intake_pad_57', obj=HE2_OilPipe([3737], [-5.7], [0.143 * real_diam_coefficient], [roughness]))
    G.add_edge('intake_pad_57', 'intake_pad_49',obj=HE2_OilPipe([3852], [-4.3], [0.309 * real_diam_coefficient], [roughness]))
    G.add_edge('PAD_57', 'intake_pad_57',obj=HE2_OilPipe([370], [-1], [0.203 * real_diam_coefficient], [roughness]))
    solver = HE2_Solver(G)
    solver.solve()
    return G

