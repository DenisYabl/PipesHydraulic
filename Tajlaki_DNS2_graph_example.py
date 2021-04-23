import networkx as nx

from Fluids.HE2_Fluid import HE2_OilWater
from GraphNodes import HE2_Vertices as vrtxs
from GraphEdges.HE2_Pipe import HE2_OilPipe
from GraphEdges.HE2_Plast import HE2_Plast
from Solver.HE2_Solver import HE2_Solver
from GraphEdges.HE2_WellPump import HE2_WellPump, create_HE2_WellPump_instance_from_dataframe
from Fluids.oil_params import oil_params
import json
from Tools.HE2_Logger import check_for_nan, getLogger
logger = getLogger(__name__)

oil_params = oil_params(Q_m3_day=500, sat_P_bar=67, plastT_C=84, gasFactor=36, oildensity_kg_m3=826,
                        waterdensity_kg_m3=1015, gasdensity_kg_m3=1, oilviscosity_Pa_s=35e-3, volumewater_percent=50, volumeoilcoeff=1.017)
fluid = HE2_OilWater(oil_params)


def build_well_graph(G, wellhead_node, pad_name, well_name, pressures, plasts, pumps, pump_curves, fluid, inclination, roughness, d_keff):
    plast = f'PAD_{pad_name}_well_{well_name}'
    zaboi = f'Zaboi_{well_name}'
    intake = f'Pump_intake_{well_name}'
    outlet = f'Pump_outlet_{well_name}'

    plast_pressure = pressures[pad_name][well_name]
    inlets = {plast:vrtxs.HE2_Source_Vertex('P', plast_pressure, fluid, 20)}

    juncs = {zaboi: vrtxs.HE2_ABC_GraphVertex(),
             intake: vrtxs.HE2_ABC_GraphVertex(),
             outlet: vrtxs.HE2_ABC_GraphVertex(),
             }

    for k, v in {**inlets, **juncs}.items():
        G.add_node(k, obj=v)

    productivity = plasts[pad_name][well_name]
    G.add_edge(plast, zaboi, obj=HE2_Plast(productivity=productivity, fluid=fluid))

    pump = pumps[pad_name][well_name]
    G.add_edge(intake, outlet,obj=create_HE2_WellPump_instance_from_dataframe(full_HPX=pump_curves, model=pump[0], fluid=fluid, frequency=pump[1]))

    nkt_dx, nkt_dy, nkt_inner_d, colon_dx, colon_dy, colon_inner_d = inclination[well_name]
    G.add_edge(outlet, wellhead_node, obj=HE2_OilPipe([nkt_dx], [nkt_dy], [nkt_inner_d * d_keff], [roughness]))
    G.add_edge(zaboi, intake, obj=HE2_OilPipe([colon_dx], [colon_dy], [colon_inner_d * d_keff], [5 * roughness]))

    return G, inlets


def build_pad_graph(G, pad_node, pad_name, well_names, pressures, plasts, pumps, pump_curves, fluid, inclination, roughness, d_keff):
    inlets = {}
    for well_name in well_names:
        wellhead = f'Wellhead_{well_name}'
        G.add_node(wellhead, obj=vrtxs.HE2_ABC_GraphVertex())
        G.add_edge(wellhead, pad_node, obj=HE2_OilPipe([100], [0], [0.125 * d_keff], [roughness]))
        G, well_inlets = build_well_graph(G, wellhead, pad_name, well_name, pressures, plasts, pumps, pump_curves, fluid, inclination, roughness, d_keff)
        inlets.update(**well_inlets)
    return G, inlets


def build_DNS2_graph_by_pad_set(use_pads, pressures, plasts, pumps=None, pump_curves=None, fluid=None,
                                roughness=0.00001, real_diam_coefficient=1, daily_debit=0, DNS_pressure=4.8):
    json_str = json.dumps(pumps)
    logger.info(f'Pumps: {json_str}')

    G = nx.DiGraph()  # Di = directed
    outlets = dict(DNS_2=vrtxs.HE2_Boundary_Vertex('P', DNS_pressure))

    #Трубопроводная система
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
