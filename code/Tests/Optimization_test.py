import networkx as nx

from Fluids.HE2_Fluid import HE2_BlackOil
from GraphNodes import HE2_Vertices as vrtxs
from GraphEdges.HE2_Pipe import HE2_OilPipe
from GraphEdges.HE2_Plast import HE2_Plast
from Solver.HE2_Solver import HE2_Solver
from GraphEdges.HE2_WellPump import HE2_WellPump, create_HE2_WellPump_instance_from_dataframe
from Tools.HE2_ABC import oil_params
import json
import colorama
from colorama import Fore, Back, Style
from Tools.HE2_Logger import check_for_nan, getLogger
import numpy as np
logger = getLogger(__name__)

class NpEncoder(json.JSONEncoder):
    def default(self, obj):
        if isinstance(obj, np.integer):
            return int(obj)
        elif isinstance(obj, np.floating):
            return float(obj)
        elif isinstance(obj, np.ndarray):
            return obj.tolist()
        else:
            return super(NpEncoder, self).default(obj)

oil_params = oil_params(sat_P_bar=67, plastT_C=84, gasFactor=36, oildensity_kg_m3=826, waterdensity_kg_m3=1015,
                        gasdensity_kg_m3=1, oilviscosity_Pa_s=35e-3, volumewater_percent=50, volumeoilcoeff=1.017)

fluid = HE2_BlackOil(oil_params)


def gimme_DNS2_inlets_outlets_Q():
    inlets_outlets_Q = dict(
                  PAD_5_well_1523=0.66,  #Куст 5
                  PAD_5_well_146=0.889,
                  PAD_5_well_142=0.442,
                  PAD_5_well_1562=0.645,

                  PAD_33_well_1385=0.53,  #Куст 33
                  PAD_33_well_736=0.494,
                  PAD_33_well_739=0.912,
                  PAD_33_well_1383=1.098,
                  PAD_33_well_738=0.8937,
                  PAD_33_well_725=1.0907,

                  PAD_34_well_731=0.912,  #Куст 34
                  PAD_34_well_196=1.459,
                  PAD_34_well_734=1.158,
                  PAD_34_well_198=0.991,
                  PAD_34_well_199=1.064,
                  PAD_34_well_197=0.6709,
                  PAD_34_well_195=0.698,
                  PAD_34_well_191=0.601,
                  PAD_34_well_729=2.99,
                  PAD_34_well_730=0.492,
                  PAD_34_well_192=1.056,
                  PAD_34_well_148=0.631,

                  PAD_39_well_3552=1.977,  # Куст 34
                  PAD_39_well_617=0.9049,
                  PAD_39_well_567=0.1068,
                  PAD_39_well_614=2.083,
                  PAD_39_well_619=0.8549,
                  PAD_39_well_609=0.61696,

                  PAD_49_well_1816=3.9381,
                  PAD_49_well_2630=0.14,
                  PAD_49_well_1815=1.414,
                  PAD_49_well_676=1.92,
                  PAD_49_well_3270=0.6378,
                  PAD_49_well_3266=1.0612,
                  PAD_49_well_1814=2.6779,
                  PAD_49_well_1817=1.8012,
                  PAD_49_well_4532=0.67042,
                  PAD_49_well_2631=1.7264,
                  PAD_49_well_677=0.98671,

                  PAD_57_well_3113=1.3827,
                  PAD_57_well_3118=0.15678,
                  PAD_57_well_3112=1.139,
                  PAD_57_well_4235=0.4738,
                  PAD_57_well_3117=2.940,
                  PAD_57_well_1493=1.833,
                  PAD_57_well_1574=1.42545,
                  PAD_57_well_1579=0.54595,
                  PAD_57_well_3116=2.62972,
                  )
    total_q = sum(inlets_outlets_Q.values())
    inlets_outlets_Q.update(DNS_2=-total_q)
    return inlets_outlets_Q

def build_DNS2_graph(pressures: dict = {}, plasts: dict = {}, pumps=None, pump_curves=None, fluid=None,
                     roughness=0.00001, real_diam_coefficient=1, DNS_daily_debit=0, DNS_pressure=4.8):

    json_str = json.dumps(pumps, cls=NpEncoder)
    logger.info(f'Pumps: {json_str}')

    # Давления в источниках
    pressure_1523 = pressures["PAD_5"]["WELL_1523"]
    pressure_146 = pressures["PAD_5"]["WELL_146"]
    pressure_142 = pressures["PAD_5"]["WELL_142"]
    pressure_1562 = pressures["PAD_5"]["WELL_1562"]

    pressure_1385 = pressures["PAD_33"]["WELL_1385"]
    pressure_736 = pressures["PAD_33"]["WELL_736"]
    pressure_739 = pressures["PAD_33"]["WELL_739"]
    pressure_1383 = pressures["PAD_33"]["WELL_1383"]
    pressure_738 = pressures["PAD_33"]["WELL_738"]
    pressure_725 = pressures["PAD_33"]["WELL_725"]

    pressure_731 = pressures["PAD_34"]["WELL_731"]
    pressure_196 = pressures["PAD_34"]["WELL_196"]
    pressure_734 = pressures["PAD_34"]["WELL_734"]
    pressure_198 = pressures["PAD_34"]["WELL_198"]
    pressure_199 = pressures["PAD_34"]["WELL_199"]
    pressure_197 = pressures["PAD_34"]["WELL_197"]
    pressure_195 = pressures["PAD_34"]["WELL_195"]
    pressure_191 = pressures["PAD_34"]["WELL_191"]
    pressure_729 = pressures["PAD_34"]["WELL_729"]
    pressure_730 = pressures["PAD_34"]["WELL_730"]
    pressure_192 = pressures["PAD_34"]["WELL_192"]
    pressure_148 = pressures["PAD_34"]["WELL_148"]

    pressure_3552 = pressures["PAD_39"]["WELL_3552"]
    pressure_617 = pressures["PAD_39"]["WELL_617"]
    pressure_567 = pressures["PAD_39"]["WELL_567"]
    pressure_614 = pressures["PAD_39"]["WELL_614"]
    pressure_619 = pressures["PAD_39"]["WELL_619"]
    pressure_609 = pressures["PAD_39"]["WELL_609"]
    
    
    pressure_1816 = pressures["PAD_49"]["WELL_1816"]
    pressure_2630 = pressures["PAD_49"]["WELL_2630"]
    pressure_1815 = pressures["PAD_49"]["WELL_1815"]
    pressure_676 = pressures["PAD_49"]["WELL_676"]
    pressure_3270 = pressures["PAD_49"]["WELL_3270"]
    pressure_3266 = pressures["PAD_49"]["WELL_3266"]
    pressure_1814 = pressures["PAD_49"]["WELL_1814"]
    pressure_1817 = pressures["PAD_49"]["WELL_1817"]
    pressure_4532 = pressures["PAD_49"]["WELL_4532"]
    pressure_2631 = pressures["PAD_49"]["WELL_2631"]
    pressure_677 = pressures["PAD_49"]["WELL_677"]

    pressure_3113 = pressures["PAD_57"]["WELL_3113"]
    pressure_3118 = pressures["PAD_57"]["WELL_3118"]
    pressure_3112 = pressures["PAD_57"]["WELL_3112"]
    pressure_4235 = pressures["PAD_57"]["WELL_4235"]
    pressure_3117 = pressures["PAD_57"]["WELL_3117"]
    pressure_1493 = pressures["PAD_57"]["WELL_1493"]
    pressure_1574 = pressures["PAD_57"]["WELL_1574"]
    pressure_1579 = pressures["PAD_57"]["WELL_1579"]
    pressure_3116 = pressures["PAD_57"]["WELL_3116"]
    
    #Источники системы
    inlets = dict(PAD_5_well_1523=vrtxs.HE2_Source_Vertex('P', pressure_1523, fluid, 20),  #Куст 5
                  PAD_5_well_146=vrtxs.HE2_Source_Vertex('P', pressure_146, fluid, 20),
                  PAD_5_well_142=vrtxs.HE2_Source_Vertex('P', pressure_142, fluid, 20),
                  PAD_5_well_1562=vrtxs.HE2_Source_Vertex('P', pressure_1562, fluid, 20),

                  PAD_33_well_1385=vrtxs.HE2_Source_Vertex('P', pressure_1385, fluid, 20),  #Куст 33
                  PAD_33_well_736=vrtxs.HE2_Source_Vertex('P', pressure_736, fluid, 20),
                  PAD_33_well_739=vrtxs.HE2_Source_Vertex('P', pressure_739, fluid, 20),
                  PAD_33_well_1383=vrtxs.HE2_Source_Vertex('P', pressure_1383, fluid, 20),
                  PAD_33_well_738=vrtxs.HE2_Source_Vertex('P', pressure_738, fluid, 20),
                  PAD_33_well_725=vrtxs.HE2_Source_Vertex('P', pressure_725, fluid, 20),


                  PAD_34_well_731=vrtxs.HE2_Source_Vertex('P', pressure_731, fluid, 20),  #Куст 34
                  PAD_34_well_196=vrtxs.HE2_Source_Vertex('P', pressure_196, fluid, 20),
                  PAD_34_well_734=vrtxs.HE2_Source_Vertex('P', pressure_734, fluid, 20),
                  PAD_34_well_198=vrtxs.HE2_Source_Vertex('P', pressure_198, fluid, 20),
                  PAD_34_well_199=vrtxs.HE2_Source_Vertex('P', pressure_199, fluid, 20),
                  PAD_34_well_197=vrtxs.HE2_Source_Vertex('P', pressure_197, fluid, 20),
                  PAD_34_well_195=vrtxs.HE2_Source_Vertex('P', pressure_195, fluid, 20),
                  PAD_34_well_191=vrtxs.HE2_Source_Vertex('P', pressure_191, fluid, 20),
                  PAD_34_well_729=vrtxs.HE2_Source_Vertex('P', pressure_729, fluid, 20),
                  PAD_34_well_730=vrtxs.HE2_Source_Vertex('P', pressure_730, fluid, 20),
                  PAD_34_well_192=vrtxs.HE2_Source_Vertex('P', pressure_192, fluid, 20),
                  PAD_34_well_148=vrtxs.HE2_Source_Vertex('P', pressure_148, fluid, 20),

                  PAD_39_well_3552=vrtxs.HE2_Source_Vertex('P', pressure_3552, fluid, 20),  # Куст 34
                  PAD_39_well_617=vrtxs.HE2_Source_Vertex('P', pressure_617, fluid, 20),  
                  PAD_39_well_567=vrtxs.HE2_Source_Vertex('P', pressure_567, fluid, 20),  
                  PAD_39_well_614=vrtxs.HE2_Source_Vertex('P', pressure_614, fluid, 20),  
                  PAD_39_well_619=vrtxs.HE2_Source_Vertex('P', pressure_619, fluid, 20),  
                  PAD_39_well_609=vrtxs.HE2_Source_Vertex('P', pressure_609, fluid, 20),

                  PAD_49_well_1816=vrtxs.HE2_Source_Vertex('P', pressure_1816, fluid, 20),
                  PAD_49_well_2630=vrtxs.HE2_Source_Vertex('P', pressure_2630, fluid, 20),
                  PAD_49_well_1815=vrtxs.HE2_Source_Vertex('P', pressure_1815, fluid, 20),
                  PAD_49_well_676=vrtxs.HE2_Source_Vertex('P', pressure_676, fluid, 20),
                  PAD_49_well_3270=vrtxs.HE2_Source_Vertex('P', pressure_3270, fluid, 20),
                  PAD_49_well_3266=vrtxs.HE2_Source_Vertex('P', pressure_3266, fluid, 20),
                  PAD_49_well_1814=vrtxs.HE2_Source_Vertex('P', pressure_1814, fluid, 20),
                  PAD_49_well_1817=vrtxs.HE2_Source_Vertex('P', pressure_1817, fluid, 20),
                  PAD_49_well_4532=vrtxs.HE2_Source_Vertex('P', pressure_4532, fluid, 20),
                  PAD_49_well_2631=vrtxs.HE2_Source_Vertex('P', pressure_2631, fluid, 20),
                  PAD_49_well_677=vrtxs.HE2_Source_Vertex('P', pressure_677, fluid, 20),

                  PAD_57_well_3113=vrtxs.HE2_Source_Vertex('P', pressure_3113, fluid, 20),
                  PAD_57_well_3118=vrtxs.HE2_Source_Vertex('P', pressure_3118, fluid, 20),
                  PAD_57_well_3112=vrtxs.HE2_Source_Vertex('P', pressure_3112, fluid, 20),
                  PAD_57_well_4235=vrtxs.HE2_Source_Vertex('P', pressure_4235, fluid, 20),
                  PAD_57_well_3117=vrtxs.HE2_Source_Vertex('P', pressure_3117, fluid, 20),
                  PAD_57_well_1493=vrtxs.HE2_Source_Vertex('P', pressure_1493, fluid, 20),
                  PAD_57_well_1574=vrtxs.HE2_Source_Vertex('P', pressure_1574, fluid, 20),
                  PAD_57_well_1579=vrtxs.HE2_Source_Vertex('P', pressure_1579, fluid, 20),
                  PAD_57_well_3116=vrtxs.HE2_Source_Vertex('P', pressure_3116, fluid, 20),
                  )



    # Узловые точки
    juncs = dict(Zaboi_1523 = vrtxs.HE2_ABC_GraphVertex(),
                 Pump_intake_1523=vrtxs.HE2_ABC_GraphVertex(),
                 Pump_outlet_1523=vrtxs.HE2_ABC_GraphVertex(),
                 Wellhead_1523=vrtxs.HE2_ABC_GraphVertex(),

                 Zaboi_146=vrtxs.HE2_ABC_GraphVertex(),
                 Pump_intake_146=vrtxs.HE2_ABC_GraphVertex(),
                 Pump_outlet_146=vrtxs.HE2_ABC_GraphVertex(),
                 Wellhead_146=vrtxs.HE2_ABC_GraphVertex(),

                 Zaboi_142=vrtxs.HE2_ABC_GraphVertex(),
                 Pump_intake_142=vrtxs.HE2_ABC_GraphVertex(),
                 Pump_outlet_142=vrtxs.HE2_ABC_GraphVertex(),
                 Wellhead_142=vrtxs.HE2_ABC_GraphVertex(),

                 Zaboi_1562=vrtxs.HE2_ABC_GraphVertex(),
                 Pump_intake_1562=vrtxs.HE2_ABC_GraphVertex(),
                 Pump_outlet_1562=vrtxs.HE2_ABC_GraphVertex(),
                 Wellhead_1562=vrtxs.HE2_ABC_GraphVertex(),

                 Zaboi_1385=vrtxs.HE2_ABC_GraphVertex(),
                 Pump_intake_1385=vrtxs.HE2_ABC_GraphVertex(),
                 Pump_outlet_1385=vrtxs.HE2_ABC_GraphVertex(),
                 Wellhead_1385=vrtxs.HE2_ABC_GraphVertex(),

                 Zaboi_736=vrtxs.HE2_ABC_GraphVertex(),
                 Pump_intake_736=vrtxs.HE2_ABC_GraphVertex(),
                 Pump_outlet_736=vrtxs.HE2_ABC_GraphVertex(),
                 Wellhead_736=vrtxs.HE2_ABC_GraphVertex(),

                 Zaboi_739=vrtxs.HE2_ABC_GraphVertex(),
                 Pump_intake_739=vrtxs.HE2_ABC_GraphVertex(),
                 Pump_outlet_739=vrtxs.HE2_ABC_GraphVertex(),
                 Wellhead_739=vrtxs.HE2_ABC_GraphVertex(),

                 Zaboi_1383=vrtxs.HE2_ABC_GraphVertex(),
                 Pump_intake_1383=vrtxs.HE2_ABC_GraphVertex(),
                 Pump_outlet_1383=vrtxs.HE2_ABC_GraphVertex(),
                 Wellhead_1383=vrtxs.HE2_ABC_GraphVertex(),

                 Zaboi_738=vrtxs.HE2_ABC_GraphVertex(),
                 Pump_intake_738=vrtxs.HE2_ABC_GraphVertex(),
                 Pump_outlet_738=vrtxs.HE2_ABC_GraphVertex(),
                 Wellhead_738=vrtxs.HE2_ABC_GraphVertex(),

                 Zaboi_725=vrtxs.HE2_ABC_GraphVertex(),
                 Pump_intake_725=vrtxs.HE2_ABC_GraphVertex(),
                 Pump_outlet_725=vrtxs.HE2_ABC_GraphVertex(),
                 Wellhead_725=vrtxs.HE2_ABC_GraphVertex(),

                 Zaboi_731=vrtxs.HE2_ABC_GraphVertex(),
                 Pump_intake_731=vrtxs.HE2_ABC_GraphVertex(),
                 Pump_outlet_731=vrtxs.HE2_ABC_GraphVertex(),
                 Wellhead_731=vrtxs.HE2_ABC_GraphVertex(),

                 Zaboi_196=vrtxs.HE2_ABC_GraphVertex(),
                 Pump_intake_196=vrtxs.HE2_ABC_GraphVertex(),
                 Pump_outlet_196=vrtxs.HE2_ABC_GraphVertex(),
                 Wellhead_196=vrtxs.HE2_ABC_GraphVertex(),

                 Zaboi_734=vrtxs.HE2_ABC_GraphVertex(),
                 Pump_intake_734=vrtxs.HE2_ABC_GraphVertex(),
                 Pump_outlet_734=vrtxs.HE2_ABC_GraphVertex(),
                 Wellhead_734=vrtxs.HE2_ABC_GraphVertex(),

                 Zaboi_198=vrtxs.HE2_ABC_GraphVertex(),
                 Pump_intake_198=vrtxs.HE2_ABC_GraphVertex(),
                 Pump_outlet_198=vrtxs.HE2_ABC_GraphVertex(),
                 Wellhead_198=vrtxs.HE2_ABC_GraphVertex(),

                 Zaboi_199=vrtxs.HE2_ABC_GraphVertex(),
                 Pump_intake_199=vrtxs.HE2_ABC_GraphVertex(),
                 Pump_outlet_199=vrtxs.HE2_ABC_GraphVertex(),
                 Wellhead_199=vrtxs.HE2_ABC_GraphVertex(),

                 Zaboi_197=vrtxs.HE2_ABC_GraphVertex(),
                 Pump_intake_197=vrtxs.HE2_ABC_GraphVertex(),
                 Pump_outlet_197=vrtxs.HE2_ABC_GraphVertex(),
                 Wellhead_197=vrtxs.HE2_ABC_GraphVertex(),

                 Zaboi_195=vrtxs.HE2_ABC_GraphVertex(),
                 Pump_intake_195=vrtxs.HE2_ABC_GraphVertex(),
                 Pump_outlet_195=vrtxs.HE2_ABC_GraphVertex(),
                 Wellhead_195=vrtxs.HE2_ABC_GraphVertex(),

                 Zaboi_191=vrtxs.HE2_ABC_GraphVertex(),
                 Pump_intake_191=vrtxs.HE2_ABC_GraphVertex(),
                 Pump_outlet_191=vrtxs.HE2_ABC_GraphVertex(),
                 Wellhead_191=vrtxs.HE2_ABC_GraphVertex(),

                 Zaboi_729=vrtxs.HE2_ABC_GraphVertex(),
                 Pump_intake_729=vrtxs.HE2_ABC_GraphVertex(),
                 Pump_outlet_729=vrtxs.HE2_ABC_GraphVertex(),
                 Wellhead_729=vrtxs.HE2_ABC_GraphVertex(),

                 Zaboi_730=vrtxs.HE2_ABC_GraphVertex(),
                 Pump_intake_730=vrtxs.HE2_ABC_GraphVertex(),
                 Pump_outlet_730=vrtxs.HE2_ABC_GraphVertex(),
                 Wellhead_730=vrtxs.HE2_ABC_GraphVertex(),

                 Zaboi_192=vrtxs.HE2_ABC_GraphVertex(),
                 Pump_intake_192=vrtxs.HE2_ABC_GraphVertex(),
                 Pump_outlet_192=vrtxs.HE2_ABC_GraphVertex(),
                 Wellhead_192=vrtxs.HE2_ABC_GraphVertex(),

                 Zaboi_148=vrtxs.HE2_ABC_GraphVertex(),
                 Pump_intake_148=vrtxs.HE2_ABC_GraphVertex(),
                 Pump_outlet_148=vrtxs.HE2_ABC_GraphVertex(),
                 Wellhead_148=vrtxs.HE2_ABC_GraphVertex(),


                 Zaboi_3552=vrtxs.HE2_ABC_GraphVertex(),
                 Pump_intake_3552=vrtxs.HE2_ABC_GraphVertex(),
                 Pump_outlet_3552=vrtxs.HE2_ABC_GraphVertex(),
                 Wellhead_3552=vrtxs.HE2_ABC_GraphVertex(),

                 Zaboi_617=vrtxs.HE2_ABC_GraphVertex(),
                 Pump_intake_617=vrtxs.HE2_ABC_GraphVertex(),
                 Pump_outlet_617=vrtxs.HE2_ABC_GraphVertex(),
                 Wellhead_617=vrtxs.HE2_ABC_GraphVertex(),

                 Zaboi_567=vrtxs.HE2_ABC_GraphVertex(),
                 Pump_intake_567=vrtxs.HE2_ABC_GraphVertex(),
                 Pump_outlet_567=vrtxs.HE2_ABC_GraphVertex(),
                 Wellhead_567=vrtxs.HE2_ABC_GraphVertex(),

                 Zaboi_614=vrtxs.HE2_ABC_GraphVertex(),
                 Pump_intake_614=vrtxs.HE2_ABC_GraphVertex(),
                 Pump_outlet_614=vrtxs.HE2_ABC_GraphVertex(),
                 Wellhead_614=vrtxs.HE2_ABC_GraphVertex(),

                 Zaboi_619=vrtxs.HE2_ABC_GraphVertex(),
                 Pump_intake_619=vrtxs.HE2_ABC_GraphVertex(),
                 Pump_outlet_619=vrtxs.HE2_ABC_GraphVertex(),
                 Wellhead_619=vrtxs.HE2_ABC_GraphVertex(),

                 Zaboi_609=vrtxs.HE2_ABC_GraphVertex(),
                 Pump_intake_609=vrtxs.HE2_ABC_GraphVertex(),
                 Pump_outlet_609=vrtxs.HE2_ABC_GraphVertex(),
                 Wellhead_609=vrtxs.HE2_ABC_GraphVertex(),

                 Zaboi_1816=vrtxs.HE2_ABC_GraphVertex(),
                 Pump_intake_1816=vrtxs.HE2_ABC_GraphVertex(),
                 Pump_outlet_1816=vrtxs.HE2_ABC_GraphVertex(),
                 Wellhead_1816=vrtxs.HE2_ABC_GraphVertex(),

                 Zaboi_2630=vrtxs.HE2_ABC_GraphVertex(),
                 Pump_intake_2630=vrtxs.HE2_ABC_GraphVertex(),
                 Pump_outlet_2630=vrtxs.HE2_ABC_GraphVertex(),
                 Wellhead_2630=vrtxs.HE2_ABC_GraphVertex(),

                 Zaboi_1815=vrtxs.HE2_ABC_GraphVertex(),
                 Pump_intake_1815=vrtxs.HE2_ABC_GraphVertex(),
                 Pump_outlet_1815=vrtxs.HE2_ABC_GraphVertex(),
                 Wellhead_1815=vrtxs.HE2_ABC_GraphVertex(),

                 Zaboi_676=vrtxs.HE2_ABC_GraphVertex(),
                 Pump_intake_676=vrtxs.HE2_ABC_GraphVertex(),
                 Pump_outlet_676=vrtxs.HE2_ABC_GraphVertex(),
                 Wellhead_676=vrtxs.HE2_ABC_GraphVertex(),

                 Zaboi_3270=vrtxs.HE2_ABC_GraphVertex(),
                 Pump_intake_3270=vrtxs.HE2_ABC_GraphVertex(),
                 Pump_outlet_3270=vrtxs.HE2_ABC_GraphVertex(),
                 Wellhead_3270=vrtxs.HE2_ABC_GraphVertex(),

                 Zaboi_3266=vrtxs.HE2_ABC_GraphVertex(),
                 Pump_intake_3266=vrtxs.HE2_ABC_GraphVertex(),
                 Pump_outlet_3266=vrtxs.HE2_ABC_GraphVertex(),
                 Wellhead_3266=vrtxs.HE2_ABC_GraphVertex(),

                 Zaboi_1814=vrtxs.HE2_ABC_GraphVertex(),
                 Pump_intake_1814=vrtxs.HE2_ABC_GraphVertex(),
                 Pump_outlet_1814=vrtxs.HE2_ABC_GraphVertex(),
                 Wellhead_1814=vrtxs.HE2_ABC_GraphVertex(),

                 Zaboi_1817=vrtxs.HE2_ABC_GraphVertex(),
                 Pump_intake_1817=vrtxs.HE2_ABC_GraphVertex(),
                 Pump_outlet_1817=vrtxs.HE2_ABC_GraphVertex(),
                 Wellhead_1817=vrtxs.HE2_ABC_GraphVertex(),

                 Zaboi_4532=vrtxs.HE2_ABC_GraphVertex(),
                 Pump_intake_4532=vrtxs.HE2_ABC_GraphVertex(),
                 Pump_outlet_4532=vrtxs.HE2_ABC_GraphVertex(),
                 Wellhead_4532=vrtxs.HE2_ABC_GraphVertex(),

                 Zaboi_2631=vrtxs.HE2_ABC_GraphVertex(),
                 Pump_intake_2631=vrtxs.HE2_ABC_GraphVertex(),
                 Pump_outlet_2631=vrtxs.HE2_ABC_GraphVertex(),
                 Wellhead_2631=vrtxs.HE2_ABC_GraphVertex(),

                 Zaboi_677=vrtxs.HE2_ABC_GraphVertex(),
                 Pump_intake_677=vrtxs.HE2_ABC_GraphVertex(),
                 Pump_outlet_677=vrtxs.HE2_ABC_GraphVertex(),
                 Wellhead_677=vrtxs.HE2_ABC_GraphVertex(),

                 Zaboi_3113=vrtxs.HE2_ABC_GraphVertex(),
                 Pump_intake_3113=vrtxs.HE2_ABC_GraphVertex(),
                 Pump_outlet_3113=vrtxs.HE2_ABC_GraphVertex(),
                 Wellhead_3113=vrtxs.HE2_ABC_GraphVertex(),

                 Zaboi_3118=vrtxs.HE2_ABC_GraphVertex(),
                 Pump_intake_3118=vrtxs.HE2_ABC_GraphVertex(),
                 Pump_outlet_3118=vrtxs.HE2_ABC_GraphVertex(),
                 Wellhead_3118=vrtxs.HE2_ABC_GraphVertex(),

                 Zaboi_3112=vrtxs.HE2_ABC_GraphVertex(),
                 Pump_intake_3112=vrtxs.HE2_ABC_GraphVertex(),
                 Pump_outlet_3112=vrtxs.HE2_ABC_GraphVertex(),
                 Wellhead_3112=vrtxs.HE2_ABC_GraphVertex(),

                 Zaboi_4235=vrtxs.HE2_ABC_GraphVertex(),
                 Pump_intake_4235=vrtxs.HE2_ABC_GraphVertex(),
                 Pump_outlet_4235=vrtxs.HE2_ABC_GraphVertex(),
                 Wellhead_4235=vrtxs.HE2_ABC_GraphVertex(),

                 Zaboi_3117=vrtxs.HE2_ABC_GraphVertex(),
                 Pump_intake_3117=vrtxs.HE2_ABC_GraphVertex(),
                 Pump_outlet_3117=vrtxs.HE2_ABC_GraphVertex(),
                 Wellhead_3117=vrtxs.HE2_ABC_GraphVertex(),

                 Zaboi_1493=vrtxs.HE2_ABC_GraphVertex(),
                 Pump_intake_1493=vrtxs.HE2_ABC_GraphVertex(),
                 Pump_outlet_1493=vrtxs.HE2_ABC_GraphVertex(),
                 Wellhead_1493=vrtxs.HE2_ABC_GraphVertex(),

                 Zaboi_1574=vrtxs.HE2_ABC_GraphVertex(),
                 Pump_intake_1574=vrtxs.HE2_ABC_GraphVertex(),
                 Pump_outlet_1574=vrtxs.HE2_ABC_GraphVertex(),
                 Wellhead_1574=vrtxs.HE2_ABC_GraphVertex(),

                 Zaboi_1579=vrtxs.HE2_ABC_GraphVertex(),
                 Pump_intake_1579=vrtxs.HE2_ABC_GraphVertex(),
                 Pump_outlet_1579=vrtxs.HE2_ABC_GraphVertex(),
                 Wellhead_1579=vrtxs.HE2_ABC_GraphVertex(),

                 Zaboi_3116=vrtxs.HE2_ABC_GraphVertex(),
                 Pump_intake_3116=vrtxs.HE2_ABC_GraphVertex(),
                 Pump_outlet_3116=vrtxs.HE2_ABC_GraphVertex(),
                 Wellhead_3116=vrtxs.HE2_ABC_GraphVertex(),

                 PAD_5=vrtxs.HE2_ABC_GraphVertex(),
                 PAD_33=vrtxs.HE2_ABC_GraphVertex(),
                 PAD_34=vrtxs.HE2_ABC_GraphVertex(),
                 PAD_39=vrtxs.HE2_ABC_GraphVertex(),
                 PAD_49=vrtxs.HE2_ABC_GraphVertex(),
                 PAD_57=vrtxs.HE2_ABC_GraphVertex(),
                 intake_pad_5=vrtxs.HE2_ABC_GraphVertex(),
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

    q = DNS_daily_debit * fluid.calc(P_bar=20, T_C=20, X_kgsec=0).CurrentLiquidDensity_kg_m3 / 86400
    outlets = dict(DNS_2=vrtxs.HE2_Boundary_Vertex('P', DNS_pressure))

    G = nx.DiGraph()  # Di = directed
    for k, v in {**inlets, **outlets, **juncs}.items():
        G.add_node(k, obj=v)


    #Продуктивности скважин
    #Куст 5
    productivity_1523 = plasts["PAD_5"]["WELL_1523"]
    productivity_146 = plasts["PAD_5"]["WELL_146"]
    productivity_142 = plasts["PAD_5"]["WELL_142"]
    productivity_1562 = plasts["PAD_5"]["WELL_1562"]

    # Куст 33
    productivity_1385= plasts["PAD_33"]["WELL_1385"]
    productivity_736 = plasts["PAD_33"]["WELL_736"]
    productivity_739 = plasts["PAD_33"]["WELL_739"]
    productivity_1383 = plasts["PAD_33"]["WELL_1383"]
    productivity_738 = plasts["PAD_33"]["WELL_738"]
    productivity_725 = plasts["PAD_33"]["WELL_725"]

    # Куст 34
    productivity_731 = plasts["PAD_34"]["WELL_731"]
    productivity_196 = plasts["PAD_34"]["WELL_196"]
    productivity_734 = plasts["PAD_34"]["WELL_734"]
    productivity_198 = plasts["PAD_34"]["WELL_198"]
    productivity_199 = plasts["PAD_34"]["WELL_199"]
    productivity_197 = plasts["PAD_34"]["WELL_197"]
    productivity_195 = plasts["PAD_34"]["WELL_195"]
    productivity_191 = plasts["PAD_34"]["WELL_191"]
    productivity_729 = plasts["PAD_34"]["WELL_729"]
    productivity_730 = plasts["PAD_34"]["WELL_730"]
    productivity_192 = plasts["PAD_34"]["WELL_192"]
    productivity_148 = plasts["PAD_34"]["WELL_148"]

    # Куст 39
    productivity_3552 = plasts["PAD_39"]["WELL_3552"]
    productivity_617 = plasts["PAD_39"]["WELL_617"]
    productivity_567 = plasts["PAD_39"]["WELL_567"]
    productivity_614 = plasts["PAD_39"]["WELL_614"]
    productivity_619 = plasts["PAD_39"]["WELL_619"]
    productivity_609 = plasts["PAD_39"]["WELL_609"]

    # Куст 39
    productivity_1816 = plasts["PAD_49"]["WELL_1816"]
    productivity_2630 = plasts["PAD_49"]["WELL_2630"]
    productivity_1815 = plasts["PAD_49"]["WELL_1815"]
    productivity_676 = plasts["PAD_49"]["WELL_676"]
    productivity_3270 = plasts["PAD_49"]["WELL_3270"]
    productivity_3266 = plasts["PAD_49"]["WELL_3266"]
    productivity_1814 = plasts["PAD_49"]["WELL_1814"]
    productivity_1817 = plasts["PAD_49"]["WELL_1817"]
    productivity_4532 = plasts["PAD_49"]["WELL_4532"]
    productivity_2631 = plasts["PAD_49"]["WELL_2631"]
    productivity_677 = plasts["PAD_49"]["WELL_677"]
    
    #Куст 57
    productivity_3113 = plasts["PAD_57"]["WELL_3113"]
    productivity_3118 = plasts["PAD_57"]["WELL_3118"]
    productivity_3112 = plasts["PAD_57"]["WELL_3112"]
    productivity_4235 = plasts["PAD_57"]["WELL_4235"]
    productivity_3117 = plasts["PAD_57"]["WELL_3117"]
    productivity_1493 = plasts["PAD_57"]["WELL_1493"]
    productivity_1574 = plasts["PAD_57"]["WELL_1574"]
    productivity_1579 = plasts["PAD_57"]["WELL_1579"]
    productivity_3116 = plasts["PAD_57"]["WELL_3116"]
    
    #Пласты скважин
    #Куст 5
    G.add_edge('PAD_5_well_1523', 'Zaboi_1523', obj=HE2_Plast(productivity=productivity_1523, fluid=fluid))
    G.add_edge('PAD_5_well_146', 'Zaboi_146', obj=HE2_Plast(productivity=productivity_146, fluid=fluid))
    G.add_edge('PAD_5_well_142', 'Zaboi_142', obj=HE2_Plast(productivity=productivity_142, fluid=fluid))
    G.add_edge('PAD_5_well_1562', 'Zaboi_1562', obj=HE2_Plast(productivity=productivity_1562, fluid=fluid))

    # Куст 33
    G.add_edge('PAD_33_well_1385', 'Zaboi_1385', obj=HE2_Plast(productivity=productivity_1385, fluid=fluid))
    G.add_edge('PAD_33_well_736', 'Zaboi_736', obj=HE2_Plast(productivity=productivity_736, fluid=fluid))
    G.add_edge('PAD_33_well_739', 'Zaboi_739', obj=HE2_Plast(productivity=productivity_739, fluid=fluid))
    G.add_edge('PAD_33_well_1383', 'Zaboi_1383', obj=HE2_Plast(productivity=productivity_1383, fluid=fluid))
    G.add_edge('PAD_33_well_738', 'Zaboi_738', obj=HE2_Plast(productivity=productivity_738, fluid=fluid))
    G.add_edge('PAD_33_well_725', 'Zaboi_725', obj=HE2_Plast(productivity=productivity_725, fluid=fluid))

    # Куст 34
    G.add_edge('PAD_34_well_731', 'Zaboi_731', obj=HE2_Plast(productivity=productivity_731, fluid=fluid))
    G.add_edge('PAD_34_well_196', 'Zaboi_196', obj=HE2_Plast(productivity=productivity_196, fluid=fluid))
    G.add_edge('PAD_34_well_734', 'Zaboi_734', obj=HE2_Plast(productivity=productivity_734, fluid=fluid))
    G.add_edge('PAD_34_well_198', 'Zaboi_198', obj=HE2_Plast(productivity=productivity_198, fluid=fluid))
    G.add_edge('PAD_34_well_199', 'Zaboi_199', obj=HE2_Plast(productivity=productivity_199, fluid=fluid))
    G.add_edge('PAD_34_well_197', 'Zaboi_197', obj=HE2_Plast(productivity=productivity_197, fluid=fluid))
    G.add_edge('PAD_34_well_195', 'Zaboi_195', obj=HE2_Plast(productivity=productivity_195, fluid=fluid))
    G.add_edge('PAD_34_well_191', 'Zaboi_191', obj=HE2_Plast(productivity=productivity_191, fluid=fluid))
    G.add_edge('PAD_34_well_729', 'Zaboi_729', obj=HE2_Plast(productivity=productivity_729, fluid=fluid))
    G.add_edge('PAD_34_well_730', 'Zaboi_730', obj=HE2_Plast(productivity=productivity_730, fluid=fluid))
    G.add_edge('PAD_34_well_192', 'Zaboi_192', obj=HE2_Plast(productivity=productivity_192, fluid=fluid))
    G.add_edge('PAD_34_well_148', 'Zaboi_148', obj=HE2_Plast(productivity=productivity_148, fluid=fluid))

    # Куст 39
    G.add_edge('PAD_39_well_3552', 'Zaboi_3552', obj=HE2_Plast(productivity=productivity_3552, fluid=fluid))
    G.add_edge('PAD_39_well_617', 'Zaboi_617', obj=HE2_Plast(productivity=productivity_617, fluid=fluid))
    G.add_edge('PAD_39_well_567', 'Zaboi_567', obj=HE2_Plast(productivity=productivity_567, fluid=fluid))
    G.add_edge('PAD_39_well_614', 'Zaboi_614', obj=HE2_Plast(productivity=productivity_614, fluid=fluid))
    G.add_edge('PAD_39_well_619', 'Zaboi_619', obj=HE2_Plast(productivity=productivity_619, fluid=fluid))
    G.add_edge('PAD_39_well_609', 'Zaboi_609', obj=HE2_Plast(productivity=productivity_609, fluid=fluid))

    # Куст 49
    G.add_edge('PAD_49_well_1816', 'Zaboi_1816', obj=HE2_Plast(productivity=productivity_1816, fluid=fluid))
    G.add_edge('PAD_49_well_2630', 'Zaboi_2630', obj=HE2_Plast(productivity=productivity_2630, fluid=fluid))
    G.add_edge('PAD_49_well_1815', 'Zaboi_1815', obj=HE2_Plast(productivity=productivity_1815, fluid=fluid))
    G.add_edge('PAD_49_well_676', 'Zaboi_676', obj=HE2_Plast(productivity=productivity_676, fluid=fluid))
    G.add_edge('PAD_49_well_3270', 'Zaboi_3270', obj=HE2_Plast(productivity=productivity_3270, fluid=fluid))
    G.add_edge('PAD_49_well_3266', 'Zaboi_3266', obj=HE2_Plast(productivity=productivity_3266, fluid=fluid))
    G.add_edge('PAD_49_well_1814', 'Zaboi_1814', obj=HE2_Plast(productivity=productivity_1814, fluid=fluid))
    G.add_edge('PAD_49_well_1817', 'Zaboi_1817', obj=HE2_Plast(productivity=productivity_1817, fluid=fluid))
    G.add_edge('PAD_49_well_4532', 'Zaboi_4532', obj=HE2_Plast(productivity=productivity_4532, fluid=fluid))
    G.add_edge('PAD_49_well_2631', 'Zaboi_2631', obj=HE2_Plast(productivity=productivity_2631, fluid=fluid))
    G.add_edge('PAD_49_well_677', 'Zaboi_677', obj=HE2_Plast(productivity=productivity_677, fluid=fluid))

    # Куст 57

    G.add_edge('PAD_57_well_3113', 'Zaboi_3113', obj=HE2_Plast(productivity=productivity_3113, fluid=fluid))
    G.add_edge('PAD_57_well_3118', 'Zaboi_3118', obj=HE2_Plast(productivity=productivity_3118, fluid=fluid))
    G.add_edge('PAD_57_well_3112', 'Zaboi_3112', obj=HE2_Plast(productivity=productivity_3112, fluid=fluid))
    G.add_edge('PAD_57_well_4235', 'Zaboi_4235', obj=HE2_Plast(productivity=productivity_4235, fluid=fluid))
    G.add_edge('PAD_57_well_3117', 'Zaboi_3117', obj=HE2_Plast(productivity=productivity_3117, fluid=fluid))
    G.add_edge('PAD_57_well_1493', 'Zaboi_1493', obj=HE2_Plast(productivity=productivity_1493, fluid=fluid))
    G.add_edge('PAD_57_well_1574', 'Zaboi_1574', obj=HE2_Plast(productivity=productivity_1574, fluid=fluid))
    G.add_edge('PAD_57_well_1579', 'Zaboi_1579', obj=HE2_Plast(productivity=productivity_1579, fluid=fluid))
    G.add_edge('PAD_57_well_3116', 'Zaboi_3116', obj=HE2_Plast(productivity=productivity_3116, fluid=fluid))
    
    #Насосы скважин

    #Куст 5
    pump_1523 = pumps["PAD_5"]["WELL_1523"]
    pump_146 = pumps["PAD_5"]["WELL_146"]
    pump_142 = pumps["PAD_5"]["WELL_142"]
    pump_1562 = pumps["PAD_5"]["WELL_1562"]

    #Куст 33
    pump_1385 = pumps["PAD_33"]["WELL_1385"]
    pump_736 = pumps["PAD_33"]["WELL_736"]
    pump_739 = pumps["PAD_33"]["WELL_739"]
    pump_1383 = pumps["PAD_33"]["WELL_1383"]
    pump_738 = pumps["PAD_33"]["WELL_738"]
    pump_725 = pumps["PAD_33"]["WELL_725"]

    # Куст 34
    pump_731 = pumps["PAD_34"]["WELL_731"]
    pump_196 = pumps["PAD_34"]["WELL_196"]
    pump_734 = pumps["PAD_34"]["WELL_734"]
    pump_198 = pumps["PAD_34"]["WELL_198"]
    pump_199 = pumps["PAD_34"]["WELL_199"]
    pump_197 = pumps["PAD_34"]["WELL_197"]
    pump_195 = pumps["PAD_34"]["WELL_195"]
    pump_191 = pumps["PAD_34"]["WELL_191"]
    pump_729 = pumps["PAD_34"]["WELL_729"]
    pump_730 = pumps["PAD_34"]["WELL_730"]
    pump_192 = pumps["PAD_34"]["WELL_192"]
    pump_148 = pumps["PAD_34"]["WELL_148"]
    
    # Куст 39
    pump_3552 = pumps["PAD_39"]["WELL_3552"]
    pump_617 = pumps["PAD_39"]["WELL_617"]
    pump_567 = pumps["PAD_39"]["WELL_567"]
    pump_614 = pumps["PAD_39"]["WELL_614"]
    pump_619 = pumps["PAD_39"]["WELL_619"]
    pump_609 = pumps["PAD_39"]["WELL_609"]

    # Куст 49
    pump_1816 = pumps["PAD_49"]["WELL_1816"]
    pump_2630 = pumps["PAD_49"]["WELL_2630"]
    pump_1815 = pumps["PAD_49"]["WELL_1815"]
    pump_676 = pumps["PAD_49"]["WELL_676"]
    pump_3270 = pumps["PAD_49"]["WELL_3270"]
    pump_3266 = pumps["PAD_49"]["WELL_3266"]
    pump_1814 = pumps["PAD_49"]["WELL_1814"]
    pump_1817 = pumps["PAD_49"]["WELL_1817"]
    pump_4532 = pumps["PAD_49"]["WELL_4532"]
    pump_2631 = pumps["PAD_49"]["WELL_2631"]
    pump_677 = pumps["PAD_49"]["WELL_677"]

    # Куст 57
    pump_3113 = pumps["PAD_57"]["WELL_3113"]
    pump_3118 = pumps["PAD_57"]["WELL_3118"]
    pump_3112 = pumps["PAD_57"]["WELL_3112"]
    pump_4235 = pumps["PAD_57"]["WELL_4235"]
    pump_3117 = pumps["PAD_57"]["WELL_3117"]
    pump_1493 = pumps["PAD_57"]["WELL_1493"]
    pump_1574 = pumps["PAD_57"]["WELL_1574"]
    pump_1579 = pumps["PAD_57"]["WELL_1579"]
    pump_3116 = pumps["PAD_57"]["WELL_3116"]
    
    #Куст 5
    G.add_edge('Pump_intake_1523', 'Pump_outlet_1523',obj=create_HE2_WellPump_instance_from_dataframe(full_HPX=pump_curves, model=pump_1523[0], fluid=fluid,
                                                             frequency=pump_1523[1]))
    G.add_edge('Pump_intake_146', 'Pump_outlet_146',obj=create_HE2_WellPump_instance_from_dataframe(full_HPX=pump_curves, model=pump_146[0], fluid=fluid,
                                                             frequency=pump_146[1]))
    G.add_edge('Pump_intake_142', 'Pump_outlet_142',obj=create_HE2_WellPump_instance_from_dataframe(full_HPX=pump_curves, model=pump_142[0], fluid=fluid,
                                                             frequency=pump_142[1]))
    G.add_edge('Pump_intake_1562', 'Pump_outlet_1562',obj=create_HE2_WellPump_instance_from_dataframe(full_HPX=pump_curves, model=pump_1562[0], fluid=fluid,
                                frequency=pump_1562[1]))
    #Куст 33
    G.add_edge('Pump_intake_1385', 'Pump_outlet_1385',obj=create_HE2_WellPump_instance_from_dataframe(full_HPX=pump_curves, model=pump_1385[0], fluid=fluid,
                                frequency=pump_1385[1]))
    G.add_edge('Pump_intake_736', 'Pump_outlet_736',obj=create_HE2_WellPump_instance_from_dataframe(full_HPX=pump_curves, model=pump_736[0], fluid=fluid,
                                frequency=pump_736[1]))
    G.add_edge('Pump_intake_739', 'Pump_outlet_739',obj=create_HE2_WellPump_instance_from_dataframe(full_HPX=pump_curves, model=pump_739[0], fluid=fluid,
                                frequency=pump_739[1]))
    G.add_edge('Pump_intake_1383', 'Pump_outlet_1383',obj=create_HE2_WellPump_instance_from_dataframe(full_HPX=pump_curves, model=pump_1383[0], fluid=fluid,
                                frequency=pump_1383[1]))
    G.add_edge('Pump_intake_738', 'Pump_outlet_738',obj=create_HE2_WellPump_instance_from_dataframe(full_HPX=pump_curves, model=pump_738[0], fluid=fluid,
                                frequency=pump_738[1]))
    G.add_edge('Pump_intake_725', 'Pump_outlet_725',obj=create_HE2_WellPump_instance_from_dataframe(full_HPX=pump_curves, model=pump_725[0], fluid=fluid,
                                frequency=pump_725[1]))

    #Куст 34
    G.add_edge('Pump_intake_731', 'Pump_outlet_731',obj=create_HE2_WellPump_instance_from_dataframe(full_HPX=pump_curves, model=pump_731[0], fluid=fluid,
                                frequency=pump_731[1]))
    G.add_edge('Pump_intake_196', 'Pump_outlet_196',obj=create_HE2_WellPump_instance_from_dataframe(full_HPX=pump_curves, model=pump_196[0], fluid=fluid,
                                frequency=pump_196[1]))
    G.add_edge('Pump_intake_734', 'Pump_outlet_734',obj=create_HE2_WellPump_instance_from_dataframe(full_HPX=pump_curves, model=pump_734[0], fluid=fluid,
                                frequency=pump_734[1]))
    G.add_edge('Pump_intake_198', 'Pump_outlet_198',obj=create_HE2_WellPump_instance_from_dataframe(full_HPX=pump_curves, model=pump_198[0], fluid=fluid,
                                frequency=pump_198[1]))
    G.add_edge('Pump_intake_199', 'Pump_outlet_199',obj=create_HE2_WellPump_instance_from_dataframe(full_HPX=pump_curves, model=pump_199[0], fluid=fluid,
                                frequency=pump_199[1]))
    G.add_edge('Pump_intake_197', 'Pump_outlet_197',obj=create_HE2_WellPump_instance_from_dataframe(full_HPX=pump_curves, model=pump_197[0], fluid=fluid,
                                frequency=pump_197[1]))
    G.add_edge('Pump_intake_195', 'Pump_outlet_195',obj=create_HE2_WellPump_instance_from_dataframe(full_HPX=pump_curves, model=pump_195[0], fluid=fluid,
                                frequency=pump_195[1]))
    G.add_edge('Pump_intake_191', 'Pump_outlet_191',obj=create_HE2_WellPump_instance_from_dataframe(full_HPX=pump_curves, model=pump_191[0], fluid=fluid,
                                frequency=pump_191[1]))
    G.add_edge('Pump_intake_729', 'Pump_outlet_729',obj=create_HE2_WellPump_instance_from_dataframe(full_HPX=pump_curves, model=pump_729[0], fluid=fluid,
                                frequency=pump_729[1]))
    G.add_edge('Pump_intake_730', 'Pump_outlet_730',obj=create_HE2_WellPump_instance_from_dataframe(full_HPX=pump_curves, model=pump_730[0], fluid=fluid,
                                frequency=pump_730[1]))
    G.add_edge('Pump_intake_192', 'Pump_outlet_192',obj=create_HE2_WellPump_instance_from_dataframe(full_HPX=pump_curves, model=pump_192[0], fluid=fluid,
                                frequency=pump_192[1]))
    G.add_edge('Pump_intake_148', 'Pump_outlet_148',obj=create_HE2_WellPump_instance_from_dataframe(full_HPX=pump_curves, model=pump_148[0], fluid=fluid,
                                frequency=pump_148[1]))

    #Куст 39
    G.add_edge('Pump_intake_3552', 'Pump_outlet_3552',obj=create_HE2_WellPump_instance_from_dataframe(full_HPX=pump_curves, model=pump_3552[0], fluid=fluid,
                                frequency=pump_3552[1]))
    G.add_edge('Pump_intake_617', 'Pump_outlet_617',obj=create_HE2_WellPump_instance_from_dataframe(full_HPX=pump_curves, model=pump_617[0], fluid=fluid,
                                frequency=pump_617[1]))
    G.add_edge('Pump_intake_567', 'Pump_outlet_567',obj=create_HE2_WellPump_instance_from_dataframe(full_HPX=pump_curves, model=pump_567[0], fluid=fluid,
                                frequency=pump_567[1]))
    G.add_edge('Pump_intake_614', 'Pump_outlet_614',obj=create_HE2_WellPump_instance_from_dataframe(full_HPX=pump_curves, model=pump_614[0], fluid=fluid,
                                frequency=pump_614[1]))
    G.add_edge('Pump_intake_619', 'Pump_outlet_619',obj=create_HE2_WellPump_instance_from_dataframe(full_HPX=pump_curves, model=pump_619[0], fluid=fluid,
                                frequency=pump_619[1]))
    G.add_edge('Pump_intake_609', 'Pump_outlet_609',obj=create_HE2_WellPump_instance_from_dataframe(full_HPX=pump_curves, model=pump_609[0], fluid=fluid,
                                frequency=pump_609[1]))
    
    # Куст 49
    G.add_edge('Pump_intake_1816', 'Pump_outlet_1816',
               obj=create_HE2_WellPump_instance_from_dataframe(full_HPX=pump_curves, model=pump_1816[0], fluid=fluid,frequency=pump_1816[1]))
    G.add_edge('Pump_intake_2630', 'Pump_outlet_2630',
               obj=create_HE2_WellPump_instance_from_dataframe(full_HPX=pump_curves, model=pump_2630[0], fluid=fluid, frequency=pump_2630[1]))
    G.add_edge('Pump_intake_1815', 'Pump_outlet_1815',
               obj=create_HE2_WellPump_instance_from_dataframe(full_HPX=pump_curves, model=pump_1815[0], fluid=fluid, frequency=pump_1815[1]))
    G.add_edge('Pump_intake_676', 'Pump_outlet_676',
               obj=create_HE2_WellPump_instance_from_dataframe(full_HPX=pump_curves, model=pump_676[0], fluid=fluid,frequency=pump_676[1]))
    G.add_edge('Pump_intake_3270', 'Pump_outlet_3270',
               obj=create_HE2_WellPump_instance_from_dataframe(full_HPX=pump_curves, model=pump_3270[0], fluid=fluid, frequency=pump_3270[1]))
    G.add_edge('Pump_intake_3266', 'Pump_outlet_3266',
               obj=create_HE2_WellPump_instance_from_dataframe(full_HPX=pump_curves, model=pump_3266[0], fluid=fluid,frequency=pump_3266[1]))
    G.add_edge('Pump_intake_1814', 'Pump_outlet_1814',
               obj=create_HE2_WellPump_instance_from_dataframe(full_HPX=pump_curves, model=pump_1814[0], fluid=fluid, frequency=pump_1814[1]))
    G.add_edge('Pump_intake_1817', 'Pump_outlet_1817',
               obj=create_HE2_WellPump_instance_from_dataframe(full_HPX=pump_curves, model=pump_1817[0], fluid=fluid, frequency=pump_1817[1]))
    G.add_edge('Pump_intake_4532', 'Pump_outlet_4532',
               obj=create_HE2_WellPump_instance_from_dataframe(full_HPX=pump_curves, model=pump_4532[0], fluid=fluid,frequency=pump_4532[1]))
    G.add_edge('Pump_intake_2631', 'Pump_outlet_2631',
               obj=create_HE2_WellPump_instance_from_dataframe(full_HPX=pump_curves, model=pump_2631[0], fluid=fluid, frequency=pump_2631[1]))
    G.add_edge('Pump_intake_677', 'Pump_outlet_677',
               obj=create_HE2_WellPump_instance_from_dataframe(full_HPX=pump_curves, model=pump_677[0], fluid=fluid, frequency=pump_677[1]))

    # Куст 57
    G.add_edge('Pump_intake_3113', 'Pump_outlet_3113',obj=create_HE2_WellPump_instance_from_dataframe(full_HPX=pump_curves, model=pump_3113[0], fluid=fluid,
                                frequency=pump_3113[1]))
    G.add_edge('Pump_intake_3118', 'Pump_outlet_3118',obj=create_HE2_WellPump_instance_from_dataframe(full_HPX=pump_curves, model=pump_3118[0], fluid=fluid,
                                frequency=pump_3118[1]))
    G.add_edge('Pump_intake_3112', 'Pump_outlet_3112',obj=create_HE2_WellPump_instance_from_dataframe(full_HPX=pump_curves, model=pump_3112[0], fluid=fluid,
                                frequency=pump_3112[1]))
    G.add_edge('Pump_intake_4235', 'Pump_outlet_4235',obj=create_HE2_WellPump_instance_from_dataframe(full_HPX=pump_curves, model=pump_4235[0], fluid=fluid,
                                frequency=pump_4235[1]))
    G.add_edge('Pump_intake_3117', 'Pump_outlet_3117',obj=create_HE2_WellPump_instance_from_dataframe(full_HPX=pump_curves, model=pump_3117[0], fluid=fluid,
                                frequency=pump_3117[1]))
    G.add_edge('Pump_intake_1493', 'Pump_outlet_1493',obj=create_HE2_WellPump_instance_from_dataframe(full_HPX=pump_curves, model=pump_1493[0], fluid=fluid,
                                frequency=pump_1493[1]))
    G.add_edge('Pump_intake_1574', 'Pump_outlet_1574',obj=create_HE2_WellPump_instance_from_dataframe(full_HPX=pump_curves, model=pump_1574[0], fluid=fluid,
                                frequency=pump_1574[1]))
    G.add_edge('Pump_intake_1579', 'Pump_outlet_1579',obj=create_HE2_WellPump_instance_from_dataframe(full_HPX=pump_curves, model=pump_1579[0], fluid=fluid,
                                frequency=pump_1579[1]))
    G.add_edge('Pump_intake_3116', 'Pump_outlet_3116',obj=create_HE2_WellPump_instance_from_dataframe(full_HPX=pump_curves, model=pump_3116[0], fluid=fluid,
                                frequency=pump_3116[1]))
    
    #Инклинометрия скважин
    #Куст 5
    G.add_edge('Pump_outlet_1523', 'Wellhead_1523', obj=HE2_OilPipe([32.6], [2561], [0.125 * real_diam_coefficient], [roughness], fluid))
    G.add_edge('Zaboi_1523', 'Pump_intake_1523',obj=HE2_OilPipe([267.3], [202.7], [0.143 * real_diam_coefficient], [5 * roughness], fluid))

    G.add_edge('Pump_outlet_146', 'Wellhead_146', obj=HE2_OilPipe([ 47.569], [2575.43], [0.125 * real_diam_coefficient], [roughness], fluid))
    G.add_edge('Zaboi_146', 'Pump_intake_146',obj=HE2_OilPipe([161.39], [210.56], [0.143 * real_diam_coefficient], [5 * roughness], fluid))


    G.add_edge('Pump_outlet_142', 'Wellhead_142', obj=HE2_OilPipe([324.89], [2511.1], [0.125 * real_diam_coefficient], [roughness], fluid))
    G.add_edge('Zaboi_142', 'Pump_intake_142',obj=HE2_OilPipe([242.24], [249.7599], [0.143 * real_diam_coefficient], [5 * roughness], fluid))

    G.add_edge('Pump_outlet_1562', 'Wellhead_1562', obj=HE2_OilPipe([549.92], [2392.08], [0.125 * real_diam_coefficient], [roughness], fluid))
    G.add_edge('Zaboi_1562', 'Pump_intake_1562',obj=HE2_OilPipe([36.4499], [325.55], [0.143 * real_diam_coefficient], [5 * roughness], fluid))

    G.add_edge('Wellhead_1562', 'PAD_5',obj=HE2_OilPipe([100], [0], [0.125 * real_diam_coefficient], [roughness], fluid))
    G.add_edge('Wellhead_1523', 'PAD_5',obj=HE2_OilPipe([100], [0], [0.125 * real_diam_coefficient], [roughness], fluid))
    G.add_edge('Wellhead_146', 'PAD_5',obj=HE2_OilPipe([100], [0], [0.125 * real_diam_coefficient], [roughness], fluid))
    G.add_edge('Wellhead_142', 'PAD_5',obj=HE2_OilPipe([100], [0], [0.125 * real_diam_coefficient], [roughness], fluid))


    #Куст 33
    G.add_edge('Pump_outlet_1385', 'Wellhead_1385', obj=HE2_OilPipe([395.8], [2139.2], [0.125 * real_diam_coefficient], [roughness], fluid))
    G.add_edge('Zaboi_1385', 'Pump_intake_1385',obj=HE2_OilPipe([80], [584], [0.143 * real_diam_coefficient], [5 * roughness], fluid))

    G.add_edge('Pump_outlet_736', 'Wellhead_736', obj=HE2_OilPipe([216], [2532], [0.125 * real_diam_coefficient], [roughness], fluid))
    G.add_edge('Zaboi_736', 'Pump_intake_736',obj=HE2_OilPipe([448], [155], [0.143 * real_diam_coefficient], [5 * roughness], fluid))

    G.add_edge('Pump_outlet_739', 'Wellhead_739', obj=HE2_OilPipe([323.45], [2356.55], [0.125 * real_diam_coefficient], [roughness], fluid))
    G.add_edge('Zaboi_739', 'Pump_intake_739',obj=HE2_OilPipe([396.35], [290.6], [0.143 * real_diam_coefficient], [5 * roughness], fluid))

    G.add_edge('Pump_outlet_1383', 'Wellhead_1383', obj=HE2_OilPipe([526], [2353.99], [0.125 * real_diam_coefficient], [roughness], fluid))
    G.add_edge('Zaboi_1383', 'Pump_intake_1383',obj=HE2_OilPipe([30.16], [377.84], [0.143 * real_diam_coefficient], [5 * roughness], fluid))

    G.add_edge('Pump_outlet_738', 'Wellhead_738', obj=HE2_OilPipe([243.83], [2128.17], [0.125 * real_diam_coefficient], [roughness], fluid))
    G.add_edge('Zaboi_738', 'Pump_intake_738',obj=HE2_OilPipe([402.45], [526.85], [0.143 * real_diam_coefficient], [5 * roughness], fluid))

    G.add_edge('Pump_outlet_725', 'Wellhead_725', obj=HE2_OilPipe([495.76], [2329.24], [0.125 * real_diam_coefficient], [roughness], fluid))
    G.add_edge('Zaboi_725', 'Pump_intake_725',obj=HE2_OilPipe([28.7], [384.3], [0.143 * real_diam_coefficient], [5 * roughness], fluid))

    G.add_edge('Wellhead_1385', 'PAD_33',obj=HE2_OilPipe([100], [0], [0.125 * real_diam_coefficient], [roughness], fluid))
    G.add_edge('Wellhead_736', 'PAD_33', obj=HE2_OilPipe([100], [0], [0.125 * real_diam_coefficient], [roughness], fluid))
    G.add_edge('Wellhead_739', 'PAD_33', obj=HE2_OilPipe([100], [0], [0.125 * real_diam_coefficient], [roughness], fluid))
    G.add_edge('Wellhead_1383', 'PAD_33', obj=HE2_OilPipe([100], [0], [0.125 * real_diam_coefficient], [roughness], fluid))
    G.add_edge('Wellhead_738', 'PAD_33', obj=HE2_OilPipe([100], [0], [0.125 * real_diam_coefficient], [roughness], fluid))
    G.add_edge('Wellhead_725', 'PAD_33', obj=HE2_OilPipe([100], [0], [0.125 * real_diam_coefficient], [roughness], fluid))

    #Куст 34
    G.add_edge('Pump_outlet_731', 'Wellhead_731', obj=HE2_OilPipe([285.7], [2414.27], [0.125 * real_diam_coefficient], [roughness], fluid))
    G.add_edge('Zaboi_731', 'Pump_intake_731',obj=HE2_OilPipe([13.5], [330], [0.143 * real_diam_coefficient], [5 * roughness], fluid))

    G.add_edge('Pump_outlet_196', 'Wellhead_196', obj=HE2_OilPipe([172.6], [2504.43], [0.125 * real_diam_coefficient], [roughness], fluid))
    G.add_edge('Zaboi_196', 'Pump_intake_196',obj=HE2_OilPipe([7.38], [245.12], [0.143 * real_diam_coefficient], [5 * roughness], fluid))

    G.add_edge('Pump_outlet_734', 'Wellhead_734', obj=HE2_OilPipe([32.21], [2362.8], [0.125 * real_diam_coefficient], [roughness], fluid))
    G.add_edge('Zaboi_734', 'Pump_intake_734',obj=HE2_OilPipe([282], [286.88], [0.143 * real_diam_coefficient], [5 * roughness], fluid))

    G.add_edge('Pump_outlet_198', 'Wellhead_198', obj=HE2_OilPipe([244], [2586], [0.125 * real_diam_coefficient], [roughness], fluid))
    G.add_edge('Zaboi_198', 'Pump_intake_198',obj=HE2_OilPipe([4.82], [169.18], [0.143 * real_diam_coefficient], [5 * roughness], fluid))

    G.add_edge('Pump_outlet_199', 'Wellhead_199', obj=HE2_OilPipe([228.11], [2520.89], [0.125 * real_diam_coefficient], [roughness], fluid))
    G.add_edge('Zaboi_199', 'Pump_intake_199',obj=HE2_OilPipe([198.34], [241.56], [0.143 * real_diam_coefficient], [5 * roughness], fluid))

    G.add_edge('Pump_outlet_197', 'Wellhead_197', obj=HE2_OilPipe([668.11], [2345.89], [0.125 * real_diam_coefficient], [roughness], fluid))
    G.add_edge('Zaboi_197', 'Pump_intake_197',obj=HE2_OilPipe([52.91], [403.09], [0.143 * real_diam_coefficient], [5 * roughness], fluid))

    G.add_edge('Pump_outlet_195', 'Wellhead_195', obj=HE2_OilPipe([610.82], [2372.96], [0.125 * real_diam_coefficient], [roughness], fluid))
    G.add_edge('Zaboi_195', 'Pump_intake_195',obj=HE2_OilPipe([50.06], [357.93], [0.143 * real_diam_coefficient], [5 * roughness], fluid))

    G.add_edge('Pump_outlet_191', 'Wellhead_191', obj=HE2_OilPipe([254.37], [2585.63], [0.125 * real_diam_coefficient], [roughness], fluid))
    G.add_edge('Zaboi_191', 'Pump_intake_191',obj=HE2_OilPipe([251.19], [166.221], [0.143 * real_diam_coefficient], [5 * roughness], fluid))

    G.add_edge('Pump_outlet_729', 'Wellhead_729', obj=HE2_OilPipe([92.78], [2345.22], [0.125 * real_diam_coefficient], [roughness], fluid))
    G.add_edge('Zaboi_729', 'Pump_intake_729',obj=HE2_OilPipe([247.95], [399.048], [0.143 * real_diam_coefficient], [5 * roughness], fluid))

    G.add_edge('Pump_outlet_730', 'Wellhead_730', obj=HE2_OilPipe([75.76], [2218.24], [0.125 * real_diam_coefficient], [roughness], fluid))
    G.add_edge('Zaboi_730', 'Pump_intake_730',obj=HE2_OilPipe([99.72], [532.78], [0.143 * real_diam_coefficient], [5 * roughness], fluid))

    G.add_edge('Pump_outlet_192', 'Wellhead_192', obj=HE2_OilPipe([337.39], [2177.61], [0.125 * real_diam_coefficient], [roughness], fluid))
    G.add_edge('Zaboi_192', 'Pump_intake_192',obj=HE2_OilPipe([255.33], [531.87], [0.143 * real_diam_coefficient], [5 * roughness], fluid))

    G.add_edge('Pump_outlet_148', 'Wellhead_148', obj=HE2_OilPipe([337.39], [2177.61], [0.125 * real_diam_coefficient], [roughness], fluid))
    G.add_edge('Zaboi_148', 'Pump_intake_148',obj=HE2_OilPipe([255.33], [531.87], [0.143 * real_diam_coefficient], [5 * roughness], fluid))


    G.add_edge('Wellhead_731', 'PAD_34',obj=HE2_OilPipe([100], [0], [0.125 * real_diam_coefficient], [roughness], fluid))
    G.add_edge('Wellhead_196', 'PAD_34',obj=HE2_OilPipe([100], [0], [0.125 * real_diam_coefficient], [roughness], fluid))
    G.add_edge('Wellhead_734', 'PAD_34',obj=HE2_OilPipe([100], [0], [0.125 * real_diam_coefficient], [roughness], fluid))
    G.add_edge('Wellhead_198', 'PAD_34',obj=HE2_OilPipe([100], [0], [0.125 * real_diam_coefficient], [roughness], fluid))
    G.add_edge('Wellhead_199', 'PAD_34',obj=HE2_OilPipe([100], [0], [0.125 * real_diam_coefficient], [roughness], fluid))
    G.add_edge('Wellhead_197', 'PAD_34',obj=HE2_OilPipe([100], [0], [0.125 * real_diam_coefficient], [roughness], fluid))
    G.add_edge('Wellhead_195', 'PAD_34',obj=HE2_OilPipe([100], [0], [0.125 * real_diam_coefficient], [roughness], fluid))
    G.add_edge('Wellhead_191', 'PAD_34',obj=HE2_OilPipe([100], [0], [0.125 * real_diam_coefficient], [roughness], fluid))
    G.add_edge('Wellhead_729', 'PAD_34',obj=HE2_OilPipe([100], [0], [0.125 * real_diam_coefficient], [roughness], fluid))
    G.add_edge('Wellhead_730', 'PAD_34',obj=HE2_OilPipe([100], [0], [0.125 * real_diam_coefficient], [roughness], fluid))
    G.add_edge('Wellhead_192', 'PAD_34',obj=HE2_OilPipe([100], [0], [0.125 * real_diam_coefficient], [roughness], fluid))
    G.add_edge('Wellhead_148', 'PAD_34',obj=HE2_OilPipe([100], [0], [0.125 * real_diam_coefficient], [roughness], fluid))

    # Куст 39
    G.add_edge('Pump_outlet_3552', 'Wellhead_3552', obj=HE2_OilPipe([130.37], [2169.63], [0.125 * real_diam_coefficient], [roughness], fluid))
    G.add_edge('Zaboi_3552', 'Pump_intake_3552',obj=HE2_OilPipe([155.83], [687.06], [0.143 * real_diam_coefficient], [5 * roughness], fluid))
    
    G.add_edge('Pump_outlet_617', 'Wellhead_617', obj=HE2_OilPipe([461.4], [2409.6], [0.125 * real_diam_coefficient], [roughness], fluid))
    G.add_edge('Zaboi_617', 'Pump_intake_617',obj=HE2_OilPipe([8], [333], [0.143 * real_diam_coefficient], [5 * roughness], fluid))
    
    G.add_edge('Pump_outlet_567', 'Wellhead_567', obj=HE2_OilPipe([674.06], [2165.94], [0.125 * real_diam_coefficient], [roughness], fluid))
    G.add_edge('Zaboi_567', 'Pump_intake_567',obj=HE2_OilPipe([112.83], [554.17], [0.143 * real_diam_coefficient], [5 * roughness], fluid))
    
    G.add_edge('Pump_outlet_614', 'Wellhead_614', obj=HE2_OilPipe([128.89], [2529.07], [0.125 * real_diam_coefficient], [roughness], fluid))
    G.add_edge('Zaboi_614', 'Pump_intake_614',obj=HE2_OilPipe([156.93], [215.51], [0.143 * real_diam_coefficient], [5 * roughness], fluid))
    
    G.add_edge('Pump_outlet_619', 'Wellhead_619', obj=HE2_OilPipe([269.27], [2343.73], [0.125 * real_diam_coefficient], [roughness], fluid))
    G.add_edge('Zaboi_619', 'Pump_intake_619',obj=HE2_OilPipe([110.95], [384.75], [0.143 * real_diam_coefficient], [5 * roughness], fluid))
    
    G.add_edge('Pump_outlet_609', 'Wellhead_609', obj=HE2_OilPipe([117.42], [2291.58], [0.125 * real_diam_coefficient], [roughness], fluid))
    G.add_edge('Zaboi_609', 'Pump_intake_609',obj=HE2_OilPipe([202.87], [535.23], [0.143 * real_diam_coefficient], [5 * roughness], fluid))

    G.add_edge('Wellhead_3552', 'PAD_39', obj=HE2_OilPipe([100], [0], [0.125 * real_diam_coefficient], [roughness], fluid))
    G.add_edge('Wellhead_617', 'PAD_39', obj=HE2_OilPipe([100], [0], [0.125 * real_diam_coefficient], [roughness], fluid))
    G.add_edge('Wellhead_567', 'PAD_39', obj=HE2_OilPipe([100], [0], [0.125 * real_diam_coefficient], [roughness], fluid))
    G.add_edge('Wellhead_614', 'PAD_39', obj=HE2_OilPipe([100], [0], [0.125 * real_diam_coefficient], [roughness], fluid))
    G.add_edge('Wellhead_619', 'PAD_39', obj=HE2_OilPipe([100], [0], [0.125 * real_diam_coefficient], [roughness], fluid))
    G.add_edge('Wellhead_609', 'PAD_39', obj=HE2_OilPipe([100], [0], [0.125 * real_diam_coefficient], [roughness], fluid))
    
    # Куст 49    

    G.add_edge('Pump_outlet_1816', 'Wellhead_1816', obj=HE2_OilPipe([61.69], [2508.31], [0.125 * real_diam_coefficient], [roughness], fluid))
    G.add_edge('Zaboi_1816', 'Pump_intake_1816',obj=HE2_OilPipe([811.92], [67.67], [0.143 * real_diam_coefficient], [5 * roughness], fluid))


    G.add_edge('Pump_outlet_2630', 'Wellhead_2630', obj=HE2_OilPipe([419.86], [2370.14], [0.125 * real_diam_coefficient], [roughness], fluid))
    G.add_edge('Zaboi_2630', 'Pump_intake_2630',obj=HE2_OilPipe([570.43], [310.38], [0.143 * real_diam_coefficient], [5 * roughness], fluid))

    G.add_edge('Pump_outlet_1815', 'Wellhead_1815', obj=HE2_OilPipe([113.77], [2516.23], [0.125 * real_diam_coefficient], [roughness], fluid))
    G.add_edge('Zaboi_1815', 'Pump_intake_1815',obj=HE2_OilPipe([482.82], [333.82], [0.143 * real_diam_coefficient], [5 * roughness], fluid))


    G.add_edge('Pump_outlet_676', 'Wellhead_676', obj=HE2_OilPipe([248.29], [2433.71], [0.125 * real_diam_coefficient], [roughness], fluid))
    G.add_edge('Zaboi_676', 'Pump_intake_676',obj=HE2_OilPipe([537.87], [336.76], [0.143 * real_diam_coefficient], [5 * roughness], fluid))


    G.add_edge('Pump_outlet_3270', 'Wellhead_3270', obj=HE2_OilPipe([310.83], [2080.17], [0.125 * real_diam_coefficient], [roughness], fluid))
    G.add_edge('Zaboi_3270', 'Pump_intake_3270',obj=HE2_OilPipe([315.89], [637.11], [0.143 * real_diam_coefficient], [5 * roughness], fluid))

    G.add_edge('Pump_outlet_3266', 'Wellhead_3266', obj=HE2_OilPipe([119.76], [2500.24], [0.125 * real_diam_coefficient], [roughness], fluid))
    G.add_edge('Zaboi_3266', 'Pump_intake_3266',obj=HE2_OilPipe([280.05], [489.3], [0.143 * real_diam_coefficient], [5 * roughness], fluid))


    G.add_edge('Pump_outlet_1814', 'Wellhead_1814', obj=HE2_OilPipe([205.63], [2455.37], [0.125 * real_diam_coefficient], [roughness], fluid))
    G.add_edge('Zaboi_1814', 'Pump_intake_1814',obj=HE2_OilPipe([186], [600.24], [0.143 * real_diam_coefficient], [5 * roughness], fluid))


    G.add_edge('Pump_outlet_1817', 'Wellhead_1817', obj=HE2_OilPipe([324.58], [2443.42], [0.125 * real_diam_coefficient], [roughness], fluid))
    G.add_edge('Zaboi_1817', 'Pump_intake_1817',obj=HE2_OilPipe([690.26], [244.53], [0.143 * real_diam_coefficient], [5 * roughness], fluid))

    G.add_edge('Pump_outlet_4532', 'Wellhead_4532', obj=HE2_OilPipe([363.75], [2482.25], [0.125 * real_diam_coefficient], [roughness], fluid))
    G.add_edge('Zaboi_4532', 'Pump_intake_4532',obj=HE2_OilPipe([403.62], [509.5], [0.143 * real_diam_coefficient], [5 * roughness], fluid))

    G.add_edge('Pump_outlet_2631', 'Wellhead_2631', obj=HE2_OilPipe([67.5], [2394.5], [0.125 * real_diam_coefficient], [roughness], fluid))
    G.add_edge('Zaboi_2631', 'Pump_intake_2631',obj=HE2_OilPipe([626.2], [278.96], [0.143 * real_diam_coefficient], [5 * roughness], fluid))


    G.add_edge('Pump_outlet_677', 'Wellhead_677', obj=HE2_OilPipe([248.2], [2351.8], [0.125 * real_diam_coefficient], [roughness], fluid))
    G.add_edge('Zaboi_677', 'Pump_intake_677',obj=HE2_OilPipe([728.15], [217.82], [0.143 * real_diam_coefficient], [5 * roughness], fluid))
    
    G.add_edge('Wellhead_1816', 'PAD_49', obj=HE2_OilPipe([100], [0], [0.125 * real_diam_coefficient], [roughness], fluid))
    G.add_edge('Wellhead_2630', 'PAD_49', obj=HE2_OilPipe([100], [0], [0.125 * real_diam_coefficient], [roughness], fluid))
    G.add_edge('Wellhead_1815', 'PAD_49', obj=HE2_OilPipe([100], [0], [0.125 * real_diam_coefficient], [roughness], fluid))
    G.add_edge('Wellhead_676', 'PAD_49', obj=HE2_OilPipe([100], [0], [0.125 * real_diam_coefficient], [roughness], fluid))
    G.add_edge('Wellhead_3270', 'PAD_49', obj=HE2_OilPipe([100], [0], [0.125 * real_diam_coefficient], [roughness], fluid))
    G.add_edge('Wellhead_3266', 'PAD_49', obj=HE2_OilPipe([100], [0], [0.125 * real_diam_coefficient], [roughness], fluid))
    G.add_edge('Wellhead_1814', 'PAD_49', obj=HE2_OilPipe([100], [0], [0.125 * real_diam_coefficient], [roughness], fluid))
    G.add_edge('Wellhead_1817', 'PAD_49', obj=HE2_OilPipe([100], [0], [0.125 * real_diam_coefficient], [roughness], fluid))
    G.add_edge('Wellhead_4532', 'PAD_49', obj=HE2_OilPipe([100], [0], [0.125 * real_diam_coefficient], [roughness], fluid))
    G.add_edge('Wellhead_2631', 'PAD_49', obj=HE2_OilPipe([100], [0], [0.125 * real_diam_coefficient], [roughness], fluid))
    G.add_edge('Wellhead_677', 'PAD_49', obj=HE2_OilPipe([100], [0], [0.125 * real_diam_coefficient], [roughness], fluid))
    # Куст 57

    G.add_edge('Pump_outlet_3113', 'Wellhead_3113', obj=HE2_OilPipe([382.13], [2367.87], [0.125 * real_diam_coefficient], [roughness], fluid))
    G.add_edge('Zaboi_3113', 'Pump_intake_3113',obj=HE2_OilPipe([652.87], [230.96], [0.143 * real_diam_coefficient], [5 * roughness], fluid))

    G.add_edge('Pump_outlet_3118', 'Wellhead_3118', obj=HE2_OilPipe([152.43], [2458.57], [0.125 * real_diam_coefficient], [roughness], fluid))
    G.add_edge('Zaboi_3118', 'Pump_intake_3118',obj=HE2_OilPipe([670.7], [109], [0.143 * real_diam_coefficient], [5 * roughness], fluid))

    G.add_edge('Pump_outlet_3112', 'Wellhead_3112', obj=HE2_OilPipe([55.74], [2481.26], [0.125 * real_diam_coefficient], [roughness], fluid))
    G.add_edge('Zaboi_3112', 'Pump_intake_3112',obj=HE2_OilPipe([753.03], [-3], [0.143 * real_diam_coefficient], [5 * roughness], fluid))

    G.add_edge('Pump_outlet_4235', 'Wellhead_4235', obj=HE2_OilPipe([351.06], [2514.94], [0.125 * real_diam_coefficient], [roughness], fluid))
    G.add_edge('Zaboi_4235', 'Pump_intake_4235',obj=HE2_OilPipe([511.21], [253.89], [0.143 * real_diam_coefficient], [5 * roughness], fluid))

    G.add_edge('Pump_outlet_3117', 'Wellhead_3117', obj=HE2_OilPipe([117.8], [2477.2], [0.125 * real_diam_coefficient], [roughness], fluid))
    G.add_edge('Zaboi_3117', 'Pump_intake_3117',obj=HE2_OilPipe([527.58], [255.96], [0.143 * real_diam_coefficient], [5 * roughness], fluid))

    G.add_edge('Pump_outlet_1493', 'Wellhead_1493', obj=HE2_OilPipe([403.96], [2442.04], [0.125 * real_diam_coefficient], [roughness], fluid))
    G.add_edge('Zaboi_1493', 'Pump_intake_1493',obj=HE2_OilPipe([409.48], [399.746], [0.143 * real_diam_coefficient], [5 * roughness], fluid))


    G.add_edge('Pump_outlet_1574', 'Wellhead_1574', obj=HE2_OilPipe([130.42], [2514.58], [0.125 * real_diam_coefficient], [roughness], fluid))
    G.add_edge('Zaboi_1574', 'Pump_intake_1574',obj=HE2_OilPipe([679.53], [-18.46], [0.143 * real_diam_coefficient], [5 * roughness], fluid))

    G.add_edge('Pump_outlet_1579', 'Wellhead_1579', obj=HE2_OilPipe([425.39], [2494.6], [0.125 * real_diam_coefficient], [roughness], fluid))
    G.add_edge('Zaboi_1579', 'Pump_intake_1579',obj=HE2_OilPipe([345.25], [281.16], [0.143 * real_diam_coefficient], [5 * roughness], fluid))
    
    G.add_edge('Pump_outlet_3116', 'Wellhead_3116', obj=HE2_OilPipe([167.57], [2453.43], [0.125 * real_diam_coefficient], [roughness], fluid))
    G.add_edge('Zaboi_3116', 'Pump_intake_3116',obj=HE2_OilPipe([715.74], [81.184], [0.143 * real_diam_coefficient], [5 * roughness], fluid))
    
    G.add_edge('Wellhead_3113', 'PAD_57', obj=HE2_OilPipe([100], [0], [0.125 * real_diam_coefficient], [roughness], fluid))
    G.add_edge('Wellhead_3118', 'PAD_57', obj=HE2_OilPipe([100], [0], [0.125 * real_diam_coefficient], [roughness], fluid))
    G.add_edge('Wellhead_3112', 'PAD_57', obj=HE2_OilPipe([100], [0], [0.125 * real_diam_coefficient], [roughness], fluid))
    G.add_edge('Wellhead_4235', 'PAD_57', obj=HE2_OilPipe([100], [0], [0.125 * real_diam_coefficient], [roughness], fluid))
    G.add_edge('Wellhead_3117', 'PAD_57', obj=HE2_OilPipe([100], [0], [0.125 * real_diam_coefficient], [roughness], fluid))
    G.add_edge('Wellhead_1493', 'PAD_57', obj=HE2_OilPipe([100], [0], [0.125 * real_diam_coefficient], [roughness], fluid))
    G.add_edge('Wellhead_1574', 'PAD_57', obj=HE2_OilPipe([100], [0], [0.125 * real_diam_coefficient], [roughness], fluid))
    G.add_edge('Wellhead_1579', 'PAD_57', obj=HE2_OilPipe([100], [0], [0.125 * real_diam_coefficient], [roughness], fluid))
    G.add_edge('Wellhead_3116', 'PAD_57', obj=HE2_OilPipe([100], [0], [0.125 * real_diam_coefficient], [roughness], fluid))
    
    
    #Трубопроводная система
    G.add_edge('PAD_5', 'intake_pad_5', obj=HE2_OilPipe([10], [-0.2], [0.143 * real_diam_coefficient], [roughness], fluid))
    G.add_edge('intake_pad_5', 'intake_pad_9_36', obj=HE2_OilPipe([785], [-0.6], [0.199 * real_diam_coefficient], [roughness], fluid))
    G.add_edge('intake_pad_9_36', 'intake_pad_35', obj=HE2_OilPipe([2326], [4.6], [0.199 * real_diam_coefficient], [roughness], fluid))
    G.add_edge('intake_pad_35', 'UDR_1',obj=HE2_OilPipe([2648], [13.6], [0.199 * real_diam_coefficient], [roughness], fluid))
    G.add_edge('UDR_1', 'DNS_2',obj=HE2_OilPipe([200], [0.2], [0.305 * real_diam_coefficient], [roughness], fluid))
    G.add_edge('PAD_33', 'intake_pad_33',obj=HE2_OilPipe([440], [-1.9], [0.139 * real_diam_coefficient], [roughness], fluid))
    G.add_edge('intake_pad_33', 'UDR_2', obj=HE2_OilPipe([4394], [-4.6], [0.199 * real_diam_coefficient], [roughness], fluid))
    G.add_edge('UDR_2', 'UDR_1', obj=HE2_OilPipe([1087], [3.4], [0.253 * real_diam_coefficient], [roughness], fluid))
    G.add_edge('PAD_34', 'intake_pad_134', obj=HE2_OilPipe([818], [0.6], [0.143 * real_diam_coefficient], [roughness], fluid))
    G.add_edge('intake_pad_134', 'UDR_2', obj=HE2_OilPipe([3344], [10.2], [0.203 * real_diam_coefficient], [roughness], fluid))
    G.add_edge('PAD_39', 'intake_pad_59', obj=HE2_OilPipe([7568], [1.6], [0.139 * real_diam_coefficient], [roughness], fluid))
    G.add_edge('intake_pad_59', 'intake_pad_59_46', obj=HE2_OilPipe([4601], [4.4], [0.143 * real_diam_coefficient], [roughness], fluid))
    G.add_edge('intake_pad_59_46', 'intake_pad_47',obj=HE2_OilPipe([2625], [18.4], [0.139 * real_diam_coefficient], [roughness], fluid))
    G.add_edge('intake_pad_47', 'intake_pad_53',obj=HE2_OilPipe([2250], [-1.3], [0.203 * real_diam_coefficient], [roughness], fluid))
    G.add_edge('intake_pad_53', 'intake_pad_41',obj=HE2_OilPipe([1391], [9], [0.203 * real_diam_coefficient], [roughness], fluid))
    G.add_edge('intake_pad_41', 'ZKL_98',obj=HE2_OilPipe([9290], [-13.5], [0.257 * real_diam_coefficient], [roughness], fluid))
    G.add_edge('ZKL_98', 'intake_pad_49',obj=HE2_OilPipe([600], [3], [0.257 * real_diam_coefficient], [roughness], fluid))
    G.add_edge('intake_pad_49', 'intake_node_33', obj=HE2_OilPipe([4001], [0], [0.410 * real_diam_coefficient], [roughness], fluid))
    G.add_edge('intake_node_33', 'intake_pad_46', obj=HE2_OilPipe([2920], [2.1], [0.410 * real_diam_coefficient], [roughness], fluid))
    G.add_edge('intake_pad_46', 'node_15',obj=HE2_OilPipe([1181], [-5.9], [0.410 * real_diam_coefficient], [roughness], fluid))
    G.add_edge('node_15', 'UDR_2',obj=HE2_OilPipe([92], [0.2], [0.199 * real_diam_coefficient], [roughness], fluid))
    G.add_edge('PAD_49', 'intake_node_49', obj=HE2_OilPipe([492], [-0.8], [0.139 * real_diam_coefficient], [roughness], fluid))
    G.add_edge('intake_node_49', 'intake_pad_57', obj=HE2_OilPipe([3737], [-5.7], [0.143 * real_diam_coefficient], [roughness], fluid))
    G.add_edge('intake_pad_57', 'intake_pad_49',obj=HE2_OilPipe([3852], [-4.3], [0.309 * real_diam_coefficient], [roughness], fluid))
    G.add_edge('PAD_57', 'intake_pad_57',obj=HE2_OilPipe([370], [-1], [0.203 * real_diam_coefficient], [roughness], fluid))

    return G, inlets, juncs, outlets


def model_DNS_2(pressures:dict = {}, plasts:dict = {}, DNS_daily_debit = 0, pumps = None, pump_curves = None,
                fluid = None, roughness = 0.00001, real_diam_coefficient = 1, DNS_pressure = 4.8):

    G, inlets, juncs, outlets = build_DNS2_graph(pressures, plasts, pumps, pump_curves, fluid, roughness,
                                                 real_diam_coefficient, DNS_daily_debit, DNS_pressure=4.8)
    # inlets_Q = gimme_DNS2_inlets_outlets_Q()
    #Создаем солвер и решаем полученный расчетный граф
    solver = HE2_Solver(G)
    # solver.set_known_Q(inlets_Q)
    solver.solve(threshold=0.25)
    return G, inlets, juncs, outlets

def cut_single_well_subgraph(G, pad_name, well):
    rez = nx.DiGraph()  # Di = directed
    plast = f'PAD_{pad_name}_well_{well}'
    zaboi = f'Zaboi_{well}'
    intake = f'Pump_intake_{well}'
    outlet = f'Pump_outlet_{well}'
    wellhead = f'Wellhead_{well}'
    pad = f'PAD_{pad_name}'
    # rez.add_nodes_from([G[plast], G[zaboi], G[intake], G[outlet], G[wellhead]])
    nodelist = [plast, zaboi, intake, outlet, wellhead, pad]
    edgelist = [(plast, zaboi), (zaboi, intake), (intake, outlet), (outlet, wellhead), (wellhead, pad)]
    rez.add_nodes_from(nodelist)
    rez.add_edges_from(edgelist)
    node_objs = {n:G.nodes[n]['obj'] for n in nodelist[:-1]}
    node_objs[pad] = vrtxs.HE2_Boundary_Vertex('P', 5)
    edge_objs = {}
    for u, v in edgelist:
        obj = G[u][v]['obj']
        edge_objs[(u, v)] = obj
    nx.set_node_attributes(rez, name='obj', values=node_objs)
    nx.set_edge_attributes(rez, name='obj', values=edge_objs)
    return rez, nodelist


# Review:
# def foo(bar = {}):
# Это плохой паттерн. Если аргумент функции по умолчанию имеет изменяемое значение (словарь, список), то вызовы функции могут приводить к затейливым спецэффектам
# Содержимое аргумента при выполнении тела функции, может быть разным, при вызове с одной и той же строкой параметров. Начинает зависеть от того с какими аргументами вызывалась функция раньше.
# Вот и PyCharm это подчеркивает
def model_DNS_2_by_parts(pressures: dict = {}, plasts: dict = {}, pumps=None, pump_curves=None, fluid=None,
                         roughness=0.00001, real_diam_coefficient=1, well_list=None, DNS_daily_debit=0):

    G, inlets, juncs, outlets = build_DNS2_graph(pressures, plasts, pumps, pump_curves, fluid, roughness,
                                                 real_diam_coefficient, DNS_daily_debit)
    wells = inlets
    if well_list:
        wells = list(set(well_list) & set(inlets))

    good_cnt, bad_cnt = 0, 0
    for well in wells:
        l = well.split('_')
        subG, well_node_list = cut_single_well_subgraph(G, l[1], l[3])
        solver = HE2_Solver(subG)
        solver.solve()
        op_result = solver.op_result
        Q = subG.nodes[well]['obj'].result['Q']
        if op_result.fun < 1e-3 and Q > 0:
            # Q > 0 means node is a source
            # print(well, ' is ok')
            good_cnt+=1
            continue

        bad_cnt +=1
        if op_result.fun > 1e-3:
            print(f'NOT SOLVED, {op_result.fun: .3f}')
        print(well, subG.nodes[well]['obj'].result)
        u, v = list(G.edges(well))[0]
        print(f'{u}-->{v}', subG[u][v]['obj'].result, '\n')

    print(good_cnt, bad_cnt, len(wells))
    return None


def model_DNS_3(daily_debit_55=0, pressure_88=0, fluid=fluid, roughness=3.5, real_diam_coefficient=0.85,
                DNS_daily_debit=0):
    inlets = dict(PAD_88=vrtxs.HE2_Source_Vertex('P', pressure_88, fluid, 20),
                  PAD_55=vrtxs.HE2_Source_Vertex('Q', daily_debit_55 * fluid.SepOilDensity * 1000 / 86400, fluid = fluid, T = 20))
                  #PAD_55=vrtxs.HE2_Source_Vertex('P', pressure_55, fluid, 20))

    juncs = dict(intake_pad_40=vrtxs.HE2_ABC_GraphVertex(),
                 intake_pad_45=vrtxs.HE2_ABC_GraphVertex(),
                 intake_pad_43=vrtxs.HE2_ABC_GraphVertex(),
                 intake_pad_44=vrtxs.HE2_ABC_GraphVertex())

    outlets = dict(DNS_3=vrtxs.HE2_Boundary_Vertex('Q', DNS_daily_debit * fluid.SepOilDensity * 1000 / 86400))

    G = nx.DiGraph()  # Di = directed
    for k, v in {**inlets, **outlets, **juncs}.items():
        G.add_node(k, obj=v)


    G.add_edge('PAD_55', 'intake_pad_40', obj=HE2_OilPipe([4871], [2.4], [0.143 * real_diam_coefficient], [roughness], fluid))
    G.add_edge('intake_pad_40', 'intake_pad_45', obj=HE2_OilPipe([3562], [5], [0.203 * real_diam_coefficient], [roughness], fluid))
    G.add_edge('PAD_88', 'intake_pad_45', obj=HE2_OilPipe([5234], [0], [0.143 * real_diam_coefficient], [roughness], fluid))
    G.add_edge('intake_pad_45', 'intake_pad_43',obj=HE2_OilPipe([3042], [-6.5], [0.203 * real_diam_coefficient], [roughness], fluid))
    G.add_edge('intake_pad_43', 'intake_pad_44',obj=HE2_OilPipe([2652], [2.1], [0.203 * real_diam_coefficient], [roughness], fluid))
    G.add_edge('intake_pad_44', 'DNS_3',obj=HE2_OilPipe([15], [-1.5], [0.203 * real_diam_coefficient], [roughness], fluid))
    solver = HE2_Solver(G)
    solver.solve()
    return G

def print_wells_pressures(G, wells):
    colorama.init()
    table_header = '                                            bottom   intake   outlet  wellhead'
    print(table_header)
    for pad_well in wells:
        l = pad_well.split('_')
        if len(l) < 4:
            continue
        pad, well = l[1], l[3]
        well_subgraph, well_nodes = cut_single_well_subgraph(G, pad, well)
        row_header = 'pad ' + f' {pad}'[-2:] + ' well ' + f'  {well}'[-4:] + ', from plast:   '
        print(row_header, end=' ')
        for n in well_nodes:
            P = G.nodes[n]['obj'].result['P_bar']
            prefix = Back.RED if P <= 1 else Style.RESET_ALL
            print(prefix + f'{P:8.3f}', end= ' ')
        print(Style.RESET_ALL + '   up to pad collector')

def print_solution(G):
    colorama.init()
    table_header = f' {"start":>20} {"P_bar":>7} {"Q kg/s":>8}   {"end":>20} {"P_bar":>7} {"Q kg/s":>8}    {"X kg/s":>7}'
    print(table_header)
    for e in G.edges:
        u, v = e
        obj = G[u][v]['obj']
        u_obj = G.nodes[u]['obj']
        v_obj = G.nodes[v]['obj']
        x = obj.result['x']
        p_u = u_obj.result['P_bar']
        pu_str = f'{p_u:8.3f}'
        if p_u <= 1:
            pu_str = Back.RED + pu_str + Style.RESET_ALL

        p_v = v_obj.result['P_bar']
        pv_str = f'{p_v:8.3f}'
        if p_v <= 1:
            pv_str = Back.RED + pv_str + Style.RESET_ALL

        q_u = u_obj.result['Q']
        q_u_str = ''
        if abs(q_u) > 1e-5:
            q_u_str = f"{q_u:8.3f}"

        q_v = v_obj.result['Q']
        q_v_str = ''
        if abs(q_v) > 1e-5:
            q_v_str = f"{q_v:8.3f}"

        row = f' {u:>20} {pu_str:>8} {q_u_str:>8}   {v:>20} {pv_str:>8} {q_v_str:>8}    {x:7.3f}{Style.RESET_ALL}'
        print(row)






