from Fluids.HE2_Fluid import HE2_BlackOil
from GraphEdges.HE2_Pipe import HE2_OilPipe
from Tools.HE2_ABC import oil_params
import matplotlib.pyplot as plt
import numpy as np
# Tailaki_oil_params = oil_params(sat_P_bar=67, plastT_C=84, gasFactor=36, oildensity_kg_m3=826,
#         waterdensity_kg_m3=1015, gasdensity_kg_m3=1, oilviscosity_Pa_s=35e-3, volumewater_percent=50, volumeoilcoeff=1.017)

op = oil_params(sat_P_bar=67, plastT_C=84, gasFactor=36, oildensity_kg_m3=826,
        waterdensity_kg_m3=1015.0, gasdensity_kg_m3=1.0, oilviscosity_Pa_s=0.035, volumewater_percent=26.830144218677827, volumeoilcoeff=1.017)

rez = []
start = 0.1
stop = 10
for X_kg_sec in np.arange(start, stop, 0.001):
        fluid = HE2_BlackOil(op)
        pipe = HE2_OilPipe([9290], [13.5], [0.2056], [1e-5], fluid)
        P_bar = 8.325713134838217
        T = 20

        P1, T1 = pipe.perform_calc_forward(P_bar, T, X_kg_sec)
        rez.append(P1)
        #print(rez)

plt.plot(np.arange(start, stop, 0.001), rez)
plt.show()
#print(P1, P2)
pass