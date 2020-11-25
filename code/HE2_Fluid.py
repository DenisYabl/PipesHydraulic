from uniflocpy.uPVT.BlackOil_model import Fluid
from HE2_ABC import HE2_ABC_Fluid
import numpy as np

class HE2_DummyWater(HE2_ABC_Fluid):
    def __init__(self):
        self.rho_wat_kgm3 = 1000
        self.mu_wat_cp = 1

    def calc(self, P_bar, T_C):
        pass


class HE2_DummyFluid(HE2_ABC_Fluid):
    def __init__(self, rho_kgm3):
        self.rho_kgm3 = rho_kgm3
        self.mu_cp = 1

    def calc(self, P_bar, T_C):
        pass

def dot_product(weigths, fluids):
    '''
    :param weigths: weigths vector
    :param fluids: fluids vector
    :return: new fluid instance, dot product weigths and fluids
    '''
    assert len(weigths) == len(fluids)
    wghts_norm = weigths / sum(weigths)
    rho_list = []
    for f in fluids:
        rho_list += [f.rho_kgm3]
    rhos = np.array(rho_list)
    rez_rho = np.dot(wghts_norm, rhos)
    rez = HE2_DummyFluid(rez_rho)
    return rez

if __name__ == '__main__':
    fl = HE2_DummyWater()
    fl.calc(100, 100)
    print(fl.rho_wat_kgm3)