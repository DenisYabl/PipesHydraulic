from Hydraulics.Properties.Mishenko import Mishenko, from_oil_params
from Tools.HE2_ABC import HE2_ABC_Fluid, oil_params
from typing import List, Tuple
import numpy as np
from Tools.HE2_Logger import getLogger
logger = getLogger(__name__)

class HE2_DummyWater(HE2_ABC_Fluid):
    def __init__(self):
        self.rho_wat_kgm3 = 1000
        self.mu_wat_cp = 1

    def calc(self, P_bar, T_C):
        pass


class HE2_BlackOil(HE2_ABC_Fluid):
    def __init__(self, oil_params: oil_params):
        # check_for_nan(**oil_params)
        self.oil_params = oil_params
        self.CurrentLiquidDensity_kg_m3 = (self.oil_params.oildensity_kg_m3 * (1 - self.oil_params.volumewater_percent / 100) +
                                                      self.oil_params.waterdensity_kg_m3 * self.oil_params.volumewater_percent / 100)


    def calc(self, P_bar, T_C, X_kgsec, IntDiameter=None):
        P_for_PVT = max(abs(P_bar), 0.75)
        calc_params = self.oil_params

        tubing = {"IntDiameter": IntDiameter} if IntDiameter else None
        temp_mishenko = from_oil_params(P_for_PVT, T_C, X_kgsec, calc_params=calc_params, tubing=tubing)
        #Side effects
        self.CurrentLiquidDensity_kg_m3 = temp_mishenko.CurrentLiquidDensity_kg_m3
        self.CurrentOilViscosity_Pa_s = temp_mishenko.CurrentOilViscosity_Pa_s
        #Return for pressure gradient calculation
        return temp_mishenko


# fieldlist = ['sat_P_bar', 'plastT_C', 'gasFactor', 'oildensity_kg_m3', 'waterdensity_kg_m3', 'gasdensity_kg_m3',
#             'oilviscosity_Pa_s', 'volumewater_percent', 'volumeoilcoeff']
# oil_params = namedtuple('oil_params', fieldlist)

def check_all_are_the_same(arr, msg):
    max_arr = max(arr)
    if min(arr) != max_arr:
        logger.error(msg)
        raise NotImplementedError
    return max_arr

def make_fluid_vectors(fluids):
    ops = [fl.oil_params for fl in fluids]
    oil_ro_vec = np.array([op.oildensity_kg_m3 for op in ops])
    wat_ro_vec = np.array([op.waterdensity_kg_m3 for op in ops])
    gas_ro_vec = np.array([op.gasdensity_kg_m3 for op in ops])
    gf_vec = np.array([op.gasFactor for op in ops])
    wc_vec = np.array([op.volumewater_percent for op in ops])*0.01
    return oil_ro_vec, wat_ro_vec, gas_ro_vec, gf_vec, wc_vec


def dot_product(xs_vec, fluids, fluid_vectors=None) -> HE2_BlackOil:

    '''
    :return: new fluid instance, dot product Xs and fluids
    '''
    # if (Xs_and_fluids is None) or (len(Xs_and_fluids) == 0):
    #     logger.error('Empty fluids list')
    #     raise ValueError

    # xs_vec = np.array([x for x, fl in Xs_and_fluids])
    # ops = [fl.oil_params for x, fl in Xs_and_fluids]

    if np.sum(xs_vec != 0) == 1:
        op = fluids[np.argmax(xs_vec)].oil_params
        return HE2_BlackOil(op)

    if fluid_vectors is None:
        fluid_vectors = make_fluid_vectors(fluids)
    oil_ro_vec, wat_ro_vec, gas_ro_vec, gf_vec, wc_vec = fluid_vectors


    # sat_P_vec = np.array([op.sat_P_bar for op in ops])
    # sat_P = check_all_are_the_same(sat_P_vec, 'dot product for fluid.saturation_pressure is not implemented')
    #
    # plast_T_vec = np.array([op.plastT_C for op in ops])
    # plast_T= check_all_are_the_same(plast_T_vec, 'dot product for fluid.plast_temperature is not implemented')
    #
    # oil_Visc_vec = np.array([op.oilviscosity_Pa_s for op in ops])
    # oil_Visc = check_all_are_the_same(oil_Visc_vec, 'dot product for fluid.oil_viscosity is not implemented')
    #
    # Volume_keff_vec = np.array([op.volumeoilcoeff for op in ops])
    # Volume_keff = check_all_are_the_same(Volume_keff_vec, 'dot product for fluid.volumeoilcoeff is not implemented')


    owg_mix_pseudo_density_vec = oil_ro_vec * (1 - wc_vec) + wat_ro_vec * wc_vec + gas_ro_vec * (1 - wc_vec) * gf_vec
    Q_owg_vec = xs_vec / owg_mix_pseudo_density_vec
    Qo_vec = (1 - wc_vec) * Q_owg_vec
    Qw_vec = wc_vec * Q_owg_vec
    Qg_vec = (1 - wc_vec) * gf_vec* Q_owg_vec
    Xo_vec = Qo_vec * oil_ro_vec
    Xw_vec = Qw_vec * wat_ro_vec
    Xg_vec = Qg_vec * gas_ro_vec
    # np.testing.assert_almost_equal(Xo_vec + Xw_vec + Xg_vec, xs_vec)

    Xo = np.sum(Xo_vec)
    Xw = np.sum(Xw_vec)
    Xg = np.sum(Xg_vec)
    Qo = np.sum(Qo_vec)
    Qw = np.sum(Qw_vec)
    Qg = np.sum(Qg_vec)
    # oil_ro = np.round(Xo / Qo, 5)
    # wat_ro = np.round(Xw / Qw, 5)
    # gas_ro = np.round(Xg / Qg, 5)
    # wc = np.round(Qw / (Qo + Qw), 5)
    # gf = np.round(Qg / Qo, 5)
    oil_ro = Xo / Qo
    wat_ro = Xw / Qw
    gas_ro = Xg / Qg
    wc = Qw / (Qo + Qw)
    gf = Qg / Qo

    # norm_keff = 1/sum(xs_vec)
    # oil_ro2 = norm_keff * np.dot(xs_vec, oil_ro_vec)
    # wat_ro2 = norm_keff * np.dot(xs_vec, wat_ro_vec)
    # gas_ro2 = norm_keff * np.dot(xs_vec, gas_ro_vec)
    # np.testing.assert_almost_equal(oil_ro, oil_ro2)
    # np.testing.assert_almost_equal(wat_ro, wat_ro2)
    # np.testing.assert_almost_equal(gas_ro, gas_ro2)

    op0 = fluids[0].oil_params
    sat_P = op0.sat_P_bar
    plast_T = op0.plastT_C
    oil_Visc = op0.oilviscosity_Pa_s
    Volume_keff = op0.volumeoilcoeff
    rez_oil_params = oil_params(sat_P, plast_T, gf, oil_ro, wat_ro, gas_ro, oil_Visc, wc*100, Volume_keff)
    rez = HE2_BlackOil(rez_oil_params)

    # if np.sum(xs_vec != 0) == 1:
    #     op1 = fluids[np.argmax(xs_vec)].oil_params
    #     op2 = rez_oil_params
    #     diff = np.array(op1) - np.array(op2)
    #     if np.linalg.norm(diff) > 1e-8:
    #         logger.warning('Something wrong with fluids')

    return rez


def gimme_dummy_oil_params(volumeWater=50):
    rez = oil_params(sat_P_bar = 66.7, plastT_C = 84, gasFactor = 39, oildensity_kg_m3 = 826,
                     waterdensity_kg_m3 = 1015, gasdensity_kg_m3 = 1, oilviscosity_Pa_s = 35e-3, volumewater_percent = volumeWater, volumeoilcoeff = 1.015)
    return rez

def gimme_dummy_BlackOil(VolumeWater=50):
    oil_params = gimme_dummy_oil_params(volumeWater=VolumeWater)
    rez = HE2_BlackOil(oil_params)
    return rez


if __name__ == '__main__':
    fl = HE2_DummyWater()
    fl.calc(100, 100)
    print(fl.rho_wat_kgm3)