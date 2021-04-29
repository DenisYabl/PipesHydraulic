import math

from Tools import HE2_ABC as abc
from Fluids.HE2_Fluid import gimme_dummy_BlackOil
import uniflocpy.uTools.uconst as uc
import numpy as np
import pandas as pd
from scipy.interpolate import interp1d
from Tools.HE2_Logger import check_for_nan, getLogger
logger = getLogger(__name__)

A_keff = 1
B_keff = 1.4

# TODO Нужно сузить класс HE2_WellPump. Не надо таскать в него датафрейм. Нужно создавать его от трех векторов Q, H, Eff
# И рядом сделать конструктор, который берет уже датафрейм, выделяет из него три вектора и создает объект HE2_WellPump

def create_HE2_WellPump_instance_from_dataframe(full_HPX:pd.DataFrame, model = "", fluid = None, frequency = 50):
    base_HPX = full_HPX[full_HPX["pumpModel"] == model].sort_values('debit')
    base_HPX = base_HPX.drop(base_HPX[base_HPX.eff == 0].index)
    p_vec = base_HPX["pressure"].values
    q_vec = base_HPX["debit"].values
    n_vec = base_HPX["power"].values
    eff_vec = base_HPX["eff"].values
    pump = HE2_WellPump(p_vec, q_vec, n_vec, eff_vec, model, fluid, frequency)
    return pump

# class HE2_WellPump(abc.HE2_ABC_GraphEdge):
#     def __init__(self, Q_P_N_table=None, fluid = HE2_DummyOil, frequency = 50):

class HE2_WellPump(abc.HE2_ABC_Pipeline, abc.HE2_ABC_GraphEdge):
    def __init__(self, p_vec, q_vec, n_vec, eff_vec, model = "", fluid = None, frequency = 50, engine_efficiency = 0.92):
        if fluid is None:
            fluid = gimme_dummy_BlackOil()
        self.engine_efficiency = engine_efficiency
        self.fluid = fluid
        self.model = model
        self.intermediate_results = []
        self.frequency = frequency
        self._printstr = self.model
        self.base_q = q_vec # self.base_* - НРХ насоса по воде
        self.base_p = p_vec
        self.base_n = n_vec
        self.base_eff = eff_vec
        self.stages_ratio = 1

        visc_approx = self.fluid.calc(30, 20, 0).CurrentOilViscosity_Pa_s
        # pseudo_vec = 1.95 * math.pow(visc_approx, 0.5) * 0.04739 * ((p_vec / 0.3048) ** 0.25739) * (((q_vec / 0.227) ** 0.5) ** 0.5)
        pseudo_vec = 1.95 * 0.04739 * visc_approx**0.5 * (p_vec / 0.3048)**0.25739 * (q_vec / 0.227)**0.25
        Cq_vec = np.polyval([0.9873, 0.009019, -0.0016233, 0.00007233, -0.0000020258, 0.000000021009][::-1], pseudo_vec)
        Cp_vec = np.polyval([1.0045, -0.002664, -0.00068292, 0.000049706, -0.0000016522, 0.000000019172][::-1], pseudo_vec)
        Ceff_vec = np.polyval([1.0522, -0.03512, -0.00090394, 0.00022218, -0.00001198, 0.00000019895][::-1], pseudo_vec)

        self.q_vec = self.base_q * Cq_vec
        self.p_vec_50hz = self.base_p * Cp_vec
        self.eff_vec = eff_vec * Ceff_vec
        self.n_vec_50hz = (self.q_vec * 9.81 * 1000 / 86400) * self.p_vec_50hz / (self.eff_vec/100)
        self.p_vec = self.p_vec_50hz * (self.frequency/50)**2 * self.stages_ratio
        self.n_vec = self.n_vec_50hz * (self.frequency/50)**2

        self.min_Q = self.q_vec.min()
        self.max_Q = self.q_vec.max()
        self.power = 0

        self.get_pressure_raise_1 = None
        self.get_pressure_raise_2 = None
        self.get_pressure_raise_3 = None
        self.make_extrapolators()


    def __str__(self):
        return self._printstr

    def perform_calc(self, P_bar, T_C, X_kgsec, unifloc_direction):
        assert unifloc_direction in [0, 1, 10, 11]
        calc_direction = 1 if unifloc_direction >= 10 else -1
        flow_direction = 1 if unifloc_direction % 10 == 1 else - 1



        if calc_direction == 1:
            return self.perform_calc_forward(P_bar, T_C, X_kgsec)
        else:
            return self.perform_calc_backward(P_bar, T_C, X_kgsec)

    def perform_calc_forward(self, P_bar, T_C, X_kgsec):
        p, t = P_bar, T_C
        fl = self.fluid.calc(P_bar, T_C, X_kgsec)
        p, t = self.calculate_pressure_differrence(p, t, X_kgsec, 1, fl)
        self.intermediate_results += [(p, t)]
        return p, t

    def perform_calc_backward(self, P_bar, T_C, X_kgsec):
        p, t = P_bar, T_C
        fl = self.fluid.calc(P_bar, T_C, X_kgsec)
        p, t = self.calculate_pressure_differrence(p, t, X_kgsec, -1, fl)
        return p, t

    def calculate_pressure_differrence(self, P_bar, T_C, X_kgsec, calc_direction, mishenko, unifloc_direction=-1):
        check_for_nan(P_bar=P_bar, T_C=T_C, X_kgsec=X_kgsec)
        #Определяем направления расчета
        liquid_debit = X_kgsec * 86400 / mishenko.CurrentLiquidDensity_kg_m3
        grav_sign, fric_sign, t_sign = self.decode_direction(X_kgsec, calc_direction, unifloc_direction)
        self.power = 0
        if liquid_debit <= 0:
            get_pressure_raise = self.get_pressure_raise_1
            get_eff = self.get_eff_1
        elif (self.min_Q < liquid_debit) and (abs(X_kgsec) * 86400 / mishenko.CurrentLiquidDensity_kg_m3 < self.max_Q) :
            get_pressure_raise = self.get_pressure_raise_2
            get_eff = self.get_eff_2
            self.power = self.n_interpolator(liquid_debit)
        else:
            get_pressure_raise = self.get_pressure_raise_3
            get_eff = self.get_eff_3

        pressure_raise = get_pressure_raise(liquid_debit) * 9.81 *  mishenko.CurrentLiquidDensity_kg_m3

        P_rez_bar = P_bar + calc_direction * uc.Pa2bar(pressure_raise)

        eff = get_eff(liquid_debit)
        #T_rez_C = T_C + calc_direction * self.calculate_temperature_raise(pressure_raise, eff, mishenko)
        T_rez_C = T_C
        check_for_nan(P_rez_bar=P_rez_bar, T_rez_C=T_rez_C)

        return P_rez_bar, T_rez_C

    def decode_direction(self, flow, calc_direction, unifloc_direction):
        '''
        :param unifloc_direction - направление расчета и потока относительно  координат.
            11 расчет и поток по координате
            10 расчет по координате, поток против
            00 расчет и поток против координаты
            00 расчет против координаты, поток по координате
            unifloc_direction перекрывает переданные flow, calc_direction
        '''
        flow_direction = np.sign(flow) if flow !=0 else 1
        if unifloc_direction in [0, 1, 10, 11]:
            calc_direction = 1 if unifloc_direction >= 10 else -1
            flow_direction = 1 if unifloc_direction % 10 == 1 else - 1

        assert calc_direction in [-1, 1]
        grav_sign = calc_direction
        fric_sign = flow_direction * calc_direction
        t_sign = calc_direction
        return grav_sign, fric_sign, t_sign

    def make_extrapolators(self):
        zero_head = self.p_vec[0]
        last_head = self.p_vec[-1]
        self.get_pressure_raise_1 = lambda x: zero_head + A_keff * abs(x) ** B_keff
        self.get_pressure_raise_2 = interp1d(self.q_vec, self.p_vec, kind="quadratic")
        self.get_pressure_raise_3 = lambda x: last_head - A_keff * (x - self.max_Q) ** B_keff

        self.get_eff_1 = lambda x: 5
        self.get_eff_2 = interp1d(self.q_vec, self.eff_vec, kind="quadratic")
        self.get_eff_3 = lambda x: 5

        self.n_interpolator = interp1d(self.q_vec, self.n_vec, kind="quadratic")

    def changeFrequency(self, new_frequency):
        self.frequency = new_frequency
        self.p_vec = self.p_vec_50hz * (self.frequency/50)**2 * self.stages_ratio
        self.n_vec = self.n_vec_50hz * (self.frequency/50)**2
        self.make_extrapolators()

    def change_stages_ratio(self, new_ratio):
        self.stages_ratio = new_ratio
        self.p_vec = self.p_vec_50hz * (self.frequency/50)**2 * self.stages_ratio
        self.n_vec = self.n_vec_50hz * (self.frequency/50)**2
        self.make_extrapolators()

    def calculate_temperature_raise(self, deltaP, eff, mishenko):
        deltaT = deltaP / (mishenko.CurrentLiquidDensity_kg_m3 * mishenko.Thermal_capacity) * (1 / (eff / 100 * self.engine_efficiency) - 1)
        return deltaT
