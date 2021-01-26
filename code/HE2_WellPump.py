import math

import HE2_ABC as abc
from HE2_Fluid import HE2_DummyWater, HE2_OilWater, HE2_DummyOil
from functools import reduce
import uniflocpy.uTools.uconst as uc
import numpy as np
import Hydraulics.Methodics.Mukherjee_Brill as mb
import pandas as pd
from scipy.interpolate import interp1d

from Hydraulics.Formulas import RecalculateNRH_coefficients


class HE2_WellPump(abc.HE2_ABC_Pipeline, abc.HE2_ABC_GraphEdge):
    def __init__(self, full_HPX:pd.DataFrame, model = "", fluid = HE2_DummyOil, IntDiameter = 0.12, frequency = 50):
        self.base_HPX = full_HPX[full_HPX["pumpModel"] == model]
        self.fluid = fluid
        self.model = model
        self.intermediate_results = []
        self.IntDiameter = IntDiameter
        self.frequency = frequency
        self._printstr = self.base_HPX.to_string()

        self.true_HPX = self.base_HPX.copy()
        self.true_HPX["pseudo"] = 1.95 * math.pow(fluid.SepOilDynamicViscosity, 0.5) * 0.04739 * ((self.true_HPX["pressure"] / 0.3048) ** 0.25739) * (((self.true_HPX["debit"] / 0.227) ** 0.5) ** 0.5)
        self.true_HPX["Cq"] = 0.9873 * (self.true_HPX["pseudo"] ** 0) + 0.009019 * (self.true_HPX["pseudo"] ** 1) - 0.0016233 * (self.true_HPX["pseudo"] ** 2) + 0.00007233 * (self.true_HPX["pseudo"] ** 3) - 0.0000020258 * (
        self.true_HPX["pseudo"] ** 4) + 0.000000021009 * (self.true_HPX["pseudo"] ** 5)
        self.true_HPX["Ch"] = 1.0045 * (self.true_HPX["pseudo"] ** 0) - 0.002664 * (self.true_HPX["pseudo"] ** 1) - 0.00068292 * (self.true_HPX["pseudo"] ** 2) + 0.000049706 * (self.true_HPX["pseudo"] ** 3) - 0.0000016522 * (
        self.true_HPX["pseudo"] ** 4) + 0.000000019172 * (self.true_HPX["pseudo"] ** 5)
        self.true_HPX["Ceff"] = 1.0522 * (self.true_HPX["pseudo"] ** 0) - 0.03512 * (self.true_HPX["pseudo"] ** 1) - 0.00090394 * (self.true_HPX["pseudo"] ** 2) + 0.00022218 * (self.true_HPX["pseudo"] ** 3) - 0.00001198 * (
        self.true_HPX["pseudo"] ** 4) + 0.00000019895 * (self.true_HPX["pseudo"] ** 5)

        self.true_HPX["debit"] = self.true_HPX["debit"] * self.true_HPX["Cq"]
        self.true_HPX["pressure"] = self.true_HPX["pressure"] * self.true_HPX["Ch"]
        self.true_HPX["eff"] = self.true_HPX["eff"] * self.true_HPX["Ceff"]
        self.true_HPX["pressure"] = self.true_HPX["pressure"] * (self.frequency / 50) ** 2
        self.true_HPX["power"] = (self.true_HPX["debit"] * self.true_HPX["pressure"] * 9.81 * 1000) / (3960 * self.true_HPX["eff"])

    def __str__(self):
        return self._printstr

    def perform_calc(self, P_bar, T_C, X_kgsec, unifloc_direction):
        assert unifloc_direction in [0, 1, 10, 11]
        calc_direction = 1 if unifloc_direction >= 10 else -1
        flow_direction = 1 if unifloc_direction % 10 == 1 else - 1



        if calc_direction == 1:
            return self.perform_calc_forward(P_bar, T_C, flow_direction * abs(X_kgsec))
        else:
            return self.perform_calc_backward(P_bar, T_C, flow_direction * abs(X_kgsec))

    def perform_calc_forward(self, P_bar, T_C, X_kgsec):
        p, t = P_bar, T_C
        p, t = self.calculate_pressure_differrence(p, t, X_kgsec, 1, self.fluid.calc(P_bar, T_C, X_kgsec, 0.12))
        self.intermediate_results += [(p, t)]
        return p, t

    def perform_calc_backward(self, P_bar, T_C, X_kgsec):
        p, t = P_bar, T_C




        p, t = self.calculate_pressure_differrence(p, t, X_kgsec, -1, self.fluid.calc(P_bar, T_C, X_kgsec, self.IntDiameter))

        self.intermediate_results += [(p, t)]
        return p, t

    def calculate_pressure_differrence(self, P_bar, T_C, X_kgsec, calc_direction, mishenko, unifloc_direction=-1):
        #Определяем направления расчета
        fric_sign, t_sign = self.decode_direction(X_kgsec, calc_direction, unifloc_direction)
        if self.true_HPX["debit"].min() <= abs(X_kgsec) * 86400 / mishenko.CurrentLiquidDensity <= self.true_HPX["debit"].max():
            get_pressure_raise = interp1d(self.true_HPX["debit"], self.true_HPX["pressure"])
        else:
            get_pressure_raise = interp1d(self.true_HPX["debit"], self.true_HPX["pressure"], fill_value="extrapolate")

        if fric_sign > 0:
            P_rez_bar = P_bar + uc.Pa2bar(get_pressure_raise(abs(X_kgsec) * 86400 / mishenko.CurrentLiquidDensity) * 9.81 *  mishenko.CurrentLiquidDensity)
            T_rez_C = T_C
        else:
            P_rez_bar = P_bar - uc.Pa2bar(get_pressure_raise(abs(X_kgsec) * 86400 / mishenko.CurrentLiquidDensity) * 9.81 *  mishenko.CurrentLiquidDensity)
            T_rez_C = T_C
        return P_rez_bar, T_rez_C

    def decode_direction(self, flow, calc_direction, unifloc_direction):
        '''
        :param unifloc_direction - направление расчета и потока относительно  координат.
            11 расчет и поток по координате
            10 расчет по координате, поток против
            00 расчет и поток против координаты
            01 расчет против координаты, поток по координате
            unifloc_direction перекрывает переданные flow, calc_direction
            grav_sign не нужен, поскольку он учитывается в Mukherjee_Brill
        '''
        flow_direction = np.sign(flow)
        if unifloc_direction in [0, 1, 10, 11]:
            calc_direction = 1 if unifloc_direction >= 10 else -1
            flow_direction = 1 if unifloc_direction % 10 == 1 else - 1

        assert calc_direction in [-1, 1]
        fric_sign = flow_direction * calc_direction
        t_sign = calc_direction
        return fric_sign, t_sign


