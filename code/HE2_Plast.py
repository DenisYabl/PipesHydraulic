import HE2_ABC as abc
from HE2_Fluid import HE2_DummyWater, HE2_OilWater, HE2_DummyOil
from functools import reduce
import uniflocpy.uTools.uconst as uc
import numpy as np
import Hydraulics.Methodics.Mukherjee_Brill as mb


class HE2_Plast(abc.HE2_ABC_Pipeline, abc.HE2_ABC_GraphEdge):
    def __init__(self, productivity = 0, fluid = HE2_DummyOil):
        self.Productivity = productivity
        self.fluid = fluid
        self.intermediate_results = []
        self._printstr = ';\n '.join(f"Productivity coefficient: {self.Productivity}")


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
        p, t = self.calculate_pressure_differrence(p, t, X_kgsec, 1)
        self.intermediate_results += [(p, t)]
        return p, t

    def perform_calc_backward(self, P_bar, T_C, X_kgsec):
        p, t = P_bar, T_C
        p, t = self.calculate_pressure_differrence(p, t, X_kgsec, -1)
        self.intermediate_results += [(p, t)]
        return p, t

    def calculate_pressure_differrence(self, P_bar, T_C, X_kgsec, calc_direction, unifloc_direction=-1):
        #Определяем направления расчета
        fric_sign, t_sign = self.decode_direction(X_kgsec, calc_direction, unifloc_direction)
        if fric_sign > 0:
            P_rez_bar =  P_bar -  ((X_kgsec * 86400 /  self.fluid.calc(P_bar, T_C, X_kgsec, 1.5).CurrentLiquidDensity) / self.Productivity)
            T_rez_C = T_C
        else:
            fl =  self.fluid
            liq = fl.calc(P_bar, T_C, X_kgsec, 1.5)
            liq_dens = liq.CurrentLiquidDensity
            P_rez_bar = P_bar + ((X_kgsec * 86400 / liq_dens) / self.Productivity)
            T_rez_C = T_C
        result_pressure = P_rez_bar
        return result_pressure, T_rez_C

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