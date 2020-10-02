import HE2_ABC as abc
from HE2_Fluid import HE2_DummyWater
import uniflocpy.uTools.uconst as uc
import numpy as np


class HE2_Pipe():
    pass

class HE2_WaterPipeSegment():
    def __init__(self, fluid=None, inner_diam_m=None, roughness_m=None, L_m=None, downhill_m=None):
        if fluid is None:
            fluid = HE2_DummyWater()
        self.fluid = fluid
        self.inner_diam_m = inner_diam_m
        self.roughness_m = roughness_m
        self.L_m = L_m
        self.downhill_m = downhill_m

    def decode_direction(self, flow, unifloc_direction):
        assert unifloc_direction == -1, 'not impl!'
        grav_sign = -1
        fric_sign = np.sign(flow)
        t_sign = 1
        return grav_sign, fric_sign, t_sign

    def calc_P_friction_gradient_Pam(self, P_bar, T_C, X_kgsec):
        assert X_kgsec >= 0
        if X_kgsec == 0:
            return 0
        # Fluid.calc will be optimized at lower level. So we will call it every time
        self.fluid.calc(P_bar, T_C)
        Rho_kgm3 = self.fluid.rho_wat_kgm3
        mu_pasec = uc.cP2pasec(self.fluid.mu_wat_cp) # dynamic viscocity
        Q_m3sec = X_kgsec / Rho_kgm3
        D_m = self.inner_diam_m
        Area_m2 = uc.pi*D_m**2/4
        V_msec = Q_m3sec / Area_m2
        Re = Rho_kgm3 * V_msec * D_m / mu_pasec
        k_m = self.roughness_m
        if Re < 2300:
            lambda_fr = 68/ Re
        else:
            lambda_fr = 0.11 * (k_m/D_m + 68.5/Re) ** 0.25
        P_fric_grad_Pam = 0.5 * lambda_fr * V_msec**2 * Rho_kgm3 / D_m
        return P_fric_grad_Pam

    def calc_T_gradient_Cm(self, P_bar, T_C, X_kgsec):
        return 0

    def perform_calc(self, P_bar, T_C, X_kgsec, unifloc_direction=-1):
        P_fric_grad_Pam = self.calc_P_friction_gradient_Pam(P_bar, T_C, abs(X_kgsec))
        dP_fric_Pa = P_fric_grad_Pam * self.L_m
        Rho_kgm3 = self.fluid.rho_wat_kgm3
        dP_gravity_Pa = Rho_kgm3 * uc.g * self.downhill_m
        grav_sign, fric_sign, t_sign = self.decode_direction(X_kgsec, unifloc_direction)
        P_drop_bar = uc.Pa2bar(grav_sign * dP_gravity_Pa + fric_sign * dP_fric_Pa)
        P_rez_bar = P_bar - P_drop_bar
        T_grad_Cm = self.calc_T_gradient_Cm(P_bar, T_C, X_kgsec)
        T_rez_C = T_C - t_sign * T_grad_Cm * self.L_m
        return P_rez_bar, T_rez_C

