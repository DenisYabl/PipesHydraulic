import HE2_ABC as abc

class HE2_MockEdge(abc.HE2_ABC_GraphEdge):
    def __init__(self, P):
        self.P = P

    def perform_calc(self, P_bar, T_C, X_kgsec, unifloc_direction=-1):
        return self.P, T_C

