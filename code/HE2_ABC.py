from abc import ABC
from abc import abstractmethod


class HE2_ABC_Fluid(ABC):
    @abstractmethod
    def calc(self, P_bar, T_C):
        pass

class HE2_ABC_GraphEdge(ABC):
    @abstractmethod
    def perform_calc(self, P_bar, T_C, X_kgsec, unifloc_direction=-1):
        pass

class HE2_ABC_GraphVertex(ABC):
    pass

class HE2_ABC_Graph(ABC):
    pass

class HE2_ABC_PipeSegment(ABC):
    pass

class HE2_ABC_Pipeline(ABC):
    pass

