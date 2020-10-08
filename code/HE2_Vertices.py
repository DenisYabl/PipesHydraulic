from HE2_ABC import HE2_ABC_GraphVertex


class HE2_Boundary_Vertex(HE2_ABC_GraphVertex):
    def __init__(self, kind, value):
        self.kind = kind
        self.value = value
        self.is_source = False


class HE2_Source_Vertex(HE2_Boundary_Vertex):
    def __init__(self, kind, value, fluid, T):
        super().__init__(kind, value)
        self.fluid = fluid
        self.T = T
        self.is_source = True