from Solver.HE2_Solver import HE2_Solver
from GraphNodes import HE2_Vertices as vrtxs

def remake_graph(G, intake_nodes):
    solver = HE2_Solver(G)
    solver.solve()
    p_dict = {n: G.nodes[n]['obj']['result'].P for n in intake_nodes}
    for n in intake_nodes:
        fluid = G.nodes[n]['obj'].fluid
        G.nodes[n]['obj'] = vrtxs.HE2_Source_Vertex('P', p_dict[n], fluid, 20)
    return G

