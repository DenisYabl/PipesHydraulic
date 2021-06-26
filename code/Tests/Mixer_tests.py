import unittest
import networkx as nx
import Fluids.HE2_MixFluids as mixer
import Fluids.HE2_MixFluids2 as mixer2
from Tools.HE2_ABC import Root
import GraphNodes.HE2_Vertices as vrtxs
import numpy as np

class TestMixer(unittest.TestCase):
    def setUp(self):
        pass

    def test_1(self):
        G = nx.DiGraph()
        nx.add_path(G, [Root, 'a', 'b', 'c', 'd'])
        nx.add_path(G, ['c', 'e', 'f', 'g'])
        for n in G.nodes:
            G.nodes[n]['obj'] = vrtxs.HE2_ABC_GraphVertex()
        G2, mappings = mixer2.make_small_graph_for_mixer(G,removeRoot=True)
        self.assertEqual(len(G2.edges), 3)
        self.assertTrue(('a', 'c') in G2.edges)
        self.assertTrue(('c', 'd') in G2.edges)
        self.assertTrue(('c', 'g') in G2.edges)

    def test_2(self):
        G = nx.DiGraph()
        nx.add_path(G, [Root, 'a', 'b', 'c', 'd'])
        nx.add_path(G, [Root, 'e', 'f', 'g'])
        nx.add_path(G, ['c', 'h', 'i', 'g'])
        for n in G.nodes:
            G.nodes[n]['obj'] = vrtxs.HE2_ABC_GraphVertex()
        G2, mappings = mixer2.make_small_graph_for_mixer(G,removeRoot=True)
        rez = set(G2.edges)
        rez_have_to_be = {('a', 'c'), ('c', 'e'), ('c', 'd')}
        self.assertEqual(rez_have_to_be, rez)

    def test_3(self):
        G = nx.DiGraph()
        nx.add_path(G, [Root, 'a', 'b', 'c', 'd'])
        nx.add_path(G, [Root, 'e', 'f', 'g'])
        # nx.add_path(G, ['c', 'h', 'i', 'g'])
        G.add_edge('c', 'h')
        G.add_edge('i', 'h')
        G.add_edge('g', 'i')
        for n in G.nodes:
            G.nodes[n]['obj'] = vrtxs.HE2_ABC_GraphVertex()
        G2, mappings = mixer2.make_small_graph_for_mixer(G,removeRoot=True)
        rez = set(G2.edges)
        rez_have_to_be = {('a', 'c'), ('e', 'c'), ('c', 'd')}
        self.assertEqual(rez_have_to_be, rez)
        fluid_em = mappings['fluid_em']
        set1 = set(fluid_em[('a', 'c')])
        set1_have_to_be = {('a', 'b'),('b', 'c')}
        self.assertEqual(set1_have_to_be, set1)
        set2 = set(fluid_em[('e', 'c')])
        set2_have_to_be = {('e', 'f'), ('f', 'g'), ('g', 'i'), ('i', 'h'), ('c', 'h')}
        self.assertEqual(set2_have_to_be, set2)
        set3 = set(fluid_em[('c', 'd')])
        set3_have_to_be = {('c', 'd')}
        self.assertEqual(set3_have_to_be, set3)
        pos_em = mappings['x_pos_em']
        set1 = set(pos_em[('a', 'c')])
        set2 = set(pos_em[('e', 'c')])
        set3 = set(pos_em[('c', 'd')])
        set1_have_to_be = {('a', 'b')}
        set2_have_to_be = {('e', 'f')}
        set3_have_to_be = {('c', 'd')}
        self.assertEqual(set1_have_to_be, set1)
        self.assertEqual(set2_have_to_be, set2)
        self.assertEqual(set3_have_to_be, set3)
        neg_em = mappings['x_neg_em']
        set1 = set(neg_em[('a', 'c')])
        set2 = set(neg_em[('e', 'c')])
        set3 = set(neg_em[('c', 'd')])
        set1_have_to_be = set()
        set2_have_to_be = set()
        set3_have_to_be = set()
        self.assertEqual(set1_have_to_be, set1)
        self.assertEqual(set2_have_to_be, set2)
        self.assertEqual(set3_have_to_be, set3)


    def test_4(self):
        G = nx.DiGraph()
        nx.add_path(G, [Root, 'a', 'b', 'c', 'd'])
        nx.add_path(G, [Root, 'e', 'f', 'g'])
        # nx.add_path(G, ['c', 'h', 'i', 'g'])
        G.add_edge('c', 'h')
        G.add_edge('i', 'h')
        G.add_edge('g', 'i')
        G.add_edge('c', 'i')
        for n in G.nodes:
            G.nodes[n]['obj'] = vrtxs.HE2_ABC_GraphVertex()
        G2, mappings = mixer2.make_small_graph_for_mixer(G,removeRoot=True)
        rez = set(G2.edges)
        rez_have_to_be = {('a', 'c'), ('e', 'c'), ('c', 'd')}
        self.assertEqual(rez_have_to_be, rez)

        fluid_em = mappings['fluid_em']
        set1 = set(fluid_em[('a', 'c')])
        set1_have_to_be = {('a', 'b'),('b', 'c')}
        self.assertEqual(set1_have_to_be, set1)
        set2 = set(fluid_em[('e', 'c')])
        set2_have_to_be = {('e', 'f'), ('f', 'g'), ('g', 'i'), ('i', 'h'), ('c', 'h'), ('c', 'i')}
        self.assertEqual(set2_have_to_be, set2)
        set3 = set(fluid_em[('c', 'd')])
        set3_have_to_be = {('c', 'd')}
        self.assertEqual(set3_have_to_be, set3)

        pos_em = mappings['x_pos_em']
        set1 = set(pos_em[('a', 'c')])
        set2 = set(pos_em[('e', 'c')])
        set3 = set(pos_em[('c', 'd')])
        set1_have_to_be = {('a', 'b')}
        set2_have_to_be = {('e', 'f')}
        set3_have_to_be = {('c', 'd')}
        self.assertEqual(set1_have_to_be, set1)
        self.assertEqual(set2_have_to_be, set2)
        self.assertEqual(set3_have_to_be, set3)
        neg_em = mappings['x_neg_em']
        set1 = set(neg_em[('a', 'c')])
        set2 = set(neg_em[('e', 'c')])
        set3 = set(neg_em[('c', 'd')])
        set1_have_to_be = set()
        set2_have_to_be = set()
        set3_have_to_be = set()
        self.assertEqual(set1_have_to_be, set1)
        self.assertEqual(set2_have_to_be, set2)
        self.assertEqual(set3_have_to_be, set3)

    def test_5(self):
        G = nx.DiGraph()
        nx.add_path(G, [Root, 'a', 'b', 'c', 'd'])
        nx.add_path(G, [Root, 'e', 'f', 'g'])
        # nx.add_path(G, ['c', 'h', 'i', 'g'])
        G.add_edge('c', 'h')
        G.add_edge('i', 'h')
        G.add_edge('g', 'i')
        G.add_edge('i', 'c')
        for n in G.nodes:
            G.nodes[n]['obj'] = vrtxs.HE2_ABC_GraphVertex()
        G2, mappings = mixer2.make_small_graph_for_mixer(G,removeRoot=True)
        rez = set(G2.edges)
        rez_have_to_be = {('a', 'c'), ('e', 'c'), ('c', 'd')}
        self.assertEqual(rez_have_to_be, rez)
        fluid_em = mappings['fluid_em']
        set1 = set(fluid_em[('a', 'c')])
        set1_have_to_be = {('a', 'b'),('b', 'c')}
        self.assertEqual(set1_have_to_be, set1)
        set2 = set(fluid_em[('e', 'c')])
        set2_have_to_be = {('e', 'f'), ('f', 'g'), ('g', 'i'), ('i', 'h'), ('c', 'h'), ('i', 'c')}
        self.assertEqual(set2_have_to_be, set2)
        set3 = set(fluid_em[('c', 'd')])
        set3_have_to_be = {('c', 'd')}
        self.assertEqual(set3_have_to_be, set3)

        pos_em = mappings['x_pos_em']
        set1 = set(pos_em[('a', 'c')])
        set2 = set(pos_em[('e', 'c')])
        set3 = set(pos_em[('c', 'd')])
        set1_have_to_be = {('a', 'b')}
        set2_have_to_be = {('e', 'f')}
        set3_have_to_be = {('c', 'd')}
        self.assertEqual(set1_have_to_be, set1)
        self.assertEqual(set2_have_to_be, set2)
        self.assertEqual(set3_have_to_be, set3)
        neg_em = mappings['x_neg_em']
        set1 = set(neg_em[('a', 'c')])
        set2 = set(neg_em[('e', 'c')])
        set3 = set(neg_em[('c', 'd')])
        set1_have_to_be = set()
        set2_have_to_be = set()
        set3_have_to_be = set()
        self.assertEqual(set1_have_to_be, set1)
        self.assertEqual(set2_have_to_be, set2)
        self.assertEqual(set3_have_to_be, set3)


    def test_6(self):
        G = nx.DiGraph()
        edges = [('a', 'b'), ('c', 'b'), ('c', 'd'), ('f', 'g'), ('g', 'h'),
                 ('i', 'c'), ('c', 'i'), ('f', 'i'), ('f', 'i')]
        G.add_edges_from(edges)
        for n in G.nodes:
            G.nodes[n]['obj'] = vrtxs.HE2_ABC_GraphVertex()
        G2, mappings = mixer2.make_small_graph_for_mixer(G,removeRoot=True)
        rez = set(G2.edges)
        rez_have_to_be = {('a', 'c'), ('c', 'd'), ('c', 'h')}
        self.assertEqual(rez_have_to_be, rez)
        fluid_em = mappings['fluid_em']
        set1 = set(fluid_em[('a', 'c')])
        set1_have_to_be = {('a', 'b'),('c', 'b')}
        self.assertEqual(set1_have_to_be, set1)
        set2 = set(fluid_em[('c', 'd')])
        set2_have_to_be = {('c', 'd')}
        self.assertEqual(set2_have_to_be, set2)
        set3 = set(fluid_em[('c', 'h')])
        set3_have_to_be = {('f', 'g'), ('g', 'h'), ('i', 'c'), ('c', 'i'), ('f', 'i'), ('f', 'i')}
        self.assertEqual(set3_have_to_be, set3)

        pos_em = mappings['x_pos_em']
        set1 = set(pos_em[('a', 'c')])
        set2 = set(pos_em[('c', 'h')])
        set3 = set(pos_em[('c', 'd')])
        set1_have_to_be = {('a', 'b')}
        set2_have_to_be = {('c', 'i')}
        set3_have_to_be = {('c', 'd')}
        self.assertEqual(set1_have_to_be, set1)
        self.assertEqual(set2_have_to_be, set2)
        self.assertEqual(set3_have_to_be, set3)

        neg_em = mappings['x_neg_em']
        set1 = set(neg_em[('a', 'c')])
        set2 = set(neg_em[('c', 'h')])
        set3 = set(neg_em[('c', 'd')])
        set1_have_to_be = set()
        set2_have_to_be = {('i', 'c')}
        set3_have_to_be = set()
        self.assertEqual(set1_have_to_be, set1)
        self.assertEqual(set2_have_to_be, set2)
        self.assertEqual(set3_have_to_be, set3)

    def test_7(self):
        G = nx.DiGraph()
        edges = [('a', 'b'), ('c', 'b'), ('c', 'd'), ('f', 'g'), ('g', 'h'),
                 ('i', 'c'), ('c', 'i'), ('f', 'i'), ('f', 'i')]
        xs = [10, -10, 5, 5, 5, 10, 15, -2, -3]
        x_dict = dict(zip(edges, xs))

        G.add_edges_from(edges)
        for n in G.nodes:
            G.nodes[n]['obj'] = vrtxs.HE2_ABC_GraphVertex()
        G2, mappings = mixer2.make_small_graph_for_mixer(G,removeRoot=True)
        x_dict2 = mixer2.build_x_dict_for_reduced_graph(G2, x_dict, **mappings)
        self.assertEqual(x_dict2[('a', 'c')], 10)
        self.assertEqual(x_dict2[('c', 'd')], 5)
        self.assertEqual(x_dict2[('c', 'h')], 5)

    def test_8(self):
        G = nx.DiGraph()
        edges = [('a', 'b'), ('c', 'b'), ('c', 'd'), ('f', 'g'), ('g', 'h'),
                 ('i', 'c'), ('c', 'i'), ('f', 'i'), ('f', 'i')]
        xs = [-15, 15, 5, -20, -20, 25, 5, 13, 7]
        x_dict = dict(zip(edges, xs))

        G.add_edges_from(edges)
        for n in G.nodes:
            G.nodes[n]['obj'] = vrtxs.HE2_ABC_GraphVertex()
        G2, mappings = mixer2.make_small_graph_for_mixer(G,removeRoot=True)
        x_dict2 = mixer2.build_x_dict_for_reduced_graph(G2, x_dict, **mappings)
        self.assertEqual(x_dict2[('a', 'c')], -15)
        self.assertEqual(x_dict2[('c', 'd')], 5)
        self.assertEqual(x_dict2[('c', 'h')], -20)

    def test_9(self):
        G = nx.DiGraph()
        edges = [('a', 'b'), ('c', 'b'), ('c', 'd'), ('f', 'g'), ('g', 'h'),
                 ('i', 'c'), ('c', 'i'), ('f', 'i'), ('f', 'i')]
        xs = [-15, 15, 5, -20, -20, 25, 5, 13, 7]
        x_dict = dict(zip(edges, xs))

        G.add_edges_from(edges)
        for n in G.nodes:
            G.nodes[n]['obj'] = vrtxs.HE2_ABC_GraphVertex()
        cocktails, srcs = mixer2.evalute_network_fluids_with_root(G, x_dict, True)

        self.assertEqual(set(srcs), {'h'})
        ethalon = np.ones(1)
        for e in edges:
            cktl = cocktails[e]
            self.assertEqual(cktl, ethalon)

    def test_10(self):
        G = nx.DiGraph()
        edges = [('a', 'b'), ('c', 'b'), ('c', 'd'), ('f', 'g'), ('g', 'h'),
                 ('i', 'c'), ('c', 'i'), ('f', 'i'), ('f', 'i')]
        xs = [15, -15, -5, 20, 20, -25, -5, -13, -7]
        x_dict = dict(zip(edges, xs))

        G.add_edges_from(edges)
        for n in G.nodes:
            G.nodes[n]['obj'] = vrtxs.HE2_ABC_GraphVertex()
        cocktails, srcs = mixer2.evalute_network_fluids_with_root(G, x_dict, True)

        self.assertEqual(set(srcs), {'a', 'd'})
        eth1 = np.array([1.0, 0])
        eth2 = np.array([0, 1.0])
        eth3 = np.array([0.75, 0.25])
        for e in edges[:2]:
            cktl = cocktails[e]
            self.assertEqual(np.linalg.norm(cktl-eth1), 0)
        for e in edges[2:3]:
            cktl = cocktails[e]
            self.assertEqual(np.linalg.norm(cktl-eth2), 0)
        for e in edges[3:]:
            cktl = cocktails[e]
            self.assertEqual(np.linalg.norm(cktl-eth3), 0)


if __name__ == "__main__":
    test = TestMixer()
    test.test_10()
