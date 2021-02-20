from matplotlib import pyplot as plt
import networkx as nx
from Tools import HE2_tools as tools
import pandas as pd

class PipelinesPlotter():
    source_col = "node_id_start"
    target_col = "node_id_end"
    node_id_col = 'node_id'
    x_col = 'x'
    y_col = 'y'
    result_P_col = 'result_P'
    well_type, junc_type, pad_type, kns_type = 2, 3, 1, 8
    node_types = [well_type, junc_type, pad_type, kns_type]
    nodes_settings = {well_type: dict(node_color='blue',  node_size=10),
                      junc_type: dict(node_color='black', node_size=5),
                      pad_type:  dict(node_color='blue',  node_size=25),
                      kns_type:  dict(node_color='red',   node_size=30)}
    node_label_settings = { well_type: {},
                            junc_type: {result_P_col: ':,%.2f'},
                            pad_type:  {node_id_col: ':,%s', result_P_col: ':,%.2f'},
                            kns_type:  {}}


    def __init__(self, input_df):
        self.input_df = input_df
        self.df_edges, self.df_nodes = tools.split_result_df_to_pipes_and_nodes(input_df)

    def make_nodes_position_dict(self):
        pos = self.df_nodes[["node_id", 'x', 'y']]
        sp = dict(zip(pos.node_id, (zip(pos.x, pos.y))))
        return sp

    def check_df(self):
        pass

    def make_graph(self):
        df = self.input_df[[PipelinesPlotter.source_col, PipelinesPlotter.target_col]]
        df.columns = ["source", "target"]
        G = nx.from_pandas_edgelist(df)
        return G

    def gimme_all_nodes_of_type(self, node_type):
        df = self.df_nodes
        mask = df.node_type == node_type
        print(sum(mask))
        df = df[mask]
        nodes = list(df.node_id)
        return nodes

    def make_labels_dict(self, type):
        sep = ' '
        node_id_col = PipelinesPlotter.node_id_col
        sets = PipelinesPlotter.node_label_settings[type]
        cols = sets.keys()
        if len(cols)==0:
            return dict()
        df = self.df_nodes.dropna()
        ids = list(df[node_id_col])
        df = df[cols]
        vals = list(df.apply(lambda rec: sep.join(map('{:,.2f}'.format, rec)), axis=1))
        rez = dict(zip(ids, vals))
        return rez


    def plot_graph(self):
        G = self.make_graph()
        fig = plt.figure(constrained_layout=True, figsize=(20, 20))
        ax = fig.add_subplot(1, 1, 1)
        ax.set_aspect('equal')
        pos = self.make_nodes_position_dict()
        ns = PipelinesPlotter.nodes_settings
        all_nodes = set()
        for k, params in ns.items():
            nodelist = self.gimme_all_nodes_of_type(k)
            s = all_nodes & set(nodelist)
            l = len(s)
            assert l==0
            all_nodes |= set(nodelist)
            lbls = self.make_labels_dict(k)
            nx.draw_networkx_nodes(G, pos=pos, nodelist=nodelist, ax=ax, **params)
            nx.draw_networkx_labels(G, pos, lbls, ax=ax)

        nx.draw_networkx_edges(G, pos, alpha=0.5, width=2, ax=ax)
        ax.tick_params(left=True, bottom=True, labelleft=True, labelbottom=True)
        limits=plt.axis('on') # turns on axis


def draw_graph(in_df, filename):
    plotter = PipelinesPlotter(in_df)
    plotter.plot_graph()
    plt.savefig(filename, format='svg')