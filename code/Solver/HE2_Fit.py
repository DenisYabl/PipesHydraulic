from DFOperations import HE2_DateframeWrapper as model
import scipy.optimize as scop
import numpy as np
import pandas as pd
import os
from Solver.HE2_Solver import HE2_Solver
from Tools.HE2_schema_maker import make_oilpipe_schema_from_OT_dataset, make_calc_df
from Tools.HE2_tools import check_solution, print_solution, print_wells_pressures, cut_single_well_subgraph
import logging
import GraphNodes.HE2_Vertices as vrtx
import matplotlib.pyplot as plt
import random
import itertools
from GraphEdges.HE2_Pipe import HE2_OilPipe
import networkx as nx
from datetime import datetime
from Tools.HE2_tools import check_solution

class HE2_PMNetwork_Model():
    def __init__(self, input_df, method=None, use_bounds=False, fit_version=None):
        # 1750023893
        # 1750024074
        # 1750040926
        # 1750026636
        # 1750037686
        # 1750028316
        # 1750037786

        df = input_df
        df.columns = model.do_upcase_columns_adhoc(df.columns)
        df.loc[df.node_id_start == 1750023893, 'start_P'] = np.nan
        df.loc[df.node_id_start == 1750024074, 'start_P'] = np.nan
        df.loc[df.node_id_start == 1750040926, 'start_P'] = np.nan
        df.loc[df.node_id_start == 1750026636, 'start_P'] = np.nan
        df.loc[df.node_id_start == 1750037686, 'start_P'] = np.nan
        df.loc[df.node_id_start == 1750028316, 'start_P'] = np.nan
        df.loc[df.node_id_start == 1750037786, 'start_P'] = np.nan
        self.input_df = df

        self.row_count = self.input_df.shape[0]
        N = self.row_count
        self.d_weights = np.ones(N)
        self.r_weights = np.ones(N)
        self.d_min = np.zeros(N) + 1e-5
        self.d_max = np.ones(N) * 1.2
        self.r_min = np.zeros(N)
        self.r_max = np.ones(N) * 100
        self.reg_d_deviation_keff = 1
        self.reg_r_deviation_keff = 1
        self.reg_d_dispersion_keff = 1
        self.reg_r_dispersion_keff = 1
        self.reg_total_keff = 1
        self.target_columns = {'start_P': 'result_start_P', 'end_P': 'result_end_P'}
        self.max_it = 100
        self.it = 0
        self.y0 = 0
        self.best_x = None
        self.best_y = 100500
        self.minimize_method = method
        self.use_bounds = use_bounds
        self.fit_version = fit_version
        self.best_rez_df = None

    def regularizator(self):
        N = self.row_count
        terms = np.zeros(4)
        terms[0] = np.linalg.norm((np.ones(N) - self.d_weights)) * self.reg_d_deviation_keff
        terms[1] = np.linalg.norm((np.ones(N) - self.r_weights)) * self.reg_r_deviation_keff
        terms[2] = np.std(self.d_weights) * self.reg_d_dispersion_keff
        terms[3] = np.std(self.r_weights) * self.reg_r_dispersion_keff
        return np.sum(terms) * self.reg_total_keff

    def apply_weights_to_model_long(self, x):
        N = self.row_count
        d_w = x[:N]
        r_w = x[N:]

        d_w = np.array([d_w, self.d_min]).max(axis=0)
        d_w = np.array([d_w, self.d_max]).min(axis=0)

        r_w = np.array([r_w, self.r_min]).max(axis=0)
        r_w = np.array([r_w, self.r_max]).min(axis=0)

        self.d_weights = d_w
        self.r_weights = r_w

        df = self.input_df.copy()
        df.D = df.D * d_w
        df.roughness = df.roughness * r_w
        return df

    def apply_weights_to_model(self, x):
        N = self.row_count
        d_w = x

        d_w = np.array([d_w, self.d_min]).max(axis=0)
        d_w = np.array([d_w, self.d_max]).min(axis=0)

        self.d_weights = d_w
        self.r_weights = np.ones(N)

        df = self.input_df.copy()
        df.D = df.D * d_w
        return df

    def evaluate_y(self, df):
        if df is None:
            return 100500100500
        y = 0
        for col_known, col_rez in self.target_columns.items():
            mask = ~df.loc[:, col_known].isna()
            known_y = df.loc[mask, col_known].values
            result_y = df.loc[mask, col_rez].values
            y += np.linalg.norm(known_y - result_y)
        return y

    def target(self, x):
        self.it += 1
        if self.it > self.max_it:
            raise Exception()
        df = self.apply_weights_to_model(x)
        rez_df = model.do_predict(df)
        y_term = self.evaluate_y(rez_df)
        reg_term = self.regularizator()
        rez = reg_term + y_term
        if (rez < self.best_y) or (self.it % 25 == 0):
            print('       ', self.it, rez)
        if rez < self.best_y:
            self.best_x = x
            self.best_y = rez
            self.best_rez_df = rez_df
        return rez

    def fit_v0(self):
        N = self.row_count
        self.it = 0
        # bnds = None
        # if self.use_bounds:
        # bnds = list(zip(self.d_min, self.d_max)) + list(zip(self.r_min, self.r_max))
        # bnds = list(zip(self.d_min, self.d_max))

        bnds = list(zip(list(np.ones(self.row_count) * 0.35), list(np.ones(self.row_count) * 1.2)))

        # baseline_y = self.target(np.ones(2*N))
        print('------------baseline_y = 240.23416053699975---------------')

        target = self.target
        # x0 = np.ones(2*N) * 0.87

        x0 = np.array([0.64338076, 0.82381463, 1., 0.96528834, 0.96071016,
                       0.95635247, 1., 0.92370242, 1., 0.38750445,
                       0.65503074, 0.95287122, 1., 1., 1.,
                       1., 1., 1., 0.97441723, 0.43580423,
                       0.87804864, 0.9425268, 0.94052543, 0.96233901, 0.9806694,
                       0.96689355, 0.4244979, 1., 0.95379705, 0.95644556,
                       0.68955511, 0.60984356, 0.96924162, 0.65517201, 0.97830807,
                       0.95069912, 0.54673536, 0.952645, 0.92537257, 0.7111141,
                       0.97352889, 0.63872855, 0.95075268, 0.92938683, 0.56179554,
                       0.73654833, 0.72238996, 0.76083504, 0.60044676, 0.91387568,
                       0.59788846, 0.74575036, 0.5212663, 0.46514458, 0.69942246,
                       0.77804001, 0.61124401, 0.80702349, 0.92731837, 0.52031628,
                       0.95585672, 0.92850279, 0.97973397, 0.91646619, 0.91254402,
                       0.91254402, 0.93168514, 0.91253429, 0.96657134, 0.98692722,
                       0.91255251, 0.91253575, 0.9456361, 0.98287691, 0.95358789,
                       0.91255343, 0.94552903])

        try:
            op_result = scop.minimize(target, x0, method=self.minimize_method, bounds=bnds, options=dict(nfev=100500))
        except:
            op_result = scop.OptimizeResult(success=True, fun=self.best_y, x=self.best_x, nfev=self.it)
        return op_result

    def fit_v1(self):
        self.it = 0

        # baseline_y = self.target(np.ones(2*N))
        print('------------baseline_y = 54.896212---------------')

        target = self.target
        x0 = np.array([0.64338076, 0.82381463, 1., 0.96528834, 0.96071016,
                       0.95635247, 1., 0.92370242, 1., 0.38750445,
                       0.65503074, 0.95287122, 1., 1., 1.,
                       1., 1., 1., 0.97441723, 0.43580423,
                       0.87804864, 0.9425268, 0.94052543, 0.96233901, 0.9806694,
                       0.96689355, 0.4244979, 1., 0.95379705, 0.95644556,
                       0.68955511, 0.60984356, 0.96924162, 0.65517201, 0.97830807,
                       0.95069912, 0.54673536, 0.952645, 0.92537257, 0.7111141,
                       0.97352889, 0.63872855, 0.95075268, 0.92938683, 0.56179554,
                       0.73654833, 0.72238996, 0.76083504, 0.60044676, 0.91387568,
                       0.59788846, 0.74575036, 0.5212663, 0.46514458, 0.69942246,
                       0.77804001, 0.61124401, 0.80702349, 0.92731837, 0.52031628,
                       0.95585672, 0.92850279, 0.97973397, 0.91646619, 0.91254402,
                       0.91254402, 0.93168514, 0.91253429, 0.96657134, 0.98692722,
                       0.91255251, 0.91253575, 0.9456361, 0.98287691, 0.95358789,
                       0.91255343, 0.94552903])

        # x0 = np.ones(self.row_count) * 0.87

        try:
            op_result = scop.basinhopping(target, x0, 100500, T=100)
        except:
            op_result = scop.OptimizeResult(success=True, fun=self.best_y, x=self.best_x, nfev=self.it)
        return op_result

    def fit_v2(self):
        self.it = 0

        # baseline_y = self.target(np.ones(2*N))
        print('------------baseline_y = 54.896212---------------')

        target = self.target
        bnds = list(zip(self.d_min, self.d_max))

        try:
            op_result = scop.differential_evolution(target, bounds=bnds, workers=2)
        except:
            op_result = scop.OptimizeResult(success=True, fun=self.best_y, x=self.best_x, nfev=self.it)
        return op_result

    def fit_v3(self):
        self.it = 0

        print('------------baseline_y = 54.896212---------------')
        print('op_result = scop.shgo(target, bounds=bnds)')

        target = self.target
        bnds = list(zip(list(np.ones(self.row_count) * 0.35), list(np.ones(self.row_count) * 1.0)))

        try:
            op_result = scop.shgo(target, bounds=bnds)
        except:
            op_result = scop.OptimizeResult(success=True, fun=self.best_y, x=self.best_x, nfev=self.it)
        return op_result

    def fit_v4(self):
        self.it = 0

        print('------------baseline_y = 54.896212---------------')
        print('op_result = scop.dual_annealing')

        target = self.target
        bnds = list(zip(list(np.ones(self.row_count) * 0.35), list(np.ones(self.row_count) * 1.0)))
        x0 = np.array([0.9368659, 0.90908377, 0.66269683, 0.84984735, 0.94433411,
                       0.94433411, 0.81021361, 0.85424077, 0.89137211, 0.84981522,
                       0.94433411, 0.8586045, 0.94433411, 0.94433411, 0.94433411,
                       0.94433411, 0.94433411, 0.94433411, 0.84980174, 0.44306014,
                       0.87813078, 0.84980174, 0.84980174, 0.84980174, 0.84980174,
                       0.93782366, 0.84459252, 0.93500233, 0.84980174, 0.84980174,
                       0.78761383, 0.59268382, 0.84980174, 0.62411044, 0.73925876,
                       0.82262774, 0.55721666, 0.73925876, 0.69122549, 0.59343754,
                       0.84980174, 0.64804504, 0.65627441, 0.84980174, 0.59175647,
                       0.70158817, 0.58099541, 0.61674218, 0.61099233, 0.7366823,
                       0.4614739, 0.66874503, 0.55128615, 0.47035048, 0.61580247,
                       0.71251471, 0.53147553, 0.64558384, 0.64558384, 0.52497384,
                       0.84980174, 0.84980174, 0.83987027, 0.71487156, 0.84980173,
                       0.84980173, 0.756754, 0.53656503, 0.84980173, 0.84980174,
                       0.84980846, 0.84979832, 0.8498047, 0.84980606, 0.84979465,
                       0.84981028, 0.84980687])

        try:
            op_result = scop.dual_annealing(target, bounds=bnds, x0=x0, maxiter=100500)
        except:
            op_result = scop.OptimizeResult(success=True, fun=self.best_y, x=self.best_x, nfev=self.it)
        return op_result

    def fit_v5(self):
        self.it = 0

        print('------------baseline_y = 54.896212---------------')
        print('op_result = scop.dual_annealing')

        target = self.target
        bnds = list(zip(list(np.ones(self.row_count) * 0.35), list(np.ones(self.row_count) * 1.1)))
        x0 = np.array([0.64338076, 0.82381463, 1., 0.96528834, 0.96071016,
                       0.95635247, 1., 0.92370242, 1., 0.38750445,
                       0.65503074, 0.95287122, 1., 1., 1.,
                       1., 1., 1., 0.97441723, 0.43580423,
                       0.87804864, 0.9425268, 0.94052543, 0.96233901, 0.9806694,
                       0.96689355, 0.4244979, 1., 0.95379705, 0.95644556,
                       0.68955511, 0.60984356, 0.96924162, 0.65517201, 0.97830807,
                       0.95069912, 0.54673536, 0.952645, 0.92537257, 0.7111141,
                       0.97352889, 0.63872855, 0.95075268, 0.92938683, 0.56179554,
                       0.73654833, 0.72238996, 0.76083504, 0.60044676, 0.91387568,
                       0.59788846, 0.74575036, 0.5212663, 0.46514458, 0.69942246,
                       0.77804001, 0.61124401, 0.80702349, 0.92731837, 0.52031628,
                       0.95585672, 0.92850279, 0.97973397, 0.91646619, 0.91254402,
                       0.91254402, 0.93168514, 0.91253429, 0.96657134, 0.98692722,
                       0.91255251, 0.91253575, 0.9456361, 0.98287691, 0.95358789,
                       0.91255343, 0.94552903])
        x0 = np.ones(self.row_count) * 0.9

        op_result = None
        for i in range(20):
            op_result = scop.dual_annealing(target, bounds=bnds, x0=x0, maxfun=1000)
            x0 = op_result.x
        return op_result

    def fit(self):
        fits = {0: self.fit_v0, 1: self.fit_v1, 2: self.fit_v2, 3: self.fit_v3, 4: self.fit_v4, 5: self.fit_v5}
        meth = fits.get(self.fit_version, self.fit_v0)
        return meth()

def heatmap(data, row_labels, col_labels, ax=None, cbar_kw={}, cbarlabel="", **kwargs):
    """
    Create a heatmap from a numpy array and two lists of labels.

    Parameters
    ----------
    data
        A 2D numpy array of shape (N, M).
    row_labels
        A list or array of length N with the labels for the rows.
    col_labels
        A list or array of length M with the labels for the columns.
    ax
        A `matplotlib.axes.Axes` instance to which the heatmap is plotted.  If
        not provided, use current axes or create a new one.  Optional.
    cbar_kw
        A dictionary with arguments to `matplotlib.Figure.colorbar`.  Optional.
    cbarlabel
        The label for the colorbar.  Optional.
    **kwargs
        All other arguments are forwarded to `imshow`.
    """

    if not ax:
        ax = plt.gca()

    # Plot the heatmap
    im = ax.imshow(data, **kwargs)

    # Create colorbar
    cbar = ax.figure.colorbar(im, ax=ax, **cbar_kw)
    cbar.ax.set_ylabel(cbarlabel, rotation=-90, va="bottom")

    # We want to show all ticks...
    ax.set_xticks(np.arange(data.shape[1]))
    ax.set_yticks(np.arange(data.shape[0]))
    # ... and label them with the respective list entries.
    ax.set_xticklabels(col_labels)
    ax.set_yticklabels(row_labels)

    # Let the horizontal axes labeling appear on top.
    ax.tick_params(top=True, bottom=False,
                   labeltop=True, labelbottom=False)

    # Rotate the tick labels and set their alignment.
    plt.setp(ax.get_xticklabels(), rotation=-30, ha="right",
             rotation_mode="anchor")

    # Turn spines off and create white grid.
    # ax.spines[:].set_visible(False)

    ax.set_xticks(np.arange(data.shape[1]+1)-.5, minor=True)
    ax.set_yticks(np.arange(data.shape[0]+1)-.5, minor=True)
    # ax.grid(which="minor", color="w", linestyle='-', linewidth=3)
    ax.tick_params(which="minor", bottom=False, left=False)

    return im, cbar


class HE2_OilGatheringNetwork_Model():
    def __init__(self, commondata_folder='', wells_folder = ''):
        # tree = os.walk(folder)
        # for fld, subfolders, files in tree:
        #     print(fld, subfolders, files)
        self.commondata_folder = commondata_folder
        self.wells_folder = wells_folder
        self.solvers = dict()
        self.calc_dfs = dict()
        self.original_dfs = dict()
        self.graphs = dict()
        self.pad_well_list = None
        self.pad_wells_dict = None
        self.fact = dict()
        self.outlayers = dict()
        self.result = dict()
        self.last_well_result = dict()
        self.well_result_before = dict()
        self.well_result_after = dict()
        self.initial_x = dict()
        self.total_target_cnt, self.not_solved = 0, 0
        self.bad_wells = []
        self.prefit_params = dict()
        self.df_file_suffixes = list(map(str, range(1, 31)))
        self.N = len(self.df_file_suffixes)
        self.original_diams = dict()
        self.last_it_count = 0
        self.ignore_watercut = False


    def gimme_original_df(self, i):
        if i in self.original_dfs:
            return self.original_dfs[i]
        folder = self.wells_folder
        suffix = self.df_file_suffixes[i]
        filename = f'{folder}/DNS2_with_wells_{suffix}.csv'
        df = pd.read_csv(filename)
        self.original_dfs[i] = df
        return df

    def gimme_calc_df(self, i):
        if i in self.calc_dfs:
            return self.calc_dfs[i]
        folder = self.wells_folder
        suffix = self.df_file_suffixes[i]
        calc_df_filename = f'{folder}/calc_df_{suffix}.csv'
        try:
            calc_df = pd.read_csv(calc_df_filename)
            self.calc_dfs[i] = calc_df
            return calc_df
        except:
            pass
        original_df = self.gimme_original_df(i)
        calc_df = make_calc_df(original_df, self.commondata_folder)
        calc_df.to_csv(calc_df_filename)
        self.calc_dfs[i] = calc_df
        return calc_df

    def gimme_graph(self, i):
        if i in self.graphs:
            return self.graphs[i]
        df = self.gimme_original_df(i)
        calc_df = self.gimme_calc_df(i)
        G, _, __ = make_oilpipe_schema_from_OT_dataset(df, self.commondata_folder, calc_df, ignore_Watercut=self.ignore_watercut)
        self.graphs[i] = G
        return self.graphs[i]

    def gimme_solver(self, i):
        if i in self.solvers:
            return self.solvers[i]

        G = self.gimme_graph(i)
        solver = HE2_Solver(G)
        self.solvers[i] = solver
        return solver

    def gimme_wells(self):
        if self.pad_wells_dict:
            return self.pad_wells_dict
        dfs = []
        for i in range(self.N):
            df = self.gimme_original_df(i)
            df = df[['juncType', 'padNum', 'wellNum']]
            df = df[df.juncType == 'oilwell']
            dfs += [df]
        df = pd.concat(dfs).drop_duplicates()[['padNum', 'wellNum']]
        raw_pw_list = list(df.to_records(index=False))
        self.pad_well_list = []
        for pad, well in raw_pw_list:
            self.pad_well_list += [(pad, int(well))]
        pad_wells_dict = dict()
        for (pad, well) in self.pad_well_list:
            wlist = pad_wells_dict.get(pad, [])
            wlist += [well]
            pad_wells_dict[pad] = wlist
        self.pad_wells_dict = pad_wells_dict
        return self.pad_wells_dict

    def grab_fact(self):
        pad_wells_dict = self.gimme_wells()
        p_zab = dict()
        p_intake = dict()
        p_head = dict()
        q_well = dict()
        freq = dict()

        cols = ['juncType', 'padNum', 'wellNum', 'zaboy_pressure','input_pressure', 'buffer_pressure', 'debit', 'frequency']
        dfs = []
        for i in range(self.N):
            df = self.gimme_original_df(i)
            df = df[cols]
            df['N'] = i
            dfs += [df]
        fact_df = pd.concat(dfs)
        for pad in pad_wells_dict:
            pad_df = fact_df[fact_df.padNum == pad]
            wells = pad_wells_dict[pad]
            for well in wells:
                df = pad_df[pad_df.wellNum == well]
                df = df.sort_values(by=['N'])
                p_zab[(pad, well)] = df.zaboy_pressure.values
                p_intake[(pad, well)] = df.input_pressure.values
                p_head[(pad, well)] = df.buffer_pressure.values
                q_well[(pad, well)] = df.debit.values
                freq[(pad, well)] = df.frequency.values

        return dict(p_zab=p_zab, p_intake=p_intake, p_head=p_head, q_well=q_well, freq=freq)

    def grab_results(self):
        pad_wells_dict = self.gimme_wells()
        p_zab = dict()
        p_intake = dict()
        p_head = dict()
        q_well = dict()
        for pad in pad_wells_dict:
            wells = pad_wells_dict[pad]
            for well in wells:
                key = (pad, well)
                zab_name = f"PAD_{pad}_WELL_{well}_zaboi"
                pump_name = f"PAD_{pad}_WELL_{well}_pump_intake"
                head_name = f"PAD_{pad}_WELL_{well}_wellhead"
                p_zab[key] = np.zeros(self.N)
                p_intake[key] = np.zeros(self.N)
                p_head[key] = np.zeros(self.N)
                q_well[key] = np.zeros(self.N)

                for i in range(self.N):
                    solver = self.gimme_solver(i)
                    G = solver.graph
                    p_zab[key][i] = G.nodes[zab_name]['obj'].result['P_bar']
                    p_intake[key][i] = G.nodes[pump_name]['obj'].result['P_bar']
                    p_head[key][i] = G.nodes[head_name]['obj'].result['P_bar']
                    obj = G[zab_name][pump_name]['obj']
                    q_well[key][i] = 86400 * obj.result['x'] / obj.result['liquid_density']
        return dict(p_zab=p_zab, p_intake=p_intake, p_head=p_head, q_well=q_well)

    def solve_em_all(self):
        self.last_it_count = 0
        success = np.zeros(self.N)
        for i in range(self.N):
            solver = self.gimme_solver(i)
            if (0, 0, i) in self.initial_x:
                solver.initial_edges_x = self.initial_x[(0, 0, i)]
            solver.solve(threshold=0.5, mix_fluids=True)
            validity = check_solution(solver.schema)
            print(validity)
            success[i] = solver.op_result.success
            self.last_it_count += solver.op_result.nfev
            self.initial_x[(0, 0, i)] = solver.edges_x.copy()
        return success

    def plot_fact_and_results(self, keys_to_plot=('head', 'intake', 'bottom', 'debit'), wells=(), pads=()):
        fig = plt.figure(constrained_layout=True, figsize=(8, 8))
        ax = fig.add_subplot(1, 1, 1)
        ax.set_title(keys_to_plot)
        ax.set_xlabel('fact')
        ax.set_ylabel('result')
        colors = plt.get_cmap('tab20c').colors
        colors = list(colors) * 3
        random.shuffle(colors)
        for i, (pad, well) in enumerate(self.pad_well_list):
            alpha = 0.2
            if well in wells or pad in pads:
                alpha = 0.8
            if 'head' in keys_to_plot:
                ax.plot(self.fact['p_head'][(pad, well)], self.result['p_head'][(pad, well)], color=colors[i], linewidth=1, alpha=alpha)
            if 'intake' in keys_to_plot:
                ax.plot(self.fact['p_intake'][(pad, well)], self.result['p_intake'][(pad, well)], color=colors[i], linewidth=1, alpha=alpha)
            if 'bottom' in keys_to_plot:
                ax.plot(self.fact['p_zab'][(pad, well)], self.result['p_zab'][(pad, well)], color=colors[i], linewidth=1, alpha=alpha)
            if 'debit' in keys_to_plot:
                ax.plot(self.fact['q_well'][(pad, well)], self.result['q_well'][(pad, well)], color=colors[i], linewidth=1, alpha=alpha)

        plt.show()

    def calc_well_score(self, pad, well):
        debit_scale = 150
        head_scale = 20
        intake_scale = 50
        bottom_scale = 80

        score = 0
        score += abs(self.fact['p_head'][(pad, well)] - self.result['p_head'][(pad, well)]) / head_scale
        score += abs(self.fact['p_intake'][(pad, well)] - self.result['p_intake'][(pad, well)]) / intake_scale
        score += abs(self.fact['p_zab'][(pad, well)] - self.result['p_zab'][(pad, well)]) / bottom_scale
        score += abs(self.fact['q_well'][(pad, well)] - self.result['q_well'][(pad, well)])  / debit_scale
        return score

    def calc_well_score2(self, pad, well, result):
        weights = dict(p_head=1, p_intake=1, p_zab=1, q_well=1)

        ol_intake = self.outlayers['p_intake'][(pad, well)]
        # ol_bottom = self.outlayers['p_zab'][(pad, well)]
        ol_debit = self.outlayers['q_well'][(pad, well)]
        ol_head = self.outlayers['p_head'][(pad, well)]

        ip_error = self.fact['p_intake'][(pad, well)] - result['p_intake'][(pad, well)]
        ip_error = np.linalg.norm(ip_error[~ol_intake])
        score = ip_error * weights['p_intake']

        # bp_error = self.fact['p_zab'][(pad, well)] - result['p_zab'][(pad, well)]
        # bp_error = np.linalg.norm(bp_error[~ol_bottom])
        # score += bp_error * weights['p_zab']

        hp_error = self.fact['p_head'][(pad, well)] - result['p_head'][(pad, well)]
        hp_error = np.linalg.norm(hp_error[~ol_head])
        score += hp_error * weights['p_head']

        q_error = self.fact['q_well'][(pad, well)] - result['q_well'][(pad, well)]
        q_error = np.linalg.norm(q_error[~ol_debit])
        score += q_error * weights['q_well']

        return score

    def prefit_all(self):
        rez = []
        self.gimme_wells()
        pad_well_list = self.pad_well_list
        for it, (pad, well) in enumerate(pad_well_list):
            if (well, pad) in self.bad_wells:
                continue
            # if not (pad, well) == ('39', 619):
            #     continue
            G0 = self.gimme_graph(0)
            nodes = [f'PAD_{pad}_WELL_{well}']
            nodes += [f'PAD_{pad}_WELL_{well}_zaboi']
            nodes += [f'PAD_{pad}_WELL_{well}_pump_intake']
            nodes += [f'PAD_{pad}_WELL_{well}_pump_outlet']
            nodes += [f'PAD_{pad}_WELL_{well}_wellhead']
            well_G, _ = cut_single_well_subgraph(G0, pad, well, nodes)
            solver = HE2_Solver(well_G)
            pump_obj = well_G[nodes[2]][nodes[3]]['obj']
            plast_obj = well_G[nodes[0]][nodes[1]]['obj']
            bounds = ((0.1, 3 * plast_obj.Productivity), (0, 1.2))
            path = []
            def target_fit(x):
                nonlocal pump_obj, plast_obj, solver, well_G, pad, well, nodes, bounds, path
                path += [x]
                productivity = x[0]
                plast_obj.Productivity = productivity
                pump_keff = x[1]
                pump_obj.change_stages_ratio(pump_keff)
                score = self.score_well(solver, well_G, pad, well, nodes)
                # print(x, score)
                return score

            x0 = np.array([plast_obj.Productivity, 1])
            target_fit(x0)
            well_res_before = self.last_well_result.copy()

            xs = np.linspace(start=bounds[0][0], stop=bounds[0][1], num=15)
            ys = np.linspace(start=bounds[1][0], stop=bounds[1][1], num=15)
            zs = np.zeros((len(xs), len(ys)))

            for i, x in enumerate(xs):
                for j, y in enumerate(ys):
                    zs[i, j] = target_fit(np.array([x, y]))

            i, j = np.unravel_index(np.argmin(zs, axis=None), zs.shape)
            x0 = np.array([xs[i], ys[j]])
            op_result = scop.minimize(target_fit, x0)
            best_x = op_result.x
            target_fit(best_x)
            well_res_after = self.last_well_result.copy()

            zs[zs==100500] = -1
            self.plot_single_well_chart(pad, well, well_res_before, well_res_after, well_G, nodes, op_result, xs, ys, zs, path)
            rez += [(well, pad, best_x, op_result.fun)]
            print(well, pad, best_x, op_result.fun)

        return rez



    def score_well(self, solver, well_G, pad, well, nodes):
        weights = dict(p_head=1, p_intake=1, p_zab=1, q_well=1)

        rez_p_i = np.ones(self.N) * 100500
        rez_p_b = np.ones(self.N) * 100500
        rez_q = np.ones(self.N) * 100500

        q0 = self.fact['q_well'][(pad, well)][0]
        known_Q = {nodes[0]: 900 * q0 / 86400}
        solver.set_known_Q(known_Q)
        fact_phs = self.fact['p_head'][(pad, well)]
        pump_obj = well_G[nodes[2]][nodes[3]]['obj']
        freqs = self.fact['freq'][(pad, well)]
        ol_freqs = self.outlayers['freq'][(pad, well)]
        freq = 50
        for i in range(self.N):
            if not ol_freqs[i]:
                freq = freqs[i]
            pump_obj.changeFrequency(freq) # Set last non-outlayer frequency

            well_G.nodes[nodes[-1]]['obj'].value = fact_phs[i]
            if (pad, well, i) in self.initial_x:
                solver.initial_edges_x = self.initial_x[(pad, well, i)]

            solver.solve(mix_fluids=True, threshold=0.5, it_limit=30)
            self.total_target_cnt += solver.op_result.nfev
            if not solver.op_result.success:
                self.not_solved += 1
                return 100500
                # continue
            self.initial_x[(pad, well, i)] = solver.edges_x.copy()

            rez_p_i[i] = well_G.nodes[nodes[2]]['obj'].result['P_bar']
            rez_p_b[i] = well_G.nodes[nodes[1]]['obj'].result['P_bar']
            rez = well_G[nodes[1]][nodes[2]]['obj'].result
            rez_q[i] = 86400 * rez['x'] / rez['liquid_density']

        self.last_well_result['p_zab'] = rez_p_b
        self.last_well_result['p_intake'] = rez_p_i
        self.last_well_result['q_well'] = rez_q

        ol_intake = self.outlayers['p_intake'][(pad, well)]
        ip_error = self.fact['p_intake'][(pad, well)] - rez_p_i
        ip_error = np.linalg.norm(ip_error[~ol_intake])
        score = ip_error * weights['p_intake']

        ol_debit = self.outlayers['q_well'][(pad, well)]
        q_error = self.fact['q_well'][(pad, well)] - rez_q
        q_error = np.linalg.norm(q_error[~ol_debit])
        score += q_error * weights['q_well']

        score += np.linalg.norm(self.fact['p_zab'][(pad, well)] - rez_p_b) * weights['p_zab']
        return score


    def plot_single_well_chart(self, pad, well, well_res_before, well_res_after, well_G, nodes, op_result, xs, ys, zs, path):
        fig = plt.figure(figsize=(8, 8))
        ax = fig.add_subplot(1, 1, 1)
        ax.set_title(f'Pad {pad} well {well} y {op_result.fun:.3f}, prod {op_result.x[0]:.3f}  pump_keff {op_result.x[1]:.3f}')
        # Xs, Ys = np.meshgrid(xs, ys)
        # Zs = zs.reshape(len(ys), len(xs))
        # im = plt.pcolormesh(Xs, Ys, Zs, cmap='plasma')
        # fig.colorbar(im, ax=ax, fraction=0.05, shrink=0.5)

        xs = np.round(xs, 3)
        ys = np.round(ys, 3)
        im, cbar = heatmap(zs, ys, xs, ax=ax, cmap="plasma")

        # i, j = np.unravel_index(np.argmin(zs, axis=None), zs.shape)
        # ax.plot(xs[i], ys[j], color='w', marker='X', markersize=15)

        ax.set_xlabel('prod')
        ax.set_ylabel('pump keff')
        plt.show()

        fig = plt.figure(constrained_layout=True, figsize=(8, 8))
        ax = fig.add_subplot(1, 1, 1)
        WC = well_G.nodes[nodes[0]]['obj'].fluid.oil_params.volumewater_percent
        ax.set_title(f'Pad {pad} well {well} WC {WC}% y {op_result.fun:.3f}, prod {op_result.x[0]:.3f}  pump_keff {op_result.x[1]:.3f}')
        # ax.set_xlabel('time')
        ax.set_ylabel('P, Q')
        colors = plt.get_cmap('Set2').colors

        xs = np.array(range(self.N))
        kwargs = dict(linewidth=1, alpha=0.9)
        mask = ~self.outlayers['p_intake'][(pad, well)]
        hndl1 = ax.plot(xs[mask], self.fact['p_intake'][(pad, well)][mask], color=colors[0], label='Fact: intake', **kwargs)
        hndl2 = ax.plot(xs, self.fact['p_zab'][(pad, well)], color=colors[1], label='Fact: zaboi', **kwargs)

        mask = ~self.outlayers['q_well'][(pad, well)]
        hndl3 = ax.plot(xs[mask], self.fact['q_well'][(pad, well)][mask], color=colors[2], label='Fact: debit', **kwargs)

        kwargs = dict(linewidth=5, alpha=1)
        mask = ~self.outlayers['freq'][(pad, well)]
        hndl6 = ax.plot(xs[mask], self.fact['freq'][(pad, well)][mask], color=colors[4], label='Fact: freq', **kwargs)

        # kwargs = dict(linewidth=1, alpha=0.9, linestyle=':')
        # hndl4 = ax.plot(xs, well_res_before['p_intake'], color=colors[0], label='Before: intake', **kwargs)
        # ax.plot(xs, well_res_before['p_zab'], color=colors[1], **kwargs)
        # ax.plot(xs, well_res_before['q_well'], color=colors[2], **kwargs)

        kwargs = dict(linewidth=1, alpha=0.9, linestyle='--')
        hndl5 = ax.plot(xs, well_res_after['p_intake'], color=colors[0], label='Predict: intake', **kwargs)
        ax.plot(xs, well_res_after['p_zab'], color=colors[1], **kwargs)
        ax.plot(xs, well_res_after['q_well'], color=colors[2], **kwargs)

        # handles = hndl1 + hndl2 + hndl3 + hndl6 + hndl4 + hndl5
        handles = hndl1 + hndl2 + hndl3 + hndl6 + hndl5
        ax.legend(handles=handles, loc='upper center')
        plt.show()

    def plot_fit_log(self, fit_log=None, filename=''):
        if not fit_log:
            if not filename:
                return
            f = open(filename)
            fit_log = []
            for s in f:
                sss = s.split(' ')
                record = tuple(map(float, sss))
                fit_log += [record]


        N = len(fit_log)
        x, y1, y2, y3 = np.zeros(N), np.zeros(N), np.zeros(N), np.zeros(N)
        prev, it_num = 0, 0
        for i, record in enumerate(fit_log):
            it_num, seconds, it_count, best_y = record
            x[i] = seconds
            y1[i] = seconds - prev
            prev = seconds
            y2[i] = it_count / self.N
            y3[i] = best_y

        fig = plt.figure(figsize=(8, 8))
        ax = fig.add_subplot(1, 1, 1)

        kwargs = dict(linewidth=0.1, alpha=0.99)
        hndl1 = ax.plot(x, y1, label='time, sec', **kwargs)
        hndl2 = ax.plot(x, y2, label='solver avg iterations', **kwargs)
        hndl3 = ax.plot(x, y3, label='Y, parrots')
        handles = hndl1 + hndl2 + hndl3
        ax.legend(handles=handles, loc='center left')
        ax.set_ylim(0, 45)
        ax.set_title(f'After {int(it_num)} iterations')
        plt.show()


    def fill_outlayers(self):
        ol_freq = dict()
        ol_intake = dict()
        ol_debit = dict()
        ol_head = dict()
        false_mask = np.zeros(self.N) == np.ones(self.N)
        for pad, well in self.pad_well_list:
            if well in self.bad_wells:
                ol_freq[(pad, well)] = false_mask
                ol_intake[(pad, well)] = false_mask
                ol_debit[(pad, well)] = false_mask
                ol_head[(pad, well)] = false_mask
                continue

            freq = self.fact['freq'][(pad, well)]
            mask1 = freq > 99
            mask2 = freq < 40
            mask3 = np.isnan(freq)
            mask = mask1 | mask2 | mask3
            ol_freq[(pad, well)] = mask

            intake = self.fact['p_intake'][(pad, well)]
            mask1 = intake > 70
            mask2 = intake < 20
            mask3 = np.isnan(intake)
            mask = mask1 | mask2 | mask3
            ol_intake[(pad, well)] = mask

            debit = self.fact['q_well'][(pad, well)]
            mask1 = debit > 1000
            mask2 = debit < 0.5
            mask3 = np.isnan(debit)
            mask = mask1 | mask2 | mask3
            ol_debit[(pad, well)] = mask

            if pad == 39:
                ol_head[(pad, well)] = false_mask
            else:
                ol_head[(pad, well)] = ~false_mask

        self.outlayers['freq'] = ol_freq
        self.outlayers['p_intake'] = ol_intake
        self.outlayers['q_well'] = ol_debit
        self.outlayers['p_head'] = ol_head


    def gimme_prefit_params(self):
        if self.prefit_params:
           return self.prefit_params

        prefit = dict()
        folder = self.folder
        filename = f'{folder}data/oilgathering_fit/prefit.csv'
        df = pd.read_csv(filename)
        lst = list(df.to_records(index=False))
        for well, pad, prod, stages in lst:
            prefit[(pad, well)] = dict(K_prod=prod, K_pump=stages)
        self.prefit_params = prefit

        return self.prefit_params

    def save_original_diams(self):
        solver = self.gimme_solver(0)
        G = solver.graph
        for edge in G.edges:
            u, v = edge
            if 'PAD' in u or 'PAD' in v:
                continue
            obj = G[u][v]['obj']
            if not isinstance(obj, HE2_OilPipe):
                continue
            ds = []
            for seg in obj.segments:
                ds += [seg.inner_diam_m]
                self.original_diams[edge] = ds
        pass


    def apply_params(self, params):
        self.gimme_wells()

        if (0, 0) in params:
            network_params = params[(0, 0)]
            diam_keff = network_params['diam_keff']
            for (u, v) in self.original_diams:
                for i in range(self.N):
                    solver = self.gimme_solver(i)
                    G = solver.graph
                    obj = G[u][v]['obj']
                    ds = self.original_diams[(u, v)]
                    for seg, d in zip(obj.segments, ds):
                        seg.inner_diam_m = d * diam_keff

        pad_well_list = self.pad_well_list
        for (pad, well) in pad_well_list:
            if not (pad, well) in params:
                continue

            prms = params[(pad, well)]
            nodes = [f'PAD_{pad}_WELL_{well}']
            nodes += [f'PAD_{pad}_WELL_{well}_zaboi']
            nodes += [f'PAD_{pad}_WELL_{well}_pump_intake']
            nodes += [f'PAD_{pad}_WELL_{well}_pump_outlet']

            ol_freq = self.outlayers['freq'][(pad, well)]

            for i in range(self.N):
                solver = self.gimme_solver(i)
                G = solver.graph
                pump_obj = G[nodes[2]][nodes[3]]['obj']
                plast_obj = G[nodes[0]][nodes[1]]['obj']
                pump_obj.change_stages_ratio(prms['K_pump'])
                plast_obj.Productivity = prms['K_prod']
                if ol_freq[i]:
                    pump_obj.changeFrequency(prms['Freq_0'])

    def dump_x(self, x, fitlog):
        f = open('x.txt', 'w')
        for key in x:
            print(*key, x[key], file=f)
        f.close()
        f = open('fitlog.txt', 'w')
        for record in fitlog:
            print(*record, file=f)
        f.close()
        self.plot_fit_log(fitlog)



    def greed_optimization(self):
        self.ignore_watercut = False
        self.N = 1
        self.gimme_wells()
        M = len(self.pad_well_list)

        def target(x):
            if x is not None:
                self.apply_params(x)
            success = self.solve_em_all()
            if not success.all():
                return 100500
            results = self.grab_results()
            score = np.zeros((M, self.N))
            for i, (pad, well) in enumerate(self.pad_well_list):
                well_score = self.calc_well_score2(pad, well, results)
                score[i] = well_score
            result = np.average(score)
            return result
        x = None
        naive_score = target(x)
        self.save_original_diams() # Not earlier, just here. Cause solver.graph is None before first solve(). And we prefer solver.graph cause it is DiGraph, not MultiDiGraph

        prefit_params = self.gimme_prefit_params()
        for key in prefit_params:
            ol_freq = self.outlayers['freq'][key]
            if ol_freq.any():
                prefit_params[key]['Freq_0'] = 50
        prefit_params[(0, 0)] = dict(diam_keff=0.85)

        x = prefit_params
        prefit_score = target(x)

        print(naive_score, prefit_score)

        print_solution(self.gimme_graph(0))

        return

        random.seed = 42
        order = [(0, 0, 'diam_keff')] * 2
        for pad, well in self.pad_well_list:
            order += [(pad, well, 'K_prod'), (pad, well, 'K_pump')]
            if 'Freq_0' in prefit_params[(pad, well)]:
                order += [(pad, well, 'Freq_0')]

        orders = []
        for i in range(10):
            random.shuffle(order)
            orders += order * 10
        order = orders

        need_dump = False
        step = 0.01
        x = prefit_params.copy()
        best_y = prefit_score
        fitlog = []

        start_time = datetime.now()
        for i, (pad, well, param) in enumerate(order):
            if i == 15:
                break
            if i%50 == 2:
                need_dump = True
            fitlog += [(i, (datetime.now()-start_time).seconds, self.last_it_count, best_y)]


            print(f'{i:#4}/{len(order)}  {pad:#2}  {well:#4}   {param:>8}', end=' ')
            x1 = x.copy()
            x1[(pad, well)][param] = x1[(pad, well)][param] * (1 + step)
            y = target(x1)
            if y >= 100500:
                print(y)

            if y < best_y:
                best_y = y
                x = x1
                print(f'{best_y:3.7}')
                if need_dump:
                    self.dump_x(x, fitlog)
                    need_dump = False
                continue

            x2 = x.copy()
            x2[(pad, well)][param] = x2[(pad, well)][param] * (1 - step)
            y = target(x2)
            if y >= 100500:
                print(y)

            if y < best_y:
                best_y = y
                x = x2
                print(f'{best_y:3.7}')
                if need_dump:
                    self.dump_x(x, fitlog)
                    need_dump = False
                continue
            print()

    def save_prefit_results_to_csv(self, rez, filename):
        df = pd.DataFrame(columns=['wellNum','padNum','K_prod','K_pump'])
        for item in rez:
            x = item[2]
            row = dict(wellNum=item[0],padNum=item[1],K_prod=x[0],K_pump=x[1])
            df = df.append(row, ignore_index=True)
        df.to_csv(filename, index=False)


# bad_wells = [(738, 33), (567, 39), (4532, 49), (2630, 49), (1579, 57), (3118, 57)]


if __name__ == '__main__':
    model = HE2_OilGatheringNetwork_Model(commondata_folder="../../CommonData/", wells_folder='../../data/fit 22-06-2021/')
    model.fact = model.grab_fact()
    # model.bad_wells = bad_wells
    model.fill_outlayers()
    # model.plot_fit_log(filename='fitlog_N5_WC_3700.txt')
    # model.greed_optimization()

    rez = model.prefit_all()
    model.save_prefit_results_to_csv(rez, 'model weights.csv')
