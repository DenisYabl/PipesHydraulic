from DFOperations import HE2_DateframeWrapper as model
import scipy.optimize as scop
import numpy as np
import pandas as pd
import os
from Solver.HE2_Solver import HE2_Solver
from Tools.HE2_schema_maker import make_oilpipe_schema_from_OT_dataset, make_calc_df
from Tools.HE2_tools import check_solution
import logging
import GraphNodes.HE2_Vertices as vrtx

class HE2_OilGatheringNetwork_Model():
    def __init__(self, folder):
        # tree = os.walk(folder)
        # for fld, subfolders, files in tree:
        #     print(fld, subfolders, files)
        self.folder = folder
        self.solvers = dict()
        self.calc_dfs = dict()
        self.original_dfs = dict()
        self.graphs = dict()
        self.N = 5

    def gimme_original_df(self, i):
        if i in self.original_dfs:
            return self.original_dfs[i]
        folder = self.folder
        filename = f'{folder}data/oilgathering_fit/DNS2_with_wells_{i}.csv'
        df = pd.read_csv(filename)
        self.original_dfs[i] = df
        return df

    def gimme_calc_df(self, i):
        if i in self.calc_dfs:
            return self.calc_dfs[i]
        folder = self.folder
        calc_df_filename = f'{folder}data/oilgathering_fit/calc_df_{i}.csv'
        try:
            calc_df = pd.read_csv(calc_df_filename)
            self.calc_dfs[i] = calc_df
            return calc_df
        except:
            pass
        original_df = self.gimme_original_df(i)
        calc_df = make_calc_df(original_df, self.folder + 'CommonData/')
        calc_df.to_csv(calc_df_filename)
        self.calc_dfs[i] = calc_df
        return calc_df


    def gimme_graph(self, i):
        if i in self.graphs:
            return self.graphs[i]
        df = self.gimme_original_df(i)
        calc_df = self.gimme_calc_df(i)
        G, _ = make_oilpipe_schema_from_OT_dataset(df, self.folder + 'CommonData/', calc_df)
        self.graphs[i] = G
        return self.graphs[i]


    def gimme_solver(self, i):
        if i in self.solvers:
            return self.solvers[i]

        G = self.gimme_graph(i)
        solver = HE2_Solver(G)
        self.solvers[i] = solver
        return solver
    
    def grab_results_to_one_dataframe(self):
        p_rez = dict()
        q_rez = dict()
        for i in range(self.N):
            G = self.gimme_graph(i)
            for n in G.nodes:
                obj = G.nodes[n]['obj']
                grab_q = isinstance(obj, vrtx.HE2_Boundary_Vertex)
                grab_q |= isinstance(obj, vrtx.HE2_Source_Vertex)
                have_to_grab = grab_q
                have_to_grab |= 'pump' in n
                have_to_grab |= 'wellhead' in n
                if not have_to_grab:
                    continue
                res = obj.result
                if grab_q:
                    q_lst = q_rez.get(n, [])
                    q_lst += [res['Q']]
                    q_rez[n] = q_lst

                p_lst = p_rez.get(n, [])
                p_lst += [res['P_bar']]
                p_rez[n] = p_lst
        for n in q_rez:
            print(n, np.round(np.array(q_rez[n]), 3))
        print('---------------------------------------------------')
        for n in p_rez:
            print(n, np.round(np.array(p_rez[n]), 3))
        pass



    def solve_em_all(self):
        for i in range(self.N):
            solver = self.gimme_solver(i)
            solver.solve(threshold=0.25)
            print(i, solver.op_result.success, solver.op_result.fun)


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
        df.loc[df.node_id_start==1750023893, 'start_P'] = np.nan
        df.loc[df.node_id_start==1750024074, 'start_P'] = np.nan
        df.loc[df.node_id_start==1750040926, 'start_P'] = np.nan
        df.loc[df.node_id_start==1750026636, 'start_P'] = np.nan
        df.loc[df.node_id_start==1750037686, 'start_P'] = np.nan
        df.loc[df.node_id_start==1750028316, 'start_P'] = np.nan
        df.loc[df.node_id_start==1750037786, 'start_P'] = np.nan
        self.input_df = df

        self.row_count = self.input_df.shape[0]
        N = self.row_count
        self.d_weights = np.ones(N)
        self.r_weights = np.ones(N)
        self.d_min = np.zeros(N)+1e-5
        self.d_max = np.ones(N)*1.2
        self.r_min = np.zeros(N)
        self.r_max = np.ones(N)*100
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
        if (rez < self.best_y) or (self.it%25==0):
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

        bnds = list(zip(list(np.ones(self.row_count)*0.35), list(np.ones(self.row_count)*1.2)))

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
        bnds = list(zip(list(np.ones(self.row_count)*0.35), list(np.ones(self.row_count)*1.0)))

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
        bnds = list(zip(list(np.ones(self.row_count)*0.35), list(np.ones(self.row_count)*1.0)))
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
        bnds = list(zip(list(np.ones(self.row_count)*0.35), list(np.ones(self.row_count)*1.1)))
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
        x0 = np.ones(self.row_count)*0.9

        op_result = None
        for i in range(20):
            op_result = scop.dual_annealing(target, bounds=bnds, x0=x0, maxfun=1000)
            x0 = op_result.x
        return op_result


    def fit(self):
        fits = {0:self.fit_v0, 1:self.fit_v1, 2:self.fit_v2, 3:self.fit_v3, 4:self.fit_v4, 5:self.fit_v5}
        meth = fits.get(self.fit_version, self.fit_v0)
        return meth()


if __name__ == '__main__':
    model = HE2_OilGatheringNetwork_Model("../../")
    model.solve_em_all()
    model.grab_results_to_one_dataframe()