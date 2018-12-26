"""
"""
import os
import pickle
import copy
import numpy as np
from scipy.interpolate import interp1d

import climada.util.plot as u_plot
from climada.util.save import save
from climada.entity import BlackMarble, ImpactFuncSet, IFTropCyclone
from climada.engine import Impact
from climada.hazard import TropCyclone, Centroids

COLORS = {'SXM': 'red', 'VIR': 'black', 'VGB': 'blue',
          'MAF': 'lightgreen', 'AIA': 'dodgerblue', 'BLM': 'darkolivegreen',
          'TCA': 'lightskyblue', 'ATG': 'orange', 'KNA': 'gold', 'NLD': 'firebrick',
          'British': 'darkblue','French': 'darkgreen', 'Dutch': 'tomato'}
""" Colors of each island """

NAMES = {'SXM': 'St Maarten', 'VIR': 'US Virgin Islands', 'VGB': 'British Virgin Islands',
         'MAF': 'St Martin', 'AIA': 'Anguilla', 'BLM': 'St Barthelemy',
         'TCA': 'Turks And Caicos Islands', 'ATG': 'Antigua And Barbuda',
         'KNA': 'St Kitts And Nevis', 'NLD': 'Saba And St Eustatius'}
""" iso codes and country names of each island """

DAMAGE_IRMA = {'SXM': 2.182E+09, 'VIR': 2.038E+09, 'VGB': 1.426E+09,
               'MAF': 1.282E+09, 'AIA': 7.915E+08, 'BLM': 7.705E+08,
               'TCA': 6.870E+08, 'ATG': 5.376E+08, 'KNA': 2.614E+08,
               'NLD': 9.311E+07}
""" damage irma of each islands """

ENGL = ['AIA', 'VGB', 'TCA']
""" british islands """

FRA = ['BLM', 'MAF']
""" french islands """

HOL = ['SXM', 'NLD']
""" dutch islands """

RET_PERIODS = np.append(np.arange(1, 10), np.arange(10, 3500, 10))
""" Return period values used in the PML """

N_SAMPLES = 100
""" Number of samples used from the distributions for uncertainty analysis """

IF_EXP = ImpactFuncSet()
IF_EM = IFTropCyclone()
IF_EM.set_emanuel_usa()
IF_EXP.add_func(IF_EM)
""" Vulnerability function of Emanuel 2011 """

def ens_if(isl_all, tc_all, data_dir):
    """ Ensembles generation changing impact functions """
    try:
        abs_path = os.path.join(data_dir, 'imp_if_samples.p')
        with open(abs_path, 'rb') as f:
            imp_if = pickle.load(f)
        print('Loaded imp_if_samples:', len(imp_if))
    except FileNotFoundError:
        imp_if = list()
        for sample_i in range(N_SAMPLES):
            v_thresh = np.random.uniform(0.9, 1.1)*25.7
            v_half = np.random.uniform(0.9, 1.1)*74.7
            scale = 1.0
            if_ens = IFTropCyclone()
            if_ens.set_emanuel_usa(v_thresh=v_thresh, v_half=v_half, scale=scale)
            if_all = ImpactFuncSet()
            if_all.add_func(if_ens)

            imp = Impact()
            imp.calc(isl_all, if_all, tc_all)
            imp_if.append(imp)
        save(os.path.join(data_dir, 'imp_if_samples.p'), imp_if)

    return imp_if

def ens_tc(hist_tracks, isl_all, data_dir):
    """ Ensembles generation changing hurricane tracks """
    try:
        abs_path = os.path.join(data_dir, 'imp_tr_samples.p')
        with open(abs_path, 'rb') as f:
            imp_tr = pickle.load(f)
        print('Loaded imp_tr_samples:', len(imp_tr))
    except FileNotFoundError:
        centr_tot = Centroids()
        centr_tot.coord = isl_all.coord
        centr_tot.id = np.arange(centr_tot.lat.size) + 1

        sample_seed = np.arange(N_SAMPLES)+100
        imp_tr = list()
        for seed in sample_seed:
            imp_tr.append(sample_tracks(seed, hist_tracks, isl_all, centr_tot))
        save(os.path.join(data_dir, 'imp_tr_samples.p'), imp_tr)

    return imp_tr

def sample_tracks(seed, hist_tracks, isl_all, centr):
    """ Sampling TC tracks impact"""
    sel_ibtracs = copy.deepcopy(hist_tracks)
    sel_ibtracs.calc_random_walk(ens_size=49, seed=seed)
    print('num tracks hist+syn:', sel_ibtracs.size)

    sel_ibtracs.equal_timestep(1, False)

    tc = TropCyclone()
    tc.set_from_tracks(sel_ibtracs, centr)

    imp = Impact()
    imp.calc(isl_all, IF_EXP, tc)
    return imp

def append_tc_centroids(tc_list):
    """ Append tropical cyclone centroids """
    tc_all = copy.deepcopy(tc_list[0])
    for tc in tc_list[1:]:
        cen_self, cen_haz = tc_all._append_haz_cent(tc.centroids, True)

        tc_all.intensity[:, cen_self] = tc.intensity[:, cen_haz]
        tc_all.fraction[:, cen_self] = tc.fraction[:, cen_haz]
        tc_all.intensity = tc_all.intensity.tocsr()
        tc_all.fraction = tc_all.fraction.tocsr()

        tc_all.check()
    return tc_all

def joint_impact(isl_all, tc_list):
    """ Compute impact of tropical cyclones with different centroids """
    tc_all = append_tc_centroids(tc_list)
    imp_all = Impact()
    imp_all.calc(isl_all, IF_EXP, tc_all)
    efc_all = imp_all.calc_freq_curve()

    return tc_all, imp_all, efc_all

def irma_efc_isl(graph, iso_order, imp_dict):
    """ Irma's imüact and exceedance frequency curve per island """
    for iso_isl in iso_order:
        imp_isl = imp_dict[iso_isl]
        efc = imp_isl.calc_freq_curve()
        irma_id = imp_isl.event_name.index('2017242N16333')
        dam_irma_all = imp_isl.at_event[irma_id]
        dam_irma = np.array([dam_irma_all])
        f_rp = interp1d(efc.impact, efc.return_per)
        rp_irma = np.array([f_rp(dam_irma_all)])

        graph.add_curve(efc.return_per, efc.impact, '-', color=COLORS[iso_isl], label=NAMES[iso_isl])
        graph.axs[graph.curr_ax].scatter(rp_irma, dam_irma, c=COLORS[iso_isl])

    return graph

def irma_efc_mainland(graph, label, expo_dict, tc_dict):
    """ Irma's imüact and exceedance frequency curve per mainland"""
    if label == 'British':
        all_str = ENGL
    elif label == 'French':
        all_str = FRA
    else:
        all_str = HOL

    isl_all = BlackMarble()
    for cntry_iso in all_str:
        isl_all.append(expo_dict[cntry_iso])

    tc_list = []
    for cntry_iso in all_str:
        tc_list.append(tc_dict[cntry_iso])

    tc_all, imp_all, efc_all = joint_impact(isl_all, tc_list)
    efc = imp_all.calc_freq_curve()
    irma_id = tc_all.get_event_id('2017242N16333')[0] - 1
    dam_irma_all = imp_all.at_event[irma_id]
    dam_irma = np.array([dam_irma_all])

    f_rp = interp1d(efc_all.impact, efc_all.return_per)
    rp_irma = np.array([f_rp(dam_irma_all)])

    graph.add_curve(efc.return_per, efc.impact, '--', color=COLORS[label], label=label)
    graph.axs[graph.curr_ax].scatter(rp_irma, dam_irma, c=COLORS[label])
    graph.axs[graph.curr_ax].set_xscale('log')
    graph.axs[graph.curr_ax].set_yscale('log')

    return graph

def filter_orig(tc_all, hist=True):
    """ Build TropCyclone instance with noly historical or synthetic events. """
    tc_filt = TropCyclone()
    tc_filt.tag = tc_all.tag
    tc_filt.units = tc_all.units
    tc_filt.centroids = tc_all.centroids

    if hist:
        filt_idx = np.argwhere(tc_all.orig).squeeze()
    else:
        filt_idx = np.argwhere(np.logical_not(tc_all.orig)).squeeze()
    tc_filt.event_id = tc_all.event_id[filt_idx]
    tc_filt.date = tc_all.date[filt_idx]
    tc_filt.orig = tc_all.orig[filt_idx]
    tc_filt.frequency = tc_all.frequency[filt_idx]
    tc_filt.intensity = tc_all.intensity[filt_idx, :]
    tc_filt.fraction = tc_all.fraction[filt_idx, :]

    for ev_idx in filt_idx:
        tc_filt.event_name.append(tc_all.event_name[ev_idx])

    return tc_filt

def prob_prob_plot(graph, obs, hist_dist, sort_est, sort_est_dist):
    """ probability-probability plot """
    graph.add_subplot('Empirical', 'Model')
    hist_est = np.zeros(hist_dist.shape)
    for val_idx, val in enumerate(obs):
        try:
            hist_est[val_idx] = sort_est_dist[np.argwhere(sort_est - val > 0).reshape(-1)[0]]
        except IndexError:
            hist_est[val_idx] = 1.0

    graph.add_curve(hist_dist, hist_est, 'ok', mfc='none', label='historical events')
    min_val = min(np.min(hist_dist), np.min(hist_est))
    max_val = max(np.max(hist_dist), np.max(hist_est))
    graph.add_curve([min_val, max_val], [min_val, max_val], 'k')

def return_levels_tr_plot(graph, hist_tracks, isl_all, data_dir):
    """ Return levels with uncertainties according to sample of tc with
    different tracks"""
    graph.add_subplot('Return Period', 'Impact')

    sam_rp = dict()
    for rp in RET_PERIODS:
        sam_rp[rp] = list()

    imp_tr_sam = ens_tc(hist_tracks, isl_all, data_dir)
    for i_imp, imp in enumerate(imp_tr_sam):
        efc_ens = imp.calc_freq_curve(RET_PERIODS)
        for (rp, rl) in zip(efc_ens.return_per, efc_ens.impact):
            sam_rp[rp].append(rl)

    mean_rl = np.zeros(len(RET_PERIODS))
    quant_25 = np.zeros(len(RET_PERIODS))
    quant_975 = np.zeros(len(RET_PERIODS))
    for idx, (rp, rl) in enumerate(sam_rp.items()):
        rl_arr = np.array(rl)
        mean_rl[idx] = np.mean(rl_arr)
        quant_25[idx] = np.quantile(rl_arr, 2.5/100)
        quant_975[idx] = np.quantile(rl_arr, 97.5/100)

    graph.add_curve(RET_PERIODS, mean_rl, 'k')
    graph.add_curve(RET_PERIODS, quant_25, 'k--', label='95% CI hazard')
    graph.add_curve(RET_PERIODS, quant_975, 'k--')
    graph.axs[graph.curr_ax].fill_between(RET_PERIODS, quant_25, quant_975)
    graph.axs[graph.curr_ax].set_xscale('log')
    graph.axs[graph.curr_ax].set_xlim(1, 3500)
    graph.axs[graph.curr_ax].set_ylim(0, 2.15e10)
    graph.axs[graph.curr_ax].legend(loc=2)

def return_levels_if_plot(graph, isl_all, tc_all, data_dir):
    """ Return levels with uncertainties according to sample of of vul
    functions with different shapes """
    graph.add_subplot('Return Period', 'Impact')

    sam_rp = dict()
    for rp in RET_PERIODS:
        sam_rp[rp] = list()

    imp_if_sam = ens_if(isl_all, tc_all, data_dir)
    for i_imp, imp in enumerate(imp_if_sam):
        efc_ens = imp.calc_freq_curve(RET_PERIODS)
        for (rp, rl) in zip(efc_ens.return_per, efc_ens.impact):
            sam_rp[rp].append(rl)

    mean_rl = np.zeros(len(RET_PERIODS))
    quant_25 = np.zeros(len(RET_PERIODS))
    quant_975 = np.zeros(len(RET_PERIODS))
    for idx, (rp, rl) in enumerate(sam_rp.items()):
        rl_arr = np.array(rl)
        mean_rl[idx] = np.mean(rl_arr)
        quant_25[idx] = np.quantile(rl_arr, 2.5/100)
        quant_975[idx] = np.quantile(rl_arr, 97.5/100)

    graph.add_curve(RET_PERIODS, mean_rl, 'k')
    graph.add_curve(RET_PERIODS, quant_25, 'k--', label='95% CI impact func')
    graph.add_curve(RET_PERIODS, quant_975, 'k--')
    graph.axs[graph.curr_ax].fill_between(RET_PERIODS, quant_25, quant_975)
    graph.axs[graph.curr_ax].set_xscale('log')
    graph.axs[graph.curr_ax].set_xlim(1, 3500)
    graph.axs[graph.curr_ax].set_ylim(0, 2.15e10)
    graph.axs[graph.curr_ax].legend(loc=2)

def quality_ensembles(hist_tracks, tc_dict, expo_dict, data_dir, cond_zero=False):
    """ Graph with 3 plots """
    # join TCs and islands and compute its impact
    isl_all = BlackMarble()
    for cntry_iso, cntry_val in expo_dict.items():
        isl_all.append(cntry_val)
    print('size exposures: ', isl_all.size)

    tc_all, _, efc_all = joint_impact(isl_all, list(tc_dict.values()))
    sort_est = efc_all.impact[::-1]
    if cond_zero:
        sort_est = sort_est[sort_est > 0]
    _, uni_idx, uni_inv, uni_cnt = np.unique(sort_est, return_index=True,
                                             return_inverse=True, return_counts=True)
    sort_est_dist = np.arange(sort_est.size)/sort_est.size
    sort_est_dist = sort_est_dist[(uni_idx + uni_cnt -1)[uni_inv]]
    print('size sort_est_dist:', sort_est_dist.size)

    # impact distributions of historical events
    tc_hist = filter_orig(tc_all, hist=True)
    imp_hist = Impact()
    imp_hist.calc(isl_all, IF_EXP, tc_hist)
    efc_hist = imp_hist.calc_freq_curve()

    obs = efc_hist.impact[::-1]
    if cond_zero:
        obs = obs[obs > 0]
    _, uni_idx, uni_inv, uni_cnt = np.unique(obs, return_index=True,
                                             return_inverse=True, return_counts=True)
    hist_dist = np.arange(obs.size)/obs.size
    hist_dist = hist_dist[(uni_idx + uni_cnt -1)[uni_inv]]
    print('size hist_dist:', hist_dist.size)

    # plots
    graph = u_plot.Graph2D(num_row=1, num_col=3)

    # probability plot
    prob_prob_plot(graph, obs, hist_dist, sort_est, sort_est_dist)

    # return levels with uncertainty plot
    return_levels_tr_plot(graph, hist_tracks, isl_all, data_dir)

    # return levels with uncertainty plot
    return_levels_if_plot(graph, isl_all, tc_all, data_dir)
    return graph

def graphs_exceedance(expo_dict, tc_dict, imp_dict):
    """ FIGURE EXCEEEDANCE FREQUENCY CURVES """
    graph = u_plot.Graph2D(num_subplots=4)
    graph.fig.set_size_inches(15, 9)
    graph.add_subplot('Return period (year)', 'Impact (USD)')
    irma_efc_isl(graph, ['VIR', 'KNA', 'ATG'], imp_dict)

    graph.add_subplot('Return period (year)', 'Impact (USD)')
    irma_efc_mainland(graph, 'British', expo_dict, tc_dict)
    irma_efc_isl(graph, ENGL, imp_dict)

    graph.add_subplot('Return period (year)', 'Impact (USD)')
    irma_efc_mainland(graph, 'French', expo_dict, tc_dict)
    irma_efc_isl(graph, FRA, imp_dict)

    graph.add_subplot('Return period (year)', 'Impact (USD)')
    irma_efc_mainland(graph, 'Dutch', expo_dict, tc_dict)
    irma_efc_isl(graph, HOL, imp_dict)

    for ax in graph.axs:
        ax.set_xscale('log')
        ax.set_yscale('log')
        ax.set_xlim(1, 3500)
        ax.set_ylim(1e4, 1e10)
        ax.legend(loc=4)

    graph.fig.subplots_adjust(hspace=0.3, wspace=0.2)
    return graph

def fig06(data_dir, fig_dir):
    """ Generate fig06. """
    abs_path = os.path.join(data_dir, 'tc_isl.p')
    with open(abs_path, 'rb') as f:
        tc_dict = pickle.load(f)
    print('Loaded tc_isl:', len(tc_dict))

    abs_path = os.path.join(data_dir, 'imp_isl.p')
    with open(abs_path, 'rb') as f:
        imp_dict = pickle.load(f)
    print('Loaded imp_isl:', len(imp_dict))

    abs_path = os.path.join(data_dir, 'exp_irma.p')
    with open(abs_path, 'rb') as f:
        expo_dict = pickle.load(f)
    print('Loaded exp_irma:', len(expo_dict))

    graph = graphs_exceedance(expo_dict, tc_dict, imp_dict)
    if fig_dir:
        graph.fig.savefig(os.path.join(fig_dir, 'fig06.png'), format='png', bbox_inches='tight')

def fig07(data_dir, fig_dir):
    """ Generate fig07. """
    abs_path = os.path.join(data_dir, 'tc_isl.p')
    with open(abs_path, 'rb') as f:
        tc_dict = pickle.load(f)
    print('Loaded tc_isl:', len(tc_dict))

    abs_path = os.path.join(data_dir, 'sel_hist.p')
    with open(abs_path, 'rb') as f:
        hist_tracks = pickle.load(f)
    print('Loaded sel_hist:', hist_tracks.size)

    abs_path = os.path.join(data_dir, 'exp_irma.p')
    with open(abs_path, 'rb') as f:
        expo_dict = pickle.load(f)
    print('Loaded exp_irma:', len(expo_dict))

    graph = quality_ensembles(hist_tracks, tc_dict, expo_dict, data_dir)

    if fig_dir:
        graph.fig.savefig(os.path.join(fig_dir, 'fig07.png'), format='png', bbox_inches='tight')
