import numpy as np
import matplotlib.pyplot as plt
from scipy.interpolate import interp1d

from climada.engine import Impact

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

def irma_efc_isl(ax, iso_order, imp_dict):
    """ Irma's imüact and exceedance frequency curve per island """
    for iso_isl in iso_order:
        imp_isl = imp_dict[iso_isl]
        efc = imp_isl.calc_freq_curve()
        irma_id = imp_isl.event_name.index('2017242N16333')
        dam_irma_all = imp_isl.at_event[irma_id]
        dam_irma = np.array([dam_irma_all])
        f_rp = interp1d(efc.impact, efc.return_per)
        rp_irma = np.array([f_rp(dam_irma_all)])

        ax.plot(efc.return_per, efc.impact, '-', color=COLORS[iso_isl], label=NAMES[iso_isl])
        ax.scatter(rp_irma, dam_irma, c=COLORS[iso_isl])

def irma_efc_mainland(ax, label, imp_dict):
    """ Irma's imüact and exceedance frequency curve per mainland"""
    if label == 'British':
        all_str = ENGL
    elif label == 'French':
        all_str = FRA
    else:
        all_str = HOL

    imp_ev_info = imp_dict[ENGL[0]]

    # add at_event: events damages for every island
    imp_all = Impact()
    imp_all.frequency = imp_ev_info.frequency
    imp_all.at_event = np.zeros((imp_dict[ENGL[0]].at_event.size,))
    for ev in range(imp_all.at_event.size):
        for isl in all_str:
            imp_all.at_event[ev] += imp_dict[isl].at_event[ev]

    efc_all = imp_all.calc_freq_curve()

    irma_id = imp_ev_info.event_name.index('2017242N16333')
    dam_irma_all = imp_all.at_event[irma_id]
    dam_irma = np.array([dam_irma_all])

    f_rp = interp1d(efc_all.impact, efc_all.return_per)
    rp_irma = np.array([f_rp(dam_irma_all)])

    ax.plot(efc_all.return_per, efc_all.impact, '--', color=COLORS[label], label=label)
    ax.scatter(rp_irma, dam_irma, c=COLORS[label])
    ax.set_xscale('log')
    ax.set_yscale('log')

def graphs_exceedance(imp_dict):
    """ FIGURE EXCEEEDANCE FREQUENCY CURVES """
    fig, ax = plt.subplots(2, 2)
    fig.set_size_inches(15, 9)
    ax[0, 0].set_xlabel('Return period (year)')
    ax[0, 0].set_ylabel('Impact (USD)')
    irma_efc_isl(ax[0, 0], ['VIR', 'KNA', 'ATG'], imp_dict)
    ax[0, 0].set_xscale('log')
    ax[0, 0].set_yscale('log')

    ax[0, 1].set_xlabel('Return period (year)')
    ax[0, 1].set_ylabel('Impact (USD)')
    irma_efc_mainland(ax[0, 1], 'British', imp_dict)
    irma_efc_isl(ax[0, 1], ENGL, imp_dict)

    ax[1, 0].set_xlabel('Return period (year)')
    ax[1, 0].set_ylabel('Impact (USD)')
    irma_efc_mainland(ax[1, 0], 'French', imp_dict)
    irma_efc_isl(ax[1, 0], FRA, imp_dict)

    ax[1, 1].set_xlabel('Return period (year)')
    ax[1, 1].set_ylabel('Impact (USD)')
    irma_efc_mainland(ax[1, 1], 'Dutch', imp_dict)
    irma_efc_isl(ax[1, 1], HOL, imp_dict)

    for axis in ax.flatten():
        axis.set_xlim(1, 3500)
        axis.set_ylim(1e4, 1e10)
        axis.legend(loc=4)

    fig.subplots_adjust(hspace=0.3, wspace=0.2)
    return fig, ax
