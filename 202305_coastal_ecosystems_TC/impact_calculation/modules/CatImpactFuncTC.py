import numpy as np
from climada.entity import ImpactFunc, ImpactFuncSet


def cat_impf_set(haz_type):
    """
    creates the impact function set based on the hazard categories and assigns 1 as their impact
    an impact function specific for the chosen hazard
    :param haz_type: abbreviation of hazard type
    :return:
    an impact function set for all hazard categories
    """
    cat_impfset = ImpactFuncSet()
    haz_cat = {
        1: [33, 42],
        2: [43, 49],
        3: [50, 58],
        4: [58, 70],
        5: [70, 200]
    }
    for cat, [i_min, i_max] in haz_cat.items():
        impf_tmp = cat_impf([i_min, i_max], haz_type)
        impf_tmp.paa[i_min: i_max] = 1
        impf_tmp.id = cat
        cat_impfset.append(impf_tmp)
    return cat_impfset


def cat_impf(intensity, haz_type):
    """
    creates a 2-step impact function, with the impact being 1 between the two thresholds
    and 0 outside of the thresholds
    :type intensity: object
    :param intensity: Four windspeeds to define the step
    :param haz_type: Abbreviation of hazard type
    :return:
    an impact function for a specific hazard category
    """
    impf = ImpactFunc()
    x1, x2 = intensity
    impf.intensity = np.array([0, x1, x1, x2, x2, x2 + 1])
    impf.mdd = np.array([0, 0, 1, 1, 0, 0])
    impf.paa = np.ones(len(impf.mdd))
    impf.haz_type = haz_type
    return impf
