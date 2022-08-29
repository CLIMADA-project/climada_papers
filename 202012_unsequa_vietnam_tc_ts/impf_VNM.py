import numpy as np
from climada.entity import IFTropCyclone, ImpactFuncSet, ImpactFunc


def impf(storm_threshold = 50, surge_threshold = 1):


    # Impact function wind-people

    ImpfSetPeople = ImpactFuncSet()

    ## Impact function surge-people

    if_ts_people = ImpactFunc()
    if_ts_people.haz_type = 'TS'
    if_ts_people.id = 3
    if_ts_people.intensity_unit = 'm'

    if_ts_people.intensity = np.arange(0, 10, 0.05)
    if_ts_people.paa = np.ones(if_ts_people.intensity.size)
    if_ts_people.mdd = np.zeros(if_ts_people.intensity.size)
    for idx in range(if_ts_people.mdd.size):
        if if_ts_people.intensity[idx] >= surge_threshold:
            if_ts_people.mdd[idx] = 1

    if_ts_people.check()

    ImpfSetPeople.append(if_ts_people)

    return ImpfSetPeople




