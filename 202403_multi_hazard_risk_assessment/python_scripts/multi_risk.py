import copy
import numpy as np
import pandas as pd
import scipy as sp
from climada.engine import Impact
from climada.util import yearsets, coordinates
import matplotlib.pyplot as plt
from climada.util.lines_polys_handler import set_imp_mat  
from sklearn.utils import shuffle


def get_indices(impact, n_samples, list_string):
    """Return the indices of the event names ordered based on the given lists of strings."""

    n_samples = int(n_samples)

    indices = {str1: np.unique([i for i, elem in enumerate(impact.event_name) if
                      (str(str1) in elem)]) for str1 in list_string}

    indices = [np.random.choice(indices[str1]) for str1 in list_string for
               n in range(n_samples)]

    return indices


def order_events_by_indices(impact, indices):
    """
    Order event names based on given strings contained in the event names.

    Parameters
    ----------
    impact: Impact
        with event_name based on the given strings
    n_events: Int
        Number of events in the output. Default: 1
    list_string : list
        A list of string based on which to order the events.
        For example climate models ['miroc5','ipsl-cm5a-lr','gfdl-esm2m','hadgem2-es']
        default is None


    Raises
    ------
    AttributeError
        If no list is providing

    Returns
    -------
    impact : Impact
        Impact yearset.

    """
    impact_ordered = Impact()
    impact_ordered.imp_mat = impact.imp_mat[indices]
    impact_ordered.event_name = [impact.event_name[index] for index in indices]
    impact_ordered.event_id = np.arange(len(impact_ordered.event_name))
    frequency = impact.frequency[indices]
    impact_ordered.frequency = frequency*(len(impact.event_id)/len(impact_ordered.event_id))
    impact_ordered.at_event = impact.at_event[indices]
    impact_ordered.aai_agg = np.median(impact_ordered.at_event)
    impact_ordered.coord_exp = impact.coord_exp
    impact_ordered.date = impact.date[indices]
    return impact_ordered


def aggregate_impact_from_event_name(impact, how='sum', exp=None):
    """
    Aggregate the impact per year to make yearsets. Maximum impact per year
    at each exposure point is exposure value if exp is not None.

    Parameters
    ----------
    impact : Impact
        Impact with an impact matrix and events with dates per year
    how : How to aggregate impacts, options are 'sum' or 'max'
    exp : Exposure
        Exposure of Impact to cap the impact value at the value of the exposure
    Raises
    ------
    AttributeError
        If impact matrix is empty.

    Returns
    -------
    impact : Impact
        Impact yearset.

    """

    imp = copy.deepcopy(impact)
    if how == 'sum':
        mask = [np.ma.make_mask(np.array(imp.event_name) == event).astype(int)
                for event in np.unique(imp.event_name)]
        mask_matrix = sp.sparse.csr_matrix(mask)
        imp_mat = mask_matrix.dot(imp.imp_mat)

    elif how == 'max':
        imp_mat = sp.sparse.csr_matrix(sp.sparse.vstack(
        [imp.imp_mat[(np.array(imp.event_name) == event).astype(bool)].max(axis=0)
         for event in np.unique(imp.event_name)]))

    if exp is not None:
        m1 = imp.imp_mat.data
        m2 = exp.gdf.value[imp.imp_mat.nonzero()[1]]
        mat = sp.sparse.csr_matrix((np.minimum(m1, m2), imp.imp_mat.indices, imp.imp_mat.indptr))
        mat.eliminate_zeros()
        imp.imp_mat = mat

    imp.frequency = np.ones(imp_mat.shape[0])/imp_mat.shape[0]
    imp = set_imp_mat(imp, imp_mat)
    imp.date = np.arange(1, imp.imp_mat.shape[0] + 1)
    imp.event_id = np.arange(1, imp.imp_mat.shape[0] + 1)
    imp.event_name = np.unique(imp.event_name)
    return imp


def combine_yearsets(impact_list, how='sum', occur_together=False, exposures=None):
    """
    Parameters
    ----------
    impact_list : list  or dict of impacts
    how : how to combine the impacts, options are 'sum', 'max' or 'min'
    exp : If the exposures are given, the impacts are caped at their value

    Returns
    -------
    imp : Impact
        Combined impact
    """

    if type(impact_list) is dict:
        impact_list = list(impact_list.values())
    imp0 = copy.deepcopy(impact_list[0])

    if how == 'sum':
        imp_mat = impact_list[0].imp_mat + impact_list[1].imp_mat # to do: adapt for more than two hazards
    elif how == 'min':
        imp_mat_min = imp0.imp_mat
        for imp in impact_list[1:]:
            imp_mat_min = imp_mat_min.minimum(imp.imp_mat)
        imp_mat = imp_mat_min

    elif how == 'max':
        imp_mat_max = imp0.imp_mat
        for imp in impact_list[1:]:
            imp_mat_max = imp_mat_max.maximum(imp.imp_mat)
        imp_mat = imp_mat_max
    else:
        raise ValueError(f"'{how}' is not a valid method. The implemented methods are sum, max or min")

    if exposures is not None:
        m1 = imp_mat.data
        m2 = exposures.gdf.value[imp_mat.nonzero()[1]]
        imp_mat = sp.sparse.csr_matrix((np.minimum(m1, m2), imp_mat.indices, imp_mat.indptr))
        imp_mat.eliminate_zeros()

    if occur_together:
        mask_list = [np.abs(impact.imp_mat.A[imp_mat.nonzero()]) == 0 for impact in impact_list]
        for mask in mask_list:
            imp_mat.data[mask] = 0
        imp_mat.eliminate_zeros()
    imp0.frequency = np.ones(len(imp0.event_id))/ len(imp0.event_id)
    imp0 = set_imp_mat(imp0, imp_mat)
    return imp0


def sample_events(impact, years, lam=1):
    #create sampling vector
    events_per_year = yearsets.sample_from_poisson(len(years), lam)
    sampling_vect = yearsets.sample_events(events_per_year, impact.frequency)
    impact_sample = yearsets.impact_from_sample(impact, years, sampling_vect)
    return impact_sample


def make_yearset_samples(impact, years, exposures=None, n_samples=1):
    [make_yearset(impact, years, exposures) for sample in range(n_samples)]


def make_yearset(impact_dict, n_years, aggregation_method=None, exposures=None):
    """
    Parameters
    ----------
    impact_dict : dict of Impact
        dictionary of impacts
    n_years : number of years
    exposures : If the exposures are given, the impacts are caped at their value

    Returns
    -------
    imp : Impact
        Combined impact
    """
    years = np.arange(1, n_years+1)
    yearset_dict = {}
    for impact in impact_dict:
        lam = np.sum(impact_dict[impact].frequency)
        lam = np.round(lam, 10)
        events_per_year = yearsets.sample_from_poisson(len(years), lam)
        sampling_vect = yearsets.sample_events(events_per_year, impact_dict[impact].frequency)
        yearset = yearsets.impact_from_sample(impact_dict[impact], years, sampling_vect)
        if yearset.imp_mat.shape[0]>len(years):
            yearset = yearsets.aggregate_impact_to_year(yearset, exp=exposures[impact], how=aggregation_method[impact])
        yearset_dict[impact] = yearset
    return yearset_dict


def order_by_event_name(impact_dict, n_samples, list_string):
    """
    Parameters
    ----------
    impact_dict : dict of Impact objects
        dictionary or list of impacts.
    n_samples : number of samples, the total number of events will be of len(list_samples_string)*n
    list_string : List of string to look for in the event names

    Returns
    -------
    ordered_impact_dict : dict of Impact objects
        Ordered impact dictionary
    """
    ordered_impact_dict = {}
    event_names_list = []
    indices_list = []
    for i, hazard in enumerate(impact_dict):
        if i ==0:
            indices = get_indices(impact_dict[hazard], n_samples,
                                  list_string=list_string)
        for n in range(i):
            if list(impact_dict[hazard].event_name) == list(event_names_list[n]):
                indices = indices_list[n]
            else:
                indices = get_indices(impact_dict[hazard], n_samples,
                                      list_string=list_string)
        event_names_list.append(impact_dict[hazard].event_name)

        indices_list.append(indices)
        ordered_impact_dict[hazard] = order_events_by_indices(impact_dict[hazard], indices_list[i])
        ordered_impact_dict[hazard].event_name = event_names_list
    return ordered_impact_dict

def corr_btw_hazards():
    return

def corr_btw_impacts(impact_dict, temporal=True,spatial=False):
    """
    Parameters
    ----------
    impact_dict : dict of Impacts
    temporal : bool, rather to consider the temporal dimension
    spatial : bool, rather to consider the spatial dimension

    Returns
    -------
    corr_df : pd.DataFrame
    """
    if temporal is True and spatial is False:
        df = pd.DataFrame.from_dict({hazard: impact_dict[hazard].at_event for hazard in impact_dict})
    if spatial is True and temporal is False:
        df = pd.DataFrame.from_dict({hazard: impact_dict[hazard].eai_exp for hazard in impact_dict})
    if spatial is True and temporal is True:
        df = pd.DataFrame.from_dict({hazard: np.array(impact_dict[hazard].imp_mat.todense().flatten())[0]
                                     for hazard in impact_dict})
    return df.corr()


def make_country_matrix(impact, countries):
    """
    Parameters
    ----------
    impact : Impact
    countries : list of countries iso3 code

    Returns
    -------
    corr_df : dict of imp_mat, keys are country iso3 codes
    """
    lat = np.array([lat for lat,lon in impact.coord_exp])
    lon = np.array([lon for lat,lon in impact.coord_exp])
    countries_num = coordinates.get_country_code(lat, lon)
    countries_num = np.array([format(num, '03d') for num in countries_num])
    country_matrices = {country: np.array([impact.imp_mat[:, countries_num == country].sum(axis=1)])
                        for country in countries}
    return country_matrices


def normalize_array(array):
    if np.sum(array)==0:
        return array
    return (array - np.min(array)) 


def make_fq_list(impact, factor):
    fq_list = []
    for n in range(10000):
        imp_sample = order_events_by_indices(impact, shuffle(np.arange(len(impact.event_name)))[0:500])
        fq = imp_sample.calc_freq_curve(np.arange(0,100))
        fq_list.append(fq.impact/factor)
    return fq_list, fq.return_per

def plot_return_period_samples(fq_list, return_per, ax, color, label, range=True, linestyle="solid", linewidth=1.5):
    ax.plot(return_per,np.median(fq_list,axis=0), color=color, label=label, linestyle=linestyle,linewidth=linewidth)
    ax.legend()
    if range:
        ax.fill_between(return_per,np.median(fq_list,axis=0),np.percentile(fq_list,q=95,axis=0), color=color, alpha=0.3)
        ax.fill_between(return_per,np.percentile(fq_list,q=5,axis=0), np.median(fq_list,axis=0),color=color, alpha=0.3)
    return ax


