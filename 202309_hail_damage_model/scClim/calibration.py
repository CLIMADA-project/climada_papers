# -*- coding: utf-8 -*-
"""
Empirical calibration functions for spatially explicit vulnerability function calibration


"""
import xarray as xr
import numpy as np
import pandas as pd
import sys #os, argparse
from scipy import sparse
import climada
from climada.hazard import Hazard
from climada import CONFIG
sys.path.append(str(CONFIG.local_data.func_dir))
import scClim as sc
from scClim.constants import CUT_OFF_DICT, INT_RANGE_DICT
import datetime as dt
import copy
import matplotlib

# %% Helper functions
def log_func(x, a, b, c, d):
    return a / (1.0 + np.exp(-c * (x - d))) + b

def fit_log_curve(xData,yData,p0):
    from scipy.optimize import curve_fit

    fittedParameters, pcov = curve_fit(log_func, xData, yData,p0)
    print('Fitted parameters:', fittedParameters)
    print()

    xrange=np.arange(20,121,1)
    modelPredictions = log_func(xData, *fittedParameters)

    absError = modelPredictions - yData

    SE = np.square(absError) # squared errors
    MSE = np.mean(SE) # mean squared errors
    RMSE = np.sqrt(MSE) # Root Mean Squared Error, RMSE
    Rsquared = 1.0 - (np.var(absError) / np.var(yData))

    modelPredictions_out = log_func(xrange, *fittedParameters)

    print('RMSE:', RMSE)
    print('R-squared:', Rsquared)

    print()

    return xrange, modelPredictions_out, Rsquared, fittedParameters

def filter_exp(exp, yearMin, yearMax):
    """filter exposure to remove buildings that are not built yet

    Args:
        exp (climada.entity.Exposures): exposure to be filtered
        yearMin (int): minimum year to be kept (usually 0, except if only subset
                                                of buildings is to be kept)
        yearMax (max): maximum year to be kept (usually the year of the hazard event)

    Returns:
        exp_out: filtered exposure
    """
    exp_out = exp.copy()
    sel = np.logical_and(exp.gdf.Baujahr >= yearMin,
                         exp.gdf.Baujahr <= yearMax)
    exp_out.gdf = exp.gdf.loc[sel,:]
    if 'centr_HL' in exp_out.gdf.columns:
        exp_out.gdf.centr_HL = exp_out.gdf.centr_HL.astype(int)
    return exp_out
    # def filter_exp_fast(exp, yearMin, yearMax): #only slightly faster!
    # sel = np.logical_and(exp.gdf.Baujahr >= yearMin,
                         # exp.gdf.Baujahr <= yearMax)
    # out_gdf = exp.gdf.where(sel).dropna()
    # if 'centr_HL' in out_gdf.columns:
        # out_gdf.centr_HL = out_gdf.centr_HL.astype(int)
    # return out_gdf

def get_exposed_assets(exposure, variable, exposure_type, variable_2 = None):
    """

    Get number and total value of exposed assets

    Parameters
    ----------
    exposure : climada.entity.exposures.base.Exposures
        exposure for which the damage function is to be calibrated
    variable : str
        hazard variable name
    exposure_type : str
        type of exposure provided. "buildings" or "agriculture"
    variable_2:

    Returns
    -------
    all_count : pandas.Dataframe
        Dataframe with hazard intensities as index and number of exposed assets as column
    all_value : pandas.Dataframe
        Dataframe with hazard intensities as index and the total value of the
        exposed assets as column
    at_centroid : pandas.Dataframae
        Dataframe with centroids as index and colums with the hazard intensity,
        the number of exposed assets (counts), and the total value

    """

    #count number of exposed assets per MESHS
    if exposure_type == 'buildings':
        all_count = exposure.gdf.groupby(variable).VersicherungsID.count().rename('count_all')
        count_centr = exposure.gdf.groupby('centr_HL').VersicherungsID.count().rename('counts')
    elif exposure_type == 'agriculture':
        all_count = exposure.gdf.groupby(variable).n_fields.sum().rename('count_all')
        count_centr = exposure.gdf.groupby('centr_HL').n_fields.sum().rename('counts')

    #get intensity at centroid
    intens_centr = exposure.gdf.groupby('centr_HL')[variable].mean().rename('intensity')

    #compute total value exposed to hail hazard per value hazard variable
    all_value = exposure.gdf.groupby(variable).value.sum().rename('value_all')

    #get total value at centroids
    value_centr = exposure.gdf.groupby('centr_HL').value.sum().rename('value')

    if variable_2 is not None and variable_2 in exposure.gdf.keys():
        intens_2_centr = exposure.gdf.groupby('centr_HL')[variable_2].max().rename('intensity_2')
        #concat intensity values, counts, and values at centroid
        at_centroid=pd.concat([intens_centr,count_centr,value_centr,intens_2_centr], axis = 1)
    else:
        at_centroid=pd.concat([intens_centr,count_centr,value_centr], axis = 1)

    return all_count, all_value, at_centroid

def add_field_area_to_damages(exp_dmg, exposure, field_area_mean):

        """
        for agriculture, assign average field area from exposure to each
        reported damage claim if damage is reported in a region without exposure,
        use Swiss average field area
        """

        #add field area column
        exp_dmg.gdf['field_area']=np.nan

        #loop over all damages and add field area
        for index in exp_dmg.gdf.index:

            #get centroids
            centroid=exp_dmg.gdf.loc[index,'centr_HL']

            #read field area from exposure
            field_area=exposure.gdf.loc[exposure.gdf[exposure.gdf['centr_HL']==centroid].index,'field_area'].values

            if field_area.size == 1:
                exp_dmg.gdf.loc[index,'field_area']=field_area[0]
            elif field_area.size == 0:
                #issue warning if damage is recorded where no exposure is present
                exp_dmg.gdf.loc[index,'field_area']=field_area_mean
                #warnings.warn('Damage reported in region without exposure. Swiss average field size used.')
            else:
                raise ValueError('more than one field area found.')

        return exp_dmg

def get_damaged_assets(damage, variable, date_now, exposure_type,
                       get_PVA = False, variable_2 = None, fraction_insured = 1):

        """ get variables specific to damage claims as a function of hazard variable

        Parameters
        ----------
        damage: climada impact
                damage claims of an event as climada impact object
        variable: str
                hazard variable name
        exposure_type: str
            type of exposure ('agriculture' or 'buildings')

        Returns
        -------
        dmg_count: pandas.Dataframe
                number of damage claims grouped by hazard intensity
        dmg_value: pandas.Dataframe
                sum of damage value grouped by hazard intensity
        affected_value: pandas.Dataframe
                total value affected grouped by hazard intensity
                (i.e. Versicherungssumme or total field area of damaged fields)
        at_centroid: pandas.Dataframe
               counts, value, and value affected by centroid
        """

        #date selection
        date_idx = damage.date == date_now.toordinal()
        dmg_1 = damage.imp_mat[date_idx,:]
        #exposure (centroid) selection
        obj_idx=dmg_1.nonzero()[1]
        #get impact at current date
        imp_now = damage.imp_mat[date_idx,obj_idx].getA1()
        #get hazard intensity at current date
        haz_intensity_now = damage.haz_intensity[date_idx,obj_idx].getA1()
        #get second hazard intensity at current date
        if hasattr(damage,'haz_intensity_2'):
            haz_intensity_2_now = damage.haz_intensity_2[date_idx,obj_idx].getA1()

        #get value affected (needs attribute aff_mat for damage impact object)
        if get_PVA == True:
            val_affected_now = damage.aff_mat[date_idx,obj_idx].getA1()

        #create Dataframe with columns damages, hazard intensity, and centroid
        if hasattr(damage,'haz_intensity_2'):
            dmg_df = pd.DataFrame({'damage': imp_now,
                                    variable: haz_intensity_now,
                                    variable_2: haz_intensity_2_now,
                                    'centr_HL':damage.centr_HL[obj_idx]})
        else:
            dmg_df = pd.DataFrame({'damage':imp_now,
                                    variable:haz_intensity_now,
                                    'centr_HL':damage.centr_HL[obj_idx]})

        #if requested, compute total value affected per intensity value and per centroid
        if get_PVA == True:
            dmg_df['value_affected'] = val_affected_now
            affected_value = dmg_df.groupby(variable).value_affected.sum().rename('value_affected')/fraction_insured
            aff_val_centr = dmg_df.groupby('centr_HL').value_affected.sum().rename('value_affected')/fraction_insured

        else:
            affected_value = None

        # get number of damage claims and damage value per intensity value and per centroid
        dmg_count = dmg_df.groupby(variable).damage.count().rename('count_dmg')/fraction_insured
        dmg_count = dmg_count.round(decimals=0) #only allow integers for count
        dmg_value = dmg_df.groupby(variable).damage.sum().rename('value_dmg')/fraction_insured
        dmg_count_centr = dmg_df.groupby('centr_HL').damage.count().rename('counts')/fraction_insured
        dmg_count_centr = dmg_count_centr.round(decimals=0) #only allow integers for count
        dmg_value_centr = dmg_df.groupby('centr_HL').damage.sum().rename('value')/fraction_insured
        intens_centr = dmg_df.groupby('centr_HL')[variable].max().rename('intensity')

        #get value of second hazard variable at centroids
        if hasattr(damage,'haz_intensity_2'):
            intens_2_centr = dmg_df.groupby('centr_HL')[variable_2].max().rename('intensity_2')
            at_centroid_list=[dmg_count_centr, dmg_value_centr, intens_centr, intens_2_centr]
        else:
            at_centroid_list=[dmg_count_centr, dmg_value_centr, intens_centr]


        #concatenate per centroid data
        #concat counts, damage value, and value affected
        if get_PVA:
            at_centroid_list.append(aff_val_centr)

        at_centroid=pd.concat(at_centroid_list, axis = 1)

        return dmg_count, dmg_value, affected_value, at_centroid

def smooth_data_running_mean(ds, variable, intensity_range, cut_off, window_size):
    """

    Parameters
    ----------
    ds : xr.dataset
        containing calibration data.
    variable : str
        Hazard variable.
    intensity_range : np.array
        Intensity range of hazard variable.
    cut_off : int
        cut off value for the smooth-empirical fit.
    window_size : int
        windowsize for rolling average.

    Returns
    -------
    ds_roll : xr.dataset
        same as ds, but with rolling average values.
    ds_roll_cut : xr.dataset
        same as ds_roll, but with values above cut_off as constant defined by
        the weighted average.

    """

    if intensity_range[0] == 0: #all values except 0, because it should not be used for rolling
        ds_roll = ds.sel({variable:slice(intensity_range[1],intensity_range[-1])})
        ds_roll = ds_roll.rolling(dim={variable:window_size}, min_periods=max(1,int( #add max(1,x) in case windowsize=1 (no rolling)
            np.floor(window_size/2))), center=True).sum()
        ds_roll = xr.concat([ds.sel({variable:0}),ds_roll],dim=variable)
    else: #include all values in rolling average
        ds_roll = ds.sel({variable:slice(intensity_range[0], intensity_range[-1])})
        ds_roll = ds_roll.rolling(dim={variable: window_size}, min_periods=max(1,int( #add max(1,x) in case windowsize=1 (no rolling)
            np.floor(window_size/2))), center=True).sum()

    #set values <1e-8 to zero. Bc of a bug (rounding errors) in xarray.rolling()
    ds_roll = ds_roll.where(ds_roll>1e-8,other=0)

    #create rolling average cut off at a threshold
    ds_roll_cut = ds_roll.copy(deep=True)
    ds_roll_cut.loc[{variable:slice(cut_off,intensity_range.max())}] = ds_roll_cut.loc[{variable:slice(
        cut_off,intensity_range.max())}].sum(dim=variable).expand_dims({variable:intensity_range[intensity_range>=cut_off]},axis=1)

    return ds_roll, ds_roll_cut

def get_unique_dates_from_impacts(damages,dates_modeled_imp=None):
    """
    Get unique dates from a measured impact (observed dates) and modeled impacts (modeled dates)'

    Parameters
    ----------
    damages : climada.engine.impact.Impact
        Impact object containing observed impacts
    dates_modeled_imp : list
        list of dates (dt.datetime) based on modeled data (e.g. a modeled impact)
        The default is None.

    Returns
    -------
    dt_dates : list
        unique list of dates (dt.datetime) that comining observed and modeled dates with impact

    """
    import datetime as dt

    # str_dates = [dt.datetime.fromordinal(d).strftime('%d%m%Y') for d in damages.date]
    dt_dates = [dt.datetime.fromordinal(d) for d in damages.date]

    if dates_modeled_imp is not None:
        dt_dates = np.unique(np.append(dt_dates, dates_modeled_imp))

    #dates must be sorted for correct exposure filtering (with 'Baujahr' column')
    return sorted(dt_dates)

def assign_intensity_to_damages(damages, haz, variable, haz_2 = None, variable_2 = None):
    """


    Parameters
    ----------
    damages : climada.Impact
        damages impact object
    haz : climada.Hazard
        hazard object
    variable : str
        variable name

    Returns
    -------
    damages : climada.Impact
        damages with new attribute haz_intensity

    """
    haz_date_idx = np.array([np.where(haz.date == date)[0][0] for date in damages.date])
    haz_loc_idx = damages.centr_HL
    damages.haz_intensity = haz.intensity[haz_date_idx,:][:,haz_loc_idx]

    if haz_2 is None:
        pass
    else:
        damages.haz_intensity_2 = haz_2.intensity[haz_date_idx,:][:,haz_loc_idx]

    return damages

def at_centroid_from_exp_and_dmg(exp_at_centroid, variable, date_now,
                                 dmg_at_centroid = None, variable_2 = None,
                                 get_PVA = False):
    """
    Parameters
    ----------
    exp_at_centroid : pd.DataFrame
        Dataframe containing exposure data at centroids
    variable : str
        hazard variable
    date_now : datetime
        current date
    dmg_at_centroid : pd.DataFrame, optional
        Dataframe containing damage data at centroids. The default is None.
    variable_2 : str, optional
        secondary hazard variable to be traced along analysis. The default is None.
    get_PVA : boolean, optional
        get PVA and MDD in output, requires a column 'value_affected' in
        dmg_at_centroids. The default is False.

    Returns
    -------
    at_centroid_dict : Dictionary
        Dictionary containing the important values at all centroids
        that have at least a damage or an overlap between
        exposure and hazard (i.e., a damage prediction)

    """

    at_centroid_dict = {}

    if dmg_at_centroid is not None:
        # reindex damages at centroid to exposure indices
        new_df = pd.DataFrame(index=np.unique(np.concatenate([dmg_at_centroid.index,
                                                              exp_at_centroid.index])))
        dmg_at_centroid=dmg_at_centroid.reindex(new_df.index)
        exp_at_centroid=exp_at_centroid.reindex(new_df.index)

        #compute PAA and MDR at centroid
        at_centroid_dict['PAA']=dmg_at_centroid['counts']/exp_at_centroid['counts']
        at_centroid_dict['MDR']=dmg_at_centroid['value']/exp_at_centroid['value']
        #get intensity at centroid
        at_centroid_dict[variable] = np.nanmax((exp_at_centroid['intensity'].replace('nan',np.nan), dmg_at_centroid['intensity'].replace('nan',np.nan)), axis=0)
        at_centroid_dict['n_exp'] = exp_at_centroid['counts']
        at_centroid_dict['n_dmgs'] = dmg_at_centroid['counts']
        at_centroid_dict['exp_val'] = exp_at_centroid['value']
        at_centroid_dict['dmg_val'] = dmg_at_centroid['value']
        at_centroid_dict['date'] = date_now

        if variable_2 is not None:
            at_centroid_dict[variable_2] = np.nanmax((exp_at_centroid['intensity_2'].replace('nan',np.nan), dmg_at_centroid['intensity_2'].replace('nan',np.nan)), axis=0)

        if get_PVA == True:
            at_centroid_dict['PVA'] = dmg_at_centroid['value_affected']/exp_at_centroid['value']
            at_centroid_dict['MDD'] = dmg_at_centroid['value']/dmg_at_centroid['value_affected']
        else:
            at_centroid_dict['PVA'] = None
            at_centroid_dict['MDD'] = None
    else:
        at_centroid_dict[variable] = exp_at_centroid['intensity'].replace('nan',np.nan)
        at_centroid_dict['n_exp'] = exp_at_centroid['counts']
        at_centroid_dict['date'] = date_now
        at_centroid_dict['n_dmgs'] = np.nan
        at_centroid_dict['dmg_val'] = np.nan

        if variable_2 is not None:
            at_centroid_dict[variable_2] = exp_at_centroid['intensity_2'].replace('nan',np.nan)

    #make dataframe
    at_centroid_df = pd.DataFrame({k: v for k, v in at_centroid_dict.items() if v is not None})
    at_centroid_df['centr_HL'] = at_centroid_df.index

    #only include data points with variable > 0 or damages > 0
    at_centroid_df['n_dmgs']=at_centroid_df['n_dmgs'].fillna(0)
    at_centroid_df['dmg_val']=at_centroid_df['dmg_val'].fillna(0)
    #at_centroid_df=at_centroid_df.loc[(at_centroid_df[variable]!=0) | (at_centroid_df['n_dmgs']!=0)]

    return at_centroid_df

# %% Main calibration functions
def empirical_calibration_per_exposure(hazard_object, exposure_object, damages,
    exposure_type = 'agriculture', variable = 'MESHS', hazard_object_2 = None,
    variable_2 = None, get_PVA = False, filter_year=None,
    dates_modeled_imp=None, roll_window_size=10, fraction_insured = 1):

    """
    main function for calibration - computes a range of variables as
    function of values of the hazard variable

    Parameters
    ----------
    hazard_object: climada.hazard
            hazard object used for calibraiton
    exposure_object: climada.entity.exposures.base.Exposures
            exposure object used for calibration
    damages: dict or climada.engine.impact.Impact
            either dictionary of exposure objects, dictionary of impact objects
            or a single impact objects hat contain damage claims
    exposure_type: str
            type of exposure ('agriculture' or 'buildings')
    variable: str
            hazard variable name
    hazard_object_2: climada.hazard
            hazard object with a secondary hazard variable that is traced trough the analysis
    variable_2: str
            secondary hazard variable name
    filter_year : tuple
            (yearMin,yearMax): to filter exposure data by year when an exposure was built
    dates_modeled_imp : np.array
        array of dt.datetime() with all dates where modelled impact is not zero
        (i.e. there is a hazard_intensity>0 somewhere)
    roll_window_size : int
        Size of rolling window. default to 10 for 10mm moving window in MESHS

    Returns
    -------

    ds: xarray.Dataset
        Dataset with additive calibration-relevant variables as function of hazard variable for each event
            count_cell: number of grid cells with non-zero exposure
            count_all: number of exposed assets
            value_all: total value of exposed assets
            count_dmg: number of affected assets
            value_dmg: total damage
            value_affected: total value affected
    df_all: pandas.Dataframe
        Dataframe which corresponds to the sum of the data in ds over all events
    ds_roll: xarray.Dataset
        ds with rolling average applied
    ds_roll_cut: xarray.Dataset
        ds_roll with cutoff for high values of hazard variable
    values_at_centroid_all: pandas.Dataframe
        A Dataframe with information on damages per gridpoint (PAA, MDD, MDR, n_dmgs),
        and exposure (n_exp), as well as hazard intensity (variable and variable_2)
    intensity_range: numpy.array
        Intensity array for hazard variable used in the calibration
    """
    ### DECLARE VARIABLES ##
    df_all = pd.DataFrame()  #output dataframe
    haz = copy.deepcopy(hazard_object) #create copy not to modify original hazard
    haz_2 = copy.deepcopy(hazard_object_2) #create copy not to modify original hazard

    exp_sel = exposure_object.copy() #create a copy of exposure not to modify original exposure


    # CHECKS
    if haz_2 is None:
        if variable_2 is None:
            pass
        else:
            print('Variable_2 set to None as no secondary hazard has been passed.')
            variable_2 = None
    else:
        if variable_2 is None:
            raise ValueError('Need to pass variable_2 to consider secondary hazard.')

    #assert no duplicate dates in hazard
    assert(len(haz.date)==len(np.unique(haz.date)))

    # ASSIGN CENTROIDS
    exp_sel.assign_centroids(haz)
    sc.assign_centroids_imp(damages,haz)

    # GET DATELISTS
    # from observed damages and modeled damages
    dt_dates=get_unique_dates_from_impacts(damages, dates_modeled_imp=dates_modeled_imp)
    # from damages only
    dt_dates_dmgs = [dt.datetime.fromordinal(date) for date in damages.date]

    # DEFINE INTENSITY RANGE
    try:
        intensity_range = INT_RANGE_DICT[variable]
    except KeyError:
        raise ValueError(f'Calibration not implemented for variable {variable}.'
                         ' Add value to INT_RANGE_DICT')

    # ASSIGN INTENSITY VALUES TO DAMAGES
    if variable == 'HKE': #round hazard intensity to steps of e.g. 50J/m2
            values, counts = np.unique(np.diff(intensity_range), return_counts=True)
            step=values[np.argmax(counts)] #get most common step size
            haz.intensity.data=np.round(haz.intensity.data/step)*step

    # assign intensities
    if hazard_object_2 is None:
        damages=assign_intensity_to_damages(damages, haz, variable)
    else:
        damages=assign_intensity_to_damages(damages, haz, variable,
                                            haz_2=haz_2, variable_2=variable_2)
    year = -np.inf #initialize year for year filtering

    # loop over dates
    for date_now in dt_dates:

        #filter exposure and recreate exp_sel once for each year
        if filter_year and date_now.year>year:
            year = date_now.year
            exp_sel = exposure_object.copy()
            exp_sel.assign_centroids(haz)
            exp_sel = filter_exp(exp_sel, filter_year[0], min(year, filter_year[1]))

        #add intensity values to exposure gdf of current date
        intensity = haz.intensity[haz.date == date_now.toordinal(), :].toarray().squeeze()
        exp_sel.gdf[variable] = intensity[exp_sel.gdf.centr_HL]
        if hazard_object_2 is not None:
            intensity_2 = haz_2.intensity[haz_2.date == date_now.toordinal(), :].toarray().squeeze()
            exp_sel.gdf[variable_2] = intensity_2[exp_sel.gdf.centr_HL]

        #for agriculture exposure compute an average field area per grid cell
        if exposure_type == 'agriculture':
                exp_sel.gdf['field_area']=exp_sel.gdf['value']/exp_sel.gdf['n_fields']
                #compute Swiss average field area to use if no field area can be defined
                #field_area_mean=exp_sel.gdf['value'].sum()/exp_sel.gdf['n_fields'].sum()

        # get counts, values ordered by intensity and centroid based information
        all_count, all_value, exp_at_centroid = \
                        get_exposed_assets(exp_sel,
                                           variable,
                                           exposure_type, variable_2 = variable_2)

        ### GET DAMAGED ASSETS

        #if damages have been reported at the current date, get them
        if date_now in list(dt_dates_dmgs):

            dmg_count, dmg_value, affected_value, dmg_at_centroid = get_damaged_assets(
                                damages,
                                variable,
                                date_now,
                                exposure_type, get_PVA,
                                variable_2 = variable_2,
                                fraction_insured = fraction_insured)

            #get at centroid data
            at_centroid_df = at_centroid_from_exp_and_dmg(exp_at_centroid, variable,
                                                          date_now, dmg_at_centroid = dmg_at_centroid,
                                                          variable_2 = variable_2, get_PVA = get_PVA)

        else: #if no damages are recorded add empty series for damage value/count

            dmg_count = pd.Series(dtype=int, name='count_dmg')
            dmg_value = pd.Series(dtype=int, name='value_dmg')
            affected_value = None

            #get at centroid data
            at_centroid_df = at_centroid_from_exp_and_dmg(exp_at_centroid, variable,
                                                          date_now, dmg_at_centroid = None,
                                                          variable_2 = variable_2, get_PVA = get_PVA)

            if get_PVA == True:
                affected_value = pd.Series(dtype=int, name='value_affected')
                # raise NotImplementedError('get_PVA with additional impact dates is not yet implemented. \
                #     All info is there->could be easily implemented')


        #concatenate at_centroid dataframe
        if 'values_at_centroid_all' not in locals():
            values_at_centroid_all = at_centroid_df
        else:
            values_at_centroid_all = pd.concat([values_at_centroid_all, at_centroid_df], ignore_index=True)

        #get number of grid cells with non-zero exposure per value of hazard variable
        vals,counts = np.unique(intensity[np.unique(exp_sel.gdf.centr_HL)],return_counts=True)
        cell_count = pd.Series(counts,index = vals).rename('count_cell')

        #create a pandas.DataFrame with all relevant data as function of hazard intensity
        df = pd.concat([cell_count,all_count,all_value,dmg_count,dmg_value,affected_value],axis = 1)

        #add a label for reported damages per day, and store full info in xarray
        if df_all.empty:
                df_all = df.fillna(0).copy()
                ds = df.reindex(intensity_range).fillna(0).to_xarray().expand_dims(dim={'date':[date_now]})
        else:
                df_all = df_all.add(df,fill_value = 0)
                ds_temp = df.reindex(intensity_range).fillna(0).to_xarray().expand_dims(dim={'date':[date_now]})
                ds = xr.concat([ds,ds_temp],dim = 'date')

    # rename the index of the output dataset
    ds = ds.rename({'index':variable})

    # SMOOTH DATA USING A RUNNING MEAN
    cut_off = CUT_OFF_DICT[variable]

    ds_roll, ds_roll_cut = smooth_data_running_mean(ds, variable, intensity_range,
                                                    cut_off, roll_window_size)

    return ds, df_all, ds_roll, ds_roll_cut, values_at_centroid_all, intensity_range

def compute_empirical_damage_functions(ds,ds_roll,ds_roll_cut,get_monotonic_fit=True):
    """
    Parameters
    ----------
    ds: xarray.Dataset
        Dataset with relevant values grouped per MESHS value
    ds_roll: xarray.Dataset
        moving average of ds
    ds_roll_cut: xarray.Dataset
        moving average of ds but with a cutoff after MESHS=60

    Returns
    --------
    df: pandas.Dataframe
        Dataframe with values computed over all events
    df_roll: pandas.Dataframe
        moving average of df
    df_roll_cut: pandas.Dataframe
        moving average of ds but with acutoff
    n_dmgs: pandas.Dataframe
        number of damage claims per intensity value
    """

    if 'value_affected' in ds.data_vars:
        use_PVA=True
    else:
        use_PVA=False

    #create dataframe including all events
    df=ds.sum(dim='date').to_dataframe()
    df_roll=ds_roll.sum(dim='date').to_dataframe()
    df_roll_cut=ds_roll_cut.sum(dim='date').to_dataframe()

    #compute measures relevant for impact function
    for dframe in [df_roll,df_roll_cut,df,ds,ds_roll,ds_roll_cut]:
        dframe['MDR'] = dframe.value_dmg/dframe.value_all
        dframe['PAA'] = dframe.count_dmg/dframe.count_all #PAA: Percent Assets Affected
        if use_PVA:
            dframe['MDD'] = dframe.value_dmg/dframe.value_affected
            dframe['PVA'] = dframe.value_affected/dframe.value_all #PVA: percent value affected
    n_dmgs=ds.sum(dim='date').count_dmg

    if get_monotonic_fit:
        df_roll_cut_fit=fit_monotonic(df_roll_cut)
        return df, df_roll, df_roll_cut, df_roll_cut_fit, n_dmgs
    else:
        return df, df_roll, df_roll_cut, n_dmgs

# %% Bootstrapping
def var_tuple(n_intensity,n_samples,variable):

    """ get tuple for mapping used to initialized xarray.Dataset
    var name: (tuple of dimension names, array-like)

    Parameters
    ---------
    n_intensity: int
        number of intensity values / length of dimension "variable"
    n_samples: int
        number of bootstrap samples / length of dimension 'b_sample'

    Returns
    --------
    tuple of dimension names and array-like
    """

    nan_arr = np.full((n_intensity, n_samples), np.nan)

    return ([variable, 'b_sample'], nan_arr)

def bootstrapping(ds, ds_roll,variable,n_samples,intensity_range,log_fit=False,
                  cut_off=60,keep_raw_values=False,do_cross_validation=False,
                  by_year=False):

    """
    Parameters
    ds: xarray.Dataset with date for each sample date
    ds_roll: ds but for rolling averages
    variable: hazard variable
    n_samples: number of bootstrap samples. If do_cross_validation=True,
                n_samples is the number of splits for cross validation
    cut_off: intensity value where function is cut off (default to 60 for MESHS=60mm)
    log_fit: boolean
        if True, fit a logistic function to the data
    keep_raw_values: boolean
        if True, also add raw values to output dataset
    do_cross_validation: boolean
        if True, do cross validation INSTEAD of bootstrapping
    n_splits: int
        number of splits for cross validation
    by_year: boolean
        if True, select cross validation splits by year

    Returns
    -------
    """
    #initialize bootstrap xarrays
    # https://elizavetalebedeva.com/bootstrapping-confidence-intervals-the-basics/

    #initialize bootstrap array with

    #get number of intensity values
    n_intensity=len(ds[variable])

    #initialize data variables
    data_vars = {'PAA': var_tuple(n_intensity,n_samples,variable),
                'MDR': var_tuple(n_intensity,n_samples,variable)}

    #determine if PVA should be calculated
    if 'value_affected' in ds.data_vars: #use PVA
        data_var_list = ['PAA', 'MDD', 'MDR','PVA']
        data_vars.update({'MDD': var_tuple(n_intensity,n_samples,variable),
                          'PVA': var_tuple(n_intensity,n_samples,variable)})

    else: #don't use PVA,MDD
        data_var_list = ['PAA', 'MDR']

    if keep_raw_values:
        #also add raw values to output dataset (e.g. count_dmg, value_all)
        data_var_list = data_var_list + ['count_cell','count_all','count_dmg',
                                         'value_all','value_dmg','value_affected']
        data_vars.update({'count_cell': var_tuple(n_intensity,n_samples,variable),
                        'count_all': var_tuple(n_intensity,n_samples,variable),
                        'count_dmg': var_tuple(n_intensity,n_samples,variable),
                        'value_all': var_tuple(n_intensity,n_samples,variable),
                        'value_dmg': var_tuple(n_intensity,n_samples,variable),
                        'value_affected': var_tuple(n_intensity,n_samples,variable)})

    ds_boot = xr.Dataset(data_vars=data_vars,
                         coords={variable: ds[variable],
                                 'b_sample': np.arange(n_samples)})

    #bootstrap array for rolling means
    ds_boot_roll = ds_boot.copy(deep=True)
    #bootstrap array for  roll cut
    ds_boot_roll_cut = ds_boot.copy(deep=True)

    #logistic fit data for each sample
    fit_data={var: pd.DataFrame(index=np.arange(n_samples), columns=['R2','params'])
              for var in data_var_list}


    if do_cross_validation:
        assert(n_samples<=10) #cross validation only for n_samples<=10

        #Create n_samples equally sized CV_splits
        if by_year:
            unique_years = np.sort(np.unique(ds.date.dt.year))
            n_years = unique_years.size
            split_size = np.floor(n_years/n_samples)
            CV_split_year = np.repeat(np.arange(n_samples), split_size)

            #pad the last CV_split to use all datapoints
            pad_length = n_years-len(CV_split_year)
            #assert that the last split is not more than 30% larger than others
            assert(pad_length/split_size<0.3)
            CV_split_year = np.pad(CV_split_year,(0,pad_length),mode='constant',
                        constant_values=max(CV_split_year))
            CV_split_dict = dict(zip(unique_years,CV_split_year))

            CV_split = np.vectorize(CV_split_dict.get)(ds.date.dt.year)
            dates_CV = pd.Series(data=CV_split,index=ds.date,name='CV_split')
        else:
            n_dates = len(ds.date)
            split_size = np.floor(n_dates/n_samples)
            CV_split = np.repeat(np.arange(n_samples), split_size)

            #pad the last CV_split to use all datapoints
            pad_length = n_dates-len(CV_split)
            CV_split = np.pad(CV_split,(0,pad_length),mode='constant',
                        constant_values=max(CV_split))
            dates_CV = pd.Series(data=CV_split,index=ds.date,name='CV_split')



    # bootstrapping n samples with replacement and with the size of the original dataset
    for num in np.arange(n_samples):

        #get random sample dates
        if do_cross_validation:

            # #get random sample dates
            # sel_samples = ds.date[num*split_size:(num+1)*split_size]
            # #get all dates except 1 split
            # sel_samples = np.setdiff1d(ds.date, sel_samples)
            sel_samples = dates_CV.index[dates_CV!=num].values

        else:
            sel_samples = np.random.choice(
                ds.date, size=len(ds.date), replace=True)


        #compute sum values for sample
        ds_now = ds.sel(date=sel_samples).sum(dim='date')
        ds_roll_now = ds_roll.sel(date=sel_samples).sum(dim='date')
        ds_roll_now_cut = ds_roll_now.copy(deep=True)
        ds_roll_now_cut.loc[{variable: slice(cut_off, intensity_range.max())}] = ds_roll_now_cut.loc[{
            variable: slice(cut_off, intensity_range.max())}].sum(dim=variable)

        #Compute relevant metrics for this sample
        for dset in [ds_now, ds_roll_now, ds_roll_now_cut]:
            dset['MDR'] = dset.value_dmg/dset.value_all
            dset['PAA'] = dset.count_dmg/dset.count_all
            if 'PVA' in data_var_list:
                dset['MDD'] = dset.value_dmg/dset.value_affected
                dset['PVA'] = dset.value_affected/dset.value_all

        #add sample metrics to output arrays
        ds_boot.loc[{'b_sample': num}] = ds_now[data_var_list]
        ds_boot_roll.loc[{'b_sample': num}] = ds_roll_now[data_var_list]
        ds_boot_roll_cut.loc[{'b_sample': num}
                             ] = ds_roll_now_cut[data_var_list]

        if log_fit == True:

            for var in data_var_list:
                xData=intensity_range[1:41]
                yData=ds_roll_now[var][1:41].values

                #initial parameter estimation
                p0=[0.2,0.075,0.2,48]

                try:
                    xrange, modelPredictions, Rsquared, parameters = fit_log_curve(xData,yData,p0)
                except:
                    Rsquared = None
                    parameters = None

                fit_data[var].loc[num,'R2']=Rsquared
                fit_data[var].loc[num,'params']=parameters
        else:
            fit_data=None

    #Rename datasets in case of cross validation
    if do_cross_validation:
        ds_boot = ds_boot.rename({'b_sample': 'CV_split'})
        ds_boot_roll = ds_boot_roll.rename({'b_sample': 'CV_split'})
        ds_boot_roll_cut = ds_boot_roll_cut.rename({'b_sample': 'CV_split'})
        if by_year:
            return (ds_boot,ds_boot_roll,ds_boot_roll_cut, fit_data,
                    dates_CV,pd.Series(CV_split_dict,name='CV_split'))
        else:
            return ds_boot,ds_boot_roll,ds_boot_roll_cut, fit_data, dates_CV
    else:
        return ds_boot,ds_boot_roll,ds_boot_roll_cut, fit_data

def fit_monotonic(df):
    """fit monotonic function to variables in df (MDR, PAA, MDD, PVA) and
    return a new dataframe with fitted values


    Parameters
    ----------
    df : pandas.Dataframe
        dataframe with values of MDR, PAA, MDD, PVA that a monotonic fit needs
        to be applied to
    Returns
    -------
    df_monotonic : pandas.Dataframe
        dataframe with monotonic fot of values in df

    """

    var_list=['MDR', 'PAA', 'MDD', 'PVA']
    data_out={}
    #logistic fit data for each sample
    for var in var_list:

        if var in df.keys():
            if df.index[0]==0:# skip_zero
                assert(np.where(df.index==0) == np.array([0]))
                x=df.index[1:]
                y=df[var][1:]
            else:
                x=df.index
                y=df[var]
            #avoid nan values
            no_nan = ~np.isnan(x) & ~np.isnan(y)
            x = x[no_nan]
            y = y[no_nan]

            monotone_fit = sc.both.smooth_monotonic(x,y,plot=False)

            data_out[var]=monotone_fit

    df_monotonic=pd.DataFrame(data_out,index=x)

    return df_monotonic

def extend_haz(haz,dates_to_extent):
    """extend hazard object with empty events for calibration

    Args:
        haz (climada.hazard): hazard object to extend
        dates_to_extent (array-like): dates to extend hazard object with
    """
    #check that none of the new dates already exists
    assert(np.all(np.isin(dates_to_extent,haz.date)==False))
    add_intensity = sparse.csr_matrix((len(dates_to_extent),haz.intensity.shape[1]),dtype=float)
    stacked_int = sparse.vstack([haz.intensity,add_intensity])

    dates_combined = np.concatenate([haz.date,dates_to_extent])
    haz2 = Hazard(
        haz_type = haz.haz_type ,
        units = haz.units,
        centroids = haz.centroids,
        event_id = np.arange(1,len(dates_combined)+1),
        frequency = np.full_like(dates_combined,haz.frequency[0],dtype=float),
        event_name = [dt.datetime.fromordinal(d).strftime("ev_%Y-%m-%d") for d in dates_combined],
        date = dates_combined,
        orig =  np.full_like(dates_combined,haz.orig[0],dtype=bool),
        intensity = stacked_int
        )
    haz2.check()
    return haz2