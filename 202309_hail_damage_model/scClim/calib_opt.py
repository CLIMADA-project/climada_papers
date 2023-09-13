"""
Impact function calibration functionalities:
    Optimization and manual calibration

based in climada.engine.calib_opt.py and extended by Timo to include
bayesian optimization
"""

import datetime as dt
import copy
import numpy as np
import pandas as pd
import warnings
from scipy import interpolate
from scipy.optimize import minimize
from itertools import combinations
import matplotlib.pyplot as plt
from climada.entity import ImpfTropCyclone, ImpactFunc


def fit_emanuel_impf_to_emp_data(emp_df,pbounds,opt_var='MDR',options=None,
                                 optimizer='Nelder-Mead',plot=True,param_start=None):
    """Fit emanuel-type impact function to empirical data

    Args:
        emp_df (pd.DataFrame): DF with empirical data, including 'MDR', 'PAA', and 'count_cell'
        pbounds (dict): dictionary of parameter bounds
        opt_var (str, optional): Variable from emp_df to fit data to. Defaults to 'MDR'.
        options (dict, optional): Additional options for the Bayesian optimizer. Defaults to None.
        optimizer (str, optional): Choice of optimizer. Defaults to 'Nelder-Mead'.
        plot (bool, optional): Whether or not to plot the data. Defaults to True.
        param_start (dict, optional): Initial parameter values. Defaults to None.

    Raises:
        ValueError: if get_emanuel_impf returns an error except for v_half <= v_thresh

    Returns:
        tuple: Parameters, optimizer object, and impact function
    """
    if not emp_df.index.name:
        raise ValueError('Careful, emp_df.index has no name. Double check if the \
                            index corresponds to the intensity unit!')

    def weighted_MSE(**param_dict):
        try:
            impf = get_emanuel_impf(**param_dict,intensity=emp_df.index.values)
        except ValueError as e:
            if 'v_half <= v_thresh' in str(e):
                #if invalid input to Emanuel_impf (v_half <= v_thresh), return zero
                return np.inf
            else:
                raise ValueError(f'Unknown Error in init_impf:{e}. Check inputs!')
        if opt_var == 'MDR':
            mdr = impf.mdd*impf.paa
            SE = np.square(mdr-emp_df.MDR)
        elif opt_var == 'PAA':
            #PAA is used implicitly, output MDD*PAA will be PAA
            paa = impf.mdd*impf.paa 
            SE = np.square(paa-emp_df.PAA)
        elif opt_var == 'MDD':
            #MDD is used implicitly, output MDD*PAA will be MDD
            mdd = impf.mdd*impf.paa
            SE = np.square(mdd-emp_df.MDD)

        SE_weigh = SE*emp_df.count_cell
        SE_weigh_noZero = SE_weigh[emp_df.index.values!=0]
        MSE = np.mean(SE_weigh_noZero)
        return MSE

    if optimizer == 'Nelder-Mead' or 'trust-constr':
        if param_start is None:
            #use mean of bounds as starting point
            param_means = [(v[0]+v[1])/2 for v in pbounds.values()]
        else: 
            param_means = param_start
        param_dict = dict(zip(pbounds.keys(),param_means))
        bounds = [(v[0],v[1]) for v in pbounds.values()]
        x0 = list(param_dict.values())

        #define function that returns the MSE, with an array x as input
        def mse_x(x):
            param_dict_temp = dict(zip(param_dict.keys(), x))
            return weighted_MSE(**param_dict_temp)
        
        if optimizer == 'trust-constr':
            print(param_dict,bounds)
            np.testing.assert_array_equal((list(pbounds.keys())[:2]),['v_thresh', 'v_half'])
            cons = ({'type': 'ineq', 'fun': lambda x:  x[1] - x[0]},#v_half-v_tresh is to be non-negative
                    {'type':'ineq', 'fun': lambda x: x[2]})         #scale is to be non-negative
            options= None
            
        elif optimizer == 'Nelder-Mead':
            cons = None
            options={'xtol': 1e-5, 'disp': True, 'maxiter': 500}
        res = minimize(mse_x, x0,
                        bounds=bounds,
                        constraints=cons,
                        method = optimizer,
                        options=options)

        optimizer = res
        param_dict_result = dict(zip(param_dict.keys(), res.x))
        print(param_dict_result)
        impf = init_impf('emanuel_HL',param_dict_result,emp_df.index.values)[0]

        if plot:
            ax=impf.plot(zorder=3)
            title = [f'{key}: {param_dict_result[key]:.2f}' if param_dict_result[key]>0.1 else 
                     f'{key}: {param_dict_result[key]:.2e}' for key in param_dict_result.keys()]
            ax.set(ylim=(0,max(impf.mdd*100)),title=title)
            #add empirical function to plot
            ax.plot(emp_df.index,emp_df[opt_var]*100,label=f'Empirical {opt_var}')
            plt.legend()

    return param_dict_result, optimizer, impf


def init_impf(impf_name_or_instance, param_dict,intensity_range, 
              df_out=pd.DataFrame(index=[0])):
    """create an ImpactFunc based on the parameters in param_dict using the
    method specified in impf_parameterisation_name and document it in df_out.

    Parameters
    ----------
    impf_name_or_instance : str or ImpactFunc
        method of impact function parameterisation e.g. 'emanuel' or an
        instance of ImpactFunc (not implemented here)
    param_dict : dict, optional
        dict of parameter_names and values
        e.g. {'v_thresh': 25.7, 'v_half': 70, 'scale': 1}
    intensity_range : array
        tuple of 3 intensity numbers along np.arange(min, max, step)
    Returns
    -------
    imp_fun : ImpactFunc
        The Impact function based on the parameterisation
    df_out : DataFrame
        Output DataFrame with headers of columns defined and with first row
        (index=0) defined with values. The impact function parameters from
        param_dict are represented here.
    """

    impact_func_final = None
    if isinstance(impf_name_or_instance, str):
        if impf_name_or_instance == 'emanuel':
            impact_func_final = ImpfTropCyclone.from_emanuel_usa(**param_dict)
            impact_func_final.haz_type = 'TC'
            impact_func_final.id = 1
            df_out['impact_function'] = impf_name_or_instance

        if impf_name_or_instance == 'emanuel_HL':
            impact_func_final = get_emanuel_impf(
                **param_dict,intensity=intensity_range,haz_type='HL')
            df_out['impact_function'] = impf_name_or_instance

        elif impf_name_or_instance == 'sigmoid_HL':
            assert('L' in param_dict.keys() and 'k' in param_dict.keys() and 
                   'x0' in param_dict.keys()) 
            impact_func_final = ImpactFunc.from_sigmoid_impf(
                **param_dict,intensity=intensity_range)#,haz_type='HL')
            impact_func_final.haz_type = 'HL'
            if intensity_range[0]==0 and not impact_func_final.mdd[0]==0:
                warnings.warn('sigmoid impact function has non-zero '
                              'impact at intensity 0. Setting impact to 0.')
                impact_func_final.mdd[0]=0
            df_out['impact_function'] = impf_name_or_instance

    elif isinstance(impf_name_or_instance, ImpactFunc):
        raise NotImplementedError('Not implemented here.')
    
    for key, val in param_dict.items():
        df_out[key] = val
    return impact_func_final, df_out

def get_emanuel_impf(v_thresh=20, v_half=60, scale=1e-3,power=3,
                    impf_id=1, intensity=np.arange(0, 110, 1), 
                    intensity_unit='mm',haz_type='HL'):
    """
    Init TC impact function using the formula of Kerry Emanuel, 2011:
    https://doi.org/10.1175/WCAS-D-11-00007.1

    Parameters
    ----------
    impf_id : int, optional
        impact function id. Default: 1
    intensity : np.array, optional
        intensity array in m/s. Default:
        5 m/s step array from 0 to 120m/s
    v_thresh : float, optional
        first shape parameter, wind speed in
        m/s below which there is no damage. Default: 25.7(Emanuel 2011)
    v_half : float, optional
        second shape parameter, wind speed in m/s
        at which 50% of max. damage is expected. Default:
        v_threshold + 49 m/s (mean value of Sealy & Strobl 2017)
    scale : float, optional
        scale parameter, linear scaling of MDD.
        0<=scale<=1. Default: 1.0
    power : int, optional
        Exponential dependence. Default to 3 (as in Emanuel (2011))

    Raises
    ------
    ValueError

    Returns
    -------
    impf : ImpfTropCyclone
        TC impact function instance based on formula by Emanuel (2011)
    """
    if v_half <= v_thresh:
        raise ValueError('Shape parameters out of range: v_half <= v_thresh.')
    if v_thresh < 0 or v_half < 0:
        raise ValueError('Negative shape parameter.')
    if scale > 1 or scale <= 0:
        raise ValueError('Scale parameter out of range.')

    impf = ImpactFunc(haz_type=haz_type, id=impf_id,intensity=intensity,
                        intensity_unit=intensity_unit,name='Emanuel-type')
    impf.paa = np.ones(intensity.shape)
    v_temp = (impf.intensity - v_thresh) / (v_half - v_thresh)
    v_temp[v_temp < 0] = 0
    impf.mdd = v_temp**power / (1 + v_temp**power)
    impf.mdd *= scale
    return impf
