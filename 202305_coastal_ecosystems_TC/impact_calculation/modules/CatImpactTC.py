import copy
import pandas as pd
from climada.engine import Impact
import scipy.sparse

def imp_per_cat(exp, impf_set, haz, path):
    """
    calculates the absolute value impacted by hazard category
    :param exp: exposure dataframe
    :param impf_set: Impact function set
    :param haz: hazard dataframe
    :param path: string file path for output csv files
    :return:
    """
    all_impact_df = pd.DataFrame()
    all_impact_df['longitude'] = exp.gdf['longitude']
    all_impact_df['latitude'] = exp.gdf['latitude']
    haz_cat = {
        1: [33, 42],
        2: [43, 49],
        3: [50, 58],
        4: [58, 70],
        5: [70, 200]
    }
    for cat in haz_cat.keys():
        exp_cop = copy.copy(exp)
        exp_cop.gdf['impf_TC'] = cat
        imp = Impact()
        imp.calc(exp_cop, impf_set, haz, save_mat=True)
        imp_mat = imp.imp_mat
        scipy.sparse.save_npz(f'{path}_pop_imp_mat{cat}.npz', imp_mat)
        imp.write_csv(f'{path}_pop_cat{cat}_impact.csv')
        all_impact_df[f'eai_exp_cat{cat}'] = imp.eai_exp
    return all_impact_df
