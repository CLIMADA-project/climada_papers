"""
Script to create calculate impact for given population exposure
Author: Sarah Hülsen
"""

from pathlib import Path

from climada.entity import Exposures
from climada.hazard import Hazard
from modules.CatImpactFuncTC import cat_impf_set
from modules.CatImpactTC import imp_per_cat


basin = 'AP'    # TC basin (e.g. AP, IO, SH, WP)
year = '2020'   # population exposure year (e.g. 2000 or 2020)
scenario = 'hist_STORM'  # current climate or SSP585 TC hazard (hist_STORM or clim_STORM_{model})

# Data paths
input_path = Path('../results/intermediate/')
output_path = Path('../results/final/')
haz_fn = f'TC_{basin}_0300as_STORM.hdf5'
pop_fn = f'exp_wp_{basin}_{year}.hdf5'

# Load hazard and exposure
haz = Hazard.from_hdf5(f'{input_path}{haz_fn}')
haz.check()
exp_pop = Exposures.from_hdf5(f'{input_path}{pop_fn}')
exp_pop.check()

# Define impact functions
impf_set = cat_impf_set('TC')

# Calculate impact
imp_cat_all = imp_per_cat(exp=exp_pop, impf_set=impf_set, haz=haz, path=output_path)
imp_cat_all.to_csv(f'{output_path}{basin}_{scenario}_STORM_{year}_wp_all_imp.csv')