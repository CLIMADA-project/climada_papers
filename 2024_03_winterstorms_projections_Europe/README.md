# Projections and uncertainties of winter windstorm damage in Europe in a changing climate
These scripts reproduce the main results of the paper:
Projections and uncertainties of winter windstorm damage in Europe in a changing climate

Luca G. Severino (1), Chahan M. Kropf (2,3) Hilla Afargan-Gerstman (1), Christopher Fairles (2), Andries Jan de Vries (4), Daniela I.V. Domeisen (4,1), and David N. Bresch (2,3)
1 Institute for Atmospheric and Climate Science, ETH Zürich, Zürich, Switzerland
2 Institute for Environmental Decisions, ETH Zürich, Zürich, Switzerland
3 Federal Office of Meteorology and Climatology MeteoSwiss, Zürich, Switzerland
4 Institute of Earth Surface Dynamics, University of Lausanne, Lausanne, Switzerland
Correspondence: Luca Severino (luca.severino@usys.ethz.ch)

## Content:
This folder contains material to reproduce the main results from the publication Projections and uncertainties of winter windstorm damage in Europe in a changing climate. The python jupyter notebook 1_load_process_data.ipynb serves to load and process the hazard data from a data archive containing CMIP6 netcdf files; the notebook 2_damage_calc.ipynb computes the damages projections; and the notebook 3_uncertainty-sensitivity_analyses.ipynb reproduces the results from the uncertainty and sensitivity analysis. The SL_bias_correction.py contains the script to bias-correct the hazard data taken from  https://github.com/samluethi/HeatMortality2021/blob/main/Data_preparation/bias_correction.py, and constants.py contains constants required for the calculations.

## Notes on hazard file generation
The section Damage Calculation of script 2_damage_calc.ipynb has been used to generate the hazard data present on the CLIMADA API: https://climada.ethz.ch/data-api/v2/docs.
The basic steps are described in the publication manuscript, currently available as a preprint at https://doi.org/10.5194/egusphere-2023-205. The daily surface wind maxima
are first bias corrected using ERA5 reanalysis as reference. Then windstorms days are detected following the algorithm:

${Stormy day}_t \iff \sum_i \{a_i|[(v_{i,t} \geq v_{i,98}) \ \& \ (v_{i,t} \geq 15)] \} \geq A_{min}$

where $a_i$ is the area of the grid cell $i$, $v_{i,t}$ is the daily sfcWindmax intensity at grid cell $i$ on the considered day $t$, $v_{i,98}$ is the 98th percentile of the daily sfcWindmax at grid cell $i$, computed over the winter periods of the historical period, and $A_{min}$ is the area threshold parameter.
Days that do not fulfill this conditions are ignored. Gridcells that for which the daily surface wind maxima intensity values are below the absolute or relative thresholds are marked as 0.


## Requirements

requires access to CMIP6 model outputs as .netcdf files (see https://wcrp-cmip.org/cmip-data-access/#esgf for more information), as well as LitPop metadata (see tutorial https://github.com/CLIMADA-project/climada_python/blob/main/doc/tutorial/climada_entity_LitPop.ipynb).

Requires Python 3.8+ and CLIMADA v4.1.1 (or later):
https://github.com/CLIMADA-project/climada_python/
Documentation: https://github.com/CLIMADA-project/climada_python/blob/master/doc/source/install.rst

## Documentation:

* Publication available as a preprint at https://doi.org/10.5194/egusphere-2023-205 .

Documentation for CLIMADA is available on Read the Docs:

* [online (recommended)](https://climada-python.readthedocs.io/en/stable/)
* [PDF file](https://buildmedia.readthedocs.org/media/pdf/climada-python/stable/climada-python.pdf)

-----
## History
created 28 March 2024

-----

www.wcr.ethz.ch