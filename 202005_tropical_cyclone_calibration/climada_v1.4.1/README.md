# Regional tropical cyclone impact functions for globally consistent risk assessments (CLIMADA v.1.4.1+)

This repository contains a folder for the code used in a research article. The code replicates the figures and explains software components implemented in [CLIMADA Python](https://github.com/CLIMADA-project/climada_python).

Eberenz, S., Lüthi, S., and Bresch, D. N.: Regional tropical cyclone impact functions for globally consistent risk assessments. Accepted for publication by the journal Natural Hazards and Earth System Sciences (NHESS). Initial submission: https://doi.org/10.5194/nhess-2020-229.       

The calibrated impact functions are part of the main branch of [CLIMADA v.1.5.1+](https://github.com/CLIMADA-project/climada_python/releases) 
in the class [climada.entity.impact_funcs.trop_cyclone.IFSTropCyclone](https://github.com/CLIMADA-project/climada_python/blob/main/climada/entity/impact_funcs/trop_cyclone.py).

This is part of the [CLIMADA_papers repository](https://github.com/CLIMADA-project/climada_papers/). An archived version of this part of CLIMADA_papers is permanently available at: [![DOI](https://zenodo.org/badge/333036354.svg)](https://zenodo.org/badge/latestdoi/333036354)

Contact: [Samuel Eberenz](mailto:samuel.eberenz@usys.ethz.ch)

## Content:

#### tc_calibration_figures_and_tables.ipynb
Jupyter notebook to reproduce the key figures and tables in the paper.
The notebook was tested for CLIMADA v1.4.1 and v1.5.1.

#### tc_calibration_main.py
Python script to replicate the calibration.

Calibration steps (set in list CALIB):
* 1:  Loading or initiating HAZARD and EXPOSURE sets (required for 2, 3, 5, 6, 7)
* 2:  With one single global impact functions, calculate EDR. Extract NRD from EM-DAT.
        Compute ratio EDR=SED/NRD. Save to CSV file. Required for CALIB 3.
* 3:  Core calibration. Loop over v_half_range to compute impact and ratio for each matched event/country. Required for CALIB 4.
* 4:  Optimization: Compute cost functions and find optimized v_half for each cost function and region,
        Save calibration results to CSV.
* 5:  Compute annual average damage (AAD) per country for comparison with GAR 2013 (calibrated and uncalibrated)
* 6:  Compute damage time series, trend and significance, as well as standard deviation of annual damage for EM-DAT and CLIMADA per country and region.

Requires:
* CLIMADA v.1.4.1 (not tested for CLIMADA 1.5)
* *tc_calibration_functions.py*, *tc_calibration_config.py*, *impact_data_stable_202006.py*, and *if_trop_cyclone_stable_202006.py*.

_Important_: You need to customize paths and parameters in *tc_calibration_config.py* before running *tc_calibration_main.py*.

Please note that the preparation of the exposure and hazard files (step CALIB=1)
and the core calibration (step CALIB=3) take a lot of resources and might need to
be run on a supercomputing cluster. Please contact the authors if you have troubles
replicating these resource intensive steps.

#### tc_calibration_config.py
Configuration script.
Customize paths and parameters before running tc_calibration_main.py,
especially:
* *DATA_DIR*: path to directory with input data.
* *TRACK_DIR*: path to directory with TC tracks downloaded from IBtRACS.
* *TRACK_FOLDER*: name of folder with subset of TC tracks to use (place in TRACK_DIR).
* *EMDAT_CSV*: full path of CSV file with data downloaded from EM-DAT.
* *CALIB*: choose which calibration steps to run.
* *v_step*: the step size for V_half in the incremental calibration. Set to a larger value to speed up CALIB=3 at the cost of precision.

#### tc_calibration_functions.py, if_trop_cyclone_stable_202006.py, impact_data_stable_202006.py
Scripts with functions required for *tc_calibration_main.py*.


## Requirements

Requires:
* Python 3.6+ environment (best to use conda for CLIMADA repository)
* _CLIMADA_ repository version 1.4.1+:
        https://wcr.ethz.ch/research/climada.html
        https://github.com/CLIMADA-project/climada_python
* TC track data from _IBTrACS_ v4, 1980-2017 (get data online or ask the author):
        https://www.ncdc.noaa.gov/ibtracs/
* _EM-DAT_ impact data for tropical cyclone 1980-2017 (get data from EM-DAT or ask the author):
        https://www.emdat.be
        https://public.emdat.be/

The calibration script was not tested for CLIMADA versions other than 1.4.1.

## Documentation:

Documentation for CLIMADA is available on Read the Docs:
* [online (recommended)](https://climada-python.readthedocs.io/en/stable/)
* [PDF file](https://buildmedia.readthedocs.org/media/pdf/climada-python/stable/climada-python.pdf)

More infromation on the methodology implemented with these specific scripts is provided in the main paper (in preparation).

If script fails, revert CLIMADA version to release v1.4.1:
* from [GitHub](https://github.com/CLIMADA-project/climada_python/releases/tag/v1.4.1)
* or from the [ETH Data Archive](http://doi.org/10.5905/ethz-1007-252)

...and download the required EM-DAT and IBTrACS data or request it form the authors.

-----

## Updates
* 2021-01-21: update of figures and tables in Jupyter Notebook after peer review and acceptance for publication in NHESS; preparation for publication with Zenodo DOI.
* 2020-07-09: initial commit [3e5696ca0c1e463c7656379eb95c4e03fb91e32a](https://github.com/CLIMADA-project/climada_papers/commit/3e5696ca0c1e463c7656379eb95c4e03fb91e32a)
-----

www.wcr.ethz.ch
