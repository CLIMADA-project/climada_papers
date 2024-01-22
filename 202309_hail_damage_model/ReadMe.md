# An open-source radar-based hail damage model for buildings and cars

The scripts here reproduce the main results of the paper:
Schmid T., Portmann R., Villiger L., Schröer K., Bresch D. N. (2023+) An open-source radar-based hail damage model for buildings and cars, Natural Hazards and Earth System Sciences, https://doi.org/10.5194/nhess-2023-158

Publication status: [accepted](https://doi.org/10.5194/nhess-2023-158)

Contact: [Timo Schmid](timo.schmid@usys.ethz.ch)

## Content

### test_notebook.ipynb
Jupyter notebook that runs through the calibration and model evaluation with a test dataset. The data is artificially created, but has the same format as the original data.
The notebook contains the relevant code from 2 scripts that are used in the calibration with the real data: *calibration_main.py* and *hail_main.py*.

To perform a calibration as in this publication, 3 main dataset are needed:
* **Hazard**: Gridded data of a natural hazard. For format see *test_data/test_meshs.nc*
* **Exposure**: Tabular data of Exposure values and coordinates.  For format see *test_data/test_exp.nc*
* **Damage**: Tabular data of reported damages **with spatial coordinates**. For format see *test_data/test_dmg.nc*

The actual data used in the paper cannot be shared as it is from private insurance companies and the Swiss national weather service. The scripts in the remaining folders reproduce the results and figures from the paper with the input data as described in the [publication](https://doi.org/10.5194/nhess-2023-158).

### ./notebooks/

Jupyter notebooks to reproduce figures that appear in the paper

### ./scripts/

Python scripts used to process data save intermediate results. In particular:
* *calibration_main.py* performs the spatially explicit calibration of vulnerability functions and saved impact function parameters and produces impact function plots.
* *hail_main.py* uses calibrated impact functions to estimate hail damages for different hazard-expsoure combinations and saves skill scores.
* *event_definition.py* performs the POH-based pre-processing of the building and car damage data
* *data_processing/grid_cantonal_data.py* collects damage and exposure data from all 4 cantons and combines it to a 1km gridded dataset.
* *data_processing/impact_to_netcdf.ipynb* saves impact data as netcdf file.

### ./scClim/

Contains functions which are called in other scripts for data pre-processing, calibration, visualizing, as well as utility functions and constants.

## Requirements
Requires:
* Python 3.9+ environment (best to use conda for CLIMADA repository)
* _CLIMADA_ repository version 3.3+:
        https://wcr.ethz.ch/research/climada.html
        https://github.com/CLIMADA-project/climada_python
* Exposure and damage data for the calibration. The hail damage data used in the paper are not public and only available within the [scClim](https://scclim.ethz.ch/) project. Calculation can be reproduced with other user-provided data as shown in *test_notebook.ipynb*.