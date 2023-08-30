# An open-source radar-based hail damage model for buildings and cars

The scripts here reproduce the main results of the paper:
Schmid T., Portmann R., Villiger L. Schröer K., Bresch D. N. An open-source radar-based hail damage model for buildings and cars, Natural Hazards and Earth System Sciences (2023+)

Publication status: submitted
Contact: [Timo Schmid](timo.schmid@usys.ethz.ch)

## Content

Notebooks/

Jupyter notebooks to reproduce figures that appear in the paper

Scripts/

Python scripts used to process data save intermediate results. In particular:
* calibration_per_exposure.py performs the spatially explicit calibration of vulnerability functions and saved impact function parameters and produces impact function plots.
* hail_main.py uses calibrated impact functions to estimate hail damages for different hazard-expsoure combinations and saves skill scores.
* data_processing/grid_cantonal_data.py collects damage and exposure data from all 4 cantons and combines it to a 1km gridded dataset.
* data_processing/impact_to_netcdf.ipynb saves impact data as netcdf file.

scClim/

Contains functions used for data pre-processing, calibration, visualizing, as well as utility functions.