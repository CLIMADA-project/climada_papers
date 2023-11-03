# Crop hail damage footprint model
These scripts can be used to reproduce the figures and analyses as well as apply the model developed in the following paper:

- **Portmann R., Schmid T., Villiger L., Bresch D., Calanca P.** Modelling crop hail damage footprints with single-polarization radar: The roles of spatial resolution, hail intensity, and cropland density, submitted to Natural Hazards and Earth System Sciences.

**Corresponding Author**: 
Raphael Portmann, Agroscope Reckenholz, Zurich, Switzerland\
E-Mail: raphael.portmann@agroscope.admin.ch or raphael.portmann@alumni.ethz.ch

Requires [**Climada 4.0.1**](https://climada-python.readthedocs.io/en/v4.0.1/misc/README.html) or higher
**python version**: 3.9.16

The scripts use hazard (hail), exposure (crops), and damage (hail_damage_crops) data that is provided with this publication via the [Climada Data-API](https://climada.ethz.ch/data-types/).


### Scripts

`run_damage_footprint_model.ipynb`
   - Notebook to apply the the hail damage footprint model and compute verification scores


`get_per_gridcell_data.py`
   - produces calibration files (*.p) in the subdirectory data/data_at_centroid/ that are used in many of the visualization scripts
   - suggested use: `python get_per_gridcell_data.py -croptypes wheat barley maize rapeseed grapevine -hazard_var MESHS -res 1 2 4 8 16 32`
   
   
`plot_skill_resolution_Figure_1.ipynb`
   - Notebook to reproduce Figure 1 of the paper

`plot_skill_metrics_single_event_Figures_2_3.ipynb`
   - Notebook to reproduce Figures 2 and 3 of the paper
   - Only works with hazard, exposure and damage as .netcdf files that are not provided via the Climada Data API. Available upon request.

`plot_skill_MESHS_threshold_Figure_4.ipynb`
   - Notebook to reproduce Figure 4 of the paper

`plot_performance_diagrams_Figure_5.py`
   - Skript to reproduce Figure 5 of the paper (suggest to use Spyder to allow for manual placement of bias labels)

`plot_skill_exposure_density_Figure_6.ipynb`
   - Skript to reproduce Figure 6 of the paper

`plot_performance_diagrams_sensitivity_Figure_7.py`
   - Skript to reproduce Figure 7 of the paper (suggest to use Spyder to allow for manual placement of bias labels)

`plot_performance_diagram.ipynb`
   - Notebook to produce a performance diagram

`plot_cropland_number_density_maps_Figure_A1.ipynb`
   - Notebook to reproduce Figure A1 of the paper 
   - Only works with hazard, exposure and damage as .netcdf files that are not provided via the Climada Data API. Available upon request.

`investigate_exposure_properties_Figures_A2_3.ipynb`
   - Notebook to plot Figures A2 and A3 and get the numbers mentionned in the paper

`plot_and_investigate_damage_exposure_mismatch.ipynb`
   - Notebook to reproduce numbers mentionned in Section 2.4 of the paper (uncertainty of random sampling and mismatch between damage and exposure data)
   - Only works with additional randomly sampled data (available upon request).

`store_gridded_data_as_exposure_and_impact_hdf5.ipynb`
   - Notebook to load and store gridded data as exposure, impact and hazard .hdf5 files

`calibration.py`
   - contains code used in `get_per_gridcell_data.py` to produce the calibration files

`utility.py`
   - contains helper code used in most files

