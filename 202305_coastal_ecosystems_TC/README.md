# Instructions

These scripts can be used to reproduce the results from the paper

"Global coastal protection benefits of ecosystems - past, present, and under climate change"

Hülsen, S.; McDonald, R.I.; Chaplin-Kramer, R.; Bresch, D.N.; Sharp, R.; Worthington, T.; Kropf, C.M.

See table below for an overview of input data. 

To do the impact calculation, first run the script `exposure_ecosystem_matching.py`, then `Impact_calculation.py` in the `impact_calculation` folder.

To perform an analysis of the impact data, run the scripts in the `impact_data_analysis` folder in the following order: `global_analysis.py`, `region_matching.py`, `regional_analysis.py`, `per_country_analysis.py`. 


| Data type                              | Description                                      | Data source                                                                                                                                                                                                                                                                                                           | Data preparation                                                                     |   |   |   |   |
|----------------------------------------|--------------------------------------------------|-----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------|--------------------------------------------------------------------------------------|---|---|---|---|
| Landcover data                         | Global land cover satellite data                 | SA CCI-LC, 2017. ESA/CCI viewer [WWW Document]. URL http://maps.elie.ucl.ac.be/CCI/viewer/download.php                                                                                                                                                                                                                | Compute coastal ecosystem protective rank using InVEST Coastal Vulnerability model   |   |   |   |   |
| Population exposure                    | Global gridded population data, 1km resolution   | WorldPop, 2018. Global 1km Population. https://doi.org/10.5258/SOTON/WP00647                                                                                                                                                                                                                                          | Clip to 10m elevation above sea level using DEM, Vectorize, Clip to TC basin extents |   |   |   |   |
| Coastal ecosystem protective rank data | Shoreline points, 450m resolution                | Natural Capital Project, 2019. InVEST [WWW Document]. Nat. Cap. Proj. URL https://naturalcapitalproject.stanford.edu/software/invest                                                                                                                                                                                  | Clip to TC basin extents                                                             |   |   |   |   |
| STORM TC hazard data                   | Synthetic tropical cyclone tracks                | Bloemendaal, N., Haigh, I.D., de Moel, H., Muis, S., Haarsma, R.J., Aerts, J.C.J.H., 2020. Generation of a global synthetic tropical cyclone hazard dataset using STORM. Sci. Data 7, 40. https://doi.org/10.1038/s41597-020-0381-2                                                                                   | Compute windfields using CLIMADA Hazard module, with 300 arc seconds resolution      |   |   |   |   |
| STORM-C TC hazard data                 | Synthetic tropical cyclone tracks SSP585 in 2050 | Bloemendaal, N., de Moel, H., Martinez, A.B., Muis, S., Haigh, I.D., van der Wiel, K., Haarsma, R.J., Ward, P.J., Roberts, M.J., Dullaart, J.C.M., Aerts, J.C.J.H., 2022. A globally consistent local-scale assessment of future tropical cyclone risk. Sci. Adv. 8, eabm8438. https://doi.org/10.1126/sciadv.abm8438 | Compute windfields using CLIMADA Hazard module, with 300 arc seconds resolution      |   |   |   |   |
