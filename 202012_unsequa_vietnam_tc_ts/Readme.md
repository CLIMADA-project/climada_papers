# Instructions

These scripts can be used to reproduce the results from the paper 

"Uncertainty and sensitivity analysis for probabilistic weather and climate risk modelling: an implementation in CLIMADA v.3.1.0"

Kropf, C. M., Ciullo, A., Otth, L., Meiler, S., Rana, A., Schmid, E., McCaughey, J. W., and Bresch, D. N.: Geoscientific Model Development (2022)

The hazard data must be downloaded from [here](https://www.research-collection.ethz.ch/handle/20.500.11850/566528) and placed in the folder `VNM_Data`.

To reproduce the uncertainty and sensitivity analysis run the scripts `unsequa_impact.py` and `unsequa_costben.py`. This example is fully self-contained and works with [CLIMADA V.3.1.0](https://zenodo.org/record/5947271) or higher. All data from the paper can be reproduced with these scripts.

The original hazard data, exposures data, and adaptation measures parametrization are taken explained in [Rana, A., Zhu, Q., Detken, A., Whalley, K., and Castet, C.: Strengthening climate-resilient development and transformation in Viet Nam, Climatic Change, 170, 4, 2022.](https://link.springer.com/article/10.1007/s10584-021-03290-y). The original scripts from this paper are available [here](https://github.com/arunranain/climada_tc_vietnam). 
