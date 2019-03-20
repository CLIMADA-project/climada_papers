These scripts reproduce the main results of paper:

Samuel Eberenz, Dario Stocker, Thomas Röösli, David N. Bresch (2019, publication in preparation). "Global LitPop: An Exposure Data Model for Disaster Risk Assessment based on Nightlight and Population Data"

Contains 3 independent scripts for Python 3.6 with CLIMADA 1.2.0:

*** litpop_maps.py ***
Compute LitPop, Lit, and Pop for metropolitan areas and plot maps.
(Section 3.1; Figures 2, A1)

*** litpop_validation.py ***
LitPop exposure data model validation + scatter and box plots
Sections 3.2, 3.3; Figures 3, 4; Tables (A1), A2, A3
! slow !

*** litpop_data.py ***
Reproduce LitPop data as in data archive:
LitPop: Global Exposure Data for Disaster Risk Assessment
DOI: 10.3929/ethz-b-000331316
https://www.research-collection.ethz.ch/handle/20.500.11850/331316

-----

Requires Python 3.6 and CLIMADA v1.2.0 (or later):
https://github.com/CLIMADA-project/climada_python/

-----
If script fails, revert CLIMADA version to release v1.2.0 (2018-03):
https://github.com/CLIMADA-project/climada_python/releases/tag/v1.2.0

