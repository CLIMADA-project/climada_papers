# Exposure data for global physical risk assessment

These scripts reproduce the main results of paper:

Eberenz, S., Stocker, D., Röösli, T., and Bresch, D. N.:
Exposure data for global physical risk assessment,
Earth Syst. Sci. Data Discuss., https://doi.org/10.5194/essd-2019-189, in review, 2019.

Contact: [Samuel Eberenz](mailto:samuel.eberenz@usys.ethz.ch)

## Content:
Four independent scripts for Python 3.6 with CLIMADA 1.3.1:

##  litpop_data.py
Reproduces LitPop based gridded asset exposure data as published in teh ETH data archive:
LitPop: Global Exposure Data for Disaster Risk Assessment
DOI: 10.3929/ethz-b-000331316
https://www.research-collection.ethz.ch/handle/20.500.11850/331316
! slow to run for all countries !

##  litpop_world_map.py
Plot global asset exposure data for 224 countries
Sections 3.1
Figure 3
Requires asset exposure data in the form of CSV as available from the ETH research repository:
https://doi.org/10.3929/ethz-b-000331316
Save CSV data in local folder ENTITY_DIR or set path in variable ENTITY_DIR before executing this script.

## litpop_evaluation.py
LitPop exposure data model evaluation for 14 countries and plotting of scatter and box plots.
Sections 2.6, 3.2, 3.3;
Figures 3, 5;
Tables (A1), A2, A3.

! slow to run for all countries !

## litpop_metropolitan_maps.py
Compute Lit^1Pop^1, Lit^1, and Pop^1 for four metropolitan areas and plot maps.
Section 3.3;
Figures 4, A1.

## Requirements

Requires Python 3.6 and CLIMADA v1.2.0 (or later):
https://github.com/CLIMADA-project/climada_python/
Documentation: https://github.com/CLIMADA-project/climada_python/blob/master/doc/source/install.rst

## Documentation:

* Publication: [doi.org/10.5194/essd-2019-189](https://doi.org/10.5194/essd-2019-189)

Documentation for CLIMADA is available on Read the Docs:

* [online (recommended)](https://climada-python.readthedocs.io/en/stable/)
* [PDF file](https://buildmedia.readthedocs.org/media/pdf/climada-python/stable/climada-python.pdf)

Documentation for the LitPop module is available on Read the Docs:

* [online](https://climada-python.readthedocs.io/en/stable/tutorial/climada_entity_LitPop.html)

If script fails, revert CLIMADA version to release v1.2.0 (2018-03):
* from [GitHub](https://github.com/CLIMADA-project/climada_python/releases/tag/v1.2.0)
* from [ETH Data Archive](http://doi.org/10.5905/ethz-1007-226)

...and download GPW v4.10 to the data/systems folder inside CLIMADA
(Download link: http://sedac.ciesin.columbia.edu/downloads/data/gpw-v4/gpw-v4-population-count-rev10/gpw-v4-population-count-rev10_2015_30_sec_tif.zip)

-----
# Update 10 February 2020:
World map added and scripts updated in line with revisions to manuscript.


# Update 25 April 2019:
For the publication and research data, v4.10 of the Gridded Population of the World (GPW) data was used.
This data needs to be downloaded manually from:
http://sedac.ciesin.columbia.edu/downloads/data/gpw-v4/gpw-v4-population-count-rev10/gpw-v4-population-count-rev10_2015_30_sec_tif.zip

Since April 2019, GPW v4.11 is available. CLIMADA is compatible to v4.11 from Release v1.2.3 onwards:
https://github.com/CLIMADA-project/climada_python/releases/tag/v1.2.3

If a GPW version different from v4.10 is used, results can deviate from the original results for certain countries.

-----

www.wcr.ethz.ch
