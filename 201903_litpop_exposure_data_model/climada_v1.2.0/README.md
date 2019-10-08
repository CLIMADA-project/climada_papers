# Exposure data for global physical risk assessment

These scripts reproduce the main results of paper:

Samuel Eberenz, Dario Stocker, Thomas Röösli, David N. Bresch (2019, in prep.). "Exposure data for global physical risk assessment"

Contact: [Samuel Eberenz](mailto:samuel.eberenz@usys.ethz.ch)

## Content:
three independent scripts for Python 3.6 with CLIMADA 1.2.0:

## litpop_maps.py
Compute LitPop, Lit, and Pop for metropolitan areas and plot maps.
(Section 3.1; Figures 2, A1)

## litpop_validation.py
LitPop exposure data model validation + scatter and box plots
Sections 3.2, 3.3; Figures 3, 4; Tables (A1), A2, A3, S1

! slow to run for all countries !

##  litpop_data.py
Reproduce LitPop data as in data archive:
LitPop: Global Exposure Data for Disaster Risk Assessment
DOI: 10.3929/ethz-b-000331316
https://www.research-collection.ethz.ch/handle/20.500.11850/331316

## Requirements

Requires Python 3.6 and CLIMADA v1.2.0 (or later):
https://github.com/CLIMADA-project/climada_python/
Documentation: https://github.com/CLIMADA-project/climada_python/blob/master/doc/source/install.rst

## Documentation:

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

# Update 25 April 2019:
For the publication and research data, v4.10 of the Gridded Population of the World (GPW) data was used.
This data needs to be downloaded manually from:
http://sedac.ciesin.columbia.edu/downloads/data/gpw-v4/gpw-v4-population-count-rev10/gpw-v4-population-count-rev10_2015_30_sec_tif.zip

Since April 2019, GPW v4.11 is available. CLIMADA is compatible to v4.11 from Release v1.2.3 onwards:
https://github.com/CLIMADA-project/climada_python/releases/tag/v1.2.3

If a GPW version different from v4.10 is used, results can deviate from the original results for certain countries.

-----

www.wcr.ethz.ch

