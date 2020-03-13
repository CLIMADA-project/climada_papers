# Comparing an insurer’s perspective on building damages with modelled damages from pan-European winter windstorm event sets: a case study from Zurich, Switzerland

These scripts reproduce the main results of the paper:

Christoph Welker (1), Thomas Röösli (2,3), David N. Bresch (2,3):
Comparing an insurer’s perspective on building damages with modelled damages from pan-European winter windstorm event sets: a case study from Zurich, Switzerland,
submitted to the journal Natural Hazards and Earth System Sciences (NHESS)

(1) GVZ Gebäudeversicherung Kanton Zürich, Zurich, Switzerland

(2) Institute for Environmental Decisions, ETH Zurich, Zurich, Switzerland

(3) Federal Office of Meteorology and Climatology MeteoSwiss, Zurich, Switzerland

Contact: [Thomas Röösli](mailto:thomas.roeoesli@usys.ethz.ch)

## Content:
one jupyter notebook and one independent script for Python 3.6 with CLIMADA 1.4.0:

##  WISC_GVZ_Zurich_Analysis.ipynb
damage model and risk assessment
Fig 1, 2, 3, 4
Tables 1 and A1

##  probabilistic_extension.py
Storm Severity Index and General Extremal Value Distribution
Selection of parameters for the creation of probabilistic events
Sections 2.2.3



## Requirements

requires the download of the "historic storm footprints" as .netcdf files from the webpage https://wisc.climate.copernicus.eu/wisc/#/help/products#footprint_section

Requires Python 3.6 and CLIMADA v1.4.0 (or later):
https://github.com/CLIMADA-project/climada_python/
Documentation: https://github.com/CLIMADA-project/climada_python/blob/master/doc/source/install.rst

## Documentation:

* Publication: now submitted (will be linked here as soon as available)

Documentation for CLIMADA is available on Read the Docs:

* [online (recommended)](https://climada-python.readthedocs.io/en/stable/)
* [PDF file](https://buildmedia.readthedocs.org/media/pdf/climada-python/stable/climada-python.pdf)

Documentation for the StromEurope module is available on Read the Docs:

* [online](https://climada-python.readthedocs.io/en/stable/tutorial/climada_hazard_StormEurope.html)

If script fails, revert CLIMADA version to release v1.4.0 (2020-03):
* from [GitHub](https://github.com/CLIMADA-project/climada_python/releases/tag/v1.4.0)
* from [ETH Data Archive](http://doi.org/10.5905/ethz-1007-226)



-----
## History
created 20 February 2020

-----

www.wcr.ethz.ch
