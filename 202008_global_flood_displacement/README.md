# Global warming and population change both heighten future risk of human displacement due to river floods

These scripts reproduce the main results of the paper:

Pui Man Kam(1), Gabriela Aznar Siguan(2), Jacob Schewe(3), Leonardo Milano(4), Justin Ginetti(5), Sven Willner(3), Jamie McCaughey(1), and David N. Bresch(1,2):
Global warming and population change both heighten future risk of human displacement due to river floods, submitted to Environmental Research Letters (ERL)

Publication status: submitted

(1) Institute for Environmental Decisions, ETH Zurich, Zurich, Switzerland

(2) Federal Office of Meteorology and Climatology MeteoSwiss, Zurich, Switzerland

(3) Potsdam Institute for Climate Impact Research, Potsdam, Germany

(4) Centre for Humanitarian Data, United Nations OCHA, The Hague, The Netherlands

(5) Internal Displacement Monitoring Centre, Geneva, Switzerland

Contact: [Pui Man Kam](mannie.kam@usys.ethz.ch)

## Content:

#### fl_cc_displacement_raster.py
Python script to compute the estimated absolute number of displacement in the climate change scenarios.
The output raster files contain the estimated number of displacement at 5km resolution globally in a particular decade, with each band (total 21) represent a combination of general circulation models (GCMs) and global hydrological models (GHMs).

#### fl_historical_displacement_raster.py
Same as fl_cc_displacement_raster.py, but compute the historical simulation from 1966-2005 (baseline period) using the constant population at base year 2000 and the historical simulation from the flood models.

## Requirements
Requires:
* Python 3.6+ environment (best to use conda for CLIMADA repository)
* _CLIMADA_ repository version 1.4.1+:
        https://wcr.ethz.ch/research/climada.html
        https://github.com/CLIMADA-project/climada_python
* flood hazard data from _ISIMIP2b_ (Frieler et al., 2017; data available online or ask the authors for access):
        https://www.isimip.org/gettingstarted/data-access/
* Spatial Population Scenarios data for SSP1 and SSP4 (Jones and O'Neil, 2016; Gao, 2017; data available from NCAR or ask the authors for access):
        http://www.cgd.ucar.edu/iam/modeling/spatial-population-scenarios.html

## Documentation:
Publication: submitted

Documentation for CLIMADA is available on Read the Docs:
* [online (recommended)](https://climada-python.readthedocs.io/en/stable/)
* [PDF file](https://buildmedia.readthedocs.org/media/pdf/climada-python/stable/climada-python.pdf)

If script fails, revert CLIMADA version to release v1.4.1:
* from [GitHub](https://github.com/CLIMADA-project/climada_python/releases/tag/v1.4.1)
* or from the [ETH Data Archive](http://doi.org/10.5905/ethz-1007-252)

## Reference
Frieler, K., Lange, S., Piontek, F., Reyer, C. P. O., Schewe, J., Warszawski, L., Zhao, F., Chini, L., Denvil, S., Emanuel, K., Geiger, T., Halladay, K., Hurtt, G., Mengel, M., Murakami, D., Ostberg, S., Popp, A., Riva, R., Stevanovic, M., Suzuki, T., Volkholz, J., Burke, E., Ciais, P., Ebi, K., Eddy, T. D., Elliott, J., Galbraith, E., Gosling, S. N., Hattermann, F., Hickler, T., Hinkel, J., Hof, C., Huber, V., Jägermeyr, J., Krysanova, V., Marcé, R., Müller Schmied, H., Mouratiadou, I., Pierson, D., Tittensor, D. P., Vautard, R., van Vliet, M., Biber, M. F., Betts, R. A., Bodirsky, B. L., Deryng, D., Frolking, S., Jones, C. D., Lotze, H. K., Lotze-Campen, H., Sahajpal, R., Thonicke, K., Tian, H., and Yamagata, Y.: Assessing the impacts of 1.5 °C global warming – simulation protocol of the Inter-Sectoral Impact Model Intercomparison Project (ISIMIP2b), Geosci. Model Dev., 10, 4321–4345, https://doi.org/10.5194/gmd-10-4321-2017, 2017.

Gao, J., 2017. Downscaling Global Spatial Population Projections from 1/8-degree to 1-km Grid Cells. NCAR Technical Note NCAR/TN-537+STR, DOI: 10.5065/D60Z721H.

Jones, B., O’Neill, B.C., 2016. Spatially explicit global population scenarios consistent with the Shared Socioeconomic Pathways. Environmental Research Letters 11, 84003. DOI:10.1088/1748-9326/11/8/084003.

## History

Created on 24 August 2020

-----

www.wcr.ethz.ch
