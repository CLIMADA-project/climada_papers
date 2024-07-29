
README 

-----------------------------------------------------------------------------------------------------------

SYSTEM REQUIREMENTS

Python (3.6+) version of CLIMADA release v1.5.1

For pre- and post-processing we also recommend Python 3

Tested on Python 3.7

Required non-standard python packages for post-processing are:

pyts.decomposition.SingularSpectrumAnalysis
pymannkendall
statsmodels
itertools
astropy.convolution

For the demo, we strongly recommend to install Jupyter-Notebook

2 INSTALLATION GUIDE

A detailed description on how to install CLIMADA is provided under

https://github.com/CLIMADA-project/climada_python

Typical installation time including all tests (~1.5h)

Post-processing (Python 3) can be done with any Python environment. 

3 DEMO

For damage generation with CLIMADA please see the RiverFlood Tutorial

https://github.com/CLIMADA-project/climada_python/blob/main/doc/tutorial/climada_hazard_RiverFlood.ipynb

For post-processing please see the DEMO_Scripts and follow the tutorial under 

https://github.com/CLIMADA-project/climada_papers/tree/master/202010_flood_attribution/Demo

Please note that only dummies are provided for observational data, as we have no rights to publish the data_sets.
Starting with the jupyter-notebook 'DemoDataAggregation.ipynb'

The output generated is an example for gaining the results presented in Fig.2-4 and Fig. SI1 and SI2
for the example of Latin America. Please note that results are only partly real,
as we only present dummies instead of real observed damages.
Outputs are small dataframes and plots for Latin America.

4 INSTRUCTIONS FOR USE
For the use of the demo data just start the jupyter-notebook 'demo_data_aggregation.ipynb' in DEMO_scripts and
follow the instructions. Only the input data set provided in DEMO_data is needed. Further input is generated.

------------------------------------------------------------------------------------------------------------
0 PREPROCESSING

The modeling of spatially explicit flood depth and fraction with CaMa-Flood and additional post-processing
can be accessed under



1 DAMAGE GENENERATION

The modeling process starts with the country based damage calculation with the impact modeling framework
CLIMADA available under:

https://github.com/CLIMADA-project/climada_python

Installation requirements and instructions are provided in the corresponding documentation.

In order to run the model, spatially explicit flooded fractions and flood depth provided by the ISIMIP2a
simulation round are required, these data are available under:

https://files.isimip.org/cama-flood/

The script "schedule_runs.py" starts a modeling run for each climate forcing-GHM combination and calls
the simulation script "schedule_sim.py". In "schedule_runs.py" all climate forcing datasets and GHMs
are defined.

The script calculates the flood damage for each country and every year, as decribed under

https://github.com/CLIMADA-project/climada_python/blob/main/doc/tutorial/climada_hazard_RiverFlood.ipynb

The output of "schedule_sim.py" are 46 .csv files containing damage-time series for each country between 1971-2010.


2 POST-PROCESSING

The entire post-processing analysis is done once on regional level and on subregional level. 
Scripts ending with '...regions.py' are used to analyse entire regions, while scripts ending with
'...subregions.py' are for the analysis of subregions. Datasets derived from scripts with the Ending
'...regions.py' have to be used as an input for Scripts with the Ending '...regions.py', similarly 
scripts with the Ending'...Subregions.py' have to be used as an input for Scripts with the Ending
'...Subregions.py'

2.1 DATA AGGREGATION

In the first step, data is aggregated to regional/subregional level and across all model-runs, so damage time-series
aggregated to model-medians for each region/subregion are the output. Aditionally, observed damages and country specific 
indicators are added and aggregated. For the aggregation the output files of the 
'schedule_sim.py' script need to be accessed by both scripts: data_aggregation_regions.py and data_aggregation_subegions.py

2.2 VULNERABILITY ASSESSMENT

The aggregated files are then used for the vulnerability assessment in 'vulnerability_adjustment_regions.py' 
and 'vulnerability_adjustment_subregions.py' , further input is not necessary. The scripts each provide a MetaData and a
TimeSeries dataset which are then used for the attribution scripts. The MetaData contains information on explained variances 
and correlation.

2.3 ATTRIBUTION ASSESSMENT

The TimeSeries output is then used as an input for the scripts 'attributionRegions.py' and 'attributionSubregions.py'.
The Scripts again produce TimeSeries and MetaData which can than used to produce the Plots 2,3 and 4.
Both data sets serve as input for the detection of teleconnections. MetaData contains climate (H), exposure (E) and 
vulnerability (V) contributions as well as modeled (I) and observed trend (N) as well as significances.

2.4 DRIVERS FOR CLIMATE-INDUCED TRENDS

The two data sets generated during the attribution assessment are used as inputs for the scripts 'teleconnections_regions.py'
and 'teleconnections_subregions.py'. Climate Oscillation Indices need to be added and are available under:
Southern Oscillation Index as a predictor for ENSO (https://www.ncdc.noaa.gov/teleconnections/enso/enso-tech.php)
Monthly data for AMO, NAO and PDO were extracted from the NOAA/Climate Prediction Center
(https://www.psl.noaa.gov/data/climateindices/list/).



3 PLOTTING

The plot scripts are named according to their Figures in the papers. Which datasets are needed to produce
the plot is indicated in the script. 
