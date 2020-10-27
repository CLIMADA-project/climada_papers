

1 DAMAGE GENENERATION

The modeling process starts with the country based damage calculation with the impact modeling framework
CLIMADA available under:

https://github.com/CLIMADA-project/climada_python

Installation requirements and instructions are provided in the corresponding documentation.

In order to drive the model, spatially explicit flooded fractions and flood depth provided by the ISIMIP2a
simulation round are required, these data are available under:

https://files.isimip.org/cama-flood/

The script "schedule_runs.py" starts a modeling run for each climate forcing-GHM combination and calls
the simulation script "schedule_sim.py". In "schedule_runs.py" all climate forcing datasets and GHMs
are defined.

The script calculates the flood damage for each country and every year, as decribed under

https://github.com/CLIMADA-project/climada_python/blob/main/doc/tutorial/climada_hazard_RiverFlood.ipynb

The output of "schedule_sim.py" are 46 .csv files containing damage-time series for each country between 1971-2010.


2 POSTPROCESSING

2.1 DATA AGGREGATION

In the first step 