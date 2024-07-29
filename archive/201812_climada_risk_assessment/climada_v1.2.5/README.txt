These scripts reproduce the main results of paper:
Aznar-Siguan, G. and Bresch, D. N.: CLIMADA v1: a global weather and climate risk assessment platform, Geosci. Model Dev., 12, 3085-3097, https://doi.org/10.5194/gmd-12-3085-2019, 2019.

Execute script_exec.py with CLIMADA version v1.2.5. 
The execution in a laptop should take less than 2 hours.

The script writes in the same folder intermediate data as well as the figures.
Consecutive executions of the scripts will be faster thanks to the loading
of intermediate results.

After the execution, the folder should have a size of 900 Mb. The results might not be exactly
the sames as the published ones, since IBTrACS is periodically correcting their data. In contrast
to older versions, this script downloads the current IBTrACS data.
