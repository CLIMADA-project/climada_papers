These scripts reproduce the main results of paper:
CLIMADA v1.4.1: Towards a globally consistent adaptation options appraisal tool; Bresch, D. N. & Aznar-Siguan, G.

First time execution: execute script_risk.py from your Python console. It will take approximately 3 hours to compute 
all the risk metrics of the region and save the results in a "results_eca" folder.
The IBTRACS data used in this publication is contained in IBTrACS.ALL.v04r00.nc. You might use this in your CLIMADA execution 
to obtain exactly the same results as in the publication.

Normal execution: the Jupyter notebook reproduce_results.ipynb loads the risk data computed in the previous step and 
performs the adaptation options analysis, reproducing the figures and results of the paper.
