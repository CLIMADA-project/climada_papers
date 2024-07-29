# Instructions

These scripts can be used to reproduce the results from the paper

"Global multi-hazard risk assessment in a changing climate"
https://doi.org/10.31223/X56W9N
Zélie Stalhandske, Carmen B. Steinmann, Simona Meiler, Inga Sauer, Thomas Vogt, David N. Bresch, Chahan M. Kropf

Contact: zelie.stalhandske@usys.ethz.ch

CLIMADA version v4.0.2

The project is separated in three main folders. The folder named python_scripts contains that scripts 
that are needed to reproduce the results, bash_scripts the scripts that were used to run the 
different python scripts on the ETH euler cluster. These are numbered in the order that they were run.
The notebooks folder contains the scripts to generate the figures used in the paper. 

The following data is used in these scripts:
1. ISIMIP flooded area and fraction from ISIMIP2b available at: https://zenodo.org/record/4627841#.Ysb_iHhBw5k
2. The TC track simulations are available for scientific purposes only and upon request from WindRiskTech (info@windrisktech.com). 
3. The heat stress hazard can be calculated based on the ESI formula using the “hurs”, “tas”, and “rsds” variables from the ISIMIP2b data https://data.isimip.org
4. LitPop population and exposure data, are available through the CLIMADA API as shown in the scripts.