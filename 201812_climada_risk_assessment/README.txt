These scripts reproduce the results of paper:
CLIMADA – a global weather and climate risk assessment platform. Aznar-Siguan, G & Bresch, D. N.

Execute first script_exec.py and then script_plots.py from this folder. The
total execution time should not exceed 15 min.

The scripts load the data in the containing folder and write in the same folder
intermediate data as well as the figures. Consecutive executions of the
scripts will be faster because of the loading of intermediate results.

The intermediate results which take longer to execute are already provided 
but will be computed if not present. imp_if_samples and imp_tr_samples
contain each the 100 samples of 34650 tropical cyclones each used in fig07, and
would need hours to be computed (consider using a cluster if you want to
reproduce it).

After the execution, the folder should have a size of 3.5 Gb.
