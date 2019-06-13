These files are the ones I'm using with `Python 3.4.8` at Michigan's computers.

# What's different with these files?

* It runs many in parallel, loaded from `initial-values.csv`.
  * How many in parallel it's determined by the variable MAX_PROCESSES in the document `simulation.py`.
  * If there's an error with the file it will ask the user (you) to manually input the necessary values.
* It saves the output into a sub-directory called "Data".
  * It can create the folder if it does not exist.
  * Before starting a new run, it will compress the contents on that folder and delete the folder.
* Now most of the initial values are calculated off of a few (to make sure there are less inconsistencies with the data).
* It now saves two CSV files, the extra one with the `runtime`, initial values and average values.
  * If there's an error in the middle of the run, it will abort it and save the current progress. The _extra data_ file will contain relevant info on the error.
* The simulations support the runs of both Normal conditions and Cystic Fibrosis.
* The data files are compressed into a .ZIP file in case they need to be downloaded one-by one.

# Units for the `initial-values.csv`

* `H`: meters.
* Concentrations (`aNa` -> `cCl`): mM ($10^{-3}$ Molar).
* Nucleotides (`ATP` -> `INO`): microM ($\mu M$ or $10^{-6}$ Molar).
* Cystic Fibrosis (`CF`): 1 (modeling cystic fibrosis) or 0 (normal conditions).
  * It now *has an effect on the simulations!*
