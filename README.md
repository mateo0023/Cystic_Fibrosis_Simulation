# What's this?

This repository hosts a model for airway epithelial cells. The model was developed with the Python programming language.

The model is heavily based on the one done by [Sandefour et. al.](https://www.ncbi.nlm.nih.gov/pubmed/28808008), with a few modifications, currently the cotransporter being used is the one proposed by Benjamin-Franklin in 1997.

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
