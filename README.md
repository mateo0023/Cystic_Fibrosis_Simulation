# What's this?

This repository hosts a model for airway epithelial cells. The model was developed with the Python programming language.

The model is heavily based on the one done by [Sandefour et. al.](https://www.ncbi.nlm.nih.gov/pubmed/28808008), with a few modifications, currently the cotransporter being used is the one proposed by Benjamin-Franklin in 1997.

# What is this branch?

This leaves only the `AirwayModel` class, removes everything else from the `.py` files. It's better when using the model for other aplications such as running the simulations with the MATLAB programming language.

# Units for the `initial-values.csv`

* `H`: meters.
* Concentrations (`aNa` -> `cCl`): mM ($10^{-3}$ Molar).
* Nucleotides (`ATP` -> `INO`): microM ($\mu M$ or $10^{-6}$ Molar).
* Cystic Fibrosis (`CF`): 1 (modeling cystic fibrosis) or 0 (normal conditions).
