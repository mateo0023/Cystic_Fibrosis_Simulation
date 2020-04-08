# What's this?

This repository hosts a model for airway epithelial cells. The model was developed with the Python programming language.

The model is heavily based on the one done by [Sandefour et. al.](https://www.ncbi.nlm.nih.gov/pubmed/28808008), with a few modifications, currently the cotransporter being used is the one proposed by Benjamin-Franklin in 1997.

# Libraries in Use

```python
multiprocessing.pool.Pool
pandas
numpy
scipy
datetime
zipfile
sys
os
```

### In `simulation.py`

```
pandas
numpy
scipy
```



# `initial-values.csv`

This file contains the values to which the simulations will be initially set.

#### Required

* Amount of Steps (`'max_steps'`) an integer value
* Time-frame (`'time_frame'`) is how much time each step represents. Since we're approximating the functions, the `time_frame` is the *precision* of the approximation. The smaller the number the better it might be. **Note that the total simulated time will be `time_frame` times `max_steps` in seconds**
* Whether it is Cystic Fibrosis or not (`'CF'`).
* ASL Height (`H`) 
* Apical and Cellular concentrations: `'aNa', 'aCl', 'aK', 'cNa', 'cCl', 'cK'`
* Nucleotide Concentrations: `'ATP', 'ADP', 'AMP', 'ADO', 'INO'`

#### Optional Values

* Membrane potentials: `'aV', 'bV'`. If not given they will be calculated by `AirwayModel.fn_voltage`

#### What about the other variables?

They will be calculated based on the given values.

### Units for the `initial-values.csv`

* `H`: meters.
* `time_frame`: seconds
* Ionic Concentrations (`aNa` -> `cCl`): mM ($10^{-3}$ Molar).
* Nucleotides (`ATP` -> `INO`): microM ($\mu M$ or $10^{-6}$ Molar).
* Voltage (`aV` and `bV`) volts.
* Cystic Fibrosis (`CF`): 1 or *'CF'* (modeling cystic fibrosis) **or** 0 or *'NL'* (normal conditions).

# `simulation.py`

This file contains two classes, `Constants` and `AirwayModel`. It is everything that the model needs to work, it is designed in this way for easier exports to other programs.

## Constants

Is basically a big dictionary for the `AirwayModel` to save it's parameter values.

The `Constants` class now calculates the parameters based on the following steady state values:

| $\frac{dN^C_{Na}}{dt} = 0$ | $\frac{dN^C_{K}}{dt}=0$ | $\frac{dN^C_{Cl}}{dt}=0$ |
| :------------------------: | :---------------------: | :----------------------: |
|    $\frac{dV_a}{dt}=0$     |                         |   $\frac{dV_b}{dt}=0$    |

The units are the same as in the paper [Sandefour et. al.](https://www.ncbi.nlm.nih.gov/pubmed/28808008)

## AirwayModel

This is the *big boy*, it has all of the functions of the model. It is meant to be run with the `AirwayModel.fn_runAll()` function

# Intended Use

## Running `main.py`

If the file `initial-values.csv` has more than one line with values, it will run the different simulations in parallel of up two 8 simulations at a time. It will automatically save files after the run is done.

## Importing `simulation.py`

1. Create an instance of a class:
   * The `init_data` can be a Dictionary, path to a csv file containing the values or left blank (user will be prompted to input the values manually)
   * `special_constants` a dictionary with the constants that will have a different value than on the paper. **Recommended only when testing**. The general naming is `where_what_ofWhat` and *concentration* is *"CONS"* because it has the same sound as it would. However, checking the file `simulation.py` for all the corresponding variable names is more recommended and probably necessary.
3. Run the function `AirwayModel.fn_runAll()`. It will complete the simulations from steps `1` through `max_steps -1`.
4. Export the data to a CSV file. This step has been made easier by declaring `self.to_csv = self.data.to_csv`  under `AirwayModel.__init__` . That means that you're calling the `to_csv` function on the `Pandas.DataFrame` object  that is the simulation.

```python
import simulation as sim

# The csv file is recomended to follow the template of the one in the repository.
simulated_data = sim.AirwayModel('initial-values.csv')
simulated_data.fn_runAll()

simulated_data.to_csv('Exported-File.csv')
```



***Read the documentation in the files for more information***