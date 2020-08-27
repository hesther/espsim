# Comparison of electrostatic potential and shape
This repository contains a small code snippet to calculate similarities of shapes and electrostatic potentials between molecules. It is based on Python3, RDKit, Numpy and Scipy. The package furthermore contains functionalities to embed (create 3D coordinates) molecules with a constrained core using RDKit functions.

## Table of Contents

- [Installation](#installation)
- [Usage](#usage)
  * [ESP similarity](#esp-similarity)
  * [Embedding, alignment and similarity](#embedding-alignment-and-similarity)
  * [Jupyter demo](#jupyter-demo)
- [Contact](#contact)

## Installation

It is easiest to install all required package inside a conda environment. If you want to use jupyter notebooks, and plot 3D representations using py3Dmol, create the following environment (here named `espsim`):

`conda create -n espsim -c conda-forge conda-forge::rdkit numpy scipy jupyter py3dmol`

If you don't need jupyter, you can create an environment with less packages:

`conda create -n espsim -c conda-forge conda-forge::rdkit numpy scipy`

Next, activate the environment:

`conda activate espsim`

(or `source activate espsim` if you have an older version of conda). If you don't want to use conda, install Python3, RDKit, Numpy and Scipy, as well as optionally Jupyter and py3Dmol in any way convenient to you.

Now, install espsim as a pip package as follows, by (1) changing the directory to wherever you saved espsim (use the correct path instead of `<pathtoespsim>`) and compiling the package (2):

1. `cd <pathtoespsim>`
2. `pip install -e .`

Then you can use `import espsim` or from `espsim import ...` in your code. To test your installation, run

`python scripts/test_imports.py`

which should print `Test passed, imports work fine.` or tell you which package/module could not be imported.

## Usage

Within a python script, you can either use only the ESP similarity routine on previously embedded, aligned RDKit molecule objects using the function `GetEspSim()` or use `EmbedAlignConstrainedScore()` to embed, align and score molecules. 

### ESP similarity

An example is provided in `scripts/test_esp_function.py`. Execute it via

`python scripts/test_esp_function.py`

The script loads two RDKit molecules from file, that have been previously embedded and aligned, and calls `GetEspSim()`, which calculates the overlap integrals of the electrostatic potentials of the two molecules. The function returns the ESP similarity. The integration is performed via fitting the `1/r` term in the Coloumb potential by three Gaussians. The overlap integral between the potentials of molecules A and B
then becomes a sum of a large number of two-center integrals of Gaussians for which an analytic solution exists (see Good et al, [doi.org/10.1021/ci00007a002](https://doi.org/10.1021/ci00007a002)). After integration, the overall similarity is calculated as the Carbo similarity (Carbo et al, [doi.org/10.1002/qua.560320412](https://doi.org/10.1002/qua.560320412) and Carbo, Arnau, Intl. J. Quantum Chem 17 (1980) 1185 (no doi available)), which is the overlap integral of A and B divided by the square root of the norms of A and B.

 The function also takes optionally the arguments `prbCharge`  and `refCharge` as input (list or 1D array of partial charges for the probe and reference molecule). If not given, RDKit's Gasteiger charges are computed.

### Embedding, alignment and similarity

If you need to embed and align the the probe and reference molecule, espsim provides a convenient way to do so. The function `EmbedAlignConstrainedScore()` takes a probe molecule, one or more reference molecules, and a core that is to be constrained (with 3D coordinates!) as input, computes constrained embeddings, compares the shape similarities of all combinations and returns both the shape and ESP similarity. The script `scripts/test_embedalignscore.py` uses this function to create 10 conformers of three probe molecules each, and compares their similiarities to 10 conformers of nine reference compounds each:

`python scripts/test_embedalignscore.py`

The function `EmbedAlignConstrainedScore()` takes several optional input arguments as well: `prbNumConfs` and `refNumConfs`, which are the number of conformers to be created for each probe and reference molecule. (The default is 10 each). The more conformers, the more accurate the results, since a better alignment may be found. Thus, the computed shape and ESP similarities are actually lower bounds to the true values (with perfectly aligned conformers). The function also takes `prbCharge` (list or 1D array of partial charges for the probe molecule) and `refCharges` (list of list, or 2D array of partial charges of all reference molecules) as input. If not specified (default), RDKit partial charges are calculated (Gasteiger charges). The option to input partial charges is especially convenient if you have already pre-computed charges, for example from quantum mechanics. The order of the inputted partial charges must be the same as the atoms in the RDKit molecules.

### Jupyter demo

The `scripts` folder also holds a short demo of espsim on Jupyter. To view it, power up a notebook (you must have jupyter and py3Dmol installed, see Installation instructions):

`jupyter notebook`

and click on `scripts` and `short_demonstration.ipynb`. The file holds more information and explanations about espsim, as well as some 3D visualizations of molecule embeddings. The 3D visualizations are especially helpful to view the results of the constrained embedding and alignment. The notebook furthermore holds more detailed information about the integration of the electrostatic potenital.

## Contact
Feel free to post questions, feedback, errors or concerns on github, or email to eheid@mit.edu.