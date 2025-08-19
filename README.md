# Comparison of electrostatic potential and shape
This repository contains a small code snippet to calculate similarities of shapes and electrostatic potentials between molecules, see [manuscript](https://doi.org/10.1021/acs.jcim.1c01535). It is based on Python3, RDKit, Numpy and Scipy. The package furthermore contains functionalities to embed (create 3D coordinates) molecules with a constrained core using RDKit functions.

## New

* On June 23, 2022 RSC-CICAG hosted an in-depth ESPsim workshop. The recording is available on [Youtube](https://www.youtube.com/watch?v=Ka08REoGYvI), and the workshop materials are available in the workshop folder. We hope this resource is valuable to you!
* December 2021: ESPsim has exciting new features: A machine-learned partial charge model, as well as aligning and scoring of molecules without a common core. We furthermore added various benchmarks, see jupyter notebooks in the benchmarks folder.

## Table of Contents

- [Installation](#installation)
- [Usage](#usage)
  * [ESP similarity](#esp-similarity)
  * [Further options](#further-options)
  * [Embedding, alignment and similarity](#embedding-alignment-and-similarity)
  * [Embedding without constrained core](#embedding-without-constrained-core)
  * [Jupyter demo](#jupyter-demo)
- [Contact](#contact)

## Installation

It is easiest to install all required packages inside a conda environment. Run

`conda env create -f environment.yml`

to install all necessary packages. Some users have reported installation problems, in that case it can be easier to recreate an existing environment: `espsim_fixed_env.yml` instead of the generic `environment.yml`. Next, activate the environment:

`conda activate espsim`

(or `source activate espsim` if you have an older version of conda). If you don't want to use conda, install Python3, RDKit, Numpy, Scipy, Scikit-Learn, Pytorch, Tqdm, as well as optionally Jupyter, py3Dmol, Psi4, RESP (for calculation of quantum mechanical partial charges), and https://github.com/hesther/chemprop-atom-bond.git (for calculation of machine-learned partial charges) in any way convenient to you.

Now, install espsim as a pip package as follows:

`pip install -e .`

Then you can use `import espsim` or from `espsim import ...` in your code. To test your installation, run

`python scripts/test_imports.py`

which should print `Test passed, imports work fine.` or tell you which package/module could not be imported.

## Usage

Within a python script, you can either use only the ESP similarity routine on previously embedded, aligned RDKit molecule objects using the function `GetEspSim()` or use `EmbedAlignConstrainedScore()` to embed, align and score molecules with a common core, or `EmbedAlignScore()` to embed, align and score molecules without a common core. 

### ESP similarity

An example is provided in `scripts/test_esp_function.py`. Execute it via

`python scripts/test_esp_function.py`

The script loads two RDKit molecules from file (mol, mol2 and sdf with optional custom charges in mol2 and sdf), that have been previously embedded and aligned, and calls `GetEspSim()`, which calculates the overlap integrals of the electrostatic potentials of the two molecules. The function returns the ESP similarity (target: 0.85). The integration is performed via fitting the `1/r` term in the Coloumb potential by three Gaussians. The overlap integral between the potentials of molecules A and B
then becomes a sum of a large number of two-center integrals of Gaussians for which an analytic solution exists (see Good et al, [doi.org/10.1021/ci00007a002](https://doi.org/10.1021/ci00007a002)). After integration, the overall similarity is calculated as the Carbo similarity (Carbo et al, [doi.org/10.1002/qua.560320412](https://doi.org/10.1002/qua.560320412) and Carbo, Arnau, Intl. J. Quantum Chem 17 (1980) 1185 (no doi available)), which is the overlap integral of A and B divided by the square root of the norms of A and B.

 The function also takes optionally the arguments `prbCharge`  and `refCharge` as input (list or 1D array of partial charges for the probe and reference molecule). If not given, RDKit's Gasteiger charges are computed.

### Further options

ESPsim provides a range of options to customize the behavior of `GetEspSim()`:

* `prbCharge`  and `refCharge` to use custom partial charges.
* `metric` determines the similarity metric used. Defaults to `metric = "carbo"` (Carbo similarity). Currently implemented: `carbo` with a range of -1 to 1 and `tanimoto` with a range of -1/3 to 1 (Tanimoto similarity, which is the overlap integral of A and B divided by the sum of the norms of A and B minus the overlap integral).
* `renormalize` is a Boolean which determines whether the obtained similarities should be rescaled to the range `[0:1]` (or to a custom range, given by `customrange`)
* `integrate` determines the integration routine. Defaults to `integrate = "gauss"` (analytic integration via Gaussian functions). Currently implemented: `gauss` and `mc` (Monte Carlo numeric integration). The main difference between the two is that integration via Gaussian functions encompasses integration over all space (including within van-der-Waals radii of atoms, and locations far away from the molecule), whereas the MC integration is local (up to a specified margin around the molecules), and excludes space within the van-der-Waals radii of atoms. The margin of the MC integration, i.e. the minimum distance between the van-der-Waals surface and the integration box can be specified via `marginMC` in Angstrom, and defaults to 10 Angstrom. Smaller integration margins yield a similarity comparison for the immediate vicinity of the molecule. The number of points of the integration per 1 cubic Angstrom can be set by `nMC`, and defaults to 1. Smaller numbers, such as 0.1 or 0.01 speed up the calculation but decrease the accuracy of the integration.
* `partialCharges` determines the partial charge distribution to use. Defaults to `partialCharges = "gasteiger"`. Currently implemented: `gasteiger`, `mmff` (Merk Molecular Force Field MMFF94), `ml` (Machine learning model adapted from [this article](https://doi.org/10.1039/D0SC04823B) - caution, only trained on neutral molecules!), or `resp` (quantum-mechanical determination of partial charges via Psi4 and fitting to the electrostatic potential under constraints (RESP)). Note that `resp` is much slower than using Gasteiger, MMFF94 or machine-learned charges. The method and basis set of the quantum mechanical calculation can be customized via `methodPsi4` (default `"scf"`, that is, Hartree-Fock) and `basisPsi4` (default `"3-21G"`). For available options, refer to Psi4. The number of electrostatic potential grid points within Psi4 can be regulated via `gridPsi4` (default `1`, that is 1 grid point per cubic Angstrom).


### Embedding, alignment and similarity

If you need to embed and align the the probe and reference molecule, espsim provides a convenient way to do so. The function `EmbedAlignConstrainedScore()` takes a probe molecule, one or more reference molecules, and a core that is to be constrained (with 3D coordinates!) as input, computes constrained embeddings, compares the shape similarities of all combinations and returns both the shape and ESP similarity. The script `scripts/test_embedalignscore.py` uses this function to create 10 conformers of three probe molecules each, and compares their similiarities to 10 conformers of nine reference compounds each:

`python scripts/test_embedalignscore.py`

The function `EmbedAlignConstrainedScore()` takes several optional input arguments as well: `prbNumConfs` and `refNumConfs`, which are the number of conformers to be created for each probe and reference molecule. (The default is 10 each). The more conformers, the more accurate the results, since a better alignment may be found. Thus, the computed shape and ESP similarities are actually lower bounds to the true values (with perfectly aligned conformers). The function also takes `prbCharge` (list or 1D array of partial charges for the probe molecule) and `refCharges` (list of list, or 2D array of partial charges of all reference molecules) as input. If not specified (default), RDKit partial charges are calculated (Gasteiger charges). The option to input partial charges is especially convenient if you have already pre-computed charges, for example from quantum mechanics. The order of the inputted partial charges must be the same as the atoms in the RDKit molecules. Further options (partial charge distribution, metric, renormalization, integration routine) can be specified analogous to the `GetEspSim()` function. Finally, you may choose to calculate ESP similarities at every alignment (not only the best alignment based on shape), which may be achieved via `getBestESP = True` (caution, this takes longer to compute).

### Embedding without constrained core

(Caution: This feature is currently experimental and not well tested yet). For molecules without common core, the function `EmbedAlignScore()` can be used instead of `EmbedAlignConstrainedScore()` (with the same input arguments except the `core` argument). Here, identificiation of the best alignment via shape scores only often leads to poor ESP similarities, so that it is often necessary to use `getBestESP = True`.

### Jupyter demo

The `scripts` folder also holds a short demo of espsim on Jupyter. To view it, power up a notebook (you must have jupyter and py3Dmol installed, see Installation instructions):

`jupyter notebook`

and click on `scripts` and `short_demonstration.ipynb`. The file holds more information and explanations about espsim, as well as some 3D visualizations of molecule embeddings. The 3D visualizations are especially helpful to view the results of the constrained embedding and alignment. The notebook furthermore holds more detailed information about the integration of the electrostatic potenital.

### Benchmarks

Visit the `benchmarks` folder to show scripts to benchmark ESP-Sim on DUD-E and [D4-rescore](https://github.com/ljmartin/d4-rescore), compare performance across different partial charges, as well as against EON.


## Contact
Feel free to post questions, feedback, errors or concerns on github, or email to esther.heid@tuwien.ac.at
