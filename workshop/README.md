# Welcome to the RSC CICAG Workshop materials
This page contains all workshop materials from the RSC CICAG Workshop June 23, 2022! If you want to participate directly, register for free [here](https://www.eventbrite.com/e/open-source-tools-for-chemistry-tickets-294585512197?). After the workshop, a full video will be available, too. 

## Materials

To follow along the example, you will either need:

* A Gmail address so you can use Google Colabs (free, no registration, no installation) OR
* A local installation of ESPsim and Jupyter


## Colabs notebook
Click on workshop_colab.ipynb and then on the "Open in Colab" Badge.

## Local installation
You can install ESPsim locally via Conda or Pip. We strongly recommend Conda or other solutions to create local environments instead of changing your system packages!

To create a conda environment and install all necessary packages, run:

* `conda create -n espsim -c rdkit -c conda-forge rdkit numpy scipy matplotlib joblib tqdm pytorch jupyter py3dmol`
* `conda activate espsim`
* `pip install git+https://github.com/hesther/chemprop-atom-bond.git`
* `pip install git+https://github.com/hesther/espsim.git`
* `pip install git+https://github.com/hesther/ehreact.git`

To instead install via pip (preferably after creating a local environment), run:

* `pip install rdkit numpy scipy matplotlib joblib tqdm pytorch jupyter py3dmol`
* `pip install git+https://github.com/hesther/chemprop-atom-bond.git`
* `pip install git+https://github.com/hesther/espsim.git`
* `pip install git+https://github.com/hesther/ehreact.git`

After installing via conda or pip, you can open a Jupyter notebook via

`jupyter notebook`

then open `workshop_jupyter.ipynb` (either download this file or the whole folder first).

