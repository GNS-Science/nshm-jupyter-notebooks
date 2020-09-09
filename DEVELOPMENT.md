# nshm-jupyter-notebooks Development

This project uses Jupyter notebooks. This note decribes setting thie project on for local development or just to run it locally.

There are no tests in this repo as it's a demonstration project. We rely on all the upstream projects and their testing regime

## Pre-requisites

 * **anaconda** >= 4.8 is required so that this project is compatible with the mybinder.org service.
    The bulk of this instruction is about setting up your conda environment for ther project
 * **git** and **github** account is needed for pulling in other https://github.com/GNS-Science projects via **environment.yml**

### Setup conda environment

```
git clone https://github.com/GNS-Science/nshm-jupyter-notebooks
cd nshm-jupyter-notebooks
conda env create
```

## Setup Jupyter notebook kernel:

```
conda activate nshm-jupyter-notebooks
python -m ipykernel install --user --name nshm-jupyter-notebooks
```

## Start up Jupyter

```
jupyter notebook
```
