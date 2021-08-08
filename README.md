
# Worm-GEM: The generic genome-scale metabolic model of _Caenorhabditis elegans_

[![GitHub version](https://badge.fury.io/gh/sysbiochalmers%2FWorm-GEM.svg)](https://badge.fury.io/gh/sysbiochalmers%2FWorm-GEM)


#### Description

This repository contains the latest version of Worm-GEM, a worm genome-scale metabolic model.


#### Citation

 > H. Wang, J. L. Robinson, P. Kocabaş, J. Gustafsson, M. Anton, P.-E. Cholley, et al. Genome-scale metabolic network reconstruction of model animals as a platform for translational research. _PNAS_ 118, e2102344118 (2021). [doi.org/10.1073/pnas.2102344118](https://doi.org/10.1073/pnas.2102344118)

#### Model Keywords

**Utilisation:** multi-omics integrative analysis, predictive simulation
**Field:** metabolic-network reconstruction  
**Type of Model:** reconstruction, curated  
**Model source:** [Human-GEM](https://doi.org/10.1126/scisignal.aaz1482)   
**Omic source:** genomics; metabolomics   
**Taxonomic name:** _Caenorhabditis elegans_  
**Taxonomy ID:** [6239](https://identifiers.org/taxonomy:6239)  
**Genome ID:** [GCA_000002985.3](https://identifiers.org/insdc.gca:GCA_000002985.3)  
**Metabolic System:** general metabolism  
**Tissue:**  
**Bioreactor:**    
**Cell type:**  
**Cell line:**  
**Condition:** generic metabolism


### Model Overview

|Taxonomy | Template Model | Reactions | Metabolites| Genes |
| ------------- |:-------------:|:-------------:|:-------------:|:-----:|
|_Caenorhabditis elegans_ |   Human-GEM |  12169  | 8158 | 1946 |


## Installation

### Required Software
* A functional MATLAB installation (MATLAB 7.3 and higher).
* The [RAVEN toolbox](https://github.com/SysBioChalmers/RAVEN).
* The [COBRA toolbox](https://github.com/opencobra/cobratoolbox) (not necessary for most functionality).


### Dependencies - Recommended Software
* The libSBML MATLAB API (version [5.13.0](https://sourceforge.net/projects/sbml/files/libsbml/5.13.0/stable/MATLAB%20interface/) is recommended).
* [Gurobi Optimizer](http://www.gurobi.com/registration/download-reg) for any simulations.


### Installation Instructions
* Clone the [master branch](https://github.com/SysBioChalmers/Worm-GEM/tree/master) of this repository, or [download the latest release](https://github.com/SysBioChalmers/Worm-GEM/releases/latest).
* Add the directory to your MATLAB path (instructions [here](https://se.mathworks.com/help/matlab/ref/addpath.html?requestedDomain=www.mathworks.com)).


## Usage

#### Loading/saving the model

`Worm-GEM.mat` (Recommended if on `master` branch)
* Load and save using the built-in MATLAB `load()` and `save()` functions.

`Worm-GEM.xml` (SBML format)
* Load using the `importModel.m` function (from [RAVEN Toolbox](https://github.com/SysBioChalmers/RAVEN))
* Save using the `exportModel.m` function (from [RAVEN Toolbox](https://github.com/SysBioChalmers/RAVEN))


## Websites

- [Metabolic Atlas](https://metabolicatlas.org/) enables visualization and exploration of Worm-GEM content.


### Contributing

Contributions are always welcome! Please read the [contributing guidelines](.github/CONTRIBUTING.md) to get started.

