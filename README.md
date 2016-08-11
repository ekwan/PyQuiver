# *PyQuiver*

## Contents
  - [Introduction](#introduction)
    - [Installation](#installation)
	- [Features](#features)
  - [Tutorial](#tutorial)
    - [Interfaces](#interfaces)
	- [.config Files](#.config-files)
	- [Input Files](#input-files)
  
## Introduction

*PyQuiver* is an open-source Python package for calculating Kinetic Isotope Effects (KIEs) and Equilibrium Isotope Effects (EIEs) based on the output of electronic structure programs. *PyQuiver* reads the Hessian matrix and geometry of a chemical system to calculate partition function ratios and thereby determine isotope effects according to established mathematical methods.

### Installation

This package is written entirely in Python 2.7 and depends only on `numpy` - a standard Python package for scientific computing. If you do not already have Python installed, you should download the appropriate [Python installer](https://www.python.org/download/releases/2.7/). To install `numpy`, you should use [`pip install numpy`](https://pip.pypa.io/en/stable/installing/).

To install *PyQuiver* simply `git clone` the repository.

### Features

* parse [`g09`](http://www.gaussian.com/g_prod/g09.htm) frequency job output files.
* calculate frequencies from the mass-weighted Hessian.
* calculate Reduced Partition Function Ratios (RPFRs) for isotopically substituted molecules.
* calculate KIEs and EIEs for isotopically substituted chemical systems (including the appropriate temperature dependent effects).
* correct KIEs for tunneling using Wigner tunneling corrections and Bell infinite parabola corrections.

## Tutorial
*PyQuiver* supports many types of inputs and offers two interfaces for general use: the command line and an [IPython Notebook](https://ipython.org/notebook.html).

### Interfaces

To run *PyQuiver* from the command line, simply move to the `src/` directory and input the following command:

`python quiver.py config_file ground_state_file transition_state_file`

This command will calculate the KIEs or EIEs associated with the isotopic substitutions specified in the configuration file.

To use the IPython Notebook interface, move to the `src/` directory and run the command `ipython notebook`. Then open the `quiver.ipynb` notebook file. To run a calculation replace the arguments for the `KIE_Calculation()` object as described in detail in the notebook.

### .config Files

Calculations performed in *PyQuiver* require a configuration file to specify the parameters (such as scaling factor and temperature) and isotopologue substitution rules.

Each configuration file is a plain text file with the following properties:
* blank lines and lines starting with "#" are ignored. All other lines must correspond to valid directives
* fields in directives are separated by spaces.

Valid configuration files have all of the following directives:
* `scaling`: a linear factor by which to scale the frequencies. Recommended value: 1.0
* `frequency_threshold`: the threshold (in units cm^-1) that defines the cutoff between the small frequencies (corresponding to translation and rotation) and the normal vibrational mode frequencies. Recommended value: 50
* `temperature`: the temperature in Kelvin at which to model the calculation.
* `reference_isoto[pomer/logue]`: the name of an isotologue to use as the reference for KIE calculations. The name "default" is specially reserved for use here. If the name "default" is specified, then the KIEs are not referenced and are calculated absolutely.
* `isoto[pomer/logue]`: the rule used for isotopic substitution. The expected fields are `name ground_state_atom_number transition_state_atom_number substitution`. The final field, `substitution` must correspond to a valid substitution weight. These weights are specified in `weights.dat`. Examples include `13C`, `18O`, and `2D`.

### Input Files

*PyQuiver* defaults to the assumption that the ground state file and transition state file are the outputs of a Gaussian09 `freq` job. To change this assumption *PyQuiver* can accept an additional command-line argument corresponding to the input file style.

Currently, *PyQuiver* can automatically read output files from the following electronic structure programs:
* Gaussian 2009. Style name: `g09`
* PyQuiver Standard. Style name: `pyquiver`

If you require *PyQuiver* to support any other program, simply email the authors with an example frequency job output associated with that program, and the parsing will be implemented as soon as possible.

The PyQuiver Standard is outlined as follows:


## Technical Details

## Authors
*PyQuiver* was written by Thayer Anderson and Eugene Kwan in the Department of Chemistry and Chemical Biology at Harvard University.

## License

