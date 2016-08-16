# *PyQuiver*

## Contents
  - [Introduction](#introduction)
    - [Features](#features)
    - [Installation](#installation)

  - [Tutorial](#tutorial)

  - [Technical Details](#technical-details)
   - [Interfaces](#interfaces)
   - [.config Files](#config-files)
   - [Input Files](#input-files)
   - [Defining Masses](#masses)
   - [Math](#math)
   - [Notes](#notes)

  - [Fine Print](#fine-print)
   - [References](#references)
   - [Authors](#authors)
   - [License](#license)
  
## Introduction

*PyQuiver* is an open-source Python program for calculating kinetic isotope effects (KIEs) and equilibrium isotope effects (EIEs) using harmonic frequencies and the Bigeleisen-Mayer equation.  *PyQuiver* requires the Cartesian Hessian matrix, which can be calculated using any electronic structure program. 

### Features

* calculate KIEs or EIEs
* automatically read frequencies from [`g09`](http://www.gaussian.com/g_prod/g09.htm) output files
* excellent performance for larger systems
* highly customizable isotopic substitutions
* arbitrary temperature
* tunnelling corrections: Wigner and Bell infinite parabola
* run via command line or simple Python API

The development of *PyQuiver* was inspired by the orignal fortran program [QUIVER](https://github.com/ekwan/quiver), written by Keith Laidig. *PyQuiver* is designed to be as compatible as possible with the QUIVER program, but to clarify some ambiguity in choice of masses, configuration files need to updated for use with *PyQuiver*. See the [Configuration](#config-files) section for detail.

### Installation

*PyQuiver* is written in pure Python 2.7.  Its only dependency is  `numpy`, a standard Python package for scientific computing that is included in virtually every Python distribution.

1. Install [Python](https://www.continuum.io/downloads) if necesary.
2. Install `numpy` if necessary: [`pip install numpy`](https://pip.pypa.io/en/stable/installing/)
3. Install [git](https://git-scm.com/downloads).  git comes pre-installed on most Linux distributions and Macs.
4. Clone the repository: `git clone https://github.com/ekwan/PyQuiver.git`

For those who do not want to deal with git, click on the green "clone or download" button on this github repository page and click on "Download ZIP" to receive an archive.

Other than downloading the source code, there is nothing to configure, compile, or fiddle with to get *PyQuiver* to run.


## Tutorial

This tutorial will walk through an example KIE calculation while explaining the components in broad terms. More detailed and sophisticated techniques are exposed in the form of the underlying Python objects, but the simple command line interface should suffice for most use cases.

All the files associated with this tutorial are available in the `test/` directory. In particular, the tutorial will reference the configuration file `claisen_demo.config` and the g09 output files `claisen_gs.out` and `claisen_ts.out` representing the ground and transition state frequency calculations, respectively.

All KIEs (and EIEs) refer to an isotopic substitution made in both the ground and transition state (or the exchanging systems for an EIE calculation). In the `claisen_demo.config` file, one line reads:

```isotopomer C1 1 1 13C```

This line provides the details of isotopic substitution necessary for *PyQuiver* to make the replacements. The keyword `isotopomer` lets *PyQuiver* know that the line will describe an isotopic substitution. The label `C1` defines the name of the isotopomer (substitution). This name has no syntactic significance - it can be any string without a space character. The numbers `1` and `1` tell *PyQuiver* which atom to replace in the ground and transition state, respectively. The final entry, `13C`, defines the weight of the isotope used in the substitution. This weight must be drawn from `src/weights.dat` file. Examples include `13C` (for Carbon-13), `2D` (for Deuterium), and `18O` for (Oxygen-18). An isotopomer need not contain only a single replacement. KIEs can be calculated for systems where the heavy isotopomer has substitutions at multiple atoms. To make a multiple replacement, the label of the isotopomer is simply repeated (as seen at the end of the example configuration file):

```isotopomer H/D 7 7 2D
isotopomer H/D 8 8 2D```

This defines an isotopomer named `H/D` that replaces hydrogens 7 and 8 in the ground and transition with deuterium.

Once the substitutions are specified in the configuration file, *PyQuiver* will read in the cartesian Hessian/second derivative matrix/force constant matrix to calculate the appropriate KIE. The Bigeleisen-Mayer method for KIE calculation relates the KIEs to the normal modes of vibration of the molecule. In particular, the frequencies of the normal modes are used to calculate partition functions for the ground and transition state, which are then divided to find the KIE. 

*PyQuiver* automates this procedure. To run a KIE calculation for the example system, move to the `src/` directory and run `quiver.py` from the command line:
```cd src/```
`python quiver.py ../test/claisen_demo.config ../test/claisen_gs.out ../test/claisen_ts.out`
This command accepts (in order) the configuration file, the ground state file, and the transition state file. When run, the command will print a summary of the configuration file used (including all isotopic substitutions) and then calculated and print the KIEs corresponding to each isotopomer. For each KIE, three numbers are printed. These numbers correspond to the uncorrected KIE and two tunneling-corrected KIEs. In general, the infinite parabola KIE is most accurate.

### Summary

The above captures the basic workflow of a *PyQuiver* calculation:
* run a frequency calculation and collect the output files.
* write a configuration file specifying the desired isotopic substitutions. (Note that there are other necessary parts of a configuration file, see the [appropriate section](#config-files) for details).
* run `python quiver.py` on the configuration, ground state, and transition state files to print the KIEs.

## Technical Details

   - [Interfaces](#interfaces)
   - [.config Files](#config-files)
   - [Input Files](#input-files)
   - [Defining Masses](#masses)
   - [Math](#math)
   - [Notes](#notes)


### Interfaces


*PyQuiver* supports many types of inputs and offers two interfaces for general use: the command line and an [IPython Notebook](https://ipython.org/notebook.html).

To run *PyQuiver* from the command line, simply move to the `src/` directory and input the following command:

```python quiver.py config_file ground_state_file transition_state_file```

The command line interface accepts some standard and some custom flags for usage: `quiver.py [-h] [-v] [-s STYLE] config gs ts`. To see more details, run `python quiver.py -h` to display a help message.


This command will calculate the KIEs or EIEs associated with the isotopic substitutions specified in the configuration file. For details, see the tutorial above.

*PyQuiver* also offers an IPython Notebook interface. This advantage of this interface is that it exposes the underlying Python objects. This allows you to run custom calculations, automate routine calculations, and inspect the internal data directly.

To use the IPython Notebook, move to the `src/` directory and run the command `ipython notebook`. Then open the `quiver.ipynb` notebook file. To run a calculation replace the arguments for the `KIE_Calculation()` object as described in detail in the notebook.

### .config Files

Calculations performed in *PyQuiver* require a configuration file to specify the parameters (such as scaling factor and temperature) and isotopologue substitution rules.

Each configuration file is a plain text file with the following properties:
* blank lines and lines starting with "#" are ignored. All other lines must correspond to valid directives
* fields in directives are separated by spaces.

Valid configuration files have all of the following directives:
* `scaling`: a linear factor by which to scale the frequencies. Recommended value: 1.0
* `frequency_threshold`: the threshold (in units cm^-1) that defines the cutoff between the small frequencies (corresponding to translation and rotation) and the normal vibrational mode frequencies. Recommended value: 50
* `temperature`: the temperature in Kelvin at which to model the calculation.
* `reference_isoto[pomer/logue]`: the name of an isotopologue to use as the reference for KIE calculations. The name "default" is specially reserved for use here. If the name "default" is specified, then the KIEs are not referenced and are calculated absolutely.
* `unsubstitued_isot[pomer/logue]`: the name of an isotopologue to replace the default mass behaviour of PyQuiver. If for some reason you wish to override the default mass assignments made by PyQuiver, then the specified isotopologue will be used for masses in both the heavy and light isotopomers in addition to any replacements made to calculate KIEs. If the name "default" is specified, then the masses are assigned using the default behaviour.
* `isoto[pomer/logue]`: the rule used for isotopic substitution. The expected fields are `name ground_state_atom_number transition_state_atom_number substitution`. The final field, `substitution` must correspond to a valid substitution weight. These weights are specified in `weights.dat`. Examples include `13C`, `18O`, and `2D`.

### Input Files

*PyQuiver* defaults to the assumption that the ground state file and transition state file are the outputs of a Gaussian09 `freq` job. To change this assumption *PyQuiver* can accept an additional command-line argument corresponding to the input file style.

Currently, *PyQuiver* can automatically read output files from the following electronic structure programs:
* Gaussian 2009. Style name: `g09`
* PyQuiver Standard. Style name: `pyquiver`

If you require *PyQuiver* to support any other program, simply email the authors with an example frequency job output associated with that program, and the parsing will be implemented as soon as possible.

The *PyQuiver* Standard is outlined as follows:
* Plain text input files.
* The first line of a ground-state or transition-state input file should have the form `NumberOfAtoms`.
* The next `n` lines, where `n` is the number of atoms specified in the first line define the geometry. Each line should be of the form `CenterNumber,AtomicNumber,XPosition,YPosition,ZPosition`. The positions should be provided in units of Angstroms. The center number simply refers to a numbered label of the atom ranging between `0` and `n-1` (inclusive).
* The next line should contain the serialized lower-right triangular Cartesian Hessian matrix defined in the usual fashion. In particular, if `H` is the Hessian matrix then `H_(3p+i,3q+j)` corresponds to taking derivatives with respect to atom `p` moving in the `i`th coordinate and atom `q` moving in the `j`th coordinate (`i` and `j` run across the three cartesian coordinates). The entries in the serialized form should be separated by commas. The serialization should occur by stringing together the rows (truncated at the main diagonal).

If input files are provided in a known format other than the *PyQuiver* standard, *PyQuiver* can dump the appropriate *PyQuiver* input files. To do this load the appropriate system (ex. `gs = System("./ground_state.g09")`) and then run `gs.dump_pyquiver_input_file()` which will create the appropriate input file at the same path with the extension `.pyq`.


## References

## Authors
*PyQuiver* was written by Thayer Anderson and Eugene Kwan at the Department of Chemistry and Chemical Biology at Harvard University.  Please email `ekwan@fas.harvard.edu` with any questions.

## License

