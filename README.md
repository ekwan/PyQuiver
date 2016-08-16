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

Picture and table references from paper for claisen reaction

This tutorial will walk through an example KIE calculation while explaining the components in broad terms. More detailed and sophisticated techniques are exposed in the form of the underlying Python objects, but the simple command line interface should suffice for most use cases.

All the files associated with this tutorial are available in the `test/` directory. In particular, the tutorial will reference the configuration file `claisen_demo.config` and the g09 output files `claisen_gs.out` and `claisen_ts.out` representing the ground and transition state frequency calculations, respectively.

All KIEs (and EIEs) refer to an isotopic substitution made in both the ground and transition state (or the exchanging systems for an EIE calculation). In the `claisen_demo.config` file, one line reads:

```
isotopomer C1 1 1 13C
```

This line provides the details of isotopic substitution necessary for *PyQuiver* to make the replacements. The keyword `isotopomer` lets *PyQuiver* know that the line will describe an isotopic substitution. The label `C1` defines the name of the isotopomer (substitution). This name has no syntactic significance - it can be any string without a space character. The numbers `1` and `1` tell *PyQuiver* which atom to replace in the ground and transition state, respectively. The final entry, `13C`, defines the weight of the isotope used in the substitution. This weight must be drawn from `src/weights.dat` file. Examples include `13C` (for Carbon-13), `2D` (for Deuterium), and `18O` for (Oxygen-18). An isotopomer need not contain only a single replacement. KIEs can be calculated for systems where the heavy isotopomer has substitutions at multiple atoms. To make a multiple replacement, the label of the isotopomer is simply repeated (as seen at the end of the example configuration file):

```
isotopomer H/D 7 7 2D
isotopomer H/D 8 8 2D
```

This defines an isotopomer named `H/D` that replaces hydrogens 7 and 8 in the ground and transition with deuterium.

Once the substitutions are specified in the configuration file, *PyQuiver* will read in the cartesian Hessian/second derivative matrix/force constant matrix to calculate the appropriate KIE. The Bigeleisen-Mayer method for KIE calculation relates the KIEs to the normal modes of vibration of the molecule. In particular, the frequencies of the normal modes are used to calculate the reduced isotopic partition functions for the ground and transition state, which are then divided to find the KIE. 

*PyQuiver* automates this procedure. To run a KIE calculation for the example system, move to the `src/` directory and run `quiver.py` from the command line:
```
cd src/
python quiver.py ../test/claisen_demo.config ../test/claisen_gs.out ../test/claisen_ts.out
```
This command accepts (in order) the configuration file, the ground state file, and the transition state file. When run, the command will print a summary of the configuration file used (including all isotopic substitutions) and then calculated and print the KIEs corresponding to each isotopomer. For each KIE, three numbers are printed. These numbers correspond to the uncorrected KIE and two tunneling-corrected KIEs.

The expected output of the claisen test case is as follows:

```
Read atomic weight data for 30 elements.

Reading configuration from ./claisen_demo.config
Reading data from claisen_gs.out... with style g09
14
Reading data from claisen_ts.out... with style g09
14
Config file: ./claisen_demo.config
Temperature: 393.0 K
Scaling: 0.961
Reference Isotopologue: reference
Frequency threshold (cm-1): 50
   Isotopologue  reference, replacement  1: replace gs atom  5  and ts atom  5  with 13C
   Isotopologue         C1, replacement  1: replace gs atom  1  and ts atom  1  with 13C
   Isotopologue         C2, replacement  1: replace gs atom  2  and ts atom  2  with 13C
   Isotopologue         C4, replacement  1: replace gs atom  4  and ts atom  4  with 13C
   Isotopologue         C5, replacement  1: replace gs atom  5  and ts atom  5  with 13C
   Isotopologue         C6, replacement  1: replace gs atom  6  and ts atom  6  with 13C
   Isotopologue        H/D, replacement  1: replace gs atom  7  and ts atom  7  with  2D
   Isotopologue        H/D, replacement  2: replace gs atom  8  and ts atom  8  with  2D
   Isotopologue         O3, replacement  1: replace gs atom  3  and ts atom  3  with 17O

=== PY-QUIVER ANALYSIS ===
Isotopologue                                              uncorrected      Widmer     infinite parabola
                                                              KIE           KIE              KIE
Isotopologue         C1                                      1.011         1.012            1.013      
Isotopologue         C2                                      1.000         1.000            1.000      
Isotopologue         C4                                      1.028         1.031            1.031      
Isotopologue         C5                                      1.000         1.000            1.000      
Isotopologue         C6                                      1.013         1.015            1.015      
Isotopologue        H/D                                      0.953         0.954            0.955      
Isotopologue         O3                                      1.017         1.018            1.019  
```


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


*PyQuiver* can be controlled from the command line or from an [IPython Notebook](https://ipython.org/notebook.html).

To run *PyQuiver* from the command line, simply move to the `src/` directory and input the following command:

```
python quiver.py config_file ground_state_file transition_state_file
```

The command line interface accepts some standard and some custom flags for usage: `quiver.py [-h] [-v] [-s STYLE] config gs ts`. To see more details, run `python quiver.py -h` to display the following help message:

```
usage: quiver.py [-h] [-v] [-s STYLE] config gs ts

A program that calculates KIEs and EIEs based on a ground and transition state
file.

positional arguments:
  config                configuration file path
  gs                    ground state file path
  ts                    transition state file path

optional arguments:
  -h, --help            show this help message and exit
  -v, --verbose         when the verbose flag is set debug information is
                        printed
  -s STYLE, --style STYLE
                        style of input files
```


This command will calculate the KIEs or EIEs associated with the isotopic substitutions specified in the configuration file. For details, see the tutorial above.

*PyQuiver* also offers an IPython Notebook interface. This advantage of this interface is that it exposes the underlying Python objects. This allows you to run custom calculations, automate routine calculations, and inspect the internal data directly.

To use the IPython Notebook, move to the `src/` directory and run the command `ipython notebook`. Then open the `quiver.ipynb` notebook file. To run a calculation replace the arguments for the `KIE_Calculation()` object as described in detail in the notebook.

### .config Files

Calculations performed in *PyQuiver* require a configuration file to specify the parameters (such as scaling factor and temperature) and isotopologue substitution rules.

Each configuration file is a plain text file with the following properties:
* blank lines and lines starting with `#` are ignored. All other lines must correspond to valid directives. (Furthermore `#` begins a comment within a line).
* fields in directives are separated by spaces.

Valid configuration files have all of the following directives:
* `scaling`: a linear factor by which to scale the frequencies. Recommended value: 1.0
* `frequency_threshold`: the threshold (in units cm^-1) that defines the cutoff between the small frequencies (corresponding to translation and rotation) and the normal vibrational mode frequencies. Typical value: 50. (Tests show that projecting out rotations and translations have no effect on the KIE).
* `temperature`: the temperature in Kelvin at which to model the calculation.
* `reference_isoto[pomer/logue]`: possible values are "default" or the name of an isotopologue. If "default" is specified, the absolute KIEs will be reported. If the name of an isotopologue is specified, all KIEs will be divided the KIE values for this isotopologue.
* `mass_override_isot[pomer/logue]`: possible values are "default" or the name of an isotopologue. If the value "default" is specified, the default masses in `weights.dat`. If the name of an isotopolgoue is given, then that isotopologue is used to replace the default mass behaviour of PyQuiver. In particular, those substitutions are made in both the ground and transition state of both the heavy and light isotopologues for all KIE calculations. For example, if for some reason you wish for Carbon 5 to have the mass of 12.5, for whatever reason, you would specify such an isotopologue as the `mass_override_isotopologue`.
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
* Example *PyQuiver* standard input files are available in the `test/` directory. The files `claisen_gs.qin` and `claisen_ts.qin` are *PyQuiver* input files corresponding to the example Claisen system discussed in the tutorial.

If input files are provided in a known format other than the *PyQuiver* standard, *PyQuiver* can dump the appropriate *PyQuiver* input files. To do this load the appropriate system (ex. `gs = System("./ground_state.g09")`) and then run `gs.dump_pyquiver_input_file()` which will create the appropriate input file at the same path with the extension `.qin`.

### Defining Masses

### Math

The math behind a quiver calculation is detailed in the PDF `doc/technical_details.pdf` generated from the TeX file `doc/technical_details.tex`.

## References

1. **Bigeleisen-Mayer theory:**
  * J. Chem. Phys.*  **1947**, *15*, 261.
  * Wolfsberg, M.  *Acc. Chem. Res.* **1972**, *5*, 225.
2. **QUIVER:**
  * Saunders, M.; Laidig, K.E. Wolfsberg, M.  *JACS* **19898*, *111*, 8989.
3. **Scaling Factors:**
  * Wong.  *Chem. Phys. Lett.* **1996**, *256*, 391-399.
  * Radom.  *J Phys. Chem.* **1996**, *100*, 16502.
4. **Tunnelling Corrections:**
  * Bell.  *Chem. Soc. Rev.*  **1974**, *3*, 513.

## Authors
*PyQuiver* was written by Thayer Anderson and Eugene Kwan at the Department of Chemistry and Chemical Biology at Harvard University.  Please email `ekwan@fas.harvard.edu` with any questions.

## License

