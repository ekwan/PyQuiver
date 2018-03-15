# *PyQuiver*

*A user-friendly program for calculating isotope effects.*

## Contents
  - [Introduction](#introduction)
    - [Features](#features)
    - [Compatibility with QUIVER](#Compatibility-with-QUIVER)
    - [Installation](#installation)
  
  - [Tutorial](#tutorial)
    - [Summary](#summary)
    - [`autoquiver.py`](#autoquiverpy)

  - [Technical Details](#technical-details)
    - [Interfaces](#interfaces)
    - [.config Files](#config-files)
    - [Input Files](#input-files)
    - [Snipping Utility](#snipping-utility)
    - [Notes](#notes)

  - [Fine Print](#fine-print)
    - [References](#references)
    - [Authors](#authors)
    - [How to Cite](#how-to-cite)
    - [License](#license)
  
## Introduction

*PyQuiver* is an open-source Python program for calculating kinetic isotope effects (KIEs) and equilibrium isotope effects (EIEs) using harmonic frequencies and the Bigeleisen-Mayer equation.  *PyQuiver* requires Cartesian Hessian matrices, which can be calculated using any electronic structure program.  *PyQuiver* works with Gaussian output files by default (Gaussian '09 and '16 are both supported).

### Features

* calculate KIEs or EIEs
* automatically read frequencies from [`Gaussian`](http://www.gaussian.com/g_prod/g09.htm) output files
* excellent performance for larger systems
* highly customizable isotopic substitutions
* arbitrary temperature
* tunnelling corrections: Wigner and Bell inverted parabola
* run via command line or simple Python API

### Compatibility with QUIVER

The development of *PyQuiver* was inspired by the [original](#ref2) Fortran program [QUIVER](https://github.com/ekwan/quiver). *PyQuiver* is designed to be as compatible as possible with the original QUIVER program, but to clarify some ambiguity in choice of masses, configuration files need to be updated for use with *PyQuiver*. See the [Configuration](#config-files) section for detail.

### Installation

*PyQuiver* is written in pure Python 2.7.  It requires a few standard Python libraries like `numpy`.

1. Install [Python](https://www.continuum.io/downloads) if necesary.  The standard Anaconda distribution will contain the necessary dependencies.
2. Install [git](https://git-scm.com/downloads).  git comes pre-installed on most Linux distributions and Macs.
3. Clone the repository: `git clone https://github.com/ekwan/PyQuiver.git`

For those who do not want to deal with `git`, click on the green "clone or download" button on this github repository page and click on "Download ZIP" to receive an archive.

Other than downloading the source code, there is nothing to configure, compile, or fiddle with to get *PyQuiver* to run.  *PyQuiver* has been tested on PC, Mac, and Linux platforms.


## Tutorial

In this tutorial, we reproduce the B3LYP/6-31G* KIE predictions for the following Claisen rearrangement reported by [Singleton](#ref5):

<img src="img/claisen_scheme.png" height=100>

<img src="img/claisen_table4.png" height=200>

(This is Table 4 in the paper.)

All the files associated with this tutorial are available in the `test/` directory. In particular, the tutorial requires the *PyQuiver* configuration file `claisen_demo.config` and the g09 output files `claisen_gs.out` and `claisen_ts.out`, representing the ground and transition state frequency calculations, respectively.

In general, all KIEs are defined as rate(light)/rate(heavy).  For example, the absolute KIE at C1 is defined as the rate of the rearrangement with carbon-12 at C1 divided by the rate with carbon-13 at C1. This definition is given by this line of the `claisen_demo.config` file:

```
isotopomer C1 1 1 13C
```

`C1` is an arbitrary label for the KIE of interest (it can be any string without a space character).  In general, we may want to calculate multiple KIEs using one configuration file.  For example, the next few lines define the KIEs at C2 and the oxygen:

```
isotopomer C2 2 2 13C
isotopomer O3 3 3 17O
```

In each case, the isotopomer definition is followed by three parameters: the atom number in the ground state, the atom number in the transition state, and the isotope to substitute with.  For example, for C1, atom 1 in the ground state and atom 1 in the transition state will be substituted with carbon-13.  In general, *PyQuiver* will try to prevent you from entering isotopomers that do not make sense.

The definition of `13C` is drawn from `src/weights.dat`:

```
carbon,6,C,12.0,12C,12.0,13C,13.00335,14C,14.0031
```

In English, this says that carbon has an atomic number of 6 and has the symbol `C`.  In all cases (`C1`, `C2`, `O`, etc.), whenever carbon appears in the "light" isotopomer, it is defined to have a mass of `12.0`.  This is called the "default mass."  If the "heavy" replacement is specified as `13C`, it is given a mass of `13.00335`.  Although a number of common replacements are already defined (e.g., `2D` (deuterium) or `18O` (oxygen-18)), it is easy to add more definitions.

Note that KIEs can be defined for multiple isotopic replacements by repeating the label of the isotopomer.  For example, these entries replace two hydrogens with two deuteriums:

```
isotopomer H/D 7 7 2D
isotopomer H/D 8 8 2D
```

Now we are ready to calculate the KIEs!  Enter in the following:

```
cd src/
python quiver.py ../test/claisen_demo.config ../test/claisen_gs.out ../test/claisen_ts.out
```

When run from the command line, *PyQuiver* expects the names (in order) of the configuration file, the ground state file, and the transition state file.  The expected output is:

```
Read atomic weight data for 30 elements.

Reading configuration from claisen_demo.config
Reading data from claisen_gs.out... with style g09
Reading data from claisen_ts.out... with style g09
Config file: claisen_demo.config
Temperature: 393.0 K
Scaling: 0.961
Reference Isotopologue: C5
Imag threshold (cm-1): 50
   Isotopologue         C5, replacement  1: replace gs atom  5  and ts atom  5  with 13C
   Isotopologue         C1, replacement  1: replace gs atom  1  and ts atom  1  with 13C
   Isotopologue         C2, replacement  1: replace gs atom  2  and ts atom  2  with 13C
   Isotopologue         C4, replacement  1: replace gs atom  4  and ts atom  4  with 13C
   Isotopologue         C6, replacement  1: replace gs atom  6  and ts atom  6  with 13C
   Isotopologue        H/D, replacement  1: replace gs atom  7  and ts atom  7  with  2D
   Isotopologue        H/D, replacement  2: replace gs atom  8  and ts atom  8  with  2D
   Isotopologue         O3, replacement  1: replace gs atom  3  and ts atom  3  with 17O

=== PyQuiver Analysis ===
Isotopologue                                              uncorrected      Wigner     inverted parabola
                                                              KIE           KIE              KIE
Isotopologue         C1                                      1.0091        1.0107          1.0110
Isotopologue         C2                                      1.0000        1.0000          1.0000
Isotopologue         O3                                      1.0153        1.0168          1.0171
Isotopologue         C4                                      1.0249        1.0276          1.0281
Isotopologue         C6                                      1.0111        1.0128          1.0131
Isotopologue        H/D                                      0.9515        0.9529          0.9532

KIEs referenced to isotopologue C5, whose absolute KIEs are:
Isotopologue         C5                                      1.0019        1.0019          1.0019
```

Note that these KIEs are *relative* to the KIE at `C5`.  This is controlled by this line of the config file:

```
reference_isotopomer C5
```

This means that all absolute KIEs will be divided by this one to give relative KIEs.  Use `none` to calculate absolute KIEs only.

These numbers agree closely with the predictions reported by Singleton.  There are small (0.001) differences that arise from roundoff errors, differing values of physical constants, slight changes in the way masses are handled, and the way the original QUIVER program handles small rotational/translational frequencies with a threshold system.  These small differences should not affect any chemical conclusions.

### Summary

The above captures the basic workflow of a *PyQuiver* calculation:

* locate ground and transition states using g09 (using the verbose `#p` flag)
* run a frequency calculations
* specify the desired isotopic substitutions in a configuration file
* run `python quiver.py` on the configuration, ground state, and transition state files

If EIEs are desired, simply replace the transition state with the equilibrium state of interest.

### `autoquiver.py`

`autoquiver` allows KIEs to be calculated over many ground state and transition state files that share a common set of desired isotopic substitutions.  (The `autoquiver` module has additional functionality when used through a Python interface (read the `src/quiver.ipynb` if interested) but the command line interface should suffice for most purposes.)

Suppose we have a directory, `auto/`, with the following files:

```
substitutions.config
gs-type1.output    ts-type1.output
gs-type2.output    ts-type2.output
gs-type3.output    ts-type3.output
gs-type4.output    ts-type4.output
```

We might want to run *PyQuiver* using the `substitutions.config` file on all pairs of ground states and transition states.  For example, these pairs may be the same calculation run at many levels of theory.  Note that it is crucial that all files have a consistent atom numbering scheme.

To accomplish this we run `autoquiver.py` as follows:

```
python src/autoquiver.py -e .output auto/ auto/substitutions.config gs ts -
```

The arguments are:

* `-e .output`: a flag to look for files with the extension `.output` as the frequency jobs for the ground and transitions states.
* `auto/`: look for files in the `auto/` directory.
* `auto/substitutions.config`: use `auto/substitutions.config` as the configuration file.
* `gs`: use the string "gs" to find ground state files. All files with the appropriate extension that contain the substring "gs" will be treated as ground state files.
* `ts`: use the string "ts" to find transition state files.
* `-`: use the field delimiter "-" to test if a ground state and transition states match. All fields after the first "-" must be identical. This means that `gs-type1.output` and `ts-type1.output` will match but `gs-type1.output` and `ts-type2.output` won't.

The output of autoquiver is directed to a number of csv files corresponding to each configuration file given. These filenames are printed when autoquiver exits.

For more information, the output of `python autoquiver.py -h` has been reproduced below:

```
usage: autoquiver.py [-h] [-v] [-s STYLE] [-e EXT]
                     config target gs_p ts_p delimiter

A program to automatically run PyQuiver on a config file and all ground state
and transition states matching certain constraints.

positional arguments:
  config                configuration file path
  target                target directory file path (where the ground and
                        transition state files live
  gs_p                  substring in ground state files
  ts_p                  substring in transition state files
  delimiter             delimiter used to match ground and transition state
                        files (all fields separated by the delimiter after the
                        first must match)

optional arguments:
  -h, --help            show this help message and exit
  -v, --verbose         when the verbose flag is set debug information is
                        printed
  -s STYLE, --style STYLE
                        style of input files
  -e EXT, --extension EXT
                        extension of input files
```

## Technical Details

   - [Interfaces](#interfaces)
   - [.config Files](#config-files)
   - [Input Files](#input-files)
   - [Snipping Utility](#snipping-utility)
   - [Defining Masses](#masses)
   - [Notes](#notes)

### Interfaces

*PyQuiver* can be controlled from the command line or its Python API.

To run *PyQuiver* from the command line, issue the following command from the `src/` directory:

```
python quiver.py config_file ground_state_file transition_state_file
```

For more details, run `python quiver.py -h` to display the following help message:

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

*PyQuiver* also has a Python API. The advantage of this interface is that it exposes the underlying Python objects. This allows you to run custom calculations, automate routine calculations, and inspect the internal data directly.

An IPython Notebook is provided as a demonstration.  To try it out, move to the `src/` directory and run the command `ipython notebook`. Then open the `quiver.ipynb` notebook file.

### .config Files

Calculations performed in *PyQuiver* require a configuration file to specify the parameters (such as scaling factor and temperature) and isotopologue substitution rules.

Each configuration file is a plain text file with the following properties:

* blank lines and lines starting with `#` are ignored
* anything after `#` within a line is assumed to be a comment
* fields in directives are separated by spaces.

Valid configuration files have all of the following directives:

* `scaling`: a linear factor by which to scale the frequencies
* `imag_threshold`: the threshold (in units cm^-1) that defines the cutoff between small and large imaginary frequencies used to distinguish EIE from KIE calculations (typical value: 50)
* `temperature`: the temperature in Kelvin
* `reference_isoto[pomer/logue]`: possible values are `none` or the name of an isotopologue. If `none` is specified, the absolute KIEs will be reported. If the name of an isotopologue is specified, all KIEs will be divided the KIE values for this isotopologue.
* `mass_override_isot[pomer/logue]`: possible values are "default" or the name of an isotopologue. If the value "default" is specified, the masses of the light isotopologue will be the defaults found in `weights.dat`. If the name of an isotopolgoue is given, then that isotopologue is then used to replace the default mass behaviour of PyQuiver at a particular atom. For example, if the isotopomer `C2` replaces carbon 2 with 13C, then specifying `mass_overide_isotopomer C2` will place carbon-13 at C2 *for every KIE calculation*.
* `isoto[pomer/logue]`: the rule used for isotopic substitution. The expected fields are `name ground_state_atom_number transition_state_atom_number substitution`. The final field, `substitution` must correspond to a valid substitution weight. These weights are specified in `weights.dat` (e.g., `13C`, `18O`, `2D`).

### Input Files

*PyQuiver* assumes that the ground state file and transition state file are the outputs of a Gaussian09 `freq` job. To change this assumption *PyQuiver* can accept an additional command-line argument corresponding to the input file style.

Currently, *PyQuiver* can automatically read output files from the following formats:
* Gaussian 2009. Style name: `g09`
* PyQuiver Standard. Style name: `pyquiver`

To specify a format other than `g09` from the command line, run with the `-s` flag. For instance:

```
python quiver.py -s pyquiver ../test/claisen_demo.config ../test/claisen_gs.qin ../test/claisen_ts.qin
```

would execute the example *PyQuiver* job on the Claisen system using the *PyQuiver* standard input files. This allows you to adapt PyQuiver to other electronic structure programs.  (We would be pleased to offer some advice on how to accomplish this.)

The *PyQuiver* Standard is a generic format for the output of an electronic structure program in plain-text outlined as follows:

* The first line of a file contains the number of atoms (*n*).
* The next *n* lines define the geometry. Each line should be of the form:

```
CenterNumber,AtomicNumber,XPosition,YPosition,ZPosition
```

The positions should be provided in units of Angstroms. The center number is the atom index (ranging from 0 and *n-1*, inclusive).
* The next line should contain the lower triangular Cartesian Hessian matrix with no line breaks. In particular, if *H* is the Hessian matrix then *H<sub>3p+i,3q+j</sub>* corresponds to taking derivatives with respect to atom *p* moving in the *i*-th coordinate and atom *q* moving in the *j*-th coordinate (*i* and *j* run across the three cartesian coordinates). The entries in the serialized form should be separated by commas. The serialization should occur by stringing together the rows (truncated at the main diagonal). For example, suppose the following is the lower-right triangular form of the Cartesian Hessian for a one atom system:

```
1.0
2.0 3.0
4.0 5.0 6.0
```

then the *PyQuiver* would expect the following line:

```
1.0,2.0,3.0,4.0,5.0,6.0
```

* Example *PyQuiver* standard input files are available in the `test/` directory. The files `claisen_gs.qin` and `claisen_ts.qin` are *PyQuiver* input files corresponding to the example Claisen system discussed in the tutorial.

If input files are provided in a known format other than the *PyQuiver* standard, *PyQuiver* can dump the appropriate *PyQuiver* input files. To do this load the appropriate system (ex. `gs = System("./ground_state.g09")`) and then run `gs.dump_pyquiver_input_file()` which will create the appropriate input file at the same path with the extension `.qin`.

### Snipping Utility

The AWK script `scripts/snip.awk` is provided to convert g09 `.out` files into `.snip` files.  These "snips" are contain only the geometry, energies, and frequencies, but retain enough of the g09 format for PyQuiver to read them. This feature is useful if you have many files and want to "compress" them to save room.

To convert a collection of `.out` files to their corresponding `.snip` files:

`awk -f snip.awk *.out` 

Note that this will overwrite any existing snips with the same name without warning.  (`awk` is a standard scripting langauge available on any platform.)

### Notes

The internal workings of *PyQuiver* are relatively simple.  Essentially, a KIE calculation involves parsing the Hessian, mass-weighting, diagonalization with `np.eigvalsh`, calculation of the reduced isotopic partition functions, substitution into the Bigeleisen-Mayer equation, and tunnelling corrections.  The tunnelling corrections are expected to work well for heavy atoms, but poorly for hydrogen, especially when considering primary KIEs.

The frequencies can optionally be scaled (see ref. 3), but this will probably not help much (the harmonic approximation is usually quite adequate).  Note that *PyQuiver* is recognizees linear vs. non-linear molecules and removes the appropriate number of small translational/rotational modes.  (The `frequency_threshold` keyword has been deprecated and is now ignored.) 

The performance of *PyQuiver* is generally excellent, even for large systems.  This is largely because of the efficiency of the `np.eighvalsh` implementation.  Note that when multiple isotopomers are calculated using the same configuration file, *PyQuiver* will recalculate the frequencies for the reference isotopomer repeatedly (i.e., once for every isotopomer).  This should not be relevant for routine use.  However, it could be avoided by using the *PyQuiver* API.

There is a known issue with getting PyQuiver to work on Cygwin systems due to a problem processing file paths correctly.  Additionally it seems to be hard to get NumPy installed properly on such systems.  We are working on a fix--please contact me if you require this urgently.  (The program works properly on all other kinds unix/linux, as far as we know).

Please note that the verbose flag should be set in the Gaussian route card.  For example:
` #p b3lyp 6-31g* freq=noraman`
This will ensure that Gaussian places the Hessian in the archive string at the end of the file.  (Not including the verbose flag will cause an error.)

Tunnelling corrections work best for heavy-atom KIEs.  H/D KIEs are more challenging (see the work of the Singleton and Truhlar groups for more details).

# Fine Print

## References

1. **Bigeleisen-Mayer theory:**
  * Bigeleisen, J.; Mayer, M.G. *J. Chem. Phys.*  **1947**, *15*, 261.
  * Wolfsberg, M.  *Acc. Chem. Res.*, **1972**, *5*, 225.
  * Wolfsberg, M. *et al.*  <span style="text-decoration:underline">Isotope Effects in the Chemical, Geological, and Bio Sciences</span>  
2. <span id="ref2">**QUIVER:**</span>
  * Saunders, M.; Laidig, K.E. Wolfsberg, M.  *J. Am. Chem. Soc.*, **1988**, *111*, 8989.
3. **Scaling Factors:**
  * Wong *et al.*  *Chem. Phys. Lett.* **1996**, *256*, 391-399.
  * Radom *et al.*  *J Phys. Chem.* **1996**, *100*, 16502.
4. **Tunnelling Corrections:**
  * Bell, R.P.  *Chem. Soc. Rev.*  **1974**, *3*, 513.
5. **Claisen Rearragenent KIEs:**
  * <span id="ref5">Meyer, M.P.; DelMonte, A.J.; Singleton, D.A.</span> *J. Am. Chem. Soc.*, **1999**, *121*, 10865-10874.

## Authors

*PyQuiver* was written by Thayer Anderson and Eugene Kwan at the Department of Chemistry and Chemical Biology at Harvard University.  Please email `ekwan@fas.harvard.edu` with any questions.  We will gladly try to help you.

## How to Cite

Anderson, T.L.; Kwan, E.E.  *PyQuiver*  **2016**, www.github.com/ekwan/PyQuiver.

## License
   
This project is licensed under the Apache License, Version 2.0. See `LICENSE.txt` for full terms and conditions.
   
*Copyright 2016 by Thayer L. Anderson and Eugene E. Kwan*
