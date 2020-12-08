## Tutorial

In this tutorial, we reproduce the B3LYP/6-31G* KIE predictions for the following Claisen rearrangement reported by Singleton ([reference 5](../DETAILS.md#references)):

<img src="../img/claisen_scheme.png" height=100>

<img src="../img/claisen_table4.png" height=200>

(This is Table 4 in the paper.)

### File Locations

All the files associated with this tutorial are available in the `tutorial/` directory.  There are separate folders for running with Gaussian, ORCA, or PyQuiver standard files.

Here, we'll focus on the use of PyQuiver with Gaussian.  This tutorial requires the *PyQuiver* configuration file `claisen_demo.config` and the g09 output files `claisen_gs.out` and `claisen_ts.out`, representing the ground and transition state frequency calculations, respectively.

If running with ORCA, use the `.hess` files instead of the output files.  If running with Quiver Standard files, use the `.qin` files.

Note that the tutorial results will differ slightly with the provided ORCA files due to subtle differences in geometry.  Thse differences are not experimentally meaningful and can be ignored.

### Getting Started

In general, all KIEs are defined as rate constant(light)/rate constant(heavy).  For example, the absolute KIE at C1 is defined as the rate of the rearrangement with carbon-12 at C1 divided by the rate with carbon-13 at C1. This definition is given by this line of the `claisen_demo.config` file:

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
python quiver.py ../tutorial/gaussian/claisen_demo.config ../tutorial/gaussian/claisen_gs.out ../tutorial/gaussian/claisen_ts.out
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

* locate ground and transition states (if using Gaussian, turn on the verbose `#p` flag)
* run a frequency calculations
* specify the desired isotopic substitutions in a configuration file
* run `python quiver.py` on the configuration, ground state, and transition state files

If EIEs are desired, simply replace the transition state with the equilibrium state of interest.
