## Technical Details

   - [Compatibility with QUIVER](#compatibility-with-quiver)
   - [How to Invoke PyQuiver](#how-to-invoke-pyquiver)
   - [`.config` Files](#config-files)
   - [Input Files](#input-files)
   - [PyQuiver Standard](#pyquiver-standard)
   - [Miscellaneous Notes](#miscellaneous-notes)
   - [References](#references)

### Compatibility with QUIVER

The development of *PyQuiver* was inspired by the [original](#ref2) Fortran program [QUIVER](https://github.com/ekwan/quiver). *PyQuiver* is designed to be as compatible as possible with the original QUIVER program, but to clarify some ambiguity in choice of masses, configuration files need to be updated for use with *PyQuiver*. See the [Configuration](#config-files) section for detail.

### How to Invoke PyQuiver

*PyQuiver* can be controlled from the command line or its Python API.

Once installed (`pip install pyquiver-kie`), run *PyQuiver* from the command line:

```
pyquiver config_file ground_state_file transition_state_file
```

For more details, run `pyquiver -h`. The help message looks like:

```
usage: pyquiver [-h] [-v] [-s STYLE] [-j JOBS] config gs ts

PyQuiver calculates KIEs and EIEs based on a ground and transition state file.

positional arguments:
  config                configuration file path
  gs                    ground state file path
  ts                    transition state file path

options:
  -h, --help            show this help message and exit
  -v, --verbose         print debug information (repeatable)
  -s STYLE, --style STYLE
                        style of input files (gaussian, orca, or pyquiver)
  -j JOBS, --jobs JOBS  number of worker threads for the isotopologue
                        computations (default 1; -1 uses all cores)
```

This command will calculate the KIEs or EIEs associated with the isotopic substitutions specified in the configuration file. For details, see the tutorial above. (The legacy `python src/quiver.py ...` invocation still works but is deprecated.)

*PyQuiver* also has a Python API that exposes the underlying objects, so you can run custom calculations, scan parameters, and inspect or export results directly:

```python
from pyquiver import Config, KIE_Calculation
config = Config.from_dict(isotopologues={"C1": [(1, 1, "13C")]},
                          temperature=393, scaling=0.9614, imag_threshold=50)
calc = KIE_Calculation(config, "gs.out", "ts.out", style="gaussian")
calc.to_dataframe()              # structured results (pandas)
calc.results["C1"].wigner        # named access
```

See the README for the full API summary (structured results, `Config.from_dict`, the `n_jobs` parameter, and `batch`). The IPython Notebook `pyquiver.ipynb` in `scripts/` demonstrates the API end to end, including the tunnelling corrections (Wigner, Bell, and Skodje-Truhlar; see ref. 4).

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
* `mass_override_isot[pomer/logue]`: possible values are "default" or the name of an isotopologue. If the value "default" is specified, the masses of the light isotopologue will be the defaults found in `weights.dat`. If the name of an isotopologue is given, then that isotopologue is used to replace the default mass behaviour of PyQuiver at a particular atom. For example, if the isotopomer `C2` replaces carbon 2 with 13C, then specifying `mass_override_isotopomer C2` will place carbon-13 at C2 *for every KIE calculation*.
* `isoto[pomer/logue]`: the rule used for isotopic substitution. The expected fields are `name ground_state_atom_number transition_state_atom_number substitution`. The final field, `substitution` must correspond to a valid substitution weight. These weights are specified in `weights.dat` (e.g., `13C`, `18O`, `2D`).

### Input Files

*PyQuiver* assumes that the ground state file and transition state file are the outputs of a Gaussian09 `freq` job. To change this assumption *PyQuiver* can accept an additional command-line argument corresponding to the input file style.

Currently, *PyQuiver* can automatically read output files from the following formats:
* Gaussian. Style name: `gaussian` (aliases: `g16`, `g09`; works with Gaussian 2009 and 2016.)
* ORCA. Style name: `orca` (Reads the `.hess` files generated by any frequency job.)
* PyQuiver Standard. Style name: `pyquiver` (aliases: `native`, `qin`.)

To specify a format other than `gaussian` from the command line, run with the `-s` flag. For instance:

```
pyquiver -s pyquiver tutorial/pyquiver/claisen_demo.config tutorial/pyquiver/claisen_gs.qin tutorial/pyquiver/claisen_ts.qin
```

would execute the example *PyQuiver* job on the Claisen system using the *PyQuiver* standard input files. This allows you to adapt PyQuiver to other electronic structure programs.  (We would be pleased to offer some advice on how to accomplish this.)

### PyQuiver Standard

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

* Example *PyQuiver* standard input files are available in the `tutorial/pyquiver/` directory. The files `claisen_gs.qin` and `claisen_ts.qin` are *PyQuiver* input files corresponding to the example Claisen system discussed in the tutorial.

If input files are provided in a known format other than the *PyQuiver* standard, *PyQuiver* can dump the appropriate *PyQuiver* input files. To do this load the appropriate system (ex. `gs = System("./ground_state.g09")`) and then run `gs.dump_pyquiver_input_file()` which will create the appropriate input file at the same path with the extension `.qin`.

### Miscellaneous Notes

The internal workings of *PyQuiver* are relatively simple.  Essentially, a KIE calculation involves parsing the Hessian, mass-weighting, diagonalization with `np.eigvalsh`, calculation of the reduced isotopic partition functions, substitution into the Bigeleisen-Mayer equation, and tunnelling corrections.  The tunnelling corrections are expected to work well for heavy atoms, but poorly for hydrogen, especially when considering primary KIEs.

The frequencies can optionally be scaled (see ref. 3), but this will probably not help much (the harmonic approximation is usually quite adequate).  Note that *PyQuiver* recognizes linear vs. non-linear molecules and removes the appropriate number of small translational/rotational modes.  (The `frequency_threshold` keyword has been deprecated and is now ignored.) 

The performance of *PyQuiver* is generally excellent, even for large systems.  This is largely because of the efficiency of the `np.eigvalsh` implementation.  When multiple isotopomers are calculated using the same configuration file, the reference isotopomer's frequencies are computed only once and cached (`Isotopologue.calculate_frequencies` short-circuits once `self.frequencies` is set), so the reference is not re-diagonalized for each isotopomer.

There is a known issue with getting PyQuiver to work on Cygwin systems due to a problem processing file paths correctly.  Additionally it seems to be hard to get NumPy installed properly on such systems.  We are working on a fix--please contact me if you require this urgently.  (The program works properly on all other kinds unix/linux, as far as we know).

Please note that the verbose flag should be set in the Gaussian route card.  For example:
` #p b3lyp 6-31g* freq=noraman`
This will ensure that Gaussian places the Hessian in the archive string at the end of the file.  (Not including the verbose flag will cause an error.)

Tunnelling corrections work best for heavy-atom KIEs.  H/D KIEs are more challenging (see the work of the Singleton and Truhlar groups for more details).

### References

1. **Bigeleisen-Mayer theory:**
  * Bigeleisen, J.; Mayer, M.G. *J. Chem. Phys.*  **1947**, *15*, 261.
  * Wolfsberg, M.  *Acc. Chem. Res.*, **1972**, *5*, 225.
  * Wolfsberg, M. *et al.*  <span style="text-decoration:underline">Isotope Effects in the Chemical, Geological, and Bio Sciences</span>  
2. <span id="ref2">**QUIVER:**</span>
  * Saunders, M.; Laidig, K.E.; Wolfsberg, M.  *J. Am. Chem. Soc.*, **1989**, *111*, 8989.
3. **Scaling Factors:**
  * Wong *et al.*  *Chem. Phys. Lett.* **1996**, *256*, 391-399.
  * Radom *et al.*  *J Phys. Chem.* **1996**, *100*, 16502.
4. **Tunnelling Corrections:**
  * Wigner, E.  *Z. Phys. Chem. B*  **1932**, *19*, 203.
  * Bell, R.P.  *Trans. Faraday Soc.*  **1959**, *55*, 1.
  * Bell, R.P.  <span style="text-decoration:underline">The Tunnel Effect in Chemistry</span>; Chapman & Hall: London, **1980**.
  * Bell, R.P.  *Chem. Soc. Rev.*  **1974**, *3*, 513.
  * Skodje, R.T.; Truhlar, D.G. *J. Phys. Chem.* **1981**, *85*, 624-626
5. **Claisen Rearrangement KIEs:**
  * <span id="ref5">Meyer, M.P.; DelMonte, A.J.; Singleton, D.A.</span> *J. Am. Chem. Soc.*, **1999**, *121*, 10865-10874.
