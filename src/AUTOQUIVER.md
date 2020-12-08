### `autoquiver.py`

`autoquiver` allows KIEs to be calculated over many ground state and transition state files that share a common set of desired isotopic substitutions.  (The `autoquiver` module has additional functionality when used through a Python interface (read the `scripts/pyquiver.ipynb` if interested) but the command line interface should suffice for most purposes.)

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
