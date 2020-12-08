### Snipping Utility

The AWK script `snip.awk` is provided to convert Gaussian `.out` files into `.snip` files.  These "snips" are contain only the geometry, energies, and frequencies, but retain enough of the Gaussian format for PyQuiver to read them. This feature is useful if you have many files and want to "compress" them to save room.

To convert a collection of `.out` files to their corresponding `.snip` files:

`awk -f snip.awk *.out` 

Note that this will overwrite any existing snips with the same name without warning.  (`awk` is a standard scripting langauge available on any platform.)
