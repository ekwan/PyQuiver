# *PyQuiver*

*A user-friendly program for calculating isotope effects.*

*PyQuiver* is an open-source Python program for calculating kinetic isotope effects (KIEs) and equilibrium isotope effects (EIEs) using harmonic frequencies and the Bigeleisen-Mayer equation.  *PyQuiver* requires Cartesian Hessian matrices, which can be calculated using any electronic structure program.  

### Features

* automatically read frequencies from [`Gaussian`](http://www.gaussian.com/g_prod/g09.htm) and [`ORCA`](https://orcaforum.kofo.mpg.de/app.php/portal) output files
* excellent performance for larger systems
* custom isotopic substitutions and arbitrary temperature
* tunnelling corrections: Wigner and Bell inverted parabola
* run via command line or simple Python API

### Installation

*PyQuiver* requires Python 3 and [`numpy`](https://numpy.org).  No other libraries are necessary.

1. Install [Python](https://www.continuum.io/downloads) if necesary.  The standard Anaconda distribution will contain the necessary dependencies.
2. Install [git](https://git-scm.com/downloads).  git comes pre-installed on most Linux distributions and Macs.
3. Clone the repository: `git clone https://github.com/ekwan/PyQuiver.git`.  (If you don't want to deal with `git`, click on the green "clone or download" button on this github repository page and click on "Download ZIP" to receive an archive.)

*PyQuiver* has been tested on PC, Mac, and Linux platforms.

## Tutorial

To learn how to use *PyQuiver*, please look at the [tutorial](tutorial/TUTORIAL.md).

## Further Reading

* [Tutorial]
* [Technical Details](DETAILS.md) (compatibility with QUIVER, how to invoke PyQuiver, `.config` files, input files, the PyQuiver Standard format, miscellaneous notes)
* [AutoQuiver](src/AUTOQUIVER.md) (how to run PyQuiver on many related files)
* [Snipping Utility](scripts/SNIP.md) (how to "snip" out only the relevant parts of a Gaussian output file for publication or archiving purposes)

## Authors

*PyQuiver* was written by Thayer Anderson and Eugene Kwan.  Please email `ekwan16@gmail.com` with any questions.  We will gladly try to help you.

We thank Gregor Giesen for writing the code to read ORCA output.  We also thank Christoph Riplinger and Giovanni Bistoni for valuable discussions.

## How to Cite

Anderson, T.L.; Kwan, E.E.  *PyQuiver*  **2020**, `www.github.com/ekwan/PyQuiver`

## License
   
This project is licensed under the Apache License, Version 2.0. See `LICENSE.txt` for full terms and conditions.
   
*Copyright 2020 by Thayer L. Anderson and Eugene E. Kwan*
