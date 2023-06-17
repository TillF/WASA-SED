!!! This repository contains the source code. For generating the executable, please use a Fortran compiler of your system.

!!! Windows-users may try the (irregularly updated) binaries provided in the [release-section](https://github.com/TillF/WASA-SED/releases).


WASA-SED Model
--------------
WASA-SED is a numerical model for simulating hydrological and sediment fluxes from meso-scale catchments [(Mueller et al., 2010)](https://www.geosci-model-dev.net/3/275/2010/gmd-3-275-2010.pdf).
This repository contains source code and documentation.

All rights under the stated license (see [```license.txt```](https://github.com/TillF/WASA-SED/blob/master/license.txt)) of the WASA-SED code are with the SESAM-Project and its successors, c/o Till Francke, Universität Potsdam, Karl-Liebknecht-Str. 24-25, 14473 Potsdam, Germany.
If you find bugs, please contact the maintainer of this repository, e.g. via the 'Issues' function.

******************************
WASA-SED is a large and complex hydrological and sediment transport model. It is not a polished point-and-click GUI tool. Extensive knowledge of its design, purpose and limitations is required in order to apply it properly. See [```license.txt```](https://github.com/TillF/WASA-SED/blob/master/license.txt) for terms of use.
******************************

Contents
--------
* docs/: Directory containg the documentation
  * CONTRIBUTING.md: Guidelines for modifying the code
  * index.md: technical documentation (also available via https://tillf.github.io/WASA-SED/)
  * variables.ods: List of variables used in the source code (incomplete)
  * tutorial_wasa_input.zip: Example input to run the model

* src/: Contains the model's Fortran 90 source code files structured in sub-directories
  * update_revision_no.\*: Script to add the revision number into the model executable when compiling the source code
  
* license.txt: terms of use

* Makefile: GNU Makefile with rules for compiling the source code and installing the model

* README.md: THIS readme in markdown language

* .gitignore: List of files / directories git shall ignore (only needed for git usage)


Installation
--------
To run the model, use the WASA.exe executable. Windows user may use the file provided in the [release-section](https://github.com/TillF/WASA-SED/releases).

Alternatively, compile the source code yourself (recommended). A Makefile (adjusted to _GNU make_, on Windows, use [MinGW](http://mingw.org/)) is provided. For more information on how to compile and install the model, type `make help` in the console within the WASA-SED main directory. 
To specify your desired compiler or adjust the compiler flags, you have to alter the Makefile. By default, it runs with the _GNU fortran_ compiler _gfortran_. The compiler flags are adapted to version 5.x. For more information, see comments in the Makefile.

Reference
---------

E.N. Mueller, A. Güntner, T. Francke, G. Mamede (2010): Modelling sediment export, retention and reservoir sedimentation in drylands with the WASA-SED model Geosci. Model Dev., 275-291, 3(1) , url: http://www.geosci-model-dev.net/3/275/2010/, doi:10.5194/gmd-3-275-2010

Auxiliary Tools
---------
[SoilDataPrep](https://github.com/TillF/SoilDataPrep) : soildata-preprocessing for WASA-SED

[lumpR](https://github.com/tpilz/lumpR) : Geodata-preprocessing for WASA-SED
