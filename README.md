!!!The executable (WASA.exe) are updated only infrequently, if at all.
To use the latest version, you need to compile from source code. !!!


WASA-SED Model
--------------
WASA-SED is a numerical model for simulation hydrological and sediment fluxes from meso-scale catchments (Mueller et al., 2010).
This repository contains documentation, sample data set and source code.
All rights under the stated license (see license.txt) of the WASA-SED code are with the SESAM-Project and successors, c/o Till Francke, Universität Potsdam, Karl-Liebknecht-Str. 24-25, 14473 Potsdam, Germany, please contact: francke@uni-potsdam.de).

******************************
WASA-SED is a large and complex hydrological and sediment transport model. It is not a polished point-and-click GUI tool. Extensive knowledge of its design, purpose, and limitations is required in order to apply it properly. See license.txt for terms of use.
******************************

Contents
--------
* doc/: Directory containg the documentation
  * coding_guidelines.txt: Guidelines for modifying the code
  * Wasa_Documentation.doc: Model technical documentation
  * variables.ods: List of variables used (incomplete)
  * tutorial_wasa_input.zip: Example input to run the model

* src/: Contains the model's Fortran 90 source code structures in sub-directories
  * update_revision_no.\*: Script to add the revision number into the model executable when compiling the source code
  
* license.txt: terms of use

* Makefile: GNU Makefile with rules for compiling the source code and installing the model

* README.md: THIS readme in markdown language

* WASA.exe: Model executable for Windows (**only infrequently updated!**)

* .gitignore: List of files / directories git shall ignore (only needed for model development)


Installation
--------
To run the model, use the WASA.exe executable (Windows; **only infrequently updated!**) or compile the source code on your own (suggested). To do so, a Makefile (adjusted to _GNU make_) is provided.

For more information on how to compile and install the model, in a command line within the WASA-SED main directory type `make help`.

To specify your desired compiler or adjust the compiler flags, you have to alter the Makefile. By default, it runs with the _GNU fortran_ compiler. The compiler flags are adapted to version 5.x. For more information, see comments in the Makefile.

Reference
---------

E. Mueller, A. Güntner, T. Francke (2010): Modelling sediment export, retention and reservoir sedimentation in drylands with the WASA-SED model Geosci. Model Dev., 275-291, 3(1) , url: http://www.geosci-model-dev.net/3/275/2010/, doi:10.5194/gmd-3-275-2010

