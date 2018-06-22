# WASA-SED

# User Manual

Eva Nora Müller, Till Francke, George Mamede and Andreas Güntner

<br>
<br>

**30.1.2018<br>
WASA-SED rev. 255**

<br>
<br>

Developed within the SESAM-Project:<br>
Sediment Export of Semi-Arid Catchment: Monitoring and Modelling 2005-2008<br>
SESAM II 2010-2014,<br>
2015-2018

<br>
<br>

Institute of Earth and Environmental Science, University of Potsdam, Potsdam,<br>
Deutsches Geoforschungszentrum Potsdam,<br>
Germany


<br>
<br>

**Recent change: snow module (“WASA-SNOW” by Erwin Rottler)**<br>
The respective documentation has not been included yet. It can be found in Rottler’s thesis (on request, to be published on https://publishup.uni-potsdam.de/opus4-ubp/home)

Updates of WASA-SED Manual (this file):
https://github.com/TillF/WASA-SED

Source-code (Fortran90) git-repository:
https://github.com/TillF/WASA-SED

Theoretical description of the WASA-SED model:

Mueller EN, Güntner A, Francke T, Mamede G (2010) Modelling sediment export, retention and reservoir sedimentation in drylands with the WASA-SED model. Geosci Model Dev 3:275–291. doi:10.5194/gmd-3-275-2010. http://www.geosci-model-dev.net/3/275/2010/

Information about the SESAM-Project:
http://uni-potsdam.de/sesam/

Contact:<br>
Till Francke<br>
Institute of Earth and Environmental Sciences<br>
University of Potsdam<br>
Karl-Liebknecht-Str. 24-25<br>
14476 Potsdam, Germany<br>
Email: till.francke@uni-potsdam.de
 
Copyrights and ownership:<br>
For the original WASA code: Andreas Güntner, Geoforschungszentrum Potsdam, Telegrafenberg, 14473 Potsdam, Germany.<br>
For all extensions of WASA-SED: SESAM-project team, University of Potsdam, 14476 Potsdam, Germany.<br>
\* Creative Commons Attribution 4.0 International License

Disclaimer:<br>
The WASA-SED program is large and complex and extensive knowledge of its design, purpose, and limitations is required in order to apply it properly. The WASA-SED and its source code is freely available under a CC4 licence (“use as you wish, don’t blame us, give credit,…”). See license.txt.
<br>
<br>
<br>

## Table of Content
- [Introduction](#heading)<br>
- [Program folders and structure](#heading)<br>
  * [Hillslope module](#sub-heading)<br>
  * [River module](#sub-heading)<br>
  * [Reservoir module](#sub-heading)<br>
  * [Input Data](#sub-heading)<br>
  * [General parameter and control files](#sub-heading)<br>
  * [Input files for the hillslope module](#sub-heading)<br>
  * [Input files for the river module](#sub-heading)<br>
  * [Input files for the reservoir module](#sub-heading)<br>
  * [Input of climate data](#sub-heading)<br>
- [Output Data](#heading)<br>
  * [Output of the hillslope module](#sub-heading)<br>
  * [Output of the river module](#sub-heading)<br>
  * [Output of the reservoir module](#sub-heading)<br>
- [Relevant Literature for the WASA-SED Model](#heading)<br>
- [Reference](#heading)<br>

## Introduction
The WASA-SED model simulates the runoff and erosion processes at the hillslope scale, the transport processes of suspended and bedload fluxes at the river scale and the retention and remobilisation processes of sediments in large reservoirs. The modelling tool enables the evaluation of management options both for sustainable land-use change scenarios to reduce erosion in the headwater catchments as well as adequate reservoir management options to lessen sedimentation in large reservoirs and reservoir networks. The model concept, its spatial discretisation and the numerical components of the hillslope, river and reservoir processes are summarised and current model applications are reviewed in Mueller et al. (2008). The hydrological routines of the model are based on the WASA model (Model for Water Availability in Semi-Arid environments), which was developed by Güntner (2002) and Güntner and Bronstert (2002, 2003) to enable the quantification of water availability in semi-arid regions. The WASA-SED model was developed within the joint Spanish-Brazilian-German research project SESAM (Sediment Export from Semi-Arid Catchments: Measurement and Modelling). The existing WASA model code has been extended to include sediment-transport routines for the three new conceptual levels of the WASA-SED model: the hillslope scale, river scale and the reservoir scale for the calculation of sedimentation. This documentation gives a short outline of the structure, computational routines and folder system of the WASA-SED code in Chapter 2, followed by a description of the input files for model parameterisation in Chapter 3 and output files for the hillslope, river and reservoir modules in Chapter 4.

## Program folders and structure
The WASA-SED model is programmed in Fortran 90 and was tested with Compaq Visual Fortran (6.6.a) and gfortran 4.x compilers on Windows gfortran-4.3 under Linux. The WASA-SED code is organised as presented in Table 1. The main program is named wasa.f90 and calls the key subroutines as summarised in Table 2.

Table 1 Folder structure of WASA code

|Folder Name | Content  |
|---|---| 
\General\ | Main program and utility routines of the WASA model  
\Hillslope\ | Hillslope routines (Overland, sub-surface flow, evapotranspiration, sediment production etc.)  
\River\	| River routines (Routing of water and sediment in the river network)
\Reservoir\ | Reservoir water and sediment modelling routines (Water balance, bed elevation change, management options)
\Input\	| Input data, contains parameter file: do.dat, 
          \Input\Hillslope, 
          \Input\River, 
          \Input\Reservoir, 
          \Input\Time_series
\Output\ | Output files of model scenarios 

A complete example input-data set is part of the Subversion repository (see top of document). Model parameterisations are available e.g. for meso-scale catchments in dryland areas of Spain and Brazil (Bengue Catchment in Spain: Mamede 2008, Ribera Salada Catchment in Spain: Mueller et al. 2008, Mueller et al. submitted to CATENA, Isábena Catchment in Spain: Francke 2009).

The original WASA code version (Güntner 2002, Güntner and Bronstert 2004) was extended within the SESAM-Project to include sediment-transport processes at the hillslope scale using various USLE-derivative approaches, a spatially distributed, semi-process-based modelling approach for the modelling of water and sediment transport through the river network and a reservoir module that computes the transport of water and sediment as well as sedimentation processes in reservoirs.



