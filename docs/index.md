# WASA-SED

# User Manual
main contributors:
Andreas Güntner (2000-2002), Eva Nora Müller (2006-2009), George Mamede(2006-2009), Till Francke (2006-...), Pedro Medeiros (2007-2008), Tobias Pilz (2016), Erwin Rottler (2017)

**30.6.2021<br>
WASA-SED rev_269**

Developed within the SESAM-Project:<br>
Sediment Export of Semi-Arid Catchment: Monitoring and Modelling 2005-2008<br>
SESAM II 2010-2014,<br>
2015-2022

Institute of Environmental Sciences and Geography, University of Potsdam, Potsdam,<br>
Deutsches Geoforschungszentrum Potsdam,<br>
Germany

Updates of WASA-SED Manual (this file):
https://github.com/TillF/WASA-SED

Source-code (Fortran90) git-repository:
https://github.com/TillF/WASA-SED

Theoretical description of the WASA-SED model:

Mueller, E.N., Güntner, A., Francke, T., Mamede, G. (2010): Modelling sediment export, retention and reservoir sedimentation in drylands with the WASA-SED model. Geosci Model Dev 3:275–291. doi:10.5194/gmd-3-275-2010. http://www.geosci-model-dev.net/3/275/2010/.

Information about the SESAM-Project:
http://uni-potsdam.de/sesam/

Contact:<br>
Till Francke<br>
Institute of Environmental Sciences and Geography<br>
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

<a name="toc"></a>
## Table of content
- [1 Introduction](#1-introduction)<br>
- [2 Program folders and structure](#2-program-folders-and-structure)<br>
  -	[2.1 Hillslope module](#2-1-hillslope-module)<br>
  - [2.2 River module](#2-2-river-module)<br>
  - [2.3 Reservoir module](#2-3-reservoir-module)<br>
- [3 Input data](#3-input-data)<br>
  - [3.1 General parameter and control files](#3-1-general-parameter-and-control-files)<br>
  - [3.2 Input files for the hillslope module](#3-2-input-files-for-the-hillslope-module)<br>
  - [3.3 Input files for the river module](#3-3-input-files-for-the-river-module)<br>
  - [3.4 Input files for the reservoir module](#3-4-input-files-for-the-reservoir-module)<br>
  - [3.5 Input of climate data](#3-5-input-of-climate-data)<br>
  - [3.6 Input/Output of state variables](#3-6-state-variables)<br>
- [4 Output data](#4-output-data)<br>
  - [4.1 Output of the hillslope module](#4-1-output-of-the-hillslope-module)<br>
  - [4.2 Output of the river module](#4-2-output-of-the-river-module)<br>
  - [4.3 Output of the reservoir module](#4-3-output-of-the-reservoir-module)<br>
- [5 References](#5-references)<br>
- [6 Further relevant literature and tools for the WASA-SED model](#6-further-relevant-literature-and-tools-for-the-wasa-sed-model)<br>


## List of figures
[Figure 1](#figure-1): WASA parameter specification file ```do.dat```.<br>
[Figure 2](#figure-2): Trapezoidal channel dimension with floodplains.

## List of tables
[Table 1](#table-1): Folder structure of WASA code.<br>
[Table 2](#table-2): Main subroutines of the main program (```wasa.f90```).<br>
[Table 3](#table-3): Main subroutines of ```hymo_all.f90``` (hydrological subroutines).<br>
[Table 4](#table-4): Summary of model input parameters.<br>
[Table 5](#table-5): Keywords and descriptions optional outfiles related to the snow routine.<br>
[Table 6](#table-6): Input data files for the hillslope component.<br>
[Table 7](#table-7): Groundwater option Gw_flag=0.<br>
[Table 8](#table-8): Groundwater option Gw_flag=1.<br>
[Table 9](#table-9): Groundwater option Gw_flag=99.<br>
[Table 10](#table-10): Input data files for the river component.<br>
[Table 11](#table-11): Input data files for the reservoir component.<br>
[Table 12](#table-12): Input/Output of state variables.<br>
[Table 13](#table-13): Output files of the hillslope module.<br>
[Table 14](#table-14): Output files of the river module.<br>
[Table 15](#table-15): Output files of the reservoir module.

## 1 Introduction

The WASA-SED model simulates the runoff and erosion processes at the hillslope scale, the transport processes of suspended and bedload fluxes at the river scale and the retention and remobilisation processes of sediments in large reservoirs. The modelling tool enables the evaluation of management options both for sustainable land-use change scenarios to reduce erosion in the headwater catchments as well as adequate reservoir management options to lessen sedimentation in large reservoirs and reservoir networks. The model concept, its spatial discretisation and the numerical components of the hillslope, river and reservoir processes are summarised and current model applications are reviewed in [Mueller et al. (2010)](#mueller-et-al-2010). The hydrological routines of the model are based on the WASA model (Model for Water Availability in Semi-Arid environments), which was developed by [Güntner (2002)](#guentner-2002) and [Güntner and Bronstert (2002,](#guentner-bronstert-2002) [2003)](#guentner-bronstert-2003) to enable the quantification of water availability in semi-arid regions. The WASA-SED model was developed within the joint Spanish-Brazilian-German research project SESAM (Sediment Export from Semi-Arid Catchments: Measurement and Modelling). The existing WASA model code has been extended to include sediment-transport routines for the three new conceptual levels of the WASA-SED model: the hillslope scale, river scale and the reservoir scale for the calculation of sedimentation. This documentation gives a short outline of the structure, computational routines and folder system of the WASA-SED code in [Chapter 2](#program-folders-and-structure), followed by a description of the input files for model parameterisation in [Chapter 3](#input-data) and output files for the hillslope, river and reservoir modules in [Chapter 4](#output-data).

\[[Table of contents](#toc)]
## 2 Program folders and structure

The WASA-SED model is programmed in Fortran 90 and was tested with Compaq Visual Fortran (6.6.a) and gfortran 4.x compilers on Windows gfortran-4.3 under Linux. The WASA-SED code is organised as presented in [Table 1](#table-1). The main program is named ```wasa.f90``` and calls the key subroutines as summarised in [Table 2](#table-2).

<a name="table-1"></a>
**Table 1:** Folder structure of WASA code.

Folder Name | Content
---|---
```/General```| Main program and utility routines of the WASA model  
```/Hillslope``` | Hillslope routines (Overland, sub-surface flow, evapotranspiration, sediment production etc.)  
```/River``` | River routines (Routing of water and sediment in the river network)
```/Reservoir``` | Reservoir water and sediment modelling routines (Water balance, bed elevation change, management options)
```/Input``` | Input data, contains parameter file: ```do.dat, /Input/Hillslope, /Input/River, /Input/Reservoir, /Input/Time_series```
```/Output``` | Output files of model scenarios 

A complete example input-data set is part of the github repository (see top of document). Model parameterisations are available e.g. for meso-scale catchments in dryland areas of Spain and Brazil (Bengue Catchment in Spain: [Mamede 2008](#mamede-2008), Ribera Salada Catchment in Spain: [Mueller et al. 2008](#mueller-et-al-2008), [Mueller et al. 2009](#mueller-et-al-2009), Isábena Catchment in Spain: [Francke 2009](#francke-2009)) and others.

The original WASA code version ([Güntner 2002](#guentner-2002), [Güntner and Bronstert 2004](#guentner-bronstert-2004)) was extended within the SESAM-Project to include sediment-transport processes at the hillslope scale using various USLE-derivative approaches, a spatially distributed, semi-process-based modelling approach for the modelling of water and sediment transport through the river network and a reservoir module that computes the transport of water and sediment as well as sedimentation processes in reservoirs. Other extensions (e.g. snow) were added later.

<a name="table-2"></a>
**Table 2:** Main subroutines of the main program (```wasa.f90```).

Routine	| Content
---|---
```readgen.f90```	| Reads parameter file: do.dat
```calcyear.f90```	| Counter for time
```calcday.f90```	| Number of days per month
```climo.f90```	| Climate calculations
```hymo_all.f90```	| Hydrological 
```routing.f90``` |	original river routing routine (unit hydrograph)
```routing_new.f90``` | river routine with Muskingum, calls reservoir routines

The following sections give some information on the computational background and the functional structure of the corresponding routines for each of the three conceptual levels: hillslope, river and reservoir.

<a name="2-1-hillslope-module"></a>
### 2.1 Hillslope module

The hillslope module comprises the modelling of the hydrological and sediment-transport processes. The hydrological modelling accounts for interception, evaporation, infiltration, surface and subsurface runoff, transpiration and ground water recharge. Details are given in [Güntner (2002](#guentner-2002), Chapter 4). The main hydrological calculations are carried out in ```hymo_all.f90``` (for daily or hourly time steps). The subroutines that are called within ```hymo_all.f90``` are summarised in [Table 3](#table-3). The temporal sequence of hydrological process modelling is summarised in [Güntner (2002](#guentner-2002), p. 36-37).

<a name="table-3"></a>
**Table 3:** Main subroutines of ```hymo_all.f90``` (hydrological subroutines).

Routine	| Content	| Subroutine	| Content of Subroutine
---|---|---|---
```readhymo.f90``` |	reads input data	| ```soildistr.f90``` | distribution functions for soil parameters
```lake.f90```	| water balance for small reservoirs	| - |	-
```soilwat.f90```	| soil water components (infiltration, percolation, runoff generation)	|```etp_max.f90, etp_soil.f90, sedi_yield.f90```	| maximal evapotranspiration, daily evapotranspiration, hillslope erosion

Sediment generation on the hillslopes in the form of soil erosion by water is modelled with four erosion equations (USLE, Onstad-Foster, MUSLE and MUST after [Williams, 1995](#williams-1995)). The following subroutines contain the main calculations for the water and sediment production for the hillslopes:

```hymo_all.f90```: general loop in daily or hourly timesteps
1.	Initialisation and reading of hillslope input data, preparation of output files
2.	For each timestep
    -	nested loop that contains the calculation for hydrological and sediment budgets for each sub-basin, landscape unit, terrain component, soil-vegetation component (```soilwat.f90```)
    -	aggregated runoff refers to the whole available area including reservoir areas; the corresponding reservoir areal fraction is substracted afterwards if reservoirs / small reservoirs are considered
    -	call module for small (diffuse, unlocated) reservoirs (```lake.f90```)
    - write daily results into output files

```soilwat.f90```: computes water and sediment balance for a terrain component in a given landscape unit of a given sub-basin
1.	For each soil vegetation component (SVCs) in the terrain component: calculate hydrological variables (infiltration, surface runoff, subsurface runoff, evapotranspiration, etc.)
2.	Re-distribute surface runoff, subsurface runoff among SVCs, allow re-infiltration and compute overall water balance for 
terrain component
3.	Estimate runoff height and peak-runoff using Manning Equation
4.	Call module for hillslope erosion (```sedi_yield.f90```)

```sedi_yield.f90```: computes sediment production and sediment export for a TC
1.	No sediment export without surface runoff out of terrain component
2.	Otherwise, compute gross erosion for current terrain component using one of four erosion equations
3.	Distribute gross erosion among particle classes according to constitution of uppermost horizons in the current terrain component
4.	Completely mix the newly generated sediment with any sediment coming from upslope terrain components
5.	If necessary, limit sediment export by transport capacity.

<a name="2-2-river-module"></a>
### 2.2 River module


The river routing of the original WASA model ([Güntner 2002](#guentner-2002)) bases on daily linear response functions ([Bronstert et al. 1999](#bronstert-et-al-1999)) similar to a triangular unit hydrograph. Its implementation does not support output in hourly resolution (only daily is produced). It was extended to include a spatially semi-distributed, semi-process-based modelling approach for the modelling of water and sediment transport through the river network. The implemented water modelling approach is similar to the routing routines from the SWAT model (Soil Water Assessment Tool, [Neitsch et al. 2002](#neitsch-et-al-2002)) model and the SWIM model (Soil Water Integrated Modelling, [Krysanova et al. 2000](#krysanova-et-al-2000)). The new water routing is based on the Muskingum kinematic wave approximation. Suspended sediment transport and bedload is modelled using the transport capacity concept. The river module can be run with variable time steps. Transmission losses through riverbed infiltration and evaporation are accounted for. The main routing calculations as well as the initialisation and reading of the river input files are carried out in ```routing.f90```. The following sub-routines are called from ```routing.f90```:

```muskingum.f90```: contains the flow calculation using the Muskingum method
1. Calculation of water volume in reach
2. Calculation of cross-section area of current flow, flow depth, wetted perimeter and hydraulic radius
3. Calculation of flow in reach with Manning Equation
4. Calculation of Muskingum coefficients
5. Calculation of discharge out of the reach and water storage in reach at end of time step

```routing_coefficients.f90```: contains the calculations for travel time and flow depth
1. Calculation of initial water storage for each river stretch
2. Calculation of channel dimensions
3. Calculation of flow and travel time at 100 % bankfull depth, 120 % bankfull depth and 10 % bankfull depth

```route_sediments.f90```: contains the calculation for suspended sediment-transport routing
1. Calculation of water volume and sediment mass in reach
2. Calculation of peak flow and peak velocity
3. Calculation of maximum sediment carrying capacity concentration
4. Comparison with current concentration
5. Calculation of net deposition and degradation
6. Calculation of sediment mass out of the reach and sediment storage in reach at end of time step

```bedload.f90```: contains the calculation for bedload transport using 5 different formulas
1.	Calculation of current width of the river
2.	Bedload formulas after Meyer-Peter and Müller (1948), Schoklitsch (1950), Smart and Jaeggi (1983), Bagnold (1956) and Rickenmann (2001), see [Mueller et al. 2008](#mueller-et-al-2008) relevant references.

<a name="2-3-reservoir-module"></a>
### 2.3 Reservoir module

The reservoir sedimentation routine was included into the WASA-SED model by [Mamede (2008)](#mamede-2008) to enable the calculation of non-uniform sediment transport along the longitudinal profile of a reservoir, of the reservoir bed changes caused by deposition/erosion processes and of reservoir management options. 

In order to perform the simulation of sediment transport in reservoirs, four important processes have to be considered: (1) reservoir water balance, (2) hydraulic calculations in the reservoir, (3) sediment transport along the longitudinal profile of the reservoir and (4) reservoir bed elevation changes. For the calculation of sediment transport in the reservoir, four different equations for the calculation of total sediment load were selected from recent literature. The reservoir bed elevation changes are calculated through the sediment balance at each cross section, taking into account three conceptual sediment layers above the original bed material. The reservoir sedimentation module is composed by the following subroutines:

```reservoir.f90```: contains the reservoir water balance
1. Calculation of reservoir level
2. Calculation of controlled and uncontrolled outflow

```semres.f90```: contains the sediment balance for each cross section
1. Calculation of sediment deposition for each cross section
2. Calculation of sediment entrainment for each cross section
3. Calculation of sediment compaction for each cross section
4. Calculation of total load transport
5. Calculation of storage capacity reduction because of the accumulated sediment
6. Calculation of effluent grain size distribution

```sedbal.f90```: contains a simplified sediment balance of the reservoir whenever its geometry is not provided (no cross-section)
1. Calculation of sediment deposition in the reservoir
2. Calculation of trapping efficiency
3. Calculation of storage capacity reduction because of the accumulated sediment
4. Calculation of effluent grain size distribution

```hydraul_res.f90```: contains the hydraulic calculations for each cross section in the reservoir
1. Calculation of mean velocity of current flow for each cross section
2. Calculation of hydraulic radius for each cross section
3. Calculation of slope of energy-grade line for each cross section
4. Calculation of top width for each cross section
5. Calculation of water depth for each cross section

```eq1_wu.f90```, ```eq2_ashida.f90```, ```eq3_tsinghua.f90``` and ```eq4_ackers.f90```: each sub-routine contains a specific sediment transport function, mentioned in Table 4
1. Calculation of fractional bed load transport for each cross section
2. Calculation of fractional suspended load transport for each cross section

```change_sec.f90```: contains the calculation of reservoir geometry changes
1. Calculation of bed elevation changes for each cross section.
```reservoir_routing.f90```: contains the calculation of level-pool reservoir routing 
1. Calculation of water routing for reservoirs with uncontrolled overflow spillways
```vert_dist.f90```: contains the calculation of vertical distribution of sediment concentration in the reservoir 
1. Calculation of sediment concentration at the reservoir outlets
```lake.f90```: contains the water balance for networks of small reservoirs 
1. Calculation of water level in the reservoirs
2. Calculation of controlled and uncontrolled outflow out the small reservoirs
```sedbal_lake.f90```: contains a simplified sediment balance for networks of small reservoirs
1. Calculation of sediment deposition in small reservoir
2. Calculation of trapping efficiency in small reservoirs
3. Calculation of storage capacity reduction because of the accumulated sediment in small reservoirs
4. Calculation of effluent grain size distribution in small reservoirs
```lake_routing.f90```: contains the calculation of level-pool routing for networks of small reservoirs
1. Calculation of water routing for small reservoirs

\[[Table of contents](#toc)]
## 3 Input data

The model runs as a Fortran console application for catchment from a few km² up to several 100,000 km²) on daily or hourly time steps. Climatic drivers are daily/hourly time series for precipitation, humidity, short-wave radiation and temperature. For model parameterisation, regional digital maps on soil associations, land-use and vegetation cover, a digital elevation model with a cell size of 100 metres (or smaller) and, optional, bathymetric surveys of the reservoirs are required. Soil data can be obtained and preprocessed with SoilDataPrep-tool(#SoilDataPrep). 
The soil, vegetation and terrain maps are processed with the lumpR tool (see below) to derive the spatial discretisation into soil-vegetation units, terrain components and landscape units. [Table 4](#table-4) summarises the input parameters for the climatic drivers and the hillslope, river and reservoir modules. The vegetation parameters may be derived with the comprehensive study of, for example, [Breuer et al. (2003)](#breuer-et-al-2003), the soil and erosion parameters with the data compilation of [FAO (1993,](#fao-1993) [2001)](#fao-2001), [Morgan (1995)](#morgan-1995), [Maidment (1993)](#maidment-1993) and [Antronico et al. (2005)](#antronico-et-al-2005).

For a semi-automated discretisation of the model domain into landscape units and terrain components, the LUMP-algorithm (Landscape Unit Mapping Program) is available ([Francke et al. 2008](#francke-et-al-2008)). This algorithm delineates areas with similar hillslope characteristics by retrieving homogeneous catenas with regard to e.g. hillslope shape, flow length and slope (provided by a digital elevation model), and additional properties such as for soil and land-use and optionally for specific model parameters such as leaf area index, albedo or soil aggregate stability. LUMP can be linked with the WASA-SED parameterisation procedure through a data base management tool, which allows to process and store digital soil, vegetation and topographical data in a coherent way and facilitates the generation of the required input files for the model. LUMP and further WASA-SED pre-processing tools have been transferred to the package lumpR ([Pilz et al. 2017](#pilz-et-al-2017)) for the free software environment for statistical computing and graphics R which is available from https://github.com/tpilz/lumpR.

The input files for general purpose, the hillslope, river and reservoir routines are explained below with details on parameter type, units, data structure including examples of the parameterisation files.

<a name="table-4"></a>
**Table 4:** Summary of model input parameters.

Type |	Model input parameter
---|---
Climate |	Daily or hourly time series on rainfall \[mm/day, mm/h], Daily time series for average short-wave radiation \[W/m2], Daily time series for humidity \[%], Daily time series for temperature \[°C]
Vegetation	| Stomata resistance \[s/m], Minimum suction \[hPa], Maximum suction \[hPa], Height \[m], Root depth \[m], LAI \[-], Albedo \[-], Manning's n of hillslope \[-], USLE C \[-]
Soil	| No. of horizons\*, Residual water content \[Vol. %], Water content at permanent wilting point \[Vol. %], Usable field capacity \[Vol. %], Saturated water content \[Vol. %], Saturated hydraulic conductivity (mm/h), Thickness \[mm], Suction at wetting front \[mm], Pore size index \[-], Bubble pressure \[cm], USLE K \[-]\*\*, Particle size distribution\*\*
Terrain and river |	Hydraulic conductivity of bedrock \[mm/d], Mean maximum depth of soil zone \[mm], Depth of river bed below terrain component \[mm], Initial depth of groundwater below surface \[mm], Storage coefficient for groundwater outflow \[day], Bankful depth of river \[m], Bankful width of river \[m], Run to rise ratio of river banks \[-], Bottom width of floodplain \[m], Run to rise ratio of floodplain side slopes \[-], River length \[km], River slope \[m/m], D<sub>50</sub> (median sediment particle size) of riverbed \[m], Manning’s n for riverbed and floodplains \[-]
Reservoir	| Longitudinal profile of reservoir \[m], Cross-section profiles of reservoir \[m], Stage-volume curves, Initial water storage and storage capacity volumes \[m<sup>3</sup>], Initial area of the reservoir \[ha], Maximal outflow through the bottom outlets \[m<sup>3</sup>/s], Manning’s roughness for reservoir bed, Depth of active layer \[m]

\* for each soil horizon, all following parameters in the column are required <br>
\** of topmost horizon

<a name="3-1-general-parameter-and-control-files"></a>
### 3.1 General parameter and control files

Four parameter files control the data input and output and some internal settings:

```do.dat``` <br>
\[can be generated with the lumpR package, manual adjustments required]

The file ```do.dat``` contains the main parameter specifications for the WASA-SED model. By default, it is located in the folder WASA\Input, but any other location may be specified as a command line argument (e.g. ```wasa.exe my_own/do.dat```).  [Figure 1](#figure-1) displays an example. Line order is important; do not add nor delete any lines. The first line of the ```do.dat``` contains the title. Line 2 and 3 specify the path for the location of WASA input and output folder. Relative paths are supported. The backslash “\” only works on Windows-platforms. The slash “/” is accepted on Windows and Unix/Linux systems. Make sure that both specified paths end with slash or backslash, respectively. Line 4 and 5 contain the start and the end year of the simulation, respectively. Line 6 and 7 contain the start and the end calendar month of the simulation, respectively. Optionally, the day of month for begin and end can be specified separated by a blank space. Line 10 contains the number of sub-basins. The number in line 9 is given by the sum of the number of terrain components in each landscape-unit of each sub-basin (e.g. if the system has only two sub-basins, sub-basin A has 1 landscape unit with 3 terrain components, sub-basin B has 2 landscape units with 1 terrain component each, then the number of combinations is 5). Line 14 specifies if the reservoir module is switched on (.t.) or is switched off (.f.). The same issue for the calculations of networks of small reservoirs in line 15. Lines 16 – 19 allow customizing the way water and sediment is (re-)distributed within and among the TCs. Line 21 allows the setting of the simulation timestep (daily / hourly). This may become obsolete in future versions by setting the timestep directly in line 30. Line 24 allows specifying a correction factor for hydraulic conductivity to account for intra-daily rainfall intensities. Optionally, this factor can also be made a function of daily rainfall by specifying two more parameters (a and b) in the same line, so that kfkorr=kfkorr0\*(a\*1/daily_precip+b+1). In line 31 the erosion and sediment-transport routines may be switched on and off. Specify the number of grain size classes you want to model in line 32. Their limits must be specified in ```part_class.dat```, if more than one class is desired. Line 33 lets you choose the hillslope erosion model to be used in WASA. Currently, this parameter is disregarded, further options can be chosen in ```erosion.ctl```. Select the model for the river routing in line 34. Possible options are: (1) UHG routing (daily resolution only, (2) Muskingum and suspended sediment, (3) Muskingum and bedload transport. Choose the sediment model in the reservoir in line 35 among 4 sediment transport equations: (1) [Wu et al. (2000)](#wu-et-al-2000); (2) [Ashida and Michiue (1973)](#ashida-michiue-1973); (3) [Yang (1973,](#yang-1973) [1984)](#yang-1984); (4) [Ackers and White (1973)](#ackers-white-1973).

The optional lines 36 and 37 allow the saving/loading of state variables (i.e. groundwater, interception and soil storages) at the end/beginning of a model run (works only if ```svc.dat``` has been specified).

Line 36 may additionally contain a second logical variable (append_output), allowing the model to append to existing output files (not yet implemented for all outputs). The consistence with the existing files is not checked! Default is .FALSE.

Line 37 may additionally contain a second logical variable (save_states_yearly), determining if the model states are saved (and overwritten) at the end of each simulation year. Default is .TRUE.

Line 38 (dosnow) defines, if the optional snow routine, implemented by [Rottler (2017)](#rottler-2017), is included. If set to .TRUE., the snow routine is active and all related calculations are performed. If set to .FALSE., the snow routine remains inactive.

The line numbers in the following template are for reference only, they must be removed in the actual file!

<a name="figure-1"></a>
```
 1   Parameter specification for the WASA Model (SESAM-Project)
 2   ..\WASA\Input\Case_study\ 	Location of model platform
 3   ..\WASA\Output\		Specification of folder for simulation output
 4   1980  //tstart (start year of simulation) 
 5   1981  //tstop (end year of simulation)
 6    1 15    //mstart (start month of simulation [optional: <space>start_day])
 7    12 15    //mstop (end month of simulation [optional: <space>end_day])
 8   10    //no. of sub-basins
 9   49    //no. of combinations of sub-basins, landscape units, terrain components
10   321   //total no. of landscape units in study area
11   515   //total no. of terrain components in study area
12   77    //total no. of soil components in study area
13   34    //total no. of vegetation units in study area
14   .f.   //doreservoir: do reservoir calculations
15   .f.   //doacudes:includes dam calculations
16   .t.   //dolattc: do latflow between TCs
17   .f.   //doalllattc: rout latflow completely to next downslope TC
18   .t.   //dolatsc: do latflow within TCs (surface runoff)
19   .t.   //dolatscsub: do latflow within TCs (subsurface runoff)
20   .f.   //dotrans: do water transpositions betwen sub-basins
21   .f.   //dohour: do hourly version (ignored, use “dt”)
22   0     //scenario: choose scenario (0:less rain (ECHAM), 1:no trend, 2:more rain (Hadley))
23   0     //krig: type of precipitation interpolation (0….)0
24   15.0  //kfkorr:  hydraulic conductivity factor (for daily model version) (kfkorr)
25   0.30  //intcf: interception capacity per unit LAI (mm)
26   0     //dointc: type of interception routine (simple bucket:0, modified bucket:1)
27   .f.   //doscale: do scaling due to rainfall interpolation ?
28   .f.   //domuncell: for muni/ezg-nocell-version, use rainfall input derived from cells ? (change kf_calib.dat !)
29   1     //sensfactor: factor for sensitivity studies
30   24    //dt: time step in [hours]
31   .t.   //dosediment
32   1     //No. of grain size classes
33   1     // type of sediment transport model at the hillslope	
34   1     //type of water / sediment model in the river: (1) UHG routing, (2) Muskingum & ss transport, (3) Muskingum & bedload modelling
35   1     //type of sediment model in the reservoir: choose sediment transport …
36   .t. .f.  //load state of storages from files (if present) at start (optional); second flag: append output
37   .f.   //save state of storages to files after simulation period (optional)
38   .t.   //dosnow: activate snow routine
```

**Figure 1:** WASA parameter specification file ```do.dat```.

```maxdim.dat``` <br>
optional \[can be generated with the lumpR package\]<br>
The file ```maxdim.dat``` serves to optimise memory management and thus improves computational performance. It is, however, optional and if not encountered, default values are assumed.

```
contains maximum dimensions of spatial units
3	//maximum no. of landscape units in a sub-basins
3	//maximum no. of terrain components in a landscape unit
4	//maximum no. of soil vegetation components in a terrain component
6	//maximum no. of horizons in a soil
2	//maximum no. transpositions between sub-basins
0	//number of cross sections at the strategic reservoirs [optional]
0	//number of points (x,y) along the cross sections [optional]
```

Example: In the given example, no more than 3 LU may occur in one sub-basin (line 2). Analogously, no LU may contain more than 3 TCs (line 3) and no TC more than 4 SVCs (line 4). The number of horizons in a soil is limited to 6 (line 5). No more than 2 transpositions between sub-basins may exist. The last two lines are optional and valid only for computation of sedimentation patterns in strategic reservoirs (assumed 200, if missing).

```part_class.dat``` <br>
optional \[can be generated with the lumpR package\]<br>
The file ```part_class.dat``` is only necessary if sediment transport in multiple particle-size classes is to be modelled. If ```part_class.dat``` is missing, sediment transport will be modelled for a single particle-size class only. Otherwise, the file defines the number and the properties of the particle sizes that will be modelled. Please note that class numbering has to be continuous, starting with 1. The particle size classes must be ordered from fine to coarse.

```
# Particle size classes to be used in sediment modelling
class_number    upper_limit[mm]
1	            0.002
2	            0.063
3	            2.0
```

Class_number: continuous numbering of particle classes <br>
Upper_limit: upper limit of particle size for the respective class \[mm\]<br>

Example: The example file describes the 3 particle-size-classes clay, silt and sand (according to German classification) with clay particles up to 0.002 mm, silt (0.002 - 0.063 mm) and sand (0.063 - 2.0 mm).

```outfiles.dat``` <br>
optional<br>
The file allows specifying, which output files are desired. Disabling unnecessary output files saves computation time and disk space. The file contains two headerlines, each following line contains a keyword, which is the filename of a possible output file (case insensitive, without the extension ```.out```). If a keyword for a certain output file is not contained in ```outfiles.dat``` the respective file is not created, any existing file of that name is deleted. Information on the content of output files can be found in the respective sections. If ```outfiles.dat``` is not found, WASA-SED creates a default set of output files. 

```
# This file describes which output files are generated
# put any character before the files you don't want to be created
//hillslope
daily_actetranspiration
#daily_potetranspiration
#daily_qhorton
daily_qin_m3s
daily_qout_m3s
daily_rain
daily_runoff
daily_sediment_production
daily_subsurface_runoff
daily_theta
daily_total_overlandflow
daily_water_sub-basin
sediment_production
water_sub-basin
deep_gw_recharge
deep_gw_discharge
tc_theta
qhorton
actetranspiration
subsurface_runoff
total_overlandflow
gw_discharge
potetranspiration
gw_loss
gw_recharge

//river
routing_response
river_degradation
river_deposition
river_flow
river_flow_dailyaverage
river_flowdepth
river_sediment_concentration
river_sediment_total
river_sediment_total_dailyaverage
river_storage
river_sediment_storage
river_susp_sediment_storage
river_velocity
river_bedload
river_infiltration
tc_surfflow
tc_sedout
lu_sedout

//reservoir
res_watbal
res_vollost
res_cav
res_hydraul
res_bedchange
res_sedbal
res_longitudunal
res_sedcomposition
lake_inflow_r
lake_outflow_r
lake_retention_r
lake_storage_r
lake_sedinflow_r
lake_sedoutflow_r
lake_sedretention_r
lake_sedimentation_r
lake_volume_r
lake_watbal
lake_sedbal
lake_inflow
lake_outflow
lake_storage
lake_retention
lake_vollost
lake_maxstorcap
lake_volume
lake_sedinflow
lake_sedoutflow
lake_sizedistoutflow
```

Example: The output files ```daily_actetranspiration.out``` and ```daily_qhorton.out``` will be created. The creation of ```daily_potetranspiration.dat``` is omitted.

The output files of the snow routine ([Rottler, 2017](#rottler-2017)) are also defined in the control file ```outfiles.dat```. These files can be generated for all state variables and mass and energy fluxes related to the snow routine. A list of possible output files with corresponding keyword is given in [Table 5](#table-5). The output is created TC-wise. For each time step, values for each TC of each LU get exported.

<a name="table-5"></a>
**Table 5**: Keywords and descriptions optional outfiles related to the snow routine.

Keyword | Description
---|---
snowEnergyCont | Snow energy content of snow cover \[kJ/m²]
snowWaterEquiv | Snow water equivalent snow cover \[m]
snowAlbedo | Snow albedo \[-]
snowCover | Snow cover \[-]
snowTemp | Snow temperature \[°C]
surfTemp | Snow surface temperature \[°C]
liquFrac | Fraction of liquid water in snowpack \[-]
fluxPrec | Precipitation flux \[m/s]
fluxFlow | Melt water flux \[m/s]
fluxSubl | Sublimation flux \[m/s]
fluxNetS | Short-wave radiation balance \[W/m²]
fluxNetL | Long-wave radiation balance \[W/m²]
fluxSoil | Soil heat flux \[W/m²]
fluxSens | Sensible heat flux \[W/m²]
stoiPrec | Conversion factor mass and energy flux precipitation \[kJ/m<sup>3</sup>]
stoiSubl | Conversion factor mass and energy flux sublimation \[kJ/m<sup>3</sup>]
stoiFlow | Conversion factor mass and energy flux melt water \[kJ/m<sup>3</sup>]
rateAlbe | Change rate snow albedo \[1/s]
precipMod | Modified precipitation signal \[mm]
cloudFrac | Cloud cover fraction \[-]
radiMod | Radiation signal corrected for aspect and slope \[W/m²]
temperaMod | Height-modified temperature signal \[°C]
rel_elevation | Relative elevation of TC to mean sub-basin \[m]

<a name="3-2-input-files-for-the-hillslope-module"></a>
### 3.2 Input files for the hillslope module

The input files for the hillslope module are located in the folder ```Input/[case_study]/Hillslope``` and are summarised in [Table 6](#table-6).

<a name="table-6"></a>
**Table 6:** Input data files for the hillslope component.

Parameter File |	Content
---|---
```hymo.dat```	|	Specification of the sub-basins
```soter.dat```	| Specification of the landscape units
```terrain.dat```	| Specification of the terrain components
```svc.dat``` |	Specifications of soil vegetation components, erosion parameters
```soil_vegetation.dat```	| Specification of soil-vegetation components
```svc_in_tc.dat``` |	Specification which SVCs are contained in each TC
```soil.dat```	| Specification of the soil properties
```horizon_particles.dat```	| Particle size distributions of soil horizons
```vegetation.dat``` |	Specification of the vegetation properties
```rainy_season.dat``` (optional) | Specification of the start and end of rainy season
```k_seasons.dat``` (optional) | Seasonality of USLE factor K
```c_seasons.dat``` (optional) | Seasonality of USLE factor C
```p_seasons.dat``` (optional) | Seasonality of USLE factor P
```coarse_seasons.dat``` (optional) | Seasonality of coarse fraction
```n_seasons.dat``` (optional)	 | Seasonality of Manning’s n
```scaling_factor.dat``` (optional)	| Scaling/calibration factors for hydraulic conductivity
```calibration.dat``` (optional)	| Calibration factors for hydraulic conductivity of soils
```../erosion.ctl``` (optional)	| Options for the erosion module
```frac_direct_gw.dat``` (optional)	| partitioning of groundwater into river and alluvia
```beta_fac_lu.dat``` (optional) |  Correction factors for beta (USLE L-factor computation)
```sdr_lu.dat``` (optional)	| LU-wise specification of sediment delivery ratio
```calib_wind.dat``` (optional)	| Calibration of wind speed (sensitive parameter for evapotranspiration)
```snow_params.ctl```(optional) | Options for the snow module
```lu2.dat```(optional) | LU-properties for the snow module
```gw_storage.stat, intercept_storage.stat, soil_moisture.stat, snow_storage.stat``` (optional) | Initialisation state variables, see [section 3.6](#3-6-state-variables)


The spatial conceptualisation of the WASA model is explained in detail in [Güntner (2002), p. 33](#guentner-2002). The following spatial modelling units are used identified:

-	*Sub-basins:* ca. 50-1000 km<sup>3</sup>, topologically referenced, defined e.g. by the location of river gauging stations, or large reservoirs with a storage capacity of more than 50x10<sup>6</sup> m<sup>3</sup> and the confluence of major rivers;

-	*Landscape units (LUs):* based on the LU concept (e.g. SOil and TERrain digital database, [FAO, 1993](#fao-1993)), i.e. structure of the landscape according to geological, topographic and soil characteristics with similarity in major landform, general lithology, soil associations and toposequences. Characterized only by their fractions in each sub-basin, not georeferenced;

-	*Terrain components (TCs):* part of a landscape unit with similarity in slope gradients, position within toposequence (e.g. highlands, slopes and valley bottoms), soil association; Characterized only by their fractions and position along hillslope in each landscape unit, not georeferenced;

-	*Soil vegetation components (SVCs):* part of a terrain component with specific combination of soil type, and vegetation/land-cover class; characterized only by their fractions in each TC, not georeferenced; 

-	*Soil profile:* descriptions of characteristic sequence of soil horizons.

The model domain is divided into sub-basins; each sub-basin has an individual ID. This ID has to be a unique number; the employed numbering scheme does not have to be continuous. The following paragraphs explain each of the respective input files.

**1)** ```hymo.dat```<br>
\[can be generated with the lumpR package\]

```
# Specification of the sub-basins and their total number, type & areal fraction of LUs
Subbasin-ID [-], Area[km**2],  nbr[-],  LU-IDs[-], areal fraction of LU[-]
49	  10    4   19   87   90  135  0.357  0.147  0.214  0.282
50	  15    2   19   87   0.827  0.173
 1	  40    3   19  103   87  0.612  0.143  0.245
44	  37    3   19   87   18  0.646  0.097  0.257
10	   5    3   19   18   87  0.483  0.173  0.344
 4	  22    3  159   19   18  0.107  0.793  0.100
15	  13    2   87  138   0.740  0.260
39	  25    5   87   18  142   56   31  0.351  0.146  0.142  0.237  0.124
 3	  31    3   87   18  136  0.416  0.471  0.113
29	  20    4   56  122   31    7  0.091  0.652  0.131  0.126                                                
```
<br>
*Subbasin-ID*: ID of sub-basin<br>
*Area*:    Area of sub-basin in \[km²\] (including reservoir areas)<br>
*nbr*:     Number of landscape units (LUs) in sub-basin<br>
*LU-IDs*:      List of LU-ids which occur in this specific sub-basin<br>
*areal fraction*:	Fraction of each LU within the sub-basin \[-]

Example: In ```do.dat```, it was specified that 10 sub-basins are to be simulated with the WASA model. Accordingly, the file ```hymo.dat``` above contains the specification of 10 sub-basins, with the map IDs 49, 50,  1, … 29. The first sub-basin has a Map ID of 49 and an area of 10 km². Within this sub-basin, four different landscape units can be identified with the LU-IDs 19, 87, 90 and 135.  The first LU (ID 19) covers an area of 35.7 % (0.357) of the total area of the sub-basin, the second one 14.7 %, the third one 21.4 % and the last one 28.2 % (total 100 %).

Important: Any sub-basin that is not listed in the file ```routing.dat``` will be ignored.

When the area of a sub-basin is set to zero, the soilwater-calculations in this basin are skipped. This may be useful when introducing dummy-basins for routing reasons.

**2)** ```soter.dat``` <br>
\[can be generated with the lumpR package\]

```
# Specification of LUs										
LU-ID [-], No._of_TC[-], TC1[-], TC2[-], TC3[-], kfsu[mm/d], length[m], meandep[mm], maxdep[mm], riverbed[mm], gwflag[0/1], gw_dist[mm], frgw_delay[day]
1	3	7	49	   11	100	601	  -1	-1	1500	0	0	1
2	1	2	100	1963.7   -1	 -1	1500	 0	   0	1
…
```
<br>
*LU_id*:		    ID of landscape units<br>
*Nb.\_of_TC*:		Number of terrain components<br>
*TC1*:		    ID of a terrain component <br>
*TC2*:		    ID of another terrain component<br>
*TC3*:			ID of a third terrain component<br>
...			    more TC-IDs according to field 2<br>
*kfsu*:			Hydraulic conductivity of bedrock \[mm/d]<br>
*length*:			Mean  slope length in LU\[m]<br>
*meandep*:		Mean maximum depth of soil zone \[mm]<br>
*maxdep*:			Maximum depth of alluvial soil zone \[mm]<br>
*riverbed*:		Depth of river bed below terrain component \[mm]<br>
*gw_flag*:		Flag for LU \[0: no groundwater, 1: with groundwater] <br>
*gw_dist*:			Initial depth of groundwater below surface \[mm] (ignored, unless gw_flag=99)<br>
*frgw_delay*:		Storage coefficient for groundwater outflow \[day]<br>

Example: The landscape unit with ID 1 has 3 terrain components with the IDs 7, 49 and 11, a hydraulic conductivity of bedrock of 100 mm/d, a mean slope length of 601 m, etc. The LU with ID 2 has only 1 terrain component with the ID-Number 2 (i.e. consisting only of one rather homogenous hillslope section, TC2 and TC3 are set to zero), a hydraulic conductivity of bedrock of 100 mm/d, a mean slope length of 1963.7 m, etc. The TCs within a LU can be listed in any order, their position in the toposequence is read from terrain.dat.
LU-ids not lister in ```hymo.dat``` are ignored.

Remarks concerning groundwater: In WASA-SED, the representation of groundwater is yet quite simplistic. The general groundwater regime is essentially controlled by the flag for groundwater (```gw_flag```) in ```soter.dat```. Currently, the three options for setting ```gw_flag``` are 0, 1 or 99 (see [Table 7](#table-7), [8](#table-8), [9](#table-9), respectively). This determines if the water that leaves the soil column is “lost” for the model (```gwflag=0```) or enters a linear storage (```gwflag=1```). ```gw_flag=99``` is an experimental modelling option which usually should not be chosen.

Please note that the depth of the bedrock may be specified in two ways: In ```soil.dat```, the bedrock is assumed beneath the deepest soil horizon, if the bedrock flag is set (```bedrock=1```). Otherwise (```bedrock=0```), its depth is taken from the ```meandep``` (```maxdep``` for alluvial soils) specification given in ```soter.dat```.

Please note that in the case of alluvial soils, for each option in Tables 1 or 2 below, ```maxdep``` as defined in ```soter.dat``` instead of ```meandep``` is used. ```maxdep``` usually should be set larger than ```meandep```. The rationale behind this is that alluvial soils in the valley bottoms, in particular in crystalline bedrock environment, usually have larger soil depths (and thus higher water storage capacity) than average soils somewhere else (at the slopes) within the landscape unit. 

IMPORTANT:
In any of the above options, ```riverbed``` is to be defined in ```soter.dat```. In WASA, only soil horizons of the lowest terrain component which are located at depths above the riverbed are allowed to exfiltrate into the river by lateral flow. Soil horizons below ```riverbed``` cannot lose water to the river, but only due to evapotranspiration or percolation to groundwater/bedrock.

<a name="table-7"></a>
**Table 7:** Groundwater option ```gw_flag=0```.

Modelling options |	Groundwater regime, internal representation of processes | gw flag (```soter.dat```) | Bed-rock (```soil.dat```) / Add. parameters (```soter.dat```)
---|---|---|---
I) Groundwater below soil zone is ignored / not relevant for surface processes |	Soil water balance is modelled within the zone above bedrock. Water percolation out of the deepest soil horizon leaves the model domain	| 0	/ Select option 1.1 or 1.2 |
I.1) No bedrock is taken into account | Water percolation out of the deepest soil horizon leaves the model domain. | 0	| Select option 1.1.1 or 1.1.2
I.1.1) Depth of soil zone is the sum of all horizons given in ```soil.dat``` | See 1.1) |	0 |	0	/ meandep=-1, maxdep=-1, riverbed
I.1.2) Depth of soil zone is meandep (or the total profile depth given in ```soil.dat``` if this is larger than meandep) | See 1.1) |	0 | 0 / Meandep, Maxdep, riverbed
I.2) Bedrock is taken into account below deepest soil horizon | Hydraulic conductivity of bedrock may affect percolation rates out of the deepest soil horizon. A saturated zone (“groundwater”) above bedrock may consequently evolve due to saturation. Groundwater percolating into the bedrock leaves the model domain. | 0 | Select option 1.2.1 or 1.2.2
I.2.1) Bedrock is given in ```soil.dat```	| See 1.2) | 0 | 1 / Kfsu, riverbed
I.2.2) If bedrock is not given in ```soil.dat```, bedrock is assumed to be in the depth defined by meandep (or below deepest horizon given in ```soil.dat``` if its depth is greater than meandep) | See 1.2) |	0 | 0 / Meandep, maxdep, kfsu, riverbed

<a name="table-8"></a>
**Table 8:** Groundwater option ```gw_flag=1```.

Modelling options | Groundwater regime, internal representation of processes | Gw flag (```soter.dat```) | Bed-rock (```soil.dat```) / Add. parameters (```soter.dat```)
---|---|---|---
II) Groundwater below soil zone is taken into account | Soil water balance is modelled within the zone above bedrock. Water percolation out of the deepest soil horizon enters linear groundwater storage | 1 | Select option 2.1 or 2.2
II.1) No bedrock is taken into account | Water percolation out of the deepest soil horizon goes directly into the groundwater storage. | 1 | 0 / meandep=-1, maxdep=-1, frgw_delay, riverbed
II.2) Bedrock is taken into account below deepest soil horizon | The hydraulic conductivity of the bedrock may affect percolation rates out of the deepest soil horizon. A saturated zone (“groundwater”) above bedrock may consequently evolve due to saturation. Percolation into the bedrock is added to linear groundwater storage. | 1 | Select option 2.2.1 or 2.2.2
II.2.1) Bedrock is given in ```soil.dat``` | See 2.2) | 1 | 1 / Kfsu, frgw_delay, riverbed
II.2.2) If bedrock is not given in soil.dat, bedrock is assumed to be in the depth defined by meandep (or below deepest horizon given in ```soil.dat``` if the profile depth is greater than meandep) | See 2.2) | 1 | 0 / Meandep, maxdep, kfsu, frgw_delay, riverbed

<a name="table-9"></a>
**Table 9:** Groundwater option ```gw_flag=99```.

Experimental option, not verified!

Modelling options | Groundwater regime, internal representation of processes | Gw flag (```soter.dat```) | Add. parameters (```soter.dat```)
---|---|---|---
III) Groundwater in soil zone. The initial depth of the groundwater surface below soil surface is defined by gw_dist in ```soter.dat```. A separate deep groundwater storage or bedrock is ignored. | Permanent groundwater in soil zone is assumed. The groundwater level may vary in time as a function of input by percolation and drainage into river. | 99 | meandep, maxdep, riverbed, gw_dist

**3)** ```terrain.dat```<br> 
\[can be generated with the lumpR package\]

```
# Specification of terrain components
TC-ID, fraction, slope [%], position [-]
1      0.65      12         1
2      1.00       2         2
…
```

*TC-ID*:  		ID of terrain component<br>
*fraction*: 		Fraction of terrain component in LU \[-]<br>
*slope*:			Slope of terrain component \[%]<br>
*position*:	    Position of terrain component in LU (e.g. 1: highlands, 2: middle slopes, 3: lowland)<br>
(*beta_fac*):	(optional) correction factor for beta (details below)<br>
(*SDR*):         (optional) TC-specific sediment delivery ratio (details below)

Example: The terrain component with ID 1 covers 65 % (0.65) of the corresponding LU. This terrain component has a slope of 12 % and is located at the hillslope top of the landscape unit (position 1). The terrain component with ID 2 covers 100 % (1.00) of the corresponding LU, has a slope of 2 % and is located at the second position within the landscape unit (position 2), i.e. downslope of another TC.<br>
For erosion modelling, ```terrain.dat``` may contain a fifth column holding correction factors for beta (see explanation section ```beta_fac_lu.dat```) on the TC-scale OR <br>
A fifth AND sixth column holding beta_fac and a sediment delivery ratio (SDR) for each TC. If either of these are given, the respective settings for the LUs are ignored.
SDRs are applied to raw erosion on TC-scale before transport capacity limitations. Normally, they should be used with USLE and without transport capacity limitation, otherwise deposition may be accounted for twice.

**4)** ```svc.dat```<br>
\[can be generated with the lumpR package\] <br>
(optional, mandatory for sediment module and saving/loading of model states)

```
# Specifications of soil vegetation components and erosion parameters
ID,	soil_id,	veg_id	musle_k[(ton acre hr)/(acre ft-ton inch)],	musle_c[-],	musle_p[-],	coarse_fraction[%],	manning_n
11	13	21	0.13	1.0	1.0	0.8	0.011
…
```

*ID*:			    unique ID for the soil-vegetation component<br>
*Soil\_id*:	   	ID of corresponding soil unit (as specified in ```soil.dat```)<br>
*Veg\_id*:			ID of corresponding vegetation component (as specified in ```vegetation.dat```)<br>
*Musle\_k*\*:		MUSLE erodibility factor \[(ton acre hr)/(acre ft-ton inch)]<br>
*Musle\_c*\*:		MUSLE crop factor <br>
*Musle\_p*\*:		MUSLE protection factor <br>
*Coarse\_fraction*\*:	fraction of soil fragments > 2 mm \[%]<br>
*Manning\_n*\*:		Manning’s n roughness coefficient<br>

\*Each of these columns can be replicated 4 times to describe seasonal dynamics of the respective parameter. In that case, the corresponding seasonality file must be created (see ```rainy_season.dat```).

Example: The SVC with the ID 11 consists of the soil with ID 13 (as specified in ```soil.dat```) and the vegetation/landuse with ID 21 (as specified in ```vegetation.dat```), resulting in the MUSLE-factors K, C, P of 0.13, 1.0 and 1.0, respectively. It contains 0.8 % coarse particles and has a surface roughness coefficient of 0.011.

**5)** ```svc_in_tc.dat``` <br>
\[can be generated with the lumpR package\]

```
# Specification of which SVCs are contained in each TC
TC-ID[-],	SVC-ID[-],	fraction[-] 
11			12			1.0
12			13			0.2
12			14			0.8
…
```

*TC-ID*:	ID of terrain component (as specified in ```terrain.dat```)<br>
*SVC-ID*:	ID of soil vegetation component (as specified in ```svc.dat```)<br>
*Fraction*:	fraction of TC that is covered by the current SVC

Note: The sum of “fraction” over a specific TC can be smaller than one as the sum over a TC plus “fraction_rocky” needs to sum up to one (see also the description for file ```soil_vegetation.dat```). If this requirement is not fulfilled, WASA-SED will issue a warning and normalize the fraction to unity automatically.

Example: The TC with the ID 11 consists only of the SVC with the ID 12. The TC with the ID 12 is covered by the SCV 12 at 20 % of its area. The remaining 80 % of TC 12 consist of SVC 14.

**6)** ```soil_vegetation.dat```<br>
\[can be generated with the lumpR package\]

```
# Specification of soil-vegetation components (links LU, terrain component, soil and vegetation properties) 
# For each block: first line Soil IDs, Second line Land use, third line fraction of SVCs in each terrain component
Subbasin-ID[-], LU-ID[-], TC-ID[-], fraction_rocky[-], nbrSVC[-], Soil-ID(30 values)[-], Vegetation-ID (30 values)[-], fraction (30 values)[-]
49	19	25	0.12	9	86	30	77	86	85	86 …	43  in total max. 30 IDs*
49	19	25	0.12	9	8002	8004	8005	8500	7203	9203 …	7203  
49	19	25	0.12	9	0.017	0.031	0.019	0.022	0.043	0.598 …	0.025 
49	19	26	0.000	2	27	27	0	0	0	0 …	0 
49	19	26	0.000	2	8002	8004	0	0	0	0 …	0 
49	19	26	0.000	2	0.024	0.042	0	0	0	0 …	0.000 
…
```
<br>

*Subbasin-ID*:	ID of sub-basin (ID), same ordering as in ```hymo.dat```<br>
*LU-ID*:  ID of corresponding LU (as determined in ```hymo.dat```)<br>
*TC-ID*:  	ID of corresponding terrain component (as determined in ```soter.dat```)<br>
*fraction\_rocky*:	   	fraction of impermeable (rock) area in each terrain component \[-]<sup>1</sup> <br>
*nbrSVC*:	 Number of soil-vegetation components (SVCs) in current TC of  sub-basin<br>
*Soil-IDs (nbrSVC values)*: 	1<sup>st</sup> row of each block: corresponding soil-IDs as defined in ```soil.dat```<br>
*Vegetation-ID (nbrSVC values)*:	2<sup>nd</sup> row of each block: corresponding vegetation-ID as defined in ```vegetation.dat``` <br>
*fraction (nbrSVC values)*: Areal fraction of SVCs in current terrain component of current sub-basin \[-]

(\*each Sub-basin, LU, Terrain Component Unit has a block of data in three lines)<br>
<sup>1</sup> fraction_rocky will probably become obsolete in future versions. If it is 0 for all TCs, the rocky fraction is determined from the fraction of soils that have with a coarse fraction of 1 in their topsoil (see ```soil.dat```). Please note that fraction_rocky and the fraction values are internally normalized to unity if the sum is greater than one. It is advisable, however, to respect this already during the pre-processing. The package lumpR is able to handle impervious areas.

Example: The combination sub-basin ID of 49, the LU-ID of 19 and the terrain component-ID 25 has a fraction of 12 % (0.12) of impermeable rock area, and 9 different soil and landuse / vegetation classes. The sum of the fraction of the impermeable area and of the areal fractions of all SVCs must equal 1.0. The first row holds the 9 different soil-IDs (86, 30, 77, 86, etc.). The second row contains the landuse / vegetation classes for the same sub-basin – LU – terrain component combination (8002, 8004, 8005, etc.). The third line holds the areal fraction of each soil-vegetation specification within each LU-terrain combination. The next three lines contain the same block of data for sub-basin 49, LU 19 but for terrain component-ID 26. <br>
The order of the sub-basins in the first column has to follow the same order of the sub-basin IDs that was used in hymo.dat (due to computational reasons), in this case 49, 50, 1, 44, etc.; otherwise an error message occurs.

**7)** ```soil.dat```<br>
\[can be generated with the lumpR package\]

```
# Specification of soil parameters
Soil-ID[-], number(horizons)[-], res[Vol-], PWP[-], FK2.5[-], FK1.8[-], nFK[-], saturated[-], depth[mm], ks[mm/d], suction[mm], pore-size-index[-], bubblepressure[cm], coarse_frag[-], shrinks[0/1], bedrock[0/1], alluvial[0/1] (1 line)\*

    1    3  0.073  0.102  0.184  0.270  0.168  0.472   200.0    504 9.037   72.9  0.3593    9.80 0.000 0  0.102  0.160  0.265  0.349  0.188  0.472   600.0    1096.288  125.5  0.2705   16.16 0.000 0  0.064  0.084  0.152  0.229  0.145  0.472   300.0  7709.267   46.1  0.3891  6.31 0.000 0 0 0 (1 line)\*

    2    4  0.070  0.094  0.164  0.240  0.146  0.452   250.0    6229.629   67.2  0.3816    9.07 0.001 0  0.075  0.102  0.176  0.253  0.151  0.453   250.0    5250.267   71.4  0.3655    9.56 0.001 0  0.076  0.106  0.184  0.262  0.156  0.458   500.0    5368.695   82.4  0.3558   10.96 0.000 0  0.066  0.087  0.155  0.232  0.145  0.462   500.0    8344.087   56.4  0.3903    7.67 0.225 0 0 1
    
    3    3  0.075  0.104  0.191  0.278  0.174  0.396   100.0    1020.368  154.2  0.3666   20.89 0.000 0  0.089  0.126  0.215  0.297  0.170  0.434   400.0    1381.206  107.1  0.3192   14.18 0.000 0  0.070  0.097  0.181  0.270  0.173  0.434   500.0    2655.635  102.0  0.3741   13.87 0.000 0 0 0
…
```

\* the data must be all in one line

*Soil-ID*:			ID of soil unit \[-]<br>
*numb* (*horizons*):	number of horizons in soil profile \[-]<br>
*res*:		residual soil water content (*horizons*) \[Vol. fraction]<br>
*PWP*:		water content at permanent wilting point (at 15000 cm suction) (*horizons*) \[Vol. fraction]<br>
*FK2.5*:		field capacity (316 hPa / pF=2.6) (*horizons*) \[Vol. fraction] <br>
*FK1.8*:		field capacity (63 hPa / pF=1.8) (*horizons*) \[Vol. fraction]<br>
*nFK*:	usable field capacity (*horizons*) \[Vol. fraction]<br>
*saturated*:	saturated water content (*horizons*) \[Vol. fraction]<br>
*depth*:	thickness of soil horizon (*horizons*) \[mm]<br>
*ks*:		saturated hydraulic conductivity (*horizons*) \[mm/d]<br>
*suction*:	suction at the wetting front (after Green & Ampt) (*horizons*) \[mm]<br>
*pore-size-index*:	pore-size index (after Green & Ampt) (*horizons*) \[-]<br>
*bubblepressure*:		bubbling pressure (=1/alpha of van Genuchten) (*horizons*) \[cm]<br>
*coarse\_frag*:		fraction of coarse fragments (*horizons*) \[Vol. fraction]<sup>1</sup><br>
*shrinks*:		flag for soil structure (*horizons*) \[0/1, currently not used, set to 0]<br>
*bedrock*:	Does bedrock occur below deepest horizon or profile \[0: no bedrock, 1: bedrock below deepest horizon]?<br>
*alluvial*:	Is this soil an alluvial soil \[0: no alluvial soil, 1: alluvial soil]?  <br>
<sup>1</sup>If set to 1 for a top horizon, these soils can contribute to fraction\_rocky of a TC (see ```soil_vegatation.dat```).

Example: The soil class with the ID 1 has 3 horizons in the soil profile. Consecutively within each row, first all horizon-specific attributes of horizon 1 are given, then those of horizon 2 and finally those of horizon 3. For example, the first horizon has a residual soil water content of 7.3 % (0.073), a water content at permanent wilting point of 10.2 % (0.102), a field capacity of 18.4 % (0.184), a field capacity \[63 %] of 27.0 % (0.270), a usable field capacity of 16.8 % (0.168), a saturated water content of 47.2 % (0.472), a thickness of 200 mm, etc… The values for the second horizon start again with a residual soil water content of 10.2 % (0.102), a water content at permanent wilting point of 16.0 % (0.160) etc. The very last two parameters of each row are the flags for bedrock and alluvium for the entire soil column. In the first line of the file dump above, the data belonging to one horizon have the same colour. All data for each individual soil ID must be given continuously in one line.

**8)** ```soil_particles.dat``` <br>
\[can be generated with the lumpR package\]

```
# Particle size distribution of topmost horizons of soils
soil_id,	part_class_id,	fraction[-]
1		1				0.3
1		2				0.4
1		3				0.3
2		1				0.1
2		2				0.2
2		3				0.7
…
```

*Soil\_id*: ID of Soil (as determined in ```soil.dat```)<br>
*part\_class_id*:		ID of particle-size class (as determined in ```part_class.dat```)<br>
*fraction*:	mass fraction of the respective particle class in the topmost horizon \[-]

Example: The topmost horizon of soil 1 contains 30 % of particles falling into the particle-size-class 1. Furthermore, it consists of 40 % of particles out of class 2 and of 30 % of particles out of class 3.

**9)** ```vegetation.dat```<br>
\[can be generated with the lumpR package\]

```
# Specification of vegetation parameters
Veg-ID, Stomata_Resistance[s/m], minsuction[hPa], maxsuction[hPa], height[m], root-depth[m], LAI [-], albedo[-] 
9101 200 10000	30000 5	 5.5  6  6  1.5 1.5  1.7 1.8	1 5.5  5.5  1.5 0.23	0.17	0.17	0.21
9102 175 10000	30000 3 3  3  3   1.25  1.25	1.25  1.25  0.4 4 4 0.8 0.25	0.18	0.18	0.23
…
```
<br>

*Veg-ID*:			ID of landuse/vegetation<br>
*Stomata_Resistance*:	stomata resistance without water stress \[s/m]<br>
*minsuction*:	suction threshold for water stress effect on resistance (begin of stomata closure) \[hPa]<br>
*maxsuction*:	suction threshold for water stress effect on resistance (total closure of stomata – wilting point) \[hPa]<br>
*height*:	Average height of vegetation canopy \[m] 4 values<br>
*root depth*:	Rooting depth of vegetation \[m] 4 values<br>
*LAI*:	Leaf area index of vegetation cover \[-] 4 values<br>
*albedo*:	Surface albedo \[-] 4 values

Example: The landuse/vegetation class with ID 9101 has a stomata resistance without water stress of 200 s/m, a minimal suction threshold of 10000 hPa and a maximal suction threshold for water stress effect on resistance of 30000 hPa. Furthermore, 4 values for vegetation height (5 m, 5.5 m, 6 m, 6 m), 4 values for root depth (1.5 m, 1.5 m, 1.7 m, 1.8 m), 4 values for leaf area index (1, 5.5, 5.5, 1.5) and 4 values for albedo (0.23, 0.17, 0.17, 0.21) are specified for this landuse/vegetation class. Do not use zero values here but very low ones instead! The four values for the four parameters (height, root depth, LAI, albedo) reflect the temporal changes of vegetation parameters as a function of seasonal changes (e.g. due to a rainy season). The first value of a set of four reflects the vegetation properties before the rainy season, the second value at the beginning of the rainy season, the third value at the end, and the fourth value after the rainy season. The exact date for these four points in time is given in the next file (```rainy_season.dat```).

**10)** ```rainy_season.dat``` (optional)<br>
```k_seasons.dat``` (optional)<br>
```c_seasons.dat``` (optional)<br>
```p_seasons.dat``` (optional)<br>
```coarse_seasons.dat``` (optional)<br>
```n_seasons.dat``` (optional)

```
# Specification of the rainy season (per year)
# for the interpolation of temporal distribution of vegetation characteristics (Rootdepth,height,lai,albedo)
subbasin_id veg_id	year	DOY1	DOY2	DOY3	DOY4
49	2	1980	10	40	175	205
49	-1	1980	20	40	175	205
50	-1	1980	10	40	175	205
-1	1981	5	35	150	180
-1	-1	-1	36	151	181
…
```

*subbasin_id*: ID of sub-basin<br>
*Year*:	Simulation year(s), each individual year has to be listed in the file<br>
*Veg\_id*: (optional column) option for vegetation-specific dynamics<br>
*DOY1 – DOY4*: day-of-year (1-365/366) for four days within the seasonal cycle, serving as support points for intraannual dynamics (e.g. start of transition phase, first day of rainy season, last day of (rainy) season, start of transition phase)

Temporal dynamics of LAI, root depth, C-factor, etc. can be specified as a piecewise linear cycle with four segments. The ordinates (y-values, e.g. minimum, maximum and transitional LAIs) are listed in ```vegetation.dat``` and ```svc.dat```. They are time-invariant.

The respective abscissa (x-values as day-of-year, DOY) are listed in ```*_seasons.dat```. They can be specified separately for each year, sub-basin and (optionally) vegetation class.

Eg. ```rainy_season.dat``` contains the points in time that serve as temporal nodes, between which the seasonal dynamics of LAI, vegetation height, root depth and albedo are linearly interpolated. The nodes are specified as julian days/DOYs in relation to the respective year. Negative values (previous) and values greater than 365/366 (next year) are allowed as long as they do not surpass the next adjacent node. Please make sure that all required rainy-seasons are specified (i.e. for all sub-basins, all simulation years), otherwise an error message occurs. If not specified otherwise, values for days BEFORE the first / AFTER the last specified node are extrapolated with constant value.

“-1” can be used as a wildcard for sub-basins, vegetation classes and years. The file is interpreted top to bottom, i.e. if multiple matches are found, the uppermost is used. Thus, a fallback specification (“all other cases”) should go into the last line (see example).
A standard (legacy) rainy_season.dat can be generated from rainfall data using the function rainy_season of the lumpR-package. If rainy_season.dat is missing or contains only the headerlines, only the first value in vegetation.dat is used (i.e. no seasonality).

The optional files ```k_seasons.dat```, ```c_seasons.dat```, ```p_seasons.dat```, ```coarse_seasons.dat```, ```n_seasons.dat``` work analogously for the USLE factors K, C, P, the coarse fraction and Manning’s n. If any of these files is found, four respective columns (instead of the default one) must be given in ```svc.dat```. Any seasonality file must reside in the input-subfolder ```hillslope```.

Example: The example file above is for two specific simulation years: 1980 and 1981. For vegetation class 2 in sub-basin 49, the rainy season in 1980 started on the 40th day (09.02.1980), and stopped on the 175th day (28.06.1980) (line 4). For all other vegetation classes, it started ten days later (DOY 20) (line 5). For sub-basin 50, the dynamics are identical for all vegetation classes (line 6).

In 1981, the same DOYs are used for all sub-basins (line 7). All other sub-basins, vegetation classes or years would use the values of the last line.

**11)** ```scaling_factor.dat```<br>
(optional)

```
Subbasin-ID	mean_kf-calib-factor
1	10 
2	4.5 
…
```

*Subbasin-ID*:		ID of sub-basin<br>
*mean\_kf-calib-factor*:	scaling factor (actually a divisor)

This file is optional and is only read if doscale (in ```do.dat```) is set to “.T.”. In this case, ```scaling_factor.dat``` is expected in the subdirectory ```Others/```.

Example: All values for saturated hydraulic conductivity during infiltration in sub-basin 1 are modified by dividing them by 10. Maximum interception is adjusted by multiplication with 1/(0.340+0.647*10).

**12)** ```calibration.dat```<br>
(optional)

```
Soil-ID,	Ksat-calib-factor
-1	20
5	0.3
…
```

*Soil-ID*:	ID of soil<br>
*Ksat-calib-factor*:	calibration factor 

This file is optional and is in the subdirectory ```Others/```.

Example: All values for saturated hydraulic conductivity of ALL soils (wildcard “-1”) are increased by factor 20. For soil 5, the values are additionally multiplied with 0.3.

**13)** ```erosion.ctl```<br>
(optional)

```
# WASA-control file for erosion and sediment transport routines;
application_scale	0	
erosion_equation	4	
ri_05_coeffs	1.911	0.807
transport_limit_mode	2	
#river module
transp_cap_a	0.016111	
transp_cap_b	1.807	
```

*Application scale*:		0: apply equations on TC-scale; 1: apply on sub-basin-scale<br>
*Erosion equation*:		erosion equation to be used 1: USLE, 2: Onstad-Foster, 3: MUSLE, 4: MUST<br>
*ri\_05\_coeffs*:	needed for USLE and OF: coefficients for estimation of maximum half-hour rainfall intensity (ri\_05) from daily/hourly rainfall data (R\_dt): ri\_05=a\*R\_dt<sup>b</sup>. <br>
*transport\_limit\_mode*:		different modes how/if transport capacity of runoff is limited: 1: no transport capacity limit; 2: transport capacity according to Everaert (1991); 3: transport capacity computed from MUSLE with maximum erodibility<br>
*transp\_cap\_a*:	empirical factor for computing suspended sediment transport capacity in river (a \* vel_peak \*\* b)<br>
*transp\_cap\_b*:		empirical factor for computing suspended sediment transport capacity in river (a \* vel_peak \*\* b)

This file and any of its entries are optional. If not present, default values are used (Application scale=0; Erosion equation=3, ri\_05\_coeffs = (a\_i30=1.1630; b\_i30=0.667981 for daily resolution;  a\_i30=1; b\_i30=1 for hourly resolution); transport\_limit\_mode=2, transp\_cap\_a=0.016111, transp\_cap\_b=1.707).

**14)** ```frac_direct_gw.dat```<br>
(optional)

```
0.7
```

This file contains a single value x that specifies the fraction of the groundwater (formed in the LUs) that is routed directly into the river (as is the default). The remaining fraction 1-x enters the lowermost TC as subsurface flow. Low values of x tend to reduce periods of very low flow in ephemeral rivers. Default x is 1.

```frac_direct_gw.dat``` resides in the root of the WASA-SED input-directory as specified in ```do.dat```.

**15)** ```beta_fac_lu.dat```<br>
(optional)

```
# specified correction factor for beta (rill/interrill ratio, used for the computation of L-factor, see Renard, 1997, pp.101
lu_id,	beta_factor
12111	2
21111	0.5
…
```

*Lu\_id*:		landscape unit ID<br>
*beta\_factor*:	factor for modification of beta

This file allows specifying a correction factor for beta for selected LUs. Beta describes the ratio of
rill/interrill erosion and is used for the computation of the L-factor (see Renard et al., 1997, pp.101, eq. 4-1 – 4-3).

Default value (for non-specified LUs) is 1. Common values are 0.5 for a low (lower yield) and 2 for a high (higher yield) rill/interrill ratio. If this correction factor is already specified for the TC-scale in ```terrain.dat```, the values in ```beta_fac_lu.dat``` are ignored.

A row with an lu_id of -1 (wild card) will set all unset LUs to the specified value.

**16)** ```sdr_lu.dat```<br>
(optional)

```
# Prespecified Sediment delivery ratio [0-1] for Landscape units		
lu_id,	sdr
11211	0.94
11311	0.94
…
```

*lu\_id*:		landscape unit ID<br>
*sdr*:		sediment delivery ratio \[0...1]

This file allows specifying a sediment delivery ratio for selected LUs. 

Default value (for non-specified LUs) is 1. Check!

If this correction factor is already specified for the TC-scale in ```terrain.dat```, the values in ```sdr_lu.dat``` are ignored. A row with an lu_id of -1 will set all unset LUs to the specified value.

Warning: Using SDR should be used without a transport capacity limitation, otherwise, deposition is considered twice.

**17)** ```calib_wind.dat```<br> 
(optional) 

This file contains a single value which will be used as static wind speed value (in m/s) within the model. If this file is not given, a value of 1 m/s is used by default. As this is a very sensitive parameter, it can be used for calibration of evapotranspiration.

**18)** ```snow_params.ctl```<br>
(optional)

The two logical parameters do_rad_corr and do_alt_corr allow controlling, whether radiation correction for aspect and slope, and height-depended temperature modifications, respectively, are applied.

```
#WASA-control file for snow routines;
a0	0.002	#empirical coefficient (m/s); linear dependence of turbulent transfer coefficient (D) in sensible heat flux: D = a0 + a1*WindSpeed
a1	0.0008	#empirical coefficient (-)  ; linear dependence of turbulent transfer coefficient (D) in sensible heat flux: D = a0 + a1*WindSpeed
kSatSnow	0.00004	#Saturated hydraulic conductivity of snow (m/s)
densDrySnow	450	#Density of dry snow (kg/m│)
specCapRet	0.05	#Capill. retention volume as fraction of solid SWE (-)
emissivitySnowMin	0.84	#Minimum snow emissivity used for old snow (-)
emissivitySnowMax	0.99	#Maximum snow emissivity used for new snow (-)
tempAir_crit	0.2	#Threshold temperature for rain-/snowfall (░C)
albedoMin	0.55	#Minimum albedo used for old snow (-)
albedoMax	0.88	#Maximum albedo used for new snow (-)
agingRate_tAirPos	0.00000111	#Aging rate for air temperatures > 0 (1/s)
agingRate_tAirNeg	0.000000462	#Aging rate for air temperatures < 0 (1/s)
soilDepth	0.1	#Depth of interacting soil layer (m)
soilDens	1300	#Density of soil (kg/m3)
soilSpecHeat	2.18	#Spec. heat capacity of soil (kJ/kg/K)
weightAirTemp	0.5	#Weighting param. for air temp. (-) in 0...1
lat	42.4	#Latitude of centre of study area
lon 0.55	#Longitude of centre of study area
do_rad_corr .TRUE. #modification of radiation with aspect and slope
do_alt_corr=.TRUE.   #modification of temperature with altitude of LU
```

**19)** ```lu2.dat```<br>
(optional)

Required when using the snow module. Hold LU-specific parameters.

```
Specification of landscape units, snow parameters		
LU-ID[id]	aspect[deg]	mean_altitude_over_subbas_mean[m]
1	104.620873988632	-169.42
2	-34.3803447238449	-111.23
3	-123.69006752598	-80.39
4	-6.84277341263094	-47.46
```

<a name="3-3-input-files-for-the-river-module"></a>
### 3.3 Input files for the river module

The input files for the river module are located in the folder ```Input\\[case_study]\River``` and are summarised in [Table 10](#table-10). Three options are available for the river routing: routing scheme 1 comprises the original river routing using time response functions, routing scheme 2 uses the Muskingum routing and suspended sediment transport and routing scheme 3 uses the Muskingum routing and bedload transport. Routing schemes 2 and 3 enable a spatially distributed representation of river stretch characteristics. Sediment-transport calculations are only possible for routing schemes 2 and 3. The flow calculations are carried out in routing order, i.e. the river stretches which are located most upstream are calculated first. The routing order is specified in ```routing.dat```. The key model input parameters for water and sediment routing are stored in an input file called ```river.dat``` that assigns each sub-basin with a specific map ID a corresponding river stretch. The input file ```response.dat``` contains the time response parameters that were used for the original version of the WASA code (routing scheme 1).

<a name="table-10"></a>
**Table 10:** Input data files for the river component.

Parameter File | Content
---|---
```routing.dat```	|	Specification of routing order (for all routing schemes)
```river.dat```	| Specification of river parameters (for routing scheme 2 and 3)
```response.dat``` |	Specification of time response parameters (for routing scheme 1)
```bedload.dat```	| Specification of bedload data (for routing scheme 3)
```subbasin_out.dat``` (optional) | pre-specification of outflow of selected sub-basins
```subbasin_outsed.dat```	(optional) | pre-specification of sediment output of selected sub-basins
```transposition.dat``` (optional)	| Specification of additional water fluxes between sub-basins
```river_storage.stat``` (optional) | Initialisation state variables, see [section 3.6](#3-6-state-variables)

<br>

**1)** ```routing.dat```

```
# Specification of routing order (flow directions)
No., Subbasin-ID(upstream), Subbasin-ID(downstream)
1	4	10
2	15	10
3	39	10
4	10	50
5	1	50
6	44	50
7	50	49
8	49	3
9	3	29
10	29	999
```

*No.*: Continuous numbering (calculation order, value is ignored) <br>
*Subbasin-ID(upstream)*: ID of sub-basin, which is located upstream of another sub-basin <br>
*Subbasin-ID(downstream)*:	ID of sub-basin, which is located downstream of the previous sub-basin

Example: This file defines the order of the calculation of the sub-basins. All sub-basins must be listed before their corresponding outlet basins, otherwise an error is issued. For example, sub-basin No. 4 is upstream of sub-basin No. 10. Sub-basins No. 15 and 39 are also upstream of No. 10. The runoff of sub-basin No. 10 flows into sub-basin No. 50 etc. The sub-basin at the outlet of the entire drainage system must drain to a sub-basin labelled 999 or 9999.

**2)** ```river.dat```

```
# Specification	of	river	parameters			
Subbasin-ID, depth(m), width(m), side ratio (m/m), bottom width of floodplain (m), side ratio floodplains (m/m), channel slope(m/m), length(km), manningn (-), manningn_floddplain (-) Ksat(mm/hr), erodibilityfactor(-),coverfactor(-), riverbedrock(-),baseflowalphafactor(days), msk_x(-), Q_spring(m3/s) 
1	1	5	2	100	4	0.006	7.4	0.02	0.05	25	0.1	1	0	0.1	0.2	4	0.1
```
<br>

*Subbasin-ID*:		ID of sub-basin <br>
*depth*:			Bankful depth of river reach \[m] <br>
*width*:			Bankful width of river reach \[m] <br>
*side ratio*:		Run to rise ratio of river banks (1/side channel slope) \[m/m] <br>
*bottom width*:		Bottom width of flood plain \[m] <br>
*side ratio floodplains*:	Run to rise ratio of floodplain side slopes (1/floodplain side slope) \[m] <br>
*slope*:			Slope of river reach \[m/m] <br>
*length*:			Length of river reach \[km] <br>
*manningn*:		Manning’s n of river reach \[-] <br>
*manningn_floodplain*:	Manning’s n of floodplain reach \[-] <br>
*Ksat*: Saturated hydraulic conductivity the river bed \[mm/h] <br>
*erodibility factor*:	River erodibility factor of river reach \[-] <br>
*cover factor*: River vegetation cover factor of river reach \[-] <br>
*riverbedrock*: Cover of solid rock in riverbed \[-] <br>
*baseflow*: baseflow alpha factor for bank storage (days) <br>
*msk\_x*: Muskingum X weighting coefficient \[-] <br>
*msk\_k*: Muskingum K storage time constant \[hours] <br>
*Q\_Spring*: Initial conditions for headwater reaches (minimum discharge) \[m<sup>3</sup>/s]

Example: The river stretch at the sub-basin with the ID of 1 has a bankful depth of 1 m, a width of 5 metres, a site ratio of 2, a bottom width of the floodplain of 100 m, a side ratio on the floodplains of 4, a channel slope of 0.006 (or 0.6 %), a length of 7.4 km, a Manning’s n of 0.02 and a Manning’s n in the floodplain of 0.05, a Ksat of 25 mm/h, an erodibility factor of 0.1, a cover factor of 1, a riverbedrock factor of 0, a baseflowalphafactor of 0.1 days, a Muskingum X coefficient of 0.2, a Muskingum K factor of 4 hours and an initial condition of 0.1 m<sup>3</sup>/s. The dimensions of the trapezoidal channels including the floodplains are depicted in [Figure 2](#figure-2). The height of the wedge at the channel bottom (enables smooth transition of low flows) is fixed to 0.1 m.
 
 <a name="figure-2"></a>
**Figure 2:** Trapezoidal channel dimension with floodplains.
![Trapezoidal channel dimension with floodplains](./Figure_2.svg)

**3)** ```response.dat```
```
# Specification of routing parameter
Subbasin-ID,lag time [d],retention[d]
49	0.5	2.0
50	1	1.5
1	4	7
…
```

*Subbasin-ID*:	ID of sub-basin<br>
*lag time*: Lag time between runoff input to sub-basin and first runoff response at its outlet in \[days]<br>
*retention*:	Retention specifies the maximum retention time in the sub-basin in \[days]

Reference is midday, partial coverage of days is considered. Autochtonous runoff (riverflow generated inside a sub-basin, not entering from upstream) is routed slightly different with zero lag time (triangular like this: /\ \_; tL\*=0, tR\*=tL+tR).
Example: The sub-basin with the ID of 49 has a lag time of 0.5 days and a retention time of 2 days (i.e. its runoff will be delayed by 0.5 day, then stretched over another 2 days). The sub-basin with the ID of 50 has a lag time of 1 day and a retention time of 1.5 days; etc.

For a detailed description of the routing process and the linear response function, see [Güntner (2002)](#guentner-2002), p. 48 and [Bronstert et al. (1999)](#bronstert-et-al-1999). The order of the sub-basins in the first column has to follow the same order of the sub-basin IDs as was used in ```hymo.dat``` (due to computational reasons); otherwise an error message occurs.

**4)** ```bedload.dat```

```
# Specification of bedload modelling parameter D50
Subbasin-ID, D50 (m)
49	0.048
50	0.048
1	0.048
…
```

*Subbasin-ID*: ID of sub-basin <br>
*D50*: median sediment particle size in the riverbed in (m)

Example: The river stretch in the sub-basin with the ID of 49 has a riverbed gradation with a D50 value of 0.048 m.

**5)** ```subbasin_out.dat``` <br>
(optional)

```
# Pre-specified mean daily river flow [m3/s] for selected sub-basins (MAP-IDs)		
Date	Timestep	Sub-basin-ID.
0	0	4
1092005	1	0.5
2092005	1	0.3
3092005	1	0.2
…
```

*Date*: Date in the format ddmmyyyy (leading zeroes are optional) <br>
*Timestep*: timestep (not interpreted in daily resolution, 1..24 for hourly resolution)<br>
*Subbasin-ID*: ID of sub-basin

This optional file allows specifying the water output of selected sub-basins. If this file is not found in the folder ```Time_series```, all sub-basins are treated regularly. Otherwise, any outflow that is specified in this file is used directly as an output of the respective sub-basin – no computations are performed within this sub-basin (evaporation, groundwater, river routing, reservoir, etc.). Consequently, no climate input needs to be specified for such subbasins. WASA reads data from this file sequentially, starting from start of simulation and every calendar year (e.g. chunks of 365 days). Warning: The subsequent rows are assumed without gaps and not checked for completeness in the time series. "-1" is regarded as "no data" and will lead to "no data" in the riverflow in all affected downstream subbasins.

Example: Sub-basin 4 has pre-specified discharge of 0.5 m<sup>3</sup>/s for 1 Sep 2005.

**6)** ```subbasin_outsed.dat``` <br> 
(optional)

```
# Pre-specified daily sediment output [t] for selected  sub-basins (MAP-IDs), mean PSD: 0.3 0.2 0.5
Date	Timestep	Sub-basin-ID.
0	0	4
1092005	1	0.5
2092005	1	0.3
3092005	1	0.2
…
```

*mean PSD*:	mean particle size distribution to be used for all records, consists of n_sed_classes fraction values that sum to 1 <br>
*Date*: Date in the format ddmmyyyy (leading zeroes are optional)<br>
*Timestep*:	timestep (not interpreted in daily resolution, 1..24 for hourly resolution) <br>
*Subbasin-ID*: ID of sub-basin

This optional file allows specifying the sediment output of selected sub-basins. WASA expects this file in the folder Time_series. If this file is not found, all sub-basins are treated regularly. Otherwise, any sediment output that is specified in this file is used directly as an output of the respective sub-basin – no sediment related computations are performed within this sub-basin. WASA reads data from this file sequentially, starting from start of simulation and every calendar year (e.g. chunks of 365 days). Warning: The subsequent rows are assumed without gaps and not checked for completeness in the time series. "-1" is regarded as "no data" and will lead to "no data" in the sediment flux in all affected downstream subbasins.

Example: Sub-basin 4 has pre-specified sediment output of 0.5 t/d for 1 Sep 2005, distributed among 3 particle size classes with the fractions 0.3, 0.2 and 0.5.

**7)** ```transposition.dat``` <br>
(optional)

```
# Water transpositions via canals or pipes between sub-basins, in order of routing scheme
Start-Subbasin-ID, Flag(reservoir/river), Flow(m3/s,) Loss(%), Destination-Subbasin-ID, Flag (reservoir/river), begin_year
        91         1      1.75      0.01        96         1      1965
        67         2       2.4       0.1        31         1      1997
        31         1         2      0.06        30         1      1993
        30         1         2      0.02        29         2      1993
…
```

*Start-Subbasin-ID*:		ID of sub-basin (source of water abstraction)<br>
*Flag(reservoir/river)*:		water abstraction from: 1 (reservoir) or 2 (river)<br>
*Flow(m<sup>3</sup>/s)*:			amount of re-routed water<br>
*Loss(%)*:			transmission loss<br>
*Destination-Subbasin-ID*:		ID of sub-basin (destination of water abstraction)<br>
*Flag (reservoir/river)*:		water diversion to: 1 (reservoir) or 2 (river)<br>
*begin\_year*:			start time of water abstraction

Currently only works for daily resolution. For suspended sediments, abstraction from reservoirs assumes zero concentration, whereas abstraction from river uses river concentration.
Abstractions are taken from the outlet of river reaches, and added to the inlet points of reaches.
This file is optional. It is only read if ```dotrans``` (in ```do.dat```) is set to ```.true.```. In this case, ```transposition.dat``` is expected in the subdirectory ```Others/```. <br>

<a name="3-4-input-files-for-the-reservoir-module"></a>
### 3.4 Input files for the reservoir module

The input files for the reservoir module are located in the folder ```Input\\[case_study]\Reservoir``` and are summarised in [Table 11](#table-11). The files listed below are required according to the simulation option defined in the file ```do.dat```. Reservoirs are considered in the model simulations if the option ```doreservoir``` is switched on. For simulations of reservoir water balance the file ```reservoir.dat``` (file 1) is required. Nevertheless, additional files can be given to improve the model results (files 2 to 6). For calculations of reservoir sediment balance, the options ```doreservoir``` and ```dosediment``` must be switched on. The reservoir sedimentation model consists of two modelling approaches, which may be applied according to reservoir size and data availability. For reservoirs with information about their geometric features (reservoir topography, stage-area and stage-volume curves) and physical properties of sediment deposits, such as deposition thickness, grain size distribution of sediment deposits and sediment densities, a detailed modelling approach to reservoir sedimentation may be applied (files 7 to 9 are required; and files 10 to 12 are used to improve model results). For reservoirs without those characteristics, a simplified modelling approach is used (file 8 is required). Networks of small reservoirs are considered in the model simulations if the option doacudes is switched on. For simulations of water and sediment routing through the reservoir networks the file 13 and 16 are required (files 14, 15 and 17 are used to improve model results).

<a name="table-11"></a>
**Table 11:** Input data files for the reservoir component.

Nr. |Parameter File	| Content
---|---|---
1 | ```reservoir.dat``` | Specification of reservoir parameters
2 | ```lateral_inflow.dat``` (optional) | Specification of lateral inflow into the sub-basin’s reservoir (multiple subbasins draining into one reservoir)
3 | ```operat_rule.dat``` (optional) | Specification of reservoir operation rule
4 | ```operat_bottom.dat``` (optional) | Specification of operation rule of reservoir bottom outlets
5 | ```cav.dat``` (optional) | Specification of stage-area and stage-volume curves of the sub-basin’s reservoir
6 | ```intake.dat``` (optional) | Specification of measured data on regulated outflow discharge through intake devices from the sub-basin’s reservoir. This file is a time series file and thus need to be in directory Time_series!
7 | ```hydraul_param.dat``` (optional)	| Specification of hydraulic parameters of the sub-basin’s reservoir
8 | ```sed.dat``` | Specification of sedimentation parameters of the sub-basin’s reservoir
9 | ```cross_sec_”ID”.dat``` (optional) | Specification of cross section geometry of the sub-basin’s reservoir (sub-basin with a specific ID)
10 | ```original_sec_”ID”.dat``` | (optional)	Specification of original cross section geometry of the sub-basin’s reservoir (sub-basin with a specific ID)
11 | ```sizedist_”ID”.dat``` (optional) | Specification of size distribution of original bed material along the cross sections of the sub-basin’s reservoir (sub-basin with a specific ID)
12 | ```main_channel.dat``` (optional) | Specification of main channel geometry of the sub-basin’s reservoir
13 | ```lake.dat``` | Specification of parameters for the reservoir size classes
14 | ```lake_maxvol.dat``` (optional) | Specification of water storage capacity for the reservoir size classes
15 | ```lake_year.dat``` (optional) | Specification of changes on the number of reservoirs in the size classes
16 | ```lake_number.dat``` | Specification of total number of reservoirs in the size classes
17 | ```lake_frarea.dat``` (optional) | Specification of runoff contributing area for the reservoir size classes
18 | ```lake_storage.stat, reservoir_storage.stat``` (optional) | Initialisation state variables, see [section 3.6](#3-6-state-variables)


**1)** ```reservoir.dat```
\[can be generated with the lumpR package\]

```
Specification of reservoir parameters
Subbasin-ID, minlevel[m], maxlevel[m], vol0([1000m**3]; unknown=-999), storcap[1000m**3], damflow[m**3/s], damq_frac[-], withdrawal[m**3/s], damyear[YYYY], maxdamarea[ha], damdead[1000m**3], damalert[1000m**3], dama[-], damb[-], qoutlet[m**3/s], fvol_bottom[-], fvol_over[-], damc[-], damd[-], elevbottom[m]
60	413.30	447.67	45213.92	91795.66	36.00	1.00	0.020	1980	718.67	4802.95	45213.92	20.935	0.716	146.84	1.00	0 	300	1.5	430
```

*Subbasin-ID*: ID of sub-basin <br>
*minlevel*: Initial minimum level in the sub-basin’s reservoir \[m]. Value varies because of the sediment accumulation <br>
*maxlevel*: Maximum water level in the sub-basin’s reservoir \[m] <br>
*vol0*: Initial volume of the sub-basin’s reservoir \[10³ m³]. If set to ```-999```, 20 % of ```storcap``` is assumed. Value varies because of the sediment accumulation <br>
*storcap*: Initial storage capacity in the sub-basin’s reservoir \[10³ m³]. Value varies because of the sediment accumulation <br>
*damflow*: Target release through the barrage’s intake devices of the sub-basin’s reservoir \[m<sup>3</sup>/s]. Contributes to intake in output file ```res_*_watbal.out```. Influenced by damq\_frac and further reduced if the reservoir’s actual storage volume is less than damalert. <br>
*damq\_frac*: Maximum fraction of damflow which is released from the sub-basin’s reservoir \[-] <br>
*withdrawal*:	Water withdrawal discharge to supply the water use sectors in the sub-basin’s reservoir \[m<sup>3</sup>/s]. Outflow discharge through the dam is not considered <br>
*damyear*: Year of construction of the dam in the sub-basin <br>
*maxdamarea*: Initial maximum area of the sub-basin’s reservoir \[ha]. Value varies because of the sediment accumulation <br>
*damdead*: Initial dead volume of the sub-basin’s reservoir \[10<sup>3</sup> m<sup>3</sup>]. Value varies because of the sediment accumulation. Influences actual values of damflow and qoutlet. <br>
*damalert*: Initial alert volume of the sub-basin’s reservoir \[10<sup>3</sup> m<sup>3</sup>]. Value varies because of the sediment accumulation. Influences damflow. Should be set to volume at the height of the barrage’s intake devices.<br>
*dama, damb*:	Parameters of the area-volume relationship in the sub-basin’s reservoir (area=dama.Vol<sup>damb</sup>) \[-]. Values of reservoir area and volume are expressed in m<sup>2</sup> and m<sup>3</sup>, respectively. These are not the Molle parameters, but can be derived thereof: dama = (1/k * (alpha * k)^(alpha/(alpha-1)))^((alpha-1)/alpha); damb = (alpha-1)/alpha (damb should be > 0 and <= 1).<br>
*qoutlet*:	Maximum outflow discharge released through the bottom outlets in the sub-basin’s reservoir \[m<sup>3</sup>/s]. Contributes to qbottom in output file ```res_*_watbal.out```. Set to zero if less than damdead or less than fvol_bottom\*storcap. <br>
*fvol\_bottom*: Fraction of storage capacity that indicates the minimum storage volume for sediment release through the bottom outlets of the sub-basin's reservoir \[-]. Influences the actual value of qoutlet. <br>
*fvol\_over*: flag to simulate the retention of reservoir overflow during spillway operation \[0 = without time delay; 1 = with time delay] <br>
*damc, damd*: Parameters of the spillway rating curve in the sub-basin’s reservoir (Qout=damc.Hv<sup>damd</sup>) \[-]. Values of water height over the spillway (Hv) and overflow discharge (Qout) are expressed in m and m<sup>3</sup>/s, respectively <br>
*elevbottom*: bottom outlet elevation of the sub-basin's reservoir \[m]. Currently not used, fill in dummy values.


Example: At the outlet point of the sub-basin with the ID 60, there is a reservoir with an initial minimum level of 413.30 m, a maximum water level of 447.67 m, an initial volume of 45,213,920 m³, an initial storage capacity of 91,795,660 m³, a target outflow discharge of 36 m³/s, a water withdrawal discharge to supply the water use sectors of 20 L/s, year of construction in 1980, an initial maximum area of 718.67 ha, an initial dead volume of 4,802,950 m³, an initial alert volume of 45,213,920 m³, an area-volume relationship with parameters dama and damb set to 20.935 and 0.716, respectively, an maximum outflow discharge through the bottom outlets of 146.84 m³, a percentage of storage capacity for the operation of bottom outlets of 100 % and a percentage of storage capacity for overflow discharge through radial gates of 80 %, a spillway rating curve with parameters damc and damd set to 300 and 1.5, respectively. Value of vol0 set to -999 means that at the beginning of the simulation period the water volume is 20 % of the storage capacity ([Güntner, 2002](#guentner-2002)). Value of damq_frac set to -999 means that the reservoir operation rule is affected by irrigation season. Thus, an additional file has to be provided, which gives the interannual variability of exploitation regime (see below the file ```operat_rule.dat```). Values of damflow may be replaced by measurements when providing the optional input file ```intake.dat``` (see file description). Value of fvol_bottom set to -999 means that an additional file must be provided with detailed information about the sediment management technique selected to routing sediment through the sub-basin’s reservoir (see below the file ```operat_bottom.dat```). The order of the sub-basins in the first column has to follow the same order of the sub-basin IDs as was used in ```hymo.dat``` (due to computational reasons); otherwise, an error message occurs. Sub-basins without outlet reservoirs must not be entered in the file.

**2)** ```lateral_inflow.dat```
(optional)

```
# Specification of lateral inflow into the sub-basin’s reservoir
Subbasin-ID, reservoir_down[-]
15	60
```

*Subbasin-ID*:		ID of sub-basin with generated runoff flowing directly into the reservoir of another sub-basin <br>
*Reservoir_down*:	ID of sub-basin with an outlet reservoir that receives lateral inflow coming from another sub-basin.

Example: This optional file allows specifying lateral inflow into the sub-basin’s reservoir (multiple subbasins draining into one reservoir). The reservoir located at the outlet point of the sub-basin 60 receives lateral inflow coming from the sub-basin 15. 
Sub-basins draining regularly into their respective "own" reservoir or downstream basin must not be entered in the file.

**3)** ```operat_rule.dat``` <br>
(optional)

```
# Specification of reservoir operation rule
Subbasin-ID, dayexplot(4 values)[-], damq_frac_season(4 values)[m**3/s]
60	59	120	212	335	0.50	0.72	0.38	0.17
```

*Subbasin-ID*: ID of sub-basin <br>
*dayexplot*:	Days of change in exploitation regime in the sub-basin's reservoir \[-]. Four days of the year have to be provided <br>
*damq\_frac_season*:	Fraction of Q90 released from the sub-basin's reservoir in different seasons in the sub-basin's reservoir \[-]

Example: This optional file allows specifying changes on the operation rule of the sub-basin’s reservoir. If this file is not found in the folder reservoir, a target value of controlled outflow discharge given in the file ```reservoir.dat``` is applied to the respective sub-basin’s reservoir. The reservoir located at the outlet point of the sub-basin with the ID 60 changes its exploitation regime on those days of the year (59, 120, 212 and 335). Such values are followed by four corresponding values of fraction of Q90 released from the sub-basin's reservoir, according to the given intervals (0.50, 0.72, 0.38 and 0.17). It means that exemplarily between the 120th and the 212th of the year 38 % of Q90 can be released. The order of the sub-basins in the first column has to follow the same order of the sub-basin IDs as was used in ```hymo.dat``` (due to computational reasons); otherwise, an error message occurs. Sub-basins without outlet reservoirs or those without data on reservoir operation rule must not be entered in the file.

**4)** ```operat_bottom.dat``` <br> 
(optional)

```
# Specification of operation rule of reservoir bottom outlets
Subbasin-ID, operat_start[-], operat_stop[-], operat_elev[m]
60	270	320	430
```

*Subbasin-ID*: ID of sub-basin <br>
*operat\_start*: Target day of year to open the bottom outlets \[-] <br>
*operat\_stop*: Target day of year to close the bottom outlets \[-] <br>
*operat\_elev*: Target water depth of the sub-basin's reservoir during the period the bottom outlets remain open \[m]

Example: This optional file allows specifying the operation rule of bottom outlets of the sub-basin’s reservoir. If this file is not found in the folder reservoir, a target value of controlled outflow discharge through the bottom outlets given in the file ```reservoir.dat``` is applied to the respective sub-basin’s reservoir. The bottom outlets of the reservoir located at the outlet point of the sub-basin with the ID 60 are opened on 120th day and closed on 212th of the year. During that period, the reservoir water level should not surpass the elevation of 430 m in order to increase flow velocity and, consequently, sediment release. The order of the sub-basins in the first column has to follow the same order of the sub-basin IDs as was used in ```hymo.dat``` (due to computational reasons); otherwise, an error message occurs. Sub-basins without outlet reservoirs or those without data on operation rule of bottom outlets must not be entered in the file.

**5)** ```cav.dat``` <br>
(optional)

```
# Specification of stage-area and stage-volume curves of the sub-basin’s reservoir
Subbasin-ID, nbr. points, 1st row: elevation [m], 2nd row: reservoir area [1000m**2], 3rd row: reservoir volume [1000m**3]
60	36	413.30	415.00	416.00	417.00	…	447.00	447.67	448.00		in totalmax 200 IDs
60	36	0.00	54.82	96.10	142.89	…	6980.25	7186.73	7288.43		in totalmax 200 IDs
60	36	0.00	35.78	111.24	231.23	…	87049.73	91795.66	94184.06		in totalmax 200 IDs
```

*Subbasin-ID*: ID of sub-basin <br>
*Nbr. points*: Number of points from the stage-area and stage-volume curves of the sub-basin's reservoir <br>
*1st row*: elevation	1st row of each sub-basin: water elevation from the stage-area and stage-volume curves of the sub-basin's reservoir \[m] <br>
*2nd row*: reservoir area	2nd row of each sub-basin: reservoir area for a given elevation at the stage-area and stage-volume curves of the sub-basin's reservoir \[103 m2] <br>
*3rd row*: reservoir volume	3rd row of each sub-basin: reservoir volume for a given elevation at the stage-area and stage-volume curves of the sub-basin's reservoir \[10<sup>3</sup> m<sup>3</sup>]

Example: This optional file allows specifying the stage-area and stage-volume curves of the sub-basin’s reservoir. If this file is not found in the folder reservoir, an area-volume relationship given in the file ```reservoir.dat``` is applied to the respective sub-basin’s reservoir. The reservoir located at the outlet point of the sub-basin with the ID 60 has 36 points at the stage-area and stage-volume curves. The first row holds 36 values of water elevation at the stage-area and stage-volume curves (413.30 m, 415.00 m, 416.00 m, 417.00, etc). The second row holds 36 values of reservoir area for the given values of elevation at the stage-area and stage-volume curves (0.00 m2, 54.82 10³m², 96.10 10³m², 142.89 10³m², etc). Finally, the third row holds 36 values of reservoir volume for the given values of elevation at the stage-area and stage-volume curves (0.00 10³m³, 35.78 10³m³, 111.14 10³m³, 231.23 10³m³, etc). The order of the sub-basins in the first column has to follow the same order of the sub-basin IDs as was used in ```hymo.dat``` (due to computational reasons); otherwise, an error message occurs. Sub-basins without outlet reservoirs or those without data on stage-area and stage-volume curves must not be entered in the file.

**6)** ```intake.dat``` <br>
(optional)

```
# Specification of controlled release through reservoir's intake devices in [m3/s]
Date,	doy,	10,	16
1012005	1	5	-999
2012005	2	5	0.1
…
```

*Date*: Date \[ddmmyyyy] (leading zeroes are optional)<br>
*Doy*:	day of year. Not used but kept for historical reasons, fill in dummy values. <br>
*Sub-basin ID*: Measured data on regulated outflow discharge through intake devices of the reservoir’s barrage in a specific sub-basin \[m³/s]

Example: This optional file allows specifying measured data on regulated outflow discharge of the sub-basin’s reservoir. If this file is not found in the folder Time_series, either a target value of regulated outflow discharge is given in the file ```reservoir.dat``` or a reservoir operation rule is provided in the file ```operat_rule.dat``` to the respective sub-basin’s reservoir. On January 1st 2005, a discharge of 5 m³/s was regulated through the intake devices of the barrage from the reservoir located at the outlet of the sub-basin 10. A regulated outflow discharge set to -999 for sub-basin 16 means that there is no measured data that day. In that case, the target value of regulated outflow discharge given in the file ```reservoir.dat``` (parameter damflow) is used. Time series of measured data on regulated outflow discharge of the sub-basin’s reservoir must be given for the whole simulation period. For those days without measurements, the value of regulated outflow discharge must be set to -999. Sub-basins without outlet reservoirs or those without measured data on regulated outflow discharge must not be entered.

**7)** ```hydraul_param.dat``` <br>
(optional)

```
# Specification of hydraulic parameters of the sub-basin’s reservoir
Subbasin-ID, nbr. cross sec, 1st row: manning [s/m**1/3], 2nd row: distance [-]
60	53	0.025	0.035	0.025	…	0.025	0.025	0.025 	in totalmax 200 IDs
60	53	209.485	199.605	162.748	…	260.775	237.29	138.492	in totalmax 200 IDs
```

*Subbasin-ID*: ID of sub-basin <br>
*nbr cross sec*: Number of cross sections in the sub-basin’s reservoir <br>
*1st row*: manning	1st row of each sub-basin: Manning's roughness for each cross section \[m<sup>-1/3</sup>/s] <br>
*2nd row*: distance	2nd row of each sub-basin: distance to the downstream cross section \[m]

Example: This optional file allows specifying hydraulic parameters for the calculation of water routing through the sub-basin’s reservoir. If this file is not found in the folder reservoir, a simplified modelling approach for the calculation of sediment balance is assumed. The reservoir located at the outlet point of the sub-basin with the ID 60 has 53 cross sections. The first row holds 53 values of Manning's roughness (0.025 m<sup>-1/3</sup>/s, 0.035 m<sup>-1/3</sup>/s, 0.025 m<sup>-1/3</sup>/s, etc). The second row holds 50 values of distance from a given cross section to the downstream cross section (209.485 m, 199.605 m, 162.748 m, etc). The order of the sub-basins in the first column has to follow the same order of the sub-basin IDs as was used in ```hymo.dat``` (due to computational reasons); otherwise an error message occurs. Sub-basins without outlet reservoirs or those without hydraulic data must not be entered in the file.

**8)** ```sed.dat```

```
# Specification of sedimentation parameters of the sub-basin’s reservoir
Subbasin-ID, dry_dens[ton/m**3], factor_actlay[-]
60	1.5	1
```

*Subbasin-ID*: ID of sub-basin <br>
*dry\_dens*: Dry bulk density of the sediment deposited in the sub-basin's reservoir \[ton/m<sup>3</sup>] <br>
*factor\_actlay*: Calibration parameter for the determination of the active layer thickness \[-]

Example: At the outlet point of the sub-basin with the ID 60 there is a reservoir with a dry bulk density of 1.5 ton/m<sup>3</sup>. The calibration parameter for the determination of the active layer thickness at that reservoir is equal to 1. It means that the default value of active layer thickness (set to 0.03 mm, derived from the simulation for the Barasona reservoir in Spain) is multiplied by a factor of 1. For the calculation of sediment balance using the simplified modelling approach, the third column with values of factor_actlay must not be entered in the file. Sub-basins without outlet reservoirs must not be entered in the file. Unknown subbasins are ignored. For unspecified reservoirs, default values of 1.5 and 1, respectively, will be used.

**9)** ```cross_sec_”ID”.dat``` <br> 
(optional)

```
# Specification of cross section geometry of the sub-basin’s reservoir
Subbasin-ID, section-ID, nbpoints, x-axis [m], y-axis[m]
60	1	8	81.18	460	119.29	455	…	319.32	460		in totalmax 200 IDs
60	2	12	60.72	460	189.24	450	…	382.93	460		in totalmax 200 IDs
…
```

*Subbasin-ID*: ID of sub-basin <br>
*Section-ID*:	ID of cross-section (continuously increasing, no gaps)<br>
*nbrpoints*: Number of points at the cross section of the sub-basin’s reservoir <br>
*x-axis*: Values at the x-axis for each point of the cross section in the sub-basin’s reservoir (from left to right, view from upstream side) \[m] <br>
*y-axis*: Values at the y-axis for each point of the cross section in the sub-basin’s reservoir (from left to right, view from upstream side) \[m]

Example: This optional file allows specifying detailed data on cross section geometry for the calculation of water routing through the sub-basin’s reservoir. If this file is not found in the folder reservoir, a simplified modelling approach for the calculation of sediment balance is assumed. The reservoir located at the outlet point of the sub-basin with the ID 60 was divided into 53 cross sections. The first row holds eight points with values at the x-axis (81.18 m, 119.29 m, etc) and y-axis (460 m, 450 m, etc) at the most upstream cross section of the of the sub-basin’s reservoir. The second row holds 12 points with values at the x-axis (60.72 m, 189.24 m, etc) and y-axis (460 m, 450 m, etc) at the next downstream cross section of the sub-basin’s reservoir. All cross sections of the sub-basin’s reservoir must be entered in the file. The value at the y-axis should be given after the value at the x-axis for a same point at the cross section. Sub-basins with data on cross section geometry must be entered in different input files (e.g. ```cross_sec_60.dat``` referred to sub-basin with ID 60). Sub-basins without outlet reservoirs or those without measured data on cross section geometry must not be entered.

**10)** ```original_sec_”ID”.dat``` <br> 
(optional)

```
# Specification of original cross section geometry of the sub-basin’s reservoir
Subbasin-ID, section-ID, nbpoints, y_original[m]
60	1	8	460	455	450	…	455	460			in totalmax 200 IDs
60	2	12	460	450	449	…	455	460			in totalmax 200 IDs
…
```

*Subbasin-ID*:	ID of sub-basin <br>
*Section-ID*:	ID of cross-section (continuously increasing, no gaps)<br>
*nbrpoints*: Number of points at the cross section of the sub-basin’s reservoir <br>
*y\_original*: Values of original bed elevation for each point of the cross section in the sub-basin’s reservoir (from left to right, view from upstream side) \[m]

Example: This optional file allows specifying detailed data on original cross section geometry for the calculation of water routing through the sub-basin’s reservoir. If this file is not found in the folder reservoir, there are two possibilities: sediment routing through the sub-basin’s reservoir is computed anyway, assuming that the original cross section geometry is the same as provided in the file ```cross_sec_”ID”.dat```; or a simplified modelling approach for the calculation of sediment balance is assumed (if the file ```cross_sec_”ID”.dat``` is not given). The reservoir located at the outlet point of the sub-basin with the ID 60 was divided into 53 cross sections. The first row holds eight values of original bed elevation for cross section 1 (460 m, 455 m, etc). The second row holds 12 values of original bed elevation for cross section 2 (460 m, 450 m, etc). All cross sections of the sub-basin’s reservoir must be entered in the file. Sub-basins with data on original cross section geometry must be entered in different input files (e.g. ```original_sec_60.dat``` referred to sub-basin with ID 60). Sub-basins without outlet reservoirs or those without measured data on original cross section geometry must not be entered.

**11)** ```sizedist_”ID”.dat``` <br> 
(optional)

```
# Specification of size distribution of original bed material along the cross sections of the sub-basin’s reservoir
Subbasin-ID, section-ID, frac_actlay[-]
60	1	0	0	0	0.8848	0.1101	0.0049	0.0002 	total number of sediment size classes
60	2	0.0307	0.0115	0.0012	0.7485	0.2044	0.0034	0.0003 	total number of sediment size classes
…
```

*Subbasin-ID*: ID of sub-basin <br>
*Section-ID*: ID of cross-section <br>
*y\_original*:	Values of sediment fraction for different size classes of the cross section in the sub-basin’s reservoir \[-]. The total number of sediment size classes is previously specified in the file ```do.dat```.

Example: This optional file allows specifying detailed data on size distribution of original bed material for the calculation of water routing through the sub-basin’s reservoir. If this file is not found in the folder reservoir, there are two possibilities: sediment routing through the sub-basin’s reservoir is computed anyway, assuming that no sediment was deposited the sub-basin’s reservoir previously; or a simplified modelling approach for the calculation of sediment balance is assumed (if the file ```cross_sec_”ID”.dat``` is not given). The reservoir located at the outlet point of the sub-basin with the ID 60 was divided into 53 cross sections. The first row holds values of sediment fraction for seven sediment size classes of cross section 1 (0, 0, etc). The second row holds values of sediment fraction for seven sediment size classes of cross section 2 (0.0307, 0.0115, etc). All cross sections of the sub-basin’s reservoir must be entered in the file. Sub-basins with data on size distribution of original bed material must be entered in different input files (e.g. ```sizedist_sec_60.dat``` referred to sub-basin with ID 60). Sub-basins without outlet reservoirs or those without measured data on size distribution of original bed material must not be entered.

**12)** ```main_channel.dat``` <br>
(optional)

```
# Specification of main channel geometry of the sub-basin’s reservoir
Subbasin-ID, nbr. cross sec, 1st row: pt1 [-], 2nd row: pt2 [-]
60	53	8	10	15	…	40	65	70 	in totalmax 200 IDs
60	53	15	17	25	…	60	80	90	in totalmax 200 IDs
```

*Subbasin-ID*: ID of sub-basin <br>
*nbr cross sec*: Number of cross sections in the sub-basin’s reservoir <br>
*1st row*: manning	1st row of each sub-basin: Point of the cross section in the sub-basin’s reservoir that identifies the beginning of main channel (from left to right, view from upstream side) \[-] <br>
*2nd row*: distance	2nd row of each sub-basin: Point of the cross section in the sub-basin’s reservoir that identifies the end of main channel (from left to right, view from upstream side) \[-]

Example: This optional file allows specifying the exact location of main channel in the cross sections of the sub-basin’s reservoir. This information is used to adjust bed profiles of cross sections, avoiding steeper slopes caused by erosion processes. If this file is not found in the folder reservoir, there are two possibilities: sediment routing through the sub-basin’s reservoir is computed anyway, disregarding the occurrence of steeper slopes; or a simplified modelling approach for the calculation of sediment balance is assumed (if the file ```cross_sec_”ID”.dat``` is not given). The reservoir located at the outlet point of the sub-basin with the ID 60 has 53 cross sections. The main channel of cross section 1 is located between the 8th and 15th points (cross section 2: located between the 10th and 17th points; etc). A value of -999 indicates unknown location of main channel for that cross section. Sub-basins without outlet reservoirs or those without data on location of main channel of cross sections must not be entered in the file.

**13)** ```lake.dat```

```
# Specification of parameters for the reservoir size classes
Reservoir_class-ID, maxlake0[m**3], lake_vol0_factor[-], lake_change[-], alpha_Molle[-], damk_Molle[-], damc_hrr[-], damd_hrr[-]
1	5000	0.2	0.10	2.7	1500	7	1.5
2	25000	0.2	0	2.7	1500	14	1.5
3	50000	0.2	0	2.7	1500	21	1.5
4	100000	0.2	0	2.7	1500	28	1.5
5	250000	0.2	0	2.7	1500	35	1.5
```

*Reservoir\_class-ID*: ID of reservoir size class <br>
*maxlake0*: Upper limit of reservoir size class in terms of volume \[m³] <br>
*lake\_vol0\_factor*: Fraction of storage capacity that indicates the initial water volume in the reservoir size classes \[-] <br>
*lake\_change*: Factor that indicates yearly variation in the number of reservoirs of the size classes \[-] <br>
*alpha\_Molle, damk_Molle*: Parameters of the area-volume relationship in the reservoir size classes (Area=alpha.k.(Vol/k)<sup>alpha/(alpha-1)</sup>) \[-]. Values of reservoir area and volume are expressed in m² and m³, respectively <br>
*damc\_hrr, damd\_hrr*: Parameters of the spillway rating curve in the reservoir size classes (Qout=damc\_hrr.Hv<sup>damd\_hrr</sup>) \[-]. Values of water height over the spillway and overflow discharges are expressed in m and m³/s, respectively

Example: The study area has a network of small reservoirs, which are grouped into five size classes according to their storage capacity (changes on the number of size classes are not available yet). The water and sediment balances of small reservoirs are computed for one hypothetical representative reservoir of mean characteristics. The size class 1 has reservoirs with storage capacity up to 5,000 m³, an initial water volume of 20% of the storage capacity, an yearly increase of 10% in the number of reservoirs for the simulation period (lake_increase set to zero means that the number of reservoirs remains constant), an area-volume relationship with parameters alpha_Molle and damk_Molle set to 2.7 and 1500, respectively, and a spillway rating curve with parameters damc_hrr and damd_hrr set to 7 and 1.5, respectively.

**14)** ```lake_maxvol.dat``` <br> 
(optional)

```
# Specification of water storage capacity for the reservoir size classes
Sub-basin-ID, maxlake[m**3] (five reservoir size classes)
60	2627.21	16591.52	0	0	0
…
```

*Subbasin-ID*: ID of sub-basin <br>
*maxlake*: Mean value of initial storage capacity of the hypothetical representative reservoirs of the size classes \[m³]. Value varies because of the sediment accumulation

Example: This optional file allows specifying data on initial storage capacity of the hypothetical representative reservoirs of the size classes. If this file is not found in the folder reservoir, the initial storage capacity is computed as a geometrical mean of the lower and upper limit of the reservoir size classes (except for class 1, computed as 50% of its upper limit). The sub-basin with the ID 60 has only two reservoir size classes with initial storage capacities of 2627.21 m³ and 16591.52 m³ (size classes 1 and 2, respectively). Therefore, there is no reservoir at the size classes 3 to 5 for that sub-basin. The order of the sub-basins in the first column has to follow the same order of the sub-basin IDs as was used in ```hymo.dat``` (due to computational reasons); otherwise, an error message occurs. Sub-basins without networks of small reservoirs must not be entered in the file.

**15)** ```lake_year.dat``` <br> 
(optional)

```
# Specification of changes on the number of reservoirs in the size classes
Year, sub-basin-ID, acudfloatyear[-] (five reservoir size classes)
2005	15	32	34	17	20	11
2005	60	111	72	22	25	59
2006	15	34	36	18	21	12
2006	60	117	76	23	26	62
2007	15	35	37	19	22	12
2007	60	122	79	24	28	65
2008	15	37	39	20	23	13
2008	60	128	83	25	29	68
```

*Year*: Year of simulation \[yyyy] <br>
*Subbasin-ID*:		ID of sub-basin
*acudfloatyear*:	Total number of reservoirs in the sub-basin and size classes for all years of simulation \[m³]

Example: This optional file allows specifying changes on the number of reservoirs in the size classes. If this file is not found in the folder reservoir, a yearly variation in the number of reservoirs for sub-basins and size classes given in the file ```lake.dat``` is assumed. In the year 2005, the sub-basin with the ID 15 has 32 reservoirs of class 1, 34 of class 2, 17 of class 3, 20 of class 4, and 11 of class 5. The order of the sub-basins in the second column has to follow the same order of the sub-basin IDs for all years of simulation as was used in ```hymo.dat``` (due to computational reasons); otherwise, an error message occurs. Sub-basins without networks of small reservoirs must not be entered in the file.

**16)** ```lake_number.dat```

```
# Specification of total number of reservoirs in the size classes
Sub-basin-ID, acud[-] (five reservoir size classes)
60	15	8	0	0	0
```

*Subbasin-ID*: ID of sub-basin <br>
*acud*: Total number of reservoirs in the size classes \[-]

Example: The sub-basin with the ID 60 has 15 and 8 reservoirs of the size classes 1 and 2, respectively. Therefore, there is no reservoir of size classes 3 to 5 for that sub-basin. The order of the sub-basins in the first column has to follow the same order of the sub-basin IDs as was used in ```hymo.dat``` (due to computational reasons); otherwise, an error message occurs. Sub-basins without networks of small reservoirs must not be entered in the file.

**17)** ```lake_frarea.dat``` <br>
(optional)

```
# Specification of runoff contributing area for the reservoir size classes
Sub-basin-ID, lakefrarea[-] (five reservoir size classes)
60	0.240	0.250	0	0	0
```

*Subbasin-ID*: ID of sub-basin <br>
*maxlake*: Fraction of sub-basin area that represents the runoff contributing area for the reservoir size classes \[-]

Example: This optional file allows specifying data on runoff contributing area for the reservoir size classes. If this file is not found in the folder reservoir, the runoff contributing area is equally divided into the five reservoir size classes (one-sixth to each class). Another sixth part is attributed to the area not-controlled by the reservoir network. The sub-basin with the ID 60 has only two reservoir size classes with a runoff contributing area covering 24% and 25% of the sub-basin area (size classes 1 and 2, respectively). Therefore, there is no reservoir of size classes 3 to 5 for that sub-basin. The order of the sub-basins in the first column has to follow the same order of the sub-basin IDs as was used in ```hymo.dat``` (due to computational reasons); otherwise, an error message occurs. Sub-basins without networks of small reservoirs must not be entered in the file.

<a name="3-5-input-of-climate-data"></a>
### 3.5 Input of climate data

The WASA model requires time series for precipitation (daily or hourly), short wave radiation, humidity and temperature (daily). The input files are located in the folder ```Input\[case_study]\Time_series``` and are summarised below.


**1)** ```rain_daily.dat``` (only needed when run in daily resolution)

```
# Daily average precipitation [mm/d] for each subbasin, ordered according to IDs				
Date,	No. of days, Subbasin-ID, Subbasin-ID, Subbasin-ID, Subbasin-ID, Subbasin-ID, Subbasin-ID, Subbasin-ID, Subbasin-ID, Subbasin-ID
0	0	49	50	1	44	10	4	15	39	3	29
01011980	1	40	40	40	40	40	40	40	40	40	40
02011980	2	40	40	40	40	40	40	40	40	40	40
03011980	3	40	40	40	40	40	40	40	40	40	40
…
```

**2)** ```rain_hourly.dat``` (only needed when run in hourly resolution)

```
hourly rainfall				
Date	No. of timestep	Subbasin-ID.		
0	0	1	2	3
13092006	20	0	0	0
13092006	21	0	0	0
13092006	32	0	0	0
13092006	23	0.2	0.3	0.1
14092006	0	0	0.1	0
14092006	1	0.4	0.2	0
…
```

**3)** ```temperature.dat```

```
# Daily average temperature (in degree Celcius) for each subbasin, ordered according to IDs
Date,	No. of days, Subbasin-ID, Subbasin-ID, Subbasin-ID, Subbasin-ID, Subbasin-ID, Subbasin-ID, Subbasin-ID, Subbasin-ID, Subbasin-ID
0	0	49	50	1	44	10	4	15	39	3	29
01011980	1	15	15	15	15	15	15	15	15	15	15
02011980	2	15	15	15	15	15	15	15	15	15	15
03011980	3	15	15	15	15	15	15	15	15	15	15
…
```

**4)** ```humidity.dat```

```
# Daily average humidity [in %] for each subbasin, ordered according to IDs
Date,	No. of days, Subbasin-ID, Subbasin-ID, Subbasin-ID, Subbasin-ID, Subbasin-ID, Subbasin-ID, Subbasin-ID, Subbasin-ID, Subbasin-ID
0	0	49	50	1	44	10	4	15	39	3	29
01011980	1	75	75	75	75	75	75	75	75	75	75
02011980	2	75	75	75	75	75	75	75	75	75	75
03011980	3	75	75	75	75	75	75	75	75	75	75
…
```

**5)** ```radiation.dat```

```
# Daily average shortwave radiation [in W/m2] for each subbasin, ordered according to IDs
Date,	No. of days, Subbasin-ID, Subbasin-ID, Subbasin-ID, Subbasin-ID, Subbasin-ID, Subbasin-ID, Subbasin-ID, Subbasin-ID, Subbasin-ID
0	0	49	50	1	44	10	4	15	39	3	29
01011980	1	260	260	260	260	260	260	260	260	260	260
02011980	2	260	260	260	260	260	260	260	260	260	260
03011980	3	260	260	260	260	260	260	260	260	260	260
…
```

*Date*: Continuous number of day month year <br>
*No. of days*: Continuous numbering <br>
*Subbasin-ID*:  ID of sub-basin, same ordering at ```hymo.dat```

Example: The four files are organised in the same manner. Here they are given for three days: 01.01.1980 until 03.01.1980. In the examples above, the time series are uniform for each sub-basin, however, it is possible to assign different time series to individual sub-basins. 

**6)** ```extraterrestrial_radiation.dat```

```
Extra-terrestrial shortwave radiation as monthly mean daily value in [W/m2]
441
447
…
```
This file specifies the extraterrestrial incoming shortwave radiation at the top of the atmosphere \[W/m2]. The values are daily averages for each month from January until December (12 values). 

<a name="3-6-state-variables"></a>
### 3.6 Input/Output of state variables (optional)
The state of (some) system variables can be stored in files. This allows their later analysis and resuming runs from them. Thus, any of these files can be generated by the model, and used as an input in a successive model run. These capabilities can be enabled in ```do.dat``` lines 36-37 (see [Figure 1](#figure-1)).
These files are optional. If present and enabled in ```do.dat```, WASA-SED loads them at model startup. Each file is first searched for in the ```Output``` directory. If not found there, WASA-SED searches in ```Input/init_conds```. If not present in either, default values are used (e.g., 100 % relative saturation for the soils, 0 for most other entities). After startup, WASA-SED will generate files with the extension ```stat_start```, that describe the status of the model in the beginning of the simulation (these should be identical to the ```stat``` when provided). New ```stat``` files are written at the the end of each simulation year.


<a name="table-12"></a>
**Table 12:** Input/Output of state variables.

Parameter File |	Content
---|---
```gw_storage.stat``` | Initialisation of storage content of ground water
```intercept_storage.stat``` (optional) | Initialisation of storage content of interception
```soil_moisture.stat``` | Initialisation of storage content of soil moisture
```gw_storage.stat``` | Ground water storage
```interflow_storage.stat``` | state of intra-hillslope subsurface flow
```snow_storage.stat``` | Snow storage
```river_storage.stat``` | river storages
```lake_storage.stat``` | small reservoirs
```reservoir_storage.stat``` | strategic reservoirs
```sediment_storage.stat'    ``` | sediment storage in riverbed (deposited)
```susp_sediment_storage.stat``` | sediment storage in riverbed (suspended)
```storage.stats``` | Summary of storages 

gw_storage.stat:
```
# Ground water storage (for analysis or model re-start)
Sub-basin,	LU,	volume_[mm],	area_[m²]
1 11111	39.51	1013437.4
…
```

intercept_storage.stat:
```
# Interception storage (for analysis or model re-start)
Sub-basin,	LU,	TC,	SVC,	storage_[mm],	area_[m²]
1	11111	5	16010	    0.00	     73491.9
…
```
soil_moisture.stat:
```
# Soil moisture status (for analysis or model re-start)
Sub-basin,	LU,	TC,	SVC,	horizon,	watercontent_[mm],	area_[m²]
1	11111	5	16010	1	   26.70	     73491.9
…
```

interflow_storage.stat:
```
# Interflow storage (for analysis or model re-start)
 Sub-basin,	LU,	TC,	horizon,	storage_[m3]
1	1	3	1	   0.000
...
```

snow_storage.stat:
```
#Snow storage (for analysis or model re-start)
Sub-basin,	LU,	horizon,	storage [m],	energy [kJ/m²],	albedo [-]
1	2	10	0.000	0.000	0.880
...
```

river_storage.stat:
The structure of ```river_storage.stat``` depends on the routing option chosen (UHG or Muskingum, see examples below).
```
# UHG routing: values of routed discharge per timestep of unit hydrograph
Sub-basin	[n_h x timestep]
101	1.286	0.649	0.074
...
```

```
# Muskingum routing: river reach volume status (for analysis or model re-start)
Sub-basin	volume[m^3]
1	2033.808
...
```

lake_storage.stat:
```
# Lake volume status (for analysis or model re-start)
Sub-basin	reservoir_size_class,	volume[m^3]
1	1	       0.00
...
```

reservoir_storage.stat:
```
Reservoir volume status (for analysis or model re-start)
Subbasin	volume[m^3]
221	 25500.000
231	 30000.000
```


\[[Table of contents](#toc)]
## 4 Output data

The location of the output folder is specified in the ```do.dat```. By default, the output folder is set to ```WASA\Output```. The parameter file ```parameter.out``` echoes the main parameter specification for the WASA model, as were given in the ```do.dat``` file.

<a name="4-1-output-of-the-hillslope-module"></a>
### 4.1 Output of the hillslope module

The output files generated by the hillslope routine are shown in [Table 13](#table-13).

<a name="table-13"></a>
**Table 13:** Output files of the hillslope module.

Output file | Content
---|---
```daily_actetranspiration.out``` | daily actual evapotranspiration \[mm/d] for all sub-basins (MAP-IDs) incl. river
```daily_potetranspiration.out```	| daily potential evapotranspiration \[mm/d] for all sub-basins 
```daily_qhorton.out``` | daily horton overland flow \[m<sup>3</sup>] for all sub-basins 
```daily_sediment_production.out``` | daily sediment production \[t] for all sub-basins and particle classes (i.e. sediment input of hillslopes into the river)
```daily_subsurface_runoff.out```	| daily total subsurface runoff \[m<sup>3</sup>/d] for all sub-basins 
```daily_theta.out```	| mean soil moisture in profile \[mm] for all sub-basins 
```daily_total_overlandflow.out``` | total overland flow \[m<sup>3</sup>] for all sub-basins 
```daily_water_sub-basin.out``` | daily water contribution into river \[m<sup>3</sup>/s] for all sub-basins 
```water_sub-basin.out``` | sub-daily contribution to river \[m<sup>3</sup>/s] for all sub-basins 
```sediment_production.out``` |	daily sediment production \[t] for all sub-basins  and particle classes
```Daily_gw_loss.out``` | daily water loss from model domain due to deep seepage in LU without GW
```deep_gw_discharge.out``` | total deep ground water discharge \[m<sup>3</sup>] for all sub-basins 
```deep_gw_recharge.out``` | total deep ground water recharge \[m<sup>3</sup>] for all sub-basins 
```actetranspiration.out``` |	Subdaily actual evapotranspiration \[mm] for all sub-basins  incl. river
```qhorton.out``` | subdaily horton overland flow \[m<sup>3</sup>] for all sub-basins 
```subsurface_runoff.out``` | subdaily total subsurface runoff \[m<sup>3</sup>/d] for all sub-basins 
```total_overlandflow.out``` | Subdaily total overland flow \[m<sup>3</sup>] for all sub-basins 
```gw_discharge.out``` | groundwater discharge \[m<sup>3</sup>/timestep] for all sub-basins 
```potetranspiration.out``` | Subdaily potential evapotranspiration \[mm/d] for all sub-basins 
```gw_loss.out``` | model loss (deep seepage) \[m<sup>3</sup>/timestep] for all sub-basins  
```gw_recharge.out``` | groundwater recharge \[m<sup>3</sup>/timestep] for all sub-basins 
```*.stat```| various files containing state of system variables at end of simulation, see ...
```snow*``` | optional output files from the snow routine, see [Table 5](#table-5)
```gw_storage.stat, intercept_storage.stat, soil_moisture.stat, snow_storage.stat``` (optional) | Initialisation state variables, see [section 3.6](#3-6-state-variables)

The output files ```daily_water_sub-basin.out```, ```sediment_production.out``` and ```water_sub-basin.out``` include the effect of the *distributed* reservoirs. All other remaining output files above contain the raw output of the hillslope module (no reservoir effects). All above-mentioned files have the same structure, as shown by the example ```daily_actetranspiration.out``` below (the subdaily output files additionally contain the timestep number in the third column):

```
# Daily actual evapotranspiration [mm/d]  for all sub-basins (MAP-IDs)
Year,    Day,      57,           15,           20,           60
1980     1         3.495         2.974         3.258         3.412
1980     2         3.504         3.076         3.424         3.436
    …
```
    
*Year*: year of simulation <br>
*Day*: day of current year of simulation <br>
*\[variable]*: respective variable for each sub-basin

Beware: “day” counts the number of days in the respective simulation year, i.e. if you start your simulation on May, 1, the number “1” refers to this day and the rest of the days in that year will be counted till 306.

```gw_storage.stat```, ```intercept_storage.stat```, ```soil_moisture.stat``` and ```storage.stats```:<br>
These files are written at the end of each simulation year, thus allowing recommencing an aborted WASA run starting from the last simulation timestep. <br>
Beware: all other output files are overwritten in this case. For file structure, see section [Input data](#input-data). ```storage.stats``` contains the overall summary of storages corresponding to the three files mentioned before.

<a name="4-2-output-of-the-river-module"></a>
### 4.2 Output of the river module

The river routine calculates the water and sediment discharge in each river stretch. The following files (see [Table 14](#table-14)) are generated as time series, if enabled, depending on the selected routing scheme.

<a name="table-14"></a>
**Table 14:** Output files of the river module.

Output file | Content
---|---
```River_flow.out``` | River discharge in m<sup>3</sup>/s
```River_storage.out``` | River storage volume in m<sup>3</sup>
```River_velocity.out``` | Flow velocity in m/s
```River_flowdepth.out```	| Flow depth in m
```River_Flow_dailyaverage.out``` | Daily averaged  flow in m<sup>3</sup>/s
```River_Sediment_total.out``` | Suspended sediment in tons/timestep
```River_Sediment_Concentration.out``` | Suspended sediment concentration in g/l
```River_Sediment_total_dailyaverage.out``` | Daily averaged sediment flux in tons/h
```River_Degradation.out``` | Degradation of sediment in riverbed in tons/stretch
```River_Deposition.out``` | Deposition of sediment in riverbed in tons/stretch
```River_Bedload.out``` | Bedload rate for 5 formulas in kg/s
```Routing_response.out``` | Linear response function for routing scheme 1
```River_Sediment_storage.out``` | Deposited sediment stored in river reach in t
```River_Susp_Sediment_storage.out``` | Suspended sediment stored in river reach in t
```River_Infiltration.out``` | Infiltration of river stretches
```river_storage.stat, sediment_storage.stat, susp_sediment_storage.stat``` (optional) | Initialisation state variables, see [section 3.6](#3-6-state-variables)


All above-mentioned files have the same structure, as shown by the example ```River_flow.out``` below:

```
# Output files for river discharge q_out (m3/s) (with MAP IDs as in hymo.dat)
Year  Day    dt   9          10             11
2009     1     1         6.313         1.797         8.922
2009     1     2         6.176         1.744         8.733
2009     1     3         4.001         1.029         5.878
    …
```
    
*Subbasin-ID*: ID of all sub-basins in the second line of the file <br>
*Timestep*: Timestep as specified in the ```do.dat``` in \[hours] <br>
*Time series*: water discharge in river stretch in m<sup>3</sup>/s

Example: After each time step (i.e. hour or day), the discharge is given for each sub-basin, e.g. Sub-basin No. 9 has a discharge of 6.313 m<sup>3</sup>/s, Sub-basin No. 10 of 1.797 m<sup>3</sup>/s and Sub-basin No. 11 of 8.922 m<sup>3</sup>/s after 1 hours.

<a name="4-3-output-of-the-reservoir-module"></a>
### 4.3 Output of the reservoir module

The reservoir module simulates the water and sediment transport through the reservoirs located in the study area. Currently, the output comprises results on water balance, hydraulic calculations, sediment transport and bed elevation changes for all reservoirs located at the outlet point of the sub-basins. The results are printed for all outlet reservoirs separately, identified by the ID of the sub-basin where it is located. Additional files are also printed for the reservoir size classes. [Table 15](#table-15) shows the generated files.

<a name="table-15"></a>
**Table 15:** Output files of the reservoir module.

Nr. |Output file | Content
---|---|---
1 | ```res_”ID”_watbal.out``` | Water balance components of outlet reservoirs 
2 | ```res_”ID”_vollost.out``` | Dead volume, alert volume, and storage capacity of outlet reservoirs 
3 | ```res_”ID”_cav.out``` | Stage-area and stage-volume curves of outlet reservoirs 
4 | ```res_”ID”_hydraul.out``` | Hydraulic components of outlet reservoirs
5 | ```res_”ID”_sec”ID”_bedchange.out``` | Bed elevation at cross sections (identified by a specific Section-ID) of outlet reservoirs
6 | ```res_”ID”_sedbal.out``` | Sediment balance components of outlet reservoirs
7 | ```res_”ID”_longitudinal.out``` | Longitudinal bed profile of outlet reservoirs
8 | ```res_”ID”_sedcomposition.out``` | Effluent grain size distribution of outlet reservoirs
9 | ```lake_inflow_r.out``` | Water inflow into the reservoir size classes<sup>1</sup>
10 | ```lake_outflow_r.out``` | Water outflow from the reservoir size classes<sup>1</sup>
11 | ```lake_retention_r.out``` | Water retention in the reservoir size classes<sup>1</sup>
12 | ```lake_storage_r.out``` | Water volume of the reservoir size classes<sup>1</sup>
13 | ```lake_sedinflow_r.out``` | Sediment inflow into the reservoir size classes<sup>1</sup>
14 | ```lake_sedoutflow_r.out``` | Sediment outflow from the reservoir size classes<sup>1</sup>
15 | ```lake_sedretention_r.out``` | Sediment retention in the reservoir size classes<sup>1</sup>
16 | ```lake_sedimentation_r.out``` | Cumulative sediment deposition in the reservoir size classes<sup>1</sup>
17 | ```lake_watbal.out``` | Water balance components of all upstream reservoirs<sup>2</sup>
18 | ```lake_sedbal.out``` | Sediment balance components of all upstream reservoirs<sup>2</sup>
19 | ```lake_inflow.out``` | Water inflow into the reservoir size classes<sup>3</sup>
20 | ```lake_outflow.out``` | Water outflow from the reservoir size classes<sup>3</sup>
22 | ```lake_retention.out``` | Water retention in the reservoir size classes<sup>3</sup>
23 | ```lake_vollost.out``` | Sediment retention in the reservoir size classes per timestep<sup>3</sup>
24 | ```lake_maxstorcap.out``` | Current (remaining) maximum sotrage capacity for reservoir size classes<sup>3</sup>
25 | ```lake_volume.out``` | Water stored in the reservoir size classes<sup>3</sup>
26 | ```lake_volume_r.out``` | Total water stored in the reservoir classes<sup>1</sup>
27 | ```lake_sedinflow.out``` | Sediment inflow into the reservoir size classes<sup>3</sup>
28 | ```lake_sedoutflow.out``` | Sediment outflow from the reservoir size classes<sup>3</sup>
29 | ```lake_sizedistoutflow.out``` | Effluent grain size distribution of the reservoir size classes <sup>4</sup>
30 | ```lake_storage.stat, reservoir_storage.stat``` (optional) | Initialisation state variables, see [section 3.6](#3-6-state-variables)

<sup>1</sup> Results are aggregated over all sub-basins grouped by reservoir size class (one value for the whole catchment and each reservoir size class) <br>
<sup>2</sup> Results are displayed for the whole catchment without distinguishing between size classes (one value for the whole catchment) <br>
<sup>3</sup> Results are displayed for each sub-basin and reservoir size class (one value for each sub-basin and reservoir size class). <br>
<sup>4</sup> Results are displayed for all sub-basins without distinguishing between size classes (one value for each sub-basin).

**1)** ```res_”ID”_watbal.out```

```
Subasin-ID	year	day	hour	qlateral(m**3/s)	inflow(m**3/s)	evap(m**3)	prec(m**3)	intake(m**3/s)	overflow(m**3/s)	qbottom(m**3/s)	qout(m**3/s)	withdrawal(m**3/s)	    elevation(m)	area(m**2)	volume(m**3
60	1980	1	1	55.04	6.12	25478	0	0		0.00	0.00	6.12	440.86	5255332.50	49625572.00
60	1980	1	2	42.01	6.12	24579	0	0		0.00	0.00	6.12	441.48	4580464.00	52922032.00
…
```

*Subbasin-ID*: ID of sub-basin <br>
*year*: Year of simulation <br>
*day*: Day of simulation <br>
*hour*: Hour of simulation <br>
*qlateral*: lateral inflow discharge into the subbasin's reservoir \[m<sup>3</sup>/s] <br>
*inflow*: Water inflow discharges into the sub-basin's reservoir \[m<sup>3</sup>/s] <br>
*evap*: evaporation from the subbasin's reservoir \[m<sup>3</sup>] <br>
*prec*: precipitation into the subbasin's reservoir \[m<sup>3</sup>] <br>
*intake*: Water outflow discharges through water intake devices in the sub-basin's reservoir \[m<sup>3</sup>/s]. Should correspond to values in  ```intake.dat```.<br>
*overflow*: Water overflow discharges in the sub-basin's reservoir \[m<sup>3</sup>/s] <br>
*qbottom*: Water outflow discharges through bottom outlets in the sub-basin's reservoir \[m<sup>3</sup>/s] <br>
*qout*: Total outflow discharges in the sub-basin's reservoir \[m<sup>3</sup>/s] <br>
*withdrawal*: withdrawal water volume to supply the water use sectors Should correspond to values in  ```reservoir.dat```.) \[m<sup>3</sup>/s] <br>
*elevation*: Reservoir level in the sub-basin's reservoir \[m] <br>
*area*: Reservoir area in the sub-basin's reservoir \[m<sup>2</sup>] <br>
*volume*: Reservoir volume in the sub-basin's reservoir \[m<sup>3</sup>]

Example: After each time step, e.g. after one day, the reservoir of the sub-basin with the ID 60 has a water inflow discharge of 55.04 m<sup>3</sup>/s, a water outflow discharge through the intake device of 6.12 m<sup>3</sup>/s, no water overflow discharge, no water outflow discharge through the bottom outlets, a total water outflow discharge of 6.12 m<sup>3</sup>/s, a water level of 440.86, a reservoir area of 5,255,332.50 m<sup>2</sup> and a reservoir volume of 49,625,572.00 m<sup>3</sup>. Currently, the model generates an output file for each reservoir considered in the simulation (e.g. ```res_60_watbal.out``` referred to sub-basin with ID 60)

**2)** ```res_”ID”_vollost.out```

```
Subbasin-ID, year, day, hour, deadvol(m**3), alertvol(m**3), storcap(m**3)
60  1980   1   1     4795484.24    45171678.11    91744848.62
60  1980   1   2     4795457.23    45171322.04    91744690.30
…
```

*Subbasin-ID*: ID of sub-basin <br>
*year*: Year of simulation <br>
*day*: Day of simulation <br>
*hour*: Hour of simulation <br>
*deadvol*: Dead volume in the sub-basin's reservoir \[m<sup>3</sup>] <br>
*alertvol*: Alert volume in the sub-basin's reservoir \[m<sup>3</sup>] <br>
*storvap*: Storage capacity in the sub-basin's reservoir \[m<sup>3</sup>]

Example: After each time step, e.g. after one day, the reservoir of the sub-basin with the ID 60 has a dead volume of 4,795,484.24 m<sup>3</sup>, an alert volume of 45,171,678.11 m<sup>3</sup>, and an alert volume of 91,744,848.62 m<sup>3</sup>. Currently, the model generates an output file for each reservoir considered in the simulation (e.g. ```res_60_watbal.out``` refers to sub-basin with ID 60).

**3)** ```res_”ID”_cav.out```

```
Subbasin-ID, year, day, hour, 1st row: elev_bat(m), 2nd row: area_bat(m**2), 3rd row: vol_bat(m**3)
60  1980   1   1         413.34        415.00         416.00	…	        447.00          447.67          448.00
60  1980   1   1           0.00      79176.34      122767.10	…	    5872791.16      6020551.37      7288430.00
60  1980   1   1           0.00      35021.39      110657.75	…   	86999642.73     91744848.62     94132893.50
…
```

*Subbasin-ID*: ID of sub-basin <br>
*year*: Year of simulation <br>
*day*: Day of simulation <br>
*hour*: Hour of simulation <br>
*1st row*: elevation	1st row: bed elevation values of the stage-area and stage-volume curves after sediment erosion/deposition in the sub-basin's reservoir \[m] <br>
*2nd row*: res. area	2nd row: reservoir area values of the stage-area and stage-volume curves after sediment erosion/deposition in the sub-basin's reservoir \[m<sup>2</sup>] <br>
*3rd row*: res. volume	3rd row: reservoir volume of the stage-area and stage-volume curves after sediment erosion/deposition in the sub-basin's reservoir \[m<sup>3</sup>]

Example: After each time step, e.g. after one day, the reservoir of the sub-basin with the ID 60 has 36 new points at the stage-area and stage-volume curves changed due to sediment erosion/deposition. The first row holds 36 values of water elevation at the stage-area-volume curve (413.34 m, 415.00 m, 416.00 m, etc). The second row holds 36 values of corresponding reservoir area (0.00 m<sup>2</sup>, 79,176.34 m<sup>2</sup>, 122,767.10 m2, etc). Finally, the third row holds 36 values of corresponding reservoir volume (0.00 m<sup>3</sup>, 35,021.39 m<sup>3</sup>, 110,657.75 m<sup>3</sup>, etc). Currently, the model generates an output file for each reservoir considered in the simulation (e.g. ```res_60_cav.out``` referred to sub-basin with ID 60).

**4)** ```res_”ID”_hydraul.out```

```
Subbasin-ID, year, day, hour, section-ID, depth_sec(m), watelev_sec(m), area_sec(m**2), topwidth_sec(m),
 energslope_sec(-), hydrad_sec(m), meanvel_sec(m/s), discharge_sec(m**3/s)
60  1980   1   1    1          1.325        448.635         34.717         40.183      0.192E-02       0.858863       1.585346      55.038498
60  1980   1   1    2          1.376        446.956         33.455         56.513      0.345E-02       0.585371       1.645158      55.038498
…
```

*Subbasin-ID*: ID of sub-basin <br>
*year*: Year of simulation <br>
*day*: Day of simulation <br>
*hour*: Hour of simulation <br>
*section-ID*: ID of cross-section in the sub-basin's reservoir <br>
*depth\_sec*: Water depth of each cross section in the sub-basin's reservoir \[m] <br>
*watelev\_sec*: Water elevation of each cross section in the sub-basin's reservoir \[m] <br>
*area\_sec*: Wetted area of each cross section in the sub-basin's reservoir \[m<sup>2</sup>] <br>
*topwidth\_sec*: Top width of each cross section in the sub-basin's reservoir \[m] <br>
*energslope\_sec*: Slope of energy-grade line of each cross section in the sub-basin's reservoir \[-] <br>
*hydrad\_sec*: Hydraulic radius of each cross section in the sub-basin's reservoir \[m] <br>
*meanvel\_sec*: Mean velocity of each cross section in the sub-basin's reservoir \[m/s] <br>
*discharge\_sec*: Discharge of each cross section in the sub-basin's reservoir \[m<sup>3</sup>/s]

Example: After each time step, e.g. after one day, the most upstream cross section (section 1) of the reservoir of the sub-basin with the ID 60 has a water depth of 1.325 m, a water elevation of 448.635 m, a wetted area of 34.717 m2, a top width of 40.183 m, a slope of energy-grade line of 0.00192, a hydraulic radius of 0.858863 m, a mean velocity of 1.585346 m/s and a discharge of 55.038498 m<sup>3</sup>/s. Currently, the model generates an output file for each reservoir considered in the simulation (e.g. ```res_60_hydraul.out``` referred to sub-basin with ID 60).

**5)** ```res_”ID”_sec”ID”_bedchange.out```

```
Subbasin-ID, section-ID, year, day, hour, nbr. points, y-axis(m)
60   1   1980   1   1   11     460.000000     451.000000     450.000000	…	450.000000     451.000000     460.000000
60   1   1980   1   2   11     460.000000     451.000000     450.000000	…	450.000000     451.000000     460.000000
…
```

*Subbasin-ID*: ID of sub-basin <br>
*section-ID*: ID of cross-section in the sub-basin's reservoir <br>
*year*: Year of simulation <br>
*day*: Day of simulation <br>
*hour*: Hour of simulation <br>
*nbr. points*: Number of points at the cross section in the sub-basin's reservoir <br>
*y-axis*: Bed elevation changes in the cross section of the reservoir (from left to right, seen from upstream) \[m]

Example: After each time step, e.g. after one day, the most upstream cross section (section 1) of the reservoir of the sub-basin with the ID 60 holds 11 values of bed elevation, changed because of either deposition or erosion processes (460.00 m, 451.00 m, 450.00 m, etc). Currently, the model generates an output file for each cross section of the sub-basin’s reservoir (e.g. ```res_60_sec1_bedchange.out``` referred to cross section 1 and sub-basin with ID 60).

**6)** ```res_”ID”_sedbal.out```

```
Subbasin-ID, year, day, hour, sed_input(ton/timestep), sed_output(ton/timestep), sedimentation(ton/timestep), cum_sedimentation(ton)
60  1980   1   1     78555.180     1799.437   76755.743      76755.743
60  1980   1   2        240.464         10.663        229.801      76985.544
…
```

*Subbasin-ID*: ID of sub-basin <br>
*year*: Year of simulation <br>
*day*: Day of simulation <br>
*hour*: Hour of simulation <br>
*sed\_input*: Sediment inflow discharges into the sub-basin's reservoir \[ton/timestep] <br>
*sed\_output*: Sediment outflow discharges out the sub-basin's reservoir \[ton/timestep] <br>
*sedimentation*: Sediment deposition rate in the sub-basin's reservoir \[ton/timestep] <br>
*cum\_sedimentation*: Cumulative sediment deposition in the sub-basin's reservoir since dam construction \[ton]

Example: After each time step, e.g. after one day, the reservoir of the sub-basin with the ID 60 has a sediment inflow discharge of 78,555.180 ton/timestep, a sediment outflow discharge of 1,799.437 ton/timestep, a sediment deposition rate of 76,755.743 ton/timestep and a cumulative sediment deposition of 76,755.743 ton/timestep. Currently, the model generates an output file for each reservoir considered in the simulation (e.g. ```res_60 _sedbal.out``` referred to sub-basin with ID 60).

**7)** ```res_”ID”_longitudinal.out```

```
Subbasin-ID, year, day, hour, nbr. sections, minelev_sec(m)
60  1980   1   1  12     447.309998     445.579987     445.239990	…	430.570007     418.160004     414.519989
60  1980   1   2  12     448.229315     445.615664     445.240689	…	430.591036     418.168831     414.525120
…
```

*Subbasin-ID*: ID of sub-basin <br>
*year*: Year of simulation <br>
*day*: Day of simulation <br>
*hour*: Hour of simulation <br>
*nbr. sections*: Number of cross sections in the sub-basin's reservoir <br>
*minelev\_sec*: Minimum elevation at the cross section of the sub-basin's reservoir \[m]

Example: After each time step, e.g. after one day, the reservoir of the sub-basin with the ID 60 has 12 values of minimum bed elevation corresponding to the 12 cross sections. They are changed by either deposition or erosion processes (447.309998 m, 445,579987 m, 445.239990 m, etc). Currently, the model generates an output file for each reservoir considered in the simulation (e.g. ```res_60 _longitudinal.out``` referred to sub-basin with ID 60).

**8)** ```res_”ID”_sedcomposition.out```

```
Subbasin-ID, year, day, hour, nbr. classes, sedcomp_outflow(-)
60  1980   1   1   3          0.999          0.001          0.000
60  1980   1   2   3          0.999          0.001          0.000
…
```

*Subbasin-ID*: ID of sub-basin <br>
*year*: Year of simulation <br>
*day*: Day of simulation <br>
*hour*: Hour of simulation <br>
*nbr. classes*: Number of sediment size classes considered in the simulation <br>
*sedcomp\_outflow*: Effluent size distribution downstream the sub-basin's reservoir \[-].The total number of sediment size classes is previously specified in the file ```do.dat```.

Example: After each time step, e.g. after one day, the reservoir of the sub-basin with the ID 60 has the following effluent size distribution for the given sediment classes (e.g. three sediment classes): 0.999, 0.001, and 0.000. Currently, the model generates an output file for each reservoir considered in the simulation (e.g. ```res_60 _sedcomposition.out``` referred to sub-basin with ID 60).

**9)** ```lake_inflow_r.out```

```
Year, day, hour, inflow_r(m**3/timestep)
1980     1     1       2224.611       1511.890         11.295          0.000          0.000
1980     2     1       2098.507       1457.613         10.831          0.000          0.000
…
```

*Year*: Year of simulation <br>
*day*: Day of simulation <br>
*hour*: Hour of simulation <br>
*nbr. classes*: Number of sediment size classes considered in the simulation <br>
*inflow\_r*: Water inflow discharges into the reservoir size classes \[m<sup>3</sup>/timestep]. Currently, the number of reservoir size classes can not be changed (total of five classes)

Example: After each time step, e.g. after one day, five values of water inflow discharges into the reservoir size classes are computed (5748.602 m<sup>3</sup>, 2409.138 m<sup>3</sup>, 31.014 m<sup>3</sup>, 0.000 m<sup>3</sup> and 0.000 m<sup>3</sup> within the timestep for the size classes 1 to 5, respectively). Results are displayed for the whole catchment after grouping them by reservoir size classes. The files 10 to 16 have the same structure, as shown by the file ```lake_ inflow_r.out``` (file 9).

**17)** ```lake_ watbal.out```

```
Year, day, hour, totallakeinflow(m**3/timestep), totallakeoutflow(m**3/timestep
 ), totallakeprecip(m**3/timestep), totallakeevap(m**3/timestep), lakevol(m**3)
1980     1     1       3747.793          0.000          0.000         53.854      38678.957
1980     2     1       3566.949          0.000          0.000         63.270      42182.633
…
```

*Year*: Year of simulation <br>
*day*: Day of simulation <br>
*hour*: Hour of simulation <br>
*totallakeinflow*: Total water inflow discharge into all upstream reservoirs of the catchment \[m<sup>3</sup>/timestep] <br>
*totallakeoutflow*: Total water outflow discharge from all upstream reservoirs of the catchment \[m<sup>3</sup>/timestep] <br>
*totallakeprecip*: Total rainfall over all upstream reservoirs of the catchment \[m<sup>3</sup>/timestep] <br>
*totallakeevap*: Total evaporation from all upstream reservoirs of the catchment \[m<sup>3</sup>/timestep] <br>
*lakevol*: Total water volume stored in all upstream reservoirs of the catchment \[m<sup>3</sup>]

Example: After each time step, e.g. after one day, a total water inflow discharge into all upstream reservoir of 3747.793 m<sup>3</sup>/timestep, no water outflow discharge, no rainfall over the reservoir areas, a total evaporation of 53.854 m<sup>3</sup>/timestep, and a total water storage of 38678.957 m<sup>3</sup> in all upstream reservoirs. Results are displayed for the whole catchment without distinguishing between size classes.

**18)** ```lake_sedbal.out```

```
Year, day, hour, totalsedinflow(ton/timestep), totalsedoutflow(ton/timestep), totalsedimentation(ton/timestep), cumsedimentation(ton)
1980     1     1        200.000          0.000          200.000          200.000
1980     2     1        100.000          50.000          50.000          250.000
…
```
 
*Year*: Year of simulation <br>
*day*: Day of simulation <br>
*hour*: Hour of simulation <br>
*totalsedinflow*: Total sediment inflow discharge into all upstream reservoirs of the catchment \[ton/timestep] <br>
*totalsedoutflow*: Total sediment outflow discharge from all upstream reservoirs of the catchment \[ton/timestep] <br>
*totalsedimentation*: Total sediment deposition in all upstream reservoirs of the catchment \[ton/timestep] <br>
*cumsedimentation*: Cumulative sediment deposition in all upstream reservoirs of the catchment \[ton]

Example: After each time step, e.g. after one day, a total sediment inflow discharge into all upstream reservoir of 200 ton/timestep, no sediment outflow discharge, a total sediment deposition of 200 ton/timestep, and a cumulative sediment deposition of 200 ton in all upstream reservoirs. Results are displayed for the whole catchment without distinguishing between size classes.

**19)** ```lake_inflow.out```

```
Year, day, hour, reservoir_class, lakeinflow(m**3/timestep)
                                              57                 15                20                60
1980     1     1     1        384.741        587.248         38.144         17.718
1980     1     1     2        411.781          0.000         70.314         40.512
1980     1     1     3          0.000          0.000          0.000         11.295
1980     1     1     4          0.000          0.000          0.000          0.000
1980     1     1     5          0.000          0.000          0.000          0.000
…
```

*Year*: Year of simulation <br>
*day*: Day of simulation <br>
*hour*: Hour of simulation <br>
*nbr. classes*: Number of sediment size classes considered in the simulation <br>
*lakeinflow*: Water inflow discharges into the reservoir size classes. Currently, the number of reservoir size classes can not be changed (total of five classes)

Example: After each time step, e.g. after one day, after one day, values of water inflow discharges into the reservoir size classes are computed for all sub-basins (e.g. for size class 1: 384.741 m<sup>3</sup>, 587.248 m<sup>3</sup>, 38.144 m<sup>3</sup>, and 17.718 m<sup>3</sup> within the timestep, for sub-basins 57, 15, 20 and 60, respectively). Results are displayed for all sub-basins after grouping them by reservoir size classes. The files 20 to 25 have the same structure, as shown by the file ```lake_ inflow_r.out``` (file 9).

**20)** ```lake_sizedistoutflow.out```

```
Year, day, hour, sediment size class, lakesizedistoutflow(-)
                                              57                 15                 20                60
1980     1     1     1              0.60              0.40              0.40             0.50
1980     1     1     2              0.30              0.30              0.40             0.40
1980     1     1     3              0.10              0.30              0.40             0.10
…
```

*Year*: Year of simulation <br>
*day*: Day of simulation <br>
*hour*: Hour of simulation <br>
*nbr. classes*: Number of sediment size classes considered in the simulation <br>
*lakeinflow*: Effluent size distribution at the sub-basin outlet after sediment routing through the reservoir cascade \[-].The total number of sediment size classes is previously specified in the file ```do.dat```.

Example: After each time step, e.g. after one day, the sediment outflow discharge at the sub-basin outlet has the following effluent size distribution for the given sediment classes (e.g. three sediment classes): fifth column displays the results of grain size distribution for the sub-basin with ID 57 (0.60, 0.30 and 0.10, for sediment classes 1 to 3, respectively). Results are displayed for all sub-basins without distinguishing between size classes. 

\[[Table of contents](#toc)]
## 5 References

<a name="ackers-white-1973"></a>
Ackers, P. and White, W.R. (1973): Sediment transport: a new approach and analysis. Proc. ASCE, Journal of the Hydraulics Division, Vol. 99, HY11, pp. 2041-2060.

<a name="antronico-et-al-2005"></a>
Antronico, L., Coscarelli, R., Terranova, O. (2005): Surface erosion assessment in two Calabrian basins (southern Italy). In: R. J. Batalla and C. Garcia (Ed.), Geomorphological Processes and Human Impacts in River Basins, IAHS, pp. 16-22.

<a name="ashida-michiue-1973"></a>
Ashida, K. and Michiue, M. (1973): Studies on bed load transport rate in alluvial streams. Trans. Japan Society of Civil Engineers, Vol. 4.

<a name="breuer-et-al-2003"></a>
Breuer, L., Eckhardt, K., Frede, H.-G. (2003): Plant parameter values for models in temperate climates, Ecological Modelling, 169: 237-293.

<a name="bronstert-et-al-1999"></a>
Bronstert, A., Güntner, A., Jaeger, A., Krol, M., and Krywkow, J. (1999): Großräumige hydrologische Parameterisierung und Modellierung als Teil der integrierten Modellierung, pp. 31-40. In N. Fohrer and P. Döll, editors, Modellierung des Wasser- und Stofftransports in großen Einzugsgebieten. Kassel University Press, Kassel.

<a name="fao-1993"></a>
FAO (1993): Global and national soils and terrain digital databases (SOTER). Procedures Manual. World Soil Resources Reports, No. 74., FAO (Food and Agriculture Organization of the United Nations), Rome, Italy.

<a name="fao-2001"></a>
FAO (2001): Global Soil and Terrain Database (WORLD-SOTER). FAO, AGL (Food and AgricultureOrganization of the United Nations, Land and Water Development Division), http://www.fao.org/ag/AGL/agll/soter.htm.

<a name="francke-et-al-2008"></a>
Francke, T., Güntner, A., Bronstert, A., Mamede, G., Müller, E. N. (2008): Automated catena-based discretisation of landscapes for the derivation of hydrological modelling units. International Journal of Geographical Information Science 22: 111-132.

<a name="francke-2009"></a>
Francke, T. (2009): Measurement and Modelling of Water and Sediment Fluxes in Meso-Scale Dryland Catchments. PhD thesis, Universität Potsdam, Germany. http://opus.kobv.de/ubp/volltexte/2009/3152/. 

<a name="guentner-2002"></a>
Güntner, A. (2002): Large-scale hydrological modelling in the semi-arid North-East of Brazil. PIK-Report No. 77. Potsdam Institute for Climate Research, Germany (http://www.pik-potsdam.de/pik_web/ publications/pik_reports/reports/reports/pr.77/pr77.pdf).

<a name="guentner-bronstert-2002"></a>
Güntner, A., Bronstert, A. (2002): Process-based modelling of large-scale water availability in a semi-arid environment: process representation and scaling issues. In G.H. Schmitz, editor, Schriftenreihe des Institutes für Abfallwirtschaft und Altlasten, Universität Dresden, Dresden, pp. 46.

<a name="guentner-bronstert-2003"></a>
Güntner, A., Bronstert, A. (2003): Large-scale hydrological modeling of a semiarid environment: model development, validation and application, In T. Gaiser, M. Krol, H. Frischkorn, and J.C.Araujo, editors, Global change and regional impacts. Springer-Verlag, Berlin.

<a name="guentner-bronstert-2004"></a>
Güntner, A., Bronstert, A. (2004): Representation of landscape variability and lateral redistribution processes for large-scale hydrological modelling in semi-arid areas, Journal of Hydrology 297: 136-161.

<a name="krysanova-et-al-2000"></a>
Krysanova, V., Wechsung, F., Arnold, J., Srinivasan, R., Williams, J. (2000): SWIM (Soil and Water Integrated Model), User Manual. PIK Report Nr. 69, pp 239.

<a name="maidment-1993"></a>
Maidment, D. R. (1993): Handbook of hydrology. MGraw-Hill, New York.

<a name="mamede-2008"></a>
Mamede, G. (2008): Reservoir sedimentation in dryland catchments: Modelling and management. PhD thesis at the University of Potsdam, Germany, published on: https://publishup.uni-potsdam.de/opus4-ubp/files/1546/mamede_diss.pdf.

<a name="morgan-1995"></a>
Morgan, R.P.C. (1995): Soil erosion and conservation Longman Group, UK Limited. 

<a name="mueller-et-al-2009"></a>
Mueller, EN., Francke, T., Batalla, RJ., Bronstert, A. (2009): Modelling the effects of land-use change on runoff and sediment yield for a meso-scale catchment in the Southern Pyrenees. Catena 79:3, 288-296. (1)

<a name="mueller-et-al-2008"></a>
Mueller, E. N., Batalla, R. J., Garcia, C., Bronstert, A. (2008): Modelling bedload rates from fine grain-size patches during small floods in a gravel-bed river. J. of Hydr. Eng. in press.

<a name="mueller-et-al-2010"></a>
Mueller, E.N., Güntner, A., Francke, T., Mamede, G. (2010): Modelling sediment export, retention and reservoir sedimentation in drylands with the WASA-SED model. Geosci Model Dev 3:275–291, doi:10.5194/gmd-3-275-2010, published on: http://www.geosci-model-dev.net/3/275/2010/.

<a name="neitsch-et-al-2002"></a>
Neitsch, S.L., Arnold, J.G., Kiniry, J.R., Williams, J.R., King, K.W. (2002): Soil and Water Assessment Tool. Theoretical Documentation, Version 2000. Published by Texas Water Resources Institute, TWRI Report TR-191.

<a name="pilz-et-al-2017"></a>
Pilz, T., T. Francke, & A. Bronstert (2017): lumpR: An R package facilitating landscape discretisation for hillslope-based hydrological models, Geosci Model Dev, 10(8), 3001–3023, doi:https://doi.org/10.5194/gmd-10-3001-2017.

<a name="rottler-2017"></a>
Rottler, E. (2017): Implementation of a snow routine into the hydrological model WASA-SED and its validation in a mountainous catchment. MSc Thesis, University Potsdam, Germany. https://doi.org/10.25932/publishup-50496.

<a name="williams-1995"></a>
Williams, J. (1995): The EPIC Model. In: Singh, V. P. (Eds.): Computer Models of Watershed Hydrology. Water Resources Publications, Highlands Ranch, CO., pp. 909-1000.

<a name="wu-et-al-2000"></a>
Wu, W., Wang, S.S.Y., Jia, Y. (2000): Nonuniform sediment transport in alluvial rivers. Journal of Hydraulic Research, Vol. 38, No. 6, pp 427-434.

<a name="yang-1973"></a>
Yang, C.T. (1973): Incipient  motion  and  sediment transport.  Journal  of  the  Hydraulic  Division, ASCE, 99, no. HY 10, pp. 1679-1704. 

<a name="yang-1984"></a>
Yang, C.T. (1984): Unit stream power equation for gravel. Journal of Hydraulic Engineering, ASCE, 110 (12), pp. 1783-1797. 

\[[Table of contents](#toc)]
## 6 Further relevant literature and tools for the WASA-SED model

Bronstert, A., Jaeger, A., Güntner, A., Hauschild, M., Döll, P., and Krol, M. (2000): Integrated modelling of water availability and water use in the semi-arid Northeast of Brazil, Physics and Chemistry of the Earth 25: 227-232.

Mamede, G.L., Bronstert, A., Araujo, J.C., Batalla, R. J., Güntner, A., Mueller, E. N., Francke, T. (2006): 1D Process-Based Modelling of Reservoir Sedimentation: a Case Study for the Barasona Reservoir in Spain. Proceedings of the International Conference on Fluvial Hydraulics, Lisbon, Vol. 2: 1585-1594.

Güntner, A. (2003): Auswirkung von Klimaänderungen auf die Wasserverfügbarkeit in Trockengebieten - Ergebnisse und Unsicherheiten am Beispiel Nordost-Brasiliens. In H.-B.Kleeberg, editor, Hydrologische Wissenschaften - Fachgemeinschaft in der ATV-DVWK, pp. 205-214.

Güntner, A., Krol, M., Araujo, J.C., Bronstert, A. (2004): Simple water balance modelling of surface reservoir systems in a large data-scarce semiarid region, Hydrological Sciences Journal 49: 901-918.

**For WASA-SED:**

Mueller, EN., Guentner, A., Francke, T., Mamede, GL. (2010): Modelling sediment export, retention and reservoir sedimentation in drylands with the WASA-SED Model.  Geoscientific Model Development 3, 275-291.

Mamede, GL. (2008): Reservoir sedimentation in dryland catchments: Modelling and management. PhD thesis, Universität Potsdam, Germany. http://opus.kobv.de/ubp/volltexte/2008/1704/.


**For WASA-SED parameterisations:**

Appel, K. (2006): Characterisation of badlands and modelling of soil erosion in the Isabena watershed, NE Spain. Unpublished MSc thesis, University of Potsdam, Germany.

Francke, T. (2009): Measurement and Modelling of Water and Sediment Fluxes in Meso-Scale Dryland Catchments. PhD thesis, Universität Potsdam, Germany. http://opus.kobv.de/ubp/volltexte/2009/3152/. 

Francke, T., Parameterisation of the Esera/Isabena Catchment, Pre-Pyrenees, Spain. SESAM Working Report, http://brandenburg.geoecology.uni-potsdam.de/projekte/sesam/publications.php.

Medeiros, PHA. (2009): Hydro-sedimentological processes and connectivity in a semiarid basin: modelling and validation in several scales. PhD thesis, Universidade Federal do Ceará, Brazil. http://www.teses.ufc.br/tde_busca/arquivo.php?codArquivo=4425.

Medeiros, P.H.A., Guentner, A., Francke, T., Mamede, GL., De Araújo, JC. (2010): Modelling spatio-temporal patterns of sediment yield and connectivity in a semi-arid catchment with the WASA-SED model. Hydrological Sciences Journal 55:4, 636-648. (1)

Mueller, EN., Francke, T., Batalla, RJ., Bronstert, A. (2009): Modelling the effects of land-use change on runoff and sediment yield for a meso-scale catchment in the Southern Pyrenees. Catena 79:3, 288-296. (1)

Mueller E. N., Batalla, R. J., Garcia, C., Bronstert, A. (2008): Modelling bedload rates from fine grain-size patches during small floods in a gravel-bed river. J. of Hydr. Eng. in press. 

**For the hydrological modules:**

Güntner, A. (2002): Large-scale hydrological modelling in the semi-arid North-East of Brazil. PIK-Report No. 77, Potsdam Institute for Climate Research, Germany.

Güntner, A. and Bronstert, A. (2004): Representation of landscape variability and lateral redistribution processes for large-scale hydrological modelling in semi-arid areas. Journal of Hydrology, 297: 136-161.

**For LUMP:**

Francke, T., Güntner, A., Bronstert, A., Mamede, G., Müller, E. N. (2008): Automated catena-based discretisation of landscapes for the derivation of hydrological modelling units. International Journal of Geographical Information Science, 22: 111-132.

Pilz, T., T. Francke, & A. Bronstert (2017): lumpR: An R package facilitating landscape discretisation for hillslope-based hydrological models, Geosci Model Dev, 10(8), 3001–3023, doi:https://doi.org/10.5194/gmd-10-3001-2017.

Pilz, T (2015): https://github.com/tpilz/LUMP.


**Auxiliary Tools:**

[SoilDataPrep](https://github.com/tillf/SoilDataPrep)<a name="SoilDataPrep"></a> : soildata-preprocessing for WASA-SED

[lumpR](https://github.com/tpilz/lumpR) : Geodata-preprocessing for WASA-SED

[MMEMO](https://github.com/TillF/MMEMO)<a name="MMEMO"></a> : assess the impact of various model enhancements, including routines for calibration. 
