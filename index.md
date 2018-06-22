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

<!-- toc -->

## Introduction
The WASA-SED model simulates the runoff and erosion processes at the hillslope scale, the transport processes of suspended and bedload fluxes at the river scale and the retention and remobilisation processes of sediments in large reservoirs. The modelling tool enables the evaluation of management options both for sustainable land-use change scenarios to reduce erosion in the headwater catchments as well as adequate reservoir management options to lessen sedimentation in large reservoirs and reservoir networks. The model concept, its spatial discretisation and the numerical components of the hillslope, river and reservoir processes are summarised and current model applications are reviewed in Mueller et al. (2008). The hydrological routines of the model are based on the WASA model (Model for Water Availability in Semi-Arid environments), which was developed by Güntner (2002) and Güntner and Bronstert (2002, 2003) to enable the quantification of water availability in semi-arid regions. The WASA-SED model was developed within the joint Spanish-Brazilian-German research project SESAM (Sediment Export from Semi-Arid Catchments: Measurement and Modelling). The existing WASA model code has been extended to include sediment-transport routines for the three new conceptual levels of the WASA-SED model: the hillslope scale, river scale and the reservoir scale for the calculation of sedimentation. This documentation gives a short outline of the structure, computational routines and folder system of the WASA-SED code in Chapter 2, followed by a description of the input files for model parameterisation in Chapter 3 and output files for the hillslope, river and reservoir modules in Chapter 4.

## Program folders and structure
The WASA-SED model is programmed in Fortran 90 and was tested with Compaq Visual Fortran (6.6.a) and gfortran 4.x compilers on Windows gfortran-4.3 under Linux. The WASA-SED code is organised as presented in Table 1. The main program is named wasa.f90 and calls the key subroutines as summarised in Table 2.

Table 1 Folder structure of WASA code

Folder Name | Content
---|---
\General\ | Main program and utility routines of the WASA model  
\Hillslope\ | Hillslope routines (Overland, sub-surface flow, evapotranspiration, sediment production etc.)  
\River\	| River routines (Routing of water and sediment in the river network)
\Reservoir\ | Reservoir water and sediment modelling routines (Water balance, bed elevation change, management options)
\Input\	| Input data, contains parameter file: do.dat, \Input\Hillslope, \Input\River, \Input\Reservoir, \Input\Time_series
\Output\ | Output files of model scenarios 

A complete example input-data set is part of the Subversion repository (see top of document). Model parameterisations are available e.g. for meso-scale catchments in dryland areas of Spain and Brazil (Bengue Catchment in Spain: Mamede 2008, Ribera Salada Catchment in Spain: Mueller et al. 2008, Mueller et al. submitted to CATENA, Isábena Catchment in Spain: Francke 2009).

The original WASA code version (Güntner 2002, Güntner and Bronstert 2004) was extended within the SESAM-Project to include sediment-transport processes at the hillslope scale using various USLE-derivative approaches, a spatially distributed, semi-process-based modelling approach for the modelling of water and sediment transport through the river network and a reservoir module that computes the transport of water and sediment as well as sedimentation processes in reservoirs.

Table 2 Main subroutines of the main program (wasa.f90)

Routine	| Content
---|---
readgen.f90	| Reads parameter file: do.dat
calcyear.f90	| Counter for time
calcday.f90	| Number of days per month
climo.f90	| Climate calculations
hymo_all.f90	| Hydrological 
routing.f90 |	Old river routing routine without sediment transport (unit hydrograph)
routing_new.f90 | New river routine with Muskingum and sediment transport, calls reservoir routines

The following sections give some information on the computational background and the functional structure of the corresponding routines for each of the three conceptual levels: hillslope, river and reservoir.

## Hillslope module
The hillslope module comprises the modelling of the hydrological and sediment-transport processes. The hydrological modelling accounts for interception, evaporation, infiltration, surface and subsurface runoff, transpiration and ground water recharge. Details are given in Güntner (2002, Chapter 4). The main hydrological calculations are carried out in hymo_all.f90 (for daily or hourly time steps). The subroutines that are called within hymo_all.f90 are summarised in Table 3. The temporal sequence of hydrological process modelling is summarised in Güntner (p. 36-37).

Table 3 Main subroutines of hymo_all.f90 (hydrological subroutines)

Routine	| Content	| Subroutine	| Content of Subroutine
---|---|---|---
readhymo.f90 |	reads input data	| soildistr.f90	| distribution functions for soil parameters
lake.f90	| water balance for small reservoirs	| - |	-
soilwat.f90	| soil water components (infiltration, percolation, runoff generation)	|etp_max.f90, etp_soil.f90, sedi_yield.f90	| maximal evapotranspiration, daily evapotranspiration, hillslope erosion

Sediment generation on the hillslopes in the form of soil erosion by water is modelled with four erosion equations (USLE, Onstad-Foster, MUSLE and MUST after Williams, 1995). The following subroutines contain the main calculations for the water and sediment production for the hillslopes:

*hymo_all.f90:* general loop in daily or hourly timesteps
1.	Initialisation and reading of hillslope input data, preparation of output files
2.	For each timestep
    -	nested loop that contains the calculation for hydrological and sediment budgets for each sub-basin, landscape unit, terrain component, soil-vegetation component (soilwat.f90)
    -	aggregated runoff refers to the whole available area including reservoir areas; the corresponding reservoir areal fraction is substracted afterwards if reservoirs / small reservoirs are considered
    -	call module for small (diffuse, unlocated) reservoirs (lake.f90)
    - write daily results into output files

*soilwat.f90:* computes water and sediment balance for a terrain component in a given landscape unit of a given sub-basin
1.	For each soil vegetation component (SVCs) in the terrain component: calculate hydrological variables (infiltration, surface runoff, subsurface runoff, evapotranspiration, etc.)
2.	Re-distribute surface runoff, subsurface runoff among SVCs, allow re-infiltration and compute overall water balance for 
terrain component
3.	Estimate runoff height and peak-runoff using Manning Equation
4.	Call module for hillslope erosion (sedi_yield.f90)

*sedi_yield.f90:* computes sediment production and sediment export for a TC
1.	No sediment export without surface runoff out of terrain component
2.	Otherwise, compute gross erosion for current terrain component using one of four erosion equations
3.	Distribute gross erosion among particle classes according to constitution of uppermost horizons in the current terrain component
4.	Completely mix the newly generated sediment with any sediment coming from upslope terrain components
5.	If necessary, limit sediment export by transport capacity.

## River module
The river routing of the original WASA model (Güntner 2002) bases on daily linear response functions (Bronstert et al. 1999) similar to a triangular unit hydrograph. Its implementation does not support output in hourly resolution (only daily is produced) and sediment transport. It was extended to include a spatially semi-distributed, semi-process-based modelling approach for the modelling of water and sediment transport through the river network. The implemented water modelling approach is similar to the routing routines from the SWAT model (Soil Water Assessment Tool, Neitsch et al. 2002) model and the SWIM model (Soil Water Integrated Modelling, Krysanova et al. 2000). The new water routing is based on the Muskingum kinematic wave approximation. Suspended sediment transport and bedload is modelled using the transport capacity concept. The river module can be run with variable time steps. Transmission losses through riverbed infiltration and evaporation are accounted for. The main routing calculations as well as the initialisation and reading of the river input files are carried out in routing.f90. The following sub-routines are called from routing.f90:

*muskingum.f90:* contains the flow calculation using the Muskingum method
1. Calculation of water volume in reach
2. Calculation of cross-section area of current flow, flow depth, wetted perimeter and hydraulic radius
3. Calculation of flow in reach with Manning Equation
4. Calculation of Muskingum coefficients
5. Calculation of discharge out of the reach and water storage in reach at end of time step

*routing_coefficients.f90:* contains the calculations for travel time and flow depth
1. Calculation of initial water storage for each river stretch
2. Calculation of channel dimensions
3. Calculation of flow and travel time at 100 % bankfull depth, 120 % bankfull depth and 10 % bankfull depth

*route_sediments.f90:* contains the calculation for suspended sediment-transport routing
1. Calculation of water volume and sediment mass in reach
2. Calculation of peak flow and peak velocity
3. Calculation of maximum sediment carrying capacity concentration
4. Comparison with current concentration
5. Calculation of net deposition and degradation
6. Calculation of sediment mass out of the reach and sediment storage in reach at end of time step

*bedload.f90:* contains the calculation for bedload transport using 5 different formulas
1.	Calculation of current width of the river
2.	Bedload formulas after Meyer-Peter & Müller (1948), Schoklitsch (1950), Smart & Jaeggi (1983), Bagnold (1956) und Rickenmann (2001), see Mueller et al. 2008 relevant references

## Reservoir module
The reservoir sedimentation routine was included into the WASA-SED model by Mamede (2008) to enable the calculation of non-uniform sediment transport along the longitudinal profile of a reservoir, of the reservoir bed changes caused by deposition/erosion processes and of reservoir management options. 
In order to perform the simulation of sediment transport in reservoirs, four important processes have to be considered: (1) reservoir water balance, (2) hydraulic calculations in the reservoir, (3) sediment transport along the longitudinal profile of the reservoir and (4) reservoir bed elevation changes. For the calculation of sediment transport in the reservoir, four different equations for the calculation of total sediment load were selected from recent literature. The reservoir bed elevation changes are calculated through the sediment balance at each cross section, taking into account three conceptual sediment layers above the original bed material. The reservoir sedimentation module is composed by the following subroutines:

*reservoir.f90:* contains the reservoir water balance
1. Calculation of reservoir level
2. Calculation of controlled and uncontrolled outflow

*semres.f90:* contains the sediment balance for each cross section
1. Calculation of sediment deposition for each cross section
2. Calculation of sediment entrainment for each cross section
3. Calculation of sediment compaction for each cross section
4. Calculation of total load transport
5. Calculation of storage capacity reduction because of the accumulated sediment
6. Calculation of effluent grain size distribution

*sedbal.f90:* contains a simplified sediment balance of the reservoir whenever its geometry is not provided (no cross-section)
1. Calculation of sediment deposition in the reservoir
2. Calculation of trapping efficiency
3. Calculation of storage capacity reduction because of the accumulated sediment
4. Calculation of effluent grain size distribution

*hydraul_res.f90:* contains the hydraulic calculations for each cross section in the reservoir
1. Calculation of mean velocity of current flow for each cross section
2. Calculation of hydraulic radius for each cross section
3. Calculation of slope of energy-grade line for each cross section
4. Calculation of top width for each cross section
5. Calculation of water depth for each cross section

*eq1_wu.f90, eq2_ashida.f90, eq3_tsinghua.f90* and *eq4_ackers.f90:* each sub-routine contains a specific sediment transport function, mentioned in Table 4
1. Calculation of fractional bed load transport for each cross section
2. Calculation of fractional suspended load transport for each cross section

*change_sec.f90:* contains the calculation of reservoir geometry changes
1. Calculation of bed elevation changes for each cross section.
*reservoir_routing.f90:* contains the calculation of level-pool reservoir routing 
1. Calculation of water routing for reservoirs with uncontrolled overflow spillways
*vert_dist.f90:* contains the calculation of vertical distribution of sediment concentration in the reservoir 
1. Calculation of sediment concentration at the reservoir outlets
*lake.f90:* contains the water balance for networks of small reservoirs 
1. Calculation of water level in the reservoirs
2. Calculation of controlled and uncontrolled outflow out the small reservoirs
*sedbal_lake.f90:* contains a simplified sediment balance for networks of small reservoirs
1. Calculation of sediment deposition in small reservoir
2. Calculation of trapping efficiency in small reservoirs
3. Calculation of storage capacity reduction because of the accumulated sediment in small reservoirs
4. Calculation of effluent grain size distribution in small reservoirs
*lake_routing.f90:* contains the calculation of level-pool routing for networks of small reservoirs
1. Calculation of water routing for small reservoirs

## Input Data
The model runs as a Fortran console application for catchment from a few km² up to several 100,000 km²) on daily or hourly time steps. Climatic drivers are daily/hourly time series for precipitation, humidity, short-wave radiation and temperature. For model parameterisation, regional digital maps on soil associations, land-use and vegetation cover, a digital elevation model with a cell size of 100 metres (or smaller) and, optional, bathymetric surveys of the reservoirs are required. The soil, vegetation and terrain maps are processed with the LUMP tool (see above) to derive the spatial discretisation into soil-vegetation units, terrain components and landscape units. Table 4 summarises the input parameters for the climatic drivers and the hillslope, river and reservoir modules. The vegetation parameters may be derived with the comprehensive study of, for example, Breuer et al. (2003), the soil and erosion parameters with the data compilation of FAO (1993, 2001), Morgan (1995), Maidment (1993) and Antronico et al. (2005).

For a semi-automated discretisation of the model domain into landscape units and terrain components, the software tool LUMP (Landscape Unit Mapping Program) is available (Francke et al. 2008). LUMP incorporates an algorithm that delineates areas with similar hillslope characteristics by retrieving homogeneous catenas with regard to e.g. hillslope shape, flow length and slope (provided by a digital elevation model), and additional properties such as for soil and land-use and optionally for specific model parameters such as leaf area index, albedo or soil aggregate stability. The LUMP tool is linked with the WASA-SED parameterisation procedure through a databank management tool, which allows to process and store digital soil, vegetation and topographical data in a coherent way and facilitates the generation of the required input files for the model. LUMP and further WASA-SED pre-processing tools have been transferred to the package lumpR for the free software environment for statistical computing and graphics R which is available from https://github.com/tpilz/lumpR.

The input files for general purpose, the hillslope, river and reservoir routines are explained below with details on parameter type, units, data structure including examples parameterisation files.

Table 4: Summary of model input parameters

Type |	Model input parameter
---|---
Climate |	Daily or hourly time series on rainfall [mm/day, mm/h], Daily time series for average short-wave radiation [W/m2], Daily time series for humidity [%], Daily time series for temperature [°C]
Vegetation	| Stomata resistance [s/m], Minimum suction [hPa], Maximum suction [hPa], Height [m], Root depth [m], LAI [-], Albedo [-], Manning's n of hillslope[-], USLE C [-]
Soil	| No. of horizons\*, Residual water content [Vol. %], Water content at permanent wilting point [Vol. %], Usable field capacity [Vol. %], Saturated water content [Vol. %], Saturated hydraulic conductivity (mm/h), Thickness [mm], Suction at wetting front [mm], Pore size index [-], Bubble pressure [cm], USLE K [-]\*\*, Particle size distribution\*\*
Terrain and river |	Hydraulic conductivity of bedrock [mm/d], Mean maximum depth of soil zone [mm], Depth of river bed below terrain component [mm], Initial depth of groundwater below surface [mm], Storage coefficient for groundwater outflow [day], Bankful depth of river [m], Bankful width of river [m], Run to rise ratio of river banks [-], Bottom width of floodplain [m], Run to rise ratio of floodplain side slopes [-], River length [km], River slope [m/m], D<sub>50</sub> (median sediment particle size) of riverbed [m], Manning’s n for riverbed and floodplains [-]
Reservoir	| Longitudinal profile of reservoir [m], Cross-section profiles of reservoir [m], Stage-volume curves, Initial water storage and storage capacity volumes [m3], Initial area of the reservoir [ha], Maximal outflow through the bottom outlets [m3/s], Manning’s roughness for reservoir bed, Depth of active layer [m]

\* for each soil horizon, all following parameters in the column are required, \** of topmost horizon

## General parameter and control files
Four parameter files control the data input and output and some internal settings:

*do.dat* \[can be generated with The LUMP package, manual completing required\]
The do.dat file is located in the folder WASA\Input and contains the main parameter specifications for the WASA-SED model. Figure 1 displays an example file for the do.dat. The first line of the do.dat contains the title. Line 2 and 3 specify the path for the location of WASA input and output folder. Relative paths are supported. The backslash “\” only works on Windows-platforms. The slash “/” is accepted on Windows and Unix/Linux systems. Make sure that both specified paths end with slash or backslash, respectively. Line 4 and 5 contain the start and the end year of the simulation, respectively. Line 6 and 7 contain the start and the end calendar month of the simulation, respectively. Optionally, the day of month for begin and end can be specified. Line 10 contains the number of sub-basins. The number in line 9 is given by the sum of the number of terrain components in each landscape-unit of each sub-basin (e.g. if the system has only two sub-basins, sub-basin A has 1 landscape unit with 3 terrain components, sub-basin B has 2 landscape units with 1 terrain component each, then the number of combinations is 5). Line 14 specifies if the reservoir module is switched on (.t.) or is switched off (.f.). The same issue for the calculations of networks of small reservoirs in line 15. Lines 16 – 19 allow customizing the way water and sediment is (re-)distributed within and among the TCs. Line 21 allows the setting of the simulation timestep (daily / hourly). This may become obsolete in future versions by setting the timestep directly in line 30. Line 24 allows specifying a correction factor for hydraulic conductivity to account for intra-daily rainfall intensities. Optionally, this factor can also be made a function of daily rainfall by specifying two more parameters (a and b) in the same line, so that kfkorr=kfkorr0\*(a\*1/daily_precip+b+1). In line 31 the erosion and sediment-transport routines may be switched on and off. Specify the number of grain size classes you want to model in line 32. Their limits must be specified in part_class.dat, if more than one class is desired. Line 33 lets you choose the hillslope erosion model to be used in WASA. Currently, this parameter is disregarded, further options can be chosen in erosion.ctl. Select the model for the river routing in line 34. Possible options are: (1) old WASA routing (daily resolution only, no sediment transport), (2) Muskingum and suspended sediment, (3) Muskingum and bedload transport. Choose the sediment model in the reservoir in line 35 among 4 sediment transport equations: (1) Wu et al. (2000); (2) Ashida & Michiue (1973); (3) Yang (1973, 1984); (4) Ackers & White (1973).

The optional lines 36 and 37 allow the saving/loading of state variables (i.e. groundwater, interception and soil storages) at the end/beginning of a model run (works only if svc.dat has been specified).

Line 36 may additionally contain a second logical variable (append_output), allowing the model to append to existing output files (not yet implemented for all outputs). The consistence with the existing files is not checked! Default is .FALSE.

Line 37 may additionally contain a second logical variable (save_states_yearly), determining if the model states are saved (and overwritten) at the end of each simulation year. Default is .TRUE.

insert here: Figure 1 -> table? figure?...





# Relevant Literature for the WASA-SED Model
**For WASA-SED:**<br>
Mueller, EN., Guentner, A., Francke, T., Mamede, GL. (2010): Modelling sediment export, retention and reservoir sedimentation in drylands with the WASA-SED Model.  Geoscientific Model Development 3, 275-291.

Mamede, GL. (2008): Reservoir sedimentation in dryland catchments: Modelling and management. PhD thesis, Universität Potsdam, Germany. http://opus.kobv.de/ubp/volltexte/2008/1704/

**For WASA-SED parameterisations:**<br>
Appel, K., 2006. Characterisation of badlands and modelling of soil erosion in the Isabena watershed, NE Spain. Unpublished MSc thesis. University of Potsdam, Germany.

Medeiros, PHA. (2009): Hydro-sedimentological processes and connectivity in a semiarid basin: modelling and validation in several scales. PhD thesis, Universidade Federal do Ceará, Brazil. http://www.teses.ufc.br/tde_busca/arquivo.php?codArquivo=4425

Medeiros, PHA., Guentner, A., Francke, T., Mamede, GL., De Araújo, JC. (2010): Modelling spatio-temporal patterns of sediment yield and connectivity in a semi-arid catchment with the WASA-SED model. Hydrological Sciences Journal 55:4, 636-648. (1)

Mueller, EN., Francke, T., Batalla, RJ., Bronstert, A. (2009): Modelling the effects of land-use change on runoff and sediment yield for a meso-scale catchment in the Southern Pyrenees. Catena 79:3, 288-296. (1)

Mueller E. N., Batalla, R. J., Garcia, C., Bronstert, A., 2008. Modelling bedload rates from fine grain-size patches during small floods in a gravel-bed river. J. of Hydr. Eng. in press 

Francke, T. (2009): Measurement and Modelling of Water and Sediment Fluxes in Meso-Scale Dryland Catchments. PhD thesis, Universität Potsdam, Germany. http://opus.kobv.de/ubp/volltexte/2009/3152/ 

**For the hydrological modules:**<br>
Güntner, A., 2002. Large-scale hydrological modelling in the semi-arid North-East of Brazil. PIK-Report No. 77, Potsdam Institute for Climate Research, Germany.
Güntner, A. and Bronstert, A., 2004. Representation of landscape variability and lateral redistribution processes for large-scale hydrological modelling in semi-arid areas, Journal of Hydrology, 297: 136-161.

**For LUMP:**<br>
Francke, T., Güntner, A., Bronstert, A., Mamede, G., Müller, E. N., 2008. Automated catena-based discretisation of landscapes for the derivation of hydrological modelling units, International Journal of Geographical Information Science, 22: 111-132.

**For LUMP package:**<br>
Pilz, T (2015): https://github.com/tpilz/LUMP

<br>
<br>

# References
Antronico, L., Coscarelli, R., Terranova, O., 2005. Surface erosion assessment in two Calabrian basins (southern Italy). In: R. J. Batalla and C. Garcia (Ed.), Geomorphological Processes and Human Impacts in River Basins, IAHS, pp. 16-22.

Appel, K., 2006. Characterisation of badlands and modelling of soil erosion in the Isabena watershed, NE Spain. Unpublished MSc thesis. University of Potsdam, Germany 

Ashida, K. and Michiue, M. 1973. “Studies on bed load transport rate in alluvial streams”, Trans. Japan Society of Civil Engineers, Vol. 4.

Ackers, P. and White, W.R. 1973. “Sediment transport: a new approach and analysis”, Proc. ASCE, Journal of the Hydraulics Division, Vol. 99, HY11, pp. 2041-2060.

Breuer, L., Eckhardt, K., Frede, H.-G., 2003. Plant parameter values for models in temperate climates, Ecological Modelling, 169: 237-293.

Bronstert, A., Güntner, A., Jaeger, A., Krol, M., and Krywkow, J. 1999. Großräumige hydrologische Parameterisierung und Modellierung als Teil der integrierten Modellierung, pp. 31-40. In N. Fohrer and P. Döll, editors, Modellierung des Wasser- und Stofftransports in großen Einzugsgebieten. Kassel University Press, Kassel

Bronstert, A., Jaeger, A., Güntner, A., Hauschild, M., Döll, P., and Krol, M. 2000. Integrated modelling of water availability and water use in the semi-arid Northeast of Brazil, Physics and Chemistry of the Earth 25: 227-232

Francke, T., Parameterisation of the Esera/Isabena Catchment, Pre-Pyrenees, Spain. SESAM Working Report, http://brandenburg.geoecology.uni-potsdam.de/projekte/sesam/publications.php 

Francke, T., Güntner, A., Bronstert, A., Mamede, G., Müller, E. N., 2008. Automated catena-based discretisation of landscapes for the derivation of hydrological modelling units. International Journal of Geographical Information Science 22: 111-132.

Francke, T., 2005. LUMP package, Manual, Auxiliary software tool to generate the input files for the hillslope module of the WASa model, SESAM working reports on http://brandenburg.geoecology.uni-potsdam.de/projekte/sesam/publications.php

FAO 1993. Global and national soils and terrain digital databases (SOTER). Procedures Manual. World Soil Resources Reports, No. 74., FAO (Food and Agriculture Organization of the United Nations), Rome, Italy.

FAO 2001. Global Soil and Terrain Database (WORLD-SOTER). FAO, AGL (Food and AgricultureOrganization of the United Nations, Land and Water Development Division), http://www.fao.org/ag/AGL/agll/soter.htm.

Güntner, A., 2002. Large-scale hydrological modelling in the semi-arid North-East of Brazil. PIK-Report No. 77. Potsdam Institute for Climate Research, Germany (http://www.pik-potsdam.de/pik_web/ publications/pik_reports/reports/reports/pr.77/pr77.pdf)

Güntner, A., Bronstert, A., 2002. Process-based modelling of large-scale water availability in a semi-arid environment: process representation and scaling issues. In G.H. Schmitz, editor, Schriftenreihe des Institutes für Abfallwirtschaft und Altlasten, Universität Dresden, Dresden, pp. 46

Güntner, A., Bronstert, A., 2003. Large-scale hydrological modeling of a semiarid environment: model development, validation and application, In T. Gaiser, M. Krol, H. Frischkorn, and J.C.Araujo, editors, Global change and regional impacts. Springer-Verlag, Berlin

Güntner, A., Bronstert, A., 2003. Large-scale hydrological modelling in the semiarid Northeast of Brazil: aspects of model sensitivity and uncertainty, In E. Servat, W. Najem, C. Leduc, and A. Shakeel, editors, Hydrology of the Mediterranean and Semi-Arid Regions. IAHS-Publication 278

Güntner, A., 2003. Auswirkung von Klimaänderungen auf die Wasserverfügbarkeit in Trockengebieten - Ergebnisse und Unsicherheiten am Beispiel Nordost-Brasiliens. In H.-B.Kleeberg, editor, Hydrologische Wissenschaften - Fachgemeinschaft in der ATV-DVWK, pp. 205-214
Güntner, A., Bronstert, A., 2004. Representation of landscape variability and lateral redistribution processes for large-scale hydrological modelling in semi-arid areas, Journal of Hydrology 297: 136-161

Güntner, A., Krol, M., Araujo, J.C., and Bronstert, A. 2004. Simple water balance modelling of surface reservoir systems in a large data-scarce semiarid region, Hydrological Sciences Journal 49: 901-918

IRTCES, 1985. Lecture notes of the training course on reservoir sedimentation. International Research of Training Center on Erosion and Sedimentation, Sediment Research Laboratory of Tsinghua University, Beijing, China.

Krysanova, V., Wechsung, F., Arnold, J., Srinivasan, R.., Williams, J., 2000. SWIM (Soil and Water Integrated Model), User Manual. PIK Report Nr. 69, pp 239.

Medeiros, PHA., Guentner, A., Francke, T., Mamede, GL., De Araújo, JC. (2010): Modelling spatio-temporal patterns of sediment yield and connectivity in a semi-arid catchment with the WASA-SED model. Hydrological Sciences Journal 55:4, 636-648. (1)

Maidment, D. R., 1993. Handbook of hydrology. MGraw-Hill, New York.

Mamede, G., 2008. Reservoir sedimentation in dryland catchments: Modelling and management. PhD thesis at the University of Potsdam, Germany, published on: urn:nbn:de:kobv:517-opus-17047.

Mamede, G.L., Bronstert, A., Araujo, J.C., Batalla, R. J., Güntner, A., Mueller, E. N., Francke, T. 2006. 1D Process-Based Modelling of Reservoir Sedimentation: a Case Study for the Barasona Reservoir in Spain. Proceedings of the International Conference on Fluvial Hydraulics, Lisbon, Vol. 2: 1585-1594.

Morgan, R.P.C., 1995. Soil erosion and conservation Longman Group, UK Limited. 

Mueller, EN., Francke, T., Batalla, RJ., Bronstert, A. (2009): Modelling the effects of land-use change on runoff and sediment yield for a meso-scale catchment in the Southern Pyrenees. Catena 79:3, 288-296. (1)

Mueller, E. N., Batalla, R. J., Garcia, C., Bronstert, A., 2008. Modelling bedload rates from fine grain-size patches during small floods in a gravel-bed river. J. of Hydr. Eng. in press

Mueller, E. N., Güntner, A., Francke, T., Mamede, G., 2008. Modelling water availability, sediment export and reservoir sedimentation in drylands with the WASA-SED Model. submitted to Geoscientific Model Development

Neitsch, S.L., Arnold, J.G., Kiniry, J.R., Williams, J.R., King, K.W., 2002. Soil and Water Assessment Tool. Theoretical Documentation, Version 2000. Published by Texas Water Resources Institute, TWRI Report TR-191

Renard, K.G., Foster, G.R., Weesies, G.A., McCool, D.K. and Yoder, D.C., 1997. Renard K, Foster G, Weesies G, McCool D, Yoder D 1997. Predicting soil loss by water: A guide to conservation planning with the Revised Universal Soil Loss Equation (RUSLE). U.S. Dep. of Agriculture, Agriculture Handbook 703.

Williams, J., 1995. The EPIC Model. In: Singh, V. P. (Eds.), Computer Models of Watershed Hydrology. Water Resources Publications, Highlands Ranch, CO., pp. 909-1000.

Wu, W., Wang, S.S.Y. and Jia, Y. 2000, “Nonuniform sediment transport in alluvial rivers”, Journal of Hydraulic Research, Vol. 38, No. 6, pp 427-434.

Yang, T.C. and Simoes, F.J.M. 2002 User’s Manual for GSTARS3 (Generalized Sediment Transport model for Alluvial River Simulation version 3.0). U.S. Department of the Interior, Bureau of Reclamation, Technical Service Center, Denver, Colorado


