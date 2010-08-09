! Till: computationally irrelevant: added program version information to parameter.out
! 2009-06-17

! Till: added output for River_Sediment_Storage.out
! 2008-11-13

! 2008-07-11
! Till: implemented optional pre-specified sediment outflow of selected subbasins

! 2008-07-03
! Till: implemented optional pre-specified outflow of selected subbasins

! 2008-01-31 
! Till: included variables for loading and saving of initial conditions
! included variables for controling evaporation calculation

! 2007-10-18
! Till: increased length of pathp variable to 160

! 2007-06-04
! Till: added flag f_tc_theta

! 2007-04-28
! Till: added flag f_tc_theta

! 2007-01-10
! Till: renamed f_deep_gw to f_deep_gw_recharge, added f_deep_gw_discharge

!Till: added variable f_deep_gw as flag for fileoutput

!
!This file contains four modules:
! common_h
! params_h
! time_h
! climo_h


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
module common_h
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
save
! logicals for inclusion of modules in calculation
CHARACTER (LEN=80) ::rev_string1,rev_string2	!Till: revision information string
LOGICAL :: doacud
LOGICAL :: doreservoir
LOGICAL :: dolattc
LOGICAL :: doalllattc
LOGICAL :: dolatsc
LOGICAL :: dolatscsub
LOGICAL :: dotrans
LOGICAL :: docellprec
LOGICAL :: dohour
LOGICAL :: doscale
LOGICAL :: domuncell
LOGICAL :: doreserv
LOGICAL :: do_pre_outflow=.TRUE.	!Till: consider pre-specified outflow time series from subbasins
LOGICAL :: do_pre_outsed=.TRUE.	!Till: consider pre-specified sediment outflow time series from subbasins

INTEGER :: ezgt
INTEGER :: zahlt
INTEGER :: nameti
INTEGER :: nmunt
INTEGER :: nroutt
INTEGER :: batch
INTEGER :: scenario
INTEGER :: krig
REAL :: sensfactor
CHARACTER (LEN=10) :: namet
CHARACTER (LEN=10) :: caset
CHARACTER (LEN=160) :: pfadn
CHARACTER (LEN=160) :: pfadp
INTEGER :: pfadi
INTEGER :: pfadj
INTEGER :: ncaset
INTEGER :: nezgt
INTEGER :: pyear
INTEGER :: dt					!time step in [hours]
LOGICAL :: dosediment			!do sediment calculation for hillslope, river and reservoir
INTEGER :: n_sed_class			!number of particles classes
INTEGER :: hill_transport		!type of hillslope sediment transport
INTEGER :: river_transport		!type of river sediment transport
INTEGER :: reservoir_transport	!type of reservoir sediment transport
INTEGER :: nt					!number of timesteps per day (computed as 24/dt)

real, allocatable :: upper_limit(:)			!upper limits of particle size class intervalls (mm)
real, allocatable :: particle_classes(:)	!upper limits of particle size classes [mm]

!flags for selecting which output files will be created
LOGICAL :: f_daily_actetranspiration,f_daily_potetranspiration, f_daily_qhorton,f_daily_qin_m3s,f_daily_qout_m3s,f_daily_rain,f_daily_runoff
LOGICAL :: f_daily_sediment_production,f_daily_subsurface_runoff,f_daily_theta,f_daily_total_overlandflow,f_daily_water_subbasin,f_routing_response,f_sediment_production,f_water_subbasin,f_deep_gw_recharge,f_deep_gw_discharge,f_daily_gw_loss,f_tc_theta
LOGICAL :: f_tc_surfflow, f_tc_sedout
LOGICAL :: f_river_degradation, f_river_deposition, f_river_flow, f_river_flow_dailyaverage, f_river_flowdepth, f_river_sediment_concentration, f_river_sediment_total, f_river_sediment_total_dailyaverage, f_river_storage, f_river_sediment_storage, f_river_velocity, f_river_bedload		 
LOGICAL :: f_actetranspiration,f_qhorton,f_subsurface_runoff,f_total_overlandflow,f_gw_discharge,f_potetranspiration,f_gw_loss,f_gw_recharge
LOGICAL :: f_res_watbal,f_res_vollost,f_res_cav,f_res_hydraul,f_res_bedchange,f_res_sedbal,f_res_longitudunal,f_res_sedcomposition
LOGICAL :: f_lake_inflow_r,f_lake_outflow_r,f_lake_retention_r,f_lake_volume_r,f_lake_sedinflow_r,f_lake_sedoutflow_r,f_lake_sedretention_r,f_lake_sedimentation_r
LOGICAL :: f_lake_watbal,f_lake_sedbal,f_lake_inflow,f_lake_outflow,f_lake_volume,f_lake_retention,f_lake_vollost,f_lake_sedinflow,f_lake_sedoutflow,f_lake_sizedistoutflow

LOGICAL :: doloadstate=.TRUE.			!load initial values before model run
LOGICAL :: dosavestate=.TRUE.			!save state of model after execution
REAL	:: default_rel_sat=1.0		!default relative saturation (initial value for all horizons, if not specified otherwise)
REAL	:: default_gw_storage=0	!default ground water storage (to be assumed for all non-specified LUs) [mm]



logical,parameter :: domean=.FALSE.	! calculate daily ETP from mean daily values (this flag is dominant over donight)
logical,parameter :: donight=.TRUE.		!if (domean==0):calculate daily ETP as sum of ETP for day and night (donight=1) or only day
real :: daily_delta_temp !daily temperature amplitude (Tmin=Tmean-daily_delta_temp; Tmax=Tmean+daily_delta_temp;)
REAL,parameter :: hours_of_daylight=12 !number of hours with sunlight per day

end module common_h



!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
module params_h
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
save
! TOTAL NUMBER OF SUB-BASINS IN STUDY AREA (e.g. 10 sub-basins)
INTEGER :: subasin
! TOTAL NUMBER OF SOTER UNITS IN STUDY AREA (e.g. 321 units)
INTEGER :: nsoter
! TOTAL NUMBER OF TERRAIN COMPONENTS IN STUDY AREA (e.g. 515 components)
INTEGER :: nterrain
! TOTAL NUMBER OF SOIL COMPONENTS IN STUDY AREA (e.g. 72 components)
INTEGER :: nsoil
! TOTAL NUMBER OF VEGETATION UNITS IN STUDY AREA (e.g. 34 units)
INTEGER :: nveg
! TOTAL NUMBER OF CELL/SOTER UNIT/TERRAIN COMPONENT COMBINATIONS (e.g. 49 combinations)
INTEGER :: sv_comb
! MAXIMUM NUMBER OF SOTER UNITS IN CELLS (7)
INTEGER :: maxsoter
!PARAMETER (maxsoter=7)
! MAXIMUM NUMBER OF TERRAIN COMPONENTS IN SOTER UNIT (3)
INTEGER :: maxterrain
!PARAMETER (maxterrain=3)
! MAXIMUM NUMBER OF SOIL COMPONENTS IN TERRAIN COMPONENTS (28)
INTEGER :: maxsoil
!PARAMETER (maxsoil=28)
! MAXIMUM NUMBER OF HORIZONS IN SOIL COMPONENTS (8)
INTEGER :: maxhori
!PARAMETER (maxhori=8)
! TOTAL NUMBER OF TRANSPOSITIONS BETWEEN Sub-basins
INTEGER :: ntrans
!PARAMETER (ntrans=2)
!common / basin_parameter / subasin, sv_comb, nsoter, nterrain, nsoil, nveg
end module params_h



!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
module time_h
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
save
! time of start/end of simulation
INTEGER :: tstart, tstop, mstart, mstop
! year
INTEGER :: t
! month
INTEGER :: m
! day (in year)
INTEGER :: d, dprev
! day since simulation start
INTEGER :: dtot
! days in year
INTEGER :: dayyear, daylastyear, dayoutsim
! days in months
INTEGER :: daymon(12)

! number of normal years (no Schaltjahr) in calculation period
INTEGER*1 nos
! number of years in calculation period
INTEGER*1 years
INTEGER :: idold,idmon,mnew
INTEGER :: mon_day(12)
!COMMON /maintime/ tstart, tstop, mstart, mstop, t, m, d,  &
!    dprev, dtot, dayyear, daylastyear, idold, idmon, mnew, nos, years
!COMMON /timedat/ daymon,daynr,monnr,mon_day
end module time_h



!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
module climo_h
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
save
! daily average temperatures (oC)
 !Allocatable      real temp(366,nmun)
real, allocatable :: temp(:,:)
! daily precipitation (mm/day)
!Allocatable      real precip(366,nmun)
real, allocatable :: precip(:,:)
! daily potential evapotranspiration (mm/day)
!Allocatable      real pet(366,nmun)
real, allocatable :: pet(:,:)
! daily mean shortwave radiation (W/m^2)
!Allocatable      real rad(366,nmun)
real, allocatable :: rad(:,:)
! incoming radiation above atmosphere (W/m²)
REAL :: radex(366)
! relative humidity (%)
!Allocatable      real rhum(366,nmun)
real, allocatable :: rhum(:,:)
! wind velocity (m/s)
!Allocatable      real wind(366,nmun)
real, allocatable :: wind(:,:)
! climate variables for hourly model version
!Allocatable      real preciph(366*24)
real, allocatable :: preciph(:,:)


end module climo_h

