! Till: separated modules to allow easier use by eclipse
! 2013-10-02 

! Till: computationally irrelevant: relocated vars
! 2012-09-21 

! Till: computationally irrelevant: minor changes to improve compiler compatibility
! 2011-04-29

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

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
module climo_h
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
save
! daily average temperatures (oC)
real, allocatable :: temp(:,:)
! daily precipitation (mm/day)
real, allocatable :: precip(:,:)
! daily potential evapotranspiration (mm/day)
real, allocatable :: pet(:,:)
! daily mean shortwave radiation (W/m^2)
real, allocatable :: rad(:,:)
! incoming radiation above atmosphere (W/m²)
REAL :: radex(366)
! relative humidity (%)
real, allocatable :: rhum(:,:)
! wind velocity (m/s)
real, allocatable :: wind(:,:)
! climate variables for hourly model version
real, allocatable :: preciph(:,:)

INTEGER :: store_day, store_timestep	!Till: day and timestep for which ETP-values were stored

end module climo_h

