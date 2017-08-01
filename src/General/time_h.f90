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

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
module time_h
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
save
! time of start/end of simulation
INTEGER :: tstart, tstop, mstart, mstop, dstart, dstop
! year
INTEGER :: t
! month
INTEGER :: mm
! day (in year)
INTEGER :: d, dprev
! day since simulation start
INTEGER :: dtot
! days in year
INTEGER :: dayyear, daylastyear, dayoutsim
! days in months
INTEGER :: daymon(12)

! number of years in calculation period
INTEGER :: years
!INTEGER :: idold,idmon,mnew
INTEGER :: mon_day(12)

contains
logical FUNCTION is_leapyear(yr) !check if this is a leap year
integer :: yr
 is_leapyear = (MOD(yr,4) == 0)  .AND. ( (MOD(yr,100) /= 0) .OR. (MOD(yr,400) == 0) )
 return
END FUNCTION is_leapyear

end module time_h


