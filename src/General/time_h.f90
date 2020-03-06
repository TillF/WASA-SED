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


