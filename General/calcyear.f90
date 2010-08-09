SUBROUTINE calcyear
 
! Code converted using TO_F90 by Alan Miller
! Date: 2005-06-30  Time: 13:50:58
use time_h

INTEGER :: i


daymon( 1) = 31
daymon( 3) = 31
daymon( 4) = 30
daymon( 5) = 31
daymon( 6) = 30
daymon( 7) = 31
daymon( 8) = 31
daymon( 9) = 30
daymon(10) = 31
daymon(11) = 30
daymon(12) = 31

IF (t == tstart) THEN
! Checks for leap year and number of days in February for first simulation year
  IF (MOD(t,4) == 0) THEN
    daymon(2) =  29
    dprev     = sum(daymon(mstart:12))
    dayyear   = sum(daymon(mstart:12))
  ELSE
    daymon(2) =  28
    dprev     = sum(daymon(mstart:12))
    dayyear   = sum(daymon(mstart:12))
  END IF
!Calculation of total days since simulation started (for first year)
  dtot=dayyear
END IF

!Ge check total number of days before the start month in the start year
IF (t == tstart .and. mstart /= 1) THEN
  dayoutsim = sum(daymon(1:mstart-1))
ELSE
  dayoutsim = 0
ENDIF 

! Checks for leap year and number of days in February for all other simulation years
IF (t /= tstart.AND.t /= tstop) THEN
  IF (MOD(t,4) == 0) THEN
    daymon(2) =  29
    dayyear   = sum(daymon(1:12))
  ELSE
    daymon(2) =  28
    dayyear   = sum(daymon(1:12))
    nos=nos+1
  END IF
  dtot=dtot+dayyear
END IF

! Calculates no. of days of previous year (as a function of leap years)
IF (MOD(t-1,4) == 0.AND.t-1 /= tstart) THEN
  daylastyear = 366
ELSE IF ((MOD(t-1,4) /= 0).AND.t-1 /= tstart) THEN
  daylastyear = 365
END IF
IF (t == tstart) THEN
  daylastyear = sum(daymon(mstart:12))
END IF

! Calculates leap year and no. of days for the last simulation year
IF (t == tstop.AND.t /= tstart) THEN
  IF (MOD(t,4) == 0) THEN
    daymon(2) =  29
    dayyear   = sum(daymon(1:mstop))
  ELSE
    daymon(2) =  28
    dayyear   = sum(daymon(1:mstop))
  END IF
  dtot=dtot+dayyear
END IF

IF (t == tstop .AND. t == tstart) THEN
  IF (MOD(t,4) == 0) THEN
    daymon(2) =  29
    dayyear   = sum(daymon(mstart:mstop))
  ELSE
    daymon(2) =  28
    dayyear   = sum(daymon(mstart:mstop))
  END IF
  dtot=dayyear
END IF

mon_day(1)= 1
DO i=2,12
  mon_day(i)=mon_day(i-1)+daymon(i-1)
END DO


m    = 1
mnew = 1
d     = 1
idold = 1
idmon = 1
years=tstop-tstart+1

RETURN
END SUBROUTINE calcyear
