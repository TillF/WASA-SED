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

dayoutsim = 0

IF (is_leapyear(t)) THEN !check for leap year
	daymon(2) =  29
ELSE
	daymon(2) =  28
END IF

mon_day(1)= 1 !compute julian day of start of each month
DO i=2,12
  mon_day(i)=mon_day(i-1)+daymon(i-1)
END DO



! Calculates no. of days of simulated days in previous year (as a function of leap years)
IF (t == tstart) THEN
  daylastyear = sum(daymon(mstart:12))
else
    IF (t-1 /= tstart) THEN
        IF ( is_leapyear(t-1) ) THEN
          daylastyear = 366
        ELSE 
          daylastyear = 365
        END IF 
    END IF      
END IF




! first simulation year
IF (t == tstart) THEN	
  dstart=min(dstart, daymon(mstart)) !correct start/stop days, when larger than days in month
  dstop =min(dstop,  daymon(mstop))

  dprev     = sum(daymon(mstart:12)) !ii: is this correct?
  dayyear   = sum(daymon(mstart:12)) !number of days to treat in current year

  dayyear	  =	dayyear   - (dstart-1) !reduce by start_day_of month
  
  dtot=dayyear !Calculation of total days since simulation started (for first year)

	!Ge check total number of days before the start month in the start year
	IF (mstart /= 1) THEN
	  dayoutsim = sum(daymon(1:mstart-1))
	ENDIF 
	
	dayoutsim = dayoutsim + (dstart-1) 

END IF


! intermediate simulation years
IF (t /= tstart .AND. t /= tstop) THEN
  dayyear   = sum(daymon(1:12))  
  dtot=dtot+dayyear !count simulation days
END IF


! last simulation year
IF (t == tstop .AND. t /= tstart) THEN
  dayyear   = sum(daymon(1:mstop))
  dayyear	  =	dayyear   - (daymon(mstop)- dstop) !reduce by end_day_of month

  dtot=dtot+dayyear !count simulation days
END IF


!first=last simulation year
IF (t == tstop .AND. t == tstart) THEN
  dayyear   = sum(daymon(mstart:mstop))
  dayyear	= dayyear    - (dstart-1) - (daymon(mstop)- dstop) !reduce by start_day_of month and end_day_of month
  dtot=dayyear
END IF



d     = 1 !set counter for days in simulation year
years=tstop-tstart+1


!dprev?		simulation day number of previous day (related to d)


RETURN
END SUBROUTINE calcyear
