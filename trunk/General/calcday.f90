! Till: computationally irrelevant: minor changes to improve compiler compatibility
! 2011-04-29

SUBROUTINE calcday

use time_h

INTEGER :: temp
! day in month given day in year (1=normal, 2=leap year)

! month given day in year (1=normal, 2=leap year)
INTEGER(kind=1) :: monnr(367,2)

DATA monnr /31* 1,28* 2,31* 3,30* 4,31* 5,30* 6,  &
    31* 7,31* 8,30* 9,31*10,30*11,31*12,2* 1,  &
    31* 1,29* 2,31* 3,30* 4,31* 5,30* 6, 31* 7,31* 8,30* 9,31*10,30*11,31*12,1* 1/


IF ( is_leapyear(t) ) THEN
  daymon(2) =  29
  temp=2
ELSE
  daymon(2) =  28
  temp=1
END IF

!idold = d
IF (t == tstart) THEN
!  idmon = daynr(d+sum(daymon(1:mstart-1)),temp)
!  m     = monnr(d+sum(daymon(1:mstart-1)),temp)
!  mnew  = monnr(d+1+sum(daymon(1:mstart-1)),temp)
ELSE
!  idmon = daynr(d,temp)
!  m     = monnr(d,temp)
!  mnew  = monnr(d,temp)
END IF

RETURN
END SUBROUTINE calcday
