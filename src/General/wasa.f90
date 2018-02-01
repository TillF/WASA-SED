PROGRAM wasa_sed

use common_h
use hymo_h
use params_h
use time_h
use routing_h
use model_state_io

IMPLICIT NONE

!CCCCCCCCCCCCCCCCCCCCC MAIN PROGRAM CCCCCCCCCCCCCCCCCCCCCCCCCCCC

CALL GETARG(1, pfadn)		!Till: try to read path to central control file (do.dat)

include 'svn_rev.f90' !Till: import revision information string

WRITE(*,'(A,A,/,A)') 'WASA model, ',trim(rev_string1),' '//trim(rev_string2)

WRITE(*,*) ':Initialization'
CALL readgen(pfadn)
t = tstart
CALL calcyear

! READ INPUT DATA (STATUS 0) AND INITIALISE ARRAYS (STATUS 0)
CALL hymo_all(0)
CALL climo(0)


!initialization of routing routines
if (river_transport.eq.1) CALL routing(0)
if (river_transport.ne.1) CALL routing_new(0)
call init_river_state  

CALL save_model_state(doloadstate, .TRUE.) !Till: do backups of state files if loaded from them, and save only summary on initial storage

!!! MAIN LOOP FOR EACH YEAR
DO t=tstart, tstop
  WRITE(*,*) tstart,tstop
  IF (t /= tstart) CALL calcyear
  WRITE(*,*) 'calculations for year ',t

! Call routines for climate time series for current year
  CALL climo(1)

! Update annual values for hydrology and agriculture at start of each simulation year
     call hymo_all(1)
	 if (river_transport.eq.1) call routing(1)
	 if (river_transport.ne.1) call routing_new(1)

!!! MAIN LOOP FOR DAILY TIME STEPS
  DO d=1,dayyear
  write(*,*) d
    CALL hymo_all(2)
 !Call routing routine
	if (river_transport.eq.1) CALL routing(2)
    if (river_transport.ne.1) CALL routing_new(2)
    dprev = d
  END DO
!!! END LOOP FOR DAILY TIME STEPS

! Generate output for hillslope
   CALL hymo_all(3)
!Generate output for river and reservoir
  if (river_transport.eq.1) CALL routing(3)
  if (river_transport.ne.1) CALL routing_new(3)

END DO


! Close climate input files
CLOSE(81)
CLOSE(82)
CLOSE(83)
CLOSE(84)

END PROGRAM wasa_sed
