SUBROUTINE routing(STATUS)
 
use lake_h
use climo_h
use common_h
use hymo_h
use params_h
use routing_h
use time_h
use reservoir_h

IMPLICIT NONE

INTEGER, INTENT(IN OUT)                  :: STATUS


! status of call (0=initialization run, 1=initialization year,
!                 2=calculation day,    3=finalization year)

INTEGER :: bat
!    id           : additional loop variable of days (total 7 days)
INTEGER :: imun,imunx,irout,irout2,irout_d,id,idummy,imeso,istate
INTEGER :: upstream, downstream
INTEGER :: itl, itr, ih, mm, imunout, iout, i, j, make
REAL :: xdum(48),check,temp2,qtemp, storcapact
character(len=1000) :: fmtstr	!string for formatting file output


! -----------------------------------------------------------------------
IF (STATUS == 0) THEN
!**  Read routing paramter
OPEN(11,FILE=pfadp(1:pfadj)// 'River/routing.dat',STATUS='old')
! upbasin: MAP ID of upstream sub-basin (MAP IDs)
! downbasin: MAP ID of downstream sub-basin (MAP IDs)
READ (11,*); READ(11,*)
DO irout=1,subasin
  READ (11,*)  idummy, upbasin(irout),downbasin(irout)
END DO
CLOSE (11)

!**  Read hydrological response and reservoir paramter
OPEN(11,FILE=pfadp(1:pfadj)// 'River/response.dat'  &
    ,STATUS='old')
! prout(i,1): lag time [d]
! prout(i,2): retention storage [d]
READ (11,*); READ(11,*)
DO i=1,subasin
  READ (11,*)  idummy,prout(i,1),  prout(i,2)
  IF (idummy /= id_subbas_extern(i)) THEN
    WRITE(*,*) 'Sub-basin-IDs in file response.dat  &
        must have the same ordering scheme as in hymo.dat'
    STOP
  END IF
END DO
CLOSE (11)

!this relates the MAP IDs (id_subbas_extern(subasin)) to the sorted CODE IDs (id_subbas_intern(subasin)),
! i.e. upbasin and downbasin are now numbered according to internal ids
!     so that in routing.dat only the MAP IDs have to be read in,
! first for the ID of upstream TEZG
DO i=1,subasin
  j=1
  DO WHILE (id_subbas_extern(j) /= upbasin(i))
    j=j+1
    IF (j > 500) THEN
      WRITE (*,*) 'upbasin(i) loop in readhymo.f'
      STOP
    END IF
  END DO
  upbasin(i)=j
END DO
! second for the ID of downstream TEZG
DO i=1,subasin
  IF (downbasin(i) /= 999.AND.downbasin(i) /= 9999) THEN
    j=1
    DO WHILE (id_subbas_extern(j) /= downbasin(i))
      j=j+1
      IF (j > 500) THEN
        WRITE (*,*) 'downsbasin(i) loop in readhymo.f'
        STOP
      END IF
    END DO
    downbasin(i)=j
  END IF
END DO

! initialize
  DO i=1,subasin
    DO id=1,dayyear+6
      qout  (id,i) = 0.
      qin   (id,i) = 0.
    END DO
  END DO
 
 
! INITIALISATION OF OUTPUT FILES

OPEN(11,FILE=pfadn(1:pfadi)//'River_Flow.out',STATUS='replace')
  if (f_river_flow) then
    WRITE (11,*) 'Output files for river discharge q_out (m3/s) (with MAP IDs as in hymo.dat)'
	write(fmtstr,'(a,i0,a)')'(2a6,',subasin,'i14)'		!generate format string
	WRITE (11,fmtstr)' year ', ' day  ', (id_subbas_extern(i), i=1,subasin)
	!WRITE (11,'(2a6,<subasin>i14)')' year ', ' day  ', (id_subbas_extern(i), i=1,subasin)
    Close (11)
  else
    close(11, status='delete') !delete any existing file, if no output is desired
  endif


!  OPEN(11,FILE=pfadn(1:pfadi)//'daily_qin_m3s.out',STATUS='replace')
!  WRITE (11,*) 'Output files for qin water_subbasin into the sub-basins  &
!      [m**3/s](with MAP IDs as in hymo.dat)'
!  WRITE (11,'(2a4,<subasin>i14)')'Year' ,' Day', (id_subbas_extern(i), i=1,subasin)
!  CLOSE (11)
  
! calculate routing response function for each sub-basin
! (given partaneter lag-time tL and retention-time tR)
  OPEN(11,FILE=pfadn(1:pfadi)//'routing_response.out' ,STATUS='replace')
  WRITE(11,'(a)') 'Output of linear response function'
  WRITE(11,'(a)')'Subasin-ID,translation [days],  &
      retention [days], uh(1,7) [-]'
  
  DO i=1,subasin
    itl = nint (prout(i,1))
    itr = nint (prout(i,2))
    DO ih =1,7
      IF (ih <= itl) THEN
        hrout(ih,i)=0.
      ELSE IF (ih > itl+itr) THEN
        hrout(ih,i)=0.
      ELSE
! Calculation of the linear response function for runoff routing in the river network
        hrout(ih,i)=1./REAL(itr*itr)* (2.*REAL(itr)-2.*REAL(ih-itl)+1.)
      END IF
    END DO
    check=0.
    DO ih =1,7
      check=check+hrout(ih,i)
    END DO
    IF (check < 0.1) THEN
      hrout(1,i) = 1.0
    END IF
    
    WRITE (11,'(I5,2f10.2,2x,7f5.2)')  &
        id_subbas_extern(i),prout(i,1),prout(i,2),(hrout(ih,i),ih=1,7)
  END DO
  
! CALL Reservoir Sedimentation and Management Modules
  IF (doreservoir) THEN
    CALL reservoir (0,upstream,idummy)
  END IF
  
END IF

! ------------------------------------------------------------------------
IF (STATUS == 1) THEN
  
! initialize ...
! ... and take qout and volact from the last 6 days of last year
  IF (t > tstart) THEN
    DO i=1,subasin
      DO id=1,6
        qout(id,i)  =qout(daylastyear+id,i)
      END DO
    END DO
    DO i=1,subasin
      DO id=7,dayyear
        qout(id,i)=0.
      END DO
    END DO
  ELSE IF (t == tstart) THEN
    DO i=1,subasin
      qout(1,i)=0.
    END DO
  END IF
  
! ... and initialize remaining values
  DO i=1,subasin
    DO id=dayyear+1, dayyear+6
      qout(id,i)=0.
    END DO
  END DO
  
! CALL Reservoir Sedimentation and Management Modules
  IF (doreservoir) THEN
    CALL reservoir (1, upstream,idummy)
  END IF
  
END IF

! ------------------------------------------------------------------------
IF (STATUS == 2) THEN
  
! ..........................................................................
!**  Transfer of water between sub-basins
!    assumption: time delay = 1 day
!    transfer variable qtemp in [m3/day]
  
  IF (dotrans) THEN
! ntrans=4 defined in params.fi
    DO i=1,ntrans
      IF (t >= y_trans(i)) THEN
!    from river
        qtemp=0.
        irout=trans_start(1,i)
        IF (trans_start(2,i) == 2) THEN
          qtemp=MIN(qout(d,irout)*86400.,q_trans(i)*86400.)
          qout(d,irout)=qout(d,irout)-qtemp/86400.
          qout(d,irout)=MAX(qout(d,irout),0.)
!    from acude (there must be a large acude in this muni !)
        ELSE IF (trans_start(2,i) == 1) THEN
          qtemp=MIN(volact(d,irout)*1.e6,q_trans(i)*86400.)
          volact(d,irout)=volact(d,irout)-qtemp/1.e6
          volact(d,irout)=MAX(volact(d,irout),0.)
        END IF
!    into river
        irout=trans_end(1,i)
        IF (trans_end(2,i) == 2) THEN
          qout(d+1,irout)=qout(d+1,irout)+qtemp*(1.-loss_trans(i))/ 86400.
!    into acude (there must be a large acude in this muni !)
        ELSE IF (trans_end(2,i) == 1) THEN
          volact(d+1,irout)=volact(d+1,irout)+ qtemp*(1.-loss_trans(i))/1.e6
        END IF
      END IF
    END DO
  END IF
  
!  water_subbasin(d,i): runoff of each sub-basin (after small reservoirs)
!  (in m**3/d), assumpiton: leaving sub-basin with time delay = 1 day
!Latest version: water_subbasin already comes in m3/s from hymo_all.f
  DO i=1,subasin
!    water_subbasin(d,i)=water_subbasin(d,i)/(3600*24)
    qout(d+1,i)=water_subbasin(d,i)+qout(d+1,i)
  END DO
! set inflow qin = zero for all sub-basin
  DO i=1,subasin
    DO ih=1,7
      qin(d+ih-1,i)=0.
    END DO
  END DO
  
!cccccccccccccccccccccccccccc
! MAIN ROUTING LOOP
!cccccccccccccccccccccccccccc
! Calculate qin and qout in Fliessbaum order (as was read in and transformed in routing.dat)
  DO i=1,subasin
    upstream=upbasin(i)  !internal code-ID for most upstream sub-basin
    downstream=downbasin(i) !internal code-ID for receiving sub-basin
    
! Route inflow from upstream sub-basins (qin) through actual sub-basin within 7 days
    DO ih=1,7
      qout(d+ih-1,upstream)=qin(d,upstream)*hrout(ih,upstream)  &
          + qout(d+ih-1,upstream)
! Transmission losses by evaporation in river
! river width for given discharge (estimated according to global
! relationship given by Leopold, 1994 (Fig 8.10) [m]
      IF (ih == 1) THEN
        IF (qout(d+ih-1,upstream) > 0.) THEN
          temp2=10.**(LOG10(qout(d+ih-1,upstream))*0.494+1.031)
        ELSE
          temp2=0.
        END IF
! Length of river stretch in sub-basin estimated as diametre of (circular)
! catchment area multiplied by sinuoisity factor of meandering stream (=1.5) [m]
        temp2=temp2*2.*((area(upstream)/3.147)**0.5)*1000.*1.5
        temp2=temp2*pet(d,upstream)/1000.
        temp2=MIN(temp2,qout(d+ih-1,upstream)*86400.)
        qloss(d,upstream)=qloss(d,upstream)+temp2
        qout(d+ih-1,upstream)=qout(d+ih-1,upstream)-temp2/86400.
        qout(d+ih-1,upstream)=MAX(0.,qout(d+ih-1,upstream))
      END IF
    END DO
! end of do-loop for summying up flows within seven days
    
! Assign outflow as inflow into the next downstream sub-basin
! and add possible inflow from other upstream sub-basins
    IF (doreservoir) THEN
      CALL reservoir (2,upstream,idummy)
    END IF

	IF (downstream /= 9999 .AND. downstream /= 999) THEN
      qin(d,downstream)=qin(d,downstream) + qout(d,upstream)
    END IF

  END DO
!  END OF MAIN ROUTING LOOP
  
  
  

END IF


! -----------------------------------------------------------------------
IF (STATUS == 3) THEN
  

! daily output of water discharge in the river for the entire year
if (f_river_flow) then
  OPEN(11,FILE=pfadn(1:pfadi)//'River_Flow.out',STATUS='old' ,POSITION='append'  )
  write(fmtstr,'(a,i0,a)')'(2i6,',subasin,'f14.3)'		!generate format string
  do j=1, dayyear
	WRITE (11,fmtstr)t, j, (qout(j,i),i=1,subasin)
	!WRITE (11,'(2i6,<subasin>f14.3)')t, j, (qout(j,i),i=1,subasin)
  enddo
  CLOSE (11)
endif
  
  IF (doreservoir) THEN
    CALL  reservoir (3,idummy,idummy)
  END IF
END IF

RETURN
END SUBROUTINE routing
