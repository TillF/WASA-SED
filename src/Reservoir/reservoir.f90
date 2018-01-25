SUBROUTINE reservoir (STATUS,upstream,res_h)

! Till: made reading of reservoir.dat, cav.dat independent of order of IDs
! memory allocation according to content of these files
! various minor code beautifications
! 2012-09-17

! Till: computationally irrelevant: outcommented unused vars
! 2012-09-14

! Till: computationally irrelevant: streamlined code and improved error message with routing.dat
! 2011-07-05

! Till: computationally irrelevant: minor changes to improve compiler compatibility
! 2011-04-29

! Code converted using TO_F90 by Alan Miller
! Date: 2005-08-23  Time: 12:56:42

use common_h
use params_h
use routing_h
use climo_h
use hymo_h
use time_h
use reservoir_h
use utils_h

! water balance for a reservoir of a large river dam

IMPLICIT NONE


INTEGER, INTENT(IN)                  :: STATUS
INTEGER, INTENT(IN OUT)                  :: upstream
INTEGER, INTENT(IN)                  :: res_h

! STATUS of CALL (0=initialization run, 1=initialization year,
!                 2=calculation day,    3=finalization year)

INTEGER :: i,id,dummy1,dummy2,cont,s,p,q,h,istate,ka,ih !,irout,imun,dummy1a,dummy2a,f
!Ge include j,nbrbat1,cont,upstream,downstream
INTEGER :: j,nbrbat1,flag_cav
REAL :: elevhelp,evaphelp2,volhelp !,elevhelp2,elevhelp3
REAL :: help,help1,help2,help3,evaphelp,areahelp,infhelp,helpout,prechelp !,helpin
!Ge actual storage capacity of large river reservoir in a certain year[10**6 m**3]
!REAL :: storcapact

REAL :: res_flag(subasin)
CHARACTER(12) :: subarea

REAL :: r_level0,r_level1,r_overflow,r_qintake,r_qbottom
REAL :: r_precip,r_etp,r_qinflow

INTEGER :: idummy,n !,npt, nbrsec1
!INTEGER :: dummy3
REAL :: dummy4,dummy5,dummy6
character(len=1000) :: fmtstr	!string for formatting file output

!*****************************************************************************
!temporary variables used to test the cascade routing scheme of the lake module
!INTEGER :: a1,b1,c1,d1,e1,k1
!REAL :: inflow_class(6),outflow_class(6),retention_class(6),volume_class(6)
!REAL :: sedinflow_class(6),sedoutflow_class(6),sedretention_class(6),sedimentation_class(6)
! -----------------------------------------------------------------------
IF (STATUS == 0) THEN

reservoir_check=0 !(0=simulation will all components; 1=simulation without hillslope and river modules)
reservoir_balance=1 !(0=inflow and outflow discharges must be provided as input file; 1=only inflow discharges must be provided as input file)
reservoir_print=0 !(0=results printed at the end of the timestep; 1=results printed at the end of the simulated year)

if (reservoir_check==0) reservoir_balance=1

! Initialize
  DO i=1,subasin
!Ge "volact" vector contains now only 365/366 values per annum
!Ge others reservoir parameters have to be inserted
    DO id=1,dayyear*nt
      volact(id,i) = 0. !actual stored volume in reservoir [10**6 m**3]
      qinflow(id,i) = 0.
      qintake(id,i) = 0.
      overflow(id,i) = 0.
    END DO
  END DO

  IF (reservoir_check == 1) THEN
    OPEN(11,FILE=pfadp(1:pfadj)// 'Hillslope/hymo.dat',STATUS='old')
    READ (11,*);READ (11,*)

    DO i=1,subasin
      READ(11,*) id_subbas_extern(i)
      id_subbas_intern(i)=id_subbas_extern(i)
    END DO
    DO i=1,subasin
      n=1
      DO WHILE (id_subbas_extern(n) /= id_subbas_intern(i))	!search until the current external ID has been found...
        n=n+1
      END DO
!Eva id_subbas_intern is sorted (1 - no. of subasin), id_subbas_extern contains the original Map IDs
      id_subbas_intern(i)=n
!write(*,*)id_subbas_intern(i),id_subbas_extern(i)
    END DO
    CLOSE (11)

! Read configuration of river system (in routing order)
	OPEN(11,FILE=pfadp(1:pfadj)// 'River/routing.dat',STATUS='old')! upbasin: MAP ID of upstream sub-basin (MAP IDs);! downbasin: MAP ID of downstream sub-basin (MAP IDs)
	READ (11,*); READ(11,*)
	DO i=1,subasin
	  READ (11,*)  idummy, upbasin(i),downbasin(i)
	END DO
	CLOSE (11)

!this relates the MAP IDs (id_subbas_extern(subasin)) to the sorted CODE IDs (id_subbas_intern(subasin)), i.e. upbasin and downbasin are now numbered according to internal ids
!so that in routing.dat only the MAP IDs have to be read in, first for the ID of upstream subasin
	!replace external with internal IDs
	DO i=1,subasin
	  !upstream basin referencing
	  ih=which1(id_subbas_extern == upbasin(i))

	  IF (ih==0) THEN
		  WRITE (*,'(A,I0,A)') 'unknown upstream subbasin ID ', upbasin(i),' in routing.dat'
		  STOP
	  else
		upbasin(i)=ih
	  END IF

	  !downstream basin referencing
	  IF (downbasin(i) == 999 .OR. downbasin(i) == 9999) cycle 	!999 and 9999 mark outlet
	  ih=which1(id_subbas_extern == downbasin(i))

	  IF (ih==0) THEN
		  WRITE (*,'(A,I0,A)') 'unknown downstream subbasin ID ', downbasin(i),' in routing.dat'
		  STOP
	  else
		downbasin(i)=ih
	  END IF
	END DO
  ENDIF

!! Read reservoir parameters
res_flag(:)=-999.
storcap(:)=0.

  OPEN(11,FILE=pfadp(1:pfadj)// 'Reservoir/reservoir.dat',IOSTAT=istate,STATUS='old')
	IF (istate/=0) THEN					!reservoir.dat not found
	  write(*,*)pfadp(1:pfadj)// 'Reservoir/reservoir.dat was not found, please provide it.'
	  stop
	ENDIF

!    ELSE
!      READ(11,*)
!	  READ(11,*)
!      READ(11,*,IOSTAT=ka) dummy1,dummy2
!      DO i=1,subasin
!        IF (dummy1==id_subbas_extern(i)) THEN
!          res_flag(i)=dummy2
!	    ENDIF
!	  ENDDO
!      DO WHILE (ka==0)
!	    READ(11,*, IOSTAT=ka) dummy1,dummy2					!read next line in file
!        DO i=1,subasin
!          IF (dummy1==id_subbas_extern(i)) THEN
!            res_flag(i)=dummy2
!		  ENDIF
!	    ENDDO
!      END DO
!	ENDIF
!  CLOSE (11)

!DO i=1,subasin
!write(*,*)id_subbas_extern(i),res_flag(i)
!enddo
!stop

  READ(11,*)
  READ(11,*)
  istate=0
  j=2 !for counting lines

  DO WHILE (.TRUE.)
	  READ (11,'(A)', IOSTAT=istate)fmtstr
	  if (istate/=0) exit !exit loop when no more line encountered
	  j=j+1
	  READ (fmtstr,*, IOSTAT=istate)dummy1

      i=which1(id_subbas_extern == dummy1)

	  IF (i==0) THEN
		  WRITE (*,'(A,I0,A)') 'WARNING: unknown upstream subbasin ID ', dummy1,' in reservoir.dat, ignored.'
		  cycle
	  END IF

	  READ (fmtstr,*,IOSTAT=istate) dummy1, minlevel(i), maxlevel(i),vol0(i),storcap(i), &
			damflow(i),damq_frac(i),withdrawal(i),damyear(i),maxdamarea(i), &
			damdead(i),damalert(i),dama(i),damb(i),qoutlet(i),fvol_bottom(i), &
			fvol_over(i),damc(i),damd(i),elevbottom(i)

	  IF (istate/=0) THEN
		  WRITE (*,'(A,i0,A,A)') 'ERROR: Format error in reservoir.dat, line',j,':', fmtstr
		  STOP
	  END IF

!Ge "storcap" and "vol0" are read in 1000m**3 and after that they are converted into 10**6 m**3  &
!Ge damarea renamed to maxdamarea  &
!Ge "damdead" is read in 1000m**3 and after that it is converted into 10**6 m**3  &

      forma_factor(i)=1.e3*storcap(i)/((maxlevel(i)-minlevel(i))**3)
      IF (vol0(i) /= -999.) THEN
        vol0(i)=vol0(i)/1.e3
      END IF
      storcap(i)=storcap(i)/1.e3
      damdead(i)=damdead(i)/1.e3
      damalert(i)=damalert(i)/1.e3
  END DO
  CLOSE (11)

! Check lateral inflow directly into the subbasins' reservoir
  latflow_res(1:subasin)=0

  OPEN(11,FILE=pfadp(1:pfadj)// 'Reservoir/lateral_inflow.dat', IOSTAT=istate,STATUS='old')
	IF (istate/=0) THEN					!lateral_inflow.dat not found
	  write(*,*)pfadp(1:pfadj)// 'Reservoir/lateral_inflow.dat was not found. Run the model anyway.'
      DO i=1,subasin
        latflow_res(i)=0
      ENDDO
    ELSE
	  READ(11,*)
	  READ(11,*)
      READ(11,*,IOSTAT=ka) dummy1
      DO i=1,subasin
        IF (dummy1==id_subbas_extern(i)) THEN
          latflow_res(i)=1
	    ENDIF
	  ENDDO
      DO WHILE (ka==0)
	    READ(11,*, IOSTAT=ka) dummy1					!read next line in file
        DO i=1,subasin
          IF (dummy1==id_subbas_extern(i)) THEN
            latflow_res(i)=1
		  ENDIF
	    ENDDO
      END DO
	ENDIF
  CLOSE (11)

  OPEN(11,FILE=pfadp(1:pfadj)// 'Reservoir/lateral_inflow.dat', IOSTAT=istate,STATUS='old')
	IF (istate==0) THEN
      READ(11,*)
	  READ(11,*)
      DO i=1,subasin
        IF (latflow_res(i)==1) THEN
          READ (11,*)dummy1,reservoir_down(i)

          IF (dummy1 /= id_subbas_extern(i)) THEN
            WRITE(*,*) 'Sub-basin-IDs in file lateral_inflow.dat must have the same ordering scheme as in hymo.dat'
            STOP
          END IF

          IF (reservoir_down(i) /= 999.AND.reservoir_down(i) /= 9999) THEN
            j=1
            DO WHILE (id_subbas_extern(j) /= reservoir_down(i))
              j=j+1
              IF (j > 500) THEN
                WRITE (*,*) 'downsbasin(i) loop in readhymo.f'
                STOP
              END IF
            END DO
            reservoir_down(i)=j
          END IF
	    ELSE
		  dummy1=id_subbas_extern(i)
		ENDIF
        IF (dummy1 /= id_subbas_extern(i)) THEN
         WRITE(*,*) 'Sub-basin-IDs in file operat_rule.dat must have the same ordering scheme as in hymo.dat'
         STOP
        END IF
     ENDDO
	ENDIF
  CLOSE (11)

  OPEN(11,FILE=pfadp(1:pfadj)// 'Reservoir/operat_rule.dat', IOSTAT=istate,STATUS='old')
  IF (istate/=0) THEN					!operat_rule.dat not found
	write(*,*)pfadp(1:pfadj)// 'Reservoir/operat_rule.dat was not found. Run the model anyway.'
    DO i=1,subasin
	  if (damq_frac(i) == -999.) then
	    write(*,*)'[operat_rule.dat must be given [or change the value of the parameter damq_frac in reservoir.dat]'
		stop
	  endif
	ENDDO
  ELSE
    READ(11,*)
    READ(11,*)
    DO i=1,subasin
	  IF (damq_frac(i) == -999.) READ (11,*)dummy1,(dayexplot(i,s),s=1,4),(damq_frac_season(i,s),s=1,4)
	  IF (damq_frac(i) /= -999.) dummy1=id_subbas_extern(i)
      IF (dummy1 /= id_subbas_extern(i)) THEN
        WRITE(*,*) 'Sub-basin-IDs in file operat_rule.dat must have the same ordering scheme as in hymo.dat'
        STOP
      END IF
	ENDDO
  ENDIF

  OPEN(11,FILE=pfadp(1:pfadj)// 'Reservoir/operat_bottom.dat', IOSTAT=istate,STATUS='old')
  IF (istate/=0) THEN					!operat_bottom.dat not found
	write(*,*)pfadp(1:pfadj)// 'Reservoir/operat_bottom.dat was not found. Run the model anyway.'
    DO i=1,subasin
	  if (fvol_bottom(i) == -999.) then
	    write(*,*)'[operat_bottom.dat must be given [or change the value of the parameter damq_frac in reservoir.dat]'
		stop
	  endif
	ENDDO
  ELSE
    READ(11,*)
    READ(11,*)
    DO i=1,subasin
	  IF (fvol_bottom(i) == -999.) READ (11,*)dummy1,operat_start(i),operat_stop(i),operat_elev(i)
	  IF (fvol_bottom(i) /= -999.) dummy1=id_subbas_extern(i)
      IF (dummy1 /= id_subbas_extern(i)) THEN
        WRITE(*,*) 'Sub-basin-IDs in file operat_bottom.dat must have the same ordering scheme as in hymo.dat'
        STOP
      END IF
	ENDDO
  ENDIF

!Ge stage-volume curves for each sub-basin
  nbrbat(1:subasin)=0

  OPEN(11,FILE=pfadp(1:pfadj)// 'Reservoir/cav.dat', IOSTAT=istate,STATUS='old')
  IF (istate/=0) THEN					!cav.dat not found
	write(*,*)pfadp(1:pfadj)// 'Reservoir/cav.dat was not found. Run the model anyway.'
	flag_cav=0
  ELSE
	flag_cav=1
	READ(11,*)
 	READ(11,*)
    READ(11,*,IOSTAT=ka) dummy1,dummy2
    DO i=1,subasin
      IF (dummy1==id_subbas_extern(i)) THEN
         nbrbat(i)=dummy2
	  ENDIF
	ENDDO
    DO WHILE (ka==0)
	  READ(11,*);READ(11,*)
	  READ(11,*, IOSTAT=ka) dummy1,dummy2					!read next line in file
      DO i=1,subasin
        IF (dummy1==id_subbas_extern(i)) THEN
          nbrbat(i)=dummy2
		ENDIF
	  ENDDO
    END DO
  ENDIF
  CLOSE(11)
  i=maxval(nbrbat) !getmaximum number of rating curve points included and allocate necessary memory
  allocate(elev_bat0(i,subasin), &
		   area_bat0(i,subasin), &
			vol_bat0(i,subasin), &
			elev_bat(i,subasin), &
		   area_bat(i,subasin), &
			vol_bat(i,subasin), &
			STAT = istate)
	elev_bat0=0
	area_bat0=0
	vol_bat0=0
	elev_bat=0
	area_bat=0
	vol_bat=0

	if (istate/=0) then
		write(*,'(A,i0,a)')'Memory allocation error (',istate,') in reservoir-module (rating curves too detailed).'
		stop
	end if


  IF (flag_cav==1) THEN
   OPEN(11,FILE=pfadp(1:pfadj)// 'Reservoir/cav.dat',STATUS='unknown')
   READ(11,*)
   READ(11,*)
   j=2

   DO WHILE (.TRUE.)
	  READ (11,'(A)', IOSTAT=istate)fmtstr

	  if (istate/=0) exit !exit loop when no more line encountered
	  j=j+1
	  READ (fmtstr,*, IOSTAT=istate) dummy1, dummy2

	  i=which1(id_subbas_extern == dummy1)

	  IF (i==0) THEN
		  WRITE (*,'(A,I0,A)') 'WARNING: unknown upstream subbasin ID ', dummy1,' in cav.dat, ignored.'
		  READ(11,*,IOSTAT=istate) !skip next tow lines
          READ(11,*,IOSTAT=istate)
          j=j+2
          cycle
	  END IF

	  READ(fmtstr,*, IOSTAT=istate) dummy1,dummy2,(elev_bat0(ka,i),ka=1,dummy2)
      IF (istate/=0) THEN
		  WRITE (*,'(A,i0,A,A)') 'Format error in cav.dat, line',j,':', fmtstr
		  STOP
	  END IF
	  READ(11,*,IOSTAT=istate) dummy1,dummy2,(area_bat0(ka,i),ka=1,dummy2)
      j=j+1
	  IF (istate/=0) THEN
		  WRITE (*,'(A,i0,A,A)') 'Format error in cav.dat, line',j,':', fmtstr
		  STOP
	  END IF
	  READ(11,*,IOSTAT=istate) dummy1,dummy2,(vol_bat0(ka,i),ka=1,dummy2)
	  j=j+1

      IF (istate/=0) THEN
		  WRITE (*,'(A,i0,A,A)') 'Format error in cav.dat, line',j,':', fmtstr
		  STOP
	  END IF

      IF (maxlevel(i) > elev_bat0(dummy2,i)) THEN
        WRITE(*,*)'ERROR subasin ',id_subbas_extern(i),  &
            'MAXIMUM RESERVOIR LEVEL VALUE IS GREATER THAN THE MAXIMUM ELEVATION AT THE STAGE-AREA-VOLUME CURVE (FILE: cav.dat)'
        STOP
      ELSE IF (minlevel(i) < elev_bat0(1,i)) THEN
        WRITE(*,*)'ERROR subasin ',id_subbas_extern(i),  &
            'MINIMUM RESERVOIR LEVEL VALUE IS LESS THAN THE MINIMUM ELEVATION AT THE STAGE-AREA-VOLUME CURVE (FILE: cav.dat)'
        STOP
      END IF

   END DO
   CLOSE(11)
  ENDIF

  DO i=1,subasin
    nbrbat1=nbrbat(i)
    IF (nbrbat(i) /= 0) THEN
      DO j=1,nbrbat1
        area_bat0(j,i)=area_bat0(j,i)*1000.
        vol_bat0(j,i)=vol_bat0(j,i)*1000.
      END DO
    END IF
  END DO

!Ge initialization of the stage-volume curves for each sub-basin (erosion/deposition process)
  DO i=1,subasin
    nbrbat1=nbrbat(i)
    IF (nbrbat(i) /= 0) THEN
      DO j=1,nbrbat1
        elev_bat(j,i)=elev_bat0(j,i)
        area_bat(j,i)=area_bat0(j,i)
        vol_bat(j,i)=vol_bat0(j,i)
      END DO
    END IF
  END DO


  DO i=1,subasin
   IF (res_flag(i) /= -999.) THEN
    nbrbat1=nbrbat(i)
    IF (nbrbat(i) /= 0) THEN
      DO j=1,nbrbat(i)-1
        IF (damdead(i)*1.e6 >= vol_bat(j,i).AND.  &
            damdead(i)*1.e6 <= vol_bat(j+1,i)) THEN
          elevdead(i)=elev_bat(j,i)+((damdead(i)*1.e6)-vol_bat  &
                (j,i))/(vol_bat(j+1,i)-vol_bat(j,i))*  &
                (elev_bat(j+1,i)-elev_bat(j,i))
	    ENDIF
        IF (damalert(i)*1.e6 >= vol_bat(j,i).AND.  &
            damalert(i)*1.e6 <= vol_bat(j+1,i)) THEN
          elevalert(i)=elev_bat(j,i)+((damalert(i)*1.e6)-vol_bat  &
                (j,i))/(vol_bat(j+1,i)-vol_bat(j,i))*  &
                (elev_bat(j+1,i)-elev_bat(j,i))
	    ENDIF
	  ENDDO
    ELSE
	  if (damdead(i) /= 0.) then
        elevhelp=((damdead(i)*1.e6)/forma_factor(i))**(1./3.)
	  else
        elevhelp=0.
	  endif
	  elevdead(i)=elevhelp+minlevel(i)
	  if (damalert(i) /= 0.) then
        elevhelp=((damalert(i)*1.e6)/forma_factor(i))**(1./3.)
	  else
        elevhelp=0.
	  endif
	  elevalert(i)=elevhelp+minlevel(i)
    ENDIF
   ENDIF
  END DO

!Ge Initialization of the parameters related to spillway overflow
  DO i=1,subasin
	IF (res_flag(i) /= -999.) THEN
	  outflow_last(i)=0.
	  volume_last(i)=0.
	  alpha_over(i)=1./(1.-damb(i))
	  k_over(i)=(dama(i)/alpha_over(i))**alpha_over(i)
!write(*,*)id_subbas_extern(i),dama(i),damb(i),k_over(i),alpha_over(i),storcap(i)*1.e6
	  hmax(i)=((storcap(i)*1.e6)/k_over(i))**(1./alpha_over(i))
!write(*,'(I6,4F12.3,F10.3,F15.1)')id_subbas_extern(i),dama(i),damb(i),k_over(i),alpha_over(i),hmax(i),storcap(i)*1.e6
	ENDIF
  ENDDO

!Ge dosediment HAS TO BE INCLUDED INTO THE do.dat AND common.fi FILES, AND READ IN readgen.fi
  IF (dosediment) THEN
    CALL semres (status,upstream)
  END IF

!Ge initialization of output files
  if (f_river_velocity) then
    OPEN(11,FILE=pfadn(1:pfadi)//'River_Velocity.out',STATUS='replace')
	WRITE (11,*) 'Output file for flow velocity in m/s (with MAP IDs as in hymo.dat)'
	write(fmtstr,'(a,i0,a)')'(3a6,',subasin,'i14)'		!generate format string
	WRITE (11,fmtstr)' Year ', ' Day  ',' dt   ',(id_subbas_extern(i), i=1,subasin)
	Close (11)
  endif

  IF (f_res_watbal) then
      DO i=1,subasin
        IF (storcap(i) /= 0.) THEN
          WRITE(subarea,*)id_subbas_extern(i)
          OPEN(11,FILE=pfadn(1:pfadi)//'res_'//trim(adjustl(subarea))//'_watbal.out',STATUS='replace')
          WRITE(11,*)'Subasin-ID, year, day, hour, qlateral(m**3/s), inflow(m**3/s), intake(m**3/s), overflow(m**3/s), qbottom(m**3/s), qout(m**3/s), elevation(m), area(m**2), volume(m**3)'
          CLOSE(11)
        ENDIF
      ENDDO
  ENDIF

!Ge initialization of output files
  IF (f_res_vollost) then
      DO i=1,subasin
        IF (storcap(i) /= 0.) THEN
          WRITE(subarea,*)id_subbas_extern(i)
          OPEN(11,FILE=pfadn(1:pfadi)//'res_'//trim(adjustl(subarea))//'_vollost.out',STATUS='replace')
          WRITE(11,*)'Subasin-ID, year, day, hour, deadvol(m**3), alertvol(m**3), storcap(m**3)'
          CLOSE(11)
        ENDIF
      ENDDO
  ENDIF

!Ge initialization of output files
  IF (f_res_cav) then
      DO i=1,subasin
       IF (nbrbat(i) /= 0) THEN
        IF (storcap(i) /= 0.) THEN
          WRITE(subarea,*)id_subbas_extern(i)
          OPEN(11,FILE=pfadn(1:pfadi)//'res_'//trim(adjustl(subarea))//'_cav.out',STATUS='replace')
        ENDIF
       ENDIF
      ENDDO
  ENDIF

!**************************************************************************************
!Ge temporary output files to test the cascade routing scheme of the lake module
!  IF (.not. doacud) THEN
!   DO i=1,subasin
!    IF (storcap(i) /= 0.) THEN
!	  OPEN(11,FILE=pfadn(1:pfadi)//'res_watbal_lake.out',STATUS='replace')
!	    WRITE(11,*)'Subasin-ID, year, day, hour, inflow_classes(m**3), outflow_classes(m**3), retention_classes(m**3), volume_classes(m**3)'
!      CLOSE(11)
!	ENDIF
!   ENDDO
!   DO i=1,subasin
!    IF (dosediment) THEN
!     IF (storcap(i) /= 0.) THEN
!	  OPEN(11,FILE=pfadn(1:pfadi)//'res_sedbal_lake.out',STATUS='replace')
!	    WRITE(11,*)'Subasin-ID, year, day, hour, sedinflow_classes(ton), sedoutflow_classes(ton), sedretention_classes(ton), sedimentation_class'
!      CLOSE(11)
!	 ENDIF
!	ENDIF
!   ENDDO
!  ELSE
!   DO i=1,subasin
!    IF (storcap(i) /= 0.) THEN
!	  OPEN(11,FILE=pfadn(1:pfadi)//'res_watbal_lake.out',STATUS='replace')
!	    WRITE(11,*)'Subasin-ID, year, day, hour, inflow_strateg(m**3), outflow_strateg(m**3), retention_strateg(m**3), volume_strateg(m**3)'
!      CLOSE(11)
!	ENDIF
!   ENDDO
!   DO i=1,subasin
!    IF (dosediment) THEN
!     IF (storcap(i) /= 0.) THEN
!	  OPEN(11,FILE=pfadn(1:pfadi)//'res_sedbal_lake.out',STATUS='replace')
!	    WRITE(11,*)'Subasin-ID, year, day, hour, sedinflow_strateg(ton), sedoutflow_strateg(ton), sedretention_strateg(ton), sedimentation_class'
!      CLOSE(11)
!	 ENDIF
!	ENDIF
!   ENDDO
!  ENDIF
!**************************************************************************************


END IF


! -----------------------------------------------------------------------
IF (STATUS == 1) THEN
! initialize ...

!Ge initialization of parameters
  DO i=1,subasin
    DO id=1,dayyear*nt
      overflow(id,i)=0.
      qintake(id,i)=0.
      qinflow(id,i)=0.
	  qbottom(id,i)=0.
	  qlateral(id,i)=0.
	  volact(1,i)=0.
    END DO
    damareaact(i)=0.
  END DO

  !IF (t > tstart) THEN !Andreas
  DO i=1,subasin !Andreas
   IF (storcap(i) > 0.) THEN
     IF (t > damyear(i) .AND. t > tstart) THEN  !Andreas
       volact(1,i)=volact(daylastyear*nt,i)
       daystorcap(1,i)=daystorcap(daylastyear*nt,i)
       daydamalert(1,i)=daydamalert(daylastyear*nt,i)
       daydamdead(1,i)=daydamdead(daylastyear*nt,i)
       daymaxdamarea(1,i)=daymaxdamarea(daylastyear*nt,i)
	   dayminlevel(1,i)=dayminlevel(daylastyear*nt,i)

     ELSE IF (t == damyear(i) .OR. (t > damyear(i) .AND. t == tstart)) THEN  !Andreas
       IF (vol0(i) /= -999.) THEN
         volact(1,i)=vol0(i)
       ELSE
         volact(1,i)=storcap(i)/5.
       END IF
       daystorcap(1,i)=storcap(i)
       daydamalert(1,i)=damalert(i)
       daydamdead(1,i)=damdead(i)
       daymaxdamarea(1,i)=maxdamarea(i)
	   dayminlevel(1,i)=minlevel(i)
     ELSE
       volact(1,i)=0.
       daystorcap(1,i)=0.
       daydamalert(1,i)=0.
       daydamdead(1,i)=0.
       daymaxdamarea(1,i)=0.
       dayminlevel(1,i)=0.
	 ENDIF
   END IF
!write(*,*)id_subbas_extern(i),volact(1,i),storcap(i),vol0(i)
  END DO !Andreas


!Ge water availability approach for reservoirs has to be included
  avail_ac(:,:)=0. !water availability
  avail_all(:,:)=0.   !water availability
  damex(:,:)=0.

!  actual storage capacity in this year (as derived from the data in ??.dat)
!  storcapact=0.
!  DO i=1,subasin
!    IF (storcap(i) > 0.) THEN
!      IF (t >= damyear(i)) THEN
!        storcapact=storcapact+storcap(i)
!      END IF
!    END IF
!  END DO


!George reservoir water surface (m**2)
  DO i=1,subasin
    IF (storcap(i) > 0. .and. t >= damyear(i)) THEN
     nbrbat1=nbrbat(i)
     IF (nbrbat(i) /= 0) THEN
      DO j=1,nbrbat(i)-1
        IF (volact(1,i)*1.e6 >= vol_bat(j,i).AND.  &
            volact(1,i)*1.e6 <= vol_bat(j+1,i)) THEN
          damareaact(i)=area_bat(j,i)+((volact(1,i)*1.e6)-vol_bat  &
                (j,i))/(vol_bat(j+1,i)-vol_bat(j,i))*  &
                (area_bat(j+1,i)-area_bat(j,i))
	    ENDIF
	  ENDDO
     ELSE
	  if (volact(1,i) /= 0.) then
		areahelp=dama(i)*((volact(1,i)*1.e6)**damb(i))
	  else
        areahelp=0.
	  endif
	  damareaact(i)=areahelp
     ENDIF
	ENDIF !Andreas
!write(*,*)id_subbas_extern(i),damareaact(i),volact(1,i)*1.e6,dama(i),damb(i)
  END DO
!  stop

!Ge read daily data on reservoir level and outflow discharges
 IF (reservoir_check == 1) THEN
   IF (reservoir_balance == 0) THEN
    DO i=1,subasin
     IF (t >= damyear(i)) THEN !Andreas
	  IF (storcap(i) > 0.) THEN
        WRITE(subarea,*)id_subbas_extern(i)
		OPEN(11,FILE=pfadp(1:pfadj)//'Reservoir/inflow_'//trim(adjustl(subarea))//'.dat', &
			STATUS='unknown')
          read(11,*)
          cont=(dtot-dayyear)*nt
          DO id=1,cont
            READ(11,*)
          END DO
          DO id=1,dayyear*nt
            READ(11,*) dummy1,dummy2,r_precip,r_etp,r_qinflow, &
				r_overflow,r_qintake,r_qbottom,r_level0,r_level1
	        damelev0(id,i)=r_level0
	        damelev1(id,i)=r_level1
            overflow(id,i)=r_overflow
            qintake(id,i)=r_qintake
	        qbottom(id,i)=r_qbottom
            res_precip(id,i)=r_precip
	        res_pet(id,i)=r_etp
	        qinflow(id,i)=r_qinflow
	      ENDDO
	    CLOSE(11)
        IF (nbrbat(i) /= 0) THEN
	      nbrbat1=nbrbat(i)
          DO id=1,dayyear
            IF (damelev0(id,i) > elev_bat0(nbrbat1,i) .or. damelev1(id,i) > elev_bat0(nbrbat1,i) ) THEN
              WRITE(*,*)'ERROR subasin ',id_subbas_extern(i),' year ',t,' day ',id, &
			     	'GIVEN VALUE OF DAILY RESERVOIR LEVEL IS GREATER THAN THE MAXIMUM RESERVOIR ELEVATION AT THE STAGE-AREA-VOLUME CURVE (FILE: cav.dat)'
			  STOP
			ELSE IF (damelev0(id,i) < elev_bat0(1,i) .or. damelev1(id,i) < elev_bat0(1,i)) THEN
			  WRITE(*,*)'ERROR subasin ',id_subbas_extern(i),' year ',t,' day ',id,  &
					'GIVEN VALUE OF DAILY RESERVOIR LEVEL IS LESS THAN THE MINIMUM ELEVATION AT THE STAGE-AREA-VOLUME CURVE (FILE: cav.dat)'
			  STOP
			END IF
		  ENDDO
		ENDIF
	  ENDIF
	 ENDIF !Andreas
	ENDDO
   ENDIF
   IF (reservoir_balance == 1) THEN
    DO i=1,subasin
     IF (t >= damyear(i)) THEN !Andreas
	  IF (storcap(i) > 0.) THEN
        WRITE(subarea,*)id_subbas_extern(i)
		OPEN(11,FILE=pfadp(1:pfadj)//'Reservoir/inflow_'//trim(adjustl(subarea))//'.dat', &
			STATUS='unknown')
          read(11,*)
          cont=(dtot-dayyear)*nt
          DO id=1,cont
            READ(11,*)
          END DO
          DO id=1,dayyear*nt
            READ(11,*) dummy1,dummy2,r_precip,r_etp,r_qinflow
            res_precip(id,i)=r_precip
	        res_pet(id,i)=r_etp
	        qinflow(id,i)=r_qinflow
	      ENDDO
	    CLOSE(11)
	  ENDIF
	 ENDIF !Andreas
	ENDDO
   ENDIF
  ELSE
   DO i=1,subasin
    IF (t >= damyear(i)) THEN !Andreas
     dummy1=0
     DO id=1,dayyear
	  if (river_transport.ne.1) then
       DO h=1,nt
	    dummy1=dummy1+1
	    res_pet(dummy1,i)=pet(id,i)/nt
       ENDDO
	  else
	   res_pet(id,i)=pet(id,i)
	  endif
	 ENDDO
    ENDIF !Andreas
   ENDDO
   DO i=1,subasin
    IF (t >= damyear(i)) THEN !Andreas
	 IF (dohour) THEN
	  IF (river_transport.ne.1) THEN
       DO id=1,dayyear*nt
	    res_precip(id,i)=preciph(id,i)
	   ENDDO
	  ENDIF
	  IF (river_transport.eq.1) THEN
       DO id=1,dayyear
	    dummy4=0.
	    DO h=1,nt
		  dummy4=dummy4+preciph(id,i)
	    ENDDO
	    res_precip(id,i)=dummy4
	   ENDDO
	  ENDIF
	 ELSE
      DO id=1,dayyear
	    res_precip(id,i)=precip(id,i)
	  ENDDO
	 ENDIF
    ENDIF !Andreas
   ENDDO
 ENDIF


! begin block Andreas
! Read daily data on (measured) regulated reservoir outflow to be considered in reservoir water balance
! data given in m**3/s, convert to m**3/d
  IF (reservoir_balance == 1) THEN
   DO i=1,subasin
    IF (t >= damyear(i)) THEN
     IF (damq_frac(i) == -888.) THEN
      WRITE(subarea,*)id_subbas_extern(i)
      OPEN(11,FILE=pfadp(1:pfadj)//'Reservoir/intake_'//trim(adjustl(subarea))//'.dat', &
			STATUS='unknown')
       read(11,*)
	   read(11,*)
       cont=(dtot-dayyear)*nt
       DO id=1,cont
         READ(11,*)
       END DO
       DO id=1,dayyear*nt
         READ(11,*) dummy1,r_qintake
		 IF (qintake(id,i) /= -999.) qintake(id,i)=r_qintake*(86400./nt)
		 IF (qintake(id,i) == -999.) qintake(id,i)=damflow(i)*damq_frac(i)
       ENDDO
      CLOSE(11)
     ENDIF
    ENDIF
   ENDDO
  ENDIF
! end block Andreas


  IF (dosediment) THEN
    CALL semres (status,upstream)
  END IF




END IF

! -----------------------------------------------------------------------
IF (STATUS == 2) THEN

! simulation timestep
  if (river_transport.ne.1) then
    hour=res_h
    step=(d-1)*nt+hour
  else
    hour=1
	step=d
  endif

! Computation of reservoir water balance
  IF (storcap(upstream) > 0.) THEN
    IF (reservoir_balance == 1) THEN
      IF (reservoir_check == 1) THEN
	    qinflow(step,upstream)=qinflow(step,upstream)*(86400./nt)
	  ELSE
        if (river_transport.eq.1)then
	      qinflow(step,upstream)=qout(step,upstream)*(86400./nt)
	    else
          qinflow(step,upstream)=(r_qout(2,upstream)+qlateral(step,upstream))*(86400./nt)
	    endif
	  ENDIF

!qinflow(step,upstream)=20*86400.
!write(*,*)upstream,id_subbas_extern(upstream)
!if(step<=10)qinflow(step,upstream)=storcap(upstream)*1.e6
!if(step>20)qinflow(step,upstream)=storcap(upstream)*1.e6
!if (step == 16 .and. t==2002) read(*,*)

      volact(step,upstream)=volact(step,upstream)*1.e6
      daymaxdamarea(step,upstream)=daymaxdamarea(step,upstream)*1.e4
      daystorcap(step,upstream)=daystorcap(step,upstream)*1.e6
      daydamalert(step,upstream)=daydamalert(step,upstream)*1.e6
      daydamdead(step,upstream)=daydamdead(step,upstream)*1.e6
      damflow(upstream)=damflow(upstream)*(86400./nt)
      damvol0(upstream)=volact(step,upstream)
	  qoutlet(upstream)=qoutlet(upstream)*(86400./nt)
	  withdrawal(upstream)=withdrawal(upstream)*(86400./nt)

!write(*,'(2I4,3F15.4)')step,id_subbas_extern(upstream),volact(step,upstream)


! 1) Calculation of the actual reservoir volume after water inflow
	  help3=volact(step,upstream)+qinflow(step,upstream)
	  help2=volact(step,upstream)
!write(*,'(2I4,4F15.4)')step,id_subbas_extern(upstream),help3,volact(step,upstream),daystorcap(step,upstream),overflow(step,upstream)
      IF (help3 > daystorcap(step,upstream)) THEN
!write(*,'(2I4,4F15.4)')step,id_subbas_extern(upstream),volact(step,upstream),help3,help,daystorcap(step,upstream)
!write(*,'(2I4,4F15.4)')step,id_subbas_extern(upstream),fvol_over(upstream)*daystorcap(step,upstream),help3,daystorcap(step,upstream)
		IF (fvol_over(upstream)==1) THEN
          help=daystorcap(step,upstream)
!write(*,'(2I4,3F15.4)')step,id_subbas_extern(upstream),overflow(step,upstream)/86400.
!write(*,'(2I4,4F15.4)')step,id_subbas_extern(upstream),volact(step,upstream),help3,help,daystorcap(step,upstream)
		ELSE
		  help=daystorcap(step,upstream)
		  overflow(step,upstream)=help3-daystorcap(step,upstream)
        END IF
	    lakeret(step,upstream)=max(0.,help-volact(step,upstream))
        volact(step,upstream)=help
	  ELSE
	    lakeret(step,upstream)=qinflow(step,upstream)
        volact(step,upstream)=help3
      END IF
!write(*,'(2I4,4F15.4)')step,id_subbas_extern(upstream),help,volact(step,upstream),daystorcap(step,upstream),overflow(step,upstream)
!if (d==31)stop
!write(*,'(2I4,4F15.3)')d,id_subbas_extern(upstream),fvol_over(upstream),daystorcap(step,upstream),volact(step,upstream)+qinflow(step,upstream),overflow(step,upstream)
!write(*,'(2I4,4F15.3)')d,id_subbas_extern(upstream),volact(step,upstream)

!write(*,'(2I4,4F15.4)')step,id_subbas_extern(upstream),volact(step,upstream),help3,help,daystorcap(step,upstream)

! 2) Substract evaporation out the reservoir
! Calculation of the actual reservoir area by interpolation (using the stage-area-volume curve)
      IF (nbrbat(upstream) /= 0) THEN
        DO j=1,nbrbat(upstream)-1
          IF (volact(step,upstream) >= vol_bat(j,upstream).AND.  &
                volact(step,upstream) <= vol_bat(j+1,upstream)) THEN
            areahelp=area_bat(j,upstream)+(volact(step,upstream)-vol_bat  &
                (j,upstream))/(vol_bat(j+1,upstream)-vol_bat(j,upstream))*  &
                (area_bat(j+1,upstream)-area_bat(j,upstream))
            elevhelp=elev_bat(j,upstream)+(volact(step,upstream)-vol_bat  &
                (j,upstream))/(vol_bat(j+1,upstream)-vol_bat(j,upstream))*  &
                (elev_bat(j+1,upstream)-elev_bat(j,upstream))
          END IF
        END DO
! Calculation of the actual reservoir area by interpolation (using the morphologic parameter alpha)
      ELSE
        elevhelp=(volact(step,upstream)/forma_factor(upstream))**(1./3.)
		elevhelp=elevhelp+dayminlevel(step,upstream)
        areahelp=dama(upstream)*(volact(step,upstream)**damb(upstream))
      END IF
      damelevact(upstream)=elevhelp
      damareaact(upstream)=areahelp
!write(*,'(2I4,4F15.3)')d,id_subbas_extern(upstream),damelevact(upstream),damareaact(upstream)


! Calculation of increase or decrease of the reservoir level due to precip and evap.
! (using the stage-area-volume curve)
      evaphelp2=0.
      IF (nbrbat(upstream) /= 0) THEN
        damelevact(upstream)=damelevact(upstream)+  &
            (res_precip(step,upstream)-res_pet(step,upstream))/1000.
        IF (damelevact(upstream) < dayminlevel(step,upstream)) THEN
          evaphelp2=dayminlevel(step,upstream)-damelevact(upstream)
          damelevact(upstream)=dayminlevel(step,upstream)
        END IF

! Check overflow due to precipitation
!George        IF (damelevact(upstream) > maxlevel(upstream)) THEN
!George          overflow(step,upstream)=overflow(step,upstream)+((damelevact(upstream)-maxlevel(upstream))*  &
!George              daymaxdamarea(step,upstream))
!George          damelevact(upstream)=maxlevel(upstream)
!George        END IF
        DO j=1,nbrbat(upstream)-1
          IF (damelevact(upstream) >= elev_bat(j,upstream).AND.  &
                damelevact(upstream) <= elev_bat(j+1,upstream)) THEN
            volhelp=vol_bat(j,upstream)+(damelevact(upstream)-elev_bat  &
                (j,upstream))/(elev_bat(j+1,upstream)-elev_bat(j,upstream))*  &
                (vol_bat(j+1,upstream)-vol_bat(j,upstream))
            areahelp=area_bat(j,upstream)+(damelevact(upstream)-elev_bat  &
                (j,upstream))/(elev_bat(j+1,upstream)-elev_bat(j,upstream))*  &
                (area_bat(j+1,upstream)-area_bat(j,upstream))
          END IF
        END DO

! Calculation of evaporation and precipitation using the truncated cone volume (m3)
! (using the morphologic parameter alpha)
        evaphelp=(areahelp+SQRT(areahelp*damareaact(upstream))+  &
            damareaact(upstream))*res_pet(step,upstream)/1000.*1./3.
        prechelp=(areahelp+SQRT(areahelp*damareaact(upstream))+  &
            damareaact(upstream))*res_precip(step,upstream)/1000.*1./3.
        infhelp=0.
        volact(step,upstream)=volhelp

      ELSE
        evaphelp=MIN(volact(step,upstream),(res_pet(step,upstream)/1000.)*areahelp)
        prechelp=(res_precip(step,upstream)/1000.)*areahelp
        infhelp=0.
        volact(step,upstream)=volact(step,upstream)+ (prechelp-evaphelp-infhelp)

! Check overflow due to precipitation
!        IF (volact(step,upstream) > daystorcap(step,upstream)) THEN
!          overflow(step,upstream)=overflow(step,upstream)+volact(step,upstream)-daystorcap(step,upstream)
!          volact(step,upstream)=daystorcap(step,upstream)
!        END IF
      END IF
!write(*,'(2I4,4F15.3)')d,id_subbas_extern(upstream),overflow(step,upstream),volact(step,upstream)
!write(*,'(2I4,3F15.4)')step,id_subbas_extern(upstream),volact(step,upstream)


! 3) reservoir evaporation in mm
      IF (nbrbat(upstream) /= 0) THEN
        evapdam(step,upstream)=MAX(0.,res_pet(step,upstream)-(evaphelp2*1.e3))
      ELSE
        evapdam(step,upstream)=evaphelp/areahelp*1.e3
      END IF
!     reservoir evaporation in Mio m**3
      etdam(step,upstream)=evaphelp/1.e6
      infdam(step,upstream)=infhelp/1E6


! 4) Check flow through the reservoir
! 4a) Water intake device
!     - dam outflow is fraction of Q90
!     - according to J.C.Araújo this fraction is estimated to be 90%,
!       and 80% for larger, strategic dams (e.g. Oros)
!     - if alert volume is reached (if given for dam), outflow is reduced
!     - no outflow if storage is below dead volume (if given for dam)
!     - for non-strategic dams outflow is demand of regular (mean precip) years
!       (not yet implemented)

! Calculation of the maximum controlled outflow discharge using a factor defined in the reservoir.dat
	  IF (damq_frac(upstream) >= 0.0) THEN !Andreas
	    helpout=damflow(upstream)*damq_frac(upstream)
! Calculation of the maximum controlled outflow discharge using a operation regime as provided in the operat_rule.dat
	  ELSE IF (damq_frac(upstream) == -999.) THEN !Andreas
	    dummy4=0.
	    do s=1,3
	      IF (dayoutsim+d >= dayexplot(upstream,s) .and. &
		      dayoutsim+d < dayexplot(upstream,s+1)) dummy4=damq_frac_season(upstream,s)
        enddo
		IF (dayoutsim+d < dayexplot(upstream,1) .or. &
		    dayoutsim+d >= dayexplot(upstream,4)) dummy4=damq_frac_season(upstream,4)
        helpout=damflow(upstream)*dummy4
!write(*,*)upstream,dayoutsim+d,(dayexplot(upstream,s),s=1,4)
!write(*,*)d,upstream,helpout,damflow(upstream),dummy4
!write(*,*)step,id_subbas_extern(upstream),helpout/86400.,dummy4,damflow(upstream)/86400.
!     use explicit (measured) daily values of regulated reservoir outflow
	  ELSE IF (damq_frac(upstream) == -888.) THEN  !Andreas
        helpout=qintake(step,upstream)                !Andreas
	  ENDIF

! Check water availability
      IF (daydamalert(step,upstream) > daydamdead(step,upstream)) THEN
        IF (volact(step,upstream) < daydamalert(step,upstream)) THEN
          IF (volact(step,upstream) > daydamdead(step,upstream)) THEN
            helpout=helpout*(volact(step,upstream)-daydamdead(step,upstream))/  &
                (daydamalert(step,upstream)-daydamdead(step,upstream))
          END IF
        END IF
      ENDIF
      IF (volact(step,upstream) < (daydamdead(step,upstream)+helpout))THEN
        helpout=volact(step,upstream)-daydamdead(step,upstream)
      END IF

      IF (volact(step,upstream) <= daydamdead(step,upstream)) helpout=0.
      IF (elevdead(upstream) <= dayminlevel(step,upstream)) helpout=0.
 !write(*,*)step,id_subbas_extern(upstream),helpout/86400.
      qintake(step,upstream)=helpout
      volact(step,upstream)=MAX(0.,volact(step,upstream)-qintake(step,upstream))

!write(*,'(2I4,3F15.4)')step,id_subbas_extern(upstream),volact(step,upstream)


! 4b) Bottom outlets
      IF (nbrbat(upstream) /= 0) THEN
        DO j=1,nbrbat(upstream)-1
          IF (operat_elev(upstream) >= elev_bat(j,upstream).AND.  &
                operat_elev(upstream) <= elev_bat(j+1,upstream)) THEN
            volhelp=vol_bat(j,upstream)+(operat_elev(upstream)-elev_bat  &
                (j,upstream))/(elev_bat(j+1,upstream)-elev_bat(j,upstream))*  &
                (vol_bat(j+1,upstream)-vol_bat(j,upstream))
         END IF
        END DO
      ELSE
        elevhelp=operat_elev(upstream)-dayminlevel(step,upstream)
        volhelp=forma_factor(upstream)*(elevhelp**3)
      END IF

	  IF (fvol_bottom(upstream) /= -999.) THEN
	   helpout=qoutlet(upstream)
       IF (volact(step,upstream) > daydamdead(step,upstream)) THEN
        IF (volact(step,upstream) < daystorcap(step,upstream)) THEN
          IF (volact(step,upstream) > fvol_bottom(upstream)*daystorcap(step,upstream)) THEN
            helpout=min(helpout,helpout*(volact(step,upstream)-(fvol_bottom(upstream)*daystorcap(step,upstream)))/ &
		             ((1.-fvol_bottom(upstream))*daystorcap(step,upstream)))
		  ELSE
		    helpout=0.
		  ENDIF
        ENDIF
       ELSE
		helpout=0.
       END IF
	  ELSE
	   IF (step >= operat_start(upstream) .and. step <= operat_stop(upstream)) THEN
	    IF (volact(step,upstream) > volhelp) THEN
		  helpout=min(qoutlet(upstream),(volact(step,upstream)-volhelp))
        ELSE
		  helpout=0.
		ENDIF
	   ELSE
	    helpout=0.
	   ENDIF
	  ENDIF

!write(*,'(3I4,3F15.4)')step,operat_start(upstream),operat_stop(upstream),fvol_bottom(upstream)

!      IF (elevbottom(upstream) <= dayminlevel(step,upstream)) helpout=0.

!write(*,*)step,id_subbas_extern(upstream),helpout,fvol_bottom(upstream),volact(step,upstream),daystorcap(step,upstream)
!write(*,'(2I4,6F12.2)')step,id_subbas_extern(upstream),helpout,fvol_bottom(upstream),volact(step,upstream),daystorcap(step,upstream)

      qbottom(step,upstream)=helpout
      volact(step,upstream)=MAX(0.,volact(step,upstream)-qbottom(step,upstream))
!write(*,'(2I4,3F15.4)')step,id_subbas_extern(upstream),volact(step,upstream)


! 4c) Withdrawal water volume to supply the water use sectors
      IF (volact(step,upstream) > .05*daystorcap(step,upstream)) THEN
        volact(step,upstream)=MAX(0.,volact(step,upstream)-withdrawal(upstream))
	  ENDIF

!write(*,'(2I4,3F15.4)')step,id_subbas_extern(upstream),volact(step,upstream)


! 5) Calculation of the overflow discharges from the reservoir
      IF (help3 > daystorcap(step,upstream) .and. fvol_over(upstream) == 1) THEN
        help=(volact(step,upstream)-daystorcap(step,upstream))+qinflow(step,upstream)
		help1=help2+help
        help=max(0.,help)
!write(*,'(2I4,6F11.1)')step,id_subbas_extern(upstream),help,help1,help2,help3,volact(step,upstream),daystorcap(step,upstream)
        IF (help1 > daystorcap(step,upstream)) THEN
!write(*,'(2I4,6F15.2)')step,id_subbas_extern(upstream),help,help2
         call reservoir_routing(upstream,help,help2)
!         help=max(0.,help-overflow(step,upstream))
         if(help2 > daystorcap(step,upstream)) lakeret(step,upstream)=lakeret(step,upstream)+(max(0.,volact(step,upstream)-help2))
         if(help2 <= daystorcap(step,upstream)) lakeret(step,upstream)=lakeret(step,upstream)+(max(0.,volact(step,upstream)-daystorcap(step,upstream)))

!write(*,'(2I4,3F15.4)')step,id_subbas_extern(upstream),help,volact(step,upstream),overflow(step,upstream)/86400.
        ELSE
		 volume_last(upstream)=0.
		 outflow_last(upstream)=0.
		 volact(step,upstream)=help1
		 overflow(step,upstream)=0.
        END IF
	  ELSE IF (volact(step,upstream) > daystorcap(step,upstream)) THEN
	    volact(step,upstream)=daystorcap(step,upstream)
		overflow(step,upstream)=overflow(step,upstream)+volact(step,upstream)-daystorcap(step,upstream)
      END IF


!write(*,'(2I4,5F15.4)')step,id_subbas_extern(upstream),qinflow(step,upstream)/86400.,overflow(step,upstream)/86400.,volact(step,upstream),lakeret(step,upstream)/86400.

!write(*,'(2I4,3F15.4)')step,id_subbas_extern(upstream),volact(step,upstream)

!if (step==20)stop
!if (step==41 .and. id_subbas_extern(upstream)==29) stop

      res_qout(step,upstream)=qintake(step,upstream)+overflow(step,upstream)+qbottom(step,upstream)



      IF (nbrbat(upstream) /= 0) THEN
        DO j=1,nbrbat(upstream)-1
          IF (volact(step,upstream) >= vol_bat(j,upstream).AND.  &
                volact(step,upstream) <= vol_bat(j+1,upstream)) THEN
            areahelp=area_bat(j,upstream)+(volact(step,upstream)-vol_bat  &
                (j,upstream))/(vol_bat(j+1,upstream)-vol_bat(j,upstream))*  &
                (area_bat(j+1,upstream)-area_bat(j,upstream))
            elevhelp=elev_bat(j,upstream)+(volact(step,upstream)-vol_bat  &
                (j,upstream))/(vol_bat(j+1,upstream)-vol_bat(j,upstream))*  &
                (elev_bat(j+1,upstream)-elev_bat(j,upstream))
          END IF
        END DO
      ELSE
        elevhelp=(volact(step,upstream)/forma_factor(upstream))**(1./3.)
		elevhelp=elevhelp+dayminlevel(step,upstream)
        areahelp=dama(upstream)*(volact(step,upstream)**damb(upstream))
      END IF
      damareaact(upstream)=areahelp
      damelevact(upstream)=elevhelp

!write(*,'(I4,4F15.4)')d,dayminlevel(step,upstream),elevhelp,forma_factor(upstream),areahelp

      res_qout(step,upstream)=res_qout(step,upstream)/(86400./nt)
      qinflow(step,upstream)=qinflow(step,upstream)/(86400./nt)
      qintake(step,upstream)=qintake(step,upstream)/(86400./nt)
      overflow(step,upstream)=overflow(step,upstream)/(86400./nt)
	  qbottom(step,upstream)=qbottom(step,upstream)/(86400./nt)


! Calculation of reservoir surface area and reservoir volume when inflow discharges, outflow discharges
! and reservoir levels are still provided in the res_"ID_SUBBAS_EXTERN"_daily.dat
    ELSE IF (reservoir_balance == 0) THEN

!      write(*,*)t,d,hour,upstream,qinflow(step,upstream)
      daymaxdamarea(step,upstream)=daymaxdamarea(step,upstream)*1.e4
      daystorcap(step,upstream)=daystorcap(step,upstream)*1.e6
      daydamalert(step,upstream)=daydamalert(step,upstream)*1.e6
      daydamdead(step,upstream)=daydamdead(step,upstream)*1.e6
      damvol0(upstream)=volact(step,upstream)

	  res_qout(step,upstream)=qintake(step,upstream)+overflow(step,upstream)+qbottom(step,upstream)
	  damelevact(upstream)=damelev1(step,upstream)

      IF (nbrbat(upstream) /= 0) THEN
        DO j=1,nbrbat(upstream)-1
          IF (damelevact(upstream) >= elev_bat(j,upstream).AND.  &
                damelevact(upstream) <= elev_bat(j+1,upstream)) THEN
            volhelp=vol_bat(j,upstream)+(damelevact(upstream)-elev_bat  &
                (j,upstream))/(elev_bat(j+1,upstream)-elev_bat(j,upstream))*  &
                (vol_bat(j+1,upstream)-vol_bat(j,upstream))
            areahelp=area_bat(j,upstream)+(damelevact(upstream)-elev_bat  &
                (j,upstream))/(elev_bat(j+1,upstream)-elev_bat(j,upstream))*  &
                (area_bat(j+1,upstream)-area_bat(j,upstream))

         END IF
        END DO
      ELSE
        elevhelp=damelevact(upstream)-dayminlevel(step,upstream)
        volhelp=forma_factor(upstream)*(elevhelp**3)
        areahelp=dama(upstream)*(volhelp**damb(upstream))
      END IF
      damareaact(upstream)=areahelp
	  volact(step,upstream)=volhelp
    ENDIF
    damex(step,upstream)=res_qout(step,upstream)+withdrawal(upstream)

! Call sediment balance sub-routine
    IF (dosediment) then
      CALL semres (status,upstream)

!write(*,'(2I4,3F15.4)')d,id_subbas_extern(upstream),decstorcap(step,upstream),dayminlevel(step,upstream),daystorcap(step,upstream)

! Calculation of storage capacity reduction
      if (decstorcap(step,upstream)/=0.) then

! CASE 1: stage-area-volume curve is not provided
       if (nbrbat(upstream) == 0) then
	    if (decstorcap(step,upstream) >= daystorcap(step,upstream) .or. daystorcap(step,upstream) == 0.) then
	     write(*,*) 'the resevoir located at the outlet point of sub-basin:',id_subbas_extern(upstream)
	     write(*,*) 'lost its total storage capacity due to sediment deposition'
		 storcap(upstream)=0.
	     daystorcap(step,upstream)=0.
		 do j=1,nbrbat(upstream)
          area_bat(j,upstream)=0.
          vol_bat(j,upstream)=0.
		 enddo
        else
         volhelp=MAX(0.,daystorcap(step,upstream)  &
			-decstorcap(step,upstream))

! Calculation of minimum reservoir elevation
         if (daystorcap(step,upstream) == 0.) then
	      dummy4=0.
	     else
	      dummy4=decstorcap(step,upstream)/daystorcap(step,upstream)
	     endif

	     elevhelp=max(0.,(maxlevel(upstream)-dayminlevel(step,upstream))*dummy4)
	     dayminlevel(step,upstream)=dayminlevel(step,upstream)+elevhelp
	     forma_factor(upstream)=volhelp/((maxlevel(upstream)-dayminlevel(step,upstream))**3)

	     daystorcap(step,upstream)=volhelp

	     elevhelp=max(0.,maxlevel(upstream)-dayminlevel(step,upstream))
         daymaxdamarea(step,upstream)=dama(upstream)*(volhelp**damb(upstream))

! Calculation of dead water volume after sediment deposition
	     if (damdead(upstream) /= 0. .and. daydamdead(step,upstream) /= 0.) then
	      elevhelp=max(0.,elevdead(upstream)-dayminlevel(step,upstream))
          daydamdead(step,upstream)=forma_factor(upstream)*(elevhelp**3.)
	     else
	      daydamdead(step,upstream)=0.
	     endif

! Calculation of alert water volume after sediment deposition
 	     if (damalert(upstream) /= 0. .and. daydamalert(step,upstream) /= 0.) then
	      elevhelp=max(0.,elevalert(upstream)-dayminlevel(step,upstream))
          daydamalert(step,upstream)=forma_factor(upstream)*(elevhelp**3.)
	     else
	      daydamalert(step,upstream)=0.
	     endif
	    endif

! CASE 2: stage-area-volume curve is provided
	   else
! CASE 2a: no detailed information of the reservoir geometry is provided (no cross section)
	    if (nbrsec(upstream) == 0) then
	     if (decstorcap(step,upstream) >= daystorcap(step,upstream) .or. daystorcap(step,upstream) == 0.) then
	      write(*,*) 'the resevoir located at the outlet point of sub-basin:',id_subbas_extern(upstream)
	      write(*,*) 'lost its total storage capacity due to sediment deposition'
		  storcap(upstream)=0.
	      daystorcap(step,upstream)=0.
		  do j=1,nbrbat(upstream)
           area_bat(j,upstream)=0.
           vol_bat(j,upstream)=0.
		  enddo
         else

	     daystorcap(step,upstream)=MAX(0.,daystorcap(step,upstream)  &
			-decstorcap(step,upstream))

! Change on stage-area-volume curve due to erosion and deposition processes
! Calculation of minimum reservoir elevation
	     elevhelp=dayminlevel(step,upstream)
		 p=0
         do j=1,nbrbat(upstream)
		  if (elev_bat(j,upstream) >= maxlevel(upstream)) then
		   dummy1=j
		    exit
		  endif
		 enddo
         do j=1,nbrbat(upstream)
	      if (elev_bat(j,upstream) < maxlevel(upstream)) then
		   dummy5=real(j)/real(dummy1)
		  else if (elev_bat(j,upstream) >= maxlevel(upstream)) then
		   dummy5=1.
		  endif
!write(*,'(3I4,F15.5,3F15.2)')step,j,dummy1,dummy5,elev_bat(j,upstream),dayminlevel(step,upstream),maxlevel(upstream)
!write(*,'(3I4,F15.5,4F15.2)')step,j,dummy1,dummy5,vol_bat(j,upstream),decstorcap(step,upstream)*dummy5,vol_bat(j,upstream)-(decstorcap(step,upstream)*dummy5)
		  dummy5=max(0.,vol_bat(j,upstream)-(decstorcap(step,upstream)*dummy5))
		  if (dummy5 == 0.)p=j
		  if (j == p+1) then
           DO q=1,nbrbat(upstream)-1
            IF (vol_bat(j,upstream)-dummy5 >= vol_bat(q,upstream).AND.  &
                vol_bat(j,upstream)-dummy5 <= vol_bat(q+1,upstream)) THEN
             elevhelp=elev_bat(q,upstream)+((vol_bat(j,upstream)-dummy5)-vol_bat  &
                (q,upstream))/(vol_bat(q+1,upstream)-vol_bat(q,upstream))*  &
                (elev_bat(q+1,upstream)-elev_bat(q,upstream))
            END IF
           END DO
		  endif
		  if(vol_bat(j,upstream)/=0.)dummy6=dummy5/vol_bat(j,upstream)
		  if(vol_bat(j,upstream)==0.)dummy6=0.
		  vol_bat(j,upstream)=dummy5
		  area_bat(j,upstream)=area_bat(j,upstream)*(dummy6**(2./3.)) ! V2/V1=h2**3/h1**3 and A2/A1=h2**2/h1**2 => A2/A1=(V2/V1)**(2/3)
	     enddo
		 dayminlevel(step,upstream)=elevhelp

         do j=1,nbrbat(upstream)
		  if (elev_bat(j,upstream) <= elevhelp) then
		   vol_bat(j,upstream)=0.
		   area_bat(j,upstream)=0.
		  endif
	     enddo

! Calculation of dead water volume after sediment deposition
         do j=1,nbrbat(upstream)-1
	      if (damdead(upstream) /= 0. .and. daydamdead(step,upstream) /= 0.) then
           if (elevdead(upstream) >= elev_bat(j,upstream).AND.  &
               elevdead(upstream) <= elev_bat(j+1,upstream)) THEN
            daydamdead(step,upstream)=vol_bat(j,upstream)+(elevdead(upstream)-elev_bat  &
                (j,upstream))/(elev_bat(j+1,upstream)-elev_bat(j,upstream))*  &
                (vol_bat(j+1,upstream)-vol_bat(j,upstream))
           endif
	      else
	       daydamdead(step,upstream)=0.
	      endif
! Calculation of alert water volume after sediment deposition
 	      if (damalert(upstream) /= 0. .and. daydamalert(step,upstream) /= 0.) then
           if (elevalert(upstream) >= elev_bat(j,upstream).AND.  &
               elevalert(upstream) <= elev_bat(j+1,upstream)) THEN
            daydamalert(step,upstream)=vol_bat(j,upstream)+(elevalert(upstream)-elev_bat  &
                (j,upstream))/(elev_bat(j+1,upstream)-elev_bat(j,upstream))*  &
                (vol_bat(j+1,upstream)-vol_bat(j,upstream))
           endif
	      else
	       daydamalert(step,upstream)=0.
	      endif
! Calculation of maximum reservoir area after sediment deposition
 	      if (daymaxdamarea(step,upstream) /= 0.) then
           if (maxlevel(upstream) >= elev_bat(j,upstream).AND.  &
               maxlevel(upstream) <= elev_bat(j+1,upstream)) THEN
            daymaxdamarea(step,upstream)=area_bat(j,upstream)+(maxlevel(upstream)-elev_bat  &
                (j,upstream))/(elev_bat(j+1,upstream)-elev_bat(j,upstream))*  &
                (area_bat(j+1,upstream)-area_bat(j,upstream))
           endif
	      else
	       daymaxdamarea(step,upstream)=0.
	      endif
	     enddo
	    endif

! CASE 2b: detailed information of the reservoir geometry is provided (cross sections)
       else
	    if (decstorcap(step,upstream) >= daystorcap(step,upstream) .or. daystorcap(step,upstream) == 0.) then
	     write(*,*) 'the resevoir located at the outlet point of sub-basin:',id_subbas_extern(upstream)
	     write(*,*) 'lost its total storage capacity due to sediment deposition'
		 storcap(upstream)=0.
	     daystorcap(step,upstream)=0.
		 do j=1,nbrbat(upstream)
          area_bat(j,upstream)=0.
          vol_bat(j,upstream)=0.
		 enddo
        else
         do j=1,nbrbat(upstream)-1
	      if (storcap(upstream) /= 0. .and. daystorcap(step,upstream) /= 0.) then
           if (maxlevel(upstream) >= elev_bat(j,upstream).AND.  &
               maxlevel(upstream) <= elev_bat(j+1,upstream)) THEN
            daystorcap(step,upstream)=vol_bat(j,upstream)+(maxlevel(upstream)-elev_bat  &
                (j,upstream))/(elev_bat(j+1,upstream)-elev_bat(j,upstream))*  &
                (vol_bat(j+1,upstream)-vol_bat(j,upstream))
           endif
	      else
	       daystorcap(step,upstream)=0.
	      endif
! Calculation of dead water volume after sediment deposition
	      if (damdead(upstream) /= 0. .and. daydamdead(step,upstream) /= 0.) then
           if (elevdead(upstream) >= elev_bat(j,upstream).AND.  &
               elevdead(upstream) <= elev_bat(j+1,upstream)) THEN
            daydamdead(step,upstream)=vol_bat(j,upstream)+(elevdead(upstream)-elev_bat  &
                (j,upstream))/(elev_bat(j+1,upstream)-elev_bat(j,upstream))*  &
                (vol_bat(j+1,upstream)-vol_bat(j,upstream))
           endif
	      else
	       daydamdead(step,upstream)=0.
	      endif
! Calculation of alert water volume after sediment deposition
 	      if (damalert(upstream) /= 0. .and. daydamalert(step,upstream) /= 0.) then
           if (elevalert(upstream) >= elev_bat(j,upstream).AND.  &
               elevalert(upstream) <= elev_bat(j+1,upstream)) THEN
            daydamalert(step,upstream)=vol_bat(j,upstream)+(elevalert(upstream)-elev_bat  &
                (j,upstream))/(elev_bat(j+1,upstream)-elev_bat(j,upstream))*  &
                (vol_bat(j+1,upstream)-vol_bat(j,upstream))
           endif
	      else
	       daydamalert(step,upstream)=0.
	      endif
! Calculation of maximum reservoir area after sediment deposition
 	      if (daymaxdamarea(step,upstream) /= 0.) then
           if (maxlevel(upstream) >= elev_bat(j,upstream).AND.  &
               maxlevel(upstream) <= elev_bat(j+1,upstream)) THEN
            daymaxdamarea(step,upstream)=area_bat(j,upstream)+(maxlevel(upstream)-elev_bat  &
                (j,upstream))/(elev_bat(j+1,upstream)-elev_bat(j,upstream))*  &
                (area_bat(j+1,upstream)-area_bat(j,upstream))
           endif
	      else
	       daymaxdamarea(step,upstream)=0.
	      endif
	     enddo
	    endif
	   endif
	  endif
     endif
    endif

!write(*,*)d,upstream,dayminlevel(step,upstream),daymaxdamarea(step,upstream)

! Print daily inflow and outflow discharges
    IF (reservoir_print == 0) THEN
     WRITE(subarea,*)id_subbas_extern(upstream)
	 IF (f_res_watbal) THEN
     OPEN(11,FILE=pfadn(1:pfadi)//'res_'//trim(adjustl(subarea))//'_watbal.out',STATUS='old',  &
		  POSITION='append')
	 WRITE(11,'(4I6,7f10.3,2f16.3)')id_subbas_extern(upstream),t,d,hour,qlateral(step,upstream),qinflow(step,upstream),  &
				qintake(step,upstream),overflow(step,upstream),qbottom(step,upstream),  &
				res_qout(step,upstream),damelevact(upstream),damareaact(upstream),volact(step,upstream)
     CLOSE(11)
	 ENDIF


! print storage losses due to reservoir sedimentation
	 IF (f_res_vollost) THEN
     OPEN(11,FILE=pfadn(1:pfadi)//'res_'//trim(adjustl(subarea))//'_vollost.out',STATUS='old',  &
		  POSITION='append')
	 WRITE(11,'(4I6,3f15.2)')id_subbas_extern(upstream),t,d,hour,daydamdead(step,upstream), &
			daydamalert(step,upstream),daystorcap(step,upstream)
     CLOSE(11)
	 ENDIF

! print temporal evolution of the stage-area-volume curve, when the initial curve is given by the user
     write(fmtstr,'(a,i0,a)')'(4I6,',nbrbat(upstream),'F15.2)'		!generate format string
     IF (nbrbat(upstream) /= 0) THEN
	  IF (f_res_cav) THEN
      OPEN(11,FILE=pfadn(1:pfadi)//'res_'//trim(adjustl(subarea))//'_cav.out',STATUS='old',  &
		    POSITION='append')
		WRITE(11,fmtstr)id_subbas_extern(upstream),t,d,hour,(elev_bat(j,upstream),j=1,nbrbat(upstream))
		WRITE(11,fmtstr)id_subbas_extern(upstream),t,d,hour,(area_bat(j,upstream),j=1,nbrbat(upstream))
		WRITE(11,fmtstr)id_subbas_extern(upstream),t,d,hour,(vol_bat(j,upstream),j=1,nbrbat(upstream))
      CLOSE(11)
	  ENDIF
	 ENDIF
    ELSE
	 daydamelevact(step,upstream)=damelevact(upstream)
	 daydamareaact(step,upstream)=damareaact(upstream)
	 DO j=1,nbrbat(upstream)
	   dayelev_bat(step,j,upstream)=elev_bat(j,upstream)
	   dayarea_bat(step,j,upstream)=area_bat(j,upstream)
	   dayvol_bat(step,j,upstream)=vol_bat(j,upstream)
	 ENDDO
	ENDIF


    volact(step,upstream)=volact(step,upstream)/1.e6
	daymaxdamarea(step,upstream)=daymaxdamarea(step,upstream)/1.e4
	daystorcap(step,upstream)=daystorcap(step,upstream)/1.e6
	daydamalert(step,upstream)=daydamalert(step,upstream)/1.e6
	daydamdead(step,upstream)=daydamdead(step,upstream)/1.e6
	damflow(upstream)=damflow(upstream)/(86400./nt)
	qoutlet(upstream)=qoutlet(upstream)/(86400./nt)
	withdrawal(upstream)=withdrawal(upstream)/(86400./nt)


	IF (step < dayyear*nt) THEN
	  volact(step+1,upstream)=volact(step,upstream)
	  daystorcap(step+1,upstream)=daystorcap(step,upstream)
	  daydamalert(step+1,upstream)=daydamalert(step,upstream)
	  daydamdead(step+1,upstream)=daydamdead(step,upstream)
	  daymaxdamarea(step+1,upstream)=daymaxdamarea(step,upstream)
	  dayminlevel(step+1,upstream)=dayminlevel(step,upstream)
	END IF

!write(*,'(2I4,3F15.4)')d,id_subbas_extern(upstream),dayminlevel(step,upstream),decstorcap(step,upstream)
!if (d==4)stop

! END of IF (storcap(upstream) > 0.) THEN
  ENDIF

! END of STATUS = 2
END IF

! -----------------------------------------------------------------------
IF (STATUS == 3) THEN

! Output files of Reservoir Modules
  IF (reservoir_print == 1) THEN
    DO i=1,subasin
      IF (storcap(i) /= 0. .and. t >= damyear(i)) THEN
        WRITE(subarea,*)id_subbas_extern(i)
		IF (f_res_watbal) THEN
	    OPEN(11,FILE=pfadn(1:pfadi)//'res_'//trim(adjustl(subarea))//'_watbal.out',STATUS='old',  &
			POSITION='append')
	    DO d=1,dayyear
	      DO ih=1,nt
		    hour=ih
            step=(d-1)*nt+hour
	        WRITE(11,'(4I6,7f10.3,2f16.3)')id_subbas_extern(i),t,d,hour,qlateral(step,i),qinflow(step,i),  &
				qintake(step,i),overflow(step,i),qbottom(step,i),res_qout(step,i),daydamelevact(step,i),daydamareaact(step,i),volact(step,i)*1.e6
		  ENDDO
		ENDDO
        CLOSE(11)
		ENDIF
		IF (f_res_vollost) THEN
		OPEN(11,FILE=pfadn(1:pfadi)//'res_'//trim(adjustl(subarea))//'_vollost.out',STATUS='old',  &
			POSITION='append')
	    DO d=1,dayyear
	      DO ih=1,nt
		    hour=ih
            step=(d-1)*nt+hour
			WRITE(11,'(4I6,3f15.2)')id_subbas_extern(i),t,d,hour,daydamdead(step,i)*1.e6, &
				daydamalert(step,i)*1.e6,daystorcap(step,i)*1.e6
		  ENDDO
		ENDDO
        CLOSE(11)
		ENDIF
	  ENDIF
      IF (storcap(i) /= 0. .and. t >= damyear(i) .and. nbrbat(i) /= 0) THEN
        WRITE(subarea,*)id_subbas_extern(i)
		IF (f_res_cav) THEN
        OPEN(11,FILE=pfadn(1:pfadi)//'res_'//trim(adjustl(subarea))//'_cav.out',STATUS='old',  &
		    POSITION='append')
	    write(fmtstr,'(a,i0,a)')'(4I6,',nbrbat(i),'F15.2)'		!generate format string
		DO d=1,dayyear
	      DO ih=1,nt
		    hour=ih
            step=(d-1)*nt+hour
			WRITE(11,fmtstr)id_subbas_extern(i),t,d,hour,(dayelev_bat(step,j,i),j=1,nbrbat(i))
			WRITE(11,fmtstr)id_subbas_extern(i),t,d,hour,(dayarea_bat(step,j,i),j=1,nbrbat(i))
			WRITE(11,fmtstr)id_subbas_extern(i),t,d,hour,(dayvol_bat(step,j,i),j=1,nbrbat(i))
!			WRITE(11,'(4I6,<nbrbat(i)>f15.2)')id_subbas_extern(i),t,d,hour,(dayelev_bat(step,j,i),j=1,nbrbat(i))
!			WRITE(11,'(4I6,<nbrbat(i)>f15.2)')id_subbas_extern(i),t,d,hour,(dayarea_bat(step,j,i),j=1,nbrbat(i))
!			WRITE(11,'(4I6,<nbrbat(i)>f15.2)')id_subbas_extern(i),t,d,hour,(dayvol_bat(step,j,i),j=1,nbrbat(i))
		  ENDDO
		ENDDO
        CLOSE(11)
		ENDIF
	  ENDIF
    ENDDO
  ENDIF

!Check simulation results for the Bengue catchment
!Grouping results on water and sediment balance of all reservoirs in the Bengue catchment into the reservoir classes
!***************************************************************************************************************
!  IF (reservoir_print == 1) THEN
!   IF (.not. doacud) THEN
!	DO d=1,dayyear
!	  DO ih=1,nt
!		hour=ih
!        step=(d-1)*nt+hour
!		inflow_class(1:6)=0.
!		outflow_class(1:6)=0.
!		retention_class(1:6)=0.
!		volume_class(1:6)=0.
!		sedinflow_class(1:6)=0.
!		sedoutflow_class(1:6)=0.
!		sedretention_class(1:6)=0.
!		sedimentation_class(1:6)=0.
!
!		a1=0
!		b1=0
!		c1=0
!		d1=0
!		e1=0
!		DO i=1,subasin
!		  IF (storcap(i) /= 0 .and. t >= damyear(i)) THEN
!	        IF (id_subbas_extern(i) /= 202 .and. id_subbas_extern(i) /= 203 .and. &
!			    id_subbas_extern(i) /= 212 .and. id_subbas_extern(i) /= 215 .and. &
!				id_subbas_extern(i) /= 219 .and. id_subbas_extern(i) /= 221 .and. &
!				id_subbas_extern(i) /= 223 .and. id_subbas_extern(i) /= 224 .and. &
!				id_subbas_extern(i) /= 227 .and. id_subbas_extern(i) /= 228 .and. &
!				id_subbas_extern(i) /= 229) THEN
!			  IF (storcap(i)*1.e6 <= 5000.) THEN
!			    a1=a1+1
!			    inflow_class(1)=inflow_class(1)+qinflow(step,i)*86400
!				outflow_class(1)=outflow_class(1)+res_qout(step,i)*86400
!				retention_class(1)=retention_class(1)+max((qinflow(step,i)-res_qout(step,i))*86400,0.)
!				volume_class(1)=volume_class(1)+volact(step,i)*1.e6
!			    if (dosediment) then
!			     sedinflow_class(1)=sedinflow_class(1)+sed_inflow(step,i)
!				 sedoutflow_class(1)=sedoutflow_class(1)+sed_outflow(step,i)
!				 sedretention_class(1)=sedretention_class(1)+sedimentation(step,i)
!				 sedimentation_class(1)=sedimentation_class(1)+daycumsed(step,i)
!			    endif
!			  ENDIF
!			  IF (storcap(i)*1.e6 > 5000. .and. storcap(i)*1.e6 < 25000.) THEN
!			    b1=b1+1
!			    inflow_class(2)=inflow_class(2)+qinflow(step,i)*86400
!				outflow_class(2)=outflow_class(2)+res_qout(step,i)*86400
!				retention_class(2)=retention_class(2)+max((qinflow(step,i)-res_qout(step,i))*86400,0.)
!				volume_class(2)=volume_class(2)+volact(step,i)*1.e6
!			    if (dosediment) then
!			     sedinflow_class(2)=sedinflow_class(2)+sed_inflow(step,i)
!				 sedoutflow_class(2)=sedoutflow_class(2)+sed_outflow(step,i)
!				 sedretention_class(2)=sedretention_class(2)+sedimentation(step,i)
!				 sedimentation_class(2)=sedimentation_class(2)+daycumsed(step,i)
!			    endif
!			  ENDIF
!			  IF (storcap(i)*1.e6 > 25000. .and. storcap(i)*1.e6 < 50000.) THEN
!			    c1=c1+1
!			    inflow_class(3)=inflow_class(3)+qinflow(step,i)*86400
!				outflow_class(3)=outflow_class(3)+res_qout(step,i)*86400
!				retention_class(3)=retention_class(3)+max((qinflow(step,i)-res_qout(step,i))*86400,0.)
!				volume_class(3)=volume_class(3)+volact(step,i)*1.e6
!			    if (dosediment) then
!			     sedinflow_class(3)=sedinflow_class(3)+sed_inflow(step,i)
!				 sedoutflow_class(3)=sedoutflow_class(3)+sed_outflow(step,i)
!				 sedretention_class(3)=sedretention_class(3)+sedimentation(step,i)
!				 sedimentation_class(3)=sedimentation_class(3)+daycumsed(step,i)
!			    endif
!			  ENDIF
!			  IF (storcap(i)*1.e6 > 50000. .and. storcap(i)*1.e6 < 100000.) THEN
!			    d1=d1+1
!			    inflow_class(4)=inflow_class(4)+qinflow(step,i)*86400
!				outflow_class(4)=outflow_class(4)+res_qout(step,i)*86400
!				retention_class(4)=retention_class(4)+max((qinflow(step,i)-res_qout(step,i))*86400,0.)
!				volume_class(4)=volume_class(4)+volact(step,i)*1.e6
!			    if (dosediment) then
!			     sedinflow_class(4)=sedinflow_class(4)+sed_inflow(step,i)
!				 sedoutflow_class(4)=sedoutflow_class(4)+sed_outflow(step,i)
!				 sedretention_class(4)=sedretention_class(4)+sedimentation(step,i)
!				 sedimentation_class(4)=sedimentation_class(4)+daycumsed(step,i)
!			    endif
!			  ENDIF
!			  IF (storcap(i)*1.e6 > 100000.) THEN
!			    e1=e1+1
!			    inflow_class(5)=inflow_class(5)+qinflow(step,i)*86400
!				outflow_class(5)=outflow_class(5)+res_qout(step,i)*86400
!				retention_class(5)=retention_class(5)+max((qinflow(step,i)-res_qout(step,i))*86400,0.)
!				volume_class(5)=volume_class(5)+volact(step,i)*1.e6
!			    if (dosediment) then
!			     sedinflow_class(5)=sedinflow_class(5)+sed_inflow(step,i)
!				 sedoutflow_class(5)=sedoutflow_class(5)+sed_outflow(step,i)
!				 sedretention_class(5)=sedretention_class(5)+sedimentation(step,i)
!				 sedimentation_class(5)=sedimentation_class(5)+daycumsed(step,i)
!			    endif
!			  ENDIF
!			ELSE
!			  inflow_class(6)=inflow_class(6)+qinflow(step,i)*86400
!			  outflow_class(6)=outflow_class(6)+res_qout(step,i)*86400
!			  retention_class(6)=retention_class(6)+max((qinflow(step,i)-res_qout(step,i))*86400,0.)
!			  volume_class(6)=volume_class(6)+volact(step,i)*1.e6
!			  if (dosediment) then
!			    sedinflow_class(6)=sedinflow_class(6)+sed_inflow(step,i)
!				sedoutflow_class(6)=sedoutflow_class(6)+sed_outflow(step,i)
!				sedretention_class(6)=sedretention_class(6)+sedimentation(step,i)
!				sedimentation_class(6)=sedimentation_class(6)+daycumsed(step,i)
!			  endif
!			ENDIF
!		  ENDIF
!		ENDDO
!        OPEN(11,FILE=pfadn(1:pfadi)//'res_watbal_lake.out',STATUS='old',  &
!				POSITION='append')
!		WRITE(11,'(3I6,24F15.3)')t,d,hour,(inflow_class(k1),k1=1,6),(outflow_class(k1),k1=1,6),(retention_class(k1),k1=1,6),(volume_class(k1),k1=1,6)
!		CLOSE(11)
!		if (dosediment) then
!         OPEN(11,FILE=pfadn(1:pfadi)//'res_sedbal_lake.out',STATUS='old',  &
!				POSITION='append')
!		 WRITE(11,'(3I6,24F15.3)')t,d,hour,(sedinflow_class(k1),k1=1,6),(sedoutflow_class(k1),k1=1,6),(sedretention_class(k1),k1=1,6),(sedimentation_class(k1),k1=1,6)
!		 CLOSE(11)
!		endif
!	  ENDDO
!	ENDDO
!   ELSE
!	DO d=1,dayyear
!	  DO ih=1,nt
!		hour=ih
!        step=(d-1)*nt+hour
!		inflow_class(6)=0.
!		outflow_class(6)=0.
!		retention_class(6)=0.
!		volume_class(6)=0.
!		sedinflow_class(6)=0.
!		sedoutflow_class(6)=0.
!		sedretention_class(6)=0.
!		sedimentation_class(6)=0.
!
!		DO i=1,subasin
!		  IF (storcap(i) /= 0 .and. t >= damyear(i)) THEN
!	        IF (id_subbas_extern(i) /= 29) THEN
!			  inflow_class(6)=inflow_class(6)+qinflow(step,i)*86400
!			  outflow_class(6)=outflow_class(6)+res_qout(step,i)*86400
!			  retention_class(6)=retention_class(6)+max((qinflow(step,i)-res_qout(step,i))*86400,0.)
!			  volume_class(6)=volume_class(6)+volact(step,i)*1.e6
!			  if (dosediment) then
!			    sedinflow_class(6)=sedinflow_class(6)+sed_inflow(step,i)
!			    sedoutflow_class(6)=sedoutflow_class(6)+sed_outflow(step,i)
!			    sedretention_class(6)=sedretention_class(6)+sedimentation(step,i)
!			    sedimentation_class(6)=sedimentation_class(6)+daycumsed(step,i)
!			  endif
!			ENDIF
!		  ENDIF
!		ENDDO
!        OPEN(11,FILE=pfadn(1:pfadi)//'res_watbal_lake.out',STATUS='old',  &
!				POSITION='append')
!		WRITE(11,'(3I6,24F15.3)')t,d,hour,inflow_class(6),outflow_class(6),retention_class(6),volume_class(6)
!		CLOSE(11)
!		if (dosediment) then
!         OPEN(11,FILE=pfadn(1:pfadi)//'res_sedbal_lake.out',STATUS='old',  &
!				POSITION='append')
!		 WRITE(11,'(3I6,24F15.3)')t,d,hour,sedinflow_class(6),sedoutflow_class(6),sedretention_class(6),sedimentation_class(6)
!		 CLOSE(11)
!		endif
!	  ENDDO
!	ENDDO
!   ENDIF
!  ENDIF
!***************************************************************************************************************

  IF (dosediment) CALL semres (status,upstream)

END IF



RETURN
END SUBROUTINE reservoir
