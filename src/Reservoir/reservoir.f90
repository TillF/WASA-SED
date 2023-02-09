SUBROUTINE reservoir (STATUS,upstream,res_h)

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

INTEGER :: i,id,dummy1,dummy2,dummy3,d4,s,p,q,h,istate,ka,ih,n
character(100) :: dummy_char
!Ge include j,nbrbat1,cont,upstream,downstream
INTEGER :: j,nbrbat1,flag_cav
REAL :: elevhelp,evaphelp2,volhelp !,elevhelp2,elevhelp3
REAL :: help,help1,help2,help3,evaphelp,areahelp,helpout,prechelp !,helpin,infhelp
!Ge actual storage capacity of large river reservoir in a certain year[10**6 m**3]
!REAL :: storcapact

CHARACTER(12) :: subarea

!REAL :: r_level0,r_level1,r_overflow,r_qbottom
!REAL :: r_precip,r_etp,r_qinflow

REAL :: dummy4,dummy5,dummy6
character(len=1000) :: fmtstr	!string for formatting file output

integer :: columnheader(1000) ! storing column heads of input files

!*****************************************************************************
!temporary variables used to test the cascade routing scheme of the lake module
!INTEGER :: a1,b1,c1,d1,e1,k1
!REAL :: inflow_class(6),outflow_class(6),retention_class(6),volume_class(6)
!REAL :: sedinflow_class(6),sedoutflow_class(6),sedretention_class(6),sedimentation_class(6)
! -----------------------------------------------------------------------
IF (STATUS == 0) THEN !begin of simulation

reservoir_check   = 0 !(0=simulation with all components; 1=simulation without hillslope and river modules)
reservoir_balance = 1 !(0=inflow and outflow discharges must be provided as input file; 1=only inflow discharges must be provided as input file)
reservoir_print   = 1 !(0=results printed at the end of the timestep; 1=results printed at the end of the simulated year)
f_intake_obs      = .false.

if (reservoir_check==0) reservoir_balance=1


  DO i=1,subasin
!Ge "volact" vector contains now only 365/366 values per annum
!Ge others reservoir parameters have to be inserted
    DO id=1,dayyear*nt
      volact(id,i) = 0. !actual stored volume in reservoir [m**3 and 10**6 m**3]
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
	  READ (11,*)  dummy3, upbasin(i),downbasin(i)
	END DO
	CLOSE (11)

!this relates the MAP IDs (id_subbas_extern(subasin)) to the sorted CODE IDs (id_subbas_intern(subasin)), i.e. upbasin and downbasin are now numbered according to internal ids
!so that in routing.dat only the MAP IDs have to be read in, first for the ID of upstream subasin
	!replace external with internal IDs
	DO i=1,subasin
	  !upstream basin referencing
	  ih=which1(id_subbas_extern == upbasin(i))

	  IF (ih==0) THEN
		  WRITE (*,'(A,I0,A)') 'ERROR: unknown upstream subbasin ID ', upbasin(i),' in routing.dat'
		  STOP
	  else
		upbasin(i)=ih
	  END IF

	  !downstream basin referencing
	  IF (downbasin(i) == 999 .OR. downbasin(i) == 9999) cycle 	!999 and 9999 mark outlet
	  ih=which1(id_subbas_extern == downbasin(i))

	  IF (ih==0) THEN
		  WRITE (*,'(A,I0,A)') 'ERROR: unknown downstream subbasin ID ', downbasin(i),' in routing.dat'
		  STOP
	  else
		downbasin(i)=ih
	  END IF
	END DO
  ENDIF

!! Read reservoir parameters
res_flag(:)=.false.
storcap(:)=0.

  OPEN(11,FILE=pfadp(1:pfadj)// 'Reservoir/reservoir.dat',IOSTAT=istate,STATUS='old')
	IF (istate/=0) THEN					!reservoir.dat not found
	  write(*,'(A)') "ERROR: ", pfadp(1:pfadj)// 'Reservoir/reservoir.dat was not found, please provide it.'
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
  n_reservoir=0
  res_index(:) = 0

  DO WHILE (.TRUE.)
	  READ (11,'(A)', IOSTAT=istate)fmtstr
	  if (istate/=0) exit !exit loop when no more line encountered
	  j=j+1
	  READ (fmtstr,*, IOSTAT=istate)dummy1

      i=which1(id_subbas_extern == dummy1) !find corresponding internal subbasin ID

	  IF (i==0) THEN
		  WRITE (*,'(A,I0,A)') 'WARNING: unknown upstream subbasin ID ', dummy1,' in reservoir.dat, ignored.'
		  cycle
      END IF 
      
	  READ (fmtstr,*,IOSTAT=istate) dummy1, minlevel(i), maxlevel(i),vol0(i),storcap(i), &
			damflow(i),damq_frac(i),withdrawal(i),damyear(i),maxdamarea(i), &
			damdead(i),damalert(i),dama(i),damb(i),qoutlet(i),fvol_bottom(i), &
			fvol_over(i),damc(i),damd(i),elevbottom(i)

      
	  IF (istate/=0) THEN
		  WRITE (*,'(A,i0,A,A)') 'ERROR: Format error in reservoir.dat, line ',j,':', fmtstr
		  STOP
	  END IF

	  ! check parameters
	  if(vol0(i) < 0 .and. vol0(i) /= -999.) then
        write(*,'(A,i3,A)') 'ERROR: Parameter vol0 in reservoir.dat is outside of plausible range (0 < vol0) for reservoir / subbasin id ', id_subbas_extern(i), '!'
        stop
	  end if
	  if(vol0(i) > storcap(i)) then
        write(*,'(A,i3,A)') 'ERROR: Parameter vol0 in reservoir.dat must not be larger than storage capacity for reservoir / subbasin id ', id_subbas_extern(i), '!'
        stop
      end if
	  if(damb(i) > .9 .or. damb(i) < 0.) then
        write(*,'(A,i3,A)') 'ERROR: Parameter damb in reservoir.dat is outside of plausible range (0 < damb <= 0.9) for reservoir / subbasin id ', id_subbas_extern(i), '!'
        stop
      end if
      if(damq_frac(i) > 1. .or. (damq_frac(i) < 0. .and. damq_frac(i) > -998.)) then
        write(*,'(A,i3,A)') 'WARNING: Parameter damq_frac in reservoir.dat is outside of plausible range (0 <= damq_frac <= 1 or eq. -999) for reservoir / subbasin id ', id_subbas_extern(i), '! During calibration this might make sense.'
	  end if
      if(fvol_bottom(i) > 1. .or. (fvol_bottom(i) < 0. .and. fvol_bottom(i) > -998.)) then
        write(*,'(A,i3,A)') 'WARNING: Parameter fvol_bottom in reservoir.dat is outside of plausible range (0 <= fvol_bottom <= 1 or eq. -999) for reservoir / subbasin id ', id_subbas_extern(i), '! During calibration this might make sense.'
	  end if
      if(fvol_over(i) > 1. .or. fvol_over(i) < 0.) then
        write(*,'(A,i3,A)') 'WARNING: Parameter fvol_over in reservoir.dat is outside of plausible range (0 <= fvol_over <= 1) for reservoir / subbasin id ', id_subbas_extern(i), '! During calibration this might make sense.'
	  end if
      if(damalert(i) < damdead(i)) then
        write(*,'(A,i3,A)') 'ERROR: Parameter damalert in reservoir.dat is less than damdead for reservoir / subbasin id ', id_subbas_extern(i), '!'
        stop
	  end if

	  ! set reservoir flag indicating that for subbasin i a reservoir exists and has been initialised
	  res_flag(i)  = .true.
	  n_reservoir  = n_reservoir+1 !count reservoirs
	  res_index(i) = n_reservoir  !note indices to corresponding entries in large arrays

!Ge "storcap" and "vol0" are read in 1000m**3 and after that they are converted into 10**6 m**3  &
!Ge damarea renamed to maxdamarea  &
!Ge "damdead" is read in 1000m**3 and after that it is converted into 10**6 m**3  &

      forma_factor(i)=1.e3*storcap(i)/((maxlevel(i)-minlevel(i))**3)
      IF (vol0(i) /= -999.) THEN
        vol0(i)=vol0(i)/1.e3 !convert in 10**6 m**3
      END IF
      storcap(i)=storcap(i)/1.e3 !convert in 10**6 m**3
      damdead(i)=damdead(i)/1.e3 !convert in 10**6 m**3
      damalert(i)=damalert(i)/1.e3 !convert in 10**6 m**3
  END DO
  CLOSE (11)

  
where (do_pre_outflow)
    res_flag(1:subasin) = .FALSE. !disable reservoirs in pre-specified basins, as their output is already given
    res_index(1:subasin) = 0
end where

!Anne & Till 2019 fix reservoir memory issue:
        !to decrease array size & only do calculations for subbasins with reservoir,  
        !moved all arrays with "subbasin" from allocate.h, line 400 ff to reservoir.f90 
        !and substituted "subbasin" by "n_reservoir";
        !plus inserted "res_index()", e.g.: dayarea_bat(step,j,i) changed to dayarea_bat(step,j,res_index(i))    
  
! allocate reservoir arrays, now that their required dimension is known 
 allocate( &
    
    corr_column_intakes(n_reservoir), &   
    reservoir_down(subasin), &    
    nbrbat(n_reservoir), &     
    dayexplot(n_reservoir,4), &
    operat_start(n_reservoir), &
    operat_stop(n_reservoir), &
    operat_elev(n_reservoir), &
    hmax(n_reservoir), &
               
    elevdead(n_reservoir), &     
    elevalert(n_reservoir), &
    damq_frac_season(n_reservoir,4), &	      
    precdam(366*nt,n_reservoir), &
    etdam(366*nt,n_reservoir), &     
    alpha_over(n_reservoir), &
    k_over(n_reservoir), &
    lakeret(366*nt,n_reservoir), &      
    qlateral(366*nt,n_reservoir), &
    qinflow(366*nt,n_reservoir), &
    overflow(366*nt,n_reservoir), &
    qintake(366*nt,n_reservoir), &
    qbottom(366*nt,n_reservoir), &   
    
    withdraw_out(366*nt,n_reservoir), &
     
    daystorcap(366*nt,n_reservoir), &
    daymaxdamarea(366*nt,n_reservoir), &    
    daydamdead(366*nt,n_reservoir), &
    daydamalert(366*nt,n_reservoir), &
    dayminlevel(366*nt,n_reservoir), &
    damelevact(n_reservoir), &    
    damvol0(n_reservoir), &
    damelev0(366*nt,n_reservoir), &
    damelev1(366*nt,n_reservoir), &
    resreach_vol(n_reservoir), &
     
    res_precip(366*nt,n_reservoir), &
    res_pet(366*nt,n_reservoir), &   
    res_qout(366*nt,n_reservoir), &      
    id_sec_extern(nxsection_res,n_reservoir), &
    nbrsec(n_reservoir), &
    npoints(nxsection_res,n_reservoir), &
	      
!    decvolact(366*nt,n_reservoir), &       !Anne variable seems to be unused
    decstorcap(366*nt,n_reservoir), &
!    decmaxdamarea(366*nt,n_reservoir), &   !Anne variable seems to be unused & is set to 0 in semres.f90
!    decdamdead(366*nt,n_reservoir), &      !Anne variable seems to be unused
!    decdamalert(366*nt,n_reservoir), &     !Anne variable seems to be unused
    manning_sec(nxsection_res,n_reservoir), &
    dist_sec(nxsection_res,n_reservoir), &
    x_sec0(npointsxsect,nxsection_res,n_reservoir), &
    y_sec0(npointsxsect,nxsection_res,n_reservoir), &
     
    sed_ret(366*nt,n_reservoir), &
    sed_overflow(366*nt,n_reservoir), &
    sed_intake(366*nt,n_reservoir), &  
    sed_bottom(366*nt,n_reservoir), &
    sed_qlateral(n_reservoir,n_sed_class), &
    sed_inflow(366*nt,n_reservoir), &
    sed_outflow(366*nt,n_reservoir), &    
    sedimentation(366*nt,n_reservoir), &    
    cum_sedimentation(n_reservoir), &      
    res_sediment_out(n_reservoir,n_sed_class), &
    frsediment_in(n_reservoir,n_sed_class), &
    frsediment_out(n_reservoir,n_sed_class), &
    dry_dens(n_reservoir), &
    factor_actlay(n_reservoir), &
    sed_flag(n_reservoir), &
    sed_routing_flag(n_reservoir), &   
    sedinflow_g(366*nt,n_reservoir,n_sed_class), &
    sedoutflow_g(366*nt,n_reservoir,n_sed_class), &
    damelev_mean(366*nt,n_reservoir), &
     
    x_minelev(nxsection_res,n_reservoir), &
    minelev_sec(nxsection_res,n_reservoir), &
    bedslope_sec(nxsection_res,n_reservoir), &
    area_sec(nxsection_res,n_reservoir), &
    resarea_sec(nxsection_res,n_reservoir), &
    resvol_sec(nxsection_res,n_reservoir), &
     
    resvol(n_reservoir), &
    topwidth_sec(nxsection_res,n_reservoir), &
    weight_sec(nxsection_res,n_reservoir), &
    discharge_sec(nxsection_res,n_reservoir), &
    depth_sec(nxsection_res,n_reservoir), &
     
    watelev_sec(nxsection_res,n_reservoir), &
    wetper_sec(nxsection_res,n_reservoir), &
    hydrad_sec(nxsection_res,n_reservoir), &
    meanvel_sec(nxsection_res,n_reservoir), &
    energslope_sec(nxsection_res,n_reservoir), &
     
    dynhead_sec(nxsection_res,n_reservoir), &
    tothead_sec(nxsection_res,n_reservoir), &
    headloss_sec(nxsection_res,n_reservoir), &
    locloss_sec(nxsection_res,n_reservoir), &
     
    calctothead_sec(nxsection_res,n_reservoir), &
    maxarea_sec(nxsection_res,n_reservoir), &
    maxelev_sec(nxsection_res,n_reservoir), &
    maxdepth_sec(nxsection_res,n_reservoir), &
    
    crdepth_sec(nxsection_res,n_reservoir), &
    crwatelev_sec(nxsection_res,n_reservoir), &
    crarea_sec(nxsection_res,n_reservoir), &
    crtopwidth_sec(nxsection_res,n_reservoir), &
    crwetper_sec(nxsection_res,n_reservoir), &     
    crslope_sec(nxsection_res,n_reservoir), &
    crvel_sec(nxsection_res,n_reservoir), &
    crhydrad_sec(nxsection_res,n_reservoir), &
    normalelev_sec(nxsection_res,n_reservoir), &
    normalarea_sec(nxsection_res,n_reservoir), &
       
    area_actlay(nxsection_res,n_reservoir), &
    area_toplay(nxsection_res,n_reservoir), &
    vol_actlay(nxsection_res,n_reservoir), &
    vol_toplay(nxsection_res,n_reservoir), &
                     
    erosion(nxsection_res,n_reservoir), &
    deposition(nxsection_res,n_reservoir), &
    retention(nxsection_res,n_reservoir), &      
    totalload(nxsection_res,n_reservoir), &	
    darea_sed(nxsection_res,n_reservoir), &
    dvol_sed(nxsection_res,n_reservoir), &

    frvol_actlay(n_sed_class,nxsection_res,n_reservoir), &
    totvol_actlay(nxsection_res,n_reservoir), &
    conc(nxsection_res,n_reservoir), &      
    area_sedim(nxsection_res,n_reservoir), &
    vol_sedim(nxsection_res,n_reservoir), &
    volbed0(n_reservoir), &
    length_plunge(n_reservoir), &
    
    cumlength_sec(nxsection_res,n_reservoir), &
    length_sec(nxsection_res,n_reservoir), &
    d50_actlay(nxsection_res,n_reservoir), &
    d90_actlay(nxsection_res,n_reservoir), &
    frsedinflow(366*nt,n_reservoir,n_sed_class), &
    frvol_actlay0(n_sed_class,nxsection_res,n_reservoir), &
    totvol_actlay0(nxsection_res,n_reservoir), &
       
    x_sec(npointsxsect,nxsection_res,n_reservoir), &
    y_sec(npointsxsect,nxsection_res,n_reservoir), &
    y_actlay(npointsxsect,nxsection_res,n_reservoir), &
    y_original(npointsxsect,nxsection_res,n_reservoir), &
    frac_actlay(n_sed_class,nxsection_res,n_reservoir), &
    frac_toplay(n_sed_class,nxsection_res,n_reservoir), &
    frac_comlay(n_sed_class,nxsection_res,n_reservoir), &
    frac_susp(n_sed_class,nxsection_res,n_reservoir), &
    partarea_actlay(npointsxsect,nxsection_res,n_reservoir), &
    partarea_toplay(npointsxsect,nxsection_res,n_reservoir), &     
     
    y_laststep(npointsxsect,nxsection_res,n_reservoir), &
    erosion_level(nxsection_res,n_reservoir), &
    pt1(nxsection_res,n_reservoir), &
    pt2(nxsection_res,n_reservoir), &
    pt3(nxsection_res,n_reservoir), &
    pt4(nxsection_res,n_reservoir), &
    pt_long0(n_reservoir), &
    pt_long(n_reservoir), &
    sideslope_pt1(nxsection_res,n_reservoir), & 
    sideslope_pt2(nxsection_res,n_reservoir), &
    slope_long(n_reservoir), &
     
    daydamelevact(366*nt,n_reservoir), & 
    daydamareaact(366*nt,n_reservoir), & 
    dayelev_bat(366*nt,nxsection_res,n_reservoir), & 
    dayarea_bat(366*nt,nxsection_res,n_reservoir), &
    dayvol_bat(366*nt,nxsection_res,n_reservoir), &
        
     STAT = istate)
    
    if (istate/=0) then
        write(*,'(A,i0,a)')'ERROR: Memory allocation error (',istate,') in reservoir-module, second allocation in reservoir.f90.'
        stop
    end if

    !initialisation
    corr_column_intakes = 0
    qinflow=0.
    qintake=0.
    overflow=0.
    
!Anne moved reservoir sediment variables from allocate_h to reservoir.f90
    
    if (dosediment) then
		    allocate( &

		      daydepth_sec(366*nt,nxsection_res,n_reservoir), &
		      daywatelev_sec(366*nt,nxsection_res,n_reservoir), &
		      dayarea_sec(366*nt,nxsection_res,n_reservoir), &
		      daytopwidth_sec(366*nt,nxsection_res,n_reservoir), &
		      dayenergslope_sec(366*nt,nxsection_res,n_reservoir), &
		      dayhydrad_sec(366*nt,nxsection_res,n_reservoir), &
		      daymeanvel_sec(366*nt,nxsection_res,n_reservoir), &
		      daydischarge_sec(366*nt,nxsection_res,n_reservoir), &
		      dayminelev_sec(366*nt,nxsection_res,n_reservoir), &
		      dayy_sec(366*nt,npointsxsect,nxsection_res,n_reservoir), &
		      daycumsed(366*nt,n_reservoir), &
		      dayfrsediment_out(366*nt,n_reservoir,n_sed_class), &
		    STAT = istate)
    
		    if (istate/=0) then
			    write(*,'(A,i0,a)')'ERROR: Memory allocation error (',istate,') in reservoir-module (sediments), second allocation in reservoir.f90.'
			    stop
		    end if
    
	    end if
    
    

! Check lateral inflow directly into the subbasins' reservoir
  latflow_res(1:subasin)=0

  OPEN(11,FILE=pfadp(1:pfadj)// 'Reservoir/lateral_inflow.dat', IOSTAT=istate,STATUS='old')
	IF (istate/=0) THEN					!lateral_inflow.dat not found
	  write(*,'(A)')'WARNING: '//pfadp(1:pfadj)// 'Reservoir/lateral_inflow.dat not found, using defaults'
    ELSE
	  READ(11,*, IOSTAT=istate)
	  READ(11,*, IOSTAT=istate)
	  h=2 !count lines
	  DO WHILE (.TRUE.)
	    READ(11,'(a)',IOSTAT=istate) fmtstr
	    if (istate /=0) exit
        h=h+1
        dummy1=GetNumberOfSubstrings(fmtstr) !Till: count number of fields/columns
        if (dummy1 /= 2) then    !incorrect number of fields in line
            write(*,'(a,i0,a)')'ERROR (lateral_inflow.dat): line ', h, ' contains more than 2 fields.'
            stop
        end if

        READ(fmtstr,*,IOSTAT=istate) dummy1, dummy2
	    if (istate/=0) then
            write(*,'(a,i0)')'ERROR (lateral_inflow.dat): Format error in line ', h
            stop
        end if

	    dummy3 = id_ext2int(dummy1, id_subbas_extern) !convert to internal ID
	    if (dummy3 == -1) then
            write(*,'(a,i0,a)')'WARNING (lateral_inflow.dat): unknown subbasin ID ', dummy1, ', skipped.'
            cycle
        end if
        latflow_res(dummy3)=1 !enable for this basin

        IF (dummy2 /= 999 .AND. dummy2 /= 9999) THEN !special marker, do not change
            d4 = id_ext2int(dummy2, id_subbas_extern) !convert to internal ID
        else
            d4 = dummy2
        end if

        if (d4 == -1) then
            write(*,'(a,i0,a)')'WARNING (lateral_inflow.dat): unknown subbasin ID ', dummy2, ', skipped.'
            latflow_res(dummy3)=0 !disable for this basin
            cycle
        end if
        reservoir_down(dummy3)=d4

!	    READ(11,*, IOSTAT=istate) dummy1					!read next line in file
!        DO i=1,subasin
!          IF (dummy1==id_subbas_extern(i)) THEN
!            latflow_res(i)=1
!		  ENDIF
!	    ENDDO

      END DO
	ENDIF
  CLOSE (11)



  OPEN(11,FILE=pfadp(1:pfadj)// 'Reservoir/operat_rule.dat', IOSTAT=istate, STATUS='old')
  IF (istate/=0) THEN					!operat_rule.dat not found
	write(*,'(A)')'WARNING: '//pfadp(1:pfadj)// 'Reservoir/operat_rule.dat not found, using defaults'
    DO i=1,subasin
	  if (damq_frac(i) == -999.) then
	    write(*,'(A)')'ERROR: operat_rule.dat must be given [or change the value of the parameter damq_frac in reservoir.dat]'
		stop
	  endif
	ENDDO
  ELSE
    READ(11,*)
    READ(11,*)
    DO i=1,subasin
	  IF (damq_frac(i) == -999.) READ (11,*)dummy1,(dayexplot(res_index(i),s),s=1,4),(damq_frac_season(res_index(i),s),s=1,4)
	  IF (damq_frac(i) /= -999.) dummy1=id_subbas_extern(i)
      IF (dummy1 /= id_subbas_extern(i)) THEN
        WRITE(*,'(A)') 'ERROR: Sub-basin-IDs in file operat_rule.dat must have the same ordering scheme as in hymo.dat'
        STOP
      END IF
	ENDDO
  ENDIF

  OPEN(11,FILE=pfadp(1:pfadj)// 'Reservoir/operat_bottom.dat', IOSTAT=istate,STATUS='old')
  IF (istate/=0) THEN					!operat_bottom.dat not found
	write(*,'(A)')'WARNING: '//pfadp(1:pfadj)// 'Reservoir/operat_bottom.dat not found, using defaults'
    DO i=1,subasin
	  if (fvol_bottom(i) == -999.) then
	    write(*,'(A)')'ERROR: operat_bottom.dat must be given [or change the value of the parameter damq_frac in reservoir.dat]'
		stop
	  endif
	ENDDO
  ELSE
    READ(11,*)
    READ(11,*)
    DO i=1,subasin
	  IF (fvol_bottom(i) == -999.) READ (11,*)dummy1,operat_start(res_index(i)),operat_stop(res_index(i)),operat_elev(res_index(i))
	  IF (fvol_bottom(i) /= -999.) dummy1=id_subbas_extern(i)
      IF (dummy1 /= id_subbas_extern(i)) THEN
        WRITE(*,'(A)') 'ERROR: Sub-basin-IDs in file operat_bottom.dat must have the same ordering scheme as in hymo.dat'
        STOP
      END IF
	ENDDO
  ENDIF

! initialise reading of intake.dat if it exists
  open(101,file=pfadp(1:pfadj)// 'Time_series/intake.dat', iostat=istate,status='old')
  if(istate /= 0) then
    write(*,'(A)') 'WARNING: '//pfadp(1:pfadj)// 'Time_series/intake.dat was not found, assumed zero.'
    close(101)
  else
    write(*,*) 'Reading controlled reservoir outflow through intake devices from file Time_series/intake.dat.'
    ! skip comment
    read(101,*)
    ! get reservoir (i.e. subbasin) ids from header line
    read(101,'(a)') fmtstr
    columnheader=0
    no_col_intake=GetNumberOfSubstrings(fmtstr)-2
    READ (fmtstr,*) dummy_char, dummy_char, (columnheader(i), i=1,no_col_intake)
    DO i=1,subasin
        if (do_pre_outflow(i)) cycle !if outflow has been pre-specified for this subbasin, just skip it
        DO j=1,size(columnheader)
            IF(columnheader(j) == id_subbas_extern(i)) THEN
                corr_column_intakes(res_index(i))= j    !for each subbasin, find position of corresponding column in input file
                f_intake_obs(i)=.true.
                exit
            END IF
        END DO
    END DO
    if(sum(corr_column_intakes) == 0) then
        write(*,*) '   File intake.dat does not contain relevant reservoir (i.e. subbasin) IDs! Running the model anyway...'
        close(101)
    else
        allocate(r_qintake(no_col_intake))
        r_qintake = 0.
        ! go to correct start line by analysing the date column
        call date_seek(101,tstart,mstart,dstart,'intake.dat')
    endif
  endif

!Ge stage-volume curves for each sub-basin
  nbrbat(1:n_reservoir)=0

  OPEN(11,FILE=pfadp(1:pfadj)// 'Reservoir/cav.dat', IOSTAT=istate,STATUS='old')
  IF (istate/=0) THEN					!cav.dat not found
	write(*,'(A)')'WARNING: '//pfadp(1:pfadj)// 'Reservoir/cav.dat not found, using defaults'
	flag_cav=0
  ELSE
	flag_cav=1
	READ(11,*)
 	READ(11,*)
    READ(11,*,IOSTAT=ka) dummy1,dummy2
    DO i=1,subasin
      IF (dummy1==id_subbas_extern(i)) THEN
         nbrbat(res_index(i))=dummy2
	  ENDIF
	ENDDO
    DO WHILE (ka==0)
	  READ(11,*);READ(11,*)
	  READ(11,*, IOSTAT=ka) dummy1,dummy2					!read next line in file
      DO i=1,subasin
        if (do_pre_outflow(i)) cycle !if outflow has been pre-specified for this subbasin, just skip it
        IF (dummy1==id_subbas_extern(i)) THEN
          nbrbat(res_index(i))=dummy2
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
		write(*,'(A,i0,a)')'ERROR: Memory allocation error (',istate,') in reservoir-module (rating curves too detailed).'
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
    if (res_index(i) /= 0) then !Anne inserted this line  
        nbrbat1=nbrbat(res_index(i))
        IF (nbrbat(res_index(i)) /= 0) THEN
          DO j=1,nbrbat1
            area_bat0(j,i)=area_bat0(j,i)*1000.
            vol_bat0(j,i)=vol_bat0(j,i)*1000.
          END DO
        END IF
     endif  !Anne 
  END DO

!Ge initialization of the stage-volume curves for each sub-basin (erosion/deposition process)
  DO i=1,subasin
    if (res_index(i) /= 0) then !Anne inserted this line  
        nbrbat1=nbrbat(res_index(i))
        IF (nbrbat(res_index(i)) /= 0) THEN
          DO j=1,nbrbat1
            elev_bat(j,i)=elev_bat0(j,i)
            area_bat(j,i)=area_bat0(j,i)
            vol_bat(j,i)=vol_bat0(j,i)
          END DO
        END IF
     endif  !Anne   
  END DO


  DO i=1,subasin
   IF (res_flag(i)) THEN
    nbrbat1=nbrbat(res_index(i))
    IF (nbrbat(res_index(i)) /= 0) THEN
      DO j=1,nbrbat(res_index(i))-1
        IF (damdead(i)*1.e6 >= vol_bat(j,i).AND.  &
            damdead(i)*1.e6 <= vol_bat(j+1,i)) THEN
          elevdead(res_index(i))=elev_bat(j,i)+((damdead(i)*1.e6)-vol_bat  &
                (j,i))/(vol_bat(j+1,i)-vol_bat(j,i))*  &
                (elev_bat(j+1,i)-elev_bat(j,i))
	    ENDIF
        IF (damalert(i)*1.e6 >= vol_bat(j,i).AND.  &
            damalert(i)*1.e6 <= vol_bat(j+1,i)) THEN
          elevalert(res_index(i))=elev_bat(j,i)+((damalert(i)*1.e6)-vol_bat  &
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
	  elevdead(res_index(i))=elevhelp+minlevel(i)
	  if (damalert(i) /= 0.) then
        elevhelp=((damalert(i)*1.e6)/forma_factor(i))**(1./3.)
	  else
        elevhelp=0.
	  endif
	  elevalert(res_index(i))=elevhelp+minlevel(i)
    ENDIF
   ENDIF
  END DO

!Ge Initialization of the parameters related to spillway overflow
  DO i=1,subasin
	IF (res_flag(i)) THEN
	  outflow_last(i)=0.
	  volume_last(i)=max(0., vol0(i) - storcap(i))
	  alpha_over(res_index(i))=1./(1.-damb(i))
	  k_over(res_index(i))=(dama(i)/alpha_over(res_index(i)))**alpha_over(res_index(i))
!write(*,*)id_subbas_extern(i),dama(i),damb(i),k_over(i),alpha_over(i),storcap(i)*1.e6
	  hmax(res_index(i))=((storcap(i)*1.e6)/k_over(res_index(i)))**(1./alpha_over(res_index(i)))
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
        IF (res_flag(i)) THEN
          WRITE(subarea,*)id_subbas_extern(i)
          OPEN(11,FILE=pfadn(1:pfadi)//'res_'//trim(adjustl(subarea))//'_watbal.out',STATUS='replace')
          WRITE(11,'(A)')'Subasin-ID'//char(9)//'year'//char(9)//'day'//char(9)//'hour'//char(9)//'qlateral(m**3/s)'//char(9)&
              //'inflow(m**3/s)'//char(9)//'evap(m**3)'//char(9)//'prec(m**3)'//char(9)//'intake(m**3/s)'//char(9)//&
              'overflow(m**3/s)'//char(9)//'qbottom(m**3/s)'//char(9)//'qout(m**3/s)'//char(9)//'withdrawal(m**3/s)'//char(9)//&
              'elevation(m)'//char(9)//'area(m**2)'//char(9)//'volume(m**3)'
          CLOSE(11)
        ENDIF
      ENDDO
  ENDIF

!Ge initialization of output files
  IF (f_res_vollost) then
      DO i=1,subasin
        IF (res_flag(i)) THEN
          WRITE(subarea,*)id_subbas_extern(i)
          OPEN(11,FILE=pfadn(1:pfadi)//'res_'//trim(adjustl(subarea))//'_vollost.out',STATUS='replace')
          WRITE(11,'(A)')'Subasin-ID, year, day, hour, deadvol(m**3), alertvol(m**3), storcap(m**3)'
          CLOSE(11)
        ENDIF
      ENDDO
  ENDIF

!Ge initialization of output files
  IF (f_res_cav) then
      DO i=1,subasin
        if (res_index(i) /= 0) then !Anne inserted this line  
           IF (nbrbat(res_index(i)) /= 0) THEN
            IF (res_flag(i)) THEN
              WRITE(subarea,*)id_subbas_extern(i)
              OPEN(11,FILE=pfadn(1:pfadi)//'res_'//trim(adjustl(subarea))//'_cav.out',STATUS='replace')
              WRITE(11,'(A)')'Subasin-ID, year, day, hour, 1st row: elev_bat(m), 2nd row: area_bat(m**2), 3rd row: vol_bat(m**3)'
            ENDIF
           ENDIF
         endif	 !Anne  
      ENDDO
  ENDIF

!**************************************************************************************
!Ge temporary output files to test the cascade routing scheme of the lake module
!  IF (.not. doacud) THEN
!   DO i=1,subasin
!    IF (storcap(i) /= 0.) THEN
!	  OPEN(11,FILE=pfadn(1:pfadi)//'res_watbal_lake.out',STATUS='replace')
!	    WRITE(11,'(A)')'Subasin-ID, year, day, hour, inflow_classes(m**3), outflow_classes(m**3), retention_classes(m**3), volume_classes(m**3)'
!      CLOSE(11)
!	ENDIF
!   ENDDO
!   DO i=1,subasin
!    IF (dosediment) THEN
!     IF (storcap(i) /= 0.) THEN
!	  OPEN(11,FILE=pfadn(1:pfadi)//'res_sedbal_lake.out',STATUS='replace')
!	    WRITE(11,'(A)')'Subasin-ID, year, day, hour, sedinflow_classes(ton), sedoutflow_classes(ton), sedretention_classes(ton), sedimentation_class'
!      CLOSE(11)
!	 ENDIF
!	ENDIF
!   ENDDO
!  ELSE
!   DO i=1,subasin
!    IF (storcap(i) /= 0.) THEN
!	  OPEN(11,FILE=pfadn(1:pfadi)//'res_watbal_lake.out',STATUS='replace')
!	    WRITE(11,'(A)')'Subasin-ID, year, day, hour, inflow_strateg(m**3), outflow_strateg(m**3), retention_strateg(m**3), volume_strateg(m**3)'
!      CLOSE(11)
!	ENDIF
!   ENDDO
!   DO i=1,subasin
!    IF (dosediment) THEN
!     IF (storcap(i) /= 0.) THEN
!	  OPEN(11,FILE=pfadn(1:pfadi)//'res_sedbal_lake.out',STATUS='replace')
!	    WRITE(11,'(A)')'Subasin-ID, year, day, hour, sedinflow_strateg(ton), sedoutflow_strateg(ton), sedretention_strateg(ton), sedimentation_class'
!      CLOSE(11)
!	 ENDIF
!	ENDIF
!   ENDDO
!  ENDIF
!**************************************************************************************


END IF


! -----------------------------------------------------------------------
IF (STATUS == 1) THEN !beginning of a new simulation year
! initialize ...

!Ge initialization of parameters
    overflow=0.
    qintake=0.
    qinflow=0.
    qbottom=0.
    qlateral=0.
 !   volact(1,:)=0. !Till: already initialized in model_state_io
    damareaact=0.

    etdam = 0.
    precdam = 0.
    res_qout = 0.
    withdraw_out = 0.
    damelevact = 0.
  

  !IF (t > tstart) THEN !Andreas
  DO i=1,subasin !Andreas
   IF (res_flag(i)) THEN
     IF (t > damyear(i) .AND. t > tstart) THEN  !continuing simulations
       volact       (1,i)            = volact       (daylastyear*nt,i)
       daystorcap   (1,res_index(i)) = daystorcap   (daylastyear*nt,res_index(i))
       daydamalert  (1,res_index(i)) = daydamalert  (daylastyear*nt,res_index(i))
       daydamdead   (1,res_index(i)) = daydamdead   (daylastyear*nt,res_index(i))
       daymaxdamarea(1,res_index(i)) = daymaxdamarea(daylastyear*nt,res_index(i))
	   dayminlevel  (1,res_index(i)) = dayminlevel  (daylastyear*nt,res_index(i))

     ELSE IF (t == damyear(i) .OR. (t > damyear(i) .AND. t == tstart)) THEN  !Andreas
       if (volact(1,i)== -1) then !only initialize those that have not been initialized by loading from stat-file
           IF (vol0(i) /= -999. .and. vol0(i) /= -9999.) THEN 
             volact(1,i)=vol0(i) !Till: initial volume [1e6 m^3]
           ELSE
             volact(1,i)=storcap(i)/5. !Till: initial volume as 20% of storage cap [1e6 m^3]
           END IF
       end if    
       daystorcap(1,res_index(i))=storcap(i)
       daydamalert(1,res_index(i))=damalert(i)
       daydamdead(1,res_index(i))=damdead(i)
       daymaxdamarea(1,res_index(i))=maxdamarea(i)
	   dayminlevel(1,res_index(i))=minlevel(i)
     ELSE  !no reservoirs or not-yet built
       volact(1,i)=0.
       daystorcap(1,res_index(i))=0.
       daydamalert(1,res_index(i))=0.
       daydamdead(1,res_index(i))=0.
       daymaxdamarea(1,res_index(i))=0.
       dayminlevel(1,res_index(i))=0.
	 ENDIF
   END IF
  END DO !Andreas


!Ge water availability approach for reservoirs has to be included
  avail_ac(:,:)=0. !water availability
  avail_all(:,:)=0.   !water availability
!  damex(:,:)=0. tp not used

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
    IF (res_flag(i) .and. t >= damyear(i)) THEN
     nbrbat1=nbrbat(res_index(i))
     IF (nbrbat(res_index(i)) /= 0) THEN
      DO j=1,nbrbat(res_index(i))-1
        IF (volact(1,i)*1.e6 >= vol_bat(j,i).AND.  &
            volact(1,i)*1.e6 <= vol_bat(j+1,i)) THEN
          damareaact(i)=area_bat(j,i)+((volact(1,i)*1.e6)-vol_bat  &
                (j,i))/(vol_bat(j+1,i)-vol_bat(j,i))*  &
                (area_bat(j+1,i)-area_bat(j,i))
	    ENDIF
	  ENDDO
     ELSE
	  if (volact(1,i) /= 0.) then
		areahelp=dama(i)*((volact(1,i)*1.e6)**damb(i)) !Till: estimate surface area [m^2]
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
!   tobias: file Reservoir/inflow_*.dat is not documented: Maybe the following commented code section can be removed?!
!   IF (reservoir_balance == 0) THEN
!    DO i=1,subasin
!     IF (t >= damyear(i)) THEN !Andreas
!	  IF (storcap(i) > 0.) THEN
!        WRITE(subarea,*)id_subbas_extern(i)
!		OPEN(11,FILE=pfadp(1:pfadj)//'Reservoir/inflow_'//trim(adjustl(subarea))//'.dat', &
!			STATUS='unknown')
!          read(11,*)
!          cont=(dtot-dayyear)*nt
!          DO id=1,cont
!            READ(11,*)
!          END DO
!          DO id=1,dayyear*nt
!            READ(11,*) dummy1,dummy2,r_precip,r_etp,r_qinflow, &
!				r_overflow,r_qintake,r_qbottom,r_level0,r_level1
!	        damelev0(id,i)=r_level0
!	        damelev1(id,i)=r_level1
!            overflow(id,i)=r_overflow
!            qintake(id,i)=r_qintake
!	        qbottom(id,i)=r_qbottom
!            res_precip(id,i)=r_precip
!	        res_pet(id,i)=r_etp
!	        qinflow(id,i)=r_qinflow
!	      ENDDO
!	    CLOSE(11)
!        IF (nbrbat(i) /= 0) THEN
!	      nbrbat1=nbrbat(i)
!          DO id=1,dayyear
!            IF (damelev0(id,i) > elev_bat0(nbrbat1,i) .or. damelev1(id,i) > elev_bat0(nbrbat1,i) ) THEN
!              WRITE(*,*)'ERROR subasin ',id_subbas_extern(i),' year ',t,' day ',id, &
!			     	'GIVEN VALUE OF DAILY RESERVOIR LEVEL IS GREATER THAN THE MAXIMUM RESERVOIR ELEVATION AT THE STAGE-AREA-VOLUME CURVE (FILE: cav.dat)'
!			  STOP
!			ELSE IF (damelev0(id,i) < elev_bat0(1,i) .or. damelev1(id,i) < elev_bat0(1,i)) THEN
!			  WRITE(*,*)'ERROR subasin ',id_subbas_extern(i),' year ',t,' day ',id,  &
!					'GIVEN VALUE OF DAILY RESERVOIR LEVEL IS LESS THAN THE MINIMUM ELEVATION AT THE STAGE-AREA-VOLUME CURVE (FILE: cav.dat)'
!			  STOP
!			END IF
!		  ENDDO
!		ENDIF
!	  ENDIF
!	 ENDIF !Andreas
!	ENDDO
!   ENDIF
!   IF (reservoir_balance == 1) THEN
!    DO i=1,subasin
!     IF (t >= damyear(i)) THEN !Andreas
!	  IF (storcap(i) > 0.) THEN
!        WRITE(subarea,*)id_subbas_extern(i)
!		OPEN(11,FILE=pfadp(1:pfadj)//'Reservoir/inflow_'//trim(adjustl(subarea))//'.dat', &
!			STATUS='unknown')
!          read(11,*)
!          cont=(dtot-dayyear)*nt
!          DO id=1,cont
!            READ(11,*)
!          END DO
!          DO id=1,dayyear*nt
!            READ(11,*) dummy1,dummy2,r_precip,r_etp,r_qinflow
!            res_precip(id,i)=r_precip
!	        res_pet(id,i)=r_etp
!	        qinflow(id,i)=r_qinflow
!	      ENDDO
!	    CLOSE(11)
!	  ENDIF
!	 ENDIF !Andreas
!	ENDDO
!   ENDIF
 ELSE ! tobias: reservoir_check == 0 which is actually the only possible value for reservoir_check as it is hard-coded
    DO i=1,subasin
    IF (t >= damyear(i)) THEN !Andreas
     dummy1=0
     if (res_index(i) /= 0) then !Anne inserted this line
            DO id=1,dayyear
	              if (river_transport.ne.1) then
                   DO h=1,nt
	                dummy1=dummy1+1
	                res_pet(dummy1,res_index(i))=pet(id,i)/nt
                   ENDDO
	              else
	               res_pet(id,res_index(i))=pet(id,i)
                  endif             
            ENDDO
     endif    !Anne
    ENDIF !Andreas
    ENDDO
 
   DO i=1,subasin
    IF (t >= damyear(i)) THEN !Andreas
	 IF (dohour) THEN
	  IF (river_transport.ne.1) THEN
          if (res_index(i) /= 0) then !Anne inserted this line
               DO id=1,dayyear*nt
	            res_precip(id,res_index(i))=preciph(id,i)
               ENDDO
          endif   !Anne  
	  ENDIF
	  IF (river_transport.eq.1) THEN
       DO id=1,dayyear
	    dummy4=0.
	    DO h=1,nt
		  dummy4=dummy4+preciph(id,i)
	    ENDDO
	    res_precip(id,res_index(i))=dummy4
	   ENDDO
	  ENDIF
     ELSE
         if (res_index(i) /= 0) then !Anne inserted this line
              DO id=1,dayyear
	            res_precip(id,res_index(i))=precip(id,i)
              ENDDO
          endif   !Anne     
	 ENDIF
    ENDIF !Andreas
   ENDDO
 ENDIF

! begin block Andreas; changed by tobias
! Read daily data on (measured) regulated reservoir outflow to be considered in reservoir water balance
! data given in m**3/s, convert to m**3/d
!  IF (reservoir_balance == 1) THEN
!   DO i=1,subasin
!    IF (t >= damyear(i)) THEN
!     IF (damq_frac(i) == -888.) THEN
!      WRITE(subarea,*)id_subbas_extern(i)
!      OPEN(11,FILE=pfadp(1:pfadj)//'Reservoir/intake_'//trim(adjustl(subarea))//'.dat', &
!			STATUS='unknown')
!       read(11,*)
!	   read(11,*)
!       cont=(dtot-dayyear)*nt
!       DO id=1,cont
!         READ(11,*)
!       END DO
!       DO id=1,dayyear*nt
!         READ(11,*) dummy1,r_qintake
!		 IF (qintake(id,i) /= -999.) qintake(id,i)=r_qintake*(86400./nt)
!		 IF (qintake(id,i) == -999.) qintake(id,i)=damflow(i)*damq_frac(i)
!       ENDDO
!      CLOSE(11)
!     ENDIF
!    ENDIF
!   ENDDO
!  ENDIF
  if(any(f_intake_obs)) then
    do id=1,dayyear*nt
        ! read data for current day from intake.dat
        read(101, *, iostat=istate) dummy1,dummy2, r_qintake !ii: read only columns that are actually needed
        IF (istate/=0) THEN
            write(*,'(A)')'ERROR: Premature end of file intake.dat.'
            stop
        END IF
        ! distribute values to qintake variable
        ! NOTE: missing observations are treated later
        do i=1,subasin
            if( f_intake_obs(i) .and. (t >= damyear(i)) ) then
                qintake(id,res_index(i)) = r_qintake(corr_column_intakes(res_index(i)))*(86400./nt)
            endif
        enddo
    enddo
  endif
! end block Andreas; changed by tobias


  IF (dosediment) THEN
    CALL semres (status,upstream)
  END IF




END IF

! -----------------------------------------------------------------------
IF (STATUS == 2) THEN !regular call during timestep


! simulation timestep
  if (river_transport.ne.1) then
    hour=res_h
    step=(d-1)*nt+hour
  else
    hour=1
	step=d
  endif

! Computation of reservoir water balance
  IF (res_flag(upstream) .and. t >= damyear(upstream)) THEN
    IF (reservoir_balance == 1) THEN
      IF (reservoir_check == 1) THEN
	    qinflow(step,res_index(upstream))=qinflow(step,res_index(upstream))*(86400./nt) !Till: convert m3/s to m3
	  ELSE
        if (river_transport.eq.1)then ! old routing
	      qinflow(step,res_index(upstream))=qout(step,upstream)*(86400./nt)
	    else ! new routing
          qinflow(step,res_index(upstream))=(r_qout(2,upstream)+qlateral(step,res_index(upstream)))*(86400./nt)
	    endif
	  ENDIF

!qinflow(step,upstream)=20*86400.
!write(*,*)upstream,id_subbas_extern(upstream)
!if(step<=10)qinflow(step,upstream)=storcap(upstream)*1.e6
!if(step>20)qinflow(step,upstream)=storcap(upstream)*1.e6
!if (step == 16 .and. t==2002) read(*,*)

      volact(step,upstream)=volact(step,upstream)*1.e6 !convert to m^3
      daymaxdamarea(step,res_index(upstream))=daymaxdamarea(step,res_index(upstream))*1.e4
      daystorcap(step,res_index(upstream))=daystorcap(step,res_index(upstream))*1.e6
      daydamalert(step,res_index(upstream))=daydamalert(step,res_index(upstream))*1.e6
      daydamdead(step,res_index(upstream))=daydamdead(step,res_index(upstream))*1.e6
      damflow(upstream)=damflow(upstream)*(86400./nt)
      damvol0(res_index(upstream))=volact(step,upstream)
	  qoutlet(upstream)=qoutlet(upstream)*(86400./nt)
	  withdrawal(upstream)=withdrawal(upstream)*(86400./nt)

!write(*,'(2I4,3F15.4)')step,id_subbas_extern(upstream),volact(step,upstream)


! 1) Calculation of the actual reservoir volume after water inflow
	  help3=volact(step,upstream)+qinflow(step,res_index(upstream)) !Till: m^3
	  help2=volact(step,upstream) !Till: m^3
!write(*,'(2I4,4F15.4)')step,id_subbas_extern(upstream),help3,volact(step,upstream),daystorcap(step,upstream),overflow(step,upstream)
      IF (help3 > daystorcap(step,res_index(upstream))) THEN
!write(*,'(2I4,4F15.4)')step,id_subbas_extern(upstream),volact(step,upstream),help3,help,daystorcap(step,upstream)
!write(*,'(2I4,4F15.4)')step,id_subbas_extern(upstream),fvol_over(upstream)*daystorcap(step,upstream),help3,daystorcap(step,upstream)
		IF (fvol_over(upstream)==1) THEN
          help=daystorcap(step,res_index(upstream))
!write(*,'(2I4,3F15.4)')step,id_subbas_extern(upstream),overflow(step,upstream)/86400.
!write(*,'(2I4,4F15.4)')step,id_subbas_extern(upstream),volact(step,upstream),help3,help,daystorcap(step,upstream)
		ELSE
		  help=daystorcap(step,res_index(upstream))
		  overflow(step,res_index(upstream))=help3-daystorcap(step,res_index(upstream)) !Till: spillover [m^3]
        END IF
	    lakeret(step,res_index(upstream))=max(0.,help-volact(step,upstream))
        volact(step,upstream)=help
	  ELSE
	    lakeret(step,res_index(upstream))=qinflow(step,res_index(upstream))
        volact(step,upstream)=help3
      END IF
!write(*,'(2I4,4F15.4)')step,id_subbas_extern(upstream),help,volact(step,upstream),daystorcap(step,upstream),overflow(step,upstream)
!if (d==31)stop
!write(*,'(2I4,4F15.3)')d,id_subbas_extern(upstream),fvol_over(upstream),daystorcap(step,upstream),volact(step,upstream)+qinflow(step,upstream),overflow(step,upstream)
!write(*,'(2I4,4F15.3)')d,id_subbas_extern(upstream),volact(step,upstream)

!write(*,'(2I4,4F15.4)')step,id_subbas_extern(upstream),volact(step,upstream),help3,help,daystorcap(step,upstream)

! 2) Substract evaporation out the reservoir
! Calculation of the actual reservoir area by interpolation (using the stage-area-volume curve)
      IF (nbrbat(res_index(upstream)) /= 0) THEN
        DO j=1,nbrbat(res_index(upstream))-1
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
		elevhelp=elevhelp+dayminlevel(step,res_index(upstream))
        areahelp=dama(upstream)*(volact(step,upstream)**damb(upstream))
      END IF
      damelevact(res_index(upstream))=elevhelp
      damareaact(upstream)=areahelp
!write(*,'(2I4,4F15.3)')d,id_subbas_extern(upstream),damelevact(upstream),damareaact(upstream)


! Calculation of increase or decrease of the reservoir level due to precip and evap.
! (using the stage-area-volume curve)
      evaphelp2=0.
      IF (nbrbat(res_index(upstream)) /= 0) THEN
        damelevact(res_index(upstream))=damelevact(res_index(upstream))+  &
            (res_precip(step,res_index(upstream))-res_pet(step,res_index(upstream)))/1000.
        IF (damelevact(res_index(upstream)) < dayminlevel(step,res_index(upstream))) THEN
          evaphelp2=dayminlevel(step,res_index(upstream))-damelevact(res_index(upstream))
          damelevact(res_index(upstream))=dayminlevel(step,res_index(upstream))
        END IF

! Check overflow due to precipitation
!George        IF (damelevact(upstream) > maxlevel(upstream)) THEN
!George          overflow(step,upstream)=overflow(step,upstream)+((damelevact(upstream)-maxlevel(upstream))*  &
!George              daymaxdamarea(step,upstream))
!George          damelevact(upstream)=maxlevel(upstream)
!George        END IF
        volhelp = -1. !flag as "not computed"
        DO j=1,nbrbat(res_index(upstream))-1 !iterate through points of CAV
          IF ((damelevact(res_index(upstream)) >= elev_bat(j,upstream).AND.  &
                damelevact(res_index(upstream)) <= elev_bat(j+1,upstream)) .OR. &
            (j == nbrbat(res_index(upstream))-1) & ! (when water stage is higher than max stage in CAV extrapolate CAV-curve
          ) THEN
            volhelp=vol_bat(j,upstream)+(damelevact(res_index(upstream))-elev_bat  &
                (j,upstream))/(elev_bat(j+1,upstream)-elev_bat(j,upstream))*  &
                (vol_bat(j+1,upstream)-vol_bat(j,upstream))
            areahelp=area_bat(j,upstream)+(damelevact(res_index(upstream))-elev_bat  &
                (j,upstream))/(elev_bat(j+1,upstream)-elev_bat(j,upstream))*  &
                (area_bat(j+1,upstream)-area_bat(j,upstream))
            exit !correct point of CAV-found, no more searching
          END IF
      END DO

      if (damelevact(res_index(upstream)) > elev_bat(nbrbat(res_index(upstream)),upstream)) then
          write(*,"(A,i0,a)")"WARNING: Water stage of reservoir ",id_subbas_extern(upstream)," exceeds CAV-curve. Curve extrapolated."
      end if


! Calculation of evaporation and precipitation using the truncated cone volume (m3)
! (using the morphologic parameter alpha)
        evaphelp=(areahelp+SQRT(areahelp*damareaact(upstream))+  &
            damareaact(upstream))*res_pet(step,res_index(upstream))/1000.*1./3.
        prechelp=(areahelp+SQRT(areahelp*damareaact(upstream))+  &
            damareaact(upstream))*res_precip(step,res_index(upstream))/1000.*1./3.
!        infhelp=0. tp TODO not used=!
        volact(step,upstream)=volhelp

      ELSE
        evaphelp=MIN(volact(step,upstream),(res_pet(step,res_index(upstream))/1000.)*areahelp) !??why zero?
        prechelp=(res_precip(step,res_index(upstream))/1000.)*areahelp
!        infhelp=0. tp TODO not used=!
        volact(step,upstream)=volact(step,upstream)+ (prechelp-evaphelp) ! -infhelp) tp TODO not used

! Check overflow due to precipitation
!        IF (volact(step,upstream) > daystorcap(step,upstream)) THEN
!          overflow(step,upstream)=overflow(step,upstream)+volact(step,upstream)-daystorcap(step,upstream)
!          volact(step,upstream)=daystorcap(step,upstream)
!        END IF
      END IF
!write(*,'(2I4,4F15.3)')d,id_subbas_extern(upstream),overflow(step,upstream),volact(step,upstream)
!write(*,'(2I4,3F15.4)')step,id_subbas_extern(upstream),volact(step,upstream)

!tp store for output in reservoir balance in m**3
       etdam(step,res_index(upstream))=evaphelp
       precdam(step,res_index(upstream))=prechelp

! 3) reservoir evaporation in mm, tp: TODO never used!
!      IF (nbrbat(upstream) /= 0) THEN
!        evapdam(step,upstream)=MAX(0.,res_pet(step,upstream)-(evaphelp2*1.e3))
!      ELSE
!        evapdam(step,upstream)=evaphelp/areahelp*1.e3
!      END IF
!!     reservoir evaporation in Mio m**3
!      etdam(step,upstream)=evaphelp/1.e6
!      infdam(step,upstream)=infhelp/1E6


! 4) Check flow through the reservoir
! 4a) Water intake device
!     - dam outflow is fraction of Q90
!     - according to J.C.Araújo this fraction is estimated to be 90%,
!       and 80% for larger, strategic dams (e.g. Oros)
!     - if alert volume is reached (if given for dam), outflow is reduced
!     - no outflow if storage is below dead volume (if given for dam)
!     - for non-strategic dams outflow is demand of regular (mean precip) years
!       (not yet implemented)

! use measured intake values if available
      if(f_intake_obs(upstream) .and. qintake(step,res_index(upstream)) > -0.5) then
        helpout=qintake(step,res_index(upstream))
      else
! Calculation of the maximum controlled outflow discharge using a factor defined in the reservoir.dat
          IF (damq_frac(upstream) >= 0.0) THEN !Andreas
            helpout=damflow(upstream)*damq_frac(upstream)
! Calculation of the maximum controlled outflow discharge using a operation regime as provided in the operat_rule.dat
          ELSE
            dummy4=0.
            do s=1,3
              IF (dayoutsim+d >= dayexplot(res_index(upstream),s) .and. &
                  dayoutsim+d < dayexplot(res_index(upstream),s+1)) dummy4=damq_frac_season(res_index(upstream),s)
            enddo
            IF (dayoutsim+d < dayexplot(res_index(upstream),1) .or. &
                dayoutsim+d >= dayexplot(res_index(upstream),4)) dummy4=damq_frac_season(res_index(upstream),4)
            helpout=damflow(upstream)*dummy4
          ENDIF
!write(*,*)upstream,dayoutsim+d,(dayexplot(upstream,s),s=1,4)
!write(*,*)d,upstream,helpout,damflow(upstream),dummy4
!write(*,*)step,id_subbas_extern(upstream),helpout/86400.,dummy4,damflow(upstream)/86400.
        endif ! measured or generic intake values

! Check water availability
      IF (daydamalert(step,res_index(upstream)) > daydamdead(step,res_index(upstream))) THEN
        IF (volact(step,upstream) < daydamalert(step,res_index(upstream))) THEN
          IF (volact(step,upstream) > daydamdead(step,res_index(upstream))) THEN
            helpout=helpout*(volact(step,upstream)-daydamdead(step,res_index(upstream)))/  &
                (daydamalert(step,res_index(upstream))-daydamdead(step,res_index(upstream)))
          END IF
        END IF
      ENDIF
      IF (volact(step,upstream) < (daydamdead(step,res_index(upstream))+helpout))THEN
        helpout=volact(step,upstream)-daydamdead(step,res_index(upstream))
      END IF

      IF (volact(step,upstream) <= daydamdead(step,res_index(upstream))) helpout=0.
      IF (elevdead(res_index(upstream)) <= dayminlevel(step,res_index(upstream))) helpout=0.
 !write(*,*)step,id_subbas_extern(upstream),helpout/86400.
      qintake(step,res_index(upstream))=helpout
      volact(step,upstream)=MAX(0.,volact(step,upstream)-qintake(step,res_index(upstream)))

!write(*,'(2I4,3F15.4)')step,id_subbas_extern(upstream),volact(step,upstream)


! 4b) Bottom outlets
      IF (nbrbat(res_index(upstream)) /= 0) THEN
        DO j=1,nbrbat(res_index(upstream))-1
          IF (operat_elev(res_index(upstream)) >= elev_bat(j,upstream).AND.  &
                operat_elev(res_index(upstream)) <= elev_bat(j+1,upstream)) THEN
            volhelp=vol_bat(j,upstream)+(operat_elev(res_index(upstream))-elev_bat  &
                (j,upstream))/(elev_bat(j+1,upstream)-elev_bat(j,upstream))*  &
                (vol_bat(j+1,upstream)-vol_bat(j,upstream))
         END IF
        END DO
      ELSE
        elevhelp=operat_elev(res_index(upstream))-dayminlevel(step,res_index(upstream))
        volhelp=forma_factor(upstream)*(elevhelp**3)
      END IF

	  IF (fvol_bottom(upstream) /= -999.) THEN
	   helpout=qoutlet(upstream)
       IF (volact(step,upstream) > daydamdead(step,res_index(upstream))) THEN
        IF (volact(step,upstream) < daystorcap(step,res_index(upstream))) THEN
          IF (volact(step,upstream) > fvol_bottom(upstream)*daystorcap(step,res_index(upstream))) THEN
            helpout=min(helpout,helpout*(volact(step,upstream)-(fvol_bottom(upstream)*daystorcap(step,res_index(upstream))))/ &
		             ((1.-fvol_bottom(upstream))*daystorcap(step,res_index(upstream))))
		  ELSE
		    helpout=0.
		  ENDIF
        ENDIF
       ELSE
		helpout=0.
       END IF
	  ELSE
	   IF (step >= operat_start(res_index(upstream)) .and. step <= operat_stop(res_index(upstream))) THEN
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

      qbottom(step,res_index(upstream))=helpout
      volact(step,upstream)=MAX(0.,volact(step,upstream)-qbottom(step,res_index(upstream)))
!write(*,'(2I4,3F15.4)')step,id_subbas_extern(upstream),volact(step,upstream)


! 4c) Withdrawal water volume to supply the water use sectors
      IF (volact(step,upstream) > .05*daystorcap(step,res_index(upstream))) THEN
        withdraw_out(step,res_index(upstream)) = withdrawal(upstream) !ii: first, we should check that enough water is really available in teh reservoir
        volact(step,upstream)=MAX(0.,volact(step,upstream)-withdrawal(upstream))
	  ELSE
        withdraw_out(step,res_index(upstream)) = 0.
      ENDIF

!write(*,'(2I4,3F15.4)')step,id_subbas_extern(upstream),volact(step,upstream)


! 5) Calculation of the overflow discharges from the reservoir
      IF (help3 > daystorcap(step,res_index(upstream)) .and. fvol_over(upstream) == 1) THEN
        ! volume in/decrease relative to storage capacity after all other water balance components were added
        help=(volact(step,upstream)-daystorcap(step,res_index(upstream)))+qinflow(step,res_index(upstream))
        ! total volume after all other water balance components were added
		help1=help2+help
        help=max(0.,help)
!write(*,'(2I4,6F11.1)')step,id_subbas_extern(upstream),help,help1,help2,help3,volact(step,upstream),daystorcap(step,upstream)
        ! total volume still larger than storage capacity -> calculate overspill
        IF (help1 > daystorcap(step,res_index(upstream))) THEN
!write(*,'(2I4,6F15.2)')step,id_subbas_extern(upstream),help,help2
         call reservoir_routing(upstream,help,help2)
!         help=max(0.,help-overflow(step,upstream))
         if(help2 > daystorcap(step,res_index(upstream))) lakeret(step,res_index(upstream))=lakeret(step,res_index(upstream))+(max(0.,volact(step,upstream)-help2))
         if(help2 <= daystorcap(step,res_index(upstream))) lakeret(step,res_index(upstream))=lakeret(step,res_index(upstream))+(max(0.,volact(step,upstream)-daystorcap(step,res_index(upstream))))

!write(*,'(2I4,3F15.4)')step,id_subbas_extern(upstream),help,volact(step,upstream),overflow(step,upstream)/86400.
        ELSE
		 volume_last(upstream)=0.
		 outflow_last(upstream)=0.
		 volact(step,upstream)=help1
		 overflow(step,res_index(upstream))=0.
        END IF
	  ELSE IF (volact(step,upstream) > daystorcap(step,res_index(upstream))) THEN
	    volact(step,upstream)=daystorcap(step,res_index(upstream))
		overflow(step,res_index(upstream))=overflow(step,res_index(upstream))+volact(step,upstream)-daystorcap(step,res_index(upstream))
      END IF

!write(*,'(2I4,5F15.4)')step,id_subbas_extern(upstream),qinflow(step,upstream)/86400.,overflow(step,upstream)/86400.,volact(step,upstream),lakeret(step,upstream)/86400.

!write(*,'(2I4,3F15.4)')step,id_subbas_extern(upstream),volact(step,upstream)

!if (step==20)stop
!if (step==41 .and. id_subbas_extern(upstream)==29) stop

      res_qout(step,res_index(upstream))=qintake(step,res_index(upstream))+overflow(step,res_index(upstream))+qbottom(step,res_index(upstream))


      IF (nbrbat(res_index(upstream)) /= 0) THEN
        DO j=1,nbrbat(res_index(upstream))-1
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
		elevhelp=elevhelp+dayminlevel(step,res_index(upstream))
        areahelp=dama(upstream)*(volact(step,upstream)**damb(upstream))
      END IF
      damareaact(upstream)=areahelp
      damelevact(res_index(upstream))=elevhelp

!write(*,'(I4,4F15.4)')d,dayminlevel(step,upstream),elevhelp,forma_factor(upstream),areahelp

      res_qout(step,res_index(upstream))=res_qout(step,res_index(upstream))/(86400./nt)
      qinflow(step,res_index(upstream))=qinflow(step,res_index(upstream))/(86400./nt)
      qintake(step,res_index(upstream))=qintake(step,res_index(upstream))/(86400./nt)
      overflow(step,res_index(upstream))=overflow(step,res_index(upstream))/(86400./nt)
	  qbottom(step,res_index(upstream))=qbottom(step,res_index(upstream))/(86400./nt)
	  withdraw_out(step,res_index(upstream))=withdraw_out(step,res_index(upstream))/(86400./nt)


! Calculation of reservoir surface area and reservoir volume when inflow discharges, outflow discharges
! and reservoir levels are still provided in the res_"ID_SUBBAS_EXTERN"_daily.dat
    ELSE IF (reservoir_balance == 0) THEN

!      write(*,*)t,d,hour,upstream,qinflow(step,upstream)
      daymaxdamarea(step,res_index(upstream))=daymaxdamarea(step,res_index(upstream))*1.e4
      daystorcap(step,res_index(upstream))=daystorcap(step,res_index(upstream))*1.e6
      daydamalert(step,res_index(upstream))=daydamalert(step,res_index(upstream))*1.e6
      daydamdead(step,res_index(upstream))=daydamdead(step,res_index(upstream))*1.e6
      damvol0(res_index(upstream))=volact(step,upstream)

	  res_qout(step,res_index(upstream))=qintake(step,res_index(upstream))+overflow(step,res_index(upstream))+qbottom(step,res_index(upstream))
	  damelevact(res_index(upstream))=damelev1(step,res_index(upstream))

      IF (nbrbat(res_index(upstream)) /= 0) THEN
        DO j=1,nbrbat(res_index(upstream))-1
          IF (damelevact(res_index(upstream)) >= elev_bat(j,upstream).AND.  &
                damelevact(res_index(upstream)) <= elev_bat(j+1,upstream)) THEN
            volhelp=vol_bat(j,upstream)+(damelevact(res_index(upstream))-elev_bat  &
                (j,upstream))/(elev_bat(j+1,upstream)-elev_bat(j,upstream))*  &
                (vol_bat(j+1,upstream)-vol_bat(j,upstream))
            areahelp=area_bat(j,upstream)+(damelevact(res_index(upstream))-elev_bat  &
                (j,upstream))/(elev_bat(j+1,upstream)-elev_bat(j,upstream))*  &
                (area_bat(j+1,upstream)-area_bat(j,upstream))

         END IF
        END DO
      ELSE
        elevhelp=damelevact(res_index(upstream))-dayminlevel(step,res_index(upstream))
        volhelp=forma_factor(upstream)*(elevhelp**3)
        areahelp=dama(upstream)*(volhelp**damb(upstream))
      END IF
      damareaact(upstream)=areahelp
	  volact(step,upstream)=volhelp
    ENDIF
!    damex(step,upstream)=res_qout(step,upstream)+withdrawal(upstream) tp not used

! Call sediment balance sub-routine
    IF (dosediment) then
      CALL semres (status,upstream)

!write(*,'(2I4,3F15.4)')d,id_subbas_extern(upstream),decstorcap(step,upstream),dayminlevel(step,upstream),daystorcap(step,upstream)

! Calculation of storage capacity reduction
      if (decstorcap(step,res_index(upstream))/=0.) then

! CASE 1: stage-area-volume curve is not provided
       if (nbrbat(res_index(upstream)) == 0) then
	    if (decstorcap(step,res_index(upstream)) >= daystorcap(step,res_index(upstream)) .or. daystorcap(step,res_index(upstream)) == 0.) then
	     write(*,'(A, i0, A)') 'Reservoir at outlet of sub-basin ', id_subbas_extern(upstream), ' completely silted up!'
		 storcap(upstream)=0.
	     daystorcap(step,res_index(upstream))=0.
		 do j=1,nbrbat(res_index(upstream))
          area_bat(j,upstream)=0.
          vol_bat(j,upstream)=0.
		 enddo
        else
         volhelp=MAX(0.,daystorcap(step,res_index(upstream))  &
			-decstorcap(step,res_index(upstream)))

! Calculation of minimum reservoir elevation
         if (daystorcap(step,res_index(upstream)) == 0.) then
	      dummy4=0.
	     else
	      dummy4=decstorcap(step,res_index(upstream))/daystorcap(step,res_index(upstream))
	     endif

	     elevhelp=max(0.,(maxlevel(upstream)-dayminlevel(step,res_index(upstream)))*dummy4)
	     dayminlevel(step,res_index(upstream))=dayminlevel(step,res_index(upstream))+elevhelp
	     forma_factor(upstream)=volhelp/((maxlevel(upstream)-dayminlevel(step,res_index(upstream)))**3)

	     daystorcap(step,res_index(upstream))=volhelp

	     elevhelp=max(0.,maxlevel(upstream)-dayminlevel(step,res_index(upstream)))
         daymaxdamarea(step,res_index(upstream))=dama(upstream)*(volhelp**damb(upstream))

! Calculation of dead water volume after sediment deposition
	     if (damdead(upstream) /= 0. .and. daydamdead(step,res_index(upstream)) /= 0.) then
	      elevhelp=max(0.,elevdead(res_index(upstream))-dayminlevel(step,res_index(upstream)))
          daydamdead(step,res_index(upstream))=forma_factor(upstream)*(elevhelp**3.)
	     else
	      daydamdead(step,res_index(upstream))=0.
	     endif

! Calculation of alert water volume after sediment deposition
 	     if (damalert(upstream) /= 0. .and. daydamalert(step,res_index(upstream)) /= 0.) then
	      elevhelp=max(0.,elevalert(res_index(upstream))-dayminlevel(step,res_index(upstream)))
          daydamalert(step,res_index(upstream))=forma_factor(upstream)*(elevhelp**3.)
	     else
	      daydamalert(step,res_index(upstream))=0.
	     endif
	    endif

! CASE 2: stage-area-volume curve is provided
	   else
! CASE 2a: no detailed information of the reservoir geometry is provided (no cross section)
	    if (nbrsec(res_index(upstream)) == 0) then
	     if (decstorcap(step,res_index(upstream)) >= daystorcap(step,res_index(upstream)) .or. daystorcap(step,res_index(upstream)) == 0.) then
	      write(*,*) 'the resevoir located at the outlet point of sub-basin:',id_subbas_extern(upstream)
	      write(*,*) 'lost its total storage capacity due to sediment deposition'
		  storcap(upstream)=0.
	      daystorcap(step,res_index(upstream))=0.
		  do j=1,nbrbat(res_index(upstream))
           area_bat(j,upstream)=0.
           vol_bat(j,upstream)=0.
		  enddo
         else

	     daystorcap(step,res_index(upstream))=MAX(0.,daystorcap(step,res_index(upstream))  &
			-decstorcap(step,res_index(upstream)))

! Change on stage-area-volume curve due to erosion and deposition processes
! Calculation of minimum reservoir elevation
	     elevhelp=dayminlevel(step,res_index(upstream))
		 p=0
         do j=1,nbrbat(res_index(upstream))
		  if (elev_bat(j,upstream) >= maxlevel(upstream)) then
		   dummy1=j
		    exit
		  endif
		 enddo
         do j=1,nbrbat(res_index(upstream))
	      if (elev_bat(j,upstream) < maxlevel(upstream)) then
		   dummy5=real(j)/real(dummy1)
		  else if (elev_bat(j,upstream) >= maxlevel(upstream)) then
		   dummy5=1.
		  endif
!write(*,'(3I4,F15.5,3F15.2)')step,j,dummy1,dummy5,elev_bat(j,upstream),dayminlevel(step,upstream),maxlevel(upstream)
!write(*,'(3I4,F15.5,4F15.2)')step,j,dummy1,dummy5,vol_bat(j,upstream),decstorcap(step,upstream)*dummy5,vol_bat(j,upstream)-(decstorcap(step,upstream)*dummy5)
		  dummy5=max(0.,vol_bat(j,upstream)-(decstorcap(step,res_index(upstream))*dummy5))
		  if (dummy5 == 0.)p=j
		  if (j == p+1) then
           DO q=1,nbrbat(res_index(upstream))-1
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
		 dayminlevel(step,res_index(upstream))=elevhelp

         do j=1,nbrbat(res_index(upstream))
		  if (elev_bat(j,upstream) <= elevhelp) then
		   vol_bat(j,upstream)=0.
		   area_bat(j,upstream)=0.
		  endif
	     enddo

! Calculation of dead water volume after sediment deposition
         do j=1,nbrbat(res_index(upstream))-1
	      if (damdead(upstream) /= 0. .and. daydamdead(step,res_index(upstream)) /= 0.) then
           if (elevdead(res_index(upstream)) >= elev_bat(j,upstream).AND.  &
               elevdead(res_index(upstream)) <= elev_bat(j+1,upstream)) THEN
            daydamdead(step,res_index(upstream))=vol_bat(j,upstream)+(elevdead(res_index(upstream))-elev_bat  &
                (j,upstream))/(elev_bat(j+1,upstream)-elev_bat(j,upstream))*  &
                (vol_bat(j+1,upstream)-vol_bat(j,upstream))
           endif
	      else
	       daydamdead(step,res_index(upstream))=0.
	      endif
! Calculation of alert water volume after sediment deposition
 	      if (damalert(upstream) /= 0. .and. daydamalert(step,res_index(upstream)) /= 0.) then
           if (elevalert(res_index(upstream)) >= elev_bat(j,upstream).AND.  &
               elevalert(res_index(upstream)) <= elev_bat(j+1,upstream)) THEN
            daydamalert(step,res_index(upstream))=vol_bat(j,upstream)+(elevalert(res_index(upstream))-elev_bat  &
                (j,upstream))/(elev_bat(j+1,upstream)-elev_bat(j,upstream))*  &
                (vol_bat(j+1,upstream)-vol_bat(j,upstream))
           endif
	      else
	       daydamalert(step,res_index(upstream))=0.
	      endif
! Calculation of maximum reservoir area after sediment deposition
 	      if (daymaxdamarea(step,res_index(upstream)) /= 0.) then
           if (maxlevel(upstream) >= elev_bat(j,upstream).AND.  &
               maxlevel(upstream) <= elev_bat(j+1,upstream)) THEN
            daymaxdamarea(step,res_index(upstream))=area_bat(j,upstream)+(maxlevel(upstream)-elev_bat  &
                (j,upstream))/(elev_bat(j+1,upstream)-elev_bat(j,upstream))*  &
                (area_bat(j+1,upstream)-area_bat(j,upstream))
           endif
	      else
	       daymaxdamarea(step,res_index(upstream))=0.
	      endif
	     enddo
	    endif

! CASE 2b: detailed information of the reservoir geometry is provided (cross sections)
       else
	    if (decstorcap(step,res_index(upstream)) >= daystorcap(step,res_index(upstream)) .or. daystorcap(step,res_index(upstream)) == 0.) then
	     write(*,*) 'the resevoir located at the outlet point of sub-basin:',id_subbas_extern(upstream)
	     write(*,*) 'lost its total storage capacity due to sediment deposition'
		 storcap(upstream)=0.
	     daystorcap(step,res_index(upstream))=0.
		 do j=1,nbrbat(res_index(upstream))
          area_bat(j,upstream)=0.
          vol_bat(j,upstream)=0.
		 enddo
        else
         do j=1,nbrbat(res_index(upstream))-1
	      if (storcap(upstream) /= 0. .and. daystorcap(step,res_index(upstream)) /= 0.) then
           if (maxlevel(upstream) >= elev_bat(j,upstream).AND.  &
               maxlevel(upstream) <= elev_bat(j+1,upstream)) THEN
            daystorcap(step,res_index(upstream))=vol_bat(j,upstream)+(maxlevel(upstream)-elev_bat  &
                (j,upstream))/(elev_bat(j+1,upstream)-elev_bat(j,upstream))*  &
                (vol_bat(j+1,upstream)-vol_bat(j,upstream))
           endif
	      else
	       daystorcap(step,res_index(upstream))=0.
	      endif
! Calculation of dead water volume after sediment deposition
	      if (damdead(upstream) /= 0. .and. daydamdead(step,res_index(upstream)) /= 0.) then
           if (elevdead(res_index(upstream)) >= elev_bat(j,upstream).AND.  &
               elevdead(res_index(upstream)) <= elev_bat(j+1,upstream)) THEN
            daydamdead(step,res_index(upstream))=vol_bat(j,upstream)+(elevdead(res_index(upstream))-elev_bat  &
                (j,upstream))/(elev_bat(j+1,upstream)-elev_bat(j,upstream))*  &
                (vol_bat(j+1,upstream)-vol_bat(j,upstream))
           endif
	      else
	       daydamdead(step,res_index(upstream))=0.
	      endif
! Calculation of alert water volume after sediment deposition
 	      if (damalert(upstream) /= 0. .and. daydamalert(step,res_index(upstream)) /= 0.) then
           if (elevalert(res_index(upstream)) >= elev_bat(j,upstream).AND.  &
               elevalert(res_index(upstream)) <= elev_bat(j+1,upstream)) THEN
            daydamalert(step,res_index(upstream))=vol_bat(j,upstream)+(elevalert(res_index(upstream))-elev_bat  &
                (j,upstream))/(elev_bat(j+1,upstream)-elev_bat(j,upstream))*  &
                (vol_bat(j+1,upstream)-vol_bat(j,upstream))
           endif
	      else
	       daydamalert(step,res_index(upstream))=0.
	      endif
! Calculation of maximum reservoir area after sediment deposition
 	      if (daymaxdamarea(step,res_index(upstream)) /= 0.) then
           if (maxlevel(upstream) >= elev_bat(j,upstream).AND.  &
               maxlevel(upstream) <= elev_bat(j+1,upstream)) THEN
            daymaxdamarea(step,res_index(upstream))=area_bat(j,upstream)+(maxlevel(upstream)-elev_bat  &
                (j,upstream))/(elev_bat(j+1,upstream)-elev_bat(j,upstream))*  &
                (area_bat(j+1,upstream)-area_bat(j,upstream))
           endif
	      else
	       daymaxdamarea(step,res_index(upstream))=0.
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
	 WRITE(11,'(4I6,2f10.3,2f13.1,6f10.3,3f14.1)')id_subbas_extern(upstream),t,d,hour,qlateral(step,res_index(upstream)),qinflow(step,res_index(upstream)), & 
                etdam(step,res_index(upstream)),precdam(step,res_index(upstream)),  &
				qintake(step,res_index(upstream)),overflow(step,res_index(upstream)),qbottom(step,res_index(upstream)),res_qout(step,res_index(upstream)), &
				withdraw_out(step,res_index(upstream)),damelevact(res_index(upstream)),damareaact(upstream),volact(step,upstream)
     CLOSE(11)
	 ENDIF


! print storage losses due to reservoir sedimentation
	 IF (f_res_vollost) THEN
     OPEN(11,FILE=pfadn(1:pfadi)//'res_'//trim(adjustl(subarea))//'_vollost.out',STATUS='old',  &
		  POSITION='append')
	 WRITE(11,'(4I6,3f15.2)')id_subbas_extern(upstream),t,d,hour,daydamdead(step,res_index(upstream)), &
			daydamalert(step,res_index(upstream)),daystorcap(step,res_index(upstream))
     CLOSE(11)
	 ENDIF

! print temporal evolution of the stage-area-volume curve, when the initial curve is given by the user
     write(fmtstr,'(a,i0,a)')'(4I6,',nbrbat(res_index(upstream)),'F15.2)'		!generate format string
     IF (nbrbat(res_index(upstream)) /= 0) THEN
	  IF (f_res_cav) THEN
      OPEN(11,FILE=pfadn(1:pfadi)//'res_'//trim(adjustl(subarea))//'_cav.out',STATUS='old',  &
		    POSITION='append')
		WRITE(11,fmtstr)id_subbas_extern(upstream),t,d,hour,(elev_bat(j,upstream),j=1,nbrbat(res_index(upstream)))
		WRITE(11,fmtstr)id_subbas_extern(upstream),t,d,hour,(area_bat(j,upstream),j=1,nbrbat(res_index(upstream)))
		WRITE(11,fmtstr)id_subbas_extern(upstream),t,d,hour,(vol_bat(j,upstream),j=1,nbrbat(res_index(upstream)))
      CLOSE(11)
	  ENDIF
	 ENDIF
    ELSE
	 daydamelevact(step,res_index(upstream))=damelevact(res_index(upstream)) !A
	 daydamareaact(step,res_index(upstream))=damareaact(upstream) !A
	 DO j=1,nbrbat(res_index(upstream))
	   dayelev_bat(step,j,res_index(upstream))=elev_bat(j,upstream) !A
	   dayarea_bat(step,j,res_index(upstream))=area_bat(j,upstream) !A
	   dayvol_bat(step,j,res_index(upstream))=vol_bat(j,upstream) !Anne: changed "upstream" to "res_index(upstream)" to save memory
	 ENDDO
	ENDIF


    volact(step,upstream)=volact(step,upstream)/1.e6 !convert to 10^6 m3
	daymaxdamarea(step,res_index(upstream))=daymaxdamarea(step,res_index(upstream))/1.e4
	daystorcap(step,res_index(upstream))=daystorcap(step,res_index(upstream))/1.e6
	daydamalert(step,res_index(upstream))=daydamalert(step,res_index(upstream))/1.e6
	daydamdead(step,res_index(upstream))=daydamdead(step,res_index(upstream))/1.e6
	damflow(upstream)=damflow(upstream)/(86400./nt)
	qoutlet(upstream)=qoutlet(upstream)/(86400./nt)
	withdrawal(upstream)=withdrawal(upstream)/(86400./nt)


	IF (step < dayyear*nt) THEN !Till: use current values to initialize the next timestep
	  volact(step+1,upstream)=volact(step,upstream)
	  daystorcap(step+1,res_index(upstream))=daystorcap(step,res_index(upstream))
	  daydamalert(step+1,res_index(upstream))=daydamalert(step,res_index(upstream))
	  daydamdead(step+1,res_index(upstream))=daydamdead(step,res_index(upstream))
	  daymaxdamarea(step+1,res_index(upstream))=daymaxdamarea(step,res_index(upstream))
	  dayminlevel(step+1,res_index(upstream))=dayminlevel(step,res_index(upstream))
	END IF

!write(*,'(2I4,3F15.4)')d,id_subbas_extern(upstream),dayminlevel(step,upstream),decstorcap(step,upstream)
!if (d==4)stop

ELSE  ! reservoir does not (yet) exist

    if (res_index(upstream)/=0) then
        if (river_transport.eq.1)then
          res_qout(step,res_index(upstream))=qout(step,upstream)
        else
        res_qout(step,res_index(upstream))=r_qout(2,upstream)+qlateral(step,res_index(upstream))
        endif
    endif
        
  ENDIF

! END of STATUS = 2 !regular call during timestep
END IF

! -----------------------------------------------------------------------
IF (STATUS == 3) THEN

! close intake.dat
if( any(f_intake_obs) .and. (t == tstop) ) then
    close(101)
    deallocate(r_qintake)
endif

! Output files of Reservoir Modules
  IF (reservoir_print == 1) THEN
    DO i=1,subasin
      IF (res_flag(i) .and. t >= damyear(i)) THEN
        WRITE(subarea,*)id_subbas_extern(i)
		IF (f_res_watbal) THEN
	    OPEN(11,FILE=pfadn(1:pfadi)//'res_'//trim(adjustl(subarea))//'_watbal.out',STATUS='old',  &
			POSITION='append')
	    DO d=1,dayyear
	      DO ih=1,nt
		    hour=ih
            step=(d-1)*nt+hour
	        WRITE(11,'(4(I6,a),2(f10.3,a),2(f13.1,a),6(f10.3,a),3(f14.1,a))')id_subbas_extern(i),char(9),&
            t,char(9),d,char(9),hour,char(9),qlateral(step,res_index(i)),char(9),qinflow(step,res_index(i)),char(9),etdam(step,res_index(i)),char(9),&
            precdam(step,res_index(i)),char(9),qintake(step,res_index(i)),char(9),overflow(step,res_index(i)),char(9),qbottom(step,res_index(i)),char(9),&
    !Anne changed, eg., daydamareaact(step,i) to daydamareaact(step,res_index(i))
            res_qout(step,res_index(i)),char(9),withdraw_out(step,res_index(i)),char(9),daydamelevact(step,res_index(i)),char(9),daydamareaact(step,res_index(i)),char(9),volact(step,i)*1.e6,char(9)
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
			WRITE(11,'(4I6,3f15.2)')id_subbas_extern(i),t,d,hour,daydamdead(step,res_index(i))*1.e6, &
				daydamalert(step,res_index(i))*1.e6,daystorcap(step,res_index(i))*1.e6
		  ENDDO
		ENDDO
        CLOSE(11)
		ENDIF
      ENDIF
      if (res_index(i) /= 0) then !Anne inserted this line
          IF (res_flag(i) .and. t >= damyear(i) .and. nbrbat(res_index(i)) /= 0) THEN
            WRITE(subarea,*)id_subbas_extern(i)
		    IF (f_res_cav) THEN
                OPEN(11,FILE=pfadn(1:pfadi)//'res_'//trim(adjustl(subarea))//'_cav.out',STATUS='old',  &
		            POSITION='append')
	            write(fmtstr,'(a,i0,a)')'(4I6,',nbrbat(res_index(i)),'F15.2)'		!generate format string
		        DO d=1,dayyear
	              DO ih=1,nt
		            hour=ih
                    step=(d-1)*nt+hour
			        WRITE(11,fmtstr)id_subbas_extern(i),t,d,hour,(dayelev_bat(step,j,res_index(i)),j=1,nbrbat(res_index(i)))
			        WRITE(11,fmtstr)id_subbas_extern(i),t,d,hour,(dayarea_bat(step,j,res_index(i)),j=1,nbrbat(res_index(i)))
			        WRITE(11,fmtstr)id_subbas_extern(i),t,d,hour,(dayvol_bat(step,j,res_index(i)),j=1,nbrbat(res_index(i))) !Anne: changed i to res_index(i) to save memory
		          ENDDO
		        ENDDO
                CLOSE(11)
		    ENDIF
          ENDIF
       endif	 !Anne   
    ENDDO
  ENDIF


  IF (dosediment) CALL semres (status,upstream)

END IF



RETURN
END SUBROUTINE reservoir
