SUBROUTINE hymo_all(STATUS)
!Till: computationally relevant in hourly version: fixed bug that led to excess runoff;
! fixed faulty summation of hortonian runoff in hourly version
! minor optimizations and code beautifying
!2009-04-01

!Till: computationally irrelevant, reformatted *water_subbasin.out (and others) with tab-separated output
!2009-03-31

!Till: unified output files, used tab as separators
!2009-03-17

!Till: fixed faulty sediment distribution to lower TCs (was: all sediment to next lower TC and river); relevant for sediment modelling
!minor optimizations, variable name clarifications
!2008-10-15

!Till: error in seasonality calculation, when simulation day was before or after last seasonal node - fixed
!2008-10-09

!Till: fixed formatting error in tc_*.dat when svc_comb was misspecified; output optimizations
!2008-10-08

!Till: fixed bug with julian day: wrong for years with endmonth/=12. Resulted in wrong vegetation dynamics
!2008-09-24

!Till: bugfix in pre-specified sediment outflow of selected subbasins
!2008-08-28

!Till: conditional saving of system state
!2008-07-29

!Till: improved TC-wise theta output
!2008-07-16

!Till: optimized and corrected TC-wise sediment output for non_January starting month 
!2008-07-15

!Till: implemented optional pre-specified sediment outflow of selected subbasins
!2008-07-11

!Till: implemented optional pre-specified outflow of selected subbasins
!2008-07-03

!Till: corrected output of gw_loss (wrong calculation)
! corrected computation of balance (for LU), computationally irrelevant
!2008-06-24

!Pedro
!TC-wise output of deposition information
!2008-02-07

!Till: optional computation of erosion at subbasin-scale
! tab separated output for daily_sediment_production.dat
!2007-11-13

!Till: corrected routing of hourly sediment flux
!2007-11-7

!Till: time-variate kfcorr
!2007-10-29

!Till: fixed bug in "redistribution of gw-flow to lowermost tc" that caused hangup
!2007-09-25

!Till: gw discharge is routed to subsurface flow of lowermost TC with specified fraction frac_direct_gw
!2007-08-21

!Till: limit groundwater discharge to maximum available stored volume - prevent negative values in gw-storage
!2007-08-13

! Till & George: fixed bug in plant parameter interpolation when starting with month/=1 
! 2007-06-21 / 2007-07-17

! Till: improved output format with tabs (for deep_gw_recharge only) 
! 2007-01-19

! Till: fixed bug in hourly version: hourly water flux did not contain direct runoff and ground water discharge
! prepared threshold value for minimum water yield necessary to produce riverflow - currently deactivated 
! 2007-01-10

! Till: fixed bug in sediment routing between TCs when doalllattc=.FALSE. (mass deficit)
! 2006-08-21

! Till: daily_water_subbasin2.out is created and contains the same as daily_water_subbasin.out, but in m**3 (easier to postprocess when low flows prevail)
! 2005-12-20

! Eva: corrected conditional calls to lake (acud = TRUE) 
! 2006-01-xx

! Till: changed addressing scheme of sediment_subbasin_t from (366*nt,subasin,n_sed_class) to (366,nt,subasin,n_sed_class)
! 2005-12-20

! Till: outcommented unnecessary vars
! added initialisation of latred
! added missing initialisation of horithact and intercept - this strongly affects hourly results!
! removed unnnecessary differences between hourly and daily version
! subdaily output files are only created in subdaily resolution
! 2005-11-03

! Till: sediment_subbasin_t contains subdaily sediment freight
! initialisation of daily values is done at the beginning of the year for entire domain instead of doing it every day and subbasin
! 2005-10-25

! Till: write selected output files only
! water_subbasin_t contains subdaily runoff
! 2005-10-24

! Till: minor optimisations and removal of unnecessary lines
! 2005-10-20

! Till: water_subbasin converted from m**3 into m**3/s 
! 2005-10-13
  
! Till: routing is no longer called from hymo, but wasa.f90
! 2005-10-12
 
! Till: fixed formatting bug in daily_sediment_production.out
! 2005-09-29

! Till: variable names clarified, comments
! 2005-09-26

! Till: variable names clarified
! 2005-08-09

! Code converted using TO_F90 by Alan Miller
! Date: 2005-06-30  Time: 13:47:07

use lake_h
use climo_h
use common_h
use hymo_h
use params_h
use routing_h
use time_h
use reservoir_h
use erosion_h
use model_state_io

IMPLICIT NONE


INTEGER, INTENT(IN)                  :: STATUS
! status of call (0=initialization run, 1=initialization year,
!                 2=calculation day,    3=finalization year)


! counters
INTEGER :: i,j,sb_counter,lu_counter,tc_counter,sc,h,hh

!(internal) id of TC-instance (unique subbas-LU-TC-combination)
INTEGER :: tcid_instance

! ids of components in work
INTEGER :: i_subbas,i_lu,id_tc_type,isoil

! julian day
INTEGER :: julian_day
! actual soil water content of terrain component
REAL :: thact, thactroot, temp2, temp3, diff, prec, temparea
REAL :: precday, prechall(24)

! area of terrain component (km²)
REAL :: tcarea
! areal fractions of terrain components in current LU
REAL :: fractemp(maxterrain)
! water flow between terrain components
REAL :: surfflow_in(maxterrain),surfflow_out(maxterrain)	!Till: overland inflow and outflow for all TCs of current LU
REAL :: sublat_in(maxterrain),sublat_out(maxterrain)

! sediment flow between terrain components
REAL,allocatable,save :: sed_in(:,:),sed_out(:,:)
!REAL,allocatable,save :: sed_in(:,:)
!REAL,pointer,save :: sed_out(:,:)

! water balance of landscape unit
REAL :: balance,balance_tc

! dummy
INTEGER :: dummy,count
REAL ::	rtemp,rtemp2

! vegetation characteristics of current day
REAL :: rootd_act(nveg), height_act(nveg)
REAL :: lai_act(nveg), alb_act(nveg)

REAL :: frac_satsu

!variables of daily version only
REAL :: gwr,deepgwr

! hourly variables
REAL :: aeth,laih,soileth,inth,horth
INTEGER :: timestep_counter

character(len=1000) :: fmtstr	!string for formatting file output

!conrad: theta output
INTEGER :: tc_counter_all						!counter for theta_tc Till: counts instances of TCs
INTEGER :: dig_sub,dig_lu,dig_tc				!digits necessary to print external IDs
!INTEGER*8, allocatable :: tc_idx(:)	!TC-IDs used for TC-wise output
CHARACTER(LEN=30),allocatable :: tc_idx(:)	!TC-IDs used for TC-wise output

!!Print hydrologic variable on TC scale. If not used, DISABLE
!!***********************************************************
!CHARACTER(20) :: year_name,day_name
!INTEGER :: k,c,dummy1,dummy2,year_print(6),day_print(6)
!!CHARACTER(LEN=12),allocatable :: print_timelabel(:)	!labels to be used for printing TC output
!
! year_print(1)=2005;day_print(1)=258!;timestep_print(1)=1	!Till: specify dates for which TC-wise output is desired
! year_print(2)=2005;day_print(2)=266!;timestep_print(2)=1
! year_print(3)=2005;day_print(3)=267!;timestep_print(3)=1
! year_print(4)=2005;day_print(4)=268!;timestep_print(4)=1
! year_print(5)=2005;day_print(5)=269!;timestep_print(5)=1
!
!
!
!!***********************************************************

!CCCCCCCCCCCCCCCCCCCC MODULE CODE CCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
IF (STATUS == 0) THEN	!Till: initialisation before first run
!** Initialization from datafile
	
  CALL readhymo

  if (do_pre_outflow) then		!if water outflow from upstream subbasins is given
		call read_pre_subbas_outflow		!read
  end if
  
!** Additional intializations
  DO i_subbas=1,subasin
!George    lakeoutflow(dprev,i_subbas)=0.
    water_subbasin(dprev,i_subbas)=0.
    soilwater(dprev,:)=0.
  END DO
  DO sb_counter=1,subasin
    deepgw(sb_counter,:)=0.
  END DO
  
  allocate(sed_in(maxval(nbrterrain),n_sed_class),sed_out(maxval(nbrterrain),n_sed_class))


  latred=0.				!set lateral subsurface runof between SVCs to zero
  horithact(:,:,:)=0.	!set soil moisture of each horizon to 0
  intercept(:,:)=0.		!set interception storage to 0
  
  if (doloadstate) then
	call init_model_state		!load initital conditions from file
  else
	!** Initialise soil moisture of each horizon (horithact)
	!   and soil moisture of each terrain component (soilwater) with fixed values ii: put this into init_soil_conds
	DO sb_counter=1,subasin
		DO lu_counter=1,nbr_lu(sb_counter)
		  i_lu=id_lu_intern(lu_counter,sb_counter)
		  IF (gw_flag(i_lu) == 0 .OR. gw_flag(i_lu) == 1) THEN
			DO tc_counter=1,nbrterrain(i_lu)
			  tcid_instance=tcallid(sb_counter,lu_counter,tc_counter)
			  DO sc=1,nbr_svc(tcid_instance)
				isoil=id_soil_intern(sc,tcid_instance)
				DO h=1,nbrhori(isoil)
				  horithact(tcid_instance,sc,h)= soilpwp(isoil,h)*horiz_thickness(tcid_instance,sc,h) !Till: set water content to wilting point
				  
				  
				END DO
    
				soilwater(dprev,tcid_instance)= soilwater(dprev,tcid_instance)+  &
					sum(horithact(tcid_instance,sc,:))* frac_svc(sc,tcid_instance)
			  END DO
			END DO
		  ELSE IF (gw_flag(i_lu) == 99) THEN	!Till: gw_flag==99: experimental flag for representing groundwater close to the surface, not documented
			DO tc_counter=1,nbrterrain(i_lu)
			  tcid_instance=tcallid(sb_counter,lu_counter,tc_counter)
			  DO sc=1,nbr_svc(tcid_instance)
				isoil=id_soil_intern(sc,tcid_instance)
				DO h=1,nbrhori(isoil)
				  temp3=sum(horiz_thickness(tcid_instance,sc,1:h))	!Till: thickness of entire soil profile
				  temp2=temp3-horiz_thickness(tcid_instance,sc,h)	!Till: thickness of all horizons except the deepest
				  IF (temp3 <= gw_dist(i_lu)) THEN					!Till: if gw table is below deepest horizon...
					horithact(tcid_instance,sc,h)= soilpwp(isoil,h)*horiz_thickness(tcid_instance,sc,h)	!Till: set initial soil moisture to pwp
				!IF (.NOT. dohour) THEN
				!      thetas(isoil,h)*horiz_thickness(tcid_instance,sc,h)
				!END IF
				  ELSE IF (temp2 >= gw_dist(i_lu)) THEN				!Till: if gw is above deepest horizon...
					horithact(tcid_instance,sc,h)= thetas(isoil,h)*horiz_thickness(tcid_instance,sc,h)	!Till: lowest horizon is saturated completely ii warum sind die darüberliegenden nicht ggf. auch gesättigt?
				  ELSE IF (temp2 < gw_dist(i_lu) .AND.  &
						temp3 > gw_dist(i_lu)) THEN					!Till: gw is within lowest horizon...
					horithact(tcid_instance,sc,h)=  &
						thetas(isoil,h)*(temp3-gw_dist(i_lu))+  &	!Till: saturated below gw level
						soilpwp(isoil,h)*(gw_dist(i_lu)-temp2)		!Till: pwp above gw level
				  END IF
				END DO
				soilwater(dprev,tcid_instance)= soilwater(dprev,tcid_instance)+  &
					sum(horithact(tcid_instance,sc,:))* frac_svc(sc,tcid_instance)
			  END DO
			END DO
		  END IF
		END DO
	END DO

  
  end if

  
  
!   initialize saturated fraction of TC
!  frac_sat(:,:)=0.5		!Till: added for quick response in runoff while using short time series
  
if (doacud)  CALL lake(0,dummy)

  
! create and open output files
! Output daily water contribution to river (m**3/s)
  OPEN(11,FILE=pfadn(1:pfadi)// 'daily_water_subbasin.out', STATUS='replace')
  IF (f_daily_water_subbasin) THEN	
	WRITE(11,'(a)') 'daily river flow [m3/s] for all sub-basins (MAP-IDs)'
	write(fmtstr,'(a,i0,a)')'(a,',subasin,'(a,i0))'		!generate format string
	WRITE(11,fmtstr)'Year'//char(9)//'Julian_Day',(char(9),id_subbas_extern(i),i=1,subasin)
	CLOSE(11)
  ELSE				!delete any existing file, if no output is desired
	CLOSE(11,status='delete')
  END IF

! Output daily water contribution to river (m**3)
  OPEN(11,FILE=pfadn(1:pfadi)// 'daily_water_subbasin2.out', STATUS='replace')
  IF (f_daily_water_subbasin) THEN	
	WRITE(11,'(a)') 'daily river flow [m3] for all sub-basins (MAP-IDs)'
	write(fmtstr,'(a,i0,a)')'(a,',subasin,'(a,i0))'		!generate format string
	WRITE(11,fmtstr)'Year'//char(9)//'Julian_Day',(char(9),id_subbas_extern(i),i=1,subasin)
	CLOSE(11)
  ELSE				!delete any existing file, if no output is desired
	CLOSE(11,status='delete')
  END IF

! Output daily actual evapotranspiration
  OPEN(11,FILE=pfadn(1:pfadi)//  &
      'daily_actetranspiration.out', STATUS='replace')
  IF (f_daily_actetranspiration) THEN	
    WRITE(11,'(a)') 'daily actual evapotranspiration [mm/d]  &
      for all sub-basins (MAP-IDs)'
	write(fmtstr,'(a,i0,a)')'(a,',subasin,'(a,i0))'		!generate format string
	WRITE(11,fmtstr)'Year'//char(9)//'Julian_Day',(char(9),id_subbas_extern(i),i=1,subasin)
	CLOSE(11)
  ELSE				!delete any existing file, if no output is desired
	CLOSE(11,status='delete')
  END IF

! Output daily potential evapotranspiration
  OPEN(11,FILE=pfadn(1:pfadi)//'daily_potetranspiration.out',  &
      STATUS='replace')
  IF (f_daily_potetranspiration) THEN	
	WRITE(11,'(a)') 'daily potential evapotranspiration [mm/d]  &
      for all sub-basins (MAP-IDs)'
	write(fmtstr,'(a,i0,a)')'(a,',subasin,'(a,i0))'		!generate format string
	WRITE(11,fmtstr)'Year'//char(9)//'Julian_Day',(char(9),id_subbas_extern(i),i=1,subasin)
	CLOSE(11)
  ELSE				!delete any existing file, if no output is desired
	CLOSE(11,status='delete')
  END IF

! Output daily soil moisture
  OPEN(11,FILE=pfadn(1:pfadi)//'daily_theta.out', STATUS='replace')
  IF (f_daily_theta) THEN	
	WRITE(11,'(a)') 'soil moisture in profile [mm] for all sub-basins (MAP-IDs)'
	write(fmtstr,'(a,i0,a)')'(a,',subasin,'(a,i0))'		!generate format string
	WRITE(11,fmtstr)'Year'//char(9)//'Julian_Day',(char(9),id_subbas_extern(i),i=1,subasin)
	CLOSE(11)
  ELSE				!delete any existing file, if no output is desired
	CLOSE(11,status='delete')
  END IF

! Output daily Hortonian Overland Flow
  OPEN(11,FILE=pfadn(1:pfadi)//'daily_qhorton.out', STATUS='replace')
  IF (f_daily_qhorton) THEN	
	WRITE(11,'(a)') 'horton overland flow [m**3] for all sub-basins (MAP-IDs)'
	write(fmtstr,'(a,i0,a)')'(a,',subasin,'(a,i0))'		!generate format string
	WRITE(11,fmtstr)'Year'//char(9)//'Julian_Day',(char(9),id_subbas_extern(i),i=1,subasin)
	CLOSE(11)
  ELSE				!delete any existing file, if no output is desired
	CLOSE(11,status='delete')
  END IF

! Output daily total overland flow
  OPEN(11,FILE=pfadn(1:pfadi)//'daily_total_overlandflow.out',  &
      STATUS='replace')
  IF (f_daily_total_overlandflow) THEN	
	WRITE(11,'(a)') 'total overland flow [m**3] for all sub-basins (MAP-IDs)'
	write(fmtstr,'(a,i0,a)')'(a,',subasin,'(a,i0))'		!generate format string
	WRITE(11,fmtstr)'Year'//char(9)//'Julian_Day',(char(9),id_subbas_extern(i),i=1,subasin)
	CLOSE(11)
  ELSE				!delete any existing file, if no output is desired
	CLOSE(11,status='delete')
  END IF


! Output deep ground water recharge (losses from model domain / into lin. storage)
  OPEN(11,FILE=pfadn(1:pfadi)//'deep_gw_recharge.out',  &
      STATUS='replace')
  IF (f_deep_gw_recharge) THEN	
	WRITE(11,'(a)') 'total deep ground water recharge [m3] for all sub-basins (MAP-IDs)'
	write(fmtstr,'(a,i0,a)')'(a4,a,a4,',subasin,'(a,i0))'		!generate format string
	WRITE(11,fmtstr)'Year', char(9),' Day', (char(9),id_subbas_extern(i),i=1,subasin)		!tab separated output
	CLOSE(11)
  ELSE				!delete any existing file, if no output is desired
	CLOSE(11,status='delete')
  END IF

  ! Output deep ground water discharge (component of water_subbasin)
  OPEN(11,FILE=pfadn(1:pfadi)//'deep_gw_discharge.out',  &
      STATUS='replace')
  IF (f_deep_gw_discharge) THEN	
	allocate (deep_gw_discharge(366,subasin))
	WRITE(11,'(a)') 'total deep ground water discharge [m3] for all sub-basins (MAP-IDs)'
	write(fmtstr,'(a,i0,a)')'(a4,a,a4,',subasin,'(a,i0))'		!generate format string
	WRITE(11,fmtstr)'Year', char(9),' Day', (char(9),id_subbas_extern(i),i=1,subasin)		!tab separated output
	CLOSE(11)
  ELSE				!delete any existing file, if no output is desired
	CLOSE(11,status='delete')
  END IF

! Output ground water loss (deep percolation in LUs with no ground water flag)
  OPEN(11,FILE=pfadn(1:pfadi)//'gw_loss.out',  &
      STATUS='replace')
  IF (f_gw_loss) THEN	
	allocate (gw_loss(366,subasin))
	WRITE(11,'(a)') 'ground water loss from model domain [m3] for all sub-basins (MAP-IDs)'
	write(fmtstr,'(a,i0,a)')'(a4,a,a4,',subasin,'(a,i0))'		!generate format string
	WRITE(11,fmtstr)'Year', char(9),' Day', (char(9),id_subbas_extern(i),i=1,subasin)		!tab separated output
	CLOSE(11)
  ELSE				!delete any existing file, if no output is desired
	CLOSE(11,status='delete')
  END IF


!     Output total subsurface runoff (m³/d)
  OPEN(11,FILE=pfadn(1:pfadi)//'daily_subsurface_runoff.out',  &
      STATUS='replace')
  IF (f_daily_subsurface_runoff) THEN	
	WRITE(11,'(a)') 'total subsurface runoff [m**3/d]  &
      for all sub-basins (MAP-IDs)'
	write(fmtstr,'(a,i0,a)')'(a,',subasin,'(a,i0))'		!generate format string
	WRITE(11,fmtstr)'Year'//char(9)//'Julian_Day',(char(9),id_subbas_extern(i),i=1,subasin)
	CLOSE(11)
  ELSE				!delete any existing file, if no output is desired
	CLOSE(11,status='delete')
  END IF

  ! Output sub-daily water flux
  OPEN(11,FILE=pfadn(1:pfadi)//'water_subbasin.out', STATUS='replace')
  IF (f_water_subbasin) THEN	!Till: if sub-daily resolution is required
	WRITE(11,'(a)') 'river flow [m3/s] for all sub-basins (MAP-IDs)'
	
	write(fmtstr,'(a,i0,a)')'(a,',subasin,'(a,i0))'		!generate format string
	WRITE(11,fmtstr)'Year'//char(9)//'Day'//char(9)//'Timestep', (char(9),id_subbas_extern(i),i=1,subasin)
	CLOSE(11)
  ELSE				!delete any existing file, if no output is desired
	CLOSE(11,status='delete')
  END IF
  
  !     Output sediment production (t)
  OPEN(11,FILE=pfadn(1:pfadi)//'daily_sediment_production.out', STATUS='replace')
  IF (f_daily_sediment_production .AND. dosediment) THEN
	WRITE(11,'(a,i2,a1)') 'total sediment production [t] for all sub-basins (MAP-IDs) and particle classes{',n_sed_class,'}'
	IF (n_sed_class==1) then
		write(fmtstr,'(a,i0,a)')'(a4,a,a4,',subasin,'(a,i14))'		!generate format string
		WRITE(11,fmtstr)'Year', char(9),'Day ', (char(9),id_subbas_extern(i),i=1,subasin)		!tab separated output
	ELSE
		write(fmtstr,'(a,i0,a)') '(a4,a,a4,',subasin*n_sed_class,'(a,i0,a,i0))'		!generate format string
		WRITE(11,fmtstr)'Year', char(9),'Day ', ((char(9),id_subbas_extern(i),':',j,j=1,n_sed_class),i=1,subasin) 
	END IF
	CLOSE(11)
  ELSE				!delete any existing file, if no output is desired
	CLOSE(11,status='delete')
  END IF
  
!     Output sub-saily sediment production (t)
  OPEN(11,FILE=pfadn(1:pfadi)//'sediment_production.out', STATUS='replace')
  IF (f_sediment_production .AND. dosediment .AND. dt<24) THEN	!only do output if file is desired, sediment modelling is enabled and model runs in sub-daily resolution
	WRITE(11,'(a,i2,a1)') 'total sediment production [t] for all sub-basins (MAP-IDs) and particle classes{',n_sed_class,'}'
	IF (n_sed_class==1) then
		write(fmtstr,'(a,i0,a)')'(A,',subasin,'i6)'		!generate format string
		WRITE(11,fmtstr)'Year	Day	Timestep',(id_subbas_extern(i),i=1,subasin) 
		!WRITE(11,'(A,<subasin>(i6))')'Year	Day	Timestep',(id_subbas_extern(i),i=1,subasin) 
	ELSE
		write(fmtstr,'(a,i0,a,i0,a)')'(A,',subasin,'(i6,a1,i1,',n_sed_class-1,'i3))'		!generate format string
		WRITE(11,fmtstr)'Year	Day	Timestep', (id_subbas_extern(i),':',(j,j=1,n_sed_class),i=1,subasin) 
		!WRITE(11,'(A,<subasin>(i6,a1,i1,<n_sed_class-1>i3))')'Year	Day	Timestep', (id_subbas_extern(i),':',(j,j=1,n_sed_class),i=1,subasin) 
	END IF
	CLOSE(11)
  ELSE				!delete any existing file, if no output is desired
	CLOSE(11,status='delete')
  END IF
  


!     debug Output for checking purposes !remove
  OPEN(11,FILE=pfadn(1:pfadi)//'debug.out', STATUS='replace')
	WRITE(11,*) 'current c-factor of lowermost TC, LU 22312.'
	debug_out=0
	CLOSE(11)
!     debug Output for checking purposes !remove
  OPEN(11,FILE=pfadn(1:pfadi)//'debug2.out', STATUS='replace')
	WRITE(11,*) 'current q	q_peak	ei	K_fac	C_fac	P_fac	LS_fac	CFRG_fac of lowermost TC, LU 31221.'
 	debug_out2=0
	CLOSE(11)

  



!Till: if any TC-specific output is desired, prepare file header
if (f_tc_theta .OR. f_tc_surfflow .OR. f_tc_sedout) then		

	!create tc identification key 
	dig_sub=max(1,ceiling(log10(1.*maxval(id_subbas_extern))))	!determine digits necessary to display the IDs
	dig_lu=max(1,ceiling(log10(1.*maxval(id_lu_extern))))
	dig_tc=max(1,ceiling(log10(1.*maxval(id_terrain_extern))))

	allocate(tc_idx(sv_comb))
	tc_idx=""

	write(fmtstr,'(a,i0,a,i0,a,i0,a,i0,a,i0,a,i0,a)'),'(I',dig_sub,'.',dig_sub,',I',dig_lu,'.',dig_lu,',I',dig_tc,'.',dig_tc,')'
	!formatstr for writing into tc_idx
	 tc_counter_all=1
	 DO i_subbas=1,subasin 
	  DO lu_counter=1,nbr_lu(i_subbas)
      i_lu=id_lu_intern(lu_counter,i_subbas)
		DO tc_counter=1,nbrterrain(i_lu)
		   !tc_idx(tc_counter_all)=(tc_counter)+(id_lu_extern(i_lu)*10**dig_tc)+(10**(dig_tc+dig_lu)*id_subbas_extern(i_subbas))
		   !write(tc_idx(tc_counter_all),fmtstr),id_subbas_extern(i_subbas),id_lu_extern(i_lu),tc_counter	!old scheme using TC-position
		   write(tc_idx(tc_counter_all),fmtstr),id_subbas_extern(i_subbas),id_lu_extern(i_lu),id_terrain_extern(id_terrain_intern(tc_counter,i_lu))	!new scheme using TC-id
		   tc_counter_all=tc_counter_all+1
		END DO
	  END DO
	 END DO
	 tc_counter_all=tc_counter_all-1
	
!	m=0		!sum up storage required
!	do i=1,length(year_print)	!compute number of timesteps that need to be saved
!		if (year_print(i)==-1) then	!save all model time
!			m=((tstop-tstart+1)*366-dayoutsim-(mstop-12)*30)*nt
!			exit
!		end if
!		if (year_print(i)==0) cycle	!save no data for this line
!	
!		if (month_print(i)==-1) then	!save all months of current year
!		
!			m=((tstop-tstart+1)*366-dayoutsim-(mstop-12)*30)*nt
!			exit
!		end if
!		if (year_print(i)==0) cycle	!save no data for this line
!
!		
!	end do


end if



!Output TC-wise soil moisture [%]
!conrad: tc-wise theta output file preparation
  OPEN(11,FILE=pfadn(1:pfadi)// 'tc_theta.out', STATUS='replace')
  IF (f_tc_theta) THEN	
  	 
	 !Till: allocate memory for TC-wise output
	 allocate(theta_tc(366,sv_comb)) 
	 theta_tc=-1 !debugging help

	!conrad: store average soil thickness for each TC (for each instance of a TC)
	allocate(meandepth_tc(subasin,maxsoter,maxterrain))
	meandepth_tc(:,:,:)=0

	 DO i_subbas=1,subasin 
	  DO lu_counter=1,nbr_lu(i_subbas)
      i_lu=id_lu_intern(lu_counter,i_subbas)
		DO tc_counter=1,nbrterrain(i_lu)
			DO i=1,nbr_svc(tc_counter)
				temp2=sum(horiz_thickness(tc_counter,i,:))	!compute depth of this soil type 
				meandepth_tc(i_subbas,lu_counter,tc_counter)=meandepth_tc(i_subbas,lu_counter,tc_counter)+temp2*frac_svc(i,tc_counter) !fraction of svc used for calculating weighted mean
			END DO
		END DO
	  END DO
	 END DO
	
	
	WRITE(11,'(A)'),'daily theta [%] for tcs in lus in sub-basins (scheme: '//REPEAT('S',dig_sub)//&
		REPEAT('L', dig_lu)//REPEAT('T', dig_tc)//') -> use with tc_plot.m in Matlab'
	!don't change headerline, needed by matlab- /R-script
	
	
	WRITE(fmtstr,'(a,i0,a,i0,a)') '(i0,a,i0,a,i0,',sv_comb,'(a,a',dig_sub+dig_lu+dig_tc,'))'
	WRITE(11,fmtstr)0,char(9),0,char(9),0,(char(9),tc_idx(i),i=1,tc_counter_all)		!tab separated output (TC indices as strings)

	CLOSE(11)
  ELSE				!delete any existing file, if no output is desired
	CLOSE(11,status='delete')
  END IF


OPEN(11,FILE=pfadn(1:pfadi)// 'tc_surfflow.out', STATUS='replace')
  IF (f_tc_surfflow) THEN	
	allocate(surfflow_tc(366,nt,sv_comb)) !Till: allocate memory for TC-wise output of surface flow
	WRITE(11,'(A)'),'surface runoff [mm] for tcs in lus in sub-basins (scheme: '//REPEAT('S',dig_sub)//&
		REPEAT('L', dig_lu)//REPEAT('T', dig_tc)//') -> use with tc_plot.m in Matlab'
	!don't change headerline, needed by matlab- /R-script
	
	WRITE(fmtstr,'(a,i0,a,i0,a)') '(i0,a,i0,a,i0,',sv_comb,'(a,a',dig_sub+dig_lu+dig_tc,'))'
	WRITE(11,fmtstr)0,char(9),0,char(9),0,(char(9),tc_idx(i),i=1,tc_counter_all)		!tab separated output (TC indices as strings)
	CLOSE(11)
  ELSE				!delete any existing file, if no output is desired
	CLOSE(11,status='delete')
  END IF

OPEN(11,FILE=pfadn(1:pfadi)// 'tc_sedout.out', STATUS='replace')
  IF (dosediment .AND. f_tc_sedout .AND. .NOT. do_musle_subbasin) THEN	
	allocate(sedout_tc(366,nt,sv_comb)) !Till: allocate memory for TC-wise output of sediment output
	WRITE(11,'(A)'),'sediment output [t/km²] for tcs in lus in sub-basins (scheme: '//REPEAT('S',dig_sub)//&
		REPEAT('L', dig_lu)//REPEAT('T', dig_tc)//') -> use with tc_plot.m in Matlab'
	!don't change headerline, needed by matlab- /R-script
	
	WRITE(fmtstr,'(a,i0,a,i0,a)') '(i0,a,i0,a,i0,',sv_comb,'(a,a',dig_sub+dig_lu+dig_tc,'))'
	WRITE(11,fmtstr)0,char(9),0,char(9),0,(char(9),tc_idx(i),i=1,tc_counter_all)		!tab separated output (TC indices as strings)
	CLOSE(11)
  ELSE				!delete any existing file, if no output is desired
	CLOSE(11,status='delete')
  END IF





if (allocated(tc_idx)) then		
	deallocate(tc_idx)		!TC-indices are no longer needed
end if


!!Print hydrologic variable on TC scale. If not used, DISABLE
!!*******************************************************************************
!!create files
!DO k=1,size(year_print)
!  WRITE(year_name,*)year_print(k)						!Till: assemble file name
!  WRITE(day_name,*)day_print(k)
!  DO c=1,12
!    IF (year_name(c:c) /= ' ') THEN
!      dummy1=c
!      EXIT
!	ENDIF
!  END DO
!  DO c=1,12
!    IF (day_name(c:c) /= ' ') THEN
!      dummy2=c
!      EXIT
!	ENDIF
!  END DO
!  OPEN(11,FILE=pfadn(1:pfadi)//'TC_'//year_name(dummy1:12)//'_'//day_name(dummy2:12)// &
!	'.out',STATUS='replace')
!   WRITE(11,*)'Subasin-ID, LU-ID, TC-ID, deposition'	
!  CLOSE(11)
!ENDDO
!!*******************************************************************************


END IF


!-------------------------------------------------------------------
IF (STATUS == 1) THEN
!  DO i_subbas=1,subasin		
!!** Initialization for one year
!!** Initialize seaflow variables
!    seaflow(i_subbas) = 0.		!Till: never used
!  END DO
  if (doacud) CALL lake(1,dummy)
!Till: set all arrays that contain (daily) data of the whole year to zero (all refer to subbasins?)
  sediment_subbasin_t=0.
  water_subbasin_t=0.

	water_subbasin=0.0		!Water contribution of each subbasin to the river per timestep
!George	lakeinflow=0.0			!river inflow
	qloss=0.0				!losses in river network by evaporation
	ovflow=0.0				!total surface runoff of subbasins for each day
	subflow=0.0				!total subsurface runoff of subbasins for each day

	!Till: these are mainly variables used for output that are not used for further computations (ii skip computation, when output is disabled)
	qgen=0.0				!river flow before acudes (small reservoirs)
!George	muniret=0.0
	gw_recharge=0.0			!groundwater recharge (percolation below root zone) !Till: into into linear storage [m3]
	!deepgw_r=0.0			!deep groundwater recharge (loss from model) 
	aet=0.0					!daily actual evapotranspiration (mm/day)
	laimun=0.0				!daily mean LAI (m²/m²)
	soilet=0.0				!daily soil evaporation
	intc=0.0				!daily interception storage evapotranspiration
	hortflow=0.0			!horton overland flow
	soilm=0.0				!soil moisture [mm]
	sofarea=0.0				!fraction of saturated area
	sediment_subbasin=0.

	soilmroot=0.0		!soil moisture first meter

	deep_gw_discharge=0.0	!ground water discharge
	gw_loss=0.0	!ground water loss from domain

	if (do_pre_outflow) then		!if water outflow from upstream subbasins is given
		call read_pre_subbas_outflow		!read
	end if

END IF


!-------------------------------------------------------------------
IF (STATUS == 2) THEN
  !** Calculation for one timestep


tc_counter_all=1		!reset TC counter
  DO i_subbas=1,subasin
	
	j=0	!Till: temporary indicator if this subbasin can be skipped because of prespecified values

	if (do_pre_outflow) then		!if water outflow from upstream subbasins is given
		if (corr_column_pre_subbas_outflow(i_subbas)>0) then			!if outflow of subbasin is prespecified
			water_subbasin (d,i_subbas) = sum(pre_subbas_outflow(d,1:nt,corr_column_pre_subbas_outflow(i_subbas)))
			water_subbasin_t(d,1:nt,i_subbas)=pre_subbas_outflow(d,1:nt,corr_column_pre_subbas_outflow(i_subbas))	
			j=1	!Till: indicate that this subbasin can be skipped because of prespecified values
			if (dosediment) then
				sediment_subbasin(d,i_subbas,:) = -1.	!Till: this is for clearly indicating missing data
				sediment_subbasin_t(d,1:nt,i_subbas,:)=-1.
			end if
		end if
	end if !do outflow
	
	if (do_pre_outsed) then		!if water outsed from upstream subbasins is given
		if (corr_column_pre_subbas_outsed(i_subbas)>0) then			!if outsed of subbasin is prespecified
			sediment_subbasin(d,i_subbas,:) = sum(pre_subbas_outsed(d,1:nt,corr_column_pre_subbas_outsed(i_subbas)))*pre_psd
			do i=1,n_sed_class
				sediment_subbasin_t(d,1:nt,i_subbas,i)=pre_subbas_outsed(d,1:nt,corr_column_pre_subbas_outsed(i_subbas))*pre_psd(i)	
			end do
			dosediment=.FALSE.		!don't compute sediment for this subbasin because it is pre-specified
		else
			dosediment=.TRUE.
		end if
	end if !do outsed
    

    if (j==1) cycle		!Till: indicator: this subbasin can be skipped because of prespecified values
    
	precday=precip(d,i_subbas)		!daily precip

    julian_day=d
	IF (t == tstart) THEN
	  julian_day=julian_day+sum(daymon(1:mstart-1))
	ENDIF 
	
	IF (MOD((t),4) == 0)  THEN	!compute julian day
	  julian_day=julian_day+1
	ENDIF


!  determine actual vegetation characteristics
!  (height, root depth, LAI, albedo) for all vegetation units
	rootd_act =calc_seasonality(i_subbas,t,julian_day,period(i_subbas,:),rootdep)	!compute root depths of current day
	height_act=calc_seasonality(i_subbas,t,julian_day,period(i_subbas,:),height)	!compute heights of current day
	lai_act   =calc_seasonality(i_subbas,t,julian_day,period(i_subbas,:),lai)	!compute LAIs of current day
	alb_act   =calc_seasonality(i_subbas,t,julian_day,period(i_subbas,:),alb)	!compute albedos of current day

	svc_k_fac_day     =calc_seasonality(i_subbas,t,julian_day,seasonality_k     (i_subbas,:),svc_k_fac)	        !compute K-factors of current day
	svc_c_fac_day     =calc_seasonality(i_subbas,t,julian_day,seasonality_c     (i_subbas,:),svc_c_fac)	        !compute c-factors of current day
	svc_p_fac_day     =calc_seasonality(i_subbas,t,julian_day,seasonality_p     (i_subbas,:),svc_p_fac)	        !compute p-factors of current day
	svc_coarse_fac_day=calc_seasonality(i_subbas,t,julian_day,seasonality_coarse(i_subbas,:),svc_coarse_fac)	!compute coarse-factors of current day
	svc_n_day         =calc_seasonality(i_subbas,t,julian_day,seasonality_n     (i_subbas,:),svc_n)	            !compute n of current day

	kfkorr_day=kfkorr*(kfkorr_a*1/precip(d,i_subbas)+kfkorr_b)	!compute kfkorr as a function of daily precipitation
	

	deepgwrsu=0.


! LOOPs through all Landscape Units of specific sub-basin
    DO lu_counter=1,nbr_lu(i_subbas)
      i_lu=id_lu_intern(lu_counter,i_subbas)
      
      balance=0.0		
      frac_satsu=0.0

	  fractemp(:)=0.
      fractemp(1:nbrterrain(i_lu))=fracterrain(id_terrain_intern(1:nbrterrain(i_lu),i_lu))	!Till: initialise auxiliary array that holds the saturated fraction of each TC in the current LU
      


!      IF (dohour) THEN			
!only for output, still left in the code, may be removed later
!		  sat_xs_overland_flow_sc(d,:)=0.
!		  sat_area_of_sc(d,:)=0.
!		  hortsc(d,:)=0.
!		  hort2sc(d,:)=0.
!		  gwrsc(d,:)=0.
!		  deepgwrsc(d,:)=0.
!		  thsc(:,d,:)=0.
!		  resistsc(d,:)=0
!		  aetsc(d,:)=0.
!		  aettotsc(d,:)=0.
!		  intetsc(d,:)=0.
!		  soiletsc(d,:)=0.
!		  nfksc(d,:)=0.
!	  END IF



! LOOPs through required number of timesteps (for current day)
      DO timestep_counter=1,nt

        qsurf_lu(lu_counter)=0.		!Till: reset variables for summing up at the LU-scale
		qsub_lu(lu_counter)=0.
		gwrsu(lu_counter)=0.
		hortsu(lu_counter)=0.
		aetsu(lu_counter)=0.
		laisu(lu_counter)=0.
		soiletsu(lu_counter)=0.
		intcsu(lu_counter)=0.
		soilmsu(lu_counter)=0.
		sedsu(lu_counter,:)=0.
		
		surfflow_in(:)=0.			!Till: reset TC-related variables
        surfflow_out(:)=0.
        sublat_in(:)=0
        sublat_out(:)=0.
		sed_in(:,:)=0.
		sed_out(:,:)=0.

!		!Till: reset TC-instance-related variables
!		horttc  (:)=0.
!		aettc   (:)=0.
!		laitc   (:)=0.
!		soilettc(:)=0.
!		intctc  (:)=0.


		
		soilmrootsu(lu_counter)=0.		
		
! LOOPs through all Terrain Components of specific Landscape Unit (of specific sub-basin)

		L_cum=0	!Pedro - zero the cumulated slope length for erosion calculations
  
        DO tc_counter=1,nbrterrain(i_lu)
          tcid_instance=tcallid(i_subbas,lu_counter,tc_counter)
          
		  id_tc_type=id_terrain_intern(tc_counter,i_lu)
		  tcarea=area(i_subbas)*frac_lu(lu_counter,i_subbas)* fracterrain(id_tc_type)

		  IF (timestep_counter == 1) THEN
			horttc   (tcid_instance)=0.				!Till: reset TC-instance-related variables - these are used to compute daily output values
			aettc    (tcid_instance)=0.
			laitc    (tcid_instance)=0.
			soilettc (tcid_instance)=0.
			intctc   (tcid_instance)=0.
			thact=soilwater(dprev,tcid_instance)

			gwrtc    (tcid_instance)=0.				!ii zero these entries at start of each timestep for entire domain 
			deepgwrtc(tcid_instance)=0.
		  END IF

		  
!!Print hydrologic variable on TC scale. If not used, DISABLE
!!**************************************************************************************************
!	      IF (precip(d,i_subbas) == 0.) deposition_TC(i_subbas,id_tc_type)=0.
!	      IF (precip(d,i_subbas) /= 0.) deposition_TC(i_subbas,id_tc_type)=1.
!!**************************************************************************************************

		  gwr=0.
		  deepgwr=0.
		  
		  IF (.NOT. dohour) THEN	!ii: join daily and hourly branch
			prec=precip(d,i_subbas)
			hh=0

	        CALL soilwat(hh,d,m,i_subbas,i_subbas,i_lu,lu_counter,tcid_instance,id_tc_type,tc_counter,thact,  &
            thactroot,surfflow_in(tc_counter),surfflow_out(tc_counter),  &
            sublat_in(tc_counter),sublat_out(tc_counter),gwr,deepgwr, horttc(tcid_instance),  &
            aettc(tcid_instance), laitc(tcid_instance),  &
            soilettc(tcid_instance), intctc(tcid_instance),  &
            prec,precday,prechall, pet(d,i_subbas),  &
            tcarea,balance_tc, rootd_act,height_act,lai_act,alb_act,sed_in(tc_counter,:),sed_out(tc_counter,:))
			balance=balance+balance_tc*fracterrain(id_tc_type)	!Till: compute water balance for LU [mm]
			
			soilwater(d,tcid_instance)=thact		!Till: save these values, because they were written to vars that will be recycled
			gwrtc(tcid_instance)=gwr			
			deepgwrtc(tcid_instance)=deepgwr	 
			
			if ( (tc_counter==nbrterrain(i_lu)) .AND. (i_lu==1)) then		!debugging output !remove
			 !debug_out(d)=aettc(tcid_instance)
			 !debug_out(d)=svc_k_fac_day(6)
			 !debug_out(d)=lai_act(1)
			 debug_out(d)=svc_c_fac_day(1)
			end if



		  ELSE
			aeth=0.
			laih=0.
			soileth=0.
			inth=0.

			prec=preciph((d-1)*24+timestep_counter,i_subbas)
			prechall(1:24)=preciph((d-1)*24+1:d*24,i_subbas)

			CALL soilwat(timestep_counter,d,m,i_subbas,i_subbas,i_lu,lu_counter,tcid_instance,id_tc_type,tc_counter,  &
			thact,thactroot,surfflow_in(tc_counter),  &
			surfflow_out(tc_counter),sublat_in(tc_counter),sublat_out(tc_counter),  &
			gwr,deepgwr,horttc(tcid_instance),aeth,laih,soileth,inth,  &
			prec,precday,prechall(1:24), pet(d,i_subbas),  &
			tcarea,balance_tc, rootd_act,height_act,lai_act,alb_act,sed_in(tc_counter,:),sed_out(tc_counter,:))
			balance=balance+balance_tc*fracterrain(id_tc_type)	!Till: compute water balance for LU [mm]
			
			gwrtc    (tcid_instance)=gwrtc    (tcid_instance)+gwr
			deepgwrtc(tcid_instance)=deepgwrtc(tcid_instance)+deepgwr	 
			aettc    (tcid_instance)=aettc    (tcid_instance)+aeth
			laitc    (tcid_instance)=laitc    (tcid_instance)+laih/24.
			soilettc (tcid_instance)=soilettc (tcid_instance)+soileth
			intctc   (tcid_instance)=intctc   (tcid_instance)+inth

			IF (timestep_counter == 24) THEN
				soilwater(d,tcid_instance)=thact
				soilmrootsu(lu_counter)=soilmrootsu(lu_counter)+ thactroot*fracterrain(id_tc_type)
			END IF

		  END IF	!dohour


		  !conrad: convert thact (in mm) into percent of soil volume using mean soil depth in tc
		  if (f_tc_theta) then
			!theta_tc(d,tc_counter_all)=thact/meandepth_tc(i_subbas,lu_counter,tc_counter) !original version, using mean theta for entire profile (instead of topmost horizon only, see below)
			
			theta_tc(d,tc_counter_all)=0
			DO i=1,nbr_svc(tcid_instance)
				rtemp=horithact(tcid_instance,i,1) / (thetas(id_soil_intern(i,tcid_instance),1)* horiz_thickness(tcid_instance,i,1))!Till: compute relative saturation of topmost horizon
				theta_tc(d,tc_counter_all)=theta_tc(d,tc_counter_all)+rtemp*frac_svc(i,tcid_instance)	!Till: compute weighted average according to fraction of SVC
			END DO	
		  end if

		  if (f_tc_surfflow) then	!Till: safe TC-wise surface runoff
			surfflow_tc(d,timestep_counter,tc_counter_all)=surfflow_out(tc_counter)/(tcarea*1e3)	
		  end if
		  
		  if (dosediment .AND. f_tc_sedout .and. .NOT. do_musle_subbasin) then	!Till: safe TC-wise sediment output
			sedout_tc(d,timestep_counter,tc_counter_all)=sum(sed_out(tc_counter,:))/tcarea	
		  end if
		  
		  tc_counter_all=tc_counter_all+1	!Till: for TC-wise output

!  lateral flow between terrain components (TCs)
!  lateral outflow is lateral inflow to deeper terrain compent
!  if dolattc is set
!  - if doalllattc is set, all outflow from one tc_counter is inflow to lower tc_counter
!  - if doalllattc is not set, outflow is subdivided between lower TCs and
!    direct runoff into river according to areal fractions of TCs in
!    landscape unit
				!Till: this implicitly accounts for preferential flow
!  this refers to surface flow only, subsurface flow is always completely
!  routed to next downstream TC if dolattc is set
        
			IF (tc_counter < nbrterrain(i_lu)) THEN

			  IF (dolattc) THEN
				sublat_in(tc_counter+1)  =sublat_in(tc_counter+1)+sublat_out(tc_counter)
				IF (doalllattc) THEN
				  surfflow_in(tc_counter+1)=surfflow_in(tc_counter+1)+surfflow_out(tc_counter)
				  hortsu(lu_counter)=hortsu(lu_counter)
				  sed_in(tc_counter+1,:)=sed_out(tc_counter,:)	!Till: all sediment leaving the upper TC enters the next downslope TC
													! prepare sed_in this for the next downslope TC
				ELSE
				  temp2=sum(fractemp(tc_counter:nbrterrain(i_lu))) !fraction of all remaining TCs in this LU 
				  DO i=tc_counter+1,nbrterrain(i_lu)				!Till: overland flow is distributed among lower TCs proportional to area
					surfflow_in(i)=surfflow_in(i)  +surfflow_out(tc_counter  )* fractemp(i)/temp2
					sed_in(i,:)   =sed_in     (i,:)+sed_out     (tc_counter,:)* fractemp(i)/temp2		
																			!Till: sediment goes directly to river and into next downslope TC (no redistribution among all more-downslope TCs as with water)
				  END DO
	!  Remaining outflow of each terrain component, which is not routed to lower TCs
	!  is surface runoff of entire landscape unit (direct dunoff to river)
				  qsurf_lu(lu_counter  )=qsurf_lu(lu_counter  )+surfflow_out(tc_counter)   * fractemp(tc_counter)/temp2
				  hortsu  (lu_counter  )=hortsu  (lu_counter  )+horttc      (tcid_instance)* fractemp(tc_counter)/temp2
				  sedsu   (lu_counter,:)=sedsu   (lu_counter,:)+sed_out     (tc_counter,:) * fractemp(tc_counter)/temp2
							!Till: remaining sediment that is not routed to the lower TC reaches the river directly
				END IF
			  ELSE IF (.NOT. dolattc) THEN		!Till: no flow between TCs, all flows leave the LU directly
				qsub_lu   (lu_counter)=qsub_lu(lu_counter)+sublat_out(tc_counter)

				qsurf_lu(lu_counter  )=qsurf_lu(lu_counter)  +surfflow_out(tc_counter)
				hortsu  (lu_counter  )=hortsu  (lu_counter)  +horttc      (tcid_instance)
				sedsu   (lu_counter,:)=sedsu   (lu_counter,:)+sed_out     (tc_counter,:)
			  END IF

	!  Outflow of lowest terrain component is added to runoff generated
	!  in landscape unit of this subbasin
			ELSE IF (tc_counter == nbrterrain(i_lu)) THEN	!lowest TC is reached, all flows leave the LU ii: join this with previous branch
			  qsurf_lu(lu_counter  )=qsurf_lu(lu_counter  )+surfflow_out(tc_counter)
			  qsub_lu (lu_counter  )=qsub_lu (lu_counter  )+sublat_out  (tc_counter)
			  sedsu   (lu_counter,:)=sedsu   (lu_counter,:)+sed_out     (tc_counter,:)
			  hortsu  (lu_counter  )=hortsu  (lu_counter  )+horttc      (tcid_instance)
			END IF
			
		
			IF (.NOT. dohour) THEN
				soilmsu(lu_counter)=soilmsu(lu_counter)+ thact*fracterrain(id_tc_type)	!compute mean soil moisture in LU  
			END IF

		
		END DO !end of calculations for each terrain component in the landscape unit

		!Till: store sub-daily values in extra array
		water_subbasin_t(d,timestep_counter,i_subbas)  = water_subbasin_t(d,timestep_counter,i_subbas) + &
			(qsurf_lu(lu_counter)+qsub_lu(lu_counter))
		!old: water_subbasin_t(d,timestep_counter,i_subbas)  = water_subbasin_t(d,timestep_counter,i_subbas) + (1.-intercepted)*(qsurf_lu(lu_counter)+qsub_lu(lu_counter))
		
		!George: disable subsurface flow (for debugging) - outcomment previous line, include following line
		!water_subbasin_t(d,timestep_counter,i_subbas)  = water_subbasin_t(d,timestep_counter,i_subbas) + (qsurf_lu(lu_counter))
		!old: water_subbasin_t(d,timestep_counter,i_subbas)  = water_subbasin_t(d,timestep_counter,i_subbas) + (1.-intercepted)*(qsurf_lu(lu_counter))
		
		
		IF (dosediment .AND. allocated(sdr) ) THEN	!apply prespecified SDR, if present
			sedsu(lu_counter,:)=sdr(i_lu)*sedsu(lu_counter,:)
			!sedsu(lu_counter,:)=0.*sedsu(lu_counter,:) !test option
		END IF
		
		sediment_subbasin_t(d,timestep_counter,i_subbas,:)=sediment_subbasin_t(d,timestep_counter,i_subbas,:)+sedsu(lu_counter,:) 
	  END DO
	  !  END of hourly calculations
      
      
	  !  calculate (daily) water balance variables at landscape unit scale

	  DO tc_counter=1,nbrterrain(i_lu)
		tcid_instance=tcallid(i_subbas,lu_counter,tc_counter)
		id_tc_type=id_terrain_intern(tc_counter,i_lu)
		
		aetsu(lu_counter)=aetsu(lu_counter)+ aettc(tcid_instance)*fracterrain(id_tc_type)
		laisu(lu_counter)=laisu(lu_counter)+ laitc(tcid_instance)*fracterrain(id_tc_type)
		soiletsu(lu_counter)=soiletsu(lu_counter)+ soilettc(tcid_instance)*fracterrain(id_tc_type)
		intcsu(lu_counter)=intcsu(lu_counter)+ intctc(tcid_instance)*fracterrain(id_tc_type)
	
		deepgwrsu(lu_counter)=deepgwrsu(lu_counter) +deepgwrtc(tcid_instance)*fracterrain(id_tc_type) !new
		gwrsu(lu_counter)=gwrsu(lu_counter)+gwrtc(tcid_instance)*fracterrain(id_tc_type)
		
		IF (dohour) THEN
			!gwrsu(lu_counter)=gwrsu(lu_counter)+gwrtc(tcid_instance)*fracterrain(id_tc_type)
			!deepgwrsu(lu_counter)=deepgwrsu(lu_counter) +deepgwrtc(tcid_instance)*fracterrain(id_tc_type)
			soilmsu(lu_counter)=soilmsu(lu_counter)+ thact*fracterrain(id_tc_type)	!Till: use last value of day to compute mean soil moisture in LU 
																					!thact ist nur von letzter TC in LU. Müsste das nicht irgendwie in die TC-Schleife mit rein?
		ELSE
			!gwrsu(lu_counter)=gwrsu(lu_counter)+gwr*fracterrain(id_tc_type)	!falsch: gwr und deepgwr müssen zwischengespeichert werden (wie in Stundenversion)
			!deepgwrsu(lu_counter)=deepgwrsu(lu_counter)+deepgwr*fracterrain(id_tc_type)
			soilmrootsu(lu_counter)=soilmrootsu(lu_counter)+ thactroot*fracterrain(id_tc_type)	!ii this is most likely wrong because thactroot is from last TC only
			frac_satsu=frac_satsu+sum(frac_sat(tcid_instance,:))* fracterrain(id_tc_type)
		END IF
	  END DO	!thru all TCs of current LU

! deep groundwater runoff of each landscape unit (daily budget)
     IF (gw_flag(i_lu) == 1 .AND. gw_dist(i_lu) <= 1.) THEN		!Till: if gw is enabled and normal case (gw_dist(i_lu)/= 99)... 
																	!Till: ii gw_dist is zero anyway, eliminate it

		tcid_instance=tcallid(i_subbas,lu_counter,nbrterrain(i_lu))	!Till: get ID of lowermost TC in current LU
        
		rtemp=deepgwrsu(lu_counter)* (1.-gw_dist(i_lu))*  &
            area(i_subbas)*frac_lu(lu_counter,i_subbas)*1.e3	!Till: compute ground water recharge [m3]
		deepgw(i_subbas,lu_counter)=deepgw(i_subbas,lu_counter)+rtemp		!Till: add gw recharge to gw storage [m3]		

		gw_recharge(d,i_subbas)= gw_recharge(d,i_subbas)+rtemp							!Till: for output of gw-recharge


		!rtemp=0. !George: disable groundwater contribution: enable this line, outcomment next line
		rtemp=min(deepgw(i_subbas,lu_counter),deepgw(i_subbas,lu_counter)/gw_delay(i_lu))		!Till: compute gw-discharge [m3], at maximum this equals the entire stored volume
		
		if ((rtemp>0) .AND. (frac_direct_gw<1)) then !Till: ground water discharge is (partially) routed into interflow of lowermost tc
			dummy=min(INT(maxval(sum(horiz_thickness(tcid_instance,:,:),2))/500.)+1,size(latred, DIM = 2)  ) !Till: compute number of horizon layers needed for storing the subsurface flow
			latred(tcid_instance,1:dummy)=latred(tcid_instance,1:dummy)+ rtemp*(1-frac_direct_gw)/dummy
			!distribute gw-fraction routed to subsurface flow of lowest TC equally among soil horizons (to be redistributed in next timestep)
		end if
		rtemp2=rtemp*frac_direct_gw	!Till: fraction directly routed to the river 

		qsub_lu(lu_counter)=qsub_lu(lu_counter)+rtemp2		!Till: add gw-discharge to total subsurface runoff
		!old water_subbasin_t(d,:,i_subbas)=water_subbasin_t(d,:,i_subbas)+(1.-intercepted)*rtemp2/nt	!Till: deep gw discharge is distributed equally among all timesteps of day
		water_subbasin_t(d,:,i_subbas)=water_subbasin_t(d,:,i_subbas)+rtemp2/nt	!Till: deep gw discharge is distributed equally among all timesteps of day
		deepgwrsu(lu_counter)=deepgwrsu(lu_counter)-deepgwrsu(lu_counter)*(1.-gw_dist(i_lu))	!Till: must be zeroed, everything still in there leaves the model domain
		deepgw(i_subbas,lu_counter)=deepgw(i_subbas,lu_counter)- rtemp				  			!Till: gw storage is reduced according to outflow (direct outflow to river and outflow to lowest TC)
		
		IF (f_deep_gw_discharge) THEN	!Till: sum up daily groundwater discharge (effective part into river) (for output only)
			deep_gw_discharge(d,i_subbas)=deep_gw_discharge(d,i_subbas)+rtemp2		!(1.-intercepted) 
		END IF

		
		!qsub_lu(lu_counter)=qsub_lu(lu_counter)+deepgw(i_subbas,lu_counter)/gw_delay(i_lu)		!Till: compute gw discharge from subbasin
        !water_subbasin_t(d,:,i_subbas)=water_subbasin_t(d,:,i_subbas)+(1.-intercepted)*deepgw(i_subbas,lu_counter)/gw_delay(i_lu)/nt	!Till: deep gw discharge is distributed equally among all timesteps of day
		!deepgwrsu(lu_counter)=deepgwrsu(lu_counter)-deepgwrsu(lu_counter)*(1.-gw_dist(i_lu))	!gw storage is reduced according to outflow (losses from model domain)
        !deepgw(i_subbas,lu_counter)=deepgw(i_subbas,lu_counter)- deepgw(i_subbas,lu_counter)/gw_delay(i_lu) !gw storage is reduced according to outflow
      END IF

!		!Till: ground water is disabled: seepage below profile is lost
!		IF (gw_flag(i_lu) == 0 ) THEN	
!
!		END IF

!  runoff generated in landscape units is simply summed to give runoff of
!  entire subbasin
!  ETP, horton flow. gwr and soil moisture is area-weighted mean (mm)
      soilm(d,i_subbas)=   soilm(d,i_subbas)+ soilmsu(lu_counter)*frac_lu(lu_counter,i_subbas)		!ii compute values for each LU, do final summing up outside loop using 'sum'
      soilmroot(d,i_subbas)=soilmroot(d,i_subbas)+  soilmrootsu(lu_counter)*frac_lu(lu_counter,i_subbas)
      aet(d,i_subbas)=     aet(d,i_subbas)+aetsu(lu_counter)* frac_lu(lu_counter,i_subbas)
      laimun(d,i_subbas)=     laimun(d,i_subbas)+laisu(lu_counter)* frac_lu(lu_counter,i_subbas)
      soilet(d,i_subbas)=  soilet(d,i_subbas)+soiletsu(lu_counter)* frac_lu(lu_counter,i_subbas)
      intc(d,i_subbas)=    intc(d,i_subbas)+intcsu(lu_counter)* frac_lu(lu_counter,i_subbas)
      hortflow(d,i_subbas)=hortflow(d,i_subbas)+hortsu(lu_counter)
      ovflow(d,i_subbas)=  ovflow(d,i_subbas)+qsurf_lu(lu_counter)
	  rtemp=ovflow(d,i_subbas)
	  subflow(d,i_subbas)= subflow(d,i_subbas)+qsub_lu(lu_counter)
      !gw_recharge(d,i_subbas)= gw_recharge(d,i_subbas)+gwrsu(lu_counter)*  frac_lu(lu_counter,i_subbas)
      !deepgw_r(d,i_subbas)=    deepgw_r(d,i_subbas)+deepgwrsu(lu_counter)* frac_lu(lu_counter,i_subbas)*area(i_subbas)*1.e3 !Till: sum up ground water recharge, convert to m3
      sofarea(d,i_subbas)= sofarea(d,i_subbas)+ frac_satsu*frac_lu(lu_counter,i_subbas)
      qgen(d,i_subbas)=    qgen(d,i_subbas)+ (qsurf_lu(lu_counter)+qsub_lu(lu_counter))
	  sediment_subbasin(d,i_subbas,:)= sediment_subbasin(d,i_subbas,:)+sedsu(lu_counter,:)
    END DO
!  this is the end of calculation for each landscape unit in current cell/subbasin
    
!  runoff generated in cells units is simply summed to give runoff of
!  entire watershed
!  ETP, horton flow. gwr and soil moisture is area-weighted mean (mm)
 
IF (dosediment .AND. do_musle_subbasin .AND. (ovflow(d,i_subbas)>0.) ) THEN	!calculate sediment yield on subbasin scale (in contrast to on TC-scale)
	!sediment_subbasin(d,i_subbas,:)=sedi_yield_subbas(subbas_id, ovflow(d,i_subbas), sed_yield_subbas)
	CALL sedi_yield_subbas(i_subbas, ovflow(d,i_subbas), sediment_subbasin(d,i_subbas,:))
	IF (nt==1)sediment_subbasin_t(d,1,i_subbas,:)=sediment_subbasin(d,i_subbas,:) !George
END IF

 
    
!George *****************************************************  
! routing between landscape units
! currently, simple summing up of runoff of all LUs
! inflow to acudes of lowlands: lakeinflow

!George if (doacud) then
!George    lakeinflow(d,i_subbas) = lakeinflow(d,i_subbas)+ intercepted*ovflow(d,i_subbas)+  &
!George        intercepted*subflow(d,i_subbas)
!George endif
!George ***************************************************** 
    
!** water_subbasin composed of components:
!    1) remaining surface flow generated in the lowlands
!    2) outflow from largest volume class of small acudes is added later
!    3) groundwater exfiltration in the lowlands downstream small acudes
	
!	!threshold value for eliminating minimum outflow
!	rtemp=(1.-intercepted)*(ovflow(d,i_subbas)+ subflow(d,i_subbas))
!	if ((rtemp>0) .AND. (rtemp <= 0.000*3600*24*area(i_subbas))) then				!Till: effective discharge only if more than 0.001 m3/s/km2 is exceeded (prevents very low flows)
!			deepgw(i_subbas,1:nbr_lu(i_subbas))=deepgw(i_subbas,1:nbr_lu(i_subbas))+ rtemp/nbr_lu(i_subbas)*frac_lu(1:nbr_lu(i_subbas),i_subbas) 
!			!if water yield is very low, it goes into the linear storage (distributed among LUs according to areal fraction) instead of the river
!					
!			water_subbasin_t(d,:,i_subbas)=0	!Till: no water reaches the river
!			water_subbasin(d,i_subbas)  = 0		
!		
!			hortflow(d,i_subbas)=0
!			ovflow(d,i_subbas)=  0
!			subflow(d,i_subbas)= 0
!			sediment_subbasin(d,i_subbas,:)= 0
!			gw_recharge(d,i_subbas)= rtemp
!			qgen(d,i_subbas)=    0
!
!		else
!			water_subbasin(d,i_subbas)  = rtemp			!normal case: all water flux components are directed to the river
!	end if
    
!	water_subbasin(d,i_subbas)  = water_subbasin(d,i_subbas)+(1.-intercepted)*(ovflow(d,i_subbas)+ subflow(d,i_subbas))
	water_subbasin(d,i_subbas)  = water_subbasin(d,i_subbas)+(ovflow(d,i_subbas)+ subflow(d,i_subbas))


!George	water_subbasin(d,i_subbas)  = (1.-intercepted)*(ovflow(d,i_subbas))
!write(*,*)water_subbasin (d,i_subbas),water_subbasin_t(d,1,i_subbas)
	
!  this is the end of calculations for each sub-basin
    
!  correct runoff volumes by fraction of actual lake surface area
!  on total watershed area
!  (which receives direct precipitation)
!  assessment of total lake surface in watershed:
!  (small lake area given in km**2, largedam area in m**2)
    
!George new values for the variables water_subbasin_t and sediment_subbasin_t are obtained
! directly in the lake.f90
    IF (doacud) THEN
      CALL lake(2,i_subbas)
	ELSE 
		IF (doreservoir) THEN
			water_subbasin (d,i_subbas)=water_subbasin(d,i_subbas)*(1-((damareaact(i_subbas)/1.e6)/area(i_subbas)))	  
		END IF
	END IF
	
	water_subbasin (d,i_subbas) = water_subbasin (d,i_subbas) / (24*3600)	!convert m**3/d into m**3/s   
	water_subbasin_t(d,:,i_subbas)=water_subbasin_t(d,:,i_subbas) / (dt*3600)	!convert m**3 into m**3/s   

  
	if (f_gw_loss) THEN
		gw_loss(d,i_subbas)=sum(deepgwrsu(1:(nbr_lu(i_subbas)))*frac_lu(1:(nbr_lu(i_subbas)),i_subbas) )*area(i_subbas)*1.e3		!Till: all deep seepage that has not been transferred to groundwater is lost from model domain, convert to m^3
	end if
	    
!!Print hydrologic variable on TC scale. If not used, DISABLE
!!************************************************************************
!	DO lu_counter=1,nbr_lu(i_subbas)
!      i_lu=id_lu_intern(lu_counter,i_subbas)
!      DO tc_counter=1,nbrterrain(i_lu)
!		id_tc_type=id_terrain_intern(tc_counter,i_lu)
!
!	    DO k=1,6
!	      if (t==year_print(k) .and. d+dayoutsim ==day_print(k)) then
!            WRITE(year_name,*)year_print(k)						!Till: assemble file name
!            WRITE(day_name,*)day_print(k)
!            DO c=1,12
!              IF (year_name(c:c) /= ' ') THEN
!                dummy1=c
!                EXIT
!			  ENDIF
!            END DO
!            DO c=1,12
!              IF (day_name(c:c) /= ' ') THEN
!                dummy2=c
!                EXIT
!			  ENDIF
!            END DO
!	        OPEN(11,FILE=pfadn(1:pfadi)//'TC_'//year_name(dummy1:12)//'_'//day_name(dummy2:12)// &
!				'.out',STATUS='old',POSITION='append')
!				write(11,'(3I10,F12.3)')id_subbas_extern(i_subbas),id_lu_extern(i_lu),id_terrain_extern(id_tc_type),&
!					deposition_TC(i_subbas,id_tc_type)
!			CLOSE(11)
!		  END IF
!	    ENDDO
!	  ENDDO
!	ENDDO
!!*************************************************************************
    
!   river runoff of each sub-basin is
!   used in subroutine "routing",incl. water balance of large reservoirs
    
!   end of loop for all sub-basins
  END DO
 
	if (do_pre_outsed) then		!if water outsed from upstream subbasins is given
		dosediment=.TRUE.		!this may have been switched off temproarily if the last subbas was a prespecified one - switch it on again
	end if
    
END IF

!----------------------------------------------------------------------
IF (STATUS == 3) THEN
  
! Output daily water contribution into river (m**3/s and m**3)
  IF (f_daily_water_subbasin) THEN	!if output file is enabled
	  OPEN(11,FILE=pfadn(1:pfadi)//  &
		  'daily_water_subbasin.out', STATUS='old',POSITION='append')
	  write(fmtstr,'(a,i0,a)') '(i0,a,i0,',subasin,'(a,f11.6))'		!generate format string 
	  DO d=1,dayyear
		WRITE (11,fmtstr)t, char(9), d, (char(9),water_subbasin(d,i),i=1,subasin)	!in m**3/s
	  END DO
	  CLOSE(11)
  	  
	  OPEN(11,FILE=pfadn(1:pfadi)//  &
		  'daily_water_subbasin2.out', STATUS='old',POSITION='append')
	  write(fmtstr,'(a,i0,a)') '(i0,a,i0,',subasin,'(a,f11.6))'		!generate format string 
	  DO d=1,dayyear
		WRITE (11,fmtstr)t, char(9), d, (char(9),water_subbasin(d,i)*3600*24,i=1,subasin) !in m**3
	  END DO
	  CLOSE(11)
  END IF


	CALL write_output(f_daily_actetranspiration,'daily_actetranspiration.out',aet)	! Output daily actual evapotranspiration
	CALL write_output(f_daily_potetranspiration,'daily_potetranspiration.out',pet)	! Output daily potential evapotranspiration
	CALL write_output(f_daily_theta,'daily_theta.out',soilm)		! Output daily soil moisture
	CALL write_output(f_daily_qhorton,'daily_qhorton.out',hortflow)		! Output daily Hortonian overland flow
	CALL write_output(f_daily_total_overlandflow,'daily_total_overlandflow.out',ovflow)			! Output daily total overland flow
	CALL write_output(f_deep_gw_recharge,'deep_gw_recharge.out',gw_recharge)			!   deep groundwater recharge (loss from modell domain / into linear storage)
	CALL write_output(f_deep_gw_discharge,'deep_gw_discharge.out',deep_gw_discharge)  !   deep groundwater discharge (component of water_subbasin)
	CALL write_output(f_gw_loss,'gw_loss.out',gw_loss)			!   groundwater losses (leaving model domain)
	CALL write_output(f_daily_subsurface_runoff,'daily_subsurface_runoff.out',subflow)		!     Output total subsurface runoff (m³/d)



! Output sub-daily water flux   
  IF (f_water_subbasin) THEN	!Till: if sub-daily resolution is required
	OPEN(11,FILE=pfadn(1:pfadi)//'water_subbasin.out', STATUS='old',POSITION='append')
	write(fmtstr,'(a,i0,a)') '(i0,a,i0,a,i0,',subasin,'(a,f11.6))'		!generate format string
	DO d=1,dayyear
		DO j=1,nt
			WRITE (11,fmtstr)t, char(9), d, char(9), j, (char(9),water_subbasin_t(d,j,i),i=1,subasin) !in m**3/s
		END DO
	END DO
	CLOSE(11)
  END IF


!     Output sediment production (t)
  IF (f_daily_sediment_production .AND. dosediment) THEN
	OPEN(11,FILE=pfadn(1:pfadi)//'daily_sediment_production.out',  &
		  STATUS='old',POSITION='append')
	write(fmtstr,'(a,i0,a)') '(i0,a,i0,',subasin*n_sed_class,'(a,f13.4))'		!generate format string
	DO d=1,dayyear
		WRITE(11,fmtstr)t, char(9),d, ((char(9),sediment_subbasin(d,i,j),j=1,n_sed_class),i=1,subasin)	!tab separated output
	END DO
	CLOSE(11)
  END IF


!     Output sub-daily sediment production (t)
  IF (f_sediment_production .AND. dosediment .AND. dt<24) THEN	!only do output if file is desired, sediment modelling is enabled and model runs in sub-daily resolution
	OPEN(11,FILE=pfadn(1:pfadi)//'sediment_production.out', STATUS='old',POSITION='append')
   	WRITE(fmtstr,'(a,i0,a,i0)')'(3i6,',n_sed_class,'(',subasin,'f13.4))'	!generate format string
	
	DO d=1,dayyear
	    DO count=1,nt
			WRITE (11,fmtstr)t, d, count, ((sediment_subbasin_t(d,count,i,j),j=1,n_sed_class),i=1,subasin)
		END DO
	END DO
	CLOSE(11)
  END IF


!  debug output !remove
  OPEN(11,FILE=pfadn(1:pfadi)//'debug.out', STATUS='old',POSITION='append')
		WRITE (11,'(3(i6,a1),f14.3)')(t,char(9), d,char(9), -1,char(9), debug_out(d),&
																				d=1,dayyear)
	CLOSE(11)

!  debug output !remove
  OPEN(11,FILE=pfadn(1:pfadi)//'debug2.out', STATUS='old',POSITION='append')
		DO d=1,dayyear
		WRITE (11,'(3(i6,a1),8(f14.3,a1))')t,char(9), d,char(9), -1,char(9), (debug_out2(d,i),char(9),i=1,8)
		END DO
	CLOSE(11)


	!conrad: output of theta [%] for each tc from theta_tc
	IF (f_tc_theta ) THEN
		OPEN(11,FILE=pfadn(1:pfadi)//'tc_theta.out', STATUS='old',POSITION='append')
		WRITE(fmtstr,'(a,i0,a)')'(i0,a,i0,a,i0,',sv_comb,'(a,f5.3))'	!generate format string
		DO d=1,dayyear
				WRITE (11,fmtstr)(t,char(9),d,char(9),count,(char(9),theta_tc(d,i),i=1,tc_counter_all-1),count=1,nt)	!tab separated output		
		END DO
		CLOSE(11)
	END IF
	!!!!!!!!!end conrad!


	IF (f_tc_surfflow ) THEN	!Till: write TC-wise surface flow
		OPEN(11,FILE=pfadn(1:pfadi)//'tc_surfflow.out', STATUS='old',POSITION='append')
		WRITE(fmtstr,'(a,i0,a)')'(i0,a,i0,a,i0,',sv_comb,'(a,f6.1))'	!generate format string
		
		DO d=1,dayyear
				WRITE (11,fmtstr)(t,char(9),d,char(9),count,(char(9),surfflow_tc(d,count,i),i=1,tc_counter_all-1),count=1,nt)	!tab separated output
		END DO
		
		CLOSE(11)
	END IF


	IF (dosediment .AND. f_tc_sedout .AND. .NOT. do_musle_subbasin) THEN	!Till: write TC-wise surface flow
		OPEN(11,FILE=pfadn(1:pfadi)//'tc_sedout.out', STATUS='old',POSITION='append')
		WRITE(fmtstr,'(a,i0,a)')'(i0,a,i0,a,i0,',sv_comb,'(a,f8.1))'	!generate format string
		
		DO d=1,dayyear
				WRITE (11,fmtstr)(t,char(9),d,char(9),count,(char(9),sedout_tc(d,count,i),i=1,tc_counter_all-1),count=1,nt)	!tab separated output
		END DO
		
		CLOSE(11)
	END IF



if (doacud) CALL lake(3,dummy)
if (dosavestate) call save_model_state !Till: saves model state (soil moisture, ground water, etc.) to file, if specified)



END IF


RETURN

900   FORMAT(i4)
999   FORMAT(3(1X,f5.0))
998   FORMAT(3(1X,f5.3))

contains
	FUNCTION calc_seasonality(subbasin,year,julian_day,node_days,param_node_array)
	!compute seasonality (value of current parameter for current timestep and subbasin) by interpolation between node1 and node2

	use utils_h
	implicit none
		real :: param_node_array(:,:)
		real :: calc_seasonality(size(param_node_array,dim=1))	!return value
		INTEGER, INTENT(IN)                  :: subbasin	!subbasin-id
		INTEGER, INTENT(IN)                  :: year,julian_day	
		integer, INTENT(IN) :: node_days(:)


		
		integer	:: i,j,k
		integer :: d		!distance between start node and current day (in days)
		integer :: d_nodes		!distance between start node and end_node (in days)
		real :: node1_value,node2_value		!parameter values at nodepoints (start and end-point of interpolation)
		integer :: i_node1,i_node2	!indices to nodes
		!integer :: i_current_year	!index to current year in seasonality_array
		

		if (size(node_days)==1) then !no seasonality for this parameter
			calc_seasonality=param_node_array(:,1)	!use single value
			return
		end if

		k=0
		do i=tstart,year-1	!compute sum of days before current simulation year
			k=k+365		!sum up total days since start of start year
			IF (MOD((i),4) == 0)  THEN	!leap year
			  k=k+1
			ENDIF
		end do
		
		i_node1=maxval(whichn(node_days<=julian_day+k,.TRUE.))	!find index to node 1
		i_node2=minval(whichn(node_days> julian_day+k,.TRUE.))	!find index to node 2
		
		if ((i_node1==0) .OR. (i_node2==0)) then  !special cases: simulation day is before or after last specified node 
			i_node1=max(i_node1,i_node2)		!set to first or last node, whichever is non-zero
			calc_seasonality(:)=param_node_array(:,1+mod(i_node1-1,4))		!just extrapolate beyond first or last node-value
		else								!normal case: real interpolation
			d_nodes=node_days(i_node2)-node_days(i_node1)	!distance between nodes
			d      =julian_day+k      -node_days(i_node1)	!distance between node1 and current day

			do i=1,size(calc_seasonality)	!loop through all entities to compute current parameter value
				node1_value=param_node_array(i,1+mod(i_node1-1,4))	!parameter values at nodepoints (start and end-point of interpolation)
				node2_value=param_node_array(i,1+mod(i_node2-1,4))
				calc_seasonality(i)=node1_value+ (node2_value-node1_value)*d/d_nodes		!linear interpolation between nodes
			end do
		end if

		return

	END FUNCTION calc_seasonality


	
	SUBROUTINE write_output(f_flag,file_name,value_array)
	! Output daily values of given array
	IMPLICIT NONE
	LOGICAL, INTENT(IN)                  :: f_flag
	CHARACTER(len=*), INTENT(IN)         :: file_name
	REAL, INTENT(IN)                  :: value_array(:,:)

	
	  IF (f_flag) THEN	!if output file is enabled
		  OPEN(11,FILE=pfadn(1:pfadi)//file_name, STATUS='old',POSITION='append')
		  
		write(fmtstr,'(a,i0,a)') '(i0,a,i0,',subasin,'(a,f14.3))'		!generate format string (daily format)
		DO d=1,dayyear
			write(11,trim(fmtstr))t,char(9),d,(char(9),value_array(d,i),i=1,subasin)
		END DO



		  CLOSE(11)
	  END IF
	END SUBROUTINE write_output

END SUBROUTINE hymo_all


