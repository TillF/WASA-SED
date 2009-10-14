SUBROUTINE readgen(path2do_dat)

!Till: computationally irrelevant: made "command-strings" in outfile.dat case-insensitive
!2009-10-14

!Till: computationally irrelevant: added program version information to parameter.out
!2009-06-17

!Till: added output for River_Sediment_Storage.out
!sediment-related output-files are all disabled, if dosediment=FALSE
!2008-11-13

!Till: if location of do.dat as argument is given, all path are interpreted relative to the location of do.dat
! 2008-05-30

!Till: optionally get location of do.dat as argument
! 2008-04-24

!Till: read in parameters for time-variate kfcorr
!2007-10-29

! 2007-10-18 Till
! increased length of pfadp variable to 160

! 2007-06-04, Till
! added flag f_tc_theta

! 2007-04-28, Till
! added flag f_tc_theta

! 2007-01-10, Till
! renamed f_deep_gw to f_deep_gw_recharge, added f_deep_gw_discharge
 

! 2005-10-24, Till
! optional reading of outfiles.dat to configure, which output files are created

! 2005-08-24, Eva
! read additional information from do.dat
 
! 2005-08-09, Till
! read optional dosediment from do.dat
! read maximum dimensions of spatial units from maxdim.dat, if available

! Code converted using TO_F90 by Alan Miller
! Date: 2005-06-30  Time: 13:47:18

use climo_h
use common_h
use hymo_h
use params_h
use time_h
use lake_h
use routing_h
use reservoir_h
use erosion_h
use utils_h

IMPLICIT NONE
CHARACTER (LEN=160) :: path2do_dat		!Till: path to central control file (do.dat)

INTEGER :: i,istate !,imun,imicro,imeso
CHARACTER (LEN=150) :: custompath
CHARACTER (30) :: dummy

if (trim(path2do_dat)=='') then
	path2do_dat='./Input/do.dat'			!Till: use default, if no command line argument was specified
	custompath=''
else
	write(*,*)'reading runtime parameters from ',path2do_dat
	i=sizeof(trim(path2do_dat))
	do while (i>0)
		if ((path2do_dat(i:i)=='/') .OR. (path2do_dat(i:i)=='\')) then	!find last slash in path
			exit
		end if 
		i=i-1
	end do
	custompath=path2do_dat(1:i)	!extract path to do.dat (without filename)
end if

! read data for simulation period from file
!OPEN(11,FILE='./Input/do.dat' ,STATUS='old')
OPEN(11,FILE=trim(path2do_dat) ,STATUS='old')
READ(11,*)
READ(11,'(a)') pfadp
READ(11,'(a)') pfadn
READ(11,*) tstart
READ(11,*) tstop
READ(11,*) mstart
READ(11,*) mstop
READ(11,*) subasin  !total no. of sub-basins
READ(11,*) sv_comb  !total no. of sub-basin / SO / TC combinations
READ(11,*) nsoter  !total no. of SOTER units in study area
READ(11,*) nterrain  !total no. of terrain components in study area
READ(11,*) nsoil  !total no. of soil components in study area
READ(11,*) nveg   !total no. of vegetation units in study area
READ(11,*) doreservoir !do large reservoir calculation (George's modules)
READ(11,*) doacud
READ(11,*) dolattc
READ(11,*) doalllattc
READ(11,*) dolatsc
READ(11,*) dolatscsub
READ(11,*) dotrans

READ(11,*) dohour
READ(11,*) scenario
READ(11,*) krig
!READ(11,*) kfkorr
READ(11,'(A)') dummy
	READ(dummy,*,IOSTAT=i)kfkorr,kfkorr_a,kfkorr_b
	IF (i/=0 .OR. kfkorr_a==0) THEN	!no parameters for variable kfcorr specified, use constant kfcorr (as in old version)
		READ(dummy,*)kfkorr
		kfkorr_a=0
		kfkorr_b=1
	END IF

READ(11,*) intcf
READ(11,*) dointc
READ(11,*) doscale
READ(11,*) domuncell
READ(11,*) sensfactor
READ(11,*) dt
READ(11,*) dosediment
READ(11,*) n_sed_class
READ(11,*) hill_transport
READ(11,*) river_transport
READ(11,*) reservoir_transport

CLOSE(11)

if (trim(custompath)/='') then		!if a custom path was specified, all paths are relative to this one
	pfadp=trim(custompath)//pfadp
	pfadn=trim(custompath)//pfadn
end if


pfadi=sizeof(trim(pfadn))	! determine length of path name for output files
pfadj=sizeof(trim(pfadp))	! determine length of basic model path name
ncaset=sizeof(trim(caset))	! determine length of case study path name


!Till: read maximum dimensions of arrays
OPEN(11,FILE=pfadp(1:pfadj)// 'maxdim.dat',IOSTAT=i,STATUS='old')	
IF (i==0) THEN
	READ(11,*)
	READ(11,*)maxsoter
	READ(11,*)maxterrain
	READ(11,*)maxsoil
	READ(11,*)maxhori
	READ(11,*)ntrans
	CLOSE(11)
else
	write(*,*)pfadp(1:pfadj)// 'maxdim.dat could not be opened. Using default dimensions specified in the source code.'
	maxsoter=7
	maxterrain=3
	maxsoil=28
	maxhori=8
	ntrans=2
END IF


!Till: read desired output files
f_daily_actetranspiration=.FALSE.	!disable all output files
f_daily_potetranspiration=.FALSE.
f_daily_qhorton=.FALSE.
f_daily_qin_m3s=.FALSE.
f_daily_qout_m3s=.FALSE.
f_daily_rain=.FALSE.
f_daily_runoff=.FALSE.
f_daily_sediment_production=.FALSE.
f_daily_subsurface_runoff=.FALSE.
f_daily_theta=.FALSE.
f_daily_total_overlandflow=.FALSE.
f_daily_water_subbasin=.FALSE.
f_routing_response=.FALSE.
f_sediment_production=.FALSE.
f_water_subbasin=.FALSE.
f_deep_gw_recharge=.FALSE.
f_deep_gw_discharge=.FALSE.
f_tc_theta=.FALSE.
f_daily_gw_loss=.FALSE.
f_river_degradation=.FALSE.
f_river_deposition=.FALSE.
f_river_flow=.FALSE. 
f_river_flow_dailyaverage=.FALSE.
f_river_flowdepth=.FALSE.
f_river_sediment_concentration =.FALSE.
f_river_sediment_total=.FALSE.
f_river_sediment_total_dailyaverage =.FALSE.
f_river_storage=.FALSE.
f_river_sediment_storage=.FALSE.
f_river_velocity=.FALSE.
f_river_bedload=.FALSE.
f_tc_surfflow=.FALSE.
f_tc_sedout=.FALSE.

f_actetranspiration=.FALSE.
f_qhorton=.FALSE.
f_subsurface_runoff=.FALSE.
f_total_overlandflow=.FALSE.
f_gw_discharge=.FALSE.
f_potetranspiration=.FALSE.
f_gw_loss=.FALSE.
f_gw_recharge=.FALSE.

OPEN(11,FILE=pfadp(1:pfadj)// 'outfiles.dat',IOSTAT=istate,STATUS='old')	
IF (istate==0) THEN
	READ(11,*,IOSTAT=istate)dummy  
	READ(11,*,IOSTAT=istate)dummy  
	READ(11,*,IOSTAT=istate)dummy  

	DO WHILE (istate==0)
		SELECT CASE (trim(locase(dummy)))		!enable/disable file output of desired results
			CASE ('daily_actetranspiration')
			  f_daily_actetranspiration=.TRUE.
			CASE ('daily_potetranspiration')
			  f_daily_potetranspiration=.TRUE.
			CASE ('daily_qhorton')
			  f_daily_qhorton=.TRUE.
			CASE ('daily_qin_m3s')
			  f_daily_qin_m3s=.TRUE.
			CASE ('daily_qout_m3s')
			  f_daily_qout_m3s=.TRUE.
			CASE ('daily_rain')
			  f_daily_rain=.TRUE.
			CASE ('daily_runoff')
			  f_daily_runoff=.TRUE.
			CASE ('daily_sediment_production')
			  f_daily_sediment_production = dosediment
			CASE ('daily_subsurface_runoff')
			  f_daily_subsurface_runoff=.TRUE.
			CASE ('daily_theta')
			  f_daily_theta=.TRUE.
			CASE ('daily_total_overlandflow')
			  f_daily_total_overlandflow=.TRUE.
			CASE ('daily_water_subbasin')
			  f_daily_water_subbasin=.TRUE.
			CASE ('routing_response')
			  f_routing_response=.TRUE.
			CASE ('sediment_production')
			  f_sediment_production=dosediment
			CASE ('water_subbasin')
			  f_water_subbasin=.TRUE.
			CASE ('deep_gw_recharge')
			  f_deep_gw_recharge=.TRUE.
			CASE ('deep_gw_discharge')
			  f_deep_gw_discharge=.TRUE.
			CASE ('f_daily_gw_loss')
			  f_daily_gw_loss=.TRUE.
			CASE ('tc_theta')
			  f_tc_theta=.TRUE.
			CASE ('river_degradation')
			  f_river_degradation=dosediment
			CASE ('river_deposition')
			  f_river_deposition=dosediment
			CASE ('river_flow')
				f_river_flow=.TRUE. 
			CASE ('river_flow_dailyaverage')
				f_river_flow_dailyaverage=.TRUE.
			CASE ('river_flowdepth')
				f_river_flowdepth=.TRUE.
			CASE ('river_sediment_concentration')
				f_river_sediment_concentration =dosediment
			CASE ('river_sediment_total')
				f_river_sediment_total=.TRUE.
			CASE ('river_sediment_total_dailyaverage')
				f_river_sediment_total_dailyaverage =dosediment
			CASE ('river_storage')
				f_river_storage=.TRUE.
			CASE ('river_sediment_storage')
				f_river_sediment_storage=dosediment
			CASE ('river_velocity')
				f_river_velocity=.TRUE.
			CASE ('river_bedload')
				f_river_bedload=dosediment
			CASE ('tc_surfflow')
				f_tc_surfflow=.TRUE.
			CASE ('tc_sedout')
				f_tc_sedout=.TRUE.

			CASE ('actetranspiration')
				f_actetranspiration=.TRUE.
			CASE ('qhorton')
				f_qhorton=.TRUE.
			CASE ('subsurface_runoff')
				f_subsurface_runoff=.TRUE.
			CASE ('total_overlandflow')
				f_total_overlandflow=.TRUE.
			CASE ('gw_discharge')
				f_gw_discharge=.TRUE.
			CASE ('potetranspiration')
				f_potetranspiration=.TRUE.
			CASE ('gw_loss')
				f_gw_loss=.TRUE.
			CASE ('gw_recharge')
				f_gw_recharge=.TRUE.


 	END SELECT
		READ(11,*,IOSTAT=istate)dummy	!try to read next line
	END DO
	CLOSE(11)
	
	f_qhorton= f_qhorton .OR. f_daily_qhorton

ELSE
	WRITE(*,*)pfadp(1:pfadj)// 'outfiles.dat could not be opened. Using default output files.'
	!these are the default output files
	f_daily_runoff=.TRUE.
	f_daily_sediment_production=dosediment
	f_daily_water_subbasin=.TRUE.
	f_water_subbasin=.TRUE.
	f_river_flow=.TRUE.
	f_river_sediment_total=dosediment
	f_river_sediment_concentration=dosediment
END IF
!end insert Till


! save settings of this run to output directory
OPEN(11,FILE=pfadn(1:pfadi)//'parameter.out', STATUS='unknown',IOSTAT=istate)	
IF (istate==0) THEN
	WRITE(11,*)
	WRITE(11,'(a)') pfadp
	WRITE(11,'(a)') pfadn
	WRITE(11,*) 'start year of simulation: ', tstart
	WRITE(11,*) 'end year of simulation: ', tstop
	WRITE(11,*) 'start month of simulation: ',mstart
	WRITE(11,*) 'end month of simulation: ',mstop
	WRITE(11,*) 'no. of sub-basin: ',subasin
	WRITE(11,*) 'no. of combinations: ',sv_comb
	WRITE(11,*) 'total no. of SOTER units: ',nsoter
	WRITE(11,*) 'total no. of terrain components: ',nterrain
	WRITE(11,*) 'total no. of soil components: ',nsoil
	WRITE(11,*) 'total no. of vegetation units: ',nveg
	WRITE(11,*) 'do reservoir calculations: ',doacud
	WRITE(11,*) dolattc
	WRITE(11,*) doalllattc
	WRITE(11,*) dolatsc
	WRITE(11,*) dolatscsub
	WRITE(11,*) dotrans
	WRITE(11,*) dohour
	WRITE(11,*) scenario
	WRITE(11,*) krig
	WRITE(11,*) kfkorr,kfkorr_a,kfkorr_b
	WRITE(11,*) intcf
	WRITE(11,*) dointc
	WRITE(11,*) doscale
	WRITE(11,*) domuncell
	WRITE(11,*) sensfactor
	write(11,*) dosediment
	write(11,*) n_sed_class
	write(11,*) hill_transport
	write(11,*) river_transport
	write(11,*) reservoir_transport
	if (dosediment) then
		WRITE(11,*) 
		WRITE(11,*) 'spatial scale for application of erosion eq : ',do_musle_subbasin
		WRITE(11,*) 'erosion eq : ',erosion_equation
	end if
	WRITE(11,*) 'WASA model, ',trim(rev_string1),'; ',trim(rev_string2)
	CLOSE(11)
else
	write(*,*)'Error: Output file ',pfadn(1:pfadi)//'parameter.out',' could not be created, aborting.'
	stop
END IF


nt = int(24/dt)	!Till: number of simulation steps per day		

INCLUDE '../General/allocat_general.var'	
INCLUDE '../Hillslope/allocat_hymo.var'
INCLUDE '../Hillslope/allocat_erosion.var'
INCLUDE '../River/allocat_routing.var'
nt = int(24/dt)	!Till: number of simulation steps per day		
if (doreservoir.or.doacud) then
	INCLUDE '../Reservoir/allocat_reservoir_lake.var'
end if


RETURN

END SUBROUTINE readgen
