SUBROUTINE petcalc
!read input file containing extra terrestrial radiation
!compute potential evaporation from water surface (applied for each subbasin and day)

! Till: computationally irrelevant: minor changes to improve compiler compatibility
! 2011-04-29

! Till: adapted day-night-calculation scheme that is also used in etp_max to have consistent results
! 2008-02-07

! Till: prevent PET from getting negative
! 2006-02-20

! Code converted using TO_F90 by Alan Miller
! Date: 2005-06-30  Time: 13:47:16
 
use climo_h
use common_h
use params_h
use time_h

! Ausgehend davon, das Temperatur, Globalstrahlung und Luftfeuchte als taegl.,
! Sub-basin abhaengige  Werte gegeben sind, also temp(id,isubbas), rad(id,isubbas),
! rhum(id,isubbas)
! Wind ist erstmal fuer alle Sub-basins und Jahre gleichgesetzt

IMPLICIT NONE

INTEGER :: isubbas,im,id,idmonstart


REAL :: rex(12)	! radiation (W/m^2)
REAL :: alpha	! albedo
REAL :: ra  , rs	! aerodynamic / stomata  resistance  (s/m)
REAL :: z0, z		! roughness length of canopy, height of measurement of wind speed  (m)
REAL :: tempd,tempn	! day and night temperature [°C]
REAL :: rshort	! short wave radiation [W/m²]


z=6.0		! free water surface
!z=10.0		! FAO-Grass reference evaporation

z0=0.00137	! free water surface
!z0=0.12	! FAO-Grass reference evaporation

alpha=0.05	! free water surface
!alpha=0.23 ! FAO-Grass reference evaporation
      

! Reads in radiation as monthly means (in contrast to rad(i,nmun) which are daily means)
OPEN(11,FILE=pfadp(1:pfadj)//'Time_series/extraterrestrial_radiation.dat',STATUS='old')	!ii: load this file only once
READ(11,*)
DO im=1, 12
  READ(11,*)rex(im)
END DO
CLOSE(11)

im=1
IF (t == tstart) THEN
  im=mstart
END IF

idmonstart=1
DO id=1,dayyear
  IF (id-idmonstart >= daymon(im)) THEN
    idmonstart = idmonstart + daymon(im)
    im         = im+1
  END IF
  radex(id)=rex(im)
  DO isubbas=1, subasin
    ra = 4.72*(LOG(z/z0))**2/(1+0.536*wind(id,isubbas))	! free water surface
	!ra = 208./wind(id,isubbas)	! FAO-Grass reference evaporation
    rs=0.												! free water surface
	!rs=69.0	! FAO-Grass reference evaporation
    
	
	if (domean) then
		pet(id,isubbas)=et_pen_mon(ra,rs,temp(id,isubbas),alpha,rhum(id,isubbas),rad(id,isubbas),radex(id))	!Till: do computation for daily mean values
	else
		!Till: do computation for day and night separately
		tempd=temp(id,isubbas)+daily_delta_temp	!Till: temperature during day
		rshort =rad(id,isubbas)*24./hours_of_daylight	!Till: radiation concentrates during daylight hours
		pet(id,isubbas)=et_pen_mon(ra,rs,tempd,alpha,rhum(id,isubbas),rshort,radex(id))*hours_of_daylight/24.
		if (donight) then		!Till: also include nighttime evaporation
			tempn=temp(d,isubbas)-daily_delta_temp	!Till: temperature during night
			pet(id,isubbas)=pet(id,isubbas) +&
						et_pen_mon(ra,rs,tempd,alpha,rhum(id,isubbas),0.,radex(id))*(1.-hours_of_daylight/24.)
		end if
	end if
	
  END DO

  
END DO


contains 
real FUNCTION et_pen_mon(ra,rs,tempr,alpha, rhum,rad, radex) 
!compute daily evapotranspiration [mm] according to Penman-Monteith
implicit none
REAL, intent(in) :: ra  , rs	! aerodynamic / stomata  resistance  (s/m)
REAL, intent(in) :: tempr		! temperature [°C]
REAL, intent(in) :: alpha		! albedo
REAL, intent(in) :: rhum		! relative humidity [%]
REAL, intent(in) :: rad			! daily mean shortwave radiation  [W/m²]
REAL, intent(in) :: radex		! incoming radiation above atmosphere [W/m²]


REAL :: e,  es, s	! vapour pressure (hPa), saturation vapour pressure (hPa) and slope of  the saturation vapour pressure curve  (hPa/K)
REAL :: nn,f	! cloud factor, cover fraction
REAL :: rnetto	! radiation (W/m^2)

REAL,parameter :: cp= 1013.	! specific (J/kg/K)
REAL,parameter :: rho= 1.18		!air density (kg/m^3)
REAL,parameter :: gamma=0.67	! psychrometric constant  (hPa/K)

    es = 6.11*EXP(17.62*tempr/(243.12+tempr))
    e  = rhum*es/100.
    s  = es*4284./((243.12+tempr)**2)
    nn = rad/radex/0.55-0.18/0.55
    nn=MAX(0.0,nn)
    nn=MIN(1.0,nn)
    f  = 0.1+0.9*nn
! Rnetto=Rshort+Rlongwave
    rnetto = -f*(0.52-0.065*SQRT(e))* 5.67E-8*(tempr+273.2)**4  &	! longwave net radiation (Brunt, 1932 in Dyck&Peschke, 1995)
        +(1.-alpha)*rad					! shortwave net radiation
    
    et_pen_mon  = (0.99*rnetto*s+rho*cp/ra*(es-e))*0.0353/ (s+gamma*(1.+rs/ra))
	
!	pmp=(s*anetto+((roh*heatc*             (es-emean)-s*rap*anettos)/ (raa+rap)))/ (s+gamma*(1.+rsp/(raa+rap)))*transn
	
	et_pen_mon = MAX(0.,et_pen_mon)	!Till: prevent PET from getting negative

END FUNCTION et_pen_mon



END SUBROUTINE petcalc


