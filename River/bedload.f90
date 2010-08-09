SUBROUTINE bedload_formulae(i)
 
!! Subroutine calculates the bedload using 5 different formulae for individual river stretches
!! Formula after: Meyer_Peter, Schoklitsch, Smart&Jaeggi, Bagnold, Rickenmann
!! For list of reference, see publication on bedload modelling Müller et al. 2006
!! Implemented: 24.07.06

!!    ~ ~ ~ OUTGOING VARIABLES ~ ~ ~
!!	bedload rate in kg/s as submerged weight

!!    ~ ~ ~ INCOMING VARIABLES ~ ~ ~
!! D50(i)			|m		   	  |D50 parameter of particle size distribution of riverbed
!! r_depth(i)		|m			  |boardful depth
!! r_depth_cur(i)   |m            |current depth of flow on day
!! r_qout(2,i)		|m3/s		  |water discharge
!! velocity(i)		|m/s		  |water flow velocity
!! bottom_width(i)  |m			  |bottom width of riverbed
!! r_sideratio(:)   |none         |change in horizontal distance per unit change in vertical distance on channel side slopes

!! LOCAL DEFINTIONS
!! q_unit           |m2/s         |unit discharge for m river width
!! q_crit			|m2/s		  |critical discharge for initiation of motion
!! shear			|kg/m*s2 or - |shear stress (real or dimensionless)
!! shear_crit		|kg/m*s2 or - |critical shear stress (real or dimensionless)
!!R_Sed				|kg/m3		  |density of sediment
!!R_Water			|kg/m3		  |density of water
!!Fr				|-			  |Froude number

use routing_h
use common_h
use time_h

IMPLICIT NONE

INTEGER i !,j
INTEGER :: R_Sed=2650, R_Water=1000
REAL :: g=9.81
real q_crit, q_unit, shear, shear_crit, Fr, width

! calculation of current width
if (r_depth_cur(i).lt.r_depth(i)) then
	width = bottom_width(i) + r_sideratio(i) * r_depth_cur(i)
else	!no bedload transport on the floodplains, therefore maximum width equivalent to boardful discharge
	width = bottom_width(i) + r_sideratio(i) * r_depth(i)
endif

! calculation of unit water discharge
	q_unit = r_qout(2,i) / width

! Meyer_Peter et al. formula:
! for reference see: Graf 1984, Graf 2001, Martin & Ham 2005, Nicholas 2000
	shear = R_water * g *r_depth_cur(i) * r_slope(i)
	shear_crit = 0.047 * (R_Sed-R_Water) * g * D50(i)
	if (shear_crit.gt.shear) then
		bedload(i,1)  = 0.
	else
		bedload(i,1) = 8 * (shear-shear_crit)**1.5/(g * (R_Sed-R_Water) * R_Water**0.5)
		bedload(i,1)  = bedload(i,1)  * R_Sed
	endif

	bedload(i,1) = bedload(i,1) * (2650-1000)/2650* width		!Unit: kg/sec over entire width

!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Schoklitsch 1950 formula:
!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! for reference see: Graf 1982, Graf 2001, Agostino & Lenzi 1999
	q_crit = 0.26 * ((R_Sed-R_Water)/R_Water)**(5./3.) * D50(i)**1.5/ (r_slope(i)**(7./6.))
	if (q_crit.gt.q_unit) then
		bedload(i,2) = 0.
	else
		bedload(i,2) = 2500 * r_slope(i)**1.5 * (q_unit - q_crit)	
	endif
!NOTE: in der original formula D40 instead of D50 is used!

	bedload(i,2) = bedload(i,2) * (2650-1000)/2650* width		!Unit: kg/sec over entire width

!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Smart and Jaeggi (1983)
! for reference see: Agostino & Lenzi 1999
!!!!!!!!!!!!!!!!!!!!!!!!!!!
	shear = r_depth_cur(i) * r_slope(i) / ((R_Sed/R_Water - 1) * D50(i))
!theorertical critical dimensionless shear stress after Shield: 
    shear_crit= 0.056*(R_Sed-R_Water)*g*D50(i)/(g*(R_Sed-R_Water)*D50(i))
	if (shear_crit.gt.shear) then
		bedload(i,3) = 0.
	else
		bedload(i,3) = 4.2 * q_unit * r_slope(i)**1.6 * (1 - shear_crit/shear) / (R_Sed/R_Water-1)
		bedload(i,3) = bedload(i,3) * R_Sed
	endif

	bedload(i,3) = bedload(i,3) * (2650-1000)/2650* width		!Unit: kg/sec over entire width

!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Bagnold (1956, Yalin 1977) 
! for reference see: Agostino & Lenzi 1999
!!!!!!!!!!!!!!!!!!!!!!!!!!!
	shear = r_depth_cur(i) * r_slope(i) / ((R_Sed/R_Water - 1) * D50(i))
!theorertical critical dimensionless shear stress after Shield: 
    shear_crit= 0.056*(R_Sed-R_Water)*g*D50(i)/(g*(R_Sed-R_Water)*D50(i))
	if (shear_crit.gt.shear) then
		bedload(i,4) = 0.
	else
		bedload(i,4) = 4.25 * shear**0.5 * (shear-shear_crit)
	bedload(i,4) = bedload(i,4) * R_Sed * ((R_sed/R_Water-1) * g * D50(i)**3)**0.5
	endif

	bedload(i,4) = bedload(i,4) * (2650-1000)/2650* width		!Unit: kg/sec over entire width

!!!!!!!!!!!!!!!!!!!
! Rickenmann 1989, 1991
! for reference see Rickenmann and Agostino & Lenzi 1999
!!!!!!!!!!!!!!!!!!!
	shear = r_depth_cur(i) * r_slope(i) / ((R_Sed/R_Water - 1) * D50(i))
!theorertical critical dimensionless shear stress after Shield: 
    shear_crit= 0.056*(R_Sed-R_Water)*g*D50(i)/(g*(R_Sed-R_Water)*D50(i))
	Fr = velocity(i) / (g * r_depth_cur(i))**0.5
	if (shear_crit.gt.shear) then
		bedload(i,5)= 0.
	else
		bedload(i,5) = 3.1*shear**0.5*(shear-shear_crit) * Fr**1.1 * (R_Sed/R_Water-1)**(-0.5)
		bedload(i,5) = 	bedload(i,5) * R_Sed * ((R_Sed/R_Water-1) * g * D50(i)**3)**0.5
	endif

	bedload(i,5) = bedload(i,5) * (2650-1000)/2650* width		!Unit: kg/sec over entire width


RETURN
END SUBROUTINE bedload_formulae




