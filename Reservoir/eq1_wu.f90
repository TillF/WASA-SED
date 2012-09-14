SUBROUTINE eq1_wu(upstream)
!Till: computationally irrelevant: outcommented unused vars
!2012-09-14 

!Till: computationally irrelevant: minor changes to improve compiler compatibility
!2011-04-29

! Code converted using TO_F90 by Alan Miller
! Date: 2005-08-23  Time: 12:57:31
 
use common_h
use time_h
use reservoir_h

IMPLICIT NONE

INTEGER, INTENT(IN OUT)                  :: upstream
INTEGER :: j,g,h !i,id,ih,irout,imun,dummy1,,factor_sed
real:: dummy2


real :: manning_sed,shear_bed,shear_sec,d50 !,accum1,accum2
real :: hidp(n_sed_class),expp(n_sed_class),hidexp(n_sed_class)
real :: wat_dens,sed_dens,shield,crshear(n_sed_class)
real :: visc !,perc_bed(n_sed_class),perc_susp(n_sed_class)
real :: transpfac_bed(n_sed_class),transpfac_susp(n_sed_class)


!Ge to include the DAILY mean temperature of the reservoir (celsius degree)
!Ge tempres=20 C (temporarily)
REAL :: tempres(step,upstream)

real :: discharge(200),depth(200),hydrad(200),energslope(200)
real :: manning(200),topwidth(200),meanvel(200)


! -----------------------------------------------------------------------

do j=1,nbrsec(upstream)
 if (j /= nbrsec(upstream)) then
    discharge(j)=(discharge_sec(j,upstream)+discharge_sec(j+1,upstream))/2.
    depth(j)=(depth_sec(j,upstream)+depth_sec(j+1,upstream))/2.
    hydrad(j)=(hydrad_sec(j,upstream)+hydrad_sec(j+1,upstream))/2.
    energslope(j)=(energslope_sec(j,upstream)+energslope_sec(j+1,upstream))/2.
    manning(j)=(manning_sec(j,upstream)+manning_sec(j+1,upstream))/2.
    topwidth(j)=(topwidth_sec(j,upstream)+topwidth_sec(j+1,upstream))/2.
    meanvel(j)=(meanvel_sec(j,upstream)+meanvel_sec(j+1,upstream))/2.
 else
    discharge(j)=discharge_sec(j,upstream)
    depth(j)=depth_sec(j,upstream)
    hydrad(j)=hydrad_sec(j,upstream)
    energslope(j)=energslope_sec(j,upstream)
    manning(j)=manning_sec(j,upstream)
    topwidth(j)=topwidth_sec(j,upstream)
    meanvel(j)=meanvel_sec(j,upstream)
 endif
enddo


! loop to calculate the sediment carrying capacity
DO j=1,nbrsec(upstream)

 if (discharge(j) /= 0.) then
  
!Ge mean diameter is calculated using the grain size distribution
  dummy2=d50_actlay(j,upstream)


! 1) BED LOAD FORMULA
  
  
  
! Gravitacional acceleration (m/s2)
! 9.807
  
! density of water and density of natural sediments (kg/m3)
  wat_dens=1.*1000.
  sed_dens=2.65*1000.
  
! bed shear stress (N/m2)
  shear_bed=9.807*wat_dens*depth(j)* energslope(j)

! shear stress on the entire cross section (N/m2)
  shear_sec=9.807*wat_dens*hydrad(j)* energslope(j)
  
! hiding probability (-)
! exposure probability (-)
  DO g=1,n_sed_class
    hidp(g)=0.
	expp(g)=0.
    DO h=1,n_sed_class
      IF (frac_actlay(g,j,upstream) /= 0.) THEN
        hidp(g)=hidp(g)+frac_actlay(h,j,upstream)* (diam(h)/(diam(g)+diam(h)))
        expp(g)=expp(g)+frac_actlay(h,j,upstream)* (diam(g)/(diam(g)+diam(h)))
      ELSE
        hidp(g)=1.
        expp(g)=1.
      END IF
    END DO
  END DO
  
! critical Shields parameter (-)
! 0.03 according to Wu et al.
  shield=0.03
  
! loop to calculate the sediment carrying capacity for bed load
  DO g=1,n_sed_class
    if(dummy2==0.) then
	  d50=diam(g)
	else
	  d50=dummy2
	endif

! Manning's roughness related to grains (s/(m**1/3))
    manning_sed=(d50**(1./6.))/20.
    
! hiding/exposure factor (-)
    hidexp(g)=(expp(g)/hidp(g))**(-0.6)
    
! critical shear stress (N/m2)
    crshear(g)=(sed_dens-wat_dens)*9.807*diam(g)*shield*hidexp(g)
    
! dimensionless transport parameter for fractional bed load yields (-)
    if ((((manning_sed/manning(j))**(3./2.))*shear_bed/crshear(g))-1. <= 0.) then
	  transpfac_bed(g)=0.
	else
      transpfac_bed(g)=0.0053*((((manning_sed/manning(j))  &
        **(3./2.))*shear_bed/crshear(g))-1.)**2.2
    endif
!    transpfac_bed(g)=MAX(0.,transpfac_bed(g))

! fractional bed sediment transport (m2/s)
    bed_frtransp(g,j)=transpfac_bed(g)*  &
        SQRT((sed_dens-wat_dens)/wat_dens*9.807*(diam(g)**3.))
    
! fractional bed sediment transport (m3/day)
    bed_frtransp(g,j)=86400.*bed_frtransp(g,j)*topwidth(j)
    
!write(*,'(3I4,6F10.6)')d,j,g,d50,manning_sed,hidexp(g),shear_bed,crshear(g),transpfac_bed(g)
  END DO

  
! 2) SUSPENDED LOAD FORMULA
  
!Ge to include the DAILY mean temperature of the reservoir (celsius degree)
!Ge tempres=20 C (temporarily)
  tempres(step,upstream)=20.
  
! kinematic viscosity (m2/s)
  visc=(1.14-0.031*(tempres(step,upstream)-15.)+  &
      0.00068*((tempres(step,upstream)-15.)**2.))*1.e-6
  
! loop to calculate the sediment carrying capacity for suspended load
  DO g=1,n_sed_class
    
! settling velocity (m/s)
    setvel(g)=SQRT((13.95*visc/diam(g))**2.+1.09*  &
        ((sed_dens-wat_dens)/wat_dens)*9.807*diam(g))-13.95*visc/diam(g)
    
! the fractional transport parameter for suspended sediment transport (-)
    if ((shear_sec/crshear(g))-1. <= 0.) then
	  transpfac_susp(g)=0.
	else
      transpfac_susp(g)=0.0000262*(((shear_sec/crshear(g))-1.)*  &
        meanvel(j)/setvel(g))**1.74
    endif
    
!    transpfac_susp(g)=MAX(0.,transpfac_susp(g))
    
    
! total suspended sediment transport (m2/s)
    susp_frtransp(g,j)=transpfac_susp(g)*  &
        SQRT((sed_dens-wat_dens)/wat_dens*9.807*(diam(g)**3.))
    
! total suspended sediment transport (m3/day)
!    susp_frtransp(g,j)=min(susp_frtransp(g,j),1000.)
    susp_frtransp(g,j)=86400.*susp_frtransp(g,j)*topwidth(j)
    
  END DO

! loop to calculate the fractional carrying capacity for total load (m3/day)
  DO g=1,n_sed_class
    fr_capacity(g,j)=susp_frtransp(g,j)+bed_frtransp(g,j)
  END DO

 else
  DO g=1,n_sed_class
    fr_capacity(g,j)=0.
  enddo
 endif

 DO g=1,n_sed_class
!write(*,'(3I4,3F20.5)')d,j,g,fr_capacity(g,j),susp_frtransp(g,j),bed_frtransp(g,j)
 END DO

!if (d==54)read(*,*)
END DO

RETURN
END SUBROUTINE eq1_wu
