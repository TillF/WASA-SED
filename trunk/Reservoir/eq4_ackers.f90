SUBROUTINE eq4_ackers(upstream)
 
! Code converted using TO_F90 by Alan Miller
! Date: 2005-08-23  Time: 12:57:31
 
use common_h
use time_h
use reservoir_h

IMPLICIT NONE

INTEGER, INTENT(IN OUT)                  :: upstream
INTEGER :: i,j,g,h,id,ih,irout,imun,dummy1
real:: dummy2


real :: manning_sed,shear_bed,shear_sec,d50,accum1,accum2
real :: hidp(n_sed_class),expp(n_sed_class),hidexp(n_sed_class)
real :: wat_dens,sed_dens,shield,crshear(n_sed_class)
real :: visc,perc_bed(n_sed_class),perc_susp(n_sed_class)
real :: transpfac_bed(n_sed_class),transpfac_susp(n_sed_class)
real :: D_gs(n_sed_class),diam_u,shear_vel

!Ge to include the DAILY mean temperature of the reservoir (celsius degree)
!Ge tempres=20 C (temporarily)
REAL :: tempres(step,upstream),n_coef(n_sed_class),m_coef(n_sed_class),k_coef(n_sed_class),Fgr_cr(n_sed_class),Fgr(n_sed_class),rel_dens

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
  
!Ge mean diameter has to be calculated from the grain size distribution
  dummy2=d50_actlay(j,upstream)



! TOTAL LOAD FORMULA
    
! Gravitacional acceleration (m/s2)
! 9.807
  
! density of water and density of natural sediments (kg/m3)
  wat_dens=1.*1000.
  sed_dens=2.65*1000.

! relative density of sediment
  rel_dens=(sed_dens-wat_dens)/wat_dens
 
! bed shear stress (N/m2)
  shear_bed=9.807*wat_dens*depth(j)* energslope(j)
  
!Ge to include the DAILY mean temperature of the reservoir (celsius degree)
!Ge tempres=20 C (temporarily)
  tempres(step,upstream)=20
  
! kinematic viscosity (m2/s)
  visc=(1.14-0.031*(tempres(step,upstream)-15.)+  &
      0.00068*((tempres(step,upstream)-15.)**2.))*1.e-6
!write(*,*)j,visc

! loop to calculate the sediment carrying capacity for bed load
  DO g=1,n_sed_class
   IF (diam(g)*1000. > 0.040) THEN
    if(dummy2==0.) then
	  d50=diam(g)
	else
	  d50=dummy2
	endif
!write(*,*)j,g,d50

! dimensionless grain size (-)
	D_gs(g)=diam(g)*((9.807*rel_dens/(visc**2))**(1./3.))
!write(*,*)j,g,D_gs(g),diam(g),rel_dens,visc

! coefficients Ackers and White Formula
    if (D_gs(g)>1. .and. D_gs(g)<60) then
	  n_coef(g)=1.-(.56*log10(D_gs(g)))
	  m_coef(g)=(9.66/D_gs(g))+1.34
	  k_coef(g)=10**(-3.53+2.86*log10(D_gs(g))-(log10(D_gs(g))**2))
	  Fgr_cr(g)=.14+(.23/sqrt(D_gs(g)))
	else if (D_gs(g)>=60.) then
	  n_coef(g)=0.
	  m_coef(g)=1.5
	  k_coef(g)=.025
	  Fgr_cr(g)=.17
	endif	
!if(g==9)write(*,*)j,g,D_gs(g),n_coef(g),m_coef(g),k_coef(g),Fgr_cr(g)
 

! shear velocity (m/s)
    shear_vel=sqrt(9.807*depth(j)*energslope(j))

! sediment mobility number (-)
    Fgr(g)=(shear_vel**n_coef(g)/(sqrt(9.807*diam(g)*rel_dens)))* &
	     (meanvel(j)/(sqrt(32.)*log10(10*depth(j)/diam(g))))**(1-n_coef(g))

!if(g==10)write(*,*)j,g,(shear_vel**n_coef(g)/(sqrt(9.807*diam(g)*rel_dens))),(meanvel(j)/(sqrt(32.)*log10(10*depth(j)/diam(g))))**(1-n_coef(g)),Fgr(g)
!if(g==10)write(*,*)j,g,depth(j),energslope(j),shear_vel,diam(g)*1000,meanvel(j)

   
! shields parameter (-)
    shield=shear_bed/(wat_dens*rel_dens*9.807*d50)

! grain sizes that requires no correction (-)
    if (shield<=.04) then
	  diam_u=1.08*d50
	else if (shield>.04 .and. shield<=.046) then
	  diam_u=(-21.6*shield+1.944)*d50
	else if (shield>.46 .and. shield<=.097) then
	  diam_u=(-9.73*shield+1.40)*d50
	else if (shield>.097) then
	  diam_u=.456*d50
	endif

!if(g==10)write(*,*)j,g,shear_bed,wat_dens,rel_dens,d50,shield,diam_u


! hiding/exposure factor (-)
    if (diam_u/=0.) then
      if (diam(g)/diam_u<=.075) then
        hidexp(g)=2.5
      else if (diam(g)/diam_u>.075 .and. diam(g)/diam_u<3.7) then
	    hidexp(g)=1/((.53*log10(diam(g)/diam_u))+1)
	  else if (diam(g)/diam_u>=3.7) then
        hidexp(g)=.769
	  endif
	endif 

!if(g==10)write(*,*)j,g,hidexp(g)
 
    
! loop to calculate the fractional carrying capacity for total load (m3/day)
    if (Fgr(g)/(Fgr_cr(g)*hidexp(g)) <= 1.) then
	  fr_capacity(g,j)=0.
	else
      fr_capacity(g,j)=k_coef(g)*meanvel(j)*diam(g)*((meanvel(j)/shear_vel)**n_coef(g))*  &
        (((Fgr(g)/(Fgr_cr(g)*hidexp(g)))-1)**m_coef(g))
    endif

!if(g==10)write(*,*)j,g,k_coef(g)*meanvel(j)*diam(g),((meanvel(j)/shear_vel)**n_coef(g)),(((Fgr(g)/(Fgr_cr(g)*hidexp(g)))-1)**m_coef(g))

!if(g==10)write(*,*)j,g,fr_capacity(g,j),topwidth(j)
    
! total suspended sediment transport (m3/day)
    fr_capacity(g,j)=min(1000.,fr_capacity(g,j))
    fr_capacity(g,j)=86400*fr_capacity(g,j)*topwidth(j)

   ELSE
    fr_capacity(g,j)=1000.
    fr_capacity(g,j)=86400*fr_capacity(g,j)*topwidth(j)
   ENDIF

   setvel(g)=SQRT((13.95*visc/diam(g))**2.+1.09*  &
        ((sed_dens-wat_dens)/wat_dens)*9.807*diam(g)) -13.95*visc/diam(g)

  END DO

  else
    DO g=1,n_sed_class
      fr_capacity(g,j)=0.
	enddo
  endif

!write(*,*)d,j,(fr_capacity(g,j),g=1,n_sed_class)
!write(*,*)d,j,fr_capacity(9,j)
  DO g=1,n_sed_class
!if (j==10)write(*,'(3I4,3F20.5)')d,j,g,fr_capacity(g,j)
  END DO
!if (j==10)stop

!if (d==54)read(*,*)
END DO

RETURN
END SUBROUTINE eq4_ackers
