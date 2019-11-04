SUBROUTINE eq2_ashida(upstream)
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
INTEGER :: j,g,h !i,id,ih,irout,imun,dummy1,factor_sed
real:: dummy2


real :: crshear_D50,d50 !manning_sed,shear_bed,shear_sec,,accum1,accum2
real :: spec_gravit,shear_vel,effshear_vel,karman,crshearvel_D50,C1
real :: crshear_vel(n_sed_class),crshear(n_sed_class),effshear(n_sed_class),shear(n_sed_class)
real :: wat_dens,sed_dens
real :: visc !,perc_bed(n_sed_class),perc_susp(n_sed_class)
real :: transpfac_bed(n_sed_class) !,transpfac_susp(n_sed_class)
real :: seta0,dummy6,dummy7,dummy3,dummy4,dummy,f_seta0,f_seta,Ca,a !,dummy5


!Ge to include the DAILY mean temperature of the reservoir (celsius degree)
!Ge tempres=20 C (temporarily)
REAL :: tempres(step,upstream)

real :: discharge(200),depth(200),hydrad(200),energslope(200)
real :: manning(200),topwidth(200),meanvel(200)


! -----------------------------------------------------------------------

do j=1,nbrsec(upstream)
  if (j /= nbrsec(upstream)) then
    discharge(j)=(discharge_sec(j,res_index(upstream))+discharge_sec(j+1,res_index(upstream)))/2.
    depth(j)=(depth_sec(j,res_index(upstream))+depth_sec(j+1,res_index(upstream)))/2.
    hydrad(j)=(hydrad_sec(j,res_index(upstream))+hydrad_sec(j+1,res_index(upstream)))/2.
    energslope(j)=(energslope_sec(j,res_index(upstream))+energslope_sec(j+1,upstream))/2.
    manning(j)=(manning_sec(j,res_index(upstream))+manning_sec(j+1,res_index(upstream)))/2.
    topwidth(j)=(topwidth_sec(j,res_index(upstream))+topwidth_sec(j+1,res_index(upstream)))/2.
    meanvel(j)=(meanvel_sec(j,res_index(upstream))+meanvel_sec(j+1,res_index(upstream)))/2.
  else
    discharge(j)=discharge_sec(j,res_index(upstream))
    depth(j)=depth_sec(j,res_index(upstream))
    hydrad(j)=hydrad_sec(j,res_index(upstream))
    energslope(j)=energslope_sec(j,res_index(upstream))
    manning(j)=manning_sec(j,res_index(upstream))
    topwidth(j)=topwidth_sec(j,res_index(upstream))
    meanvel(j)=meanvel_sec(j,res_index(upstream))
  endif
enddo

! loop to calculate the sediment carrying capacity
DO j=1,nbrsec(upstream)

 if (discharge(j) /= 0.) then
  
!Ge mean diameter has to be calculated from the grain size distribution
  dummy2=d50_actlay(j,res_index(upstream))

! 1) BED LOAD FORMULA
  
  
  
! Gravitacional acceleration (m/s2)
! 9.807
  
! density of water and density of natural sediments (kg/m3)
  wat_dens=1.*1000.
  sed_dens=2.65*1000.

! specific gravity of sediments in water
  spec_gravit=(sed_dens-wat_dens)/wat_dens
  
! critical shear stress for D50
  crshear_D50=.05  
    
! loop to calculate the sediment carrying capacity for bed load
  DO g=1,n_sed_class
    if(dummy2==0.) then
	  d50=diam(g)
	else
	  d50=dummy2
	endif

! critical shear stress for D50
    crshearvel_D50=sqrt(crshear_D50*spec_gravit*9.807*d50)

! shear velocity (m/s)
    shear_vel=sqrt(9.807*hydrad(j)*energslope(j))

! nondimensional shear stress ti
	shear(g)=(shear_vel**2)/(spec_gravit*9.807*diam(g))

! effective shear velocity
    effshear_vel=meanvel(j)/(5.75*log10(hydrad(j)/D50)/(1.+(2.*shear(g)))+6.)

! critical shear velocity Uci (m/s)
    if (diam(g)/D50 < .4) then
	  crshear_vel(g)=crshearvel_D50*sqrt(.85)
	else
	  crshear_vel(g)=crshearvel_D50*(log10(19.)/log10(19.*diam(g)/D50))
	endif

! nondimensional critical shear stress tci
    crshear(g)=(crshear_vel(g)**2)/(spec_gravit*9.807*diam(g))

! nondimensional effective shear stress tei
    effshear(g)=(effshear_vel**2)/(spec_gravit*9.807*diam(g))


! dimensionless transport parameter for fractional bed load yields (-)
    if ((1.-(crshear(g)/shear(g))) <= 0.) then
	  transpfac_bed(g)=0.
	else
      transpfac_bed(g)=1.-(crshear(g)/shear(g))
    endif

! fractional bed sediment transport (m2/s)
    if (shear(g)/=0.) then
      bed_frtransp(g,j)=17.*effshear_vel*diam(g)*effshear(g)*(transpfac_bed(g))* &
        (1.-sqrt((crshear(g)/shear(g))))
	else
	  bed_frtransp(g,j)=0.
	endif
	bed_frtransp(g,j)=max(0.,bed_frtransp(g,j))

! fractional bed sediment transport (m3/day)
    bed_frtransp(g,j)=86400.*bed_frtransp(g,j)* topwidth(j)
    
!write(*,'(3I4,6F10.6)')d,j,g,d50,crshear(g),transpfac_bed(g)
  END DO

  
! 2) SUSPENDED LOAD FORMULA
  
!Ge to include the DAILY mean temperature of the reservoir (celsius degree)
!Ge tempres=20 C (temporarily)
  tempres(step,upstream)=20.
  
! kinematic viscosity (m2/s)
  visc=(1.14-0.031*(tempres(step,upstream)-15.)+  &
      0.00068*((tempres(step,upstream)-15.)**2.))*1.e-6
  
! Von Kàrmàn constant (-)
	karman=0.412

! loop to calculate the sediment carrying capacity for suspended load
  DO g=1,n_sed_class
   IF (diam(g)*1000. > 0.040) THEN
    
! settling velocity (m/s)
    setvel(g)=SQRT((13.95*visc/diam(g))**2.+1.09*  &
        ((sed_dens-wat_dens)/wat_dens)*9.807*diam(g)) -13.95*visc/diam(g)
    
! C1 parameter
    C1=6.*setvel(g)/(karman*shear_vel*depth(j))

! seta initial parameter
	seta0=setvel(g)/(.75*shear_vel)

! f_seta function of the seta initial parameter
    f_seta0=(1./sqrt(2.*3.14))*exp(-.5*(seta0**2.))

! f_seta function of the seta parameter
    dummy6=seta0
	dummy7=f_seta0
	h=0
	f_seta=0.
	dummy=0.
!	dummy5=0.
!    do while (dummy>.0001 .or. dummy==0 .or. dummy5<dummy)
    do while (dummy6<5.)
!      if (h/=0) dummy5=dummy
	  if (dummy6<.001) dummy3=dummy6*10.
	  if (dummy6>=.001) dummy3=dummy6+.05
	  dummy4=(1./sqrt(2.*3.14))*exp(-.5*(dummy3**2.))
	  if(dummy4==0.) then
	    f_seta=0.
	    exit
	  endif
	  dummy=(dummy7+dummy4)*(dummy3-dummy6)/2.
	  dummy6=dummy3
	  dummy7=dummy4
!	  if (dummy<.0001 .and. dummy5<dummy)exit
	  h=h+1
	  f_seta=f_seta+dummy
!write(*,'(3I4,F5.2,7ES9.1)')j,g,h,dummy3,dummy4,dummy5,dummy6,dummy7,dummy,f_seta
!if (j==51)stop
!if (h==100)stop
	enddo

  
! concentration at a=0,05h (-) 
!    Ca=.025*((f_seta0/seta0)-f_seta)
    if (seta0/=0.) Ca=.025*((f_seta0/seta0)-f_seta)
	Ca=max(Ca,0.)
!write(*,*)d,j,g,ca
!if (Ca<0.)write(*,*)d,j,g,ca
! depth a


    a=.05*depth(j)
    
! total suspended sediment transport (m2/s)
    if (Ca/=0.) then
      susp_frtransp(g,j)=Ca*meanvel(j)*(exp(-C1*a)-exp(-C1*depth(j)))*exp(C1*a)/C1
    else
	  susp_frtransp(g,j)=0.
	endif

    susp_frtransp(g,j)=min(susp_frtransp(g,j),1000.)

!write(*,'(2I4,I5,F10.2,6ES11.3)')j,g,h,seta0,f_seta0,dummy,f_seta,Ca,susp_frtransp(g,j)

! constrained by bed shear stress  (-)
!    if ((1-(crshear(g)/shear(g))) <= 0.) susp_frtransp(g,j)=0.

   ELSE
    susp_frtransp(g,j)=1000.

    setvel(g)=SQRT((13.95*visc/diam(g))**2.+1.09*  &
        ((sed_dens-wat_dens)/wat_dens)*9.807*diam(g)) -13.95*visc/diam(g)
   ENDIF
   


! total suspended sediment transport (m3/day)
   susp_frtransp(g,j)=86400.*susp_frtransp(g,j)* topwidth(j)

    
    
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

!if (d==1 .and. j==1)stop
!if (d==2)stop
END DO

RETURN
END SUBROUTINE eq2_ashida
