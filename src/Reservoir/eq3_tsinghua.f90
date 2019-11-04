SUBROUTINE eq3_tsinghua(upstream)

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
INTEGER :: j,g,factor_sed !i,h,id,ih,irout,imun,dummy1,
!real:: dummy2


real :: wat_dens,sed_dens
real :: visc


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
    energslope(j)=(energslope_sec(j,res_index(upstream))+energslope_sec(j+1,res_index(upstream)))/2.
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

! TOTAL LOAD FORMULA

! loop to calculate the sediment carrying capacity
DO j=1,nbrsec(upstream)

  if (discharge(j) /= 0.) then
    DO g=1,n_sed_class
      if(diam(g)*1000.<=0.01) factor_sed=1600
      if(diam(g)*1000.>0.01 .and. diam(g)*1000.<=0.1) factor_sed=650
      if(diam(g)*1000.>0.1  .and. diam(g)*1000.<=2.)  factor_sed=300
      if(diam(g)*1000.>2.) factor_sed=180
	enddo  

! total suspended sediment transport (m3/day)
    DO g=1,n_sed_class
      fr_capacity(g,j)=factor_sed*(discharge(j)**1.6)*(energslope(j)**1.2)/(topwidth(j)*.6)
      fr_capacity(g,j)=fr_capacity(g,j)*86400./1.65
    END DO

  else
    DO g=1,n_sed_class
      fr_capacity(g,j)=0.
	enddo
  endif

! Gravitacional acceleration (m/s2)
! 9.807
  
! density of water and density of natural sediments (kg/m3)
  wat_dens=1.*1000.
  sed_dens=2.65*1000.


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
  enddo

  DO g=1,n_sed_class
!write(*,'(3I4,3F20.5)')d,j,g,fr_capacity(g,j),susp_frtransp(g,j),bed_frtransp(g,j)
  END DO

!if (d==54)read(*,*)
END DO

RETURN
END SUBROUTINE eq3_tsinghua
