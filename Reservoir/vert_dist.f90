SUBROUTINE vert_dist(g,j,upstream,fbottom,fintake,fover)
 
! Till: computationally irrelevant: minor changes to improve compiler compatibility
! 2011-04-29

! Code converted using TO_F90 by Alan Miller
! Date: 2005-08-23  Time: 12:57:31
 
use common_h
use time_h
use reservoir_h


implicit none

INTEGER, INTENT(IN OUT)                  :: upstream
INTEGER:: g,j
INTEGER:: i,numb_interv
REAL:: fbottom,fintake,fover
REAL:: a,dep,z
REAL:: factor,interv,depth_bottom,depth_intake,depth_over
REAL:: wat_dens,sed_dens,visc,shear_vel,karman,tempres(d,upstream)
REAL:: retention_s,time_s,suspension_0

!Ge to include the DAILY mean temperature of the reservoir (celsius degree)
!Ge tempres=20 C (temporarily)
tempres(d,upstream)=20.

! initialization
a=2.*diam(g)
depth_bottom=a
depth_intake=max(0.,elevdead(upstream)-dayminlevel(step,upstream))
depth_over=max(0.,(1.-a)*(maxlevel(upstream)-dayminlevel(step,upstream)))

! density of water and density of natural sediments (kg/m3)
wat_dens=1.*1000.
sed_dens=2.65*1000.
  
! kinematic viscosity (m2/s)
visc=(1.14-0.031*(tempres(d,upstream)-15.)+  &
      0.00068*((tempres(d,upstream)-15.)**2.))*1.e-6

! settling velocity (m/s)
setvel(g)=SQRT((13.95*visc/diam(g))**2.+1.09*  &
        ((sed_dens-wat_dens)/wat_dens)*9.807*diam(g))-13.95*visc/diam(g)
!write(*,*)upper_limit(g),diam(g)*1000.,setvel(g)
!if (upper_limit(g) <= .0025) setvel(g)=0.00000456
if (upper_limit(g) <= .0040) setvel(g)=0.0000117

if (discharge_sec(j,upstream) /= 0.) then

! shear velocity (m/s)
  shear_vel=sqrt(9.807*hydrad_sec(j,upstream)*energslope_sec(j,upstream))
!  shear_vel=sqrt(9.807*hydrad_sec(j,upstream)*max(energslope_sec(j,upstream),.0001))
!write(*,'(2I4,8E10.2)')j,g,a,hydrad_sec(j,upstream),energslope_sec(j,upstream),shear_vel

! Von Kàrmàn constant (-)
  karman=0.412

  dep=max(0.,maxlevel(upstream)-dayminlevel(step,upstream))
  z=setvel(g)/(karman*shear_vel)


  if (depth_bottom /= 0.) fbottom=(((dep-depth_bottom)/depth_bottom)*(a/(dep-a)))**z
  if (depth_intake /= 0.) fintake=(((dep-depth_intake)/depth_intake)*(a/(dep-a)))**z
  if (depth_over /= 0.) fover=(((dep-depth_over)/depth_over)*(a/(dep-a)))**z

else
  fbottom=0.
  fintake=0.
  fover=0.
endif

!write(*,'(2I4,8F10.3)')j,g,depth_bottom,depth_intake,depth_over,fbottom,fintake,fover
!write(*,'(2I4,8E12.4)')j,g,a,dep,z,depth_over,shear_vel,setvel(g)
!if (step==2)stop
END SUBROUTINE vert_dist

