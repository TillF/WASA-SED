SUBROUTINE reservoir_routing(upstream,help,help2)

! Code converted using TO_F90 by Alan Miller
! Date: 2005-08-23  Time: 12:57:31

use common_h
use time_h
use reservoir_h
use params_h

IMPLICIT NONE

INTEGER, INTENT(IN OUT)               :: upstream
REAL, INTENT(IN OUT)                  :: help
REAL, INTENT(IN OUT)                  :: help2

INTEGER :: i,n !j,k,,npts
INTEGER :: ninterac !,nbrbat1

REAL :: par_c,par_d
REAL :: y1,y2,vol1,vol2,dvol !out1,out2,elev
REAL :: frtime,inflow_t,outflow_mean
REAL :: inflow,outflow
REAL :: dummy1,dummy2,error !,dummy3,dummy4,dummy5
REAL :: volmax,volstep,time !,areamax
REAL :: hmax0,vmax0,par_k,par_alpha
! -----------------------------------------------------------------------

! Initialization
inflow=help/(86400./nt)
outflow=outflow_last(upstream)
volmax=daystorcap(step,res_index(upstream)) !storage capacity in the subbasin's reservoir [m**3]
volstep=help2
par_c=damc(upstream)
par_d=damd(upstream)
par_k=k_over(res_index(upstream))
par_alpha=alpha_over(res_index(upstream))
vmax0=(storcap(upstream)*1.e6)
hmax0=hmax(res_index(upstream))

!write(*,'(I6,4F7.3,F13.1,2F8.3,F13.1)')upstream,par_c,par_d,par_k,par_alpha,hmax0,vmax0

!stop

! if this is the initial timestep, volume_last might be > 0 even if outflow_last = 0
! because there is no information about outflow_last at initial timestep but volume_last might be > 0 depending on parameter vol0
IF (outflow==0. .and. dtot > 1) THEN
  vol2=0.
ELSE
  vol2=volume_last(upstream)
ENDIF

!write(*,'(10F12.1)') vol2,volume_last(upstream),volstep,volmax,inflow,outflow

! initial condition from previous timestep
y2=max((((vol2+vmax0)/par_k)**(1./par_alpha))-hmax0,0.)

!write(*,'(10F20.6)') y2,vol2,volstep,volmax,inflow,outflow



IF (vol2 > 0.) THEN
  time=(86400./nt)
ELSE
  time=(86400./nt)-((volmax-volstep)/inflow)
ENDIF

!write(*,'(10F20.6)') time,volmax,volstep,inflow



!write(*,'(2I6,10F20.6)') i,n,y2,elev,volact(step,upstream),inflow,outflow

ninterac=60
frtime=0.
error=1.
dummy1=outflow
dummy2=y2
outflow_mean=0.

! loop over sub-steps
DO i=1,ninterac
  frtime=frtime+(time/ninterac)
  inflow_t=inflow*(time/ninterac)

!write(*,'(I6,10F20.6)') i,frtime,inflow_t

  y1=y2
  vol1=max((par_k*((hmax0+y1)**par_alpha))-vmax0,0.)

!write(*,'(I6,10F13.3)') i,y1,y2,vol1,vol2,volmax

  ! predictor-corrector until convergence (?)
  n=0
  DO WHILE (ABS(error) > 0.001)
    n=n+1
	dvol=inflow_t-((outflow+dummy1)*(time/ninterac)/2.)
	vol2=max(vol1+dvol,0.)
!write(*,'(5F12.3)')frtime,inflow,outflow,dummy1,dvol
	y2=max((((vol2+vmax0)/par_k)**(1./par_alpha))-hmax0,0.)
	dummy1=par_c*(y2**par_d) ! outflow m3/s

!if(upstream==29)write(*,'(I4,8F9.1)')upstream,frtime,y1,y2,vol1,vol2,dummy1,hmax0,vmax0
! boundary condition for constant inflow
!******************************************************************************************
    IF (dummy1*(time/ninterac)>inflow_t*(time/ninterac)+vol1 .or. y2==0. .or. n>=10) THEN !tp corrected units, Feb 2018
	 dummy1=inflow
     y2=0.
	 vol2=0.
	 error=0.
!write(*,*)frtime,inflow,dummy1
!write(*,'(2I6,10F20.6)') i,n,frtime,inflow,inflow_t,y1,vol1,dvol,vol2,y2,dummy1,error
!if(upstream==29)write(*,'(I4,F10.1,F10.3,3F12.3,2F8.3)')upstream,frtime,inflow,inflow_t,vol1,vol2,y2,dummy1
	 exit
!else
!if(upstream==29)write(*,'(I4,F10.1,F10.3,3F12.3,2F8.3)')upstream,frtime,inflow,inflow_t,vol1,vol2,y2,dummy1
	ENDIF
!******************************************************************************************
	error=y2-dummy2
	dummy2=y2
!write(*,'(2I6,10F20.6)') i,n,frtime,inflow,inflow_t,y1,vol1,dvol,vol2,y2,dummy1,error
  ENDDO
!write(*,*)frtime,inflow,dummy1
  outflow=dummy1
  outflow_mean=outflow_mean+outflow*(time/ninterac) ! cumulated outflow over sub-steps in m3
  error=1.

  outflow_last(upstream)=outflow
ENDDO

overflow(step,res_index(upstream))=outflow_mean
volume_last(upstream)=vol2
volact(step,upstream)=volmax+vol2 ! current reservoir volume in m3

!write(*,'(I4,10F15.3)')upstream,frtime,overflow(step,upstream)/86400.,outflow_last(upstream),volact(step,upstream)

RETURN
END SUBROUTINE reservoir_routing
