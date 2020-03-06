SUBROUTINE lake_routing(muni,k,help,help2,help3)

! Code converted using TO_F90 by Alan Miller
! Date: 2005-08-23  Time: 12:57:31


use lake_h
use common_h
use time_h
use reservoir_h
use params_h

IMPLICIT NONE

INTEGER, INTENT(IN)               :: muni
INTEGER, INTENT(IN)               :: k
REAL, INTENT(IN)                  :: help
REAL, INTENT(IN)                  :: help2
REAL, INTENT(OUT)                  :: help3

INTEGER :: i,n !,j,npts
INTEGER :: ninterac !,nbrbat1

REAL :: par_c,par_d !,par_a,par_b
REAL :: y1,y2,vol1,vol2,dvol !,elev,out1,out2
REAL :: frtime,inflow_t,inflow_mean,outflow_mean
REAL :: inflow,outflow
REAL :: dummy1,dummy2,error !,dummy3,dummy4,dummy5
REAL :: volmax,volstep,time !,areamax
REAL :: hmax0,vmax0,par_k,par_alpha
REAL :: dummy
! -----------------------------------------------------------------------

! Initialization
! Initialization
if (acud(muni,k)>0.) then
 inflow=help/(86400./nt)
 volstep=help2

 inflow=inflow/acud(muni,k)
 volstep=volstep/acud(muni,k)

 outflow=outflow_last_hrr(muni,k)
 volmax=maxlake(muni,k)
 par_c=damc_hrr(k)
 par_d=damd_hrr(k)
 par_k=damk_Molle(k)
 par_alpha=alpha_Molle(k)
 vmax0=maxlakesub0(muni,k)
 hmax0=hmax_hrr(muni,k)

!write(*,'(2I6,2F7.3,F13.1,2F8.3,F13.1)')muni,k,par_c,par_d,par_k,par_alpha,hmax0,vmax0

!stop


 IF (outflow==0.) THEN
  vol2=max(0.,volstep-volmax)
 ELSE
  vol2=volume_last_hrr(muni,k)
 ENDIF

 dummy=vol2
!write(*,'(10F20.6)') vol2,volume_last(upstream),volstep,volmax,inflow,outflow

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
 inflow_mean=0.
 outflow_mean=0.


 DO i=1,ninterac
  frtime=frtime+(time/ninterac)
  inflow_t=inflow*(time/ninterac)

!write(*,'(I6,10F20.6)') i,frtime,inflow_t

  y1=y2
  vol1=max((par_k*((hmax0+y1)**par_alpha))-vmax0,0.)

!write(*,'(I6,10F13.3)') i,y1,y2,vol1,vol2,volmax

  n=0
  DO WHILE (ABS(error) > 0.001)
    n=n+1
	dvol=inflow_t-((outflow+dummy1)*(time/ninterac)/2.)
	vol2=max(vol1+dvol,0.)
!write(*,'(5F12.3)')frtime,inflow,outflow,dummy1,dvol
	y2=max((((vol2+vmax0)/par_k)**(1./par_alpha))-hmax0,0.)
	dummy1=par_c*(y2**par_d)

!if(upstream==29)write(*,'(I4,8F9.1)')upstream,frtime,y1,y2,vol1,vol2,dummy1,hmax0,vmax0
! boundary condition for constant inflow
!******************************************************************************************
    IF (dummy1*(time/ninterac)>inflow_t*(time/ninterac)+vol1 .or. y2==0. .or. n>=10) THEN
	 dummy1=inflow
     y2=0.
	 vol2=0.
	 error=0.
!write(*,*)frtime,inflow,dummy1
!write(*,'(2I6,10F20.6)') i,n,frtime,inflow,inflow_t,y1,vol1,dvol,vol2,y2,dummy1,error
!if(muni==2 .and. k==1)write(*,'(2I4,F10.1,F10.3,3F10.1,2F8.3)')muni,k,frtime,inflow,inflow_t,vol1,vol2,y2,dummy1
	 exit
!else
!if(muni==2 .and. k==1)write(*,'(2I4,F10.1,F10.3,3F10.1,2F8.3)')muni,k,frtime,inflow,inflow_t,vol1,vol2,y2,dummy1
	ENDIF
!******************************************************************************************
	error=y2-dummy2
	dummy2=y2
!if(muni==2 .and. k==1)write(*,'(2I4,F10.1,F10.3,3F10.1,2F8.3)') muni,k,frtime,inflow,inflow_t,vol1,vol2,y2,dummy1
  ENDDO
!write(*,*)frtime,inflow,dummy1
  outflow=dummy1
  outflow_mean=outflow_mean+outflow*(time/ninterac)
  inflow_mean=inflow_mean+inflow*(time/ninterac)
  error=1.


  outflow_last_hrr(muni,k)=outflow
 ENDDO
! help3=outflow_mean*acud(muni,k)
 help3=(inflow_mean+(dummy-vol2))*acud(muni,k)
!if(muni==2 .and. k==1) write(*,*)inflow_mean,dummy,vol2,help3
 volume_last_hrr(muni,k)=vol2
 lakewater(step,muni,k)=(volmax+vol2)*acud(muni,k)
 lakeoutflow_hrr(step,muni,k)=(inflow_mean+(dummy-vol2))*acud(muni,k)
! lakeoutflow_hrr(step,muni,k)=outflow_mean

!write(*,'(2I4,10F10.1)')muni,k,frtime,inflow,help3/86400.,volume_last_hrr(muni,k),lakewater(step,muni,k)
ELSE
 help3=help
 volume_last_hrr(muni,k)=0.
 outflow_last_hrr(muni,k)=0.
 lakewater(step,muni,k)=0.
ENDIF
RETURN
END SUBROUTINE lake_routing
