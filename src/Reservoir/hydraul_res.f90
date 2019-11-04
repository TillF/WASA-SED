SUBROUTINE hydraul_res(upstream)
!Till: computationally irrelevant: outcommented unused vars
!2012-09-14

!Till: computationally irrelevant: minor changes to improve compiler compatibility
!2011-04-29

! Code converted using TO_F90 by Alan Miller
! Date: 2005-08-23  Time: 12:57:41

use common_h
use routing_h
use time_h
use reservoir_h

IMPLICIT NONE

INTEGER, INTENT(IN OUT)                  :: upstream

INTEGER :: j,n,p,k,l,dummy1,dummy2,dummy3,q !i,id,ih,irout,imun,nbrsec1,
INTEGER :: npt
INTEGER :: a,b,e,f !,c,c2
INTEGER :: interv

REAL :: damelev,TAN,dist0_sec,discharge_calc !,side_slope
!REAL :: area1,elev1
REAL :: loclosscoef,coef,weight !maxelev,
REAL ::  lowlim,toplim,dep,error,error1
REAL ::  area_part(200,200) !,centroid(200,200)
REAL :: flow_regime,normaldepth(200)

!temporary
!INTEGER :: dummy(200)

REAL :: dummy4


! -----------------------------------------------------------------------

! initialization of hydraulic parameters
resvol(upstream)=0.
DO j=1,nbrsec(upstream)
  resvol_sec(j,res_index(upstream))=0.
  bedslope_sec(j,res_index(upstream))=0.
  area_sec(j,res_index(upstream))=0.
  resarea_sec(j,res_index(upstream))=0.
  hydrad_sec(j,res_index(upstream))=0.
  meanvel_sec(j,res_index(upstream))=0.
  energslope_sec(j,res_index(upstream))=0.
  topwidth_sec(j,res_index(upstream))=0.
  depth_sec(j,res_index(upstream))=0.
  discharge_sec(j,res_index(upstream))=0.
  weight_sec(j,res_index(upstream))=0.
  wetper_sec(j,res_index(upstream))=0.
  meanvel_sec(j,res_index(upstream))=0. 
  energslope_sec(j,res_index(upstream))=0. 
  dynhead_sec(j,res_index(upstream))=0.
  tothead_sec(j,res_index(upstream))=0.
  headloss_sec(j,res_index(upstream))=0.
  locloss_sec(j,res_index(upstream))=0.
  calctothead_sec(j,res_index(upstream))=0.
  maxarea_sec(j,res_index(upstream))=0.
  maxelev_sec(j,res_index(upstream))=0.
  maxdepth_sec(j,res_index(upstream))=0.
  crdepth_sec(j,res_index(upstream))=0.
  crwatelev_sec(j,res_index(upstream))=0.
  crarea_sec(j,res_index(upstream))=0.
  crtopwidth_sec(j,res_index(upstream))=0.
  crwetper_sec(j,res_index(upstream))=0.
END DO

! Mean reservoir level (m)
IF (reservoir_balance == 1) THEN
  damelev0(step,upstream)= -1. !flag as "not set"
  DO l=1,(nbrbat(upstream)-1)
    IF (damvol0(upstream) >= vol_bat(l,upstream).AND.  &
        damvol0(upstream) <= vol_bat(l+1,upstream)) THEN
      damelev0(step,upstream)=elev_bat(l,upstream)+(damvol0(upstream)-  &
        vol_bat(l,upstream))/(vol_bat(l+1,upstream)- vol_bat(l,upstream))*  &
        (elev_bat(l+1,upstream)-elev_bat(l,upstream))
    END IF
  END DO
  if (damelev0(step,upstream) == -1.) then !flag as "not set"
    write(*,*)"ERROR: actual volume beyond range of cav.dat"
    !damelev0(step,upstream)= 1.
    stop
  end if
ENDIF

damelev=(damelevact(upstream)+damelev0(step,upstream))/2.
damelev_mean(step,res_index(upstream))=damelev
!write(*,*)damelevact(upstream),damelev0(step,upstream),damelev_mean(step,upstream)

! Minimum elevation for each cross_section (m)
DO j=1,nbrsec(upstream)
  minelev_sec(j,res_index(upstream))=y_sec(1,j,res_index(upstream))
  npt=npoints(j,res_index(upstream))
  DO m=2,npt
    IF (minelev_sec(j,res_index(upstream)) > y_sec(m,j,res_index(upstream))) THEN
      minelev_sec(j,res_index(upstream))=y_sec(m,j,res_index(upstream))
	  x_minelev(j,res_index(upstream))=x_sec(m,j,res_index(upstream))
    END IF
  END DO
END DO

! Bed slope of each cross section
DO j=1,nbrsec(upstream)
  IF (j /= nbrsec(upstream)) THEN
    IF (minelev_sec(j,res_index(upstream))-minelev_sec(j+1,res_index(upstream)) <= 0.) THEN
      bedslope_sec(j,res_index(upstream))=0.001
    ELSE
      bedslope_sec(j,res_index(upstream))=(minelev_sec(j,res_index(upstream))-  &
          minelev_sec(j+1,res_index(upstream)))/dist_sec(j,res_index(upstream))
	  bedslope_sec(j,res_index(upstream))=max(bedslope_sec(j,res_index(upstream)),0.001)
    END IF
  ELSE
    bedslope_sec(j,res_index(upstream))=bedslope_sec(j-1,res_index(upstream))
  END IF
END DO

! Computation of maximum wetted surface of cross sections
DO j=1,nbrsec(upstream)
  maxarea_sec(j,res_index(upstream))=0.
  npt=npoints(j,res_index(upstream))
  maxelev_sec(j,res_index(upstream))=MIN(y_sec(1,j,res_index(upstream)), y_sec(npt,j,res_index(upstream)))
  maxdepth_sec(j,res_index(upstream))=maxelev_sec(j,res_index(upstream))- minelev_sec(j,res_index(upstream))

  DO m=2,npt
    TAN=ABS(y_sec(m,j,res_index(upstream))-y_sec(m-1,j,res_index(upstream)))/  &
        (x_sec(m,j,res_index(upstream))-x_sec(m-1,j,res_index(upstream)))
    IF (maxelev_sec(j,res_index(upstream)) >= y_sec(m-1,j,res_index(upstream))  &
          .AND.maxelev_sec(j,res_index(upstream)) >= y_sec(m,j,res_index(upstream))) THEN
      maxarea_sec(j,res_index(upstream))=maxarea_sec(j,res_index(upstream))+  &
          (x_sec(m,j,res_index(upstream))-x_sec(m-1,j,res_index(upstream)))  &
          *(maxelev_sec(j,res_index(upstream))-  &
          (y_sec(m,j,res_index(upstream))+y_sec(m-1,j,res_index(upstream)))/2.)
    ELSE IF (maxelev_sec(j,res_index(upstream)) < y_sec(m-1,j,res_index(upstream))  &
          .AND.maxelev_sec(j,res_index(upstream)) >= y_sec(m,j,res_index(upstream))) THEN
      maxarea_sec(j,res_index(upstream))=maxarea_sec(j,res_index(upstream))+  &
          (((maxelev_sec(j,res_index(upstream))-y_sec(m,j,res_index(upstream)))**2.)/ (2.*TAN))
    ELSE IF (maxelev_sec(j,res_index(upstream)) >= y_sec(m-1,j,res_index(upstream))  &
          .AND.maxelev_sec(j,res_index(upstream)) < y_sec(m,j,res_index(upstream))) THEN
      maxarea_sec(j,res_index(upstream))=maxarea_sec(j,res_index(upstream))+  &
          (((maxelev_sec(j,res_index(upstream))-y_sec(m-1,j,res_index(upstream)))**2.)/ (2.*TAN))
    END IF
  END DO
END DO

! Check if the cross sections are completely sedimented
DO j=1,nbrsec(upstream)
  if (maxarea_sec(j,res_index(upstream)) == 0.) then
    WRITE(*,*)'ERROR: cross section ',j
    WRITE(*,*)'is completely sedimented'
    STOP
  END IF
END DO



! Check if the cross section belongs to either the river or the reservoir
n=0

! Computation of the water discharge along the longitudinal profile of the reservoir
! from a weighted averaged value using the river inflow and the reservoir outflow discharges
DO j=1,nbrsec(upstream)
  IF (damelev >= minelev_sec(j,res_index(upstream))) THEN
    IF (damelev > maxelev_sec(j,res_index(upstream))) THEN
      WRITE(*,*)'ERROR: the water discharge overflows the cross section geometry ',j
      STOP
    END IF

! RESERVOIR SUBREACH
    n=n+1
! Cross section area
    resarea_sec(j,res_index(upstream))=0.
    npt=npoints(j,res_index(upstream))
    DO m=2,npt
      TAN=ABS(y_sec(m,j,res_index(upstream))-y_sec(m-1,j,res_index(upstream)))/  &
          (x_sec(m,j,res_index(upstream))-x_sec(m-1,j,res_index(upstream)))
      IF (damelev >= y_sec(m-1,j,res_index(upstream)).AND.damelev >=  &
            y_sec(m,j,res_index(upstream))) THEN
        resarea_sec(j,res_index(upstream))=resarea_sec(j,res_index(upstream))+  &
            (x_sec(m,j,res_index(upstream))-x_sec(m-1,j,res_index(upstream)))*(damelev-  &
            (y_sec(m,j,res_index(upstream))+y_sec(m-1,j,res_index(upstream)))/2.)
      ELSE IF (damelev < y_sec(m-1,j,res_index(upstream)).AND.damelev >=  &
            y_sec(m,j,res_index(upstream))) THEN
        resarea_sec(j,res_index(upstream))=resarea_sec(j,res_index(upstream))+  &
            (((damelev-y_sec(m,j,res_index(upstream)))**2.)/(2.*TAN))
      ELSE IF (damelev >= y_sec(m-1,j,res_index(upstream)).AND.damelev <  &
            y_sec(m,j,res_index(upstream))) THEN
        resarea_sec(j,res_index(upstream))=resarea_sec(j,res_index(upstream))+  &
            (((damelev-y_sec(m-1,j,res_index(upstream)))**2.)/(2.*TAN))
      END IF
    END DO

! Reservoir's volume represented by two consecutive cross section
    IF (n == 1 .and. j /= 1) THEN
      dist0_sec=(damelev-minelev_sec(j,res_index(upstream)))/ bedslope_sec(j-1,res_index(upstream))
      resvol_sec(j,res_index(upstream))=resarea_sec(j,res_index(upstream))*dist0_sec/3.
    ELSE IF (n == 1 .and. j == 1) THEN
      dist0_sec=(damelev-minelev_sec(j,res_index(upstream)))/ bedslope_sec(j,res_index(upstream))
      resvol_sec(j,res_index(upstream))=resarea_sec(j,res_index(upstream))*dist0_sec/3.
    ELSE
      resvol_sec(j,res_index(upstream))=(dist_sec(j-1,res_index(upstream))/3.)*  &
          (resarea_sec(j,res_index(upstream))+resarea_sec(j-1,res_index(upstream))+  &
          SQRT(resarea_sec(j,res_index(upstream))*resarea_sec(j-1,res_index(upstream))))
    END IF
    resvol(upstream)=resvol(upstream)+resvol_sec(j,res_index(upstream))

  ELSE
! RIVER SUBREACH
    resarea_sec(j,res_index(upstream))=0.
  END IF
END DO


! Reservoir subreach must not include sections with minimum elevation above the reservoir level
q=0
DO j=1,nbrsec(upstream)
  if (resarea_sec(j,res_index(upstream)) == 0.) q=j
ENDDO

if (q /= 0) then
  DO j=1,q
    dummy4=resarea_sec(j,res_index(upstream))
    resarea_sec(j,res_index(upstream))=0.
  enddo
endif

! Weighted averaged value
weight=0.
DO j=1,nbrsec(upstream)
  IF (damelev >= minelev_sec(j,res_index(upstream))) THEN
    weight=weight+resvol_sec(j,res_index(upstream))/resvol(upstream)
    weight_sec(j,res_index(upstream))=weight
    IF (j == nbrsec(upstream)) THEN
      weight_sec(j,res_index(upstream))=1.
    END IF
  END IF
END DO

! Discharge at the cross section j (m3/s)
DO j=1,nbrsec(upstream)
  IF (resarea_sec(j,res_index(upstream)) /= 0.) THEN
! reservoir subreach
    discharge_sec(j,res_index(upstream))=qinflow(step,upstream)-  &
        (qinflow(step,upstream)-res_qout(step,res_index(upstream)))*weight_sec(j,res_index(upstream))
  ELSE
! river subreach
    discharge_sec(j,res_index(upstream))=(qinflow(step,upstream)+res_qout(step,res_index(upstream)))/2.
  END IF
END DO

! 1) calculation of the normal depth at cross section by trial-and-error

!write(*,*)'Step (1)'

! Find the amount of cross section located in the river subreach
k=0

DO j=1,nbrsec(upstream)
  npt=npoints(j,res_index(upstream))

  if (j==1) depth_sec(j,res_index(upstream))=1.
  if (j/=1) depth_sec(j,res_index(upstream))=depth_sec(j-1,res_index(upstream))

  watelev_sec(j,res_index(upstream))=depth_sec(j,res_index(upstream))+ minelev_sec(j,res_index(upstream))

! trial-and-error
  error=1.
  error1=0.
  toplim=0.
  lowlim=0.
  interv=0
  dep=0.

  a=0
  dummy3=0
  IF(discharge_sec(j,res_index(upstream)) > 0.) THEN
    DO WHILE (ABS(error) > 0.01)
      a=a+1
!dummy(j)=a

! hydraulic parameters at the first cross section
      area_sec(j,res_index(upstream))=0.
      topwidth_sec(j,res_index(upstream))=0.
      wetper_sec(j,res_index(upstream))=0.

      IF (watelev_sec(j,res_index(upstream)) > maxelev_sec(j,res_index(upstream))) THEN
        watelev_sec(j,res_index(upstream))=maxelev_sec(j,res_index(upstream))
        depth_sec(j,res_index(upstream))=maxdepth_sec(j,res_index(upstream))
      END IF

      DO m=2,npt
        TAN=ABS(y_sec(m,j,res_index(upstream))-y_sec(m-1,j,res_index(upstream)))/  &
            (x_sec(m,j,res_index(upstream))-x_sec(m-1,j,res_index(upstream)))
        IF (watelev_sec(j,res_index(upstream)) >= y_sec(m-1,j,res_index(upstream)).AND.  &
              watelev_sec(j,res_index(upstream)) >= y_sec(m,j,res_index(upstream))) THEN
          area_sec(j,res_index(upstream))=area_sec(j,res_index(upstream))+  &
              (x_sec(m,j,res_index(upstream))-x_sec(m-1,j,res_index(upstream)))*  &
              (watelev_sec(j,res_index(upstream))-  &
              (y_sec(m,j,res_index(upstream))+y_sec(m-1,j,res_index(upstream)))/2.)
          topwidth_sec(j,res_index(upstream))=topwidth_sec(j,res_index(upstream))+  &
              (x_sec(m,j,res_index(upstream))-x_sec(m-1,j,res_index(upstream)))
          wetper_sec(j,res_index(upstream))=wetper_sec(j,res_index(upstream))+  &
              SQRT((x_sec(m,j,res_index(upstream))-x_sec(m-1,j,res_index(upstream)))**2.  &
              +(y_sec(m,j,res_index(upstream))-y_sec(m-1,j,res_index(upstream)))**2.)
        ELSE IF (watelev_sec(j,res_index(upstream)) < y_sec(m-1,j,res_index(upstream))  &
              .AND.watelev_sec(j,res_index(upstream)) >= y_sec(m,j,res_index(upstream))) THEN
          area_sec(j,res_index(upstream))=area_sec(j,res_index(upstream))+  &
              (((watelev_sec(j,res_index(upstream))-y_sec(m,j,res_index(upstream)))**2.)/ (2.*TAN))
          topwidth_sec(j,res_index(upstream))=topwidth_sec(j,res_index(upstream))+  &
              ((x_sec(m,j,res_index(upstream))-x_sec(m-1,j,res_index(upstream)))*  &
              (watelev_sec(j,res_index(upstream))-y_sec(m,j,res_index(upstream)))/  &
              (y_sec(m-1,j,res_index(upstream))-y_sec(m,j,res_index(upstream))))
          wetper_sec(j,res_index(upstream))=wetper_sec(j,res_index(upstream))+  &
              SQRT(((watelev_sec(j,res_index(upstream))-y_sec(m,j,res_index(upstream)))/TAN)  &
              **2+(watelev_sec(j,res_index(upstream))-y_sec(m,j,res_index(upstream)))**2.)
        ELSE IF (watelev_sec(j,res_index(upstream)) >= y_sec(m-1,j,res_index(upstream))  &
              .AND.watelev_sec(j,res_index(upstream)) < y_sec(m,j,res_index(upstream))) THEN
          area_sec(j,res_index(upstream))=area_sec(j,res_index(upstream))+  &
              (((watelev_sec(j,res_index(upstream))-y_sec(m-1,j,res_index(upstream)))**2.)/ (2.*TAN))
          topwidth_sec(j,res_index(upstream))=topwidth_sec(j,res_index(upstream))+  &
              ((x_sec(m,j,res_index(upstream))-x_sec(m-1,j,res_index(upstream)))*  &
              (watelev_sec(j,res_index(upstream))-y_sec(m-1,j,res_index(upstream)))/  &
              (y_sec(m,j,res_index(upstream))-y_sec(m-1,j,res_index(upstream))))
          wetper_sec(j,res_index(upstream))=wetper_sec(j,res_index(upstream))+  &
              SQRT(((watelev_sec(j,res_index(upstream))-y_sec(m-1,j,res_index(upstream)))/TAN)  &
              **2+(watelev_sec(j,res_index(upstream))-y_sec(m-1,j,res_index(upstream)))**2.)
        END IF
      END DO


! check if there is neither inflow nor outflow at the reservoir
      hydrad_sec(j,res_index(upstream))=area_sec(j,res_index(upstream))/ wetper_sec(j,res_index(upstream))
      discharge_calc=(1./(manning_sec(j,res_index(upstream))))*  &
          SQRT(bedslope_sec(j,res_index(upstream)))*hydrad_sec(j,res_index(upstream))**(2./3.)  &
          *area_sec(j,res_index(upstream))

      error=discharge_calc-discharge_sec(j,res_index(upstream))

      IF (watelev_sec(j,res_index(upstream)) == maxelev_sec(j,res_index(upstream)).AND.  &
            discharge_calc < discharge_sec(j,res_index(upstream)) .AND. a/=1) THEN
        WRITE(*,*)'ERROR: the water discharge overflows the cross section ',j
        STOP
      END IF

      if (dummy3==1) then
	    error=0.00001
        watelev_sec(j,res_index(upstream))=minelev_sec(j,res_index(upstream))+  &
              depth_sec(j,res_index(upstream))
		normalelev_sec(j,res_index(upstream))=watelev_sec(j,res_index(upstream))
        exit
	  endif
      IF (ABS(error) > 0.01) THEN
        IF (interv == 0) THEN
          IF (error1 >= 0. .AND.error > 0.) THEN
            dep=depth_sec(j,res_index(upstream))
            depth_sec(j,res_index(upstream))=.80*depth_sec(j,res_index(upstream))
          ELSE IF (error1 <= 0. .AND. error < 0.) THEN
            dep=depth_sec(j,res_index(upstream))
            depth_sec(j,res_index(upstream))=1.20*depth_sec(j,res_index(upstream))
          ELSE
            interv=1
          END IF
        END IF

        IF (interv /= 0) THEN
          IF (error1 < 0. .AND. error < 0.) THEN
            lowlim=MAX(depth_sec(j,res_index(upstream)),dep)
          ELSE IF (error1 > 0. .AND. error > 0.) THEN
            toplim=MIN(depth_sec(j,res_index(upstream)),dep)
          ELSE
            lowlim=MIN(depth_sec(j,res_index(upstream)),dep)
            toplim=MAX(depth_sec(j,res_index(upstream)),dep)
          END IF
          dep=depth_sec(j,res_index(upstream))
          depth_sec(j,res_index(upstream))=(lowlim+toplim)/2.
        END IF

        error1=error

		IF (toplim-lowlim < 1e-4 .and. toplim/=0. .and. lowlim/=0.) then
		  dep=depth_sec(j,res_index(upstream))
          depth_sec(j,res_index(upstream))=(lowlim+toplim)/2.
		  dummy3=1
!		  write(*,*)'warning! Trial-and-error calculation at the reservoir routing'
!		  write(*,*)'stopped for an approximation error of ',error1
!		  write(*,*)'hidraulic calculation (step 1) and cross section',j
		ENDIF

        watelev_sec(j,res_index(upstream))=minelev_sec(j,res_index(upstream))+  &
              depth_sec(j,res_index(upstream))
		normalelev_sec(j,res_index(upstream))=watelev_sec(j,res_index(upstream))
		normalarea_sec(j,res_index(upstream))=area_sec(j,res_index(upstream))
      END IF

!if(j==1)write(*,'(2I4,11F10.5,I4)')j,a,depth_sec(j,upstream),dep, &
!     watelev_sec(j,upstream),area_sec(j,upstream),wetper_sec(j,upstream),hydrad_sec(j,upstream), &
!     discharge_calc,discharge_sec(j,upstream), &
!     error,lowlim,toplim,interv

    END DO
    IF (resarea_sec(j,res_index(upstream))==0.) k=k+1
    IF (resarea_sec(j,res_index(upstream))/=0.) then
      if (watelev_sec(j,res_index(upstream))>damelev) k=k+1
	ENDIF
	normaldepth(j)=depth_sec(j,res_index(upstream))

  END IF
END DO

! 2) calculation of the critical depth at cross section by trial-and-error

! the first approximation to calculate the water depth (triangular cross section)
!write(*,*)'Step (2)'

DO j=1,nbrsec(upstream)
  if(j==k+1 .and. k/=nbrsec(upstream))exit

  npt=npoints(j,res_index(upstream))
  crdepth_sec(j,res_index(upstream))=depth_sec(j,res_index(upstream))
  crwatelev_sec(j,res_index(upstream))=watelev_sec(j,res_index(upstream))

! trial-and-error
  error=1.
  error1=0.
  toplim=0.
  lowlim=0.
  interv=0
  dep=0.

! check if there is neither inflow nor outflow at the reservoir
  b=0
  dummy3=0
  IF (discharge_sec(j,res_index(upstream)) > 0.) THEN
    DO WHILE (ABS(error) > 0.01)
      b=b+1
!dummy(j)=b

! hydraulic parameters at the first cross section
      crarea_sec(j,res_index(upstream))=0.
      crtopwidth_sec(j,res_index(upstream))=0.
      crwetper_sec(j,res_index(upstream))=0.

      IF (watelev_sec(j,res_index(upstream)) > maxelev_sec(j,res_index(upstream))) THEN
        crwatelev_sec(j,res_index(upstream))=maxelev_sec(j,res_index(upstream))
        crdepth_sec(j,res_index(upstream))=maxdepth_sec(j,res_index(upstream))
      END IF

      DO m=2,npt
        TAN=ABS(y_sec(m,j,res_index(upstream))-y_sec(m-1,j,res_index(upstream)))/  &
            (x_sec(m,j,res_index(upstream))-x_sec(m-1,j,res_index(upstream)))
        IF (crwatelev_sec(j,res_index(upstream)) >= y_sec(m-1,j,res_index(upstream)).AND.  &
              crwatelev_sec(j,res_index(upstream)) >= y_sec(m,j,res_index(upstream))) THEN
          crarea_sec(j,res_index(upstream))=crarea_sec(j,res_index(upstream))+  &
              (x_sec(m,j,res_index(upstream))-x_sec(m-1,j,res_index(upstream)))*  &
              (crwatelev_sec(j,res_index(upstream))-  &
              (y_sec(m,j,res_index(upstream))+y_sec(m-1,j,res_index(upstream)))/2.)
          crtopwidth_sec(j,res_index(upstream))=crtopwidth_sec(j,res_index(upstream))+  &
              (x_sec(m,j,res_index(upstream))-x_sec(m-1,j,res_index(upstream)))
          crwetper_sec(j,res_index(upstream))=crwetper_sec(j,res_index(upstream))+  &
              SQRT((x_sec(m,j,res_index(upstream))-x_sec(m-1,j,res_index(upstream)))**2.  &
              +(y_sec(m,j,res_index(upstream))-y_sec(m-1,j,res_index(upstream)))**2.)
        ELSE IF (crwatelev_sec(j,res_index(upstream)) < y_sec(m-1,j,res_index(upstream))  &
              .AND.crwatelev_sec(j,res_index(upstream)) >= y_sec(m,j,res_index(upstream))) THEN
          crarea_sec(j,res_index(upstream))=crarea_sec(j,res_index(upstream))+  &
              (((crwatelev_sec(j,res_index(upstream))-y_sec(m,j,res_index(upstream)))**2.)/ (2.*TAN))
          crtopwidth_sec(j,res_index(upstream))=crtopwidth_sec(j,res_index(upstream))+  &
              ((x_sec(m,j,res_index(upstream))-x_sec(m-1,j,res_index(upstream)))*  &
              (crwatelev_sec(j,res_index(upstream))-y_sec(m,j,res_index(upstream)))/  &
              (y_sec(m-1,j,res_index(upstream))-y_sec(m,j,res_index(upstream))))
          crwetper_sec(j,res_index(upstream))=crwetper_sec(j,res_index(upstream))+  &
              SQRT(((crwatelev_sec(j,res_index(upstream))-y_sec(m,j,res_index(upstream)))/TAN)  &
              **2+(crwatelev_sec(j,res_index(upstream))-y_sec(m,j,res_index(upstream)))**2.)
        ELSE IF (crwatelev_sec(j,res_index(upstream)) >= y_sec(m-1,j,res_index(upstream))  &
              .AND.crwatelev_sec(j,res_index(upstream)) < y_sec(m,j,res_index(upstream))) THEN
          crarea_sec(j,res_index(upstream))=crarea_sec(j,res_index(upstream))+  &
              (((crwatelev_sec(j,res_index(upstream))-y_sec(m-1,j,res_index(upstream)))**2.)/ (2.*TAN))
          crtopwidth_sec(j,res_index(upstream))=crtopwidth_sec(j,res_index(upstream))+  &
              ((x_sec(m,j,res_index(upstream))-x_sec(m-1,j,res_index(upstream)))*  &
              (crwatelev_sec(j,res_index(upstream))-y_sec(m-1,j,res_index(upstream)))/  &
              (y_sec(m,j,res_index(upstream))-y_sec(m-1,j,res_index(upstream))))
          crwetper_sec(j,res_index(upstream))=crwetper_sec(j,res_index(upstream))+  &
              SQRT(((crwatelev_sec(j,res_index(upstream))-y_sec(m-1,j,res_index(upstream)))/TAN)  &
              **2+(crwatelev_sec(j,res_index(upstream))-y_sec(m-1,j,res_index(upstream)))**2.)
        END IF
      END DO

      discharge_calc=SQRT(9.807*(crarea_sec(j,res_index(upstream))**3.)/  &
          crtopwidth_sec(j,res_index(upstream)))


      error=discharge_calc-discharge_sec(j,res_index(upstream))


      IF (crwatelev_sec(j,res_index(upstream)) == maxelev_sec(j,res_index(upstream)).AND.  &
            discharge_calc < discharge_sec(j,res_index(upstream))) THEN
        WRITE(*,*)'ERROR: the water discharge overflows the cross section ',j
        STOP
      END IF

      if (dummy3==1) then
	    error=0.00001
        crwatelev_sec(j,res_index(upstream))=minelev_sec(j,res_index(upstream))+  &
              crdepth_sec(j,res_index(upstream))
	    crvel_sec(j,res_index(upstream))=discharge_sec(j,res_index(upstream))/crarea_sec(j,res_index(upstream))
	    crhydrad_sec(j,res_index(upstream))=crarea_sec(j,res_index(upstream))/crwetper_sec(j,res_index(upstream))
		crslope_sec(j,res_index(upstream))=((manning_sec(j,res_index(upstream))*crvel_sec(j,res_index(upstream)))**2)/ &
			  (crhydrad_sec(j,res_index(upstream))**(4./3.))
        exit
	  endif


      IF (ABS(error) > 0.01) THEN
        IF (interv == 0) THEN
          IF (error1 >= 0. .AND. error > 0.) THEN
            dep=crdepth_sec(j,res_index(upstream))
            crdepth_sec(j,res_index(upstream))=.80*crdepth_sec(j,res_index(upstream))
          ELSE IF (error1 <= 0. .AND. error < 0.) THEN
            dep=crdepth_sec(j,res_index(upstream))
            crdepth_sec(j,res_index(upstream))=1.20*crdepth_sec(j,res_index(upstream))
          ELSE
            interv=1
          END IF
        END IF

        IF (interv /= 0) THEN
          IF (error1 < 0. .AND. error < 0.) THEN
            lowlim=MAX(crdepth_sec(j,res_index(upstream)),dep)
          ELSE IF (error1 > 0. .AND. error > 0.) THEN
            toplim=MIN(crdepth_sec(j,res_index(upstream)),dep)
          ELSE
            lowlim=MIN(crdepth_sec(j,res_index(upstream)),dep)
            toplim=MAX(crdepth_sec(j,res_index(upstream)),dep)
          END IF
          dep=crdepth_sec(j,res_index(upstream))
          crdepth_sec(j,res_index(upstream))=(lowlim+toplim)/2.
        END IF

        error1=error

		IF (toplim-lowlim<1e-4 .and. toplim/=0. .and. lowlim/=0.) then
		  dep=crdepth_sec(j,res_index(upstream))
          crdepth_sec(j,res_index(upstream))=(lowlim+toplim)/2.
          dummy3=1
!		  write(*,*)'warning! Trial-and-error calculation at the reservoir routing'
!		  write(*,*)'stopped for an approximation error of ',error1
!		  write(*,*)'hidraulic calculation (step 2) and cross section',j
		ENDIF

        crwatelev_sec(j,res_index(upstream))=minelev_sec(j,res_index(upstream))+  &
              crdepth_sec(j,res_index(upstream))
	    crvel_sec(j,res_index(upstream))=discharge_sec(j,res_index(upstream))/crarea_sec(j,res_index(upstream))
	    crhydrad_sec(j,res_index(upstream))=crarea_sec(j,res_index(upstream))/crwetper_sec(j,res_index(upstream))
		crslope_sec(j,res_index(upstream))=((manning_sec(j,res_index(upstream))*crvel_sec(j,res_index(upstream)))**2)/ &
			  (crhydrad_sec(j,res_index(upstream))**(4./3.))
      END IF


   END DO
!write(*,'(2I4,7F10.5)')j,dummy(j),crdepth_sec(j,upstream),dep, &
!     crwatelev_sec(j,upstream),crarea_sec(j,upstream),crwetper_sec(j,upstream),hydrad_sec(j,upstream),error

    if(crwatelev_sec(j,res_index(upstream))>maxelev_sec(j,res_index(upstream))) then
      WRITE(*,*)'ERROR: the water discharge overflows the cross section ',j
      STOP
    END IF

  END IF


END DO


! 3) routing


!write(*,*)'Step (3.1)'
DO j=1,nbrsec(upstream)

! 3.1) computation of hydraulic parameters if the cross section j belongs to the reservoir subreach
  IF (resarea_sec(j,res_index(upstream)) /= 0.) THEN

    npt=npoints(j,res_index(upstream))

	  area_sec(j,res_index(upstream))=resarea_sec(j,res_index(upstream))
      depth_sec(j,res_index(upstream))=damelev-minelev_sec(j,res_index(upstream))
	  watelev_sec(j,res_index(upstream))=damelev

	  topwidth_sec(j,res_index(upstream))=0.
	  wetper_sec(j,res_index(upstream))=0.

      DO m=2,npt
        TAN=ABS(y_sec(m,j,res_index(upstream))-y_sec(m-1,j,res_index(upstream)))/  &
            (x_sec(m,j,res_index(upstream))-x_sec(m-1,j,res_index(upstream)))
        IF (watelev_sec(j,res_index(upstream)) >= y_sec(m-1,j,res_index(upstream)).AND.  &
              watelev_sec(j,res_index(upstream)) >= y_sec(m,j,res_index(upstream))) THEN
          topwidth_sec(j,res_index(upstream))=topwidth_sec(j,res_index(upstream))+  &
              (x_sec(m,j,res_index(upstream))-x_sec(m-1,j,res_index(upstream)))
          wetper_sec(j,res_index(upstream))=wetper_sec(j,res_index(upstream))+  &
              SQRT((x_sec(m,j,res_index(upstream))-x_sec(m-1,j,res_index(upstream)))**2.  &
              +(y_sec(m,j,res_index(upstream))-y_sec(m-1,j,res_index(upstream)))**2.)
        ELSE IF (watelev_sec(j,res_index(upstream)) < y_sec(m-1,j,res_index(upstream))  &
              .AND.watelev_sec(j,res_index(upstream)) >= y_sec(m,j,res_index(upstream))) THEN
          topwidth_sec(j,res_index(upstream))=topwidth_sec(j,res_index(upstream))+  &
              ((x_sec(m,j,res_index(upstream))-x_sec(m-1,j,res_index(upstream)))*  &
              (watelev_sec(j,res_index(upstream))-y_sec(m,j,res_index(upstream)))/  &
              (y_sec(m-1,j,res_index(upstream))-y_sec(m,j,res_index(upstream))))
          wetper_sec(j,res_index(upstream))=wetper_sec(j,res_index(upstream))+  &
              SQRT(((watelev_sec(j,res_index(upstream))-y_sec(m,j,res_index(upstream)))/TAN)  &
              **2.+(watelev_sec(j,res_index(upstream))-y_sec(m,j,res_index(upstream)))**2.)
        ELSE IF (watelev_sec(j,res_index(upstream)) >= y_sec(m-1,j,res_index(upstream))  &
              .AND.watelev_sec(j,res_index(upstream)) < y_sec(m,j,res_index(upstream))) THEN
          topwidth_sec(j,res_index(upstream))=topwidth_sec(j,res_index(upstream))+  &
              ((x_sec(m,j,res_index(upstream))-x_sec(m-1,j,res_index(upstream)))*  &
              (watelev_sec(j,res_index(upstream))-y_sec(m-1,j,res_index(upstream)))/  &
              (y_sec(m,j,res_index(upstream))-y_sec(m-1,j,res_index(upstream))))
          wetper_sec(j,res_index(upstream))=wetper_sec(j,res_index(upstream))+  &
              SQRT(((watelev_sec(j,res_index(upstream))-y_sec(m-1,j,res_index(upstream)))/TAN)  &
              **2.+(watelev_sec(j,res_index(upstream))-y_sec(m-1,j,res_index(upstream)))**2.)
        END IF
      END DO
      hydrad_sec(j,res_index(upstream))=area_sec(j,res_index(upstream))/ wetper_sec(j,res_index(upstream))
      meanvel_sec(j,res_index(upstream))=discharge_sec(j,res_index(upstream))/ area_sec(j,res_index(upstream))
      dynhead_sec(j,res_index(upstream))=(meanvel_sec(j,res_index(upstream))**2.) /(2.*9.807)
      tothead_sec(j,res_index(upstream))=watelev_sec(j,res_index(upstream))+ dynhead_sec(j,res_index(upstream))
      energslope_sec(j,res_index(upstream))=(meanvel_sec(j,res_index(upstream))**2.)*  &
          (manning_sec(j,res_index(upstream))**2.)/ (hydrad_sec(j,res_index(upstream))**(4./3.))
      calctothead_sec(j,res_index(upstream))=tothead_sec(j,res_index(upstream))
  END IF

!if (t.eq.1986.and.d.eq.8.and.j.eq.24) then
!write(*,*)t,d,j,depth_sec(j,upstream), &
!     watelev_sec(j,upstream),area_sec(j,upstream), &
!     topwidth_sec(j,upstream),wetper_sec(j,upstream), &
!	 hydrad_sec(j,upstream),discharge_sec(j,upstream),meanvel_sec(j,upstream),&
!	 dynhead_sec(j,upstream),tothead_sec(j,upstream)
!endif
!write(*,*)j,topwidth_sec(j,upstream)
END DO

! 3.2) routing computation in the river subreach

dummy2=0
DO p=1,k
  if (bedslope_sec(k+1-p,res_index(upstream)) <= crslope_sec(k+1-p,res_index(upstream))) then !M1
!    dummy2=dummy2+1
  endif
ENDDO
!flow_regime=real(dummy2)/real(k)
flow_regime=1.

IF (k /= 0) THEN


! 3.2a) for the case M1 (mild slope at the river subreach)
 if (flow_regime>=0.5) then
! for the case, which the reservoir is not completely empty (simulation affected by the reservoir subreach)

  IF (k == nbrsec(upstream)) THEN
    depth_sec(k,res_index(upstream))=crdepth_sec(k,res_index(upstream))
	watelev_sec(k,res_index(upstream))=crwatelev_sec(k,res_index(upstream))
	topwidth_sec(k,res_index(upstream))=crtopwidth_sec(k,res_index(upstream))
    wetper_sec(k,res_index(upstream))=crwetper_sec(k,res_index(upstream))
    area_sec(k,res_index(upstream))=crarea_sec(k,res_index(upstream))
	hydrad_sec(k,res_index(upstream))=crarea_sec(k,res_index(upstream))/ crwetper_sec(k,res_index(upstream))
    meanvel_sec(k,res_index(upstream))=discharge_sec(k,res_index(upstream))/ crarea_sec(k,res_index(upstream))
    dynhead_sec(k,res_index(upstream))=(meanvel_sec(k,res_index(upstream))**2.) /(2.*9.807)
    tothead_sec(k,res_index(upstream))=crwatelev_sec(k,res_index(upstream))+ dynhead_sec(k,res_index(upstream))
    energslope_sec(k,res_index(upstream))=(meanvel_sec(k,res_index(upstream))**2.)*  &
		(manning_sec(k,res_index(upstream))**2.)/ (hydrad_sec(k,res_index(upstream))**(4./3.))
    calctothead_sec(k,res_index(upstream))=tothead_sec(k,upstream)
	k=k-1
  endif

!write(*,*)'Step (3.2a)'
  DO p=1,k
    npt=npoints(k+1-p,res_index(upstream))

! trial-and-error
    error=1.
    error1=0.
    toplim=0.
    lowlim=0.
    interv=0
    dep=0.

    e=0
    dummy3=0
    IF (discharge_sec(k+1-p,res_index(upstream)) > 0.) THEN
      DO WHILE (ABS(error) > 0.001)
        e=e+1
!dummy(k+1-p)=e

! hydraulic parameters at the cross section j
        area_sec(k+1-p,res_index(upstream))=0.
        topwidth_sec(k+1-p,res_index(upstream))=0.
        wetper_sec(k+1-p,res_index(upstream))=0.


        IF (watelev_sec(k+1-p,res_index(upstream)) > maxelev_sec(k+1-p,res_index(upstream))) THEN
          watelev_sec(k+1-p,res_index(upstream))=maxelev_sec(k+1-p,res_index(upstream))
          depth_sec(k+1-p,res_index(upstream))=maxdepth_sec(k+1-p,res_index(upstream))
        END IF

        DO m=1,npt
		  area_part(m,k+1-p)=0.
		enddo

        DO m=2,npt
          TAN=ABS(y_sec(m,k+1-p,res_index(upstream))-y_sec(m-1,k+1-p,res_index(upstream)))/  &
              (x_sec(m,k+1-p,res_index(upstream))-x_sec(m-1,k+1-p,res_index(upstream)))
          IF (watelev_sec(k+1-p,res_index(upstream)) >= y_sec(m-1,k+1-p,res_index(upstream))  &
                .AND.watelev_sec(k+1-p,res_index(upstream)) >= y_sec(m,k+1-p,res_index(upstream))) THEN
			area_part(m,k+1-p)=(x_sec(m,k+1-p,res_index(upstream))-x_sec(m-1,k+1-p,res_index(upstream)))*  &
                (watelev_sec(k+1-p,res_index(upstream))-  &
                (y_sec(m,k+1-p,res_index(upstream))+y_sec(m-1,k+1-p,res_index(upstream)))/2.)
            area_sec(k+1-p,res_index(upstream))=area_sec(k+1-p,res_index(upstream))+area_part(m,k+1-p)
            topwidth_sec(k+1-p,res_index(upstream))=topwidth_sec(k+1-p,res_index(upstream))+  &
                (x_sec(m,k+1-p,res_index(upstream))-x_sec(m-1,k+1-p,res_index(upstream)))
            wetper_sec(k+1-p,res_index(upstream))=wetper_sec(k+1-p,res_index(upstream))+  &
                SQRT((x_sec(m,k+1-p,res_index(upstream))-x_sec(m-1,k+1-p,res_index(upstream)))**2.  &
                +(y_sec(m,k+1-p,res_index(upstream))-y_sec(m-1,k+1-p,res_index(upstream)))**2.)
          ELSE IF (watelev_sec(k+1-p,res_index(upstream)) <  &
                y_sec(m-1,k+1-p,res_index(upstream)).AND.watelev_sec(k+1-p,res_index(upstream))  &
                 >= y_sec(m,k+1-p,res_index(upstream))) THEN
            area_part(m,k+1-p)=(((watelev_sec(k+1-p,res_index(upstream))-  &
                y_sec(m,k+1-p,res_index(upstream)))**2.)/(2.*TAN))
            area_sec(k+1-p,res_index(upstream))=area_sec(k+1-p,res_index(upstream))+area_part(m,k+1-p)
            topwidth_sec(k+1-p,res_index(upstream))=topwidth_sec(k+1-p,res_index(upstream))+  &
                ((x_sec(m,k+1-p,res_index(upstream))-x_sec(m-1,k+1-p,res_index(upstream)))*  &
                (watelev_sec(k+1-p,res_index(upstream))-y_sec(m,k+1-p,res_index(upstream)))/  &
                (y_sec(m-1,k+1-p,res_index(upstream))-y_sec(m,k+1-p,res_index(upstream))))
            wetper_sec(k+1-p,res_index(upstream))=wetper_sec(k+1-p,res_index(upstream))+  &
                SQRT(((watelev_sec(k+1-p,res_index(upstream))-y_sec(m,k+1-p,res_index(upstream)))  &
                /TAN)**2.+(watelev_sec(k+1-p,res_index(upstream))-  &
                y_sec(m,k+1-p,res_index(upstream)))**2.)
          ELSE IF (watelev_sec(k+1-p,res_index(upstream)) >=  &
                y_sec(m-1,k+1-p,res_index(upstream)).AND.watelev_sec(k+1-p,res_index(upstream))  &
                 < y_sec(m,k+1-p,res_index(upstream))) THEN
            area_part(m,k+1-p)=(((watelev_sec(k+1-p,res_index(upstream))-  &
                y_sec(m-1,k+1-p,res_index(upstream)))**2.)/(2.*TAN))
            area_sec(k+1-p,res_index(upstream))=area_sec(k+1-p,res_index(upstream))+area_part(m,k+1-p)
            topwidth_sec(k+1-p,res_index(upstream))=topwidth_sec(k+1-p,res_index(upstream))+  &
                ((x_sec(m,k+1-p,res_index(upstream))-x_sec(m-1,k+1-p,res_index(upstream)))*  &
                (watelev_sec(k+1-p,res_index(upstream))-y_sec(m-1,k+1-p,res_index(upstream)))/  &
                (y_sec(m,k+1-p,res_index(upstream))-y_sec(m-1,k+1-p,res_index(upstream))))
            wetper_sec(k+1-p,res_index(upstream))=wetper_sec(k+1-p,res_index(upstream))+  &
                SQRT(((watelev_sec(k+1-p,res_index(upstream))-  &
                y_sec(m-1,k+1-p,res_index(upstream)))/TAN)**2.+  &
                (watelev_sec(k+1-p,res_index(upstream))- y_sec(m-1,k+1-p,res_index(upstream)))**2.)
          END IF
        END DO

        if (wetper_sec(k+1-p,res_index(upstream)) == 0. .OR. area_sec(k+1-p,res_index(upstream)) ==0.) then
            !no water in cross section
            hydrad_sec(k+1-p,res_index(upstream))= 1.
            meanvel_sec(k+1-p,res_index(upstream))= 0.
            energslope_sec(k+1-p,res_index(upstream)) = 1e-5 !use very small value as a quick and dirty fix
        else !regular case
            hydrad_sec(k+1-p,res_index(upstream))=area_sec(k+1-p,res_index(upstream))/  &
                wetper_sec(k+1-p,res_index(upstream))
            meanvel_sec(k+1-p,res_index(upstream))=discharge_sec(k+1-p,res_index(upstream))/  &
                area_sec(k+1-p,res_index(upstream))
            energslope_sec(k+1-p,res_index(upstream))=(meanvel_sec(k+1-p,res_index(upstream))  &
                **2.)*(manning_sec(k+1-p,res_index(upstream))**2.)  &
                /(hydrad_sec(k+1-p,res_index(upstream))**(4./3.))
         end if
         dynhead_sec(k+1-p,res_index(upstream))=(meanvel_sec(k+1-p,res_index(upstream))**2.)  &
                /(2.*9.807)
         tothead_sec(k+1-p,res_index(upstream))=watelev_sec(k+1-p,res_index(upstream))+  &
                dynhead_sec(k+1-p,res_index(upstream))

		if (bedslope_sec(k+1-p,res_index(upstream)) < crslope_sec(k+1-p,res_index(upstream))) then
		  if (watelev_sec(k+1-p,res_index(upstream)) >= normalelev_sec(k+1-p,res_index(upstream))) then !M1
		    headloss_sec(k+1-p,res_index(upstream))=dist_sec(k+1-p,res_index(upstream))*((2.*discharge_sec(k+1-p,res_index(upstream))/((discharge_sec(k+1-p,res_index(upstream))/sqrt(energslope_sec(k+1-p,res_index(upstream))))+ &
				(discharge_sec(k+2-p,res_index(upstream))/sqrt(energslope_sec(k+2-p,res_index(upstream))))))**2.)
dummy1=11
		  else if (watelev_sec(k+1-p,res_index(upstream)) < crwatelev_sec(k+1-p,res_index(upstream))) then !M3
			headloss_sec(k+1-p,res_index(upstream))=(1./2.)*dist_sec(k+1-p,res_index(upstream))*  &
				(energslope_sec(k+1-p,res_index(upstream))+ energslope_sec(k+2-p,res_index(upstream)))
dummy1=13
		  else !M2
			headloss_sec(k+1-p,res_index(upstream))=dist_sec(k+1-p,res_index(upstream))*sqrt(energslope_sec(k+1-p,res_index(upstream))*energslope_sec(k+2-p,res_index(upstream)))
dummy1=12
          endif
		else
		  if (watelev_sec(k+1-p,res_index(upstream)) >= crwatelev_sec(k+1-p,res_index(upstream))) then !S1
			headloss_sec(k+1-p,res_index(upstream))=dist_sec(k+1-p,res_index(upstream))*sqrt(energslope_sec(k+1-p,res_index(upstream))*energslope_sec(k+2-p,res_index(upstream)))
dummy1=21
		  else if (watelev_sec(k+1-p,res_index(upstream)) < normalelev_sec(k+1-p,res_index(upstream))) then !S3
			headloss_sec(k+1-p,res_index(upstream))=(1./2.)*dist_sec(k+1-p,res_index(upstream))*  &
				(energslope_sec(k+1-p,res_index(upstream))+ energslope_sec(k+2-p,res_index(upstream)))
dummy1=23
		  else !S2
		    headloss_sec(k+1-p,res_index(upstream))=dist_sec(k+1-p,res_index(upstream))*((2.*discharge_sec(k+1-p,res_index(upstream))/((discharge_sec(k+1-p,res_index(upstream))/sqrt(energslope_sec(k+1-p,res_index(upstream))))+ &
				(discharge_sec(k+2-p,res_index(upstream))/sqrt(energslope_sec(k+2-p,res_index(upstream))))))**2.)
dummy1=22
          endif
		endif

! local head-loss coefficient (dimensionless)
        IF (area_sec(k+1-p,res_index(upstream)) >= 1.2*area_sec(k+2-p,res_index(upstream))) THEN
          loclosscoef=0.1
        ELSEIF (area_sec(k+1-p,res_index(upstream)) <= 0.8*area_sec(k+2-p,res_index(upstream))) THEN
		  loclosscoef=0.3
        ELSE
          loclosscoef=0.
        END IF

        if ( area_sec(k+1-p,res_index(upstream)) == 0. .or. area_sec(k+2-p,res_index(upstream)) == 0.) then
            locloss_sec(k+1-p,res_index(upstream)) = 0. !Till: quick fix for dry section
        else
            !regular case
        locloss_sec(k+1-p,res_index(upstream))=(loclosscoef/(2.*9.807))*  &
            (discharge_sec(k+1-p,res_index(upstream))**2.)*  &
            ABS((1./(area_sec(k+1-p,res_index(upstream))**2.))  &
            -(1./(area_sec(k+2-p,res_index(upstream))**2.)))
        end if

        coef=1.

        calctothead_sec(k+1-p,res_index(upstream))=tothead_sec(k+2-p,res_index(upstream))+coef*(headloss_sec(k+1-p,res_index(upstream))  &
            +locloss_sec(k+1-p,res_index(upstream)))

        error=tothead_sec(k+1-p,res_index(upstream))-calctothead_sec(k+1-p,res_index(upstream))

		if (e>20 .and. ABS(error1)>0.1) exit
		if (e>30) exit
        IF (watelev_sec(k+1-p,res_index(upstream)) == maxelev_sec(k+1-p,res_index(upstream))  &
              .AND.calctothead_sec(k+1-p,res_index(upstream)) <  &
              tothead_sec(k+1-p,res_index(upstream)) &
			  .and. lowlim == 0. .and. toplim == 0.) exit
        IF (watelev_sec(k+1-p,res_index(upstream)) == maxelev_sec(k+1-p,res_index(upstream))  &
              .AND.calctothead_sec(k+1-p,res_index(upstream)) <  &
              tothead_sec(k+1-p,res_index(upstream))) THEN
          WRITE(*,*)'ERROR: the water discharge overflows the cross section ', k+1-p
          STOP
        END IF

        if (dummy3==1) then
	      error=0.00001
          watelev_sec(k+1-p,res_index(upstream))=minelev_sec(k+1-p,res_index(upstream))+  &
              depth_sec(k+1-p,res_index(upstream))
          exit
	    endif

        IF (ABS(error) > 0.001) THEN
          IF (interv == 0) THEN
            IF (error1 >= 0. .AND. error > 0.) THEN
              dep=depth_sec(k+1-p,res_index(upstream))
              depth_sec(k+1-p,res_index(upstream))=.80*depth_sec(k+1-p,res_index(upstream))
            ELSE IF (error1 <= 0.AND.error < 0) THEN
              dep=depth_sec(k+1-p,res_index(upstream))
              depth_sec(k+1-p,res_index(upstream))=1.20*depth_sec(k+1-p,res_index(upstream))
            ELSE
              interv=1
            END IF
          END IF

          IF (interv /= 0) THEN
            IF (error1 < 0.AND.error < 0) THEN
              lowlim=MAX(depth_sec(k+1-p,res_index(upstream)),dep)
            ELSE IF (error1 > 0.AND.error > 0) THEN
              toplim=MIN(depth_sec(k+1-p,res_index(upstream)),dep)
            ELSE
              lowlim=MIN(depth_sec(k+1-p,res_index(upstream)),dep)
              toplim=MAX(depth_sec(k+1-p,res_index(upstream)),dep)
            END IF
            dep=depth_sec(k+1-p,res_index(upstream))
            depth_sec(k+1-p,res_index(upstream))=(lowlim+toplim)/2.
          END IF

          error1=error

		  IF (toplim-lowlim<1e-5 .and. toplim/=0. .and. lowlim/=0.) then
		    dep=depth_sec(k+1-p,res_index(upstream))
            depth_sec(k+1-p,res_index(upstream))=(lowlim+toplim)/2.
            dummy3=1
!			write(*,*)'warning! Trial-and-error calculation at the reservoir routing'
!			write(*,*)'stopped for an approximation error of ',error1
!		    write(*,*)'hidraulic calculation (step 3.2a) and cross section',k+1-p
		  ENDIF

          watelev_sec(k+1-p,res_index(upstream))=minelev_sec(k+1-p,res_index(upstream))+  &
              depth_sec(k+1-p,res_index(upstream))

        END IF

!if (error==0.00001) then
!write(*,'(2I4,9F10.5,I4,3F10.5)')k+1-p,dummy1,depth_sec(k+1-p,upstream),dep, &
!     watelev_sec(k+1-p,upstream),area_sec(k+1-p,upstream), &
!     calctothead_sec(k+1-p,upstream),tothead_sec(k+1-p,upstream), &
!     error,lowlim,toplim,interv,headloss_sec(k+1-p,upstream), &
!     locloss_sec(k+1-p,upstream),energslope_sec(k+1-p,upstream)
!endif

      END DO
	  IF (depth_sec(k+1-p,res_index(upstream))<crdepth_sec(k+1-p,res_index(upstream)) .or. ABS(error1)>0.1 .or. e>30) then
	    depth_sec(k+1-p,res_index(upstream))=crdepth_sec(k+1-p,res_index(upstream))
		watelev_sec(k+1-p,res_index(upstream))=crwatelev_sec(k+1-p,res_index(upstream))
	    topwidth_sec(k+1-p,res_index(upstream))=crtopwidth_sec(k+1-p,res_index(upstream))
        wetper_sec(k+1-p,res_index(upstream))=crwetper_sec(k+1-p,res_index(upstream))
        area_sec(k+1-p,res_index(upstream))=crarea_sec(k+1-p,res_index(upstream))
        hydrad_sec(k+1-p,res_index(upstream))=crarea_sec(k+1-p,res_index(upstream))/  &
				crwetper_sec(k+1-p,res_index(upstream))
		meanvel_sec(k+1-p,res_index(upstream))=discharge_sec(k+1-p,res_index(upstream))/  &
				crarea_sec(k+1-p,res_index(upstream))
		dynhead_sec(k+1-p,res_index(upstream))=(meanvel_sec(k+1-p,res_index(upstream))**2.)  &
				/(2.*9.807)
		tothead_sec(k+1-p,res_index(upstream))=watelev_sec(k+1-p,res_index(upstream))+  &
				dynhead_sec(k+1-p,res_index(upstream))
		energslope_sec(k+1-p,res_index(upstream))=(meanvel_sec(k+1-p,res_index(upstream))  &
				**2.)*(manning_sec(k+1-p,res_index(upstream))**2.)  &
				/(hydrad_sec(k+1-p,res_index(upstream))**(4./3.))
!write(*,*)'warning! forced critical flow regime at cross section',k+1-p
	  endif

      if(watelev_sec(k+1-p,res_index(upstream))>maxelev_sec(k+1-p,res_index(upstream))) then
        WRITE(*,*)'ERROR: the water discharge overflows the cross section ',k+1-p
        STOP
      END IF

    END IF


  END DO
!  stop

! 3.2b) for the case S1 (steep slope at the river subreach)
 else
! for the case, which the reservoir is completely empty (simulation is not affected by the reservoir subreach)
  depth_sec(1,res_index(upstream))=crdepth_sec(1,res_index(upstream))
  watelev_sec(1,res_index(upstream))=crwatelev_sec(1,res_index(upstream))
  topwidth_sec(1,res_index(upstream))=crtopwidth_sec(1,res_index(upstream))
  wetper_sec(1,res_index(upstream))=crwetper_sec(1,res_index(upstream))
  area_sec(1,res_index(upstream))=crarea_sec(1,res_index(upstream))
  hydrad_sec(1,res_index(upstream))=crarea_sec(1,res_index(upstream))/ crwetper_sec(1,res_index(upstream))
  meanvel_sec(1,res_index(upstream))=discharge_sec(1,res_index(upstream))/ crarea_sec(1,res_index(upstream))
  dynhead_sec(1,res_index(upstream))=(meanvel_sec(1,res_index(upstream))**2.) /(2.*9.807)
  tothead_sec(1,res_index(upstream))=crwatelev_sec(1,res_index(upstream))+ dynhead_sec(1,res_index(upstream))
  energslope_sec(1,res_index(upstream))=(meanvel_sec(1,res_index(upstream))**2.)*  &
		(manning_sec(1,res_index(upstream))**2.)/ (hydrad_sec(1,res_index(upstream))**(4./3.))
  calctothead_sec(1,res_index(upstream))=tothead_sec(1,upstream)

  if (k>1) then
!write(*,*)'Step (3.2b)'
   DO p=2,k
    npt=npoints(p,res_index(upstream))

! trial-and-error
    error=1.
    error1=0.
    toplim=0.
    lowlim=0.
    interv=0
    dep=0.

    f=0
    dummy3=0
    IF (discharge_sec(p,res_index(upstream)) > 0.) THEN
      DO WHILE (ABS(error) > 0.001)
        f=f+1
!dummy(p)=e

! hydraulic parameters at the cross section j
        area_sec(p,res_index(upstream))=0.
        topwidth_sec(p,res_index(upstream))=0.
        wetper_sec(p,res_index(upstream))=0.


        IF (watelev_sec(p,res_index(upstream)) > maxelev_sec(p,res_index(upstream))) THEN
          watelev_sec(p,res_index(upstream))=maxelev_sec(p,res_index(upstream))
          depth_sec(p,res_index(upstream))=maxdepth_sec(p,res_index(upstream))
        END IF

        DO m=2,npt
          TAN=ABS(y_sec(m,p,res_index(upstream))-y_sec(m-1,p,res_index(upstream)))/  &
              (x_sec(m,p,res_index(upstream))-x_sec(m-1,p,res_index(upstream)))
          IF (watelev_sec(p,res_index(upstream)) >= y_sec(m-1,p,res_index(upstream))  &
                .AND.watelev_sec(p,res_index(upstream)) >= y_sec(m,p,res_index(upstream))) THEN
            area_part(m,p)=(x_sec(m,p,res_index(upstream))-x_sec(m-1,p,res_index(upstream)))*  &
                (watelev_sec(p,res_index(upstream))-  &
                (y_sec(m,p,res_index(upstream))+y_sec(m-1,p,res_index(upstream)))/2.)
			area_sec(p,res_index(upstream))=area_sec(p,res_index(upstream))+area_part(m,p)
            topwidth_sec(p,res_index(upstream))=topwidth_sec(p,res_index(upstream))+  &
                (x_sec(m,p,res_index(upstream))-x_sec(m-1,p,res_index(upstream)))
            wetper_sec(p,res_index(upstream))=wetper_sec(p,res_index(upstream))+  &
                SQRT((x_sec(m,p,res_index(upstream))-x_sec(m-1,p,res_index(upstream)))**2.  &
                +(y_sec(m,p,res_index(upstream))-y_sec(m-1,p,res_index(upstream)))**2.)
          ELSE IF (watelev_sec(p,res_index(upstream)) <  &
                y_sec(m-1,p,res_index(upstream)).AND.watelev_sec(p,res_index(upstream))  &
                 >= y_sec(m,p,res_index(upstream))) THEN
            area_part(m,p)=(((watelev_sec(p,res_index(upstream))-  &
                y_sec(m,p,res_index(upstream)))**2.)/(2.*TAN))
			area_sec(p,res_index(upstream))=area_sec(p,res_index(upstream))+area_part(m,p)
            topwidth_sec(p,res_index(upstream))=topwidth_sec(p,res_index(upstream))+  &
                ((x_sec(m,p,res_index(upstream))-x_sec(m-1,p,res_index(upstream)))*  &
                (watelev_sec(p,res_index(upstream))-y_sec(m,p,res_index(upstream)))/  &
                (y_sec(m-1,p,res_index(upstream))-y_sec(m,p,res_index(upstream))))
            wetper_sec(p,res_index(upstream))=wetper_sec(p,res_index(upstream))+  &
                SQRT(((watelev_sec(p,res_index(upstream))-y_sec(m,p,res_index(upstream)))  &
                /TAN)**2.+(watelev_sec(p,res_index(upstream))-  &
                y_sec(m,p,res_index(upstream)))**2.)
          ELSE IF (watelev_sec(p,res_index(upstream)) >=  &
                y_sec(m-1,p,res_index(upstream)).AND.watelev_sec(p,res_index(upstream))  &
                 < y_sec(m,p,res_index(upstream))) THEN
            area_part(m,p)=(((watelev_sec(p,res_index(upstream))-  &
                y_sec(m-1,p,res_index(upstream)))**2.)/(2.*TAN))
			area_sec(p,res_index(upstream))=area_sec(p,res_index(upstream))+area_part(m,p)
            topwidth_sec(p,res_index(upstream))=topwidth_sec(p,res_index(upstream))+  &
                ((x_sec(m,p,res_index(upstream))-x_sec(m-1,p,res_index(upstream)))*  &
                (watelev_sec(p,res_index(upstream))-y_sec(m-1,p,res_index(upstream)))/  &
                (y_sec(m,p,res_index(upstream))-y_sec(m-1,p,res_index(upstream))))
            wetper_sec(p,res_index(upstream))=wetper_sec(p,res_index(upstream))+  &
                SQRT(((watelev_sec(p,res_index(upstream))-  &
                y_sec(m-1,p,res_index(upstream)))/TAN)**2.+  &
                (watelev_sec(p,res_index(upstream))- y_sec(m-1,p,res_index(upstream)))**2.)
          END IF
        END DO

        hydrad_sec(p,res_index(upstream))=area_sec(p,res_index(upstream))/  &
            wetper_sec(p,res_index(upstream))
        meanvel_sec(p,res_index(upstream))=discharge_sec(p,res_index(upstream))/  &
            area_sec(p,res_index(upstream))
        dynhead_sec(p,res_index(upstream))=(meanvel_sec(p,res_index(upstream))**2.)  &
            /(2.*9.807)
        tothead_sec(p,res_index(upstream))=watelev_sec(p,res_index(upstream))+  &
            dynhead_sec(p,res_index(upstream))
        energslope_sec(p,res_index(upstream))=(meanvel_sec(p,res_index(upstream))  &
            **2.)*(manning_sec(p,res_index(upstream))**2.)  &
            /(hydrad_sec(p,res_index(upstream))**(4./3.))

		if (bedslope_sec(p,res_index(upstream)) < crslope_sec(p,res_index(upstream))) then
		  if (watelev_sec(p,res_index(upstream)) >= normalelev_sec(p,res_index(upstream))) then !M1
		    headloss_sec(p,res_index(upstream))=dist_sec(p,res_index(upstream))*((2.*discharge_sec(p,res_index(upstream))/((discharge_sec(p,res_index(upstream))/sqrt(energslope_sec(p,res_index(upstream))))+ &
				(discharge_sec(p-1,res_index(upstream))/sqrt(energslope_sec(p-1,res_index(upstream))))))**2.)
dummy1=11
		  else if (watelev_sec(p,res_index(upstream)) < crwatelev_sec(p,res_index(upstream))) then !M3
			headloss_sec(p,res_index(upstream))=(1./2.)*dist_sec(p,res_index(upstream))*  &
				(energslope_sec(p,res_index(upstream))+ energslope_sec(p-1,res_index(upstream)))
dummy1=13
		  else !M2
			headloss_sec(p,res_index(upstream))=dist_sec(p,res_index(upstream))*sqrt(energslope_sec(p,res_index(upstream))*energslope_sec(p-1,res_index(upstream)))
dummy1=12
          endif
		else
		  if (watelev_sec(p,res_index(upstream)) >= crwatelev_sec(p,res_index(upstream))) then !S1
			headloss_sec(p,res_index(upstream))=dist_sec(p,res_index(upstream))*sqrt(energslope_sec(p,res_index(upstream))*energslope_sec(p-1,res_index(upstream)))
dummy1=21
		  else if (watelev_sec(p,res_index(upstream)) < normalelev_sec(p,res_index(upstream))) then !S3
			headloss_sec(p,res_index(upstream))=(1./2.)*dist_sec(p,res_index(upstream))*  &
				(energslope_sec(p,res_index(upstream))+ energslope_sec(p-1,res_index(upstream)))
dummy1=23
		  else !S2
		    headloss_sec(p,res_index(upstream))=dist_sec(p,res_index(upstream))*((2.*discharge_sec(p,res_index(upstream))/((discharge_sec(p,res_index(upstream))/sqrt(energslope_sec(p,res_index(upstream))))+ &
				(discharge_sec(p-1,res_index(upstream))/sqrt(energslope_sec(p-1,res_index(upstream))))))**2.)
dummy1=22
          endif
		endif

! local head-loss coefficient (dimensionless)
        IF (area_sec(p,res_index(upstream)) >= 1.2*area_sec(p-1,res_index(upstream))) THEN
          loclosscoef=0.05
        ELSEIF (area_sec(p,res_index(upstream)) <= 0.8*area_sec(p-1,res_index(upstream))) THEN
		  loclosscoef=0.1
        ELSE
          loclosscoef=0.
        END IF

        locloss_sec(p,res_index(upstream))=(loclosscoef/(2.*9.807))*  &
            (discharge_sec(p,res_index(upstream))**2.)*  &
            ABS((1./(area_sec(p,res_index(upstream))**2.))  &
            -(1./(area_sec(p-1,res_index(upstream))**2.)))

        coef=-1.

        calctothead_sec(p,res_index(upstream))=tothead_sec(p-1,upstream)+coef*(headloss_sec(p,res_index(upstream))  &
            +locloss_sec(p,res_index(upstream)))

        error=tothead_sec(p,res_index(upstream))-calctothead_sec(p,res_index(upstream))

		if (e>20 .and. ABS(error1)>0.1) exit
		if (e>30) exit
        IF (watelev_sec(p,res_index(upstream)) == maxelev_sec(p,res_index(upstream))  &
              .AND.calctothead_sec(p,res_index(upstream)) <  &
              tothead_sec(p,res_index(upstream)) &
			  .and. lowlim == 0. .and. toplim == 0.) exit
        IF (watelev_sec(p,res_index(upstream)) == maxelev_sec(p,res_index(upstream))  &
              .AND.calctothead_sec(p,res_index(upstream)) <  &
              tothead_sec(p,res_index(upstream))) THEN
          WRITE(*,*)'ERROR: the water discharge overflows the cross section ', p
          STOP
        END IF

        if (dummy3==1) then
	      error=0.00001
          watelev_sec(p,res_index(upstream))=minelev_sec(p,res_index(upstream))+  &
              depth_sec(p,res_index(upstream))
          exit
	    endif


        IF (ABS(error) > 0.001) THEN
          IF (interv == 0) THEN
            IF (error1 >= 0. .AND. error > 0.) THEN
              dep=depth_sec(p,res_index(upstream))
              depth_sec(p,res_index(upstream))=.80*depth_sec(p,res_index(upstream))
            ELSE IF (error1 <= 0.AND.error < 0) THEN
              dep=depth_sec(p,res_index(upstream))
              depth_sec(p,res_index(upstream))=1.20*depth_sec(p,res_index(upstream))
            ELSE
              interv=1
            END IF
          END IF

          IF (interv /= 0) THEN
            IF (error1 < 0.AND.error < 0) THEN
              lowlim=MAX(depth_sec(p,res_index(upstream)),dep)
            ELSE IF (error1 > 0.AND.error > 0) THEN
              toplim=MIN(depth_sec(p,res_index(upstream)),dep)
            ELSE
              lowlim=MIN(depth_sec(p,res_index(upstream)),dep)
              toplim=MAX(depth_sec(p,res_index(upstream)),dep)
            END IF
            dep=depth_sec(p,res_index(upstream))
            depth_sec(p,res_index(upstream))=(lowlim+toplim)/2.
          END IF

          error1=error

		  IF (toplim-lowlim <1e-5 .and. toplim/=0. .and. lowlim/=0.) then
		    dep=depth_sec(p,res_index(upstream))
            depth_sec(p,res_index(upstream))=(lowlim+toplim)/2.
            dummy3=1
!			write(*,*)'warning! Trial-and-error calculation at the reservoir routing'
!			write(*,*)'stopped for an approximation error of ',error1
!		    write(*,*)'hidraulic calculation (step 3.2a) and cross section',p
		  ENDIF

          watelev_sec(p,res_index(upstream))=minelev_sec(p,res_index(upstream))+  &
              depth_sec(p,res_index(upstream))
        END IF

!if (error==0.00001) then
!write(*,'(2I4,9F10.5,I4,2F10.5)')p,dummy1,depth_sec(p,upstream),dep, &
!     watelev_sec(p,upstream),area_sec(p,upstream), &
!     calctothead_sec(p,upstream),tothead_sec(p,upstream), &
!     error,lowlim,toplim,interv,headloss_sec(p,upstream), &
!     locloss_sec(p,upstream)
!endif

      END DO

	  IF (depth_sec(p,res_index(upstream))<crdepth_sec(p,res_index(upstream)) .or. ABS(error1)>0.1 .or. e>30) then
	    depth_sec(p,res_index(upstream))=crdepth_sec(p,res_index(upstream))
		watelev_sec(p,res_index(upstream))=crwatelev_sec(p,res_index(upstream))
	    topwidth_sec(p,res_index(upstream))=crtopwidth_sec(p,res_index(upstream))
        wetper_sec(p,res_index(upstream))=crwetper_sec(p,res_index(upstream))
        area_sec(p,res_index(upstream))=crarea_sec(p,res_index(upstream))
        hydrad_sec(p,res_index(upstream))=crarea_sec(p,res_index(upstream))/  &
				crwetper_sec(p,res_index(upstream))
		meanvel_sec(p,res_index(upstream))=discharge_sec(p,res_index(upstream))/  &
				crarea_sec(p,res_index(upstream))
		dynhead_sec(p,res_index(upstream))=(meanvel_sec(p,res_index(upstream))**2.)  &
				/(2.*9.807)
		tothead_sec(p,res_index(upstream))=watelev_sec(p,res_index(upstream))+  &
				dynhead_sec(p,res_index(upstream))
		energslope_sec(p,res_index(upstream))=(meanvel_sec(p,res_index(upstream))  &
				**2.)*(manning_sec(p,res_index(upstream))**2.)  &
				/(hydrad_sec(p,res_index(upstream))**(4./3.))
!write(*,*)'warning! forced critical flow regime at the cross section',p
	  endif

      if(watelev_sec(p,res_index(upstream))>maxelev_sec(p,res_index(upstream))) then
        WRITE(*,*)'ERROR: the water discharge overflows the cross section ',p
        STOP
      END IF

    END IF


   END DO

  endif

 endif


END IF

!DO j=1,nbrsec(upstream)
!write(*,'(2I4,3F10.6)')j,k,crdepth_sec(j,upstream),normaldepth(j),depth_sec(j,upstream)
!ENDDO

END SUBROUTINE hydraul_res
