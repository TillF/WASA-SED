SUBROUTINE change_sec (j1,upstream)
! Code converted using TO_F90 by Alan Miller
! Date: 2005-08-23  Time: 12:57:15
 
use common_h
use time_h
use reservoir_h

IMPLICIT NONE

INTEGER, INTENT(IN OUT)                  :: upstream

INTEGER :: j1,g !,npt

!REAL :: max_dheight,min_dheight,dummy,dist_bank1,dist_bank2
!REAL :: max_adjust
REAL ::elev
!REAL :: predom_grain,perc_grain

!REAL :: dummy1,dummy2
real :: dummy3,dummy4 !,damelev, dummy,xmin_tw,
real :: wfactor(200),dummy5
! -----------------------------------------------------------------------

! calculation of reservoir bed elevation changes

! identification of the sediment balance process (erosion/deposition)

    
! 1) erosion
!dummy1=0

IF (dvol_sed(j1,res_index(upstream)) < 0.) THEN
  elev=erosion_level(j1,res_index(upstream))
!  elev=maxelev_sec(j1,upstream)
!  elev=watelev_sec(j1,upstream)
!  elev=min(minelev_sec(j1,upstream)+(depth_sec(j1,upstream)*1.5),maxelev_sec(j1,upstream))
!  elev=watelev_sec(j1,upstream)+((maxelev_sec(j1,upstream)-watelev_sec(j1,upstream))/20.)
!if (id_sec_extern(j1,upstream)==43) elev=max(elev,y_sec(11,j1,upstream))      
!if (id_sec_extern(j1,upstream)==47) elev=max(elev,y_sec(18,j1,upstream))      
!if (id_sec_extern(j1,upstream)==47)elev=max(elev,435.)
!  elev=min(watelev_sec(j1,upstream)+5.,maxelev_sec(j1,upstream))
!  elev=min(watelev_sec(j1,upstream)+1.5,watelev_sec(j1,upstream)+((maxelev_sec(j1,upstream)-watelev_sec(j1,upstream))/5.))
  dummy5=0.
  DO m=1,npoints(j1,res_index(upstream))
    if(y_sec(m,j1,res_index(upstream))<elev) then
	  dummy5=dummy5+(y_sec(m,j1,res_index(upstream))-y_actlay(m,j1,res_index(upstream)))
	endif
  enddo

  DO m=1,npoints(j1,res_index(upstream))
    if(y_sec(m,j1,res_index(upstream))<elev .and. dummy5/=0.) then
	  wfactor(m)=(y_sec(m,j1,res_index(upstream))-y_actlay(m,j1,res_index(upstream)))/dummy5
	else
	  wfactor(m)=0.
    endif
  enddo

  dummy3=0.
  DO m=2,npoints(j1,res_index(upstream))-1
    IF (geom(m,j1) == 1)dummy3=dummy3+(((x_sec(m,j1,res_index(upstream))-x_sec(m-1,j1,res_index(upstream)))*wfactor(m)*(1./2.))+ &
							((x_sec(m+1,j1,res_index(upstream))-x_sec(m,j1,res_index(upstream)))*wfactor(m)*(3./8.))+ &
							((x_sec(m+1,j1,res_index(upstream))-x_sec(m,j1,res_index(upstream)))*wfactor(m+1)*(1./8.)))
    IF (geom(m,j1) == 2)dummy3=dummy3+((x_sec(m,j1,res_index(upstream))-x_sec(m-1,j1,res_index(upstream)))*wfactor(m)*(3./8.))+ &
							((x_sec(m,j1,res_index(upstream))-x_sec(m-1,j1,res_index(upstream)))*wfactor(m-1)*(1./8.))+ &
							((x_sec(m+1,j1,res_index(upstream))-x_sec(m,j1,res_index(upstream)))*wfactor(m)*(3./8.))+ &
							((x_sec(m+1,j1,res_index(upstream))-x_sec(m,j1,res_index(upstream)))*wfactor(m+1)*(1./8.))
    IF (geom(m,j1) == 3)dummy3=dummy3+((x_sec(m,j1,res_index(upstream))-x_sec(m-1,j1,res_index(upstream)))*wfactor(m)*(3./8.))+ &
							((x_sec(m,j1,res_index(upstream))-x_sec(m-1,j1,res_index(upstream)))*wfactor(m-1)*(1./8.))+ &
							((x_sec(m+1,j1,res_index(upstream))-x_sec(m,j1,res_index(upstream)))*wfactor(m)*(1./2.))
    IF (geom(m,j1) == 4)dummy3=dummy3+((x_sec(m+1,j1,res_index(upstream))-x_sec(m-1,j1,res_index(upstream)))*wfactor(m)*(1./2.))
!if(j1<30)write(*,'(3I4,4F10.2)')j1,m,geom(m,j1),dummy3
  enddo
  
  dummy4=darea_sed(j1,res_index(upstream))/dummy3
  DO m=2,npoints(j1,res_index(upstream))-1
	y_laststep(m,j1,res_index(upstream))=y_sec(m,j1,res_index(upstream))
	IF (geom(m,j1) /= 0)y_sec(m,j1,res_index(upstream))=y_sec(m,j1,res_index(upstream))-dummy4*wfactor(m)	
!if(j1<20)write(*,'(2I4,3I3,6F12.6)')t,d,j1,m,geom(m,j1),y_sec(m,j1,upstream),y_laststep(m,j1,upstream),dummy4*wfactor(m),dummy4,wfactor(m)
!write(*,'(2I4,3I3,6F12.6)')t,d,j1,m,geom(m,j1),y_sec(m,j1,upstream),y_laststep(m,j1,upstream),dummy4*wfactor(m),dummy4,wfactor(m)
  ENDDO

!if(d==2)stop

! 2) deposition
ELSE IF (dvol_sed(j1,res_index(upstream)) > 0.) THEN

  DO m=1,npoints(j1,res_index(upstream))
    if(y_sec(m,j1,res_index(upstream))<watelev_sec(j1,res_index(upstream))) then
!	  dummy5=(watelev_sec(j1,upstream)-y_sec(m,j1,upstream))/(watelev_sec(j1,upstream)-minelev_sec(j1,upstream))
!	  wfactor(m)=1.-((1.-(dummy5))**2.9)
	  wfactor(m)=(watelev_sec(j1,res_index(upstream))-y_sec(m,j1,res_index(upstream)))/(watelev_sec(j1,res_index(upstream))-minelev_sec(j1,res_index(upstream)))
	else
!	  dummy5=0.
	  wfactor(m)=0.
    endif
!if(j1<20)write(*,'(2I4,3I3,2F15.6)')t,d,j1,m,geom(m,j1),wfactor(m),dummy5
  enddo
    
  dummy3=0.
  DO m=2,npoints(j1,res_index(upstream))-1
    IF (geom(m,j1) == 1)dummy3=dummy3+(((x_sec(m,j1,res_index(upstream))-x_sec(m-1,j1,res_index(upstream)))*wfactor(m)*(1./2.))+ &
							((x_sec(m+1,j1,res_index(upstream))-x_sec(m,j1,res_index(upstream)))*wfactor(m)*(3./8.))+ &
							((x_sec(m+1,j1,res_index(upstream))-x_sec(m,j1,res_index(upstream)))*wfactor(m+1)*(1./8.)))
    IF (geom(m,j1) == 2)dummy3=dummy3+((x_sec(m,j1,res_index(upstream))-x_sec(m-1,j1,res_index(upstream)))*wfactor(m)*(3./8.))+ &
							((x_sec(m,j1,res_index(upstream))-x_sec(m-1,j1,res_index(upstream)))*wfactor(m-1)*(1./8.))+ &
							((x_sec(m+1,j1,res_index(upstream))-x_sec(m,j1,res_index(upstream)))*wfactor(m)*(3./8.))+ &
							((x_sec(m+1,j1,res_index(upstream))-x_sec(m,j1,res_index(upstream)))*wfactor(m+1)*(1./8.))
    IF (geom(m,j1) == 3)dummy3=dummy3+((x_sec(m,j1,res_index(upstream))-x_sec(m-1,j1,res_index(upstream)))*wfactor(m)*(3./8.))+ &
							((x_sec(m,j1,res_index(upstream))-x_sec(m-1,j1,res_index(upstream)))*wfactor(m-1)*(1./8.))+ &
							((x_sec(m+1,j1,res_index(upstream))-x_sec(m,j1,res_index(upstream)))*wfactor(m)*(1./2.))
    IF (geom(m,j1) == 4)dummy3=dummy3+((x_sec(m+1,j1,res_index(upstream))-x_sec(m-1,j1,res_index(upstream)))*wfactor(m)*(1./2.))
!if(j1<30)write(*,'(3I4,4F10.2)')j1,m,geom(m,j1),dummy3
  enddo
  
  dummy4=darea_sed(j1,res_index(upstream))/dummy3
  DO m=2,npoints(j1,res_index(upstream))-1
	y_laststep(m,j1,res_index(upstream))=y_sec(m,j1,res_index(upstream))
	IF (geom(m,j1) /= 0)y_sec(m,j1,res_index(upstream))=y_sec(m,j1,res_index(upstream))+dummy4*wfactor(m)	
!if(j1<20)write(*,'(2I4,3I3,6F12.6)')t,d,j1,m,geom(m,j1),y_sec(m,j1,upstream),y_laststep(m,j1,upstream),dummy4*wfactor(m),dummy4,wfactor(m)
  ENDDO

!if(d==2)stop

END IF

!if(d==41)stop




DO g=1,n_sed_class
  if (frac_actlay(g,j1,res_index(upstream)) /= 0.) then
!	dummy1=diam(g)     
  endif
enddo

!Ge User-specified positive multiplication factor for defining the thickness of active layer
!Ge constant=10/20/30
!dummy2=20.*dummy1*24.



!DO m=1,npoints(j1,upstream)
!  if (x_sec(m,j1,upstream) == x_minelev(j1,upstream)) then
!    y_actlay(m,j1,upstream)=max(y_sec(m,j1,upstream)-dummy,y_original(m,j1,upstream))
!	dummy=y_actlay(m,j1,upstream)
!  endif
!enddo
DO m=1,npoints(j1,res_index(upstream))
!  if (x_sec(m,j1,upstream) /= x_minelev(j1,upstream)) then
!    y_actlay(m,j1,upstream)=max(y_sec(m,j1,upstream)-dummy2,y_original(m,j1,upstream))
!  endif
!if(j1==11)write(*,'(2I4,3F15.8)')j1,m,y_original(m,j1,upstream),y_actlay(m,j1,upstream),y_sec(m,j1,upstream)  
ENDDO


IF (dvol_sed(j1,res_index(upstream)) < 0.) THEN

if (j1/=nbrsec(res_index(upstream))) then
!write(*,'(I4,4F10.2,A10)')j1,minelev_sec(j1,upstream),dheight_sed(j1,upstream),minelev_sec(j1,upstream)-dheight_sed(j1,upstream),minelev_sec(j1+1,upstream)+(mean_slope*dist_sec(j1,upstream)),'eros'
else
!write(*,'(I4,4F10.2,A10)')j1,minelev_sec(j1,upstream),dheight_sed(j1,upstream),minelev_sec(j1,upstream)-dheight_sed(j1,upstream),minelev_sec(j1-1,upstream)-(mean_slope*dist_sec(j1-1,upstream)),'eros'
endif

else

if (j1/=nbrsec(res_index(upstream))) then
!write(*,'(I4,4F10.2,A10)')j1,minelev_sec(j1,upstream),dheight_sed(j1,upstream),minelev_sec(j1,upstream)+dheight_sed(j1,upstream),minelev_sec(j1+1,upstream)+(mean_slope*dist_sec(j1,upstream)),'depo'
else
!write(*,'(I4,4F10.2,A10)')j1,minelev_sec(j1,upstream),dheight_sed(j1,upstream),minelev_sec(j1,upstream)+dheight_sed(j1,upstream),minelev_sec(j1-1,upstream)-(mean_slope*dist_sec(j1-1,upstream)),'depo'
endif

endif

DO m=1,npoints(j1,res_index(upstream))
!write(*,*)j1,y_sec(m,j1,upstream),y_laststep(m,j1,upstream)
ENDDO

!if (d==7)stop


RETURN
END SUBROUTINE change_sec
