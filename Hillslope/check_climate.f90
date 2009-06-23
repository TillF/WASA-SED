SUBROUTINE check_climate
!check read climate data for implausible data

!Till: computationally irrelevant: improved error reporting in faulty climate data
!2009-06-22

!Till: fixed faulty check of hourly rainfall
!2009-04-06
 
!Till: added check of hourly rainfall
!2009-03-11
 
!Till: removed erroneous error messages that could occur when starting midyear
!2008-04-11

!Till: modified error messages
!2008-01-21


use climo_h
use common_h
!use params_h
use hymo_h
use time_h

IMPLICIT NONE


INTEGER :: ii,jj
REAL :: store=1

!check temperature data
DO jj=1,min(dayyear,size(temp,DIM=1))
	DO ii=1,size(temp,DIM=2)
		if ((temp(jj,ii)>50) .OR. (temp(jj,ii)<-50)) then
			if (jj/=366) then	!leap year values might be missing, don't alert
				write(*,'(A,i0,a,i0,a,i0,a,f7.1,a)')'Temperature of subbasin ',id_subbas_extern(ii),', day ', jj+dayoutsim,', year ',t,' is implausible (',temp(jj,ii),'). using previous day.'
			end if
			temp(jj,ii)=store	!implausible value, use the previous day
		else
			store=temp(jj,ii)
		end if
	end do
END DO


!check humidity data
DO jj=1,min(dayyear,size(rhum,DIM=1))
	DO ii=1,size(rhum,DIM=2)
		if ((rhum(jj,ii)>100) .OR. (rhum(jj,ii)<0)) then
			if (jj/=366) then	!leap year values might be missing, don't alert
				write(*,'(A,i0,a,i0,a,i0,a,f7.1,a)')'Humidity of subbasin ',id_subbas_extern(ii),', day ', jj+dayoutsim,', year '&
					,t,' is implausible (',rhum(jj,ii),'). using previous day.'
			end if
			rhum(jj,ii)=store	!implausible value, use the previous day
		else
			store=rhum(jj,ii)
		end if
	end do
END DO


!check radiation data
DO jj=1,min(dayyear,size(rad,DIM=1))
	DO ii=1,size(rad,DIM=2)
		if ((rad(jj,ii)>700) .OR. (rad(jj,ii)<1)) then
			if (jj/=366) then	!leap year values might be missing, don't alert
				write(*,'(A,i0,a,i0,a,i0,a,f7.1,a)')'Radiation of subbasin ',id_subbas_extern(ii),', day ', jj+dayoutsim,', year '&
					,t,' is implausible (',rad(jj,ii),'). using previous day.'
			end if
			rad(jj,ii)=store	!implausible value, use the previous day
		else
			store=rad(jj,ii)
		end if
	end do
END DO

!check hourly precip data
if (dohour) then
	DO jj=1,min(dayyear*nt,size(preciph,DIM=1))
		DO ii=1,size(preciph,DIM=2)
			if ((preciph(jj,ii)>50) .OR. (preciph(jj,ii)<0)) then
				if (jj>=365*nt) then	!leap year values might be missing, don't alert
					write(*,'(A,i0,a,i0,a,i0,a,i0,a,f7.1,a)')'Hourly precipitation data of subbasin ',id_subbas_extern(ii),', timestep ', mod(jj,nt),', day ', jj/nt+1,', year '&
						,t,' is implausible(',preciph(jj,ii),'). using previous time step.'
				end if
				preciph(jj,ii)=store	!implausible value, use the previous time step
			else
				store=preciph(jj,ii)	!last plausible value encountered
			end if
		end do
	END DO
end if


!check precip data
DO jj=1,min(dayyear,size(precip,DIM=1))
	DO ii=1,size(precip,DIM=2)
		if ((precip(jj,ii)>50*24) .OR. (precip(jj,ii)<0)) then
			if (jj/=366) then	!leap year values might be missing, don't alert
				write(*,'(A,i0,a,i0,a,i0,a,f7.1,a)')'Precipitation data of subbasin ',id_subbas_extern(ii),', day ', jj,', year '&
					,t,' is implausible(',precip(jj,ii),'). using previous day.'
			end if
			precip(jj,ii)=store	!implausible value, use the previous day
		else
			store=precip(jj,ii)
		end if
	end do
END DO



RETURN
END SUBROUTINE check_climate
