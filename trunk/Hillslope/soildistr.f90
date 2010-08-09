SUBROUTINE soildistr(theta,out)
 
! Code converted using TO_F90 by Alan Miller
! Date: 2005-06-30  Time: 13:53:36

use common_h
use hymo_h
use params_h

!** subroutine creates array with
!** distribution functions of soil parameters
!**
!** - one distribution for each soil component (SVC)
!** - linear interpolation between points of distribution

IMPLICIT NONE


REAL, INTENT(IN)                         :: theta(sv_comb,maxsoil)
REAL, INTENT(OUT)                        :: out(5,2,sv_comb,maxsoil)

!  theta:   input soil parameters for SVC (Vol%)
!  out:     values of distribution funtion of soil parameter
!           fraction of SVC versus  storage volume (mm)



REAL :: tempx,var1,var2

INTEGER :: k,i


!  variability around given value of SVC (first interval)
!  (reference to soil characteristic in mm)
var1=0.05
!  variability around given value of SVC (second interval)
!  (reference to soil characteristic in mm)
var2=0.1

!** Loop for each soil component
DO k=1,sv_comb
  DO i=1,nbr_svc(k)
    tempx=theta(k,i)
    
    out(1,1,k,i)=tempx-var2*tempx
    out(2,1,k,i)=tempx-var1*tempx
    out(3,1,k,i)=tempx
    out(4,1,k,i)=tempx+var1*tempx
    out(5,1,k,i)=tempx+var2*tempx
    
    out(1,2,k,i)=0.0
    out(2,2,k,i)=0.1
    out(3,2,k,i)=0.5
    out(4,2,k,i)=0.9
    out(5,2,k,i)=1.0
    
!  end of loop for all soil components
  END DO
END DO

RETURN
END SUBROUTINE soildistr

