! Till: computationally irrelevant: outcommented unused variables
! 2012-09-14
!
! Till: computationally irrelevant: minor changes to improve compiler compatibility
! 2011-04-29
!
!The module file contains two modules: reservoir_h and lake_h
!

MODULE numint

PRIVATE
PUBLIC qtrap

INTEGER(4),PARAMETER:: NPAR_ARTH= 16, NPAR2_ARTH= 8

CONTAINS

FUNCTION arth_d(first,increment,n)
IMPLICIT NONE
REAL, INTENT(IN) :: first,increment
INTEGER(4), INTENT(IN) :: n
REAL, DIMENSION(n) :: arth_d
INTEGER(4) :: k,k2
REAL :: temp
if (n > 0) arth_d(1)=first
if (n <= NPAR_ARTH) then
  do k=2,n
    arth_d(k)=arth_d(k-1)+increment
  end do
else
  do k=2,NPAR2_ARTH
    arth_d(k)=arth_d(k-1)+increment
  end do
  temp=increment*NPAR2_ARTH
  k=NPAR2_ARTH
  do
    if (k >= n) exit
    k2=k+k
    arth_d(k+1:min(k2,n))=temp+arth_d(1:min(k,n-k))
    temp=temp+temp
    k=k2
  end do
end if
END FUNCTION arth_d


SUBROUTINE trapzd(func,par,a,b,s,n)
IMPLICIT NONE
REAL, INTENT(IN) :: a,b
REAL, DIMENSION(:), INTENT(IN):: par
REAL, INTENT(INOUT) :: s
INTEGER(4), INTENT(IN) :: n
INTERFACE
  FUNCTION func(x,par)
    REAL, DIMENSION(:), INTENT(IN) :: x
    REAL, DIMENSION(:), INTENT(IN) :: par
    REAL, DIMENSION(size(x)) :: func
  END FUNCTION func
END INTERFACE
!This routine computes the nth stage of refinement of an extended trapezoidal rule. func is
!input as the name of the function to be integrated between limits a and b, also input. When
!called with n=1, the routine returns as s the crudest estimate of .b
!a f(x)dx. Subsequent
!calls with n=2,3,... (in that sequential order) will improve the accuracy of s by adding 2n-2
!additional interior points. s should not be modified between sequential calls.
REAL :: del,fsum
INTEGER(4) :: it
if (n == 1) then
  s=0.5*(b-a)*sum(func( (/ a,b /),par(:) ))
else
  it=2**(n-2)
  del=(b-a)/real(it)                                 !This is the spacing of the points to be added.
  fsum=sum(func(arth_d(a+0.5*del,del,it),par(:)))
  s=0.5*(s+del*fsum)                               !This replaces s by its refined value.
end if
END SUBROUTINE trapzd


FUNCTION qtrap(func,par,a,b)
IMPLICIT NONE
REAL, INTENT(IN) :: a,b
REAL, DIMENSION(:), INTENT(IN):: par
REAL :: qtrap
INTERFACE
  FUNCTION func(x,par)
    REAL, DIMENSION(:), INTENT(IN) :: x
    REAL, DIMENSION(:), INTENT(IN) :: par
    REAL, DIMENSION(size(x)) :: func
  END FUNCTION func
END INTERFACE
INTEGER(4), PARAMETER :: JMAX=30
REAL, PARAMETER :: EPS=1.0e-1
!Returns the integral of the function func from a to b. The parameter EPS should be set to
!the desired fractional accuracy and JMAX so that 2 to the power JMAX-1 is the maximum
!allowed number of steps. Integration is performed by the trapezoidal rule.
REAL :: olds
INTEGER(4) :: j
olds= 0.0 !Initial value of olds is arbitrary.
do j=1,JMAX
  call trapzd(func,par(:),a,b,qtrap,j)
  if (j > 5) then !Avoid spurious early convergence.
    if (abs(qtrap-olds) < EPS*abs(olds) .or. (qtrap == 0.0 .and. olds == 0.0)) RETURN
  end if
  olds=qtrap
end do
write(*,*)"Too many steps in 'qtrap'"
stop 1
END FUNCTION qtrap

END MODULE NumInt
