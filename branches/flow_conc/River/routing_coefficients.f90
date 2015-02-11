SUBROUTINE routing_coefficients(i, STATUS, flow, r_area, p)
! Till: computationally irrelevant: minor changes to improve compiler compatibility
! 2011-04-29

!Till: computationally relevant for very initial phase of simulation: fixed bug in iterative estimation of depth from discharge
!2009-12-18

!Till: computationally relevant: major changes, bugfixes and optimzation
!2009-12-17

!Till: computationally irrelevant: increased iteration counter for determination of river cross section
!2009-12-11


!!    this subroutine computes travel time coefficients for the Muskingum Routing
!!    along the main channel

!!    ~ ~ ~ INCOMING VARIABLES ~ ~ ~
!!    name        |units         |definition
!!    ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!!    r_depth(:)     |m             |average depth of main channel
!!    r_length(:)    |km            |length of main channel
!!    manning(:)     |none          |Manning's "n" value for the main channel
!!    r_slope(:)     |m/m           |average slope of main channel
!!    r_width(:)     |m             |average width of main channel
!!    ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~

!!    ~ ~ ~ OUTGOING VARIABLES ~ ~ ~
!!    name        |units         |definition
!!    ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!!    r_sideratio(:)   |none          |change in horizontal distance per unit
!!                               |change in vertical distance on channel side slopes
!!    area_bankful(:)          |m^2           |cross-sectional area of flow at bankfull
!!    bottom_width(:)    |m             |bottom width of main channel depth
!!    phi5(i)    |m^3/s         |flow rate when reach is at bankfull depth
!!    phi10(i)   |hr            |storage time constant for reach at bankfull depth (ratio of storage to discharge)
!!    phi13(i)   |hr            |storage time constant for reach at 0.1 bankfull depth (low flow) (ratio of storage to discharge)

!!    ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~

!!    ~ ~ ~ LOCAL DEFINITIONS ~ ~ ~
!!    name        |units         |definition
!!    ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!!    phi8        |m/s           |average velocity when reach is at bankfull depth
!!    phi9        |m/s           |wave celerity when reach is at bankfull depth
!!    phi11       |m/s           |average velocity when reach is at 0.1 bankfull depth (low flow)
!!    phi12       |m/s           |wave celerity when reach is at 0.1 bankfull depth (low flow)
!!    aa          |none          |area/area=1 (used to calculate velocity with
!!                               |Manning's equation)
!!    a           |m^2           |cross-sectional area of channel
!!    bottom_width(i)           |m             |bottom width of channel
!!    dep         |m             |depth of flow
!!    fps         |none          |change in horizontal distance per unit
!!                               |change in vertical distance on floodplain side
!!                               |slopes; always set to 4 (slope=1/4)
!!    jj          |none          |counter
!!    p           |m             |wetting perimeter
!!    qq1         |m^3/s         |flow rate for a specified depth
!!    rh          |m             |hydraulic radius of channel
!!    tt1         |km s/m        |time coefficient for specified depth
!!    tt2         |km s/m        |time coefficient for bankfull depth
!!    ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~

!!    ~ ~ ~ SUBROUTINES/FUNCTIONS CALLED ~ ~ ~
!!    Intrinsic: Sqrt
!!    SWAT: Qman

!!    ~ ~ ~ ~ ~ ~ END SPECIFICATIONS ~ ~ ~ ~ ~ ~

use routing_h
use common_h
use hymo_h
use model_state_io !Jose Miguel: in order to be able to use init_river_conds

implicit none
INTEGER, INTENT(IN) :: i, STATUS
REAL, INTENT(OUT) :: r_area, p, flow

!REAL :: fps, qq1, tt1, tt2, aa, phi8, phi9, phi11, phi12, sed_con
REAL :: q_bankful100, dep
REAL::  d_it !, q_it, percent, error,
REAL :: vol, s1, s2 


s1 = r_sideratio(i)	
s2 = r_sideratio_fp(i)


flow =-1000.
r_area=-100.
p=-1000.
r_depth_cur(i)=-1000.


! -----------------------------------------------------------------------
IF (STATUS == 1) THEN

	! Calculation of bottom width
	dep = r_depth(i)
	bottom_width(i) = r_width(i) - 2. * dep * r_sideratio(i)

	!! check if bottom width (bottom_width(i)) is < 0
	IF (bottom_width(i) <= 0.) THEN
	  write(*,*)"WARNING: to low river width for depth / side slopes combination at ",i,"th stretch, side slopes adjusted."
	  bottom_width(i) = .5 * r_width(i)
	  r_sideratio(i) = (r_width(i) - bottom_width(i)) / (2. * dep)
	END IF

!! compute flow and travel time at bankfull depth
	CALL calc_q_a_p(r_depth(i), q_bankful100, area_bankful(i), p)

	!!determine the correct flow area a for the given input discharge r_qin
	! for flow in the main channel
	if (r_qin(2,i) > 0.)  then
		d_it =calc_d(r_qin(2,i))
		CALL calc_q_a_p(d_it, flow, r_area, p)
	else
		r_area=0.
	endif

	if (flow < 0.) r_area=0.

	!Calculation of initial water volume for each river stretch [m3]
 IF (.NOT. doloadstate) then
	 if (r_qin(2,i) == 0.) then
		r_storage(i) = 0. 
	 else
		r_storage(i) =   r_area * r_length(i) * 1000.
		sed_storage(i,1) = 0.
	 endif
 END IF
	 !IF (doloadstate) CALL init_river_conds(trim(pfadn)//'river_storage.stat')!Jose Miguel: load river storage of previous run of the model
     !Till: outcommented previous line as initialisation has been done before already
ENDIF

IF (STATUS == 3) THEN
	vol = r_qin(2,i) * dt * 3600. !Till: fill up storage with inflow - just to allow computations when storage was empty
ELSE
	vol = r_storage(i)
END IF



IF (STATUS == 1  .OR.STATUS == 2 .OR. STATUS == 3) THEN
! adjust water depth to current volume in reach
	!! calculate volume of water in reach

	!! calculate cross-sectional area of flow, Equation 23.2.3
	r_area = vol / (r_length(i) * 1000.)

	!! calculate depth of flow, Equation 23.2.4 or 23.2.5
	if (vol <= 0.) then
		r_depth_cur(i) = 0.
	elseif (r_area <= area_bankful(i)) THEN
	  r_depth_cur(i) = SQRT( r_area                    / s1 + bottom_width(i) * bottom_width(i) / (4. * s1 * s1)) - bottom_width(i) / (2. * s1) !Till: chose only positve solution of quadratic equation
	  IF (r_depth_cur(i) < 0.) r_depth_cur(i) = 0.	!Till: this should never occur
	ELSE
	  r_depth_cur(i) = SQRT((r_area - area_bankful(i)) / s2 + (r_width_fp(i)*r_width_fp(i))     / (4. * s2 * s2)) - r_width_fp(i)   / (2. * s2) !Till: chose only positve solution of quadratic equation
	  IF (r_depth_cur(i) < 0.) r_depth_cur(i) = 0.	!Till: this should never occur
	  r_depth_cur(i) = r_depth_cur(i) + r_depth(i)	
	END IF
END IF

! -----------------------------------------------------------------------
IF (STATUS >= 2) THEN
	!CALCULATION OF FLOW PARAMETERS FOR CONTINUOUS, PERENNIAL FLOW
	call calc_q_a_p(r_depth_cur(i), flow, r_area, p)
ENDIF


! -----------------------------------------------------------------------
IF (STATUS == 33) THEN	!obsolete

	!CALCULATION OF FLOW PARAMETERS FOR DISCONTINUOUS, EPHEMERAL FLOW

   !determine the flow area a for the given input discharge r_qin
	if (r_qin(2,i) /= 0.) then
		d_it =calc_d(r_qin(2,i))
		call calc_q_a_p(r_depth_cur(i), flow, r_area, p)
	else
  			r_area=0.
			p=0.
			flow=0.
			d_it=0.
	endif


	r_depth_cur(i)=d_it			!Till: set current depth if river			!!why was r_depth_cur not set for STATUS=3?


	flow = r_qin(2,i)			!Till: then why is q_it computed above?
	

ENDIF


!! Calculation of flow velocity [m/s]
if (r_area == 0.) then
	velocity(i)=0.
else
  velocity(i)= flow/r_area
endif



contains    
REAL FUNCTION calc_q(dep)
!compute discharge for a given depth dep in the current subreach indexed with i
!i:reach_index
!dep: reach_index
! Calculation of composite Manning factor taking into account the roughness values of the floodplains (for derivation of composite Manning: siehe Wasserbauskriptum)
implicit none

real, intent(in) :: dep	!water depth, for which the discharge is to be calculated
real :: d_fp = 0.	!depth of water in flood-plain
real :: q_ch, q_fp, a_ch, a_fp, p_fp, p_ch 

	if (dep<0.) then
		write(*,*)"ERROR: negative river depth in routing_coefficients.f90"
		stop
	end if
	p_ch = bottom_width(i) + 2. * min(dep,r_depth(i)) * SQRT(1. + s1 * s1)	!wetted perimeter in channel
	a_ch = (bottom_width(i) + r_sideratio(i) * dep) * dep						!cross-section area
	q_ch = a_ch * (a_ch/p_ch) ** 0.6666 * SQRT(r_slope(i))/manning(i)				!Till: discharge according to Manning's equation

	if (dep <= r_depth(i)) then			
		a_fp = 0.	!no flow in floodplain
		p_fp = 0.
		q_fp = 0.
	else
		d_fp = dep - r_depth(i)
		a_fp = ((r_width_fp(i)-r_width(i)) + s2 * d_fp) * d_fp
		p_fp =  (r_width_fp(i)-r_width(i)) + 2. * d_fp * SQRT(1. + s2 * s2)		!flow on floodplains
		q_fp = a_fp * (a_fp/p_fp) ** 0.6666 * SQRT(r_slope(i))/manning_fp(i)				!Till: discharge according to Manning's equation (flood plain)
	end if
		
	calc_q  = q_ch + q_fp
	
END FUNCTION calc_q

SUBROUTINE calc_q_a_p(dep, q, a, p)
!compute discharge q. cross-section-area a and wetted perimeter p for a given depth dep in the current subreach indexed with i
!i:reach_index
!dep: reach_index
! Calculation of composite Manning factor taking into account the roughness values of the floodplains (for derivation of composite Manning: siehe Wasserbauskriptum)
implicit none
real, intent(in) :: dep	!water depth, for which the discharge is to be calculated
real, intent(out) :: q, a, p
real :: d_fp = 0.	!depth of water in flood-plain
real :: q_ch, q_fp, a_ch, a_fp, p_fp, p_ch

	if (dep<0.) then
		write(*,*)"ERROR: negative river depth in routing_coefficients.f90"
		stop
	end if
	
	if (dep == 0.) then
		a = 0.
		p = 0.
		q = 0.
		return
	end if

	p_ch = bottom_width(i) + 2. * min(dep,r_depth(i)) * SQRT(1. + s1 * s1)	!wetted perimeter in channel
	a_ch = (bottom_width(i) + r_sideratio(i) * dep) * dep						!cross-section area
	q_ch = a_ch * (a_ch/p_ch) ** 0.6666 * SQRT(r_slope(i))/manning(i)				!Till: discharge according to Manning's equation

	if (dep <= r_depth(i)) then			
		a_fp = 0.	!no flow in floodplain
		p_fp = 0.
		q_fp = 0.
	else
		d_fp = dep - r_depth(i)
		a_fp = ((r_width_fp(i)-r_width(i)) + s2 * d_fp) * d_fp
		p_fp =  (r_width_fp(i)-r_width(i)) + 2. * d_fp * SQRT(1. + s2 * s2)		!flow on floodplains
		q_fp = a_fp * (a_fp/p_fp) ** 0.6666 * SQRT(r_slope(i))/manning_fp(i)				!Till: discharge according to Manning's equation (flood plain)
	end if
		
	a = a_ch + a_fp		!simple summation of components
	p = p_ch + p_fp
	q = q_ch + q_fp
END SUBROUTINE calc_q_a_p

REAL FUNCTION calc_d(q)
!compute depth dep for a given discharge q in the current subreach indexed with i
!see http://nptel.iitm.ac.in/courses/IIT-MADRAS/Hydraulics/pdfs/Unit12/12_1b.pdf
!i:reach_index
!q: discharge
implicit none
real, intent(in) :: q	!discharge, for which the water depth is to be calculated
real :: dd, df	!differentials
real :: d_est0, d_est1, q_est0, q_est1, f0, f1, error_tolerance=0.01	!Till: max. relative error indicating convergence
integer :: j, max_iter=50			!Till: max number of iterations

	q_bankful100 = calc_q(r_depth(i))		!compute bankful discharge	!ii: compute only once and store in array
	q_est0       = q_bankful100
	d_est0       = r_depth(i)
	f0 = q_est0 - q

	if (q < q_bankful100) then
		d_est1=r_depth(i)/2.	!Till: initial estimate	for small discharge
	else
		d_est1=r_depth(i)*1.5	!Till: initial estimate	for high discharge
	end if
	q_est1 = calc_q(d_est1)
	!f1 = q_est1 - q

	do j=1,max_iter
		f1 = q_est1-q
		df= f1     - f0
		dd= d_est1 - d_est0  

		d_est0 = d_est1		!prepare next iteration step
		q_est0 = q_est1
		f0	   = f1

		d_est1 = max(0.,d_est1 - f1 / (df/dd))
		q_est1 = calc_q(d_est1)
		if (abs(q-q_est1)/q <= error_tolerance) exit	!Till: converged
	end do

    if (j>max_iter) write(*,*)"WARNING: iteration limit reached in  routing_coefficients.f90"
	calc_d = d_est1
END FUNCTION calc_d


END SUBROUTINE routing_coefficients


   
