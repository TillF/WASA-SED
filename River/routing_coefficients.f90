SUBROUTINE routing_coefficients(i, STATUS, flow, r_area, p)
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
!!    b           |m             |bottom width of channel
!!    d           |m             |depth of flow
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

INTEGER, INTENT(IN OUT) :: i, STATUS
INTEGER :: jj, k
REAL :: fps, d, b, a, qq1, rh, tt1, tt2, aa, phi8, phi9, phi11, phi12, sed_con
REAL :: q_bankful100, a_100
REAL:: percent, q_it, error, d_it, t
REAL :: vol, s1, s2, r_area, p, flow, p_bed, p_fp, manning_composite


! -----------------------------------------------------------------------
IF (STATUS == 1) THEN

! Calculation of bottom width
s1 = r_sideratio(i)	
s2 = r_sideratio_fp(i)
d = r_depth(i)
bottom_width(i) = r_width(i) - 2. * d * r_sideratio(i)
b = r_width(i) - 2. * d * r_sideratio(i)

!! check if bottom width (b) is < 0
IF (b <= 0.) THEN
  b = 0.
  r_sideratio(i) = 0.
  b = .5 * r_width(i)
  r_sideratio(i) = (r_width(i) - b) / (2. * d)
END IF

!! compute flow and travel time at bankfull depth
p = b + 2. * d * SQRT(r_sideratio(i) * r_sideratio(i) + 1.)
a_100 = b * d + r_sideratio(i) * d * d
rh = a_100 / p
area_bankful(i) = a_100
q_bankful100 = a_100 * rh ** 0.6666 * Sqrt(r_slope(i))/manning(i)


!!Iteration to determine the correct flow area a for the given input discharge r_qin
t=0
! for flow in the main channel
if (r_qin(2,i).eq.0) then

elseif (r_qin(2,i).lt.q_bankful100) then
 error=1
 percent=0.01
!loop until change between given and iterated discharge becomes reasonably small
   do while (error.gt.0)
    d_it = (1 - percent) * r_depth(i)
    p = b + 2. * d_it * SQRT(r_sideratio(i) * r_sideratio(i) + 1.)
    a = b * d_it + r_sideratio(i) * d_it * d_it
    rh = a / p
    q_it = a * rh ** 0.6666 * Sqrt(r_slope(i))/manning(i)
    error = q_it - r_qin(2,i)
    percent = percent + 0.01
	t=t+1
    if (t.gt.99.or.q_it.lt.0) then
     write(*,*) 'Problems in routing_coefficients.f90: Bankful flow is used for initial river storage'
	 a=a_100
     exit
    endif
   enddo
!for flow in the main channel and the floodplains
elseif (r_qin(2,i).gt.q_bankful100) then
 error=-1
 percent=0.01
   do while (error.lt.0)
    d_it = percent * r_depth(i)    !=(1. + percent) * r_depth(i) - r_depth(i)	!iterative estimation of depth
    p = b + 2. * r_depth(i) * SQRT(1. + s1 * s1) + (r_width_fp(i)-r_width(i)) &	!Till: perimeter
    + 2. * d_it * SQRT(1 + s2 * s2)
! Calculation of composite Manning factor taking into account the roughness values of the floodplains (for derivation of composite Manning: siehe Wasserbauskriptum)
	p_bed = b + 2. * r_depth(i) * SQRT(1. + s1 * s1)
	p_fp = (r_width_fp(i)-r_width(i)) + 2. * d_it * SQRT(1 + s2 * s2)
	manning_composite = ((p_bed * manning(i)**1.5 + p_fp * manning_fp(i)**1.5)/p)**(0.666)	
  	a = a_100 + r_width_fp(i) * d_it + s2 * d_it * d_it							!Till: cross-section area
    rh = a / p																	!Till: hydraulic radius
    q_it = a * rh ** 0.6666 * SQRT(r_slope(i))/manning_composite				!Till: discharge according to Manning's equation
    error = q_it - r_qin(2,i)
    percent = percent + 0.01
    t=t+1
    if (t.gt.10000) then
     write(*,*) 'Problems in routing_coefficients.f90: floodplain calculations'
     stop
    endif
   enddo
endif

if (q_it.lt.0) a=0

!Calculation of initial water volume for each river stretch [m3]
if (r_qin(2,i).eq.0) then
  r_storage(i) = 0. 
else
  r_storage(i) =   a * r_length(i) * 1000
  sed_storage(i,1) = 0.
endif

!! Calculation of flow velocity [m/s]
if (a.eq.0) then
	velocity(i)=0.
else
  velocity(i)= q_it/a
endif

ENDIF
! -----------------------------------------------------------------------
IF (STATUS == 2) THEN
!CALCULATION OF FLOW PARAMETERS FOR CONTINUOUS, PERENNIAL FLOW

!! calculate volume of water in reach
vol = r_storage(i)

!! calculate cross-sectional area of flow, Equation 23.2.3
r_area = vol / (r_length(i) * 1000.)

!! calculate depth of flow, Equation 23.2.4 or 23.2.5
s1 = r_sideratio(i)	
s2 = r_sideratio_fp(i)
if (vol.lt.0.or.vol.eq.0) then
	r_depth_cur(i) = 0.
elseif (r_area <= area_bankful(i)) THEN
  r_depth_cur(i) = SQRT(r_area / s1 + bottom_width(i) * bottom_width(i) / (4. * s1 * s1)) - bottom_width(i) / (2. * s1)
  IF (r_depth_cur(i) < 0.) r_depth_cur(i) = 0.
ELSE
  r_depth_cur(i) = SQRT((r_area - area_bankful(i)) / s2 + (r_width_fp(i)*r_width_fp(i)) / (4*s2*s2)) - r_width_fp(i) / (2*s2)
  IF (r_depth_cur(i) < 0.) r_depth_cur(i) = 0.
  r_depth_cur(i) = r_depth_cur(i) + r_depth(i)
END IF

!! calculate wetted perimeter
p = 0.
if (vol.lt.0.or.vol.eq.0) then
  p = 0.
elseIF (r_depth_cur(i) <= r_depth(i)) THEN
  p = bottom_width(i) + 2. * r_depth_cur(i) * SQRT(1. + s1 * s1)
ELSE
  p = bottom_width(i) + 2. * r_depth(i) * SQRT(1. + s1 * s1) + (r_width_fp(i)-r_width(i)) &
      + 2. * (r_depth_cur(i) - r_depth(i)) * SQRT(1+s2*s2)
END IF

!! calculate hydraulic radius
rh = 0.
if (vol.lt.0.or.vol.eq.0) then
  rh = 0.
elseIF (p > 0.01) THEN
  rh = r_area / p
ELSE
  rh = 0.
END IF

!! Calculation of flow in reach with Manning Equation (m^3/s)
if (r_depth_cur(i).lt.r_depth(i).or.r_depth_cur(i).eq.r_depth(i)) then
  flow = r_area * rh ** 0.6666 * Sqrt(r_slope(i))/manning(i)
elseif (r_depth_cur(i).gt.r_depth(i)) then
! Calculation of composite Manning factor taking into account the roughness values of the floodplains 
  p_bed = b + 2. * r_depth(i) * SQRT(1. + s1 * s1)
  p_fp = (r_width_fp(i)-r_width(i)) + 2. * d_it * SQRT(1 + s2 * s2)
  manning_composite = ((p_bed * manning(i)**1.5 + p_fp * manning_fp(i)**1.5)/p)**(0.666)
  flow = r_area * rh ** 0.6666 * Sqrt(r_slope(i))/manning_composite
endif


ENDIF
! -----------------------------------------------------------------------
IF (STATUS == 3) THEN

!CALCULATION OF FLOW PARAMETERS FOR DISCONTINUOUS, EPHEMERAL FLOW

! Calculation of bottom width
s1 = r_sideratio(i)	
s2 = r_sideratio_fp(i)
d = r_depth(i)
b = r_width(i) - 2. * d * r_sideratio(i)

!! check if bottom width (b) is < 0
IF (b <= 0.) THEN
  b = 0.
  r_sideratio(i) = 0.
  b = .5 * r_width(i)
  r_sideratio(i) = (r_width(i) - b) / (2. * d)
END IF

!! compute flow and travel time at bankfull depth
p = b + 2. * d * SQRT(r_sideratio(i) * r_sideratio(i) + 1.)
a_100 = b * d + r_sideratio(i) * d * d
rh = a_100 / p
area_bankful(i) = a_100
q_bankful100 = a_100 * rh ** 0.6666 * Sqrt(r_slope(i))/manning(i)

!!Iteration to determine the correct flow area a for the given input discharge r_qin
t=0
! for flow in the main channel
if (r_qin(2,i).eq.0) then

elseif (r_qin(2,i).lt.q_bankful100) then
 error=1.
 percent=5.e-2
!loop until change between given and iterated discharge becomes reasonably small
   do while (error.gt.0)
    d_it = (1. - percent) * r_depth(i)
    p = b + 2. * d_it * SQRT(r_sideratio(i) * r_sideratio(i) + 1.)
    a = b * d_it + r_sideratio(i) * d_it * d_it
    rh = a / p
    q_it = a * rh ** 0.6666 * Sqrt(r_slope(i))/manning(i)
    error = q_it - r_qin(2,i)
    percent = percent + 5.e-2
	t=t+1
    if (t.eq.19.or.q_it.lt.0) then
	 if (r_qin(2,i).lt.q_it.or.error.lt.0) then
!       write(*,*) 'Problems in routing_coefficients.f90: Very small flow is used for beginning of ephermeal flow', i
	 else
!       write(*,*) 'Problems in routing_coefficients.f90: Bankful flow is used for beginning of ephermeal flow', i
	   a=a_100
     endif
     exit
    endif
   enddo
!for flow in the main channel and the floodplains
elseif (r_qin(2,i).gt.q_bankful100) then
 error=-1.
 percent=5.e-2
   do while (error.lt.0)
    d_it = percent * r_depth(i)    !=(1. + percent) * r_depth(i) - r_depth(i)
    p = b + 2. * r_depth(i) * SQRT(1. + s1 * s1) + (r_width_fp(i)-r_width(i)) &
    + 2. * d_it * SQRT(1 + s2 * s2)
	! Calculation of composite Manning factor taking into account the roughness values of the floodplains 
	p_bed = b + 2. * r_depth(i) * SQRT(1. + s1 * s1)
	p_fp = (r_width_fp(i)-r_width(i)) + 2. * d_it * SQRT(1 + s2 * s2)
	manning_composite = ((p_bed * manning(i)**1.5 + p_fp * manning_fp(i)**1.5)/p)**(0.666)
    a = a_100 + r_width_fp(i) * d_it + s2 * d_it * d_it
    rh = a / p
    q_it = a * rh ** 0.6666 * Sqrt(r_slope(i))/manning_composite
    error = q_it - r_qin(2,i)
    percent = percent + 5.e-2
    t=t+1

    if (t.gt.1000) then
     write(*,*) 'Problems in routing_coefficients.f90: Bankful flow is used for beginning of ephermeal flow', i
     stop
    endif
   enddo
endif


if (q_it.lt.0) a=0

!! Calculation of flow velocity [m/s]
if (a.eq.0) then
	velocity(i)=0.
else
  velocity(i)= q_it/a
endif

flow = r_qin(2,i)
r_area = a

ENDIF


RETURN
END SUBROUTINE routing_coefficients
