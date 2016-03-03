SUBROUTINE muskingum (i, flow, r_area,h)

! Till: computationally irrelevant: minor changes to improve compiler compatibility
! 2011-04-29

!Till: re-modified transmission losses, various optimizations 
!2009-12-17
 
!Till: modified transmission losses (affect rout instead of storage in reach); infiltration is restricted to bottom width, not full wetted perimeter
!2008-10-30

!! subroutine routes a daily or hourly water flow through a reach using the Muskingum method

!!    ~ ~ ~ INCOMING VARIABLES ~ ~ ~
!! r_sideratio(:)   |none        |change in horizontal distance per unit change in vertical distance on channel side slopes

!!    ~ ~ ~ OUTGOING VARIABLES ~ ~ ~
!! flow             |m^3/s       |average flow on day in reach
!! r_Vin(1,i)		|m3			 |volume of water into reach on current day at beginning of time step
!! r_Vin(2,i)		|m3			 |volume of water entering reach on day at end of the time step
!! r_Vout(1,i)		|m3			 |volume of flow out of reach on current day at beginning of time step
!! r_Vout(2,i)		|m3			 |water leaving reach on day at end of time step
!! r_ksat(i)        |mm/hr       |effective hydraulic conductivity of main channel alluvium
!! r_infil          |m^3 H2O     |transmission losses from reach in time step

!! LOCAL DEFINTIONS
!! vol              |m^3 H2O     |volume of water in reach at beginning of day
!! r_area           |m^2         |cross-sectional area of flow
!! c                |none        |change in horizontal distance per unit change in vertical distance on channel side slopes; always set to 2 (slope=1/2)
!! r_depth_cur(i)   |m           |depth of flow on day
!! rttime           |hr          |reach travel time
!! tbase            |none        |flow duration (fraction of 24 hr)
!! topw             |m           |top width of main channel
!! r_evp            |m^3 H2O     |evaporation from reach per timestep


use routing_h
use time_h
use hymo_h
use common_h
use climo_h

IMPLICIT NONE
INTEGER, INTENT(IN) :: i !internal subbasin-ID
INTEGER, INTENT(IN) :: h !subdaily timestep
integer :: dummy
real ::  r_area,  p !wetted perimeter
real ::  rttime, topw, flow, r_evp,r_infil, dummy2 !,vol, c, rh, tbase, s1, s2
real :: c0, c1, c2, c3, yy, total_water, total_losses !, Fr

!! Initialise water and sediment storage in each reach   
if (t == tstart .and. d == 1 .and. h == 1) then ! Till: is this necessary (already covered by clause below)? 
    call routing_coefficients(i,1, dummy2, dummy2, dummy2)
endif
!-------------------------------------------------------------

r_qout(2,i) = -1. !flag for "not yet computed"

if (r_storage(i) == 0.) then
    ! Calculation of discharge coefficients for ephemeral rivers
    if (r_qin(2,i) > 1.e-3) then
      call routing_coefficients (i,3,flow,r_area,p)
    ! Calculation of flow time [h]
      rttime = r_length(i)*1000./(velocity(i)*3600.) !Till: this uses the velocity of the previous timestep (?), ii
      if (rttime > dt) then
          r_qout(2,i) =      dt / rttime * r_qin(2,i)      !runoff gets stored instead of passing thru 
          r_storage(i)= (1- dt / rttime) * r_qin(2,i)*3600.*dt
          velocity(i) = velocity(i) * dt / rttime !gradually decrease velocity to avoid oscillations
      endif
    else
      r_qout(2,i) = 0. 
      r_storage(i) = r_qin(2,i)*3600.*dt 
      velocity(i) = velocity(i) * 0.5 !gradually decrease velocity to avoid oscillations
    endif
else
    ! Calculation of discharge coefficients for perennial rivers
  call routing_coefficients (i,2,flow,r_area,p)
end if
!------------------------------------------------------------

if (r_qout(2,i) == -1.) then !Till: do Muskingum routing unless already treated as ephemeral above
    !! Compute coefficients
    yy = dt / msk_k(i)
    c0 = yy  + 2. * (1. - msk_x(i))
    c1 = (yy + 2. * msk_x(i))  / c0
    c2 = (yy - 2. * msk_x(i))  / c0
    !c3 = (2. * (1. - msk_x(i)) - yy) / c0 
    c3 = 1 - c1 - c2   !equivalent to line above, but faster and numerically more stable

    !! Compute new outflow r_qout2
    IF (t == tstart .AND. d == 1 .and. h == 1) THEN
      r_qout(2,i) = r_qin(2,i) 
    ELSE
      r_qout(2,i) = c1 * r_qin(1,i) + c2 * r_qin(2,i) + c3 * r_qout(1,i)
      r_qout(2,i) = min (r_qout(2,i), r_storage(i)/(3600.*dt) + r_qin(2,i)) !Till: not more than the inflow and the storage can flow out of the reach	
    END IF
    IF (r_qout(2,i) < 0.) r_qout(2,i) = 0.

    !! Calculate flow velocity [m/s]
    if (r_area > 0.) then
      velocity(i)= flow/ r_area
    else
        velocity(i) = velocity(i) * 0.5 !gradually decrease velocity to avoid oscillations
    endif

end if !muskingum




!! Calculate travel time in [h]
!IF (flow > 1.e-4) THEN
!  rttime = r_length(i) * r_area / (3.6 * flow)
!END IF

!! Calculate transmission losses via riverbed infiltration
!r_infil = 0.
!! Calculate only if groundwater contributions are zero 
!if (gw_recharge(d,i) == 0..or.subflow(d,i) == 0.) then
!  r_infil=0.
!else
!  if (r_qout(2,i) > 1.e-2) then
    r_infil = r_ksat(i)* 1e-3 * dt * (r_length(i)* 1000.) * p ![m3] mm/h * h * km
    
!  END IF
!endif

!! Calculate evaporation of river stretch (in m3)
r_evp = 0.
IF (r_qout(2,i) > 1.e-3) THEN
  IF (r_depth_cur(i) <= r_depth(i)) THEN
    topw = bottom_width(i) + 2. * r_depth_cur(i) * r_sideratio(i)	! width of channel at water level [m]
  ELSE
    topw = r_width_fp(i) + 2. * (r_depth_cur(i) - r_depth(i)) * r_sideratio_fp(i)	! width of channel at water level [m]
  END IF
  r_evp = (pet(d,i)/24.)* 1e-3 * dt * (r_length(i)*1000.) * topw	!river evaporation [m³]
  if (r_evp < 0.) r_evp = 0.
END IF


total_water =  r_storage(i) + r_qout(2,i) *dt*3600. !total amount of water available in this timestep [m3]
total_losses = r_infil + r_evp

if (total_water == 0.) then !no water in reach
        r_infil = 0.
        r_evp   = 0.
        total_losses = 0.
end if

if (total_water < total_losses) then !if there is less water in the channel to fulfill seepage and evaporative demand...
        r_infil = r_infil * total_water / total_losses   !...reduce both proportionally
        r_evp   = r_evp   * total_water / total_losses  
        total_losses = total_water                        !reduce total losses
end if
if (f_river_infiltration) river_infiltration_t(d, h, i) = r_infil !store for output
if (f_actetranspiration)       aet_t(d, h, i) = aet_t(d, h, i) + r_evp / (area(i)*1e6) !store for output
if (f_daily_actetranspiration) aet  (d,    i) = aet(d,      i) + r_evp / (area(i)*1e6) !store for output


!update amount of water in channel at end of the time step
   r_storage(i) = r_storage(i) + (r_qin(2,i) - r_qout(2,i))*dt*3600.	
   r_storage(i) = max(0.,r_storage(i))	!shouldn't happen, but sometimes still does due to rounding errors

 
if (total_losses > 0) then
! distribute losses (infiltration, evap) between storage and outflow
    r_storage(i) = r_storage(i) - total_losses * r_storage(i)          / total_water ![m3]
    r_qout(2,i)  = r_qout(2,i)  - total_losses * r_qout(2,i)           / total_water ![m3/s]
     !prevent negative values (rather safe than sorry)
    r_storage(i) = max(0.,r_storage(i))	
    r_qout(2,i)  = max(0.,r_qout(2,i))	

end if !total losses>0    

!!Calculation of Froude Number
! Fr=velocity(i)/(sqrt(9.81 * r_depth_cur(i)))


RETURN
END SUBROUTINE muskingum
