SUBROUTINE route_sediments(i)

! Till: computationally irrelevant: removed unused parameters flow, rarea
! 2011-05-05

! Till: computationally irrelevant: minor changes to improve compiler compatibility
! 2011-04-29

!!    ~ ~ ~ PURPOSE ~ ~ ~
!!    this subroutine routes sediment through the river system;
!!    deposition is based on method velocity and degradation on stream

!!    ~ ~ ~ INCOMING VARIABLES ~ ~ ~
!!    name        |units         |definition
!!    ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!!    r_cover(:)   |none          |channel cover factor (0.0-1.0)
!!                               |0 channel is completely protected from
!!                               |  erosion by cover
!!                               |1 no vegetative cover on channel
!!    r_depth(:)     |m             |average depth of main channel
!!    r_efactor(:)  |none          |channel erodibility factor (0.0-1.0)
!!                               |0 non-erosive channel
!!                               |1 no resistance to erosion
!!    r_length(:)    |km            |initial length of main channel
!!    r_width(:)   |m             |average width of main channel
!!    ideg        |none          |channel degredation code
!!                               |0: do not compute channel degradation
!!                               |1: compute channel degredation (downcutting
!!                               |   and widening)
!!    phi(5,:)    |m^3/s         |flow rate when reach is at bankfull depth
!!    prf         |none          |Peak rate adjustment factor for sediment
!!                               |routing in the channel. Allows impact of
!!                               |peak flow rate on sediment routing and
!!                               |channel reshaping to be taken into account
!!    r_depth_cur(i)      |m             |depth of flow on day
!!    flow        |m^3/s         |average flow on day in reach
!!    sed_storage(:)    |metric tons   |amount of sediment stored in reach
!!    spcon       |none          |linear parameter for calculating sediment 0.0001-0.01
!!                               |reentrained in channel sediment routing
!!    spexp       |none          |exponent parameter for calculating sediment 1.0-1.5 
!!                               |reentrained in channel sediment routing
!!	  r_efactor(i) |none          |0.0-1.0    |channel cover factor: 0.0=channel is completely protected from erosion by cover; 1.0=no vegetative cover on channel

!!    varoute(3,:)|metric tons   |sediment
!!    ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~

!!    ~ ~ ~ OUTGOING VARIABLES ~ ~ ~
!!    name        |units         |definition
!!    ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!!    peakr       |m^3/s         |peak runoff rate in channel
!!    sed_storage(:)    |metric tons   |amount of sediment stored in reach
!!    sediment_out(:)      |metric tons   |sediment transported out of channel
!!                               |during time step
!!    ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~

!!    ~ ~ ~ LOCAL DEFINITIONS ~ ~ ~
!!    name        |units         |definition
!!    ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!!    river_degradation |metric tons   |sediment reentrained in water by channel degradation
!!    river_deposition  |metric tons   |sediment deposited on river bottom
!!    depdeg      |m             |depth of degradation/deposition from original
!!    depnet      |metric tons   |
!!    dot         |
!!    i           |none          |reach number
!!    volume        |m^3 H2O       |water in reach during time step
!!    vel_peak    |m/s           |peak flow velocity in reach
!!	  r_storage_previous  | m3	| river storage on previous time step
!!    ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~

!!    ~ ~ ~ ~ ~ ~ END SPECIFICATIONS ~ ~ ~ ~ ~ ~

use routing_h
use common_h
use time_h

IMPLICIT NONE
INTEGER, INTENT(IN OUT) :: i
INTEGER :: det, k
REAL ::  vel_peak,  conc_in, conc_max, depnet !peakr,peakflow,
REAL :: volume, sed_mass !flow, r_area
REAL ::   r_storage_previous !depdeg,dot,


!if no or only very little water flow in reach, no sediment calculation is performed
if (r_storage(i).eq.0.) then
 sediment_out(i,:) = 0.
 sed_storage(i,:) = 0.
 return
endif

IF (r_storage(i) .eq. 0. .and. r_qin(2,i) .lt. 1.e-1) then
 sediment_out(i,:) = 0.
 sed_storage(i,:) = 0.
 return
ENDIF

!! initialize water in reach during time step [m3]
r_storage_previous=r_storage(i)+(r_qout(2,i)-r_qin(2,i))*dt*3600.
volume = r_storage_previous+r_qin(2,i)*3600.*dt

!! check if it is an ephemeral river (as set in Muskingum.f90) that starts to flow
if (r_qout(1,i) .eq. 0. .and. r_storage(i) .eq. (r_qin(2,i)+3600.*dt)) then
  sediment_out(i,:) = 0.
  sed_storage(i,:) = 0.
  volume = r_storage(i)
  return
endif
!! do not perform sediment routing if no water in reach
IF (volume <= 0.01) then 
 sediment_out(i,:) = 0.
 sed_storage(i,:) = 0.
 RETURN
ENDIF
!! initialize sediment mass in reach during time step [tons]

prf= 1.
spcon(:)=  0.016111		!0.0001-0.01  
spexp (:)= 1.707			!1 - 1.5
det = 24/dt	!Amount of simulation step per day



!! Calculation of flow velocity [m/s]
  IF (velocity(i) < .010) THEN
    vel_peak = 0.01
    write(*,*) 'very small flow velocity for sediment transport in sub-basin: ', i
  ELSE
    vel_peak = prf * velocity(i)
  END IF
  IF (vel_peak > 5.) then
    vel_peak = 5.
    write(*,*) 'very large flow velocity for sediment transport in sub-basin ', i
  endif
  

!Loop through all sediment classes
do k=1, n_sed_class

  sed_mass = sediment_in(i,k)+ sed_storage(i,k)


! Calculation of current and maximum sediment carrying capacity concentration [ton/m3, kg/l]
  conc_in = sed_mass / volume
  conc_max = spcon(k) * vel_peak ** spexp(k)
  ! conc_max = conc_in

! Calculation of net amount of sediment deposited [tons]
  depnet = volume * (conc_max - conc_in)
  IF (depnet > 0.) THEN
! Calculation of degradation [tons]
! Extension to incorporate detachment limitations (e.g. rocky riverbeds, where nothing can be degradeded, r_rock: percentage cover of riverbed with rock)
! only allow degradation if discharge larger than 0.1 m3/s
	if (r_qout(1,i).gt.0.1) then   
	   river_degradation(i,k) = depnet * r_efactor(i) * r_cover(i) * (1. - r_rock(i))
	else
       river_degradation (i,k) = 0.
	endif

    river_deposition(i,k) = 0.
!  Check if sediment is available in temporary sediment storage of riverbed
!	if (river_degradation(i,k).gt.riverbed_storage(i,k)) river_degradation(i,k)=riverbed_storage(i,k) 
! Calculation of deposition [tons]
  ELSE
    river_deposition(i,k) = -depnet
    river_degradation(i,k) = 0.
  END IF

! Calculation of temporary storage of sediments in the riverbed of the stretch [tons]
!  riverbed_storage(i,k)=riverbed_storage(i,k) + river_deposition(i,k) - river_degradation(i,k)
!  if (riverbed_storage(i,k).lt.0.) riverbed_storage(i,k) = 0.

! Calculation of sediment balance [tons] (assuming instaneous mixing)
  sed_storage(i,k) = sed_storage(i,k)+ sediment_in(i,k) + river_degradation(i,k) - river_deposition(i,k)
!Calculation of sediment leaving the reach [tons]
  sediment_out(i,k) = sed_storage(i,k)/ volume * (r_qout(2,i)*3600.*dt)
  IF (sediment_out(i,k) < 0.) sediment_out(i,k) = 0.

  sed_storage(i,k) = sed_storage(i,k) - sediment_out(i,k)
  IF (sed_storage(i,k) < 0.) sed_storage(i,k) = 0.


! end of loop for all particle size classes
ENDDO

RETURN
END SUBROUTINE route_sediments
