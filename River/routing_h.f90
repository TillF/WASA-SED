module routing_h
save

!ROUTING ARRAYS
! Arrays used for the routing of q
integer, allocatable :: upbasin(:), downbasin(:)
! qin: runoff into the Sub-basin [m**3/s]
!Allocatable      real qin (372,nmun)
real, allocatable :: qin(:,:)
! qout: runoff out of the Sub-basin [m**3/s]
!Allocatable      real qout (372,nmun)
real, allocatable :: qout(:,:)
!river inflow (m3/s) at time t and t+1 (2,subasin)
real, allocatable :: r_qin(:,:)
!river outflow (m3/s) at time t and t+1 (2,subasin)
real, allocatable :: r_qout(:,:)
!storage volume (m3) in river reach
real, allocatable :: r_storage (:)
!array stores runoff from hillslopes (m3/s) and erosion (tons), if it is generated outside WASA
real, allocatable :: runoff(:,:)
real, allocatable :: sediment(:,:,:)

! ORIGINAL WASA ROUTING SCHEME	
! routing parameter (lag time and retention time) (days)
!Allocatable      real prout(nmun,2)
real, allocatable :: prout (:,:)
! hrout(imun,7): unit hydr. ordinate [1/d] of routing function
!Allocatable     real hrout(7,nmun)
real, allocatable :: hrout(:,:)

! RIVER CHARACTERISTICS
! river channel width (bankful width) (m)
real, allocatable ::	r_width(:)
! river channel depth (m)
real, allocatable :: r_depth(:)
! current river depth at calculation step (m)
real, allocatable :: r_depth_cur(:)
! river slope(m/m)
real, allocatable :: r_slope(:)
!river side ratio: change in horizontal distance per unit change in vertical distance on channel side slopes (m/m)
real, allocatable :: r_sideratio(:)
! river bottom width of floodplain (m)
real, allocatable :: r_width_fp(:)
! river side ratio for the floodplain (m/m)
real, allocatable :: r_sideratio_fp(:)
! river length (km)
real, allocatable :: r_length(:)
! manning's n(-)
real, allocatable :: manning(:)
! manning's n for the floodplains (-)
real, allocatable :: manning_fp(:)
! effective hydraulic conduvtivity Ksat(mm/hr)
real, allocatable :: r_ksat(:)
! river erodibility factor (-)
real, allocatable :: r_efactor(:)
! river cover factor (-)
real, allocatable :: r_cover(:)
! percentage solid rock in river bed (-)
real, allocatable :: r_rock(:)
! baseflow alpha factor for bank storage (days)
real, allocatable :: r_alpha(:)
! weighting factor controlling relative importance of inflow rate and outflow rate in determining storage on reach
real, allocatable :: msk_x(:)
!storage time constant for the reach on current day (hours)
real, allocatable :: msk_k(:)
! Inflow of spring water at the start of a river
real, allocatable :: Q_spring(:)


!MUSKINGUM variables
!!  cross-sectional area of flow at bankfull (m^2)
real, allocatable :: area_bankful(:)
!!  bottom width of main channel depth (m)
real, allocatable :: bottom_width(:)
!!    phi5(i)    |m^3/s         |flow rate when reach is at bankfull depth
real, allocatable :: phi5(:)
!!    phi10(i)   |hr            |storage time constant for reach at bankfull depth (ratio of storage to discharge)
real, allocatable :: phi10(:)
!!    phi13(i)   |hr            |storage time constant for reach at 0.1 bankfull depth (low flow) (ratio of storage to discharge)
real, allocatable :: phi13(:)
! flow velocity (m/s)
real, allocatable :: velocity(:)


!SEDIMENT ROUTING VARIABLES
!! amount of sediment stored in river reach for each particle size class k (tons/timestep)
real, allocatable :: sed_storage(:,:)
!! sediment entering reach for each particle size class k (tons/timestep)
real, allocatable :: sediment_in(:,:)
!! sediment leaving reach for each particle size class k (tons/timestep)
real, allocatable :: sediment_out(:,:)
!! total sediment leaving reach (sum of all size classes) (tons/timestep)
real, allocatable :: qsediment_t(:)
!! daily total sediment output (tons/timestep) (372,subasin)
real, allocatable :: qsediment(:,:)
!! sediment concentration leaving reach (sum of all size classes) (tons/timestep)
real, allocatable :: r_sediment_concentration(:)
!! linear parameter for calculating sediment reentrained in channel sediment routing (-)
real, allocatable :: spcon(:)
!! exponent parameter for calculating sediment reentrained in channel sediment routing (-)
real, allocatable :: spexp(:)
!! peak flow rate adjustment factor for sediment routing in the channel. Allows impact of peak flow rate on sediment routing and channel reshaping to be taken into account
real prf
!! sediments deposited on river bottom in tons/(timestep and river stretch)
real, allocatable :: river_deposition(:,:)
!! sediment reentrained in water by channel degradation in tons/(timestep and river stretch)
real, allocatable :: river_degradation(:,:)
!! temporary storage of sediment on top of the riverbed in river reach [tons]
real, allocatable :: riverbed_storage(:,:) !!! why 2 dimensions?


!BEDLOAD MODELLING VARIABLES
! bedload(subasin,5) rate in (kg/s) as submerged weight for five bedload formula
real, allocatable :: bedload(:,:)
! D50 parameers in (m)
real, allocatable :: D50(:)


! VARIABLES FOR WATER TRANSPOSITIONS BETWEEN MUNIs OR EZGs
! Muni-ID of start Sub-basin of transposition
! and flag indicating source of water (1:from acude, 2:from river)
INTEGER :: trans_start(2,2)
! Muni-ID of destination Sub-basin of transposition
! and flag indicating source (1:to acude, 2:to river)
INTEGER :: trans_end(2,2)
! flow rate of transposition (m**3/s)
REAL :: q_trans(2)
! transport losses of flow rate (% of input flow rate)
REAL :: loss_trans(2)
! year of begin of operation of transposition
INTEGER :: y_trans(2)




end module routing_h