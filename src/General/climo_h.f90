!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
module climo_h
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
save
! daily average temperatures (oC)
real, allocatable :: temp(:,:)
! daily precipitation (mm/day)
real, allocatable :: precip(:,:)
! daily potential evapotranspiration (mm/day)
real, allocatable :: pet(:,:)
! daily mean shortwave radiation (W/m^2)
real, allocatable :: rad(:,:)
! incoming radiation above atmosphere (W/m²)
REAL :: radex(366)
! relative humidity (%)
real, allocatable :: rhum(:,:)
! wind velocity (m/s)
real, allocatable :: wind(:,:)
! climate variables for hourly model version
real, allocatable :: preciph(:,:)

INTEGER :: store_day, store_timestep	!Till: day and timestep for which ETP-values were stored

end module climo_h

