! Till: computationally irrelevant: outcommented unused variables
! 2012-09-14
!
! Till: computationally irrelevant: minor changes to improve compiler compatibility
! 2011-04-29
!
!The module file contains two modules: reservoir_h and lake_h
!

module lake_h
save

! Logical for reading or calculation ofthe number of reservoir per size class per year 
LOGICAL :: doacudyear
! number of small reservoirs per volume class (5 volume classes)
! rounded to nearest integer
!Allocatable      real acud(nmun,5)
real, allocatable :: acud(:,:) 
! number of acudes per volume class (5 volume classes)
! not rounded, for calculation of growth per year
!Allocatable       real acudfloat(nmun,5)
real, allocatable :: acudfloat(:,:)
! Andreas 
! number of acudes per volume class (5 volume classes) and year
! specific values for each year, read from small_rservoirs_year.dat
! calculation of growth per year is not used if data are in the file
!Allocatable       real acudfloatyear(nmun,5,200)
real, allocatable :: acudfloatyear(:,:,:)
! initial maximum storage capacity in one acude of each volume classes (m3)
real :: maxlake0(5)
! initial maximum storage capacity in one acude of each volume classes (m3)
real, allocatable :: maxlakesub0(:,:)
! water stored initially (m3)
!Allocatable       real lakewater0(nmun,5)
real, allocatable :: lakewater0(:,:)
! trend in number of acudes for each volume class (+-nbr/year)
!Allocatable       real laketrend(nmun,5)
real, allocatable :: laketrend(:,:)
! total storage capacity of sub-basin in small acudes
!Allocatable       real totalacud(nmun)
real, allocatable :: totalacud(:)
! maximum capacity of water storage (m3)
!Allocatable       real maxlakewater(nmun,5)
real, allocatable ::  maxlakewater(:,:)
! fraction of storage capacity per volume class on total storage
!Allocatable       real acudfraction(nmun,5)
real, allocatable :: acudfraction(:,:) 
! fraction of generated river runoff being routed through small acudes
real :: intercepted
! water stored (m3)
!Allocatable      real lakewater(366*nt,nmun,5)
real, allocatable :: lakewater(:,:,:)
! water stored (m3)
!Allocatable      real laketot(366*nt,nmun)
real, allocatable :: laketot(:,:)
! relative water storage (-)
!Allocatable      real laketotfrac(366*nt,nmun)
real, allocatable :: laketotfrac(:,:)
! evaporation from lakes, all classes (m3)
!Allocatable      real lakeevap(366*nt,nmun,5)
real, allocatable :: lakeevap(:,:,:)
! precipitation on lakes, all classes (m3)
!Allocatable      real lakeprec(366*nt,nmun,5)
real, allocatable :: lakeprec(:,:,:)
! consumptive water use (effective extraction) from lakes, all classes (m3)
!Allocatable      real lakeex(366*nt,nmun,5)
real, allocatable :: lakeex(:,:,:)
! area of water storages (km2)
!Allocatable      real lakearea(nmun,5)
real, allocatable :: lakearea(:,:)
! Percentage  of upper limit of reservoir size class that represents the maximum storage capacity 
! of the hypothetical representative reservoir of each class in subbasin(-)
real, allocatable :: maxlake_factor(:,:)
! Percentage  of maximum volume per reservoir size class (-)
real :: lake_vol0_factor(5)
! Change on the number of reservoir per size class (-)
real :: lake_change(5)
! river inflow (m3)
!Allocatable       real lakeinflow(366*nt,subasin)
real, allocatable :: lakeinflow(:,:)
! river outflow (m3)
!Allocatable       real lakeoutflow(366*nt,subasin)
real, allocatable ::  lakeoutflow(:,:)
!Allocatable       real muniret(366*nt,subasin)
real, allocatable :: muniret(:,:)
!alpha factor for area-volume relationship (Molle,1989a)[-]
real :: alpha_Molle(5)
!K factor for area-volume relationship (Molle,1989a)[-]
real :: damk_Molle(5)
!Qout=c*Hv**d rating curve of the spillway, where Hv is the water height above the spillway of hypothetical representative reservoir from each class
real :: damc_hrr(5)
real :: damd_hrr(5)
! maximum storage capacity of the hypothetical representative reservoir from each class (m**3)
real, allocatable :: maxlake(:,:)
! maximum storage capacity of the hypothetical representative reservoir from each class (m**3) time-dependent variable
real, allocatable :: maxstorcap_hrr(:,:,:)
! inflow discharge into the hypothetical representative reservoir from each class (m**3)
real, allocatable :: lakeinflow_hrr(:,:,:)
! outflow discharge from the hypothetical representative reservoir from each class (m**3)
real, allocatable :: lakeoutflow_hrr(:,:,:)
! water retention in the hypothetical representative reservoir from each class (m**3)
real, allocatable :: lakeretention_hrr(:,:,:)
! storage volume in the hypothetical representative reservoir from each class (m**3)
real, allocatable :: lakewater_hrr(:,:,:)
!Spillway overflow on the day before in the hypothetical representative reservoir from each class (m3/timestep)
real, allocatable :: outflow_last_hrr(:,:)
!Maximum water depth of the hypothetical representative reservoir from each class used for the overflow calculation (m)
real, allocatable :: hmax_hrr(:,:)
!Water volume above the spillway elevation on the day before in the hypothetical representative reservoir from each class (m**3/timestep)
real, allocatable :: volume_last_hrr(:,:)
!fraction of contributing area of a reservoir class within the sub-basin  (-)
real, allocatable :: lakefrarea(:,:)
!fraction of sub-basin area not controlled by reservoirs (-)
real, allocatable :: subfrarea(:)
!fraction of sub-basin area not controlled by reservoirs receiving outflow discharge from reservoir classes of smaller storage capacity (-)
real, allocatable :: subfrout(:)
!fraction of contributing area of a reservoir class receiving outflow discharge from reservoir classes of smaller storage capacity (-)
real, allocatable :: lakefrout(:,:)
!fraction of sub-basin area receiving outflow discharge from reservoir classes of smaller storage capacity (-)
real, allocatable :: lakecumfrout(:,:)
!generated runoff within the subbasin (m**3)
real, allocatable :: lakerunoff(:,:)

!sediment balance of small reservoirs
! sediment inflow into the small reservoirs (ton)
real, allocatable :: lakesedin(:,:)
! sediment outflow without retention in the small reservoirs (ton)
real, allocatable :: lakesedout(:,:)
! percentage of grain size g that flows into the small reservoirs (-)
real, allocatable :: lakefracsedin(:,:,:)
! percentage of grain size g that flows out of the sub-basin(-)
real, allocatable :: lakefracsedout(:,:,:)
! sediment inflow into the hypothetical representative reservoir from each class (ton)
real, allocatable :: lakesedin_hrr(:,:,:)
! sediment outflow from the hypothetical representative reservoir from each class (ton)
real, allocatable :: lakesedout_hrr(:,:,:)
! percentage of grain size g that flows into hypothetical representative reservoir from each class (-)
real, allocatable :: lakefracsedin_hrr(:,:,:,:)
! percentage of grain size g that flows out of hypothetical representative reservoir from each class(-)
real, allocatable :: lakefracsedout_hrr(:,:,:,:)
! sediment deposition at the hypothetical representative reservoir from each class(-)
real, allocatable :: lakedep_hrr(:,:,:)
! cumulative sediment deposition at the hypothetical representative reservoir from each class(-)
real, allocatable :: lakecumdep_hrr(:,:)
! storage capacity reduction of the hypothetical representative reservoir from each class(-)
real, allocatable :: lake_vollost(:,:,:)
! cumulative sediment deposition at the hypothetical representative reservoir from each class(-)
real, allocatable :: cumseddep(:)
! cumulative sediment deposition at the reservoir class(-)
real, allocatable :: lake_cumseddep(:,:,:)

end module lake_h
