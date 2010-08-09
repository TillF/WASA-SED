!
!The module file contains two modules: reservoir_h and lake_h
!

module reservoir_h
save

! Options to run the WASA Model without river and hillslope modules
INTEGER :: reservoir_balance	!option to read (0) or calculate (1) outflow discharges and reservoir levels 
INTEGER :: reservoir_check		!option to run the WASA Model with river and hillslope modules (0) or without them (1)
INTEGER :: reservoir_print		!option to print the output files within the timestep (0) or once a year (1)

!Reservoir water balance arrays
!simulation timestep
integer :: step
!hour (on day)
integer :: hour

!flag to calculate the ratio between the reservoir volumes given in the file cav.dat and that derived from the cross sections  [0 = initial value; 1 = after first calculation]
integer, allocatable :: fcav(:)
!flag to simulate reservoir routing during spillway overflow of the sub-basin's reservoir [0 = without time delay; 1 = with time delay]
integer, allocatable :: fvol_over(:)
!flag to identify the ocurrence of lateral inflow into the dowstream reservoir (0 = no and 1 = yes)
integer, allocatable :: latflow_res(:)
!arrays used for the routing of lateral inflow into the dowstream reservoir
integer, allocatable :: reservoir_down(:)
!year of construction of the dam in the sub-basin
integer, allocatable :: damyear(:)
!number of points from the stage-area-volume curve
integer, allocatable :: nbrbat(:)
!day of change in exploitation regime in the sub-basin's reservoir [-]
integer, allocatable :: dayexplot(:,:)
!time of start of the sediment management measures [m**3]
integer, allocatable :: operat_start(:)
!time of stop of the sediment management measures [m**3]
integer, allocatable :: operat_stop(:)
!maximum reservoir level for the sediment management measures [m**3]
integer, allocatable :: operat_elev(:)
!maximum water depth used for the overflow calculation (m)
real, allocatable :: hmax(:)
!water volume above the spillway elevation on the day before (m3/timestep)
real, allocatable :: volume_last(:)
!spillway overflow on the day before (m3/timestep)
real, allocatable :: outflow_last(:)
!effective water extraction (consumptive use) from the Sub-basin's reservoir [m**3]
real, allocatable :: damex(:,:)
!initial storage capacity in the sub-basin's reservoir [10**3 m**3]
real, allocatable :: storcap(:)
!initial dead volume of the sub-basin's reservoir [10**3 m**3]
real, allocatable :: damdead(:)
!elevation corresponding to dead volume of the sub-basin's reservoir [m]
real, allocatable :: elevdead(:)
!initial alert volume of the sub-basin's reservoir [10**3 m**3]
real, allocatable :: damalert(:)
!elevation corresponding to alert volume of the sub-basin's reservoir [m]
real, allocatable :: elevalert(:)
!bottom outlet elevation of the sub-basin's reservoir [m]
real, allocatable :: elevbottom(:)
!base-flow out of the sub-basin's reservoir (90% garanteed) [m**3/s]
real, allocatable :: damflow(:)
!percentage of Q90 released from the sub-basin's reservoir in regular years [-]
real, allocatable :: damq_frac(:)
!percentage of Q90 released from the sub-basin's reservoir in different seasons in the sub-basin's reservoir [-]
real, allocatable :: damq_frac_season(:,:)
!actual stored volume in the subbasin's reservoir [10**6 m**3]
real, allocatable :: volact(:,:)
!evaporation from the subbasin's reservoir [mm]
real, allocatable :: evapdam(:,:)
!evaporation from the subbasin's reservoir [m**3]
real, allocatable :: etdam(:,:)
!infiltration losses from the subbasin's reservoir [m**3]
real, allocatable :: infdam(:,:)
!initial maximum reservoir area in subbasin [ha]
real, allocatable :: maxdamarea(:)
!actual reservoir area [m**2]
real, allocatable :: damareaact(:)
!volume=k*Hv**alpha (Volume/heigth) relationship parameters
real, allocatable :: alpha_over(:)
real, allocatable :: k_over(:)
!area=a*Vol**b (Volume/area) relationship parameters
real, allocatable :: dama(:)
real, allocatable :: damb(:)
!Qout=c*Hv**d rating curve of the spillway, where Hv is the water height above the spillway
real, allocatable :: damc(:)
real, allocatable :: damd(:)
!volume=ff*h**3 (Volume/area) relationship parameters
real, allocatable :: forma_factor(:)
!maximum outflow released through the bottom outlets of the sub-basin's reservoir (90% garanteed) [m**3/s]
real, allocatable :: qoutlet(:)
!percentage of storage capacity that indicates the minimum storage volume for sediment release through the bottom outlets of the sub-basin's reservoir [-]
real, allocatable :: fvol_bottom(:)
!withdrawal water volume to supply the water use sectors (outflow discharge through the dam is not considered) [m**3/s] 
real, allocatable :: withdrawal(:)
!maximum retention in the subbasin's reservoir [m**3]
real, allocatable :: lakeret(:,:)
!water elevation from the initial stage-area-volume curve [m]
real, allocatable :: elev_bat0(:,:)
!reservoir area for a given elevation at the initial stage-area-volume curve [10**3 m**2]
real, allocatable :: area_bat0(:,:)
!reservoir volume for a given elevation at the initial stage-area-volume curve [10**3 m**3]
real, allocatable :: vol_bat0(:,:)
!water elevation from the actual stage-area-volume curve [m]
real, allocatable :: elev_bat(:,:)
!reservoir area for a given elevation at the actual stage-area-volume curve [10**3 m**2]
real, allocatable :: area_bat(:,:)
!reservoir volume for a given elevation at the actual stage-area-volume curve [10**3 m**3]
real, allocatable :: vol_bat(:,:)
!initial volume of the sub-basin's reservoir [10**3 m**3]
real, allocatable :: vol0(:)
!lateral inflow discharge into the subbasin's reservoir [m**3/s]
real, allocatable :: qlateral(:,:)
!inflow discharge into the subbasin's reservoir [m**3/s]
real, allocatable :: qinflow(:,:)
!overflow discharge out the sub-basin's reservoir [m**3/s]
real, allocatable :: overflow(:,:)
!controlled outflow by intake device in the sub-basin's reservoir [m**3/s]
real, allocatable :: qintake(:,:)
!controlled outflow relased through the bottom outlets of the sub-basin's reservoir [m**3/s]
real, allocatable :: qbottom(:,:)
!storage capacity in the subbasin's reservoir [10**6 m**3]
real, allocatable :: daystorcap(:,:)
!maximum reservoir area in the subbasin's reservoir [10**4 m**2]
real, allocatable :: daymaxdamarea(:,:)
!dead volume of the sub-basin's reservoir [10**3 m**3]
real, allocatable :: daydamdead(:,:)
!alert volume of the sub-basin's reservoir [10**3 m**3]
real, allocatable :: daydamalert(:,:)
!minimum level of the sub-basin's reservoir [mm]
real, allocatable :: dayminlevel(:,:)
!actual reservoir level in the subbasin [m]
real, allocatable :: damelevact(:)
!maximum water level in the sub-basin's reservoir [m]
real, allocatable :: maxlevel(:)
!initial minimum level in the sub-basin's reservoir [m]
real, allocatable :: minlevel(:)
!reservoir volume at the beginning of the simulation timestep in the sub-basin [10**6 m**3]
real, allocatable :: damvol0(:)
!reservoir level at the beginning of the simulation timestep in the sub-basin [m]
real, allocatable :: damelev0(:,:)
!reservoir level at the end of the simulation timestep in the sub-basin  [m]
real, allocatable :: damelev1(:,:)
!water volume of the reservoir sub-reach by summing up all cross sections' volume [m**3]
real, allocatable :: resreach_vol(:)
!precipitation over the reservoir area (mm)
real, allocatable :: res_precip(:,:)
!potential evapotranspiration (mm)
real, allocatable :: res_pet(:,:)
!total outflow discharge in the sub-basin's reservoir (m**3/s)
real, allocatable :: res_qout(:,:)

!Reservoir sedimentation (semres) arrays
!general parameters
!external IDs of cross sections in the sub-basin's reservoir
integer, allocatable :: id_sec_extern(:,:)
!number of cross sections in the sub-basin's reservoir
integer, allocatable :: nbrsec(:)
!number of points at the cross section of the reservoir in the sub-basin
integer, allocatable :: npoints(:,:)
!number of cases for the calculation of bed geometry changes between two adjacent points 
!at the cross section of the reservoir in the sub-basin
integer, allocatable :: geom(:,:)
!decrease of stored volume in the subbasin's reservoir [m**3]
real, allocatable :: decvolact(:,:)
!decrease of storage capacity in the subbasin's reservoir [m**3]
real, allocatable :: decstorcap(:,:)
!decrease of maximum reservoir area in subbasins [m**2]
real, allocatable :: decmaxdamarea(:,:)
!decrease of dead volume in the subbasin's reservoir [m**3]
real, allocatable :: decdamdead(:,:)
!decrease of alert volume in the subbasin's reservoir [m**3]
real, allocatable :: decdamalert(:,:)
!Manning's roughness for each cross section [m**(-1/3).s]
real, allocatable :: manning_sec(:,:)
!distance to the downstream cross section [m]
real, allocatable :: dist_sec(:,:)
!values at the x-axis for each point in the cross section of the reservoir
!(from left to right, seen from upstream) [m]
real, allocatable :: x_sec0(:,:,:)
!values at the y-axis for each point in the cross section of the reservoir
!(from left to right, seen from upstream) [m]
real, allocatable :: y_sec0(:,:,:)
!suspended sediment in the subbasin's reservoir [ton/timestep]
real, allocatable :: sed_susp(:,:)
!sediment retention in the subbasin's reservoir [ton/timestep]
real, allocatable :: sed_ret(:,:)
!volume of overflowed sediment out the subbasin's reservoir [ton/timestep]
real, allocatable :: sed_overflow(:,:)
!volume of released sediment through the water intake device of the subbasin's reservoir through the water intake [ton/timestep] 
real, allocatable :: sed_intake(:,:)
!volume of released sediment through the bottom outlets of the subbasin's reservoir through the water intake [ton/timestep] 
real, allocatable :: sed_bottom(:,:)
!lateral sediment inflow into the dowstream reservoirs for each particle size class g [ton/timestep]
real, allocatable :: sed_qlateral(:,:)
!sediment inflow discharge into the subbasin's reservoir [ton/timestep]
real, allocatable :: sed_inflow(:,:)
!sediment outflow discharge from the subbasin's reservoir [ton/timestep]
real, allocatable :: sed_outflow(:,:)
!mean diameter of particles considered in the sediment routing in the subbasin's reservoir [m]
real, allocatable :: diam(:)
!sedimentation in the sub-basin's reservoir [ton]
real, allocatable :: sedimentation(:,:)
!cumulative sedimentation in the sub-basin's reservoir [ton]
real, allocatable :: cum_sedimentation(:)
!minimum sediment concentration in the sub-basin's reservoir [mg/l]
real, allocatable :: min_conc(:)
!initial sediment concentration in the sub-basin's reservoir [mg/l]
real, allocatable :: sed_conc0(:)
!wet bulk density for sediments in the sub-basin's reservoir [ton/m**3]
real, allocatable :: wet_dens(:)
!! sediment release from reservoir for each particle size class k (tons/timestep)
real, allocatable :: res_sediment_out(:,:)
!grain size distribution of the sediment inflow into the subbasin's reservoir [-]
real, allocatable :: frsediment_in(:,:)
!grain size distribution of the sediment outflow from the subbasin's reservoir [-]
real, allocatable :: frsediment_out(:,:)
!dry bulk density of the sediment deposited in the subbasin's reservoir [ton/m**3]
real, allocatable :: dry_dens(:)
!calibration parameter for the determination of the active layer thickness [-]
real, allocatable :: factor_actlay(:)
!flag to control changes on sideslope at cross sections in the sub-basin's reservoir  [0 = changes on sideslope is not controlled; 1 = changes on sideslope is controlled avoiding steeper slopes by erosion processes]
integer, allocatable :: sed_flag(:)
!flag to simulate sediment routing by using of flushing technique in the sub-basin's reservoir  [0 = without flushing scenario; 1 = with flushing scenario]
integer, allocatable :: sed_routing_flag(:)
!amount of correlations between water inflow discharge and incoming sediment into the subbasin's reservoir [-]
integer, allocatable :: numbeqs(:)
!water discharge that represents the upper limit of applicability of each correlation 
!between water inflow discharge and incoming sediment into the subbasin's reservoir [-]
real, allocatable :: Q_max(:,:)
!Qsed=a*Q**b (water input x sediment input) relationship parameters [-]
real, allocatable :: param_a(:,:)
real, allocatable :: param_b(:,:)
!number of size distribution of the incoming sediment into the subbasin's reservoir [-]
integer, allocatable :: nbsizedist(:)
!water discharge value with a given size distribution of the incoming sediment into the subbasin's reservoir [-]
real, allocatable :: Q_refer(:,:)
!computed size distribution of the incoming sediment into the subbasin's reservoir [-]
real, allocatable :: perc(:,:,:)
!sediment inflow into the subbasin's reservoir related to grain size g [ton(timestep]
real, allocatable :: sedinflow_g(:,:,:)
!sediment inflow into the subbasin's reservoir related to grain size g [ton(timestep]
real, allocatable :: sedoutflow_g(:,:,:)

!hydraulic calculations
!mean reservoir level at the simulation timestep in the sub-basin used for the hydraulic calculation [m]
real, allocatable :: damelev_mean(:,:)
!location of the minimum elevation at the x-axis of each cross section in the sub-basin's reservoir [m]
real, allocatable :: x_minelev(:,:)
!minimum elevation of each cross section in the sub-basin's reservoir [m]
real, allocatable :: minelev_sec(:,:)
!bed slope of each cross section in the sub-basin's reservoir [m/m]
real, allocatable :: bedslope_sec(:,:)
!wetted area of each cross section in the sub-basin's reservoir [m**2]
real, allocatable :: area_sec(:,:)
!wetted area of each cross section at the reservoir sub-reach in the sub-basin [m**2]
real, allocatable :: resarea_sec(:,:)
!top width of each cross section in the sub-basin's reservoir [m]
real, allocatable :: topwidth_sec(:,:)
!storage volume of the sub-basin's reservoir by summing up all cross sections' volume [m**3]
real, allocatable :: resvol(:)
!reservoir's volume represented by two adjacent cross section [m**3]
real, allocatable :: resvol_sec(:,:)
!weighting parameter for each cross section in the sub-basin's reservoir [-]
real, allocatable :: weight_sec(:,:)
!discharge of each cross section in the sub-basin's reservoir [m**3/s]
real, allocatable :: discharge_sec(:,:)
!water depth of each cross section in the sub-basin's reservoir [m]
real, allocatable :: depth_sec(:,:)
!water elevation of each cross section in the sub-basin's reservoir [m]
real, allocatable :: watelev_sec(:,:)
!wetted perimeter of each cross section in the sub-basin's reservoir [m]
real, allocatable :: wetper_sec(:,:)
!hydraulic radius of each cross section in the sub-basin's reservoir [m]
real, allocatable :: hydrad_sec(:,:)
!mean velocity of each cross section in the sub-basin's reservoir [m/s]
real, allocatable :: meanvel_sec(:,:)
!slope of energy-grade line of each cross section in the sub-basin's reservoir [-]
real, allocatable :: energslope_sec(:,:)
!dynamic head of each cross section in the sub-basin's reservoir [m]
real, allocatable :: dynhead_sec(:,:)
!total head of each cross section in the sub-basin's reservoir [m]
real, allocatable :: tothead_sec(:,:)
!head loss between two adjacent cross sections in the sub-basin's reservoir due friction[m]
real, allocatable :: headloss_sec(:,:)
!local head loss between two adjacent cross sections in the sub-basin's reservoir [m]
real, allocatable :: locloss_sec(:,:)
!calculated total head of each cross section in the sub-basin's reservoir [m]
real, allocatable :: calctothead_sec(:,:)
!maximum area of each cross section in the sub-basin's reservoir [m**2]
real, allocatable :: maxarea_sec(:,:)
!maximum elevation of each cross section in the sub-basin's reservoir [m
real, allocatable :: maxelev_sec(:,:)
!maximum depth of each cross section in the sub-basin's reservoir [m]
real, allocatable :: maxdepth_sec(:,:)
!critical depth of each cross section in the sub-basin's reservoir [m]
real, allocatable :: crdepth_sec(:,:)
!critical elevation of each cross section in the sub-basin's reservoir [m]
real, allocatable :: crwatelev_sec(:,:)
!critical area of each cross section in the sub-basin's reservoir [m**2]
real, allocatable :: crarea_sec(:,:)
!critical top width of each cross section in the sub-basin's reservoir [m]
real, allocatable :: crtopwidth_sec(:,:)
!critical wetted perimeter of each cross section in the sub-basin's reservoir [m]
real, allocatable :: crwetper_sec(:,:)
!critical bed slope of each cross section in the sub-basin's reservoir [-]
real, allocatable :: crslope_sec(:,:)
!critical bed slope of each cross section in the sub-basin's reservoir [-]
real, allocatable :: crvel_sec(:,:)
!critical hydraulic radius of each cross section in the sub-basin's reservoir [m]
real, allocatable :: crhydrad_sec(:,:)
!normal water elevation of each cross section in the sub-basin's reservoir [m]
real, allocatable :: normalelev_sec(:,:)
!normal wetted area of each cross section in the sub-basin's reservoir [m]
real, allocatable :: normalarea_sec(:,:)

!sediment transport
!settling velocity [m/s]
real, allocatable :: setvel(:)
!closest point above the water line at the left site (seen from upstream) of the cross section in the sub-basin's reservoir [m]
real, allocatable :: point1_bank(:,:)
!closest point above the water line at the right site (seen from upstream) of the cross section in the sub-basin's reservoir [m]
real, allocatable :: point2_bank(:,:)
!first point below the water line (seen from upstream) of the cross section in the sub-basin's reservoir [m**2]
real, allocatable :: point1_sub(:,:)
!last point below the water line (seen from upstream) of the cross section in the sub-basin's reservoir [m**2]
real, allocatable :: point2_sub(:,:)
!active layer area for each cross section in the sub-basin's reservoir [m**2]
real, allocatable :: area_actlay(:,:)
!top layer area for each cross section in the sub-basin's reservoir [m**2]
real, allocatable :: area_toplay(:,:)
!active layer volume for each cross section in the sub-basin's reservoir [m**3]
real, allocatable :: vol_actlay(:,:)
!top layer volume for each cross section in the sub-basin's reservoir [m**3]
real, allocatable :: vol_toplay(:,:)
!fractional sediment availability at the active layer for each cross section in the sub-basin's reservoir [ton]
real, allocatable :: frsedavailab(:,:)
!fractional sediment erosion at the active layer for each cross section in the sub-basin's reservoir [ton]
real, allocatable :: frerosion(:,:)
!fractional sediment deposition at the active layer for each cross section in the sub-basin's reservoir [ton]
real, allocatable :: frdeposition(:,:)
!fractional deposition of the material coming from the dam for each cross section in the sub-basin's reservoir [ton]
real, allocatable :: frretention(:,:)
!fractional suspention of the material retained at the top layer for each cross section in the sub-basin's reservoir [ton]
real, allocatable :: frsuspension(:,:)
!fractional bed load transport for each cross section in the sub-basin's reservoir [ton]
real, allocatable :: frbed_discharge(:,:)
!fractional suspended load transport for each cross section in the sub-basin's reservoir [ton]
real, allocatable :: frsusp_discharge(:,:)
!fractional sediment transport for each cross section in the sub-basin's reservoir [ton]
real, allocatable :: frtotal_discharge(:,:)
!total sediment erosion at the active layer for each cross section in the sub-basin's reservoir [ton]
real, allocatable :: erosion(:,:)
!total sediment deposition at the active layer for each cross section in the sub-basin's reservoir [ton]
real, allocatable :: deposition(:,:)
!total deposition of the material coming from the dam at the active layer for each cross section in the sub-basin's reservoir [ton]
real, allocatable :: retention(:,:)
!total suspention of the material retained at the top layer for each cross section in the sub-basin's reservoir [ton]
real, allocatable :: suspension(:,:)
!total sediment tranport at the active layer for each cross section in the sub-basin's reservoir [ton]
real, allocatable :: totalload(:,:)
!fractional bed carrying capacity (ton)
real, allocatable :: bed_frtransp(:,:)
!fractional suspended carrying capacity (ton)
real, allocatable :: susp_frtransp(:,:)
!fractional carrying capacity for total load (ton)
real, allocatable :: fr_capacity(:,:)
!depth change of bed sediment due to deposition or scour for each cross section in the sub-basin's reservoir [m]
real, allocatable :: dheight_sed(:,:)
!change in the area of bed sediment due to deposition or scour for each cross section in the sub-basin's reservoir [m**2]
real, allocatable :: darea_sed(:,:)
!change in the volume of bed sediment due to deposition or scour for each cross section in the sub-basin's reservoir [m**3]
real, allocatable :: dvol_sed(:,:)
!fractional sediment volume at the active layer for each cross section in the sub-basin's reservoir [m**3]
real, allocatable :: frvol_actlay(:,:,:)
!actual sediment volume at the active layer for each cross section in the sub-basin's reservoir [m**3]
real, allocatable :: totvol_actlay(:,:)
!sediment concentration for each cross section in the sub-basin's reservoir [g/l]
real, allocatable :: conc(:,:)
!fractional sediment volumetric concentration for each cross section in the sub-basin's reservoir [-]
real, allocatable :: frconc(:,:)
!sedimentation area for each cross section in the sub-basin's reservoir [m**2]
real, allocatable :: area_sedim(:,:)
!sedimentation volume for each cross section in the sub-basin's reservoir [m**3]
real, allocatable :: vol_sedim(:,:)
!initial sediment deposition in the sub-basin's reservoir [m**3]
real, allocatable :: volbed0(:)
!distance of the plunge point to the dam [m]
real, allocatable :: length_plunge(:)
!distance to the dam [m]
real, allocatable :: cumlength_sec(:,:)
!length of the reach represented by each cross section [m]
real, allocatable :: length_sec(:,:)
!mean diameter D50 related to the bed material at the active layer [m]
real, allocatable :: d50_actlay(:,:)
!D90 related to the bed material at the active layer [m]
real, allocatable :: d90_actlay(:,:)
!sediment entering sub-basin's reservoir for each particle size class k [tons/timestep]
real, allocatable :: frsedinflow(:,:,:)
!internal variables (temporarily described here) 
real, allocatable :: frvol_actlay0(:,:,:)
real, allocatable :: totvol_actlay0(:,:)


!reservoir bed elevation changes
!bed changes at the x-axis for each point in the cross section of the reservoir
!(from left to right, seen from upstream) [m]
real, allocatable :: x_sec(:,:,:)
!bed changes at the y-axis for each point in the cross section of the reservoir
!(from left to right, seen from upstream) [m]
real, allocatable :: y_sec(:,:,:)
!active layer elevation for each cross section in the sub-basin's reservoir [m]
real, allocatable :: y_actlay(:,:,:)
!active layer elevation for each cross section in the sub-basin's reservoir [m]
real, allocatable :: y_original(:,:,:)
!grain size distribution at the active layer for each cross section in the sub-basin's reservoir [-]
real, allocatable :: frac_actlay(:,:,:)
!grain size distribution at the top layer for each cross section in the sub-basin's reservoir [-]
real, allocatable :: frac_toplay(:,:,:)
!grain size distribution at the inactive layer for each cross section in the sub-basin's reservoir [-]
real, allocatable :: frac_comlay(:,:,:)
!grain size distribution at the top layer that remains in suspension for each cross section in the sub-basin's reservoir [-]
real, allocatable :: frac_susp(:,:,:)
!partial area of the active layer between two adjacent points for each cross section in the sub-basin's reservoir [m**2]
real, allocatable :: partarea_actlay(:,:,:)
!partial area of the top layer between two adjacent points for each cross section in the sub-basin's reservoir [m**2]
real, allocatable :: partarea_toplay(:,:,:)
!weighting factor for the calculation of changes at the active layer elevation for each cross section in the sub-basin's reservoir [-]
real, allocatable :: weightfac_actlay(:,:,:)
!weighting factor for the calculation of changes at the top layer elevation for each cross section in the sub-basin's reservoir [-]
real, allocatable :: weightfac_toplay(:,:,:)
!bed elevation at the beginning of the actual simulation timestep for each cross section in the sub-basin's reservoir [m]
real, allocatable :: y_laststep(:,:,:)
!water level used to attenuate elevation changes caused by erosion [m] 
real, allocatable :: erosion_level(:,:)
!Point of the cross section in the sub-basin's reservoir that identifies the beginning of main channel (from left to right, view from upstream side) [-]
integer, allocatable :: pt1(:,:)
!Point of the cross section in the sub-basin's reservoir that identifies the end of main channel (from left to right, view from upstream side) [-]
integer, allocatable :: pt2(:,:)
!Point located before pt1 that defines the minimum reach available to be eroded (from left to right, view from upstream side) [-] (used to minimize sediment erosion concentrated on the main channel)
integer, allocatable :: pt3(:,:)
!Point located after pt2 that defines the minimum reach available to be eroded  (from left to right, view from upstream side) [-] (used to minimize sediment erosion concentrated on the main channel)
integer, allocatable :: pt4(:,:)
!Initial location of the delta plunge point along the reservoir cross sections [-] (location varies acording to sedimento erosion/deposition in the reservoi)
integer, allocatable :: pt_long0(:)
!Location of the delta plunge point along the reservoir cross sections [-] 
integer, allocatable :: pt_long(:)
!Side slope at the beginning of main channel defined by pt1 [m/m]
real, allocatable :: sideslope_pt1(:,:)
!Side slope at the end of main channel defined by pt2 [m/m]
real, allocatable :: sideslope_pt2(:,:)
!Bed slope at the end of the delta plunge point in the sub-basin's reservoir [m/m]
real, allocatable :: slope_long(:)

!printable variables
!reservoir water balance
real, allocatable :: daydamelevact(:,:)			!described above as damelevact
real, allocatable :: daydamareaact(:,:)			!described above as damareaact
real, allocatable :: dayelev_bat(:,:,:)			!described above as elev_bat
real, allocatable :: dayarea_bat(:,:,:)			!described above as area_bat
real, allocatable :: dayvol_bat(:,:,:)			!described above as vol_bat
!reservoir sediment balance
real, allocatable :: daydepth_sec(:,:,:)		!described above as depth_sec
real, allocatable :: daywatelev_sec(:,:,:)		!described above as watelev_sec
real, allocatable :: dayarea_sec(:,:,:)			!described above as area_sec
real, allocatable :: daytopwidth_sec(:,:,:)		!described above as topwidth_sec
real, allocatable :: dayenergslope_sec(:,:,:)	!described above as energslope_sec
real, allocatable :: dayhydrad_sec(:,:,:)		!described above as hydrad_sec
real, allocatable :: daymeanvel_sec(:,:,:)		!described above as meanvel_sec
real, allocatable :: daydischarge_sec(:,:,:)	!described above as discharge_sec
real, allocatable :: dayminelev_sec(:,:,:)		!described above as minelev_sec
real, allocatable :: dayy_sec(:,:,:,:)			!described above as y_sec
real, allocatable :: daycumsed(:,:)				!described above as cum_sedimentation
real, allocatable :: dayfrsediment_out(:,:,:)	!described above as frsediment_out

end module reservoir_h


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

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
  del=(b-a)/dble(it)                                 !This is the spacing of the points to be added.
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
REAL, PARAMETER :: EPS=1.0d-1
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