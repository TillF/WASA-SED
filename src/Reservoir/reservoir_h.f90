! Till: computationally irrelevant: outcommented unused variables
! 2012-09-14
!
! Till: computationally irrelevant: minor changes to improve compiler compatibility
! 2011-04-29
!
!The module file contains two modules: reservoir_h and lake_h
!

module reservoir_h
save

INTEGER :: nxsection_res=200 !number of cross sections at the strategic reservoirs
INTEGER :: npointsxsect=200  !number of points (x,y) along the cross sections



! Options to run the WASA Model without river and hillslope modules
INTEGER :: reservoir_balance	!option to read (0) or calculate (1) outflow discharges and reservoir levels
INTEGER :: reservoir_check		!option to run the WASA Model with river and hillslope modules (0) or without them (1)
INTEGER :: reservoir_print		!option to print the output files within the timestep (0) or once a year (1)
!
!!Reservoir water balance arrays
!!simulation timestep
integer :: step
!hour (on day)
integer :: hour
! temporary array for reading of data from intake.dat
real, allocatable :: r_qintake(:)
! number of columns in intake.dat
integer :: no_col_intake
! for each subasin position in intake.dat
integer, allocatable :: corr_column_intakes(:)
! does intake.dat exist for specific subasin (and contain information)?
logical, allocatable :: f_intake_obs(:)
! flag indicating if a specific subbasin contains a strategic reservoir; FALSE: no reservoir; TRUE: contains a reservoir
logical, allocatable :: res_flag(:)
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
!!maximum water depth used for the overflow calculation (m)
real, allocatable :: hmax(:)
!water volume above the spillway elevation on the day before (m3/timestep)
real, allocatable :: volume_last(:)
!spillway overflow on the day before (m3/timestep)
real, allocatable :: outflow_last(:)
!tp not used
!effective water extraction (consumptive use) from the Sub-basin's reservoir [m**3]
!real, allocatable :: damex(:,:)
!initial storage capacity in the sub-basin's reservoir [read as 10**3 m**3, later converted into 10**6 m**3]
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
!actual stored volume in the subbasin's reservoir [m^3 and 10^6 m**3]
real, allocatable :: volact(:,:)
!precipitation into the subbasin's reservoir [m**3]
real, allocatable :: precdam(:,:)
!evaporation from the subbasin's reservoir [m**3]
real, allocatable :: etdam(:,:)
! tp TODO never used
!!evaporation from the subbasin's reservoir [mm]
!real, allocatable :: evapdam(:,:)
!!infiltration losses from the subbasin's reservoir [m**3]
!real, allocatable :: infdam(:,:)
!initial maximum reservoir area in subbasin [ha]
real, allocatable :: maxdamarea(:)
!actual reservoir area [m**2]
real, allocatable :: damareaact(:)
!!volume=k*Hv**alpha (Volume/heigth) relationship parameters
real, allocatable :: alpha_over(:)
real, allocatable :: k_over(:)
!!area=a*Vol**b (Volume/area) relationship parameters
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
!inflow discharge into the subbasin's reservoir [m**3/s and m³]
real, allocatable :: qinflow(:,:)
!overflow discharge out the sub-basin's reservoir [m**3/s]
real, allocatable :: overflow(:,:)
!controlled outflow by intake device in the sub-basin's reservoir [m**3/s]
real, allocatable :: qintake(:,:)
!controlled outflow relased through the bottom outlets of the sub-basin's reservoir [m**3/s]
real, allocatable :: qbottom(:,:)
!actual withdrawal from the sub-basin's reservoir (e.g. for irrigation; in the model not further used) [m**3/s]
real, allocatable :: withdraw_out(:,:)
!storage capacity in the subbasin's reservoir [m**3 and 10**6 m**3]
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
!
!!Reservoir sedimentation (semres) arrays
!!general parameters
!external IDs of cross sections in the sub-basin's reservoir
integer, allocatable :: id_sec_extern(:,:)
!number of cross sections in the sub-basin's reservoir
integer, allocatable :: nbrsec(:)
!number of points at the cross section of the reservoir in the sub-basin
integer, allocatable :: npoints(:,:)
!!number of cases for the calculation of bed geometry changes between two adjacent points
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
!!values at the x-axis for each point in the cross section of the reservoir
!(from left to right, seen from upstream) [m]
real, allocatable :: x_sec0(:,:,:)
!!values at the y-axis for each point in the cross section of the reservoir
!(from left to right, seen from upstream) [m]
real, allocatable :: y_sec0(:,:,:)
!!suspended sediment in the subbasin's reservoir [ton/timestep]
!real, allocatable :: sed_susp(:,:)
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
!!minimum sediment concentration in the sub-basin's reservoir [mg/l]
!real, allocatable :: min_conc(:)
!!initial sediment concentration in the sub-basin's reservoir [mg/l]
!real, allocatable :: sed_conc0(:)
!!wet bulk density for sediments in the sub-basin's reservoir [ton/m**3]
!real, allocatable :: wet_dens(:)
!! sediment release from reservoir for each particle size class k (tons/timestep)
real, allocatable :: res_sediment_out(:,:)
!!grain size distribution of the sediment inflow into the subbasin's reservoir [-]
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
!!amount of correlations between water inflow discharge and incoming sediment into the subbasin's reservoir [-]
!integer, allocatable :: numbeqs(:)
!!water discharge that represents the upper limit of applicability of each correlation
!!between water inflow discharge and incoming sediment into the subbasin's reservoir [-]
!real, allocatable :: Q_max(:,:)
!!Qsed=a*Q**b (water input x sediment input) relationship parameters [-]
!real, allocatable :: param_a(:,:)
!real, allocatable :: param_b(:,:)
!!number of size distribution of the incoming sediment into the subbasin's reservoir [-]
!integer, allocatable :: nbsizedist(:)
!!water discharge value with a given size distribution of the incoming sediment into the subbasin's reservoir [-]
!real, allocatable :: Q_refer(:,:)
!!computed size distribution of the incoming sediment into the subbasin's reservoir [-]
!real, allocatable :: perc(:,:,:)
!sediment inflow into the subbasin's reservoir related to grain size g [ton(timestep]
real, allocatable :: sedinflow_g(:,:,:)
!sediment inflow into the subbasin's reservoir related to grain size g [ton(timestep]
real, allocatable :: sedoutflow_g(:,:,:)
!
!!hydraulic calculations
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
!
!!sediment transport
!settling velocity [m/s]
real, allocatable :: setvel(:)
!!closest point above the water line at the left site (seen from upstream) of the cross section in the sub-basin's reservoir [m]
!!real, allocatable :: point1_bank(:,:)
!!closest point above the water line at the right site (seen from upstream) of the cross section in the sub-basin's reservoir [m]
!real, allocatable :: point2_bank(:,:)
!!first point below the water line (seen from upstream) of the cross section in the sub-basin's reservoir [m**2]
!real, allocatable :: point1_sub(:,:)
!!last point below the water line (seen from upstream) of the cross section in the sub-basin's reservoir [m**2]
!real, allocatable :: point2_sub(:,:)
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
!!fractional suspention of the material retained at the top layer for each cross section in the sub-basin's reservoir [ton]
!real, allocatable :: frsuspension(:,:)
!!fractional bed load transport for each cross section in the sub-basin's reservoir [ton]
!real, allocatable :: frbed_discharge(:,:)
!!fractional suspended load transport for each cross section in the sub-basin's reservoir [ton]
!real, allocatable :: frsusp_discharge(:,:)
!fractional sediment transport for each cross section in the sub-basin's reservoir [ton]
real, allocatable :: frtotal_discharge(:,:)
!total sediment erosion at the active layer for each cross section in the sub-basin's reservoir [ton]
real, allocatable :: erosion(:,:)
!total sediment deposition at the active layer for each cross section in the sub-basin's reservoir [ton]
real, allocatable :: deposition(:,:)
!total deposition of the material coming from the dam at the active layer for each cross section in the sub-basin's reservoir [ton]
real, allocatable :: retention(:,:)
!!total suspention of the material retained at the top layer for each cross section in the sub-basin's reservoir [ton]
!real, allocatable :: suspension(:,:)
!total sediment tranport at the active layer for each cross section in the sub-basin's reservoir [ton]
real, allocatable :: totalload(:,:)
!fractional bed carrying capacity (ton)
real, allocatable :: bed_frtransp(:,:)
!fractional suspended carrying capacity (ton)
real, allocatable :: susp_frtransp(:,:)
!fractional carrying capacity for total load (ton)
real, allocatable :: fr_capacity(:,:)
!!depth change of bed sediment due to deposition or scour for each cross section in the sub-basin's reservoir [m]
!real, allocatable :: dheight_sed(:,:)
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
!!internal variables (temporarily described here)
real, allocatable :: frvol_actlay0(:,:,:)
real, allocatable :: totvol_actlay0(:,:)
!
!
!!reservoir bed elevation changes
!!bed changes at the x-axis for each point in the cross section of the reservoir
!(from left to right, seen from upstream) [m]
real, allocatable :: x_sec(:,:,:)
!X-coordinate of bed changes for each point in the cross section of the reservoir (from left to right, seen from upstream)
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
!!weighting factor for the calculation of changes at the active layer elevation for each cross section in the sub-basin's reservoir [-]
!real, allocatable :: weightfac_actlay(:,:,:)
!!weighting factor for the calculation of changes at the top layer elevation for each cross section in the sub-basin's reservoir [-]
!real, allocatable :: weightfac_toplay(:,:,:)
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
!
!!printable variables
!!reservoir water balance
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

integer :: m
end module reservoir_h
