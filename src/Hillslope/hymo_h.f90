module hymo_h
    use common_h
    save

    ! Code converted using TO_F90 by Alan Miller
    ! Date: 2005-06-30  Time: 13:59:19

    ! WATERSHED / Sub-basin PARAMETERS

    ! IDs of cells in watershed (changed to model internal number)
    !Allocatable      integer id_subbas_intern(subasin)
    ! area of sub-basin (square kilometer)
    !      real area(subasin)

    real, allocatable :: area(:)

    ! SUBBASIN / GRID CELL PARAMETERS
    ! IDs of subbasins (external and internal IDs)
    !Allocatable integer id_subbas_extern(subasin)
    integer, allocatable ::  id_subbas_extern(:),id_subbas_intern(:)
    ! number of landscape units within each subbasin
    !Allocatable     integer nbr_lu(subasin)
    integer, allocatable ::  nbr_lu(:)
    ! IDs of landscape units in cell
    !Allocatable       integer id_lu_intern(maxsoter,subasin)
    integer, allocatable ::  id_lu_intern(:,:)
    ! fraction of LU in subbasin
    !Allocatable       real frac_lu(maxsoter,subasin)
    real, allocatable :: frac_lu(:,:)
    ! scaling factor for hydraulic conductivity in infiltration routine
    !Allocatable       real kfkorrc(subasin)
    real, allocatable :: kfkorrc(:)
    ! scaled interception capacity (at the scale of cells for interpolated daily rainfall)
    !Allocatable       real intcfc(subasin)
    real, allocatable :: intcfc(:)


    ! Landscape UNIT PARAMETERS
    ! IDs of LUs
    !Allocatable      integer id_lu_extern(nsoter)
    integer, allocatable :: id_lu_extern(:)
    ! number of terrain components in landscape unit
    !Allocatable      integer nbrterrain(nsoter)
    integer, allocatable ::  nbrterrain(:)
    ! IDs of terrain component in landscape unit
    !Allocatable      integer id_terrain_intern(maxterrain,nsoter)
    integer, allocatable ::   id_terrain_intern(:,:)
    ! position of landscape unit in landscape

    !!Allocatable      integer possoter(nsoter)					!Till: is this used at all?
    !integer, allocatable ::   possoter(:)

    ! hydr. conductivity of bedrock in landscape unit
    !Allocatable      real kfsu(nsoter)
    real, allocatable :: kfsu(:)
    ! mean slope length in landscape unit
    !Allocatable      real slength(nsoter)
    real, allocatable ::  slength(:)
    ! mean maximum depth of soil zone in SU
    !Allocatable      real meandep(nsoter)
    real, allocatable ::   meandep(:)
    ! maximum depth of alluvial soil zone
    !Allocatable      real maxdep(nsoter)
    real, allocatable ::    maxdep(:)
    ! depth of river bed below terrain surface (mm)
    !Allocatable      real riverbed(nsoter)
    real, allocatable ::    riverbed(:)
    ! flag for landscape unit with groundwater table
    !Allocatable      integer gw_flag(nsoter)
    integer, allocatable ::   gw_flag(:)
    ! initial depth of groundwater below surface (mm)
    !Allocatable      real gw_dist(nsoter)
    real, allocatable ::     gw_dist(:)
    ! storage coefficient for groundwater outflow (days)
    !Allocatable      real gw_delay(nsoter)
    real, allocatable ::      gw_delay(:)
    ! ordered positon of terrain component in landscape unit
    ! (highest first)
    !Allocatable      integer orderterrain(maxterrain,nsoter)
    integer, allocatable :: orderterrain(:,:)

    ! TERRAIN COMPONENT PARAMETERS
    ! IDs of terrain components
    !Allocatable      integer id_terrain_extern(nterrain)
    integer, allocatable :: id_terrain_extern(:)
    ! slope of terrain component (%)
    !Allocatable     real slope(nterrain)
    real, allocatable :: slope(:)
    ! fraction of terrain component in landscape unit
    !Allocatable      real fracterrain(nterrain)
    real, allocatable ::  fracterrain(:)
    ! positon of terrain component in landscape unit
    !Allocatable      real posterrain(nterrain)
    real, allocatable ::  posterrain(:)
    ! number of Soil-Vegetation components (SVCs) in each TC of each Sub-basin and landscape Unit
    !Allocatable      integer nbr_svc(nmunsutc)
    integer, pointer :: nbr_svc(:)
    ! IDs of soils of each SVC in each TC of each Sub-basin and landscape Unit
    !Allocatable      integer id_soil_intern(maxsoil,nmunsutc)
    integer, pointer :: id_soil_intern(:,:)
    ! IDs of landuse units of each SVC in each TC of each Sub-basin and landscape Unit
    !Allocatable      integer id_veg_intern(maxsoil,nmunsutc)
    integer, pointer ::  id_veg_intern(:,:)
    ! fractions of SVCs in each TC
    !Allocatable      real frac_scv(maxsoil,nmunsutc)
    real, pointer :: frac_svc(:,:)
    ! fraction of impermeable (rock) area in each TC
    !Allocatable      real rocky(nmunsutc)
    real, pointer :: rocky(:)
    ! IDs of all Subbasin-LU-TC-combinations
    !Allocatable      real tcallid(subasin,maxsoter,maxterrain)
    integer, allocatable :: tcallid(:,:,:)
    !fraction of sheetflow that gets concentrated within the TC (and vice-versa)
    real, allocatable :: frac_diff2conc(:), frac_conc2diff(:)


    ! SOIL-VEGETATION COMPONENT PARAMETERS
    ! Till thickness of horizons (lowest horizon may be different from that of profiles)
    !Allocatable      real horiz_thickness(nmunsutc,maxsoil,maxhori)
    real, pointer :: horiz_thickness(:,:,:)
    ! lowest horizon to which roots go down
    !Allocatable      integer svcrooth(nmunsutc,maxsoil)
    integer, pointer ::  svcrooth(:,:)
    ! flag indicating bedrock or no bedrock (1/0)
    !Allocatable      integer svcbedr(nmunsutc,maxsoil)
    integer, pointer ::   svcbedr(:,:)
    ! water content at permanent wilting point (VOL%) (variable depending on
    ! plant characteristics, for each Soil-Vegetation unit)
    !Allocatable      real pwpsc(nmunsutc,maxsoil,maxhori)
    real, pointer ::  pwpsc(:,:,:)
    ! distribution of saturated soil water among horizons
    !Allocatable      real horiths(nmunsutc,maxsoil,maxhori)
    real, pointer ::  horiths(:,:,:)
    ! saturated water content integrated for profile (mm)
    !Allocatable      real thsprof(nmunsutc,maxsoil)
    real, pointer ::   thsprof(:,:)
    ! saturated soil water content distribution (mm)
    !Allocatable      real tctheta_s(5,2,ntcinst,maxsoil)
    real, pointer ::   tctheta_s(:,:,:,:)

    ! SOIL PARAMETERS
    ! IDs of soil components
    !Allocatable      integer id_soil_extern(nsoil)
    integer, allocatable :: id_soil_extern(:)
    ! number of horizons in soil
    !Allocatable      integer nbrhori(nsoil)
    integer, allocatable ::nbrhori(:)
    ! residual soil water content (VOL%)
    !Allocatable      real thetar(nsoil,maxhori)
    real, allocatable :: thetar(:,:)
    ! water content at permanent wilting point (VOL fraction) (standard 15000 cm suction)
    !Allocatable      real soilpwp(nsoil,maxhori)
    real, allocatable :: soilpwp(:,:)
    ! saturated soil water content (VOL%)
    !Allocatable      real thetas(nsoil,maxhori)
    real, allocatable ::  thetas(:,:)
    ! usable field capacity (nFK) (VOL%)
    !Allocatable      real soilnfk(nsoil,maxhori)
    real, allocatable ::   soilnfk(:,:)
    !
    ! field capacity (FK) (VOL%)
    !two FC values Saugspannung 63 hPa (also pF=1.8, FK1.8), and 316hPa (pF=2.6, FK2.5)
    ! within the WAVES-Project, it was not clear which value would be more appropriate,
    ! so both were tested.
    !Allocatable      real soilfc(nsoil,maxhori)
    real, allocatable :: soilfc(:,:)
    ! field capacity (FK63) (VOL%)
    !Allocatable      real soilfc63(nsoil,maxhori)
    real, allocatable ::  soilfc63(:,:)

    ! does bedrock occur below deepest horizon of profile ? (0/1)
    !Allocatable      integer bedrock(nsoil)
    integer, allocatable :: bedrock(:)
    ! is this an alluvial soil ? (0/1)
    integer, allocatable :: alluvial_flag(:)
    ! saturated hydraulic conductivity (mm/day)
    !Allocatable      real k_sat(nsoil,maxhori)
    real, allocatable :: k_sat(:,:)
    ! suction at the wetting front
    !Allocatable      real saug(nsoil,maxhori)
    real, allocatable ::  saug(:,:)
    ! Brooks-Corey pore-size index (-) (transformed into Van-Genuchten m-parameter)
    !Allocatable      real poresz(nsoil,maxhori)
    real, allocatable ::   poresz(:,:)
    ! Van-Genuchten n-parameter ( = pore-size index +1.) (-)
    !Allocatable      real porem(nsoil,maxhori)
    real, allocatable ::porem(:,:)
    ! Bubbling pressure (cm)
    !Allocatable      real bubble(nsoil,maxhori)
    real, allocatable ::bubble(:,:)
    ! Fraction of coarse fragments (-)
    !Allocatable      real coarse(nsoil,maxhori)
    real, allocatable :: coarse(:,:)
    ! Flag for soil structure (-)
    !Allocatable      real shrink(nsoil,maxhori)
    real, allocatable ::  shrink(:,:)


    !-------------------------------------------------------
    ! VEGETATION PARAMETER
    ! ID of vegetation unit
    !Allocatable      integer id_veg_extern(nveg)
    integer, pointer :: id_veg_extern(:)
    ! height (m) (monthly)
    !Allocatable       real height(nveg,4)
    real, pointer ::height(:,:)
    ! root depth (mm) (monthly)
    !Allocatable       real rootdep(nveg,4)
    real, pointer :: rootdep(:,:)
    ! leaf area index (-) (monthly)
    !Allocatable       real lai(nveg,4)
    real, pointer :: lai(:,:)
    ! albedo (-) (monthly)
    !Allocatable       real alb(nveg,4)
    real, pointer ::  alb(:,:)
    ! stomata resistance without water stress (s/m)
    !Allocatable       real resist(nveg)
    real, allocatable ::   resist(:)
    ! suction threshold for water stress effect on resistance
    ! (begin of stomata closure)
    !Allocatable       real wstressmin(nveg)
    real, allocatable ::    wstressmin(:)
    ! suction threshold for water stress effect on resistance
    ! (total closure of stomata - wilting point)
    !Allocatable       real wstressmax(nveg)
    real, allocatable ::    wstressmax(:)
    ! four key points in time for temporal distribution of vegetation
    ! characteristics within year
    ! (begin/end of rainy period)
    ! specific for each Sub-basin and year
    !      integer period(4,subasin,200)
    integer, pointer :: period(:,:)	! four key nodes in time for temporal vegetation dynamics within year (index: subbasin,(1:4)*simulation_year,)


    ! daily mean LAI (m²/m²)
    !Allocatable      real laimun(366,subasin)
    real, allocatable :: laimun (:,:)
    !Allocatable      real lai_c(366,subasin)
    REAL, allocatable :: laisu(:)
    !Allocatable      real laitc(nmunsutc)
    real, pointer ::    laitc(:)

    !-------------------------------------------------------
    ! WATER BALANCE VARIABLES

    ! WATER BALANCE VARIABLES      MUNICIP/WATERSHED SCALE
    ! soil moisture (mm)
    !Allocatable      real soilm(366,subasin)
    real, allocatable ::  soilm(:,:)
    ! soil moisture first meter (mm)
    !Allocatable       real soilmroot(366,subasin)
    real, allocatable ::  soilmroot(:,:)
    ! daily actual evapotranspiration (mm/day)
    !Allocatable       real aet(366,subasin)
    real, allocatable ::  aet(:,:)

    real, pointer ::  aet_t(:,:,:)				!actual evapotrans for each day and timestep(366,nt,subasin)
    real, pointer ::  hortflow_t(:,:,:)			!horton flow for each day and timestep(366,nt,subasin)
    real, pointer ::  subflow_t(:,:,:)	!subsurface runoff for each day and timestep(366,nt,subasin)
    real, pointer ::  ovflow_t(:,:,:) !total overland flow for each day and timestep(366,nt,subasin)
    real, pointer ::  deep_gw_discharge_t(:,:,:)		!groundwater discharge for each day and timestep(366,nt,subasin)
    real, pointer ::  pet_t(:,:,:)				!potential evapotrans for each day and timestep(366,nt,subasin)
    real, pointer ::  deep_gw_recharge_t(:,:,:)		!groundwater rescharge for each day and timestep(366,nt,subasin)
    real, pointer ::  gw_loss_t(:,:,:)				!ground water loss (deep percolation in LUs with no ground water flag) for each day and timestep(366,nt,subasin)
    real, pointer ::  river_infiltration_t(:,:,:)	!infiltration into riverbed, loss from model domain
    real, pointer ::  riverflow_t(:,:,:)	!Till: flow in the river in m3/s for each day and timestep(366,nt,subasin)
    real, pointer :: irri_supply_record(:,:,:)   !Paul irrigation water that each subbasin receives for each day and timestep(366,nt,subasin)
    real, pointer :: irri_abstraction_record(:,:,:)  !Paul irrigation water that is abstracted from each subbasin for each day and timestep(366,nt,subasin)

    ! daily soil evaporation (mm/day)
    !Allocatable       real soilet(366,subasin)
    real, allocatable :: soilet(:,:)
    ! daily interception storage evapotranspiration (mm/day)
    !Allocatable       real intc(366,subasin)
    real, allocatable ::  intc(:,:)
    ! total surface runoff (m³/d)
    !Allocatable       real ovflow(366,subasin)
    real, allocatable ::  ovflow(:,:)
    ! horton overland flow (m³/d)
    !Allocatable      real hortflow(366,subasin)
    real, allocatable :: hortflow(:,:)
    ! total subsurface runoff (m³/d)
    !Allocatable       real subflow(366,subasin)
    real, allocatable ::  subflow(:,:)

    real, allocatable ::   deep_gw_discharge(:,:)	!groundwater discharge into river (366,subasin)
    real, allocatable ::   gw_loss(:,:)				!ground water loss (deep percolation in LUs with no ground water flag)(366,subasin)


    ! groundwater recharge (percolation below root zone)
    !Allocatable       real gw_recharge(366,subasin)
    real, allocatable ::   gw_recharge(:,:)
    !! deep groundwater recharge (loss from model / into lin. GW storage)
    !!Allocatable       real deepgw_r(366,subasin)
    !real, allocatable ::    deepgw_r(:,:)
    ! river flow before acudes (m3)
    !Allocatable       real qgen(366,subasin)
    real, allocatable :: qgen(:,:)
    !Till: contribution of each subbasin to the river in m3/s for each day (366,subasin)
    real, allocatable ::  water_subbasin(:,:)
    !Till: contribution of each subbasin to the river in m3/s for each day and timestep(366,nt,subasin)
    real, allocatable ::  water_subbasin_t(:,:,:)

    ! losses in river network by evaporation
    !Allocatable       real qloss(366,subasin)
    real, allocatable ::   qloss(:,:)
    ! fraction of saturated area (-)
    !Allocatable       real sofarea(366,subasin)
    real, allocatable ::  sofarea(:,:)

    !REAL :: seaflow(366)	! river discharge to ocean (m3)	!Till: never used
    ! water availability
    !Allocatable       real avail_all(366,subasin),avail_ac(366,subasin)
    real, allocatable :: avail_all(:,:), avail_ac(:,:)
    !Till: pre-specified outflow of subbasins (optionally read from file) [m³/s]
    real, allocatable ::  pre_subbas_outflow(:,:,:)
    !Till: holds corresponding columns of input files to be related to internal numbering of subbasins
    integer, pointer :: corr_column_pre_subbas_outflow(:)



    ! WATER BALANCE VARIABLES      LANDSCAPE UNIT SCALE
    ! soil moisture (mm)
    REAL, allocatable  :: soilmsu(:)
    ! soil moisture first meter (mm)
    REAL, allocatable  :: soilmrootsu(:)
    ! daily actual evapotranspiration (mm/day)
    REAL, allocatable  :: aetsu(:)
    ! daily soil evaporation (mm/day)
    REAL, allocatable  :: soiletsu(:)
    ! daily interception storage evaporation (mm/day)
    REAL, allocatable  :: intcsu(:)
    ! total surface runoff (m³/d) for all LUs of current subbasin
    REAL, allocatable  :: qsurf_lu(:)
    ! total subsurface runoff (m³/d) for all LUs of current subbasin
    REAL, allocatable  :: qsub_lu(:)
    ! groundwater recharge (percolation from root zone)
    REAL, allocatable  :: gwrsu(:)
    ! deep groundwater recharge (loss from model / into lin. GW storage)
    REAL, allocatable  :: deepgwrsu(:)
    ! horton overland flow
    REAL, allocatable  :: hortsu(:)
    ! actual deep groundwater storage
    !Allocatable       real deepgw(subasin,:)
    real, allocatable ::    deepgw(:,:)


    ! WATER BALANCE VARIABLES      TERRAIN COMPONENT SCALE
    ! average soil moisture of every terrain component (mm)
    !Allocatable      real soilwater(366,nmunsutc)
    real, pointer :: soilwater(:,:)
    
    ! soil moisture in every horizon of each SVC in each TC in each LU in each subbasin (mm)
    !Allocatable      real horithact(nmunsutc,maxsoil,maxhori)
    real, pointer ::  horithact(:,:,:)
    ! lateral subsurface runoff to be redistributed between SVCs (m**3)
    !Allocatable      real latred(nmunsutc,maxhori*3)
    real, allocatable ::   latred(:,:)
    ! interception storage of each SVC in each TC (mm)
    !Allocatable      real intercept(nmunsutc,maxsoil)
    real, pointer ::   intercept(:,:)
    ! interception evaporation daily, proportionally distributed among hours (mm)
    !AllocatableREAL :: intcept_mem(maxterrain,maxsoil,24)
    real, allocatable :: intcept_mem (:,:,:)
    ! reduction of evaporation due to interception evaporation for hourly version (-)
    !Allocatalbe REAL :: aet_red_mem(maxterrain,maxsoil)
    real, allocatable :: aet_red_mem(:,:)

    ! saturated fraction of each SVC in each TC (-)
    !Allocatable      real frac_sat(nmunsutc,maxsoil)
    real, pointer ::   frac_sat(:,:)
    ! average real evapotranspiration of every terrain component (mm)
    !Allocatable      real aettc(nmunsutc)
    real, pointer :: aettc(:)
    ! average soil evaporation of each terrain component (mm)
    !Allocatable      real soilettc(nmunsutc)
    real, pointer :: soilettc(:)
    ! average interception evaporation of every terrain component (mm)
    !Allocatable      real intctc(nmunsutc)
    real, pointer :: intctc(:)
    ! horton overland flow
    !Allocatable      real horttc(nmunsutc)
    real, pointer :: horttc(:)
    ! groundwater recharge (percolation)
    !Allocatable      real gwrtc(nmunsutc)
    real, pointer ::  gwrtc(:)
    ! deep groundwater recharge
    !Allocatable      real deepgwrtc(nmunsutc)
    real, pointer ::   deepgwrtc(:)
    ! actual transpiration of each SVC (only plants)
    REAL , allocatable :: aet1sc(:,:)
    ! soil evaporation of each SVC
    REAL , allocatable :: soilet1sc(:,:)


    !Conrad/Till: TC-wise output of soil moisture, surface flow and sediment output
    real, allocatable :: meandepth_tc(:,:,:)!mean soil depths of TC instances
    real, allocatable :: theta_tc(:,:,:)		!theta of tcs [%]			(day,timestep,ntcinst)
    real, allocatable :: surfflow_tc(:,:,:)		!surface runoff of TCs [mm] (day,timestep,ntcinst)
    real, allocatable :: sedout_tc(:,:,:)		!sediment output of TCs [t/km2] (day,timestep,ntcinst)

    real, allocatable :: sedout_lu(:,:,:)       !sediment output of TCs [t/km2] (day,timestep,ntcinst)


    !-----------------------------------------------------
    ! IRRIGATION VARIABLES
    REAL, allocatable :: frac_irr_sub(:)  ! Anteil der bewässerten Flächen innerhalb der Subbasins
    INTEGER, allocatable :: svc_irr(:)  ! irrigation variable svc_dat
    INTEGER, allocatable :: sub_source(:)           !for reading irri.dat
    INTEGER, allocatable :: sub_receiver(:)         !for reading irri.dat
    CHARACTER(len=12), allocatable :: irri_rule(:)  !for reading irri.dat
    REAL, allocatable :: irri_rate_res(:,:,:)               !for reading irri.dat, rates for reservoir abstraction
    REAL, allocatable :: irri_rate_lake(:,:,:)               !for reading irri.dat, rates for lake abstraction
    REAL, allocatable :: irri_rate_riv(:,:,:)               !for reading irri.dat, rates for river abstraction
    REAL, allocatable :: irri_rate_gw(:,:,:)               !for reading irri.dat, rates for groundwater abstraction
    REAL, allocatable :: irri_rate_ext(:,:)               !for reading irri.dat, rates for external abstraction
    REAL, allocatable :: cwd_gw(:,:,:)                     ! for storing coefficients for crop water demand
    REAL, allocatable :: cwd_ext(:,:)                ! for storing coefficients for crop water demand
    REAL, allocatable :: cwd_riv(:,:,:)                     ! for storing coefficients for crop water demand
    REAL, allocatable :: cwd_res(:,:,:)                     ! for storing coefficients for crop water demand
    REAL, allocatable :: cwd_lake(:,:,:)                     ! for storing coefficients for crop water demand
    CHARACTER(len=11), allocatable :: irri_source(:)    !for reading irri.dat
    INTEGER ::  nbr_irri_records !total number of vavlid irrigation records in irri.dat
    REAL, allocatable :: irri_supply(:) !stores the amout of irrigation water each subbasin recieves for each timestep  !allocated with dimension subasin in readhymo
    REAL, allocatable :: irri_abstraction(:) !stores the amout of irrigation water that is taken from each subbasin  for each timestep  !allocated with dimension subasin in readhymo
    INTEGER, pointer :: seasonality_irri(:,:)
    !INTEGER, pointer :: seasonality_irri_res(:,:)   !Paul 09.11.2020 First version with different seasonality depending on irri_source. For this option were 5 seasons.dat files necesssary. Simplified to one seasonality per sreceiver basin
    !INTEGER, pointer :: seasonality_irri_lake(:,:)
    !INTEGER, pointer :: seasonality_irri_riv(:,:)
    !INTEGER, pointer :: seasonality_irri_gw(:,:)
    !INTEGER, pointer :: seasonality_irri_ext(:,:)
    REAL, allocatable :: loss_gw(:)
    REAL, allocatable :: loss_riv(:)
    REAL, allocatable :: loss_res(:)
    REAL, allocatable :: loss_lake(:)
    REAL, allocatable :: loss_ext(:)
    REAL :: test_supply             ! Test variables: DELETE!
    REAL :: testPAUL
    REAL :: Paul
    INTEGER :: testcounter


    !Till: these are all output variables that are currently not used
    !! horton overland flow of each SVC
    !REAL , allocatable :: hortsc(:,:)
    !! saturated area of each SVC, relative to total TC area
    !REAL, allocatable  :: sat_area_of_sc(:,:)
    !! saturation excess overland flow of each SVC
    !REAL , allocatable :: sat_xs_overland_flow_sc(:,:)
    !! horton overland flow 2 of each SVC
    !REAL , allocatable :: hort2sc(:,:)
    !! groundwater recharge in each SVC
    !REAL , allocatable :: gwrsc(:,:)
    !! groundwater recharge in each SVC
    !REAL , allocatable :: deepgwrsc(:,:)
    !! soil moisture each SVC
    !REAL, allocatable  :: thsc(:,:,:)
    !! canopy resistance of each SVC
    !REAL, allocatable :: resistsc(:,:)
    !! actual transpiration of each SVC (only plants)
    !REAL , allocatable :: aetsc(:,:)
    !! total evapotranspiration of each SVC
    !REAL , allocatable :: aettotsc(:,:)
    !! interception evaporation of each SVC
    !REAL , allocatable :: intetsc(:,:)
    !! soil evaporation of each SVC
    !REAL , allocatable :: soiletsc(:,:)
    !! actual available field capacity of each SVC
    !REAL , allocatable :: nfksc(:,:)



    ! AGGREGATED VALUES IN TIME
    ! monthly average dummy array
    !      real mondummy(12,subasin)
    ! monthly average soil water content (mm)
    !      real monsoilw(12,subasin)
    ! annual average soil water content (mm)
    !      real annsoilw(subasin)
    ! monthly average soil water content (mm)
    !      real monsoilwroot(12,subasin)
    ! annual average soil water content (mm)
    !      real annsoilwroot(subasin)
    ! monthly average real evapotranspiration (mm)
    !      real monatp(12,subasin)
    ! annual average real evapotranspiration (mm)
    !      real annatp(subasin)
    ! annual average soil evaporation (mm)
    !      real annsoilet(subasin)
    ! monthly average soil evaporation (mm)
    !      real monsoilet(12,subasin)
    ! annual average LAI
    !      real annlai(subasin)
    ! monthly average LAI
    !      real monlai(12,subasin)
    ! monthly average interception (mm)
    !      real monintc(12,subasin)
    ! annual average interception (mm)
    !      real annintc(subasin)
    ! monthly average groundwater recharge (mm)
    !      real mongwr(12,subasin)
    ! annual average groundwater recharge (mm)
    !      real anngwr(subasin)
    ! monthly average deep groundwater recharge (mm)
    !      real mondeepgwr(12,subasin)
    ! annual average deep groundwater recharge (mm)
    !      real anndeepgwr(subasin)
    ! monthly average river runoff (m3/s) (before acudes)
    !      real monqgen(12,subasin)
    ! annual average river runoff (m3/s) (before acudes)
    !      real annqgen(subasin)
    ! monthly accumulated retention in reservoirs (m3)
    !      real monret(12,subasin)
    ! annualy accumulated retention in reservoirs (m3)
    !      real annret(subasin)
    ! annual horton overland flow (mm)
    !      real annhort(subasin)
    ! annual subsurface runoff (mm)
    !      real annsublat(subasin)

    !      real aveatp(subasin),aveintc(subasin),avesoilw(subasin)
    !      real avesoilwroot(subasin)
    !      real avegwr(subasin),aveqgen(subasin),avesoilet(subasin)
    !      real avelai(subasin),avehort(subasin),avesublat(subasin)
    !      real avedeepgwr(subasin)


    ! AGGREGATED VALUES IN SPACE
    ! data at mesoregion and state level
    !      real mesosoilw(2)  ,statesoilw(2)
    !      real mesosoilwroot(2)  ,statesoilwroot(2)
    !      real mesoatp(2)    ,stateatp(2)
    !      real mesolai(2)    ,statelai(2)
    !      real mesosoilet(2) ,statesoilet(2)
    !      real mesointc(2)   ,stateintc(2)
    !      real mesogwr(2)    ,stategwr(2)
    !      real mesodeepgwr(2),statedeepgwr(2)
    !      real mesoqgen(2) ,stateqgen(2)



    ! CALIBRATION PARAMETER
    ! calibration factor of hydraulic conductivity on daily basis
    REAL :: kfkorr
    REAL :: kfkorr_a,kfkorr_b	!coefficients for time-variant kfkorr (kfkorr=kfkorr*(kfkorr_a*1/daily_precip+kfkorr_b)
    REAL :: kfkorr_day			!kfkorr for current subbasin and day

    ! interception capacity per unit LAI (mm)
    REAL :: intcf
    ! type of interception routine (simple or modified bucket approach)
    INTEGER :: dointc

    ! fraction of ground water discharge routed directly into river (instead of subsurface flow into lowermost TC)
    REAL 	:: frac_direct_gw


 !   real, allocatable :: debug_out(:)		!for debugging purposes !remove
 !   real, allocatable :: debug_out2(:,:)		!for debugging purposes !remove
    integer :: debug_flag=0

contains

    FUNCTION allocate_hourly_array(f_flag)
        !allocate array for hourly output
        use params_h
        use common_h

        implicit none
        REAL, pointer :: allocate_hourly_array(:,:,:)
        LOGICAL, INTENT(IN)                  :: f_flag
        !CHARACTER(len=*), INTENT(IN)         :: name
        integer :: istate
        if (f_flag) then
            allocate(allocate_hourly_array(366,nt,subasin),STAT = istate)
            if (istate/=0) then
                write(*,'(A,i0,a)')'ERROR: Memory allocation error (',istate,') in hymo-module. Try disabling some hourly output.'
                stop
            end if
        else
            nullify(allocate_hourly_array)
        end if
    END FUNCTION allocate_hourly_array


	FUNCTION calc_seasonality2(subbas_id, year, julian_day ,seasonality_array, support_values)
        !replaces calc_seasonality
		!compute seasonality (value of current parameter for current timestep and subbasin) by interpolation between node_n and node_n+1

        use utils_h
        implicit none
		INTEGER, INTENT(IN) :: subbas_id, year, julian_day
        INTEGER, INTENT(IN) :: seasonality_array(:,:) !seasonality values as read from file
		REAL, INTENT(IN) :: support_values(:,:) !real values for n classes and 4 DOYs to be interpolated

        !real, target :: calc_seasonality2(size(support_values,dim=1))    !return value: a single value for each class (e.g. vegetation)
        !Error	in readhymo:	 error #6678: When the target is an expression it must deliver a pointer result.	E:\till\uni\wasa\wasa_svn_comp\Hillslope\readhymo.f90	1695
        !real, pointer :: calc_seasonality2(size(support_values,dim=1))    !return value: a single value for each class (e.g. vegetation)
        !error: ALLOCATABLE or POINTER attribute dictates a deferred-shape-array   [CALC_SEASONALITY2]
        real :: calc_seasonality2(size(support_values,dim=1))    !failure for some, return value: a single value for each class (e.g. vegetation)
        
        !REAL, allocatable :: calc_seasonality2(:) !fails during allocation
        !real, pointer :: calc_seasonality2(:)    !return value: a single value for each class (e.g. vegetation)

        integer    :: k, irow, search_year, istate
        integer :: d        !distance between start node and current day (in days)
        integer :: d_nodes        !distance between start node and end_node (in days)
        real :: node1_value, node2_value        !parameter values at nodepoints (start and end-point of interpolation)
        integer :: i_node1, i_node2    !indices to relevant nodes in seasonality_array
        integer :: doy_node1, doy_node2    !corresponding DOYs
        integer :: i_matchrow1, i_matchrow2    !index to matching in seasonality_array
        integer :: dy !number of days in the current year (365 or 366 in leap years)


        !allocate(calc_seasonality2(size(support_values,dim=1)), STAT = istate)
        !if (istate/=0) then
        !    write(*,'(A,i0,a)')'ERROR: Memory allocation error (',istate,') in hymo-module, seasonality computation. Try reducing nmber of SVCs or vegetation classes or disable seasonality.'
        !    stop
        !end if

        !handling of leap years
        IF (MOD((year),4) == 0)  THEN
            dy=366
        ELSE
            dy=365
        ENDIF


        if (size(seasonality_array)==1) then !no seasonality for this parameter
            calc_seasonality2=support_values(:,1)    !use single value
            return
        end if
		calc_seasonality2(:) = tiny(calc_seasonality2(1)) !flag for "not set" - initial value for all entities

		DO irow=1, size(support_values,dim=1) !do loop for all rows (i.e. all vegetation classes)
			!find matching row in seasonality array for current entity

			i_node2 = 0 !not set
			i_matchrow1 = &
			which1(	 (seasonality_array(:,1) == subbas_id .OR. seasonality_array(:,1) == -1) .AND. &
					 (seasonality_array(:,2) == irow      .OR. seasonality_array(:,2) == -1) .AND. &
					 (seasonality_array(:,3) == year      .OR. seasonality_array(:,3) == -1), nowarning=.TRUE.)

			if (i_matchrow1 == 0) cycle  !no matching row found
			i_matchrow2 = i_matchrow1    !default: other node is also in the same year (row)
   			do k=4,7
				if (seasonality_array(i_matchrow1, k) > julian_day) exit
			end do

            i_node1 = k - 4

			if (i_node1 == 0) then !current DOY is BEFORE first node, lookup other node in previous year
				search_year = -1 !search in the previous year for the other node
			elseif (i_node1 == 4) then !current DOY is AFTER last node, lookup other node in next year
				search_year =  1 !search in the next year for the other node
			else
				search_year = 0 !other node is still in the same year
			end if

			if (search_year /=0) then !search in the other year for the other node (interpolation over year break)
				i_matchrow2 = &
					which1(	 (seasonality_array(:,1) == subbas_id .OR. seasonality_array(:,1) == -1) .AND. &
							 (seasonality_array(:,2) == irow      .OR. seasonality_array(:,2) == -1) .AND. &
							 (seasonality_array(:,3) == year+search_year    .OR. seasonality_array(:,3) == -1), nowarning=.TRUE.)
			end if
			if (i_node2 == 0) i_node2 = MOD(i_node1,4) + 1   !only modify i_node2, if it has not been set before
			if (i_matchrow2 == 0) then  !no matching row found...
					if (search_year == 1) then
						i_node1 = 4  !extrapolate last value
					else
					    i_node1 = 1	!extrapolate first value
					end if
					calc_seasonality2(irow) = support_values(irow, i_node1)
					cycle
			end if
			if (i_node1 == 0) then !swap node_1 and node_2 (interpolation over year break)
				k           = i_matchrow2
				i_matchrow2 = i_matchrow1
				i_matchrow1 = k
				i_node1 = 4
			end if
			doy_node1 = seasonality_array(i_matchrow1, 3 + i_node1)
			doy_node2 = seasonality_array(i_matchrow2, 3 + i_node2)
			if (i_node1 == 4) then
				if (doy_node1 > 0) doy_node1 = doy_node1 - dy !force a negative value, as we are looking at the previous year
				if (doy_node2 < 0) doy_node2 = doy_node2 + dy !force a positive value, as we are looking at the next year
            end if
			d       = julian_day - doy_node1       !distance between start node and current day (in days)
			if (d >= dy) d = d - dy !if current day is after all nodes
            d_nodes = doy_node2  - doy_node1       !distance between start node and end_node (in days)

			if (d < 0 .OR. d_nodes < 0) then !error (presumably in input file)
				calc_seasonality2(irow) = tiny(calc_seasonality2(1)) !flag for "not set" - initial value for all entities
				cycle
			end if

			node1_value = support_values(irow, i_node1)    !parameter values at nodepoints (start and end-point of interpolation)
			node2_value = support_values(irow, i_node2)
			calc_seasonality2(irow) = node1_value+ (node2_value-node1_value) * real(d)/d_nodes        !linear interpolation between nodes
        END DO  !end loop for all rows (i.e. all vegetation classes

        return

    END FUNCTION calc_seasonality2

    FUNCTION soildistr()

    use common_h
    use params_h

    !** subroutine creates array with
    !** distribution functions of soil parameters
    !**
    !** - one distribution for each soil component (SVC)
    !** - linear interpolation between points of distribution

    IMPLICIT NONE


    !REAL, INTENT(IN)                         :: thsprof(ntcinst,maxsoil)
    !REAL, INTENT(soildistr)                        :: soildistr(5,2,ntcinst,maxsoil)
    !REAL, INTENT(soildistr)                        :: soildistr(:,:,:,:)
    REAL, pointer :: soildistr(:,:,:,:)

    !  thsprof:   input soil parameters for SVC (Vol%)
    !  soildistr:     values of distribution function of soil parameter
    !           fraction of SVC versus storage volume (mm)



    REAL :: tempx,var1,var2

    INTEGER :: k,i

    allocate(soildistr(5,2,size(thsprof, dim=1), size(thsprof, dim=2)),STAT = i)
    if (i/=0) then
        write(*,'(A,i0,a)')'ERROR: Memory allocation error (',i,') in soil-distr-module. Try disabling some hourly output.'
        stop
    end if

    !  variability around given value of SVC (first interval)
    !  (reference to soil characteristic in mm)
    var1=0.05
    !  variability around given value of SVC (second interval)
    !  (reference to soil characteristic in mm)
    var2=0.1

    !Till: saturated soil water content distribution. For each soil, the maximum water storage is modified by -10,-5,0,5,10 % (5 values). The second dimension just holds these percentages (why?)

    !** Loop for each soil component
    DO k=1,ntcinst
      DO i=1,nbr_svc(k)
        tempx=thsprof(k,i)

        if (k==5 .AND. i==3) then
            var2=0.1
        end if

        soildistr(1,1,k,i)=tempx-var2*tempx
        soildistr(2,1,k,i)=tempx-var1*tempx
        soildistr(3,1,k,i)=tempx
        soildistr(4,1,k,i)=tempx+var1*tempx
        soildistr(5,1,k,i)=tempx+var2*tempx

        soildistr(1,2,k,i)=0.0
        soildistr(2,2,k,i)=0.1
        soildistr(3,2,k,i)=0.5
        soildistr(4,2,k,i)=0.9
        soildistr(5,2,k,i)=1.0

    !  end of loop for all soil components
      END DO
    END DO

    RETURN
    END FUNCTION soildistr



end module hymo_h
