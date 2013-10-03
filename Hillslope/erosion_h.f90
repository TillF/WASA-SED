! Till: added specifiyable SDR for TCs
! 2013-10-02

! Pedro: added variables for TC_wise output
! 2009-06-03

! Till: include optional specification for LU-based correction of beta/m (exponent for L-factor calculation)
! 2008-10-22

! Till: added specifyable SDR for LUs
! 2008-09-09

! Till: implemented variables for optional pre-specified sediment outflow of selected subbasins
! 2008-07-11

! Till: included coefficients for estimation of maximum half-hour rainfall intensity 
! 2008-07-10

! Pedro
! added variables
! 2008-02-07

! Till: changed addressing scheme of sediment_subbasin_t from (366*nt,subasin,n_sed_class) to (366,nt,subasin,n_sed_class)
! 2005-12-20


module erosion_h



    INTEGER  :: erosion_equation            !erosion equation to be used:1: USLE, 2: Onstad-Foster, 3:    MUSLE, 4: MUST
    LOGICAL :: do_musle_subbasin        !1: compute erosion on subbasin scale; 0: compute erosion on TC-scale

    !Las Paules
    REAL    :: a_i30=1.1630            ! coefficients for estimation of maximum half-hour rainfall intensity (ri_05) from daily rainfall data (R_day) (default, taken from erosion.ctl, if present)
    REAL    :: b_i30=0.667981            ! ri_05=a*R_day^b


    INTEGER :: transport_limit_mode=2    !different modes how/if transport capacity of runoff is limited (default, taken from erosion.ctl, if present)
        !(1)                        !no transport capacity limit
        !(2)                        !transport capacity according to Everaert (1991)
        !(3)                        !transport capacity computed from MUSLE with maximum erodibility

                                

    TYPE t_svc2
        INTEGER :: svc_id
        REAL    :: fraction
    END TYPE t_svc2

    TYPE t_svc_pointer_array
        type (t_svc2), pointer :: p(:)
    END TYPE t_svc_pointer_array

       


    !Till: Hillslope erosion
    !real, allocatable :: musle_fac(:,:)    !contains MUSLE-factors (K C P CFRG ManningsN) for each SVC (LS comes from the TC itself)
    real, pointer :: svc_k_fac(:,:)    !contains MUSLE-factor K for each SVC and seasonal node(LS comes from the TC itself)
    real, pointer :: svc_c_fac(:,:)    !contains MUSLE-factor C for each SVC and seasonal node
    real, pointer :: svc_p_fac(:,:)    !contains MUSLE-factor P for each SVC and seasonal node
    real, pointer :: svc_coarse_fac(:,:)    !contains MUSLE-factor CFRG  for each SVC and seasonal node
    real, pointer :: svc_n(:,:)    !contains ManningsN for each SVC and seasonal node

    integer, allocatable :: id_svc_extern(:)    !external ids of soil vegetation components
    integer, allocatable :: svc_soil_veg(:,:)    !contains internal IDs of soils and vegetation units belonging to each SVC
    type (t_svc_pointer_array), allocatable :: tc_contains_svc2(:)            !contains which SVCs are contained in each TC (indexed by TC-type)

    ! sediment flux for each landscape unit and particle size class [t] in the current subbasin
    REAL, allocatable  :: sedsu(:,:)


    !sediment contribution of each subbasin per timestep (analogous to water_subbasin), for each particle class [t/timestep]
    real, allocatable ::  sediment_subbasin(:,:,:)

    !sediment contribution of each subbasin per timestep(366*nt,subasin,n_sed_class) (analogous to water_subbasin), for each particle class [t/timestep]
    !alt sedi real, allocatable ::  sediment_subbasin_t(:,:,:) !sediment contribution of each subbasin per timestep (analogous to water_subbasin), for each particle class [t/timestep]
    !sediment contribution of each subbasin per timestep(366,nt,subasin,n_sed_class) (analogous to water_subbasin), for each particle class [t/timestep]
    real, allocatable ::  sediment_subbasin_t(:,:,:,:) !sediment contribution of each subbasin per timestep (analogous to water_subbasin), for each particle class [t/timestep]

    real, allocatable :: soil_particles(:,:)        !fraction of particle classes for each soil (uppermost horizon)


    ! TOTAL NUMBER OF SVCs IN STUDY AREA
    INTEGER :: nsvc


    !auxiliary vars for erosivity computation (computation at subbasin scale)
    REAL :: q_peak, q_ov, v_ov, t_conc, alpha_t_conc, alpha_05,manning_n,h_depth         !variables used for calculation of peak runoff rate

    REAL :: L_fac, S_fac, LS_fac, K_fac, C_fac, P_fac, CFRG_fac, L_cum        !(M)USLE erosion factors
    REAL :: m_ls                        !auxiliary variable for computing LS-factor
    REAL :: L_slp                    !slope length of TC

    real ::    R_d            !daily rainfall [mm]
    real ::    r_p            !peak rainfall rate [mm/h]
    real ::    R_05        !maximum amount of rain in 30 min [mm]
    real ::    ri_05        !maximum 0.5-h rainfall intensity [mm/h]
    real ::    q_peak_mm_h    !peak runoff as before in [mm/h]
    real :: ei            !USLE energy factor ["proper units" (Williams, 1995]
    real :: earea        !area of erosive unit [m2]

    real, allocatable :: runoff_TC(:,:)
    real, allocatable :: sed_yield_TC(:,:)
    real, allocatable :: deposition_TC(:,:)
    real, allocatable :: area_TC(:,:)
    real, allocatable :: cum_erosion_TC(:,:)
    real, allocatable :: cum_deposition_TC(:,:)



    real, allocatable ::  pre_subbas_outsed(:,:,:)    !Till: pre-specified outflow of sediments of subbasins (optionally read from file) [m³/s]
    integer, pointer :: corr_column_pre_subbas_outsed(:)    !Till: holds corresponding columns of input files to be related to internal numbering of subbasins
    real, allocatable :: pre_psd(:)                        !Till: mean particle size distribution of pre-specified sediment outflow

    real, allocatable :: sdr_lu(:)                        !Till: optionally specified sediment delivery ratio for landscape units on LU-scale
    real, allocatable :: sdr_tc(:)                      !Till: optionally specified sediment delivery ratio for landscape units on TC-scale

    real, allocatable :: beta_fac(:)                !Till: optionally specified correction factor for beta on LU-scale(rill/interrill ratio, used for the computation of L-factor, see Renard, 1997, pp.101
    real, allocatable :: beta_fac_tc(:)                !Till: optionally specified correction factor for beta on TC-scale(rill/interrill ratio, used for the computation of L-factor, see Renard, 1997, pp.101



    integer, pointer :: seasonality_k(:,:)    ! four key nodes in time for temporal dynamics of K-factor within year (index: subbasin,(1:4)*simulation_year,)
    REAL, pointer :: svc_k_fac_day(:)        !erosion factor dynamics of current day and subbasin

    integer, pointer :: seasonality_c(:,:)    ! four key nodes in time for temporal dynamics of C-factor within year (index: subbasin,(1:4)*simulation_year,)
    REAL, pointer :: svc_c_fac_day(:)        !erosion factor dynamics of current day and subbasin

    integer, pointer :: seasonality_p(:,:)    ! four key nodes in time for temporal dynamics of P-factor within year (index: subbasin,(1:4)*simulation_year,)
    REAL, pointer :: svc_p_fac_day(:)        !erosion factor dynamics of current day and subbasin

    integer, pointer :: seasonality_coarse(:,:)    ! four key nodes in time for temporal dynamics of coarse fraction-factor within year (index: subbasin,(1:4)*simulation_year,)
    REAL, pointer :: svc_coarse_fac_day(:)        !erosion factor dynamics of current day and subbasin

    integer, pointer :: seasonality_n(:,:)    ! four key nodes in time for temporal dynamics of Mannings n within year (index: subbasin,(1:4)*simulation_year,)
    REAL, pointer :: svc_n_day(:)        !erosion factor dynamics of current day and subbasin



end module erosion_h
