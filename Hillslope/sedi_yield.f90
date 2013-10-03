SUBROUTINE sedi_yield(d, subbas_id, lu_id, tc_type_id, q_in, q_out, q_peak_in, v_ov_in, sed_in, timestep, tc_area, sed_yield)

    ! hillslope erosion module for WASA
    ! to be called by soilwat.f90, Till Francke (till@comets.de)

    ! Till: runoff during no-rain days caused math error (log(0)) - fixed
    ! 2012-05-30

    ! Till: computationally irrelevant: removed unused parameter h
    ! 2011-05-05

    ! Till: computationally relevant: fixed overestimation of transport capacity according to Everaert due to integer division
    ! minor changes to improve compiler compatibility
    ! 2011-04-29

    ! Pedro: added variables for TC_wise output
    ! 2009-06-03

    ! Pedro: fixed stream power calculation - computationally relevant for sediment estimations
    ! 2009-06-03

    ! Till: include optional specification for LU-based correction of beta/m (exponent for L-factor calculation)
    ! 2008-10-22

    ! Till: switched back from cumulative to ordinary L-factor
    ! 2008-10-15

    ! Till: cleanup of storage structure for SVCs
    ! 2008-09-11

    ! 2008-07-15
    ! Till: clarified code

    ! 2008-07-10
    ! Till: outsourced coefficients for estimation of maximum half-hour rainfall intensity
    ! minor correction for transport capacity limit=off

    ! 2008-01-23
    ! Version with the modifications proposed by Pedro

    ! 2007-11-07
    ! Till: corrected MUSLE (SWAT-version had wrong units) to version by Krysanova et al. in Summer&Walling (2002) which is equivalent to Williams, 1995

    ! 2005-10-10
    ! Till: d50 no more calculated internally but used from global array

    ! 2005-09-29
    ! Till: structure changed from FUNCTION to SUBROUTINE: variables are passed by reference anyway and SUBROUTINE allows easier handling of array return values
    ! sediment yield is returned in fractions (MUSLE-result distributed according to mean particle size fractions in TC)

    ! 2005-08-10
    ! Till: fixed minor bugs with unit conversions

    ! 2005-07-08
    ! Till: adapt to new WASA-structure and headers

    ! 2005-06-23:
    ! Till: calculate sediment yield of a TC based on MUSLE:
    ! use the mean of flow coming into the TC and flow going out of TC as average flow


    use hymo_h
    use erosion_h
    use common_h
    use climo_h
    !use utils_h    !only for output of debug2 (otherwise remove)

    IMPLICIT NONE

    INTEGER, PARAMETER :: routing_mode=2 !different modes how incoming sediment from upslope is treated
        !(1)                        !sediment from upslope TC is stored completely, no further transport
        !(2)                        !sediment from upslope TC is transferred completely

    INTEGER, INTENT(IN):: d            !julian day of current day
    INTEGER, INTENT(IN):: subbas_id    !ID of subbasin currently treated (internal numbering scheme)
    INTEGER, INTENT(IN):: lu_id        !ID of LU currently treated TC belongs to (internal numbering scheme)
    INTEGER, INTENT(IN):: tc_type_id        !ID of TC-type that is currently treated (internal numbering scheme)
    REAL, INTENT(IN):: q_in            !overland flow entering the TC from above [m**3]
    REAL, INTENT(IN):: q_out        !overland flow leaving the TC downslope [m**3]
    REAL, INTENT(IN):: q_peak_in    !peak overland flow rate  [m**3/s]
    REAL, INTENT(IN):: sed_in(1:n_sed_class)        !sediment entering the TC from uplope [tons/timestep] for each particle size class [currently not considered yet]
    REAL, INTENT(IN):: timestep        !timestep which the given flows are related to [h]
    REAL, INTENT(IN) :: tc_area        !area of TC [km**2]
    REAL, INTENT(IN) :: v_ov_in        !overland flow velocity [m/s]
    REAL, INTENT(OUT) :: sed_yield(1:n_sed_class)                !sediment yield [tons/timestep] (usually applied daily) for each particle size class

    REAL :: q                        !mean overland flow during timestep [m H2O]
    REAL :: r, r2                    !temporary real variable

    REAL :: mean_particle(1:n_sed_class)    !auxiliary variable to compute the mean particle size distribution of current TC

    INTEGER :: i,svc_id
    REAL :: d50                        !median diameter of soil particles

    REAL :: gsize, accum1, accum2, dummy4, dummy15    !temporary auxiliary variables for computing D50
    INTEGER :: dummy14                !temporary auxiliary variable for computing D50

    REAL :: trans_cap(1:n_sed_class)        !tranport capacities
    REAL :: ef_str_pwr,tcap_everaert        !variable for computing transport capacity, transport capacity according to Everaert (1991)



    IF (q_out==0.) then        !no surface flow leaving this TC...
        sed_yield=0.            !...no sediment export
        return
    END IF

    q  =q_in                    !q_in is not used (yet), this dummy assigment prevents compiler warnings
    !q = (q_in+q_out)/2/(tc_area*1e6)    !mean overland flow during timestep [m H2O]
    q = q_out/(tc_area*1e6)        !overland flow during timestep [m H2O]

    !ii: ganze Schleife einmal berechnen für jede TC am Anfang des Programms
    K_fac        =0.
    C_fac        =0.
    P_fac        =0.
    CFRG_fac     =0.
    mean_particle=0.
    r=            0.        !for summing up fractions (actually, these should sum up to 1, but we'll do so just in case)



    !i=0
    !temp_svc=>tc_contains_svc(35)%p
    !DO WHILE (associated(temp_svc))            !average erosion-factors throughout TC, weighted by areal fraction of SVC in TC
    !    i=i+1                                !SVC counter in current TC
    !    temp_svc=>temp_svc%next                !treat next SVC
    !END DO
    !write(*,*),i



    DO i=1,size(tc_contains_svc2(tc_type_id)%p)
        svc_id=tc_contains_svc2(tc_type_id)%p(i)%svc_id                !get id of current SVC to be treated
        r2=    tc_contains_svc2(tc_type_id)%p(i)%fraction                    !get fraction of SVC
        K_fac=K_fac            +(svc_k_fac_day(svc_id))*r2            !erodibility factor that was read from svc.dat
        C_fac=C_fac            +(svc_c_fac_day(svc_id))*r2            !crop factor that was read from svc.dat
        P_fac=P_fac            +(svc_p_fac_day(svc_id))*r2            !practice factor that was read from svc.dat
        CFRG_fac=CFRG_fac   +exp(-0.03*svc_coarse_fac_day(svc_id))*r2        !coarse fragment factor that was read from svc.dat
                                            !ii: convert musle_fac(*,4) to exp(musle_fac(*,4)) at beginning of program
        mean_particle=mean_particle+soil_particles(svc_soil_veg(svc_id,1),:)*r2        !sum up all fractions for particle classes

        r=r+r2                                !sum up fractions
    end do


    K_fac=K_fac/r
    C_fac=C_fac/r
    P_fac=P_fac/r
    CFRG_fac=CFRG_fac/r
    mean_particle=mean_particle/r

    L_cum=0.    !disable cumulated L
    !Pedro - LS factor computed in a cumulative way according to Haan et al., 1994 (p. 262)
    L_cum=L_cum+L_slp        !acumulated slope length for erosion calculations (L_cum set to zero in file hymo_all.f90)

    m_ls = 0.3*(slope(tc_type_id)/100.)/((slope(tc_type_id)/100.)+exp(-1.47-(61.09*(slope(tc_type_id)/100.))))+0.2    !according to Williams (1995)

    r=1.
    IF (allocated(beta_fac)   ) r=beta_fac   (lu_id)                !use LU-prespecified correction factor for beta/m (rill/interill ratio)
    IF (allocated(beta_fac_tc)) r=beta_fac_tc(tc_type_id)        !use TC-prespecified correction factor for beta/m (rill/interill ratio), override LU value
    IF (r/=1.) THEN    !apply prespecified correction factor for beta/m (rill/interill ratio)
        !m_ls = 2*m_ls/(m_ls+1.)  !TEST option
        m_ls = r*m_ls/(-m_ls+1.+r*m_ls) !(Haan 1994, p. table 8.6; Renard et al. 1997)
    END IF


    S_fac=(slope(tc_type_id)/100.)*(65.41*(slope(tc_type_id)/100.)+4.56)+0.065
    L_fac=(L_cum**(m_ls+1.)-(L_cum-L_slp)**(m_ls+1.))/(L_slp*22.1**m_ls)
    LS_fac=L_fac*S_fac



    earea=tc_area*1e6    !area of erosive unit [m2]
    q_peak=q_peak_in
    v_ov=v_ov_in
    q_ov=(q_out)/(dt*3600./kfkorr_day)/(earea/L_slp)        !compute average overland flow rate [m**3/s] on a  1-m-strip


    !compute maximum half-hour rain: fraction of daily precipitation that falls within 30 min (needed for estimation of peak flow and rainfall erosivity)
    !IF (dt<=1) THEN        !if high resolution precipitation is available
    !    alpha_05=(maxval(preciph((d-1)*24+1:d*24,subbas_id))/precip(d,subbas_id))/(dt*2)    !estimate the maximum rainfall rate based on given precipitation data
    !ELSE
    !    IF (kfkorr_a==0) THEN
    !        alpha_05=0.5*kfkorr/24    !time invariant kfkorr
    !    ELSE
    !        alpha_05=min(1.0,0.5*(kfkorr_day/kfkorr)/24)    !time variant kfkorr: only the the precipitation-dependent part of the term is used (the kfkorr factor is considered a calibration factor for Ksat and not responsible for sub-daily rainfall dynamics)
    !    END IF
    !END IF


    IF ( (erosion_equation==1) .OR. (erosion_equation==2))  THEN
        !compute USLE-rainfall energy factor for USLE (1) and Onstad-Foster (2)
        IF (dt<=1) THEN    !if high resolution precipitation is available
            R_d=-1.        !to do: add equation here for subdaily rainfall
            r_p=-1.
        ELSE
            R_d=precip(d,subbas_id)    !daily rainfall [mm]
            !r_p=-2*R_d*log(1-min(alpha_05,0.99))    !peak rainfall rate [mm/h] Williams, 1995; eq. 25.132
            !R_05=alpha_05*R_d                        !maximum amount of rain in 30 min [mm]; Williams, 1995; 25.131
            !ri_05=R_05/0.5                            !maximum 0.5-h rainfall intensity [mm/h]
            !rainfall intensities based on kfkorr showed low values, resulting in low erosion as well
        
            if (R_d==0.) then
                ei=0.    !do computations only when there is rainfall
            else
                ri_05=a_i30*(R_d**b_i30)
                ri_05=min(ri_05,2.*R_d)                    !maximum possible intensity
                r_p=-2.*R_d*log(1-min((ri_05/2./R_d),0.99))
            end if

        END IF

        if (R_d /=0.) ei=R_d*(12.1+8.9*(log10(r_p)-0.434))*ri_05/1000.    !USLE-energy factor in the "proper" units according to Williams, 1995 in Singh,1995, p.934,25.128
    ELSE
        ei=-1.        !just for debugging
    END IF


    IF ((erosion_equation==2) .OR. (erosion_equation==3) .OR. (erosion_equation==4)) THEN
        !compute surface runoff and peak flow for Onstad-Foster (2), MUSLE (3), MUST (4)
        !given as parameter q_surf=q_out/earea/1e6/1000        !surface runoff [mm H2O]
        !given as parameter q_peak = alpha_t_conc*q_surf*earea/1e6 / (3.6*t_conc)        !(6.3.20) estimation of peak runoff rate [m**3/s]
        q_peak_mm_h=q_peak*3.6/(earea/1e6)    !convert peak runoff from [m³/s] to [mm/h]
    ELSE
        q_peak=-1.        !just for debugging
        q_peak_mm_h=-1.
    END IF

    SELECT CASE (erosion_equation)
        CASE (1)    !USLE
            r = earea*1e-4 * ei*(K_fac*C_fac*P_fac*LS_fac*CFRG_fac)
            !yield [t] according to USLE as given by Williams,1995 in Singh, 1995, p. 933
        CASE (2)    !Onstad-Foster
            r = earea*1e-4 * (0.646*ei+0.45*((q*1000.)*q_peak_mm_h)**0.33)*(K_fac*C_fac*P_fac*LS_fac*CFRG_fac)
            !yield [t] according to Onstad-Forster as given by Williams,1995 in Singh, 1995, p. 933
        CASE (3)    !MUSLE
            !ii area kürzt sich raus mit folgender Gleichung
            r = 11.8*(q*q_peak*earea)**0.56*(K_fac*C_fac*P_fac*LS_fac*CFRG_fac)
            !yield [t] according to MUSLE as given by Summer&Walling (equivalent to Williams,1995) !different units as in SWAT-manual!
        CASE (4)    !MUST
            r = earea*1e-4 * 2.5*((q*1000.)*q_peak_mm_h)**0.5*(K_fac*C_fac*P_fac*LS_fac*CFRG_fac)
            !yield [t] according to MUST as given by Williams,1995 in Singh, 1995, p. 933
        CASE DEFAULT
            WRITE (*,*) "Erosion equation ",erosion_equation," not defined, please check manual"
    END SELECT

    if (allocated(sdr_tc))    r = sdr_tc(tc_type_id) * r !use pre-specified TC-wise SDR

    sed_yield=r*mean_particle        !overall yield is distributed among size classes according to their percentage, no selective removal of finer fractions


    SELECT CASE (routing_mode)        !different modes how incoming sediment from upslope is treated
        CASE (1)        !sediment from upslope TC is stored completely

        CASE (2)        !sediment from upslope TC is transferred completely
            sed_yield = sed_yield + sed_in
    END SELECT


    SELECT CASE (transport_limit_mode)        !different modes how/if transport capacity of runoff is limited
        CASE (1)        !no limit

        CASE (2)        !transport capacity according to Everaert (1991)

            !d50 calculations for each terrain component
            accum2=0.
            DO i=1,n_sed_class
                accum2=accum2+mean_particle(i)
            ENDDO

            IF (accum2==0.) THEN
                d50=0.
            ELSE
                gsize=.5
                accum1=0.
                accum2=0.
                dummy14=0
                dummy15=0.
                DO i=1,n_sed_class
                    IF (accum2 <= gsize) THEN
                        dummy14=dummy14+1
                        accum1=accum2
                        accum2=accum2+mean_particle(i)
                    ELSE
                        EXIT
                    END IF
                ENDDO
                IF (dummy14 > 1) THEN
                    dummy15=10.**(log10(upper_limit(dummy14-1))+((log10(upper_limit(dummy14))-&
                        log10(upper_limit(dummy14-1)))*(gsize-accum1)/(accum2-accum1)))
                ELSE
                    dummy15=10.**(log10(upper_limit(dummy14)/2.)+((log10(upper_limit(dummy14))-&
                        log10(upper_limit(dummy14)/2.))*(gsize-accum1)/(accum2-accum1)))
                END IF
                d50=dummy15    !d50 em mm
            END IF

            !transport capacity
            ef_str_pwr=((1000.*1e3*9.81*q_ov*slope(tc_type_id)/100.)**1.5)/((q*100.)**(2./3.))

            IF (d50>0.15) THEN
                tcap_everaert=(4.6 *1.e-7*ef_str_pwr**1.75*(d50*1000.)**(-0.56))*(earea/L_slp*100)/1.e6*86400. !(ton)
            ELSE
                tcap_everaert=(1.74*1.e-6*ef_str_pwr**1.07*(d50*1000.)**  0.47 )*(earea/L_slp*100)/1.e6*86400. !(ton)
            END IF

            trans_cap(:)=tcap_everaert*mean_particle(:)

        CASE (3)        !Transport capacity from max MUSLE-erosion
            !trans_cap=sed_yield*0.5/K_fac        !sediment generation on current TC with maximum erodibility (K=0.5, Williams, 1995)
            trans_cap=r*0.5/K_fac        !sediment generation on current TC with maximum erodibility (K=0.5, Williams, 1995)

    END SELECT


    dummy4=sum(sed_yield)        !for comparing to yield after the application of transport capacity limits


    !deposition_TC(subbas_id,tc_type_id)=sum(sed_yield)        !Till: temporary: TC-wise output of sediment yield

    !limitation of sediment yield by transport capacity

    IF (transport_limit_mode/=1) THEN        !current sediment load is limited by transport capacity
        sed_yield=min(sed_yield,trans_cap)
    END IF


    !Pedro: TC-wise outputs
    runoff_TC(subbas_id,tc_type_id)=q*1000.
    sed_yield_TC(subbas_id,tc_type_id)=sum(sed_yield)/tc_area
    deposition_TC(subbas_id,tc_type_id)=1.-(sum(sed_yield)/(r+sum(sed_in)))

    area_TC(subbas_id,tc_type_id)=tc_area
    cum_erosion_TC(subbas_id,tc_type_id)=cum_erosion_TC(subbas_id,tc_type_id)+r
    cum_deposition_TC(subbas_id,tc_type_id)=cum_deposition_TC(subbas_id,tc_type_id)+(r+sum(sed_in)-sum(sed_yield))


    !!deposition_TC(subbas_id,tc_type_id)=q_out/tc_area
    !IF (dummy4/=0.) THEN                    !Till: for TC-wise output of deposition
    !    deposition_TC(subbas_id,tc_type_id)=1-sum(sed_yield)/dummy4    !relative deposition
    !!    deposition_TC(subbas_id,tc_type_id)=(dummy4-sum(sed_yield))/tc_area    !absolute deposition per area
    !ELSE
    !    deposition_TC(subbas_id,tc_type_id)=0.
    !ENDIF
    !

    return
end SUBROUTINE sedi_yield
