module snow_h

    use common_h
    use snow_params

    save

    !SNOW MODULE

    real, allocatable :: snowEnergyCont(:,:,:)         !Snow energy content (kJ/m2)
    real, allocatable :: snowWaterEquiv(:,:,:)         !Snow water equivalent (m)
    real, allocatable :: snowAlbedo(:,:,:)             !Albedo (-)

    real, allocatable :: rel_elevation(:)              !Relative elevation of TC above foot of toposequence/LU (i.e. river) [m]

    real, allocatable :: snowTemp(:,:,:)               !Mean temperatur of the snow pack [°C]
    real, allocatable :: surfTemp(:,:,:)               !Snow surface temperature [°C]
    real, allocatable :: liquFrac(:,:,:)               !Fraction of liquid water (mass water / (mass water + mass ice)); Unit: Dimensionless, range 0...1
    real, allocatable :: fluxPrec(:,:,:)               !Precipitation mass flux [m/s]
    real, allocatable :: fluxSubl(:,:,:)               !Sublimation mass flux [m/s]
    real, allocatable :: fluxFlow(:,:,:)               !Meltwater flux [m/s]
    real, allocatable :: fluxNetS(:,:,:)               !Short-wave radiation balance [W/m²]
    real, allocatable :: fluxNetL(:,:,:)               !Long-wave radiation balance [W/m²]
    real, allocatable :: fluxSoil(:,:,:)               !Soil heat flux [W/m²]
    real, allocatable :: fluxSens(:,:,:)               !Sensible heat flux [W/m²]
    real, allocatable :: stoiPrec(:,:,:)               !Conversion of precipitaion mass flux (m/s) to energy flux (kJ/m2/s); Unit of result: kJ/m3
    real, allocatable :: stoiSubl(:,:,:)               !Conversion of sublimation mass flux (m/s) to energy flux (kJ/kg/K); Unit of result: kJ/m3
    real, allocatable :: stoiFlow(:,:,:)               !Conversion of meltwater loss mass flux (m/s) to energy flux (kJ/m2/s); Unit of result: kJ/m3
    real, allocatable :: rateAlbe(:,:,:)               !Change rate of albedo [1/s]

    !Variables to collect output of snow_compute
    real              :: snowEnergyCont_new
    real              :: snowWaterEquiv_new
    real              :: albedo_new
    real              :: precip_new

    real              :: TEMP_MEAN
    real              :: TEMP_SURF
    real              :: LIQU_FRAC
    real              :: flux_M_prec
    real              :: flux_M_subl
    real              :: flux_M_flow
    real              :: flux_R_netS
    real              :: flux_R_netL
    real              :: flux_R_soil
    real              :: flux_R_sens
    real              :: stoi_f_prec
    real              :: stoi_f_subl
    real              :: stoi_f_flow
    real              :: rate_G_alb

   !Computations

    CALL snow_compute(precip(i), temp(i), radia(i), airpress(i), relhumi(i), windspeed(i), cloudcover(i), &
                      snowEnergyCont(i), snowWaterEquiv(i), albedo(i), snowEnergyCont_new, snowWaterEquiv_new, albedo_new, &
                      precip_new, TEMP_MEAN, TEMP_SURF, LIQU_FRAC, flux_M_prec, flux_M_subl, flux_M_flow, flux_R_netS, &
                      flux_R_netL, flux_R_soil, flux_R_sens, stoi_f_prec, stoi_f_subl, stoi_f_flow, rate_G_alb)

       snowEnergyCont(i+1)      =      snowEnergyCont_new
       snowWaterEquiv(i+1)      =      snowWaterEquiv_new
       albedo(i+1)              =      albedo_new

       precip(i)                =      precip_new !Precipitation modified by snow accumulation and snow melt

       snowTemp(i) = TEMP_MEAN
       surfTemp(i) = TEMP_SURF
       liquFrac(i) = LIQU_FRAC
       fluxPrec(i) = flux_M_prec
       fluxSubl(i) = flux_M_subl
       fluxFlow(i) = flux_M_flow
       fluxNetS(i) = flux_R_netS
       fluxNetL(i) = flux_R_netL
       fluxSoil(i) = flux_R_soil
       fluxSens(i) = flux_R_sens
       stoiPrec(i) = stoi_f_prec
       stoiSubl(i) = stoi_f_subl
       stoiFlow(i) = stoi_f_flow
       rateAlbe(i) = rate_G_alb

    END DO


contains

    !Modification of meteo-drivers according to time and location
    SUBROUTINE snow_prepare_input(hh, day, i_subbas2, lu_counter2, tc_counter2,
        prec_mod, temp_mod, rad_mod)
        IMPLICIT NONE
        INTEGER, INTENT(IN)                  :: hh, day, i_subbas2, lu_counter2, tc_counter2
        REAL, INTENT(INOUT)             :: prec_mod, temp_mod, rad_mod
        
        !simple, lapse-rate-based modification of meteo-drivers
        prec_mod = prec_mod + lapse_prec * rel_elevation(tcallid(sb_counter, lu_counter, tc_counter))
        temp_mod = temp_mod + lapse_temp * rel_elevation(tcallid(sb_counter, lu_counter, tc_counter))
        rad_mod  = rad_mod  + lapse_rad  * rel_elevation(tcallid(sb_counter, lu_counter, tc_counter))
        
    END SUBROUTINE snow_prepare_input

end module snow_h
