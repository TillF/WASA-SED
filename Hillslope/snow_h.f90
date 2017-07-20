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


contains

    !Modification of meteo-drivers according to time and location
    SUBROUTINE snow_prepare_input(hh, day, i_subbas2, lu_counter2, tc_counter2,
        prec_mod, temp_mod, rad_mod)
        
        implicit none

        integer, intent(IN)                  :: hh, day, i_subbas2, lu_counter2, tc_counter2
        real,    intent(INOUT)               :: prec_mod, temp_mod, rad_mod

        real                                 :: lapse_prec = 0.
        real                                 :: lapse_temp = 0.
        real                                 :: lapse_rad  = 0.

        !simple, lapse-rate-based modification of meteo-drivers
        prec_mod = prec_mod + lapse_prec * rel_elevation(tcallid(sb_counter, lu_counter, tc_counter))
        temp_mod = temp_mod + lapse_temp * rel_elevation(tcallid(sb_counter, lu_counter, tc_counter))
        rad_mod  = rad_mod  + lapse_rad  * rel_elevation(tcallid(sb_counter, lu_counter, tc_counter))
        
    END SUBROUTINE snow_prepare_input

end module snow_h
