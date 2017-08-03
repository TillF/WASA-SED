module snow_h

    use common_h
    use snow_params
    use hymo_h

    save

    !SNOW MODULE

    real, pointer :: snowEnergyCont(:,:,:)         !Snow energy content (kJ/m2)
    real, pointer :: snowWaterEquiv(:,:,:)         !Snow water equivalent (m)
    real, pointer :: snowAlbedo(:,:,:)             !Albedo (-)

    real, pointer :: rel_elevation(:)              !Relative elevation of TC above foot of toposequence/LU (i.e. river) [m]

    real, pointer :: snowTemp(:,:,:)               !Mean temperatur of the snow pack [°C]
    real, pointer :: surfTemp(:,:,:)               !Snow surface temperature [°C]
    real, pointer :: liquFrac(:,:,:)               !Fraction of liquid water (mass water / (mass water + mass ice)); Unit: Dimensionless, range 0...1
    real, pointer :: fluxPrec(:,:,:)               !Precipitation mass flux [m/s]
    real, pointer :: fluxSubl(:,:,:)               !Sublimation mass flux [m/s]
    real, pointer :: fluxFlow(:,:,:)               !Meltwater flux [m/s]
    real, pointer :: fluxNetS(:,:,:)               !Short-wave radiation balance [W/m²]
    real, pointer :: fluxNetL(:,:,:)               !Long-wave radiation balance [W/m²]
    real, pointer :: fluxSoil(:,:,:)               !Soil heat flux [W/m²]
    real, pointer :: fluxSens(:,:,:)               !Sensible heat flux [W/m²]
    real, pointer :: stoiPrec(:,:,:)               !Conversion of precipitaion mass flux (m/s) to energy flux (kJ/m2/s); Unit of result: kJ/m3
    real, pointer :: stoiSubl(:,:,:)               !Conversion of sublimation mass flux (m/s) to energy flux (kJ/kg/K); Unit of result: kJ/m3
    real, pointer :: stoiFlow(:,:,:)               !Conversion of meltwater loss mass flux (m/s) to energy flux (kJ/m2/s); Unit of result: kJ/m3
    real, pointer :: rateAlbe(:,:,:)               !Change rate of albedo [1/s]



contains

    !Modification of meteo-drivers according to time and location
    SUBROUTINE snow_prepare_input(hh, day, sb_counter, lu_counter2, tc_counter2, prec_mod, temp_mod, rad_mod)
        
        implicit none

        integer, intent(IN)                  :: hh, day, sb_counter, lu_counter2, tc_counter2
        real,    intent(INOUT)               :: prec_mod, temp_mod, rad_mod

        real                                 :: lapse_prec = 0.
        real                                 :: lapse_temp = 0.
        real                                 :: lapse_rad  = 0.

        !simple, lapse-rate-based modification of meteo-drivers
        prec_mod = prec_mod + lapse_prec * rel_elevation(tcallid(sb_counter, lu_counter2, tc_counter2))
        temp_mod = temp_mod + lapse_temp * rel_elevation(tcallid(sb_counter, lu_counter2, tc_counter2))
        rad_mod  = rad_mod  + lapse_rad  * rel_elevation(tcallid(sb_counter, lu_counter2, tc_counter2))
        !Here modification radiation with time, aspect,...
    END SUBROUTINE snow_prepare_input

end module snow_h
