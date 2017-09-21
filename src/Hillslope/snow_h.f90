module snow_h

    use common_h
    use hymo_h
    use climo_h

    save

    !SNOW MODULE

    real :: precipSeconds                          !Length of referenceInterval (seconds)
    real :: a0                                     !Empirical coeff. (m/s)
    real :: a1                                     !Empirical coeff. (-)
    real :: kSatSnow                               !Saturated hydraulic conductivity of snow (m/s)
    real :: densDrySnow                            !Density of dry snow (kg/m3)
    real :: specCapRet                             !Capill. retent. vol as fraction of solid SWE (-)
    real :: emissivitySnowMin                      !Minimum snow emissivity used for old snow (-)
    real :: emissivitySnowMax                      !Maximum snow emissivity used for new snow (-)
    real :: tempAir_crit                           !Threshold temp. for rain-/snowfall (°C)
    real :: albedoMin                              !Minimum albedo used for old snow (-)
    real :: albedoMax                              !Maximum albedo used for new snow (-)
    real :: agingRate_tAirPos                      !Aging rate for air temperatures > 0 (1/s)
    real :: agingRate_tAirNeg                      !Aging rate for air temperatures < 0 (1/s)
    real :: soilDepth                              !Depth of interacting soil layer (m)
    real :: soilDens                               !Density of soil (kg/m3)
    real :: soilSpecHeat                           !Spec. heat capacity of soil (kJ/kg/K)
    real :: weightAirTemp                          !Weighting param. for air temp. (-) in 0...1

    real, pointer :: snowEnergyCont(:,:,:)         !Snow energy content [kJ/m²]
    real, pointer :: snowWaterEquiv(:,:,:)         !Snow water equivalent [m]
    real, pointer :: snowAlbedo(:,:,:)             !Albedo [-]

    real, pointer :: snowCover(:,:,:)              !Snow cover [-]
    real, pointer :: precipMod(:,:,:)              !Precipitation signal modified by snow module [mm]
    real, pointer :: cloudFrac(:,:,:)              !Cloud fraction [-]
    real, pointer :: precipBal(:,:,:)              !Precipitation input signal TC-wise for balance correction
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
    SUBROUTINE snow_prepare_input(hh, day, sb_counter, lu_counter2, tc_counter2, prec_mod, temp_mod, rad_mod, cloudFraction, precipBalance)
        
        implicit none

        integer, intent(IN)                  :: hh, day, sb_counter, lu_counter2, tc_counter2
        real,    intent(INOUT)               :: prec_mod, temp_mod, rad_mod
        real,    intent(INOUT)               :: cloudFraction !cloudiness fraction
        real,    intent(INOUT)               :: precipBalance !Precipiation input signal for balance check

        real                                 :: lapse_prec = 0.
        real                                 :: lapse_temp = -0.8/100
        real                                 :: lapse_rad  = 0.

        !Calculations based on approach included by developer Andreas Güntner (see etp_max.f90)
        !Calculations according to Shuttleworth (1992) Handbook of Hydrology, Chapter 4
        !IMPORTANT: nn does not equal cloudFrac => nn = 1-cloudFrac
        cloudFraction = 1 - (rad_mod/radex(day)/0.55-0.18/0.55)
        cloudFraction = MAX(0.0,cloudFraction)
        cloudFraction = MIN(1.0,cloudFraction)

        !Angstrom coefficients:
        !0.18 (fraction of extratesetrial radiation on overcast day)
        !0.55+0.18 (fraction of extraterestrial radiation on clear days)

        !simple, lapse-rate-based modification of meteo-drivers
        prec_mod = prec_mod + lapse_prec * rel_elevation(tcallid(sb_counter, lu_counter2, tc_counter2))
        temp_mod = temp_mod + lapse_temp * rel_elevation(tcallid(sb_counter, lu_counter2, tc_counter2))
        rad_mod  = rad_mod  + lapse_rad  * rel_elevation(tcallid(sb_counter, lu_counter2, tc_counter2))

        precipBalance = prec_mod

        !Here modification radiation with time, aspect,...
        !...


    END SUBROUTINE snow_prepare_input

end module snow_h
