
subroutine snow_compute(precipSumMM, tempAir, shortRad, pressAir, relHumid, windSpeed, cloudCoverage, &
                        snowEnergyCont, snowWaterEquiv, albedo, snowEnergyCont_new, snowWaterEquiv_new, albedo_new, &
                        TEMP_MEAN, TEMP_SURF, LIQU_FRAC, flux_M_prec, flux_M_subl, flux_M_flow, flux_R_netS, &
                        flux_R_netL, flux_R_soil, flux_R_sens, stoi_f_prec, stoi_f_subl, stoi_f_flow, rate_G_alb)

    use snow_params

    implicit none

    REAL, INTENT(IN OUT)  ::      precipSumMM             !Precipitation sum (mm/referenceInterval)
    REAL, INTENT(IN)      ::      tempAir                 !Air temperature (°C)
    REAL, INTENT(IN)      ::      shortRad                !Incoming short wave radiation (W/m2)
    REAL, INTENT(IN)      ::      pressAir                !Air pressure (hPa)
    REAL, INTENT(IN)      ::      relHumid                !Relative humidity (%)
    REAL, INTENT(IN)      ::      windSpeed               !Wind speed (m/s)
    REAL, INTENT(IN)      ::      cloudCoverage           !Cloud cover (0 = clear sky, 1 = fully covered)

    REAL, INTENT(IN)      ::      snowEnergyCont          !Snow energy content (kJ/m2)
    REAL, INTENT(IN)      ::      snowWaterEquiv          !Snow water equivalent (m)
    REAL, INTENT(IN)      ::      albedo                  !Albedo (-)

    REAL, INTENT(OUT)     ::      snowEnergyCont_new      !Snow energy content (kJ/m2)
    REAL, INTENT(OUT)     ::      snowWaterEquiv_new      !Snow water equivalent (m)
    REAL, INTENT(OUT)     ::      albedo_new              !Albedo (-)

    !REAL, INTENT(OUT)     ::      precip_new              !Precipitation modified by snow accumulation and snow melt [mm/referenceInterval]

    REAL, DIMENSION(1:5)  ::      ddt_states              !array to collect results of calculations

    REAL, DIMENSION(1:14) ::      debug                   !array to collect output for debugging purposes and in-depth analyses of the behaviour

    REAL, INTENT(OUT)     ::      TEMP_MEAN               !Mean temperatur of the snow pack [°C]
    REAL, INTENT(OUT)     ::      TEMP_SURF               !Snow surface temperature [°C]
    REAL, INTENT(OUT)     ::      LIQU_FRAC               !Fraction of liquid water (mass water / (mass water + mass ice)); Unit: Dimensionless, range 0...1
    REAL, INTENT(OUT)     ::      flux_M_prec             !Precipitation mass flux [m/s]
    REAL, INTENT(OUT)     ::      flux_M_subl             !Sublimation mass flux [m/s]
    REAL, INTENT(OUT)     ::      flux_M_flow             !Meltwater flux [m/s]
    REAL, INTENT(OUT)     ::      flux_R_netS             !Short-wave radiation balance [W/m²]
    REAL, INTENT(OUT)     ::      flux_R_netL             !Long-wave radiation balance [W/m²]
    REAL, INTENT(OUT)     ::      flux_R_soil             !Soil heat flux [W/m²]
    REAL, INTENT(OUT)     ::      flux_R_sens             !Sensible heat flux [W/m²]
    REAL, INTENT(OUT)     ::      stoi_f_prec             !Conversion of precipitaion mass flux (m/s) to energy flux (kJ/m2/s); Unit of result: kJ/m3
    REAL, INTENT(OUT)     ::      stoi_f_subl             !Conversion of sublimation mass flux (m/s) to energy flux (kJ/kg/K); Unit of result: kJ/m3
    REAL, INTENT(OUT)     ::      stoi_f_flow             !Conversion of meltwater loss mass flux (m/s) to energy flux (kJ/m2/s); Unit of result: kJ/m3
    REAL, INTENT(OUT)     ::      rate_G_alb              !Change rate of albedo [1/s]

    !Integrate and update states
    ddt_states(1:5) = snowModel_derivs(precipSumMM, shortRad, tempAir, pressAir, relHumid, windSpeed, cloudCoverage, &
                                       precipSeconds, a0, a1, kSatSnow, densDrySnow, SpecCapRet, emissivitySnowMin, &
                                       emissivitySnowMax, tempAir_crit, albedoMin, albedoMax, agingRate_tAirPos, &
                                       agingRate_tAirNeg, soilDepth, soilDens, soilSpecHeat, weightAirTemp, &
                                       snowEnergyCont, snowWaterEquiv, albedo)


    !Debug function for in-depth analyses of model behaviour
    debug(1:14) = snowModel_debug(precipSumMM, shortRad, tempAir, pressAir, relHumid, windSpeed, cloudCoverage, &
                                  precipSeconds, a0, a1, kSatSnow, densDrySnow, SpecCapRet, emissivitySnowMin, &
                                  emissivitySnowMax, tempAir_crit, albedoMin, albedoMax, agingRate_tAirPos, &
                                  agingRate_tAirNeg, soilDepth, soilDens, soilSpecHeat, weightAirTemp, &
                                  snowEnergyCont, snowWaterEquiv, albedo)


       snowEnergyCont_new     =     snowEnergyCont   +   ddt_states(1) * precipSeconds
       snowWaterEquiv_new     =     snowWaterEquiv   +   ddt_states(2) * precipSeconds
       albedo_new             =     albedo           +   ddt_states(3) * precipSeconds



    !Correct if SWE would become < 0
    if(ddt_states(2) < 0 .and. ABS(ddt_states(2))*precipSeconds > snowWaterEquiv)   then
      snowWaterEquiv_new = 0.
       snowEnergyCont_new = 0.
       albedo_new         = albedoMax
    end if


    !Special case rain on snow
    if(snowWaterEquiv_new > 0. .or. tempAir < tempAir_crit) then !if snow cover/fall

       if(tempAir > tempAir_crit) then !case rain on snow; to include energy input due to rainfall re-run with sec-new

       ddt_states(1:5) = snowModel_derivs(precipSumMM, shortRad, tempAir, pressAir, relHumid, windSpeed, cloudCoverage, &
                                          precipSeconds, a0, a1, kSatSnow, densDrySnow, SpecCapRet, emissivitySnowMin, &
                                          emissivitySnowMax, tempAir_crit, albedoMin, albedoMax, agingRate_tAirPos, &
                                          agingRate_tAirNeg, soilDepth, soilDens, soilSpecHeat, weightAirTemp, &
                                          snowEnergyCont_new, snowWaterEquiv, albedo)

       snowEnergyCont_new     =     snowEnergyCont   +   ddt_states(1) * precipSeconds
       snowWaterEquiv_new     =     snowWaterEquiv   +   ddt_states(2) * precipSeconds
       albedo_new             =     albedo           +   ddt_states(3) * precipSeconds

          !Correct if SWE would become < 0
          if(ddt_states(2) < 0 .and. ABS(ddt_states(2))*precipSeconds > snowWaterEquiv)   then
             snowWaterEquiv_new = 0.
             snowEnergyCont_new = 0.
             albedo_new         = albedoMax
          end if

       end if


       precipSumMM = min((ddt_states(4) * 1000 * precipSeconds), (snowWaterEquiv*1000 + precipSumMM)) ! mm/referenceInterval
       !in case rain on snow rain is swe melted + precipitation current time step

    end if


      TEMP_MEAN       =     debug(1)
      TEMP_SURF       =     debug(2)
      LIQU_FRAC       =     debug(3)
      flux_M_prec     =     debug(4)
      flux_M_subl     =     debug(5)
      flux_M_flow     =     debug(6)
      flux_R_netS     =     debug(7)
      flux_R_netL     =     debug(8)
      flux_R_soil     =     debug(9)
      flux_R_sens     =     debug(10)
      stoi_f_prec     =     debug(11)
      stoi_f_subl     =     debug(12)
      stoi_f_flow     =     debug(13)
      rate_G_alb      =     debug(14)

    contains

    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    ! Definition of stoichiomietry factors
    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


    !Conversion of precipitaion mass flux (m/s) to energy flux (kJ/m2/s)
    !Unit of result: kJ/m3
    function f_prec(tempAir,tempAir_crit)
        real :: tempAir
        real :: tempAir_crit
        real :: f_prec

        if (tempAir > tempAir_crit) then
        ! 4180. = 1000. * 4.18 = densityWater [kg/m3] *specHeatWater [kJ/kg/K]
        !333.5E+03 = 1000. * 333.5 = densityWater [kg/m3] * meltHeatIce [kJ/kg]
        f_prec = 4180. * max(tempAir,0.) + 333.5E+03
        else
        !2090. = 1000. * 2.09 = densityWater [kg/m3] * specHeatIce [kJ/kg/K]
        f_prec = 2090. * min(tempAir,0.)
        endif

    end function f_prec


    !Conversion of sublimation mass flux (m/s) to energy flux (kJ/kg/K)
    !Unit of result: kJ/m3
    !2837.E+03 = 1000. * 2837. = densityWater [kg/m3] * subHeatIce [kJ/kg]
    function f_subl()
        real :: f_subl

        f_subl = 2837.E+03

    end function f_subl


    !Conversion of meltwater loss mass flux (m/s) to energy flux (kJ/m2/s)
    !Unit of result: kJ/m3
    !333.5E+03 = 1000. * 333.5 = densityWater [kg/m3] * meltHeatIce [kJ/kg]
    function f_flow()
        real :: f_flow

        f_flow = 333.5E+03

    end function f_flow



    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    !General relations for estimation of derived variables
    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


    !Fraction of liquid water (mass water / (mass water + mass ice))
    !Unit of result: Dimensionless, range 0...1
    !333500. = 333.5 * 1000 = meltHeatOfIce (kJ/kg) * densWater (kg/m3)
    !The result is bounded to the range of 0...1
    !If there is no snow, a fraction of 1 is returned.
    function snowLiquidFrac(snowEnergyCont, snowWaterEquiv)
        real :: snowEnergyCont
        real :: snowWaterEquiv
        real :: snowLiquidFrac

        if (snowWaterEquiv > 0.) then
        snowLiquidFrac = min(1. , max(0., snowEnergyCont / (snowWaterEquiv * 333500.)))
        else
        snowLiquidFrac=1.
        endif

    end function snowLiquidFrac


    !Mean temperatur of the snow pack
    !Unit of result; °C (Range: -Inf...0)
    function snowTemp_mean(snowEnergyCont, snowWaterEquiv, soilDepth, soilDens, soilSpecHeat)
       real :: snowEnergyCont
       real :: snowWaterEquiv
       real :: soilDepth
       real :: soilDens
       real :: soilSpecHeat
       real :: snowTemp_mean

       if (snowWaterEquiv > 0.) then
          !If the snow pack is free of liquid water
          if (snowEnergyCont < 0.) then
          !2090. = 1000. * 2.09 = WaterDensity (kg/m3) * specHeatCapIce (kJ/kg/K)
              snowTemp_mean = snowEnergyCont / (snowWaterEquiv * 2090. + soilDepth * soilDens * soilSpecHeat)
          else
              snowTemp_mean = 0.
          !Note: Temperature for the case where all water is liquid is not computed
          !Note: If there is no snow, a temperature of zero is returned
          endif
       else
       snowTemp_mean = 0.
       endif

    end function snowTemp_mean



    !Snow surface temperature
    !Unit of result: °C (Range: -Inf...0)
    function snowTemp_surf(tempSnow_mean, tempAir, weightAirTemp)
       real :: tempSnow_mean
       real :: tempAir
       real :: weightAirTemp
       real :: snowTemp_surf

       if (tempSnow_mean < 0.) then
          snowTemp_surf = min(0. , (1.-weightAirTemp) * tempSnow_mean + weightAirTemp * tempAir)
       else
          snowTemp_surf = 0.
       endif

    end function snowTemp_Surf


    !--------------------------------------------------------------------------------------------
    !Part1: Radiation flux rates
    !--------------------------------------------------------------------------------------------

    !Short-wave radiation balance
    !Unit of result: W/m2
    function R_netS(snowWaterEquiv, shortRad, albedo)
       real :: snowWaterEquiv
       real :: shortRad
       real :: albedo
       real :: R_netS

        if (snowWaterEquiv > 0.) then
           R_netS = shortRad * (1 - albedo)
        else
           R_netS = 0.
        endif

     end function R_netS


    !Long-wave radiation balance
    !Unit of result: W/m2
    function R_netL(tempSnow_surf, emissivitySnowMin, emissivitySnowMax, tempAir, relHumid, cloudCoverage, &
                    albedo, albedoMin, albedoMax, snowWaterEquiv)
       real :: tempSnow_surf
       real :: emissivitySnowMin
       real :: emissivitySnowMax
       real :: tempAir
       real :: relHumid
       real :: cloudCoverage
       real :: albedo
       real :: albedoMin
       real :: albedoMax
       real :: snowWaterEquiv
       real :: R_outL
       real :: R_inL_clear
       real :: R_inL_cloud
       real :: R_netL

       if (snowWaterEquiv > 0.) then
          !Outgoing part
          !Note: The snow emissivity decreases with the age of the snow cover (Dyck &
          !      Peschke, 1995). We use the dynamically computed albedo as an indicator
          !      for the age of the snow pack. The constant 0.001 is to prevent a Div0
          !      in case of a constant albedo (i.e. albedoMax = albedoMin).
          !5.67E-08 = Stefan-Boltzmann constant (W/m2/K4)

          R_outL = (emissivitySnowMin + (emissivitySnowMax - emissivitySnowMin) &
                   * (albedo-albedoMin+0.001) / (albedoMax-albedoMin+0.001)) &
                   * 5.67E-08 * ((tempSnow_surf + 273.15)**4.)

          !Incoming part 1 - Clear sky emission
          !This is the Stefan-Boltzmann equation with cleaer sky emissivity estimated
          !from vaporPressure using the empirical Brunt-formula
          R_inL_clear = (0.51 + 0.066 * sqrt(vapPress_overIce(tempAir,relHumid))) * 5.67E-08 * ((tempAir + 273.15)**4.)

          !Incoming part 2 - Cloud emission
          !This is the Stefan-Boltzmann equation with the emissivity set to 1 (clouds
          !treated as black body). The cloud temperature is approximated by the dewpoint temperature

          R_inL_cloud = 1.0 * 5.67E-08 * (dewTemp_overIce(vapPress_overIce(tempAir,relHumid)) + 273.15)**4.

          R_netL = (1. - cloudCoverage) * R_inL_clear + cloudCoverage * R_inL_cloud - R_outL

       else
          R_netL = 0.
       endif

    end function R_netL


    !Soil heat flux
    !Unit of result: W/m2
    function R_soil()
       real :: R_soil

       R_soil = 0.

    end function R_soil


    !Sensible heat flux
    !Unit of result: W/m2
    function R_sens(tempSnow_surf, tempAir, pressAir, windSpeed, a0, a1, snowWaterEquiv)
       real :: tempSnow_surf
       real :: tempAir
       real :: pressAir
       real :: windSpeed
       real :: a0
       real :: a1
       real :: snowWaterEquiv
       real :: R_sens

       if (snowWaterEquiv > 0.) then
          R_sens = (a0 +a1 *windSpeed) * densityDryAir(tempAir,pressAir) * 1005. &
          * (tempAir - tempSnow_surf)
       else
          R_sens = 0.
       endif
    end function R_sens

    !--------------------------------------------------------------------------------------------
    !Part2: Mass flux rates
    !--------------------------------------------------------------------------------------------


    !Precipitation mass flux
    !Unit of result: m/s
    function M_prec(precipSumMM, precipSeconds)
       real :: precipSumMM
       real :: precipSeconds
       real :: M_prec

       M_prec = precipSumMM /1000. / precipSeconds

    end function M_prec


    !Sublimation mass flux
    !Unit of result: m/s
    function M_subl(tempSnow_surf, tempAir, pressAir, relHumid, windSpeed, a0, a1, snowWaterEquiv)
       real :: tempSnow_surf
       real :: tempAir
       real :: pressAir
       real :: relHumid
       real :: windSpeed
       real :: a0
       real :: a1
       real :: snowWaterEquiv
       real :: M_subl


       if (snowWaterEquiv > 0.) then
       !(a0 + a1 * windSpeed) = Empirical estimatin of transfer coefficient (m/s)
       !1000. = Density of water (kg/m3)
       !100. = Relative humidity at saturation (%)
       !Negative flux rates are set to zero (no freezing air moisture)
          M_subl = max(0., (a0 + a1 * windSpeed) * densityDryAir(tempAir,pressAir) / 1000. &
                   * (specificHumidity(pressAir, vapPress_overIce(tempSnow_surf,100.)) &
                   -  specificHumidity(pressAir, vapPress_overIce(tempAir,relHumid))))
       else
          M_subl = 0.
       endif
    end function M_subl


    !Meltwater flux
    !Unit of result: m/s
    function M_flow(snowLiquidFrac, kSatSnow, densDrySnow, specCapRet, snowWaterEquiv)
       real :: snowLiquidFrac
       real :: kSatSnow
       real :: densDrySnow
       real :: specCapRet
       real :: snowWaterEquiv
       real :: M_flow
       real :: rss

       if(snowWaterEquiv > 0.) then
          !Relative saturation (-)
          !Don´t allow 100% liquid water as this will cause a division by zero

          rss= (min(0.99, snowLiquidFrac) / (1.-min(0.99, snowLiquidFrac)) - specCapRet) &
          / (1000./densDrySnow - 1000./922. - specCapRet)
          !1000. = Density of water (kg/m3)
          !922 = Density of ice (kg/m3)
          !Fix negative values of the relative saturation (-)
          M_flow = kSatSnow * max(0., rss)**3.

       else
          M_flow = 0.
       endif
    end function M_flow



    !--------------------------------------------------------------------------------------------
    !Part3: Other rates
    !--------------------------------------------------------------------------------------------

    !Change rate of albedo
    !Unit of results: 1/s
    function G_alb(albedo, precipSumMM, precipSeconds, tempAir, tempAir_crit, albedoMin, albedoMax, &
                   agingRate_tAirPos, agingRate_tAirNeg, snowWaterEquiv)
       real :: albedo
       real :: precipSumMM
       real :: precipSeconds
       real :: tempAir
       real :: tempAir_crit
       real :: albedoMin
       real :: albedoMax
       real :: agingRate_tAirPos
       real :: agingRate_tAirNeg
       real :: snowWaterEquiv
       real :: G_alb

       if(snowWaterEquiv > 0.) then
       !Surface renewal if snow falls
       !Time of renewal set to the reference interval to keep the change rate reasonably
       !small in compariosn to other change rates (prevents stiff ODE system)
       !Complete renewal is assumed for a snowfall of 10mm (as in UTAH Energy Balance Model) and
       !a partly renewal is assumed for less precipitation
          if((precipSumMM > 0.) .and. (tempAir < tempAir_crit)) then
             G_alb = (albedoMax - albedo) / precipSeconds * min(1., precipSumMM/10.) !This is positive
             !Surface aging
          else
             if (tempAir >= 0.) then
                G_alb = agingRate_tAirPos * (albedoMin - albedo) !This is negative
             else
                G_alb = agingRate_tAirNeg * (albedoMin - albedo) !This is negative

             endif
          endif
       else
          G_alb=0.
       endif
    end function G_alb



    !--------------------------------------------------------------------------------------------
    !Meteorological basics
    !--------------------------------------------------------------------------------------------


    !Saturation vapor pressure as a function of temperature
    !Source: Dyck & Peschke (1995), p. 58, Eq. 4.6
    !Input: Temperature (°C)
    !Output: Result in (hPa)

    function satVapPress_overWater(temp)
        real :: temp
        real :: satVapPress_overWater
        satVapPress_overWater = 6.11 * 10.**(7.5*temp/(237.3+temp))
    end function satVapPress_overWater

    function satVapPress_overIce(temp)
        real :: temp
        real :: satVapPress_overIce
        satVapPress_overIce = 6.11 * 10.**(9.5*temp/(265.5+temp))
    end function satVapPress_overIce

    !Vapor pressure as a funciton of temperature and relative humidity
    !Source: Dycka & Peschke (1995), p.65, Eq. 4.12 together with the equation to compute
    !the saturation vapor pressure
    !Input: Temperature (degree celsius)
    !Input: Relative humidity (%)
    !Output: Result in (hPa)

    function vapPress_overWater(temp, relhum)
        real :: temp
        real :: relhum
        real :: vapPress_overWater
        vapPress_overWater = satVapPress_overWater(temp) * relhum / 100.
    end function vapPress_overWater


    function vapPress_overIce(temp, relhum)
        real :: temp
        real :: relhum
        real :: vapPress_overIce
        vapPress_overIce = satVapPress_overIce(temp) * relhum / 100.
    end function vapPress_overIce


    !Dew point temperature as a function of temperature and relative humidity
    !Source: See explanation in documentation of the snow model
    !Input: Vapor Pressure (hPa)
    !Output: Result in (hPa)

    function dewTemp_overWater(vapPress)
        real :: vapPress
        real :: dewTemp_overWater
        dewTemp_overWater = 237.3 * log10(vapPress / 6.11) / (7.5 - log10(vapPress / 6.11))
    end function dewTemp_overWater


    function dewTemp_overIce(vapPress)
        real :: vapPress
        real :: dewTemp_overIce
        dewTemp_overIce = 265.5 * log10(vapPress / 6.11) / (9.5 - log10(vapPress / 6.11))
    end function dewTemp_overIce



    !Density of dry air
    !Note: The error when using this equation for moist air is low as under normal atmospheric condition
    !the vapor pressure is small when compared ot the air pressure
    !Source: Dyck & Peschke (1995), p.60, Eq.4.9
    !Input: Temperature (°C)
    !Output: Air pressure (hPa)

    function densityDryAir(temp,airpress)
        real :: temp
        real :: airpress
        real :: densityDryAir
        densityDryAir = airpress * 0.1 / (0.287 * (273.15 + temp))
    end function densityDryAir

    !Specific Humidity
    !Source: Dycke & Peschke (1995), p.60, Eq. 4.10
    !Input: Air pressure (hPa)
    !Input: Vapor pressure (hPa)
    !Output: Result is dimensionless

    function specificHumidity(pressAir,vapPress)
        real :: pressAir
        real :: vapPress
        real :: specificHumidity
        specificHumidity = 0.622 * vapPress / (pressAir - 0.378 * vapPress)
    end function specificHumidity



    !--------------------------------------------------------------------------------------------
    !Derivates of the snow model´s state varibles with respect to time
    !--------------------------------------------------------------------------------------------

    function snowModel_derivs(precipSumMM, shortRad, tempAir, pressAir, relHumid, windSpeed, cloudCoverage, &
                              precipSeconds, a0, a1, kSatSnow, densDrySnow, SpecCapRet, emissivitySnowMin, &
                              emissivitySnowMax, tempAir_crit, albedoMin, albedoMax, agingRate_tAirPos, &
                              agingRate_tAirNeg, soilDepth, soilDens, soilSpecHeat, weightAirTemp, &
                              snowEnergyCont, snowWaterEquiv, albedo) result(res)
        !Inputs
        !Part 1: Forcings
        real :: precipSumMM
        real :: shortRad
        real :: tempAir
        real :: pressAir
        real :: relHumid
        real :: windSpeed
        real :: cloudCoverage

        !Part 2: Parameters
        real :: precipSeconds
        real :: a0
        real :: a1
        real :: kSatSnow
        real :: densDrySnow
        real :: SpecCapRet
        real :: emissivitySnowMin
        real :: emissivitySnowMax
        real :: tempAir_crit
        real :: albedoMin
        real :: albedoMax
        real :: agingRate_tAirPos
        real :: agingRate_tAirNeg
        real :: soilDepth
        real :: soilDens
        real :: soilSpecHeat
        real :: weightAirTemp

        !Part 3: States
        real :: snowEnergyCont
        real :: snowWaterEquiv
        real :: albedo

        !Derived Variables and rate expressions
        real :: TEMP_MEAN
        real :: TEMP_SURF
        real :: LIQU_FRAC
        real :: M_P
        real :: M_S
        real :: M_F

        !Outputs
        real :: ddt_sec
        real :: ddt_swe
        real :: ddt_alb
        real :: flux_melt
        real :: flux_subl
        real, dimension(:), allocatable :: res

        !Derived variables
        TEMP_MEAN = snowTemp_mean(snowEnergyCont, snowWaterEquiv, soilDepth, soilDens, soilSpecHeat)
        TEMP_SURF = snowTemp_surf(TEMP_MEAN, tempAir, weightAirTemp)
        LIQU_FRAC = snowLiquidFrac(snowEnergyCont, snowWaterEquiv)

        !Rate expressions used multiple times
        M_P = M_prec(precipSumMM, precipSeconds)

        if(snowWaterEquiv < 0.0001 .and. tempAir > tempAir_crit) then !if no snow cover present and precipitation liquid, no addition to swe
        M_P = 0.
        end if

        M_S = M_subl(TEMP_SURF, tempAir, pressAir, relHumid, windSpeed, a0, a1, snowWaterEquiv)
        M_F = M_flow(LIQU_FRAC, kSatSnow, densDrySnow, specCapRet, snowWaterEquiv)

        !Computation of derivatives
        ddt_sec = 0.001 * (R_netS(snowWaterEquiv, shortRad, albedo) + &
                           R_netL(TEMP_SURF, emissivitySnowMin, emissivitySnowMax, tempAir, relHumid, &
                           cloudCoverage, albedo, albedoMin, albedoMax, snowWaterEquiv) + &
                           R_soil() + &
                           R_sens(TEMP_SURF, tempAir, pressAir, windSpeed, a0, a1, snowWaterEquiv)) + &
                           f_prec(tempAir,tempAir_crit) * M_P - &
                           f_subl() * M_S - &
                           f_flow() * M_F
        ddt_swe = M_P - M_S - M_F
        ddt_alb = G_alb(albedo, precipSumMM, precipSeconds, tempAir, tempAir_crit, albedoMin, albedoMax, &
                        agingRate_tAirPos, agingRate_tAirNeg, snowWaterEquiv)
        flux_melt = M_F
        flux_subl = M_S

        res = (/ ddt_sec, ddt_swe, ddt_alb, flux_melt, flux_subl/)

     end function snowModel_derivs



     !--------------------------------------------------------------------------------------------
     !Function to return the fundamental variables of the snow model
     !(for debugging purposes and in-depth analyses of the behaviour)
     !--------------------------------------------------------------------------------------------

     function snowModel_debug(precipSumMM, shortRad, tempAir, pressAir, relHumid, windSpeed, cloudCoverage, &
                              precipSeconds, a0, a1, kSatSnow, densDrySnow, SpecCapRet, emissivitySnowMin, &
                              emissivitySnowMax, tempAir_crit, albedoMin, albedoMax, agingRate_tAirPos, &
                              agingRate_tAirNeg, soilDepth, soilDens, soilSpecHeat, weightAirTemp, &
                              snowEnergyCont, snowWaterEquiv, albedo) result(res)
        !Inputs
        !Part 1: Forcings
        real :: precipSumMM
        real :: shortRad
        real :: tempAir
        real :: pressAir
        real :: relHumid
        real :: windSpeed
        real :: cloudCoverage

        !Part 2: Parameters
        real :: precipSeconds
        real :: a0
        real :: a1
        real :: kSatSnow
        real :: densDrySnow
        real :: SpecCapRet
        real :: emissivitySnowMin
        real :: emissivitySnowMax
        real :: tempAir_crit
        real :: albedoMin
        real :: albedoMax
        real :: agingRate_tAirPos
        real :: agingRate_tAirNeg
        real :: soilDepth
        real :: soilDens
        real :: soilSpecHeat
        real :: weightAirTemp

        !Part 3: States
        real :: snowEnergyCont
        real :: snowWaterEquiv
        real :: albedo

        !Outputs
        real :: TEMP_MEAN
        real :: TEMP_SURF
        real :: LIQU_FRAC
        real :: flux_R_netS
        real :: flux_R_netL
        real :: flux_R_soil
        real :: flux_R_sens
        real :: stoi_f_prec
        real :: stoi_f_subl
        real :: stoi_f_flow
        real :: flux_M_prec
        real :: flux_M_subl
        real :: flux_M_flow
        real :: rate_G_alb
        real, dimension(:), allocatable :: res

        !Derived variables
        TEMP_MEAN= snowTemp_mean(snowEnergyCont, snowWaterEquiv, soilDepth, soilDens, soilSpecHeat)
        TEMP_SURF= snowTemp_surf(TEMP_MEAN, tempAir, weightAirTemp)
        LIQU_FRAC= snowLiquidFrac(snowEnergyCont, snowWaterEquiv)

        !Mass fluxes
        flux_M_prec= M_prec(precipSumMM, precipSeconds)
        flux_M_subl= M_subl(TEMP_SURF, tempAir, pressAir, relHumid, windSpeed, a0, a1, snowWaterEquiv)
        flux_M_flow= M_flow(LIQU_FRAC, kSatSnow, densDrySnow, specCapRet, snowWaterEquiv)

        !Radiation fluxes
        flux_R_netS= R_netS(shortRad, albedo, snowWaterEquiv)
        flux_R_netL= R_netL(TEMP_SURF, emissivitySnowMin, emissivitySnowMax, tempAir, relHumid, &
                            cloudCoverage, albedo, albedoMin, albedoMax, snowWaterEquiv)
        flux_R_soil= R_soil()
        flux_R_sens= R_sens(TEMP_SURF, tempAir, pressAir, windSpeed, a0, a1, snowWaterEquiv)

        !Stochiometry factors
        stoi_f_prec= f_prec(tempAir, tempAir_crit)
        stoi_f_subl= f_subl()
        stoi_f_flow= f_flow()

        !Other rates
        rate_G_alb= G_alb(albedo, precipSumMM, precipSeconds, tempAir, tempAir_crit, albedoMin, &
                          albedoMax, agingRate_tAirPos, agingRate_tAirNeg, snowWaterEquiv)

        res = (/TEMP_MEAN, TEMP_SURF, LIQU_FRAC, flux_M_prec, flux_M_subl, flux_M_flow, &
                flux_R_netS, flux_R_netL, flux_R_soil, flux_R_sens, &
                stoi_f_prec, stoi_f_subl, stoi_f_flow, rate_G_alb/)

     end function snowModel_debug

end subroutine snow_compute
