
subroutine snow_compute(precipSumMM, tempAir, shortRad, pressAir, relHumid, windSpeed, cloudCoverage, &
                        snowEnergyContent, snowWaterEquivalent, albedo, snowEnergyCont_new, snowWaterEquiv_new, albedo_new, &
                        snowCov, TEMP_MEAN, TEMP_SURF, LIQU_FRAC, flux_M_prec, flux_M_subl, flux_M_flow, flux_R_netS, &
                        flux_R_netL, flux_R_soil, flux_R_sens, stoi_f_prec, stoi_f_subl, stoi_f_flow, rate_G_alb, &
                        precipModif, cloudFraction)

    use snow_h

    implicit none

    REAL, INTENT(IN)      ::      precipSumMM             !Precipitation sum (mm/referenceInterval)

    REAL, INTENT(IN), DIMENSION(1:24)  ::  shortRad       !Incoming short wave radiation (W/m2)
    REAL, INTENT(IN), DIMENSION(1:24)  ::  tempAir        !Air temperature (°C)

    REAL, INTENT(IN)      ::      pressAir                !Air pressure (hPa)
    REAL, INTENT(IN)      ::      relHumid                !Relative humidity (%)
    REAL, INTENT(IN)      ::      windSpeed               !Wind speed (m/s)
    REAL, INTENT(IN)      ::      cloudCoverage           !Cloud cover (0 = clear sky, 1 = fully covered)

    REAL, INTENT(IN)      ::      snowEnergyContent       !Snow energy content (kJ/m2)
    REAL, INTENT(IN)      ::      snowWaterEquivalent     !Snow water equivalent (m)
    REAL, INTENT(IN)      ::      albedo                  !Albedo (-)

    REAL, INTENT(OUT)     ::      snowEnergyCont_new      !Snow energy content (kJ/m2)
    REAL, INTENT(OUT)     ::      snowWaterEquiv_new      !Snow water equivalent (m)
    REAL, INTENT(OUT)     ::      albedo_new              !Albedo (-)

    REAL, INTENT(OUT)     ::      snowCov                 !Snow cover (-)
    REAL, INTENT(OUT)     ::      precipModif             !Precipitation modified by snow module [mm] (collect for debugging purposes)
    REAL, INTENT(OUT)     ::      cloudFraction           !cloud fraction caclulated based on radiation [-] (collect for debugging purposes)

    REAL, DIMENSION(1:17) ::      snowResults             !array to collect results of calculations

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

    INTEGER               ::     i                        !counter for do loop
    INTEGER               ::     n                        !amount intermediate timesteps


    !Do full calculations with intermediate steps only if snow or rain
    if(snowWaterEquivalent > 0. .OR.  precipSumMM > 0.)then


       snowEnergyCont_new     =     snowEnergyContent
       snowWaterEquiv_new     =     snowWaterEquivalent
       albedo_new             =     albedo
       flux_M_flow            =     0.
       flux_M_subl            =     0.
       flux_M_prec            =     0.
       TEMP_MEAN              =     0.
       TEMP_SURF              =     0.
       LIQU_FRAC              =     0.
       flux_R_netS            =     0.
       flux_R_netL            =     0.
       flux_R_soil            =     0.
       flux_R_sens            =     0.
       stoi_f_prec            =     0.
       stoi_f_subl            =     0.
       stoi_f_flow            =     0.
       rate_G_alb             =     0.

       if(precipSeconds >= 84600.)then
          n = 24
       else
          n=1
       end if

       DO i=1,n !Daily resolution data too coarse; causing instabilities in melting process

          snowResults(1:17) = f_snowModel(precipSumMM/n, shortRad(i), tempAir(i), pressAir, relHumid, windSpeed, cloudCoverage, &
                                          precipSeconds/n, a0, a1, kSatSnow, densDrySnow, SpecCapRet, emissivitySnowMin, &
                                          emissivitySnowMax, tempAir_crit, albedoMin, albedoMax, agingRate_tAirPos, &
                                          agingRate_tAirNeg, soilDepth, soilDens, soilSpecHeat, weightAirTemp, &
                                          snowEnergyCont_new, snowWaterEquiv_new, albedo_new)

          snowEnergyCont_new     =     snowEnergyCont_new   +   snowResults(1) * precipSeconds/n
          snowWaterEquiv_new     =     snowWaterEquiv_new   +   snowResults(2) * precipSeconds/n
          albedo_new             =     albedo_new           +   snowResults(3) * precipSeconds/n


          !Correct if SWE is 0
          if(snowWaterEquiv_new <=  0.)   then
             snowWaterEquiv_new  =  0.
             snowEnergyCont_new  =  0.
             albedo_new          =  albedoMax
          end if

          flux_M_flow            =     flux_M_flow          +   snowResults(4)/n
          flux_M_subl            =     flux_M_subl          +   snowResults(5)/n
          flux_M_prec            =     flux_M_prec          +   snowResults(6)/n
          TEMP_MEAN              =     TEMP_MEAN            +   snowResults(7)/n
          TEMP_SURF              =     TEMP_SURF            +   snowResults(8)/n
          LIQU_FRAC              =     LIQU_FRAC            +   snowResults(9)/n
          flux_R_netS            =     flux_R_netS          +   snowResults(10)/n
          flux_R_netL            =     flux_R_netL          +   snowResults(11)/n
          flux_R_soil            =     flux_R_soil          +   snowResults(12)/n
          flux_R_sens            =     flux_R_sens          +   snowResults(13)/n
          stoi_f_prec            =     stoi_f_prec          +   snowResults(14)/n
          stoi_f_subl            =     stoi_f_subl          +   snowResults(15)/n
          stoi_f_flow            =     stoi_f_flow          +   snowResults(16)/n
          rate_G_alb             =     rate_G_alb           +   snowResults(17)/n

       END DO

       if(snowWaterEquiv_new > 0.) then
         precipModif = min((flux_M_flow * 1000 * precipSeconds), (snowWaterEquivalent*1000 + precipSumMM)) ! mm/referenceInterval
       end if

       if(snowWaterEquiv_new <= 0.) then
         precipModif = snowWaterEquivalent*1000 + precipSumMM
       end if

       snowCov        =  snowDepl(snowWaterEquiv_new)
       cloudFraction  =  cloudCoverage

    else
       snowResults(1:17) = f_snowModel(precipSumMM, sum(shortRad), tempAir(1), pressAir, relHumid, windSpeed, cloudCoverage, &
                                       precipSeconds, a0, a1, kSatSnow, densDrySnow, SpecCapRet, emissivitySnowMin, &
                                       emissivitySnowMax, tempAir_crit, albedoMin, albedoMax, agingRate_tAirPos, &
                                       agingRate_tAirNeg, soilDepth, soilDens, soilSpecHeat, weightAirTemp, &
                                       snowEnergyCont_new, snowWaterEquiv_new, albedo_new)
                                       !Calculations here only to get values for stochiometric factors and cloud coverage
       snowEnergyCont_new     =     0.
       snowWaterEquiv_new     =     0.
       albedo_new             =     albedoMax

       flux_M_flow            =     0.
       flux_M_subl            =     0.
       flux_M_prec            =     0.
       TEMP_MEAN              =     -9999.
       TEMP_SURF              =     -9999.
       LIQU_FRAC              =     1.
       flux_R_netS            =     0.
       flux_R_netL            =     0.
       flux_R_soil            =     0.
       flux_R_sens            =     0.
       stoi_f_prec            =     snowResults(14)
       stoi_f_subl            =     snowResults(15)
       stoi_f_flow            =     snowResults(16)
       rate_G_alb             =     0.

       snowCov                =     0.
       cloudFraction          =     cloudCoverage

       precipModif = precipSumMM !No snow cover present; rain not modified; when executing this line precipSumMM always 0. Necessary to  overwrite value of previous year of this time step

    end if

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
    function snowLiquidFrac(snowEnergyContent, snowWaterEquivalent)
        real :: snowEnergyContent
        real :: snowWaterEquivalent
        real :: snowLiquidFrac

        if (snowWaterEquivalent > 0.) then
        snowLiquidFrac = min(1. , max(0., snowEnergyContent / (snowWaterEquivalent * 333500.)))
        else
        snowLiquidFrac=1.
        endif

    end function snowLiquidFrac


    !Mean temperatur of the snow pack
    !Unit of result; °C (Range: -Inf...0)
    function snowTemp_mean(snowEnergyContent, snowWaterEquivalent, soilDepth, soilDens, soilSpecHeat)
       real :: snowEnergyContent
       real :: snowWaterEquivalent
       real :: soilDepth
       real :: soilDens
       real :: soilSpecHeat
       real :: snowTemp_mean

       if (snowWaterEquivalent > 0.) then
          !If the snow pack is free of liquid water
          if (snowEnergyContent < 0.) then
          !2090. = 1000. * 2.09 = WaterDensity (kg/m3) * specHeatCapIce (kJ/kg/K)
              snowTemp_mean = snowEnergyContent / (snowWaterEquivalent * 2090. + soilDepth * soilDens * soilSpecHeat)
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
    function R_netS(snowWaterEquivalent, shortRad, albedo)
       real :: snowWaterEquivalent
       real :: shortRad
       real :: albedo
       real :: R_netS

        if (snowWaterEquivalent > 0.) then
           R_netS = shortRad * (1 - albedo)
        else
           R_netS = 0.
        endif

     end function R_netS


    !Long-wave radiation balance
    !Unit of result: W/m2
    function R_netL(tempSnow_surf, emissivitySnowMin, emissivitySnowMax, tempAir, relHumid, cloudCoverage, &
                    albedo, albedoMin, albedoMax, snowWaterEquivalent)
       real :: tempSnow_surf
       real :: emissivitySnowMin
       real :: emissivitySnowMax
       real :: tempAir
       real :: relHumid
       real :: cloudCoverage
       real :: albedo
       real :: albedoMin
       real :: albedoMax
       real :: snowWaterEquivalent
       real :: R_outL
       real :: R_inL_clear
       real :: R_inL_cloud
       real :: R_netL

       if (snowWaterEquivalent > 0.) then
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
    function R_sens(tempSnow_surf, tempAir, pressAir, windSpeed, a0, a1, snowWaterEquivalent)
       real :: tempSnow_surf
       real :: tempAir
       real :: pressAir
       real :: windSpeed
       real :: a0
       real :: a1
       real :: snowWaterEquivalent
       real :: R_sens

       if (snowWaterEquivalent > 0.) then
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
    function M_subl(tempSnow_surf, tempAir, pressAir, relHumid, windSpeed, a0, a1, snowWaterEquivalent)
       real :: tempSnow_surf
       real :: tempAir
       real :: pressAir
       real :: relHumid
       real :: windSpeed
       real :: a0
       real :: a1
       real :: snowWaterEquivalent
       real :: M_subl

       if (snowWaterEquivalent > 0.) then
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
    function M_flow(snowLiquidFrac, kSatSnow, densDrySnow, specCapRet, snowWaterEquivalent)
       real :: snowLiquidFrac
       real :: kSatSnow
       real :: densDrySnow
       real :: specCapRet
       real :: snowWaterEquivalent
       real :: M_flow
       real :: rss

       if(snowWaterEquivalent > 0.) then
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
                   agingRate_tAirPos, agingRate_tAirNeg, snowWaterEquivalent)
       real :: albedo
       real :: precipSumMM
       real :: precipSeconds
       real :: tempAir
       real :: tempAir_crit
       real :: albedoMin
       real :: albedoMax
       real :: agingRate_tAirPos
       real :: agingRate_tAirNeg
       real :: snowWaterEquivalent
       real :: G_alb

       if(snowWaterEquivalent > 0.) then
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


    !Snow depletion curve
    !Function to describe the variability of snow water equivalent within a TC

    function snowDepl(snowWaterEquivalent)
        !no depeltion curve included yet; for now all or nothing
        real :: snowWaterEquivalent
        real :: snowDepl

        if(snowWaterEquivalent > snowFracThresh) then
           snowDepl = 1.
        else
           snowDepl = 0.
        end if
    end function snowDepl


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
    !Function to calculate derivates of the snow model´s state variables with respect to time
    !and to optionally return fundamental variables of the snow model (for debugging purposes and
    !in depth analysis of model behaviour
    !--------------------------------------------------------------------------------------------

    function f_snowModel(precipSumMM, shortRad, tempAir, pressAir, relHumid, windSpeed, cloudCoverage, &
                              precipSeconds, a0, a1, kSatSnow, densDrySnow, SpecCapRet, emissivitySnowMin, &
                              emissivitySnowMax, tempAir_crit, albedoMin, albedoMax, agingRate_tAirPos, &
                              agingRate_tAirNeg, soilDepth, soilDens, soilSpecHeat, weightAirTemp, &
                              snowEnergyContent, snowWaterEquivalent, albedo) result(snowModelRes)
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
        real :: snowEnergyContent
        real :: snowWaterEquivalent
        real :: albedo


        !Outputs
        real :: TEMP_MEAN
        real :: TEMP_SURF
        real :: LIQU_FRAC

        real :: flux_M_prec
        real :: flux_M_subl
        real :: flux_M_flow

        real :: flux_R_netS
        real :: flux_R_netL
        real :: flux_R_soil
        real :: flux_R_sens

        real :: stoi_f_prec
        real :: stoi_f_subl
        real :: stoi_f_flow

        real :: rate_G_alb
        real :: ddt_sec
        real :: ddt_swe
        real :: ddt_alb
        real, dimension(1:17) :: snowModelRes

        !Derived variables
        TEMP_MEAN = snowTemp_mean(snowEnergyContent, snowWaterEquivalent, soilDepth, soilDens, soilSpecHeat)
        TEMP_SURF = snowTemp_surf(TEMP_MEAN, tempAir, weightAirTemp)
        LIQU_FRAC = snowLiquidFrac(snowEnergyContent, snowWaterEquivalent)

        !Mass fluxes
        flux_M_prec = M_prec(precipSumMM, precipSeconds)
        flux_M_subl = M_subl(TEMP_SURF, tempAir, pressAir, relHumid, windSpeed, a0, a1, snowWaterEquivalent)
        flux_M_flow = M_flow(LIQU_FRAC, kSatSnow, densDrySnow, specCapRet, snowWaterEquivalent)

        !if no snow cover present and precipitation liquid, no addition to swe
        if(snowWaterEquivalent <= 0.0 .and. tempAir > tempAir_crit) then
           flux_M_prec = 0.
        end if

        !no more melt outflow than in snow cover
        if(flux_M_flow*precipSeconds > snowWaterEquivalent)then
           flux_M_flow = snowWaterEquivalent / precipSeconds
        endif


        !Radiation fluxes
        flux_R_netS = R_netS(snowWaterEquivalent, shortRad, albedo)
        flux_R_netL = R_netL(TEMP_SURF, emissivitySnowMin, emissivitySnowMax, tempAir, relHumid, &
                            cloudCoverage, albedo, albedoMin, albedoMax, snowWaterEquivalent)
        flux_R_soil = R_soil()
        flux_R_sens = R_sens(TEMP_SURF, tempAir, pressAir, windSpeed, a0, a1, snowWaterEquivalent)

        !Stochiometry factors
        stoi_f_prec = f_prec(tempAir, tempAir_crit)
        stoi_f_subl = f_subl()
        stoi_f_flow = f_flow()

        !Other rates
        rate_G_alb = G_alb(albedo, precipSumMM, precipSeconds, tempAir, tempAir_crit, albedoMin, &
                          albedoMax, agingRate_tAirPos, agingRate_tAirNeg, snowWaterEquivalent)

        !Computation of derivatives
        ddt_sec = 0.001 * (R_netS(snowWaterEquivalent, shortRad, albedo) + &
                           R_netL(TEMP_SURF, emissivitySnowMin, emissivitySnowMax, tempAir, relHumid, &
                           cloudCoverage, albedo, albedoMin, albedoMax, snowWaterEquivalent) + &
                           R_soil() + &
                           R_sens(TEMP_SURF, tempAir, pressAir, windSpeed, a0, a1, snowWaterEquivalent)) + &
                           f_prec(tempAir,tempAir_crit) * flux_M_prec - &
                           f_subl() * flux_M_subl - &
                           f_flow() * flux_M_flow

        ddt_swe = flux_M_prec - flux_M_subl - flux_M_flow

        ddt_alb = G_alb(albedo, precipSumMM, precipSeconds, tempAir, tempAir_crit, albedoMin, albedoMax, &
                        agingRate_tAirPos, agingRate_tAirNeg, snowWaterEquivalent)

        !Collect output
        snowModelRes = (/ ddt_sec,     ddt_swe,     ddt_alb,     flux_M_flow, flux_M_subl, flux_M_prec, &
                          TEMP_MEAN,   TEMP_SURF,   LIQU_FRAC,   flux_R_netS, flux_R_netL, flux_R_soil, &
                          flux_R_sens, stoi_f_prec, stoi_f_subl, stoi_f_flow, rate_G_alb/)

     end function f_snowModel


end subroutine snow_compute
