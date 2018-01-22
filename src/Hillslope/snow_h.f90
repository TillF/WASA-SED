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
    real :: lat,lon                                !Latitude / Longitude of centre of study area
    real :: radiation_tc(24)                       !hourly radiation for current TC corrected for slope and aspect
    real :: temperature_tc(24)                     !hourly temperature from daily temperature for current TC modified according to altitude with simple daily dynamic
    real :: tempLaps                               !Temperature lapse rate for modification depending on elevation of TC (°C/m)
    real :: tempAmplitude                          !Temperature amplitude to simulate daily cycle (°C)
    real :: tempMaxOffset                          !Offset of daily temperature maximum from 12:00 (h)
    real :: snowFracThresh                         !Threshold to determine when TC snow covered (m)

    real, pointer :: snowEnergyCont(:,:,:)         !Snow energy content [kJ/m²]
    real, pointer :: snowWaterEquiv(:,:,:)         !Snow water equivalent [m]
    real, pointer :: snowAlbedo(:,:,:)             !Albedo [-]

    real, pointer :: snowCover(:,:,:)              !Snow cover [-]
    real, pointer :: precipMod(:,:,:)              !Precipitation signal modified by snow module [mm]
    real, pointer :: radiMod(:,:,:)                !Radiation signal corrected for aspect and slope [W/m²]
    real, pointer :: temperaMod(:,:,:)             !Height-modified temperature signal [°C]
    real, pointer :: cloudFrac(:,:,:)              !Cloud fraction [-]
    real, pointer :: rel_elevation(:)              !Relative elevation of TC to mean subbasin [m]

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
    real, pointer :: lu_aspect(:)                  !mean aspect of a LU in radiants (south=0°, north=180° and west=90°, east=-90)
    real, pointer :: lu_alt(:)                     !mean LU altitude over subbasin outlet (m)
    
    !Indices for optional output arrays
    !When output activated only with dimensions (1,1,1) allocated
    !Depending on logical factor of outfiles.dat those indicy-arrays modified with function fLog() (in snow_h)
    integer,dimension(3) :: radiModIndices
    integer,dimension(3) :: temperaModIndices
    integer,dimension(3) :: cloudFracIndices
    integer,dimension(3) :: snowTempIndices
    integer,dimension(3) :: surfTempIndices
    integer,dimension(3) :: liquFracIndices
    integer,dimension(3) :: fluxPrecIndices
    integer,dimension(3) :: fluxSublIndices
    integer,dimension(3) :: fluxFlowIndices
    integer,dimension(3) :: fluxNetSIndices
    integer,dimension(3) :: fluxNetLIndices
    integer,dimension(3) :: fluxSoilIndices
    integer,dimension(3) :: fluxSensIndices
    integer,dimension(3) :: stoiPrecIndices
    integer,dimension(3) :: stoiSublIndices
    integer,dimension(3) :: stoiFlowIndices
    integer,dimension(3) :: rateAlbeIndices

contains

    !Modification of meteo-drivers according to time and location
    SUBROUTINE snow_prepare_input(hh, day, sb_counter, lu_counter2, tc_counter2, prec_mod, temp_mod, rad_mod, cloudFraction, &
                                  lapse_temp, tempAmpli, tempMaxOff)
        
        implicit none

        integer, intent(IN)                  :: hh, day, sb_counter, lu_counter2, tc_counter2
        real,    intent(INOUT)               :: prec_mod, temp_mod, rad_mod
        real,    intent(INOUT)               :: cloudFraction   !cloudiness fraction
        real,    intent(IN)                  :: lapse_temp      !Temperature lapse rate (°C/m)
        real,    intent(IN)                  :: tempAmpli       !Daily temperature amplitude (°C)
        real,    intent(IN)                  :: tempMaxOff      !Offset daily temperature maximum from 12:00 (h)
        real                                 :: lapse_prec = 0. !Lapse rate for precipitation modification according to elevation
        
        real                                 :: slope_rad, aspect  !slope, aspect of TC [rad]
        integer                              :: lu_id,i
        real,parameter                       :: S_C = 1366.944 !Solar constant [W/(m^2)]
        real                                 :: E_0     !Eccentricity of Earth orbit
        real                                 :: Gamma   !day angle [rad]
        real                                 :: delta   !Declination [rad] (the angular position of the sun at solar noon with respect to the plane of the equator
        real                                 :: omega_s !Local sunrise hour angle for a horizontal surface [rad]
        real                                 :: H0      !Daily extraterrestrial radiation [MJ/m2/d]
        real                                 :: K_t     !Ratio of global to extraterrestrial radiation [-]
        real                                 :: K_r     !Ratio of diffuse radiation to global radiation for a horizontal surface [-]
        real                                 :: f_beta  !Slope reduction factor [-]
        real, DIMENSION(1:24)                :: omega1  !Hour angle [rad]
        real, DIMENSION(1:24)                :: i_hor   !Angle of incidence onto horizontal surface [rad]
        real, DIMENSION(1:24)                :: i_slope !Angle of incidence onto inclinated surface [rad]
        real, DIMENSION(1:24)                :: I0      !Hourly extraterrestrial radiation [W/m²]
        real, DIMENSION(1:24)                :: G_h     !Hourly global radiation onto horizontal plane, rescaled from daily measurement [W/m²]
        real, DIMENSION(1:24)                :: G_B     !Direct radiation onto inclinated plane [W/m²]
        real, DIMENSION(1:24)                :: G_D     !Diffuse irradiance on inclinated surface [W/m²]


        !Modification precipitaion using simple, lapse-rate-based approach
           prec_mod = prec_mod + lapse_prec * rel_elevation(tcallid(sb_counter, lu_counter2, tc_counter2))

        !Modification temperature using simple, lapse-rate-based approach, if activated
        if (do_alt_corr) then
            if(dohour) then!if hourly run only modification according to altitude
               temp_mod = temp_mod + lapse_temp * rel_elevation(tcallid(sb_counter, lu_counter2, tc_counter2))
               temperature_tc(1:24) = temp_mod !all 24 values equal; in snow compute only first taken when hourly

            else !if daily altitudinal modification + daily cycle
               temp_mod = temp_mod + lapse_temp * rel_elevation(tcallid(sb_counter, lu_counter2, tc_counter2))

               DO i=1,24
                  temperature_tc(i) = temp_mod + sin(-pi/2 + (i-tempMaxOff)*(abs(-pi/2)+abs(3*pi/2))/24) * tempAmpli/2
               END DO

            end if
        else
          temperature_tc(1:24) = temp_mod !no modification all 24 values original daily mean
        end if
        
        !Modification of radiation according to aspect and slope if activated
           !References:
           !Tian et al. 2001: Estimating solar radiation on slopes of arbitrary aspec (https://doi.org/10.1016/S0168-1923(01)00245-3)
           !Maleki et al. 2017:Estimation of Hourly, Daily and Monthly Global Solar Radiation on Inclined Surfaces: Models Re-Visited  (doi:10.3390/en10010134)

           !ii: many of these caculations include computationally-demanding trigonometric functions
           !wherever possible, these should be done once at the start of the model run, or once per day and stored for re-use
        
           if (dohour)then
              write(*,*)'Radiation modification for aspect and slope for hourly resolution runs not yet integrated!'
           end if


           if (hh==1) then !if this is the first hour, compute values for the entire day (speedup)
              lu_id = id_lu_intern(lu_counter2, sb_counter) !get LU-id

              slope_rad = atan(slope(id_terrain_intern(tc_counter2, lu_id))) !slope angle [radiant]
              aspect = lu_aspect(lu_id) !aspect [rad]; (south=0, north=+/-pi and west=pi/2, east=-pi/2); in lu2.dat given in [°], converted to [rad] during read in
  
              E_0 = 1+0.033*cos(2.*julian_day/365) !eccentricity of Earth orbit; Tian, A.2, Beckman in Maleki et al. 2017
              Gamma = (2*pi*(julian_day-1)) / 365  !day angle [rad], Tian A.4
    
              ! [rad] declination (the angular position of the sun at solar noon with respect to the plane of the equator
              delta = (0.006918-0.399912*cos(Gamma)+0.070257*sin(Gamma)-0.002697*cos(2*Gamma) +0.000907*sin(2*Gamma)-0.002697*cos(3*Gamma)+0.00148*sin(3*Gamma))  !Tian, A.3
 
              omega_s = acos(-tan(lat)*tan(delta)) !local sunrise hour angle for a horizontal surface [rad], (Tian, A.5)

              !daily extraterrestrial radiation [MJ/m2/d] (Maleki 2017), p.20
              H0 =  (24*3.6/pi) * S_C * E_0 * (sin(lat)*sin(delta)*omega_s + cos(lat)*cos(delta)*sin(omega_s))/1000


              !hourly extraterrestrial radiation [W/m²] (Maleki 2017, p.4)

              !computation of solar time - omitted here, assuming that this is sufficiently presented by local time
              !lt1 = 1:24 !begin of hourly time step
              !B = 2*pi* (julian_day-81)/365
              !et = 9.87*sin(2*B)-7.53*cos(B)-1.5*sin(B) !equation of time (Tasdemiroglu 1988, wrong in Maleki (2017), taken from Bakirci, 2009 (eq.8) or Jamil & Khan, 2014(eq. 8)
              !L_s = 0 !standard meridian for local zone
              !L_l = 0 !longitude of location
              !st1 = lt1 + et/60+4/60*(L_s-L_l) !local solar time
              !omega1 = 15*(12-st1) * pi/180   !hour angle [rad]; angular distance between the observer’s meridian and the meridian whose plane contains the sun (Maleki 2017, p.3)

              do i=0,23
                omega1(i+1) = 15*(12-i) * pi/180   !hour angle [rad]; angular distance between the observer’s meridian and the meridian whose plane contains the sun (Maleki 2017, p.3)
              end do
    
              !angle of incidence onto horizontal surface [rad], Allen, 2006, eq.3
              !ii: should only be computed once per day for whole catchment
              i_hor = acos( &
                 sin(delta)*sin(lat)*1 &
                 - 0 &
                 + cos(delta)*cos(lat)*1*cos(-omega1) &
                 +0 &
                 +0 &
                )
      
              !angle of incidence onto inclinated surface [rad], Allen, 2006, eq.3. (hourly angle omega1 with minus, as calculations base on approach from Maleki with moring/afternoon angle defined with opposite sign
              !ii: should re-use terms of i_hor
               i_slope = acos( &
                 sin(delta)*sin(lat)*cos(slope_rad) &
               - sin(delta)*cos(lat)*sin(slope_rad)*cos(aspect) &
               + cos(delta)*cos(lat)*cos(slope_rad)*cos(-omega1) &
               + cos(delta)*sin(lat)*sin(slope_rad)*cos(aspect)*cos(-omega1) &
               + cos(delta)*sin(aspect)*sin(slope_rad)*sin(-omega1) &
               )
  
              I0 = S_C * E_0 * cos (i_hor)    !hourly extraterrestrial radiation [W/m²]
              where(abs(i_hor)>= pi/2) I0 = 0 !no radiation before/after sunrise/sunset
  
              !H0 = sum(I0) *3600/1e6  !daily extraterrestrial radiation [MJ/m²/d]

              K_t = rad_mod*3600*24/1e6 / H0 !ratio of global to extraterrestrial radiation (an inverse index of cloudiness)
              ! "clearness index" M_t in Maleki et al, 2017
              if (K_t >= 1) then
                  write(*,'(A,i0,a)')'Warning: Radiation for subbasin ', id_subbas_extern(sb_counter),' is higher than the computed extraterrestrial radiation. Truncated.'
                  K_t = 0.99
              end if

         if (do_rad_corr) then !Calculations needed before, as K_t used as measure for cloudiness

              G_h = I0 * K_t     !hourly global radiation onto horizontal plane, rescaled from daily measurement [W/m²]
  
    
              !K_r = a*K_t +aspect !ratio of diffuse radiation to global radiation for a horizontal surface [-], eq. 1 in Tian
  
              !Tian et al give no values for the coefficents. We consulted Maleki et al and selected Reindl et al, because it was linear and covered whole range of K_t
              if(K_t<=0.3) then
                  K_r = 1.02-0.248*K_t 
              else
                if (K_t < 0.78) then 
                    K_r = 1.45-1.67 * K_t 
                else
                  K_r = 0.147
                end if  
              end if
              K_r = min(1., max(0., K_r)) !restrain to 0..1
  
 
              G_B     = S_C * E_0 * K_t * (1-K_r)*cos (i_slope) !direct radiation onto inclinated plane [W/m²] Hay & McKay, eq. 1 [W/m²]
              where ((abs(i_slope)>= pi/2 .OR. abs(i_hor)>= pi/2 )) G_B = 0.     !no direct radiation before/after sunrise/sunset
   
              f_beta = 1-slope_rad/pi !slope reduction factor [-], Tian
              G_D = G_h * (f_beta*K_r + 0.2 * (1-f_beta)) !diffuse irradiance on inclinated surface
              radiation_tc(:) = G_B + G_D  !total irradiance onto inclinated surface
              rad_mod = sum(radiation_tc)/24

           end if !end first hour

        else !If radiation correction inactiv, distribute daily value equally over 24 time steps

           radiation_tc(:)=rad_mod/24
      
        end if !do_radCorr (radiation correction)

        !Calculations cloud fraction
        !Old method based on approach included Andreas Güntner (see etp_max.f90)
        !According to Shuttleworth (1992) Handbook of Hydrology, Chapter 4
        !IMPORTANT: nn does not equal cloudFrac => nn = 1-cloudFrac
           !cloudFraction = 1 - (rad_mod/radex(day)/0.55-0.18/0.55)
           !cloudFraction = MAX(0.0,cloudFraction)
           !cloudFraction = MIN(1.0,cloudFraction)
           !Angstrom coefficients:
           !0.18 (fraction of extratesetrial radiation on overcast day)
           !0.55+0.18 (fraction of extraterestrial radiation on clear days)
        !New approach using calculated daily extratesetrial radiation
          cloudFraction = 1 - K_t



    END SUBROUTINE snow_prepare_input

    !Function to addapt indicis of optional output arrays
    !When output not activated arrays only with dimensions (1,1,1) allocated and overwritten each run
    function fLog(factorLogic, day, hour, tc)

       logical                 :: factorLogic !Logical factor read from outfiles.dat
       integer                 :: day         !Current day
       integer                 :: hour        !Current timestep
       integer                 :: tc          !Currtent TC
       integer,dimension(3)    :: fLog        !Array to collect new indices

       if(.not. factorLogic)then !if optional output not activated, indices (1,1,1)
          fLog(:)=1
       else                      !if activated keep current indices and collect values
          fLog(1) = day
          fLog(2) = hour
          fLog(3) = tc
       end if

    end function fLog



end module snow_h
