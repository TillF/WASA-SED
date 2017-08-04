module snow_params

    !Parameters for snow routine

    real :: precipSeconds=86400.           !Length of referenceInterval (seconds)
    real :: a0=0.002                       !Empirical coeff. (m/s)
    real :: a1=0.0008                      !Empirical coeff. (m/s)
    real :: kSatSnow = 0.00004             !Saturated hydraulic conductivity of snow (m/s)
    real :: densDrySnow=450                !Density of dry snow (kg/m3)
    real :: specCapRet=0.05                !Capill. retent. vol as fraction of solid SWE (-)
    real :: emissivitySnowMin=0.84          !Minimum snow emissivity used for old snow (-)
    real :: emissivitySnowMax=0.99          !Maximum snow emissivity used for new snow (-)
    real :: tempAir_crit=0.2               !Threshold temp. for rain-/snowfall (°C)
    real :: albedoMin=0.55                 !Minimum albedo used for old snow (-)
    real :: albedoMax=0.88                 !Maximum albedo used for new snow (-)
    real :: agingRate_tAirPos=0.12/86400   !Aging rate for air temperatures > 0 (1/s)
    real :: agingRate_tAirNeg=0.05/86400   !Aging rate for air temperatures < 0 (1/s)
    real :: soilDepth=0.1                  !Depth of interacting soil layer (m)
    real :: soilDens=1300.                 !Density of soil (kg/m3)
    real :: soilSpecHeat=2.18              !Spec. heat capacity of soil (kJ/kg/K)
    real :: weightAirTemp=0.5              !Weighting param. for air temp. (-) in 0...1

end module snow_params
