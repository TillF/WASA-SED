SUBROUTINE readgen(path2do_dat)

    ! Code converted using TO_F90 by Alan Miller
    ! Date: 2005-06-30  Time: 13:47:18

    use allocate_h
    use climo_h
    use common_h
    use hymo_h
    use params_h
    use time_h
    use lake_h
    use routing_h
    use reservoir_h
    use erosion_h
    use utils_h
    use snow_h

    IMPLICIT NONE
    CHARACTER (LEN=160) :: path2do_dat		!Till: path to central control file (do.dat)

    INTEGER :: i,istate !,imun,imicro,imeso
    CHARACTER (LEN=150) :: custompath
    CHARACTER (LEN=250) :: dummy
    CHARACTER (LEN=150) :: dummy2


    if (trim(path2do_dat)=='') then
        path2do_dat='./Input/do.dat'    !Till: use default, if no command line argument was specified
        custompath=''
    else
        write(*,*)'reading runtime parameters from ',path2do_dat
        i=len_trim(path2do_dat)

        do while (i>0)
            if ((path2do_dat(i:i)=='/') .OR. (path2do_dat(i:i)==achar(92))) then	!find last slash or backslash in path
                exit
            end if
            i=i-1
        end do
        custompath=path2do_dat(1:i)	!extract path to do.dat (without filename)
    end if

    ! read run-time parameters
    OPEN(11,FILE=trim(path2do_dat) ,IOSTAT=istate,STATUS='old')
    IF (istate/=0) THEN
        write(*,'(A)')'Error: Control file '//trim(path2do_dat)//' could not be opened.'
        stop
    END IF
	

    READ(11,*)
    READ(11,'(a)') pfadp
    READ(11,'(a)') pfadn
    READ(11,*) tstart
    READ(11,*) tstop

    READ(11,'(A)') dummy !READ mstart (optional: dstart)
    READ(dummy,*,IOSTAT=i)mstart,dstart
    IF (i/=0 .OR. dstart==0.) THEN	!no dstart specified, assume 1
        READ(dummy,*)mstart
        dstart=1
    END IF

    READ(11,'(A)') dummy !READ mstop (optional: dstop)
    READ(dummy,*,IOSTAT=i)mstop,dstop
    IF (i/=0 .OR. dstop==0.) THEN	!no dstop specified, assume 31
        READ(dummy,*)mstop
        dstop=31
    END IF

!check time specifications    
    dummy = ""
    if (mstart<1 .OR. mstart>12) dummy=       ' ERROR: Invalid specification of start month.'
    if (mstop <1 .OR. mstop >12) dummy=trim(dummy)//' ERROR: Invalid specification of end month.'
    if (dstart<1 .OR. dstart>31) dummy=trim(dummy)//' ERROR: Invalid specification of start day.' 
    if (dstop<1  .OR. dstop >31) dummy=trim(dummy)//' ERROR: Invalid specification of end day.'
    if ( (tstart  > tstop) .OR. &
        ((tstart == tstop) .AND. (mstart  > mstop)) .OR. &
        ((tstart == tstop) .AND. (mstart == mstop) .and. (dstart > dstop) )) dummy = trim(dummy)//' ERROR: Simulation start must be before simulation end.'
    
    if (dummy /="") then
        write(*,*)trim(dummy)
        stop
    end if

    READ(11,*) subasin  !total no. of sub-basins
    READ(11,*) ntcinst  !total no. of sub-basin / SO / TC combinations
    READ(11,*) nsoter  !total no. of SOTER units / LUs in study area
    READ(11,*) nterrain  !total no. of terrain components in study area
    READ(11,*) nsoil  !total no. of soil components in study area
    READ(11,*) nveg   !total no. of vegetation units in study area
    READ(11,*) doreservoir !do large reservoir calculation (George's modules)
    READ(11,*) doacud
    READ(11,*) dolattc
    READ(11,*) doalllattc
    READ(11,*) dolatsc
    READ(11,*) dolatscsub
    READ(11,*) dotrans

    READ(11,*) dohour !overwritten by value of dt
    READ(11,*) scenario
    READ(11,*) krig

    READ(11,'(A)') dummy !READ kfkorr (optional: kfkorr_a and kfkorr_b)
    READ(dummy,*,IOSTAT=i)kfkorr,kfkorr_a,kfkorr_b
    IF (i/=0 .OR. kfkorr_a==0.) THEN	!no parameters for variable kfcorr specified, use constant kfcorr (as in old version)
        READ(dummy,*)kfkorr
        kfkorr_a=0.
        kfkorr_b=1.
    END IF

    READ(11,*) intcf
    READ(11,*) dointc
    READ(11,*) doscale
    READ(11,*) domuncell
    READ(11,*) sensfactor
    READ(11,*) dt
    READ(11,*) dosediment
    READ(11,*) n_sed_class
    READ(11,*)  !line ignored
    READ(11,*) river_transport
    READ(11,*) reservoir_transport
    
    READ(11,'(A)', IOSTAT=i) dummy !READ doloadstate
    IF (i==0) then  
        READ(dummy,*,IOSTAT=i) doloadstate, append_output !try to read doloadstate AND append_output
        IF (i/=0 ) READ(dummy,*,IOSTAT=i) doloadstate !read doloadstate only
    END IF
    
    !READ(11,*, IOSTAT=i) dosavestate
    READ(11,'(A)', IOSTAT=i) dummy !READ dosavestate
    IF (i==0) then  
        READ(dummy,*,IOSTAT=i) dosavestate, save_states_yearly !try to read dosavestate AND save_states_yearly
        IF (i/=0 ) READ(dummy,*,IOSTAT=i) dosavestate !read dosavestate only
    END IF
    
    
    READ(11,*, IOSTAT=i) dosnow !ii: rather use existence of snow input files as indicator

    CLOSE(11)

    if (dt/=1 .AND. dt /=24) then
		write(*,*)'ERROR: dt must be set to 1 or 24 (do.dat).'
        stop
    end if
	dohour = dt == 1
	
	if (trim(custompath)/='') then		!if a custom path was specified, all paths are relative to this one
        pfadp=trim(custompath)//pfadp
        pfadn=trim(custompath)//pfadn
    end if


    pfadi=len_trim(pfadn)	! determine length of path name for output files
    pfadj=len_trim(pfadp)	! determine length of basic model path name
    ncaset=len_trim(caset)	! determine length of case study path name


    !Till: read maximum dimensions of arrays
    OPEN(11,FILE=pfadp(1:pfadj)// 'maxdim.dat',IOSTAT=i,STATUS='old')
    IF (i==0) THEN
        READ(11,*)
        READ(11,*)maxsoter
        READ(11,*)maxterrain
        READ(11,*)maxsoil
        READ(11,*)maxhori
        READ(11,*)ntrans
        READ(11,*, IOSTAT=i)nxsection_res
        READ(11,*, IOSTAT=i)npointsxsect
        CLOSE(11)
    else
        write(*,*)pfadp(1:pfadj)// 'maxdim.dat could not be opened. Using default dimensions specified in the source code.'
        maxsoter=7
        maxterrain=3
        maxsoil=28
        maxhori=8
        ntrans=2
        nxsection_res=200
        npointsxsect=200
    END IF

!read general erosion parameters
    if (dosediment) THEN
        allocate( spcon(n_sed_class), spexp(n_sed_class))
        
        spcon(:)=  0.016111		!0.0001-0.01  default values
        spexp (:)= 1.707			!1 - 1.5
        erosion_equation=0
        OPEN(11,FILE=pfadp(1:pfadj)// 'erosion.ctl',IOSTAT=istate,STATUS='old')
        IF (istate==0) THEN
            READ(11,'(a)',IOSTAT=istate)dummy
            do while (istate==0)
                READ(dummy,*,IOSTAT=istate)dummy2
                SELECT CASE (trim(dummy2))
                    CASE ('application_scale')
                        READ(dummy,*,IOSTAT=istate) dummy2, do_musle_subbasin
                    CASE ('erosion_equation')
                        READ(dummy,*,IOSTAT=istate) dummy2, erosion_equation
                    CASE ('ri_05_coeffs')
                        READ(dummy,*,IOSTAT=istate) dummy2, a_i30,b_i30
                    CASE ('transport_limit_mode')
                        READ(dummy,*,IOSTAT=istate) dummy2, transport_limit_mode
                    CASE ('transp_cap_a')
                        READ(dummy,*,IOSTAT=istate) dummy2, spcon(1)
                        spcon(:)=spcon(1)
                    CASE ('transp_cap_b')
                        READ(dummy,*,IOSTAT=istate) dummy2, spexp(1)
                        spexp(:)=spexp(1)
                END SELECT
                READ(11,'(a)',IOSTAT=istate)dummy
            end do
            CLOSE(11)
            IF ((erosion_equation>4.) .OR. (erosion_equation<0.)) THEN
                write(*,*)'WARNING: erosion_equation was outside [1..4], assumed to be 3 (MUSLE).'
                erosion_equation=3            !default erosion equation to be used: MUSLE
            END IF
        ELSE        !erosion.ctl not found
            write(*,*)'WARNING: erosion.ctl not found, using defaults.'
            erosion_equation=3            !default erosion equation to be used: MUSLE
            do_musle_subbasin=.FALSE.            !default 0: compute erosion on TC-scale
            transport_limit_mode=2 !transport capacity according to Everaert (1991)
        END IF
        
        !(taken from erosion.ctl, if present): default coefficients for estimation of maximum half-hour rainfall intensity (ri_05) from daily rainfall data (R_day) 
        if (a_i30==-1) then !scaling coefficients for half-hour-intensity not set, use defaults (ri_05=a*R_dt^b)
            if (dt==24) then
                a_i30=1.1630         !default coefficients for estimation of maximum half-hour rainfall intensity (ri_05) from daily rainfall data (R_day)    
                b_i30=0.667981            
            else
                a_i30=1.            
                b_i30=1.
            end if
        end if
    END IF !dosediments



    !Parameter for snow routine
    !Read, if present. Else using default values.
    if (dosnow /= 0) THEN

    precipSeconds     =    dt * 3600.
    !Read parameters for snow routine
     OPEN(12, file=pfadp(1:pfadj)// 'Hillslope/snow_params.ctl',IOSTAT=istate, STATUS='old')
         IF (istate==0) THEN
            READ(12,'(a)',IOSTAT=istate)dummy
            do while (istate==0)
                READ(12,'(a)',IOSTAT=istate)dummy
                READ(dummy,*)dummy2
                SELECT CASE (trim(dummy2))
                    CASE ('a0')
                        READ(dummy,*) dummy2, a0
                    CASE ('a1')
                        READ(dummy,*) dummy2, a1
                    CASE ('kSatSnow')
                        READ(dummy,*) dummy2, kSatSnow
                    CASE ('densDrySnow')
                        READ(dummy,*) dummy2, densDrySnow
                    CASE ('specCapRet')
                        READ(dummy,*) dummy2, specCapRet
                    CASE ('emissivitySnowMin')
                        READ(dummy,*) dummy2, emissivitySnowMin
                    CASE ('emissivitySnowMax')
                        READ(dummy,*) dummy2, emissivitySnowMax
                    CASE ('tempAir_crit')
                        READ(dummy,*) dummy2, tempAir_crit
                    CASE ('albedoMin')
                        READ(dummy,*) dummy2, albedoMin
                    CASE ('albedoMax')
                        READ(dummy,*) dummy2, albedoMax
                    CASE ('agingRate_tAirPos')
                        READ(dummy,*) dummy2, agingRate_tAirPos
                    CASE ('agingRate_tAirNeg')
                        READ(dummy,*) dummy2, agingRate_tAirNeg
                    CASE ('soilDepth')
                        READ(dummy,*) dummy2, soilDepth
                    CASE ('soilDens')
                        READ(dummy,*) dummy2, soilDens
                    CASE ('soilSpecHeat')
                        READ(dummy,*) dummy2, soilSpecHeat
                    CASE ('weightAirTemp')
                        READ(dummy,*) dummy2, weightAirTemp
                    CASE ('lat')
                        READ(dummy,*) dummy2, lat
                        lat = lat *pi/180
                    CASE ('lon')
                        READ(dummy,*) dummy2, lon
                        lon = lon *pi/180
                    CASE ('do_rad_corr')
                        READ(dummy,*) dummy2, do_rad_corr
                    CASE ('do_alt_corr')
                        READ(dummy,*) dummy2, do_alt_corr
                    CASE ('tempLaps')
                        READ(dummy,*) dummy2, tempLaps
                    CASE ('tempAmplitude')
                        READ(dummy,*) dummy2, tempAmplitude
                    CASE ('tempMaxOffset')
                        READ(dummy,*) dummy2, tempMaxOffset
                    CASE ('snowFracThresh')
                        READ(dummy,*) dummy2, snowFracThresh

                END SELECT
            end do
     CLOSE(12)

        ELSE        !snow.ctl not found
            write(*,*)'WARNING: snow_params.ctl not found, using defaults.'

                !Default values parameters for snow routine
                a0=0.002                       !Empirical coeff. (m/s)
                a1=0.0008                      !Empirical coeff. (-)
                kSatSnow = 0.00004             !Saturated hydraulic conductivity of snow (m/s)
                densDrySnow=450                !Density of dry snow (kg/m3)
                specCapRet=0.05                !Capill. retent. vol as fraction of solid SWE (-)
                emissivitySnowMin=0.84         !Minimum snow emissivity used for old snow (-)
                emissivitySnowMax=0.99         !Maximum snow emissivity used for new snow (-)
                tempAir_crit=0.2               !Threshold temp. for rain-/snowfall (°C)
                albedoMin=0.55                 !Minimum albedo used for old snow (-)
                albedoMax=0.88                 !Maximum albedo used for new snow (-)
                agingRate_tAirPos=0.00000111   !Aging rate for air temperatures > 0 (1/s)
                agingRate_tAirNeg=0.000000462  !Aging rate for air temperatures < 0 (1/s)
                soilDepth=0.1                  !Depth of interacting soil layer (m)
                soilDens=1300.                 !Density of soil (kg/m3)
                soilSpecHeat=2.18              !Spec. heat capacity of soil (kJ/kg/K)
                weightAirTemp=0.5              !Weighting param. for air temp. (-) in 0...1
                lat = 42.4                     !Latitude of centre of study area
                lon = 0.55                     !Longitude of centre of study area
                do_rad_corr = .TRUE.           !modification of radiation with aspect and slope
                do_alt_corr = .TRUE.           !modification of temperature with altitude of LU
                tempLaps = -0.006              !Temperature lapse rate for modification depending on elevation of TC (°C/m)
                tempAmplitude = 8              !Temperature amplitude to simulate daily cycle (°C])
                tempMaxOffset = 2              !Offset of daily temperature maximum from 12:00 (h)
                snowFracThresh = 0.02          !Threshold to determine when TC snow covered (m)

        END IF
    END IF !dosnow


!Till: read list of desired output files
    f_daily_actetranspiration=.FALSE.	!disable all output files
    f_daily_potetranspiration=.FALSE.
    f_daily_qhorton=.FALSE.
    f_daily_qin_m3s=.FALSE.
    f_daily_qout_m3s=.FALSE.
    f_daily_rain=.FALSE.
    f_daily_runoff=.FALSE.
    f_daily_sediment_production=.FALSE.
    f_daily_subsurface_runoff=.FALSE.
    f_daily_theta=.FALSE.
    f_daily_total_overlandflow=.FALSE.
    f_daily_water_subbasin=.FALSE.
    f_routing_response=.FALSE.
    f_sediment_production=.FALSE.
    f_water_subbasin=.FALSE.
    f_deep_gw_recharge=.FALSE.
    f_deep_gw_discharge=.FALSE.
    f_tc_theta=.FALSE.
    f_daily_gw_loss=.FALSE.
    f_river_degradation=.FALSE.
    f_river_deposition=.FALSE.
    f_river_flow=.FALSE.
    f_river_flow_dailyaverage=.FALSE.
    f_river_flowdepth=.FALSE.
    f_river_sediment_concentration =.FALSE.
    f_river_sediment_total=.FALSE.
    f_river_sediment_total_dailyaverage =.FALSE.
    f_river_storage=.FALSE.
    f_river_sediment_storage=.FALSE.
    f_river_susp_sediment_storage=.FALSE.
    f_river_velocity=.FALSE.
    f_river_bedload=.FALSE.
    f_river_infiltration=.FALSE.
    f_tc_surfflow=.FALSE.
    f_tc_sedout=.FALSE.
    f_lu_sedout=.FALSE.

    f_actetranspiration=.FALSE.
    f_qhorton=.FALSE.
    f_subsurface_runoff=.FALSE.
    f_total_overlandflow=.FALSE.
    f_gw_discharge=.FALSE.
    f_potetranspiration=.FALSE.
    f_gw_loss=.FALSE.
    f_gw_recharge=.FALSE.

    f_snowEnergyCont=.FALSE.
    f_snowWaterEquiv=.FALSE.
    f_snowAlbedo=.FALSE.
    f_snowCover=.FALSE.
    f_snowTemp=.FALSE.
    f_surfTemp=.FALSE.
    f_liquFrac=.FALSE.
    f_fluxPrec=.FALSE.
    f_fluxSubl=.FALSE.
    f_fluxFlow=.FALSE.
    f_fluxNetS=.FALSE.
    f_fluxNetL=.FALSE.
    f_fluxSoil=.FALSE.
    f_fluxSens=.FALSE.
    f_stoiPrec=.FALSE.
    f_stoiSubl=.FALSE.
    f_stoiFlow=.FALSE.
    f_rateAlbe=.FALSE.
    f_precipMod=.FALSE.
    f_cloudFrac=.FALSE.

    f_res_watbal=.FALSE.
    f_res_vollost=.FALSE.
    f_res_cav=.FALSE.
    f_res_hydraul=.FALSE.
    f_res_bedchange=.FALSE.
    f_res_sedbal=.FALSE.
    f_res_longitudunal=.FALSE.
    f_res_sedcomposition=.FALSE.
    f_lake_inflow_r=.FALSE.
    f_lake_outflow_r=.FALSE.
    f_lake_retention_r=.FALSE.
    f_lake_volume_r=.FALSE.
    f_lake_sedinflow_r=.FALSE.
    f_lake_sedoutflow_r=.FALSE.
    f_lake_sedretention_r=.FALSE.
    f_lake_sedimentation_r=.FALSE.
    f_lake_watbal=.FALSE.
    f_lake_sedbal=.FALSE.
    f_lake_inflow=.FALSE.
    f_lake_outflow=.FALSE.
    f_lake_volume=.FALSE.
    f_lake_retention=.FALSE.
    f_lake_vollost=.FALSE.
    f_lake_sedinflow=.FALSE.
    f_lake_sedoutflow=.FALSE.
    f_lake_sizedistoutflow=.FALSE.

    OPEN(11,FILE=pfadp(1:pfadj)// 'outfiles.dat',IOSTAT=istate,STATUS='old')
    IF (istate==0) THEN
        READ(11,*,IOSTAT=istate)dummy
        READ(11,*,IOSTAT=istate)dummy

        DO WHILE (istate==0)
            SELECT CASE (trim(locase(dummy)))		!enable/disable file output of desired results
                CASE ('daily_actetranspiration')
                    f_daily_actetranspiration=.TRUE.
                CASE ('daily_potetranspiration')
                    f_daily_potetranspiration=.TRUE.
                CASE ('daily_qhorton')
                    f_daily_qhorton=.TRUE.
                CASE ('daily_qin_m3s')
                    f_daily_qin_m3s=.TRUE.
                CASE ('daily_qout_m3s')
                    f_daily_qout_m3s=.TRUE.
                CASE ('daily_rain')
                    f_daily_rain=.TRUE.
                CASE ('daily_runoff')
                    f_daily_runoff=.TRUE.
                CASE ('daily_sediment_production')
                    f_daily_sediment_production = dosediment
                CASE ('daily_subsurface_runoff')
                    f_daily_subsurface_runoff=.TRUE.
                CASE ('daily_theta')
                    f_daily_theta=.TRUE.
                CASE ('daily_total_overlandflow')
                    f_daily_total_overlandflow=.TRUE.
                CASE ('daily_water_subbasin')
                    f_daily_water_subbasin=.TRUE.
                CASE ('routing_response')
                    f_routing_response=.TRUE.
                CASE ('sediment_production')
                    f_sediment_production=dosediment
                CASE ('water_subbasin')
                    f_water_subbasin=.TRUE.
                CASE ('deep_gw_recharge')
                    f_deep_gw_recharge=.TRUE.
                CASE ('deep_gw_discharge')
                    f_deep_gw_discharge=.TRUE.
                CASE ('f_daily_gw_loss')
                    f_daily_gw_loss=.TRUE.
                CASE ('tc_theta')
                    f_tc_theta=.TRUE.

                CASE ('snowenergycont')
                    f_snowEnergyCont=.TRUE. .AND. (dosnow /= 0)
                CASE ('snowwaterequiv')
                    f_snowWaterEquiv=.TRUE. .AND. (dosnow /= 0)
                CASE ('snowalbedo')
                    f_snowAlbedo=.TRUE. .AND. (dosnow /= 0)
                CASE ('snowcover')
                    f_snowCover=.TRUE. .AND. (dosnow /= 0)
                CASE ('snowtemp')
                    f_snowTemp=.TRUE. .AND. (dosnow /= 0)
                CASE ('surftemp')
                    f_surfTemp=.TRUE. .AND. (dosnow /= 0)
                CASE ('liqufrac')
                    f_liquFrac=.TRUE. .AND. (dosnow /= 0)
                CASE ('fluxprec')
                    f_fluxPrec=.TRUE. .AND. (dosnow /= 0)
                CASE ('fluxsubl')
                    f_fluxSubl=.TRUE. .AND. (dosnow /= 0)
                CASE ('fluxflow')
                    f_fluxFlow=.TRUE. .AND. (dosnow /= 0)
                CASE ('fluxnets')
                    f_fluxNetS=.TRUE. .AND. (dosnow /= 0)
                CASE ('fluxnetl')
                    f_fluxNetL=.TRUE. .AND. (dosnow /= 0)
                CASE ('fluxsoil')
                    f_fluxSoil=.TRUE. .AND. (dosnow /= 0)
                CASE ('fluxsens')
                    f_fluxSens=.TRUE. .AND. (dosnow /= 0)
                CASE ('stoiprec')
                    f_stoiPrec=.TRUE. .AND. (dosnow /= 0)
                CASE ('stoisubl')
                    f_stoiSubl=.TRUE. .AND. (dosnow /= 0)
                CASE ('stoiflow')
                    f_stoiFlow=.TRUE. .AND. (dosnow /= 0)
                CASE ('ratealbe')
                    f_rateAlbe=.TRUE. .AND. (dosnow /= 0)
                CASE ('precipmod')
                    f_precipMod=.TRUE. .AND. (dosnow /= 0)
                CASE ('cloudfrac')
                    f_cloudFrac=.TRUE. .AND. (dosnow /= 0)
                CASE ('radimod')
                    f_radiMod=.TRUE. .AND. (dosnow /= 0)
                CASE ('temperamod')
                    f_temperaMod=.TRUE. .AND. (dosnow /= 0)
                CASE ('rel_elevation')
                    f_rel_elevation=.TRUE. .AND. (dosnow /= 0)

                CASE ('river_degradation')
                    f_river_degradation=dosediment
                CASE ('river_deposition')
                    f_river_deposition=dosediment
                CASE ('river_flow')
                    f_river_flow=.TRUE.
                CASE ('river_flow_dailyaverage')
                    f_river_flow_dailyaverage=.TRUE.
                CASE ('river_flowdepth')
                    f_river_flowdepth=.TRUE.
                CASE ('river_sediment_concentration')
                    f_river_sediment_concentration =dosediment
                CASE ('river_sediment_total')
                    f_river_sediment_total=dosediment .AND. river_transport /= 1
                CASE ('river_sediment_total_dailyaverage')
                    f_river_sediment_total_dailyaverage =dosediment
                CASE ('river_storage')
                    f_river_storage=.TRUE.
                CASE ('river_sediment_storage')
                    f_river_sediment_storage=dosediment
                CASE ('river_susp_sediment_storage')
                    f_river_susp_sediment_storage=dosediment
                CASE ('river_velocity')
                    f_river_velocity=.TRUE.
                CASE ('river_bedload')
                    f_river_bedload=dosediment
                CASE ('river_infiltration')
                    f_river_infiltration=.TRUE. .AND. river_transport /= 1
                CASE ('tc_surfflow')
                    f_tc_surfflow=.TRUE.
                CASE ('tc_sedout')
                    f_tc_sedout=dosediment
                CASE ('lu_sedout')
                    f_lu_sedout=dosediment

                CASE ('actetranspiration')
                    f_actetranspiration=.TRUE.
                CASE ('qhorton')
                    f_qhorton=.TRUE.
                CASE ('subsurface_runoff')
                    f_subsurface_runoff=.TRUE.
                CASE ('total_overlandflow')
                    f_total_overlandflow=.TRUE.
                CASE ('gw_discharge')
                    f_gw_discharge=.TRUE.
                CASE ('potetranspiration')
                    f_potetranspiration=.TRUE.
                CASE ('gw_loss')
                    f_gw_loss=.TRUE.
                CASE ('gw_recharge')
                    f_gw_recharge=.TRUE.

                CASE ('res_watbal')
                    f_res_watbal=.TRUE.
                CASE ('res_vollost')
                    f_res_vollost=.TRUE.
                CASE ('res_cav')
                    f_res_cav=.TRUE.
                CASE ('res_hydraul')
                    f_res_hydraul=.TRUE.
                CASE ('res_bedchange')
                    f_res_bedchange=dosediment
                CASE ('res_sedbal')
                    f_res_sedbal=dosediment
                CASE ('res_longitudunal')
                    f_res_longitudunal=.TRUE.
                CASE ('res_sedcomposition')
                    f_res_sedcomposition=dosediment
                CASE ('lake_inflow_r')
                    f_lake_inflow_r=.TRUE.
                CASE ('lake_outflow_r')
                    f_lake_outflow_r=.TRUE.
                CASE ('lake_retention_r')
                    f_lake_retention_r=.TRUE.
                CASE ('lake_volume_r')
                    f_lake_volume_r=.TRUE.
                CASE ('lake_sedinflow_r')
                    f_lake_sedinflow_r=dosediment
                CASE ('lake_sedoutflow_r')
                    f_lake_sedoutflow_r=dosediment
                CASE ('lake_sedretention_r')
                    f_lake_sedretention_r=dosediment
                CASE ('lake_sedimentation_r')
                    f_lake_sedimentation_r=dosediment
                CASE ('lake_watbal')
                    f_lake_watbal=.TRUE.
                CASE ('lake_sedbal')
                    f_lake_sedbal=dosediment
                CASE ('lake_inflow')
                    f_lake_inflow=.TRUE.
                CASE ('lake_outflow')
                    f_lake_outflow=.TRUE.
                CASE ('lake_volume')
                    f_lake_volume=.TRUE.
                CASE ('lake_retention')
                    f_lake_retention=.TRUE.
                CASE ('lake_vollost')
                    f_lake_vollost=.TRUE.
                CASE ('lake_sedinflow')
                    f_lake_sedinflow=dosediment
                CASE ('lake_sedoutflow')
                    f_lake_sedoutflow=dosediment
                CASE ('lake_sizedistoutflow')
                    f_lake_sizedistoutflow=dosediment


            END SELECT
            READ(11,*,IOSTAT=istate)dummy	!try to read next line
        END DO
        CLOSE(11)
	
        f_qhorton= f_qhorton .OR. f_daily_qhorton

    ELSE
        WRITE(*,*)pfadp(1:pfadj)// 'outfiles.dat could not be opened. Using default output files.'
        !these are the default output files
        f_daily_runoff=.TRUE.
        f_daily_sediment_production=dosediment
        f_daily_water_subbasin=.TRUE.
        f_water_subbasin=.TRUE.
        f_river_flow=.TRUE.
        f_river_sediment_total=dosediment
        f_river_sediment_concentration=dosediment
    END IF
    !end insert Till


    
    
!allocate necessary memory
    nt = int(24/dt)	!Till: number of simulation steps per day

    call allocate_hymo()

    INCLUDE '../Hillslope/allocat_erosion.var'
    call allocate_routing()

        nt = int(24/dt)	!Till: number of simulation steps per day !ii: can be removed?
    !INCLUDE '../Reservoir/allocat_reservoir_lake.var'
    call allocate_reservoir()


    
    
    
! save summary of settings to output directory
    OPEN(11,FILE=pfadn(1:pfadi)//'parameter.out', STATUS='unknown', IOSTAT=istate)
    IF (istate/=0) THEN !output dir not found, try creating using system command
        write(*,*)'Warning: Output path ',pfadn(1:pfadi),' not found, trying to create...'
        dummy=pfadn(1:pfadi)
        !try create with *nix command
            do i=1,pfadi
                if (dummy(i:i)=="\") dummy(i:i)="/" !using slashes as delimiter (*nix)
            end do

        CALL system('mkdir '//dummy, istate)
	
        OPEN(11,FILE=pfadn(1:pfadi)//'parameter.out', STATUS='unknown',IOSTAT=istate) !try again to open
        if (istate/=0) then !try create with Windows command
            do i=1,pfadi
                if (dummy(i:i)=="/") dummy(i:i)="\" !using backslashes as delimiter (windows)
            end do
            CALL system('mkdir '//dummy)
		
            OPEN(11,FILE=pfadn(1:pfadi)//'parameter.out', STATUS='unknown',IOSTAT=istate)
            if (istate/=0) then
                write(*,*)'Error: Output file ',pfadn(1:pfadi)//'parameter.out',' could not be created, aborting.'
                stop
            end if
        end if
    end if
    
    WRITE(11,*)
    WRITE(11,'(a)') pfadp
    WRITE(11,'(a)') pfadn
    WRITE(11,*) 'start year of simulation: ', tstart
    WRITE(11,*) 'end year of simulation: ', tstop
    WRITE(11,*) 'start of simulation (month, day): ',mstart, dstart
    WRITE(11,*) 'end of simulation (month, day): ',mstop, dstop
    WRITE(11,*) 'no. of sub-basin: ',subasin
    WRITE(11,*) 'no. of combinations: ',ntcinst
    WRITE(11,*) 'total no. of SOTER units: ',nsoter
    WRITE(11,*) 'total no. of terrain components: ',nterrain
    WRITE(11,*) 'total no. of soil components: ',nsoil
    WRITE(11,*) 'total no. of vegetation units: ',nveg
    WRITE(11,*) 'do reservoir calculations: ',doacud
    WRITE(11,*) dolattc
    WRITE(11,*) doalllattc
    WRITE(11,*) dolatsc
    WRITE(11,*) dolatscsub
    WRITE(11,*) dotrans
    WRITE(11,*) dohour
    WRITE(11,*) scenario
    WRITE(11,*) krig
    WRITE(11,*) kfkorr,kfkorr_a,kfkorr_b
    WRITE(11,*) intcf
    WRITE(11,*) dointc
    WRITE(11,*) doscale
    WRITE(11,*) domuncell
    WRITE(11,*) sensfactor
    write(11,*) dosediment
    write(11,*) n_sed_class
    write(11,*) "[not used]"
    write(11,*) river_transport
    write(11,*) reservoir_transport
    if (dosediment) then
        WRITE(11,*)
        WRITE(11,*) 'spatial scale for application of erosion eq : ',do_musle_subbasin
        WRITE(11,*) 'erosion eq : ',erosion_equation
    end if
    WRITE(11,*) 'WASA model, ',trim(rev_string1),'; ',trim(rev_string2)
    CLOSE(11)
    

    RETURN

END SUBROUTINE readgen
