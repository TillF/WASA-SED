    SUBROUTINE soilwat(hh,day,month,i_subbas2,i_ce,i_lu,lu_counter2,tcid_instance2,id_tc_type2,tc_counter2,  &
        thact,thactroot,q_surf_in,q_surf_out,q_sub_in,q_sub_out,qgw,deepqgw,  &
        hortf,tcaet,tclai,  &
        tcsoilet,tcintc,prec_in,precday_in,prechall2_in,petday,  &
        tcarea2,bal, rootd_act,height_act,lai_act,alb_act,sed_in_tc,sed_out_tc, q_rill_in)


    ! Code converted using TO_F90 by Alan Miller
    ! Date: 2005-06-22  Time: 14:51:47

    ! (0) Initialization
    ! (1.1) Update of soil moisture due to lateral inflow of subsurface flow from other SVCs of this TC
    ! (1.2) Update of soil moisture due to lateral inflow of subsurface flow from higher TC
    ! (11) INTERCEPTION
    ! (2)   SURFACE RUNOFF
    !     (2.0) on rock outcrops: overland flow
    !     (2.1) on saturated part of TC: saturation excess overland flow occurs
    !     (2.2) on non-saturated part of TC: infiltration excess overland flow may occur
    ! (3) EVAPOTRANSPIRATION
    ! (4) Percolation and lateral flow
    ! (5) Capillary rise
    ! (6) End of calculations for terrain component


    use common_h
    use params_h
    use hymo_h
    use erosion_h
    use utils_h
    use climo_h
    use snow_h

    IMPLICIT NONE


    !real :: sedi_yield						!Till: declared function

    INTEGER, INTENT(IN)                  :: hh
    INTEGER, INTENT(IN)                  :: day
    INTEGER, INTENT(IN)                  :: month
    INTEGER, INTENT(IN)                  :: i_subbas2		!internal representation of i_subbas2in (id of subbasin)
    INTEGER, INTENT(IN)                  :: i_ce			!internal representation of i_subbas2in (id of subbasin)
    INTEGER, INTENT(IN)                  :: i_lu			!internal representation of isot (id of LU)
    INTEGER, INTENT(IN)                  :: lu_counter2	    !internal representation of lu_counter
    INTEGER, INTENT(IN)                  :: tcid_instance2		!(internal) id of TC-instance (unique subbas-LU-TC-combination)
    INTEGER, INTENT(IN)                  :: id_tc_type2			!ID of TC type
    INTEGER, INTENT(IN)                  :: tc_counter2
    REAL, INTENT(IN OUT)                 :: thact			!internal representation of thact ((average) actual water content of TC [mm])
    REAL, INTENT(OUT)                    :: thactroot
    REAL, INTENT(IN)                     :: q_surf_in !sheet flow entering TC from uphill
    REAL, INTENT(OUT)                    :: q_surf_out
    REAL, INTENT(IN)                     :: q_sub_in			!internal representation of sublat_in
    REAL, INTENT(OUT)                    :: q_sub_out			!internal representation of sublat_out
    REAL, INTENT(OUT)                    :: qgw
    REAL, INTENT(OUT)                    :: deepqgw
    REAL, INTENT(OUT)                    :: hortf
    REAL, INTENT(OUT)                    :: tcaet
    REAL, INTENT(OUT)                    :: tclai
    REAL, INTENT(OUT)                    :: tcsoilet
    REAL, INTENT(OUT)                    :: tcintc
    REAL, INTENT(IN)                     :: prec_in     !input precipitation of current time step
    REAL, INTENT(IN)                     :: precday_in  !input daily precipitation
    REAL, INTENT(IN)                     :: prechall2_in(24) !input excerpt of preciph (rainfall) valid for 24 hours of the current subbas
    REAL, INTENT(IN)                     :: petday
    REAL, INTENT(IN)                     :: tcarea2
    REAL, INTENT(OUT)                    :: bal		!Till: water balance of current TC
    REAL, INTENT(IN)                     :: rootd_act(nveg)
    REAL, INTENT(IN)                     :: height_act(nveg)
    REAL, INTENT(IN)                     :: lai_act(nveg)
    REAL, INTENT(IN)                     :: alb_act(nveg)
    REAL, INTENT(IN)			         :: sed_in_tc(n_sed_class)
    REAL, INTENT(OUT)				     :: sed_out_tc(n_sed_class)
    REAL, INTENT(IN)                     :: q_rill_in !rill flow entering TC from uphill

    REAL                                 :: prec            !precipitation current time step potentially modified in snow module
    REAL                                 :: precday         !daily precipitation potentially modified in snow module
    REAL                                 :: prechall2(24)   !excerpt of preciph (24h of current subbas) potentially modified in snow module
    REAL                                 :: temperature     !temperature of current time step potentially modified in snow module
    REAL                                 :: radiation       !radiation of current time step potentially modified in snow module
    REAL                                 :: cloudFraction   !cloudiness fraction
    REAL                                 :: airPress        !Air pressure of current time step

    logical :: isnan
    !** Soil water model for terrain component (TC)

    !** Main variables:
    !**   thact:  actual water content of TC (mm)
    !**   thfree: water content above field capacity (mm)
    !**   q_surf_in:    surface runoff from upper TC (m**3)
    !**   q_surf_out:    surface runoff towards lower TC (m**3)
    !**   q_sub_in:    lateral subsurface flow from upper TC (m**3)
    !**   q_sub_out:    lateral subsurface flow towards lower TC (m**3)
    !**   qgw:    percolation below root zone (mm)
    !**   deepqgw:deep groundwater recharge (leaving profile)(mm)
    !**   hortf:  hortonian overland flow (m**3)
    !**   k_sat:saturated hydraulic conductivity (mm/day)
    !**   slope:  slope of terrain component (%)
    !**   saug:   suction at the wetting front (mm)
    !**   prec:   precipitation (mm)
    !**   tcarea2:   area of terrain component (km**2)
    !**   na:     soil pore volume to be filled by infiltrating water (Vol%/100)
    !**   ts:     time to surface saturation by infiltration (day)
    !**   inf:    cumulative inflitration during time-step (mm)




    REAL :: thact1,watbal
    REAL :: evapveg,evaps,evapt
    REAL :: rockintc,scaet, satdep
    REAL :: etpmax(maxsoil)
    REAL :: evapsfrac,rsfinday !,evapvegfrac
    REAL :: ERR,inf,infalt,infsatt,test,INPUT,def,facw,infh,dt_per_day !,infall
    REAL :: tempx,temp2,temp3,temp4, temp5,tempna,tempalt,templat,kftemp
    REAL :: remainlat(size(latred,2))
    REAL :: na(maxsoil,maxhori)
    REAL :: tsh(maxsoil)
    REAL :: zfrac(maxsoil)	!Till: fraction of empty pore space (saturation deficit) [Vol. fraction]
    REAL :: thfree(maxsoil,maxhori)
    REAL :: etpfrac,aktnfk(maxhori)
    REAL :: l2(maxhori),l2eff(maxhori)
    REAL :: percol(maxhori),infup,fillup,remain,tshup
    REAL :: conduns(maxhori),qsurf(maxsoil)
    REAL :: qmerk(maxsoil),qmerk2(maxsoil),qotemp
    REAL :: horitemp(maxhori),inputrem(maxsoil)
    REAL :: macro(maxsoil,maxhori),namic,infmac(maxhori)
    REAL :: suction(maxhori),vangen(maxhori)
    REAL :: resf3(maxhori),etpvert(maxhori)
    REAL :: frac_old,tempth,percolmac(maxhori)
    REAL :: intc_evap(maxsoil),intcred(maxsoil),maxintc,aetred(maxsoil)
    REAL :: hfrac(maxsoil),merkalluv(maxsoil)

    REAL :: allalluv


    INTEGER :: nbrrooth(maxsoil)
    INTEGER :: i,it,j,h,soilid,n_iter
    INTEGER :: horton,hsat(maxsoil),testi,hmerk,testi2,lath

    REAL :: q_surf								!surface runoff [mm H2O]
    REAL :: q_rill_out, q_surf_out2				!modified surface fluxes to pass to sedi_yield
    REAL :: r
    REAL :: exp_x !temporary variable to check argument for exp() for reasonable value range (within -15...15 to avoid Inf als result)

    ! REAL :: dsi, Aq, Q_li, l2_test1   !these variables are only needed for derivation of lateral subsurface flow equations


    INTEGER :: doshrink


    ! character(len=80):: errstring

    if (tcarea2 == 0.) then !for dummy basins, don't do calculations
        !set all output arguments to zero (otherwise, they may end up some as undefined value)
        thactroot = 0.
        q_surf_out = 0.
        q_sub_out = 0.
        qgw = 0.
        deepqgw = 0.
        hortf = 0.
        tcaet = 0.
        tclai = 0.
        tcsoilet = 0.
        tcintc = 0.
        bal = 0.
        sed_out_tc(:) = 0.
        return   
    end if    
    
    !** ------------------------------------------------------
    !** (0) Initialization
    !** Soil water content at end of last timestep (mm)
    !   thact is average soil moisture of TC (mm).
    !   watbal checks water balance of TC for each timestep (mm)

    !for debugging - remove
    !DO i=1,nbr_svc(tcid_instance2)
    !	DO h=1,nbrhori(id_soil_intern(i,tcid_instance2))
    !		IF (horithact(tcid_instance2,i,h) - thetas(id_soil_intern(i,tcid_instance2),h)* horiz_thickness(tcid_instance2,i,h)> 0.1 ) THEN	!Till: if water content exceeds saturation...
    !			write(*,*)"start: Oversaturated horizon in TC/svc/horizon ",tcid_instance2,i,h
    !			!call pause1
    !		END IF
    !	END DO
    !END DO
    !for debugging - remove end


    remainlat=0.
    IF (dohour) THEN
        dt_per_day=24.
        doshrink=0

    ELSE
        !Till: these are all output variables that are currently not used
        !  hortsc(day,:)=0.		! horton overland flow of each SVC
        !  sat_area_of_sc(day,:)=0.		! saturated area each SVC, relative to total TC area
        !  sat_xs_overland_flow_sc(day,:)=0.		! saturation excess overland flow of each SVC
        !  hort2sc(day,:)=0.		! horton overland flow 2 of each SVC
        !
        !  deepgwrsc(day,:)=0.	!Till: deep groundwater recharge in each SVC
        !  thsc(tc_counter2,day,:)=0.	! soil moisture of each SVC
        !  resistsc(day,:)=0.
        !  aetsc(day,:)=0.		! actual transpiration of each SVC (only plants)
        !  aettotsc(day,:)=0.	! total evapotranspiration of each SVC
        !  intetsc(day,:)=0.		! interception evaporation of each SVC
        !  soiletsc(day,:)=0.	! soil evaporation of each SVC
        !  nfksc(day,:)=0.		! actual available field capacity of each SVC

        dt_per_day=1.
        doshrink=0			!ii: shrinkage is currently disabled in the model
    END IF

    q_surf_out=0.
    q_sub_out=0.
    qsurf(:)=0.
    thact1=thact
    watbal=0.
    !** ------------------------------------------------------------------
    !** (1.1) Update of soil moisture due to lateral inflow of
    !   subsurface flow from other soil-vegetation components (SVCs) of
    !   this terrain component (TC)
    !   (subsurface flow produced in timestep before)

    !   difference between upslope and lowest terrain component:
    !   in lowest TC, first alluvial soils, if existent, are refilled by lateral
    !   subsurface inflow
    !
    ! Till, all comments "?":
    !		Input from
    !		- lateral subssurface flow of last timestep (latred)
    !		Affects / returns
    !		- saturated soil water content distribution (mm) (tcthetas)
    !		- saturated soil water content (VOL%)	(tcthetas)
    !		- saturated fraction of each SVC in each TC (-) (frac_sat)



    IF (dolatscsub) THEN										!Till: if lateral subsurface runoff is to be computed...
        templat=sum(latred(tcid_instance2,:))					!lateral subsurface runoff to be redistributed between SVCs (m**3)
    ELSE
        templat=0.
    END IF

    IF (templat > 0.) THEN					!Till: is there subsurface runoff to be redistributed?

        !return flow from rocky areas...
        !...onto other SVCs
        IF (dolatsc .AND. rocky(tcid_instance2)>0. ) THEN
            DO i=1,nbr_svc(tcid_instance2)
                qsurf(i)=qsurf(i)+templat*rocky(tcid_instance2)*frac_svc(i,tcid_instance2)
            END DO

            !fix 2
            tempx = templat*rocky(tcid_instance2)*(1-rocky(tcid_instance2)) !amount of subsurface flow to convert to return flow
            latred(tcid_instance2,:) = latred(tcid_instance2,:) * (1-tempx/templat) !Till: reduce pending lateral inflow
            templat = templat - tempx
        END IF

        !...onto rocky areas -> surface runoff
        watbal=watbal+templat/(tcarea2*1.e3)	!ii tcarea gleich am Anfang mit 1000 multiplizieren
        IF (nbr_svc(tcid_instance2) == 0) THEN	!Till: if whole TC consists of rock...
            tempx = 1.				!Till: subsurface inflow becomes returnflow and leaves the TC as surface runoff
        ELSE IF (.NOT. dolatsc) THEN				!Till: if no lateral flowdistribution...
            tempx = rocky(tcid_instance2)	!Till: surface runoff is increased according to fraction of rocky parts
        ELSE
            tempx = rocky(tcid_instance2)*rocky(tcid_instance2) !Till: surface runoff is increased according to fraction of rocky parts, but reduced due to lateral redistribution
        END IF
        q_surf_out = q_surf_out + tempx * templat
        !fix 1
        latred(tcid_instance2,:) = latred(tcid_instance2,:) * (1-tempx) !Till: reduce pending lateral inflow !debugadded - remove?!
        templat = templat * (1-tempx)

        IF (tc_counter2 == nbrterrain(i_lu)) THEN	!Till: treat the very lowest TC in a different way: treat alluvial soils
            merkalluv(:)=0. !fraction of SVC covered by alluvial soil (i.e. 0 for non-alluvials, otherwise fraction of SVC)
            ! check if alluvial soils exist and sum up their fractions !ii: this is a static property and could be computed only once
            merkalluv(1:nbr_svc(tcid_instance2)) = &
                alluvial_flag(id_soil_intern(1:nbr_svc(tcid_instance2), tcid_instance2)) *&
                frac_svc     (               1:nbr_svc(tcid_instance2), tcid_instance2)

            ! total fraction of alluvial soils in current TC
            allalluv=sum(merkalluv(1:nbr_svc(tcid_instance2)))

            ! if any alluvial soils occur in TC
            IF (allalluv > 0.) THEN
                !  maximum lateral inflow into this SVC
                DO i=1,nbr_svc(tcid_instance2)
                    IF (merkalluv(i) > 0.) THEN	!Till: if the current SVC has an alluvial soil, treat it
                        soilid=id_soil_intern(i,tcid_instance2)
                        remainlat(:)=0.
                        testi=0

                        !  maximum lateral inflow into this SVC
                        h=size(remainlat)
                        remainlat(1:h)=(latred(tcid_instance2, 1:h)* merkalluv(i)/allalluv)/ (tcarea2*1.e3*frac_svc(i,tcid_instance2))
                        !Till: all potential subsurface flow is distributed among alluvial SVCs [mm] (still "remains to be distributed")

                        !  distribute inflow among horizons
                        DO h=1,nbrhori(soilid)
                            lath=min(size(remainlat),INT(sum(horiz_thickness(tcid_instance2,i,1:h))/500.))	!Till: compute number of "exchange horizons" to be used from remainlat to be fed into the current profile
                            IF (lath > 0) THEN
                                temp2=sum(remainlat(1:lath))			!Till: max amount of water that is available to be routed into this horizon [mm]
                                !ii: better simply use an extra var "excessfromhigherhorizons" or add the excess to remainlat(lath+1)
                                !  update soil moisture of this horizon
                                !  maximum inflow may be limited by hydraulic conductivity
                                IF (temp2 > 0.) THEN
                                    !Till: compute maximum inflow according to "exchange surface"
                                    !Andreas: has to be in unit mm, correct calculation similar to lateral outflow calculation
                                    temp3=k_sat(soilid,h)/dt_per_day*(1.-coarse(soilid,h))* &           ! Q_li
                                    horiz_thickness(tcid_instance2,i,h)*slope(id_tc_type2)/100.  &
                                        /  slength(i_lu) / fracterrain(id_tc_type2)/1000.

                                    temp3=MIN(temp3,temp2)					!Till: actual amount of water going into current horizon
                                    horithact(tcid_instance2,i,h)=horithact(tcid_instance2,i,h)+temp3 !mm

                                    remainlat(1:lath)=0.
                                    IF (temp3 < temp2) THEN					!Till: if inflow is limited by conductivity, the rest is stored for further redistribution
                                        remainlat(lath)=temp2-temp3
                                    END IF
                                    !  if horizon is saturated, remaining lateral inflow may be added to
                                    !  lower horizon
                                    temp5=thetas(soilid,h)*horiz_thickness(tcid_instance2,i,h)	!Till: compute saturated water content [mm]

                                    IF (horithact(tcid_instance2,i,h) > temp5) THEN !Till: if water content exceeds saturation...
                                        remainlat(lath)=remainlat(lath)+ (horithact(tcid_instance2,i,h)- temp5)
                                        horithact(tcid_instance2,i,h)=temp5
                                    END IF
                                END IF
                            END IF
                            ! enddo for all horizons
                        END DO

                        ! check if remaining lateral subsurface inflow
                        ! may then be distributed also in higher horizons which are below river bed
                        ! repeat procedure for all deeper horizons
                        ! remaining latflow is river runoff
                        temp2=sum(remainlat(:))
						IF (temp2 > 0.) THEN
                            lath=maxhori !default: use all horizons
                            DO h=1,maxhori !reduce number of exchange horizons according to depth of riverbed !ii: use "nbrhori(soilid)" instead?
                                IF (sum(horiz_thickness(tcid_instance2,i,1:h)) > riverbed(i_lu)) THEN
                                    lath=h !Till: wrong, when total soil profile is shallower than river: then, lowest horizon is still treated (but shouldn't be). Not (yet) fixed for legacy reasons.
                                    exit
                                END IF
                            END DO

							DO h=lath,nbrhori(soilid)
                                temp2=sum(remainlat(:))
                                !  update soil moisture of this horizon
                                !  maximum inflow may be limited by hydraulic conductivity
                                IF (temp2 > 0.) THEN

									!Till: compute maximum inflow according to "exchange surface"
                                    !Andreas: has to be in unit mm, correct calculation similar to lateral outflow calculation
                                    temp3=k_sat(soilid,h)/dt_per_day*(1.-coarse(soilid,h))* &           ! Q_li
                                    horiz_thickness(tcid_instance2,i,h)*slope(id_tc_type2)/100.  &
                                        /  slength(i_lu) / fracterrain(id_tc_type2)/1000.

                                    temp3=MIN(temp3,temp2)

                                    horithact(tcid_instance2,i,h)=horithact(tcid_instance2,i,h)+temp3

                                    remainlat(:)=0.

									IF (temp3 < temp2) THEN
                                        remainlat(1)=temp2-temp3
                                    END IF

                                    !  if horizon is saturated, remaining lateral inflow may be added to
                                    !  lower horizon
                                    temp5=thetas(soilid,h)*horiz_thickness(tcid_instance2,i,h)	!Till: compute saturated water content [mm]

                                    IF (horithact(tcid_instance2,i,h) > temp5) THEN !Till: if water content exceeds saturation...
                                        remainlat(1)=remainlat(1)+ (horithact(tcid_instance2,i,h)- temp5)
                                        horithact(tcid_instance2,i,h)=temp5
                                    END IF

                                END IF

                                ! enddo for all horizons
                            END DO

                        END IF

                        !  if not all lateral inflow to this SVC can be distributed among the horizons,
                        !  the remaining inflow is assumed to become subsurface runoff from this
                        !  terrain component
                        !  (this is e.g. also the case if the current soil profile is very shallow)
                        temp2=sum(remainlat(:))
						IF (temp2 > 0.) THEN
                            q_sub_out=q_sub_out+temp2*tcarea2*1.e3*frac_svc(i,tcid_instance2)	!Till: excess water is temporarily assumed subsurface runoff from TC [m3]
                            watbal=watbal-temp2*frac_svc(i,tcid_instance2)						!Till: ATTENTION: this water is redistributed among non-alluvial SVCs later
                        END IF


                        !  Calculate new fraction of saturation of current SVC
                        !  frac_sat value is relative to TC, not to SVC !!
                        tempth=sum(horithact(tcid_instance2,i,:)) !sum up total water content in entire profile [mm]
                        IF (tempth <= tctheta_s(1,1,tcid_instance2,i)) THEN
                        ! no part of SVC is saturated
                            frac_sat(tcid_instance2,i)=0.
                        ELSE
                        ! soil component is partly or completely saturated
                            !ii: Till: the following branch seems to serve for simulationg distributions of saturation
                            ! is is, however, disabled by having testi=0. Why?
                            !same for more instances down below
                            tempalt=0.
                            testi=0

                            DO j=1,4
                                IF (testi == 0) THEN

                                    tempx=tempalt+(tctheta_s(j+1,1,tcid_instance2,i)-  &
                                        tctheta_s(j,1,tcid_instance2,i))* (tctheta_s(j+1,2,tcid_instance2,i)-  &
                                        tctheta_s(j,2,tcid_instance2,i))/2.+ (tctheta_s(j+1,1,tcid_instance2,i)-  &
                                        tctheta_s(j,1,tcid_instance2,i))* (tctheta_s(5,2,tcid_instance2,i)-  &
                                        (tctheta_s(j+1,2,tcid_instance2,i)- tctheta_s(1,2,tcid_instance2,i)))

                                    IF (tempth > tempalt .AND. tempth < tempx) THEN

                                        frac_sat(tcid_instance2,i)=((tctheta_s(j,2,tcid_instance2,i)-  &
                                            tctheta_s(1,2,tcid_instance2,i))+  &
                                            (tempth-tempalt)/(tempx-tempalt)*  &
                                            (tctheta_s(j+1,2,tcid_instance2,i)- tctheta_s(j,2,tcid_instance2,i)))*  &
                                            frac_svc(i,tcid_instance2)
                                        testi=1
                                    END IF
                                    tempalt=tempx
                                END IF
                            END DO
                        END IF ! end of update saturated fraction of SVC

                    END IF	! endif alluvial SVC

                END DO	! enddo over all SVCs


                if (q_sub_out > 0) then	!Till: if there are also non-alluvial soils AND not all subsurface inflow has been absorbed by the alluvials before
                    if ((allalluv <= 1.)) then
						!Till: distribute among exchange horizons with the same proportions as before the alluvial influence
						latred(tcid_instance2,:)=latred(tcid_instance2,:)/sum(latred(tcid_instance2,:))*q_sub_out
						watbal=templat/(tcarea2*1.e3)							!reset water balance
						q_sub_out=0.											!Till: no surface runoff yet - may result later after treatment of non-alluvial SVCs
                    else
						temp5=3 !should not occur
					end if
                else
                    latred(tcid_instance2,:)=0.								!Till: nothing more to do
                end if

            END IF ! endif alluvial soils occur in this most downslope TC
        ELSE
            allalluv=0. !set alluvial fraction to 0 for non-lowermost TCs !ii probably obsolete
        END IF ! end if this is the lowermost TC


        IF ((sum(latred(tcid_instance2,:)) > 0.) ) THEN !Till: if there is excess water not infiltrated into alluvial soils

            DO i=1,nbr_svc(tcid_instance2) !regular, non-alluvial SVCs

                soilid=id_soil_intern(i,tcid_instance2)
                remainlat(:)=0.
                testi=0

                IF ((tc_counter2 == nbrterrain(i_lu)) .AND. (alluvial_flag(soilid)==1)) THEN	!Till: check if this soil is an alluvial soil
                    cycle		!Till: alluvial soils already have been treated, so skip them
                END IF


                !  maximum lateral inflow into this SVC
                h=size(remainlat)
                remainlat(1:h)=(latred(tcid_instance2,1:h)*frac_svc(i,tcid_instance2)/(1.-allalluv)) / & !Till: distribute subsurface flow according to fraction among non-alluvial soils
                (tcarea2*1.e3*frac_svc(i,tcid_instance2))

                !  check each horizon of SVC if there is lateral inflow (from soil depth
                !  above this horizon)
                DO h=1,nbrhori(soilid)
                    lath=INT(sum(horiz_thickness(tcid_instance2,i,1:h))/500.)
                    IF (lath > 0) THEN
                        temp2=sum(remainlat(1:lath))

                        !  update soil moisture of this horizon
                        !  maximum inflow may be limited by hydraulic conductivity
                        IF (temp2 > 0.) THEN
                            !Till: compute maximum inflow according to "exchange surface"
                            IF (tc_counter2 == nbrterrain(i_lu)) THEN	!Till: treat the very lowest TC in a different way
                                !ii: this different treatment is most likely an error. I assume Andreas' bugfix (below) has only be applied to non-lowermost TCs (before merging the two code branches into the current form). Check this!
                                temp3=(k_sat(soilid,h)/dt_per_day*(1.-coarse(soilid,h))*  &
                                    horiz_thickness(tcid_instance2,i,h)*slope(id_tc_type2)/100)* &
                                    area(i_subbas2)*frac_lu(lu_counter2,i_subbas2) * frac_svc(i,tcid_instance2) /&
                                    slength(i_lu)		!subsurface inflow into horizon(s) admitted due to conductivity
                            ELSE
                                !Andreas: has to be in unit mm, correct calculation similar to lateral outflow calculation
                                temp3=k_sat(soilid,h)/dt_per_day*(1.-coarse(soilid,h))* &           ! Q_li
                                horiz_thickness(tcid_instance2,i,h)*slope(id_tc_type2)/100.  &
                                    /  slength(i_lu) / fracterrain(id_tc_type2)/1000.
                            END IF


                            temp3=MIN(temp3,temp2)
                            horithact(tcid_instance2,i,h)=horithact(tcid_instance2,i,h)+temp3

                            remainlat(1:lath)=0.
                            IF (temp3 < temp2) THEN
                                remainlat(lath)=temp2-temp3
                            END IF

                           !  if horizon is saturated, remaining lateral inflow may be added to
                           !  lower horizon
                            temp5=thetas(soilid,h)*horiz_thickness(tcid_instance2,i,h)	!Till: compute saturated water content [mm]

                            IF (horithact(tcid_instance2,i,h) > temp5) THEN !Till: if water content exceeds saturation...
                                remainlat(lath)=remainlat(lath)+ (horithact(tcid_instance2,i,h)- temp5)
                                horithact(tcid_instance2,i,h)=temp5
                            END IF

                    END IF
                END IF
                !  enddo horizons
            END DO

            !  if not all lateral inflow to this SVC can be distributed among the horizons,
            !  the remaining inflow is assumed to become subsurface runoff from this
            !  terrain component
            !  (this is e.g. also the case if the current soil profile is very shallow)
            temp2=sum(remainlat(:))
            IF (temp2 > 0.) THEN
                q_sub_out=q_sub_out+temp2*tcarea2*1.e3*frac_svc(i,tcid_instance2)
                watbal=watbal-temp2*frac_svc(i,tcid_instance2)
            END IF


            !  Calculate new fraction of saturation of current SVC
            !  frac_sat value is relative to TC, not to SVC !!

            tempth=sum(horithact(tcid_instance2,i,:))

            ! no part of SVC is saturated

            IF (tempth <= tctheta_s(1,1,tcid_instance2,i)) THEN
                frac_sat(tcid_instance2,i)=0.

                ! soil component is partly or completely saturated
            ELSE

                tempalt=0.
                testi=0

                DO j=1,4
                    IF (testi == 0) THEN

                        tempx=tempalt+(tctheta_s(j+1,1,tcid_instance2,i)-  &
                            tctheta_s(j,1,tcid_instance2,i))* (tctheta_s(j+1,2,tcid_instance2,i)-  &
                            tctheta_s(j,2,tcid_instance2,i))/2.+ (tctheta_s(j+1,1,tcid_instance2,i)-  &
                            tctheta_s(j,1,tcid_instance2,i))* (tctheta_s(5,2,tcid_instance2,i)-  &
                            (tctheta_s(j+1,2,tcid_instance2,i)- tctheta_s(1,2,tcid_instance2,i)))

                        IF (tempth > tempalt .AND. tempth < tempx) THEN

                            frac_sat(tcid_instance2,i)=((tctheta_s(j,2,tcid_instance2,i)-  &
                                tctheta_s(1,2,tcid_instance2,i))+ (tempth-tempalt)/(tempx-tempalt)*  &
                                (tctheta_s(j+1,2,tcid_instance2,i)- tctheta_s(j,2,tcid_instance2,i)))*  &
                                frac_svc(i,tcid_instance2)
                            testi=1
                        END IF
                        tempalt=tempx
                    END IF
                END DO
            END IF
            ! end of update saturated fraction of SVC

            ! end of loop for all non-alluvial SVCs
        END DO
    END IF      ! end if lateral inflow > 0 (after alluvials)

    END IF		! end if (lateral inflow templat > 0
    latred(tcid_instance2,:)=0.				!Till: all lateral subsurface flow has been treated




    !** ------------------------------------------------------------------
    !** (1.2) Update of soil moisture due to lateral inflow of
    !   subsurface flow from higher TC.
    !   If a part of TC is already saturated, lateral flow into this
    !   area gives directly surface runoff (return flow)
    !
    ! Till, all comments "?":
    !		Input from
    !		- lateral subsurface flow from higher TC (q_sub_in)
    !		Affects / returns
    !		- saturated soil water content distribution (mm) (tcthetas)
    !		- saturated soil water content (VOL%)	(tcthetas)
    !		- saturated fraction of each SVC in each TC (-) (frac_sat)
    !		- surface runoff (qsurf)

    IF (q_sub_in /= 0.) THEN				!Till: lateral inflow present

        watbal=watbal+q_sub_in/(tcarea2*1.e3)
        IF (nbr_svc(tcid_instance2) == 0) THEN	!Till: if whole TC consists of rock...
            q_surf_out=q_surf_out+q_sub_in
        ELSE IF (.NOT. dolatsc) THEN	!Till: if no lateral flowdistribution...
            q_surf_out=q_surf_out+q_sub_in*rocky(tcid_instance2)
        ELSE
            q_surf_out=q_surf_out+q_sub_in*rocky(tcid_instance2)*rocky(tcid_instance2)
        END IF

        DO i=1,nbr_svc(tcid_instance2)

            soilid=id_soil_intern(i,tcid_instance2)
            h=1
            remain=0.
            testi=0
            frac_old=frac_sat(tcid_instance2,i)

            !  inflow into SVCs by return flow of rocky areas
            IF (dolatsc) THEN
                qsurf(i)=qsurf(i)+q_sub_in*rocky(tcid_instance2)*frac_svc(i,tcid_instance2)
            END IF

            !  return flow from saturated areas (which are saturated at the
            !  beginning of this timestep)
            qotemp=q_sub_in*frac_sat(tcid_instance2,i)	!Till: potential return flow originating from current SVC
            IF (dolatsc) THEN
                q_surf_out=q_surf_out+qotemp*(frac_svc(i,tcid_instance2)+rocky(tcid_instance2))	!Till: the fraction of the current SVC and the rocky fraction directly lead to runoff
                DO j=1,nbr_svc(tcid_instance2)
                    IF (j /= i) THEN
                        qsurf(j)=qsurf(j)+qotemp*frac_svc(j,tcid_instance2)	!Till: return flow redistributed among other SVCs; leads to surface runoff, according to their degree saturation
                    END IF
                END DO
            ELSE
                q_surf_out=q_surf_out+qotemp
            END IF

            !  return flow from SVC if storage capacity is exceeded in this
            !  timestep by lateral inflow

            remain=(q_sub_in*frac_svc(i,tcid_instance2)-qotemp)/ (tcarea2*1.e3*frac_svc(i,tcid_instance2)) !Till: amount of subsurface flow that needs to be distributed in the current SVC (that has not been diverted to other SVCs already) [mm]


            tempx = sum(horithact(tcid_instance2,i,:)  - thetas(soilid,:) * horiz_thickness(tcid_instance2,i,:) ) + remain !Till: compute amount of excess water [mm]

            IF (tempx > 0.) THEN
                tempth=thsprof(tcid_instance2,i)				 !Till: set new water content of profile to saturation
                frac_sat(tcid_instance2,i)=frac_svc(i,tcid_instance2)

                IF (dolatsc) THEN
                    q_surf_out=q_surf_out+tempx*tcarea2*1.e3*frac_svc(i,tcid_instance2)*  &
                        (frac_svc(i,tcid_instance2)+rocky(tcid_instance2))
                    DO j=1,nbr_svc(tcid_instance2)
                        IF (j /= i) THEN
                            qsurf(j)=qsurf(j)+tempx*tcarea2*1.e3*frac_svc(i,tcid_instance2)*  &
                                frac_svc(j,tcid_instance2)
                        END IF
                    END DO
                ELSE
                    q_surf_out=q_surf_out+tempx*tcarea2*1.e3*frac_svc(i,tcid_instance2)
                END IF

                horithact(tcid_instance2,i, 1:nbrhori(soilid))=thetas(soilid,1:nbrhori(soilid))* horiz_thickness(tcid_instance2,i,1:nbrhori(soilid))	!Till: whole profile is saturated

            ELSE

                !  distribution of lateral inflow among horizons if profile
                !  is not completely saturated

                tempth = tempx + thsprof(tcid_instance2,i) !Till: new total water content in soil profile [mm]
                horitemp(1:maxhori)=horiths(tcid_instance2,i,1:maxhori) !Till: get relative distribution of maximum water storage among horizons of a profile

                DO WHILE (remain > 0.001)

                    !  distribute according to fraction of thetasat [mm] of each horizon
                    !  on total thetasat of profile [mm] (given by horiths)

                    horithact(tcid_instance2,i,1:nbrhori(soilid)) = horithact(tcid_instance2,i,1:nbrhori(soilid))+ remain*horitemp(1:nbrhori(soilid))

                    remain=0.

                    !  check if any horizon is more than saturated

                    DO h=1,nbrhori(soilid)
                        temp2=thetas(soilid,h)*horiz_thickness(tcid_instance2,i,h)	!max water content of horizon [mm]
                        IF (horithact(tcid_instance2,i,h) > temp2) THEN
                            remain=remain+horithact(tcid_instance2,i,h)- temp2
                            horithact(tcid_instance2,i,h)=temp2
                            horitemp(h)=0.
                        END IF
                    END DO

                    !  if there is remaining lateral inflow, calculate new horitemp,
                    !  excluding saturated horizons (horitemp = 0)
                    !  and repeat entire procedure until remain gets zero.
                    IF (remain > 0.001) THEN
                        tempx=sum(horitemp(:))
                        DO j=1,nbrhori(soilid)
                            horitemp(j)=horitemp(j)*(1./tempx)
                        END DO
                    END IF
                    !  end of loop remain > 0
                END DO

                !  end if-loop if SVC is saturated
            END IF


            !** calculate new fraction of saturation of current SVC
            !  frac_sat value is relative to TC, not to SVC !!

            ! no part of SVC is saturated

            IF (tempth <= tctheta_s(1,1,tcid_instance2,i)) THEN
                frac_sat(tcid_instance2,i)=0.

                ! soil component is partly or completely saturated
            ELSE

                tempalt=0.
                testi=0

                DO j=1,4
                    IF (testi == 0) THEN

                        tempx=tempalt+(tctheta_s(j+1,1,tcid_instance2,i)- tctheta_s(j,1,tcid_instance2,i))*  &
                            (tctheta_s(j+1,2,tcid_instance2,i)- tctheta_s(j,2,tcid_instance2,i))/2.+  &
                            (tctheta_s(j+1,1,tcid_instance2,i)- tctheta_s(j,1,tcid_instance2,i))*  &
                            (tctheta_s(5,2,tcid_instance2,i)- (tctheta_s(j+1,2,tcid_instance2,i)-  &
                            tctheta_s(1,2,tcid_instance2,i)))

                        IF (tempth > tempalt .AND. tempth < tempx) THEN

                            frac_sat(tcid_instance2,i)=((tctheta_s(j,2,tcid_instance2,i)-  &
                                tctheta_s(1,2,tcid_instance2,i))+ (tempth-tempalt)/(tempx-tempalt)*  &
                                (tctheta_s(j+1,2,tcid_instance2,i)- tctheta_s(j,2,tcid_instance2,i)))*  &
                                frac_svc(i,tcid_instance2)
                            testi=1
                        END IF
                        tempalt=tempx
                    END IF
                END DO
            END IF

            !  close loops for all SVCs and if lateral inflow > 0
        END DO
    END IF

    !only for output
    !DO i=1,nbr_svc(tcid_instance2)	!Till: enlarge saturated area fraction of each svc
    !  sat_area_of_sc(day,i)=sat_area_of_sc(day,i)+frac_sat(tcid_instance2,i)/dt_per_day
    !END DO
    !

    !thsc(tc_counter2,day,1:nbr_svc(tcid_instance2))=sum(horithact(tcid_instance2,1:nbr_svc(tcid_instance2),:))
    thact=sum(sum(horithact(tcid_instance2,1:nbr_svc(tcid_instance2),:), dim=2)*frac_svc(1:nbr_svc(tcid_instance2),tcid_instance2)	)




    !!for debugging - remove
    !DO i=1,nbr_svc(tcid_instance2)
    !	DO h=1,nbrhori(id_soil_intern(i,tcid_instance2))
    !		IF (horithact(tcid_instance2,i,h) - thetas(id_soil_intern(i,tcid_instance2),h)* horiz_thickness(tcid_instance2,i,h)> 0.1 ) THEN	!Till: if water content exceeds saturation...
    !			write(*,*)"after qsubin: Oversaturated horizon in TC/svc/horizon ",tcid_instance2,i,h
    !			call pause1
    !		END IF
    !	END DO
    !END DO
    !!for debugging - remove end



    prec          =   prec_in !these assignments are necessary even if snow module not active; in steps after snow module usage of prec
    prechall2     =   prechall2_in
    precday       =   precday_in

    !** -------------------------------------------------------------------------
    !SNOW MODULE
    !Physically-based simulations based on energy balance method of ECHSE (Eco-hydrological Simulation Environment)

    if(dosnow > 0) then

        temperature   =   temp(day,i_subbas2)
        radiation     =   rad( day,i_subbas2)
        airPress      =   1000.

        !Determination indices for optional arrays
        !When optional output not activated array only with dimensions (1,1,1) allocated
        !Check whether activated using logical factor read in from outfiles.dat
        !Function fLog() in snow_h.f90
        radiModIndices(:)    = fLog(f_snowTemp,  day, hh, tcid_instance2)
        temperaModIndices(:) = fLog(f_surfTemp,  day, hh, tcid_instance2)
        cloudFracIndices(:)  = fLog(f_cloudFrac, day, hh, tcid_instance2)
        snowTempIndices(:)   = fLog(f_snowTemp,  day, hh, tcid_instance2)
        surfTempIndices(:)   = fLog(f_surfTemp,  day, hh, tcid_instance2)
        liquFracIndices(:)   = fLog(f_liquFrac,  day, hh, tcid_instance2)
        fluxPrecIndices(:)   = fLog(f_fluxPrec,  day, hh, tcid_instance2)
        fluxSublIndices(:)   = fLog(f_fluxSubl,  day, hh, tcid_instance2)
        fluxFlowIndices(:)   = fLog(f_fluxFlow,  day, hh, tcid_instance2)
        fluxNetSIndices(:)   = fLog(f_fluxNetS,  day, hh, tcid_instance2)
        fluxNetLIndices(:)   = fLog(f_fluxNetL,  day, hh, tcid_instance2)
        fluxSoilIndices(:)   = fLog(f_fluxSoil,  day, hh, tcid_instance2)
        fluxSensIndices(:)   = fLog(f_fluxSens,  day, hh, tcid_instance2)
        stoiPrecIndices(:)   = fLog(f_stoiPrec,  day, hh, tcid_instance2)
        stoiSublIndices(:)   = fLog(f_stoiSubl,  day, hh, tcid_instance2)
        stoiFlowIndices(:)   = fLog(f_stoiFlow,  day, hh, tcid_instance2)
        rateAlbeIndices(:)   = fLog(f_rateAlbe,  day, hh, tcid_instance2)

       !Subroutine to modify meteo-drivers according to location
       !Preparation before feeding the snow model

       call snow_prepare_input(hh, day, i_subbas2, lu_counter2, tc_counter2, prec, temperature, radiation, &
                               cloudFraction, tempLaps, tempAmplitude, tempMaxOffset)

       !Collect modified radiation(aspect, slope) and temperature (elevation) signal
       radiMod(radiModIndices(1),radiModIndices(2),radiModIndices(3))             = radiation
       temperaMod(temperaModIndices(1),temperaModIndices(2),temperaModIndices(3)) = temperature

       !Wind currently not read from input file; assumed constant wind = 1; see climo.f90
       !Air pressure constant; for now set to 1000 hPa
       if(dohour) then
         call snow_compute(prec, temperature_tc, radiation_tc, airPress, rhum(day,i_subbas2), wind(day,i_subbas2), cloudFraction, &
                         snowEnergyCont(max(1,day-1),max(1,hh-1),tcid_instance2), snowWaterEquiv(max(1,day-1), max(1,hh-1),tcid_instance2), &
                         snowAlbedo(max(1,day-1),max(1,hh-1),tcid_instance2), snowEnergyCont(day, hh, tcid_instance2), snowWaterEquiv(day, hh, tcid_instance2), &
                         snowAlbedo(day, hh, tcid_instance2), snowCover(day, hh, tcid_instance2), &
                         snowTemp(snowTempIndices(1),snowTempIndices(2),snowTempIndices(3)), surfTemp(surfTempIndices(1),surfTempIndices(2),surfTempIndices(3)), &
                         liquFrac(liquFracIndices(1),liquFracIndices(2),liquFracIndices(3)), fluxPrec(fluxPrecIndices(1),fluxPrecIndices(2),fluxPrecIndices(3)), &
                         fluxSubl(fluxSublIndices(1),fluxSublIndices(2),fluxSublIndices(3)), fluxFlow(fluxFlowIndices(1),fluxFlowIndices(2),fluxFlowIndices(3)), &
                         fluxNetS(fluxNetSIndices(1),fluxNetSIndices(2),fluxNetSIndices(3)), fluxNetL(fluxNetLIndices(1),fluxNetLIndices(2),fluxNetLIndices(3)), &
                         fluxSoil(fluxSoilIndices(1),fluxSoilIndices(2),fluxSoilIndices(3)), fluxSens(fluxSensIndices(1),fluxSensIndices(2),fluxSensIndices(3)), &
                         stoiPrec(stoiPrecIndices(1),stoiPrecIndices(2),stoiPrecIndices(3)), stoiSubl(stoiSublIndices(1),stoiSublIndices(2),stoiSublIndices(3)), &
                         stoiFlow(stoiFlowIndices(1),stoiFlowIndices(2),stoiFlowIndices(3)), rateAlbe(rateAlbeIndices(1),rateAlbeIndices(2),rateAlbeIndices(3)), &
                         precipMod(day, hh, tcid_instance2), cloudFrac(cloudFracIndices(1),cloudFracIndices(2),cloudFracIndices(3)))

         prec           = precipMod(day, hh, tcid_instance2) !Further calculations with modified precipitation signal
         prechall2(hh)  = precipMod(day, hh, tcid_instance2) !Further calculations with modified precipitation signal

         if(hh ==24) then !to get into the next day; have start value
            snowEnergyCont(day+1, 1, tcid_instance2) = snowEnergyCont(day, 24, tcid_instance2)
            snowWaterEquiv(day+1, 1, tcid_instance2) = snowWaterEquiv(day, 24, tcid_instance2)
            snowAlbedo    (day+1, 1, tcid_instance2) = snowAlbedo    (day, 24, tcid_instance2)
         end if

       else !daily

         call snow_compute(prec, temperature_tc, radiation_tc, airPress, rhum(day,i_subbas2), wind(day,i_subbas2), cloudFraction, &
                         snowEnergyCont(max(1,day-1),max(1,hh-1),tcid_instance2), snowWaterEquiv(max(1,day-1), max(1,hh-1),tcid_instance2), &
                         snowAlbedo(max(1,day-1),max(1,hh-1),tcid_instance2), snowEnergyCont(day, hh, tcid_instance2), snowWaterEquiv(day, hh, tcid_instance2), &
                         snowAlbedo(day, hh, tcid_instance2), snowCover(day, hh, tcid_instance2), &
                         snowTemp(snowTempIndices(1),snowTempIndices(2),snowTempIndices(3)), surfTemp(surfTempIndices(1),surfTempIndices(2),surfTempIndices(3)), &
                         liquFrac(liquFracIndices(1),liquFracIndices(2),liquFracIndices(3)), fluxPrec(fluxPrecIndices(1),fluxPrecIndices(2),fluxPrecIndices(3)), &
                         fluxSubl(fluxSublIndices(1),fluxSublIndices(2),fluxSublIndices(3)), fluxFlow(fluxFlowIndices(1),fluxFlowIndices(2),fluxFlowIndices(3)), &
                         fluxNetS(fluxNetSIndices(1),fluxNetSIndices(2),fluxNetSIndices(3)), fluxNetL(fluxNetLIndices(1),fluxNetLIndices(2),fluxNetLIndices(3)), &
                         fluxSoil(fluxSoilIndices(1),fluxSoilIndices(2),fluxSoilIndices(3)), fluxSens(fluxSensIndices(1),fluxSensIndices(2),fluxSensIndices(3)), &
                         stoiPrec(stoiPrecIndices(1),stoiPrecIndices(2),stoiPrecIndices(3)), stoiSubl(stoiSublIndices(1),stoiSublIndices(2),stoiSublIndices(3)), &
                         stoiFlow(stoiFlowIndices(1),stoiFlowIndices(2),stoiFlowIndices(3)), rateAlbe(rateAlbeIndices(1),rateAlbeIndices(2),rateAlbeIndices(3)), &
                         precipMod(day, hh, tcid_instance2), cloudFrac(cloudFracIndices(1),cloudFracIndices(2),cloudFracIndices(3)))

          !Correct if SWE is 0
          if(snowWaterEquiv(day,hh,tcid_instance2) <=  0.)   then
             snowWaterEquiv(day,hh,tcid_instance2)  =  0.
             snowEnergyCont(day,hh,tcid_instance2)  =  0.
             snowAlbedo(day,hh,tcid_instance2)      =  albedoMax
          end if

          prec    = precipMod(day, hh, tcid_instance2) !Further calculations with modified precipitation signal
          precday = precipMod(day, hh, tcid_instance2) !Further calculations with modified precipitation signal

       end if

    end if


    !** -------------------------------------------------------------------------
    !** (11) INTERCEPTION
    !        approach similar to WASIM-ETH by Schulla (p.50f)
    !        one-bucket model
    !        storage capacity per unit LAI (intcf) is set to 0.25 mm
    !        maximum storage capacity of canopy: maxintc
    !        ! calculation approach depends on wether LAI refers to entire area
    !          or only to vegetation cover area !
    !        intc_evap: evaporation from intc storage (intercept)
    !        (update of TC water balance is done at end of evapotranspiration
    !         for non-rocky areas)
    !         reduction of plant and soil evapotranspiration as function
    !         of fraction of intercep evap to potential evapotranspiration etp_max
    !         etp_max is calculated in separate subroutine for rs=0, and rs,soil=0
    !
    ! Till, all comments "?":
    !		Input from
    !		- precipitation [mm] (prec, precday)
    !		- current interception storage of each SVC in each TC (mm) (intercept)
    !		Affects / returns
    !		- current interception storage of each SVC in each TC (mm) (intercept)
    !		- reduction term for potential soil evaporation and transpiration (aetred)
    !		- Water balance for TC (watbal)

    tcintc=0.
    rockintc=0.
    tclai=0.
    intc_evap(:)=0.
    intcred(:)=0.
    aetred(:)=0.
    etpmax(:)=0.

    !  initial interception storage volume
    tempx=sum(intercept(tcid_instance2,1:nbr_svc(tcid_instance2))*frac_svc(1:nbr_svc(tcid_instance2),tcid_instance2))	!Till: sum up amount of intercepted water in the entire TC [mm]


    ! only rock surface
    IF (nbr_svc(tcid_instance2) == 0) THEN	!Till: if whole TC consists of rock...
        IF (precday > 0.) THEN
            tcintc=MIN(intcfc(i_ce),precday)	!Till: fill interception storage according to rainfall and capacity, intercepted water evaporates completely
            !Wird in der Stundenvariante dann jede Stunde der komplette interzeptionsspeicher geleert?
            rockintc=tcintc
            tclai=0.
            watbal=watbal-tcintc
        END IF

    ELSE

        DO j=1,nbr_svc(tcid_instance2)
            ! determine maximum (potential) evapotranspiration of SVC
            CALL etp_max(i_subbas2,id_veg_intern(j,tcid_instance2),etpmax(j), height_act,lai_act,alb_act,day,hh)


            ! mean LAI of TC
            tclai=tclai+lai_act(id_veg_intern(j,tcid_instance2))*frac_svc(j,tcid_instance2)

            IF (dointc == 0) THEN	! interception, simple bucket model

                ! hourly version of interception calculated for aggregated daily rainfall
                ! distributed among hours according to rainfall fraction of each
                ! hour on total daily rainfall (Till: because interception capacity is usually related to daily rainfall)
                ! stored for later hours of day in intcept_mem
                ! similarly, reduction of plant aet is stored in aet_red_mem
                IF (dohour) THEN	!Till: hourly version
                    IF (hh == 1) THEN
                        intcept_mem(tc_counter2,j,:)=0.	!ii Till: not needed, overwritten below
                        aet_red_mem(tc_counter2,j)=0.
                        IF (precday > 0.) THEN
                            maxintc=lai_act(id_veg_intern(j,tcid_instance2))*intcfc(i_ce)		!Till: compute maximum possible interception (daily value)
                            intcred(j)=MIN(precday,maxintc-intercept(tcid_instance2,j))	!Till: compute amount that is intercepted in this timestep (daily value)
                            intercept(tcid_instance2,j)=intercept(tcid_instance2,j)+ intcred(j)	!Till: update information on status of interception storage
                        END IF

                        IF (intercept(tcid_instance2,j) > 0.) THEN						!Till: if there is intercepted water...
                            !             intc_evap(j)=min(intercept(tcid_instance2,j),petday)
                            intc_evap(j)=MIN(intercept(tcid_instance2,j),etpmax(j))		!Till: compute interception loss by evaporation
                            intercept(tcid_instance2,j)=intercept(tcid_instance2,j)- intc_evap(j)	!Till: reduce interception storage by evaporated amount
                            aetred(j)=intc_evap(j)/max(etpmax(j),0.00001)						!Till: compute evaporation reduction term
                            aet_red_mem(tc_counter2,j)=aetred(j)							!Till: remember this reduction term on hourly basis
                            tcintc=tcintc+intc_evap(j)*frac_svc(j,tcid_instance2)			!Till: sum up all interception losses for entire TC

                            ! distribute among hourly timesteps								!ii das soll dann auch mal mit anderen Zeitschritten funktionieren
                            if (precday==0.) then						!no precip this day - distribute evaporation of leftover interception equally
                                intcept_mem=intc_evap(j)/24.
                            else
                                intcept_mem(tc_counter2,j, 1:24)=intc_evap(j)* prechall2(1:24)/precday	!Till: distribute according to rainfall amount that fell
                            end if

                            intcred(j)=intcept_mem(tc_counter2,j,hh)
                            intc_evap(j)=intcept_mem(tc_counter2,j,hh)
                        END IF
                    ELSE IF (hh > 1) THEN
                        intc_evap(j)=intcept_mem(tc_counter2,j,hh)			!Till: evaporation from interception storage (interception loss)
                        intcred(j)=intcept_mem(tc_counter2,j,hh)				!Till: Precip reduction due to interception in the current SVC and timestep
                        aetred(j)=aet_red_mem(tc_counter2,j)
                    END IF
                ELSE	! if not (dohour) - daily version
                    IF (prec > 0.) THEN
                        !ii: folgende 13 Zeilen scheinen gleich in beiden Schleifen bis auf prec/precday, auslagern
                        maxintc=lai_act(id_veg_intern(j,tcid_instance2))*intcfc(i_ce)			!Till: compute maximum possible interception
                        ! if only partial vegetation cover
                        !            maxintc=lai_act(id_veg_intern(j,tcid_instance2))*intcf*
                        !     .            (1.-0.7**lai_act(id_veg_intern(j,tcid_instance2)))+
                        !     .            intcf*(0.7**lai_act(id_veg_intern(j,tcid_instance2)))
                        intcred(j)=MIN(prec,maxintc-intercept(tcid_instance2,j))				!Till: compute amount that is intercepted in this timestep
                        intercept(tcid_instance2,j)=intercept(tcid_instance2,j)+ intcred(j)	!Till: update information on status of interception storage
                    END IF
                    IF (intercept(tcid_instance2,j) > 0.) THEN								!Till: if there is intercepted water...
                        !             intc_evap(j)=min(intercept(tcid_instance2,j),petday)
                        intc_evap(j)=MIN(intercept(tcid_instance2,j),etpmax(j))				!Till: compute interception loss by evaporation
                        intercept(tcid_instance2,j)=intercept(tcid_instance2,j)- intc_evap(j)	!Till: reduce interception storage by evaporated amount
                        aetred(j)=intc_evap(j)/max(etpmax(j),0.00001)										!Till: compute evaporation reduction term
                        tcintc=tcintc+intc_evap(j)*frac_svc(j,tcid_instance2)					!Till: sum up all interception losses for entire TC
                    END IF
                END IF

            ELSE IF (dointc == 1) THEN
                ! interception, modified bucket model
                ! (not yet updated to etp_max evaporation)
                maxintc=lai_act(id_veg_intern(j,tcid_instance2))*intcfc(i_ce)
                IF (prec > 0.) THEN
                    intcred(j)=maxintc*(1.-EXP(-0.8*prec))- intercept(tcid_instance2,j)
                    !            intcred(j)=min(prec,maxintc-intercept(tcid_instance2,j))
                    intercept(tcid_instance2,j)=intercept(tcid_instance2,j)+ intcred(j)
                    !        IF (intercept(tcid_instance2,j) < 0.) THEN
                    !          WRITE(*,*) intercept(tcid_instance2,j)
                    !        END IF
                END IF
                IF (intercept(tcid_instance2,j) > 0.) THEN
                    IF (dohour) THEN
                        IF (hh <= 6 .OR. hh >= 19) THEN
                            intc_evap(j)=0.
                        ELSE
                            intc_evap(j)=MIN(intercept(tcid_instance2,j), petday/(dt_per_day/2.)*  &
                                ((intercept(tcid_instance2,j)/maxintc)**(2./3.)))
                        END IF
                    ELSE
                        intc_evap(j)=MIN(intercept(tcid_instance2,j),petday*  &
                            ((intercept(tcid_instance2,j)/maxintc)**(2./3.)))
                    END IF
                    intercept(tcid_instance2,j)=intercept(tcid_instance2,j)- intc_evap(j)
                    tcintc=tcintc+intc_evap(j)*frac_svc(j,tcid_instance2)
                END IF
            END IF	!Till: end (modified bucket)

        END DO	!Till: end (thru all SVCs)

        IF (rocky(tcid_instance2) > 0.) THEN	!Till: treat rocky part of this TC
            IF (prec > 0.) THEN
                tcintc=tcintc+MIN(intcfc(i_ce),prec)*rocky(tcid_instance2)
                rockintc=MIN(intcfc(i_ce),prec)*rocky(tcid_instance2)
                watbal=watbal-MIN(intcfc(i_ce),prec)*rocky(tcid_instance2)
            END IF
        END IF

    END IF	!Till: end (if whole TC consists not of rock)

    !  final interception storage volume

    temp3=sum(intercept(tcid_instance2,1:nbr_svc(tcid_instance2))*frac_svc(1:nbr_svc(tcid_instance2),tcid_instance2)) !Till: sum up how much water is in the interception storage of the TC [mm]


    watbal=watbal-(temp3-tempx)							!Till: modify the water balance by the difference between storage before and after calculation [mm]



    !** -------------------------------------------------------------------------
    !** (2)   SURFACE RUNOFF
    !   (2.0) on rock outcrops:            overland flow
    !   (2.1) on saturated part of TC:     saturation excess overland flow occurs
    !   (2.2) on non-saturated part of TC: infiltration excess overland flow may occur


    !** Total TC input (mm)
    !** Input by precipitation and by surface runoff from higher TC (mm)
    !
    ! Till, all comments "?":
    !		Input from
    !		- precipitation [mm] (prec)
    !		- overland flow from higher TC (q_surf_in)
    !		Affects / returns
    !		-
    !		-
    !		-


    INPUT=prec+q_surf_in/(tcarea2*1.e3)		!Till: precipitation and overland flow from higher TCs [mm]
    watbal=watbal+INPUT	![mm]


    !** .....................................................................
    !**  (2.0) Runoff from impermeable surfaces
    !          input reduced by evaporation from surface storage

    IF (rocky(tcid_instance2) > 0.) THEN
        IF (INPUT-MIN(prec,intcfc(i_ce)) > 0.) THEN
            IF (nbr_svc(tcid_instance2) == 0) THEN						!Till: if whole TC consists of rock...
                q_surf_out=q_surf_out+(INPUT-MIN(prec,intcfc(i_ce)))*tcarea2*1.e3	!Till: ...everything except some interception runs off
            ELSE IF (.NOT. dolatsc) THEN						!Till: if no lateral flowdistribution...
                q_surf_out=q_surf_out+(INPUT-MIN(prec,intcfc(i_ce)))*tcarea2*1.e3*rocky(tcid_instance2)!...the rocky fraction gets exactly according to its share
            ELSE
                q_surf_out=q_surf_out+(INPUT-MIN(prec,intcfc(i_ce)))*tcarea2*1.e3*rocky(tcid_instance2)*rocky(tcid_instance2)
            END IF
            IF (dolatsc) THEN									!ii: use intcred(i) instead of computing interception anew?
                DO j=1,nbr_svc(tcid_instance2)
                    qsurf(j)=qsurf(j)+(INPUT-MIN(prec,intcfc(i_ce)))*tcarea2*1.e3*  &		!Till: to each SVC the runoff produced by the rocky fraction is assigned according to its fraction
                    rocky(tcid_instance2)*frac_svc(j,tcid_instance2)
                END DO
            END IF
        END IF
    END IF


    !** .....................................................................
    !** (2.1) Saturation excess runoff
    !         Input from
    !         - precipitation exceeding retention in interception storage (prec-intcred)
    !         - surface runoff from higher TC (q_surf_in)
    !         - surface runoff from rocky area of this TC (qsurf)
    !         - return flow from other SVCs of this TC (qsurf)
    !Till:
    !		Affects / returns
    !		- overland flow to lower TC / river (q_surf_out)
    !		- saturation excess overland flow of each SVC (sat_xs_overland_flow_sc)
    !		- cumulated "excess" water that is considered for surface runoff and infiltration (inputrem)
    !		- ...


    inputrem(:)=0.
    DO i=1,nbr_svc(tcid_instance2) !Till: for all SVC in the current TC, compute saturation excess runoff

        INPUT=prec-intcred(i)+ q_surf_in/(tcarea2*1.e3)+qsurf(i)/(tcarea2*1.e3*frac_svc(i,tcid_instance2))
        !Till: precipitation-intercepted+surface flow from above+surface(return) flow from within (mm)
        IF (INPUT > 0.) THEN

            IF (frac_sat(tcid_instance2,i) > 0.) THEN	!Till: if there is a saturated area in this SVC...

                IF (dolatsc) THEN					!Till: if lateral redistribution is considered...
                    q_surf_out=q_surf_out+INPUT*tcarea2*1.e3*frac_sat(tcid_instance2,i)*(frac_svc(i,tcid_instance2)+rocky(tcid_instance2))
                    inputrem(i)=inputrem(i)+INPUT*tcarea2*1.e3*  (frac_svc(i,tcid_instance2)-frac_sat(tcid_instance2,i))
                    DO j=1,nbr_svc(tcid_instance2)
                        IF (j /= i) THEN
                            inputrem(j)=inputrem(j)+INPUT*tcarea2*1.e3*frac_sat(tcid_instance2,i)*frac_svc(j,tcid_instance2)
                        END IF
                    END DO
                ELSE								!no lateral redistribution...
                    q_surf_out=q_surf_out+INPUT*tcarea2*1.e3*frac_sat(tcid_instance2,i)	!Till: everything that falls onto saturated area becomes runoff (m**3)
                    inputrem(i)=inputrem(i)+INPUT*tcarea2*1.e3*(frac_svc(i,tcid_instance2)-frac_sat(tcid_instance2,i)) !this water remains, doesn't drain instantly (mm*3)
                END IF

            ELSE								!Till: no saturated area yet....
                inputrem(i)=inputrem(i)+INPUT*tcarea2*1.e3*frac_svc(i,tcid_instance2)	!Till: ...this is the remaining water for this SVC (m**3)
            END IF
            !only for output
            !sat_xs_overland_flow_sc(day,i)=sat_xs_overland_flow_sc(day,i)+ INPUT*(frac_sat(tcid_instance2,i)/frac_svc(i,tcid_instance2))	!Till: sum up saturation excess overland flow of this TC
        END IF	!Till: end of INPUT>0
    END DO


    !**..............................................................
    !** (2.2a) Shrinkages in clay horizons and macroporosity
    IF (doshrink == 1) THEN
        DO i=1,nbr_svc(tcid_instance2)
            soilid=id_soil_intern(i,tcid_instance2)		!Till: get internal soil-id of SVC that is to be treated
            macro(i,:)=0.
            DO h=1,nbrhori(soilid)
                IF (shrink(soilid,h) == 1.0) THEN

                    tempx=((1./(1.+(3100.0/bubble(soilid,h))**  &
                        (poresz(soilid,h)+1.)))** porem(soilid,h))*  &
                        (thetas(soilid,h)/(1.-coarse(soilid,h))-  &
                        thetar(soilid,h)/(1.-coarse(soilid,h)))  &
                        +thetar(soilid,h)/(1.-coarse(soilid,h))
                    tempx=tempx*horiz_thickness(tcid_instance2,i,h)*(1.-coarse(soilid,h))
                    IF (horithact(tcid_instance2,i,h)< tempx) THEN
                        macro(i,h)=((tempx-horithact(tcid_instance2,i,h))/  &
                            (tempx-thetar(soilid,h)))*0.5* thetas(soilid,h)*  &
                            horiz_thickness(tcid_instance2,i,h)
                    END IF
                ELSE IF (shrink(soilid,h) == 0.0) THEN
                    macro(i,h)=0.01*horiz_thickness(tcid_instance2,i,h)
                    !           else
                    !            macro(i,h)=0.0025*horiz_thickness(tcid_instance2,i,h)
                END IF
            END DO
        END DO
    END IF

    !!for debugging - remove
    !DO k=1,nbr_svc(tcid_instance2)
    !	DO h=1,nbrhori(id_soil_intern(k,tcid_instance2))
    !		iF (horithact(tcid_instance2,k,h) - thetas(id_soil_intern(k,tcid_instance2),h)* horiz_thickness(tcid_instance2,k,h)> 0.1 ) THEN	!Till: if water content exceeds saturation...
    !			write(*,*)"vor loop(",n_iter,"): Oversaturated horizon in TC/svc/horizon ",tcid_instance2,k,h
    !			call pause1
    !		END iF
    !	END DO
    !END DO
    !!for debugging - remove end



    !** .....................................................................
    !** (2.2) Inflitration excess overland flow
    !         Green-Ampt-approach according to WASIM-ETH (J.Schulla / Peschke)
    !         extended by a multi-layer approach

    !         Input from
    !         - precipitation exceeding retention in interception storage
    !         - surface runoff from higher TC
    !         - surface runoff from rocky area of this TC
    !         - return flow from other SVCs of this TC
    !         - saturation excess runoff from SVCs of this TC
    !         -> all included in inputrem(1:nbr_svc)
    !Till:
    !		Affects / returns
    !		- overland flow to lower TC / river (q_surf_out)
    !		- (average) actual water content of TC (thact)
    !		- soil moisture in every horizon of each SVC in each TC in each LU in each subbasin (horithact)
    !		- ...



    qmerk(:)=0. !Till: for remembering infiltration excess between iterations
    qmerk2(:)=0.
    hortf=0.
    macro=0.	!Doris: inserted 10.12.2010 macro has to be initialised, otherwise namic is not calculated and
    !		the refillable porosity is not calcultated correctly (na(i,h) and namic)

    DO n_iter=1,2


        !!for debugging - remove
        !DO k=1,nbr_svc(tcid_instance2)
        !	DO h=1,nbrhori(id_soil_intern(k,tcid_instance2))
        !		IF (horithact(tcid_instance2,k,h) - thetas(id_soil_intern(k,tcid_instance2),h)* horiz_thickness(tcid_instance2,k,h)> 0.1 ) THEN	!Till: if water content exceeds saturation...
        !			write(*,*)"in loop(",n_iter,"): Oversaturated horizon in TC/svc/horizon ",tcid_instance2,k,h
        !			call pause1
        !		END iF
        !	END DO
        !END DO
        !!for debugging - remove end



        !*loop 1:.....................................................................
        !   (2.2a)
        !   if lateral flow between SVCs in this terrain component should be
        !   taken into account, in a first loop the expected infiltration
        !   excess overland flow of each SVC (according to the input of inputrem)
        !   is assessed and its contribution to input to other SVCs determined.
        !   In a second loop, this combined input, which may be larger than inputrem,
        !   is used to calculate real infiltration and soil moisture change of each SVC

        !   flow input (in m**3)(inputrem) is distributed over total SVC area
        !   (not exluding saturated part) to give input intensity (mm)

        !*loop 2:.....................................................................
        !** (2.2b) Infiltration excess overland flow, final calculation (second iteration)
        !         taking into account lateral inflow
        !         due to infiltration excess from other SVCs (qmerk)



        IF ((n_iter==1) .AND. (.NOT. dolatsc)) CYCLE !Till: no first iteration if disabled


        na(:,:)=0.
        DO i=1,nbr_svc(tcid_instance2)
            !for debugging - remove
            !	DO k=1,nbr_svc(tcid_instance2)
            !		DO h=1,nbrhori(id_soil_intern(k,tcid_instance2))
            !			IF (horithact(tcid_instance2,k,h) - thetas(id_soil_intern(k,tcid_instance2),h)* horiz_thickness(tcid_instance2,k,h)> 0.1 ) THEN	!Till: if water content exceeds saturation...
            !				write(*,*)"in svcloop(",n_iter,"): Oversaturated horizon in TC/svc/horizon ",tcid_instance2,k,h
            !				call pause1
            !			END iF
            !		END DO
            !	END DO
            !	!for debugging - remove end


            soilid=id_soil_intern(i,tcid_instance2)
            infh=0.
            zfrac(i)=0.

            hsat(i)=0
            infmac(:)=0.
            tsh(i)=0.
            kftemp=0.

            ! inputrem is in m**3, flow onto not-saturated areas of SVC
            IF (n_iter/=1) THEN !Till: second iteration
                inputrem(i)=inputrem(i)+qmerk(i)
            END IF
            INPUT=inputrem(i)/(tcarea2*1.e3*frac_svc(i,tcid_instance2))	!Till: converting remaining water input to water height [mm]


            !for debugging - remove
            !	DO k=1,nbr_svc(tcid_instance2)
            !		DO h=1,nbrhori(id_soil_intern(k,tcid_instance2))
            !			IF (horithact(tcid_instance2,k,h) - thetas(id_soil_intern(k,tcid_instance2),h)* horiz_thickness(tcid_instance2,k,h)> 0.1 ) THEN	!Till: if water content exceeds saturation...
            !				write(*,*)"input>0(",n_iter,"): Oversaturated horizon in TC/svc/horizon ",tcid_instance2,k,h
            !				call pause1
            !			END iF
            !		END DO
            !	END DO
            !	!for debugging - remove end
            !


            IF (INPUT > 0.) THEN

                !   soil component is saturated
                !   if dolatsc (lateral redistribution) there might be input from other SVCs
                IF (frac_sat(tcid_instance2,i) == frac_svc(i,tcid_instance2)) THEN
                    zfrac(i)=0.
                    IF (n_iter==1) THEN !Till first iteration
                        DO j=1,nbr_svc(tcid_instance2)
                            IF (j /= i) THEN
                                temp3=INPUT*tcarea2*1.e3*frac_svc(i,tcid_instance2)*frac_svc(j,tcid_instance2)
                                qmerk(j)=qmerk(j)+temp3
                                qmerk2(i)=qmerk2(i)+temp3
                            END IF
                        END DO
                    ELSE !Till: second iteration
                        IF (dolatsc) THEN
                            q_surf_out=q_surf_out+INPUT*tcarea2*1.e3*frac_svc(i,tcid_instance2)-qmerk2(i)
                            IF (dolattc .AND. tc_counter2 == nbrterrain(i_lu)) THEN
                                !c              hortf=hortf+(input*tcarea2*1.E3*frac_svc(i,tcid_instance2)-qmerk2(i))/
                                !c     .                    (area*1.E3)
                            END IF
                            !only for output
                            !        hortsc(day,i)=hortsc(day,i)+ (INPUT*tcarea2*1.e3*frac_svc(i,tcid_instance2)-qmerk2(i))/ (tcarea2*1.e3*frac_svc(i,tcid_instance2))

                            !   (else should not occur ??)
                        ELSE
                            q_surf_out=q_surf_out+INPUT*tcarea2*1.e3*frac_svc(i,tcid_instance2)
                            IF (dolattc .AND. tc_counter2 == nbrterrain(i_lu)) THEN
                                !c              hortf=hortf+input*frac_svc(i,tcid_instance2)
                            END IF
                            !only for output
                            !        hortsc(day,i)=hortsc(day,i)+INPUT
                        END IF
                    END IF !Till: branching between iterations


                    !   soil component is partly or not saturated
                ELSE
                    zfrac(i)=frac_svc(i,tcid_instance2)-frac_sat(tcid_instance2,i) !Till: fraction of non-saturation of current SVC

                    DO h=1,nbrhori(soilid)
                        na(i,h)=thetas(soilid,h)- horithact(tcid_instance2,i,h)/horiz_thickness(tcid_instance2,i,h) !Till: compute refillable porosity for all horizons [-]
                        !   horizon is saturated
                        IF (na(i,h) <= 0.0001) THEN
                            na(i,h)=0.0
                            hsat(i)=h				!Till: determine number of topmost saturated horizon (?)
                            exit
                        END IF
                    END DO
                END IF	!if (frac_sat(tcid_instance2,i) == frac_svc(i,tcid_instance2))


                !** if rainfall intensity plus surface runoff from higher TC
                !   exceeds hydr. conductivity, saturation may occur.
                !   tests if infiltration excess occurs for each SVC (index i)
                !   and horizon



                !for debugging - remove
                !	DO k=1,nbr_svc(tcid_instance2)
                !		DO h=1,nbrhori(id_soil_intern(k,tcid_instance2))
                !			IF (horithact(tcid_instance2,k,h) - thetas(id_soil_intern(k,tcid_instance2),h)* horiz_thickness(tcid_instance2,k,h)> 0.1 ) THEN	!Till: if water content exceeds saturation...
                !				write(*,*)"zfrac(i)>0(",n_iter,"): Oversaturated horizon in TC/svc/horizon ",tcid_instance2,k,h
                !				call pause1
                !			END iF
                !		END DO
                !	END DO
                !	!for debugging - remove end



                !   test only if soil is not saturated

                IF (zfrac(i) /= 0.) THEN
                    horton=0
                    fillup=0.
                    tshup=0.

                    !   repeat tests for all horizons until one with inf excess occurs
                    h=1
                    kftemp=k_sat(soilid,1)/dt_per_day*(1.-coarse(soilid,1))/  &
                        (kfkorrc(i_ce)*kfkorr_day)					!Till: compute kf in topmost horizon of current SVC ii: constant, compute only once for all soils
                    DO WHILE (horton == 0 .AND. h <= nbrhori(soilid)  &		!Till: do test for all unsaturated horizons
                        .AND. h /= hsat(i))

                        !   do test only if minimum saturated conductivity of all horizons
                        !   down to current horizon is is below input rate
                        !   if not uppermost horizon, calculation only if input exceeds
                        !   refillable porosity of higher horizons

                        IF (h > 1) THEN
                            fillup=fillup+na(i,h-1)*horiz_thickness(tcid_instance2,i,h-1)		!Till: reduction of input due to refilling of the horizon above [mm]
                            tshup=tshup+(na(i,h-1)*horiz_thickness(tcid_instance2,i,h-1))/INPUT !Till: relative reduction of input ? [-]
                        END IF


                        IF (INPUT > fillup .AND. tshup < 1.) THEN

                            kftemp=MIN(kftemp,k_sat(soilid,h)/dt_per_day*(1.-coarse(soilid,h))/  &
                                (kfkorrc(i_ce)*kfkorr_day))					!Till: compute kf in current horizon of current SVC ii: constant, compute only once for all soils
                            IF (INPUT > kftemp) THEN						!Till: why compared to INPUT, not INPUT-fillup?

                                !   effective refillable porosity in micropores is smaller then na
                                !   calculated above when taking into account shrinkages and makropores

                                IF (macro(i,h) > 0. .AND. shrink(soilid,h) == 1.0) THEN
                                    namic=thetas(soilid,h)-(macro(i,h)+  &
                                        horithact(tcid_instance2,i,h))/horiz_thickness(tcid_instance2,i,h)
                                ELSE IF (macro(i,h) == 0.) THEN
                                    namic=na(i,h)
                                END IF

                                !  depth of saturation front at time of surface saturation taking
                                !  wetting front suction of current horizon
                                satdep=saug(soilid,h)/(INPUT/kftemp-1.0)

                                !  if saturation front depth is larger than soil horizon depth:
                                !  test next horizon
                                !  else test for infiltration excess in this horizon
                                !  and infiltration routine ends with this horizon
                                IF (satdep <= sum(horiz_thickness(tcid_instance2,i,1:h))) THEN

                                    !  time of saturation
                                    tsh(i)=(namic*satdep)/INPUT

                                    !   if tsh is smaller than time step -> inf excess saturation of horizon
                                    IF (tsh(i) <= 1.0) THEN

                                        !   calculate cumulative infiltration (infh) for this time step
                                        !   iteration for current horizon
                                        !   total amount of infiltrated water is sum of
                                        !     - infiltration into horizon until tsh (before surface saturation)
                                        !     - infiltration into horizon after tsh (after  surface saturation)

                                        horton=h
                                        infsatt=tsh(i)*INPUT
                                        temp3=1.0-tsh(i)
                                        temp4=namic*saug(soilid,h)
                                        it=0
                                        ERR=1.0
                                        infalt=infsatt+temp3*INPUT/2.0
                                        DO WHILE ((it < 12) .AND. (ERR > 0.05))
                                            infh=kftemp*temp3 + temp4* LOG((infalt+temp4)/  &
                                                (infsatt+temp4))+infsatt
                                            ERR=ABS(infh-infalt)
                                            infalt=infh
                                            it=it+1
                                        END DO

                                        !  check the amount of infiltrating water in this horizon (infh)
                                        !  (may not exceed its refillable porosity)
                                        !  and for entire profile (infh+fillup=infall)
                                        !             infh=min(infh,namic*horiz_thickness(tcid_instance2,i,h))
                                        !             infall=infh+fillup
                                        !             infall=min(infall,input)
                                        !             infh=infall-fillup

                                        infh=MIN(infh,INPUT)
                                        tempx=fillup+namic*horiz_thickness(tcid_instance2,i,h)	!Till: maximum refillable porosity until current horizon [mm]
                                        IF (infh > tempx) THEN	!Till: maximum refillable porosity should not be exceeded
                                            infh=tempx
                                        END IF

                                        !  remaining water from above
                                        remain=INPUT-infh		!Till: water that has neither been redistributed nor infiltrated [mm]

                                        IF (n_iter==1) THEN !Till first iteration
                                        ELSE !Till: second iteration
                                            !  check if remain is smaller than what has been estimated as
                                            !  inf.excess runoff in preliminary infiltration estimation in 2.2a
                                            !  (this may occur in some special cases, might also be inaccuracy of
                                            !   iterative routine above)
                                            !  in order to not confuse water balance with other SVCs, this previous
                                            !  estimate has to be at least the value of remain

                                            test=qmerk2(i)/(tcarea2*1.e3*frac_svc(i,tcid_instance2))	!Till: "remain" of previous iteration, convert m3 to mm
                                            IF (test > remain) THEN
                                                !	               write(*,*) 'test,remain',tcid_instance2,i,test,remain
                                                !               infall=infall-(test-remain)
                                                !infh=infh-(test-remain)  !Till: reduzier
                                                infh=infh-(test-remain)  !Till: reduce infiltration to ensure sufficient contribution to redistributable water in the entire TC
                                                remain=test
                                                IF (infh < 0.) THEN		!Till: infiltration couldn't be reduced enough - this will lead to a water balance problem. Probably another iteration should be conducted?
                                                    WRITE(*,*) 'inf excess problem, tcid,i',tcid_instance2,i,infh
                                                    !call pause1
                                                END IF
                                            END IF
                                        END IF !Till: branching between iterations


                                        !  check if infiltration into macropores occurs
                                        !  this can occur if current horizon has macropores
                                        !  soil column's macropores are being filled starting with lowest
                                        !  horizon with macropores

                                        IF (macro(i,h) > 0.) THEN
                                            IF (remain > infh) THEN
                                                hmerk=0
                                                testi2=1
                                                j=h
                                                DO WHILE (testi2 == 1 .AND. j <= nbrhori(soilid))
                                                    IF (macro(i,j) > 0.001) THEN
                                                        hmerk=j
                                                    ELSE IF (macro(i,j) <= 0.001) THEN
                                                        testi2=0
                                                    END IF
                                                    j=j+1
                                                END DO
                                                IF (hmerk >= h) THEN
                                                    DO j=hmerk,h,-1
                                                        IF (macro(i,j) > 0. .AND. tempx > 0.0001) THEN
                                                            infmac(j)=MIN(tempx,macro(i,j))
                                                            tempx=tempx-infmac(j)
                                                        END IF
                                                    END DO
                                                END IF
                                            END IF
                                        END IF



                                        !  if remain is still larger than what can infiltrate into
                                        !  micro- and macropores -> surface runoff
                                        IF (remain > 0.) THEN
                                            IF (n_iter==1) THEN !Till first iteration
                                                DO j=1,nbr_svc(tcid_instance2)
                                                    IF (j /= i) THEN
                                                        temp3=remain*tcarea2*1.e3*frac_svc(i,tcid_instance2)*frac_svc(j,tcid_instance2)
                                                        qmerk(j)=qmerk(j)+temp3	!Till: assign this part of the excess runoff to other SVCs for the next iteration [m³]
                                                        qmerk2(i)=qmerk2(i)+temp3 !Till: sum up the redistributed amount [m³]
                                                    END IF
                                                END DO
                                            ELSE !Till: second iteration
                                                IF (dolatsc) THEN
                                                    q_surf_out=q_surf_out + remain*tcarea2*1.e3*frac_svc(i,tcid_instance2)-qmerk2(i)
                                                    hortf=         hortf  +(remain*tcarea2*1.e3*frac_svc(i,tcid_instance2)-qmerk2(i))
                                                    !only for output
                                                    !hortsc(day,i)=hortsc(day,i)+(remain*tcarea2*1.e3*frac_svc(i,tcid_instance2)-qmerk2(i))/(tcarea2*1.e3*frac_svc(i,tcid_instance2))
                                                ELSE
                                                    q_surf_out = q_surf_out + remain*tcarea2*1.e3*frac_svc(i,tcid_instance2)
                                                    hortf      = hortf      + remain*tcarea2*1.e3*frac_svc(i,tcid_instance2)
                                                    !only for output
                                                    !hortsc(day,i)=hortsc(day,i)+remain
                                                END IF


                                                !   updating soil moisture distribution among horizons
                                                !   within soil class

                                                !   horizons above inf.excess horizon
                                                DO j=1,h
                                                    tempx=MIN(na(i,j)*horiz_thickness(tcid_instance2,i,j),infh)
                                                    horithact(tcid_instance2,i,j)=horithact(tcid_instance2,i,j)+tempx

                                                    infh=infh-tempx
                                                END DO
                                                IF (infh > 0.001) THEN
                                                    WRITE(*,*) 'infall zu gross,tcod2,i',tcid_instance2,i,infh
                                                    !						  call pause1
                                                END IF

                                                !   all horizons if macropore infiltration occured
                                                IF (sum(infmac(:)) > 0.) THEN
                                                    DO j=1,nbrhori(soilid)
                                                        horithact(tcid_instance2,i,j)= horithact(tcid_instance2,i,j)+infmac(j)
                                                    END DO
                                                END IF
                                            END IF !Till: branching between iterations
                                        END IF

                                        !**  ts is larger than time step -> no saturation excess in this horizon
                                        !    (infiltrating water should then not be enough for saturation of this horizon
                                        !     -> hopefully fillup < input for next horizon)
                                    ELSE
                                        h=h+1
                                    END IF

                                    !**  saturation depth is larger than soil depth down to this horizon
                                ELSE
                                    h=h+1
                                END IF

                                !**   if input intensity is smaller than conductivity
                                !**   -> no infiltration excess in this horizon
                            ELSE
                                h=h+1
                            END IF

                            !**   input is larger than refillable porosity of higher horizons,
                            !**   but time for refill of upper horizons exceeds timestep
                            !**   assumption is the that after complete fill upper horizons
                            !**   inf excess surface runoff occurs (horton horizon is h-1)
                            !**   (no macropore flow in this case yet included)
                        ELSE IF (INPUT > fillup .AND. tshup > 1.) THEN
                            horton=h-1
                            remain=INPUT-fillup
                            IF (remain > 0.) THEN
                                DO j=1,nbr_svc(tcid_instance2)
                                    IF (j /= i) THEN
                                        temp3=remain*tcarea2*1.e3*frac_svc(i,tcid_instance2)*frac_svc(j,tcid_instance2)
                                        qmerk(j)=qmerk(j)+temp3
                                        qmerk2(i)=qmerk2(i)+temp3
                                    END IF
                                END DO
                            END IF

                            IF (n_iter==1) THEN !Till first iteration
                            ELSE !Till: second iteration

                                !  check if remain is smaller than what has been estimated as
                                !  inf.excess runoff in preliminary infiltration estimation in 2.2a


                                !  (this may occur in some special cases, might also be inaccuracy of
                                !   iterative routine above)
                                !  in order to not confuse water balance with other SVCs, this previous
                                !  estimate has to be at least the value of remain

                                test=qmerk2(i)/(tcarea2*1.e3*frac_svc(i,tcid_instance2))
                                IF (test > remain) THEN
                                    WRITE(*,*) 'inf excess problem remin,tcid,i',tcid_instance2,i, test,remain
                                    !            call pause1
                                END IF

                                !  check if infiltration into macropores occurs
                                !  this can occur if current horizon has macropores
                                !  soil column's macropores are being filled starting with lowest
                                !  horizon with macropores

                                IF (macro(i,h) > 0.) THEN
                                    IF (remain > 0.) THEN
                                        hmerk=0
                                        testi2=1
                                        j=h
                                        DO WHILE (testi2 == 1 .AND. j <= nbrhori(soilid))
                                            IF (macro(i,j) > 0.001) THEN
                                                hmerk=j
                                            ELSE IF (macro(i,j) <= 0.001) THEN
                                                testi2=0
                                            END IF
                                            j=j+1
                                        END DO
                                        IF (hmerk >= h) THEN
                                            DO j=hmerk,h,-1
                                                IF (macro(i,j) > 0. .AND. tempx > 0.0001) THEN
                                                    infmac(j)=MIN(tempx,macro(i,j))
                                                    tempx=tempx-infmac(j)
                                                END IF
                                            END DO
                                        END IF
                                    END IF
                                END IF


                                !  if remain is still larger than what can infiltrate into
                                !  micro- and macropores -> surface runoff
                                IF (remain > 0.) THEN
                                    IF (dolatsc) THEN
                                        q_surf_out=q_surf_out+ remain*tcarea2*1.e3*frac_svc(i,tcid_instance2)-qmerk2(i)
                                        hortf     =hortf     +(remain*tcarea2*1.e3*frac_svc(i,tcid_instance2)-qmerk2(i))
                                        !hortsc(day,i)=hortsc(day,i)+(remain*tcarea2*1.e3*frac_svc(i,tcid_instance2)-qmerk2(i))/(tcarea2*1.e3*frac_svc(i,tcid_instance2))
                                    ELSE
                                        q_surf_out=q_surf_out+remain*tcarea2*1.e3*frac_svc(i,tcid_instance2)
                                        hortf     =hortf     +remain*tcarea2*1.e3*frac_svc(i,tcid_instance2)
                                        !only for output
                                        !hortsc(day,i)=hortsc(day,i)+remain
                                    END IF
                                END IF

                                !   updating soil moisture distribution among horizons
                                !   within soil class

                                !   horizons down to that with infiltration excess (which is
                                !   in this case completely filled) and no water to lower horizons
                                DO j=1,horton
                                    horithact(tcid_instance2,i,j)=thetas(soilid,j)*horiz_thickness(tcid_instance2,i,j)

                                END DO

                                !   all horizons if macropore infiltration occured
                                IF (sum(infmac(:)) > 0.) THEN
                                    DO j=1,nbrhori(soilid)
                                        horithact(tcid_instance2,i,j)= horithact(tcid_instance2,i,j)+infmac(j)

                                    END DO
                                END IF

                            END IF !Till: branching between iterations


                            !**   input is smaller than refillable porosity of higher horizons
                        ELSE
                            h=99
                        END IF

                        !    end of loop for horizons
                    END DO


                    !!for debugging - remove
                    !DO k=1,nbr_svc(tcid_instance2)
                    !	DO h=1,nbrhori(id_soil_intern(k,tcid_instance2))
                    !		IF (horithact(tcid_instance2,k,h) - thetas(id_soil_intern(k,tcid_instance2),h)* horiz_thickness(tcid_instance2,k,h)> 0.1 ) THEN	!Till: if water content exceeds saturation...
                    !			write(*,*)"vor horton=0(",n_iter,"): Oversaturated horizon in TC/svc/horizon ",tcid_instance2,k,h
                    !			call pause1
                    !		END iF
                    !	END DO
                    !END DO
                    !!for debugging - remove end
                    !



                    !   do the following only if no infiltration excess occured in any horizon
                    IF (horton == 0) THEN
                        infup=INPUT
                        IF (n_iter==1) THEN !Till first iteration
                            horitemp(1:maxhori)=horithact(tcid_instance2,i,1:maxhori)		!Till: get current soil moisture of each horizon of current soil in current TC
                            !   higher horizons that are (partly) filled with water
                            j=1
                            DO WHILE (j <= nbrhori(soilid) .AND.  &		!Till: fill up all non-saturated horizons from the top, until all water is distributed or all horizons are full
                                infup > 0.001 .AND. j /= hsat(i))
                                tempx=MIN(infup,na(i,j)*horiz_thickness(tcid_instance2,i,j))		!Till: maximum amount of water that can go into current horizon
                                !horitemp(j)=horitemp(j) +tempx										!possibly obsolete, horitemp is never used later !Till: increase soil water content
                                infup=infup-tempx													!Till: decrease amount of water available for lower horizons
                                j=j+1
                            END DO

                            !   update soil moisture of terrain component with infiltrated water
                            !   of this soil class

                            tempx=inputrem(i)			!Till: water available for infiltration and runoff[m3] = INPUT [mm]
                            inf=INPUT-infup			!Till: water infiltrated during timestep [mm]

                            tempx=tempx-inf*tcarea2*1.e3*frac_svc(i,tcid_instance2)	!Till: remaining water not infiltrated during this iteration step
                            IF (INPUT-inf > 0.001) THEN
                                DO j=1,nbr_svc(tcid_instance2)		!Till: distribute any remaining water among other SVCs
                                    IF (j /= i) THEN
                                        temp3=tempx*frac_svc(j,tcid_instance2)	!Till: distribute according to areal fraction of SVCs
                                        qmerk(j)=qmerk(j)+temp3					!Till: increase surplus for other SVCs
                                        qmerk2(i)=qmerk2(i)+temp3				!Till: increase amount of water redirected away from this SVC
                                    END IF
                                END DO
                            END IF
                        ELSE !Till: second iteration


                            !   higher horizons that are (partly) filled with water
                            j=1
                            DO WHILE (j <= nbrhori(soilid) .AND.  &
                                infup > 0.001 .AND. j /= hsat(i))
                            tempx=MIN(infup,na(i,j)*horiz_thickness(tcid_instance2,i,j))		!Till: maximum amount of water that can go into current horizon
                            horithact(tcid_instance2,i,j)= horithact(tcid_instance2,i,j) + tempx	!Till: increase soil water content
                            infup=infup-tempx													!Till: decrease amount of water available for lower horizons
                            j=j+1
                            END DO


                            !				!for debugging - remove
                            !				DO k=1,nbr_svc(tcid_instance2)
                            !					DO h=1,nbrhori(id_soil_intern(k,tcid_instance2))
                            !						IF (horithact(tcid_instance2,k,h) - thetas(id_soil_intern(k,tcid_instance2),h)* horiz_thickness(tcid_instance2,k,h)> 0.1 ) THEN	!Till: if water content exceeds saturation...
                            !							write(*,*)"update sm(",n_iter,"): Oversaturated horizon in TC/svc/horizon ",tcid_instance2,k,h
                            !							call pause1
                            !						END IF
                            !					END DO
                            !				END DO
                            !				!for debugging - remove end


                            !   update soil moisture of terrain component with infiltrated water
                            !   of this soil class


                            inf=INPUT-infup				!Till: water infiltrated during timestep [mm]
                            inputrem(i)=inputrem(i)-inf*tcarea2*1.e3*frac_svc(i,tcid_instance2)		!Till: remaining water not infiltrated during this iteration step

                            IF (INPUT-inf > 0.001) THEN
                                IF (dolatsc) THEN
                                    q_surf_out=q_surf_out+inputrem(i)-qmerk2(i)
                                    if (q_surf_out<0.) then
                                        write(*,*) 'das ist mist'
                                    end if
                                    !only for output
                                    !hort2sc(day,i)=hort2sc(day,i)+(inputrem(i)-qmerk2(i))/(tcarea2*1.e3*frac_svc(i,tcid_instance2))
                                    IF (dolattc .AND. tc_counter2 == nbrterrain(i_lu)) THEN
                                        !c             hortf=hortf+(inputrem(i)-qmerk2(i))/(tcarea2*1.E3)
                                    END IF

                                ELSE
                                    q_surf_out=q_surf_out+inputrem(i)
                                    !only for output
                                    !hort2sc(day,i)=hort2sc(day,i)+inputrem(i)/ (tcarea2*1.e3*frac_svc(i,tcid_instance2))
                                    IF (dolattc .AND. tc_counter2 == nbrterrain(i_lu)) THEN
                                        !c             hortf=hortf+inputrem(i)/(tcarea2*1.E3)
                                    END IF
                                END IF	!if dolat
                            END IF	!IF (INPUT-inf > 0.001)


                            !  determine new saturated fraction of this SVC after infiltration
                            !  not required here, as satfrac is not needed in evaporation and
                            !  percolation routines
                            !  new satfrac is determined at the end of soilwat


                        END IF !Till: branching between iterations
                    END IF !(horton == 0)
                END IF	!    end of question for unsaturated soil component
            END IF	!    end of question if input > 0

            !	!for debugging - remove
            !	DO k=1,nbr_svc(tcid_instance2)
            !		DO h=1,nbrhori(id_soil_intern(k,tcid_instance2))
            !			IF (horithact(tcid_instance2,k,h) - thetas(id_soil_intern(k,tcid_instance2),h)* horiz_thickness(tcid_instance2,k,h)> 0.1 ) THEN	!Till: if water content exceeds saturation...
            !				write(*,*)"end loop svc(",n_iter,"): Oversaturated horizon in TC/svc/horizon ",tcid_instance2,k,h
            !				call pause1
            !			END IF
            !		END DO
            !	END DO
            !	!for debugging - remove end


        END DO	!    end of loop for all soil components


        !!for debugging - remove
        !DO i=1,nbr_svc(tcid_instance2)
        !	DO h=1,nbrhori(id_soil_intern(i,tcid_instance2))
        !		IF (horithact(tcid_instance2,i,h) - thetas(id_soil_intern(i,tcid_instance2),h)* horiz_thickness(tcid_instance2,i,h)> 0.1 ) THEN	!Till: if water content exceeds saturation...
        !			write(*,*)"after infiltration excess(",n_iter,"): Oversaturated horizon in TC/svc/horizon ",tcid_instance2,i,h
        !			call pause1
        !		END IF
        !	END DO
        !END DO
        !!for debugging - remove end



    END DO !Till: loop over two iterations



    !1:nbr_svc(tcid_instance2)
    !thsc(tc_counter2,day,1:nbr_svc(tcid_instance2))=thsc(tc_counter2,day,1:nbr_svc(tcid_instance2))+sum(horithact(tcid_instance2,1:nbr_svc(tcid_instance2),:))/dt_per_day
    thact=sum( sum(horithact(tcid_instance2,1:nbr_svc(tcid_instance2),:), dim=2) * frac_svc(1:nbr_svc(tcid_instance2),tcid_instance2) )









    !** -------------------------------------------------------------------------
    !** (3) EVAPOTRANSPIRATION

    !   Dual-layer (canopy+soil surface) and dual-source (vegetation+bare soil)
    !   approach according to Brenner & Incoln (1997)
    !   reduces to the classical Shuttleworth & Wallace approach (1985)
    !   if fractional vegetation cover is 1

    !   canopy resistance is function of given stomatal resistance and LAI
    !   canopy resistance is modified as function of water suction in root zone,
    !   water vapour pressure deficit (Hannan et al.)

    !   soil surface resistance is function of soil moisture in uppermost horizon
    !   (according to  )

    !   separate calculation for night and daytime (Schulla)

    tcaet=0.


    !   compare soil depth versus root depth
    !   determine number of horizons within root zone (nbrrooth)
    !   and fraction of depth of lowest horizon to be reached by roots (hfrac)

    nbrrooth(:)=0
    hfrac(:)=0.
    DO i=1,nbr_svc(tcid_instance2)
        soilid=id_soil_intern(i,tcid_instance2)
        tempx=0.                                         !Till: for summing up soil layer thickness

        temp3=rootd_act(id_veg_intern(i,tcid_instance2)) !Till: current seasonal root depth
        if (temp3 ==0) then
            nbrrooth(i)=1
            hfrac(i)=0.01 !Till: prevent later division by 0
        else
            h=1
            DO WHILE (nbrrooth(i) == 0 .AND. h <= nbrhori(soilid))
                tempx=tempx+horiz_thickness(tcid_instance2,i,h)
                IF (tempx >= temp3) THEN
                    nbrrooth(i)=h
                    IF (h > 1) THEN
                        hfrac(i)=(temp3-sum(horiz_thickness(tcid_instance2,i,1:h-1)))/ horiz_thickness(tcid_instance2,i,h)
                    ELSE
                        hfrac(i)=temp3/horiz_thickness(tcid_instance2,i,h)
                    END IF
                ELSE IF (h == nbrhori(soilid)) THEN
                    nbrrooth(i)=h
                    hfrac(i)=1.
                END IF
                h=h+1
            END DO

            !    ansonsten: mist:
            IF (nbrrooth(i) == 0.) THEN
                WRITE(*,*) 'no roots ?, tcid_instance2,i',tcid_instance2,i
                !call pause1
            END IF

            IF (hfrac(i)==0.) THEN
                WRITE(*,*) 'calculated rooted fraction is 0, set to 0.001, tcid_instance2,i',tcid_instance2,i
                hfrac(i)=0.01
            END IF
        end if

    END DO


    !  for all Soil-Vegetation-components
    DO i=1,nbr_svc(tcid_instance2)
        soilid=id_soil_intern(i,tcid_instance2)
        etpfrac=frac_svc(i,tcid_instance2)
        aktnfk(:)=0.
        vangen(:)=0.
        suction(:)=0.
        resf3(:)=0.

        !  calculate water potential (suction) for each horizon of root zone (hPa)
        DO h=1,nbrrooth(i)
            IF (horithact(tcid_instance2,i,h)/horiz_thickness(tcid_instance2,i,h) > pwpsc(tcid_instance2,i,h)) THEN

                vangen(h)=(horithact(tcid_instance2,i,h)/ horiz_thickness(tcid_instance2,i,h)-  &
                    thetar(soilid,h))/ (thetas(soilid,h)-thetar(soilid,h))	!Till: compute relative saturation
                vangen(h)=MAX(vangen(h),0.)
                vangen(h)=MIN(vangen(h),1.)									!Till: relative saturation must be within 0..1

                if (vangen(h)==1.) then										!Till: complete saturation...
                    suction(h)=0.												!Till: suction is 0 (this produces NaN otherwise)
                else
                    suction(h)=((1./(vangen(h)** (1./porem(soilid,h))) -1.)**(1./(poresz(soilid,h)+1.)))/ (1./bubble(soilid,h))
                end if



                ! numeric overflow test values -- old obsolete check replaced by method below
                !      temp3=(thetar(soilid,h)+0.001- thetar(soilid,h))/  &
                !          (thetas(soilid,h)- thetar(soilid,h))
                !      temp3=MAX(temp3,0.)
                !      temp4=(((1./(temp3** (1./porem(soilid,h)))  &
                !          -1.)**(1./(poresz(soilid,h)+1.)))/ (1./bubble(soilid,h)))



                IF (suction(h) > huge(suction(h))-1) THEN
                    suction(h)=huge(suction(h))-1.		!Till: prevents infinite suction values
                END IF


                !  soil water stress factor for water uptake by plants
                !  in each horizon (0 <= resf3 <= 1)
                !  linear with respect to suction between critical suction
                !  (start of stomata closure) and
                !  wilting point (total stomata closure)
                resf3(h)=1.-(suction(h)-wstressmin(id_veg_intern(i,tcid_instance2)))/  &
                    (wstressmax(id_veg_intern(i,tcid_instance2))- wstressmin(id_veg_intern(i,tcid_instance2)))
                IF (resf3(h) < 0.01) THEN
                    resf3(h)=0.01
                END IF
                IF (resf3(h) > 1.) THEN
                    resf3(h)=1.
                END IF

                !  if horizon is at PWP
            ELSE
                resf3(h)=0.01
            END IF
        END DO

        !   horizon-depth weighted mean of stress parameter
        !   (harmonic mean will result in more severe restriction)
        temp3=sum( resf3(1:nbrrooth(i)-1)*horiz_thickness(tcid_instance2,i,1:nbrrooth(i)-1) )
        temp4=sum(                        horiz_thickness(tcid_instance2,i,1:nbrrooth(i)-1) )


        h=nbrrooth(i)
        temp3=temp3+resf3(h)*horiz_thickness(tcid_instance2,i,h)*hfrac(i)
        temp4=temp4+horiz_thickness(tcid_instance2,i,h)*hfrac(i)

        def=temp3/temp4

        if (isNaN(def)) then
            write(*,*)'NAN value produced! (2)'
            def=1.				!this happens if a vegetation class with "no vegetation" is encountered
            !ii: fix this more elegantly
        end if

        !   use stress parameter of wettest horizon only
        !       def=maxval(resf3)

        !  soil surface resistance (Taylor, 2000)
        !       facw=2.4*(horithact(tcid_instance2,i,1)/
        !     .           horiz_thickness(tcid_instance2,i,1))**(-1.*1.9)
        !  soil surface resistance (mean Domingo et al. 1999)
        !  (based on grav. water content (bulk density = 1.5)
        if (horithact(tcid_instance2,i,1)<=0.001) then
            facw=huge(facw)				!Till: if soil moisture is very low, set resistance very high (otherwise this caused problems when horithact was 0)
        else
            facw=26.*(horithact(tcid_instance2,i,1)/ horiz_thickness(tcid_instance2,i,1)/1.5)**(-1.)
        end if

        if (isnan(facw)) then
            write(*,*)'NAN value produced (3)'
        end if

        !   actual nFK of each horizon in root zone
        !   (lowest horizon is confined to root zone part)
        DO h=1,nbrrooth(i)
            tempx=MIN(horithact(tcid_instance2,i,h), soilfc63(soilid,h)*  &
                horiz_thickness(tcid_instance2,i,h))
            tempx=tempx-pwpsc(tcid_instance2,i,h)* horiz_thickness(tcid_instance2,i,h)
            tempx=MAX(0.0,tempx)
            IF (h == nbrrooth(i)) THEN
                tempx=tempx*hfrac(i)
            END IF
            aktnfk(h)=tempx
        END DO

        !  subroutine for calculation of actual evapotranspiration
        evapveg=0.
        evaps=0.


        IF (dohour) THEN

            !	evapt=-7.
            !	evapveg=-7.
            !	evaps=-7.
            CALL etp_soil_hour(hh,i_subbas2,id_veg_intern(i,tcid_instance2),def,facw,  &
                evapt,evapveg,evaps, height_act,lai_act,alb_act,rsfinday)
            !only for output
            !resistsc(day,i)=resistsc(day,i)+rsfinday/dt_per_day
        ELSE
            CALL etp_soil(i_subbas2,id_veg_intern(i,tcid_instance2),def,facw, evapt,evapveg,evaps,  &
                height_act,lai_act,alb_act,rsfinday)
            !only for output
            !resistsc(day,i)=resistsc(day,i)+rsfinday/dt_per_day
        END IF



        if (isnan(evaps)) then
            write(*,*)'NAN value produced! (4)'
            !	write(*,*)i_subbas2,id_veg_intern(i,tcid_instance2),def,facw, evapt,evapveg,evaps,  &
            !	height_act,lai_act,alb_act,rsfinday

        end if


        ! sometimes there may be some net negative allday transpiration (dew deposit)
        ! if daytime transpiration is 0 due to no active vegetation
        ! allday transpiration is set to 0 then
        IF (evapveg < 0.) THEN
            evapveg=0.
        END IF
        IF (evaps < 0.) THEN
            evaps=0.
        END IF
        !evaps=0
        evapt=evaps+evapveg

        if (isnan(evapt)) then
            evaps=1.
            write(*,*)'ERROR: problems with evaps1'
            stop
        end if


        IF (evapt > 0.) THEN
            !  soil and vegetation fractions of total evap

            evapsfrac=evaps/(evaps+evapveg)	!Till: fraction of soil evaporation in relation to total ET

            !  soil and vegetation evaporation plus interception evaporation
            !  may in maximum be potential evaporation of this day
            !  (new version)

            evapt=evapt*(1.-aetred(i))		!Till: reduce total evapotranspiration by the ratio of evap from interception and Epot

            !	IF (evapt+intc_evap(i) > etpmax(i)) THEN
            !		evapt=etpmax(i)-intc_evap(i)
            !
            !		if(isnan(evapt)) then		!for debugging only
            !			evapt=etpmax(i)-intc_evap(i)
            !		end if
            !
            !		if(isnan(evaps)) then		!for debugging only
            !			evapt=etpmax(i)-intc_evap(i)
            !		end if
            !    END IF

            evaps  =evapt*evapsfrac	!Till: reduced total ET is redistributed, retaining the original ration between EP and TR
            evapveg=evapt-evaps

            if (isnan(evaps)) then
                evaps=1.
                write(*,*)'ERROR: problems with evaps'
                stop
            end if


        END IF

        !  transpiration of current soil component
        !  (must not be larger than actual nfk)
        scaet=MIN(evapveg,sum(aktnfk(:)))
        evapveg=scaet

        !  ETP from each horizon
        !  fraction of total ETP from each horizon according
        !  to relative distribution of aktnfk (in mm ?) between horizons

        tempx=sum(aktnfk(:))
        !only for output
        !nfksc(day,i)=nfksc(day,i)+tempx/dt_per_day
        etpvert(:)=0.
        temp3=0.
        IF (tempx > 0.0 .AND. scaet > 0.) THEN
            !       if (tempx > 0.0) then
            !  weighting factors
            DO h=1,nbrrooth(i)
                IF (aktnfk(h) == 0.) THEN
                    etpvert(h)=0.
                ELSE IF (suction(h) <=  &
                    wstressmin(id_veg_intern(i,tcid_instance2))) THEN
                etpvert(h)=horiz_thickness(tcid_instance2,i,h)
                IF (h == nbrrooth(i)) THEN
                    etpvert(h)=horiz_thickness(tcid_instance2,i,h)*hfrac(i)
                END IF
                ELSE
                    etpvert(h)=horiz_thickness(tcid_instance2,i,h)*  &
                        (suction(h)-wstressmin(id_veg_intern(i,tcid_instance2)))/  &
                        (wstressmax(id_veg_intern(i,tcid_instance2))- wstressmin(id_veg_intern(i,tcid_instance2)))
                    IF (h == nbrrooth(i)) THEN
                        etpvert(h)=horiz_thickness(tcid_instance2,i,h)*hfrac(i)*  &
                            (suction(h)-wstressmin(id_veg_intern(i,tcid_instance2)))/  &
                            (wstressmax(id_veg_intern(i,tcid_instance2))- wstressmin(id_veg_intern(i,tcid_instance2)))

                    END IF
                END IF
                temp3=temp3+etpvert(h)
            END DO

            !  update soil moisture of all horizons in root zone
            !  by plant transpiration

            !  available field capacity of a horizon might be slightly smaller than
            !  expected evap attributed to this horizon by above weighting
            !  -> then evapveg of lower horizon is taken to be higher by this amount
            !  if problem occurs for lowest horizon
            !  -> start again with uppermost horizon
            tempx=0.
            temp4=scaet
            DO h=1,nbrrooth(i)
                IF (aktnfk(h) > 0.) THEN
                    horithact(tcid_instance2,i,h)=horithact(tcid_instance2,i,h)  &
                        -MIN(evapveg*etpvert(h)/temp3+tempx,aktnfk(h))


                    temp4=temp4- MIN(evapveg*etpvert(h)/temp3+tempx,aktnfk(h))
                    IF (aktnfk(h) < evapveg*etpvert(h)/temp3+tempx) THEN
                        tempx=evapveg*etpvert(h)/temp3+tempx-aktnfk(h)
                    ELSE
                        tempx=0.
                    END IF
                END IF
            END DO
            IF (tempx > 0.) THEN
                IF (h == nbrrooth(i)+1) THEN
                    DO h=1,nbrrooth(i)-1
                        IF (aktnfk(h) > 0.) THEN
                            IF (tempx > 0.001) THEN
                                horithact(tcid_instance2,i,h)=horithact(tcid_instance2,i,h) -MIN(tempx,aktnfk(h))

                                temp4=temp4- MIN(tempx,aktnfk(h))
                                IF (aktnfk(h) < tempx) THEN
                                    tempx=tempx-aktnfk(h)
                                ELSE
                                    tempx=0.
                                END IF
                            END IF
                        END IF
                    END DO
                END IF
                IF (tempx > 0.001) THEN
                    WRITE(*,*) 'evapveg  > nfk-availability !'
                    WRITE(*,*) tempx,h,maxhori
                    !call pause1
                    scaet=scaet-tempx
                END IF
            END IF
            IF (temp4-scaet > 0.001) THEN
                WRITE(*,*) 'temp4>scaet!!'
                !call pause1
            END IF
        END IF

        if (isnan(evaps)) then
            evaps=1.
        end if

        !  update soil moisture of uppermost horizon by evaporation from soil surface
        IF (evaps > 0) THEN
            evaps=MIN(horithact(tcid_instance2,i,1)- thetar(soilid,1)*  &
                horiz_thickness(tcid_instance2,i,1),evaps)
            evaps=MAX(0.,evaps)
            horithact(tcid_instance2,i,1)=horithact(tcid_instance2,i,1)-evaps
        END IF


        !  contribution of this SVC to evapotranspiration of entire TC
        !  and end of loop for all Soil-Veg-components
        tcaet=tcaet+(scaet+intc_evap(i)+evaps)*etpfrac
        tcsoilet=tcsoilet+evaps*etpfrac

        !only for output
        !aetsc(day,i)=aetsc(day,i)+scaet
        !intetsc(day,i)=intetsc(day,i)+intc_evap(i)
        !soiletsc(day,i)=soiletsc(day,i)+evaps
        !aettotsc(day,i)=aettotsc(day,i)+(scaet+intc_evap(i)+evaps)


        !  determine new saturated fraction of this SVC after evapotranspiration
        !  frac_sat value is relative to TC, not to SVC !!

        frac_old=frac_sat(tcid_instance2,i)
        tempth=sum(horithact(tcid_instance2,i,:))
        testi=0

        IF (tempth <= tctheta_s(1,1,tcid_instance2,i)) THEN
            frac_sat(tcid_instance2,i)=0.
        ELSE
            testi=0
            tempalt=tctheta_s(1,1,tcid_instance2,i)
            DO j=1,4
                IF (testi == 0) THEN
                    tempx=tempalt+(tctheta_s(j+1,1,tcid_instance2,i)- tctheta_s(j,1,tcid_instance2,i))*  &
                        (tctheta_s(j+1,2,tcid_instance2,i)- tctheta_s(j,2,tcid_instance2,i))/2.+  &
                        (tctheta_s(j+1,1,tcid_instance2,i)- tctheta_s(j,1,tcid_instance2,i))*  &
                        (tctheta_s(5,2,tcid_instance2,i)- (tctheta_s(j+1,2,tcid_instance2,i)-  &
                        tctheta_s(1,2,tcid_instance2,i)))

                    IF (tempth > tempalt .AND. tempth < tempx) THEN
                        frac_sat(tcid_instance2,i)=((tctheta_s(j,2,tcid_instance2,i)-  &
                            tctheta_s(1,2,tcid_instance2,i))+ (tempth-tempalt)/(tempx-tempalt)*  &
                            (tctheta_s(j+1,2,tcid_instance2,i)- tctheta_s(j,2,tcid_instance2,i)))  &
                            *frac_svc(i,tcid_instance2)
                        testi=100
                    END IF
                    tempalt=tempx
                END IF
            END DO

            !  soil component is completely saturated
            IF (testi == 0) THEN
                frac_sat(tcid_instance2,i)=frac_svc(i,tcid_instance2)
            END IF
        END IF

        !  end loop of all soil-veg-components
    END DO


    !  new soil moisture of terrain component
    watbal=watbal-tcaet

    !  add evaporation from rocky surfaces to evapotranspiration of TC
    !  (has been taken into account for watbal already in interception calculation)
    tcaet=tcaet+rockintc



    !** -------------------------------------------------------
    !** (4) Percolation and lateral flow


    !** free drainage water above field capacity yields
    !   percolation to deeper soil horizon,
    !   vertical recharge to deep groundwater and
    !   lateral subsurface flow to lower terrain component

    qgw=0.
    deepqgw=0.
    thfree(:,:)=0.

    !** do for all soil components

    DO i=1,nbr_svc(tcid_instance2)

        soilid=id_soil_intern(i,tcid_instance2)
        l2(:)=0.0
        l2eff(:)=0.0
        percol(:)=0.
        percolmac(:)=0.
        vangen(:)=0.
        conduns(:)=0.
        suction(:)=0.

        !** shrinkage
        macro(i,:)=0.
        IF (doshrink == 1) THEN
            DO h=1,nbrhori(soilid)
                IF (shrink(soilid,h) == 1.0) THEN

                    tempx=((1./(1.+(3100.0/bubble(soilid,h))**  &
                        (poresz(soilid,h)+1.)))** porem(soilid,h))*  &
                        (thetas(soilid,h)/(1.-coarse(soilid,h))-  &
                        thetar(soilid,h)/(1.-coarse(soilid,h)))  &
                        +thetar(soilid,h)/(1.-coarse(soilid,h))
                    tempx=tempx*horiz_thickness(tcid_instance2,i,h)*(1.-coarse(soilid,h))
                    IF (horithact(tcid_instance2,i,h) < tempx) THEN
                        macro(i,h)=((tempx-horithact(tcid_instance2,i,h))/  &
                            (tempx-thetar(soilid,h)))*0.5* thetas(soilid,h)*  &
                            horiz_thickness(tcid_instance2,i,h)
                    END IF
                ELSE IF (shrink(soilid,h) == 0.0) THEN
                    macro(i,h)=0.01*horiz_thickness(tcid_instance2,i,h)
                    !           else
                    !            macro(i,h)=0.0025*horiz_thickness(tcid_instance2,i,h)
                END IF
            END DO
        END IF


        !**   do for all soil horizons

        DO h=1,nbrhori(soilid)

            thfree(i,h)=horithact(tcid_instance2,i,h)- soilfc(soilid,h)*  &		!Till: compute amount of water above field capacity [mm]
            horiz_thickness(tcid_instance2,i,h)

            IF (thfree(i,h) > 0.0) THEN

                !  calculate unsaturated conductivity of horizon
                !  (calculation of vangen: coarse fragment reduction equals out)

                vangen(h)=(horithact(tcid_instance2,i,h)/ horiz_thickness(tcid_instance2,i,h)-  &
                    thetar(soilid,h))/ (thetas(soilid,h)-  &
                    thetar(soilid,h))				!Till: compute relative saturation
                vangen(h)=MAX(vangen(h),0.)
                vangen(h)=MIN(vangen(h),1.)
                conduns(h)=(k_sat(soilid,h)/dt_per_day*(vangen(h)**0.5)*  &
                    (1.-(1.-(vangen(h)** (1./porem(soilid,h)))  &
                    )**porem(soilid,h))**2)* (1.-coarse(soilid,h))

                ! catch unrealistic conduns value due to possibly unrealistic soil parameters
                ! would lead to floating point exception during percolation calculation
                if(conduns(h) < 1e-12) then
                    write(*,*) 'ERROR: Calculated unsaturated conductivity is zero!'
                    write(*,*) 'Check your soil parameters for soil id ', soilid, ', horizon ', h
                    stop
                end if

                exp_x = -1./(thfree(i,h)/conduns(h))


                IF (h < nbrhori(soilid)) THEN
                    !          percol(h)=k_sat(soilid,h)/dt_per_day*(1.-coarse(soilid,h))

                    IF (exp_x < -15) THEN !Avoid floating point exception
                        percol(h) = thfree(i,h)
                    ELSE
                        percol(h)=thfree(i,h)*(1.-EXP(exp_x))	!Till: compute percolation [mm]
                    END IF

                    !  if thact is greater than field capacity and lower horizons
                    !  have skrinkages or macropores, than rapid percolation
                    !  macropores are filled from below (from deepest horizon with macropores)
                    !  if there is a continous macropore space throughout all horizons, starting
                    !  with the current horizon
                    !  macropores can be filled twice every timestep
                    !  (first time during infiltration, second time during percolation)
                    !  macro is reduced if flow into macropores occured for one horizon
                    !  (available macropores volume is smaller for higher horizons)
                    percolmac(:)=0.
                    IF (macro(i,h+1) > 0.) THEN
                        IF (thfree(i,h) > 0.) THEN
                            hmerk=h
                            testi2=1
                            j=h+1
                            DO WHILE (testi2 == 1 .AND.  &
                                j <= nbrhori(soilid))
                            IF (macro(i,j) > 0.001) THEN
                                hmerk=j
                            ELSE IF (macro(i,j) <= 0.001) THEN
                                testi2=0
                            END IF
                            j=j+1
                            END DO
                            IF (hmerk > h) THEN
                                tempx=thfree(i,h)
                                DO j=hmerk,h+1,-1
                                    IF (macro(i,j) > 0. .AND. tempx > 0.0001) THEN
                                        percolmac(j)=MIN(tempx,macro(i,j))
                                        tempx=tempx-percolmac(j)
                                        macro(i,j)=macro(i,j)-percolmac(j)
                                    END IF
                                END DO
                            END IF
                        END IF
                    END IF

                ELSE IF (h == nbrhori(soilid)) THEN				!Till: deepest horizon
                    IF (gw_flag(i_lu) == 0 .OR. gw_flag(i_lu) == 1) THEN
                        IF (exp_x < -15) THEN !Avoid floating point exception
                            percol(h) = thfree(i,h)
                        ELSE
                            percol(h)=thfree(i,h)*(1.-EXP(exp_x))	!Till: compute percolation [mm]
                        END IF
                        IF (svcbedr(tcid_instance2,i) == 1) THEN		!Till: if there is bedrock...
                            percol(h)=MIN(percol(h),kfsu(i_lu)/dt_per_day) !Till: ...percolation is limited by bedrock conductivity
                        END IF
                    ELSE IF (gw_flag(i_lu) == 99) THEN				!Till: experimental groundwater configuration (undocumented)
                        percol(h)=0.									!Till: no percolation at all
                    END IF
                END IF

                !  lateral flow out of horizon as function of slope gradient of TC, ksat,
                !  slope length and outflow area, which is height of horizon multiplied by
                !  downslope contour length of Soil-Veg-component, leading to
                !  downslope terrain component or to river
                !  slope of terrain component is given in % (100% = 45 degree)
                !  gradient for vertical flow is 1, then if slope=45degree: gradient=0.5
                !  mean slope length in mm
                ! (area of SVC equals out, see notes)


                !  lateral outflow to river for lowest terrain component
                !  only if saturated horizon is above level of river bed
                !  (calibration parameter, given in soter.dat)

                ! (transmission losses only if lateral subsurface inflow ?)



                IF (thfree(i,h) > 0.) THEN !Till: if free water is available, try draining to riverbed

                    ! Calculate temp3 (mm) which is the water content above field capacity

                    !depth (thickness) of saturated zone in this horizon
                    temp3=thfree(i,h)/(thetas(soilid,h)-soilfc(soilid,h))


                    !Calculate temp2 (-), which is the saturated fraction of the horizon above the riverbed
                    !for the lowermost TC and otherwise 1
                    !- For the lowest terrain component lateral flow out of the TC is subsurface runoff to
                    !  the river. This lateral outflow to the river only occurs for the saturated part above
                    !  the river bed.
                    !- For all other TCs lateral flow is lateral inflow to the next lower TC and calculated over
                    !  the entire horizon.
                    !
                    IF (tc_counter2 == nbrterrain(i_lu)) THEN		!Till: if this TC is the most downslope of the entire LU
                        !depth of lower limit of horizon below ground
                        tempx=sum(horiz_thickness(tcid_instance2,i,1:h))
                        !  depth of water table below ground surface
                        temp4=tempx-temp3

                        !  horizon is completely above river bed -> all goes to river
                        IF (tempx <= riverbed(i_lu)) THEN
                            temp2=1.

                            !  saturated zone is completely below river bed -> no runoff
                        ELSE IF (temp4 >= riverbed(i_lu)) THEN
                            temp2=0.

                            !  horizon limit is below, saturated zone above river bed
                        ELSE IF (tempx > riverbed(i_lu) .AND. temp4 < riverbed(i_lu)) THEN
                            temp2=(riverbed(i_lu)-temp4)/temp3		!Till: drain saturated fraction of horizon that is above riverbed
                        ELSE
                            WRITE(*,*) 'ERROR: groundwater level not valid'
                            WRITE(*,*) i,tcid_instance2,soilid
                            STOP
                        END IF

                    ELSE IF (tc_counter2 < nbrterrain(i_lu)) THEN	!Doris, Andreas: if this TC is not the most downslope of the entire LU
                        temp2=1.
                    END IF


                    !Till: compute maximum inflow according to "exchange surface" (revised)
                    !Doris, Till, Andreas: calculation results in wrong unit, l2 should result in mm

                    !!Eq 4.46 saturated depth in m
                    !dsi = temp3*temp2 /1000.
                    !
                    !!Eq 4.45 cross section of lateral flow in m2
                    !Aq  = frac_svc(i,tcid_instance2)*area(i_subbas2)*frac_lu(lu_counter2,i_subbas2)*1000000  &
                    !	  /slength(i_lu)  &
                    !	  *dsi
                    !
                    !!Eq 4.44 lateral flow in m3/day
                    !Q_li= Aq * (k_sat(soilid,h)/1000)/dt_per_day*(1.-coarse(soilid,h))  &
                    !	   * slope(id_tc_type2)/100
                    !
                    !
                    !!lateral flow in mm (convert area from km2 to m2, then convert flow from m to mm)
                    !l2_test1=Q_li/(area(i_subbas2)*frac_lu(lu_counter2,i_subbas2) * frac_svc(i,tcid_instance2) * fracterrain(id_tc_type2)*1000000)*1000


                    !! above lines in one step:
                    !l2(h)=frac_svc(i,tcid_instance2)*area(i_subbas2)*frac_lu(lu_counter2,i_subbas2)*1000000  /  slength(i_lu)  &    ! Aq
                    !	  * temp3*temp2 /1000	&											   ! dsi
                    !	  * (k_sat(soilid,h)/1000)/dt_per_day*(1.-coarse(soilid,h)) * slope(id_tc_type2)/100.  &		   ! Q_li
                    !	  / (area(i_subbas2)*frac_lu(lu_counter2,i_subbas2) * frac_svc(i,tcid_instance2) * fracterrain(id_tc_type2)*1000000)  *  1000  !/area


                    ! above lines in one step and shorter
                    l2(h)= k_sat(soilid,h)/dt_per_day*(1.-coarse(soilid,h))  &		   ! Q_li
                    * temp3*temp2* slope(id_tc_type2)/100. 	&
                        /  slength(i_lu) / fracterrain(id_tc_type2)/1000.


                    !Original equation from Andreas: equal to new corrected equation with extra factor /2.
                    !l2_test2=k_sat(soilid,h)/dt_per_day*(1.-coarse(soilid,h))*  &
                    !		 thfree(i,h)/(thetas(soilid,h)-soilfc(soilid,h)) *temp2  &
                    !		/(slength(i_lu)*fracterrain(id_tc_type2)*1000.) *(slope(id_tc_type2)/100./2.)

                ELSE

                    !l2(h)=k_sat(soilid,h)/dt_per_day*(1.-coarse(soilid,h))*  &	!ii obsolete, since thfree is 0 anyway in this case, just set l2(h)=0





                    !	thfree(i,h)/(thetas(soilid,h)-soilfc(soilid,h))  &
                    !	/(slength(i_lu)*fracterrain(id_tc_type2)*1000.) *(slope(id_tc_type2)/100./2.)
                    !	!temp2=temp2

                    l2(h)=0. !Till: no free draining water

                END IF	! if (thfree(i,h) > 0.)



                ! maximum possible outflow of horizon is defined by its
                ! available moisture content above FC
                test=thfree(i,h)			!Till: this is available
                tempx=percol(h)+l2(h)		!Till: this can potentially drain
                IF (tempx > test) THEN	!Till: percolation and lateral outflow are rescaled to be (in sum) thfree at maximum
                    percol(h)=test*percol(h)/tempx
                    l2(h)=test*l2(h)/tempx
                END IF

                ! if lateral redistribution between SVCs occurs, effective lateral
                ! subsurface flow is function of areal fraction of this SVC in TC;
                ! for the most downslope TC, lateral outflow of TC only after complete
                ! redistribution among SVCs

                !     Andreas/Till 2011-03-01: This part of the code has changed, previously redistribution was only performed
                !               in the case of no groundwater. In the previous version there was additional code for specific
                !               tests in Brazil.
                IF (dolatscsub) THEN
                    lath=INT(sum(horiz_thickness(tcid_instance2,i,1:h))/500.)+1							!Till: compute number of "exchange" horizon
                    latred(tcid_instance2,lath)=latred(tcid_instance2,lath)+&							!Till: compute amount of water to go into exchange horizon [m3]
                    l2(h)*(1.-frac_svc(i,tcid_instance2))* tcarea2*frac_svc(i,tcid_instance2)*1.e3
                    l2eff(h)=l2(h)*frac_svc(i,tcid_instance2)											!Till: assign subsurface flow to river or lower TC
                ELSE
                    l2eff(1:maxhori)=l2(1:maxhori)
                END IF

                ! lateral subsurface outflow of SVC add to outflow of TC
                q_sub_out=q_sub_out+l2eff(h)*tcarea2*frac_svc(i,tcid_instance2)*1.e3

                !  update soil moisture of deeper horizons by macroporeflow
                horithact(tcid_instance2,i, 1:nbrhori(soilid))= horithact(tcid_instance2,i, 1:nbrhori(soilid))+percolmac(1:nbrhori(soilid))

                !  update soil moisture of deeper horizon
                !  check if refillable porosity is larger than percolation
                IF (h < nbrhori(soilid)) THEN
                    tempna=thetas(soilid,h+1)*horiz_thickness(tcid_instance2,i,h+1)- horithact(tcid_instance2,i,h+1)	!Till: compute refillable porosity
                    horithact(tcid_instance2,i,h+1)=horithact(tcid_instance2,i,h+1)+ MIN(percol(h),tempna)				!Till: increase moisture of lower horizon
                    !  update soil moisture of current horizon due to percolation and
                    !  lateral flow
                    horithact(tcid_instance2,i,h)=horithact(tcid_instance2,i,h)-  &
                        (MIN(percol(h),tempna)+l2(h)+sum(percolmac(:)))			!Till: decrease moisture of current horizon due to percolation an lateral flow. why is l2 (instead of l2eff) used here?


                    !  check if current horizon is lowest of root zone
                    !  then percolation is considered as groundwater recharge
                    !  (only for water flow study, not important for water balance of SVC,
                    !   as this water does not leave the profile)
                    IF (h == svcrooth(tcid_instance2,i)) THEN
                        qgw=qgw+MIN(percol(h),tempna)*frac_svc(i,tcid_instance2)
                    END IF

                    !  if no deeper horizon, water leaves the root zone
                    !  (groundwater recharge)
                    !  (output of gwr in mm)
                    !  if landscape unit with groundwater body, percolation gets river runoff,
                    !  independent of position of current terrain component
                ELSE IF (h == nbrhori(soilid)) THEN
                    deepqgw=deepqgw+percol(h)*frac_svc(i,tcid_instance2)	!Till: increase groundwater recharge [mm]
                    watbal=watbal-percol(h)*frac_svc(i,tcid_instance2)
                    horithact(tcid_instance2,i,h)=horithact(tcid_instance2,i,h)-percol(h)-l2(h) !Till: decrease moisture of current horizon due to percolation an lateral flow


                    !only for output
                    !deepgwrsc(day,i)=deepgwrsc(day,i)+percol(h)
                    IF (h == svcrooth(tcid_instance2,i)) THEN
                        qgw=qgw+percol(h)*frac_svc(i,tcid_instance2)
                    END IF
                END IF

                !  end of question for free draininage water
            END IF
            !  end of loop for all horizons
        END DO


        !  update soil moisture of entire terrain component  due to
        !  lateral outflow of current soil-veg component
        watbal=watbal-sum(l2(:))*frac_svc(i,tcid_instance2)

        !  determine new saturated fraction of this SVC after infiltration
        !  and percolation/lateral outflow
        !  frac_sat value is relative to TC, not to SVC !!

        frac_old=frac_sat(tcid_instance2,i)
        tempth=sum(horithact(tcid_instance2,i,:))
        testi=0

        IF (tempth <= tctheta_s(1,1,tcid_instance2,i)) THEN
            frac_sat(tcid_instance2,i)=0.
        ELSE
            testi=0
            tempalt=tctheta_s(1,1,tcid_instance2,i)
            DO j=1,4
                IF (testi == 0) THEN
                    tempx=tempalt+(tctheta_s(j+1,1,tcid_instance2,i)- tctheta_s(j,1,tcid_instance2,i))*  &
                        (tctheta_s(j+1,2,tcid_instance2,i)- tctheta_s(j,2,tcid_instance2,i))/2.+  &
                        (tctheta_s(j+1,1,tcid_instance2,i)- tctheta_s(j,1,tcid_instance2,i))*  &
                        (tctheta_s(5,2,tcid_instance2,i)- (tctheta_s(j+1,2,tcid_instance2,i)-  &
                        tctheta_s(1,2,tcid_instance2,i)))

                    IF (tempth > tempalt .AND. tempth < tempx) THEN
                        frac_sat(tcid_instance2,i)=((tctheta_s(j,2,tcid_instance2,i)-  &
                            tctheta_s(1,2,tcid_instance2,i))+ (tempth-tempalt)/(tempx-tempalt)*  &
                            (tctheta_s(j+1,2,tcid_instance2,i)- tctheta_s(j,2,tcid_instance2,i)))  &
                            *frac_svc(i,tcid_instance2)
                        testi=100
                    END IF
                    tempalt=tempx
                END IF
            END DO

            !  soil component is completely saturated
            IF (testi == 0) THEN
                frac_sat(tcid_instance2,i)=frac_svc(i,tcid_instance2)
            END IF
        END IF

        !  end of loop for all soil components
    END DO


    !only for output
    !thsc(tc_counter2,day,1:nbr_svc(tcid_instance2))=sum(horithact(tcid_instance2,1:nbr_svc(tcid_instance2),:))
    thact=sum(sum(horithact(tcid_instance2,1:nbr_svc(tcid_instance2),:), dim=2)*frac_svc(1:nbr_svc(tcid_instance2),tcid_instance2))




    !** -------------------------------------------------------
    !** (5) Capillary rise


    !  capillary rise occurs if water content of one horizon is below FC
    !  and water content of the next horizon is larger than FC (field capacity)

    thfree(:,:)=0.

    !** do for all soil components
    DO i=1,nbr_svc(tcid_instance2)

        soilid=id_soil_intern(i,tcid_instance2)
        percol(:)=0.
        vangen(:)=0.
        conduns(:)=0.
        suction(:)=0.


        !**   do for all soil horizons:
        DO h=1,nbrhori(soilid)

            thfree(i,h)=horithact(tcid_instance2,i,h)- soilfc(soilid,h)*  &
                horiz_thickness(tcid_instance2,i,h)

            IF (thfree(i,h) > 0.0) THEN
                !  calculate unsaturated conductivity of horizon
                !  (calculation of vangen: coarse fragment reduction equals out)

                vangen(h)=(horithact(tcid_instance2,i,h)/ horiz_thickness(tcid_instance2,i,h)-  &
                    thetar(soilid,h))/ (thetas(soilid,h)-  &
                    thetar(soilid,h))
                vangen(h)=MAX(vangen(h),0.)
                vangen(h)=MIN(vangen(h),1.)
                conduns(h)=(k_sat(soilid,h)/dt_per_day*(vangen(h)**0.5)*  &
                    (1.-(1.-(vangen(h)** (1./porem(soilid,h)))  &
                    )**porem(soilid,h))**2)* (1.-coarse(soilid,h))
            END IF
        END DO

        !**   do for all soil horizons
        DO h=1,nbrhori(soilid)-1
            IF (thfree(i,h) < 0.) THEN
                IF (thfree(i,h+1) > 0.) THEN

                    exp_x = -1./(thfree(i,h+1)/conduns(h+1))
                    IF (exp_x < -15) THEN !Avoid floating point exception
                        percol(h)=-1.*thfree(i,h+1)
                    ELSE
                        percol(h)=-1.*thfree(i,h+1) * (1.-EXP(exp_x))
                    END IF

                    !   maximum inflow rate is defined by remaining deficit until FC
                    percol(h)=MAX(percol(h),thfree(i,h))

                    !  update soil moisture of both horizons (percol is negative !)
                    horithact(tcid_instance2,i,h+1)=horithact(tcid_instance2,i,h+1)+percol(h)

                    horithact(tcid_instance2,i,h)=horithact(tcid_instance2,i,h)-percol(h)

                END IF
            END IF
        END DO

        !  end of loop for all soil components
    END DO



    !!for debugging - remove
    !DO k=1,nbr_svc(tcid_instance2)
    !	DO h=1,nbrhori(id_soil_intern(k,tcid_instance2))
    !		iF (horithact(tcid_instance2,k,h) - thetas(id_soil_intern(k,tcid_instance2),h)* horiz_thickness(tcid_instance2,k,h)> 0.1 ) THEN	!Till: if water content exceeds saturation...
    !			write(*,*)"nach percol(",n_iter,"): Oversaturated horizon in TC/svc/horizon ",tcid_instance2,k,h
    !			call pause1
    !		END iF
    !	END DO
    !END DO
    !!for debugging - remove end


    !** ------------------------------------------------------------------------
    !** (6) End of calculations for terrain component

    !   new soil moisture of TC


    !only for output
    !thsc(tc_counter2,day,1:nbr_svc(tcid_instance2))=sum(horithact(tcid_instance2,1:nbr_svc(tcid_instance2),:))
    thact=sum(sum(horithact(tcid_instance2,1:nbr_svc(tcid_instance2),:), dim=2)*frac_svc(1:nbr_svc(tcid_instance2),tcid_instance2))



    !   soil moisture in uppermost meter of the profile
    !   exceeding wilting point

    thactroot=0.
    DO i=1,nbr_svc(tcid_instance2)
        soilid=id_soil_intern(i,tcid_instance2)
        temp2=0.
        temp3=0.
        h=1
        DO WHILE (temp2 <= 1000.0 .AND. h <= nbrhori(soilid))
            temp2=temp2+horiz_thickness(tcid_instance2,i,h)
            temp4=horithact(tcid_instance2,i,h)-pwpsc(tcid_instance2,i,h)* horiz_thickness(tcid_instance2,i,h)
            temp4=MAX(temp4,0.)
            IF (temp2 <= 1000.0) THEN
                temp3=temp3+temp4
            ELSE
                temp3=temp3+temp4*(temp2-1000.)/horiz_thickness(tcid_instance2,i,h)
            END IF
            h=h+1
        END DO
        !only for output
        !thsc(tc_counter2,day,i)=temp3
        thactroot=thactroot+temp3*frac_svc(i,tcid_instance2)
    END DO

    !   plant avail. soil moisture in actual rooted soil zone

    !DO i=1,nbr_svc(tcid_instance2)
    !  soilid=id_soil_intern(i,tcid_instance2)
    !  temp2=0.
    !  temp3=0.
    !  temptest=rootd_act(id_veg_intern(i,tcid_instance2))
    !  h=1
    !  DO WHILE (temp2 <= temptest .AND. h <= nbrhori(soilid))
    !    temp2=temp2+horiz_thickness(tcid_instance2,i,h)
    !    temp4=horithact(tcid_instance2,i,h)-pwpsc(tcid_instance2,i,h)* horiz_thickness(tcid_instance2,i,h)
    !    temp4=MAX(temp4,0.)
    !    IF (temp2 <= temptest) THEN
    !      temp3=temp3+temp4
    !    ELSE
    !      temp3=temp3+temp4*(temp2-temptest)/horiz_thickness(tcid_instance2,i,h)
    !    END IF
    !    h=h+1
    !  END DO
    !!only for output
    !!  thsc(tc_counter2,day,i)=temp3
    !END DO


    !    check water balance of timestep
    watbal=watbal-q_surf_out/(tcarea2*1.e3)
    watbal=watbal+(thact1-thact)
    !IF (watbal**2 > 0.001) THEN
    !         write(*,*) 'WATBAL +- delta_theta',watbal,(thact1-thact)
    !         write(*,*) i_subbas2,i_lu,lu_counter2,tcid_instance2,tc_counter2,nbr_svc(tcid_instance2)
    !         write(*,*) 'tcid_instance2,frac_svc',tcid_instance2,frac_svc(:,tcid_instance2)
    !         write(*,*) 'fracsat',frac_sat(tcid_instance2,:)
    !         !call pause1
    !END IF
    !IF (q_surf_out/(tcarea2*1.e3) < -0.001) THEN
    !         write(*,*) 'WATBAL +- delta_theta',watbal,(thact1-thact)
    !         write(*,*) i_subbas2,i_lu,lu_counter2,tcid_instance2,tc_counter2,nbr_svc(tcid_instance2)
    !         write(*,*) 'tcid_instance2,frac_svc',tcid_instance2,frac_svc(:,tcid_instance2)
    !         write(*,*) 'fracsat',frac_sat(tcid_instance2,:)
    !         write(*,*) 'q_surf_out (m**3)',q_surf_out
    !         write(*,*) 'q_surf_out (mm), endI',q_surf_out/(tcarea2*1.e3)
    !         !call pause1
    !END IF
    ! write(*,*)watbal
    ! errstring= GetCalcError()
    ! write(*,*)trim(errstring)
    ! !call pause1


    bal=watbal

    if (dosediment  .AND. .NOT. do_musle_subbasin .AND. (q_surf_out + q_rill_in > 0.)) then !if hillslope erosion is to be computed and there is runoff
        !the following calculation in the peak runoff rate q_peak uses the equations given in the SWAT Theoretical Documentation, pp.105, 2002

        q_rill_out      = q_rill_in
        q_surf_out2     = q_surf_out

        IF (allocated(frac_diff2conc)) then !Till: redistribute sheet/rillflow (temporarily). This is later done again in hymo_all again. Could be done more elegantly only once (ii)
            tempx = q_surf_out * frac_diff2conc(id_tc_type2) !amount of sheetflow that will be channelized
            temp2 = q_rill_in  * frac_conc2diff(id_tc_type2) !amount of concentrated flow that will be redistributed to sheetflow
            q_surf_out2  = q_surf_out2 - tempx + temp2      !sheet flow
            q_rill_out   = q_rill_out  + tempx - temp2      !rill flow
            !beta = q_sheet_out / q_rill_out !interrill/rill-ratio
        END IF

        q_surf_out2     = q_surf_out2 + q_rill_out !qfix Sum up total surface runoff. The distinction is not important any more

        q_surf=q_surf_out2/(tcarea2*1e3)		!surface runoff [mm H2O]

        !ii do this only once for all models at the beginning of model run
        manning_n=0.
        r=0.		!for summing up fractions (actually, these should sum up to 1, but we'll do so just in case)

        DO i=1,size(tc_contains_svc2(id_tc_type2)%p)
            j      =tc_contains_svc2(id_tc_type2)%p(i)%svc_id				!get id of current SVC to be treated
            temp2  =tc_contains_svc2(id_tc_type2)%p(i)%fraction					!get fraction of SVC
            manning_n=manning_n+(svc_n_day(j))*temp2			!average Manning-factors throughout TC, weighted by fraction
            r=r+temp2								!sum up fractions
        end do


        manning_n=manning_n/r


        !debugcheck(i_subbas2,1)=debugcheck(i_subbas2,1)+manning_n
        !debugcheck(i_subbas2,2)=debugcheck(i_subbas2,2)+1

        r=atan(slope(id_tc_type2)/100.)				!convert slope [%] to radiant, avoid multiple computation
        L_slp=slength(i_lu)*fracterrain(id_tc_type2)	 !absolute slope length of TC [m] (projected length, Haan 1994, p.261) (was:.../ cos(r))

        !q_ov=(q_surf_in+q_surf_out)/2./(dt*3600./kfkorr_day)/(tcarea2*1e6/L_slp)		!compute average overland flow rate [m**3/s] on a  1-m-strip (no consieration of rillflow, delete)

        q_ov=(q_surf_in+q_surf_out+q_rill_in+q_rill_out)/2. &
            /(dt*3600./kfkorr_day)/(tcarea2*1e6/L_slp)		!compute average overland flow rate [m**3/s] on a  1-m-strip qfix


        v_ov=(q_ov**0.4)*((slope(id_tc_type2)/100.)**0.3)/manning_n**0.6				!overland flow velocity [m/s] (6.3.4)



        !ii: L_slp pass this to sedi_yield
        t_conc=L_slp/(3600.*v_ov)								!compute time of concentration [h] (6.3.3)

        !compute maximum half_hour-rain: fraction of daily precipitation that falls within 30 min
        if (dt<=1) then	!if high resolution precipitation is available
            if (precday==0) then
                alpha_05=1. !overland flow without rain: may happen in small amounts, thus alpha_05 is not really important
            else
                alpha_05=(maxval(prechall2)/precday)/(dt*2)	!estimate the maximum rainfall rate based on given precipitation data
            end if
        else
            if (kfkorr_a==0.) then
                alpha_05=0.5*kfkorr/24.	!time invariant kfkorr
            else
                !based on kfcorr: time variant kfkorr: only the the precipitation-dependent part of the term is used
                !(the "constant" part of kfkorr factor is considered a calibration factor for Ksat and not
                ! not representing sub-daily rainfall dynamics )
                alpha_05=min(1.0,0.5*(kfkorr_day/kfkorr)/24.)

                !based on USLE Ri50-coefficients - doubles sediment yield in Isábena
                !alpha_05=0.5*a_i30*(precday**(b_i30-1))
            end if
            alpha_05=min(alpha_05,1.0)	!at maximum, all daily rainfall fell within 30 min

        end if

        if (alpha_05==1.0) then
            alpha_t_conc=1.0
        else
            alpha_t_conc=1-exp(2.*t_conc*log(1.-alpha_05))				!(6.3.19) compute fraction of rain falling during the time of concentration
        end if

        q_peak = alpha_t_conc*q_surf*tcarea2 / (3.6*t_conc)			! (6.3.20) estimation of peak runoff rate [m**3/s]


        !old one-particle class call
        !CALL sedi_yield(i_lu, id_tc_type2, q_surf_in, q_surf_out, q_peak, sed_in_tc, dt, tcarea2, r)
        !write(*,'(A,f8.4)') "sed_yield: ", r
        !sed_out_tc(:)=r/n_sed_class	!currently, the particle classes are not treated seperately yet

        CALL sedi_yield(prec, i_subbas2, i_lu, id_tc_type2, tcid_instance2, q_surf_in, q_surf_out2, q_peak, v_ov, sed_in_tc, dt, tcarea2,sed_out_tc)
    else
        sed_out_tc(:)=0.
    end if !end (do sediment)


    RETURN
    END SUBROUTINE soilwat




    logical function isnan(a)
    real :: a

    if (a.ne.a) then
        isnan = .true.
    else
        isnan = .false.
    end if

    return
    end
