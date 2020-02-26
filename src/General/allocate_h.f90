!allocate dynamic arrays
module allocate_h

    contains

subroutine allocate_general()

    use params_h
    use climo_h
    use common_h

    implicit none

    INTEGER :: istate
   	allocate( &
     temp(366,subasin), & !ab r
	 precip(366,subasin), & !ab r
	 preciph(366*24,subasin), & !abc r
	 pet(366,subasin), & !ab r
	 rad(366,subasin), & !ab r
	 rhum(366,subasin), & !ab r
	 wind(366,subasin), & !ab r
	 do_pre_outflow(subasin), & !a r
      STAT = istate)

	if (istate/=0) then
		write(*,'(A,i0,a)')'ERROR: Memory allocation error (',istate,') in general-module: '
		stop
	end if
    do_pre_outflow=.true.

    contains
    subroutine mem_usage() !estimate required memory and identify most important saving options
        integer(kind=4) :: total_required !kbytes required
        integer(kind=4) a, ab, abc !auxiliary variables for computing similarly structured variables

        if (.not. do_mem_diagnostics) return
        a = subasin
        ab = a * 366
        abc = ab * 24
        total_required = nint(  (1.0 * a + 6 * ab + 1 * abc) * SIZEOF(1.1) / 1024) !size in kBytes
        write(*,*) 'Memory requirements general [kB]: ', total_required
        write(*,*) ' influence rates:'
        write(*,*) '   number of subbasins : ', nint ((1.0 + 6*ab/a + 1*abc/a)* SIZEOF(1.1) /a )

    end subroutine mem_usage

end subroutine allocate_general


subroutine allocate_hymo()

use hymo_h
use params_h

implicit none

INTEGER :: istate

	 allocate( &
     area(subasin), &
	 id_subbas_extern(subasin), &
	 id_subbas_intern(subasin), &
	 nbr_lu(subasin), &
	 id_lu_intern(maxsoter,subasin), &
	 tcallid(subasin, maxsoter, maxterrain), &
	 frac_lu(maxsoter,subasin), &
	 kfkorrc(subasin), &
	 intcfc(subasin), &
	 !period(4,subasin,tstop-tstart+1), &
	 laimun(366,subasin), &
	 soilm(366,subasin), &
	 soilmroot(366,subasin), &
	 aet(366,subasin), &
	 soilet(366,subasin), &
	 intc(366,subasin), &
	 ovflow(366,subasin), &
	 hortflow(366,subasin), &
	 subflow(366,subasin), &
	 gw_recharge(366,subasin), &
	 deepgw(subasin,maxsoter), &
	 !deepgw_r(366,subasin), &
	 qgen(366,subasin), &
	 water_subbasin_t(366,nt,subasin), &
	 water_subbasin(366,subasin), &
	 qloss(366,subasin), &
	 sofarea(366,subasin), &
	 avail_all(366,subasin), &
	 avail_ac(366,subasin), &

	 nbr_svc(ntcinst), &

	 rocky(ntcinst), &
	 laitc(ntcinst), &
	 aettc(ntcinst), &
	 soilettc(ntcinst), &
	 intctc(ntcinst), &
	 horttc(ntcinst), &
	 gwrtc(ntcinst), &
	 deepgwrtc(ntcinst), &

	 id_soil_intern(maxsoil,ntcinst), &
	 id_veg_intern(maxsoil,ntcinst), &
	 svcrooth(ntcinst,maxsoil), &
	 svcbedr(ntcinst,maxsoil), &

	 frac_svc(maxsoil,ntcinst), &
	 soilwater(366,ntcinst), &

	 intercept(ntcinst,maxsoil), &
	 frac_sat(ntcinst,maxsoil), &
	 thsprof(ntcinst,maxsoil), &

	 horiz_thickness(ntcinst,maxsoil,maxhori), &
	 pwpsc(ntcinst,maxsoil,maxhori), &
	 horiths(ntcinst,maxsoil,maxhori), &
	 horithact(ntcinst,maxsoil,maxhori), &

	 tctheta_s(5,2,ntcinst,maxsoil), & !allocated in soil_distr

!	 latred(ntcinst,maxhori*3), &		!Till: latred is allocated in readhymo.f90 to adjust to maximum number of exchange horizons required

	 id_lu_extern(nsoter), &
	 nbrterrain(nsoter), &
	 id_terrain_intern(maxterrain,nsoter), &
	 gw_flag(nsoter), &
	 kfsu(nsoter), &
	 slength(nsoter), &
	 meandep(nsoter), &
	 maxdep(nsoter), &
	 riverbed(nsoter), &
	 gw_dist(nsoter), &
	 gw_delay(nsoter), &
	 orderterrain(maxterrain,nsoter), &

	 id_terrain_extern(nterrain), &
	 slope(nterrain), &
	 fracterrain(nterrain), &
	 posterrain(nterrain), &

	 id_soil_extern(nsoil), &
	 nbrhori(nsoil), &
	 bedrock(nsoil), &
	 alluvial_flag(nsoil), &

	 thetar(nsoil,maxhori), &
	 soilpwp(nsoil,maxhori), &
	 thetas(nsoil,maxhori), &
	 soilnfk(nsoil,maxhori), &
	 soilfc(nsoil,maxhori), &
	 soilfc63(nsoil,maxhori), &
	 k_sat(nsoil,maxhori), &
	 saug(nsoil,maxhori), &
	 poresz(nsoil,maxhori), &
	 porem(nsoil,maxhori), &
	 bubble(nsoil,maxhori), &
	 coarse(nsoil,maxhori), &
	 shrink(nsoil,maxhori), &

	 id_veg_extern(nveg), &
	 height(nveg,4), &
	 rootdep(nveg,4), &
	 lai(nveg,4), &
	 alb(nveg,4), &
	 resist(nveg), &
	 wstressmin(nveg), &
	 wstressmax(nveg), &

 	 intcept_mem(maxterrain,maxsoil,24), &
	 aet_red_mem(maxterrain,maxsoil), &

	 laisu(maxsoter), &
	 soilmsu(maxsoter), &
	 soilmrootsu(maxsoter), &
	 aetsu(maxsoter), &
	 soiletsu(maxsoter), &
	 intcsu(maxsoter), &
	 qsurf_lu(maxsoter), &
	 qsub_lu(maxsoter), &
	 gwrsu(maxsoter), &
	 deepgwrsu(maxsoter), &
	 hortsu(maxsoter), &
!	 resistsc(366,maxsoil), &



!	 sat_xs_overland_flow_sc(366,maxsoil), &
!     sat_area_of_sc(366,maxsoil), &
!     hortsc(366,maxsoil), &
!     hort2sc(366,maxsoil), &
!     gwrsc(366,maxsoil), &
!     deepgwrsc(366,maxsoil), &
!     allocate(thsc(3,366,maxsoil), &
!     nfksc(366,maxsoil), &
!     aetsc(366,maxsoil), &
!     soiletsc(366,maxsoil), &
     aet1sc(366,maxsoil), &
     soilet1sc(366,maxsoil), &
!     intetsc(366,maxsoil), &
!     aettotsc(366,maxsoil), &

!	 debug_out(366), & !remove
!	 debug_out2(366,8), & !remove

	 STAT = istate)

	if (istate/=0) then
		write(*,'(A,i0,a)')'ERROR: Memory allocation error (',istate,') in hymo-module: '
		stop
	end if

	aet_t               => allocate_hourly_array(f_actetranspiration)
	hortflow_t          => allocate_hourly_array(f_qhorton)
	deep_gw_discharge_t => allocate_hourly_array(f_gw_discharge)
	pet_t               => allocate_hourly_array(f_potetranspiration)
	deep_gw_recharge_t  => allocate_hourly_array(f_gw_recharge)
	gw_loss_t           => allocate_hourly_array(f_gw_loss)
	river_infiltration_t=> allocate_hourly_array(f_river_infiltration)
	riverflow_t         => allocate_hourly_array(f_river_flow)

    subflow_t          => allocate_hourly_array(f_subsurface_runoff .OR. .TRUE.)	!Till: currently, these are needed in any case
	ovflow_t           => allocate_hourly_array(f_total_overlandflow .OR. .TRUE.)

END subroutine allocate_hymo


subroutine allocate_routing()

use hymo_h
use params_h
use routing_h

implicit none

INTEGER :: istate


!ii: do allocation specific to routing mode to save memory
	 allocate ( &
     upbasin(subasin), &
	 downbasin(subasin), &
!	 qin(372,subasin), &
! 	 qout(372,subasin), &
	 prout (subasin,2), &
	 r_qin(2,subasin), &
	 r_qout(2,subasin), &
	 r_storage(subasin), &
!	 qsediment(372,subasin), &

!	 runoff(8784,subasin), &

	 r_width(subasin), &
     r_depth(subasin), &
	 r_depth_cur(subasin), &
	 r_sideratio(subasin), &
	 r_width_fp(subasin), &
	 r_sideratio_fp(subasin), &
     r_slope(subasin), &
     r_length(subasin), &
     manning(subasin), &
	 manning_fp(subasin), &
     r_ksat(subasin), &
     r_efactor(subasin), &
     r_cover(subasin), &
	 r_rock(subasin), &
     r_alpha(subasin), &
	 msk_x(subasin), &
	 msk_k(subasin), &
	 Q_spring(subasin), &

	 area_bankful(subasin), &
	 area_loflo(subasin), &
	 q_bankful100(subasin), &
	 bottom_width(subasin), &


	 phi5(subasin), &
	 phi10(subasin), &
	 phi13(subasin), &
	 velocity(subasin), &


	 sed_storage(subasin,n_sed_class), &
	 sediment_in(subasin,n_sed_class), &
	 sediment_out(subasin,n_sed_class), &
	 qsediment2_t(366,nt,subasin), &
	 r_sediment_concentration(subasin), &
!	 sediment(8784,subasin,n_sed_class), & !use dynamic margins!
	 river_deposition(subasin,n_sed_class), &
	 river_degradation(subasin,n_sed_class), &
	 riverbed_storage(subasin,n_sed_class), &

	 D50(subasin), &
	 bedload(subasin,5), &

	 STAT = istate)

	if (istate/=0) then
		write(*,'(A,i0,a)')'Memory allocation error (',istate,') in river-module.'
	end if

	if (f_river_sediment_storage) allocate (river_sediment_storage_t(366,nt,subasin),  STAT = istate)
	if (istate/=0) write(*,'(A,i0,a)')'Memory allocation error (',istate,') in river-module.'

	if (f_river_susp_sediment_storage) allocate (river_susp_sediment_storage_t(366,nt,subasin),  STAT = istate)
	if (istate/=0) write(*,'(A,i0,a)')'Memory allocation error (',istate,') in river-module.'

    if (ntrans>0) then
        allocate(trans_start(2, ntrans), STAT = istate)
        allocate(trans_end  (2, ntrans), STAT = istate)
        if (istate/=0) write(*,'(A,i0,a)')'Memory allocation error (',istate,') in river-module.'
    end if    

END subroutine allocate_routing

subroutine allocate_reservoir()

    use hymo_h
    use params_h
    use lake_h
    use time_h
    use reservoir_h

    implicit none

    INTEGER :: istate

    IF (doacud) THEN
	    allocate( &
	    acud(subasin,5), &
	    acudfloat(subasin,5), &
	    acudfloatyear(subasin,5,tstop-tstart+1), &
	      maxlakesub0(subasin,5), &
	    !  lakewater0(subasin,5), &
	      laketrend(subasin,5), &
	      maxlakewater(subasin,5), &
	      lakewater(366*nt,subasin,5), &
	      laketot(366*nt,subasin), &
	      laketotfrac(366*nt,subasin), &
	      lakeevap(366*nt,subasin,5), &
	      lakeprec(366*nt,subasin,5), &
	      lakeex(366*nt,subasin,5), &
	      lakearea(subasin,5), &
	      maxlake_factor(subasin,5), &
	      lakeinflow(366*nt,subasin), &
	      lakeoutflow(366*nt,subasin), &
	      muniret(366*nt,subasin), &
          maxlake(subasin,5), &
	      maxstorcap_hrr(366*nt,subasin,5), &
	      lakeinflow_hrr(366*nt,subasin,5), &
	      lakeoutflow_hrr(366*nt,subasin,5), &
	      lakeretention_hrr(366*nt,subasin,5), &
	      lakewater_hrr(366*nt,subasin,5), &
          outflow_last_hrr(subasin,5), &
          hmax_hrr(subasin,5), &
          volume_last_hrr(subasin,5), &
          lakefrarea(subasin,5), &
          subfrarea(subasin), &
          subfrout(subasin), &
          lakefrout(subasin,5), &
          lakecumfrout(subasin,5), &
          lakerunoff(366*nt,subasin), &

	    STAT = istate)

	    if (istate/=0) then
		    write(*,'(A,i0,a)')'ERROR: Memory allocation error (',istate,') in acude-module (small reservoirs).'
		    stop
	    end if


	    if (dosediment) then
		    allocate( &
		      lakesedin(366*nt,subasin), &
		      lakesedout(366*nt,subasin), &
		      lakefracsedin(366*nt,subasin,n_sed_class), &
		      lakefracsedout(366*nt,subasin,n_sed_class), &
		      lakesedin_hrr(366*nt,subasin,5), &
		      lakesedout_hrr(366*nt,subasin,5), &
		      lakefracsedin_hrr(366*nt,subasin,5,n_sed_class), &
		      lakefracsedout_hrr(366*nt,subasin,5,n_sed_class), &
		      lakedep_hrr(366*nt,subasin,5), &
		      lakecumdep_hrr(subasin,5), &
		      lake_vollost(366*nt,subasin,5), &
		      cumseddep(subasin), &
		      lake_cumseddep(366*nt,subasin,5), &

		    STAT = istate)
		    if (istate/=0) then
			    write(*,'(A,i0,a)')'Memory allocation error (',istate,') in acude-module (small reservoirs, sediment).'
			    stop
		    end if
	    end if

    END IF




    IF (doreservoir) THEN
	    allocate( &
            
 !Anne & Till 2019 fix reservoir memory issue:
        !to decrease array size & only do calculations for subbasins with reservoir,  
        !moved all arrays with "subbasin" to reservoir.f90 and substituted by "n_reservoir";
        !plus inserted "res_index()", e.g.: dayarea_bat(step,j,i) changed to dayarea_bat(step,j,res_index(i))   
  
  !The following variables could not be moved/changed, because they are allocated after READ-operation from file in reservoir.f90
  !	        minlevel(i), maxlevel(i),vol0(i),storcap(i), &
  !			damflow(i),damq_frac(i),withdrawal(i),damyear(i),maxdamarea(i), &
  !			damdead(i),damalert(i),dama(i),damb(i),qoutlet(i),fvol_bottom(i), &
  !			fvol_over(i),damc(i),damd(i),elevbottom(i)
  !         forma_factor(i)
  !         damareaact(i), volact(id,i)
  !         volume_last(i), outflow_last(i)   !for lake.f90       
                                
	      res_index(subasin), &
          f_intake_obs(subasin), &      !flag
          res_flag(subasin), &          !flag
   	      fcav(subasin), &              !flag
	      fvol_over(subasin), &         !flag
          latflow_res(subasin), &       !flag
 	      damyear(subasin), &
          volume_last(subasin), &
          outflow_last(subasin), &
!	      damex(366*nt,subasin), & tp not used
	      storcap(subasin), &
	      damdead(subasin), &
	      damalert(subasin), &
	      elevbottom(subasin), &
	      damflow(subasin), &
	      damq_frac(subasin), &
          volact(366*nt,subasin), &
            
! tp TODO never used
!	      evapdam(366*nt,subasin), &
!	      infdam(366*nt,subasin), &
	      maxdamarea(subasin), &
	      damareaact(subasin), &
	      dama(subasin), &
	      damb(subasin), &
	      damc(subasin), &
	      damd(subasin), &
	      forma_factor(subasin), &
	      qoutlet(subasin), &
	      fvol_bottom(subasin), &
	      withdrawal(subasin), &

    !	  elev_bat0(maxreservoir_b,subasin), &
    !	  area_bat0(maxreservoir_b,subasin), &
    !	  vol_bat0(maxreservoir_b,subasin), &
    !	  elev_bat(maxreservoir_b,subasin), &
    !	  area_bat(maxreservoir_b,subasin), &
    !	  vol_bat(maxreservoir_b,subasin), &
	      vol0(subasin), &
	      maxlevel(subasin), &
	      minlevel(subasin), &
            
	      geom(npointsxsect,nxsection_res), &           
    !	  sed_susp(366*nt,subasin), &           
	      diam(n_sed_class), &
    
    !	  min_conc(subasin), &
    !	  sed_conc0(subasin), &
    !	  wet_dens(subasin), &                      
    !	  numbeqs(subasin), &
    !	  Q_max(10,subasin), &
    !	  param_a(10,subasin), &
    !	  param_b(10,subasin), &
    !	  nbsizedist(subasin), &
    !	  Q_refer(10,subasin), &          	          
	      setvel(n_sed_class), &
                        
	      frsedavailab(n_sed_class,nxsection_res), &
	      frerosion(n_sed_class,nxsection_res), &
	      frdeposition(n_sed_class,nxsection_res), &
	      frretention(n_sed_class,nxsection_res), &           
    !	  frsuspension(n_sed_class,nxsection_res), &
    !	  frbed_discharge(n_sed_class,nxsection_res), &
    !	  frsusp_discharge(n_sed_class,nxsection_res), &           
	      frtotal_discharge(n_sed_class,nxsection_res), &        
    !	  suspension(nxsection_res,subasin), &
            
	      bed_frtransp(n_sed_class,nxsection_res), &
	      susp_frtransp(n_sed_class,nxsection_res), &
	      fr_capacity(n_sed_class,nxsection_res), &            
    !	  dheight_sed(nxsection_res,subasin), &
           
	      frconc(n_sed_class,nxsection_res), &                        
    !	  weightfac_actlay(npointsxsect,nxsection_res,subasin), &
    !	  weightfac_toplay(npointsxsect,nxsection_res,subasin), &
	      
	     STAT = istate)

	    if (istate/=0) then
		    write(*,'(A,i0,a)')'ERROR: Memory allocation error (',istate,') in reservoir-module.'
		    stop
	    end if


!Anne moved reservoir sediment variables from allocate_h to reservoir.f90
!A	    if (dosediment) then
!A		    allocate( &

!A		      daydepth_sec(366*nt,nxsection_res,subasin), &
!A		      daywatelev_sec(366*nt,nxsection_res,subasin), &
!A		      dayarea_sec(366*nt,nxsection_res,subasin), &
!A		      daytopwidth_sec(366*nt,nxsection_res,subasin), &
!A		      dayenergslope_sec(366*nt,nxsection_res,subasin), &
!A		      dayhydrad_sec(366*nt,nxsection_res,subasin), &
!A		      daymeanvel_sec(366*nt,nxsection_res,subasin), &
!A		      daydischarge_sec(366*nt,nxsection_res,subasin), &
!A		      dayminelev_sec(366*nt,nxsection_res,subasin), &
!A		      dayy_sec(366*nt,npointsxsect,nxsection_res,subasin), &
!A		      daycumsed(366*nt,subasin), &
!A		      dayfrsediment_out(366*nt,subasin,n_sed_class), &
!A		    STAT = istate)

!A		    if (istate/=0) then
!A			    write(*,'(A,i0,a)')'ERROR: Memory allocation error (',istate,') in reservoir-module (sediments).'
!A			    stop
!A		    end if

!A	    end if



    END IF

end subroutine allocate_reservoir

!subroutine allocate_erosion()
!
! ! variables are allocated in the end of readhymo, when all dimensions are known
!
!end subroutine allocate_erosion

subroutine allocate_snow()
!Allocate variables related to the snow module

use snow_h
use params_h

implicit none

integer :: idummy, istate, i_lu, tc_counter
real :: temp2

    idummy = 0 !for summing up IO-error flags

    allocate(snowEnergyCont(366,nt,ntcinst), &
             snowWaterEquiv(366,nt,ntcinst), &
             snowAlbedo    (366,nt,ntcinst), &
             snowCover     (366,nt,ntcinst), &
             precipMod     (366,nt,ntcinst), &
             STAT = istate)
    idummy = istate + idummy

    !Initial values states
    snowEnergyCont(1:366,1:nt,:)     =     0.
    snowWaterEquiv(1:366,1:nt,:)     =     0.
    snowAlbedo(1:366,1:nt,:)         =     albedoMax

    allocate(rel_elevation(nterrain), STAT = istate)
    idummy = istate + idummy

    if (f_radiMod) then
       allocate(radiMod(366,nt,ntcinst), STAT = istate)
    else
       allocate(radiMod(1,1,1), STAT = istate)
    end if
    idummy = istate + idummy

    if (f_temperaMod) then
       allocate(temperaMod(366,nt,ntcinst), STAT = istate)
    else
       allocate(temperaMod(1,1,1), STAT = istate)
    end if
    idummy = istate + idummy

    if (f_cloudFrac) then
       allocate(cloudFrac(366,nt,ntcinst), STAT = istate)
    else
       allocate(cloudFrac(1,1,1), STAT = istate)
    end if
    idummy = istate + idummy

    if(f_snowTemp) then
       allocate(snowTemp(366,nt,ntcinst), STAT = istate)
    else
       allocate(snowTemp(1,1,1), STAT = istate)
    end if
    idummy = istate + idummy

    if (f_surfTemp) then
       allocate(surfTemp(366,nt,ntcinst), STAT = istate)
    else
       allocate(surfTemp(1,1,1), STAT = istate)
    end if
    idummy = istate + idummy

    if (f_liquFrac) then
       allocate(liquFrac(366,nt,ntcinst), STAT = istate)
    else
       allocate(liquFrac(1,1,1), STAT = istate)
    end if
    idummy = istate + idummy

    if (f_fluxPrec) then
       allocate(fluxPrec(366,nt,ntcinst), STAT = istate)
    else
       allocate(fluxPrec(1,1,1), STAT = istate)
    end if
    idummy = istate + idummy

    if (f_fluxSubl) then
       allocate(fluxSubl(366,nt,ntcinst), STAT = istate)
    else
       allocate(fluxSubl(1,1,1), STAT = istate)
    end if
    idummy = istate + idummy

    if (f_fluxFlow) then
       allocate(fluxFlow(366,nt,ntcinst), STAT = istate)
    else
       allocate(fluxFlow(1,1,1), STAT = istate)
    end if
    idummy = istate + idummy

    if (f_fluxNetS) then
       allocate(fluxNetS(366,nt,ntcinst), STAT = istate)
    else
       allocate(fluxNetS(1,1,1), STAT = istate)
    end if
    idummy = istate + idummy

    if (f_fluxNetL) then
       allocate(fluxNetL(366,nt,ntcinst), STAT = istate)
    else
       allocate(fluxNetL(1,1,1), STAT = istate)
    end if
    idummy = istate + idummy

    if (f_fluxSoil) then
       allocate(fluxSoil(366,nt,ntcinst), STAT = istate)
    else
       allocate(fluxSoil(1,1,1), STAT = istate)
    end if
    idummy = istate + idummy

    if (f_fluxSens) then
       allocate(fluxSens(366,nt,ntcinst), STAT = istate)
    else
       allocate(fluxSens(1,1,1), STAT = istate)
    end if
    idummy = istate + idummy

    if (f_stoiPrec) then
       allocate(stoiPrec(366,nt,ntcinst), STAT = istate)
    else
       allocate(stoiPrec(1,1,1), STAT = istate)
    end if
    idummy = istate + idummy

    if (f_stoiSubl) then
       allocate(stoiSubl(366,nt,ntcinst), STAT = istate)
    else
       allocate(stoiSubl(1,1,1), STAT = istate)
    end if
    idummy = istate + idummy

    if (f_stoiFlow) then
       allocate(stoiFlow(366,nt,ntcinst), STAT = istate)
    else
       allocate(stoiFlow(1,1,1), STAT = istate)
    end if
    idummy = istate + idummy

    if (f_rateAlbe) then
       allocate(rateAlbe(366,nt,ntcinst), STAT = istate)
    else
       allocate(rateAlbe(1,1,1), STAT = istate)
    end if
    idummy = istate + idummy

    if (do_rad_corr .OR. do_alt_corr) then
        allocate(lu_aspect(nsoter), STAT = istate)
        allocate(lu_alt(nsoter), STAT = istate)
    end if
    idummy = istate + idummy

    if (idummy/=0) then
        write(*,'(A,i0,i0,a)')'ERROR: Memory allocation error (',istate,'/',idummy,') in hymo-module (snow): '
        stop
    end if



    !Computation relative elevation of TC centre above foot of toposequence/LU (i.e. river) [m]
    DO i_lu=1,nsoter
       temp2 = 0. !cumulative elevation of downslope TCs
       DO tc_counter=1,nbrterrain(i_lu)
           !Relative elevation of center TC
           !Note: slope(tc_counter) in %
           rel_elevation(id_terrain_intern(tc_counter,i_lu)) = temp2 + slength(i_lu)*fracterrain(tc_counter)*slope(tc_counter)/100/2
           temp2 = temp2 + slength(i_lu)*fracterrain(tc_counter)*slope(tc_counter)/100 !cumulate y-extent TCs
       END DO

       DO tc_counter=1,nbrterrain(i_lu)
           rel_elevation(id_terrain_intern(tc_counter,i_lu)) = rel_elevation(id_terrain_intern(tc_counter,i_lu)) - (temp2 / 2.) !"normalize" to middle elevation of the LU
           !Relative elevation of LU not available at this point; read in later; added in readhymo after read in of lu2.dat
       END DO

    END DO


end subroutine allocate_snow


end module allocate_h
