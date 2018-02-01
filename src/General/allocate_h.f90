!allocate dynamic arrays
module allocate_h

contains
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
		write(*,'(A,i0,a)')'Memory allocation error (',istate,') in hymo-module: ' 
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

    
END MODULE allocate_h