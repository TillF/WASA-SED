SUBROUTINE sedi_yield_subbas(subbas_id, q_out, sed_yield_subbas)

! Till: computationally irrelevant: minor changes to improve compiler compatibility
! 2011-04-29

! Pedro: computationally relevant correction on q_ov by kfkorrday (overland flow with duration according to rainfall, and not during 24 h)
! 2009-06-03

! Till: cleanup of storage structure for SVCs
! 2008-09-11

! 2008-07-10
! Till: outsourced coefficients for estimation of maximum half-hour rainfall intensity 

! 2008-01-23
! Version with the modifications proposed by Pedro

! 2007-11-13:
! Till: calculate sediment yield of a subbasin based on MUSLE


use hymo_h
use erosion_h
use common_h
use climo_h
use time_h

IMPLICIT NONE


INTEGER, INTENT(IN):: subbas_id		!ID of subbasin currently treated internal numbering scheme)
REAL, INTENT(IN):: q_out			!overland flow leaving the subbasin [m**3]

REAL, INTENT(OUT) :: sed_yield_subbas(1:n_sed_class)		!sediment yield [tons/timestep] (usually applied daily) for each particle size class

REAL :: q		!mean overland flow during timestep [m H2O]
REAL :: frac_svc_x,frac_tc,frac_lu_x				!areal fraction of SVC, TC, LU within the respective higher hierarchy
REAL :: frac_svc_x_sum,frac_tc_sum,frac_lu_x_sum	!sum of areal fraction of SVC, TC, LU within the respective higher hierarchy
REAL :: K_tc, C_tc, P_tc, CFRG_tc, L_tc_cum, S_tc, L_tc, LS_tc, K_lu, C_lu, P_lu, CFRG_lu, LS_lu, K_subbas, C_subbas, P_subbas, CFRG_subbas, LS_subbas		!MUSLE erosion factors on TC, LU and subbasin scale
REAL :: mean_particle_tc(1:n_sed_class),mean_particle_lu(1:n_sed_class),mean_particle_subbas(1:n_sed_class)		!auxiliary variables to compute the mean particle size distribution of current TC, LU and subbasin, respectively
REAL :: manning_n_tc,manning_n_lu,manning_n_subbas	!auxiliary variable to compute manning_n
REAL :: L_slp_subbas			!"averaged" slope length in subbasin [m]
REAL :: slope_lu, slope_subbas	!"averaged" slope in LU and subbasin [%], respectively
INTEGER:: id_tc_type			!ID of TC that is currently treated (internal numbering scheme of TC-types, not instances)
INTEGER:: tc_counter,lu_counter	!counter for treating all TC in current LU, and LUs in current subbasin
INTEGER:: svc_id, lu_id,i			!ID of SVC, LU currently treated
REAL :: r						!temporary real variable

!REAL :: q_peak,q_ov,v_ov,t_conc,alpha_05,alpha_t_conc  		!variables used for calculation of peak runoff rate
REAL :: q_surf					!surface runoff [mm H2O]


q = q_out/(area(subbas_id)*1e6)	!mean overland flow during timestep [m H2O]

!ii: ganze Schleife einmal berechnen für jedes subbasin am Anfang des Programms
slope_subbas =0.
K_subbas     =0.
C_subbas     =0.
P_subbas     =0.
CFRG_subbas  =0.
mean_particle_subbas =0.
manning_n_subbas     =0.
LS_subbas            =0.
L_slp_subbas         =0.
frac_lu_x_sum=0.		!for summing up fractions (actually, these should sum up to 1, but we'll do so just in case)


DO lu_counter=1,nbr_lu(subbas_id)		!loop over all LUs
	lu_id=id_lu_intern(lu_counter,subbas_id)
	frac_lu_x=frac_lu(lu_counter,subbas_id)		!get fraction of LU
	frac_lu_x_sum=frac_lu_x_sum+frac_lu_x		!sum up fractions

	slope_lu=0.
	K_lu=0.
	C_lu=0.
	P_lu=0.
	CFRG_lu=0.
	mean_particle_lu=0.
	manning_n_lu=0.
	LS_lu=0.

	frac_tc_sum=0.		!for summing up fractions (actually, these should sum up to 1, but we'll do so just in case)
	L_tc_cum=0.			!acumulated L factor for erosion calculations

	DO tc_counter=1,nbrterrain(lu_id)		!loop over all TCs
		id_tc_type=id_terrain_intern(tc_counter,lu_id)	!get ID of TC-type
		frac_tc=fracterrain(id_tc_type)					!get fraction of TC
		frac_tc_sum=frac_tc_sum+frac_tc					!sum up fractions

		K_tc=0.
		C_tc=0.
		P_tc=0.
		CFRG_tc=0.
		mean_particle_tc=0.
		manning_n_tc=0.
		frac_svc_x_sum=0.		!for summing up fractions (actually, these should sum up to 1, but we'll do so just in case)

		

		DO i=1,size(tc_contains_svc2(id_tc_type)%p)
			frac_svc_x=    tc_contains_svc2(id_tc_type)%p(i)%fraction					!get fraction of SVC
			frac_svc_x_sum=frac_svc_x_sum+frac_svc_x			!sum up fractions
			svc_id=tc_contains_svc2(id_tc_type)%p(i)%svc_id				!get id of current SVC to be treated
			
			K_tc=K_tc+svc_k_fac_day(svc_id)*frac_svc_x			!erodibility factor that was read from svc.dat with seasonal variability
			C_tc=C_tc+svc_c_fac_day(svc_id)*frac_svc_x			!crop factor that was read from svc.dat with seasonal variability
			P_tc=P_tc+svc_p_fac_day(svc_id)*frac_svc_x			!practice factor that was read from svc.dat
			CFRG_tc=CFRG_tc+exp(-0.03*svc_coarse_fac_day(svc_id))*frac_svc_x	!coarse fragment factor that was read from svc.dat
												!ii: convert musle_fac(*,4) to exp(musle_fac(*,4)) at beginning of program
			mean_particle_tc=mean_particle_tc+soil_particles(svc_soil_veg(svc_id,1),:)*frac_svc_x	!sum up all fractions for particle classes
			manning_n_tc=manning_n_tc+(svc_n_day(svc_id))*frac_svc_x		!average Manning-factors throughout TC, weighted by fraction
		end do



		K_tc=K_tc/frac_svc_x_sum		!normalize summed up values -> final values for TC
		C_tc=C_tc/frac_svc_x_sum
		P_tc=P_tc/frac_svc_x_sum
		CFRG_tc=CFRG_tc/frac_svc_x_sum
		mean_particle_tc=mean_particle_tc/frac_svc_x_sum
		manning_n_tc=manning_n_tc/frac_svc_x_sum

		!Pedro - LS factor computed in a cumulative way according to Haan et al., 1994 (p. 262)
		L_tc_cum=L_tc_cum+slength(lu_id)*frac_tc
		m_ls=0.3*(slope(id_tc_type)/100.)/((slope(id_tc_type)/100.)+exp(-1.47-(61.09*(slope(id_tc_type)/100.))))+0.2	! according to Williams (1995)
		S_tc=(slope(id_tc_type)/100.)*(65.41*(slope(id_tc_type)/100.)+4.56)+0.065
		L_tc=((L_tc_cum**(m_ls+1.))-(max((L_tc_cum-(slength(lu_id)*frac_tc)),0.)**(m_ls+1.)))/(slength(lu_id)*frac_tc*(22.1**m_ls))
		LS_tc=L_tc*S_tc

		slope_lu=slope_lu+slope(id_tc_type)*frac_tc
		K_lu=K_lu+K_tc*frac_tc
		C_lu=C_lu+C_tc*frac_tc
		P_lu=P_lu+P_tc*frac_tc
		CFRG_lu=CFRG_lu+CFRG_tc*frac_tc
		mean_particle_lu=mean_particle_lu+mean_particle_tc*frac_tc
		manning_n_lu=manning_n_lu+manning_n_tc*frac_tc
		LS_lu=LS_lu+LS_tc*frac_tc
	END DO	!loop over all TCs

	slope_lu=slope_lu/frac_tc_sum		!normalize summed up values -> final values for LU
	K_lu=K_lu/frac_tc_sum
	C_lu=C_lu/frac_tc_sum
	P_lu=P_lu/frac_tc_sum
	CFRG_lu=CFRG_lu/frac_tc_sum
	mean_particle_lu=mean_particle_lu/frac_tc_sum
	manning_n_lu=manning_n_lu/frac_tc_sum
	LS_lu=LS_lu/frac_tc_sum

	slope_subbas=slope_subbas+slope_lu*frac_lu_x
	K_subbas=K_subbas+K_lu*frac_lu_x
	C_subbas=C_subbas+C_lu*frac_lu_x
	P_subbas=P_subbas+P_lu*frac_lu_x
	CFRG_subbas=CFRG_subbas+CFRG_lu*frac_lu_x
	mean_particle_subbas=mean_particle_subbas+mean_particle_lu*frac_lu_x
	manning_n_subbas=manning_n_subbas+manning_n_lu*frac_lu_x
	LS_subbas=LS_subbas+LS_lu*frac_lu_x
	L_slp_subbas=L_slp_subbas+slength(lu_id)*frac_lu_x
END DO	!loop over all LUs

slope_subbas=slope_subbas/frac_lu_x_sum			!normalize summed up values -> final values for subbasin
K_subbas=K_subbas/frac_lu_x_sum
C_subbas=C_subbas/frac_lu_x_sum
P_subbas=P_subbas/frac_lu_x_sum
CFRG_subbas=CFRG_subbas/frac_lu_x_sum
mean_particle_subbas=mean_particle_subbas/frac_lu_x_sum
manning_n_subbas=manning_n_subbas/frac_lu_x_sum
LS_subbas=LS_subbas/frac_lu_x_sum
L_slp_subbas=L_slp_subbas/frac_lu_x_sum


K_fac=K_subbas	!Till: this is done just to match the style below which is also used in sedi_yield.f90
C_fac=C_subbas
P_fac=P_subbas
CFRG_fac=CFRG_subbas
manning_n=manning_n_subbas
LS_fac=LS_subbas
L_slp=L_slp_subbas

q_ov=q_out/(dt*3600./kfkorr_day)/(area(subbas_id)*1e6/L_slp_subbas)				!compute average overland flow rate [m**3/s] on a 1-m-strip
v_ov=(q_ov**0.4)*((slope_subbas/100.)**0.3)/manning_n_subbas**0.6	!overland flow velocity [m/s] (6.3.4)
t_conc=L_slp_subbas/(3600.*v_ov)										!compute time of concentration [h] (6.3.3)

earea=area(subbas_id)*1e6	!area of erosive unit [m2]


!compute maximum half-hour rain: fraction of daily precipitation that falls within 30 min (needed for estimation of peak flow and rainfall erosivity)
!IF (dt<=1) THEN		!if high resolution precipitation is available
!	alpha_05=(maxval(preciph((d-1)*24+1:d*24,subbas_id))/precip(d,subbas_id))/(dt*2)	!estimate the maximum rainfall rate based on given precipitation data
!ELSE
!	IF (kfkorr_a==0) THEN
!		alpha_05=0.5*kfkorr/24	!time invariant kfkorr
!	ELSE
!		alpha_05=min(1.0,0.5*(kfkorr_day/kfkorr)/24)	!time variant kfkorr: only the the precipitation-dependent part of the term is used (the kfkorr factor is considered a calibration factor for Ksat and not responsible for sub-daily rainfall dynamics)
!	END IF
!END IF


IF ((erosion_equation==1) .OR. (erosion_equation==2))  THEN
	!compute USLE-rainfall energy factor for USLE (1) and Onstad-Foster (2)
	IF (dt<=1) THEN	!if high resolution precipitation is available
		R_d=-1.		!to do: add equation here for subdaily rainfall
		r_p=-1.
	ELSE
		R_d=precip(d,subbas_id)	!daily rainfall [mm]
		!r_p=-2*R_d*log(1-min(alpha_05,0.99))	!peak rainfall rate [mm/h] Williams, 1995; eq. 25.132
		!R_05=alpha_05*R_d						!maximum amount of rain in 30 min [mm]; Williams, 1995; 25.131
		!ri_05=R_05/0.5							!maximum 0.5-h rainfall intensity [mm/h]

		!rainfall intensities based on kfkorr showed low values, resulting in low erosion as well
		ri_05=a_i30*(R_d**b_i30)
		r_p=-2.*R_d*log(1-min((ri_05/2./R_d),0.99))

	END IF

	ei=R_d*(12.1+8.9*(log10(r_p)-0.434))*ri_05/1000.	!USLE-energy factor in the "proper" units according to Williams, 1995 in Singh,1995, p.934,25.128
ELSE
	ei=-1.		!just for debugging
END IF


IF ((erosion_equation==2) .OR. (erosion_equation==3) .OR. (erosion_equation==4)) THEN
R_d=precip(d,subbas_id)	!daily rainfall [mm]
ri_05=a_i30*(R_d**b_i30)
!alpha_05=ri_05/2/R_d
alpha_05=min(1.0,0.5*(kfkorr_day/kfkorr)/24.)

	IF (alpha_05==1.0) THEN
		alpha_t_conc=1.0
	ELSE
		alpha_t_conc=1.-exp(2.*t_conc*log(1.-min(alpha_05,0.99)))		!(6.3.19) compute fraction of rain falling during the time of concentration
	END IF

	!compute surface runoff and peak flow for Onstad-Foster (2), MUSLE (3), MUST (4)
	q_surf=(q_out/earea)*1000.		!surface runoff [mm H2O]
	q_peak=alpha_t_conc*q_surf*earea/1e6/(3.6*t_conc)	!(6.3.20) estimation of peak runoff rate [m**3/s]
	q_peak_mm_h=q_peak*3.6/(earea/1e6)		!peak runoff as before in [mm/h]
ELSE
	q_peak=-1.	!just for debugging
	q_peak_mm_h=-1.
END IF


SELECT CASE (erosion_equation)
	CASE (1)	!USLE
		r = earea*1e-4 * ei*(K_fac*C_fac*P_fac*LS_fac*CFRG_fac)
		!yield [t] according to USLE as given by Williams,1995 in Singh, 1995, p. 933
	CASE (2)	!Onstad-Foster
		r = earea*1e-4 * (0.646*ei+0.45*((q*1000.)*q_peak_mm_h)**0.33)*(K_fac*C_fac*P_fac*LS_fac*CFRG_fac)
		!yield [t] according to Onstad-Forster as given by Williams,1995 in Singh, 1995, p. 933
	CASE (3)	!MUSLE
		!ii area kürzt sich raus mit folgender Gleichung
		r = 11.8*(q*q_peak*earea)**0.56*(K_fac*C_fac*P_fac*LS_fac*CFRG_fac)
		!yield [t] according to MUSLE as given by Summer&Walling (equivalent to Williams,1995) !different units as in SWAT-manual!
	CASE (4)	!MUST
		r = earea*1e-4 * 2.5*((q*1000.)*q_peak_mm_h)**0.5*(K_fac*C_fac*P_fac*LS_fac*CFRG_fac)
		!yield [t] according to MUST as given by Williams,1995 in Singh, 1995, p. 933
	CASE DEFAULT
	   WRITE (*,*) "Erosion equation ",erosion_equation," not defined, please check manual"
END SELECT

sed_yield_subbas = r*mean_particle_subbas		!overall yield is distributed among size classes according to their percentage, no selective removal of finer fractions

return
end SUBROUTINE sedi_yield_subbas
