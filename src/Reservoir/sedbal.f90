SUBROUTINE sedbal(upstream)

!Till: computationally irrelevant: outcommented unused vars
!2012-09-14

! Till: computationally irrelevant: minor changes to improve compiler compatibility
! 2011-04-29

! Code converted using TO_F90 by Alan Miller
! Date: 2005-08-23  Time: 12:57:31

use common_h
use time_h
use reservoir_h
use routing_h
use hymo_h

IMPLICIT NONE

INTEGER, INTENT(IN OUT)                  :: upstream
INTEGER :: k,g,l,nbrsteps !i,
real:: dummy

real :: wat_dens,sed_dens,delta,grav,visc
real :: overflow_rate,diam_equiv,frac_finer,trap_eff
real :: frretent,frrelease !frcum,
real :: cumfrsed_in(n_sed_class),cumfrsed_out(n_sed_class),lower_limit(n_sed_class)
real :: par_a,par_b,par_c
real :: interv(200),setveloc(200),upper_diam(200),lower_diam(200),mean_diam(200),frfiner(200)
real :: frfiner_out(200),cumfrrelease(200)
real :: par_x(200),par_y(200) !,par_z(200)

!Ge to include the DAILY mean temperature of the reservoir (celsius degree)
!Ge tempres=20 C (temporarily)
REAL :: tempres(step,upstream)

! -----------------------------------------------------------------------

if (res_qout(step,res_index(upstream))/=0.) then
!Ge to include the DAILY mean temperature of the reservoir (celsius degree)
!Ge tempres=20 C (temporarily)
  tempres(step,upstream)=20.

! density of water and density of natural sediments (kg/m3)
  wat_dens=1.*1000.
  sed_dens=2.65*1000.

! specific gravity of sediments in water (-)
  delta=(sed_dens-wat_dens)/wat_dens

! gravitacional acceleration (m2/s)
  grav=9.807

! kinematic viscosity (m2/s)
  visc=(1.14-0.031*(tempres(step,upstream)-15.)+  &
      0.00068*((tempres(step,upstream)-15.)**2.))*1.e-6

! Overflow rate  Vc (m/s)
  overflow_rate=res_qout(step,res_index(upstream))/damareaact(upstream)
!  overflow_rate=0.00032

! calculation of the corresponding equivalent diameter (mm)
  par_a=1.09*(delta)*grav
  par_b=-(overflow_rate**2)
  par_c=-27.9*overflow_rate*visc

  diam_equiv=1000.*(-par_b+sqrt(((par_b)**2.)-4.*par_a*par_c))/(2.*par_a)
!if (t==2004 .and. step<100) write(*,'(2I4,10F15.6)') step,id_subbas_extern(upstream),res_qout(step,upstream),damareaact(upstream),overflow_rate
!if (t==2004 .and. step<100) write(*,'(2I4,10F15.6)') step,id_subbas_extern(upstream),par_a,par_b,par_c,diam_equiv
!if (t==2004 .and. step<100) pause
! affluent grain size distribution (-)
  dummy=0.
  DO g=1,n_sed_class
    dummy=dummy+frsediment_in(res_index(upstream),g)
    cumfrsed_in(g)=dummy
!write(*,*)g,cumfrsed_in(g)
  enddo

! lower limits of particle size class intervalls (mm)
  lower_limit(1)=upper_limit(1)/2.
  DO g=2,n_sed_class
    lower_limit(g)=upper_limit(g-1)
  enddo
!write(*,*)g,(lower_limit(g),g=1,n_sed_class)

! fraction of particles smaller than equivalent diameter (-)
  DO g=1,n_sed_class
    if (diam_equiv>=lower_limit(g) .and. diam_equiv<=upper_limit(g) .and. g>1) then
      frac_finer=cumfrsed_in(g-1)+(cumfrsed_in(g)-cumfrsed_in(g-1))* &
		(log10(diam_equiv)-log10(lower_limit(g)))/ &
		(log10(upper_limit(g))-log10(lower_limit(g)))
	  exit
	else if (diam_equiv>=lower_limit(1) .and. diam_equiv<=upper_limit(1)) then
      frac_finer=cumfrsed_in(g)* &
		(log10(diam_equiv)-log10(lower_limit(g)))/ &
		(log10(upper_limit(g))-log10(lower_limit(g)))
	  exit
	else if (diam_equiv<lower_limit(1)) then
	  frac_finer=0.
	else if (diam_equiv>upper_limit(n_sed_class)) then
	  frac_finer=1.
	endif
  enddo

  DO g=1,n_sed_class
!write(*,'(I3,5F10.6)')g,frac_finer,cumfrsed_in(g),diam_equiv,lower_limit(g),upper_limit(g)
  enddo
!write(*,*)frac_finer,diam_equiv,lower_limit(1)

  if (frac_finer == 1. .or. frac_finer == 0.) then
    trap_eff=1.-frac_finer
	if (frac_finer == 1.) then
      DO g=1,n_sed_class
	    frsediment_out(res_index(upstream),g)=frsediment_in(res_index(upstream),g)
	  enddo
	else
      DO g=1,n_sed_class
	    frsediment_out(res_index(upstream),g)=0.
	  enddo
    endif
  else
    dummy=0.
    DO g=1,n_sed_class
      if (diam_equiv>=lower_limit(g)) then
	    dummy=dummy+1.
	  else
	    exit
	  endif
	enddo
    nbrsteps=int(10.*dummy)
	lower_diam(1)=lower_limit(1)
!	interv=frac_finer/nbrsteps
    k=0
    DO g=1,n_sed_class
	  if (g>1) then
	    if (frac_finer>=cumfrsed_in(g)) &
		  dummy=(cumfrsed_in(g)-cumfrsed_in(g-1))/10.
	    if (frac_finer<cumfrsed_in(g) .and. frac_finer>cumfrsed_in(g-1)) &
		  dummy=(frac_finer-cumfrsed_in(g-1))/10.
	    if (frac_finer<cumfrsed_in(g) .and. frac_finer<cumfrsed_in(g-1)) exit
	  else
		if (frac_finer>=cumfrsed_in(g)) &
		  dummy=cumfrsed_in(g)/10.
	    if (frac_finer<cumfrsed_in(g)) &
		  dummy=frac_finer/10.
	  endif

      do l=1,10
	    k=k+1
	    interv(k)=dummy
	  enddo
	enddo
    DO k=1,nbrsteps
	  if (k>1) frfiner(k)=frfiner(k-1)+interv(k)
	  if (k==1) frfiner(k)=interv(k)
	enddo
    DO k=1,nbrsteps
!	  frfiner(k)=k*interv
      DO g=1,n_sed_class
        if (frfiner(k)<=cumfrsed_in(1)) then
          upper_diam(k)=10.**(log10(lower_limit(g))+((log10(upper_limit(g))-log10(lower_limit(g)))* &
			frfiner(k)/cumfrsed_in(g)))
		  exit
		else if (frfiner(k)<=cumfrsed_in(g)) then
          upper_diam(k)=10.**(log10(lower_limit(g))+((log10(upper_limit(g))-log10(lower_limit(g)))* &
			(frfiner(k)-cumfrsed_in(g-1))/(cumfrsed_in(g)-cumfrsed_in(g-1))))
		  exit
	    endif
	  enddo
!write(*,'(I3,4F10.6)')k,interv(k),frfiner(k),upper_diam(k)
	enddo

    DO k=2,nbrsteps
	  lower_diam(k)=upper_diam(k-1)
	enddo

	frretent=0.
	frrelease=0.
    DO k=1,nbrsteps
	  mean_diam(k)=sqrt(lower_diam(k)*upper_diam(k))
! settling velocity (m/s)
      setveloc(k)=SQRT((13.95*visc/(mean_diam(k)/1000.))**2.+(1.09*  &
        (delta)*9.807*(mean_diam(k)/1000.)))-(13.95*visc/(mean_diam(k)/1000.))
	  par_x(k)=(setveloc(k)/overflow_rate)*interv(k)
	  par_y(k)=(1.-(setveloc(k)/overflow_rate))*interv(k)
!	  par_x(k)=(setveloc(k)/overflow_rate)*interv
!	  par_y(k)=(1-(setveloc(k)/overflow_rate))*interv
      frretent=frretent+par_x(k)
      frrelease=frrelease+par_y(k)
!write(*,'(I3,6F12.8)')k,mean_diam(k),lower_diam(k),upper_diam(k),setveloc(k),par_x(k),par_y(k)
!write(*,'(I3,6F12.8)')k,mean_diam(k),setveloc(k),par_x(k),par_y(k),frretent,frrelease
	enddo

	trap_eff=(1.-frac_finer)+frretent

	dummy=0.
    DO k=1,nbrsteps
	  frfiner_out(k)=par_y(k)/frrelease
	  dummy=dummy+frfiner_out(k)
	  cumfrrelease(k)=dummy
!write(*,'(I3,4F12.8)')k,upper_diam(k),frfiner_out(k),cumfrrelease(k),trap_eff
	enddo

    DO g=1,n_sed_class
	  if (upper_limit(g)>upper_diam(nbrsteps)) then
	    cumfrsed_out(g)=1.
	  else
	    do k=1,nbrsteps
	      if (upper_limit(g)>=lower_diam(k) .and. upper_limit(g)<=upper_diam(k) .and. k>1) then
            cumfrsed_out(g)=cumfrrelease(k-1)+(cumfrrelease(k)-cumfrrelease(k-1))* &
				(log10(upper_limit(g))-log10(lower_diam(k)))/ &
				(log10(upper_diam(k))-log10(lower_diam(k)))
	        exit
		  else if (upper_limit(g)>=lower_diam(1) .and. upper_limit(g)<=upper_diam(1)) then
            cumfrsed_out(g)=cumfrrelease(k)* &
				(log10(upper_limit(g))-log10(lower_diam(k)))/ &
				(log10(upper_diam(k))-log10(lower_diam(k)))
	        exit
		  endif
		enddo
	  endif
	  if (g==1) frsediment_out(res_index(upstream),g)=cumfrsed_out(g)
	  if (g>1) frsediment_out(res_index(upstream),g)=cumfrsed_out(g)-cumfrsed_out(g-1)

!write(*,'(I3,2F12.8)')g,frsediment_out(upstream,g),cumfrsed_out(g)
	enddo
  endif
  sed_outflow(step,res_index(upstream))=sed_inflow(step,res_index(upstream))*(1.-trap_eff)
else
  sed_outflow(step,res_index(upstream))=0.
endif


if (t==tstart .and. step==1) cum_sedimentation(upstream)=0.
IF (t>tstart .and. t== damyear(upstream) .and. step==1) cum_sedimentation(upstream)=0.
sedimentation(step,res_index(upstream))=max(sed_inflow(step,res_index(upstream))-sed_outflow(step,res_index(upstream)),0.)
cum_sedimentation(upstream)=cum_sedimentation(upstream)+sedimentation(step,res_index(upstream))
!sed_susp(step,upstream)=0.


! storage capacity reduction (m3)
decstorcap(step,res_index(upstream))=sedimentation(step,res_index(upstream))/dry_dens(upstream)


!write(*,'(3I4,2F15.1)')t,d,upstream,sedimentation(step,upstream),cum_sedimentation(upstream)

!write(*,'(3I4,4F15.1)')t,d,upstream,sed_inflow(step,upstream),sed_outflow(step,upstream),sedimentation(step,upstream),sed_susp(step,upstream)

!if (d==2)stop

RETURN
END SUBROUTINE sedbal

