SUBROUTINE sedbal_lake(muni,k)
 
! Till: computationally irrelevant: minor changes to improve compiler compatibility
! 2011-04-29

! Code converted using TO_F90 by Alan Miller
! Date: 2005-08-23  Time: 12:57:31
 
use common_h
use time_h
use reservoir_h
use lake_h
use routing_h

IMPLICIT NONE

INTEGER, INTENT(IN OUT)                  :: muni
INTEGER, INTENT(IN OUT)                  :: k
INTEGER :: i,n,g,l,nbrsteps
real:: dummy

real :: wat_dens,sed_dens,delta,grav,visc
real :: overflow_rate,diam_equiv,frac_finer,trap_eff
real :: frcum,frretent,frrelease
real :: cumfrsed_in(n_sed_class),cumfrsed_out(n_sed_class),lower_limit(n_sed_class)
real :: par_a,par_b,par_c
real :: interv(200),setveloc(200),upper_diam(200),lower_diam(200),mean_diam(200),frfiner(200)
real :: frfiner_out(200),cumfrrelease(200)
real :: par_x(200),par_y(200),par_z(200)

!Ge to include the DAILY mean temperature of the reservoir (celsius degree)
!Ge tempres=20 C (temporarily)
REAL :: tempres(step,muni)

! -----------------------------------------------------------------------
if (lakesedin_hrr(step,muni,k)/=0. .and. lakearea(muni,k)/=0.) then
!Ge to include the DAILY mean temperature of the reservoir (celsius degree)
!Ge tempres=20 C (temporarily)
  tempres(step,muni)=20.

! density of water and density of natural sediments (kg/m3)
  wat_dens=1.*1000.
  sed_dens=2.65*1000.

! specific gravity of sediments in water (-)
  delta=(sed_dens-wat_dens)/wat_dens

! gravitacional acceleration (m2/s)
  grav=9.807

! kinematic viscosity (m2/s)
  visc=(1.14-0.031*(tempres(step,muni)-15.)+  &
      0.00068*((tempres(step,muni)-15.)**2.))*1.e-6

! Overflow rate  Vc (m/s)
  overflow_rate=(lakeoutflow_hrr(step,muni,k)/(86400./nt))/(lakearea(muni,k)*1.e6)
!  overflow_rate=0.00032

! calculation of the corresponding equivalent diameter (mm) 
  par_a=1.09*(delta)*grav
  par_b=-(overflow_rate**2)
  par_c=-27.9*overflow_rate*visc

  diam_equiv=1000.*(-par_b+sqrt(((par_b)**2.)-4.*par_a*par_c))/(2.*par_a)

! affluent grain size distribution (-)
  dummy=0.
  DO g=1,n_sed_class
    dummy=dummy+lakefracsedin_hrr(step,muni,k,g)
    cumfrsed_in(g)=dummy
!write(*,'(3I4,3F10.3)')muni,k,g,cumfrsed_in(g)
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
	    lakefracsedout_hrr(step,muni,k,g)=lakefracsedin_hrr(step,muni,k,g)
	  enddo
	else
      DO g=1,n_sed_class
	    lakefracsedout_hrr(step,muni,k,g)=0.
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
    nbrsteps=10.*dummy
	lower_diam(1)=lower_limit(1)
!	interv=frac_finer/nbrsteps
    n=0
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
	    n=n+1
	    interv(n)=dummy
	  enddo
	enddo
    DO n=1,nbrsteps
	  if (n>1) frfiner(n)=frfiner(n-1)+interv(n)
	  if (n==1) frfiner(n)=interv(n)
	enddo
    DO n=1,nbrsteps
!	  frfiner(n)=n*interv
      DO g=1,n_sed_class
        if (frfiner(n)<=cumfrsed_in(1)) then
          upper_diam(n)=10.**(log10(lower_limit(g))+((log10(upper_limit(g))-log10(lower_limit(g)))* &
			frfiner(n)/cumfrsed_in(g)))
		  exit
		else if (frfiner(n)<=cumfrsed_in(g)) then
          upper_diam(n)=10.**(log10(lower_limit(g))+((log10(upper_limit(g))-log10(lower_limit(g)))* &
			(frfiner(n)-cumfrsed_in(g-1))/(cumfrsed_in(g)-cumfrsed_in(g-1))))
		  exit
	    endif
	  enddo
!write(*,'(I3,4F10.6)')n,interv(n),frfiner(n),upper_diam(n)
	enddo

    DO n=2,nbrsteps
	  lower_diam(n)=upper_diam(n-1)
	enddo

	frretent=0.
	frrelease=0.
    DO n=1,nbrsteps
	  mean_diam(n)=sqrt(lower_diam(n)*upper_diam(n))
! settling velocity (m/s)
      setveloc(n)=SQRT((13.95*visc/(mean_diam(n)/1000.))**2.+(1.09*  &
        (delta)*9.807*(mean_diam(n)/1000.)))-(13.95*visc/(mean_diam(n)/1000.))
	  par_x(n)=(setveloc(n)/overflow_rate)*interv(n)
	  par_y(n)=(1.-(setveloc(n)/overflow_rate))*interv(n)
!	  par_x(n)=(setveloc(n)/overflow_rate)*interv
!	  par_y(n)=(1-(setveloc(n)/overflow_rate))*interv
      frretent=frretent+par_x(n)
      frrelease=frrelease+par_y(n)
!write(*,'(I3,6F12.8)')n,mean_diam(n),lower_diam(n),upper_diam(n),setveloc(n),par_x(n),par_y(n)
!write(*,'(I3,6F12.8)')n,mean_diam(n),setveloc(n),par_x(n),par_y(n),frretent,frrelease
	enddo

	trap_eff=(1.-frac_finer)+frretent

	dummy=0.
    DO n=1,nbrsteps
	  frfiner_out(n)=par_y(n)/frrelease
	  dummy=dummy+frfiner_out(n)
	  cumfrrelease(n)=dummy
!write(*,'(I3,4F12.8)')n,upper_diam(n),frfiner_out(n),cumfrrelease(n),trap_eff
	enddo

    DO g=1,n_sed_class
	  if (upper_limit(g)>upper_diam(nbrsteps)) then
	    cumfrsed_out(g)=1.
	  else
	    do n=1,nbrsteps
	      if (upper_limit(g)>=lower_diam(n) .and. upper_limit(g)<=upper_diam(n) .and. n>1) then
            cumfrsed_out(g)=cumfrrelease(n-1)+(cumfrrelease(n)-cumfrrelease(n-1))* &
				(log10(upper_limit(g))-log10(lower_diam(n)))/ &
				(log10(upper_diam(n))-log10(lower_diam(n)))
	        exit
		  else if (upper_limit(g)>=lower_diam(1) .and. upper_limit(g)<=upper_diam(1)) then
            cumfrsed_out(g)=cumfrrelease(n)* &
				(log10(upper_limit(g))-log10(lower_diam(n)))/ &
				(log10(upper_diam(n))-log10(lower_diam(n)))
	        exit
		  endif
		enddo
	  endif
	  if (g==1) lakefracsedout_hrr(step,muni,k,g)=cumfrsed_out(g)     
	  if (g>1) lakefracsedout_hrr(step,muni,k,g)=cumfrsed_out(g)-cumfrsed_out(g-1)    

!write(*,'(3I4,3F10.3)')muni,k,g,lakefracsedout_hrr(step,muni,k,g),cumfrsed_out(g)
!write(*,'(3I4,3F10.3)')muni,k,g,cumfrsed_in(g),cumfrsed_out(g)
!write(*,'(3I4,3F10.3)')muni,k,g,lakefracsedin_hrr(step,muni,k,g),lakefracsedout_hrr(step,muni,k,g)
	enddo		    
  endif
  lakesedout_hrr(step,muni,k)=lakesedin_hrr(step,muni,k)*(1.-trap_eff)
else
  lakesedout_hrr(step,muni,k)=0.
endif

!write(*,'(3I4,5F12.2)')step,muni,k,acud(muni,k),lakeoutflow_hrr(step,muni,k),lakesedin_hrr(step,muni,k),lakesedout_hrr(step,muni,k),trap_eff

if (t==tstart .and. step==1) lakecumdep_hrr(muni,k)=0.
lakedep_hrr(step,muni,k)=max(lakesedin_hrr(step,muni,k)-lakesedout_hrr(step,muni,k),0.)
lakecumdep_hrr(muni,k)=lakecumdep_hrr(muni,k)+lakedep_hrr(step,muni,k)


! storage capacity reduction (m3)
lake_vollost(step,muni,k)=lakedep_hrr(step,muni,k)/1.5
!lake_vollost(step,muni,k)=lake_vollost(step,muni,k)*acud(muni,k)


!write(*,'(3I4,2F15.1)')t,d,upstream,sedimentation(step,upstream),cum_sedimentation(upstream)

!write(*,'(3I4,4F15.1)')t,d,upstream,sed_inflow(step,upstream),sed_outflow(step,upstream),sedimentation(step,upstream),sed_susp(step,upstream)

!if (d==2)stop

RETURN
END SUBROUTINE sedbal_lake

