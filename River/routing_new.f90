SUBROUTINE routing_new(STATUS)
!Till: computationally irrelevant: quick-fix of severe performance loss when using output of river_flow.dat in hourly version
!2009-04-02

!Till: added output for River_Sediment_Storage.out
!2008-11-13

!(George:) various changes concerning routing of reservoir fluxes
!2008-10-07 

!Till: file and header of river_sediment_total was not created when river_transport=3
!2008-07-11

!! Routing of water and sediment fluxes through the river system
!!
!! Status 0: reads input files and initialises run, calls subroutine reservoir (status 0)
!! Status 1: initialises year, calls subroutine reservoir (status 1) 
!! Status 2: main calculations, calls subroutine reservoir (status 2)
!! Status 3: finalises year, calls subroutine reservoir (status 3)
!! Parameter file: routing_h.f90

!!	Q_spring	(m3/s)		spring inflow at the start of a river
 
use lake_h
use climo_h
use common_h
use hymo_h
use params_h
use routing_h
use time_h
use reservoir_h
use erosion_h

IMPLICIT NONE

INTEGER, INTENT(IN OUT) :: STATUS
!INTEGER :: bat
INTEGER :: idummy !,imun,imunx,irout,irout2,irout_d,id,imeso,istate ! id: additional loop variable of days (total 7 days)
INTEGER :: upstream, downstream
INTEGER :: i, j, h,k !itl, itr, ih, mm, imunout, iout,  make
!REAL :: xdum(48),check,temp2,qtemp, storcapact, con_sed
REAL :: flow, r_area !, sediment_temp(24), temp_rain(366), dummy
Real :: temp_water(17), temp_sediment(17)
character(len=1000) :: fmtstr	!string for formatting file output
Real :: r_sediment_storage(subasin)		!sediment storage in reach [t]

! -----------------------------------------------------------------------
IF (STATUS == 0) THEN

! READ INPUT FILES
! Read hydrological response and river paramters
 OPEN(11,FILE=pfadp(1:pfadj)// 'River/river.dat'  &
    ,STATUS='old')
  READ (11,*); READ(11,*)
  DO i=1,subasin
   READ (11,*) idummy, r_depth(i),r_width(i), r_sideratio(i),r_width_fp(i),r_sideratio_fp(i), &
   r_slope(i), r_length(i),manning(i), manning_fp(i),r_ksat(i),r_efactor(i),r_cover(i),r_rock(i),r_alpha(i), &
   msk_x(i), msk_k(i),Q_spring(i)
     IF (idummy /= id_subbas_extern(i)) THEN
     WRITE(*,*) 'Sub-basin-IDs in file river.dat must have the same ordering scheme as in hymo.dat'
     STOP 
   END IF
  END DO
  CLOSE (11)

! READ INPUT FILES
! Read configuration of river system (in routing order)
OPEN(11,FILE=pfadp(1:pfadj)// 'River/routing.dat',STATUS='old')! upbasin: MAP ID of upstream sub-basin (MAP IDs);! downbasin: MAP ID of downstream sub-basin (MAP IDs)
READ (11,*); READ(11,*)
DO i=1,subasin
  READ (11,*)  idummy, upbasin(i),downbasin(i)
END DO
CLOSE (11)


! if bedload modelleing is switched on
if(river_transport.eq.3) then
  OPEN(11,FILE=pfadp(1:pfadj)// 'River/bedload.dat',STATUS='old')!
  read(11,*); read(11,*)
  DO i=1,subasin
	READ (11,*) idummy, D50(i)
	IF (idummy /= id_subbas_extern(i)) THEN
		WRITE(*,*) 'Sub-basin-IDs in file river.dat must have the same ordering scheme as in hymo.dat'
		STOP 
	END IF
  END DO
endif

!this relates the MAP IDs (id_subbas_extern(subasin)) to the sorted CODE IDs (id_subbas_intern(subasin)), i.e. upbasin and downbasin are now numbered according to internal ids
!so that in routing.dat only the MAP IDs have to be read in, first for the ID of upstream subasin
DO i=1,subasin
  j=1
  DO WHILE (id_subbas_extern(j) /= upbasin(i))
    j=j+1
    IF (j > 1000) THEN
      WRITE (*,*) 'upbasin(i) loop in routing_new.f'
      STOP
    END IF
  END DO
  upbasin(i)=j
END DO
! second for the ID of downstream subasin
DO i=1,subasin
  IF (downbasin(i) /= 999.AND.downbasin(i) /= 9999) THEN
    j=1
    DO WHILE (id_subbas_extern(j) /= downbasin(i))
      j=j+1
      IF (j > 1000) THEN
        WRITE (*,*) 'downbasin(i) loop in routing_new.f'
        STOP
      END IF
    END DO
    downbasin(i)=j
  END IF
END DO


! INITIALISATION OF RESERVOIR MODULE
!George (status,upstream,h) instead (status,upstream)
  IF (doreservoir) CALL reservoir (status,upstream,h)
 

! INITIALISATION OF OUTPUT FILES
  OPEN(111,FILE=pfadn(1:pfadi)//'River_Flow.out',STATUS='replace')
  if (f_river_flow) then
    WRITE (111,*) 'Output files for river discharge q_out (m3/s) (with MAP IDs as in hymo.dat)'
	    
	if (dohour) then
		write(fmtstr,'(a,i0,a)')'(3a6,',subasin,'i14)'		!generate format string
		WRITE (111,fmtstr)' Year ', ' Day  ','  dt  ', (id_subbas_extern(i), i=1,subasin)
	else
   		write(fmtstr,'(a,i0,a)')'(2a6,',subasin,'i14)'		!generate format string
		WRITE (111,fmtstr)' Year ', ' Day  ', (id_subbas_extern(i), i=1,subasin)
		Close (111)		!Till: hourly version leaves the file open to save time
	endif
	
  else
    close(111, status='delete') !delete any existing file, if no output is desired
  endif


if ((river_transport.eq.2)) then
  OPEN(11,FILE=pfadn(1:pfadi)//'River_Sediment_total.out',STATUS='replace')
  if(f_river_sediment_total) then
	WRITE (11,*) 'Output file for sediment mass in ton/timestep (with MAP IDs as in hymo.dat)'
	if (dohour) then
		write(fmtstr,'(a,i0,a)')'(3a6,',subasin,'i14)'		!generate format string
		WRITE (11,fmtstr)' Year ', ' Day  ','  dt  ', (id_subbas_extern(i), i=1,subasin)
		!WRITE (11,'(3a6,<subasin>i14)')' Year ', ' Day  ','  dt  ', (id_subbas_extern(i), i=1,subasin)
	else
		write(fmtstr,'(a,i0,a)')'(2a6,',subasin,'i14)'		!generate format string
		WRITE (11,fmtstr)' Year ', ' Day  ', (id_subbas_extern(i), i=1,subasin)
		!WRITE (11,'(2a6,<subasin>i14)')' Year ', ' Day  ', (id_subbas_extern(i), i=1,subasin)
	endif
	Close (11)
  else
    close(11, status='delete') !delete any existing file, if no output is desired
  endif

  OPEN(11,FILE=pfadn(1:pfadi)//'River_Sediment_Concentration.out',STATUS='replace')
  if (f_river_sediment_concentration) then
	WRITE (11,*) 'Output file for sediment concentration (g/l) (with MAP IDs as in hymo.dat)'
	write(fmtstr,'(a,i0,a)')'(3a6,',subasin,'i14)'		!generate format string
	WRITE (11,fmtstr)' Year ', ' Day  ','  dt  ', (id_subbas_extern(i), i=1,subasin)
	!WRITE (11,'(3a6,<subasin>i14)')' Year ', ' Day  ','  dt  ', (id_subbas_extern(i), i=1,subasin)
	Close (11)
  else
    close(11, status='delete') !delete any existing file, if no output is desired
  endif

elseif (river_transport.eq.3) then
  if (f_river_bedload) then
	OPEN(11,FILE=pfadn(1:pfadi)//'River_bedload.out',STATUS='replace')
	WRITE (11,*) 'Output file for river bedload rate (kg/s) as submerged weight(with MAP IDs as in hymo.dat)'
	WRITE (11,*) ' Year ', ' Day  ',' dt   ',&
		'Meyer_Peter, Schoklitsch, Smart&Jaeggi, Bagnold, Rickenmann for each subasin in successive columns'
	Close (11)
  else
    close(11, status='delete') !delete any existing file, if no output is desired
  endif
endif

  OPEN(11,FILE=pfadn(1:pfadi)//'River_Velocity.out',STATUS='replace')
  if (f_river_velocity) then
	WRITE (11,*) 'Output file for flow velocity in m/s (with MAP IDs as in hymo.dat)'
	write(fmtstr,'(a,i0,a)')'(3a6,',subasin,'i14)'		!generate format string
	WRITE (11,fmtstr)' Year ', ' Day  ',' dt   ',(id_subbas_extern(i), i=1,subasin)
	Close (11)
  else
    close(11, status='delete') !delete any existing file, if no output is desired
  endif

  OPEN(11,FILE=pfadn(1:pfadi)//'River_Flowdepth.out',STATUS='replace')
  if (f_river_flowdepth) then
	WRITE (11,*) 'Output file for flow depth in m (with MAP IDs as in hymo.dat)'
	write(fmtstr,'(a,i0,a)')'(3a6,',subasin,'i14)'		!generate format string
	WRITE (11,fmtstr)' Year ', ' Day  ','  dt  ', (id_subbas_extern(i), i=1,subasin)
	Close (11)
  else
    close(11, status='delete') !delete any existing file, if no output is desired
  endif

  OPEN(11,FILE=pfadn(1:pfadi)//'River_Storage.out',STATUS='replace')
  if (f_river_storage) then
	WRITE (11,*) 'Output file for river water storage in m3 (with MAP IDs as in hymo.dat)'
	write(fmtstr,'(a,i0,a)')'(3a6,',subasin,'i14)'		!generate format string
	WRITE (11,fmtstr)' Year ', ' Day  ','  dt  ', (id_subbas_extern(i), i=1,subasin)
	Close (11)
   else
    close(11, status='delete') !delete any existing file, if no output is desired
  endif

  OPEN(11,FILE=pfadn(1:pfadi)//'River_Sediment_Storage.out',STATUS='replace')
  if (f_river_sediment_storage) then
	WRITE (11,*) 'Output file for river sediment storage in t (with MAP IDs as in hymo.dat)'
	write(fmtstr,'(a,i0,a)')'(5a,',subasin,'(a,i14))'		!generate format string
	WRITE (11,fmtstr)' Year ',char(9), ' Day  ',char(9),'  dt  ', (char(9),id_subbas_extern(i), i=1,subasin)

	Close (11)
   else
    close(11, status='delete') !delete any existing file, if no output is desired
  endif

  OPEN(11,FILE=pfadn(1:pfadi)//'River_Deposition.out',STATUS='replace')
  if (f_river_deposition) then
	WRITE (11,*) 'Output file for deposition of sediments in the riverbed in tons/river stretch (with MAP IDs as in hymo.dat)'
	write(fmtstr,'(a,i0,a)')'(3a6,',subasin,'i14)'		!generate format string
	WRITE (11,fmtstr)' Year ', ' Day  ','  dt  ', (id_subbas_extern(i), i=1,subasin)
	Close (11)
  else
    close(11, status='delete') !delete any existing file, if no output is desired
  endif

  OPEN(11,FILE=pfadn(1:pfadi)//'River_Degradation.out',STATUS='replace')
  if (f_river_degradation) then
	WRITE (11,*) 'Output file for erosion of sediments in the riverbed in tons/river stretch (with MAP IDs as in hymo.dat)'
	write(fmtstr,'(a,i0,a)')'(3a6,',subasin,'i14)'		!generate format string
	WRITE (11,fmtstr)' Year ', ' Day  ','  dt  ', (id_subbas_extern(i), i=1,subasin)
	Close (11)
  else
    close(11, status='delete') !delete any existing file, if no output is desired
  endif


  OPEN(11,FILE=pfadn(1:pfadi)//'River_Flow_dailyaverage.out',STATUS='replace')
  if (f_river_flow_dailyaverage) then
	WRITE (11,*) 'Output files for river discharge q_out (m3/s) (with MAP IDs as in hymo.dat)'
	write(fmtstr,'(a,i0,a)')'(3a6,',subasin,'i14)'		!generate format string
	WRITE (11,fmtstr)' Year ', ' Day  ','  dt  ', (id_subbas_extern(i), i=1,subasin)
	Close (11)
  else
    close(11, status='delete') !delete any existing file, if no output is desired
  endif

  OPEN(11,FILE=pfadn(1:pfadi)//'River_Sediment_total_dailyaverage.out',STATUS='replace')
  if (f_river_sediment_total_dailyaverage) then
	WRITE (11,*) 'Output file for sediment mass in ton/timestep (with MAP IDs as in hymo.dat)'
	write(fmtstr,'(a,i0,a)')'(3a6,',subasin,'i14)'		!generate format string
	WRITE (11,fmtstr)' Year ', ' Day  ','  dt  ', (id_subbas_extern(i), i=1,subasin)
	Close (11)
  else
    close(11, status='delete') !delete any existing file, if no output is desired
  endif

END IF
! ------------------------------------------------------------------------
IF (STATUS == 1) THEN
  
! Initialisation for first year and day
 IF (t == tstart) THEN
  DO i=1,subasin
   do j=1,2
	  r_qout(j,i)=0.
	  r_qin(j,i)=0.
   enddo
   r_depth_cur(i)=0.
   if(river_transport.eq.2) then
     sediment_total(i) = 0.
     do k=1,n_sed_class
      sediment_in(i,k)=0.
	  sediment_out(i,k)=0.
	  sed_storage(i,k)=0.
	  river_deposition(i,k)=0.
	  river_degradation(i,k)=0.
	  riverbed_storage(i,k)=0.
     enddo
	elseif(river_transport.eq.3) then
	  bedload(i,:) = 0.
	endif
  enddo

  do i=1,subasin 
   r_qin(1,i) = Q_spring(i)	!spring water, start of a river
   r_qin(2,i) = Q_spring(i) !spring water, start of a river
   r_qout(1,i)= Q_spring(i) !guestimation of the outflow discharge
  END DO

 ENDIF

  
! CALL Reservoir Sedimentation and Management Modules
!George (status,upstream,h) instead (status,upstream)
 IF (doreservoir) CALL reservoir (status,upstream,h)

END IF

! ------------------------------------------------------------------------
IF (STATUS == 2) THEN

! Calculation of inflow and outflow in stream flow order
DO h=1,nt					
  DO i=1,subasin

  !Water (m3/s) and sediment (ton/dt) fluxes from the hillslope module are added to river stretch
	r_qin(2,i) = water_subbasin_t(d,h,i) + Q_Spring(i)
	if(dosediment.and.river_transport.eq.2) then
		do k= 1,n_sed_class
				sediment_in(i,k) = sediment_subbasin_t(d,h,i,k)
		enddo
    endif !(dosediment)
  enddo	!(i=1, subasin)

!George (initialization of the variable that calculates lateral water inflow into the reservoirs)
!*************************************************************************
  IF (doreservoir) THEN
   if (dosediment) then
    DO i=1,subasin
	  do k= 1,n_sed_class
	    sed_qlateral(i,k)=0.
	  enddo
    ENDDO
   endif
  ENDIF
!*************************************************************************

  DO i=1,subasin
   upstream=upbasin(i)  !internal code-ID for most upstream sub-basin
   downstream=downbasin(i) !internal code-ID for receiving sub-basin

   if (do_pre_outflow .AND. (corr_column_pre_subbas_outflow(i)>0)) then		!Till: if water outflow from upstream subbasins is given
     r_qout(2,upstream)=water_subbasin_t(d,h,upstream)
   else
     call muskingum (upstream, flow, r_area,h)								!normal routing
   end if 

   if (do_pre_outsed .AND. (corr_column_pre_subbas_outsed(i)>0)) then		!Till: if water outsed from upstream subbasins is given and if outsed of subbasin is prespecified
	do k=1,n_sed_class
		sediment_out(upstream,k)=sediment_subbasin_t(d,h,upstream,k)
	end do
   else																		!normal sediment routing
	if (dosediment) then
     if (river_transport.eq.2) call route_sediments(upstream, flow, r_area)
     if(river_transport.eq.3) call bedload_formulae(upstream)
	endif
   end if !do outsed

   
!George Loop was changed in April 2007 after including lateral inflow into the dowstream reservoirs
!*************************************************************************************
!   IF (downstream /= 9999 .AND. downstream /= 999) THEN !George
!     r_qin(2,downstream)=r_qin(2,downstream) + r_qout(2,upstream) !George
!	 if(dosediment.and.river_transport.eq.2) then !George
!	   do k=1,n_sed_class !George
!		 sediment_in(downstream, k)=sediment_in(downstream,k) + sediment_out(upstream,k) !George
!	   enddo !George
!	 endif !George
!   END IF !George

   IF (doreservoir) THEN
!George calculation of the simulation timestep "step"
     step=(d-1)*nt+h
     IF (storcap(upstream) > 0. .and. t >= damyear(upstream)) THEN
!George (status,upstream,h) instead (status,upstream)
       CALL reservoir (status,upstream,h)
	 ENDIF
   ENDIF

   IF (doreservoir) THEN
     IF (latflow_res(upstream) == 0) THEN
       IF (storcap(upstream) > 0. .and. t >= damyear(upstream)) THEN
         IF (downstream /= 9999 .AND. downstream /= 999) THEN
!George res_qout(step,upstream) instead qout(d,upstream)
!write(*,*)step,id_subbas_extern(upstream),"case 1a"
           r_qin(2,downstream)=r_qin(2,downstream)+res_qout(step,upstream)
	       if(dosediment) then
	         do k=1,n_sed_class
	           sediment_in(downstream,k)=sediment_in(downstream,k)+res_sediment_out(upstream,k)
	         enddo
	       endif
         END IF
	   ELSE
         IF (downstream /= 9999 .AND. downstream /= 999) THEN !George
!write(*,*)step,id_subbas_extern(upstream),"case 1b"
           r_qin(2,downstream)=r_qin(2,downstream) + r_qout(2,upstream) !George
	       if(dosediment.and.river_transport.eq.2) then !George
	         do k=1,n_sed_class !George
		       sediment_in(downstream, k)=sediment_in(downstream,k) + sediment_out(upstream,k) !George
	         enddo !George
	       endif !George
         END IF !George
	   ENDIF
	 ELSE IF (latflow_res(upstream) == 1) THEN
       IF (storcap(upstream) > 0. .and. t >= damyear(upstream)) THEN
         IF (storcap(reservoir_down(upstream)) > 0. .and. t >= damyear(reservoir_down(upstream))) THEN
!write(*,*)step,id_subbas_extern(upstream),"case 2a"
           qlateral(step,reservoir_down(upstream))=qlateral(step,reservoir_down(upstream))+res_qout(step,upstream)
	       if(dosediment) then
	         do k=1,n_sed_class
	           sed_qlateral(reservoir_down(upstream),k)=sed_qlateral(reservoir_down(upstream),k)+res_sediment_out(upstream,k)
	         enddo
	       endif
		 ELSE
           IF (downstream /= 9999 .AND. downstream /= 999) THEN !George
!write(*,*)step,id_subbas_extern(upstream),"case 2b"
             r_qin(2,downstream)=r_qin(2,downstream) + res_qout(step,upstream) !George
	         if(dosediment.and.river_transport.eq.2) then !George
	           do k=1,n_sed_class !George
		         sediment_in(downstream, k)=sediment_in(downstream,k) + res_sediment_out(upstream,k) !George
	           enddo !George
	         endif !George
           END IF !George
		 ENDIF
	   ELSE
         IF (storcap(reservoir_down(upstream)) > 0. .and. t >= damyear(reservoir_down(upstream))) THEN
!write(*,*)step,id_subbas_extern(upstream),"case 3a"
           qlateral(step,reservoir_down(upstream))=qlateral(step,reservoir_down(upstream))+r_qout(2,upstream)
	       if(dosediment) then
	         do k=1,n_sed_class
	           sed_qlateral(reservoir_down(upstream),k)=sed_qlateral(reservoir_down(upstream),k)+sediment_out(upstream,k)
	         enddo
	       endif
		 ELSE
           IF (downstream /= 9999 .AND. downstream /= 999) THEN !George
!write(*,*)step,id_subbas_extern(upstream),"case 3b"
             r_qin(2,downstream)=r_qin(2,downstream) + r_qout(2,upstream) !George
	         if(dosediment.and.river_transport.eq.2) then !George
	           do k=1,n_sed_class !George
		         sediment_in(downstream, k)=sediment_in(downstream,k) + sediment_out(upstream,k) !George
	           enddo !George
	         endif !George
           END IF !George
		 ENDIF
	   ENDIF
	 ENDIF
   ELSE
     IF (downstream /= 9999 .AND. downstream /= 999) THEN !George
!write(*,*)step,id_subbas_extern(upstream),"case 4"
       r_qin(2,downstream)=r_qin(2,downstream) + r_qout(2,upstream) !George
	   if(dosediment.and.river_transport.eq.2) then !George
	     do k=1,n_sed_class !George
		   sediment_in(downstream, k)=sediment_in(downstream,k) + sediment_out(upstream,k) !George
	     enddo !George
	   endif !George
     END IF !George
   ENDIF ! doreservoir
!do k=1,n_sed_class
!dummy=dummy+sediment_subbasin_t(d,h,i,k)
!enddo
!write(*,*)d,id_subbas_extern(upstream),dummy,(sediment_out(upstream,k),k=1,n_sed_class)
!if (d==2)stop

   if (dosediment) then
	 do k=1,n_sed_class !George
	  IF (doreservoir) THEN
       sedinflow_g(step,upstream,k)=sediment_in(upstream,k)+sed_qlateral(upstream,k)
       IF (storcap(upstream) > 0. .and. t >= damyear(upstream)) THEN
         sedoutflow_g(step,upstream,k)=res_sediment_out(upstream,k)
	   ELSE
         sedoutflow_g(step,upstream,k)=sedinflow_g(step,upstream,k)
	   ENDIF
	  ENDIF
	 enddo
   endif
!*************************************************************************

!for daily averages
if (dohour) then
  temp_water(i) = temp_water(i) + r_qout(2,i)
  temp_sediment(i) = temp_sediment(i) + sediment_out(i,1)
endif


	
  END DO ! i=1,subasin

	r_sediment_storage=sum(sed_storage,dim=2)  !Till: sum up sediment storage over all particle classes for all subbasins

! add up all sediment size classes to obtain total sediment mass
  do i = 1, subasin
    do k = 1, n_sed_class
     sediment_total(i) = sediment_total(i) + sediment_out(i,k)
    enddo
  enddo

! calculate total sediment concentration
   do i = 1, subasin
     if (r_qout(2,i).eq. 0.) then
        r_sediment_concentration(i) = 0.
     else
        r_sediment_concentration(i) = sediment_total(i)/(r_qout(2,i)*3.6*dt)
     endif
  enddo
 
 if (dohour) then
	if (f_river_flow) then	!daily output in routing_new(3) to decrease simulation time
!		OPEN(11,FILE=pfadn(1:pfadi)//'River_Flow.out',STATUS='old' ,POSITION='append'  )
		write(fmtstr,'(a,i0,a)')'(3i6,',subasin,'f14.3)'		!generate format string
		WRITE (111,fmtstr)t,d,h, (r_qout(2,i),i=1,subasin)
!		CLOSE (11)
	endif
 endif

if (river_transport.eq.2) then !daily output in routing_new(3) to decrease simulation time
	if (dohour) then
	  if (f_river_sediment_total) then
		OPEN(11,FILE=pfadn(1:pfadi)//'River_Sediment_total.out',STATUS='old' ,POSITION='append'  )
		write(fmtstr,'(a,i0,a)')'(3i6,',subasin,'f14.3)'		!generate format string
		WRITE (11,fmtstr)t,d,h, (sediment_total(i),i=1,subasin)
		CLOSE (11) 
	  endif
	endif

  if (f_river_sediment_concentration) then
	OPEN(11,FILE=pfadn(1:pfadi)//'River_Sediment_Concentration.out',STATUS='old' ,POSITION='append'  )
	write(fmtstr,'(a,i0,a)')'(3i6,',subasin,'f14.3)'		!generate format string
	WRITE (11,fmtstr)t,d,h, (r_sediment_concentration(i),i=1,subasin)
	CLOSE (11) 
  endif
elseif(river_transport.eq.3) then
   if(f_river_bedload) then	
	OPEN(11,FILE=pfadn(1:pfadi)//'River_Bedload.out',STATUS='old' ,POSITION='append'  )
	write(fmtstr,'(a,i0,a)')'(3i6,',5*subasin,'f14.3)'		!generate format string
	WRITE (11,fmtstr)t,d,h, (bedload(i,1), bedload(i,2),bedload(i,3), bedload(i,4), bedload(i,5), i=1,subasin)
	CLOSE (11) 
   endif
endif

  if(f_river_velocity) then
	OPEN(11,FILE=pfadn(1:pfadi)//'River_Velocity.out',STATUS='old' ,POSITION='append'  )
	write(fmtstr,'(a,i0,a)')'(3i6,',subasin,'f14.3)'		!generate format string
	WRITE (11,fmtstr)t, d,h, (velocity(i),i=1,subasin)
	CLOSE (11) 
  endif

  if(f_river_storage) then
	OPEN(11,FILE=pfadn(1:pfadi)//'River_Storage.out',STATUS='old' ,POSITION='append'  )
	write(fmtstr,'(a,i0,a)')'(3i6,',subasin,'f14.3)'		!generate format string
	WRITE (11,fmtstr)t, d,h, (r_storage(i),i=1,subasin)
	CLOSE (11) 
  endif

  if(f_river_sediment_storage) then
	OPEN(11,FILE=pfadn(1:pfadi)//'River_Sediment_Storage.out',STATUS='old' ,POSITION='append'  )
	write(fmtstr,'(a,i0,a)')'(i0,a,i0,a,i0,',subasin,'(a,f0.1))'		!generate format string
	WRITE (11,fmtstr)t,char(9),d,char(9),h,(char(9),r_sediment_storage(i),i=1,subasin)
	CLOSE (11) 
  endif

  if(f_river_flowdepth) then
	OPEN(11,FILE=pfadn(1:pfadi)//'River_Flowdepth.out',STATUS='old' ,POSITION='append'  )
	write(fmtstr,'(a,i0,a)')'(3i6,',subasin,'f14.3)'		!generate format string
	WRITE (11,fmtstr)t, d,h, (r_depth_cur(i),i=1,subasin)
	CLOSE (11) 
  endif

  if(f_river_deposition) then
	OPEN(11,FILE=pfadn(1:pfadi)//'River_Deposition.out',STATUS='old' ,POSITION='append'  )
	write(fmtstr,'(a,i0,a)')'(3i6,',subasin,'f14.3)'		!generate format string
	WRITE (11,fmtstr)t, d,h, (river_deposition(i,1),i=1,subasin)
	CLOSE (11)
  endif

  if(f_river_degradation) then
	OPEN(11,FILE=pfadn(1:pfadi)//'River_Degradation.out',STATUS='old' ,POSITION='append'  )
	write(fmtstr,'(a,i0,a)')'(3i6,',subasin,'f14.3)'		!generate format string
	WRITE (11,fmtstr)t, d,h, (river_degradation(i,1),i=1,subasin)
	CLOSE (11)
  endif

 

 ! for daily averages
  if (dohour.and.h.eq.24) then
    if (f_river_flow_dailyaverage) then
		OPEN(11,FILE=pfadn(1:pfadi)//'River_Flow_dailyaverage.out',STATUS='old' ,POSITION='append')
		write(fmtstr,'(a,i0,a)')'(3i6,',subasin,'f14.3)'		!generate format string
		WRITE (11,fmtstr)t,d,h, ((temp_water(i)/24.),i=1,subasin)
		Close (11)
	endif

	if (f_river_sediment_total_dailyaverage) then
		OPEN(11,FILE=pfadn(1:pfadi)//'River_Sediment_total_dailyaverage.out',STATUS='old' ,POSITION='append')
		write(fmtstr,'(a,i0,a)')'(3i6,',2*subasin,'f14.3)'		!generate format string
		WRITE (11,fmtstr)t,d,h, (temp_sediment(i)/24.,i=1,subasin)
		Close (11)
	endif

    temp_water(:) = 0.
	temp_sediment(:) = 0.
  endif

! Update flows for next timestep 
  DO i=1,subasin
    qout(d,i)= r_qout(2,i)
	qsediment(d,i)=sediment_total(i)
	r_qin(1,i) = r_qin(2,i)			
	r_qout(1,i)= r_qout(2,i)	
	r_qin(2,i)=0.
	r_qout(2,i)=0.
	sediment_total(i) = 0.
	do k=1,n_sed_class
	  sediment_in(i,k)=0.
	enddo
  ENDDO
 
  
ENDDO	! end of day or hourly loop
  
 




END IF


! -----------------------------------------------------------------------
IF (STATUS == 3) THEN
  IF (doreservoir) THEN
    CALL  reservoir (status,upstream,h)
  END IF


! daily output of water and discharge discharge in the river for entire year
	if (.not.dohour) then
		if (f_river_flow) then
			OPEN(11,FILE=pfadn(1:pfadi)//'River_Flow.out',STATUS='old' ,POSITION='append'  )
			write(fmtstr,'(a,i0,a)')'(2i6,',subasin,'f14.3)'		!generate format string
			do j=1, dayyear
				WRITE (11,fmtstr)t, j, (qout(j,i),i=1,subasin)
			enddo
			CLOSE (11)
		endif

!		if (river_transport.eq.2) then 
!			if(f_river_sediment_total) then
!				OPEN(11,FILE=pfadn(1:pfadi)//'River_Sediment_total.out',STATUS='old' ,POSITION='append'  )
!				write(fmtstr,'(a,i0,a,i0,a)')'(',dayyear,'(2(i0,a),',subasin,'(f14.3,a)))'		!generate format string
!				WRITE (11,fmtstr)(t,char(9),j,char(9),(qsediment(j,i),char(9),i=1,subasin),j=1, dayyear)
!				CLOSE (11) 
!			endif
!		endif
	endif

END IF




RETURN
END SUBROUTINE routing_new

