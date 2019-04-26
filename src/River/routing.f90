SUBROUTINE routing(STATUS)

use lake_h
use climo_h
use common_h
use hymo_h
use params_h
use routing_h
use time_h
use reservoir_h
use utils_h
use erosion_h

IMPLICIT NONE

INTEGER, INTENT(IN)                  :: STATUS


! status of call (0=initialization run, 1=initialization year,
!                 2=calculation day,    3=finalization year)

!INTEGER :: bat

INTEGER :: irout,idummy,id !,imun,imunx,irout2,irout_d,imeso,istate
INTEGER :: upstream, downstream
INTEGER :: itl, ih, i, j, istate, h !, mm, imunout, iout, make
REAL :: temp2, temp3, temp4, qtemp, qtemp_sed(n_sed_class) !, x, b, y, hi  !,xdum(48),storcapact
character(len=1000) :: fmtstr	!string for formatting file output


! -----------------------------------------------------------------------
IF (STATUS == 0) THEN

!**  Read hydrological response and reservoir paramter
OPEN(11,FILE=pfadp(1:pfadj)// 'River/response.dat', STATUS='old', IOSTAT=istate)
IF (istate/=0) THEN
		write(*,*)'Error (response.dat): File not found'
		stop
 END IF

prout=0.
READ (11,*, IOSTAT=istate); READ(11,*, IOSTAT=istate)
h=2
i=0 !count treated subbasins

DO WHILE (i<subasin)
  h=h+1
  READ (11,*, IOSTAT=istate)  idummy,temp2,  temp3
  IF (istate/=0) THEN
		write(*,'(a, i0)')'Error (response.dat): Format error or unexpected end in line ',h
		stop
  END IF

  j=which1(idummy == id_subbas_extern) !relate to external IDs from routing.dat

  if (j==0) then
	write(*,'(a, i0, a)')'Warning (response.dat): Unknown subbasin ',idummy,', skipped.'
  else
	prout(j,1)=temp2
	prout(j,2)=temp3
    if (temp2 + temp3 > 200) then
        write(*,'(a, i0, a)')'Warning (response.dat): Subbasin ',idummy,': lag+retention must be < 200, rescaled.'
        prout(j,1)=prout(j,1) * 200/(temp2+temp3)
	    prout(j,2)=prout(j,2) * 200/(temp2+temp3)
    end if     
	i=i+1
  end if

END DO
CLOSE (11)

!this relates the external IDs (id_subbas_extern(subasin)) to the sorted CODE IDs (id_subbas_intern(subasin)),
! i.e. upbasin and downbasin are now numbered according to internal ids
!     so that in routing.dat only the MAP IDs have to be read in,

!replace external with internal IDs
DO i=1,subasin
  !upstream basin referencing
  ih=which1(id_subbas_extern == upbasin(i))

  IF (ih==0) THEN
      WRITE (*,'(A,I0,A)') 'ERROR: unknown upstream subbasin ID ', upbasin(i),' in routing.dat'
      STOP
  else
	upbasin(i)=ih
  END IF

  !downstream basin referencing
  IF (downbasin(i) == 999 .OR. downbasin(i) == 9999) cycle 	!999 and 9999 mark outlet
  ih=which1(id_subbas_extern == downbasin(i))

  IF (ih==0) THEN
      WRITE (*,'(A,I0,A)') 'ERROR: unknown downstream subbasin ID ', downbasin(i),' in routing.dat'
      STOP
  else
	downbasin(i)=ih
  END IF
END DO



! INITIALISATION OF OUTPUT FILES

OPEN(11,FILE=pfadn(1:pfadi)//'River_Flow.out',STATUS='replace')
  if (f_river_flow) then
    WRITE (11,*) 'Output files for river discharge q_out (m3/s) (with MAP IDs as in hymo.dat)'
	write(fmtstr,'(a,i0,a)')'(2a6,',subasin,'i14)'		!generate format string
	WRITE (11,fmtstr)' year ', ' day  ', (id_subbas_extern(i), i=1,subasin)
	!WRITE (11,'(2a6,<subasin>i14)')' year ', ' day  ', (id_subbas_extern(i), i=1,subasin)
    Close (11)
  else
    close(11, status='delete') !delete any existing file, if no output is desired
  endif

! calculate external routing response function for each sub-basin (triangular like this: _|\_ ) - used for riverflow ENTERING the subbasin from upstream
! (given parameter lag-time tL and retention-time tR)
  !ii: the dimensions must be computed differently when delta_t is not 24 hours!
  allocate( hrout       (maxval(ceiling(1+sum(prout, dim=2))) ,subasin)) !allocate memory for triangular unit hydrograph (allochtonous runoff)
! routing of autochtonous runoff  (triangular like this: /\_ ) - used for riverflow generated WITHIN the basin (tL*=0, tR*=tL+tR)
  allocate( hrout_intern(maxval(ceiling(1+sum(prout, dim=2))) ,subasin)) !allocate memory for triangular unit hydrograph (autochtonous runoff)
  hrout=0.
  hrout_intern=0.
 
  !allocate arrays for in- and outflow into/out of subbasins, as their length needs to accomodate hrout, too
  allocate( qout(366 + size(hrout,dim=1), subasin))
  allocate( qin (subasin))
  qout(:,:)=0.
  qin (:)=0.

  if (dosediment) then
    allocate( qsediment2_t(366 + size(hrout,dim=1), nt, subasin))
    allocate( qin_sed (subasin))
    qsediment2_t(:,:,:)=0.
    qin_sed (:)=0.
  end if  


    DO i=1,subasin
        itl = ceiling (prout(i,1)) !index to position in hrout where peak will be located
        itl = max(1,itl) !prevent 0 index
        j   = ceiling   ( prout(i,1)+prout(i,2) ) !                   (integer index to end of triangle)

    if(itl > 1) then
      do ih = 1,(itl-1) !ascending part 
        hrout_intern(ih,i) = (ih-0.5) / prout(i,1)
        hrout       (ih,i) = 0
      end do
      temp2 = hrout_intern(itl-1,i)
    end if 
  
    !peak interval
    if (prout(i,1) /= 0.) then
        temp2 = (itl-1) / prout(i,1) !value AT interval border before itl
        temp4 = 1 -(itl - prout(i,1))  !fraction of rising limb in interval of peak !!r
        hrout_intern(itl,i) =   (temp2 + 1)/2 * temp4
     else   
        temp4 = 0  
        hrout_intern(itl,i) =   1
     end if
     
    temp3 = 1- temp4 - max(0., (itl - (prout(i,1) + prout(i,2)))) !fraction of falling limb in interval of peak 

    if (prout(i,1) /= 0.) then
        temp2 = 1 - temp3 / prout(i,2) !value AT interval border after itl (or end of triangle) 
        hrout       (itl,i) =   (1 + temp2)/2 * temp3 
    else
        hrout       (itl,i) =   1
    end if
    
    hrout_intern(itl,i) =  hrout_intern(itl,i)  + hrout       (itl,i)

    !recession part - do fully covered intervals only
    if (itl+1 < prout(i,1) + prout(i,2)) then
      do ih = (itl+1),(j-1) !recession part 
        temp4 = 1 - (ih-0.5-prout(i,1)) / prout(i,2)
        hrout_intern(ih,i) = temp4  !recession parts of internal and external UHG are identical
        hrout       (ih,i) = temp4
      end do
    end if

    !remaining part of triangle that ends midway in interval
    if (prout(i,2) /= 0.) then
        !temp4 = j - (prout(i,1) + prout(i,2)) !fraction of interval covered by triangle tail
        !strange effect in release mode: 7.63 and 20.37 produce wrong result 
        temp4 =  - (prout(i,1) + prout(i,2)-j) !fraction of interval covered by triangle tail
        
        if (temp4 == 0) temp4=1. !special case: triangle ends exactly at interval border
    
        temp2 = 1 - (j-1-prout(i,1)) / prout(i,2) !value AT interval border before j
         
        temp3 =   (temp2 + 0)/2 * temp4
     
        hrout_intern(j,i) = temp3  !recession parts of internal and external UHG are identical
        hrout       (j,i) = temp3
    end if
            
    hrout(:,i)        = hrout(:,i)        / sum(hrout(:,i))          !normalize response function
    hrout_intern(:,i) = hrout_intern(:,i) / sum(hrout_intern(:,i))   !normalize response function

    if (sum(hrout(:,i))==0 .OR. sum(hrout_intern(:,i))==0 .OR. any(hrout(:,i)<0) .OR. any(hrout_intern(:,i)<0)) then
        write(*,"(A)") "Error when computing response functions: Please send response.dat to the developers."
        stop
    end if

  END DO

  if (f_routing_response) then
      OPEN(11,FILE=pfadn(1:pfadi)//'routing_response.out' ,STATUS='replace')
      WRITE(11,'(a)') 'Output of linear response function'
      WRITE(11,'(a,i0,a)')'Subasin-ID,translation [days], retention [days], uh(1,',size(hrout,dim=1),') [-]'
      DO i=1,subasin
          write(fmtstr,'(a,i0,a)')'(I5,2f10.2,2x,',size(hrout,dim=1),'f5.2)'		!generate format string
	        WRITE (11,fmtstr)  &
            id_subbas_extern(i),prout(i,1),prout(i,2),(hrout(ih,i),ih=1, size(hrout,dim=1))
      END DO

      OPEN(11,FILE=pfadn(1:pfadi)//'routing_response_intern.out' ,STATUS='replace')
      WRITE(11,'(a)') 'Output of linear response function (for routing of autochtonous fluxes)'
      WRITE(11,'(a,i0,a)')'Subasin-ID,translation [days], retention [days], uh(1,',size(hrout,dim=1),') [-]'
      DO i=1,subasin
          write(fmtstr,'(a,i0,a)')'(I5,2f10.2,2x,',size(hrout,dim=1),'f5.2)'		!generate format string
	        WRITE (11,fmtstr)  &
            id_subbas_extern(i),prout(i,1),prout(i,2),(hrout_intern(ih,i),ih=1, size(hrout_intern,dim=1))
      END DO
  end if


! CALL Reservoir Sedimentation and Management Modules
  IF (doreservoir) THEN
    CALL reservoir (0,upstream,idummy)
  END IF

END IF

! ------------------------------------------------------------------------
IF (STATUS == 1) THEN   !beginning of a new simulation year

! initialize ...
! ... and take qout and volact from the last nn days of last year (nn:length of unit hydrograph)

  IF (t > tstart) THEN
      DO id=1,size(hrout,dim=1)-1   !shift routed riverflow that reaches beyond boundary of year to the beginning of new year
        qout(id,1:subasin)  = qout(daylastyear+id,1:subasin)
      END DO
      
      if (dosediment) then
        DO id=1,size(hrout,dim=1)-1   !shift routed sediment flow that reaches beyond boundary of year to the beginning of new year
            qsediment2_t(id,1,1:subasin)  = qsediment2_t(daylastyear+id,1,1:subasin)
        END DO
      end if
      
  END IF

! ... and initialize remaining values
  qout(size(hrout,dim=1):size(qout,dim=1),1:subasin)=0. !reset the rest of the year to 0
  if (dosediment) then
      qsediment2_t(size(hrout,dim=1):size(qout,dim=1),:,1:subasin)=0. !reset the rest of the year to 0
  end if    

! CALL Reservoir Sedimentation and Management Modules
  IF (doreservoir) THEN
    CALL reservoir (1, upstream,idummy)
  END IF

END IF

! ------------------------------------------------------------------------
IF (STATUS == 2) THEN !regular call during timestep

! ..........................................................................
!**  Transfer of water between sub-basins
!    assumption: time delay = 1 day (?)
!    transfer variable qtemp in [m3/day]
!ii: needs to be refined for subdaily resolution
    
  IF (dotrans) THEN
    DO i=1,ntrans
      IF (t >= y_trans(i)) THEN
        qtemp=0.
        qtemp_sed(:)=0.
        irout=trans_start(1,i)
        IF (trans_start(2,i) == 2) THEN ! if transposition starts in river 
          qtemp=MIN(qout(d,irout)*86400.,q_trans(i)*86400.) !take as much water as specified, but no more than is available [m³]
          qout(d,irout)=qout(d,irout)-qtemp/86400.   !reduce start of transposition by amount abstracted
          qout(d,irout)=MAX(qout(d,irout),0.)   !prevent negative values (is this necessary?)
          if (dosediment) then
            !qtemp_sed = qsediment(d,irout) * qtemp/86400. / qout(d,irout)  ![t/d] assume homogenous mixing, abstract the same fraction of sediment as water
            qtemp_sed = sediment_out(irout,:) * qtemp/86400. / qout(d,irout)  ! assume homogenous mixing, abstract the same fraction of sediment as water for every particle size class [t/d]
            
            sediment_out(irout,:) =  sediment_out(irout,:) -  qtemp_sed  ! assume homogenous mixing, abstract the same fraction of sediment as water for every particle size class [t/d]
        
            !qsediment2_t(d,1,irout) = qsediment2_t(d,1,irout) - sum(qtemp_sed)  !abstract sediment from source (total over all particle size classes)
            !qsediment2_t(d,1,irout)=MAX(qsediment2_t(d,1,irout),0.)   !prevent negative values (is this necessary?)
          end if    
!    from acude (there must be a large acude in this muni !)
        ELSE IF (trans_start(2,i) == 1) THEN
          qtemp=MIN(volact(d,irout)*1.e6,q_trans(i)*86400.)
          volact(d,irout)=volact(d,irout)-qtemp/1.e6
          volact(d,irout)=MAX(volact(d,irout),0.)
          if (dosediment) then
            qtemp_sed(:) = 0.  ![t/d] assume no sediment is transfered when water is abstracted from reservoir
          end if
        END IF
        irout=trans_end(1,i)
        IF (trans_end(2,i) == 2) THEN !    into river
          qout(d+1,irout)=qout(d+1,irout)+qtemp*(1.-loss_trans(i))/ 86400.
          if (dosediment) then
              sediment_in(irout,:) = sediment_in(irout,:) + qtemp_sed
              qsediment2_t(d+1,1,irout) = qsediment2_t(d+1,1,irout)+ sum(qtemp_sed)
          end if     
        ELSE IF (trans_end(2,i) == 1) THEN !    into acude (there must be a large acude in this muni !)
          volact(d+1,irout)=volact(d+1,irout)+ qtemp*(1.-loss_trans(i))/1.e6
          if (dosediment) then
            sediment_in(irout,:) = sediment_in(irout,:) + qtemp_sed !particle size specific
            sed_inflow(d+1,irout) = sed_inflow(d+1,irout) + sum(qtemp_sed) !sum over all particle size fractions
          end if  
        END IF
      END IF
    END DO
  END IF

!  water_subbasin(d,i): autochtonous runoff of each sub-basin (after small reservoirs)
  !distribute autochtonous runoff according to "reduced" response function, because it doesn't travel all the way through the subbasin but only half the distance on average
  ! Effectively, the unit hydrograph (hrout) is shrunk by half to account for less translation and retention of the autochtonous runoff compared to the runoff entering from upstream ("ishft" equals division by two, but is faster)

 DO i=1,subasin
     if (do_pre_outflow(i))  then        !if water outflow from current subbasins is NOT given
        qout(d,i)=qout(d,i) + water_subbasin(d,i)
     else
         DO ih=1,size(hrout_intern,dim=1)
             j =  d+ih-1 !index for qout to write to.
             qout(j,i)=qout(j,i) + water_subbasin(d,i)*hrout_intern(ih,i)    !m3/s
             if (dosediment) then
                 qsediment2_t(j,1,i)=qsediment2_t(j,1,i) + sum(sediment_subbasin_t(d,1,i,:))*hrout_intern(ih,i)    !m3/s
             end if
         END Do
     end if
  END DO


!cccccccccccccccccccccccccccc
! MAIN ROUTING LOOP
!cccccccccccccccccccccccccccc
! Calculate qin and qout in order of routing scheme (as was read in and transformed from routing.dat)
  qin=0.  ! set inflow qin = zero for all sub-basin
  if (dosediment) then
    qin_sed=0.  ! set inflow qin = zero for all sub-basin
  end if  

  DO i=1,subasin
    upstream=upbasin(i)  !internal code-ID for most upstream sub-basin (should usually just be i)
    downstream=downbasin(i) !internal code-ID for receiving sub-basin

    if (.NOT. do_pre_outflow(i))  then        !if water outflow from current subbasins is NOT given

! Route inflow from upstream sub-basins (qin) through current sub-basin within nn days
        if (qin(upstream) == -1) then !nodata in (prespecified) outflow of subbasin
            qout(d:(d-1+size(hrout,dim=1)),upstream) = -1 !set entire effected period to "no data"
            cycle
        end if
        
        DO ih=1,size(hrout,dim=1)
          if (qout(d+ih-1,upstream) /= -1) then !only do computation if this timestep is not affected by prior nodata
            qout(d+ih-1,upstream)=qin(upstream)*hrout(ih,upstream)  &
              + qout(d+ih-1,upstream)
          end if 
          
          if (dosediment) then
              if (qsediment2_t(d+ih-1,1,upstream) /= -1) then !only do computation if this timestep is not affected by prior nodata
                qsediment2_t(d+ih-1,1,upstream)=qin_sed(upstream)*hrout(ih,upstream)  &
                  + qsediment2_t(d+ih-1,1,upstream)
               end if 
           end if
              
        if (qout(d+ih-1,upstream) == -1) cycle !no data, don't do further calculations
! Transmission losses by evaporation in river
! river width for given discharge (estimated according to global
! relationship given by Leopold, 1994 (Fig 8.10) [m]
          IF (ih == 1) THEN
            IF (qout(d+ih-1,upstream) > 0.) THEN
              temp2=10.**(LOG10(qout(d+ih-1,upstream))*0.494+1.031)
            ELSE
              temp2=0.
            END IF
! Length of river stretch in sub-basin estimated as diametre of (circular)
! catchment area multiplied by sinuoisity factor of meandering stream (=1.5) [m]
            temp2=temp2*2.*((area(upstream)/3.147)**0.5)*1000.*1.5
            temp2=temp2*pet(d,upstream)/1000.
            temp2=MIN(temp2,qout(d+ih-1,upstream)*86400.)
            qloss(d,upstream)=qloss(d,upstream)+temp2
            qout(d+ih-1,upstream)=qout(d+ih-1,upstream)-temp2/86400.
            qout(d+ih-1,upstream)=MAX(0.,qout(d+ih-1,upstream))
          END IF
        END DO
! end of do-loop for summing up flows within nn days

! Assign outflow as inflow into the next downstream sub-basin
! and add possible inflow from other upstream sub-basins
        IF (doreservoir .AND. qout(d,upstream)/= -1) THEN
          CALL reservoir (2,upstream,idummy)
          qout(d,upstream) = res_qout(d,upstream) ! replace qout with reservoir outflow
        END IF

    END IF ! water outflow from current subbasin given?

	IF (downstream /= 9999 .AND. downstream /= 999) THEN
        if (qout(d,upstream)/= -1 .AND. qin(downstream)/=1 ) then
            qin(downstream)=qin(downstream) + qout(d,upstream) !the downstream subbasin receives what has been computed for the upstream basin
        else
            qin(downstream) = -1 !mark as "no data"
        end if
      if (dosediment) then
        if (qsediment2_t(d,1,upstream)/= -1 .AND. qin_sed(downstream)/=1 ) then
         qin_sed(downstream)=qin_sed(downstream) + qsediment2_t(d,1,upstream) 
        else
         qin_sed(downstream)= -1  !mark as "no data" 
        end if 
      end if
    END IF

  END DO !1,subasin
!  END OF MAIN ROUTING LOOP




END IF


! -----------------------------------------------------------------------
IF (STATUS == 3) THEN


! daily output of water discharge in the river for the entire year
if (f_river_flow) then
  OPEN(11,FILE=pfadn(1:pfadi)//'River_Flow.out',STATUS='old' ,POSITION='append'  )
  write(fmtstr,'(a,i0,a)')'(2i6,',subasin,'f14.3)'		!generate format string
  do j=1, dayyear
	WRITE (11,fmtstr)t, j, (qout(j,i),i=1,subasin)
	!WRITE (11,'(2i6,<subasin>f14.3)')t, j, (qout(j,i),i=1,subasin)
  enddo
  CLOSE (11)
endif

  IF (doreservoir) THEN
    CALL  reservoir (3,idummy,idummy)
  END IF
END IF

RETURN
END SUBROUTINE routing
