SUBROUTINE lake(STATUS,muni)
 
! Code converted using TO_F90 by Alan Miller
! Date: 2005-06-30  Time: 13:47:13
 
use lake_h
use climo_h
use common_h
use hymo_h
use params_h
use time_h
use reservoir_h
use erosion_h
use utils_h

IMPLICIT NONE


INTEGER, INTENT(IN)                      :: STATUS
INTEGER, INTENT(IN)                     :: muni

INTEGER :: imun,i,j,ii,ie,mm,k,dummy1,id,ih,d_laststep,g,loop,istate,ka,idummy2(20),dummy,dummy5(200)
REAL :: xx(5),frac,frac_hrr(5)
REAL :: lakeinflow_r(5),lakeoutflow_r(5),lakeretention_r(5),lakevolume_r(5),dummy9
REAL :: lakesedinflow_r(5),lakesedoutflow_r(5),lakesedretention_r(5),lakesedimentation_r(5)
!  evaporation form acudes
REAL :: evap(5),rain(5)
INTEGER :: overflow_delay
REAL :: temparea,dummy2(5),dummy3(5,200),dummy4,dummy6(5),dummy7
REAL :: help,help1,help2,help3,delta_vol(5),class
!*************************************************************************************
REAL :: totallakeinflow,totallakeoutflow,totallakeprec,totallakeevap,totallakearea,totallakevol
REAL :: totalsedinflow,totalsedoutflow,totalsedimentation,cumsedimentation
REAL :: totalrunoff,directrunoff,wateryield,totalrunoff2,totalarea
CHARACTER(12) :: subarea
REAL :: cumarea,cumvolume,cumrunoff,cumrunoff2,cuminflow
!*************************************************************************************
character(len=1000) :: fmtstr	!string for formatting file output


!REAL :: dummy7(24),dummy8(24,100)
!** -------------------------------------------------------------------
IF (STATUS == 0) THEN
  
write(*,*)doacudyear
! Read of small reservoir characteristics
! Read number of reservoir per classes
  OPEN(11,FILE=pfadp(1:pfadj)// 'Reservoir/lake.dat',STATUS='old')
  READ (11,*)
  DO k=1,5
!    READ(11,*) dummy1,maxlake0(k),maxlake_factor(k),lake_vol0_factor(k),lake_increase(k),alpha_Molle(k),damk_Molle(k),damc_hrr(k),damd_hrr(k)
    READ(11,*) dummy1,maxlake0(k),lake_vol0_factor(k),lake_increase(k),alpha_Molle(k),damk_Molle(k),damc_hrr(k),damd_hrr(k)
  END DO
  CLOSE(11)

! Read maximum fraction of volume of each hypothetical representative reservoir of class k (m3)
  OPEN(11,FILE=pfadp(1:pfadj)// 'Reservoir/lake_maxfrvol.dat',IOSTAT=istate,STATUS='old')
	IF (istate/=0) THEN					!lake_frarea.dat not found
	  write(*,*)pfadp(1:pfadj)// 'Reservoir/lake_maxfrvol.dat was not found. Run the model anyway.'
      DO imun=1,subasin
	    DO k=1,5
		  maxlake_factor(imun,k)=.5
	    ENDDO
	    dummy5(imun)=id_subbas_extern(imun)
	  ENDDO
	ELSE
      READ(11,*)
      DO imun=1,subasin
	    DO k=1,5
		  maxlake_factor(imun,k)=0
	    ENDDO
	    dummy5(imun)=id_subbas_extern(imun)
	  ENDDO
      READ(11,*,IOSTAT=ka)dummy1,(dummy6(k),k=1,5)
      DO imun=1,subasin
        IF (dummy1==id_subbas_extern(imun)) THEN
	      DO k=1,5
		    maxlake_factor(imun,k)=dummy6(k)
		  ENDDO
		  dummy5(imun)=dummy1
	    ENDIF
	  ENDDO
	  DO WHILE (ka==0)
        READ(11,*,IOSTAT=ka) dummy1,(dummy6(k),k=1,5)
	    dummy5(imun)=dummy1
        DO imun=1,subasin
	      IF (dummy1==id_subbas_extern(imun)) THEN
	        DO k=1,5
		      maxlake_factor(imun,k)=dummy6(k)
		    ENDDO
		    dummy5(imun)=dummy1
		  ENDIF
	    ENDDO
	  ENDDO
	ENDIF
    DO imun=1,subasin
!write(*,'(2I4,10F10.4)')id_subbas_extern(imun),dummy5(imun),(maxlake_factor(imun,k),k=1,5)
      IF (dummy5(imun) /= id_subbas_extern(imun)) THEN
        WRITE(*,*) 'Sub-basin-IDs in file lake_frarea.dat must  &
				have the same ordering scheme as in hymo.dat'
        STOP
      END IF
    END DO
  CLOSE(11)
!stop

!   Andreas block begin
!** read number of small reservoirs for each year and subbasin
!   if data are in the file (index in line 3 (doacudyear) is .true.) then
!   the exact values will be used instead of the function 
!   of temporal variation of acudes number (which is used if only small_reservoirs.dat is read)
  acudfloatyear(:,:,:)=-99
  OPEN(11,FILE=pfadp(1:pfadj)// 'Reservoir/lake_year.dat', IOSTAT=istate,STATUS='old')
	IF (istate/=0) THEN					!lake_year.dat not found
	  write(*,*)pfadp(1:pfadj)// 'Reservoir/lake_year.dat was not found. Run the model anyway.'
	  doacudyear=.FALSE.
	ELSE
	  doacudyear=.TRUE.
	  DO t=tstart,tstop
       DO imun=1,subasin
	    DO k=1,5
		  acudfloatyear(imun,k,t-tstart+1)=0
	    ENDDO
	   ENDDO
	  ENDDO
	  READ(11,*)
	  READ(11,*)
	  DO t=tstart,tstop
       DO imun=1,subasin
        READ(11,*) dummy, dummy1,(acudfloatyear(imun,k,t-tstart+1),k=1,5) 
       ENDDO 
      END DO
	ENDIF
  CLOSE(11)
!   Andreas block end

!DO t=tstart,tstop
!  DO imun=1,subasin
!    write(*,*)t,doacudyear,imun,(acudfloatyear(imun,k,t-tstart+1),k=1,5)
!  ENDDO
!ENDDO
!stop

! Read number of reservoir per classes at the beginning of the simulation period
  IF (.NOT. doacudyear) THEN					!lake_year.dat not found
   OPEN(11,FILE=pfadp(1:pfadj)// 'Reservoir/lake_number.dat',STATUS='old')
    READ(11,*)
    DO imun=1,subasin
	  DO k=1,5
		acud(imun,k)=0
	  ENDDO
	  dummy5(imun)=id_subbas_extern(imun)
	ENDDO
    READ(11,*,IOSTAT=ka)dummy1,(dummy6(k),k=1,5)
    DO imun=1,subasin
      IF (dummy1==id_subbas_extern(imun)) THEN
	    DO k=1,5
		  acud(imun,k)=dummy6(k)
		ENDDO
		dummy5(imun)=dummy1
	  ENDIF
	ENDDO
	DO WHILE (ka==0)
      READ(11,*,IOSTAT=ka) dummy1,(dummy6(k),k=1,5)
	  dummy5(imun)=dummy1
      DO imun=1,subasin
	    IF (dummy1==id_subbas_extern(imun)) THEN
	      DO k=1,5
		    acud(imun,k)=dummy6(k)
		  ENDDO
		  dummy5(imun)=dummy1
		ENDIF
	  ENDDO
	ENDDO
    DO imun=1,subasin
!write(*,'(2I4,10F7.2)')id_subbas_extern(imun),dummy5(imun),(acud(imun,k),k=1,5)
      IF (dummy5(imun) /= id_subbas_extern(imun)) THEN
        WRITE(*,*) 'Sub-basin-IDs in file lake_number.dat must  &
				have the same ordering scheme as in hymo.dat'
        STOP
      END IF
    END DO
   CLOSE(11)
  ENDIF

! Read ratio between the runoff contributing area of the reservoir class and the total runoff contributing area of the sub-basin (-)
  OPEN(11,FILE=pfadp(1:pfadj)// 'Reservoir/lake_frarea.dat',IOSTAT=istate,STATUS='old')
	IF (istate/=0) THEN					!lake_frarea.dat not found
	  write(*,*)pfadp(1:pfadj)// 'Reservoir/lake_frarea.dat was not found. Run the model anyway.'
      DO imun=1,subasin
	    DO k=1,5
		  lakefrarea(imun,k)=.166
	    ENDDO
	    dummy5(imun)=id_subbas_extern(imun)
	  ENDDO
	ELSE
      READ(11,*)
      DO imun=1,subasin
	    DO k=1,5
		  lakefrarea(imun,k)=0
	    ENDDO
	    dummy5(imun)=id_subbas_extern(imun)
	  ENDDO
      READ(11,*,IOSTAT=ka)dummy1,(dummy6(k),k=1,5)
      DO imun=1,subasin
        IF (dummy1==id_subbas_extern(imun)) THEN
	      DO k=1,5
		    lakefrarea(imun,k)=dummy6(k)
		  ENDDO
		  dummy5(imun)=dummy1
	    ENDIF
	  ENDDO
	  DO WHILE (ka==0)
        READ(11,*,IOSTAT=ka) dummy1,(dummy6(k),k=1,5)
	    dummy5(imun)=dummy1
        DO imun=1,subasin
	      IF (dummy1==id_subbas_extern(imun)) THEN
	        DO k=1,5
		      lakefrarea(imun,k)=dummy6(k)
		    ENDDO
		    dummy5(imun)=dummy1
		  ENDIF
	    ENDDO
	  ENDDO
	ENDIF
    DO imun=1,subasin
!write(*,'(2I4,10F10.4)')id_subbas_extern(imun),dummy5(imun),(lakefrarea(imun,k),k=1,5)
      IF (dummy5(imun) /= id_subbas_extern(imun)) THEN
        WRITE(*,*) 'Sub-basin-IDs in file lake_frarea.dat must  &
				have the same ordering scheme as in hymo.dat'
        STOP
      END IF
    END DO
  CLOSE(11)

!Distribution of overflow discharges into the classes of larger reservoirs
! - fraction of the total sub-basin area not controlled by small reservoirs (-)
  DO imun=1,subasin
    dummy7=0.
    DO k=1,5
	  dummy7=dummy7+lakefrarea(imun,k)
	ENDDO
	subfrarea(imun)=max(0.,1.-dummy7)
	subfrout(imun)=subfrarea(imun)
!write(*,'(I4,10F10.4)')id_subbas_extern(imun),dummy7,subfrarea(imun)
  ENDDO

  DO imun=1,subasin
    DO k=1,5
	  lakefrout(imun,k)=lakefrarea(imun,k)
	ENDDO
  ENDDO

! - fraction of outflow discharge from the reservoir class j that flows into the reservoir class r of larger storage volume (-)
  DO imun=1,subasin
    DO k=1,4
	  dummy7=0.
	  DO j=1,k
	    dummy7=dummy7+lakefrarea(imun,j)
      ENDDO
	  lakecumfrout(imun,k)=1.-dummy7
!write(*,'(I4,10F10.4)')k,dummy7,lakefrout(imun,k)
	ENDDO
    lakecumfrout(imun,5)=subfrarea(imun)
!write(*,'(I4,10F10.4)')id_subbas_extern(imun),(lakefrout(imun,k),k=1,5)
  ENDDO
!stop

! for small catchments the cascade routing scheme is not valid
  DO imun=1,subasin
    IF (area(imun)<1000.) then
      DO k=1,5
!		lakefrout(imun,k)=0.
!	    lakecumfrout(imun,k)=1.
	  ENDDO
!	  subfrout(imun)=1.
	ENDIF
!write(*,'(I4,F10.2,5F7.4)')id_subbas_extern(imun),area(imun),(lakefrarea(imun,k),k=1,5)
!write(*,'(I4,F10.2,5F7.4)')id_subbas_extern(imun),subfrout(imun),(lakefrout(imun,k),k=1,5)
  ENDDO
	  

!George Initialization of the maximum storage capacity of each reservoir class (m**3)
  DO imun=1,subasin
    maxlake(imun,1:5)=maxlake0(1:5)*maxlake_factor(imun,1:5)
!write(*,*)imun,(maxlake0(k),k=1,5)
!write(*,*)imun,(maxlake_factor(imun,k),k=1,5)
!write(*,*)imun,(maxlake(imun,k),k=1,5)
    maxlakesub0(imun,1:5)=maxlake0(1:5)*maxlake_factor(imun,1:5)
  ENDDO

!  dummy7=0.
!  DO imun=1,subasin
!    DO k=1,5
!      dummy7=dummy7+maxlake(imun,k)*acud(imun,k)
!write(*,*)imun,maxlake(imun,k)*acud(imun,k),dummy7
!	ENDDO
!  ENDDO
!stop

!george  maxlake0(:)=maxlake0(:)*maxlake_factor(:)
!write(*,*)imun,(maxlake0(k),k=1,5)

!** Acudes Initalisierung (focus area or sub-basins)
!George  maxlake(1:5) = (/5.e4,5.e5,2.e6,6.5E6,25.e6/)

!** Modelling unit Initalisierung
  DO imun=1,subasin
    lakewater0(imun,1:5)=lake_vol0_factor(1:5)*maxlake(imun,1:5)
    totalacud(imun)=sum(maxlake(imun,1:5)*acud(imun,1:5))
    acudfraction(imun,1:5)=(maxlake(imun,1:5)*acud(imun,1:5))/totalacud(imun)
!write(*,*)imun,(lakewater0(imun,k),k=1,5)
!write(*,*)imun,totalacud(imun)
!George    laketrend(imun,1:5)=(/0.,0.,0.,0.,0./)
  END DO

!George  Estimation of small reservoirs' areas depending on volume (km**2)
  DO imun=1,subasin
    DO k=1,5
      lakearea(imun,k)=(alpha_Molle(k)*damk_Molle(k)*((lakewater0(imun,k)  &
			/damk_Molle(k))**((alpha_Molle(k)-1.)/alpha_Molle(k))))/1.e6
	ENDDO
!write(*,'(I4,5F15.3)')imun,(lakearea(imun,k)*acud(imun,k)*1.e6,k=1,5)
  ENDDO

  IF (.NOT. doacudyear) THEN
    acudfloat(:,:)=acud(:,:)
  ELSE IF (doacudyear) THEN  !begin block Andreas
   DO imun=1,subasin
    acud(imun,1:5)=nint(acudfloatyear(imun,1:5,1))*1.
   END DO
  END IF !end block Andreas

  DO imun=1,subasin
    outflow_last_hrr(imun,1:5)=0.
    volume_last_hrr(imun,1:5)=0.
	cumseddep(imun)=0.
  ENDDO

  DO imun=1,subasin
    hmax_hrr(imun,1:5)=((maxlake(imun,1:5))/damk_Molle(1:5))**(1./alpha_Molle(1:5))
    lakewater(1,imun,1:5) = lakewater0(imun,1:5)*acud(imun,1:5)
!write(*,'(I4,5F15.3)')imun,(lakewater0(imun,k),k=1,5)
!write(*,'(I4,5F15.3)')imun,(lakewater(1,imun,k),k=1,5)
!write(*,'(I4,5F15.3)')imun,(acud(imun,k),k=1,5)
!dummy7=0.
!do k=1,5
!dummy7=dummy7+lakewater(1,imun,k)
!enddo
!write(*,'(I4,6F12.3)')imun,dummy7,(lakewater(1,imun,k),k=1,5)
  END DO
!stop

  cumsedimentation=0.
  totalrunoff2=0.


!Ge initialization of output files
  OPEN(11,FILE=pfadn(1:pfadi)//'lake_inflow_r.out',STATUS='replace')
	WRITE(11,*)'Year, day, hour, inflow_r(m**3)'
  CLOSE(11)

  OPEN(11,FILE=pfadn(1:pfadi)//'lake_outflow_r.out',STATUS='replace')
	WRITE(11,*)'Year, day, hour, outflow_r(m**3)'
  CLOSE(11)

  OPEN(11,FILE=pfadn(1:pfadi)//'lake_retention_r.out',STATUS='replace')
	WRITE(11,*)'Year, day, hour, retention_r(m**3)'
  CLOSE(11)

  OPEN(11,FILE=pfadn(1:pfadi)//'lake_volume_r.out',STATUS='replace')
	WRITE(11,*)'Year, day, hour, volume_r(m**3)'
  CLOSE(11)

  IF (dosediment) then
   OPEN(11,FILE=pfadn(1:pfadi)//'lake_sedinflow_r.out',STATUS='replace')
	WRITE(11,*)'Year, day, hour, sedinflow_r(m**3)'
   CLOSE(11)

   OPEN(11,FILE=pfadn(1:pfadi)//'lake_sedoutflow_r.out',STATUS='replace')
	WRITE(11,*)'Year, day, hour, sedoutflow_r(m**3)'
   CLOSE(11)

   OPEN(11,FILE=pfadn(1:pfadi)//'lake_sedretention_r.out',STATUS='replace')
	WRITE(11,*)'Year, day, hour, sedretention_r(m**3)'
   CLOSE(11)

   OPEN(11,FILE=pfadn(1:pfadi)//'lake_sedimentation_r.out',STATUS='replace')
	WRITE(11,*)'Year, day, hour, sedretention_r(m**3)'
   CLOSE(11)
  ENDIF

  OPEN(11,FILE=pfadn(1:pfadi)//'lake_watbal.out',STATUS='replace')
	WRITE(11,*)'Year, day, hour, totallakeinflow(m**3), totallakeoutflow(m**3), &
				totallakeprecip(m**3), totallakeevap(m**3), lakevol(m**3)'
  CLOSE(11)

  IF (dosediment) then
   OPEN(11,FILE=pfadn(1:pfadi)//'lake_sedbal.out',STATUS='replace')
	WRITE(11,*)'Year, day, hour, totalsedinflow(ton), totalsedoutflow(ton), &
				totalsedimentation(ton), cumsedimentation(ton)'
   CLOSE(11)
  ENDIF

  OPEN(11,FILE=pfadn(1:pfadi)//'lake_inflow.out',STATUS='replace')
	WRITE(11,*)'Year, day, hour, reservoir_class, lakeinflow(m**3)'
    write(fmtstr,'(a,i0,a)')'(A24,',subasin,'I15)'		!generate format string	    
	WRITE(11,fmtstr)'                  ', (id_subbas_extern(imun),imun=1,subasin)
	!WRITE(11,'(A24,<subasin>I15)')'                  ', (id_subbas_extern(imun),imun=1,subasin)
  CLOSE(11)

  OPEN(11,FILE=pfadn(1:pfadi)//'lake_outflow.out',STATUS='replace')
	WRITE(11,*)'Year, day, hour, reservoir_class, lakeoutflow(m**3)'
	WRITE(11,fmtstr)'                  ', (id_subbas_extern(imun),imun=1,subasin)
	!WRITE(11,'(A24,<subasin>I15)')'                  ', (id_subbas_extern(imun),imun=1,subasin)
  CLOSE(11)

  OPEN(11,FILE=pfadn(1:pfadi)//'lake_volume.out',STATUS='replace')
	WRITE(11,*)'Year, day, hour, reservoir_class, lakevolume(m**3)'
	WRITE(11,fmtstr)'                  ', (id_subbas_extern(imun),imun=1,subasin)
	!WRITE(11,'(A24,<subasin>I15)')'                  ', (id_subbas_extern(imun),imun=1,subasin)
  CLOSE(11)

  OPEN(11,FILE=pfadn(1:pfadi)//'lake_retention.out',STATUS='replace')
	WRITE(11,*)'Year, day, hour, reservoir_class, lakeretention(m**3)'
	WRITE(11,fmtstr)'                  ', (id_subbas_extern(imun),imun=1,subasin)
	!WRITE(11,'(A24,<subasin>I15)')'                  ', (id_subbas_extern(imun),imun=1,subasin)
  CLOSE(11)

  OPEN(11,FILE=pfadn(1:pfadi)//'lake_vollost.out',STATUS='replace')
	WRITE(11,*)'Year, day, hour, reservoir_class, lakevollost(m**3)'
	WRITE(11,fmtstr)'                  ', (id_subbas_extern(imun),imun=1,subasin)
	!WRITE(11,'(A24,<subasin>I15)')'                  ', (id_subbas_extern(imun),imun=1,subasin)
  CLOSE(11)

  IF (dosediment) then
   OPEN(11,FILE=pfadn(1:pfadi)//'lake_sedinflow.out',STATUS='replace')
	WRITE(11,*)'Year, day, hour, reservoir_class, lakesedinflow(ton)'
	WRITE(11,fmtstr)'                  ', (id_subbas_extern(imun),imun=1,subasin)
	!WRITE(11,'(A24,<subasin>I15)')'                  ', (id_subbas_extern(imun),imun=1,subasin)
   CLOSE(11)

   OPEN(11,FILE=pfadn(1:pfadi)//'lake_sedoutflow.out',STATUS='replace')
	WRITE(11,*)'Year, day, hour, reservoir_class, lakesedoutflow(ton)'
	WRITE(11,fmtstr)'                  ', (id_subbas_extern(imun),imun=1,subasin)
	!WRITE(11,'(A24,<subasin>I15)')'                  ', (id_subbas_extern(imun),imun=1,subasin)
   CLOSE(11)

   OPEN(11,FILE=pfadn(1:pfadi)//'lake_sizedistoutflow.out',STATUS='replace')
	WRITE(11,*)'Year, day, hour, sediment size class, lakesizedistoutflow(m**3)'
	WRITE(11,fmtstr)'                  ', (id_subbas_extern(imun),imun=1,subasin)
	!WRITE(11,'(A24,<subasin>I15)')'                  ', (id_subbas_extern(imun),imun=1,subasin)
   CLOSE(11)
  ENDIF

  OPEN(11,FILE=pfadn(1:pfadi)//'lake_checkwatbal.out',STATUS='replace')
	WRITE(11,*)'Year, day, hour, totalrunoff(m**3), directrunoff(m**3), totallakeinflow(m**3), totallakeoutflow(m**3), &
				wateryield(m**3)'
  CLOSE(11)

 
END IF


!** ----------------------------------------------------------------------
IF (STATUS == 1) THEN
  
!   trend of construction of new reservoirs increases storage capacity
  
!*****************************************************************
! OLD VERSION
!  BAZIN (1992) found an exponential increase in the number of
!  acudes in Taua between 1900-1992 of 4.92% / year.
!  calculation of increase based on not rounded numbers (acudfloat)
!  water balance calculations with rounded numbers (acud)
!George  IF (t < 1992 .AND. t /= tstart) THEN
!*****************************************************************

  IF (.NOT. doacudyear) THEN !Andreas
    IF (t /= tstart) THEN
      DO imun=1,subasin
        acudfloat(imun,:)=acudfloat(imun,:)+acudfloat(imun,:)*lake_increase(:)
        acud(imun,:)=nint(acudfloat(imun,:))*1.
	  ENDDO
    END IF
!  if the number of acudes is given lake_year.dat for each year
  ELSE IF (doacudyear) THEN  !begin block Andreas
    DO imun=1,subasin
      acud(imun,1:5)=nint(acudfloatyear(imun,1:5,t-tstart+1))*1.
    END DO
  END IF !end block Andreas

  
!** Initialization for actual year of calculation
  DO imun=1,subasin
    DO id=1,dayyear*nt
	  lakeoutflow(id,imun) = 0.
	  lakesedout(id,imun)=0.
	  DO g=1,n_sed_class
		lakefracsedout(id,imun,g)=0.
	  ENDDO
    END DO
	damareaact(imun)=0. !for subbasin without strategic reservoir
  END DO				!for subbasin with strategic reservoir, surface area is updated in the reservoir module

!  DO imun=1,subasin
!    maxlakewater(imun,1:5)=maxlake(imun,1:5)*acud(imun,1:5)
!  END DO

END IF


!** -----------------------------------------------------------------------
IF (STATUS == 2) THEN
  
!   Lakes are filled with surface flow (lakeinflow) and are emptied
!   by evaporation
  
!  water balance of acudes, water is distributed between volume classes
!  it is assumed that smaller volume class acudes are lying above
!  larger acudes in the basin, i.e. if total storage capacity of
!  one volume class is filled, additional inflow is added to next
!  larger volume class
  
!  modification 06.06.01:
!  this routing scheme may be too strict (percentage of total runoff
!  which is captured by small reservoirs is too high), as not all dams
!  of a smaller volume class will are situated upstream of a dam of the next
!  larger volume class:
!  outflow of one volume class is sub-divided into e.g. five parts (in the
!  case of the smallest class):
!  one fifth goes to each larger volume class and one-fifth is directly runoff
!  from the catchment

!if(d==13)stop

!George 
! routing between landscape units
! currently, simple summing up of runoff of all LUs
! inflow to acudes of lowlands: lakeinflow

! actual maximum storage capacity changed by sediment deposition
  maxlakewater(muni,1:5)=maxlake(muni,1:5)*acud(muni,1:5)
!write(*,'(I4,5F15.2)')muni,(maxlakewater(muni,k),k=1,5)

! check water and sediment inflow discharge into the small reservoirs
!George revised runoff contributing area by subtraction of the small reservoirs' areas
  DO ih=1,nt					
    hour=ih
    step=(d-1)*nt+hour
!**********************************************************************************
!dummy7(ih)=water_subbasin_t(d,ih,muni)/(1.-intercepted)
!**********************************************************************************
!write(*,'(4I5,4F12.3)')muni,t,d,ih,water_subbasin_t(d,ih,muni),water_subbasin(d,muni)

!	lakeinflow(step,muni)=water_subbasin_t(d,ih,muni)*(intercepted/(1.-intercepted))
	lakeinflow(step,muni)=water_subbasin_t(d,ih,muni)
	totalrunoff2=totalrunoff2+lakeinflow(step,muni)

!water_subbasin_t(d,ih,muni)=50.
cumrunoff=cumrunoff+water_subbasin_t(d,ih,muni)
if(muni==subasin .and. ih==nt) then
!write(*,'(2I6,4F15.3)')step,muni,cumrunoff,water_subbasin_t(d,ih,muni)
cumrunoff=0.
endif
!if(muni==subasin .and. ih==nt) write(*,'(2I5,F20.3)')muni,step,totalrunoff2
!write(*,'(4I5,4F12.3)')muni,t,d,ih,lakeinflow(step,muni),water_subbasin_t(d,ih,muni),water_subbasin(d,muni),water_subbasin_t(d,ih,muni)*(intercepted/(1.-intercepted))
    temparea=0.
    DO k=1,5
      temparea=temparea+lakearea(muni,k)*acud(muni,k)
!write(*,'(3I6,4F10.3)')step,muni,k,temparea,lakearea(muni,k),acud(muni,k)
    END DO
    temparea=temparea+damareaact(muni)/1.e6
cumarea=cumarea+temparea
if(muni==subasin .and. ih==nt) then
!write(*,*)step,cumarea*1.e6
cumarea=0.
endif
!write(*,'(2I6,4F10.3)')step,muni,damareaact(muni)

!write(*,'(3I6,3F15.3)')t,d,muni,temparea,damareaact(muni)/1.e6
    temparea=MAX(0.,1.-temparea/area(muni))
!write(*,'(3I6,3F15.3)')t,d,muni,temparea
!George revised water inflow into the small reservoirs 
!    IF (.NOT. dohour) THEN
!      lakeinflow(d,muni)=lakeinflow(d,muni)*temparea
!    ELSE
    lakeinflow(step,muni)=lakeinflow(step,muni)*temparea
 	water_subbasin_t(d,ih,muni)=water_subbasin_t(d,ih,muni)*temparea	!reduce sub-daily runoff
!    ENDIF
!dummy2=water_subbasin(d,muni)
    if (ih==nt) water_subbasin(d,muni)=water_subbasin(d,muni)*temparea			!reduce daily runoff

!write(*,'(I4,4F15.3)')step,id_subbas_extern(muni),(damareaact(muni)/1.e6)+(lakearea(muni,k)*acud(muni,k)),area(muni),water_subbasin(d,muni)

	lakerunoff(step,muni)=water_subbasin_t(d,ih,muni)

cumrunoff2=cumrunoff2+water_subbasin_t(d,ih,muni)
if(muni==subasin .and. ih==nt) then
!write(*,'(2I6,4F15.3)')step,muni,cumrunoff2,water_subbasin_t(d,ih,muni)
cumrunoff2=0.
endif

cuminflow=cuminflow+water_subbasin_t(d,ih,muni)*(1.-subfrarea(muni))
if(muni==subasin .and. ih==nt) then
!write(*,'(2I6,4F15.3)')step,muni,cuminflow,water_subbasin_t(d,ih,muni)*(1.-subfrarea(muni)),subfrarea(muni)
cuminflow=0.
endif

!write(*,'(2I6,4F15.3)')step,

  ENDDO
!write(*,'(3I6,3F15.3)')t,d,muni,temparea,dummy2,water_subbasin(d,muni)

!write(*,'(4I6,4F15.3)')muni,t,d,ih,temparea,lakeinflow(step,muni),water_subbasin(d,muni)
!write(*,'(4I6,5F11.3)')muni,t,d,ih,(lakewater(d_laststep,muni,k),k=1,5)


  IF (dosediment) then
    DO ih=1,nt					
      hour=ih
      step=(d-1)*nt+hour
	  lakesedin(step,muni)=0.
	  DO g=1,n_sed_class
!**********************************************************************************
!dummy8(ih,g)=sediment_subbasin_t(d,ih,muni,g)
!**********************************************************************************
!write(*,'(3I5,4F12.3)')d,muni,g,sediment_subbasin_t(d,ih,muni,g)
!	    sediment_subbasin_t(d,ih,muni,g)=sediment_subbasin_t(d,ih,muni,g)*(1.-intercepted)
	    sediment_subbasin_t(d,ih,muni,g)=sediment_subbasin_t(d,ih,muni,g)
!write(*,'(2I5,4F12.3)')d,muni,sediment_subbasin_t(d,ih,muni,g),intercepted
	    lakesedin(step,muni)=lakesedin(step,muni)+sediment_subbasin_t(d,ih,muni,g)
!write(*,'(2I5,4F12.3)')d,muni,sediment_subbasin_t(d,ih,muni,g),intercepted,lakesedin(step,muni)
	  ENDDO
	  IF (lakesedin(step,muni)>0.) THEN
	    DO g=1,n_sed_class
	      lakefracsedin(step,muni,g)=sediment_subbasin_t(d,ih,muni,g)/lakesedin(step,muni)
!write(*,'(2I5,4F12.3)')d,muni,lakefracsedin(step,muni,g)
	    ENDDO
	  ELSE
	    DO g=1,n_sed_class
	      lakefracsedin(step,muni,g)=0.
!write(*,'(2I5,4F12.3)')d,muni,lakefracsedin(step,muni,g)
	    ENDDO
	  ENDIF

!write(*,'(2I5,4F12.3)')d,muni,lakesedin(step,muni),(intercepted/(1.-intercepted))
	  lakesedin(step,muni)=lakesedin(step,muni)
!	  lakesedin(step,muni)=lakesedin(step,muni)*(intercepted/(1.-intercepted))
!write(*,'(2I5,4F12.3)')d,muni,lakesedin(step,muni)
      temparea=0.
      DO k=1,5
        temparea=temparea+lakearea(muni,k)*acud(muni,k)
!write(*,'(4I6,3F15.3)')t,d,muni,k,temparea,lakearea(muni,k),acud(muni,k)
      END DO
      temparea=temparea+damareaact(muni)/1.e6
!write(*,'(3I6,3F15.3)')t,d,muni,temparea,damareaact(muni)/1.e6
      temparea=MAX(0.,1.-temparea/area(muni))
!George revised water inflow into the small reservoirs 
      lakesedin(step,muni)=lakesedin(step,muni)*temparea
	  DO g=1,n_sed_class
 	    sediment_subbasin_t(d,ih,muni,g)=sediment_subbasin_t(d,ih,muni,g)*temparea	!reduce sub-daily runoff
!write(*,'(3I5,4F12.3)')d,muni,g,sediment_subbasin_t(d,ih,muni,g),lakesedin(step,muni),temparea
        if (ih==nt) then
!		  sediment_subbasin(d,muni,g)=sediment_subbasin(d,muni,g)*(1.-intercepted)
		  sediment_subbasin(d,muni,g)=sediment_subbasin(d,muni,g)
		  sediment_subbasin(d,muni,g)=sediment_subbasin(d,muni,g)*temparea
		endif
      END DO
    ENDDO 
  ENDIF


!George calculation of actual water volume and outflow discharges 
  DO ih=1,nt
    hour=ih
    step=(d-1)*nt+hour
	if (step == 1 .and. t==tstart) d_laststep=1
	if (step == 1 .and. t/=tstart) d_laststep=(daylastyear)*nt
	if (step /= 1) d_laststep=step-1

!**********************************************************************************
!if(step<5)water_subbasin_t(d,ih,muni)=1500000.
!if(step<5)lakeinflow(step,muni)=1500000.			
!if(step<5)sediment_subbasin_t(d,ih,muni,:)=150.
!if(step<5)lakesedin(step,muni)=150.		
!**********************************************************************************

!George water balance calculation, not considering outflow discharges
!George delta_vol(k) is defined as either increment or decrement of inflow discharges into the reservoir class k

!  Estimate area of open water storage depending on volume (km**2)
!  formula by Molle (1989) (for acudes < 2 Mio m**3)
!  alpha,mean=2.7 ; K häufig = 1000.
    DO k=1,5
      IF (acud(muni,k) > 0) THEN
        lakearea(muni,k)=(alpha_Molle(k)*damk_Molle(k)*((lakewater(d_laststep,muni,k)/acud(muni,k)  &
			/damk_Molle(k))**((alpha_Molle(k)-1.)/alpha_Molle(k))))/1.e6
	  ELSE
	    lakearea(muni,k)=0.
      END IF
!write(*,*)muni,k,step,d_laststep,lakearea(muni,k),lakewater(d_laststep,muni,k)
    END DO
  
!  Calculate direct evaporation losses (m**3)
    DO k=1,5
	  IF (dohour) THEN	
        evap(k)=MIN(lakewater(d_laststep,muni,k),(pet(d,muni)/nt)*lakearea(muni,k)*  &
			1.e3*acud(muni,k))
	  ELSE
        evap(k)=MIN(lakewater(d_laststep,muni,k),pet(step,muni)*lakearea(muni,k)*  &
			1.e3*acud(muni,k))
	  ENDIF
    END DO
  
!  Infiltration losses (according to Molle, 1989:
!  mean additional fraction of evaporation of 34%)
!  study was made for acudes up to about 2Mio m**3 storage volume
!  -> losses are currently taken only for smaller classes 1-3
!  ! to be improved !
! 1-5 only for taua

    DO k=1,5
      if (maxlake(muni,k)>2.e6) evap(k)=evap(k)+evap(k)*0.34
    END DO
    lakeevap(step,muni,1:5)=evap(1:5)
  
!  Input by direct rainfall on lake surface (m**3)
    DO k=1,5
	  IF (dohour) THEN	
        rain(k)=preciph(step,muni)*lakearea(muni,k)* 1.e3*acud(muni,k)
	  ELSE
        rain(k)=precip(d,muni)*lakearea(muni,k)* 1.e3*acud(muni,k)
	  ENDIF
    END DO
	lakeprec(step,muni,1:5)=rain(1:5)
  
!  water balance of lake
    delta_vol(1:5)=rain(1:5)-evap(1:5)

    lakeoutflow(step,muni)=0.
	dummy2(1:5)=0.
	overflow_delay=0

    DO k=1,4
	  lakeinflow_hrr(step,muni,k)=dummy2(k)+lakeinflow(step,muni)*lakefrarea(muni,k)
!write(*,'(3I5,6F10.1)')step,muni,k,lakeinflow_hrr(step,muni,k),dummy2(k),lakeinflow(step,muni)*lakefrarea(muni,k)
	  help2=lakewater(d_laststep,muni,k)
	  help1=max(0.,help2+lakeinflow(step,muni)*lakefrarea(muni,k)+delta_vol(k))
	  help=max(0.,lakeinflow(step,muni)*lakefrarea(muni,k)+delta_vol(k))
!write(*,'(5I5,6F10.1)')d,d_laststep,daylastyear,muni,k,help,help1,help2,lakeinflow(step,muni)/5.,delta_vol(k)
!write(*,'(5I5,4F12.3)')d,d_laststep,daylastyear,muni,k,lakewater(d_laststep,muni,k)
!write(*,'(5I5,4F12.3)')muni,t,d,d_laststep,k,lakewater(d_laststep,muni,k),lakeinflow(step,muni)/5.,lakewater(step,muni,k),maxlakewater(muni,k)
!write(*,*)lakewater(step,muni,k),maxlakewater(muni,k)
      IF (help1 > maxlakewater(muni,k)) THEN
!  Andreas' version of redistribution
        if (acud(muni,k) /= 0.) then
		  if (overflow_delay == 1) then
		    call lake_routing(muni,k,help,help2,help3)
		  else
		    lakewater(step,muni,k)=maxlakewater(muni,k)
			lakeoutflow_hrr(step,muni,k)=max(0.,help1-maxlakewater(muni,k))
		    help3=lakeoutflow_hrr(step,muni,k)
		  endif
		endif
        if (acud(muni,k) == 0.) then
		  volume_last_hrr(muni,k)=0.
		  lakewater(step,muni,k)=0.
		  lakeoutflow_hrr(step,muni,k)=lakeinflow_hrr(step,muni,k)
		  help3=lakeoutflow_hrr(step,muni,k)
		endif

!if(muni==2)write(*,'(5I5,5F10.1)')d,d_laststep,daylastyear,muni,k,help,help1,help2,help3,maxlakewater(muni,k)
        DO j=k+1,5
		  frac=help3*(lakefrout(muni,j)/lakecumfrout(muni,k))
          delta_vol(j)=delta_vol(j)+frac
		  dummy2(j)=dummy2(j)+frac
        END DO
        lakeoutflow(step,muni)=lakeoutflow(step,muni)+help3*(subfrout(muni)/lakecumfrout(muni,k))
	  ELSE
!	    lakeoutflow(step,muni)=0.
		lakeoutflow_hrr(step,muni,k)=0.
		volume_last_hrr(muni,k)=0.
		outflow_last_hrr(muni,k)=0.
		lakewater(step,muni,k)=help1
      END IF
!write(*,'(2I4,5F10.1)')muni,k,lakeinflow(step,muni)/5.,help,lakeoutflow(step,muni)
    END DO
	lakeinflow_hrr(step,muni,5)=dummy2(5)+lakeinflow(step,muni)*lakefrarea(muni,5)
	help2=lakewater(d_laststep,muni,5)
	help1=max(0.,help2+lakeinflow(step,muni)*lakefrarea(muni,5)+delta_vol(5))
	help=max(0.,lakeinflow(step,muni)*lakefrarea(muni,5)+delta_vol(5))
    IF (help1 > maxlakewater(muni,5)) THEN
      if (acud(muni,5) /= 0.) then
		if (overflow_delay == 1) then
	      call lake_routing(muni,5,help,help2,help3)
		else
		  lakewater(step,muni,5)=maxlakewater(muni,5)
		  lakeoutflow_hrr(step,muni,5)=max(0.,help1-maxlakewater(muni,5))
		  help3=lakeoutflow_hrr(step,muni,5)
		endif
	  endif
      if (acud(muni,5) == 0.) then
		volume_last_hrr(muni,5)=0.
		lakewater(step,muni,5)=0.
		lakeoutflow_hrr(step,muni,5)=lakeinflow_hrr(step,muni,5)
		help3=lakeoutflow_hrr(step,muni,5)
	  endif

      lakeoutflow(step,muni)=lakeoutflow(step,muni)+help3
    ELSE
!      lakeoutflow(step,muni)=0.
	  lakeoutflow_hrr(step,muni,5)=0.
	  volume_last_hrr(muni,5)=0.
	  outflow_last_hrr(muni,5)=0.
	  lakewater(step,muni,5)=help1
    END IF
!k=5
!write(*,'(2I4,5F10.1)')muni,k,lakeinflow(step,muni)/5.,help,lakeoutflow(step,muni)
!write(*,'(4I5,5F11.2)')muni,t,d,d_laststep,(lakewater(step,muni,k),k=1,5)
!if (step==10) stop
!write(*,'(4I5,5F11.2)')muni,t,d,d_laststep,(lakewater(d_laststep,muni,k),k=1,5)
!do k=1,5
!write(*,'(4I5,5F11.2)')muni,t,d,k,lakewater(step,muni,k),maxlakewater(muni,k)
!enddo

    DO k=1,5
	  if (acud(muni,k)>0.) then
	    lakeinflow_hrr(step,muni,k)=lakeinflow_hrr(step,muni,k)/acud(muni,k)
	    lakeoutflow_hrr(step,muni,k)=lakeoutflow_hrr(step,muni,k)/acud(muni,k)
	  else
	    lakeinflow_hrr(step,muni,k)=0.
		lakeoutflow_hrr(step,muni,k)=0.
	  endif
	  lakeretention_hrr(step,muni,k)=max(0.,lakeinflow_hrr(step,muni,k)-lakeoutflow_hrr(step,muni,k))
!write(*,'(3I5,6F10.1)')step,muni,k,acud(muni,k),lakewater(step,muni,k),lakeinflow(step,muni),lakeinflow_hrr(step,muni,k),lakeoutflow_hrr(step,muni,k),lakeoutflow(step,muni)
	ENDDO
!if(d==1 .and. muni==5)write(*,'(3I5,6F10.1)')step,muni,k,(acud(muni,k),k=1,5)
!if(muni==5)write(*,'(3I5,6F10.1)')step,muni,k,(lakewater(step,muni,k),k=1,5)
!if(muni==5)write(*,'(3I5,6F10.1)')step,muni,k,(lakeinflow_hrr(step,muni,k),k=1,5)
!if(muni==5)write(*,'(3I5,6F10.1)')step,muni,k,(lakeoutflow_hrr(step,muni,k),k=1,5)
!stop

  
!  water inflow into the small reservoirs (m**3)
    lakeinflow(step,muni)=lakeinflow(step,muni)*(1.-subfrarea(muni))


    dummy7=0.
    DO k=1,5
!write(*,'(3I5,6F10.1)')step,muni,k,acud(muni,k),lakewater(step,muni,k),maxlakewater(muni,k)
	  dummy7=dummy7+lakeinflow_hrr(step,muni,k)*acud(muni,k)
!	  dummy7=dummy7+lakeoutflow_hrr(step,muni,k)*acud(muni,k)
	ENDDO
!write(*,'(3I5,6F10.1)')step,muni,k,dummy7,lakeinflow(step,muni)
!write(*,'(3I5,6F10.1)')step,muni,k,dummy7,lakeinflow(step,muni)
!if (step==43)stop

!  total retention of incoming runoff in all acudes in this timestep (m**3)
    muniret(step,muni)=MAX(0.,(lakeinflow(step,muni)-  &
		lakeoutflow(step,muni)))
  
!  Estimate area of open water storage depending on volume (km**2)
!  formula by Molle (1989) (for acudes < 2 Mio m**3)
!  alpha,mean=2.7 ; K häufig = 1000.
    DO k=1,5
      IF (acud(muni,k) > 0) THEN
        lakearea(muni,k)=(alpha_Molle(k)*damk_Molle(k)*((lakewater(step,muni,k)/acud(muni,k)  &
			/damk_Molle(k))**((alpha_Molle(k)-1.)/alpha_Molle(k))))/1.e6
	  ELSE
	    lakearea(muni,k)=0.
      END IF
    END DO
  

!  total amount of water storaged in sub-basin (m**3)
    laketot(step,muni)=sum(lakewater(step,muni,1:5))

!if(muni==2)write(*,'(I5,6F10.1)')muni,lakeinflow(step,muni),lakeoutflow(step,muni),laketot(step,muni)

!if (d==10)stop
  
!  water storage relative to storage capacity
    DO k=1,5
      IF (maxlakewater(muni,k) > 0.) THEN
        xx(k)=lakewater(step,muni,k)/maxlakewater(muni,k)
      END IF
    END DO
    IF (sum(maxlakewater(muni,1:5)) > 0.) THEN
      laketotfrac(step,muni)=laketot(step,muni)/ sum(maxlakewater(muni,1:5))
    END IF

! calculation of storage volume in the hypothetical representative reservoir from each class (m**3)
    DO k=1,5
	  if (acud(muni,k)>0.) then
	    lakewater_hrr(step,muni,k)=lakewater(step,muni,k)/acud(muni,k)
	  else
	    lakewater_hrr(step,muni,k)=0.
        lakewater(step,muni,k)=0.
	  endif
!write(*,'(3I5,5F12.3)')step,muni,k,acud(muni,k),lakewater(step,muni,k),lakewater_hrr(step,muni,k),lakewater(d_laststep,muni,k)
!write(*,'(3I5,6F10.1)')step,muni,k,acud(muni,k),lakewater(d_laststep,muni,k),lakewater(step,muni,k),lakewater_hrr(step,muni,k),lakeinflow(step,muni),maxlakewater(muni,k)
	ENDDO

!if(t==2002)read(*,*)

	
!	IF (muni==subasin .and. ih==nt) THEN
!      OPEN(11,FILE=pfadn(1:pfadi)//'lake_watbal.out',STATUS='old',  &
!		 POSITION='append')
!	    WRITE(11,'(3I6,6f15.3)')t,d,hour,totallakeinflow,totallakeoutflow,totallakeprecip,&
!			totallakeevap,totallakevol
!     CLOSE(11)

!	  totallakeinflow=0.
!	  totallakeoutflow=0.
!	  totallakeprecip=0.
!	  totallakeevap=0.
!	  totallakevol=0.
!	ENDIF

    IF (dosediment) then
	  lakesedin_hrr(step,muni,1:5)=0.
      lakesedout(step,muni)=0.
	  dummy2(1:5)=0.
	  DO g=1,n_sed_class
	    dummy3(1:5,g)=0.
	  ENDDO
      DO k=1,4
	    lakesedin_hrr(step,muni,k)=dummy2(k)+lakesedin(step,muni)*lakefrarea(muni,k)
	    DO g=1,n_sed_class
	      IF (lakesedin_hrr(step,muni,k) > 0.) THEN
		    lakefracsedin_hrr(step,muni,k,g)=(dummy3(k,g)+ &
				((lakesedin(step,muni)*lakefrarea(muni,k))*lakefracsedin(step,muni,g)))/ &
		        lakesedin_hrr(step,muni,k)
		  ELSE
		    lakefracsedin_hrr(step,muni,k,g)=0.
		  ENDIF
	    ENDDO
        IF (lakeoutflow_hrr(step,muni,k) > 0. .and. acud(muni,k) /= 0.) THEN
		  lakesedin_hrr(step,muni,k)=lakesedin_hrr(step,muni,k)/acud(muni,k)
		  call sedbal_lake(muni,k)
		  lakesedin_hrr(step,muni,k)=lakesedin_hrr(step,muni,k)*acud(muni,k)
		  lakesedout_hrr(step,muni,k)=lakesedout_hrr(step,muni,k)*acud(muni,k)
          DO j=k+1,5
            frac=lakesedout_hrr(step,muni,k)*(lakefrout(muni,j)/lakecumfrout(muni,k))
		    dummy2(j)=dummy2(j)+frac
	        DO g=1,n_sed_class
		      frac_hrr(g)=(lakefracsedout_hrr(step,muni,k,g)*lakesedout_hrr(step,muni,k))*(lakefrout(muni,j)/lakecumfrout(muni,k))
			  dummy3(j,g)=dummy3(j,g)+frac_hrr(g)
            END DO		    
          END DO
		  dummy4=lakesedout(step,muni)
          lakesedout(step,muni)=lakesedout(step,muni)+(lakesedout_hrr(step,muni,k)*(subfrout(muni)/lakecumfrout(muni,k)))
	      DO g=1,n_sed_class
	        IF (lakesedout(step,muni) > 0.) THEN
		      lakefracsedout(step,muni,g)=((dummy4*lakefracsedout(step,muni,g))+ &
					frac_hrr(g))/lakesedout(step,muni)
		    ELSE
		      lakefracsedout(step,muni,g)=0.
		    ENDIF
		  ENDDO
	    ELSE IF (acud(muni,k) == 0.) THEN
		  lakesedout_hrr(step,muni,k)=lakesedin_hrr(step,muni,k)
	      DO g=1,n_sed_class
		    lakefracsedout_hrr(step,muni,k,g)=lakefracsedin_hrr(step,muni,k,g)
		  ENDDO
          DO j=k+1,5
            frac=(lakesedout_hrr(step,muni,k))*(lakefrout(muni,j)/lakecumfrout(muni,k))
		    dummy2(j)=dummy2(j)+frac
	        DO g=1,n_sed_class
		      frac_hrr(g)=(lakefracsedout_hrr(step,muni,k,g)*lakesedout_hrr(step,muni,k))*(lakefrout(muni,j)/lakecumfrout(muni,k))
			  dummy3(j,g)=dummy3(j,g)+frac_hrr(g)
            END DO		    
          END DO
		  dummy4=lakesedout(step,muni)
          lakesedout(step,muni)=lakesedout(step,muni)+(lakesedout_hrr(step,muni,k)*(subfrout(muni)/lakecumfrout(muni,k)))
	      DO g=1,n_sed_class
	        lakefracsedout_hrr(step,muni,k,g)=lakefracsedin_hrr(step,muni,k,g)
	        IF (lakesedout(step,muni) > 0.) THEN
		      lakefracsedout(step,muni,g)=((dummy4*lakefracsedout(step,muni,g))+ &
					frac_hrr(g))/lakesedout(step,muni)
		    ELSE
		      lakefracsedout(step,muni,g)=0.
		    ENDIF
	      ENDDO
		  lake_vollost(step,muni,k)=0.
	    ELSE
		  lakesedout_hrr(step,muni,k)=0.
	      DO g=1,n_sed_class
	        lakefracsedout_hrr(step,muni,k,g)=0.
	      ENDDO
		  lake_vollost(step,muni,k)=(lakesedin_hrr(step,muni,k)/acud(muni,k))/1.5
        END IF
      END DO
	  lakesedin_hrr(step,muni,5)=dummy2(5)+lakesedin(step,muni)*lakefrarea(muni,5)
	  DO g=1,n_sed_class
	    IF (lakesedin_hrr(step,muni,5) > 0.) THEN
	      lakefracsedin_hrr(step,muni,5,g)=(dummy3(5,g)+ &
				((lakesedin(step,muni)*lakefrarea(muni,5))*lakefracsedin(step,muni,g)))/ &
		        lakesedin_hrr(step,muni,5)
		ELSE
		  lakefracsedin_hrr(step,muni,5,g)=0.
		ENDIF
	  ENDDO
      IF (lakeoutflow_hrr(step,muni,5) > 0. .and. acud(muni,5) /= 0.) THEN
		lakesedin_hrr(step,muni,5)=lakesedin_hrr(step,muni,5)/acud(muni,5)
		call sedbal_lake(muni,5)
		lakesedin_hrr(step,muni,5)=lakesedin_hrr(step,muni,5)*acud(muni,5)
		lakesedout_hrr(step,muni,5)=lakesedout_hrr(step,muni,5)*acud(muni,5)
		dummy4=lakesedout(step,muni)
		lakesedout(step,muni)=lakesedout(step,muni)+lakesedout_hrr(step,muni,5)
		DO g=1,n_sed_class
	      IF (lakesedout(step,muni) > 0.) THEN
		    lakefracsedout(step,muni,g)=((dummy4*lakefracsedout(step,muni,g))+ &
					(lakesedout_hrr(step,muni,5)*lakefracsedout_hrr(step,muni,k,g)))/ &
					lakesedout(step,muni)
		  ELSE
		    lakefracsedout(step,muni,g)=0.
		  ENDIF
		ENDDO
	  ELSE IF (acud(muni,5) == 0.) THEN
		lakesedout_hrr(step,muni,5)=lakesedin_hrr(step,muni,5)
		dummy4=lakesedout(step,muni)
		lakesedout(step,muni)=lakesedout(step,muni)+lakesedout_hrr(step,muni,5)
	    DO g=1,n_sed_class
	      lakefracsedout_hrr(step,muni,5,g)=lakefracsedin_hrr(step,muni,5,g)
	      IF (lakesedout(step,muni) > 0.) THEN
		    lakefracsedout(step,muni,g)=((dummy4*lakefracsedout(step,muni,g))+ &
					(lakesedout_hrr(step,muni,5)*lakefracsedout_hrr(step,muni,k,g)))/ &
					lakesedout(step,muni)
		  ELSE
		    lakefracsedout(step,muni,g)=0.
		  ENDIF
	    ENDDO
		lake_vollost(step,muni,5)=0.
	  ELSE
		lakesedout_hrr(step,muni,5)=0.
	    DO g=1,n_sed_class
	      lakefracsedout_hrr(step,muni,5,g)=0.
	    ENDDO
		lake_vollost(step,muni,5)=(lakesedin_hrr(step,muni,5)/acud(muni,5))/1.5
	  ENDIF
      DO k=1,5
	    maxlake(muni,k)=max(0.,maxlake(muni,k)-lake_vollost(step,muni,k))
!write(*,'(3I4,5F15.3)')step,muni,k,lake_vollost(step,muni,k),maxlake(muni,k)
	  ENDDO
	ENDIF

    DO k=1,5
      maxstorcap_hrr(step,muni,k)=maxlake(muni,k)
	ENDDO
!if(muni==2)write(*,'(I5,6F10.1)')muni,cumseddep(muni)

!    DO k=1,5
!      cumseddep(muni)=cumseddep(muni)+(lake_vollost(step,muni,k)*acud(muni,k)*1.5)
!	ENDDO

!if(muni==2)write(*,'(I5,6F10.1)')muni,lakesedin(step,muni),lakesedout(step,muni),sum(lake_vollost(step,muni,1:5)*acud(muni,1:5)*1.5),cumseddep(muni)
	
!DO k=1,5
!write(*,'(3I5,5F12.3)')step,muni,k,acud(muni,k),lakesedin_hrr(step,muni,k),lakesedout_hrr(step,muni,k),lakesedin(step,muni),lakesedout(step,muni)
!DO g=1,n_sed_class		    
!if(g==1)write(*,'(4I4,4F10.3)')step,muni,k,g,frac_hrr(g),dummy3(k,g),lakesedin_hrr(step,muni,k)*lakefracsedout_hrr(step,muni,k,g),lakesedin_hrr(step,muni,k)
!if(g==1)write(*,'(4I4,5F10.3)')step,muni,k,g,acud(muni,k),frac_hrr(g),dummy3(k,g),lakesedin_hrr(step,muni,k)*lakefracsedin_hrr(step,muni,k,g),lakesedout_hrr(step,muni,k)*lakefracsedout_hrr(step,muni,k,g)
!if(g==1)write(*,'(3I5,5F12.3)')step,muni,k,acud(muni,k),lakefracsedin_hrr(step,muni,k,g),lakefracsedout_hrr(step,muni,k,g)
!if(g==1)write(*,'(3I5,5F12.3)')step,muni,k,acud(muni,k),dummy3(k,g),lakesedin(step,muni)/5.,lakefracsedin(step,muni,g),lakesedin_hrr(step,muni,k)
!if (lakesedin_hrr(step,muni,k)*lakefracsedin_hrr(step,muni,k,g) < lakesedout_hrr(step,muni,k)*lakefracsedout_hrr(step,muni,k,g)) stop
!END DO
!ENDDO


    IF (dosediment) then
     DO k=1,5
	  if (acud(muni,k)>0.) then
	    lakesedin_hrr(step,muni,k)=lakesedin_hrr(step,muni,k)/acud(muni,k)
	    lakesedout_hrr(step,muni,k)=lakesedout_hrr(step,muni,k)/acud(muni,k)
	  else
	    lakesedin_hrr(step,muni,k)=0.
		lakesedout_hrr(step,muni,k)=0.
	  endif
	 ENDDO
	ENDIF


! The lake module returns values of water outflow and sediment outflow of each sub-basin
! after the passage of the cascade routing scheme
    totalacud(muni)=sum(maxlake(muni,1:5)*acud(muni,1:5))
    acudfraction(muni,1:5)=(maxlake(muni,1:5)*acud(muni,1:5))/totalacud(muni)

!	water_subbasin_t(d,ih,muni)=water_subbasin_t(d,ih,muni)+ &
!			lakeoutflow(step,muni)
	water_subbasin_t(d,ih,muni)=water_subbasin_t(d,ih,muni)*subfrarea(muni)+ &
			lakeoutflow(step,muni)
!*******************************************************
!George temporary (discuss with Eva and Till)
!	water_subbasin(d,muni)=water_subbasin(d,muni)+lakeoutflow(step,muni)
	water_subbasin(d,muni)=water_subbasin(d,muni)*subfrarea(muni)+lakeoutflow(step,muni)
!*******************************************************


    IF (dosediment) then
	 DO g=1,n_sed_class
!	  sediment_subbasin_t(d,ih,muni,g)=sediment_subbasin_t(d,ih,muni,g)+lakesedout(step,muni)*lakefracsedout(step,muni,g)
	  sediment_subbasin_t(d,ih,muni,g)=sediment_subbasin_t(d,ih,muni,g)*subfrarea(muni)+lakesedout(step,muni)*lakefracsedout(step,muni,g)
!*******************************************************
!George temporary (discuss with Eva and Till)
	  sediment_subbasin(d,muni,g)=sediment_subbasin(d,muni,g)*subfrarea(muni)+lakesedout(step,muni)*lakefracsedout(step,muni,g)
	 ENDDO
	ENDIF
!*******************************************************

!write(*,'(2I4,4F15.6)')step,muni,dummy7(ih),(dummy8(ih,g),g=1,n_sed_class)
!write(*,'(2I4,4F15.6)')step,muni,water_subbasin_t(d,ih,muni),(sediment_subbasin_t(d,ih,muni,g),g=1,n_sed_class)
!write(*,'(2I4,4F15.6)')step,muni,lakesedin(step,muni),lakesedout(step,muni)

!  sediment inflow into the small reservoirs (m**3)
    IF (dosediment) then
	 lakesedin(step,muni)=lakesedin(step,muni)*(1.-subfrarea(muni))
     DO k=1,5
	  lake_cumseddep(step,muni,k)=max(0.,1.5*((maxlakesub0(muni,k)*acud(muni,k))-maxlakewater(muni,k)))
	 ENDDO
	ENDIF


!	totalsedinflow=totalsedinflow+lakesedin(step,muni)
!	totalsedoutflow=totalsedoutflow+lakesedout(step,muni)
!    DO k=1,5
!	  totalsedimentation=totalsedimentation+(lake_vollost(step,muni,k)*acud(muni,k)*1.5)
!	  cumsedimentation=cumsedimentation+(lake_vollost(step,muni,k)*acud(muni,k)*1.5)
!	ENDDO

	
!	IF (muni==subasin .and. ih==nt) THEN
!      OPEN(11,FILE=pfadn(1:pfadi)//'lake_sedbal.out',STATUS='old',  &
!		 POSITION='append')
!	    WRITE(11,'(3I6,4f15.3)')t,d,hour,totalsedinflow,totalsedoutflow,totalsedimentation,cumsedimentation
!      CLOSE(11)
!	  totalsedinflow=0.
!	  totalsedoutflow=0.
!	  totalsedimentation=0.
!	ENDIF

  ENDDO
END IF

!--------------------------------------------------------------------
IF (STATUS == 3) THEN
  
!** Obtain monthly and annual values
  
!        call average(laketot,annlaketot,monlaketot,1.)
!        call average(lakewater(:,:,1),annlaketot1,monlaketot1,1.)
!        call average(lakewater(:,:,2),annlaketot2,monlaketot2,1.)
!        call average(lakewater(:,:,3),annlaketot3,monlaketot3,1.)
!        call average(lakewater(:,:,4),annlaketot4,monlaketot4,1.)
!        call average(lakewater(:,:,5),annlaketot5,monlaketot5,1.)
!        call total  (lakeevap(:,:,1),annlakeevap1,monlakeevap1,1.)
!        call total  (lakeevap(:,:,2),annlakeevap2,monlakeevap2,1.)
!        call total  (lakeevap(:,:,3),annlakeevap3,monlakeevap3,1.)
!        call total  (lakeevap(:,:,4),annlakeevap4,monlakeevap4,1.)
!        call total  (lakeevap(:,:,5),annlakeevap5,monlakeevap5,1.)
!        call average(laketotfrac,annlakefrac,monlakefrac,1.)
  
!Eva output files deleted

    OPEN(11,FILE=pfadn(1:pfadi)//'lake_checkwatbal.out',STATUS='old',  &
		 POSITION='append')
	  DO d=1,dayyear
	    DO ih=1,nt
		  hour=ih
          step=(d-1)*nt+hour
		  totalrunoff=0.
		  totallakeinflow=0.
		  totallakeoutflow=0.
		  directrunoff=0.
		  DO imun=1,subasin
		    totalrunoff=totalrunoff+lakerunoff(step,imun)
			directrunoff=directrunoff+lakerunoff(step,imun)*subfrarea(imun)
		    totallakeinflow=totallakeinflow+lakeinflow(step,imun)
		    totallakeoutflow=totallakeoutflow+lakeoutflow(step,imun)
		  ENDDO
		  wateryield=totallakeoutflow+directrunoff
	      WRITE(11,'(3I6,6f15.3)')t,d,hour,totalrunoff,directrunoff,totallakeinflow,totallakeoutflow,wateryield
		ENDDO
	  ENDDO
    CLOSE(11)

    OPEN(11,FILE=pfadn(1:pfadi)//'lake_inflow_r.out',STATUS='old',  &
		 POSITION='append')
	  DO d=1,dayyear
	    DO ih=1,nt
		  hour=ih
          step=(d-1)*nt+hour
		  DO k=1,5
		    lakeinflow_r(k)=0.
		    DO imun=1,subasin
              lakeinflow_r(k)=lakeinflow_r(k)+(lakeinflow_hrr(step,imun,k)*acud(imun,k))
			ENDDO
		  ENDDO
	      WRITE(11,'(3I6,5f15.3)')t,d,hour,(lakeinflow_r(k),k=1,5)
		ENDDO
	  ENDDO
    CLOSE(11)

    OPEN(11,FILE=pfadn(1:pfadi)//'lake_outflow_r.out',STATUS='old',  &
		 POSITION='append')
	  DO d=1,dayyear
	    DO ih=1,nt
		  hour=ih
          step=(d-1)*nt+hour
		  DO k=1,5
		    lakeoutflow_r(k)=0.
		    DO imun=1,subasin
              lakeoutflow_r(k)=lakeoutflow_r(k)+(lakeoutflow_hrr(step,imun,k)*acud(imun,k))
			ENDDO
		  ENDDO
	      WRITE(11,'(3I6,5f15.3)')t,d,hour,(lakeoutflow_r(k),k=1,5)
		ENDDO
	  ENDDO
    CLOSE(11)

    OPEN(11,FILE=pfadn(1:pfadi)//'lake_retention_r.out',STATUS='old',  &
		 POSITION='append')
	  DO d=1,dayyear
	    DO ih=1,nt
		  hour=ih
          step=(d-1)*nt+hour
		  DO k=1,5
		    dummy9=0.
		    DO imun=1,subasin
			  dummy9=dummy9+max(0.,(lakeinflow_hrr(step,imun,k)-lakeoutflow_hrr(step,imun,k))*acud(imun,k))
			ENDDO
            lakeretention_r(k)=dummy9
		  ENDDO
	      WRITE(11,'(3I6,5f15.3)')t,d,hour,(lakeretention_r(k),k=1,5)
		ENDDO
	  ENDDO
    CLOSE(11)

    OPEN(11,FILE=pfadn(1:pfadi)//'lake_volume_r.out',STATUS='old',  &
		 POSITION='append')
	  DO d=1,dayyear
	    DO ih=1,nt
		  hour=ih
          step=(d-1)*nt+hour
		  DO k=1,5
		    lakevolume_r(k)=0.
		    DO imun=1,subasin
              lakevolume_r(k)=lakevolume_r(k)+lakewater(step,imun,k)
			ENDDO
		  ENDDO
	      WRITE(11,'(3I6,5f15.3)')t,d,hour,(lakevolume_r(k),k=1,5)
		ENDDO
	  ENDDO
    CLOSE(11)

    IF (dosediment) then
     OPEN(11,FILE=pfadn(1:pfadi)//'lake_sedinflow_r.out',STATUS='old',  &
		 POSITION='append')
	  DO d=1,dayyear
	    DO ih=1,nt
		  hour=ih
          step=(d-1)*nt+hour
		  DO k=1,5
		    lakesedinflow_r(k)=0.
		    DO imun=1,subasin
              lakesedinflow_r(k)=lakesedinflow_r(k)+(lakesedin_hrr(step,imun,k)*acud(imun,k))
			ENDDO
		  ENDDO
	      WRITE(11,'(3I6,5f15.3)')t,d,hour,(lakesedinflow_r(k),k=1,5)
		ENDDO
	  ENDDO
     CLOSE(11)

     OPEN(11,FILE=pfadn(1:pfadi)//'lake_sedoutflow_r.out',STATUS='old',  &
		 POSITION='append')
	  DO d=1,dayyear
	    DO ih=1,nt
		  hour=ih
          step=(d-1)*nt+hour
		  DO k=1,5
		    lakesedoutflow_r(k)=0.
		    DO imun=1,subasin
              lakesedoutflow_r(k)=lakesedoutflow_r(k)+(lakesedout_hrr(step,imun,k)*acud(imun,k))
			ENDDO
		  ENDDO
	      WRITE(11,'(3I6,5f15.3)')t,d,hour,(lakesedoutflow_r(k),k=1,5)
		ENDDO
	  ENDDO
     CLOSE(11)

     OPEN(11,FILE=pfadn(1:pfadi)//'lake_sedretention_r.out',STATUS='old',  &
		 POSITION='append')
	  DO d=1,dayyear
	    DO ih=1,nt
		  hour=ih
          step=(d-1)*nt+hour
		  DO k=1,5
		    dummy9=0.
		    DO imun=1,subasin
			  dummy9=dummy9+max(0.,(lakesedin_hrr(step,imun,k)-lakesedout_hrr(step,imun,k))*acud(imun,k))
			ENDDO
            lakesedretention_r(k)=dummy9
		  ENDDO
	      WRITE(11,'(3I6,5f15.3)')t,d,hour,(lakesedretention_r(k),k=1,5)
		ENDDO
	  ENDDO
     CLOSE(11)

     OPEN(11,FILE=pfadn(1:pfadi)//'lake_sedimentation_r.out',STATUS='old',  &
		 POSITION='append')
	  DO d=1,dayyear
	    DO ih=1,nt
		  hour=ih
          step=(d-1)*nt+hour
		  DO k=1,5
		    lakesedimentation_r(k)=0.
		    DO imun=1,subasin
              lakesedimentation_r(k)=lakesedimentation_r(k)+lake_cumseddep(step,imun,k)
			ENDDO
		  ENDDO
	      WRITE(11,'(3I6,5f15.3)')t,d,hour,(lakesedimentation_r(k),k=1,5)
		ENDDO
	  ENDDO
     CLOSE(11)
	ENDIF

    OPEN(11,FILE=pfadn(1:pfadi)//'lake_watbal.out',STATUS='old',  &
		 POSITION='append')
	  DO d=1,dayyear
	    DO ih=1,nt
		  hour=ih
          step=(d-1)*nt+hour
		  totallakevol=0.
		  totallakeinflow=0.
		  totallakeoutflow=0.
		  totallakeprec=0.
		  totallakeevap=0.
		  DO imun=1,subasin
		    totallakevol=totallakevol+laketot(step,imun)
		    totallakeinflow=totallakeinflow+lakeinflow(step,imun)
		    totallakeoutflow=totallakeoutflow+lakeoutflow(step,imun)
		    DO k=1,5
			  totallakeprec=totallakeprec+lakeprec(step,imun,k)
			  totallakeevap=totallakeevap+lakeevap(step,imun,k)
			ENDDO
		  ENDDO
	      WRITE(11,'(3I6,6f15.3)')t,d,hour,totallakeinflow,totallakeoutflow,totallakeprec,&
			totallakeevap,totallakevol
		ENDDO
	  ENDDO
    CLOSE(11)

    IF (dosediment) then
     OPEN(11,FILE=pfadn(1:pfadi)//'lake_sedbal.out',STATUS='old',  &
		 POSITION='append')
	  DO d=1,dayyear
	    DO ih=1,nt
		  hour=ih
          step=(d-1)*nt+hour
		  totalsedinflow=0.
		  totalsedoutflow=0.
		  totalsedimentation=0.
		  DO imun=1,subasin
		    totalsedinflow=totalsedinflow+lakesedin(step,imun)
		    totalsedoutflow=totalsedoutflow+lakesedout(step,imun)
		    DO k=1,5
			  totalsedimentation=totalsedimentation+(lake_vollost(step,imun,k)*acud(imun,k)*1.5)
			  cumsedimentation=cumsedimentation+(lake_vollost(step,imun,k)*acud(imun,k)*1.5)
			ENDDO
		  ENDDO
		  WRITE(11,'(3I6,4f15.3)')t,d,hour,totalsedinflow,totalsedoutflow,totalsedimentation,cumsedimentation
		ENDDO
	  ENDDO
     CLOSE(11)
	ENDIF

    OPEN(11,FILE=pfadn(1:pfadi)//'lake_inflow.out',STATUS='old',  &
		 POSITION='append')
      write(fmtstr,'(a,i0,a)')'(4I6,',subasin,'F15.3)'		!generate format string	    

	  DO d=1,dayyear
	    DO ih=1,nt
		  hour=ih
          step=(d-1)*nt+hour
		  DO k=1,5
	        WRITE(11,fmtstr)t,d,hour,k,(lakeinflow_hrr(step,imun,k),imun=1,subasin)
			!WRITE(11,'(4I6,<subasin>f15.3)')t,d,hour,k,(lakeinflow_hrr(step,imun,k),imun=1,subasin)
		  ENDDO
		ENDDO
	  ENDDO
    CLOSE(11)

    OPEN(11,FILE=pfadn(1:pfadi)//'lake_outflow.out',STATUS='old',  &
		 POSITION='append')
	  DO d=1,dayyear
	    DO ih=1,nt
		  hour=ih
          step=(d-1)*nt+hour
		  DO k=1,5
	        WRITE(11,fmtstr)t,d,hour,k,(lakeoutflow_hrr(step,imun,k),imun=1,subasin)
		    !WRITE(11,'(4I6,<subasin>f15.3)')t,d,hour,k,(lakeoutflow_hrr(step,imun,k),imun=1,subasin)
		  ENDDO
		ENDDO
	  ENDDO
    CLOSE(11)

    OPEN(11,FILE=pfadn(1:pfadi)//'lake_volume.out',STATUS='old',  &
		 POSITION='append')
	  DO d=1,dayyear
	    DO ih=1,nt
		  hour=ih
          step=(d-1)*nt+hour
		  DO k=1,5
	        WRITE(11,fmtstr)t,d,hour,k,(lakewater_hrr(step,imun,k),imun=1,subasin)
			!WRITE(11,'(4I6,<subasin>f15.3)')t,d,hour,k,(lakewater_hrr(step,imun,k),imun=1,subasin)
		  ENDDO
		ENDDO
	  ENDDO
    CLOSE(11)

    OPEN(11,FILE=pfadn(1:pfadi)//'lake_retention.out',STATUS='old',  &
		 POSITION='append')
	  DO d=1,dayyear
	    DO ih=1,nt
		  hour=ih
          step=(d-1)*nt+hour
		  DO k=1,5
	        WRITE(11,fmtstr)t,d,hour,k,(lakeretention_hrr(step,imun,k),imun=1,subasin)
			!WRITE(11,'(4I6,<subasin>f15.3)')t,d,hour,k,(lakeretention_hrr(step,imun,k),imun=1,subasin)
		  ENDDO
		ENDDO
	  ENDDO
    CLOSE(11)

    OPEN(11,FILE=pfadn(1:pfadi)//'lake_vollost.out',STATUS='old',  &
		 POSITION='append')
	  DO d=1,dayyear
	    DO ih=1,nt
		  hour=ih
          step=(d-1)*nt+hour
		  DO k=1,5
	        WRITE(11,fmtstr)t,d,hour,k,(maxstorcap_hrr(step,imun,k),imun=1,subasin)
			!WRITE(11,'(4I6,<subasin>f15.3)')t,d,hour,k,(maxstorcap_hrr(step,imun,k),imun=1,subasin)
		  ENDDO
		ENDDO
	  ENDDO
    CLOSE(11)

    IF (dosediment) then
     OPEN(11,FILE=pfadn(1:pfadi)//'lake_sedinflow.out',STATUS='old',  &
		 POSITION='append')
	  DO d=1,dayyear
	    DO ih=1,nt
		  hour=ih
          step=(d-1)*nt+hour
		  DO k=1,5
	        WRITE(11,fmtstr)t,d,hour,k,(lakesedin_hrr(step,imun,k),imun=1,subasin)
			!WRITE(11,'(4I6,<subasin>f15.3)')t,d,hour,k,(lakesedin_hrr(step,imun,k),imun=1,subasin)
		  ENDDO
		ENDDO
	  ENDDO
     CLOSE(11)

     OPEN(11,FILE=pfadn(1:pfadi)//'lake_sedoutflow.out',STATUS='old',  &
		 POSITION='append')
	  DO d=1,dayyear
	    DO ih=1,nt
		  hour=ih
          step=(d-1)*nt+hour
		  DO k=1,5
	        WRITE(11,fmtstr)t,d,hour,k,(lakesedout_hrr(step,imun,k),imun=1,subasin)
			!WRITE(11,'(4I6,<subasin>f15.3)')t,d,hour,k,(lakesedout_hrr(step,imun,k),imun=1,subasin)
		  ENDDO
		ENDDO
	  ENDDO
     CLOSE(11)
    
     OPEN(11,FILE=pfadn(1:pfadi)//'lake_sizedistoutflow.out',STATUS='old',  &
		 POSITION='append')
	  DO d=1,dayyear
	    DO ih=1,nt
		  hour=ih
          step=(d-1)*nt+hour
		  DO g=1,n_sed_class
	        WRITE(11,fmtstr)t,d,hour,g,(lakefracsedout(step,imun,g),imun=1,subasin)
			!WRITE(11,'(4I6,<subasin>f15.6)')t,d,hour,g,(lakefracsedout(step,imun,g),imun=1,subasin)
		  ENDDO
		ENDDO
	  ENDDO
     CLOSE(11)
	ENDIF
    
    
!        do  j=1,12
!        call areasum(monlaketot(j,1:subasin),mesolaketot(j,:),
!     .               statelaketot(j,:),1)
!        enddo
  
END IF
!--------------------------------------------------------------------

RETURN

900   FORMAT(i4)
999   FORMAT(3(1X,f5.0))
998   FORMAT(3(1X,f5.3))
800   FORMAT(3(1X,f7.0))

END SUBROUTINE lake
