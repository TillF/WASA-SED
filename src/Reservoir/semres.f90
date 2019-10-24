SUBROUTINE semres (STATUS,upstream)

! Till: computationally irrelevant: outcommented unused vars
! 2012-09-14

! Till: computationally irrelevant: minor changes to improve compiler compatibility
! 2011-04-29

! Code converted using TO_F90 by Alan Miller
! Date: 2005-08-23  Time: 12:56:59

use common_h
use params_h
use routing_h
use hymo_h
use time_h
use reservoir_h

IMPLICIT NONE

INTEGER, INTENT(IN)                  :: STATUS
INTEGER, INTENT(IN OUT)                  :: upstream



! STATUS of CALL (0=initialization run, 1=initialization year,
!                 2=calculation day,    3=finalization year)

INTEGER :: i,j,j1,id,ih,dummy1,dummy2,nbrsec1,k,b,p,cont,istate,dummy1a,ka !irout,imun,f,q,
INTEGER :: g, npt,nbrbat1

CHARACTER(20) :: subarea,section
INTEGER :: c,dummy14

! Computation of fractional erosion for each cross section (fractional suspended load transport)
!Ge the weighting factor to include incoming sediment into the carrying capacity has to be read
!Ge in do.dat
REAL :: weigth_factor
!real :: sediment_outflow

real :: conc_inflow
REAL :: recov


real :: loadincoming !TAN,
!real :: predom_grain(200),perc_grain(200) !frsedinflow(366*nt,subasin,n_sed_class)
!real :: coarsest_grain
real :: thickness_act(200)
real :: dummy3
real :: elev
real :: x_p1,x_p2,side_p1,side_p2 !,y_p1,y_p2,factor_x
integer :: p1,p2

real :: dummy4,dummy5,dummy6,dummy7
!INTEGER :: dummy8
real :: weight,retent(200) !discharge_mod(200),
!real :: perc1(10),perc2(10),discharg1,discharg2

real :: discharge(200),topwidth(200),area_mean(200)
!real :: dep_fac
real :: factor_over,factor_intake,factor_bottom
real :: factor_over2,factor_intake2,factor_bottom2
CHARACTER (LEN=200) :: cdummy
REAL :: temp1 !,temp2
real :: areasec(200),widthsec(200),areasec2(200),widthsec2(200)
real :: elevhelp,areahelp(200),volhelp(200),areahelp2(200),volhelp2(200)
real :: dummy9 !,dummy12,dummy13
real :: accum1,accum2,dummy15,gsize(2)
character(len=1000) :: fmtstr	!string for formatting file output



! -----------------------------------------------------------------------
IF (STATUS == 0) THEN

  DO i=1,subasin
    DO g=1,n_sed_class
      sed_qlateral(i,g)=0.
	ENDDO
  ENDDO
  DO i=1,subasin
    DO id=1,dayyear*nt
      decvolact(id,i)=0.
      decstorcap(id,i)=0.
      decdamalert(id,i)=0.
      decdamdead(id,i)=0.
      decmaxdamarea(id,i)=0.
    END DO
	cum_sedimentation(i)=0.
	fcav(i)=0
  END DO

!** read particle size classes
  IF (reservoir_check == 1) THEN
	OPEN(11,FILE=pfadp(1:pfadj)// 'part_class.dat', IOSTAT=k,STATUS='old')
	IF (k/=0) THEN					!part_class.dat not found
		write(*,*)'WARNING: '//pfadp(1:pfadj)// 'part_class.dat not found, using defaults (one size class)'
        n_sed_class=1				!treat one particle size class only
		allocate(upper_limit(n_sed_class))
		upper_limit(1)=2.		!set upper limit of sediment that is considered to 2 mm
	else							!part_class.dat sucessfully opened
		READ(11,*)
		READ(11,*)
		i=0							!for counting number of particle size classes
		DO WHILE (.TRUE.)
			READ(11,'(a)', IOSTAT=k) cdummy
			IF (k/=0) THEN			!no further line
				n_sed_class=i		!store total number of particle size classes that have been read
				exit				!exit loop
			END IF
			i=i+1					!increase record counter
		END DO

		allocate(upper_limit(n_sed_class))	!allocate memory for particle size classes to be read
		allocate(particle_classes(n_sed_class))	!allocate memory for particle size classes to be read
		rewind(11)
		READ(11,*)
		READ(11,*)

		DO i=1,n_sed_class
			READ(11,*,  IOSTAT=k) j, temp1
			IF (k/=0) THEN			!no further line
				n_sed_class=i-1		!correct total number of particle size classes that have been read
				exit				!exit loop
			END IF
			upper_limit(i)=temp1		!store upper limit of particle size class
		END DO
		CLOSE(11)

		!convert class limits to "mean" diameter
		DO i=2,n_sed_class
			particle_classes(i)=sqrt(upper_limit(i)*upper_limit(i-1))		!use geometric mean as a mean diameter for each class
		END DO
		particle_classes(1)=upper_limit(1)/sqrt(2.)	!lower limit of the particle size class is assumed to be half of the upper limit (geometric mean is then applied)
!particle_classes(1)=0.020
!write(*,*)(particle_classes(i),i=1,n_sed_class)
	END IF
  ENDIF

! Read geometric parameters of the reservoir
  DO i=1,subasin
	nbrsec(i)=0
  ENDDO
  OPEN(11,FILE=pfadp(1:pfadj)// 'Reservoir/hydraul_param.dat', IOSTAT=istate,STATUS='old')
	IF (istate/=0) THEN					!hydraul_param.dat not found
	  write(*,*)'WARNING: '//pfadp(1:pfadj)// 'hydraul_param.dat not found, using defaults'
    ELSE
	  READ(11,*);READ(11,*)
      READ(11,*,IOSTAT=ka) dummy1,dummy2
      DO i=1,subasin
        IF (dummy1==id_subbas_extern(i)) THEN
          nbrsec(i)=dummy2
	    ENDIF
	  ENDDO
      DO WHILE (ka==0)
	    READ(11,*, IOSTAT=ka) dummy1,dummy2					!read next line in file
        DO i=1,subasin
          IF (dummy1==id_subbas_extern(i)) THEN
		    dummy1a=i
            nbrsec(i)=dummy2
		  ENDIF
	    ENDDO
      END DO
	ENDIF
  CLOSE (11)

  OPEN(11,FILE=pfadp(1:pfadj)// 'Reservoir/hydraul_param.dat', IOSTAT=istate,STATUS='old')
	IF (istate==0) THEN
      READ(11,*);READ(11,*)
      DO i=1,subasin
        nbrsec1=nbrsec(i)
        IF (nbrbat(i) /= 0 .AND. nbrsec(i) /= 0) THEN
          READ(11,*) dummy1,dummy2,(manning_sec(j,i),j=1,nbrsec1)
          READ(11,*) dummy1,dummy2,(dist_sec(j,i),j=1,nbrsec1)
        ELSE IF (nbrbat(i) == 0.AND.nbrsec(i) /= 0) THEN
          WRITE(*,*)'ERROR - if sections are defined, the file cav.dat must be given'
          WRITE(*,*)'subasin:',id_subbas_extern(i)
          STOP
        ELSE IF (nbrsec(i) == 0) THEN
          dummy1=id_subbas_extern(i)
        END IF
        IF (dummy1 /= id_subbas_extern(i)) THEN
          WRITE(*,*) 'ERROR: Sub-basin-IDs in file hydraul_param.dat must have the same ordering scheme as in hymo.dat'
          STOP
        END IF
      END DO
	ENDIF
  CLOSE(11)

!Distance from the cross section to the dam (m)
  DO i=1,subasin
    nbrsec1=nbrsec(i)
    IF (nbrsec(i) /= 0) THEN
      DO j=1,nbrsec1
	    cumlength_sec(j,i)=0.
	    DO j1=j,nbrsec1
		  cumlength_sec(j,i)=cumlength_sec(j,i)+dist_sec(j1,i)
		ENDDO
!write(*,'(2I4,2F15.3)')id_subbas_extern(i),j,dist_sec(j,i),cumlength_sec(j,i)
	  ENDDO
	ENDIF
  ENDDO

!Length of the reach represented by each cross section (m)
  DO i=1,subasin
    nbrsec1=nbrsec(i)
    IF (nbrsec(i) /= 0) THEN
      DO j=1,nbrsec1
        IF (j == 1) THEN
          length_sec(j,i)=dist_sec(j,i)
        ELSE IF (j == nbrsec(i)) THEN
          length_sec(j,i)=dist_sec(j,i)+(dist_sec(j-1,i)/2.)
        ELSE
          length_sec(j,i)=(dist_sec(j,i)+dist_sec(j-1,i))/2.
        END IF
      END DO
    endif
  enddo

! Read sedimentological parameters
  OPEN(11,FILE=pfadp(1:pfadj)// 'Reservoir/sed.dat',STATUS='unknown')
  READ(11,*);READ(11,*)
  DO i=1,subasin
    IF (storcap(i) > 0.) THEN
     IF (nbrsec(i) == 0) read(11,*) dummy1,dry_dens(i)
     IF (nbrsec(i) /= 0) read(11,*) dummy1,dry_dens(i),factor_actlay(i)
!     IF (nbrsec(i) /= 0) read(11,*) dummy1,dry_dens(i),sed_flag(i)
	ENDIF
    IF (storcap(i) == 0.) dummy1=id_subbas_extern(i)
    IF (dummy1 /= id_subbas_extern(i)) THEN
      WRITE(*,*) 'ERROR: Sub-basin-IDs in file sed.dat must have the same ordering scheme as in hymo.dat'
      STOP
    END IF
!write(*,*) dummy1,dry_dens(i),factor_actlay(i)
!write(*,*)param_a(i,1),param_b(i,1)
!write(*,*)param_a(i,2),param_b(i,2)
  END DO
  CLOSE(11)

!  OPEN(11,FILE=pfadp(1:pfadj)// 'Reservoir/sed.dat',STATUS='unknown')
!  READ(11,*)
!  DO i=1,subasin
!    IF (nbrsec(i) /= 0) THEN
!	  IF (sed_flag(i)==0) read(11,*) dummy1,dummy3,dummy2,factor_actlay(i)
!	  IF (sed_flag(i)==1) read(11,*) dummy1,dummy3,dummy2,factor_actlay(i),factor_bottom(i),pt_long0(i)
!	ENDIF
!write(*,*) dummy1,pt_long0(i)
!  END DO
!  CLOSE(11)
!stop

! cross sections geometry / initial bed elevation
  g=0 !flag for indicating file coherence
  DO i=1,subasin
    IF (nbrsec(i) /= 0) THEN
      WRITE(subarea,*)id_subbas_extern(i)
	  OPEN(11,FILE=pfadp(1:pfadj)//'Reservoir/cross_sec_'//trim(adjustl(subarea))//'.dat',STATUS='unknown')
		READ(11,*);READ(11,*)
        nbrsec1=nbrsec(i)
        DO j=1,nbrsec1
          READ(11,*, IOSTAT=istate) dummy1,dummy2,npoints(j,i),  &
            (x_sec0(m,j,i),y_sec0(m,j,i),m=1,npoints(j,i))
            IF (istate/=0) THEN
                write(*,"(A)")"ERROR: Premature end of file in cross_sec_"//trim(adjustl(subarea))//".dat. Check specs in hydraul_param.dat."
                stop
            end if
            do m=1,npoints(j,i)-1 !check for increasing x-coordinates
               if (x_sec0(m,j,i) >= x_sec0(m+1,j,i)) then
                    write(*,'(A,i0,A,i0,A,f6.1,A)')"ERROR: x-coordinates in reservoir cross-section must be increasing (line ", j+2,", point ", m,", x=", x_sec0(m,j,i), ",). "
                    g=1 !indicate error
               end if
            end do

		  id_sec_extern(j,i)=dummy2
!write(*,*)dummy1,dummy2,npoints(j,i),(x_sec0(m,j,i),y_sec0(m,j,i),m=1,npoints(j,i))
        END DO
      CLOSE(11)
      if (g == 1) then
            write(*,*)"Increase floating point precision or decrease resolution."
            stop
       end if
	ENDIF
  ENDDO

  DO i=1,subasin
   IF (sed_flag(i)==1) THEN
    IF (nbrsec(i) /= 0) THEN
	  pt_long0(res_index(i))=0			!pt_long0(i) should be read in the file sed.dat (not implemented)
	  pt_long(res_index(i))=pt_long0(res_index(i)) !Anne changed i to res_index(i)
      nbrsec1=nbrsec(i)
	  dummy3=0.
	  dummy4=0.
	  dummy1=pt_long0(res_index(i))
	  length_plunge(i)=cumlength_sec(dummy1,i)
	  IF (nbrsec1>1 .and. dummy1/=nbrsec1) then
	    dummy3=y_sec0(1,dummy1,i)
	    dummy4=y_sec0(1,dummy1+1,i)
!write(*,*)nbrsec1,dummy1,dummy3,dummy4
        DO m=2,npoints(dummy1,i)
          IF (dummy3 > y_sec0(m,dummy1,i)) dummy3=y_sec0(m,dummy1,i)
          IF (dummy4 > y_sec0(m,dummy1+1,i)) dummy4=y_sec0(m,dummy1+1,i)
		ENDDO
	    slope_long(res_index(i))=min((dummy3-dummy4)/dist_sec(dummy1,i),0.005)  !A
	  ENDIF
	  IF (nbrsec1>1 .and. dummy1==nbrsec1) then
	    dummy3=y_sec0(1,dummy1-1,i)
	    dummy4=y_sec0(1,dummy1,i)
        DO m=2,npoints(dummy1,i)
          IF (dummy3 > y_sec0(m,dummy1-1,i)) dummy3=y_sec0(m,dummy1-1,i)
          IF (dummy4 > y_sec0(m,dummy1,i)) dummy4=y_sec0(m,dummy1,i)
		ENDDO
	    slope_long(res_index(i))=min((dummy3-dummy4)/dist_sec(dummy1-1,i),0.005) !A
	  ENDIF
!write(*,*) dummy1,dummy1,dummy3,dummy4,slope_long(i)
     ENDIF
	ENDIF
  ENDDO
!stop

! initialization of the cross section geometry / actual bed elevation
  DO i=1,subasin
    nbrsec1=nbrsec(i)
    IF (nbrsec(i) /= 0) THEN
      DO j=1,nbrsec1
        DO m=1,npoints(j,i)
          x_sec(m,j,res_index(i))=x_sec0(m,j,i)
          y_sec(m,j,res_index(i))=y_sec0(m,j,i)
        END DO
!write(*,*)j,npoints(j,i),(x_sec0(m,j,i),y_sec0(m,j,i),m=1,npoints(j,i))
      END DO
    END IF
  END DO

!Ge cross sections geometry / original bed elevation
!Ge it represents the sediment layer thickness below the initial bed elevation
!Ge y_original(m,j,i) <= y_sec0(m,j,i)
  DO i=1,subasin
    IF (nbrsec(i) /= 0) THEN
      WRITE(subarea,*)id_subbas_extern(i)
      OPEN(11,FILE=pfadp(1:pfadj)// 'Reservoir/original_sec_'//trim(adjustl(subarea))//'.dat', IOSTAT=istate,STATUS='old')
        IF (istate/=0) THEN					!original_sec_"subbasin".dat not found
          write(*,*)'WARNING: '//pfadp(1:pfadj)// 'original_sec.dat not found, using defaults'
          nbrsec1=nbrsec(i)
          DO j=1,nbrsec1
            DO m=1,npoints(j,i)
              y_original(m,j,res_index(i))=y_sec0(m,j,i)
			ENDDO
		  ENDDO
		  dummy2=0
	    ELSE
	      READ(11,*);READ(11,*)
          DO j=1,nbrsec1
		    READ(11,*)dummy1,dummy2,(y_original(m,j,res_index(i)),m=1,npoints(j,i))
		  ENDDO
		  dummy2=1
	    ENDIF
	  CLOSE(11)
	  IF (dummy2==0) THEN
        DO j=1,nbrsec1
          DO g=1,n_sed_class
            frac_actlay(g,j,res_index(i))=0.
		  ENDDO
		ENDDO
	  ELSE
        OPEN(11,FILE=pfadp(1:pfadj)// 'Reservoir/sizedist_'//trim(adjustl(subarea))//'.dat', STATUS='unknown')
	      READ(11,*);READ(11,*)
          DO j=1,nbrsec1
            READ(11,*)dummy1,dummy2,(frac_actlay(g,j,res_index(i)),g=1,n_sed_class)
		  ENDDO
		CLOSE(11)
	  ENDIF
	ENDIF
  ENDDO
  DO i=1,subasin
    IF (nbrsec(i) /= 0) THEN
	  nbrsec1=nbrsec(i)
      DO j=1,nbrsec1
        DO m=1,npoints(j,i)
!write(*,*)j,m,y_original(m,j,i),y_sec0(m,j,i)
	      IF (y_original(m,j,res_index(i)) > y_sec0(m,j,i)) THEN
		    write(*,*)'ERROR - original bed elevation over the initial bed elevation:'
			write(*,*)'    subbasin=',id_subbas_extern(i),'section=',id_sec_extern(j,i)
		    stop
		  ENDIF
		ENDDO
      ENDDO
	ENDIF
  ENDDO

! Initial deposition area of each cross section (m2)
  DO i=1,subasin
    nbrsec1=nbrsec(i)
    IF (nbrsec(i) /= 0) THEN
      DO j=1,nbrsec1
	    area_sedim(j,i)=0.
        DO m=2,npoints(j,i)
          area_sedim(j,i)=area_sedim(j,i)+  &
			(x_sec(m,j,res_index(i))-x_sec(m-1,j,res_index(i)))*  &
 			((y_sec(m-1,j,res_index(i))-y_original(m-1,j,res_index(i)))+ &
 			(y_sec(m,j,res_index(i))-y_original(m,j,res_index(i))))/2.
!write(*,*)m,y_sec(m,j,i),y_original(m,j,i),area_sedim(j,i)
	    ENDDO
	  ENDDO
	ENDIF
  ENDDO

! Initial deposition volume of each cros section (m3)
  DO i=1,subasin
    nbrsec1=nbrsec(i)
    IF (nbrsec(i) /= 0) THEN
      DO j=1,nbrsec1
        vol_sedim(j,i)=area_sedim(j,i)*length_sec(j,i)
!write(*,'(2I4,5F15.1)')id_subbas_extern(i),j,area_sedim(j,i),length_sec(j,i),vol_sedim(j,i)

      ENDDO
	ENDIF
  ENDDO

! Initial sediment volume of the sub-basins' reservoir (m3)
  DO i=1,subasin
    volbed0(i)=0.
    nbrsec1=nbrsec(i)
    IF (nbrsec(i) /= 0) THEN
      DO j=1,nbrsec1
        volbed0(i)=volbed0(i)+vol_sedim(j,i)
!write(*,'(2I4,5F18.3)')id_subbas_extern(i),j,vol_sedim(j,i),volbed0(i)
      ENDDO
	ENDIF
  ENDDO

  DO i=1,subasin
    nbrsec1=nbrsec(i)
    IF (nbrsec(i) /= 0) THEN
      DO j=1,nbrsec1
        npt=npoints(j,i)
        DO m=1,npt
          y_actlay(m,j,res_index(i))=y_sec0(m,j,i)
        END DO
      END DO
    END IF
  END DO
  DO i=1,subasin
    nbrsec1=nbrsec(i)
    IF (nbrsec(i) /= 0) THEN
      DO j=1,nbrsec1
	    totvol_actlay0(j,i)=0.
        DO g=1,n_sed_class
          frac_comlay(g,j,res_index(i))=frac_actlay(g,j,res_index(i))
          frac_toplay(g,j,res_index(i))=0.
          frac_susp(g,j,res_index(i))=0.
	      frvol_actlay(g,j,i)=frac_actlay(g,j,res_index(i))*vol_sedim(j,i)
	      frvol_actlay0(g,j,i)=frvol_actlay(g,j,i)
          totvol_actlay0(j,i)=totvol_actlay0(j,i)+frvol_actlay0(g,j,i)
!write(*,'(2I4,5F15.4)')id_subbas_extern(i),j,vol_sedim(j,i),frac_actlay(g,j,i),frvol_actlay(g,j,i),totvol_actlay(j,i)
!write(*,'(2I4,5F15.4)')id_subbas_extern(i),j,totvol_actlay0(j,i),frvol_actlay0(g,j,i)
        END DO
!write(*,'(2I4,5F15.4)')id_subbas_extern(i),j,totvol_actlay(j,i),dummy4,volbed0(i)
      END DO
    END IF
  END DO

!Identification of main channel in the cross sections
  OPEN(11,FILE=pfadp(1:pfadj)// 'Reservoir/main_channel.dat', IOSTAT=istate,STATUS='old')
	IF (istate/=0) THEN					!main_channel.dat not found
      write(*,*)'WARNING: '//pfadp(1:pfadj)// 'main_channel.dat not found, using defaults'
      DO i=1,subasin
	    sed_flag(i)=0 !0 = changes on sideslope is not controlled
	  ENDDO
    ELSE
      READ(11,*);READ(11,*)
      DO i=1,subasin
 	    sed_flag(i)=1 !changes on sideslope is controlled avoiding steeper slopes by erosion processes
        nbrsec1=nbrsec(i)
        IF (nbrbat(i) /= 0 .AND. nbrsec(i) /= 0) THEN
          READ(11,*) dummy1,dummy2,(pt1(j,res_index(i)),j=1,nbrsec1)
          READ(11,*) dummy1,dummy2,(pt2(j,res_index(i)),j=1,nbrsec1)    !A
        ELSE IF (nbrbat(i) == 0.AND.nbrsec(i) /= 0) THEN
          WRITE(*,*)'ERROR - if sections are defined, the file cav.dat must be given'
          WRITE(*,*)'subasin:',id_subbas_extern(i)
          STOP
        ELSE IF (nbrsec(i) == 0) THEN
          dummy1=id_subbas_extern(i)
        END IF
        IF (dummy1 /= id_subbas_extern(i)) THEN
          WRITE(*,*) 'ERROR: Sub-basin-IDs in file main_channel.dat must have the same ordering scheme as in hymo.dat'
          STOP
        END IF
      END DO
	ENDIF
  CLOSE(11)

  DO i=1,subasin
    nbrsec1=nbrsec(i)
    IF (sed_flag(i)==1) THEN
     IF (nbrsec(i) /= 0) THEN
      DO j=1,nbrsec1
	    npt=npoints(j,i)
		p1=pt1(j,res_index(i))
		p2=pt2(j,res_index(i))
	    if (p1/=-999) then
		  sideslope_pt1(j,res_index(i))=abs(y_sec(p1,j,res_index(i))-y_sec(p1+1,j,res_index(i)))/abs(x_sec(p1,j,res_index(i))-x_sec(p1+1,j,res_index(i)))   !A changed i to res_index(i)
		  sideslope_pt1(j,res_index(i))=max(sideslope_pt1(j,res_index(i)),.05)                                  !A
!		  elev=y_sec(p1,j,i)+2.
		  elev=y_sec(p1,j,res_index(i))+(min(y_sec(1,j,res_index(i)),y_sec(npt,j,res_index(i)))-y_sec(p1,j,res_index(i)))/2.
          DO m=2,npoints(j,i)-1
            if (y_sec(m,j,res_index(i)) < elev) then
		      pt3(j,res_index(i))=m-1
		      exit
		    endif
	      ENDDO
		else
		  pt3(j,res_index(i))=-999
		endif
	    if (p2/=-999) then
	      sideslope_pt2(j,res_index(i))=abs(y_sec(p2,j,res_index(i))-y_sec(p2-1,j,res_index(i)))/abs(x_sec(p2,j,res_index(i))-x_sec(p2-1,j,res_index(i)))   !A
		  sideslope_pt2(j,res_index(i))=max(sideslope_pt1(j,res_index(i)),.05)                                  !A
!		  elev=y_sec(p2,j,i)+2.
		  elev=y_sec(p2,j,res_index(i))+(min(y_sec(1,j,res_index(i)),y_sec(npt,j,res_index(i)))-y_sec(p2,j,res_index(i)))/2.
          DO m=1,npoints(j,i)
            if (y_sec(m,j,res_index(i)) < elev) then
		      pt4(j,res_index(i))=m+1
		    endif
	      ENDDO
		else
		  pt4(j,res_index(i))=-999
		endif
!write(*,'(6I5,8F15.4)')id_subbas_extern(i),j,pt1(j,i),pt2(j,i),pt3(j,i),pt4(j,i),sideslope_pt1(j,i),sideslope_pt2(j,i),elev
      END DO
     END IF
    END IF
  END DO
!stop

! stage-area and stage-vulume curves given in the file cav.dat is disregarded. Values derived from cross section are used instead
  DO i=1,subasin
    IF (nbrsec(i) /= 0) THEN
	 nbrbat1=nbrbat(i)
     DO b=1,nbrbat1
      elevhelp=elev_bat(b,i)
      DO j=1,nbrsec(i)
        areasec2(j)=0.
        widthsec2(j)=0.

        npt=npoints(j,i)
        DO m=2,npt
		  IF (elevhelp >= y_sec(m-1,j,res_index(i))  &
			.AND.elevhelp >= y_sec(m,j,res_index(i))) THEN
			areasec2(j)=areasec2(j)+(x_sec(m,j,res_index(i))-x_sec(m-1,j,res_index(i)))  &
				*(elevhelp-(y_sec(m,j,res_index(i))+y_sec(m-1,j,res_index(i)))/2.)
			widthsec2(j)=widthsec2(j)+  &
				(x_sec(m,j,res_index(i))-x_sec(m-1,j,res_index(i)))
		  ELSE IF (elevhelp < y_sec(m-1,j,res_index(i))  &
			.AND.elevhelp >= y_sec(m,j,res_index(i))) THEN
			areasec2(j)=areasec2(j)+(((elevhelp-y_sec(m,j,res_index(i)))**2.)/ &
				(2.*ABS(y_sec(m,j,res_index(i))-y_sec(m-1,j,res_index(i)))/  &
				(x_sec(m,j,res_index(i))-x_sec(m-1,j,res_index(i)))))
			widthsec2(j)=widthsec2(j)+  &
				((x_sec(m,j,res_index(i))-x_sec(m-1,j,res_index(i)))*  &
				(elevhelp-y_sec(m,j,res_index(i)))/  &
				(y_sec(m-1,j,res_index(i))-y_sec(m,j,res_index(i))))
		  ELSE IF (elevhelp >= y_sec(m-1,j,res_index(i))  &
			.AND.elevhelp < y_sec(m,j,res_index(i))) THEN
			areasec2(j)=areasec2(j)+(((elevhelp-y_sec(m-1,j,res_index(i)))**2.)/ &
				(2.*ABS(y_sec(m,j,res_index(i))-y_sec(m-1,j,res_index(i)))/  &
				(x_sec(m,j,res_index(i))-x_sec(m-1,j,res_index(i)))))
			widthsec2(j)=widthsec2(j)+  &
				((x_sec(m,j,res_index(i))-x_sec(m-1,j,res_index(i)))*  &
				(elevhelp-y_sec(m-1,j,res_index(i)))/  &
				(y_sec(m,j,res_index(i))-y_sec(m-1,j,res_index(i))))
		  END IF
	    END DO
	   END DO

       volhelp2(b)=0.
       areahelp2(b)=0.

       DO j=1,nbrsec(i)
        IF (j /= nbrsec(i)) THEN
	      volhelp2(b)=volhelp2(b)+(areasec2(j)+areasec2(j+1))*dist_sec(j,i)/2.
	      areahelp2(b)=areahelp2(b)+(widthsec2(j)+widthsec2(j+1))*dist_sec(j,i)/2.
        ELSE
	      volhelp2(b)=volhelp2(b)+areasec2(j)*dist_sec(j,i)/2.
	      areahelp2(b)=areahelp2(b)+widthsec2(j)*dist_sec(j,i)/2.
        END IF
!write(*,*)j,areasec2(j),dist_sec(j,i),areahelp2(b),volhelp2(b)
	   END DO
!write(*,*)j,areahelp2(b),volhelp2(b)
     ENDDO
     DO b=1,nbrbat1
	   IF (elev_bat(b,i) >= maxlevel(i)) THEN
!write(*,*)b,elev_bat(b,i),maxlevel(i)
         vol_bat(b,i)=vol_bat0(b,i)
         area_bat(b,i)=area_bat0(b,i)
	   ELSE
!write(*,*)b,areahelp2(b),volhelp2(b)
	     vol_bat(b,i)=volhelp2(b)
	     area_bat(b,i)=areahelp2(b)
	   ENDIF
!write(*,*)b,area_bat(b,i),vol_bat(b,i),elev_bat(b,i),maxlevel(i)
	 ENDDO
	ENDIF
  ENDDO


!Ge initialization of output files
!Ge check if could be read in do.dat to print either general results or detailed results
  DO i=1,subasin
    IF (nbrsec(i) /= 0) THEN
      WRITE(subarea,*)id_subbas_extern(i)
	  OPEN(11,FILE=pfadn(1:pfadi)//'res_'//trim(adjustl(subarea))//'_hydraul.out',STATUS='replace')
      IF (f_res_hydraul) then
	    WRITE(11,*)'Subasin-ID, year, day, hour, section-ID, depth_sec(m), watelev_sec(m), area_sec(m**2), topwidth_sec(m), energslope_sec(-), hydrad_sec(m), meanvel_sec(m/s), discharge_sec(m**3/s)'
        CLOSE(11)
	  ELSE
        CLOSE(11, status='delete') !delete any existing file, if no output is desired
	  ENDIF
	ENDIF
  ENDDO

  DO i=1,subasin
    IF (nbrsec(i) /= 0) THEN
      WRITE(subarea,*)id_subbas_extern(i)

      DO j=1,nbrsec(i)
        WRITE(section,*)id_sec_extern(j,i)
	    OPEN(11,FILE=pfadn(1:pfadi)//'res_'//trim(adjustl(subarea))//'_sec'//trim(adjustl(section))// &
			'_bedchange.out',STATUS='replace')
        IF (f_res_hydraul) then
	      WRITE(11,*)'Subasin-ID, section-ID, year, day, hour, nbr. points, y-axis(m)'
          CLOSE(11)
	    ELSE
          CLOSE(11, status='delete') !delete any existing file, if no output is desired
	    ENDIF
	  ENDDO
	ENDIF
  ENDDO

  DO i=1,subasin
   IF (storcap(i) /= 0.) THEN
    WRITE(subarea,*)id_subbas_extern(i)
	OPEN(11,FILE=pfadn(1:pfadi)//'res_'//trim(adjustl(subarea))//'_sedbal.out',STATUS='replace')
    IF (f_res_sedbal) then
	  WRITE(11,*)'Subasin-ID, year, day, hour, sed_input(ton/timestep), sed_output(ton/timestep), sedimentation(ton/timestep), cum_sedimentation(ton)'
      CLOSE(11)
    ELSE
      CLOSE(11, status='delete') !delete any existing file, if no output is desired
    ENDIF
   ENDIF
  ENDDO

  DO i=1,subasin
   IF (storcap(i) /= 0.) THEN
    IF (nbrsec(i) /= 0) THEN
      WRITE(subarea,*)id_subbas_extern(i)
	  OPEN(11,FILE=pfadn(1:pfadi)//'res_'//trim(adjustl(subarea))//'_longitudunal.out',STATUS='replace')
      IF (f_res_longitudunal) then
	    WRITE(11,*)'Subasin-ID, year, day, hour, nbr. sections, minelev_sec(m)'
	  ELSE
        CLOSE(11, status='delete') !delete any existing file, if no output is desired
	  ENDIF
	ENDIF
   ENDIF
  ENDDO

  DO i=1,subasin
   IF (storcap(i) /= 0.) THEN
    WRITE(subarea,*)id_subbas_extern(i)
	OPEN(11,FILE=pfadn(1:pfadi)//'res_'//trim(adjustl(subarea))//'_sedcomposition.out',STATUS='replace')
    IF (f_res_sedcomposition) then
	  WRITE(11,*)'Subasin-ID, year, day, hour, nbr. classes, sedcomp_outflow(-)'
      CLOSE(11)
    ELSE
      CLOSE(11, status='delete') !delete any existing file, if no output is desired
    ENDIF
   ENDIF
  ENDDO

!Ge output file for check some parameters (temporary)
!***************************************************
!OPEN(11,FILE=pfadn(1:pfadi)//'check.out', &
!		STATUS='replace')
!write(11,*)'check paramaters'
!CLOSE(11)

!OPEN(11,FILE=pfadn(1:pfadi)//'check1.out', &
!		STATUS='replace')
!write(11,*)'check paramaters'
!CLOSE(11)

!OPEN(11,FILE=pfadn(1:pfadi)//'check2.out', &
!		STATUS='replace')
!write(11,*)'check paramaters'
!CLOSE(11)

!OPEN(11,FILE=pfadn(1:pfadi)//'check3.out', &
!		STATUS='replace')
!write(11,*)'check paramaters'
!CLOSE(11)
!***************************************************


! end of status 0
END IF


! -----------------------------------------------------------------------
IF (STATUS == 1) THEN
! initialize ...

! initialization of parameters
  DO i=1,subasin
    DO id=1,dayyear*nt
!      sed_susp(id,i)=0.
      sed_ret(id,i)=0.
      sed_overflow(id,i)=0.
      sed_intake(id,i)=0.
      sed_inflow(id,i)=0.
    END DO
  END DO

!Ge read daily data on reservoir level and outflow discharges
  IF (reservoir_check == 1) THEN
    DO i=1,subasin
      IF (t >= damyear(i)) THEN !Andreas
	    IF (storcap(i) > 0.) THEN
          WRITE(subarea,*)id_subbas_extern(i)
	      OPEN(11,FILE=pfadp(1:pfadj)//'Reservoir/frsedinflow_'//trim(adjustl(subarea))//'.dat', &
			STATUS='unknown')
          read(11,*)
          cont=dtot*nt
          DO id=1,cont
            READ(11,*)
          END DO
          DO id=1,dayyear*nt
            READ(11,*) dummy1,dummy2,(frsedinflow(id,i,g),g=1,n_sed_class)
!write(*,'(2I4,<n_sed_class>F10.4)')dummy1,dummy2,(frsedinflow(id,i,g),g=1,n_sed_class)
!write(*,*)dummy1,dummy2,(frsedinflow(id,i,g),g=1,n_sed_class)
!if (id==70 .and. id_subbas_extern(i)==21)stop
	      ENDDO
	      CLOSE(11)
        ENDIF
      ENDIF
	ENDDO
  ENDIF

! end of status 1
END IF


! -----------------------------------------------------------------------
IF (STATUS == 2) THEN

! sediment inflow into the reservoir is the generate sediment flow
! from the rainfall-runoff processes, calculated using the WASA model

  IF (reservoir_check == 1) THEN
    DO g=1,n_sed_class
      sediment_in(upstream,g)=frsedinflow(step,upstream,g)
!write(*,*)step,upstream,frsedinflow(step,upstream,g),sediment_in(upstream,g)
	ENDDO
  ENDIF
!write(*,'(2I4,<n_sed_class>F10.4)')step,upstream,(sediment_in(upstream,g),g=1,n_sed_class)
!if (id==70)stop


! Determination of total sediment inflow (ton/timestep) and inflow sediment concentration (g/l)
  do g=1,n_sed_class
	sed_inflow(step,upstream)=sed_inflow(step,upstream)+sediment_in(upstream,g)+sed_qlateral(upstream,g)
!write(*,'(2I4,3F15.3)')id_subbas_extern(upstream),g,sediment_in(upstream,g),sed_inflow(step,upstream),sed_qlateral(upstream,g)
!write(*,*)upstream,g,sediment_in(upstream,g),sed_inflow(step,upstream),sed_qlateral(upstream,g)
  enddo

!****************************************************************
!sed_inflow(step,upstream)=0.
!sediment_in(upstream,1)=576000.
!  do g=1,n_sed_class
!sediment_in(upstream,1)=576000.
!sediment_in(upstream,g)=2*(1./particle_classes(g))**2.
!sed_inflow(step,upstream)=sed_inflow(step,upstream)+sediment_in(upstream,g)
!  enddo
!write(*,'(F13.1,<n_sed_class>F12.1)')sed_inflow(step,upstream),(sediment_in(upstream,g),g=1,n_sed_class)
!stop
!****************************************************************

if (qinflow(step,upstream)==0.) then
    conc_inflow = 0. 
else
    conc_inflow=sed_inflow(step,upstream)/(qinflow(step,upstream)*(86400./nt)*.001)
end if


!write(*,*)step,upstream,(sediment_in(upstream,g),g=1,n_sed_class)

  if (sed_inflow(step,upstream) > 0.) then
	do g=1,n_sed_class
	  frsediment_in(upstream,g)=(sediment_in(upstream,g)+sed_qlateral(upstream,g))/sed_inflow(step,upstream)
	enddo
  else
	do g=1,n_sed_class
	  frsediment_in(upstream,g)=0.
	enddo
  endif

!****************************************************************
!    diam(1)=0.0014/1000.  !medium clay			(0.0010 - 0.0020)
!    diam(2)=0.0028/1000.  !coarse clay			(0.0020 - 0.0039)
!    diam(3)=0.0055/1000.  !very fine silt		(0.0039 - 0.0078)
!    diam(4)=0.0110/1000.  !fine silt			(0.0078 - 0.0156)
!    diam(5)=0.0221/1000.  !medium silt			(0.0156 - 0.0313)
!    diam(6)=0.0442/1000.  !coarse silt			(0.0313 - 0.0625)
!    diam(7)=0.088/1000.  !very fine sand		(0.0625 - 0.1250)
!    diam(8)=0.350/1000.  !medium sand			(0.1250 - 0.5000)
!    diam(9)=1.400/1000.  !very coarse sand		(0.5000 - 2.0000)
!    diam(10)=2.800/1000.  !very fine gravel	(2.0000 - 4.0000)

! values from field measurements
!perc1(1)=0.2998
!perc1(2)=0.2377
!perc1(3)=0.2221
!perc1(4)=0.1726
!perc1(5)=0.0641
!perc1(6)=0.0035
!perc1(7)=0.0001
!perc1(8)=0.0001
!perc1(9)=0.0000
!perc1(10)=0.0000

!perc2(1)=0.116
!perc2(2)=0.148
!perc2(3)=0.153
!perc2(4)=0.153
!perc2(5)=0.154
!perc2(6)=0.112
!perc2(7)=0.027
!perc2(8)=0.055
!perc2(9)=0.055
!perc2(10)=0.027
!****************************************************************

! sediment diameter (m)
  do g=1,n_sed_class
	diam(g)=particle_classes(g)/1000.
  enddo


! 1) calculation of the sediment balance in the reservoir without taking into account
!    the longitudinal sediment transport processes (using pond performance modelling)

  IF (nbrsec(upstream) == 0) THEN
!write(*,'(2I4,10F15.6)') step,id_subbas_extern(upstream),sed_inflow(step,upstream)
    IF (sed_inflow(step,upstream) /= 0.) then
      call sedbal(upstream)
	ELSE
	  sed_outflow(step,upstream)=0.
	  sedimentation(step,upstream)=0.
      if (t==tstart .and. step==1) cum_sedimentation(upstream)=0.
      IF (t>tstart .and. t== damyear(upstream) .and. step==1) cum_sedimentation(upstream)=0.

      cum_sedimentation(upstream)=cum_sedimentation(upstream)+sedimentation(step,upstream)
	  decstorcap(step,upstream)=0.
	  do g=1,n_sed_class
		frsediment_out(upstream,g)=0.
	  enddo
	ENDIF


! 2) calculation of the sediment balance in the reservoir taking into account
!    the longitudinal sediment transport processes
  ELSE
! Determination of grain size distribution at the top layer of the most upstream
! cross section provided by size distribution of sediment input
	do g=1,n_sed_class
	  if (sed_inflow(step,upstream) > 0.) then
	    frac_toplay(g,1,res_index(upstream))=frsediment_in(upstream,g)
	  else
	    frac_toplay(g,1,res_index(upstream))=0.
	  endif
	enddo

! 2a) determination of hydraulic parameters
    CALL hydraul_res(upstream)


! 2b) sediment carrying capacity

!Ge Calculation of the sediment diameters D50 and D90 (m)
!Ge mean diameter is calculated using the grain size distribution
    DO j=1,nbrsec(upstream)
	 accum2=0.
     DO g=1,n_sed_class
	   accum2=accum2+frac_actlay(g,j,res_index(upstream))
	 ENDDO
!write(*,'(2I4,<n_sed_class>F7.3)')upstream,j,(frac_actlay(g,j,upstream),g=1,n_sed_class)

     IF (accum2==0.) THEN
	  d50_actlay(j,upstream)=0.
	  d90_actlay(j,upstream)=0.
	 ELSE
      gsize(1)=.5
	  gsize(2)=.9
	  DO c=1,2
        accum1=0.
        accum2=0.
        dummy14=0
        dummy15=0.
        DO g=1,n_sed_class
          if (accum2 <= gsize(c)) then
            dummy14=dummy14+1
            accum1=accum2
            accum2=accum2+frac_actlay(g,j,res_index(upstream))
!write(*,'(4I4,4F14.6)')upstream,j,g,dummy14,accum1,accum2
!read(*,*)
		  else
		    exit
		  endif
		ENDDO
        IF (dummy14 > 1) then
          dummy15=10.**(log10(upper_limit(dummy14-1))+((log10(upper_limit(dummy14))-log10(upper_limit(dummy14-1)))*(gsize(c)-accum1)/(accum2-accum1)))
          dummy15=dummy15/1000.
	    ELSE
	      IF (c == 1) dummy15=diam(1)
	      IF (c == 2) then
		    dummy15=10.**(log10(upper_limit(dummy14)/2.)+((log10(upper_limit(dummy14))-log10(upper_limit(dummy14)/2.))*(gsize(c)-accum1)/(accum2-accum1)))
		    dummy15=dummy15/1000.
		  ENDIF
		ENDIF
	    IF (c == 1) d50_actlay(j,upstream)=dummy15
	    IF (c == 2) d90_actlay(j,upstream)=dummy15
!IF (c == 1)write(*,'(3I4,4F14.6)')upstream,j,dummy14,accum1,accum2,d50_actlay(j,upstream)*1000.
!IF (c == 2)write(*,'(3I4,4F14.6)')upstream,j,dummy14,accum1,accum2,d90_actlay(j,upstream)*1000.
	  ENDDO
	 ENDIF
!write(*,'(3I4,6F14.6)')upstream,j,g,d50_actlay(j,upstream)*1000.,d90_actlay(j,upstream)*1000.
	ENDDO
!if (step==2)stop

! Call of the sub-routine to calculate the sediment carrying capacity, as  defined in the do.dat
    if (reservoir_transport == 1) call eq1_wu(upstream)
    if (reservoir_transport == 2) call eq2_ashida(upstream)
    if (reservoir_transport == 3) call eq3_tsinghua(upstream)
    if (reservoir_transport == 4) call eq4_ackers(upstream)

!DO j=1,nbrsec(upstream)
!write(*,'(I6,10F14.1)')j,(fr_capacity(g,j),g=1,n_sed_class)
!ENDDO
!if (step==1)stop

! Identification of flushing operation
!	p=0
!    DO j=1,nbrsec(upstream)
!      if (watelev_sec(j,upstream) == damelev_mean(step,upstream)) then
!	    p=j
!		exit
!	  endif
!	ENDDO

    IF (sed_routing_flag(upstream)==1) THEN
	 if (p == 0) then
!	if (qbottom(step,upstream) /= 0.) then
      dummy4=0.
	  DO j=pt_long(res_index(upstream)),pt_long0(res_index(upstream))
	    dummy4=dummy4+dist_sec(j,upstream)
	  ENDDO
!write(*,*)dummy4

      dummy3=(minelev_sec(pt_long(res_index(upstream)),upstream)-minelev_sec(pt_long0(res_index(upstream))+1,upstream))/dummy4
!write(*,*)cumlength_sec(pt_long(upstream),upstream),cumlength_sec(pt_long0(upstream)+1,upstream),dummy4,dummy3 !(cumlength_sec(pt_long(upstream),upstream)-cumlength_sec(pt_long0(upstream)+1,upstream)),dist_sec(pt_long0(upstream),upstream)
      IF (dummy3 > slope_long(res_index(upstream))) THEN   !A
	    pt_long(res_index(upstream))=max(pt_long(res_index(upstream))-1,42)		!the section 42 of the Barasona reservoir was assumed to be the upstream limit (variable pt_long_min(upstream) should be read in the file sed.dat
        length_plunge(upstream)=cumlength_sec(max(pt_long(res_index(upstream)),1),upstream)
	  ENDIF
!write(*,'(2I6,F10.2,2F10.6,4F10.2)')pt_long0(upstream),pt_long(upstream),length_plunge(upstream),dummy3,slope_long(upstream),minelev_sec(pt_long0(upstream),upstream),minelev_sec(pt_long0(upstream)+1,upstream),cumlength_sec(pt_long(upstream)&
!,upstream)-cumlength_sec(pt_long0(upstream)+1,upstream)
	 ENDIF
	ENDIF
!    length_plunge(upstream)=.4*cumlength_sec(1,upstream)
!write(*,*)p

! Calculation of sediment volume in the active layer
    DO j=1,nbrsec(upstream)
      area_actlay(j,upstream)=0.
      area_toplay(j,upstream)=0.

      DO m=1,npoints(j,upstream)
	    partarea_actlay(m,j,res_index(upstream))=0.     !A changed upstream to res_index(upstream)
	    partarea_toplay(m,j,res_index(upstream))=0.
	  ENDDO

!	  elev=maxelev_sec(j,upstream)
	  elev=watelev_sec(j,upstream)
!	  elev=min(minelev_sec(j,upstream)+(depth_sec(j,upstream)*1.5),maxelev_sec(j,upstream))
!	  elev=watelev_sec(j,upstream)+((maxelev_sec(j,upstream)-watelev_sec(j,upstream))/20.)
!if (id_sec_extern(j,upstream)==43) elev=max(elev,y_sec(11,j,upstream))
!if (id_sec_extern(j,upstream)==47) elev=max(elev,y_sec(18,j,upstream))
!if (id_sec_extern(j,upstream)==47)elev=max(elev,435.)
!	  elev=min(watelev_sec(j,upstream)+1.5,watelev_sec(j,upstream)+((maxelev_sec(j,upstream)-watelev_sec(j,upstream))/5.))
!	  elev=watelev_sec(j,upstream)
	  if (sed_flag(upstream) == 0) then
        p1=-999
        p2=-999
	  else if (sed_flag(upstream) == 1) then
        p1=pt1(j,res_index(upstream))
        p2=pt2(j,res_index(upstream))
	  endif

	  npt=npoints(j,upstream)

!write(*,'(3I6,6F10.3)')id_sec_extern(j,upstream),p1,p2,elev,watelev_sec(j,upstream)

	  if (p1/=-999 .and. y_sec(p1,j,res_index(upstream))>watelev_sec(j,upstream)) then !main channel was defined for each cross section
        x_p1=x_sec(p1,j,res_index(upstream))
	    side_p1=abs(y_sec(p1,j,res_index(upstream))-y_sec(p1+1,j,res_index(upstream)))/abs(x_sec(p1,j,res_index(upstream))-x_sec(p1+1,j,res_index(upstream)))
	    if (side_p1>sideslope_pt1(j,res_index(upstream))) then      !A changed upstream to res_index(upstream)
		  p1=max(pt3(j,res_index(upstream)),p1-1)
		  x_p1=x_sec(p1,j,res_index(upstream))
	      pt1(j,res_index(upstream))=p1
		endif
		elev=max(y_sec(p1,j,res_index(upstream)),elev)
      else
        DO m=2,npoints(j,upstream)-1
          if (y_sec(m,j,res_index(upstream)) < elev) then
		    p1=m-1
		    exit
		  endif
	    ENDDO
	  endif
	  if (p2/=-999 .and. y_sec(p2,j,res_index(upstream))>watelev_sec(j,upstream)) then !main channel was defined for each cross section
        x_p2=x_sec(p2,j,res_index(upstream))
	    side_p2=abs(y_sec(p2,j,res_index(upstream))-y_sec(p2-1,j,res_index(upstream)))/abs(x_sec(p2,j,res_index(upstream))-x_sec(p2-1,j,res_index(upstream)))
	    if (side_p2>sideslope_pt2(j,res_index(upstream))) then              !A changed upstream to res_index(upstream)
		  p2=min(pt4(j,res_index(upstream)),p2+1)
		  x_p2=x_sec(p2,j,res_index(upstream))
	      pt2(j,res_index(upstream))=p2
		endif
		elev=max(y_sec(p2,j,res_index(upstream)),elev)
	  else
        DO m=1,npoints(j,upstream)
          if (y_sec(m,j,res_index(upstream)) < elev) then
		    p2=m+1
		  endif
	    ENDDO
	  endif

      erosion_level(j,res_index(upstream))=elev

!write(*,'(3I6,6F10.3)')id_sec_extern(j,upstream),p1,p2,elev,watelev_sec(j,upstream)
!if(step==2)stop

!write(*,'(I6,6F10.3)')id_sec_extern(j,upstream),elev,y_sec(p1-1,j,upstream),y_sec(p2+1,j,upstream),y_p1,y_p2
!write(*,'(I6,6F10.3)')id_sec_extern(j,upstream),topwidth_sec(j,upstream),x_p2-x_p1
!write(*,'(I6,6F10.3)')id_sec_extern(j,upstream),x_p1,elev,x_sec(p1-1,j,upstream),x_sec(p1,j,upstream),y_sec(p1-1,j,upstream),y_sec(p1,j,upstream)
!write(*,'(I6,6F10.3)')id_sec_extern(j,upstream),x_p2,elev,x_sec(p2,j,upstream),x_sec(p2+1,j,upstream),y_sec(p2,j,upstream),y_sec(p2+1,j,upstream)
!write(*,'(3I6,6F10.3)')id_sec_extern(j,upstream),p1,p2,elev,x_sec(p1,j,upstream),x_sec(p2,j,upstream),y_sec(p1,j,upstream),y_sec(p2,j,upstream)
!if (step==2)stop

!write(*,*)id_sec_extern(j,upstream),elev,watelev_sec(j,upstream)

      DO m=1,npoints(j,upstream)
        IF (elev > y_sec(m,j,res_index(upstream))) THEN
          geom(m,j)=1
		ELSE
		  geom(m,j)=0
		ENDIF
	  ENDDO
      DO m=2,npoints(j,upstream)-1
		IF (geom(m-1,j) == 1 .AND. geom(m,j) == 1 .AND. geom(m+1,j) == 1 &
			.OR. geom(m-1,j) == 2  .AND. geom(m,j) == 1 .AND. geom(m+1,j) == 1) THEN
		  geom(m,j)=2
		ELSE IF (geom(m-1,j) == 1  .AND. geom(m,j) == 1 .AND. geom(m+1,j) == 0 &
			.OR. geom(m-1,j) == 2  .AND. geom(m,j) == 1 .AND. geom(m+1,j) == 0) THEN
		  geom(m,j)=3
		ENDIF
        IF (elev < y_sec(m-1,j,res_index(upstream))  &
              .AND. elev > y_sec(m,j,res_index(upstream))  &
			  .AND. elev < y_sec(m+1,j,res_index(upstream))) THEN
		  geom(m,j)=4
		ENDIF
	  ENDDO

! Maximum sediment diameter of each cross section (m)
!      coarsest_grain=0.
!      DO g=1,n_sed_class
!        if (frac_actlay(g,j,upstream) /= 0.) then
!		  coarsest_grain=diam(g)
!		endif
!	  enddo

! Determination of active layer thickness using a constant of proportionality (m)	!hourly active layer thickness = 1 * factor_actlay(upstream)*d90_actlay(j,upstream)
!factor_actlay(upstream)=5.
!      thickness_act(j)=factor_actlay(upstream)*d90_actlay(j,upstream)*(24./nt)		!nt is the number of simulation steps per day
!	  if (p /= 0) then
	  if (sed_routing_flag(upstream) == 0) then
        thickness_act(j)=.03*factor_actlay(upstream)*(24./nt)	!0,03 value of active layer thickness along the reservoir derived from the simulation without sediment management technique for the Barasona reservoir
	  else if (sed_routing_flag(upstream) == 1) then
	    if (length_plunge(upstream) > cumlength_sec(j,upstream)) then !0,25 value of active layer thickness close to the dam derived from the simulation with sediment management technique for the Barasona reservoir
          thickness_act(j)=max(.03*factor_actlay(upstream),(.25**factor_actlay(upstream))*(length_plunge(upstream)-cumlength_sec(j,upstream))/length_plunge(upstream))
		else
          thickness_act(j)=.03*factor_actlay(upstream)
		endif
	  endif

!write(*,'(2I4,4F14.3)')upstream,j,length_plunge(upstream),cumlength_sec(j,upstream),thickness_act(j),max(.05,.30*(length_plunge(upstream)-cumlength_sec(j,upstream))/length_plunge(upstream))
!write(*,'(2I4,4F14.6)')upstream,j,factor_actlay(upstream),d90_actlay(j,upstream),thickness_act(j),(24./nt),nt

!Ge erosion height using a normalized distance along the wetted perimeter
      DO m=1,npoints(j,upstream)
	    if(x_sec(m,j,res_index(upstream))>x_sec(p1,j,res_index(upstream)) .and. x_sec(m,j,res_index(upstream))<=x_minelev(j,upstream)) then
		  dummy5=(x_sec(m,j,res_index(upstream))-x_sec(p1,j,res_index(upstream)))/(x_minelev(j,upstream)-x_sec(p1,j,res_index(upstream)))
		else if(x_sec(m,j,res_index(upstream))<x_sec(p2,j,res_index(upstream)) .and. x_sec(m,j,res_index(upstream))>x_minelev(j,upstream)) then
		  dummy5=(x_sec(p2,j,res_index(upstream))-x_sec(m,j,res_index(upstream)))/(x_sec(p2,j,res_index(upstream))-x_minelev(j,upstream))
		else
		  dummy5=.0
		endif
		dummy6=thickness_act(j)*(1.-((1.-dummy5)**2.9))
!dummy6=dummy(j)
        y_actlay(m,j,res_index(upstream))=max(y_sec(m,j,res_index(upstream))-dummy6,y_original(m,j,res_index(upstream)))
!if(j<20)write(*,'(4I4,5F12.6)')t,d,j,m,y_actlay(m,j,upstream),y_sec(m,j,upstream),dummy6,(1.-((1.-(dummy5))**2.9)),dummy5
!if(d==2 .and. j==10)stop
      ENDDO

! Erosion height using a normalized distance along the wetted perimeter
!      DO m=1,npoints(j,upstream)
!        if (y_sec(m,j,upstream)<elev) then
!	      dummy5=(elev-y_sec(m,j,upstream))/(elev-minelev_sec(j,upstream))
!	    else
!          dummy5=0.
!		endif
!		dummy6=thickness_act(j)*(1.-((1.-(dummy5))**2.9))
!        dummy6=thickness_act(j)
!        dummy6=thickness_act(j)*dummy5
!        y_actlay(m,j,upstream)=max(y_sec(m,j,upstream)-dummy6,y_original(m,j,upstream))
!write(*,*)j,m,dummy6
!      ENDDO

      DO m=2,npoints(j,upstream)-1
        IF (geom(m,j) == 1) THEN
          partarea_actlay(m,j,res_index(upstream))=(y_sec(m,j,res_index(upstream))-y_actlay(m,j,res_index(upstream)))* &
		      (x_sec(m,j,res_index(upstream))-x_sec(m-1,j,res_index(upstream)))/2.
          partarea_toplay(m,j,res_index(upstream))=(((elev-y_sec(m,j,res_index(upstream)))**2.)/ &
		      (2.*ABS(y_sec(m,j,res_index(upstream))-y_sec(m-1,j,res_index(upstream)))/(x_sec(m,j,res_index(upstream))-x_sec(m-1,j,res_index(upstream)))))
	      area_actlay(j,upstream)=area_actlay(j,upstream)+partarea_actlay(m,j,res_index(upstream))
	    ELSE IF (geom(m,j) == 2) THEN
          partarea_actlay(m,j,res_index(upstream))=((y_sec(m,j,res_index(upstream))-y_actlay(m,j,res_index(upstream)))+ &
			  (y_sec(m-1,j,res_index(upstream))-y_actlay(m-1,j,res_index(upstream))))*  &
              (x_sec(m,j,res_index(upstream))-x_sec(m-1,j,res_index(upstream)))/2.
          partarea_toplay(m,j,res_index(upstream))=(x_sec(m,j,res_index(upstream))-x_sec(m-1,j,res_index(upstream)))  &
              *(2.*elev-(y_sec(m,j,res_index(upstream))+y_sec(m-1,j,res_index(upstream))))/2.
          area_actlay(j,upstream)=area_actlay(j,upstream)+partarea_actlay(m,j,res_index(upstream))
	    ELSE IF (geom(m,j) == 3) THEN
          partarea_actlay(m,j,res_index(upstream))=((y_sec(m,j,res_index(upstream))-y_actlay(m,j,res_index(upstream)))+ &
			  (y_sec(m-1,j,res_index(upstream))-y_actlay(m-1,j,res_index(upstream))))*  &
              (x_sec(m,j,res_index(upstream))-x_sec(m-1,j,res_index(upstream)))/2.+ &
			  (y_sec(m,j,res_index(upstream))-y_actlay(m,j,res_index(upstream)))* &
		      (x_sec(m+1,j,res_index(upstream))-x_sec(m,j,res_index(upstream)))/2.
          partarea_toplay(m,j,res_index(upstream))=(x_sec(m,j,res_index(upstream))-x_sec(m-1,j,res_index(upstream)))  &
              *(2.*elev- (y_sec(m,j,res_index(upstream))+y_sec(m-1,j,res_index(upstream))))/2.+ &
			  (((elev-y_sec(m,j,res_index(upstream)))**2.)/ &
			  (2.*ABS(y_sec(m+1,j,res_index(upstream))-y_sec(m,j,res_index(upstream)))/(x_sec(m+1,j,res_index(upstream))-x_sec(m,j,res_index(upstream)))))
          area_actlay(j,upstream)=area_actlay(j,upstream)+partarea_actlay(m,j,res_index(upstream))
	    ELSE IF (geom(m,j) == 4) THEN
          partarea_actlay(m,j,res_index(upstream))=(y_sec(m,j,res_index(upstream))-y_actlay(m,j,res_index(upstream)))* &
		      (x_sec(m+1,j,res_index(upstream))-x_sec(m-1,j,res_index(upstream)))/2.
          partarea_toplay(m,j,res_index(upstream))=(((elev-y_sec(m,j,res_index(upstream)))**2.)/ &
			  (2.*ABS(y_sec(m,j,res_index(upstream))-y_sec(m-1,j,res_index(upstream)))/(x_sec(m,j,res_index(upstream))-x_sec(m-1,j,res_index(upstream)))))+ &
              (((elev-y_sec(m,j,res_index(upstream)))**2.)/ &
			  (2.*ABS(y_sec(m+1,j,res_index(upstream))-y_sec(m,j,res_index(upstream)))/(x_sec(m+1,j,res_index(upstream))-x_sec(m,j,res_index(upstream)))))
          area_actlay(j,upstream)=area_actlay(j,upstream)+partarea_actlay(m,j,res_index(upstream))
		ENDIF
	  ENDDO

	  DO m=1,npoints(j,upstream)
	    IF (area_sec(j,upstream) /= 0.) THEN
!		  weightfac_toplay(m,j,upstream)=partarea_toplay(m,j,upstream)/area_sec(j,upstream)
	    ELSE
!		  weightfac_toplay(m,j,upstream)=0.
	    ENDIF

	    IF (area_actlay(j,upstream) /= 0.) THEN
!	      weightfac_actlay(m,j,upstream)=partarea_actlay(m,j,upstream)/area_actlay(j,upstream)
	    ELSE
!	      weightfac_actlay(m,j,upstream)=0.
	    ENDIF
	  enddo

	ENDDO


    DO j=1,nbrsec(upstream)
      area_toplay(j,upstream)=area_sec(j,upstream)
	ENDDO

! Active layer volume (m3)
    DO j=1,nbrsec(upstream)
      vol_actlay(j,upstream)=area_actlay(j,upstream)*length_sec(j,upstream)
      vol_toplay(j,upstream)=area_toplay(j,upstream)*length_sec(j,upstream)
    END DO

! Determination of water volume of the reservoir reach by summing up all cross sections' volume
! that belong to it [m**3]
    resreach_vol(upstream)=0.
	p=0
!george    DO j=1,nbrsec(upstream)
!george      if (watelev_sec(j,upstream) /= damelev_mean(step,upstream)) p=j
!george	ENDDO
!george	if (p /=0 ) then
      DO j=1,nbrsec(upstream)
	    resreach_vol(upstream)=resreach_vol(upstream)+vol_toplay(j,upstream)
!george	    if (j > p) resreach_vol(upstream)=resreach_vol(upstream)+vol_toplay(j,upstream)
!write(*,'(2I4,4F15.3)')j,p,watelev_sec(j,upstream),damelev_mean(step,upstream),vol_toplay(j,upstream),resreach_vol(upstream)
	  enddo
!george	else
!george	  resreach_vol(upstream)=0.
!george	endif

! Fractional sediment availability at the active layer (ton)
    DO j=1,nbrsec(upstream)
      DO g=1,n_sed_class
        frsedavailab(g,j)=vol_actlay(j,upstream)*frac_actlay(g,j,res_index(upstream))*dry_dens(upstream)
!if(j<10)write(*,'(2I4,10F15.2)')j,g,frsedavailab(g,j),vol_actlay(j,upstream),frac_actlay(g,j,upstream),dry_dens(upstream)
      END DO
    END DO
    DO g=1,n_sed_class
	  retent(g)=0.
	enddo


! 2c) Sediment balance for each cross section (ton)

! Calculation of hydraulic parameter at the outlet of sub-reaches represented by each cross section
    do j=1,nbrsec(upstream)
      if (j /= nbrsec(upstream)) then
        discharge(j)=(discharge_sec(j,upstream)+discharge_sec(j+1,upstream))/2.
        topwidth(j)=(topwidth_sec(j,upstream)+topwidth_sec(j+1,upstream))/2.
		area_mean(j)=(area_sec(j,upstream)+area_sec(j+1,upstream))/2.
      else
        discharge(j)=discharge_sec(j,upstream)
        topwidth(j)=topwidth_sec(j,upstream)
		area_mean(j)=area_sec(j,upstream)
      endif
    enddo

	j=nbrsec(upstream)
!    dummy9=discharge_sec(j,upstream)/vol_toplay(j,upstream)
	dummy9=1.
!write(*,'(2I4,2F11.1,F11.6)')d,upstream,discharge_sec(j,upstream),vol_toplay(j,upstream),dummy9
    DO j=1,nbrsec(upstream)
      totalload(j,upstream)=0.

!	  call ssc_function(j,upstream,dummy9)

!XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
! Calculation of a representative discharge of each cross section based on area fraction
!      if (watelev_sec(j,upstream)>damelev_mean(step,upstream)) then
!	    discharge_mod(j)=discharge(j)
!	  else
!	    dummy4=discharge(j)/1.
!		if(dummy4/=0.) then
!		  if (area_mean(j)/dummy4>1.) discharge_mod(j)= &
!				discharge(j)* (area_mean(j)/dummy4)
!		  if (area_mean(j)/dummy4<=1.) &
!		        discharge_mod(j)=discharge(j)
!        else
!		  discharge_mod(j)=0.
!        endif
!      endif

!write(*,'(2I4,3F8.3,4F11.3)')d,j,watelev_sec(j,upstream),damelev_mean(step,upstream), &
! discharge(j),discharge_mod(j),area_mean(j), &
! normalarea(j),area_mean(j)/normalarea(j)

!write(*,'(2I4,3F8.3,4F11.3)')d,j,watelev_sec(j,upstream),damelev_mean(step,upstream), &
! discharge(j),discharge_mod(j),area_mean(j), &
! dummy4,area_mean(j)/dummy4
!XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
	  dummy5=0.
	  dummy6=0.
	  dummy7=0.
      DO g=1,n_sed_class
        frerosion(g,j)=0.
        frdeposition(g,j)=0.
        frretention(g,j)=0.
		frtotal_discharge(g,j)=0.



!New approach
!*************************************************************************
! Computation of fractional erosion for each cross section
        IF (j > 1) THEN
		  frac_toplay(g,j,res_index(upstream))=frac_toplay(g,j-1,res_index(upstream))
		ENDIF

! The weighting factor to include incoming sediment into the carrying capacity has to be read in do.dat
! temporarily selected as 0.9 according to the literature
        weigth_factor=0.9
		fr_capacity(g,j)=fr_capacity(g,j)*dry_dens(upstream)
        IF (fr_capacity(g,j) > 0.) THEN
          fr_capacity(g,j)=fr_capacity(g,j)*(weigth_factor*frac_actlay(g,j,res_index(upstream))+ &
			   (1.-weigth_factor)*frac_toplay(g,j,res_index(upstream)))
        ELSE
          fr_capacity(g,j)=0.
        END IF
!write(*,*)j,g,fr_capacity(g,j),frac_actlay(g,j,upstream),frac_toplay(g,j,upstream),(weigth_factor*frac_actlay(g,j,upstream)+ &
!			   (1.-weigth_factor)*frac_toplay(g,j,upstream))

! Conversion (ton/timestep to g/l)
!george		fr_capacity(g,j)=1000*fr_capacity(g,j)/(discharge(j)*(86400./nt))
!        IF (upper_limit(g) <= .0078) setvel(g)=0.000044
!S7        IF (upper_limit(g) <= .0040) setvel(g)=0.000012
!        IF (upper_limit(g) <= .0055) setvel(g)=0.000022
!write(*,'(2I4,2F15.6)')j,g,upper_limit(g),setvel(g)*1000
        IF (j == 1) THEN
          loadincoming=sed_inflow(step,upstream)*frac_toplay(g,j,res_index(upstream))
        ELSE
          loadincoming=frtotal_discharge(g,j-1)
        END IF
!write(*,*)j,g,fr_capacity(g,j),loadincoming
        if (discharge(j) /= 0.) then
	      IF (j == 1) THEN
!write(*,*)j,g,loadincoming
		    if (conc_inflow*frac_toplay(g,j,res_index(upstream))-fr_capacity(g,j)>=0.)recov=.25
		    if (conc_inflow*frac_toplay(g,j,res_index(upstream))-fr_capacity(g,j)<0.)recov=1.
!george		    frconc(g,j)=fr_capacity(g,j)+((conc_inflow*frac_toplay(g,j,upstream))-fr_capacity(g,j))*  &
!george					exp(-recov*setvel(g)*length_sec(j,upstream)*topwidth(j)/discharge(j))
!write(*,*)j,g,fr_capacity(g,j),loadincoming,recov,setvel(g),length_sec(j,upstream),topwidth(j),discharge(j)
		    frconc(g,j)=fr_capacity(g,j)+(loadincoming-fr_capacity(g,j))*  & !george
				exp(-recov*setvel(g)*length_sec(j,upstream)*topwidth(j)/discharge(j)) !george
		  else
		    if (frconc(g,j-1)-fr_capacity(g,j)>=0.)recov=.25
		    if (frconc(g,j-1)-fr_capacity(g,j)<0.)recov=1.
!write(*,*)j,g,fr_capacity(g,j),loadincoming,recov,setvel(g),length_sec(j,upstream),topwidth(j),discharge(j)
		    frconc(g,j)=fr_capacity(g,j)+(frconc(g,j-1)-fr_capacity(g,j))* &
				exp(-recov*setvel(g)*length_sec(j,upstream)*topwidth(j)/discharge(j))
		  endif
		else
		  frconc(g,j)=0.
		endif
!write(*,*)j,g,loadincoming,frconc(g,j)
!if(loadincoming<0. .or. frconc(g,j)<0.)stop

!write(*,'(2I4,2F10.1,2F10.6)')j,g,fr_capacity(g,j),frconc(g,j),exp(-recov*setvel(g)*length_sec(j,upstream)*topwidth(j)/discharge(j)),setvel(g)*length_sec(j,upstream)*topwidth(j)/discharge(j)

! Conversion (g/l to ton/timestep)
!george		fr_capacity(g,j)=frconc(g,j)*discharge(j)*(86400./nt)/1000
        fr_capacity(g,j)=frconc(g,j) !george

! Computation of remobilization of sediment at the cross section (fractional suspended load transport)
        IF (fr_capacity(g,j) > loadincoming) THEN
          frerosion(g,j)=MIN(frsedavailab(g,j),  &
              fr_capacity(g,j)-loadincoming)
		  frtotal_discharge(g,j)=loadincoming+frerosion(g,j)
        ELSE
          frerosion(g,j)=0.
        END IF

!if(j<10)write(*,*)step,j,g,frerosion(g,j),frsedavailab(g,j)


!        frsedavailab(g,j)=frsedavailab(g,j)-frerosion(g,j)
!        frsedavailab(g,j)=MAX(frsedavailab(g,j),0.)

! Computation of deposition at the cross section (fractional suspended load transport)
        IF (fr_capacity(g,j) < loadincoming .and. j/=nbrsec(upstream)) THEN
		  frtotal_discharge(g,j)=fr_capacity(g,j)
          frdeposition(g,j)=MAX(loadincoming-fr_capacity(g,j),0.)
        ELSE IF (fr_capacity(g,j) < loadincoming .and. j==nbrsec(upstream)) THEN
		  IF (fr_capacity(g,j) /= 0.) THEN
		    frtotal_discharge(g,j)=fr_capacity(g,j)
            frdeposition(g,j)=MAX(loadincoming-fr_capacity(g,j),0.)
		  ELSE
		  ENDIF
        ELSE
	      frdeposition(g,j)=0.
        END IF

		IF (fr_capacity(g,j) == loadincoming) frtotal_discharge(g,j)=loadincoming


!if (j>=52)write(*,*)step,j,g,loadincoming,frdeposition(g,j)
!george		dummy4=frtotal_discharge(g,j)
!george        frtotal_discharge(g,j)=frtotal_discharge(g,j)*dummy9
!george		if (frtotal_discharge(g,j) < dummy4) frdeposition(g,j)=frdeposition(g,j)+max(dummy4-frtotal_discharge(g,j),0.)



!write(*,'(2I4,5F13.3)')j,g,frconc(g,j),frtotal_discharge(g,j),frsedavailab(g,j),frerosion(g,j),loadincoming+frerosion(g,j)
!if (j==10)stop

! Calculation of sediment load in each cross section (g/l)
        if (discharge(j)/=0.) then
!george		  frconc(g,j)=frtotal_discharge(g,j)*1000./(discharge(j)*(86400./nt))
!write(*,'(2I4,5F20.6)')j,g,frconc(g,j),frtotal_discharge(g,j)
		  frconc(g,j)=frtotal_discharge(g,j) !george
!		  frconc(g,j)=frtotal_discharge(g,j)*1650/(discharge(j)*(86400./nt))

          if (j==nbrsec(upstream)) then
!temp		    frconc(g,j)=frtotal_discharge(g,j)*1000/(discharge_mod(j)*(86400./nt))
!temp		    frtotal_discharge(g,j)=frconc(g,j)*(discharge(j)*(86400./nt))/1000.
!temp            if (frtotal_discharge(g,j) >= loadincoming) THEN
!temp              retent(g)=0.
!temp			else
!temp              retent(g)=max(loadincoming-frtotal_discharge(g,j),0.)
!temp			endif
            if (resreach_vol(upstream) /= 0.) then
			  if (frtotal_discharge(g,j) /= 0.) then
!                call ssc_function(g,j,upstream,dummy9)
			    call vert_dist(g,j,upstream,factor_bottom,factor_intake,factor_over)
			  else
			    factor_bottom=0.
			    factor_intake=0.
			    factor_over=0.
				dummy9=0.
			  endif
			  if (qbottom(step,upstream)/=0.) factor_intake=max(factor_intake,.7)	!to account for resuspension of sediment by opening of the bottom outlets
!			    factor_bottom=1.
!			    factor_intake=1.
!			    factor_over=1.

!write(*,'(2I4,5F10.3)')j,g,factor_bottom,factor_intake,factor_over,dummy9

			  dummy4=frtotal_discharge(g,j)
              frtotal_discharge(g,j)=frtotal_discharge(g,j)*dummy9
			  if (frtotal_discharge(g,j) < dummy4) retent(g)=max(dummy4-frtotal_discharge(g,j),0.)
			  if (res_qout(step,upstream) /= 0.)then
 			    dummy5=factor_bottom*(qbottom(step,upstream)/res_qout(step,upstream))
			    dummy6=factor_intake*(qintake(step,upstream)/res_qout(step,upstream))
			    dummy7=factor_over*(overflow(step,upstream)/res_qout(step,upstream))
			  else
 			    dummy5=0.
 			    dummy6=0.
 			    dummy7=0.
			  endif
			  dummy4=frtotal_discharge(g,j)
              frtotal_discharge(g,j)=frtotal_discharge(g,j)* &
				(dummy5+dummy6+dummy7)
              if (frtotal_discharge(g,j) >= dummy4) THEN
!2007                retent(g)=0.
			  else
                retent(g)=retent(g)+max(dummy4-frtotal_discharge(g,j),0.)
			  endif
			  factor_bottom2=factor_bottom2+dummy5*dummy4
			  factor_intake2=factor_intake2+dummy6*dummy4
			  factor_over2=factor_over2+dummy7*dummy4
			else
			  retent(g)=0.
			  factor_bottom2=factor_bottom2+(qbottom(step,upstream)/res_qout(step,upstream))*frtotal_discharge(g,j)
			  factor_intake2=factor_intake2+(qintake(step,upstream)/res_qout(step,upstream))*frtotal_discharge(g,j)
			  factor_over2=factor_over2+(overflow(step,upstream)/res_qout(step,upstream))*frtotal_discharge(g,j)
			endif

!write(*,'(2I4,3F11.6,3F12.3)')j,g,dummy5,dummy6,dummy7,factor_bottom2,factor_intake2,factor_over2
!temp		    frconc(g,j)=frtotal_discharge(g,j)*1000/(discharge(j)*(86400./nt))
!write(*,'(I4,4F15.3)')j,factor_over2+factor_intake2+factor_bottom2,dummy4,frtotal_discharge(g,j),retent(g)
!write(*,'(2I4,5F10.3)')j,g,loadincoming,frdeposition(g,j),frerosion(g,j),frtotal_discharge(g,j),retent(g)
		  endif
	      totalload(j,upstream)=totalload(j,upstream)+frtotal_discharge(g,j)
		else
          if (j==nbrsec(upstream)) then
		    retent(g)=loadincoming
			frdeposition(g,j)=0.
!george		    retent(g)=0.
		  endif
		  frconc(g,j)=0.
	      totalload(j,upstream)=0.
		endif
!if(j<20)write(*,'(2I4,6F11.3)')j,g,frconc(g,j),frtotal_discharge(g,j),frsedavailab(g,j),frerosion(g,j),frdeposition(g,j),retent(g)

      END DO

	  IF (j==nbrsec(upstream)) THEN
!write(*,'(2I4,6F15.6)')j,g,factor_bottom2,factor_intake2,factor_over2,factor_bottom2+factor_intake2+factor_over2,factor_intake2/(factor_bottom2+factor_intake2+factor_over2)
        IF (totalload(j,upstream) /= 0.) THEN
		  dummy4=factor_bottom2+factor_intake2+factor_over2
	      factor_bottom2=factor_bottom2/dummy4
	      factor_intake2=factor_intake2/dummy4
	      factor_over2=factor_over2/dummy4
	    ELSE !Till: compute partitioning of sediment outflow (?)
	      if (res_qout(step,upstream) ==0) then
	        factor_bottom2=0.
            factor_intake2=0.
            factor_over2  =0.
          else
	        factor_bottom2=qbottom(step,upstream)/res_qout(step,upstream)
	        factor_intake2=qintake(step,upstream)/res_qout(step,upstream)
	        factor_over2=overflow(step,upstream)/res_qout(step,upstream)
	      end if
	    ENDIF
!write(*,'(2I4,6F12.6)')j,g,factor_bottom2,factor_intake2,factor_over2,factor_over2+factor_intake2+factor_bottom2
	  ENDIF

      DO g=1,n_sed_class
        IF (totalload(j,upstream) /= 0.) THEN
	      frac_toplay(g,j,res_index(upstream))=frtotal_discharge(g,j)/totalload(j,upstream)
	    ELSE
	      frac_toplay(g,j,res_index(upstream))=0.
	    ENDIF
!		write(*,*)j,g,frac_toplay(g,j,upstream),frtotal_discharge(g,j),totalload(j,upstream)
	  ENDDO
    END DO

!if(step==290)stop
!if (step==2)stop

! Distribution of the material retained near the dam along the reservoir reach
! using a weighting factor based on the water volume represented by each cross section
    weight=0.
    dummy5=0.
    DO j=1,nbrsec(upstream)
      if (j > p .and. resreach_vol(upstream) /= 0.) then
		if (j /= nbrsec(upstream)) then
          weight=vol_toplay(j,upstream)/resreach_vol(upstream)
		else
          weight=max(0.,1.-dummy5)
		endif
        DO g=1,n_sed_class
		  frretention(g,j)=weight*retent(g)
		enddo
	  else
		weight=0.
        DO g=1,n_sed_class
		  frretention(g,j)=0.
		enddo
      END IF
      dummy5=dummy5+weight
!write(*,'(I4,4F15.3)')j,weight,dummy5,vol_toplay(j,upstream),resreach_vol(upstream)
!write(*,'(I4,<n_sed_class>F8.0)')j,(frretention(g,j),g=1,n_sed_class)
    END DO


! Calculation of suspended sediment concentration for each cross section (g/l)
    DO j=1,nbrsec(upstream)
	  if (discharge(j) /= 0.) then
!temp        if (j/=nbrsec(upstream)) then
!temp	      conc(j,upstream)=(totalload(j,upstream)*1000)/(discharge_mod(j)*(86400./nt))
!temp		else
!temp	      conc(j,upstream)=(totalload(j,upstream)*1000)/(discharge(j)*(86400./nt))
!temp		endif
        conc(j,upstream)=(totalload(j,upstream)*1000.)/(discharge(j)*(86400./nt))
	  else
	    conc(j,upstream)=0.
	  endif
    END DO
!*****************************************************************************

! Calculation of the actual sediment volume at the active layer for each cross section
    DO j=1,nbrsec(upstream)
	  totvol_actlay(j,upstream)=0.
      DO g=1,n_sed_class
	    frvol_actlay(g,j,upstream)=frvol_actlay(g,j,upstream)+ &
			(frdeposition(g,j)-frerosion(g,j)+frretention(g,j))/dry_dens(upstream)
		frvol_actlay(g,j,upstream)=MAX(frvol_actlay(g,j,upstream),0.)
        totvol_actlay(j,upstream)=totvol_actlay(j,upstream)+frvol_actlay(g,j,upstream)
!write(*,'(2I4,5F15.4)')id_subbas_extern(upstream),j,frvol_actlay(g,j,upstream),totvol_actlay(j,upstream)
      END DO
    END DO

! Computation of the actual bed composition
    DO j=1,nbrsec(upstream)
      DO g=1,n_sed_class
        IF (totvol_actlay(j,upstream) /= 0.) THEN
          frac_actlay(g,j,res_index(upstream))=frvol_actlay(g,j,upstream)/totvol_actlay(j,upstream)
        ELSE
          frac_actlay(g,j,res_index(upstream))=0.
        END IF
      END DO
    END DO

! Computation of total erosion and total deposition for each cross section (ton/timestep)
    DO j=1,nbrsec(upstream)
      erosion(j,upstream)=0.
      deposition(j,upstream)=0.
	  retention(j,upstream)=0.
      DO g=1,n_sed_class
        erosion(j,upstream)=erosion(j,upstream)+frerosion(g,j)
        deposition(j,upstream)=deposition(j,upstream)+ frdeposition(g,j)
		retention(j,upstream)=retention(j,upstream)+frretention(g,j)
!write(*,'(2I4,6F11.3)')j,g,frsedavailab(g,j),frerosion(g,j),erosion(j,upstream),deposition(j,upstream),retention(j,upstream)
      END DO
    END DO

! Computation of sediment volume variation at the reservoir bed for each cross section
    DO j=1,nbrsec(upstream)
	  dvol_sed(j,upstream)=(deposition(j,upstream)-erosion(j,upstream)+retention(j,upstream))/dry_dens(upstream)
!write(*,'(I4,4F15.3)')j,deposition(j,upstream),erosion(j,upstream),retention(j,upstream),dvol_sed(j,upstream)
    END DO
!if(step==290)stop
! Computation of sediment area variation at the reservoir bed for each cross section
    DO j=1,nbrsec(upstream)
      darea_sed(j,upstream)=abs(dvol_sed(j,upstream)/  &
            length_sec(j,upstream))
    END DO

! Daily sediment retention in the subbasin's reservoir [ton/timestep]
    DO j=1,nbrsec(upstream)
      sed_ret(step,upstream)=sed_ret(step,upstream)+retention(j,upstream)
	enddo

! Cross section geometry before bed elevtion change
    DO j=1,nbrsec(upstream)
      DO m=1,npoints(j,upstream)
        y_laststep(m,j,res_index(upstream))=y_sec(m,j,res_index(upstream))
	  ENDDO
	ENDDO

! 2d) Calculation of reservoir bed elevation changes
    DO j=1,nbrsec(upstream)
      IF (dvol_sed(j,upstream) /= 0.) THEN
	    CALL change_sec(j,upstream)
	  ENDIF
	ENDDO

! Deposition area of each cross section (m2)
!2010    DO j=1,nbrsec(upstream)
!2010	  area_sedim(j,upstream)=0.
!2010      DO m=2,npoints(j,upstream)
!2010        area_sedim(j,upstream)=area_sedim(j,upstream)+  &
!2010          (x_sec(m,j,upstream)-x_sec(m-1,j,upstream))*  &
!2010 		  ((y_sec(m-1,j,upstream)-y_original(m-1,j,upstream))+ &
!2010 		  (y_sec(m,j,upstream)-y_original(m,j,upstream)))/2.
!2010	  enddo
!2010	enddo

! Calculation of sediment volume at the active layer of each cros section
!2010    DO j=1,nbrsec(upstream)
!2010      vol_sedim(j,upstream)=area_sedim(j,upstream)*length_sec(j,upstream)

! Check sediment volume discrepancy between values computed by summing up sediment deposition from each
! grain size and computed by comparison between original bed elevation and current bed elevation
!      dummy4=totvol_actlay(j,upstream)/length_sec(j,upstream)
!      if(abs(dummy4-area_sedim(j,upstream))>.001) then
!        write(*,*)'******************************************************************'
!        write(*,'(I4,4F15.6)')j,dummy4,area_sedim(j,upstream),darea_sed(j,upstream),dvol_sed(j,upstream)
!        write(*,*)'******************************************************************'
!      endif
!2010    END DO

! Sediment outflow (ton/time step)
    j=nbrsec(upstream)
	if (totalload(j,upstream) /= 0.) then
	  sed_overflow(step,upstream)=totalload(j,upstream)*factor_over2
	  sed_intake(step,upstream)=totalload(j,upstream)*factor_intake2
	  sed_bottom(step,upstream)=totalload(j,upstream)*factor_bottom2
	  sed_outflow(step,upstream)=totalload(j,upstream)
!temp	  sed_overflow(step,upstream)=max(0.,totalload(j,upstream)*overflow(step,upstream)/(res_qout(step,upstream)))
!temp	  sed_intake(step,upstream)=max(0.,totalload(j,upstream)*qintake(step,upstream)/(res_qout(step,upstream)))
!temp	  sed_bottom(step,upstream)=max(0.,totalload(j,upstream)*qbottom(step,upstream)/(res_qout(step,upstream)))
!temp	  sed_outflow(step,upstream)=sed_overflow(step,upstream)+sed_intake(step,upstream)+sed_bottom(step,upstream)
	else
	  sed_overflow(step,upstream)=0.
	  sed_intake(step,upstream)=0.
	  sed_bottom(step,upstream)=0.
	  sed_outflow(step,upstream)=0.
	endif

! ***********************************************************
! Compaction of the deposited material has to be calculated
! ***********************************************************

! Calculation of daily sedimentation and cumulative sedimentation (m3)
!2010    cum_sedimentation(upstream)=0.
	sedimentation(step,upstream)=0.
    DO j=1,nbrsec(upstream)
!2010      cum_sedimentation(upstream)=cum_sedimentation(upstream)+vol_sedim(j,upstream)
	  sedimentation(step,upstream)=sedimentation(step,upstream)+dvol_sed(j,upstream)
	ENDDO
!write(*,*)cum_sedimentation(upstream)
! Revised cumulative sedimentation, i.e. sedimentation above the initial bed elevation (m3)
!2010    cum_sedimentation(upstream)=cum_sedimentation(upstream)-volbed0(upstream)
!write(*,'(I4,3F13.1)')upstream,cum_sedimentation(upstream),volbed0(upstream)
! Storage capacity reduction (m3)
    decstorcap(step,upstream)=sedimentation(step,upstream)

! Conversion (m3 to ton)
!2010    cum_sedimentation(upstream)=cum_sedimentation(upstream)*dry_dens(upstream)
	sedimentation(step,upstream)=sedimentation(step,upstream)*dry_dens(upstream)
	if (step==1 .and. t==tstart) cum_sedimentation(upstream)=cum_sedimentation(upstream)+volbed0(upstream) !2010
	if (t>tstart .and. t== damyear(upstream) .and. step==1) cum_sedimentation(upstream)=cum_sedimentation(upstream)+volbed0(upstream) !2010
    cum_sedimentation(upstream)=cum_sedimentation(upstream)+sedimentation(step,upstream) !2010
!write(*,'(2F18.3)')cum_sedimentation(upstream),sedimentation(step,upstream)
!if (step==30)stop

! Calculation of total sediment release (ton)
	j=nbrsec(upstream)
    DO g=1,n_sed_class
	  frsediment_out(upstream,g)=frac_toplay(g,j,res_index(upstream))
	enddo


! Change on stage-area-volume curve due to erosion and deposition processes
    dummy9=0.
	nbrbat1=nbrbat(upstream)
	if (decstorcap(step,upstream) < daystorcap(step,upstream) .and. decstorcap(step,upstream)/=0.) then
     DO b=1,nbrbat1
      elevhelp=elev_bat(b,upstream)
      DO j=1,nbrsec(upstream)
        areasec(j)=0.
        widthsec(j)=0.
        areasec2(j)=0.
        widthsec2(j)=0.

        npt=npoints(j,upstream)
        DO m=2,npt
		  IF (elevhelp >= y_sec(m-1,j,res_index(upstream))  &
			.AND.elevhelp >= y_sec(m,j,res_index(upstream))) THEN
		    areasec(j)=areasec(j)+(x_sec(m,j,res_index(upstream))-x_sec(m-1,j,res_index(upstream)))  &
				*(elevhelp-(y_sec(m,j,res_index(upstream))+y_sec(m-1,j,res_index(upstream)))/2.)
			widthsec(j)=widthsec(j)+  &
				(x_sec(m,j,res_index(upstream))-x_sec(m-1,j,res_index(upstream)))
			IF (fcav(upstream)==0) THEN
			 areasec2(j)=areasec2(j)+(x_sec(m,j,res_index(upstream))-x_sec(m-1,j,res_index(upstream)))  &
				*(elevhelp-(y_laststep(m,j,res_index(upstream))+y_laststep(m-1,j,res_index(upstream)))/2.)
			 widthsec2(j)=widthsec2(j)+  &
				(x_sec(m,j,res_index(upstream))-x_sec(m-1,j,res_index(upstream)))
			ENDIF
		  ELSE IF (elevhelp < y_sec(m-1,j,res_index(upstream))  &
			.AND.elevhelp >= y_sec(m,j,res_index(upstream))) THEN
			areasec(j)=areasec(j)+(((elevhelp-y_sec(m,j,res_index(upstream)))**2.)/ &
				(2.*ABS(y_sec(m,j,res_index(upstream))-y_sec(m-1,j,res_index(upstream)))/  &
				(x_sec(m,j,res_index(upstream))-x_sec(m-1,j,res_index(upstream)))))
			widthsec(j)=widthsec(j)+  &
				((x_sec(m,j,res_index(upstream))-x_sec(m-1,j,res_index(upstream)))*  &
				(elevhelp-y_sec(m,j,res_index(upstream)))/  &
				(y_sec(m-1,j,res_index(upstream))-y_sec(m,j,res_index(upstream))))
			IF (fcav(upstream)==0) THEN
			 areasec2(j)=areasec2(j)+(((elevhelp-y_laststep(m,j,res_index(upstream)))**2.)/ &
				(2.*ABS(y_laststep(m,j,res_index(upstream))-y_laststep(m-1,j,res_index(upstream)))/  &
				(x_sec(m,j,res_index(upstream))-x_sec(m-1,j,res_index(upstream)))))
			 widthsec2(j)=widthsec2(j)+  &
				((x_sec(m,j,res_index(upstream))-x_sec(m-1,j,res_index(upstream)))*  &
				(elevhelp-y_laststep(m,j,res_index(upstream)))/  &
				(y_laststep(m-1,j,res_index(upstream))-y_laststep(m,j,res_index(upstream))))
			ENDIF
		  ELSE IF (elevhelp >= y_sec(m-1,j,res_index(upstream))  &
			.AND.elevhelp < y_sec(m,j,res_index(upstream))) THEN
			areasec(j)=areasec(j)+(((elevhelp-y_sec(m-1,j,res_index(upstream)))**2.)/ &
				(2.*ABS(y_sec(m,j,res_index(upstream))-y_sec(m-1,j,res_index(upstream)))/  &
				(x_sec(m,j,res_index(upstream))-x_sec(m-1,j,res_index(upstream)))))
			widthsec(j)=widthsec(j)+  &
				((x_sec(m,j,res_index(upstream))-x_sec(m-1,j,res_index(upstream)))*  &
				(elevhelp-y_sec(m-1,j,res_index(upstream)))/  &
				(y_sec(m,j,res_index(upstream))-y_sec(m-1,j,res_index(upstream))))
			IF (fcav(upstream)==0) THEN
			 areasec2(j)=areasec2(j)+(((elevhelp-y_laststep(m-1,j,res_index(upstream)))**2.)/ &
				(2.*ABS(y_laststep(m,j,res_index(upstream))-y_laststep(m-1,j,res_index(upstream)))/  &
				(x_sec(m,j,res_index(upstream))-x_sec(m-1,j,res_index(upstream)))))
			 widthsec2(j)=widthsec2(j)+  &
				((x_sec(m,j,res_index(upstream))-x_sec(m-1,j,res_index(upstream)))*  &
				(elevhelp-y_laststep(m-1,j,res_index(upstream)))/  &
				(y_laststep(m,j,res_index(upstream))-y_laststep(m-1,j,res_index(upstream))))
			ENDIF
		  END IF
	    END DO
	   END DO

       volhelp(b)=0.
       areahelp(b)=0.
       volhelp2(b)=0.
       areahelp2(b)=0.

       DO j=1,nbrsec(upstream)
        IF (j /= nbrsec(upstream)) THEN
	      volhelp(b)=volhelp(b)+(areasec(j)+areasec(j+1))*dist_sec(j,upstream)/2.
	      areahelp(b)=areahelp(b)+(widthsec(j)+widthsec(j+1))*dist_sec(j,upstream)/2.
	      IF (fcav(upstream)==0) volhelp2(b)=volhelp2(b)+(areasec2(j)+areasec2(j+1))*dist_sec(j,upstream)/2.
	      IF (fcav(upstream)==0) areahelp2(b)=areahelp2(b)+(widthsec2(j)+widthsec2(j+1))*dist_sec(j,upstream)/2.
        ELSE
	      volhelp(b)=volhelp(b)+areasec(j)*dist_sec(j,upstream)/2.
	      areahelp(b)=areahelp(b)+widthsec(j)*dist_sec(j,upstream)/2.
	      IF (fcav(upstream)==0) volhelp2(b)=volhelp2(b)+areasec2(j)*dist_sec(j,upstream)/2.
	      IF (fcav(upstream)==0) areahelp2(b)=areahelp2(b)+widthsec2(j)*dist_sec(j,upstream)/2.
        END IF
	   END DO
     ENDDO

! Calculation of the ratio between the reservoir volume calculated according to the cross sections and that value given in the file cav.dat
	 nbrbat1=nbrbat(upstream)
     DO b=1,nbrbat1
	   if (elev_bat(b,upstream) >= maxlevel(upstream)) THEN
!	     if (step==1 .and. t==tstart) dummy12=vol_bat(b,upstream)/volhelp2(b)
!	     if (t>tstart .and. t== damyear(upstream) .and. step==1) dummy12=vol_bat(b,upstream)/volhelp2(b)
!	     if (step==1 .and. t==tstart) dummy13=area_bat(b,upstream)/areahelp2(b)
!	     if (t>tstart .and. t== damyear(upstream) .and. step==1) dummy13=area_bat(b,upstream)/areahelp2(b)
		 IF (fcav(upstream)==0) THEN
!2010		   dummy12=vol_bat(b,upstream)/volhelp2(b)
!2010		   dummy13=area_bat(b,upstream)/areahelp2(b)
!2010		   fcav(upstream)=1
		   exit
		 ENDIF
	   endif
!write(*,*)b,nbrbat1,elev_bat(b,upstream),maxlevel(upstream)
!write(*,'(2I4,3F13.1,2F13.3)')step,b,volhelp2(b),volhelp(b),volhelp2(b)-volhelp(b),dummy12,dummy13
!write(*,'(2I4,3F13.1,2F13.3)')step,b,volhelp2(b),volhelp(b),areahelp2(b)-areahelp(b),dummy12,dummy13
!write(*,'(4I4,5F13.3)')d,hour,b,fcav(upstream),volhelp2(b),volhelp(b),vol_bat(b,upstream),dummy12,dummy13
!if(d>199 .and. hour==24)write(*,'(3I4,F12.8,4F13.3,F9.4)')d,hour,b,dummy10(b),volhelp2(b),volhelp(b),vol_bat(b,upstream)
!write(*,'(3I4,F10.3,F12.8,5F13.3)')d,hour,b,elev_bat(b,upstream),dummy10(b),volhelp(b),volhelp2(b),vol_bat(b,upstream)
     ENDDO
!if (d==95)stop
!if (d==2)stop
!write(*,*)dummy12,dummy13,fcav(upstream)

     DO b=1,nbrbat1
	   if (elev_bat(b,upstream) >= maxlevel(upstream) .and. volhelp(b) /= 0.) THEN
!2010         dummy9=(vol_bat(b,upstream)-decstorcap(step,upstream))/(dummy12*volhelp(b))
!write(*,'(2I4,3F13.1,2F13.8)')step,b,vol_bat(b,upstream),decstorcap(step,upstream),volhelp(b),dummy12,dummy9
		 exit
	   ENDIF
	 ENDDO

     DO b=1,nbrbat(upstream)
	   if (elev_bat(b,upstream) >= maxlevel(upstream)) THEN
	     vol_bat(b,upstream)=max(vol_bat(b,upstream)-decstorcap(step,upstream),0.)
	   else
!	     dummy3=vol_bat(b,upstream)
!write(*,'(2I4,3F13.1,2F13.3)')step,b,volhelp(b),dummy12*volhelp(b),vol_bat(b,upstream),dummy12,dummy9
!2010	     vol_bat(b,upstream)=dummy12*volhelp(b)*dummy9
!2010	     area_bat(b,upstream)=dummy13*areahelp(b)
	     vol_bat(b,upstream)=volhelp(b)
	     area_bat(b,upstream)=areahelp(b)
	   endif
!write(*,'(2I4,3F13.1,2F13.3)')step,b,volhelp(b),dummy12*volhelp(b),vol_bat(b,upstream),dummy12,dummy9
!write(*,'(2I4,4F13.1,2F13.3)')step,b,vol_bat(b,upstream),dummy3,vol_bat(b,upstream)-dummy3,decstorcap(step,upstream)
	 ENDDO
!if (step==21)stop

! Minimum elevation for each cross_section (m)
     j=nbrsec(upstream)
     dummy5=minelev_sec(j,upstream)
	 minelev_sec(j,upstream)=y_sec(1,j,res_index(upstream))
	 npt=npoints(j,upstream)
	 DO m=2,npt
	   IF (minelev_sec(j,upstream) > y_sec(m,j,res_index(upstream))) THEN
		 minelev_sec(j,upstream)=y_sec(m,j,res_index(upstream))
	   END IF
	 END DO
	 dummy4=minelev_sec(j,upstream)-dummy5
	 dayminlevel(step,upstream)=min(dayminlevel(step,upstream)+dummy4,maxlevel(upstream))

!write(*,'(2I4,3F15.3)')d,p,dummy4,dummy12,dayminlevel(step,upstream)
!if (step==126)stop

	 if (dayminlevel(step,upstream) < elev_bat(1,upstream)) elev_bat(1,upstream)=dayminlevel(step,upstream)
	 if (dayminlevel(step,upstream) > elev_bat(nbrbat(upstream),upstream)) elev_bat(nbrbat(upstream),upstream)=dayminlevel(step,upstream)

     DO b=1,nbrbat(upstream)
	   if (elev_bat(b,upstream) < dayminlevel(step,upstream)) then
	     vol_bat(b,upstream)=0.
		 area_bat(b,upstream)=0.
	   endif
	 enddo


	ENDIF


	IF (decstorcap(step,upstream) >= daystorcap(step,upstream) .or.&
	    cum_sedimentation(upstream)/dry_dens(upstream) >= storcap(upstream)*1.e6) then
	 storcap(upstream)=0.
	 daystorcap(step,upstream)=0.

     DO b=1,nbrbat(upstream)
	  vol_bat(b,upstream)=0.
	  area_bat(b,upstream)=0.
	 ENDDO
	ENDIF



!Ge output file for check some parameters (temporary)
!***************************************************
!OPEN(11,FILE=pfadn(1:pfadi)//'check.out', &
! 		STATUS='old',POSITION='append')
! DO j=1,nbrsec(upstream)
! write(11,'(3I4,7F20.6)')t,d,id_sec_extern(j,upstream),conc(j,upstream),totalload(j,upstream),dvol_sed(j,upstream),area_actlay(j,upstream),thickness_act(j),discharge(j),vol_toplay(j,upstream)
! write(11,'(4I4,20F20.6)')t,d,id_sec_extern(j,upstream),pt_long(upstream),totvol_actlay(j,upstream),totvol_actlay0(j,upstream),totvol_actlay(j,upstream)-totvol_actlay0(j,upstream)
!  enddo
! CLOSE(11)
! OPEN(11,FILE=pfadn(1:pfadi)//'check1.out', &
! 		STATUS='old',POSITION='append')

!write(fmtstr,'(a,i0,a)')'(3I4,',n_sed_class,'F15.6)'		!generate format string

!Ge write(11,'(2I4,<nbrsec(upstream)>F15.6)')t,d,(dvol_sed(j,upstream)*dry_dens(upstream),j=1,nbrsec(upstream))
! DO j=1,nbrsec(upstream)
	 !write(11,'(3I4,<n_sed_class>F20.6)')t,d,id_sec_extern(j,upstream),(frvol_actlay(g,j,upstream)-frvol_actlay0(g,j,upstream),g=1,n_sed_class)
!	 write(11,fmtstr)t,d,id_sec_extern(j,upstream),(frvol_actlay(g,j,upstream)-frvol_actlay0(g,j,upstream),g=1,n_sed_class)

	!Ge if (j==nbrsec(upstream)) then
	!write(11,'(2I4,3F20.6)')t,d,discharge_sec(j,upstream),discharge_mod(j),vol_toplay(j,upstream)
	!Ge endif
	!write(11,'(3I4,<n_sed_class>F20.6)')t,d,id_sec_extern(j,upstream),(frac_actlay(g,j,upstream),g=1,n_sed_class)
	!write(11,'(3I4,<n_sed_class>F10.4)')t,d,id_sec_extern(j,upstream),(frac_toplay(g,j,upstream),g=1,n_sed_class)
! enddo
! CLOSE(11)


! OPEN(11,FILE=pfadn(1:pfadi)//'check2.out', &
! 		STATUS='old',POSITION='append')
! write(fmtstr,'(a,i0,a)')'(3I4,',n_sed_class,'F20.3)'		!generate format string

! DO j=1,nbrsec(upstream)
!	write(11,fmtstr)t,d,id_sec_extern(j,upstream),(fr_capacity(g,j),g=1,n_sed_class)
	!write(11,'(3I4,<n_sed_class>F20.3)')t,d,id_sec_extern(j,upstream),(fr_capacity(g,j),g=1,n_sed_class)
! enddo
! DO j=1,nbrsec(upstream)
! write(11,'(3I4,<n_sed_class>F15.6)')t,d,id_sec_extern(j,upstream),(frtotal_discharge(g,j),g=1,n_sed_class)
! enddo
! CLOSE(11)

! OPEN(11,FILE=pfadn(1:pfadi)//'check3.out', &
! 		STATUS='old',POSITION='append')
! 		STATUS='unknown')
! if (t==tstop .and. d==dayyear) then
!if (t==1986 .and. d==5) then
!Ge write(11,*)'1993 computed'



! DO j=1,nbrsec(upstream)
!	 npt=npoints(j,upstream)
!	 write(fmtstr,'(a,i0,a)')'(2I4,',npt,'F15.6)'		!generate format string
!	 write(11,'(2I4,<npt>F15.6)')id_sec_extern(j,upstream),npt,(x_sec(m,j,upstream),m=1,npt)
!	 write(11,'(2I4,<npt>F15.6)')id_sec_extern(j,upstream),npt,(y_sec(m,j,upstream),m=1,npt)
!	 write(11,fmtstr)id_sec_extern(j,upstream),npt,(x_sec(m,j,upstream),m=1,npt)
!	 write(11,fmtstr)id_sec_extern(j,upstream),npt,(y_sec(m,j,upstream),m=1,npt)
	!ge write(11,'(3I4,<n_sed_class>F15.6)')t,d,id_sec_extern(j,upstream),(frvol_actlay(g,j,upstream),g=1,n_sed_class)
! enddo
!stop
! endif

!Ge CLOSE(11)
!***************************************************


    IF (reservoir_print == 0) THEN
! Print results on hydraulic calculations
     WRITE(subarea,*)id_subbas_extern(upstream)
	 IF (f_res_hydraul) THEN
     OPEN(11,FILE=pfadn(1:pfadi)//'res_'//trim(adjustl(subarea))//'_hydraul.out',STATUS='old',  &
		  POSITION='append')
      DO j=1,nbrsec(upstream)
        WRITE(11,'(5I6,4F15.3,E15.3E2,3F15.6)')id_subbas_extern(upstream),t,d,hour,id_sec_extern(j,upstream),  &
			depth_sec(j,upstream),watelev_sec(j,upstream),  &
			area_sec(j,upstream),topwidth_sec(j,upstream),  &
			energslope_sec(j,upstream),hydrad_sec(j,upstream),  &
			meanvel_sec(j,upstream),discharge_sec(j,upstream)
	  END DO
	 ENDIF
     CLOSE(11)

! Print results on bed elevation change of each cross section of the sub-basin's reservoir
     DO j=1,nbrsec(upstream)
	  npt=npoints(j,upstream)
      WRITE(section,*)id_sec_extern(j,upstream)
	  IF (f_res_bedchange) THEN
	  OPEN(11,FILE=pfadn(1:pfadi)//'res_'//trim(adjustl(subarea))//'_sec'//trim(adjustl(section))// &
			'_bedchange.out',STATUS='old',POSITION='append')
         write(fmtstr,'(a,i0,a)')'(6I6,',npt,'F15.6)'		!generate format string
		 WRITE(11,fmtstr)id_subbas_extern(upstream),id_sec_extern(j,upstream),t,d,hour,npt,  &
			(y_sec(m,j,res_index(upstream)),m=1,npoints(j,upstream))
		 !WRITE(11,'(6I6,<npt>F15.6)')id_subbas_extern(upstream),id_sec_extern(j,upstream),t,d,hour,npt,  &
		!	(y_sec(m,j,upstream),m=1,npoints(j,upstream))
      CLOSE(11)
	 ENDIF
	 ENDDO

! Print results on bed elevation change of the longitudinal profile of the sub-basin's reservoir
     nbrsec1=nbrsec(upstream)
     IF (f_res_longitudunal) THEN
	 OPEN(11,FILE=pfadn(1:pfadi)//'res_'//trim(adjustl(subarea))//'_longitudunal.out', &
		STATUS='old',POSITION='append')
       write(fmtstr,'(a,i0,a)')'(5I6,',nbrsec1,'F15.6)'		!generate format string
	   WRITE(11,fmtstr)id_subbas_extern(upstream),t,d,hour,nbrsec1, &
			(minelev_sec(j,upstream),j=1,nbrsec1)
	   !WRITE(11,'(5I6,<nbrsec1>F15.6)')id_subbas_extern(upstream),t,d,hour,nbrsec1, &
		!	(minelev_sec(j,upstream),j=1,nbrsec1)
     CLOSE(11)
	 ENDIF

	ELSE
      DO j=1,nbrsec(upstream)
	  	daydepth_sec(step,j,res_index(upstream))=depth_sec(j,upstream)
	  	daywatelev_sec(step,j,res_index(upstream))=watelev_sec(j,upstream)
	  	dayarea_sec(step,j,res_index(upstream))=area_sec(j,upstream)
	  	daytopwidth_sec(step,j,res_index(upstream))=topwidth_sec(j,upstream)
	  	dayenergslope_sec(step,j,res_index(upstream))=energslope_sec(j,upstream)
	  	dayhydrad_sec(step,j,res_index(upstream))=hydrad_sec(j,upstream)
	  	daymeanvel_sec(step,j,res_index(upstream))=meanvel_sec(j,upstream)
	  	daydischarge_sec(step,j,res_index(upstream))=discharge_sec(j,upstream)
		dayminelev_sec(step,j,res_index(upstream))=minelev_sec(j,upstream)
		DO m=1,npoints(j,upstream)
		  dayy_sec(step,m,j,res_index(upstream))=y_sec(m,j,res_index(upstream))
		ENDDO
	  ENDDO
	ENDIF

  END IF !Case 2 (with cross sections)

! Calculation of fractional sediment transport to the next dowstream sub-basin (ton/timestep)
  DO g=1,n_sed_class
    res_sediment_out(upstream,g)=frsediment_out(upstream,g)*sed_outflow(step,upstream)
  enddo


! Print results on bed elevation change of the longitudinal profile of the sub-basin's reservoir
  IF (reservoir_print == 0) THEN
   WRITE(subarea,*)id_subbas_extern(upstream)
   j=nbrsec(upstream)
   IF (f_res_sedbal) THEN
   OPEN(11,FILE=pfadn(1:pfadi)//'res_'//trim(adjustl(subarea))//'_sedbal.out',STATUS='old',POSITION='append')
    WRITE(11,'(4I6,4F15.3)')id_subbas_extern(upstream),t,d,hour,sed_inflow(step,upstream),sed_outflow(step,upstream), &
		sedimentation(step,upstream),cum_sedimentation(upstream)
   CLOSE(11)
   ENDIF

! Print results on sediment composition released out the sub-basin's reservoir
   IF (f_res_sedcomposition) THEN
   OPEN(11,FILE=pfadn(1:pfadi)//'res_'//trim(adjustl(subarea))//'_sedcomposition.out',STATUS='old',POSITION='append')
    write(fmtstr,'(a,i0,a)')'(5I6,',n_sed_class,'F15.3)'		!generate format string
	WRITE(11,fmtstr)id_subbas_extern(upstream),t,d,hour,n_sed_class, &
		(frsediment_out(upstream,g),g=1,n_sed_class)
	!WRITE(11,'(5I6,<n_sed_class>F15.3)')id_subbas_extern(upstream),t,d,hour,n_sed_class, &
	!	(frsediment_out(upstream,g),g=1,n_sed_class)
   CLOSE(11)
   ENDIF
  ELSE
   daycumsed(step,res_index(upstream))=cum_sedimentation(upstream)
   DO g=1,n_sed_class
     dayfrsediment_out(step,res_index(upstream),g)=frsediment_out(upstream,g)
   ENDDO
  ENDIF


!Ge end of status 2
END IF
! -----------------------------------------------------------------------
IF (STATUS == 3) THEN

! Output files of Reservoir Modules
  IF (reservoir_print == 1) THEN
    DO i=1,subasin
     IF (storcap(i) /= 0. .and. t >= damyear(i)) THEN
      WRITE(subarea,*)id_subbas_extern(i)
      IF (nbrsec(i) /= 0) THEN
        IF (f_res_hydraul) THEN
		OPEN(11,FILE=pfadn(1:pfadi)//'res_'//trim(adjustl(subarea))//'_hydraul.out',STATUS='old',  &
			POSITION='append')
	    DO d=1,dayyear
	      DO ih=1,nt
		    hour=ih
            step=(d-1)*nt+hour
	        DO j=1,nbrsec(i)
			  WRITE(11,'(5I6,4F15.3,E15.3E2,3F15.6)')id_subbas_extern(i),t,d,hour,id_sec_extern(j,i),  &
				daydepth_sec(step,j,res_index(i)),daywatelev_sec(step,j,res_index(i)),  &
				dayarea_sec(step,j,res_index(i)),daytopwidth_sec(step,j,res_index(i)),  &
				dayenergslope_sec(step,j,res_index(i)),dayhydrad_sec(step,j,res_index(i)),  &
				daymeanvel_sec(step,j,res_index(i)),daydischarge_sec(step,j,res_index(i))
			ENDDO
		  ENDDO
		ENDDO
        CLOSE(11)
		ENDIF
        DO j=1,nbrsec(i)
	      npt=npoints(j,i)
		  WRITE(section,*)id_sec_extern(j,i)
          IF (f_res_bedchange) THEN
  		  OPEN(11,FILE=pfadn(1:pfadi)//'res_'//trim(adjustl(subarea))//'_sec'//trim(adjustl(section))// &
			'_bedchange.out',STATUS='old',POSITION='append')
		  write(fmtstr,'(a,i0,a)')'(6I6,',npt,'F15.6)'		!generate format string
		  DO d=1,dayyear
	        DO ih=1,nt
		      hour=ih
              step=(d-1)*nt+hour
			  WRITE(11,fmtstr)id_subbas_extern(i),id_sec_extern(j,i),t,d,hour,npt,  &
				(dayy_sec(step,m,j,res_index(i)),m=1,npoints(j,i))
		    ENDDO
		  ENDDO
          CLOSE(11)
		  ENDIF
	    ENDDO
        IF (f_res_longitudunal) THEN
		OPEN(11,FILE=pfadn(1:pfadi)//'res_'//trim(adjustl(subarea))//'_longitudunal.out', &
			STATUS='old',POSITION='append')
		write(fmtstr,'(a,i0,a)')'(5I6,',nbrsec(i),'F15.6)'		!generate format string
		DO d=1,dayyear
	      DO ih=1,nt
		    hour=ih
            step=(d-1)*nt+hour
			WRITE(11,fmtstr)id_subbas_extern(i),t,d,hour,nbrsec1, &
				(dayminelev_sec(step,j,res_index(i)),j=1,nbrsec1)
		  ENDDO
		ENDDO
        CLOSE(11)
		ENDIF
        IF (f_res_sedbal) THEN
        OPEN(11,FILE=pfadn(1:pfadi)//'res_'//trim(adjustl(subarea))//'_sedbal.out',STATUS='old',POSITION='append')
	    DO d=1,dayyear
	      DO ih=1,nt
		    hour=ih
            step=(d-1)*nt+hour
			WRITE(11,'(4I6,4F15.3)')id_subbas_extern(i),t,d,hour,sed_inflow(step,i),sed_outflow(step,i), &
				sedimentation(step,i),daycumsed(step,res_index(i))
		  ENDDO
		ENDDO
        CLOSE(11)
		ENDIF
        IF (f_res_sedcomposition) THEN
        OPEN(11,FILE=pfadn(1:pfadi)//'res_'//trim(adjustl(subarea))//'_sedcomposition.out',STATUS='old',POSITION='append')
		write(fmtstr,'(a,i0,a)')'(5I6,',n_sed_class,'F15.6)'		!generate format string
		DO d=1,dayyear
	      DO ih=1,nt
		    hour=ih
            step=(d-1)*nt+hour
			WRITE(11,fmtstr)id_subbas_extern(i),t,d,hour,n_sed_class, &
				(dayfrsediment_out(step,res_index(i),g),g=1,n_sed_class)
!			WRITE(11,'(5I6,<n_sed_class>F15.3)')id_subbas_extern(i),t,d,hour,n_sed_class, &
!				(dayfrsediment_out(step,i,g),g=1,n_sed_class)
		  ENDDO
		ENDDO
        CLOSE(11)
		ENDIF
	  ELSE
        IF (f_res_sedbal) THEN
        OPEN(11,FILE=pfadn(1:pfadi)//'res_'//trim(adjustl(subarea))//'_sedbal.out',STATUS='old',POSITION='append')
	    DO d=1,dayyear
	      DO ih=1,nt
		    hour=ih
            step=(d-1)*nt+hour
			WRITE(11,'(4I6,4F15.3)')id_subbas_extern(i),t,d,hour,sed_inflow(step,i),sed_outflow(step,i), &
				sedimentation(step,i),daycumsed(step,res_index(i))
		  ENDDO
		ENDDO
        CLOSE(11)
		ENDIF
        IF (f_res_sedcomposition) THEN
        OPEN(11,FILE=pfadn(1:pfadi)//'res_'//trim(adjustl(subarea))//'_sedcomposition.out',STATUS='old',POSITION='append')
	    write(fmtstr,'(a,i0,a)')'(5I6,',n_sed_class,'F15.6)'		!generate format string
		DO d=1,dayyear
	      DO ih=1,nt
		    hour=ih
            step=(d-1)*nt+hour
			WRITE(11,fmtstr)id_subbas_extern(i),t,d,hour,n_sed_class, &
				(dayfrsediment_out(step,res_index(i),g),g=1,n_sed_class)
!			WRITE(11,'(5I6,<n_sed_class>F15.3)')id_subbas_extern(i),t,d,hour,n_sed_class, &
!				(dayfrsediment_out(step,i,g),g=1,n_sed_class)
		  ENDDO
		ENDDO
        CLOSE(11)
		ENDIF
	  ENDIF
	 ENDIF
    ENDDO
  ENDIF
END IF



RETURN
END SUBROUTINE semres
