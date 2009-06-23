module model_state_io
!contains subroutines for saving and laoding model state (soil water, ground water, interception)

!Till: computationally irrelevant: added missing CLOSE for ic_conds_file (caused error in linux); variable renaming
!2009-06-22

!Till: save summary on storage contents on start of model in storage.stats_start
!2009-03-27

!Till: fixed more init problems with dummy basins (prespecified outflow)
!2008-10-09

!Till: fixed problems with dummy basins (prespecified outflow)
!2008-08-27

!Till: fixed faulty initialisation whith default values that could lead to saturation >1
! check for invalidly high soil moisture contents and correct these, if necessary
!2008-07-28

!Till: storage values were summed up incorrectly in storage.stats
!2008-06-26

!Till: when initial condition files are missing, throw warning but continue with default values 
!2008-05-30

!Till: fixed bugs in saving and loading concerning confusion of tc_instance and id_tc_type 
!2008-02-11

!Till: created 
!2008-02-01

use common_h

contains
subroutine init_model_state		!load initial conditions
	if (.TRUE.) then
		call init_soil_conds(trim(pfadn)//'soil_moisture.stat',default_rel_sat)	!Till: load initial status of soil moisture
		call init_gw_conds(trim(pfadn)//'gw_storage.stat',default_gw_storage)	!Till: load initial status of gw storage
	else
		call init_soil_conds('',default_rel_sat)	!Till: load initial status of soil moisture
		call init_gw_conds('',default_gw_storage)	!Till: load initial status of gw storage
	end if
	CALL save_all_conds('','','',trim(pfadn)//'storage.stats_start')		!Till: save only summary on initial storage
end subroutine init_model_state

subroutine save_model_state		!save model state variables
	call save_all_conds(trim(pfadn)//'soil_moisture.stat',trim(pfadn)//'gw_storage.stat',trim(pfadn)//'intercept_storage.stat',trim(pfadn)//'storage.stats')	!Till: save status 
end subroutine save_model_state

		  

subroutine save_all_conds(soil_conds_file, gw_conds_file, ic_conds_file, summary_file)
!store current conditions of soil moisture, ground water and interception in the specified files
	use hymo_h
	use params_h
	use erosion_h
	use utils_h

	implicit none

	character(len=*),intent(in):: soil_conds_file, gw_conds_file, ic_conds_file,summary_file		!files to save to
	
	INTEGER :: i,j,sb_counter,lu_counter,tc_counter,svc_counter,h	! counters
	INTEGER :: i_subbas,i_lu,id_tc_type,i_svc,i_soil,i_veg		! ids of components in work
	INTEGER :: tcid_instance	!(internal) id of TC-instance (unique subbas-LU-TC-combination)
	REAL	:: total_storage_soil, total_storage_gw, total_storage_intercept 	!total amount of water stored  [m3]
	REAL	:: lu_area, svc_area	!area of current lu/svc [m3]
	INTEGER	::	soil_file_hdle, gw_file_hdle, intercept_file_hdle	!file handles to output files

	total_storage_soil=0
	total_storage_gw=0
	total_storage_intercept=0 
	
	if (trim(soil_conds_file)=='') then		!don't do anything if an empty filename is specified
		soil_file_hdle=0
	else 
		soil_file_hdle=11
		OPEN(soil_file_hdle,FILE=soil_conds_file, STATUS='replace')
		WRITE(soil_file_hdle,'(a)') 'soil moisture status (for analysis or model re-start)'
		WRITE(soil_file_hdle,*)'Subbasin', char(9),'LU', char(9),'TC' , char(9),'SVC' , char(9),'horizon', char(9),'watercontent_[mm]', char(9),'area_[m²]'		!tab separated output
	end if

	if (trim(gw_conds_file)=='') then		!don't do anything if an empty filename is specified
		gw_file_hdle=0
	else
		gw_file_hdle=12
		OPEN(gw_file_hdle,FILE=gw_conds_file, STATUS='replace')
		WRITE(gw_file_hdle,'(a)') 'ground water storage (for analysis or model re-start)'
		WRITE(gw_file_hdle,*)'Subbasin', char(9),'LU', char(9),'volume_[mm]', char(9),'area_[m²]'		!tab separated output
	end if

	if (trim(ic_conds_file)=='') then		!don't do anything if an empty filename is specified
		intercept_file_hdle=0
	else
		intercept_file_hdle=13
		OPEN(intercept_file_hdle,FILE=ic_conds_file, STATUS='replace')
		WRITE(intercept_file_hdle,'(a)') 'interception storage (for analysis or model re-start)'
		WRITE(intercept_file_hdle,*)'Subbasin', char(9),'LU', char(9),'TC' , char(9),'SVC' , char(9),'storage_[mm]', char(9),'area_[m²]'		!tab separated output
	end if

	


	DO sb_counter=1,subasin
		DO lu_counter=1,nbr_lu(sb_counter)
			i_lu=id_lu_intern(lu_counter,sb_counter)
			lu_area=area(sb_counter)*frac_lu(lu_counter,sb_counter)*1e6
			if (gw_file_hdle/=0) then
				WRITE(gw_file_hdle,'(2(I0,A1),F8.2,A1,F12.1)') id_subbas_extern(sb_counter), char(9),id_lu_extern(i_lu), char(9),&
					deepgw(sb_counter,lu_counter)/lu_area*1e3, char(9),area(sb_counter)*frac_lu(lu_counter,sb_counter)*1e6	!tab separated output
			end if
			total_storage_gw=total_storage_gw+deepgw(sb_counter,lu_counter) !sum up total storage

			DO tc_counter=1,nbrterrain(i_lu)
				tcid_instance=tcallid(sb_counter,lu_counter,tc_counter)	!id of TC instance
				if (tcid_instance==-1) cycle							!this may happen if this is merely a dummy basin with prespecified outflow
				id_tc_type=id_terrain_intern(tc_counter,i_lu)			!id of TC type
				
				DO svc_counter=1,nbr_svc(tcid_instance)
					i_soil=id_soil_intern(svc_counter,tcid_instance)		!internal id of soil type
					i_veg=  id_veg_intern(svc_counter,tcid_instance)
					i_svc=which1(svc_soil_veg(:,1)==i_soil .AND. svc_soil_veg(:,2)==i_veg)			!get internal id of svc-type

					svc_area=lu_area*fracterrain(id_terrain_intern(tc_counter,i_lu))*frac_svc(svc_counter,tcid_instance)
						
					if (intercept_file_hdle/=0) then
							! ##
							WRITE(intercept_file_hdle,'(4(I0,A1),F8.2,A1,F12.1)') id_subbas_extern(sb_counter), char(9),id_lu_extern(i_lu),char(9),&
								id_terrain_extern(id_tc_type), char(9),id_svc_extern(i_svc), char(9),intercept(tcid_instance,svc_counter),&
								char(9),	svc_area	!tab separated output
					end if
					total_storage_intercept=total_storage_intercept+intercept(tcid_instance,svc_counter)*1e-3*svc_area !sum up total storage
					

					DO h=1,nbrhori(i_soil)
						if (soil_file_hdle/=0) then
							WRITE(soil_file_hdle,'(5(I0,A1),F8.2,A1,F12.1)') id_subbas_extern(sb_counter),&
								char(9),id_lu_extern(i_lu), char(9),id_terrain_extern(id_tc_type),&
								char(9),id_svc_extern(i_svc), char(9),h, char(9),horithact(tcid_instance,svc_counter,h), &
								char(9),	svc_area	!tab separated output
						end if
						total_storage_soil=total_storage_soil+horithact(tcid_instance,svc_counter,h)*1e-3*svc_area !sum up total storage
					END DO	!loop horizons
				END DO	!loop SVC
			END DO	!loop TCs
		END DO	!loop LUs
	END DO	!loop subbasins
	CLOSE(soil_file_hdle, iostat=i_lu)	!close output files
	CLOSE(gw_file_hdle, iostat=i_lu)
	CLOSE(intercept_file_hdle, iostat=i_lu)

	
	OPEN(11,FILE=summary_file, STATUS='replace')		!write to summary file
		WRITE(11,*)'total water storage in catchment after model run [m3]'
		WRITE(11,*)'soil_storage', char(9),total_storage_soil
		WRITE(11,*)'gw_storage', char(9),total_storage_gw
		WRITE(11,*)'interception_storage', char(9),total_storage_intercept
	CLOSE(11)
end subroutine save_all_conds
		  




subroutine init_soil_conds(soil_conds_file,default_rel_sat)

!load soil moistures information from file soil_conds_file
	use hymo_h
	use params_h
	use utils_h
	use erosion_h
	implicit none

	character(len=*),intent(in):: soil_conds_file		!file to load from
	real :: default_rel_sat								!default relative saturation (to be assumed for all non-specified horizons) [0: theta_r; 1:theta_s]
	INTEGER :: i,line,errors,sb_counter,lu_counter,tc_counter,svc_counter,h	! counters
	INTEGER :: i_subbasx,i_lux,i_tcx,i_svcx, i_soilx		! external ids of components in work
	INTEGER :: i_subbas,i_lu,i_tc,i_svc, i_soil, i_veg,id_tc_type		! internal ids of components in work
	INTEGER :: lu_instance,tcid_instance,soil_instance	!(internal) id of LU,TC,soil-instance (unique subbas-LU-TC-soil-combination)
	REAL	:: horithact_temp, x
	INTEGER	:: file_read=0
	character(len=160) :: error_msg=''


	i=0
	OPEN(11,FILE=soil_conds_file,STATUS='old',action='read',  IOSTAT=i)	!check existence of file
	if (i/=0) then
		write(*,'(a,a,a)')'WARNING: Soil moisture file ''',trim(soil_conds_file),''' not found, using defaults.'
		call pause1
	else
		CLOSE(11)
	end if
	
	
	horithact=-9999.					!mark all horizons as "not (yet) initialised"
	if (trim(soil_conds_file)/='' .AND. i==0) then		!load values from file
		write(*,'(a,a,a)')'Inititalize soil moisture from file ''',trim(soil_conds_file),''''
		
		OPEN(11,FILE=soil_conds_file,STATUS='old',action='read')
		!OPEN(11,FILE=soil_conds_file,STATUS='old',action='read',readonly,shared)
		READ(11,*); READ (11,*)	!skip header lines
		line=2
		errors=0

		do while (.TRUE.)		!read whole file
			IF (len(trim(error_msg))/=0) THEN	!print error message, if occured
				if (errors==0) then !print heading at before first error
					write(*,'(A,/,6a12)')' Entities not found in current domain (ignored):','Line','subbasin','LU','TC','SVC','horizon'
				end if
				write(*,*)trim(error_msg)
				error_msg='' !reset error message
				errors=errors+1
			END IF
			
			READ(11,*,  IOSTAT=i) i_subbasx,i_lux,i_tcx,i_svcx,h,horithact_temp, x
			IF (i/=0) THEN	!no further line
				exit		!exit loop
			END IF
			line=line+1
			i_subbas=id_ext2int(i_subbasx,id_subbas_extern)	!convert to internal subbas id
			if (i_subbas==-1) then
				write(error_msg,'(2i12)')line,i_subbasx
				cycle	!proceed with next line
			end if
			i_lu=id_ext2int(i_lux,id_lu_extern)	!convert to internal lu id
			if (i_lu==-1) then
				write(error_msg,'(i12,1a12,i12)')line,'-',i_lux
				cycle	!proceed with next line
			end if
			i_tc=id_ext2int(i_tcx,id_terrain_extern)	!convert to internal tc-type id
			if (i_tc==-1) then
				write(error_msg,'(i12,2a12,i12)')line,'-','-',i_tcx
				cycle	!proceed with next line
			end if
			i_svc=id_ext2int(i_svcx,id_svc_extern)	!convert to internal svc-type id
			if (i_svc==-1) then
				write(error_msg,'(i12,3a12,i12)')line,'-','-','-',i_svcx
				cycle	!proceed with next line
			end if
			lu_counter=id_ext2int(i_lu,id_lu_intern(:,i_subbas))	!convert to position/index of lu instance in current subbasin
			if (lu_counter==-1) then
				write(error_msg,'(3i12)')line,i_subbasx,i_lux
				cycle	!proceed with next line
			end if

			tc_counter=id_ext2int(i_tc,id_terrain_intern(:,i_lu))	!convert to position/index of tc instance in current lu
			if (tc_counter==-1) then
				write(error_msg,'(4i12)')line,i_subbasx,i_lux,i_tcx
				cycle	!proceed with next line
			end if

			tcid_instance=tcallid(i_subbas,lu_counter,tc_counter)	!get the ID of the TC instance
			if (tcid_instance==-1) cycle							!this may happen if this is merely a dummy basin with prespecified outflow

			svc_counter=which1(id_soil_intern(:,tcid_instance)==svc_soil_veg(i_svc,1) .AND. &
			id_veg_intern(:,tcid_instance)==svc_soil_veg(i_svc,2)) !convert to position/index of svc instance in current tc

			!svc_counter=id_ext2int(i_tc,id_terrain_intern(:,i_tc))	!convert to position/index of svc instance in current tc
			if (svc_counter==-1 .OR. svc_counter==0) then
				write(error_msg,'(5i12)')line,i_subbasx,i_lux,i_tcx,i_svcx
				cycle	!proceed with next line
			end if


			if (horithact_temp> thetas(id_soil_intern(svc_counter,tcid_instance),h)* horiz_thickness(tcid_instance,svc_counter,h)*1.01) then			!exceeds sotrage capacity
				if (errors==0) then	!produce header before first warning only
					write(*,'(A,/,5a12)')' water content of following horizons exceed porosity. Corrected to thetaS.','subbasin','LU','TC','SVC','horizon'
				end if
				errors=errors+1
				write(*,'(5i12)')i_subbasx, i_lux, i_tcx,i_svcx, h			!issue warning
				horithact_temp = thetas(id_soil_intern(svc_counter,tcid_instance),h)* horiz_thickness(tcid_instance,svc_counter,h)
			end if

			horithact(tcid_instance,svc_counter,h)=horithact_temp			!set soil water content


		END DO
		file_read=1
		CLOSE(11)
		if (errors>0) then
			!call pause1 
		end if
	end if

	errors=0

	DO sb_counter=1,subasin			!check, if all relevant horizons have been initialized, if not, use default values
		DO lu_counter=1,nbr_lu(sb_counter)
			i_lu=id_lu_intern(lu_counter,sb_counter)
			DO tc_counter=1,nbrterrain(i_lu)
				tcid_instance=tcallid(sb_counter,lu_counter,tc_counter) !id of TC instance
				if (tcid_instance==-1) cycle							!this may happen if this is merely a dummy basin with prespecified outflow
				id_tc_type=id_terrain_intern(tc_counter,i_lu)			!id of TC type
				DO svc_counter=1,nbr_svc(tcid_instance)
					i_soil=id_soil_intern(svc_counter,tcid_instance)		!internal id of soil type
					i_veg=  id_veg_intern(svc_counter,tcid_instance)
					i_svc=which1(svc_soil_veg(:,1)==i_soil .AND. svc_soil_veg(:,2)==i_veg)			!get internal id of svc-type


					DO h=1,nbrhori(i_soil)
						if (horithact(tcid_instance,svc_counter,h)==-9999) then			!not yet set?
							if (file_read==1) then						!but this should have been done before
								if (errors==0) then	!produce header before first warning only
									write(*,'(A,f4.2,a,/,5a12)')' Following entities not initalised, using defaults (rel.saturation=',default_rel_sat,'):',&
									'subbasin','LU','TC','SVC','horizon'
								end if
								errors=errors+1
								write(*,'(5i12)')id_subbas_extern(sb_counter), id_lu_extern(i_lu), id_terrain_extern(id_tc_type),id_svc_extern(i_svc), h			!issue warning
								!write(*,5i12)'Subbasin ',id_subbas_extern(sb_counter), ', LU ',id_lu_extern(i_lu), ', TC ',id_terrain_extern(id_tc_type),&
								!', SVC ',id_svc_extern(i_svc), ', horiz ',h, ' not found, setting to ',default_rel_sat,' rel. sat.'			!issue warning
							end if
							
							horithact(tcid_instance,svc_counter,h)=(thetar(i_soil,h)+&
							(thetas(i_soil,h)-thetar(i_soil,h))*default_rel_sat)*horiz_thickness(tcid_instance,svc_counter,h)		!set to default relative saturation
						end if
					END DO
					horithact(tcid_instance,svc_counter,(nbrhori(i_soil)+1):maxhori)=0		!set water content of irrelevant horizons values to 0
				END DO
			END DO
		END DO
	END DO


!!for debugging - remove
!DO tcid_instance=1,30
!	DO i=1,nbr_svc(tcid_instance)
!		DO h=1,nbrhori(id_soil_intern(i,tcid_instance)	)
!			IF (horithact(tcid_instance,i,h) - thetas(id_soil_intern(i,tcid_instance)	,h)* horiz_thickness(tcid_instance,i,h)> 0.1 ) THEN	!Till: if water content exceeds saturation...
!				write(*,*)"Init: Oversaturated horizon in TC/svc/horizon ",tcid_instance,i,h
!				call pause1
!			END IF
!		END DO
!	END DO
!END DO
!!for debugging - remove
	

	if (errors>0) then
		call pause1 
	end if	

end subroutine init_soil_conds




subroutine init_gw_conds(gw_conds_file,default_gw_storage)

!load gw conditions from file gw_conds_file
	use hymo_h
	use params_h
	use utils_h
	use erosion_h
	implicit none

	character(len=*),intent(in):: gw_conds_file		!file to load from
	real :: default_gw_storage								!default ground water storage (to be assumed for all non-specified LUs) [mm]
	INTEGER :: i,line,errors,lu_counter	! counters
	INTEGER :: i_subbasx,i_lux		! external ids of components in work
	INTEGER :: i_subbas,i_lu		! internal ids of components in work
	INTEGER :: lu_instance	!(internal) id of LU
	REAL	:: gwvol_temp, x
	INTEGER	:: file_read=0
	character(len=160) :: error_msg=''
	REAL	:: lu_area	!area of current lu [m3]

	
	i=0
	OPEN(11,FILE=gw_conds_file,STATUS='old',action='read',  IOSTAT=i)	!check existence of file
	if (i/=0) then
		write(*,'(a,a,a)')'WARNING: GW storage file ''',trim(gw_conds_file),''' not found, using defaults.'
		call pause1
	else
		CLOSE(11)
	end if
	
	
	deepgw=-9999.					!mark all gw storages as "not (yet) read"
	if (trim(gw_conds_file)/='' .AND. i==0) then		!load values from file
		write(*,'(a,a,a)')'Inititalize ground water storage from file ''',trim(gw_conds_file),''''
		
		
		OPEN(11,FILE=gw_conds_file,STATUS='old',action='read')
		!OPEN(11,FILE=gw_conds_file,STATUS='old',action='read',readonly,shared)	
		READ(11,*); READ (11,*)	!skip header lines
		line=2
		errors=0

		do while (.TRUE.)		!read whole file
			IF (len(trim(error_msg))/=0) THEN	!print error message, if occured
				if (errors==0) then
					write(*,'(A,/,3a12)')' Entities not found in current domain (ignored):','Line','subbasin','LU'
				end if
				write(*,*)trim(error_msg)
				error_msg=''
				errors=errors+1
			END IF
			

			READ(11,*,  IOSTAT=i) i_subbasx,i_lux,gwvol_temp, x
			IF (i/=0) THEN	!no further line
				exit		!exit loop
			END IF
			line=line+1
			i_subbas=id_ext2int(i_subbasx,id_subbas_extern)	!convert to internal subbas id
			if (i_subbas==-1) then
				write(error_msg,'(2i12)')line,i_subbasx
				cycle	!proceed with next line
			end if
			i_lu=id_ext2int(i_lux,id_lu_extern)	!convert to internal lu id
			if (i_lu==-1) then
				write(error_msg,'(i12,1a12,i12)')line,'-',i_lux
				cycle	!proceed with next line
			end if


			lu_counter=id_ext2int(i_lu,id_lu_intern(:,i_subbas))	!convert to position/index of lu instance in current subbasin
			if (lu_counter==-1) then
				write(error_msg,'(3i12)')line,i_subbasx,i_lux
				cycle	!proceed with next line
			end if

			lu_area=area(i_subbas)*frac_lu(lu_counter,i_subbas)*1e6
			
			deepgw(i_subbas,lu_counter)=gwvol_temp/1000*lu_area			!set gw water content [m3]
		END DO
		file_read=1
		CLOSE(11)
		if (errors>0) then
			call pause1 
		end if
	end if

	errors=0
	DO i_subbas=1,subasin			!check, if all relevant LUs have been initialized, if not, use default values
		DO lu_counter=1,nbr_lu(i_subbas)
			i_lu=id_lu_intern(lu_counter,i_subbas)
			lu_area=area(i_subbas)*frac_lu(lu_counter,i_subbas)*1e6

			if (deepgw(i_subbas,lu_counter)==-9999) then			!not yet set?
				if (file_read==1) then						!but this should have been done before
					if (errors==0) then	!produce header before first warning only
						write(*,'(A,f5.1,a,/,2a12)')' Following entities not initalised, using defaults (storage[mm]=',&
						default_gw_storage,'):','subbasin','LU'
					end if
					errors=errors+1
					write(*,'(5i12)')id_subbas_extern(i_subbas), id_lu_extern(i_lu)			!issue warning
				end if
				deepgw(i_subbas,lu_counter)=default_gw_storage/1000*lu_area		!set to default gw storage
			end if
			
		END DO
		deepgw(i_subbas,(nbr_lu(i_subbas)+1):maxsoter)=0		!set water content of irrelevant storages to 0
	END DO
	if (errors>0) then
		call pause1 
	end if	

end subroutine init_gw_conds



end module model_state_io