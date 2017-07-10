module model_state_io
    !contains subroutines for saving and loading model state (soil water, ground water, interception)#

    !Till: dirty temporary fix to prevent crash at startup
    ! code formatting fro increased compatibility (long lines, indent, removed tabs)
    !2013-10-02

    !Jose Miguel: replaced "call pause1" for "return", which takes the execution of the subprogram to the program of an upper level
    !2013-01-14

    !Jose Miguel: added the initialization of river_storage from the river_storage.stat file.
    !2013-01-14
 
    !Jose Miguel: added code to save additional state variable r_storage, which represents the volume stored in the river chanel at the end of the model run. It is stored in river_storage.stat
    !2012-12-05

    !Jose Miguel: added code to save additional state variable lakewater_hrr, which represents the volume stored in the lake classes at the end of the model run. It is stored in lake_volume.stat
    !2012-12-05

    !Till: computationally irrelevant: minor changes in output details, improved error handling
    !2012-06-14

    !Till: computationally irrelevant: minor changes to improve compiler compatibility
    !2011-04-29

    !Till: computationally irrelevant: removed subroutine arguments for init_soil_conds and init_gw_conds that were global anyway (created confusion with gfortran 4.4)
    !2011-03-23

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
    subroutine init_model_state        !load initial conditions
            call init_soil_conds(trim(pfadn)//'soil_moisture.stat')    !Till: load initial status of soil moisture
            call init_gw_conds(trim(pfadn)//'gw_storage.stat')    !Till: load initial status of gw storage
            call init_lake_conds(trim(pfadn)//'lake_volume.stat')    !Jose Miguel: load initial status of lake storage
            call init_river_conds(trim(pfadn)//'river_storage.stat')    !Jose Miguel: load initial status of river storage
            if (dosediment) then
               call init_sediment_conds(trim(pfadn)//'sediment_storage.stat')    !Jose Miguel: load initial status of sediment storage
               call init_susp_sediment_conds(trim(pfadn)//'susp_sediment_storage.stat')    !Jose Miguel: load initial status of sediment storage
            endif
    end subroutine init_model_state

    subroutine save_model_state(backup_files, start)        !save all model state variables, optionally backup older files
    implicit none
    logical, intent(in) :: backup_files
    logical, intent(in) :: start
    
        if (.not. dosavestate) return    
        if (backup_files) then
            !keep files with initial conditions, save only summary on initial storages
            call rename(trim(pfadn)//'soil_moisture.stat'     , trim(pfadn)//'soil_moisture.stat_start')
            call rename(trim(pfadn)//'gw_storage.stat'        ,trim(pfadn)//'gw_storage.stat_start')
            call rename(trim(pfadn)//'intercept_storage.stat' ,trim(pfadn)//'intercept_storage.stat_start')
            call rename(trim(pfadn)//'lake_volume.stat'       ,trim(pfadn)//'lake_volume.stat_start')
            call rename(trim(pfadn)//'river_storage.stat'     ,trim(pfadn)//'river_storage.stat_start')
            if (dosediment) then
                call rename(trim(pfadn)//'sediment_storage.stat'     ,trim(pfadn)//'sediment_storage.stat_start')
                call rename(trim(pfadn)//'susp_sediment_storage.stat',trim(pfadn)//'susp_sediment_storage.stat_start')
            end if    
            CALL save_all_conds('','','','','','','',trim(pfadn)//'storage.stats', start)        !Till: save only summary on initial storage
        else    
            if (dosediment) then
                call save_all_conds(trim(pfadn)//'soil_moisture.stat',trim(pfadn)//'gw_storage.stat',trim(pfadn)//'intercept_storage.stat',&
                trim(pfadn)//'lake_volume.stat',trim(pfadn)//'river_storage.stat',&
                trim(pfadn)//'sediment_storage.stat',trim(pfadn)//'susp_sediment_storage.stat',trim(pfadn)//'storage.stats', start)    !Till: save status
            else
               call save_all_conds(trim(pfadn)//'soil_moisture.stat',trim(pfadn)//'gw_storage.stat',trim(pfadn)//'intercept_storage.stat',&
                trim(pfadn)//'lake_volume.stat',trim(pfadn)//'river_storage.stat',&
                '','',trim(pfadn)//'storage.stats', start)
            endif
        end if    
    end subroutine save_model_state



    subroutine save_all_conds(soil_conds_file, gw_conds_file, ic_conds_file, lake_conds_file, river_conds_file, sediment_conds_file, susp_sediment_conds_file, summary_file, start)
        !store current conditions of soil moisture, ground water and interception in the specified files
        use hymo_h
        use params_h
        use erosion_h
        use utils_h
        !Jose Miguel: the following modules are also necessary, because we are now saving state variables that are declared inside these modules
        use routing_h
        use lake_h
        use reservoir_h
        use time_h
        use common_h

        implicit none

        character(len=*),intent(in):: soil_conds_file, gw_conds_file, ic_conds_file,lake_conds_file,river_conds_file,sediment_conds_file,susp_sediment_conds_file,summary_file        !files to save to
        logical,intent(in), optional:: start

        INTEGER :: sb_counter,lu_counter,tc_counter,svc_counter,h,acud_class, k, tt, digits    ! counters
        INTEGER :: i_lu,id_tc_type,i_svc,i_soil,i_veg        ! ids of components in work
        INTEGER :: tcid_instance    !(internal) id of TC-instance (unique subbas-LU-TC-combination)
        REAL    :: total_storage_soil, total_storage_gw, total_storage_intercept, total_storage_lake(5), total_storage_river !, total_storage_sediment,total_storage_suspsediment    !total amount of water stored  [m3]
        REAL    :: lu_area, svc_area    !area of current lu/svc [m3]
        INTEGER    ::    soil_file_hdle, gw_file_hdle, intercept_file_hdle, lake_file_hdle, river_file_hdle,sediment_file_hdle,susp_sediment_file_hdle    !file handles to output files
		character(len=1000) :: fmtstr    !string for formatting file output
        character(len=6) :: suffix

        suffix=""
        if (present(start)) then
            if (start) suffix="_start"
        end if
        total_storage_soil=0.
        total_storage_gw=0.
        total_storage_intercept=0.
        total_storage_lake=0.
        total_storage_river=0.
!        total_storage_sediment=0.

!write file headers
        if (trim(soil_conds_file)=='') then        !don't do anything if an empty filename is specified
            soil_file_hdle=0
        else
            soil_file_hdle=11
            OPEN(soil_file_hdle,FILE=soil_conds_file//trim(suffix), STATUS='replace')
            WRITE(soil_file_hdle,'(a)') 'soil moisture status (for analysis or model re-start)'
            WRITE(soil_file_hdle,*)'Subbasin', char(9),'LU', char(9),'TC' , char(9),'SVC' , char(9),'horizon', char(9),&
                'watercontent_[mm]', char(9),'area_[m²]'        !tab separated output
        end if

        if (trim(gw_conds_file)=='') then        !don't do anything if an empty filename is specified
            gw_file_hdle=0
        else
            gw_file_hdle=12
            OPEN(gw_file_hdle,FILE=gw_conds_file//trim(suffix), STATUS='replace')
            WRITE(gw_file_hdle,'(a)') 'ground water storage (for analysis or model re-start)'
            WRITE(gw_file_hdle,*)'Subbasin', char(9),'LU', char(9),'volume_[mm]', char(9),'area_[m²]'        !tab separated output
        end if

        if (trim(ic_conds_file)=='') then        !don't do anything if an empty filename is specified
            intercept_file_hdle=0
        else
            intercept_file_hdle=13
            OPEN(intercept_file_hdle,FILE=ic_conds_file//trim(suffix), STATUS='replace')
            WRITE(intercept_file_hdle,'(a)') 'interception storage (for analysis or model re-start)'
            WRITE(intercept_file_hdle,*)'Subbasin', char(9),'LU', char(9),'TC' , char(9),'SVC' , char(9),&
                'storage_[mm]', char(9),'area_[m²]'        !tab separated output
        end if

        if (trim(river_conds_file)=='' .OR. (river_transport == 1)) then        !don't do anything if an empty filename is specified
            river_file_hdle=0
        else
            river_file_hdle=14
            OPEN(river_file_hdle,FILE=river_conds_file//trim(suffix), STATUS='replace')
            WRITE(river_file_hdle,'(a)') 'river reach volume status (for analysis or model re-start)'
            WRITE(river_file_hdle,*)'Subbasin', char(9),'volume[m^3]' !tab separated output
        endif

        if (trim(lake_conds_file)=='') then        !don't do anything if an empty filename is specified
            lake_file_hdle=0
        else
            lake_file_hdle=15
            OPEN(lake_file_hdle,FILE=lake_conds_file//trim(suffix), STATUS='replace')
            if (doacud) then
                WRITE(lake_file_hdle,'(a)') 'Lake volume status (for analysis or model re-start)'
                WRITE(lake_file_hdle,*)'Subbasin', char(9),'reservoir_size_class',char(9),'volume[m^3]' !tab separated output
            else
                CLOSE(lake_file_hdle,status='delete')
            end if
        end if


        if (trim(sediment_conds_file)=='') then        !don't do anything if an empty filename is specified
            sediment_file_hdle=0
        else
            sediment_file_hdle=16
            OPEN(sediment_file_hdle,FILE=sediment_conds_file//trim(suffix), STATUS='replace')
            WRITE(sediment_file_hdle,'(a)') 'river reach deposition sediment weight status (for analysis or model re-start)'
            WRITE(sediment_file_hdle,'(A)')'Subbasin'//char(9)//'particle_size_class'//char(9)//'mass[t]' !tab separated output
        endif
        
        if (trim(susp_sediment_conds_file)=='') then        !don't do anything if an empty filename is specified
            susp_sediment_file_hdle=0
        else
            susp_sediment_file_hdle=17
            OPEN(susp_sediment_file_hdle,FILE=susp_sediment_conds_file//trim(suffix), STATUS='replace')
            WRITE(susp_sediment_file_hdle,'(a)') 'river reach suspended sediment weight status (for analysis or model re-start)'
            WRITE(susp_sediment_file_hdle,'(A)')'Subbasin'//char(9)//'particle_size_class'//char(9)//'mass[t]' !tab separated output
        endif
        
!write file contents        
        !Jose Miguel: write river storage state to .stat file .
        digits=floor(log10(max(1.0,maxval(r_storage))))+1    !Till: number of pre-decimal digits required
        if (digits<10) then
            write(fmtstr,'(a,i0,a,i0,a)') '(I0,A1,F',min(11,digits+4),'.',min(3,11-digits-1),'))'        !generate format string
        else       
            fmtstr='(I0,A1,E12.5)' !for large numbers, use exponential notation
        end if
        DO sb_counter=1,subasin !we wrap subbasin loop around each single entity to save time in generating format string
            if (river_file_hdle/=0) then
                WRITE(river_file_hdle,trim(fmtstr))id_subbas_extern(sb_counter), char(9),r_storage(sb_counter) !tab separated output
            endif
            total_storage_river=total_storage_river+r_storage(sb_counter) !sum up total storage
        END DO
        
        if (dosediment) then
            !write riverbed sediment storage
            !generate format string
            if (sediment_file_hdle/=0) then
                digits=floor(log10(max(1.0,maxval(riverbed_storage))))+1    !Till: number of pre-decimal digits required
                if (digits<10) then
                    write(fmtstr,'(a,i0,a,i0,a)') '(I0,A1,I0,A1,F',min(11,digits+4),'.',min(3,11-digits-1),'))'        !generate format string
                else       
                    fmtstr='(I0,A1,I0,A1,E12.5)' !for large numbers, use exponential notation
                end if
            END IF    
            DO sb_counter=1,subasin !we wrap subbasin loop around each single entity to save time in generating format string
                if (sediment_file_hdle/=0) then
                    do k=1, n_sed_class
                            WRITE(sediment_file_hdle,fmtstr)id_subbas_extern(sb_counter), char(9), k, char(9) ,riverbed_storage(sb_counter,k) !print each sediment class
                    enddo
                endif
                !total_storage_sediment=total_storage_sediment+sum(riverbed_storage(sb_counter,:)) !sum up total storage
            END DO
        
            !write suspended sediment storage
            if (susp_sediment_file_hdle/=0) then
                digits=floor(log10(max(1.0,maxval(sed_storage))))+1    !Till: number of pre-decimal digits required
                if (digits<10) then
                  write(fmtstr,'(a,i0,a,i0,a)') '(I0,A1,I0,A1,F',min(11,digits+4),'.',min(3,11-digits-1),'))'        !generate format string
                else       
                   fmtstr='(I0,A1,I0,A1,E12.5)' !for large numbers, use exponential notation
                end if
            end if
            DO sb_counter=1,subasin !we wrap subbasin loop around each single entity to save time in generating format string
                if (susp_sediment_file_hdle/=0) then
                    do k=1, n_sed_class
                         WRITE(susp_sediment_file_hdle,fmtstr)id_subbas_extern(sb_counter), char(9), k, char(9), sed_storage(sb_counter,k) !print each sediment class
                    enddo
                endif
                !total_storage_suspsediment=total_storage_suspsediment+sum(sed_storage(sb_counter,:)) !sum up total storage            
            END DO    
        end if !dosediment
                
        
        IF (doacud) THEN
            tt = (d-2)*nt+hour !index to last valid value
            if (tt<1) tt=1 !Till: dirty fix to prevent crash at start up. José, please check this
            DO sb_counter=1,subasin
                DO acud_class=1,5
				    if (lake_file_hdle/=0) then
					    WRITE(lake_file_hdle,'(I0,A1,I0,A1,F8.2)') id_subbas_extern(sb_counter), char(9),acud_class,char(9),&
						    lakewater_hrr(tt,sb_counter,acud_class)
				    endif
				    total_storage_lake(acud_class) = total_storage_lake(acud_class)+lakewater_hrr(tt,sb_counter,acud_class) * acud(sb_counter,acud_class) !sum up total storage
                ENDDO
            END DO
        END IF !small reservoirs
       
        DO sb_counter=1,subasin
            DO lu_counter=1,nbr_lu(sb_counter) !groundwater, intercept and soil storages inside
                i_lu=id_lu_intern(lu_counter,sb_counter)
                lu_area=area(sb_counter)*frac_lu(lu_counter,sb_counter)*1e6
                if (gw_file_hdle/=0) then
                    WRITE(gw_file_hdle,'(2(I0,A1),F8.2,A1,F12.1)') id_subbas_extern(sb_counter), char(9),id_lu_extern(i_lu),&
                        char(9),deepgw(sb_counter,lu_counter)/lu_area*1e3,&
                        char(9),area(sb_counter)*frac_lu(lu_counter,sb_counter)*1e6    !tab separated output
                endif
                total_storage_gw=total_storage_gw+deepgw(sb_counter,lu_counter) !sum up total storage

                DO tc_counter=1,nbrterrain(i_lu) !intercept and soil storages inside
                    tcid_instance=tcallid(sb_counter,lu_counter,tc_counter)    !id of TC instance
                    if (tcid_instance==-1) cycle                            !this may happen if this is merely a dummy basin with prespecified outflow
                    id_tc_type=id_terrain_intern(tc_counter,i_lu)            !id of TC type

                    DO svc_counter=1,nbr_svc(tcid_instance) !intercept and soil storages inside
                        i_soil=id_soil_intern(svc_counter,tcid_instance)        !internal id of soil type
                        i_veg=  id_veg_intern(svc_counter,tcid_instance)
                        i_svc=which1(svc_soil_veg(:,1)==i_soil .AND. svc_soil_veg(:,2)==i_veg)            !get internal id of svc-type

                        svc_area=lu_area*fracterrain(id_terrain_intern(tc_counter,i_lu))*frac_svc(svc_counter,tcid_instance)

                        if (intercept_file_hdle/=0) then
                            ! ##
                            WRITE(intercept_file_hdle,'(4(I0,A1),F8.2,A1,F12.1)') id_subbas_extern(sb_counter), char(9),&
                                id_lu_extern(i_lu),char(9),id_terrain_extern(id_tc_type), char(9),&
                                id_svc_extern(i_svc), char(9),intercept(tcid_instance,svc_counter),&
                                char(9),    svc_area    !tab separated output
                        endif
                        total_storage_intercept=total_storage_intercept+intercept(tcid_instance,svc_counter)*1e-3*svc_area !sum up total storage


                        DO h=1,nbrhori(i_soil) !soil storages inside
                            if (soil_file_hdle/=0) then
                                WRITE(soil_file_hdle,'(5(I0,A1),F8.2,A1,F12.1)') id_subbas_extern(sb_counter),&
                                    char(9),id_lu_extern(i_lu), char(9),id_terrain_extern(id_tc_type),&
                                    char(9),id_svc_extern(i_svc), char(9),h, char(9),horithact(tcid_instance,svc_counter,h), &
                                    char(9),    svc_area    !tab separated output
                            endif
                            total_storage_soil=total_storage_soil+horithact(tcid_instance,svc_counter,h)*1e-3*svc_area !sum up total storage
                        ENDDO    !loop horizons
                    ENDDO    !loop SVC
                ENDDO    !loop TCs
            ENDDO    !loop LUs
        ENDDO    !loop subbasins
        
    
        CLOSE(soil_file_hdle, iostat=i_lu)    !close output files
        CLOSE(gw_file_hdle, iostat=i_lu)
        CLOSE(intercept_file_hdle, iostat=i_lu)
        CLOSE(lake_file_hdle, iostat=i_lu)
        CLOSE(river_file_hdle, iostat=i_lu)    !close output files
        CLOSE(sediment_file_hdle, iostat=i_lu)    !close output files
        CLOSE(susp_sediment_file_hdle, iostat=i_lu)    !close output files
        

        OPEN(11,FILE=summary_file//trim(suffix), STATUS='replace')        !write to summary file
        WRITE(11,*)'total water storage in catchment after model run [m3]'
        WRITE(11,*)'soil_storage', char(9),total_storage_soil
        WRITE(11,*)'gw_storage', char(9),total_storage_gw
        WRITE(11,*)'interception_storage', char(9),total_storage_intercept

        DO acud_class=1,5
            WRITE(11,'(A,I0,A,F12.1)')'lake_storage', acud_class, char(9),total_storage_lake(acud_class)
        ENDDO
        WRITE(11,*)'river_storage', char(9),total_storage_river

        CLOSE(11)
    end subroutine save_all_conds





    subroutine init_soil_conds(soil_conds_file)

        !load soil moistures information from file soil_conds_file
        use hymo_h
        use params_h
        use utils_h
        use erosion_h
        implicit none

        character(len=*),intent(in):: soil_conds_file        !file to load from
        INTEGER :: i,line,errors,sb_counter,lu_counter,tc_counter,svc_counter,h    ! counters
        INTEGER :: i_subbasx,i_lux,i_tcx,i_svcx !, i_soilx        ! external ids of components in work
        INTEGER :: i_subbas,i_lu,i_tc,i_svc, i_soil, i_veg,id_tc_type        ! internal ids of components in work
        INTEGER :: tcid_instance !,lu_instance, !,soil_instance    !(internal) id of LU,TC,soil-instance (unique subbas-LU-TC-soil-combination)
        REAL    :: horithact_temp, x
        INTEGER    :: file_read=0
        character(len=160) :: linestr='', error_msg=''


        i=0
        OPEN(11,FILE=soil_conds_file,STATUS='old',action='read',  IOSTAT=i)    !check existence of file
        if (i/=0) then
            write(*,'(a,a,a)')'WARNING: Soil moisture file ''',trim(soil_conds_file),''' not found, using defaults.'
            CLOSE(11)
			return
        end if


        horithact=-9999.                    !mark all horizons as "not (yet) initialised"
        if (trim(soil_conds_file)/='' .AND. i==0) then        !load values from file
            write(*,'(a,a,a)')'Initialize soil moisture from file ''',trim(soil_conds_file),''''

            READ(11,*); READ (11,*)    !skip header lines
            line=2
            errors=0

            do while (.TRUE.)        !read whole file
                IF (len(trim(error_msg))/=0) THEN    !print error message, if occured
                    if (errors==0) then !print heading at before first error
                        write(*,'(A,/,6a12)')' Entities not found in current domain (ignored):','Line','subbasin','LU','TC','SVC',&
                            'horizon'
                    end if
                    write(*,*)trim(error_msg)
                    error_msg='' !reset error message
                    errors=errors+1
                END IF

                READ(11,'(A)',  IOSTAT=i) linestr

                IF (i==24 .OR. i==-1) THEN    !end of file
                    exit        !exit loop
                END IF
                line=line+1

                READ(linestr,*,  IOSTAT=i) i_subbasx,i_lux,i_tcx,i_svcx,h,horithact_temp, x
                IF (i/=0) THEN    !format error
                    write(error_msg,'(i12,1a12,a)')line,'-',trim(linestr)
                    cycle    !proceed with next line
                END IF

                i_subbas=id_ext2int(i_subbasx,id_subbas_extern)    !convert to internal subbas id
                if (i_subbas==-1) then
                    write(error_msg,'(2i12)')line,i_subbasx
                    cycle    !proceed with next line
                end if
                i_lu=id_ext2int(i_lux,id_lu_extern)    !convert to internal lu id
                if (i_lu==-1) then
                    write(error_msg,'(i12,1a12,i12)')line,'-',i_lux
                    cycle    !proceed with next line
                end if
                i_tc=id_ext2int(i_tcx,id_terrain_extern)    !convert to internal tc-type id
                if (i_tc==-1) then
                    write(error_msg,'(i12,2a12,i12)')line,'-','-',i_tcx
                    cycle    !proceed with next line
                end if
                i_svc=id_ext2int(i_svcx,id_svc_extern)    !convert to internal svc-type id
                if (i_svc==-1) then
                    write(error_msg,'(i12,3a12,i12)')line,'-','-','-',i_svcx
                    cycle    !proceed with next line
                end if
                lu_counter=id_ext2int(i_lu,id_lu_intern(:,i_subbas))    !convert to position/index of lu instance in current subbasin
                if (lu_counter==-1) then
                    write(error_msg,'(3i12)')line,i_subbasx,i_lux
                    cycle    !proceed with next line
                end if

                tc_counter=id_ext2int(i_tc,id_terrain_intern(:,i_lu))    !convert to position/index of tc instance in current lu
                if (tc_counter==-1) then
                    write(error_msg,'(4i12)')line,i_subbasx,i_lux,i_tcx
                    cycle    !proceed with next line
                end if

                tcid_instance=tcallid(i_subbas,lu_counter,tc_counter)    !get the ID of the TC instance
                if (tcid_instance==-1) cycle                            !this may happen if this is merely a dummy basin with prespecified outflow

                svc_counter=which1(id_soil_intern(:,tcid_instance)==svc_soil_veg(i_svc,1) .AND. &
                    id_veg_intern(:,tcid_instance)==svc_soil_veg(i_svc,2)) !convert to position/index of svc instance in current tc

                !svc_counter=id_ext2int(i_tc,id_terrain_intern(:,i_tc))    !convert to position/index of svc instance in current tc
                if (svc_counter==-1 .OR. svc_counter==0) then
                    write(error_msg,'(5i12)')line,i_subbasx,i_lux,i_tcx,i_svcx
                    cycle    !proceed with next line
                end if


                if (horithact_temp > thetas(id_soil_intern(svc_counter,tcid_instance),h)*&
                    horiz_thickness(tcid_instance,svc_counter,h)*1.01) then            !exceeds sotrage capacity
                    if (errors==0) then    !produce header before first warning only
                        write(*,'(A,/,5a12)')' water content of following horizons exceed porosity. Corrected to thetaS.',&
                            'subbasin','LU','TC','SVC','horizon'
                    end if
                    errors=errors+1
                    write(*,'(5i12)')i_subbasx, i_lux, i_tcx,i_svcx, h            !issue warning
                    horithact_temp = thetas(id_soil_intern(svc_counter,tcid_instance),h) * &
                        horiz_thickness(tcid_instance,svc_counter,h)
                end if

                horithact(tcid_instance,svc_counter,h)=horithact_temp            !set soil water content


            END DO
            file_read=1
            CLOSE(11)
!            if (errors>0) then
!                return
!            end if
        end if

        errors=0

        DO sb_counter=1,subasin            !check, if all relevant horizons have been initialized, if not, use default values
            DO lu_counter=1,nbr_lu(sb_counter)
                i_lu=id_lu_intern(lu_counter,sb_counter)
                DO tc_counter=1,nbrterrain(i_lu)
                    tcid_instance=tcallid(sb_counter,lu_counter,tc_counter) !id of TC instance
                    if (tcid_instance==-1) cycle                            !this may happen if this is merely a dummy basin with prespecified outflow
                    id_tc_type=id_terrain_intern(tc_counter,i_lu)            !id of TC type
                    DO svc_counter=1,nbr_svc(tcid_instance)
                        i_soil=id_soil_intern(svc_counter,tcid_instance)        !internal id of soil type
                        i_veg=  id_veg_intern(svc_counter,tcid_instance)
                        i_svc=which1(svc_soil_veg(:,1)==i_soil .AND. svc_soil_veg(:,2)==i_veg)            !get internal id of svc-type


                        DO h=1,nbrhori(i_soil)
                            if (horithact(tcid_instance,svc_counter,h)==-9999.) then            !not yet set?
                                if (file_read==1) then                        !but this should have been done before
                                    if (errors==0) then    !produce header before first warning only
                                        write(*,'(A,f4.2,a,/,5a12)')' Following entities not initialised, using defaults '// &
                                            '(rel.saturation=', default_rel_sat,'):','subbasin','LU','TC','SVC','horizon'
                                    end if
                                    errors=errors+1
                                    write(*,'(5i12)')id_subbas_extern(sb_counter), id_lu_extern(i_lu),&
                                        id_terrain_extern(id_tc_type), id_svc_extern(i_svc), h            !issue warning
                                   !write(*,5i12)'Subbasin ',id_subbas_extern(sb_counter), ', LU ',id_lu_extern(i_lu), ', TC ',id_terrain_extern(id_tc_type),&
                                   !', SVC ',id_svc_extern(i_svc), ', horiz ',h, ' not found, setting to ',default_rel_sat,' rel. sat.'            !issue warning
                                end if

                                horithact(tcid_instance,svc_counter,h)=(thetar(i_soil,h)+&
                                    (thetas(i_soil,h)-thetar(i_soil,h))*default_rel_sat)*horiz_thickness(tcid_instance,svc_counter,h)        !set to default relative saturation
                            end if
                        END DO
                        horithact(tcid_instance,svc_counter,(nbrhori(i_soil)+1):maxhori)=0.        !set water content of irrelevant horizons values to 0
                    END DO
                END DO
            END DO
        END DO


        if (errors>0) then
            return
        end if

    end subroutine init_soil_conds




    subroutine init_gw_conds(gw_conds_file)

        !load gw conditions from file gw_conds_file
        use hymo_h
        use params_h
        use utils_h
        use erosion_h
        implicit none

        character(len=*),intent(in):: gw_conds_file        !file to load from

        INTEGER :: i,line,errors,lu_counter    ! counters
        INTEGER :: i_subbasx,i_lux        ! external ids of components in work
        INTEGER :: i_subbas,i_lu        ! internal ids of components in work
        !    INTEGER :: lu_instance    !(internal) id of LU
        REAL    :: gwvol_temp, x
        INTEGER    :: file_read=0
        character(len=160) :: linestr='', error_msg=''
        REAL    :: lu_area    !area of current lu [m3]


        i=0
        OPEN(11,FILE=gw_conds_file,STATUS='old',action='read',  IOSTAT=i)    !check existence of file
        if (i/=0) then
            write(*,'(a,a,a)')'WARNING: GW storage file ''',trim(gw_conds_file),''' not found, using defaults.'
            CLOSE(11)
			return
        end if


        deepgw=-9999.                    !mark all gw storages as "not (yet) read"
        if (trim(gw_conds_file)/='' .AND. i==0) then        !load values from file
            write(*,'(a,a,a)')'Initialize ground water storage from file ''',trim(gw_conds_file),''''

            READ(11,*); READ (11,*)    !skip header lines
            line=2
            errors=0

            do while (.TRUE.)        !read whole file
                IF (len(trim(error_msg))/=0) THEN    !print error message, if occured
                    if (errors==0) then
                        write(*,'(A,/,3a12)')' Entities not found in current domain (ignored):','Line','subbasin','LU'
                    end if
                    write(*,*)trim(error_msg)
                    error_msg=''
                    errors=errors+1
                END IF

                READ(11,'(A)',  IOSTAT=i) linestr

                IF (i==24 .OR. i==-1) THEN    !end of file
                    exit        !exit loop
                END IF
                line=line+1

                READ(linestr,*,  IOSTAT=i) i_subbasx,i_lux,gwvol_temp, x
                IF (i/=0) THEN    !format error
                    write(error_msg,'(i12,1a12,a)')line,'-',trim(linestr)
                    cycle    !proceed with next line
                END IF

                i_subbas=id_ext2int(i_subbasx,id_subbas_extern)    !convert to internal subbas id
                if (i_subbas==-1) then
                    write(error_msg,'(2i12)')line,i_subbasx
                    cycle    !proceed with next line
                end if
                i_lu=id_ext2int(i_lux,id_lu_extern)    !convert to internal lu id
                if (i_lu==-1) then
                    write(error_msg,'(i12,1a12,i12)')line,'-',i_lux
                    cycle    !proceed with next line
                end if


                lu_counter=id_ext2int(i_lu,id_lu_intern(:,i_subbas))    !convert to position/index of lu instance in current subbasin
                if (lu_counter==-1) then
                    write(error_msg,'(3i12)')line,i_subbasx,i_lux
                    cycle    !proceed with next line
                end if

                lu_area=area(i_subbas)*frac_lu(lu_counter,i_subbas)*1e6

                deepgw(i_subbas,lu_counter)=gwvol_temp/1000.*lu_area            !set gw water content [m3]
            END DO
            file_read=1
            CLOSE(11)
            if (errors>0) then
                return
            end if
        end if

        errors=0
        DO i_subbas=1,subasin            !check, if all relevant LUs have been initialized, if not, use default values
            DO lu_counter=1,nbr_lu(i_subbas)
                i_lu=id_lu_intern(lu_counter,i_subbas)
                lu_area=area(i_subbas)*frac_lu(lu_counter,i_subbas)*1e6

                if (deepgw(i_subbas,lu_counter)==-9999.) then            !not yet set?
                    if (file_read==1) then                        !but this should have been done before
                        if (errors==0) then    !produce header before first warning only
                            write(*,'(A,f5.1,a,/,2a12)')' Following entities not initialised, using defaults (gw_storage[mm]=',&
                                default_gw_storage,'):','subbasin','LU'
                        end if
                        errors=errors+1
                        write(*,'(5i12)')id_subbas_extern(i_subbas), id_lu_extern(i_lu)            !issue warning
                    end if
                    deepgw(i_subbas,lu_counter)=default_gw_storage/1000.*lu_area        !set to default gw storage
                end if

            END DO
            deepgw(i_subbas,(nbr_lu(i_subbas)+1):maxsoter)=0.        !set water content of irrelevant storages to 0
        END DO
        if (errors>0) then
            return
        end if
        close(11)
    end subroutine init_gw_conds


    subroutine init_river_conds(river_conds_file)
        use routing_h
        use params_h
        use utils_h
    	use hymo_h

        character(len=*),intent(in):: river_conds_file        !file to load from
        integer :: subbas_id, iostatus, i
        real :: dummy1

        if (river_transport == 2) then
            OPEN(11,FILE=river_conds_file,STATUS='old',action='read', IOSTAT=i)    !check existence of file
            if (i/=0) then
                write(*,'(a,a,a)')'WARNING: River storage file ''',trim(river_conds_file),''' not found, using defaults.'
                CLOSE(11)
			    return
            end if

            write(*,'(a,a,a)')'Initialize river storage from file ''',trim(river_conds_file),'''.'
    
            !read 2 header lines into buffer
            READ(11,*); READ(11,*)
		    r_storage(:)=-1. !indicator for "not read"
        
            DO WHILE (.TRUE.) 
                READ(11,*,IOSTAT=iostatus) i, dummy1
			    IF (iostatus /=0) exit 
            
                subbas_id = id_ext2int(i, id_subbas_extern) !convert external to internal id
			    if (subbas_id < 1 .OR. subbas_id > subasin) then
				    WRITE(*,'(a,i0,a)') 'WARNING: unknown subbasin ', i,' in river_storage.stat, ignored.'
                    cycle 
                end if
            
                if (dummy1 < 0.) then
				    WRITE(*,'(a,i0,a)') 'Error: negative value for subbasin ', i,' in river_storage.stat.'
                    stop 
			    end if
       
                r_storage(subbas_id)=dummy1        !add the previous storage to the river reach additionally to potential volume from spring or runoff contribution.
            END DO
            close(11)
        
            if (count(r_storage==-1.) > 0) then  
                WRITE(*,'(A)') 'WARNING: could not read initial river storage from river_storage.stat for the following subbasins, assumed 0:'
                DO subbas_id=1,subasin
                    if (r_storage(i)==-1.) WRITE(*,'(i0)') subbas_id
                    r_storage(i)=0.
                END DO
            end if
        else
           r_storage(:)=0.  !old unit-hydrograph routing, set Muskingum storages to 0
        end if !river_transport == 2    
        
    end subroutine init_river_conds

    
   subroutine init_sediment_conds(sediment_conds_file)
        use routing_h
        use params_h
        use utils_h
	    use hymo_h
        implicit none

        character(len=*),intent(in):: sediment_conds_file        !file to load from
        integer :: subbas_id, iostatus, i, k
        real :: dummy1

        riverbed_storage(:,:)=-1. !indicator for "not read"
        OPEN(11,FILE=sediment_conds_file,STATUS='old',action='read', IOSTAT=i)    !check existence of file
        if (i/=0) then
            write(*,'(a,a,a)')'WARNING: Sediment storage file ''',trim(sediment_conds_file),''' not found, using defaults.'
            CLOSE(11)
			return
        end if

        write(*,'(a,a,a)')'Initialize sediment storage from file ''',trim(sediment_conds_file),'''.'
    
        !read 2 header lines into buffer
        READ(11,*); READ(11,*)
        
        DO WHILE (.TRUE.) 
	        READ(11,*,IOSTAT=iostatus) i, k, dummy1
            IF (iostatus == -1) exit !end of file
		    IF (iostatus /= 0) THEN
		        WRITE(*,'(a,a,a)') 'WARNING: format error in sediment_storage.stat, line skipped, assumed 0.'
            ENDIF

            subbas_id = id_ext2int(i, id_subbas_extern) !convert external to internal id
			if (subbas_id < 1 .OR. subbas_id > subasin) then
				WRITE(*,'(a,i0,a)') 'WARNING: unknown subbasin ',i,' in sediment_storage.stat, ignored.'
                cycle 
            end if
                
            if (k < 1 .OR. k > n_sed_class) then
				WRITE(*,'(a,i0,a)') 'WARNING: unknown particle size class ',k,' in sediment_storage.stat, ignored.'
                cycle 
			end if
	        riverbed_storage(subbas_id,k)=dummy1
        ENDDO
        close(11)
        
        if (count(riverbed_storage==-1.) > 0) then  
            WRITE(*,'(A)') 'WARNING: could not read initial river sediment storage from sediment_storage.stat for the following subbasins, assumed 0:'
            DO subbas_id=1,subasin
                if (count(riverbed_storage(subbas_id,:)==-1)> 0) WRITE(*,'(i0)') subbas_id
            END DO
            where(riverbed_storage==-1) riverbed_storage=0.
        end if

   end subroutine init_sediment_conds

       
   subroutine init_susp_sediment_conds(susp_sediment_conds_file)
        use routing_h
        use params_h
        use utils_h
	    use hymo_h
        implicit none

        character(len=*),intent(in):: susp_sediment_conds_file        !file to load from
        integer :: subbas_id, iostatus, i, k
        real :: dummy1

        sed_storage(:,:)=-1. !indicator for "not read"
        OPEN(11,FILE=susp_sediment_conds_file,STATUS='old',action='read', IOSTAT=i)    !check existence of file
        if (i/=0) then
            write(*,'(a,a,a)')'WARNING: Sediment storage file ''',trim(susp_sediment_conds_file),''' not found, using defaults.'
            CLOSE(11)
			return
        end if

        write(*,'(a,a,a)')'Initialize sediment storage from file ''',trim(susp_sediment_conds_file),'''.'
    
        !read 2 header lines into buffer
        READ(11,*); READ(11,*)
        
        DO WHILE (.TRUE.) 
	        READ(11,*,IOSTAT=iostatus) i, k, dummy1
            IF (iostatus == -1) exit !end of file
		    IF (iostatus /= 0) THEN
		        WRITE(*,'(a,a,a)') 'WARNING: format error in susp_sediment_storage.stat, line skipped, assumed 0.'
            ENDIF

            subbas_id = id_ext2int(i, id_subbas_extern) !convert external to internal id
			if (subbas_id < 1 .OR. subbas_id > subasin) then
				WRITE(*,'(a,i0,a)') 'WARNING: unknown subbasin ',i,' in susp_sediment_storage.stat, ignored.'
                cycle 
            end if
                
            if (k < 1 .OR. k > n_sed_class) then
				WRITE(*,'(a,i0,a)') 'WARNING: unknown particle size class ',k,' in susp_sediment_storage.stat, ignored.'
                cycle 
			end if
	        sed_storage(subbas_id,k)=dummy1
        ENDDO
        close(11)
        
        if (count(sed_storage==-1.) > 0) then  
            WRITE(*,'(A)') 'WARNING: could not read initial river sediment storage from susp_sediment_storage.stat for the following subbasins, assumed 0:'
            DO subbas_id=1,subasin
				if (count(sed_storage(subbas_id,:)==-1) > 0) WRITE(*,'(i0)') subbas_id
            END DO
            where(sed_storage==-1) sed_storage=0.
        end if

    end subroutine init_susp_sediment_conds

   
    subroutine init_lake_conds(lake_conds_file)
        !the variable that should be initialized here is the lakewater0 which is used in lake.f90 and declared in reservoir_lake_h.f90
        use lake_h
        use time_h
        use params_h
        use utils_h
	    use hymo_h
        
        implicit none

        character(len=*),intent(in):: lake_conds_file        !file to load from
        integer :: sb_counter, acud_class, iostatus, i, k, subbas_id
        real :: dummy1
    
        if (.not. doacud) then !don't try to load file if reservoirs have been disabled anyway
            return
        end if

        lakewater_hrr(1,:,:) = -1 !for detecting uninitialized values later
        
        OPEN(11,FILE=lake_conds_file,STATUS='old',action='read',  IOSTAT=i)    !check existence of file
        if (i/=0) then
            write(*,'(a,a,a)')'WARNING: Lake storage file ''',trim(lake_conds_file),''' not found, using defaults.'
            CLOSE(11)
            return
        end if

        write(*,'(a,a,a)')'Initialize lake storage from file ''',trim(lake_conds_file),'''.'
    
        READ(11,*); READ(11,*)!read 2 header lines into buffer

        DO WHILE (.TRUE.) 
	        READ(11, *, IOSTAT=iostatus) i, k, dummy1
            
            IF (iostatus == -1) exit !end of file
		    IF (iostatus /= 0) THEN
		        WRITE(*,'(a,a,a)') 'WARNING: format error in '//trim(lake_conds_file)//', line skipped.'
            ENDIF
        
            subbas_id = id_ext2int(i, id_subbas_extern) !convert external to internal id
			if (subbas_id < 1 .OR. subbas_id > subasin) then
				WRITE(*,'(a,i0,a)') 'WARNING: unknown subbasin ',i,' in '//trim(lake_conds_file)//', ignored.'
                cycle 
            end if
                
            if (k < 1 .OR. k > 5) then
				WRITE(*,'(a,i0,a)') 'WARNING: unknown reservoir class ',k,' in '//trim(lake_conds_file)//', ignored.'
                cycle 
			end if
            lakewater_hrr(1,subbas_id,k) = dummy1
        ENDDO
        close(11)
        
        DO sb_counter=1,subasin
            DO acud_class=1,5
                IF (lakewater_hrr(1,sb_counter,acud_class) < 0.) then
                    WRITE(*,'(a,a,a)') 'WARNING: Problem with state variable file ''',trim(lake_conds_file),&
                    '''. No specification for subbasin ', id_subbas_extern(sb_counter),&
                    ', reservoir size class ',acud_class,' found. Using fraction specified in lake.dat'
                END IF    
            ENDDO
        END DO  
        
        
        
    end subroutine init_lake_conds

end module model_state_io
