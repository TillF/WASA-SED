module model_state_io
    !contains subroutines for saving and loading model state (soil water, ground water, interception)


    use common_h
    !character(len=max_path_length) :: source_dir !source directory, either "input" or "output"
    integer :: file_handle !handle to currently treated input file, either in  "input/init_conds" or "output"


contains
    function get_file_handle(input_dir, output_dir, file_to_load, file_descr)
    !search for requested file in input and output dir and return the first match
        implicit none
        character(len=*), intent(in) :: input_dir, output_dir        !input and output directories
        character(len=*), intent(in):: file_to_load        !input and output directories
        !character(len=max_path_length) :: get_file_handle !chosen directory string
        character(len=max_path_length) :: src_dir !chosen directory string
        integer :: get_file_handle !handle to currently treated input file
        character(len=*) :: file_descr !description for file name
        integer :: iostate
        iostate=0

        get_file_handle = 11 !default file handle to use
        src_dir = trim(output_dir) !try output directory first
        OPEN(get_file_handle,FILE=trim(src_dir)//file_to_load, STATUS='old', action='read', IOSTAT=iostate)    !check existence of file
        if (iostate/=0) then
            src_dir = trim(input_dir)//'init_conds/' !nothing found? try input dir
            OPEN(get_file_handle,FILE=trim(src_dir)//file_to_load, STATUS='old', action='read', IOSTAT=iostate)    !check existence of file
            if (iostate/=0) then
                get_file_handle =0
                write(*,'(a,a,a)')' WARNING: '//file_descr//' file ''',trim(file_to_load),''' not found, using defaults.'
                return
            end if
        end if
        write(*,'(a,a,a)')' ... '//file_descr//' from file ''',trim(src_dir)//file_to_load,''''

    end function


    subroutine init_hillslope_state        !load initial conditions for hillslopes and reservoirs
            write(*,"(A)") "Initialize hillslope entities..."
            !if (.not. doloadstate) return   !do not load files, if disabled
            call init_soil_conds      ('soil_moisture.stat')    !Till: load initial status of soil moisture
            call init_intercept_conds ('intercept_storage.stat')    !Till: load initial status of gw storage
            call init_gw_conds        ('gw_storage.stat')    !Till: load initial status of gw storage
            call init_interflow_conds ('interflow_storage.stat')    !Till: load initial status of interflow storage
            call init_snow_conds      ('snow_storage.stat')    !Till: load initial status of snow storage
    end subroutine init_hillslope_state

   subroutine init_river_state        !load initial conditions for riverscape
            if (.not. doloadstate) return   !do not load files, if disabled
            write(*,"(A)") "Initialize river entities..."
            call init_river_conds('river_storage.stat')    !Jose Miguel: load initial status of river storage
            if (dosediment) then
               call init_sediment_conds     ('sediment_storage.stat')    !Jose Miguel: load initial status of deposited sediment storage
               call init_susp_sediment_conds('susp_sediment_storage.stat')    !Jose Miguel: load initial status of suspended sediment storage
            endif
   end subroutine init_river_state
   
    subroutine init_reservoir_lake_state        !load initial conditions for small and large reservoirs
        if (.not. doloadstate) return   !do not load files, if disabled
        write(*,"(A)") "Initialize reservoir entities..."
        call init_reservoir_conds ('reservoir_storage.stat')    !load initial status of reservoir storage
        call init_lake_conds      ('lake_storage.stat')    !Jose Miguel: load initial status of lake storage
    end subroutine init_reservoir_lake_state

    subroutine save_model_state(backup_files, start)        !save all model state variables, optionally backup older files
    implicit none
    logical, intent(in) :: backup_files
    logical, intent(in) :: start

        if (.not. dosavestate) return
            !keep files with initial conditions (if existing), save new summary on initial storages
        CALL save_all_conds(&
        rename_or_return('soil_moisture.stat'        , backup_files),&
        rename_or_return('gw_storage.stat'           , backup_files),&
        rename_or_return('intercept_storage.stat'    , backup_files),&
        rename_or_return('interflow_storage.stat'    , backup_files),&
        rename_or_return('snow_storage.stat'         , backup_files),&
        rename_or_return('lake_storage.stat'          , backup_files),&
        rename_or_return('reservoir_storage.stat'     , backup_files),&
        rename_or_return('river_storage.stat'        , backup_files),&
        rename_or_return('sediment_storage.stat'     , backup_files),&
        rename_or_return('susp_sediment_storage.stat', backup_files),&
        trim(pfadn)//'storage.stats', start)        !Till: save only summary on initial storage

    contains
    function rename_or_return(src, backup)
    !creates a backup of the existing file, if backup=TRUE and the file exists and returns '' (so save_all_conds doesn't treat this file).
    !Otherwise, returns file src
        implicit none
        character(*), intent(in) :: src
        LOGICAL, intent(in) :: backup
        character(len=len(src)+len(trim(pfadn))) :: rename_or_return
        LOGICAL :: file_exists

        rename_or_return = trim(pfadn)//src !default: return name of file
        if (.not. backup) return !no backup requested

        INQUIRE(FILE=src, EXIST=file_exists)
        if (.not. file_exists) return !file to backup not found, so return its name

        call rename(trim(pfadn)//src     , trim(pfadn)//src//'_start') !rename file
        rename_or_return = '' !no further wrting of this file is required
    end function
    end subroutine save_model_state



    subroutine save_all_conds(soil_conds_file, gw_conds_file, ic_conds_file, interflow_conds_file, snow_conds_file, lake_conds_file, reservoir_conds_file, river_conds_file, sediment_conds_file, susp_sediment_conds_file, summary_file, start)
        !store current conditions of soil moisture, ground water, interception, snow, lakes, reservoirs, river(water + sediment) in the specified files
        use hymo_h
        use params_h
        use erosion_h
        use utils_h
        use routing_h
        use lake_h
        use reservoir_h
        use time_h
        use common_h
        use snow_h

        implicit none

        character(len=*),intent(in):: soil_conds_file, gw_conds_file, ic_conds_file,interflow_conds_file, snow_conds_file, lake_conds_file, reservoir_conds_file, river_conds_file,sediment_conds_file,susp_sediment_conds_file,summary_file        !files to save to
        logical,intent(in), optional:: start

        INTEGER :: sb_counter,lu_counter,tc_counter,svc_counter,h,acud_class, i, k, tt, digits    ! counters
        INTEGER :: i_lu,id_tc_type,i_svc,i_soil,i_veg        ! ids of components in work
        INTEGER :: tcid_instance    !(internal) id of TC-instance (unique subbas-LU-TC-combination)
        REAL    :: total_storage_soil, total_storage_gw, total_storage_intercept, total_storage_lake(5), total_storage_river, total_storage_interflow, total_storage_snow, total_storage_reservoir  !total amount of water stored  [m3]
        !, total_storage_sediment,total_storage_suspsediment    
        REAL    :: lu_area, svc_area, rtemp    !area of current lu/svc [m3]
        INTEGER    ::    soil_file_hdle, gw_file_hdle, intercept_file_hdle, lake_file_hdle, reservoir_file_hdle, river_file_hdle,sediment_file_hdle,susp_sediment_file_hdle, interflow_file_hdle, snow_file_hdle    !file handles to output files
		character(len=1000) :: fmtstr, fmtstr2    !string for formatting file output
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
        total_storage_interflow=0.
        total_storage_snow=0.
        total_storage_reservoir=0.
!        total_storage_sediment=0.
!        total_storage_suspsediment=0.

!write file headers
        if (trim(soil_conds_file)=='') then        !don't do anything if an empty filename is specified
            soil_file_hdle=0
        else
            soil_file_hdle=11
            OPEN(soil_file_hdle,FILE=soil_conds_file//trim(suffix), STATUS='replace')
            WRITE(soil_file_hdle,'(a)') 'soil moisture status (for analysis or model re-start)'
            WRITE(soil_file_hdle,'(13a)')'Subbasin', char(9),'LU', char(9),'TC' , char(9),'SVC' , char(9),'horizon', char(9),&
                'watercontent_[mm]', char(9),'area_[m²]'        !tab separated output
        end if

        if (trim(gw_conds_file)=='') then        !don't do anything if an empty filename is specified
            gw_file_hdle=0
        else
            gw_file_hdle=12
            OPEN(gw_file_hdle,FILE=gw_conds_file//trim(suffix), STATUS='replace')
            WRITE(gw_file_hdle,'(a)') 'ground water storage (for analysis or model re-start)'
            WRITE(gw_file_hdle,'(7a)')'Subbasin', char(9),'LU', char(9),'volume_[mm]', char(9),'area_[m²]'        !tab separated output
        end if

        if (trim(ic_conds_file)=='') then        !don't do anything if an empty filename is specified
            intercept_file_hdle=0
        else
            intercept_file_hdle=13
            OPEN(intercept_file_hdle,FILE=ic_conds_file//trim(suffix), STATUS='replace')
            WRITE(intercept_file_hdle,'(a)') 'interception storage (for analysis or model re-start)'
            WRITE(intercept_file_hdle,'(11a)')'Subbasin', char(9),'LU', char(9),'TC' , char(9),'SVC' , char(9),&
                'storage_[mm]', char(9),'area_[m²]'        !tab separated output
        end if

        if (trim(lake_conds_file)=='') then        !don't do anything if an empty filename is specified
            lake_file_hdle=0
        else
            lake_file_hdle=15
            OPEN(lake_file_hdle,FILE=lake_conds_file//trim(suffix), STATUS='replace')
            if (doacud) then
                WRITE(lake_file_hdle,'(a)') 'Lake volume status (for analysis or model re-start)'
                WRITE(lake_file_hdle,'(a)')'Subbasin'//char(9)//'reservoir_size_class'//char(9)//'volume[m^3]' !tab separated output
            else
                CLOSE(lake_file_hdle,status='delete')
            end if
        end if


        if (trim(sediment_conds_file)=='' .or. .not. dosediment .or. (river_transport == 1) ) then        !don't do anything if an empty filename is specified
            sediment_file_hdle=0
        else
            sediment_file_hdle=16
            OPEN(sediment_file_hdle,FILE=sediment_conds_file//trim(suffix), STATUS='replace')
            WRITE(sediment_file_hdle,'(a)') 'river reach deposition sediment weight status (for analysis or model re-start)'
            WRITE(sediment_file_hdle,'(A)')'Subbasin'//char(9)//'particle_size_class'//char(9)//'mass[t]' !tab separated output
        endif

        if (trim(susp_sediment_conds_file)=='' .or. .not. dosediment) then        !don't do anything if an empty filename is specified
            susp_sediment_file_hdle=0
        else
            susp_sediment_file_hdle=17
            OPEN(susp_sediment_file_hdle,FILE=susp_sediment_conds_file//trim(suffix), STATUS='replace')
            WRITE(susp_sediment_file_hdle,'(a)') 'river reach suspended sediment weight status (for analysis or model re-start)'
            WRITE(susp_sediment_file_hdle,'(A)')'Subbasin'//char(9)//'particle_size_class'//char(9)//'mass[t]' !tab separated output
        endif
        
        if (trim(reservoir_conds_file)=='') then        !don't do anything if an empty filename is specified
            reservoir_file_hdle=0
        else
            reservoir_file_hdle=18
            OPEN(reservoir_file_hdle,FILE=reservoir_conds_file//trim(suffix), STATUS='replace')
            if (doreservoir) then
                WRITE(reservoir_file_hdle,'(a)') 'Reservoir volume status (for analysis or model re-start)'
                WRITE(reservoir_file_hdle,'(a)')'Subbasin'//char(9)//'volume[m^3]' !tab separated output
                
                tt = (d-2)*nt+hour !index to last valid value
                if (tt<1) tt=1 !Till: dirty fix to prevent crash at start up. José, please check this
                
                digits=ceiling(log10(max(1.0,maxval(abs(volact(tt,:)))*1.e6)))+2    !Till: number of pre-decimal digits required
                if (digits<10) then
                    write(fmtstr,'(A1,i0,a1,i0)') 'F',min(11,digits+4),'.',min(3,11-digits-1)        !generate format string
                else
                    fmtstr='E12.5' !for large numbers, use exponential notation
                end if
                write(fmtstr2,*) '(I0,',1,'(A1,',trim(fmtstr),'))'        !generate format string

                DO sb_counter=1,subasin 
                    IF (res_flag(sb_counter)) then
                        WRITE(reservoir_file_hdle,trim(fmtstr2))id_subbas_extern(sb_counter), char(9), volact(tt,sb_counter)*1.e6 !tab separated output
                    end if
                END DO
                total_storage_reservoir=sum(volact(tt,:)) *1.e6
                CLOSE(reservoir_file_hdle, iostat=i_lu)    !close output file
            else
                CLOSE(reservoir_file_hdle,status='delete')
            end if
        end if


!write file contents
        !write river storage state to .stat file .
        if (trim(river_conds_file)=='') then        !don't do anything if an empty filename is specified
            river_file_hdle=0
        else
            river_file_hdle=14
            OPEN(river_file_hdle,FILE=river_conds_file//trim(suffix), STATUS='replace')
            if (river_transport == 1) then !UHG routing
                WRITE(river_file_hdle,'(a)') 'UHG routing: values of routed discharge per timestep of unit hydrograph'
                WRITE(river_file_hdle,'(a)')'Subbasin'// char(9)//'[n_h x timestep]' !tab separated output
            else !Muskingum routing
                WRITE(river_file_hdle,'(a)') 'Muskingum routing: river reach volume status (for analysis or model re-start)'
                WRITE(river_file_hdle,'(a)')'Subbasin'//char(9)//'volume[m^3]' !tab separated output
            end if

            if (river_transport == 1) then !UHG routing
                tt = size(hrout,dim=1)-1 !length of UHG minus 1
                digits=ceiling(log10(max(1.0,maxval(qout(d:d+tt,1:subasin)))))+1    !Till: number of pre-decimal digits required
                if (digits<10) then
                    write(fmtstr,'(A1,i0,a1,i0)') 'F',min(11,digits+4),'.',min(3,11-digits-1)        !generate format string
                else
                    fmtstr='E12.5' !for large numbers, use exponential notation
                end if
                write(fmtstr2,*) '(I0,',tt,'(A1,',trim(fmtstr),'))'        !generate format string
                DO sb_counter=1,subasin !we wrap subbasin loop around each single entity to save time in generating format string
                    WRITE(river_file_hdle,trim(fmtstr2))id_subbas_extern(sb_counter), (char(9),qout(d+k-1,sb_counter), k=1, tt) !tab separated output
                END DO
                total_storage_river=sum(qout(d:d+tt-1,1:subasin)) * 3600 * dt !sum up total storage, convert m3/s to m3
            else !Muskingum routing
                digits=ceiling(log10(max(1.0,maxval(r_storage))))+1    !Till: number of pre-decimal digits required
                if (digits<10) then
                    write(fmtstr,'(a,i0,a,i0,a)') '(I0,A1,F',min(11,digits+4),'.',min(3,11-digits-1),')'        !generate format string
                else
                    fmtstr='(I0,A1,E12.5)' !for large numbers, use exponential notation
                end if

                DO sb_counter=1,subasin !we wrap subbasin loop around each single entity to save time in generating format string
                    WRITE(river_file_hdle,trim(fmtstr))id_subbas_extern(sb_counter), char(9),r_storage(sb_counter) !tab separated output
                END DO
                total_storage_river=sum(r_storage(1:subasin)) !sum up total storage
            end if !routing options

            CLOSE(river_file_hdle, iostat=i_lu)    !close output file
        endif


        !write interflow storage state to .stat file
        if (trim(river_conds_file)=='') then        !don't do anything if an empty filename is specified
            interflow_file_hdle=0
        else
            interflow_file_hdle=14
            OPEN(interflow_file_hdle,FILE=interflow_conds_file//trim(suffix), STATUS='replace')
            WRITE(interflow_file_hdle,'(a)') 'interflow storage (for analysis or model re-start)'
            WRITE(interflow_file_hdle,'(9a)')'Subbasin', char(9),'LU', char(9),'TC' , char(9),'horizon' , char(9),&
                'storage_[m3]'        !tab separated output

            digits=ceiling(log10(max(1.0,maxval(latred))))+1    !Till: number of pre-decimal digits required
            if (digits<10) then
                write(fmtstr,'(A1,i0,a1,i0)') 'F',min(11,digits+4),'.',min(3,11-digits-1)        !generate format string
            else
                fmtstr='E12.5' !for large numbers, use exponential notation
            end if
            tt = size(latred,dim=2) !number of interchange horizons in latred
            write(fmtstr2,*) '(4(I0,A1),',trim(fmtstr),')'        !generate format string
            DO sb_counter=1,subasin
                DO lu_counter=1,nbr_lu(sb_counter)
                    i_lu=id_lu_intern(lu_counter,sb_counter)
                    DO tc_counter=1,nbrterrain(i_lu) !intercept and soil storages inside
                        tcid_instance=tcallid(sb_counter,lu_counter,tc_counter)    !id of TC instance
                        if (tcid_instance==-1) cycle                            !this may happen if this is merely a dummy basin with prespecified outflow
                        id_tc_type=id_terrain_intern(tc_counter,i_lu)            !id of TC type
                        DO  h=1,tt
                            WRITE(interflow_file_hdle,fmtstr2) id_subbas_extern(sb_counter), char(9),&
                                    id_lu_extern(i_lu),char(9),id_terrain_extern(id_tc_type), &
                                    char(9), h, char(9), latred(tcid_instance,h)   !tab separated output
                        END DO
                    ENDDO    !loop TCs
                ENDDO    !loop LUs
            ENDDO    !loop subbasins
            total_storage_interflow = sum(latred) !sum up total storage
            CLOSE(interflow_file_hdle, iostat=i_lu)    !close output file
        endif !interflow output

        !write snow storage state to .stat file
        if (trim(snow_conds_file)=='' .OR. .not. dosnow) then        !don't do anything if an empty filename is specified
            snow_file_hdle=0
        else
            snow_file_hdle=19
            OPEN(snow_file_hdle,FILE=snow_conds_file//trim(suffix), STATUS='replace')
            WRITE(snow_file_hdle,'(a)') 'snow storage (for analysis or model re-start)'
            WRITE(snow_file_hdle,'(11a)')'Subbasin', char(9),'LU', char(9), 'TC' , char(9),&
                'storage [m]', char(9), 'energy [kJ/m²]', char(9), 'albedo [-]'         !tab separated output

            digits=ceiling(log10(max(1.0, maxval(snowWaterEquiv), maxval(snowEnergyCont), maxval(snowAlbedo))))+1    !Till: number of pre-decimal digits required
            if (digits<10) then
                write(fmtstr,'(A1,i0,a1,i0)') 'F',min(11,digits+4),'.',min(3,11-digits-1)        !generate format string
            else
                fmtstr='E12.5' !for large numbers, use exponential notation
            end if

            write(fmtstr2,*) '(2(I0,A1),I0,3(A1,',trim(fmtstr),'))'        !generate format string
            tt = max(2, d) -1  !so it works both at the beginning and end of simulation year/period
            DO sb_counter=1,subasin
                DO lu_counter=1,nbr_lu(sb_counter)
                    i_lu=id_lu_intern(lu_counter,sb_counter)
                    lu_area=area(sb_counter)*frac_lu(lu_counter,sb_counter)*1e6
                    DO tc_counter=1,nbrterrain(i_lu) !intercept and soil storages inside
                        tcid_instance=tcallid(sb_counter,lu_counter,tc_counter)    !id of TC instance
                        if (tcid_instance==-1) cycle                            !this may happen if this is merely a dummy basin with prespecified outflow
                        id_tc_type=id_terrain_intern(tc_counter,i_lu)            !id of TC type
                            WRITE(snow_file_hdle,fmtstr2) id_subbas_extern(sb_counter), char(9),&
                                    id_lu_extern(i_lu),char(9),id_terrain_extern(id_tc_type), &
                                    char(9), snowWaterEquiv(tt, nt, tcid_instance), &
                                    char(9), snowEnergyCont(tt, nt, tcid_instance), &
                                    char(9), snowAlbedo    (tt, nt, tcid_instance)
                        total_storage_snow = total_storage_snow + snowWaterEquiv(tt, nt, tcid_instance)*lu_area*fracterrain(id_terrain_intern(tc_counter,i_lu)) !sum up total storage [m3]
                    ENDDO    !loop TCs
                ENDDO    !loop LUs
            ENDDO    !loop subbasins

            CLOSE(snow_file_hdle, iostat=i_lu)    !close output file
        endif !snow output


        if (dosediment) then
            !write riverbed sediment storage
            !generate format string
            if (sediment_file_hdle/=0) then
                digits=ceiling(log10(max(1.0,maxval(riverbed_storage))))+1    !Till: number of pre-decimal digits required
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
            CLOSE(sediment_file_hdle, iostat=i_lu)    !close output file

            !write suspended sediment storage
            if (susp_sediment_file_hdle/=0) then
                digits=ceiling(log10(max(1.0,maxval(sed_storage))))+1    !Till: number of pre-decimal digits required
                if (digits<10) then
                  write(fmtstr,'(a,i0,a,i0,a)') '(I0,A1,I0,A1,F',min(11,digits+4),'.',min(3,11-digits-1),'))'        !generate format string
                else
                   fmtstr='(I0,A1,I0,A1,E12.5)' !for large numbers, use exponential notation
                end if
            end if
            DO sb_counter=1,subasin !we wrap subbasin loop around each single entity to save time in generating format string
                if (susp_sediment_file_hdle/=0) then
                    if (river_transport == 1) then !UHG routing
                        tt = size(hrout,dim=1)-1 !length of UHG minus 1
                        digits=ceiling(log10(max(1.0,maxval(qsediment2_t(d:d+tt,1,1:subasin)))))+1    !Till: number of pre-decimal digits required
                        if (digits<10) then
                            write(fmtstr,'(A1,i0,a1,i0)') 'F',min(11,digits+4),'.',min(3,11-digits-1)        !generate format string
                        else
                            fmtstr='E12.5' !for large numbers, use exponential notation
                        end if
                        write(fmtstr2,*) '(I0,A1,I0,',tt,'(A1,',trim(fmtstr),'))'        !generate format string
                        
                        do k=1, n_sed_class
                             WRITE(susp_sediment_file_hdle,fmtstr2)id_subbas_extern(sb_counter), char(9), k, (char(9), qsediment2_t(d+i-1,1,sb_counter)/n_sed_class, i=1, tt) !print each sediment class
                        enddo
                    else !Muskingum routing
                        do k=1, n_sed_class
                             WRITE(susp_sediment_file_hdle,fmtstr)id_subbas_extern(sb_counter), char(9), k, char(9), sed_storage(sb_counter,k) !print each sediment class
                        enddo
                    end if    
                endif
                !total_storage_suspsediment=total_storage_suspsediment+sum(sed_storage(sb_counter,:)) !sum up total storage
            END DO
            CLOSE(susp_sediment_file_hdle, iostat=i_lu)    !close output file
        end if !dosediment


        IF (doacud) THEN
            tt = (d-2)*nt+hour !index to last valid value
            if (tt<1) tt=1 !Till: dirty fix to prevent crash at start up. José, please check this
            DO sb_counter=1,subasin
                DO acud_class=1,5
				    if (lake_file_hdle/=0) then
					    WRITE(lake_file_hdle,'(I0,A1,I0,A1,F13.2)') id_subbas_extern(sb_counter), char(9),acud_class,char(9),&
						    lakewater_hrr(tt,sb_counter,acud_class)
				    endif
				    total_storage_lake(acud_class) = total_storage_lake(acud_class)+lakewater_hrr(tt,sb_counter,acud_class) * acud(sb_counter,acud_class) !sum up total storage
                ENDDO
            END DO
        END IF !small reservoirs
        CLOSE(lake_file_hdle, iostat=i_lu)


!save groundwater, intercept and soil storages
        DO sb_counter=1,subasin
            DO lu_counter=1,nbr_lu(sb_counter)
                i_lu=id_lu_intern(lu_counter,sb_counter)
                lu_area=area(sb_counter)*frac_lu(lu_counter,sb_counter)*1e6
                if (gw_file_hdle/=0) then
                    if (lu_area /= 0.) then
                        rtemp = deepgw(sb_counter,lu_counter)/lu_area*1e3 
                    else
                        rtemp = 0. !in case of dummy subbasins
                    end if    
                    WRITE(gw_file_hdle,'(2(I0,A1),F8.2,A1,F0.0)') id_subbas_extern(sb_counter), char(9),id_lu_extern(i_lu),&
                        char(9), rtemp,&
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



        OPEN(11,FILE=summary_file//trim(suffix), STATUS='replace')        !write to summary file
        WRITE(11,'(A)')'total water storage in catchment after model run [m3]'
        WRITE(11,'(A,F16.1)')'soil_storage'        // char(9), total_storage_soil
        WRITE(11,'(A,F16.1)')'gw_storage'          // char(9), total_storage_gw
        WRITE(11,'(A,F16.1)')'interception_storage'// char(9), total_storage_intercept
        WRITE(11,'(A,F16.1)')'interflow_storage'   // char(9), total_storage_interflow
        WRITE(11,'(A,F16.1)')'snow_storage'        // char(9), total_storage_snow

        DO acud_class=1,5
            WRITE(11,'(A,I0,A,F16.1)')'lake_storage', acud_class, char(9),total_storage_lake(acud_class)
        ENDDO
        WRITE(11,'(A,F16.1)')'river_storage'       // char(9), total_storage_river
        WRITE(11,'(A,F16.1)')'reservoir_storage'   // char(9), total_storage_reservoir

        CLOSE(11)
    end subroutine save_all_conds

    subroutine init_soil_conds(soil_conds_file)

        !load soil moistures information from file soil_conds_file
        use hymo_h
        use params_h
        use utils_h
        use erosion_h
        implicit none

        character(len=*), intent(in):: soil_conds_file        !file to load from
        INTEGER :: i,line,errors,sb_counter,lu_counter,tc_counter,svc_counter,h    ! counters
        INTEGER :: i_subbasx,i_lux,i_tcx,i_svcx !, i_soilx        ! external ids of components in work
        INTEGER :: i_subbas,i_lu,i_tc,i_svc, i_soil, i_veg,id_tc_type        ! internal ids of components in work
        INTEGER :: tcid_instance !,lu_instance, !,soil_instance    !(internal) id of LU,TC,soil-instance (unique subbas-LU-TC-soil-combination)
        REAL    :: horithact_temp, x
        INTEGER    :: file_read=0
        character(len=160) :: linestr='', error_msg=''

        horithact=-9999.                    !mark all horizons as "not (yet) initialised"

        file_handle = get_file_handle(pfadp, pfadn, soil_conds_file, "soil moisture") !pick source directory

        if (file_handle /= 0) then        !load values from file

            READ(file_handle,*, IOSTAT=i); READ (11,*, IOSTAT=i)    !skip header lines
            line=2
            errors=0

            do while (.TRUE.)        !read whole file
                IF (len(trim(error_msg))/=0) THEN    !print error message, if occured
                    if (errors==0) then !print heading at before first error
                        write(*,'(A,/,6a12)')'  Entities not found in current domain (ignored):','Line','subbasin','LU','TC','SVC',&
                            'horizon'
                    end if
                    write(*,*)trim(error_msg)
                    error_msg='' !reset error message
                    errors=errors+1
                END IF

                READ(file_handle,'(A)',  IOSTAT=i) linestr

                IF (i==24 .OR. i==-1) THEN    !end of file
                    exit        !exit loop
                END IF
                line=line+1

                READ(linestr,*,  IOSTAT=i) i_subbasx, i_lux, i_tcx, i_svcx, h ,horithact_temp, x
                IF (i/=0 .OR. isnan(0. + i_subbasx+i_lux+i_tcx+i_svcx+h+horithact_temp+x) ) THEN    !format error
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
                    horiz_thickness(tcid_instance,svc_counter,h)*1.01) then            !exceeds storage capacity
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
            CLOSE(file_handle)
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
                                        write(*,'(A,f5.2,a,/,5a12)')' Following entities not initialised, using defaults '// &
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

  subroutine init_interflow_conds(interflow_conds_file)

        !load interflow information from file interflow_conds_file
        use hymo_h
        use params_h
        use utils_h
        implicit none

        character(len=*),intent(in):: interflow_conds_file        !file to load from
        INTEGER :: i,line,errors,sb_counter,lu_counter,tc_counter,h    ! counters
        INTEGER :: i_subbasx,i_lux,i_tcx       ! external ids of components in work
        INTEGER :: i_subbas,i_lu,i_tc,  id_tc_type, tt        ! internal ids of components in work
        INTEGER :: tcid_instance     !(internal) id of LU,TC,soil-instance (unique subbas-LU-TC-soil-combination)
        REAL    :: horithact_temp
        INTEGER    :: file_read=0
        character(len=160) :: linestr='', error_msg=''

        tt = size(latred,dim=2) !number of interchange horizons in latred
        latred=-9999.                    !mark all interchange horizons as "not (yet) initialised"

        file_handle = get_file_handle(pfadp, pfadn, interflow_conds_file, "interflow") !pick source directory

        if (file_handle /= 0) then        !load values from file
            READ(file_handle,*); READ (file_handle,*)    !skip header lines
            line=2
            errors=0

            do while (.TRUE.)        !read whole file
                IF (len(trim(error_msg))/=0) THEN    !print error message, if occured
                    if (errors==0) then !print heading at before first error
                        write(*,'(A,/,5a12)')' Entities not found in current domain (ignored):','Line','subbasin','LU','TC','horizon'
                    end if
                    write(*,*)trim(error_msg)
                    error_msg='' !reset error message
                    errors=errors+1
                END IF

                READ(file_handle,'(A)',  IOSTAT=i) linestr

                IF (i==24 .OR. i==-1) THEN    !end of file
                    exit        !exit loop
                END IF
                line=line+1

                READ(linestr,*,  IOSTAT=i) i_subbasx,i_lux,i_tcx,h,horithact_temp
                IF (i/=0  .OR. isnan(0. + i_subbasx+i_lux+i_tcx+h+horithact_temp)) THEN    !format error
                    write(*,'(A,i0,A2,a)')' Format error in line ', line,': ',trim(linestr)
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

                if (h> tt) then
                    write(error_msg,'(5i12)')line,i_subbasx,i_lux,i_tcx,h
                    cycle    !proceed with next line
                end if

                if (horithact_temp < 0) then    !ignore negative values
                    write(error_msg,'(5i12)')line,i_subbasx,i_lux,i_tcx,h
                    cycle    !proceed with next line
                end if

                latred(tcid_instance,h)=horithact_temp            !set interchange horizon water content

            END DO
            file_read=1
            CLOSE(file_handle)
        end if

        errors=0

        DO sb_counter=1,subasin            !check, if all relevant interchange horizons have been initialized, if not, use default values
            DO lu_counter=1,nbr_lu(sb_counter)
                i_lu=id_lu_intern(lu_counter,sb_counter)
                DO tc_counter=1,nbrterrain(i_lu)
                    tcid_instance=tcallid(sb_counter,lu_counter,tc_counter) !id of TC instance
                    if (tcid_instance==-1) cycle                            !this may happen if this is merely a dummy basin with prespecified outflow
                    id_tc_type=id_terrain_intern(tc_counter,i_lu)            !id of TC type

                    DO h=1,tt
                        if (latred(tcid_instance,h)==-9999.) then            !not yet set?
                            if (file_read==1) then                        !but this should have been done before
                                if (errors==0) then    !produce header before first warning only
                                    write(*,'(a,/,4a12)')' Following entities not initialised, using defaults (content=0):', &
                                    'subbasin','LU','TC','horizon'
                                end if
                                errors=errors+1
                                write(*,'(4i12)')id_subbas_extern(sb_counter), id_lu_extern(i_lu),&
                                    id_terrain_extern(id_tc_type), h            !issue warning
                            end if

                            latred(tcid_instance,h)=0.  !resume to default (0)
                        end if
                    END DO
                END DO
            END DO
        END DO
        where (latred==-9999.) latred=0.    !remove unchecked* "not-read" markers to avoid later confusion (*: may happen for subbasins with prespecified outflow)


end subroutine init_interflow_conds

  subroutine init_snow_conds(snow_conds_file)

        !load snow information from file snow_conds_file
        use hymo_h
        use utils_h
        use snow_h
        use params_h
        implicit none

        character(len=*),intent(in):: snow_conds_file        !file to load from
        INTEGER :: i,line,errors,sb_counter,lu_counter,tc_counter    ! counters
        INTEGER :: i_subbasx,i_lux,i_tcx       ! external ids of components in work
        INTEGER :: i_subbas,i_lu,i_tc,  id_tc_type        ! internal ids of components in work
        INTEGER :: tcid_instance     !(internal) id of LU,TC,soil-instance (unique subbas-LU-TC-soil-combination)
        REAL    :: snowweqv_temp,snowetemp,albedo_temp
        INTEGER    :: file_read=0
        character(len=160) :: linestr='', error_msg=''

        if (.not. dosnow) return
        snowWaterEquiv=-9999.                    !mark all snow storages as "not (yet) initialised"
        file_handle = get_file_handle(pfadp, pfadn, snow_conds_file, "snow") !pick source directory

        if (file_handle /= 0) then        !load values from file
            write(*,'(a,a,a)')' ... snow from file ''',trim(snow_conds_file),''''

            READ(file_handle,*); READ (file_handle,*)    !skip header lines
            line=2
            errors=0

            do while (.TRUE.)        !read whole file
                IF (len(trim(error_msg))/=0) THEN    !print error message, if occured
                    if (errors==0) then !print heading at before first error
                        write(*,'(A,/,5a12)')' Entities not found in current domain (ignored):','Line','subbasin','LU','TC'
                    end if
                    write(*,*)trim(error_msg)
                    error_msg='' !reset error message
                    errors=errors+1
                END IF

                READ(file_handle,'(A)',  IOSTAT=i) linestr

                IF (i==24 .OR. i==-1) THEN    !end of file
                    exit        !exit loop
                END IF
                line=line+1

                READ(linestr,*,  IOSTAT=i) i_subbasx,i_lux,i_tcx,snowweqv_temp,snowetemp,albedo_temp
                IF (i/=0 .OR. isnan(0. + i_subbasx+i_lux+i_tcx+snowweqv_temp+snowetemp+albedo_temp)) THEN    !format error
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

                snowWaterEquiv(1,1,tcid_instance) = snowweqv_temp
                snowEnergyCont(1,1,tcid_instance) = snowetemp
                snowAlbedo    (1,1,tcid_instance) = albedo_temp

            END DO
            file_read=1
            CLOSE(file_handle)
        end if

        errors=0

        DO sb_counter=1,subasin            !check, if all relevant interchange horizons have been initialized, if not, use default values
            DO lu_counter=1,nbr_lu(sb_counter)
                i_lu=id_lu_intern(lu_counter,sb_counter)
                DO tc_counter=1,nbrterrain(i_lu)
                    tcid_instance=tcallid(sb_counter,lu_counter,tc_counter) !id of TC instance
                    if (tcid_instance==-1) cycle                            !this may happen if this is merely a dummy basin with prespecified outflow
                    id_tc_type=id_terrain_intern(tc_counter,i_lu)            !id of TC type

                    if (snowWaterEquiv(1,1,tcid_instance)==-9999.) then            !not yet set?
                        if (file_read==1) then                        !but this should have been done before
                            if (errors==0) then    !produce header before first warning only
                                write(*,'(A,/,4a12)')' Following entities not initialised, using defaults '// &
                                    '(content=0):','subbasin','LU','TC'
                            end if
                            errors=errors+1
                            write(*,'(4i12)')id_subbas_extern(sb_counter), id_lu_extern(i_lu),&
                                id_terrain_extern(id_tc_type)            !issue warning
                        end if
                        snowWaterEquiv(1,1,tcid_instance) = 0.
                        snowEnergyCont(1,1,tcid_instance) = 0.
                        snowAlbedo    (1,1,tcid_instance) = 0.

                        !resume to default (0)
                    end if
                END DO
            END DO
        END DO


    end subroutine init_snow_conds


    subroutine init_intercept_conds(intercept_conds_file)

        !load intercept state from file intercept_conds_file
        use hymo_h
        use params_h
        use utils_h
        use erosion_h
        implicit none
        !
        character(len=*), intent(in):: intercept_conds_file        !file to load from
        INTEGER :: i,line,errors,sb_counter,lu_counter,tc_counter,svc_counter !    ! counters
        INTEGER :: i_subbasx,i_lux,i_tcx,i_svcx         ! external ids of components in work
        INTEGER :: i_subbas,i_lu,i_tc,i_svc, i_soil, i_veg,id_tc_type        ! internal ids of components in work
        INTEGER :: tcid_instance    !(internal) id of LU,TC,soil-instance (unique subbas-LU-TC-soil-combination)
        REAL    :: int_temp, x
        INTEGER    :: file_read=0
        character(len=160) :: error_msg='', linestr=''

        intercept(:,:)=-9999.                    !mark all SVC-int-storages as "not (yet) initialised"

        file_handle = get_file_handle(pfadp, pfadn, intercept_conds_file, "interception") !pick source directory

        if (file_handle /= 0) then        !load values from file
            READ(file_handle,*); READ (file_handle,*)    !skip header lines
            line=2
            errors=0

            do while (.TRUE.)        !read whole file
                IF (len(trim(error_msg))/=0) THEN    !print error message, if occured
                    if (errors==0) then !print heading at before first error
                        write(*,'(A,/,5a12)')' Entities not found in current domain (ignored):','Line','subbasin','LU','TC','SVC'
                    end if
                    write(*,*)trim(error_msg)
                    error_msg='' !reset error message
                    errors=errors+1
                END IF

                READ(file_handle,'(A)',  IOSTAT=i) linestr

                IF (i==24 .OR. i==-1) THEN    !end of file
                    exit        !exit loop
                END IF
                line=line+1

                READ(linestr,*,  IOSTAT=i) i_subbasx,i_lux,i_tcx,i_svcx,int_temp, x
                IF (i/=0 .OR. isnan(0. + i_subbasx+i_lux+i_tcx+i_svcx+int_temp+ x)) THEN    !format error
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


                if (int_temp < 0.) then            !negative interception
                    if (errors==0) then    !produce header before first warning only
                        write(*,'(A,/,4a12)')' Interception of following SVCs is negative. Corrected to 0.',&
                            'subbasin','LU','TC','SVC'
                    end if
                    errors=errors+1
                    write(*,'(4i12)')i_subbasx, i_lux, i_tcx,i_svcx            !issue warning
                    int_temp = 0.
                end if

                intercept(tcid_instance,svc_counter)=int_temp            !set interception value


            END DO
            file_read=1
            CLOSE(file_handle)
        end if

        !check, if all relevant SVCs have been initialized, if not, use default values
        errors=0
        DO sb_counter=1,subasin
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

                        if (intercept(tcid_instance,svc_counter)==-9999.) then            !not yet set?
                            if (file_read==1) then                        !but this should have been done before
                                if (errors==0) then    !produce header before first warning only
                                    write(*,'(A,f5.2,a,/,4a12)')' Following entities not initialised, using defaults '// &
                                        '(intercept=', 0.,'):','subbasin','LU','TC','SVC'
                                end if
                                errors=errors+1
                                write(*,'(4i12)')id_subbas_extern(sb_counter), id_lu_extern(i_lu),&
                                    id_terrain_extern(id_tc_type), id_svc_extern(i_svc)            !issue warning
                            end if

                            intercept(tcid_instance,svc_counter)=0.        !set to default
                        end if

                        intercept(tcid_instance, (nbr_svc(tcid_instance)+1):maxsoil)=0.        !set water content of irrelevant horizons values to 0
                    END DO
                END DO
            END DO
        END DO

        if (errors>0) then
            return
        end if

    end subroutine init_intercept_conds

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

        deepgw=-9999.                    !mark all gw storages as "not (yet) read"

        file_handle = get_file_handle(pfadp, pfadn, gw_conds_file, "GW storage") !pick source directory

        if (file_handle /= 0) then        !load values from file
            READ(file_handle,*); READ (file_handle,*)    !skip header lines
            line=2
            errors=0

            do while (.TRUE.)        !read whole file
                IF (len(trim(error_msg))/=0) THEN    !print error message, if occured
                    if (errors==0) then
                        write(*,'(A,/,3a12)')' Format errors or entities not found in current domain (ignored):','Line','subbasin','LU'
                    end if
                    write(*,*)trim(error_msg)
                    error_msg=''
                    errors=errors+1
                END IF

                READ(file_handle,'(A)',  IOSTAT=i) linestr

                IF (i==24 .OR. i==-1) THEN    !end of file
                    exit        !exit loop
                END IF
                line=line+1

                READ(linestr,*,  IOSTAT=i) i_subbasx,i_lux,gwvol_temp, x
                IF (i/=0 .OR. isnan(0. + i_subbasx+i_lux + gwvol_temp + x)) THEN    !format error
                    write(error_msg,'(i12,1a2,a)')line,' ',trim(linestr)
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

                if (lu_area /= 0) then
                    deepgw(i_subbas,lu_counter)=gwvol_temp/1000.*lu_area            !set gw water content [m3]
                else
                    deepgw(i_subbas,lu_counter)=0. !may be the case for dummy subbasins
                end if    
            END DO
            file_read=1
            CLOSE(file_handle)
            !if (errors>0) then
            !    return
            !end if
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
        CLOSE(file_handle)
    end subroutine init_gw_conds


    subroutine init_river_conds(river_conds_file)
        use routing_h
        use params_h
        use utils_h
    	use hymo_h
        implicit none

        character(len=*),intent(in):: river_conds_file        !file to load from
        integer :: subbas_id, iostatus, i, j, tt, k, line
        real :: dummy1
        real, pointer :: array_ptr(:)
        character(len=1000) :: linestr

        file_handle = get_file_handle(pfadp, pfadn, river_conds_file, "river storage") !pick source directory

        if (file_handle == 0) then        !load values from file
		    return
        end if

        !verify that this river file belongs to the selected routing mode
        READ(file_handle,'(A)') linestr;
        if ( (river_transport == 1 .AND. linestr(1:3) /= "UHG") .OR. &
                ( ((river_transport == 2) .or. (river_transport == 3)) .AND. linestr(1:3) /= "Mus") )            then
            write(*,'(a,a,a)')' WARNING: River storage file ''',trim(river_conds_file),''' does not match selected routing mode. File ignored, using defaults.'
            CLOSE(11)
		    return
        end if

        !skip next header line
        READ(file_handle,*)

	    if (river_transport == 1) then
            tt = size(hrout,dim=1)-1 !length of UHG minus 1
            qout                   = 0.
            qout(1,1:subasin)      =-1. !indicator for "not read"
            
            !check that the current file matches the specs of the current UHG
            READ(file_handle,'(A)',IOSTAT=iostatus) linestr
             IF (iostatus /=0) then
                write(*,'(a,a,a)')'WARNING: format error in river storage file ''',trim(river_conds_file),', line ',line,', using defaults.'             
             else
               line = line+2
               i = GetNumberOfSubstrings(linestr) - 1
               if (i > tt .OR. i<1) write(*,'(a,a,a,i0,a,i0)')' WARNING: River storage file ''',trim(river_conds_file),''': expecting ',tt,' ordinates, found ',i,'. Truncated.'
                backspace(11) !rewind line just read
             end if   
        end if
        if ( ((river_transport == 2) .or. (river_transport == 3))) r_storage(1:subasin)=-1. !indicator for "not read"

        line = 1 !line counter
        DO WHILE (.TRUE.)
            READ(file_handle,'(A)',IOSTAT=iostatus) linestr
            IF (iostatus /=0) exit
            line = line+2
           
            READ(linestr,*,IOSTAT=iostatus) i, dummy1
            
            if (iostatus /=0 .or. isnan(dummy1)) then
			    WRITE(*,'(a,i0,a)') ' WARNING: format error in line ', line,' of river_storage.stat, line ignored.'
                cycle
            end if

            subbas_id = id_ext2int(i, id_subbas_extern) !convert external to internal id
		    if (subbas_id < 1 .OR. subbas_id > subasin) then
			    WRITE(*,'(a,i0,a,i0,a)') ' WARNING: unknown subbasin ', i,' in line ', line, ' of river_storage.stat, ignored.'
                cycle
            end if

            if (dummy1 < 0.) then
			    WRITE(*,'(a,i0,a,i0,a)') 'ERROR: negative value for subbasin ', i,' in line ', line, ' of river_storage.stat.'
                stop
            end if

            if (river_transport == 1) then
                READ(linestr,*,IOSTAT=iostatus) k, (qout(j,subbas_id), j=1,tt)
                if (iostatus /= 0) WRITE(*,'(a,i0,a,i0,a)') 'WARNING: length of saved UHG not matching for subbasin ',i,' in line ', line, ' of river_storage.stat; truncated/padded.'
            end if    
            if (((river_transport == 2) .or. (river_transport == 3))) r_storage(subbas_id)=dummy1
        END DO
        close(11)

        !check for completeness
        if (river_transport == 1) array_ptr => qout(1,1:subasin)
        if (((river_transport == 2) .or. (river_transport == 3))) array_ptr => r_storage

        if (count(array_ptr==-1.) > 0) then
            WRITE(*,'(A)') ' WARNING: could not read initial river storage from river_storage.stat for the following subbasins, assumed 0:'
            DO subbas_id=1,subasin
                if (array_ptr(subbas_id)==-1.) then
                    WRITE(*,'(i0)') id_subbas_extern(subbas_id)
                    array_ptr(subbas_id)=0.
                end if
            END DO
        end if

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

        if (river_transport == 1) then !UHG routing does not consider deposited sediment so far
            return
        end if
        
        riverbed_storage(:,:)=-1. !indicator for "not read"

        file_handle = get_file_handle(pfadp, pfadn, sediment_conds_file, "river sediment storage") !pick source directory

        if (file_handle == 0) then        !load values from file
            riverbed_storage(:,:)=0. !set to default
            return
        end if

        !read 2 header lines into buffer
        READ(file_handle,*); READ(file_handle,*)

        DO WHILE (.TRUE.)
	        READ(file_handle,*,IOSTAT=iostatus) i, k, dummy1
            IF (iostatus == -1) exit !end of file
		    IF (iostatus /= 0) THEN
		        WRITE(*,'(a,a,a)') ' WARNING: format error in sediment_storage.stat, line skipped, assumed 0.'
            ENDIF

            subbas_id = id_ext2int(i, id_subbas_extern) !convert external to internal id
			if (subbas_id < 1 .OR. subbas_id > subasin) then
				WRITE(*,'(a,i0,a)') ' WARNING: unknown subbasin ',i,' in sediment_storage.stat, ignored.'
                cycle
            end if

            if (k < 1 .OR. k > n_sed_class) then
				WRITE(*,'(a,i0,a)') ' WARNING: unknown particle size class ',k,' in sediment_storage.stat, ignored.'
                cycle
			end if
	        riverbed_storage(subbas_id,k)=dummy1
        ENDDO
        close(11)

        if (count(riverbed_storage==-1.) > 0) then
            WRITE(*,'(A)') ' WARNING: could not read initial river sediment storage from sediment_storage.stat for the following subbasins, assumed 0:'
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

        file_handle = get_file_handle(pfadp, pfadn, susp_sediment_conds_file, "river susp. sediment storage") !pick source directory

        if (file_handle == 0) then        !load values from file
            sed_storage(:,:) = 0. !default value
            return
        end if

        !read 2 header lines into buffer
        READ(file_handle,*); READ(file_handle,*)

        DO WHILE (.TRUE.)
	        READ(file_handle,*,IOSTAT=iostatus) i, k, dummy1
            IF (iostatus == -1) exit !end of file
		    IF (iostatus /= 0) THEN
		        WRITE(*,'(a,a,a)') ' WARNING: format error in susp_sediment_storage.stat, line skipped, assumed 0.'
            ENDIF

            subbas_id = id_ext2int(i, id_subbas_extern) !convert external to internal id
			if (subbas_id < 1 .OR. subbas_id > subasin) then
				WRITE(*,'(a,i0,a)') ' WARNING: unknown subbasin ',i,' in susp_sediment_storage.stat, ignored.'
                cycle
            end if

            if (k < 1 .OR. k > n_sed_class) then
				WRITE(*,'(a,i0,a)') ' WARNING: unknown particle size class ',k,' in susp_sediment_storage.stat, ignored.'
                cycle
			end if

	        if (isnan(dummy1)) then
                WRITE(*,'(a,i0,a)') ' WARNING: NaN in susp_sediment_storage.stat, ignored.'
                cycle
            end if

	        sed_storage(subbas_id,k)=dummy1
        ENDDO
        close(11)

        if (count(sed_storage==-1.) > 0) then
            WRITE(*,'(A)') ' WARNING: could not read initial river sediment storage from susp_sediment_storage.stat for the following subbasins, assumed 0:'
            DO subbas_id=1,subasin
				if (count(sed_storage(subbas_id,:)==-1) > 0) WRITE(*,'(i0)') subbas_id
            END DO
            where(sed_storage==-1) sed_storage=0.
        end if

    end subroutine init_susp_sediment_conds


    subroutine init_lake_conds(lake_conds_file)
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
        file_handle = get_file_handle(pfadp, pfadn, lake_conds_file, "small reservoirs") !pick source directory

        if (file_handle /= 0) then        !load values from file

            READ(file_handle,*, IOSTAT=iostatus); READ(file_handle,*, IOSTAT=iostatus)!read 2 header lines

            DO WHILE (.TRUE.)
                READ(file_handle, *, IOSTAT=iostatus) i, k, dummy1

                IF (iostatus == -1) exit !end of file
                IF (iostatus /= 0) THEN
                    WRITE(*,'(a,a,a)') ' WARNING: format error in '//trim(lake_conds_file)//', line skipped.'
                ENDIF

                subbas_id = id_ext2int(i, id_subbas_extern) !convert external to internal id
                if (subbas_id < 1 .OR. subbas_id > subasin) then
                    WRITE(*,'(a,i0,a)') ' WARNING: unknown subbasin ',i,' in '//trim(lake_conds_file)//', ignored.'
                    cycle
                end if

                if (k < 1 .OR. k > 5) then
                    WRITE(*,'(a,i0,a)') ' WARNING: unknown reservoir class ',k,' in '//trim(lake_conds_file)//', ignored.'
                    cycle
                end if
                lakewater_hrr(1,subbas_id,k) = dummy1
            ENDDO
            close(11)

            DO sb_counter=1,subasin
                DO acud_class=1,5
                    IF (lakewater_hrr(1,sb_counter,acud_class) < 0.) then
                        WRITE(*,'(a,i0,a,i0,a)') ' WARNING: No specification for subbasin ',&
                         id_subbas_extern(sb_counter), ', reservoir size class ',&
                        acud_class,' found in '''//trim(lake_conds_file)//'''. Using fraction specified in lake.dat'
                    END IF
                ENDDO
            END DO
        end if !initialisation via files
        CALL lake(0, 0) !this is the old default initialisation. It also initialises some dependend variables
    end subroutine init_lake_conds

     subroutine init_reservoir_conds(reservoir_conds_file)
        use reservoir_h
        use time_h
        use params_h
        use utils_h
	    use hymo_h

        implicit none

        character(len=*),intent(in):: reservoir_conds_file        !file to load from
        integer :: sb_counter, iostatus, i, subbas_id
        logical :: reservoir_read(subasin)
        real :: dummy1

        if (.not. doreservoir .OR.& !don't try to load file if reservoirs have been disabled anyway
            .not. doloadstate) return   !do not load files, if disabled
        
        where (res_flag)
            reservoir_read = .false. !mark all reservoirs a "not read"
        elsewhere
            reservoir_read = .true. !for detecting uninitialized reservoirs later: non-existing reservoirs don't need to be red, though
        end where

        file_handle = get_file_handle(pfadp, pfadn, reservoir_conds_file, "reservoirs") !pick source directory

        if (file_handle == 0) then        !load values from file
            return
        end if

        READ(file_handle,*, IOSTAT=iostatus); READ(file_handle,*, IOSTAT=iostatus)!read 2 header lines

        DO WHILE (.TRUE.)
	        READ(file_handle, *, IOSTAT=iostatus) i, dummy1

            IF (iostatus == -1) exit !end of file
		    IF (iostatus /= 0) THEN
		        WRITE(*,'(a,a,a)') ' WARNING: format error in '//trim(reservoir_conds_file)//', line skipped.'
            ENDIF

            subbas_id = id_ext2int(i, id_subbas_extern) !convert external to internal id
			if (subbas_id < 1 .OR. subbas_id > subasin) then
				WRITE(*,'(a,i0,a)') ' WARNING: unknown subbasin ',i,' in '//trim(reservoir_conds_file)//', ignored.'
                cycle
            end if

            volact(1,subbas_id) = dummy1 / 1e6 !internally used in [10^6 m3]
            reservoir_read(subbas_id) = .true. !mark as "storage read"
        ENDDO
        close(11)

        DO sb_counter=1,subasin
             IF (.not. reservoir_read(sb_counter) ) THEN !reservoir has not been read before, but should be present
                    WRITE(*,'(a,i0,a)') ' WARNING: No specification for subbasin ',&
                    id_subbas_extern(sb_counter),' found in '''// trim(reservoir_conds_file)//'''. Using values specified in reservoir.dat'
            END IF
        END DO
    end subroutine init_reservoir_conds

    
end module model_state_io
