    SUBROUTINE readhymo

    ! Code converted using TO_F90 by Alan Miller
    ! Date: 2005-06-30  Time: 13:47:20

    use allocate_h
    use lake_h
    use common_h
    use hymo_h
    use params_h
    use time_h
    use reservoir_h
    use routing_h
    use utils_h
    use erosion_h
    use climo_h
    use snow_h


    IMPLICIT NONE

    INTEGER :: idummy,i,j,c,k,i_min,n,h,ii
    INTEGER :: dummy1, loop
    INTEGER ::id_sub_int,id_lu_int,id_lu_ext,id_tc_int,id_soil_int,test
    INTEGER :: idummy2(20),tausch,istate
    integer, pointer :: unique_vals(:)
    REAL :: temp1,temp2, maxthickness
    REAL :: sortier(maxterrain),tauschr
    CHARACTER (LEN=8000) :: cdummy

    !----------------------------------------------
    INTEGER :: PAUL ! for testing GetNumberOfSubstrings(cdummy)
    REAL :: frac_irr_tc
    REAL :: frac_irr_lu
    REAL :: frac_svc_x
    INTEGER :: sub_area_irr(subasin)     ! Für AREA

    INTEGER :: tcid_instance    !(internal) id of TC-instance (unique subbas-LU-TC-combination)
    INTEGER :: soilid            !internal id of current soil
    INTEGER :: sb_counter,lu_counter,tc_counter,svc_counter    ! counters
    INTEGER :: i_lu,id_tc_type, check1, check2, check3 !, svc_id,i_subbas,i_svc,i_soil,i_veg        ! ids of components in work

    INTEGER :: lu_temp(maxsoter)        !auxiliary arrays for reading
    REAL    :: frac_lu_temp(maxsoter)

    !temporary var for reading thickness of soil horizons (mm) (nsoil,maxhori)
    real, allocatable ::   temp_hori_thick(:,:)

    INTEGER :: nbr_svc2(nterrain)    !number of SVCs for each TC-type
    character(len=1000) :: fmtstr    !string for formatting file output

    ! temporary variables for reading landuse characteristics for all
    ! Subbasin-LandscapeUnit-TerrainComponent-combinations
    INTEGER :: luse_subbas(ntcinst),luse_lu(ntcinst),luse_tc(ntcinst)

    INTEGER :: i_soil,i_veg        ! ids of components in work

    REAL    :: wind_t


    !Till: Read routing.dat, which determines which of the given subbasins are to be modelled
    OPEN(11,FILE=pfadp(1:pfadj)// 'River/routing.dat',STATUS='old',IOSTAT=istate)    ! upbasin: MAP ID of upstream sub-basin (MAP IDs);! downbasin: MAP ID of downstream sub-basin (MAP IDs)
    IF (istate/=0) THEN
        write(*,'(a)')'Error: Input file '//pfadp(1:pfadj)// 'River/routing.dat'//' not found, aborting.'
        stop
    END IF


    READ (11,*); READ(11,*)
    istate=0
    upbasin=0
    downbasin=0
    i=0
    DO WHILE (istate==0)
        READ (11,'(A)',IOSTAT=istate)  cdummy
        if (trim(cdummy)=='' .OR. istate/=0) cycle
        i=i+1
        if (i > subasin) then
            write(*,'(A,i0,a,i0,a)')'ERROR: do.dat specifies ',subasin,' subbasins, routing.dat contains more. Correct this.'
            stop
        end if
        READ (cdummy,*,IOSTAT=istate)  idummy, upbasin(i),downbasin(i)
        if (istate /= 0) then
            write(*,'(a)')'ERROR: Format error in routing.dat.'
            stop
        end if
    END DO

    CLOSE (11)
    if (i/=subasin) then        !Till: correct number of subbasins to be modelled
        write(*,'(A,i0,a,i0,a)')'WARNING: do.dat specifies ',subasin,' subbasins, routing.dat contains only ',i,'. Only these will be modelled.'
        subasin=i !ii: check, if unused memory can be freed (see new_real_array3)
        call pause1
    end if

    !check integrity of routing order
    IF (downbasin(subasin)/=999 .AND. downbasin(subasin)/=9999)  THEN
        WRITE (*,'(A,I0,A)') 'ERROR in routing.dat: Last line must contain the outlet (downstream ID must be 999)'
        STOP
    end if
    DO i=1,subasin
        idummy=which1(upbasin(i)==upbasin(i+1:subasin))
        IF (idummy/=0) THEN
            WRITE (*,'(A,I0,A)') 'ERROR in routing.dat: subbasin ', upbasin(i),' cannot drain to more than one basin.'
            STOP
        end if

        if (i==subasin) exit !don't do the next check for the last line
        idummy=which1(upbasin(i)==downbasin(i+1:subasin))
        IF (idummy/=0) THEN
            WRITE (*,'(A,I0,A,I0,A)') 'ERROR in routing.dat: subbasin ', upbasin(i),' receives input from ', upbasin(i+idummy), ". Relocate the former to after the latter in the routing scheme."
            STOP
        end if
    END DO


    call allocate_general()

    !** read subbasin parameters
    OPEN(11,FILE=pfadp(1:pfadj)// 'Hillslope/hymo.dat',STATUS='old',IOSTAT=istate)
    IF (istate/=0) THEN
        write(*,'(a)')'Error: Input file '//pfadp(1:pfadj)// 'Hillslope/hymo.dat'//' not found, aborting.'
        stop
    END IF


    READ(11,*)
    READ(11,*)
    i=1
    istate=0
    id_subbas_extern=upbasin
    id_subbas_intern=0
    id_lu_intern=0
    c=0 !count successfully read subbasins
    h=3 !count lines

    DO WHILE ((istate==0) .AND. (c<subasin))
        cdummy=''
        READ(11,'(a)',IOSTAT=istate) cdummy
        if (istate/=0 .AND. c<subasin-1) then    !(premature) end of file
            write(*,'(a,i0,a)')'ERROR (hymo.dat): At least ',subasin,' lines (#subbasins) expected'
            stop
        end if

        dummy1=GetNumberOfSubstrings(cdummy) !Till: count number of fields/columns
        if (dummy1-3 > 2*maxsoter) then    !too many fields in line
            write(*,'(a,i0,a,i0,a,i0,a)')'ERROR (hymo.dat): line ',h,' contains more (',dummy1,') than the expected 3 + 2 * ',maxsoter,' fields (maxdim.dat).'
            stop
        end if

        lu_temp = 0
        READ(cdummy,*,IOSTAT=istate) n,temp1, k,  (lu_temp(j),j=1,k), (frac_lu_temp(j),j=1,k)
        if (istate/=0) then    !format error
            write(*,'(a,i0)')'ERROR (hymo.dat): Format error in line ',h
            stop
        end if

        i=which1(upbasin==n)
        if (i==0) then    !Till: the current subbasin was not contained in routing.dat, skip it
            WRITE(*,'(a, I0, a, I0, a)') 'Subbasin ',n,' not listed in routing.dat, ignored.'
        else
            area(i)=temp1
            nbr_lu(i)=k
            id_lu_intern(1:k,i)=lu_temp     (1:k)
            frac_lu     (1:k,i)=frac_lu_temp(1:k)
            id_subbas_intern(i)=id_subbas_extern(i)
            c=c+1 !count successfully read subbasins
        end if
        h=h+1 !count lines
    END DO
    CLOSE(11)

    if (c<subasin) then    !less subbasins found than expected
        write(*,'(a,i0,a,i0)')'ERROR: hymo.dat: expected ',subasin,' subbasins, found ',c
        stop
    end if

    if (count(id_subbas_intern(1:subasin)==0) >0) then !check if there are subbasins read from routing.dat that were not found in hymo.dat
        write(*,'(a)')'ERROR: The following subbasins have been listed in routing.dat, but are missing in hymo.dat:'
        do i=1,subasin
            if (id_subbas_intern(i)==0) write(*,*)id_subbas_extern(i) 
        end do    
        stop
    end if



    !** read Landscape units parameters
    OPEN(11,FILE=pfadp(1:pfadj)// 'Hillslope/soter.dat',STATUS='old')
    READ(11,*); READ (11,*)
    h=2

    id_terrain_intern=0
    id_terrain_extern=0

    kfsu(:) = -1. !flag as "not read"
    i=1
    k=0 !count "unused" but specified LUs
    DO WHILE (.TRUE.)
        cdummy=''
        READ(11,'(a)',IOSTAT=istate) cdummy
        if (istate/=0 .AND. len(trim(cdummy))==0) then
            if (h<nsoter+2) then    !premature end of file
                write(*,'(a,i0,a)')'ERROR (soter.dat): At least ',nsoter,' lines (#LUs) expected'
                stop
            end if
            exit
        end if
        h=h+1 !count lines

        dummy1=GetNumberOfSubstrings(cdummy) !Till: count number of fields/columns
        if (dummy1 > 2+maxterrain+8) then    !too many fields in line
            write(*,'(a,i0,a,i0,a,i0,a)')'ERROR (soter.dat): line ',h,' contains more (',dummy1,') than the max expected 2+',maxterrain,'+8 fields (maxdim.dat).'
            stop
        end if

        READ(cdummy,*,IOSTAT=istate) dummy1
        if (istate/=0) then    !format error
            write(*,'(a,i0)')'ERROR (soter.dat): Format error in line ',h
            stop
        end if

        if (.NOT. any(id_lu_intern==dummy1)) then
            !WRITE(*,'(a, I0, a, I0, a)') 'LU ',dummy1,' not listed in hymo.dat, ignored.'
            !k=k+1
            !cycle
        end if

        ! IDs are read into id_terrain_intern - conversion to internal id is done later
        READ(cdummy,*,IOSTAT=istate) id_lu_extern(i),nbrterrain(i), (id_terrain_intern(j,i),j=1,nbrterrain(i)),  &
            kfsu(i),slength(i), meandep(i),maxdep(i),riverbed(i),  &
            gw_flag(i),gw_dist(i),gw_delay(i)
        if (istate/=0) then    !format error
            write(*,'(a,i0)')'ERROR (soter.dat): Format error in line ',h
            stop
        end if

        if (slength(i)<=0) then    !faulty specs
            write(*,'(a,i0,a)')'WARNING (soter.dat): Slope length cannot be <= 0 in line ',h,', corrected to 0.1'
            slength(i)=0.1
        end if

        if (riverbed(i)<=0) then    !faulty specs
            write(*,'(a,i0,a)')'WARNING (soter.dat): Riverbed depth cannot be <= 0 in line ',h,', corrected to 1000'
            riverbed(i)=1000
        end if


        if (gw_flag(i)/=99) then        !Till: gw_dist is only used if gw_flag is 99. To avoid side effects, gw_dist is et to 0
            gw_dist(i)=0.
        end if
        i=i+1
    END DO
    CLOSE(11)
    
    !plausibility checks
    
    do i=1,nsoter-k
        if (nbrterrain(i) < 0.) then
            write(*,'(a,i0,a,i0,a)')'ERROR (soter.dat): line ', i+2,': number of TCs (', nbrterrain(i),') out of range.'
            stop
        end if
        if (kfsu(i) > 1000. .OR. kfsu(i)<0.) then
            write(*,'(a,i0,a,f0.0,a)')'ERROR (soter.dat): line ', i+2,': bedrock hydraulic conductivity (', kfsu(i),') out of range.'
            stop
        end if
        if (slength(i) > 20000. .OR. slength(i)<0.) then
            write(*,'(a,i0,a,f0.0,a)')'ERROR (soter.dat): line ', i+2,': LU length (', slength(i),') out of range.'
            stop
        end if
        if ((meandep(i) > 100000. .OR. meandep(i)<0.) .AND. meandep(i)/=-1.) then
            write(*,'(a,i0,a,f0.0,a)')'ERROR (soter.dat): line ', i+2,': mean soil depth (', meandep(i),') out of range.'
            stop
        end if
        if ((maxdep(i) > 100000. .OR. maxdep(i)<0.) .AND. maxdep(i)/=-1.) then
            write(*,'(a,i0,a,f0.0,a)')'ERROR (soter.dat): line ', i+2,': max soil depth (', maxdep(i),') out of range.'
            stop
        end if
        
        if (riverbed(i) > 100000. .OR. riverbed(i)<0.) then
            write(*,'(a,i0,a,f0.0,a)')'ERROR (soter.dat): line ', i+2,': depth to riverbed (', riverbed(i),') out of range.'
            stop
        end if
        
        if (gw_delay(i) > 100000. .OR. gw_delay(i)<0.) then
            write(*,'(a,i0,a,f0.0,a)')'ERROR (soter.dat): line ', i+2,': GW-delay (', gw_delay(i),') out of range.'
            stop
        end if
    end do
    
    gw_delay=max(1.,gw_delay)    !Till: gw_delay cannot be less than 1 (=all GW leaves LU immediately)

    !** read terrain component parameters
    OPEN(11,FILE=pfadp(1:pfadj)// 'Hillslope/terrain.dat',STATUS='old')
    READ(11,*);READ(11,*)
    READ(11,'(a)') cdummy
    dummy1=GetNumberOfSubstrings(cdummy) !Till: count number of fields/columns
    rewind(11)
    READ(11,*);READ(11,*)

    !ii: reading like this may be problematic with input files with more or less than the xpected nterrain lines
    if (dummy1>4) then !more than 4 columns
        allocate(beta_fac_tc(nterrain))
        beta_fac_tc=1.        !default value, no beta correction
    end if
    if (dummy1>5) then !6 columns, containing beta correction factor and TC-wise SDR
        allocate(sdr_tc(nterrain))
        sdr_tc=1.        !default value, no SDR correction
    end if

    if (dummy1 == 8) then !8 columns, containing beta correction factor, TC-wise SDR and concentration factors
        allocate(frac_diff2conc(nterrain))
        allocate(frac_conc2diff(nterrain))
        frac_diff2conc = 0. !default values, same as dolattc = doalllattc = TRUE
        frac_conc2diff = 0.
    end if


    if (dummy1==4) then !4 columns, old version
        READ(11,*,IOSTAT=istate) (id_terrain_extern(i),fracterrain(i), slope(i),posterrain(i),i=1,nterrain)
    elseif (dummy1==5) then !5 columns, containing beta correction factor
        READ(11,*,IOSTAT=istate) (id_terrain_extern(i),fracterrain(i), slope(i),posterrain(i),beta_fac_tc(i),i=1,nterrain)
    elseif (dummy1==6) then !6 columns, containing beta correction factor and TC-wise SDR
        READ(11,*,IOSTAT=istate) (id_terrain_extern(i),fracterrain(i), slope(i),posterrain(i),beta_fac_tc(i),sdr_tc(i),i=1,nterrain)
    elseif (dummy1==8) then !8 columns, containing beta correction factor, TC-wise SDR and concentration factors
        READ(11,*,IOSTAT=istate) (id_terrain_extern(i),fracterrain(i), slope(i),posterrain(i),beta_fac_tc(i),sdr_tc(i),frac_diff2conc(i),frac_conc2diff(i),i=1,nterrain)
    end if

    if (istate/=0) then    !format error
        write(*,'(a,i0,a)')'ERROR (terrain.dat): Format error in line ',i+2,". Check formatting and number of TCs specified in do.dat"
        stop
    end if

    READ(11,'(a)',IOSTAT=istate) cdummy
    if (istate==0 .AND. cdummy/="") then    !excess lines
        write(*,'(a,i0,a)')'WARNING (terrain.dat): Found more than the expected ',nterrain," data lines. File truncated. Check number of TCs specified in do.dat"
    end if

    !free memory containing correction factors, if sediment is disabled or all set to 1
    if (allocated(beta_fac_tc)) then
        if ((.NOT. dosediment) .OR. all(beta_fac_tc==1.)) deallocate(beta_fac_tc)
    end if
    if ( allocated(sdr_tc))     then
        if ((.NOT. dosediment) .OR. all(sdr_tc==1.))       deallocate(sdr_tc)
    end if

    CLOSE(11)


    !** read landuse characteristics for all terrain components
    !   (soil-vegetation components)
    OPEN(11,FILE=pfadp(1:pfadj)//  &
        'Hillslope/soil_vegetation.dat',STATUS='old')
    READ(11,*); READ(11,*);READ(11,*)
    h=4
    luse_subbas=0    !initial values
    luse_lu=0
    luse_tc=0
    rocky=0.
    nbr_svc=0
    frac_svc=0.
    nbr_svc=0
    id_soil_intern=0
    id_veg_intern =0

    i=1 !count instances that have been read
    DO WHILE (.TRUE.) !read till end of file
        READ(11,'(a)',IOSTAT=istate) cdummy
        if (k /= 0 .or. istate /=0) then
            if ((k == 7)  ) then    !inconsistent line scheme
                write(*,'(a,i0,a,i0,a)')'ERROR (soil_vegetation.dat): Inconsistent IDs in block lines ', h-3,' - ', h-1,'. Check IDs and 3 headerlines.'
                stop
            endif    
            if ((h-4 < ntcinst)  ) then    !premature end of file
                write(*,'(a,i0,a)')'ERROR (soil_vegetation.dat): ', ntcinst,' x 3 lines (#TC-LU-SUBBAS-combinations) + 3 headerlines expected'
                stop
            else
                if (i-1/=ntcinst) then    !less entities read than expected
                    write(*,'(a,i0,a,i0,a)')'WARNING (soil_vegetation.dat): ',i-1,' instead of the expected ',ntcinst,' TCs read.'
                    ntcinst=i-1 !correct value
                    ! resize arrays
                    nbr_svc =>new_int_array1(nbr_svc, ntcinst, 1)
                    rocky=>     new_real_array1(rocky, ntcinst,1)
                    laitc=>     new_real_array1(laitc, ntcinst,1)
                    aettc=>     new_real_array1(aettc, ntcinst,1)
                    soilettc=>  new_real_array1(soilettc, ntcinst,1)
                    intctc=>    new_real_array1(intctc, ntcinst,1)
                    horttc=>    new_real_array1(horttc, ntcinst,1)
                    gwrtc=>     new_real_array1(gwrtc, ntcinst,1)
                    deepgwrtc=> new_real_array1(deepgwrtc, ntcinst,1)

                    id_soil_intern=> new_int_array2(id_soil_intern, ntcinst,2)
                    id_veg_intern => new_int_array2(id_veg_intern, ntcinst,2)
                    svcrooth =>      new_int_array2(svcrooth, ntcinst,1)
                    svcbedr  =>      new_int_array2(svcbedr, ntcinst,1)

                    frac_svc  => new_real_array2(frac_svc, ntcinst,2)
                    soilwater => new_real_array2(soilwater,ntcinst,2)

                    intercept=> new_real_array2(intercept, ntcinst,1)
                    frac_sat => new_real_array2(frac_sat, ntcinst,1)

                    horiz_thickness=> new_real_array3(horiz_thickness,ntcinst,1)
                    pwpsc          => new_real_array3(pwpsc          ,ntcinst,1)
                    horiths        => new_real_array3(horiths        ,ntcinst,1)
                    horithact      => new_real_array3(horithact      ,ntcinst,1)
                end if
                exit !enough lines read, abort loop
            end if
        end if

        if (i-1 > ntcinst  ) then    !more lines than expected
            write(*,'(a,i0,a)')'ERROR (soil_vegetation.dat): ', ntcinst,' x 3 lines (#TC-LU-SUBBAS-combinations) expected, found more'
            stop
        end if


        dummy1=GetNumberOfSubstrings(cdummy)-5 !Till: count number of fields (ie SVCs) specified for this combination
        if (dummy1 > maxsoil) then
            write (*,'(a,i0,a,i0,a)')'ERROR: Line ',h,' in soil_vegetation.dat contains ', dummy1,' soil types - more than specified in maxdim.dat or assumed by default.'
            stop
        end if
        k=0
        READ(cdummy,*,IOSTAT=istate) luse_subbas(i) !read subbasin ID only
        h=h+1 !count lines
        if (which1(luse_subbas(i) == id_subbas_extern) == 0) then
            !write(*,'(a,i0)')'WARNING (soil_vegetation.dat): Unknown subbasin ',luse_subbas(i),' in line ',h,', ignored.'
            READ(11,'(a)',IOSTAT=istate) cdummy !skip next two lines
            READ(11,'(a)',IOSTAT=istate) cdummy
            h=h+2 !count lines
        else
            READ(cdummy,*,IOSTAT=istate) luse_subbas(i),luse_lu(i),luse_tc(i),  &
                rocky(i),nbr_svc(i),(id_soil_intern(c,i),c=1,dummy1)    !!ii: nbr_svc, rocky, id_soil_intern sind feststehende Parameter für einen TC-Typ und sollten nur einmal pro TC-typ gespeichert werden
            k=k+istate
            READ(11,*,IOSTAT=istate) check1, check2, check3,  &
                rocky(i),nbr_svc(i),(id_veg_intern(c,i),c=1,dummy1)
            k=k+istate
            if (check1 /= luse_subbas(i) .OR. check2 /= luse_lu(i) .OR. check3 /= luse_tc(i)) k = 7 !flag an error if there is an inconsistency of IDs in first three cols
            
            h=h+1 !count lines
   
            READ(11,*,IOSTAT=istate) check1, check2, check3, &
                rocky(i),nbr_svc(i),(frac_svc(c,i),c=1,dummy1)
            k=k+istate
            if (check1 /= luse_subbas(i) .OR. check2 /= luse_lu(i) .OR. check3 /= luse_tc(i)) k = 7 !flag an error if there is an inconsistency of IDs in first three cols
            
            h=h+1 !count lines
            i=i+1 !count instances that have been read
        end if

        if (istate/=0) then    !premature end of file
            write(*,'(a,i0)')'ERROR (soil_vegetation.dat): Format error or unexpected end in line ',h
            stop
        end if

    END DO

    cdummy=''
    READ(11,'(a)', IOSTAT=istate) cdummy
    if (trim(cdummy)/='') then
        write(*,'(a,i0,a)')'ERROR: in soil_vegetation.dat: more than the expected 3 x ',ntcinst,', lines found.'
        stop
    end if

    CLOSE(11)



    !** read soil properties
    allocate (temp_hori_thick(nsoil,maxhori))    !ii: gleich in hori_thick einlesen
    temp_hori_thick(:,:)=0.

    OPEN(11,FILE=pfadp(1:pfadj)// 'Hillslope/soil.dat',STATUS='old')
    READ(11,*)
    READ(11,*)
    READ(11,'(a)') cdummy
    c=MOD(GetNumberOfSubstrings(cdummy),13) !Till: count number of fields to decide if this is an older input file version
    if (c==3) then
        write(*,'(a)')'WARNING: old version of soil.dat without alluvial_flag. Soils 1-3 assumed to be alluvial.'
    elseif (c/=4) then
        write(*,'(a)')'ERROR: unknown version of soil.dat'
        stop
    end if

    rewind(11)
    READ(11,*)
    READ(11,*)
    h=3

    alluvial_flag=0


    shrink = 0.
    DO j=1,nsoil


        READ(11,'(a)',IOSTAT=istate) cdummy
        if ((istate/=0).OR. (trim(cdummy)=='')) then
            write(*,'(a,i0,a,i0,a)')'ERROR: in soil.dat: Expected ',nsoil,', found ',j-1,' soil types.'
            stop
        end if

        dummy1=GetNumberOfSubstrings(cdummy) !Till: count number of fields

        if (dummy1 > 2+maxhori*13+1+(c-3)) then
            write (*,'(A,i0,a,i0,a)')'ERROR: Line ',h,' in soil.dat contains ',(dummy1-3-(c-3))/13,' horizons - more than specified in maxdim.dat or assumed by default.'
            stop
        end if

        if (c==3) then    !old format
            READ(cdummy,*,IOSTAT=istate) id_soil_extern(j),nbrhori(j), (thetar(j,k),soilpwp(j,k),  &
                soilfc(j,k),soilfc63(j,k), soilnfk(j,k),thetas(j,k),  &
                temp_hori_thick(j,k), k_sat(j,k),saug(j,k),poresz(j,k),  &
                bubble(j,k),coarse(j,k),shrink(j,k),k=1,nbrhori(j)), bedrock(j)
            if (id_soil_extern(j)==1 .OR. id_soil_extern(j)==2 .OR. id_soil_extern(j)==3) then
                alluvial_flag(j)=1    !setting that was implied in old code
            end if
        else
            READ(cdummy,*,IOSTAT=istate) id_soil_extern(j),nbrhori(j), (thetar(j,k),soilpwp(j,k),  &
                soilfc(j,k),soilfc63(j,k), soilnfk(j,k),thetas(j,k),  &
                temp_hori_thick(j,k), k_sat(j,k),saug(j,k),poresz(j,k),  &
                bubble(j,k),coarse(j,k),shrink(j,k),k=1,nbrhori(j)), bedrock(j), alluvial_flag(j)
        end if

        if (istate/=0) then
            write(*,'(a, i0)')'ERROR (soil.dat): Format error in , line ',h
            stop
        end if
        h=h+1

        DO k=1,nbrhori(j)
            if (coarse(j,k)>1.) then
                write(*,'(a,i0)')'ERROR: coarse fraction cannot be larger than 1. Check soil.dat in line ', j+2
                stop
            end if
        END DO

    END DO

    shrink = 0.    !Till shrinkage currently disabled (buggy, "macro" can become huge and lead to crashes in soilwat 2.2a)
    CLOSE(11)





    !** read vegetation characteristics
    OPEN(11,FILE=pfadp(1:pfadj)// 'Hillslope/vegetation.dat',STATUS='old')
    READ(11,*); READ (11,*)
    j=1
    h=3 !line counter
    id_veg_extern=0
    resist=0
    wstressmin=0
    wstressmax=0
    height=0
    rootdep=0
    lai=0
    alb=0
    DO WHILE (.TRUE.)

        cdummy=''
        READ(11,'(a)',IOSTAT=istate) cdummy

        if (istate == -1 ) then
!            if ((trim(cdummy)=='') .and. j < nveg) then
!                write(*,'(a,i0,a,i0)')'ERROR: vegetation.dat, expected ', nveg, ' valid entries, found ',j
!                stop
!            else
                exit
!            end if
        end if
        
         if ((trim(cdummy)/='') .and. j > nveg) then
            write(*,'(a,i0,a)')'WARNING: found more than the expected ', nveg, ' entries in vegetation.dat, ignored.'
            exit
         END if

        READ(cdummy,*, IOSTAT=istate) id_veg_extern(j),resist(j),wstressmin(j),wstressmax(j),  &
            (height(j,i),i=1,4),(rootdep(j,i),i=1,4), (lai(j,i),i=1,4),(alb(j,i),i=1,4)

        if (istate /= 0) then
            write(*,'(a,i0)')'ERROR: vegetation.dat, format error in line ',h
            stop
        end if

        h=h+1


        if (.not. (any (id_veg_intern == id_veg_extern(j)))) then !check if the vegetation-ID occurred in soil_vegetation.dat
            !if (  size(pack(id_veg_intern, id_veg_extern(j) == id_veg_intern(:,:))) == 0  ) then
            write(*,'(a,i0,a,i0,a)')'WARNING: unused vegetation-id ',id_veg_extern(j),' in vegetation.dat, line ',h-1
            cycle !ii: we should not cycle this here, otherwise later errors in reading SVCs may occur (?)
        end if

        if (wstressmin(j) >= wstressmax(j)) then
            write(*,'(a,i0,a)')'ERROR: vegetation.dat, line ', h,': wstressmin must be < wstressmax'
            stop
        end if

        j=j+1

    END DO
    CLOSE(11)
    
    
    if (j < nveg) then !less 
        unique_vals => unique2d(id_veg_intern) !get unique IDs loaded from soil_vegetation.dat
        do i=1, size(unique_vals)
           if (unique_vals(i)==0) cycle
            if (.not. (any (id_veg_extern== unique_vals(i)))) then !check if the all vegetation-ID from soil_vegetation.dat have been found
              write(*,'(a,i0,a,i0,a,i0)')'ERROR: vegetation.dat, expected ', nveg, ' valid entries, found ',j,', e.g. missing class ',unique_vals(i)
              stop
           end if    
        end do
        nveg = j !decrease actual number of vegetation classes
        id_veg_extern => new_int_array1(id_veg_extern, nveg, 1) !shrink array size
        rootdep    => new_real_array2(rootdep, nveg, 1) 
        lai        => new_real_array2(lai,     nveg, 1)  
        alb        => new_real_array2(alb,     nveg, 1)   
        height     => new_real_array2(height,  nveg, 1)  
        
    end if
    
    
    rootdep(:,:)=rootdep(:,:)*1000.
    !!** read key points for temporal distribution of vegetation
    !!   characteristics within year (end/begin of rainy season)
    period => seasonality_array2('rainy_season.dat')
    CALL check_seasonality(period, "rainy_season.dat", id_veg_extern)  !check completeness of seasonality file


    if (doloadstate .OR. dosavestate .OR. dosediment) then
        seasonality_k     => seasonality_array2('k_seasons.dat')    !read seasonality of K-factor
        seasonality_c     => seasonality_array2('c_seasons.dat')    !read seasonality of C-factor

        seasonality_p     => seasonality_array2('p_seasons.dat'     )    !read seasonality of P-factor
        seasonality_coarse=> seasonality_array2('coarse_seasons.dat')    !read seasonality of coarse fraction factor
        seasonality_n     => seasonality_array2('n_seasons.dat'     )    !read seasonality of Manning's n


        !** read SVC information (numbering scheme, erosion properties)
        OPEN(11,FILE=pfadp(1:pfadj)// 'Hillslope/svc.dat', IOSTAT=istate,STATUS='old')
        IF (istate/=0) THEN
            write(*,'(a)') "ERROR: "//pfadp(1:pfadj)// 'Hillslope/svc.dat could not be opened. Supply this file or disable sediment modelling or loading/saving states. Aborting.'
            stop
        END IF

        READ(11,*)  !----2 Zeilen überspringen, dann wird die Anzahl der SVC ermittelt, also Im Prinzip die Anzahl der Zeilen
        READ(11,*)

        i=0            !for counting SVC records
        DO WHILE (.TRUE.)
            READ(11,'(a)', IOSTAT=k) cdummy
            IF (k/=0) THEN    !no further line
                nsvc=i        !store total number of SVCs that have been read
                exit        !exit loop
            END IF
            i=i+1            !increase record counter
        END DO

        k=0    !count number of columns to be expected in file

         !------------------------------------------------

        READ(11,'(a)', IOSTAT=istate) cdummy

        PAUL = GetNumberOfSubstrings(cdummy)  ! Feststellen der Spalten von Svc.dat


        !--------------------------------------------------------------------------------------------

        allocate(id_svc_extern(nsvc))    !allocate memory for SVC-IDs to be read


        j=(tstop-tstart)+1        !number of years to simulate

        if (SIZE(seasonality_k,dim=2) == 1) then
            allocate (svc_k_fac(nsvc, 1))
            svc_k_fac_day=>svc_k_fac(:,1)        !no dynamics, daily values remain static as read from file
        else
            allocate (svc_k_fac(nsvc, 4))
            allocate(svc_k_fac_day(nsvc))    !allocate memory for temporal dynamics (current day and subbasin)
        end if
        k=k+SIZE(svc_k_fac,dim=2)
        svc_k_fac=0.

        if (SIZE(seasonality_c,dim=2) == 1) then
            allocate (svc_c_fac(nsvc, 1))
            svc_c_fac_day=>svc_c_fac(:,1)        !no dynamics, daily values remain static as read from file
        else
            allocate (svc_c_fac(nsvc, 4))
            allocate(svc_c_fac_day(nsvc))    !allocate memory for temporal dynamics (current day and subbasin)
        end if
        k=k+SIZE(svc_c_fac,dim=2)
        svc_c_fac=0.

        if (SIZE(seasonality_p,dim=2) == 1) then
            allocate (svc_p_fac(nsvc, 1))
            svc_p_fac_day=>svc_p_fac(:,1)        !no dynamics, daily values remain static as read from file
        else
            allocate (svc_p_fac(nsvc, 4))
            allocate(svc_p_fac_day(nsvc))    !allocate memory for temporal dynamics (current day and subbasin)
        end if
        k=k+SIZE(svc_p_fac,dim=2)
        svc_p_fac=0.

        if (SIZE(seasonality_coarse,dim=2) == 1) then
            allocate (svc_coarse_fac(nsvc, 1))
            svc_coarse_fac_day=>svc_coarse_fac(:,1)        !no dynamics, daily values remain static as read from file
        else
            allocate (svc_coarse_fac(nsvc, 4))
            allocate(svc_coarse_fac_day(nsvc))    !allocate memory for temporal dynamics (current day and subbasin)
        end if
        k=k+SIZE(svc_coarse_fac,dim=2)
        svc_coarse_fac=0.

        if (SIZE(seasonality_n,dim=2) == 1) then
            allocate (svc_n(nsvc, 1))
            svc_n_day=>svc_n(:,1)        !no dynamics, daily values remain static as read from file
        else
            allocate (svc_n(nsvc, 4))
            allocate(svc_n_day(nsvc))    !allocate memory for temporal dynamics (current day and subbasin)
        end if
        k=k+SIZE(svc_n,dim=2)
        svc_n=0.

        allocate(svc_soil_veg(nsvc,2))        !allocate memory for SVCs (soil-vegetation relation) to be read

        rewind(11)
        READ(11,*)
        READ(11,*)

        k=k+3    !number of fields plus columns for 3 IDs
        h=2 !for counting lines
        i=1
        write(fmtstr,*)'(3i, ', SIZE(seasonality_k,dim=2),'F, ',SIZE(seasonality_c,dim=2),'F, ',SIZE(seasonality_p,dim=2),'F, ',SIZE(seasonality_coarse,dim=2),'F, ',SIZE(seasonality_n,dim=2),'F)'    !generate format string according to number of columns to be expected

        !-------------------------------------------------
        IF (PAUL /= k) THEN
        allocate (svc_irr(nsvc))  !svc_irr vektor anlegen. Wird in erosion.dat deklariert -> ändern in hymo_h
        svc_irr = 0
        allocate (frac_irr_sub(subasin))
        frac_irr_sub = 0.

        END IF
        !-------------------------------------------------------

        DO WHILE (.TRUE.)
            READ(11,'(a)', IOSTAT=istate) cdummy
            h=h+1
            IF (istate==-1) THEN    !no further line
                nsvc=i-1    !correct total number of SVCs that have been read
                !ii: check, if unused memory can be freed (see new_real_array3)
                exit        !exit loop
            END IF
            if (trim(cdummy)=='') cycle    !skip blank lines

            !------------------------------------ Falls es mehr als 8 Spalten gibt, lies svc_irr ein
            IF (PAUL == k) THEN
                READ(cdummy,*, IOSTAT=istate) id_svc_extern(i), j, n, svc_k_fac(i,:), svc_c_fac(i,:), svc_p_fac(i,:), svc_coarse_fac(i,:), svc_n(i,:)
            ELSE
                READ(cdummy,*, IOSTAT=istate) id_svc_extern(i), j, n, svc_k_fac(i,:), svc_c_fac(i,:), svc_p_fac(i,:), svc_coarse_fac(i,:), svc_n(i,:), svc_irr(i)
            END IF

            IF (istate /= 0 .OR. (PAUL /=k .AND. PAUL /= k + 1)) THEN    !format error
                write(*,'(a,i0,a)')'ERROR: svc.dat, line ',h,': unexpected number of fields. Check *_seasons.dat'
                stop
            END IF

            svc_soil_veg(i,1)=id_ext2int(j, id_soil_extern)    !store internal soil id of this SVC
            IF (svc_soil_veg(i,1)==-1) THEN    !specified soil not found
                write(*,'(a,i0,a,i0,a)')'ERROR in svc.dat, line ',h,': could not find soil ', j,' in soil.dat'
                cycle
                !stop
            END IF

            svc_soil_veg(i,2)=id_ext2int(n, id_veg_extern)     !store internal vegetation id of this SVC
            IF (svc_soil_veg(i,2)==-1) THEN    !specified vegetation not found
                write(*,'(a,i0,a,i0,a)')'ERROR in svc.dat, line ',h,': could not find vegetation ', n,' in vegetation.dat'
                cycle
                !stop
            END IF
            i=i+1
        END DO
        CLOSE(11)
        
		if (SIZE(seasonality_k,dim=2) /= 1) then
            CALL check_seasonality(seasonality_k, "k_seasons.dat", id_svc_extern)  !check completeness of seasonality file
        end if
        if (SIZE(seasonality_c,dim=2) /= 1) then
            CALL check_seasonality(seasonality_c, "c_seasons.dat", id_svc_extern)  !check completeness of seasonality file
        end if
        if (SIZE(seasonality_p,dim=2) /= 1) then
            CALL check_seasonality(seasonality_p, "p_seasons.dat", id_svc_extern)  !check completeness of seasonality file
        end if
        if (SIZE(seasonality_coarse,dim=2) /= 1) then
            CALL check_seasonality(seasonality_coarse, "coarse_seasons.dat", id_svc_extern)  !check completeness of seasonality file
        end if
        if (SIZE(seasonality_n,dim=2) /= 1) then
            CALL check_seasonality(seasonality_n, "n_seasons.dat", id_svc_extern)  !check completeness of seasonality file
        end if
    end if


    CALL check_seasonality_superfluous(seasonality_k,      "k_seasons.dat",      id_svc_extern)  !check for obsolete entries in seasonality file
    CALL check_seasonality_superfluous(seasonality_c,      "c_seasons.dat",      id_svc_extern)  !check for obsolete entries in seasonality file
    CALL check_seasonality_superfluous(seasonality_p,      "p_seasons.dat",      id_svc_extern)  !check for obsolete entries in seasonality file
    CALL check_seasonality_superfluous(seasonality_coarse, "coarse_seasons.dat", id_svc_extern)  !check for obsolete entries in seasonality file
    CALL check_seasonality_superfluous(seasonality_n,      "n_seasons.dat",      id_svc_extern)  !check for obsolete entries in seasonality file


    if (dosediment)    then                !Till: if erosion is to be modelled, read these files

        !** read particle size classes
        OPEN(11,FILE=pfadp(1:pfadj)// 'part_class.dat', IOSTAT=istate,STATUS='old')
        IF (istate/=0) THEN                    !part_class.dat not found
            write(*,'(a)')pfadp(1:pfadj)// 'part_class.dat could not be opened. Sediment treated as one size class.'
            n_sed_class=1                !treat one particle size class only
            allocate(upper_limit(n_sed_class))    !allocate memory for particle size classes to be read
            allocate(particle_classes(n_sed_class))    !allocate memory for particle size classes to be read
            upper_limit(1)=2.        !set upper limit of sediment that is considered to 2 mm
            particle_classes(1)=sqrt(2.*0.064)        !set mean particle size to something in the fine sand range
        else                            !part_class.dat sucessfully opened
            READ(11,*)
            READ(11,*)
            i=0                            !for counting number of particle size classes
            DO WHILE (.TRUE.)
                READ(11,'(a)', IOSTAT=istate) cdummy
                IF (istate/=0) THEN            !no further line
                    n_sed_class=i        !store total number of particle size classes that have been read
                    exit                !exit loop
                END IF
                i=i+1                    !increase record counter
            END DO

            allocate(upper_limit(n_sed_class))    !allocate memory for particle size classes to be read
            allocate(particle_classes(n_sed_class))    !allocate memory for particle size classes to be read
            rewind(11)
            READ(11,*)
            READ(11,*)

            DO i=1,n_sed_class
                READ(11,*,  IOSTAT=istate) j, temp1
                IF (istate/=0) THEN            !no further line
                    n_sed_class=i-1        !correct total number of particle size classes that have been read
                    !ii: check, if unused memory can be freed (see new_real_array3)
                    exit                !exit loop
                END IF
                upper_limit(i)=temp1        !store upper limit of particle size class
            END DO
            CLOSE(11)

            !convert class limits to "mean" diameter
            DO i=2,n_sed_class
                particle_classes(i)=sqrt(upper_limit(i)*upper_limit(i-1))        !use geometric mean as a mean diameter for each class
            END DO
            particle_classes(1)=upper_limit(1)/sqrt(2.)    !use arithmetic mean for first class (geometric won't work here)

        END IF

        !** read particle size information for uppermost soil horizons
        k=0
        IF (n_sed_class > 1) THEN    !if several particle classes are to be treated, try to open file
            OPEN(11,FILE=pfadp(1:pfadj)// 'Hillslope/soil_particles.dat', IOSTAT=istate,STATUS='old')
        END IF
        IF (istate/=0) THEN    !file could not be opened
            write(*,'(a)')pfadp(1:pfadj)// 'soil_particles.dat could not be opened. Evenly distributed particle-size fractions assumed.'
        END IF
        IF ((istate/=0).OR.(n_sed_class == 1)) THEN    !only 1 size class or file could not be opened
            allocate(soil_particles(nsoil, n_sed_class))    !allocate memory for particle-size distributions of soil surfaces
            soil_particles=1./n_sed_class                    !evenly distribute particle classes for all soils
        ELSE
            READ(11,*)
            READ(11,*)

            allocate(soil_particles(nsoil, n_sed_class))    !allocate memory for each soil and particle class
            soil_particles(:,:)=-1.                            !to detect missing values later

            loop=3                                            !line counter
            READ(11,*,  IOSTAT=istate) i, ii, temp1
            DO WHILE (istate==0)
                j = id_ext2int(i, id_soil_extern)            !convert external to internal ID
                IF (j==-1) THEN                                        !found unknown subbasin id
                    WRITE(*,'(a, I0, a, I0, a)') 'soil_particles.dat, line ', loop, ': soil-ID ', i, ' not found in soil.dat, ignored.'
                ELSEIF (ii<1 .OR. ii>n_sed_class) THEN
                    WRITE(*,'(a, I0, a, I0, a)') 'soil_particles.dat, line ', loop, ': particle-class ', ii,&
                        ' not defined in part_class.dat, ignored.'
                ELSE
                    soil_particles(j,ii)=temp1                !store read values in proper position
                ENDIF

                READ(11,*,  IOSTAT=istate) i, ii, temp1
                loop=loop+1                                    !increase line counter
            END DO
            CLOSE(11)
            DO j=1,n_sed_class  !check completeness
                DO i=1,nsoil
                    IF (soil_particles(i,j)==-1.) then        !found soil/particle class for that no fraction has been read
                        WRITE(*,'(a, I3, a, I3)') 'ERROR: soil_particles.dat: Soil ', id_soil_extern(i), ' is lacking particle-size data for class ', j
                        STOP
                    END IF
                END DO
            END DO
        END IF







    END IF

    IF (dosediment .OR. doirrigation) THEN
        !** read SVC relations towards TCs  if do sediment or do irrigation
        OPEN(11,FILE=pfadp(1:pfadj)// 'Hillslope/svc_in_tc.dat', IOSTAT=istate,STATUS='old')
        IF (istate/=0) THEN
            write(*,'(a)') "ERROR: ", pfadp(1:pfadj)// 'Hillslope/svc_in_tc.dat could not be opened. Aborting.'
            stop
        END IF
        READ(11,*)
        READ(11,*)

        nbr_svc2=0
        ii=-1
        DO WHILE (.TRUE.)        !count SVCs in each TC
            READ(11,*,  IOSTAT=istate) i, j, temp1
            IF (istate/=0) THEN    !no further line
                exit        !exit loop
            END IF
            k=id_ext2int(i, id_terrain_extern)    !convert to internal TC-type id
            IF (k==-1) THEN    !ID not found
                if (i /=ii) write(*,'(A,i0,A)')'WARNING: TC-ID ', i,' in svc_in_tc.dat not found in terrain.dat. Ignored.'
                ii=i !prevent repeated error messages of same TC
                cycle
            END IF
            nbr_svc2(k)=nbr_svc2(k)+1        !increase number SVCs for this TC
        END DO

        allocate(tc_contains_svc2(nterrain))                    !allocate space for each TC-type
        DO i=1,nterrain
            allocate(tc_contains_svc2(i)%p(nbr_svc2(i)))    !allocate necessary number of SVCs for each TC-type
        END DO

        rewind(11)
        READ(11,*)
        READ(11,*)

        nbr_svc2=0    !re-used as a pointer to the last SVC-entry of each TC
        DO WHILE (.TRUE.)
            READ(11,*,  IOSTAT=istate) i, j, temp1
            IF (istate/=0) THEN    !no further line
                exit        !exit loop
            END IF
            k=id_ext2int(i, id_terrain_extern)    !convert to internal id
            IF (k==-1) cycle!ID not found

            n=id_ext2int(j, id_svc_extern)    !convert to internal id

            nbr_svc2(k)=nbr_svc2(k)+1
            tc_contains_svc2(k)%p(nbr_svc2(k))%svc_id=n    !store internal id
            tc_contains_svc2(k)%p(nbr_svc2(k))%fraction=temp1    !store SVC-fraction

        END DO
        CLOSE(11)
    END IF ! (dosediment .OR. doirrigation)



    ! ** read cell-based scaling factor, see Guentner (2002), p. 67
    !Till: the scaling factor described there refers to Kfkorr,
    ! but this here affects both Kfkorr and interception capacity
    IF (doscale) THEN
        OPEN(11,FILE=pfadp(1:pfadj)// 'Others/scaling_factor.dat', IOSTAT=istate, STATUS='old')
        IF (istate==0) THEN                    !scaling_factor.dat found
            READ(11,*)
            DO WHILE (.TRUE.)
                READ(11,'(a)',IOSTAT=istate)cdummy    !try to read next line
                if (istate/=0) exit

                READ(cdummy,*)dummy1,temp1
                if (dummy1==-1) then
                    kfkorrc(:) = temp1 !universal calibration factor for all subbasins
                    intcfc(:)  = intcf/(0.340+0.647*kfkorrc(:))  !Till: as in the original code: interception capacity is also modified. I dunno why.
                    cycle                         !go to next line
                end if
                i=id_ext2int(dummy1,id_subbas_extern)    !convert external to internal ID
                if (i==-1) then
                    write(*,'(a,i0,a)')'Unknown subbasin-ID ',dummy1,' in scaling_factor.dat, ignored.'
                    cycle
                end if
                kfkorrc(i) = temp1                !modify kfkorrc of specified subbasin
                intcfc(i)  =intcf/(0.340+0.647*kfkorrc(i)) !Till: as in the original code: interception capacity is also modified. I dunno why.
            END DO
            CLOSE(11)
        END IF
    ELSE
        kfkorrc(:)=1.
        intcfc(:)=intcf
    END IF


    !Till: read calibration factors
    OPEN(11,FILE=pfadp(1:pfadj)// 'Others/calibration.dat', IOSTAT=istate, STATUS='old')
    IF (istate==0) THEN                    !calibration.dat found
        READ(11,*)
        DO WHILE (.TRUE.)
            READ(11,'(a)',IOSTAT=istate)cdummy    !try to read next line
            if (istate/=0) exit

            READ(cdummy,*)dummy1,temp1
            if (dummy1==-1) then
                k_sat(:,:)=k_sat(:,:)*temp1     !universal calibration factor for all soils and horizons
                cycle                         !go to next line
            end if
            i=id_ext2int(dummy1,id_soil_extern)    !convert external to internal ID
            if (i==-1) then
                write(*,'(a,i0,a)')'ERROR: Unknown soil-ID ',dummy1,' in calibration.dat'
                stop
            end if
            k_sat(i,:)=k_sat(i,:)*temp1                !modify Ksat of specified soil

        END DO
        CLOSE(11)
    END IF


    !** make arrays which point to the (model internal) subscript of
    !   soter/terrain/soil/vegetation unit if the (external) ID of this unit
    !   is given.
    !   this transformation essentially corresponds to the exchange of the
    !   external ID of a unit to the subscript it has in the model according
    !   to the order when reading the .dat file.

    !   IDs of subbasins
    DO i=1,subasin
        n=1
        DO WHILE (id_subbas_extern(n) /= id_subbas_intern(i))    !search until the current external ID has been found...
            n=n+1
            if (n>subasin) then
                write(*,*)'ERROR: Subbasin-ID mismatch.'
                stop
            end if
        END DO
        id_subbas_intern(i)=n
    END DO

    DO i=1,subasin
        j=1
        DO WHILE (j <= nbr_lu(i))
            n=1
            DO WHILE ((n<=nsoter) .AND. (id_lu_extern(n) /= id_lu_intern(j,i)))     !search until the current external ID has been found...
                n=n+1
            END DO

            IF (n > nsoter) THEN    !reference to undefined LU
                WRITE(*,'(a, I0,a, I0, a)') 'ERROR: LU ', id_lu_intern(j,i), ' (Subbasin ', id_subbas_extern(i), ') not found in soter.dat.'
                STOP
            END IF
            id_lu_intern(j,i)=n                                !replace external ID with internal ID
            j=j+1
        END DO
    END DO


    !   IDs of terrain components: convert external to internal
    DO i=1,nsoter
        DO j=1, nbrterrain(i)
            n=which1(id_terrain_extern == id_terrain_intern(j,i))
            IF (n == 0) THEN    !reference to undefined TC
                WRITE(*,'(a, I0,a, I0, a)') 'WARNING: TC ', id_terrain_intern(j,i), ' (Subbasin ', id_subbas_extern(i), ') in soter.dat not found in terrain.dat. Toposequence truncated.'
                nbrterrain(i)=j-1
                exit
            END IF
            id_terrain_intern(j,i)=n                                    !replace internal ID with external ID
        END DO
    END DO

    !check fractions of TCs in LUs and normalize, if necessary
    DO i_lu = 1,nsoter
        temp1 = 0.
        DO tc_counter=1,nbrterrain(i_lu)
            !tcid_instance=tcallid(sb_counter,lu_counter,tc_counter)    !id of TC instance
            id_tc_type=id_terrain_intern(tc_counter,i_lu)            !id of TC type
            temp1 = temp1 + fracterrain(id_tc_type) !current sum of fractions
        END DO
        if (temp1>1.05 .OR. temp1<0.95) then !correct fractions
            write(*,'(A,i0,a,f5.2,a)')'WARNING: Sum of TC-fractions for LU ', id_lu_extern(i_lu),' was ',temp1,', now normalized to 1.'
            DO tc_counter=1,nbrterrain(i_lu)
                id_tc_type=id_terrain_intern(tc_counter,i_lu)            !id of TC type
                fracterrain(id_tc_type) = fracterrain(id_tc_type) /temp1   !correct fractions
            END DO
        end if
    END DO


    ! insert internal numbers of soils into
    ! Soil-Vegetation component arrays
    DO i=1,ntcinst
        DO j=1, nbr_svc(i)
            n=1
            DO WHILE ((n<=nsoil) .AND. (id_soil_extern(n) /= id_soil_intern(j,i)))            !search until the current external ID has been found...
                n=n+1
            END DO

            IF (n > nsoil) THEN    !reference to undefined SOIL
                WRITE(*,'(a, I0, a)') 'ERROR: Soil ', id_soil_intern(j,i), ' not found in soil.dat.'
                STOP
            END IF

            id_soil_intern(j,i)=n                                        !Till: replace external ID with internal ID
        END DO
    END DO

    ! insert internal numbers of vegetation units into
    ! Soil-Vegetation component arrays
    DO i=1,ntcinst
        DO j=1, nbr_svc(i)
            DO n=1,nveg
                if (id_veg_extern(n) == id_veg_intern(j,i)) exit           !search until the current external ID has been found...
            END DO

            if (n>nveg) then
                write(*,'(a, I0, a, I0, a)')'ERROR: Vegetation class ',id_veg_intern(j,i),' not found in vegetation.dat.'
                stop
            end if
            id_veg_intern(j,i)=n                                        !replace internal ID with external ID
        END DO
    END DO


    !** sort terrain components within Soter unit
    !** according to upstream-downstream relationship
    !** (the order of terrain-IDs is changed: ID of uppermost TC first in array)

    DO k=1,nsoter
        j=1
        n=0
        DO WHILE (j <= nbrterrain(k))
            n=n+1
            sortier(j)=posterrain(id_terrain_intern(j,k))
            j=j+1
        END DO
        DO i=1,n-1
            i_min=i
            DO c=i+1,n
                IF (sortier(i_min) > sortier(c)) i_min=c
            END DO
            tauschr=sortier(i_min)
            sortier(i_min)=sortier(i)
            sortier(i)=tauschr
            tausch=id_terrain_intern(i_min,k)
            id_terrain_intern(i_min,k)=id_terrain_intern(i,k)
            id_terrain_intern(i,k)=tausch
        END DO
    END DO


    tcallid=-1
    DO n=1,ntcinst  !sort in contents of soil_vegation.dat and check for excess entries
        id_sub_int = id_ext2int(luse_subbas(n), id_subbas_extern)
        id_lu_int  = id_ext2int(    luse_lu(n), id_lu_extern)
        id_tc_int  = id_ext2int(    luse_tc(n), id_terrain_extern)

        i=which1(id_lu_intern(:,id_sub_int)==id_lu_int) !get "position" of current LU in current subbasin
        if (i==0) then
            WRITE(*,'(a,i0,a,i0,a,i0)')'Error in soil_vegetation.dat, line ',(n-1)*3+4,': LU ',luse_lu(n),&
                ' not expected for subbasin ',luse_subbas(n)
            STOP
        END IF

        j=which1(id_terrain_intern(:,id_lu_int)==id_tc_int) !get "position" of current TC in current LU
        if (j==0) then
            WRITE(*,'(a,i0,a,i0,a,i0)')'Error in soil_vegetation.dat, line ',(n-1)*3+4,': TC ',luse_tc(n),&
                ' not expected for LU ',luse_lu(n)
            STOP
        END IF

        tcallid(id_sub_int,i,j) = n
    END DO

    if (allocated(svc_soil_veg)) then
        DO sb_counter=1,subasin            !check, if all relevant SVCs have been specified
            DO lu_counter=1,nbr_lu(sb_counter)
                i_lu=id_lu_intern(lu_counter,sb_counter)
                DO tc_counter=1,nbrterrain(i_lu)
                    tcid_instance=tcallid(sb_counter,lu_counter,tc_counter) !id of TC instance
                    id_tc_type=id_terrain_intern(tc_counter,i_lu)            !id of TC type
                    if (tcid_instance==-1) cycle                            !this may happen if this is merely a dummy basin with prespecified outflow
                    DO svc_counter=1,nbr_svc(tcid_instance)
                        i_soil=id_soil_intern(svc_counter,tcid_instance)        !internal id of soil type
                        i_veg=  id_veg_intern(svc_counter,tcid_instance)
                        i=which1(svc_soil_veg(:,1)==i_soil .AND. svc_soil_veg(:,2)==i_veg)
                        IF (i==0) THEN    !soil-vegetation combination not specified as an SVC
                            write(*,'(a,i0,a,i0)')'ERROR in svc.dat: could not find soil-vegetation combination ', id_soil_extern(i_soil),' :', id_veg_extern(i_veg)
                            write(*,'(a,i0,a,i0,a,i0,a,i0,a)')' found for subbasin ', id_subbas_extern(sb_counter),', LU ', id_lu_extern(i_lu), ', TC ', id_terrain_extern(id_tc_type), ' (position ', tc_counter,')'
                            stop
                        END IF
                    END DO
                END DO
            END DO
        END DO
    end if





    DO id_sub_int=1,subasin !check contents of soil_vegation.dat and for missing entries
        k = count(tcallid(id_sub_int,:,1)/=-1)
        if (k /= nbr_lu(id_sub_int)) then
            WRITE(*,'(a,i0,a,i0,a,i0)')'Error in soil_vegetation.dat: Expecting ',nbr_lu(id_sub_int),&
                ' LU(s) for subbasin ',id_subbas_extern(id_sub_int),', found ',k
            STOP
        END IF

        DO j=1,nbr_lu(id_sub_int)
            id_lu_int=id_lu_intern(j,id_sub_int)

            k=count(tcallid(id_sub_int,j,:)/=-1)
            if (k /= nbrterrain(id_lu_int)) then
                WRITE(*,'(a,i0,a,i0)')'Error in soil_vegetation.dat: Expecting ',nbrterrain(id_lu_int),&
                    ' TCs for LU ',id_lu_extern(id_lu_int),',found ',k
                STOP
            END IF
        END DO
    END DO


    !----------------------------------------------------------------------------------
    ! calculation of the fractions within each subbasin that have the irrigation flags. -> irrigated fraction of each subbasin

    IF (allocated(svc_irr)) THEN

        DO sb_counter=1, subasin !Loop over all Subbasins

            frac_irr_sub(sb_counter) = 0.

            DO lu_counter=1,nbr_lu(sb_counter)  !Loop over all LU's
                i_lu=id_lu_intern(lu_counter,sb_counter)

               frac_irr_lu = 0.

                DO tc_counter=1,nbrterrain(i_lu)  ! Loop over all TC's
                    tcid_instance=tcallid(sb_counter,lu_counter,tc_counter)    !id of TC instance
                    if (tcid_instance==-1) cycle                            !this may happen if this is merely a dummy basin with prespecified outflow
                    id_tc_type=id_terrain_intern(tc_counter,i_lu)            !id of TC type

                    frac_irr_tc = 0.

                    DO svc_counter=1,nbr_svc(tcid_instance)    !check all SVCs of the current TC, Loop over all SVC's
                       IF (svc_irr(tc_contains_svc2(id_tc_type)%p(svc_counter)%svc_id) == 1 ) THEN
                           frac_svc_x=    tc_contains_svc2(id_tc_type)%p(svc_counter)%fraction
                           frac_irr_tc = frac_irr_tc + frac_svc_x
                       END IF
                    END DO

                    frac_irr_lu = frac_irr_lu + fracterrain(id_tc_type) * frac_irr_tc
                END DO

                frac_irr_sub(sb_counter) = frac_irr_sub(sb_counter) + frac_irr_lu * frac_lu(lu_counter, sb_counter)

            END DO !LU-Loop

        END DO     ! Subbasin Loop
    END IF ! if irrigation is on


    !-------------------------------------------


    !** define van-genuchten parameter from pore-size index for
    !   calculation of unsaturated hydraulic conductivity
    DO i=1,nsoil
        DO h=1,nbrhori(i)
            porem(i,h)=poresz(i,h)/(poresz(i,h)+1.)
        END DO
    END DO

   !Till: check consistency between SVC-fractions in soil_vegetation.dat and svc_in_tc.dat
    if (allocated(tc_contains_svc2)) then !check consistency between SVC-fractions in soil_vegetation.dat and svc_in_tc.dat
        DO sb_counter=1,subasin
            DO lu_counter=1,nbr_lu(sb_counter)
                i_lu=id_lu_intern(lu_counter,sb_counter)
                DO tc_counter=1,nbrterrain(i_lu)
                    tcid_instance=tcallid(sb_counter,lu_counter,tc_counter)    !id of TC instance
                    if (tcid_instance==-1) cycle                            !this may happen if this is merely a dummy basin with prespecified outflow
                    id_tc_type=id_terrain_intern(tc_counter,i_lu)            !id of TC type

                    do svc_counter=1,nbr_svc(tcid_instance)    !check all SVCs of the current TC
                        frac_svc_x=    tc_contains_svc2(id_tc_type)%p(svc_counter)%fraction
                        if (abs(frac_svc_x - frac_svc(svc_counter,tcid_instance)) > 0.01) then
                                write(*,'(A,i0,a,i0,a)')'WARNING: Inconsistent SVC-fractions in soil_vegetation.dat and svc_in_tc.dat for TC ', id_terrain_extern(id_tc_type),', SVC at position ', svc_counter,'. Using the latter. '
                                frac_svc(svc_counter,tcid_instance) = tc_contains_svc2(id_tc_type)%p(svc_counter)%fraction  !we correct this way, as svc_in_tc.dat (if present) has a structure similar to db, which is more consistent
                        end if
                    end do

                END DO
            END DO
        END DO
    end if


    !Till: determine rocky_fraction from areal coverage of SVCs with a soil that has 100% coarse fragments in the top layer
    if (sum(rocky)==0.0) then        !this will only be done if rocky has not been specified in soil_vegetation.dat explicitely

        !    !remove rocky SVCs from tc_contains_svc2 - should not be done, they influence areal USLE factors
        !    DO sb_counter=1,subasin
        !        DO lu_counter=1,nbr_lu(sb_counter)
        !            i_lu=id_lu_intern(lu_counter,sb_counter)
        !            DO tc_counter=1,nbrterrain(i_lu)
        !                tcid_instance=tcallid(sb_counter,lu_counter,tc_counter)    !id of TC instance
        !                if (tcid_instance==-1) cycle                            !this may happen if this is merely a dummy basin with prespecified outflow
        !                id_tc_type=id_terrain_intern(tc_counter,i_lu)            !id of TC type
        !                if (id_tc_type<0) cycle                            !already treated
        !                temp1=0
        !                i=size(tc_contains_svc2(id_tc_type)%p)        !number of SVCs in current TC
        !                DO svc_counter=1,i    !loop over all SVCs in TC
        !                    soilid=id_soil_intern(svc_counter,tcid_instance)        !get soil-ID of current SVC
        !                    if (soilid==0 ) then                        !check if the topmost horizon consists of coarse fragments only
        !                        soilid=-1
        !                    end if
        !
        !                    if (soilid==-1 .OR. coarse(soilid,1)==1) then                        !check if the topmost horizon consists of coarse fragments only
        !                        tc_contains_svc2(id_tc_type)%p(svc_counter:i-1)=tc_contains_svc2(id_tc_type)%p(svc_counter+1:i)    !remove entry
        !                        tc_contains_svc2(id_tc_type)%p(svc_counter)%fraction=0.    !just to make sure is further disregarded
        !                    end if
        !                    temp1=temp1+ tc_contains_svc2(id_tc_type)%p(svc_counter)%fraction                    !sum up fraction of remaining SVCs
        !                END DO !loop SVCs
        !                tc_contains_svc2(id_tc_type)%p(:)%fraction=tc_contains_svc2(id_tc_type)%p(:)%fraction/temp1    !normalize fractions
        !                where (id_terrain_intern==id_tc_type)
        !                    id_terrain_intern=-id_tc_type        !mark these TC-types as treated by making ID negativd
        !                endwhere
        !
        !            end do
        !        end do
        !    end do
        !    id_terrain_intern=abs(id_terrain_intern)    !restore real IDs (no negatives)
        !


        DO tcid_instance=1,ntcinst    !go thru all instances of TCs
            svc_counter=1
            do while (svc_counter<=nbr_svc(tcid_instance))    !check all SVCs of the current TC
                soilid=id_soil_intern(svc_counter,tcid_instance)        !get soil-ID
                if (coarse(soilid,1)==1.) then                        !check if the topmost horizon consists of coarse fragments only
                    rocky(tcid_instance)=rocky(tcid_instance)+frac_svc(svc_counter,tcid_instance)    !increase rocky fraction of this TC by the fraction of this SVC
                    !remove SVC from this TC
                    frac_svc(svc_counter:nbr_svc(tcid_instance)-1,tcid_instance)=frac_svc(svc_counter+1:nbr_svc(tcid_instance),tcid_instance)
                    frac_svc(nbr_svc(tcid_instance),tcid_instance)=0.    !just to make sure its one less
                    id_soil_intern(svc_counter:nbr_svc(tcid_instance)-1,tcid_instance)=id_soil_intern(svc_counter+1:nbr_svc(tcid_instance),tcid_instance)
                    id_soil_intern(nbr_svc(tcid_instance),tcid_instance)=-1    !just to make sure its one less
                    id_veg_intern(svc_counter:nbr_svc(tcid_instance)-1,tcid_instance)=id_veg_intern(svc_counter+1:nbr_svc(tcid_instance),tcid_instance)
                    id_veg_intern(nbr_svc(tcid_instance),tcid_instance)=-1    !just to make sure its one less
                    nbr_svc(tcid_instance)=nbr_svc(tcid_instance)-1        !this TC has on SVC less
                else
                    svc_counter=svc_counter+1            !next SVC
                end if
            END DO    !end loop thru all SVCs

        END DO


        !    do id_tc_type=1,size(tc_contains_svc2)    !loop over all TC-types
        !        i=size(tc_contains_svc2(id_tc_type)%p)        !number of SVCs in current TC
        !        temp1=0
        !        DO svc_counter=1,i    !loop over all SVCs in TC
        !            svc_id=tc_contains_svc2(id_tc_type)%p(svc_counter)%svc_id                !get id of current SVC to be treated
        !
        !            tc_instance=0
        !
        !            DO tc_counter=1,maxterrain
        !                do i_lu=1:nsoter
        !                    if (id_terrain_intern(tc_counter,i_lu)==id_tc_type) then
        !                        tcid_instance=tcallid(i_subbas,lu_counter,tc_counter)
        !          id_tc_type=id_terrain_intern(tc_counter,i_lu)
        !
        !
        !
        !                        exit
        !                    end if
        !                end do
        !            end do
        !
        !            soilid=id_soil_intern(svc_counter,tc_instance)        !get soil-ID of current SVC
        !            if (soilid==0 ) then                        !check if the topmost horizon consists of coarse fragments only
        !                soilid=-1
        !            end if
        !
        !            if (soilid==-1 .OR. coarse(soilid,1)==1) then                        !check if the topmost horizon consists of coarse fragments only
        !                tc_contains_svc2(id_tc_type)%p(svc_counter:i-1)=tc_contains_svc2(id_tc_type)%p(svc_counter+1:i)    !remove entry
        !                tc_contains_svc2(id_tc_type)%p(svc_counter)%fraction=0.    !just to make sure is further disregarded
        !            end if
        !            temp1=temp1+ tc_contains_svc2(id_tc_type)%p(svc_counter)%fraction                    !sum up fraction of remaining SVCs
        !        END DO !loop SVCs
        !        tc_contains_svc2(id_tc_type)%p(:)%fraction=tc_contains_svc2(id_tc_type)%p(:)%fraction/temp1    !normalize fractions
        !    end do
        !
    end if





    !Till: make sure, rocky fraction and svc fraction in each SVC sum up to 1
    DO sb_counter=1,subasin
        DO lu_counter=1,nbr_lu(sb_counter)
            i_lu=id_lu_intern(lu_counter,sb_counter)
            DO tc_counter=1,nbrterrain(i_lu)
                tcid_instance=tcallid(sb_counter,lu_counter,tc_counter)    !id of TC instance
                if (tcid_instance==-1) cycle                            !this may happen if this is merely a dummy basin with prespecified outflow
                !normalize
                temp1= sum(frac_svc(:,tcid_instance)) + rocky(tcid_instance) !current sum of fractions
                if (temp1>1.05 .OR. temp1<0.95) then
                    id_tc_type=id_terrain_intern(tc_counter,i_lu)            !id of TC type
                    write(*,'(A,f5.2,a,i0)')'WARNING: Sum of SVC- and rocky fraction was ',temp1,', now normalized to 1. (TC ',&
                        id_terrain_extern(id_tc_type),')'
                end if
                rocky(tcid_instance)     = rocky(tcid_instance)     / temp1    !correct fractions
                frac_svc(:,tcid_instance)= frac_svc(:,tcid_instance)/ temp1
            END DO
        END DO
    END DO



    ! .......................................................................
    !  for all Soil-Vegetation-components (SVC):

    !  compare (1) root depth of landuse class in this SVC
    !          (2) maximum depth of soil zone in actual landscape unit of this SVC
    !              (depth to bedrock)
    !          (3) given soil depth of characteristic soil profile
    !          (4) calibration factor for profile depth down to bedrock

    !  modify data according to these rules:
    !   (1) if given soil profile has bedrock below deepest horizon:
    !       no modifications
    !   (2) if depth of bedrock is given for LU, but not for profile:
    !       extend depth of lowest soil horizon to depth of bedrock of LU
    !       (if depth of profile is larger than depth to bedrock for LU:
    !        profile depth is not modified, bedrock is assumed directly below
    !        deepest horizon)
    !   (3) if maximum root depth > total soildepth of given profile
    !       and no depth of bedrock is given :
    !       extend depth of lowest horizon to root depth

    !  for each SVC is saved:
    !   (1)  the depth of the lowest horizon
    !   (2)  flag indicating if bedrock occurs below deepest horizon or not (1/0)
    horiz_thickness(:,:,:)=0.
    maxthickness=0.    !determine maximum soil thickness
    DO n=1,ntcinst
        DO j=1,nbr_svc(n)
            test=luse_subbas(n)
            id_soil_int=id_soil_intern(j,n)

            h=nbrhori(id_soil_int)
            DO i=1,h
                horiz_thickness(n,j,i)=temp_hori_thick(id_soil_int,i) !Till: horizon thicknesses are stored for each instance of a soil
            END DO

            !  determine lowest horizon to which roots go down
            temp1=maxval(rootdep(id_veg_intern(j,n),:))
            svcrooth(n,j)=99
            temp2=0.
            i=1
            DO WHILE (i <= h .AND. svcrooth(n,j) == 99)
                temp2=temp2+horiz_thickness(n,j,i)
                IF (temp2 >= temp1) THEN
                    svcrooth(n,j)=i
                END IF
                i=i+1
            END DO


            id_lu_ext=luse_lu(n)
            id_lu_int=1
            DO WHILE (id_lu_extern(id_lu_int) /= id_lu_ext)
                id_lu_int=id_lu_int+1
            END DO

            !  if no permanent groundwater table (normal case)
            IF (gw_flag(id_lu_int) == 0 .OR. gw_flag(id_lu_int) == 1) THEN

                !  bedrock is given for soil profile
                IF (bedrock(id_soil_int) == 1) THEN
                    svcbedr(n,j)=1
                END IF

                !  no bedrock is given for soil profile or this is an alluvial soil (is always extended to maxdepth, regardless of bedrock flag)
                IF (bedrock(id_soil_int) == 0 .OR. (alluvial_flag(id_soil_int)==1)) THEN

                    !id_soil_extern(id_soil_int)
                    IF (alluvial_flag(id_soil_int)==1) THEN
                        temp1=maxdep(id_lu_int)    !Till: alluvial soils reach at least to maxdepth of LU
                    ELSE
                        temp1=meandep(id_lu_int)!Till: other soils are extended to mean soil depth (if no bedrock is specified)
                    END IF

                    !  if bedrock depth is given for LU...
                    IF (temp1>0.) THEN    !...extent lowest horizon, if necessary

                        IF (alluvial_flag(id_soil_int)==1) THEN
                            svcbedr(n,j)=bedrock(id_soil_int)    !alluvial soils are always extended but keep their original bedrock flag
                        ELSE
                            svcbedr(n,j)=1
                        END IF


                        IF (sum(temp_hori_thick(id_soil_int,:)) < temp1) THEN !Till: expand deepest horizon to bedrock, if necessary
                            IF (h > 1) THEN
                                horiz_thickness(n,j,h)=temp1-sum(temp_hori_thick(id_soil_int,1:h-1))
                            ELSE
                                horiz_thickness(n,j,h)=temp1
                            END IF
                        END IF
                        !  no bedrock (nor for profile nor for LU) given
                    ELSE
                        svcbedr(n,j)=0                                !Till: expand deepest horizon to rooting depth, if necessary
                        temp1=maxval(rootdep(id_veg_intern(j,n),:))    !get maximum rooting depth

                        IF (sum(temp_hori_thick(id_soil_int,:)) < temp1) THEN    !if rooting depth is greater than specified horizons
                            IF (h > 1) THEN                                !enlarge deepest horizon to maximum rooting depth
                                horiz_thickness(n,j,h)=temp1-sum(temp_hori_thick(id_soil_int,1:h-1))
                            ELSE
                                horiz_thickness(n,j,h)=temp1
                            END IF
                        END IF
                    END IF
                END IF

                !  landscape unit with permanent groundwater table
                !  extend lowest horizon to the depth of the riverbed if this is deeper
            ELSE IF (gw_flag(id_lu_int) == 99) THEN
                IF (h > 1) THEN
                    temp1=riverbed(id_lu_int)-sum(temp_hori_thick(id_soil_int,1:h))
                    IF (temp1 > 0.) THEN
                        horiz_thickness(n,j,h)=riverbed(id_lu_int)-sum(temp_hori_thick(id_soil_int,1:h-1))
                    END IF
                ELSE
                    temp1=riverbed(id_lu_int)-temp_hori_thick(id_soil_int,h)
                    IF (temp1 > 0.) THEN
                        horiz_thickness(n,j,h)=riverbed(id_lu_int)
                    END IF
                END IF
            END IF
        END DO

        maxthickness=max(maxthickness, maxval(sum(horiz_thickness(n,:,:),2) )  )



    END DO
    !allocate(latred(ntcinst,maxhori*4)) !Till: latred is now allocated in readhymo.f90 to adjust to maximum number of exchange horizons required
    !maxthickness=maxval(sum(horiz_thickness, 3))    !Till: compute maximum total thickness among soil profiles - crashes with large arrays, so we use the loop above


    allocate(latred(ntcinst,ceiling(maxthickness/500.)+1)) !Till: latred is now allocated in readhymo.f90 to adjust to maximum number of exchange horizons required

    ! ..................................................................


    !  for all Soil-Vegetation-components:
    !  calculate water content at permanent wilting point which depends
    !  on plant characteristics (not fixed at 15000 cm suction)
    !  Van-Genuchten formula
    !  but pwp must not be smaller than residual water content

    DO i=1,ntcinst
        DO j=1,nbr_svc(i)
            idummy=id_soil_intern(j,i)
            DO h=1,nbrhori(idummy)
                temp1=poresz(idummy,h)+1.

                if (bubble(idummy,h)<=0.) then
                    write(*,'(a,i0,a,i0)')'ERROR: Non-positive bubble pressure point for soil ',id_soil_extern(idummy), ', horizon ',h
                    stop
                end if

                pwpsc(i,j,h)=((1./(1.+(wstressmax(id_veg_intern(j,i))/  &
                    bubble(idummy,h))**temp1))** porem(idummy,h))*  &
                    (thetas(idummy,h)- thetar(idummy,h))+  &
                    thetar(idummy,h)

                pwpsc(i,j,h)=pwpsc(i,j,h)*(1.-coarse(idummy,h))
                IF (pwpsc(i,j,h) < thetar(idummy,h)*(1.-coarse(idummy,h))) THEN
                    WRITE(*,*) 'lower PWP',idummy,h,i
                END IF
            END DO
        END DO
    END DO


    !** establish arrays with distributions of soil characteristics between
    !** horizons of each soil component
    !   horiths is used in soilwat.f90 for lateral flow input (check if needed !)

    DO i=1,ntcinst
        DO j=1,nbr_svc(i)
            id_soil_int=id_soil_intern(j,i)

            temp2=sum(horiz_thickness(i,j,1:nbrhori(id_soil_int))*thetas(id_soil_int,1:nbrhori(id_soil_int))*(1.-coarse(id_soil_int,1:nbrhori(id_soil_int))))    !Till: sum up overall potential storage of this profile

            horiths(i,j,1:nbrhori(id_soil_int))=thetas(id_soil_int,1:nbrhori(id_soil_int))*horiz_thickness(i,j,1:nbrhori(id_soil_int))* (1.-coarse(id_soil_int,1:nbrhori(id_soil_int)))/temp2    !Till: relative distribution of maximum water storage among horizons of a profile
        END DO
    END DO


    !  correct water content characteristics of all horizons
    !  for coarse fragments
    DO i=1,nsoil
        DO h=1,nbrhori(i)
            soilnfk(i,h) = soilnfk(i,h)*(1.-coarse(i,h))
            thetas(i,h)  = thetas(i,h)* (1.-coarse(i,h))
            thetar(i,h)  = thetar(i,h)* (1.-coarse(i,h))
            soilpwp(i,h) = soilpwp(i,h)*(1.-coarse(i,h))
            soilfc(i,h)  = soilfc(i,h)* (1.-coarse(i,h))
            soilfc63(i,h)= soilfc63(i,h)*(1.-coarse(i,h))
        END DO
    END DO


    !** summary values of all horizons of each soil-veg-component
    DO i=1,ntcinst
        DO j=1,nbr_svc(i)
            id_soil_int=id_soil_intern(j,i)
            thsprof(i,j)=0.
            DO h=1,nbrhori(id_soil_int)
                thsprof(i,j)=thsprof(i,j)+thetas(id_soil_int,h)*horiz_thickness(i,j,h)
            END DO
        END DO
    END DO

    !** establish values/arrays of distribution functions
    !** for soil properties for soil components
    !   Saturated water content (corrected for coarse fragments)
    tctheta_s => soildistr()



    !George read of small_reservoirs.dat was removed


    !IF (doacud) THEN
    !!   fraction of generated runoff routed through small acudes
    !!   as there are five volume classes it is assumed that
    !!   each volume class receives 1/6 of generated runoff directly
    !!   and 1/6 is not captured by any small acude
    !
    !    intercepted=0.83
    !ELSE
    !    intercepted=0.0        !Till: if the small reservoirs are disabled, don't use route anything to them
    !END IF


    !George modelling unit initalisierung related to lake module was removed


    ! reminder: reservoir related parameters

    !IF (.NOT. doacud) THEN    !Till: outcommented, have not been allocated yet anyway
    !  storcap(:)=0.
    !  damflow(:)=0.
    !END IF


    ! ..........................................................................
    ! For TRANSPOSITIONS ONLY

    !**  Read data on water transpositions between sub-basins
    IF (dotrans) THEN
        OPEN(11,FILE=pfadp(1:pfadj)// 'Others/transposition.dat',STATUS='old')
        READ(11,*);READ(11,*)
        DO i=1,ntrans
            READ(11,*) trans_start(1,i),trans_start(2,i), q_trans(i),loss_trans(i), trans_end(1,i),trans_end(2,i),y_trans(i)
            ! Use Map IDs in the previous expressions
        END DO
        CLOSE(11)

        !Eva this relates the MAP IDs (id_subbas_extern(subasin)) to the sorted CODE IDs (id_subbas_intern(subasin)),
        DO i=1,ntrans
            j=1
            DO WHILE (id_subbas_extern(j) /= trans_start(1,i))
                j=j+1
                IF (j > 500) THEN
                    WRITE (*,*) 'trans_start(1,i) loop in readhymo.f'
                END IF
            END DO
            trans_start(1,i)=j !ii: do this more elegantly
        END DO
        DO i=1,ntrans
            j=1
            DO WHILE (id_subbas_extern(j) /= trans_end(1,i))
                j=j+1
                IF (j > 500) THEN
                    WRITE (*,*) 'trans_end(1,i) loop in readhymo.f'
                END IF
            END DO
            trans_end(1,i)=j !ii: do this more elegantly
        END DO

    END IF





    OPEN(11,FILE=pfadp(1:pfadj)// 'Hillslope/frac_direct_gw.dat',IOSTAT=istate,STATUS='old')
    IF (istate==0) THEN
        READ(11,*,IOSTAT=istate)frac_direct_gw
        CLOSE(11)
        IF ((frac_direct_gw>1.).OR.(frac_direct_gw<0.)) THEN
            write(*,'(a)')'WARNING: frac_direct_gw was outside [0..1], assumed to be 1.'
            frac_direct_gw=1.
        END IF
    ELSE
        frac_direct_gw=1.
    END IF




    ! Calibration factor for static wind
    OPEN(11,FILE=pfadp(1:pfadj)// 'Others/calib_wind.dat',IOSTAT=istate,STATUS='old')
    IF (istate==0) THEN
        READ(11,*,IOSTAT=istate) wind_t
        wind = wind_t
        CLOSE(11)
    ELSE
        wind = 1.
    END IF




    if (dosediment) then
        !Till: eg for scenanario with badland remediation
        OPEN(11,FILE=pfadp(1:pfadj)// 'Hillslope/sdr_lu.dat',IOSTAT=istate,STATUS='old')
        IF (istate==0) THEN
            if (allocated(sdr_tc)) then
                write(*,'(a)')'WARNING: SDRs for TCs found in terrain.dat, contents of sdr_lu.dat ignored.'
            else
                allocate(sdr_lu(nsoter))
                sdr_lu=1.    !default value: all sediment is delivered
                READ(11,*,IOSTAT=istate) !skip header
                READ(11,*,IOSTAT=istate)
                do while (istate==0)
                    READ(11,*,IOSTAT=istate)id_lu_ext,temp1
                    if (id_lu_ext==-1) then
                        where(sdr_lu == 1.) sdr_lu=temp1  !set all entries that have not been set before
                    else
                        i_lu=id_ext2int(id_lu_ext,id_lu_extern)    !convert to internal lu id
                        if (i_lu==-1) then
                            write(*,'(a,i0,a)')'Unknown LU-ID ',id_lu_ext,' in sdr_lu.dat, ignored.'
                            cycle    !proceed with next line
                        end if
                        sdr_lu(i_lu)=temp1
                    end if
                end do
            end if
            CLOSE(11)
        END IF

        if ((allocated(sdr_tc) .OR. allocated(sdr_tc)) .AND. erosion_equation /=1) then
            write(*,'(a)')'WARNING: Usage of pre-specified SDRs (terrain.dat or sdr_lu.dat) and not using USLE may not be appriate.'
        end if

        !Till: eg implementing highly erodible surface such as badlands
        OPEN(11,FILE=pfadp(1:pfadj)// 'Hillslope/beta_fac_lu.dat',IOSTAT=istate,STATUS='old')
        IF (istate==0) THEN
            if (allocated(beta_fac_tc)) then
                write(*,'(a)')'WARNING: Beta values for TCs found in terrain.dat, contents of beta_fac_lu.dat ignored.'
            else
                allocate(beta_fac(nsoter))
                beta_fac=1.    !default value: no correction of beta
                READ(11,*,IOSTAT=istate) !skip header
                READ(11,*,IOSTAT=istate)
                do while (istate==0)
                    READ(11,*,IOSTAT=istate)id_lu_ext,temp1
                    if (id_lu_ext==-1) then
                        where(sdr_lu == 1.) beta_fac=temp1  !set all entries that have not been set before
                    else
                        i_lu=id_ext2int(id_lu_ext,id_lu_extern)    !convert to internal lu id
                        if (i_lu==-1) then
                            write(*,'(a,i0,a)')'Unknown LU-ID ',id_lu_ext,' in beta_fac.dat, ignored.'
                            cycle    !proceed with next line
                        end if
                        beta_fac(i_lu)=temp1
                    end if
                end do
            end if
            CLOSE(11)
        END IF
END IF !dosediments

if (dosnow) then
    call allocate_snow()
     !** read Landscape units parameters related to snow
        OPEN(11,FILE=pfadp(1:pfadj)// 'Hillslope/lu2.dat',STATUS='old')
        READ(11,*,IOSTAT=istate); READ (11,*,IOSTAT=istate)
        h=2

        lu_aspect(:)=-999.
        lu_alt(:)   =-999.

        do while (istate==0)
            READ(11,'(a)',IOSTAT=istate) cdummy
            h=h+1 !count lines
            IF (istate/=-0) exit
            READ(cdummy,*,IOSTAT=istate) j !read subbasin ID only
            h=h+1 !count lines
		    k=id_ext2int(j, id_lu_extern)    !get internal id of this LU
                IF (k==-1) THEN    !specified LU not found
                    write(*,'(a,i0,a,i0,a)')'lu2.dat, line ',h,': Unknown LU ', j,' (not in soter.dat).'
                    cycle
                END IF

            READ(cdummy,*,IOSTAT=istate) j, lu_aspect(k), lu_alt(k)
            if (istate/=0) then    !format error
                write(*,'(a,i0)')'ERROR (lu2.dat): Format error in line ',h
                stop
            end if
        END DO
        CLOSE(11)

        DO i=1,nsoter !check completeness
            if (lu_aspect(i)==-999 .OR. lu_alt(i)==-999) then
                WRITE(*,'(a, I0, a)') 'ERROR: lu2.dat: LU ', id_lu_extern(i), ' needs data for aspect and/or altitude.'
                STOP
            END IF
        END DO

        lu_aspect = lu_aspect * pi/180 !convert degree to radiants

        DO i_lu=1,nsoter  !Add relative LU elevation to rel_elevation() for all TCs of all LUs
           DO tc_counter=1,nbrterrain(i_lu)
              rel_elevation(id_terrain_intern(tc_counter,i_lu)) = rel_elevation(id_terrain_intern(tc_counter,i_lu)) + lu_alt(i_lu)
           END DO
        END DO

end if ! do_snow


    !Till: allocation part - these variables are allocated here, because their dimension is not known before
    !!Print hydrologic variable on TC scale. If not used, DISABLE
    !!************************************************************************
    !    allocate (runoff_TC(subasin,nterrain))
    !    allocate (area_TC(subasin,nterrain))

    !Till: ii: these shouldn't be allocated if dosediment=FALSE
    allocate (sedsu(maxsoter,n_sed_class))
    allocate (sediment_subbasin(366,subasin,n_sed_class))
    allocate (sediment_subbasin_t(366,nt,subasin,n_sed_class))
    !!Print hydrologic variable on TC scale. If not used, DISABLE
    !!***********************************************************
    !    allocate (sed_yield_TC(subasin,nterrain))
    !    allocate (deposition_TC(subasin,nterrain))
    !    allocate (cum_erosion_TC(subasin,nterrain))
    !    allocate (cum_deposition_TC(subasin,nterrain))


    contains


    FUNCTION seasonality_array2(inputfile_name) result(node_days)        !read seasonality information (replaces seasonality_array)
    use params_h
    use hymo_h

    implicit none
    character(len=*), INTENT(IN)                  :: inputfile_name    !name of input file

    integer, pointer :: node_days(:,:), n2(:,:)
    integer :: i,j,k, dummy1, dummy2, dummy3, fields, irow

    OPEN(11,FILE=pfadp(1:pfadj)// 'Hillslope/'//inputfile_name,STATUS='old',IOSTAT=istate)

    i=0
    if (istate==0) then        !seasonality file found
        k = 0
        DO WHILE (k==0) !count lines in file (maximum of memory required)
            READ(11,*,IOSTAT=k)
            i=i+1
        END DO
        rewind(11) !back to start of file
        i=i-3-1 !substract headerlines and last unsuccessful line
    end if
        
    if (i <= 0) then    
        allocate(node_days(1,1))
        write(*,'(a)') 'WARNING: '//trim(inputfile_name)//' not found or empty, using defaults'
        return
    else
        allocate(node_days(i, 3+4)) !accomodate all lines; per line: subbasin-id, (veg-id/SVC-id,), year, 4 doys      !!, last_line indicator (for calc_seasonality)
        node_days(:,:)=0
        READ(11,*); READ(11,*); READ(11,*)    !read headerlines
        READ(11,'(a)') cdummy
        fields = GetNumberOfSubstrings(cdummy) !Till: count number of fields/columns

        BACKSPACE (11)		!rewind line just read

        loop = 3 !line position in file
        irow = 0 !counter for valid lines that have been read

        DO WHILE (.TRUE.)

            dummy2 = -1 !default for "all vegetation/SVC classes"
            if (fields == 6) then
                READ(11,*,IOSTAT=k) dummy1,         dummy3,(idummy2(i),i=1,4) !old format without veg/SVC-id
            else
                READ(11,*,IOSTAT=k) dummy1, dummy2, dummy3,(idummy2(i),i=1,4)
            end if
            if (k /=0) exit !no more lines

            loop = loop + 1 !increase line counter

            IF (dummy3 /= -1 .AND. ((dummy3 > tstop) .OR. (dummy3 < tstart))) cycle        !found specification for a year that is out of the simulation range

            if (dummy1 /=-1) then !wildcard for "all subbasins"
                j = id_ext2int(dummy1, id_subbas_extern)            !convert external to internal ID
                IF (j == -1) THEN                                        !found unknown subbasin id
                    WRITE(*,'(a, I0, a, I0, a)') inputfile_name//', line ', loop, ': Sub-basin-ID ', dummy1, ' not found in hymo.dat, ignored.'
                    cycle
                end if
                dummy1 = j
            end if

            !validity of veg/SVC-dis cannot be checked here yet, because they are not yet available
            irow = irow + 1 !increase counter for valid rows
            node_days(irow, :) = [dummy1, dummy2, dummy3, idummy2(1:4)]        !store read values in proper position
        END DO

        CLOSE(11)

        if (size(node_days,1) > irow) then         !less memory needed than expected, free it
            if (irow == 0) then
                allocate(n2(1, 1))   !no valid seasonality found
                write(*,'(a)') 'WARNING: no valid seasonal data found in '//trim(inputfile_name)//', using defaults.'
            else
                allocate(n2(irow, 3+4))   !accomodate all lines; per line: subbasin-id, (veg-id,), year, 4 doys, last_line indicator (for calc_seasonality)
            end if
            n2 = node_days(1:irow,:)
            deallocate(node_days)
            node_days => n2           !set return value to adjusted array
        end if

    end if !end if file present


    END FUNCTION seasonality_array2


    SUBROUTINE check_seasonality(seasonality_array, inputfile_name, external_ids)        !check completeness in seasonality data

    use hymo_h
    implicit none
    character(len=*), INTENT(IN)                  :: inputfile_name    !name of input file
    integer, pointer :: seasonality_array(:,:)
    integer, INTENT(IN) :: external_ids(:) !array holding external IDs of entities to check (e.g. SVCs, vegetation classes, ...)
    integer  :: missing(size(external_ids, dim=1)), n_missing = 0, k
    real  :: dummy(size(external_ids, dim=1), 4)
    !real, pointer :: test_ar(:)
    real :: test_ar(max(nveg, nsvc))


    dummy = .0

    DO i = tstart, tstop !check completeness
        DO j=1,subasin
            test_ar=0.
            !test_ar => calc_seasonality2(j, i, 1, seasonality_array, dummy)            !compute seasonal value of first day (actual values are not important, therefore we simply use dummy)
            test_ar(1:size(external_ids, dim=1)) = calc_seasonality2(j, i, 1, seasonality_array, dummy)            !compute seasonal value of first day (actual values are not important, therefore we simply use dummy)

            IF (any(test_ar == tiny(test_ar))) then
                IF (all(test_ar == tiny(test_ar))) then
                    WRITE(*,'(a, I0, a, I0)') 'ERROR: '//inputfile_name//': Sub-basin ', id_subbas_extern(j), &
                        ' lacks seasonality data for begin of simulation year ', i
                    STOP
                END IF
                missing = 0

                DO k=1, size(test_ar)
                    IF (test_ar(k) == tiny(test_ar)) then
                        n_missing = n_missing+1
                        missing(n_missing) = k
                    end if
                END DO

                IF (missing(1) /= 0) then        !found subbasin/veg_id for which no seasonality data has been read
                    write(fmtstr,*)'(a, I0, a, ', n_missing,'(I0,A), a, I0)' !generate format string
                    WRITE(*,fmtstr) inputfile_name//': Sub-basin ', id_subbas_extern(j),', SVC/vegetation-ID(s) ', (external_ids(missing(k)),", ", k=1, n_missing),&
                        ' lack(s) seasonality data for for start of simulation year ', i
                    STOP
                END IF
            END IF
        END DO
    END DO
    END SUBROUTINE check_seasonality

    SUBROUTINE check_seasonality_superfluous(seasonality_array, inputfile_name, external_ids)        !check for superfluous entries in seasonality data
    use hymo_h
    implicit none
    character(len=*), INTENT(IN)                  :: inputfile_name    !name of input file
    integer, pointer :: seasonality_array(:,:), s2(:,:)
    integer, INTENT(IN) :: external_ids(:) !array holding external IDs of entities to check (e.g. SVCs, vegetation classes, ...)
    integer :: i, j, dummy2
    integer :: obsolete_lines(size(seasonality_array, dim=1))

    if (size(seasonality_array, dim=2) == 1) return !no seasonality to check
    obsolete_lines = 0
    do i=1, size(seasonality_array, dim=1) !through all "lines" of seasonality file
        dummy2 = seasonality_array(i, 2)
        if (dummy2 /=-1) then !wildcard for "all vegetation" / SVCs
            j = id_ext2int(dummy2, external_ids)            !convert external to internal ID
            IF (j == -1) THEN                                        !found unknown id
                WRITE(*,'(a, I0, a, I0, a)') inputfile_name//': vegetation/SVC-ID ', dummy2, ' not found in vegetation/svc.dat, ignored.'
                obsolete_lines(i) = 1 !mark line as "obsolete"
            end if
        end if
    end do !through all "lines" of seasonality file

    if (sum(obsolete_lines) > 0) then         !some obsolete lines: remove and free memory
        allocate(s2(  size(seasonality_array, dim=1) - sum(obsolete_lines), &
            size(seasonality_array, dim=2)))    !allocate new space for cleansed seasonality

        j=1
        do i=1,size(obsolete_lines)
            if (obsolete_lines(i) == 0) then
                s2(j,:) = seasonality_array(obsolete_lines(i), :)     
                j=j+1
            end if
        end do    

        deallocate(seasonality_array)
        seasonality_array => s2           !point to cleansed array
    end if

    END SUBROUTINE check_seasonality_superfluous


    END SUBROUTINE readhymo




    SUBROUTINE read_pre_subbas_outflow    !Till: reads predefined subbasin outflow (water and sediment)
    ! status of call (0=initialization, 1=start of year)
    use common_h
    use hymo_h
    use time_h
    use utils_h
    use erosion_h
    use params_h
    IMPLICIT NONE

    INTEGER :: i,j,dc,istate
    INTEGER   :: columnheader(2*subasin+10)    !Till: for storing column headings of input files
    CHARACTER(len=1000) :: linedummy
    CHARACTER(len=30) :: dstr    !Till: dummies for reading input header
    INTEGER,save  :: no_columns(2)        !number of columns of input files for the input file

    !first call of function: memory allocation and seeking the required file position
    IF (any(do_pre_outflow) .AND. .NOT. allocated(pre_subbas_outflow)) THEN
        do_pre_outflow=.FALSE.
        OPEN(91,FILE=pfadp(1:pfadj)// '/Time_series/subbasin_out.dat',IOSTAT=i,STATUS='old')
        IF (i==0) THEN
            write(*,'(a)', advance='no')'Reading pre-specified subbasin outflow from subbasin_out.dat...'

            READ(91,*)        !skip headerlines
            READ(91,*)

            READ(91,'(a)') linedummy
            columnheader=0
            no_columns(1)=GetNumberOfSubstrings(linedummy)-2    !Till: count number of columns
            IF (no_columns(1) < 1) then        !nothing found, issue warning
                WRITE(*,'(A)') 'WARNING: subbasin_out.dat contains no valid header columns, check format'
            END IF
            allocate (pre_subbas_outflow(366,nt,no_columns(1)))
            pre_subbas_outflow(:,:,:)=-2.

            READ (linedummy,*) dstr, dstr, (columnheader(i), i=1,no_columns(1))    !Till: extract column headers
            corr_column_pre_subbas_outflow=>set_corr_column(columnheader, 'subbasin_out.dat')
            WHERE(corr_column_pre_subbas_outflow/=0)
                nbr_lu=0        !Till: set number of LUs for prespecified subbasins to 0
                do_pre_outflow=.TRUE.
            ENDWHERE
            if (sum(corr_column_pre_subbas_outflow)==0) then
                write(*,'(a)')'no relevant subbasins found.' 
                deallocate(pre_subbas_outflow)
            else
                write(*,'(a)', advance='no')'found subbasin(s): ' 
                DO i=1, subasin
                    if (corr_column_pre_subbas_outflow(i) /=0) write(*,"(i0, A)", advance='no') id_subbas_extern(i),", "
                END DO
                write(*,'(a)')'.' !end line of output
                call date_seek(91,tstart,mstart, dstart, 'subbasin_out.dat')    !set internal filepointer to correct line in file
            end if
        END IF

        OPEN(11,FILE=pfadn(1:pfadi)//'parameter.out', STATUS='old',POSITION='append', IOSTAT=istate) !write to log file
        write(11,*) "using prespecified outflow:",       allocated(pre_subbas_outflow)
        close(11)
    END IF

    IF (.NOT. dosediment) THEN
        do_pre_outsed=.FALSE.
    ELSE
        IF (do_pre_outsed .AND. .NOT. allocated(pre_subbas_outsed)) THEN    !first call of function
            do_pre_outsed=.FALSE.
            OPEN(92,FILE=pfadp(1:pfadj)// '/Time_series/subbasin_outsed.dat',IOSTAT=i,STATUS='old')
            IF (i==0) THEN
                write(*,'(a)', advance='no')'Reading pre-specified subbasin sediment output from subbasin_outsed.dat...'

                allocate (pre_psd(n_sed_class))    !allocate storage for mean particle size distribution

                !READ(92,*)        !skip headerlines
                READ(92,'(a)') linedummy    !read first line
                DO i=len_trim(linedummy),1,-1
                    if (linedummy(i:i)==':') then
                        linedummy=linedummy(i+1:len_trim(linedummy))    !extract what is behind :
                        exit
                    end if
                END DO

                READ (linedummy,*,IOSTAT=i) pre_psd

                if (abs(1.-sum(pre_psd))>0.05) then
                    stop 'ERROR: trouble reading mean PSD from subbasin_outsed.dat. Not specified, format error or sum<>1'
                end if

                READ(92,*)
                READ(92,'(a)') linedummy
                columnheader=0
                no_columns(2)=GetNumberOfSubstrings(linedummy)-2    !Till: count number of columns
                allocate (pre_subbas_outsed(366,nt,no_columns(2)))
                pre_subbas_outsed(:,:,:)=-2.

                READ (linedummy,*) dstr, dstr, (columnheader(i), i=1,no_columns(2))    !Till: extract column headers
                corr_column_pre_subbas_outsed=>set_corr_column(columnheader, 'subbasin_outsed.dat')
                if (sum(corr_column_pre_subbas_outsed)==0) then
                    write(*,'(a)')'no relevant subbasins found.' 
                    deallocate(pre_subbas_outsed)
                else
                    do_pre_outsed=.TRUE.
                    write(*,'(a)', advance='no')'found subbasin(s): ' 
                    DO i=1, subasin
                        if (corr_column_pre_subbas_outsed(i) /=0) write(*,"(i0, A)", advance='no') id_subbas_extern(i),", "
                    END DO
                    write(*,'(a)')'.' !end line of output
                    call date_seek(92,tstart,mstart,dstart, 'subbasin_outsed.dat')    !set internal filepointer to correct line in file
                end if
            END IF
            OPEN(11,FILE=pfadn(1:pfadi)//'parameter.out', STATUS='old', POSITION='append', IOSTAT=istate) !write to log file
            write(11,*) "using prespecified sediment flux:", allocated(pre_subbas_outsed)
            close(11)
        END IF !do_pre_outsed
    END IF !dosediment


    !subsequent calls of function: actual reading of data
    IF (any(do_pre_outflow)) then
        if (t/=tstart .OR. pre_subbas_outflow(1,1,1)==-2.) THEN    !this clumsy check ensures that we do not read further when this function is called twice before the actual simulation starts
            !read outflow
            pre_subbas_outflow(:,:,:)=-1.
            READ(91, *, IOSTAT=istate) ((dstr,dstr,(pre_subbas_outflow(dc,j,i),i=1,no_columns(1)),j=1,nt),dc=1,dayyear)        !Reads in daily time series for pre-specified outflow from upstream basins
            IF (istate/=0) THEN
                write(*,'(a)')'ERROR: Premature end of file subbasin_out.dat.'
                stop
            END IF

            WHERE(pre_subbas_outflow < 0.)
                pre_subbas_outflow=-1.        !Till: mask negative values
            ENDWHERE
        END IF

    END IF


    IF (allocated(pre_subbas_outsed)) THEN
        IF (do_pre_outsed .AND. (t/=tstart .OR. pre_subbas_outsed(1,1,1)==-2.)) THEN    !this clumsy check ensures that we do not read further when this function is called twice before the actual simulation starts
            !read sediment outflow
            pre_subbas_outsed(:,:,:)=-1.
            if (do_pre_outsed) then
                READ(92, *, IOSTAT=istate) ((dstr,dstr,(pre_subbas_outsed(dc,j,i),i=1,no_columns(2)),j=1,nt),dc=1,dayyear)        !Reads in daily time series for pre-specified outflow from upstream basins
                IF (istate/=0) THEN
                    write(*,'(a)')'ERROR: Premature end of file subbasin_outsed.dat.'
                    stop
                END IF

                WHERE(pre_subbas_outsed < 0.)
                    pre_subbas_outsed=-1000.        !Till: mask negative values
                ENDWHERE
            end if
        END IF
    END IF


    RETURN


    contains
    FUNCTION set_corr_column(input_header,inputfile_name)
    !provide "pointer" to to match order of subbasins to that in hymo.dat
    use params_h
    use hymo_h

    implicit none
    integer, pointer :: set_corr_column(:)
    INTEGER, INTENT(IN)                  :: input_header(:)    !order of subbasins in input file
    character(len=*), INTENT(IN)                  :: inputfile_name    !name of input file
    integer    :: i,j

    allocate(set_corr_column(subasin))
    set_corr_column=0

    DO i=1,subasin
        DO j=1,size(input_header)
            IF(input_header(j) == id_subbas_extern(i)) THEN
                set_corr_column(i)= j    !Till: for each subbasin, find position of corresponding column in input file
                exit
            END IF

        END DO

    END DO
    END FUNCTION set_corr_column



    END SUBROUTINE read_pre_subbas_outflow





