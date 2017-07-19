!Till: allocation of variables related to snow module
	 
    idummy = 0 !for summing up IO-error flags
	allocate(
	snowEnergyCont(366,nt,ntcinst), &
    snowWaterEquiv(366,nt,ntcinst), &
    snowAlbedo(366,nt,ntcinst), &
	STAT = istate)
	idummy = istate + idummy
    
	allocate(rel_elevation(nterrain), STAT = istate)
    idummy = istate + idummy
    
    if (f_snowalbedo) allocate(snowTemperature(366,nt,ntcinst), STAT = istate)
    idummy = istate + idummy
    
    
    if (idummy/=0) then
		write(*,'(A,i0,i0,a)')'Memory allocation error (',istate,'/',idummy,') in hymo-module (snow): ' 
		stop
    end if
    
    ! compute relative elevation of TC centre above foot of toposequence/LU (i.e. river) [m]
    DO i_lu=1,nsoter
        temp2 = 0. !cumulative elevation of downslope TCs
        DO tc_counter=1,nbrterrain(i_lu) 
           rel_elevation(id_terrain_intern(tc_counter,i_lu) = fracterrain(tc_counter)*slope(tc_counter)*slope(tc_counter)/2 + &   !centre of TC plus ...
                                                            temp2  ! cumulative downslope elevation
           !temp1 = temp1 + slength(i_lu)*fracterrain(tc_counter) ! cumulate downslope x-extent
           temp2 = temp2 + slength(i_lu)*fracterrain(tc_counter)*slope(tc_counter) ! cumulate downslope x-extent
       END DO
    END DO
       