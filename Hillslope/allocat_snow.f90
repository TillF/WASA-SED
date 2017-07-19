!Allocate variables related to the snow module

    idummy = 0 !for summing up IO-error flags

    allocate(
       snowEnergyCont(366,nt,ntcinst), &
       snowWaterEquiv(366,nt,ntcinst), &
       snowAlbedo(366,nt,ntcinst), &
       STAT = istate)
       idummy = istate + idummy

    !Initial values states
    snowEnergyCont(1,1,)     =     0.
    snowWaterEquiv(1,1,)     =     0.
    albedo(1,1,)             =     albedoMax

    allocate(rel_elevation(nterrain), STAT = istate)
       idummy = istate + idummy

    if (f_snowTemp) allocate(snowTemp(366,nt,ntcinst), STAT = istate)
    idummy = istate + idummy
    if (f_surfTemp) allocate(surfTemp(366,nt,ntcinst), STAT = istate)
    idummy = istate + idummy
    if (f_liquFrac) allocate(liquFrac(366,nt,ntcinst), STAT = istate)
    idummy = istate + idummy
    if (f_fluxPrec) allocate(fluxPrec(366,nt,ntcinst), STAT = istate)
    idummy = istate + idummy
    if (f_fluxSubl) allocate(fluxSubl(366,nt,ntcinst), STAT = istate)
    idummy = istate + idummy
    if (f_fluxFlow) allocate(fluxFlow(366,nt,ntcinst), STAT = istate)
    idummy = istate + idummy
    if (f_fluxNetS) allocate(fluxNetS(366,nt,ntcinst), STAT = istate)
    idummy = istate + idummy
    if (f_fluxNetL) allocate(fluxNetL(366,nt,ntcinst), STAT = istate)
    idummy = istate + idummy
    if (f_fluxSoil) allocate(fluxSoil(366,nt,ntcinst), STAT = istate)
    idummy = istate + idummy
    if (f_fluxSens) allocate(fluxSens(366,nt,ntcinst), STAT = istate)
    idummy = istate + idummy
    if (f_stoiPrec) allocate(stoiPrec(366,nt,ntcinst), STAT = istate)
    idummy = istate + idummy
    if (f_stoiSubl) allocate(stoiSubl(366,nt,ntcinst), STAT = istate)
    idummy = istate + idummy
    if (f_stoiFlow) allocate(stoiFlow(366,nt,ntcinst), STAT = istate)
    idummy = istate + idummy
    if (f_rateAlbe) allocate(rateAlbe(366,nt,ntcinst), STAT = istate)
    idummy = istate + idummy

    if (idummy/=0) then
        write(*,'(A,i0,i0,a)')'Memory allocation error (',istate,'/',idummy,') in hymo-module (snow): '
        stop
    end if





    !compute relative elevation of TC centre above foot of toposequence/LU (i.e. river) [m]
    DO i_lu=1,nsoter
        temp2 = 0. !cumulative elevation of downslope TCs
        DO tc_counter=1,nbrterrain(i_lu) 
           rel_elevation(id_terrain_intern(tc_counter,i_lu) = fracterrain(tc_counter)*slope(tc_counter)*slope(tc_counter)/2 + &   !centre of TC plus ...
                                                            temp2  ! cumulative downslope elevation
           !temp1 = temp1 + slength(i_lu)*fracterrain(tc_counter) ! cumulate downslope x-extent
           temp2 = temp2 + slength(i_lu)*fracterrain(tc_counter)*slope(tc_counter) ! cumulate downslope x-extent
       END DO
    END DO
       
