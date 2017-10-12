!Allocate variables related to the snow module

    idummy = 0 !for summing up IO-error flags

    allocate(snowEnergyCont(366,nt,ntcinst), &
             snowWaterEquiv(366,nt,ntcinst), &
             snowAlbedo    (366,nt,ntcinst), &
             snowCover     (366,nt,ntcinst), &
             precipMod     (366,nt,ntcinst), &
             STAT = istate)
    idummy = istate + idummy

    !Initial values states
    snowEnergyCont(1:366,1:nt,:)     =     0.
    snowWaterEquiv(1:366,1:nt,:)     =     0.
    snowAlbedo(1:366,1:nt,:)         =     albedoMax

    allocate(rel_elevation(nterrain), STAT = istate)
    idummy = istate + idummy

    if (f_radiMod) then
       allocate(radiMod(366,nt,ntcinst), STAT = istate)
    else
       allocate(radiMod(1,1,1), STAT = istate)
    end if
    idummy = istate + idummy

    if (f_temperaMod) then
       allocate(temperaMod(366,nt,ntcinst), STAT = istate)
    else
       allocate(temperaMod(1,1,1), STAT = istate)
    end if
    idummy = istate + idummy

    if (f_cloudFrac) then
       allocate(cloudFrac(366,nt,ntcinst), STAT = istate)
    else
       allocate(cloudFrac(1,1,1), STAT = istate)
    end if
    idummy = istate + idummy

    if(f_snowTemp) then
       allocate(snowTemp(366,nt,ntcinst), STAT = istate)
    else
       allocate(snowTemp(1,1,1), STAT = istate)
    end if
    idummy = istate + idummy

    if (f_surfTemp) then
       allocate(surfTemp(366,nt,ntcinst), STAT = istate)
    else
       allocate(surfTemp(1,1,1), STAT = istate)
    end if    
    idummy = istate + idummy

    if (f_liquFrac) then
       allocate(liquFrac(366,nt,ntcinst), STAT = istate)
    else
       allocate(liquFrac(1,1,1), STAT = istate)
    end if
    idummy = istate + idummy

    if (f_fluxPrec) then 
       allocate(fluxPrec(366,nt,ntcinst), STAT = istate)
    else
       allocate(fluxPrec(1,1,1), STAT = istate)
    end if
    idummy = istate + idummy

    if (f_fluxSubl) then 
       allocate(fluxSubl(366,nt,ntcinst), STAT = istate)
    else
       allocate(fluxSubl(1,1,1), STAT = istate)
    end if
    idummy = istate + idummy

    if (f_fluxFlow) then 
       allocate(fluxFlow(366,nt,ntcinst), STAT = istate)
    else
       allocate(fluxFlow(1,1,1), STAT = istate)
    end if
    idummy = istate + idummy
    
    if (f_fluxNetS) then 
       allocate(fluxNetS(366,nt,ntcinst), STAT = istate)
    else
       allocate(fluxNetS(1,1,1), STAT = istate)
    end if
    idummy = istate + idummy

    if (f_fluxNetL) then 
       allocate(fluxNetL(366,nt,ntcinst), STAT = istate)
    else
       allocate(fluxNetL(1,1,1), STAT = istate)
    end if
    idummy = istate + idummy

    if (f_fluxSoil) then
       allocate(fluxSoil(366,nt,ntcinst), STAT = istate)
    else
       allocate(fluxSoil(1,1,1), STAT = istate)
    end if
    idummy = istate + idummy

    if (f_fluxSens) then 
       allocate(fluxSens(366,nt,ntcinst), STAT = istate)
    else
       allocate(fluxSens(1,1,1), STAT = istate)
    end if
    idummy = istate + idummy

    if (f_stoiPrec) then
       allocate(stoiPrec(366,nt,ntcinst), STAT = istate)
    else
       allocate(stoiPrec(1,1,1), STAT = istate)
    end if
    idummy = istate + idummy
    
    if (f_stoiSubl) then
       allocate(stoiSubl(366,nt,ntcinst), STAT = istate)
    else
       allocate(stoiSubl(1,1,1), STAT = istate)
    end if
    idummy = istate + idummy

    if (f_stoiFlow) then
       allocate(stoiFlow(366,nt,ntcinst), STAT = istate)
    else
       allocate(stoiFlow(1,1,1), STAT = istate)
    end if
    idummy = istate + idummy

    if (f_rateAlbe) then
       allocate(rateAlbe(366,nt,ntcinst), STAT = istate)
    else
       allocate(rateAlbe(1,1,1), STAT = istate)
    end if
    idummy = istate + idummy

    if (do_rad_corr .OR. do_alt_corr) then
        allocate(lu_aspect(nsoter), STAT = istate)
        allocate(lu_alt(nsoter), STAT = istate)
    end if
    idummy = istate + idummy
    
    if (idummy/=0) then
        write(*,'(A,i0,i0,a)')'Memory allocation error (',istate,'/',idummy,') in hymo-module (snow): '
        stop
    end if



    !Computation relative elevation of TC centre above foot of toposequence/LU (i.e. river) [m]
    DO i_lu=1,nsoter
       temp2 = 0. !cumulative elevation of downslope TCs
       DO tc_counter=1,nbrterrain(i_lu)
           !Relative elevation of center TC
           !Note: slope(tc_counter) in %
           rel_elevation(id_terrain_intern(tc_counter,i_lu)) = temp2 + slength(i_lu)*fracterrain(tc_counter)*slope(tc_counter)/100/2
           temp2 = temp2 + slength(i_lu)*fracterrain(tc_counter)*slope(tc_counter)/100 !cumulate y-extent TCs
       END DO

       DO tc_counter=1,nbrterrain(i_lu)
           rel_elevation(id_terrain_intern(tc_counter,i_lu)) = rel_elevation(id_terrain_intern(tc_counter,i_lu)) - (temp2 / 2.) !"normalize" to middle elevation of the LU
           !Relative elevation of LU not available at this point; read in later; added in readhymo after read in of lu2.dat
       END DO

    END DO
       
