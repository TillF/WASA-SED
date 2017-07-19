module snow_h
    use common_h
    save

    ! SNOW MODULE
    real, allocatable :: snowEnergyCont(:,:,:)
    real, allocatable :: snowWaterEquiv(:,:,:)
    real, allocatable :: snowAlbedo(:,:,:)
    
    real, allocatable :: snowTemperature(:,:,:)
    !relative elevation of TC above foot of toposequence/LU (i.e. river) [m]
    real, allocatable :: rel_elevation(:)

contains
!modify meteo-drivers accoding to time and location
    SUBROUTINE snow_prepare_input(hh, day, i_subbas2, lu_counter2, tc_counter2,
        prec_mod, temp_mod, rad_mod)
        IMPLICIT NONE
        INTEGER, INTENT(IN)                  :: hh, day, i_subbas2, lu_counter2, tc_counter2
        REAL, INTENT(INOUT)             :: prec_mod, temp_mod, rad_mod
        
        !simple, lapse-rate-based modification of meteo-drivers
        prec_mod = prec_mod + lapse_prec * rel_elevation(tcallid(sb_counter, lu_counter, tc_counter))
        temp_mod = temp_mod + lapse_temp * rel_elevation(tcallid(sb_counter, lu_counter, tc_counter))
        rad_mod  = rad_mod  + lapse_rad  * rel_elevation(tcallid(sb_counter, lu_counter, tc_counter))
        
    END SUBROUTINE snow_prepare_input

end module snow_h
