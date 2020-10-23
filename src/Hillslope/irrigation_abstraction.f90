SUBROUTINE irrigation_abstraction

use lake_h
use climo_h
use common_h
use hymo_h
use params_h
use routing_h
use time_h
use reservoir_h
use utils_h

    IMPLICIT NONE
    !counters
    INTEGER :: sb_counter, i, sufficient_LUs
    REAL :: gw_abstraction

    !-- Abstraction from groundwater, Loop through all subbasins and check if they're included in irri.dat
    DO sb_counter = 1, subasin ! Subbasin Loop
        gw_abstraction = 0
        sufficient_LUs = 0

        DO i=1, nbr_irri_records !Find all groundwater extractions in current subbasin and sum them up
            IF (sub_source(i) == sb_counter .AND. irri_source(i) == "groundwater") THEN
                gw_abstraction = gw_abstraction + irri_rate(i)
            END IF
        END DO

        IF (gw_abstraction > 0 ) THEN
            DO i=1, maxsoter            ! count all available LUs in current subbasin that have enough groundwater equaling the total abstraction demand
                IF (deepgw(i,sb_counter) > gw_abstraction) THEN
                    sufficient_LUs = sufficient_LUs + 1
                END IF
            END DO

            IF ( sufficient_LUs == 0 ) THEN   !Warning if abstraction demand can't be met by any LU within current subbasin
                WRITE(*,'(a,I0,a)') 'WARNING: Not enough groundwater in Subbasin ',id_subbas_extern(sb_counter), ' to meet abstraction rate. No extraction for this timestep.'
            END IF

            DO i=1, maxsoter            ! take abstraction equally divided out of all LUs that have enough storage
                IF (deepgw(i,sb_counter) > gw_abstraction) THEN
                    deepgw(i,sb_counter) = deepgw(i,sb_counter) - (gw_abstraction / sufficient_LUs) ! Numerischer Fehler? Statt 900 werden nur 896 abgezogen
                END IF
            END DO
        END IF

    END DO ! Subbasin Loop



















END SUBROUTINE
