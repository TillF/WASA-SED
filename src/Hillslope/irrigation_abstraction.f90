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
        gw_abstraction = 0 !Till: besser mit "0.0" (Real-Zahl) initialisieren, sonst Fehlermeldung bei manchen Compilern
        sufficient_LUs = 0

        DO i=1, nbr_irri_records !Find all groundwater extractions FROM current subbasin and sum them up
		!vielleicht besser mit PACK alle Array-Indizes holen, die wir brauchen -> "current_receivers"
		!gw_abstraction_requested = sum(irri_rate(current_receivers))
            IF (sub_source(i) == sb_counter .AND. irri_source(i) == "groundwater") THEN
                gw_abstraction = gw_abstraction + irri_rate(i)
            END IF
        END DO

        IF (gw_abstraction > 0. ) THEN
			!gw_abstraction_available = sum(deepgw(1:nbr_lu(sb_counter),sb_counter)
            DO i=1, nbr_lu(sb_counter)
            ! count all available LUs in current subbasin that have enough groundwater equaling the total abstraction demand
                IF (deepgw(i,sb_counter) > gw_abstraction) THEN !Till: Warum muss JEDE EINZELNE LU soviel GW bieten, wie gefordert ist? Besser: insgesamt (Summe) muss genügend GW da sein.
				! wenn nicht genügend vorhanden, dann gw_abstraction auf verfügbare Menge reduzieren:
				! gw_abstraction_requested = gw_abstraction_available
                    sufficient_LUs = sufficient_LUs + 1
                END IF
            END DO

            IF ( sufficient_LUs == 0 ) THEN   !Warning if abstraction demand can't be met by any LU within current subbasin
                WRITE(*,'(a,I0,a)') 'WARNING: Not enough groundwater in Subbasin ',id_subbas_extern(sb_counter), ' to meet abstraction rate. No extraction for this timestep.'
            END IF

            !Till: ich würde die Entnahme nicht glaichmäßih über alle LUs verteilen, sondern hinsichtlich der verfügbaren Speicher aufteilen (i.e. aus großen Speichern wird auch viel entnommen, das dürfte in der Praxis auch oft so sein:
			!deepgw(1:nbr_lu(sb_counter),sb_counter) = deepgw(1:nbr_lu(sb_counter),sb_counter) - &
			!deepgw(1:nbr_lu(sb_counter),sb_counter)/gw_abstraction_available * gw_abstraction_requested
			
			DO i=1, nbr_lu(sb_counter)            ! take abstraction equally divided out of all LUs that have enough storage
                IF (deepgw(i,sb_counter) > gw_abstraction) THEN
                    deepgw(i,sb_counter) = deepgw(i,sb_counter) - (gw_abstraction / sufficient_LUs) ! Numerischer Fehler? Statt 900 werden nur 896 abgezogen
                END IF
            END DO
			
			
        END IF

    END DO ! Subbasin Loop



















END SUBROUTINE
