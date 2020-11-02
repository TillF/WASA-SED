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
    INTEGER :: sb_counter, i
    REAL :: abstraction_requested, abstraction_available
    INTEGER, allocatable :: receiver_basins(:)


    irri_supply = 0.0

    !-- Abstraction from groundwater, Loop through all subbasins and check if they're included in irri.dat-----------------------
    DO sb_counter = 1, subasin ! Subbasin Loop
        abstraction_requested = 0.0
        abstraction_available = 0.0

       ! get the indices of irri.dat for all subbasins that receive irrigation water from current subbasin
        receiver_basins = pack([(i, i=1, nbr_irri_records)], MASK=sub_source == sb_counter .AND. irri_source=="groundwater")
        abstraction_requested = sum(irri_rate(receiver_basins)) ! get total amount of abstracted water in current subbasin

        IF (abstraction_requested== 0.0) THEN
            cycle
        END IF

        IF (abstraction_requested > 0. ) THEN
			abstraction_available = sum(deepgw(1:nbr_lu(sb_counter),sb_counter)) !total available deep groundwater

                IF (abstraction_available < abstraction_requested) THEN
				    abstraction_requested = abstraction_available  !If theres not enough water to fulfill demand, use all available water
				    WRITE(*,'(a,I0,a)') 'WARNING: Not enough deep groundwater in Subbasin ',id_subbas_extern(sb_counter), ' to meet abstraction rate. All available (deep) groundwater will be used.'
                END IF
        END IF

        !  Extraction of groundwater from all the LU's in Subbasin proportionally to the amount of groundwater stored in each LU. (Full LU's give a lot of water, empty ones just a bit)
	    deepgw(1:nbr_lu(sb_counter),sb_counter) = deepgw(1:nbr_lu(sb_counter),sb_counter) - deepgw(1:nbr_lu(sb_counter),sb_counter)/abstraction_available * abstraction_requested

        !Write extracted water for each reciever basin but skip external reciever basins (code 9999) since this water just disappears outside the model
        receiver_basins = pack([(i, i=1, nbr_irri_records)], MASK=sub_source == sb_counter .AND. irri_source=="groundwater" .AND. sub_receiver /= 9999)
        abstraction_requested = sum(irri_rate(receiver_basins)) ! get total amount of abstracted water in current subbasin minus the amount that will go to external

		irri_supply(sub_receiver(receiver_basins)) = irri_supply(sub_receiver(receiver_basins)) + abstraction_requested * irri_rate(receiver_basins)/sum(irri_rate(receiver_basins)) !This vector contains the total irrigation ammount for every subbasin
    END DO ! Subbasin Loop

    !--Add water from external sources (9999) to irri_supply ----------------------------------------------

    receiver_basins = pack([(i, i=1, nbr_irri_records)], MASK=sub_source == 9999) !get indices of all external sources
    abstraction_requested = sum(irri_rate(receiver_basins))  !total amount of water coming from external sources

    IF (sum(irri_rate(receiver_basins)) > 0.0) THEN
        irri_supply(sub_receiver(receiver_basins)) = irri_supply(sub_receiver(receiver_basins)) + abstraction_requested * irri_rate(receiver_basins)/sum(irri_rate(receiver_basins))
    END IF








 !Vektor mit allen Entnahmen (sub_source) analog zu irri_supply für spätere log.file. evt. eine für jede Art (groundwater/river/reservoir)

! qout für Entnahme river





















END SUBROUTINE
