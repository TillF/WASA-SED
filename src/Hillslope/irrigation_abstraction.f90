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
    INTEGER :: sb_counter, i     !counters
    REAL :: abstraction_requested, abstraction_available, all_request
    INTEGER, allocatable :: receiver_basins(:)

    irri_supply = 0.0
    irri_abstraction = 0.0

    !-- Abstraction from groundwater, Loop through all subbasins and check if they're included in irri.dat-----------------------
    DO sb_counter = 1, subasin ! Subbasin Loop for groundwater abstraction
        abstraction_requested = 0.0
        abstraction_available = 0.0


       ! get the indices of irri.dat for all subbasins that receive irrigation water from current subbasin
        receiver_basins = pack([(i, i=1, nbr_irri_records)], MASK =sub_source == sb_counter .AND. irri_source == "groundwater")
        abstraction_requested = sum(irri_rate(receiver_basins)) ! get total amount of abstracted water in current subbasin
        all_request = abstraction_requested !needed if demand is higher than availibilty to calculate the right fractions

        IF (abstraction_requested== 0.0) THEN
            cycle
        END IF

        IF (abstraction_requested > 0. ) THEN
			abstraction_available = sum(deepgw(1:nbr_lu(sb_counter),sb_counter)) !total available deep groundwater in current subbasin

		    IF (abstraction_available == 0.0) THEN !On first day the storage might be empty. Skip this day
		        cycle
		    END IF

            IF (abstraction_available < abstraction_requested) THEN
			    abstraction_requested = abstraction_available  !If theres not enough water to fulfill demand, use all available water
			WRITE(*,'(a,I0,a)') 'WARNING: Not enough deep groundwater in Subbasin ',id_subbas_extern(sb_counter), ' to meet abstraction rate. All available (deep) groundwater will be used.'
            END IF
        END IF

        irri_abstraction(sb_counter) = irri_abstraction(sb_counter) + abstraction_requested  !Write extracted water from current subbasin !

        !  Extraction of groundwater from all the LU's in Subbasin proportionally to the amount of groundwater stored in each LU. (Full LU's give a lot of water, empty ones just a bit)
	    deepgw(1:nbr_lu(sb_counter),sb_counter) = deepgw(1:nbr_lu(sb_counter),sb_counter) - deepgw(1:nbr_lu(sb_counter),sb_counter)/abstraction_available * abstraction_requested

        !Write extracted water for each receiver basin but skip external receiver basins (code 9999) since this water just disappears outside the model
        receiver_basins = pack([(i, i=1, nbr_irri_records)], MASK=sub_source == sb_counter .AND. irri_source=="groundwater" .AND. sub_receiver /= 9999)
        irri_supply(sub_receiver(receiver_basins)) = irri_supply(sub_receiver(receiver_basins)) + irri_rate(receiver_basins)/all_request * abstraction_requested  !This vector contains the total irrigation ammount for every subbasin

    END DO ! Subbasin Loop for groundwater abstraction

    !--Add water from external sources (9999) to irri_supply ----------------------------------------------

    receiver_basins = pack([(i, i=1, nbr_irri_records)], MASK=sub_source == 9999) !get indices of all external sources

    IF (sum(irri_rate(receiver_basins)) > 0.0) THEN
        irri_supply(sub_receiver(receiver_basins)) = irri_supply(sub_receiver(receiver_basins)) + irri_rate(receiver_basins)
    END IF

    !-- Abstraction from river, Loop through all subbasins and check if they're included in irri.dat-----------------------

    DO sb_counter = 1, subasin ! Subbasin Loop for river
        abstraction_requested = 0.0
        abstraction_available = 0.0

       ! get the indices of irri.dat for all subbasins that receive irrigation water from current subbasin
        receiver_basins = pack([(i, i=1, nbr_irri_records)], MASK= sub_source == sb_counter .AND. irri_source=="river")
        abstraction_requested = sum(irri_rate(receiver_basins)) ! get total amount of abstracted water in current subbasin
        all_request = abstraction_requested !needed if demand is higher than availibilty to calculate the right fractions

        IF (abstraction_requested== 0.0) THEN
            cycle
        END IF

        IF (abstraction_requested > 0. ) THEN
            abstraction_available = qout(1,sb_counter) !river water in first timestep

            IF (abstraction_available == 0.0) THEN !On first day the storage might be empty. Skip this day
                cycle
            END IF

            IF (abstraction_available < abstraction_requested) THEN
                abstraction_requested = abstraction_available  !If theres not enough water to fulfill demand, use all available water from river flow
                WRITE(*,'(a,I0,a)') 'WARNING: Not enough river flow in Subbasin ',id_subbas_extern(sb_counter), ' to meet abstraction rate. All available river flow will be used.'
            END IF
        END IF

        irri_abstraction(sb_counter) = irri_abstraction(sb_counter) + abstraction_requested  !Write extracted water from current subbasin !

        !  Abstraction of river water from first timestep of river routing
        qout(1,sb_counter) = qout(1,sb_counter) - abstraction_requested

        !Write abstracted water for each receiver basin but skip external receiver basins (code 9999) since this water just disappears outside the model. Indexing irri_supply with 9999  would crash the program
        receiver_basins = pack([(i, i=1, nbr_irri_records)], MASK=sub_source == sb_counter .AND. irri_source=="river" .AND. sub_receiver /= 9999)
        irri_supply(sub_receiver(receiver_basins)) = irri_supply(sub_receiver(receiver_basins)) + irri_rate(receiver_basins)/all_request * abstraction_requested  !This vector contains the total irrigation ammount for every subbasin

    END DO ! Subbasin Loop for river abstraction





     !irri_supply und irri_abstraction in irri_abstraction_record und irri_receiver_record schreiben -> Ausgabedatei um Zwischenstände abzuspeichern

     irri_abstraction_record(d,nt,1:subasin) = irri_abstraction
     irri_supply_record(d,nt,1:subasin) = irri_supply



         ! doacudes weitere Option


















END SUBROUTINE
