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
	REAL :: rate(subasin + 1)
    INTEGER, allocatable :: receiver_basins(:), sub_exists(:)


    irri_supply = 0.0
    irri_abstraction = 0.0
    rate = 0.0

    !-- Abstraction from groundwater, Loop through all subbasins and check if they're included in irri.dat-----------------------
	DO  sb_counter = 1, subasin ! Subbasin Loop for groundwater abstraction
        abstraction_requested = 0.0
        abstraction_available = 0.0

        sub_exists = pack(sub_receiver, MASK=sub_source == sb_counter .AND. irri_source == "groundwater") !If
        IF (SIZE(sub_exists) == 0) THEN
            cycle
        END IF

        ! get the indices of irri.dat for all subbasins that receive irrigation water from current subbasin and are irrigated seasonal
        receiver_basins = pack([(i, i=1, nbr_irri_records)], MASK =sub_source == sb_counter .AND. irri_source == "groundwater" .AND. irri_rule == "seasonal")

        ! get total amount of seasonal and fixed abstracted water in current subbasin
        IF (SIZE(receiver_basins) /= 0) THEN
            rate = calc_seasonality2(sub_receiver(receiver_basins(1)), t, d, seasonality_irri_gw, irri_rate_gw) !this function needs one seasonally irrigated subbasinId as first argument and then strangely computes the values for all the seasonal and fixed receivers
            abstraction_requested = sum(rate)
        ELSE IF (SIZE(receiver_basins) == 0) THEN
            receiver_basins = pack([(i, i=1, nbr_irri_records)], MASK =sub_source == sb_counter .AND. irri_source == "groundwater" .AND. irri_rule == "fixed" .AND. sub_receiver /= 9999)
            rate(sub_receiver(receiver_basins)) = irri_rate_gw(sub_receiver(receiver_basins),1)
            abstraction_requested = sum(irri_rate_gw(sub_receiver(receiver_basins),1)) + irri_rate_gw(subasin + 1,1) ! <- can't index with external. External is always the last row in irri_rate arrays
        END IF

        !Einfügen abstraction requested von soil moisture thresh

        all_request = abstraction_requested !add water from sm_thresh to total request

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

        irri_supply = irri_supply + rate(1:subasin)/all_request * abstraction_requested * loss_gw  !This vector contains the total irrigation ammount for every subbasin

    END DO ! Subbasin Loop for groundwater abstraction

    !--Add water from external sources (9999) to irri_supply ----------------------------------------------
    rate = 0.0
    receiver_basins = pack([(i, i=1, nbr_irri_records)], MASK=sub_source == 9999) !get indices of all external sources

    IF (SIZE(receiver_basins) > 0)  THEN !If there are external sources
        receiver_basins = pack([(i, i=1, nbr_irri_records)], MASK=sub_source == 9999 .AND. irri_rule == "seasonal") !get indices of all external sources

        IF (SIZE(receiver_basins) /= 0) THEN
            rate = calc_seasonality2(sub_receiver(receiver_basins(1)), t, d, seasonality_irri_ext, irri_rate_ext) !this function needs one seasonally irrigated subbasinId as first argument and then strangely computes the values for all the seasonal and fixed receivers
            irri_supply = irri_supply + rate(1:subasin) * loss_ext
        ELSE IF (SIZE(receiver_basins) == 0) THEN
            receiver_basins = pack([(i, i=1, nbr_irri_records)], MASK =sub_source == sb_counter .AND. irri_source == "groundwater" .AND. irri_rule == "fixed" .AND. sub_receiver /= 9999)
            rate(sub_receiver(receiver_basins)) = irri_rate_ext(sub_receiver(receiver_basins),1)
            irri_supply = irri_supply + rate(1:subasin) * loss_ext
        END IF
    END IF

    !-- Abstraction from river, Loop through all subbasins and check if they're included in irri.dat-----------------------
    rate = 0.0

    DO sb_counter = 1, subasin ! Subbasin Loop for river
        abstraction_requested = 0.0
        abstraction_available = 0.0

        sub_exists = pack(sub_receiver, MASK=sub_source == sb_counter .AND. irri_source == "river")
        IF (SIZE(sub_exists) == 0) THEN
            cycle
        END IF

        ! get the indices of irri.dat for all subbasins that receive irrigation water from current subbasin and are irrigated seasonal
        receiver_basins = pack([(i, i=1, nbr_irri_records)], MASK =sub_source == sb_counter .AND. irri_source == "river" .AND. irri_rule == "seasonal")

        ! get total amount of seasonal and fixed abstracted water in current subbasin
        IF (SIZE(receiver_basins) /= 0) THEN
            rate = calc_seasonality2(sub_receiver(receiver_basins(1)), t, d, seasonality_irri_riv, irri_rate_riv) !this function needs one seasonally irrigated subbasinId as first argument and then strangely computes the values for all the seasonal and fixed receivers
            abstraction_requested = sum(rate)
        ELSE IF (SIZE(receiver_basins) == 0) THEN
            receiver_basins = pack([(i, i=1, nbr_irri_records)], MASK =sub_source == sb_counter .AND. irri_source == "river" .AND. irri_rule == "fixed" .AND. sub_receiver /= 9999)
            rate(sub_receiver(receiver_basins)) = irri_rate_riv(sub_receiver(receiver_basins),1)
            abstraction_requested = sum(irri_rate_riv(sub_receiver(receiver_basins),1)) + irri_rate_riv(subasin + 1,1) ! <- can't index with external. External is always the last row in irri_rate arrays
        END IF

        !Einfügen abstraction requested von soil moisture thresh

        all_request = abstraction_requested !add water from sm_thresh to total request

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
        irri_supply = irri_supply + rate(1:subasin)/all_request * abstraction_requested * loss_riv !This vector contains the total irrigation ammount for every subbasin

    END DO ! Subbasin Loop for river abstraction





     !irri_supply und irri_abstraction in irri_abstraction_record und irri_receiver_record schreiben -> Ausgabedatei um Zwischenstände abzuspeichern

     irri_abstraction_record(d,nt,1:subasin) = irri_abstraction
     irri_supply_record(d,nt,1:subasin) = irri_supply



         ! doacudes weitere Option


















END SUBROUTINE
