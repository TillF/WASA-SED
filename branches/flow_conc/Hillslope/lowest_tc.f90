  
  
  IF (templat > 0.) THEN					!Till: is there subsurface runoff to be redistributed?
    watbal=watbal+templat/(tcarea2*1.e3)	!ii tcarea gleich am Anfang mit 1000 multiplizieren
    IF (nbr_svc(tcid_instance2) == 0) THEN	!Till: if whole TC consists of rock...
      q_surf_out=q_surf_out+templat				!Till: subsurface inflow becomes returnflow and leaves the TC as surface runoff
    ELSE IF (.NOT. dolatsc) THEN				!Till: if no lateral flowdistribution...
      q_surf_out=q_surf_out+templat*rocky(tcid_instance2)	!Till: surface runoff is increased according to fraction of rocky parts
    ELSE
      q_surf_out=q_surf_out+templat*rocky(tcid_instance2)*rocky(tcid_instance2) !Till: surface runoff is increased according to fraction of rocky parts, but reduced due to lateral redistribution
    END IF
    
IF (tc_counter2 == nbrterrain(i_lu)) THEN	!Till: treat the very lowest TC in a different way    
! check if alluvial soils exist and sum up their fractions !ii: this is a static property and could by computed once
    merkalluv(:)=0. !fraction of SVC covered by alluvial soil (i.e. 0 for non-alluvials, otherwise fraction of SVC)
    
    merkalluv(1:nbr_svc(tcid_instance2)) = &
        alluvial_flag(id_soil_intern(1:nbr_svc(tcid_instance2), tcid_instance2)) *&
        frac_svc     (               1:nbr_svc(tcid_instance2), tcid_instance2)
        
! total fraction of alluvial soils in current TC
    allalluv=sum(merkalluv(1:nbr_svc(tcid_instance2))) 
    
! if any alluvial soils occur in TC
    !debugremove: this eliminates the excess!
    IF (allalluv > 0.) THEN 
      
!  maximum lateral inflow into this SVC
      DO i=1,nbr_svc(tcid_instance2)
        IF (merkalluv(i) > 0.) THEN	!Till: if the current SVC has an alluvial soil, treat it
		  soilid=id_soil_intern(i,tcid_instance2)	
          remainlat(:)=0.
          testi=0
          
!  maximum lateral inflow into this SVC
          h=size(remainlat)
		  remainlat(1:h)=(latred(tcid_instance2, 1:h)* merkalluv(i)/allalluv)/ (tcarea2*1.e3*frac_svc(i,tcid_instance2)) !ii vectorize
          				!Till: all potential subsurface flow is distributed among alluvial SVCs [mm] (still "remains to be distributed")
                    
!  distribute inflow among horizons
          DO h=1,nbrhori(soilid)
            lath=min(size(remainlat),INT(sum(horiz_thickness(tcid_instance2,i,1:h))/500.))	!Till: compute number of "exchange horizons" to be used from remainlat to be fed into the current profile
            IF (lath > 0) THEN
              temp2=sum(remainlat(1:lath))			!Till: max amount of water that is available to be routed into this horizon [mm]
              
!  update soil moisture of this horizon
!  maximum inflow may be limited by hydraulic conductivity
              IF (temp2 > 0.) THEN
				!Till: compute maximum inflow according to "exchange surface"
				!Andreas: has to be in unit mm, correct calculation similar to lateral outflow calculation
			    temp3=k_sat(soilid,h)/dt_per_day*(1.-coarse(soilid,h))* &           ! Q_li 
				    horiz_thickness(tcid_instance2,i,h)*slope(id_tc_type2)/100.  &
					/  slength(i_lu) / fracterrain(id_tc_type2)/1000.  

                temp3=MIN(temp3,temp2)					!Till: actual amount of water going into current horizon
                horithact(tcid_instance2,i,h)=horithact(tcid_instance2,i,h)+temp3

				remainlat(1:lath)=0.
                IF (temp3 < temp2) THEN					!Till: if inflow is limited by conductivity, the rest is stored for further redistribution
                  remainlat(lath)=temp2-temp3
                END IF
!  if horizon is saturated, remaining lateral inflow may be added to
!  lower horizon
                temp5=horithact(tcid_instance2,i,h) - thetas(soilid,h)* horiz_thickness(tcid_instance2,i,h)	!Till: compute if watercontent exceeds saturation
					  
				IF (temp5 > 0. ) THEN	!Till: if water content exceeds saturation...
				  remainlat(lath)=remainlat(lath)+ temp5				!Till: excess cannot flow into this horizon
                  horithact(tcid_instance2,i,h)=thetas(soilid,h)* horiz_thickness(tcid_instance2,i,h)	!Till: horizon is saturated
                END IF
              END IF
            END IF
! enddo for all horizons
          END DO
          
! check if remaining lateral subsurface inflow
! may then be distributed also in higher horizons which are below river bed
! repeat procedure for all deeper horizons
! remaining latflow is river runoff
          temp2=sum(remainlat(:))
          IF (temp2 > 0.) THEN
            
            lath=maxhori !default: use all horizons
            DO h=1,maxhori !reduce number of exchange horizons according to depth of riverbed
              IF (sum(horiz_thickness(tcid_instance2,i,1:h)) > riverbed(i_lu)) THEN
                lath=h
				exit
              END IF
            END DO
!			DO WHILE (lath == 0)						!ii: Schleife restrukturieren: for cycle exit continue 
!              IF (sum(horiz_thickness(tcid_instance2,i,1:h)) > riverbed(i_lu)) THEN
!                lath=h
!              END IF
!              h=h+1
!              IF (h >= maxhori) THEN
!                lath=maxhori
!              END IF
!            END DO
            DO h=lath,nbrhori(soilid)
              temp2=sum(remainlat(:))
              
!  update soil moisture of this horizon
!  maximum inflow may be limited by hydraulic conductivity
              IF (temp2 > 0.) THEN

				!Till: compute maximum inflow according to "exchange surface"
				!Andreas: has to be in unit mm, correct calculation similar to lateral outflow calculation
			    temp3=k_sat(soilid,h)/dt_per_day*(1.-coarse(soilid,h))* &           ! Q_li 
				    horiz_thickness(tcid_instance2,i,h)*slope(id_tc_type2)/100.  &
					/  slength(i_lu) / fracterrain(id_tc_type2)/1000.  

				temp3=MIN(temp3,temp2)
                horithact(tcid_instance2,i,h)=horithact(tcid_instance2,i,h)+temp3

				remainlat(:)=0.
                IF (temp3 < temp2) THEN
                  remainlat(1)=temp2-temp3
                END IF
                
!  if horizon is saturated, remaining lateral inflow may be added to
!  lower horizon
                temp2=thetas(soilid,h)*horiz_thickness(tcid_instance2,i,h)	!Till: compute saturated water content [mm]
				
				IF (horithact(tcid_instance2,i,h) > temp2) THEN
                  remainlat(1)=remainlat(1)+ (horithact(tcid_instance2,i,h)- temp2)
                  horithact(tcid_instance2,i,h)=temp2
                END IF
              END IF
              
! enddo for all horizons
            END DO
          END IF
          
!  if not all lateral inflow to this SVC can be distributed among the horizons,
!  the remaining inflow is assumed to become subsurface runoff from this
!  terrain component
!  (this is e.g. also the case if the current soil profile is very shallow)
          temp2=sum(remainlat(:))
          IF (temp2 > 0.) THEN
            q_sub_out=q_sub_out+temp2*tcarea2*frac_svc(i,tcid_instance2)*1.e3	!Till: excess water is temporarily assumed subsurface runoff from TC [m3]
            watbal=watbal-temp2*frac_svc(i,tcid_instance2)						!Till: ATTENTION: this water is redistributed among non-alluvial SVCs later
          END IF
          
          
!  Calculate new fraction of saturation of current SVC
!  frac_sat value is relative to TC, not to SVC !!
          
          tempth=sum(horithact(tcid_instance2,i,:))
          
! no part of SVC is saturated
          
          IF (tempth <= tctheta_s(1,1,tcid_instance2,i)) THEN
            frac_sat(tcid_instance2,i)=0.
            
! soil component is partly or completely saturated
          ELSE
            
            tempalt=0.
            testi=0
            
            DO j=1,4
              IF (testi == 0) THEN
                
                tempx=tempalt+(tctheta_s(j+1,1,tcid_instance2,i)-  &
                    tctheta_s(j,1,tcid_instance2,i))* (tctheta_s(j+1,2,tcid_instance2,i)-  &
                    tctheta_s(j,2,tcid_instance2,i))/2.+ (tctheta_s(j+1,1,tcid_instance2,i)-  &
                    tctheta_s(j,1,tcid_instance2,i))* (tctheta_s(5,2,tcid_instance2,i)-  &
                    (tctheta_s(j+1,2,tcid_instance2,i)- tctheta_s(1,2,tcid_instance2,i)))
                
                IF (tempth > tempalt .AND. tempth < tempx) THEN
                  
                  frac_sat(tcid_instance2,i)=((tctheta_s(j,2,tcid_instance2,i)-  &
                      tctheta_s(1,2,tcid_instance2,i))+  &
                      (tempth-tempalt)/(tempx-tempalt)*  &
                      (tctheta_s(j+1,2,tcid_instance2,i)- tctheta_s(j,2,tcid_instance2,i)))*  &
                      frac_svc(i,tcid_instance2)
                  testi=1
                END IF
                tempalt=tempx
              END IF
            END DO
          END IF ! end of update saturated fraction of SVC

        END IF	! endif alluvial SVC

      END DO	! enddo over all SVCs


	  if ((allalluv <= 1.) .AND. (q_sub_out > 0.)) then		!Till: if there are also non-alluvial soils AND not all subsurface inflow has been absorbed by the alluvials before
	    !Till: distribute among exchange horizons with the same proportions as before the alluvial influence
		latred(tcid_instance2,:)=latred(tcid_instance2,:)/sum(latred(tcid_instance2,:))*q_sub_out
		watbal=templat/(tcarea2*1.e3)							!reset water balance
		q_sub_out=0.											!Till: no surface runoff yet - may result later after treatment of non-alluvial SVCs
	  else
		latred(tcid_instance2,:)=0.								!Till: nothing more to do
	  end if

	END IF ! endif alluvial soils occur in this most downslope TC
ELSE
 allalluv=0. !set alluvial fraction to 0 for non-lowermost TCs
END IF ! end if this is the lowermost TC
	

      IF ((sum(latred(tcid_instance2,:)) > 0.) .OR. & !Till: if there is excess water not infiltrated into alluvial soils
             (tc_counter2 /= nbrterrain(i_lu))) THEN   !ii: this condition is probably obsolete, check
     

		     DO i=1,nbr_svc(tcid_instance2) !regular, non-alluvial SVCs
        
			soilid=id_soil_intern(i,tcid_instance2)
			remainlat(:)=0.
			testi=0

			IF ((tc_counter2 == nbrterrain(i_lu)) .AND. (alluvial_flag(soilid)==1)) THEN	!Till: check if this soil is an alluvial soil
				cycle		!Till: alluvial soils already have been treated, so skip them
			END IF
        
     
	!  maximum lateral inflow into this SVC
			DO h=1,size(remainlat)
			  remainlat(h)=(latred(tcid_instance2,h)*frac_svc(i,tcid_instance2)/(1.-allalluv)) / & !Till: distribute subsurface flow according to fraction among non-alluvial soils
				  (tcarea2*1.e3*frac_svc(i,tcid_instance2))
			END DO
        
        
	!  check each horizon of SVC if there is lateral inflow (from soil depth
	!  above this horizon)
			DO h=1,nbrhori(soilid)
			  lath=INT(sum(horiz_thickness(tcid_instance2,i,1:h))/500.)
			  IF (lath > 0) THEN
				temp2=sum(remainlat(1:lath))
            
	!  update soil moisture of this horizon
	!  maximum inflow may be limited by hydraulic conductivity
				IF (temp2 > 0.) THEN
				  !Till: compute maximum inflow according to "exchange surface"
                                  IF (tc_counter2 == nbrterrain(i_lu)) THEN	!Till: treat the very lowest TC in a different way
                                       !ii: this different treatment is most likely an error. I assume Andreas' bugfix has only be applied to non-lowermost TCs. Check this!
                                      temp3=(k_sat(soilid,h)/dt_per_day*(1.-coarse(soilid,h))*  &
					horiz_thickness(tcid_instance2,i,h)*slope(id_tc_type2)/100)* &
					area(i_subbas2)*frac_lu(lu_counter2,i_subbas2) * frac_svc(i,tcid_instance2) /&
					slength(i_lu)		!subsurface inflow into horizon(s) admitted due to conductivity
                                  ELSE
                                      !Andreas: has to be in unit mm, correct calculation similar to lateral outflow calculation
                                      temp3=k_sat(soilid,h)/dt_per_day*(1.-coarse(soilid,h))* &           ! Q_li 
                                      horiz_thickness(tcid_instance2,i,h)*slope(id_tc_type2)/100.  &
                                      /  slength(i_lu) / fracterrain(id_tc_type2)/1000.  
                                  END IF


				  temp3=MIN(temp3,temp2)
				  horithact(tcid_instance2,i,h)=horithact(tcid_instance2,i,h)+temp3

				  remainlat(1:lath)=0.
				  IF (temp3 < temp2) THEN
					remainlat(lath)=temp2-temp3
				  END IF
              
	!  if horizon is saturated, remaining lateral inflow may be added to
	!  lower horizon
				  IF (horithact(tcid_instance2,i,h) > thetas(soilid,h)*  &
						horiz_thickness(tcid_instance2,i,h)) THEN
					remainlat(lath)=remainlat(lath)+ (horithact(tcid_instance2,i,h)-  &
						thetas(soilid,h)*horiz_thickness(tcid_instance2,i,h))
					horithact(tcid_instance2,i,h)=thetas(soilid,h)* horiz_thickness(tcid_instance2,i,h)
					

				  END IF
				END IF
			  END IF
	!  enddo horizons
			END DO
        
	!  if not all lateral inflow to this SVC can be distributed among the horizons,
	!  the remaining inflow is assumed to become subsurface runoff from this
	!  terrain component
	!  (this is e.g. also the case if the current soil profile is very shallow)
			temp2=sum(remainlat(:))
			IF (temp2 > 0.) THEN
			  q_sub_out=q_sub_out+temp2*tcarea2*frac_svc(i,tcid_instance2)*1.e3
			  watbal=watbal-temp2*frac_svc(i,tcid_instance2)
			END IF
        
        
	!  Calculate new fraction of saturation of current SVC
	!  frac_sat value is relative to TC, not to SVC !!
        
			tempth=sum(horithact(tcid_instance2,i,:))
        
	! no part of SVC is saturated
        
			IF (tempth <= tctheta_s(1,1,tcid_instance2,i)) THEN
			  frac_sat(tcid_instance2,i)=0.
          
	! soil component is partly or completely saturated
			ELSE
          
			  tempalt=0.
			  testi=0
          
			  DO j=1,4
				IF (testi == 0) THEN
              
				  tempx=tempalt+(tctheta_s(j+1,1,tcid_instance2,i)-  &
					  tctheta_s(j,1,tcid_instance2,i))* (tctheta_s(j+1,2,tcid_instance2,i)-  &
					  tctheta_s(j,2,tcid_instance2,i))/2.+ (tctheta_s(j+1,1,tcid_instance2,i)-  &
					  tctheta_s(j,1,tcid_instance2,i))* (tctheta_s(5,2,tcid_instance2,i)-  &
					  (tctheta_s(j+1,2,tcid_instance2,i)- tctheta_s(1,2,tcid_instance2,i)))
              
				  IF (tempth > tempalt .AND. tempth < tempx) THEN
                
					frac_sat(tcid_instance2,i)=((tctheta_s(j,2,tcid_instance2,i)-  &
						tctheta_s(1,2,tcid_instance2,i))+ (tempth-tempalt)/(tempx-tempalt)*  &
						(tctheta_s(j+1,2,tcid_instance2,i)- tctheta_s(j,2,tcid_instance2,i)))*  &
						frac_svc(i,tcid_instance2)
					testi=1
				  END IF
				  tempalt=tempx
				END IF
			  END DO
			END IF
	! end of update saturated fraction of SVC
        
	! end of loop for all non-alluvial SVCs
		  END DO

		END IF      ! end if lateral inflow > 0 (after alluvials)

  END IF		! end if (lateral inflow templat > 0
  latred(tcid_instance2,:)=0.				!Till: all lateral subsurface flow has been treated
