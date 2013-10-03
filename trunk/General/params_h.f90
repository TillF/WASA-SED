! Till: separated modules to allow easier use by eclipse
! 2013-10-02 

! Till: computationally irrelevant: relocated vars
! 2012-09-21 

! Till: computationally irrelevant: minor changes to improve compiler compatibility
! 2011-04-29

! Till: computationally irrelevant: added program version information to parameter.out
! 2009-06-17

! Till: added output for River_Sediment_Storage.out
! 2008-11-13

! 2008-07-11
! Till: implemented optional pre-specified sediment outflow of selected subbasins

! 2008-07-03
! Till: implemented optional pre-specified outflow of selected subbasins

! 2008-01-31 
! Till: included variables for loading and saving of initial conditions
! included variables for controling evaporation calculation

! 2007-10-18
! Till: increased length of pathp variable to 160

! 2007-06-04
! Till: added flag f_tc_theta

! 2007-04-28
! Till: added flag f_tc_theta

! 2007-01-10
! Till: renamed f_deep_gw to f_deep_gw_recharge, added f_deep_gw_discharge

!Till: added variable f_deep_gw as flag for fileoutput

!

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
module params_h
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
save
! TOTAL NUMBER OF SUB-BASINS IN STUDY AREA (e.g. 10 sub-basins)
INTEGER :: subasin
! TOTAL NUMBER OF SOTER UNITS IN STUDY AREA (e.g. 321 units)
INTEGER :: nsoter
! TOTAL NUMBER OF TERRAIN COMPONENTS IN STUDY AREA (e.g. 515 components)
INTEGER :: nterrain
! TOTAL NUMBER OF SOIL COMPONENTS IN STUDY AREA (e.g. 72 components)
INTEGER :: nsoil
! TOTAL NUMBER OF VEGETATION UNITS IN STUDY AREA (e.g. 34 units)
INTEGER :: nveg
! TOTAL NUMBER OF CELL/SOTER UNIT/TERRAIN COMPONENT COMBINATIONS (e.g. 49 combinations)
INTEGER :: ntcinst
! MAXIMUM NUMBER OF SOTER UNITS IN CELLS (7)
INTEGER :: maxsoter
!PARAMETER (maxsoter=7)
! MAXIMUM NUMBER OF TERRAIN COMPONENTS IN SOTER UNIT (3)
INTEGER :: maxterrain
!PARAMETER (maxterrain=3)
! MAXIMUM NUMBER OF SOIL COMPONENTS IN TERRAIN COMPONENTS (28)
INTEGER :: maxsoil
!PARAMETER (maxsoil=28)
! MAXIMUM NUMBER OF HORIZONS IN SOIL COMPONENTS (8)
INTEGER :: maxhori
!PARAMETER (maxhori=8)
! TOTAL NUMBER OF TRANSPOSITIONS BETWEEN Sub-basins

INTEGER :: ntrans
!PARAMETER (ntrans=2)
!common / basin_parameter / subasin, ntcinst, nsoter, nterrain, nsoil, nveg
end module params_h
