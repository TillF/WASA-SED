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
! MAXIMUM NUMBER OF TERRAIN COMPONENTS IN SOTER UNIT (3)
INTEGER :: maxterrain
! MAXIMUM NUMBER OF SOIL COMPONENTS IN TERRAIN COMPONENTS (28)
INTEGER :: maxsoil
! MAXIMUM NUMBER OF HORIZONS IN SOIL COMPONENTS (8)
INTEGER :: maxhori
! TOTAL NUMBER OF TRANSPOSITIONS BETWEEN Sub-basins
INTEGER :: ntrans

INTEGER, parameter :: output_int=366 !interval of file output [days] (still to be implemented fully)
end module params_h
