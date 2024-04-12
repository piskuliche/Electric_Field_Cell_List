! ****************************************************************************************************
! CellList Module
! 
! This module contains the subroutines that are used to create the cell list and the linked list
! that are used to calculate the internal distance of the particles. The module also contains the
! function that calculates the distance between two particles.
!
! ****************************************************************************************************


MODULE CellList

    IMPLICIT NONE

    CONTAINS

    SUBROUTINE TestMap(i, cell_map)
        IMPLICIT NONE
        INTEGER :: i
        INTEGER, DIMENSION(:,:) :: cell_map

        INTEGER l, m, n

        DO l = -2,2
            WRITE(*,*) "l *", i, l, i+l+2
            WRITE(*,*) cell_map(1,i+l+2)
        ENDDO
        WRITE(*,*) "**"
        DO m=1, 34
            WRITE(*,*) cell_map(1,m)
        ENDDO 

    END SUBROUTINE TestMap

    INCLUDE "construct_grid.f90"
    INCLUDE "linked_list_build.f90"
    INCLUDE "internal_distance.f90"
    INCLUDE "functions.f90"


END MODULE CellList