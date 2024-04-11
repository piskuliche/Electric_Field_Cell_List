MODULE CellList

    IMPLICIT NONE

    CONTAINS

    SUBROUTINE Construct_Grid(box_lengths, r_cutoff, cell_nbins, cell_map)
    ! This subroutine constructs the cell_nbins array.
    ! 
    ! The way this works is that if you have a box length of 30.0 angstroms
    ! and the cutoff is 7 angstroms. 30/7 = 4.29 bins, which rounds down to 4.0 bins.
    !
    ! 
    !
    ! Parameters:
    ! -----------
    ! box_lengths :: REAL, DIMENSION(3)
    !       The box dimensions in angstroms
    ! r_cutoff :: REAL
    !       The width of the bins in each direction.
    ! cell_nbins :: INTEGER, DIMENSION(3)
    !       The number of bins in each direction.
    !
    !
        IMPLICIT NONE
        REAL,    DIMENSION(3),      INTENT(IN) :: box_lengths
        REAL,                       INTENT(IN) :: r_cutoff
        INTEGER, DIMENSION(3),      INTENT(OUT):: cell_nbins
        INTEGER, DIMENSION(:,:), ALLOCATABLE,   INTENT(OUT) :: cell_map

        INTEGER :: i, j, max_nbins

        DO i=1, 3
            cell_nbins(i) = FLOOR(box_lengths(i)/(r_cutoff/2.0))
        ENDDO 
        max_nbins = MAXVAL(cell_nbins)+4

        IF (.NOT. ALLOCATED(cell_map)) THEN
            ALLOCATE(cell_map(3, max_nbins))
        ENDIF

        cell_map=0
        DO i = 1, 3
            cell_map(i,1)  = cell_nbins(i) - 1
            cell_map(i,2)  = cell_nbins(i)
            DO j=1, cell_nbins(i)
                cell_map(i,j+2) = j
            ENDDO
            cell_map(i,cell_nbins(i)+3) = 1
            cell_map(i,cell_nbins(i)+4) = 2
        ENDDO

        !CALL TestMap(29, cell_map)

    END SUBROUTINE Construct_Grid

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

    SUBROUTINE Linked_List_Build(positions,  box_lengths, r_cutoff, cell_nbins, ll_head, cell_list, cell_map)
    ! This subroutine builds a linked list based on the cell dimensions. It assigns all
    ! the positions to locations on the grid, and then puts these into the cell list.
    ! 
    ! Parameters:
    ! -----------
    ! positions :: REAL, DIMENSION(:,:)
    !       Coordinates of all the atoms for the build
    ! cell_nbins ::  INTEGER, DIMENSION(:)
    !       The # of bins in each direction
    ! box_lengths :: REAL, DIMENSION(3)
    !       The box dimensions in each direction
    ! ll_head :: INTEGER, DIMENSION(:,:,:)
    !       The location of the head of the linked list for each
    !       bin.
    ! cell_list :: INTEGER, DIMENSION(:)
    !       The linked list.
    !
    !

        IMPLICIT NONE

        REAL, DIMENSION(:,:),       INTENT(IN) :: positions
        REAL, DIMENSION(3),         INTENT(IN) :: box_lengths
        REAL,                       INTENT(IN) :: r_cutoff
        INTEGER, DIMENSION(3),      INTENT(IN) :: cell_nbins
        INTEGER, DIMENSION(:,:,:),  INTENT(OUT) :: ll_head
        INTEGER, DIMENSION(:),      INTENT(OUT) :: cell_list
        INTEGER, DIMENSION(:,:),    INTENT(OUT) :: cell_map

        INTEGER :: i
        INTEGER, DIMENSION(3) :: location_on_grid
        

        ll_head = 0; cell_list = 0

        DO i=1, SIZE(positions,1)
            ! Assign the location on the grid
            location_on_grid(:) = Assign_To_Cell_Grid(positions(i,:), cell_nbins, box_lengths)
            ! Put that information into the cell_list
            cell_list(i) = ll_head(location_on_grid(1), location_on_grid(2), location_on_grid(3))
            ! Put the head of that location at i.
            ll_head(location_on_grid(1), location_on_grid(2), location_on_grid(3)) = i
        ENDDO

    END SUBROUTINE Linked_List_Build

    INCLUDE "internal_distance.f90"
    INCLUDE "functions.f90"


END MODULE CellList