
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