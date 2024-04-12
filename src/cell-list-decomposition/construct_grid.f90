! **************************************************************************************************** 
! This subroutine constructs a grid of cells for the calculation.
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
! cell_map :: INTEGER, DIMENSION(:,:), ALLOCATABLE
!       The map of the cells. The first index is the direction (1, 2, or 3) and
!         the second index is the cell number.
!
! ****************************************************************************************************

SUBROUTINE Construct_Grid(box_lengths, r_cutoff, cell_nbins, cell_map)

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

END SUBROUTINE Construct_Grid