
FUNCTION PBC_Wrap(position, box_lengths) RESULT(wrapped)
! This function wraps coordinates back into the box
    IMPLICIT NONE
    REAL, INTENT(IN) :: position, box_lengths
    REAL :: wrapped


    wrapped = 0.0


    IF(position < 0) THEN
        wrapped = position + box_lengths*CEILING(ABS(position)/box_lengths)
    ELSE IF(position > box_lengths) THEN
        wrapped = position - box_lengths*CEILING(ABS(position)/box_lengths)
    ELSE
        wrapped = position
    ENDIF

END FUNCTION PBC_Wrap


FUNCTION Assign_To_Cell_Grid(position, cell_nbins, box_lengths) RESULT(location_on_grid)
! This function takes position, a number of bins in each direction, and box lengths and then
! Assigns them to a specific bin.
! 
!        --- --- ---
!       | 1 | 2 | 3 |
!        --- --- ---
!       | 4 | 5 | 6 |
!        --- --- ---
!       | 7 | 8 | 9 |
!        --- --- ---
!       
!       |---L=10----|
! 
!       In this case, the point (1,1.5) would fall into cell 1
!       whereas the point (4,9) would fall into cell 6 
!       (note that in reality, cell 6 would be (3,2) in the 2d notation)
! 
! Parameters:
! -----------
! position : REAL, DIMENSION(3)
!       The x,y,z coordinates
! cell_nbins: INTEGER, DIMENSION(3)
!       The # of bins in each direction
! box_lengths: REAL, DIMENSION(3)
!       The box dimensions
! 
! Returns:
! --------
! location_on_grid : INTEGER, DIMENSION(3)
!       The bin that the positions belong to.
!
    IMPLICIT NONE

    REAL, DIMENSION(3), INTENT(IN) :: position, box_lengths
    INTEGER, DIMENSION(3), INTENT(IN) :: cell_nbins

    INTEGER :: i
    INTEGER, DIMENSION(3) :: location_on_grid
    REAL, DIMENSION(3) :: temp_pos

    location_on_grid = 0

    DO i=1,3
        temp_pos(i) = PBC_Wrap(position(i), box_lengths(i))
        location_on_grid(i) = MIN(cell_nbins(i), MAX(1, CEILING(temp_pos(i)/box_lengths(i)*REAL(cell_nbins(i)))))
    ENDDO
END FUNCTION Assign_To_Cell_Grid


FUNCTION Electric_Field_Component(charge, distance_vector, r_distance) RESULT(EField)
    IMPLICIT NONE

    REAL,               INTENT(IN) :: charge, r_distance
    REAL, DIMENSION(3), INTENT(IN) :: distance_vector
    REAL, DIMENSION(3) :: Efield
    
    Efield = charge * distance_vector / (r_distance**3.)

END FUNCTION Electric_Field_Component

FUNCTION Field_Contribution(q, ra, vector, dist) RESULT (efield)
! *********************************************************************
! This subroutine calculates the contribution to the electric field
! from a single atom
!
! Inputs:
!   - q: Charge of the atom
!   - ra: Coordinates of the atom
!   - vector: Vector pointing from the atom to the OH bond
!   - dist: Distance between the atom and the OH bond
!
! Outputs:
!   - efield: Electric field contribution from the atom
!
! *********************************************************************
    
    IMPLICIT NONE

    REAL, INTENT(IN) :: q, dist
    REAL, DIMENSION(3), INTENT(IN) :: ra, vector
    REAL, DIMENSION(3) :: efield

    INTEGER :: k
    efield = 0
    DO k=1,3
        efield(k) = efield(k) + q * (ra(k) - vector(k))/(dist**3)
    ENDDO

END Function



SUBROUTINE PBC_DIST(rA, rB, box_lengths, dist_vec, dist)
    IMPLICIT NONE
    REAL, DIMENSION(3) :: rA, rB, box_lengths
    REAL, DIMENSION(3) :: dist_vec
    REAL :: dist

    dist_vec = rA(:) - rB(:)
    dist_vec = dist_vec - box_lengths*ANINT(dist_vec/box_lengths)
    dist = SQRT(SUM(dist_vec))
END SUBROUTINE

