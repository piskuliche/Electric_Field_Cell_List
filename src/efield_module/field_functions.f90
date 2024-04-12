FUNCTION OH_Vector(ra, rb) RESULT(eOH)
! *********************************************************************
! This subroutine calculates the OH vector, and normalizes it
!
! Inputs:
!   - ra: Coordinates of the H atom
!   - rb: Coordinates of the O atom
!
! Outputs:
!   - eOH: OH vector
!
! *********************************************************************

    IMPLICIT NONE

    REAL, DIMENSION(3), INTENT(IN) :: ra, rb
    REAL, DIMENSION(3) :: eOH

    REAL :: norm
    INTEGER :: k

    eOH = 0.0
    norm = 0.0
    DO k=1,3
        eOH(k) = ra(k) - rb(k)
        norm = norm + eOH(k)**2
    ENDDO ! k
    norm = SQRT(norm)

    DO k=1,3
        eOH(k) = eOH(k) / norm
    ENDDO ! k

END FUNCTION OH_Vector

SUBROUTINE BondVector(rA, rB, eOH)
    REAL, DIMENSION(:,:) :: rA, rB
    REAL, DIMENSION(:,:) :: eOH

    INTEGER :: i

    eOH = rA - rB
    DO i=1, SIZE(eOH, DIM=1)
        eOH(i,:) = eOH(i,:) / NORM2(eOH(i,:))
    ENDDO
END SUBROUTINE BondVector


FUNCTION PBC_Dist(input_dr, box_lengths) RESULT(output_dr)

    IMPLICIT NONE
    REAL, DIMENSION(3), INTENT(IN) :: input_dr, box_lengths
    REAL, DIMENSION(3) :: output_dr

    output_dr(:) = input_dr(:) - box_lengths(:)*ANINT(input_dr(:)/box_lengths(:))

END FUNCTION PBC_Dist

Function PBC_dr(pos1, pos2, box_lengths) RESULT(output)
    IMPLICIT NONE
    REAL, DIMENSION(3) :: pos1, pos2, box_lengths, dr
    REAL, DIMENSION(4) :: output

    dr = 0.0; output = 0.0
    dr = pos1 - pos2
    dr = dr - box_lengths*ANINT(dr/box_lengths)

    output(4) = SUM(dr**2)
    output(1:3) = dr

END FUNCTION PBC_dr

Function PBC_dr_new(pos1, pos2, box_lengths) RESULT(output)
    IMPLICIT NONE
    REAL, DIMENSION(3) :: pos1, pos2, box_lengths, dr
    REAL, DIMENSION(4) :: output

    INTEGER :: k
    REAL :: dist

    dr = 0.0; output = 0.0
    dist = 0
    DO k=1,3
        dr(k) = pos1(k) - box_lengths(k)*ANINT((pos1(k)-pos2(k))/box_lengths(k))
        dist = dist + (pos2(k)-dr(k))**2.
    ENDDO

    output(4) = dist
    output(1:3) = dr

END FUNCTION PBC_dr_new
