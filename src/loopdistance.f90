MODULE DistanceModule
    IMPLICIT NONE
    INTEGER, parameter :: max_bins = 500
    CONTAINS

    SUBROUTINE LoopDistance(positions, r_cutoff, box_lengths, distance, local_count)
        REAL, DIMENSION(:,:), INTENT(IN) :: positions
        REAL,                 INTENT(IN) :: r_cutoff
        REAL, DIMENSION(3),   INTENT(IN) :: box_lengths
        REAL, DIMENSION(:,:), INTENT(OUT) :: distance
        INTEGER, DIMENSION(:), INTENT(OUT) :: local_count

        REAL, DIMENSION(3) :: dr

        INTEGER :: i, j
        INTEGER :: number_of_atoms
        INTEGER :: count

        count = 0
        distance = 0.0
        local_count = 0

        number_of_atoms = SIZE(positions, 1)

        DO i = 1, number_of_atoms
            DO j = 1, number_of_atoms

                dr = positions(i,:) - positions(j,:)
                dr = dr - box_lengths(:)*ANINT(dr(:)/box_lengths(:))
                distance(i, j) = SQRT(SUM(dr(:)**2))
                IF (distance(i,j) < r_cutoff) THEN
                    count = count + 1
                    local_count(i) = local_count(i) + 1
                ENDIF
            END DO
        END DO
        WRITE(*,*) count
        WRITE(*,*) MAXVAL(distance)
    END SUBROUTINE LoopDistance

END MODULE DistanceModule