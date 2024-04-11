MODULE CoordinatesModule
    IMPLICIT NONE
    CONTAINS

    SUBROUTINE Coordinates(natoms, box, positions)
        IMPLICIT NONE
        INTEGER, INTENT(IN) :: natoms
        REAL, DIMENSION(3), INTENT(IN) :: box
        REAL, DIMENSION(:,:), INTENT(INOUT) :: positions
        INTEGER :: i, j

        ! Randomize Positions in Box
        DO i = 1, natoms
            DO j = 1, 3
                positions(i,j) = box(j) * RAND()
            END DO
        END DO
    END SUBROUTINE Coordinates

    SUBROUTINE WriteXYZ(natoms, positions, filename)
        IMPLICIT NONE
        INTEGER, INTENT(IN) :: natoms
        REAL, DIMENSION(:,:), INTENT(IN) :: positions
        CHARACTER(LEN=*), INTENT(IN) :: filename
        INTEGER :: i, j

        OPEN(UNIT=10, FILE=trim(filename), STATUS='UNKNOWN')
        WRITE(10,*) natoms
        WRITE(10,*) "xyz file"

        DO i = 1, natoms
            WRITE(10,*) i, positions(i,1), positions(i,2), positions(i,3)
        END DO

        CLOSE(10)
    END SUBROUTINE WriteXYZ

    SUBROUTINE ReadXYZ(natoms, positions, filename)
        IMPLICIT NONE
        INTEGER, INTENT(OUT) :: natoms
        REAL, DIMENSION(:,:), ALLOCATABLE, INTENT(OUT) :: positions
        CHARACTER(LEN=*), INTENT(IN) :: filename
        INTEGER :: i, j

        OPEN(UNIT=10, FILE=trim(filename), STATUS='OLD')
        READ(10,*) natoms
        READ(10,*)

        ALLOCATE(positions(natoms,3))

        DO i = 1, natoms
            READ(10,*) j, positions(i,1), positions(i,2), positions(i,3)
        END DO

        CLOSE(10)
    END SUBROUTINE ReadXYZ


END MODULE CoordinatesModule

