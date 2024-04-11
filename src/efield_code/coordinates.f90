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

    SUBROUTINE Read_XYZ_Frame(xyz_unit, rO, r1, r2, rmol)
    ! ************************************************************************************
    ! This subroutine reads a frame from the trajectory file in the XYZ format
    ! In:
    !  -nmoltypes: number of molecule types
    !  -nmols: array of number of molecules of each type
    !  -natoms: array of number of atoms of each type
    !  -which_is_wat: index of the water molecule type
    !  -box_lengths: box length
    ! Out:
    !  -rO: array of oxygen positions
    !  -r1: array of hydrogen 1 positions
    !  -r2: array of hydrogen 2 positions
    !  -rmol: array of all the other atoms positions
    ! ************************************************************************************
        USE InputData
        IMPLICIT NONE

        INTEGER, INTENT(IN) :: xyz_unit
        REAL, DIMENSION(:,:), INTENT(INOUT) :: rO, r1, r2
        REAL, DIMENSION(:,:,:,:), INTENT(INOUT) :: rmol

        INTEGER :: i, type, jatom,k
        CHARACTER(LEN=10) :: ctmp


        DO i=1,2
            READ(xyz_unit,*) 
        END DO

        DO type=1, nmoltypes
            IF (type == which_is_water) THEN
                DO i=1, nmols(type)
                    ! Reading rule for water
                    read(xyz_unit,*) ctmp, (rO(i,k), k=1,3)
                    read(xyz_unit,*) ctmp, (r1(i,k), k=1,3)
                    read(xyz_unit,*) ctmp, (r2(i,k), k=1,3)
                    ! Make the molecule whole
                    DO k=1, 3
                        r1(i,k) = r1(i,k) - box_lengths(k)*anint((r1(i,k)-rO(i,k))/box_lengths(k))
                        r2(i,k) = r2(i,k) - box_lengths(k)*anint((r2(i,k)-rO(i,k))/box_lengths(k))
                    ENDDO ! k
                ENDDO ! i
            ELSE ! (type == which_is_wat)
                ! Reading rule for not water
                DO i=1, nmols(type)
                    DO jatom=1, natoms(type)
                        read(xyz_unit,*) ctmp, (rmol(type,i,jatom,k), k=1,3)
                        ! Make the molecule whole
                        IF (jatom > 1) THEN
                            DO k=1, 3
                                rmol(type,i,jatom,k) = rmol(type,i,jatom,k) &
                                & - box_lengths(k)*anint((rmol(type,i,jatom,k)-rmol(type,i,1,k))/box_lengths(k))
                            ENDDO ! k
                        END IF ! (jatom > 1)
                    END DO !jatom
                ENDDO ! i
            END IF ! type == which_is_wat
        END DO ! type
    END SUBROUTINE Read_XYZ_Frame


    SUBROUTINE Read_XYZ_Next(xyz_unit, positions)
    ! ************************************************************************************
    ! This subroutine reads a frame from the trajectory file in the XYZ format
    ! In:
    !  -nmoltypes: number of molecule types
    !  -nmols: array of number of molecules of each type
    !  -natoms: array of number of atoms of each type
    !  -which_is_wat: index of the water molecule type
    !  -box_lengths: box length
    ! Out:
    !  -positions: array of oxygen positions
    ! 
    ! ************************************************************************************
        USE InputData
        IMPLICIT NONE

        INTEGER, INTENT(IN) :: xyz_unit
        REAL, DIMENSION(:,:) :: positions
        
        INTEGER :: i, k
        INTEGER :: type, jatom
        INTEGER :: count, start_count
        CHARACTER(LEN=10) :: ctmp


        DO i=1,2
            READ(xyz_unit,*) 
        END DO

        count = 1
        DO type=1, nmoltypes
            IF (type == which_is_water) THEN
                wstart = count
            ENDIF
            DO i=1, nmols(type)
                DO jatom=1, natoms(type)
                    read(xyz_unit,*) ctmp, (positions(count,k), k=1,3)
                    ! Make the molecule whole
                    IF (jatom > 1) THEN
                        DO k=1, 3
                            positions(count,k) = positions(count,k) &
                            & - box_lengths(k)*anint((positions(count,k)-positions(start_count,k))/box_lengths(k))
                        ENDDO ! k
                    ELSE
                        start_count = count
                    END IF ! (jatom > 1)
                    count = count + 1
                END DO !jatom
            ENDDO ! i
            IF (type == which_is_water) THEN
                wend = count - 1
            ENDIF
        END DO ! type
    END SUBROUTINE Read_XYZ_Next

INCLUDE 'field_functions.f90'

END MODULE CoordinatesModule

