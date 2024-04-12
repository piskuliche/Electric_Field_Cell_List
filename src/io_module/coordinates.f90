! ************************************************************************************
! CoordinatesModule - Module for handling coordinates
! This module contains subroutines for reading and writing coordinates in the XYZ format.
! 
! Subroutines:
!  - Coordinates: Randomizes the positions of the atoms in the box
!  - WriteXYZ: Writes the positions of the atoms in the XYZ format
!  - Read_XYZ_Frame: Reads a frame from the trajectory file in the XYZ format
!  - Read_XYZ_Next: Reads the next frame from the trajectory file in the XYZ format
!
! Copyright (C) 2024, Zeke Piskulich, Boston University
! ************************************************************************************

SUBROUTINE Coordinates(natoms, box, positions)
! ************************************************************************************
! This subroutine randomizes the positions of the atoms in the box. This is not used
! in the actual simulation, but is useful for generating initial configurations for 
! testing purposes.
!
! Parameters:
! -----------
! natoms: integer
!       Number of atoms in the system
! box: real, dimension(3)
!       Box dimensions in the x, y, and z directions
! positions: real, dimension(:,:)
!       Array containing the positions of the atoms. First dimension is the number of
!       atoms. The second dimension is the dimensionality of the system.(3)
!
! ************************************************************************************
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
! ************************************************************************************
! This subroutine writes the positions of the atoms in the XYZ format
!
! Parameters:
! -----------
! natoms: integer
!       Number of atoms in the system
! positions: real, dimension(:,:)
!       Array containing the positions of the atoms. First dimension is the number of
!       atoms. The second dimension is the dimensionality of the system.(3)
! filename: character(len=*)
!       Name of the file to write the positions to.
!
! ************************************************************************************
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

SUBROUTINE Read_XYZ_Next(xyz_unit, positions)
! ************************************************************************************
! This subroutine reads a frame from the trajectory file in the XYZ format
! 
! Parameters:
! -----------
! xyz_unit: integer
!       Unit number of the trajectory file, must be an integer. The file must be open 
!       before calling this subroutine.
! positions: real, dimension(:,:)
!       Array for the positions to be stored by the read operation. First dimension is 
!       the number of atoms. The second dimension is dimensionality of the system.(3)
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
                ! This next segment applies PBC to make molecules remain whole.
                IF (jatom > 1) THEN
                    DO k=1, 3
                        positions(count,k) = positions(count,k) &
                        & - box_lengths(k)*anint((positions(count,k)-positions(start_count,k))/box_lengths(k))
                    ENDDO 
                ELSE
                    start_count = count
                END IF 
                count = count + 1
            END DO 
        ENDDO 
        ! Stores the ending location of the water molecules.
        IF (type == which_is_water) THEN
            wend = count - 1
        ENDIF
    END DO 
END SUBROUTINE Read_XYZ_Next


