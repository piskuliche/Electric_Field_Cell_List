
SUBROUTINE Read_Input
! ************************************************************************************
! This subroutine reads the input file for the field calculation.
! 
! It also reads the molecular input files for each molecule type.
!
!
! ************************************************************************************
    USE InputData
    IMPLICIT NONE

    INTEGER :: i, nwater

    nmols = 0; natoms = 0
    chrgs = 0.0; box_lengths = 0.0
    
    ! Read the field input file
    OPEN(10, file='field_input.in', status='old')
        READ(10,*)
        READ(10,*) nmoltypes, nconfig, which_is_water
        READ(10,*) 
        READ(10,*) r_cutoff, box_lengths(1), box_lengths(2), box_lengths(3)
        READ(10,*)
        IF (nmoltypes > 10) STOP 'Too many molecule types'
        READ(10,*) (molnames(i), i=1, nmoltypes) ! molecule names
        READ(10,*)
        READ(10,*) (nmols(i), i=1, nmoltypes) ! number of mols
    CLOSE(10)

    DO i=1, nmoltypes
        CALL Read_Molecule(i, molnames(i), chrgs(i,:), atps(i,:), natoms(i))
        IF ( i == which_is_water) THEN
            nwater = nmols(i)
        END IF
    ENDDO

END SUBROUTINE Read_Input

SUBROUTINE BuildArrays(charges, molids, atypes)
    USE InputData
    IMPLICIT NONE
    INTEGER,    DIMENSION(:) :: molids, atypes
    REAL,       DIMENSION(:) :: charges

    INTEGER :: atom_count, mol_count
    INTEGER :: i, j, k

    INTEGER :: max_types

    charges = 0.0; molids=0; atypes=0
    atom_count = 1; mol_count = 1

    max_types = 0
    DO i=1, nmoltypes
    DO j=1, nmols(i)
        DO k=1, natoms(i)
            charges(atom_count) =   chrgs(i,k)
            molids(atom_count)  =   mol_count
            atypes(atom_count)  =   atps(i,k) + max_types
            atom_count = atom_count + 1
        ENDDO
    mol_count = mol_count + 1
    ENDDO
    max_types = MAXVAL(atypes)
    ENDDO
END SUBROUTINE 

SUBROUTINE DisplayInput
    USE InputData
    IMPLICIT NONE
    INTEGER :: i
    WRITE(*,*) '***********************************'
    WRITE(*,*) 'nmoltypes, nconfig, which_is_water'
    WRITE(*,*) nmoltypes, nconfig, which_is_water
    WRITE(*,*) 'r_cutoff, box_lengths(:)'
    WRITE(*,*) r_cutoff, box_lengths(1), box_lengths(2), box_lengths(3)
    WRITE(*,*) 'molnames(:)'
    IF (nmoltypes > 10) STOP 'Too many molecule types'
    WRITE(*,*) (molnames(i), i=1, nmoltypes) ! molecule names
    WRITE(*,*) 'numberofmols(:)'
    WRITE(*,*) (nmols(i), i=1, nmoltypes) ! number of mols

    WRITE(*,*) 'Charges for moltype 1'
    WRITE(*,*) (chrgs(1,i), i=1,3)
    WRITE(*,*) '***********************************'

END SUBROUTINE DisplayInput

SUBROUTINE Read_Molecule(imol, molname, q, atype, natoms)
! ************************************************************************************
! This subroutine reads the input file for a molecule
!
! Input:
!  -imol: index of the molecule
!  -molname: name of the molecule
!
! Output:
!  -q: array of chrgs of each atom
!  -natoms: number of atoms
! ************************************************************************************
    IMPLICIT NONE

    CHARACTER(LEN=10), INTENT(IN) :: molname
    INTEGER, INTENT(IN) :: imol
    INTEGER, INTENT(OUT) :: natoms
    REAL, DIMENSION(2000) :: q
    INTEGER, DIMENSION(2000) :: atype

    INTEGER :: i, ifile

    ifile = imol + 50

    atype = 0

    OPEN(ifile, file=trim(molname)//".in", status='old')
        READ(ifile,*) natoms
        DO i=1, natoms
            READ(ifile,*) q(i), atype(i)
        END DO
    CLOSE(ifile)

END SUBROUTINE Read_Molecule





