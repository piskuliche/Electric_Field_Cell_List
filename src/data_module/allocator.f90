! **************************************************************************************************
! This module contains the subroutines that allocate and deallocate the arrays used in the program.
! 
! Copyright (C) 2024 Zeke Piskulich
! **************************************************************************************************

MODULE AllocHandler
    IMPLICIT NONE
    CONTAINS

    SUBROUTINE AllocateArrays(positions, dot1, dot2, eOH1, eOH2, z0, &
                            & oo_count, charges, molids, atypes, &
                            & electric_field1, electric_field2)
    ! **************************************************************************************************
    ! This subroutine allocates the arrays used in the program.
    !
    ! Parameters:
    ! -----------
    ! positions : REAL, ALLOCATABLE, DIMENSION(:,:,:)
    !     The positions of the atoms in the system.
    ! dot1 : REAL, ALLOCATABLE, DIMENSION(:,:)
    !     The dot product of the electric field and the vector pointing from the oxygen to the hydrogen.
    ! dot2 : REAL, ALLOCATABLE, DIMENSION(:,:)
    !     The dot product of the electric field and the vector pointing from the oxygen to the hydrogen.
    ! eOH1 : REAL, ALLOCATABLE, DIMENSION(:,:,:)
    !     The first OH bond vector.
    ! eOH2 : REAL, ALLOCATABLE, DIMENSION(:,:,:)
    !     The second OH bond vector
    ! z0 : REAL, ALLOCATABLE, DIMENSION(:,:)
    !     The z-coordinate of the oxygen atom.
    ! oo_count : INTEGER, ALLOCATABLE, DIMENSION(:)
    !     The number of water molecules in each configuration.
    ! charges : REAL, ALLOCATABLE, DIMENSION(:)
    !     The charges of the atoms in the system.
    ! molids : INTEGER, ALLOCATABLE, DIMENSION(:)
    !     The molecule id of each atom.
    ! atypes : INTEGER, ALLOCATABLE, DIMENSION(:)
    !     The atom type of each atom.
    ! electric_field1 : REAL, ALLOCATABLE, DIMENSION(:,:,:)
    !     The electric field at each atom.
    ! electric_field2 : REAL, ALLOCATABLE, DIMENSION(:,:,:)
    !     The electric field at each atom.
    !
    ! **************************************************************************************************
        USE InputData
        IMPLICIT NONE
        INTEGER, ALLOCATABLE,   DIMENSION(:)      :: oo_count, molids, atypes
        REAL, ALLOCATABLE,      DIMENSION(:)      :: charges
        REAL, ALLOCATABLE,      DIMENSION(:,:)    :: dot1, dot2, z0
        REAL, ALLOCATABLE,      DIMENSION(:,:,:)  :: eOH1, eOH2, positions, electric_field1, electric_field2

        INTEGER :: i
        INTEGER :: atom_count
        
        atom_count = 0
        DO i=1, nmoltypes
            atom_count = atom_count + nmols(i)*natoms(i)
        ENDDO
        WRITE(*,*) "There are ", atom_count, " atoms."

        ALLOCATE(positions(atom_count, ndim, max_config))

        ! Output Allocations
        ALLOCATE(dot1(nmols(which_is_water),    nconfig))
        ALLOCATE(dot2(nmols(which_is_water),    nconfig))
        ALLOCATE(z0(nmols(which_is_water),      nconfig))
        ALLOCATE(eOH1(nmols(which_is_water),    ndim,   nconfig))
        ALLOCATE(eOH2(nmols(which_is_water),    ndim,   nconfig))
        ALLOCATE(oo_count(nmols(which_is_water)))
        ALLOCATE(charges(atom_count))
        ALLOCATE(molids(atom_count))
        ALLOCATE(atypes(atom_count))
        ALLOCATE(electric_field1(atom_count,3,nconfig))
        ALLOCATE(electric_field2(atom_count,3,nconfig))


        WRITE(*,*) "** Arrays have been allocated **"
        WRITE(*,*) "positions shape: ", SHAPE(positions)
        WRITE(*,*) "dot1 shape: ", SHAPE(dot1)
        WRITE(*,*) "dot2 shape: ", SHAPE(dot2)
        WRITE(*,*) "z0 shape: ", SHAPE(z0)
        WRITE(*,*) "eOH1 shape: ", SHAPE(eOH1)
        WRITE(*,*) "eOH2 shape: ", SHAPE(eOH2)
    END SUBROUTINE AllocateArrays
    
    SUBROUTINE DeallocateArrays(positions, dot1, dot2, eOH1, eOH2, z0, &
                            & oo_count, charges, molids, atypes, &
                            & electric_field1, electric_field2)
    ! **************************************************************************************************
    ! This subroutine deallocates the arrays used in the program.
    !
    ! Parameters:
    ! -----------
    ! positions : REAL, ALLOCATABLE, DIMENSION(:,:,:)
    !     The positions of the atoms in the system.
    ! dot1 : REAL, ALLOCATABLE, DIMENSION(:,:)
    !     The dot product of the electric field and the vector pointing from the oxygen to the hydrogen.
    ! dot2 : REAL, ALLOCATABLE, DIMENSION(:,:)
    !     The dot product of the electric field and the vector pointing from the oxygen to the hydrogen.
    ! eOH1 : REAL, ALLOCATABLE, DIMENSION(:,:,:)
    !     The first OH bond vector.
    ! eOH2 : REAL, ALLOCATABLE, DIMENSION(:,:,:)
    !     The second OH bond vector
    ! z0 : REAL, ALLOCATABLE, DIMENSION(:,:)
    !     The z-coordinate of the oxygen atom.
    ! oo_count : INTEGER, ALLOCATABLE, DIMENSION(:)
    !     The number of water molecules in each configuration.
    ! charges : REAL, ALLOCATABLE, DIMENSION(:)
    !     The charges of the atoms in the system.
    ! molids : INTEGER, ALLOCATABLE, DIMENSION(:)
    !     The molecule id of each atom.
    ! atypes : INTEGER, ALLOCATABLE, DIMENSION(:)
    !     The atom type of each atom.
    ! electric_field1 : REAL, ALLOCATABLE, DIMENSION(:,:,:)
    !     The electric field at each atom.
    ! electric_field2 : REAL, ALLOCATABLE, DIMENSION(:,:,:)
    !     The electric field at each atom.
    !
    ! **************************************************************************************************
        IMPLICIT NONE
        INTEGER, ALLOCATABLE,   DIMENSION(:)      :: oo_count, molids, atypes
        REAL, ALLOCATABLE,      DIMENSION(:)      :: charges
        REAL, ALLOCATABLE,      DIMENSION(:,:)    :: dot1, dot2, z0
        REAL, ALLOCATABLE,      DIMENSION(:,:,:)  :: eOH1, eOH2, positions, electric_field1, electric_field2

        IF (ALLOCATED(positions)) THEN
            DEALLOCATE(positions)
        END IF
        IF (ALLOCATED(dot1)) THEN
            DEALLOCATE(dot1)
        END IF
        IF (ALLOCATED(dot2)) THEN
            DEALLOCATE(dot2)
        END IF
        IF (ALLOCATED(z0)) THEN
            DEALLOCATE(z0)
        END IF
        IF (ALLOCATED(eOH1)) THEN
            DEALLOCATE(eOH1)
        END IF
        IF (ALLOCATED(eOH2)) THEN
            DEALLOCATE(eOH2)
        END IF
        IF (ALLOCATED(oo_count)) THEN
            DEALLOCATE(oo_count)
        END IF
        IF (ALLOCATED(charges)) THEN
            DEALLOCATE(charges)
        END IF
        IF (ALLOCATED(molids)) THEN
            DEALLOCATE(molids)
        END IF
        IF (ALLOCATED(atypes)) THEN
            DEALLOCATE(atypes)
        END IF
        IF (ALLOCATED(electric_field1)) THEN
            DEALLOCATE(electric_field1)
        END IF
        IF (ALLOCATED(electric_field2)) THEN
            DEALLOCATE(electric_field2)
        END IF
    END SUBROUTINE DeallocateArrays

END MODULE AllocHandler