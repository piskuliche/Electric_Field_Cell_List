MODULE AllocHandler
    IMPLICIT NONE
    CONTAINS

    SUBROUTINE AllocateArrays(positions, dot1, dot2, eOH1, eOH2, z0)
        USE InputData
        IMPLICIT NONE
        REAL, ALLOCATABLE, DIMENSION(:,:,:) :: positions
            
        REAL, ALLOCATABLE, DIMENSION(:,:)    :: dot1, dot2
        REAL, ALLOCATABLE, DIMENSION(:,:,:)  :: eOH1, eOH2
        REAL, ALLOCATABLE, DIMENSION(:,:)    :: z0

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


        WRITE(*,*) "** Arrays have been allocated **"
        WRITE(*,*) "positions shape: ", SHAPE(positions)
        WRITE(*,*) "dot1 shape: ", SHAPE(dot1)
        WRITE(*,*) "dot2 shape: ", SHAPE(dot2)
        WRITE(*,*) "z0 shape: ", SHAPE(z0)
        WRITE(*,*) "eOH1 shape: ", SHAPE(eOH1)
        WRITE(*,*) "eOH2 shape: ", SHAPE(eOH2)
    END SUBROUTINE AllocateArrays
    
    SUBROUTINE DeallocateArrays(positions, dot1, dot2, eOH1, eOH2, z0)
        IMPLICIT NONE
        REAL, DIMENSION(:,:,:), ALLOCATABLE :: positions

        REAL, ALLOCATABLE, DIMENSION(:,:)    :: dot1, dot2
        REAL, ALLOCATABLE, DIMENSION(:,:,:)  :: eOH1, eOH2
        REAL, ALLOCATABLE, DIMENSION(:,:)    :: z0

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
    END SUBROUTINE DeallocateArrays

END MODULE AllocHandler