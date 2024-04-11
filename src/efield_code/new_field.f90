PROGRAM ElectricField
    USE CellList
    USE InputData
    USE AllocHandler
    USE CoordinatesModule
    USE EfieldCalc
    USE InputModule
    USE Units
    USE HDF5_IO
    IMPLICIT NONE

    INTEGER :: i, k
    INTEGER :: chunk, number_of_chunks
    INTEGER :: z, iconfig

    REAL, ALLOCATABLE, DIMENSION(:,:,:) :: positions


    REAL, ALLOCATABLE, DIMENSION(:,:)    :: dot1, dot2
    REAL, ALLOCATABLE, DIMENSION(:,:,:)  :: eOH1, eOH2
    REAL, ALLOCATABLE, DIMENSION(:,:)    :: z0
    REAL,              DIMENSION(3)      :: ef1_tmp, ef2_tmp

    REAL, ALLOCATABLE, DIMENSION(:)      :: charges
    INTEGER, ALLOCATABLE, DIMENSION(:)   :: molids, atypes

    INTEGER, ALLOCATABLE, DIMENSION(:) :: jcontributes
    REAL, DIMENSION(3) :: ef_component
    REAL, ALLOCATABLE, DIMENSION(:) :: ef_components
    REAL, ALLOCATABLE, DIMENSION(:,:,:) :: electric_field1, electric_field2

    INTEGER, ALLOCATABLE, DIMENSION(:)   :: oo_count
    INTEGER :: molecule_type, j

    REAL :: dist
    REAL, DIMENSION(3) :: dist_vec

    REAL :: tstart, tend

    ! File Opens
    OPEN(12, FILE='traj.xyz', STATUS='OLD', ACTION='READ')

    ! Read Input 
    CALL Read_Input
    CALL DisplayInput


    ! Allocates the Position Arrays
    CALL AllocateArrays(positions, dot1, dot2, eOH1, eOH2, z0)
    ALLOCATE(oo_count(nmols(which_is_water)))
    number_of_atoms = SIZE(positions,1)
    ALLOCATE(charges(number_of_atoms))
    ALLOCATE(molids(number_of_atoms))
    ALLOCATE(atypes(number_of_atoms))
    ALLOCATE(electric_field1(number_of_atoms,3,nconfig))
    ALLOCATE(electric_field2(number_of_atoms,3,nconfig))

    CALL BuildArrays(charges, molids, atypes)

    WRITE(*,*) "ALL ARRAYS ALLOCATED"

    number_of_chunks = CEILING(REAL(nconfig)/REAL(max_config))

    CALL CPU_TIME(tstart)

    DO chunk=1, number_of_chunks
        WRITE(*,*) "Starting chunk ", chunk, " of ", number_of_chunks
        ! Read the Frame
        positions = 0.0
        DO z=1, max_config
            IF ((chunk-1)*max_config + z > nconfig) then
                EXIT
            ENDIF
            CALL Read_XYZ_Next(12, positions(:,:,z))
        ENDDO ! z 

        ! Calculate OH bond vectors
        DO z=1, MIN(nconfig, max_config)
            iconfig = (chunk-1)*max_config + z
            
            CALL BondVector(positions(wstart+1:wend:3,:,z), positions(wstart:wend:3,:,z), eOH1(:,:,iconfig))
            CALL BondVector(positions(wstart+2:wend:3,:,z), positions(wstart:wend:3,:,z), eOH2(:,:,iconfig))
            
            electric_field1(:,:, iconfig) = 0.0; electric_field2(:,:,iconfig) = 0.0
            CALL FieldCalc(positions(:,:,z), charges, molids, atypes, electric_field1(:,:,iconfig), electric_field2(:,:,iconfig))
        ENDDO

        DO z=1, MIN(nconfig, max_config)
            iconfig = (chunk-1)*max_config + z
            DO i=1, nmols(which_is_water)
                electric_field1(i,:,iconfig) = electric_field1(i,:,iconfig) * angperau**2.
                electric_field2(i,:,iconfig) = electric_field2(i,:,iconfig) * angperau**2.
                dot1(i, iconfig) = DOT_PRODUCT(eOH1(i,:,iconfig), electric_field1(i,:,iconfig))
                dot2(i, iconfig) = DOT_PRODUCT(eOH2(i,:,iconfig), electric_field2(i,:,iconfig))
                z0(i, iconfig) = positions(wstart + (i-1)*3,3, z)
            ENDDO
        ENDDO

    ENDDO 

    CALL CPU_TIME(tend)

    WRITE(*,*) "Main Loop Execution Time", tend-tstart

    CALL CPU_TIME(tstart)
    WRITE(*,*) shape(dot1)
    CALL WRITE_HD5F(dot1, dot2, eOH1, eOH2, z0, nmols(which_is_water), nconfig)
    CALL CPU_TIME(tend)
    WRITE(*,*) "File Write Time", tend-tstart

    ! Write HDF5 Archive
    CALL DeallocateArrays(positions,dot1,dot2,eOH1,eOH2,z0)
    DEALLOCATE(oo_count)
    ! File Closing
    CLOSE(12)

END PROGRAM ElectricField

