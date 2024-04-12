PROGRAM ElectricField
    USE CellList
    USE InputData
    USE AllocHandler
    USE Efield_Module
    USE IO_Module
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

    ! Open the trajectory file to unit 12.
    OPEN(12, FILE='traj.xyz', STATUS='OLD', ACTION='READ')

    ! Read Input 
    CALL Read_Input
    CALL DisplayInput


    ! Allocates the allocatable arrays. All of these should be deallocated
    ! at the end of the program.
    CALL AllocateArrays(positions, dot1, dot2, eOH1, eOH2, z0, &
                        & oo_count, charges, molids, atypes, &
                        & electric_field1, electric_field2)
    number_of_atoms = SIZE(positions,1)

    ! Build the charges, molids, and atypes arrays based on the 
    ! input data. These arrays are used to distinguish between molecules
    ! and atom types. The charges array is used to calculate the electric
    ! field at each atom.
    CALL BuildArrays(charges, molids, atypes)

    WRITE(*,*) "ALL ARRAYS ALLOCATED"

    ! Calculate the number of chunks that will be read from the trajectory file.
    ! This is done so that the program can be run on a machine with limited memory.
    number_of_chunks = CEILING(REAL(nconfig)/REAL(max_config))

    CALL CPU_TIME(tstart)

    DO chunk=1, number_of_chunks
        WRITE(*,*) "Starting chunk ", chunk, " of ", number_of_chunks
        ! I. Read the next max_config configurations from the trajectory file.
        positions = 0.0
        DO z=1, max_config
            IF ((chunk-1)*max_config + z > nconfig) then
                EXIT
            ENDIF
            CALL Read_XYZ_Next(12, positions(:,:,z))
        ENDDO ! z 

        ! II. Calculate the OH bond vectors for the water molecules and the electric field.
        DO z=1, MIN(nconfig, max_config)
            iconfig = (chunk-1)*max_config + z
            
            ! II.A - Calculate the OH bond vectors for the water molecules.
            CALL BondVector(positions(wstart+1:wend:3,:,z), positions(wstart:wend:3,:,z), eOH1(:,:,iconfig))
            CALL BondVector(positions(wstart+2:wend:3,:,z), positions(wstart:wend:3,:,z), eOH2(:,:,iconfig))
            
            ! II.B - Calculate the electric field at each atom.
            electric_field1(:,:, iconfig) = 0.0; electric_field2(:,:,iconfig) = 0.0
            CALL FieldCalc(positions(:,:,z), charges, molids, atypes, electric_field1(:,:,iconfig), electric_field2(:,:,iconfig))
        ENDDO

        ! III. Calculate the dot product of the OH bond vectors and the electric field.
        DO z=1, MIN(nconfig, max_config)
            iconfig = (chunk-1)*max_config + z
            DO i=1, nmols(which_is_water)
                ! III.A - Convert the electric field to atomic units.
                electric_field1(i,:,iconfig) = electric_field1(i,:,iconfig) * angperau**2.
                electric_field2(i,:,iconfig) = electric_field2(i,:,iconfig) * angperau**2.

                ! III.B - Calculate the dot product of the OH bond vectors and the electric field.
                dot1(i, iconfig) = DOT_PRODUCT(eOH1(i,:,iconfig), electric_field1(i,:,iconfig))
                dot2(i, iconfig) = DOT_PRODUCT(eOH2(i,:,iconfig), electric_field2(i,:,iconfig))

                ! III.C - Store the z-coordinate of the water molecule.
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

    ! Deallocate the allocatable arrays that were allocated at the top of the program.
    CALL DeallocateArrays(positions,dot1,dot2,eOH1,eOH2,z0 &
                          & ,oo_count,charges,molids,atypes &
                          & ,electric_field1,electric_field2)

    ! Close the Trajectory File
    CLOSE(12)

END PROGRAM ElectricField

