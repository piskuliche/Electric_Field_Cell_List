PROGRAM ElectricField
    USE CoordinatesModule
    USE DistanceModule
    USE CellList
    IMPLICIT NONE

    INTEGER :: number_of_atoms 
    REAL, DIMENSION(3) :: box
    REAL, ALLOCATABLE, DIMENSION(:,:) :: positions, distance
    INTEGER, ALLOCATABLE, DIMENSION(:) :: local_loop, local_cell

    INTEGER :: i

    REAL :: tstart, tend
    REAL :: r_cutoff

    r_cutoff = 2.0

    !number_of_atoms = 10000
    
    CALL ReadXYZ(number_of_atoms, positions, 'coordinates.xyz')

    !ALLOCATE(positions(number_of_atoms,3))
    ALLOCATE(distance(number_of_atoms,number_of_atoms))
    
    WRITE(*,*) 'Electric Field Calculation'
    box = 30.0

    positions = 0.0


    CALL Coordinates(number_of_atoms, box, positions)

    !CALL WriteXYZ(number_of_atoms, positions, 'coordinates.xyz')


    ALLOCATE(local_loop(number_of_atoms))
    ALLOCATE(local_cell(number_of_atoms))

    CALL cpu_time(tstart)
    CALL LoopDistance(positions, r_cutoff, box, distance, local_loop)
    CALL cpu_time(tend)
    WRITE(*,*) "Iteration Time Loop ", tend - tstart, " seconds"

    CALL cpu_time(tstart)
    CALL Cell_Distance(positions, positions, box, r_cutoff, local_cell)
    CALL cpu_time(tend)
    WRITE(*,*) "Iteration Time Cell ", tend - tstart, " seconds"

    DO i=1, number_of_atoms
        !IF (local_loop(i) - local_cell(i) .ne. 0) THEN
        !    WRITE(*,*) "Calculation error", i
        !    WRITE(*,*) local_loop(i)
        !    WRITE(*,*) local_cell(i)
        !    WRITE(*,*) positions(i,:)
        !ENDIF 
    ENDDO

    

    DEALLOCATE(positions)

END PROGRAM ElectricField