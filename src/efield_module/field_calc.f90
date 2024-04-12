! **************************************************************************************************
! This is the main subroutine for calculating the electric field at each atom in the system.
! The electric field is calculated by summing the contributions from all other atoms in the system.
! 
! The electric field is calculated by first building a cell list, and then calculating the electric
! field at each atom by summing the contributions from all other atoms in the system.
!
! Parameters:
!   positions:  The positions of all atoms in the system.
!   charges:    The charges of all atoms in the system.
!   molids:     The molecule IDs of all atoms in the system.
!   atypes:     The atom types of all atoms in the system.
!   electric_field1: The electric field at each atom in the system.
!   electric_field2: The electric field at each atom in the system.
!
! This uses the following subroutines:
!   Construct_Grid:  Builds the cell list. (From the CellList module)
!   Linked_List_Build:  Builds the linked list. (From the CellList module)
!   WaterListCalc:  Calculates the electric field for water molecules.
!   OtherListCalc:  Calculates the electric field for other molecules.
!
! Copyright (C) 2024 Zeke Piskulich
! **************************************************************************************************


SUBROUTINE FieldCalc(positions, charges, molids, atypes, electric_field1, electric_field2)
    USE InputData
    USE CellList
    IMPLICIT NONE
    ! Input Arguments
    REAL, DIMENSION(:,:),   INTENT(IN) :: positions
    REAL, DIMENSION(:),     INTENT(IN) :: charges
    INTEGER, DIMENSION(:),  INTENT(IN) :: molids, atypes
    
    ! Cell List Declarations
    INTEGER, DIMENSION(3)                    :: cell_nbins
    INTEGER, DIMENSION(:,:,:),  ALLOCATABLE  :: ll_head
    INTEGER, DIMENSION(:,:),    ALLOCATABLE  :: cell_map  
    INTEGER, DIMENSION(:),      ALLOCATABLE  :: cell_list 

    REAL, DIMENSION(:,:)                     :: electric_field1, electric_field2

    ! Other Variables
    INTEGER :: i, j, k, l, m, n
    INTEGER :: moltype

    REAL :: r_cutoff_sq



    electric_field1 = 0.0; electric_field2 = 0.0


    ! ****** Build the Cell List ********
    ALLOCATE(cell_list(SIZE(positions,1)))
    CALL Construct_Grid(box_lengths, r_cutoff, cell_nbins, cell_map)
    ALLOCATE(ll_head(cell_nbins(1), cell_nbins(2), cell_nbins(3)))
    ll_head = 0; cell_list = 0
    CALL Linked_List_Build(positions, box_lengths, r_cutoff, cell_nbins, ll_head, cell_list, cell_map)

    ! ****** Get Distances *******
    electric_field1 = 0.0; electric_field2 = 0.0
    DO moltype=1, nmoltypes
        IF (moltype==which_is_water) THEN
            ! First OH Calculation
            CALL WaterListCalc(positions, charges, molids, atypes, 1, &
                        & cell_list, cell_nbins, ll_head, cell_map, &
                        & electric_field1)
            ! Second OH Calculation
            CALL WaterListCalc(positions, charges, molids, atypes, 2, &
                        & cell_list, cell_nbins, ll_head, cell_map, &
                        & electric_field2)
        ELSE
            CALL OtherListCalc(positions, charges, molids, atypes, 1, &
                        & cell_list, cell_nbins, ll_head, cell_map, &
                        & electric_field1)
            CALL OtherListCalc(positions, charges, molids, atypes, 2, &
                        & cell_list, cell_nbins, ll_head, cell_map, &
                        & electric_field2)
        ENDIF
    ENDDO

    DEALLOCATE(cell_map)
    DEALLOCATE(cell_list); DEALLOCATE(ll_head)

END SUBROUTINE FieldCalc