
MODULE EfieldCalc
    IMPLICIT NONE
    CONTAINS

    SUBROUTINE FieldCalc(positions, charges, molids, atypes, electric_field1, electric_field2)
        USE InputData
        USE CellList
        USE Units
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

    END SUBROUTINE


    SUBROUTINE WaterListCalc(positions, charges, molids, atypes, hi, &
                        & cell_list, cell_nbins, ll_head, cell_map, &
                        electric_field)
        USE InputData
        USE CellList
        uSE CoordinatesModule
        USE Units
        USE, INTRINSIC :: IEEE_ARITHMETIC
        IMPLICIT NONE
        REAL,    DIMENSION(:,:),    INTENT(IN) :: positions
        REAL,    DIMENSION(:),      INTENT(IN) :: charges
        INTEGER, DIMENSION(:),      INTENT(IN) :: molids, atypes, cell_list, cell_nbins
        INTEGER, DIMENSION(:,:),    INTENT(IN) :: cell_map
        INTEGER, DIMENSION(:,:,:),  INTENT(IN) :: ll_head
        INTEGER,                    INTENT(IN) :: hi

        REAL, DIMENSION(:,:),       INTENT(INOUT) :: electric_field

        INTEGER :: hid, l, m, n, ii, jj, kk, jhead
        INTEGER :: count
        INTEGER, DIMENSION(3) :: hid_grid
        REAL, DIMENSION(4) :: tmp_dr
        REAL :: r_cutoff_sq

        INTEGER :: wcount
        
        r_cutoff_sq = r_cutoff*r_cutoff 
        count = 1
        DO hid=wstart+hi, wend, 3
            wcount = 0
            hid_grid = Assign_To_Cell_Grid(positions(hid,:), cell_nbins, box_lengths)
            DO l = -2, 2
            DO m = -2, 2
            DO n = -2, 2
                ii = cell_map(1, hid_grid(1) + l + 2)
                jj = cell_map(2, hid_grid(2) + m + 2)
                kk = cell_map(3, hid_grid(3) + n + 2)
                
                jhead = ll_head(ii, jj, kk)
                DO WHILE (jhead /= 0)
                    IF ((atypes(jhead) == atypes(wstart)) .and. (molids(hid) /= molids(jhead))) THEN ! ONLY go for oxygen atoms
                        tmp_dr = PBC_dr_new(positions(jhead,:), positions(hid,:), box_lengths)
                        IF (tmp_dr(4) < r_cutoff_sq) THEN
                            wcount = wcount + 1
                            ! Add contribution to the electric field from o_atom
                            electric_field(count,:) = electric_field(count,:) + &
                                    &  Field_Contribution(charges(jhead), positions(hid,:), &
                                        & tmp_dr(1:3), SQRT(tmp_dr(4)))

                            ! Add contribution to the electric field from h_atom 1
                            tmp_dr = PBC_dr_new(positions(jhead+1,:), positions(hid,:), box_lengths)
                            electric_field(count,:) = electric_field(count,:) &
                                    & + Field_Contribution(charges(jhead+1), positions(hid,:), tmp_dr(1:3), SQRT(tmp_dr(4)))

                            ! Add contribution to the electric field from h_atom 2
                            tmp_dr = PBC_dr_new(positions(jhead+2,:), positions(hid,:), box_lengths)
                            electric_field(count,:) = electric_field(count,:) &
                                    & + Field_Contribution(charges(jhead+2), positions(hid,:), tmp_dr(1:3), SQRT(tmp_dr(4)))
                        ENDIF
                    ENDIF
                    jhead = cell_list(jhead)
                ENDDO
            ENDDO
            ENDDO
            ENDDO
            count = count + 1
        ENDDO

    END SUBROUTINE

    SUBROUTINE OtherListCalc(positions, charges, molids, atypes, hi, &
                        & cell_list, cell_nbins, ll_head, cell_map, &
                        electric_field)
        USE InputData
        USE CellList
        uSE CoordinatesModule
        USE, INTRINSIC :: IEEE_ARITHMETIC
        IMPLICIT NONE
        REAL,    DIMENSION(:,:),    INTENT(IN) :: positions
        REAL,    DIMENSION(:),      INTENT(IN) :: charges
        INTEGER, DIMENSION(:),      INTENT(IN) :: molids, atypes, cell_list, cell_nbins
        INTEGER, DIMENSION(:,:),    INTENT(IN) :: cell_map
        INTEGER, DIMENSION(:,:,:),  INTENT(IN) :: ll_head
        INTEGER,                    INTENT(IN) :: hi

        REAL, DIMENSION(:,:),       INTENT(INOUT) :: electric_field

        INTEGER :: hid, l, m, n, ii, jj, kk, jhead
        INTEGER :: count, wcount
        INTEGER, DIMENSION(3) :: hid_grid
        REAL, DIMENSION(4) :: tmp_dr
        REAL :: r_cutoff_sq
        r_cutoff_sq = r_cutoff*r_cutoff 
        count = 1
        DO hid=wstart+hi, wend, 3
            wcount = 0
            hid_grid = Assign_To_Cell_Grid(positions(hid,:), cell_nbins, box_lengths)
            DO l = -2, 2
            DO m = -2, 2
            DO n = -2, 2
                ii = cell_map(1, hid_grid(1) + l + 2)
                jj = cell_map(2, hid_grid(2) + m + 2)
                kk = cell_map(3, hid_grid(3) + n + 2)
                
                jhead = ll_head(ii, jj, kk)
                DO WHILE (jhead /= 0)
                    IF ((molids(hid) /= molids(jhead)) .and. (jhead > wend .or. jhead < wstart)) THEN ! ONLY go for oxygen atoms
                        tmp_dr = PBC_dr_new(positions(jhead,:), positions(hid,:), box_lengths)
                        IF (tmp_dr(4) < r_cutoff_sq) THEN
                            wcount = wcount + 1
                            ! Add contribution to the electric field
                            electric_field(count,:) = electric_field(count,:) + &
                                    &  Field_Contribution(charges(jhead), positions(hid,:), &
                                        & tmp_dr(1:3), SQRT(tmp_dr(4)))
                        ENDIF
                    ENDIF
                    jhead = cell_list(jhead)
                ENDDO
            ENDDO
            ENDDO
            ENDDO
            count = count + 1
        ENDDO

    END SUBROUTINE OtherListCalc
END MODULE EfieldCalc