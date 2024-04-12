! ****************************************************************************************************
! This subroutine calculates the electric field at each grid point due to the charges in the system. 
! This part works for the charges that do not belong to the water molecules.
! 
! Parameters:
!   positions: 2D array containing the positions of all the atoms in the system
!   charges: 1D array containing the charges of all the atoms in the system
!   molids: 1D array containing the molecule ids of all the atoms in the system
!   atypes: 1D array containing the atom types of all the atoms in the system
!   hi: integer containing the index of the hydrogen atom of interest in the system
!   cell_list: 1D array containing the linked list of atoms in each cell
!   cell_nbins: 1D array containing the number of cells in each direction
!   ll_head: 3D array containing the head of the linked list for each cell
!   cell_map: 2D array containing the mapping of the cells to the grid
!   electric_field: 2D array containing the electric field at each grid point
!
! Copyright 2024, Zeke Piskulich
! ****************************************************************************************************


SUBROUTINE OtherListCalc(positions, charges, molids, atypes, hi, &
                    & cell_list, cell_nbins, ll_head, cell_map, &
                    electric_field)
    USE InputData
    USE CellList
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