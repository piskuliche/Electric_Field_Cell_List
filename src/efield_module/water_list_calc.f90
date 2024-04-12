! ****************************************************************************************************
! WaterListCalc: This subroutine calculates the electric field at the position of the hydrogen atoms
!               of the water molecules. The electric field is calculated by summing the contribution
!               from the oxygen atom and the two hydrogen atoms of the water molecule. The electric
!               field is calculated using the formula:
!               E = q * (r - r0) / |r - r0|^3
!               where q is the charge of the atom, r is the position of the atom, r0 is the position
!               of the hydrogen atom and |r - r0| is the distance between the two atoms.
!               The electric field is calculated for all the water molecules in the system.
! Parameters:
!       positions: 2D array containing the positions of all the atoms in the system
!       charges: 1D array containing the charges of all the atoms in the system
!       molids: 1D array containing the molecule id of all the atoms in the system
!       atypes: 1D array containing the atom type of all the atoms in the system
!       hi: Integer containing the index of the hydrogen atom in the system
!       cell_list: 1D array containing the linked list of atoms in each cell
!       cell_nbins: 1D array containing the number of bins in each dimension of the cell grid
!       ll_head: 3D array containing the head of the linked list of atoms in each cell
!       cell_map: 2D array containing the mapping of the cell grid to the cell list
!       electric_field: 2D array containing the electric field at the position of the hydrogen atoms
!
! Copyright 2024 (C) Zeke Piskulich, Boston University
! ****************************************************************************************************

SUBROUTINE WaterListCalc(positions, charges, molids, atypes, hi, &
                    & cell_list, cell_nbins, ll_head, cell_map, &
                    electric_field)
    USE InputData
    USE CellList
    uSE IO_Module
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
    ! Setup the square of the cutoff distance
    r_cutoff_sq = r_cutoff*r_cutoff

    count = 1
    DO hid=wstart+hi, wend, 3
        wcount = 0
        ! I. Assign the position of hid to a cell on the grid.
        hid_grid = Assign_To_Cell_Grid(positions(hid,:), cell_nbins, box_lengths)
        ! II. Loop over the neighboring cells
        DO l = -2, 2
        DO m = -2, 2
        DO n = -2, 2
            ii = cell_map(1, hid_grid(1) + l + 2)
            jj = cell_map(2, hid_grid(2) + m + 2)
            kk = cell_map(3, hid_grid(3) + n + 2)
            
            ! III. Loop over the atoms in the neighboring cells
            jhead = ll_head(ii, jj, kk)
            DO WHILE (jhead /= 0)
                ! The section of code in this if statment calculates the electric field on hid.
                ! This works by first checking that the atom is an oxygen atom and that the atom is not
                ! in the same molecule as hid. If these conditions are met, the electric field is calculated
                ! by summing the contribution from the oxygen atom and the two hydrogen atoms of the water
                ! molecule.
                IF ((atypes(jhead) == atypes(wstart)) .and. (molids(hid) /= molids(jhead))) THEN 
                    tmp_dr = PBC_dr_new(positions(jhead,:), positions(hid,:), box_lengths)                                  ! Get the OO Distance
                    IF (tmp_dr(4) < r_cutoff_sq) THEN                                                                       ! If in Cutoff Distance, Calculate rest of values
                        wcount = wcount + 1
                        ! Add contribution to the electric field from o_atom
                        electric_field(count,:) = electric_field(count,:) + &                                               ! Get O contribution to the electric field
                                &  Field_Contribution(charges(jhead), positions(hid,:), &
                                    & tmp_dr(1:3), SQRT(tmp_dr(4)))

                        ! Add contribution to the electric field from h_atom 1
                        tmp_dr = PBC_dr_new(positions(jhead+1,:), positions(hid,:), box_lengths)                            ! Calculate the HH distance
                        electric_field(count,:) = electric_field(count,:) &                                                 ! Calculate the HH contribution to the electric field
                                & + Field_Contribution(charges(jhead+1), positions(hid,:), tmp_dr(1:3), SQRT(tmp_dr(4)))

                        ! Add contribution to the electric field from h_atom 2
                        tmp_dr = PBC_dr_new(positions(jhead+2,:), positions(hid,:), box_lengths)                            ! Calculate the HH distance
                        electric_field(count,:) = electric_field(count,:) &                                                 ! Calculate the HH contribution to the electric field
                                & + Field_Contribution(charges(jhead+2), positions(hid,:), tmp_dr(1:3), SQRT(tmp_dr(4)))
                    ENDIF
                ENDIF
                jhead = cell_list(jhead)        ! Move to the Next Atom in the Cell
            ENDDO
        ENDDO
        ENDDO
        ENDDO
        count = count + 1   ! Increment the count
    ENDDO

END SUBROUTINE