SUBROUTINE Cell_Internal_Distance(ihead_init, jhead_init, &
                                & icell_list, jcell_list, &
                                & ipositions, jpositions, &
                                box_lengths, r_cutoff_sq, inner_count, local_count)
! This function calcluates the distance between atoms within two cells.
!
! Parameters:
! -----------
! ihead_init : INTEGER
!       The head of the outer list
! jhead_init : INTEGER
!       The head of the inner list
! cell_list : INTEGER, DIMENSION(:)
!       The linked-list for the cells.
! positions : REAL, DIMENSION(:,:)
!       The positions list.
! box_lengths : REAL, DIMENSION(3)
!       Array of box positions in angstroms.
! rc_sq : REAL
!       Cutoff for the distance calculation
! dr : REAL, DIMENSION(:)
!       Array of distances.
! dr_ids_1 : INTEGER, DIMENSION(:)
!       List of ids of atoms in dr
! dr_ids_2 : INTEGER, DIMENSION(:)
!       List of ids of atoms in dr
!
    INTEGER,                  INTENT(IN)  :: ihead_init, jhead_init
    INTEGER, DIMENSION(:),    INTENT(IN)  :: icell_list, jcell_list
    REAL,    DIMENSION(:,:),  INTENT(IN)  :: ipositions, jpositions
    REAL,    DIMENSION(3),    INTENT(IN)  :: box_lengths
    REAL,                     INTENT(IN)  :: r_cutoff_sq
    INTEGER,                  INTENT(OUT) :: inner_count
    INTEGER, DIMENSION(:),     INTENT(INOUT) :: local_count
    
    REAL, DIMENSION(3) :: distance_vector

    INTEGER :: ihead, jhead
    REAL    :: distance, distance_sq

    inner_count = 0
    ihead = ihead_init
    DO WHILE (ihead /= 0)
        jhead = jhead_init
        DO WHILE(jhead /= 0)
            distance_vector(:) = ipositions(ihead,:) - jpositions(jhead,:)
            distance_vector(:) = distance_vector(:) - box_lengths(:)*ANINT(distance_vector(:)/box_lengths(:))
            distance_sq = SUM((distance_vector(:))**2.)

            IF (distance_sq < r_cutoff_sq) THEN
                inner_count = inner_count + 1
                local_count(ihead) = local_count(ihead) + 1.0
            ENDIF
            jhead = jcell_list(jhead)
        ENDDO
        ihead = icell_list(ihead)
    ENDDO

END SUBROUTINE Cell_Internal_Distance

SUBROUTINE Cell_Distance(ipositions, jpositions, box_lengths, r_cutoff, local_count)
!
    IMPLICIT NONE
    REAL, DIMENSION(:,:),       INTENT(IN) :: ipositions, jpositions
    REAL, DIMENSION(3),         INTENT(IN) :: box_lengths
    REAL,                       INTENT(IN) :: r_cutoff
    
    INTEGER, DIMENSION(3)                  :: cell_nbins
    INTEGER, DIMENSION(:),     ALLOCATABLE :: icell_list, jcell_list ! Dimensions of # positions
    INTEGER, DIMENSION(:,:),   ALLOCATABLE :: cell_map  ! Dimensions of # bins
    INTEGER, DIMENSION(:,:,:), ALLOCATABLE ::  ill_head, jll_head ! Dimensions
    INTEGER, DIMENSION(:), INTENT(OUT) :: local_count


    INTEGER :: i, j, k, l, m, n, ii, jj, kk
    INTEGER :: ihead, jhead
    INTEGER :: max_nbins
    INTEGER :: count, tmpcount

    REAL :: r_cutoff_sq

    r_cutoff_sq = r_cutoff*r_cutoff

    local_count = 0
    ! Part 1: Construct Grid, and allocate arrays.
    ALLOCATE(icell_list(SIZE(ipositions,1)))
    ALLOCATE(jcell_list(SIZE(jpositions,1)))
    icell_list = 0; jcell_list = 0

    CALL Construct_Grid(box_lengths, r_cutoff, cell_nbins, cell_map)

    ALLOCATE(ill_head(cell_nbins(1), cell_nbins(2), cell_nbins(3)))
    ALLOCATE(jll_head(cell_nbins(1), cell_nbins(2), cell_nbins(3)))
    ill_head = 0; jll_head = 0

    ! Part 2: Build linked list
    CALL Linked_List_Build(ipositions, box_lengths, r_cutoff, cell_nbins, ill_head, icell_list, cell_map)
    CALL Linked_List_Build(jpositions, box_lengths, r_cutoff, cell_nbins, jll_head, jcell_list, cell_map)

    count = 0

    ! Part 3: Calcualte Distances
    DO i=1, cell_nbins(1)
    DO j=1, cell_nbins(2)
    DO k=1, cell_nbins(3)
        DO l=-2,2
        DO m=-2,2
        DO n=-2,2
            ii = cell_map(1, i + l + 2)
            jj = cell_map(2, j + m + 2)
            kk = cell_map(3, k + n + 2)
            
            ihead = ill_head(i,  j,  k)
            jhead = jll_head(ii, jj, kk)

            CALL Cell_Internal_Distance(ihead, jhead, &
                                    &   icell_list, jcell_list, &
                                    &   ipositions, jpositions, &
                                    &   box_lengths, r_cutoff_sq, tmpcount, local_count)
            count = count + tmpcount

        ENDDO !l
        ENDDO !m
        ENDDO !n
    ENDDO !j
    ENDDO !j
    ENDDO !k

    DEALLOCATE(cell_map)
    DEALLOCATE(icell_list); DEALLOCATE(jcell_list)
    DEALLOCATE(ill_head); DEALLOCATE(jll_head)

END SUBROUTINE Cell_Distance

SUBROUTINE EF_List(iposition, jpositions, jcharges, box_lengths, r_cutoff, &
                &  cell_list, ll_head, cell_nbins, cell_map, jcontributes, electric_field)
    IMPLICIT NONE
    REAL, DIMENSION(3),         INTENT(IN) :: iposition
    REAL, DIMENSION(:,:),       INTENT(IN) :: jpositions
    REAL,                       INTENT(IN) :: jcharges
    REAL, DIMENSION(3),         INTENT(IN) :: box_lengths
    REAL,                       INTENT(IN) :: r_cutoff
    INTEGER, DIMENSION(3),      INTENT(IN) :: cell_nbins
    INTEGER, DIMENSION(:),     ALLOCATABLE, INTENT(IN) :: cell_list ! Dimensions of # positions
    INTEGER, DIMENSION(:,:),   ALLOCATABLE, INTENT(IN) :: cell_map  ! Dimensions of # bins
    INTEGER, DIMENSION(:,:,:), ALLOCATABLE, INTENT(IN) :: ll_head ! Dimensions


    INTEGER :: l, m, n, ii, jj, kk
    INTEGER :: head, count, tmpcount
    REAL :: r_cutoff_sq
    REAL, DIMENSION(3) :: distance_vector
    REAL :: r_distance, distance_sq
    
    INTEGER,    DIMENSION(3)    :: location_on_grid
    INTEGER,    DIMENSION(:)    :: jcontributes
    REAL,       DIMENSION(:)    :: electric_field




    r_cutoff_sq = r_cutoff*r_cutoff
    location_on_grid = Assign_To_Cell_Grid(iposition, cell_nbins, box_lengths)

    DO l=-2,2
    DO m=-2,2
    DO n=-2,2
        ii = cell_map(1, location_on_grid(1) + l + 2)
        jj = cell_map(2, location_on_grid(2) + m + 2)
        kk = cell_map(3, location_on_grid(3) + n + 2)
        
        head = ll_head(ii, jj, kk)

        DO WHILE(head /= 0)
            distance_vector(:) = iposition(:) - jpositions(head,:)
            distance_vector(:) = distance_vector(:) - box_lengths(:)*ANINT(distance_vector(:)/box_lengths(:))
            distance_sq = SUM((distance_vector(:))**2.)

            IF (distance_sq < r_cutoff_sq) THEN
                jcontributes(head) = 1
                electric_field(:) = electric_field(:) + Electric_Field_Component(jcharges, distance_vector, r_distance)
            ENDIF
            head = cell_list(head)
        END DO

        count = count + tmpcount
    ENDDO !l
    ENDDO !m
    ENDDO !n

END SUBROUTINE EF_LIST

