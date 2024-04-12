! ****************************************************************************************************
!  This module stores the data for the program.
! 
! Copyright (C) 2024 Zeke Piskulich
! ****************************************************************************************************


MODULE InputData

    IMPLICIT NONE
    INTEGER, PARAMETER :: ndim = 3
    INTEGER, PARAMETER :: max_config = 100
    INTEGER :: nconfig
    INTEGER :: nmoltypes, which_is_water
    REAL    :: r_cutoff

    INTEGER,    DIMENSION(10) :: nmols, natoms
    REAL,       DIMENSION(10,2000) :: chrgs
    INTEGER,    DIMENSION(10,2000) :: atps
    REAL,       DIMENSION(3) :: box_lengths
    CHARACTER(LEN=10), DIMENSION(10) :: molnames

    INTEGER :: wstart, wend, number_of_atoms

    REAL, PARAMETER :: angperau = 0.52917721092d0

END MODULE

MODULE TimingData
    
    IMPLICIT NONE
    REAL :: read_time, write_time, calc_time, total_time

END MODULE TimingData