MODULE InputData

    IMPLICIT NONE
    INTEGER, PARAMETER :: ndim = 3
    INTEGER, PARAMETER :: max_config = 100
    INTEGER :: nconfig
    INTEGER :: nmoltypes, which_is_water
    INTEGER :: nsamples ! Deprecate
    REAL    :: r_cutoff

    INTEGER,    DIMENSION(10) :: nmols, natoms
    REAL,       DIMENSION(10,2000) :: chrgs
    INTEGER,    DIMENSION(10,2000) :: atps
    REAL,       DIMENSION(3) :: box_lengths
    CHARACTER(LEN=10), DIMENSION(10) :: molnames

    INTEGER :: wstart, wend, number_of_atoms

END MODULE

MODULE Units
    IMPLICIT NONE
    REAL, PARAMETER :: angperau = 0.52917721092d0
END MODULE UNITS

