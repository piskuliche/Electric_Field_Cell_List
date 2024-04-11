# Electric_Field_Cell_List

[![License](https://img.shields.io/badge/license-MIT-blue.svg)](https://opensource.org/licenses/MIT)

## Description

This repository contains the code for calculating electric field using a cell list data structure. The code is written in fortran and is designed to be efficient and scalable for large-scale simulations.

## Features

- Efficient calculation of electric field using a cell list data structure
- Support for OH electric field calculations of water
- Support for other molecules in an extensible way
- Writes output in an efficient, HDF5 binary archive format

## Installation

1. Clone the repository:

    ```bash
    git clone https://github.com/your-username/Electric_Field_Cell_List.git
    ```

2. Make a Build directory

    ```bash
    mkdir build
    cd build
    ```

3. Run CMAKE

    ```bash
    cmake ..
    make
    ```

## Usage

This code is run by using the executable (new_efield) alongside a trajectory.

Right now, the default trajectory format and filename are xyz, and traj.xyz. There isn't a way to change this right now, but eventually this could be added as an input feature.

The code takes a few other files:

1) field_input.in

    ```
    # Number of Molecules, Number of Configurations, which_is_water
    2 100 2
    # cutoff distance, L(1), L(2), L(3)
    8.000   80.22729   80.22729   94.99541
    # molecule 1, molecule 2
    popc water
    # nmols(1), nmols(2)
    200 12305
    # DEPRECATED FEATURE
    20
    ```

2) foo.in: a molecular input file. Format looks something like this.

    ```
    3
    -0.834  1
     0.417  2
     0.417  2
    ```

    The first line is the number of atoms in a single molecule of type foo (in this case, water - this would be water.in)

    The next lines are the charge and a user-defined atom type. This should start from 1 and work upwards. If you have more than one type of molecule, then those can all number from 1 for the atom types (a line in the code adds to a max_type in order to avoid problems. Just don't set this to zero.)
## Contributing

Contributions are welcome! If you find any issues or have suggestions for improvements, please open an issue or submit a pull request.

## License

This project is licensed under the [MIT License](LICENSE).

