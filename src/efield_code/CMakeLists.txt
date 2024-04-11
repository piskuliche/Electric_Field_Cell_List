add_executable(new_efield new_field.f90 coordinates.f90 read_input.f90 data_mods.f90 allocator.f90 field_functions.f90 field_linked_list.f90 read_write_hdf5.f90)


target_link_libraries(new_efield cell_list ${HDF5_LIBRARIES})