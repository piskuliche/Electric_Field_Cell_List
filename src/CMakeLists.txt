add_executable(get_efield electric_field.f90 coordinates.f90 loopdistance.f90)

target_link_libraries(get_efield cell_list)

add_subdirectory(cell-list-decomposition)
add_subdirectory(efield_code)