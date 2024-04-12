add_library(io_module io_module.f90)

target_link_libraries(io_module data_module ${HDF5_LIBRARIES})