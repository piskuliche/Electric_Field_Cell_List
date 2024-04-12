add_library(data_module data_mods.f90)
add_library(alloc_handler allocator.f90)

target_link_libraries(alloc_handler data_module)