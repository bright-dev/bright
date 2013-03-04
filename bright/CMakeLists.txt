# Python extension modules
include_directories(${BRIGHT_INCLUDE_DIRS})

# Print include dir
get_property(inc_dirs DIRECTORY PROPERTY INCLUDE_DIRECTORIES)
message("-- C_INCLUDE_PATH for ${CMAKE_CURRENT_SOURCE_DIR}: ${inc_dirs}")

# type converters
set_source_files_properties(${PROJECT_SOURCE_DIR}/bright/typeconverters.pyx
                            PROPERTIES CYTHON_IS_CXX TRUE)
cython_add_module(typeconverters typeconverters.pyx)
target_link_libraries(typeconverters pyne_nucname pyne_material)

# bright config
set_source_files_properties(${PROJECT_SOURCE_DIR}/bright/bright_config.pyx
                            PROPERTIES CYTHON_IS_CXX TRUE)
cython_add_module(bright_config bright_config.pyx)
target_link_libraries(bright_config bright_bright)

# fccomp
set_source_files_properties(${PROJECT_SOURCE_DIR}/bright/fccomp.pyx
                            PROPERTIES CYTHON_IS_CXX TRUE)
cython_add_module(fccomp fccomp.pyx)
target_link_libraries(fccomp bright_fccomp)

# enrichment
set_source_files_properties(${PROJECT_SOURCE_DIR}/bright/enrichment.pyx
                            PROPERTIES CYTHON_IS_CXX TRUE)
cython_add_module(enrichment enrichment.pyx)
target_link_libraries(enrichment bright_enrichment)

# reprocess
set_source_files_properties(${PROJECT_SOURCE_DIR}/bright/reprocess.pyx
                            PROPERTIES CYTHON_IS_CXX TRUE)
cython_add_module(reprocess reprocess.pyx)
target_link_libraries(reprocess bright_reprocess)

# storage
set_source_files_properties(${PROJECT_SOURCE_DIR}/bright/storage.pyx
                            PROPERTIES CYTHON_IS_CXX TRUE)
cython_add_module(storage storage.pyx)
target_link_libraries(storage bright_storage)

# reactor_parameters
set_source_files_properties(${PROJECT_SOURCE_DIR}/bright/reactor_parameters.pyx
                            PROPERTIES CYTHON_IS_CXX TRUE)
cython_add_module(reactor_parameters reactor_parameters.pyx)
target_link_libraries(reactor_parameters bright_reactor_parameters)

# fluence_point
set_source_files_properties(${PROJECT_SOURCE_DIR}/bright/fluence_point.pyx
                            PROPERTIES CYTHON_IS_CXX TRUE)
cython_add_module(fluence_point fluence_point.pyx)
target_link_libraries(fluence_point bright_fluence_point)

# reactor1g
set_source_files_properties(${PROJECT_SOURCE_DIR}/bright/reactor1g.pyx
                            PROPERTIES CYTHON_IS_CXX TRUE)
cython_add_module(reactor1g reactor1g.pyx)
target_link_libraries(reactor1g bright_reactor1g)

# light_water_reactor1g
set_source_files_properties(${PROJECT_SOURCE_DIR}/bright/light_water_reactor1g.pyx
                            PROPERTIES CYTHON_IS_CXX TRUE)
cython_add_module(light_water_reactor1g light_water_reactor1g.pyx)
target_link_libraries(light_water_reactor1g bright_light_water_reactor1g)

# fast_reactor1g
set_source_files_properties(${PROJECT_SOURCE_DIR}/bright/fast_reactor1g.pyx
                            PROPERTIES CYTHON_IS_CXX TRUE)
cython_add_module(fast_reactor1g fast_reactor1g.pyx)
target_link_libraries(fast_reactor1g bright_fast_reactor1g)

# fuel_fabrication
set_source_files_properties(${PROJECT_SOURCE_DIR}/bright/fuel_fabrication.pyx
                            PROPERTIES CYTHON_IS_CXX TRUE)
cython_add_module(fuel_fabrication fuel_fabrication.pyx)
target_link_libraries(fuel_fabrication bright_fuel_fabrication)

# reactormg
set_source_files_properties(${PROJECT_SOURCE_DIR}/bright/reactormg.pyx
                            PROPERTIES CYTHON_IS_CXX TRUE)
cython_add_module(reactormg reactormg.pyx)
target_link_libraries(reactormg bright_reactormg)

