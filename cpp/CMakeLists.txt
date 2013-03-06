set(BRIGHT_INCLUDE_DIRS ${CMAKE_CURRENT_SOURCE_DIR} ${BRIGHT_INCLUDE_DIRS} 
    PARENT_SCOPE)

set(PYNE_LINK_LIBS pyne pyne_nucname pyne_data pyne_material)

include_directories(${CMAKE_CURRENT_SOURCE_DIR})

# bright
add_library(bright_bright bright.cpp)
set_target_properties(bright_bright PROPERTIES 
                      LIBRARY_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}/bright/lib")
target_link_libraries(bright_bright ${PYNE_LINK_LIBS} hdf5 hdf5_hl hdf5_cpp hdf5_hl_cpp)
install_lib(bright_bright)

# fccomp
add_library(bright_fccomp fccomp.cpp)
set_target_properties(bright_fccomp PROPERTIES 
                      LIBRARY_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}/bright/lib")
target_link_libraries(bright_fccomp bright_bright ${PYNE_LINK_LIBS})
install_lib(bright_fccomp)

# enrichment
add_library(bright_enrichment_parameters enrichment_parameters.cpp)
set_target_properties(bright_enrichment_parameters PROPERTIES 
                      LIBRARY_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}/bright/lib")
target_link_libraries(bright_enrichment_parameters bright_bright ${PYNE_LINK_LIBS})
install_lib(bright_enrichment_parameters)

# enrichment
add_library(bright_enrichment bright_enrichment.cpp)
set_target_properties(bright_enrichment PROPERTIES 
                      LIBRARY_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}/bright/lib")
target_link_libraries(bright_enrichment bright_bright bright_fccomp 
                      bright_enrichment_parameters ${PYNE_LINK_LIBS})
install_lib(bright_enrichment)

# reprocess
add_library(bright_reprocess reprocess.cpp)
set_target_properties(bright_reprocess PROPERTIES 
                      LIBRARY_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}/bright/lib")
target_link_libraries(bright_reprocess bright_bright bright_fccomp ${PYNE_LINK_LIBS})
install_lib(bright_reprocess)

# storage
add_library(bright_storage storage.cpp)
set_target_properties(bright_storage PROPERTIES 
                      LIBRARY_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}/bright/lib")
target_link_libraries(bright_storage bright_bright bright_fccomp ${PYNE_LINK_LIBS})
install_lib(bright_storage)

# reactor parameters
add_library(bright_reactor_parameters reactor_parameters.cpp)
set_target_properties(bright_reactor_parameters PROPERTIES 
                      LIBRARY_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}/bright/lib")
install_lib(bright_reactor_parameters)

# fluence point
add_library(bright_fluence_point fluence_point.cpp)
set_target_properties(bright_fluence_point PROPERTIES 
                      LIBRARY_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}/bright/lib")
install_lib(bright_fluence_point)

# reactor1g
add_library(bright_reactor1g reactor1g.cpp)
set_target_properties(bright_reactor1g PROPERTIES 
                      LIBRARY_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}/bright/lib")
target_link_libraries(bright_reactor1g bright_bright bright_fccomp 
                      bright_reactor_parameters bright_fluence_point ${PYNE_LINK_LIBS})
install_lib(bright_reactor1g)

# light water reactor1g
add_library(bright_light_water_reactor1g light_water_reactor1g.cpp)
set_target_properties(bright_light_water_reactor1g PROPERTIES 
                      LIBRARY_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}/bright/lib")
target_link_libraries(bright_light_water_reactor1g bright_bright bright_fccomp 
                      bright_reactor_parameters bright_fluence_point 
                      bright_reactor1g ${PYNE_LINK_LIBS})
install_lib(bright_light_water_reactor1g)

# fast reactor1g
add_library(bright_fast_reactor1g fast_reactor1g.cpp)
set_target_properties(bright_fast_reactor1g PROPERTIES 
                      LIBRARY_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}/bright/lib")
target_link_libraries(bright_fast_reactor1g bright_bright bright_fccomp 
                      bright_reactor_parameters bright_fluence_point 
                      bright_reactor1g ${PYNE_LINK_LIBS})
install_lib(bright_fast_reactor1g)

# fuel fabrication
add_library(bright_fuel_fabrication fuel_fabrication.cpp)
set_target_properties(bright_fuel_fabrication PROPERTIES 
                      LIBRARY_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}/bright/lib")
target_link_libraries(bright_fuel_fabrication bright_bright bright_fccomp 
                      bright_reactor_parameters bright_fluence_point 
                      bright_reactor1g ${PYNE_LINK_LIBS})
install_lib(bright_fuel_fabrication)

# reactormg
add_library(bright_reactormg reactormg.cpp)
set_target_properties(bright_reactormg PROPERTIES 
                      LIBRARY_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}/bright/lib")
target_link_libraries(bright_reactormg bright_bright bright_fccomp 
                      bright_reactor_parameters bright_fluence_point ${PYNE_LINK_LIBS})
install_lib(bright_reactormg)


# Print include dir
get_property(inc_dirs DIRECTORY PROPERTY INCLUDE_DIRECTORIES)
message("-- C_INCLUDE_PATH for ${CMAKE_CURRENT_SOURCE_DIR}: ${inc_dirs}")
