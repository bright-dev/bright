set(BRIGHT_INCLUDE_DIRS ${CMAKE_CURRENT_SOURCE_DIR} ${BRIGHT_INCLUDE_DIRS} 
    PARENT_SCOPE)

include_directories(${CMAKE_CURRENT_SOURCE_DIR})

# bright
add_library(bright
    utils.cpp
    fccomp.cpp
    enrichment_parameters.cpp
    bright_enrichment.cpp
    reprocess.cpp
    storage.cpp
    reactor_parameters.cpp
    fluence_point.cpp
    reactor1g.cpp
    light_water_reactor1g.cpp
    fast_reactor1g.cpp
    fuel_fabrication.cpp
    reactormg.cpp
    )
target_link_libraries(bright pyne hdf5 hdf5_hl hdf5_cpp hdf5_hl_cpp)
install_lib(bright)

# Print include dir
get_property(inc_dirs DIRECTORY PROPERTY INCLUDE_DIRECTORIES)
message("-- C_INCLUDE_PATH for ${CMAKE_CURRENT_SOURCE_DIR}: ${inc_dirs}")
