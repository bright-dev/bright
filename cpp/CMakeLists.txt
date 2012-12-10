set(BRIGHT_INCLUDE_DIRS ${CMAKE_CURRENT_SOURCE_DIR} ${BRIGHT_INCLUDE_DIRS} 
    PARENT_SCOPE)

set(PYNE_LINK_LIBS pyne pyne_nucname pyne_data pyne_material)

include_directories(${CMAKE_CURRENT_SOURCE_DIR})

# bright
add_library(bright_bright bright.cpp)
set_target_properties(bright_bright PROPERTIES 
                      LIBRARY_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}/bright/lib")
target_link_libraries(bright_bright ${PYNE_LINK_LIBS})
install_lib(bright_bright)



# Print include dir
get_property(inc_dirs DIRECTORY PROPERTY INCLUDE_DIRECTORIES)
message("-- C_INCLUDE_PATH for ${CMAKE_CURRENT_SOURCE_DIR}: ${inc_dirs}")

