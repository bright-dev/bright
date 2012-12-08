# Python extension modules
include_directories(${BRIGHT_INCLUDE_DIRS})

# jsoncpp
set_source_files_properties(${PROJECT_SOURCE_DIR}/bright/bright_config.pyx
                            PROPERTIES CYTHON_IS_CXX TRUE)
cython_add_module(bright_config bright_config.pyx)
target_link_libraries(bright_config.pyx bright_bright)
