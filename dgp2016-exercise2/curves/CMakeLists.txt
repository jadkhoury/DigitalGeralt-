get_filename_component(EXERCISENAME ${CMAKE_CURRENT_LIST_DIR} NAME)

file(GLOB_RECURSE SOURCES "*.cpp")
file(GLOB_RECURSE HEADERS "*.h")
file(GLOB_RECURSE SHADERS "*.glsl")

add_executable (smoothing curves.cpp)
target_link_libraries(smoothing ${LIBRARIES})

add_executable (reconstruction reconstruction.cpp)
target_link_libraries(reconstruction ${LIBRARIES})