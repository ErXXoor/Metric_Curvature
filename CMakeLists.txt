cmake_minimum_required(VERSION 3.2 FATAL_ERROR)
project(MetricCurvature)

list(PREPEND CMAKE_MODULE_PATH ${CMAKE_CURRENT_SOURCE_DIR}/cmake)
include(libigl)
igl_include(glfw)
igl_include(copyleft comiso)
add_subdirectory(src)