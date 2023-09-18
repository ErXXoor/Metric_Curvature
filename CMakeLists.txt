cmake_minimum_required(VERSION 3.2 FATAL_ERROR)
project(MetricCurvature)

set(CMAKE_CXX_STANDARD 17)
list(PREPEND CMAKE_MODULE_PATH ${CMAKE_CURRENT_SOURCE_DIR}/cmake)
include(libigl)
igl_include(glfw)

add_subdirectory(src)