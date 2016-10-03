cmake_minimum_required(VERSION 3.0)

project(LLEA Fortran)
enable_language(Fortran)
enable_testing()

include(cmake/flags.cmake)
include(deps/deps.cmake)
include(src/src.cmake)

add_library(llea SHARED ${LLEA_SRCS})
install(TARGETS llea LIBRARY DESTINATION lib)

include(test/test.cmake)
