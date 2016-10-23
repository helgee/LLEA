cmake_minimum_required(VERSION 2.8)

project(LLEA Fortran)
enable_language(Fortran)
enable_testing()

include(cmake/data.cmake)
include(cmake/flags.cmake)
include(deps/deps.cmake)
include(src/src.cmake)

add_library(llea SHARED ${LLEA_SRCS})
if(CMAKE_VERSION VERSION_GREATER 2.8.12)
    set_target_properties(llea PROPERTIES MACOSX_RPATH ON)
endif()
install(TARGETS llea LIBRARY DESTINATION lib)

include(test/test.cmake)
