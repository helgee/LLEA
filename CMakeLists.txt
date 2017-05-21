cmake_minimum_required(VERSION 2.8)

project(LLEA Fortran)
enable_language(Fortran)
enable_testing()

# if(WIN32)
	# set(CMAKE_WINDOWS_EXPORT_ALL_SYMBOLS true CACHE BOOL "blob")
# endif()

include(cmake/data.cmake)
include(cmake/flags.cmake)
include(deps/deps.cmake)
include(src/src.cmake)

add_library(llea ${LLEA_SRCS})
target_link_libraries(llea dierckx)
if(CMAKE_VERSION VERSION_GREATER 2.8.12)
    set_target_properties(llea dierckx PROPERTIES MACOSX_RPATH ON)
endif()
install(
	TARGETS llea
	LIBRARY DESTINATION lib
	ARCHIVE DESTINATION lib
)

include(test/test.cmake)
