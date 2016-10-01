cmake_minimum_required(VERSION 2.6)

project(LLEA Fortran)
enable_language(Fortran)
enable_testing()

include(cmake/flags.cmake)
include(deps/deps.cmake)
include(src/src.cmake)
include(test/test.cmake)

# target_link_libraries(imdex
#     snopt
#     orderpack
#     FoX_common
#     FoX_dom
#     FoX_fsys
#     FoX_sax
#     FoX_utils
#     FoX_wxml)
# install(TARGETS imdex snopt orderpack
#         RUNTIME DESTINATION bin
#         LIBRARY DESTINATION lib
#         ARCHIVE DESTINATION lib/static)
