set(SPK_430 "${PROJECT_BINARY_DIR}/data/de430.bsp")
set(SPK_405 "${PROJECT_BINARY_DIR}/data/de405.bsp")
set(TEST_430 "${PROJECT_BINARY_DIR}/data/testpo.430")
set(TEST_405 "${PROJECT_BINARY_DIR}/data/testpo.405")

if(NOT EXISTS ${SPK_430})
    file(DOWNLOAD "http://naif.jpl.nasa.gov/pub/naif/generic_kernels/spk/planets/de430.bsp" ${SPK_430} SHOW_PROGRESS)
endif()
if(NOT EXISTS ${SPK_405})
    file(DOWNLOAD "http://naif.jpl.nasa.gov/pub/naif/generic_kernels/spk/planets/a_old_versions/de405.bsp" ${SPK_405} SHOW_PROGRESS)
endif()
if(NOT EXISTS ${TEST_430})
    file(DOWNLOAD "ftp://ssd.jpl.nasa.gov/pub/eph/planets/test-data/430/testpo.430" ${TEST_430} SHOW_PROGRESS)
endif()
if(NOT EXISTS ${TEST_405})
    file(DOWNLOAD "ftp://ssd.jpl.nasa.gov/pub/eph/planets/test-data/testpo.405" ${TEST_405} SHOW_PROGRESS)
endif()

add_executable(test_assertions test/test_assertions.f90)
target_link_libraries(test_assertions llea)
add_test(TEST_ASSERTIONS test_assertions)

add_executable(test_exceptions test/test_exceptions.f90)
target_link_libraries(test_exceptions llea)
add_test(TEST_EXCEPTIONS test_exceptions)

add_executable(test_ephemeris test/test_ephemeris.f90)
target_link_libraries(test_ephemeris llea)
add_test(TEST_EPHEMERIS test_ephemeris)

add_executable(test_epochs test/test_epochs.f90)
target_link_libraries(test_epochs llea)
add_test(TEST_EPOCHS test_epochs)

add_executable(test_math test/test_math.f90)
target_link_libraries(test_math llea)
add_test(TEST_MATH test_math)

add_executable(test_rotations test/test_rotations.f90)
target_link_libraries(test_rotations llea)
add_test(TEST_ROTATIONS test_rotations)

add_executable(test_iau test/test_iau.f90)
target_link_libraries(test_iau llea)
add_test(TEST_IAU test_iau)

add_executable(test_states test/test_states.f90)
target_link_libraries(test_states llea)
add_test(TEST_STATES test_states)
