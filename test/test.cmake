add_executable(test_assertions test/test_assertions.f90)
target_link_libraries(test_assertions llea)
add_test(TEST_ASSERTIONS test_assertions)

add_executable(test_exceptions test/test_exceptions.f90)
target_link_libraries(test_exceptions llea)
add_test(TEST_EXCEPTIONS test_exceptions)

add_executable(test_ephemerides test/test_ephemerides.f90)
target_link_libraries(test_ephemerides llea)
add_test(TEST_EPHEMERIDES test_ephemerides)

add_executable(test_epochs test/test_epochs.f90)
target_link_libraries(test_epochs llea)
add_test(TEST_EPOCHS test_epochs)

add_executable(test_math test/test_math.f90)
target_link_libraries(test_math llea)
add_test(TEST_MATH test_math)

add_executable(test_propagators test/test_propagators.f90)
target_link_libraries(test_propagators llea)
add_test(TEST_PROPAGATORS test_propagators)

add_executable(test_rotations test/test_rotations.f90)
target_link_libraries(test_rotations llea)
add_test(TEST_ROTATIONS test_rotations)

add_executable(test_iau test/test_iau.f90)
target_link_libraries(test_iau llea)
add_test(TEST_IAU test_iau)

add_executable(test_states test/test_states.f90)
target_link_libraries(test_states llea)
add_test(TEST_STATES test_states)

add_executable(test_trajectories test/test_trajectories.f90)
target_link_libraries(test_trajectories llea)
add_test(TEST_TRAJECTORIES test_trajectories)
