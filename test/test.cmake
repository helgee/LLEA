add_executable(test_assertions test/test_assertions.f90)
target_link_libraries(test_assertions llea)
add_test(TEST_ASSERTIONS test_assertions)

add_executable(test_exceptions test/test_exceptions.f90)
target_link_libraries(test_exceptions llea)
add_test(TEST_EXCEPTIONS test_exceptions)

add_executable(test_epochs test/test_epochs.f90)
target_link_libraries(test_epochs llea)
add_test(TEST_EPOCHS test_epochs)

add_executable(test_math test/test_math.f90)
target_link_libraries(test_math llea)
add_test(TEST_MATH test_math)

add_executable(test_rotations test/test_rotations.f90)
target_link_libraries(test_rotations llea)
add_test(TEST_ROTATIONS test_rotations)
