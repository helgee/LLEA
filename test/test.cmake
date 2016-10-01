add_executable(test_assertions test/test_assertions.f90)
target_link_libraries(test_assertions llea)
add_test(TEST_ASSERTIONS test_assertions)

add_executable(test_exceptions test/test_exceptions.f90)
target_link_libraries(test_exceptions llea)
add_test(TEST_EXCEPTIONS test_exceptions)
#
# add_executable(test_integrator test/test_integrator.f90)
# target_link_libraries(test_integrator imdex)
# add_test(TEST_INTEGRATOR test_integrator)
#
# add_executable(test_parallel_integrator test/test_parallel_integrator.f90)
# target_link_libraries(test_parallel_integrator imdex)
# add_test(TEST_PARALLEL_INTEGRATOR test_parallel_integrator)
#
# add_executable(test_time test/test_time.f90)
# target_link_libraries(test_time imdex)
# add_test(TEST_TIME test_time)
#
# add_executable(test_optimizer test/test_optimizer.f90)
# target_link_libraries(test_optimizer imdex)
# add_test(TEST_OPTIMIZER test_optimizer)
#
# add_executable(test_arrays test/test_arrays.f90)
# target_link_libraries(test_arrays imdex)
# add_test(TEST_ARRAYS test_arrays)
#
# add_executable(test_de test/test_de.f90)
# target_link_libraries(test_de imdex)
# add_test(TEST_DE test_de)
#
# add_executable(test_trajectories test/test_trajectories.f90)
# target_link_libraries(test_trajectories imdex)
# add_test(TEST_TRAJECTORIES test_trajectories)
#
# add_executable(test_util test/test_util.f90)
# target_link_libraries(test_util imdex)
# add_test(TEST_UTIL test_util)
