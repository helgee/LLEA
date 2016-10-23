set(SPK_430 "${PROJECT_SOURCE_DIR}/data/de430.bsp")
set(SPK_405 "${PROJECT_SOURCE_DIR}/data/de405.bsp")
set(TEST_430 "${PROJECT_SOURCE_DIR}/data/testpo.430")
set(TEST_405 "${PROJECT_SOURCE_DIR}/data/testpo.405")

if(NOT EXISTS ${SPK_430})
    file(DOWNLOAD "http://naif.jpl.nasa.gov/pub/naif/generic_kernels/spk/planets/de430.bsp"
        ${SPK_430} SHOW_PROGRESS)
endif()
if(NOT EXISTS ${SPK_405})
    file(DOWNLOAD "http://naif.jpl.nasa.gov/pub/naif/generic_kernels/spk/planets/a_old_versions/de405.bsp"
        ${SPK_405} SHOW_PROGRESS)
endif()
if(NOT EXISTS ${TEST_430})
    file(DOWNLOAD "ftp://ssd.jpl.nasa.gov/pub/eph/planets/test-data/430/testpo.430"
        ${TEST_430} SHOW_PROGRESS)
endif()
if(NOT EXISTS ${TEST_405})
    file(DOWNLOAD "ftp://ssd.jpl.nasa.gov/pub/eph/planets/test-data/testpo.405"
        ${TEST_405} SHOW_PROGRESS)
endif()
