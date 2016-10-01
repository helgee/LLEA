set(LLEA_SRCS
    src/assertions.f90
    src/bodies.f90
    src/epochs.f90
    src/exceptions.f90
    src/math.f90
    src/states.f90
    src/types.f90
    src/util.f90)

add_library(llea ${LLEA_SRCS})
