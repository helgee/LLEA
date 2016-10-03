module llea

use iso_c_binding, only: c_int

use assertions
use bodies
use constants
use epochs
use exceptions
use math
use rotations
use states
use types
use util

implicit none

contains

function init_llea() result(code)
    integer :: code
    call init_constants
    code = 0
end function init_llea

function init_llea_c() result(code) bind(c)
    integer(c_int) :: code

    code = int(init_llea(), c_int)
end function init_llea_c

end module llea
