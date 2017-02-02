!
! Copyright (c) 2016 Helge Eichhorn
!
! This Source Code Form is subject to the terms of the Mozilla Public
! License, v. 2.0. If a copy of the MPL was not distributed with this
! file, You can obtain one at http://mozilla.org/MPL/2.0/.
!
module llea

use iso_c_binding, only: c_int

use assertions
use bodies
use constants
use containers
use dopri
use ephemerides
use epochs
use events
use exceptions
use forces
use integrators
use interfaces
use libration
use math
use naif
use propagators
use rotations
use splines
use states
use trajectories
use types
use util

implicit none

contains

function init_llea(ephempath) result(code)
    character(len=*), intent(in), optional :: ephempath
    integer :: code
    call init_constants
    if (present(ephempath)) then
        call init_ephemeris(path=ephempath)
    else
        call init_ephemeris
    end if
    code = 0
end function init_llea

function init_llea_c() result(code) bind(c)
    integer(c_int) :: code

    code = int(init_llea(), c_int)
end function init_llea_c

end module llea
