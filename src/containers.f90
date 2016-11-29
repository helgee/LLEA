!
! Copyright (c) 2016 Helge Eichhorn
!
! This Source Code Form is subject to the terms of the Mozilla Public
! License, v. 2.0. If a copy of the MPL was not distributed with this
! file, You can obtain one at http://mozilla.org/MPL/2.0/.
!
module containers

use bodies, only: body
use states, only: state, framelen
use types, only: dp

implicit none

type parameters
    type(state) :: s0
    character(len=framelen) :: frame
    type(body) :: center
    real(dp), dimension(:), allocatable :: con
    integer, dimension(:), allocatable :: icomp
end type parameters

interface parameters
    module procedure parameters_init
end interface parameters

contains

function parameters_init(s0, frame, center) result(p)
    type(state), intent(in) :: s0
    character(len=framelen), intent(in) :: frame
    type(body), intent(in) :: center
    type(parameters) :: p

    p%s0 = s0
    p%frame = frame
    p%center = center
end function parameters_init

end module containers
