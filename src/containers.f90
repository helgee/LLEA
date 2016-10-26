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
    real(dp), dimension(:,:), allocatable :: con
    integer, dimension(:,:), allocatable :: icomp
end type parameters

end module containers
