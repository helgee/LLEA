!
! Copyright (c) 2016 Helge Eichhorn
!
! This Source Code Form is subject to the terms of the Mozilla Public
! License, v. 2.0. If a copy of the MPL was not distributed with this
! file, You can obtain one at http://mozilla.org/MPL/2.0/.
!
! This file incorporates work covered by the following copyright and
! permission notice:
!
module trajectories

use epochs, only: epochdelta, operator (+)
use states, only: state
use types, only: dp

implicit none

integer, parameter :: fieldlen = 12

type listnode
    type(listnode), pointer :: next => null()
    real(dp), dimension(:), allocatable :: val
end type listnode

type trajectory
    type(state) :: initial_state
    type(state) :: final_state
    real(dp), dimension(:), allocatable :: tindex
    character(len=fieldlen), dimension(:), allocatable :: fields
    real(dp), dimension(:,:), allocatable :: vectors
end type trajectory

interface trajectory
    module procedure new_trajectory_array
end interface trajectory

contains

function new_trajectory_array(s0, t, arr, fields) result(tra)
    type(state), intent(in) :: s0
    real(dp), dimension(:), intent(in) :: t
    real(dp), dimension(:,:), intent(in) :: arr
    character(len=fieldlen), dimension(:), intent(in), optional :: fields
    type(trajectory) :: tra

    integer :: n

    n = size(t)
    tra%fields = ["x ", "y ", "z ", "vx", "vy", "vz"]
    if (present(fields)) tra%fields = fields
    tra%initial_state = s0
    tra%vectors = arr
    tra%tindex = t
    tra%final_state = state(s0%ep + epochdelta(seconds=t(n)), arr(n, :), s0%frame, s0%center)
end function new_trajectory_array

end module trajectories
