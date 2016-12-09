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

type tranode
    type(tranode), pointer :: next => null()
    real(dp) :: t
    real(dp), dimension(:), allocatable :: y
end type tranode

type trajectory
    type(state) :: initial_state
    type(state) :: final_state
    character(len=fieldlen), dimension(:), allocatable :: fields
    real(dp), dimension(:), allocatable :: tindex
    real(dp), dimension(:,:), allocatable :: vectors
end type trajectory

interface trajectory
    module procedure new_trajectory_array
    ! Gfortran chokes on this
    ! module procedure new_trajectory_tranode
end interface trajectory

interface tranode
    module procedure init_tranode
end interface tranode

private

public :: tranode, trajectory, length, add_node, save_trajectory

contains

function init_tranode(t, y) result(node)
    real(dp), intent(in) :: t
    real(dp), dimension(:), intent(in) :: y
    type(tranode) :: node
    node%t = t
    node%y = y
end function init_tranode

function new_trajectory_array(s0, t, arr, fields) result(tra)
    type(state), intent(in) :: s0
    real(dp), dimension(:), intent(in) :: t
    real(dp), dimension(:,:), intent(in) :: arr
    character(len=*), dimension(:), intent(in), optional :: fields
    type(trajectory) :: tra

    integer :: n

    n = size(t)
    call setfields(tra, ["x ", "y ", "z ", "vx", "vy", "vz"])
    if (present(fields)) then
        if (allocated(tra%fields)) deallocate(tra%fields)
        call setfields(tra, fields)
    end if
    tra%initial_state = s0
    tra%vectors = arr
    tra%tindex = t
    tra%final_state = state(s0%ep + epochdelta(seconds=t(n)), arr(n, :), s0%frame, s0%center)
end function new_trajectory_array

function save_trajectory(s0, node, fields, delete) result(tra)
    type(state), intent(in) :: s0
    type(tranode), intent(in), target :: node
    character(len=*), dimension(:), intent(in), optional :: fields
    logical, intent(in), optional :: delete
    type(trajectory) :: tra

    integer :: n
    integer :: i
    logical :: delete_
    type(tranode), pointer :: current

    delete_ = .true.
    if (present(delete)) delete_ = delete

    n = length(node)
    call setfields(tra, ["x ", "y ", "z ", "vx", "vy", "vz"])
    if (present(fields)) then
        call setfields(tra, fields)
    end if

    allocate(tra%tindex(n))
    allocate(tra%vectors(n, size(node%y)))
    tra%tindex(1) = node%t
    tra%vectors(1,:) = node%y

    i = 1
    current => node%next
    do
        if (associated(current)) then
            i = i + 1
            tra%tindex(i) = current%t
            tra%vectors(i,:) = current%y
            current => current%next
            if (delete_) then
                deallocate(current)
                nullify(current)
            end if
        else
            exit
        end if
    end do

    tra%initial_state = s0
    tra%final_state = state(s0%ep + epochdelta(seconds=tra%tindex(n)), tra%vectors(n, :), s0%frame, s0%center)
end function save_trajectory

subroutine setfields(tra, fields)
    type(trajectory), intent(inout) :: tra
    character(len=*), dimension(:), intent(in) :: fields

    integer :: i
    integer :: n

    n = size(fields)
    if (allocated(tra%fields)) deallocate(tra%fields)
    allocate(tra%fields(n))
    do i = 1, n
        if (len(fields(i)) < fieldlen) then
            tra%fields(i) = fields(i) // repeat(" ", fieldlen - len(fields(i)))
        else
            tra%fields(i) = fields(i)
        end if
    end do
end subroutine setfields

function length(nd) result(len)
    type(tranode), intent(in) :: nd
    integer :: len

    type(tranode), pointer :: current

    current => nd%next
    len = 1
    do while (associated(current))
        len = len + 1
        current => current%next
    end do
end function length

subroutine add_node(nd, t, y, head)
    type(tranode), intent(in), target :: nd
    real(dp), intent(in) :: t
    real(dp), dimension(:), intent(in) :: y
    type(tranode), intent(inout), pointer, optional :: head

    type(tranode), pointer :: current

    current => nd
    do while (associated(current%next))
        current => current%next
    end do

    allocate(current%next)
    current%next%t = t
    current%next%y = y
    if (present(head)) head => current%next
end subroutine add_node

end module trajectories
