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

use epochs, only: epochdelta, operator (+), seconds
use exceptions
use splines, only: spline1d, evaluate
use states, only: state
use types, only: dp

implicit none

private

public :: tranode, trajectory, len_dirty, add_node, isdirty, save_trajectory, getfield, &
    state_trajectory_epd, state_trajectory_dp, interpolate_scalar

integer, parameter :: fieldlen = 12

type tranode
    type(tranode), pointer :: next => null()
    real(dp) :: t
    real(dp), dimension(:), allocatable :: y
end type tranode

type trajectory
    type(state) :: initial_state
    type(state) :: final_state
    type(tranode), pointer :: head => null()
    type(tranode), pointer :: tail => null()
    character(len=fieldlen), dimension(:), allocatable :: fields
    real(dp), dimension(:), allocatable :: t
    real(dp), dimension(:,:), allocatable :: y
    type(spline1d), dimension(:), allocatable :: spl
end type trajectory

interface tranode
    module procedure init_tranode
end interface tranode

interface trajectory
    module procedure init_trajectory_array
    module procedure init_trajectory
end interface trajectory

contains

function init_trajectory(s0, fields) result(tra)
    type(state), intent(in) :: s0
    character(len=*), dimension(:), intent(in), optional :: fields
    type(trajectory) :: tra

    tra%initial_state = s0

    call setfields(tra, ["x ", "y ", "z ", "vx", "vy", "vz"])
    if (present(fields)) then
        if (allocated(tra%fields)) deallocate(tra%fields)
        call setfields(tra, fields)
    end if
end function init_trajectory

function init_tranode(t, y) result(node)
    real(dp), intent(in) :: t
    real(dp), dimension(:), intent(in) :: y
    type(tranode) :: node
    node%t = t
    node%y = y
end function init_tranode

subroutine getfield(tra, field, res, err)
    type(trajectory), intent(in) :: tra
    character(len=*), intent(in) :: field
    real(dp), dimension(:), intent(out), allocatable :: res
    type(exception), intent(inout), optional :: err

    type(exception) :: err_
    integer :: ind
    integer :: i

    if (.not.allocated(tra%y)) then
        err_ = error("Trajectory is empty.", "getfield", __FILE__, __LINE__)
        if (present(err)) then
            err = err_
            return
        else
            call raise(err_)
        end if
    end if

    ind = 0
    do i = 1, size(tra%fields)
        if (tra%fields(i) == field) then
            ind = i
            exit
        end if
    end do
    if (ind == 0) then
        err_ = error("Unknown field: "//trim(field), "getfield", __FILE__, __LINE__)
        if (present(err)) then
            err = err_
            return
        else
            call raise(err_)
        end if
    end if

    res = tra%y(ind,:)
end subroutine getfield

function init_trajectory_array(s0, t, arr, fields) result(tra)
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
    tra%y = arr
    tra%t = t
    tra%final_state = state(s0%ep + epochdelta(seconds=t(n)), arr(:,n), s0%frame, s0%center)
    call generate_splines(tra)
end function init_trajectory_array

subroutine save_trajectory(tra)
    type(trajectory), intent(inout) :: tra

    integer :: n
    integer :: m
    integer :: i
    type(tranode), pointer :: current

    if (.not.isdirty(tra)) return

    n = len_dirty(tra)
    m = size(tra%head%y)

    allocate(tra%t(n))
    allocate(tra%y(m,n))
    tra%t(1) = tra%head%t
    tra%y(:,1) = tra%head%y
    i = 1

    current => tra%head%next
    do while (associated(current))
        i = i + 1
        tra%t(i) = current%t
        tra%y(:,i) = current%y
        current => current%next
    end do
    call delete_node(tra%head)
    tra%final_state = state(tra%initial_state%ep + epochdelta(seconds=tra%t(n)), tra%y(:,n), &
        tra%initial_state%frame, tra%initial_state%center)
    call generate_splines(tra)
end subroutine save_trajectory

subroutine generate_splines(tra)
    type(trajectory), intent(inout) :: tra

    integer :: n
    integer :: i

    n = size(tra%fields)
    allocate(tra%spl(n))

    do i = 1, n
        tra%spl(i) = spline1d(tra%t, tra%y(i,:), bc="extrapolate")
    end do
end subroutine generate_splines

function interpolate_scalar(tra, t, n, extrapolate, err) result(y)
    type(trajectory), intent(in) :: tra
    real(dp), intent(in) :: t
    integer, intent(in), optional :: n
    logical, intent(in), optional :: extrapolate
    type(exception), intent(inout), optional :: err
    real(dp), dimension(:), allocatable :: y

    integer :: n_
    integer :: i
    logical :: extrapolate_
    type(exception) :: err_

    n_ = size(tra%fields)
    if (present(n)) then
        if (n < n_) n_ = n
    end if
    allocate(y(n_))

    extrapolate_ = .false.
    if (present(extrapolate)) extrapolate_ = extrapolate

    if ((t < tra%t(1).or.t > tra%t(size(tra%t))).and..not.extrapolate_) then
        err_ = error("t is outside range.", "interpolate_scalar", __FILE__, __LINE__)
        if (present(err)) then
            err = err_
            return
        else
            call raise(err_)
        end if
    end if

    do i = 1, n_
        y(i) = evaluate(tra%spl(i), t)
    end do
end function interpolate_scalar

function state_trajectory_dp(tra, t, extrapolate, err) result(s)
    type(trajectory), intent(in) :: tra
    real(dp), intent(in) :: t
    logical, intent(in), optional :: extrapolate
    type(exception), intent(inout), optional :: err
    type(state) :: s

    logical :: extrapolate_
    type(exception) :: err_
    real(dp), dimension(:), allocatable :: y

    extrapolate_ = .false.
    if (present(extrapolate)) extrapolate_ = extrapolate

    y = interpolate_scalar(tra, t, extrapolate=extrapolate_, err=err_)
    if (iserror(err_)) then
        call catch(err_, "state_trajectory", __FILE__, __LINE__)
        if (present(err)) then
            err = err_
            return
        else
            call raise(err_)
        end if
    end if
    s = state(tra%initial_state%ep + epochdelta(seconds=t), y, tra%initial_state%frame, tra%initial_state%center)
end function state_trajectory_dp

function state_trajectory_epd(tra, epd, extrapolate, err) result(s)
    type(trajectory), intent(in) :: tra
    type(epochdelta), intent(in) :: epd
    logical, intent(in), optional :: extrapolate
    type(exception), intent(inout), optional :: err
    type(state) :: s

    logical :: extrapolate_
    type(exception) :: err_
    real(dp), dimension(:), allocatable :: y
    real(dp) :: t

    extrapolate_ = .false.
    if (present(extrapolate)) extrapolate_ = extrapolate

    t = seconds(epd)
    y = interpolate_scalar(tra, t, extrapolate=extrapolate_, err=err_)
    if (iserror(err_)) then
        call catch(err_, "state_trajectory", __FILE__, __LINE__)
        if (present(err)) then
            err = err_
            return
        else
            call raise(err_)
        end if
    end if
    s = state(tra%initial_state%ep + epd, y, tra%initial_state%frame, tra%initial_state%center)
end function state_trajectory_epd

recursive subroutine delete_node(node)
    type(tranode), pointer :: node

    if (associated(node)) then
        call delete_node(node%next)
        deallocate(node)
        nullify(node)
    end if
end subroutine delete_node

function isdirty(tra) result(res)
    type(trajectory), intent(in) :: tra
    logical :: res
    res = associated(tra%head)
end function isdirty

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

function len_dirty(tra) result(len)
    type(trajectory), intent(in) :: tra
    integer :: len

    type(tranode), pointer :: current

    if (.not.isdirty(tra)) then
        len = 0
        return
    else
        current => tra%head%next
        len = 1
        do while (associated(current))
            len = len + 1
            current => current%next
        end do
    end if
end function len_dirty

subroutine add_node(tra, t, y, err)
    type(trajectory), intent(inout), target :: tra
    real(dp), intent(in) :: t
    real(dp), dimension(:), intent(in) :: y
    type(exception), intent(inout), optional :: err

    type(exception) :: err_

    if (size(y) /= size(tra%fields)) then
        err_ = error("Wrong number of fields.", "add_node", __FILE__, __LINE__)
        if (present(err)) then
            err = err_
            return
        else
            call raise(err_)
        end if
    end if

    if (.not.associated(tra%head)) then
        allocate(tra%head, source=tranode(t, y))
        tra%tail => tra%head
    else
        allocate(tra%tail%next, source=tranode(t, y))
        tra%tail => tra%tail%next
    end if
end subroutine add_node

end module trajectories
