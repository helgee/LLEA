!
! Copyright (c) 2016 Helge Eichhorn
!
! This Source Code Form is subject to the terms of the Mozilla Public
! License, v. 2.0. If a copy of the MPL was not distributed with this
! file, You can obtain one at http://mozilla.org/MPL/2.0/.
!
module forces

use bodies, only: body
use containers, only: parameters
use ephemerides, only: ephem, getposition
use epochs, only: epoch, epochdelta, operator(+), epochdelta_
use math, only: norm
use states, only: state
use types, only: dp

implicit none

type, abstract :: model
contains
    procedure(eval), deferred :: eval
    procedure(update), deferred :: update
end type model

interface
    function eval(this, t, y, p) result(f)
        import :: model, dp, parameters
        class(model), intent(in) :: this
        real(dp), intent(in) :: t
        real(dp), dimension(:), intent(in) :: y
        type(parameters), intent(in) :: p
        real(dp), dimension(:), allocatable :: f
    end function eval
    subroutine update(this, f, t, y, p)
        import :: model, dp, parameters
        class(model), intent(in) :: this
        real(dp), dimension(:), intent(inout) :: f
        real(dp), intent(in) :: t
        real(dp), dimension(:), intent(in) :: y
        type(parameters), intent(in) :: p
    end subroutine update
end interface

type, extends(model), abstract :: gravity
end type gravity

type, extends(model), abstract :: drag
end type drag

type, extends(gravity) :: uniformgravity
contains
    procedure :: eval => uniformgravity_eval
    procedure :: update => uniformgravity_update
end type uniformgravity

type, extends(gravity) :: j2gravity
contains
    procedure :: eval => j2gravity_eval
    procedure :: update => j2gravity_update
end type j2gravity

type, extends(model) :: thirdbody
    type(body), dimension(:), allocatable :: bodies
contains
    procedure :: eval => thirdbody_eval
    procedure :: update => thirdbody_update
end type thirdbody

contains

function uniformgravity_eval(this, t, y, p) result(f)
    class(uniformgravity), intent(in) :: this
    real(dp), intent(in) :: t
    real(dp), dimension(:), intent(in) :: y
    type(parameters), intent(in) :: p
    real(dp), dimension(:), allocatable :: f

    allocate(f(6))
    call this%update(f, t, y, p)
end function uniformgravity_eval

subroutine uniformgravity_update(this, f, t, y, p)
    class(uniformgravity), intent(in) :: this
    real(dp), dimension(:), intent(inout) :: f
    real(dp), intent(in) :: t
    real(dp), dimension(:), intent(in) :: y
    type(parameters), intent(in) :: p

    real(dp) :: r
    real(dp) :: r3

    r = norm(y(1:3))
    r3 = r * r * r
    f(1:3) = f(1:3) + y(4:6)

    f(4:6) = f(4:6) - p%center%mu * y(1:3) / r3
end subroutine uniformgravity_update

function j2gravity_eval(this, t, y, p) result(f)
    class(j2gravity), intent(in) :: this
    real(dp), intent(in) :: t
    real(dp), dimension(:), intent(in) :: y
    type(parameters), intent(in) :: p
    real(dp), dimension(:), allocatable :: f

    allocate(f(6))
    call this%update(f, t, y, p)
end function j2gravity_eval

subroutine j2gravity_update(this, f, t, y, p)
    class(j2gravity), intent(in) :: this
    real(dp), dimension(:), intent(inout) :: f
    real(dp), intent(in) :: t
    real(dp), dimension(:), intent(in) :: y
    type(parameters), intent(in) :: p

    real(dp) :: r
    real(dp) :: r2
    real(dp) :: r3
    real(dp) :: pj

    r = norm(y(1:3))
    r2 = r * r
    r3 = r * r * r
    f(1:3) = f(1:3) + y(4:6)

    pj = -3._dp / 2._dp * p%center%mu * p%center%j2 * p%center%radii(1)**2 / (r2*r3)
    f(4:5) = f(4:5) - p%center%mu * y(1:2) / r3 + pj * y(1:2) * (1._dp - 5._dp * y(3)*y(3) / r2)
    f(6) = f(6) - p%center%mu * y(3) / r3 + pj * y(3) * (3._dp - 5._dp * y(3)*y(3) / r2)
end subroutine j2gravity_update

function thirdbody_eval(this, t, y, p) result(f)
    class(thirdbody), intent(in) :: this
    real(dp), intent(in) :: t
    real(dp), dimension(:), intent(in) :: y
    type(parameters), intent(in) :: p
    real(dp), dimension(:), allocatable :: f

    allocate(f(6))
    call this%update(f, t, y, p)
end function thirdbody_eval

subroutine thirdbody_update(this, f, t, y, p)
    class(thirdbody), intent(in) :: this
    real(dp), dimension(:), intent(inout) :: f
    real(dp), intent(in) :: t
    real(dp), dimension(:), intent(in) :: y
    type(parameters), intent(in) :: p

    integer :: i
    real(dp) :: mu
    real(dp), dimension(3) :: rc3
    real(dp), dimension(3) :: rs3
    type(epoch) :: ep

    ep = p%s0%ep + epochdelta_(seconds=t)
    if (allocated(this%bodies)) then
        do i = 1, size(this%bodies)
            mu = this%bodies(i)%mu
            rc3 = getposition(ephem, ep, from=p%center%id, to=this%bodies(i)%id)
            rs3 = rc3 - y(1:3)
            f(4:6) = f(4:6) + mu * (rs3 / norm(rs3)**3 - rc3 / norm(rc3)**3)
        end do
    end if
end subroutine thirdbody_update

end module forces
