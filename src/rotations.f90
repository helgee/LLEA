!
! Copyright (c) 2016 Helge Eichhorn
!
! This Source Code Form is subject to the terms of the Mozilla Public
! License, v. 2.0. If a copy of the MPL was not distributed with this
! file, You can obtain one at http://mozilla.org/MPL/2.0/.
!
module rotations

use types, only: dp
use exceptions

implicit none

interface rotationmatrix
    module procedure rotationmatrix_single
    module procedure rotationmatrix_euler
end interface rotationmatrix

interface ratematrix
    module procedure ratematrix_single
    module procedure ratematrix_euler
end interface ratematrix

contains

function checkaxes(axes) result(err)
    character(len=3), intent(in) :: axes
    type(exception) :: err

    integer :: i

    do i = 1, 3
        select case (axes(i:i))
        case ("1", "2", "3", "x", "y", "z", "X", "Y", "Z")
            continue
        case default
            err = error("Invalid rotation axis: '"//axes(i:i)//"'.", "checkaxes", &
                __FILE__, __LINE__)
            return
        end select
    end do

    if ((axes(1:1) == axes(2:2)).or.(axes(2:2) == axes(3:3))) then
        err = error("Subsequent rotations around the same axis are meaningless.", "checkaxes", &
            __FILE__, __LINE__)
        return
    end if
end function checkaxes

function rotationmatrix_euler(axes, angle1, angle2, angle3, err) result(m)
    character(len=3), intent(in) :: axes
    real(dp), intent(in) :: angle1
    real(dp), intent(in) :: angle2
    real(dp), intent(in) :: angle3
    type(exception), intent(inout), optional :: err

    type(exception) :: err_
    real(dp), dimension(3,3) :: m
    real(dp), dimension(3,3) :: m1
    real(dp), dimension(3,3) :: m2
    real(dp), dimension(3,3) :: m3

    err_ = checkaxes(axes)
    if (iserror(err_)) then
        if (present(err)) then
            err = err_
        else
            call raise(err_)
        end if
    end if

    m1 = rotationmatrix_single(axes(1:1), angle1)
    m2 = rotationmatrix_single(axes(2:2), angle2)
    m3 = rotationmatrix_single(axes(3:3), angle3)
    m = matmul(m3, matmul(m2, m1))
end function rotationmatrix_euler

function ratematrix_euler(axes, angle1, rate1, angle2, rate2, angle3, rate3, err) result(m)
    character(len=3), intent(in) :: axes
    real(dp), intent(in) :: angle1
    real(dp), intent(in) :: rate1
    real(dp), intent(in) :: angle2
    real(dp), intent(in) :: rate2
    real(dp), intent(in) :: angle3
    real(dp), intent(in) :: rate3
    type(exception), intent(inout), optional :: err

    type(exception) :: err_
    real(dp), dimension(3,3) :: m
    real(dp), dimension(3,3) :: ma1
    real(dp), dimension(3,3) :: ma2
    real(dp), dimension(3,3) :: ma3
    real(dp), dimension(3,3) :: mr1
    real(dp), dimension(3,3) :: mr2
    real(dp), dimension(3,3) :: mr3

    err_ = checkaxes(axes)
    if (iserror(err_)) then
        if (present(err)) then
            err = err_
        else
            call raise(err_)
        end if
    end if

    ma1 = rotationmatrix_single(axes(1:1), angle1)
    ma2 = rotationmatrix_single(axes(2:2), angle2)
    ma3 = rotationmatrix_single(axes(3:3), angle3)
    mr1 = ratematrix_single(axes(1:1), angle1, rate1)
    mr2 = ratematrix_single(axes(2:2), angle2, rate2)
    mr3 = ratematrix_single(axes(3:3), angle3, rate3)
    m = matmul(mr3, matmul(ma2, ma1)) + matmul(ma3, matmul(mr2, ma1)) + matmul(ma3, matmul(ma2, mr1))
end function ratematrix_euler

! Function: rot
!   Generates a rotation matrix.
!
! Parameters:
!   ang - Rotation angle in radians.
!   axis - Rotation axis (1,2,3).
!
! Returns:
!   m - 3x3 Rotation matrix
function rotationmatrix_single(axis, angle, err) result(m)
    character(len=1), intent(in) :: axis
    real(dp), intent(in) :: angle
    type(exception), intent(inout), optional :: err
    real(dp), dimension(3,3) :: m

    real(dp) :: sina
    real(dp) :: cosa
    type(exception) :: err_

    m = 0._dp
    sina = sin(angle)
    cosa = cos(angle)

    select case (axis)
    case ("1", "x", "X")
        m(1,1) = 1._dp
        m(2,2) = cosa
        m(2,3) = sina
        m(3,2) = -sina
        m(3,3) = cosa
    case ("2", "y", "Y")
        m(1,1) = cosa
        m(1,3) = -sina
        m(2,2) = 1._dp
        m(3,1) = sina
        m(3,3) = cosa
    case("3", "z", "Z")
        m(1,1) = cosa
        m(1,2) = sina
        m(2,1) = -sina
        m(2,2) = cosa
        m(3,3) = 1._dp
    case default
        err_ = error("Invalid rotation axis: '"//axis//"'.", "rotationmatrix_single", &
            __FILE__, __LINE__)
        if (present(err)) then
            err = err_
            return
        else
            call raise(err_)
        end if
    end select
end function rotationmatrix_single

! Function: rotd
!   Generates a derivative rotation matrix.
!
! Parameters:
!   ang - Rotation angle in radians.
!   dang - Rotational velocity in radians/sec.
!   axis - Rotation axis (1,2,3).
!
! Returns:
!   m - 3x3 Derivative rotation matrix
function ratematrix_single(axis, angle, rate, err) result(m)
    character(len=1), intent(in) :: axis
    real(dp), intent(in) :: angle
    real(dp), intent(in) :: rate
    type(exception), intent(inout), optional :: err
    real(dp), dimension(3,3) :: m

    real(dp) :: sina
    real(dp) :: cosa
    type(exception) :: err_

    m = 0._dp
    sina = sin(angle)
    cosa = cos(angle)

    select case (axis)
    case ("1", "x", "X")
        m(2,2) = -rate * sina
        m(2,3) = rate * cosa
        m(3,2) = -rate * cosa
        m(3,3) = -rate * sina
    case ("2", "y", "Y")
        m(1,1) = -rate * sina
        m(1,3) = rate * cosa
        m(3,1) = -rate * cosa
        m(3,3) = -rate * sina
    case("3", "z", "Z")
        m(1,1) = -rate * sina
        m(1,2) = rate * cosa
        m(2,1) = -rate * cosa
        m(2,2) = -rate * sina
    case default
        err_ = error("Invalid rotation axis: '"//axis//"'.", "ratematrix_single", &
            __FILE__, __LINE__)
        if (present(err)) then
            err = err_
            return
        else
            call raise(err_)
        end if
    end select
end function ratematrix_single

end module rotations
