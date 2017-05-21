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
!   Copyright (c) 2001, 2002 Enthought, Inc.
!   All rights reserved.
!
!   Copyright (c) 2003-2012 SciPy Developers.
!   All rights reserved.
!
!   Copyright (c) 2014 Kyle Barbary.
!   All rights reserved.
!
!   Redistribution and use in source and binary forms, with or without
!   modification, are permitted provided that the following conditions are met:
!
!     a. Redistributions of source code must retain the above copyright notice,
!        this list of conditions and the following disclaimer.
!     b. Redistributions in binary form must reproduce the above copyright
!        notice, this list of conditions and the following disclaimer in the
!        documentation and/or other materials provided with the distribution.
!     c. Neither the name of Enthought nor the names of the SciPy Developers
!        may be used to endorse or promote products derived from this software
!        without specific prior written permission.
!
!
!   THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
!   AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
!   IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
!   ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDERS OR CONTRIBUTORS
!   BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY,
!   OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF
!   SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS
!   INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN
!   CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE)
!   ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF
!   THE POSSIBILITY OF SUCH DAMAGE.
!
module splines

use exceptions
use types, only: dp
use util, only: resize

implicit none

private

public :: spline1d, evaluate, knots, coeffs, residual, spline1d_

type spline1d
    real(dp), dimension(:), allocatable :: t
    real(dp), dimension(:), allocatable :: c
    integer :: k
    integer :: bc
    real(dp) :: fp
end type spline1d

interface spline1d_
    module procedure init_spline1d
end interface spline1d_

interface evaluate
    module procedure eval_1d_scalar
    module procedure eval_1d_vector
end interface evaluate

interface knots
    module procedure knots_1d
end interface knots

interface coeffs
    module procedure coeffs_1d
end interface coeffs

interface residual
    module procedure residual_1d
end interface residual

interface
    subroutine splev(t, n, c, k, x, y, m, e, ier)
        import :: dp
        integer, intent(in) :: n
        integer, intent(in) :: m
        real(dp), dimension(n), intent(in) :: t
        real(dp), dimension(n), intent(in) :: c
        integer, intent(in) :: k
        real(dp), dimension(m), intent(in) :: x
        real(dp), dimension(m), intent(out) :: y
        integer, intent(in) :: e
        integer, intent(out) :: ier
    end subroutine splev
end interface

contains

function init_spline1d(x, y, w, k, s, bc, err) result(spl)
    real(dp), dimension(:), intent(in) :: x
    real(dp), dimension(:), intent(in) :: y
    real(dp), dimension(:), intent(in), optional :: w
    integer, intent(in), optional :: k
    real(dp), intent(in), optional :: s
    character(len=*), intent(in), optional :: bc
    type(exception), intent(inout), optional :: err
    type(spline1d) :: spl

    character(len=11) :: bc_
    integer :: bccode
    integer :: ier
    integer :: k_
    integer :: lwrk
    integer :: m
    integer :: n
    integer :: nest
    real(dp) :: fp
    real(dp) :: s_
    real(dp), dimension(:), allocatable :: c
    real(dp), dimension(:), allocatable :: t
    real(dp), dimension(:), allocatable :: w_
    type(exception) :: err_
    real(dp), dimension(:), allocatable :: wrk
    integer, dimension(:), allocatable :: iwrk

    m = size(x)
    if (present(w)) then
        w_ = w
    else
        allocate(w_(m))
        w_ = 1._dp
    end if
    k_ = 3
    if (present(k)) k_ = k
    s_ = 0._dp
    if (present(s)) s_ = s
    bc_ = "nearest"
    if (present(bc)) bc_ = bc

    if (size(y) /= m) err_ = error("Length of x and y must match.", "init_spline1d", __FILE__, __LINE__)
    if (present(w)) then
        if (size(w) /= m) err_ = error("Length of x and w must match.", "init_spline1d", __FILE__, __LINE__)
    end if
    if (m <= k_) err_ = error("Length of x must be greater than k.", "init_spline1d", __FILE__, __LINE__)
    if (iserror(err_)) then
        if (present(err)) then
            err = err_
            return
        else
            call raise(err_)
        end if
    end if

    nest = m + k_ + 1
    allocate(t(nest))
    allocate(c(nest))

    lwrk = m * (k_ + 1) + nest * (7 + 3 * k_)
    allocate(wrk(lwrk))
    allocate(iwrk(nest))

    call curfit(0, m, x, y, w_, x(1), x(m), k_, s_, nest, n, t, c, fp, wrk, lwrk, iwrk, ier)
    if (ier > 0) then
        err_ = error("s is too small or other invalid input.", "init_spline1d", __FILE__, __LINE__)
        if (present(err)) then
            err = err_
            return
        else
            call raise(err_)
        end if
    end if
    ! Resize output arrays
    call resize(t, n)
    call resize(c, n-k_-1)
    bccode = translate_bc(bc_, err_)
    spl = spline1d(t, c, k_, bccode, fp)
    if (iserror(err_)) then
        call catch(err_, "init_spline1d", __FILE__, __LINE__)
        if (present(err)) then
            err = err_
            return
        else
            call raise(err_)
        end if
    end if
end function init_spline1d

function knots_1d(spl) result(res)
    type(spline1d), intent(in) :: spl
    real(dp), dimension(:), allocatable :: res

    res = spl%t(spl%k+1:size(spl%t)-spl%k)
end function knots_1d

function coeffs_1d(spl) result(res)
    type(spline1d), intent(in) :: spl
    real(dp), dimension(:), allocatable :: res

    res = spl%c(:size(spl%c)-spl%k+1)
end function coeffs_1d

function residual_1d(spl) result(res)
    type(spline1d), intent(in) :: spl
    real(dp) :: res

    res = spl%fp
end function residual_1d

function eval_1d_scalar(spl, x, err) result(y)
    type(spline1d), intent(in) :: spl
    real(dp), intent(in) :: x
    type(exception), intent(inout), optional :: err
    real(dp) :: y

    type(exception) :: err_
    real(dp), dimension(1) :: x_
    real(dp), dimension(1) :: y_
    integer :: ier

    y = 0._dp
    x_ = x
    call splev(spl%t, size(spl%t), spl%c, spl%k, x_, y_, 1, spl%bc, ier)
    if (ier /= 0) then
        err_ = error(trim(eval_1d_message(ier)), "eval_1d_scalar", __FILE__, __LINE__)
        if (present(err)) then
            err = err_
            return
        else
            call raise(err_)
        end if
    end if
    y = y_(1)
end function eval_1d_scalar

function eval_1d_vector(spl, x, err) result(y)
    type(spline1d), intent(in) :: spl
    real(dp), dimension(:), intent(in) :: x
    type(exception), intent(inout), optional :: err
    real(dp), dimension(:), allocatable :: y

    type(exception) :: err_
    integer :: ier
    integer :: m

    m = size(x)
    allocate(y(m))

    call splev(spl%t, size(spl%t), spl%c, spl%k, x, y, m, spl%bc, ier)
    if (ier /= 0) then
        err_ = error(trim(eval_1d_message(ier)), "eval_1d_vector", __FILE__, __LINE__)
        if (present(err)) then
            err = err_
            return
        else
            call raise(err_)
        end if
    end if
end function eval_1d_vector

function translate_bc(bc, err) result(res)
    character(len=*), intent(in) :: bc
    type(exception), intent(inout) :: err
    integer :: res

    res = -1
    select case (bc)
    case ("extrapolate")
        res = 0
    case ("zero")
        res = 1
    case ("error")
        res = 2
    case ("nearest")
        res = 3
    case default
        err = error("Unknown boundary condition: "//bc, "translate_bc", __FILE__, __LINE__)
    end select
end function translate_bc

function eval_1d_message(ier) result(msg)
    integer, intent(in) :: ier
    character(len=25) :: msg

    select case (ier)
    case (10)
        msg = "Invalid input data."
    case (1)
        msg = "Input point out of range."
    end select
end function eval_1d_message

end module splines
