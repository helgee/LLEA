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
!   The MIT License (MIT)
!
!   Copyright (c) 2012-2015 Juan Luis Cano Rodríguez
!
!   Permission is hereby granted, free of charge, to any person obtaining a copy
!   of this software and associated documentation files (the "Software"), to deal
!   in the Software without restriction, including without limitation the rights
!   to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
!   copies of the Software, and to permit persons to whom the Software is
!   furnished to do so, subject to the following conditions:
!
!   The above copyright notice and this permission notice shall be included in all
!   copies or substantial portions of the Software.
!
!   THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
!   IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
!   FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
!   AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
!   LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
!   OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
!   SOFTWARE.
!
module states

use types, only: dp
use bodies, only: body, iaumatrix
use constants, only: earth
use epochs, only: epoch
use exceptions
use math, only: pi, cross, isapprox, mod2pi, norm
use rotations, only: rotationmatrix

implicit none

integer, parameter :: framelen = 8

type state
    type(epoch) :: ep
    real(dp), dimension(6) :: rv = 0._dp
    character(len=framelen) :: frame = "GCRF"
    type(body) :: center
end type state

interface keplerian
    module procedure keplerian_vectors
    module procedure keplerian_state
end interface keplerian

interface state
    module procedure state_init
end interface state

private

public :: state, rotate_inplace, rotate, cartesian, keplerian, state_from_elements, framelen

contains

function state_init(ep, rv, frame, center) result(s)
    type(epoch), intent(in) :: ep
    real(dp), dimension(:), intent(in) :: rv
    character(len=framelen), intent(in), optional :: frame
    type(body), intent(in), optional :: center
    type(state) :: s

    s%ep = ep
    s%rv = rv
    if (present(frame)) s%frame = frame
    s%center = earth
    if (present(center)) s%center = center
end function state_init

subroutine rotate_inplace(s, to, err)
    type(state), intent(inout) :: s
    character(len=*), intent(in) :: to
    type(exception), intent(inout), optional :: err

    type(exception) :: err_
    real(dp), dimension(3,3) :: m
    real(dp), dimension(3,3) :: dm
    real(dp), dimension(6,6) :: rot

    rot = 0._dp

    select case (s%frame)
    case ("GCRF")
        select case (to)
        case ("IAU")
            call iaumatrix(s%center, s%ep, m, dm)
            rot(1:3,1:3) = m
            rot(4:6,4:6) = m
            rot(4:6,1:3) = dm
            s%rv = matmul(rot, s%rv)
        end select
    case ("IAU")
        select case (to)
        case ("GCRF")
            call iaumatrix(s%center, s%ep, m, dm)
            m = transpose(m)
            dm = transpose(dm)
            rot(1:3,1:3) = m
            rot(4:6,4:6) = m
            rot(4:6,1:3) = dm
            s%rv = matmul(rot, s%rv)
        end select
    case default
        err_ = error("Unknown target frame: "//to, "rotate_inplace", __FILE__, __LINE__)
        if (present(err)) then
            err = err_
            return
        else
            call raise(err_)
        end if
    end select
    s%frame = to
end subroutine rotate_inplace

function rotate(s, to, err) result(s1)
    type(state), intent(in) :: s
    character(len=*), intent(in) :: to
    type(exception), intent(inout), optional :: err

    type(exception) :: err_
    type(state) :: s1

    s1 = s
    call rotate_inplace(s1, to, err_)
    if (iserror(err_)) then
        call catch(err_, "rotate", __FILE__, __LINE__)
        if (present(err)) then
            err = err_
        else
            call raise(err_)
        end if
    end if
end function rotate

pure function keplerian_vectors(rv, mu) result(ele)
    real(dp), dimension(:), intent(in) :: rv
    real(dp), intent(in) :: mu
    real(dp), dimension(6) :: ele

    real(dp) :: rm, vm2, hm, nm, xi
    real(dp), dimension(3) :: h, n, e
    logical :: equatorial
    logical :: circular

    rm = norm(rv(1:3))
    vm2 = sum(rv(4:6)**2)
    h = cross(rv(1:3),rv(4:6))
    hm = norm(h)
    n = cross([0._dp, 0._dp, 1._dp], h)
    nm = norm(n)
    xi = vm2 / 2._dp - mu / rm
    e = ((vm2 - mu / rm) * rv(1:3) - rv(4:6) &
        * dot_product(rv(1:3), rv(4:6))) / mu
    ele(2) = norm(e)
    ele(3) = acos(h(3) / hm)

    equatorial = isapprox(abs(ele(3)), 0._dp)
    circular = isapprox(ele(2), 0._dp)

    if (circular) then
        ele(1) = hm**2 / mu
    else
        ele(1) = -mu / (2._dp * xi)
    end if

    if (equatorial .and. .not.circular) then
        ele(4) = 0._dp
        ele(5) = mod2pi(atan2(e(2), e(1)))
        ele(6) = mod2pi(atan2(dot_product(h, cross(e, rv(1:3))) / hm, &
            dot_product(rv(1:3), e)))
    else if (.not.equatorial .and. circular) then
        ele(4) = mod2pi(atan2(n(2), n(1)))
        ele(5) = 0._dp
        ele(6) = mod2pi(atan2(dot_product(rv(1:3), cross(h, n)) / hm, &
            dot_product(rv(1:3), n)))
    else if (equatorial .and. circular) then
        ele(4) = 0._dp
        ele(5) = 0._dp
        ele(6) = mod2pi(atan2(rv(2), rv(1)))
    else
        ele(4) = mod2pi(atan2(n(2), n(1)))
        ele(5) = mod2pi(atan2(dot_product(e, cross(h, n)) / hm, &
            dot_product(e, n)))
        ele(6) = mod2pi(atan2(dot_product(rv(1:3), cross(h, e)) / hm, &
            dot_product(rv(1:3), e)))
    end if
end function keplerian_vectors

pure function keplerian_state(s) result(ele)
    type(state), intent(in) :: s
    double precision, dimension(6) :: ele

    ele = keplerian_vectors(s%rv, s%center%mu)
end function keplerian_state

pure function perifocal(p, ecc, ano, mu) result(rvp)
    real(dp), intent(in) :: p
    real(dp), intent(in) :: ecc
    real(dp), intent(in) :: ano
    real(dp), intent(in) :: mu

    real(dp), dimension(6) :: rvp

    rvp(1:3) = [p * cos(ano) / (1 + ecc * cos(ano)), p * sin(ano) / (1 + ecc * cos(ano)), 0._dp]
    rvp(4:6) = [-sqrt(mu / p) * sin(ano), sqrt(mu / p) * (ecc + cos(ano)), 0._dp]
end function perifocal

function cartesian(ele, mu) result(rv)
    real(dp), dimension(:), intent(in) :: ele
    real(dp), intent(in) :: mu

    real(dp), dimension(6) :: rv
    real(dp), dimension(6) :: rvp
    real(dp), dimension(6,6) :: m
    real(dp), dimension(3,3) :: m1
    real(dp) :: p

    m = 0._dp

    if (isapprox(ele(2), 0._dp)) then
        p = ele(1)
    else
        p = ele(1) * (1 - ele(2)**2)
    end if

    rvp = perifocal(p, ele(2), ele(6), mu)
    m1 = rotationmatrix("313", -ele(5), -ele(3), -ele(4))
    m(1:3,1:3) = m1
    m(4:6,4:6) = m1

    rv = matmul(m, rvp)
end function cartesian

function state_from_elements(ep, ele, frame, center) result(s)
    type(epoch), intent(in) :: ep
    real(dp), dimension(6), intent(in) :: ele
    character(len=*), intent(in), optional :: frame
    type(body), intent(in), optional :: center

    type(state) :: s
    character(len=framelen) :: frame_
    type(body) :: center_
    real(dp), dimension(6) :: rv

    frame_ = "GCRF"
    if (present(frame)) frame_ = frame
    center_ = earth
    if (present(center)) center_ = center

    rv = cartesian(ele, center%mu)
    s%ep = ep
    s%rv = rv
    if (present(frame)) s%frame = frame
    s%center = earth
    if (present(center)) s%center = center
end function state_from_elements

end module states
