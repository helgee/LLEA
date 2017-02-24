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
!   Copyright (c) 2003-2016 SciPy Developers.
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
module math

use, intrinsic :: iso_c_binding, only: c_ptr, c_null_ptr
use types, only: dp
use exceptions

implicit none

private

public :: pih, pi, twopi, deg, rad, &
    eps, cross, unitmat, mod2pi, &
    vecang, polcart, cartpol, norm, &
    getsign, binsearch, linspace, &
    interp, findroot, isapprox, &
    issorted, greatcircle, isin, &
    deg2rad, rad2deg, cot

! Constants: Mathematical constants
!
!   eps - Machine epsilon for double precision.
!   pih - pi/2
!   pi - pi
!   twopi - 2*pi
!   deg - Conversion factor from radians to degrees: 180/pi
!   rad - Conversion factor from degrees to radians: pi/180
real(dp), parameter :: eps = epsilon(0.0_dp)
real(dp), parameter :: pih = 3.1415926535897931_dp / 2
real(dp), parameter :: pi = 3.1415926535897931_dp
real(dp), parameter :: twopi = 2*pi
real(dp), parameter :: deg = 180._dp / pi
real(dp), parameter :: rad = pi / 180._dp

interface isapprox
    module procedure isapprox_scalar
    module procedure isapprox_vector
    module procedure isapprox_matrix
end interface isapprox

interface norm
    module procedure norm_vector
    module procedure norm_matrix
end interface norm

interface isin
    module procedure isin_vector_dp
    module procedure isin_matrix_dp
    module procedure isin_vector_int
    module procedure isin_matrix_int
end interface isin

interface findroot
    module procedure findroot_brent
    module procedure findroot_steffensen
end interface findroot

abstract interface
    function func(x, rpar, ipar) result(res)
        import :: dp
        real(dp), intent(in) :: x
        real(dp), dimension(:), intent(in) :: rpar
        integer, dimension(:), intent(in) :: ipar
        real(dp) :: res
    end function func
end interface

contains

elemental function cot(arg) result(res)
    real(dp), intent(in) :: arg
    real(dp) :: res

    res = 1._dp / tan(arg)
end function cot

elemental function deg2rad(d) result(r)
    real(dp), intent(in) :: d
    real(dp) :: r

    r = d * rad
end function deg2rad

elemental function rad2deg(r) result(d)
    real(dp), intent(in) :: r
    real(dp) :: d

    d = r * deg
end function rad2deg

! Function: isapprox
!
! Compares two double precision numbers.
!
! Parameters:
!   a - The first double precision number.
!   b - The second double presicion number.
!   rtol - Relative tolerance (optional, default: 1e-9).
!   atol - Absolute tolerance (optional, default: 0.0).
!
! Returns:
!   True if the numbers are approximately equal.

pure function isapprox_scalar(x, y, rtol, atol) result(res)
    real(dp), intent(in) :: x
    real(dp), intent(in) :: y
    real(dp), intent(in), optional :: rtol
    real(dp), intent(in), optional :: atol

    logical :: res
    real(dp) :: rtol_
    real(dp) :: atol_

    rtol_ = sqrt(eps)
    if (present(rtol)) rtol_ = rtol
    atol_ = 0._dp
    if (present(atol)) atol_ = atol

    res = abs(x - y) <= atol_ + rtol_ * max(abs(x), abs(y))
end function isapprox_scalar

pure function isapprox_vector(x, y, rtol, atol) result(res)
    real(dp), dimension(:), intent(in) :: x
    real(dp), dimension(:), intent(in) :: y
    real(dp), intent(in), optional :: rtol
    real(dp), intent(in), optional :: atol

    logical :: res
    real(dp) :: rtol_
    real(dp) :: atol_

    rtol_ = sqrt(eps)
    if (present(rtol)) rtol_ = rtol
    atol_ = 0._dp
    if (present(atol)) atol_ = atol

    res = norm(x - y) <= atol_ + rtol_ * max(norm(x), norm(y))
end function isapprox_vector

pure function isapprox_matrix(x, y, rtol, atol) result(res)
    real(dp), dimension(:,:), intent(in) :: x
    real(dp), dimension(:,:), intent(in) :: y
    real(dp), intent(in), optional :: rtol
    real(dp), intent(in), optional :: atol

    logical :: res
    real(dp) :: rtol_
    real(dp) :: atol_

    rtol_ = sqrt(eps)
    if (present(rtol)) rtol_ = rtol
    atol_ = 0._dp
    if (present(atol)) atol_ = atol

    res = norm(x - y) <= atol_ + rtol_ * max(norm(x), norm(y))
end function isapprox_matrix

pure function isin_vector_dp(el, arr) result(res)
    real(dp), intent(in) :: el
    real(dp), dimension(:), intent(in) :: arr
    logical :: res

    integer :: i

    res = .false.

    do i = 1, size(arr)
        if (isapprox(el, arr(i))) then
            res = .true.
            exit
        end if
    end do
end function isin_vector_dp

pure function isin_matrix_dp(el, arr) result(res)
    real(dp), intent(in) :: el
    real(dp), dimension(:,:), intent(in) :: arr
    logical :: res

    integer :: i
    integer :: j
    integer, dimension(2) :: shp

    shp = shape(arr)
    res = .false.

    do i = 1, shp(1)
        do j = 1, shp(2)
            if (isapprox(el, arr(i,j))) then
                res = .true.
                exit
            end if
        end do
    end do
end function isin_matrix_dp

pure function isin_vector_int(el, arr) result(res)
    integer, intent(in) :: el
    integer, dimension(:), intent(in) :: arr
    logical :: res

    res = any(el == arr)
end function isin_vector_int

pure function isin_matrix_int(el, arr) result(res)
    integer, intent(in) :: el
    integer, dimension(:,:), intent(in) :: arr
    logical :: res

    res = any(el == arr)
end function isin_matrix_int

! Function: cross
!   Calculates cross product.
!
! Parameters:
!   v1 - First vector
!   v2 - Second vector
!
! Returns:
!   prod - Resulting vector
pure function cross(a, b)
  real(dp), dimension(3) :: cross
  real(dp), dimension(:), intent(in) :: a, b

  cross(1) = a(2) * b(3) - a(3) * b(2)
  cross(2) = a(3) * b(1) - a(1) * b(3)
  cross(3) = a(1) * b(2) - a(2) * b(1)
end function cross

! Function: unitmat
!   Generates a unity matrix.
!
! Parameters:
!   n - Matrix rank
!
! Returns:
!   unity - Unity matrix of rank n
pure function unitmat(n) result(unity)
    integer, intent(in) :: n
    real(dp) :: unity(n, n)
    integer :: i
    unity = 0._dp
    do i = 1, n
        unity(i,i) = 1
    end do
end function unitmat

! Function: mod2pi
!   Truncates angle to be within 0 and 2*pi.
!
! Parameters:
!   angin - Angle in radians
!
! Returns:
!   ang - Truncated angle in radians
elemental function mod2pi(angin) result(ang)
    real(dp), intent(in) :: angin
    real(dp) :: ang
    ang = mod(angin, twopi)
    if (ang < 0._dp) then
        ang = ang + twopi
    end if
end function mod2pi

! Function: vecang
!   Calculates the angle between two vectors.
!
! Parameters:
!   v1 - First vector
!   v2 - Second vector
!
! Returns:
!   ang - Angle in radians
pure function vecang(v1, v2) result(ang)
    real(dp), intent(in) :: v1(:)
    real(dp), intent(in) :: v2(:)
    real(dp) :: ang
    real(dp) :: v1norm
    real(dp) :: v2norm

    v1norm = norm(v1)
    v2norm = norm(v2)

    if (v1norm < eps) then
        ang = 0._dp
        return
    elseif(v2norm < eps) then
        ang = 0._dp
        return
    end if

    ang = dot_product(v1, v2)
    ang = ang/v1norm/v2norm
    ang = min(ang, 1._dp)
    ang = max(ang, -1._dp)
    ang = acos(ang)
end function vecang

! Function: polcart
!   Converts from polar to cartesian coordinates.
!
! Parameters:
!   p - Vector in polar coordinates
!
! Returns:
!   x - Vector in cartesian coordinates
pure function polcart(p) result(x)
    real(dp), intent(in) :: p(:)
    real(dp) :: x(3)
    real(dp) :: pxy

    x(3) = sin(p(3)) * p(1)
    pxy = cos(p(3)) * p(1)
    x(1) = pxy * cos(p(2))
    x(2) = pxy * sin(p(2))
end function polcart

! Function: cartpol
!   Converts from cartesian to polar coordinates.
!
! Parameters:
!   x - Vector in cartesian coordinates
!
! Returns:
!   p - array Vector in polar coordinates
pure function cartpol(x) result(p)
    real(dp), intent(in) :: x(:)
    real(dp) :: p(3)

    p(1) = norm(x)
    if (p(1) < eps) then
        p(2) = 0._dp
        p(3) = 0._dp
    else
        p(2) = atan2(x(2), x(1))
        p(3) = asin(x(3) / p(1))
    end if
end function cartpol

! Function: norm
!   Calculates vector norm.
!
! Parameters:
!   x - Vector
!
! Returns:
!   l - Vector norm
pure function norm_vector(x) result(l)
    real(dp), dimension(:), intent(in) :: x
    real(dp) :: l

    l = sqrt(sum(x**2))
end function norm_vector

pure function norm_matrix(x) result(l)
    real(dp), dimension(:,:), intent(in) :: x
    real(dp) :: l

    l = sqrt(sum(x**2))
end function norm_matrix

! Function: getsign
!   Calculates sign
!
! Parameters:
!   x - Scalar
!
! Returns:
!   l - Sign : neg = -1, 0 = 0, pos = +1
pure function getsign(x) result(l)
    real(dp), intent(in) :: x
    real(dp) :: l

    if (isapprox(x, 0._dp)) then
        l = 0
    else
        l = x/abs(x)
    end if

end function getsign

! Function: binsearch
!   Binary search within a sorted array.
!
! Parameters:
!   val - Search parameter
!   x - Array sorted in ascending/descending order
!   desc - Flag for descending order
!
! Returns:
!   ind - Index of the first element larger/smaller than val
function binsearch(val, x, desc, err) result(ind)
    real(dp), intent(in) :: val
    real(dp), intent(in), dimension(:) :: x
    logical, intent(in), optional :: desc
    type(exception), intent(inout), optional :: err
    integer :: ind
    integer :: imin
    integer :: imax
    integer :: imid
    integer :: n
    logical :: desc_
    type(exception) :: err_

    n = size(x)
    imin = 1
    imax = n
    ind = 0

    desc_ = .false.
    if (present(desc)) desc_ = desc

    if (.not.issorted(x, desc_)) then
        err_ = error("Array is not monotonically "//&
            "increasing/decreasing.", "binsearch", __FILE__, __LINE__)
        if (present(err)) then
            err = err_
            return
        else
            call raise(err_)
        end if
    end if

    if (desc_) then
        if (val < x(n)) then
            ind = n
            return
        end if
    else
        if (val > x(n)) then
            ind = n
            return
        end if
    end if

    do while (imin < imax)
        imid = imin + (imax - imin) / 2
        if (desc_) then
            if (val <= x(imid)) then
                imin = imid + 1
            else
                imax = imid
            end if
        else
            if (val >= x(imid)) then
                imin = imid + 1
            else
                imax = imid
            end if
        end if
    end do

    ind = imin
end function binsearch

! Function: interp
!   Linear and zero-order interpolation.
!
! Parameters:
!   x - Intermediate point
!   xt - Vector of support points.
!   yt - Vector of data values.
!   kind - Interpolation method ("zero": zero-order ascending, "descnd":
!   zero-order descending, "linear": linear interpolation)
!   left - Extrapolation value for x <= xt(1) (optional, default: yt(1))
!   right - Extrapolation value for x >= xt(end) (optional, default: yt(end))
!
! Returns:
!   y - Interpolated value
function interp(x, xt, yt, kind, left, right, err) result(y)
    real(dp), intent(in) :: x
    real(dp), dimension(:), intent(in) :: xt
    real(dp), dimension(:), intent(in) :: yt
    character(len=*), intent(in) :: kind
    real(dp), intent(in), optional :: left
    real(dp), intent(in), optional :: right
    type(exception), intent(out), optional :: err
    real(dp) :: y

    type(exception) :: err_
    integer :: ind
    integer :: n
    real(dp) :: x0
    real(dp) :: x1
    real(dp) :: y0
    real(dp) :: y1
    real(dp) :: left_
    real(dp) :: right_

    y = 0._dp

    n = size(xt)

    if (present(left)) then
        left_ = left
    else
        left_ = yt(1)
    end if

    if (present(right)) then
        right_ = right
    else
        right_ = yt(n)
    end if

    select case (kind)
    case ("zero")
        if (x <= xt(1)) then
            y = left_
            return
        else if (x >= xt(n)) then
            y = right_
            return
        end if
        ind = binsearch(x, xt, err=err_)
        if (iserror(err_)) then
            call catch(err_, "interp", __FILE__, __LINE__)
            if (present(err)) then
                err = err_
            else
                call raise(err_)
            end if
        end if
        y = yt(ind-1)
    case ("descnd")
        if (x >= xt(1)) then
            y = left_
            return
        else if (x <= xt(n)) then
            y = right_
            return
        end if
        ind = binsearch(x, xt, desc=.true., err=err_)
        if (iserror(err_)) then
            call catch(err_, "interp", __FILE__, __LINE__)
            if (present(err)) then
                err = err_
            else
                call raise(err_)
            end if
        end if
        y = yt(ind-1)
    case ("linear")
        if (x <= xt(1)) then
            x0 = xt(1)
            x1 = xt(2)
            y0 = yt(1)
            y1 = yt(2)
        else if (x >= xt(n)) then
            x0 = xt(n-1)
            x1 = xt(n)
            y0 = yt(n-1)
            y1 = yt(n)
        else
            ind = binsearch(x, xt, err=err_)
            if (iserror(err_)) then
                call catch(err_, "interp", __FILE__, __LINE__)
                if (present(err)) then
                    err = err_
                else
                    call raise(err_)
                end if
            end if
            x0 = xt(ind-1)
            x1 = xt(ind)
            y0 = yt(ind-1)
            y1 = yt(ind)
        end if

        y = y0 + (y1 - y0)/(x1 - x0) * (x - x0)
    end select
end function interp

! Function: findroot
!   Find the root of a function within a given interval using Brent's method.
!
! Source:
!   https://github.com/scipy/scipy/blob/v0.12.0/scipy/optimize/Zeros/brentq.c
!
!   Fortran90 port
!
! Parameters:
!   fun - User function: y = fun(x, rpar, ipar)
!   xa - Lower bound of the interval
!   xb - Upper bound of the interval
!   rpar - Real valued function parameters
!   ipar - Integer valued function parameters
!   xtol - Absolute tolerance (optional, default: 1e-12)
!   rtol - Relative tolerance (optional, default: 1.48e-8)
!   max_iter - Maximum number of iterations (optional, default: 100)
!
! Returns:
!   root - Root of the function.
function findroot_brent(f, xa, xb, rpar, ipar, xtol, rtol, max_iter, err)&
        result(root)
    procedure(func) :: f
    real(dp), intent(in) :: xa
    real(dp), intent(in) :: xb
    real(dp), dimension(:), intent(in), optional :: rpar
    integer, dimension(:), intent(in), optional :: ipar
    real(dp), intent(in), optional :: xtol
    real(dp), intent(in), optional :: rtol
    integer, intent(in), optional :: max_iter
    type(exception), intent(inout), optional :: err

    real(dp) :: root

    type(exception) :: err_

    real(dp) :: fpre
    real(dp) :: fcur
    real(dp) :: fblk
    real(dp) :: xpre
    real(dp) :: xcur
    real(dp) :: xblk
    real(dp) :: xtol_
    real(dp) :: rtol_
    real(dp) :: spre
    real(dp) :: scur
    real(dp) :: sbis
    real(dp) :: stry
    real(dp) :: dpre
    real(dp) :: dblk
    real(dp) :: tol
    real(dp), dimension(:), allocatable :: rpar_
    integer, dimension(:), allocatable :: ipar_

    integer :: max_iter_
    integer :: i

    if (present(rpar)) then
        rpar_ = rpar
    else
        allocate(rpar_(1))
        rpar_ = 0._dp
    end if

    if (present(ipar)) then
        ipar_ = ipar
    else
        allocate(ipar_(1))
        ipar_ = 0
    end if

    i = 0
    fblk = 0._dp
    xblk = 0._dp
    xpre = xa
    xcur = xb
    spre = 0._dp
    scur = 0._dp
    root = 0._dp

    max_iter_ = 100
    if (present(max_iter)) max_iter_ = max_iter

    xtol_ = 1e-6_dp
    if (present(xtol)) xtol_ = xtol

    rtol_ = sqrt(eps)
    if (present(rtol)) rtol_ = rtol

    fpre = f(xpre, rpar_, ipar_)
    fcur = f(xcur, rpar_, ipar_)

    if (fpre * fcur > 0._dp) then
        err_ = error("Root not in bracket.", "findroot", __FILE__, __LINE__)
        if (present(err)) then
            err = err_
        else
            call raise(err_)
        end if
    end if

    if (isapprox(fpre, 0._dp)) then
        root = xpre
        return
    end if

    if (isapprox(fcur, 0._dp)) then
        root = xcur
        return
    end if

    do while (i < max_iter_)
        if (fpre * fcur < 0._dp) then
            xblk = xpre
            fblk = fpre
            spre = xcur - xpre
            scur = xcur - xpre
        end if

        if (abs(fblk) < abs(fcur)) then
            xpre = xcur
            xcur = xblk
            xblk = xpre
            fpre = fcur
            fcur = fblk
            fblk = fpre
        end if

        tol = xtol_ + rtol_ * abs(xcur)
        sbis = (xblk - xcur) / 2._dp

        if (isapprox(fcur, 0._dp).or.(abs(sbis) < tol)) then
            root = xcur
            return
        end if

        if ((abs(spre) > tol).and.(abs(fcur) < abs(fpre))) then
            if (isapprox(xpre, xblk)) then
                stry = -fcur*(xcur - xpre)/(fcur - fpre)
            else
                dpre = (fpre - fcur)/(xpre - xcur)
                dblk = (fblk - fcur)/(xblk - xcur)
                stry = -fcur*(fblk*dblk - fpre*dpre)&
                    /(dblk*dpre*(fblk - fpre))
            end if

            if (2*abs(stry) < min(abs(spre), 3*abs(sbis) - tol)) then
                spre = scur
                scur = stry
            else
                spre = sbis
                scur = sbis
            end if
        else
            spre = sbis
            scur = sbis
        end if

        xpre = xcur
        fpre = fcur

        if (abs(scur) > tol) then
            xcur = xcur + scur
        else
            if (sbis > 0._dp) xcur = xcur + tol
            if (sbis < 0._dp) xcur = xcur - tol
        end if

        fcur = f(xcur, rpar_, ipar_)

        i = i + 1
    end do
end function findroot_brent

function findroot_steffensen(fun, x, rpar, ipar, rtol, atol, maxiter, err) result(root)
    procedure(func) :: fun
    real(dp), intent(in) :: x
    real(dp), dimension(:), intent(in), optional :: rpar
    integer, dimension(:), intent(in), optional :: ipar
    real(dp), intent(in), optional :: rtol
    real(dp), intent(in), optional :: atol
    integer, intent(in), optional :: maxiter
    type(exception), intent(inout), optional :: err
    real(dp) :: root

    real(dp) :: f
    real(dp) :: f1
    real(dp) :: f2
    real(dp) :: atol_
    real(dp) :: rtol_
    real(dp), dimension(:), allocatable :: rpar_
    integer, dimension(:), allocatable :: ipar_
    logical :: converged
    type(exception) :: err_
    integer :: i

    integer :: maxiter_

    maxiter_ = 1000
    if (present(maxiter)) maxiter_ = maxiter

    atol_ = 0._dp
    if (present(atol)) atol_ = atol

    rtol_ = sqrt(eps)
    if (present(rtol)) rtol_ = rtol

    if (present(rpar)) then
        rpar_ = rpar
    else
        allocate(rpar_(1))
        rpar_ = 0._dp
    end if

    if (present(ipar)) then
        ipar_ = ipar
    else
        allocate(ipar_(1))
        ipar_ = 0
    end if

    root = x
    converged = .false.

    do i = 1, maxiter_
        f1 = root + fun(root, rpar_, ipar_)
        f2 = f1 + fun(f1, rpar_, ipar_)
        f = root - (f1 - root)**2 / (f2 - 2 * f1 + root)
        converged = isapprox(f, root, rtol=rtol_, atol=atol_)
        if (converged) exit
        root = f
    end do

    if (.not.converged) then
        err_ = error("Failed to converge.", "findroot_steffensen", __FILE__, __LINE__)
        if (present(err)) then
            err = err_
            return
        else
            call raise(err_)
        end if
    end if
end function findroot_steffensen

! Function: issorted
!   Tests whether an array is monotonically increasing or decreasing.
!
! Parameters:
!   v - Array
!   desc - Flag for descending order (optional, default: false)
!
! Returns:
!   True if the array is sorted.
function issorted(v, desc) result(sorted)
    real(dp), dimension(:), intent(in) :: v
    logical, intent(in), optional :: desc
    logical :: sorted
    logical :: desc_

    integer :: i

    desc_ = .false.
    if (present(desc)) desc_ = desc

    sorted = .true.
    do i = 2, size(v)
        if (.not.(desc_)) then
            if (v(i-1) >= v(i)) then
               sorted = .false.
               return
            end if
        else
            if (v(i-1) <= v(i)) then
               sorted = .false.
               return
            end if
        end if
    end do
end function issorted

function greatcircle(lat1, lat2, lon1, lon2, r) result(d)
    real(dp), intent(in) :: lat1
    real(dp), intent(in) :: lat2
    real(dp), intent(in) :: lon1
    real(dp), intent(in) :: lon2
    real(dp), intent(in) :: r
    real(dp) :: d
    real(dp) :: dlon
    real(dp) :: dsigma
    real(dp) :: a
    real(dp) :: b

    dlon = abs(lon2 - lon1)
    a = sqrt((cos(lat2)*sin(dlon))**2 + (cos(lat1)*sin(lat2) - sin(lat1)*cos(lat2)*cos(dlon))**2)
    b = sin(lat1)*sin(lat2) + cos(lat1)*cos(lat2)*cos(dlon)
    dsigma = atan2(a, b)
    d = r*dsigma
end function greatcircle

function linspace(x1, xn, n) result(res)
    real(dp), intent(in) :: x1
    real(dp), intent(in) :: xn
    integer, intent(in) :: n
    integer :: i

    real(dp), dimension(n) :: res

    do i = 1, n-1
        res(i) = x1 + (xn - x1) * real(i-1, dp) / real(n-1, dp)
    end do
    res(size(res)) = xn
end function linspace

end module math
