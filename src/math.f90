module math

use types, only: dp
use exceptions

implicit none

private

public :: pih, pi, twopi, deg, rad,&
    eps, cross, unitmat, mod2pi,&
    vecang, polcart, cartpol, norm,&
    getsign, binsearch,&
    interp, findroot, isclose,&
    issorted, greatcircle, deg2rad, rad2deg

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

contains

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

! Function: is_close
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

logical function isclose(a, b, rtol, atol)
    real(dp), intent(in) :: a
    real(dp), intent(in) :: b
    real(dp), intent(in), optional :: rtol
    real(dp), intent(in), optional :: atol

    real(dp) :: rtol_
    real(dp) :: atol_
    real(dp) :: diff

    rtol_ = 1e-9_dp
    if (present(rtol)) rtol_ = rtol
    atol_ = 0._dp
    if (present(atol)) atol_ = atol

    diff = abs(b - a)

    isclose = (((diff <= abs(rtol_*b)).or.(diff <= abs(rtol_*a)))&
        .or.(diff <= atol_))
end function isclose

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
pure function mod2pi(angin) result(ang)
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
pure function norm(x) result(l)
    real(dp), intent(in) :: x(:)
    real(dp) :: l

    l = sqrt(sum(x**2))
end function norm

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

    if (x == 0) then
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

    n = size(x)
    imin = 1
    imax = n

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
function findroot(fun, xa, xb, rpar, ipar, xtol, rtol, max_iter, err)&
        result(root)
    interface
        real(dp) function fun(x, rpar, ipar)
            use types, only: dp
            real(dp), intent(in) :: x
            real(dp), dimension(:), intent(in) :: rpar
            integer, dimension(:), intent(in) :: ipar
        end function fun
    end interface
    real(dp), intent(in) :: xa
    real(dp), intent(in) :: xb
    real(dp), dimension(:), intent(in) :: rpar
    integer, dimension(:), intent(in) :: ipar
    real(dp), intent(in), optional :: xtol
    real(dp), intent(in), optional :: rtol
    integer, intent(in), optional :: max_iter
    type(exception), intent(out), optional :: err

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

    integer :: max_iter_
    integer :: i

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

    fpre = fun(xpre, rpar, ipar)
    fcur = fun(xcur, rpar, ipar)

    if (fpre * fcur > 0._dp) then
        err_ = error("Root not in bracket.", "findroot", __FILE__, __LINE__)
        if (present(err)) then
            err = err_
        else
            call raise(err_)
        end if
    end if

    if (fpre == 0._dp) then
        root = xpre
        return
    end if

    if (fcur == 0._dp) then
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

        if ((fcur == 0._dp).or.(abs(sbis) < tol)) then
            root = xcur
            return
        end if

        if ((abs(spre) > tol).and.(abs(fcur) < abs(fpre))) then
            if (xpre == xblk) then
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

        fcur = fun(xcur, rpar, ipar)

        i = i + 1
    end do
end function findroot

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

end module math
