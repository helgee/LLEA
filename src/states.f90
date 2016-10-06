module states

use types, only: dp
use bodies, only: body, iaumatrix
use epochs, only: epoch
use exceptions
use math, only: pi, cross

implicit none

integer, parameter :: framelen = 8

type state
    type(epoch) :: ep
    real(dp), dimension(6) :: rv
    character(len=framelen) :: frame
    type(body) :: center
end type state

interface keplerian
    module procedure keplerian_vectors
    module procedure keplerian_state
end interface keplerian

private

public :: state, rotate_inplace, rotate

contains

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

! TODO: Port Julia implementation
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

pure function keplerian_vectors(r, v, mu) result(ele)
    double precision, dimension(:), intent(in) :: r
    double precision, dimension(:), intent(in) :: v
    double precision, intent(in) :: mu
    double precision, dimension(6) :: ele

    double precision :: r_mag, v_mag, h_mag, n_mag, xi
    double precision, dimension(3) :: h, n, k, e

    r_mag = norm2(r)
    v_mag = norm2(v)
    h = cross(r,v)
    h_mag = norm2(h)
    k = [0d0, 0d0, 1d0]
    n = cross(k, h)
    n_mag = norm2(n)
    xi = v_mag**2/2 - mu/r_mag
    e = ((v_mag**2 - mu/r_mag)*r - v*dot_product(r,v))/mu
    ele(2) = norm2(e)
    if (ele(2) /= 1.0) then
        ele(1) = -mu/(2*xi)
    else
        ele(1) = h_mag**2/mu
    end if
    ele(3) = acos(h(3)/h_mag)
    ele(4) = acos(n(1)/n_mag)
    ele(5) = acos(dot_product(n,e)/(ele(2)*n_mag))
    ele(6) = acos(dot_product(e,r)/(ele(2)*r_mag))
    if (n(2) < 0) then
        ele(4) = 2*pi - ele(4)
    end if
    if (e(3) < 0) then
        ele(5) = 2*pi - ele(5)
    end if
    if (dot_product(r,v) < 0) then
        ele(6) = 2*pi - ele(6)
    end if
end function keplerian_vectors

pure function keplerian_state(s) result(ele)
    type(state), intent(in) :: s
    double precision, dimension(6) :: ele

    ele = keplerian_vectors(s%rv(1:3), s%rv(4:6), s%center%mu)
end function keplerian_state

end module states
