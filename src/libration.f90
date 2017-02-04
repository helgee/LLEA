module libration

use bodies, only: body
use ephemerides, only: ephem, getstate
use epochs, only: epoch
use exceptions
use math, only: findroot, cross, norm
use types, only: dp

implicit none

private

public :: librationdist, gcrftolib, libtogcrf

contains

function l1dist(x, rpar, ipar) result(res)
    real(dp), intent(in) :: x
    real(dp), dimension(:), intent(in) :: rpar
    integer, dimension(:), intent(in) :: ipar
    real(dp) :: res
    real(dp) :: mu

    mu = rpar(1)
    res = x**3._dp - mu * (1._dp - x)**2._dp / (3._dp - 2._dp * mu - x * (3._dp - mu - x))
end function l1dist

function l2dist(x, rpar, ipar) result(res)
    real(dp), intent(in) :: x
    real(dp), dimension(:), intent(in) :: rpar
    integer, dimension(:), intent(in) :: ipar
    real(dp) :: res
    real(dp) :: mu

    mu = rpar(1)
    res = x**3._dp - mu * (1._dp + x)**2._dp / (3._dp - 2._dp * mu + x * (3._dp - mu + x))
end function l2dist

function l3dist(x, rpar, ipar) result(res)
    real(dp), intent(in) :: x
    real(dp), dimension(:), intent(in) :: rpar
    integer, dimension(:), intent(in) :: ipar
    real(dp) :: res
    real(dp) :: mu

    mu = rpar(1)
    res = x**3._dp - (1._dp - mu) * (1._dp + x)**2._dp / (1._dp + 2._dp * mu + x * (2._dp + mu + x))
end function l3dist

function librationdist(mu_prime, mu_sec, point, err) result(res)
    real(dp), intent(in) :: mu_prime
    real(dp), intent(in) :: mu_sec
    character(len=*), intent(in) :: point
    type(exception), intent(inout), optional :: err
    real(dp) :: res

    real(dp), dimension(1) :: rpar
    real(dp) :: mu_lib
    real(dp) :: initial
    type(exception) :: err_

    res = 0._dp
    mu_lib = 1._dp / (1._dp + mu_prime / mu_sec)
    rpar = mu_lib

    select case (point)
    case ("L1")
        initial = (mu_lib / 3._dp / (1._dp - mu_lib))**(1._dp/3._dp)
        res = findroot(l1dist, initial, rpar=rpar)
    case ("L2")
        initial = (mu_lib / 3._dp / (1._dp - mu_lib))**(1._dp/3._dp)
        res = findroot(l2dist, initial, rpar=rpar)
    case ("L3")
        initial = 1._dp - (7._dp * mu_lib / 12._dp)
        res = findroot(l3dist, initial, rpar=rpar)
    case default
        err_ = error("Libration point name must be 'L1', 'L2', or 'L3'.", "librationdist", __FILE__, __LINE__)
        if (present(err)) then
            err = err_
        else
            call raise(err_)
        end if
    end select
end function librationdist

subroutine rotationparams(rv, m, rm, omega)
    real(dp), dimension(:), intent(in) :: rv
    real(dp), dimension(:,:), intent(out) :: m
    real(dp), intent(out) :: rm
    real(dp), intent(out) :: omega

    real(dp), dimension(3) :: omegav
    real(dp) :: om

    m = 0._dp
    rm = norm(rv(1:3))
    m(:,1) = rv(1:3) / rm
    omegav = cross(rv(1:3), rv(4:6))
    om = norm(omegav)
    m(:,3) = omegav / om
    m(:,2) = cross(m(:,3), m(:,1))
    omega = om / rm**2
end subroutine rotationparams

function gcrftolib(rv, ep, primary, secondary, point, err) result(lib)
    real(dp), dimension(:), intent(in) :: rv
    type(epoch), intent(in) :: ep
    type(body), intent(in) :: primary
    type(body), intent(in) :: secondary
    character(len=*), intent(in) :: point
    type(exception), intent(inout), optional :: err
    real(dp), dimension(6) :: lib

    type(exception) :: err_
    real(dp) :: gam
    real(dp), dimension(6) :: sec
    real(dp), dimension(3,3) :: m
    real(dp) :: rm
    real(dp) :: omega

    gam = librationdist(primary%mu, secondary%mu, point, err=err_)
    if (iserror(err_)) then
        call catch(err_, "gcrftolib", __FILE__, __LINE__)
        if (present(err)) then
            err = err_
            return
        else
            call raise(err_)
        end if
    end if
    sec = ephem%state(ep, secondary%id, primary%id/100, err_)
    if (iserror(err_)) then
        call catch(err_, "gcrftolib", __FILE__, __LINE__)
        if (present(err)) then
            err = err_
            return
        else
            call raise(err_)
        end if
    end if
    lib = rv - sec * (1 + gam)

    call rotationparams(sec, m, rm, omega)
    lib(1:3) = matmul(transpose(m), lib(1:3))
    lib(4:6) = matmul(transpose(m), lib(4:6)) - [-omega * lib(2), omega * lib(1), 0._dp]
    lib(1:3) = lib(1:3) / rm
    lib(4:6) = lib(4:6) / (rm * omega)
end function gcrftolib

function libtogcrf(lib, ep, primary, secondary, point, err) result(rv)
    real(dp), dimension(:), intent(in) :: lib
    type(epoch), intent(in) :: ep
    type(body), intent(in) :: primary
    type(body), intent(in) :: secondary
    character(len=*), intent(in) :: point
    type(exception), intent(inout), optional :: err
    real(dp), dimension(6) :: res

    type(exception) :: err_
    real(dp) :: gam
    real(dp), dimension(6) :: sec
    real(dp), dimension(6) :: rv
    real(dp), dimension(3,3) :: m
    real(dp) :: rm
    real(dp) :: omega

    gam = librationdist(primary%mu, secondary%mu, point, err=err_)
    if (iserror(err_)) then
        call catch(err_, "libtogcrf", __FILE__, __LINE__)
        if (present(err)) then
            err = err_
            return
        else
            call raise(err_)
        end if
    end if
    sec = ephem%state(ep, secondary%id, primary%id/100, err_)
    if (iserror(err_)) then
        call catch(err_, "libtogcrf", __FILE__, __LINE__)
        if (present(err)) then
            err = err_
            return
        else
            call raise(err_)
        end if
    end if

    call rotationparams(sec, m, rm, omega)
    rv = lib
    rv(4:6) = rv(4:6) + [-lib(2), lib(1), 0._dp]
    rv(1:3) = matmul(m, rv(1:3))
    rv(4:6) = matmul(m, rv(4:6))
    rv(1:3) = rv(1:3) * rm
    rv(4:6) = rv(4:6) * omega * rm
    rv = rv + sec * (1 + gam)
end function libtogcrf

end module libration
