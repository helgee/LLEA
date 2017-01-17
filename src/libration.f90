module libration

use exceptions
use math, only: findroot
use types, only: dp

implicit none

private

public :: librationdist

contains

function l1dist(x, rpar) result(res)
    real(dp), intent(in) :: x
    real(dp), dimension(:), intent(in) :: rpar
    real(dp) :: res
    real(dp) :: mu

    mu = rpar(1)
    res = x**3._dp - mu * (1._dp - x)**2._dp / (3._dp - 2._dp * mu - x * (3._dp - mu - x))
end function l1dist

function l2dist(x, rpar) result(res)
    real(dp), intent(in) :: x
    real(dp), dimension(:), intent(in) :: rpar
    real(dp) :: res
    real(dp) :: mu

    mu = rpar(1)
    res = x**3._dp - mu * (1._dp + x)**2._dp / (3._dp - 2._dp * mu + x * (3._dp - mu + x))
end function l2dist

function l3dist(x, rpar) result(res)
    real(dp), intent(in) :: x
    real(dp), dimension(:), intent(in) :: rpar
    real(dp) :: res
    real(dp) :: mu

    mu = rpar(1)
    res = x**3._dp - (1._dp - mu) * (1._dp + x)**2._dp / (1._dp + 2._dp * mu + x * (2._dp + mu + x))
end function l3dist

function librationdist(mu_prime, mu_sec, point) result(res)
    real(dp), intent(in) :: mu_prime
    real(dp), intent(in) :: mu_sec
    character(len=*), intent(in) :: point
    real(dp) :: res

    real(dp), dimension(1) :: rpar
    real(dp) :: mu_lib
    real(dp) :: initial
    type(exception) :: err

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
        err = error("Libration point name must be 'L1', 'L2', or 'L3'.", "librationdist", __FILE__, __LINE__)
        call raise(err)
    end select
end function librationdist

end module libration
