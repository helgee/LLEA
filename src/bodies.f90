module bodies

use types, only: dp
use epochs, only: epoch, seconds_per_century, seconds, seconds_per_day
use math, only: mod2pi, pih

implicit none

private

public :: body, rightascension, rightascensionrate, declination, declinationrate, rotationangle, rotationrate

type body
    real(dp) :: mu
    real(dp) :: j2
    character(len=6) :: bodytype
    character(len=7) :: parent
    real(dp), dimension(3) :: radii
    integer :: id
    real(dp), dimension(3) :: ra
    real(dp), dimension(3) :: dec
    real(dp), dimension(3) :: ww
    real(dp), dimension(:), allocatable :: a
    real(dp), dimension(:), allocatable :: d
    real(dp), dimension(:), allocatable :: w
    real(dp), dimension(:), allocatable :: theta0
    real(dp), dimension(:), allocatable :: theta1
end type body

contains

pure function rightascension(b, ep) result(angle)
    type(body), intent(in) :: b
    type(epoch), intent(in) :: ep

    real(dp) :: angle
    real(dp) :: t

    t = seconds(ep)
    angle = mod2pi(b%ra(1) + b%ra(2) * t / seconds_per_century + b%ra(3) * t**2 / seconds_per_century**2 &
        + sum(b%a * sin(b%theta0 + b%theta1 * t / seconds_per_century)))
end function rightascension

pure function declination(b, ep) result(angle)
    type(body), intent(in) :: b
    type(epoch), intent(in) :: ep

    real(dp) :: angle
    real(dp) :: t

    t = seconds(ep)
    angle = mod2pi(b%dec(1) + b%dec(2) * t / seconds_per_century + b%dec(3) * t**2 / seconds_per_century**2 &
        + sum(b%d * cos(b%theta0 + b%theta1 * t / seconds_per_century)))
end function declination

pure function rotationangle(b, ep) result(angle)
    type(body), intent(in) :: b
    type(epoch), intent(in) :: ep

    real(dp) :: angle
    real(dp) :: t

    t = seconds(ep)
    angle = mod2pi(b%ww(1) + b%ww(2) * t / seconds_per_day + b%ww(3) * t**2 / seconds_per_day**2 &
        + sum(b%w * sin(b%theta0 + b%theta1 * t / seconds_per_century)))
end function rotationangle

pure function rightascensionrate(b, ep) result(rate)
    type(body), intent(in) :: b
    type(epoch), intent(in) :: ep

    real(dp) :: rate
    real(dp) :: t

    t = seconds(ep)
    rate = b%ra(2) / seconds_per_century + 2 * b%ra(3) * t / seconds_per_century**2 &
        + sum(b%a * b%theta1 / seconds_per_century * cos(b%theta0 + b%theta1 * t / seconds_per_century))
end function rightascensionrate

pure function declinationrate(b, ep) result(rate)
    type(body), intent(in) :: b
    type(epoch), intent(in) :: ep

    real(dp) :: rate
    real(dp) :: t

    t = seconds(ep)
    rate = b%dec(2) / seconds_per_century + 2 * b%dec(3) * t / seconds_per_century**2 &
        + sum(b%d * b%theta1 / seconds_per_century * sin(b%theta0 + b%theta1 * t / seconds_per_century))
end function declinationrate

pure function rotationrate(b, ep) result(rate)
    type(body), intent(in) :: b
    type(epoch), intent(in) :: ep

    real(dp) :: rate
    real(dp) :: t

    t = seconds(ep)
    rate = b%ww(2) / seconds_per_day + 2 * b%ww(3) * t / seconds_per_day**2 &
        + sum(b%w * b%theta1 / seconds_per_century * cos(b%theta0 + b%theta1 * t / seconds_per_century))
end function rotationrate

subroutine iaumatrix(b, ep, m, dm)
    type(body), intent(in) :: b
    type(epoch), intent(in) :: ep
    real(dp), dimension(3,3), intent(out) :: m
    real(dp), dimension(3,3), intent(out) :: dm

    real(dp) :: alpha
    real(dp) :: dalpha
    real(dp) :: delta
    real(dp) :: ddelta
    real(dp) :: omega
    real(dp) :: domega
    real(dp) :: phi
    real(dp) :: xi

    alpha = rightascension(b, ep)
    dalpha = rightascensionrate(b, ep)
    delta = declination(b, ep)
    ddelta = declinationrate(b, ep)
    omega = rotationangle(b, ep)
    domega = rotationrate(b, ep)
    phi = alpha + pih
    xi = pih - delta

    m = rotationmatrix("313", phi, xi, omega)
    dm = ratematrix("313", phi, dalpha, xi, -ddelta, omega, domega)
end subroutine iaumatrix

end module bodies
