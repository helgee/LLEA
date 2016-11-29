!
! Copyright (c) 2016 Helge Eichhorn
!
! This Source Code Form is subject to the terms of the Mozilla Public
! License, v. 2.0. If a copy of the MPL was not distributed with this
! file, You can obtain one at http://mozilla.org/MPL/2.0/.
!
module events

use containers, only: parameters
use exceptions
use math, only: eps, isapprox
use types, only: dp

implicit none

private

public :: event, update, findevent

type, abstract :: update
contains
    procedure(abstract_apply), deferred :: apply
end type update

type, abstract :: detector
contains
    procedure(abstract_haspassed), deferred :: haspassed
    procedure(abstract_detect), deferred :: detect
    procedure :: find => findevent
end type detector

type event
    logical :: done = .false.
    real(dp), dimension(:), allocatable :: t
    class(detector), allocatable :: det
    class(update), allocatable :: up
end type event

abstract interface
    subroutine abstract_apply(this, y, p)
        import :: update, dp, parameters
        class(update), intent(in) :: this
        real(dp), dimension(:), intent(inout) :: y
        type(parameters), intent(inout) :: p
    end subroutine abstract_apply

    function abstract_haspassed(this, told, t, yold, y, p) result(res)
        import :: detector, dp, parameters
        class(detector), intent(in) :: this
        real(dp), intent(in) :: told
        real(dp), intent(in) :: t
        real(dp), dimension(:), intent(in) :: yold
        real(dp), dimension(:), intent(in) :: y
        type(parameters), intent(in) :: p
        logical :: res
    end function abstract_haspassed

    function abstract_detect(this, t, p) result(res)
        import :: detector, dp, parameters
        class(detector), intent(in) :: this
        real(dp), intent(in) :: t
        type(parameters), intent(in) :: p
    end function abstract_detect
end interface

contains

function findevent(this, xa, xb, params, xtol, rtol, max_iter, err)&
        result(root)
    class(detector), intent(in) :: this
    real(dp), intent(in) :: xa
    real(dp), intent(in) :: xb
    type(parameters), intent(in) :: params
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

    fpre = this%detect(xpre, params)
    fcur = this%detect(xcur, params)

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

        fcur = this%detect(xcur, params)

        i = i + 1
    end do
end function findevent

end module events
