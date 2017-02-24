!
! Copyright (c) 2016 Helge Eichhorn
!
! This Source Code Form is subject to the terms of the Mozilla Public
! License, v. 2.0. If a copy of the MPL was not distributed with this
! file, You can obtain one at http://mozilla.org/MPL/2.0/.
!
module integrators

use, intrinsic :: iso_c_binding, only: c_ptr, c_null_ptr
use dopri, only: dop853, dopri5, contd8, contd5, soldummy, dopmessage
use exceptions
use math, only: eps
use types, only: dp

implicit none

private

public :: integrate

contains

! Subroutine: integrate
!   Adaptive step-size integrator with dense output.
!
! Parameters:
!   fcn - User subroutine: fcn(n, x, y, f, rpar, ipar)
!   y - Initial state (on entry) and final state (upon return)
!   t - Initial time (on entry) and final time (upon return)
!   tend - Stop time
!   rpar - Real-valued parameters
!   ipar - Integer-valued parameters
!   maxstep - Maximum step size (optional)
!   stat - Exit code (optional)
!   rtol - Relative tolerance (optional, default: 1e-6)
!   atol - Absolute tolerance (optional, default: 1.48e-8)
subroutine integrate(integrator, fcn, y, t, tend, tnk,&
        solout, maxstep, nstiff, nsteps, rtol, atol, err)
    character(len=*), intent(in) :: integrator
    real(dp), dimension(:), intent(inout) :: y
    real(dp), intent(inout) :: t
    real(dp), intent(in) :: tend
    type(c_ptr), intent(in), optional :: tnk
    real(dp), intent(in), optional :: maxstep
    integer, intent(in), optional :: nstiff
    integer, intent(in), optional :: nsteps
    real(dp), dimension(:), intent(in), optional :: rtol
    real(dp), dimension(:), intent(in), optional :: atol
    type(exception), intent(inout), optional :: err
    interface
        subroutine fcn(n, t, y, f, tnk)
            import :: dp, c_ptr
            integer, intent(in) :: n
            real(dp), intent(in) :: t
            real(dp), dimension(n), intent(in) :: y
            real(dp), dimension(n), intent(out) :: f
            type(c_ptr), intent(in) :: tnk
        end subroutine fcn
        subroutine solout(nr, told, t, y, n, con, icomp,&
                            nd, tnk, irtrn, tout)
            import :: dp, c_ptr
            integer, intent(in) :: n
            integer, intent(in) :: nr
            integer, intent(in) :: nd
            integer, intent(inout) :: irtrn
            integer, dimension(nd), intent(in) :: icomp
            real(dp), intent(in) :: told
            real(dp), intent(inout) :: t
            real(dp), dimension(n), intent(inout) :: y
            real(dp), dimension(8*nd), intent(in) :: con
            real(dp), intent(inout) :: tout
            type(c_ptr), intent(in) :: tnk
        end subroutine solout
    end interface
    optional :: solout

    integer :: lwork
    integer :: liwork
    integer :: iout
    integer :: idid
    integer :: itol
    integer, dimension(:), allocatable :: iwork
    integer :: n
    real(dp), dimension(:), allocatable :: rtol_
    real(dp), dimension(:), allocatable :: atol_
    real(dp), dimension(:), allocatable :: work
    real(dp), dimension(:), allocatable :: y0
    type(c_ptr) :: tnk_
    type(exception) :: err_

    tnk_ = c_null_ptr
    if (present(tnk)) tnk_ = tnk

    select case (integrator)
    case ("dop853", "dopri5")
        n = size(y)
        lwork = 11*n + 8*n + 21
        liwork = n + 21
        allocate(work(lwork))
        allocate(iwork(liwork))
        allocate(y0(n))
        allocate(rtol_(n))
        allocate(atol_(n))

        if (present(rtol)) then
            rtol_ = rtol
        else
            rtol_ = 1e-6_dp
        end if
        if (present(atol)) then
            atol_ = atol
        else
            atol_ = sqrt(eps)
        end if

        ! Perform dense output after every successful step.
        iout = 2
        ! Tolerances are given in vector form.
        itol = 1
        ! Set all parameters to default values.
        iwork = 0
        work = 0._dp

        if (present(maxstep)) work(6) = maxstep
        if (present(nstiff)) iwork(4) = nstiff
        if (present(nsteps)) iwork(1) = nsteps

        ! Do not print verbose error messages.
        iwork(3) = -1
        ! Number of components for dense output.
        iwork(5) = n

        y0 = y
        if (present(solout)) then
            call dop853(n, fcn, t, y0, tend, rtol_, atol_,&
                itol, solout, iout, work, lwork, iwork,&
                liwork, tnk_, idid)
        else
            call dop853(n, fcn, t, y0, tend, rtol_, atol_,&
                itol, soldummy, iout, work, lwork, iwork,&
                liwork, tnk_, idid)
        end if
        y = y0

        if (idid < 0) then
            err_ = error(dopmessage(idid), "integrate", __FILE__, __LINE__)
            if (present(err)) then
                err = err_
                return
            else
                call raise(err_)
            end if
        end if
    end select
end subroutine integrate

end module integrators
