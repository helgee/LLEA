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
module propagators

use, intrinsic :: iso_c_binding, only: c_ptr, c_f_pointer, c_loc
use bodies, only: body
use constants, only: earth
use containers, only: parameters
use dopri, only: densestate
use epochs, only: epochdelta, seconds, operator (+)
use events, only: findevent, event
use exceptions
use forces, only: model, gravity, drag, thirdbody, uniformgravity
use integrators, only: integrate
use math, only: eps, isapprox, norm, pih, cross, cot, linspace, findroot, &
    isin
use states, only: state, framelen
use trajectories, only: trajectory, add_node, save_trajectory
use types, only: dp

implicit none

type, abstract :: propagator
contains
    procedure(abstract_trajectory), deferred :: trajectory
    procedure(abstract_state), deferred :: state
end type propagator

abstract interface
    function abstract_trajectory(p, s0, epd, err) result(tra)
        import :: propagator, state, epochdelta, exception, trajectory
        class(propagator), intent(inout), target :: p
        type(state), intent(in) :: s0
        type(epochdelta), intent(in) :: epd
        type(exception), intent(inout), optional :: err
        type(trajectory) :: tra
    end function abstract_trajectory

    function abstract_state(p, s0, epd, err) result(s)
        import :: propagator, state, epochdelta, exception, trajectory
        class(propagator), intent(inout), target :: p
        type(state), intent(in) :: s0
        type(epochdelta), intent(in) :: epd
        type(exception), intent(inout), optional :: err
        type(state) :: s
    end function abstract_state
end interface

type, extends(propagator) :: kepler
    integer :: iterations = 50
    integer :: points = 100
    real(dp) :: rtol = sqrt(eps)
contains
    procedure :: trajectory => kepler_trajectory
    procedure :: state => kepler_state
end type kepler

interface kepler
    module procedure kepler_init
end interface kepler

type, extends(propagator) :: ode
    character(len=framelen) :: frame = "GCRF"
    character(len=6) :: integrator = "dop853"
    type(body) :: center
    real(dp) :: maxstep = 0._dp
    integer :: numstep = 100000
    class(gravity), allocatable :: gravity
    class(drag), allocatable :: drag
    class(model), allocatable :: thirdbody
    logical :: savetrajectory = .false.
    type(event), dimension(:), allocatable :: events
    type(parameters) :: parameters
contains
    procedure :: trajectory => ode_trajectory
    procedure :: state => ode_state
end type ode

interface ode
    module procedure :: ode_init
end interface ode

private

public :: kepler, solve_kepler, getstate, gettrajectory, ode

contains

function ode_init(frame, integrator, center, maxstep, numstep, gravmodel, dragmodel, tbmodel, events) result(p)
    character(len=*), intent(in), optional :: frame
    character(len=*), intent(in), optional :: integrator
    type(body), intent(in), optional :: center
    real(dp), intent(in), optional :: maxstep
    integer, intent(in), optional :: numstep
    class(gravity), intent(in), optional :: gravmodel
    class(drag), intent(in), optional :: dragmodel
    class(model), intent(in), optional :: tbmodel
    type(event), dimension(:), intent(in), optional :: events
    type(ode) :: p

    type(uniformgravity) :: grav

    if (present(frame)) p%frame = frame
    if (present(integrator)) p%integrator = integrator
    p%center = earth
    if (present(center)) p%center = center
    if (present(maxstep)) p%maxstep = maxstep
    if (present(numstep)) p%numstep = numstep
    grav = uniformgravity()
    allocate(p%gravity, source=grav)
    if (present(gravmodel)) allocate(p%gravity, source=gravmodel)
    if (present(dragmodel)) allocate(p%drag, source=dragmodel)
    if (present(tbmodel)) allocate(p%thirdbody, source=tbmodel)
    if (present(events)) p%events = events
end function ode_init

subroutine rhs(n, t, y, f, tnk)
    integer, intent(in) :: n
    real(dp), intent(in) :: t
    real(dp), dimension(n), intent(in) :: y
    real(dp), dimension(n), intent(out) :: f
    type(c_ptr), intent(in) :: tnk

    type(ode), pointer :: p
    call c_f_pointer(tnk, p)

    f = 0._dp
    call p%gravity%update(f, t, y, p%parameters)
    if (allocated(p%thirdbody)) then
        call p%thirdbody%update(f, t, y, p%parameters)
    end if
    if (allocated(p%drag)) then
        call p%drag%update(f, t, y, p%parameters)
    end if
end subroutine rhs

subroutine callback(nr, told, t, y, n, con, icomp,&
                    nd, tnk, irtrn, xout)
    integer, intent(in) :: n
    integer, intent(in) :: nr
    integer, intent(in) :: nd
    integer, intent(inout) :: irtrn
    integer, dimension(nd), intent(in) :: icomp
    real(dp), intent(in) :: told
    real(dp), intent(inout) :: t
    real(dp), dimension(n), intent(inout) :: y
    real(dp), dimension(8*nd), intent(in) :: con
    real(dp), intent(inout) :: xout
    type(c_ptr), intent(in) :: tnk

    type(ode), pointer :: p
    type(parameters), pointer :: params
    type(event), pointer :: evt
    integer :: i
    logical :: firststep
    real(dp), dimension(n) :: yevt
    type(exception) :: err_
    real(dp) :: tevt
    logical :: passed

    ! Flag that signals the integrator to abort the integration
    integer, parameter :: abort = -1
    ! Flag that signals an altered numerical solution do the integrator
    integer, parameter :: altered = 2

    firststep = isapprox(told, t)
    xout = 0._dp
    call c_f_pointer(tnk, p)
    params => p%parameters
    if (allocated(p%events)) then
        p%parameters%con = con
        p%parameters%icomp = icomp
        p%parameters%nd = nd
        do i = 1, size(p%events)
            evt => p%events(i)
            ! Skip discontinuities that have already been applied
            if (evt%done) cycle

            passed = .false.

            if (allocated(evt%det) .and. .not.firststep) then
                ! Handle events with a detector
                passed = evt%det%haspassed(t, told, p%parameters)
                if (passed) then
                    tevt = evt%det%find(t, told, p%parameters)
                    evt%t = tevt
                end if
            else if (allocated(evt%t)) then
                ! Handle events with a fixed time
                if (isapprox(t, evt%t)) then
                    passed = .true.
                    tevt = t
                else if (told < evt%t .and. t > evt%t) then
                    passed = .true.
                    tevt = evt%t
                end if
            end if

            if (.not.passed) cycle

            if (passed .and. evt%abort) then
                y = densestate(tevt, n, con, icomp)
                t = tevt
                irtrn = abort
                return
            end if

            if (allocated(evt%up)) then
                ! Activate the state vector updater for discontinuities
                evt%done = .true.
                yevt = densestate(tevt, n, con, icomp)
                call evt%up%apply(yevt, p%parameters)
                y = yevt
                t = tevt
                irtrn = altered
            else
                ! Only log the time of normal events
                if (allocated(evt%tlog) .and. .not.isin(tevt, evt%tlog)) then
                    evt%tlog = [evt%tlog, tevt]
                else
                    evt%tlog = [tevt]
                end if
            end if
        end do
    end if

    if (allocated(params%trajectory)) then
        call add_node(params%trajectory, t, y)
    end if
end subroutine callback

function getstate(s0, dt, p, err) result(s)
    type(state), intent(in) :: s0
    class(*), intent(in) :: dt
    class(propagator), intent(inout) :: p
    type(exception), intent(out), optional :: err
    type(state) :: s

    type(exception) :: err_

    select type (p)
    class is (propagator)
        select type (dt)
        type is (real(dp))
            s = p%state(s0, epochdelta(seconds=dt), err_)
        type is (epochdelta)
            s = p%state(s0, dt, err_)
        end select
    end select
    if (iserror(err_)) then
        call catch(err_, "getstate", __FILE__, __LINE__)
        if (present(err)) then
            err = err_
            return
        else
            call raise(err_)
        end if
    end if
end function getstate

function gettrajectory(s0, dt, p, err) result(tra)
    type(state), intent(in) :: s0
    class(*), intent(in) :: dt
    class(propagator), intent(inout) :: p
    type(exception), intent(out), optional :: err
    type(trajectory) :: tra

    type(exception) :: err_

    select type (p)
    class is (propagator)
        select type (dt)
        type is (real(dp))
            tra = p%trajectory(s0, epochdelta(seconds=dt), err_)
        type is (epochdelta)
            tra = p%trajectory(s0, dt, err_)
        end select
    end select
    if (iserror(err_)) then
        call catch(err_, "gettrajectory", __FILE__, __LINE__)
        if (present(err)) then
            err = err_
            return
        else
            call raise(err_)
        end if
    end if
end function gettrajectory

function ode_trajectory(p, s0, epd, err) result(tra)
    class(ode), intent(inout), target :: p
    type(state), intent(in) :: s0
    type(epochdelta), intent(in) :: epd
    type(exception), intent(inout), optional :: err
    type(trajectory) :: tra

    type(exception) :: err_
    type(ode), pointer :: p_
    type(c_ptr) :: tnk
    real(dp), dimension(6) :: rv
    real(dp) :: t
    real(dp) :: tend

    p%parameters = parameters(s0, p%frame, p%center)
    allocate(p%parameters%trajectory, source=trajectory(s0))
    p_ => null()
    ! Gfortran workaround
    select type (p)
    type is (ode)
        p_ => p
    end select

    tnk = c_loc(p_)
    rv = s0%rv
    t = 0._dp
    tend = seconds(epd)
    call integrate(p%integrator, rhs, rv, t, tend, tnk, maxstep=p%maxstep, err=err_, &
        solout=callback)
    if (iserror(err_)) then
        call catch(err_, "ode_trajectory", __FILE__, __LINE__)
        if (present(err)) then
            err = err_
            return
        else
            call raise(err_)
        end if
    end if
    call save_trajectory(p%parameters%trajectory)
    tra = p%parameters%trajectory
end function ode_trajectory

function ode_state(p, s0, epd, err) result(s1)
    class(ode), intent(inout), target :: p
    type(state), intent(in) :: s0
    type(epochdelta), intent(in) :: epd
    type(exception), intent(inout), optional :: err
    type(state) :: s1

    type(exception) :: err_
    type(ode), pointer :: p_
    type(c_ptr) :: tnk
    real(dp), dimension(6) :: rv
    real(dp) :: t
    real(dp) :: tend

    p%parameters = parameters(s0, p%frame, p%center)
    p_ => null()
    ! Gfortran workaround
    select type (p)
    type is (ode)
        p_ => p
    end select

    tnk = c_loc(p_)
    rv = s0%rv
    t = 0._dp
    tend = seconds(epd)
    call integrate(p%integrator, rhs, rv, t, tend, tnk, maxstep=p%maxstep, err=err_, &
        solout=callback)
    if (iserror(err_)) then
        call catch(err_, "ode_state", __FILE__, __LINE__)
        if (present(err)) then
            err = err_
            return
        else
            call raise(err_)
        end if
    end if
    s1 = state(s0%ep + epochdelta(seconds=t), rv, s0%frame, s0%center)
end function ode_state

function kepler_trajectory(p, s0, epd, err) result(tra)
    class(kepler), intent(inout), target :: p
    type(state), intent(in) :: s0
    type(epochdelta), intent(in) :: epd
    type(exception), intent(inout), optional :: err
    type(trajectory) :: tra

    integer :: i
    real(dp), dimension(p%points) :: times
    real(dp), dimension(:,:), allocatable :: vectors
    type(exception) :: err_

    times = linspace(0._dp, seconds(epd), p%points)
    allocate(vectors(6,p%points))
    do i = 1, p%points
        vectors(:,i) = solve_kepler(s0%center%mu, s0%rv, times(i), p%iterations, &
            p%rtol, err_)
        if (iserror(err_)) then
            call catch(err_, "kepler_trajectory", __FILE__, __LINE__)
            if (present(err)) then
                err = err_
                return
            else
                call raise(err_)
            end if
        end if
    end do
    tra = trajectory(s0, times, vectors)
end function kepler_trajectory

function kepler_state(p, s0, epd, err) result(s1)
    class(kepler), intent(inout), target :: p
    type(state), intent(in) :: s0
    type(epochdelta), intent(in) :: epd
    type(exception), intent(inout), optional :: err

    type(state) :: s1
    type(exception) :: err_
    real(dp), dimension(6) :: rv1

    rv1 = solve_kepler(s0%center%mu, s0%rv, seconds(epd), p%iterations, p%rtol, err_)
    if (iserror(err_)) then
        call catch(err_, "kepler_state", __FILE__, __LINE__)
        if (present(err)) then
            err = err_
            return
        else
            call raise(err_)
        end if
    end if
    s1 = state(s0%ep + epd, rv1, s0%frame, s0%center)
end function kepler_state

function kepler_init(iterations, points, rtol) result(kep)
    integer, intent(in), optional :: iterations
    integer, intent(in), optional :: points
    real(dp), intent(in), optional :: rtol

    type(kepler) :: kep

    if (present(iterations)) kep%iterations = iterations
    if (present(points)) kep%points = points
    if (present(rtol)) kep%rtol = rtol
end function kepler_init

function solve_kepler(mu, rv, dt, iterations, rtol, err) result(rv1)
    real(dp), intent(in) :: mu
    real(dp), dimension(:), intent(in) :: rv
    real(dp), intent(in) :: dt
    integer, intent(in), optional :: iterations
    real(dp), intent(in), optional :: rtol
    type(exception), intent(inout), optional :: err

    real(dp), dimension(6) :: rv1
    type(exception) :: err_
    integer :: iterations_
    real(dp) :: rtol_
    real(dp) :: rm
    real(dp) :: rdotv
    real(dp) :: alpha
    real(dp) :: xi1
    real(dp) :: xi2
    real(dp), dimension(3) :: h
    real(dp) :: hm
    real(dp) :: p
    real(dp) :: s
    real(dp) :: w
    real(dp) :: a
    integer :: counter
    logical :: converged
    real(dp) :: xi
    real(dp) :: c2
    real(dp) :: c3
    real(dp) :: r
    real(dp) :: psi
    real(dp) :: deltat
    real(dp) :: f
    real(dp) :: g
    real(dp) :: fdot
    real(dp) :: gdot

    iterations_ = 50
    if (present(iterations)) iterations_ = iterations
    rtol_ = sqrt(eps)
    if (present(rtol)) rtol_ = rtol

    if (abs(dt) < rtol_) then
        rv1 = rv
        return
    end if

    rm = norm(rv(1:3))
    rdotv = dot_product(rv(1:3), rv(4:6))
    alpha = -dot_product(rv(4:6), rv(4:6)) / mu + 2._dp / rm

    ! Elliptic orbit
    if (alpha > rtol_) then
        if (isapprox(alpha, 1._dp)) then
            xi1 = sqrt(mu) * dt * alpha * 0.97_dp
        else
            xi1 = sqrt(mu) * dt * alpha
        end if
    ! Parabolic orbit
    else if (abs(alpha) < rtol_) then
        h = cross(rv(1:3), rv(4:6))
        hm = norm(h)
        p = hm * hm / mu
        s = 0.5_dp * (pih - atan(3._dp * sqrt(mu / p**3) * dt))
        w = atan(tan(s)**(1._dp/3._dp))
        xi1 = sqrt(p) * (2._dp * cot(2._dp * w))
        alpha = 0._dp
    ! Hyperbolic orbit
    else
        a = 1._dp / alpha
        xi1 = sign(1._dp, dt) * sqrt(-a) * log(-2._dp * mu * dt &
            / (a * (rdotv + sign(1._dp, dt) * sqrt(-mu * alpha) * (1._dp - rm * alpha))))
    end if

    counter = 0
    converged = .false.
    xi = 0._dp
    xi2 = 0._dp
    c2 = 0._dp
    c3 = 0._dp
    r = 0._dp
    psi = 0._dp

    do while (counter < iterations_)
        counter = counter + 1
        xi = xi1
        xi2 = xi * xi
        psi = xi2 * alpha
        c2 = stumpffc2(psi)
        c3 = stumpffc3(psi)
        r = xi2 * c2 + rdotv / sqrt(mu) * xi * (1 - psi * c3) + rm * (1 - psi * c2)
        deltat = xi2 * xi * c3 + rdotv / sqrt(mu) * xi2 * c2 + rm * xi * (1 - psi * c3)
        xi1 = xi + (dt * sqrt(mu) - deltat) / r
        if (abs(xi - xi1) < rtol_) then
            converged = .true.
            exit
        end if
    end do

    if (.not.converged) then
        err_ = error("Kepler solver did not converge.", "solve_kepler", __FILE__, __LINE__)
        if (present(err)) then
            err = err_
            return
        else
            call raise(err_)
        end if
    end if

    f = 1._dp - xi2 / rm * c2
    g = dt - xi2 * xi / sqrt(mu) * c3
    fdot = sqrt(mu) / (r * rm) * xi * (psi * c3 - 1)
    gdot = 1 - xi2 / r * c2
    rv1 = [f * rv(1:3) + g * rv(4:6), fdot * rv(1:3) + gdot * rv(4:6)]
end function solve_kepler

function stumpffc2(psi) result(res)
    real(dp), intent(in) :: psi
    real(dp) :: res

    real(dp) :: eps
    real(dp) :: delta
    integer :: k

    eps = 1._dp
    if (psi > eps) then
        res = (1._dp - cos(sqrt(psi))) / psi
    else if (psi < -eps) then
        res = (cosh(sqrt(-psi)) - 1._dp) / (-psi)
    else
        res = 1._dp / 2._dp
        delta = (-psi) / gamma(2._dp + 2._dp + 1._dp)
        k = 1
        do while (.not.isapprox(res + delta, res))
            res = res + delta
            k = k + 1
            delta = (-psi)**k / gamma(2._dp*k + 2._dp + 1._dp)
        end do
    end if
end function stumpffc2

function stumpffc3(psi) result(res)
    real(dp), intent(in) :: psi
    real(dp) :: res

    real(dp) :: eps
    real(dp) :: delta
    integer :: k

    eps = 1._dp
    if (psi > eps) then
        res = (sqrt(psi) - sin(sqrt(psi))) / (psi * sqrt(psi))
    else if (psi < -eps) then
        res = (sinh(sqrt(-psi)) - sqrt(-psi)) / (-psi * sqrt(-psi))
    else
        res = 1._dp / 6._dp
        delta = (-psi) / gamma(2._dp + 3._dp + 1._dp)
        k = 1
        do while (.not.isapprox(res + delta, res))
            res = res + delta
            k = k + 1
            delta = (-psi)**k / gamma(2._dp*k + 3._dp + 1._dp)
        end do
    end if
end function stumpffc3

end module propagators
