!
! Copyright (c) 2016 Helge Eichhorn
!
! This Source Code Form is subject to the terms of the Mozilla Public
! License, v. 2.0. If a copy of the MPL was not distributed with this
! file, You can obtain one at http://mozilla.org/MPL/2.0/.
!
module epochs

use types, only: dp

implicit none

private

public :: epoch, epochdelta, juliandate, jd2000, &
    operator (+), operator (-), seconds_per_day, mjd2000, days, seconds, calendardate, isostring, &
    datetime, to_datetime, days_per_century, seconds_per_century, centuries
!DEC$ ATTRIBUTES DLLEXPORT :: epoch, epochdelta, juliandate, jd2000, seconds_per_day, mjd2000, days
!DEC$ ATTRIBUTES DLLEXPORT :: seconds, calendardate, isostring, datetime, to_datetime, days_per_century
!DEC$ ATTRIBUTES DLLEXPORT :: seconds_per_century, centuries

real(dp), parameter :: seconds_per_day = 86400._dp
real(dp), parameter :: days_per_century = 36525._dp
real(dp), parameter :: seconds_per_century = days_per_century * seconds_per_day
real(dp), parameter :: mjd2000 = 2451544.5_dp
real(dp), parameter :: j2000 = 2451545._dp

type datetime
    integer :: year
    integer :: month
    integer :: day
    integer :: hour
    integer :: minute
    real(dp) :: seconds
end type datetime

type epoch
    real(dp) :: jd = 0._dp
    real(dp) :: jd1 = 0._dp
end type epoch

interface epoch
    module procedure new_epoch_dummy
    module procedure new_epoch_init
    module procedure new_epoch_calendar
end interface epoch

interface epochdelta
    module procedure new_epochdelta
end interface epochdelta

type epochdelta
    real(dp) :: deltajd = 0._dp
    real(dp) :: deltajd1 = 0._dp
end type epochdelta

interface operator (+)
    module procedure apply_epochdelta
end interface

interface operator (-)
    module procedure subtract_epochs
end interface

interface isostring
    module procedure isostring_datetime
    module procedure isostring_epoch
end interface isostring

interface seconds
    module procedure seconds_epochdelta
    module procedure seconds_epoch
end interface seconds

interface days
    module procedure days_epochdelta
    module procedure days_epoch
end interface days

contains

pure function new_epoch_dummy() result(ep)
    type(epoch) :: ep

    ep = new_epoch_calendar(2000, 1, 1, 12)
end function new_epoch_dummy

pure function new_epoch_init(jd, jd1) result(ep)
    real(dp), intent(in) :: jd
    real(dp), intent(in), optional :: jd1

    type(epoch) :: ep

    ep%jd = jd
    if (present(jd1)) ep%jd1 = jd1
end function new_epoch_init

pure function new_epoch_calendar(year, month, day, hour, minute, seconds) result(ep)
    integer, intent(in) :: year
    integer, intent(in) :: month
    integer, intent(in) :: day
    integer, intent(in), optional :: hour
    integer, intent(in), optional :: minute
    real(dp), intent(in), optional :: seconds

    type(epoch) :: ep
    integer :: hour_
    integer :: minute_
    real(dp) :: seconds_
    real(dp), dimension(2) :: jd

    hour_ = 0
    if (present(hour)) hour_ = hour
    minute_ = 0
    if (present(minute)) minute_ = minute
    seconds_ = 0._dp
    if (present(seconds)) seconds_ = seconds

    jd = juliandate(year, month, day, hour_, minute_, seconds_)
    ep%jd = jd(1)
    ep%jd1 = jd(2)
end function new_epoch_calendar

pure function centuries(ep, base) result(c)
    type(epoch), intent(in) :: ep
    real(dp), intent(in), optional :: base

    real(dp) :: c
    real(dp) :: base_

    base_ = j2000
    if (present(base)) base_ = base

    c = (ep%jd + ep%jd1 - base_) / days_per_century
end function centuries

pure function days_epoch(ep, base) result(d)
    type(epoch), intent(in) :: ep
    real(dp), intent(in), optional :: base

    real(dp) :: d
    real(dp) :: base_

    base_ = j2000
    if (present(base)) base_ = base

    d = ep%jd + ep%jd1 - base_
end function days_epoch

pure function seconds_epoch(ep, base) result(s)
    type(epoch), intent(in) :: ep
    real(dp), intent(in), optional :: base

    real(dp) :: s
    real(dp) :: base_

    base_ = j2000
    if (present(base)) base_ = base

    s = (ep%jd + ep%jd1 - base_) * seconds_per_day
end function seconds_epoch

pure function jd2000(ep) result(jd)
    type(epoch), intent(in) :: ep

    real(dp) :: jd

    jd = ep%jd + ep%jd1 - mjd2000
end function jd2000

pure function new_epochdelta(days, seconds) result(epd)
    real(dp), intent(in), optional :: days
    real(dp), intent(in), optional :: seconds

    type(epochdelta) :: epd

    if (present(days)) epd%deltajd = days
    if (present(seconds)) epd%deltajd1 = seconds / seconds_per_day
end function new_epochdelta

pure function apply_epochdelta(ep, epd) result(ep1)
    type(epoch), intent(in) :: ep
    type(epochdelta), intent(in) :: epd

    type(epoch) :: ep1

    ep1%jd = ep%jd + epd%deltajd
    ep1%jd1 = ep%jd1 + epd%deltajd1
end function apply_epochdelta

pure function subtract_epochs(ep1, ep2) result(epd)
    type(epoch), intent(in) :: ep1
    type(epoch), intent(in) :: ep2

    type(epochdelta) :: epd

    epd%deltajd = ep1%jd - ep2%jd
    epd%deltajd1 = ep1%jd1 - ep2%jd1
end function subtract_epochs

pure function days_epochdelta(epd) result(d)
    type(epochdelta), intent(in) :: epd

    real(dp) :: d

    d = epd%deltajd + epd%deltajd1
end function days_epochdelta

pure function seconds_epochdelta(epd) result(s)
    type(epochdelta), intent(in) :: epd

    real(dp) :: s

    s = (epd%deltajd + epd%deltajd1) * seconds_per_day
end function seconds_epochdelta

pure function juliandate(year, month, day, hour, minute, seconds) result(jd)
    integer, intent(in) :: year
    integer, intent(in) :: month
    integer, intent(in) :: day
    integer, intent(in) :: hour
    integer, intent(in) :: minute
    real(dp), intent(in) :: seconds

    real(dp), dimension(2) :: jd
    integer :: b
    integer :: y
    integer :: m

    if (month < 3) then
        y = year - 1
        m = month + 12
    else
        y = year
        m = month
    end if

    b = 2 - y / 100 + y / 400
    jd(1) = floor(365.25_dp * (y + 4716)) + floor(30.6001 * (m + 1)) + day + b - 1524.5_dp
    jd(2) = ((seconds / 60 + minute) / 60 + hour) / 24
end function juliandate

pure function calendardate(jd, jd1) result(dt)
    real(dp), intent(in) :: jd
    real(dp), intent(in), optional :: jd1

    type(datetime) :: dt
    integer :: z
    integer :: a
    integer :: b
    integer :: c
    integer :: d
    integer :: e
    integer :: alpha
    real(dp) :: f
    real(dp) :: date
    real(dp) :: t
    real(dp) :: hmod
    real(dp) :: jd1_

    jd1_ = 0._dp
    if (present(jd1)) jd1_ = jd1

    date = jd + jd1_ + 0.5_dp

    z = int(date)
    f = date - z
    if (z < 2299161) then
        a = z
    else
        alpha = int((z - 1867216.25_dp) / 36524.25_dp)
        a = z + 1 + alpha - alpha / 4
    end if
    b = a + 1524
    c = int((b - 122.1) / 365.25)
    d = int(365.25 * c)
    e = int((b - d) / 30.6001)
    dt%day = int(b - d - int(30.6001 * e) + f)
    if (e < 14) then
        dt%month = e - 1
    else
        dt%month = e - 13
    end if
    if (dt%month > 2) then
        dt%year = c - 4716
    else
        dt%year = c - 4715
    end if
    t = f * seconds_per_day
    dt%hour = int(t / 3600._dp)
    hmod = mod(t, 3600._dp)
    dt%seconds = mod(hmod, 60._dp)
    dt%minute = int((hmod - dt%seconds) / 60._dp)
end function calendardate

pure function isostring_datetime(dt) result(str)
    type(datetime), intent(in) :: dt

    character(len=23) :: str
    character(len=57) :: fmt

    fmt =  "(i4.4, a, i2.2, a, i2.2, a, i2.2, a, i2.2, a, i2.2, f0.3)"
    write(str, fmt) dt%year, "-", dt%month , "-", dt%day, "T", dt%hour, ":", dt%minute, ":",&
        int(dt%seconds), dt%seconds - int(dt%seconds)
end function isostring_datetime

pure function isostring_epoch(ep) result(str)
    type(epoch), intent(in) :: ep

    character(len=23) :: str

    str = isostring_datetime(to_datetime(ep))
end function isostring_epoch

pure function to_datetime(ep) result(dt)
    type(epoch), intent(in) :: ep

    type(datetime) :: dt

    dt = calendardate(ep%jd, ep%jd1)
end function to_datetime

end module epochs
