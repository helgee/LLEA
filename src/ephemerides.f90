!
! Copyright (c) 2016 Helge Eichhorn
!
! This Source Code Form is subject to the terms of the Mozilla Public
! License, v. 2.0. If a copy of the MPL was not distributed with this
! file, You can obtain one at http://mozilla.org/MPL/2.0/.
!
module ephemerides

use epochs, only: seconds_per_day, epoch
use exceptions
use naif, only: naifid
use types, only: dp
use util, only: joinpath, sep, projectdir

implicit none

type, abstract :: ephemeris
contains
    procedure(abstract_position), deferred :: position
    procedure(abstract_velocity), deferred :: velocity
    procedure(abstract_state), deferred :: state
end type ephemeris

abstract interface
    function abstract_position(this, ep, to, from, err) result(r)
        import :: ephemeris, epoch, dp, exception
        class(ephemeris), intent(inout) :: this
        type(epoch), intent(in) :: ep
        integer, intent(in) :: to
        integer, intent(in), optional :: from
        type(exception), intent(inout), optional :: err
        real(dp), dimension(3) :: r
    end function abstract_position

    function abstract_velocity(this, ep, to, from, err) result(v)
        import :: ephemeris, epoch, dp, exception
        class(ephemeris), intent(inout) :: this
        type(epoch), intent(in) :: ep
        integer, intent(in) :: to
        integer, intent(in), optional :: from
        type(exception), intent(inout), optional :: err
        real(dp), dimension(3) :: v
    end function abstract_velocity

    function abstract_state(this, ep, to, from, err) result(s)
        import :: ephemeris, epoch, dp, exception
        class(ephemeris), intent(inout) :: this
        type(epoch), intent(in) :: ep
        integer, intent(in) :: to
        integer, intent(in), optional :: from
        type(exception), intent(inout), optional :: err
        real(dp), dimension(6) :: s
    end function abstract_state
end interface

type daf
    character(len=1024) :: path = ""
    character(len=13) :: endianness = ""
    character(len=8) :: id = ""
    character(len=60) :: name = ""
    integer :: nd = -1
    integer :: ni = -1
    integer :: first = -1
    integer :: last = -1
    integer :: ss = -1
    integer :: nc = -1
end type daf

interface daf
    module procedure new_daf
end interface daf

type spksegment
    real(dp) :: firstsec = -1._dp
    real(dp) :: lastsec = -1._dp
    real(dp) :: firstdate = -1._dp
    real(dp) :: lastdate = -1._dp
    integer :: targ = -1
    integer :: origin = -1
    integer :: frame = -1
    integer :: spktype = -1
    integer :: firstaddr = -1
    integer :: lastaddr = -1
    real(dp) :: initialsecond = -1._dp
    real(dp) :: intervall = -1._dp
    integer :: recordsize = -1
    integer :: nrecords = -1
    real(dp), dimension(:,:), allocatable :: cache
    integer :: cachedrecord = -1
end type spksegment

type, extends(ephemeris) :: jplephem
    type(daf) :: daffile
    type(spksegment), dimension(:), allocatable :: segments
contains
    procedure :: coefficients => coefficients
    procedure :: getsegnum => getsegnum
    procedure :: position => jplephem_position
    procedure :: position_lowlevel => jplephem_position_lowlevel
    procedure :: velocity => jplephem_velocity
    procedure :: velocity_lowlevel => jplephem_velocity_lowlevel
    procedure :: state => jplephem_state
    procedure :: state_lowlevel => jplephem_state_lowlevel
end type jplephem

interface jplephem
    module procedure new_jplephem
end interface jplephem

interface getposition
    module procedure getposition_jplephem_name
    module procedure getposition_jplephem_id
end interface getposition

interface getvelocity
    module procedure getvelocity_jplephem_name
    module procedure getvelocity_jplephem_id
end interface getvelocity

interface getstate
    module procedure getstate_jplephem_name
    module procedure getstate_jplephem_id
end interface getstate

class(ephemeris), allocatable :: ephem

private

public :: daf, jplephem, ephemeris, ephem, naifid, init_ephemeris, getpath, &
    getposition, getvelocity, getstate

contains

function getposition_jplephem_name(eph, ep, to, from, err) result(r)
    type(jplephem), intent(inout) :: eph
    type(epoch), intent(in) :: ep
    character(len=*), intent(in) :: to
    character(len=*), intent(in), optional :: from
    type(exception), intent(inout), optional :: err
    real(dp), dimension(3) :: r

    type(exception) :: err_
    integer :: to_
    integer :: from_

    to_ = naifid(to, err_)
    if (iserror(err_)) then
        call catch(err_, "getposition_jplephem_name", __FILE__, __LINE__)
        if (present(err)) then
            err = err_
            return
        else
            call raise(err_)
        end if
    end if
    from_ = 0
    if (present(from)) then
        from_ = naifid(from, err_)
        if (iserror(err_)) then
            call catch(err_, "getposition_jplephem_name", __FILE__, __LINE__)
            if (present(err)) then
                err = err_
                return
            else
                call raise(err_)
            end if
        end if
    end if
    r = eph%position(ep, to_, from_, err_)
    if (iserror(err_)) then
        call catch(err_, "getposition_jplephem_name", __FILE__, __LINE__)
        if (present(err)) then
            err = err_
            return
        else
            call raise(err_)
        end if
    end if
end function getposition_jplephem_name

function getposition_jplephem_id(eph, ep, to, from, err) result(r)
    type(jplephem), intent(inout) :: eph
    type(epoch), intent(in) :: ep
    integer, intent(in) :: to
    integer, intent(in), optional :: from
    type(exception), intent(inout), optional :: err
    real(dp), dimension(3) :: r

    type(exception) :: err_
    integer :: from_

    from_ = 0
    if (present(from)) from_ = from
    r = eph%position(ep, to, from_, err_)
    if (iserror(err_)) then
        call catch(err_, "getposition_jplephem_id", __FILE__, __LINE__)
        if (present(err)) then
            err = err_
            return
        else
            call raise(err_)
        end if
    end if
end function getposition_jplephem_id

function getvelocity_jplephem_name(eph, ep, to, from, err) result(r)
    type(jplephem), intent(inout) :: eph
    type(epoch), intent(in) :: ep
    character(len=*), intent(in) :: to
    character(len=*), intent(in), optional :: from
    type(exception), intent(inout), optional :: err
    real(dp), dimension(3) :: r

    type(exception) :: err_
    integer :: to_
    integer :: from_

    to_ = naifid(to, err_)
    if (iserror(err_)) then
        call catch(err_, "getvelocity_jplephem_name", __FILE__, __LINE__)
        if (present(err)) then
            err = err_
            return
        else
            call raise(err_)
        end if
    end if
    from_ = 0
    if (present(from)) then
        from_ = naifid(from, err_)
        if (iserror(err_)) then
            call catch(err_, "getvelocity_jplephem_name", __FILE__, __LINE__)
            if (present(err)) then
                err = err_
                return
            else
                call raise(err_)
            end if
        end if
    end if
    r = eph%velocity(ep, to_, from_, err_)
    if (iserror(err_)) then
        call catch(err_, "getvelocity_jplephem_name", __FILE__, __LINE__)
        if (present(err)) then
            err = err_
            return
        else
            call raise(err_)
        end if
    end if
end function getvelocity_jplephem_name

function getvelocity_jplephem_id(eph, ep, to, from, err) result(r)
    type(jplephem), intent(inout) :: eph
    type(epoch), intent(in) :: ep
    integer, intent(in) :: to
    integer, intent(in), optional :: from
    type(exception), intent(inout), optional :: err
    real(dp), dimension(3) :: r

    type(exception) :: err_
    integer :: from_

    from_ = 0
    if (present(from)) from_ = from
    r = eph%velocity(ep, to, from_, err_)
    if (iserror(err_)) then
        call catch(err_, "getvelocity_jplephem_id", __FILE__, __LINE__)
        if (present(err)) then
            err = err_
            return
        else
            call raise(err_)
        end if
    end if
end function getvelocity_jplephem_id

function getstate_jplephem_name(eph, ep, to, from, err) result(r)
    type(jplephem), intent(inout) :: eph
    type(epoch), intent(in) :: ep
    character(len=*), intent(in) :: to
    character(len=*), intent(in), optional :: from
    type(exception), intent(inout), optional :: err
    real(dp), dimension(6) :: r

    type(exception) :: err_
    integer :: to_
    integer :: from_

    to_ = naifid(to, err_)
    if (iserror(err_)) then
        call catch(err_, "getstate_jplephem_name", __FILE__, __LINE__)
        if (present(err)) then
            err = err_
            return
        else
            call raise(err_)
        end if
    end if
    from_ = 0
    if (present(from)) then
        from_ = naifid(from, err_)
        if (iserror(err_)) then
            call catch(err_, "getstate_jplephem_name", __FILE__, __LINE__)
            if (present(err)) then
                err = err_
                return
            else
                call raise(err_)
            end if
        end if
    end if
    r = eph%state(ep, to_, from_, err_)
    if (iserror(err_)) then
        call catch(err_, "getstate_jplephem_name", __FILE__, __LINE__)
        if (present(err)) then
            err = err_
            return
        else
            call raise(err_)
        end if
    end if
end function getstate_jplephem_name

function getstate_jplephem_id(eph, ep, to, from, err) result(r)
    type(jplephem), intent(inout) :: eph
    type(epoch), intent(in) :: ep
    integer, intent(in) :: to
    integer, intent(in), optional :: from
    type(exception), intent(inout), optional :: err
    real(dp), dimension(6) :: r

    type(exception) :: err_
    integer :: from_

    from_ = 0
    if (present(from)) from_ = from
    r = eph%state(ep, to, from_, err_)
    if (iserror(err_)) then
        call catch(err_, "getstate_jplephem_id", __FILE__, __LINE__)
        if (present(err)) then
            err = err_
            return
        else
            call raise(err_)
        end if
    end if
end function getstate_jplephem_id

subroutine init_ephemeris(denum)
    character(len=3), intent(in), optional :: denum

    character(len=3) :: denum_

    denum_ = "430"
    if (present(denum)) denum_ = denum
    allocate(ephem, source=jplephem(joinpath(projectdir("data"), "de"//denum//".bsp")))
end subroutine init_ephemeris

function new_jplephem(path, err) result(eph)
    character(len=*), intent(in) :: path
    type(exception), intent(inout), optional :: err
    type(jplephem), target :: eph

    type(exception) :: err_
    integer :: i
    integer :: u
    integer :: nsegments
    integer :: nsum
    integer :: next
    integer :: prev
    integer :: pos
    real(dp) :: nsum_dp
    real(dp) :: next_dp
    real(dp) :: prev_dp
    real(dp) :: recordsize
    real(dp) :: nrecords
    type(spksegment), pointer :: seg

    seg => null()

    eph%daffile = daf(path, err_)
    if (iserror(err_)) then
        call catch(err_, "new_spk", __FILE__, __LINE__)
        if (present(err)) then
            err = err_
            return
        else
            call raise(err_)
        end if
    end if

    nsegments = 0
    next = eph%daffile%first
    do while (next /= 0)
        open(newunit=u, file=eph%daffile%path, status="old", form="unformatted", access="stream", &
            action="read", convert=eph%daffile%endianness)
        pos = 1024 * (next - 1) + 1
        read(u,pos=pos) next_dp, prev_dp, nsum_dp
        close(u)
        next = int(next_dp)
        prev = int(prev_dp)
        nsum = int(nsum_dp)
        nsegments = nsegments + nsum
    end do
    allocate(eph%segments(nsegments))

    nsegments = 0
    next = eph%daffile%first
    do while (next /= 0)
        open(newunit=u, file=eph%daffile%path, status="old", form="unformatted", access="stream", &
            action="read", convert=eph%daffile%endianness)
        pos = 1024 * (next - 1) + 1
        read(u,pos=pos) next_dp, prev_dp, nsum_dp
        close(u)
        nsum = int(nsum_dp)
        open(newunit=u, file=eph%daffile%path, status="old", form="unformatted", access="stream", &
            action="read", convert=eph%daffile%endianness)
        do i = 1, nsum
            nsegments = nsegments + 1
            seg => eph%segments(nsegments)
            pos = 1024 * (next - 1) + 25 + dp * eph%daffile%ss * (i-1)
            read(u,pos=pos) seg%firstsec, seg%lastsec, seg%targ, seg%origin, seg%frame, seg%spktype, &
                seg%firstaddr, seg%lastaddr
            seg%firstdate = jd(seg%firstsec)
            seg%lastdate = jd(seg%lastsec)
            pos = seg%lastaddr * dp - 4 * dp + 1
            read(u,pos=pos) seg%initialsecond, seg%intervall, recordsize, nrecords
            seg%recordsize = int(recordsize)
            seg%nrecords = int(nrecords)
        end do
        close(u)
        next = int(next_dp)
        prev = int(prev_dp)
    end do
end function new_jplephem

subroutine coefficients(eph, segnum, tdb, tdb2, x, twotc, err)
    class(jplephem), intent(inout), target :: eph
    integer, intent(in) :: segnum
    real(dp), intent(in) :: tdb
    real(dp), intent(in) :: tdb2
    real(dp), dimension(:), intent(out), allocatable :: x
    real(dp), intent(out) :: twotc
    type(exception), intent(inout), optional :: err

    type(exception) :: err_
    type(spksegment), pointer :: seg
    integer :: components
    real(dp) :: secs
    integer :: recordnum
    real(dp) :: frac
    real(dp) :: tc
    integer :: order
    integer :: pos
    integer :: i
    integer :: u
    real(dp), dimension(:), allocatable :: c

    seg => eph%segments(segnum)
    if ((tdb + tdb2 < seg%firstdate).or.(tdb + tdb2 > seg%lastdate)) then
        err_ = error("Date is out of range of the ephemeris.", "coefficients", __FILE__, __LINE__, &
            id="OutOfRangeError")
        if (present(err)) then
            err = err_
            return
        else
            call raise(err_)
        end if
    end if

    components = 3
    order = (seg%recordsize - 2) / 3
    allocate(c(order*components))
    secs = (seconds(tdb) - seg%initialsecond) + tdb2 * seconds_per_day
    recordnum = floor(secs/seg%intervall)
    frac = mod(secs, seg%intervall)
    if (recordnum == seg%nrecords) then
        recordnum = recordnum - 1
        frac = seg%intervall
    end if
    if (recordnum /= seg%cachedrecord) then
        pos = dp * (seg%firstaddr + seg%recordsize * recordnum + 1) + 1
        open(newunit=u, file=eph%daffile%path, status="old", form="unformatted", access="stream", &
            action="read", convert=eph%daffile%endianness)
        read(u,pos=pos) c
        close(u)
        if (allocated(seg%cache)) deallocate(seg%cache)
        allocate(seg%cache(components, order))
        seg%cache = transpose(reshape(c, [order, components]))
        seg%cachedrecord = recordnum
    end if
    allocate(x(order))
    tc = 2._dp * frac / seg%intervall - 1._dp
    x(1) = 1._dp
    x(2) = tc
    twotc = tc + tc
    do i = 3, order
        x(i) = twotc * x(i-1) - x(i-2)
    end do
end subroutine coefficients

recursive subroutine getpath(path, to)
    integer, dimension(:), intent(inout), allocatable :: path
    integer, intent(in) :: to

    integer :: current

    current = path(size(path))
    if (current == to) then
        return
    elseif (isbody(current)) then
        path = [path, current/100]
        call getpath(path, to)
    elseif (isbarycenter(current).and.ischild(current, to)) then
        path = [path, to]
        call getpath(path, to)
    elseif (isbarycenter(current).and..not.ischild(current, to)) then
        path = [path, 0]
        call getpath(path, to)
    elseif (current == 0.and.isbarycenter(to)) then
        path = [path, to]
        call getpath(path, to)
    elseif (current == 0.and.isbody(to)) then
        path = [path, to/100]
        call getpath(path, to)
    end if
end subroutine getpath

function isbarycenter(id) result(res)
    integer, intent(in) :: id
    logical :: res
    res = id < 100 .and. id > 0
end function isbarycenter

function isbody(id) result(res)
    integer, intent(in) :: id
    logical :: res
    res = id > 100
end function isbody

function ischild(parent, id) result(res)
    integer, intent(in) :: parent
    integer, intent(in) :: id
    logical :: res
    res = parent == id / 100
end function ischild

function jplephem_position(this, ep, to, from, err) result(r)
    class(jplephem), intent(inout) :: this
    type(epoch), intent(in) :: ep
    integer, intent(in) :: to
    integer, intent(in), optional :: from
    type(exception), intent(inout), optional :: err
    real(dp), dimension(3) :: r

    type(exception) :: err_
    integer :: from_
    integer :: i
    integer, dimension(:), allocatable :: path
    integer :: segnum

    from_ = 0
    if (present(from)) then
        from_ = from
    end if
    path = [from_]
    call getpath(path, to)

    r = 0._dp
    do i = 1, size(path) - 1
        segnum = this%getsegnum(path(i), path(i+1))
        r = r + sign(1, segnum) * &
            this%position_lowlevel(ep%jd, ep%jd1, abs(segnum), err_)
        if (iserror(err_)) then
            call catch(err_, "jplephem_position", __FILE__, __LINE__)
            if (present(err)) then
                err = err_
                return
            else
                call raise(err_)
            end if
        end if
    end do
end function jplephem_position

function jplephem_position_lowlevel(this, tdb, tdb2, segnum, err) result(r)
    class(jplephem), intent(inout) :: this
    real(dp), intent(in) :: tdb
    real(dp), intent(in) :: tdb2
    integer, intent(in) :: segnum
    type(exception), intent(inout), optional :: err
    real(dp), dimension(3) :: r

    real(dp), dimension(:), allocatable :: x
    real(dp) :: twotc
    type(exception) :: err_

    call this%coefficients(segnum, tdb, tdb2, x, twotc, err_)
    if (iserror(err_)) then
        call catch(err_, "jplephem_position_lowlevel", __FILE__, __LINE__)
        if (present(err)) then
            err = err_
            return
        else
            call raise(err_)
        end if
    end if
    r = matmul(this%segments(segnum)%cache, x)
end function jplephem_position_lowlevel

function jplephem_velocity(this, ep, to, from, err) result(v)
    class(jplephem), intent(inout) :: this
    type(epoch), intent(in) :: ep
    integer, intent(in) :: to
    integer, intent(in), optional :: from
    type(exception), intent(inout), optional :: err
    real(dp), dimension(3) :: v

    type(exception) :: err_
    integer :: from_
    integer :: i
    integer, dimension(:), allocatable :: path
    integer :: segnum

    from_ = 0
    if (present(from)) then
        from_ = from
    end if
    path = [from_]
    call getpath(path, to)

    v = 0._dp
    do i = 1, size(path) - 1
        segnum = this%getsegnum(path(i), path(i+1))
        v = v + sign(1, segnum) * &
            this%velocity_lowlevel(ep%jd, ep%jd1, abs(segnum), err_)
        if (iserror(err_)) then
            call catch(err_, "jplephem_position", __FILE__, __LINE__)
            if (present(err)) then
                err = err_
                return
            else
                call raise(err_)
            end if
        end if
    end do
end function jplephem_velocity

function jplephem_velocity_lowlevel(this, tdb, tdb2, segnum, err) result(v)
    class(jplephem), intent(inout) :: this
    real(dp), intent(in) :: tdb
    real(dp), intent(in) :: tdb2
    integer, intent(in) :: segnum
    type(exception), intent(inout), optional :: err
    real(dp), dimension(3) :: v

    real(dp), dimension(:), allocatable :: x
    real(dp) :: twotc
    real(dp), dimension(:), allocatable :: t
    integer :: i
    type(exception) :: err_

    call this%coefficients(segnum, tdb, tdb2, x, twotc, err_)
    if (iserror(err_)) then
        call catch(err_, "jplephem_velocity_lowlevel", __FILE__, __LINE__)
        if (present(err)) then
            err = err_
            return
        else
            call raise(err_)
        end if
    end if

    allocate(t(size(x)))
    t(1) = 0._dp
    t(2) = 1._dp
    if (size(t) > 2) then
        t(3) = twotc + twotc
        do i = 4, size(t)
            t(i) = twotc * t(i-1) - t(i-2) + x(i-1) + x(i-1)
        end do
    end if
    t = 2 * t / this%segments(segnum)%intervall
    v = matmul(this%segments(segnum)%cache, t)
end function jplephem_velocity_lowlevel

function jplephem_state(this, ep, to, from, err) result(s)
    class(jplephem), intent(inout) :: this
    type(epoch), intent(in) :: ep
    integer, intent(in) :: to
    integer, intent(in), optional :: from
    type(exception), intent(inout), optional :: err
    real(dp), dimension(6) :: s

    type(exception) :: err_
    integer :: from_
    integer :: i
    integer, dimension(:), allocatable :: path
    integer :: segnum

    from_ = 0
    if (present(from)) then
        from_ = from
    end if
    path = [from_]
    call getpath(path, to)

    s = 0._dp
    do i = 1, size(path) - 1
        segnum = this%getsegnum(path(i), path(i+1))
        s = s + sign(1, segnum) * &
            this%state_lowlevel(ep%jd, ep%jd1, abs(segnum), err_)
        if (iserror(err_)) then
            call catch(err_, "jplephem_position", __FILE__, __LINE__)
            if (present(err)) then
                err = err_
                return
            else
                call raise(err_)
            end if
        end if
    end do
end function jplephem_state

function jplephem_state_lowlevel(this, tdb, tdb2, segnum, err) result(s)
    class(jplephem), intent(inout) :: this
    real(dp), intent(in) :: tdb
    real(dp), intent(in) :: tdb2
    integer, intent(in) :: segnum
    type(exception), intent(inout), optional :: err
    real(dp), dimension(6) :: s

    real(dp), dimension(:), allocatable :: x
    real(dp) :: twotc
    real(dp), dimension(:), allocatable :: t
    integer :: i
    type(exception) :: err_

    call this%coefficients(segnum, tdb, tdb2, x, twotc, err_)
    if (iserror(err_)) then
        call catch(err_, "jplephem_state_lowlevel", __FILE__, __LINE__)
        if (present(err)) then
            err = err_
            return
        else
            call raise(err_)
        end if
    end if

    allocate(t(size(x)))
    t(1) = 0._dp
    t(2) = 1._dp
    if (size(t) > 2) then
        t(3) = twotc + twotc
        do i = 4, size(t)
            t(i) = twotc * t(i-1) - t(i-2) + x(i-1) + x(i-1)
        end do
    end if
    t = 2 * t / this%segments(segnum)%intervall
    s(1:3) = matmul(this%segments(segnum)%cache, x)
    s(4:6) = matmul(this%segments(segnum)%cache, t)
end function jplephem_state_lowlevel

function jd(sec) result(res)
    real(dp), intent(in) :: sec
    real(dp) :: res
    res = 2451545 + sec / seconds_per_day
end function jd

function seconds(jd) result(res)
    real(dp), intent(in) :: jd
    real(dp) :: res
    res = (jd - 2451545) * seconds_per_day
end function seconds

function getsegnum(this, origin, targ, err) result(segnum)
    class(jplephem), intent(in) :: this
    integer, intent(in) :: origin
    integer, intent(in) :: targ
    type(exception), intent(inout), optional :: err
    integer :: segnum

    integer :: i
    type(exception) :: err_
    character(len=8) :: tstr
    character(len=8) :: ostr

    segnum = -999

    do i = 1, size(this%segments)
        if (this%segments(i)%targ == targ .and. this%segments(i)%origin == origin) then
            segnum = i
            return
        elseif (this%segments(i)%targ == origin .and. this%segments(i)%origin == targ) then
            segnum = -i
            return
        end if
    end do
    write(tstr,"(i8)") targ
    write(ostr,"(i8)") origin
    err_ = error("No segment with origin="//trim(adjustl(ostr))//&
        " and targ="//trim(adjustl(tstr))//" available.", &
        "getsegnum", __FILE__, __LINE__)
    if (present(err)) then
        err = err_
        return
    else
        call raise(err_)
    end if
end function getsegnum

function new_daf(path, err) result(d)
    character(len=*), intent(in) :: path
    type(exception), intent(inout), optional :: err
    type(daf) :: d

    type(exception) :: err_
    integer :: u
    logical :: isfile

    inquire(file=path, exist=isfile)
    if (.not.isfile) then
        err_ = error("File not found: "//trim(path), "new_daf", __FILE__, __LINE__)
        if (present(err)) then
            err = err_
            return
        else
            call raise(err_)
        end if
    end if
    d%path = path

    d%endianness = check_endianness(path, err_)
    if (iserror(err_)) then
        call catch(err_, "new_daf", __FILE__, __LINE__)
        if (present(err)) then
            err = err_
            return
        else
            call raise(err_)
        end if
    end if

    open(newunit=u, file=path, status="old", form="unformatted", access="direct", recl=84, &
        action="read", convert=d%endianness)
    read(u,rec=1) d%id, d%nd, d%ni, d%name, d%first, d%last
    close(u)
    d%ss = d%nd + (d%ni + 1) / 2
    d%nc = 8 * d%ss
end function new_daf

function check_endianness(path, err) result(endianness)
    character(len=*), intent(in) :: path
    type(exception), intent(inout), optional :: err
    character(len=13) :: endianness

    type(exception) :: err_
    character(len=8) :: id
    character(len=76) :: tmp
    character(len=8) :: endian
    integer :: nd
    integer :: u

    open(newunit=u, file=trim(path), status="old", form="unformatted", access="direct", recl=96, action="read")
    read(u,rec=1) id, nd, tmp, endian
    close(u)

    if (id == "NAIF/DAF") then
        if (nd == 2) then
            endianness = "little_endian"
        else
            open(newunit=u, file=path, status="old", form="unformatted", access="direct", recl=12, &
                action="read", convert="big_endian")
            read(u,rec=1) id, nd
            close(u)
            if (nd == 2) then
                endianness = "big_endian"
            else
                err_ = error("Endianess could not be detected.", "check_endianness", __FILE__, __LINE__)
                if (present(err)) then
                    err = err_
                    return
                else
                    call raise(err_)
                end if
            end if
        end if
    else
        if (endian == "LTL-IEEE") then
            endianness = "little_endian"
        else
            endianness = "big_endian"
        end if
    end if
end function check_endianness

end module ephemerides
