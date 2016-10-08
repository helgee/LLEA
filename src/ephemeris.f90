module ephemeris

use types, only: dp
use exceptions
use epochs, only: seconds_per_day

implicit none

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

type segment
    real(dp) :: firstsec = -1._dp
    real(dp) :: lastsec = -1._dp
    real(dp) :: firstdate = -1._dp
    real(dp) :: lastdate = -1._dp
    integer :: targ = -1
    integer :: center = -1
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
end type segment

type spk
    type(daf) :: daffile
    type(segment), dimension(:), allocatable :: segments
end type spk

interface spk
    module procedure new_spk
end interface spk

interface getposition
    module procedure getposition_lowlevel
    module procedure getposition_tdb
end interface getposition

interface getvelocity
    module procedure getvelocity_lowlevel
    module procedure getvelocity_tdb
end interface getvelocity

interface getstate
    module procedure getstate_tdb
end interface getstate

private

public :: daf, spk, getposition, getvelocity, getstate

contains

function new_spk(path, err) result(s)
    character(len=*), intent(in) :: path
    type(exception), intent(inout), optional :: err
    type(spk), target :: s

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
    type(segment), pointer :: seg

    seg => null()

    s%daffile = daf(path, err_)
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
    next = s%daffile%first
    do while (next /= 0)
        open(newunit=u, file=s%daffile%path, status="old", form="unformatted", access="direct", recl=1024, &
            action="read", convert=s%daffile%endianness)
        read(u,rec=next) next_dp, prev_dp, nsum_dp
        close(u)
        next = int(next_dp)
        prev = int(prev_dp)
        nsum = int(nsum_dp)
        nsegments = nsegments + nsum
    end do
    allocate(s%segments(nsegments))

    nsegments = 0
    next = s%daffile%first
    do while (next /= 0)
        open(newunit=u, file=s%daffile%path, status="old", form="unformatted", access="direct", recl=1024, &
            action="read", convert=s%daffile%endianness)
        read(u,rec=next) next_dp, prev_dp, nsum_dp
        close(u)
        open(newunit=u, file=s%daffile%path, status="old", form="unformatted", access="stream", &
            action="read", convert=s%daffile%endianness)
        do i = 1, nsum
            nsegments = nsegments + 1
            seg => s%segments(nsegments)
            pos = 1024 * (next - 1) + 25 + dp * s%daffile%ss * (i-1)
            read(u,pos=pos) seg%firstsec, seg%lastsec, seg%targ, seg%center, seg%frame, seg%spktype, &
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
        nsum = int(nsum_dp)
    end do
end function new_spk

subroutine coefficients(s, segnum, tdb, tdb2, x, twotc, err)
    type(spk), intent(inout), target :: s
    integer, intent(in) :: segnum
    real(dp), intent(in) :: tdb
    real(dp), intent(in) :: tdb2
    real(dp), dimension(:), intent(out), allocatable :: x
    real(dp), intent(out) :: twotc
    type(exception), intent(inout), optional :: err

    type(exception) :: err_
    type(segment), pointer :: seg
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

    seg => s%segments(segnum)
    if ((tdb + tdb2 < seg%firstdate).or.(tdb + tdb2 > seg%lastdate)) then
        err_ = error("Date is out of range of the ephemeris.", "coefficients", __FILE__, __LINE__)
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
        open(newunit=u, file=s%daffile%path, status="old", form="unformatted", access="stream", &
            action="read", convert=s%daffile%endianness)
        read(u,pos=pos) c
        close(u)
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

function getposition_lowlevel(c, x) result(r)
    real(dp), dimension(:,:), intent(in) :: c
    real(dp), dimension(:), intent(in) :: x
    real(dp), dimension(3) :: r
    r = matmul(c, x)
end function getposition_lowlevel

function getposition_tdb(s, segnum, tdb, tdb2, err) result(r)
    type(spk), intent(inout) :: s
    integer, intent(in) :: segnum
    real(dp), intent(in) :: tdb
    real(dp), intent(in) :: tdb2
    type(exception), intent(inout), optional :: err
    real(dp), dimension(3) :: r

    real(dp), dimension(:), allocatable :: x
    real(dp) :: twotc
    type(exception) :: err_

    call coefficients(s, segnum, tdb, tdb2, x, twotc, err_)
    if (iserror(err_)) then
        call catch(err_, "getposition_tdb", __FILE__, __LINE__)
        if (present(err)) then
            err = err_
            return
        else
            call raise(err_)
        end if
    end if
    r = getposition_lowlevel(s%segments(segnum)%cache, x)
end function getposition_tdb

function getvelocity_lowlevel(seg, x, twotc) result(v)
    type(segment), intent(in) :: seg
    real(dp), dimension(:), intent(in) :: x
    real(dp), intent(in) :: twotc
    real(dp), dimension(3) :: v

    real(dp), dimension(size(x)) :: t
    integer :: i

    t(1) = 0._dp
    t(2) = 1._dp
    if (size(t) > 2) then
        t(3) = twotc + twotc
        do i = 4, size(t)
            t(i) = twotc * t(i-1) - t(i-2) + x(i-1) + x(i-1)
        end do
    end if
    t = 2 * t / seg%intervall
    v = matmul(seg%cache, t)
end function getvelocity_lowlevel

function getvelocity_tdb(s, segnum, tdb, tdb2, err) result(v)
    type(spk), intent(inout) :: s
    integer, intent(in) :: segnum
    real(dp), intent(in) :: tdb
    real(dp), intent(in) :: tdb2
    type(exception), intent(inout), optional :: err
    real(dp), dimension(3) :: v

    real(dp), dimension(:), allocatable :: x
    real(dp) :: twotc
    type(exception) :: err_

    call coefficients(s, segnum, tdb, tdb2, x, twotc, err_)
    if (iserror(err_)) then
        call catch(err_, "getvelocity_tdb", __FILE__, __LINE__)
        if (present(err)) then
            err = err_
            return
        else
            call raise(err_)
        end if
    end if
    v = getvelocity_lowlevel(s%segments(segnum), x, twotc)
end function getvelocity_tdb

function getstate_tdb(s, segnum, tdb, tdb2, err) result(st)
    type(spk), intent(inout) :: s
    integer, intent(in) :: segnum
    real(dp), intent(in) :: tdb
    real(dp), intent(in) :: tdb2
    type(exception), intent(inout), optional :: err
    real(dp), dimension(6) :: st

    real(dp), dimension(:), allocatable :: x
    real(dp) :: twotc
    type(exception) :: err_

    call coefficients(s, segnum, tdb, tdb2, x, twotc, err_)
    if (iserror(err_)) then
        call catch(err_, "getvelocity_tdb", __FILE__, __LINE__)
        if (present(err)) then
            err = err_
            return
        else
            call raise(err_)
        end if
    end if
    st(1:3) = getposition_lowlevel(s%segments(segnum)%cache, x)
    st(4:6) = getvelocity_lowlevel(s%segments(segnum), x, twotc)
end function getstate_tdb

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

end module ephemeris
