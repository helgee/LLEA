program testephemeris

use assertions
use constants, only: au
use ephemerides
use epochs
use exceptions
use types, only: dp
use util, only: joinpath, projectdir

implicit none

character(len=:), allocatable :: path
integer, dimension(:), allocatable :: ipath
type(epoch) :: ep
type(daf) :: d
type(jplephem) :: eph
real(dp), dimension(6) :: ref405
real(dp), dimension(6) :: ref430
real(dp), dimension(6) :: res6
real(dp), dimension(3) :: res3

ep = epoch(2000, 1, 1, 12)
ref405 = [-2.0529325063213453e7_dp,-6.032395669044396e7_dp,-3.013084448833336e7_dp,37.00430439938116_dp,-8.541375986872962_dp,-8.398373128367817_dp]
ref430 = [-2.052932489501512e7_dp,-6.0323956764362626e7_dp,-3.0130843855883405e7_dp,37.00430445042317_dp,-8.541376874554961_dp,-8.398372276759316_dp]

call assert_equal(naifid("ssb"), 0, __LINE__)

ipath = [301]
call getpath(ipath, 399)
call assert_equal(ipath, [301, 3, 399], __LINE__)
ipath = [301]
call getpath(ipath, 0)
call assert_equal(ipath, [301, 3, 0], __LINE__)
ipath = [301]
call getpath(ipath, 1)
call assert_equal(ipath, [301, 3, 0, 1], __LINE__)
ipath = [301]
call getpath(ipath, 199)
call assert_equal(ipath, [301, 3, 0, 1, 199], __LINE__)

path = joinpath(projectdir("data"), "de430.bsp")
d = daf(path)
call assert_equal(d%id, "DAF/SPK", __LINE__)
call assert_equal(d%endianness, "little_endian", __LINE__)
eph = jplephem(path)

res3 = getposition(eph, ep, 1)
call assert_almost_equal(res3, ref430(1:3), __LINE__)
res3 = getposition(eph, ep, 1, from=0)
call assert_almost_equal(res3, ref430(1:3), __LINE__)
res3 = getposition(eph, ep, "mercury barycenter", from="ssb")
call assert_almost_equal(res3, ref430(1:3), __LINE__)
res3 = getvelocity(eph, ep, 1)
call assert_almost_equal(res3, ref430(4:6), __LINE__)
res3 = getvelocity(eph, ep, 1, from=0)
call assert_almost_equal(res3, ref430(4:6), __LINE__)
res3 = getvelocity(eph, ep, "mercury barycenter", from="ssb")
call assert_almost_equal(res3, ref430(4:6), __LINE__)
res6 = getstate(eph, ep, 1)
call assert_almost_equal(res6, ref430, __LINE__)
res6 = getstate(eph, ep, 1, from=0)
call assert_almost_equal(res6, ref430, __LINE__)
res6 = getstate(eph, ep, "mercury barycenter", from="ssb")
call assert_almost_equal(res6, ref430, __LINE__)

path = joinpath(projectdir("data"), "testpo.430")
call jpltest(eph, path)

path = joinpath(projectdir("data"), "de405.bsp")
d = daf(path)
call assert_equal(d%id, "NAIF/DAF", __LINE__)
call assert_equal(d%endianness, "big_endian", __LINE__)
eph = jplephem(path)
res6 = eph%state(ep, 1)
call assert_almost_equal(res6, ref405, __LINE__)

path = joinpath(projectdir("data"), "testpo.405")
call jpltest(eph, path)

contains

subroutine jpltest(ephem, path)
    type(jplephem), intent(inout) :: ephem
    character(len=*), intent(in) :: path

    integer :: u
    character(len=100) :: tmp
    integer :: de
    character(len=10) :: date
    real(dp) :: jed
    integer :: t
    integer :: c
    integer :: x
    real(dp) :: val
    integer :: stat
    real(dp), dimension(6) :: ts
    real(dp), dimension(6) :: cs
    real(dp), dimension(6) :: st
    type(exception) :: err
    type(epoch) :: ep

    stat = 0
    tmp = ""

    open(newunit=u, file=path, status="old")
    do while (tmp /= "EOT")
        read(u, *) tmp
    end do
    do while (stat == 0)
        read(u, *, iostat=stat) de, date, jed, t, c, x, val
        write(*, *) de, date, jed, t, c, x, val
        if ((t == 14).or.(t == 15)) cycle

        ep = epoch(jed)
        ts = teststate(ephem, t, ep, err)
        if (iserror(err).and.hasid(err, "OutOfRangeError")) then
            call reset(err)
            cycle
        end if
        cs = teststate(ephem, c, ep, err)
        if (iserror(err).and.hasid(err, "OutOfRangeError")) then
            call reset(err)
            cycle
        end if

        ts(4:6) = ts(4:6) * 86400
        cs(4:6) = cs(4:6) * 86400
        st = (ts - cs) / au

        call assert_almost_equal(st(x), val, __LINE__, atol=1e-13_dp)
    end do
    close(u)
end subroutine jpltest

function teststate(ephem, targ, tdb, err) result(st)
    type(jplephem), intent(inout) :: ephem
    integer, intent(in) :: targ
    type(epoch), intent(in) :: tdb
    type(exception) :: err

    real(dp), dimension(6) :: st
    real(dp), dimension(6) :: st1
    real(dp), dimension(6) :: st2

    if (targ == 3) then
        st1 = getstate(ephem, tdb, 3, 0, err=err)
        if (iserror(err)) return
        st2 = getstate(ephem, tdb, 399, 3, err=err)
        if (iserror(err)) return
        st = st1 + st2
    else if (targ == 10) then
        st1 = getstate(ephem, tdb, 3, 0, err=err)
        if (iserror(err)) return
        st2 = getstate(ephem, tdb, 301, 3, err=err)
        if (iserror(err)) return
        st = st1 + st2
    else if (targ == 11) then
        st = getstate(ephem, tdb, 10, 0, err=err)
        if (iserror(err)) return
    else if (targ == 12) then
        st = 0._dp
    else if (targ == 13) then
        st = getstate(ephem, tdb, 3, 0, err=err)
        if (iserror(err)) return
    else
        st = getstate(ephem, tdb, targ, 0, err=err)
        if (iserror(err)) return
    end if
end function teststate

end program testephemeris
