program testepochs

use assertions
use epochs
use types, only: dp

implicit none

type(epoch) :: ep0
type(epoch) :: ep1
type(epoch) :: ep2
type(epochdelta) :: epd

! Julian day conversions
call assert_almost_equal(sum(juliandate(2000, 1, 1, 12, 0, 0._dp)), 2451545._dp, __LINE__)
call assert_almost_equal(sum(juliandate(1999, 1, 1, 0, 0, 0._dp)), 2451179.5_dp, __LINE__)
call assert_almost_equal(sum(juliandate(1987, 1, 27, 0, 0, 0._dp)), 2446822.5_dp, __LINE__)
call assert_almost_equal(sum(juliandate(1987, 6, 19, 12, 0, 0._dp)), 2446966.0_dp, __LINE__)
call assert_almost_equal(sum(juliandate(1988, 1, 27, 0, 0, 0._dp)), 2447187.5_dp, __LINE__)
call assert_almost_equal(sum(juliandate(1988, 6, 19, 12, 0, 0._dp)), 2447332.0_dp, __LINE__)
call assert_almost_equal(sum(juliandate(1900, 1, 1, 0, 0, 0._dp)), 2415020.5_dp, __LINE__)
call assert_almost_equal(sum(juliandate(1600, 1, 1, 0, 0, 0._dp)), 2305447.5_dp, __LINE__)
call assert_almost_equal(sum(juliandate(1600, 12, 31, 0, 0, 0._dp)), 2305812.5_dp, __LINE__)

! epoch and epochdelta types
ep0 = new_epoch(2000, 1, 1)
ep1 = new_epoch(2000, 1, 2)
call assert_almost_equal(jd2000(ep0), 0._dp, __LINE__)
call assert_almost_equal(jd2000(ep1), 1._dp, __LINE__)
epd = new_epochdelta(days=1._dp)
call assert_almost_equal(jd2000(ep0 + epd), jd2000(ep1), __LINE__)
epd = ep1 - ep0
call assert_almost_equal(days(epd), 1._dp, __LINE__)
call assert_almost_equal(seconds(epd), seconds_per_day, __LINE__)
ep2 = new_epoch(2000, 1, 1, hour=12)
call assert_almost_equal(jd2000(ep2), 0.5_dp, __LINE__)
ep2 = new_epoch(2000, 1, 1, minute=720)
call assert_almost_equal(jd2000(ep2), 0.5_dp, __LINE__)
ep2 = new_epoch(2000, 1, 1, seconds=43200._dp)
call assert_almost_equal(jd2000(ep2), 0.5_dp, __LINE__)

! isostring and calendardate
call assert_equal(isostring(ep0), "2000-01-01T00:00:00.000", __LINE__)
call assert_equal(isostring(ep1), "2000-01-02T00:00:00.000", __LINE__)
call assert_equal(isostring(ep2), "2000-01-01T12:00:00.000", __LINE__)
epd = new_epochdelta(seconds=7.1234_dp)
call assert_equal(isostring(ep2+epd), "2000-01-01T12:00:07.123", __LINE__)
epd = new_epochdelta(days=1.5_dp)
call assert_equal(isostring(ep2+epd), "2000-01-03T00:00:00.000", __LINE__)
epd = new_epochdelta(days=3.25_dp, seconds=7.1234_dp)
call assert_equal(isostring(ep2+epd), "2000-01-04T18:00:07.123", __LINE__)
epd = new_epochdelta(days=-1.5_dp)
call assert_equal(isostring(ep2+epd), "1999-12-31T00:00:00.000", __LINE__)

end program testepochs
