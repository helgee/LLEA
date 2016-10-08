program testephemeris

use assertions
use ephemeris
use epochs
use types, only: dp

implicit none

character(len=14) :: path
integer :: u
integer :: stat
type(epoch) :: ep
type(daf) :: d
type(spk) :: s
real(dp), dimension(6) :: ref405
real(dp), dimension(6) :: ref430
real(dp), dimension(6) :: res405
real(dp), dimension(6) :: res430

ep = epoch(2000, 1, 1, 12)
ref405 = [-2.0529325063213453e7_dp,-6.032395669044396e7_dp,-3.013084448833336e7_dp,37.00430439938116_dp,-8.541375986872962_dp,-8.398373128367817_dp]
ref430 = [-2.052932489501512e7_dp,-6.0323956764362626e7_dp,-3.0130843855883405e7_dp,37.00430445042317_dp,-8.541376874554961_dp,-8.398372276759316_dp]
open(newunit=u, file="data/de430.bsp", status="old", form="unformatted", access="direct", recl=1024, iostat=stat)
call assert_equal(stat, 0, __LINE__)
close(u)

path = "data/de405.bsp"
d = daf(path)
call assert_equal(d%id, "NAIF/DAF", __LINE__)
call assert_equal(d%endianness, "big_endian", __LINE__)
s = spk(path)
res405 = getstate(s, 1, ep%jd, ep%jd1)
call assert_almost_equal(res405, ref405, __LINE__)

path = "data/de430.bsp"
d = daf(path)
call assert_equal(d%id, "DAF/SPK", __LINE__)
call assert_equal(d%endianness, "little_endian", __LINE__)
s = spk(path)
res430 = getstate(s, 1, ep%jd, ep%jd1)
call assert_almost_equal(res430, ref430, __LINE__)

end program testephemeris
