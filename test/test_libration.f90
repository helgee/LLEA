program testlibration

use assertions
use constants
use epochs, only: epoch
use ephemerides, only: init_ephemeris
use libration
use types, only: dp

implicit none

real(dp) :: l1ref
real(dp) :: l2ref
real(dp) :: l3ref
real(dp) :: l1dist
real(dp) :: l2dist
real(dp) :: l3dist
real(dp), dimension(6) :: ref
real(dp), dimension(6) :: rv
real(dp), dimension(6) :: lib
real(dp), dimension(6) :: eci
type(epoch) :: ep

call init_constants
call init_ephemeris

l1ref = 0.15093428265127021_dp
l1dist = librationdist(earth%mu, moon%mu, "L1")
call assert_almost_equal(l1ref, l1dist, __LINE__)
l2ref = 0.16783274367911588_dp
l2dist = librationdist(earth%mu, moon%mu, "L2")
call assert_almost_equal(l2ref, l2dist, __LINE__)
l3ref = 0.99291206108848129_dp
l3dist = librationdist(earth%mu, moon%mu, "L3")
call assert_almost_equal(l3ref, l3dist, __LINE__)

ref = [-1.1717170062386328_dp, -1.1212220608364672E-003_dp, 1.6674369835388907E-003_dp, &
    -1.6420893636515004_dp, -0.45661886467983320_dp, 0.68484198823320419_dp]
rv = 1._dp
rv(1:3) = rv(1:3) * 1000._dp
ep = epoch(2000, 1, 1)
lib = gcrftolib(rv, ep, earth, moon, "L2")
call assert_almost_equal(lib, ref, __LINE__)
eci = libtogcrf(lib, ep, earth, moon, "L2")
call assert_almost_equal(eci, rv, __LINE__)

lib = gcrftolib(rv, ep, earth, moon, "L2", normalise=.false.)
eci = libtogcrf(lib, ep, earth, moon, "L2", normalise=.false.)
call assert_almost_equal(eci, rv, __LINE__)

end program testlibration
