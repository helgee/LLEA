program testlibration

use assertions
use constants
use libration
use types, only: dp

implicit none

real(dp) :: l1ref
real(dp) :: l2ref
real(dp) :: l3ref
real(dp) :: l1dist
real(dp) :: l2dist
real(dp) :: l3dist

call init_constants

l1ref = 0.15093428265127021_dp
l1dist = librationdist(earth%mu, moon%mu, "L1")
! call assert_almost_equal(l1ref, l1dist, __LINE__)
l2ref = 0.16783274367911588_dp
l2dist = librationdist(earth%mu, moon%mu, "L2")
call assert_almost_equal(l2ref, l2dist, __LINE__)
l3ref = 0.99291206108848129_dp
l3dist = librationdist(earth%mu, moon%mu, "L3")
call assert_almost_equal(l3ref, l3dist, __LINE__)

end program testlibration
