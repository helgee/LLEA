program testiau

use assertions
use bodies
use constants
use epochs
use math, only: pih
use rotations
use types, only: dp

implicit none

type(epoch) :: ep
real(dp) :: phi
real(dp) :: xi
real(dp) :: omega
real(dp) :: domega
real(dp) :: alpha
real(dp) :: dalpha
real(dp) :: delta
real(dp) :: ddelta
real(dp), dimension(3,3) :: iau
real(dp), dimension(3) :: res

call init_constants
ep = new_epoch(2000, 1, 1, 12) + new_epochdelta(seconds=1000._dp)
alpha = rightascension(moon, ep)
call assert_almost_equal(alpha, 4.657546533178125_dp, __LINE__)
dalpha = rightascensionrate(moon, ep)
call assert_almost_equal(dalpha, 4.556005792252131e-10_dp, __LINE__)
delta = declination(moon, ep)
call assert_almost_equal(delta, 1.1456521072552015_dp, __LINE__)
ddelta = declinationrate(moon, ep)
call assert_almost_equal(ddelta, 1.2597678243559196e-9_dp, __LINE__)
omega = rotationangle(moon, ep)
call assert_almost_equal(omega, 0.7216545147226898_dp, __LINE__)
domega = rotationrate(moon, ep)
call assert_almost_equal(domega, 2.6615172281048365e-6_dp, __LINE__)
phi = pih + alpha
call assert_almost_equal(phi, 6.228342859973021_dp, __LINE__)
xi = pih - delta
call assert_almost_equal(xi, 0.42514421953969506_dp, __LINE__)

call assert_almost_equal(rotationmatrix("3", phi), &
    reshape([0.9984965298803808_dp,0.054814959790533614_dp,0.0_dp,-0.054814959790533614_dp, &
    0.9984965298803808_dp,0.0_dp,0.0_dp,0.0_dp,1.0_dp], [3,3]), __LINE__)
call assert_almost_equal(rotationmatrix("1", xi), &
    reshape([1.0_dp,0.0_dp,0.0_dp,0.0_dp,0.9109792587664173_dp,-0.4124521670416935_dp,0.0_dp, &
    0.4124521670416935_dp,0.9109792587664173_dp], [3,3]), __LINE__)
call assert_almost_equal(rotationmatrix("3", omega), &
    reshape([0.7507137389875419_dp,-0.6606276425455907_dp,0.0_dp,0.6606276425455907_dp, &
    0.7507137389875419_dp,0.0_dp,0.0_dp,0.0_dp,1.0_dp], [3,3]), __LINE__)
call assert_almost_equal(ratematrix("3", phi, dalpha), &
    reshape([2.497372743077388e-11_dp,-4.5491559736786684e-10_dp,0.0_dp,4.5491559736786684e-10_dp, &
    2.497372743077388e-11_dp,0.0_dp,0.0_dp,0.0_dp,0.0_dp], [3,3]), __LINE__)
call assert_almost_equal(ratematrix("1", xi, -ddelta), &
    reshape([0.0_dp,0.0_dp,0.0_dp,0.0_dp,5.195939691249986e-10_dp,1.1476223588495377e-9_dp,0.0_dp, &
    -1.1476223588495377e-9_dp,5.195939691249986e-10_dp], [3,3]), __LINE__)
call assert_almost_equal(ratematrix("3", omega, domega), &
    reshape([-1.7582718519973733e-6_dp,-1.9980375496903402e-6_dp,0.0_dp,1.9980375496903402e-6_dp, &
    -1.7582718519973733e-6_dp,0.0_dp,0.0_dp,0.0_dp,0.0_dp], [3,3]), __LINE__)
call assert_almost_equal(rotationmatrix("313", phi, xi, omega), &
    reshape([0.7825736971759526_dp,-0.622147299281006_dp,-0.022608548951908884_dp, &
    0.5597629212090675_dp,0.7190687230215297_dp,-0.4118320575327742_dp, &
    0.2724773027755742_dp,0.30963350847338394_dp,0.9109792587664173_dp], [3,3]), __LINE__)
iau = ratematrix("313", phi, dalpha, xi, -ddelta, omega, domega)
call assert_almost_equal(iau, &
    reshape([-1.6560919680853432e-6_dp,-2.0831396039519395e-6_dp,2.50537797410498e-10_dp,1.9145130776034144e-6_dp, &
    -1.4897126292438954e-6_dp,1.1355964749264694e-9_dp,8.233367661469968e-7_dp,-7.26064571476686e-7_dp,5.195939691249986e-10_dp], &
    [3,3]), __LINE__)

end program testiau
