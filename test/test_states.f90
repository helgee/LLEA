program teststates

use assertions
use bodies
use constants
use ephemerides, only: init_ephemeris
use epochs
use math, only: mod2pi
use types, only: dp
use states

implicit none

type(epoch) :: ep
real(dp), dimension(6,6) :: ref
real(dp), dimension(3,3) :: m
real(dp), dimension(3,3) :: dm
real(dp), dimension(6) :: refvec = [1._dp, 2._dp, 3._dp, 4._dp, 5._dp, 6._dp]
type(state) :: s
type(state) :: s1
real(dp) :: mu
real(dp), dimension(6) :: rvexp
real(dp), dimension(6) :: elexp
real(dp), dimension(6) :: rv
real(dp), dimension(6) :: el

ep = epoch_(2000, 1, 1, 12) + epochdelta_(seconds=1000._dp)
call init_constants
call init_ephemeris

! Reference data from WebGeocalc (http://wgc.jpl.nasa.gov/)
ref = transpose(reshape([ &
    0.93161546_dp, &
    -0.27117375_dp, &
    -0.24198643_dp, &
    0.00000000E+00_dp, &
    0.00000000E+00_dp, &
    0.00000000E+00_dp, &
    0.35177121_dp, &
    0.84016564_dp, &
    0.41276956_dp, &
    0.00000000E+00_dp, &
    0.00000000E+00_dp, &
    0.00000000E+00_dp, &
    0.09137641_dp, &
    -0.46966636_dp, &
    0.87810242_dp, &
    0.00000000E+00_dp, &
    0.00000000E+00_dp, &
    0.00000000E+00_dp, &
    4.36141153E-07_dp, &
    1.04167363E-06_dp, &
    5.11769580E-07_dp, &
    0.93161546_dp, &
    -0.27117375_dp, &
    -0.24198643_dp, &
    -1.15505708E-06_dp, &
    3.36212903E-07_dp, &
    3.00025285E-07_dp, &
    0.35177121_dp, &
    0.84016564_dp, &
    0.41276956_dp, &
    -8.06548746E-14_dp, &
    -3.99346680E-14_dp, &
    -1.29666162E-14_dp, &
    0.09137641_dp, &
    -0.46966636_dp, &
    0.87810242_dp &
], [6, 6]))
call iaumatrix(mercury, ep, m, dm)
call assert_almost_equal(m, ref(1:3, 1:3), __LINE__)
call assert_almost_equal(dm, ref(4:6, 1:3), __LINE__)

ref = transpose(reshape([ &
    -0.95473270_dp, &
    0.26677448_dp, &
    0.13159348_dp, &
    0.00000000E+00_dp, &
    0.00000000E+00_dp, &
    0.00000000E+00_dp, &
    -0.29687729_dp, &
    -0.88233398_dp, &
    -0.36517205_dp, &
    0.00000000E+00_dp, &
    0.00000000E+00_dp, &
    0.00000000E+00_dp, &
    0.01869081_dp, &
    -0.38770881_dp, &
    0.92159239_dp, &
    0.00000000E+00_dp, &
    0.00000000E+00_dp, &
    0.00000000E+00_dp, &
    8.88390262E-08_dp, &
    2.64033981E-07_dp, &
    1.09275888E-07_dp, &
    -0.95473270_dp, &
    0.26677448_dp, &
    0.13159348_dp, &
    -2.85698932E-07_dp, &
    7.98309152E-08_dp, &
    3.93786831E-08_dp, &
    -0.29687729_dp, &
    -0.88233398_dp, &
    -0.36517205_dp, &
    0.00000000E+00_dp, &
    0.00000000E+00_dp, &
    0.00000000E+00_dp, &
    0.01869081_dp, &
    -0.38770881_dp, &
    0.92159239_dp &
], [6, 6]))
call iaumatrix(venus, ep, m, dm)
call assert_almost_equal(m, ref(1:3, 1:3), __LINE__)
call assert_almost_equal(dm, ref(4:6, 1:3), __LINE__)

ref = transpose(reshape([ &
    0.24742306_dp, &
    -0.96890755_dp, &
    -7.62199719E-10_dp, &
    0.00000000E+00_dp, &
    0.00000000E+00_dp, &
    0.00000000E+00_dp, &
    0.96890755_dp, &
    0.24742306_dp, &
    -2.98477054E-09_dp, &
    0.00000000E+00_dp, &
    0.00000000E+00_dp, &
    0.00000000E+00_dp, &
    3.08055248E-09_dp, &
    -1.09209406E-17_dp, &
    1.00000000_dp, &
    0.00000000E+00_dp, &
    0.00000000E+00_dp, &
    0.00000000E+00_dp, &
    7.06538526E-05_dp, &
    1.80423738E-05_dp, &
    -9.79852602E-13_dp, &
    0.24742306_dp, &
    -0.96890755_dp, &
    -7.62199719E-10_dp, &
    -1.80423738E-05_dp, &
    7.06538526E-05_dp, &
    -2.92918995E-12_dp, &
    0.96890755_dp, &
    0.24742306_dp, &
    -2.98477054E-09_dp, &
    3.08055237E-12_dp, &
    -2.18418808E-20_dp, &
    -9.48980323E-21_dp, &
    3.08055248E-09_dp, &
    -1.09209406E-17_dp, &
    1.00000000_dp &
], [6, 6]))
call iaumatrix(earth, ep, m, dm)
call assert_almost_equal(m, ref(1:3, 1:3), __LINE__)
call assert_almost_equal(dm, ref(4:6, 1:3), __LINE__)

ref = transpose(reshape([ &
    -0.66608963_dp, &
    -0.74583619_dp, &
    -0.00727954_dp, &
    0.00000000E+00_dp, &
    0.00000000E+00_dp, &
    0.00000000E+00_dp, &
    0.59771816_dp, &
    -0.52791987_dp, &
    -0.60335198_dp, &
    0.00000000E+00_dp, &
    0.00000000E+00_dp, &
    0.00000000E+00_dp, &
    0.44615873_dp, &
    -0.40623761_dp, &
    0.79744178_dp, &
    0.00000000E+00_dp, &
    0.00000000E+00_dp, &
    0.00000000E+00_dp, &
    4.23675661E-05_dp, &
    -3.74201109E-05_dp, &
    -4.27669044E-05_dp, &
    -0.66608963_dp, &
    -0.74583619_dp, &
    -0.00727954_dp, &
    4.72138849E-05_dp, &
    5.28664951E-05_dp, &
    5.15989714E-07_dp, &
    0.59771816_dp, &
    -0.52791987_dp, &
    -0.60335198_dp, &
    -3.97806628E-14_dp, &
    -4.42633886E-13_dp, &
    -2.03232447E-13_dp, &
    0.44615873_dp, &
    -0.40623761_dp, &
    0.79744178_dp &
], [6, 6]))
call iaumatrix(mars, ep, m, dm)
call assert_almost_equal(m, ref(1:3, 1:3), __LINE__)
call assert_almost_equal(dm, ref(4:6, 1:3), __LINE__)

ref = transpose(reshape([ &
    0.39505525_dp, &
    -0.83169183_dp, &
    -0.39015388_dp, &
    0.00000000E+00_dp, &
    0.00000000E+00_dp, &
    0.00000000E+00_dp, &
    0.91854138_dp, &
    0.35086343_dp, &
    0.18214441_dp, &
    0.00000000E+00_dp, &
    0.00000000E+00_dp, &
    0.00000000E+00_dp, &
    -0.01459729_dp, &
    -0.43032959_dp, &
    0.90255380_dp, &
    0.00000000E+00_dp, &
    0.00000000E+00_dp, &
    0.00000000E+00_dp, &
    0.00016153_dp, &
    6.17004691E-05_dp, &
    3.20306845E-05_dp, &
    0.39505525_dp, &
    -0.83169183_dp, &
    -0.39015388_dp, &
    -6.94717442E-05_dp, &
    0.00014626_dp, &
    6.86098223E-05_dp, &
    0.91854138_dp, &
    0.35086343_dp, &
    0.18214441_dp, &
    -5.69989759E-14_dp, &
    4.64700233E-15_dp, &
    1.29378658E-15_dp, &
    -0.01459729_dp, &
    -0.43032959_dp, &
    0.9025538_dp &
], [6, 6]))
call iaumatrix(jupiter, ep, m, dm)
call assert_almost_equal(m, ref(1:3, 1:3), __LINE__)
call assert_almost_equal(dm, ref(4:6, 1:3), __LINE__, atol=1e-8_dp)

ref = transpose(reshape([ &
    -0.99620412_dp, &
    0.02275182_dp, &
    0.08402210_dp, &
    0.00000000E+00_dp, &
    0.00000000E+00_dp, &
    0.00000000E+00_dp, &
    -0.01645380_dp, &
    -0.99705510_dp, &
    0.07490256_dp, &
    0.00000000E+00_dp, &
    0.00000000E+00_dp, &
    0.00000000E+00_dp, &
    0.08547883_dp, &
    0.07323576_dp, &
    0.99364475_dp, &
    0.00000000E+00_dp, &
    0.00000000E+00_dp, &
    0.00000000E+00_dp, &
    -2.69488600E-06_dp, &
    -0.00016330_dp, &
    1.22679155E-05_dp, &
    -0.99620412_dp, &
    0.02275182_dp, &
    0.08402210_dp, &
    0.00016316_dp, &
    -3.72640641E-06_dp, &
    -1.37615586E-05_dp, &
    -0.01645380_dp, &
    -0.99705510_dp, &
    0.07490256_dp, &
    3.12743406E-14_dp, &
    -2.71699092E-15_dp, &
    -2.49013867E-15_dp, &
    0.08547883_dp, &
    0.07323576_dp, &
    0.99364475_dp &
], [6, 6]))
call iaumatrix(saturn, ep, m, dm)
call assert_almost_equal(m, ref(1:3, 1:3), __LINE__)
call assert_almost_equal(dm, ref(4:6, 1:3), __LINE__, atol=1e-8_dp)

ref = transpose(reshape([ &
    -0.91000094_dp, &
    0.28785233_dp, &
    -0.29839458_dp, &
    0.00000000E+00_dp, &
    0.00000000E+00_dp, &
    0.00000000E+00_dp, &
    0.35630671_dp, &
    0.17494968_dp, &
    -0.91784429_dp, &
    0.00000000E+00_dp, &
    0.00000000E+00_dp, &
    0.00000000E+00_dp, &
    -0.21199958_dp, &
    -0.94155916_dp, &
    -0.26176809_dp, &
    0.00000000E+00_dp, &
    0.00000000E+00_dp, &
    0.00000000E+00_dp, &
    -3.60714918E-05_dp, &
    -1.77114146E-05_dp, &
    9.29199819E-05_dp, &
    -0.91000094_dp, &
    0.28785233_dp, &
    -0.29839458_dp, &
    -9.21259430E-05_dp, &
    2.91413627E-05_dp, &
    -3.02086310E-05_dp, &
    0.35630671_dp, &
    0.17494968_dp, &
    -0.91784429_dp, &
    0.00000000E+00_dp, &
    0.00000000E+00_dp, &
    0.00000000E+00_dp, &
    -0.21199958_dp, &
    -0.94155916_dp, &
    -0.26176809_dp &
], [6, 6]))
call iaumatrix(uranus, ep, m, dm)
call assert_almost_equal(m, ref(1:3, 1:3), __LINE__)
call assert_almost_equal(dm, ref(4:6, 1:3), __LINE__, atol=1e-8_dp)

ref = transpose(reshape([ &
    0.16782169_dp, &
    -0.67394650_dp, &
    -0.71946646_dp, &
    0.00000000E+00_dp, &
    0.00000000E+00_dp, &
    0.00000000E+00_dp, &
    0.91829122_dp, &
    0.37233151_dp, &
    -0.13457519_dp, &
    0.00000000E+00_dp, &
    0.00000000E+00_dp, &
    0.00000000E+00_dp, &
    0.35857651_dp, &
    -0.63809510_dp, &
    0.68136446_dp, &
    0.00000000E+00_dp, &
    0.00000000E+00_dp, &
    0.00000000E+00_dp, &
    9.94860666E-05_dp, &
    4.03377447E-05_dp, &
    -1.45796407E-05_dp, &
    0.16782169_dp, &
    -0.67394650_dp, &
    -0.71946646_dp, &
    -1.81815093E-05_dp, &
    7.30141884E-05_dp, &
    7.79457377E-05_dp, &
    0.91829122_dp, &
    0.37233151_dp, &
    -0.13457519_dp, &
    2.28629820E-12_dp, &
    1.20926573E-12_dp, &
    -7.07202958E-14_dp, &
    0.35857651_dp, &
    -0.63809510_dp, &
    0.68136446_dp &
], [6, 6]))
call iaumatrix(neptune, ep, m, dm)
call assert_almost_equal(m, ref(1:3, 1:3), __LINE__)
call assert_almost_equal(dm, ref(4:6, 1:3), __LINE__, atol=1e-8_dp)

ref = transpose(reshape([ &
    0.78257370_dp, &
    0.55976292_dp, &
    0.27247730_dp, &
    0.00000000E+00_dp, &
    0.00000000E+00_dp, &
    0.00000000E+00_dp, &
    -0.62214730_dp, &
    0.71906872_dp, &
    0.30963351_dp, &
    0.00000000E+00_dp, &
    0.00000000E+00_dp, &
    0.00000000E+00_dp, &
    -0.02260855_dp, &
    -0.41183206_dp, &
    0.91097926_dp, &
    0.00000000E+00_dp, &
    0.00000000E+00_dp, &
    0.00000000E+00_dp, &
    -1.65612960E-06_dp, &
    1.91382763E-06_dp, &
    8.24852982E-07_dp, &
    0.78257370_dp, &
    0.55976292_dp, &
    0.27247730_dp, &
    -2.08318236E-06_dp, &
    -1.49049154E-06_dp, &
    -7.24341598E-07_dp, &
    -0.62214730_dp, &
    0.71906872_dp, &
    0.30963351_dp, &
    1.24731231E-10_dp, &
    -1.15606661E-09_dp, &
    -5.19534660E-10_dp, &
    -0.02260855_dp, &
    -0.41183206_dp, &
    0.91097926_dp &
], [6, 6]))
call iaumatrix(moon, ep, m, dm)
call assert_almost_equal(m, ref(1:3, 1:3), __LINE__)
call assert_almost_equal(dm, ref(4:6, 1:3), __LINE__, atol=1e-8_dp)

s = state(ep, [1000._dp, 1000._dp, 1000._dp, 1._dp, 1._dp, 1._dp], "GCRF", moon)
refvec = [1614.8139211605944_dp,406.5549322139076_dp,476.5386522817342_dp,1.6158956790362593_dp,0.4022560154092351_dp,0.4765405580099757_dp]
s1 = rotate(s, "IAU")
call assert_almost_equal(s1%rv, refvec, __LINE__)
call rotate_inplace(s1, "GCRF")
call assert_almost_equal(s1%rv, s%rv, __LINE__)

! Example from Vallado, reference values from Orekit
mu = 3.986004415e5_dp
rvexp = [6524.834_dp, 6862.875_dp, 6448.296_dp, 4.901327_dp, 5.533756_dp, 1.976341_dp]
elexp = [36127.33776397479_dp, 0.970200843222015_dp, 1.4857581497380534_dp, &
    -2.2728703614097303_dp, -0.026581796620510634_dp, 2.5678767183357816_dp]
elexp(3:) = mod2pi(elexp(3:))
el = keplerian(rvexp, mu)
call assert_almost_equal(el, elexp, __LINE__)
rv = cartesian(el, mu)
call assert_almost_equal(rv, rvexp, __LINE__)

ep = epoch_(2000, 1, 1)
s = state_(ep, [1000._dp, 1000._dp, 1000._dp, 1._dp, 0._dp, 0._dp])
rvexp = [-2.3320498362289613e8_dp, 1.3379366576629393e8_dp, 6.365110638632974e7_dp, &
    -30.134836845536128_dp, -28.69614433942715_dp, -13.006923220167772_dp]
s1 = wrt(s, mars)
call assert_almost_equal(s1%rv, rvexp, __LINE__)

end program teststates
