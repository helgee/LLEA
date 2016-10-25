module constants

use types, only: dp
use bodies, only: body
use epochs, only: days_per_century
use math, only: deg2rad

implicit none

! planets
type(body) :: mercury
type(body) :: venus
type(body) :: earth
type(body) :: mars
type(body) :: jupiter
type(body) :: saturn
type(body) :: uranus
type(body) :: neptune

! satellites
type(body) :: moon
!DEC$ ATTRIBUTES DLLEXPORT :: mercury, venus, earth, mars, jupiter, saturn, uranus, neptune, moon, au

real(dp), parameter :: au = 0.149597870700000000e+09_dp

private

public :: mercury, venus, earth, mars, jupiter, saturn, uranus, neptune, moon, au, init_constants

contains

subroutine init_constants
    mercury%name = "Mercury"
    mercury%parent = "Sun"
    mercury%bodytype = "Planet"
    mercury%id = 199
    mercury%mu = 2.2032e4_dp
    mercury%j2 = 2.027e-4_dp
    mercury%radii = [2439.7_dp, 2439.7_dp, 2439.7_dp]
    mercury%ra = deg2rad([281.0097_dp, -0.0328_dp, 0._dp])
    mercury%dec = deg2rad([61.4143_dp, -0.0049_dp, 0._dp])
    mercury%ww = deg2rad([329.5469_dp, 6.1385025_dp, 0._dp])
    allocate(mercury%a(5))
    mercury%a = 0._dp
    allocate(mercury%d(5))
    mercury%d = 0._dp
    mercury%w = deg2rad([0.00993822_dp, -0.00104581_dp, -0.00010280_dp, -0.00002364_dp, -0.00000532_dp])
    mercury%theta0 = deg2rad([174.791086_dp, 349.582171_dp, 164.373257_dp, 339.164343_dp, 153.955429_dp])
    mercury%theta1 = deg2rad([4.092335_dp, 8.184670_dp, 12.277005_dp, 16.369340_dp, 20.461675_dp]) &
        * days_per_century

    venus%name = "Venus"
    venus%parent = "Sun"
    venus%bodytype = "Planet"
    venus%id = 299
    venus%mu = 3.24859e5_dp
    venus%j2 = 6e-5_dp
    venus%radii = [6051.8_dp, 6051.8_dp, 6051.8_dp]
    venus%ra = deg2rad([272.76_dp, 0._dp, 0._dp])
    venus%dec = deg2rad([67.16_dp, 0._dp, 0._dp])
    venus%ww = deg2rad([160.20_dp, -1.4813688_dp, 0._dp])
    allocate(venus%a(1))
    allocate(venus%d(1))
    allocate(venus%w(1))
    allocate(venus%theta0(1))
    allocate(venus%theta1(1))
    venus%a = 0._dp
    venus%d = 0._dp
    venus%w = 0._dp
    venus%theta0 = 0._dp
    venus%theta1 = 0._dp

    earth%name = "earth"
    earth%parent = "sun"
    earth%bodytype = "planet"
    earth%id = 399
    earth%mu = 3.986004418e5_dp
    earth%j2 = 1.08262668e-3_dp
    earth%radii = [6371.0084_dp, 6378.1366_dp, 6356.7519_dp]
    earth%ra = deg2rad([0._dp, -0.641_dp, 0._dp])
    earth%dec = deg2rad([90._dp, -0.557_dp, 0._dp])
    earth%ww = deg2rad([190.147_dp, 360.9856235_dp, 0._dp])
    allocate(earth%a(1))
    allocate(earth%d(1))
    allocate(earth%w(1))
    allocate(earth%theta0(1))
    allocate(earth%theta1(1))
    earth%a = 0._dp
    earth%d = 0._dp
    earth%w = 0._dp
    earth%theta0 = 0._dp
    earth%theta1 = 0._dp

    mars%name = "Mars"
    mars%parent = "Sun"
    mars%bodytype = "Planet"
    mars%id = 4
    mars%mu = 4.282837e4_dp
    mars%j2 = 1.964e-3_dp
    mars%radii = [3389.5_dp, 3396.19_dp, 3376.20_dp]
    mars%ra = deg2rad([317.68143_dp, -0.1061_dp, 0._dp])
    mars%dec = deg2rad([52.88650_dp, -0.0609_dp, 0._dp])
    mars%ww = deg2rad([176.630_dp, 350.89198226_dp, 0._dp])
    allocate(mars%a(1))
    allocate(mars%d(1))
    allocate(mars%w(1))
    allocate(mars%theta0(1))
    allocate(mars%theta1(1))
    mars%a = 0._dp
    mars%d = 0._dp
    mars%w = 0._dp
    mars%theta0 = 0._dp
    mars%theta1 = 0._dp

    jupiter%name = "Jupiter"
    jupiter%parent = "Sun"
    jupiter%bodytype = "Planet"
    jupiter%id = 5
    jupiter%mu = 1.26686534e8_dp
    jupiter%j2 = 1.475e-2_dp
    jupiter%radii = [69911._dp, 71492._dp, 66854._dp]
    jupiter%ra = deg2rad([268.056595_dp, -0.006499_dp, 0._dp])
    jupiter%dec = deg2rad([64.495303_dp, 0.002413_dp, 0._dp])
    jupiter%ww = deg2rad([284.95_dp, 870.536_dp, 0._dp])
    allocate(jupiter%w(5))
    jupiter%w = 0._dp
    jupiter%a = deg2rad([0.000117_dp, 0.000938_dp, 0.001432_dp, 0.00003_dp, 0.002150_dp])
    jupiter%d = deg2rad([0.00005_dp, 0.000404_dp, 0.000617_dp, -0.000013_dp, 0.000926_dp])
    jupiter%theta0 = deg2rad([99.360714_dp, 175.895369_dp, 300.323162_dp, 114.012305_dp, 49.511251_dp])
    jupiter%theta1 = deg2rad([4850.4046_dp, 1191.9605_dp, 262.5475_dp, 6070.2476_dp, 64.3_dp])

    saturn%name = "Saturn"
    saturn%parent = "Sun"
    saturn%bodytype = "Planet"
    saturn%id = 6
    saturn%mu = 3.7931187e7_dp
    saturn%j2 = 1.645e-2_dp
    saturn%radii = [58232._dp, 60268._dp, 54364._dp]
    saturn%ra = deg2rad([40.589_dp, -0.036_dp, 0._dp])
    saturn%dec = deg2rad([83.537_dp, -0.004_dp, 0._dp])
    saturn%ww = deg2rad([38.90_dp, 810.7939024_dp, 0._dp])
    allocate(saturn%a(1))
    allocate(saturn%d(1))
    allocate(saturn%w(1))
    allocate(saturn%theta0(1))
    allocate(saturn%theta1(1))
    saturn%a = 0._dp
    saturn%d = 0._dp
    saturn%w = 0._dp
    saturn%theta0 = 0._dp
    saturn%theta1 = 0._dp

    uranus%name = "Uranus"
    uranus%parent = "Sun"
    uranus%bodytype = "Planet"
    uranus%id = 7
    uranus%mu = 5.793939e6_dp
    uranus%j2 = 1.2e-2_dp
    uranus%radii = [25362._dp, 25559._dp, 24973._dp]
    uranus%ra = deg2rad([257.311_dp, 0._dp, 0._dp])
    uranus%dec = deg2rad([-15.175_dp, 0._dp, 0._dp])
    uranus%ww = deg2rad([203.81_dp, -501.1600928_dp, 0._dp])
    allocate(uranus%a(1))
    allocate(uranus%d(1))
    allocate(uranus%w(1))
    allocate(uranus%theta0(1))
    allocate(uranus%theta1(1))
    uranus%a = 0._dp
    uranus%d = 0._dp
    uranus%w = 0._dp
    uranus%theta0 = 0._dp
    uranus%theta1 = 0._dp

    neptune%name = "Neptune"
    neptune%parent = "Sun"
    neptune%bodytype = "Planet"
    neptune%id = 8
    neptune%mu = 6.836529e6_dp
    neptune%j2 = 4e-3_dp
    neptune%radii = [24622._dp, 24764._dp, 24341._dp]
    neptune%ra = deg2rad([299.36_dp, 0._dp, 0._dp])
    neptune%dec = deg2rad([43.46_dp, 0._dp, 0._dp])
    neptune%ww = deg2rad([253.18_dp, 536.3128492_dp, 0._dp])
    neptune%a = deg2rad([0.7_dp])
    neptune%d = deg2rad([-0.51_dp])
    neptune%w = deg2rad([-0.48_dp])
    neptune%theta0 = deg2rad([357.85_dp])
    neptune%theta1 = deg2rad([52.316_dp])

    moon%name = "Moon"
    moon%parent = "Earth"
    moon%bodytype = "Satellite"
    moon%id = 301
    moon%mu = 4902.80007_dp
    moon%j2 = 203.21568e-6_dp
    moon%radii = [1737.4_dp, 1737.4_dp, 1737.4_dp]
    moon%ra = deg2rad([269.9949_dp, 0.0031_dp, 0._dp])
    moon%dec = deg2rad([66.5392_dp, 0.01300_dp, 0._dp])
    moon%ww = deg2rad([38.3213_dp, 13.17635815_dp, -1.4e-12_dp])
    moon%a = deg2rad([-3.8787_dp, -0.1204_dp, 0.07_dp, -0.0172_dp, 0.0_dp, &
        0.0072_dp, 0.0_dp, 0.0_dp, 0.0_dp, -0.0052_dp, 0.0_dp, 0.0_dp, 0.0043_dp])
    moon%d = deg2rad([1.5419_dp, 0.0239_dp, -0.0278_dp, 0.0068_dp, 0.0_dp, &
        -0.0029_dp, 0.0009_dp, 0.0_dp, 0.0_dp, 0.0008_dp, 0.0_dp, 0.0_dp, -0.0009_dp])
    moon%w = deg2rad([3.5610_dp, 0.1208_dp, -0.0642_dp, 0.0158_dp, 0.0252_dp, &
        -0.0066_dp, -0.0047_dp, -0.0046_dp, 0.0028_dp, 0.0052_dp, 0.004_dp, 0.0019_dp, -0.0044_dp])
    moon%theta0 = deg2rad([125.045_dp, 250.089_dp, 260.008_dp, 176.625_dp, &
        357.529_dp, 311.589_dp, 134.963_dp, 276.617_dp, 34.226_dp, 15.134_dp, &
        119.743_dp, 239.961_dp, 25.053_dp])
    moon%theta1 = deg2rad([-0.0529921_dp, -0.1059842_dp, 13.0120009_dp, &
        13.3407154_dp, 0.9856003_dp, 26.4057084_dp, 13.0649930_dp, 0.3287146_dp, &
        1.7484877_dp, -0.1589763_dp, 0.0036096_dp, 0.1643573_dp, 12.9590088_dp]) * days_per_century
end subroutine init_constants

end module constants
