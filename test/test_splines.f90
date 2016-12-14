program testsplines

use assertions
use exceptions
use math, only: linspace
use splines
use types, only: dp

implicit none

type(spline1d) :: spl
integer :: i
real(dp), dimension(:), allocatable :: x
real(dp), dimension(100) :: xp
real(dp), dimension(100) :: xp_zeros
real(dp), dimension(100) :: xp_clip
real(dp), dimension(:), allocatable :: y
real(dp), dimension(:), allocatable :: yi
real(dp), dimension(:), allocatable :: w
real(dp), dimension(:), allocatable :: desired
real(dp), dimension(:), allocatable :: actual
real(dp) :: ys
type(exception) :: err

x = [1._dp, 2._dp, 3._dp]
y = [0._dp, 2._dp, 4._dp]
spl = spline1d(x, y, k=1, s=real(size(x), kind=dp))
yi = knots(spl)
call assert_almost_equal(yi, [1._dp, 3._dp], __LINE__)
yi = coeffs(spl)
call assert_almost_equal(yi, [0._dp, 4._dp], __LINE__)
ys = residual(spl)
call assert_almost_equal(ys, 0._dp, __LINE__)

yi = evaluate(spl, [1._dp, 1.5_dp, 2._dp])
call assert_almost_equal(yi, [0._dp, 1._dp, 2._dp], __LINE__)
ys = evaluate(spl, 1.5_dp)
call assert_almost_equal(ys, 1._dp, __LINE__)

x = [-1._dp, -0.65016502_dp, -0.58856235_dp, -0.26903553_dp, -0.17370892_dp, &
     -0.10011001_dp, 0._dp, 0.10011001_dp, 0.17370892_dp, 0.26903553_dp, 0.58856235_dp, &
     0.65016502_dp, 1._dp]
y = [1._dp,0.62928599_dp, 0.5797223_dp, 0.39965815_dp, 0.36322694_dp, 0.3508061_dp, &
     0.35214793_dp, 0.3508061_dp, 0.36322694_dp, 0.39965815_dp, 0.5797223_dp, &
     0.62928599_dp, 1._dp]
w = [1.00000000e+12_dp, 6.88875973e+02_dp, 4.89314737e+02_dp, 4.26864807e+02_dp, &
     6.07746770e+02_dp, 4.51341444e+02_dp, 3.17480210e+02_dp, 4.51341444e+02_dp, &
     6.07746770e+02_dp, 4.26864807e+02_dp, 4.89314737e+02_dp, 6.88875973e+02_dp, &
     1.00000000e+12_dp]
spl = spline1d(x, y, w=w, s=real(size(x), kind=dp))
desired = [0.35100374_dp, 0.51715855_dp, 0.87789547_dp, 0.98719344_dp]
actual = evaluate(spl, [0.1_dp, 0.5_dp, 0.9_dp, 0.99_dp])
call assert_almost_equal(actual, desired, __LINE__)

spl = spline1d(x, y, bc="unknown", err=err)
call assert_raises("Unknown boundary condition: unknown", err, __LINE__)

x = [0._dp, 1._dp, 2._dp, 3._dp, 4._dp]
y = x**3
xp = linspace(-8._dp, 13._dp, 100)
do i = 1, 100
    if (xp(i) < 0._dp) then
        xp_zeros(i) = 0._dp
        xp_clip(i) = 0._dp
    elseif (xp(i) >= 0._dp .and. xp(i) <= 4._dp) then
        xp_zeros(i) = xp(i)
        xp_clip(i) = xp(i)
    elseif (xp(i) > 4._dp) then
        xp_zeros(i) = 0._dp
        xp_clip(i) = 4._dp
    end if
end do

spl = spline1d(x, y)
yi = evaluate(spl, xp)
call assert_almost_equal(yi, xp_clip**3, __LINE__)
spl = spline1d(x, y, bc="extrapolate")
yi = evaluate(spl, xp)
call assert_almost_equal(yi, xp**3, __LINE__)
spl = spline1d(x, y, bc="zero")
yi = evaluate(spl, xp)
call assert_almost_equal(yi, xp_zeros**3, __LINE__)
spl = spline1d(x, y, bc="error")
yi = evaluate(spl, xp, err=err)
call assert_raises("Input point out of range.", err, __LINE__)

end program testsplines
