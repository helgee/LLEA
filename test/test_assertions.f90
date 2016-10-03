program test_assertions

use types, only: dp
use math, only: isapprox
use assertions

implicit none

real(dp) :: s1 = 0._dp
real(dp) :: s2 = 1._dp
real(dp), dimension(3) :: v1 = [0._dp, 0._dp, 0._dp]
real(dp), dimension(3) :: v2 = [1._dp, 1._dp, 1._dp]

call assert(isapprox(s1, s1), __LINE__)
call assert_false(isapprox(s1, s2), __LINE__)
call assert_equal(s1, s1, __LINE__)
call assert_not_equal(s1, s2, __LINE__)
call assert_equal(v1, v1, __LINE__)
call assert_not_equal(v1, v2, __LINE__)
call assert_almost_equal(s2, s2+epsilon(0._dp), __LINE__)
call assert_almost_equal(v2, v2+epsilon(0._dp), __LINE__)
end program test_assertions
