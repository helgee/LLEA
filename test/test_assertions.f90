program test_assertions

use types, only: dp
use math, only: isclose
use assertions

implicit none

real(dp) :: s1 = 0._dp
real(dp) :: s2 = 1._dp
real(dp), dimension(3) :: v1 = [0._dp, 0._dp, 0._dp]
real(dp), dimension(3) :: v2 = [1._dp, 1._dp, 1._dp]

call assert(isclose(s1, s1), __LINE__)
call assert_false(s1 == s2, __LINE__)
call assert(v1 == v1, __LINE__)
call assert_false(v1 == v2, __LINE__)
call assert_equal(s1, s1, __LINE__)
call assert_not_equal(s1, s2, __LINE__)
call assert_equal(v1, v1, __LINE__)
call assert_not_equal(v1, v2, __LINE__)
call assert_almost_equal(s1, s1+epsilon(0._dp), __LINE__)
call assert_almost_equal(v1, v1+epsilon(0._dp), __LINE__)
end program test_assertions
