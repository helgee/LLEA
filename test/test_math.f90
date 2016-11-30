program testmath

use assertions
use math
use types, only: dp

implicit none

call assert_almost_equal(deg2rad(180._dp), pi, __LINE__)
call assert_almost_equal(deg2rad([180._dp, 180._dp]), [pi, pi], __LINE__)
call assert_almost_equal(rad2deg(pi), 180._dp, __LINE__)
call assert_almost_equal(rad2deg([pi, pi]), [180._dp, 180._dp], __LINE__)

call assert(isin(3._dp, [1._dp, 2._dp, 3._dp]), __LINE__)
call assert(isin(3, [1, 2, 3]), __LINE__)
call assert_false(isin(4, [1, 2, 3]), __LINE__)
call assert(isin(3._dp, reshape([1._dp, 2._dp, 3._dp, 4._dp], [2,2])), __LINE__)
call assert(isin(3, reshape([1, 2, 3, 4], [2,2])), __LINE__)

end program testmath
