program testmath

use assertions
use math
use types, only: dp

implicit none

call assert_almost_equal(deg2rad(180._dp), pi, __LINE__)
call assert_almost_equal(deg2rad([180._dp, 180._dp]), [pi, pi], __LINE__)
call assert_almost_equal(rad2deg(pi), 180._dp, __LINE__)
call assert_almost_equal(rad2deg([pi, pi]), [180._dp, 180._dp], __LINE__)

end program testmath
