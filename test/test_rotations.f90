program testrotations

use assertions
use exceptions
use rotations
use types, only: dp

implicit none

type(exception) :: err
real(dp), dimension(3,3) :: m

m = rotationmatrix("J", 0._dp, err)
call assert_raises("Invalid rotation axis: 'J'.", err, __LINE__)
call reset(err)
m = ratematrix("J", 0._dp, 0._dp, err)
call assert_raises("Invalid rotation axis: 'J'.", err, __LINE__)
call reset(err)
err = checkaxes("JXZ")
call assert_raises("Invalid rotation axis: 'J'.", err, __LINE__)
err = checkaxes("XJZ")
call assert_raises("Invalid rotation axis: 'J'.", err, __LINE__)
err = checkaxes("XZJ")
call assert_raises("Invalid rotation axis: 'J'.", err, __LINE__)
err = checkaxes("XZZ")
call assert_raises("Subsequent rotations around the same axis are meaningless.", err, __LINE__)
err = checkaxes("XXZ")
call assert_raises("Subsequent rotations around the same axis are meaningless.", err, __LINE__)

end program testrotations
