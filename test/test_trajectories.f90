program testtrajectories

use assertions
use epochs, only: epoch
use exceptions
use states, only: state
use trajectories
use types, only: dp

implicit none

type(exception) :: err
type(trajectory) :: tra
type(state) :: s0
real(dp), dimension(:), allocatable :: x

s0 = state(epoch(), [0._dp])
tra = trajectory(s0, ["x"])

call getfield(tra, "x", x, err)
call assert_raises("Trajectory is empty.", err, __LINE__)
call assert_false(isdirty(tra), __LINE__)
call assert_equal(len_dirty(tra), 0, __LINE__)
call add_node(tra, 1._dp, [1._dp])
call add_node(tra, 1._dp, [1._dp, 1._dp], err)
call assert_raises("Wrong number of fields.", err, __LINE__)
call assert(isdirty(tra), __LINE__)
call assert_equal(len_dirty(tra), 1, __LINE__)
call add_node(tra, 2._dp, [2._dp])
call assert_equal(len_dirty(tra), 2, __LINE__)
call add_node(tra, 3._dp, [3._dp])
call assert_equal(len_dirty(tra), 3, __LINE__)
call save_trajectory(tra)
call assert_almost_equal(tra%t, [1._dp, 2._dp, 3._dp], __LINE__)
call assert_almost_equal(tra%y(1,:), [1._dp, 2._dp, 3._dp], __LINE__)
call assert_false(isdirty(tra), __LINE__)
call getfield(tra, "x", x)
call assert_almost_equal(x, [1._dp, 2._dp, 3._dp], __LINE__)
call getfield(tra, "y", x, err)
call assert_raises("Unknown field: y", err, __LINE__)
end program testtrajectories
