program testtrajectories

use assertions
use epochs, only: epoch
use states, only: state
use trajectories
use types, only: dp

implicit none

type(trajectory) :: tra
type(state) :: s0

s0 = state(epoch(), [0._dp])

call assert_false(isdirty(tra), __LINE__)
call assert_equal(length(tra), 0, __LINE__)
call add_node(tra, 1._dp, [1._dp])
call assert(isdirty(tra), __LINE__)
call assert_equal(length(tra), 1, __LINE__)
call add_node(tra, 2._dp, [2._dp])
call assert_equal(length(tra), 2, __LINE__)
call add_node(tra, 3._dp, [3._dp])
call assert_equal(length(tra), 3, __LINE__)
call save_trajectory(tra)
call assert_almost_equal(tra%t, [1._dp, 2._dp, 3._dp], __LINE__)
call assert_almost_equal(tra%y(1,:), [1._dp, 2._dp, 3._dp], __LINE__)
call assert_false(isdirty(tra), __LINE__)
end program testtrajectories
