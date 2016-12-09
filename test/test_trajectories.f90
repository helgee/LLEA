program testtrajectories

use assertions
use epochs, only: epoch
use states, only: state
use trajectories
use types, only: dp

implicit none

type(tranode) :: node
type(tranode), pointer :: head
type(trajectory) :: tra
type(state) :: s0

s0 = state(epoch(), [0._dp])

node = tranode(1._dp, [1._dp])
call assert_equal(length(node), 1, __LINE__)
call add_node(node, 2._dp, [2._dp])
call assert_equal(length(node), 2, __LINE__)
call add_node(node, 3._dp, [3._dp], head)
call add_node(head, 4._dp, [4._dp], head)
call assert_equal(length(node), 4, __LINE__)
tra = save_trajectory(s0, node, ["x"])
! call assert_almost_equal(tra%tindex, [1._dp, 2._dp, 3._dp, 4._dp], __LINE__)
call assert_almost_equal(tra%vectors(:, 1), [1._dp, 2._dp, 3._dp, 4._dp], __LINE__)
end program testtrajectories
