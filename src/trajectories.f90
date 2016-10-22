module trajectories

use states, only: state
use types, only: dp

implicit none

integer, parameter :: fieldlen = 12

type listnode
    type(listnode), pointer :: next => null()
    real(dp), dimension(:), allocatable :: val
end type listnode

type trajectory
    type(state) :: initial_state
    type(state) :: final_state
    character(len=fieldlen), dimension(:), allocatable :: fields
    real(dp), dimension(:,:), allocatable :: vectors
end type trajectory

end module trajectories
