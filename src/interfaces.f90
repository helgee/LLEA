module interfaces

use states, only: init_state, state
use propagators, only: state_propagator, trajectory_propagator
use trajectories, only: trajectory, init_trajectory, init_trajectory_array

implicit none

private

public :: state, trajectory

interface state
    module procedure init_state
    module procedure state_propagator
end interface state

interface trajectory
    module procedure init_trajectory_array
    module procedure init_trajectory
    module procedure trajectory_propagator
end interface trajectory

end module interfaces
