module interfaces

use propagators, only: state_propagator, trajectory_propagator
use trajectories, only: state_trajectory_dp, state_trajectory_epd

implicit none

private

public :: getstate, gettrajectory

interface getstate
    module procedure state_propagator
    module procedure state_trajectory_dp
    module procedure state_trajectory_epd
end interface getstate

interface gettrajectory
    module procedure trajectory_propagator
end interface gettrajectory

end module interfaces
