module interfaces

use propagators, only: state_propagator, trajectory_propagator

implicit none

private

public :: getstate, gettrajectory

interface getstate
    module procedure state_propagator
end interface getstate

interface gettrajectory
    module procedure trajectory_propagator
end interface gettrajectory

end module interfaces
