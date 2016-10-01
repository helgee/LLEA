module states

use types, only: dp
use bodies, only: body
use epochs, only: epoch

implicit none

type state
    type(epoch) :: t
    real(dp), dimension(3) :: r
    real(dp), dimension(3) :: v
    character(len=8) :: frame
    type(body) :: center
end type state

end module states
