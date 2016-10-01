module epochs

use types, only: dp

implicit none

type epoch
    real(dp) :: jd
    real(dp) :: jd1
    character(len=3) :: scale
end type epoch

end module epochs
