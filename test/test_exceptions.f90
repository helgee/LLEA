program testexceptions

use exceptions
use assertions

implicit none

type(exception) :: err

! Works:
! call faulty2()

call faulty1(err)
call catch(err, "main", __FILE__, __LINE__)

! Works:
! call raise(err)

call assert_equal(err%line, [14, 27, 35], __LINE__)

contains

subroutine faulty1(err)
    type(exception), intent(inout) :: err

    call faulty2(err)
    call catch(err, "faulty1", __FILE__, __LINE__)
end subroutine faulty1

subroutine faulty2(err)
    type(exception), intent(inout), optional :: err

    type(exception) :: err_

    err_ = error("Everything is broken.", "faulty2", __FILE__, __LINE__)
    if (present(err)) then
        err = err_
        return
    else
        call raise(err_)
    end if
end subroutine faulty2

end program testexceptions
