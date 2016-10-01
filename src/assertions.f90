module assertions

use iso_fortran_env, only: error_unit
use types, only: dp
use exceptions, only: exception
use util, only: stop_error
use math, only: eps

implicit none

private

public :: assert, assert_false, assert_equal, assert_not_equal,&
    assert_almost_equal, assert_raises

interface assert_equal
    module procedure assert_equal_character
    module procedure assert_equal_scalar_int
    module procedure assert_equal_vector_int
    module procedure assert_equal_scalar_dp
    module procedure assert_equal_vector_dp
    module procedure assert_equal_matrix_dp
end interface

interface assert_not_equal
    module procedure assert_not_equal_scalar_dp
    module procedure assert_not_equal_vector_dp
end interface

interface assert_almost_equal
    module procedure assert_almost_equal_scalar
    module procedure assert_almost_equal_vector
end interface

interface assert
    module procedure assert_scalar
    module procedure assert_vector
end interface

interface assert_false
    module procedure assert_false_scalar
    module procedure assert_false_vector
end interface

contains

subroutine assert_raises(message, err, line)
    character(len=*), intent(in) :: message
    type(exception), intent(in) :: err
    integer, intent(in) :: line

    character(len=20) :: c

    if (trim(message) /= trim(err%message)) then
        c = " Expected message:"
        write(error_unit,"(a,a)") c, trim(message)
        c = " Actual message:"
        write(error_unit,"(a,a)") c, trim(err%message)
        write(error_unit,*) "Line:", line
        call stop_error("Assertion failed.")
    end if
end subroutine assert_raises

subroutine assert_equal_character(a, b, line)
    character(len=*), intent(in) :: a
    character(len=*), intent(in) :: b
    integer, intent(in) :: line
    if (a /= b) then
        write(error_unit,*) a, " != ", b
        write(error_unit,*) "Line:", line
        call stop_error("Assertion failed.")
    end if
end subroutine assert_equal_character

subroutine assert_equal_scalar_int(a, b, line)
    integer, intent(in) :: a
    integer, intent(in) :: b
    integer, intent(in) :: line
    if (a /= b) then
        write(error_unit,*) a, "!=", b
        write(error_unit,*) "Line:", line
        call stop_error("Assertion failed.")
    end if
end subroutine assert_equal_scalar_int

subroutine assert_equal_vector_int(a, b, line)
    integer, dimension(:), intent(in) :: a
    integer, dimension(:), intent(in) :: b
    integer, intent(in) :: line
    integer :: i
    if (size(a) /= size(b)) then
        write(error_unit,*) "size(a)=",size(a),"!= size(b)=",size(b)
        write(error_unit,*) "Line:", line
        call stop_error("Assertion failed.")
    end if
    do i=1,size(a)
        call assert_equal_scalar_int(a(i), b(i), line)
    end do
end subroutine assert_equal_vector_int

subroutine assert_equal_scalar_dp(a, b, line)
    real(dp), intent(in) :: a
    real(dp), intent(in) :: b
    integer, intent(in) :: line
    if (a /= b) then
        write(error_unit,*) a, "!=", b
        write(error_unit,*) "Line:", line
        call stop_error("Assertion failed.")
    end if
end subroutine assert_equal_scalar_dp

subroutine assert_equal_vector_dp(a, b, line)
    real(dp), dimension(:), intent(in) :: a
    real(dp), dimension(:), intent(in) :: b
    integer, intent(in) :: line
    integer :: i
    if (size(a) /= size(b)) then
        write(error_unit,*) "size(a)=",size(a),"!= size(b)=",size(b)
        write(error_unit,*) "Line:", line
        call stop_error("Assertion failed.")
    end if
    do i=1,size(a)
        call assert_equal_scalar_dp(a(i), b(i), line)
    end do
end subroutine assert_equal_vector_dp

subroutine assert_equal_matrix_dp(a, b, line)
    real(dp), dimension(:,:), intent(in) :: a
    real(dp), dimension(:,:), intent(in) :: b
    integer, intent(in) :: line
    integer, dimension(2) :: shp
    integer :: i
    integer :: j
    if (size(a) /= size(b)) then
        write(error_unit,*) "size(a)=",size(a),"!= size(b)=",size(b)
        write(error_unit,*) "Line:", line
        call stop_error("Assertion failed.")
    end if
    shp = shape(a)
    do i=1,shp(2)
        do j=1,shp(1)
            call assert_equal_scalar_dp(a(j, i), b(j, i), line)
        end do
    end do
end subroutine assert_equal_matrix_dp

subroutine assert_not_equal_scalar_dp(a, b, line)
    real(dp), intent(in) :: a
    real(dp), intent(in) :: b
    integer, intent(in) :: line
    if (a == b) then
        write(error_unit,*) a, "==", b
        write(error_unit,*) "Line:", line
        call stop_error("Assertion failed.")
    end if
end subroutine assert_not_equal_scalar_dp

subroutine assert_not_equal_vector_dp(a, b, line)
    real(dp), dimension(:), intent(in) :: a
    real(dp), dimension(:), intent(in) :: b
    integer, intent(in) :: line
    integer :: i
    if (size(a) /= size(b)) then
        write(error_unit,*) "size(a)=",size(a),"!= size(b)=",size(b)
        write(error_unit,*) "Line:", line
        call stop_error("Assertion failed.")
    end if
    do i=1,size(a)
        call assert_not_equal_scalar_dp(a(i), b(i), line)
    end do
end subroutine assert_not_equal_vector_dp

subroutine assert_scalar(a, line)
    logical, intent(in) :: a
    integer, intent(in) :: line
    if (a .neqv. .true.) then
        write(error_unit,*) "Line:", line
        call stop_error("Assertion failed.")
    end if
end subroutine assert_scalar

subroutine assert_false_scalar(a, line)
    logical, intent(in) :: a
    integer, intent(in) :: line
    if (a .neqv. .false.) then
        write(error_unit,*) "Line:", line
        call stop_error("Assertion failed.")
    end if
end subroutine assert_false_scalar

subroutine assert_vector(a, line)
    logical, dimension(:), intent(in) :: a
    integer, intent(in) :: line
    if (any(a .neqv. .true.)) then
        write(error_unit,*) "Line:", line
        call stop_error("Assertion failed.")
    end if
end subroutine assert_vector

subroutine assert_false_vector(a, line)
    logical, dimension(:), intent(in) :: a
    integer, intent(in) :: line
    if (any(a .neqv. .false.)) then
        write(error_unit,*) "Line:", line
        call stop_error("Assertion failed.")
    end if
end subroutine assert_false_vector

subroutine assert_almost_equal_scalar(a, b, line, eps)
    ! Source: http://floating-point-gui.de/errors/comparison/
    real(dp), intent(in) :: a
    real(dp), intent(in) :: b
    integer, intent(in) :: line
    real(dp), intent(in), optional :: eps
    real(dp) :: eps_
    real(dp) :: abs_a
    real(dp) :: abs_b
    real(dp) :: diff
    real(dp) :: min_normal

    eps_ = sqrt(epsilon(0._dp))
    if (present(eps)) then
        eps_ = eps
    end if

    min_normal = tiny(0._dp)
    abs_a = abs(a)
    abs_b = abs(b)
    diff = abs(a-b)
    if (a == b) then
        return
    ! TODO: This does not seem to work.
    ! else if ((a == 0._dp) .or. (b == 0._dp) .or. (diff < min_normal)) then
    !     if (diff > (eps_*min_normal)) then
    !         write(error_unit,*) "a=",a,"!~= b=",b
    !         call stop_error("Assertion failed.")
    !     end if
    else
        if (diff > eps_) then
            write(error_unit,*) "abs(a=",abs_a,"-b=",b,") > eps=",eps_
            write(error_unit,*) "Line:", line
            call stop_error("Assertion failed.")
        end if
        ! if ((diff / (abs_a + abs_b)) > eps_) then
        !     write(error_unit,*) "diff=",diff,"/(abs(a)=",abs_a,"+abs(b)=",abs_b,") > eps=",eps_
        !     call stop_error("Assertion failed.")
        ! end if
    end if
end subroutine assert_almost_equal_scalar

subroutine assert_almost_equal_vector(a, b, line, eps)
    real(dp), dimension(:), intent(in) :: a
    real(dp), dimension(:), intent(in) :: b
    integer, intent(in) :: line
    real(dp), intent(in), optional :: eps
    integer :: i
    real(dp) :: eps_

    eps_ = sqrt(epsilon(0._dp))
    if (present(eps)) then
        eps_ = eps
    end if

    if (size(a) /= size(b)) then
        write(error_unit,*) "size(a)=",size(a),"!= size(b)=",size(b)
        write(error_unit,*) "Line:", line
        call stop_error("Assertion failed.")
    end if
    do i=1,size(a)
        call assert_almost_equal_scalar(a(i), b(i), line, eps_)
    end do
end subroutine assert_almost_equal_vector

end module assertions
