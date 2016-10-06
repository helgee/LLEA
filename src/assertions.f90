module assertions

use iso_fortran_env, only: error_unit
use types, only: dp
use exceptions, only: exception
use util, only: stop_error
use math, only: isapprox, eps

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
    module procedure assert_almost_equal_matrix
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

subroutine assert_almost_equal_scalar(a, b, line, rtol, atol)
    real(dp), intent(in) :: a
    real(dp), intent(in) :: b
    integer, intent(in) :: line
    real(dp), intent(in), optional :: rtol
    real(dp), intent(in), optional :: atol

    real(dp) :: atol_
    real(dp) :: rtol_

    rtol_ = sqrt(eps)
    if (present(rtol)) rtol_ = rtol
    atol_ = 0._dp
    if (present(atol)) atol_ = atol

    if (.not.isapprox(a, b, rtol_, atol_)) then
        write(error_unit,*) "abs(a=",a,"-b=",b,") > tolerance"
        write(error_unit,*) "Line:", line
        call stop_error("Assertion failed.")
    end if
end subroutine assert_almost_equal_scalar

subroutine assert_almost_equal_vector(a, b, line, rtol, atol, elementwise)
    real(dp), dimension(:), intent(in) :: a
    real(dp), dimension(:), intent(in) :: b
    integer, intent(in) :: line
    real(dp), intent(in), optional :: rtol
    real(dp), intent(in), optional :: atol
    logical, intent(in), optional :: elementwise

    real(dp) :: atol_
    real(dp) :: rtol_
    logical :: elementwise_
    integer :: i

    rtol_ = sqrt(eps)
    if (present(rtol)) rtol_ = rtol
    atol_ = 0._dp
    if (present(atol)) atol_ = atol
    elementwise_ = .false.
    if (present(elementwise)) elementwise_ = elementwise

    if (size(a) /= size(b)) then
        write(error_unit,*) "size(a)=",size(a),"!= size(b)=",size(b)
        write(error_unit,*) "Line:", line
        call stop_error("Assertion failed.")
    end if

    if (.not.elementwise_) then
        if (.not.isapprox(a, b, rtol_, atol_)) then
            write(error_unit,*) "a = ", a
            write(error_unit,*) "b = ", b
            write(error_unit,*) "norm(a - b) > tolerance"
            write(error_unit,*) "Line:", line
            call stop_error("Assertion failed.")
        end if
    else
        do i=1,size(a)
            call assert_almost_equal_scalar(a(i), b(i), line, rtol_, atol_)
        end do
    end if
end subroutine assert_almost_equal_vector

subroutine assert_almost_equal_matrix(a, b, line, atol, rtol, elementwise)
    real(dp), dimension(:,:), intent(in) :: a
    real(dp), dimension(:,:), intent(in) :: b
    integer, intent(in) :: line
    real(dp), intent(in), optional :: atol
    real(dp), intent(in), optional :: rtol
    logical, intent(in), optional :: elementwise
    integer, dimension(2) :: shp

    real(dp) :: atol_
    real(dp) :: rtol_
    logical :: elementwise_
    integer :: i
    integer :: j

    atol_ = 0._dp
    if (present(atol)) atol_ = atol
    rtol_ = sqrt(eps)
    if (present(rtol)) rtol_ = rtol
    elementwise_ = .false.
    if (present(elementwise)) elementwise_ = elementwise

    if (size(a) /= size(b)) then
        write(error_unit,*) "size(a)=",size(a),"!= size(b)=",size(b)
        write(error_unit,*) "Line:", line
        call stop_error("Assertion failed.")
    end if

    if (.not.elementwise_) then
        if (.not.isapprox(a, b, rtol_, atol_)) then
            write(error_unit,*) "a = ", a
            write(error_unit,*) "b = ", b
            write(error_unit,*) "norm(a - b) > tolerance"
            write(error_unit,*) "Line:", line
            call stop_error("Assertion failed.")
        end if
    else
        shp = shape(a)
        do i=1,shp(2)
            do j=1,shp(1)
                call assert_almost_equal_scalar(a(j, i), b(j, i), line, atol_, rtol_)
            end do
        end do
    end if
end subroutine assert_almost_equal_matrix

end module assertions
