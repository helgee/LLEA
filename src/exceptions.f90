!
! Copyright (c) 2016: Helge Eichhorn
!
! This Source Code Form is subject to the terms of the Mozilla Public
! License, v. 2.0. If a copy of the MPL was not distributed with this
! file, You can obtain one at http://mozilla.org/MPL/2.0/.
!
module exceptions

use iso_fortran_env, only: error_unit

implicit none

integer, parameter :: messagelen = 80
integer, parameter :: flen = 20
integer, parameter :: filelen = 65
integer, parameter :: llen = 6

type exception
    character(len=messagelen) :: message = ""
    character(len=flen), dimension(:), allocatable :: func
    character(len=filelen), dimension(:), allocatable :: file
    integer, dimension(:), allocatable :: line
end type exception

private

public :: exception, error, catch, raise, iserror

contains

logical function iserror(err)
    type(exception), intent(in) :: err
    iserror = len_trim(err%message) /= 0
end function iserror

function error(message, func, file, line) result(err)
    character(len=*), intent(in) :: message
    character(len=*), intent(in) :: func
    character(len=*), intent(in) :: file
    integer, intent(in) :: line
    type(exception) :: err

    allocate(err%func(1))
    allocate(err%file(1))
    allocate(err%line(1))
    err%message = message
    err%func = func
    err%file = file
    err%line = line
end function error

subroutine catch(err, func, file, line)
    type(exception), intent(inout) :: err
    character(len=*), intent(in) :: func
    character(len=*), intent(in) :: file
    integer, intent(in) :: line

    err%func = [func // repeat(" ", flen - len(func)), err%func]
    err%file = [file // repeat(" ", filelen - len(file)), err%file]
    err%line = [line, err%line]
end subroutine catch

subroutine raise(err)
    type(exception), intent(in) :: err

    integer :: i
    character(len=2) :: ffmt
    character(len=2) :: filefmt
    character(len=1) :: lfmt
    character(len=flen) :: routine
    character(len=filelen) :: file
    character(len=llen) :: line
    character(len=16) :: fmt

    routine = "Routine"
    file = "File"
    line = "Line"

    write(ffmt,"(i2.2)") flen
    write(filefmt,"(i2.2)") filelen
    write(lfmt,"(i1.1)") llen
    write(fmt,"(7(a))") "(a", ffmt, ",a", filefmt, ",a", lfmt, ")"

    write(error_unit,fmt) routine, file, line
    do i = 1, size(err%func)
        write(line,"(i6)") err%line(i)
        write(error_unit, fmt) err%func(i), err%file(i), adjustl(line)
    end do
    write(error_unit,*) "ERROR: "//trim(err%message)
    stop 1
end subroutine raise

end module exceptions
