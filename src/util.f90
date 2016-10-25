!
! Copyright (c) 2016 Helge Eichhorn
!
! This Source Code Form is subject to the terms of the Mozilla Public
! License, v. 2.0. If a copy of the MPL was not distributed with this
! file, You can obtain one at http://mozilla.org/MPL/2.0/.
!
module util

use, intrinsic :: iso_fortran_env, only: error_unit
use types, only: dp

implicit none

#if __unix__
character(len=1), parameter :: sep='/'
#elif _WIN32
character(len=1), parameter :: sep='\'
#else
character(len=1), parameter :: sep='/'
#endif

interface joinpath
    module procedure joinpath_1
    module procedure joinpath_2
    module procedure joinpath_3
end interface joinpath

private

public :: stop_error, printbool, joinpath, sep, mode_open,&
    timestamp_file, uppercase, projectdir

contains

function projectdir(sub) result(p)
    character(len=*), intent(in), optional :: sub
    character(len=:), allocatable :: p

    integer :: ind

    p = __FILE__
    ind = index(p(1:index(p, sep, .true.) - 1), sep, .true.) - 1
    if (present(sub)) then
        p = joinpath(p(1:ind), sub)
    end if
end function projectdir

function uppercase(str) result(upper)
    character(len=*), intent(in) :: str
    character(len=len(str)) :: upper
    integer :: j
    do j = 1,len(str)
        if(str(j:j) >= "a" .and. str(j:j) <= "z") then
           upper(j:j) = achar(iachar(str(j:j)) - 32)
        else
           upper(j:j) = str(j:j)
        end if
    end do
end function uppercase

function timestamp_file(file, ext) result(f)
    character(len=*), intent(in) :: file
    character(len=*), intent(in), optional :: ext
    character(len=:), allocatable :: f

    integer :: i
    integer :: l
    logical :: exists
    character(len=8) :: date
    character(len=10) :: time
    character(len=14) :: timestamp
    character(len=:), allocatable :: root
    character(len=:), allocatable :: ext_
    character(len=1) :: ind

    l = index(trim(file), '.', .true.)
    root = file(:l-1)
    if (present(ext)) then
        ext_ = ext
    else
        ext_ = file(l+1:len_trim(file))
    end if
    call date_and_time(date, time)
    timestamp = date//time(:6)

    l = len(root) + len(ext_) + 16
    allocate(character(len=l) :: f)
    f = root//'_'//timestamp//'.'//ext_
    inquire(file=f, exist=exists)
    if (exists) then
        i = 1
        l = l + 2
        do
            write(ind,'(i1.1)') i
            deallocate(f)
            allocate(character(len=l) :: f)
            f = root//'_'//timestamp//'_'//ind//'.'//ext_
            inquire(file=f, exist=exists)
            if (.not.exists) exit
            i = i + 1
        end do
    end if
end function timestamp_file

function mode_open(filename, mode) result(u)
    character(len=*), intent(in) :: filename
    character(len=1), intent(in) :: mode
    integer :: u

    if (mode == "w") then
        open(newunit=u, file=filename, status="replace", action="write")
    else if (mode == "a") then
        open(newunit=u, file=filename, status="old", position="append")
    else if (mode == "r") then
        open(newunit=u, file=filename, status="old", action="read")
    end if
end function mode_open

logical function wrong_sep(s)
    character(len=1), intent(in) :: s
    if (s == '/' .or. s == '\') then
        wrong_sep = s /= sep
    else
        wrong_sep = .false.
    end if
end function wrong_sep

function clean_sep(s0) result(s)
    character(len=*), intent(in) :: s0
    character(len=:), allocatable :: s
    integer :: i

    s = trim(s0)
    do i = 1, len(s)
        if (wrong_sep(s(i:i))) s(i:i) = sep
    end do
end function clean_sep

function len_sep(s) result(l)
    character(len=*), intent(in) :: s
    integer :: l

    l = len_trim(s)
    if (s(l:l) == sep) then
        l = l - 1
    end if
end function len_sep

subroutine insert_path(p, s, start)
    character(len=*), intent(inout) :: p
    character(len=*), intent(in) :: s
    integer, intent(inout) :: start

    integer :: l

    l = len_sep(s)
    p(start:start+l-1) = s(:l)
    start = start + l
    if (start < len(p)) then
        p(start:start) = sep
        start = start + 1
    end if
end subroutine insert_path

function joinpath_1(part0, part1) result(p)
    character(len=*), intent(in) :: part0
    character(len=*), intent(in) :: part1
    character(len=:), allocatable :: p

    character(len=:), allocatable :: part0_
    character(len=:), allocatable :: part1_
    integer :: n
    integer :: m

    part0_ = clean_sep(part0)
    part1_ = clean_sep(part1)

    n = len_sep(part0_) + len_sep(part1_) + 1

    allocate(character(len=n) :: p)
    m = 1
    call insert_path(p, part0_, m)
    call insert_path(p, part1_, m)
end function joinpath_1

function joinpath_2(part0, part1, part2) result(p)
    character(len=*), intent(in) :: part0
    character(len=*), intent(in) :: part1
    character(len=*), intent(in) :: part2
    character(len=:), allocatable :: p

    integer :: n
    integer :: m
    character(len=:), allocatable :: part0_
    character(len=:), allocatable :: part1_
    character(len=:), allocatable :: part2_

    part0_ = clean_sep(part0)
    part1_ = clean_sep(part1)
    part2_ = clean_sep(part2)

    n = len_sep(part0_) + len_sep(part1_)&
       + len_sep(part2_) + 2

    allocate(character(len=n) :: p)
    m = 1
    call insert_path(p, part0_, m)
    call insert_path(p, part1_, m)
    call insert_path(p, part2_, m)
end function joinpath_2

function joinpath_3(part0, part1, part2, part3) result(p)
    character(len=*), intent(in) :: part0
    character(len=*), intent(in) :: part1
    character(len=*), intent(in) :: part2
    character(len=*), intent(in) :: part3
    character(len=:), allocatable :: p

    integer :: n
    integer :: m
    character(len=:), allocatable :: part0_
    character(len=:), allocatable :: part1_
    character(len=:), allocatable :: part2_
    character(len=:), allocatable :: part3_

    part0_ = clean_sep(part0)
    part1_ = clean_sep(part1)
    part2_ = clean_sep(part2)
    part3_ = clean_sep(part3)

    n = len_sep(part0_) + len_sep(part1_)&
       + len_sep(part2_) + len_sep(part3_) + 3

    allocate(character(len=n) :: p)
    m = 1
    call insert_path(p, part0_, m)
    call insert_path(p, part1_, m)
    call insert_path(p, part2_, m)
    call insert_path(p, part3_, m)
end function joinpath_3

! Function: printbool
!   Print the string true/false instead of T/F.
!
! Parameter:
!   bool - Boolean variable
function printbool(bool) result(str)
    logical, intent(in) :: bool
    character(len=7) :: str
    if (bool) then
        str = " true"
    else
        str = " false"
    end if
end function printbool

! Subroutine: stop_error
!   Print an error string to stderr and stop the program with exitcode 1.
!
! Parameter:
!   str - Error string
subroutine stop_error(str, file, line)
    character(len=*), intent(in) :: str
    character(len=*), intent(in), optional :: file
    integer, intent(in), optional :: line
    if (present(file)) write(error_unit,*) "File:", file
    if (present(line)) write(error_unit,*) "Line:", line
    write(error_unit,*) str
    stop 1
end subroutine stop_error

! Subroutine: warning
!   Print an error string to stderr.
!
! Parameter:
!   str - Error string
subroutine warning(str)
    character(len=*), intent(in) :: str
    write(error_unit,*) str
end subroutine warning

end module util
