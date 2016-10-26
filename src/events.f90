!
! Copyright (c) 2016 Helge Eichhorn
!
! This Source Code Form is subject to the terms of the Mozilla Public
! License, v. 2.0. If a copy of the MPL was not distributed with this
! file, You can obtain one at http://mozilla.org/MPL/2.0/.
!
module events

use containers, only: parameters
use types, only: dp

implicit none

private

public :: event, update, discontinuity

type, abstract :: event
    real(dp), allocatable :: t
contains
    procedure(abstract_haspassed), deferred :: haspassed
    procedure(abstract_detect), deferred :: detect
end type event

abstract interface
    function abstract_haspassed(this, told, t, yold, y, p) result(res)
        import :: event, dp, parameters
        class(event), intent(in) :: this
        real(dp), intent(in) :: told
        real(dp), intent(in) :: t
        real(dp), intent(in) :: yold
        real(dp), intent(in) :: y
        type(parameters), intent(in) :: p
        logical :: res
    end function abstract_haspassed

    function abstract_detect(this, t, p) result(res)
        import :: event, dp, parameters
        class(event), intent(inout) :: this
        real(dp), intent(in) :: t
        type(parameters), intent(in) :: p
    end function abstract_detect
end interface

type, abstract :: update
end type update

type discontinuity
    class(event), allocatable :: event
    class(update), allocatable :: update
end type discontinuity


end module events
