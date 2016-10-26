!
! Copyright (c) 2016 Helge Eichhorn
!
! This Source Code Form is subject to the terms of the Mozilla Public
! License, v. 2.0. If a copy of the MPL was not distributed with this
! file, You can obtain one at http://mozilla.org/MPL/2.0/.
!
module events

use types, only: dp

implicit none

private

public :: event, update, discontinuity

type, abstract :: event
end type event

type, abstract :: update
end type update

type discontinuity
    class(event), allocatable :: event
    class(update), allocatable :: update
end type discontinuity


end module events
