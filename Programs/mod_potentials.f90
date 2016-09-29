! This module contains the potentials that we will use
module mod_potentials
real(8),parameter:: L0 = 20.d0 ! fm
real(8),parameter:: V0 = 20.0 ! 10.d0 !200.d0 ! MeV

contains


! The square well box, potential
real(8) function Vbox(x)
real(8)::x

Vbox = 0.d0

if((x >= -L0).and.(x <= L0)) then
Vbox = -V0
end if

!Vbox = V0

!if((x >= -L0).and.(x <= L0)) then
!Vbox = 0.d0
!end if


end function


end module
