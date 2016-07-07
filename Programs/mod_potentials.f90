! This module contains the potentials that we will use
module mod_potentials
real(8),parameter:: L0 = 10.d0 ! fm

contains


! The square well box, potential
real(8) function Vbox(x)
real(8)::x
real(8),parameter:: V0 = 200.d0 ! MeV

Vbox = 0.d0

if((x >= -L0).and.(x <= L0)) then
Vbox = -V0
end if

end function


end module
