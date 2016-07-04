! This module contains the potentials that we will use
module mod_potentials
real(8),parameter:: boxLength =1.d0 ! fm

contains


! The square well box, potential
real(8) function Vbox(x)
real(8)::x
real(8),parameter:: V0 = -10.d0 ! MeV

Vbox = V0

if(x.ge.1.d0) then
Vbox = 0.d0
end if

end function


end module
