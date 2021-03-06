! Here we define the HO-basis functions

module HO_basis
use mod_parameters


contains

! This is the basis function in the 1D-Harmonic Oscillator Basis
real(8) function Rn(n,x)
integer::n
real(8)::x
real(8)::s1,s2

! This is the older implementation
!s1 = Norm(n)*HermitePoly(n,dsqrt(2.d0*v)*x)

! This works for very Large Nmax
s2 = HermitePolyOrtho(n,dsqrt(2.d0*v)*x)*((2.d0*v)**(0.25d0))

Rn = s2

end function

! Normalization
! (1/Sqrt[2^n * Gamma[n+1]])*((2*v/pi)^(1/4))
real(8) function Norm(n)
integer::n
real(8):: logN

logN = 0.25d0*dlog(2.d0*v/pi)-dfloat(n)*0.5d0*dlog(2.d0)-0.5d0*DLGAMA(dfloat(n+1))
Norm =  Dexp(logN)

end function


! The Hermite Polynomials
! Multiplied by the Exponential Factor Exp[-0.5*x*x]
real(8) function HermitePoly(n,x)
implicit none
integer::n,ni
real(8)::x
real(8)::um1, um2,up

um2 = 1.d0*dExp(-0.5d0*x**2)
um1 = 2.d0*x*dExp(-0.5d0*x**2)

if(n==0) then
HermitePoly = 1.d0*dExp(-0.5d0*x**2)
else if(n==1) then
HermitePoly = 2.d0*x*dExp(-0.5d0*x**2)

else

	do ni =2,n
	up = 2.d0*x*um1-2.d0*dfloat(ni-1)*um2
	um2 = um1
	um1 = up 
	end do

HermitePoly = up
end if


end function


! The Othonormal Hermite Polynomials generated by a recursion
! Multiplied by the Exponential weight
real(8) function HermitePolyOrtho(n,x)
implicit none
integer::n,ni
real(8)::x
real(8)::um0, um1,up

um1 = 0.d0
um0 = (pi**(-0.25d0))*dExp(-0.5d0*x*x)

if(n.eq.0) then 
up = um0

else
	do ni =0,n-1
	up = x*dsqrt(2.d0/dfloat(ni+1))*um0-dsqrt(dfloat(ni)/dfloat(ni+1))*um1
	um1 = um0
	um0 = up 
	end do
end if

HermitePolyOrtho = up


end function


end module
