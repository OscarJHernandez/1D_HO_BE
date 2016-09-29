module mod_excitation_operator
use mod_potentials
use mod_matrix_elements
use mod_parameters
use HO_basis
use mod_potentials
implicit none

real(8),allocatable::TransMatrix(:,:) ! This is the matrix that causes the G.S to be excited by an operator
real(8),allocatable:: gamma0(:) ! |gamma0> = T[N,N]|0>

complex(8),allocatable::Hcomplex(:,:)
complex(8),allocatable:: psi_tilde(:) ! The Lorentz vector

contains

! This module calculates the transition operator from initial to final state
subroutine computeTransMatrix(p)
implicit none
integer::i,j
real(8)::p

allocate(TransMatrix(0:Nmax,0:Nmax))

TransMatrix(:,:) =0.d0

do i =0,Nmax
	do j =0,Nmax
	TransMatrix(i,j) = Operator_n1n2(i,j,p)
	end do
end do


end subroutine

! |k > = TransMatrix[Nmax,Nmax]|0 >
subroutine excite_OpPsi0(p)
implicit none
integer ::i,j
real(8)::s,p

allocate(gamma0(0:Nmax))

gamma0(:) = 0.d0

	do j=0,Nmax
	
	   s=0.d0
	   
	   do i=0,Nmax
	   s = s+Operator_n1n2(i,j,p)*psi0(i)
	   end do
	   
	   gamma0(j) = s
	end do



end subroutine

! <psi_tilde|psi_tilde>*gamma^2/pi
real(8) function build_H_matrix_complex(W,gam)
implicit none
Complex(8),allocatable:: G(:,:)
real(8)::W,gam
integer::i,j
real(8)::s,pi

allocate(G(0:Nmax,0:Nmax))
allocate(Hcomplex(0:Nmax,0:Nmax))
allocate(psi_tilde(0:Nmax))

do i=0,Nmax
do j=0,Nmax
    if(i==j) then
       G(i,j) = (Eig(i+1)-E0-W- DCMPLX(0.d0, gam))**(-1)
    end if
end do
end do

Hcomplex = MATMUL(Udagger,MATMUL(G,U))

psi_tilde = MATMUL(Hcomplex,psi0)

! Calculate the psi_tilde norm
s = 0.d0
do i =0, Nmax
s = s+REALPART(psi_tilde(i)*dCONJG(psi_tilde(i)))
end do

pi = datan(1.d0)*4.d0


build_H_matrix_complex =s*(gam**2)*(1.d0/pi)

deallocate(G)
deallocate(Hcomplex)
deallocate(psi_tilde)

end function


end module
