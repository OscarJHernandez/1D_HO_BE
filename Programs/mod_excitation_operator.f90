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

G(:,:)=0.d0

do i=0,Nmax
do j=0,Nmax
    if(i==j) then
       G(i,j) = ((Eig(i+1)-E0-W)- DCMPLX(0.d0, gam))**(-1)
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


build_H_matrix_complex =s*(gam)*(1.d0/pi)

deallocate(G)
deallocate(Hcomplex)
deallocate(psi_tilde)

end function

! This function integrates ! <psi_tilde|psi_tilde>*gamma^2/pi
! int_{0}^{1} f(x)dx 
! int_{a}^{b} f(x) dx = \int_{0}^{1}f(m*x+a) dx
! n=1,Nmax+1
real(8) function integrate_response(n,gam,de)
implicit none
integer:: n ! the eigen value of interest
real(8):: gam ! The gamma parameter 
real(8):: de ! The epsilon paramter for the integration
real(8)::a,b,y0
integer::i
real(8)::s,yi,m,wi

!real(8),allocatable::da(:),db(:),dxx(:),dw(:),e(:)
!real*8,parameter:: depsma=1.0d-18
!allocate(da(NquadPosition),dbX(NquadPosition),dx(NquadPosition),dwX(NquadPosition),eX(NquadPosition))
a = (Eig(n)-Eig(1))-de
b = (Eig(n)-Eig(1))+de
y0 = 0.5d0*(b+a)
m = 0.5d0*(b-a)

!print *, 'w0',(Eig(n)-Eig(1))

s=0.d0
do i=1, NquadPosition
yi = dx(i)
wi = m*yi+y0
s= s+build_H_matrix_complex(wi,gam)*dwX(i)*((b-a)*0.5d0)
end do

integrate_response = s
end function


end module
