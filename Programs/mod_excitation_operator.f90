module mod_excitation_operator
use mod_potentials
use mod_matrix_elements
use mod_parameters
use HO_basis
use mod_potentials
implicit none

real(8),allocatable::TransMatrix(:,:) ! This is the matrix that causes the G.S to be excited by an operator
real(8),allocatable:: gamma0(:) ! |gamma0> = T[N,N]|0>

!real(8),allocatable::Hgauss(:,:)
!complex(8),allocatable::Hcomplex(:,:)
complex(8),allocatable:: psi_tilde(:) ! The Lorentz vector
real(8),allocatable:: psi_gauss(:) ! The Gaussian vector

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

allocate(gamma0(Nmax+1))

gamma0(:) = 0.d0

	do j=1,Nmax+1
	
	   s=0.d0
	   
	   do i=1,Nmax+1
	   s = s+Operator_n1n2(i,j,p)*psi0(i)
	   end do
	   
	   gamma0(j) = s
	end do


end subroutine

! Calculates the overlap integral |<f|x^p|0>|^2
real(8) function calculate_Overlap(n)
implicit none
integer::n
integer::i,j 
real(8)::s

s = 0.d0
do i=1,Nmax+1
s = s + U(i,n)*gamma0(i)  
end do

s=s**2

calculate_Overlap = s

end function


! <psi_tilde|psi_tilde>*gamma^2/pi
real(8) function build_H_matrix_complex(W,gam)
implicit none
Complex(8),allocatable:: G(:,:)
real(8)::W,gam
integer::i,j
real(8)::s,pi
Complex(8),allocatable:: Hcomplex(:,:)

allocate(G(Nmax+1,Nmax+1))
allocate(Hcomplex(Nmax+1,Nmax+1))
allocate(psi_tilde(Nmax+1))

G(:,:)=0.d0

do i=1,Nmax+1
do j=1,Nmax+1
    if(i==j) then
      ! G(i,j) = ((Eig(i+1)-E0-W)- DCMPLX(0.d0, gam))**(-1)
       G(i,j) = ((Eig(i)-W)- DCMPLX(0.d0, gam))**(-1)
    end if
end do
end do

! initialize Hcomplex
Hcomplex(:,:) = 0.d0

!========================================
! H = UD U_dagger
!Hcomplex = MATMUL(U,MATMUL(G,Udagger))
!========================================

! This is the complex matrix, now multiply against oper|0>
Hcomplex = MATMUL(U,MATMUL(G,Udagger))

psi_tilde = MATMUL(Hcomplex,gamma0)

! Calculate the psi_tilde norm
s = 0.d0
do i =1, Nmax+1
s = s+REALPART(psi_tilde(i)*dCONJG(psi_tilde(i)))
!s = s+psi_tilde(i)**2 !dCONJG(psi_tilde(i))
end do

!print *, '<psi tilde| psi tilde>', s


!call exit()

pi = datan(1.d0)*4.d0


build_H_matrix_complex =s*(gam)*(1.d0/pi)

deallocate(G)
deallocate(Hcomplex)
deallocate(psi_tilde)

end function


! <psi_tilde|psi_tilde>*gamma^2/pi
! Gaussian kernel
real(8) function build_H_matrix_gauss(W,sigma)
implicit none
real(8),allocatable:: G(:,:)
real(8)::W,sigma
integer::i,j
real(8)::s,pi
real(8),allocatable:: Hgauss(:,:)

allocate(G(Nmax+1,Nmax+1))
allocate(Hgauss(Nmax+1,Nmax+1))
allocate(psi_gauss(Nmax+1))

G(:,:)=0.d0

do i=1,Nmax+1
do j=1,Nmax+1
    if(i==j) then
      ! G(i,j) = dexp(-0.25d0*(1.d0/(sigma**2))*(Eig(i+1)-E0-W)**2)
       G(i,j) = dexp(-0.25d0*(1.d0/(sigma**2))*(Eig(i)-W)**2)
    end if
end do
end do

Hgauss = MATMUL(U,MATMUL(G,Udagger))
psi_gauss = MATMUL(Hgauss,gamma0) ! |psi'> = Exp{} Oper|0>

! Calculate the psi_tilde norm
s = 0.d0
do i =1, Nmax+1
s = s+psi_gauss(i)**2
end do

pi = datan(1.d0)*4.d0


build_H_matrix_gauss = s*(1.d0/dsqrt(2.d0*pi*sigma**2))

deallocate(G)
deallocate(Hgauss)
deallocate(psi_gauss)

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
!a = (Eig(n)-Eig(1))-de
!b = (Eig(n)-Eig(1))+de
a = Eig(n)-de
b = Eig(n)+de
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


! This function integrates ! <psi_tilde|psi_tilde>*gamma^2/pi
! int_{0}^{1} f(x)dx 
! int_{a}^{b} f(x) dx = \int_{0}^{1}f(m*x+a) dx
! n=1,Nmax+1
real(8) function integrate_response_gauss(n,sigma,de)
implicit none
integer:: n ! the eigen value of interest
real(8):: sigma ! The sigma parameter 
real(8):: de ! The epsilon paramter for the integration
real(8)::a,b,y0
integer::i
real(8)::s,yi,m,wi

!a = (Eig(n)-Eig(1))-de
!b = (Eig(n)-Eig(1))+de
a = Eig(n)-de
b = Eig(n)+de
y0 = 0.5d0*(b+a)
m = 0.5d0*(b-a)

s=0.d0
do i=1, NquadPosition
yi = dx(i)
wi = m*yi+y0
s= s+build_H_matrix_gauss(wi,sigma)*dwX(i)*((b-a)*0.5d0)
end do

integrate_response_gauss = s
end function

end module
