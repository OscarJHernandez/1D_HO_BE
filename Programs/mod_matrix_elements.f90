module mod_matrix_elements
use mod_potentials
use mod_parameters
use HO_basis
use mod_potentials
implicit none


!real(8),allocatable::daP(:),dbP(:),dp(:),dwP(:),eP(:)
real(8),allocatable::daX(:),dbX(:),dx(:),dwX(:),eX(:)
real(8),allocatable::daY(:),dbY(:),dy(:),dwY(:),eY(:)

!integer,parameter:: ipolyMomentum=2 !  Legendre polynomial from 0 to 1 
integer,parameter:: ipolyPosition=1 !  Legendre polynomial from -1 to 1 
integer,parameter:: ipoly_Y=2 !  Legendre polynomial from 0 to 1 
real*8,parameter:: depsma=1.0d-18
real*8:: al,ierr,iderr,dbe
real(8),allocatable:: Rn_xi(:,:) ! (n,x)
real(8),allocatable:: H(:,:) ! The Hamiltonian Matrix
real(8),allocatable:: U(:,:) ! This matrix contains the Eigenvalues of H
real(8),allocatable:: Udagger(:,:) ! The Transpose of the U matrix
real(8),allocatable:: Eig(:) ! The Eigenvalues of te Matrix
real(8),allocatable:: E0 ! The ground state energy
real(8),allocatable:: psi0(:) ! The ground state Eigenvector
integer,parameter:: steps =1000  !! This is the number of steps used for the wavefunction


contains


! Initialize the matrix to store the Basis function at mesh-points
subroutine initMesh()
implicit none
integer:: i,ni
real(8)::xi

! Allocate the Gauss-Legendre quadrature routine
    !allocate(daP(NquadMomentum),dbP(NquadMomentum),dp(NquadMomentum),dwP(NquadMomentum),eP(NquadMomentum))
    allocate(daX(NquadPosition),dbX(NquadPosition),dx(NquadPosition),dwX(NquadPosition),eX(NquadPosition))
    allocate(daY(NquadPosition),dbY(NquadPosition),dY(NquadPosition),dwY(NquadPosition),eY(NquadPosition))
    
    al = 0.d0
    
    call drecur(NquadPosition,ipoly_Y,al,dbe,daY,dbY,iderr)
    call dgauss(NquadPosition,daY,dbY,depsma,dy,dwY,ierr,eY)
    
    call drecur(NquadPosition,ipolyPosition,al,dbe,daX,dbX,iderr)
    call dgauss(NquadPosition,daX,dbX,depsma,dx,dwX,ierr,eX)

allocate(Rn_xi(0:Nmax,NquadPosition))

	do i=1, NquadPosition
		do ni=0, Nmax
			xi = Xmax*dx(i)
			!xi = dx(i)
			Rn_xi(ni,i) = Rn(ni,xi)
		end do
	end do



end subroutine


! The Kinetic Energy matrix element
! <n1| T |n2>
real(8) function KineticEnergy(n1,n2)
implicit none
integer::n1,n2
real(8)::s1,s2,s3

s1=0.d0
s2=0.d0
s3=0.d0

if(n1.eq.n2) then
s1 =0.5d0*hw*(dfloat(n1)+0.5d0)
end if

if(n1.eq.(n2-2)) then
s2 = -0.25d0*hw*dsqrt(dfloat(n2*(n2-1)))
end if

if(n1.eq.(n2+2)) then
s3 = -0.25d0*hw*dsqrt(dfloat((n2+1)*(n2+2)))
end if


KineticEnergy = s1+s2+s3

end function


! This is the potential matrix element
! <n1|V|n2>
real(8) function Vn1n2(n1,n2)
implicit none
integer::n1,n2
real(8)::f,xi
real(8)::s
integer::i

s =0.d0
do i=1,NquadPosition
xi = Xmax*dx(i)


f = Rn_xi(n1,i)*Rn_xi(n2,i)*Vbox(xi)

! The harmonic Potential for Testing purposes
!f = Rn_xi(n1,i)*Rn_xi(n2,i)*hw*v*xi*xi
s = s +f*dwX(i)*Xmax
end do

!print *, n1,n2,s
Vn1n2 = s
end function

! Calculates <n1|x**p|n2> 
real(8) function Operator_n1n2(n1,n2,p)
implicit none
integer::n1,n2
real(8)::p
real(8)::s,xi,f
integer::i

s=0.d0

do i=1,NquadPosition
xi = Xmax*dx(i)

f = Rn_xi(n1,i)*Rn_xi(n2,i)*(xi**p)

s = s +f*dwX(i)*Xmax
end do

Operator_n1n2 = s
end function




! A subroutine meant for testing the subroutine
subroutine testHermiteElements()
implicit none
integer::n1,n2,i
real(8)::s,xi,f


do n1=0,Nmax
do n2=0,Nmax
	s =0.d0
	do i=1,NquadPosition
	xi = Xmax*dx(i)
	!print *, i, xi
    !xi = dx(i)
	!f = HermitePoly(n1,xi)*HermitePoly(n2,xi)!dExp(-1.d0*xi*xi)!*xi**2
	f =  Rn(n1,xi)*Rn(n2,xi)
	
	f = HermitePolyOrtho(n1,xi)*HermitePolyOrtho(n2,xi)*dExp(-1.d0*xi*xi)

	!f =  HermitePoly(n1,xi)
	!f = dexp(-1.d0*xi*xi)
	!f = xi*dexp(-xi*xi)
	!f= (4.d0*xi*xi-2.d0)*dexp(-1.d0*xi*xi)
	!!f=1.d0
	s = s +f*dwX(i)*Xmax
	end do
	
	print *, n1,n2,s

end do
end do



end subroutine


subroutine build_H_matrix()
implicit none
integer::n1,n2

allocate(H(Nmax+1,Nmax+1))
allocate(Eig(Nmax+1))

! Now we fill the Matrix
do n1=0,Nmax
	do n2=0,Nmax
		H(n1+1,n2+1) = KineticEnergy(n1,n2)+Vn1n2(n1,n2)
	end do
end do



end subroutine



SUBROUTINE diagonalize(M,n,W)
INTEGER :: info, lwork, lwmax = 2000000
integer::n
REAL*8 :: M(n,n)
real*8::W(n)
REAL*8, DIMENSION(:), ALLOCATABLE :: work
INTEGER, DIMENSION(:), ALLOCATABLE :: IPIV
!REAL*8, DIMENSION(:), ALLOCATABLE :: W !Vector storing eigenvalues
ALLOCATE(work(lwmax))
ALLOCATE(IPIV(n))
!ALLOCATE(W(n))

! Optimizing lwork
  LWORK = -1
  CALL dsyev( 'Vectors', 'Upper', n, M, n, W, work, lwork, info )
  lwork = min( lwmax, int( work( 1 ) ) )
  CALL dsyev( 'Vectors', 'Upper', n, M, n, W, work, lwork, info )

  IF( info.GT.0 ) THEN
  
  print *, 'Algorithm Failed to Converge'
  
  else if(info.LT.0) THEN
  
  print *, 'the ith argument had an illegal Value', INFO*(-1)
  
  END IF

DEALLOCATE(IPIV,work)
RETURN
END SUBROUTINE diagonalize

! This subroutine will diagnalize our array
subroutine solveSchrodinger()
implicit none
integer:: i

allocate(U(0:Nmax,0:Nmax))
U= H(:,:) ! We make a deep copy of the Hamiltonian Matrix

! now we transpose u
allocate(Udagger(0:Nmax,0:Nmax))
Udagger = Transpose(U)


! Now we diagonalize the Hamiltonian, and store Eigenvectors in U
call diagonalize(U,Nmax+1,Eig)

!=================================================================
! Store the ground state eigenvector
allocate(psi0(0:Nmax))
psi0(:)=0.d0

do i=0,Nmax
psi0(i) = U(i,0)
end do

E0 = Eig(1)
!=================================================================

do i=1,10
Print *, i,Eig(i)
end do

Print *,''

do i=1,10
Print *, i, dsqrt(2.d0*Eig(i)/(hw))
end do

end subroutine

! This will evaluate our ith eigenfunction Psi(x) at x,
real(8) function psi(x,i)
implicit none
real(8)::x,s
integer::i,j

s=0.d0

do j=0,Nmax
s = s+U(j,i)*Rn(j,x)
end do

psi = s

end function


! This subroutine will then generate the wave functions for examination purposes
! Phi(x,Ei) = Sum_{i} C_i Phi(n,x)
! Generates the wavefunction in the region from a to b
subroutine generateWavefunction(a,b,ni)
implicit none
integer::ni ! The ith eigenvector
real(8)::a,b,xi
integer::i
!real(8):: Matrix(0:Nmax,0:Nmax)
character(50):: f_out,f_out2,x1
real(8)::dh

character(len=8) :: fmt ! format descriptor
fmt = '(I3.3)' ! an integer of width 5 with zeros at the left
write(x1,fmt) ni ! converting integer to string using a 'internal file'

f_out='EigenVec_'//trim(x1)//'.txt'
open(unit = 2, file = f_out)


!steps = integer((b-a)/dh)
dh = abs(b-a)/dfloat(steps)

do i=0,steps
xi = a+dfloat(i)*dh
write(2,*), xi,psi(xi,ni)
end do

close(2)

! Print off the eigenvalues of the Eigenvectors
f_out2 = 'EigenVectors.txt'
open(unit=3, file = f_out2)
do i=1,Nmax+1
write(3,*) i,Eig(i)
end do
close(3)

end subroutine




end module
