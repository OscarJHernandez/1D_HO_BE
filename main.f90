program main
use mod_parameters
use mod_matrix_elements
use mod_potentials
use mod_excitation_operator
implicit none
integer::i
real(8)::s

! Read the input file
call readInput()

print *, 'Nmax', Nmax
print *, 'hw', hw
print *, 'v',v
print *, 'Xmax',Xmax
print *, 'V0', V0
print *, 'L0', L0
print *, 'Nquad', NquadPosition

! Initialize the Legendre quadrature mesh
call initMesh()


! Construct the Hamiltonian Matrix
call build_H_matrix()

!do i=1,Nmax+1
!do j=1,Nmax+1
!Print *, i,j,H(i,j)
!end do
!end do

call solveSchrodinger()

call excite_OpPsi0(1.d0)

s=0.d0
do i=0,Nmax
s = s+psi0(i)**2
end do

print *,'psi0', s


s=0.d0
do i=0,Nmax
s = s+gamma0(i)**2
end do

print *,'gamma0', s

!================================
! Lorentz vector
do i=1,4*Nmax+1
s= Eig(2)-Eig(1)
s = s+0.005d0*dfloat(i)
print *,s,build_H_matrix_complex(s,0.01d0)
end do
!print *,psi_tilde(0)
!print *,psi_tilde(10)
!=================================



! Here we print off the first few eigen_vectors
do i=0,9
call generateWavefunction(-2000.d0,2000.d0,i)
end do
! Diagonalize the Matrix
!call diagonalize_H_matrix()
!call testHermiteElements()

!call print_Eigen_Values()


end program
