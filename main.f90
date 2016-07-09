program main
use mod_parameters
use mod_matrix_elements
use mod_potentials
implicit none
integer::i

! Read the input file
call readInput()

print *, 'Nmax', Nmax
print *, 'hw', hw
print *, 'v',v
print *, 'Xmax',Xmax
print *, 'V0', V0
print *, 'L0', L0

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

! Here we print off the first eigen_vector
do i=0,9
call generateWavefunction(-2000.d0,2000.d0,i)
end do
! Diagonalize the Matrix
!call diagonalize_H_matrix()
!call testHermiteElements()

!call print_Eigen_Values()


end program
