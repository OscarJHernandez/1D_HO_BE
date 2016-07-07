program main
use mod_parameters
use mod_matrix_elements
implicit none
integer::i,j

! Read the input file
call readInput()

print *, 'Nmax', Nmax
print *, 'hw', hw
print *, 'v',v

! Initialize the Legendre quadrature mesh
call initMesh()


! Construct the Hamiltonian Matrix
call build_H_matrix()

!do i=1,Nmax+1
!do j=1,Nmax+1
!Print *, i,j,H(i,j)
!end do
!end do

call diagonalize(H,Nmax+1,Eig)

do i=1,Nmax+1
Print *, i,Eig(i)
end do

! Diagonalize the Matrix
!call diagonalize_H_matrix()
!call testHermiteElements()

!call print_Eigen_Values()


end program
