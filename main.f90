program main
use mod_parameters
use mod_matrix_elements

! Read the input file
call readInput()

print *, 'Nmax', Nmax
print *, 'hw', hw
print *, 'v',v

! Initialize the Legendre quadrature mesh
call initMesh()


! Construct the Hamiltonian Matrix
call build_H_matrix()

! Diagonalize the Matrix
!call diagonalize_H_matrix()
call testHermiteElements()

!call print_Eigen_Values()


end program
