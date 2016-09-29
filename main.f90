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

call excite_OpPsi0(0.d0)

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
open(unit = 1, file = "response_1.txt")

do i=1,200*Nmax+1
!s= Eig(2)-Eig(1)
s = 0.0005d0*dfloat(i)
write(1,*) s, build_H_matrix_complex(s,0.01d0)

end do

close(1)
!=====================================================================
print *,'overlap,1',calculate_Overlap(0)
print *,'overlap,2',calculate_Overlap(1)
print *,'overlap,3',calculate_Overlap(2)

!====================================================================
! Integrate the response around First contiuum eigenstate
print *, 'w0',(Eig(1)-Eig(1))
print *, integrate_response(1,0.1d0,0.1d0)
print *, integrate_response(1,0.01d0,0.1d0)
print *, integrate_response(1,0.001d0,0.1d0)

print *, ''
print *, integrate_response(1,0.001d0,1.d0)
print *, integrate_response(1,0.001d0,0.5d0)
print *, integrate_response(1,0.001d0,0.25d0)
print *, integrate_response(1,0.001d0,0.125d0)

print *, ''
print *, integrate_response(1,0.0005d0,0.5d0)
print *, integrate_response(1,0.0005d0,0.25d0)
print *, integrate_response(1,0.0005d0,0.125d0)
print *, integrate_response(1,0.0005d0,0.1d0)
!=====================================================================
open(unit = 1, file = "response_2.txt")
do i=1,200*Nmax+1
!s= Eig(2)-Eig(1)
s = 0.0005d0*dfloat(i)
write(1,*) s, build_H_matrix_complex(s,0.001d0)
end do

close(1)
!=====================================================================

! Here we print off the first few eigen_vectors
do i=0,9
call generateWavefunction(-2000.d0,2000.d0,i)
end do


!print *,psi_tilde(0)
!print *,psi_tilde(10)
!=================================


! Diagonalize the Matrix
!call diagonalize_H_matrix()
!call testHermiteElements()

!call print_Eigen_Values()


end program
