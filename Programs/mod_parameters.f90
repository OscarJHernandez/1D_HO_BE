module mod_parameters
real(8),parameter:: mn=1.d0 ! Units of MeV
real(8),parameter:: hbarc =197.327d0
real(8),parameter:: pi = 4.d0*datan(1.d0) ! The definition of Pi
real(8):: hw
real(8):: v
real(8):: Xmax
integer:: Nmax ! Max number of basis functions
integer:: NquadPosition ! Number of quadrature points in x
integer:: NquadMomentum ! Number of quadrature points in p
integer::st

contains

! This subroutine reads in the input from a file
! Must read Nmax,hw
! also compute v
subroutine readInput()
implicit none
character(50)::word,title

title = 'input.txt'
open(unit=1,file=title, iostat = st)
read(1,*) word,hw
read(1,*) word,Nmax
read(1,*) word,Xmax
read(1,*) word,NquadPosition
read(1,*) word,NquadMomentum
close(1)


! Initialize v
v = (0.5d0*hw*mn)/(hbarc**2)


end subroutine



end module
