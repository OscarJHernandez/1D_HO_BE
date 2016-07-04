! Here we define the HO-basis functions

module HO_basis
use mod_parameters


contains


real(8) function Rn(n,x)
integer::n
real(8)::x

Rn = Norm(n)*HermitePoly(n,x)

end function

! Normalization
!
real(8) function Norm(n)
integer::n

Norm = 1

end function


! The Hermite Polynomials
real(8) function HermitePoly(n,x)
integer::n,ni
real(8)::x
real(8)::um1, um2,up1
real(8),allocatable::pl(:)
real(8),allocatable::dpl(:)
integer,parameter::kfn = 4


allocate(pl(0:n+4),dpl(0:n+4))

!um2 = 1.d0
!um1 = 2.d0*x

!if(n==0) then
!HermitePoly = 1.d0
!else if(n==1) then
!HermitePoly = 2.d0*x

!else

!	do ni =2,n
!	up = 2.d0*x*um1-2.d0*dfloat(ni-1)*um2
!	um2 = um1
!	um1 = up 
!	end do

!HermitePoly = up
!end if


call othpl( kf, n+4, x, pl, dpl )

HermitePoly = pl(n)

deallocate(pl,dpl)

end function



end module
