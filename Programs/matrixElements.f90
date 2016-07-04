! This modules will contain all of the pieces needed to generate the matrix elements of an Operator
! < 1 || O || 2 >
module matrixElements
use deuteronMEC_params
implicit none

real(8),allocatable::daP(:),dbP(:),dp(:),dwP(:),eP(:)
real(8),allocatable::daR(:),dbR(:),dr(:),dwR(:),eR(:)
integer,parameter:: ipolyMomentum=2 !  Legendre polynomial from 0 to 1 
integer,parameter:: ipolyPosition=2 !  Legendre polynomial from 0 to 1 
real*8,parameter:: depsma=1.0d-18
real*8:: al,ierr,iderr,dbe

! Here we define some constants that will be used for our Matrix elements
real(8)::Mdelta
real(8)::Fpi
real(8)::mpi
real(8)::Mnucleon
real(8)::gA
real(8)::amp0
real(8)::ampc
real(8)::ff

real(8):: fconst ! The f-coupling constant used by Friar

! Here we have the Momentum-cutoff of the spectral fun
real(8):: Lambda

! Here we have arrays that store calculated points
real(8),allocatable:: LagExpArray(:,:,:) ! The array that stores the associated Laguerre Exp(-2*v*x*x)L(n,l+1/2,x) polynomials
real(8),allocatable:: SphericalBessel(:,:) ! Array that stores the spherical Bessel functions evaluated at quadrature points
real(8),allocatable:: PhiFunction(:,:,:,:) ! Array that stores the Phi functions evaluated at the quadrature points
real(8),allocatable:: RadialSeagullME(:,:,:,:,:) ! Array that stores the radial ME for the seagull term
real(8),parameter::  dh =0.00001d0 ! The small parameter used to define derivatives


contains 


! Here we determine which set of constants to use
subroutine init_ME_params()
implicit none
     
    print *,'constants',constants 
    
	! Friar's Constants
	if(constants.eq.0) then
	Mdelta = 1232.d0 ! [MeV] The mass of the delta isobar resonance
	mpi = hbarc*0.70728055d0 ! To be used for comparison with Friar
	Mnucleon = 938.9265d0! The average of the proton and neutron mass for use with Friar
	gA = 1.27d0  ! The nucleon axial decay constant (used for Friar and Sauri Constants)
	ff = 0.079d0
    Fpi = (mpi*gA)/(dsqrt(4.d0*pi*ff))
	
    ! Saori's Constants
	else if(constants.eq.1) then
	Mdelta = 1232.d0 ! [MeV] The mass of the delta isobar resonance
	Fpi = 184.6d0 ! [MeV] The Pion decay amplitude taken from Phys Rev C 78,064002 (2008)
	amp0 =134.9766d0  ! MeV            ,
	ampc =139.5702d0   ! MeV
	mpi = (amp0+2.d0*ampc)/3.d0 ! MeV,
	Mnucleon =mp 
	gA = 1.27d0  ! The nucleon axial decay constant (used for Friar and Sauri Constants)
    
    ! Here we employ Arenhovels constants
    else if(constants.eq.2) then
    mpi = 0.6995d0*hbarc ! MeV
    Mnucleon = 938.9055d0 ! MeV
    gA = 1.27d0 ! This is the same for Friar and Sauri constants
    ff = 0.08d0
    Fpi = (mpi*gA)/(dsqrt(4.d0*pi*ff))
    end if
    

    ! Here we initialize the cut-off
    Lambda = cutOff
    
print *, 'The cut-off [MeV]: ', Lambda
print *,'mpi',mpi
print *, 'Initialized'

end subroutine

! Initialize all of the guassian quadrature points and weights for momentum and position
subroutine gauleg_init()
implicit none

call init_ME_params()


    allocate(daP(NquadMomentum),dbP(NquadMomentum),dp(NquadMomentum),dwP(NquadMomentum),eP(NquadMomentum))
    allocate(daR(NquadPosition),dbR(NquadPosition),dr(NquadPosition),dwR(NquadPosition),eR(NquadPosition))
    
    al = 0.d0
    
    call drecur(NquadMomentum,ipolyMomentum,al,dbe,daP,dbP,iderr)
    call dgauss(NquadMomentum,daP,dbP,depsma,dp,dwP,ierr,eP)
    
    call drecur(NquadPosition,ipolyPosition,al,dbe,daR,dbR,iderr)
    call dgauss(NquadPosition,daR,dbR,depsma,dr,dwR,ierr,eR)
    
    ! Factorials are also initialized, this will be used by angular reduced ME
    call factorial

end subroutine


! This subroutine initializes all of the Norm(n,l)*Laguerre*Exp(-0.5*v*r*r)*r**(L+1) functions in the basis
subroutine initExpLaguerre()
implicit none
integer::ni,li,i
real(8)::xi,alphai,ri

! Determine the Max order of n an L
nPolyMax = Nmax/2 ! max order of N
lPolyMax = (jjMax+2)/2 ! Max order of L
print *,'nPolyMax,lpolyMax', nPolyMax,lPolyMax


allocate(LagExpArray(0:nPolyMax,0:lPolyMax,NquadPosition))

do ni=0,nPolyMax
    do li=0,lPolyMax
        do i=1,NquadPosition
            ri = rMax*dr(i)
            xi = 2.d0*v*ri*ri
            alphai = 0.5d0+dfloat(li)
            LagExpArray(ni,li,i) = Norm(ni,li,v)*LagExp(ni,alphai,xi)*(ri**(li+1))
        end do
    end do
end do


end subroutine


! This function tests the orthognonality of the Basis functions, with respect to n1 and n2
real(8) function radialMeTest(n1,l1,n2,l2)
implicit none
integer::n1,n2,l1,l2
real(8)::s
integer::i

s =0.d0

do i=1,NquadPosition
s = s+LagExpArray(n1,l1,i)*LagExpArray(n2,l2,i)*dwR(i)
end do

s = rMax*s

radialMeTest =s

end function



! Calculates the radial matrix elements of the operator r^p 
! Note, we must divide the operator by r/2 
real(8) function radial_matrixElement_rp(n1,l1,n2,l2,k)
implicit none
integer::n1,n2,l1,l2,k
real(8)::ri
real(8)::s
integer::i

s=0.d0

do i =1,NquadPosition
ri = rMax*dr(i)
s = s+LagExpArray(n1,l1,i)*LagExpArray(n2,l2,i)*dwR(i)*((ri/2.d0)**k)
end do

s=rMax*s

radial_matrixElement_rp =s

end function





! This function evaluates the associated Laguerre polynomial with an exponential at a given point.
! This function behaves quite well even for large values of x. The exponential helps keep the polynomial from
! overflow problems.
! LagExp = Exp(-1/2*x)*Lag[n,alpha,x]
real*8 function LagExp(n,alpha,x)
implicit none
integer,intent(in):: n
real*8,intent(in):: alpha
real*8,intent(in):: x
real*8::L,Lm,Lp
integer:: i 

! Initialize the values
Lm = 1.d0*dexp(-0.5d0*x)
L = (1.d0+ alpha -x)*dexp(-0.5d0*x)
Lp = 0.d0


    if(n.eq.0) then
        Lp = Lm
    else if (n.eq.1) then
        Lp = L
    else
    
        do i=1,n-1
        Lp = ((2.d0*dfloat(i)+1.d0+alpha-x)/(dfloat(i)+1.d0)*L - ((dfloat(i)+alpha)/(dfloat(i)+1.d0))*Lm)
        Lm = L
        L = Lp
        end do
    
    end if
    
    LagExp = Lp

end function

! This is the Normalization coefficients for the harmonic Oscillator
! This function computes the normalization for the Harmonic Oscillator wave function
! N(n,l,v) = sqrt( Gamma(n+1)/Gamma(n+l+3/2) * 2 *(2*v)^(l+3/2))
real*8 function Norm(n,l,v)
implicit none
integer:: n,l,i
integer:: Nm
real*8::v
real*8:: frac
real*8:: pi


pi = datan(1.d0)*4.d0


frac = 1.d0/dsqrt(pi)

do i=1,n
    frac =frac*(dfloat(i)/(dfloat(i)-0.5d0))
end do

Nm = n+l+1

do i=n+1,Nm
    frac = frac*(1.d0/(dfloat(i) -0.5d0))
end do

frac = frac*2.d0*((2.d0*v)**(dfloat(l)+1.5d0))
Norm = dsqrt(frac)

return
end function Norm



! This is the reduced Matrix Element: <(l2 s)j2||(Y^p)_q||(l1 s)j1>
real*8 function reduced_M_element_1(ll1,ss1,jj1,ll2,ss2,jj2,pp)
integer::ll1,jj1,jj2,ll2,ss1,ss2,pp
real*8:: l1,l2,j1,j2,p,s
real*8,external::CG ! Generates Clebsh-gordan Coefficients
real*8,external::SIXJ ! Generates the six-j symbols
!real*8,external::reduced_M_element_2
logical,external::TRIANG
real*8::factor1
integer::mm1,mm2,mm3

mm1=0
mm2=0
mm3=0

!call factorial


    if(ss1.eq.ss2) then
    
       !if(TRIANG(jj1,jj2,pp).and.TRIANG(ll1,ss1,jj1).and.TRIANG(ll2,ss2,jj2)) then
        l1=DBLE(ll1)*0.5d0
        l2=DBLE(ll2)*0.5d0
        j1=DBLE(jj1)*0.5d0
        j2=DBLE(jj2)*0.5d0
        s=DBLE(ss1)*0.5d0
        p=DBLE(pp)*0.5d0

        factor1=  (1-2*mod(ll2/2+ss1/2+jj1/2+pp/2,2))*dsqrt((2.d0*j1+1.d0)*(2.d0*j2+1.d0))

        reduced_M_element_1 = factor1*SIXJ(ll2,jj2,ss1,jj1,ll1,pp)*reduced_M_element_Y(ll2,pp,ll1)
    else
        reduced_M_element_1 =0.d0
    end if

return
end function reduced_M_element_1 


! This is the reduced Matrix Element <l2||y^p ||l1>
real*8 function reduced_M_element_Y(ll1,pp,ll2)
integer::ll2,pp,ll1
real*8,external::CG ! Generates Clebsh-gordan Coefficients
real(8),external::THREEJ
real*8:: l1,l2,p,mm1,mm2,mm3
!call factorial

l1=DBLE(ll1)*0.5d0
l2=DBLE(ll2)*0.5d0
p=DBLE(pp)*0.5d0

mm1=0
mm2=0
mm3=0

!reduced_M_element_Y = (1-2*mod(pp/2,2))*dsqrt(((2.d0*l2+1.d0)*(2.d0*p+1.d0))/(4.d0*Pi))*CG(ll2,mm1,pp,mm2,ll1,mm3)
reduced_M_element_Y = ((-1)**(l2))*dsqrt(dfloat(ll1+1))*dsqrt(dfloat(ll2+1))*dsqrt(dfloat(pp+1))*(1.d0/dsqrt(4.d0*pi))
reduced_M_element_Y = reduced_M_element_Y*THREEJ(ll2,pp,ll1,0,0,0)

end function reduced_M_element_Y

! This function computes <s1 l1 j1|| S || s0 l0 j0 >, where S is the total spin operators for the proton 
! and neutron respectively. 
real*8 function reduced_M_element_S(ll1,ss1,jj1,ll0,ss0,jj0)
implicit none    
integer::ll1,ss1,jj1
integer:: ll0,ss0,jj0
real*8:: factor1,factor2,s1
real*8,external::CG ! Generates Clebsh-gordan Coefficients
real*8,external::SIXJ ! Generates the six-j symbols
real(8)::factor0

!First we call Factorial, to generate factorial numbers:
!call factorial

if(ll1.eq.ll0) then
factor0 = 1.d0
else
factor0 = 0.d0
end if

factor1 = (1-2*mod( ll0/2+ss0/2+jj1/2+1,2))*dsqrt((DBLE(jj0)+1.d0)*(DBLE(jj1)+1.d0))*SIXJ(ss1,jj1,ll0,jj0,ss0,2)

    if(ss1.eq.ss0) then

    s1 = DBLE(ss1)*0.5d0
    factor2 = dsqrt(s1*(s1+1.d0)*(2.d0*s1+1))

    else
    factor2=0.d0
    end if


reduced_M_element_S = factor0*factor1*factor2

end function reduced_M_element_S

! This function computes <s1 l1 j1|| L || s0 l0 j0 >, where L is the relative angular momentum operator for the proton 
! and neutron. 
real*8 function reduced_M_element_L(ll1,ss1,jj1,ll0,ss0,jj0)
implicit none    
integer::ll1,ss1,jj1
integer:: ll0,ss0,jj0
real*8:: factor1,factor2,s1,l1
real*8,external::CG ! Generates Clebsh-gordan Coefficients
real*8,external::SIXJ ! Generates the six-j symbols

!First we call Factorial, to generate factorial numbers:
!call factorial
    
    if(ss1.eq.ss0) then
    factor1 = (1-2*mod( ll1/2+ss0/2+jj0/2+1,2))*dsqrt((DBLE(jj0)+1.d0)*(DBLE(jj1)+1.d0))*SIXJ(ll1,jj1,ss0,jj0,ll0,2)
    else
    factor1=0.d0
    end if

    if(ll1.eq.ll0) then
    l1 =DBLE(ll1)/2.d0
    factor2 = dsqrt(l1*(l1+1.d0)*(2*l1+1))
    else
    factor2=0.d0
    end if


reduced_M_element_L = factor1*factor2

end function reduced_M_element_L

! This function computes <s1 l1 j1|| s^p || s0 l0 j0 >, where s^p, s^n are the spin operators for the proton 
! and neutron respectively. This is the M1 magnetic transition operator
real*8 function reduced_M_element_sp(ll1,ss1,jj1,ll0,ss0,jj0)    
implicit none
integer::ll1,ss1,jj1
integer:: ll0,ss0,jj0
real*8:: factor1,factor2,factor3
real*8,external::CG ! Generates Clebsh-gordan Coefficients
real*8,external::SIXJ ! Generates the six-j symbols

!First we call Factorial, to generate factorial numbers:
!call factorial

!factor1 = ((-1)**( ll0/2 + jj1/2 +1))*sqrt((DBLE(jj0)+1.d0)*(DBLE(jj1)+1.d0))*SIXJ(ss1,jj1,ll0,jj0,ss0,2)
factor1 = (1-2*mod( ll0/2 + jj1/2 +1,2))*dsqrt((DBLE(jj0)+1.d0)*(DBLE(jj1)+1.d0))*SIXJ(ss1,jj1,ll0,jj0,ss0,2)
factor2 = dsqrt( (DBLE(ss0)+1.d0)*(DBLE(ss1)+1.d0) )*SIXJ(1,ss1,1,ss0,1,2)
factor3 = dsqrt(3.d0/2.d0)


reduced_M_element_sp = factor1*factor2*factor3

end function reduced_M_element_sp


! This function computes <s1 l1 j1|| s^n || s0 l0 j0 >, where s^p, s^n are the spin operators for the proton 
! and neutron respectively. This is the M1 magnetic transition operator
real*8 function reduced_M_element_sn(ll1,ss1,jj1,ll0,ss0,jj0)    
implicit none
integer::ll1,ss1,jj1
integer:: ll0,ss0,jj0
real*8:: factor1,factor2,factor3
real*8,external::CG ! Generates Clebsh-gordan Coefficients
real*8,external::SIXJ ! Generates the six-j symbols

!First we call Factorial, to generate factorial numbers:
!call factorial

!factor1 = ((-1)**( ll0/2 + jj1/2 + ss0/2 + ss1/2 +1))*dsqrt((DBLE(jj0)+1.d0)*(DBLE(jj1)+1.d0))*SIXJ(ss1,jj1,ll0,jj0,ss0,2)
factor1 = (1-2*MOD( ll0/2 + jj1/2 + ss0/2 + ss1/2 +1,2))*dsqrt((DBLE(jj0)+1.d0)*(DBLE(jj1)+1.d0))*SIXJ(ss1,jj1,ll0,jj0,ss0,2)
factor2 = dsqrt( (DBLE(ss0)+1.d0)*(DBLE(ss1)+1.d0) )*SIXJ(1,ss1,1,ss0,1,2)
factor3 = dsqrt(3.d0/2.d0)


reduced_M_element_sn = factor1*factor2*factor3

end function reduced_M_element_sn

!! < ||sp -sn||>
real*8 function reduced_M_element_spMsn(ll1,ss1,jj1,ll0,ss0,jj0)    
implicit none
integer::ll1,ss1,jj1
integer:: ll0,ss0,jj0
real*8:: factor1,factor2,factor3
real*8,external::CG ! Generates Clebsh-gordan Coefficients
real*8,external::SIXJ ! Generates the six-j symbols

!First we call Factorial, to generate factorial numbers:
!call factorial

factor1 = (1-2*mod( ll0/2+ss0/2+jj1/2+1,2))*dsqrt((DBLE(jj0)+1.d0)*(DBLE(jj1)+1.d0))*SIXJ(ss1,jj1,ll0,jj0,ss0,2)
factor2 = dsqrt( (DBLE(ss0)+1.d0)*(DBLE(ss1)+1.d0) )*SIXJ(1,ss1,1,ss0,1,2)*( (-1)**(ss0/2) -(-1)**(ss1/2) )
factor3 = dsqrt(3.d0/2.d0)

reduced_M_element_spMsn = factor1*factor2*factor3

end function reduced_M_element_spMsn   

! M1 operator
! M1 = ((gp+gn)/2 S+ 1/2*L+(gp-gn)/2*(sp-sn))*(1/(2.d0*mp))
real(8) function M1_operator(n1,ll1,ss1,jj1,n0,ll0,ss0,jj0,t,mt,t0,mt0)
implicit none
integer::ll1,ss1,jj1
integer:: ll0,ss0,jj0
integer::n1,n0
integer:: t0,mt,t,mt0
integer::tt0,tt
real(8)::s,factor1

if(n0.eq.n1) then
factor1 = 1.d0
else
factor1 = 0.d0
end if


s = 0.d0
! The Isoscalar operators
if((t.eq.0).and.(mt.eq.0)) then
s = 0.5d0*(gp+gn)*reduced_M_element_S(ll1,ss1,jj1,ll0,ss0,jj0)
s = s+0.5d0*reduced_M_element_L(ll1,ss1,jj1,ll0,ss0,jj0)

! The isovector operator
else if((t.eq.1).and.(mt.eq.0)) then
s = 0.5d0*(gp-gn)*reduced_M_element_spMsn(ll1,ss1,jj1,ll0,ss0,jj0)
end if

M1_operator = s*(1.d0/(2.d0*Mnucleon))*factor1

end function




! Reduced matrix elements of Sigma^{k}_{12} for the deuteron
! <S0(1/2,1/2)||Sigma^{k}_{12} || S(1/2,1/2)>
! <2|| Sigma ||1>
real(8) function reduced_ME_Sigma12_S(ss1,k,ss2)
implicit none
integer::k,kk
integer::ss1,ss2
real(8)::factor1,factor2,factor3
real(8),external::ANINEJ

kk=2*k

factor1 = 6.d0
factor2 = dsqrt(dfloat(ss1+1))*dsqrt(dfloat(ss2+1))*dsqrt(dfloat(kk+1))
factor3 = ANINEJ(1,1,2,1,1,2,ss2,ss1,kk)

reduced_ME_Sigma12_S = factor1*factor2*factor3

end function

! This is just for the M1 MEC operator
! The reduces matrix element <J0(S0,L0)||[Sigma^{k1}_{12} x Y^k2 ]^k|| J(S,L)>
! <2||Omega||1>
real(8) function reduced_ME_Omega(l1,jj1,ss1,l2,jj2,ss2,k1,k2,k)
implicit none
integer::l1,l2,jj1,jj2,ss1,ss2
integer::ll1,ll2
integer::k1,k2,k
integer::kk1,kk2,kk
real(8)::factor1,factor2,factor3,factor4
real(8),external::ANINEJ

kk1=2*k1
kk2=2*k2
kk=2*k
ll1=2*l1
ll2=2*l2
factor1 = dsqrt(dfloat(kk+1))*dsqrt(dfloat(jj1+1))*dsqrt(dfloat(jj2+1))
factor2 = ANINEJ(ss2,ss1,kk1,ll2,ll1,kk2,jj2,jj1,kk)
factor3 = reduced_ME_Sigma12_S(ss1,k1,ss2)
factor4 = reduced_M_element_Y(ll1,kk2,ll2)


reduced_ME_Omega = factor1*factor2*factor3*factor4

end function

! Reduced matrix elements of Sigma^{k}_{12} for the deuteron
! <J0(S0,L0)||Sigma^{k}_{12} || J(S,L)>
! <1 || Operator ||2>
real(8) function reduced_ME_Sigma12_J(l1,jj1,ss1,k,l2,jj2,ss2)
implicit none
integer::l1,l2
integer::ll1,ll2
integer::jj1,jj2,ss1,ss2
integer::j1,j2
integer::k,kk
integer::s1,s2
real(8)::factor1,factor2,factor3,factor4
real(8),external::ANINEJ
real(8),external::SIXJ

kk=2*k
s1 = ss1/2
j1=jj1/2
j2=jj2/2
ll1=2*l1
ll2=2*l2

factor1 = (-1)**(l1+s1+j2+k)
if(l1.eq.l2) then
factor2 = dsqrt(dfloat(jj1+1))*dsqrt(dfloat(jj2+1))
else
factor2 = 0.d0
end if

factor3 = SIXJ(ss2,jj2,ll1,jj1,ss1,kk)
factor4 = reduced_ME_Sigma12_S(ss1,k,ss2)

reduced_ME_Sigma12_J= factor1*factor2*factor3*factor4

end function

! This is the isospin matrix element for the deuteron
! <T_0 ||[tau1 x tau2 ]_0|| T>
! <2|| Operator ||1>
real(8) function matrixElement_MEC_Isospin(t1,t2) 
implicit none
integer::t1,t2,tt1,tt2
real(8)::factor1
real(8),external::ANINEJ

tt1=2*t1
tt2=2*t2

factor1 = 6.d0*dsqrt(3.d0)*dsqrt(dfloat(2*t1+1))*dsqrt(dfloat(2*t2+1))
factor1 = factor1*ANINEJ(1,1,2,1,1,2,tt2,tt1,2)

matrixElement_MEC_Isospin = factor1

end function matrixElement_MEC_Isospin

! This is the full isospin matrix element
! <T2 M2 ||[tau^1_1 x tau^1_2]_0 || T1, M1>
real(8) function full_matrixElement_MEC_Isospin(t1,m1,t2,m2)
implicit none
integer::t1,m1,t2,m2
integer::tt1,tt2,mm1,mm2
real(8)::s
real(8),external::CG

! Using Wigner Eckart theorem
tt1 = 2*t1
tt2 = 2*t2
mm1=2*m1
mm2=2*m2
s = CG(2,0,tt1,mm1,tt2,mm2)*((-1)**(1-t1+t2))*matrixElement_MEC_Isospin(t1,t2)*(1.d0/(dsqrt(dfloat(2*t2+1))))
!s = CG(2,0,tt1,0,tt2,0)*((-1)**(1-t1+t2))*matrixElement_MEC_Isospin(t1,t2)*(1.d0/(dsqrt(dfloat(2*t2+1))))
full_matrixElement_MEC_Isospin =s

! The analytical result should be: (and delta_m1,0 and delta_m2,0)
!full_matrixElement_MEC_Isospin=dsqrt(2.d0)*((-1)**(t2+1))

end function full_matrixElement_MEC_Isospin




end module

