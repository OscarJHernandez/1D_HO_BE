c      program math
c      implicit real*8 (a-h,o-z)
c      parameter(n=100)
c      !real*8 :: wg(n),zre(n)
c      real*8 :: da(n),db(n),dx(n),dw(n),e(n)
c      real*8 :: s
c      integer:: i 
      
c      ipoly=2
c      depsma=1.0d-18       
c      al=0.0d0
c      call drecur(n,ipoly,al,dbe,da,db,iderr)
c      call dgauss(n,da,db,depsma,dx,dw,ierr,e)
c      !write(*,*) "zeros"
c      !write(*,*) dx
c      !write(*,*) "weights"
c      !write(*,*) dw
      
c      s = 0.d0
c      do i=1,n
c      s = s+ dw(i)*(1.d0/(1+dx(i)*dx(i)))
c      end do
      
c      print *, s
      
      
      
c      end program math


      double precision function d1mach(i)
c extracted from  a package of routines, called ORTHPOL, for generating 
c orthogonal polynomials and Gauss-type quadrature rules developed
c by Walter Gautschi. A description of the underlying methods can be
c found in a companion paper published in ``ACM Transactions on
c Mathematical Software''.
c ALGORITHM 726, COLLECTED ALGORITHMS FROM ACM.
c THIS WORK PUBLISHED IN TRANSACTIONS ON MATHEMATICAL SOFTWARE,
c VOL. 20, NO. 1, MARCH, 1994, PP. 21-62.
c adapted by Jeremy Dohet-Eraly (29/08/2014)


c
c  Double-precision machine constants
c
c  d1mach( 1) = b**(emin-1), the smallest positive magnitude.
c
c  d1mach( 2) = b**emax*(1 - b**(-t)), the largest magnitude.
c
c  d1mach( 3) = b**(-t), the smallest relative spacing.
c
c  d1mach( 4) = b**(1-t), the largest relative spacing.
c
c  d1mach( 5) = log10(b)
c
      integer sc
      double precision dmach(5),x
      data sc/987/

      dmach(1)=tiny(x)
      dmach(2)=huge(x)
      dmach(3)=epsilon(x)/radix(x)
      dmach(4)=epsilon (x)
      dmach(5)=log(1.0d0*radix(x))
      if (sc .ne. 987) stop 779
c  ***  issue stop 778 if all data statements are obviously wrong...
      if (dmach(4) .ge. 1.0d0) stop 778
      if (i .lt. 1  .or.  i .gt. 5) goto 999
      d1mach = dmach(i)
      return
  999 write(*,1999) i
 1999 format(' d1mach - i out of bounds',i10)
      stop
      end



c ##################################
c #                                #
c # 2. CLASSICAL WEIGHT FUNCTIONS  #
c #                                #
c ##################################

c extracted from  a package of routines, called ORTHPOL, for generating 
c orthogonal polynomials and Gauss-type quadrature rules developed
c by Walter Gautschi. A description of the underlying methods can be
c found in a companion paper published in ``ACM Transactions on
c Mathematical Software''.
c ALGORITHM 726, COLLECTED ALGORITHMS FROM ACM.
c THIS WORK PUBLISHED IN TRANSACTIONS ON MATHEMATICAL SOFTWARE,
c VOL. 20, NO. 1, MARCH, 1994, PP. 21-62.
c adapted by Jeremy Dohet-Eraly (29/08/2014)

c This subroutine generates the coefficients  a(k),b(k), k=0,1,...,n-1,
c in the recurrence relation
c
c       p(k+1)(x)=(x-a(k))*p(k)(x)-b(k)*p(k-1)(x),
c                            k=0,1,...,n-1,
c
c       p(-1)(x)=0,  p(0)(x)=1,
c
c for some classical (monic) orthogonal polynomials, and sets  b(0)
c equal to the total mass of the weight distribution. The results are
c stored in the arrays  a,b,  which hold, respectively, the coefficients
c a(k-1),b(k-1), k=1,2,...,n.
c
c       Input:  n - - the number of recursion coefficients desired
c               ipoly-integer identifying the polynomial as follows:
c                     1=Legendre polynomial on (-1,1)
c                     2=Legendre polynomial on (0,1)
c                     3=Chebyshev polynomial of the first kind
c                     4=Chebyshev polynomial of the second kind
c                     5=Jacobi polynomial with parameters  al=-.5,be=.5
c                     6=Jacobi polynomial with parameters  al,be
c                     7=generalized Laguerre polynomial with
c                       parameter  al
c                     8=Hermite polynomial
c               al,be-input parameters for Jacobi and generalized
c                     Laguerre polynomials
c
c       Output: a,b - arrays containing, respectively, the recursion
c                     coefficients  a(k-1),b(k-1), k=1,2,...,n.
c               ierr -an error flag, equal to  0  on normal return, 
c                     equal to  1  if  al  or  be  are out of range 
c                     when  ipoly=6  or  ipoly=7, equal to  2  if  b(0) 
c                     overflows when  ipoly=6  or  ipoly=7, equal to  3 
c                     if  n  is out of range, and equal to  4  if  ipoly
c                     is not an admissible integer. In the case  ierr=2,
c                     the coefficient  b(0)  is set equal to the largest
c                     machine-representable number.
c
c The subroutine calls for the function subroutine  d1mach.
c
      subroutine drecur(n,ipoly,dal,dbe,da,db,iderr)
c
c This is a double-precision version of the routine  recur.
c
      external dgamma
      double precision dal,dbe,da,db,dlmach,d1mach,dkm1,dalpbe,dt,
     *dlga,dal2,dbe2,dgamma
      dimension da(n),db(n)
      if(n.lt.1) then
        iderr=3
        return
      end if
      dlmach=dlog(d1mach(2))
      iderr=0
      do 10 k=1,n
        da(k)=0.d0
   10 continue
      if(ipoly.eq.1) then
        db(1)=2.d0
        if (n.eq.1) return
        do 20 k=2,n
          dkm1=dble(k-1)
          db(k)=1.d0/(4.d0-1.d0/(dkm1*dkm1))
   20   continue
        return
      else if(ipoly.eq.2) then
        da(1)=.5d0
        db(1)=1.d0
        if(n.eq.1) return
        do 30 k=2,n
          da(k)=.5d0
          dkm1=dble(k-1)
          db(k)=.25d0/(4.d0-1.d0/(dkm1*dkm1))
   30   continue
        return
      else if(ipoly.eq.3) then
        db(1)=4.d0*datan(1.d0)
        if(n.eq.1) return
        db(2)=.5d0
        if(n.eq.2) return
        do 40 k=3,n
          db(k)=.25d0
   40   continue
        return
      else if(ipoly.eq.4) then
        db(1)=2.d0*datan(1.d0)
        if(n.eq.1) return
        do 50 k=2,n
          db(k)=.25d0
   50   continue
        return
      else if(ipoly.eq.5) then
        db(1)=4.d0*datan(1.d0)
        da(1)=.5d0
        if(n.eq.1) return
        do 60 k=2,n
          db(k)=.25d0
   60   continue
        return
      else if(ipoly.eq.6) then
        if(dal.le.-1.d0 .or. dbe.le.-1.d0) then
          iderr=1
          return
        else
          dalpbe=dal+dbe
          da(1)=(dbe-dal)/(dalpbe+2.d0)
          dt=(dalpbe+1.d0)*dlog(2.d0)+dlga(dal+1.d0)+dlga(dbe+1.d0)-
     *      dlga(dalpbe+2.d0)
          if(dt.gt.dlmach) then
            iderr=2
            db(1)=d1mach(2)
          else
            db(1)=dexp(dt)
          end if
          if(n.eq.1) return
          dal2=dal*dal
          dbe2=dbe*dbe
          da(2)=(dbe2-dal2)/((dalpbe+2.d0)*(dalpbe+4.d0))
          db(2)=4.d0*(dal+1.d0)*(dbe+1.d0)/((dalpbe+3.d0)*(dalpbe+
     *      2.d0)**2)
          if(n.eq.2) return
          do 70 k=3,n
            dkm1=dble(k-1)
            da(k)=.25d0*(dbe2-dal2)/(dkm1*dkm1*(1.d0+.5d0*dalpbe/dkm1)
     *        *(1.d0+.5d0*(dalpbe+2.d0)/dkm1))
            db(k)=.25d0*(1.d0+dal/dkm1)*(1.d0+dbe/dkm1)*(1.d0+dalpbe/
     *        dkm1)/((1.d0+.5d0*(dalpbe+1.d0)/dkm1)*(1.d0+.5d0*(dalpbe
     *      -1.d0)/dkm1)*(1.d0+.5d0*dalpbe/dkm1)**2)
   70     continue
          return
        end if
      else if(ipoly.eq.7) then
        if(dal.le.-1.d0) then
          iderr=1
          return
        else
          da(1)=dal+1.d0
          db(1)=dgamma(dal+1.d0,iderr)
          if(iderr.eq.2) db(1)=d1mach(2)
          if(n.eq.1) return
          do 80 k=2,n
            dkm1=dble(k-1)
            da(k)=2.d0*dkm1+dal+1.d0
            db(k)=dkm1*(dkm1+dal)
   80     continue
          return
        end if
      else if(ipoly.eq.8) then
        db(1)=dsqrt(4.d0*datan(1.d0))
        if(n.eq.1) return
        do 90 k=2,n
          db(k)=.5d0*dble(k-1)
   90   continue
        return
      else
        iderr=4
      end if
      end

      double precision function dlga(dx)
      double precision dbnum,dbden,dx,d1mach,dc,dp,dy,dt,ds
      dimension dbnum(8),dbden(8)
c 
c This routine evaluates the logarithm of the gamma function by a
c combination of recurrence and asymptotic approximation.
c
c The entries in the next data statement are the numerators and
c denominators, respectively, of the quantities B[16]/(16*15),
c B[14]/(14*13),..., B[2]/(2*1), where B[2n] are the Bernoulli
c numbers.
c
      data dbnum/-3.617d3,1.d0,-6.91d2,1.d0,-1.d0,1.d0,-1.d0,1.d0/,
     *     dbden/1.224d5,1.56d2,3.6036d5,1.188d3,1.68d3,1.26d3,3.6d2,
     *1.2d1/
c
c The quantity  dprec  in the next statement is the number of decimal
c digits carried in double-precision floating-point arithmetic.
c
      dprec=-alog10(sngl(d1mach(3)))
      dc=.5d0*dlog(8.d0*datan(1.d0))
      dp=1.d0
      dy=dx
      y=sngl(dy)
c
c The quantity  y0  below is the threshold value beyond which asymptotic
c evaluation gives sufficient accuracy; see Eq. 6.1.42 in M. Abramowitz
c and I.A. Stegun,``Handbook of Mathematical Functions''. The constants 
c are .12118868... = ln(10)/19 and .05390522... = ln(|B[20]|/190)/19.
c
      y0=exp(.121189*dprec+.053905)
   10 if(y.gt.y0) goto 20
      dp=dy*dp
      dy=dy+1.d0
      y=sngl(dy)
      goto 10
   20 dt=1.d0/(dy*dy)
c
c The right-hand side of the next assignment statement is B[18]/(18*17).
c
      ds=4.3867d4/2.44188d5
      do 30 i=1,8
        ds=dt*ds+dbnum(i)/dbden(i)
   30 continue
      dlga=(dy-.5d0)*dlog(dy)-dy+dc+ds/dy-dlog(dp)
      return
      end

      double precision function dgamma(dx,iderr)
c
c This evaluates the gamma function for real positive  dx, using the
c function subroutine  dlga.
c
      double precision dx,dlmach,d1mach,dt,dlga
      dlmach=dlog(d1mach(2))
      iderr=0
      dt=dlga(dx)
      if(dt.ge.dlmach) then
        iderr=2
        dgamma=d1mach(2)
        return
      else
        dgamma=dexp(dt)
        return
      end if
      end

C-END-OF-FILE

c ##############################
c #                            #
c # 3. MOMENT-RELATED METHODS  #
c #                            #
c ##############################


c
      subroutine dcheb(n,da,db,dnu,dalpha,dbeta,ds,iderr,ds0,ds1,ds2)
c
c This is a double-precision version of the routine  cheb.
c
      double precision da,db,dnu,dalpha,dbeta,ds,ds0,ds1,ds2,dtiny,
     *d1mach,dhuge
      dimension da(*),db(*),dnu(*),dalpha(n),dbeta(n),ds(n),
     *ds0(*),ds1(*),ds2(*)
c
c The arrays  da,db  are assumed to have dimension  2*n-1, the arrays
c dnu,ds0,ds1,ds2  dimension  2*n.
c
      nd=2*n
      dtiny=10.d0*d1mach(1)
      dhuge=.1d0*d1mach(2)
      iderr=0
      if(dabs(dnu(1)).lt.dtiny) then
        iderr=1
        return
      end if
      if(n.lt.1) then
        iderr=2
        return
      end if
      dalpha(1)=da(1)+dnu(2)/dnu(1)
      dbeta(1)=dnu(1)
      if(n.eq.1) return
      ds(1)=dnu(1)
      do 10 l=1,nd
        ds0(l)=0.d0
        ds1(l)=dnu(l)
   10 continue
      do 40 k=2,n
        lk=nd-k+1
        do 20 l=k,lk
          ds2(l)=ds1(l+1)-(dalpha(k-1)-da(l))*ds1(l)-dbeta(k-1)*ds0(l)
     *      +db(l)*ds1(l-1)
        if(l.eq.k) ds(k)=ds2(k)
   20   continue
        if(dabs(ds(k)).lt.dtiny) then
          iderr=-(k-1)
          return
        else if(dabs(ds(k)).gt.dhuge) then
          iderr=k-1
          return
        end if
        dalpha(k)=da(k)+(ds2(k+1)/ds2(k))-(ds1(k)/ds1(k-1))
        dbeta(k)=ds2(k)/ds1(k-1)
        do 30 l=k,lk
          ds0(l)=ds1(l)
          ds1(l)=ds2(l)
   30   continue
   40 continue
      return
      end

C-END-OF-FILE

c ######################################################################
c #                                                                    #
c # 4. STIELTJES, ORTHOGONAL REDUCTION, AND DISCRETIZATION PROCEDURES  #
c #                                                                    #
c ######################################################################


cat <<C-END-OF-FILE > dsti.f
c
c
      subroutine dsti(n,ncap,dx,dw,dalpha,dbeta,ierr,dp0,dp1,dp2)
c
c This is a double-precision version of the routine  sti.
c
      double precision dx,dw,dalpha,dbeta,dp0,dp1,dp2,dtiny,d1mach,
     *dhuge,dsum0,dsum1,dsum2,dt
      dimension dx(ncap),dw(ncap),dalpha(n),dbeta(n),dp0(ncap),
     *dp1(ncap),dp2(ncap)
      dtiny=10.d0*d1mach(1)
      dhuge=.1d0*d1mach(2)
      ierr=0
      if(n.le.0 .or. n.gt.ncap) then
        ierr=1
        return
      end if
      nm1=n-1
      dsum0=0.d0
      dsum1=0.d0
      do 10 m=1,ncap
        dsum0=dsum0+dw(m)
        dsum1=dsum1+dw(m)*dx(m)
   10 continue
      dalpha(1)=dsum1/dsum0
      dbeta(1)=dsum0
      if(n.eq.1) return
      do 20 m=1,ncap
        dp1(m)=0.d0
        dp2(m)=1.d0
   20 continue
      do 40 k=1,nm1
        dsum1=0.d0
        dsum2=0.d0
        do 30 m=1,ncap
          if(dw(m).eq.0.d0) goto 30
          dp0(m)=dp1(m)
          dp1(m)=dp2(m)
          dp2(m)=(dx(m)-dalpha(k))*dp1(m)-dbeta(k)*dp0(m)
          if(dabs(dp2(m)).gt.dhuge .or. dabs(dsum2).gt.dhuge) then
            ierr=k
            return
          end if
          dt=dw(m)*dp2(m)*dp2(m)
          dsum1=dsum1+dt
          dsum2=dsum2+dt*dx(m)
   30   continue
        if(dabs(dsum1).lt.dtiny) then
          ierr=-k
          return
        end if
        dalpha(k+1)=dsum2/dsum1
        dbeta(k+1)=dsum1/dsum0
        dsum0=dsum1
   40 continue
      return
      end

C-END-OF-FILE

cat <<C-END-OF-FILE > lancz.f
c
c
 
      subroutine dlancz(n,ncap,dx,dw,dalpha,dbeta,ierr,dp0,dp1)
c
c This is a double-precision version of the routine  lancz.
c
      double precision dx(ncap),dw(ncap),dalpha(n),dbeta(n),
     *dp0(ncap),dp1(ncap),dpi,dgam,dsig,dt,dxlam,drho,dtmp,
     *dtsig,dtk
      if(n.le.0 .or. n.gt.ncap) then
        ierr=1
        return
      else
        ierr=0
      end if
      do 10 i=1,ncap
        dp0(i)=dx(i)
        dp1(i)=0.d0
   10 continue
      dp1(1)=dw(1)
      do 30 i=1,ncap-1
        dpi=dw(i+1)
        dgam=1.d0
        dsig=0.d0
        dt=0.d0
        dxlam=dx(i+1)
        do 20 k=1,i+1
          drho=dp1(k)+dpi
          dtmp=dgam*drho
          dtsig=dsig
          if(drho.le.0.d0) then
            dgam=1.d0
            dsig=0.d0
          else
            dgam=dp1(k)/drho
            dsig=dpi/drho
          end if
          dtk=dsig*(dp0(k)-dxlam)-dgam*dt
          dp0(k)=dp0(k)-(dtk-dt)
          dt=dtk
          if(dsig.le.0.d0) then
            dpi=dtsig*dp1(k)
          else
            dpi=(dt**2)/dsig
          end if
          dtsig=dsig
          dp1(k)=dtmp
   20   continue
   30 continue
      do 40 k=1,n
        dalpha(k)=dp0(k)
        dbeta(k)=dp1(k)
   40 continue
      return
      end 

C-END-OF-FILE

c ###############################
c #                             #
c # 5. MODIFICATION ALGORITHMS  #
c #                             #
c ###############################


      subroutine dchri(n,iopt,da,db,dx,dy,dhr,dhi,dalpha,dbeta,ierr)
c
c This is a double-precision version of the routine  chri.
c
      double precision da,db,dx,dy,dhr,dhi,dalpha,dbeta,deps,d1mach,
     *de,dq,ds,dt,deio,dd,der,dei,deroo,deioo,dso,dero,deoo,deo,du,dc,
     *dc0,dgam,dcm1,dp2
      dimension da(*),db(*),dalpha(n),dbeta(n)
c
c The arrays  da,db  are assumed to have dimension  n+1.
c
      deps=5.d0*d1mach(3)
      ierr=0
      if(n.lt.2) then
        ierr=1
        return
      end if
      if(iopt.eq.1) then
        de=0.d0
        do 10 k=1,n
          dq=da(k)-de-dx
          dbeta(k)=dq*de
          de=db(k+1)/dq
          dalpha(k)=dx+dq+de
   10   continue
        dbeta(1)=db(1)*(da(1)-dx)
        return
      else if(iopt.eq.2) then
        ds=dx-da(1)
        dt=dy
        deio=0.d0
        do 20 k=1,n
          dd=ds*ds+dt*dt
          der=-db(k+1)*ds/dd
          dei=db(k+1)*dt/dd
          ds=dx+der-da(k+1)
          dt=dy+dei
          dalpha(k)=dx+dt*der/dei-ds*dei/dt
          dbeta(k)=dt*deio*(1.d0+(der/dei)**2)
          deio=dei
   20   continue
        dbeta(1)=db(1)*(db(2)+(da(1)-dx)**2+dy*dy)
        return
      else if(iopt.eq.3) then
        dt=dy
        deio=0.d0
        do 30 k=1,n
          dei=db(k+1)/dt
          dt=dy+dei
          dalpha(k)=0.d0
          dbeta(k)=dt*deio
          deio=dei
   30   continue
        dbeta(1)=db(1)*(db(2)+dy*dy)
        return
      else if(iopt.eq.4) then
        dalpha(1)=dx-db(1)/dhr
        dbeta(1)=-dhr
        dq=-db(1)/dhr
        do 40 k=2,n
          de=da(k-1)-dx-dq
          dbeta(k)=dq*de
          dq=db(k)/de
          dalpha(k)=dq+de+dx
   40   continue
        return
      else if(iopt.eq.5) then
        nm1=n-1
        dd=dhr*dhr+dhi*dhi
        deroo=da(1)-dx+db(1)*dhr/dd
        deioo=-db(1)*dhi/dd-dy
        dalpha(1)=dx+dhr*dy/dhi
        dbeta(1)=-dhi/dy
        dalpha(2)=dx-db(1)*dhi*deroo/(dd*deioo)+dhr*deioo/dhi
        dbeta(2)=dy*deioo*(1.d0+(dhr/dhi)**2)
        if(n.eq.2) return
        dso=db(2)/(deroo**2+deioo**2)
        dero=da(2)-dx-dso*deroo
        deio=dso*deioo-dy
        dalpha(3)=dx+deroo*deio/deioo+dso*deioo*dero/deio
        dbeta(3)=-db(1)*dhi*deio*(1.d0+(deroo/deioo)**2)/dd
        if(n.eq.3) return
        do 50 k=3,nm1
          ds=db(k)/(dero**2+deio**2)
          der=da(k)-dx-ds*dero
          dei=ds*deio-dy
          dalpha(k+1)=dx+dero*dei/deio+ds*deio*der/dei
          dbeta(k+1)=dso*deioo*dei*(1.d0+(dero/deio)**2)
          deroo=dero
          deioo=deio
          dero=der
          deio=dei
          dso=ds
   50   continue
        return
      else if(iopt.eq.6) then
        nm1=n-1
        deoo=-db(1)/dhi-dy
        deo=db(2)/deoo-dy
        dalpha(1)=0.d0
        dbeta(1)=-dhi/dy
        dalpha(2)=0.d0
        dbeta(2)=dy*deoo
        if(n.eq.2) return
        dalpha(3)=0.d0
        dbeta(3)=-db(1)*deo/dhi
        if(n.eq.3) return
        do 60 k=3,nm1
          de=db(k)/deo-dy
          dbeta(k+1)=db(k-1)*de/deoo
          dalpha(k+1)=0.d0
          deoo=deo
          deo=de
   60   continue
        return
      else if(iopt.eq.7) then
        du=0.d0
        dc=1.d0
        dc0=0.d0
        do 70 k=1,n
          dgam=da(k)-dx-du
          dcm1=dc0
          dc0=dc
          if(dabs(dc0).gt.deps) then
            dp2=(dgam**2)/dc0
          else
            dp2=dcm1*db(k)
          end if
          if(k.gt.1) dbeta(k)=ds*(dp2+db(k+1))
          ds=db(k+1)/(dp2+db(k+1))
          dc=dp2/(dp2+db(k+1))
          du=ds*(dgam+da(k+1)-dx)
          dalpha(k)=dgam+du+dx
   70   continue
        dbeta(1)=db(1)*(db(2)+(dx-da(1))**2)
        return
      else
        ierr=2
        return
      end if
      end

C-END-OF-FILE


c
      subroutine dknum(n,nu0,numax,dx,dy,deps,da,db,drhor,drhoi,nu,
     *ierr,droldr,droldi)
c
c This is a double-precision version of the routine  knum.
c
      double precision dx,dy,deps,da(numax),db(numax),drhor(*),
     *drhoi(*),droldr(*),droldi(*),drr,dri,dden,dt
c
c The arrays  drhor,drhoi,droldr,droldi  are assumed to have
c dimension  n+1.
c
      ierr=0
      np1=n+1
      if(nu0.gt.numax) then
        ierr=nu0
        return
      end if
      if(nu0.lt.np1) nu0=np1
      nu=nu0-5
      do 10 k=1,np1
        drhor(k)=0.d0
        drhoi(k)=0.d0
   10 continue
   20 nu=nu+5
      if(nu.gt.numax) then
        ierr=numax
        goto 60
      end if
      do 30 k=1,np1
        droldr(k)=drhor(k)
        droldi(k)=drhoi(k)
   30 continue
      drr=0.d0
      dri=0.d0
      do 40 j=1,nu
        j1=nu-j+1
        dden=(dx-da(j1)-drr)**2+(dy-dri)**2
        drr=db(j1)*(dx-da(j1)-drr)/dden
        dri=-db(j1)*(dy-dri)/dden
        if(j1.le.np1) then
          drhor(j1)=drr
          drhoi(j1)=dri
        end if
   40 continue
      do 50 k=1,np1
C Following statement replaced -- authors remark
c
c       if((drhor(k)-droldr(k))**2+(drhoi(k)-droldi(k))**2.gt.
c    *    deps*(drhor(k)**2+drhoi(k)**2)) goto 20
        if((drhor(k)-droldr(k))**2+(drhoi(k)-droldi(k))**2.gt.
     *    (deps**2)*(drhor(k)**2+drhoi(k)**2)) goto 20
   50 continue
   60 if(n.eq.0) return
      do 70 k=2,np1
        dt=drhor(k)*drhor(k-1)-drhoi(k)*drhoi(k-1)
        drhoi(k)=drhor(k)*drhoi(k-1)+drhoi(k)*drhor(k-1)
        drhor(k)=dt
   70 continue
      return
      end 

C-END-OF-FILE


c
      subroutine dkern(n,nu0,numax,dx,dy,deps,da,db,dkerr,dkeri,
     *  nu,ierr,droldr,droldi)
c
c This is a double-precision version of the routine  kern.
c
      double precision dx,dy,deps,da(numax),db(numax),dkerr(*),
     *  dkeri(*),droldr(*),droldi(*),dp0r,dp0i,dpr,dpi,dpm1r,
     *  dpm1i,dden,dt
c
c The arrays  dkerr,dkeri,droldr,droldi  are assumed to have
c dimension  n+1.
c
      call dknum(n,nu0,numax,dx,dy,deps,da,db,dkerr,dkeri,nu,ierr,
     *  droldr,droldi)
      if(ierr.ne.0) return
      dp0r=0.d0
      dp0i=0.d0
      dpr=1.d0
      dpi=0.d0
      do 10 k=1,n
        dpm1r=dp0r
        dpm1i=dp0i
        dp0r=dpr
        dp0i=dpi
        dpr=(dx-da(k))*dp0r-dy*dp0i-db(k)*dpm1r
        dpi=(dx-da(k))*dp0i+dy*dp0r-db(k)*dpm1i
        dden=dpr**2+dpi**2
        dt=(dkerr(k+1)*dpr+dkeri(k+1)*dpi)/dden
        dkeri(k+1)=(dkeri(k+1)*dpr-dkerr(k+1)*dpi)/dden
        dkerr(k+1)=dt
   10 continue
      return
      end

C-END-OF-FILE

cat <<C-END-OF-FILE > gchri.f
c
c

cat <<C-END-OF-FILE > dgchri.f
c
c
      subroutine dgchri(n,iopt,nu0,numax,deps,da,db,dx,dy,dalpha,dbeta,
     *nu,ierr,ierrc,dnu,drhor,drhoi,droldr,droldi,ds,ds0,ds1,ds2)
c
c This is a double-precision version of the routine  gchri.
c
      double precision deps,da(numax),db(numax),dx,dy,dalpha(n),
     *dbeta(n),dnu(*),drhor(*),drhoi(*),droldr(*),droldi(*),
     *ds(n),ds0(*),ds1(*),ds2(*)
c
c The arrays  dnu,drhor,drhoi,droldr,droldi,ds0,ds1,ds2  are assumed
c to have dimension  2*n.
c
      if(n.lt.1) then
        ierr=-1
        return
      end if
      ierr=0
      nd=2*n
      ndm1=nd-1
      if(iopt.eq.1) then
        call dknum(ndm1,nu0,numax,dx,dy,deps,da,db,drhor,drhoi,nu,
     *    ierr,droldr,droldi)
        do 10 k=1,nd
          dnu(k)=-drhor(k)
   10   continue
        call dcheb(n,da,db,dnu,dalpha,dbeta,ds,ierrc,ds0,ds1,ds2)
        return
      else if(iopt.eq.2) then
        dy=dabs(dy)
        call dknum(ndm1,nu0,numax,dx,dy,deps,da,db,drhor,drhoi,nu,
     *    ierr,droldr,droldi)
        do 20 k=1,nd
          dnu(k)=-drhoi(k)/dy
   20   continue
        call dcheb(n,da,db,dnu,dalpha,dbeta,ds,ierrc,ds0,ds1,ds2)
        return
      else
        ierr=1
        return
      end if
      end

C-END-OF-FILE

c ###################################
c #                                 #
c # 6. GAUSS-TYPE QUADRATURE RULES  #
c #                                 #
c ###################################

      subroutine dgauss(n,dalpha,dbeta,deps,dzero,dweigh,ierr,de)
c
c This is a double-precision version of the routine  gauss.
c
      double precision dalpha,dbeta,deps,dzero,dweigh,de,dp,dg,dr,
     *ds,dc,df,db
      dimension dalpha(n),dbeta(n),dzero(n),dweigh(n),de(n)
      if(n.lt.1) then
        ierr=-1
        return
      end if
      ierr=0
      dzero(1)=dalpha(1)
      if(dbeta(1).lt.0.d0) then
        ierr=-2
        return
      end if
      dweigh(1)=dbeta(1)
      if (n.eq.1) return
      dweigh(1)=1.d0
      de(n)=0.d0
      do 100 k=2,n
        dzero(k)=dalpha(k)
        if(dbeta(k).lt.0.d0) then
          ierr=-2
          return
        end if
        de(k-1)=dsqrt(dbeta(k))
        dweigh(k)=0.d0
  100 continue
      do 240 l=1,n
        j=0
  105   do 110 m=l,n
          if(m.eq.n) goto 120
          if(dabs(de(m)).le.deps*(dabs(dzero(m))+dabs(dzero(m+1)))) 
     *      goto 120
  110   continue
  120   dp=dzero(l)
        if(m.eq.l) goto 240
        if(j.eq.30) goto 400
        j=j+1
        dg=(dzero(l+1)-dp)/(2.d0*de(l))
        dr=dsqrt(dg*dg+1.d0)
        dg=dzero(m)-dp+de(l)/(dg+dsign(dr,dg))
        ds=1.d0
        dc=1.d0
        dp=0.d0
        mml=m-l
        do 200 ii=1,mml
          i=m-ii
          df=ds*de(i)
          db=dc*de(i)
          if(dabs(df).lt.dabs(dg)) goto 150
          dc=dg/df
          dr=dsqrt(dc*dc+1.d0)
          de(i+1)=df*dr
          ds=1.d0/dr
          dc=dc*ds
          goto 160
  150     ds=df/dg
          dr=dsqrt(ds*ds+1.d0)
          de(i+1)=dg*dr
          dc=1.d0/dr
          ds=ds*dc
  160     dg=dzero(i+1)-dp
          dr=(dzero(i)-dg)*ds+2.d0*dc*db
          dp=ds*dr
          dzero(i+1)=dg+dp
          dg=dc*dr-db
          df=dweigh(i+1)
          dweigh(i+1)=ds*dweigh(i)+dc*df
          dweigh(i)=dc*dweigh(i)-ds*df
  200   continue
        dzero(l)=dzero(l)-dp
        de(l)=dg
        de(m)=0.d0
        goto 105
  240 continue
      do 300 ii=2,n
        i=ii-1
        k=i
        dp=dzero(i)
        do 260 j=ii,n
          if(dzero(j).ge.dp) goto 260
          k=j
          dp=dzero(j)
  260   continue
        if(k.eq.i) goto 300
        dzero(k)=dzero(i)
        dzero(i)=dp
        dp=dweigh(i)
        dweigh(i)=dweigh(k)
        dweigh(k)=dp
  300 continue
      do 310 k=1,n
        dweigh(k)=dbeta(1)*dweigh(k)*dweigh(k)
  310 continue
      return
  400 ierr=l
      return
      end

C-END-OF-FILE


      subroutine dradau(n,dalpha,dbeta,dend,dzero,dweigh,ierr,de,
     *da,db)
c
c This is a double-precision version of the routine  radau.
c
      double precision dend,depsma,dp0,dp1,dpm1,dalpha(*),dbeta(*),
     *dzero(*),dweigh(*),de(*),da(*),db(*),d1mach
c
c The arrays  dalpha,dbeta,dzero,dweigh,de,da,db  are assumed to have
c dimension  n+1.
c
      depsma=d1mach(3)
c
c depsma is the machine double precision.
c
      np1=n+1
      do 10 k=1,np1
        da(k)=dalpha(k)
        db(k)=dbeta(k)
   10 continue
      dp0=0.d0
      dp1=1.d0
      do 20 k=1,n
        dpm1=dp0
        dp0=dp1
        dp1=(dend-da(k))*dp0-db(k)*dpm1
   20 continue
      da(np1)=dend-db(np1)*dp0/dp1
      call dgauss(np1,da,db,depsma,dzero,dweigh,ierr,de)
      return
      end

C-END-OF-FILE


cat <<C-END-OF-FILE > dlob.f
c
c

      subroutine dlob(n,dalpha,dbeta,dleft,dright,dzero,dweigh,
     *ierr,de,da,db)
c
c This is a double-precision version of the routine  lob.
c
      double precision dleft,dright,depsma,dp0l,dp0r,dp1l,dp1r,dpm1l,
     *dpm1r,ddet,dalpha(*),dbeta(*),dzero(*),dweigh(*),de(*),da(*),
     *db(*),d1mach
c
c The arrays  dalpha,dbeta,dzero,dweigh,de,da,db  are assumed to have
c dimension  n+2.
c
      depsma=d1mach(3)
c
c depsma is the machine double precision.
c
      np1=n+1
      np2=n+2
      do 10 k=1,np2
        da(k)=dalpha(k)
        db(k)=dbeta(k)
   10 continue
      dp0l=0.d0
      dp0r=0.d0
      dp1l=1.d0
      dp1r=1.d0
      do 20 k=1,np1
        dpm1l=dp0l
        dp0l=dp1l
        dpm1r=dp0r
        dp0r=dp1r
        dp1l=(dleft-da(k))*dp0l-db(k)*dpm1l
        dp1r=(dright-da(k))*dp0r-db(k)*dpm1r
   20 continue
      ddet=dp1l*dp0r-dp1r*dp0l
      da(np2)=(dleft*dp1l*dp0r-dright*dp1r*dp0l)/ddet
      db(np2)=(dright-dleft)*dp1l*dp1r/ddet
      call dgauss(np2,da,db,depsma,dzero,dweigh,ierr,de)
      return
      end

C-END-OF-FILE
