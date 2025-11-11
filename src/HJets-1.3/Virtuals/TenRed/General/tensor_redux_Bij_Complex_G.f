CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C                tens_red2 = 2-point tensors
C     Francisco Campanario
C    Email: Francam@particle.uni-karlsruhe.de
C Date: 31/10/2008
C  Tensor reduction a la Passarino-Vetlman for the tensor
C In this file, we should put as an input the B0 scalar invariants, 
C the 2 A0
C Notation: Example BR. the "R" means the real part of B0
C                   BI, the "I" means the imaginary part of B0
C It speeds up the code. 
C OUTPUT:  is given in B0r,B0i,Bijr,BijI. "r":  real part, "I" imaginary part
C Massive case
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
       subroutine tens_red2_Complex(m0sq,m1sq,p1sq,A0r2,A0I2,A0r1,A0I1, 
c     #                     B0,!Cij,
     & B0r,B0I,Bijr,BijI)
      implicit none
      real * 8  p1sq
c      complex*16 B0
c
c  determine the Passarino-Veltman tensor decomposition for the three-point
c  tensor integrals
c
c                       d^4k           1;  k_mu;  
c B0;B_mu=Int   -----------------------------------------
c                     (2pi)^4    [k^2-m0^2][(k+p1)^2-m1^2]
c
c with
c
c    B_mu = p1_mu B1
c
c
C INPUT:  p1sq,              external invariants: p1^2<=>0
C         A0, B0                   3 scalar integrals; the 2 A0 are, 
c                                  in PV notation:
c         A0(1) =A0(m0)  A_0 function with subtraction of 
c         A0(2) =A0(m1)  divergent term

c OUTPUT: B1          form factors in the tensor integrals

      real*8 deter,Ideter
      real*8 A0r2,A0r1,Bijr
      real*8 A0I2,A0I1,BijI
      real*8 B0r 
      real*8 B0I 
      real*8 eps,inv2
      Parameter(eps=1d-7)
      complex*16 m0sq,m1sq!,m0,m1,
      complex*16 mdif,msum,mprud,BijrT,DerB0rT,B0rt
      complex*16 A0r1T,A0r2T

      deter = 2.d0*p1sq

      Ideter=1.d0/deter
      inv2=1d0/2d0
c      m0sq=m0*m0
c      m1sq=m1*m1
      
      B0rT=DCMPLX(B0R,B0I)
      A0r1T=DCMPLX(A0r1,A0I1)
      A0r2T=DCMPLX(A0r2,A0I2)



      If(abs(p1sq).ge.eps) then
         BijrT=((-(p1sq-m1sq)-m0sq)*B0rT+A0r1T-A0r2T)*Ideter

         else  ! p1sq =0
            mdif= m0sq-m1sq
           If(abs(mdif).ge.eps) then  ! m0=!=m1
            msum=m0sq+m1sq
            mprud=m0sq*m1sq
                If(abs(mprud).gt.eps*1d-2) then  ! m0sq>0 & m1sq>0
                  DerB0RT=inv2*msum/(mdif*mdif)+mprud/(mdif*mdif*mdif)*log(m1sq/m0sq)
                else   ! m0sq=0 or m1sq=0
                  DerB0RT=inv2*msum/(mdif*mdif)
                endif
               BijRT=-inv2*(B0rT+mdif*DerB0rT)

           else  ! m0=m1
               BijrT=-inv2*B0rT
           endif
       endif

       BijR=Dble(BijrT)
       BijI=DImag(BijrT)

       return
      end




