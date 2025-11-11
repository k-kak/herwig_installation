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
      subroutine tens_red2_new_Re_ComDiv(m0,m1,p1sq,A0r2,A0I2,A0r1,A0I1, 
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

      real*8 deter,Ideter,DerB0R
      real*8 A0r2,A0r1,Bijr
      real*8 A0I2,A0I1,BijI
      real*8 B0r 
      real*8 B0I 
      real*8 eps
      Parameter(eps=1d-07)
      real*8 m0,m1,m0sq,m1sq

   
      
       Bijr=-1d0/2d0*B0r
       BijI=-1d0/2d0*B0I

       

      end




