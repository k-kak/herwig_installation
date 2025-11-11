       subroutine tens_red3_Complex_G(m0sq,m1sq,m2sq,p1sq,p2sq,s12,
     & B0r_23,B0r_13,B0r_12, 
     & B0I_23,B0I_13,B0I_12, 
     & B1r_23,B1r_13,B1r_12, 
     & B1I_23,B1I_13,B1I_12, 
     & C0,!Cij,
     & C0r,C0I,Cijr,CijI)
C
C
CB0_23,B0_13,B0_12, 
c       subroutine ten_red_LUdecom_G(m0,m1,m2,p1sq,p2sq,s12,B0_23,B0_13,B0_12, 
c     &                     C0,C0r,C0I,Cijr,CijI)
C                tens_red3 = 3-point tensors
C     Francisco Campanario
C    Email: Francam@particle.uni-karlsruhe.de
C Date: 25/02/2010
c
c  determine the Passarino-Veltman tensor decomposition for the three-point
c  tensor integrals
c
c                       d^4k           1;  k_mu;   k_mu k_nu
c C0;C_mu;C_munu =Int -------  -----------------------------------------
c                     (2pi)^4    [k^2-m^2][(k+p1)^2-m^2][(k+p1+p2)^2-m^2] 
c
c with
c
c    C_mu = p1_mu C11  +  p_2_mu C12
c
c  C_munu = p1_mu p1_nu C21 + p2_mu p2_nu C22 + 
c           (p1_mu p2_nu + p1_nu p2_mu) C23  -  g_munu C24
c
c  for notation see Passarino&Veltman, NP B160 (1979) 151 and my notes
c
C INPUT:  p1sq, p2sq, s12          external invariants: p1^2, p2^2, s12=(p1+p2)^2
C         B0, C0                   4 scalar integrals; the 3 B0 are, 
c                                  in PV notation:
c         B0(1) = B0(1,2) = B0(p1)  B_0 function with subtraction of 
c         B0(2) = B0(2,3) = B0(p2)  divergent term
c         B0(3) = B0(1,3) = B0(s12)
c
c OUTPUT: Cij(n,m) = C_nm          form factors in the tensor integrals
c          n=1,2,3,4; n=1,2        a la PV
c

      implicit none
      real * 8  p1sq, p2sq, s12
      complex*16  C0
      real*8 det
      real*8 B1r_12,B1r_13,B1r_23,Cijr(4,2)
      real*8 B1I_12,B1I_13,B1I_23,CijI(4,2)
      real*8 B0r_23, B0r_13, B0r_12, C0r
      real*8 B0I_23, B0I_13, B0I_12, C0I 
      real*8 z11,z12,z21,z22,iz11,iz22
c
      complex*16 r1, r2r1,p1p2,r1r0
      complex*16 Rr(2),RI(2),PRr(2),PRI(2)
      complex*16 m0sq,m1sq,m2sq!m1,m0,m2,
      complex*16 B1r_12T,B1r_13T,B1r_23T,CijrT(4,2)
      complex*16 B0r_23T, B0r_13T, B0r_12T
c
	  real*8 deter,detAbs
	  Logical Singular
	  COMMON Singular
	  
      p1p2 = (s12 - p1sq - p2sq)*0.5d0

c      m0sq=m0*m0
c      m1sq=m1*m1
c      m2sq=m2*m2

      r1 = p1sq -m1sq
      r1r0= r1 + m0sq
      r2r1 = (s12-m2sq) - r1

      deter = abs(2.d0*(p1sq*p2sq - p1p2*p1p2))
      detAbs=  abs(2.d0*(abs(p1sq*p2sq)+abs(p1p2*p1p2)))
      
      If( (deter/detAbs).lt.1d-6) Singular=.true.
  
      C0r=Dble(C0)
      C0I=DImag(C0)

      B0r_13T=DCMPLX(B0r_13,B0I_13)
      B0r_23T=DCMPLX(B0r_23,B0I_23)
      B0r_12T=DCMPLX(B0r_12,B0I_12)

      B1r_13T=DCMPLX(B1r_13,B1I_13)
      B1r_23T=DCMPLX(B1r_23,B1I_23)
      B1r_12T=DCMPLX(B1r_12,B1I_12)

    
      If(abs(p1sq).gt.abs(p1p2)) then
          z11=2d0*p1sq
          iz11=1d0/z11
          z12=2d0*p1p2 
          z21=z12*iz11
          z22=2d0*p2sq-z12*z21
          iz22=1d0/z22
c          iorder(1)=1
c          iorder(2)=2
          det=z11*z22 
c 1-2
       PRr(1) = (B0r_13T - B0r_23T - C0*r1r0)
       PRr(2) = (B0r_12T - B0r_13T - C0*r2r1)
   
       Rr(2)=(PRr(2)-z21*PRr(1))*iz22 
       Rr(1)=(PRr(1)-z12*Rr(2))*iz11 

       CijrT(1,1) = Rr(1)
       CijrT(2,1) = Rr(2)

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
c C00
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

      CijrT(4,2) = ( B0r_23T +2d0*m0sq*C0 + CijrT(1,1)*r1r0 +CijrT(2,1)*r2r1 +1.d0)*0.25d0
     
       PRr(1) = (B1r_13T + B0r_23T - CijrT(1,1)*r1r0 - CijrT(4,2)*2.d0)
       PRr(2) = (B1r_12T - B1r_13T - CijrT(1,1)*r2r1)

       Rr(2)=(PRr(2)-z21*PRr(1))*iz22 
       Rr(1)=(PRr(1)-z12*Rr(2))*iz11 

       CijrT(1,2) =Rr(1)
       CijrT(3,2) =Rr(2)


c 4-6
       PRr(1) = (  B1r_13T - B1r_23T - CijrT(2,1)*r1r0)
       PRr(2) = (- B1r_13T          -CijrT(2,1)*r2r1 -CijrT(4,2)*2.d0)
      
       Rr(2)=(PRr(2)-z21*PRr(1))*iz22 
c      Rr(1)=(PRr(1)-z12*Rr(2))*iz11 

       CijrT(2,2) = Rr(2)
c      CijrT(3,2) = Rr(1)

       
        Cijr(1,1) = Dble(CijrT(1,1))
        Cijr(2,1) = Dble(CijrT(2,1))

        CijI(1,1) =DImag(CijrT(1,1))
        CijI(2,1) =DImag(CijrT(2,1))


        Cijr(1,2) = Dble(CijrT(1,2))
        Cijr(2,2) = Dble(CijrT(2,2))
        Cijr(3,2) = Dble(CijrT(3,2))
        Cijr(4,2) = Dble(CijrT(4,2))


        CijI(1,2) =DImag(CijrT(1,2))
        CijI(2,2) =DImag(CijrT(2,2))
        CijI(3,2) =DImag(CijrT(3,2))
        CijI(4,2) =DImag(CijrT(4,2))



       return
ccccccccccccccccc
cccccccccccccccc
ccccccccccccccccc
       else
		  z11=2d0*p1p2
		  iz11=1d0/z11
		  z21=2d0*p1sq*iz11
		  z12=2d0*p2sq
		  z22=z11-z12*z21
		  iz22=1d0/z22
c          iorder(1)=2
c          iorder(2)=1
          det=-z11*z22

c 1-2
       PRr(1) = (B0r_13T - B0r_23T - C0*r1r0)
       PRr(2) = (B0r_12T - B0r_13T - C0*r2r1)
   
       Rr(2)=(PRr(1)-z21*PRr(2))*iz22 
       Rr(1)=(PRr(2)-z12*Rr(2))*iz11 
 
       CijrT(1,1) = Rr(1)
       CijrT(2,1) = Rr(2)
 

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
c C00
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

      CijrT(4,2) = ( B0r_23T +2d0*m0sq*C0 + CijrT(1,1)*r1r0 +CijrT(2,1)*r2r1 +1.d0)*0.25d0

c 3-5
       PRr(1) = (B1r_13T + B0r_23T - CijrT(1,1)*r1r0 - CijrT(4,2)*2.d0)
       PRr(2) = (B1r_12T - B1r_13T - CijrT(1,1)*r2r1)
       
       Rr(2)=(PRr(1)-z21*PRr(2))*iz22 
       Rr(1)=(PRr(2)-z12*Rr(2))*iz11 
       
       CijrT(1,2) =Rr(1)
       CijrT(3,2) =Rr(2)
       


c 4-6
       PRr(1) = (  B1r_13T - B1r_23T - CijrT(2,1)*r1r0)
       PRr(2) = (- B1r_13T          -CijrT(2,1)*r2r1 -CijrT(4,2)*2.d0)
       
       Rr(2)=(PRr(1)-z21*PRr(2))*iz22 
c      Rr(1)=(PRr(2)-z12*Rr(2))*iz11 

       CijrT(2,2) = Rr(2)
c      CijrT(3,2) = Rr(1)

        Cijr(1,1) = Dble(CijrT(1,1))
        Cijr(2,1) = Dble(CijrT(2,1))

        CijI(1,1) =DImag(CijrT(1,1))
        CijI(2,1) =DImag(CijrT(2,1))


        Cijr(1,2) = Dble(CijrT(1,2))
        Cijr(2,2) = Dble(CijrT(2,2))
        Cijr(3,2) = Dble(CijrT(3,2))
        Cijr(4,2) = Dble(CijrT(4,2))


        CijI(1,2) =DImag(CijrT(1,2))
        CijI(2,2) =DImag(CijrT(2,2))
        CijI(3,2) =DImag(CijrT(3,2))
        CijI(4,2) =DImag(CijrT(4,2))


       endif
	return
	end

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
