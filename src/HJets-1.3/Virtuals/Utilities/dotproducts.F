****************************************************************************
*  file for dot-products and tensor contractions
****************************************************************************
*   This file contains functions for dot-products of real and complex 
*   4-vectors and routines for contractions with the Levi-Civita 
*   epsilon-tensor and arbitrary tensors.
c***************************************************************************
c   LIST OF ALL FUNCTIONS AND SUBROUTINES IN THIS FILE:
c
c   double precision FUNCTION dot0p(p1,p2)
c   double complex   FUNCTION dotqj(j1)
c   real*8           FUNCTION mjj2(p1,p2)
c   real*8           FUNCTION qjj2(p1,p2)
c   real*8           FUNCTION mass2(p)
c   double precision FUNCTION psumsq(p1,p2)
c   double complex   FUNCTION contract_Tjj(T,j1,j2)
c   complex*16       FUNCTION contract_TJJr(T,j1,j2)
c   real*8           FUNCTION contract_Trjjr(T,j1,j2)
c   SUBROUTINE EPSRRR(RES,A,B,C)
c   SUBROUTINE EPSCRR(RES,A,B,C)
c   SUBROUTINE EPSCCR(RES,A,B,C)
c   SUBROUTINE EPSCCC(RES,A,B,C)
c   SUBROUTINE EPSRR (T,A,B)
c   SUBROUTINE contract_CPO(J1cJ2,J1,J2,k1,k2)
c   SUBROUTINE contract_CPE(J1cJ2,J1,J2,k1,k2)
c   SUBROUTINE contract_T1j(T,jc,jout)
c   SUBROUTINE contract_T2j(T,jc,jout)
c   SUBROUTINE contract_Tr1j(T,jc,jout)
c   SUBROUTINE contract_T1jr(T,jr,jout)
c   SUBROUTINE contract_Tr2j(T,jc,jout)
c   SUBROUTINE contract_T2jr(T,jr,jout)
c   SUBROUTINE contract_Tr3j(T,jc,i,Tout)
c   SUBROUTINE contract_Tr3jr(T,jr,i,Tout)
c   SUBROUTINE contract_T3j(T,jc,i,Tout)
c   SUBROUTINE contract_T3jr(T,jr,i,Tout)
c   SUBROUTINE contract_Tr4j(T,jc,i,Tout)
c   SUBROUTINE contract_T4j(T,jc,i,Tout)
c   SUBROUTINE contract_Tr4jr(T,jr,i,Tout)
c   double complex   FUNCTION EPSCRRR(a,b,c,d)
c   double complex   FUNCTION EPSCCRR(a,b,c,d)
c   double precision FUNCTION EPSRRRR(a,b,c,d)
c   double complex   FUNCTION EPSCRRMU(a,b,c,mu)
c   double precision FUNCTION EPSRRRMU(a,b,c,mu)
c   double complex   FUNCTION EPSCRMUNU(a,b,mu,nu)
c   double precision FUNCTION EPSRRMUNU(a,b,mu,nu)
c   double precision FUNCTION EPSRMUNURHO(a,mu,nu,rho)
c   double real*8    FUNCTION EPSMUNURHOSIGMA(mu,nu,rho,sigma)
c
c***************************************************************************
*   Last modified: 06.07.2006

c **************************************************************************
c***************************************************************************
      function EPSRRRR(a,b,c,d)
c***************************************************************************      
c EPSRRRR calculates the the value of a 4 dimensional totally antisymmetric tensor
c EPS-mu-nu-rho-sigma which is contracted with 4 real momenta
c A(mu),B(nu),C(rho),D(sigma)      
      IMPLICIT NONE      
      double precision a(0:3),b(0:3),c(0:3),d(0:3),epsrrrr
      epsrrrr= a(0)*b(1)*c(2)*d(3)-a(0)*b(1)*c(3)*d(2)-a(0)*b(2)*c(1)*d(3)
     &        +a(0)*b(2)*c(3)*d(1)+a(0)*b(3)*c(1)*d(2)-a(0)*b(3)*c(2)*d(1)
     &        -a(1)*b(0)*c(2)*d(3)+a(1)*b(0)*c(3)*d(2)+a(1)*b(2)*c(0)*d(3)
     &        -a(1)*b(2)*c(3)*d(0)-a(1)*b(3)*c(0)*d(2)+a(1)*b(3)*c(2)*d(0)
     &        +a(2)*b(0)*c(1)*d(3)-a(2)*b(0)*c(3)*d(1)-a(2)*b(1)*c(0)*d(3)
     &        +a(2)*b(1)*c(3)*d(0)+a(2)*b(3)*c(0)*d(1)-a(2)*b(3)*c(1)*d(0)
     &        -a(3)*b(0)*c(1)*d(2)+a(3)*b(0)*c(2)*d(1)+a(3)*b(1)*c(0)*d(2)
     &        -a(3)*b(1)*c(2)*d(0)-a(3)*b(2)*c(0)*d(1)+a(3)*b(2)*c(1)*d(0)
      END
c***********
c***************************************************************************
      function EPSCRRR(a,b,c,d)
c****************************************************************************      
c EPSCRRR calculates the the value of a 4 dimensional totally antisymmetric tensor
c EPS-mu-nu-rho-sigma which is contracted with 1 complex current A(mu) and 3 real momenta
c B(nu),C(rho),D(sigma)      
      IMPLICIT NONE
      double complex a(0:3),epscrrr
      double precision b(0:3),c(0:3),d(0:3)
      epscrrr= a(0)*b(1)*c(2)*d(3)-a(0)*b(1)*c(3)*d(2)-a(0)*b(2)*c(1)*d(3)
     &        +a(0)*b(2)*c(3)*d(1)+a(0)*b(3)*c(1)*d(2)-a(0)*b(3)*c(2)*d(1)
     &        -a(1)*b(0)*c(2)*d(3)+a(1)*b(0)*c(3)*d(2)+a(1)*b(2)*c(0)*d(3)
     &        -a(1)*b(2)*c(3)*d(0)-a(1)*b(3)*c(0)*d(2)+a(1)*b(3)*c(2)*d(0)
     &        +a(2)*b(0)*c(1)*d(3)-a(2)*b(0)*c(3)*d(1)-a(2)*b(1)*c(0)*d(3)
     &        +a(2)*b(1)*c(3)*d(0)+a(2)*b(3)*c(0)*d(1)-a(2)*b(3)*c(1)*d(0)
     &        -a(3)*b(0)*c(1)*d(2)+a(3)*b(0)*c(2)*d(1)+a(3)*b(1)*c(0)*d(2)
     &        -a(3)*b(1)*c(2)*d(0)-a(3)*b(2)*c(0)*d(1)+a(3)*b(2)*c(1)*d(0)
      END
c***************************************************************************


c***************************************************************************
      function EPSCCRR(a,b,c,d)
c***************************************************************************      
c EPSCCRR calculates the the value of a 4 dimensional totally antisymmetric tensor
c EPS-mu-nu-rho-sigma which is contracted with 2 complex currents A(mu) and B(nu) and 2 real c momenta C(rho) and D(sigma)      
      IMPLICIT NONE
      double complex a(0:3),b(0:3),epsccrr
      double precision c(0:3),d(0:3)
      epsccrr= a(0)*b(1)*c(2)*d(3)-a(0)*b(1)*c(3)*d(2)-a(0)*b(2)*c(1)*d(3)
     &        +a(0)*b(2)*c(3)*d(1)+a(0)*b(3)*c(1)*d(2)-a(0)*b(3)*c(2)*d(1)
     &        -a(1)*b(0)*c(2)*d(3)+a(1)*b(0)*c(3)*d(2)+a(1)*b(2)*c(0)*d(3)
     &        -a(1)*b(2)*c(3)*d(0)-a(1)*b(3)*c(0)*d(2)+a(1)*b(3)*c(2)*d(0)
     &        +a(2)*b(0)*c(1)*d(3)-a(2)*b(0)*c(3)*d(1)-a(2)*b(1)*c(0)*d(3)
     &        +a(2)*b(1)*c(3)*d(0)+a(2)*b(3)*c(0)*d(1)-a(2)*b(3)*c(1)*d(0)
     &        -a(3)*b(0)*c(1)*d(2)+a(3)*b(0)*c(2)*d(1)+a(3)*b(1)*c(0)*d(2)
     &        -a(3)*b(1)*c(2)*d(0)-a(3)*b(2)*c(0)*d(1)+a(3)*b(2)*c(1)*d(0)
      END      
c***************************************************************************      
c***************************************************************************
      function EPSCCCR(a,b,c,d)
c***************************************************************************      
c EPSCCRR calculates the the value of a 4 dimensional totally antisymmetric tensor
c EPS-mu-nu-rho-sigma which is contracted with 2 complex currents A(mu) and B(nu) and 2 real c momenta C(rho) and D(sigma)      
      IMPLICIT NONE
      double complex a(0:3),b(0:3),c(0:3),epscccr
      double precision d(0:3)
      epscccr= a(0)*b(1)*c(2)*d(3)-a(0)*b(1)*c(3)*d(2)-a(0)*b(2)*c(1)*d(3)
     &        +a(0)*b(2)*c(3)*d(1)+a(0)*b(3)*c(1)*d(2)-a(0)*b(3)*c(2)*d(1)
     &        -a(1)*b(0)*c(2)*d(3)+a(1)*b(0)*c(3)*d(2)+a(1)*b(2)*c(0)*d(3)
     &        -a(1)*b(2)*c(3)*d(0)-a(1)*b(3)*c(0)*d(2)+a(1)*b(3)*c(2)*d(0)
     &        +a(2)*b(0)*c(1)*d(3)-a(2)*b(0)*c(3)*d(1)-a(2)*b(1)*c(0)*d(3)
     &        +a(2)*b(1)*c(3)*d(0)+a(2)*b(3)*c(0)*d(1)-a(2)*b(3)*c(1)*d(0)
     &        -a(3)*b(0)*c(1)*d(2)+a(3)*b(0)*c(2)*d(1)+a(3)*b(1)*c(0)*d(2)
     &        -a(3)*b(1)*c(2)*d(0)-a(3)*b(2)*c(0)*d(1)+a(3)*b(2)*c(1)*d(0)
      END      
c********  
c***************************************************************************
      function EPSCCCC(a,b,c,d)
c***************************************************************************      
c EPSCCRR calculates the the value of a 4 dimensional totally antisymmetric tensor
c EPS-mu-nu-rho-sigma which is contracted with 2 complex currents A(mu) and B(nu) and 2 real c momenta C(rho) and D(sigma)      
      IMPLICIT NONE
      double complex a(0:3),b(0:3),c(0:3),d(0:3),epscccc
      epscccc= a(0)*b(1)*c(2)*d(3)-a(0)*b(1)*c(3)*d(2)-a(0)*b(2)*c(1)*d(3)
     &        +a(0)*b(2)*c(3)*d(1)+a(0)*b(3)*c(1)*d(2)-a(0)*b(3)*c(2)*d(1)
     &        -a(1)*b(0)*c(2)*d(3)+a(1)*b(0)*c(3)*d(2)+a(1)*b(2)*c(0)*d(3)
     &        -a(1)*b(2)*c(3)*d(0)-a(1)*b(3)*c(0)*d(2)+a(1)*b(3)*c(2)*d(0)
     &        +a(2)*b(0)*c(1)*d(3)-a(2)*b(0)*c(3)*d(1)-a(2)*b(1)*c(0)*d(3)
     &        +a(2)*b(1)*c(3)*d(0)+a(2)*b(3)*c(0)*d(1)-a(2)*b(3)*c(1)*d(0)
     &        -a(3)*b(0)*c(1)*d(2)+a(3)*b(0)*c(2)*d(1)+a(3)*b(1)*c(0)*d(2)
     &        -a(3)*b(1)*c(2)*d(0)-a(3)*b(2)*c(0)*d(1)+a(3)*b(2)*c(1)*d(0)
      END      
c********   
c**********
 

