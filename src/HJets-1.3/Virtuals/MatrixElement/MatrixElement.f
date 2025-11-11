c This file contain the functions required to evaluate all the matrix elements 
c that appears in the hexlinexxxx, penlinexxx subroutines and in the boxline subroutines.
c They are extension of the function SC1 in the file brakets.f
c They evaluates matrix elements in the quark line with up to 7 insertion of gamma
c matrices. The sufix added to SC are: "1","3","5","7" number of gamma matrices
c "r" the position for a slash(p), and "c" the position of 
c a gamma matrix that has to be contracted with a external current of polarization vector.
c Alpha is the helicity of the spinors.
c Contain also the funtion randini, and fillmomenta subroutines which generate 
c randon momenta for pp-> WWW and alike processes. The sufix added to fillmomenta "1","3","4","5" indicates the number of massles particles,execpt for historical reasons "1" which set p1^2 and p2^2 to zero.
c Author: Francisco Campanario
c Date: 13/2/2008
c Modified:25/11/2009
c Include SC0,SC2cc,SC2cr,SC2rc,SC2rr,SC3crc,SC4cccc,SC4ccrc,SC4crrc
c SC5crcrc
c Modified:11/12/2009
C Include subroutine SCP2,SCP1,SCP2R,SCP2C,SCM2,SCM1,SCM2R,SCM2C
c Include Massive spinors
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC


      FUNCTION SC0(CHII,CHIF,ALPHA)
      IMPLICIT NONE 
      COMPLEX*16 SC0
      COMPLEX*16  CHII(2), CHIF(2)
      integer alpha
      SC0  = CHII(1)*CHIF(1) + CHII(2)*CHIF(2)
      RETURN
      END

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC


CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC


      FUNCTION SC1c(CHII,A1,CHIF,ALPHA)
      IMPLICIT NONE 
      COMPLEX*16 ZI, SC1c
      PARAMETER ( ZI=(0D0,1D0) )
      INTEGER  ALPHA, ALP, N, I
      COMPLEX*16  CHII(2), CHIF(2), CHIAUX(2), CHIDUM, ASLASH(2,2)
      COMPLEX*16  A1(0:3), AUX(0:3)
   
C
     
      DO I = 0,3
        AUX(I) = A1(I)
      ENDDO
      
      CHIAUX(1) = CHII(1)
      CHIAUX(2) = CHII(2)
      ALP = ALPHA
C
      
         IF (ALP.GT.0) THEN
            ASLASH(1,1) = AUX(0) - AUX(3)
            ASLASH(1,2) = -AUX(1) + ZI*AUX(2)
            ASLASH(2,1) = -AUX(1) - ZI*AUX(2)
            ASLASH(2,2) = AUX(0) + AUX(3)
         ELSE
            ASLASH(1,1) = AUX(0) + AUX(3)
            ASLASH(1,2) = AUX(1) - ZI*AUX(2)
            ASLASH(2,1) = AUX(1) + ZI*AUX(2)
            ASLASH(2,2) = AUX(0) - AUX(3)
         ENDIF
         CHIDUM = CHIAUX(1)*ASLASH(1,1) + CHIAUX(2)*ASLASH(2,1)
         CHIAUX(2) = CHIAUX(1)*ASLASH(1,2) + CHIAUX(2)*ASLASH(2,2)
         CHIAUX(1) = CHIDUM
         ALP = - ALP

C
      SC1c  = CHIAUX(1) * CHIF(1) + CHIAUX(2) * CHIF(2)
      RETURN
      END

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC



CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC


      FUNCTION SC1r(CHII,A1,CHIF,ALPHA)
      IMPLICIT NONE 
      COMPLEX*16 ZI, SC1r
      PARAMETER ( ZI=(0D0,1D0) )
      INTEGER  ALPHA, ALP, N, I
      COMPLEX*16  CHII(2), CHIF(2), CHIAUX(2), CHIDUM, ASLASH(2,2)
      COMPLEX*16  AUX(0:3)
      REAL*8      A1(0:3)   
C
     
      DO I = 0,3
        AUX(I) = A1(I)
      ENDDO
      CHIAUX(1) = CHII(1)
      CHIAUX(2) = CHII(2)
      ALP = ALPHA
C
      
         IF (ALP.GT.0) THEN
            ASLASH(1,1) = AUX(0) - AUX(3)
            ASLASH(1,2) = -AUX(1) + ZI*AUX(2)
            ASLASH(2,1) = -AUX(1) - ZI*AUX(2)
            ASLASH(2,2) = AUX(0) + AUX(3)
         ELSE
            ASLASH(1,1) = AUX(0) + AUX(3)
            ASLASH(1,2) = AUX(1) - ZI*AUX(2)
            ASLASH(2,1) = AUX(1) + ZI*AUX(2)
            ASLASH(2,2) = AUX(0) - AUX(3)
         ENDIF
         CHIDUM = CHIAUX(1)*ASLASH(1,1) + CHIAUX(2)*ASLASH(2,1)
         CHIAUX(2) = CHIAUX(1)*ASLASH(1,2) + CHIAUX(2)*ASLASH(2,2)
         CHIAUX(1) = CHIDUM
         ALP = - ALP

C
      SC1r  = CHIAUX(1) * CHIF(1) + CHIAUX(2) * CHIF(2)
      RETURN
      END



CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
C                         2
C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC




      FUNCTION SC2cc(CHII,A1,A2,CHIF,ALPHA)
      IMPLICIT NONE 
      COMPLEX*16 ZI, SC2cc
      PARAMETER ( ZI=(0D0,1D0) )
      INTEGER  ALPHA, ALP, N, I
      COMPLEX*16  CHII(2), CHIF(2), CHIAUX(2), CHIDUM, ASLASH(2,2)
      COMPLEX*16  AUX(0:3,2),A1(0:3),A2(0:3)
C
      N = 2
      DO I = 0,3
         AUX(I,2) =A2(I)
         AUX(I,1) =A1(I)
      ENDDO
      CHIAUX(1) = CHII(1)
      CHIAUX(2) = CHII(2)
      ALP = ALPHA
C
      
      DO 30 I = 1,N
         IF (ALP.GT.0) THEN
            ASLASH(1,1) = AUX(0,I) - AUX(3,I)
            ASLASH(1,2) = -AUX(1,I) + ZI*AUX(2,I)
            ASLASH(2,1) = -AUX(1,I) - ZI*AUX(2,I)
            ASLASH(2,2) = AUX(0,I) + AUX(3,I)
         ELSE
            ASLASH(1,1) = AUX(0,I) + AUX(3,I)
            ASLASH(1,2) = AUX(1,I) - ZI*AUX(2,I)
            ASLASH(2,1) = AUX(1,I) + ZI*AUX(2,I)
            ASLASH(2,2) = AUX(0,I) - AUX(3,I)
         ENDIF
         CHIDUM = CHIAUX(1)*ASLASH(1,1) + CHIAUX(2)*ASLASH(2,1)
         CHIAUX(2) = CHIAUX(1)*ASLASH(1,2) + CHIAUX(2)*ASLASH(2,2)
         CHIAUX(1) = CHIDUM
         ALP = - ALP
 30   CONTINUE
C
      SC2cc  = CHIAUX(1) * CHIF(1) + CHIAUX(2) * CHIF(2)
      RETURN
      END





CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
C                         2
C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC


      subroutine  SC3P1Mr(M,vM,A,uM,FM)
      IMPLICIT NONE 
      COMPLEX*16 ZI
      PARAMETER (ZI=(0D0,1D0) )
      COMPLEX*16  VM(2),UM(2),FM(0:3,0:3)
      COMPLEX*16 uM1vM1,uM1vM2,uM2vM1,uM2vM2
      COMPLEX*16 aux1,aux2,aux3,aux4
      COMPLEX*16 aux5,aux6,aux7,aux8
      COMPLEX*16 g0,g1,g2,g3
      real*8 m,A(0:3)
      

      uM1vM1= uM(1)*vM(1)
      uM1vM2= uM(1)*vM(2)
      uM2vM1= uM(2)*vM(1)
      uM2vM2= uM(2)*vM(2)
      
      IF(A(3).gt.0d0)then
      g3=A(0)+A(3)
      g0=m*m/g3
      else
      g0=A(0)-A(3)
c      g1=-A(1)+ZI*A(2)
c      g2=-(A(1)+ZI*A(2))
      g3=m*m/g0
      endif

      aux1=g0*uM1vM1 
      aux2=g3*uM2vM2
      aux3=g0*uM2vM1
      aux4=g3*uM1vM2
      aux5=g3*uM2vM1
      aux6=g0*uM1vM2
      aux7=g3*uM1vM1
      aux8=g0*uM2vM2

      FM(0,0)=aux1+aux2
      FM(3,3)=FM(0,0)
      FM(0,1)=aux3+aux4
      FM(3,2)=-ZI*FM(0,1)
      FM(3,1)=aux3-aux4
      FM(0,2)=-ZI*FM(3,1)
      FM(0,3)=aux1-aux2
      FM(3,0)=FM(0,3)
      FM(1,0)=aux5+aux6
      FM(2,3)=ZI*FM(1,0)
      FM(1,1)=aux7+aux8
      FM(2,2)=FM(1,1)
      FM(1,2)=ZI*(aux7-aux8)
      FM(2,1)=-FM(1,2)
      FM(1,3)=aux6-aux5
      FM(2,0)=ZI*FM(1,3)
 
      RETURN
      END


      subroutine  SC3P1Pr(M,vP,A,uP,FP)
      IMPLICIT NONE 
      COMPLEX*16 ZI
      PARAMETER (ZI=(0D0,1D0) )
      COMPLEX*16  VP(2),UP(2),FP(0:3,0:3)
      COMPLEX*16 uP1vP1,uP1vP2,uP2vP1,uP2vP2
      COMPLEX*16 aux1,aux2,aux3,aux4
      COMPLEX*16 aux5,aux6,aux7,aux8
      COMPLEX*16 g0,g3
      real*8 m,A(0:3)

      uP1vP1= uP(1)*vP(1)
      uP1vP2= uP(1)*vP(2)
      uP2vP1= uP(2)*vP(1)
      uP2vP2= uP(2)*vP(2)
      
      IF(A(3).gt.0d0) then
      g3=A(0)+A(3)
      g0=m*m/g3
      else
      g0=A(0)-A(3)
c      g1=-A(1)+ZI*A(2)
c      g2=-(A(1)+ZI*A(2))
      g3=m*m/g0
      endif

      aux1=g3*uP1vP1 
      aux2=g0*uP2vP2
      aux3=g3*uP2vP1
      aux4=g0*uP1vP2
      aux5=g0*uP2vP1
      aux6=g3*uP1vP2
      aux7=g0*uP1vP1
      aux8=g3*uP2vP2



      FP(0,0)=aux1+aux2
      FP(3,3)=FP(0,0)
      FP(0,1)=-(aux3+aux4)
      FP(3,2)=ZI*FP(0,1)
      FP(3,1)=aux3-aux4
      FP(0,2)=ZI*FP(3,1)
      FP(0,3)=aux2-aux1
      FP(3,0)=FP(0,3)
      FP(1,0)=-(aux5+aux6)
      FP(2,3)=-ZI*FP(1,0)
      FP(1,1)=aux7+aux8
      FP(2,2)=FP(1,1)
      FP(1,2)=ZI*(aux7-aux8)
      FP(2,1)=-FP(1,2)
      FP(1,3)=aux6-aux5
      FP(2,0)=-ZI*FP(1,3) 


      RETURN
      END














      subroutine  SC1M(vM,uM,FM)
      IMPLICIT NONE 
      COMPLEX*16 ZI
      PARAMETER (ZI=(0D0,1D0) )
      COMPLEX*16  VM(2),UM(2),FM(0:3)
      COMPLEX*16 uM1vM1,uM1vM2,uM2vM1,uM2vM2
 
      uM1vM1= uM(1)*vM(1)
      uM1vM2= uM(1)*vM(2)
      uM2vM1= uM(2)*vM(1)
      uM2vM2= uM(2)*vM(2)

      FM(0)=uM1vM1 + uM2vM2
      FM(1)=uM2vM1 + uM1vM2
      FM(2)=ZI*(uM1vM2 - uM2vM1)
      FM(3)=uM1vM1 - uM2vM2
   
      RETURN
      END

      subroutine  SC1P(vP,uP,FP)
      IMPLICIT NONE 
      COMPLEX*16 ZI
      PARAMETER (ZI=(0D0,1D0) )
      COMPLEX*16  VP(2),UP(2),FP(0:3)
      COMPLEX*16 uP1vP1,uP1vP2,uP2vP1,uP2vP2
 
      uP1vP1= uP(1)*vP(1)
      uP1vP2= uP(1)*vP(2)
      uP2vP1= uP(2)*vP(1)
      uP2vP2= uP(2)*vP(2)

      FP(0)=uP1vP1 + uP2vP2
      FP(1)=-uP2vP1 - uP1vP2
      FP(2)=ZI*(uP2vP1-uP1vP2)
      FP(3)=uP2vP2-uP1vP1
   
      RETURN
      END


      subroutine  SC2M(vP,uM,FM)
      IMPLICIT NONE 
      COMPLEX*16 ZI
      PARAMETER (ZI=(0D0,1D0) )
      COMPLEX*16  VP(2),UM(2),FM(0:3,0:3)
      COMPLEX*16 uM1vP1,uM1vP2,uM2vP1,uM2vP2
 
      uM1vP1= uM(1)*vP(1)
      uM1vP2= uM(1)*vP(2)
      uM2vP1= uM(2)*vP(1)
      uM2vP2= uM(2)*vP(2)

      FM(0,0)=uM1vP1 + uM2vP2
      FM(1,1)=-FM(0,0)
      FM(2,2)=FM(1,1)
      FM(3,3)=FM(2,2)
      FM(0,1)=uM2vP1 + uM1vP2
      FM(1,0)=-FM(0,1)
      FM(0,2)=ZI*(uM1vP2-uM2vP1)
      FM(2,0)=-FM(0,2)
      FM(0,3)=uM1vP1-uM2vP2
      FM(3,0)=-FM(0,3)
      FM(1,2)=ZI*(uM2vP2-uM1vP1)
      FM(2,1)=-FM(1,2)
      FM(1,3)=uM2vP1-uM1vP2
      FM(3,1)=-FM(1,3)
      FM(2,3)=-ZI*(uM2vP1+uM1vP2)
      FM(3,2)=-FM(2,3)
      RETURN
      END


      subroutine  SC2P(vM,uP,FP)
      IMPLICIT NONE 
      COMPLEX*16 ZI
      PARAMETER (ZI=(0D0,1D0) )
      COMPLEX*16  VM(2),UP(2),FP(0:3,0:3)
      COMPLEX*16 uP1vM1,uP1vM2,uP2vM1,uP2vM2
 
      uP1vM1= uP(1)*vM(1)
      uP1vM2= uP(1)*vM(2)
      uP2vM1= uP(2)*vM(1)
      uP2vM2= uP(2)*vM(2)

      FP(0,0)=uP1vM1 + uP2vM2
      FP(1,1)=-FP(0,0)
      FP(2,2)=FP(1,1)
      FP(3,3)=FP(2,2)
      FP(0,1)=-(uP2vM1 + uP1vM2)
      FP(1,0)=-FP(0,1)
      FP(0,2)=-ZI*(uP1vM2-uP2vM1)
      FP(2,0)=-FP(0,2)
      FP(0,3)=-(uP1vM1-uP2vM2)
      FP(3,0)=-FP(0,3)
      FP(1,2)=ZI*(uP2vM2-uP1vM1)
      FP(2,1)=-FP(1,2)
      FP(1,3)=(uP2vM1-uP1vM2)
      FP(3,1)=-FP(1,3)
      FP(2,3)=-ZI*(uP2vM1+uP1vM2)
      FP(3,2)=-FP(2,3)
 
      RETURN
      END


c  FVD[v1].GAD
      subroutine  SC21Mc(vP,A,uM,FM)
      IMPLICIT NONE 
      COMPLEX*16 ZI
      PARAMETER (ZI=(0D0,1D0) )
      COMPLEX*16  VP(2),UM(2),FM(0:3),A(0:3)
      COMPLEX*16 uM1vP1,uM1vP2,uM2vP1,uM2vP2
      COMPLEX*16 aux1,aux2,aux3,aux4
      COMPLEX*16 g0,g1,g2,g3
      

      uM1vP1= uM(1)*vP(1)
      uM1vP2= uM(1)*vP(2)
      uM2vP1= uM(2)*vP(1)
      uM2vP2= uM(2)*vP(2)
      
      g0=A(0)-A(3)
      g1=-A(1)+ZI*A(2)
      g2=-(A(1)+ZI*A(2))
      g3=A(0)+A(3)

      aux1= g3*uM1vP1-g1*uM2vP1
      aux2=-g2*uM1vP2+g0*uM2vP2
      aux3= g2*uM1vP1-g0*uM2vP1
      aux4=-g3*uM1vP2+g1*uM2vP2

      FM(0)=aux1+aux2
      FM(1)=aux3+aux4
      FM(2)=ZI*(-aux3+aux4)
      FM(3)=-aux1+aux2

      RETURN
      END
      subroutine  SC21MR(vP,A,uM,FM)
      IMPLICIT NONE 
      COMPLEX*16 ZI
      PARAMETER (ZI=(0D0,1D0) )
      COMPLEX*16  VP(2),UM(2),FM(0:3)
      COMPLEX*16 uM1vP1,uM1vP2,uM2vP1,uM2vP2
      COMPLEX*16 aux1,aux2,aux3,aux4
      COMPLEX*16 g0,g1,g2,g3
      REAL*8 A(0:3)

      uM1vP1= uM(1)*vP(1)
      uM1vP2= uM(1)*vP(2)
      uM2vP1= uM(2)*vP(1)
      uM2vP2= uM(2)*vP(2)
      
      g0=A(0)-A(3)
      g1=-A(1)+ZI*A(2)
      g2=-(A(1)+ZI*A(2))
      g3=A(0)+A(3)

      aux1= g3*uM1vP1-g1*uM2vP1
      aux2=-g2*uM1vP2+g0*uM2vP2
      aux3= g2*uM1vP1-g0*uM2vP1
      aux4=-g3*uM1vP2+g1*uM2vP2

      FM(0)=aux1+aux2
      FM(1)=aux3+aux4
      FM(2)=ZI*(-aux3+aux4)
      FM(3)=-aux1+aux2

      RETURN
      END


c  GAD.FVD[v1]
      subroutine  SC2Mc(vP,A,uM,FM)
      IMPLICIT NONE 
      COMPLEX*16 ZI
      PARAMETER (ZI=(0D0,1D0) )
      COMPLEX*16  VP(2),UM(2),FM(0:3),A(0:3)
      COMPLEX*16 uM1vP1,uM1vP2,uM2vP1,uM2vP2
      COMPLEX*16 aux1,aux2,aux3,aux4
      COMPLEX*16 g0,g1,g2,g3
      

      uM1vP1= uM(1)*vP(1)
      uM1vP2= uM(1)*vP(2)
      uM2vP1= uM(2)*vP(1)
      uM2vP2= uM(2)*vP(2)
      
      g0=A(0)-A(3)
      g1=-A(1)+ZI*A(2)
      g2=-(A(1)+ZI*A(2))
      g3=A(0)+A(3)

      aux1= g0*uM1vP1+g2*uM1vP2
      aux2= g1*uM2vP1+g3*uM2vP2
      aux3= g1*uM1vP1+g3*uM1vP2
      aux4= g0*uM2vP1+g2*uM2vP2

      FM(0)=aux1+aux2
      FM(1)=aux3+aux4
      FM(2)=ZI*(aux3-aux4)
      FM(3)=aux1-aux2

      RETURN
      END

      subroutine  SC2Mr(vP,A,uM,FM)
      IMPLICIT NONE 
      COMPLEX*16 ZI
      PARAMETER (ZI=(0D0,1D0) )
      COMPLEX*16  VP(2),UM(2),FM(0:3)
      COMPLEX*16 uM1vP1,uM1vP2,uM2vP1,uM2vP2
      COMPLEX*16 aux1,aux2,aux3,aux4
      COMPLEX*16 g0,g1,g2,g3
      real*8 A(0:3)
      

      uM1vP1= uM(1)*vP(1)
      uM1vP2= uM(1)*vP(2)
      uM2vP1= uM(2)*vP(1)
      uM2vP2= uM(2)*vP(2)
      
      g0=A(0)-A(3)
      g1=-A(1)+ZI*A(2)
      g2=-(A(1)+ZI*A(2))
      g3=A(0)+A(3)

      aux1= g0*uM1vP1+g2*uM1vP2
      aux2= g1*uM2vP1+g3*uM2vP2
      aux3= g1*uM1vP1+g3*uM1vP2
      aux4= g0*uM2vP1+g2*uM2vP2

      FM(0)=aux1+aux2
      FM(1)=aux3+aux4
      FM(2)=ZI*(aux3-aux4)
      FM(3)=aux1-aux2

      RETURN
      END







c  GAD.FVD[v1]
      subroutine  SC2MrM(M,vP,A,uM,FM)
      IMPLICIT NONE 
      COMPLEX*16 ZI
      PARAMETER (ZI=(0D0,1D0) )
      COMPLEX*16  VP(2),UM(2),FM(0:3)
      COMPLEX*16 uM1vP1,uM1vP2,uM2vP1,uM2vP2
      COMPLEX*16 aux1,aux2,aux3,aux4
      COMPLEX*16 g0,g1,g2,g3
      real*8 A(0:3),m
      

      uM1vP1= uM(1)*vP(1)
      uM1vP2= uM(1)*vP(2)
      uM2vP1= uM(2)*vP(1)
      uM2vP2= uM(2)*vP(2)
      
c      g0=A(0)-A(3)
      g1=-A(1)+ZI*A(2)
      g2=-(A(1)+ZI*A(2))
c      g3=A(0)+A(3)

      IF(A(3).gt.0d0) then
      g3=A(0)+A(3)
      g0=m*m/g3
      else
      g0=A(0)-A(3)
c      g1=-A(1)+ZI*A(2)
c      g2=-(A(1)+ZI*A(2))
      g3=m*m/g0
      endif



      aux1= g0*uM1vP1+g2*uM1vP2
      aux2= g1*uM2vP1+g3*uM2vP2
      aux3= g1*uM1vP1+g3*uM1vP2
      aux4= g0*uM2vP1+g2*uM2vP2

      FM(0)=aux1+aux2
      FM(1)=aux3+aux4
      FM(2)=ZI*(aux3-aux4)
      FM(3)=aux1-aux2

      RETURN
      END










cccccccccccccccccccccccccccccccccccccccccccccccc


      subroutine  SC21Pc(vM,A,uP,FP)
      IMPLICIT NONE 
      COMPLEX*16 ZI
      PARAMETER (ZI=(0D0,1D0) )
      COMPLEX*16  VM(2),UP(2),FP(0:3),A(0:3)
      COMPLEX*16 uP1vM1,uP1vM2,uP2vM1,uP2vM2
      COMPLEX*16 aux1,aux2,aux3,aux4
      COMPLEX*16 g0,g1,g2,g3
      

      uP1vM1= uP(1)*vM(1)
      uP1vM2= uP(1)*vM(2)
      uP2vM1= uP(2)*vM(1)
      uP2vM2= uP(2)*vM(2)
      
      g0=A(0)-A(3)
      g1=-A(1)+ZI*A(2)
      g2=-(A(1)+ZI*A(2))
      g3=A(0)+A(3)

      aux1= g0*uP1vM1+g1*uP2vM1
      aux2=+g2*uP1vM2+g3*uP2vM2
      aux3= g2*uP1vM1+g3*uP2vM1
      aux4= g0*uP1vM2+g1*uP2vM2

      FP(0)=aux1+aux2
      FP(1)=aux3+aux4
      FP(2)=ZI*(-aux3+aux4)
      FP(3)=aux1-aux2

      RETURN
      END

      subroutine  SC21Pr(vM,A,uP,FP)
      IMPLICIT NONE 
      COMPLEX*16 ZI
      PARAMETER (ZI=(0D0,1D0) )
      COMPLEX*16  VM(2),UP(2),FP(0:3)
      COMPLEX*16 uP1vM1,uP1vM2,uP2vM1,uP2vM2
      COMPLEX*16 aux1,aux2,aux3,aux4
      COMPLEX*16 g0,g1,g2,g3
      REAL*8 A(0:3)

      uP1vM1= uP(1)*vM(1)
      uP1vM2= uP(1)*vM(2)
      uP2vM1= uP(2)*vM(1)
      uP2vM2= uP(2)*vM(2)
      
      g0=A(0)-A(3)
      g1=-A(1)+ZI*A(2)
      g2=-(A(1)+ZI*A(2))
      g3=A(0)+A(3)

      aux1= g0*uP1vM1+g1*uP2vM1
      aux2=+g2*uP1vM2+g3*uP2vM2
      aux3= g2*uP1vM1+g3*uP2vM1
      aux4= g0*uP1vM2+g1*uP2vM2

      FP(0)=aux1+aux2
      FP(1)=aux3+aux4
      FP(2)=ZI*(-aux3+aux4)
      FP(3)=aux1-aux2

      RETURN
      END

cccccccccccccccccccccccccccccccccccccccccccccccc


      subroutine  SC2Pc(vM,A,uP,FP)
      IMPLICIT NONE 
      COMPLEX*16 ZI
      PARAMETER (ZI=(0D0,1D0) )
      COMPLEX*16  VM(2),UP(2),FP(0:3),A(0:3)
      COMPLEX*16 uP1vM1,uP1vM2,uP2vM1,uP2vM2
      COMPLEX*16 aux1,aux2,aux3,aux4
      COMPLEX*16 g0,g1,g2,g3
      

      uP1vM1= uP(1)*vM(1)
      uP1vM2= uP(1)*vM(2)
      uP2vM1= uP(2)*vM(1)
      uP2vM2= uP(2)*vM(2)
      
      g0=A(0)-A(3)
      g1=-A(1)+ZI*A(2)
      g2=-(A(1)+ZI*A(2))
      g3=A(0)+A(3)

      aux1= g3*uP1vM1-g2*uP1vM2
      aux2=-g1*uP2vM1+g0*uP2vM2
      aux3= g1*uP1vM1-g0*uP1vM2
      aux4=-g3*uP2vM1+g2*uP2vM2

      FP(0)=aux1+aux2
      FP(1)=aux3+aux4
      FP(2)=ZI*(aux3-aux4)
      FP(3)=-aux1+aux2

      RETURN
      END

      subroutine  SC2Pr(vM,A,uP,FP)
      IMPLICIT NONE 
      COMPLEX*16 ZI
      PARAMETER (ZI=(0D0,1D0) )
      COMPLEX*16  VM(2),UP(2),FP(0:3)
      COMPLEX*16 uP1vM1,uP1vM2,uP2vM1,uP2vM2
      COMPLEX*16 aux1,aux2,aux3,aux4
      COMPLEX*16 g0,g1,g2,g3
      REAL*8 A(0:3)

      uP1vM1= uP(1)*vM(1)
      uP1vM2= uP(1)*vM(2)
      uP2vM1= uP(2)*vM(1)
      uP2vM2= uP(2)*vM(2)
      
      g0=A(0)-A(3)
      g1=-A(1)+ZI*A(2)
      g2=-(A(1)+ZI*A(2))
      g3=A(0)+A(3)

 
      aux1= g3*uP1vM1-g2*uP1vM2
      aux2=-g1*uP2vM1+g0*uP2vM2
      aux3= g1*uP1vM1-g0*uP1vM2
      aux4=-g3*uP2vM1+g2*uP2vM2

      FP(0)=aux1+aux2
      FP(1)=aux3+aux4
      FP(2)=ZI*(aux3-aux4)
      FP(3)=-aux1+aux2

      RETURN
      END



cccccccccccccccccccccccccccccccccccccccccccccccc

      subroutine  SC2PrM(M,vM,A,uP,FP)
      IMPLICIT NONE 
      COMPLEX*16 ZI
      PARAMETER (ZI=(0D0,1D0) )
      COMPLEX*16  VM(2),UP(2),FP(0:3)
      COMPLEX*16 uP1vM1,uP1vM2,uP2vM1,uP2vM2
      COMPLEX*16 aux1,aux2,aux3,aux4
      COMPLEX*16 g0,g1,g2,g3
      REAL*8 A(0:3),M

      uP1vM1= uP(1)*vM(1)
      uP1vM2= uP(1)*vM(2)
      uP2vM1= uP(2)*vM(1)
      uP2vM2= uP(2)*vM(2)
      
c      g0=A(0)-A(3)
      g1=-A(1)+ZI*A(2)
      g2=-(A(1)+ZI*A(2))
c      g3=A(0)+A(3)

       IF(A(3).gt.0d0) then
      g3=A(0)+A(3)
      g0=m*m/g3
      else
      g0=A(0)-A(3)
c      g1=-A(1)+ZI*A(2)
c      g2=-(A(1)+ZI*A(2))
      g3=m*m/g0
      endif

      aux1= g3*uP1vM1-g2*uP1vM2
      aux2=-g1*uP2vM1+g0*uP2vM2
      aux3= g1*uP1vM1-g0*uP1vM2
      aux4=-g3*uP2vM1+g2*uP2vM2

      FP(0)=aux1+aux2
      FP(1)=aux3+aux4
      FP(2)=ZI*(aux3-aux4)
      FP(3)=-aux1+aux2

      RETURN
      END












CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC


      FUNCTION SC2cr(CHII,A1,A2,CHIF,ALPHA)
      IMPLICIT NONE 
      COMPLEX*16 ZI, SC2cr
      PARAMETER ( ZI=(0D0,1D0) )
      INTEGER  ALPHA, ALP, N, I
      COMPLEX*16  CHII(2), CHIF(2), CHIAUX(2), CHIDUM, ASLASH(2,2)
      COMPLEX*16  AUX(0:3,2),A1(0:3)
      REAL*8   A2(0:3)


C
      N = 2
      DO I = 0,3
         AUX(I,2) =A2(I)
         AUX(I,1) =A1(I)
      ENDDO
      CHIAUX(1) = CHII(1)
      CHIAUX(2) = CHII(2)
      ALP = ALPHA
C
      
      DO 30 I = 1,N
         IF (ALP.GT.0) THEN
            ASLASH(1,1) = AUX(0,I) - AUX(3,I)
            ASLASH(1,2) = -AUX(1,I) + ZI*AUX(2,I)
            ASLASH(2,1) = -AUX(1,I) - ZI*AUX(2,I)
            ASLASH(2,2) = AUX(0,I) + AUX(3,I)
         ELSE
            ASLASH(1,1) = AUX(0,I) + AUX(3,I)
            ASLASH(1,2) = AUX(1,I) - ZI*AUX(2,I)
            ASLASH(2,1) = AUX(1,I) + ZI*AUX(2,I)
            ASLASH(2,2) = AUX(0,I) - AUX(3,I)
         ENDIF
         CHIDUM = CHIAUX(1)*ASLASH(1,1) + CHIAUX(2)*ASLASH(2,1)
         CHIAUX(2) = CHIAUX(1)*ASLASH(1,2) + CHIAUX(2)*ASLASH(2,2)
         CHIAUX(1) = CHIDUM
         ALP = - ALP
 30   CONTINUE
C
      SC2cr  = CHIAUX(1) * CHIF(1) + CHIAUX(2) * CHIF(2)
      RETURN
      END

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

      FUNCTION SC2rc(CHII,A1,A2,CHIF,ALPHA)
      IMPLICIT NONE 
      COMPLEX*16 ZI, SC2rc
      PARAMETER ( ZI=(0D0,1D0) )
      INTEGER  ALPHA, ALP, N, I
      COMPLEX*16  CHII(2), CHIF(2), CHIAUX(2), CHIDUM, ASLASH(2,2)
      COMPLEX*16  AUX(0:3,2), A2(0:3)
      REAL*8      A1(0:3)


C
      N = 2
      DO I = 0,3
         AUX(I,2) =A2(I)
         AUX(I,1) =A1(I)
      ENDDO
      CHIAUX(1) = CHII(1)
      CHIAUX(2) = CHII(2)
      ALP = ALPHA
C
      
      DO 30 I = 1,N
         IF (ALP.GT.0) THEN
            ASLASH(1,1) = AUX(0,I) - AUX(3,I)
            ASLASH(1,2) = -AUX(1,I) + ZI*AUX(2,I)
            ASLASH(2,1) = -AUX(1,I) - ZI*AUX(2,I)
            ASLASH(2,2) = AUX(0,I) + AUX(3,I)
         ELSE
            ASLASH(1,1) = AUX(0,I) + AUX(3,I)
            ASLASH(1,2) = AUX(1,I) - ZI*AUX(2,I)
            ASLASH(2,1) = AUX(1,I) + ZI*AUX(2,I)
            ASLASH(2,2) = AUX(0,I) - AUX(3,I)
         ENDIF
         CHIDUM = CHIAUX(1)*ASLASH(1,1) + CHIAUX(2)*ASLASH(2,1)
         CHIAUX(2) = CHIAUX(1)*ASLASH(1,2) + CHIAUX(2)*ASLASH(2,2)
         CHIAUX(1) = CHIDUM
         ALP = - ALP
 30   CONTINUE
C
      SC2rc  = CHIAUX(1) * CHIF(1) + CHIAUX(2) * CHIF(2)
      RETURN
      END


CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      FUNCTION SC2rr(CHII,A1,A2,CHIF,ALPHA)
      IMPLICIT NONE 
      COMPLEX*16 ZI, SC2rr
      PARAMETER ( ZI=(0D0,1D0) )
      INTEGER  ALPHA, ALP, N, I
      COMPLEX*16  CHII(2), CHIF(2), CHIAUX(2), CHIDUM, ASLASH(2,2)
      COMPLEX*16  AUX(0:3,2)
      REAL*8      A1(0:3), A2(0:3)


C
      N = 2
      DO I = 0,3
         AUX(I,2) =A2(I)
         AUX(I,1) =A1(I)
      ENDDO
      CHIAUX(1) = CHII(1)
      CHIAUX(2) = CHII(2)
      ALP = ALPHA
C
      
      DO 30 I = 1,N
         IF (ALP.GT.0) THEN
            ASLASH(1,1) = AUX(0,I) - AUX(3,I)
            ASLASH(1,2) = -AUX(1,I) + ZI*AUX(2,I)
            ASLASH(2,1) = -AUX(1,I) - ZI*AUX(2,I)
            ASLASH(2,2) = AUX(0,I) + AUX(3,I)
         ELSE
            ASLASH(1,1) = AUX(0,I) + AUX(3,I)
            ASLASH(1,2) = AUX(1,I) - ZI*AUX(2,I)
            ASLASH(2,1) = AUX(1,I) + ZI*AUX(2,I)
            ASLASH(2,2) = AUX(0,I) - AUX(3,I)
         ENDIF
         CHIDUM = CHIAUX(1)*ASLASH(1,1) + CHIAUX(2)*ASLASH(2,1)
         CHIAUX(2) = CHIAUX(1)*ASLASH(1,2) + CHIAUX(2)*ASLASH(2,2)
         CHIAUX(1) = CHIDUM
         ALP = - ALP
 30   CONTINUE
C
      SC2rr  = CHIAUX(1) * CHIF(1) + CHIAUX(2) * CHIF(2)
      RETURN
      END


CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C                             3
C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC     




      FUNCTION SC3ccc(CHII,A1,A2,A3,CHIF,ALPHA)
      IMPLICIT NONE 
      COMPLEX*16 ZI, SC3ccc
      PARAMETER ( ZI=(0D0,1D0) )
      INTEGER  ALPHA, ALP, N, I
      COMPLEX*16  CHII(2), CHIF(2), CHIAUX(2), CHIDUM, ASLASH(2,2)
      COMPLEX*16  A1(0:3), A2(0:3), A3(0:3), AUX(0:3,3)
     
C
      N = 3
      DO I = 0,3
         AUX(I,3) = A3(I)
         AUX(I,2) = A2(I)
         AUX(I,1) = A1(I)
      ENDDO
      CHIAUX(1) = CHII(1)
      CHIAUX(2) = CHII(2)
      ALP = ALPHA
C
      DO 30 I = 1,N
         IF (ALP.GT.0) THEN
            ASLASH(1,1) = AUX(0,I) - AUX(3,I)
            ASLASH(1,2) = -AUX(1,I) + ZI*AUX(2,I)
            ASLASH(2,1) = -AUX(1,I) - ZI*AUX(2,I)
            ASLASH(2,2) = AUX(0,I) + AUX(3,I)
         ELSE
            ASLASH(1,1) = AUX(0,I) + AUX(3,I)
            ASLASH(1,2) = AUX(1,I) - ZI*AUX(2,I)
            ASLASH(2,1) = AUX(1,I) + ZI*AUX(2,I)
            ASLASH(2,2) = AUX(0,I) - AUX(3,I)
         ENDIF
         CHIDUM = CHIAUX(1)*ASLASH(1,1) + CHIAUX(2)*ASLASH(2,1)
         CHIAUX(2) = CHIAUX(1)*ASLASH(1,2) + CHIAUX(2)*ASLASH(2,2)
         CHIAUX(1) = CHIDUM
         ALP = - ALP
 30   CONTINUE
C
      SC3ccc  = CHIAUX(1) * CHIF(1) + CHIAUX(2) * CHIF(2)
      RETURN
      END

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC



      FUNCTION SC3ccr(CHII,A1,A2,A3,CHIF,ALPHA)
      IMPLICIT NONE 
      COMPLEX*16 ZI, SC3ccr
      PARAMETER ( ZI=(0D0,1D0) )
      INTEGER  ALPHA, ALP, N, I
      COMPLEX*16  CHII(2), CHIF(2), CHIAUX(2), CHIDUM, ASLASH(2,2)
      COMPLEX*16  A2(0:3), A1(0:3), AUX(0:3,3)
      REAL*8  A3(0:3)
C
      N = 3
      DO I = 0,3
         AUX(I,3) = A3(I)
         AUX(I,2) = A2(I)
         AUX(I,1) = A1(I)
      ENDDO
      CHIAUX(1) = CHII(1)
      CHIAUX(2) = CHII(2)
      ALP = ALPHA
C
      DO 30 I = 1,N
         IF (ALP.GT.0) THEN
            ASLASH(1,1) = AUX(0,I) - AUX(3,I)
            ASLASH(1,2) = -AUX(1,I) + ZI*AUX(2,I)
            ASLASH(2,1) = -AUX(1,I) - ZI*AUX(2,I)
            ASLASH(2,2) = AUX(0,I) + AUX(3,I)
         ELSE
            ASLASH(1,1) = AUX(0,I) + AUX(3,I)
            ASLASH(1,2) = AUX(1,I) - ZI*AUX(2,I)
            ASLASH(2,1) = AUX(1,I) + ZI*AUX(2,I)
            ASLASH(2,2) = AUX(0,I) - AUX(3,I)
         ENDIF
         CHIDUM = CHIAUX(1)*ASLASH(1,1) + CHIAUX(2)*ASLASH(2,1)
         CHIAUX(2) = CHIAUX(1)*ASLASH(1,2) + CHIAUX(2)*ASLASH(2,2)
         CHIAUX(1) = CHIDUM
         ALP = - ALP
 30   CONTINUE
C
      SC3ccr  = CHIAUX(1) * CHIF(1) + CHIAUX(2) * CHIF(2)
      RETURN
      END



CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC



      FUNCTION SC3crc(CHII,A1,A2,A3,CHIF,ALPHA)
      IMPLICIT NONE 
      COMPLEX*16 ZI, SC3crc
      PARAMETER ( ZI=(0D0,1D0) )
      INTEGER  ALPHA, ALP, N, I
      COMPLEX*16  CHII(2), CHIF(2), CHIAUX(2), CHIDUM, ASLASH(2,2)
      COMPLEX*16  A1(0:3), AUX(0:3,3),A3(0:3)
      REAL*8  A2(0:3)
C
      N = 3
      DO I = 0,3
         AUX(I,3) = A3(I)
         AUX(I,2) = A2(I)
         AUX(I,1) = A1(I)
      ENDDO
      CHIAUX(1) = CHII(1)
      CHIAUX(2) = CHII(2)
      ALP = ALPHA
C
      DO 30 I = 1,N
         IF (ALP.GT.0) THEN
            ASLASH(1,1) = AUX(0,I) - AUX(3,I)
            ASLASH(1,2) = -AUX(1,I) + ZI*AUX(2,I)
            ASLASH(2,1) = -AUX(1,I) - ZI*AUX(2,I)
            ASLASH(2,2) = AUX(0,I) + AUX(3,I)
         ELSE
            ASLASH(1,1) = AUX(0,I) + AUX(3,I)
            ASLASH(1,2) = AUX(1,I) - ZI*AUX(2,I)
            ASLASH(2,1) = AUX(1,I) + ZI*AUX(2,I)
            ASLASH(2,2) = AUX(0,I) - AUX(3,I)
         ENDIF
         CHIDUM = CHIAUX(1)*ASLASH(1,1) + CHIAUX(2)*ASLASH(2,1)
         CHIAUX(2) = CHIAUX(1)*ASLASH(1,2) + CHIAUX(2)*ASLASH(2,2)
         CHIAUX(1) = CHIDUM
         ALP = - ALP
 30   CONTINUE
C
      SC3crc  = CHIAUX(1) * CHIF(1) + CHIAUX(2) * CHIF(2)
      RETURN
      END



CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC



      FUNCTION SC3crr(CHII,A1,A2,A3,CHIF,ALPHA)
      IMPLICIT NONE 
      COMPLEX*16 ZI, SC3crr
      PARAMETER ( ZI=(0D0,1D0) )
      INTEGER  ALPHA, ALP, N, I
      COMPLEX*16  CHII(2), CHIF(2), CHIAUX(2), CHIDUM, ASLASH(2,2)
      COMPLEX*16  A1(0:3), AUX(0:3,3)
      REAL*8  A2(0:3),A3(0:3)
C
      N = 3
      DO I = 0,3
         AUX(I,3) = A3(I)
         AUX(I,2) = A2(I)
         AUX(I,1) = A1(I)
      ENDDO
      CHIAUX(1) = CHII(1)
      CHIAUX(2) = CHII(2)
      ALP = ALPHA
C
      DO 30 I = 1,N
         IF (ALP.GT.0) THEN
            ASLASH(1,1) = AUX(0,I) - AUX(3,I)
            ASLASH(1,2) = -AUX(1,I) + ZI*AUX(2,I)
            ASLASH(2,1) = -AUX(1,I) - ZI*AUX(2,I)
            ASLASH(2,2) = AUX(0,I) + AUX(3,I)
         ELSE
            ASLASH(1,1) = AUX(0,I) + AUX(3,I)
            ASLASH(1,2) = AUX(1,I) - ZI*AUX(2,I)
            ASLASH(2,1) = AUX(1,I) + ZI*AUX(2,I)
            ASLASH(2,2) = AUX(0,I) - AUX(3,I)
         ENDIF
         CHIDUM = CHIAUX(1)*ASLASH(1,1) + CHIAUX(2)*ASLASH(2,1)
         CHIAUX(2) = CHIAUX(1)*ASLASH(1,2) + CHIAUX(2)*ASLASH(2,2)
         CHIAUX(1) = CHIDUM
         ALP = - ALP
 30   CONTINUE
C
      SC3crr  = CHIAUX(1) * CHIF(1) + CHIAUX(2) * CHIF(2)
      RETURN
      END


CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC




      FUNCTION SC3rcc(CHII,A1,A2,A3,CHIF,ALPHA)
      IMPLICIT NONE 
      COMPLEX*16 ZI, SC3rcc
      PARAMETER ( ZI=(0D0,1D0) )
      INTEGER  ALPHA, ALP, N, I
      COMPLEX*16  CHII(2), CHIF(2), CHIAUX(2), CHIDUM, ASLASH(2,2)
      COMPLEX*16  A2(0:3), A3(0:3), AUX(0:3,3)
      REAL*8  A1(0:3)
C
      N = 3
      DO I = 0,3
         AUX(I,3) = A3(I)
         AUX(I,2) = A2(I)
         AUX(I,1) = A1(I)
      ENDDO
      CHIAUX(1) = CHII(1)
      CHIAUX(2) = CHII(2)
      ALP = ALPHA
C
      DO 30 I = 1,N
         IF (ALP.GT.0) THEN
            ASLASH(1,1) = AUX(0,I) - AUX(3,I)
            ASLASH(1,2) = -AUX(1,I) + ZI*AUX(2,I)
            ASLASH(2,1) = -AUX(1,I) - ZI*AUX(2,I)
            ASLASH(2,2) = AUX(0,I) + AUX(3,I)
         ELSE
            ASLASH(1,1) = AUX(0,I) + AUX(3,I)
            ASLASH(1,2) = AUX(1,I) - ZI*AUX(2,I)
            ASLASH(2,1) = AUX(1,I) + ZI*AUX(2,I)
            ASLASH(2,2) = AUX(0,I) - AUX(3,I)
         ENDIF
         CHIDUM = CHIAUX(1)*ASLASH(1,1) + CHIAUX(2)*ASLASH(2,1)
         CHIAUX(2) = CHIAUX(1)*ASLASH(1,2) + CHIAUX(2)*ASLASH(2,2)
         CHIAUX(1) = CHIDUM
         ALP = - ALP
 30   CONTINUE
C
      SC3rcc  = CHIAUX(1) * CHIF(1) + CHIAUX(2) * CHIF(2)
      RETURN
      END

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC



CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC



      FUNCTION SC3rrc(CHII,A1,A2,A3,CHIF,ALPHA)
      IMPLICIT NONE 
      COMPLEX*16 ZI, SC3rrc
      PARAMETER ( ZI=(0D0,1D0) )
      INTEGER  ALPHA, ALP, N, I
      COMPLEX*16  CHII(2), CHIF(2), CHIAUX(2), CHIDUM, ASLASH(2,2)
      COMPLEX*16  A3(0:3), AUX(0:3,3)
      REAL*8      A1(0:3), A2(0:3)
C
      N = 3
      DO I = 0,3
         AUX(I,3) = A3(I)
         AUX(I,2) = A2(I)
         AUX(I,1) = A1(I)
      ENDDO
      CHIAUX(1) = CHII(1)
      CHIAUX(2) = CHII(2)
      ALP = ALPHA
C
      DO 30 I = 1,N
         IF (ALP.GT.0) THEN
            ASLASH(1,1) = AUX(0,I) - AUX(3,I)
            ASLASH(1,2) = -AUX(1,I) + ZI*AUX(2,I)
            ASLASH(2,1) = -AUX(1,I) - ZI*AUX(2,I)
            ASLASH(2,2) = AUX(0,I) + AUX(3,I)
         ELSE
            ASLASH(1,1) = AUX(0,I) + AUX(3,I)
            ASLASH(1,2) = AUX(1,I) - ZI*AUX(2,I)
            ASLASH(2,1) = AUX(1,I) + ZI*AUX(2,I)
            ASLASH(2,2) = AUX(0,I) - AUX(3,I)
         ENDIF
         CHIDUM = CHIAUX(1)*ASLASH(1,1) + CHIAUX(2)*ASLASH(2,1)
         CHIAUX(2) = CHIAUX(1)*ASLASH(1,2) + CHIAUX(2)*ASLASH(2,2)
         CHIAUX(1) = CHIDUM
         ALP = - ALP
 30   CONTINUE
C
      SC3rrc  = CHIAUX(1) * CHIF(1) + CHIAUX(2) * CHIF(2)
      RETURN
      END


CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC




CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      

      FUNCTION SC3rrr(CHII,A1,A2,A3,CHIF,ALPHA)
      IMPLICIT NONE 
      COMPLEX*16 ZI, SC3rrr
      PARAMETER ( ZI=(0D0,1D0) )
      INTEGER  ALPHA, ALP, N, I
      COMPLEX*16  CHII(2), CHIF(2), CHIAUX(2), CHIDUM, ASLASH(2,2)
      COMPLEX*16  AUX(0:3,3)
      REAL*8      A1(0:3), A2(0:3), A3(0:3)

c      print*, CHII(1),A1(1),A2(1),A3(1),CHIF(1),ALPHA
C
      N = 3
      DO I = 0,3
         AUX(I,3) =A3(I)
         AUX(I,2) =A2(I)
         AUX(I,1) =A1(I)
      ENDDO
      CHIAUX(1) = CHII(1)
      CHIAUX(2) = CHII(2)
      ALP = ALPHA
C
      
      DO 30 I = 1,N
         IF (ALP.GT.0) THEN
            ASLASH(1,1) = AUX(0,I) - AUX(3,I)
            ASLASH(1,2) = -AUX(1,I) + ZI*AUX(2,I)
            ASLASH(2,1) = -AUX(1,I) - ZI*AUX(2,I)
            ASLASH(2,2) = AUX(0,I) + AUX(3,I)
         ELSE
            ASLASH(1,1) = AUX(0,I) + AUX(3,I)
            ASLASH(1,2) = AUX(1,I) - ZI*AUX(2,I)
            ASLASH(2,1) = AUX(1,I) + ZI*AUX(2,I)
            ASLASH(2,2) = AUX(0,I) - AUX(3,I)
         ENDIF
         CHIDUM = CHIAUX(1)*ASLASH(1,1) + CHIAUX(2)*ASLASH(2,1)
         CHIAUX(2) = CHIAUX(1)*ASLASH(1,2) + CHIAUX(2)*ASLASH(2,2)
         CHIAUX(1) = CHIDUM
         ALP = - ALP
 30   CONTINUE
C
      SC3rrr  = CHIAUX(1) * CHIF(1) + CHIAUX(2) * CHIF(2)
      RETURN
      END

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
CCCC                 4
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC


CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      
      FUNCTION SC4cccc(CHII,A1,A2,A3,A4,CHIF,ALPHA)
      IMPLICIT NONE 
      COMPLEX*16 ZI, SC4cccc
      PARAMETER ( ZI=(0D0,1D0) )
      INTEGER  ALPHA, ALP, N, I
      COMPLEX*16  CHII(2), CHIF(2), CHIAUX(2), CHIDUM, ASLASH(2,2)
      COMPLEX*16  A1(0:3), A2(0:3), A3(0:3), A4(0:3), AUX(0:3,4)
     
C
      N = 4
      DO I = 0,3
         AUX(I,4) = A4(I)
         AUX(I,3) = A3(I)
         AUX(I,2) = A2(I)
         AUX(I,1) = A1(I)
      ENDDO
      CHIAUX(1) = CHII(1)
      CHIAUX(2) = CHII(2)
      ALP = ALPHA
C
      DO 30 I = 1,N
         IF (ALP.GT.0) THEN
            ASLASH(1,1) = AUX(0,I) - AUX(3,I)
            ASLASH(1,2) = -AUX(1,I) + ZI*AUX(2,I)
            ASLASH(2,1) = -AUX(1,I) - ZI*AUX(2,I)
            ASLASH(2,2) = AUX(0,I) + AUX(3,I)
         ELSE
            ASLASH(1,1) = AUX(0,I) + AUX(3,I)
            ASLASH(1,2) = AUX(1,I) - ZI*AUX(2,I)
            ASLASH(2,1) = AUX(1,I) + ZI*AUX(2,I)
            ASLASH(2,2) = AUX(0,I) - AUX(3,I)
         ENDIF
         CHIDUM = CHIAUX(1)*ASLASH(1,1) + CHIAUX(2)*ASLASH(2,1)
         CHIAUX(2) = CHIAUX(1)*ASLASH(1,2) + CHIAUX(2)*ASLASH(2,2)
         CHIAUX(1) = CHIDUM
         ALP = - ALP
 30   CONTINUE
C
      SC4cccc  = CHIAUX(1) * CHIF(1) + CHIAUX(2) * CHIF(2)
      RETURN
      END
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      
      FUNCTION SC4ccrc(CHII,A1,A2,A3,A4,CHIF,ALPHA)
      IMPLICIT NONE 
      COMPLEX*16 ZI, SC4ccrc
      PARAMETER ( ZI=(0D0,1D0) )
      INTEGER  ALPHA, ALP, N, I
      COMPLEX*16  CHII(2), CHIF(2), CHIAUX(2), CHIDUM, ASLASH(2,2)
      COMPLEX*16  A1(0:3), A2(0:3), A4(0:3), AUX(0:3,4)
      REAL*8 A3(0:3)
C
      N = 4
      DO I = 0,3
         AUX(I,4) = A4(I)
         AUX(I,3) = A3(I)
         AUX(I,2) = A2(I)
         AUX(I,1) = A1(I)
      ENDDO
      CHIAUX(1) = CHII(1)
      CHIAUX(2) = CHII(2)
      ALP = ALPHA
C
      DO 30 I = 1,N
         IF (ALP.GT.0) THEN
            ASLASH(1,1) = AUX(0,I) - AUX(3,I)
            ASLASH(1,2) = -AUX(1,I) + ZI*AUX(2,I)
            ASLASH(2,1) = -AUX(1,I) - ZI*AUX(2,I)
            ASLASH(2,2) = AUX(0,I) + AUX(3,I)
         ELSE
            ASLASH(1,1) = AUX(0,I) + AUX(3,I)
            ASLASH(1,2) = AUX(1,I) - ZI*AUX(2,I)
            ASLASH(2,1) = AUX(1,I) + ZI*AUX(2,I)
            ASLASH(2,2) = AUX(0,I) - AUX(3,I)
         ENDIF
         CHIDUM = CHIAUX(1)*ASLASH(1,1) + CHIAUX(2)*ASLASH(2,1)
         CHIAUX(2) = CHIAUX(1)*ASLASH(1,2) + CHIAUX(2)*ASLASH(2,2)
         CHIAUX(1) = CHIDUM
         ALP = - ALP
 30   CONTINUE
C
      SC4ccrc  = CHIAUX(1) * CHIF(1) + CHIAUX(2) * CHIF(2)
      RETURN
      END



CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      
      FUNCTION SC4crrc(CHII,A1,A2,A3,A4,CHIF,ALPHA)
      IMPLICIT NONE 
      COMPLEX*16 ZI, SC4crrc
      PARAMETER ( ZI=(0D0,1D0) )
      INTEGER  ALPHA, ALP, N, I
      COMPLEX*16  CHII(2), CHIF(2), CHIAUX(2), CHIDUM, ASLASH(2,2)
      COMPLEX*16  A1(0:3), A4(0:3), AUX(0:3,4)
      REAL*8 A2(0:3), A3(0:3)
C
      N = 4
      DO I = 0,3
         AUX(I,4) = A4(I)
         AUX(I,3) = A3(I)
         AUX(I,2) = A2(I)
         AUX(I,1) = A1(I)
      ENDDO
      CHIAUX(1) = CHII(1)
      CHIAUX(2) = CHII(2)
      ALP = ALPHA
C
      DO 30 I = 1,N
         IF (ALP.GT.0) THEN
            ASLASH(1,1) = AUX(0,I) - AUX(3,I)
            ASLASH(1,2) = -AUX(1,I) + ZI*AUX(2,I)
            ASLASH(2,1) = -AUX(1,I) - ZI*AUX(2,I)
            ASLASH(2,2) = AUX(0,I) + AUX(3,I)
         ELSE
            ASLASH(1,1) = AUX(0,I) + AUX(3,I)
            ASLASH(1,2) = AUX(1,I) - ZI*AUX(2,I)
            ASLASH(2,1) = AUX(1,I) + ZI*AUX(2,I)
            ASLASH(2,2) = AUX(0,I) - AUX(3,I)
         ENDIF
         CHIDUM = CHIAUX(1)*ASLASH(1,1) + CHIAUX(2)*ASLASH(2,1)
         CHIAUX(2) = CHIAUX(1)*ASLASH(1,2) + CHIAUX(2)*ASLASH(2,2)
         CHIAUX(1) = CHIDUM
         ALP = - ALP
 30   CONTINUE
C
      SC4crrc  = CHIAUX(1) * CHIF(1) + CHIAUX(2) * CHIF(2)
      RETURN
      END







CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
CCCC                 5
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      
      FUNCTION SC5ccccc(CHII,A1,A2,A3,A4,A5,CHIF,ALPHA)
      IMPLICIT NONE 
      COMPLEX*16 ZI, SC5ccccc
      PARAMETER ( ZI=(0D0,1D0) )
      INTEGER  ALPHA, ALP, N, I
      COMPLEX*16  CHII(2), CHIF(2), CHIAUX(2), CHIDUM, ASLASH(2,2)
      COMPLEX*16  A1(0:3), A2(0:3), A3(0:3), A4(0:3),  A5(0:3), AUX(0:3,5)

     
C
      N = 5
      DO I = 0,3
         AUX(I,5) = A5(I)
         AUX(I,4) = A4(I)
         AUX(I,3) = A3(I)
         AUX(I,2) = A2(I)
         AUX(I,1) = A1(I)
      ENDDO
      CHIAUX(1) = CHII(1)
      CHIAUX(2) = CHII(2)
      ALP = ALPHA
C
      DO 30 I = 1,N
         IF (ALP.GT.0) THEN
            ASLASH(1,1) = AUX(0,I) - AUX(3,I)
            ASLASH(1,2) = -AUX(1,I) + ZI*AUX(2,I)
            ASLASH(2,1) = -AUX(1,I) - ZI*AUX(2,I)
            ASLASH(2,2) = AUX(0,I) + AUX(3,I)
         ELSE
            ASLASH(1,1) = AUX(0,I) + AUX(3,I)
            ASLASH(1,2) = AUX(1,I) - ZI*AUX(2,I)
            ASLASH(2,1) = AUX(1,I) + ZI*AUX(2,I)
            ASLASH(2,2) = AUX(0,I) - AUX(3,I)
         ENDIF
         CHIDUM = CHIAUX(1)*ASLASH(1,1) + CHIAUX(2)*ASLASH(2,1)
         CHIAUX(2) = CHIAUX(1)*ASLASH(1,2) + CHIAUX(2)*ASLASH(2,2)
         CHIAUX(1) = CHIDUM
         ALP = - ALP
 30   CONTINUE
C
      SC5ccccc  = CHIAUX(1) * CHIF(1) + CHIAUX(2) * CHIF(2)
      RETURN
      END


CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      
      FUNCTION SC5ccccr(CHII,A1,A2,A3,A4,A5,CHIF,ALPHA)
      IMPLICIT NONE 
      COMPLEX*16 ZI, SC5ccccr
      PARAMETER ( ZI=(0D0,1D0) )
      INTEGER  ALPHA, ALP, N, I
      COMPLEX*16  CHII(2), CHIF(2), CHIAUX(2), CHIDUM, ASLASH(2,2)
      COMPLEX*16  A1(0:3), A2(0:3), A3(0:3), A4(0:3), AUX(0:3,5)
      REAL*8   A5(0:3)
     
C
      N = 5
      DO I = 0,3
         AUX(I,5) = A5(I)
         AUX(I,4) = A4(I)
         AUX(I,3) = A3(I)
         AUX(I,2) = A2(I)
         AUX(I,1) = A1(I)
      ENDDO
      CHIAUX(1) = CHII(1)
      CHIAUX(2) = CHII(2)
      ALP = ALPHA
C
      DO 30 I = 1,N
         IF (ALP.GT.0) THEN
            ASLASH(1,1) = AUX(0,I) - AUX(3,I)
            ASLASH(1,2) = -AUX(1,I) + ZI*AUX(2,I)
            ASLASH(2,1) = -AUX(1,I) - ZI*AUX(2,I)
            ASLASH(2,2) = AUX(0,I) + AUX(3,I)
         ELSE
            ASLASH(1,1) = AUX(0,I) + AUX(3,I)
            ASLASH(1,2) = AUX(1,I) - ZI*AUX(2,I)
            ASLASH(2,1) = AUX(1,I) + ZI*AUX(2,I)
            ASLASH(2,2) = AUX(0,I) - AUX(3,I)
         ENDIF
         CHIDUM = CHIAUX(1)*ASLASH(1,1) + CHIAUX(2)*ASLASH(2,1)
         CHIAUX(2) = CHIAUX(1)*ASLASH(1,2) + CHIAUX(2)*ASLASH(2,2)
         CHIAUX(1) = CHIDUM
         ALP = - ALP
 30   CONTINUE
C
      SC5ccccr  = CHIAUX(1) * CHIF(1) + CHIAUX(2) * CHIF(2)
      RETURN
      END



CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      
      FUNCTION SC5cccrc(CHII,A1,A2,A3,A4,A5,CHIF,ALPHA)
      IMPLICIT NONE 
      COMPLEX*16 ZI, SC5cccrc
      PARAMETER ( ZI=(0D0,1D0) )
      INTEGER  ALPHA, ALP, N, I
      COMPLEX*16  CHII(2), CHIF(2), CHIAUX(2), CHIDUM, ASLASH(2,2)
      COMPLEX*16  A1(0:3), A2(0:3),A3(0:3), A5(0:3) ,AUX(0:3,5)
      REAL*8    A4(0:3)
     
C
      N = 5
      DO I = 0,3
         AUX(I,5) = A5(I)
         AUX(I,4) = A4(I)
         AUX(I,3) = A3(I)
         AUX(I,2) = A2(I)
         AUX(I,1) = A1(I)
      ENDDO
      CHIAUX(1) = CHII(1)
      CHIAUX(2) = CHII(2)
      ALP = ALPHA
C
      DO 30 I = 1,N
         IF (ALP.GT.0) THEN
            ASLASH(1,1) = AUX(0,I) - AUX(3,I)
            ASLASH(1,2) = -AUX(1,I) + ZI*AUX(2,I)
            ASLASH(2,1) = -AUX(1,I) - ZI*AUX(2,I)
            ASLASH(2,2) = AUX(0,I) + AUX(3,I)
         ELSE
            ASLASH(1,1) = AUX(0,I) + AUX(3,I)
            ASLASH(1,2) = AUX(1,I) - ZI*AUX(2,I)
            ASLASH(2,1) = AUX(1,I) + ZI*AUX(2,I)
            ASLASH(2,2) = AUX(0,I) - AUX(3,I)
         ENDIF
         CHIDUM = CHIAUX(1)*ASLASH(1,1) + CHIAUX(2)*ASLASH(2,1)
         CHIAUX(2) = CHIAUX(1)*ASLASH(1,2) + CHIAUX(2)*ASLASH(2,2)
         CHIAUX(1) = CHIDUM
         ALP = - ALP
 30   CONTINUE
C
      SC5cccrc  = CHIAUX(1) * CHIF(1) + CHIAUX(2) * CHIF(2)
      RETURN
      END























CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      
      FUNCTION SC5cccrr(CHII,A1,A2,A3,A4,A5,CHIF,ALPHA)
      IMPLICIT NONE 
      COMPLEX*16 ZI, SC5cccrr
      PARAMETER ( ZI=(0D0,1D0) )
      INTEGER  ALPHA, ALP, N, I
      COMPLEX*16  CHII(2), CHIF(2), CHIAUX(2), CHIDUM, ASLASH(2,2)
      COMPLEX*16  A1(0:3), A2(0:3), A3(0:3),  AUX(0:3,5)
      REAL*8   A4(0:3),A5(0:3)
     
C
      N = 5
      DO I = 0,3
         AUX(I,5) = A5(I)
         AUX(I,4) = A4(I)
         AUX(I,3) = A3(I)
         AUX(I,2) = A2(I)
         AUX(I,1) = A1(I)
      ENDDO
      CHIAUX(1) = CHII(1)
      CHIAUX(2) = CHII(2)
      ALP = ALPHA
C
      DO 30 I = 1,N
         IF (ALP.GT.0) THEN
            ASLASH(1,1) = AUX(0,I) - AUX(3,I)
            ASLASH(1,2) = -AUX(1,I) + ZI*AUX(2,I)
            ASLASH(2,1) = -AUX(1,I) - ZI*AUX(2,I)
            ASLASH(2,2) = AUX(0,I) + AUX(3,I)
         ELSE
            ASLASH(1,1) = AUX(0,I) + AUX(3,I)
            ASLASH(1,2) = AUX(1,I) - ZI*AUX(2,I)
            ASLASH(2,1) = AUX(1,I) + ZI*AUX(2,I)
            ASLASH(2,2) = AUX(0,I) - AUX(3,I)
         ENDIF
         CHIDUM = CHIAUX(1)*ASLASH(1,1) + CHIAUX(2)*ASLASH(2,1)
         CHIAUX(2) = CHIAUX(1)*ASLASH(1,2) + CHIAUX(2)*ASLASH(2,2)
         CHIAUX(1) = CHIDUM
         ALP = - ALP
 30   CONTINUE
C
      SC5cccrr  = CHIAUX(1) * CHIF(1) + CHIAUX(2) * CHIF(2)
      RETURN
      END






CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      
      FUNCTION SC5ccrrc(CHII,A1,A2,A3,A4,A5,CHIF,ALPHA)
      IMPLICIT NONE 
      COMPLEX*16 ZI, SC5ccrrc
      PARAMETER ( ZI=(0D0,1D0) )
      INTEGER  ALPHA, ALP, N, I
      COMPLEX*16  CHII(2), CHIF(2), CHIAUX(2), CHIDUM, ASLASH(2,2)
      COMPLEX*16  A1(0:3), A2(0:3), A5(0:3) ,AUX(0:3,5)
      REAL*8    A3(0:3),A4(0:3)
     
C
      N = 5
      DO I = 0,3
         AUX(I,5) = A5(I)
         AUX(I,4) = A4(I)
         AUX(I,3) = A3(I)
         AUX(I,2) = A2(I)
         AUX(I,1) = A1(I)
      ENDDO
      CHIAUX(1) = CHII(1)
      CHIAUX(2) = CHII(2)
      ALP = ALPHA
C
      DO 30 I = 1,N
         IF (ALP.GT.0) THEN
            ASLASH(1,1) = AUX(0,I) - AUX(3,I)
            ASLASH(1,2) = -AUX(1,I) + ZI*AUX(2,I)
            ASLASH(2,1) = -AUX(1,I) - ZI*AUX(2,I)
            ASLASH(2,2) = AUX(0,I) + AUX(3,I)
         ELSE
            ASLASH(1,1) = AUX(0,I) + AUX(3,I)
            ASLASH(1,2) = AUX(1,I) - ZI*AUX(2,I)
            ASLASH(2,1) = AUX(1,I) + ZI*AUX(2,I)
            ASLASH(2,2) = AUX(0,I) - AUX(3,I)
         ENDIF
         CHIDUM = CHIAUX(1)*ASLASH(1,1) + CHIAUX(2)*ASLASH(2,1)
         CHIAUX(2) = CHIAUX(1)*ASLASH(1,2) + CHIAUX(2)*ASLASH(2,2)
         CHIAUX(1) = CHIDUM
         ALP = - ALP
 30   CONTINUE
C
      SC5ccrrc  = CHIAUX(1) * CHIF(1) + CHIAUX(2) * CHIF(2)
      RETURN
      END












CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      
      FUNCTION SC5ccrrr(CHII,A1,A2,A3,A4,A5,CHIF,ALPHA)
      IMPLICIT NONE 
      COMPLEX*16 ZI, SC5ccrrr
      PARAMETER ( ZI=(0D0,1D0) )
      INTEGER  ALPHA, ALP, N, I
      COMPLEX*16  CHII(2), CHIF(2), CHIAUX(2), CHIDUM, ASLASH(2,2)
      COMPLEX*16  A1(0:3), A2(0:3),  AUX(0:3,5)
      REAL*8    A3(0:3),A4(0:3),A5(0:3)
     
C
      N = 5
      DO I = 0,3
         AUX(I,5) = A5(I)
         AUX(I,4) = A4(I)
         AUX(I,3) = A3(I)
         AUX(I,2) = A2(I)
         AUX(I,1) = A1(I)
      ENDDO
      CHIAUX(1) = CHII(1)
      CHIAUX(2) = CHII(2)
      ALP = ALPHA
C
      DO 30 I = 1,N
         IF (ALP.GT.0) THEN
            ASLASH(1,1) = AUX(0,I) - AUX(3,I)
            ASLASH(1,2) = -AUX(1,I) + ZI*AUX(2,I)
            ASLASH(2,1) = -AUX(1,I) - ZI*AUX(2,I)
            ASLASH(2,2) = AUX(0,I) + AUX(3,I)
         ELSE
            ASLASH(1,1) = AUX(0,I) + AUX(3,I)
            ASLASH(1,2) = AUX(1,I) - ZI*AUX(2,I)
            ASLASH(2,1) = AUX(1,I) + ZI*AUX(2,I)
            ASLASH(2,2) = AUX(0,I) - AUX(3,I)
         ENDIF
         CHIDUM = CHIAUX(1)*ASLASH(1,1) + CHIAUX(2)*ASLASH(2,1)
         CHIAUX(2) = CHIAUX(1)*ASLASH(1,2) + CHIAUX(2)*ASLASH(2,2)
         CHIAUX(1) = CHIDUM
         ALP = - ALP
 30   CONTINUE
C
      SC5ccrrr  = CHIAUX(1) * CHIF(1) + CHIAUX(2) * CHIF(2)
      RETURN
      END
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      
      FUNCTION SC5crccc(CHII,A1,A2,A3,A4,A5,CHIF,ALPHA)
      IMPLICIT NONE 
      COMPLEX*16 ZI, SC5crccc
      PARAMETER ( ZI=(0D0,1D0) )
      INTEGER  ALPHA, ALP, N, I
      COMPLEX*16  CHII(2), CHIF(2), CHIAUX(2), CHIDUM, ASLASH(2,2)
      COMPLEX*16  A1(0:3), A3(0:3), A5(0:3), A4(0:3), AUX(0:3,5)
      REAL*8   A2(0:3)
     
C
      N = 5
      DO I = 0,3
         AUX(I,5) = A5(I)
         AUX(I,4) = A4(I)
         AUX(I,3) = A3(I)
         AUX(I,2) = A2(I)
         AUX(I,1) = A1(I)
      ENDDO
      CHIAUX(1) = CHII(1)
      CHIAUX(2) = CHII(2)
      ALP = ALPHA
C
      DO 30 I = 1,N
         IF (ALP.GT.0) THEN
            ASLASH(1,1) = AUX(0,I) - AUX(3,I)
            ASLASH(1,2) = -AUX(1,I) + ZI*AUX(2,I)
            ASLASH(2,1) = -AUX(1,I) - ZI*AUX(2,I)
            ASLASH(2,2) = AUX(0,I) + AUX(3,I)
         ELSE
            ASLASH(1,1) = AUX(0,I) + AUX(3,I)
            ASLASH(1,2) = AUX(1,I) - ZI*AUX(2,I)
            ASLASH(2,1) = AUX(1,I) + ZI*AUX(2,I)
            ASLASH(2,2) = AUX(0,I) - AUX(3,I)
         ENDIF
         CHIDUM = CHIAUX(1)*ASLASH(1,1) + CHIAUX(2)*ASLASH(2,1)
         CHIAUX(2) = CHIAUX(1)*ASLASH(1,2) + CHIAUX(2)*ASLASH(2,2)
         CHIAUX(1) = CHIDUM
         ALP = - ALP
 30   CONTINUE
C
      SC5crccc  = CHIAUX(1) * CHIF(1) + CHIAUX(2) * CHIF(2)
      RETURN
      END

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      
      FUNCTION SC5crcrc(CHII,A1,A2,A3,A4,A5,CHIF,ALPHA)
      IMPLICIT NONE 
      COMPLEX*16 ZI, SC5crcrc
      PARAMETER ( ZI=(0D0,1D0) )
      INTEGER  ALPHA, ALP, N, I
      COMPLEX*16  CHII(2), CHIF(2), CHIAUX(2), CHIDUM, ASLASH(2,2)
      COMPLEX*16  A1(0:3), A3(0:3), A5(0:3), AUX(0:3,5)
      REAL*8   A2(0:3), A4(0:3)
     
C
      N = 5
      DO I = 0,3
         AUX(I,5) = A5(I)
         AUX(I,4) = A4(I)
         AUX(I,3) = A3(I)
         AUX(I,2) = A2(I)
         AUX(I,1) = A1(I)
      ENDDO
      CHIAUX(1) = CHII(1)
      CHIAUX(2) = CHII(2)
      ALP = ALPHA
C
      DO 30 I = 1,N
         IF (ALP.GT.0) THEN
            ASLASH(1,1) = AUX(0,I) - AUX(3,I)
            ASLASH(1,2) = -AUX(1,I) + ZI*AUX(2,I)
            ASLASH(2,1) = -AUX(1,I) - ZI*AUX(2,I)
            ASLASH(2,2) = AUX(0,I) + AUX(3,I)
         ELSE
            ASLASH(1,1) = AUX(0,I) + AUX(3,I)
            ASLASH(1,2) = AUX(1,I) - ZI*AUX(2,I)
            ASLASH(2,1) = AUX(1,I) + ZI*AUX(2,I)
            ASLASH(2,2) = AUX(0,I) - AUX(3,I)
         ENDIF
         CHIDUM = CHIAUX(1)*ASLASH(1,1) + CHIAUX(2)*ASLASH(2,1)
         CHIAUX(2) = CHIAUX(1)*ASLASH(1,2) + CHIAUX(2)*ASLASH(2,2)
         CHIAUX(1) = CHIDUM
         ALP = - ALP
 30   CONTINUE
C
      SC5crcrc  = CHIAUX(1) * CHIF(1) + CHIAUX(2) * CHIF(2)
      RETURN
      END









CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      
      FUNCTION SC5crrrc(CHII,A1,A2,A3,A4,A5,CHIF,ALPHA)
      IMPLICIT NONE 
      COMPLEX*16 ZI, SC5crrrc
      PARAMETER ( ZI=(0D0,1D0) )
      INTEGER  ALPHA, ALP, N, I
      COMPLEX*16  CHII(2), CHIF(2), CHIAUX(2), CHIDUM, ASLASH(2,2)
      COMPLEX*16  A1(0:3), A5(0:3) ,AUX(0:3,5)
      REAL*8    A2(0:3), A3(0:3),A4(0:3)
     
C
      N = 5
      DO I = 0,3
         AUX(I,5) = A5(I)
         AUX(I,4) = A4(I)
         AUX(I,3) = A3(I)
         AUX(I,2) = A2(I)
         AUX(I,1) = A1(I)
      ENDDO
      CHIAUX(1) = CHII(1)
      CHIAUX(2) = CHII(2)
      ALP = ALPHA
C
      DO 30 I = 1,N
         IF (ALP.GT.0) THEN
            ASLASH(1,1) = AUX(0,I) - AUX(3,I)
            ASLASH(1,2) = -AUX(1,I) + ZI*AUX(2,I)
            ASLASH(2,1) = -AUX(1,I) - ZI*AUX(2,I)
            ASLASH(2,2) = AUX(0,I) + AUX(3,I)
         ELSE
            ASLASH(1,1) = AUX(0,I) + AUX(3,I)
            ASLASH(1,2) = AUX(1,I) - ZI*AUX(2,I)
            ASLASH(2,1) = AUX(1,I) + ZI*AUX(2,I)
            ASLASH(2,2) = AUX(0,I) - AUX(3,I)
         ENDIF
         CHIDUM = CHIAUX(1)*ASLASH(1,1) + CHIAUX(2)*ASLASH(2,1)
         CHIAUX(2) = CHIAUX(1)*ASLASH(1,2) + CHIAUX(2)*ASLASH(2,2)
         CHIAUX(1) = CHIDUM
         ALP = - ALP
 30   CONTINUE
C
      SC5crrrc  = CHIAUX(1) * CHIF(1) + CHIAUX(2) * CHIF(2)
      RETURN
      END




















CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      
      FUNCTION SC5rrccc(CHII,A1,A2,A3,A4,A5,CHIF,ALPHA)
      IMPLICIT NONE 
      COMPLEX*16 ZI, SC5rrccc
      PARAMETER ( ZI=(0D0,1D0) )
      INTEGER  ALPHA, ALP, N, I
      COMPLEX*16  CHII(2), CHIF(2), CHIAUX(2), CHIDUM, ASLASH(2,2)
      COMPLEX*16  A3(0:3), A4(0:3), A5(0:3), AUX(0:3,5)
      REAL*8      A1(0:3), A2(0:3)
     
C
      N = 5
      DO I = 0,3
         AUX(I,5) = A5(I)
         AUX(I,4) = A4(I)
         AUX(I,3) = A3(I)
         AUX(I,2) = A2(I)
         AUX(I,1) = A1(I)
      ENDDO
      CHIAUX(1) = CHII(1)
      CHIAUX(2) = CHII(2)
      ALP = ALPHA
C
      DO 30 I = 1,N
         IF (ALP.GT.0) THEN
            ASLASH(1,1) = AUX(0,I) - AUX(3,I)
            ASLASH(1,2) = -AUX(1,I) + ZI*AUX(2,I)
            ASLASH(2,1) = -AUX(1,I) - ZI*AUX(2,I)
            ASLASH(2,2) = AUX(0,I) + AUX(3,I)
         ELSE
            ASLASH(1,1) = AUX(0,I) + AUX(3,I)
            ASLASH(1,2) = AUX(1,I) - ZI*AUX(2,I)
            ASLASH(2,1) = AUX(1,I) + ZI*AUX(2,I)
            ASLASH(2,2) = AUX(0,I) - AUX(3,I)
         ENDIF
         CHIDUM = CHIAUX(1)*ASLASH(1,1) + CHIAUX(2)*ASLASH(2,1)
         CHIAUX(2) = CHIAUX(1)*ASLASH(1,2) + CHIAUX(2)*ASLASH(2,2)
         CHIAUX(1) = CHIDUM
         ALP = - ALP
 30   CONTINUE
C
      SC5rrccc  = CHIAUX(1) * CHIF(1) + CHIAUX(2) * CHIF(2)
      RETURN
      END


CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      
      FUNCTION SC5rrcrc(CHII,A1,A2,A3,A4,A5,CHIF,ALPHA)
      IMPLICIT NONE 
      COMPLEX*16 ZI, SC5rrcrc
      PARAMETER ( ZI=(0D0,1D0) )
      INTEGER  ALPHA, ALP, N, I
      COMPLEX*16  CHII(2), CHIF(2), CHIAUX(2), CHIDUM, ASLASH(2,2)
      COMPLEX*16  A3(0:3), A5(0:3), AUX(0:3,5)
      REAL*8      A1(0:3), A2(0:3), A4(0:3)
     
C
      N = 5
      DO I = 0,3
         AUX(I,5) = A5(I)
         AUX(I,4) = A4(I)
         AUX(I,3) = A3(I)
         AUX(I,2) = A2(I)
         AUX(I,1) = A1(I)
      ENDDO
      CHIAUX(1) = CHII(1)
      CHIAUX(2) = CHII(2)
      ALP = ALPHA
C
      DO 30 I = 1,N
         IF (ALP.GT.0) THEN
            ASLASH(1,1) = AUX(0,I) - AUX(3,I)
            ASLASH(1,2) = -AUX(1,I) + ZI*AUX(2,I)
            ASLASH(2,1) = -AUX(1,I) - ZI*AUX(2,I)
            ASLASH(2,2) = AUX(0,I) + AUX(3,I)
         ELSE
            ASLASH(1,1) = AUX(0,I) + AUX(3,I)
            ASLASH(1,2) = AUX(1,I) - ZI*AUX(2,I)
            ASLASH(2,1) = AUX(1,I) + ZI*AUX(2,I)
            ASLASH(2,2) = AUX(0,I) - AUX(3,I)
         ENDIF
         CHIDUM = CHIAUX(1)*ASLASH(1,1) + CHIAUX(2)*ASLASH(2,1)
         CHIAUX(2) = CHIAUX(1)*ASLASH(1,2) + CHIAUX(2)*ASLASH(2,2)
         CHIAUX(1) = CHIDUM
         ALP = - ALP
 30   CONTINUE
C
      SC5rrcrc  = CHIAUX(1) * CHIF(1) + CHIAUX(2) * CHIF(2)
      RETURN
      END











CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      
      FUNCTION SC7ccccrrr(CHII,A1,A2,A3,A4,A5,A6,A7,CHIF,ALPHA)
      IMPLICIT NONE 
      COMPLEX*16 ZI, SC7ccccrrr
      PARAMETER ( ZI=(0D0,1D0) )
      INTEGER  ALPHA, ALP, N, I
      COMPLEX*16  CHII(2), CHIF(2), CHIAUX(2), CHIDUM, ASLASH(2,2)
      COMPLEX*16  A1(0:3), A2(0:3), A3(0:3), A4(0:3), AUX(0:3,7)
      REAL*8      A5(0:3), A6(0:3), A7(0:3)
     
C
      N = 7
      DO I = 0,3
         AUX(I,7) = A7(I)
         AUX(I,6) = A6(I)
         AUX(I,5) = A5(I)
         AUX(I,4) = A4(I)
         AUX(I,3) = A3(I)
         AUX(I,2) = A2(I)
         AUX(I,1) = A1(I)
      ENDDO
      CHIAUX(1) = CHII(1)
      CHIAUX(2) = CHII(2)
      ALP = ALPHA
C
      DO 30 I = 1,N
         IF (ALP.GT.0) THEN
            ASLASH(1,1) = AUX(0,I) - AUX(3,I)
            ASLASH(1,2) = -AUX(1,I) + ZI*AUX(2,I)
            ASLASH(2,1) = -AUX(1,I) - ZI*AUX(2,I)
            ASLASH(2,2) = AUX(0,I) + AUX(3,I)
         ELSE
            ASLASH(1,1) = AUX(0,I) + AUX(3,I)
            ASLASH(1,2) = AUX(1,I) - ZI*AUX(2,I)
            ASLASH(2,1) = AUX(1,I) + ZI*AUX(2,I)
            ASLASH(2,2) = AUX(0,I) - AUX(3,I)
         ENDIF
         CHIDUM = CHIAUX(1)*ASLASH(1,1) + CHIAUX(2)*ASLASH(2,1)
         CHIAUX(2) = CHIAUX(1)*ASLASH(1,2) + CHIAUX(2)*ASLASH(2,2)
         CHIAUX(1) = CHIDUM
         ALP = - ALP
 30   CONTINUE
C
      SC7ccccrrr  = CHIAUX(1) * CHIF(1) + CHIAUX(2) * CHIF(2)
      RETURN
      END


CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      
      FUNCTION SC7ccrrccc(CHII,A1,A2,A3,A4,A5,A6,A7,CHIF,ALPHA)
      IMPLICIT NONE 
      COMPLEX*16 ZI, SC7ccrrccc
      PARAMETER ( ZI=(0D0,1D0) )
      INTEGER  ALPHA, ALP, N, I
      COMPLEX*16  CHII(2), CHIF(2), CHIAUX(2), CHIDUM, ASLASH(2,2)
      COMPLEX*16  A1(0:3), A2(0:3), A5(0:3), A6(0:3), A7(0:3), AUX(0:3,7)
      REAL*8      A3(0:3), A4(0:3)
     
C
      N = 7
      DO I = 0,3
         AUX(I,7) = A7(I)
         AUX(I,6) = A6(I)
         AUX(I,5) = A5(I)
         AUX(I,4) = A4(I)
         AUX(I,3) = A3(I)
         AUX(I,2) = A2(I)
         AUX(I,1) = A1(I)
      ENDDO
      CHIAUX(1) = CHII(1)
      CHIAUX(2) = CHII(2)
      ALP = ALPHA
C
      DO 30 I = 1,N
         IF (ALP.GT.0) THEN
            ASLASH(1,1) = AUX(0,I) - AUX(3,I)
            ASLASH(1,2) = -AUX(1,I) + ZI*AUX(2,I)
            ASLASH(2,1) = -AUX(1,I) - ZI*AUX(2,I)
            ASLASH(2,2) = AUX(0,I) + AUX(3,I)
         ELSE
            ASLASH(1,1) = AUX(0,I) + AUX(3,I)
            ASLASH(1,2) = AUX(1,I) - ZI*AUX(2,I)
            ASLASH(2,1) = AUX(1,I) + ZI*AUX(2,I)
            ASLASH(2,2) = AUX(0,I) - AUX(3,I)
         ENDIF
         CHIDUM = CHIAUX(1)*ASLASH(1,1) + CHIAUX(2)*ASLASH(2,1)
         CHIAUX(2) = CHIAUX(1)*ASLASH(1,2) + CHIAUX(2)*ASLASH(2,2)
         CHIAUX(1) = CHIDUM
         ALP = - ALP
 30   CONTINUE
C
      SC7ccrrccc  = CHIAUX(1) * CHIF(1) + CHIAUX(2) * CHIF(2)
      RETURN
      END



CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      
      FUNCTION SC7ccrrcrc(CHII,A1,A2,A3,A4,A5,A6,A7,CHIF,ALPHA)
      IMPLICIT NONE 
      COMPLEX*16 ZI, SC7ccrrcrc
      PARAMETER ( ZI=(0D0,1D0) )
      INTEGER  ALPHA, ALP, N, I
      COMPLEX*16  CHII(2), CHIF(2), CHIAUX(2), CHIDUM, ASLASH(2,2)
      COMPLEX*16  A1(0:3), A2(0:3), A5(0:3), A7(0:3), AUX(0:3,7)
      REAL*8      A3(0:3), A4(0:3), A6(0:3)
     
C
      N = 7
      DO I = 0,3
         AUX(I,7) = A7(I)
         AUX(I,6) = A6(I)
         AUX(I,5) = A5(I)
         AUX(I,4) = A4(I)
         AUX(I,3) = A3(I)
         AUX(I,2) = A2(I)
         AUX(I,1) = A1(I)
      ENDDO
      CHIAUX(1) = CHII(1)
      CHIAUX(2) = CHII(2)
      ALP = ALPHA
C
      DO 30 I = 1,N
         IF (ALP.GT.0) THEN
            ASLASH(1,1) = AUX(0,I) - AUX(3,I)
            ASLASH(1,2) = -AUX(1,I) + ZI*AUX(2,I)
            ASLASH(2,1) = -AUX(1,I) - ZI*AUX(2,I)
            ASLASH(2,2) = AUX(0,I) + AUX(3,I)
         ELSE
            ASLASH(1,1) = AUX(0,I) + AUX(3,I)
            ASLASH(1,2) = AUX(1,I) - ZI*AUX(2,I)
            ASLASH(2,1) = AUX(1,I) + ZI*AUX(2,I)
            ASLASH(2,2) = AUX(0,I) - AUX(3,I)
         ENDIF
         CHIDUM = CHIAUX(1)*ASLASH(1,1) + CHIAUX(2)*ASLASH(2,1)
         CHIAUX(2) = CHIAUX(1)*ASLASH(1,2) + CHIAUX(2)*ASLASH(2,2)
         CHIAUX(1) = CHIDUM
         ALP = - ALP
 30   CONTINUE
C
      SC7ccrrcrc  = CHIAUX(1) * CHIF(1) + CHIAUX(2) * CHIF(2)
      RETURN
      END











CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

      SUBROUTINE PSI0M( NF,PBAR,SIGN,PSI )
      IMPLICIT NONE
      INTEGER  IF, NF, SIGN(NF)
      DOUBLE COMPLEX  PSI(2,-1:1,NF)
      DOUBLE PRECISION  PBAR(0:3,NF), PABS, PBARX, PBARY, PBARZ, PAPZ
      DOUBLE PRECISION  NORMAL, NX, NY, NZ, EPS
      PARAMETER ( EPS=1.d-30 )
CC
      DO IF = 1,NF
        PABS  = PBAR(0,IF)
        PBARX = PBAR(1,IF)
        IF ( MOD(IF,2).EQ.1 ) THEN
          PBARY =  PBAR(2,IF)
        ELSE
          PBARY = -PBAR(2,IF)
        ENDIF
        PBARZ = PBAR(3,IF)
        IF ( PBARZ.GT.0.d0 ) THEN
           PAPZ = PABS + PBARZ
        ELSE
           PAPZ = (PBARX**2 + PBARY**2)/(PABS - PBARZ)
        ENDIF
c
c treat the pbar along negative z-axis case first
c
        IF ( PAPZ.LE.EPS*PABS ) THEN
           NORMAL = SIGN(IF)*SQRT(2.d0*PBAR(0,IF))
c          Print*, 'Normal',PBAR(0,IF)
c          Print*, 'Normal',SQRT(2.d0*PBAR(0,IF))
c           Print*, 'Normal',Normal
           PSI(1,-1,IF) = -NORMAL 
           PSI(2,-1,IF) = 0.d0
           PSI(1,1,IF)  = 0.d0
           PSI(2,1,IF)  = NORMAL
        ELSE
c
c and now the general case
c
           NORMAL = SIGN(IF)/SQRT(PAPZ)
           NX = NORMAL*PBARX
           NY = NORMAL*PBARY
           NZ = NORMAL*PAPZ
           PSI(1,-1,IF) = DCMPLX(-NX,NY)
           PSI(2,-1,IF) = NZ
           PSI(1,1,IF)  = NZ
           PSI(2,1,IF)  = DCMPLX(NX,NY)
        ENDIF
      ENDDO
ccc
      END



      SUBROUTINE PSI0M_MAS(P,SIGN,PSI,m,hel,bar)
      IMPLICIT NONE
      INTEGER  SIGN,hel
      DOUBLE COMPLEX  PSI(4),ZI,u1(2),u2(2),EPHI,invEPHI,xp1,xp2
      DOUBLE PRECISION  p(0:3),m,EPS,norm1,norm2,cth2,sth2
      DOUBLE PRECISION pt2,pz2,TH,pi
      PARAMETER ( EPS=1.d-16,ZI=(0d0,1d0),pi=3.141592653589793d0)
      LOGICAL bar

      pt2=(p(1)*p(1)+ p(2)*p(2))
      pz2=p(3)*p(3)
c      print*, 'Sqrt(pt2/pz2)',Sqrt(pt2/pz2)
      IF (Sqrt(pt2/pz2).LT.eps) THEN
c         print*, 'here0'
         IF(P(3).gE.0D0) THEN 
c         print*, 'herem01'
          TH=0D0
          cth2=1d0
          sth2=0d0
         ELSE
c         print*, 'herem02'
          TH=pi
          cth2=0d0
          sth2=-1d0
         endif
          EPHI=1D0
          norm1=sqrt(abs(p(0))+abs(p(3)))
          norm2=m/norm1
c          print*, 'm',m
c        print*,'norm1*norm2',norm1*norm2
c         print*,'norm1',norm1
c         print*,'norm2',norm2
c          print*,'norm1*norm2',norm1*norm2
      ELSE
c          print*, 'here01'
          norm1=sqrt(p(0)+sqrt(pt2+pz2))
          norm2=m/norm1
c          print*,'norm1*norm2',norm1*norm2
	  if (p(3).gt.0d0) then
          th=atan(sqrt(pt2)/p(3))
	  else
	  th=pi-atan(sqrt(pt2)/(-p(3))) 
	  endif
c	  print*, 'th',th
          ephi=(p(1)+ ZI*p(2))/sqrt(pt2)
          cth2=cos(th/2d0)
c	  print*,'cth2',cth2
          sth2=sin(th/2d0)
      ENDIF
      If (Hel.eq.1) then
          if (sign.eq.1) then
c           print*, 'here11'
            xp1=cth2
            xp2=ephi*sth2
            u1(1)=norm2*xp1
            u1(2)=norm2*xp2
            u2(1)=norm1*xp1
            u2(2)=norm1*xp2
         else
c           print*, 'here12'
            invephi=1d0/ephi
            xp1=-invephi*sth2
            xp2=cth2
            u1(1)=-norm1*xp1
            u1(2)=-norm1*xp2
            u2(1)=norm2*xp1
            u2(2)=norm2*xp2            
         endif
      else
         IF (SIGN.EQ.1) then
c          print*, 'here21'
c          print*, 'ephi',ephi

            invephi=1d0/ephi
            xp1=-invephi*sth2
            xp2=cth2
c            print*, 'xp1',xp1
c           print*, 'xp2',xp2
c         print*,'norm1*norm2',norm1*norm2
c         print*,'norm1',norm1
c         print*,'norm2',norm2
            u1(1)=norm1*xp1
            u1(2)=norm1*xp2
            u2(1)=norm2*xp1
            u2(2)=norm2*xp2
         else
c          print*, 'here22'
            xp1=cth2
            xp2=ephi*sth2
            u1(1)=norm2*xp1
            u1(2)=norm2*xp2
            u2(1)=-norm1*xp1
            u2(2)=-norm1*xp2
         endif   
      endif   
      if (bar) then
c          print*, 'here31'

      PSI(1)=Conjg(u2(1))
      PSI(2)=Conjg(u2(2))
      PSI(3)=Conjg(u1(1))
      PSI(4)=Conjg(u1(2))
      else
c          print*, 'here32'
      PSI(1)=u1(1)
      PSI(2)=u1(2)
      PSI(3)=u2(1)
      PSI(4)=u2(2)
      endif
      return
      end
C************************************  CURR  *************************
C
C  CURR calculates the current 
C
C        PSIBAR gamm^mu PSI
C 
C  of Hagiwara, Zeppenfeld, Nucl.Phys.B313 (1989) 560, eq. 2.20 for 
C  the two possible polarizations of the Weyl spinors PSIBAR and PSI
C
C  INPUT:
C  ------
C
C    PSIBAR(2,-1:1)  two Weyl spinors for 2 helicity states each
C    PSI(2,-1:1)
C
C  OUTPUT:
C  -------
C
C  J(0:3,-1:1)       the current as defined in eq.2.27 for two possible
C                    helicity combinations for massless fermions
C
      subroutine curr( sigmax,psibar,psi,j )
      double complex  j(0:3,-1:1), psibar(2,-1:1), psi(2,-1:1), zj2
      double complex  z1, z2, z3, z4
      integer  sigmax,sig
cc
      do sig = -1,sigmax,2
         z1 = psibar(1,sig) * psi(1,sig)
         z2 = psibar(2,sig) * psi(2,sig)
         z3 = psibar(1,sig) * psi(2,sig)
         z4 = psibar(2,sig) * psi(1,sig)
         j(0,sig) = z1 + z2
         if (sig.eq.-1) then
            j(1,sig) = -(z3+z4)
            zj2 = z3-z4
            j(2,sig) = dcmplx(-dimag(zj2),dreal(zj2))
            j(3,sig) = z2-z1
         else
            j(1,sig) = z3+z4
            zj2 = z4-z3
            j(2,sig) = dcmplx(-dimag(zj2),dreal(zj2))
            j(3,sig) = z1-z2
         endif
      enddo
ccc
      end
C


