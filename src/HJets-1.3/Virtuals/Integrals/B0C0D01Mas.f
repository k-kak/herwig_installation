c      complex*16 function A0fin(M,musq)
c      complex*16 function A0finDiv(M,i)
c      complex*16 function B0t1(M1,s,musq)
c      complex*16 function B0tMDiv(M1,s,i)
c      complex*16 function C0fin1M(m,q1sq,q2sq,Psq)
c      complex*16 function I3point(m,q1sq,q2sq,Psq)
c      complex*16 function C0finG2(M1,M2,M3,s1,s2,s3,musq)
c        FUNCTION D04(P1,P2,P3,P4,P12,P23,M1,M2,M3,M4)
c        FUNCTION ETA(C1,C2)                                            
c        FUNCTION ETAS(Y,R,RS)                                            
c        FUNCTION SQE(A,B,C)                                            
c      function E01M(m,p1,p2,p3,p4,p1p2,p1p3,p1p4,p2p3,p2p4,p3p4,
c      complex*16 function B0tM1(m,qsq) 
c       complex*16 function cdilog(z)    
c






c --------------------------------------
c   A0finG2=A0finG(M,musq)                                 
c --------------------------------------
      complex*16 function A0fin(M,musq)
      implicit none
      real*8 pi
      complex*16 Ipi, Ieps 
      parameter(pi=3.14159265358979324d0,
     & Ipi=(0d0,3.14159265358979324d0),Ieps=(0d0,1d-38))

      double precision M,musq
      A0fin = M*M*(-dLog(M*M/musq)+1.d0)
      end

c --------------------------------------
c   A0finG2=A0finG(M,musq)                                 
c --------------------------------------
      complex*16 function A0finDiv(M,i)
      implicit none
      double precision M,musq
      integer i
      if(i.eq.1) then
      A0finDiv = M*M
      else
      A0finDiv=0d0
      endif
      end


c --------------------------------------
c   B0finG3=B0finG(M1,M1,s,musq)                                 
c --------------------------------------
      complex*16 function B0t1(M1,s,musq)
      implicit none
      real*8 pi
      complex*16 Ipi, Ieps 
      parameter(pi=3.14159265358979324d0,
     & Ipi=(0d0,3.14159265358979324d0),Ieps=(0d0,1d-38))
      double precision M1,s,musq,msq
      complex*16 b,g1,g2

      msq = M1*M1
      if(abs(s).gt.1d-6) then
      b = Sqrt(s*s-4.d0*s*(msq-Ieps))
      g1 = .5d0*(s+b)/s
      g2 = .5d0*(s-b)/s

      B0t1 = -Log((s-Ieps)/musq) 
     &        + (g1*Log((g1-1.d0)/g1)-Log(g1-1))
     &        + (g2*Log((g2-1.d0)/g2)-Log(g2-1)) + 2.d0 


      else
      B0t1  = -dLog(msq/musq)
      
      endif


      end



c --------------------------------------
c   B0finG3=B0finG(M1,M1,s,musq)                                 
c --------------------------------------
      complex*16 function B0tMDiv(M1,s,i)
       double precision M1,s,musq,msq
      complex*16 b,g1,g2
      integer i
      if(i.eq.1) then
      B0tMDiv = 1.d0
      else
      B0tMDiv = 0d0
      endif
      end




c------------- C(m,q1^2,q2^2,P^2)  3-point fuction ---------------------------
c
      complex*16 function C0fin1M(m,q1sq,q2sq,Psq)
      implicit none
      double precision m, q1sq, q2sq, Psq,eps 
      double complex I3point
      External I3point
      parameter(eps=1d-07 )

c Author: Francisco Campanario
c Date: 15 08 2008
c(,,+)
      If (Psq.gt.eps) then
c(-,,+)
      If (q1sq.lt.-eps) then
      C0fin1M=I3point(m,Psq,q2sq,q1sq) 
      Return   
c(,-,+)
      else if(q2sq.lt.-eps) then
      C0fin1M=I3point(m,q1sq,Psq,q2sq)
      Return  
      else 
c(+,+,+)
c      Print*,3
C      C0fin1M=-Dble(I3point(m,-q1sq,-Psq,-q2sq))+(0,1)*DIMAG(I3point(m,-q1sq,-Psq,-q2sq))   
      C0fin1M=I3point(m,q1sq,q2sq,Psq)
      Return   
      endif
c(,,-)
      elseif(Psq.lt.-eps) then
      C0fin1M=I3point(m,q1sq,q2sq,Psq)
c      Print*,4
      Return  
      else
c(-,,0)
      If (q1sq.lt.-eps) then
      C0fin1M=I3point(m,Psq,q2sq,q1sq)   
      Return  
c(,-,0)
      else if(q2sq.lt.-eps) then
      C0fin1M=I3point(m,q1sq,Psq,q2sq)
      Return  
      else 
c(+,+,0)
c      C0fin1M=-Dble(I3point(m,-q1sq,-Psq,-q2sq))+(0,1)*DIMAG(I3point(m,-q1sq,-Psq,-q2sq)) 
      C0fin1M=I3point(m,q1sq,q2sq,Psq)
      Return  
      endif
      endif
      end


c
c------------- C(m,q1^2,q2^2,P^2)  3-point fuction ---------------------------
c
      complex*16 function I3point(m,q1sq,q2sq,Psq)
      implicit none
      double precision m, q1sq, q2sq, psq 

c evaluate scalar 3-point function for equal masses m on propagators 
c  
c  I3 = 1/(i*pi^2) * Int d^4k [k^2-m^2]^-1 [(k+q1)^2-m^2]^-1 [(k-q2)^2-m^2]^-1
c
c     = i/pi^2 C('t Hooft,Veltman) = -C0(Passarino,Veltman)
c
c in terms of 12 dilogarithms. P = q1+q2 is assumed to be space-like. q1 
c and q2 may be space-, time- or lightlike
c
c	Dieter Zeppenfeld, <dieter@pheno.physics.wisc.edu>
c	Initial version:  2000 November 12
c	Last modified:    
c  

      double precision lambda, rtlam, q1q2, q1p, q2p, fmsq
      double precision rp(3), rm(3), eps, precise, pi

      complex*16 sp(3), sm(3), fmt2, ieps, fac, z, di, li2l, li2
      external li2, li2l
      parameter (eps=1d-18, ieps = eps*(0,1), precise=1d-3 )
      parameter (pi=3.14159 26535 89793d0, fac=-1d0/2d0)

c determine dot products of input momenta
      q1q2 = 0.5d0*(psq - q1sq - q2sq)
      q1p  = 0.5d0*(psq + q1sq - q2sq)
      q2p  = 0.5d0*(psq - q1sq + q2sq)
      lambda = q1q2**2-q1sq*q2sq
      
      
      if (lambda.le.0d0) then
c        write(*,*) ' singular lambda in 3-point function. Reset to 0 '
         I3point = 0
         return
      else
         rtlam = sqrt(lambda)
      endif
      
      
      
c determine factors for call of spence functions
      fmsq = 4*m**2
      fmt2 = fmsq - ieps
      di = q1sq*q2sq*psq + fmt2*lambda 
      if ( abs(dreal(di)).gt.eps )then
         di = 1d0/di
      else
c         print*,' singular point called in I3point '
         I3point = 0
         return
c         stop
      endif

      rp(1) = q2p + rtlam
      rm(1) = q2p - rtlam
      if (abs(rm(1)).lt.abs(rp(1))*precise ) then
         rm(1) = q2sq*psq/rp(1)
      elseif (abs(rp(1)).lt.abs(rm(1))*precise ) then
         rp(1) = q2sq*psq/rm(1)
      endif

      rp(2) = q1q2 + rtlam
      rm(2) = q1q2 - rtlam
      if (abs(rm(2)).lt.abs(rp(2))*precise ) then
         rm(2) = q1sq*q2sq/rp(2)
      elseif (abs(rp(2)).lt.abs(rm(2))*precise ) then
         rp(2) = q1sq*q2sq/rm(2)
      endif

      rp(3) = q1p + rtlam
      rm(3) = q1p - rtlam
      if (abs(rm(3)).lt.abs(rp(3))*precise ) then
         rm(3) = q1sq*psq/rp(3)
      elseif (abs(rp(3)).lt.abs(rm(3))*precise ) then
         rp(3) = q1sq*psq/rm(3)
      endif

      if (q1sq.ne.0d0) then
         z = rtlam*sqrt((q1sq-fmt2)*q1sq) 
         sp(1) = di*( q2p*q1sq + z )
         sm(1) = di*( q2p*q1sq - z )
      else
         sp(1) = 0
         sm(1) = 0
      endif

      if (psq.ne.0d0) then
         z = rtlam*sqrt((psq-fmt2)*psq) 
         sp(2) = di*( q1q2*psq + z )
         sm(2) = di*( q1q2*psq - z )
      else
         sp(2) = 0
         sm(2) = 0
      endif

      if (q2sq.ne.0d0) then
         z = rtlam*sqrt((q2sq-fmt2)*q2sq) 
         sp(3) = di*( q1p*q2sq + z )
         sm(3) = di*( q1p*q2sq - z )
      else
         sp(3) = 0
         sm(3) = 0
      endif

      z = 0
      if (q1sq.lt.0d0 .or. q1sq.gt.fmsq) then
         z = z + li2(rm(1)*sp(1)) + li2(rm(1)*sm(1)) - 
     &           li2(rp(1)*sm(1)) - li2(rp(1)*sp(1))
      elseif (q1sq.gt.0d0) then
         z = z + 2*dreal( li2(rm(1)*sp(1)) ) - 
     &           2*dreal( li2(rp(1)*sp(1)) )
      endif
      if (psq.lt.0d0 .or. psq.gt.fmsq) then
         z = z - li2(rm(2)*sp(2)) - li2(rm(2)*sm(2)) + 
     &           li2(rp(2)*sm(2)) + li2(rp(2)*sp(2))
      elseif (psq.gt.0d0) then
         z = z - 2*dreal( li2(rm(2)*sp(2)) ) + 
     &           2*dreal( li2(rp(2)*sp(2)) )
      endif
      if (q2sq.lt.0d0 .or. q2sq.gt.fmsq) then
         z = z + li2(rm(3)*sp(3)) + li2(rm(3)*sm(3)) - 
     &           li2(rp(3)*sm(3)) - li2(rp(3)*sp(3))
      elseif (q2sq.gt.0d0) then
         z = z + 2*dreal( li2(rm(3)*sp(3)) ) - 
     &           2*dreal( li2(rp(3)*sp(3)) )
      endif

      I3point = z*fac/rtlam
      return
      end
c







c --------------------------------------
c   C0finG2=C0finG full finite
c   only s3 can be zero
c --------------------------------------
      complex*16 function C0finG2(M1,M2,M3,s1,s2,s3,musq)
      implicit none
      real*8 pi
      complex*16 Ipi, Ieps 
      parameter(pi=3.14159265358979324d0,
     & Ipi=(0d0,3.14159265358979324d0),Ieps=(0d0,1d-38))
      double precision M1,M2,M3,s1,s2,s3,musq,m1sq,m2sq,m3sq
      complex*16 x(3), ys(3,2), a, b, c, d, e, f,
     &                 alfa, deno, sqr,ce,ad,abc,de
      complex*16  cdilog, x1, x2, C0
      external    cdilog
      integer i,j
      REAL*8 eps
      parameter (eps=1d-9)
     
      m1sq = M1*M1
      m2sq = M2*M2
      m3sq = M3*M3
c     ---
      a = s3
      b = s2
      c = s1-s2-s3
      d = m3sq-m1sq-s3
      e = m2sq-m3sq+s3-s1
      f = m1sq-Ieps
      ce=m2sq-m3sq-s2
      ad=m3sq-m1sq
      abc=s1
      de=m2sq-m1sq-s1      
c     ---
      alfa = (-c + Sqrt(c*c-4.d0*b*a-2.d0*Ieps*(s1+s2+s3)))/(2.d0*b)
      deno = c + 2.d0*alfa*b

      x(1) = - (d+2.d0*a+(ce)*alfa)/deno
      x(2) = - (d+e*alfa)/((1.d0-alfa)*deno)
      x(3) = (d+e*alfa)/(alfa*deno)

      sqr = Sqrt((ce)*(ce)-4.d0*b*(ad+f))
      ys(1,1) = (-(ce)+sqr)/(2.d0*b)
      ys(1,2) = (-(ce)-sqr)/(2.d0*b)


  
      sqr = Sqrt((de)*(de)-4.d0*f*(abc))
      ys(2,1) = (-(de)+sqr)/(2.d0*(abc))
      ys(2,2) = (-(de)-sqr)/(2.d0*(abc))
  

      
      C0 = (0.d0,0.d0)

      if (abs((a*a)).gt.eps) then 
 

      sqr = Sqrt(d*d-4.d0*a*f)
      ys(3,1) = (-d+sqr)/(2.d0*a)
      ys(3,2) = (-d-sqr)/(2.d0*a)
      do i = 1,3
         do j = 1,2
            x1 = x(i)/(x(i)-ys(i,j))
            x2 = (x(i)-1.d0)/(x(i)-ys(i,j))
            C0 = C0 + ((-1.d0)**i)*(cdilog(x1) -  cdilog(x2))
         enddo
      enddo 
      else
      do i = 1,2
         do j = 1,2
           x1 = x(i)/(x(i)-ys(i,j))
            x2 = (x(i)-1.d0)/(x(i)-ys(i,j))
            C0 = C0 + ((-1.d0)**i)*(cdilog(x1) -  cdilog(x2))

         enddo
      enddo 
      endif
      C0finG2 = C0/deno
      end














************************************************************************
c        FUNCTION D04(P1,P2,P3,P4,P12,P23,M1,M2,M3,M4)
        FUNCTION D04(P1t,P2t,P3t,P4t,P12t,P23t,M1t,M2t,M3t,M4t)
************************************************************************
*  SCALAR 4-POINT FUNCTION WITH AT LEAST ONE MASS ZERO                 *
*  P1,P2,P3,P4 = SQUARED EXTERNAL MOMENTA			       *
*  P12 = (p1+p2)**2,  P23 = (p2+p3)**2				       *
*----------------------------------------------------------------------*
*  2.1.92  SD	         					       *
*  rearrangement to exploit massless external momenta   14.3.01  DZ    *
*  Modified: Michael Kubocz                                            *
*  Interception of NANs e.g. caused by log(0) etc. (see below)         *
************************************************************************
        IMPLICIT REAL*8 (A-Z)
	REAL*8 M(4),P(4,4),K(4,4)
        real*8 pi,eps,eps1
        real*8 im1,im2
        real*8 m1,m2,m3,m4
        real*8 m1t,m2t,m3t,m4t
        real*8 m02,m12,m22,m32,m42
        real*8 mm0,mm1,mm2,mm3,mm4 
        real*8 p1,p2,p3,p4,p12,p23
        real*8 p1t,p2t,p3t,p4t,p12t,p23t
        real*8 q0,q1,q2,q3,q4,q00,q12,q23
	COMPLEX*16 A1,A2,A3,A4,SWAP
	COMPLEX*16 SS(4), XX(2), X(2,4),RS(4,4)
	COMPLEX*16 S0(4),XX0(2),X0(2,4), R(4,4),G(2)
        COMPLEX*16 D04,LI2,ETA,SQE,ETAS
	COMPLEX*16 AA,BB,CC,DD,IEPS,H,HH,L1,L2,L3,L4
	INTEGER I,J, i1,i2,i3,i4
        complex*16 D1,D2,D3,D4

        EXTERNAL LI2,ETA,ETAS,SQE
 
        D1=DCMPLX(0d0,0d0)
        D2=DCMPLX(0d0,0d0)
        D3=DCMPLX(0d0,0d0)
        D4=DCMPLX(0d0,0d0)
        eps1=1d-7

        if(abs(P1t).le.eps1) then
           P1=0d0
        else
           P1=P1t
        endif
        if(abs(P2t).le.eps1) then
           P2=0d0
        else
           P2=P2t
        endif
        if(abs(P3t).le.eps1) then
           P3=0d0
        else
           P3=P3t
        endif
        if(abs(P4t).le.eps1) then
           P4=0d0
        else
           P4=P4t
        endif
        if(abs(P12t).le.eps1) then
           P12=0d0
        else
           P12=P12t
        endif
        if(abs(P23t).le.eps1) then
           P23=0d0
        else
           P23=P23t
        endif
        if(abs(M1t).le.eps1) then
           M1=0d0
        else
           M1=M1t
        endif
        if(abs(M2t).le.eps1) then
           M2=0d0
        else
           M2=M2t
        endif
        if(abs(M3t).le.eps1) then
           M3=0d0
        else
           M3=M3t
        endif
        if(abs(M4t).le.eps1) then
           M4=0d0
        else
           M4=M4t
        endif
c FC %        print*,"D04arg ",p1t,p2t,p3t,p4t,p12t,p23t,m1t,m2t,m3t,m4t
c FC %        print*,"D04arg ",p1,p2,p3,p4,p12,p23,m1,m2,m3,m4
        MM1=M1
        MM2=M2
        MM3=M3
        MM4=M4
        M12=M1*M1
        M22=M2*M2
        M32=M3*M3
        M42=M4*M4
        Q1=P1
        Q2=P2
        Q3=P3
	Q4=P4
        Q12=P12
        Q23=P23
c FC %       print*, "Equality", MM1*MM2*MM3*MM4
c FC %        print*, "Equality", MM1*MM2*MM3*MM4.NE.0D0

C	IS AT LEAST ONE MASS ZERO ???
	IF (MM1*MM2*MM3*MM4.NE.0D0) GOTO 130

C	PERMUTATE UNTIL MM3=0D0
	GOTO 20
10	CONTINUE
	MM0=MM1
	MM1=MM2
	MM2=MM3
	MM3=MM4
	MM4=MM0
	M02=M12
	M12=M22
	M22=M32
	M32=M42
	M42=M02
	Q00=Q12
	Q12=Q23
	Q23=Q00
	Q0=Q1
	Q1=Q2
	Q2=Q3
	Q3=Q4
	Q4=Q0
20	IF (MM3.NE.0D0) GOTO 10
C	ONLY MM3 IS ZERO
	IF (MM1*MM2*MM4.NE.0D0) GOTO 30
C	ONLY MM3 AND MM4 ARE ZERO ==> 3->2, 4->3...
	IF ((MM1*MM2.NE.0D0).AND.(MM4.EQ.0D0)) GOTO 10
C	ONLY MM2 AND MM3 ARE ZERO
	IF ((MM1*MM4.NE.0D0).AND.(MM2.EQ.0D0)) GOTO 40
	WRITE(*,*)'CASE OF THIS SPECIAL D0-FUNCTION NOT IMPLEMENTED!'
	STOP

C	****** NO MASS EQUAL TO ZERO ******
130	CONTINUE
	EPS=1D-18
	IEPS=DCMPLX(0D0,EPS)
c check for massless external momentum: excellent candidate for p13,
c leading to r13 >=1 and real.
c$$$        if (q1*q2*q3*q4.eq.0d0) then !org
        if (abs(q1*q2*q3*q4).le.eps1) then
c$$$           if (q2.eq.0d0) then
           if (abs(q2).le.eps1) then
              I1 = 2
              I2 = 3
              I3 = 1
              I4 = 4
c$$$           elseif (q1.eq.0d0) then
           elseif (abs(q1).le.eps1) then
              I1 = 1
              I2 = 3
              I3 = 2
              I4 = 4
c$$$           elseif (q3.eq.0d0) then
           elseif (abs(q3).le.eps1) then
              I1 = 2
              I2 = 4
              I3 = 1
              I4 = 3
           else
              I1 = 1
              I2 = 4
              I3 = 2
              I4 = 3
           endif
           M(i1)=MM1
           M(i2)=MM2
           M(i3)=MM3
           M(i4)=MM4
           P(i1,i2)=Q1
           P(i3,i2)=Q2
           P(i3,i4)=Q3
           P(i1,i4)=Q4
           P(i1,i3)=Q12
           P(i3,i1)=Q12
           P(i2,i4)=Q23
           P(i4,i2)=Q23
	ELSEIF( ABS((MM1**2+MM3**2-Q12)/MM1/MM3).LT.2D0 ) THEN
C	R13 WOULD BE NOT REAL. -> PERMUTATION! -> R(2,4) IS NOT REAL.
	   M(1)=MM2
	   M(2)=MM3
	   M(3)=MM4
	   M(4)=MM1
	   P(1,2)=Q2
	   P(1,3)=Q23
	   P(1,4)=Q1
	   P(2,3)=Q3
	   P(2,4)=Q12
	   P(3,4)=Q4
	ELSE
C	R(1,3) IS REAL.
	   M(1)=MM1
	   M(2)=MM2
	   M(3)=MM3
	   M(4)=MM4
	   P(1,2)=Q1
	   P(1,3)=Q12
	   P(1,4)=Q4
	   P(2,3)=Q2
	   P(2,4)=Q23
	   P(3,4)=Q3
	ENDIF

	DO 11 J=2,4
	DO 11 I=1,J-1
	K(I,J)=(M(I)**2+M(J)**2-P(I,J))/M(I)/M(J)
	R(I,J) =SQE(DCMPLX(1D0,0D0),DCMPLX(-K(I,J),0D0),
     *	            DCMPLX(1D0,0D0))
        IF( DBLE(K(I,J)).LT.-2D0 ) THEN
c        IF( DIMAG(R(I,J)).EQ.0D0 ) THEN
	   RS(I,J)=SQE(DCMPLX(1D0,0D0),DCMPLX(-K(I,J),EPS),
     *	               DCMPLX(1D0,0D0))
	ELSE
	   RS(I,J)=R(I,J)
	ENDIF
11	CONTINUE

	SS(1)=RS(1,2)
	SS(2)=RS(2,3)
	SS(3)=RS(3,4)
	SS(4)=RS(1,4)
	S0(1)=R(1,2)
	S0(2)=R(2,3)
	S0(3)=R(3,4)
	S0(4)=R(1,4)
	AA=K(3,4)/R(2,4)+R(1,3)*K(1,2)-K(1,4)*R(1,3)/R(2,4)-K(2,3)
	BB=(R(2,4)-1D0/R(2,4))*(R(1,3)-1D0/R(1,3))
     *		+K(1,2)*K(3,4)-K(1,4)*K(2,3)
	CC=K(1,2)/R(1,3)+R(2,4)*K(3,4)-K(1,4)*R(2,4)/R(1,3)-K(2,3)
	DD=K(2,3)-R(1,3)*K(1,2)-R(2,4)*K(3,4)+R(1,3)*R(2,4)*K(1,4)
	XX(1)=SQE(AA,BB,CC+IEPS*DD)
	XX(2)=(CC+IEPS*DD)/AA/XX(1)
	XX0(1)=SQE(AA,BB,CC)
	XX0(2)=CC/AA/XX0(1)
c	IF (ABS(DREAL(XX0(1)-XX(2))).LT.ABS(DREAL(XX0(1)-XX(1)))) THEN
	IF (ABS(XX0(1)-XX(2)).LT.ABS(XX0(1)-XX(1))) THEN
	  SWAP  =XX0(1)
	  XX0(1)=XX0(2)
	  XX0(2)=SWAP
	ENDIF

	DO 12 I=1,2
	G(I)  =SIGN( 1D0,DREAL(AA*(XX(I)-XX(3-I))) )
	 X(I,1)= XX(I)/R(2,4)
	X0(I,1)=XX0(I)/R(2,4)
	 X(I,2)= XX(I)/R(2,4)*R(1,3)
	X0(I,2)=XX0(I)/R(2,4)*R(1,3)
	 X(I,3)= XX(I)*R(1,3)
	X0(I,3)=XX0(I)*R(1,3)
	 X(I,4)= XX(I)
	X0(I,4)=XX0(I)
12	CONTINUE

	D04 = DCMPLX(0D0,0D0)
	DO 13 I=1,2
	DO 13 J=1,4
	A1 = 1D0+X0(I,J)*S0(J) + ABS(1D0+X0(I,J)*S0(J))*IEPS*
     *				  SIGN(1D0,DIMAG(X(I,J)*SS(J)))
	A2 = 1D0+X0(I,J)/S0(J) + ABS(1D0+X0(I,J)/S0(J))*IEPS*
     *				  SIGN(1D0,DIMAG(X(I,J)/SS(J)))
c  org. code:
c$$$           D04 = D04 + (-1D0)**(I+J)*(
c$$$     &          LI2(A1)+ETA(-X(I,J),SS(J))*LOG(A1)
c$$$     &          +LI2(A2)+ETA(-X(I,J),1D0/SS(J))*LOG(A2))

        if(abs(ETA(-X(I,J),SS(J))).ne.0d0) then
           D1=ETA(-X(I,J),SS(J))*LOG(A1)
        else
           D1=DCMPLX(0d0,0d0)
        endif
        if(abs(ETA(-X(I,J),1D0/SS(J))).ne.0d0) then
           D2=ETA(-X(I,J),1D0/SS(J))*LOG(A2)
        else
           D2=DCMPLX(0d0,0d0)
        endif
        D04=D04+(-1D0)**(I+J)*(LI2(A1)+LI2(A2)+D1+D2)

c FC %        Print*, 'abs(ETA(-X(I,J),SS(J)))',abs(ETA(-X(I,J),SS(J)))
c FC %        Print*,  'A1',A1
c FC %        Print*, 'abs(ETA(-X(I,J),1D0/SS(J)))',abs(ETA(-X(I,J),1D0/SS(J)))
c FC %        print*,  'a2',a2
c FC %        Print*, 'D04',D04
c   The enquiry avoids occurrence of NANs causing by LOG(A1) for A1=0 
c   and LOG(A2) for A2=0. At that points also ETA(-X(I,J),1D0/SS(J)) 
c   or ETA(-X(I,J),SS(J) are 0. (Michael Kubocz)

13	CONTINUE

c        print*,'DIMAG(R(1,3))',DIMAG(R(1,3))
c	IF( DIMAG(R(1,3)).EQ.0D0 ) THEN !org (makes troubles in squark pentagons)
	IF( abs(DIMAG(R(1,3))).le.eps1 ) THEN
	DO 14 I=1,2
	   A1 = (K(1,3)-2D0*R(1,3))/XX0(I)
     *		      -R(1,3)*K(1,4)+K(3,4)
     	   A2 = ((K(2,4)-2D0*R(2,4))*R(1,3)*XX0(I)
     *		      -R(2,4)*K(3,4)+K(2,3))/DD
	   A3 = (K(1,3)-2D0*R(1,3))*R(2,4)/XX0(I)
     *		      -R(1,3)*K(1,2)+K(2,3)
	   A4 = ((K(2,4)-2D0*R(2,4))*XX0(I)
     *		      -R(2,4)*K(1,4)+K(1,2))/DD
	   L1 = LOG( A1-ABS(A1)*IEPS )
     	   L2 = LOG( A2+ABS(A2)*IEPS*G(I)*SIGN(1D0,DREAL(R(1,3))
     *				        	  *DIMAG(RS(2,4))) )
	   L3 = LOG( A3-ABS(A3)*IEPS )
	   L4 = LOG( A4+ABS(A4)*IEPS*G(I)*SIGN(1D0,DIMAG(RS(2,4))) )

c org. code:
c$$$	   D04 = D04 
c$$$     &         + (3D0-2D0*I)*(
c$$$     *		       ETAS( -XX(I),R(1,3),RS(1,3) )
c$$$     *		          *( LOG(R(1,3)*XX(I)) + L1 + L2 )
c$$$     *		     + ETAS( -XX(I),1D0/R(2,4),1D0/RS(2,4) )
c$$$     *		          *( LOG(XX(I)/R(2,4)) + L3 + L4 )
c$$$     *		     - ( ETAS( -XX(I),R(1,3)/R(2,4),RS(1,3)/RS(2,4) )
c$$$     *		       + ETA( RS(1,3),1D0/RS(2,4) )                  )
c$$$     *		        *( LOG(XX(I)*R(1,3)/R(2,4)) + L3 + L2 )
c$$$     *	  	     + ETA( RS(1,3),1D0/RS(2,4) )
c$$$     *		       *ETAS(-XX(I),-R(1,3)/R(2,4),-RS(1,3)/RS(2,4))   )

            if(abs(ETAS(-XX(I),R(1,3),RS(1,3))).ne.0d0) then
               D1=ETAS(-XX(I),R(1,3),RS(1,3))*(LOG(R(1,3)*XX(I))+L1+L2)
            else
               D1=DCMPLX(0d0,0d0)
            endif
            if(abs(ETAS(-XX(I),1D0/R(2,4),1D0/RS(2,4))).ne.0d0) then
               D2=ETAS(-XX(I),1D0/R(2,4),1D0/RS(2,4))*(LOG(XX(I)/R(2,4))
     &           +L3+L4)
            else
               D2=DCMPLX(0d0,0d0)
            endif
            if((abs(ETAS(-XX(I),R(1,3)/R(2,4),RS(1,3)/RS(2,4))).ne.0d0)
     &          .or.(abs(ETA(RS(1,3),1D0/RS(2,4))).ne.0d0)) then 
               D3=-(ETAS(-XX(I),R(1,3)/R(2,4),RS(1,3)/RS(2,4))
     &              +ETA(RS(1,3),1D0/RS(2,4)))*(LOG(XX(I)*R(1,3)/R(2,4))
     $              +L3+L2)
            else
               D3=DCMPLX(0d0,0d0)
            endif
            D4=ETA(RS(1,3),1D0/RS(2,4))*ETAS(-XX(I),-R(1,3)/R(2,4),
     $           -RS(1,3)/RS(2,4))
            D04=D04+(3D0-2D0*I)*(D1+D2+D3+D4)
c   The enquiry avoids occurrence of NANs causing by LOG(0). At that points 
c   also ETA(...) and ETAS(...) are 0. (Michael Kubocz)

14	CONTINUE
	ELSE
	DO 15 I=1,2
	   L1 = LOG( R(2,4)/XX0(I)+XX0(I)/R(2,4)+K(1,2)
     *		     -XX0(I)/R(2,4)*EPS*BB*G(I) )
	   L2 = LOG( R(1,3)*XX0(I)+1D0/XX0(I)/R(1,3)+K(3,4)
     *		     -XX0(I)*R(1,3)*EPS*BB*G(I) )
	   L3 = LOG( R(1,3)/R(2,4)*XX0(I)+R(2,4)/XX0(I)/R(1,3)+K(2,3)
     *		     -XX0(I)*R(1,3)/R(2,4)*EPS*BB*G(I) )

	   D04 = D04 
     &          + (3D0-2D0*I)*(
     *		     ETA(-XX(I),1D0/R(2,4))
     *		      *( LOG(XX(I)/R(2,4)) + L1 )
     *		    +ETA(-XX(I),R(1,3))
     *		      *( LOG(R(1,3)*XX(I)) + L2 )
     *		    -( ETA(-XX(I),R(1,3)/R(2,4))
     *		      +ETA(R(1,3),1D0/R(2,4)) )
     *		      *( LOG(XX(I)*R(1,3)/R(2,4)) + L3 )
     *	  	    +ETA(R(1,3),1D0/R(2,4))
     *		      *ETA(-XX(I),-R(1,3)/R(2,4))
     *		      *(1D0-G(I)*SIGN(1D0,DREAL(BB)))	    )
15	CONTINUE
	ENDIF

	D04 = D04/M(1)/M(2)/M(3)/M(4)/AA/(XX(1)-XX(2))
	RETURN


C--->	***************** SPEZIELL ( --> T.SACK-PROMOTION )
C	D1=Q12-M12
C	D2=Q2 -M22
C	D3=Q3 -M42
C	IF ((D1*D2.LE.0D0).OR.(D2*D3.LE.0D0)) THEN
C	   WRITE(*,*) 'THE CASE OF DIFFERENT SIGNS OF THE D1,D2,D3'
C	   WRITE(*,*) 'IN D04(...) IS NOT IMPLEMENTED !!!'
C	   STOP
C	ENDIF
C	NM1=ABS(MM1/D1)
C	NM2=ABS(MM2/D2)
C	NM3=ABS(MM4/D3)
C	NP1=Q2/D2**2+Q12/D1**2+(Q1-Q2-Q12)/D1/D2
C	NP2=Q2/D2**2+ Q3/D3**2+(Q23-Q2-Q3)/D2/D3
C	NP3=Q3/D3**2+Q12/D1**2+(Q4-Q3-Q12)/D1/D3
C	D04=C04(NP1,NP2,NP3,NM1,NM2,NM3)/D1/D2/D3

C	*************** ALLGEMEIN


C	****** ONLY MM3 IS ZERO ******
30	CONTINUE
	EPS=1D-17
	IEPS=DCMPLX(0D0,EPS)
	M(1)=MM1
	M(2)=MM2
	M(3)=10D0
	M(4)=MM4
	P(1,2)=Q1
	P(1,3)=Q12
	P(1,4)=Q4
	P(2,3)=Q2
	P(2,4)=Q23
	P(3,4)=Q3
	DO 1 J=2,4
	DO 1 I=1,J-1
	K(I,J)=(M(I)**2+M(J)**2-P(I,J))/M(I)/M(J)
	IF (I.EQ.3) K(I,J)=K(I,J)-M(I)/M(J)
	IF (J.EQ.3) K(I,J)=K(I,J)-M(J)/M(I)
	R(I,J) =SQE(DCMPLX(1D0,0D0),DCMPLX(-K(I,J),0D0),
     *	            DCMPLX(1D0,0D0))
	IF( DIMAG(R(I,J)).EQ.0D0 ) THEN
	   RS(I,J)=SQE(DCMPLX(1D0,0D0),DCMPLX(-K(I,J),EPS),
     *	               DCMPLX(1D0,0D0))
	ELSE
	   RS(I,J)=R(I,J)
	ENDIF
1	CONTINUE
	SS(1)=RS(1,2)
	SS(2)=RS(2,3)
	SS(3)=RS(3,4)
	SS(4)=RS(1,4)
	AA=K(3,4)/R(2,4)-K(2,3)
	BB=K(1,3)*(1D0/R(2,4)-R(2,4))+K(1,2)*K(3,4)-K(1,4)*K(2,3)
	CC=K(1,2)*K(1,3)-K(1,3)*K(1,4)*R(2,4)+R(2,4)*K(3,4)-K(2,3)
	DD=K(2,3)-R(2,4)*K(3,4)
	XX(1)=SQE(AA,BB,CC+IEPS*DD)
	XX(2)=(CC+IEPS*DD)/AA/XX(1)
	DO 2 I=1,2
	X(I,1)=XX(I)/R(2,4)
	X(I,2)=XX(I)/R(2,4)*R(1,3)
	X(I,3)=XX(I)*R(1,3)
	X(I,4)=XX(I)
2	CONTINUE
	D04 = DCMPLX(0D0,0D0)
	DO 3 I=1,2
	D04 = D04 + (2D0*I-3D0)*(
     *		LI2(1D0+SS(4)*X(I,4))
     *	       -LI2(1D0+SS(1)*X(I,1))
     *	       +LI2(1D0+X(I,4)/SS(4))
     *	       -LI2(1D0+X(I,1)/SS(1))
     *	       +ETA(-X(I,4),SS(4))*LOG(1D0+SS(4)*X(I,4))
     *	       -ETA(-X(I,1),SS(1))*LOG(1D0+SS(1)*X(I,1))
     *	       +ETA(-X(I,4),1D0/SS(4))*LOG(1D0+X(I,4)/SS(4))
     *	       -ETA(-X(I,1),1D0/SS(1))*LOG(1D0+X(I,1)/SS(1))
     *	       -LI2(1D0+X(I,4)*(K(3,4)-IEPS)/(K(1,3)-IEPS))
     *	       +LI2(1D0+X(I,1)*(K(2,3)-IEPS)/(K(1,3)-IEPS))
     *	       -ETA(-X(I,4),(K(3,4)-IEPS)/(K(1,3)-IEPS))
     *	           *LOG(1D0+X(I,4)*(K(3,4)-IEPS)/(K(1,3)-IEPS))
     *	       +ETA(-X(I,1),(K(2,3)-IEPS)/(K(1,3)-IEPS))
     *	           *LOG(1D0+X(I,1)*(K(2,3)-IEPS)/(K(1,3)-IEPS))   )
	IF (DIMAG(R(2,4)).NE.0D0) THEN
	   H=ETA(-1D0/XX(I),R(2,4))
	ELSE
	   H=DCMPLX(0D0,0D0)
	   IF (DREAL(R(2,4)).LT.0D0) THEN
	      HH=-1D0/XX(I)
	      IM1=DIMAG(HH)
	      IM2=DIMAG(RS(2,4))
              pi = 4.D0*datan(1.D0)
	      IF ((IM1.GT.0D0).AND.(IM2.GT.0D0)) THEN
	         H=-DCMPLX(0D0,2D0*PI)
	      ENDIF
	      IF ((IM1.LT.0D0).AND.(IM2.LT.0D0)) THEN
	         H=+DCMPLX(0D0,2D0*PI)
	      ENDIF
	   ENDIF
	ENDIF
	D04 = D04 + (2D0*I-3D0)*
     *	          H*( LOG( (K(1,2)-R(2,4)*K(1,4)
     *			  +XX(I)*(1D0/R(2,4)-R(2,4)))/DD )
     *		     +LOG(K(1,3)-IEPS) )
3	CONTINUE
	D04 = D04/M(1)/M(2)/M(3)/M(4)/AA/(XX(1)-XX(2))
	RETURN

C	****** ONLY MM2 AND MM3 ARE ZERO ******
40	CONTINUE
	EPS=1D-17
	IEPS=DCMPLX(0D0,EPS)

	M(1)=MM1
	M(2)=10D0
	M(3)=10D0
	M(4)=MM4
	P(1,2)=Q1
	P(1,3)=Q12
	P(1,4)=Q4
	P(2,3)=Q2
	P(2,4)=Q23
	P(3,4)=Q3
	DO 4 J=2,4
	DO 4 I=1,J-1
	K(I,J)=(M(I)**2+M(J)**2-P(I,J))/M(I)/M(J)
	IF (I.EQ.2) K(I,J)=K(I,J)-M(I)/M(J)
	IF (J.EQ.2) K(I,J)=K(I,J)-M(J)/M(I)
	IF (I.EQ.3) K(I,J)=K(I,J)-M(I)/M(J)
	IF (J.EQ.3) K(I,J)=K(I,J)-M(J)/M(I)
	R(I,J) =SQE(DCMPLX(1D0,0D0),DCMPLX(-K(I,J),0D0),
     *	            DCMPLX(1D0,0D0))
	IF( DIMAG(R(I,J)).EQ.0D0 ) THEN
	   RS(I,J)=SQE(DCMPLX(1D0,0D0),DCMPLX(-K(I,J),EPS),
     *	               DCMPLX(1D0,0D0))
	ELSE
	   RS(I,J)=R(I,J)
	ENDIF
4	CONTINUE
	SS(1)=RS(1,2)
	SS(2)=RS(2,3)
	SS(3)=RS(3,4)
	SS(4)=RS(1,4)
	AA=K(2,4)*K(3,4)-K(2,3)
	BB=K(1,3)*K(2,4)+K(1,2)*K(3,4)-K(1,4)*K(2,3)
	CC=K(1,2)*K(1,3)-K(2,3)
	DD=K(2,3)
	XX(1)=SQE(AA,BB,CC+IEPS*DD)
	XX(2)=(CC+IEPS*DD)/AA/XX(1)
	DO 5 I=1,2
	X(I,1)=XX(I)/R(2,4)
	X(I,2)=XX(I)/R(2,4)*R(1,3)
	X(I,3)=XX(I)*R(1,3)
	X(I,4)=XX(I)
5	CONTINUE
	D04 = DCMPLX(0D0,0D0)
	DO 6 I=1,2
	D04 = D04 + (2D0*I-3D0)*(
     *		LI2(1D0+SS(4)*X(I,4))
     *	       +LI2(1D0+X(I,4)/SS(4))
     *	       +ETA(-X(I,4),SS(4))*LOG(1D0+SS(4)*X(I,4))
     *	       +ETA(-X(I,4),1D0/SS(4))*LOG(1D0+X(I,4)/SS(4))
     *	       -LI2(1D0+XX(I)*(K(3,4)-IEPS)/(K(1,3)-IEPS))
     *	       -LI2(1D0+XX(I)*(K(2,4)-IEPS)/(K(1,2)-IEPS))
     *	       -ETA(-XX(I),(K(3,4)-IEPS)/(K(1,3)-IEPS))
     *	           *LOG(1D0+XX(I)*(K(3,4)-IEPS)/(K(1,3)-IEPS))
     *	       -ETA(-XX(I),(K(2,4)-IEPS)/(K(1,2)-IEPS))
     *	           *LOG(1D0+XX(I)*(K(2,4)-IEPS)/(K(1,2)-IEPS))
     *	       +LOG(-XX(I))*( LOG(K(1,2)-IEPS)
     *			     +LOG(K(1,3)-IEPS)-LOG(K(2,3)-IEPS) ) )
6	CONTINUE
	D04 = D04/M(1)/M(2)/M(3)/M(4)/AA/(XX(1)-XX(2))

	RETURN

	END











c FC %************************************************************************
c FC %        FUNCTION D04(P1,P2,P3,P4,P12,P23,M1,M2,M3,M4)
c FC %************************************************************************
c FC %*  SCALAR 4-POINT FUNCTION WITH AT LEAST ONE MASS ZERO                 *
c FC %*  P1,P2,P3,P4 = SQUARED EXTERNAL MOMENTA			       *
c FC %*  P12 = (p1+p2)**2,  P23 = (p2+p3)**2				       *
c FC %*----------------------------------------------------------------------*
c FC %*  2.1.92  SD	         					       *
c FC %*  rearrangement to exploit massless external momenta   14.3.01  DZ    *  
c FC %************************************************************************
c FC %        IMPLICIT REAL*8 (A-Z)
c FC %	REAL*8 M(4),P(4,4),K(4,4)
c FC %        real*8 pi,eps
c FC %        real*8 im1,im2
c FC %        real*8 m1,m2,m3,m4
c FC %        real*8 m02,m12,m22,m32,m42
c FC %        real*8 mm0,mm1,mm2,mm3,mm4 
c FC %        real*8 p1,p2,p3,p4,p12,p23
c FC %        real*8 q0,q1,q2,q3,q4,q00,q12,q23
c FC %	COMPLEX*16 A1,A2,A3,A4,SWAP
c FC %	COMPLEX*16 SS(4), XX(2), X(2,4),RS(4,4)
c FC %	COMPLEX*16 S0(4),XX0(2),X0(2,4), R(4,4),G(2)
c FC %        COMPLEX*16 D04,LI2,ETA,SQE,ETAS
c FC %	COMPLEX*16 AA,BB,CC,DD,IEPS,H,HH,L1,L2,L3,L4
c FC %	INTEGER I,J, i1,i2,i3,i4
c FC %
c FC %        MM1=M1
c FC %        MM2=M2
c FC %        MM3=M3
c FC %        MM4=M4
c FC %        M12=M1*M1
c FC %        M22=M2*M2
c FC %        M32=M3*M3
c FC %        M42=M4*M4
c FC %        Q1=P1
c FC %        Q2=P2
c FC %        Q3=P3
c FC %	Q4=P4
c FC %        Q12=P12
c FC %        Q23=P23
c FC %
c FC %C	IS AT LEAST ONE MASS ZERO ???
c FC %	IF (MM1*MM2*MM3*MM4.NE.0D0) GOTO 130
c FC %
c FC %C	PERMUTATE UNTIL MM3=0D0
c FC %	GOTO 20
c FC %10	CONTINUE
c FC %	MM0=MM1
c FC %	MM1=MM2
c FC %	MM2=MM3
c FC %	MM3=MM4
c FC %	MM4=MM0
c FC %	M02=M12
c FC %	M12=M22
c FC %	M22=M32
c FC %	M32=M42
c FC %	M42=M02
c FC %	Q00=Q12
c FC %	Q12=Q23
c FC %	Q23=Q00
c FC %	Q0=Q1
c FC %	Q1=Q2
c FC %	Q2=Q3
c FC %	Q3=Q4
c FC %	Q4=Q0
c FC %20	IF (MM3.NE.0D0) GOTO 10
c FC %C	ONLY MM3 IS ZERO
c FC %	IF (MM1*MM2*MM4.NE.0D0) GOTO 30
c FC %C	ONLY MM3 AND MM4 ARE ZERO ==> 3->2, 4->3...
c FC %	IF ((MM1*MM2.NE.0D0).AND.(MM4.EQ.0D0)) GOTO 10
c FC %C	ONLY MM2 AND MM3 ARE ZERO
c FC %	IF ((MM1*MM4.NE.0D0).AND.(MM2.EQ.0D0)) GOTO 40
c FC %
c FC %C	WRITE(*,*)'CASE OF THIS SPECIAL D0-FUNCTION NOT IMPLEMENTED!'
c FC %C	STOP
c FC %
c FC %C	****** NO MASS EQUAL TO ZERO ******
c FC %130	CONTINUE
c FC %	EPS=1D-18
c FC %	IEPS=DCMPLX(0D0,EPS)
c FC %c check for massless external momentum: excellent candidate for p13,
c FC %c leading to r13 >=1 and real.
c FC %        if (q1*q2*q3*q4.eq.0d0) then
c FC %           if (q2.eq.0d0) then
c FC %              I1 = 2
c FC %              I2 = 3
c FC %              I3 = 1
c FC %              I4 = 4
c FC %           elseif (q1.eq.0d0) then
c FC %              I1 = 1
c FC %              I2 = 3
c FC %              I3 = 2
c FC %              I4 = 4
c FC %           elseif (q3.eq.0d0) then
c FC %              I1 = 2
c FC %              I2 = 4
c FC %              I3 = 1
c FC %              I4 = 3
c FC %           else
c FC %              I1 = 1
c FC %              I2 = 4
c FC %              I3 = 2
c FC %              I4 = 3
c FC %           endif
c FC %           M(i1)=MM1
c FC %           M(i2)=MM2
c FC %           M(i3)=MM3
c FC %           M(i4)=MM4
c FC %           P(i1,i2)=Q1
c FC %           P(i3,i2)=Q2
c FC %           P(i3,i4)=Q3
c FC %           P(i1,i4)=Q4
c FC %           P(i1,i3)=Q12
c FC %           P(i3,i1)=Q12
c FC %           P(i2,i4)=Q23
c FC %           P(i4,i2)=Q23
c FC %	ELSEIF( ABS((MM1**2+MM3**2-Q12)/MM1/MM3).LT.2D0 ) THEN
c FC %C	R13 WOULD BE NOT REAL. -> PERMUTATION! -> R(2,4) IS NOT REAL.
c FC %	   M(1)=MM2
c FC %	   M(2)=MM3
c FC %	   M(3)=MM4
c FC %	   M(4)=MM1
c FC %	   P(1,2)=Q2
c FC %	   P(1,3)=Q23
c FC %	   P(1,4)=Q1
c FC %	   P(2,3)=Q3
c FC %	   P(2,4)=Q12
c FC %	   P(3,4)=Q4
c FC %	ELSE
c FC %C	R(1,3) IS REAL.
c FC %	   M(1)=MM1
c FC %	   M(2)=MM2
c FC %	   M(3)=MM3
c FC %	   M(4)=MM4
c FC %	   P(1,2)=Q1
c FC %	   P(1,3)=Q12
c FC %	   P(1,4)=Q4
c FC %	   P(2,3)=Q2
c FC %	   P(2,4)=Q23
c FC %	   P(3,4)=Q3
c FC %	ENDIF
c FC %
c FC %	DO 11 J=2,4
c FC %	DO 11 I=1,J-1
c FC %	K(I,J)=(M(I)**2+M(J)**2-P(I,J))/M(I)/M(J)
c FC %	R(I,J) =SQE(DCMPLX(1D0,0D0),DCMPLX(-K(I,J),0D0),
c FC %     *	            DCMPLX(1D0,0D0))
c FC %        IF( DBLE(K(I,J)).LT.-2D0 ) THEN
c FC %c        IF( DIMAG(R(I,J)).EQ.0D0 ) THEN
c FC %	   RS(I,J)=SQE(DCMPLX(1D0,0D0),DCMPLX(-K(I,J),EPS),
c FC %     *	               DCMPLX(1D0,0D0))
c FC %	ELSE
c FC %	   RS(I,J)=R(I,J)
c FC %	ENDIF
c FC %11	CONTINUE
c FC %
c FC %	SS(1)=RS(1,2)
c FC %	SS(2)=RS(2,3)
c FC %	SS(3)=RS(3,4)
c FC %	SS(4)=RS(1,4)
c FC %	S0(1)=R(1,2)
c FC %	S0(2)=R(2,3)
c FC %	S0(3)=R(3,4)
c FC %	S0(4)=R(1,4)
c FC %	AA=K(3,4)/R(2,4)+R(1,3)*K(1,2)-K(1,4)*R(1,3)/R(2,4)-K(2,3)
c FC %	BB=(R(2,4)-1D0/R(2,4))*(R(1,3)-1D0/R(1,3))
c FC %     *		+K(1,2)*K(3,4)-K(1,4)*K(2,3)
c FC %	CC=K(1,2)/R(1,3)+R(2,4)*K(3,4)-K(1,4)*R(2,4)/R(1,3)-K(2,3)
c FC %	DD=K(2,3)-R(1,3)*K(1,2)-R(2,4)*K(3,4)+R(1,3)*R(2,4)*K(1,4)
c FC %	XX(1)=SQE(AA,BB,CC+IEPS*DD)
c FC %	XX(2)=(CC+IEPS*DD)/AA/XX(1)
c FC %	XX0(1)=SQE(AA,BB,CC)
c FC %	XX0(2)=CC/AA/XX0(1)
c FC %c	IF (ABS(DREAL(XX0(1)-XX(2))).LT.ABS(DREAL(XX0(1)-XX(1)))) THEN
c FC %	IF (ABS(XX0(1)-XX(2)).LT.ABS(XX0(1)-XX(1))) THEN
c FC %	  SWAP  =XX0(1)
c FC %	  XX0(1)=XX0(2)
c FC %	  XX0(2)=SWAP
c FC %	ENDIF
c FC %
c FC %	DO 12 I=1,2
c FC %	G(I)  =SIGN( 1D0,DREAL(AA*(XX(I)-XX(3-I))) )
c FC %	 X(I,1)= XX(I)/R(2,4)
c FC %	X0(I,1)=XX0(I)/R(2,4)
c FC %	 X(I,2)= XX(I)/R(2,4)*R(1,3)
c FC %	X0(I,2)=XX0(I)/R(2,4)*R(1,3)
c FC %	 X(I,3)= XX(I)*R(1,3)
c FC %	X0(I,3)=XX0(I)*R(1,3)
c FC %	 X(I,4)= XX(I)
c FC %	X0(I,4)=XX0(I)
c FC %12	CONTINUE
c FC %
c FC %	D04 = DCMPLX(0D0,0D0)
c FC %	DO 13 I=1,2
c FC %	DO 13 J=1,4
c FC %	A1 = 1D0+X0(I,J)*S0(J) + ABS(1D0+X0(I,J)*S0(J))*IEPS*
c FC %     *				  SIGN(1D0,DIMAG(X(I,J)*SS(J)))
c FC %	A2 = 1D0+X0(I,J)/S0(J) + ABS(1D0+X0(I,J)/S0(J))*IEPS*
c FC %     *				  SIGN(1D0,DIMAG(X(I,J)/SS(J)))
c FC %	D04 = D04 + (-1D0)**(I+J)*(
c FC %     *		LI2(A1)+ETA(-X(I,J),SS(J))*LOG(A1)
c FC %     *	       +LI2(A2)+ETA(-X(I,J),1D0/SS(J))*LOG(A2)     )
c FC %13	CONTINUE
c FC %
c FC %	IF( DIMAG(R(1,3)).EQ.0D0 ) THEN
c FC %	DO 14 I=1,2
c FC %	   A1 = (K(1,3)-2D0*R(1,3))/XX0(I)
c FC %     *		      -R(1,3)*K(1,4)+K(3,4)
c FC %     	   A2 = ((K(2,4)-2D0*R(2,4))*R(1,3)*XX0(I)
c FC %     *		      -R(2,4)*K(3,4)+K(2,3))/DD
c FC %	   A3 = (K(1,3)-2D0*R(1,3))*R(2,4)/XX0(I)
c FC %     *		      -R(1,3)*K(1,2)+K(2,3)
c FC %	   A4 = ((K(2,4)-2D0*R(2,4))*XX0(I)
c FC %     *		      -R(2,4)*K(1,4)+K(1,2))/DD
c FC %	   L1 = LOG( A1-ABS(A1)*IEPS )
c FC %     	   L2 = LOG( A2+ABS(A2)*IEPS*G(I)*SIGN(1D0,DREAL(R(1,3))
c FC %     *				        	  *DIMAG(RS(2,4))) )
c FC %	   L3 = LOG( A3-ABS(A3)*IEPS )
c FC %	   L4 = LOG( A4+ABS(A4)*IEPS*G(I)*SIGN(1D0,DIMAG(RS(2,4))) )
c FC %
c FC %	   D04 = D04 
c FC %     &         + (3D0-2D0*I)*(
c FC %     *		       ETAS( -XX(I),R(1,3),RS(1,3) )
c FC %     *		          *( LOG(R(1,3)*XX(I)) + L1 + L2 )
c FC %     *		     + ETAS( -XX(I),1D0/R(2,4),1D0/RS(2,4) )
c FC %     *		          *( LOG(XX(I)/R(2,4)) + L3 + L4 )
c FC %     *		     - ( ETAS( -XX(I),R(1,3)/R(2,4),RS(1,3)/RS(2,4) )
c FC %     *		       + ETA( RS(1,3),1D0/RS(2,4) )                  )
c FC %     *		        *( LOG(XX(I)*R(1,3)/R(2,4)) + L3 + L2 )
c FC %     *	  	     + ETA( RS(1,3),1D0/RS(2,4) )
c FC %     *		       *ETAS(-XX(I),-R(1,3)/R(2,4),-RS(1,3)/RS(2,4))   )
c FC %14	CONTINUE
c FC %	ELSE
c FC %	DO 15 I=1,2
c FC %	   L1 = LOG( R(2,4)/XX0(I)+XX0(I)/R(2,4)+K(1,2)
c FC %     *		     -XX0(I)/R(2,4)*EPS*BB*G(I) )
c FC %	   L2 = LOG( R(1,3)*XX0(I)+1D0/XX0(I)/R(1,3)+K(3,4)
c FC %     *		     -XX0(I)*R(1,3)*EPS*BB*G(I) )
c FC %	   L3 = LOG( R(1,3)/R(2,4)*XX0(I)+R(2,4)/XX0(I)/R(1,3)+K(2,3)
c FC %     *		     -XX0(I)*R(1,3)/R(2,4)*EPS*BB*G(I) )
c FC %
c FC %	   D04 = D04 
c FC %     &          + (3D0-2D0*I)*(
c FC %     *		     ETA(-XX(I),1D0/R(2,4))
c FC %     *		      *( LOG(XX(I)/R(2,4)) + L1 )
c FC %     *		    +ETA(-XX(I),R(1,3))
c FC %     *		      *( LOG(R(1,3)*XX(I)) + L2 )
c FC %     *		    -( ETA(-XX(I),R(1,3)/R(2,4))
c FC %     *		      +ETA(R(1,3),1D0/R(2,4)) )
c FC %     *		      *( LOG(XX(I)*R(1,3)/R(2,4)) + L3 )
c FC %     *	  	    +ETA(R(1,3),1D0/R(2,4))
c FC %     *		      *ETA(-XX(I),-R(1,3)/R(2,4))
c FC %     *		      *(1D0-G(I)*SIGN(1D0,DREAL(BB)))	    )
c FC %15	CONTINUE
c FC %	ENDIF
c FC %
c FC %	D04 = D04/M(1)/M(2)/M(3)/M(4)/AA/(XX(1)-XX(2))
c FC %	RETURN
c FC %
c FC %
c FC %C--->	***************** SPEZIELL ( --> T.SACK-PROMOTION )
c FC %C	D1=Q12-M12
c FC %C	D2=Q2 -M22
c FC %C	D3=Q3 -M42
c FC %C	IF ((D1*D2.LE.0D0).OR.(D2*D3.LE.0D0)) THEN
c FC %C	   WRITE(*,*) 'THE CASE OF DIFFERENT SIGNS OF THE D1,D2,D3'
c FC %C	   WRITE(*,*) 'IN D04(...) IS NOT IMPLEMENTED !!!'
c FC %C	   STOP
c FC %C	ENDIF
c FC %C	NM1=ABS(MM1/D1)
c FC %C	NM2=ABS(MM2/D2)
c FC %C	NM3=ABS(MM4/D3)
c FC %C	NP1=Q2/D2**2+Q12/D1**2+(Q1-Q2-Q12)/D1/D2
c FC %C	NP2=Q2/D2**2+ Q3/D3**2+(Q23-Q2-Q3)/D2/D3
c FC %C	NP3=Q3/D3**2+Q12/D1**2+(Q4-Q3-Q12)/D1/D3
c FC %C	D04=C04(NP1,NP2,NP3,NM1,NM2,NM3)/D1/D2/D3
c FC %
c FC %C	*************** ALLGEMEIN
c FC %
c FC %
c FC %C	****** ONLY MM3 IS ZERO ******
c FC %30	CONTINUE
c FC %	EPS=1D-17
c FC %	IEPS=DCMPLX(0D0,EPS)
c FC %	M(1)=MM1
c FC %	M(2)=MM2
c FC %	M(3)=10D0
c FC %	M(4)=MM4
c FC %	P(1,2)=Q1
c FC %	P(1,3)=Q12
c FC %	P(1,4)=Q4
c FC %	P(2,3)=Q2
c FC %	P(2,4)=Q23
c FC %	P(3,4)=Q3
c FC %	DO 1 J=2,4
c FC %	DO 1 I=1,J-1
c FC %	K(I,J)=(M(I)**2+M(J)**2-P(I,J))/M(I)/M(J)
c FC %	IF (I.EQ.3) K(I,J)=K(I,J)-M(I)/M(J)
c FC %	IF (J.EQ.3) K(I,J)=K(I,J)-M(J)/M(I)
c FC %	R(I,J) =SQE(DCMPLX(1D0,0D0),DCMPLX(-K(I,J),0D0),
c FC %     *	            DCMPLX(1D0,0D0))
c FC %	IF( DIMAG(R(I,J)).EQ.0D0 ) THEN
c FC %	   RS(I,J)=SQE(DCMPLX(1D0,0D0),DCMPLX(-K(I,J),EPS),
c FC %     *	               DCMPLX(1D0,0D0))
c FC %	ELSE
c FC %	   RS(I,J)=R(I,J)
c FC %	ENDIF
c FC %1	CONTINUE
c FC %	SS(1)=RS(1,2)
c FC %	SS(2)=RS(2,3)
c FC %	SS(3)=RS(3,4)
c FC %	SS(4)=RS(1,4)
c FC %	AA=K(3,4)/R(2,4)-K(2,3)
c FC %	BB=K(1,3)*(1D0/R(2,4)-R(2,4))+K(1,2)*K(3,4)-K(1,4)*K(2,3)
c FC %	CC=K(1,2)*K(1,3)-K(1,3)*K(1,4)*R(2,4)+R(2,4)*K(3,4)-K(2,3)
c FC %	DD=K(2,3)-R(2,4)*K(3,4)
c FC %	XX(1)=SQE(AA,BB,CC+IEPS*DD)
c FC %	XX(2)=(CC+IEPS*DD)/AA/XX(1)
c FC %	DO 2 I=1,2
c FC %	X(I,1)=XX(I)/R(2,4)
c FC %	X(I,2)=XX(I)/R(2,4)*R(1,3)
c FC %	X(I,3)=XX(I)*R(1,3)
c FC %	X(I,4)=XX(I)
c FC %2	CONTINUE
c FC %	D04 = DCMPLX(0D0,0D0)
c FC %	DO 3 I=1,2
c FC %	D04 = D04 + (2D0*I-3D0)*(
c FC %     *		LI2(1D0+SS(4)*X(I,4))
c FC %     *	       -LI2(1D0+SS(1)*X(I,1))
c FC %     *	       +LI2(1D0+X(I,4)/SS(4))
c FC %     *	       -LI2(1D0+X(I,1)/SS(1))
c FC %     *	       +ETA(-X(I,4),SS(4))*LOG(1D0+SS(4)*X(I,4))
c FC %     *	       -ETA(-X(I,1),SS(1))*LOG(1D0+SS(1)*X(I,1))
c FC %     *	       +ETA(-X(I,4),1D0/SS(4))*LOG(1D0+X(I,4)/SS(4))
c FC %     *	       -ETA(-X(I,1),1D0/SS(1))*LOG(1D0+X(I,1)/SS(1))
c FC %     *	       -LI2(1D0+X(I,4)*(K(3,4)-IEPS)/(K(1,3)-IEPS))
c FC %     *	       +LI2(1D0+X(I,1)*(K(2,3)-IEPS)/(K(1,3)-IEPS))
c FC %     *	       -ETA(-X(I,4),(K(3,4)-IEPS)/(K(1,3)-IEPS))
c FC %     *	           *LOG(1D0+X(I,4)*(K(3,4)-IEPS)/(K(1,3)-IEPS))
c FC %     *	       +ETA(-X(I,1),(K(2,3)-IEPS)/(K(1,3)-IEPS))
c FC %     *	           *LOG(1D0+X(I,1)*(K(2,3)-IEPS)/(K(1,3)-IEPS))   )
c FC %	IF (DIMAG(R(2,4)).NE.0D0) THEN
c FC %	   H=ETA(-1D0/XX(I),R(2,4))
c FC %	ELSE
c FC %	   H=DCMPLX(0D0,0D0)
c FC %	   IF (DREAL(R(2,4)).LT.0D0) THEN
c FC %	      HH=-1D0/XX(I)
c FC %	      IM1=DIMAG(HH)
c FC %	      IM2=DIMAG(RS(2,4))
c FC %              pi = 4.D0*datan(1.D0)
c FC %	      IF ((IM1.GT.0D0).AND.(IM2.GT.0D0)) THEN
c FC %	         H=-DCMPLX(0D0,2D0*PI)
c FC %	      ENDIF
c FC %	      IF ((IM1.LT.0D0).AND.(IM2.LT.0D0)) THEN
c FC %	         H=+DCMPLX(0D0,2D0*PI)
c FC %	      ENDIF
c FC %	   ENDIF
c FC %	ENDIF
c FC %	D04 = D04 + (2D0*I-3D0)*
c FC %     *	          H*( LOG( (K(1,2)-R(2,4)*K(1,4)
c FC %     *			  +XX(I)*(1D0/R(2,4)-R(2,4)))/DD )
c FC %     *		     +LOG(K(1,3)-IEPS) )
c FC %3	CONTINUE
c FC %	D04 = D04/M(1)/M(2)/M(3)/M(4)/AA/(XX(1)-XX(2))
c FC %	RETURN
c FC %
c FC %C	****** ONLY MM2 AND MM3 ARE ZERO ******
c FC %40	CONTINUE
c FC %	EPS=1D-17
c FC %	IEPS=DCMPLX(0D0,EPS)
c FC %
c FC %	M(1)=MM1
c FC %	M(2)=10D0
c FC %	M(3)=10D0
c FC %	M(4)=MM4
c FC %	P(1,2)=Q1
c FC %	P(1,3)=Q12
c FC %	P(1,4)=Q4
c FC %	P(2,3)=Q2
c FC %	P(2,4)=Q23
c FC %	P(3,4)=Q3
c FC %	DO 4 J=2,4
c FC %	DO 4 I=1,J-1
c FC %	K(I,J)=(M(I)**2+M(J)**2-P(I,J))/M(I)/M(J)
c FC %	IF (I.EQ.2) K(I,J)=K(I,J)-M(I)/M(J)
c FC %	IF (J.EQ.2) K(I,J)=K(I,J)-M(J)/M(I)
c FC %	IF (I.EQ.3) K(I,J)=K(I,J)-M(I)/M(J)
c FC %	IF (J.EQ.3) K(I,J)=K(I,J)-M(J)/M(I)
c FC %	R(I,J) =SQE(DCMPLX(1D0,0D0),DCMPLX(-K(I,J),0D0),
c FC %     *	            DCMPLX(1D0,0D0))
c FC %	IF( DIMAG(R(I,J)).EQ.0D0 ) THEN
c FC %	   RS(I,J)=SQE(DCMPLX(1D0,0D0),DCMPLX(-K(I,J),EPS),
c FC %     *	               DCMPLX(1D0,0D0))
c FC %	ELSE
c FC %	   RS(I,J)=R(I,J)
c FC %	ENDIF
c FC %4	CONTINUE
c FC %	SS(1)=RS(1,2)
c FC %	SS(2)=RS(2,3)
c FC %	SS(3)=RS(3,4)
c FC %	SS(4)=RS(1,4)
c FC %	AA=K(2,4)*K(3,4)-K(2,3)
c FC %	BB=K(1,3)*K(2,4)+K(1,2)*K(3,4)-K(1,4)*K(2,3)
c FC %	CC=K(1,2)*K(1,3)-K(2,3)
c FC %	DD=K(2,3)
c FC %	XX(1)=SQE(AA,BB,CC+IEPS*DD)
c FC %	XX(2)=(CC+IEPS*DD)/AA/XX(1)
c FC %	DO 5 I=1,2
c FC %	X(I,1)=XX(I)/R(2,4)
c FC %	X(I,2)=XX(I)/R(2,4)*R(1,3)
c FC %	X(I,3)=XX(I)*R(1,3)
c FC %	X(I,4)=XX(I)
c FC %5	CONTINUE
c FC %	D04 = DCMPLX(0D0,0D0)
c FC %	DO 6 I=1,2
c FC %	D04 = D04 + (2D0*I-3D0)*(
c FC %     *		LI2(1D0+SS(4)*X(I,4))
c FC %     *	       +LI2(1D0+X(I,4)/SS(4))
c FC %     *	       +ETA(-X(I,4),SS(4))*LOG(1D0+SS(4)*X(I,4))
c FC %     *	       +ETA(-X(I,4),1D0/SS(4))*LOG(1D0+X(I,4)/SS(4))
c FC %     *	       -LI2(1D0+XX(I)*(K(3,4)-IEPS)/(K(1,3)-IEPS))
c FC %     *	       -LI2(1D0+XX(I)*(K(2,4)-IEPS)/(K(1,2)-IEPS))
c FC %     *	       -ETA(-XX(I),(K(3,4)-IEPS)/(K(1,3)-IEPS))
c FC %     *	           *LOG(1D0+XX(I)*(K(3,4)-IEPS)/(K(1,3)-IEPS))
c FC %     *	       -ETA(-XX(I),(K(2,4)-IEPS)/(K(1,2)-IEPS))
c FC %     *	           *LOG(1D0+XX(I)*(K(2,4)-IEPS)/(K(1,2)-IEPS))
c FC %     *	       +LOG(-XX(I))*( LOG(K(1,2)-IEPS)
c FC %     *			     +LOG(K(1,3)-IEPS)-LOG(K(2,3)-IEPS) ) )
c FC %6	CONTINUE
c FC %	D04 = D04/M(1)/M(2)/M(3)/M(4)/AA/(XX(1)-XX(2))
c FC %
c FC %	RETURN
c FC %
c FC %	END



***********************************************************************
        FUNCTION ETA(C1,C2)                                            
***********************************************************************
*       COMPLEX ETA-FUNKTION                                           
*---------------------------------------------------------------------*
*       8.06.90    ANSGAR DENNER                                       
***********************************************************************
        IMPLICIT   LOGICAL(A-Z)                                        
        COMPLEX*16 ETA,C1,C2                                           
        REAL*8     PI,IM1,IM2,IM12                                     
                                                                       
        PI     = 4D0*DATAN(1D0)                                        
        IM1    = DIMAG(C1)                                             
        IM2    = DIMAG(C2)                                             
        IM12   = DIMAG(C1*C2)                                          
                                                                       
        IF(IM1.LT.0D0.AND.IM2.LT.0D0.AND.IM12.GT.0D0) THEN             
            ETA = DCMPLX(0D0,2D0*PI)                                   
        ELSE IF (IM1.GT.0D0.AND.IM2.GT.0D0.AND.IM12.LT.0D0) THEN       
            ETA = DCMPLX(0D0,-2D0*PI)                                  
        ELSE                                                           
            ETA = DCMPLX(0D0)                                          
        END IF                                                         
        END                                                            
***********************************************************************
        FUNCTION ETAS(Y,R,RS)                                            
***********************************************************************
*       MODIFIED ETA-FUNKTION                                           
*---------------------------------------------------------------------*
*       18.1.94   SD                                       
***********************************************************************
        IMPLICIT   LOGICAL(A-Z)                                        
        COMPLEX*16 ETA,ETAS,Y,R,RS
        REAL*8     PI,IMY,IMRS
                                                                       
        PI     = 4D0*DATAN(1D0)                                        

	IF( DIMAG(R).NE.0D0 ) THEN
	    ETAS = ETA(Y,R)
	ELSE	    
	    IF( DREAL(R).GT.0D0 ) THEN
		ETAS = DCMPLX(0D0,0D0)
	    ELSE
	 	IMY  = DIMAG(Y)
		IMRS = DIMAG(RS)
		ETAS = 2D0*DCMPLX(0D0,PI)*(
     *			(1D0+SIGN(1D0,-IMY))*(1D0+SIGN(1D0,-IMRS))-
     *			(1D0+SIGN(1D0, IMY))*(1D0+SIGN(1D0, IMRS))
     *					  )/4D0
	    ENDIF
	ENDIF
        END  


c FC %***********************************************************************
c FC %        FUNCTION ETAS(Y,R,RS)                                            
c FC %***********************************************************************
c FC %*       MODIFIED ETA-FUNKTION                                           
c FC %*---------------------------------------------------------------------*
c FC %*       18.1.94   SD                                       
c FC %***********************************************************************
c FC %        IMPLICIT   LOGICAL(A-Z)                                        
c FC %        COMPLEX*16 ETA,ETAS,Y,R,RS
c FC %        REAL*8     PI,IMY,IMRS
c FC %                                                                       
c FC %        PI     = 4D0*DATAN(1D0)                                        
c FC %
c FC %	IF( DIMAG(R).NE.0D0 ) THEN
c FC %	    ETAS = ETA(Y,R)
c FC %	ELSE	    
c FC %	    IF( DREAL(R).GT.0D0 ) THEN
c FC %		ETAS = DCMPLX(0D0,0D0)
c FC %	    ELSE
c FC %	 	IMY  = DIMAG(Y)
c FC %		IMRS = DIMAG(RS)
c FC %		ETAS = 2D0*DCMPLX(0D0,PI)*(
c FC %     *			(1D0+SIGN(1D0,-IMY))*(1D0+SIGN(1D0,-IMRS))-
c FC %     *			(1D0+SIGN(1D0, IMY))*(1D0+SIGN(1D0, IMRS))
c FC %     *					  )/4D0
c FC %	    ENDIF
c FC %	ENDIF
c FC %        END                                                            

***********************************************************************
        FUNCTION SQE(A,B,C)                                            
***********************************************************************
*       SOLUTION OF QUADRATIC EQUATION				      *
*---------------------------------------------------------------------*
*       13.1.92  SD						      *
***********************************************************************
        IMPLICIT REAL*8 (A-Z)                                        
        COMPLEX*16 A,B,C,SQE,X1,X2

	X1=(-B+SQRT(B**2-4D0*A*C))/2D0/A
	X2=(-B-SQRT(B**2-4D0*A*C))/2D0/A

	IF (ABS(X1).GT.ABS(X2)) THEN
	   SQE=X1
	ELSE
	   SQE=X2
	ENDIF

        END                                                            


ccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
cccccccccccccccccccccccccccccccccccccccccccccccccccccccc
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
cccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c
c
c
c
      function E01M(m,p1,p2,p3,p4,p1p2,p1p3,p1p4,p2p3,p2p4,p3p4,
     #     D0_2345,D0_1345,D0_1245,D0_1235,D0_1234)
      implicit none
      complex * 16 E01M
      real * 8 p1,p2,p3,p4,p1p2,p1p3,p1p4,p2p3,p2p4,p3p4,m,msq
      real * 8  tempX5,ItempX5,tempX4(0:4)
      complex * 16 D0_2345,D0_1345,D0_1245,D0_1235,D0_1234

       msq=m*m

       tempX5=-8*(p1**2*p2p4**2*p3+4*p1*p1p2*p2p4**2*p3+4*p1p2**2*p2p4*
     -   *2*p3+2*p1*p1p3*p2p4**2*p3+4*p1p2*p1p3*p2p4**2*p3+p1*p2p4**2
     -   *p3**2+p1p4**2*p2*p3*(p2+2*p2p3+p3)+2*p1*p1p3*p2*p2p4*p3p4+4
     -   *p1p2*p1p3*p2*p2p4*p3p4+4*p1p3**2*p2*p2p4*p3p4-2*p1**2*p2p3*
     -   p2p4*p3p4-8*p1*p1p2*p2p3*p2p4*p3p4-8*p1p2**2*p2p3*p2p4*p3p4-
     -   4*p1*p1p3*p2p3*p2p4*p3p4-8*p1p2*p1p3*p2p3*p2p4*p3p4+2*p1*p1p
     -   2*p2p4*p3*p3p4+4*p1p2**2*p2p4*p3*p3p4+4*p1p2*p1p3*p2p4*p3*p3
     -   p4-2*p1*p2*p2p4*p3*p3p4-4*p1*p2p3*p2p4*p3*p3p4+p1**2*p2*p3p4
     -   **2+2*p1*p1p2*p2*p3p4**2+4*p1*p1p3*p2*p3p4**2+4*p1p2*p1p3*p2
     -   *p3p4**2+4*p1p3**2*p2*p3p4**2+p1*p2**2*p3p4**2-4*p1*p1p2*p2p
     -   3*p3p4**2-8*p1p2**2*p2p3*p3p4**2-8*p1p2*p1p3*p2p3*p3p4**2+4*
     -   p1*p2*p2p3*p3p4**2+4*p1*p2p3**2*p3p4**2+(p1**2*(p2p3**2-p2*p
     -   3)+p1p3**2*p2*(p2+2*p2p4+4*p3p4)+p1*(-(p3*(p2**2+2*p2p3*p2p4
     -   +p2*(2*(p2p3+p2p4)+p3)))+2*p1p2*(2*p2p3**2+(-p2+p2p4)*p3+p2p
     -   3*(p3-p3p4))+2*(p2p3*(p2+2*p2p3)-p2*p3)*p3p4)+p1p2**2*(4*p2p
     -   3**2+4*p2p3*(p3-p3p4)+p3*(4*p2p4+p3+2*p3p4))-2*p1p3*(p1*(p2p
     -   3*p2p4+p2*(p2p3+p3-p3p4))+p1p2*(2*p2p3*p2p4-p2p4*p3+p2*(2*p2
     -   p3+p3-p3p4)+4*p2p3*p3p4)))*p4+(p1p3**2*p2-2*p1p2*p1p3*p2p3+p
     -   1p2**2*p3+p1*(p2p3**2-p2*p3))*p4**2+4*msq*(p1p3**2*p2p4**2-p
     -   1*p2p4**2*p3+p1p4**2*(p2p3**2-p2*p3)-2*p1p2*p1p3*p2p4*p3p4+2
     -   *p1*p2p3*p2p4*p3p4+p1p2**2*p3p4**2-p1*p2*p3p4**2-2*p1p4*(p1p
     -   3*p2p3*p2p4-p1p2*p2p4*p3-p1p3*p2*p3p4+p1p2*p2p3*p3p4)-(p1p3*
     -   *2*p2-2*p1p2*p1p3*p2p3+p1*p2p3**2+p1p2**2*p3-p1*p2*p3)*p4)-2
     -   *p1p4*(p1p3*p2*(p2p4*p3+p2*p3p4+(p2p3+p3)*(2*p3p4+p4))+p1*(p
     -   2*(-(p2p3*p3p4)+p3*(p2p4+p3p4+p4))+p2p3*(p2p4*p3-p2p3*(2*p3p
     -   4+p4)))+p1p2*(p2*(-2*p2p3*p3p4+p3*(2*p2p4+p3p4+p4))-(2*p2p3+
     -   p3)*(-(p2p4*p3)+p2p3*(2*p3p4+p4)))))

      ItempX5=1d0/tempX5


       tempX4(0)=-8*(p1+2*p1p2+p1p3)*p2p4**2*p3-8*p2p4*(p1p3*(p2-2*p2p
     -   3)-2*(p1+2*p1p2)*p2p3+p1p2*p3)*p3p4-8*((p1+p1p2+2*p1p3)*p2-2
     -   *p1p2*p2p3)*p3p4**2+8*(p1*(-p2p3**2+p2*p3)+p1p3*(p2p3*p2p4+p
     -   2*(p2p3+p3-p3p4))+p1p2*(-2*p2p3**2+p2*p3-p2p3*p3-p2p4*p3+p2p
     -   3*p3p4))*p4+8*p1p4*(p2*(-(p2p3*p3p4)+p3*(p2p4+p3p4+p4))+p2p3
     -   *(p2p4*p3-p2p3*(2*p3p4+p4)))
       tempX4(1)=-8*(p1p4**2*(p2+p2p3)*p3-p1*p2p4**2*p3-2*p1p2*p2p4**2
     -   *p3-p1p3*p2p4**2*p3+p1*p1p3*p2p4*p3p4+2*p1p2*p1p3*p2p4*p3p4+
     -   2*p1p3**2*p2p4*p3p4-p1p3*p2*p2p4*p3p4+2*p1*p2p3*p2p4*p3p4+4*
     -   p1p2*p2p3*p2p4*p3p4+2*p1p3*p2p3*p2p4*p3p4-p1*p2p4*p3*p3p4-p1
     -   p2*p2p4*p3*p3p4-p1*p1p2*p3p4**2-2*p1p2**2*p3p4**2-2*p1p2*p1p
     -   3*p3p4**2-p1p2*p2*p3p4**2-2*p1p3*p2*p3p4**2+2*p1*p2p3*p3p4**
     -   2+2*p1p2*p2p3*p3p4**2+(p1p3**2*(p2+p2p4)-p1*(p2p3**2+(-p1p2+
     -   p2p4)*p3+p2p3*(p3-p3p4))+p1p3*(-(p1*p2p3)-2*p1p2*p2p3+p2*p2p
     -   3+p2p3*p2p4+p1p2*p3+p2*p3-(p1p2+p2)*p3p4)+p1p2*(-2*p2p3**2+(
     -   2*p1p2+p2-p2p4)*p3+p2p3*(-p3+p3p4)))*p4+p1p4*(-(p1*p2p4*p3)-
     -   2*p1p2*p2p4*p3-p1p3*p2p4*p3+p2*p2p4*p3+p2p3*p2p4*p3-2*p1p3*p
     -   2*p3p4+p1*p2p3*p3p4+2*p1p2*p2p3*p3p4-2*p1p3*p2p3*p3p4-p2*p2p
     -   3*p3p4-2*p2p3**2*p3p4+p1p2*p3*p3p4+p2*p3*p3p4+(-(p2p3*(p1p3+
     -   p2p3))+(p1p2+p2)*p3)*p4))
       tempX4(2)=8*(p1p4**2*p2p3*(p2+2*p2p3+p3)+p1*(p2p4+p3p4)*(p1p3*p
     -   2p4-p2p4*p3+(-p1p2+p2+2*p2p3)*p3p4)-p1*(p2*p2p3+2*p2p3**2+p1
     -   p3*(p2+p2p3)+p2p3*p2p4+p2p3*p3+p2p4*p3-p1p2*(p2p3+p3)-(p2+p2
     -   p3)*p3p4)*p4+(p1p2+p1p3)*(2*p1p3*p2p4*(p2p4+p3p4)-2*p1p2*p3p
     -   4*(p2p4+p3p4)+p1p3*(-p2+p2p4)*p4+p1p2*(2*p2p3+p3-p3p4)*p4)-p
     -   1p4*(p1*p2p4*(p2p3+p3)-p1*(p2+p2p3)*p3p4+p1p2*(2*p2p3*p2p4-(
     -   p2+p3)*p3p4-(p2p3+p3)*p4)+p1p3*(p2p4*p3+p2*(p2p4+p4)+p2p3*(4
     -   *p2p4+2*p3p4+p4))))
       tempX4(3)=-8*(-(p1*p1p2*p2p4*p3)-2*p1p2**2*p2p4*p3+p1*p2*p2p4*p
     -   3+p1*p2p3*p2p4*p3+p1*p2p4**2*p3+p1p4**2*p2*(p2p3+p3)+p1*p1p2
     -   *p2p3*p3p4+2*p1p2**2*p2p3*p3p4-p1*p2*p2p3*p3p4-2*p1*p2p3**2*
     -   p3p4-p1*p1p2*p2p4*p3p4-2*p1p2**2*p2p4*p3p4+p1*p2*p2p4*p3p4-2
     -   *p1*p2p3*p2p4*p3p4-p1p2**2*p3*p3p4+p1*p2*p3*p3p4-2*p1p2**2*p
     -   3p4**2+2*p1*p2*p3p4**2+(-(p1*p2p3*(-p1p2+p2+p2p3+p2p4))+p1p2
     -   **2*(2*p2p3-p3p4)+p1*p2*p3p4)*p4-p1p3**2*p2*(p2p4+2*p3p4+p4)
     -   +p1p4*(p1*(-(p2p3*(p2p3+p2p4))+p2*(p3+p3p4))+p1p2*(-2*p2p3**
     -   2-2*p2p4*p3+p2*(p3+p3p4)+p2p3*(-2*p2p4-p3+2*p3p4+p4)))+p1p3*
     -   (p1p4*p2*(p2p3-p2p4+p3-2*p3p4-p4)+p1*(p2p4*(p2p3+p2p4)-p2*(p
     -   3p4+p4))+p1p2*(-(p2*(p3p4+p4))+2*p2p3*(p2p4+2*p3p4+p4)+p2p4*
     -   (2*p2p4-p3+2*p3p4+p4))))
       tempX4(4)=8*(-(p1p3**2*p2*(p2p4+2*p3p4+p4))+p1*(-(p1p4*p2p3**2)
     -   +p1p4*p2*p3-p1p2*p2p4*p3+p2*p2p4*p3+p2p3*p2p4*p3+p1p2*p2p3*p
     -   3p4-p2*p2p3*p3p4-2*p2p3**2*p3p4+p2*p3*p3p4-p2p3**2*p4+p2*p3*
     -   p4)+p1p2*(p1p4*(p2*p3-p2p3*(2*p2p3+p3))-p1p2*(-2*p2p3*p3p4+p
     -   3*(2*p2p4+p3p4+p4)))+p1p3*(p1p4*p2*(p2p3+p3)+p1*(p2p3*p2p4-p
     -   2*p3p4)+p1p2*(-(p2p4*p3)-p2*p3p4+2*p2p3*(p2p4+2*p3p4+p4))))










  

       E01M=(D0_2345*tempX4(0)+D0_1345*tempX4(1)+D0_1245*tempX4(2)
     - +D0_1235*tempX4(3)+D0_1234*tempX4(4))*ItempX5    



      end

    



c$$$c --------------------------------------
c$$$c   B0finG1=B0finG(0,0,s,musq)                                 
c$$$c --------------------------------------
c$$$      complex*16 function B0finG1(s,musq)
c$$$      implicit none
c$$$      real*8 pi
c$$$      complex*16 Ipi, Ieps 
c$$$      parameter(pi=3.14159265358979324d0,
c$$$     & Ipi=(0d0,3.14159265358979324d0),Ieps=(0d0,1d-38))
c$$$      double precision s,musq,b,c
c$$$
c$$$      if (s.gt.0.d0) then 
c$$$         B0finG1 = -(dLog(s/musq)-Ipi)+2.d0
c$$$      else
c$$$         B0finG1 = -dLog(-s/musq)+2.d0
c$$$      endif
c$$$      end
c$$$
c$$$c --------------------------------------
c$$$c   B0finG2=B0finG(M1,0,s,musq)                                 
c$$$c --------------------------------------
c$$$      complex*16 function B0finG2(M1,s,musq)
c$$$      implicit none
c$$$      real*8 pi
c$$$      complex*16 Ipi, Ieps 
c$$$      parameter(pi=3.14159265358979324d0,
c$$$     & Ipi=(0d0,3.14159265358979324d0),Ieps=(0d0,1d-38))
c$$$      double precision M1,s,musq,msq
c$$$
c$$$      msq = M1*M1
c$$$      if (s.gt.msq) then 
c$$$         B0finG2 = -dLog(msq/musq)
c$$$     &          +(msq-s)/s*(dLog((s-msq)/msq)-Ipi)+2.d0
c$$$      else
c$$$         B0finG2 = -dLog(msq/musq)
c$$$     &          +(msq-s)/s*dLog(-(s-msq)/msq)+2.d0
c$$$      endif 
c$$$      end
c$$$c
c--------------  B0tM(m,qsq): regularized 2-point function --------------
c
      complex*16 function B0tM1(m,qsq) 
      implicit none
      double precision m, qsq, qsqn

c evaluate scalar 2-point function for equal masses m on propagators 
c  
c    B0 = Int d^4k [k^2-m^2]^-1 [(k-q)^2-m^2]^-1 
c
c Subtracting the divergent piece, 1/eps - gamma + log(4pi mu^2/m^2),
c one obtains the modified scalar 2-point function B_0~ which is evaluated 
c here
c
c   B0tM(m,q^2) = - int_0^1 dx log[ 1 - q^2/(m^2-i eps) x(1-x) ]
c
c	Dieter Zeppenfeld, <dieter@pheno.physics.wisc.edu>
c	Initial version:  2000 April 7
c	Last modified:    2000 November 12
c  

      double precision phi, beta, srt, lnfac, re, im
      double precision eps, pi
      parameter (pi=3.14159 26535 89793d0)
      parameter (eps=5d-4)  ! limit of q^2/m^2 << 1 approximation 

      qsqn = qsq/m**2
      if ( qsqn.lt.-eps ) then
         srt = sqrt(1d0-4d0/qsqn)
         lnfac = log( (srt-1d0)/(srt+1d0) )
         B0tM1 = 2d0 + srt*lnfac
      elseif ( abs(qsqn).le.eps )then
         B0tM1 = qsqn/6d0* ( 1d0+qsqn*0.1d0*(1d0+qsqn/7d0 *
     &        ( 1d0+qsqn/7d0*( 1d0+2d0/11d0*qsqn ) ) ) )
      elseif (qsqn.lt.4d0) then
         srt = sqrt(4d0/qsqn-1d0)
         phi = atan(1d0/srt)
         B0tM1 = 2d0 - 2d0*srt*phi
      elseif (qsqn.eq.4d0) then
         B0tM1 = 2d0
      else
         beta = sqrt(1d0-4d0/qsqn)
         lnfac = log( (1d0-beta)/(1d0+beta) )
         re = 2d0 + beta*lnfac
         im = pi*beta
         B0tM1 = cmplx( re, im )
      endif
      return
      end








c ------------------------------------------------------------------------
c ---- complex dilogarithm -----------------------------------------------
c ------------------------------------------------------------------------
       complex*16 function cdilog(z)                                    
       implicit none                                   
       complex*16 z,zl,coef,dilog1,u,caux                               
       real*8 pi,sign                                                   
       integer n,i 
       pi=3.141592653589793238462643d0                                  
       zl=z                                                             
       dilog1=dcmplx(pi**2/6.d0)                                        
       if(dreal(zl).eq.1.and.dimag(zl).eq.0.) then                      
          cdilog=dilog1                                                    
          return                                                           
       else if (cdabs(zl).lt.1.d-2) then   
          n=-40./dlog(cdabs(zl))                                           
          caux=(0.d0,0.d0)                                                 
          do i=1,n                                                         
             caux=caux+zl**i/dble(i**2)                                       
          enddo                                                           
          cdilog=caux                                                      
          return                                                           
       else if(cdabs(zl).lt.1.) then                                    
          sign=1.d0                                                       
          coef=dcmplx(dble(0.))                                           
       else                                                            
          coef=-cdlog(-zl)**2/2.d0-dilog1                                 
          sign=-1.d0                                                      
          zl=1.d0/zl                                                      
       endif                                                          
       if(dreal(zl).gt.0.5) then                   
          coef=coef+sign*(dilog1-cdlog(zl)*cdlog(1.d0-zl))                
          sign=-sign                                                      
          zl=1.d0-zl                                                      
       else   
       endif  
       u=-cdlog(1.d0-zl)                                               
       cdilog=u-u**2/4.d0+u**3/36.d0-u**5/3600.d0+u**7/211680.d0       
     &  -u**9/10886400.d0+u**11*5.d0/2634508800.d0                     
       cdilog=cdilog-u**13*691.d0/2730.d0/6227020800.d0                
       cdilog=cdilog+u**15*7.d0/6.d0/1.307674368d12                    
       cdilog=cdilog-u**17*3617.d0/510.d0/3.5568742810d14              
       cdilog=cdilog+u**19*43867.d0/798.d0/1.2164510041d17              
       cdilog=cdilog-u**21*174611.d0/330.d0/5.1090942172d19            
       cdilog=sign*cdilog+coef                                         
       return                                                          
       end                                                             
