    
c
c-------------  spence function  li2(z) ------------------------
c
      complex*16 function li2(zin)
      implicit none
      complex*16 zin, z, u, u2, unpo, ans, zext
      double precision r, r2, r2n, fac
c
c determine the value of the dilogarithm 
c
c    li2(z) = - int_0^1  log(1-zt)/t dt  with cut along the positive 
c                                        real axis, z>1
c
c	Dieter Zeppenfeld, <dieter@pheno.physics.wisc.edu>
c	Initial version:  2000 November 6
c	Last modified:    2000 November 12
c  
      integer i
      double precision c0,c1,c2,c4,c6,c8,c10,c12,c14,c16,c18,c20,c22
      double precision b0,b1,b2,b4,b6,b8,b10,b12,b14,b16,b18,b20,b22
      double precision d0,d1,d2,d4,d6,d8,d10,d12,d14,d16,d18,d20,d22
      parameter (b0=1d0,            d0 =1d0,      c0= b0/d0)
      parameter (b1=-1d0/2d0,       d1 =d0*2d0,   c1= b1/d1)
      parameter (b2= 1d0/6d0,       d2 =d1*3d0,   c2= b2/d2)
      parameter (b4=-1d0/30d0,      d4 =d2*20d0,  c4= b4/d4)
      parameter (b6=1d0/42d0,       d6 =d4*42d0,  c6= b6/d6)
      parameter (b8=-1d0/30d0,      d8 =d6*72d0,  c8= b8/d8)
      parameter (b10=5d0/66d0,      d10=d8*110d0, c10=b10/d10)
      parameter (b12=-691d0/2730d0, d12=d10*156d0,c12=b12/d12)
      parameter (b14=7d0/6d0,       d14=d12*210d0,c14=b14/d14)
      parameter (b16=-3617d0/510d0, d16=d14*272d0,c16=b16/d16)
      parameter (b18=43867d0/798d0, d18=d16*342d0,c18=b18/d18)
      parameter (b20=-174611d0/330d0,d20=d18*420d0,c20=b20/d20)
      parameter (b22=854513d0/138d0,d22=d20*506d0,c22=b22/d22)
      double precision eps, epst, pi, pi2o6
      parameter (eps=1d-16, epst=1d-3)
      parameter (pi=3.14159 26535 89793, pi2o6=pi**2/6)
c
c debug information
      logical ldebug
      parameter (ldebug=.false.)

      z = zin
c      print*,' li2 call with z = ',z
      u = z**2
      r2 = dreal(z)**2+dimag(z)**2 
      if (r2.lt.eps) then
         li2 = z + u/4d0
         return
      elseif (r2.lt.epst) then
         ans = z + u/4
         do i = 3,11
            u = u*z
            ans = ans + u/i**2
         enddo
         li2 = ans
         return
      endif
      if (dreal(z).ge.1d0 .and. dimag(z).eq.0 ) then
         z = z + (0d0,1d0)*eps
      endif
c
c use z-->1/z and z--> 1-z mappings of the spence function to restrict 
c agument to unit circle in the complex plane with Re(z) <= 0.5
c
      zext = (0d0,0d0)
      fac = 1
      if (r2.gt.1d0) then     ! map z ---> 1/z
         fac = -fac
         zext = -pi2o6 - 0.5d0*(log(-z))**2
         z = 1/z
      endif
      if (dreal(z).gt.0.5d0) then     ! map new z ---> 1-z
         zext = zext + fac*(pi2o6-log(z)*log(1-z))
         fac = -fac
         z = 1-z
      endif
c
c now use t = 1 - exp(-u) mapping to write Li(z) in terms of Bernoulli 
c numbers
c
      u = - log(1-z)
      r2 = abs(u)**2
      u2 = u*u
      ans = u*(c0 + u*(c1+c2*u))
      r2n = r2*r2       !r^4

      unpo = u2*u2*u
      ans = ans + c4*unpo

      unpo = unpo*u2
      ans = ans + c6*unpo

      r = r2n*r2n       !r^8
      unpo = unpo*u2
      ans = ans + c8*unpo

      r2n = r*r2        !r^10 
      if ((r2n*c10).gt.eps) then
         unpo = unpo*u2
         ans = ans + c10*unpo
      else
         li2 = fac * ans + zext
         if (ldebug) print*,' exit li2s at n=8 '
         return
      endif

      unpo = unpo*u2
      ans = ans + c12*unpo

      unpo = unpo*u2
      ans = ans + c14*unpo

      unpo = unpo*u2
      ans = ans + c16*unpo

      r2n = r2n*r
      if ((r2n*c18).gt.eps) then
         unpo = unpo*u2
         ans = ans + c18*unpo
      else
         li2 = fac * ans + zext
         if (ldebug) print*,' exit li2s at n=16 '
         return
      endif

      unpo = unpo*u2
      ans = ans + c20*unpo

      unpo = unpo*u2
      ans = ans + c22*unpo

      li2 = fac * ans + zext
      end










c$$$
c$$$  DOUBLE PRECISION FUNCTION li2(X)                                       
c$$$                                                                                
c$$$      DOUBLE PRECISION X,Y,T,S,A,PI3,PI6,ZERO,ONE,HALF,MALF,MONE,MTWO           
c$$$      DOUBLE PRECISION C(0:18),H,ALFA,B0,B1,B2,B3,B4,B5,B6,B7,B8,B9,B10
c$$$      DOUBLE PRECISION B11,B12,B13,B14,B15,B16,B17,B18                                  
c$$$                                                                                
c$$$      DATA ZERO /0.0D0/, ONE /1.0D0/                                            
c$$$      DATA HALF /0.5D0/, MALF /-0.5D0/, MONE /-1.0D0/, MTWO /-2.0D0/            
c$$$      DATA PI3 /3.28986 81336 96453D0/, PI6 /1.64493 40668 48226D0/             
c$$$                                                                                
c$$$      DATA C( 0) / 0.42996 69356 08136 97D0/                                     
c$$$      DATA C( 1) / 0.40975 98753 30771 05D0/                                     
c$$$      DATA C( 2) /-0.01858 84366 50145 92D0/                                     
c$$$      DATA C( 3) / 0.00145 75108 40622 68D0/                                     
c$$$      DATA C( 4) /-0.00014 30418 44423 40D0/                                     
c$$$      DATA C( 5) / 0.00001 58841 55418 80D0/                                     
c$$$      DATA C( 6) /-0.00000 19078 49593 87D0/                                     
c$$$      DATA C( 7) / 0.00000 02419 51808 54D0/                                     
c$$$      DATA C( 8) /-0.00000 00319 33412 74D0/                                     
c$$$      DATA C( 9) / 0.00000 00043 45450 63D0/                                     
c$$$      DATA C(10) /-0.00000 00006 05784 80D0/                                     
c$$$      DATA C(11) / 0.00000 00000 86120 98D0/                                     
c$$$      DATA C(12) /-0.00000 00000 12443 32D0/                                     
c$$$      DATA C(13) / 0.00000 00000 01822 56D0/                                     
c$$$      DATA C(14) /-0.00000 00000 00270 07D0/                                     
c$$$      DATA C(15) / 0.00000 00000 00040 42D0/                                     
c$$$      DATA C(16) /-0.00000 00000 00006 10D0/                                     
c$$$      DATA C(17) / 0.00000 00000 00000 93D0/                                     
c$$$      DATA C(18) /-0.00000 00000 00000 14D0/                                     
c$$$c      DATA C(19) /-0.00000 00000 00000 02D0/                                                                               
c$$$      IF(X .EQ. ONE) THEN                                                       
c$$$       DILOG=PI6                                                               
c$$$       RETURN                                                                   
c$$$      ELSE IF(X .EQ. MONE) THEN                                                 
c$$$       DILOG=MALF*PI6                                                          
c$$$       RETURN                                                                   
c$$$      END IF                                                                    
c$$$      T=-X                                                                      
c$$$      IF(T .LE. MTWO) THEN                                                      
c$$$       Y=MONE/(ONE+T)                                                           
c$$$       S=ONE                                                                    
c$$$       A=-PI3+HALF*(LOG(-T)**2-LOG(ONE+ONE/T)**2)                               
c$$$      ELSE IF(T .LT. MONE) THEN                                                 
c$$$       Y=MONE-T                                                                 
c$$$       S=MONE                                                                   
c$$$       A=LOG(-T)                                                                
c$$$       A=-PI6+A*(A+LOG(ONE+ONE/T))                                              
c$$$      ELSE IF(T .LE. MALF) THEN                                                 
c$$$       Y=(MONE-T)/T                                                             
c$$$       S=ONE                                                                    
c$$$       A=LOG(-T)                                                                
c$$$       A=-PI6+A*(MALF*A+LOG(ONE+T))                                             
c$$$      ELSE IF(T .LT. ZERO) THEN                                                 
c$$$       Y=-T/(ONE+T)                                                             
c$$$       S=MONE                                                                   
c$$$       A=HALF*LOG(ONE+T)**2                                                     
c$$$      ELSE IF(T .LE. ONE) THEN                                                  
c$$$       Y=T                                                                      
c$$$       S=ONE                                                                    
c$$$       A=ZERO                                                                   
c$$$      ELSE                                                                      
c$$$       Y=ONE/T                                                                  
c$$$       S=MONE                                                                   
c$$$       A=PI6+HALF*LOG(T)**2                                                     
c$$$      END IF                                                                    
c$$$                                                                                
c$$$      H=Y+Y-ONE                                                                 
c$$$      ALFA=H+H                                                                                                                                           
c$$$      B18=C(18)
c$$$      B17=C(17)+ALFA*B18
c$$$      B16=C(16)+ALFA*B17-B18
c$$$      B15=C(15)+ALFA*B16-B17
c$$$      B14=C(14)+ALFA*B15-B16
c$$$      B13=C(13)+ALFA*B14-B15 
c$$$      B12=C(12)+ALFA*B13-B14 
c$$$      B11=C(11)+ALFA*B12-B13 
c$$$      B10=C(10)+ALFA*B11-B12 
c$$$      B9=C(9)  +ALFA*B10-B11
c$$$      B8=C(8)  +ALFA*B9-B10
c$$$      B7=C(7)  +ALFA*B8-B9 
c$$$      B6=C(6)  +ALFA*B7-B8 
c$$$      B5=C(5)  +ALFA*B6-B7 
c$$$      B4=C(4)  +ALFA*B5-B6 
c$$$      B3=C(3)  +ALFA*B4-B5 
c$$$      B2=C(2)  +ALFA*B3-B4 
c$$$      B1=C(1)  +ALFA*B2-B3 
c$$$      B0=C(0)  +ALFA*B1-B2 
c$$$
c$$$
c$$$      DILOG=-(S*(B0-H*B1)+A)    
c$$$C      write(*,'(I2,E28.19)') I,DILOG
c$$$C      enddo
c$$$                                               
c$$$      RETURN                                                                    
c$$$      END                                                                       
