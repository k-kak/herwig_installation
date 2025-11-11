      DOUBLE PRECISION FUNCTION DILOG(X)                                       
                                                                                
      DOUBLE PRECISION X,Y,T,S,A,PI3,PI6,ZERO,ONE,HALF,MALF,MONE,MTWO           
      DOUBLE PRECISION C(0:18),H,ALFA,B0,B1,B2,B3,B4,B5,B6,B7,B8,B9,B10
      DOUBLE PRECISION B11,B12,B13,B14,B15,B16,B17,B18                                  
                                                                                
      DATA ZERO /0.0D0/, ONE /1.0D0/                                            
      DATA HALF /0.5D0/, MALF /-0.5D0/, MONE /-1.0D0/, MTWO /-2.0D0/            
      DATA PI3 /3.28986 81336 96453D0/, PI6 /1.64493 40668 48226D0/             
                                                                                
      DATA C( 0) / 0.42996 69356 08136 97D0/                                     
      DATA C( 1) / 0.40975 98753 30771 05D0/                                     
      DATA C( 2) /-0.01858 84366 50145 92D0/                                     
      DATA C( 3) / 0.00145 75108 40622 68D0/                                     
      DATA C( 4) /-0.00014 30418 44423 40D0/                                     
      DATA C( 5) / 0.00001 58841 55418 80D0/                                     
      DATA C( 6) /-0.00000 19078 49593 87D0/                                     
      DATA C( 7) / 0.00000 02419 51808 54D0/                                     
      DATA C( 8) /-0.00000 00319 33412 74D0/                                     
      DATA C( 9) / 0.00000 00043 45450 63D0/                                     
      DATA C(10) /-0.00000 00006 05784 80D0/                                     
      DATA C(11) / 0.00000 00000 86120 98D0/                                     
      DATA C(12) /-0.00000 00000 12443 32D0/                                     
      DATA C(13) / 0.00000 00000 01822 56D0/                                     
      DATA C(14) /-0.00000 00000 00270 07D0/                                     
      DATA C(15) / 0.00000 00000 00040 42D0/                                     
      DATA C(16) /-0.00000 00000 00006 10D0/                                     
      DATA C(17) / 0.00000 00000 00000 93D0/                                     
      DATA C(18) /-0.00000 00000 00000 14D0/                                     
c      DATA C(19) /-0.00000 00000 00000 02D0/                                                                               
      IF(X .EQ. ONE) THEN                                                       
       DILOG=PI6                                                               
       RETURN                                                                   
      ELSE IF(X .EQ. MONE) THEN                                                 
       DILOG=MALF*PI6                                                          
       RETURN                                                                   
      END IF                                                                    
      T=-X                                                                      
      IF(T .LE. MTWO) THEN                                                      
       Y=MONE/(ONE+T)                                                           
       S=ONE                                                                    
       A=-PI3+HALF*(LOG(-T)**2-LOG(ONE+ONE/T)**2)                               
      ELSE IF(T .LT. MONE) THEN                                                 
       Y=MONE-T                                                                 
       S=MONE                                                                   
       A=LOG(-T)                                                                
       A=-PI6+A*(A+LOG(ONE+ONE/T))                                              
      ELSE IF(T .LE. MALF) THEN                                                 
       Y=(MONE-T)/T                                                             
       S=ONE                                                                    
       A=LOG(-T)                                                                
       A=-PI6+A*(MALF*A+LOG(ONE+T))                                             
      ELSE IF(T .LT. ZERO) THEN                                                 
       Y=-T/(ONE+T)                                                             
       S=MONE                                                                   
       A=HALF*LOG(ONE+T)**2                                                     
      ELSE IF(T .LE. ONE) THEN                                                  
       Y=T                                                                      
       S=ONE                                                                    
       A=ZERO                                                                   
      ELSE                                                                      
       Y=ONE/T                                                                  
       S=MONE                                                                   
       A=PI6+HALF*LOG(T)**2                                                     
      END IF                                                                    
                                                                                
      H=Y+Y-ONE                                                                 
      ALFA=H+H                                                                                                                                           
      B18=C(18)
      B17=C(17)+ALFA*B18
      B16=C(16)+ALFA*B17-B18
      B15=C(15)+ALFA*B16-B17
      B14=C(14)+ALFA*B15-B16
      B13=C(13)+ALFA*B14-B15 
      B12=C(12)+ALFA*B13-B14 
      B11=C(11)+ALFA*B12-B13 
      B10=C(10)+ALFA*B11-B12 
      B9=C(9)  +ALFA*B10-B11
      B8=C(8)  +ALFA*B9-B10
      B7=C(7)  +ALFA*B8-B9 
      B6=C(6)  +ALFA*B7-B8 
      B5=C(5)  +ALFA*B6-B7 
      B4=C(4)  +ALFA*B5-B6 
      B3=C(3)  +ALFA*B4-B5 
      B2=C(2)  +ALFA*B3-B4 
      B1=C(1)  +ALFA*B2-B3 
      B0=C(0)  +ALFA*B1-B2 


      DILOG=-(S*(B0-H*B1)+A)    
C      write(*,'(I2,E28.19)') I,DILOG
C      enddo
                                               
      RETURN                                                                    
      END                                                                       
