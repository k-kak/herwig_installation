        subroutine tens_red4_Complex_G(m0sq,m1sq,m2sq,m3sq,p1sq,p2sq,p3sq,p1p2,p1p3,p2p3, 
     &                     C0rr_234, C0rr_134, C0rr_124, C0rr_123,   
     &                     Cijrr_234, Cijrr_134, Cijrr_124, Cijrr_123,             
     &                     C0I_234, C0I_134, C0I_124, C0I_123,   
     &                     CijI_234, CijI_134, CijI_124, CijI_123,   
     &                     D0, D0rr, D0I,Dijrr,DijI)
C                tens_red3 = 3-point tensors
C     Francisco Campanario
C    Email: Francam@particle.uni-karlsruhe.de
C Date: 25/02/2010
c
c
c  determine the Passarino-Veltman tensor decomposition for the four-point
c  tensor integrals
c
c                                          d^4k
c   D0; D_mu; D_mu,nu; D_mu,nu,rho =  Int ------ 
c                                         (2pi)^4
c
c              1;  k_mu;   k_mu k_nu; k_mu k_nu k_rho
c      -------------------------------------------------------------------
c         [k^2-m^2][(k+p1)^2-m^2][(k+p1+p2)^2-m^2][(k+p1+p2+p3)^2-m^2]
c with
c
c   D_mu = p1_mu D11  +  p2_mu D12  +  p3_mu D13
c
c   D_munu = p1_mu p2_nu D21 + p2_mu p2_nu D22 + ...
c
c  for notation see Passarino&Veltman, NP B160 (1979) 151 
c
C INPUT:  psq, pq,...                        kinematics invariants
C         C0_123 = C0(1,2,3) = C0(p1,p2)     scalar three point 
C         C0_124 = C0(1,2,4) = C0(p1,p2+p3)  functions in PV notation
C         C0_134 = C0(1,3,4) = C0(p1+p2,p3)
C         C0_234 = C0(2,3,4) = C0(p2,p3)
C         Cij_123(n,m) = C_nm(1,2,3) ....    higher C_nm form factors
C                                            as in tens_red3
c         D0 = D0(p,q,l)                  scalar four point function
c
c OUTPUT: Dij(n,m) = D_nm                    form factors in the tensor 
c                                            integrals a la PV
c         nm = 11, 21, 31                    ff"s for D_mu
c         nm = 12, 22, 32, 42, 52, 62, 72    ff"s for D_munu
c         nm = 13, 23, 33, ..., 93, 103, 113, 123  ff"s for D_mu,nu,rho
c
        implicit none
        real * 8 p1sq, p2sq, p3sq, p1p2, p1p3,p2p3
        complex*16 D0
        real*8  det,Iv,Ivv
        real*8 z11,z12,z13,z21,z22,z23,z31,z32,z33
        real*8 iz11,iz22,iz33,piv
        real*8 C0rr_234, C0rr_134, C0rr_124, C0rr_123
        real*8 Cijrr_234(4,2), Cijrr_134(4,2),Cijrr_124(4,2), Cijrr_123(4,2)
        real*8 C0I_234, C0I_134, C0I_124, C0I_123
        real*8 CijI_234(4,2), CijI_134(4,2), CijI_124(4,2), CijI_123(4,2)
        real*8 D0rr, Dijrr(13,3)
        real*8 D0I, DijI(13,3)
c New variables to compare       
        Integer indexC(2),i,j,indexD(3)
        real*8 s23
c Added for the General case
        complex*16 m0sq,m1sq,m2sq,m3sq,r1r0,r2r1, r3r2!,m0,m1,m2,m3
c
        complex*16 PRr(1:3)
        complex*16 Rr(1:3)
        complex*16 C0r_234, C0r_134, C0r_124, C0r_123
        complex*16 Cijr_234(4,2),Cijr_134(4,2),Cijr_124(4,2), Cijr_123(4,2)
        complex*16 D0r, Dijr(13,3)

        indexC(1)=2
        indexC(2)=4

        indexD(1)=13
        indexD(2)=7
        indexD(2)=13

c        m0sq=m0*m0
c        m1sq=m1*m1
c        m2sq=m2*m2
c        m3sq=m3*m3
  

        s23=p2sq+p3sq+2d0*p2p3
  
        r1r0 = p1sq-m1sq+m0sq
        r2r1 = p2sq+2.d0*p1p2-m2sq+m1sq
        r3r2 = 2.d0*(p1p3+p2p3)+p3sq-m3sq+m2sq
    
       det =-2.d0*(-2.d0*p1p2*p1p3*p2p3 +p1p3*p1p3*p2sq + p1p2*p1p2*p3sq 
     -  + p1sq*(p2p3*p2p3 - p2sq*p3sq))


       D0rr=Dble(D0)
       D0I=DImag(D0)
  
       D0r=D0

       do j=1,2
        do i=1,indexC(j)
          Cijr_123(i,j)=DCMPLX(Cijrr_123(i,j),Ciji_123(i,j))
          Cijr_124(i,j)=DCMPLX(Cijrr_124(i,j),Ciji_124(i,j))
          Cijr_134(i,j)=DCMPLX(Cijrr_134(i,j),Ciji_134(i,j))
          Cijr_234(i,j)=DCMPLX(Cijrr_234(i,j),Ciji_234(i,j))
        enddo   
       enddo  
          C0r_123=DCMPLX(C0rr_123,C0i_123)
          C0r_124=DCMPLX(C0rr_124,C0i_124)
          C0r_134=DCMPLX(C0rr_134,C0i_134)
          C0r_234=DCMPLX(C0rr_234,C0i_234)

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC      
	if(abs(p1sq).gt.abs(p1p2)) then
		if(abs(p1sq).gt.abs(p1p3)) then
			z11=2d0*p1sq
			iz11=1d0/z11
			z12=2d0*p1p2
			z13=2d0*p1p3
c			
			z21=z12*iz11
			z31=z13*iz11
c	
			z22=2d0*p2sq-z21*z12
			z32=2d0*p2p3-z31*z12
c
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC	
			if(abs(z32).gt.abs(z22)) then	
c				iorder(1)=1
c				iorder(2)=3
c                                iorder(3)=2
				piv=z21
				z21=z31
				z31=piv
				piv=z22
				z22=z32
				iz22=1d0/z22
				z32=piv
c
				z23=2d0*p3sq
				z33=2d0*p2p3
c				det=-det
				z32=z32*iz22
	            z23=z23-z21*z13
	            z33=z33-z31*z13-z32*z23
	            iz33=1d0/z33
        	  det=-z11*z22*z33
c             print*, "det1",det



       PRr(1) = (C0r_134 - C0r_234 - D0r*r1r0   )
       PRr(2) = (C0r_124 - C0r_134 - D0r*r2r1 )
       PRr(3) = (C0r_123 - C0r_124 - D0r*r3r2 )


           Rr(3)=(PRr(2)-z32*(PRr(3)-z21*PRr(1))-z31*PRr(1))*iz33
           Rr(2)=(PRr(3)-z21*PRr(1) - z23*Rr(3))*iz22
           Rr(1)=(PRr(1)-  z12*Rr(2)-  z13*Rr(3))*iz11 	
       

       Dijr(1,1) = Rr(1)
       Dijr(2,1) = Rr(2)
       Dijr(3,1) = Rr(3)

CCCCCC CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
c D00  
CCCCCC CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
               
      Dijr(7,2) = (C0r_234+ 2d0*m0sq*D0r + Dijr(1,1)*r1r0 + Dijr(2,1)*r2r1
     &     +Dijr(3,1)*r3r2)*0.5d0


c FC %CCCCCC CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
c FC %CCCCCC CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
c FC %c 30-3 1-32
c FC %       
c FC %c      
c FC %       PRr(1) = (Cijr_134(1,1) + C0r_234 - Dijr(1,1)*r1r0 - Dijr(7,2)*2.d0)
c FC %       PRr(2) = (Cijr_124(1,1) - Cijr_134(1,1) - Dijr(1,1)*r2r1        )
c FC %       PRr(3) = (Cijr_123(1,1) - Cijr_124(1,1) - Dijr(1,1)*r3r2        )     
c FC %c       
c FC %           Rr(3)=(PRr(2)-z32*(PRr(3)-z21*PRr(1))-z31*PRr(1))*iz33
c FC %           Rr(2)=(PRr(3)-z21*PRr(1) - z23*Rr(3))*iz22
c FC %           Rr(1)=(PRr(1)-  z12*Rr(2)-  z13*Rr(3))*iz11 
c FC %       
c FC %       Dijr(1,2) = Rr(1)
c FC %       Dijr(4,2) = Rr(2)
c FC %       Dijr(5,2) = Rr(3)
c FC %       
c FC %c 33-3 4-35
c FC %       PRr(1) = (Cijr_134(1,1) -Cijr_234(1,1) - Dijr(2,1)*r1r0          )
c FC %       PRr(2)=(Cijr_124(2,1)-Cijr_134(1,1)-Dijr(2,1)*r2r1-Dijr(7,2)*2.d0)
c FC %       PRr(3) = (Cijr_123(2,1) -Cijr_124(2,1) - Dijr(2,1)*r3r2        )
c FC %       
c FC %       Rr(3)=(PRr(2)-z32*(PRr(3)-z21*PRr(1))-z31*PRr(1))*iz33
c FC %           Rr(2)=(PRr(3)-z21*PRr(1) - z23*Rr(3))*iz22
c FC %c          Rr(1)=(PRr(1)-  z12*Rr(2)-  z13*Rr(3))*iz11 
c FC %       
c FC %c       Dijr(4,2) = Rr(1)
c FC %       Dijr(2,2) = Rr(2)
c FC %       Dijr(6,2) = Rr(3)
c FC %       
c FC %c 36-3 7-38
c FC %       
c FC %       PRr(1) = (Cijr_134(2,1) - Cijr_234(2,1) - Dijr(3,1)*r1r0      )
c FC %       PRr(2) = (Cijr_124(2,1) - Cijr_134(2,1) - Dijr(3,1)*r2r1    )
c FC %       PRr(3) = (       - Cijr_124(2,1) - Dijr(3,1)*r3r2 -Dijr(7,2)*2.d0)
c FC %       
c FC %       Rr(3)=(PRr(2)-z32*(PRr(3)-z21*PRr(1))-z31*PRr(1))*iz33
c FC %c          Rr(2)=(PRr(3)-z21*PRr(1) - z23*Rr(3))*iz22
c FC %c          Rr(1)=(PRr(1)-  z12*Rr(2)-  z13*Rr(3))*iz11 
c FC %       
c FC %c       Dijr(5,2) = Rr(1)
c FC %c       Dijr(6,2) = Rr(2)
c FC %       Dijr(3,2) = Rr(3)
c FC %    
c FC %CCCCCC CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
c FC %CCCCCC CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
c FC %        
c FC %       Dijr(11,3)=0.25d0*(-C0r_234+2d0*m0sq*Dijr(1,1)+Dijr(1,2)*r1r0+Dijr(4,2)*r2r1
c FC %     &       +Dijr(5,2)*r3r2)
c FC %       
c FC %       Dijr(12,3)=0.25d0*(Cijr_234(1,1)+2d0*m0sq*Dijr(2,1)+Dijr(4,2)*r1r0+Dijr(2,2)*r2r1
c FC %     &                   +Dijr(6,2)*r3r2)
c FC %       Dijr(13,3)=0.25d0*(Cijr_234(2,1)+2d0*m0sq*Dijr(3,1)+Dijr(5,2)*r1r0 + Dijr(6,2)*r2r1
c FC %     &                    + Dijr(3,2)*r3r2)
c FC %       
c FC %     
c FC %c 41-4 2-43
c FC %       
c FC %       PRr(1) =(Cijr_134(1,2) - C0r_234 - Dijr(1,2)*r1r0  - Dijr(11,3)*4.d0)
c FC %       PRr(2) =(Cijr_124(1,2) - Cijr_134(1,2) - Dijr(1,2)*r2r1       )
c FC %       PRr(3) =(Cijr_123(1,2) - Cijr_124(1,2) - Dijr(1,2)*r3r2       )
c FC %       
c FC %       Rr(3)=(PRr(2)-z32*(PRr(3)-z21*PRr(1))-z31*PRr(1))*iz33
c FC %           Rr(2)=(PRr(3)-z21*PRr(1) - z23*Rr(3))*iz22
c FC %           Rr(1)=(PRr(1)-  z12*Rr(2)-  z13*Rr(3))*iz11 
c FC %       
c FC %       Dijr(1,3) = Rr(1)
c FC %       Dijr(4,3) = Rr(2)
c FC %       Dijr(5,3) = Rr(3)
c FC %       
c FC %c 50-5 1-52
c FC %       
c FC %       PRr(1) =  (Cijr_134(1,2) - Cijr_234(1,2) - Dijr(2,2)*r1r0)
c FC %       PRr(2) =  (Cijr_124(2,2) - Cijr_134(1,2) - Dijr(2,2)*r2r1
c FC %     &       -Dijr(12,3)*4.d0)
c FC %       PRr(3) =  (Cijr_123(2,2) - Cijr_124(2,2) - Dijr(2,2)*r3r2)
c FC %       
c FC %       Rr(3)=(PRr(2)-z32*(PRr(3)-z21*PRr(1))-z31*PRr(1))*iz33
c FC %           Rr(2)=(PRr(3)-z21*PRr(1) - z23*Rr(3))*iz22
c FC %           Rr(1)=(PRr(1)-  z12*Rr(2)-  z13*Rr(3))*iz11 
c FC %       
c FC %       Dijr(6,3) = Rr(1)
c FC %       Dijr(2,3) = Rr(2)
c FC %       Dijr(8,3) = Rr(3)
c FC %       
c FC %c      
c FC %c 56-5 7-58
c FC %c      
c FC %c      
c FC %       PRr(1) = (Cijr_134(2,2) - Cijr_234(2,2) - Dijr(3,2)*r1r0)
c FC %       PRr(2) = (Cijr_124(2,2) - Cijr_134(2,2) - Dijr(3,2)*r2r1)
c FC %       PRr(3) = (     -Cijr_124(2,2) - Dijr(3,2)*r3r2 - Dijr(13,3)*4.d0)
c FC %c       
c FC %       
c FC %       Rr(3)=(PRr(2)-z32*(PRr(3)-z21*PRr(1))-z31*PRr(1))*iz33
c FC %           Rr(2)=(PRr(3)-z21*PRr(1) - z23*Rr(3))*iz22
c FC %           Rr(1)=(PRr(1)-  z12*Rr(2)-  z13*Rr(3))*iz11 
c FC %       
c FC %       Dijr(7,3) = Rr(1)
c FC %       Dijr(9,3) = Rr(2)
c FC %       Dijr(3,3) = Rr(3)
c FC %c      
c FC %c 44-4 5-46
c FC %       PRr(1) = (Cijr_134(1,2) + Cijr_234(1,1) - Dijr(4,2)*r1r0 
c FC %     &      -Dijr(12,3)*2.d0)
c FC %       PRr(2) = (Cijr_124(3,2) - Cijr_134(1,2) - Dijr(4,2)*r2r1
c FC %     &      - Dijr(11,3)*2.d0)
c FC %       PRr(3) = (Cijr_123(3,2) - Cijr_124(3,2)  -Dijr(4,2)*r3r2)
c FC %c      
c FC %       
c FC %       Rr(3)=(PRr(2)-z32*(PRr(3)-z21*PRr(1))-z31*PRr(1))*iz33
c FC %c          Rr(2)=(PRr(3)-z21*PRr(1) - z23*Rr(3))*iz22
c FC %c          Rr(1)=(PRr(1)-  z12*Rr(2)-  z13*Rr(3))*iz11 
c FC %       
c FC %C       Dijr(4,3) = Rr(1)
c FC %C       Dijr(6,3) = Rr(2)
c FC %       Dijr(10,3) = Rr(3)
c FC %c      
c FC %c      
c FC %c      
c FC %c to d efine D00ij and D0000 functions!
c FC %c In P V notation we have:
c FC %c Dij( 1,7)=D416, Dij(1,8)=D417, Dij(1,9)=D418
c FC %c Dij( 1,10)=D419, Dij(1,11)=D420, Dij(1,12)=D421, D(1,13)=D422
c FC %c      
c FC %       Iv=1.d0/6.d0
c FC %c      
c FC %c 
c FC %
c FC %      Dijr(7,1)=(C0r_234+2*m0sq*Dijr(1,2)+Dijr(1,3)*r1r0+Dijr(4,3)*r2r1+Dijr(5,3)*r3r2 )*Iv
c FC %c      
c FC %      Dijr(8,1)=(Cijr_234(1,2)+2*m0sq*Dijr(2,2)+Dijr(6,3)*r1r0+Dijr(2,3)*r2r1
c FC %     - +Dijr(8,3)*r3r2)*Iv
c FC %c
c FC %      Dijr(9,1)=(Cijr_234(2,2)+2*m0sq*Dijr(3,2)+Dijr(7,3)*r1r0+Dijr(9,3)*r2r1
c FC %     - +Dijr(3,3)*r3r2)*Iv
c FC %c
c FC %      Dijr(10,1)=(-Cijr_234(1,1)+2*m0sq*Dijr(4,2)+Dijr(4,3)*r1r0+Dijr(6,3)*r2r1
c FC %     - +Dijr(10,3)*r3r2)*Iv
c FC %c      
c FC %      Dijr(11,1)=(-Cijr_234(2,1)+2*m0sq*Dijr(5,2)+Dijr(5,3)*r1r0+Dijr(10,3)*r2r1
c FC %     - +Dijr(7,3)*r3r2)*Iv
c FC %c
c FC %      Dijr(12,1)=(Cijr_234(3,2)+2*m0sq*Dijr(6,2)+Dijr(10,3)*r1r0+Dijr(8,3)*r2r1
c FC %     - +Dijr(9,3)*r3r2)*Iv
c FC %c
c FC %      Dijr(13,1)=(Cijr_234(4,2)+2*m0sq*Dijr(7,2)+Dijr(11,3)*r1r0+Dijr(12,3)*r2r1
c FC %     -      + Dijr(13,3)*r3r2+Iv)*Iv        
c FC %
c	RETURN
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC				
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCc
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC				
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCc
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC	
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC				
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCc
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC				
			else
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCc
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCc
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC			
c				iorder(1)=1
c				iorder(2)=2
c				iorder(3)=3
c 

				iz22=1d0/z22
                                z23=2d0*p2p3
				z33=2d0*p3sq
				z32=z32*iz22
				z23=z23-z21*z13
				z33=z33-z31*z13-z32*z23
				iz33=1d0/z33 
	  det=z11*z22*z33
c         print*, "det2",det

 


       PRr(1) = (C0r_134 - C0r_234 - D0r*r1r0   )
       PRr(2) = (C0r_124 - C0r_134 - D0r*r2r1 )
       PRr(3) = (C0r_123 - C0r_124 - D0r*r3r2 )


           Rr(3)=(PRr(3)-z32*(PRr(2)-z21*PRr(1))-z31*PRr(1))*iz33
           Rr(2)=(PRr(2)-z21*PRr(1) - z23*Rr(3))*iz22
           Rr(1)=(PRr(1)-  z12*Rr(2)-  z13*Rr(3))*iz11 	
       

       Dijr(1,1) = Rr(1)
       Dijr(2,1) = Rr(2)
       Dijr(3,1) = Rr(3)
c      
CCCCCC CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
c D00  
CCCCCC CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
       
      Dijr(7,2) = (C0r_234+ 2d0*m0sq*D0r + Dijr(1,1)*r1r0 + Dijr(2,1)*r2r1
     &     +Dijr(3,1)*r3r2)*0.5d0

c FC %CCCCCC CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
c FC %CCCCCC CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
c FC %c 30-3 1-32
c FC %       
c FC %c      
c FC %       PRr(1) = (Cijr_134(1,1) + C0r_234 - Dijr(1,1)*r1r0 - Dijr(7,2)*2.d0)
c FC %       PRr(2) = (Cijr_124(1,1) - Cijr_134(1,1) - Dijr(1,1)*r2r1        )
c FC %       PRr(3) = (Cijr_123(1,1) - Cijr_124(1,1) - Dijr(1,1)*r3r2        )     
c FC %c       
c FC %           Rr(3)=(PRr(3)-z32*(PRr(2)-z21*PRr(1))-z31*PRr(1))*iz33
c FC %           Rr(2)=(PRr(2)-z21*PRr(1) - z23*Rr(3))*iz22
c FC %           Rr(1)=(PRr(1)-  z12*Rr(2)-  z13*Rr(3))*iz11 
c FC %       
c FC %       Dijr(1,2) = Rr(1)
c FC %       Dijr(4,2) = Rr(2)
c FC %       Dijr(5,2) = Rr(3)
c FC %c       
c FC %c 33-3 4-35
c FC %       PRr(1) = (Cijr_134(1,1) -Cijr_234(1,1) - Dijr(2,1)*r1r0          )
c FC %       PRr(2)=(Cijr_124(2,1)-Cijr_134(1,1)-Dijr(2,1)*r2r1-Dijr(7,2)*2.d0)
c FC %       PRr(3) = (Cijr_123(2,1) -Cijr_124(2,1) - Dijr(2,1)*r3r2        )
c FC %       
c FC %       Rr(3)=(PRr(3)-z32*(PRr(2)-z21*PRr(1))-z31*PRr(1))*iz33
c FC %           Rr(2)=(PRr(2)-z21*PRr(1) - z23*Rr(3))*iz22
c FC %c          Rr(1)=(PRr(1)-  z12*Rr(2)-  z13*Rr(3))*iz11 
c FC %       
c FC %c       Dijr(4,2) = Rr(1)
c FC %       Dijr(2,2) = Rr(2)
c FC %       Dijr(6,2) = Rr(3)
c FC %       
c FC %c 36-3 7-38
c FC %       
c FC %       PRr(1) = (Cijr_134(2,1) - Cijr_234(2,1) - Dijr(3,1)*r1r0      )
c FC %       PRr(2) = (Cijr_124(2,1) - Cijr_134(2,1) - Dijr(3,1)*r2r1    )
c FC %       PRr(3) = (       - Cijr_124(2,1) - Dijr(3,1)*r3r2 -Dijr(7,2)*2.d0)
c FC %       
c FC %       Rr(3)=(PRr(3)-z32*(PRr(2)-z21*PRr(1))-z31*PRr(1))*iz33
c FC %c          Rr(2)=(PRr(2)-z21*PRr(1) - z23*Rr(3))*iz22
c FC %c          Rr(1)=(PRr(1)-  z12*Rr(2)-  z13*Rr(3))*iz11 
c FC %       
c FC %c       Dijr(5,2) = Rr(1)
c FC %c       Dijr(6,2) = Rr(2)
c FC %       Dijr(3,2) = Rr(3)
c FC %       
c FC %
c FC %CCCCCC CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
c FC %CCCCCC CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
c FC %
c FC %       Dijr(11,3)=0.25d0*(-C0r_234+2d0*m0sq*Dijr(1,1)+Dijr(1,2)*r1r0+Dijr(4,2)*r2r1
c FC %     &       +Dijr(5,2)*r3r2)
c FC %       
c FC %       Dijr(12,3)=0.25d0*(Cijr_234(1,1)+2d0*m0sq*Dijr(2,1)+Dijr(4,2)*r1r0+Dijr(2,2)*r2r1
c FC %     &                   +Dijr(6,2)*r3r2)
c FC %       Dijr(13,3)=0.25d0*(Cijr_234(2,1)+2d0*m0sq*Dijr(3,1)+Dijr(5,2)*r1r0 + Dijr(6,2)*r2r1
c FC %     &                    + Dijr(3,2)*r3r2)
c FC %       
c FC %
c FC %c 41-4 2-43
c FC %       
c FC %       PRr(1) =(Cijr_134(1,2) - C0r_234 - Dijr(1,2)*r1r0  - Dijr(11,3)*4.d0)
c FC %       PRr(2) =(Cijr_124(1,2) - Cijr_134(1,2) - Dijr(1,2)*r2r1       )
c FC %       PRr(3) =(Cijr_123(1,2) - Cijr_124(1,2) - Dijr(1,2)*r3r2       )
c FC %       
c FC %       Rr(3)=(PRr(3)-z32*(PRr(2)-z21*PRr(1))-z31*PRr(1))*iz33
c FC %           Rr(2)=(PRr(2)-z21*PRr(1) - z23*Rr(3))*iz22
c FC %           Rr(1)=(PRr(1)-  z12*Rr(2)-  z13*Rr(3))*iz11 
c FC %       
c FC %       Dijr(1,3) = Rr(1)
c FC %       Dijr(4,3) = Rr(2)
c FC %       Dijr(5,3) = Rr(3)
c FC %       
c FC %
c FC %c 50-5 1-52
c FC %       
c FC %       PRr(1) =  (Cijr_134(1,2) - Cijr_234(1,2) - Dijr(2,2)*r1r0)
c FC %       PRr(2) =  (Cijr_124(2,2) - Cijr_134(1,2) - Dijr(2,2)*r2r1
c FC %     &       -Dijr(12,3)*4.d0)
c FC %       PRr(3) =  (Cijr_123(2,2) - Cijr_124(2,2) - Dijr(2,2)*r3r2)
c FC %       
c FC %       Rr(3)=(PRr(3)-z32*(PRr(2)-z21*PRr(1))-z31*PRr(1))*iz33
c FC %           Rr(2)=(PRr(2)-z21*PRr(1) - z23*Rr(3))*iz22
c FC %           Rr(1)=(PRr(1)-  z12*Rr(2)-  z13*Rr(3))*iz11 
c FC %       
c FC %       Dijr(6,3) = Rr(1)
c FC %       Dijr(2,3) = Rr(2)
c FC %       Dijr(8,3) = Rr(3)
c FC %       
c FC %
c FC %       
c FC %c      
c FC %c 56-5 7-58
c FC %c      
c FC %c      
c FC %       PRr(1) = (Cijr_134(2,2) - Cijr_234(2,2) - Dijr(3,2)*r1r0)
c FC %       PRr(2) = (Cijr_124(2,2) - Cijr_134(2,2) - Dijr(3,2)*r2r1)
c FC %       PRr(3) = (     -Cijr_124(2,2) - Dijr(3,2)*r3r2 - Dijr(13,3)*4.d0)
c FC %c       
c FC %       
c FC %       Rr(3)=(PRr(3)-z32*(PRr(2)-z21*PRr(1))-z31*PRr(1))*iz33
c FC %           Rr(2)=(PRr(2)-z21*PRr(1) - z23*Rr(3))*iz22
c FC %           Rr(1)=(PRr(1)-  z12*Rr(2)-  z13*Rr(3))*iz11 
c FC %       
c FC %       Dijr(7,3) = Rr(1)
c FC %       Dijr(9,3) = Rr(2)
c FC %       Dijr(3,3) = Rr(3)
c FC %c      
c FC %
c FC %c 44-4 5-46
c FC %       PRr(1) = (Cijr_134(1,2) + Cijr_234(1,1) - Dijr(4,2)*r1r0 
c FC %     &      -Dijr(12,3)*2.d0)
c FC %       PRr(2) = (Cijr_124(3,2) - Cijr_134(1,2) - Dijr(4,2)*r2r1
c FC %     &      - Dijr(11,3)*2.d0)
c FC %       PRr(3) = (Cijr_123(3,2) - Cijr_124(3,2)  -Dijr(4,2)*r3r2)
c FC %c      
c FC %       
c FC %       Rr(3)=(PRr(3)-z32*(PRr(2)-z21*PRr(1))-z31*PRr(1))*iz33
c FC %c          Rr(2)=(PRr(2)-z21*PRr(1) - z23*Rr(3))*iz22
c FC %c          Rr(1)=(PRr(1)-  z12*Rr(2)-  z13*Rr(3))*iz11 
c FC %       
c FC %C       Dijr(4,3) = Rr(1)
c FC %C       Dijr(6,3) = Rr(2)
c FC %       Dijr(10,3) = Rr(3)
c FC %c      
c FC %c      
c FC %
c FC %c to d efine D00ij and D0000 functions!
c FC %c In P V notation we have:
c FC %c Dij( 1,7)=D416, Dij(1,8)=D417, Dij(1,9)=D418
c FC %c Dij( 1,10)=D419, Dij(1,11)=D420, Dij(1,12)=D421, D(1,13)=D422
c FC %c      
c FC %       Iv=1.d0/6.d0
c FC %c      
c FC %c 
c FC %
c FC %      Dijr(7,1)=(C0r_234+2*m0sq*Dijr(1,2)+Dijr(1,3)*r1r0+Dijr(4,3)*r2r1+Dijr(5,3)*r3r2 )*Iv
c FC %c      
c FC %      Dijr(8,1)=(Cijr_234(1,2)+2*m0sq*Dijr(2,2)+Dijr(6,3)*r1r0+Dijr(2,3)*r2r1
c FC %     - +Dijr(8,3)*r3r2)*Iv
c FC %c
c FC %      Dijr(9,1)=(Cijr_234(2,2)+2*m0sq*Dijr(3,2)+Dijr(7,3)*r1r0+Dijr(9,3)*r2r1
c FC %     - +Dijr(3,3)*r3r2)*Iv
c FC %c
c FC %      Dijr(10,1)=(-Cijr_234(1,1)+2*m0sq*Dijr(4,2)+Dijr(4,3)*r1r0+Dijr(6,3)*r2r1
c FC %     - +Dijr(10,3)*r3r2)*Iv
c FC %c      
c FC %      Dijr(11,1)=(-Cijr_234(2,1)+2*m0sq*Dijr(5,2)+Dijr(5,3)*r1r0+Dijr(10,3)*r2r1
c FC %     - +Dijr(7,3)*r3r2)*Iv
c FC %c
c FC %      Dijr(12,1)=(Cijr_234(3,2)+2*m0sq*Dijr(6,2)+Dijr(10,3)*r1r0+Dijr(8,3)*r2r1
c FC %     - +Dijr(9,3)*r3r2)*Iv
c FC %c
c FC %      Dijr(13,1)=(Cijr_234(4,2)+2*m0sq*Dijr(7,2)+Dijr(11,3)*r1r0+Dijr(12,3)*r2r1
c FC %     -      + Dijr(13,3)*r3r2+Iv)*Iv        
c FC %c
c FC %
c          RETURN
         		endif    			
CCCCCC CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
CCCCCC CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
CCCCCC CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCc				
CCCCCC CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
CCCCCC CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
CCCCCC CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCc				
CCCCCC CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
CCCCCC CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
CCCCCC CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCc				
CCCCCC CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
CCCCCC CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
CCCCCC CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCc				
              	else   ! p1p3 is the pivot
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

         		z11=2d0*p1p3
         		iz11=1d0/z11
         		z12=2d0*p2p3
         		z13=2d0*p3sq
c      
         		z21=2d0*p1p2*iz11
         		z31=2d0*p1sq*iz11
         		
         		z22=2d0*p2sq-z21*z12
         		z32=2d0*p1p2-z31*z12
c        
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
c
         		if(abs(z32).gt.abs(z22)) then	
c        			iorder(1)=2
c        			iorder(2)=3
c                               iorder(3)=1
         			piv=z21
         			z21=z31
         			z31=piv
         			piv=z22
         			z22=z32
         			z32=piv
c      
         			iz22=1d0/z22
         			z23=2d0*p1p3
         			z33=2d0*p2p3
         			z32=z32*iz22
         			z23=z23-z21*z13
         			z33=z33-z31*z13-z32*z23
         			iz33=1d0/z33
           det=z11*z22*z33
c               print*, "det3",det





       PRr(1) = (C0r_134 - C0r_234 - D0r*r1r0   )
       PRr(2) = (C0r_124 - C0r_134 - D0r*r2r1 )
       PRr(3) = (C0r_123 - C0r_124 - D0r*r3r2 )


           Rr(3)=(PRr(2)-z32*(PRr(1)-z21*PRr(3))-z31*PRr(3))*iz33
           Rr(2)=(PRr(1)-z21*PRr(3) - z23*Rr(3))*iz22
           Rr(1)=(PRr(3)-  z12*Rr(2)-  z13*Rr(3))*iz11 	
       

       Dijr(1,1) = Rr(1)
       Dijr(2,1) = Rr(2)
       Dijr(3,1) = Rr(3)
c      
CCCCCC CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
c D00  
CCCCCC CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      Dijr(7,2) = (C0r_234+ 2d0*m0sq*D0r + Dijr(1,1)*r1r0 + Dijr(2,1)*r2r1
     &     +Dijr(3,1)*r3r2)*0.5d0

     
c FC %CCCCCC CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
c FC %CCCCCC CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
c FC %c 30-3 1-32
c FC %       
c FC %c      
c FC %       PRr(1) = (Cijr_134(1,1) + C0r_234 - Dijr(1,1)*r1r0 - Dijr(7,2)*2.d0)
c FC %       PRr(2) = (Cijr_124(1,1) - Cijr_134(1,1) - Dijr(1,1)*r2r1        )
c FC %       PRr(3) = (Cijr_123(1,1) - Cijr_124(1,1) - Dijr(1,1)*r3r2        )     
c FC %c       
c FC %           Rr(3)=(PRr(2)-z32*(PRr(1)-z21*PRr(3))-z31*PRr(3))*iz33
c FC %           Rr(2)=(PRr(1)-z21*PRr(3) - z23*Rr(3))*iz22
c FC %           Rr(1)=(PRr(3)-  z12*Rr(2)-  z13*Rr(3))*iz11 
c FC %       
c FC %       Dijr(1,2) = Rr(1)
c FC %       Dijr(4,2) = Rr(2)
c FC %       Dijr(5,2) = Rr(3)
c FC %c       
c FC %       
c FC %c 33-3 4-35
c FC %       PRr(1) = (Cijr_134(1,1) -Cijr_234(1,1) - Dijr(2,1)*r1r0          )
c FC %       PRr(2)=(Cijr_124(2,1)-Cijr_134(1,1)-Dijr(2,1)*r2r1-Dijr(7,2)*2.d0)
c FC %       PRr(3) = (Cijr_123(2,1) -Cijr_124(2,1) - Dijr(2,1)*r3r2        )
c FC %       
c FC %       Rr(3)=(PRr(2)-z32*(PRr(1)-z21*PRr(3))-z31*PRr(3))*iz33
c FC %           Rr(2)=(PRr(1)-z21*PRr(3) - z23*Rr(3))*iz22
c FC %c          Rr(1)=(PRr(3)-  z12*Rr(2)-  z13*Rr(3))*iz11 
c FC %       
c FC %c       Dijr(4,2) = Rr(1)
c FC %       Dijr(2,2) = Rr(2)
c FC %       Dijr(6,2) = Rr(3)
c FC %       
c FC %   
c FC %c 36-3 7-38
c FC %       
c FC %       PRr(1) = (Cijr_134(2,1) - Cijr_234(2,1) - Dijr(3,1)*r1r0      )
c FC %       PRr(2) = (Cijr_124(2,1) - Cijr_134(2,1) - Dijr(3,1)*r2r1    )
c FC %       PRr(3) = (       - Cijr_124(2,1) - Dijr(3,1)*r3r2 -Dijr(7,2)*2.d0)
c FC %       
c FC %       Rr(3)=(PRr(2)-z32*(PRr(1)-z21*PRr(3))-z31*PRr(3))*iz33
c FC %c          Rr(2)=(PRr(1)-z21*PRr(3) - z23*Rr(3))*iz22
c FC %c          Rr(1)=(PRr(3)-  z12*Rr(2)-  z13*Rr(3))*iz11 
c FC %       
c FC %c       Dijr(5,2) = Rr(1)
c FC %c       Dijr(6,2) = Rr(2)
c FC %       Dijr(3,2) = Rr(3)
c FC %       
c FC %  
c FC %       
c FC %CCCCCC CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
c FC %CCCCCC CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
c FC %       
c FC %
c FC %       Dijr(11,3)=0.25d0*(-C0r_234+2d0*m0sq*Dijr(1,1)+Dijr(1,2)*r1r0+Dijr(4,2)*r2r1
c FC %     &       +Dijr(5,2)*r3r2)
c FC %       
c FC %       Dijr(12,3)=0.25d0*(Cijr_234(1,1)+2d0*m0sq*Dijr(2,1)+Dijr(4,2)*r1r0+Dijr(2,2)*r2r1
c FC %     &                   +Dijr(6,2)*r3r2)
c FC %       Dijr(13,3)=0.25d0*(Cijr_234(2,1)+2d0*m0sq*Dijr(3,1)+Dijr(5,2)*r1r0 + Dijr(6,2)*r2r1
c FC %     &                    + Dijr(3,2)*r3r2)
c FC %       
c FC %      
c FC %c 41-4 2-43
c FC %       
c FC %       PRr(1) =(Cijr_134(1,2) - C0r_234 - Dijr(1,2)*r1r0  - Dijr(11,3)*4.d0)
c FC %       PRr(2) =(Cijr_124(1,2) - Cijr_134(1,2) - Dijr(1,2)*r2r1       )
c FC %       PRr(3) =(Cijr_123(1,2) - Cijr_124(1,2) - Dijr(1,2)*r3r2       )
c FC %       
c FC %       Rr(3)=(PRr(2)-z32*(PRr(1)-z21*PRr(3))-z31*PRr(3))*iz33
c FC %           Rr(2)=(PRr(1)-z21*PRr(3) - z23*Rr(3))*iz22
c FC %           Rr(1)=(PRr(3)-  z12*Rr(2)-  z13*Rr(3))*iz11 
c FC %       
c FC %       Dijr(1,3) = Rr(1)
c FC %       Dijr(4,3) = Rr(2)
c FC %       Dijr(5,3) = Rr(3)
c FC %       
c FC %      
c FC %c 50-5 1-52
c FC %       
c FC %       PRr(1) =  (Cijr_134(1,2) - Cijr_234(1,2) - Dijr(2,2)*r1r0)
c FC %       PRr(2) =  (Cijr_124(2,2) - Cijr_134(1,2) - Dijr(2,2)*r2r1
c FC %     &       -Dijr(12,3)*4.d0)
c FC %       PRr(3) =  (Cijr_123(2,2) - Cijr_124(2,2) - Dijr(2,2)*r3r2)
c FC %       
c FC %       Rr(3)=(PRr(2)-z32*(PRr(1)-z21*PRr(3))-z31*PRr(3))*iz33
c FC %           Rr(2)=(PRr(1)-z21*PRr(3) - z23*Rr(3))*iz22
c FC %           Rr(1)=(PRr(3)-  z12*Rr(2)-  z13*Rr(3))*iz11 
c FC %       
c FC %       Dijr(6,3) = Rr(1)
c FC %       Dijr(2,3) = Rr(2)
c FC %       Dijr(8,3) = Rr(3)
c FC %       
c FC %
c FC %       
c FC %c      
c FC %c 56-5 7-58
c FC %c      
c FC %c      
c FC %       PRr(1) = (Cijr_134(2,2) - Cijr_234(2,2) - Dijr(3,2)*r1r0)
c FC %       PRr(2) = (Cijr_124(2,2) - Cijr_134(2,2) - Dijr(3,2)*r2r1)
c FC %       PRr(3) = (     -Cijr_124(2,2) - Dijr(3,2)*r3r2 - Dijr(13,3)*4.d0)
c FC %c       
c FC %       
c FC %       Rr(3)=(PRr(2)-z32*(PRr(1)-z21*PRr(3))-z31*PRr(3))*iz33
c FC %           Rr(2)=(PRr(1)-z21*PRr(3) - z23*Rr(3))*iz22
c FC %           Rr(1)=(PRr(3)-  z12*Rr(2)-  z13*Rr(3))*iz11 
c FC %       
c FC %       Dijr(7,3) = Rr(1)
c FC %       Dijr(9,3) = Rr(2)
c FC %       Dijr(3,3) = Rr(3)
c FC %c      
c FC %
c FC %c 44-4 5-46
c FC %       PRr(1) = (Cijr_134(1,2) + Cijr_234(1,1) - Dijr(4,2)*r1r0 
c FC %     &      -Dijr(12,3)*2.d0)
c FC %       PRr(2) = (Cijr_124(3,2) - Cijr_134(1,2) - Dijr(4,2)*r2r1
c FC %     &      - Dijr(11,3)*2.d0)
c FC %       PRr(3) = (Cijr_123(3,2) - Cijr_124(3,2)  -Dijr(4,2)*r3r2)
c FC %c      
c FC %       
c FC %       Rr(3)=(PRr(2)-z32*(PRr(1)-z21*PRr(3))-z31*PRr(3))*iz33
c FC %c          Rr(2)=(PRr(1)-z21*PRr(3) - z23*Rr(3))*iz22
c FC %c          Rr(1)=(PRr(3)-  z12*Rr(2)-  z13*Rr(3))*iz11 
c FC %       
c FC %C       Dijr(4,3) = Rr(1)
c FC %C       Dijr(6,3) = Rr(2)
c FC %       Dijr(10,3) = Rr(3)
c FC %c      
c FC %c      
c FC %
c FC %
c FC %c to d efine D00ij and D0000 functions!
c FC %c In P V notation we have:
c FC %c Dij( 1,7)=D416, Dij(1,8)=D417, Dij(1,9)=D418
c FC %c Dij( 1,10)=D419, Dij(1,11)=D420, Dij(1,12)=D421, D(1,13)=D422
c FC %c      
c FC %       Iv=1.d0/6.d0
c FC %c      
c FC %c 
c FC %
c FC %      Dijr(7,1)=(C0r_234+2*m0sq*Dijr(1,2)+Dijr(1,3)*r1r0+Dijr(4,3)*r2r1+Dijr(5,3)*r3r2 )*Iv
c FC %c      
c FC %      Dijr(8,1)=(Cijr_234(1,2)+2*m0sq*Dijr(2,2)+Dijr(6,3)*r1r0+Dijr(2,3)*r2r1
c FC %     - +Dijr(8,3)*r3r2)*Iv
c FC %c
c FC %      Dijr(9,1)=(Cijr_234(2,2)+2*m0sq*Dijr(3,2)+Dijr(7,3)*r1r0+Dijr(9,3)*r2r1
c FC %     - +Dijr(3,3)*r3r2)*Iv
c FC %c
c FC %      Dijr(10,1)=(-Cijr_234(1,1)+2*m0sq*Dijr(4,2)+Dijr(4,3)*r1r0+Dijr(6,3)*r2r1
c FC %     - +Dijr(10,3)*r3r2)*Iv
c FC %c      
c FC %      Dijr(11,1)=(-Cijr_234(2,1)+2*m0sq*Dijr(5,2)+Dijr(5,3)*r1r0+Dijr(10,3)*r2r1
c FC %     - +Dijr(7,3)*r3r2)*Iv
c FC %c
c FC %      Dijr(12,1)=(Cijr_234(3,2)+2*m0sq*Dijr(6,2)+Dijr(10,3)*r1r0+Dijr(8,3)*r2r1
c FC %     - +Dijr(9,3)*r3r2)*Iv
c FC %c
c FC %      Dijr(13,1)=(Cijr_234(4,2)+2*m0sq*Dijr(7,2)+Dijr(11,3)*r1r0+Dijr(12,3)*r2r1
c FC %     -      + Dijr(13,3)*r3r2+Iv)*Iv        
c FC %c
c FC %




c             RETURN
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
                       else
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
c        		iorder(1)=3
c        		iorder(2)=2
c                       iorder(3)=1
         	        iz22=1d0/z22	
                        z23=2d0*p2p3
                        z33=2d0*p1p3
                        z32=z32*iz22
         		z23=z23-z21*z13
         		z33=z33-z31*z13-z32*z23
         		iz33=1d0/z33
           det=-z11*z22*z33
c               print*, "det4",det         


       PRr(1) = (C0r_134 - C0r_234 - D0r*r1r0   )
       PRr(2) = (C0r_124 - C0r_134 - D0r*r2r1 )
       PRr(3) = (C0r_123 - C0r_124 - D0r*r3r2 )


           Rr(3)=(PRr(1)-z32*(PRr(2)-z21*PRr(3))-z31*PRr(3))*iz33
           Rr(2)=(PRr(2)-z21*PRr(3) - z23*Rr(3))*iz22
           Rr(1)=(PRr(3)-  z12*Rr(2)-  z13*Rr(3))*iz11 	
       

       Dijr(1,1) = Rr(1)
       Dijr(2,1) = Rr(2)
       Dijr(3,1) = Rr(3)
c      

CCCCCC CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
c D00  
CCCCCC CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      Dijr(7,2) = (C0r_234+ 2d0*m0sq*D0r + Dijr(1,1)*r1r0 + Dijr(2,1)*r2r1
     &     +Dijr(3,1)*r3r2)*0.5d0


c FC %CCCCCC CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
c FC %CCCCCC CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
c FC %c 30-3 1-32
c FC %       
c FC %c      
c FC %       PRr(1) = (Cijr_134(1,1) + C0r_234 - Dijr(1,1)*r1r0 - Dijr(7,2)*2.d0)
c FC %       PRr(2) = (Cijr_124(1,1) - Cijr_134(1,1) - Dijr(1,1)*r2r1        )
c FC %       PRr(3) = (Cijr_123(1,1) - Cijr_124(1,1) - Dijr(1,1)*r3r2        )     
c FC %c       
c FC %           Rr(3)=(PRr(1)-z32*(PRr(2)-z21*PRr(3))-z31*PRr(3))*iz33
c FC %           Rr(2)=(PRr(2)-z21*PRr(3) - z23*Rr(3))*iz22
c FC %           Rr(1)=(PRr(3)-  z12*Rr(2)-  z13*Rr(3))*iz11 
c FC %       
c FC %       Dijr(1,2) = Rr(1)
c FC %       Dijr(4,2) = Rr(2)
c FC %       Dijr(5,2) = Rr(3)
c FC %c       
c FC %
c FC %       
c FC %c 33-3 4-35
c FC %       PRr(1) = (Cijr_134(1,1) -Cijr_234(1,1) - Dijr(2,1)*r1r0          )
c FC %       PRr(2)=(Cijr_124(2,1)-Cijr_134(1,1)-Dijr(2,1)*r2r1-Dijr(7,2)*2.d0)
c FC %       PRr(3) = (Cijr_123(2,1) -Cijr_124(2,1) - Dijr(2,1)*r3r2        )
c FC %       
c FC %       Rr(3)=(PRr(1)-z32*(PRr(2)-z21*PRr(3))-z31*PRr(3))*iz33
c FC %           Rr(2)=(PRr(2)-z21*PRr(3) - z23*Rr(3))*iz22
c FC %c          Rr(1)=(PRr(3)-  z12*Rr(2)-  z13*Rr(3))*iz11 
c FC %       
c FC %c       Dijr(4,2) = Rr(1)
c FC %       Dijr(2,2) = Rr(2)
c FC %       Dijr(6,2) = Rr(3)
c FC %       
c FC %
c FC %c 36-3 7-38
c FC %       
c FC %       PRr(1) = (Cijr_134(2,1) - Cijr_234(2,1) - Dijr(3,1)*r1r0      )
c FC %       PRr(2) = (Cijr_124(2,1) - Cijr_134(2,1) - Dijr(3,1)*r2r1    )
c FC %       PRr(3) = (       - Cijr_124(2,1) - Dijr(3,1)*r3r2 -Dijr(7,2)*2.d0)
c FC %       
c FC %       Rr(3)=(PRr(1)-z32*(PRr(2)-z21*PRr(3))-z31*PRr(3))*iz33
c FC %c          Rr(2)=(PRr(2)-z21*PRr(3) - z23*Rr(3))*iz22
c FC %c          Rr(1)=(PRr(3)-  z12*Rr(2)-  z13*Rr(3))*iz11 
c FC %       
c FC %c       Dijr(5,2) = Rr(1)
c FC %c       Dijr(6,2) = Rr(2)
c FC %       Dijr(3,2) = Rr(3)
c FC %       
c FC %
c FC %       
c FC %CCCCCC CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
c FC %CCCCCC CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
c FC %
c FC %       Dijr(11,3)=0.25d0*(-C0r_234+2d0*m0sq*Dijr(1,1)+Dijr(1,2)*r1r0+Dijr(4,2)*r2r1
c FC %     &       +Dijr(5,2)*r3r2)
c FC %       
c FC %       Dijr(12,3)=0.25d0*(Cijr_234(1,1)+2d0*m0sq*Dijr(2,1)+Dijr(4,2)*r1r0+Dijr(2,2)*r2r1
c FC %     &                   +Dijr(6,2)*r3r2)
c FC %       Dijr(13,3)=0.25d0*(Cijr_234(2,1)+2d0*m0sq*Dijr(3,1)+Dijr(5,2)*r1r0 + Dijr(6,2)*r2r1
c FC %     &                    + Dijr(3,2)*r3r2)
c FC %       
c FC %
c FC %
c FC %       
c FC %       
c FC %c 41-4 2-43
c FC %       
c FC %       PRr(1) =(Cijr_134(1,2) - C0r_234 - Dijr(1,2)*r1r0  - Dijr(11,3)*4.d0)
c FC %       PRr(2) =(Cijr_124(1,2) - Cijr_134(1,2) - Dijr(1,2)*r2r1       )
c FC %       PRr(3) =(Cijr_123(1,2) - Cijr_124(1,2) - Dijr(1,2)*r3r2       )
c FC %       
c FC %       Rr(3)=(PRr(1)-z32*(PRr(2)-z21*PRr(3))-z31*PRr(3))*iz33
c FC %           Rr(2)=(PRr(2)-z21*PRr(3) - z23*Rr(3))*iz22
c FC %           Rr(1)=(PRr(3)-  z12*Rr(2)-  z13*Rr(3))*iz11 
c FC %       
c FC %       Dijr(1,3) = Rr(1)
c FC %       Dijr(4,3) = Rr(2)
c FC %       Dijr(5,3) = Rr(3)
c FC %       
c FC %
c FC %c 50-5 1-52
c FC %       
c FC %       PRr(1) =  (Cijr_134(1,2) - Cijr_234(1,2) - Dijr(2,2)*r1r0)
c FC %       PRr(2) =  (Cijr_124(2,2) - Cijr_134(1,2) - Dijr(2,2)*r2r1
c FC %     &       -Dijr(12,3)*4.d0)
c FC %       PRr(3) =  (Cijr_123(2,2) - Cijr_124(2,2) - Dijr(2,2)*r3r2)
c FC %       
c FC %       Rr(3)=(PRr(1)-z32*(PRr(2)-z21*PRr(3))-z31*PRr(3))*iz33
c FC %           Rr(2)=(PRr(2)-z21*PRr(3) - z23*Rr(3))*iz22
c FC %           Rr(1)=(PRr(3)-  z12*Rr(2)-  z13*Rr(3))*iz11 
c FC %       
c FC %       Dijr(6,3) = Rr(1)
c FC %       Dijr(2,3) = Rr(2)
c FC %       Dijr(8,3) = Rr(3)
c FC %       
c FC %
c FC %       
c FC %c      
c FC %c 56-5 7-58
c FC %c      
c FC %c      
c FC %       PRr(1) = (Cijr_134(2,2) - Cijr_234(2,2) - Dijr(3,2)*r1r0)
c FC %       PRr(2) = (Cijr_124(2,2) - Cijr_134(2,2) - Dijr(3,2)*r2r1)
c FC %       PRr(3) = (     -Cijr_124(2,2) - Dijr(3,2)*r3r2 - Dijr(13,3)*4.d0)
c FC %c       
c FC %       
c FC %       Rr(3)=(PRr(1)-z32*(PRr(2)-z21*PRr(3))-z31*PRr(3))*iz33
c FC %           Rr(2)=(PRr(2)-z21*PRr(3) - z23*Rr(3))*iz22
c FC %           Rr(1)=(PRr(3)-  z12*Rr(2)-  z13*Rr(3))*iz11 
c FC %       
c FC %       Dijr(7,3) = Rr(1)
c FC %       Dijr(9,3) = Rr(2)
c FC %       Dijr(3,3) = Rr(3)
c FC %c      
c FC %
c FC %c 44-4 5-46
c FC %       PRr(1) = (Cijr_134(1,2) + Cijr_234(1,1) - Dijr(4,2)*r1r0 
c FC %     &      -Dijr(12,3)*2.d0)
c FC %       PRr(2) = (Cijr_124(3,2) - Cijr_134(1,2) - Dijr(4,2)*r2r1
c FC %     &      - Dijr(11,3)*2.d0)
c FC %       PRr(3) = (Cijr_123(3,2) - Cijr_124(3,2)  -Dijr(4,2)*r3r2)
c FC %c      
c FC %       
c FC %       Rr(3)=(PRr(1)-z32*(PRr(2)-z21*PRr(3))-z31*PRr(3))*iz33
c FC %c          Rr(2)=(PRr(2)-z21*PRr(3) - z23*Rr(3))*iz22
c FC %c          Rr(1)=(PRr(3)-  z12*Rr(2)-  z13*Rr(3))*iz11 
c FC %       
c FC %C       Dijr(4,3) = Rr(1)
c FC %C       Dijr(6,3) = Rr(2)
c FC %       Dijr(10,3) = Rr(3)
c FC %c      
c FC %c      
c FC %
c FC %
c FC %
c FC %c to d efine D00ij and D0000 functions!
c FC %c In P V notation we have:
c FC %c Dij( 1,7)=D416, Dij(1,8)=D417, Dij(1,9)=D418
c FC %c Dij( 1,10)=D419, Dij(1,11)=D420, Dij(1,12)=D421, D(1,13)=D422
c FC %c      
c FC %       Iv=1.d0/6.d0
c FC %c      
c FC %c 
c FC %
c FC %      Dijr(7,1)=(C0r_234+2*m0sq*Dijr(1,2)+Dijr(1,3)*r1r0+Dijr(4,3)*r2r1+Dijr(5,3)*r3r2 )*Iv
c FC %c      
c FC %      Dijr(8,1)=(Cijr_234(1,2)+2*m0sq*Dijr(2,2)+Dijr(6,3)*r1r0+Dijr(2,3)*r2r1
c FC %     - +Dijr(8,3)*r3r2)*Iv
c FC %c
c FC %      Dijr(9,1)=(Cijr_234(2,2)+2*m0sq*Dijr(3,2)+Dijr(7,3)*r1r0+Dijr(9,3)*r2r1
c FC %     - +Dijr(3,3)*r3r2)*Iv
c FC %c
c FC %      Dijr(10,1)=(-Cijr_234(1,1)+2*m0sq*Dijr(4,2)+Dijr(4,3)*r1r0+Dijr(6,3)*r2r1
c FC %     - +Dijr(10,3)*r3r2)*Iv
c FC %c      
c FC %      Dijr(11,1)=(-Cijr_234(2,1)+2*m0sq*Dijr(5,2)+Dijr(5,3)*r1r0+Dijr(10,3)*r2r1
c FC %     - +Dijr(7,3)*r3r2)*Iv
c FC %c
c FC %      Dijr(12,1)=(Cijr_234(3,2)+2*m0sq*Dijr(6,2)+Dijr(10,3)*r1r0+Dijr(8,3)*r2r1
c FC %     - +Dijr(9,3)*r3r2)*Iv
c FC %c
c FC %      Dijr(13,1)=(Cijr_234(4,2)+2*m0sq*Dijr(7,2)+Dijr(11,3)*r1r0+Dijr(12,3)*r2r1
c FC %     -      + Dijr(13,3)*r3r2+Iv)*Iv        
c FC %c
c FC %


c      RETURN
         		endif
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC         		
         	endif  ! p1sq or p1p3 is the pivot
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
         else  ! p1p2 or p1p3 is the pivot
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
             if(abs(p1p2).gt.abs(p1p3)) then
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
         		z11=2d0*p1p2
         		iz11=1d0/z11
         		z12=2d0*p2sq
         		z13=2d0*p2p3
             
         		z21=2d0*p1sq*iz11
         		z31=2d0*p1p3*iz11
         		
         		z22=z11-z21*z12
         		z32=z13-z31*z12
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC       
        		if(abs(z32).gt.abs(z22)) then	
c        			iorder(1)=3
c        			iorder(2)=1
c                               iorder(3)=2
         			piv=z21
         			z21=z31
         			z31=piv
         			piv=z22
         			z22=z32
         			z32=piv
c      
         			iz22=1d0/z22 
         			z23=2d0*p3sq
         			z33=2d0*p1p3
         			z32=z32*iz22
         			z23=z23-z21*z13
         			z33=z33-z31*z13-z32*z23
           			iz33=1d0/z33 
           det=z11*z22*z33
c           print*, "det5",det





       PRr(1) = (C0r_134 - C0r_234 - D0r*r1r0   )
       PRr(2) = (C0r_124 - C0r_134 - D0r*r2r1 )
       PRr(3) = (C0r_123 - C0r_124 - D0r*r3r2 )


           Rr(3)=(PRr(1)-z32*(PRr(3)-z21*PRr(2))-z31*PRr(2))*iz33
           Rr(2)=(PRr(3)-z21*PRr(2) - z23*Rr(3))*iz22
           Rr(1)=(PRr(2)-  z12*Rr(2)-  z13*Rr(3))*iz11 	
       

       Dijr(1,1) = Rr(1)
       Dijr(2,1) = Rr(2)
       Dijr(3,1) = Rr(3)
c      

CCCCCC CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
c D00  
CCCCCC CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      Dijr(7,2) = (C0r_234+ 2d0*m0sq*D0r + Dijr(1,1)*r1r0 + Dijr(2,1)*r2r1
     &     +Dijr(3,1)*r3r2)*0.5d0


       
c FC %CCCCCC CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
c FC %CCCCCC CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
c FC %c 30-3 1-32
c FC %       
c FC %c      
c FC %       PRr(1) = (Cijr_134(1,1) + C0r_234 - Dijr(1,1)*r1r0 - Dijr(7,2)*2.d0)
c FC %       PRr(2) = (Cijr_124(1,1) - Cijr_134(1,1) - Dijr(1,1)*r2r1        )
c FC %       PRr(3) = (Cijr_123(1,1) - Cijr_124(1,1) - Dijr(1,1)*r3r2        )     
c FC %c       
c FC %           Rr(3)=(PRr(1)-z32*(PRr(3)-z21*PRr(2))-z31*PRr(2))*iz33
c FC %           Rr(2)=(PRr(3)-z21*PRr(2) - z23*Rr(3))*iz22
c FC %           Rr(1)=(PRr(2)-  z12*Rr(2)-  z13*Rr(3))*iz11 
c FC %       
c FC %       Dijr(1,2) = Rr(1)
c FC %       Dijr(4,2) = Rr(2)
c FC %       Dijr(5,2) = Rr(3)
c FC %c       
c FC %
c FC %       
c FC %c 33-3 4-35
c FC %       PRr(1) = (Cijr_134(1,1) -Cijr_234(1,1) - Dijr(2,1)*r1r0          )
c FC %       PRr(2)=(Cijr_124(2,1)-Cijr_134(1,1)-Dijr(2,1)*r2r1-Dijr(7,2)*2.d0)
c FC %       PRr(3) = (Cijr_123(2,1) -Cijr_124(2,1) - Dijr(2,1)*r3r2        )
c FC %       
c FC %       Rr(3)=(PRr(1)-z32*(PRr(3)-z21*PRr(2))-z31*PRr(2))*iz33
c FC %           Rr(2)=(PRr(3)-z21*PRr(2) - z23*Rr(3))*iz22
c FC %c          Rr(1)=(PRr(2)-  z12*Rr(2)-  z13*Rr(3))*iz11 
c FC %       
c FC %c       Dijr(4,2) = Rr(1)
c FC %       Dijr(2,2) = Rr(2)
c FC %       Dijr(6,2) = Rr(3)
c FC %       
c FC %
c FC %       
c FC %c 36-3 7-38
c FC %       
c FC %       PRr(1) = (Cijr_134(2,1) - Cijr_234(2,1) - Dijr(3,1)*r1r0      )
c FC %       PRr(2) = (Cijr_124(2,1) - Cijr_134(2,1) - Dijr(3,1)*r2r1    )
c FC %       PRr(3) = (       - Cijr_124(2,1) - Dijr(3,1)*r3r2 -Dijr(7,2)*2.d0)
c FC %       
c FC %       Rr(3)=(PRr(1)-z32*(PRr(3)-z21*PRr(2))-z31*PRr(2))*iz33
c FC %c          Rr(2)=(PRr(3)-z21*PRr(2) - z23*Rr(3))*iz22
c FC %c          Rr(1)=(PRr(2)-  z12*Rr(2)-  z13*Rr(3))*iz11 
c FC %       
c FC %c       Dijr(5,2) = Rr(1)
c FC %c       Dijr(6,2) = Rr(2)
c FC %       Dijr(3,2) = Rr(3)
c FC %       
c FC %
c FC %       
c FC %       
c FC %       
c FC %CCCCCC CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
c FC %CCCCCC CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
c FC %
c FC %       Dijr(11,3)=0.25d0*(-C0r_234+2d0*m0sq*Dijr(1,1)+Dijr(1,2)*r1r0+Dijr(4,2)*r2r1
c FC %     &       +Dijr(5,2)*r3r2)
c FC %       
c FC %       Dijr(12,3)=0.25d0*(Cijr_234(1,1)+2d0*m0sq*Dijr(2,1)+Dijr(4,2)*r1r0+Dijr(2,2)*r2r1
c FC %     &                   +Dijr(6,2)*r3r2)
c FC %       Dijr(13,3)=0.25d0*(Cijr_234(2,1)+2d0*m0sq*Dijr(3,1)+Dijr(5,2)*r1r0 + Dijr(6,2)*r2r1
c FC %     &                    + Dijr(3,2)*r3r2)
c FC %       
c FC %
c FC %       
c FC %   
c FC %       
c FC %       
c FC %c 41-4 2-43
c FC %       
c FC %       PRr(1) =(Cijr_134(1,2) - C0r_234 - Dijr(1,2)*r1r0  - Dijr(11,3)*4.d0)
c FC %       PRr(2) =(Cijr_124(1,2) - Cijr_134(1,2) - Dijr(1,2)*r2r1       )
c FC %       PRr(3) =(Cijr_123(1,2) - Cijr_124(1,2) - Dijr(1,2)*r3r2       )
c FC %       
c FC %       Rr(3)=(PRr(1)-z32*(PRr(3)-z21*PRr(2))-z31*PRr(2))*iz33
c FC %           Rr(2)=(PRr(3)-z21*PRr(2) - z23*Rr(3))*iz22
c FC %           Rr(1)=(PRr(2)-  z12*Rr(2)-  z13*Rr(3))*iz11 
c FC %       
c FC %       Dijr(1,3) = Rr(1)
c FC %       Dijr(4,3) = Rr(2)
c FC %       Dijr(5,3) = Rr(3)
c FC %       
c FC % 
c FC %       
c FC %c 50-5 1-52
c FC %       
c FC %       PRr(1) =  (Cijr_134(1,2) - Cijr_234(1,2) - Dijr(2,2)*r1r0)
c FC %       PRr(2) =  (Cijr_124(2,2) - Cijr_134(1,2) - Dijr(2,2)*r2r1
c FC %     &       -Dijr(12,3)*4.d0)
c FC %       PRr(3) =  (Cijr_123(2,2) - Cijr_124(2,2) - Dijr(2,2)*r3r2)
c FC %       
c FC %       Rr(3)=(PRr(1)-z32*(PRr(3)-z21*PRr(2))-z31*PRr(2))*iz33
c FC %           Rr(2)=(PRr(3)-z21*PRr(2) - z23*Rr(3))*iz22
c FC %           Rr(1)=(PRr(2)-  z12*Rr(2)-  z13*Rr(3))*iz11 
c FC %       
c FC %       Dijr(6,3) = Rr(1)
c FC %       Dijr(2,3) = Rr(2)
c FC %       Dijr(8,3) = Rr(3)
c FC %       
c FC %
c FC %       
c FC %c      
c FC %c 56-5 7-58
c FC %c      
c FC %c      
c FC %       PRr(1) = (Cijr_134(2,2) - Cijr_234(2,2) - Dijr(3,2)*r1r0)
c FC %       PRr(2) = (Cijr_124(2,2) - Cijr_134(2,2) - Dijr(3,2)*r2r1)
c FC %       PRr(3) = (     -Cijr_124(2,2) - Dijr(3,2)*r3r2 - Dijr(13,3)*4.d0)
c FC %c       
c FC %       
c FC %       Rr(3)=(PRr(1)-z32*(PRr(3)-z21*PRr(2))-z31*PRr(2))*iz33
c FC %           Rr(2)=(PRr(3)-z21*PRr(2) - z23*Rr(3))*iz22
c FC %           Rr(1)=(PRr(2)-  z12*Rr(2)-  z13*Rr(3))*iz11 
c FC %       
c FC %       Dijr(7,3) = Rr(1)
c FC %       Dijr(9,3) = Rr(2)
c FC %       Dijr(3,3) = Rr(3)
c FC %c      
c FC %
c FC %c 44-4 5-46
c FC %       PRr(1) = (Cijr_134(1,2) + Cijr_234(1,1) - Dijr(4,2)*r1r0 
c FC %     &      -Dijr(12,3)*2.d0)
c FC %       PRr(2) = (Cijr_124(3,2) - Cijr_134(1,2) - Dijr(4,2)*r2r1
c FC %     &      - Dijr(11,3)*2.d0)
c FC %       PRr(3) = (Cijr_123(3,2) - Cijr_124(3,2)  -Dijr(4,2)*r3r2)
c FC %c      
c FC %       
c FC %       Rr(3)=(PRr(1)-z32*(PRr(3)-z21*PRr(2))-z31*PRr(2))*iz33
c FC %c          Rr(2)=(PRr(3)-z21*PRr(2) - z23*Rr(3))*iz22
c FC %c          Rr(1)=(PRr(2)-  z12*Rr(2)-  z13*Rr(3))*iz11 
c FC %       
c FC %C       Dijr(4,3) = Rr(1)
c FC %C       Dijr(6,3) = Rr(2)
c FC %       Dijr(10,3) = Rr(3)
c FC %c      
c FC %c      
c FC %
c FC %
c FC %c to d efine D00ij and D0000 functions!
c FC %c In P V notation we have:
c FC %c Dij( 1,7)=D416, Dij(1,8)=D417, Dij(1,9)=D418
c FC %c Dij( 1,10)=D419, Dij(1,11)=D420, Dij(1,12)=D421, D(1,13)=D422
c FC %c      
c FC %       Iv=1.d0/6.d0
c FC %c      
c FC %c 
c FC %
c FC %      Dijr(7,1)=(C0r_234+2*m0sq*Dijr(1,2)+Dijr(1,3)*r1r0+Dijr(4,3)*r2r1+Dijr(5,3)*r3r2 )*Iv
c FC %c      
c FC %      Dijr(8,1)=(Cijr_234(1,2)+2*m0sq*Dijr(2,2)+Dijr(6,3)*r1r0+Dijr(2,3)*r2r1
c FC %     - +Dijr(8,3)*r3r2)*Iv
c FC %c
c FC %      Dijr(9,1)=(Cijr_234(2,2)+2*m0sq*Dijr(3,2)+Dijr(7,3)*r1r0+Dijr(9,3)*r2r1
c FC %     - +Dijr(3,3)*r3r2)*Iv
c FC %c
c FC %      Dijr(10,1)=(-Cijr_234(1,1)+2*m0sq*Dijr(4,2)+Dijr(4,3)*r1r0+Dijr(6,3)*r2r1
c FC %     - +Dijr(10,3)*r3r2)*Iv
c FC %c      
c FC %      Dijr(11,1)=(-Cijr_234(2,1)+2*m0sq*Dijr(5,2)+Dijr(5,3)*r1r0+Dijr(10,3)*r2r1
c FC %     - +Dijr(7,3)*r3r2)*Iv
c FC %c
c FC %      Dijr(12,1)=(Cijr_234(3,2)+2*m0sq*Dijr(6,2)+Dijr(10,3)*r1r0+Dijr(8,3)*r2r1
c FC %     - +Dijr(9,3)*r3r2)*Iv
c FC %c
c FC %      Dijr(13,1)=(Cijr_234(4,2)+2*m0sq*Dijr(7,2)+Dijr(11,3)*r1r0+Dijr(12,3)*r2r1
c FC %     -      + Dijr(13,3)*r3r2+Iv)*Iv        
c FC %c
c FC %


c             	RETURN
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC   
              		else
c        			iorder(1)=2
c        			iorder(2)=1
c        			iorder(3)=3
c      
         			z23=2d0*p1p3
         			z33=2d0*p3sq
         			iz22=1d0/z22
c        			det=-det
         			z32=z32*iz22
         			z23=z23-z21*z13
         			z33=z33-z31*z13-z32*z23
         			iz33=1d0/z33
            det=-z11*z22*z33
         	
c           print*, "det6",det
       
      




       PRr(1) = (C0r_134 - C0r_234 - D0r*r1r0   )
       PRr(2) = (C0r_124 - C0r_134 - D0r*r2r1 )
       PRr(3) = (C0r_123 - C0r_124 - D0r*r3r2 )


           Rr(3)=(PRr(3)-z32*(PRr(1)-z21*PRr(2))-z31*PRr(2))*iz33
           Rr(2)=(PRr(1)-z21*PRr(2) - z23*Rr(3))*iz22
           Rr(1)=(PRr(2)-  z12*Rr(2)-  z13*Rr(3))*iz11 	
       

       Dijr(1,1) = Rr(1)
       Dijr(2,1) = Rr(2)
       Dijr(3,1) = Rr(3)
c      

CCCCCC CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
c D00  
CCCCCC CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      Dijr(7,2) = (C0r_234+ 2d0*m0sq*D0r + Dijr(1,1)*r1r0 + Dijr(2,1)*r2r1
     &     +Dijr(3,1)*r3r2)*0.5d0



CCCCCC CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
c FC %CCCCCC CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
c FC %c 30-3 1-32
c FC %       
c FC %c      
c FC %       PRr(1) = (Cijr_134(1,1) + C0r_234 - Dijr(1,1)*r1r0 - Dijr(7,2)*2.d0)
c FC %       PRr(2) = (Cijr_124(1,1) - Cijr_134(1,1) - Dijr(1,1)*r2r1        )
c FC %       PRr(3) = (Cijr_123(1,1) - Cijr_124(1,1) - Dijr(1,1)*r3r2        )     
c FC %c       
c FC %           Rr(3)=(PRr(3)-z32*(PRr(1)-z21*PRr(2))-z31*PRr(2))*iz33
c FC %           Rr(2)=(PRr(1)-z21*PRr(2) - z23*Rr(3))*iz22
c FC %           Rr(1)=(PRr(2)-  z12*Rr(2)-  z13*Rr(3))*iz11 
c FC %       
c FC %       Dijr(1,2) = Rr(1)
c FC %       Dijr(4,2) = Rr(2)
c FC %       Dijr(5,2) = Rr(3)
c FC %c       
c FC %
c FC %       
c FC %c 33-3 4-35
c FC %       PRr(1) = (Cijr_134(1,1) -Cijr_234(1,1) - Dijr(2,1)*r1r0          )
c FC %       PRr(2)=(Cijr_124(2,1)-Cijr_134(1,1)-Dijr(2,1)*r2r1-Dijr(7,2)*2.d0)
c FC %       PRr(3) = (Cijr_123(2,1) -Cijr_124(2,1) - Dijr(2,1)*r3r2        )
c FC %       
c FC %       Rr(3)=(PRr(3)-z32*(PRr(1)-z21*PRr(2))-z31*PRr(2))*iz33
c FC %           Rr(2)=(PRr(1)-z21*PRr(2) - z23*Rr(3))*iz22
c FC %c          Rr(1)=(PRr(2)-  z12*Rr(2)-  z13*Rr(3))*iz11 
c FC %       
c FC %c       Dijr(4,2) = Rr(1)
c FC %       Dijr(2,2) = Rr(2)
c FC %       Dijr(6,2) = Rr(3)
c FC %       
c FC %
c FC %       
c FC %c 36-3 7-38
c FC %       
c FC %       PRr(1) = (Cijr_134(2,1) - Cijr_234(2,1) - Dijr(3,1)*r1r0      )
c FC %       PRr(2) = (Cijr_124(2,1) - Cijr_134(2,1) - Dijr(3,1)*r2r1    )
c FC %       PRr(3) = (       - Cijr_124(2,1) - Dijr(3,1)*r3r2 -Dijr(7,2)*2.d0)
c FC %       
c FC %       Rr(3)=(PRr(3)-z32*(PRr(1)-z21*PRr(2))-z31*PRr(2))*iz33
c FC %c          Rr(2)=(PRr(1)-z21*PRr(2) - z23*Rr(3))*iz22
c FC %c          Rr(1)=(PRr(2)-  z12*Rr(2)-  z13*Rr(3))*iz11 
c FC %       
c FC %c       Dijr(5,2) = Rr(1)
c FC %c       Dijr(6,2) = Rr(2)
c FC %       Dijr(3,2) = Rr(3)
c FC %       
c FC %
c FC %       
c FC %CCCCCC CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
c FC %CCCCCC CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
c FC %       
c FC %       Dijr(11,3)=0.25d0*(-C0r_234+2d0*m0sq*Dijr(1,1)+Dijr(1,2)*r1r0+Dijr(4,2)*r2r1
c FC %     &       +Dijr(5,2)*r3r2)
c FC %       
c FC %       Dijr(12,3)=0.25d0*(Cijr_234(1,1)+2d0*m0sq*Dijr(2,1)+Dijr(4,2)*r1r0+Dijr(2,2)*r2r1
c FC %     &                   +Dijr(6,2)*r3r2)
c FC %       Dijr(13,3)=0.25d0*(Cijr_234(2,1)+2d0*m0sq*Dijr(3,1)+Dijr(5,2)*r1r0 + Dijr(6,2)*r2r1
c FC %     &                    + Dijr(3,2)*r3r2)
c FC %       
c FC %
c FC %     
c FC %       
c FC %       
c FC %c 41-4 2-43
c FC %       
c FC %       PRr(1) =(Cijr_134(1,2) - C0r_234 - Dijr(1,2)*r1r0  - Dijr(11,3)*4.d0)
c FC %       PRr(2) =(Cijr_124(1,2) - Cijr_134(1,2) - Dijr(1,2)*r2r1       )
c FC %       PRr(3) =(Cijr_123(1,2) - Cijr_124(1,2) - Dijr(1,2)*r3r2       )
c FC %       
c FC %       Rr(3)=(PRr(3)-z32*(PRr(1)-z21*PRr(2))-z31*PRr(2))*iz33
c FC %           Rr(2)=(PRr(1)-z21*PRr(2) - z23*Rr(3))*iz22
c FC %           Rr(1)=(PRr(2)-  z12*Rr(2)-  z13*Rr(3))*iz11 
c FC %       
c FC %       Dijr(1,3) = Rr(1)
c FC %       Dijr(4,3) = Rr(2)
c FC %       Dijr(5,3) = Rr(3)
c FC %       
c FC %
c FC %       
c FC %c 50-5 1-52
c FC %       
c FC %       PRr(1) =  (Cijr_134(1,2) - Cijr_234(1,2) - Dijr(2,2)*r1r0)
c FC %       PRr(2) =  (Cijr_124(2,2) - Cijr_134(1,2) - Dijr(2,2)*r2r1
c FC %     &       -Dijr(12,3)*4.d0)
c FC %       PRr(3) =  (Cijr_123(2,2) - Cijr_124(2,2) - Dijr(2,2)*r3r2)
c FC %       
c FC %       Rr(3)=(PRr(3)-z32*(PRr(1)-z21*PRr(2))-z31*PRr(2))*iz33
c FC %           Rr(2)=(PRr(1)-z21*PRr(2) - z23*Rr(3))*iz22
c FC %           Rr(1)=(PRr(2)-  z12*Rr(2)-  z13*Rr(3))*iz11 
c FC %       
c FC %       Dijr(6,3) = Rr(1)
c FC %       Dijr(2,3) = Rr(2)
c FC %       Dijr(8,3) = Rr(3)
c FC %       
c FC %
c FC %       
c FC %c      
c FC %c 56-5 7-58
c FC %c      
c FC %c      
c FC %       PRr(1) = (Cijr_134(2,2) - Cijr_234(2,2) - Dijr(3,2)*r1r0)
c FC %       PRr(2) = (Cijr_124(2,2) - Cijr_134(2,2) - Dijr(3,2)*r2r1)
c FC %       PRr(3) = (     -Cijr_124(2,2) - Dijr(3,2)*r3r2 - Dijr(13,3)*4.d0)
c FC %c       
c FC %       
c FC %       Rr(3)=(PRr(3)-z32*(PRr(1)-z21*PRr(2))-z31*PRr(2))*iz33
c FC %           Rr(2)=(PRr(1)-z21*PRr(2) - z23*Rr(3))*iz22
c FC %           Rr(1)=(PRr(2)-  z12*Rr(2)-  z13*Rr(3))*iz11 
c FC %       
c FC %       Dijr(7,3) = Rr(1)
c FC %       Dijr(9,3) = Rr(2)
c FC %       Dijr(3,3) = Rr(3)
c FC %c      
c FC %
c FC %c 44-4 5-46
c FC %       PRr(1) = (Cijr_134(1,2) + Cijr_234(1,1) - Dijr(4,2)*r1r0 
c FC %     &      -Dijr(12,3)*2.d0)
c FC %       PRr(2) = (Cijr_124(3,2) - Cijr_134(1,2) - Dijr(4,2)*r2r1
c FC %     &      - Dijr(11,3)*2.d0)
c FC %       PRr(3) = (Cijr_123(3,2) - Cijr_124(3,2)  -Dijr(4,2)*r3r2)
c FC %c      
c FC %       
c FC %       Rr(3)=(PRr(3)-z32*(PRr(1)-z21*PRr(2))-z31*PRr(2))*iz33
c FC %c          Rr(2)=(PRr(1)-z21*PRr(2) - z23*Rr(3))*iz22
c FC %c          Rr(1)=(PRr(2)-  z12*Rr(2)-  z13*Rr(3))*iz11 
c FC %       
c FC %C       Dijr(4,3) = Rr(1)
c FC %C       Dijr(6,3) = Rr(2)
c FC %       Dijr(10,3) = Rr(3)
c FC %c      
c FC %c      
c FC %
c FC %
c FC %c to d efine D00ij and D0000 functions!
c FC %c In P V notation we have:
c FC %c Dij( 1,7)=D416, Dij(1,8)=D417, Dij(1,9)=D418
c FC %c Dij( 1,10)=D419, Dij(1,11)=D420, Dij(1,12)=D421, D(1,13)=D422
c FC %c      
c FC %       Iv=1.d0/6.d0
c FC %c      
c FC %c 
c FC %
c FC %      Dijr(7,1)=(C0r_234+2*m0sq*Dijr(1,2)+Dijr(1,3)*r1r0+Dijr(4,3)*r2r1+Dijr(5,3)*r3r2 )*Iv
c FC %c      
c FC %      Dijr(8,1)=(Cijr_234(1,2)+2*m0sq*Dijr(2,2)+Dijr(6,3)*r1r0+Dijr(2,3)*r2r1
c FC %     - +Dijr(8,3)*r3r2)*Iv
c FC %c
c FC %      Dijr(9,1)=(Cijr_234(2,2)+2*m0sq*Dijr(3,2)+Dijr(7,3)*r1r0+Dijr(9,3)*r2r1
c FC %     - +Dijr(3,3)*r3r2)*Iv
c FC %c
c FC %      Dijr(10,1)=(-Cijr_234(1,1)+2*m0sq*Dijr(4,2)+Dijr(4,3)*r1r0+Dijr(6,3)*r2r1
c FC %     - +Dijr(10,3)*r3r2)*Iv
c FC %c      
c FC %      Dijr(11,1)=(-Cijr_234(2,1)+2*m0sq*Dijr(5,2)+Dijr(5,3)*r1r0+Dijr(10,3)*r2r1
c FC %     - +Dijr(7,3)*r3r2)*Iv
c FC %c
c FC %      Dijr(12,1)=(Cijr_234(3,2)+2*m0sq*Dijr(6,2)+Dijr(10,3)*r1r0+Dijr(8,3)*r2r1
c FC %     - +Dijr(9,3)*r3r2)*Iv
c FC %c
c FC %      Dijr(13,1)=(Cijr_234(4,2)+2*m0sq*Dijr(7,2)+Dijr(11,3)*r1r0+Dijr(12,3)*r2r1
c FC %     -      + Dijr(13,3)*r3r2+Iv)*Iv        
c FC %c
c FC %



c      RETURN	
	      		endif	
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC         			
              else
         		z11=2d0*p1p3
         		iz11=1d0/z11
         		z12=2d0*p2p3
         		z13=2d0*p3sq
c      
         		z21=2d0*p1p2*iz11
         		z31=2d0*p1sq*iz11
         		
         		z22=2d0*p2sq-z21*z12
         		z32=2d0*p1p2-z31*z12
c  
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC               
         		if(abs(z32).gt.abs(z22)) then	
c        			iorder(1)=2
c        			iorder(2)=3
c                               iorder(3)=1
         			piv=z21
         			z21=z31
         			z31=piv
         			piv=z22
         			z22=z32
         			z32=piv
c      
         			iz22=1d0/z22
         			z23=2d0*p1p3
         			z33=2d0*p2p3
         			z32=z32*iz22
         			z23=z23-z21*z13
         			z33=z33-z31*z13-z32*z23
         			iz33=1d0/z33
           det=z11*z22*z33
c           print*, "det7",det

 


       PRr(1) = (C0r_134 - C0r_234 - D0r*r1r0   )
       PRr(2) = (C0r_124 - C0r_134 - D0r*r2r1 )
       PRr(3) = (C0r_123 - C0r_124 - D0r*r3r2 )


           Rr(3)=(PRr(2)-z32*(PRr(1)-z21*PRr(3))-z31*PRr(3))*iz33
           Rr(2)=(PRr(1)-z21*PRr(3) - z23*Rr(3))*iz22
           Rr(1)=(PRr(3)-  z12*Rr(2)-  z13*Rr(3))*iz11 	
       

       Dijr(1,1) = Rr(1)
       Dijr(2,1) = Rr(2)
       Dijr(3,1) = Rr(3)
c      

CCCCCC CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
c D00  
CCCCCC CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      Dijr(7,2) = (C0r_234+ 2d0*m0sq*D0r + Dijr(1,1)*r1r0 + Dijr(2,1)*r2r1
     &     +Dijr(3,1)*r3r2)*0.5d0


CCCCCC CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
CCCCCC CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
c FC %c 30-3 1-32
c FC %       
c FC %c      
c FC %       PRr(1) = (Cijr_134(1,1) + C0r_234 - Dijr(1,1)*r1r0 - Dijr(7,2)*2.d0)
c FC %       PRr(2) = (Cijr_124(1,1) - Cijr_134(1,1) - Dijr(1,1)*r2r1        )
c FC %       PRr(3) = (Cijr_123(1,1) - Cijr_124(1,1) - Dijr(1,1)*r3r2        )     
c FC %c       
c FC %           Rr(3)=(PRr(2)-z32*(PRr(1)-z21*PRr(3))-z31*PRr(3))*iz33
c FC %           Rr(2)=(PRr(1)-z21*PRr(3) - z23*Rr(3))*iz22
c FC %           Rr(1)=(PRr(3)-  z12*Rr(2)-  z13*Rr(3))*iz11 
c FC %       
c FC %       Dijr(1,2) = Rr(1)
c FC %       Dijr(4,2) = Rr(2)
c FC %       Dijr(5,2) = Rr(3)
c FC %c       
c FC %
c FC %       
c FC %c 33-3 4-35
c FC %       PRr(1) = (Cijr_134(1,1) -Cijr_234(1,1) - Dijr(2,1)*r1r0          )
c FC %       PRr(2)=(Cijr_124(2,1)-Cijr_134(1,1)-Dijr(2,1)*r2r1-Dijr(7,2)*2.d0)
c FC %       PRr(3) = (Cijr_123(2,1) -Cijr_124(2,1) - Dijr(2,1)*r3r2        )
c FC %       
c FC %       Rr(3)=(PRr(2)-z32*(PRr(1)-z21*PRr(3))-z31*PRr(3))*iz33
c FC %           Rr(2)=(PRr(1)-z21*PRr(3) - z23*Rr(3))*iz22
c FC %c          Rr(1)=(PRr(3)-  z12*Rr(2)-  z13*Rr(3))*iz11 
c FC %       
c FC %c       Dijr(4,2) = Rr(1)
c FC %       Dijr(2,2) = Rr(2)
c FC %       Dijr(6,2) = Rr(3)
c FC %       
c FC %
c FC %       
c FC %c 36-3 7-38
c FC %       
c FC %       PRr(1) = (Cijr_134(2,1) - Cijr_234(2,1) - Dijr(3,1)*r1r0      )
c FC %       PRr(2) = (Cijr_124(2,1) - Cijr_134(2,1) - Dijr(3,1)*r2r1    )
c FC %       PRr(3) = (       - Cijr_124(2,1) - Dijr(3,1)*r3r2 -Dijr(7,2)*2.d0)
c FC %       
c FC %       Rr(3)=(PRr(2)-z32*(PRr(1)-z21*PRr(3))-z31*PRr(3))*iz33
c FC %c          Rr(2)=(PRr(1)-z21*PRr(3) - z23*Rr(3))*iz22
c FC %c          Rr(1)=(PRr(3)-  z12*Rr(2)-  z13*Rr(3))*iz11 
c FC %       
c FC %c       Dijr(5,2) = Rr(1)
c FC %c       Dijr(6,2) = Rr(2)
c FC %       Dijr(3,2) = Rr(3)
c FC %       
c FC %
c FC %       
c FC %CCCCCC CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
c FC %CCCCCC CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
c FC %       
c FC %       Dijr(11,3)=0.25d0*(-C0r_234+2d0*m0sq*Dijr(1,1)+Dijr(1,2)*r1r0+Dijr(4,2)*r2r1
c FC %     &       +Dijr(5,2)*r3r2)
c FC %       
c FC %       Dijr(12,3)=0.25d0*(Cijr_234(1,1)+2d0*m0sq*Dijr(2,1)+Dijr(4,2)*r1r0+Dijr(2,2)*r2r1
c FC %     &                   +Dijr(6,2)*r3r2)
c FC %       Dijr(13,3)=0.25d0*(Cijr_234(2,1)+2d0*m0sq*Dijr(3,1)+Dijr(5,2)*r1r0 + Dijr(6,2)*r2r1
c FC %     &                    + Dijr(3,2)*r3r2)
c FC %       
c FC %
c FC %      
c FC %       
c FC %       
c FC %c 41-4 2-43
c FC %       
c FC %       PRr(1) =(Cijr_134(1,2) - C0r_234 - Dijr(1,2)*r1r0  - Dijr(11,3)*4.d0)
c FC %       PRr(2) =(Cijr_124(1,2) - Cijr_134(1,2) - Dijr(1,2)*r2r1       )
c FC %       PRr(3) =(Cijr_123(1,2) - Cijr_124(1,2) - Dijr(1,2)*r3r2       )
c FC %       
c FC %       Rr(3)=(PRr(2)-z32*(PRr(1)-z21*PRr(3))-z31*PRr(3))*iz33
c FC %           Rr(2)=(PRr(1)-z21*PRr(3) - z23*Rr(3))*iz22
c FC %           Rr(1)=(PRr(3)-  z12*Rr(2)-  z13*Rr(3))*iz11 
c FC %       
c FC %       Dijr(1,3) = Rr(1)
c FC %       Dijr(4,3) = Rr(2)
c FC %       Dijr(5,3) = Rr(3)
c FC %       
c FC %
c FC %       
c FC %c 50-5 1-52
c FC %       
c FC %       PRr(1) =  (Cijr_134(1,2) - Cijr_234(1,2) - Dijr(2,2)*r1r0)
c FC %       PRr(2) =  (Cijr_124(2,2) - Cijr_134(1,2) - Dijr(2,2)*r2r1
c FC %     &       -Dijr(12,3)*4.d0)
c FC %       PRr(3) =  (Cijr_123(2,2) - Cijr_124(2,2) - Dijr(2,2)*r3r2)
c FC %       
c FC %       Rr(3)=(PRr(2)-z32*(PRr(1)-z21*PRr(3))-z31*PRr(3))*iz33
c FC %           Rr(2)=(PRr(1)-z21*PRr(3) - z23*Rr(3))*iz22
c FC %           Rr(1)=(PRr(3)-  z12*Rr(2)-  z13*Rr(3))*iz11 
c FC %       
c FC %       Dijr(6,3) = Rr(1)
c FC %       Dijr(2,3) = Rr(2)
c FC %       Dijr(8,3) = Rr(3)
c FC %       
c FC % 
c FC %       
c FC %c      
c FC %c 56-5 7-58
c FC %c      
c FC %c      
c FC %       PRr(1) = (Cijr_134(2,2) - Cijr_234(2,2) - Dijr(3,2)*r1r0)
c FC %       PRr(2) = (Cijr_124(2,2) - Cijr_134(2,2) - Dijr(3,2)*r2r1)
c FC %       PRr(3) = (     -Cijr_124(2,2) - Dijr(3,2)*r3r2 - Dijr(13,3)*4.d0)
c FC %c       
c FC %       
c FC %       Rr(3)=(PRr(2)-z32*(PRr(1)-z21*PRr(3))-z31*PRr(3))*iz33
c FC %           Rr(2)=(PRr(1)-z21*PRr(3) - z23*Rr(3))*iz22
c FC %           Rr(1)=(PRr(3)-  z12*Rr(2)-  z13*Rr(3))*iz11 
c FC %       
c FC %       Dijr(7,3) = Rr(1)
c FC %       Dijr(9,3) = Rr(2)
c FC %       Dijr(3,3) = Rr(3)
c FC %c      
c FC %
c FC %c 44-4 5-46
c FC %       PRr(1) = (Cijr_134(1,2) + Cijr_234(1,1) - Dijr(4,2)*r1r0 
c FC %     &      -Dijr(12,3)*2.d0)
c FC %       PRr(2) = (Cijr_124(3,2) - Cijr_134(1,2) - Dijr(4,2)*r2r1
c FC %     &      - Dijr(11,3)*2.d0)
c FC %       PRr(3) = (Cijr_123(3,2) - Cijr_124(3,2)  -Dijr(4,2)*r3r2)
c FC %c      
c FC %       
c FC %       Rr(3)=(PRr(2)-z32*(PRr(1)-z21*PRr(3))-z31*PRr(3))*iz33
c FC %c          Rr(2)=(PRr(1)-z21*PRr(3) - z23*Rr(3))*iz22
c FC %c          Rr(1)=(PRr(3)-  z12*Rr(2)-  z13*Rr(3))*iz11 
c FC %       
c FC %C       Dijr(4,3) = Rr(1)
c FC %C       Dijr(6,3) = Rr(2)
c FC %       Dijr(10,3) = Rr(3)
c FC %c      
c FC %c      
c FC % 
c FC %
c FC %c to d efine D00ij and D0000 functions!
c FC %c In P V notation we have:
c FC %c Dij( 1,7)=D416, Dij(1,8)=D417, Dij(1,9)=D418
c FC %c Dij( 1,10)=D419, Dij(1,11)=D420, Dij(1,12)=D421, D(1,13)=D422
c FC %c      
c FC %       Iv=1.d0/6.d0
c FC %c      
c FC %c 
c FC %
c FC %      Dijr(7,1)=(C0r_234+2*m0sq*Dijr(1,2)+Dijr(1,3)*r1r0+Dijr(4,3)*r2r1+Dijr(5,3)*r3r2 )*Iv
c FC %c      
c FC %      Dijr(8,1)=(Cijr_234(1,2)+2*m0sq*Dijr(2,2)+Dijr(6,3)*r1r0+Dijr(2,3)*r2r1
c FC %     - +Dijr(8,3)*r3r2)*Iv
c FC %c
c FC %      Dijr(9,1)=(Cijr_234(2,2)+2*m0sq*Dijr(3,2)+Dijr(7,3)*r1r0+Dijr(9,3)*r2r1
c FC %     - +Dijr(3,3)*r3r2)*Iv
c FC %c
c FC %      Dijr(10,1)=(-Cijr_234(1,1)+2*m0sq*Dijr(4,2)+Dijr(4,3)*r1r0+Dijr(6,3)*r2r1
c FC %     - +Dijr(10,3)*r3r2)*Iv
c FC %c      
c FC %      Dijr(11,1)=(-Cijr_234(2,1)+2*m0sq*Dijr(5,2)+Dijr(5,3)*r1r0+Dijr(10,3)*r2r1
c FC %     - +Dijr(7,3)*r3r2)*Iv
c FC %c
c FC %      Dijr(12,1)=(Cijr_234(3,2)+2*m0sq*Dijr(6,2)+Dijr(10,3)*r1r0+Dijr(8,3)*r2r1
c FC %     - +Dijr(9,3)*r3r2)*Iv
c FC %c
c FC %      Dijr(13,1)=(Cijr_234(4,2)+2*m0sq*Dijr(7,2)+Dijr(11,3)*r1r0+Dijr(12,3)*r2r1
c FC %     -      + Dijr(13,3)*r3r2+Iv)*Iv        
c FC %c
c FC %
c              RETURN
       
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC        
         		else
c        			iorder(1)=3
c        			iorder(2)=2
c        			iorder(3)=1
c      
         			z23=2d0*p2p3
         			z33=2d0*p1p3
c        			det=-det
       
         		  iz22=1d0/z22	
           	  	  z32=z32*iz22
           	  	  z23=z23-z21*z13
           	  	  z33=z33-z31*z13-z32*z23
           	  	  iz33=1d0/z33
                   det=-z11*z22*z33
c         	 print*, "det8",det
       	      

       PRr(1) = (C0r_134 - C0r_234 - D0r*r1r0   )
       PRr(2) = (C0r_124 - C0r_134 - D0r*r2r1 )
       PRr(3) = (C0r_123 - C0r_124 - D0r*r3r2 )


           Rr(3)=(PRr(1)-z32*(PRr(2)-z21*PRr(3))-z31*PRr(3))*iz33
           Rr(2)=(PRr(2)-z21*PRr(3) - z23*Rr(3))*iz22
           Rr(1)=(PRr(3)-  z12*Rr(2)-  z13*Rr(3))*iz11 	
       

       Dijr(1,1) = Rr(1)
       Dijr(2,1) = Rr(2)
       Dijr(3,1) = Rr(3)
c      

CCCCCC CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
c D00  
CCCCCC CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      Dijr(7,2) = (C0r_234+ 2d0*m0sq*D0r + Dijr(1,1)*r1r0 + Dijr(2,1)*r2r1
     &     +Dijr(3,1)*r3r2)*0.5d0


CCCCCC CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
CCCCCC CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
c FC %c 30-3 1-32
c FC %       
c FC %c      
c FC %       PRr(1) = (Cijr_134(1,1) + C0r_234 - Dijr(1,1)*r1r0 - Dijr(7,2)*2.d0)
c FC %       PRr(2) = (Cijr_124(1,1) - Cijr_134(1,1) - Dijr(1,1)*r2r1        )
c FC %       PRr(3) = (Cijr_123(1,1) - Cijr_124(1,1) - Dijr(1,1)*r3r2        )     
c FC %c       
c FC %           Rr(3)=(PRr(1)-z32*(PRr(2)-z21*PRr(3))-z31*PRr(3))*iz33
c FC %           Rr(2)=(PRr(2)-z21*PRr(3) - z23*Rr(3))*iz22
c FC %           Rr(1)=(PRr(3)-  z12*Rr(2)-  z13*Rr(3))*iz11 
c FC %       
c FC %       Dijr(1,2) = Rr(1)
c FC %       Dijr(4,2) = Rr(2)
c FC %       Dijr(5,2) = Rr(3)
c FC %c       
c FC %
c FC %       
c FC %c 33-3 4-35
c FC %       PRr(1) = (Cijr_134(1,1) -Cijr_234(1,1) - Dijr(2,1)*r1r0          )
c FC %       PRr(2)=(Cijr_124(2,1)-Cijr_134(1,1)-Dijr(2,1)*r2r1-Dijr(7,2)*2.d0)
c FC %       PRr(3) = (Cijr_123(2,1) -Cijr_124(2,1) - Dijr(2,1)*r3r2        )
c FC %       
c FC %       Rr(3)=(PRr(1)-z32*(PRr(2)-z21*PRr(3))-z31*PRr(3))*iz33
c FC %           Rr(2)=(PRr(2)-z21*PRr(3) - z23*Rr(3))*iz22
c FC %c          Rr(1)=(PRr(3)-  z12*Rr(2)-  z13*Rr(3))*iz11 
c FC %       
c FC %c       Dijr(4,2) = Rr(1)
c FC %       Dijr(2,2) = Rr(2)
c FC %       Dijr(6,2) = Rr(3)
c FC %       
c FC %
c FC %       
c FC %c 36-3 7-38
c FC %       
c FC %       PRr(1) = (Cijr_134(2,1) - Cijr_234(2,1) - Dijr(3,1)*r1r0      )
c FC %       PRr(2) = (Cijr_124(2,1) - Cijr_134(2,1) - Dijr(3,1)*r2r1    )
c FC %       PRr(3) = (       - Cijr_124(2,1) - Dijr(3,1)*r3r2 -Dijr(7,2)*2.d0)
c FC %       
c FC %       Rr(3)=(PRr(1)-z32*(PRr(2)-z21*PRr(3))-z31*PRr(3))*iz33
c FC %c          Rr(2)=(PRr(2)-z21*PRr(3) - z23*Rr(3))*iz22
c FC %c          Rr(1)=(PRr(3)-  z12*Rr(2)-  z13*Rr(3))*iz11 
c FC %       
c FC %c       Dijr(5,2) = Rr(1)
c FC %c       Dijr(6,2) = Rr(2)
c FC %       Dijr(3,2) = Rr(3)
c FC %       
c FC % 
c FC %       
c FC %       
c FC %CCCCCC CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
c FC %CCCCCC CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
c FC %
c FC %       Dijr(11,3)=0.25d0*(-C0r_234+2d0*m0sq*Dijr(1,1)+Dijr(1,2)*r1r0+Dijr(4,2)*r2r1
c FC %     &       +Dijr(5,2)*r3r2)
c FC %       
c FC %       Dijr(12,3)=0.25d0*(Cijr_234(1,1)+2d0*m0sq*Dijr(2,1)+Dijr(4,2)*r1r0+Dijr(2,2)*r2r1
c FC %     &                   +Dijr(6,2)*r3r2)
c FC %       Dijr(13,3)=0.25d0*(Cijr_234(2,1)+2d0*m0sq*Dijr(3,1)+Dijr(5,2)*r1r0 + Dijr(6,2)*r2r1
c FC %     &                    + Dijr(3,2)*r3r2)
c FC %       
c FC % 
c FC %       
c FC %      
c FC %       
c FC %       
c FC %c 41-4 2-43
c FC %       
c FC %       PRr(1) =(Cijr_134(1,2) - C0r_234 - Dijr(1,2)*r1r0  - Dijr(11,3)*4.d0)
c FC %       PRr(2) =(Cijr_124(1,2) - Cijr_134(1,2) - Dijr(1,2)*r2r1       )
c FC %       PRr(3) =(Cijr_123(1,2) - Cijr_124(1,2) - Dijr(1,2)*r3r2       )
c FC %       
c FC %       Rr(3)=(PRr(1)-z32*(PRr(2)-z21*PRr(3))-z31*PRr(3))*iz33
c FC %           Rr(2)=(PRr(2)-z21*PRr(3) - z23*Rr(3))*iz22
c FC %           Rr(1)=(PRr(3)-  z12*Rr(2)-  z13*Rr(3))*iz11 
c FC %       
c FC %       Dijr(1,3) = Rr(1)
c FC %       Dijr(4,3) = Rr(2)
c FC %       Dijr(5,3) = Rr(3)
c FC %       
c FC %  
c FC %       
c FC %c 50-5 1-52
c FC %       
c FC %       PRr(1) =  (Cijr_134(1,2) - Cijr_234(1,2) - Dijr(2,2)*r1r0)
c FC %       PRr(2) =  (Cijr_124(2,2) - Cijr_134(1,2) - Dijr(2,2)*r2r1
c FC %     &       -Dijr(12,3)*4.d0)
c FC %       PRr(3) =  (Cijr_123(2,2) - Cijr_124(2,2) - Dijr(2,2)*r3r2)
c FC %       
c FC %       Rr(3)=(PRr(1)-z32*(PRr(2)-z21*PRr(3))-z31*PRr(3))*iz33
c FC %           Rr(2)=(PRr(2)-z21*PRr(3) - z23*Rr(3))*iz22
c FC %           Rr(1)=(PRr(3)-  z12*Rr(2)-  z13*Rr(3))*iz11 
c FC %       
c FC %       Dijr(6,3) = Rr(1)
c FC %       Dijr(2,3) = Rr(2)
c FC %       Dijr(8,3) = Rr(3)
c FC %       
c FC %
c FC %       
c FC %c      
c FC %c 56-5 7-58
c FC %c      
c FC %c      
c FC %       PRr(1) = (Cijr_134(2,2) - Cijr_234(2,2) - Dijr(3,2)*r1r0)
c FC %       PRr(2) = (Cijr_124(2,2) - Cijr_134(2,2) - Dijr(3,2)*r2r1)
c FC %       PRr(3) = (     -Cijr_124(2,2) - Dijr(3,2)*r3r2 - Dijr(13,3)*4.d0)
c FC %c       
c FC %       
c FC %       Rr(3)=(PRr(1)-z32*(PRr(2)-z21*PRr(3))-z31*PRr(3))*iz33
c FC %           Rr(2)=(PRr(2)-z21*PRr(3) - z23*Rr(3))*iz22
c FC %           Rr(1)=(PRr(3)-  z12*Rr(2)-  z13*Rr(3))*iz11 
c FC %       
c FC %       Dijr(7,3) = Rr(1)
c FC %       Dijr(9,3) = Rr(2)
c FC %       Dijr(3,3) = Rr(3)
c FC %c      
c FC %     
c FC %c 44-4 5-46
c FC %       PRr(1) = (Cijr_134(1,2) + Cijr_234(1,1) - Dijr(4,2)*r1r0 
c FC %     &      -Dijr(12,3)*2.d0)
c FC %       PRr(2) = (Cijr_124(3,2) - Cijr_134(1,2) - Dijr(4,2)*r2r1
c FC %     &      - Dijr(11,3)*2.d0)
c FC %       PRr(3) = (Cijr_123(3,2) - Cijr_124(3,2)  -Dijr(4,2)*r3r2)
c FC %c      
c FC %       
c FC %       Rr(3)=(PRr(1)-z32*(PRr(2)-z21*PRr(3))-z31*PRr(3))*iz33
c FC %c          Rr(2)=(PRr(2)-z21*PRr(3) - z23*Rr(3))*iz22
c FC %c          Rr(1)=(PRr(3)-  z12*Rr(2)-  z13*Rr(3))*iz11 
c FC %       
c FC %C       Dijr(4,3) = Rr(1)
c FC %C       Dijr(6,3) = Rr(2)
c FC %       Dijr(10,3) = Rr(3)
c FC %c      
c FC %c      
c FC %     
c FC %
c FC %
c FC %c to d efine D00ij and D0000 functions!
c FC %c In P V notation we have:
c FC %c Dij( 1,7)=D416, Dij(1,8)=D417, Dij(1,9)=D418
c FC %c Dij( 1,10)=D419, Dij(1,11)=D420, Dij(1,12)=D421, D(1,13)=D422
c FC %c      
c FC %       Iv=1.d0/6.d0
c FC %c      
c FC %c 
c FC %
c FC %      Dijr(7,1)=(C0r_234+2*m0sq*Dijr(1,2)+Dijr(1,3)*r1r0+Dijr(4,3)*r2r1+Dijr(5,3)*r3r2 )*Iv
c FC %c      
c FC %      Dijr(8,1)=(Cijr_234(1,2)+2*m0sq*Dijr(2,2)+Dijr(6,3)*r1r0+Dijr(2,3)*r2r1
c FC %     - +Dijr(8,3)*r3r2)*Iv
c FC %c
c FC %      Dijr(9,1)=(Cijr_234(2,2)+2*m0sq*Dijr(3,2)+Dijr(7,3)*r1r0+Dijr(9,3)*r2r1
c FC %     - +Dijr(3,3)*r3r2)*Iv
c FC %c
c FC %      Dijr(10,1)=(-Cijr_234(1,1)+2*m0sq*Dijr(4,2)+Dijr(4,3)*r1r0+Dijr(6,3)*r2r1
c FC %     - +Dijr(10,3)*r3r2)*Iv
c FC %c      
c FC %      Dijr(11,1)=(-Cijr_234(2,1)+2*m0sq*Dijr(5,2)+Dijr(5,3)*r1r0+Dijr(10,3)*r2r1
c FC %     - +Dijr(7,3)*r3r2)*Iv
c FC %c
c FC %      Dijr(12,1)=(Cijr_234(3,2)+2*m0sq*Dijr(6,2)+Dijr(10,3)*r1r0+Dijr(8,3)*r2r1
c FC %     - +Dijr(9,3)*r3r2)*Iv
c FC %c
c FC %      Dijr(13,1)=(Cijr_234(4,2)+2*m0sq*Dijr(7,2)+Dijr(11,3)*r1r0+Dijr(12,3)*r2r1
c FC %     -      + Dijr(13,3)*r3r2+Iv)*Iv        
c FC %c
c FC %     


c       RETURN

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC   
			endif
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC   
     	    endif
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC   
     	endif


       do j=1,3
        do i=1,indexD(j)
          Dijrr(i,j)=Dble(Dijr(i,j))  
          DijI(i,j)=DImag(Dijr(i,j))  
        enddo   
       enddo 


        return
        end

