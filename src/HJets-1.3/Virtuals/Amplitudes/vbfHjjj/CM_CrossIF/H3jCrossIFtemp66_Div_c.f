c ************************************************************************************
c ************************************************************************************
c Author: Francisco Campanario
c E-mail: francam@particle.uni-karlsruhe.de
c Date: 18/01/2010
c Modified:8/7/2011
c ************************************************************************************
c determine the  finite virtual corrections of 
c                     p_5              
c                     H              
c                     H              
c                     H              
c psi(p6) ---->-----$$$$$--V----$$$--->---   psi(p4)
c                   |              |              
c                   |              | 
c                   |              |             
c                   |              |$$$$gluon$$$$$$$$  p3, mu_p3 
c                   |              |             
c                   |              |              
c                   |              |             
c barpsi(p1)-->-$$$gluon$$$$$$------>---   bar_psi(p2)
c                                               
c Note: To make it shorter in the promgram: mu_p3,...->mup3,... 
c Notation of External momenta: p1+p2+p3+p4+p5+p6=0 
c mu_p3, should be think as external current 
c alpha is the helicity of the initial spinor 
c musq is the renormalization scale energy  
c comp: integer value.The first time called with p1...p5, comp=1
c ATTENTION: ONLY!!!If you have to call the subroutine consecutively with the same arguments
c(p1,p2,p3,p4,p5). Then, comp=-1 (it safes 4000 lines of code) 
c This applies when you have for examples the same diagram for an off-shell photon
c and a Z boson. The differences are in the coupling and  the part that depends on the
c polarization vector that are calculated at the end of this program.
c************************************************************************************
c************************************************************************************
       subroutine H3jCrossIF66_Div_c(M1sq,p1,p2,p3,p4,p5,p6,psi_p1,barpsi_p6,ps
     &   i_p2,barpsi_p4,mup3,musq,comp,result,resultborn,Div)
c************************************************************************************
c************************************************************************************
c    Declaration of variables 
c************************************************************************************
c************************************************************************************
c************************************************************************************
c************************************************************************************
c************************************************************************************
c************************************************************************************
c************************************************************************************
c************************************************************************************
       IMPLICIT NONE
      Complex*16 SMB(95),Fa(40),F(203)
       Complex*16 SMB0(16),SMB11(0:3),SMB12(0:3),SMB13(0:3),SMB14(0:3)
      Complex*16 SMB1(38)
      Real*8 FI(203),FR(203)
      Complex*16 psi_p1_P(2),barpsi_p6_P(2),psi_p2_P(2),barpsi_p4_P(2) 
       Complex*16 psi_p1_M(2),barpsi_p6_M(2),psi_p2_M(2),barpsi_p4_M(2)
       Complex*16 psi_p1(4),barpsi_p6(4),psi_p2(4),barpsi_p4(4)
       Real*8 P(61)
       Complex*16  SC1c,SC1r, SC3crr,SC3ccr,SC3crc,SC3ccc,SC3rrc,SC3rrr, SC5ccrrr,SC5ccrrc
       EXTERNAL    SC1c,SC1r, SC3crr,SC3ccr,SC3crc,SC3ccc,SC3rrc,SC3rrr, SC5ccrrr,SC5ccrrc 
       Real*8 delta
       External  delta
       Complex*16 ten2cc
       External  ten2cc
       Complex*16 v1(0:3),v2(0:3),v3(0:3) 
       Integer alpha_1,alpha_2,i,j,k,comp
c************************************************************************************
c************************************************************************************
c************************************************************************************
c************************************************************************************
c    Declaration of variables 
c************************************************************************************
c************************************************************************************
c************************************************************************************
c************************************************************************************
      Real*8   p1(0:3),p2(0:3),p3(0:3),p4(0:3),p5(0:3),p6(0:3)
       Complex*16   mup3(0:3),mup4(0:3)
       Complex*16   p1mup3, p1mup4, p2mup3, p2mup4, p3mup3, p3mup4, 
     -          p4mup3, p4mup4, p5mup3, p5mup4, p6mup3, p6mup4
       Complex*16   mup3mup4
       Real*8 dotrr
       Complex*16 A0finG,B0finG,C0finG,D0finG,E0finG,F0finG
       EXTERNAL dotrr,A0finG,B0finG,C0finG,D0finG,E0finG,F0finG
        Real*8   p1sq, p1p2, p1p3, p1p4, p1p5,p1p6 
       Real*8   p2sq, p2p3, p2p4, p2p5 
       Real*8   p3sq, p3p4, p3p5 
       Real*8   p4sq, p4p5 
       Real*8   p5sq, p5p6 
       Real*8   p6sq 
       Real*8   s12, s13, s14, s15, s16 
       Real*8   s23, s24, s25,s26 
       Real*8   s34, s35,s36 
       Real*8   s45,s46 
       Real*8   s56 
       Real*8   s123,s234,s345   
       Real*8   p2p6,p3p6,p4p6  
       Complex*16  A01,A02,A03,A04,A05,A06  
       Real*8  A01R,A02R,A03R,A04R,A05R,A06R 
       Real*8  A01I,A02I,A03I,A04I,A05I,A06I   
       Complex*16  B012,B013,B014,B015,B016 
       Complex*16  B023,B024,B025,B026 
       Complex*16  B034,B035,B036 
       Complex*16  B045,B046 
       Complex*16  B056   
       Real*8  B012R,B013R,B014R,B015R,B016R 
       Real*8  B023R,B024R,B025R,B026R 
       Real*8  B034R,B035R,B036R 
       Real*8  B045R,B046R 
       Real*8  B056R  
       Real*8  B012I,B013I,B014I,B015I,B016I 
       Real*8  B023I,B024I,B025I,B026I 
       Real*8  B034I,B035I,B036I 
       Real*8  B045I,B046I 
       Real*8  B056I     
       Real*8  Bij12R,Bij13R,Bij14R,Bij15R,Bij16R 
       Real*8  Bij23R,Bij24R,Bij25R,Bij26R 
       Real*8  Bij34R,Bij35R,Bij36R 
       Real*8  Bij45R,Bij46R 
       Real*8  Bij56R  
       Real*8  Bij12I,Bij13I,Bij14I,Bij15I,Bij16I 
       Real*8  Bij23I,Bij24I,Bij25I,Bij26I 
       Real*8  Bij34I,Bij35I,Bij36I 
       Real*8  Bij45I,Bij46I 
       Real*8  Bij56I       
       Complex*16 C0123,C0124,C0125,C0126 
       Complex*16 C0134,C0135,C0136 
       Complex*16 C0145,C0146 
       Complex*16 C0156 
       Complex*16 C0234,C0235,C0236 
       Complex*16 C0245,C0246 
       Complex*16 C0256 
       Complex*16 C0345,C0346 
       Complex*16 C0356 
       Complex*16 C0456   
       Real*8 C0123R,C0124R,C0125R,C0126R 
       Real*8 C0134R,C0135R,C0136R 
       Real*8 C0145R,C0146R 
       Real*8 C0156R 
       Real*8 C0234R,C0235R,C0236R 
       Real*8 C0245R,C0246R 
       Real*8 C0256R 
       Real*8 C0345R,C0346R 
       Real*8 C0356R 
       Real*8 C0456R    
       Real*8 C0123I,C0124I,C0125I,C0126I 
       Real*8 C0134I,C0135I,C0136I 
       Real*8 C0145I,C0146I 
       Real*8 C0156I 
       Real*8 C0234I,C0235I,C0236I 
       Real*8 C0245I,C0246I 
       Real*8 C0256I 
       Real*8 C0345I,C0346I 
       Real*8 C0356I 
       Real*8 C0456I      
       Real*8 C123R(4,2),C124R(4,2),C125R(4,2),C126R(4,2) 
       Real*8 C134R(4,2),C135R(4,2),C136R(4,2) 
       Real*8 C145R(4,2),C146R(4,2) 
       Real*8 C156R(4,2) 
       Real*8 C234R(4,2),C235R(4,2),C236R(4,2) 
       Real*8 C245R(4,2),C246R(4,2) 
       Real*8 C256R(4,2) 
       Real*8 C345R(4,2),C346R(4,2) 
       Real*8 C356R(4,2) 
       Real*8 C456R(4,2)  
       Real*8 C123I(4,2),C124I(4,2),C125I(4,2),C126I(4,2) 
       Real*8 C134I(4,2),C135I(4,2),C136I(4,2) 
       Real*8 C145I(4,2),C146I(4,2) 
       Real*8 C156I(4,2) 
       Real*8 C234I(4,2),C235I(4,2),C236I(4,2) 
       Real*8 C245I(4,2),C246I(4,2) 
       Real*8 C256I(4,2) 
       Real*8 C345I(4,2),C346I(4,2) 
       Real*8 C356I(4,2) 
       Real*8 C456I(4,2)  
       Complex*16  D01234,D01235,D01236 
       Complex*16 D01245,D01246 
       Complex*16 D01256 
       Complex*16 D01345,D01346 
       Complex*16 D01356 
       Complex*16 D01456 
       Complex*16 D02345,D02346 
       Complex*16 D02356 
       Complex*16 D02456 
       Complex*16 D03456    
       Real*8 D01234R,D01235R,D01236R 
       Real*8 D01245R,D01246R 
       Real*8 D01256R 
       Real*8 D01345R,D01346R 
       Real*8 D01356R 
       Real*8 D01456R 
       Real*8 D02345R,D02346R 
       Real*8 D02356R 
       Real*8 D02456R 
       Real*8 D03456R  
       Real*8 D01234I,D01235I,D01236I 
       Real*8 D01245I,D01246I 
       Real*8 D01256I 
       Real*8 D01345I,D01346I 
       Real*8 D01356I 
       Real*8 D01456I 
       Real*8 D02345I,D02346I 
       Real*8 D02356I 
       Real*8 D02456I 
       Real*8 D03456I   
       Real*8 D1234R(13,3),D1235R(13,3),D1236R(13,3) 
       Real*8 D1245R(13,3),D1246R(13,3) 
       Real*8 D1256R(13,3) 
       Real*8 D1345R(13,3),D1346R(13,3) 
       Real*8 D1356R(13,3) 
       Real*8 D1456R(13,3) 
       Real*8 D2345R(13,3),D2346R(13,3) 
       Real*8 D2356R(13,3) 
       Real*8 D2456R(13,3) 
       Real*8 D3456R(13,3)   
       Real*8 D1234I(13,3),D1235I(13,3),D1236I(13,3) 
       Real*8 D1245I(13,3),D1246I(13,3) 
       Real*8 D1256I(13,3) 
       Real*8 D1345I(13,3),D1346I(13,3) 
       Real*8 D1356I(13,3) 
       Real*8 D1456I(13,3) 
       Real*8 D2345I(13,3),D2346I(13,3) 
       Real*8 D2356I(13,3) 
       Real*8 D2456I(13,3) 
       Real*8 D3456I(13,3)   
       Complex*16 E012345,E012346,E012356,E012456,E013456,E023456 
       Real*8 E012345R,E012346R,E012356R,E012456R,E013456R,E023456R 
       Real*8 E012345I,E012346I,E012356I,E012456I,E013456I,E023456I  
       Real*8 E12345R(46,4),E12346R(46,4),E12356R(46,4) 
       Real*8 E12456R(46,4),E13456R(46,4),E23456R(46,4)  
       Real*8 E12345I(46,4),E12346I(46,4),E12356I(46,4) 
       Real*8 E12456I(46,4),E13456I(46,4),E23456I(46,4)  
       Complex*16 F0123456 
       Real*8 F0123456R 
       Real*8 F0123456I 
       Real*8 F123456R(166,5) 
       Real*8 F123456I(166,5) 
       Complex*16 F123456(166,5) 
       Logical PrintB,PrintC,PrintD,PrintE,PrintF 
      Complex*16 dotrc,dotcc,result(3),resultn,resultborn
       Real*8 musq
      EXTERNAL   dotrc,dotcc
      Integer alpha
       COMMON/H3jCrossIFFaFunctions/Fa
       COMMON/H3jCrossIFFhlFunctions/F
      Save/H3jCrossIFFhlFunctions/
       COMMON/H3jCrossIFPFunctions/P
      Save/H3jCrossIFPFunctions/
       COMMON/H3jCrossIFInvariants/p1sq,p2sq,p3sq,p4sq,p5sq,p6sq,s12,s
     &   23,s34,s45,s56,s16,s123,s234,s345
       COMMON/H3jCrossIFFVALFunctions/F0123456R,F123456R,F0123456I,F12
     &   3456I
       COMMON/H3jCrossIFEVALFunctions/ E012345R,E12345R,E012345I,E1234
     &   5I, E012346R,E12346R,E012346I,E12346I, E012356R,E12356R,E012
     &   356I,E12356I, E012456R,E12456R,E012456I,E12456I, E013456R,E1
     &   3456R,E013456I,E13456I, E023456R,E23456R,E023456I,E23456I
       COMMON/H3jCrossIFDVALFunctions/ D01234R,D1234R,D01234I,D1234I, 
     &   D01235R,D1235R,D01235I,D1235I, D01236R,D1236R,D01236I,D1236I
     &   , D01245R,D1245R,D01245I,D1245I, D01246R,D1246R,D01246I,D124
     &   6I, D01256R,D1256R,D01256I,D1256I, D01345R,D1345R,D01345I,D1
     &   345I, D01346R,D1346R,D01346I,D1346I, D01356R,D1356R,D01356I,
     &   D1356I, D01456R,D1456R,D01456I,D1456I, D02345R,D2345R,D02345
     &   I,D2345I, D02346R,D2346R,D02346I,D2346I, D02356R,D2356R,D023
     &   56I,D2356I, D02456R,D2456R,D02456I,D2456I, D03456R,D3456R,D0
     &   3456I,D3456I
       COMMON/H3jCrossIFCVALFunctions/ C0123R,C123R,C0123I,C123I, C012
     &   4R,C124R,C0124I,C124I, C0125R,C125R,C0125I,C125I, C0126R,C12
     &   6R,C0126I,C126I, C0134R,C134R,C0134I,C134I, C0135R,C135R,C01
     &   35I,C135I, C0136R,C136R,C0136I,C136I, C0145R,C145R,C0145I,C1
     &   45I, C0146R,C146R,C0146I,C146I, C0156R,C156R,C0156I,C156I, C
     &   0234R,C234R,C0234I,C234I, C0235R,C235R,C0235I,C235I, C0236R,
     &   C236R,C0236I,C236I, C0245R,C245R,C0245I,C245I, C0246R,C246R,
     &   C0246I,C246I, C0256R,C256R,C0256I,C256I, C0345R,C345R,C0345I
     &   ,C345I, C0346R,C346R,C0346I,C346I, C0356R,C356R,C0356I,C356I
     &   , C0456R,C456R,C0456I,C456I
       COMMON/H3jCrossIFBVALFunctions/ B012R,B012I, B013R,B013I, B014R
     &   ,B014I, B015R,B015I, B016R,B016I, B023R,B023I, B024R,B024I, 
     &   B025R,B025I, B026R,B026I, B034R,B034I, B035R,B035I, B036R,B0
     &   36I, B045R,B045I, B046R,B046I, B056R,B056I
       Integer ngluon, posgluon,Div
       complex*16 m1sq




       if(comp.gt.0) then
       call FgetCrossIF_Div_c(M1sq,p1,p2,p3,p4,p5,p6,musq,Div)



c************************************************************************************
c************************************************************************************
c************************************************************************************
c       Definition of the F,P functions:Independent of the currents    
c************************************************************************************
c************************************************************************************
c************************************************************************************
       call H3jCrossIFFFhl1(F(1))
       call H3jCrossIFFFhl2(F(41))
       call H3jCrossIFFFhl3(F(81))
       call H3jCrossIFFFhl4(F(121))
       call H3jCrossIFFFhl5(F(161))
c************************************************************************************
c************************************************************************************
c************************************************************************************
       endif  
c               PART THAT DEPENDS ON THE EXTERNAL CURRENT
c************************************************************************************
c************************************************************************************
c************************************************************************************
       p1mup3 = dotrc(p1,mup3)
       p1mup4 = dotrc(p1,mup4)
       p2mup3 = dotrc(p2,mup3)
       p2mup4 = dotrc(p2,mup4)
       p3mup3 = dotrc(p3,mup3)
       p3mup4 = dotrc(p3,mup4)
       p4mup3 = dotrc(p4,mup3)
       p4mup4 = dotrc(p4,mup4)
       p5mup3 = dotrc(p5,mup3)
       p5mup4 = dotrc(p5,mup4)
       p6mup3 = dotrc(p6,mup3)
       p6mup4 = dotrc(p6,mup4)
       mup3mup4 = dotcc(mup3,mup4)
c************** Calling the Fa functions**********************************************************************
c************************************************************************************
c************************************************************************************
       call H3jCrossIFFa1(p1mup3,p1mup4,p2mup3,p2mup4,p3mup3,p3mup4,p4
     &   mup3,p4mup4,p5mup3,p5mup4,p6mup3,p6mup4,mup3mup4,Fa(1))
       call H3jCrossIFFa2(p1mup3,p1mup4,p2mup3,p2mup4,p3mup3,p3mup4,p4
     &   mup3,p4mup4,p5mup3,p5mup4,p6mup3,p6mup4,mup3mup4,Fa(21))
c************************************************************************************
c************************************************************************************
c************************************************************************************
c       Definition of the Matrix Element  
c************************************************************************************
c************************************************************************************
c************************************************************************************
c************************************************************************************
c************************************************************************************
c************************************************************************************
c************************************************************************************
c************************************************************************************
       barpsi_p6_P(1)=barpsi_p6(1)
       barpsi_p6_P(2)=barpsi_p6(2)
       barpsi_p6_M(1)=barpsi_p6(3)
       barpsi_p6_M(2)=barpsi_p6(4)
       barpsi_p4_P(1)=barpsi_p4(1)
       barpsi_p4_P(2)=barpsi_p4(2)
       barpsi_p4_M(1)=barpsi_p4(3)
       barpsi_p4_M(2)=barpsi_p4(4)
       psi_p1_M(1)=psi_p1(1)
       psi_p1_M(2)=psi_p1(2)
       psi_p1_P(1)=psi_p1(3)
       psi_p1_P(2)=psi_p1(4)
       psi_p2_M(1)=psi_p2(1)
       psi_p2_M(2)=psi_p2(2)
       psi_p2_P(1)=psi_p2(3)
       psi_p2_P(2)=psi_p2(4)
c************************************************************************************
c************************************************************************************
       SMB0(1) = SC1c(barpsi_p4_P,mup3,psi_p2_P,+1)
       SMB0(2) = SC1r(barpsi_p6_P,p2,psi_p1_P,+1)
       SMB0(3) = SC1r(barpsi_p6_P,p4,psi_p1_P,+1)
       SMB0(4) = SC1r(barpsi_p6_P,p5,psi_p1_P,+1)
       SMB0(5) = SC3rrr(barpsi_p6_P,p5,p4,p2,psi_p1_P,+1)
       SMB0(6) = SC1r(barpsi_p4_P,p5,psi_p2_P,+1)
       SMB0(7) = SC1r(barpsi_p4_P,p1,psi_p2_P,+1)
       SMB0(8) = SC1c(barpsi_p6_P,mup3,psi_p1_P,+1)
       SMB0(9) = SC3crr(barpsi_p6_P,mup3,p4,p2,psi_p1_P,+1)
       SMB0(10) = SC3crr(barpsi_p6_P,mup3,p5,p2,psi_p1_P,+1)
       SMB0(11) = SC3crr(barpsi_p6_P,mup3,p5,p4,psi_p1_P,+1)
       SMB0(12) = SC1r(barpsi_p4_P,p6,psi_p2_P,+1)
       SMB0(13) = SC3rrr(barpsi_p4_P,p6,p5,p1,psi_p2_P,+1)
       SMB0(14) = SC3crr(barpsi_p4_P,mup3,p6,p1,psi_p2_P,+1)
       SMB0(15) = SC3crr(barpsi_p4_P,mup3,p5,p1,psi_p2_P,+1)
       SMB0(16) = SC3crr(barpsi_p4_P,mup3,p6,p5,psi_p2_P,+1)
c************************************************************************************
c************************************************************************************
      do i=0,3
        v1(0)=delta(i,0)
        v1(1)=delta(i,1)
        v1(2)=delta(i,2)
        v1(3)=delta(i,3)
       SMB11(i) = SC1c(barpsi_p4_P,v1,psi_p2_P,+1)
       SMB12(i) = SC1c(barpsi_p6_P,v1,psi_p1_P,+1)
       SMB13(i) = SC3crc(barpsi_p4_P,mup3,p5,v1,psi_p2_P,+1)
       SMB14(i) = SC3crc(barpsi_p6_P,mup3,p5,v1,psi_p1_P,+1)
      enddo
c************************************************************************************
c************************************************************************************
c************************************************************************************
c************************************************************************************
c************************************************************************************
c************************************************************************************
       SMB1(1) = SMB0(1)*SMB0(2)
       SMB1(2) = SMB0(1)*SMB0(3)
       SMB1(3) = SMB0(1)*SMB0(4)
       SMB1(4) = SMB0(1)*SMB0(5)
       SMB1(7) = SMB0(2)*SMB0(6)
       SMB1(8) = SMB0(7)*SMB0(8)
       SMB1(9) = SMB0(2)*SMB0(7)
       SMB1(10) = SMB0(3)*SMB0(7)
       SMB1(11) = SMB0(4)*SMB0(7)
       SMB1(12) = SMB0(7)*SMB0(9)
       SMB1(13) = SMB0(7)*SMB0(10)
       SMB1(14) = SMB0(7)*SMB0(11)
       SMB1(15) = SMB0(6)*SMB0(8)
       SMB1(16) = SMB0(3)*SMB0(6)
       SMB1(17) = SMB0(4)*SMB0(6)
       SMB1(18) = SMB0(6)*SMB0(9)
       SMB1(19) = SMB0(6)*SMB0(10)
       SMB1(20) = SMB0(6)*SMB0(11)
       SMB1(21) = SMB0(8)*SMB0(12)
       SMB1(22) = SMB0(2)*SMB0(12)
       SMB1(23) = SMB0(3)*SMB0(12)
       SMB1(24) = SMB0(4)*SMB0(12)
       SMB1(25) = SMB0(9)*SMB0(12)
       SMB1(26) = SMB0(10)*SMB0(12)
       SMB1(27) = SMB0(11)*SMB0(12)
       SMB1(28) = SMB0(8)*SMB0(13)
       SMB1(29) = SMB0(4)*SMB0(14)
       SMB1(30) = SMB0(2)*SMB0(15)
       SMB1(31) = SMB0(2)*SMB0(14)
       SMB1(32) = SMB0(2)*SMB0(16)
       SMB1(33) = SMB0(3)*SMB0(15)
       SMB1(34) = SMB0(3)*SMB0(14)
       SMB1(35) = SMB0(3)*SMB0(16)
       SMB1(36) = SMB0(4)*SMB0(15)
       SMB1(37) = SMB0(4)*SMB0(16)
c************************************************************************************
       SMB1(5) = dotcc(SMB11,SMB12)
       SMB1(6) = dotcc(SMB12,SMB13)
       SMB1(38) = dotcc(SMB11,SMB14)
c************************************************************************************
c************************************************************************************
       SMB(1) = SMB1(8)
       SMB(2) = SMB1(15)
       SMB(3) = SMB1(21)
       SMB(4) = SMB1(2)
       SMB(5) = SMB1(3)
       SMB(6) = SMB1(1)
       SMB(7) = SMB1(5)
       SMB(8) = 2*p1mup3*SMB1(5)-2*SMB1(8)
       SMB(9) = SMB1(6)
       SMB(10) = 0
       SMB(11) = SMB1(9)
       SMB(12) = SMB1(10)
       SMB(13) = SMB1(11)
       SMB(14) = SMB1(7)
       SMB(15) = SMB1(16)
       SMB(16) = SMB1(17)
       SMB(17) = SMB1(22)
       SMB(18) = SMB1(23)
       SMB(19) = SMB1(24)
       SMB(20) = SMB1(30)
       SMB(21) = SMB1(31)
       SMB(22) = SMB1(33)
       SMB(23) = SMB1(34)
       SMB(24) = SMB1(36)
       SMB(25) = SMB1(29)
       SMB(26) = SMB1(28)
       SMB(27) = SMB1(35)
       SMB(28) = SMB1(37)
       SMB(29) = SMB1(32)
       SMB(30) = 0
       SMB(31) = 2*p5mup3*SMB1(5)-SMB1(6)
       SMB(32) = -2*SMB1(1)+2*p2mup3*SMB1(5)
       SMB(33) = SMB1(12)
       SMB(34) = SMB1(13)
       SMB(35) = SMB1(18)
       SMB(36) = SMB1(19)
       SMB(37) = SMB1(25)
       SMB(38) = SMB1(26)
       SMB(39) = -((s16-s234+s56)*SMB1(5))-2*SMB1(11)
       SMB(40) = s16*SMB1(5)
       SMB(41) = (-p5sq+s56)*SMB1(5)
       SMB(42) = 4*p1mup3*SMB1(5)
       SMB(43) = 2*(SMB1(6)+SMB1(38))
       SMB(44) = 2*((p5sq-s45)*SMB1(8)+SMB1(14))
       SMB(45) = -2*p5mup3*(s16-s234+s56)*SMB1(5)+(s16-s234+s56)*SMB1(
     &   6)-2*p5sq*SMB1(8)
       SMB(46) = 4*p6mup3*SMB1(5)
       SMB(47) = 2*(-p5sq-s123+s45+s56)*SMB1(8)+4*p6mup3*SMB1(10)
       SMB(48) = 2*p5mup3*s16*SMB1(5)-s16*SMB1(6)+2*(p5sq-s56)*SMB1(8)
     &   +4*p6mup3*SMB1(11)
       SMB(49) = (2*p6mup3*(-p5sq+s45)-2*p5mup3*(p5sq+s123-s45-s56))*S
     &   MB1(5)+(p5sq+s123-s45-s56)*SMB1(6)
       SMB(50) = 2*p5sq*p6mup3*SMB1(5)
       SMB(51) = 2*((s16-s234+s56)*SMB1(1)-p2mup3*(s16-s234+s56)*SMB1(
     &   5)-(s16-s234+s34-s345)*SMB1(8)+SMB1(13))
       SMB(52) = -2*s16*SMB1(1)+2*p2mup3*s16*SMB1(5)+2*(s12+s16-s345)*
     &   SMB1(8)+4*p6mup3*SMB1(9)
       SMB(53) = 2*(p5sq-s56)*SMB1(1)+(2*p5mup3*(s12+s16-s345)+2*p2mup
     &   3*(-p5sq+s56))*SMB1(5)-(s12+s16-s345)*SMB1(6)+4*p6mup3*SMB1(
     &   7)
       SMB(54) = 4*p2mup3*SMB1(5)
       SMB(55) = 8*p5mup3*SMB1(5)-2*(SMB1(6)+SMB1(38))
       SMB(56) = 4*p4mup3*SMB1(5)
       SMB(57) = (s16-s234+s34-s345)*SMB1(5)-2*SMB1(7)
       SMB(58) = 2*((s16-s234+s56)*SMB1(1)+p1mup3*(s16-s234+s34-s345)*
     &   SMB1(5)+(-s16+s234-s34+s345)*SMB1(8)+SMB1(30))
       SMB(59) = -2*p5sq*SMB1(1)+(s16-s234+s34-s345)*(2*p5mup3*SMB1(5)
     &   -SMB1(38))
       SMB(60) = -2*SMB1(32)
       SMB(61) = (p3sq-s23+s234-s34)*SMB1(5)
       SMB(62) = 2*((s123-s23+s234-s56)*SMB1(1)+p1mup3*(p3sq-s23+s234-
     &   s34)*SMB1(5)+(-p3sq+s23-s234+s34)*SMB1(8)+2*p4mup3*SMB1(9))
       SMB(63) = 2*(p5sq-s45)*SMB1(1)+2*p5mup3*(p3sq-s23+s234-s34)*SMB
     &   1(5)+4*p4mup3*SMB1(7)+(-p3sq+s23-s234+s34)*SMB1(38)
       SMB(64) = 2*(-p5sq-s123+s45+s56)*SMB1(1)+4*p4mup3*SMB1(22)
       SMB(65) = SMB1(4)
       SMB(66) = 0
       SMB(67) = SMB1(14)
       SMB(68) = SMB1(20)
       SMB(69) = SMB1(27)
       SMB(70) = -2*p5mup3*(s123-s23+s234-s56)*SMB1(5)-4*p4mup3*SMB1(1
     &   1)+(s123-s23+s234-s56)*SMB1(38)
       SMB(71) = -2*(p4mup3*p5sq+p5mup3*(p5sq-s45))*SMB1(5)+(p5sq-s45)
     &   *SMB1(38)
       SMB(72) = 2*p4mup3*(p5sq-s56)*SMB1(5)+(p5sq+s123-s45-s56)*(2*p5
     &   mup3*SMB1(5)-SMB1(38))
       SMB(73) = -2*SMB1(29)
       SMB(74) = 16*SMB1(5)
       SMB(75) = 2*s12*SMB1(5)
       SMB(76) = -2*(s16-s234+s56)*SMB1(5)
       SMB(77) = 2*(s16-s234+s34-s345)*SMB1(5)
       SMB(78) = 4*p5sq*SMB1(5)
       SMB(79) = -2*(s12+s16-s345)*SMB1(5)
       SMB(80) = 2*(-p5sq+s56)*SMB1(5)
       SMB(81) = 32*(p1mup3*SMB1(5)-SMB1(8))
       SMB(82) = 32*p5mup3*SMB1(5)-16*SMB1(38)
       SMB(83) = 0
       SMB(84) = 4*p5mup3*s12*SMB1(5)-2*s12*SMB1(38)
       SMB(85) = 0
       SMB(86) = 0
       SMB(87) = 2*(s16-s234+s56)*(-2*p5mup3*SMB1(5)+SMB1(38))
       SMB(88) = 0
       SMB(89) = 0
       SMB(90) = 2*(-s123+s23-s234+s56)*SMB1(5)
       SMB(91) = 2*(-p5sq+s45)*SMB1(5)
       SMB(92) = 2*(p5sq+s123-s45-s56)*SMB1(5)
       SMB(93) = 2*(s123-s23+s234-s56)*(-2*p5mup3*SMB1(5)+SMB1(38))
       SMB(94) = 0
       SMB(95) = 0
c       Print*," SMB(1) ", SMB(1)
c       Print*," SMB(2) ", SMB(2)
c       Print*," SMB(3) ", SMB(3)
c       Print*," SMB(4) ", SMB(4)
c       Print*," SMB(5) ", SMB(5)
c       Print*," SMB(6) ", SMB(6)
c       Print*," SMB(7) ", SMB(7)
c       Print*," SMB(8) ", SMB(8)
c       Print*," SMB(9) ", SMB(9)
c       Print*," SMB(10) ", SMB(10)
c       Print*," SMB(11) ", SMB(11)
c       Print*," SMB(12) ", SMB(12)
c       Print*," SMB(13) ", SMB(13)
c       Print*," SMB(14) ", SMB(14)
c       Print*," SMB(15) ", SMB(15)
c       Print*," SMB(16) ", SMB(16)
c       Print*," SMB(17) ", SMB(17)
c       Print*," SMB(18) ", SMB(18)
c       Print*," SMB(19) ", SMB(19)
c       Print*," SMB(20) ", SMB(20)
c       Print*," SMB(21) ", SMB(21)
c       Print*," SMB(22) ", SMB(22)
c       Print*," SMB(23) ", SMB(23)
c       Print*," SMB(24) ", SMB(24)
c       Print*," SMB(25) ", SMB(25)
c       Print*," SMB(26) ", SMB(26)
c       Print*," SMB(27) ", SMB(27)
c       Print*," SMB(28) ", SMB(28)
c       Print*," SMB(29) ", SMB(29)
c       Print*," SMB(30) ", SMB(30)
c       Print*," SMB(31) ", SMB(31)
c       Print*," SMB(32) ", SMB(32)
c       Print*," SMB(33) ", SMB(33)
c       Print*," SMB(34) ", SMB(34)
c       Print*," SMB(35) ", SMB(35)
c       Print*," SMB(36) ", SMB(36)
c       Print*," SMB(37) ", SMB(37)
c       Print*," SMB(38) ", SMB(38)
c       Print*," SMB(39) ", SMB(39)
c       Print*," SMB(40) ", SMB(40)
c       Print*," SMB(41) ", SMB(41)
c       Print*," SMB(42) ", SMB(42)
c       Print*," SMB(43) ", SMB(43)
c       Print*," SMB(44) ", SMB(44)
c       Print*," SMB(45) ", SMB(45)
c       Print*," SMB(46) ", SMB(46)
c       Print*," SMB(47) ", SMB(47)
c       Print*," SMB(48) ", SMB(48)
c       Print*," SMB(49) ", SMB(49)
c       Print*," SMB(50) ", SMB(50)
c       Print*," SMB(51) ", SMB(51)
c       Print*," SMB(52) ", SMB(52)
c       Print*," SMB(53) ", SMB(53)
c       Print*," SMB(54) ", SMB(54)
c       Print*," SMB(55) ", SMB(55)
c       Print*," SMB(56) ", SMB(56)
c       Print*," SMB(57) ", SMB(57)
c       Print*," SMB(58) ", SMB(58)
c       Print*," SMB(59) ", SMB(59)
c       Print*," SMB(60) ", SMB(60)
c       Print*," SMB(61) ", SMB(61)
c       Print*," SMB(62) ", SMB(62)
c       Print*," SMB(63) ", SMB(63)
c       Print*," SMB(64) ", SMB(64)
c       Print*," SMB(65) ", SMB(65)
c       Print*," SMB(66) ", SMB(66)
c       Print*," SMB(67) ", SMB(67)
c       Print*," SMB(68) ", SMB(68)
c       Print*," SMB(69) ", SMB(69)
c       Print*," SMB(70) ", SMB(70)
c       Print*," SMB(71) ", SMB(71)
c       Print*," SMB(72) ", SMB(72)
c       Print*," SMB(73) ", SMB(73)
c       Print*," SMB(74) ", SMB(74)
c       Print*," SMB(75) ", SMB(75)
c       Print*," SMB(76) ", SMB(76)
c       Print*," SMB(77) ", SMB(77)
c       Print*," SMB(78) ", SMB(78)
c       Print*," SMB(79) ", SMB(79)
c       Print*," SMB(80) ", SMB(80)
c       Print*," SMB(81) ", SMB(81)
c       Print*," SMB(82) ", SMB(82)
c       Print*," SMB(83) ", SMB(83)
c       Print*," SMB(84) ", SMB(84)
c       Print*," SMB(85) ", SMB(85)
c       Print*," SMB(86) ", SMB(86)
c       Print*," SMB(87) ", SMB(87)
c       Print*," SMB(88) ", SMB(88)
c       Print*," SMB(89) ", SMB(89)
c       Print*," SMB(90) ", SMB(90)
c       Print*," SMB(91) ", SMB(91)
c       Print*," SMB(92) ", SMB(92)
c       Print*," SMB(93) ", SMB(93)
c       Print*," SMB(94) ", SMB(94)
c       Print*," SMB(95) ", SMB(95)
c************************************************************************************
c************************************************************************************
c************************************************************************************
c************************************************************************************
c       Amplitude                         
c************************************************************************************
c************************************************************************************
c************************************************************************************
c************************************************************************************


c The Finite virtual piece should be multiplied for (-1)  since 
c I have multiplied by (-I) to get the F's and k's without (I) factor
c . The factorization from the B_ij is Fact=(I/(4\[Pi])^2 (4 \[Pi])^Eps Gamma[1+Eps] (musq)^(-Eps))
c  c So, I*I=(-1)!!!
       result(1) = (F(3)*SMB(8)+8*(F(6)*SMB(20)+F(7)*SMB(21))+4*(s12*F
     &   (1)*SMB(5)-F(2)*SMB(6)-F(8)*SMB(24)+F(9)*SMB(25)-F(10)*SMB(2
     &   9))+F(13)*SMB(54)+F(14)*SMB(55)-2*(F(4)*SMB(9)+s12*F(5)*SMB(
     &   10)-F(11)*SMB(31)-F(12)*SMB(32)+F(15)*SMB(58)-F(16)*SMB(59)-
     &   F(17)*SMB(60)+F(18)*SMB(73))-F(201)*(SMB(81)+SMB(82)+SMB(83)
     &   )+F(19)*SMB(84)+F(20)*SMB(85)+F(21)*SMB(86)+F(22)*SMB(87)+F(
     &   23)*SMB(88)+F(24)*SMB(89))/s34
       result(2) = (F(26)*SMB(2)+F(27)*SMB(3)+F(28)*SMB(4)+F(29)*SMB(5
     &   )+Fa(1)*SMB(7)+F(34)*SMB(8)+F(35)*SMB(9)+F(36)*SMB(10)-8*(Fa
     &   (2)*SMB(12)+Fa(3)*SMB(13)+Fa(4)*SMB(15)+Fa(5)*SMB(16)-Fa(6)*
     &   SMB(18)+Fa(7)*SMB(19))+F(52)*SMB(30)+F(53)*SMB(31)+Fa(9)*SMB
     &   (40)+F(61)*SMB(55)+F(62)*SMB(56)-4*(F(25)*SMB(1)-F(50)*SMB(2
     &   2)-F(38)*SMB(23)-F(40)*SMB(25)+F(51)*SMB(26)-F(43)*SMB(27)-F
     &   (45)*SMB(28)+Fa(8)*SMB(39)-Fa(10)*SMB(41)+F(58)*SMB(49)+F(59
     &   )*SMB(50)-Fa(11)*SMB(66)+F(48)*SMB(69))-F(202)*(6*(SMB(42)+S
     &   MB(43)+SMB(46))-2*Fa(12)*SMB(74)+SMB(81)+SMB(82)+SMB(83))+F(
     &   39)*(4*(SMB(24)+SMB(45))+SMB(87))+F(60)*(4*SMB(48)+SMB(88))+
     &   F(67)*SMB(89)+2*(F(41)*(2*SMB(67)+SMB(70))+F(65)*(2*SMB(68)+
     &   SMB(71))+F(66)*SMB(72)-F(54)*SMB(73)-Fa(13)*SMB(76)+Fa(14)*S
     &   MB(78)+Fa(15)*SMB(80)-Fa(16)*SMB(90)+Fa(17)*SMB(91)+Fa(18)*S
     &   MB(92))+F(55)*(4*SMB(44)+SMB(93))+F(57)*(4*SMB(47)+SMB(94))+
     &   F(68)*SMB(95))/s23
       result(3) = Fa(19)*SMB(7)+F(80)*SMB(8)+8*(Fa(20)*SMB(11)-Fa(21)
     &   *SMB(12)-Fa(22)*SMB(13)+Fa(23)*SMB(14)-Fa(24)*SMB(17))+F(138
     &   )*SMB(54)+F(139)*SMB(55)+F(140)*SMB(56)+4*(F(69)*SMB(1)+s12*
     &   (F(70)*SMB(2)+F(71)*SMB(3)+F(72)*SMB(4)+F(73)*SMB(5))-F(74)*
     &   SMB(6)-F(83)*SMB(21)+F(105)*SMB(22)+F(106)*SMB(23)+F(107)*(S
     &   MB(25)+SMB(26))-F(108)*SMB(29)-F(109)*SMB(30)-F(110)*SMB(32)
     &   +F(84)*SMB(33)+F(104)*(-SMB(20)+SMB(34))+F(111)*SMB(35)+F(11
     &   3)*SMB(37)-Fa(25)*SMB(39)+Fa(26)*SMB(40)-Fa(27)*SMB(41)-Fa(2
     &   8)*SMB(57)-Fa(29)*SMB(61)+F(152)*SMB(65)+Fa(30)*SMB(66))-F(2
     &   03)*(8*SMB(31)+SMB(81)+SMB(82)+SMB(83))+F(188)*SMB(84)+F(189
     &   )*SMB(85)+F(190)*SMB(86)+F(129)*(2*SMB(45)+SMB(87))+F(132)*(
     &   2*SMB(48)+SMB(88))+F(134)*(2*SMB(50)+SMB(89))+2*(F(81)*SMB(9
     &   )+F(82)*SMB(10)-F(126)*SMB(42)-F(127)*SMB(43)-F(130)*SMB(46)
     &   -F(135)*SMB(51)-F(136)*SMB(52)-F(137)*SMB(53)+F(143)*SMB(58)
     &   +F(112)*(2*SMB(36)+SMB(59))+F(114)*(2*SMB(38)+SMB(60))+F(149
     &   )*SMB(62)+F(150)*SMB(63)+F(151)*SMB(64)+F(157)*SMB(70)-F(158
     &   )*SMB(71)-F(159)*SMB(72)-F(160)*SMB(73)+Fa(31)*SMB(74)+Fa(32
     &   )*SMB(75)-Fa(33)*SMB(76)+Fa(34)*SMB(77)-Fa(35)*SMB(78)-Fa(36
     &   )*SMB(79)-Fa(37)*SMB(80)-Fa(38)*SMB(90)-Fa(39)*SMB(91)-Fa(40
     &   )*SMB(92))+F(128)*(2*SMB(44)+SMB(93))+F(131)*(2*SMB(47)+SMB(
     &   94))+F(133)*(2*SMB(49)+SMB(95))
       ngluon=1
       If (ngluon.eq.0) then
       resultn=-(result(1)+result(2)+result(3))
       elseif (ngluon.eq.1) then
       result(1)=-(result(1))
       result(2)=-(result(2))
       result(3)=-(result(3))
       else
       Write(*,*) "Error: The position of the gluon is badly indicated
     &   . Look to the heading for explanation" 
       endif
       resultborn = -((2*SMB(6)+SMB(8)+SMB(9)+SMB(10))/((M1SQ-s16)*(M1
     &   SQ-s234)*s34))-(2*SMB(1)+2*SMB(2)+2*SMB(3)-2*SMB(4)-2*SMB(5)
     &   -2*(p1mup3+p5mup3+p6mup3)*SMB(7)+SMB(8)+SMB(9)+SMB(10))/((M1
     &   SQ-s16)*s23*(M1SQ-s234))
c************************************************************************************
c************************************************************************************
c************************************************************************************
       Return
       End
