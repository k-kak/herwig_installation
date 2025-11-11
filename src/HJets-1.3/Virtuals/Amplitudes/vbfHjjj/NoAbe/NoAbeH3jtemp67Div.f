c ************************************************************************************
c ************************************************************************************
c Author: Francisco Campanario
c E-mail: francam@particle.uni-karlsruhe.de
c Date: 18/01/2010
c Modified:6/2/2012
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
c                   |              | 
c                   |              |             
c                   |              |              
c                   |              |             
c barpsi(p1)-->-$$$gluon$$$$$$------>---   bar_psi(p3)
c                     $             
c                     g              
c                     l             
c                     u              
c                     o             
c                     n              
c                     $             
c                     p_2              
c                                               
c Note: To make it shorter in the promgram: mu_p2,...->mup2,... 
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
       subroutine NoAbeH3j67Div_c(M1sq,p1,p2,p3,p4,p5,p6,barpsi_p1,psi_p6,barp
     &   si_p3,psi_p4,mup2,musq,comp,result,resultborn,Div)
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
      Complex*16 SMB(72),Fa(22),F(172)
       Complex*16 SMB0(16),SMB11(0:3),SMB12(0:3),SMB13(0:3),SMB14(0:3)
      Complex*16 SMB1(33)
      Real*8 FI(172),FR(172)
      Complex*16 barpsi_p1_P(2),psi_p6_P(2),barpsi_p3_P(2),psi_p4_P(2) 
       Complex*16 barpsi_p1_M(2),psi_p6_M(2),barpsi_p3_M(2),psi_p4_M(2)
       Complex*16 barpsi_p1(4),psi_p6(4),barpsi_p3(4),psi_p4(4)
       Real*8 P(116)
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
       Complex*16   mup2(0:3)
       Complex*16   p1mup2, p2mup2, p3mup2, p4mup2, p5mup2, p6mup2
       Complex*16   mup3mup4
       Real*8 dotrr
       Complex*16 A0finG_c,B0finG_c,C0finG_c,D0finG_c,E0finG_c,F0finG_c
       EXTERNAL dotrr,A0finG_c,B0finG_c,C0finG_c,D0finG_c,E0finG_c,F0finG_c
       Complex*16 B0finG,C0finG,D0finG,E0finG,F0finG
       EXTERNAL B0finG,C0finG,D0finG,E0finG,F0finG
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
      Complex*16 dotrc,dotcc,result(2),resultn,resultborn
       Real*8 musq
      EXTERNAL   dotrc,dotcc
      Integer alpha
       COMMON/NoAbeH3jFaFunctions/Fa
       COMMON/NoAbeH3jFhlFunctions/F
      Save/NoAbeH3jFhlFunctions/
       COMMON/NoAbeH3jPFunctions/P
      Save/NoAbeH3jPFunctions/
       COMMON/NoAbeH3jInvariants/p1sq,p2sq,p3sq,p4sq,p5sq,p6sq,s12,s23
     &   ,s34,s45,s56,s16,s123,s234,s345
       COMMON/NoAbeH3jFVALFunctions/F0123456R,F123456R,F0123456I,F1234
     &   56I
       COMMON/NoAbeH3jEVALFunctions/ E012345R,E12345R,E012345I,E12345I
     &   , E012346R,E12346R,E012346I,E12346I, E012356R,E12356R,E01235
     &   6I,E12356I, E012456R,E12456R,E012456I,E12456I, E013456R,E134
     &   56R,E013456I,E13456I, E023456R,E23456R,E023456I,E23456I
       COMMON/NoAbeH3jDVALFunctions/ D01234R,D1234R,D01234I,D1234I, D0
     &   1235R,D1235R,D01235I,D1235I, D01236R,D1236R,D01236I,D1236I, 
     &   D01245R,D1245R,D01245I,D1245I, D01246R,D1246R,D01246I,D1246I
     &   , D01256R,D1256R,D01256I,D1256I, D01345R,D1345R,D01345I,D134
     &   5I, D01346R,D1346R,D01346I,D1346I, D01356R,D1356R,D01356I,D1
     &   356I, D01456R,D1456R,D01456I,D1456I, D02345R,D2345R,D02345I,
     &   D2345I, D02346R,D2346R,D02346I,D2346I, D02356R,D2356R,D02356
     &   I,D2356I, D02456R,D2456R,D02456I,D2456I, D03456R,D3456R,D034
     &   56I,D3456I
       COMMON/NoAbeH3jCVALFunctions/ C0123R,C123R,C0123I,C123I, C0124R
     &   ,C124R,C0124I,C124I, C0125R,C125R,C0125I,C125I, C0126R,C126R
     &   ,C0126I,C126I, C0134R,C134R,C0134I,C134I, C0135R,C135R,C0135
     &   I,C135I, C0136R,C136R,C0136I,C136I, C0145R,C145R,C0145I,C145
     &   I, C0146R,C146R,C0146I,C146I, C0156R,C156R,C0156I,C156I, C02
     &   34R,C234R,C0234I,C234I, C0235R,C235R,C0235I,C235I, C0236R,C2
     &   36R,C0236I,C236I, C0245R,C245R,C0245I,C245I, C0246R,C246R,C0
     &   246I,C246I, C0256R,C256R,C0256I,C256I, C0345R,C345R,C0345I,C
     &   345I, C0346R,C346R,C0346I,C346I, C0356R,C356R,C0356I,C356I, 
     &   C0456R,C456R,C0456I,C456I
       COMMON/NoAbeH3jBVALFunctions/ B012R,B012I, B013R,B013I, B014R,B
     &   014I, B015R,B015I, B016R,B016I, B023R,B023I, B024R,B024I, B0
     &   25R,B025I, B026R,B026I, B034R,B034I, B035R,B035I, B036R,B036
     &   I, B045R,B045I, B046R,B046I, B056R,B056I
       Integer ngluon, posgluon,Div
       COMPLEX*16 m1sq,m0
       real*8 m
c Old routines
       logical RealRu
       Parameter (RealRu=.false.)

c       real*8 m
        m0=DCMPLX(0d0,0d0)
c       m1=m

c************************************************************************************
c************************************************************************************
c************************************************************************************
c************************************************************************************
c************************************************************************************
c************************************************************************************
c************************************************************************************
c************************************************************************************
c       Definition of the scalar products. Not inlcueded the contraction of the
c       moments with the external currents  
c************************************************************************************
c************************************************************************************
c************************************************************************************
c************************************************************************************
c************************************************************************************
c************************************************************************************
c************************************************************************************
       p1sq = dotrr(p1,p1)
       p1p2 = dotrr(p1,p2)
       p1p3 = dotrr(p1,p3)
       p1p4 = dotrr(p1,p4)
       p1p5 = dotrr(p1,p5)
       p1p6 = dotrr(p1,p6)
       p2sq = dotrr(p2,p2)
       p2p3 = dotrr(p2,p3)
       p2p4 = dotrr(p2,p4)
       p2p5 = dotrr(p2,p5)
       p2p6 = dotrr(p2,p6)
       p3sq = dotrr(p3,p3)
       p3p4 = dotrr(p3,p4)
       p3p5 = dotrr(p3,p5)
       p3p6 = dotrr(p3,p6)
       p4sq = dotrr(p4,p4)
       p4p5 = dotrr(p4,p5)
       p4p6 = dotrr(p4,p6)
       p5sq = dotrr(p5,p5)
       p5p6 = dotrr(p5,p6)
       p6sq = dotrr(p6,p6)
       s12 = (p1sq +p2sq+ 2*p1p2) 
       s13 = (p1sq +p3sq+ 2*p1p3) 
       s14 = (p1sq +p4sq+ 2*p1p4) 
       s15 = (p1sq +p5sq+ 2*p1p5) 
       s16 = (p1sq +p6sq+ 2*p1p6) 
       s23 = (p2sq +p3sq+ 2*p2p3) 
       s24 = (p2sq +p4sq+ 2*p2p4) 
       s25 = (p2sq +p5sq+ 2*p2p5) 
       s26 = (p2sq +p6sq+ 2*p2p6) 
       s34 = (p3sq +p4sq+ 2*p3p4) 
       s35 = (p3sq +p5sq+ 2*p3p5) 
       s36 = (p3sq +p6sq+ 2*p3p6) 
       s45 = (p4sq +p5sq+ 2*p4p5) 
       s46 = (p4sq +p6sq+ 2*p4p6) 
       s56 = (p5sq +p6sq+ 2*p5p6) 
c       Write(*,'(a5,F20.10)')," p1sq ", p1sq 
c       Write(*,'(a5,F20.10)')," p1p2 ", p1p2
c       Write(*,'(a5,F20.10)')," p1p3 ", p1p3
c       Write(*,'(a5,F20.10)')," p1p4 ", p1p4
c       Write(*,'(a5,F20.10)')," p1p5 ", p1p5
c       Write(*,'(a5,F20.10)')," p1p6 ", p1p6
c       Write(*,'(a5,F20.10)')," p2sq ", p2sq 
c       Write(*,'(a5,F20.10)')," p2p3 ", p2p3
c       Write(*,'(a5,F20.10)')," p2p4 ", p2p4
c       Write(*,'(a5,F20.10)')," p2p5 ", p2p5
c       Write(*,'(a5,F20.10)')," p2p6 ", p2p6
c       Write(*,'(a5,F20.10)')," p3sq ", p3sq 
c       Write(*,'(a5,F20.10)')," p3p4 ", p3p4
c       Write(*,'(a5,F20.10)')," p3p5 ", p3p5
c       Write(*,'(a5,F20.10)')," p3p6 ", p3p6
c       Write(*,'(a5,F20.10)')," p4sq ", p4sq 
c       Write(*,'(a5,F20.10)')," p4p5 ", p4p5
c       Write(*,'(a5,F20.10)')," p4p6 ", p4p6
c       Write(*,'(a5,F20.10)')," p5sq ", p5sq 
c       Write(*,'(a5,F20.10)')," p5p6 ", p5p6
c       Write(*,'(a5,F20.10)')," p6sq ", p6sq 
      s123=p1sq+p2sq+p3sq+2*(p1p2+p1p3+p2p3) 
      s234=p2sq+p3sq+p4sq+2*(p2p3+p2p4+p3p4) 
      s345=p3sq+p4sq+p5sq+2*(p3p4+p3p5+p4p5)
      PrintB=.False. 
      PrintC=.False. 
      PrintD=.False. 
      PrintE=.False. 
      PrintF=.False.
c************************************************************************************
c************************************************************************************
c************************************************************************************
c************************************************************************************
        if (comp.gt.0) then    
c    Calling C_ij,D_ij,E_ij,F_ij Functions    

      call FgetNoAbe_Div_c(M1SQ,p1,p2,p3,p4,p5,p6,musq,Div)

c************************************************************************************
c************************************************************************************
c************************************************************************************
c       Definition of the F,P functions:Independent of the currents    
c************************************************************************************
c************************************************************************************
c************************************************************************************
       call NoAbeH3jFFhl1(F(1))
       call NoAbeH3jFFhl2(F(35))
       call NoAbeH3jFFhl3(F(69))
       call NoAbeH3jFFhl4(F(103))
       call NoAbeH3jFFhl5(F(137))
c************************************************************************************
c************************************************************************************
c************************************************************************************
       endif  
c               PART THAT DEPENDS ON THE EXTERNAL CURRENT
c************************************************************************************
c************************************************************************************
c************************************************************************************
       p1mup2 = dotrc(p1,mup2)
       p2mup2 = dotrc(p2,mup2)
       p3mup2 = dotrc(p3,mup2)
       p4mup2 = dotrc(p4,mup2)
       p5mup2 = dotrc(p5,mup2)
       p6mup2 = dotrc(p6,mup2)
c************** Calling the Fa functions**********************************************************************
c************************************************************************************
c************************************************************************************
       call NoAbeH3jFa1(p1mup2,p2mup2,p3mup2,p4mup2,p5mup2,p6mup2,Fa(1
     &   ))
       call NoAbeH3jFa2(p1mup2,p2mup2,p3mup2,p4mup2,p5mup2,p6mup2,Fa(1
     &   2))
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
       psi_p6_M(1)=psi_p6(1)
       psi_p6_M(2)=psi_p6(2)
       psi_p6_P(1)=psi_p6(3)
       psi_p6_P(2)=psi_p6(4)
       psi_p4_M(1)=psi_p4(1)
       psi_p4_M(2)=psi_p4(2)
       psi_p4_P(1)=psi_p4(3)
       psi_p4_P(2)=psi_p4(4)
       barpsi_p1_P(1)=barpsi_p1(1)
       barpsi_p1_P(2)=barpsi_p1(2)
       barpsi_p1_M(1)=barpsi_p1(3)
       barpsi_p1_M(2)=barpsi_p1(4)
       barpsi_p3_P(1)=barpsi_p3(1)
       barpsi_p3_P(2)=barpsi_p3(2)
       barpsi_p3_M(1)=barpsi_p3(3)
       barpsi_p3_M(2)=barpsi_p3(4)
c************************************************************************************
c************************************************************************************
       SMB0(1) = SC1c(barpsi_p1_P,mup2,psi_p6_P,+1)
       SMB0(2) = SC1r(barpsi_p3_M,p1,psi_p4_M,-1)
       SMB0(3) = SC1r(barpsi_p3_M,p5,psi_p4_M,-1)
       SMB0(4) = SC1r(barpsi_p3_M,p6,psi_p4_M,-1)
       SMB0(5) = SC3rrr(barpsi_p3_M,p6,p5,p1,psi_p4_M,-1)
       SMB0(6) = SC1r(barpsi_p1_P,p5,psi_p6_P,+1)
       SMB0(7) = SC1c(barpsi_p3_M,mup2,psi_p4_M,-1)
       SMB0(8) = SC1r(barpsi_p1_P,p3,psi_p6_P,+1)
       SMB0(9) = SC3crr(barpsi_p3_M,mup2,p6,p1,psi_p4_M,-1)
       SMB0(10) = SC3crr(barpsi_p3_M,mup2,p6,p5,psi_p4_M,-1)
       SMB0(11) = SC1r(barpsi_p1_P,p4,psi_p6_P,+1)
       SMB0(12) = SC3rrr(barpsi_p1_P,p5,p4,p3,psi_p6_P,+1)
       SMB0(13) = SC3crr(barpsi_p1_P,mup2,p4,p3,psi_p6_P,+1)
       SMB0(14) = SC3crr(barpsi_p1_P,mup2,p5,p4,psi_p6_P,+1)
       SMB0(15) = SC3crr(barpsi_p1_P,mup2,p5,p3,psi_p6_P,+1)
       SMB0(16) = SC3crr(barpsi_p3_M,mup2,p5,p1,psi_p4_M,-1)
c************************************************************************************
c************************************************************************************
      do i=0,3
        v1(0)=delta(i,0)
        v1(1)=delta(i,1)
        v1(2)=delta(i,2)
        v1(3)=delta(i,3)
       SMB11(i) = SC1c(barpsi_p1_P,v1,psi_p6_P,+1)
       SMB12(i) = SC1c(barpsi_p3_M,v1,psi_p4_M,-1)
       SMB13(i) = SC3crc(barpsi_p1_P,mup2,p5,v1,psi_p6_P,+1)
       SMB14(i) = SC3crc(barpsi_p3_M,mup2,p5,v1,psi_p4_M,-1)
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
       SMB1(6) = SMB0(6)*SMB0(7)
       SMB1(8) = SMB0(2)*SMB0(6)
       SMB1(9) = SMB0(7)*SMB0(8)
       SMB1(10) = SMB0(4)*SMB0(8)
       SMB1(11) = SMB0(8)*SMB0(9)
       SMB1(12) = SMB0(8)*SMB0(10)
       SMB1(13) = SMB0(7)*SMB0(11)
       SMB1(14) = SMB0(2)*SMB0(11)
       SMB1(15) = SMB0(3)*SMB0(11)
       SMB1(16) = SMB0(4)*SMB0(11)
       SMB1(17) = SMB0(9)*SMB0(11)
       SMB1(18) = SMB0(10)*SMB0(11)
       SMB1(19) = SMB0(4)*SMB0(6)
       SMB1(20) = SMB0(6)*SMB0(9)
       SMB1(21) = SMB0(6)*SMB0(10)
       SMB1(22) = SMB0(7)*SMB0(12)
       SMB1(23) = SMB0(3)*SMB0(8)
       SMB1(24) = SMB0(2)*SMB0(13)
       SMB1(25) = SMB0(2)*SMB0(14)
       SMB1(26) = SMB0(3)*SMB0(13)
       SMB1(27) = SMB0(3)*SMB0(14)
       SMB1(28) = SMB0(4)*SMB0(13)
       SMB1(29) = SMB0(4)*SMB0(14)
       SMB1(30) = SMB0(2)*SMB0(15)
       SMB1(32) = SMB0(3)*SMB0(6)
       SMB1(33) = SMB0(8)*SMB0(16)
c************************************************************************************
       SMB1(5) = dotcc(SMB11,SMB12)
       SMB1(7) = dotcc(SMB12,SMB13)
       SMB1(31) = dotcc(SMB11,SMB14)
c************************************************************************************
c************************************************************************************
       SMB(1) = SMB1(5)
       SMB(2) = SMB1(3)
       SMB(3) = 2*p3mup2*SMB1(5)-2*SMB1(9)
       SMB(4) = 0
       SMB(5) = SMB1(7)
       SMB(6) = SMB1(10)
       SMB(7) = s34*SMB1(5)
       SMB(8) = -((s34-s345+s45)*SMB1(5))-2*SMB1(23)
       SMB(9) = -2*SMB1(1)+2*p1mup2*SMB1(5)
       SMB(10) = -2*SMB1(2)+2*SMB1(6)+SMB1(7)
       SMB(11) = 0
       SMB(12) = SMB1(13)
       SMB(13) = SMB1(9)
       SMB(14) = SMB1(11)
       SMB(15) = SMB1(12)
       SMB(16) = SMB1(16)
       SMB(17) = SMB1(17)
       SMB(18) = SMB1(18)
       SMB(19) = SMB1(6)
       SMB(20) = SMB1(19)
       SMB(21) = SMB1(20)
       SMB(22) = SMB1(21)
       SMB(23) = SMB1(22)
       SMB(24) = 0
       SMB(25) = 0
       SMB(26) = 0
       SMB(27) = -2*s34*(SMB1(1)-p1mup2*SMB1(5))
       SMB(28) = s34*(-2*SMB1(2)+2*SMB1(6)+SMB1(7))
       SMB(29) = 0
       SMB(30) = 2*(2*SMB1(6)+SMB1(7)-SMB1(31))
       SMB(31) = 2*(s34-s345+s45)*SMB1(1)+2*(p2sq-s12+s123-s23)*SMB1(2
     &   )-2*p1mup2*(s34-s345+s45)*SMB1(5)+(-p2sq+s12-s123+s23)*SMB1(
     &   7)+2*(s16-s234+s56)*SMB1(9)+2*SMB1(33)
       SMB(32) = 2*p3mup2*p5sq*SMB1(5)-2*(s34-s345+s45)*SMB1(6)-2*p5sq
     &   *SMB1(9)-4*p3mup2*SMB1(32)
       SMB(33) = 2*(s12-s123-s345+s45)*SMB1(2)+2*p3mup2*(-p5sq+s56)*SM
     &   B1(5)+(-s12+s123+s345-s45)*SMB1(7)-2*SMB1(12)-4*p3mup2*SMB1(
     &   19)
       SMB(34) = 0
       SMB(35) = 0
       SMB(36) = 0
       SMB(37) = s16*SMB1(5)
       SMB(38) = (-p5sq+s56)*SMB1(5)
       SMB(39) = 4*SMB1(5)
       SMB(40) = SMB1(14)
       SMB(41) = 0
       SMB(42) = SMB1(15)
       SMB(43) = 0
       SMB(44) = 0
       SMB(45) = -((s16-s234+s56)*SMB1(5))-2*SMB1(8)
       SMB(46) = 0
       SMB(47) = 4*SMB1(32)
       SMB(48) = 0
       SMB(49) = 0
       SMB(50) = 0
       SMB(51) = SMB1(1)
       SMB(52) = SMB1(2)
       SMB(53) = SMB1(4)
       SMB(54) = SMB1(24)
       SMB(55) = SMB1(25)
       SMB(56) = SMB1(26)
       SMB(57) = SMB1(27)
       SMB(58) = SMB1(28)
       SMB(59) = SMB1(29)
       SMB(60) = 0
       SMB(61) = 2*(2*SMB1(2)-SMB1(7)+SMB1(31))
       SMB(62) = 0
       SMB(63) = 2*(s34-s345+s45)*SMB1(1)-2*p3mup2*(s16-s234+s56)*SMB1
     &   (5)+2*(p2sq-s12+s123-s23)*SMB1(6)+2*(s16-s234+s56)*SMB1(9)+2
     &   *SMB1(30)+(-p2sq+s12-s123+s23)*SMB1(31)
       SMB(64) = 2*s16*(p3mup2*SMB1(5)-SMB1(9))
       SMB(65) = 2*(p5sq-s56)*(-(p3mup2*SMB1(5))+SMB1(9))
       SMB(66) = 2*(p5sq-s45)*SMB1(1)+2*p1mup2*(-p5sq+s45)*SMB1(5)-2*(
     &   s123-s23+s234-s56)*SMB1(6)-4*p1mup2*SMB1(15)+2*SMB1(25)+(s12
     &   3-s23+s234-s56)*SMB1(31)
       SMB(67) = 0
       SMB(68) = 0
       SMB(69) = -2*p5sq*SMB1(1)-2*(s16-s234+s56)*SMB1(2)+2*p1mup2*p5s
     &   q*SMB1(5)-4*p1mup2*SMB1(32)
       SMB(70) = s16*(2*SMB1(2)-2*SMB1(6)+SMB1(31))
       SMB(71) = -((p5sq-s56)*(2*SMB1(2)-2*SMB1(6)+SMB1(31)))
       SMB(72) = 0
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
       result(1) = Fa(1)*SMB(1)+F(7)*SMB(3)+F(8)*SMB(4)+F(9)*SMB(5)+F(
     &   25)*SMB(9)+F(26)*SMB(10)+F(27)*SMB(11)+8*Fa(5)*SMB(16)+Fa(7)
     &   *SMB(24)+Fa(8)*SMB(37)+Fa(10)*SMB(39)+Fa(12)*SMB(41)+Fa(14)*
     &   SMB(43)+Fa(15)*SMB(44)+Fa(16)*SMB(45)+Fa(17)*SMB(46)+Fa(18)*
     &   SMB(47)+Fa(19)*SMB(48)+Fa(20)*SMB(49)+Fa(21)*SMB(50)+4*(Fa(2
     &   )*SMB(6)+F(29)*SMB(13)-F(30)*SMB(14)-F(31)*SMB(15)-F(33)*SMB
     &   (17)-F(34)*SMB(18)-F(35)*SMB(19)-Fa(6)*SMB(20)-F(41)*SMB(21)
     &   -F(42)*SMB(22)+F(43)*SMB(23)+Fa(11)*SMB(40)-Fa(13)*SMB(42)+F
     &   (127)*SMB(51)+F(128)*SMB(52)+F(129)*SMB(53)+F(130)*SMB(54)-F
     &   (131)*SMB(55)-F(132)*SMB(56)+F(133)*SMB(57)-F(32)*SMB(58)+F(
     &   134)*SMB(59))+2*(F(6)*SMB(2)+Fa(3)*SMB(7)+Fa(4)*SMB(8)+F(28)
     &   *SMB(12)+F(49)*SMB(25)+F(50)*SMB(26)-F(51)*SMB(27)-F(52)*SMB
     &   (28)-F(53)*SMB(29)+F(54)*SMB(30)-F(55)*SMB(31)-F(56)*SMB(32)
     &   -F(57)*SMB(33)-F(58)*SMB(34)-F(59)*SMB(35)-F(60)*SMB(36)+Fa(
     &   9)*SMB(38)-F(135)*SMB(60)-F(136)*SMB(61)-F(137)*SMB(62)+F(13
     &   8)*SMB(63)+F(139)*SMB(64)+F(140)*SMB(65)+F(141)*SMB(66)+F(14
     &   2)*SMB(67)+F(143)*SMB(68)+F(144)*SMB(69)+F(145)*SMB(70)+F(14
     &   6)*SMB(71))+Fa(22)*SMB(72)
       result(2) = F(152)*SMB(1)+F(155)*SMB(8)+F(158)*SMB(24)-2*(F(153
     &   )*SMB(6)-F(154)*SMB(7)+F(156)*SMB(16)+F(157)*SMB(20)-F(159)*
     &   SMB(37))+F(160)*SMB(38)+F(161)*SMB(39)+F(162)*SMB(40)+F(163)
     &   *SMB(41)+F(164)*SMB(42)+F(165)*SMB(43)+F(166)*SMB(44)+F(167)
     &   *SMB(45)+F(168)*SMB(46)+F(169)*SMB(47)+F(170)*SMB(48)+F(171)
     &   *SMB(49)+F(172)*SMB(50)
       ngluon=1
       If (ngluon.eq.0) then
       resultn=(0,-1d0)*(result(1))
       elseif (ngluon.eq.1) then
       result(1)=(0,-1d0)*(result(1))
       result(2)=(0,-1d0)*(result(2))
       else
       Write(*,*) "Error: The position of the gluon is badly indicated
     &   . Look to the heading for explanation" 
       endif
       resultborn = -((2*SMB(2)+SMB(3)+SMB(4)+SMB(5))/(s12*(M1SQ-s34)*(
     &   M1SQ-s345)))-(SMB(9)+SMB(10)+SMB(11)+2*SMB(12))/((M1SQ-s16)*s2
     &   3*(M1SQ-s234))-(-2*(p3mup2+p4mup2+p5mup2)*SMB(1)+SMB(3)+SMB(4
     &   )+SMB(5)+2*(SMB(12)+SMB(13)+SMB(19)-SMB(51)-SMB(52)))/((M1SQ-
     &   s34)*(M1SQ-s345)*(p2sq-s12-s16+s345))-(-2*(p1mup2+p5mup2+p6mu
     &   p2)*SMB(1)+2*SMB(2)+SMB(9)+SMB(10)+SMB(11)-2*SMB(13)-2*SMB(1
     &   9)+2*SMB(51)+2*SMB(52))/((M1SQ-s16)*(M1SQ-s234)*(p2sq-s23+s234
     &   -s34))
c************************************************************************************
c************************************************************************************
c************************************************************************************
       Return
       End
