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
       subroutine H3jCrossF76(M,p1,p2,p3,p4,p5,p6,barpsi_p1,psi_p6,psi
     &   _p2,barpsi_p4,mup3,musq,comp,result,resultborn)
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
      Complex*16 SMB(95),Fa(44),F(214)
       Complex*16 SMB0(16),SMB11(0:3),SMB12(0:3),SMB13(0:3),SMB14(0:3)
      Complex*16 SMB1(38)
      Real*8 FI(214),FR(214)
      Complex*16 barpsi_p1_P(2),psi_p6_P(2),psi_p2_P(2),barpsi_p4_P(2) 
       Complex*16 barpsi_p1_M(2),psi_p6_M(2),psi_p2_M(2),barpsi_p4_M(2)
       Complex*16 barpsi_p1(4),psi_p6(4),psi_p2(4),barpsi_p4(4)
       Real*8 P(175)
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
       Real*8 musq,M,m1sq
      EXTERNAL   dotrc,dotcc
      Integer alpha
       COMMON/H3jCrossFFaFunctions/Fa
       COMMON/H3jCrossFFhlFunctions/F
      Save/H3jCrossFFhlFunctions/
       COMMON/H3jCrossFPFunctions/P
      Save/H3jCrossFPFunctions/
       COMMON/H3jCrossFInvariants/p1sq,p2sq,p3sq,p4sq,p5sq,p6sq,s12,s2
     &   3,s34,s45,s56,s16,s123,s234,s345
       COMMON/H3jCrossFFVALFunctions/F0123456R,F123456R,F0123456I,F123
     &   456I
       COMMON/H3jCrossFEVALFunctions/ E012345R,E12345R,E012345I,E12345
     &   I, E012346R,E12346R,E012346I,E12346I, E012356R,E12356R,E0123
     &   56I,E12356I, E012456R,E12456R,E012456I,E12456I, E013456R,E13
     &   456R,E013456I,E13456I, E023456R,E23456R,E023456I,E23456I
       COMMON/H3jCrossFDVALFunctions/ D01234R,D1234R,D01234I,D1234I, D
     &   01235R,D1235R,D01235I,D1235I, D01236R,D1236R,D01236I,D1236I,
     &    D01245R,D1245R,D01245I,D1245I, D01246R,D1246R,D01246I,D1246
     &   I, D01256R,D1256R,D01256I,D1256I, D01345R,D1345R,D01345I,D13
     &   45I, D01346R,D1346R,D01346I,D1346I, D01356R,D1356R,D01356I,D
     &   1356I, D01456R,D1456R,D01456I,D1456I, D02345R,D2345R,D02345I
     &   ,D2345I, D02346R,D2346R,D02346I,D2346I, D02356R,D2356R,D0235
     &   6I,D2356I, D02456R,D2456R,D02456I,D2456I, D03456R,D3456R,D03
     &   456I,D3456I
       COMMON/H3jCrossFCVALFunctions/ C0123R,C123R,C0123I,C123I, C0124
     &   R,C124R,C0124I,C124I, C0125R,C125R,C0125I,C125I, C0126R,C126
     &   R,C0126I,C126I, C0134R,C134R,C0134I,C134I, C0135R,C135R,C013
     &   5I,C135I, C0136R,C136R,C0136I,C136I, C0145R,C145R,C0145I,C14
     &   5I, C0146R,C146R,C0146I,C146I, C0156R,C156R,C0156I,C156I, C0
     &   234R,C234R,C0234I,C234I, C0235R,C235R,C0235I,C235I, C0236R,C
     &   236R,C0236I,C236I, C0245R,C245R,C0245I,C245I, C0246R,C246R,C
     &   0246I,C246I, C0256R,C256R,C0256I,C256I, C0345R,C345R,C0345I,
     &   C345I, C0346R,C346R,C0346I,C346I, C0356R,C356R,C0356I,C356I,
     &    C0456R,C456R,C0456I,C456I
       COMMON/H3jCrossFBVALFunctions/ B012R,B012I, B013R,B013I, B014R,
     &   B014I, B015R,B015I, B016R,B016I, B023R,B023I, B024R,B024I, B
     &   025R,B025I, B026R,B026I, B034R,B034I, B035R,B035I, B036R,B03
     &   6I, B045R,B045I, B046R,B046I, B056R,B056I
       Integer ngluon, posgluon,Div
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
c************************************************************************************
c************************************************************************************
      A01=A0finG(0d0,musq)
      A01R=Dble(A01)
      A01I=DImag(A01)
      A02=A0finG(0d0,musq)
      A02R=Dble(A02)
      A02I=DImag(A02)
      A03=A0finG(0d0,musq)
      A03R=Dble(A03)
      A03I=DImag(A03)
      A04=A0finG(0d0,musq)
      A04R=Dble(A04)
      A04I=DImag(A04)
      A05=A0finG(M,musq)
      A05R=Dble(A05)
      A05I=DImag(A05)
      A06=A0finG(M,musq)
      A06R=Dble(A06)
      A06I=DImag(A06)

      B012=B0finG(0d0,0d0,p1sq,musq)
      B012R=Dble(B012)
      B012I=DImag(B012)
      call tens_red2_new_Re_Com(0d0,0d0,p1sq,A02R,A02I,A01R,A01I,B012R,B
     &   012I,Bij12R,Bij12I)
      B023=B0finG(0d0,0d0,p2sq,musq)
      B023R=Dble(B023)
      B023I=DImag(B023)
      call tens_red2_new_Re_Com(0d0,0d0,p2sq,A03R,A03I,A02R,A02I,B023R,B
     &   023I,Bij23R,Bij23I)
      B034=B0finG(0d0,0d0,p3sq,musq)
      B034R=Dble(B034)
      B034I=DImag(B034)
      call tens_red2_new_Re_Com(0d0,0d0,p3sq,A04R,A04I,A03R,A03I,B034R,B
     &   034I,Bij34R,Bij34I)
      B045=B0finG(0d0,M,p4sq,musq)
      B045R=Dble(B045)
      B045I=DImag(B045)
      call tens_red2_new_Re_Com(0d0,M,p4sq,A05R,A05I,A04R,A04I,B045R,B
     &   045I,Bij45R,Bij45I)
      B056=B0finG(M,M,p5sq,musq)
      B056R=Dble(B056)
      B056I=DImag(B056)
      call tens_red2_new_Re_Com(M,M,p5sq,A06R,A06I,A05R,A05I,B056R,B
     &   056I,Bij56R,Bij56I)
      B016=B0finG(0d0,M,p6sq,musq)
      B016R=Dble(B016)
      B016I=DImag(B016)
      call tens_red2_new_Re_Com(0d0,M,p6sq,A06R,A06I,A01R,A01I,B016R,B
     &   016I,Bij16R,Bij16I)
      B013=B0finG(0d0,0d0,s12,musq)
      B013R=Dble(B013)
      B013I=DImag(B013)
      call tens_red2_new_Re_Com(0d0,0d0,s12,A03R,A03I,A01R,A01I,B013R,B
     &   013I,Bij13R,Bij13I)
      B014=B0finG(0d0,0d0,s123,musq)
      B014R=Dble(B014)
      B014I=DImag(B014)
      call tens_red2_new_Re_Com(0d0,0d0,s123,A04R,A04I,A01R,A01I,B014R,B
     &   014I,Bij14R,Bij14I)
      B026=B0finG(0d0,M,s16,musq)
      B026R=Dble(B026)
      B026I=DImag(B026)
      call tens_red2_new_Re_Com(0d0,M,s16,A06R,A06I,A02R,A02I,B026R,B
     &   026I,Bij26R,Bij26I)
      B024=B0finG(0d0,0d0,s23,musq)
      B024R=Dble(B024)
      B024I=DImag(B024)
      call tens_red2_new_Re_Com(0d0,0d0,s23,A04R,A04I,A02R,A02I,B024R,B
     &   024I,Bij24R,Bij24I)
      B025=B0finG(0d0,M,s234,musq)
      B025R=Dble(B025)
      B025I=DImag(B025)
      call tens_red2_new_Re_Com(0d0,M,s234,A05R,A05I,A02R,A02I,B025R,B
     &   025I,Bij25R,Bij25I)
      B035=B0finG(0d0,M,s34,musq)
      B035R=Dble(B035)
      B035I=DImag(B035)
      call tens_red2_new_Re_Com(0d0,M,s34,A05R,A05I,A03R,A03I,B035R,B
     &   035I,Bij35R,Bij35I)
      B036=B0finG(0d0,M,s345,musq)
      B036R=Dble(B036)
      B036I=DImag(B036)
      call tens_red2_new_Re_Com(0d0,M,s345,A06R,A06I,A03R,A03I,B036R,B
     &   036I,Bij36R,Bij36I)
      B046=B0finG(0d0,M,s45,musq)
      B046R=Dble(B046)
      B046I=DImag(B046)
      call tens_red2_new_Re_Com(0d0,M,s45,A06R,A06I,A04R,A04I,B046R,B
     &   046I,Bij46R,Bij46I)
      B015=B0finG(0d0,M,s56,musq)
      B015R=Dble(B015)
      B015I=DImag(B015)
      call tens_red2_new_Re_Com(0d0,M,s56,A05R,A05I,A01R,A01I,B015R,B
     &   015I,Bij15R,Bij15I)

      C0123=C0finG(0d0,0d0,0d0,p1sq,p2sq,s12,musq)
      call tens_red3_new_Re_Com_G(0d0,0d0,0d0,p1sq,p2sq,s12,
     & B023R,B013R,B012R,B023I,B013I,B012I,Bij23R,Bij13R,Bij12R,Bij23I,Bij13I,Bij12I,
     & C0123,C0123R,C0123I,C123R,C123I)
      C0126=C0finG(0d0,0d0,M,p1sq,s16,p6sq,musq)
      call tens_red3_new_Re_Com_G(0d0,0d0,M,p1sq,s16,p6sq,
     & B026R,B016R,B012R,B026I,B016I,B012I,Bij26R,Bij16R,Bij12R,Bij26I,Bij16I,Bij12I,
     & C0126,C0126R,C0126I,C126R,C126I)
      C0124=C0finG(0d0,0d0,0d0,p1sq,s23,s123,musq)
      call tens_red3_new_Re_Com_G(0d0,0d0,0d0,p1sq,s23,s123,
     & B024R,B014R,B012R,B024I,B014I,B012I,Bij24R,Bij14R,Bij12R,Bij24I,Bij14I,Bij12I,
     & C0124,C0124R,C0124I,C124R,C124I)
      C0125=C0finG(0d0,0d0,M,p1sq,s234,s56,musq)
      call tens_red3_new_Re_Com_G(0d0,0d0,M,p1sq,s234,s56,
     & B025R,B015R,B012R,B025I,B015I,B012I,Bij25R,Bij15R,Bij12R,Bij25I,Bij15I,Bij12I,
     & C0125,C0125R,C0125I,C125R,C125I)
      C0234=C0finG(0d0,0d0,0d0,p2sq,p3sq,s23,musq)
      call tens_red3_new_Re_Com_G(0d0,0d0,0d0,p2sq,p3sq,s23,
     & B034R,B024R,B023R,B034I,B024I,B023I,Bij34R,Bij24R,Bij23R,Bij34I,Bij24I,Bij23I,
     & C0234,C0234R,C0234I,C234R,C234I)
      C0235=C0finG(0d0,0d0,M,p2sq,s34,s234,musq)
      call tens_red3_new_Re_Com_G(0d0,0d0,M,p2sq,s34,s234,
     & B035R,B025R,B023R,B035I,B025I,B023I,Bij35R,Bij25R,Bij23R,Bij35I,Bij25I,Bij23I,
     & C0235,C0235R,C0235I,C235R,C235I)
      C0236=C0finG(0d0,0d0,M,p2sq,s345,s16,musq)
      call tens_red3_new_Re_Com_G(0d0,0d0,M,p2sq,s345,s16,
     & B036R,B026R,B023R,B036I,B026I,B023I,Bij36R,Bij26R,Bij23R,Bij36I,Bij26I,Bij23I,
     & C0236,C0236R,C0236I,C236R,C236I)
      C0345=C0finG(0d0,0d0,M,p3sq,p4sq,s34,musq)
      call tens_red3_new_Re_Com_G(0d0,0d0,M,p3sq,p4sq,s34,
     & B045R,B035R,B034R,B045I,B035I,B034I,Bij45R,Bij35R,Bij34R,Bij45I,Bij35I,Bij34I,
     & C0345,C0345R,C0345I,C345R,C345I)
      C0346=C0finG(0d0,0d0,M,p3sq,s45,s345,musq)
      call tens_red3_new_Re_Com_G(0d0,0d0,M,p3sq,s45,s345,
     & B046R,B036R,B034R,B046I,B036I,B034I,Bij46R,Bij36R,Bij34R,Bij46I,Bij36I,Bij34I,
     & C0346,C0346R,C0346I,C346R,C346I)
      C0456=C0finG(0d0,M,M,p4sq,p5sq,s45,musq)
      call tens_red3_new_Re_Com_G(0d0,M,M,p4sq,p5sq,s45,
     & B056R,B046R,B045R,B056I,B046I,B045I,Bij56R,Bij46R,Bij45R,Bij56I,Bij46I,Bij45I,
     & C0456,C0456R,C0456I,C456R,C456I)
      C0134=C0finG(0d0,0d0,0d0,s12,p3sq,s123,musq)
      call tens_red3_new_Re_Com_G(0d0,0d0,0d0,s12,p3sq,s123,
     & B034R,B014R,B013R,B034I,B014I,B013I,Bij34R,Bij14R,Bij13R,Bij34I,Bij14I,Bij13I,
     & C0134,C0134R,C0134I,C134R,C134I)
      C0135=C0finG(0d0,0d0,M,s12,s34,s56,musq)
      call tens_red3_new_Re_Com_G(0d0,0d0,M,s12,s34,s56,
     & B035R,B015R,B013R,B035I,B015I,B013I,Bij35R,Bij15R,Bij13R,Bij35I,Bij15I,Bij13I,
     & C0135,C0135R,C0135I,C135R,C135I)
      C0136=C0finG(0d0,0d0,M,s12,s345,p6sq,musq)
      call tens_red3_new_Re_Com_G(0d0,0d0,M,s12,s345,p6sq,
     & B036R,B016R,B013R,B036I,B016I,B013I,Bij36R,Bij16R,Bij13R,Bij36I,Bij16I,Bij13I,
     & C0136,C0136R,C0136I,C136R,C136I)
      C0145=C0finG(0d0,0d0,M,s123,p4sq,s56,musq)
      call tens_red3_new_Re_Com_G(0d0,0d0,M,s123,p4sq,s56,
     & B045R,B015R,B014R,B045I,B015I,B014I,Bij45R,Bij15R,Bij14R,Bij45I,Bij15I,Bij14I,
     & C0145,C0145R,C0145I,C145R,C145I)
      C0146=C0finG(0d0,0d0,M,s123,s45,p6sq,musq)
      call tens_red3_new_Re_Com_G(0d0,0d0,M,s123,s45,p6sq,
     & B046R,B016R,B014R,B046I,B016I,B014I,Bij46R,Bij16R,Bij14R,Bij46I,Bij16I,Bij14I,
     & C0146,C0146R,C0146I,C146R,C146I)
      C0245=C0finG(0d0,0d0,M,s23,p4sq,s234,musq)
      call tens_red3_new_Re_Com_G(0d0,0d0,M,s23,p4sq,s234,
     & B045R,B025R,B024R,B045I,B025I,B024I,Bij45R,Bij25R,Bij24R,Bij45I,Bij25I,Bij24I,
     & C0245,C0245R,C0245I,C245R,C245I)
      C0246=C0finG(0d0,0d0,M,s23,s45,s16,musq)
      call tens_red3_new_Re_Com_G(0d0,0d0,M,s23,s45,s16,
     & B046R,B026R,B024R,B046I,B026I,B024I,Bij46R,Bij26R,Bij24R,Bij46I,Bij26I,Bij24I,
     & C0246,C0246R,C0246I,C246R,C246I)
      C0256=C0finG(0d0,M,M,s234,p5sq,s16,musq)
      call tens_red3_new_Re_Com_G(0d0,M,M,s234,p5sq,s16,
     & B056R,B026R,B025R,B056I,B026I,B025I,Bij56R,Bij26R,Bij25R,Bij56I,Bij26I,Bij25I,
     & C0256,C0256R,C0256I,C256R,C256I)
      C0356=C0finG(0d0,M,M,s34,p5sq,s345,musq)
      call tens_red3_new_Re_Com_G(0d0,M,M,s34,p5sq,s345,
     & B056R,B036R,B035R,B056I,B036I,B035I,Bij56R,Bij36R,Bij35R,Bij56I,Bij36I,Bij35I,
     & C0356,C0356R,C0356I,C356R,C356I)
      C0156=C0finG(0d0,M,M,s56,p5sq,p6sq,musq)
      call tens_red3_new_Re_Com_G(0d0,M,M,s56,p5sq,p6sq,
     & B056R,B016R,B015R,B056I,B016I,B015I,Bij56R,Bij16R,Bij15R,Bij56I,Bij16I,Bij15I,
     & C0156,C0156R,C0156I,C156R,C156I)

      D01234=D0finG(0d0,0d0,0d0,0d0,s12,s23,p1sq,p2sq,p3sq,s123,musq)
      call tens_red4_new_Re_Com_G(0d0,0d0,0d0,0d0,p1sq,p2sq,p3sq,p1p2,p1p3,p2p3,C0234R,
     &   C0134R,C0124R,C0123R,C234R,C134R,C124R,C123R,C0234I,C0134I,C
     &   0124I,C0123I,C234I,C134I,C124I,C123I,D01234,D01234R,D01234I,
     &   D1234R,D1234I)
      D01235=D0finG(0d0,0d0,0d0,M,s12,s234,p1sq,p2sq,s34,s56,musq)
      call tens_red4_new_Re_Com_G(0d0,0d0,0d0,M,p1sq,p2sq,s34,p1p2,p1p3+p1p4,p2p3+p2p
     &   4,C0235R,C0135R,C0125R,C0123R,C235R,C135R,C125R,C123R,C0235I
     &   ,C0135I,C0125I,C0123I,C235I,C135I,C125I,C123I,D01235,D01235R
     &   ,D01235I,D1235R,D1235I)
      D01236=D0finG(0d0,0d0,0d0,M,s12,s16,p1sq,p2sq,s345,p6sq,musq)
      call tens_red4_new_Re_Com_G(0d0,0d0,0d0,M,p1sq,p2sq,s345,p1p2,p1p3+p1p4+p1p5,p2
     &   p3+p2p4+p2p5,C0236R,C0136R,C0126R,C0123R,C236R,C136R,C126R,C
     &   123R,C0236I,C0136I,C0126I,C0123I,C236I,C136I,C126I,C123I,D01
     &   236,D01236R,D01236I,D1236R,D1236I)
      D01245=D0finG(0d0,0d0,0d0,M,s123,s234,p1sq,s23,p4sq,s56,musq)
      call tens_red4_new_Re_Com_G(0d0,0d0,0d0,M,p1sq,s23,p4sq,p1p2+p1p3,p1p4,p2p4+p3p
     &   4,C0245R,C0145R,C0125R,C0124R,C245R,C145R,C125R,C124R,C0245I
     &   ,C0145I,C0125I,C0124I,C245I,C145I,C125I,C124I,D01245,D01245R
     &   ,D01245I,D1245R,D1245I)
      D01246=D0finG(0d0,0d0,0d0,M,s123,s16,p1sq,s23,s45,p6sq,musq)
      call tens_red4_new_Re_Com_G(0d0,0d0,0d0,M,p1sq,s23,s45,p1p2+p1p3,p1p4+p1p5,p2p4
     &   +p2p5+p3p4+p3p5,C0246R,C0146R,C0126R,C0124R,C246R,C146R,C126
     &   R,C124R,C0246I,C0146I,C0126I,C0124I,C246I,C146I,C126I,C124I,
     &   D01246,D01246R,D01246I,D1246R,D1246I)
      D01256=D0finG(0d0,0d0,M,M,s56,s16,p1sq,s234,p5sq,p6sq,musq)
      call tens_red4_new_Re_Com_G(0d0,0d0,M,M,p1sq,s234,p5sq,p1p2+p1p3+p1p4,p1p5,p2
     &   p5+p3p5+p4p5,C0256R,C0156R,C0126R,C0125R,C256R,C156R,C126R,C
     &   125R,C0256I,C0156I,C0126I,C0125I,C256I,C156I,C126I,C125I,D01
     &   256,D01256R,D01256I,D1256R,D1256I)
      D02345=D0finG(0d0,0d0,0d0,M,s23,s34,p2sq,p3sq,p4sq,s234,musq)
      call tens_red4_new_Re_Com_G(0d0,0d0,0d0,M,p2sq,p3sq,p4sq,p2p3,p2p4,p3p4,C0345R,
     &   C0245R,C0235R,C0234R,C345R,C245R,C235R,C234R,C0345I,C0245I,C
     &   0235I,C0234I,C345I,C245I,C235I,C234I,D02345,D02345R,D02345I,
     &   D2345R,D2345I)
      D02346=D0finG(0d0,0d0,0d0,M,s23,s345,p2sq,p3sq,s45,s16,musq)
      call tens_red4_new_Re_Com_G(0d0,0d0,0d0,M,p2sq,p3sq,s45,p2p3,p2p4+p2p5,p3p4+p3p
     &   5,C0346R,C0246R,C0236R,C0234R,C346R,C246R,C236R,C234R,C0346I
     &   ,C0246I,C0236I,C0234I,C346I,C246I,C236I,C234I,D02346,D02346R
     &   ,D02346I,D2346R,D2346I)
      D02356=D0finG(0d0,0d0,M,M,s234,s345,p2sq,s34,p5sq,s16,musq)
      call tens_red4_new_Re_Com_G(0d0,0d0,M,M,p2sq,s34,p5sq,p2p3+p2p4,p2p5,p3p5+p4p
     &   5,C0356R,C0256R,C0236R,C0235R,C356R,C256R,C236R,C235R,C0356I
     &   ,C0256I,C0236I,C0235I,C356I,C256I,C236I,C235I,D02356,D02356R
     &   ,D02356I,D2356R,D2356I)
      D03456=D0finG(0d0,0d0,M,M,s34,s45,p3sq,p4sq,p5sq,s345,musq)
      call tens_red4_new_Re_Com_G(0d0,0d0,M,M,p3sq,p4sq,p5sq,p3p4,p3p5,p4p5,C0456R,
     &   C0356R,C0346R,C0345R,C456R,C356R,C346R,C345R,C0456I,C0356I,C
     &   0346I,C0345I,C456I,C356I,C346I,C345I,D03456,D03456R,D03456I,
     &   D3456R,D3456I)
      D01345=D0finG(0d0,0d0,0d0,M,s123,s34,s12,p3sq,p4sq,s56,musq)
      call tens_red4_new_Re_Com_G(0d0,0d0,0d0,M,s12,p3sq,p4sq,p1p3+p2p3,p1p4+p2p4,p3p
     &   4,C0345R,C0145R,C0135R,C0134R,C345R,C145R,C135R,C134R,C0345I
     &   ,C0145I,C0135I,C0134I,C345I,C145I,C135I,C134I,D01345,D01345R
     &   ,D01345I,D1345R,D1345I)
      D01346=D0finG(0d0,0d0,0d0,M,s123,s345,s12,p3sq,s45,p6sq,musq)
      call tens_red4_new_Re_Com_G(0d0,0d0,0d0,M,s12,p3sq,s45,p1p3+p2p3,p1p4+p1p5+p2p4
     &   +p2p5,p3p4+p3p5,C0346R,C0146R,C0136R,C0134R,C346R,C146R,C136
     &   R,C134R,C0346I,C0146I,C0136I,C0134I,C346I,C146I,C136I,C134I,
     &   D01346,D01346R,D01346I,D1346R,D1346I)
      D01356=D0finG(0d0,0d0,M,M,s56,s345,s12,s34,p5sq,p6sq,musq)
      call tens_red4_new_Re_Com_G(0d0,0d0,M,M,s12,s34,p5sq,p1p3+p1p4+p2p3+p2p4,p1p5
     &   +p2p5,p3p5+p4p5,C0356R,C0156R,C0136R,C0135R,C356R,C156R,C136
     &   R,C135R,C0356I,C0156I,C0136I,C0135I,C356I,C156I,C136I,C135I,
     &   D01356,D01356R,D01356I,D1356R,D1356I)
      D01456=D0finG(0d0,0d0,M,M,s56,s45,s123,p4sq,p5sq,p6sq,musq)
      call tens_red4_new_Re_Com_G(0d0,0d0,M,M,s123,p4sq,p5sq,p1p4+p2p4+p3p4,p1p5+p2
     &   p5+p3p5,p4p5,C0456R,C0156R,C0146R,C0145R,C456R,C156R,C146R,C
     &   145R,C0456I,C0156I,C0146I,C0145I,C456I,C156I,C146I,C145I,D01
     &   456,D01456R,D01456I,D1456R,D1456I)
      D02456=D0finG(0d0,0d0,M,M,s234,s45,s23,p4sq,p5sq,s16,musq)
      call tens_red4_new_Re_Com_G(0d0,0d0,M,M,s23,p4sq,p5sq,p2p4+p3p4,p2p5+p3p5,p4p
     &   5,C0456R,C0256R,C0246R,C0245R,C456R,C256R,C246R,C245R,C0456I
     &   ,C0256I,C0246I,C0245I,C456I,C256I,C246I,C245I,D02456,D02456R
     &   ,D02456I,D2456R,D2456I)

      E012345=E0finG(0d0,0d0,0d0,0d0,M,p1sq,p2sq,p3sq,p4sq,s56,s12,s23,s34,s123,s234,D02
     &   345,D01345,D01245,D01235,D01234)
      E012345R=Dble(E012345)
      E012345I=DImag(E012345)
      call tens_red5_new_Re_Com_G(0d0,0d0,0d0,0d0,M,p1sq,p2sq,p3sq,p4sq,p1p2,p1p3,p1p4,p2
     &   p3,p2p4,p3p4,D02345R,D01345R,D01245R,D01235R,D01234R,D2345R,
     &   D1345R,D1245R,D1235R,D1234R,D02345I,D01345I,D01245I,D01235I,
     &   D01234I,D2345I,D1345I,D1245I,D1235I,D1234I,E12345R,E12345I)
      E012346=E0finG(0d0,0d0,0d0,0d0,M,p1sq,p2sq,p3sq,s45,p6sq,s12,s23,s345,s123,s16,D02
     &   346,D01346,D01246,D01236,D01234)
      E012346R=Dble(E012346)
      E012346I=DImag(E012346)
      call tens_red5_new_Re_Com_G(0d0,0d0,0d0,0d0,M,p1sq,p2sq,p3sq,s45,p1p2,p1p3,p1p4+p1p
     &   5,p2p3,p2p4+p2p5,p3p4+p3p5,D02346R,D01346R,D01246R,D01236R,D
     &   01234R,D2346R,D1346R,D1246R,D1236R,D1234R,D02346I,D01346I,D0
     &   1246I,D01236I,D01234I,D2346I,D1346I,D1246I,D1236I,D1234I,E12
     &   346R,E12346I)
      E012356=E0finG(0d0,0d0,0d0,M,M,p1sq,p2sq,s34,p5sq,p6sq,s12,s234,s345,s56,s16,D02
     &   356,D01356,D01256,D01236,D01235)
      E012356R=Dble(E012356)
      E012356I=DImag(E012356)
      call tens_red5_new_Re_Com_G(0d0,0d0,0d0,M,M,p1sq,p2sq,s34,p5sq,p1p2,p1p3+p1p4,p1p
     &   5,p2p3+p2p4,p2p5,p3p5+p4p5,D02356R,D01356R,D01256R,D01236R,D
     &   01235R,D2356R,D1356R,D1256R,D1236R,D1235R,D02356I,D01356I,D0
     &   1256I,D01236I,D01235I,D2356I,D1356I,D1256I,D1236I,D1235I,E12
     &   356R,E12356I)
      E012456=E0finG(0d0,0d0,0d0,M,M,p1sq,s23,p4sq,p5sq,p6sq,s123,s234,s45,s56,s16,D02
     &   456,D01456,D01256,D01246,D01245)
      E012456R=Dble(E012456)
      E012456I=DImag(E012456)
      call tens_red5_new_Re_Com_G(0d0,0d0,0d0,M,M,p1sq,s23,p4sq,p5sq,p1p2+p1p3,p1p4,p1p
     &   5,p2p4+p3p4,p2p5+p3p5,p4p5,D02456R,D01456R,D01256R,D01246R,D
     &   01245R,D2456R,D1456R,D1256R,D1246R,D1245R,D02456I,D01456I,D0
     &   1256I,D01246I,D01245I,D2456I,D1456I,D1256I,D1246I,D1245I,E12
     &   456R,E12456I)
      E023456=E0finG(0d0,0d0,0d0,M,M,p2sq,p3sq,p4sq,p5sq,s16,s23,s34,s45,s234,s345,D03
     &   456,D02456,D02356,D02346,D02345)
      E023456R=Dble(E023456)
      E023456I=DImag(E023456)
      call tens_red5_new_Re_Com_G(0d0,0d0,0d0,M,M,p2sq,p3sq,p4sq,p5sq,p2p3,p2p4,p2p5,p3
     &   p4,p3p5,p4p5,D03456R,D02456R,D02356R,D02346R,D02345R,D3456R,
     &   D2456R,D2356R,D2346R,D2345R,D03456I,D02456I,D02356I,D02346I,
     &   D02345I,D3456I,D2456I,D2356I,D2346I,D2345I,E23456R,E23456I)
      E013456=E0finG(0d0,0d0,0d0,M,M,s12,p3sq,p4sq,p5sq,p6sq,s123,s34,s45,s56,s345,D03
     &  456,D01456,D01356,D01346,D01345)
      E013456R=Dble(E013456)
      E013456I=DImag(E013456)
      call tens_red5_new_Re_Com_G(0d0,0d0,0d0,M,M,s12,p3sq,p4sq,p5sq,p1p3+p2p3,p1p4+p2p
     &   4,p1p5+p2p5,p3p4,p3p5,p4p5,D03456R,D01456R,D01356R,D01346R,D
     &   01345R,D3456R,D1456R,D1356R,D1346R,D1345R,D03456I,D01456I,D0
     &   1356I,D01346I,D01345I,D3456I,D1456I,D1356I,D1346I,D1345I,E13
     &   456R,E13456I)

      F0123456=F0finG(0d0,0d0,0d0,0d0,M,M,p1sq,p2sq,p3sq,p4sq,p5sq,p6sq,s12,s23,s34,s45,s5
     &   6,s16,s123,s234,s345,E023456,E013456,E012456,E012356,E012346
     &   ,E012345)
      F0123456R=Dble(F0123456)
      F0123456I=DImag(F0123456)
      call tens_red6_new_Re_Com_G(0d0,0d0,0d0,0d0,M,M,p1sq,p2sq,p3sq,p4sq,p5sq,p1p2,p1p3,p1
     &   p4,p1p5,p2p3,p2p4,p2p5,p3p4,p3p5,p4p5,E023456R,E013456R,E012
     &   456R,E012356R,E012346R,E012345R,E23456R,E13456R,E12456R,E123
     &   56R,E12346R,E12345R,E023456I,E013456I,E012456I,E012356I,E012
     &   346I,E012345I,E23456I,E13456I,E12456I,E12356I,E12346I,E12345
     &   I,F123456R,F123456I,F123456)
c************************************************************************************
c************************************************************************************
c************************************************************************************
c       Definition of the F,P functions:Independent of the currents    
c************************************************************************************
c************************************************************************************
c************************************************************************************
       call H3jCrossFFFhl1(F(1))
       call H3jCrossFFFhl2(F(43))
       call H3jCrossFFFhl3(F(85))
       call H3jCrossFFFhl4(F(127))
       call H3jCrossFFFhl5(F(169))
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
       call H3jCrossFFa1(p1mup3,p1mup4,p2mup3,p2mup4,p3mup3,p3mup4,p4m
     &   up3,p4mup4,p5mup3,p5mup4,p6mup3,p6mup4,mup3mup4,Fa(1))
       call H3jCrossFFa2(p1mup3,p1mup4,p2mup3,p2mup4,p3mup3,p3mup4,p4m
     &   up3,p4mup4,p5mup3,p5mup4,p6mup3,p6mup4,mup3mup4,Fa(23))
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
       barpsi_p4_P(1)=barpsi_p4(1)
       barpsi_p4_P(2)=barpsi_p4(2)
       barpsi_p4_M(1)=barpsi_p4(3)
       barpsi_p4_M(2)=barpsi_p4(4)
       barpsi_p1_P(1)=barpsi_p1(1)
       barpsi_p1_P(2)=barpsi_p1(2)
       barpsi_p1_M(1)=barpsi_p1(3)
       barpsi_p1_M(2)=barpsi_p1(4)
       psi_p2_M(1)=psi_p2(1)
       psi_p2_M(2)=psi_p2(2)
       psi_p2_P(1)=psi_p2(3)
       psi_p2_P(2)=psi_p2(4)
c************************************************************************************
c************************************************************************************
       SMB0(1) = SC1c(barpsi_p1_M,mup3,psi_p6_M,-1)
       SMB0(2) = SC1r(barpsi_p4_P,p1,psi_p2_P,+1)
       SMB0(3) = SC1r(barpsi_p4_P,p5,psi_p2_P,+1)
       SMB0(4) = SC1r(barpsi_p4_P,p6,psi_p2_P,+1)
       SMB0(5) = SC3rrr(barpsi_p4_P,p6,p5,p1,psi_p2_P,+1)
       SMB0(6) = SC1r(barpsi_p1_M,p5,psi_p6_M,-1)
       SMB0(7) = SC1c(barpsi_p4_P,mup3,psi_p2_P,+1)
       SMB0(8) = SC3crr(barpsi_p4_P,mup3,p6,p1,psi_p2_P,+1)
       SMB0(9) = SC1r(barpsi_p1_M,p2,psi_p6_M,-1)
       SMB0(10) = SC3crr(barpsi_p4_P,mup3,p5,p1,psi_p2_P,+1)
       SMB0(11) = SC3crr(barpsi_p4_P,mup3,p6,p5,psi_p2_P,+1)
       SMB0(12) = SC1r(barpsi_p1_M,p4,psi_p6_M,-1)
       SMB0(13) = SC3rrr(barpsi_p1_M,p5,p4,p2,psi_p6_M,-1)
       SMB0(14) = SC3crr(barpsi_p1_M,mup3,p4,p2,psi_p6_M,-1)
       SMB0(15) = SC3crr(barpsi_p1_M,mup3,p5,p2,psi_p6_M,-1)
       SMB0(16) = SC3crr(barpsi_p1_M,mup3,p5,p4,psi_p6_M,-1)
c************************************************************************************
c************************************************************************************
      do i=0,3
        v1(0)=delta(i,0)
        v1(1)=delta(i,1)
        v1(2)=delta(i,2)
        v1(3)=delta(i,3)
       SMB11(i) = SC1c(barpsi_p1_M,v1,psi_p6_M,-1)
       SMB12(i) = SC1c(barpsi_p4_P,v1,psi_p2_P,+1)
       SMB13(i) = SC3crc(barpsi_p1_M,mup3,p5,v1,psi_p6_M,-1)
       SMB14(i) = SC3crc(barpsi_p4_P,mup3,p5,v1,psi_p2_P,+1)
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
       SMB1(9) = SMB0(6)*SMB0(8)
       SMB1(10) = SMB0(7)*SMB0(9)
       SMB1(11) = SMB0(2)*SMB0(9)
       SMB1(12) = SMB0(3)*SMB0(9)
       SMB1(13) = SMB0(4)*SMB0(9)
       SMB1(14) = SMB0(9)*SMB0(10)
       SMB1(15) = SMB0(8)*SMB0(9)
       SMB1(16) = SMB0(9)*SMB0(11)
       SMB1(17) = SMB0(7)*SMB0(12)
       SMB1(18) = SMB0(2)*SMB0(12)
       SMB1(19) = SMB0(3)*SMB0(12)
       SMB1(20) = SMB0(4)*SMB0(12)
       SMB1(21) = SMB0(10)*SMB0(12)
       SMB1(22) = SMB0(8)*SMB0(12)
       SMB1(23) = SMB0(11)*SMB0(12)
       SMB1(24) = SMB0(3)*SMB0(6)
       SMB1(25) = SMB0(4)*SMB0(6)
       SMB1(26) = SMB0(6)*SMB0(10)
       SMB1(27) = SMB0(6)*SMB0(11)
       SMB1(28) = SMB0(7)*SMB0(13)
       SMB1(29) = SMB0(2)*SMB0(14)
       SMB1(30) = SMB0(2)*SMB0(15)
       SMB1(31) = SMB0(2)*SMB0(16)
       SMB1(32) = SMB0(3)*SMB0(14)
       SMB1(33) = SMB0(3)*SMB0(15)
       SMB1(34) = SMB0(3)*SMB0(16)
       SMB1(35) = SMB0(4)*SMB0(14)
       SMB1(36) = SMB0(4)*SMB0(15)
       SMB1(37) = SMB0(4)*SMB0(16)
c************************************************************************************
       SMB1(5) = dotcc(SMB11,SMB12)
       SMB1(7) = dotcc(SMB12,SMB13)
       SMB1(38) = dotcc(SMB11,SMB14)
c************************************************************************************
c************************************************************************************
       SMB(1) = SMB1(17)
       SMB(2) = SMB1(6)
       SMB(3) = SMB1(1)
       SMB(4) = SMB1(2)
       SMB(5) = SMB1(3)
       SMB(6) = SMB1(5)
       SMB(7) = 2*p4mup3*SMB1(5)-2*SMB1(17)
       SMB(8) = SMB1(7)
       SMB(9) = SMB1(10)
       SMB(10) = 0
       SMB(11) = SMB1(11)
       SMB(12) = SMB1(12)
       SMB(13) = SMB1(13)
       SMB(14) = SMB1(18)
       SMB(15) = SMB1(19)
       SMB(16) = SMB1(20)
       SMB(17) = SMB1(8)
       SMB(18) = SMB1(24)
       SMB(19) = SMB1(25)
       SMB(20) = SMB1(29)
       SMB(21) = SMB1(30)
       SMB(22) = SMB1(32)
       SMB(23) = SMB1(33)
       SMB(24) = SMB1(35)
       SMB(25) = SMB1(36)
       SMB(26) = SMB1(28)
       SMB(27) = SMB1(31)
       SMB(28) = SMB1(34)
       SMB(29) = SMB1(37)
       SMB(30) = 0
       SMB(31) = 0
       SMB(32) = -2*SMB1(2)+2*SMB1(6)+SMB1(7)
       SMB(33) = 2*(2*SMB1(6)+SMB1(7)-SMB1(38))
       SMB(34) = 0
       SMB(35) = -2*SMB1(1)+2*p1mup3*SMB1(5)
       SMB(36) = SMB1(15)
       SMB(37) = SMB1(16)
       SMB(38) = SMB1(22)
       SMB(39) = SMB1(23)
       SMB(40) = SMB1(9)
       SMB(41) = SMB1(27)
       SMB(42) = 0
       SMB(43) = 0
       SMB(44) = (-p5sq+s45)*SMB1(5)-2*SMB1(19)
       SMB(45) = 0
       SMB(46) = 0
       SMB(47) = 0
       SMB(48) = 0
       SMB(49) = 0
       SMB(50) = 0
       SMB(51) = 2*(p5sq-s45)*SMB1(1)-2*(s123-s23+s234-s56)*SMB1(2)+2*
     &   p1mup3*(-p5sq+s45)*SMB1(5)+(s123-s23+s234-s56)*SMB1(7)+2*(s1
     &   6-s234+s56)*SMB1(17)+2*SMB1(21)
       SMB(52) = 2*p4mup3*p5sq*SMB1(5)+2*(-p5sq+s45)*SMB1(6)-2*p5sq*SM
     &   B1(17)-4*p4mup3*SMB1(24)
       SMB(53) = 2*(p5sq+s123-s45-s56)*SMB1(2)+2*p4mup3*(-p5sq+s56)*SM
     &   B1(5)+(-p5sq-s123+s45+s56)*SMB1(7)-2*SMB1(23)-4*p4mup3*SMB1(
     &   25)
       SMB(54) = SMB1(14)
       SMB(55) = SMB1(21)
       SMB(56) = SMB1(26)
       SMB(57) = SMB1(4)
       SMB(58) = -((s16-s234+s56)*SMB1(5))-2*SMB1(8)
       SMB(59) = s16*SMB1(5)
       SMB(60) = (-p5sq+s56)*SMB1(5)
       SMB(61) = 0
       SMB(62) = 2*(2*SMB1(2)-SMB1(7)+SMB1(38))
       SMB(63) = 0
       SMB(64) = 2*(p5sq-s45)*SMB1(1)-2*p4mup3*(s16-s234+s56)*SMB1(5)-
     &   2*(s123-s23+s234-s56)*SMB1(6)+2*(s16-s234+s56)*SMB1(17)+2*SM
     &   B1(31)+(s123-s23+s234-s56)*SMB1(38)
       SMB(65) = 2*s16*(p4mup3*SMB1(5)-SMB1(17))
       SMB(66) = 2*(p5sq-s56)*(-(p4mup3*SMB1(5))+SMB1(17))
       SMB(67) = -2*p5sq*SMB1(1)-2*(s16-s234+s56)*SMB1(2)+2*p1mup3*p5s
     &   q*SMB1(5)-4*p1mup3*SMB1(24)
       SMB(68) = s16*(2*SMB1(2)-2*SMB1(6)+SMB1(38))
       SMB(69) = -((p5sq-s56)*(2*SMB1(2)-2*SMB1(6)+SMB1(38)))
       SMB(70) = -2*(s16-s234+s34-s345)*SMB1(1)+2*p1mup3*(s16-s234+s34
     &   -s345)*SMB1(5)+2*s12*SMB1(6)-4*p1mup3*SMB1(12)+2*SMB1(30)-s1
     &   2*SMB1(38)
       SMB(71) = 0
       SMB(72) = 0
       SMB(73) = -2*SMB1(9)
       SMB(74) = 4*SMB1(5)
       SMB(75) = -8*SMB1(1)+8*p1mup3*SMB1(5)
       SMB(76) = -8*SMB1(2)+8*SMB1(6)+4*SMB1(7)
       SMB(77) = 0
       SMB(78) = 0
       SMB(79) = 0
       SMB(80) = 0
       SMB(81) = 0
       SMB(82) = 0
       SMB(83) = 0
       SMB(84) = 0
       SMB(85) = 4*SMB1(24)
       SMB(86) = 0
       SMB(87) = 0
       SMB(88) = 0
       SMB(89) = 4*SMB1(27)
       SMB(90) = 0
       SMB(91) = 0
       SMB(92) = 0
       SMB(93) = 0
       SMB(94) = 0
       SMB(95) = 4*(p5sq+s123-s45-s56)*SMB1(6)-8*p4mup3*SMB1(25)
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
       result(1) = (F(5)*SMB(30)+F(6)*SMB(31)+F(8)*SMB(33)-8*F(10)*SMB
     &   (36)+4*(F(1)*SMB(2)+F(3)*SMB(9)+F(11)*SMB(37)-F(12)*SMB(40)+
     &   F(17)*SMB(54)+F(18)*SMB(56))+2*(F(2)*SMB(8)+F(4)*SMB(10)+F(7
     &   )*SMB(32)+F(9)*SMB(35)-F(14)*SMB(48)+F(15)*SMB(49)+F(16)*SMB
     &   (50)+F(19)*SMB(73))-F(212)*(SMB(75)+SMB(76)+SMB(77))+F(20)*S
     &   MB(81)+F(21)*SMB(82)+F(22)*SMB(83)+F(23)*SMB(87)+F(24)*SMB(8
     &   8)+F(13)*(-4*SMB(41)+SMB(89)))/s34
       result(2) = (F(25)*SMB(1)+F(26)*SMB(2)+F(27)*SMB(3)+F(28)*SMB(4
     &   )+F(29)*SMB(5)+Fa(1)*SMB(6)+F(34)*SMB(7)+F(35)*SMB(8)-8*(Fa(
     &   2)*SMB(14)-Fa(3)*SMB(15)+Fa(4)*SMB(16)+Fa(5)*SMB(17)+Fa(6)*S
     &   MB(18)+Fa(7)*SMB(19))+F(49)*SMB(30)+F(50)*SMB(32)+F(51)*SMB(
     &   33)+F(52)*SMB(34)+F(53)*SMB(35)+F(58)*SMB(57)-4*(F(47)*SMB(2
     &   9)-F(37)*SMB(38)-F(54)*SMB(39)-F(43)*SMB(40)+F(45)*SMB(41)-F
     &   a(8)*SMB(44)-F(38)*SMB(55)+Fa(9)*SMB(58)+Fa(10)*SMB(59)-Fa(1
     &   1)*SMB(60)+F(63)*SMB(66)+F(64)*SMB(69))-F(213)*(6*(SMB(61)+S
     &   MB(62)+SMB(63))-2*Fa(12)*SMB(74)+SMB(75)+SMB(76)+SMB(77))+F(
     &   42)*(4*(SMB(56)+SMB(67))+SMB(87))+F(65)*(4*SMB(68)+SMB(88))+
     &   F(66)*SMB(89)+2*(F(44)*(2*SMB(27)+SMB(51))+F(48)*(2*SMB(28)+
     &   SMB(52))+F(57)*SMB(53)+F(60)*SMB(73)-Fa(13)*SMB(84)+Fa(14)*S
     &   MB(85)-Fa(15)*SMB(86)-Fa(16)*SMB(90)+Fa(17)*SMB(91)+Fa(18)*S
     &   MB(92))+F(59)*(4*SMB(64)+SMB(93))+F(61)*(4*SMB(65)+SMB(94))+
     &   F(67)*SMB(95))/s23
       result(3) = Fa(19)*SMB(6)-8*(Fa(20)*SMB(11)+Fa(21)*SMB(12)+Fa(2
     &   2)*SMB(13)-Fa(23)*SMB(14)-Fa(24)*SMB(15)-Fa(25)*SMB(16)-Fa(2
     &   6)*SMB(17)-Fa(27)*SMB(18)-Fa(28)*SMB(19))+F(125)*SMB(31)+F(1
     &   27)*SMB(33)+F(128)*SMB(34)+Fa(33)*SMB(59)+4*(F(68)*SMB(1)+F(
     &   69)*SMB(2)+F(70)*SMB(3)+F(71)*SMB(4)+F(72)*SMB(5)-F(78)*SMB(
     &   7)+F(79)*SMB(9)-F(80)*SMB(10)+F(120)*SMB(20)+F(121)*SMB(21)+
     &   F(122)*SMB(22)+F(103)*SMB(24)+F(123)*SMB(26)-F(83)*SMB(36)+F
     &   (87)*SMB(37)+F(130)*SMB(38)-F(118)*SMB(39)-Fa(29)*SMB(42)-Fa
     &   (30)*SMB(43)+Fa(31)*SMB(44)+F(148)*SMB(55)+F(131)*(SMB(40)+S
     &   MB(57))-Fa(32)*SMB(58)+Fa(34)*SMB(60))-F(214)*(8*SMB(8)+SMB(
     &   75)+SMB(76)+SMB(77))+F(189)*SMB(81)+F(190)*SMB(82)+F(191)*SM
     &   B(83)+F(149)*(4*SMB(56)+2*SMB(67)+SMB(87))+F(117)*(2*SMB(68)
     &   +SMB(88))+F(119)*(2*SMB(69)+SMB(89))+2*(F(124)*SMB(30)+F(126
     &   )*SMB(32)+F(129)*SMB(35)+F(81)*SMB(45)+F(85)*SMB(46)+F(143)*
     &   SMB(47)+F(144)*SMB(48)+F(86)*(2*SMB(23)+SMB(49))+F(116)*(2*S
     &   MB(25)+SMB(50))+F(145)*SMB(51)-F(146)*SMB(52)-F(147)*SMB(53)
     &   -F(162)*SMB(61)-F(163)*SMB(62)-F(164)*SMB(63)-F(84)*(4*SMB(5
     &   4)+SMB(70))-F(166)*SMB(71)-F(167)*SMB(72)+F(168)*SMB(73)+Fa(
     &   35)*SMB(74)+Fa(36)*SMB(78)-Fa(37)*SMB(79)-Fa(38)*SMB(80)-Fa(
     &   39)*SMB(84)-Fa(40)*SMB(85)-Fa(41)*SMB(86)-Fa(42)*SMB(90)-Fa(
     &   43)*SMB(91)-Fa(44)*SMB(92))+F(165)*(2*SMB(64)+SMB(93))+F(104
     &   )*(2*SMB(65)+SMB(94))+F(105)*(2*SMB(66)+SMB(95))
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
       m1sq=m*m
       resultborn = -((-2*SMB(1)-2*SMB(2)+2*SMB(3)+2*SMB(4)+2*SMB(5)-2
     &   *(p1mup3+p5mup3+p6mup3)*SMB(6)+SMB(30)+SMB(32)+SMB(35))/((M1
     &   SQ-s16)*s23*(M1SQ-s234)))-(2*SMB(9)+SMB(30)+SMB(32)+SMB(35))
     &   /((M1SQ-s16)*(M1SQ-s234)*s34)
c************************************************************************************
c************************************************************************************
c************************************************************************************
       Return
       End
