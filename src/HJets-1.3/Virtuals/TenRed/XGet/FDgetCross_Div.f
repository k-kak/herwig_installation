       subroutine FgetCross_Div(M,p1,p2,p3,p4,p5,p6,musq,id)
       IMPLICIT NONE
       integer id
       Real*8   p1(0:3),p2(0:3),p3(0:3),p4(0:3),p5(0:3),p6(0:3)
       Complex*16 A0finGDiv,B0finGDiv,C0finGDiv,D0finGDiv,E0finG,F0finG
       EXTERNAL A0finGDiv,B0finGDiv,C0finGDiv,D0finGDiv,E0finG,F0finG
       real*8   dotrr
       external dotrr
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
       COMMON/H3jcrossInvariants/p1sq,p2sq,p3sq,p4sq,p5sq,p6sq,s12,s23,s34,
     &   s45,s56,s16,s123,s234,s345
       COMMON/H3jcrossFVALFunctions/F0123456R,F123456R,F0123456I,F123456I
       COMMON/H3jcrossEVALFunctions/ E012345R,E12345R,E012345I,E12345I, E01
     &   2346R,E12346R,E012346I,E12346I, E012356R,E12356R,E012356I,E1
     &   2356I, E012456R,E12456R,E012456I,E12456I, E013456R,E13456R,E
     &   013456I,E13456I, E023456R,E23456R,E023456I,E23456I
       COMMON/H3jcrossDVALFunctions/ D01234R,D1234R,D01234I,D1234I, D01235R
     &   ,D1235R,D01235I,D1235I, D01236R,D1236R,D01236I,D1236I, D0124
     &   5R,D1245R,D01245I,D1245I, D01246R,D1246R,D01246I,D1246I, D01
     &   256R,D1256R,D01256I,D1256I, D01345R,D1345R,D01345I,D1345I, D
     &   01346R,D1346R,D01346I,D1346I, D01356R,D1356R,D01356I,D1356I,
     &    D01456R,D1456R,D01456I,D1456I, D02345R,D2345R,D02345I,D2345
     &   I, D02346R,D2346R,D02346I,D2346I, D02356R,D2356R,D02356I,D23
     &   56I, D02456R,D2456R,D02456I,D2456I, D03456R,D3456R,D03456I,D
     &   3456I
       COMMON/H3jcrossCVALFunctions/ C0123R,C123R,C0123I,C123I, C0124R,C124
     &   R,C0124I,C124I, C0125R,C125R,C0125I,C125I, C0126R,C126R,C012
     &   6I,C126I, C0134R,C134R,C0134I,C134I, C0135R,C135R,C0135I,C13
     &   5I, C0136R,C136R,C0136I,C136I, C0145R,C145R,C0145I,C145I, C0
     &   146R,C146R,C0146I,C146I, C0156R,C156R,C0156I,C156I, C0234R,C
     &   234R,C0234I,C234I, C0235R,C235R,C0235I,C235I, C0236R,C236R,C
     &   0236I,C236I, C0245R,C245R,C0245I,C245I, C0246R,C246R,C0246I,
     &   C246I, C0256R,C256R,C0256I,C256I, C0345R,C345R,C0345I,C345I,
     &    C0346R,C346R,C0346I,C346I, C0356R,C356R,C0356I,C356I, C0456
     &   R,C456R,C0456I,C456I
       COMMON/H3jcrossBVALFunctions/ B012R,B012I, B013R,B013I, B014R,B014I,
     &    B015R,B015I, B016R,B016I, B023R,B023I, B024R,B024I, B025R,B
     &   025I, B026R,B026I, B034R,B034I, B035R,B035I, B036R,B036I, B0
     &   45R,B045I, B046R,B046I, B056R,B056I
       
       real*8 M,musq


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

c    Calling C_ij,D_ij,E_ij,F_ij Functions    
c************************************************************************************
c************************************************************************************
      A01=A0finGDiv(0d0,musq,id)
      A01R=Dble(A01)
      A01I=DImag(A01)
      A02=A0finGDiv(0d0,musq,id)
      A02R=Dble(A02)
      A02I=DImag(A02)
      A03=A0finGDiv(0d0,musq,id)
      A03R=Dble(A03)
      A03I=DImag(A03)
      A04=A0finGDiv(0d0,musq,id)
      A04R=Dble(A04)
      A04I=DImag(A04)
      A05=A0finGDiv(M,musq,id)
      A05R=Dble(A05)
      A05I=DImag(A05)
      A06=A0finGDiv(M,musq,id)
      A06R=Dble(A06)
      A06I=DImag(A06)

      B012=B0finGDiv(0d0,0d0,p1sq,musq,id)
      B012R=Dble(B012)
      B012I=DImag(B012)
      call tens_red2_new_Re_ComDiv(0d0,0d0,p1sq,A02R,A02I,A01R,A01I,B012R,B
     &   012I,Bij12R,Bij12I)
      B023=B0finGDiv(0d0,0d0,p2sq,musq,id)
      B023R=Dble(B023)
      B023I=DImag(B023)
      call tens_red2_new_Re_ComDiv(0d0,0d0,p2sq,A03R,A03I,A02R,A02I,B023R,B
     &   023I,Bij23R,Bij23I)
      B034=B0finGDiv(0d0,0d0,p3sq,musq,id)
      B034R=Dble(B034)
      B034I=DImag(B034)
      call tens_red2_new_Re_ComDiv(0d0,0d0,p3sq,A04R,A04I,A03R,A03I,B034R,B
     &   034I,Bij34R,Bij34I)
      B045=B0finGDiv(0d0,M,p4sq,musq,id)
      B045R=Dble(B045)
      B045I=DImag(B045)
      call tens_red2_new_Re_ComDiv(0d0,M,p4sq,A05R,A05I,A04R,A04I,B045R,B
     &   045I,Bij45R,Bij45I)
      B056=B0finGDiv(M,M,p5sq,musq,id)
      B056R=Dble(B056)
      B056I=DImag(B056)
      call tens_red2_new_Re_ComDiv(M,M,p5sq,A06R,A06I,A05R,A05I,B056R,B
     &   056I,Bij56R,Bij56I)
      B016=B0finGDiv(0d0,M,p6sq,musq,id)
      B016R=Dble(B016)
      B016I=DImag(B016)
      call tens_red2_new_Re_ComDiv(0d0,M,p6sq,A06R,A06I,A01R,A01I,B016R,B
     &   016I,Bij16R,Bij16I)
      B013=B0finGDiv(0d0,0d0,s12,musq,id)
      B013R=Dble(B013)
      B013I=DImag(B013)
      call tens_red2_new_Re_ComDiv(0d0,0d0,s12,A03R,A03I,A01R,A01I,B013R,B
     &   013I,Bij13R,Bij13I)
      B014=B0finGDiv(0d0,0d0,s123,musq,id)
      B014R=Dble(B014)
      B014I=DImag(B014)
      call tens_red2_new_Re_ComDiv(0d0,0d0,s123,A04R,A04I,A01R,A01I,B014R,B
     &   014I,Bij14R,Bij14I)
      B026=B0finGDiv(0d0,M,s16,musq,id)
      B026R=Dble(B026)
      B026I=DImag(B026)
      call tens_red2_new_Re_ComDiv(0d0,M,s16,A06R,A06I,A02R,A02I,B026R,B
     &   026I,Bij26R,Bij26I)
      B024=B0finGDiv(0d0,0d0,s23,musq,id)
      B024R=Dble(B024)
      B024I=DImag(B024)
      call tens_red2_new_Re_ComDiv(0d0,0d0,s23,A04R,A04I,A02R,A02I,B024R,B
     &   024I,Bij24R,Bij24I)
      B025=B0finGDiv(0d0,M,s234,musq,id)
      B025R=Dble(B025)
      B025I=DImag(B025)
      call tens_red2_new_Re_ComDiv(0d0,M,s234,A05R,A05I,A02R,A02I,B025R,B
     &   025I,Bij25R,Bij25I)
      B035=B0finGDiv(0d0,M,s34,musq,id)
      B035R=Dble(B035)
      B035I=DImag(B035)
      call tens_red2_new_Re_ComDiv(0d0,M,s34,A05R,A05I,A03R,A03I,B035R,B
     &   035I,Bij35R,Bij35I)
      B036=B0finGDiv(0d0,M,s345,musq,id)
      B036R=Dble(B036)
      B036I=DImag(B036)
      call tens_red2_new_Re_ComDiv(0d0,M,s345,A06R,A06I,A03R,A03I,B036R,B
     &   036I,Bij36R,Bij36I)
      B046=B0finGDiv(0d0,M,s45,musq,id)
      B046R=Dble(B046)
      B046I=DImag(B046)
      call tens_red2_new_Re_ComDiv(0d0,M,s45,A06R,A06I,A04R,A04I,B046R,B
     &   046I,Bij46R,Bij46I)
      B015=B0finGDiv(0d0,M,s56,musq,id)
      B015R=Dble(B015)
      B015I=DImag(B015)
      call tens_red2_new_Re_ComDiv(0d0,M,s56,A05R,A05I,A01R,A01I,B015R,B
     &   015I,Bij15R,Bij15I)

      C0123=C0finGDiv(0d0,0d0,0d0,p1sq,p2sq,s12,musq,id)
      call tens_red3_new_Re_Com_GDiv(0d0,0d0,0d0,p1sq,p2sq,s12,
     & B023R,B013R,B012R,B023I,B013I,B012I,Bij23R,Bij13R,Bij12R,Bij23I,Bij13I,Bij12I,
     & C0123,C0123R,C0123I,C123R,C123I)
      C0126=C0finGDiv(0d0,0d0,M,p1sq,s16,p6sq,musq,id)
      call tens_red3_new_Re_Com_GDiv(0d0,0d0,M,p1sq,s16,p6sq,
     & B026R,B016R,B012R,B026I,B016I,B012I,Bij26R,Bij16R,Bij12R,Bij26I,Bij16I,Bij12I,
     & C0126,C0126R,C0126I,C126R,C126I)
      C0124=C0finGDiv(0d0,0d0,0d0,p1sq,s23,s123,musq,id)
      call tens_red3_new_Re_Com_GDiv(0d0,0d0,0d0,p1sq,s23,s123,
     & B024R,B014R,B012R,B024I,B014I,B012I,Bij24R,Bij14R,Bij12R,Bij24I,Bij14I,Bij12I,
     & C0124,C0124R,C0124I,C124R,C124I)
      C0125=C0finGDiv(0d0,0d0,M,p1sq,s234,s56,musq,id)
      call tens_red3_new_Re_Com_GDiv(0d0,0d0,M,p1sq,s234,s56,
     & B025R,B015R,B012R,B025I,B015I,B012I,Bij25R,Bij15R,Bij12R,Bij25I,Bij15I,Bij12I,
     & C0125,C0125R,C0125I,C125R,C125I)
      C0234=C0finGDiv(0d0,0d0,0d0,p2sq,p3sq,s23,musq,id)
      call tens_red3_new_Re_Com_GDiv(0d0,0d0,0d0,p2sq,p3sq,s23,
     & B034R,B024R,B023R,B034I,B024I,B023I,Bij34R,Bij24R,Bij23R,Bij34I,Bij24I,Bij23I,
     & C0234,C0234R,C0234I,C234R,C234I)
      C0235=C0finGDiv(0d0,0d0,M,p2sq,s34,s234,musq,id)
      call tens_red3_new_Re_Com_GDiv(0d0,0d0,M,p2sq,s34,s234,
     & B035R,B025R,B023R,B035I,B025I,B023I,Bij35R,Bij25R,Bij23R,Bij35I,Bij25I,Bij23I,
     & C0235,C0235R,C0235I,C235R,C235I)
      C0236=C0finGDiv(0d0,0d0,M,p2sq,s345,s16,musq,id)
      call tens_red3_new_Re_Com_GDiv(0d0,0d0,M,p2sq,s345,s16,
     & B036R,B026R,B023R,B036I,B026I,B023I,Bij36R,Bij26R,Bij23R,Bij36I,Bij26I,Bij23I,
     & C0236,C0236R,C0236I,C236R,C236I)
      C0345=C0finGDiv(0d0,0d0,M,p3sq,p4sq,s34,musq,id)
      call tens_red3_new_Re_Com_GDiv(0d0,0d0,M,p3sq,p4sq,s34,
     & B045R,B035R,B034R,B045I,B035I,B034I,Bij45R,Bij35R,Bij34R,Bij45I,Bij35I,Bij34I,
     & C0345,C0345R,C0345I,C345R,C345I)
      C0346=C0finGDiv(0d0,0d0,M,p3sq,s45,s345,musq,id)
      call tens_red3_new_Re_Com_GDiv(0d0,0d0,M,p3sq,s45,s345,
     & B046R,B036R,B034R,B046I,B036I,B034I,Bij46R,Bij36R,Bij34R,Bij46I,Bij36I,Bij34I,
     & C0346,C0346R,C0346I,C346R,C346I)
      C0456=C0finGDiv(0d0,M,M,p4sq,p5sq,s45,musq,id)
      call tens_red3_new_Re_Com_GDiv(0d0,M,M,p4sq,p5sq,s45,
     & B056R,B046R,B045R,B056I,B046I,B045I,Bij56R,Bij46R,Bij45R,Bij56I,Bij46I,Bij45I,
     & C0456,C0456R,C0456I,C456R,C456I)
      C0134=C0finGDiv(0d0,0d0,0d0,s12,p3sq,s123,musq,id)
      call tens_red3_new_Re_Com_GDiv(0d0,0d0,0d0,s12,p3sq,s123,
     & B034R,B014R,B013R,B034I,B014I,B013I,Bij34R,Bij14R,Bij13R,Bij34I,Bij14I,Bij13I,
     & C0134,C0134R,C0134I,C134R,C134I)
      C0135=C0finGDiv(0d0,0d0,M,s12,s34,s56,musq,id)
      call tens_red3_new_Re_Com_GDiv(0d0,0d0,M,s12,s34,s56,
     & B035R,B015R,B013R,B035I,B015I,B013I,Bij35R,Bij15R,Bij13R,Bij35I,Bij15I,Bij13I,
     & C0135,C0135R,C0135I,C135R,C135I)
      C0136=C0finGDiv(0d0,0d0,M,s12,s345,p6sq,musq,id)
      call tens_red3_new_Re_Com_GDiv(0d0,0d0,M,s12,s345,p6sq,
     & B036R,B016R,B013R,B036I,B016I,B013I,Bij36R,Bij16R,Bij13R,Bij36I,Bij16I,Bij13I,
     & C0136,C0136R,C0136I,C136R,C136I)
      C0145=C0finGDiv(0d0,0d0,M,s123,p4sq,s56,musq,id)
      call tens_red3_new_Re_Com_GDiv(0d0,0d0,M,s123,p4sq,s56,
     & B045R,B015R,B014R,B045I,B015I,B014I,Bij45R,Bij15R,Bij14R,Bij45I,Bij15I,Bij14I,
     & C0145,C0145R,C0145I,C145R,C145I)
      C0146=C0finGDiv(0d0,0d0,M,s123,s45,p6sq,musq,id)
      call tens_red3_new_Re_Com_GDiv(0d0,0d0,M,s123,s45,p6sq,
     & B046R,B016R,B014R,B046I,B016I,B014I,Bij46R,Bij16R,Bij14R,Bij46I,Bij16I,Bij14I,
     & C0146,C0146R,C0146I,C146R,C146I)
      C0245=C0finGDiv(0d0,0d0,M,s23,p4sq,s234,musq,id)
      call tens_red3_new_Re_Com_GDiv(0d0,0d0,M,s23,p4sq,s234,
     & B045R,B025R,B024R,B045I,B025I,B024I,Bij45R,Bij25R,Bij24R,Bij45I,Bij25I,Bij24I,
     & C0245,C0245R,C0245I,C245R,C245I)
      C0246=C0finGDiv(0d0,0d0,M,s23,s45,s16,musq,id)
      call tens_red3_new_Re_Com_GDiv(0d0,0d0,M,s23,s45,s16,
     & B046R,B026R,B024R,B046I,B026I,B024I,Bij46R,Bij26R,Bij24R,Bij46I,Bij26I,Bij24I,
     & C0246,C0246R,C0246I,C246R,C246I)
      C0256=C0finGDiv(0d0,M,M,s234,p5sq,s16,musq,id)
      call tens_red3_new_Re_Com_GDiv(0d0,M,M,s234,p5sq,s16,
     & B056R,B026R,B025R,B056I,B026I,B025I,Bij56R,Bij26R,Bij25R,Bij56I,Bij26I,Bij25I,
     & C0256,C0256R,C0256I,C256R,C256I)
      C0356=C0finGDiv(0d0,M,M,s34,p5sq,s345,musq,id)
      call tens_red3_new_Re_Com_GDiv(0d0,M,M,s34,p5sq,s345,
     & B056R,B036R,B035R,B056I,B036I,B035I,Bij56R,Bij36R,Bij35R,Bij56I,Bij36I,Bij35I,
     & C0356,C0356R,C0356I,C356R,C356I)
      C0156=C0finGDiv(0d0,M,M,s56,p5sq,p6sq,musq,id)
      call tens_red3_new_Re_Com_GDiv(0d0,M,M,s56,p5sq,p6sq,
     & B056R,B016R,B015R,B056I,B016I,B015I,Bij56R,Bij16R,Bij15R,Bij56I,Bij16I,Bij15I,
     & C0156,C0156R,C0156I,C156R,C156I)

      D01234=D0finGDiv(0d0,0d0,0d0,0d0,s12,s23,p1sq,p2sq,p3sq,s123,musq,id)
      call tens_red4_new_Re_Com_GDiv(0d0,0d0,0d0,0d0,p1sq,p2sq,p3sq,p1p2,p1p3,p2p3,C0234R,
     &   C0134R,C0124R,C0123R,C234R,C134R,C124R,C123R,C0234I,C0134I,C
     &   0124I,C0123I,C234I,C134I,C124I,C123I,D01234,D01234R,D01234I,
     &   D1234R,D1234I)
      D01235=D0finGDiv(0d0,0d0,0d0,M,s12,s234,p1sq,p2sq,s34,s56,musq,id)
      call tens_red4_new_Re_Com_GDiv(0d0,0d0,0d0,M,p1sq,p2sq,s34,p1p2,p1p3+p1p4,p2p3+p2p
     &   4,C0235R,C0135R,C0125R,C0123R,C235R,C135R,C125R,C123R,C0235I
     &   ,C0135I,C0125I,C0123I,C235I,C135I,C125I,C123I,D01235,D01235R
     &   ,D01235I,D1235R,D1235I)
      D01236=D0finGDiv(0d0,0d0,0d0,M,s12,s16,p1sq,p2sq,s345,p6sq,musq,id)
      call tens_red4_new_Re_Com_GDiv(0d0,0d0,0d0,M,p1sq,p2sq,s345,p1p2,p1p3+p1p4+p1p5,p2
     &   p3+p2p4+p2p5,C0236R,C0136R,C0126R,C0123R,C236R,C136R,C126R,C
     &   123R,C0236I,C0136I,C0126I,C0123I,C236I,C136I,C126I,C123I,D01
     &   236,D01236R,D01236I,D1236R,D1236I)
      D01245=D0finGDiv(0d0,0d0,0d0,M,s123,s234,p1sq,s23,p4sq,s56,musq,id)
      call tens_red4_new_Re_Com_GDiv(0d0,0d0,0d0,M,p1sq,s23,p4sq,p1p2+p1p3,p1p4,p2p4+p3p
     &   4,C0245R,C0145R,C0125R,C0124R,C245R,C145R,C125R,C124R,C0245I
     &   ,C0145I,C0125I,C0124I,C245I,C145I,C125I,C124I,D01245,D01245R
     &   ,D01245I,D1245R,D1245I)
      D01246=D0finGDiv(0d0,0d0,0d0,M,s123,s16,p1sq,s23,s45,p6sq,musq,id)
      call tens_red4_new_Re_Com_GDiv(0d0,0d0,0d0,M,p1sq,s23,s45,p1p2+p1p3,p1p4+p1p5,p2p4
     &   +p2p5+p3p4+p3p5,C0246R,C0146R,C0126R,C0124R,C246R,C146R,C126
     &   R,C124R,C0246I,C0146I,C0126I,C0124I,C246I,C146I,C126I,C124I,
     &   D01246,D01246R,D01246I,D1246R,D1246I)
      D01256=D0finGDiv(0d0,0d0,M,M,s56,s16,p1sq,s234,p5sq,p6sq,musq,id)
      call tens_red4_new_Re_Com_GDiv(0d0,0d0,M,M,p1sq,s234,p5sq,p1p2+p1p3+p1p4,p1p5,p2
     &   p5+p3p5+p4p5,C0256R,C0156R,C0126R,C0125R,C256R,C156R,C126R,C
     &   125R,C0256I,C0156I,C0126I,C0125I,C256I,C156I,C126I,C125I,D01
     &   256,D01256R,D01256I,D1256R,D1256I)
      D02345=D0finGDiv(0d0,0d0,0d0,M,s23,s34,p2sq,p3sq,p4sq,s234,musq,id)
      call tens_red4_new_Re_Com_GDiv(0d0,0d0,0d0,M,p2sq,p3sq,p4sq,p2p3,p2p4,p3p4,C0345R,
     &   C0245R,C0235R,C0234R,C345R,C245R,C235R,C234R,C0345I,C0245I,C
     &   0235I,C0234I,C345I,C245I,C235I,C234I,D02345,D02345R,D02345I,
     &   D2345R,D2345I)
      D02346=D0finGDiv(0d0,0d0,0d0,M,s23,s345,p2sq,p3sq,s45,s16,musq,id)
      call tens_red4_new_Re_Com_GDiv(0d0,0d0,0d0,M,p2sq,p3sq,s45,p2p3,p2p4+p2p5,p3p4+p3p
     &   5,C0346R,C0246R,C0236R,C0234R,C346R,C246R,C236R,C234R,C0346I
     &   ,C0246I,C0236I,C0234I,C346I,C246I,C236I,C234I,D02346,D02346R
     &   ,D02346I,D2346R,D2346I)
      D02356=D0finGDiv(0d0,0d0,M,M,s234,s345,p2sq,s34,p5sq,s16,musq,id)
      call tens_red4_new_Re_Com_GDiv(0d0,0d0,M,M,p2sq,s34,p5sq,p2p3+p2p4,p2p5,p3p5+p4p
     &   5,C0356R,C0256R,C0236R,C0235R,C356R,C256R,C236R,C235R,C0356I
     &   ,C0256I,C0236I,C0235I,C356I,C256I,C236I,C235I,D02356,D02356R
     &   ,D02356I,D2356R,D2356I)
      D03456=D0finGDiv(0d0,0d0,M,M,s34,s45,p3sq,p4sq,p5sq,s345,musq,id)
      call tens_red4_new_Re_Com_GDiv(0d0,0d0,M,M,p3sq,p4sq,p5sq,p3p4,p3p5,p4p5,C0456R,
     &   C0356R,C0346R,C0345R,C456R,C356R,C346R,C345R,C0456I,C0356I,C
     &   0346I,C0345I,C456I,C356I,C346I,C345I,D03456,D03456R,D03456I,
     &   D3456R,D3456I)
      D01345=D0finGDiv(0d0,0d0,0d0,M,s123,s34,s12,p3sq,p4sq,s56,musq,id)
      call tens_red4_new_Re_Com_GDiv(0d0,0d0,0d0,M,s12,p3sq,p4sq,p1p3+p2p3,p1p4+p2p4,p3p
     &   4,C0345R,C0145R,C0135R,C0134R,C345R,C145R,C135R,C134R,C0345I
     &   ,C0145I,C0135I,C0134I,C345I,C145I,C135I,C134I,D01345,D01345R
     &   ,D01345I,D1345R,D1345I)
      D01346=D0finGDiv(0d0,0d0,0d0,M,s123,s345,s12,p3sq,s45,p6sq,musq,id)
      call tens_red4_new_Re_Com_GDiv(0d0,0d0,0d0,M,s12,p3sq,s45,p1p3+p2p3,p1p4+p1p5+p2p4
     &   +p2p5,p3p4+p3p5,C0346R,C0146R,C0136R,C0134R,C346R,C146R,C136
     &   R,C134R,C0346I,C0146I,C0136I,C0134I,C346I,C146I,C136I,C134I,
     &   D01346,D01346R,D01346I,D1346R,D1346I)
      D01356=D0finGDiv(0d0,0d0,M,M,s56,s345,s12,s34,p5sq,p6sq,musq,id)
      call tens_red4_new_Re_Com_GDiv(0d0,0d0,M,M,s12,s34,p5sq,p1p3+p1p4+p2p3+p2p4,p1p5
     &   +p2p5,p3p5+p4p5,C0356R,C0156R,C0136R,C0135R,C356R,C156R,C136
     &   R,C135R,C0356I,C0156I,C0136I,C0135I,C356I,C156I,C136I,C135I,
     &   D01356,D01356R,D01356I,D1356R,D1356I)
      D01456=D0finGDiv(0d0,0d0,M,M,s56,s45,s123,p4sq,p5sq,p6sq,musq,id)
      call tens_red4_new_Re_Com_GDiv(0d0,0d0,M,M,s123,p4sq,p5sq,p1p4+p2p4+p3p4,p1p5+p2
     &   p5+p3p5,p4p5,C0456R,C0156R,C0146R,C0145R,C456R,C156R,C146R,C
     &   145R,C0456I,C0156I,C0146I,C0145I,C456I,C156I,C146I,C145I,D01
     &   456,D01456R,D01456I,D1456R,D1456I)
      D02456=D0finGDiv(0d0,0d0,M,M,s234,s45,s23,p4sq,p5sq,s16,musq,id)
      call tens_red4_new_Re_Com_GDiv(0d0,0d0,M,M,s23,p4sq,p5sq,p2p4+p3p4,p2p5+p3p5,p4p
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
c***************************************



      return
      end
