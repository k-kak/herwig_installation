c ************************************************************************************
c ************************************************************************************
c Author: Francisco Campanario
c E-mail: francam@particle.uni-karlsruhe.de
c Date: 18/01/2010
c Modified:12/6/2012
c ************************************************************************************
c determine the  finite virtual corrections of 
c psi(p6) ---->--$$$$photon$$$$$$--->---   psi(p5)
c                   |              |              
c                   |              |$$$$$$$$$$$$  p4,mu_p4 
c                   |              |             
c                   |              |$$$$$$$$$$$$  p3, mu_p3 
c                   |              |             
c                   |              |              
c                   |              |             
c barpsi(p1)-->-$$$photon$$$$$$------>---   bar_psi(p2)
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
       subroutine HjjCrossIF66(M1SQ,p1,p2,p3,p4,p5,psi_p1,barpsi_p5,ps
     &   i_p2,barpsi_p3,musq,comp,ngluon,posgluon,result,resultborn)
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
      Complex*16 SMB(10),F(10)
       Complex*16 SMB0(4),SMB11(0:3),SMB12(0:3)
      Complex*16 SMB1(3)
      Real*8 FI(10),FR(10)
      Complex*16 psi_p1_P(2),barpsi_p5_P(2),psi_p2_P(2),barpsi_p3_P(2) 
       Complex*16 psi_p1_M(2),barpsi_p5_M(2),psi_p2_M(2),barpsi_p3_M(2)
       Complex*16 psi_p1(4),barpsi_p5(4),psi_p2(4),barpsi_p3(4)
       Real*8 P(2)
       Complex*16  SC1c,SC1r, SC3crr,SC3ccr,SC3crc,SC3ccc,SC3rrc,SC3rrr
       EXTERNAL    SC1c,SC1r, SC3crr,SC3ccr,SC3crc,SC3ccc,SC3rrc,SC3rrr 
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
      Real*8   p1(0:3),p2(0:3),p3(0:3),p4(0:3),p5(0:3)
       Real*8 dotrr
       Complex*16 A0finG_c,B0finG_c,C0finG_c,D0finG_c,E0finG_c,F0finG_c
       EXTERNAL dotrr,A0finG_c,B0finG_c,C0finG_c,D0finG_c,E0finG_c,F0finG_c
       Complex*16 B0finG,C0finG,D0finG,E0finG,F0finG
       EXTERNAL B0finG,C0finG,D0finG,E0finG,F0finG
        Real*8   p1sq, p1p2, p1p3, p1p4, p1p5 
       Real*8   p2sq, p2p3, p2p4, p2p5 
       Real*8   p3sq, p3p4, p3p5 
       Real*8   p4sq, p4p5 
       Real*8   p5sq 
       Real*8   s12, s13, s14, s15 
       Real*8   s23, s24, s25 
       Real*8   s34, s35 
       Real*8   s45 
       Complex*16  B012,B013,B014,B015 
       Complex*16  B023,B024,B025 
       Complex*16  B034,B035 
       Complex*16  B045 
       Real*8  B012R,B013R,B014R,B015R 
       Real*8  B023R,B024R,B025R 
       Real*8  B034R,B035R 
       Real*8  B045R 
       Real*8  B012I,B013I,B014I,B015I 
       Real*8  B023I,B024I,B025I 
       Real*8  B034I,B035I 
       Real*8  B045I 
       Complex*16 C0123,C0124,C0125 
       Complex*16 C0134,C0135 
       Complex*16 C0145 
       Complex*16 C0234,C0235 
       Complex*16 C0245 
       Complex*16 C0345 
       Real*8 C0123R,C0124R,C0125R 
       Real*8 C0134R,C0135R 
       Real*8 C0145R 
       Real*8 C0234R,C0235R 
       Real*8 C0245R 
       Real*8 C0345R 
       Real*8 C0123I,C0124I,C0125I 
       Real*8 C0134I,C0135I 
       Real*8 C0145I 
       Real*8 C0234I,C0235I 
       Real*8 C0245I 
       Real*8 C0345I 
       Real*8 Cij123R(4,2),Cij124R(4,2),Cij125R(4,2) 
       Real*8 Cij134R(4,2),Cij135R(4,2) 
       Real*8 Cij145R(4,2) 
       Real*8 Cij234R(4,2),Cij235R(4,2) 
       Real*8 Cij245R(4,2) 
       Real*8 Cij345R(4,2) 
       Real*8 Cij123I(4,2),Cij124I(4,2),Cij125I(4,2) 
       Real*8 Cij134I(4,2),Cij135I(4,2) 
       Real*8 Cij145I(4,2) 
       Real*8 Cij234I(4,2),Cij235I(4,2) 
       Real*8 Cij245I(4,2) 
       Real*8 Cij345I(4,2) 
       Complex*16  D01234,D01235 
       Complex*16 D01245 
       Complex*16 D01345 
       Complex*16 D02345 
       Real*8 D01234R,D01235R 
       Real*8 D01245R 
       Real*8 D01345R 
       Real*8 D02345R 
       Real*8 D01234I,D01235I 
       Real*8 D01245I 
       Real*8 D01345I 
       Real*8 D02345I 
       Real*8 Dij1234R(13,3),Dij1235R(13,3) 
       Real*8 Dij1245R(13,3) 
       Real*8 Dij1345R(13,3) 
       Real*8 Dij2345R(13,3) 
       Real*8 Dij1234I(13,3),Dij1235I(13,3) 
       Real*8 Dij1245I(13,3) 
       Real*8 Dij1345I(13,3) 
       Real*8 Dij2345I(13,3) 
       Complex*16 EE0 
       Real*8 EE0R 
       Real*8 EE0I 
       Real*8 EijR(46,4) 
       Real*8 EijI(46,4) 
       Real*8 Invs23MU,Invs34MU 
       Logical PrintB,PrintC,PrintD,PrintE 
      Complex*16 dotrc,dotcc,result(3),resultn,resultborn
       Real*8 musq
      EXTERNAL   dotrc,dotcc
      Integer alpha
       COMMON/HjjCrossIFFhlFunctions/F
      Save/HjjCrossIFFhlFunctions/
       COMMON/HjjCrossIFPFunctions/P
      Save/HjjCrossIFPFunctions/
       COMMON/HjjCrossIFInvariants/p1sq,p2sq,p3sq,p4sq,p5sq,s12,s23,s3
     &   4,s45,s15
       COMMON/HjjCrossIFEVALFunctions/ EE0R,EijR,EE0I,EijI 
       COMMON/HjjCrossIFDVALFunctions/ D01234R,Dij1234R,D01234I,Dij123
     &   4I, D01235R,Dij1235R,D01235I,Dij1235I, D01245R,Dij1245R,D012
     &   45I,Dij1245I, D01345R,Dij1345R,D01345I,Dij1345I, D02345R,Dij
     &   2345R,D02345I,Dij2345I 
       COMMON/HjjCrossIFCVALFunctions/ C0123R,Cij123R,C0123I,Cij123I, 
     &   C0124R,Cij124R,C0124I,Cij124I, C0125R,Cij125R,C0125I,Cij125I
     &   , C0134R,Cij134R,C0134I,Cij134I, C0135R,Cij135R,C0135I,Cij13
     &   5I, C0145R,Cij145R,C0145I,Cij145I, C0234R,Cij234R,C0234I,Cij
     &   234I, C0235R,Cij235R,C0235I,Cij235I, C0245R,Cij245R,C0245I,C
     &   ij245I, C0345R,Cij345R,C0345I,Cij345I 
       COMMON/HjjCrossIFBVALFunctions/ B012R,B012I, B013R,B013I, B014R
     &   ,B014I, B015R,B015I, B023R,B023I, B024R,B024I, B025R,B025I, 
     &   B034R,B034I, B035R,B035I, B045R,B045I
       Integer ngluon, posgluon,Div
       Complex*16 m1sq, m0
       Real*8 m
cc Old routines
       Logical RealRu
       Parameter (RealRu=.false.)
       m0=DCMPLX(0d0,0d0)
c************************************************************************************
c************************************************************************************
c************************************************************************************
c       Definition of the scalar products. Not inlcueded the contraction of the
c       moments with the external currents  
c************************************************************************************
c************************************************************************************
c************************************************************************************
c************************************************************************************
       p1sq = dotrr(p1,p1)
       p1p2 = dotrr(p1,p2)
       p1p3 = dotrr(p1,p3)
       p1p4 = dotrr(p1,p4)
       p1p5 = dotrr(p1,p5)
       p2sq = dotrr(p2,p2)
       p2p3 = dotrr(p2,p3)
       p2p4 = dotrr(p2,p4)
       p2p5 = dotrr(p2,p5)
       p3sq = dotrr(p3,p3)
       p3p4 = dotrr(p3,p4)
       p3p5 = dotrr(p3,p5)
       p4sq = dotrr(p4,p4)
       p4p5 = dotrr(p4,p5)
       p5sq = dotrr(p5,p5)
       s12 = (p1sq +p2sq+ 2*p1p2) 
       s13 = (p1sq +p3sq+ 2*p1p3) 
       s14 = (p1sq +p4sq+ 2*p1p4) 
       s15 = (p1sq +p5sq+ 2*p1p5) 
       s23 = (p2sq +p3sq+ 2*p2p3) 
       s24 = (p2sq +p4sq+ 2*p2p4) 
       s25 = (p2sq +p5sq+ 2*p2p5) 
       s34 = (p3sq +p4sq+ 2*p3p4) 
       s35 = (p3sq +p5sq+ 2*p3p5) 
       s45 = (p4sq +p5sq+ 2*p4p5) 
c       Write(*,'(a5,F20.10)')," p1sq ", p1sq 
c       Write(*,'(a5,F20.10)')," p1p2 ", p1p2
c       Write(*,'(a5,F20.10)')," p1p3 ", p1p3
c       Write(*,'(a5,F20.10)')," p1p4 ", p1p4
c       Write(*,'(a5,F20.10)')," p1p5 ", p1p5
c       Write(*,'(a5,F20.10)')," p2sq ", p2sq 
c       Write(*,'(a5,F20.10)')," p2p3 ", p2p3
c       Write(*,'(a5,F20.10)')," p2p4 ", p2p4
c       Write(*,'(a5,F20.10)')," p2p5 ", p2p5
c       Write(*,'(a5,F20.10)')," p3sq ", p3sq 
c       Write(*,'(a5,F20.10)')," p3p4 ", p3p4
c       Write(*,'(a5,F20.10)')," p3p5 ", p3p5
c       Write(*,'(a5,F20.10)')," p4sq ", p4sq 
c       Write(*,'(a5,F20.10)')," p4p5 ", p4p5
c       Write(*,'(a5,F20.10)')," p5sq ", p5sq 
      Invs23MU=1d0/s23
        Invs34MU=1d0/s34
      PrintB=.False. 
      PrintC=.False. 
      PrintD=.False. 
      PrintE=.False.
c************************************************************************************
c************************************************************************************
        if (comp.gt.0) then    
c    Calling C_ij,D_ij Functions    
c************************************************************************************
c************************************************************************************
c************************************************************************************
c************************************************************************************
c************************************************************************************
       if(RealRu) then
           m=Dble(sqrt(m1sq))
       B012=B0finG(0d0,0d0,p1sq,musq)
       B012R=Dble(B012)
       B012I=DImag(B012)
       B023=B0finG(0d0,0d0,p2sq,musq)
       B023R=Dble(B023)
       B023I=DImag(B023)
       B034=B0finG(0d0,M,p3sq,musq)
       B034R=Dble(B034)
       B034I=DImag(B034)
       B045=B0finG(M,M,p4sq,musq)
       B045R=Dble(B045)
       B045I=DImag(B045)
       B013=B0finG(0d0,0d0,s12,musq)
       B013R=Dble(B013)
       B013I=DImag(B013)
       B014=B0finG(0d0,M,s45,musq)
       B014R=Dble(B014)
       B014I=DImag(B014)
       B024=B0finG(0d0,M,s23,musq)
       B024R=Dble(B024)
       B024I=DImag(B024)
       B025=B0finG(0d0,M,s15,musq)
       B025R=Dble(B025)
       B025I=DImag(B025)
       B035=B0finG(0d0,M,s34,musq)
       B035R=Dble(B035)
       B035I=DImag(B035)
       B015=B0finG(0d0,M,p5sq,musq)
       B015R=Dble(B015)
       B015I=DImag(B015)
        C0123=C0finG(0d0,0d0,0d0,p1sq,p2sq,s12,musq)
       C0123R=Dble(C0123)
       C0123I=DImag(C0123)
       C0124=C0finG(0d0,0d0,M,p1sq,s23,s45,musq)
       C0124R=Dble(C0124)
       C0124I=DImag(C0124)
       C0125=C0finG(0d0,0d0,M,p1sq,s15,p5sq,musq)
       C0125R=Dble(C0125)
       C0125I=DImag(C0125)
       C0234=C0finG(0d0,0d0,M,p2sq,p3sq,s23,musq)
       C0234R=Dble(C0234)
       C0234I=DImag(C0234)
       C0235=C0finG(0d0,0d0,M,p2sq,s34,s15,musq)
       C0235R=Dble(C0235)
       C0235I=DImag(C0235)
       C0345=C0finG(0d0,M,M,p3sq,p4sq,s34,musq)
       C0345R=Dble(C0345)
       C0345I=DImag(C0345)
       C0134=C0finG(0d0,0d0,M,s12,p3sq,s45,musq)
       C0134R=Dble(C0134)
       C0134I=DImag(C0134)
       C0135=C0finG(0d0,0d0,M,s12,s34,p5sq,musq)
       C0135R=Dble(C0135)
       C0135I=DImag(C0135)
       C0145=C0finG(0d0,M,M,s45,p4sq,p5sq,musq)
       C0145R=Dble(C0145)
       C0145I=DImag(C0145)
       C0245=C0finG(0d0,M,M,s23,p4sq,s15,musq)
       C0245R=Dble(C0245)
       C0245I=DImag(C0245)
        D01234=D0finG(0d0,0d0,0d0,M,s12,s23,p1sq,p2sq,p3sq,s45,musq)
       call tens_red4_new_Re_Com_G(0d0,0d0,0d0,M,p1sq,p2sq,p3sq,p1p2,p1p3,p2p3,C0234R,
     &   C0134R,C0124R,C0123R,Cij234R,Cij134R,Cij124R,Cij123R,C0234I,C0134I,C
     &    0124I,C0123I,Cij234I,Cij134I,Cij124I,Cij123I,D01234,D01234R,D01234I,
     &   Dij1234R,Dij1234I)
       D01235=D0finG(0d0,0d0,0d0,M,s12,s15,p1sq,p2sq,s34,p5sq,musq)
       call tens_red4_new_Re_Com_G(0d0,0d0,0d0,M,p1sq,p2sq,s34,p1p2,p1p3+p1p4,p2p3+p2p
     &   4,C0235R,C0135R,C0125R,C0123R,Cij235R,Cij135R,Cij125R,Cij123R,C0235I
     &   ,C0135I,C0125I,C0123I,Cij235I,Cij135I,Cij125I,Cij123I,D01235,D01235R
     &   ,D01235I,Dij1235R,Dij1235I)
       D01245=D0finG(0d0,0d0,M,M,s45,s15,p1sq,s23,p4sq,p5sq,musq)
       call tens_red4_new_Re_Com_G(0d0,0d0,M,M,p1sq,s23,p4sq,p1p2+p1p3,p1p4,p2p4+p3p
     &   4,C0245R,C0145R,C0125R,C0124R,Cij245R,Cij145R,Cij125R,Cij124R,C0245I
     &   ,C0145I,C0125I,C0124I,Cij245I,Cij145I,Cij125I,Cij124I,D01245,D01245R
     &   ,D01245I,Dij1245R,Dij1245I)
       D02345=D0finG(0d0,0d0,M,M,s23,s34,p2sq,p3sq,p4sq,s15,musq)
       call tens_red4_new_Re_Com_G(0d0,0d0,M,M,p2sq,p3sq,p4sq,p2p3,p2p4,p3p4,C0345R,
     &   C0245R,C0235R,C0234R,Cij345R,Cij245R,Cij235R,Cij234R,C0345I,C0245I,C
     &   0235I,C0234I,Cij345I,Cij245I,Cij235I,Cij234I,D02345,D02345R,D02345I,
     &   Dij2345R,Dij2345I)
       D01345=D0finG(0d0,0d0,M,M,s45,s34,s12,p3sq,p4sq,p5sq,musq)
       call tens_red4_new_Re_Com_G(0d0,0d0,M,M,s12,p3sq,p4sq,p1p3+p2p3,p1p4+p2p4,p3p
     &   4,C0345R,C0145R,C0135R,C0134R,Cij345R,Cij145R,Cij135R,Cij134R,C0345I
     &   ,C0145I,C0135I,C0134I,Cij345I,Cij145I,Cij135I,Cij134I,D01345,D01345R
     &   ,D01345I,Dij1345R,Dij1345I)
         EE0=E0finG(0d0,0d0,0d0,M,M,p1sq,p2sq,p3sq,p4sq,p5sq,s12,s23,s34,s45,s15,D02
     &   345,D01345,D01245,D01235,D01234)
       EE0R=Dble(EE0)
       EE0I=DImag(EE0)
       call tens_red5_new_Re_Com_G(0d0,0d0,0d0,M,M,p1sq,p2sq,p3sq,p4sq,p1p2,p1p3,p1p4,p2
     &   p3,p2p4,p3p4,D02345R,D01345R,D01245R,D01235R,D01234R,Dij2345R,
     &   Dij1345R,Dij1245R,Dij1235R,Dij1234R,D02345I,D01345I,D01245I,D01235I,
     &   D01234I,Dij2345I,Dij1345I,Dij1245I,Dij1235I,Dij1234I,EijR,EijI)
            else
c************************************************************************************
c************************************************************************************
       B012=B0finG_c(m0,m0,p1sq,musq)
       B012R=Dble(B012)
       B012I=DImag(B012)
       B023=B0finG_c(m0,m0,p2sq,musq)
       B023R=Dble(B023)
       B023I=DImag(B023)
       B034=B0finG_c(m0,M1SQ,p3sq,musq)
       B034R=Dble(B034)
       B034I=DImag(B034)
       B045=B0finG_c(M1SQ,M1SQ,p4sq,musq)
       B045R=Dble(B045)
       B045I=DImag(B045)
       B013=B0finG_c(m0,m0,s12,musq)
       B013R=Dble(B013)
       B013I=DImag(B013)
       B014=B0finG_c(m0,M1SQ,s45,musq)
       B014R=Dble(B014)
       B014I=DImag(B014)
       B024=B0finG_c(m0,M1SQ,s23,musq)
       B024R=Dble(B024)
       B024I=DImag(B024)
       B025=B0finG_c(m0,M1SQ,s15,musq)
       B025R=Dble(B025)
       B025I=DImag(B025)
       B035=B0finG_c(m0,M1SQ,s34,musq)
       B035R=Dble(B035)
       B035I=DImag(B035)
       B015=B0finG_c(m0,M1SQ,p5sq,musq)
       B015R=Dble(B015)
       B015I=DImag(B015)
        C0123=C0finG_c(m0,m0,m0,p1sq,p2sq,s12,musq)
       C0123R=Dble(C0123)
       C0123I=DImag(C0123)
       C0124=C0finG_c(m0,m0,M1SQ,p1sq,s23,s45,musq)
       C0124R=Dble(C0124)
       C0124I=DImag(C0124)
       C0125=C0finG_c(m0,m0,M1SQ,p1sq,s15,p5sq,musq)
       C0125R=Dble(C0125)
       C0125I=DImag(C0125)
       C0234=C0finG_c(m0,m0,M1SQ,p2sq,p3sq,s23,musq)
       C0234R=Dble(C0234)
       C0234I=DImag(C0234)
       C0235=C0finG_c(m0,m0,M1SQ,p2sq,s34,s15,musq)
       C0235R=Dble(C0235)
       C0235I=DImag(C0235)
       C0345=C0finG_c(m0,M1SQ,M1SQ,p3sq,p4sq,s34,musq)
       C0345R=Dble(C0345)
       C0345I=DImag(C0345)
       C0134=C0finG_c(m0,m0,M1SQ,s12,p3sq,s45,musq)
       C0134R=Dble(C0134)
       C0134I=DImag(C0134)
       C0135=C0finG_c(m0,m0,M1SQ,s12,s34,p5sq,musq)
       C0135R=Dble(C0135)
       C0135I=DImag(C0135)
       C0145=C0finG_c(m0,M1SQ,M1SQ,s45,p4sq,p5sq,musq)
       C0145R=Dble(C0145)
       C0145I=DImag(C0145)
       C0245=C0finG_c(m0,M1SQ,M1SQ,s23,p4sq,s15,musq)
       C0245R=Dble(C0245)
       C0245I=DImag(C0245)
        D01234=D0finG_c(m0,m0,m0,M1SQ,s12,s23,p1sq,p2sq,p3sq,s45,musq)
       call tens_red4_Complex_G(m0,m0,m0,M1SQ,p1sq,p2sq,p3sq,p1p2,p1p3,p2p3,C0234R,
     &   C0134R,C0124R,C0123R,Cij234R,Cij134R,Cij124R,Cij123R,C0234I,C0134I,C
     &   0124I,C0123I,Cij234I,Cij134I,Cij124I,Cij123I,D01234,D01234R,D01234I,
     &   Dij1234R,Dij1234I)
       D01235=D0finG_c(m0,m0,m0,M1SQ,s12,s15,p1sq,p2sq,s34,p5sq,musq)
       call tens_red4_Complex_G(m0,m0,m0,M1SQ,p1sq,p2sq,s34,p1p2,p1p3+p1p4,p2p3+p2p
     &   4,C0235R,C0135R,C0125R,C0123R,Cij235R,Cij135R,Cij125R,Cij123R,C0235I
     &   ,C0135I,C0125I,C0123I,Cij235I,Cij135I,Cij125I,Cij123I,D01235,D01235R
     &   ,D01235I,Dij1235R,Dij1235I)
       D01245=D0finG_c(m0,m0,M1SQ,M1SQ,s45,s15,p1sq,s23,p4sq,p5sq,musq)
       call tens_red4_Complex_G(m0,m0,M1SQ,M1SQ,p1sq,s23,p4sq,p1p2+p1p3,p1p4,p2p4+p3p
     &   4,C0245R,C0145R,C0125R,C0124R,Cij245R,Cij145R,Cij125R,Cij124R,C0245I
     &   ,C0145I,C0125I,C0124I,Cij245I,Cij145I,Cij125I,Cij124I,D01245,D01245R
     &   ,D01245I,Dij1245R,Dij1245I)
       D02345=D0finG_c(m0,m0,M1SQ,M1SQ,s23,s34,p2sq,p3sq,p4sq,s15,musq)
       call tens_red4_Complex_G(m0,m0,M1SQ,M1SQ,p2sq,p3sq,p4sq,p2p3,p2p4,p3p4,C0345R,
     &   C0245R,C0235R,C0234R,Cij345R,Cij245R,Cij235R,Cij234R,C0345I,C0245I,C
     &   0235I,C0234I,Cij345I,Cij245I,Cij235I,Cij234I,D02345,D02345R,D02345I,
     &   Dij2345R,Dij2345I)
       D01345=D0finG_c(m0,m0,M1SQ,M1SQ,s45,s34,s12,p3sq,p4sq,p5sq,musq)
       call tens_red4_Complex_G(m0,m0,M1SQ,M1SQ,s12,p3sq,p4sq,p1p3+p2p3,p1p4+p2p4,p3p
     &   4,C0345R,C0145R,C0135R,C0134R,Cij345R,Cij145R,Cij135R,Cij134R,C0345I
     &   ,C0145I,C0135I,C0134I,Cij345I,Cij145I,Cij135I,Cij134I,D01345,D01345R
     &   ,D01345I,Dij1345R,Dij1345I)
        EE0=E0finG_c(m0,m0,m0,M1SQ,M1SQ,p1sq,p2sq,p3sq,p4sq,p5sq,s12,s23,s34,s45,s15,D02
     &   345,D01345,D01245,D01235,D01234)
       EE0R=Dble(EE0)
       EE0I=DImag(EE0)
       call tens_red5_Complex_G(m0,m0,m0,M1SQ,M1SQ,p1sq,p2sq,p3sq,p4sq,p1p2,p1p3,p1p4,p2
     &   p3,p2p4,p3p4,D02345R,D01345R,D01245R,D01235R,D01234R,Dij2345R,
     &   Dij1345R,Dij1245R,Dij1235R,Dij1234R,D02345I,D01345I,D01245I,D01235I,
     &   D01234I,Dij2345I,Dij1345I,Dij1245I,Dij1235I,Dij1234I,EijR,EijI)
       endif
c       Definition of the F,P functions:Independent of the currents    
c************************************************************************************
c************************************************************************************
c************************************************************************************
       call HjjCrossIFFFhl1(F(1))
       call HjjCrossIFFFhl2(F(6))
c************************************************************************************
c************************************************************************************
c************************************************************************************
       endif  
c               PART THAT DEPENDS ON THE EXTERNAL CURRENT
c************************************************************************************
c************************************************************************************
c************************************************************************************
c************** Calling the Fa functions**********************************************************************
c************************************************************************************
c************************************************************************************
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
       psi_p1_M(1)=psi_p1(1)
       psi_p1_M(2)=psi_p1(2)
       psi_p1_P(1)=psi_p1(3)
       psi_p1_P(2)=psi_p1(4)
       psi_p2_M(1)=psi_p2(1)
       psi_p2_M(2)=psi_p2(2)
       psi_p2_P(1)=psi_p2(3)
       psi_p2_P(2)=psi_p2(4)
       barpsi_p5_P(1)=barpsi_p5(1)
       barpsi_p5_P(2)=barpsi_p5(2)
       barpsi_p5_M(1)=barpsi_p5(3)
       barpsi_p5_M(2)=barpsi_p5(4)
       barpsi_p3_P(1)=barpsi_p3(1)
       barpsi_p3_P(2)=barpsi_p3(2)
       barpsi_p3_M(1)=barpsi_p3(3)
       barpsi_p3_M(2)=barpsi_p3(4)
c************************************************************************************
c************************************************************************************
       SMB0(1) = SC1r(barpsi_p3_P,p1,psi_p2_P,+1)
       SMB0(2) = SC1r(barpsi_p5_P,p3,psi_p1_P,+1)
       SMB0(3) = SC1r(barpsi_p3_P,p5,psi_p2_P,+1)
       SMB0(4) = SC1r(barpsi_p5_P,p2,psi_p1_P,+1)
c************************************************************************************
c************************************************************************************
      do i=0,3
        v1(0)=delta(i,0)
        v1(1)=delta(i,1)
        v1(2)=delta(i,2)
        v1(3)=delta(i,3)
       SMB11(i) = SC1c(barpsi_p3_P,v1,psi_p2_P,+1)
       SMB12(i) = SC1c(barpsi_p5_P,v1,psi_p1_P,+1)
      enddo
c************************************************************************************
c************************************************************************************
c************************************************************************************
c************************************************************************************
c************************************************************************************
c************************************************************************************
       SMB1(2) = SMB0(1)*SMB0(2)
       SMB1(3) = SMB0(3)*SMB0(4)
c************************************************************************************
       SMB1(1) = dotcc(SMB11,SMB12)
c************************************************************************************
c************************************************************************************
       SMB(1) = SMB1(2)
       SMB(2) = SMB1(3)
       SMB(3) = SMB1(1)
       SMB(4) = s15*SMB1(1)
       SMB(5) = s23*SMB1(1)
       SMB(6) = 16*SMB1(1)
       SMB(7) = 2*s12*SMB1(1)
       SMB(8) = -2*(s12+s23-s45)*SMB1(1)
       SMB(9) = -2*(s12+s15-s34)*SMB1(1)
       SMB(10) = 2*(p4sq+s12-s34-s45)*SMB1(1)
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
c The born factor is 1. The virtual (-I)
c . I have multiplied by (-I) to eliminated (i) factors.
c The factorization from the B_ij is Fact=(I/(4\[Pi])^2 (4 \[Pi])^Eps Gamma[1+Eps] (musq)^(-Eps))
c  c So, I*(I)=(-1)!!!
       result(1) = 4*(F(1)*SMB(1)+F(2)*SMB(2))+2*F(3)*SMB(3)+F(4)*SMB(
     &   4)+F(5)*SMB(5)+F(6)*SMB(6)+F(7)*SMB(7)+F(8)*SMB(8)+F(9)*SMB(
     &   9)+F(10)*SMB(10)
       result(1) =-result(1)
       resultborn = SMB(3)/((s15-M1SQ)*(s23-M1SQ))
c************************************************************************************
c************************************************************************************
c************************************************************************************
       Return
       End
