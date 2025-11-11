       subroutine HjjCrossIFFFhl1(F1)
       IMPLICIT NONE
       Real*8 P(2),FRe(1:5),FIm(1:5)
       Complex*16 F1(1:5)
       Real*8   p1sq 
       Real*8   p2sq 
       Real*8   p3sq 
       Real*8   p4sq 
       Real*8   p5sq 
       Real*8   s12, s15 
       Real*8   s23 
       Real*8   s34 
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
       Logical PrintB,PrintC,PrintD,PrintE 
       COMMON/HjjCrossIFInvariants/p1sq,p2sq,p3sq,p4sq,p5sq,s12,s23,s3
     &   4,s45,s15
       COMMON/HjjCrossIFPFunctions/P
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
c       Definition of the F,P functions:Independent of the currents    
c************************************************************************************
c************************************************************************************
c************************************************************************************
c************************************************************************************
c************************************************************************************
       FRe(1) = EijR(4,2)+EijR(6,2)-EijR(7,2)-EijR(10,2)
       FIm(1) = EijI(4,2)+EijI(6,2)-EijI(7,2)-EijI(10,2)
       F1(1)=DCMPLX(FRe(1),FIm(1))
       FRe(2) = -EijR(4,1)+EijR(4,2)-EijR(9,2)
       FIm(2) = -EijI(4,1)+EijI(4,2)-EijI(9,2)
       F1(2)=DCMPLX(FRe(2),FIm(2))
       P(1) = p4sq+2*s12-s34-s45
       P(2) = -p4sq-s12+s34+s45
       FRe(3) = s12*(EijR(1,1)-EijR(4,1)+EijR(5,2)-EijR(7,2)-EijR(9,2)
     &   )+EijR(4,2)*P(1)+EijR(10,2)*P(2)
       FIm(3) = s12*(EijI(1,1)-EijI(4,1)+EijI(5,2)-EijI(7,2)-EijI(9,2)
     &   )+EijI(4,2)*P(1)+EijI(10,2)*P(2)
       F1(3)=DCMPLX(FRe(3),FIm(3))
       FRe(4) = -2*(EijR(4,1)-2*(EijR(4,2)-EijR(7,2)))
       FIm(4) = -2*(EijI(4,1)-2*(EijI(4,2)-EijI(7,2)))
       F1(4)=DCMPLX(FRe(4),FIm(4))
       FRe(5) = 2*(EijR(3,1)-EijR(4,1)+2*(EijR(4,2)+EijR(8,2)-EijR(9,2
     &   )-EijR(10,2)))
       FIm(5) = 2*(EijI(3,1)-EijI(4,1)+2*(EijI(4,2)+EijI(8,2)-EijI(9,2
     &   )-EijI(10,2)))
       F1(5)=DCMPLX(FRe(5),FIm(5))
c************************************************************************************
c************************************************************************************
c************************************************************************************
c************************************************************************************
c************************************************************************************
c************************************************************************************
       Return
       End
