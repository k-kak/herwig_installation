c************** Calling the Fa functions*************************
       subroutine H3jCrossFa1(p1mup3,p1mup4,p2mup3,p2mup4,p3mup3,p3mup
     &   4,p4mup3,p4mup4,p5mup3,p5mup4,p6mup3,p6mup4,mup3mup4,Fa1)
       IMPLICIT NONE
      Complex*16   p1mup2, p1mup3, p1mup4, p1mup6, p2mup2, p2mup3, 
     -          p2mup4, p2mup6, p3mup2, p3mup3, p3mup4, p3mup6, 
     -          p4mup2, p4mup3, p4mup4, p4mup6, p5mup2, p5mup3, 
     -          p5mup4, p5mup6, p6mup2, p6mup3, p6mup4, p6mup6
       Complex*16   mup2mup3, mup2mup4, mup2mup6, mup3mup4, mup3mup6, 
     -          mup4mup6
        common/H3jCrossFhlFunctions/F
       COMMON/H3jCrossInvariants/p1sq,p2sq,p3sq,p4sq,p5sq,p6sq,s12,s23
     &   ,s34,s45,s56,s16,s123,s234,s345
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
      Complex*16 Fa(44),F(218)
      Real*8 P(111) 
        COMMON/H3jCrossPFunctions/P
       Complex*16 Fa1(1:22)
       COMMON/H3jCrossFaFunctions/Fa
       Fa1(1) = p1mup3*F(24)-p5mup3*F(25)+p6mup3*F(26)-4*p2mup3*F(27)
       Fa1(2) = 16*p6mup3*F(34)-8*p5mup3*F(35)
       Fa1(3) = p5mup3*F(35)-p2mup3*F(36)+p6mup3*F(37)
       Fa1(4) = p1mup3*F(38)-p5mup3*F(39)+p6mup3*F(40)
       Fa1(5) = p2mup3*F(39)-p1mup3*F(41)+p6mup3*F(42)
       Fa1(6) = p1mup3*F(43)-p5mup3*F(44)
       Fa1(7) = -(p5mup3*F(42))+p2mup3*F(45)+p1mup3*F(46)
       Fa1(8) = p2mup3*F(51)-p5mup3*F(52)+p6mup3*F(53)
       Fa1(9) = 4*p2mup3*F(54)+p5mup3*F(55)
       Fa1(10) = p5mup3*F(42)-p1mup3*F(53)+p2mup3*F(56)
       Fa1(11) = p5mup3*F(39)-p1mup3*F(61)+p6mup3*F(62)
       Fa1(12) = p1mup3+p5mup3+p6mup3
       Fa1(13) = p5mup3*F(52)+p6mup3*F(57)
       Fa1(14) = -(p6mup3*F(42))+p1mup3*F(52)
       Fa1(15) = p5mup3*F(42)+p1mup3*F(57)
       Fa1(16) = p5mup3*F(51)+p6mup3*F(54)
       Fa1(17) = p1mup3*F(51)+p6mup3*F(56)
       Fa1(18) = p1mup3*F(54)-p5mup3*F(56)
       Fa1(19) = p1mup3*F(68)+p5mup3*F(70)-4*(p2mup3*F(69)-p6mup3*F(71
     &   ))+p4mup3*F(72)
       Fa1(20) = p5mup3*F(80)+p6mup3*F(81)+p4mup3*F(82)
       Fa1(21) = p2mup3*F(83)+p6mup3*F(84)-p5mup3*F(85)+p4mup3*F(86)+p
     &   1mup3*F(87)
       Fa1(22) = -(p2mup3*F(88))+p5mup3*F(89)+p6mup3*F(90)+p4mup3*F(91
     &   )
       Return
       End
