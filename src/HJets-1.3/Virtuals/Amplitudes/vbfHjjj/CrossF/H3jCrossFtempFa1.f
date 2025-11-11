c************** Calling the Fa functions*************************
       subroutine H3jCrossFFa1(p1mup3,p1mup4,p2mup3,p2mup4,p3mup3,p3mu
     &   p4,p4mup3,p4mup4,p5mup3,p5mup4,p6mup3,p6mup4,mup3mup4,Fa1)
       IMPLICIT NONE
      Complex*16   p1mup2, p1mup3, p1mup4, p1mup6, p2mup2, p2mup3, 
     -          p2mup4, p2mup6, p3mup2, p3mup3, p3mup4, p3mup6, 
     -          p4mup2, p4mup3, p4mup4, p4mup6, p5mup2, p5mup3, 
     -          p5mup4, p5mup6, p6mup2, p6mup3, p6mup4, p6mup6
       Complex*16   mup2mup3, mup2mup4, mup2mup6, mup3mup4, mup3mup6, 
     -          mup4mup6
        common/H3jCrossFFhlFunctions/F
       COMMON/H3jCrossFInvariants/p1sq,p2sq,p3sq,p4sq,p5sq,p6sq,s12,s2
     &   3,s34,s45,s56,s16,s123,s234,s345
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
      Complex*16 Fa(44),F(214)
      Real*8 P(175) 
        COMMON/H3jCrossFPFunctions/P
       Complex*16 Fa1(1:22)
       COMMON/H3jCrossFFaFunctions/Fa
       Fa1(1) = p5mup3*F(30)+p6mup3*F(31)+p1mup3*F(32)+p4mup3*F(33)
       Fa1(2) = p5mup3*F(36)+p6mup3*F(37)
       Fa1(3) = p1mup3*F(38)+p5mup3*F(39)+p6mup3*F(40)
       Fa1(4) = -2*p1mup3*F(41)+p5mup3*F(211)
       Fa1(5) = p5mup3*F(42)+p6mup3*F(43)-p4mup3*F(44)
       Fa1(6) = p4mup3*F(39)-p1mup3*F(42)-p6mup3*F(45)
       Fa1(7) = -2*p1mup3*F(46)+p4mup3*F(47)+p5mup3*F(211)
       Fa1(8) = p5mup3*F(39)-p1mup3*F(55)+p6mup3*F(56)
       Fa1(9) = p5mup3*F(42)+p4mup3*F(59)+p6mup3*F(60)
       Fa1(10) = p4mup3*F(61)-p5mup3*F(62)
       Fa1(11) = -(p1mup3*F(60))+p4mup3*F(63)+p5mup3*F(64)
       Fa1(12) = p1mup3+p5mup3+p6mup3
       Fa1(13) = p5mup3*F(42)+p6mup3*F(65)
       Fa1(14) = p1mup3*F(42)+p6mup3*F(64)
       Fa1(15) = p5mup3*F(64)-p1mup3*F(65)
       Fa1(16) = p5mup3*F(59)+p6mup3*F(61)
       Fa1(17) = p1mup3*F(59)+p6mup3*F(63)
       Fa1(18) = p1mup3*F(61)-p5mup3*F(63)
       Fa1(19) = p1mup3*F(75)+4*(p5mup3*F(73)-p6mup3*F(74)+p2mup3*F(76
     &   )-p4mup3*F(77))
       Fa1(20) = p4mup3*F(81)-p5mup3*F(82)-p6mup3*F(83)
       Fa1(21) = p1mup3*F(84)+p4mup3*F(85)+p5mup3*F(86)+p6mup3*F(87)
       Fa1(22) = p1mup3*F(88)+p4mup3*F(89)-p5mup3*F(90)+p2mup3*F(91)-p
     &   6mup3*F(92)
       Return
       End
