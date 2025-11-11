c************** Calling the Fa functions*************************
       subroutine H3jCrossIFFa2(p1mup3,p1mup4,p2mup3,p2mup4,p3mup3,p3m
     &   up4,p4mup3,p4mup4,p5mup3,p5mup4,p6mup3,p6mup4,mup3mup4,Fa2)
       IMPLICIT NONE
      Complex*16   p1mup2, p1mup3, p1mup4, p1mup6, p2mup2, p2mup3, 
     -          p2mup4, p2mup6, p3mup2, p3mup3, p3mup4, p3mup6, 
     -          p4mup2, p4mup3, p4mup4, p4mup6, p5mup2, p5mup3, 
     -          p5mup4, p5mup6, p6mup2, p6mup3, p6mup4, p6mup6
       Complex*16   mup2mup3, mup2mup4, mup2mup6, mup3mup4, mup3mup6, 
     -          mup4mup6
        common/H3jCrossIFFhlFunctions/F
       COMMON/H3jCrossIFInvariants/p1sq,p2sq,p3sq,p4sq,p5sq,p6sq,s12,s
     &   23,s34,s45,s56,s16,s123,s234,s345
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
      Complex*16 Fa(40),F(203)
      Real*8 P(61) 
        COMMON/H3jCrossIFPFunctions/P
       Complex*16 Fa2(21:40)
       COMMON/H3jCrossIFFaFunctions/Fa
       Fa2(21) = p5mup3*F(85)+p6mup3*F(86)+p4mup3*F(87)-p2mup3*F(88)+p
     &   1mup3*F(89)
       Fa2(22) = p4mup3*F(85)+p5mup3*F(90)+p6mup3*F(91)-p2mup3*F(92)+p
     &   1mup3*F(93)
       Fa2(23) = -(p5mup3*F(94))-p4mup3*F(95)+p6mup3*F(96)+p2mup3*F(97
     &   )+p1mup3*F(98)
       Fa2(24) = p5mup3*F(99)+p4mup3*F(100)-p2mup3*F(101)+p1mup3*F(102
     &   )+p6mup3*F(103)
       Fa2(25) = p4mup3*F(85)+p5mup3*F(90)-p2mup3*F(92)+p1mup3*F(93)-p
     &   6mup3*F(115)
       Fa2(26) = p5mup3*F(116)+p2mup3*F(117)-p4mup3*F(118)+p1mup3*F(11
     &   9)+p6mup3*F(120)
       Fa2(27) = p4mup3*F(121)-p5mup3*F(122)-p6mup3*F(123)+p1mup3*F(12
     &   4)+p2mup3*F(125)
       Fa2(28) = p5mup3*F(94)-p2mup3*F(97)-p1mup3*F(98)+p6mup3*F(141)+
     &   p4mup3*F(142)
       Fa2(29) = p5mup3*F(144)+p6mup3*F(145)-p2mup3*F(146)+p4mup3*F(14
     &   7)-p1mup3*F(148)
       Fa2(30) = p2mup3*F(142)+p5mup3*F(153)+p6mup3*F(154)+p4mup3*F(15
     &   5)+p1mup3*F(156)
       Fa2(31) = p5mup3*F(161)+p6mup3*F(162)-p1mup3*F(163)-p2mup3*F(16
     &   4)+p4mup3*F(165)
       Fa2(32) = p4mup3*F(88)+p5mup3*F(166)+p6mup3*F(167)+p1mup3*F(168
     &   )+p2mup3*F(169)
       Fa2(33) = p4mup3*F(156)+p5mup3*F(170)+p6mup3*F(171)-p2mup3*F(17
     &   2)-p1mup3*F(173)
       Fa2(34) = p1mup3*F(92)-p5mup3*F(174)-p6mup3*F(175)-p4mup3*F(176
     &   )+p2mup3*F(177)
       Fa2(35) = p4mup3*F(153)-p5mup3*F(178)-p6mup3*F(179)+p1mup3*F(18
     &   0)+p2mup3*F(181)
       Fa2(36) = p5mup3*F(182)+p6mup3*F(183)+p4mup3*F(184)-p1mup3*F(18
     &   5)-p2mup3*F(186)
       Fa2(37) = -(p5mup3*F(122))+p1mup3*F(124)+p2mup3*F(125)+p4mup3*F
     &   (154)+p6mup3*F(187)
       Fa2(38) = p5mup3*F(191)+p6mup3*F(192)-p2mup3*F(193)-p1mup3*F(19
     &   4)+p4mup3*F(195)
       Fa2(39) = p2mup3*F(142)+p5mup3*F(153)+p6mup3*F(154)+p1mup3*F(15
     &   6)-p4mup3*F(196)
       Fa2(40) = p5mup3*F(121)+p6mup3*F(197)-p4mup3*F(198)+p1mup3*F(19
     &   9)+p2mup3*F(200)
       Return
       End
