c************** Calling the Fa functions*************************
       subroutine H3jFa2(p1mup3,p1mup4,p2mup3,p2mup4,p3mup3,p3mup4,p4m
     &   up3,p4mup4,p5mup3,p5mup4,p6mup3,p6mup4,mup3mup4,Fa2)
       IMPLICIT NONE
      Complex*16   p1mup2, p1mup3, p1mup4, p1mup6, p2mup2, p2mup3, 
     -          p2mup4, p2mup6, p3mup2, p3mup3, p3mup4, p3mup6, 
     -          p4mup2, p4mup3, p4mup4, p4mup6, p5mup2, p5mup3, 
     -          p5mup4, p5mup6, p6mup2, p6mup3, p6mup4, p6mup6
       Complex*16   mup2mup3, mup2mup4, mup2mup6, mup3mup4, mup3mup6, 
     -          mup4mup6
        common/H3jFhlFunctions/F
       COMMON/H3jInvariants/p1sq,p2sq,p3sq,p4sq,p5sq,p6sq,s12,s23,s34,
     &   s45,s56,s16,s123,s234,s345
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
      Complex*16 Fa(40),F(206)
      Real*8 P(169) 
        COMMON/H3jPFunctions/P
       Complex*16 Fa2(21:40)
       COMMON/H3jFaFunctions/Fa
       Fa2(21) = -(p2mup3*F(89))+p4mup3*F(90)+p5mup3*F(91)-p6mup3*F(92
     &   )+p1mup3*F(93)
       Fa2(22) = p1mup3*F(94)+p2mup3*F(95)+p4mup3*F(96)+p5mup3*F(97)+p
     &   6mup3*F(98)
       Fa2(23) = p2mup3*F(99)+p1mup3*F(100)
       Fa2(24) = p2mup3*F(85)+p1mup3*F(104)+p4mup3*F(105)+p5mup3*F(106
     &   )-p6mup3*F(107)
       Fa2(25) = p6mup3*F(114)+4*(p1mup3*F(111)-p4mup3*F(112)-p5mup3*F
     &   (113)+p2mup3*F(115))
       Fa2(26) = p6mup3*F(85)-p1mup3*F(116)-p2mup3*F(117)+p4mup3*F(118
     &   )+p5mup3*F(119)
       Fa2(27) = p4mup3*F(96)+p5mup3*F(97)+p6mup3*F(98)+p2mup3*F(118)+
     &   p1mup3*F(120)
       Fa2(28) = p4mup3*F(137)+p5mup3*F(138)-p6mup3*F(139)-p2mup3*F(14
     &   0)+p1mup3*F(141)
       Fa2(29) = p2mup3*F(142)-p4mup3*F(143)-p5mup3*F(144)+p6mup3*F(14
     &   5)+p1mup3*F(146)
       Fa2(30) = p2mup3*F(85)+p4mup3*F(105)-p6mup3*F(107)+p1mup3*F(147
     &   )-p5mup3*F(148)
       Fa2(31) = -(p1mup3*F(161))-p2mup3*F(162)+p4mup3*F(163)+p5mup3*F
     &   (164)+p6mup3*F(206)
       Fa2(32) = -(p2mup3*F(165))-p1mup3*F(166)+p4mup3*F(167)+p5mup3*F
     &   (168)+p6mup3*F(169)
       Fa2(33) = p6mup3*F(105)+p2mup3*F(118)+p1mup3*F(170)+p4mup3*F(17
     &   1)+p5mup3*F(172)
       Fa2(34) = p4mup3*F(174)+2*(p2mup3*F(84)+p1mup3*F(173)+p5mup3*F(
     &   175)+p6mup3*F(176))
       Fa2(35) = -(p2mup3*F(116))+p4mup3*F(168)-p1mup3*F(177)+p5mup3*F
     &   (178)+p6mup3*F(179)
       Fa2(36) = p6mup3*F(106)+p2mup3*F(119)+p4mup3*F(172)+p1mup3*F(18
     &   0)+p5mup3*F(181)
       Fa2(37) = p5mup3*F(183)+2*(p2mup3*F(85)+p4mup3*F(175)+p1mup3*F(
     &   182)+p6mup3*F(184))
       Fa2(38) = p5mup3*F(140)+p1mup3*F(185)-p2mup3*F(186)+p4mup3*F(18
     &   7)+p6mup3*F(188)
       Fa2(39) = -(p1mup3*F(189))+p4mup3*F(190)+p5mup3*F(191)+p6mup3*F
     &   (192)-p2mup3*F(193)
       Fa2(40) = -(p1mup3*F(194))+p4mup3*F(195)+p5mup3*F(196)+p6mup3*F
     &   (197)-p2mup3*F(198)
       Return
       End
