c************** Calling the Fa functions*************************
       subroutine H3jCrossFa2(p1mup3,p1mup4,p2mup3,p2mup4,p3mup3,p3mup
     &   4,p4mup3,p4mup4,p5mup3,p5mup4,p6mup3,p6mup4,mup3mup4,Fa2)
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
       Complex*16 Fa2(23:44)
       COMMON/H3jCrossFaFunctions/Fa
       Fa2(23) = p1mup3*F(92)-p2mup3*F(93)+p4mup3*F(94)-p6mup3*F(95)+p
     &   5mup3*F(96)
       Fa2(24) = p1mup3*F(97)-p2mup3*F(98)-p5mup3*F(99)-p6mup3*F(100)
       Fa2(25) = p2mup3*F(101)+p5mup3*F(102)-p6mup3*F(103)+p4mup3*F(10
     &   4)+p1mup3*F(105)
       Fa2(26) = p1mup3*F(106)-p2mup3*F(107)+p4mup3*F(108)+p5mup3*F(10
     &   9)+p6mup3*F(110)
       Fa2(27) = p1mup3*F(111)-p2mup3*F(112)-p5mup3*F(113)
       Fa2(28) = p2mup3*F(114)+p1mup3*F(115)+p4mup3*F(116)-p5mup3*F(11
     &   7)-p6mup3*F(118)
       Fa2(29) = p2mup3*F(133)-p5mup3*F(134)-p4mup3*F(135)-p6mup3*F(13
     &   6)+p1mup3*F(137)
       Fa2(30) = 4*(p2mup3*F(138)-p5mup3*F(139)-p4mup3*F(140)+p1mup3*F
     &   (141))+p6mup3*F(142)
       Fa2(31) = p2mup3*F(109)-p5mup3*F(117)-p6mup3*F(118)+p1mup3*F(14
     &   3)+p4mup3*F(144)
       Fa2(32) = p6mup3*F(156)+4*(p1mup3*F(154)-p5mup3*F(155)-p4mup3*F
     &   (157)+p2mup3*F(158))
       Fa2(33) = p2mup3*F(93)-p5mup3*F(96)+p1mup3*F(159)-p6mup3*F(160)
     &   -p4mup3*F(161)
       Fa2(34) = p5mup3*F(104)+p2mup3*F(161)+p1mup3*F(165)+p4mup3*F(16
     &   6)+p6mup3*F(167)
       Fa2(35) = p1mup3*F(171)+p2mup3*F(172)-p5mup3*F(173)-p4mup3*F(17
     &   4)-p6mup3*F(218)
       Fa2(36) = -(p2mup3*F(175))-p1mup3*F(176)+p5mup3*F(177)+p6mup3*F
     &   (178)+p4mup3*F(179)
       Fa2(37) = p5mup3*F(105)-p2mup3*F(159)+p4mup3*F(177)-p1mup3*F(18
     &   0)+p6mup3*F(181)
       Fa2(38) = p2mup3*F(161)+p1mup3*F(182)+p5mup3*F(183)+p6mup3*F(18
     &   4)+p4mup3*F(185)
       Fa2(39) = p2mup3*F(96)+p5mup3*F(102)+p4mup3*F(183)+p1mup3*F(186
     &   )+p6mup3*F(187)
       Fa2(40) = p4mup3*F(191)+2*(p1mup3*F(188)+p5mup3*F(189)+p6mup3*F
     &   (190)+p2mup3*F(192))
       Fa2(41) = p5mup3*F(194)+2*(p2mup3*F(160)+p4mup3*F(189)+p1mup3*F
     &   (193)+p6mup3*F(195))
       Fa2(42) = p1mup3*F(196)-p2mup3*F(197)+p5mup3*F(198)+p6mup3*F(19
     &   9)+p4mup3*F(200)
       Fa2(43) = -(p1mup3*F(201))+p5mup3*F(202)+p6mup3*F(203)+p4mup3*F
     &   (204)-p2mup3*F(205)
       Fa2(44) = -(p1mup3*F(206))+p5mup3*F(207)+p6mup3*F(208)+p4mup3*F
     &   (209)-p2mup3*F(210)
       Return
       End
