c************** Calling the Fa functions*************************
       subroutine H3jCrossFFa2(p1mup3,p1mup4,p2mup3,p2mup4,p3mup3,p3mu
     &   p4,p4mup3,p4mup4,p5mup3,p5mup4,p6mup3,p6mup4,mup3mup4,Fa2)
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
       Complex*16 Fa2(23:44)
       COMMON/H3jCrossFFaFunctions/Fa
       Fa2(23) = -(p5mup3*F(93))-p6mup3*F(94)+p2mup3*F(95)+p4mup3*F(96
     &   )+p1mup3*F(97)
       Fa2(24) = p2mup3*F(98)+p1mup3*F(99)+p4mup3*F(100)+p5mup3*F(101)
     &   +p6mup3*F(102)
       Fa2(25) = p2mup3*F(103)+p1mup3*F(104)+p5mup3*F(105)
       Fa2(26) = -(p5mup3*F(106))-p6mup3*F(107)+p4mup3*F(108)+p2mup3*F
     &   (109)+p1mup3*F(110)
       Fa2(27) = p2mup3*F(111)+p1mup3*F(112)+p4mup3*F(113)+p5mup3*F(11
     &   4)+p6mup3*F(115)
       Fa2(28) = p2mup3*F(116)+p1mup3*F(117)-p4mup3*F(118)+p5mup3*F(11
     &   9)
       Fa2(29) = p5mup3*F(132)+p6mup3*F(133)+p4mup3*F(134)-p2mup3*F(13
     &   5)-p1mup3*F(136)
       Fa2(30) = p5mup3*F(111)+p6mup3*F(137)-p2mup3*F(138)-p1mup3*F(13
     &   9)+p4mup3*F(140)
       Fa2(31) = p4mup3*F(100)+p5mup3*F(101)+p2mup3*F(140)+p6mup3*F(14
     &   1)+p1mup3*F(142)
       Fa2(32) = p4mup3*F(93)+p5mup3*F(106)-p1mup3*F(110)-p2mup3*F(150
     &   )+p6mup3*F(151)
       Fa2(33) = p6mup3*F(155)+4*(p2mup3*F(152)-p4mup3*F(153)+p5mup3*F
     &   (154)+p1mup3*F(156))
       Fa2(34) = p2mup3*F(157)+p4mup3*F(158)-p5mup3*F(159)-p6mup3*F(16
     &   0)+p1mup3*F(161)
       Fa2(35) = p5mup3*F(169)+p6mup3*F(170)-p1mup3*F(171)-p2mup3*F(17
     &   2)+p4mup3*F(173)
       Fa2(36) = p5mup3*F(174)+p6mup3*F(175)+p2mup3*F(176)+p1mup3*F(17
     &   7)+p4mup3*F(178)
       Fa2(37) = p4mup3*F(179)+p5mup3*F(180)+p6mup3*F(181)-p1mup3*F(18
     &   2)-p2mup3*F(183)
       Fa2(38) = p4mup3*F(184)+p5mup3*F(185)+p6mup3*F(186)-p1mup3*F(18
     &   7)-p2mup3*F(188)
       Fa2(39) = p5mup3*F(112)+p4mup3*F(142)+p6mup3*F(192)-p2mup3*F(19
     &   3)-p1mup3*F(194)
       Fa2(40) = p4mup3*F(101)-p5mup3*F(195)-p6mup3*F(196)+p1mup3*F(19
     &   7)+p2mup3*F(198)
       Fa2(41) = p4mup3*F(141)+p1mup3*F(161)-p5mup3*F(199)+p6mup3*F(20
     &   0)+p2mup3*F(201)
       Fa2(42) = p5mup3*F(99)+p6mup3*F(202)-p2mup3*F(203)-p1mup3*F(204
     &   )+p4mup3*F(205)
       Fa2(43) = p5mup3*F(101)+p2mup3*F(140)+p6mup3*F(141)+p1mup3*F(14
     &   2)-p4mup3*F(206)
       Fa2(44) = p5mup3*F(102)-p4mup3*F(207)+p6mup3*F(208)+p1mup3*F(20
     &   9)+p2mup3*F(210)
       Return
       End
