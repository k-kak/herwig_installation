c************** Calling the Fa functions*************************
       subroutine NoAbeH3jCrossFFa2(p1mup2,p2mup2,p3mup2,p4mup2,p5mup2
     &   ,p6mup2,Fa2)
       IMPLICIT NONE
      Complex*16   p1mup2, p1mup3, p1mup4, p1mup6, p2mup2, p2mup3, 
     -          p2mup4, p2mup6, p3mup2, p3mup3, p3mup4, p3mup6, 
     -          p4mup2, p4mup3, p4mup4, p4mup6, p5mup2, p5mup3, 
     -          p5mup4, p5mup6, p6mup2, p6mup3, p6mup4, p6mup6
       Complex*16   mup2mup3, mup2mup4, mup2mup6, mup3mup4, mup3mup6, 
     -          mup4mup6
        common/NoAbeH3jCrossFFhlFunctions/F
       COMMON/NoAbeH3jCrossFInvariants/p1sq,p2sq,p3sq,p4sq,p5sq,p6sq,s
     &   12,s23,s34,s45,s56,s16,s123,s234,s345
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
      Complex*16 Fa(26),F(187)
      Real*8 P(111) 
        COMMON/NoAbeH3jCrossFPFunctions/P
       Complex*16 Fa2(14:26)
       COMMON/NoAbeH3jCrossFFaFunctions/Fa
       Fa2(14) = p3mup2*F(89)+p4mup2*F(90)+p5mup2*F(91)+2*(p6mup2*F(92
     &   )+p1mup2*F(93))
       Fa2(15) = p3mup2*F(94)+p4mup2*F(95)+p5mup2*F(96)+p6mup2*F(97)-p
     &   1mup2*F(98)
       Fa2(16) = p3mup2*F(99)+p4mup2*F(100)+p5mup2*F(101)+p6mup2*F(102
     &   )+p1mup2*F(103)
       Fa2(17) = p3mup2*F(104)+p4mup2*F(105)+p5mup2*F(106)+p6mup2*F(10
     &   7)+p1mup2*F(108)
       Fa2(18) = p3mup2*F(109)+p4mup2*F(110)+p5mup2*F(111)+p6mup2*F(11
     &   2)+p1mup2*F(113)
       Fa2(19) = p3mup2*F(114)+p4mup2*F(115)+p5mup2*F(116)+p6mup2*F(11
     &   7)+p1mup2*F(118)
       Fa2(20) = p5mup2*F(121)-2*(p1mup2*F(38)+p3mup2*F(119)+p4mup2*F(
     &   120)+p6mup2*F(122))
       Fa2(21) = p3mup2*F(123)+p4mup2*F(124)+p5mup2*F(125)+p6mup2*F(12
     &   6)+p1mup2*F(127)
       Fa2(22) = p3mup2*F(128)+p4mup2*F(129)+p5mup2*F(130)+p6mup2*F(13
     &   1)+p1mup2*F(132)
       Fa2(23) = p3mup2*F(133)+p4mup2*F(134)+p5mup2*F(135)+p6mup2*F(13
     &   6)+p1mup2*F(137)
       Fa2(24) = p3mup2*F(138)+p4mup2*F(139)+p5mup2*F(140)+p6mup2*F(14
     &   1)+p1mup2*F(142)
       Fa2(25) = p3mup2*F(143)+p4mup2*F(144)+p5mup2*F(145)+p6mup2*F(14
     &   6)+p1mup2*F(147)
       Fa2(26) = p3mup2*F(160)+p4mup2*F(161)+p5mup2*F(162)+p6mup2*F(16
     &   3)+p1mup2*F(164)
       Return
       End
