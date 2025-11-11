c************** Calling the Fa functions*************************
       subroutine NoAbeH3jCrossFa2(p1mup2,p2mup2,p3mup2,p4mup2,p5mup2,
     &   p6mup2,Fa2)
       IMPLICIT NONE
      Complex*16   p1mup2, p1mup3, p1mup4, p1mup6, p2mup2, p2mup3, 
     -          p2mup4, p2mup6, p3mup2, p3mup3, p3mup4, p3mup6, 
     -          p4mup2, p4mup3, p4mup4, p4mup6, p5mup2, p5mup3, 
     -          p5mup4, p5mup6, p6mup2, p6mup3, p6mup4, p6mup6
       Complex*16   mup2mup3, mup2mup4, mup2mup6, mup3mup4, mup3mup6, 
     -          mup4mup6
        common/NoAbeH3jCrossFhlFunctions/F
       COMMON/NoAbeH3jCrossInvariants/p1sq,p2sq,p3sq,p4sq,p5sq,p6sq,s1
     &   2,s23,s34,s45,s56,s16,s123,s234,s345
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
      Complex*16 Fa(26),F(186)
      Real*8 P(91) 
        COMMON/NoAbeH3jCrossPFunctions/P
       Complex*16 Fa2(14:26)
       COMMON/NoAbeH3jCrossFaFunctions/Fa
       Fa2(14) = p4mup2*F(77)+2*(p1mup2*F(76)-p5mup2*F(78)-p6mup2*F(79
     &   )-p3mup2*F(80))
       Fa2(15) = p1mup2*F(81)+p4mup2*F(82)+p5mup2*F(83)+p6mup2*F(84)+p
     &   3mup2*F(85)
       Fa2(16) = p1mup2*F(86)+p4mup2*F(87)+p5mup2*F(88)+p6mup2*F(89)+p
     &   3mup2*F(90)
       Fa2(17) = p1mup2*F(91)+p4mup2*F(92)+p5mup2*F(93)+p6mup2*F(94)+p
     &   3mup2*F(95)
       Fa2(18) = p1mup2*F(96)+p4mup2*F(97)+p5mup2*F(98)+p6mup2*F(99)+p
     &   3mup2*F(100)
       Fa2(19) = p1mup2*F(101)+p4mup2*F(102)+p5mup2*F(103)+p6mup2*F(10
     &   4)+p3mup2*F(105)
       Fa2(20) = p1mup2*F(106)+p4mup2*F(107)+p5mup2*F(108)+p6mup2*F(10
     &   9)+p3mup2*F(110)
       Fa2(21) = p1mup2*F(111)+p4mup2*F(112)+p5mup2*F(113)+p6mup2*F(11
     &   4)+p3mup2*F(115)
       Fa2(22) = p1mup2*F(116)+p4mup2*F(117)+p5mup2*F(118)+p6mup2*F(11
     &   9)+p3mup2*F(120)
       Fa2(23) = p5mup2*F(37)-p1mup2*F(121)+p4mup2*F(122)+p6mup2*F(123
     &   )-p3mup2*F(124)
       Fa2(24) = p1mup2*F(125)+p4mup2*F(126)+p5mup2*F(127)+p6mup2*F(12
     &   8)+p3mup2*F(129)
       Fa2(25) = p1mup2*F(130)+p4mup2*F(131)+p5mup2*F(132)+p6mup2*F(13
     &   3)+p3mup2*F(134)
       Fa2(26) = p1mup2*F(135)+p4mup2*F(136)+p5mup2*F(137)+p6mup2*F(13
     &   8)+p3mup2*F(139)
       Return
       End
