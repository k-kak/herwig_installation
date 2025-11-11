c************** Calling the Fa functions*************************
       subroutine NoAbeH3jFa2(p1mup2,p2mup2,p3mup2,p4mup2,p5mup2,p6mup
     &   2,Fa2)
       IMPLICIT NONE
      Complex*16   p1mup2, p1mup3, p1mup4, p1mup6, p2mup2, p2mup3, 
     -          p2mup4, p2mup6, p3mup2, p3mup3, p3mup4, p3mup6, 
     -          p4mup2, p4mup3, p4mup4, p4mup6, p5mup2, p5mup3, 
     -          p5mup4, p5mup6, p6mup2, p6mup3, p6mup4, p6mup6
       Complex*16   mup2mup3, mup2mup4, mup2mup6, mup3mup4, mup3mup6, 
     -          mup4mup6
        common/NoAbeH3jFhlFunctions/F
       COMMON/NoAbeH3jInvariants/p1sq,p2sq,p3sq,p4sq,p5sq,p6sq,s12,s23
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
      Complex*16 Fa(22),F(172)
      Real*8 P(116) 
        COMMON/NoAbeH3jPFunctions/P
       Complex*16 Fa2(12:22)
       COMMON/NoAbeH3jFaFunctions/Fa
       Fa2(12) = p3mup2*F(79)+p4mup2*F(80)+p5mup2*F(81)+p6mup2*F(82)+p
     &   1mup2*F(83)
       Fa2(13) = p5mup2*F(46)+p6mup2*F(47)+p3mup2*F(84)+p4mup2*F(85)-p
     &   1mup2*F(86)
       Fa2(14) = p3mup2*F(87)+p4mup2*F(88)+p5mup2*F(89)+p6mup2*F(90)+p
     &   1mup2*F(91)
       Fa2(15) = p3mup2*F(92)+p4mup2*F(93)+p5mup2*F(94)+p6mup2*F(95)+p
     &   1mup2*F(96)
       Fa2(16) = p5mup2*F(99)-2*(p3mup2*F(97)+p4mup2*F(98)+p6mup2*F(10
     &   0)+p1mup2*F(101))
       Fa2(17) = p3mup2*F(102)+p4mup2*F(103)+p5mup2*F(104)+p6mup2*F(10
     &   5)+p1mup2*F(106)
       Fa2(18) = p3mup2*F(107)+p4mup2*F(108)+p5mup2*F(109)+p6mup2*F(11
     &   0)+p1mup2*F(111)
       Fa2(19) = p3mup2*F(112)+p4mup2*F(113)+p5mup2*F(114)+p6mup2*F(11
     &   5)+p1mup2*F(116)
       Fa2(20) = p3mup2*F(117)+p4mup2*F(118)+p5mup2*F(119)+p6mup2*F(12
     &   0)+p1mup2*F(121)
       Fa2(21) = p3mup2*F(122)+p4mup2*F(123)+p5mup2*F(124)+p6mup2*F(12
     &   5)+p1mup2*F(126)
       Fa2(22) = p3mup2*F(147)+p4mup2*F(148)+p5mup2*F(149)+p6mup2*F(15
     &   0)+p1mup2*F(151)
       Return
       End
