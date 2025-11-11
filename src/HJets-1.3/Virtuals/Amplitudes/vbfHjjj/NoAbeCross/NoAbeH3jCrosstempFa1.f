c************** Calling the Fa functions*************************
       subroutine NoAbeH3jCrossFa1(p1mup2,p2mup2,p3mup2,p4mup2,p5mup2,
     &   p6mup2,Fa1)
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
       Complex*16 Fa1(1:13)
       COMMON/NoAbeH3jCrossFaFunctions/Fa
       Fa1(1) = p1mup2*F(1)+p4mup2*F(2)+p5mup2*F(3)+2*p6mup2*F(4)+p3mu
     &   p2*F(5)
       Fa1(2) = p1mup2*F(9)+p4mup2*F(10)+p5mup2*F(11)+p6mup2*F(12)+p3m
     &   up2*F(13)
       Fa1(3) = p4mup2*F(14)+p5mup2*F(15)+p6mup2*F(16)+p3mup2*F(17)
       Fa1(4) = p1mup2*F(18)-p3mup2*F(19)-p5mup2*F(20)-p6mup2*F(21)
       Fa1(5) = p1mup2*F(22)+p4mup2*F(23)+p5mup2*F(24)+p6mup2*F(25)+p3
     &   mup2*F(26)
       Fa1(6) = p1mup2*F(27)-p3mup2*F(28)-p5mup2*F(29)
       Fa1(7) = p1mup2*F(30)+p4mup2*F(31)+p5mup2*F(32)+p6mup2*F(33)+p3
     &   mup2*F(34)
       Fa1(8) = p1mup2*F(35)-p4mup2*F(36)-p5mup2*F(37)-p6mup2*F(38)-p3
     &   mup2*F(39)
       Fa1(9) = p4mup2*F(41)+p5mup2*F(42)+p6mup2*F(43)+4*(p1mup2*F(40)
     &   +p3mup2*F(44))
       Fa1(10) = p4mup2*F(45)+p5mup2*F(46)+p6mup2*F(47)
       Fa1(11) = p1mup2*F(52)+p4mup2*F(53)+p5mup2*F(54)+p6mup2*F(55)+p
     &   3mup2*F(56)
       Fa1(12) = p1mup2*F(57)+p4mup2*F(58)+p5mup2*F(59)+p6mup2*F(60)+p
     &   3mup2*F(61)
       Fa1(13) = p5mup2*F(32)+p6mup2*F(33)-p1mup2*F(70)+p4mup2*F(71)+p
     &   3mup2*F(72)
       Return
       End
