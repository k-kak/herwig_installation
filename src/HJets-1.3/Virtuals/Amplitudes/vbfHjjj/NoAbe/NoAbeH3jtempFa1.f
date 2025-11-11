c************** Calling the Fa functions*************************
       subroutine NoAbeH3jFa1(p1mup2,p2mup2,p3mup2,p4mup2,p5mup2,p6mup
     &   2,Fa1)
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
       Complex*16 Fa1(1:11)
       COMMON/NoAbeH3jFaFunctions/Fa
       Fa1(1) = p3mup2*F(1)+p5mup2*F(3)-2*(p4mup2*F(2)+p6mup2*F(4))+p1
     &   mup2*F(5)
       Fa1(2) = p3mup2*F(10)+p4mup2*F(11)+p5mup2*F(12)+p6mup2*F(13)+p1
     &   mup2*F(14)
       Fa1(3) = p3mup2*F(15)+p4mup2*F(16)+p5mup2*F(17)+p6mup2*F(18)+p1
     &   mup2*F(19)
       Fa1(4) = p3mup2*F(20)-p4mup2*F(21)-p5mup2*F(22)-p6mup2*F(23)+p1
     &   mup2*F(24)
       Fa1(5) = p3mup2*F(32)+p1mup2*F(33)
       Fa1(6) = p3mup2*F(36)-p4mup2*F(37)+p5mup2*F(38)+p6mup2*F(39)-p1
     &   mup2*F(40)
       Fa1(7) = p4mup2*F(45)-2*(p3mup2*F(44)+p5mup2*F(46)+p6mup2*F(47)
     &   -p1mup2*F(48))
       Fa1(8) = p3mup2*F(61)+p4mup2*F(62)+p5mup2*F(63)+2*(p6mup2*F(64)
     &   +p1mup2*F(65))
       Fa1(9) = p5mup2*F(38)+p6mup2*F(39)+p3mup2*F(66)+p4mup2*F(67)-p1
     &   mup2*F(68)
       Fa1(10) = p3mup2*F(69)+p4mup2*F(70)+p5mup2*F(71)+p6mup2*F(72)+p
     &   1mup2*F(73)
       Fa1(11) = p3mup2*F(74)+p4mup2*F(75)+p5mup2*F(76)+p6mup2*F(77)+p
     &   1mup2*F(78)
       Return
       End
