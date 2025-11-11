       subroutine HjjCrossIF66T_c(M1SQ,p1,p2,p3,p4,p5,psi_p1,barpsi_p5,ps
     &   i_p2,barpsi_p3,musq,comp,result_F,resultborn,Div)
       implicit none
       complex*16 M1SQ
       real*8 p1(0:3),p2(0:3),p3(0:3),p4(0:3),p5(0:3)
       Complex*16 psi_p1(4),barpsi_p5(4),psi_p2(4),barpsi_p3(4)
       real*8 musq
       integer comp, Div
       complex*16 result_F, resultborn
       complex*16 result(3)

       if (Div.eq.0) then 

       call HjjCrossIF66(M1SQ,p1,p2,p3,p4,p5,psi_p1,barpsi_p5,ps
     &   i_p2,barpsi_p3,musq,comp,0,0,result,resultborn)

          else

       call HjjCrossIF66Div(M1SQ,p1,p2,p3,p4,p5,psi_p1,barpsi_p5,ps
     &   i_p2,barpsi_p3,musq,comp,0,0,result,resultborn,Div)

       endif
       result_F=result(1)   
       end




       subroutine HjjCrossIF67T_c(M1SQ,p1,p2,p3,p4,p5,psi_p1,barpsi_p5,ps
     &   i_p2,barpsi_p3,musq,comp,result_F,resultborn,Div)
       implicit none
       complex*16 M1SQ
       real*8 p1(0:3),p2(0:3),p3(0:3),p4(0:3),p5(0:3)
       Complex*16 psi_p1(4),barpsi_p5(4),psi_p2(4),barpsi_p3(4)
       real*8 musq
       integer comp, Div
       complex*16 result_F, resultborn
       complex*16 result(3)

       if (Div.eq.0) then 

       call HjjCrossIF67(M1SQ,p1,p2,p3,p4,p5,psi_p1,barpsi_p5,ps
     &   i_p2,barpsi_p3,musq,comp,0,0,result,resultborn)

          else

       call HjjCrossIF67Div(M1SQ,p1,p2,p3,p4,p5,psi_p1,barpsi_p5,ps
     &   i_p2,barpsi_p3,musq,comp,0,0,result,resultborn,Div)

       endif
       result_F=result(1)   
       end


       subroutine HjjCrossIF76T_c(M1SQ,p1,p2,p3,p4,p5,psi_p1,barpsi_p5,ps
     &   i_p2,barpsi_p3,musq,comp,result_F,resultborn,Div)
       implicit none
       complex*16 M1SQ
       real*8 p1(0:3),p2(0:3),p3(0:3),p4(0:3),p5(0:3)
       Complex*16 psi_p1(4),barpsi_p5(4),psi_p2(4),barpsi_p3(4)
       real*8 musq
       integer comp, Div
       complex*16 result_F, resultborn
       complex*16 result(3)

       if (Div.eq.0) then 

       call HjjCrossIF76(M1SQ,p1,p2,p3,p4,p5,psi_p1,barpsi_p5,ps
     &   i_p2,barpsi_p3,musq,comp,0,0,result,resultborn)

          else

       call HjjCrossIF76Div(M1SQ,p1,p2,p3,p4,p5,psi_p1,barpsi_p5,ps
     &   i_p2,barpsi_p3,musq,comp,0,0,result,resultborn,Div)

       endif
       result_F=result(1)   
       end





       subroutine HjjCrossIF77T_c(M1SQ,p1,p2,p3,p4,p5,psi_p1,barpsi_p5,ps
     &   i_p2,barpsi_p3,musq,comp,result_F,resultborn,Div)
       implicit none
       complex*16 M1SQ
       real*8 p1(0:3),p2(0:3),p3(0:3),p4(0:3),p5(0:3)
       Complex*16 psi_p1(4),barpsi_p5(4),psi_p2(4),barpsi_p3(4)
       real*8 musq
       integer comp, Div
       complex*16 result_F, resultborn
       complex*16 result(3)

       if (Div.eq.0) then 

       call HjjCrossIF77(M1SQ,p1,p2,p3,p4,p5,psi_p1,barpsi_p5,ps
     &   i_p2,barpsi_p3,musq,comp,0,0,result,resultborn)

          else

       call HjjCrossIF77Div(M1SQ,p1,p2,p3,p4,p5,psi_p1,barpsi_p5,ps
     &   i_p2,barpsi_p3,musq,comp,0,0,result,resultborn,Div)

       endif
       result_F=result(1)   
       end

