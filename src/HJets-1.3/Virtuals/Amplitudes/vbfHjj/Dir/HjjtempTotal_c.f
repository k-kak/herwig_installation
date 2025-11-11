
       subroutine Hjj66T_c(M1SQ,p1,p2,p3,p4,p5,barpsi_p1,psi_p5,barpsi_p2
     &   ,psi_p3,musq,comp,result_F,resultborn,Div)
       implicit none
       complex*16 M1SQ
       real*8 p1(0:3),p2(0:3),p3(0:3),p4(0:3),p5(0:3)
       Complex*16 barpsi_p1(4),psi_p5(4),barpsi_p2(4),psi_p3(4)
       real*8 musq
       integer comp, Div
       complex*16 result_F, resultborn
       complex*16 result(3)

       if (Div.eq.0) then 

       call  Hjj66(M1SQ,p1,p2,p3,p4,p5,barpsi_p1,psi_p5,barpsi_p2
     &   ,psi_p3,musq,comp,0,0,result,resultborn)

       else
       call  Hjj66Div(M1SQ,p1,p2,p3,p4,p5,barpsi_p1,psi_p5,barpsi_p2
     &   ,psi_p3,musq,comp,0,0,result,resultborn,Div)

       endif

       result_F=result(1)


       end




       subroutine Hjj67T_c(M1SQ,p1,p2,p3,p4,p5,barpsi_p1,psi_p5,barpsi_p2
     &   ,psi_p3,musq,comp,result_F,resultborn,Div)
       implicit none
       complex*16 M1SQ
       real*8 p1(0:3),p2(0:3),p3(0:3),p4(0:3),p5(0:3)
       Complex*16 barpsi_p1(4),psi_p5(4),barpsi_p2(4),psi_p3(4)
       real*8 musq
       integer comp, Div
       complex*16 result_F, resultborn
       complex*16 result(3)

       if (Div.eq.0) then 

       call  Hjj67(M1SQ,p1,p2,p3,p4,p5,barpsi_p1,psi_p5,barpsi_p2
     &   ,psi_p3,musq,comp,0,0,result,resultborn)

       else
       call  Hjj67Div(M1SQ,p1,p2,p3,p4,p5,barpsi_p1,psi_p5,barpsi_p2
     &   ,psi_p3,musq,comp,0,0,result,resultborn,Div)

       endif

       result_F=result(1)


       end



       subroutine Hjj76T_c(M1SQ,p1,p2,p3,p4,p5,barpsi_p1,psi_p5,barpsi_p2
     &   ,psi_p3,musq,comp,result_F,resultborn,Div)
       implicit none
       complex*16 M1SQ
       real*8 p1(0:3),p2(0:3),p3(0:3),p4(0:3),p5(0:3)
       Complex*16 barpsi_p1(4),psi_p5(4),barpsi_p2(4),psi_p3(4)
       real*8 musq
       integer comp, Div
       complex*16 result_F, resultborn
       complex*16 result(3)

       if (Div.eq.0) then 

       call  Hjj76(M1SQ,p1,p2,p3,p4,p5,barpsi_p1,psi_p5,barpsi_p2
     &   ,psi_p3,musq,comp,0,0,result,resultborn)

       else
       call  Hjj76Div(M1SQ,p1,p2,p3,p4,p5,barpsi_p1,psi_p5,barpsi_p2
     &   ,psi_p3,musq,comp,0,0,result,resultborn,Div)

       endif

       result_F=result(1)


       end


       subroutine Hjj77T_c(M1SQ,p1,p2,p3,p4,p5,barpsi_p1,psi_p5,barpsi_p2
     &   ,psi_p3,musq,comp,result_F,resultborn,Div)
       implicit none
       complex*16 M1SQ
       real*8 p1(0:3),p2(0:3),p3(0:3),p4(0:3),p5(0:3)
       Complex*16 barpsi_p1(4),psi_p5(4),barpsi_p2(4),psi_p3(4)
       real*8 musq
       integer comp, Div
       complex*16 result_F, resultborn
       complex*16 result(3)

       if (Div.eq.0) then 

       call  Hjj77(M1SQ,p1,p2,p3,p4,p5,barpsi_p1,psi_p5,barpsi_p2
     &   ,psi_p3,musq,comp,0,0,result,resultborn)

       else
       call  Hjj77Div(M1SQ,p1,p2,p3,p4,p5,barpsi_p1,psi_p5,barpsi_p2
     &   ,psi_p3,musq,comp,0,0,result,resultborn,Div)

       endif

       result_F=result(1)


       end










