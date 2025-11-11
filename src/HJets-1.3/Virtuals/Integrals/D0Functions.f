c      function D0fin(s,t,p1sq,p2sq,p3sq,p4sq,musq)
c      function D00m_fin(s,t,musq)
c      function D01m_fin(s,t,m4sq,musq)
c      function D02m_fin(s,t,m3sq,m4sq,musq)
c      function D02m_fine(s,t,m1sq,m3sq,musq)
c      function D03m_fin(s,t,m2sq,m3sq,m4sq,musq)
c      function D03m_fin_BDK(s,t,m2sq,m3sq,m4sq,musq)
c      function D04m_fin(s,t,m1sq,m2sq,m3sq,m4sq,musq)
c      function B0fin(p1sq,musq)
c      function C0fin(p1sq,p2sq,p3sq,musq)
c      function E0(psq, qsq, lsq, tsq,
c      function sgnratio(a,b)
c      function im_part(a,b)
c      function dilog_1_m_aob(a,b)
c      function eta3(im_a,im_b,im_ab)
c      function theta(x)
c      function sgn(x)
c      function delta(x,y)
c




      function D0fin(s,t,p1sq,p2sq,p3sq,p4sq,musq)
      implicit none
      complex * 16 D0fin
      real * 8 s,t,musq,p1sq,p2sq,p3sq,p4sq
      Complex*16 D01m_fin,D02m_fin,D02m_fine,D03m_fin,D04m_fin,D00m_fin
      External D01m_fin,D02m_fin,D02m_fine,D03m_fin,D04m_fin,D00m_fin
      real * 8 tiny
      parameter (tiny=1d-6)
c D0fin(s,t,p1sq,p2sq,p3sq,p4sq,musq)
c Author: Francisco Campanario
c Date:2008-04-16
c Modified: 2009-01-14(Including s=0 or t=0)
c I have used Carlo Oleari D02m,D03m,D04m
c Bug solve in analytical continuation of D03m
c Need  my  D00m,D02me,D01m functions
c
c                     k
c     p1 ->-----------<-------------<-- p4
c             |                 |
c             |                 |
c             |                 |  
c     p2 ->-------------------------<-- p3
c
c
c    s=(p1+p2)^2
c    t=(p2+p3)^2
c    musq = mu^2 reference dimensional scale
c 
c  int d^dk/(2 pi)^d 1/(k^2)/(k+p1)^2/(k+p1+p2)^2/(k+p1+p2+q1)^2= 
c          N_ep * D0fin(s,t,q1^2,q2^2)
c          N_ep = i/(4 pi)^2 (4 pi)^ep Gamma(1+ep) (musq)^(-ep)
c   the external momenta could be whaterver value      
      if(abs(p1sq).lt.tiny) then
          if(abs(p2sq).lt.tiny)then
               if(abs(p3sq).lt.tiny)then
                   if(abs(p4sq).lt.tiny) then
                      if(abs(s).lt.tiny) then
                         if(abs(t).lt.tiny) then   
c  D(0,0,0,0,0,0) 
                             D0fin=0d0
                               RETURN
                         else    
c  D(0,t,0,0,0,0)                         
                            D0fin=0d0
                               RETURN
                         endif  !t if
                      else
                         if(abs(t).lt.tiny) then
c  D(s,0,0,0,0,0)
                            D0fin=0d0
                               RETURN
                         else      
c  D(s,t,0,0,0,0)                                                      
                             D0fin=D00m_fin(s,t,musq)
                             RETURN
                         endif  !t if
                       endif    !s if   
                   else   ! p4sq
                      if(abs(s).lt.tiny) then
                         if(abs(t).lt.tiny) then   
c  D(0,0,0,0,0,p4sq) 
                             D0fin=0d0
                               RETURN
                         else    
c  D(0,t,0,0,0,p4sq)                         
                            D0fin=0d0
                               RETURN
                         endif  !t if
                      else
                         if(abs(t).lt.tiny) then
c  D(s,0,0,0,0,p4sq)
                            D0fin=0d0
                               RETURN
                         else      
c D(s,t,0,0,0,p4sq) 
                            D0fin=D01m_fin(s,t,p4sq,musq)
                              RETURN

                         endif  !t if
                       endif    !s if   
                   endif        !m4sq if
               else             !m3sq if
                   if(abs(p4sq).lt.tiny) then
                      if(abs(s).lt.tiny) then
                         if(abs(t).lt.tiny) then   
c  D(0,0,0,0,p3sq,0) 
                             D0fin=0d0
                               RETURN
                         else    
c  D(0,t,0,0,p3sq,0)                         
                            D0fin=0d0
                               RETURN
                         endif  !t if
                      else
                         if(abs(t).lt.tiny) then
c  D(s,0,0,0,p3sq,0)
                            D0fin=0d0
                               RETURN
                         else      
c  D(s,t,0,0,p3sq,0)                                                      
                             D0fin=D01m_fin(s,t,p3sq,musq) 
                             RETURN
                         endif  !t if
                       endif     
                   else   ! m4sq
                      if(abs(s).lt.tiny) then
                         if(abs(t).lt.tiny) then   
c  D(0,0,0,0,p3sq,p4sq) 
                             D0fin=0d0
                               RETURN
                         else    
c  D(0,t,0,0,p3sq,p4sq)                         
                            D0fin=0d0
                               RETURN
                         endif  !t if
                      else
                         if(abs(t).lt.tiny) then
c  D(s,0,0,0,p3sq,p4sq)
                            D0fin=0d0
                               RETURN
                         else      
c  D(s,t,0,0,p3sq,p4sq)                                                      
                         D0fin=D02m_fin(s,t,p3sq,p4sq,musq)     
                             RETURN
                         endif  !t if
                       endif    !s if   
                   endif  ! p4sq 
               endif      ! p3sq
          else          ! p2sq
                     if(abs(p3sq).lt.tiny)then
                           if(abs(p4sq).lt.tiny) then
                              if(abs(s).lt.tiny) then
                                 if(abs(t).lt.tiny) then   
c  D(0,0,0,p2sq,0,0) 
                                   D0fin=0d0
                                   RETURN
                                 else    
c  D(0,t,0,p2sq,0,0)                         
                                   D0fin=0d0
                                   RETURN
                                 endif  !t if
                              else
                                 if(abs(t).lt.tiny) then
c  D(s,0,0,p2sq,0,0)
                                   D0fin=0d0
                                   RETURN
                                 else      
c  D(s,t,0,p2sq,0,0)                                                      
                                    D0fin=D01m_fin(s,t,p2sq,musq)
                                    RETURN
                                 endif  !t if
                              endif    !s if   
                           else   !p4sq
                              if(abs(s).lt.tiny) then
                                 if(abs(t).lt.tiny) then   
c  D(0,0,0,p2sq,0,p4sq) 
                                   D0fin=D00m_fin(p2sq,p4sq,musq)
                                   RETURN
                                 else    
c  D(0,t,0,p2sq,0,p4sq)                         
                                   D0fin=D01m_fin(p2sq,p4sq,t,musq)
                                   RETURN
                                 endif  !t if
                              else
                                 if(abs(t).lt.tiny) then
c  D(s,0,0,p2sq,0,p4sq)
                                   D0fin=D01m_fin(p2sq,p4sq,s,musq)
                                   RETURN
                                 else      
c  D(s,t,0,p2sq,0,p4sq)                                                      
                                     D0fin=D02m_fine(s,t,p2sq,p4sq,musq) 
                                    RETURN
                                 endif  !t if
                              endif    !s if   
                           endif  ! p4sq
                     else    !p3sq
                            if(abs(p4sq).lt.tiny) then
                              if(abs(s).lt.tiny) then
                                 if(abs(t).lt.tiny) then   
c  D(0,0,0,p2sq,p3sq,0) 
                                   D0fin=0d0
                                   RETURN
                                 else    
c  D(0,t,0,p2sq,p3sq,0)                         
                                   D0fin=0d0
                                   RETURN
                                 endif  !t if
                              else
                                 if(abs(t).lt.tiny) then
c  D(s,0,0,p2sq,p3sq,0)
                                   D0fin=0d0
                                   RETURN
                                 else      
c  D(s,t,0,p2sq,p3sq,0)                                                      
                                    D0fin=D02m_fin(t,s,p2sq,p3sq,musq) 
                                    RETURN
                                 endif  !t if
                              endif    !s if   
                             else  ! p4sq
                              if(abs(s).lt.tiny) then
                                 if(abs(t).lt.tiny) then   
c  D(0,0,0,p2sq,p3sq,p4sq) 
                                   D0fin=D01m_fin(p2sq,p4sq,p3sq,musq)
                                   RETURN
                                 else    
c  D(0,t,0,p2sq,p3sq,p4sq)                         
                                   D0fin=D02m_fin(p2sq,p4sq,p3sq,t,musq)
                                   RETURN
                                 endif  !t if
                              else
                                 if(abs(t).lt.tiny) then
c  D(s,0,0,p2sq,p3sq,p4sq)
                                   D0fin=D02m_fin(p4sq,p2sq,s,p3sq,musq)
                                   RETURN
                                 else      
c  D(s,t,0,p2sq,p3sq,p4sq)                                                      
                                   D0fin=D03m_fin(s,t,p2sq,p3sq,p4sq,musq) 
                                   RETURN
                                 endif  !t if
                              endif    !s if   
                             endif  ! p4sq
                      endif !p3sq
           endif   !p2sq
      else    ! p1sq
           if(abs(p2sq).lt.tiny)then
               if(abs(p3sq).lt.tiny)then
                   if(abs(p4sq).lt.tiny) then


                              if(abs(s).lt.tiny) then
                                 if(abs(t).lt.tiny) then   
c  D(0,0,p1sq,0,0,0) 
                                   D0fin=0d0
                                   RETURN
                                 else    
c  D(0,t,0,p1sq,0,0)                         
                                   D0fin=0d0
                                   RETURN
                                 endif  !t if
                              else
                                 if(abs(t).lt.tiny) then
c  D(s,0,0,p1sq,0,0)
                                   D0fin=0d0
                                   RETURN
                                 else      
c  D(s,t,0,p1sq,0,0)                                                      

                             D0fin=D01m_fin(s,t,p1sq,musq)
                             RETURN
                                 endif  !t if
                              endif    !s if   
                   else  ! p4sq
                              if(abs(s).lt.tiny) then
                                 if(abs(t).lt.tiny) then   
c  D(0,0,p1sq,0,0,p4sq) 
                                   D0fin=0d0
                                   RETURN
                                 else    
c  D(0,t,p1sq,0,0,p4sq)                         
                                   D0fin=0d0
                                   RETURN
                                 endif  !t if
                              else
                                 if(abs(t).lt.tiny) then
c  D(s,0,p1sq,0,0,p4sq)
                                   D0fin=0d0
                                   RETURN
                                 else      
c  D(s,t,p1sq,0,0,p4sq)                                                      
                             D0fin=D02m_fin(t,s,p4sq,p1sq,musq) 
                             RETURN
                                 endif  !t if
                              endif    !s if   
                   endif !p4sq
               else  !p3sq
                   if(abs(p4sq).lt.tiny) then
                              if(abs(s).lt.tiny) then
                                 if(abs(t).lt.tiny) then   
c  D(0,0,p1sq,0,p3sq,0) 
                                   D0fin=D00m_fin(p1sq,p3sq,musq)
                                   RETURN
                                 else    
c  D(0,t,p1sq,0,p3sq,0)                         
                                   D0fin=D01m_fin(p1sq,p3sq,t,musq)
                                   RETURN
                                 endif  !t if
                              else
                                 if(abs(t).lt.tiny) then
c  D(s,0,p1sq,0,p3sq,0)
                                   D0fin=D01m_fin(p1sq,p3sq,s,musq)
                                   RETURN
                                 else      
c  D(s,t,p1sq,0,p3sq,0)                                                      
                             D0fin=D02m_fine(s,t,p1sq,p3sq,musq) 
                             RETURN
                                 endif  !t if
                              endif    !s if   
                   else   !p4sq
                              if(abs(s).lt.tiny) then
                                 if(abs(t).lt.tiny) then   
c  D(0,0,p1sq,0,p3sq,p4sq) 
                                   D0fin=D01m_fin(p1sq,p3sq,p4sq,musq)
                                   RETURN
                                 else    
c  D(0,t,p1sq,0,p3sq,p4sq)                         
                                   D0fin=D02m_fin(p1sq,p3sq,t,p4sq,musq)
                                   RETURN
                                 endif  !t if
                              else
                                 if(abs(t).lt.tiny) then
c  D(s,0,p1sq,0,p3sq,p4sq)
                                   D0fin=D02m_fin(p3sq,p1sq,s,p4sq,musq)
                                   RETURN
                                 else      
c  D(s,t,p1sq,0,p3sq,p4sq)                                                      
                             D0fin=D03m_fin(t,s,p3sq,p4sq,p1sq,musq)
                            RETURN
                                 endif  !t if
                              endif    !s if   
                   endif  !p4sq
               endif   !p3sq
          else   !p2sq
               if(abs(p3sq).lt.tiny)then
                   if(abs(p4sq).lt.tiny) then
                              if(abs(s).lt.tiny) then
                                 if(abs(t).lt.tiny) then   
c  D(0,0,p1sq,p2sq,0,0) 
                                   D0fin=0d0
                                   RETURN
                                 else    
c  D(0,t,p1sq,p2sq,0,0)                         
                                   D0fin=0d0
                                   RETURN
                                 endif  !t if
                              else
                                 if(abs(t).lt.tiny) then
c  D(s,0,p1sq,p2sq,0,0)
                                   D0fin=0d0
                                   RETURN
                                 else      
c  D(s,t,p1sq,p2sq,0,0) 
                             D0fin=D02m_fin(s,t,p1sq,p2sq,musq)
                             RETURN
                                 endif  !t if
                              endif    !s if   
                   else !p4sq
                              if(abs(s).lt.tiny) then
                                 if(abs(t).lt.tiny) then   
c  D(0,0,p1sq,p2sq,0,p4sq) 
                                   D0fin=D01m_fin(p2sq,p4sq,p1sq,musq)
                                   RETURN
                                 else    
c  D(0,t,p1sq,p2sq,0,p4sq)                         
                                   D0fin=D02m_fin(p4sq,p2sq,t,p1sq,musq)
                                   RETURN
                                 endif  !t if
                              else
                                 if(abs(t).lt.tiny) then
c  D(s,0,p1sq,p2sq,0,p4sq)
                                   D0fin=D02m_fin(p2sq,p4sq,s,p1sq,musq)
                                   RETURN
                                 else      
c  D(s,t,p1sq,p2sq,0,p4sq)                                                      
                             D0fin=D03m_fin(s,t,p4sq,p1sq,p2sq,musq) 
                                   RETURN
                                 endif  !t if
                              endif    !s if   
                   endif !p4sq
               else   !p3sq
                   if(abs(p4sq).lt.tiny) then

                              if(abs(s).lt.tiny) then
                                 if(abs(t).lt.tiny) then   
c  D(0,0,p1sq,p2sq,p3sq,0) 
                                   D0fin=D01m_fin(p1sq,p3sq,p2sq,musq)
                                   RETURN
                                 else    
c  D(0,t,p1sq,p2sq,p3sq,0)                         
                                   D0fin=D02m_fin(p3sq,p1sq,p2sq,t,musq)
                                   RETURN
                                 endif  !t if
                              else
                                 if(abs(t).lt.tiny) then
c  D(s,0,p1sq,p2sq,p3sq,0)
                                   D0fin=D02m_fin(p1sq,p3sq,p2sq,s,musq)
                                   RETURN
                                 else      
c  D(s,t,p1sq,p2sq,p3sq,0) 
                             D0fin=D03m_fin(t,s,p1sq,p2sq,p3sq,musq) 
                             RETURN
                                 endif  !t if
                              endif    !s if   
                   else  !p4sq
                              if(abs(s).lt.tiny) then
                                 if(abs(t).lt.tiny) then   
c  D(0,0,p1sq,p2sq,p3sq,p4sq) 
                                   D0fin=D02m_fine(p2sq,p4sq,p1sq,p3sq,musq)
                                   RETURN
                                 else    
c  D(0,t,p1sq,p2sq,p3sq,p4sq)                         
                                   D0fin=D03m_fin(p4sq,p2sq,p3sq,t,p1sq,musq)
                                   RETURN
                                 endif  !t if
                              else
                                 if(abs(t).lt.tiny) then
c  D(s,0,p1sq,p2sq,p3sq,p4sq)
                                   D0fin=D03m_fin(p4sq,p2sq,p1sq,s,p3sq,musq)
                                   RETURN
                                 else      
c  D(s,t,p1sq,p2sq,p3sq,p4sq)                                                      
                             D0fin=D04m_fin(s,t,p1sq,p2sq,p3sq,p4sq,musq)
                             RETURN
                                 endif  !t if
                              endif    !s if   
                   endif !p4sq
               endif  !p3sq
           endif  ! p2sq
         endif  !p1sq
         end








CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
CCCCCCCCCCC                SCALAR INTEGRALS                CCCCCCCCCC
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC




c    ***********   D00m_fin(s,t,musq)   *************
c
c Scalar box with MASSLESS PROPAGATORS and with external kinematics: 
c
c                     k
c     p1 ->-----------<-------------<-- p4
c             |                 |
c             |                 |
c             |                 |  
c     p2 ->-------------------------<-- p3
c
c
c    s = (p1+p2)^2
c    t = (p2+p3)^2
c    p1^2 = 0, p2^2 = 0, p3^2 = 0, p4^2 = 0
c    musq = mu^2 = reference dimensional scale 
c 
c  int d^dk/(2 pi)^d 1/(k^2)/(k+p1)^2/(k+p1+p2)^2/(k+p1+p2+p3)^2 = 
c          N_ep * D00m_fin(s,t,musq);
c          N_ep = i/(4 pi)^2 (4 pi)^ep Gamma(1+ep) (musq)^(-ep)

      function D00m_fin(s,t,musq)
      implicit none
      complex * 16 D00m_fin
      real * 8 s,t,m4sq,musq
      real * 8 ms,mt,mm4sq
      complex * 16 lnms,lnmt,lnmm4sq,lnsot
      complex * 16 li2arg1,li2arg2
      real * 8 arg1,arg2
      complex * 16 ipi
      parameter (ipi=(0d0,3.14159265358979323846264338328d0))
      real * 8 pi,pi2, pi2o3t5
      parameter (pi=3.14159265358979323846264338328d0,
     #     pi2 = 9.86960440108935861883449099988d0,
     #   pi2o3t5=16.4493406684822643647241516664602519d0 )

      complex * 16 ris
      real * 8 prefactor, theta, dilog, im_part
      
      if (musq.lt.0d0) then
         write(*,*) 
     #        'POSSIBLE ERROR IN D00m_fin: SCALE MUSQ LESS THAN ZERO!!'
      endif

      prefactor = 1d0/(s*t)
      ms = -s/musq
      mt = -t/musq
      
      
      lnms = log(abs(ms)) - ipi*theta(-ms)
      lnmt = log(abs(mt)) - ipi*theta(-mt)
      lnsot = log(abs(s/t)) + ipi * im_part(s,t)

      ris = lnms**2+lnmt**2-lnsot**2-pi2o3t5

      D00m_fin = prefactor * ris
      end


c    ***********   D01m_fin(s,t,m4sq,musq)   *************
c
c Scalar box with MASSLESS PROPAGATORS and with external kinematics: 
c
c                     k
c     p1 ->-----------<-------------<-- p4
c             |                 |
c             |                 |
c             |                 |  
c     p2 ->-------------------------<-- p3
c
c
c    s = (p1+p2)^2
c    t = (p2+p3)^2
c    p1^2 = 0, p2^2 = 0, p3^2 = 0, p4^2 = m4sq <>0
c    musq = mu^2 = reference dimensional scale 
c 
c  int d^dk/(2 pi)^d 1/(k^2)/(k+p1)^2/(k+p1+p2)^2/(k+p1+p2+p3)^2 = 
c          N_ep * D01m_fin(s,t,m3sq,m4sq,musq);
c          N_ep = i/(4 pi)^2 (4 pi)^ep Gamma(1+ep) (musq)^(-ep)

      function D01m_fin(s,t,m4sq,musq)
      implicit none
      complex * 16 D01m_fin
      real * 8 s,t,m4sq,musq
      real * 8 ms,mt,mm4sq
      complex * 16 lnms,lnmt,lnmm4sq,lnsot
      complex * 16 li2arg1,li2arg2
      real * 8 arg1,arg2
      complex * 16 ipi
      parameter (ipi=(0d0,3.14159265358979323846264338328d0))
      real * 8 pi,pi2, pi2o3t2
      parameter (pi=3.14159265358979323846264338328d0,
     #     pi2 = 9.86960440108935861883449099988d0,
     #   pi2o3t2= 6.57973626739290574588966066658410076d0 )

      complex * 16 ris
      real * 8 prefactor, theta, dilog, im_part
      
      if (musq.lt.0d0) then
         write(*,*) 
     #        'POSSIBLE ERROR IN D01m_fin: SCALE MUSQ LESS THAN ZERO!!'
      endif

      prefactor = 1d0/(s*t)
      ms = -s/musq
      mt = -t/musq
      mm4sq = - m4sq/musq
      
      lnms = log(abs(ms)) - ipi*theta(-ms)
      lnmt = log(abs(mt)) - ipi*theta(-mt)
      lnmm4sq = log(abs(mm4sq)) - ipi*theta(-mm4sq)
      lnsot = log(abs(s/t)) + ipi * im_part(s,t)

      arg1 = 1-m4sq/s
      arg2 = 1-m4sq/t

c     (m4sq/s.lt.0d0)
      if (arg1.gt.1d0) then
         li2arg1 = -dilog(1d0/arg1) - log(arg1)**2/2+pi2/3 
     #        - ipi*log(arg1)*im_part(m4sq,s)
      else
         li2arg1 = dilog(arg1)
      endif

c     (m4sq/t.lt.0d0)
      if (arg2.gt.1d0)  then
         li2arg2 = -dilog(1d0/arg2) - log(arg2)**2/2+pi2/3 
     #        - ipi*log(arg2)*im_part(m4sq,t)
      else
         li2arg2 = dilog(arg2)
      endif

      ris = -lnmm4sq**2+lnms**2+lnmt**2-lnsot**2-2d0*(li2arg1+li2arg2)
     #   -pi2o3t2

      D01m_fin = prefactor * ris
      end









CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
CCCCCCCCCCC                SCALAR INTEGRALS                CCCCCCCCCC
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

c    ***********   D02m_fin(s,t,m3sq,m4sq,musq)   *************
c
c Scalar box with MASSLESS PROPAGATORS and with external kinematics: 
c
c                     k
c     p1 ->-----------<-------------<-- p4
c             |                 |
c             |                 |
c             |                 |  
c     p2 ->-------------------------<-- p3
c
c
c    s = (p1+p2)^2
c    t = (p2+p3)^2
c    p1^2 = 0, p2^2 = 0, p3^2 = m3sq <>0, p4^2 = m4sq <>0
c    musq = mu^2 = reference dimensional scale 
c 
c  int d^dk/(2 pi)^d 1/(k^2)/(k+p1)^2/(k+p1+p2)^2/(k+p1+p2+p3)^2 = 
c          N_ep * D02m_fin(s,t,m3sq,m4sq,musq);
c          N_ep = i/(4 pi)^2 (4 pi)^ep Gamma(1+ep) (musq)^(-ep)

      function D02m_fin(s,t,m3sq,m4sq,musq)
      implicit none
      complex * 16 D02m_fin
      real * 8 s,t,m3sq,m4sq,musq
      real * 8 ms,mt,mm3sq,mm4sq
      complex * 16 lnms,lnmt,lnmm3sq,lnmm4sq,lnsot
      complex * 16 li2arg1,li2arg2
      real * 8 arg1,arg2
      complex * 16 ipi
      parameter (ipi=(0d0,3.14159265358979323846264338328d0))
      real * 8 pi,pi2
      parameter (pi=3.14159265358979323846264338328d0,
     #     pi2 = 9.86960440108935861883449099988d0)

      complex * 16 ris
      real * 8 prefactor, theta, dilog, im_part
      
      if (musq.lt.0d0) then
         write(*,*) 
     #        'POSSIBLE ERROR IN D02m_fin: SCALE MUSQ LESS THAN ZERO!!'
      endif

      prefactor = 1d0/(s*t)
      ms = -s/musq
      mt = -t/musq
      mm3sq = - m3sq/musq
      mm4sq = - m4sq/musq
      
      lnms = log(abs(ms)) - ipi*theta(-ms)
      lnmt = log(abs(mt)) - ipi*theta(-mt)
      lnmm3sq = log(abs(mm3sq)) - ipi*theta(-mm3sq)
      lnmm4sq = log(abs(mm4sq)) - ipi*theta(-mm4sq)
      lnsot = log(abs(s/t)) + ipi * im_part(s,t)

      arg1 = 1-m3sq/t
      arg2 = 1-m4sq/t

c     (m3sq/t.lt.0d0)
      if (arg1.gt.1d0) then
         li2arg1 = -dilog(1d0/arg1) - log(arg1)**2/2+pi2/3 
     #        - ipi*log(arg1)*im_part(m3sq,t)
      else
         li2arg1 = dilog(arg1)
      endif

c     (m4sq/t.lt.0d0)
      if (arg2.gt.1d0)  then
         li2arg2 = -dilog(1d0/arg2) - log(arg2)**2/2+pi2/3 
     #        - ipi*log(arg2)*im_part(m4sq,t)
      else
         li2arg2 = dilog(arg2)
      endif

      ris = lnms**2/2+lnmt**2-lnmm3sq**2/2-lnmm4sq**2/2
     #     +lnmm3sq*lnmm4sq-(lnmm4sq+lnmm3sq-lnms)*lnms
     #     -2*li2arg1-2*li2arg2-lnsot**2 - pi2/6
        
      D02m_fin = prefactor * ris
      end


c    ***********   D02m_fine(s,t,m1sq,m3sq,musq)   *************
c
c Scalar box with MASSLESS PROPAGATORS and with external kinematics: 
c
c                     k
c     p1 ->-----------<-------------<-- p4
c             |                 |
c             |                 |
c             |                 |  
c     p2 ->-------------------------<-- p3
c
c
c    s = (p1+p2)^2
c    t = (p2+p3)^2
c    p1^2 = m1sq, p2^2 = 0, p3^2 = m3sq , p4^2 =0
c    musq = mu^2 = reference dimensional scale 
c 
c  int d^dk/(2 pi)^d 1/(k^2)/(k+p1)^2/(k+p1+p2)^2/(k+p1+p2+p3)^2 = 
c          N_ep * D02m_fin3(s,t,m1sq,m3sq,musq);
c          N_ep = i/(4 pi)^2 (4 pi)^ep Gamma(1+ep) (musq)^(-ep)

      function D02m_fine(s,t,m1sq,m3sq,musq)
      implicit none
      complex * 16 D02m_fine
      real * 8 s,t,m1sq,m3sq,musq
      real * 8 ms,mt,mm1sq,mm3sq
      complex * 16 lnms,lnmt,lnmm1sq,lnmm3sq,lnsot,lnarg5
      complex * 16 li2arg1,li2arg2,li2arg3,li2arg4,li2arg5
      real * 8 arg1,arg2,arg3,arg4,arg5
      complex * 16 ipi
      parameter (ipi=(0d0,3.14159265358979323846264338328d0))
      real * 8 pi,pi2
      parameter (pi=3.14159265358979323846264338328d0,
     #     pi2 = 9.86960440108935861883449099988d0)

      complex * 16 ris
      real * 8 prefactor, theta, dilog, im_part, sgn
      real * 8 im_a, im_b,im_ab, eta3


      if (musq.lt.0d0) then
         write(*,*) 
     #        'POSSIBLE ERROR IN D02m_fine: SCALE MUSQ LESS THAN ZERO!!'
      endif

      prefactor = 1d0/(m1sq*m3sq- s*t)
      ms = -s/musq
      mt = -t/musq
      mm1sq = - m1sq/musq
      mm3sq = - m3sq/musq
      
      lnms = log(abs(ms)) - ipi*theta(-ms)
      lnmt = log(abs(mt)) - ipi*theta(-mt)
      lnmm1sq = log(abs(mm1sq)) - ipi*theta(-mm1sq)
      lnmm3sq = log(abs(mm3sq)) - ipi*theta(-mm3sq)
      lnsot = log(abs(s/t)) + ipi*im_part(s,t)
c      write(*,*) "lnsot", lnsot
 
      arg1 = 1-m1sq/t
      arg2 = 1-m3sq/t
      arg3 = 1-m1sq/s
      arg4 = 1-m3sq/s
      arg5= 1-(m1sq*m3sq)/(s*t)

c$$$      write(*,*) "arg1", arg1
c$$$      write(*,*) "arg2", arg2
c$$$      write(*,*) "arg3", arg3
c$$$      write(*,*) "arg4", arg4
c$$$      write(*,*) "arg5", arg5

c     (m1sq/t.lt.0d0)
      if (arg1.gt.1d0) then
         li2arg1 = -dilog(1d0/arg1) - log(arg1)**2/2+pi2/3 
     #        - ipi*log(arg1)*im_part(m1sq,t)
      else
         li2arg1 = dilog(arg1)
      endif
     

c     (m3sq/t.lt.0d0)
      if (arg2.gt.1d0)  then
         li2arg2 = -dilog(1d0/arg2) - log(arg2)**2/2+pi2/3 
     #        - ipi*log(arg2)*im_part(m3sq,t)
      else
         li2arg2 = dilog(arg2)
      endif
     


c     (m1sq/s.lt.0d0)
      if (arg3.gt.1d0) then
         li2arg3 = -dilog(1d0/arg3) - log(arg3)**2/2+pi2/3 
     #        - ipi*log(arg3)*im_part(m1sq,s)
      else
         li2arg3 = dilog(arg3)
      endif

c     (m3sq/s.lt.0d0)
      if (arg4.gt.1d0)  then
         li2arg4 = -dilog(1d0/arg4) - log(arg4)**2/2+pi2/3 
     #        - ipi*log(arg4)*im_part(m3sq,s)
      else
         li2arg4 = dilog(arg4)
      endif

      im_a = m1sq/s*(1d0/m1sq-1d0/s)
      im_b = m3sq/t*(1d0/m3sq-1d0/t)
      im_ab =m1sq*m3sq/(s*t)*(1d0/m1sq-1d0/s+1d0/m3sq-1d0/t)


       

c     m1sq*m3sq/(s*t)  < 0
      if (arg5.gt.1d0) then
         li2arg5 = -dilog(1d0/arg5) - log(arg5)**2/2+pi2/3 
     #        - ipi*log(arg5)*sgn(im_ab)
      else
         li2arg5 = dilog(arg5)
      endif

      if (arg5.lt.0d0) then
         lnarg5 = log(abs(arg5)) - ipi*sgn(im_ab)
      else
         lnarg5 = log(arg5)
      endif
      
c      li2arg5 = li2arg5  + 2*ipi*eta3(im_a,im_b,im_ab)*lnarg5


c$$$      write(*,*) "li2arg1", li2arg1
c$$$      write(*,*) "li2arg2", li2arg2
c$$$      write(*,*) "li2arg3", li2arg3
c$$$      write(*,*) "li2arg4", li2arg4
c$$$      write(*,*) "li2arg5", li2arg5
c$$$      write(*,*) "********************"
c$$$      write(*,*) "eta3(im_a,im_b,im_ab)",  eta3(im_a,im_b,im_ab)
c$$$      write(*,*) "********************"
c$$$      write(*,*) "lnms", lnms
c$$$      write(*,*) "lnmt", lnmt
c$$$      write(*,*) "lnmms1sq", lnmm1sq
c$$$      write(*,*) "lnmms3sq", lnmm3sq
c$$$      write(*,*) "lnsot", lnsot
 
      ris = lnmm1sq**2+lnmm3sq**2-lnms**2-lnmt**2+lnsot**2+
     #  2.d0*(li2arg1+li2arg2+li2arg3+li2arg4-li2arg5)   

c      write(*,*) "ris", ris

      D02m_fine = prefactor*ris
      end




c    ***********   D03m_fin(s,t,m2sq,m3sq,m4sq,musq)   *************
c
c Scalar box with MASSLESS PROPAGATORS and with external kinematics: 
c
c
c                     k
c     p1 ->-----------<-------------<-- p4
c             |                 |
c             |                 |
c             |                 |  
c     p2 ->-------------------------<-- p3
c
c
c    s = (p1+p2)^2
c    t = (p2+p3)^2
c    p1^2 = 0, p2^2 = m2sq <>0, p3^2 = m3sq <>0, p4^2 = m4sq <>0
c    musq = mu^2 = reference dimensional scale 
c 
c  int d^dk/(2 pi)^d 1/(k^2)/(k+p1)^2/(k+p1+p2)^2/(k+p1+p2+p3)^2 =  
c          N_ep * D03m_fin(s,t,m2sq,m3sq,m4sq,musq);
c          N_ep = i/(4 pi)^2 (4 pi)^ep Gamma(1+ep) (musq)^(-ep)

      function D03m_fin(s,t,m2sq,m3sq,m4sq,musq)
      implicit none
      complex * 16 D03m_fin
      real * 8 s,t,m2sq,m3sq,m4sq,musq
      complex * 16 D03m_fin_BDK,D03m_fin_DN
      D03m_fin = D03m_fin_BDK(s,t,m2sq,m3sq,m4sq,musq)
c      D03m_fin = D03m_fin_DN(s,t,m2sq,m3sq,m4sq,musq)
      
      end




      function D03m_fin_BDK(s,t,m2sq,m3sq,m4sq,musq)
      implicit none
      complex * 16 D03m_fin_BDK
      real * 8 s,t,m2sq,m3sq,m4sq,musq
      real * 8 ms,mt,mm2sq,mm3sq,mm4sq,sot
      complex * 16 lnms,lnmt,lnmm2sq,lnmm3sq,lnmm4sq,lnsot,lnarg3
      complex * 16 li2arg1,li2arg2,li2arg3
      real * 8 arg1,arg2,arg3
      complex * 16 ipi
      parameter (ipi=(0d0,3.14159265358979323846264338328d0))
      real * 8 pi,pi2
      parameter (pi=3.14159265358979323846264338328d0,
     #     pi2 = 9.86960440108935861883449099988d0)

      complex * 16 ris
      real * 8 prefactor, theta, dilog, sgn, im_part
      real * 8 im_a, im_b,im_ab, eta3

      if (musq.lt.0d0) then
         write(*,*) 
     #        'POSSIBLE ERROR IN D03m_fin: SCALE MUSQ LESS THAN ZERO!!'
      endif


c  1/(ms*mt-mm2s*mm4s)*
c  (1/2*ln(ms)^2+1/2*ln(mt)^2-1/2*ln(mm2s)^2-1/2*ln(mm4s)^2-ln(ms/mt)^2
c  +ln(mm2s)*ln(mm3s)-(ln(mm3s)+ln(mm2s)-ln(mt))*ln(mt)
c  +ln(mm3s)*ln(mm4s)-(ln(mm4s)+ln(mm3s)-ln(ms))*ln(ms)
c  -2*li2(1-mm2s/ms)
c  -2*li2(1-mm4s/mt)
c  +2*li2(1-mm2s*mm4s/ms/mt))

      prefactor = 1d0/(s*t-m2sq*m4sq)
      ms = -s/musq
      mt = -t/musq
      mm2sq = - m2sq/musq
      mm3sq = - m3sq/musq
      mm4sq = - m4sq/musq
      sot = s/t
      
      lnms = log(abs(ms)) - ipi*theta(-ms)
      lnmt = log(abs(mt)) - ipi*theta(-mt)
      lnmm2sq = log(abs(mm2sq)) - ipi*theta(-mm2sq)
      lnmm3sq = log(abs(mm3sq)) - ipi*theta(-mm3sq)
      lnmm4sq = log(abs(mm4sq)) - ipi*theta(-mm4sq)
      lnsot = log(abs(sot)) + ipi * im_part(s,t)

     
      arg1 = 1-m2sq/s
      arg2 = 1-m4sq/t
      arg3 = 1-m2sq*m4sq/s/t
      

c     m2sq/s < 0 
      if (arg1.gt.1d0) then
         li2arg1 = -dilog(1d0/arg1) - log(arg1)**2/2+pi2/3 
     #        - ipi*log(arg1)*im_part(m2sq,s)
      else
         li2arg1 = dilog(arg1)
      endif

c     m4sq/t < 0
      if (arg2.gt.1d0) then
         li2arg2 = -dilog(1d0/arg2) - log(arg2)**2/2+pi2/3 
     #        - ipi*log(arg2)*im_part(m4sq,t)
      else
         li2arg2 = dilog(arg2)
      endif


      im_a = m2sq/s*(1d0/m2sq-1d0/s)
      im_b = m4sq/t*(1d0/m4sq-1d0/t)
      im_ab = m2sq*m4sq/s/t*(1d0/m2sq-1d0/s+1d0/m4sq-1d0/t)


c     m2sq*m4sq/(s*t)  < 0
      if (arg3.gt.1d0) then
         li2arg3 = -dilog(1d0/arg3) - log(arg3)**2/2+pi2/3 
     #        - ipi*log(arg3)*sgn(im_ab)
      else
         li2arg3 = dilog(arg3)
      endif

      if (arg3.lt.0d0) then
         lnarg3 = log(abs(arg3)) - ipi*sgn(im_ab)
      else
         lnarg3 = log(arg3)
      endif
c Not checked, this line, but without it, the results does not agree with
C QCD looptools      

      li2arg3 = li2arg3  + 2*ipi*eta3(im_a,im_b,im_ab)*lnarg3

      
      ris = lnms**2/2+lnmt**2/2-lnmm2sq**2/2-lnmm4sq**2/2
     #     +lnmm2sq*lnmm3sq
     #     -(lnmm3sq+lnmm2sq-lnmt)*lnmt+lnmm3sq*lnmm4sq
     #     -(lnmm4sq+lnmm3sq-lnms)*lnms      
      
      ris = ris - 2*li2arg1 - 2*li2arg2 + 2*li2arg3 
     #     - lnsot**2
        
      D03m_fin_BDK = prefactor * ris


      end





c    ***********   D04m_fin(s,t,m1sq,m2sq,m3sq,m4sq,musq)   *************
c
c Scalar box with MASSLESS PROPAGATORS and with external kinematics: 
c
c
c                     k
c     p1 ->-----------<-------------<-- p4
c             |                 |
c             |                 |
c             |                 |  
c     p2 ->-------------------------<-- p3
c
c
c    s = (p1+p2)^2
c    t = (p2+p3)^2
c    p1^2 = m1sq <>0, p2^2 = m2sq <>0, p3^2 = m3sq <>0, p4^2 = m4sq <>0
c    musq = mu^2 = reference dimensional scale  (NOT used here, since 
c                 this box is finite 
c 
c  int d^dk/(2 pi)^d 1/(k^2)/(k+p1)^2/(k+p1+p2)^2/(k+p1+p2+p3)^2 =  
c          N_ep * D04m_fin(s,t,m1sq,m2sq,m3sq,m4sq,musq);
c          N_ep = i/(4 pi)^2 (4 pi)^ep Gamma(1+ep) (musq)^(-ep)
c

c	m2 ... arbitrary, but positive -> set m2 = 1d0


c
c
c     D0 with 4 nonzero external masses
c     and zero internal propagators      
c     use Eq.(41) from Denner, Nierste, Scharf		     
c

      function D04m_fin(s,t,m1sq,m2sq,m3sq,m4sq,musq)
      implicit none
      complex * 16 D04m_fin
      real * 8 s,t,m1sq,m2sq,m3sq,m4sq,musq

      complex * 16 ipi
      parameter (ipi=(0d0,3.14159265358979323846264338328d0))
      real * 8 pi,pi2
      parameter (pi=3.14159265358979323846264338328d0,
     #     pi2 = 9.86960440108935861883449099988d0)

      integer k

      real * 8  dilog, eta3, theta, sgn      
      complex * 16 ln_mx(2),lnk12,lnk13,lnk14,lnk23,total

      real * 8 m2,k12,k13,k14,k23,k24,k34
      real * 8 a,b,c,d,sqrtdelta,x(2)
      real * 8 k34ok13,im_k34ok13,k24ok12,im_k24ok12
      real * 8 im_x(2),im_k34ok13_x(2),im_k24ok12_x(2),arg,logarg
      complex * 16 eta_k34ok13(2),eta_k24ok12(2)
      complex * 16 lik34ok13(2), logk34ok13(2)
      complex * 16 lik24ok12(2), logk24ok12(2)
      complex * 16 D04m_fin_real

      m2 = 1d0                  ! > 0 ALWAYS!
      
      k12 = -m1sq/m2
      k13 = -s/m2
      k14 = -m4sq/m2
      k23 = -m2sq/m2
      k24 = -t/m2
      k34 = -m3sq/m2
      

      a = k24*k34
      b = k13*k24+k12*k34-k14*k23
      c = k12*k13
      d = k23
      

      if ((b**2-4*a*c).lt.0d0) then
c         write(*,*) 'ERROR: delta less than zero!!'
c         write(*,*) 'RETURN ZERO'
c use BDK version of D04m, since in this case it's REAL!!!
c         D04m_fin = D04m_fin_real(s,t,m1sq,m2sq,m3sq,m4sq,musq)
         D04m_fin = 0d0
c         write(*,*) 's,t,m1sq,m2sq,m3sq,m4sq'
c         write(*,*) s,t,m1sq,m2sq,m3sq,m4sq
         return         
c         stop
      else
         sqrtdelta = sqrt(b**2-4*a*c)
      endif


      x(1) = (-b + sqrtdelta)/(2d0*a)
      x(2) = (-b - sqrtdelta)/(2d0*a)

c      write(*,*) 'x(1) x(2)',x(1),x(2)

      k34ok13 = k34/k13
      im_k34ok13 = (1d0-m3sq/s)/s
      
      k24ok12 = k24/k12
      im_k24ok12 = (1d0-t/m1sq)/m1sq

      im_x(1) = m2sq/sqrtdelta
      im_x(2) = -im_x(1)


      do k=1,2
         im_k34ok13_x(k) = im_k34ok13*x(k) + k34ok13*im_x(k)
         im_k24ok12_x(k) = im_k24ok12*x(k) + k24ok12*im_x(k)
         eta_k34ok13(k) = 
     #        2*ipi*eta3(-im_x(k),im_k34ok13,-im_k34ok13_x(k))
         eta_k24ok12(k) = 2*ipi*
     #        eta3(-im_x(k),im_k24ok12,-im_k24ok12_x(k))
         ln_mx(k) = log(abs(x(k))) - theta(x(k))*sgn(im_x(k))*ipi
      enddo

      lnk12 = log(abs(k12))-ipi*theta(-k12)
      lnk13 = log(abs(k13))-ipi*theta(-k13)
      lnk14 = log(abs(k14))-ipi*theta(-k14)
      lnk23 = log(abs(k23))-ipi*theta(-k23)

c      write(*,*) lnk12,lnk13,lnk14,lnk23

      do k=1,2
         arg = 1d0 + k34ok13*x(k)
         if (arg.gt.1d0) then
c     complex dilog
            logarg = log(arg)
            lik34ok13(k) = -dilog(1/arg) - 0.5d0 * logarg**2 + pi2/3
     #           + ipi * sgn(im_k34ok13_x(k)) * logarg
         else
            lik34ok13(k) = dilog(arg)
         endif
         logk34ok13(k) = log(abs(arg)) 
     #        + ipi*theta(-arg) * sgn(im_k34ok13_x(k))
c         write(*,*) lik34ok13(k),logk34ok13(k)
      enddo


      do k=1,2
         arg = 1d0 + k24ok12*x(k)
         if (arg.gt.1d0) then
c     complex dilog
            logarg = log(arg)
            lik24ok12(k) = -dilog(1/arg) - 0.5d0 * logarg**2 + pi2/3
     #           + ipi * sgn(im_k24ok12_x(k)) * logarg
         else
            lik24ok12(k) = dilog(arg)
         endif
         logk24ok12(k) = log(abs(arg)) 
     #        + ipi*theta(-arg) * sgn(im_k24ok12_x(k))
      enddo



c      write(*,*) '=========>', 
c     #     (lnk12+lnk13-lnk14-lnk23)*(-ln_mx(1)+ln_mx(2))

      total = 0d0
      do k=1,2
         total = total + (-1)**k*( - 0.5d0 *ln_mx(k)**2 
     #        - lik34ok13(k) - eta_k34ok13(k) * logk34ok13(k)
     #        - lik24ok12(k) - eta_k24ok12(k) * logk24ok12(k))
      enddo
      
      total = total + (lnk12+lnk13-lnk14-lnk23)*(-ln_mx(1)+ln_mx(2))

      D04m_fin = total/a/m2**2/(x(1)-x(2))
      
      end





c   Bubble correction: finite part.           
c
c    musq = mu^2 e' una scala dimensionale esterna!!
c 
c  int d^dk/(2 pi)^d 1/(k^2)/(k+p1)^2 = 
c          N_ep * B0fin(p1sq,musq)
c          N_ep = i/(4 pi)^2 (4 pi)^ep Gamma(1+ep) (musq)^(-ep)

      function B0fin(p1sq,musq)
      implicit none
      complex * 16 B0fin
      real * 8 p1sq,musq
c     ris(1) = finite part, ris(2) = coeff of 1/ep, 
c     ris(3) = coeff of 1/ep^2      
      complex * 16 ris(3)
      complex * 16 l1
      complex * 16 I
      parameter (I=(0,1))
      double precision pi
      parameter (pi=3.141592653589793238462643383279502884197D0)
      complex * 16 ipi
      parameter (ipi=(0,3.141592653589793238462643383279502884197D0))
      integer j,offshellleg


      complex * 16 B0
      logical debugB0C0D0
      common/debug_B0C0D0/debugB0C0D0
      real * 8 tiny
      parameter (tiny=1d-6)
      
c      if (debugB0C0D0) then
c         B0fin = B0(p1sq,musq)
c         return
c      endif

      if (musq.lt.0) then
         write(*,*) 'ERROR in B0fin: mu^2 MUST be a positive number'
         stop
      endif

c      if (p1sq.eq.0) then
      if (abs(p1sq).lt.tiny) then
c BUBBLE WITH ZERO EXTERNAL INVARIANTS
c         write(*,*) 'Warning: B0fin called with external'//
c     #       ' invariant equal to zero'
         ris(3) = 0.d0
         ris(2) = 0.d0
         ris(1) = 0.d0
         B0fin = ris(1)
         return
      endif
         
      l1 =  log(abs(p1sq/musq))
      if (p1sq.gt.0) then
         l1 = l1 -ipi
      endif
      ris(3) = 0d0
      ris(2) = 1d0
      ris(1) = 2-l1
      B0fin = ris(1)
      end




c     C0fin(p1sq,p2sq,p3sq,musq)
c
c             
c              /|===== p3^2
c             / | 
c            /  |
c           /   |
c p1^2 =====    | 
c           \   |
c            \  |
c             \ |
c              \|===== p2^2
c             
c
c    musq = mu^2 e' una scala dimensionale esterna!!
c 
c  int d^dk/(2 pi)^d 1/(k^2)/(k+p1)^2/(k+p1+p2)^2 = 
c          N_ep * C0fin(p1sq,p2sq,p3sq,musq)
c          N_ep = i/(4 pi)^2 (4 pi)^ep Gamma(1+ep) (musq)^(-ep)
      function C0fin(p1sq,p2sq,p3sq,musq)
      implicit none
      complex * 16 C0fin
      real * 8 p1sq,p2sq,p3sq,musq
c     ris(1) = finite part, ris(2) = coeff of 1/ep, 
c     ris(3) = coeff of 1/ep^2      
      complex * 16 ris(3)
      real * 8 qsq(3),tmp(3)
      real * 8 dilog
      complex * 16 lr,l1,l2,lr2,lr3
      complex * 16 I
      parameter (I=(0,1d0))
      double precision pi
      parameter (pi=3.141592653589793238462643383279502884197D0)
      complex * 16 ipi
      parameter (ipi=(0,3.141592653589793238462643383279502884197D0))
      integer j,offshellleg,imax,ii
      complex * 16 C03
      real * 8 r3,r2,detsq,det,x,y,max,lomx,lomy
c      complex*16 C0t

      real * 8 tiny
      parameter (tiny=1d-7)
      

      complex * 16 C0
      logical debugB0C0D0
      common/debug_B0C0D0/debugB0C0D0

c      if (debugB0C0D0) then
c         C0fin = C0(p1sq,p2sq,p3sq,musq)
c         return
c      endif


      if (musq.lt.0) then
         write(*,*) 'ERROR in C0fin: mu^2 MUST be a positive number'
         stop
      endif
      
c      c0fin = C0t(p1sq,p2sq,p3sq,musq)
c      return
      
c      write(*,*) 'C0 with args',p1sq,p2sq,p3sq,musq

      offshellleg = 0
      do j=1,3 
         qsq(j) = 0.d0
      enddo
c      if (p1sq.ne.0) then
      if (abs(p1sq).gt.tiny) then
         offshellleg = offshellleg +1
         qsq(offshellleg) = p1sq
      endif
c      if (p2sq.ne.0) then
      if (abs(p2sq).gt.tiny) then
         offshellleg = offshellleg +1
         qsq(offshellleg) = p2sq
      endif
c      if (p3sq.ne.0) then
      if (abs(p3sq).gt.tiny) then
         offshellleg = offshellleg +1
         qsq(offshellleg) = p3sq
      endif

      if (offshellleg.eq.1) then
c     TRIANGLE WITH ONLY ONE EXTERNAL INVARIANT      
         l1 =  log(abs(qsq(1)/musq))
         if (qsq(1).lt.0) then
c     do nothing
         else
            l1 = l1 -ipi
         endif
         ris(3) = 1.d0
         ris(2) = -l1
         ris(1) = 1.d0/2*l1**2 - pi**2/6         
         do j=1,3
            ris(j) = ris(j)/qsq(1)
         enddo

      elseif (offshellleg.eq.2) then
c TRIANGLE WITH TWO EXTERNAL INVARIANTS  
         lr = log(abs(qsq(2)/qsq(1)))
         l1 = log(abs(qsq(1)/musq))
         l2 = log(abs(qsq(2)/musq))
         if ((qsq(1).lt.0).and.(qsq(2).lt.0)) then
c     do nothing            
         elseif ((qsq(1).gt.0).and.(qsq(2).lt.0)) then
            lr = lr + ipi
            l1 = l1 - ipi
         elseif ((qsq(1).lt.0).and.(qsq(2).gt.0)) then
            lr = lr - ipi
            l2 = l2 - ipi
         elseif ((qsq(1).gt.0).and.(qsq(2).gt.0)) then
            l1 = l1 - ipi
            l2 = l2 - ipi
         endif
         ris(3) = 0.d0
         ris(2) = lr
         ris(1) = 1.d0/2*(l1**2 - l2**2)
         do j=1,3
            ris(j) = ris(j)/(qsq(1)-qsq(2))
         enddo
      elseif (offshellleg.eq.3) then         
c TRIANGLE WITH THREE EXTERNAL INVARIANTS           
         ris(3) = 0.d0
         ris(2) = 0.d0

c     ris(1) = C03(qsq(1),qsq(2),qsq(3),0d0,0d0,0d0)
         
         
c order the qsq(i) with absolute max value in the first position
         max = 0d0
         do ii=1,3
            if (abs(qsq(ii)).gt.max) then
               max = abs(qsq(ii))
               imax = ii
            endif
         enddo
         tmp(1) = qsq(imax)

         j = 2
         do ii=1,3
            if (ii.ne.imax) then               
               tmp(j) = qsq(ii)
               j = j+1
            endif
         enddo
         do ii=1,3
            qsq(ii) = tmp(ii)
         enddo

c         write(*,*) 'qsq(i)  ============> ', qsq

         r3 = qsq(3)/qsq(1)
         r2 = qsq(2)/qsq(1)

c         if (sqrt(r2)+sqrt(r3).gt.1d0) then 
c            write(*,*) 'NOT YET IMPLEMENTED'
c            stop
c         endif
         
         detsq = (1-r2-r3)**2 - 4*r2*r3
         if (detsq.lt.0d0) then
c            write(*,*) qsq,detsq
c            write(*,*) r2,r3
c            write(*,*) 'WARNING: this case has NOT yet been implemented',detsq
c            write(*,*) 'RETURN 0 from C0fin function'
            C0fin = 0d0
            detsq = 0
            RETURN
         endif
         
         det = sqrt(detsq)
         x = 1.d0/2/r2*(r3+r2-1+det)
         y = 1.d0/2/r3*(r3+r2-1+det)

c         write(*,*) 'x, y ==> ',x,y

         if ((x.gt.1d0).or.(y.gt.1d0)) then             
c            write(*,*) 'ERROR in C0fin: x and/or y have values 
c     #           bigger than one', x, y
c            stop
            C0fin = 0d0
            detsq = 0
            RETURN
         endif
        
         lomx = log(1-x)
         lomy = log(1-y)
         lr2 = log(abs(r2))
         lr3 = log(abs(r3))
         if (r2.lt.0d0) then
            if (qsq(1).lt.0) then
               lr2 = lr2 - ipi
            else
               lr2 = lr2 + ipi
            endif
         endif
         if (r3.lt.0d0) then
            if (qsq(1).lt.0) then
               lr3 = lr3 - ipi
            else
               lr3 = lr3 + ipi
            endif
         endif

         
         
         ris(1) = 1/qsq(1)*(1-x)*(1-y)/(1-x*y)*(2*dilog(x)+2*dilog(y)+
     #        (lomx+lomy)**2+2*lr2*lomy+2*lr3*lomx+lr3*lr2+Pi**2/3)
 

      else
c TRIANGLE WITH ZERO EXTERNAL INVARIANTS
c         write(*,*) 'Warning: C0fin called with external '//
c     #        'invariants equal to zero'
c         write(*,*) 'dot prods',p1sq,p2sq,p3sq,musq
         ris(3) = 0.d0
         ris(2) = 0.d0
         ris(1) = 0.d0
      endif

      C0fin = ris(1)
      end





















ccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
cccccccccccccccccccccccccccccccccccccccccccccccccccccccc
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
cccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c
c
c
c
      function E0(psq, qsq, lsq, tsq,
     #     pq, pl, pt, ql, qt, lt, 
     #     D0_2345,D0_1345,D0_1245,D0_1235,D0_1234)
      implicit none
      complex * 16 E0
      real * 8 psq, qsq, lsq, tsq, pq, pl, pt, ql, qt, lt
      real * 8 q1s, q2s, ppp, pq1, ppq1, pq2, ppq2, q1q2, pps
      real * 8 det, coeff(1:5)
      complex * 16 D0_2345,D0_1345,D0_1245,D0_1235,D0_1234, J(1:5)
      integer i
      real * 8 tiny
      logical lpr
      common /e0print/lpr

      tiny = 1d-4

      if (abs(psq).gt.tiny) then
         write(*,*) 'p^2 =',psq
         write(*,*) 'E0 computed only with p^2 and pp^2 = 0'
         stop
      endif

      q1s  = qsq
      q2s  = lsq
      ppp  = psq+pq+pl+pt
      pq1  = pq
      ppq1 = pq + qsq+ql+qt
      pq2  = pl
      ppq2 = pl + ql + lsq + lt
      q1q2 = ql
      pps = lsq+2*lt+2*ql+2*pl+tsq+2*qt+2*pt+qsq+2*pq+psq
      if (abs(pps).gt.tiny) then
         write(*,*) 'pp^2 =',pps
         write(*,*) 'E0 computed only with p^2 and pp^2 = 0'
         lpr = .true.
c         stop
      endif


      det = -8*ppp*(4*q1s*pq2**2*ppq1+q1s**2*q2s*ppp-4*pq1*pq2*q1s*ppq2-
     #4*pq1*pq2*ppq1*q2s+4*q1s*pq2*q1q2*ppp+4*q1s*pq2*q2s*ppp-4*pq1*q2s*
     #q1q2*ppp+4*q1s*q2s*pq1*ppp-8*pq1*pq2*ppq1*q1q2-8*pq1*q1q2*pq2*ppp+
     #4*q1s*pq1*q1q2*ppq2+2*q1s*pq1*q2s*ppq2-8*pq1*q2s*ppq1*q1q2+2*q2s*p
     #q2*q1s*ppq1+4*q1s*pq2*ppq1*q1q2+2*q1s*q1q2*q2s*ppp-2*q1s**2*pq2*pp
     #q2+4*pq1**2*q2s*ppq2+8*pq1**2*q1q2*ppq2+4*q1s*pq2**2*ppp-8*q1q2**2
     #*pq1*ppp+4*pq1**2*q2s*ppp+q1s*q2s**2*ppp-8*q1q2**2*pq1*ppq1-2*pq1*
     #q2s**2*ppq1)

      
c      write(*,*)  'DETERMINANTE in E0',det

      coeff(1)=(-16*pq1**2*q2s*ppp-8*q1s*pq1*ppq2**2+16*q1q2**2*ppq1*ppp
     #-8*pq1**2*q2s*ppq2-16*pq1**2*q1q2*ppq2-16*q1q2*pq2*ppq1**2-8*q1s*p
     #q2**2*ppq1-16*q1s*pq2**2*ppp+16*q1q2**2*pq1*ppp+16*q1q2**2*ppp**2-
     #8*q2s*pq2*ppq1**2-16*q1s*q2s*ppp**2+8*pq1*pq2*q1s*ppq2+8*pq1*pq2*p
     #pq1*q2s-8*q1s*pq2*q1q2*ppp-8*q1s*pq2*q2s*ppp+8*pq1*q2s*q1q2*ppp-8*
     #q1s*q2s*pq1*ppp+8*q1s*pq2*ppq1*ppq2-8*q1s*q1q2*ppq2*ppp-8*q1s*q2s*
     #ppq1*ppp-8*q1s*q2s*ppq2*ppp+8*pq1*q2s*ppq1*ppq2+8*q1q2*q2s*ppq1*pp
     #p+16*pq1*pq2*ppq1*q1q2+16*q1s*pq2*ppq2*ppp+16*pq1*q2s*ppq1*ppp+32*
     #pq1*q1q2*pq2*ppp-16*pq1*q1q2*ppq2*ppp-16*q1q2*pq2*ppq1*ppp+16*pq1*
     #q1q2*ppq1*ppq2)/det


      coeff(2)=(16*pq1*q2s*ppp**2+8*q2s*pq2*ppq1*ppp+16*pq2**2*ppq1**2-1
     #6*pq1*pq2*ppq2*ppp-16*q1q2*pq2*ppp**2+16*pq1**2*ppq2**2-32*pq1*pq2
     #*ppq1*ppq2+16*pq1*q2s*ppq1*ppp-8*q1q2*q2s*ppp**2+8*q1s*q2s*ppp**2+
     #16*pq2**2*ppq1*ppp+8*pq1*q2s*ppq2*ppp-16*q1q2**2*ppp**2+8*q1s*pq1*
     #ppq2**2-16*q1q2**2*ppq1*ppp+8*q2s*pq2*ppq1**2-8*q1s*pq2*ppq1*ppq2-
     #16*pq1*q1q2*ppq1*ppq2+8*q1s*q1q2*ppq2*ppp+16*q1q2*pq2*ppq1**2+8*q1
     #s*q2s*ppq1*ppp+8*q1s*q2s*ppq2*ppp-8*pq1*q2s*ppq1*ppq2-8*q1q2*q2s*p
     #pq1*ppp)/det


      coeff(3)=8*ppp*(-2*pq1*q1q2*ppp+2*q1s*pq2*ppq2-2*pq1*q2s*ppq1-4*pq
     #1*q1q2*ppq1-2*pq1*q1q2*ppq2-2*pq1*pq2*ppq1+q1s*q1q2*ppp+2*pq1**2*p
     #pq2+2*q1s*pq2*ppp+q1s*pq1*ppq2+q1s*pq2*ppq1-pq1*q2s*ppq2-q2s*pq2*p
     #pq1+2*pq1*pq2*ppq2+2*q1q2**2*ppp+q1q2*q2s*ppp-2*pq2**2*ppq1-2*pq1*
     #q2s*ppp-2*q1q2*pq2*ppq1+2*q1q2*pq2*ppp)/det


      coeff(4)=-8*ppp*(-2*pq1*q1q2*ppp-2*pq1*q2s*ppq1-4*pq1*q1q2*ppq1-2*
     #pq1*pq2*ppq1+q1s*q2s*ppp+q1s*q1q2*ppp+2*pq1**2*ppq2+2*q1s*pq2*ppp+
     #q1s*pq1*ppq2+q1s*pq2*ppq1)/det


      coeff(5)=(16*pq1**2*q1q2*ppq2+8*pq1**2*q2s*ppq2+16*q1s*pq2**2*ppp-
     #8*pq1*pq2*q1s*ppq2-16*pq1*pq2*ppq1*q1q2-8*pq1*pq2*ppq1*q2s+8*q1s*p
     #q2*q1q2*ppp+8*q1s*pq2*q2s*ppp-8*pq1*q2s*q1q2*ppp-32*pq1*q1q2*pq2*p
     #pp-16*q1q2**2*pq1*ppp+8*q1s*pq2**2*ppq1+16*pq1**2*q2s*ppp+8*q1s*q2
     #s*pq1*ppp)/det


c D0(q1,q2,pp-p-q1-q2) 
      J(1) = D0_2345
c D0(p+q1,q2,pp-p-q1-q2)
      J(2) = D0_1345
c D0(p,q1+q2,pp-p-q1-q2)
      J(3) = D0_1245
c D0(p,q1,pp-p-q1)
      J(4) = D0_1235
c D0(p,q1,q2)
      J(5) = D0_1234
      
      E0 = 0d0
      do i=1,5
         E0 = E0 + J(i)*coeff(i)
      enddo
      end







    

c this function compute the imaginary part of a/b, where:
c (a+I*ep)/(b+I*ep) = a/b + I*ep*a/b*(1/a-1/b)
c  and return the sign of this imaginary part
      function sgnratio(a,b)
      implicit none
      real * 8 a,b,sgnratio
      if ((1d0/b-a/b**2).gt.0d0) then
         sgnratio = +1d0
      else
         sgnratio = -1d0          
      endif

      end


c     imaginary part of a/b  
c     IF a/b > 0 then returns 0
c     IF a/b < 0 then
c        if a < 0 then return +1
c        else return -1
      function im_part(a,b)
      implicit none
      real * 8 a,b,im_part
      if ((a/b).gt.0d0) then
         im_part = 0d0
      else
         if (a.lt.0d0) then
            im_part = +1d0
         else
            im_part = -1d0
         endif
      endif
      
      end


c     compute  Li(1-a/b) 
c     IF a/b < 0 then, it attaches a POSITIVE +I*ep to a or b, 
c     according to which of them is negative
      function dilog_1_m_aob(a,b)
      implicit none
      complex * 16 dilog_1_m_aob
      real * 8 a,b
      real * 8 arg
      complex * 16 ipi
      parameter (ipi=(0d0,3.14159265358979323846264338328d0))
      real * 8 pi,pi2
      parameter (pi=3.14159265358979323846264338328d0,
     #     pi2 = 9.86960440108935861883449099988d0)
      real * 8  dilog, sign

      arg = 1d0-a/b
      if (arg.gt.1d0) then
         if (a.lt.0d0) then
            sign = +1d0
         else
            sign = -1d0
         endif
         dilog_1_m_aob = -dilog(1d0/arg) - log(arg)**2/2+pi2/3 
     #        - ipi*log(arg)*sign
      else
         dilog_1_m_aob = dilog(arg)
      endif
      
      end


      function eta3(im_a,im_b,im_ab)
      implicit none
      real * 8 eta3,im_a,im_b,im_ab
      real * 8 theta
      eta3 = theta(-im_a)*theta(-im_b)*theta(im_ab)
     #     - theta(im_a)*theta(im_b)*theta(-im_ab)
      end

      function theta(x)
      implicit none
      real * 8 theta,x
      if (x.gt.0d0) then
         theta = 1d0
      else
         theta = 0d0
      endif
      end

      function sgn(x)
      implicit none
      real * 8 sgn,x
      if (x.gt.0d0) then
         sgn = +1d0
      else
         sgn = -1d0
      endif
      end



      function delta(x,y)
      implicit none
      real * 8 delta
      integer x,y
      if (x.eq.y) then
         delta = 1d0
      else
         delta = 0d0
      endif
      end





c$$$c E0fin
c$$$c Author: Francisco Campanario
c$$$c Date:2008-04-16
c$$$c Generalized the old one by Carlo Oleari
c$$$C This is valid for general kinematics pisq<=>0
c$$$
c$$$      function E0fin(p1s,p2s,p3s,p4s,p5s,s12,s23,s34,s45,s15, 
c$$$     #     D02345,D01345,D01245,D01235,D01234)
c$$$      implicit none
c$$$      complex * 16 E0fin,D02345,D01345,D01245,D01235,D01234
c$$$      real * 8 p1s,p2s,p3s,p4s,p5s,s12,s23,s34,s45,s15
c$$$      real * 8 p1s2,p2s2,p3s2,p4s2,p5s2,s12s,s23s,s34s,s45s,s15s
c$$$      real * 8 d,x1,x2,x3,x4,x5
c$$$
c$$$      p1s2=p1s*p1s
c$$$      p2s2=p2s*p2s
c$$$      p3s2=p3s*p3s
c$$$      p4s2=p4s*p4s
c$$$      p5s2=p5s*p5s
c$$$      s12s=s12*s12
c$$$      s23s=s23*s23
c$$$      s34s=s34*s34
c$$$      s45s=s45*s45
c$$$      s15s=s15*s15
c$$$
c$$$       d=2*(p1s2*p3s*p4s*s34+p2s2*p4s*p5s*s45+p2s*(p3s*p5s*(p5s*s23-s1
c$$$     -   5*s45)+s34*s45*(-(p5s*s23)+s15*s45)-p4s*s12*(p5s*s23+s15*s45
c$$$     -   ))+s12*(p3s*s15*(-(p5s*s23)+s15*s45)+s23*(p4s*s12*s15+p5s*s2
c$$$     -   3*s34-s15*s34*s45))-p1s*(-(p3s2*p5s*s15)+s23*s34*(p4s*s12-s3
c$$$     -   4*s45)+p2s*p4s*(p3s*p5s-p4s*s12+s34*s45)+p3s*(p4s*s12*s15+p5
c$$$     -   s*s23*s34+s15*s34*s45)))
c$$$       x1=-(p3s2*p5s*s15)+p3s*p4s*s12*s15+p3s2*s15s-p3s*s12*s15s+p3s*p
c$$$     -   5s*s15*s23+p3s*s12*s15*s23-2*p4s*s12*s15*s23+p3s*p5s*s23*s34
c$$$     -   +p4s*s12*s23*s34-2*p3s*s15*s23*s34+s12*s15*s23*s34-p5s*s23s*
c$$$     -   s34-s12*s23s*s34+p1s*(-(p3s2*s15)+s23*(p4s-s34)*s34+p2s*p4s*
c$$$     -   (p3s-p4s+s34)+p3s*(p4s*(s15-2*s34)+(s15+s23)*s34))+s23s*s34s
c$$$     -   +p2s2*p4s*(p4s-p5s-s45)-p3s*s15s*s45+p3s*s15*s34*s45+s15*s23
c$$$     -   *s34*s45-s23*s34s*s45+p2s*(-(p4s2*s12)+p3s*(p4s*(p5s-2*s15)+
c$$$     -   p5s*s15-2*p5s*s23+s15*s45)+s34*(p5s*s23+(-2*s15+s23)*s45)+p4
c$$$     -   s*(p5s*s23+s12*(s15+s23)-2*s23*s34+s15*s45+s34*s45))
c$$$       x2=p3s2*p5s2-2*p3s*p4s*p5s*s12+p4s2*s12s-p3s2*p5s*s15+p3s*p4s*s
c$$$     -   12*s15+p3s*p5s*s12*s15-p4s*s12s*s15-p3s*p5s2*s23+p3s*p5s*s12
c$$$     -   *s23+p4s*p5s*s12*s23-p4s*s12s*s23+p3s*p5s*s23*s34+p4s*s12*s2
c$$$     -   3*s34-2*p5s*s12*s23*s34+p3s*p5s*s15*s45-2*p3s*s12*s15*s45+p4
c$$$     -   s*s12*s15*s45-2*p3s*p5s*s34*s45-2*p4s*s12*s34*s45+p3s*s15*s3
c$$$     -   4*s45+s12*s15*s34*s45+p5s*s23*s34*s45+s12*s23*s34*s45-s23*s3
c$$$     -   4s*s45+p2s*(-(p4s2*s12)+s34*(p5s-s45)*s45+p3s*p5s*(p4s-p5s+s
c$$$     -   45)+p4s*(p5s*(s12-2*s45)+(s12+s34)*s45))+p1s*(-(p3s2*p5s)-(p
c$$$     -   4s-s34)*(p4s*s12-s34*s45)+p3s*(p4s*(p5s+s12-2*s34)+s34*(p5s+
c$$$     -   s45)))-s15*s34*s45s+s34s*s45s
c$$$       x3=-(p3s*p5s2*s23)+p4s*p5s*s12*s23+p3s*p5s*s15*s23-2*p4s*s12*s1
c$$$     -   5*s23+p5s*s12*s15*s23+p5s2*s23s-p5s*s12*s23s+p1s2*p4s*(-p3s+
c$$$     -   p4s-s34)-p5s*s23s*s34+p3s*p5s*s15*s45+p4s*s12*s15*s45-p3s*s1
c$$$     -   5s*s45-s12*s15s*s45-2*p5s*s15*s23*s45+s12*s15*s23*s45+p5s*s2
c$$$     -   3*s34*s45+s15*s23*s34*s45+p1s*(-(p4s2*s12)+p4s*s12*s15-2*p4s
c$$$     -   *p5s*s23+p4s*s12*s23+p4s*s23*s34+p5s*s23*s34-2*p4s*s15*s45+p
c$$$     -   4s*s34*s45+s15*s34*s45-2*s23*s34*s45+p2s*p4s*(-p4s+p5s+s45)+
c$$$     -   p3s*(-2*p5s*s15+p4s*(p5s+s15)+p5s*s23+s15*s45))+p2s*(-((p5s-
c$$$     -   s45)*(p5s*s23-s15*s45))+p4s*(p5s*(s23-2*s45)+s15*s45))+s15s*
c$$$     -   s45s-s15*s34*s45s
c$$$       x4=p1s2*s34*(-p3s-p4s+s34)+p2s2*p5s*(-p4s+p5s-s45)+p2s*(-2*p5s*
c$$$     -   s12*s15+p3s*p5s*(-p5s+s15)+p4s*s12*(p5s+s15)-p5s2*s23+p5s*s1
c$$$     -   2*s23+p5s*s23*s34+p5s*s15*s45+s12*s15*s45+p5s*s34*s45-2*s15*
c$$$     -   s34*s45)+s12*(-(p4s*s12*s15)+p3s*(p5s-s15)*s15+s12*s15s+p5s*
c$$$     -   s15*s23-s12*s15*s23-2*p5s*s23*s34+s15*s23*s34-s15s*s45+s15*s
c$$$     -   34*s45)+p1s*(p4s*s12*s15+p4s*s12*s34-2*s12*s15*s34+p5s*s23*s
c$$$     -   34+s12*s23*s34+p3s*(s15*(s12+s34)+p5s*(-2*s15+s34))-s23*s34s
c$$$     -   +s15*s34*s45-s34s*s45+p2s*(p3s*p5s+p4s*(p5s-2*s12+s34)+s34*(
c$$$     -   -2*p5s+s45)))
c$$$       x5=p1s2*p3s*(p3s-p4s-s34)+p2s2*s45*(-p4s-p5s+s45)+s12*(p3s*(p5s
c$$$     -   *s23+s15*(s23-2*s45))+s23*(-(p4s*s12)-s12*s15-p5s*s23+s12*s2
c$$$     -   3-s23*s34+s15*s45+s34*s45))+p1s*(-(p3s2*(p5s+s15))+s23*(p4s*
c$$$     -   s12+s12*s34-2*s34*s45)+p2s*(-2*p4s*s12+p3s*(p4s+p5s-2*s45)+p
c$$$     -   4s*s45+s34*s45)+p3s*(p4s*s12+s12*s15+p5s*s23-2*s12*s23+s23*s
c$$$     -   34+s15*s45+s34*s45))+p2s*(p5s*s12*s23+s12*s15*s45+p5s*s23*s4
c$$$     -   5-2*s12*s23*s45+s23*s34*s45+p4s*s12*(s23+s45)+p3s*(s15*s45+p
c$$$     -   5s*(-2*s23+s45))-s15*s45s-s34*s45s)
c$$$
c$$$       E0fin=-(x1*D02345+x2*D01345+x3*D01245+x4*D01235+x5*D01234)/d
c$$$     
c$$$      end
c$$$
c$$$
c$$$
c$$$
c$$$
c$$$
c$$$











c$$$c D0fin(s,t,p1sq,p2sq,p3sq,p4sq,musq)
c$$$c Author: Francisco Campanario
c$$$c Date:2008-04-16
c$$$c I have used Carlo Oleari D02m,D03m,D04m
c$$$c Bug solve in analytical continuation of D03m
c$$$c Need  my  D00m,D02me,D01m functions
c$$$c
c$$$c                     k
c$$$c     p1 ->-----------<-------------<-- p4
c$$$c             |                 |
c$$$c             |                 |
c$$$c             |                 |  
c$$$c     p2 ->-------------------------<-- p3
c$$$c
c$$$c
c$$$c    s=(p1+p2)^2
c$$$c    t=(p2+p3)^2
c$$$c    musq = mu^2 reference dimensional scale
c$$$c 
c$$$c  int d^dk/(2 pi)^d 1/(k^2)/(k+p1)^2/(k+p1+p2)^2/(k+p1+p2+q1)^2= 
c$$$c          N_ep * D0fin(s,t,q1^2,q2^2)
c$$$c          N_ep = i/(4 pi)^2 (4 pi)^ep Gamma(1+ep) (musq)^(-ep)
c$$$c   the external momenta could be whaterver value
c$$$
c$$$
c$$$      function D0fin(s,t,p1sq,p2sq,p3sq,p4sq,musq)
c$$$      implicit none
c$$$      complex * 16 D0fin
c$$$      real * 8 s,t,musq,p1sq,p2sq,p3sq,p4sq
c$$$      Complex*16 D01m_fin,D02m_fin,D02m_fine,D03m_fin,D04m_fin,D00m_fin
c$$$      External D01m_fin,D02m_fin,D02m_fine,D03m_fin,D04m_fin,D00m_fin
c$$$      real * 8 tiny
c$$$      parameter (tiny=1d-6)
c$$$      
c$$$      if(abs(p1sq).lt.tiny) then
c$$$          if(abs(p2sq).lt.tiny)then
c$$$               if(abs(p3sq).lt.tiny)then
c$$$                   if(abs(p4sq).lt.tiny) then
c$$$c  D(s,t,0,0,0,0)                                                      
c$$$                             D0fin=D00m_fin(s,t,musq)
c$$$                             RETURN
c$$$                   else
c$$$c D(s,t,0,0,0,p4sq) 
c$$$                            D0fin=D01m_fin(s,t,p4sq,musq)
c$$$                              RETURN
c$$$                   endif 
c$$$               else 
c$$$                   if(abs(p4sq).lt.tiny) then
c$$$c D(s,t,0,0,p3sq,0)
c$$$                             D0fin=D01m_fin(s,t,p3sq,musq)
c$$$                             RETURN 
c$$$                   else
c$$$c D(s,t,0,0,p3sq,p4sq)
c$$$                             D0fin=D02m_fin(s,t,p3sq,p4sq,musq)
c$$$                             RETURN
c$$$                   endif
c$$$               endif
c$$$          else
c$$$                     if(abs(p3sq).lt.tiny)then
c$$$                           if(abs(p4sq).lt.tiny) then
c$$$c D(s,t,0,p2sq,0,0)
c$$$                                    D0fin=D01m_fin(s,t,p2sq,musq)
c$$$                             RETURN
c$$$                           else
c$$$c D(s,t,0,p2sq,0,p4sq)
c$$$                                     D0fin=D02m_fine(s,t,p2sq,p4sq,musq) 
c$$$                             RETURN
c$$$                           endif 
c$$$                     else 
c$$$                            if(abs(p4sq).lt.tiny) then
c$$$c D(s,t,0,p2sq,p3sq,0)
c$$$                                      D0fin=D02m_fin(t,s,p2sq,p3sq,musq) 
c$$$                             RETURN
c$$$                             else 
c$$$c D(s,t,0,p2sq,p3sq,p4sq)
c$$$                                      D0fin=D03m_fin(s,t,p2sq,p3sq,p4sq,musq) 
c$$$                             endif
c$$$                      endif
c$$$           endif 
c$$$      else
c$$$           if(abs(p2sq).lt.tiny)then
c$$$               if(abs(p3sq).lt.tiny)then
c$$$                   if(abs(p4sq).lt.tiny) then
c$$$c  D(s,t,p1sq,0,0,0) 
c$$$                             D0fin=D01m_fin(s,t,p1sq,musq)
c$$$                             RETURN
c$$$                   else
c$$$c  D(s,t,p1sq,0,0,p4sq) 
c$$$                             D0fin=D02m_fin(t,s,p4sq,p1sq,musq) 
c$$$                             RETURN
c$$$                   endif 
c$$$               else 
c$$$                   if(abs(p4sq).lt.tiny) then
c$$$c  D(s,t,p1sq,0,p3sq,0) 
c$$$                             D0fin=D02m_fine(s,t,p1sq,p3sq,musq) 
c$$$                             RETURN
c$$$                   else
c$$$c  D(s,t,p1sq,0,p3sq,p4sq) 
c$$$                             D0fin=D03m_fin(t,s,p3sq,p4sq,p1sq,musq)
c$$$                            RETURN
c$$$                   endif
c$$$               endif
c$$$          else
c$$$               if(abs(p3sq).lt.tiny)then
c$$$                   if(abs(p4sq).lt.tiny) then
c$$$c  D(s,t,p1sq,p2sq,0,0) 
c$$$                             D0fin=D02m_fin(s,t,p1sq,p2sq,musq)
c$$$                             RETURN
c$$$                   else
c$$$c  D(s,t,p1sq,p2sq,0,p4sq) 
c$$$                             D0fin=D03m_fin(s,t,p4sq,p1sq,p2sq,musq) 
c$$$                             RETURN
c$$$                   endif 
c$$$               else 
c$$$                   if(abs(p4sq).lt.tiny) then
c$$$c  D(s,t,p1sq,p2sq,p3sq,0) 
c$$$                             D0fin=D03m_fin(t,s,p1sq,p2sq,p3sq,musq) 
c$$$                             RETURN
c$$$                   else
c$$$c  D(s,t,p1sq,p2sq,p3sq,p4sq) 
c$$$                             D0fin=D04m_fin(s,t,p1sq,p2sq,p3sq,p4sq,musq)
c$$$                             RETURN
c$$$                   endif
c$$$               endif
c$$$           endif
c$$$         endif
c$$$         end
c$$$
c$$$

