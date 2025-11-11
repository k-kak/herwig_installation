c
c  <I(eps)> for q(pa) + q(pb) -> q(p1)+q(p2) + H
c
      subroutine ComputeColorCorrsHjj(cT1T2,cTaT1,cTaT2,cTaTb,cTbT1,
     $     cTbT2,msq,born)
      implicit none 
      double precision cT1T2,cTaT1,cTaT2,cTaTb,cTbT1,cTbT2,msq,CA,CF

      double complex born(2),m1,m2

      m1 = born(1)
      m2 = born(2)

      CA = 3.0d0
      CF = 4.0d0/3.0d0

      cT1T2 = -(CA*CF*m2*DConjg(m1)) - CA*CF*m1*DConjg(m2)
      
      cTaT1 = CA*CF*m2*DConjg(m1) + m1*(-(CA**2*CF*DConjg(m1)) + 
     $     CA*CF*DConjg(m2))
      cTaT2 = CA*CF*m1*DConjg(m2) + m2*(CA*CF*DConjg(m1) - 
     $     CA**2*CF*DConjg(m2))

      cTaTb = -(CA*CF*m2*DConjg(m1)) - CA*CF*m1*DConjg(m2)

      cTbT1 = CA*CF*m1*DConjg(m2) + m2*(CA*CF*DConjg(m1) - 
     $     CA**2*CF*DConjg(m2))
      cTbT2 = CA*CF*m2*DConjg(m1) + m1*(-(CA**2*CF*DConjg(m1)) + 
     $     CA*CF*DConjg(m2))

      msq = m1*(CA**2*Dconjg(m1) - CA*Dconjg(m2)) + 
     $     m2*(-(CA*Dconjg(m1)) + CA**2*Dconjg(m2))
      

      end
c
c

c 
      subroutine CSIepsHjj(pa,pb,p1,p2,born,mu2,result,div)
      implicit none 
      double precision cT1T2, cTaT1, cTaT2, 
     $     cTaTb, cTbT1, cTbT2,msq
      double complex born(2)

      integer div
      double precision pa(0:3),pb(0:3),p1(0:3),p2(0:3)
      double precision mu2
      double complex result

      call ComputeColorCorrsHjj(cT1T2,cTaT1,cTaT2,cTaTb,cTbT1,
     $     cTbT2,msq,born)
      
      call CSIeps1Hjj(cT1T2,cTaT1,cTaT2,cTaTb,cTbT1,
     $     cTbT2,msq,
     $     pa,pb,p1,p2,mu2,result,div)

      

      end

      subroutine CSIeps1Hjj(cT1T2,cTaT1,cTaT2,cTaTb,cTbT1,
     $     cTbT2,msq,
     $     pa,pb,p1,p2,mu2,result,div)
      implicit none

      double precision cT1T2,cTaT1,cTaT2, 
     $     cTaTb, cTbT1, cTbT2,msq

      

      integer div
      double precision pa(0:3),pb(0:3),p1(0:3),p2(0:3)
      double precision mu2
      double complex result
    
      double precision CA,CF,gammaQ,gammaG,KG,KQ,Pi,NF,TR
      double precision ga,gb,g1,g2,g3
      parameter(Pi=3.141592653589793d0)
      double precision sij
      double complex MyLog
      external MyLog
      external sij


      CA = 3.0d0
      CF = 4.0d0/3.0d0
      NF = 5.0d0 ! 5 flavors ?
      TR = 1.0d0/2.0d0

      gammaQ = 3.0d0/2.0d0*CF
      gammaG = 11.0d0/6.0d0*CA-2.0d0/3.0d0*TR*NF
C      KG = (67.0d0/18.0d0 - Pi**2/6.0d0)*CA-10.0d0/9.0d0*TR*NF
C      KQ = (7.0d0/2.0d0 - Pi**2/6.0d0)*CF


      ga = gammaQ
      gb = gammaQ
      g1 = gammaQ
      g2 = gammaQ

c      m1 = born(1)
c      m2 = born(2)
c      m3 = born(3)
c      m4 = born(4)

c      mu2 = mu2/2.0d0 ! rescale there is bug


c      print*,'div=',div
c      print*,'mu=',dsqrt(mu2)

c finite part
      if(div.eq.2) then
         result = 2*(cT1T2 + cTaT1 + cTaT2 + 
     $         cTaTb + cTbT1 + cTbT2 )*dcmplx(1.0d0,0.0d0)
c         print*,"result(2)=",result
      elseif(div.eq.1) then
c         
c     change to the V convention 
c
         result =  ((-g1 - g2 - gA - gB)*msq**dcmplx(1.0d0,0.0d0) + 
     $        2*(cT1T2*MyLog(mu2/Sij(p1,p2)) + 
     $        cTaT1*MyLog(mu2/Sij(p1,pa)) + 
     $        cTbT1*MyLog(mu2/Sij(p1,pb)) + 
     $        cTaT2*MyLog(mu2/Sij(p2,pa)) + 
     $        cTbT2*MyLog(mu2/Sij(p2,pb)) +  
     $        cTaTb*MyLog(mu2/Sij(pa,pb))))

c         print*,"result(1)=",result
c         print*,"sij(pa,p3)/mu2=",sij(pa,p3)/mu2

      else
         result = 0.0d0
      endif

      end



      subroutine ConvertVirtHjj(pa,pb,p1,p2,born,Qscale2,
     $     reslt,div)
      implicit none
      integer div
      double complex reslt
      double complex born(2)
      double precision cT1T2, cTaT1, cTaT2, 
     $     cTaTb, cTbT1, cTbT2, msq

      double precision pa(0:3),pb(0:3),p1(0:3),p2(0:3)
      double precision Qscale2

      double precision sij
      double complex MyLog
      external MyLog
      external sij

       call ComputeColorCorrsHjj(cT1T2,cTaT1,cTaT2,cTaTb,cTbT1,
     $     cTbT2,msq,born)

     


      if(div.eq.0) then ! assume q2 = mur2
         reslt = -MyLog(Qscale2/sij(p1,p2))**2*cT1T2 -   
     $        MyLog(Qscale2/sij(pa,p1))**2*cTaT1 -
     $        MyLog(Qscale2/sij(pa,p2))**2*cTaT2 -   
     $        MyLog(Qscale2/sij(pa,pb))**2*cTaTb -
     $        MyLog(Qscale2/sij(pb,p1))**2*cTbT1 -
     $        MyLog(Qscale2/sij(pb,p2))**2*cTbT2 
   
         
      elseif(div.eq.1) then

         reslt = - 2.0d0*(cT1T2*MyLog(Qscale2/sij(p1,p2))  +
     $        cTaT1*MyLog(Qscale2/sij(pa,p1))  +
     $        cTaT2*MyLog(Qscale2/sij(pa,p2))  +
     $        cTaTb*MyLog(Qscale2/sij(pa,pb))  +
     $        cTbT1*MyLog(Qscale2/sij(pb,p1))  +
     $        cTbT2*MyLog(Qscale2/sij(pb,p2))) 
     

      elseif(div.eq.2) then
         reslt = 0.0d0
      else
         stop

      endif
      
      end

c
C 
C     This subroutine computes <I(epsilon)> 
C     A factor -alpha_{s}(mu^2)/(2 pi) * (4*pi)^{epsilon}/Gamma[1-epsilon]
C     has been factored out. 
C     q(pa) q(pb) -> q(p1) q(p2) g(p3) H

      subroutine ComputeColorCorrs(cT1T2, cT1T3, cT2T3, cTaT1, cTaT2, 
     $     cTaT3, cTaTb, cTbT1, cTbT2, cTbT3, msq,born)

      implicit none 

      double precision cT1T2, cT1T3, cT2T3, cTaT1, cTaT2, cTaT3, 
     $     cTaTb, cTbT1, cTbT2, cTbT3, msq,CA,CF

      double complex born(4),mLO1,mLO2,mLO3,mLO4
      logical lremovediag

      lremovediag = .false.

      mLO1 = born(1)
      mLO2 = born(2)
      mLO3 = born(3)
      mLO4 = born(4)

      CA = 3.0d0
      CF = 4.0d0/3.0d0

    


      cTaTb=mLO4*(-(CF*DCONJG(mLO1))/2. + CA*CF**2*DCONJG(mLO2) + 
     $     (CA*CF*DCONJG(mLO3))/2.) + mLO1*((CA*CF*DCONJG(mLO2))/2. + 
     $     CA*CF**2*DCONJG(mLO3) - (CF*DCONJG(mLO4))/2.) + 
     $     mLO3*(CA*CF**2*DCONJG(mLO1) - (CF*DCONJG(mLO2))/2. + 
     $     (CA*CF*DCONJG(mLO4))/2.) + mLO2*((CA*CF*DCONJG(mLO1))/2. - 
     $     (CF*DCONJG(mLO3))/2. + CA*CF**2*DCONJG(mLO4))

c      cTbTa = cTaTb

      cTaT1=mLO4*((CF*DCONJG(mLO1))/2. - CA*CF**2*DCONJG(mLO2) - 
     $     (CA*CF*DCONJG(mLO3))/2.) + mLO1*((CA*CF*DCONJG(mLO1))/2. + 
     $     (CF*DCONJG(mLO3))/2. + (CF*DCONJG(mLO4))/2.) + 
     $     mLO3*((CF*DCONJG(mLO1))/2. - CA*CF**2*DCONJG(mLO2) - 
     $     (CA*CF*DCONJG(mLO4))/2.) + 
     $     mLO2*(-(CA**2*CF**2*DCONJG(mLO2)) - CA*CF**2*DCONJG(mLO3) - 
     $     CA*CF**2*DCONJG(mLO4))

c      cT1Ta = cTaT1

      cTaT2=mLO3*(-(CA*CF**2*DCONJG(mLO1)) - CA*CF**2*DCONJG(mLO2) - 
     $     CA**2*CF**2*DCONJG(mLO3)) + mLO2*(-(CA*CF*DCONJG(mLO1))/2. - 
     $     CA*CF**2*DCONJG(mLO3) + (CF*DCONJG(mLO4))/2.) + 
     $     mLO1*(-(CA*CF*DCONJG(mLO2))/2. - CA*CF**2*DCONJG(mLO3) + 
     $     (CF*DCONJG(mLO4))/2.) + mLO4*((CF*DCONJG(mLO1))/2. + 
     $     (CF*DCONJG(mLO2))/2. + (CA*CF*DCONJG(mLO4))/2.)

c      cT2Ta = cTaT2

      cTaT3=mLO3*(-(CA**2*CF*DCONJG(mLO1))/2. + 
     $     (CA**2*CF*DCONJG(mLO2))/2.) + 
     $     mLO1*(-(CA**3*CF*DCONJG(mLO1))/2. - 
     $     (CA**2*CF*DCONJG(mLO3))/2. - (CA**2*CF*DCONJG(mLO4))/2.) + 
     $     mLO2*((CA**2*CF*DCONJG(mLO3))/2. - 
     $     (CA**2*CF*DCONJG(mLO4))/2.) + 
     $     mLO4*(-(CA**2*CF*DCONJG(mLO1))/2. - (CA**2*CF*DCONJG(mLO2))/2. - 
     $     (CA**3*CF*DCONJG(mLO4))/2.)

      
c      cT3Ta = cTaT3

      cTbT1=mLO3*((CF*DCONJG(mLO1))/2. + (CF*DCONJG(mLO2))/2. + 
     $     (CA*CF*DCONJG(mLO3))/2.) + mLO2*(-(CA*CF*DCONJG(mLO1))/2. + 
     $     (CF*DCONJG(mLO3))/2. - CA*CF**2*DCONJG(mLO4)) + 
     $     mLO1*(-(CA*CF*DCONJG(mLO2))/2. + (CF*DCONJG(mLO3))/2. - 
     $     CA*CF**2*DCONJG(mLO4)) + mLO4*(-(CA*CF**2*DCONJG(mLO1)) - 
     $     CA*CF**2*DCONJG(mLO2) - CA**2*CF**2*DCONJG(mLO4))


c      cT1Tb = cTbT1

      cTbT2=mLO4*(-(CA*CF**2*DCONJG(mLO1)) + (CF*DCONJG(mLO2))/2. - 
     $     (CA*CF*DCONJG(mLO3))/2.) + mLO2*((CA*CF*DCONJG(mLO2))/2. + 
     $     (CF*DCONJG(mLO3))/2. + (CF*DCONJG(mLO4))/2.) + 
     $     mLO3*(-(CA*CF**2*DCONJG(mLO1)) + (CF*DCONJG(mLO2))/2. - 
     $     (CA*CF*DCONJG(mLO4))/2.) + 
     $     mLO1*(-(CA**2*CF**2*DCONJG(mLO1)) - CA*CF**2*DCONJG(mLO3) - 
     $     CA*CF**2*DCONJG(mLO4))

c      cT2Tb = cTbT2

      cTbT3=mLO4*((CA**2*CF*DCONJG(mLO1))/2. - 
     $     (CA**2*CF*DCONJG(mLO2))/2.) + 
     $     mLO3*(-(CA**2*CF*DCONJG(mLO1))/2. - 
     $     (CA**2*CF*DCONJG(mLO2))/2. - (CA**3*CF*DCONJG(mLO3))/2.) + 
     $     mLO2*(-(CA**3*CF*DCONJG(mLO2))/2. - (
     $     CA**2*CF*DCONJG(mLO3))/2. - (CA**2*CF*DCONJG(mLO4))/2.) + 
     $     mLO1*(-(CA**2*CF*DCONJG(mLO3))/2. + (CA**2*CF*DCONJG(mLO4))/2.)


c      cT3Tb = cTbT3

      cT1T2=mLO4*(CA*CF**2*DCONJG(mLO1) - (CF*DCONJG(mLO2))/2. + 
     $     (CA*CF*DCONJG(mLO3))/2.) + mLO2*((CA*CF*DCONJG(mLO1))/2. + 
     $     CA*CF**2*DCONJG(mLO3) - (CF*DCONJG(mLO4))/2.) + 
     $     mLO3*(-(CF*DCONJG(mLO1))/2. + CA*CF**2*DCONJG(mLO2) + 
     $     (CA*CF*DCONJG(mLO4))/2.) + mLO1*((CA*CF*DCONJG(mLO2))/2. - 
     $     (CF*DCONJG(mLO3))/2. + CA*CF**2*DCONJG(mLO4))

c      cT2T1 = cT1T2

      cT1T3=mLO4*(-(CA**2*CF*DCONJG(mLO1))/2. + 
     $     (CA**2*CF*DCONJG(mLO2))/2.) + 
     $     mLO3*(-(CA**2*CF*DCONJG(mLO1))/2. - 
     $     (CA**2*CF*DCONJG(mLO2))/2. - (CA**3*CF*DCONJG(mLO3))/2.) + 
     $     mLO1*(-(CA**3*CF*DCONJG(mLO1))/2. - 
     $     (CA**2*CF*DCONJG(mLO3))/2. - (CA**2*CF*DCONJG(mLO4))/2.) + 
     $     mLO2*(-(CA**2*CF*DCONJG(mLO3))/2. + (CA**2*CF*DCONJG(mLO4))/2.)

c      cT3T1 = cT1T3

      cT2T3=mLO3*((CA**2*CF*DCONJG(mLO1))/2. - 
     $     (CA**2*CF*DCONJG(mLO2))/2.) + 
     $     mLO2*(-(CA**3*CF*DCONJG(mLO2))/2. - 
     $     (CA**2*CF*DCONJG(mLO3))/2. - (CA**2*CF*DCONJG(mLO4))/2.) + 
     $     mLO1*((CA**2*CF*DCONJG(mLO3))/2. - 
     $     (CA**2*CF*DCONJG(mLO4))/2.) + 
     $     mLO4*(-(CA**2*CF*DCONJG(mLO1))/2. - 
     $     (CA**2*CF*DCONJG(mLO2))/2. - (CA**3*CF*DCONJG(mLO4))/2.)

c      cT3T2 = cT2T3

        msq = CA*CF*((CA*mLO1 + mLO3 + mLO4)*DCONJG(mLO1) + 
     $     (CA*mLO2 + mLO3 + mLO4)*DCONJG(mLO2) + 
     $    (mLO1 + mLO2 + CA*mLO3)*DCONJG(mLO3) + 
     $     (mLO1 + mLO2 + CA*mLO4)*DCONJG(mLO4))
c       

     

c     set these to zero if pentagons are off
c     These cancle

c      print*,"cTaTb =", cTaTb
c      print*,"cT1T2 =", cT1T2
c      print*,"cT2Ta =", cTaT2
c      print*,"cT1Tb =", cTbT1

c      print*,"cTaT1 =", cTaT1
c      print*,"cTbT2 =", cTbT2
c      print*,"cT1T3 =", cT1T3
c      print*,"cT2T3 =", cT2T3
c      print*,"cTaT3=",cTaT3
c      print*,"cTbT3=",cTbT3
     
c     For checking pure boxes: 
c      cTaTb = 0.0d0
c      cTbTa = 0.0d0
c      cTbT1 = 0.0d0 
c      cTbT2 = 0.0d0 
c      cTaT2 = 0.0d0 
cc      cTaT2 = 0.0d0 
c      cT1T2 = 0.0d0
cc      cT1T2 = 0.0d0

c      cTaT2 = 0.0d0
c      cTbT1 = 0.0d0
c      cTaT3 = 0.0d0
c      cTbT3 = 0.0d0
c      cT1T3 = 0.0d0
c      cT2T3 = 0.0d0
c      msq = 0.0d0


      end


      subroutine ConvertVirt(pa,pb,p1,p2,p3,born,Qscale2,
     $     reslt,div)
      implicit none
      integer div
      double complex reslt
      double complex born(4)
      double precision cT1T2, cT1T3, cT2T3, cTaT1, cTaT2, cTaT3, 
     $     cTaTb, cTbT1, cTbT2, cTbT3,msq

      double precision pa(0:3),pb(0:3),p1(0:3),p2(0:3),p3(0:3)
      double precision Qscale2

      double precision sij
      double complex MyLog
      external MyLog
      external sij


      call ComputeColorCorrs(cT1T2, cT1T3, cT2T3, cTaT1, cTaT2, cTaT3, 
     $     cTaTb, cTbT1, cTbT2, cTbT3,msq,born)


      if(div.eq.0) then ! assume q2 = mur2
         reslt = -MyLog(Qscale2/sij(p1,p2))**2*cT1T2 -
     $        MyLog(Qscale2/sij(p1,p3))**2*cT1T3 -
     $        MyLog(Qscale2/sij(p2,p3))**2*cT2T3 -
     $        MyLog(Qscale2/sij(pa,p1))**2*cTaT1 -
     $        MyLog(Qscale2/sij(pa,p2))**2*cTaT2 -
     $        MyLog(Qscale2/sij(pa,p3))**2*cTaT3 -
     $        MyLog(Qscale2/sij(pa,pb))**2*cTaTb -
     $        MyLog(Qscale2/sij(pb,p1))**2*cTbT1 -
     $        MyLog(Qscale2/sij(pb,p2))**2*cTbT2 -
     $        MyLog(Qscale2/sij(pb,p3))**2*cTbT3 
         
      elseif(div.eq.1) then

         reslt = - 2.0d0*(cT1T2*MyLog(Qscale2/sij(p1,p2))  +
     $        cT1T3*MyLog(Qscale2/sij(p1,p3))  +
     $        cT2T3*MyLog(Qscale2/sij(p2,p3))  +
     $        cTaT1*MyLog(Qscale2/sij(pa,p1))  +
     $        cTaT2*MyLog(Qscale2/sij(pa,p2))  +
     $        cTaT3*MyLog(Qscale2/sij(pa,p3))  +
     $        cTaTb*MyLog(Qscale2/sij(pa,pb))  +
     $        cTbT1*MyLog(Qscale2/sij(pb,p1))  +
     $        cTbT2*MyLog(Qscale2/sij(pb,p2))  +
     $        cTbT3*MyLog(Qscale2/sij(pb,p3)) )
         

      elseif(div.eq.2) then
         reslt = 0.0d0
      else
         stop

      endif
      
      end

      subroutine CSIeps(pa,pb,p1,p2,p3,born,mu2,NF,result,div)
      implicit none 
      double precision NF ! number of flavors
      double precision cT1T2, cT1T3, cT2T3, cTaT1, cTaT2, cTaT3, 
     $     cTaTb, cTbT1, cTbT2, cTbT3,msq
      double complex born(4)

      integer div
      double precision pa(0:3),pb(0:3),p1(0:3),p2(0:3),p3(0:3)
      double precision mu2
      double complex result

      call ComputeColorCorrs(cT1T2, cT1T3, cT2T3, cTaT1, cTaT2, cTaT3, 
     $     cTaTb, cTbT1, cTbT2, cTbT3,msq,born)
      
      call CSIeps1(cT1T2, cT1T3, cT2T3, cTaT1, cTaT2, cTaT3, 
     $     cTaTb, cTbT1, cTbT2, cTbT3,msq,
     $     pa,pb,p1,p2,p3,mu2,NF,result,div)

      

      end

      subroutine CSIeps1(cT1T2, cT1T3, cT2T3, cTaT1, cTaT2, cTaT3, 
     $     cTaTb, cTbT1, cTbT2, cTbT3,msq,
     $     pa,pb,p1,p2,p3,mu2,NF,result,div)
      implicit none
      
      double precision NF ! number of flavors
      double precision cT1T2, cT1T3, cT2T3, cTaT1, cTaT2, cTaT3, 
     $     cTaTb, cTbT1, cTbT2, cTbT3,msq

      

      integer div
      double precision pa(0:3),pb(0:3),p1(0:3),p2(0:3),p3(0:3)
      double precision mu2
      double complex result
    
      double precision CA,CF,gammaQ,gammaG,KG,KQ,Pi,TR
      double precision ga,gb,g1,g2,g3
      parameter(Pi=3.141592653589793d0)
      double precision sij
      double complex MyLog
      external MyLog
      external sij


      CA = 3.0d0
      CF = 4.0d0/3.0d0
    
      TR = 1.0d0/2.0d0

      gammaQ = 3.0d0/2.0d0*CF
      gammaG = 11.0d0/6.0d0*CA-2.0d0/3.0d0*TR*NF
C      KG = (67.0d0/18.0d0 - Pi**2/6.0d0)*CA-10.0d0/9.0d0*TR*NF
C      KQ = (7.0d0/2.0d0 - Pi**2/6.0d0)*CF


      ga = gammaQ
      gb = gammaQ
      g1 = gammaQ
      g2 = gammaQ
      g3 = gammaG

c      m1 = born(1)
c      m2 = born(2)
c      m3 = born(3)
c      m4 = born(4)

c      mu2 = mu2/2.0d0 ! rescale there is bug


c      print*,'div=',div
c      print*,'mu=',dsqrt(mu2)

c finite part
      if(div.eq.2) then
         result = 2.0d0*(cT1T2 + cT1T3 + cT2T3 + cTaT1 + cTaT2 + 
     $        cTaT3 + cTaTb + cTbT1 + cTbT2 + cTbT3)*dcmplx(1.0d0,0.0d0)
c         print*,"result(2)=",result
      elseif(div.eq.1) then
c         
c     change to the V convention 
c
         result =  ((-g1 - g2 - g3 - gA - gB)*msq**dcmplx(1.0d0,0.0d0) + 
     $        2.0d0*(cT1T2*MyLog(mu2/Sij(p1,p2)) + 
     $        cT1T3*MyLog(mu2/Sij(p1,p3)) + 
     $        cTaT1*MyLog(mu2/Sij(p1,pa)) + 
     $        cTbT1*MyLog(mu2/Sij(p1,pb)) + 
     $        cT2T3*MyLog(mu2/Sij(p2,p3)) + 
     $        cTaT2*MyLog(mu2/Sij(p2,pa)) + 
     $        cTbT2*MyLog(mu2/Sij(p2,pb)) + 
     $        cTaT3*MyLog(mu2/Sij(p3,pa)) + 
     $        cTbT3*MyLog(mu2/Sij(p3,pb)) + 
     $        cTaTb*MyLog(mu2/Sij(pa,pb))))

c         print*,"result(1)=",result
c         print*,"sij(pa,p3)/mu2=",sij(pa,p3)/mu2

      else
         result = 0.0d0
      endif

      end

      double precision function sij(a,b)
      implicit none 
      double precision a(0:3),b(0:3)!,dotrr
      !external dotrr
!not sure about this
      sij = -2.0d0*
     $     (a(0)*b(0) - 
     $     a(1)*b(1) - 
     $     a(2)*b(2) - 
     $     a(3)*b(3))

      return
      end

