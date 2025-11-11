c      subroutine qqhqqvbf(PBAR1,PBAR2,PBAR3,PBAR4,PBAR5,
c     $     PSIGN,CL1,CR1,CL3,CR3,
c     $     bosonMass, bosonwidth,IsVaw,divMax,GHVV,scale,born,virt)
      subroutine qqhqqvbf(PBAR1,PBAR2,PBAR3,PBAR4,PBAR5,
     $     PSIGN,ReCL1,ReCR1,ReCL3,ReCR3,ImCL1,ImCR1,ImCL3,ImCR3,
     $     bosonMass,bosonwidth,IsVaw,divMax,ReGHVV,ImGHVV,
     $     scale,born,virt)

c
c     q1(p1) + q3(p3) -> q2(p2) + q4(p4) + h(p5)
c
      implicit none
      logical lconvert,ldebug,ldebug1,ldebug2
      parameter(lconvert=.false.,ldebug=.false.,ldebug1=.false.)
      parameter(ldebug2=.false.)

      double precision PBAR1(0:3),PBAR2(0:3),PBAR3(0:3),PBAR4(0:3),
     $     PBAR5(0:3)
      double precision pa(0:3),pb(0:3),p1(0:3),p2(0:3)
      double precision born,virt(0:2),trivirt(0:3)
      double precision Ieps(0:2)
      double complex temp 
      double precision PBAR(0:3,5),P(0:3,5)
      integer isVaW(2)          ! 1 for W and 0 for Z
      double complex CL1(2),CR1(2),CL3(2),CR3(2)
      double precision ReCL1(2),ReCR1(2),ReCL3(2),ReCR3(2)
      double precision ImCL1(2),ImCR1(2),ImCL3(2),ImCR3(2)

      double complex CLR1(2,-1:1),CLR3(2,-1:1)

      
c     maybe set up couplings in a more general fashion
      double precision ReGHVV(2),ImGHVV(2),scale
      double complex GHVV(2)
      double precision bosonMass(2),bosonwidth(2)
      integer divMax
      integer isig1,isig3,div,iclass,mu,I,icol,isig2,isig4
      integer psign(5),fsign(4)
      double complex PSI(2,-1:1,4)
      double complex prop21,prop23,prop43,prop41
      double precision p21(0:4),p23(0:4),p41(0:4),p43(0:4)
      double complex M2_c(2)
      double precision massV2,widthV
      double complex VPROP
      external VPROP
      double complex j21(0:3,-1:1),j43(0:3,-1:1),
     $     j23(0:3,-1:1),j41(0:3,-1:1)
      double complex mat2143b(-1:1,-1:1),mat4123b(-1:1,-1:1)
      double complex tri21(0:2),tri23(0:3),
     $     tri41(0:3),tri43(0:3)
      double complex PENT2143(-1:1,-1:1,4,0:2)
      double complex PENT4123(-1:1,-1:1,4,0:2)
      double complex BORN2143(-1:1,-1:1,4)
      double complex BORN4123(-1:1,-1:1,4)
      complex*16 dotcc
      external dotcc
      double complex czero
      parameter (czero = (0.0d0,0.0d0))
      double complex jamp_virt(-1:1,-1:1,-1:1,-1:1,2,0:2),mv(2),mb(2)
      double complex jamp_born(-1:1,-1:1,-1:1,-1:1,2),
     $     amp_virt(-1:1,-1:1,-1:1,-1:1,4,0:2)

      double precision CA,CF,TR,NC,gammaQ
      parameter(CA = 3.0d0,NC=3.0d0,CF=4.0d0/3.0d0,TR=1.0d0/2.0d0)
      real*8 pi
      parameter(pi=3.14159265358979323846264338328d0)
      double precision maxCoup2143,maxCoup4123

c      divMax = 2
      
c fill coupling with real and imaginary parts
      do I=1,2
         CL3(I) = dcmplx(ReCL3(I),ImCL3(I))
         CL1(I) = dcmplx(ReCL1(I),ImCL1(I))

         CR3(I) = dcmplx(ReCR3(I),ImCR3(I))
         CR1(I) = dcmplx(ReCR1(I),ImCR1(I))
         
         GHVV(I) = dcmplx(ReGHVV(I),ImGHVV(I))
      enddo

c      print*,"test=",abs(ghvv(1)*ghvv(2))

      do I=1,2
         CLR3(I,-1) = CL3(I)
         CLR3(I,1) = CR3(I)

         CLR1(I,-1) = CL1(I)
         CLR1(I,1) = CR1(I)
      enddo
      if(ldebug2) then
         print*,"-----" 
         print*,"CL1(1),CR1(1),bosonMass(1)=",cl1(1),cr1(1),bosonMass(1)
         print*,"CL3(1),CR3(1)=",cl3(1),cr3(1)
         
         print*,"CL1(2),CR1(2),bosonMass(2)=",cl1(2),cr1(2),bosonMass(2)
         print*,"CL3(2),CR3(2)=",cl3(2),cr3(2)
         print*,"-----"
      endif
      do isig1=-1,1
         do isig3=-1,1
            mat2143b(isig1,isig3) = czero
            mat4123b(isig1,isig3) = czero
            do div=0,2
               do iclass=1,4
                  PENT2143(isig1,isig3,iclass,div)=czero
                  PENT4123(isig1,isig3,iclass,div)=czero
               enddo
            enddo
         enddo
      enddo
c
c     Define diagrammatic momenta
c     HZ style
c     
      do mu=0,3
         PBAR(mu,1) = PBAR1(mu)
         PBAR(mu,2) = PBAR2(mu)
         PBAR(mu,3) = PBAR3(mu)
         PBAR(mu,4) = PBAR4(mu)
         PBAR(mu,5) = PBAR5(mu) ! higgs

         do I=1,5
            P(mu,I) = dble(psign(I))*PBAR(mu,I)
         enddo
      enddo
         
      do I=1,4 ! fermions
         fsign(I)=psign(I)
      enddo
      
c
c     set up born matrix element
c
c     build 2-spinors
c
      call PSI0M(4,PBAR,FSIGN,PSI)
c
c     construct j43,j21 currents
c     
      call curr(1,psi(1,-1,2),psi(1,-1,1),j21)
      call curr(1,psi(1,-1,4),psi(1,-1,3),j43)

c      
c     construct j41,j23 currents
c     
      call curr(1,psi(1,-1,2),psi(1,-1,3),j23)
      call curr(1,psi(1,-1,4),psi(1,-1,1),j41)
c      
c     build propagators
c
      do mu=0,3
         p21(mu) = p(mu,2) - p(mu,1)
         p43(mu) = p(mu,4) - p(mu,3)
      enddo

      p21(4) = p21(0)**2-p21(1)**2-p21(2)**2-p21(3)**2
      p43(4) = p43(0)**2-p43(1)**2-p43(2)**2-p43(3)**2
      
c      
c     Build weak boson propagators
c
c     
      massV2 = bosonMass(1)**2
      widthV = bosonWidth(1) ! add to arguments
      
      prop21 = VPROP(p21(4),massV2,widthV)
      prop43  = VPROP(p43(4),massV2,widthV)

c            
c     build propagators
c
      do mu=0,3
         p41(mu) = p(mu,4) - p(mu,1)
         p23(mu) = p(mu,2) - p(mu,3)
      enddo

      p41(4) = p41(0)**2-p41(1)**2-p41(2)**2-p41(3)**2
      p23(4) = p23(0)**2-p23(1)**2-p23(2)**2-p23(3)**2
      
c      
c     Build weak boson propagators
c
c     
      massV2 = bosonMass(2)**2
      widthV = bosonWidth(2) ! add to arguments
      
      prop41 = VPROP(p41(4),massV2,widthV)
      prop23  = VPROP(p23(4),massV2,widthV)
c
c     
c
      do isig1=-1,1,2
         do isig3=-1,1,2
            mat2143b(isig1,isig3) = dotcc(j21(0,isig1),j43(0,isig3))*
     $           prop21*prop43

            mat4123b(isig1,isig3) = dotcc(j41(0,isig1),j23(0,isig3))*
     $           prop41*prop23
            
         enddo
      enddo
c
c     QCD triangles
c
      call VBF_TRI(tri43,p43(4),scale)
      call VBF_TRI(tri21,p21(4),scale)

      call VBF_TRI(tri41,p41(4),scale)
      call VBF_TRI(tri23,p23(4),scale)
      
c    
c QCD pentagons
c

      M2_c(1)=dcmplx(bosonMass(1)**2,-bosonMass(1)*bosonWidth(1))
      M2_c(2)=dcmplx(bosonMass(2)**2,-bosonMass(2)*bosonWidth(2))
c      print*,"bosonMass=",bosonMass
      
c For speeding up code. We do not need pentagons if the couplings are zero.
c If both ghvv(1) and ghvv(2) is non-zero we compute the pentagons.


      if(abs(GHVV(1))*abs(GHVV(2)).gt.0d0) then
         call PENT_VBF(PENT2143,BORN2143,M2_c(1),
     $        p(0,1),PSI(1,-1,1),p(0,2),PSI(1,-1,2),
     $        p(0,3),PSI(1,-1,3),p(0,4),PSI(1,-1,4),
     $        p(0,5),scale,IsVaW(1),divmax)
      endif
c     
c     
c

      if(ldebug1) then
      do isig1=-1,1,2
         do isig3=-1,1,2
            print*,"isig1=",isig1,"isig3=",isig3
            do iclass=1,4
               print*,"from PENT_VBF:born2143(",iclass,")=",
     $              born2143(isig1,isig3,iclass)
            enddo
            print*,"born2143=",mat2143b(isig1,isig3)
         enddo
      enddo
      endif

      if(abs(GHVV(2))*abs(GHVV(1)).gt.0d0) then
         call PENT_VBF(PENT4123,BORN4123,M2_c(2),
     $        p(0,1),PSI(1,-1,1),p(0,4),PSI(1,-1,4),
     $        p(0,3),PSI(1,-1,3),p(0,2),PSI(1,-1,2),
     $        p(0,5),scale,IsVaW(2),divmax)
      endif
c
c     
c
      if(ldebug1) then
      do isig1=-1,1,2
         do isig3=-1,1,2
            print*,"isig1=",isig1,"isig3=",isig3
            do iclass=1,4
               print*,"from PENT_VBF:born4123(",iclass,")=",
     $              born4123(isig1,isig3,iclass)
            enddo
            print*,"born4123=",mat4123b(isig1,isig3)
         enddo
      enddo
      endif
c
c     initialize jamp_virt and jamp_born
c
c
      do isig1=-1,1
         do isig2=-1,1
            do isig3=-1,1
               do isig4=-1,1
                  do icol=1,2
                     jamp_born(isig1,isig2,isig3,isig4,icol)= czero
                     do div=0,2
                        jamp_virt(isig1,isig2,isig3,isig4,icol,div)= 
     $                       czero
                     enddo
                  enddo
                  do icol=1,4
                     do div=0,2
                        amp_virt(isig1,isig2,isig3,isig4,icol,div)= 
     $                       czero
                     enddo
                  enddo
                  
               enddo
            enddo
         enddo
      enddo
c
c     born color sub-amplitudes
c
      do isig1=-1,1,2
         do isig3=-1,1,2
            if(abs(CLR3(1,isig3)*CLR1(1,isig1)).gt.0d0) then
               isig2=isig1
               isig4=isig3
               
               jamp_born(isig1,isig2,isig3,isig4,1) = 
     $              GHVV(1)*CLR3(1,isig3)*CLR1(1,isig1)*
     $              mat2143b(isig1,isig3)

            else
                isig2=isig1
                isig4=isig3
               
               jamp_born(isig1,isig2,isig3,isig4,1) = czero

            endif

            if(abs(CLR3(2,isig3)*CLR1(2,isig1)).gt.0d0) then
               isig2=isig3 ! 
               isig4=isig1 !
c     
               jamp_born(isig1,isig2,isig3,isig4,2) = 
     $              GHVV(2)*CLR3(2,isig3)*CLR1(2,isig1)*
     $              mat4123b(isig1,isig3)
            else
               isig2=isig3 !
               isig4=isig1 !
c     
               jamp_born(isig1,isig2,isig3,isig4,2) = czero

            endif

         enddo
      enddo
c
c     Virtual color sub-amplitudes 
c
c      
c     NC = 3 
c     CF 
c     We have used the Fierz identity on Gell-Mann matrices. 
      do isig1=-1,1,2
         do isig3=-1,1,2
            if(abs(CLR3(1,isig3)*CLR1(1,isig1)).gt.0d0) then
               isig2=isig1
               isig4=isig3
c
               do div=0,divmax
                  amp_virt(isig1,isig2,isig3,isig4,1,div) = ! delta_i2 i1 delta_i4 i3
     $                 GHVV(1)*CLR3(1,isig3)*CLR1(1,isig1)*
     $                 mat2143b(isig1,isig3)*CF*(tri21(div)+tri43(div))
                  
                  amp_virt(isig1,isig2,isig3,isig4,3,div) = ! t^{a}_{i2 i1} t^{a}_{i4 i3}
     $                 GHVV(1)*CLR3(1,isig3)*CLR1(1,isig1)*
     $                 (PENT2143(isig1,isig3,1,div)+
     $                 PENT2143(isig1,isig3,2,div)+
     $                 PENT2143(isig1,isig3,3,div)+
     $                 PENT2143(isig1,isig3,4,div)) 

                  
               enddo

            else
               isig2=isig1
               isig4=isig3
               
               do div=0,divmax
                  amp_virt(isig1,isig2,isig3,isig4,1,div) = czero

                  amp_virt(isig1,isig2,isig3,isig4,3,div) = czero
               enddo

               
            endif
            
            if(abs(CLR3(2,isig3)*CLR1(2,isig1)).gt.0d0) then
               isig2=isig3
               isig4=isig1
               do div=0,divmax
                  amp_virt(isig1,isig2,isig3,isig4,2,div) = ! -delta_{i2 i3} delta_{i4 i1}
     $                 GHVV(2)*CLR3(2,isig3)*CLR1(2,isig1)*
     $                 mat4123b(isig1,isig3)*CF*(tri41(div)+tri23(div))

                  amp_virt(isig1,isig2,isig3,isig4,4,div) = ! -t^{a}_{i2 i3} t^{a}_{i4 i1}
     $                  GHVV(2)*CLR3(2,isig3)*CLR1(2,isig1)*
     $                 (PENT4123(isig1,isig3,1,div)+
     $                 PENT4123(isig1,isig3,2,div)+
     $                 PENT4123(isig1,isig3,3,div)+
     $                 PENT4123(isig1,isig3,4,div))         
                  
               enddo

            else
               isig4=isig1
               isig2=isig3
               
               do div=0,divmax
                  amp_virt(isig1,isig2,isig3,isig4,2,div) = czero

                  amp_virt(isig1,isig2,isig3,isig4,4,div) = czero
               enddo

            endif
            
         enddo
      enddo
c
c     now use fierz identity 
c

      do div=0,divmax
         do isig1=-1,1,2
            do isig2=-1,1,2
               do isig3=-1,1,2
                  do isig4=-1,1,2
                     
                     jamp_virt(isig1,isig2,isig3,isig4,1,div) = ! detla_i2_i1 delta_i4_i3 
     $                    amp_virt(isig1,isig2,isig3,isig4,1,div) - 
     $                    0.5d0/NC *  
     $                    amp_virt(isig1,isig2,isig3,isig4,3,div) -
     $                    0.5d0 * amp_virt(isig1,isig2,isig3,isig4,4,div)
                     
                     
                     jamp_virt(isig1,isig2,isig3,isig4,2,div) = ! -delta_i4_i1 delta_i2_i3
     $                    amp_virt(isig1,isig2,isig3,isig4,2,div) - 
     $                    0.5d0/NC *  
     $                    amp_virt(isig1,isig2,isig3,isig4,4,div) -
     $                    0.5d0 * amp_virt(isig1,isig2,isig3,isig4,3,div)
                     
                     
                     
                  enddo
               enddo
            enddo
         enddo
      enddo
c
c
c     Compute |M_B|^2 summed over colors and helicities 
c     Compute 2 Re[M_V^{\dagger} M_B] summed over colors and helicities
c
      born = 0.0d0

c       isig1=1
c       isig2=1
c       isig3=1
c       isig4=1

      do isig1=-1,1,2
         do isig2=-1,1,2
            do isig3=-1,1,2
               do isig4=-1,1,2
                  
                  mb(1) = jamp_born(isig1,isig2,isig3,isig4,1)
                  mb(2) = jamp_born(isig1,isig2,isig3,isig4,2)

                  born = born + NC*NC*(dreal(mb(1))**2+dimag(mb(1))**2 + 
     $                  (dreal(mb(2))**2+dimag(mb(2))**2)) -  
     $                  NC*2.0d0*dreal(dconjg(mb(1))*mb(2))
                
                  
               enddo
            enddo
         enddo
      enddo

      


      do div=0,divmax
         virt(div) = 0.0d0
         trivirt(div)=0.0d0
         
c         print*,"isig",isig1,isig2,isig3,isig4
         do isig1=-1,1,2
            do isig2=-1,1,2
               do isig3=-1,1,2
                  do isig4=-1,1,2
                     ! temp variables 
                     mb(1) = jamp_born(isig1,isig2,isig3,isig4,1)
                     mb(2) = jamp_born(isig1,isig2,isig3,isig4,2)

                     mv(1) = jamp_virt(isig1,isig2,isig3,isig4,1,div)/
     $                    (4.0d0*pi)
                     mv(2) = jamp_virt(isig1,isig2,isig3,isig4,2,div)/
     $                    (4.0d0*pi)
                     
                     virt(div) = virt(div) +2.0d0*dreal(dconjg(mb(1))*
     $                    (mv(1)*NC*NC-mv(2)*NC) +dconjg(mb(2))*
     $                    (mv(2)*NC*NC-mv(1)*NC))
                         
                     

c                     virt(div) = virt(div) + 2.0d0*dreal(NC*NC*(dconjg(mv(1))*mb(1) + 
c     $                    dconjg(mv(2))*mb(2)) - NC*(dconjg(mv(1))*mb(2) + 
c     $                    dconjg(mv(2))*mb(1)))  

                     mv(1) = amp_virt(isig1,isig2,isig3,isig4,1,div)/
     $                    (4.0d0*pi) ! only triangles
                     mv(2) = amp_virt(isig1,isig2,isig3,isig4,2,div)/
     $                    (4.0d0*pi) ! only triangles


                     trivirt(div) = trivirt(div) +2.0d0*dreal(dconjg(mb(1))*
     $                    (mv(1)*NC*NC-mv(2)*NC) +dconjg(mb(2))*
     $                    (mv(2)*NC*NC-mv(1)*NC))

c                     trivirt(div) = trivirt(div) + 2.0d0*dreal(NC*NC*(dconjg(mv(1))*mb(1) + 
c     $                    dconjg(mv(2))*mb(2)) - NC*(dconjg(mv(1))*mb(2) + 
c     $                    dconjg(mv(2))*mb(1)))  
c                     
                  enddo
               enddo
            enddo
         enddo
      enddo
      
      call ConvertToCDRHjj(virt,virt,born,1.0d0)

c
c     check IR poles
c
      do mu=0,3
         pa(mu)=-1.0d0*p(mu,1)! quark
         pb(mu)=-1.0d0*p(mu,3)! quark
         p1(mu)=p(mu,2)! quark
         p2(mu)=p(mu,4)! quark
      enddo

      do div=0,divmax
         Ieps(div) = 0.0d0
         do isig1=-1,1,2
            do isig2=-1,1,2
               do isig3=-1,1,2
                  do isig4=-1,1,2
                     mb(1) = jamp_born(isig1,isig2,isig3,isig4,1)
                     mb(2) = jamp_born(isig1,isig2,isig3,isig4,2)

                     if(lconvert) then

                        call ConvertVirtHjj(pa,pb,p1,p2,mb,scale,
     $                       temp,div)
                        virt(div) = virt(div) + dreal(temp)/(2.0d0*pi)
                     endif
                     
                     
                     call CSIepsHjj(pa,pb,p1,p2,mb,scale,temp,div)
                     
                     Ieps(div) = Ieps(div) + dreal(temp)/(2.0d0*pi)
                     
                  enddo 
               enddo
            enddo
         enddo
         
      enddo
      
      if(ldebug) then
      gammaQ = 3.0d0/2.0d0*CF
      print*,"virt(1)/(-sum gamma_i)",virt(1)/
     $     (-4.0d0*gammaQ)/born*2.0d0*pi



      do div=0,2
         if(div.eq.0) then 
            print*,"finite terms"
            print*,"Virt(0)/born=",
     $        virt(div)/born,"Tri(0)/born=",
     $        trivirt(div)/born,"born=",born

            print*,"Virt: percent diff",100.0d0*(trivirt(div)-virt(div))/virt(div)
            
         else
         print*,"Ieps(",div,")/born",Ieps(div)/born,"Virt/born=",
     $        virt(div)/born,"Virt/Ieps=",
     $        virt(div)/Ieps(div),"born=",born

          print*,"Tri/Ieps=",
     $        trivirt(div)/Ieps(div),"born=",born
         endif
      enddo

      print*,"ratio=",Ieps(2)/virt(2)
      print*,"ratio=",Ieps(1)/virt(1)

      endif
c reweight output

      do div=0,2
         virt(div) = virt(div)*(2.0d0*pi)
      enddo

      end
c
c
c           Convert result from DR to CDR.
c
      subroutine ConvertToCDRHjj(virtCDR,virtDR,born,alphas)
      implicit none
      double precision alphas,virtCDR(0:2),virtDR(0:2),born
      double precision CF,CA,PI
      parameter(PI=3.14159265358979323846264338328d0)


      CA = 3.0d0
      CF = 4.0d0/3.0d0

      virtCDR(2) = virtDR(2)
      virtCDR(1) = virtDR(1)
      virtCDR(0) = virtDR(0) - alphas/(2.0d0*pi)*
     $     born*(2.0d0*CF)

      end
      
c
c
      subroutine PENT_VBF(PENT,BORN,M,p1,PSI1,p2,PSI2,p3,PSI3,
     $     p4,PSI4,
     $     p5,musq,IsVaWtemp,divmax)
      
      
C     
c     p1(mu) PSI1 q1 
c     p2(mu) PSI2 q2bar
c     p3(mu) PSI3 q3
c     p4(mu) PSI4 q4bar
c     p5(mu) higgs

      implicit none
      integer IsVaWtemp
      integer divmax
      double complex Pent(-1:1,-1:1,4,0:2)
  
      integer div,comp,mu
      double precision musq
      double complex M 

      double complex psi_k5(4),psi_k3(4),barpsi_k2(4),barpsi_k1(4)
      double precision k1(0:3),k2(0:3),k3(0:3),k4(0:3),k5(0:3)
      double precision p1(0:3),p2(0:3),p3(0:3),p4(0:3),p5(0:3)
      double complex PSI1(2,-1:1),PSI2(2,-1:1),PSI3(2,-1:1),PSI4(2,-1:1) 
      double complex Born(-1:1,-1:1,4)
   
 

C     define inflowing momentum (k1+k2+ ... + k5 = 0)
C     interface to Paco's code
      do mu=0,3
         k1(mu) = -p2(mu)! psibar 
         k2(mu) = -p4(mu)! psibar
         k3(mu) = p3(mu) ! psi
         k4(mu) = -p5(mu)! higgs
         k5(mu) = p1(mu) ! psi
      enddo

C
C calculate bar and ket or psi and psibar spinors
C psi is a 4-spinor which will be constructed frow two 2-spinors.
C psibar is a 4-spinor which will be constructed from two 2-spinors.
C
C     define psibar_k1,psibar_k2, psi_k3,psi_k5
C
C build 4-spinors
      barpsi_k1(1) = PSI2(1,1)!P
      barpsi_k1(2) = PSI2(2,1)!P
      barpsi_k1(3) = PSI2(1,-1)!M
      barpsi_k1(4) = PSI2(2,-1)!M


      barpsi_k2(1) = PSI4(1,1)!P
      barpsi_k2(2) = PSI4(2,1)!P
      barpsi_k2(3) = PSI4(1,-1)!M
      barpsi_k2(4) = PSI4(2,-1)!M


      psi_k5(1) = PSI1(1,-1)!M
      psi_k5(2) = PSI1(2,-1)!M
      psi_k5(3) = PSI1(1,1)!P
      psi_k5(4) = PSI1(2,1)!P


      psi_k3(1) = PSI3(1,-1)!M
      psi_k3(2) = PSI3(2,-1)!M
      psi_k3(3) = PSI3(1,1)!P
      psi_k3(4) = PSI3(2,1)!P
      

      do div=0,divmax
c FC only computed once for each divergences
            comp=1
c Fc finish
c            comp=1
            call Hjj77T_c(M,k1,k2,k3,k4,k5,barpsi_k1,psi_k5,barpsi_k2,
     $           psi_k3,musq,comp,
     $           Pent(-1,-1,1,div),
     $           Born(-1,-1,1),div)
            comp=-1
c Fix the rest now.


            if(IsVaWtemp.eq.0) then  

               call Hjj67T_c(M,k1,k2,k3,k4,k5,barpsi_k1,psi_k5,barpsi_k2,
     $           psi_k3,musq,comp,
     $           Pent(1,-1,1,div),
     $           Born(1,-1,1),div)

               call Hjj76T_c(M,k1,k2,k3,k4,k5,barpsi_k1,psi_k5,barpsi_k2,
     $           psi_k3,musq,comp,
     $           Pent(-1,1,1,div),
     $           Born(-1,1,1),div)

               call Hjj66T_c(M,k1,k2,k3,k4,k5,barpsi_k1,psi_k5,barpsi_k2,
     $           psi_k3,musq,comp,
     $           Pent(1,1,1,div),
     $           Born(1,1,1),div)
            
            
            endif
          
       enddo


         
c Cross 

       do div=0,divmax
c     FC only computed once for each divergences
          comp=1
c     Fc finish
c     comp=1
          call HjjCross77T_c(M,k5,k2,k3,k4,k1,psi_k5,barpsi_k1,barpsi_k2,
     $         psi_k3,musq,comp,
     $         Pent(-1,-1,2,div),
     $         Born(-1,-1,2),div)
          comp=-1
          
          
            if(IsVaWtemp.eq.0) then  
               
               call HjjCross67T_c(M,k5,k2,k3,k4,k1,psi_k5,barpsi_k1,barpsi_k2,
     $              psi_k3,musq,comp,
     $              Pent(1,-1,2,div),
     $              Born(1,-1,2),div)
               
               call HjjCross76T_c(M,k5,k2,k3,k4,k1,psi_k5,barpsi_k1,barpsi_k2,
     $           psi_k3,musq,comp,
     $              Pent(-1,1,2,div),
     $              Born(-1,1,2),div)
               
               call HjjCross66T_c(M,k5,k2,k3,k4,k1,psi_k5,barpsi_k1,barpsi_k2,
     $              psi_k3,musq,comp,
     $           Pent(1,1,2,div),
     $              Born(1,1,2),div)
               
            endif
         enddo
         
c CrossF
         
         do div=0,divmax
c FC only computed once for each divergences
            comp=1
c Fc finish
c            comp=1
            call HjjCrossF77T_c(M,k1,k3,k2,k4,k5,barpsi_k1,psi_k5,psi_k3,
     $           barpsi_k2,musq,comp,
     $           Pent(-1,-1,3,div),
     $           Born(-1,-1,3),div)
            comp=-1


            if(IsVaWtemp.eq.0) then  

              call HjjCrossF67T_c(M,k1,k3,k2,k4,k5,barpsi_k1,psi_k5,psi_k3,
     $           barpsi_k2,musq,comp,
     $           Pent(1,-1,3,div),
     $           Born(1,-1,3),div)

              call HjjCrossF76T_c(M,k1,k3,k2,k4,k5,barpsi_k1,psi_k5,psi_k3,
     $           barpsi_k2,musq,comp,
     $           Pent(-1,1,3,div),
     $           Born(-1,1,3),div)

              call HjjCrossF66T_c(M,k1,k3,k2,k4,k5,barpsi_k1,psi_k5,psi_k3,
     $           barpsi_k2,musq,comp,
     $           Pent(1,1,3,div),
     $           Born(1,1,3),div)
            
            
            endif
          
       enddo

c CrossIF

        do div=0,divmax
c FC only computed once for each divergences
            comp=1
c Fc finish
c            comp=1
            call HjjCrossIF77T_c(M,k5,k3,k2,k4,k1,psi_k5,barpsi_k1,psi_k3,
     $           barpsi_k2,musq,comp,
     $           Pent(-1,-1,4,div),
     $           Born(-1,-1,4),div)
            comp=-1


            if(IsVaWtemp.eq.0) then  

              call HjjCrossIF67T_c(M,k5,k3,k2,k4,k1,psi_k5,barpsi_k1,psi_k3,
     $           barpsi_k2,musq,comp,
     $           Pent(1,-1,4,div),
     $           Born(1,-1,4),div)

              call HjjCrossIF76T_c(M,k5,k3,k2,k4,k1,psi_k5,barpsi_k1,psi_k3,
     $           barpsi_k2,musq,comp,
     $           Pent(-1,1,4,div),
     $           Born(-1,1,4),div)

              call HjjCrossIF66T_c(M,k5,k3,k2,k4,k1,psi_k5,barpsi_k1,psi_k3,
     $           barpsi_k2,musq,comp,
     $           Pent(1,1,4,div),
     $           Born(1,1,4),div)
            
            
            endif
          
       enddo


         
      end

