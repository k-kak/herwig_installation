C     This subroutine shall compute all types of VBF graphs.
C
C     
C     virtual matrix element for q1(p1) + q3(p3) -> q2(p2) +  q4(p4) + g(p5) + H(p6)
c                               
C     neutral currents: g[Zccbar]_{L,R}=g[Zuubar]_{L,R}
C
C     helicity
C     define a color sub amplitude
C
C      This is an alternate way of calling the amplitudes

c      subroutine qqhqqgvbf(PBAR1,PBAR2,PBAR3,PBAR4,PBAR5,PBAR6,
c     $     PSIGN,CL1,CR1,CL3,CR3,
c     $     Mass,Width,IsVaW,divMax,gaugegood,gaugebad,GHVV,scale,NF,
c     $     born,virt)
      subroutine qqhqqgvbf(PBAR1,PBAR2,PBAR3,PBAR4,PBAR5,PBAR6,
     $     PSIGN,ReCL1,ReCR1,ReCL3,ReCR3,ImCL1,ImCR1,ImCL3,ImCR3,
     $     Mass,Width,IsVaW,divMax,gaugegood,gaugebad,gaugethres,
     $     gaugeacc,
     $     ReGHVV,ImGHVV,
     $     scale,NF,
     $     born,virt)
      implicit none

      double precision gaugethres
      double precision PBAR1(0:3),PBAR2(0:3),PBAR3(0:3),PBAR4(0:3),
     $     PBAR5(0:3),PBAR6(0:3),PBAR(0:3,6)
      double precision ReCL1(2),ReCR1(2)! left-right couplings
      double precision ReCL3(2),ReCR3(2) 
      double precision ImCL1(2),ImCR1(2) ! left-right couplings
      double precision ImCL3(2),ImCR3(2)

      double complex  CL1(2),CR1(2) ! left-right couplings
      double complex CL3(2),CR3(2) 

      integer PSIGN(6),FLAV(4)
      double precision born,virt(0:2),scale

      double complex CLR(4,2,-1:1)

      double precision Width(2),Mass(2)
c      double precision finvirt,eps1virt,eps2virt
      double precision ReGHVV(2),ImGHVV(2)
      double complex GHVV(2)
c      double precision GHVV1,GHVV2
      double precision NF
      integer i
      integer IsVaW(2)! 1 for W and 0 for Z
      integer divMax!0 finite,1 finte+1/Eps, 2 finite+1/eps+1/eps^2
      integer Gaugegood(3)
      integer Gaugebad(3)
      integer Gaugetemp(3,2)
      DOUBLE PRECISION gaugeacc
      FLAV(1)=1
      FLAV(2)=2 ! dummy
      FLAV(3)=3
      FLAV(4)=4 ! dummy
c
c
c     We assume their are not flavor changing processes.
c fill coupling with real and imaginary parts
      do I=1,2
         CL3(I) = dcmplx(ReCL3(I),ImCL3(I))
         CL1(I) = dcmplx(ReCL1(I),ImCL1(I))
         CR3(I) = dcmplx(ReCR3(I),ImCR3(I))
         CR1(I) = dcmplx(ReCR1(I),ImCR1(I))

         GHVV(I) = dcmplx(ReGHVV(I),ImGHVV(I))
      enddo
c
c q1 --------- CL1,CR1 ---- 
c            !
c            ! V
c
c
c
c     Left/right couplings 
c     of q1 -> V q.

      CLR(1,1,-1) = CL1(1)
      CLR(1,2,-1) = CL1(2)
      CLR(1,1,1) = CR1(1)
      CLR(1,2,1) = CR1(2)

c     Left/right couplings 
c     of q3 -> V q.

      CLR(3,1,-1) = CL3(1)
      CLR(3,2,-1) = CL3(2)
      CLR(3,1,1) = CR3(1)
      CLR(3,2,1) = CR3(2)

      

      do i=0,3
         PBAR(i,1) = PBAR1(i)
         PBAR(i,2) = PBAR2(i)
         PBAR(i,3) = PBAR3(i)
         PBAR(i,4) = PBAR4(i)
         PBAR(i,5) = PBAR5(i) ! gluon
         PBAR(i,6) = PBAR6(i) ! higgs
      enddo
      

      call matqqhqqgVBF(.false.,PBAR,PSIGN,FLAV,Mass,Width,IsVaW,divMax,
     $   gaugetemp,gaugethres,gaugeacc,CLR,GHVV,scale,NF,born,virt)
      

      do i=1,3
         gaugegood(i)=gaugetemp(i,1)
         gaugebad(i)=gaugetemp(i,2)
      enddo
      

      end
C
      subroutine matqqhqqgVBF(lcheckgauge,PBar,psign,flavor,bosonMass,
     $     bosonWidth,IsVaW,divMax,gaugebad,gaugethres1,GAUGEACCURACY1,
     $     CLR,GHVV,scale,NF,
     $     born,virt)
      implicit none 
      integer IsVaW(2) ! 1 for W and 0 for Z
      integer divMax!0 finite,1 finte+1/Eps, 2 finite+1/eps+1/eps^2
      logical lcheckgauge,ldebug
      parameter(ldebug=.false.)
      double precision scale, alphas ! muR^2
      integer flavor(4)! flavor of quarks 1--4
      double precision bosonMass(2) ! mass of exchange weak boson in direct and exchange diagrams
      double complex GHVV(2) ! coupling of Higgs to gauge bosons for set 1 and set 2 Feynman diagrams


      double precision ColorSumLO, ColorSumNLO
      double complex VPROP
      external ColorSumLO
      external ColorSumNLO
      external VPROP

      logical compPentagons,compBoxes,lconvert,ldebug1
      parameter(ldebug1=.false.)
      parameter(compPentagons=.true.,compBoxes=.true.,lconvert=.false.)

      double complex CLR(4,2,-1:1) ! CLR(flavor index, gauge boson index, helicity index) 

      DOUBLE COMPLEX  PSI(2,-1:1,4)
C
      double complex jamp_virt(-1:1,-1:1,-1:1,-1:1,1:2,1:14,0:2)

      double complex amp_virt(-1:1,-1:1,-1:1,-1:1,1:2,1:4,0:2)
      double complex jamp_born(-1:1,-1:1,-1:1,-1:1,1:2,1:4)

      double complex prop21,prop21g,prop23,prop23g,prop43,prop43g,prop41
      double complex prop41g

      double complex test,test1
      double precision Ieps(0:2)
      double complex temp
      double precision pa(0:3),pb(0:3)

      double complex j21(0:3,-1:1),j43(0:3,-1:1),
     $     j23(0:3,-1:1),j41(0:3,-1:1)

      double complex trifac21(0:2),trifac23(0:3),
     $     trifac41(0:3),trifac43(0:3)

      double complex tri21(0:2),tri23(0:3),
     $     tri41(0:3),tri43(0:3)


      double complex M2_c(2)

      double complex mat43(5,-1:1,-1:1,2,3,0:2), 
     $     mat43g(2,-1:1,-1:1,2,3,0:2), 
     $     mat43gb(2,-1:1,-1:1,2,3,0:2), 
     $     mat43b(-1:1,-1:1,2,3,0:2)
c
      double complex mat21(5,-1:1,-1:1,2,3,0:2), 
     $     mat21g(2,-1:1,-1:1,2,3,0:2), 
     $     mat21gb(2,-1:1,-1:1,2,3,0:2), 
     $     mat21b(-1:1,-1:1,2,3,0:2)

      double complex mat41(5,-1:1,-1:1,2,3,0:2), 
     $     mat41g(2,-1:1,-1:1,2,3,0:2), 
     $     mat41gb(2,-1:1,-1:1,2,3,0:2), 
     $     mat41b(-1:1,-1:1,2,3,0:2)
c
      double complex mat23(5,-1:1,-1:1,2,3,0:2), 
     $     mat23g(2,-1:1,-1:1,2,3,0:2), 
     $     mat23gb(2,-1:1,-1:1,2,3,0:2), 
     $     mat23b(-1:1,-1:1,2,3,0:2)

      double complex NOABEmat21(-1:1,-1:1,2,0:2), 
     $     NOABEmat41(-1:1,-1:1,2,0:2), NOABEmat23(-1:1,-1:1,2,0:2), 
     $     NOABEmat43(-1:1,-1:1,2,0:2) 
       double complex NOABEmat21g(2,-1:1,-1:1,2,0:2), 
     $     NOABEmat41g(2,-1:1,-1:1,2,0:2),
     $     NOABEmat23g(2,-1:1,-1:1,2,0:2), 
     $     NOABEmat43g(2,-1:1,-1:1,2,0:2)

      double complex eps(0:3,2)
      double precision epsR(0:3,2)
      integer isig,j
      double precision q(0:4),pk(0:4,4)
      complex*16 braket(2,-1:1,4,2), 
     1        jh1(0:3,-1:1), jh2(0:3,-1:1), e21(0:3,-1:1,2),
     2        e43(0:3,-1:1,2), e23(0:3,-1:1,2),e41(0:3,-1:1,2)

      double complex PentHex21(3,-1:1,-1:1,2,4,0:2)
      double complex PentHex43(3,-1:1,-1:1,2,4,0:2)
      double complex PentHex23(3,-1:1,-1:1,2,4,0:2)
      double complex PentHex41(3,-1:1,-1:1,2,4,0:2)

      double complex born21(-1:1,-1:1,2,4),born43(-1:1,-1:1,2,4)
      double complex born41(-1:1,-1:1,2,4),born23(-1:1,-1:1,2,4)

      double complex HEX2143(2,-1:1,-1:1,2,4,0:2)
      double complex HEX4123(2,-1:1,-1:1,2,4,0:2)
      double complex BORN2143(-1:1,-1:1,2,4)
      double complex BORN4123(-1:1,-1:1,2,4)
c
c     mat43g(1,isig1,isig3,gluon_hel,1,div)
c     array store amplitude piece for a give helicity
      
      integer psign(6),FSIGN(4),mu,I,gluon_hel
      double precision PBAR(0:3,6)
      double precision P(0:3,6)
      double precision p1(0:3),p2(0:3),p3(0:3),p4(0:3),pg(0:3)
      double precision p21g(0:4),p23g(0:4),p23(0:4),p43(0:4),p21(0:4)
      double precision p43g(0:4),p41g(0:4),p41(0:4)
c
      double precision virt(0:2),born,boxvirt(0:2)
      double precision massV2,widthV
      double precision bosonWidth(2)
      double complex matLO(4),matNLO(4)
      integer icol
c
cFC   integer variable control the repetition of
cFC  calculation of the virtual integrals for the 
cFC  same momenta
      integer comp,div,isig1,isig3,isig2,isig4,iclass
      double complex Ic,czero
      parameter(Ic=(0.0d0,1.0d0),czero=(0.0D0,0.0D0))
      double precision SgnF,SgnF2
      parameter (SgnF=-1.0) ! Sign factor due to exchange of external fermion lines
      
cFC
      double precision fac1,fac2
      double precision CA,CF
      parameter(CA = 3.0d0,CF=4.0d0/3.0d0)

c      parameter(CA = 0.0d0,CF=4.0d0/3.0d0)

      complex*16 dotcc
      external dotcc
c FC
      real*8 pi
      parameter(pi=3.14159265358979323846264338328d0)
      real*8 ptemp12(0:3),ptemp12sq
      complex*16 ctetemp(0:2),lns
      real*8 theta
      external theta
      logical VER
      parameter(VER=.false.)
      logical Wardbox1,Wardbox2,wardbox3,wardbox4
      logical Wardhex1,wardhex2,wardhex3,wardhex4
      logical WardhexN1, wardhexN2
      integer Gaugebad(3,2) ! 1 box, 2 hex, 3, hexNo Abe, 
                            ! second argument 1 good 2 bad
      integer j1,j2
      double precision nf
      double precision Gaugethres,gaugethres1,gaugeAccuracy(2)
      double precision GAUGEACCURACY1
      common/Gauge/gaugethres,gaugeAccuracy
      gaugethres=gaugethres1
      gaugeAccuracy(1)=0d0
      gaugeAccuracy(2)=0d0
      GAUGEACCURACY1=0D0
      do j1=1,3
         do j2=1,2
            Gaugebad(j1,j2)=0
         enddo
      enddo      

cFC

      SgnF2 = 1.0d0

c      CF = 4.0d0/3.0d0

c     To begin we will take the color sum and sum over helicities 
c     virt(i) = 2 Re[born^* virt]
c     virt(2) is 1/eps^2,virt(1) is 1/esp^1 and virt(0) is the finite piece   
C     born = |born|^2
C
C     pentagon and hexagon 




         do isig1=-1,1
            do isig3=-1,1
               do gluon_hel=1,2
                  
                  do div=0,2
                  NOABEmat21(isig1,isig3,gluon_hel,div) = czero
                  NOABEmat23(isig1,isig3,gluon_hel,div) = czero
                  NOABEmat41(isig1,isig3,gluon_hel,div) = czero
                  NOABEmat23(isig1,isig3,gluon_hel,div) = czero

                  enddo

                  do div =0,2
                     do j=1,3
                        do i=1,2
                           do iclass = 1,4
                              PentHex21(j,isig1,isig3,i,iclass,div)=czero
                              PentHex43(j,isig1,isig3,i,iclass,div)=czero
                              PentHex23(j,isig1,isig3,i,iclass,div)=czero
                              PentHex41(j,isig1,isig3,i,iclass,div)=czero
                           enddo
                        enddo
                     enddo
                     
                     do iclass=1,4
                        do i=1,2
                           HEX2143(2,isig1,isig3,i,iclass,div)=czero
                           HEX4123(2,isig1,isig3,i,iclass,div)=czero
                           
                           HEX2143(1,isig1,isig3,i,iclass,div)=czero
                           HEX4123(1,isig1,isig3,i,iclass,div)=czero
                        enddo
                     enddo
                  enddo
                  

                  do j=1,3
                  do div=0,2
                     
                     do i=1,5
                        mat43(i,isig1,isig3,gluon_hel,j,div) = czero
                        mat41(i,isig1,isig3,gluon_hel,j,div) = czero
                        mat23(i,isig1,isig3,gluon_hel,j,div) = czero
                        mat21(i,isig1,isig3,gluon_hel,j,div) = czero
                     enddo
                     do i=1,2
                        mat43g(i,isig1,isig3,gluon_hel,j,div) = czero
                        mat43gb(i,isig1,isig3,gluon_hel,j,div) = czero 

                        mat23g(i,isig1,isig3,gluon_hel,j,div) = czero
                        mat23gb(i,isig1,isig3,gluon_hel,j,div) = czero 

                        mat41g(i,isig1,isig3,gluon_hel,j,div) = czero
                        mat41gb(i,isig1,isig3,gluon_hel,j,div) = czero 
                        
                        mat21g(i,isig1,isig3,gluon_hel,j,div) = czero
                        mat21gb(i,isig1,isig3,gluon_hel,j,div) = czero 
                        
                     enddo
                     mat43b(isig1,isig3,gluon_hel,j,div) = czero
                     mat23b(isig1,isig3,gluon_hel,j,div) = czero
                     mat21b(isig1,isig3,gluon_hel,j,div) = czero 
                     mat41b(isig1,isig3,gluon_hel,j,div) = czero
                  enddo
               enddo
            enddo
         enddo
      enddo
C
C     Define diagrammatic momenta
C     HZ style labeling convention

c      print*,"inside ucHucg"
      
      do mu=0,3
         do I=1,6
            P(mu,I) = dble(psign(I))*PBAR(mu,I)
c            print*,"p(",mu,I,")=",p(mu,I)
         enddo

         q(mu) = p(mu,5)
      enddo
      q(4) = 0.0d0 ! massless gluon

c      do mu=0,3
c         print*,"sum mom = 0?",p(mu,1)+p(mu,3)-p(mu,2)-
c     $        p(mu,4)-p(mu,5)-p(mu,6)
c      enddo

C
      do I =1,4
         fsign(I) = psign(I)
      enddo
c
C     build 2-spinors
      call PSI0M( 4,PBAR,FSIGN,PSI )
c
c
c      alternate build of two-spinors
c
c      call bra(PBAR(0,2),PSI(1,-1,2))
c      call bra(PBAR(0,4),PSI(1,-1,4))

c      call ket(PBAR(0,1),PSI(1,-1,1))
c      call ket(PBAR(0,3),PSI(1,-1,3))
      
c
c      print*,"PSI=",PSI
C
C     construct currents J43^mu and J21^mu
C
      call curr(1,psi(1,-1,2),psi(1,-1,1),j21)
      call curr(1,psi(1,-1,4),psi(1,-1,3),j43)     
C     construct currents J41^mu and J23^mu
C
      call curr(1,psi(1,-1,2),psi(1,-1,3),j23)
      call curr(1,psi(1,-1,4),psi(1,-1,1),j41)


c gluon polarisation

c      do i=1,2
c         call polvec(pbar(0,5),i,epsR(0,i)) 
c      enddo

      do i=1,2
         call polvec(pbar(0,5),i,epsR(0,i))
         do mu=0,3
            if(lcheckgauge) then
               eps(mu,i) = dcmplx(pbar(mu,5),0d0) ! dcmplx(epsR(mu,i),0.0d0)
            else
               eps(mu,i) = dcmplx(epsR(mu,i),0.0d0)
            endif
            
         enddo

         if(ldebug) then

            do isig = -1,1,2
               call ket2r(psi(1,isig,1),.true.,p(0,1),isig,q,epsR(0,i),
     1              braket(1,isig,1,i),pk(0,1))
               call bra2r(psi(1,isig,2),.true.,p(0,2),isig,q,epsR(0,i),
     1              braket(1,isig,2,i),pk(0,2))
               call ket2r(psi(1,isig,3),.true.,p(0,3),isig,q,epsR(0,i),
     1              braket(1,isig,3,i),pk(0,3))
               call bra2r(psi(1,isig,4),.true.,p(0,4),isig,q,epsR(0,i),
     1              braket(1,isig,4,i),pk(0,4))
            enddo
         endif
      enddo
      if(ldebug) then
      do i = 1,2
         call curr(1,psi(1,-1,2),braket(1,-1,1,i),jh1)
         call curr(1,braket(1,-1,2,i),psi(1,-1,1),jh2)
         do isig=-1,1,2
            do mu = 0,3
               e21(mu,isig,i) = jh1(mu,isig) + jh2(mu,isig)
            enddo
         enddo
         
         call curr(1,psi(1,-1,4),braket(1,-1,3,i),jh1)
         call curr(1,braket(1,-1,4,i),psi(1,-1,3),jh2)
         do isig = -1,1,2
            do mu = 0,3
               e43(mu,isig,i) = jh1(mu,isig) + jh2(mu,isig)
            enddo
         enddo


         call curr(1,psi(1,-1,4),braket(1,-1,1,i),jh1)
         call curr(1,braket(1,-1,4,i),psi(1,-1,1),jh2)
         do isig = -1,1,2
            do mu = 0,3
               e41(mu,isig,i) = jh1(mu,isig) + jh2(mu,isig)
            enddo
         enddo

         call curr(1,psi(1,-1,2),braket(1,-1,3,i),jh1)
         call curr(1,braket(1,-1,2,i),psi(1,-1,3),jh2)
         do isig=-1,1,2
            do mu = 0,3
               e23(mu,isig,i) = jh1(mu,isig) + jh2(mu,isig)
            enddo
         enddo
         
       
         
      enddo
      endif


c     
c
c
      
c      print*,"inside ucHucg 1"
      
C     map momenta to Paco's notation
C     Perhap's perform this map in VBF_BOX?
      do mu=0,3
         p1(mu) = p(mu,1)
         p2(mu) = -1.0d0*p(mu,2)
         p3(mu) = p(mu,3)
         p4(mu) = -1.0d0*p(mu,4)
         pg(mu) = -1.0d0*p(mu,5) ! gluon
      enddo

c
c
c     map momentum to Paco's notation
c     Need to loop over isig1 isig2
c     gluon_hel is the gluon helicity
c     isig1 is the helicity of q1
c     isig2 is the helicity of q2
c     isig3 is the helicity of q3
c     div is 0,1,2 for eps^0 eps^-1, eps^-2 coeffiecients
c

c     Last index in mat() will be the graph id

      do mu=0,3
         p21g(mu) = -(p1(mu)+p2(mu)+pg(mu))
         p21(mu) = -(p1(mu)+p2(mu))

         p43g(mu) = -(p3(mu)+p4(mu)+pg(mu))
         p43(mu) = -(p3(mu)+p4(mu))      
      enddo
c   
      p21(4) = p21(0)**2-p21(1)**2-p21(2)**2-p21(3)**2
      p43(4) = p43(0)**2-p43(1)**2-p43(2)**2-p43(3)**2
      p21g(4) = p21g(0)**2-p21g(1)**2-p21g(2)**2-p21g(3)**2
      p43g(4) = p43g(0)**2-p43g(1)**2-p43g(2)**2-p43g(3)**2 
c
    
c
c     Build weak boson propagators
c
c     
      massV2 = bosonMass(1)**2
      widthV = bosonWidth(1) ! add to arguments
      
      prop21 = VPROP(p21(4),massV2,widthV)
      prop21g  = VPROP(p21g(4),massV2,widthV)
      prop43  = VPROP(p43(4),massV2,widthV)
      prop43g  = VPROP(p43g(4),massV2,widthV)
      
c
c     emission off <2|1> line
c Something here?
c
cFC  
      
      call BOX_VBF(psi(1,-1,2),psi(1,-1,1),p1,p2,p21g,pg,eps,j43,IsVaW(1),divmax,scale,
     $     mat21g,mat21,mat21gb,mat21b,NOABEmat21,NOABEmat21g,wardbox1)

      if(wardbox1) then
       gaugebad(1,1)=gaugebad(1,1)+1
       else
       gaugebad(1,2)=gaugebad(1,2)+1
      endif

      if(ldebug) then

       do isig1=-1,1,2
         do isig3=-1,1,2
            do gluon_hel=1,2
               

               print*,"NO ABE ratio=",((mat21(2,isig1,isig3,gluon_hel,1,2)+
     $              mat21(2,isig1,isig3,gluon_hel,2,2))*CA*(-1d0/2d0)+
     $              NOABEmat21(isig1,isig3,gluon_hel,2)*CA)/(-1d0*CA*
     $              (mat21b(isig1,isig3,gluon_hel,1,0)+
     $              mat21b(isig1,isig3,gluon_hel,2,0)))

               test = prop43*prop21g*(mat21b(isig1,isig3,gluon_hel,1,0)+
     $              mat21b(isig1,isig3,gluon_hel,2,0))

               print*,"box_mat21b(",isig1,",",isig3,",",gluon_hel,")=",
     $              test

    
c               test = mat23b(isig1,isig3,gluon_hel,1,0)+
c     $              mat23b(isig1,isig3,gluon_hel,2,0)

             
                test1 = (1.0d0)*prop43*prop21g*
     $              dotcc(e21(0,isig1,gluon_hel),j43(0,isig3))


c               test1 = prop41g*prop23*
c     $              dotcc(e41(0,isig1,gluon_hel),j23(0,isig3)) + 
c     $               prop23g*prop41*
c     $              dotcc(e23(0,isig3,gluon_hel),j41(0,isig1))

               print*,"ratio=",abs(test)/abs(test1)
               
               print*,"box_mat21b_p(",isig1,",",isig3,",",gluon_hel,")=",
     $              test1
            print*,"====="
            enddo
         enddo
      enddo
      endif

c
c     Here we compute a current trifac43(0:2) 
c
     
      call VBF_TRI(tri43,p43(4),scale)


CFC VERT
      if(VER) then
c FC 
       do mu=0,3
          ptemp12(mu)=-p1(mu)-p2(mu)
       enddo   
       ptemp12sq=ptemp12(0)**2-ptemp12(1)**2-ptemp12(2)**2-ptemp12(3)**2
 
      
      lns=log(Abs(ptemp12sq))-(0,1)*pi*theta(ptemp12sq)-log(scale)
c      print*, "p2sq 2",ptemp12sq

      ctetemp(0)=-lns*lns+3d0*lns -8d0 + pi*pi/3d0
      ctetemp(1)=2*lns-3d0
      ctetemp(2)=-2d0

      do mu=0,2
      call Verline(p1,ptemp12,p2,psi(1,-1,2),psi(1,-1,1),eps(0,1),-1,scale,
     &   1,trifac43(1),trifac43(2),mu)

      print*,"fact VERTEX",trifac43(1)
c      print*,"fact born",trifac43(2)
      print*, "ratio 1",trifac43(1)/(trifac43(2)*ctetemp(mu))
      print*, "ratio 2",tri21(mu)/(ctetemp(mu))
      print*

      enddo

      stop
      endif


C     emission off the <4|3> line

     

      call BOX_VBF(psi(1,-1,4),psi(1,-1,3),p3,p4,p43g,pg,eps,j21,IsVaW(1),divmax,scale,
     $     mat43g,mat43,mat43gb,mat43b,NOABEmat43,NOABEmat43g,wardbox2)

      if(wardbox2) then
       gaugebad(1,1)=gaugebad(1,1)+1
       else
       gaugebad(1,2)=gaugebad(1,2)+1
      endif



      if(ldebug) then 

       do isig1=-1,1,2
         do isig3=-1,1,2
            do gluon_hel=1,2

               print*,"NO ABE ratio=",((mat43(2,isig1,isig3,gluon_hel,1,2)+
     $              mat43(2,isig1,isig3,gluon_hel,2,2))*CA*(-1d0/2d0)+
     $              NOABEmat43(isig1,isig3,gluon_hel,2)*CA)/(-1d0*CA*
     $              (mat43b(isig1,isig3,gluon_hel,1,0)+
     $              mat43b(isig1,isig3,gluon_hel,2,0)))



               test = prop43g*prop21*(mat43b(isig3,isig1,gluon_hel,1,0)+
     $              mat43b(isig3,isig1,gluon_hel,2,0))

               print*,"box_mat43b(",isig1,",",isig3,",",gluon_hel,")=",
     $              test

    
c               test = mat23b(isig1,isig3,gluon_hel,1,0)+
c     $              mat23b(isig1,isig3,gluon_hel,2,0)

             
                test1 = (1.0d0)*prop43g*prop21*
     $              dotcc(e43(0,isig3,gluon_hel),j21(0,isig1))


c               test1 = prop41g*prop23*
c     $              dotcc(e41(0,isig1,gluon_hel),j23(0,isig3)) + 
c     $               prop23g*prop41*
c     $              dotcc(e23(0,isig3,gluon_hel),j41(0,isig1))

               print*,"ratio=",test/test1
               
               print*,"box_mat43b_p(",isig1,",",isig3,",",gluon_hel,")=",
     $              test1
            print*,"====="
            enddo
         enddo
      enddo

      endif

c      
c     Here we compute a current trifac21(0:2) 
c     
      
      call VBF_TRI(tri21,p21(4),scale)
     


c      print*,"inside ucHucg 2"
c
cccccccccccccccccccccccccccccc
c
C     Only needed for indentical particles in the 2 and 4 

      do mu=0,3
         p23g(mu) = -(p3(mu)+p2(mu)+pg(mu))
         p23(mu) = -(p3(mu)+p2(mu))

         p41g(mu) = -(p1(mu)+p4(mu)+pg(mu))
         p41(mu) = -(p1(mu)+p4(mu))
      enddo


      p23(4) = p23(0)**2-p23(1)**2-p23(2)**2-p23(3)**2
      p41(4) = p41(0)**2-p41(1)**2-p41(2)**2-p41(3)**2
      p23g(4) = p23g(0)**2-p23g(1)**2-p23g(2)**2-p23g(3)**2
      p41g(4) = p41g(0)**2-p41g(1)**2-p41g(2)**2-p41g(3)**2  

      massV2 = bosonMass(2)**2
      widthV = bosonWidth(2) ! add to arguments
      
      prop41 = VPROP(p41(4),massV2,widthV)
      prop41g = VPROP(p41g(4),massV2,widthV)
      prop23 = VPROP(p23(4),massV2,widthV)
      prop23g = VPROP(p23g(4),massV2,widthV)
       
C     Only needed for indentical particles in the 2 and 4 

C     emission off the <2|3> line 

      call BOX_VBF(psi(1,-1,2),psi(1,-1,3),p3,p2,p23g,pg,eps,j41,IsVaW(2),divmax,scale,
     $     mat23g,mat23,mat23gb,mat23b,NOABEmat23,NOABEmat23g,wardbox3)

      if(wardbox3) then
       gaugebad(1,1)=gaugebad(1,1)+1
       else
       gaugebad(1,2)=gaugebad(1,2)+1
      endif



      if(ldebug) then
      do isig1=-1,1,2
         do isig3=-1,1,2
            do gluon_hel=1,2

               print*,"NO ABE ratio=",((mat23(2,isig1,isig3,gluon_hel,1,2)+
     $              mat23(2,isig1,isig3,gluon_hel,2,2))*CA*(-1d0/2d0)+
     $              NOABEmat23(isig1,isig3,gluon_hel,2)*CA)/(-1d0*CA*
     $              (mat23b(isig1,isig3,gluon_hel,1,0)+
     $              mat23b(isig1,isig3,gluon_hel,2,0)))


               
               test = prop23g*prop41*(mat23b(isig3,isig1,gluon_hel,1,0)+
     $              mat23b(isig3,isig1,gluon_hel,2,0))

               print*,"box_mat23b(",isig1,",",isig3,",",gluon_hel,")=",
     $              test

    
c               test = mat23b(isig1,isig3,gluon_hel,1,0)+
c     $              mat23b(isig1,isig3,gluon_hel,2,0)

             
                test1 = prop23g*prop41*
     $              dotcc(e23(0,isig3,gluon_hel),j41(0,isig1))


c               test1 = prop41g*prop23*
c     $              dotcc(e41(0,isig1,gluon_hel),j23(0,isig3)) + 
c     $               prop23g*prop41*
c     $              dotcc(e23(0,isig3,gluon_hel),j41(0,isig1))

               print*,"ratio=",test/test1
               
               print*,"box_mat23b_p(",isig1,",",isig3,",",gluon_hel,")=",
     $              test1
            print*,"====="
            enddo
         enddo
      enddo
      endif

c     Here we compute a current trifac41(0:2) 
c     
      
      call VBF_TRI(tri41,p41(4),scale)

c      print*,mat23,mat23b

c     emission off the <4|1> line

      call BOX_VBF(psi(1,-1,4),psi(1,-1,1),p1,p4,p41g,pg,eps,j23,IsVaW(2),divmax,scale,
     $     mat41g,mat41,mat41gb,mat41b,NOABEmat41,NOABEmat41g,wardbox4)

      if(wardbox4) then
       gaugebad(1,1)=gaugebad(1,1)+1
       else
       gaugebad(1,2)=gaugebad(1,2)+1
      endif



      if(ldebug) then
       do isig1=-1,1,2
         do isig3=-1,1,2
            do gluon_hel=1,2


               print*,"NO ABE ratio=",((mat41(2,isig1,isig3,gluon_hel,1,2)+
     $              mat41(2,isig1,isig3,gluon_hel,2,2))*CA*(-1d0/2d0)+
     $              NOABEmat41(isig1,isig3,gluon_hel,2)*CA)/(-1d0*CA*
     $              (mat41b(isig1,isig3,gluon_hel,1,0)+
     $              mat41b(isig1,isig3,gluon_hel,2,0)))


               test = prop41g*prop23*(mat41b(isig1,isig3,gluon_hel,1,0)+
     $              mat41b(isig1,isig3,gluon_hel,2,0))

               print*,"box_mat41b(",isig1,",",isig3,",",gluon_hel,")=",
     $              test

    
c               test = mat23b(isig1,isig3,gluon_hel,1,0)+
c     $              mat23b(isig1,isig3,gluon_hel,2,0)

             
                test1 = prop41g*prop23*
     $              dotcc(e41(0,isig1,gluon_hel),j23(0,isig3))


c               test1 = prop41g*prop23*
c     $              dotcc(e41(0,isig1,gluon_hel),j23(0,isig3)) + 
c     $               prop23g*prop41*
c     $              dotcc(e23(0,isig3,gluon_hel),j41(0,isig1))

               print*,"ratio=", test/test1
               
               print*,"box_mat41b_p(",isig1,",",isig3,",",gluon_hel,")=",
     $              test1
            print*,"====="
            enddo
         enddo
      enddo
      endif

           
c     Here we compute a current trifac23(0:2) 
c     
      
      call VBF_TRI(tri23,p23(4),scale)


      
c      print*,"inside ucHucg 3"


C     Assembly color sub-amplitudes
C     This convention can be adjusted later.

C     color index 1: T(2,b,1) * delta(4,3)
C                 2: T(4,b,3) * delta(2,1)
C                 3: T(2,b,3) * delta(4,1)
C                 4: T(4,b,1) * delta(2,3)
C                 5: T(4,a,3) * T(2,b,a,1) 
C                 6: T(4,a,3) * T(2,a,b,1)
C                 7: T(2,a,1) * T(4,b,a,3)
C                 8: T(2,a,1) * T(4,a,b,3)
C                 9: T(2,a,3) * T(4,b,a,1)
C                10: T(2,a,3) * T(4,a,b,1)
C                11: T(4,a,1) * T(2,b,a,3)
c                12: T(4,a,1) * T(2,a,b,3)
C
C
C    

c     initialize jamp_virt, jamp_born

      do isig1=-1,1
         do isig2=-1,1
            do isig3=-1,1
               do isig4=-1,1
                  do gluon_hel = 1,2
                     do icol=1,4
                        jamp_born(isig1,isig2,isig3,isig4,
     $                       gluon_hel,icol) = 
     $                       czero
                        do div=0,2
                        amp_virt(isig1,isig2,isig3,isig4,
     $                       gluon_hel,icol,div) = 
     $                      czero
                        enddo
                     enddo
                     do div=0,2
                        do icol=1,14
                           jamp_virt(isig1,isig2,isig3,isig4,
     $                          gluon_hel,icol,div) = 
     $                          czero
                        enddo  
                     enddo
                  enddo
               enddo
            enddo
         enddo
      enddo
C     Only compute is coupling is nonzero
c
c


      do isig1=-1,1,2
         do isig3 = -1,1,2
            do gluon_hel = 1,2
               do div=0,2
                  if(compBoxes) then
                  if(abs(CLR(flavor(3),1,isig3)*
     $                 CLR(flavor(1),1,isig1)).gt.0d0) 
     $                 then

                     isig2 = isig1
                     isig4 = isig3

                     jamp_virt(isig1,isig2,isig3,isig4,gluon_hel,1,div) = 
     $                    prop21g*prop43*GHVV(1)*
     $                    CLR(flavor(3),1,isig3)*
     $                    CLR(flavor(1),1,isig1)*                 
     $                    (CF*(mat21(1,isig1,isig3,gluon_hel,1,div) + 
     $                    mat21(1,isig1,isig3,gluon_hel,2,div)) +
     $                    (CF-1.0d0/2.0d0*CA)*
     $                    (mat21(2,isig1,isig3,gluon_hel,1,div) + 
     $                    mat21(2,isig1,isig3,gluon_hel,2,div)) +
     $                    CA*NOABEmat21(isig1,isig3,gluon_hel,div) +
     $                    CF*tri43(div)*
     $                    (mat21b(isig1,isig3,gluon_hel,1,0) + 
     $                    mat21b(isig1,isig3,gluon_hel,2,0)))
c               
                  jamp_virt(isig1,isig2,isig3,isig4,gluon_hel,2,div) = 
     $                 prop43g*prop21*GHVV(1)*
     $                 CLR(flavor(3),1,isig3)*
     $                 CLR(flavor(1),1,isig1)*
     $                 (CF*(mat43(1,isig3,isig1,gluon_hel,1,div) + 
     $                 mat43(1,isig3,isig1,gluon_hel,2,div)) +
     $                 (CF-1.0d0/2.0d0*CA)*
     $                 (mat43(2,isig3,isig1,gluon_hel,1,div) + 
     $                 mat43(2,isig3,isig1,gluon_hel,2,div)) +
     $                 CA*NOABEmat43(isig3,isig1,gluon_hel,div) +
     $                 CF*tri21(div)*
     $                    (mat43b(isig3,isig1,gluon_hel,1,0) + 
     $                      mat43b(isig3,isig1,gluon_hel,2,0)))


               endif

C     Only need to W - Z interference and identical quarks in final state
               if(abs(CLR(flavor(3),2,isig3)*
     $              CLR(flavor(1),2,isig1)).gt.0d0) 
     $                 then
                  isig2 = isig3
                  isig4 = isig1

                   jamp_virt(isig1,isig2,isig3,isig4,gluon_hel,3,div) = 
     $                 SgnF*prop23g*prop41*GHVV(2)*
     $                 CLR(flavor(3),2,isig3)*
     $                 CLR(flavor(1),2,isig1)*                 
     $                 (CF*(mat23(1,isig3,isig1,gluon_hel,1,div) + 
     $                 mat23(1,isig3,isig1,gluon_hel,2,div)) +
     $                 (CF-1.0d0/2.0d0*CA)*
     $                 (mat23(2,isig3,isig1,gluon_hel,1,div) + 
     $                 mat23(2,isig3,isig1,gluon_hel,2,div)) +
     $                 CA*NOABEmat23(isig3,isig1,gluon_hel,div) +
     $                 CF*tri41(div)*
     $                    (mat23b(isig3,isig1,gluon_hel,1,0) + 
     $                      mat23b(isig3,isig1,gluon_hel,2,0)))

                  jamp_virt(isig1,isig2,isig3,isig4,gluon_hel,4,div) = 
     $                 SgnF*prop41g*prop23*GHVV(2)*
     $                 CLR(flavor(3),2,isig3)*
     $                 CLR(flavor(1),2,isig1)*
     $                 (CF*(mat41(1,isig1,isig3,gluon_hel,1,div) + 
     $                 mat41(1,isig1,isig3,gluon_hel,2,div)) +
     $                 (CF-1.0d0/2.0d0*CA)*
     $                 (mat41(2,isig1,isig3,gluon_hel,1,div) + 
     $                 mat41(2,isig1,isig3,gluon_hel,2,div)) +
     $                 CA*NOABEmat41(isig1,isig3,gluon_hel,div) +
     $                 CF*tri23(div)*
     $                    (mat41b(isig1,isig3,gluon_hel,1,0) + 
     $                      mat41b(isig1,isig3,gluon_hel,2,0)))



               endif
               endif
            enddo               ! loop over div
C     born 
               
               if(abs(CLR(flavor(3),1,isig3)*
     $              CLR(flavor(1),1,isig1)).gt.0d0) 
     $              then
                  isig2 = isig1
                  isig4 = isig3

                  jamp_born(isig1,isig2,isig3,isig4,gluon_hel,1) =  
     $                 prop21g*prop43*GHVV(1)*
     $                 CLR(flavor(3),1,isig3)*
     $                 CLR(flavor(1),1,isig1)*                 
     $                 (mat21b(isig1,isig3,gluon_hel,1,0) + 
     $                 mat21b(isig1,isig3,gluon_hel,2,0))

c     alternate version
c                  jamp_born(isig1,isig3,gluon_hel,1) = prop21g*prop43*
c     $                 CLR(flavor(3),1,isig3)*CLR(flavor(1),1,isig1)*
c     $                 dotcc(e21(0,isig1,gluon_hel),j43(0,isig3))
                  
                  
                  jamp_born(isig1,isig2,isig3,isig4,gluon_hel,2) = 
     $                 prop43g*prop21*GHVV(1)*
     $                 CLR(flavor(3),1,isig3)*
     $                 CLR(flavor(1),1,isig1)*                 
     $                 (mat43b(isig3,isig1,gluon_hel,1,0) + 
     $                 mat43b(isig3,isig1,gluon_hel,2,0))

c                  jamp_born(isig1,isig3,gluon_hel,2) = prop21*prop43g*
c     $                 CLR(flavor(3),1,isig3)*CLR(flavor(1),1,isig1)*
c     $                 dotcc(e43(0,isig3,gluon_hel),j21(0,isig1))
                  

               endif
C     only needed for same generations/quarks

               if(abs(CLR(flavor(3),2,isig3)*
     $              CLR(flavor(1),2,isig1)).gt.0d0) 
     $              then
                  isig2 = isig3
                  isig4 = isig1

                  jamp_born(isig1,isig2,isig3,isig4,gluon_hel,3) =  
     $              SgnF*prop23g*prop41*GHVV(2)*
     $                 CLR(flavor(3),2,isig3)*
     $                 CLR(flavor(1),2,isig1)*                 
     $                 (mat23b(isig3,isig1,gluon_hel,1,0) + 
     $                 mat23b(isig3,isig1,gluon_hel,2,0))

c                   jamp_born(isig1,isig3,gluon_hel,3) = 
c     $                 SgnF*prop23g*prop41*
c     $                 CLR(flavor(3),2,isig3)*CLR(flavor(1),2,isig1)*
c     $                 dotcc(e23(0,isig3,gluon_hel),j41(0,isig1))
                  
                  
                  jamp_born(isig1,isig2,isig3,isig4,gluon_hel,4) = 
     $                 SgnF*prop41g*prop23*GHVV(2)*
     $                 CLR(flavor(3),2,isig3)*
     $                 CLR(flavor(1),2,isig1)*                 
     $                 (mat41b(isig1,isig3,gluon_hel,1,0) + 
     $              mat41b(isig1,isig3,gluon_hel,2,0))


c                  jamp_born(isig1,isig3,gluon_hel,4) = 
c     $                 SgnF*prop23*prop41g*
c     $                 CLR(flavor(3),2,isig3)*CLR(flavor(1),2,isig1)*
c     $                 dotcc(e41(0,isig1,gluon_hel),j23(0,isig3))
                  
               endif
               

             
            enddo
         enddo
      enddo
C     
     
c      print*,"inside ucHucg 4"


      if(compPentagons) then
CCCC  Pentagon and hexagons are computed here.

C     
C     call Paco's routines to compute hexagons and pentagon diagrams

     
cFC                !M! complex! 
cFC   FeyCalc notation P-=7 P+=6!!!

CTF   The first 7 indicates the helicity of k1 and the second 7 indicates the helicity of k4.
C
C     PentHex(48,isig1,isig3,gluon_hel,iclass,div) (only needed 1-3) store the amplitude piece for 
C     particular set of helicities
C     PentHexn(isig1,isig3,gluon_hel,iclass,div) is the result of 0 gluons.
C     Need to loop over gluon polarization and div
C
C     iclass is 1 for direct, 2 for cross, 3 for crossF, and 4 for crossIF 

c     emission off 21 

c     Fill complex mass parameter for hexagon and pentagons
c
      M2_c(1) = dcmplx(bosonMass(1)**2,-bosonMass(1)*bosonWidth(1))
      M2_c(2) = dcmplx(bosonMass(2)**2,-bosonMass(2)*bosonWidth(2))

c      print*, "M2_c(2)",M2_c(2)

c     Shall add born piece for testing purposes.
c     born21(isig1,isig3,gluon_hel)

      call PENTHEX_VBF(PentHex21,born21,M2_c(1),
     $     p(0,1),PSI(1,-1,1),p(0,2),PSI(1,-1,2),
     $     p(0,3),PSI(1,-1,3),p(0,4),PSI(1,-1,4),
     $     p(0,5),eps,p(0,6),scale,IsVaW(1),divmax,wardhex1)

      if(wardhex1) then
       gaugebad(2,1)=gaugebad(2,1)+1
       else
       gaugebad(2,2)=gaugebad(2,2)+1
      endif

      if(ldebug) then
      print*,"testing pentagons ---"

      do iclass=1,4
         print*,"iclass=",iclass
      do isig1=-1,1,2
         do isig3=-1,1,2
            do gluon_hel=1,2
               
               print*,"born21(",isig1,",",isig3,",",gluon_hel,iclass,")=",
     $              born21(isig1,isig3,gluon_hel,iclass)

               test = born21(isig1,isig3,gluon_hel,iclass)

               test1 =
     $              prop21g*prop43*
     $              dotcc(e21(0,isig1,gluon_hel),j43(0,isig3))

c               test1 = born21(isig1,isig3,gluon_hel,1)

               print*,"ratio=",test/test1
               
               print*,"born21_p(",isig1,",",isig3,",",gluon_hel,iclass,")=",
     $              test1
            print*,"====="
            enddo
         enddo
      enddo
      enddo
      endif

c      stop
C    emission off 43 

c      print*,"testing pentagons emission off 43 line"
C
c      ! YOU MUST SWAP isig1 and isig3 !

      call PENTHEX_VBF(PentHex43,born43,M2_c(1),
     $     p(0,3),PSI(1,-1,3),p(0,4),PSI(1,-1,4),
     $     p(0,1),PSI(1,-1,1),p(0,2),PSI(1,-1,2),
     $     p(0,5),eps,p(0,6),scale,IsVaW(1),divmax,wardhex2)

      if(wardhex2) then
       gaugebad(2,1)=gaugebad(2,1)+1
       else
       gaugebad(2,2)=gaugebad(2,2)+1
      endif

      if(ldebug) then
      do iclass=1,4
         print*,"iclass=",iclass

       do isig1=-1,1,2
         do isig3=-1,1,2
            do gluon_hel=1,2
               
            print*,"born43(",isig3,",",isig1,",",gluon_hel,iclass,")=",
     $              born43(isig3,isig1,gluon_hel,iclass)

            test = born43(isig3,isig1,gluon_hel,iclass)


            test1=prop21*prop43g*
     $           dotcc(e43(0,isig3,gluon_hel),j21(0,isig1))

c            test1 = born43(isig1,isig3,gluon_hel,1)

            print*,"ratio=",test/test1
      
            print*,"born43_p(",isig3,",",isig1,",",gluon_hel,iclass,")=",
     $          test1
            print*,"====="
         enddo
         enddo
         enddo
      enddo
      endif
c
c      stop
C     emission off 41 

c      print*,"emission off the 41 line"

      call PENTHEX_VBF(PentHex41,born41,M2_c(2),
     $     p(0,1),PSI(1,-1,1),p(0,4),PSI(1,-1,4),
     $     p(0,3),PSI(1,-1,3),p(0,2),PSI(1,-1,2),
     $     p(0,5),eps,p(0,6),scale,IsVaW(2),divmax,wardhex3)

      if(wardhex3) then
       gaugebad(2,1)=gaugebad(2,1)+1
       else
       gaugebad(2,2)=gaugebad(2,2)+1
      endif
C
c$$$       do isig1=-1,1,2
c$$$         do isig3=-1,1,2
c$$$            
c$$$            print*,"born41=",born41(isig1,isig3,1,1)
c$$$   
c$$$         enddo
c$$$      enddo

      if(ldebug) then
      do iclass = 1,4
         print*,"iclass =",iclass
         do isig1=-1,1,2
            do isig3=-1,1,2
               do gluon_hel=1,2
                  print*,"born41(",isig1,",",isig3,",",gluon_hel,")=",
     $                 born41(isig1,isig3,gluon_hel,iclass)

                  test = born41(isig1,isig3,gluon_hel,iclass)
                  
                  test1 =
     $                 prop41g*prop23*
     $                 dotcc(e41(0,isig1,gluon_hel),j23(0,isig3))
                  
c     test1 = born41(isig1,isig3,gluon_hel,1)
                  
                  print*,"ratio=",test/test1
               
                  print*,"born41_p(",isig1,",",isig3,",",gluon_hel,")=",
     $                 test1
                  print*,"====="
               enddo
            enddo
         enddo
      enddo
c      stop

C     only needed for indentical flavors/same generation
C     Need master formula
C
C     emission off 23 line

      print*,"emission off the 23 line"

      endif

      call PENTHEX_VBF(PentHex23,born23,M2_c(2),
     $     p(0,3),PSI(1,-1,3),p(0,2),PSI(1,-1,2),
     $     p(0,1),PSI(1,-1,1),p(0,4),PSI(1,-1,4),
     $     p(0,5),eps,p(0,6),scale,IsVaW(2),divmax,wardhex4)


      if(wardhex4) then
       gaugebad(2,1)=gaugebad(2,1)+1
       else
       gaugebad(2,2)=gaugebad(2,2)+1
      endif

c     you must swap isig1 and isig3 !
       

      if(ldebug) then
      do iclass = 1,4
         print*,"iclass=",iclass
         do isig1=-1,1,2
            do isig3=-1,1,2
               do gluon_hel=1,2
                  print*,"born23(",isig3,",",isig1,",",gluon_hel,")=",
     $                 born23(isig3,isig1,gluon_hel,iclass)
                  
                  test = born23(isig3,isig1,gluon_hel,iclass)
                  
                  test1 =
     $                 prop23g*prop41*
     $                 dotcc(e23(0,isig3,gluon_hel),j41(0,isig1))
                  
c     test1 = born23(isig1,isig3,gluon_hel,1)
                  
                  print*,"ratio=",test/test1
                  
                  print*,"born23_p(",isig3,",",isig1,",",gluon_hel,")=",
     $                 test1
                  print*,"====="
               enddo
            enddo
         enddo
      enddo
      endif
c      stop

c      print*,"inside ucHucg 5"
c
c
c    1 ---->----- 2       1 ---->----4
c
c    3 ---->------4       3----->----2 
c

C
c     Need to compute non-abelian hexagon 
C     Only needed for same gen quarks or same flavor quarks
C
c     Check this color factor
c
c
      call NOABE_HEX_VBF(Hex2143,Born2143,M2_c(1),
     $     p(0,1),PSI(1,-1,1),p(0,2),PSI(1,-1,2),
     $     p(0,3),PSI(1,-1,3),p(0,4),PSI(1,-1,4),
     $     p(0,5),eps,p(0,6),scale,IsVaW(1),divmax,wardhexN1)

      if(wardhexN1) then
       gaugebad(3,1)=gaugebad(3,1)+1
       else
       gaugebad(3,2)=gaugebad(3,2)+1
      endif


      if(ldebug) then
c      print*,"HEX2143g=",HEX2143(2,

      print*,"testing non-abelian hexagons"
c      iclass = 2    
c      print*,"testing iclass=",iclass
      do iclass = 1,4
         print*,"iclass=",iclass
       do isig1=-1,1,2
         do isig3=-1,1,2
            do gluon_hel=1,2
c               print*,"born2143(",isig1,",",isig3,",",gluon_hel,")=",
c     $              born2143(isig1,isig3,gluon_hel,iclass)

               do div = 0,2
                print*,"HEX ---- iclass = ",iclass
                print*,"HEX ---- sig1=",isig1,"isig3=",isig3,
     $                           "gluon_hel=",gluon_hel,"div=",div
                print*,"HEX2143g=",HEX2143(2,isig1,isig3,gluon_hel,iclass,div)
                print*,"HEX2143 =",HEX2143(1,isig1,isig3,gluon_hel,iclass,div)
                print*,"HEX ratio =", HEX2143(2,isig1,isig3,gluon_hel,iclass,div)
     $          /HEX2143(1,isig1,isig3,gluon_hel,iclass,div)
                print*,"HEX ----"
             enddo

                print*,"born2143_1(",isig1,",",isig3,",",gluon_hel,")=",
     $              born2143(isig1,isig3,gluon_hel,1)

               print*,"born2143_2(",isig1,",",isig3,",",gluon_hel,")=",
     $              born2143(isig1,isig3,gluon_hel,2)

               print*,"born2143_3(",isig1,",",isig3,",",gluon_hel,")=",
     $              born2143(isig1,isig3,gluon_hel,3)

                print*,"born2143_4(",isig1,",",isig3,",",gluon_hel,")=",
     $              born2143(isig1,isig3,gluon_hel,4)


c               test = born2143(isig1,isig3,gluon_hel,iclass)*(-1.0d0)

               if(iclass.eq.1 .or. iclass.eq.3) then
                   test = born2143(isig1,isig3,gluon_hel,iclass)*(-1.0d0)
                elseif(iclass.eq.2 .or. iclass.eq.4) then
                   test = born2143(isig3,isig1,gluon_hel,iclass)*(-1.0d0)
                endif

             
               test1 = prop21g*prop43*
     $              dotcc(e21(0,isig1,gluon_hel),j43(0,isig3)) + 
     $               prop43g*prop21*
     $              dotcc(e43(0,isig3,gluon_hel),j21(0,isig1))

               print*,"ratio=",test/test1
               
               print*,"born2143_p(",isig1,",",isig3,",",gluon_hel,")=",
     $              test1
            print*,"====="
            enddo
         enddo
      enddo
      enddo
      endif
c      stop

c     Only needed for same gens/identical quarks

c
       call NOABE_HEX_VBF(Hex4123,Born4123,M2_c(2),
     $     p(0,1),PSI(1,-1,1),p(0,4),PSI(1,-1,4),
     $     p(0,3),PSI(1,-1,3),p(0,2),PSI(1,-1,2),
     $     p(0,5),eps,p(0,6),scale,IsVaW(2),divmax,wardhexN2) 


      if(wardhexN2) then
       gaugebad(3,1)=gaugebad(3,1)+1
       else
       gaugebad(3,2)=gaugebad(3,2)+1
      endif


       if(ldebug) then
       do iclass = 1,4
c       iclass = 1
          print*,"iclass=",iclass
        do isig1=-1,1,2
         do isig3=-1,1,2
            do gluon_hel=1,2

               do div = 0,2
                print*,"HEX ---- iclass = ",iclass
                print*,"HEX ---- sig1=",isig1,"isig3=",isig3,
     $       "gluon_hel=",
     $  gluon_hel,"div=",div
                print*,"HEX4123g=",HEX4123(2,isig1,isig3,gluon_hel,iclass,div)
                print*,"HEX4123 =",HEX4123(1,isig1,isig3,gluon_hel,iclass,div)
                print*,"HEX ratio =", HEX4123(2,isig1,isig3,gluon_hel,iclass,div)
     $       /HEX4123(1,isig1,isig3,gluon_hel,iclass,div)
                print*,"HEX ----"
             enddo
               print*,"born4123_1(",isig1,",",isig3,",",gluon_hel,")=",
     $              born4123(isig1,isig3,gluon_hel,1)

               print*,"born4123_2(",isig1,",",isig3,",",gluon_hel,")=",
     $              born4123(isig1,isig3,gluon_hel,2)

               print*,"born4123_3(",isig1,",",isig3,",",gluon_hel,")=",
     $              born4123(isig1,isig3,gluon_hel,3)

                print*,"born4123_4(",isig1,",",isig3,",",gluon_hel,")=",
     $              born4123(isig1,isig3,gluon_hel,4)

c               test = born4123(isig1,isig3,gluon_hel,iclass)*
c     $               SgnF*(-1.0d0)

                if(iclass.eq.1 .or. iclass.eq.3) then
                   test = born4123(isig1,isig3,gluon_hel,iclass)*(-1.0d0)
     $                  *SgnF
                elseif(iclass.eq.2 .or. iclass.eq.4) then
                   test = born4123(isig3,isig1,gluon_hel,iclass)*(-1.0d0)
     $                  *SgnF
                endif
c swap isig1 and isig3 ?
             
               test1 = SgnF*(prop41g*prop23*
     $              dotcc(e41(0,isig1,gluon_hel),j23(0,isig3)) + 
     $               prop23g*prop41*
     $              dotcc(e23(0,isig3,gluon_hel),j41(0,isig1)))

               print*,"ratio=",test/test1
               
               print*,"born4123_p(",isig1,",",isig3,",",gluon_hel,")=",
     $              test1
            print*,"====="
            enddo
         enddo
      enddo
      enddo
      endif

c
c     
c      fac2 = 1d0

       do div = 0,2
         do isig1=-1,1,2
            do isig3 = -1,1,2
               do gluon_hel = 1,2
                  
                   if(abs(CLR(flavor(3),1,isig3)*
     $                 CLR(flavor(1),1,isig1)).gt.0d0) 
     $                 then
                      isig2 = isig1
                      isig4 = isig3
                  jamp_virt(isig1,isig2,isig3,isig4,gluon_hel,5,div) = 
     $                 (-1.0d0)*CLR(flavor(3),1,isig3)*GHVV(1)*
     $                 CLR(flavor(1),1,isig1)*
     $                 (PentHex21(3,isig1,isig3,gluon_hel,1,div) + 
     $                  PentHex21(2,isig1,isig3,gluon_hel,1,div) +
     $                  PentHex21(3,isig1,isig3,gluon_hel,2,div) +
     $                  PentHex21(2,isig1,isig3,gluon_hel,2,div) +
     $                  PentHex21(2,isig1,isig3,gluon_hel,3,div) +
     $                  PentHex21(2,isig1,isig3,gluon_hel,4,div))
                     
                  jamp_virt(isig1,isig2,isig3,isig4,gluon_hel,6,div) = 
     $                 (-1.0d0)*CLR(flavor(3),1,isig3)*GHVV(1)*
     $                 CLR(flavor(1),1,isig1)* 
     $                 (PentHex21(1,isig1,isig3,gluon_hel,1,div) + 
     $                  PentHex21(1,isig1,isig3,gluon_hel,2,div) +
     $                  PentHex21(3,isig1,isig3,gluon_hel,3,div) +
     $                  PentHex21(1,isig1,isig3,gluon_hel,3,div) +
     $                  PentHex21(3,isig1,isig3,gluon_hel,4,div) +
     $                  PentHex21(1,isig1,isig3,gluon_hel,4,div))
     $                
c

                    
c TF edit
                  jamp_virt(isig1,isig2,isig3,isig4,gluon_hel,7,div) = 
     $                 CLR(flavor(3),1,isig3)*GHVV(1)*
     $                 CLR(flavor(1),1,isig1)*(-1.0d0)*
     $                 (PentHex43(3,isig3,isig1,gluon_hel,1,div) + 
     $                  PentHex43(2,isig3,isig1,gluon_hel,1,div) +
     $                  PentHex43(3,isig3,isig1,gluon_hel,2,div) +
     $                  PentHex43(2,isig3,isig1,gluon_hel,2,div) +
     $                  PentHex43(2,isig3,isig1,gluon_hel,3,div) +
     $                  PentHex43(2,isig3,isig1,gluon_hel,4,div)) 
                      
                  jamp_virt(isig1,isig2,isig3,isig4,gluon_hel,8,div) = 
     $                 CLR(flavor(3),1,isig3)*GHVV(1)*
     $                 CLR(flavor(1),1,isig1)*(-1.0d0)*
     $                 (PentHex43(1,isig3,isig1,gluon_hel,1,div) + 
     $                  PentHex43(1,isig3,isig1,gluon_hel,2,div) +
     $                  PentHex43(3,isig3,isig1,gluon_hel,3,div) +
     $                  PentHex43(1,isig3,isig1,gluon_hel,3,div) +
     $                  PentHex43(3,isig3,isig1,gluon_hel,4,div) +
     $                  PentHex43(1,isig3,isig1,gluon_hel,4,div))

C     Put non-abelian piece here

                  jamp_virt(isig1,isig2,isig3,isig4,gluon_hel,13,div) = 
     $                 CLR(flavor(3),1,isig3)*GHVV(1)*
     $                 CLR(flavor(1),1,isig1)*(-1.0d0)*
     $                 (HEX2143(1,isig1,isig3,gluon_hel,1,div) + 
     $                   HEX2143(1,isig3,isig1,gluon_hel,2,div) +
     $                   HEX2143(1,isig1,isig3,gluon_hel,3,div) +
     $                   HEX2143(1,isig3,isig1,gluon_hel,4,div))


     
               endif
    
cc only needed for same generations/same quarks

                if(abs(CLR(flavor(3),2,isig3)*
     $                 CLR(flavor(1),2,isig1)).gt.0d0) 
     $                 then
                   isig2 = isig3
                   isig4 = isig1

                  jamp_virt(isig1,isig2,isig3,isig4,gluon_hel,9,div) = 
     $                 SgnF*CLR(flavor(3),2,isig3)*GHVV(2)*
     $                 CLR(flavor(1),2,isig1)*(-1.0D0)*
     $                 (PentHex41(3,isig1,isig3,gluon_hel,1,div) + 
     $                  PentHex41(2,isig1,isig3,gluon_hel,1,div) +
     $                  PentHex41(3,isig1,isig3,gluon_hel,2,div) +
     $                  PentHex41(2,isig1,isig3,gluon_hel,2,div) +
     $                  PentHex41(2,isig1,isig3,gluon_hel,3,div) +
     $                  PentHex41(2,isig1,isig3,gluon_hel,4,div)) 
                     
                  jamp_virt(isig1,isig2,isig3,isig4,gluon_hel,10,div) = 
     $                 SgnF*CLR(flavor(3),2,isig3)*GHVV(2)*
     $                 CLR(flavor(1),2,isig1)*(-1.0D0)* 
     $                 (PentHex41(1,isig1,isig3,gluon_hel,1,div) + 
     $                  PentHex41(1,isig1,isig3,gluon_hel,2,div) +
     $                  PentHex41(3,isig1,isig3,gluon_hel,3,div) +
     $                  PentHex41(1,isig1,isig3,gluon_hel,3,div) +
     $                  PentHex41(3,isig1,isig3,gluon_hel,4,div) +
     $                  PentHex41(1,isig1,isig3,gluon_hel,4,div))

                  jamp_virt(isig1,isig2,isig3,isig4,gluon_hel,11,div) = 
     $                 SgnF*CLR(flavor(3),2,isig3)*GHVV(2)*
     $                 CLR(flavor(1),2,isig1)*(-1.0D0)*
     $                 (PentHex23(3,isig3,isig1,gluon_hel,1,div) + 
     $                  PentHex23(2,isig3,isig1,gluon_hel,1,div) +
     $                  PentHex23(3,isig3,isig1,gluon_hel,2,div) +
     $                  PentHex23(2,isig3,isig1,gluon_hel,2,div) +
     $                  PentHex23(2,isig3,isig1,gluon_hel,3,div) +
     $                  PentHex23(2,isig3,isig1,gluon_hel,4,div))
                      
                  jamp_virt(isig1,isig2,isig3,isig4,gluon_hel,12,div) = 
     $                 SgnF*CLR(flavor(3),2,isig3)*GHVV(2)*
     $                 CLR(flavor(1),2,isig1)* (-1.0D0)*
     $                 (PentHex23(1,isig3,isig1,gluon_hel,1,div) + 
     $                  PentHex23(1,isig3,isig1,gluon_hel,2,div) +
     $                  PentHex23(3,isig3,isig1,gluon_hel,3,div) +
     $                  PentHex23(1,isig3,isig1,gluon_hel,3,div) +
     $                  PentHex23(3,isig3,isig1,gluon_hel,4,div) +
     $                  PentHex23(1,isig3,isig1,gluon_hel,4,div))

                  



                    jamp_virt(isig1,isig2,isig3,isig4,gluon_hel,14,div) =
     $                   SgnF*CLR(flavor(3),2,isig3)*GHVV(2)*
     $                   CLR(flavor(1),2,isig1)* (-1.0D0)*
     $                   (HEX4123(1,isig1,isig3,gluon_hel,1,div) + 
     $                   HEX4123(1,isig3,isig1,gluon_hel,2,div) +
     $                   HEX4123(1,isig1,isig3,gluon_hel,3,div) +
     $                   HEX4123(1,isig3,isig1,gluon_hel,4,div))


                  endif
               enddo
            enddo
         enddo
      enddo

      endif

      fac2 = 1.0d0


c     new color basis: I apply Fierz identity to reduce to c1,c2,c3,c4

      do div=0,2
         do isig1=-1,1,2
            do isig2=-1,1,2
               do isig3=-1,1,2
                  do isig4=-1,1,2
                     do i=1,2   ! gluon
 
                        
                        amp_virt(isig1,isig2,isig3,isig4,i,1,div) = 
     $                       jamp_virt(isig1,isig2,isig3,isig4,i,1,div) +  
     $                       fac2*((jamp_virt(isig1,isig2,isig3,isig4,i,9,div)/2.0d0 + 
     $                       jamp_virt(isig1,isig2,isig3,isig4,i,12,div)/2.0d0 - 
     $                       jamp_virt(isig1,isig2,isig3,isig4,i,5,div)/6.0d0 - 
     $                       jamp_virt(isig1,isig2,isig3,isig4,i,6,div)/6.0d0) -                       
     $                       (Ic/2.0d0*jamp_virt(isig1,isig2,isig3,isig4,i,14,div)))
                        
                       

                        amp_virt(isig1,isig2,isig3,isig4,i,2,div) = 
     $                       jamp_virt(isig1,isig2,isig3,isig4,i,2,div) +  
     $                       fac2*((jamp_virt(isig1,isig2,isig3,isig4,i,11,div)/2.0d0 + 
     $                       jamp_virt(isig1,isig2,isig3,isig4,i,10,div)/2.0d0 - 
     $                       jamp_virt(isig1,isig2,isig3,isig4,i,7,div)/6.0d0 - 
     $                       jamp_virt(isig1,isig2,isig3,isig4,i,8,div)/6.0d0) +                      
     $                       (Ic/2.0d0*jamp_virt(isig1,isig2,isig3,isig4,i,14,div)))

                        amp_virt(isig1,isig2,isig3,isig4,i,3,div) = 
     $                       jamp_virt(isig1,isig2,isig3,isig4,i,3,div) +  
     $                       fac2*((jamp_virt(isig1,isig2,isig3,isig4,i,6,div)/2.0d0 + 
     $                       jamp_virt(isig1,isig2,isig3,isig4,i,7,div)/2.0d0 - 
     $                       jamp_virt(isig1,isig2,isig3,isig4,i,11,div)/6.0d0 - 
     $                       jamp_virt(isig1,isig2,isig3,isig4,i,12,div)/6.0d0)+                      
     $                       (Ic/2.0d0*jamp_virt(isig1,isig2,isig3,isig4,i,13,div)))
 
                        amp_virt(isig1,isig2,isig3,isig4,i,4,div) = 
     $                       jamp_virt(isig1,isig2,isig3,isig4,i,4,div) +  
     $                       fac2*((jamp_virt(isig1,isig2,isig3,isig4,i,8,div)/2.0d0 + 
     $                       jamp_virt(isig1,isig2,isig3,isig4,i,5,div)/2.0d0 - 
     $                       jamp_virt(isig1,isig2,isig3,isig4,i,9,div)/6.0d0 - 
     $                       jamp_virt(isig1,isig2,isig3,isig4,i,10,div)/6.0d0)-                       
     $                       (Ic/2.0d0*jamp_virt(isig1,isig2,isig3,isig4,i,13,div)))
                        
ccc

                        if(ldebug) then

                        if(div.eq.0) then
                            print*,"amp_virt1(",isig1,isig2,isig3,isig4,i,div,")",
     $                       amp_virt(isig1,isig2,isig3,isig4,i,1,div)

                        print*,"amp_virt2(",isig1,isig2,isig3,isig4,i,div,")",
     $                       amp_virt(isig1,isig2,isig3,isig4,i,2,div)
                        
                        print*,"amp_virt3(",isig1,isig2,isig3,isig4,i,div,")",
     $                       amp_virt(isig1,isig2,isig3,isig4,i,3,div)
                        
                        print*,"amp_virt4(",isig1,isig2,isig3,isig4,i,div,")",
     $                       amp_virt(isig1,isig2,isig3,isig4,i,4,div)
                        endif
                       
                        endif
                        
                     enddo
                  enddo
               enddo
            enddo
         enddo
      enddo
                     
c
c
      born = 0.0d0 

      do div=0,2
         boxvirt(div) = 0.0d0
         virt(div) = 0.0d0
      enddo


C    Perfom sum over colors and helicities here

C     Perform color sum (2 Re[virt^{dagger} born] for a specific helicity
C     Perfom sum over helicities

C     Virtual piece 

      alphas = 1.0d0

      do div=0,2
         virt(div) = 0.0d0
         
         do gluon_hel = 1,2
c     gluon_hel = 2
            do isig1=-1,1,2
               do isig3 = -1,1,2
                  do isig2=-1,1,2
                     do isig4 = -1,1,2
c      gluon_hel = 1
c      isig1 = -1
c      isig2 = -1
c      isig3 = -1
c      isig4 = -1
                        
                  do icol=1,4
                     matNLO(icol) = 
     $                    amp_virt(isig1,isig2,isig3,isig4,
     $                    gluon_hel,icol,div)/
     $                    (4.0d0*pi)*alphas
                  enddo
                  
                  do icol=1,4
                     matLO(icol) = jamp_born(isig1,isig2,
     $                    isig3,isig4,gluon_hel,icol)
                  enddo
c      matLO(3) = czero
c      matNLO(3) = czero
c      matLO(4) = czero
c      matNLO(4) = czero
c      matNLO(4) = czero
c     matNLO(1) = jamp_virt(isig1,isig3,gluon_hel,1,div)
c     do icol = 5,12
c                     matNLO(icol)=0.0d0
c     enddo
                  
                  virt(div) = virt(div) + 
     $                 ColorSumNLO(matNLO,matLO)
                  
               enddo
            enddo
         enddo
      enddo
      enddo
      enddo


      do div=0,2
         boxvirt(div) = 0.0d0
         
         do gluon_hel = 1,2

            do isig1=-1,1,2
               do isig3 = -1,1,2
                  do isig2=-1,1,2
                     do isig4 = -1,1,2
        
                  do icol=1,4
                     matNLO(icol) = 
     $                    jamp_virt(isig1,isig2,isig3,isig4,
     $                    gluon_hel,icol,div)/
     $                    (4.0d0*pi)*alphas
                  enddo

                    
                  do icol=1,4
                     matLO(icol) = jamp_born(isig1,isig2,
     $                    isig3,isig4,gluon_hel,icol)
                  enddo

                  
                  boxvirt(div) = boxvirt(div) + 
     $                 ColorSumNLO(matNLO,matLO)
                  
               enddo
            enddo
         enddo
      enddo
      enddo
      enddo

C
C     Include QCD renormalization
C


C Born piece

       
      born = 0.0d0
      do gluon_hel = 1,2
c         gluon_hel = 2
         do isig1=-1,1,2
            do isig3 = -1,1,2
               do isig2=-1,1,2
                  do isig4 = -1,1,2
c         isig1 = 1
c         isig2 = 1
c         isig3 = 1
c         isig4 = 1
c               print*,"amp_born=",jamp_born(isig1,isig3,gluon_hel,1)
               do icol=1,4
                  matLO(icol) = jamp_born(isig1,isig2,isig3,isig4,
     $                 gluon_hel,icol)
               enddo
c               matLO(2) = czero
c               matLO(3) = czero
c               matLO(4) = czero
               born = born + 
     $              ColorSumLO(matLO)
              
            enddo
         enddo
      enddo
      enddo
      enddo
C
C     QCD renormalization 
C
      if(compboxes) then
         call QCDRenorm(virt,virt,born,alphas,nf)
         IF(ldebug) then
         call QCDRenorm(boxvirt,boxvirt,born,alphas,nf)
         endif
      endif

C     Return virt(0:2)= 2Re[born*virt] where the index is div=0 (finite), 1 (eps^-1), 2 (eps^-2)
c
c     Convert to CDR 

      call ConvertToCDR(virt,virt,born,alphas)
         IF(ldebug) then
      call ConvertToCDR(boxvirt,boxvirt,born,alphas)
      endif

   
C
C     Physical momentum 

      if(lconvert .or. ldebug1) then
     
      do mu =0,3
         pa(mu) = -1.0d0*p(mu,1) !quark
         pb(mu) = -1.0d0*p(mu,3) !quark
         p1(mu) = p(mu,2) !quark
         p2(mu) = p(mu,4) !quark
         p3(mu) = p(mu,5) !gluon
      enddo


      do div=0,2
         Ieps(div)=0.0d0
      do gluon_hel = 1,2
c     gluon_hel = 2
            
         do isig1=-1,1,2
            do isig3 = -1,1,2
               do isig2=-1,1,2
                  do isig4 = -1,1,2
c     
                     do icol=1,4
                        matLO(icol) = jamp_born(isig1,isig2,
     $                       isig3,isig4,gluon_hel,icol)
                     enddo
c     
                     if(ldebug1) then
                        call CSIeps(pa,pb,p1,p2,p3,matLO,scale,NF,temp,div)
                        Ieps(div) = Ieps(div) + dreal(temp)*alphas/2.0d0/pi
                     endif
                     
c     Convert to Simon's convention (see note in /terrance of svn):
c$$$  if(lconvert) then
c$$$  call ConvertVirt(pa,pb,p1,p2,p3,matLO,scale,
c$$$  $                          temp,div)
c$$$  
c$$$  virt(div) = virt(div)+dreal(temp)*alphas/2.0d0/pi
c$$$  endif
                  enddo
               enddo
            enddo
         enddo
      enddo
      enddo
      endif
     
       IF(ldebug1) then     

c$$$      do div=0,2
c$$$         if(div.eq.0) then 
c$$$            print*,"finite terms"
c$$$            print*,"Virt(0)/born=",
c$$$     $        virt(div)/born,
c$$$     $      "box/born=",boxvirt(div)/born,"born=",born
c$$$            
c$$$         else
c$$$         print*,"Ieps(",div,")/born",Ieps(div)/born,"Virt/born=",
c$$$     $        virt(div)/born,"Virt/Ieps=",
c$$$     $        virt(div)/Ieps(div),"box/Ieps=",boxvirt(div)/Ieps(div),
c$$$     $        "box/born=",boxvirt(div)/born,"born=",born
c$$$         endif
c$$$      enddo

      do div=0,2
         if(div.eq.0) then 
            print*,"Virt(0)/born=",
     $        virt(div)/born
            print*,"born=",born
            
         else
            print*,"Ieps(",div,")/born",Ieps(div)/born
            print*,"Virt/born(",div,")=",virt(div)/born
            print*,"Virt/Ieps(",div,")=",virt(div)/Ieps(div)
         endif
      enddo

      ENDIF

      

      
      do div=0,2
         virt(div) = 2.0d0*pi*virt(div) ! pull out factor of 2pi
         IF(ldebug) then
            boxvirt(div) = 2.0d0*pi*boxvirt(div) 
         ENDIF
      enddo
      
      GAUGEACCURACY1=GAUGEACCURACY(2)

      end


      end ! last line of routine

ccccccc
c
c
c     QCD correction for triangles 
C     A factor of CF will be included outside this routine !
C
Cccccccc
      subroutine VBF_TRI(virt_fac,q2,muRsq)
  
      implicit none
      double complex virt_fac(0:2)
      double precision q2 ! virtuality of gauge boson
      double precision muRsq ! renorm. scale
      double precision minusq2
      double complex MyLOG
      external MyLOG
      double precision PI
      parameter(Pi=3.14159265358979323846264338328d0)
      logical lDimReg
      parameter(lDimReg=.true.)
      
      
C      
C     A factor of Alphas/(4 PI) (4 PI)^{\epsilon}/Gamma[1-\epsilon] has been factored out!
C
      minusq2 = -1.0d0*q2 ! -q2

      virt_fac(2) = dcmplx(-2.0d0,0.0d0) ! 1/eps^2
      virt_fac(1) = dcmplx(-3.0d0,0.0d0) - 2.0d0*MyLOG(muRsq/minusq2) ! 1/eps^1

      if(lDimReg) then
c     regularization: DR      
         virt_fac(0) = dcmplx(-7.0d0,0.0d0) -3.0d0*MyLOG(muRsq/minusq2) - 
     $        MyLOG(muRsq/minusq2)**2 ! finite piece
      else
c     regularization: CDR
         virt_fac(0) = dcmplx(-8.0d0,0.0d0) -3.0d0*MyLOG(muRsq/minusq2) - 
     $        MyLOG(muRsq/minusq2)**2 ! finite piece
      endif




      end
c
c     My natural log
      double complex function MyLOG(x)
      implicit none 
      double precision x
      double precision PI
      parameter(Pi=3.14159265358979323846264338328d0)

      if(x.lt.0d0) then
         MyLOG = dcmplx(DLOG(abs(x)),PI)
      else
         MyLOG = dcmplx(DLOG(x),0.0d0)
      endif


      return
      end

c
c
c
c     Convert result from DR to CDR.
c
      subroutine ConvertToCDR(virtCDR,virtDR,born,alphas)
      implicit none
      double precision alphas,virtCDR(0:2),virtDR(0:2),born
      double precision CF,CA,PI
      parameter(PI=3.14159265358979323846264338328d0)


      CA = 3.0d0
      CF = 4.0d0/3.0d0

      virtCDR(2) = virtDR(2)
      virtCDR(1) = virtDR(1)
      virtCDR(0) = virtDR(0) - alphas/(2.0d0*pi)*
     $     born*(2.0d0*CF + CA/6.0d0)

      end

CCCC
C
C     QCD renormalization factor

      subroutine QCDRenorm(virtBARE,virtRENORM,born,alphas,nf)
      implicit none 
      double precision alphas
      double precision virtBARE(0:2),virtRENORM(0:2),born
      double precision NF,CA,TR ! number of flavors, Casmir in adjoint repres, 
      parameter(!NF=5.0d0,
     $  CA=3.0d0,TR=1.0d0/2.0d0)

      double precision BETA0Tilde,BETA0,Pi
c      parameter(BETA0Tilde=0.0d0) ! MSbar scheme
      parameter(Pi=3.14159265358979323846264338328d0)

      logical lDRED
      parameter(lDRED=.true.)
C      
C     A factor of (4 PI)^{\epsilon}/Gamma[1-\epsilon] has been factore out!
C     This needs to be the case for the BARE 1-loop matrix element in order to 
C     cancel the UV divergence.

      BETA0 = (11.0d0*CA-4.0d0*TR*NF)/6.0d0 ! \beta_{0}

      virtRENORM(2) = virtBARE(2) ! 1/eps^2
      virtRENORM(1) = virtBARE(1) - BETA0 * born * 
     $     alphas/(2.0d0*pi)    !1/eps
      if(lDRED) then 
         BETA0TILDE = -CA/6.0d0
      else
         BETA0TILDE = 0d0 ! use this for MSbar alpha_{s}
      endif
      
       virtRENORM(0) = virtBARE(0) - BETA0Tilde * born * 
     $     alphas/(2.0d0*pi)    ! eps^0
      end




C
C     Computes 2 Re[matNLO^{dagger} matLO] summed over colors
C
      double precision function ColorSumNLO(matNLO,matLO)
      implicit none 
      double complex matNLO(4),matLO(4),result
      double complex mNLO1,mNLO2,mNLO3,mNLO4
      double complex mLO1,mLO2,mLO3,mLO4
      double precision CA,CF

      CA = 3.0d0
      CF = 4.0d0/3.0d0

      mNLO1 = matNLO(1)
      mNLO2 = matNLO(2)
      mNLO3 = matNLO(3)
      mNLO4 = matNLO(4)

      mLO1 = matLO(1)
      mLO2 = matLO(2)

      mLO3 = matLO(3)
      mLO4 = matLO(4)
      

      result =  CA*CF*((CA*mNLO1 + mNLO3 + mNLO4)*DConjg(mLO1) + 
     $     (CA*mNLO2 + mNLO3 + mNLO4)*DConjg(mLO2) + 
     $     mNLO1*DConjg(mLO3) + mNLO2*DConjg(mLO3) + 
     $    CA*mNLO3*DConjg(mLO3) + mNLO1*DConjg(mLO4) + 
     $     mNLO2*DConjg(mLO4) + CA*mNLO4*DConjg(mLO4) + 
     $     CA*mLO1*DConjg(mNLO1) + mLO3*DConjg(mNLO1) + 
     $    mLO4*DConjg(mNLO1) + CA*mLO2*DConjg(mNLO2) + 
     $     mLO3*DConjg(mNLO2) + mLO4*DConjg(mNLO2) + 
     $     mLO1*DConjg(mNLO3) + mLO2*DConjg(mNLO3) + 
     $    CA*mLO3*DConjg(mNLO3) + mLO1*DConjg(mNLO4) + 
     $     mLO2*DConjg(mNLO4) + CA*mLO4*DConjg(mNLO4))
      
      ColorSumNLO = dreal(result) ! 2 Re[MB*MV^{dagger}]


      return

      end

C
C
C     Computes 2 Re[matNLO^{dagger} matLO] summed over colors
C
      double precision function ColorSumLO(matLO)
      implicit none 
      double complex matLO(4),result
      double precision CF,CA
      double complex mLO1,mLO2,mLO3,mLO4

      CA = 3.0d0
      CF = 4.0d0/3.0d0

      
      mLO1 = matLO(1)
      mLO2 = matLO(2)

      mLO3 = matLO(3)
      mLO4 = matLO(4)

      result = CA*CF*((CA*mLO1 + mLO3 + mLO4)*DConjg(mLO1) + (CA*mLO2 + 
     $     mLO3 + mLO4)*DConjg(mLO2) + mLO1*DConjg(mLO3) + 
     $     mLO2*DConjg(mLO3) + 
     $    CA*mLO3*DConjg(mLO3) + mLO1*DConjg(mLO4) + 
     $     mLO2*DConjg(mLO4) + CA*mLO4*DConjg(mLO4))


      
      ColorSumLO = dreal(result) ! |M_B|^2 

c      ColorSumLO = dimag(matLO(1))**2+dreal(matLO(1))**2 + 
c     $     dimag(matLO(2))**2+dreal(matLO(2))**2
      return

      end


      double complex function VPROP(p2,m2,width)
c     We will be working a complex mass scheme.
C     
      implicit none
      double precision p2 ! mom squared
      double precision m2 ! mass squared
      double precision width ! width

c      if (p2.le.0d0) then
c         VPROP = 1/(p2-m2 )
c      else
         VPROP = 1.0d0/dcmplx(p2-m2,dsqrt(m2)*width)

c         VPROP = -1.0d0/((p2-m2)**2 + m2*width**2)*
c     $        dcmplx(m2-p2,dsqrt(m2)*width)
c      endif

      return
      end


      subroutine NOABE_HEX_VBF(Hex,Born,M,p1,PSI1,p2,PSI2,p3,PSI3,p4,PSI4,p5,
     $     eps,p6,musq,IsVaWtemp,divmax,Ward)
  
c     p1(mu) PSI1 q1 
c     p2(mu) PSI2 q2bar
c     p3(mu) PSI3 q3
c     p4(mu) PSI4 q4bar
c     p5(mu) gluon
c     p6(mu) higgs

      implicit none
      integer IsVaWtemp
      integer divmax
      double complex Hex(2,-1:1,-1:1,2,4,0:2)
      double complex Born(-1:1,-1:1,2,4)
      double complex Hexg(2,-1:1,-1:1,2,4,0:2)
      double complex Borng(-1:1,-1:1,2,4)
      integer div,comp,gluon_hel,mu
      double precision musq
      double complex M 
      double complex eps(0:3,2),MyEps(0:3)

      double complex psi_k6(4),psi_k4(4),barpsi_k2(4),barpsi_k1(4)
      double precision k1(0:3),k2(0:3),k3(0:3),k4(0:3),k5(0:3),k6(0:3)
      double precision p1(0:3),p2(0:3),p3(0:3),p4(0:3),p5(0:3),p6(0:3)
      double complex PSI1(2,-1:1),PSI2(2,-1:1),PSI3(2,-1:1),PSI4(2,-1:1)
      double precision temp,temp1
      logical Ward
      integer isig1,isig3,j1,j2
      double precision Gaugethres,gaugeAccuracy(2)
      common/Gauge/gaugethres,gaugeAccuracy
      ward=.false.
c
      do mu=0,3
         k1(mu) = -p4(mu)       ! psibar 
         k2(mu) = -p2(mu)       ! psibar
         k3(mu) = -p5(mu)       ! gluon
         k4(mu) = p1(mu)        ! psi
         k5(mu) = -p6(mu)       ! higgs
         k6(mu) = p3(mu)        ! psi
         MyEps(mu)= k3(mu)
      enddo

c     
      barpsi_k1(1) = PSI4(1,1)  !P
      barpsi_k1(2) = PSI4(2,1)!P
      barpsi_k1(3) = PSI4(1,-1)!M
      barpsi_k1(4) = PSI4(2,-1)!M


      barpsi_k2(1) = PSI2(1,1)!P
      barpsi_k2(2) = PSI2(2,1)!P
      barpsi_k2(3) = PSI2(1,-1)!M
      barpsi_k2(4) = PSI2(2,-1)!M


      psi_k6(1) = PSI3(1,-1)!M
      psi_k6(2) = PSI3(2,-1)!M
      psi_k6(3) = PSI3(1,1)!P
      psi_k6(4) = PSI3(2,1)!P


      psi_k4(1) = PSI1(1,-1)!M
      psi_k4(2) = PSI1(2,-1)!M
      psi_k4(3) = PSI1(1,1)!P
      psi_k4(4) = PSI1(2,1)!P


      do div=0,divMax
         comp = 1
         do gluon_hel = 1,2
            call NoAbeH3j77T_c(M,k1,k3,k2,k4,k5,k6,barpsi_k1,psi_k6,barpsi_k2,
     $           psi_k4,eps(0,gluon_hel),musq,comp,
     $           Hex(1,-1,-1,gluon_hel,1,div), 
     $           Born(-1,-1,gluon_hel,1),div)

            
c For Gauge Test
           if(div.eq.0) then
            if(gluon_hel.eq.1) then
 
            call NoAbeH3j77T_c(M,k1,k3,k2,k4,k5,k6,barpsi_k1,psi_k6,barpsi_k2,
     $           psi_k4,Myeps,musq,comp,
     $           Hexg(1,-1,-1,1,1,0), 
     $           Borng(-1,-1,1,1),0)
           
            temp= abs(Hexg(1,-1,-1,1,1,0)/Hex(2,-1,-1,1,1,0)-1d0)
            
c            print*,"temp",temp
c            stop
           endif
            endif
c Finish gauge test


            comp = -1

            if(IsVaWtemp.eq.0) then

            call NoAbeH3j67T_c(M,k1,k3,k2,k4,k5,k6,barpsi_k1,psi_k6,barpsi_k2,
     $           psi_k4,eps(0,gluon_hel),musq,comp,
     $           Hex(1,-1,1,gluon_hel,1,div),
     $           Born(-1,1,gluon_hel,1),div)
            
            call NoAbeH3j76T_c(M,k1,k3,k2,k4,k5,k6,barpsi_k1,psi_k6,barpsi_k2,
     $           psi_k4,eps(0,gluon_hel),musq,comp,
     $           Hex(1,1,-1,gluon_hel,1,div),
     $           Born(1,-1,gluon_hel,1),div)
            
            call NoAbeH3j66T_c(M,k1,k3,k2,k4,k5,k6,barpsi_k1,psi_k6,barpsi_k2,
     $           psi_k4,eps(0,gluon_hel),musq,comp,
     $           Hex(1,1,1,gluon_hel,1,div),
     $           Born(1,1,gluon_hel,1),div)

            endif
         enddo
      enddo

      
      do div=0,divmax
         comp = 1
         do gluon_hel = 1,2
            call NoAbeH3jCross77T_c(M,k6,k3,k2,k4,k5,k1,psi_k6,barpsi_k1,barpsi_k2,
     $           psi_k4,eps(0,gluon_hel),musq,comp,
     $           Hex(1,-1,-1,gluon_hel,2,div),
     $           Born(-1,-1,gluon_hel,2),div)

            comp = -1

c For Gauge Test
           if(div.eq.0) then
            if(gluon_hel.eq.1) then
 
            call NoAbeH3jCross77T_c(M,k6,k3,k2,k4,k5,k1,psi_k6,barpsi_k1,barpsi_k2,
     $           psi_k4,MyEps,musq,comp,
     $           Hexg(1,-1,-1,1,2,0),
     $           Borng(-1,-1,1,2),0)

            temp1= abs(Hexg(1,-1,-1,1,2,0)/Hex(2,-1,-1,1,2,0)-1d0)
            if(temp1.gt.temp) temp=temp1 
c            print*,"tempCross",temp
c            stop
           endif
            endif
c Finish gauge test



         if(IsVaWtemp.eq.0) then

            call NoAbeH3jCross67T_c(M,k6,k3,k2,k4,k5,k1,psi_k6,barpsi_k1,barpsi_k2,
     $           psi_k4,eps(0,gluon_hel),musq,comp,
     $           Hex(1,-1,1,gluon_hel,2,div),
     $           Born(-1,1,gluon_hel,2),div)
            
            call NoAbeH3jCross76T_c(M,k6,k3,k2,k4,k5,k1,psi_k6,barpsi_k1,barpsi_k2,
     $           psi_k4,eps(0,gluon_hel),musq,comp,
     $           Hex(1,1,-1,gluon_hel,2,div),
     $           Born(1,-1,gluon_hel,2),div)
            
            call NoAbeH3jCross66T_c(M,k6,k3,k2,k4,k5,k1,psi_k6,barpsi_k1,barpsi_k2,
     $           psi_k4,eps(0,gluon_hel),musq,comp,
     $           Hex(1,1,1,gluon_hel,2,div),
     $           Born(1,1,gluon_hel,2),div)

            endif

         enddo
      enddo


      do div=0,divmax
         comp = 1
         do gluon_hel = 1,2
            call NoAbeH3jCrossF77T_c(M,k1,k3,k4,k2,k5,k6,barpsi_k1,psi_k6,psi_k4,
     $           barpsi_k2,eps(0,gluon_hel),musq,comp,
     $           Hex(1,-1,-1,gluon_hel,3,div),
     $           Born(-1,-1,gluon_hel,3),div)

            comp = -1

c For Gauge Test
           if(div.eq.0) then
            if(gluon_hel.eq.1) then
 
            call  NoAbeH3jCrossF77T_c(M,k1,k3,k4,k2,k5,k6,barpsi_k1,psi_k6,psi_k4,
     $           barpsi_k2,MyEps,musq,comp,
     $           Hexg(1,-1,-1,1,3,0),
     $           Borng(-1,-1,1,3),0)

            temp1= abs(Hexg(1,-1,-1,1,3,0)/Hex(2,-1,-1,1,3,0)-1d0)
            if(temp1.gt.temp) temp=temp1 
c            print*,"tempCrossF",temp
c            stop
           endif
            endif
c Finish gauge test



            if(IsVaWtemp.eq.0) then
            
            call NoAbeH3jCrossF67T_c(M,k1,k3,k4,k2,k5,k6,barpsi_k1,psi_k6,psi_k4,
     $           barpsi_k2,eps(0,gluon_hel),musq,comp,
     $           Hex(1,-1,1,gluon_hel,3,div),
     $           Born(-1,1,gluon_hel,3),div)
            
            call NoAbeH3jCrossF76T_c(M,k1,k3,k4,k2,k5,k6,barpsi_k1,psi_k6,psi_k4,
     $           barpsi_k2,eps(0,gluon_hel),musq,comp,
     $           Hex(1,1,-1,gluon_hel,3,div),
     $           Born(1,-1,gluon_hel,3),div)
            
            call NoAbeH3jCrossF66T_c(M,k1,k3,k4,k2,k5,k6,barpsi_k1,psi_k6,psi_k4,
     $           barpsi_k2,eps(0,gluon_hel),musq,comp,
     $           Hex(1,1,1,gluon_hel,3,div),
     $           Born(1,1,gluon_hel,3),div)

            endif
         enddo
      enddo

     

      do div=0,divmax
         comp = 1
         do gluon_hel = 1,2
            call NoAbeH3jCrossIF77T_c(M,k6,k3,k4,k2,k5,k1,psi_k6,barpsi_k1,psi_k4,
     $           barpsi_k2,eps(0,gluon_hel),musq,comp,
     $           Hex(1,-1,-1,gluon_hel,4,div),
     $           Born(-1,-1,gluon_hel,4),div)

            comp = -1

c For Gauge Test
           if(div.eq.0) then
            if(gluon_hel.eq.1) then
 
            call  NoAbeH3jCrossIF77T_c(M,k6,k3,k4,k2,k5,k1,psi_k6,barpsi_k1,psi_k4,
     $           barpsi_k2,MyEps,musq,comp,
     $           Hexg(1,-1,-1,gluon_hel,4,div),
     $           Borng(-1,-1,gluon_hel,4),div)

            temp1= abs(Hexg(1,-1,-1,1,4,0)/Hex(2,-1,-1,1,4,0)-1d0)
            if(temp1.gt.temp) temp=temp1 
c            print*,"tempCrossIF",temp
c            stop
           endif
            endif
c Finish gauge test


          if(IsVaWtemp.eq.0) then      

            call NoAbeH3jCrossIF67T_c(M,k6,k3,k4,k2,k5,k1,psi_k6,barpsi_k1,psi_k4,
     $           barpsi_k2,eps(0,gluon_hel),musq,comp,
     $           Hex(1,-1,1,gluon_hel,4,div),
     $           Born(-1,1,gluon_hel,4),div)
            
            call NoAbeH3jCrossIF76T_c(M,k6,k3,k4,k2,k5,k1,psi_k6,barpsi_k1,psi_k4,
     $           barpsi_k2,eps(0,gluon_hel),musq,comp,
     $           Hex(1,1,-1,gluon_hel,4,div),
     $           Born(1,-1,gluon_hel,4),div)
            
            call NoAbeH3jCrossIF66T_c(M,k6,k3,k4,k2,k5,k1,psi_k6,barpsi_k1,psi_k4,
     $           barpsi_k2,eps(0,gluon_hel),musq,comp,
     $           Hex(1,1,1,gluon_hel,4,div),
     $           Born(1,1,gluon_hel,4),div)

            endif

         enddo
      enddo

      if(temp.lt.gaugethres) then
       Ward=.true.

       gaugeAccuracy(1)=temp

      if(gaugeAccuracy(1).gt.gaugeAccuracy(2)) then
         gaugeAccuracy(2)=gaugeAccuracy(1)
      endif       

c$$$      else
c$$$
c$$$      do div=0,divmax
c$$$           do gluon_hel = 1,2
c$$$              do isig1=-1,1,2
c$$$                 do isig3=-1,1,2
c$$$                    do j1=1,4                       
c$$$                   Hex(1,isig1,isig3,gluon_hel,j1,div)=0d0
c$$$                   Hex(2,isig1,isig3,gluon_hel,j1,div)=0d0
c$$$                   Born(isig1,isig3,gluon_hel,j1)=0d0
c$$$                    enddo
c$$$               enddo
c$$$              enddo
c$$$         enddo
c$$$      enddo
      endif



c      print*,"temp",temp, Ward
       
      end



      subroutine PENTHEX_VBF(PentHex,BORN,M,p1,PSI1,p2,PSI2,p3,PSI3,
     $     p4,PSI4,
     $     p5,eps,p6,musq,IsVaWtemp,divmax,Ward)

cFC
cFC      subroutine PENTHEX_VBF(PentHex,M,p1,PSI1,p2,PSI2,p3,PSI3,p4,PSI4,
cFC     $     p5,eps,p6,musq,Zboson)
CFC finish 
C
C     Emission off 21 line 
C
c     p1(mu) PSI1 q1 
c     p2(mu) PSI2 q2bar
c     p3(mu) PSI3 q3
c     p4(mu) PSI4 q4bar
c     p5(mu) gluon
c     p6(mu) higgs

      implicit none
      integer IsVaWtemp
      integer divmax
      double complex PentHex(3,-1:1,-1:1,2,4,0:2)
      double complex PentHexg(3,-1:1,-1:1,2,4,0:2)
      double complex Borng(-1:1,-1:1,2,4)
      integer div,comp,mu,gluon_hel
      double precision musq
      double complex M 
      double complex eps(0:3,2)

      double complex psi_k6(4),psi_k4(4),barpsi_k2(4),barpsi_k1(4)
      double precision k1(0:3),k2(0:3),k3(0:3),k4(0:3),k5(0:3),k6(0:3)
      double precision p1(0:3),p2(0:3),p3(0:3),p4(0:3),p5(0:3),p6(0:3)
      double complex PSI1(2,-1:1),PSI2(2,-1:1),PSI3(2,-1:1),PSI4(2,-1:1) 
      double complex Born(-1:1,-1:1,2,4)
      double complex MyEps(0:3)
      double precision temp,temp1
      logical Ward
      integer j1,j2,isig1,isig3
      double precision Gaugethres,gaugeAccuracy(2)
      common/Gauge/gaugethres,gaugeAccuracy
      ward=.false.

C     define inflowing momentum (k1+k2+ ... + k6 = 0)
C     interface to Paco's code
      do mu=0,3
         k1(mu) = -p4(mu)! psibar 
         k2(mu) = -p2(mu)! psibar
         k3(mu) = -p5(mu)! gluon
         k4(mu) = p1(mu) ! psi
         k5(mu) = -p6(mu)! higgs
         k6(mu) = p3(mu) ! psi
         MyEps(mu) = k3(mu)
      enddo

C
C calculate bar and ket or psi and psibar spinors
C psi is a 4-spinor which will be constructed frow two 2-spinors.
C psibar is a 4-spinor which will be constructed from two 2-spinors.
C
C     define psibar_k1,psibar_k2, psi_k4,psi_k6
C
C build 4-spinors
      barpsi_k1(1) = PSI4(1,1)!P
      barpsi_k1(2) = PSI4(2,1)!P
      barpsi_k1(3) = PSI4(1,-1)!M
      barpsi_k1(4) = PSI4(2,-1)!M


      barpsi_k2(1) = PSI2(1,1)!P
      barpsi_k2(2) = PSI2(2,1)!P
      barpsi_k2(3) = PSI2(1,-1)!M
      barpsi_k2(4) = PSI2(2,-1)!M


      psi_k6(1) = PSI3(1,-1)!M
      psi_k6(2) = PSI3(2,-1)!M
      psi_k6(3) = PSI3(1,1)!P
      psi_k6(4) = PSI3(2,1)!P


      psi_k4(1) = PSI1(1,-1)!M
      psi_k4(2) = PSI1(2,-1)!M
      psi_k4(3) = PSI1(1,1)!P
      psi_k4(4) = PSI1(2,1)!P

c      print*, " inside PENT HEX ABE 2" 

      do div=0,divmax
c FC only computed once for each divergences
            comp=1
c Fc finish
         do gluon_hel=1,2
c            comp=1
            call H3j77T_c(M,k1,k2,k3,k4,k5,k6,barpsi_k1,psi_k6,barpsi_k2,
     $           psi_k4,eps(0,gluon_hel),musq,comp,
     $           PentHex(1,-1,-1,gluon_hel,1,div),
     $           Born(-1,-1,gluon_hel,1),Div)
            comp=-1

c For Gauge Test
           if(div.eq.0) then
            if(gluon_hel.eq.1) then
            call H3j77T_c(M,k1,k2,k3,k4,k5,k6,barpsi_k1,psi_k6,barpsi_k2,
     $           psi_k4,Myeps,musq,comp,
     $           PentHexg(1,-1,-1,1,1,0),
     $           Borng(-1,-1,1,1),0)

            call HexAux(PentHexg(1,-1,-1,1,1,0),temp)
c            print*, "temp", temp

           endif
            endif
c Finish gauge test



            if(IsVaWtemp.eq.0) then         
            
            call H3j76T_c(M,k1,k2,k3,k4,k5,k6,barpsi_k1,psi_k6,barpsi_k2,
     &           psi_k4,eps(0,gluon_hel),musq,comp,
     $           PentHex(1,1,-1,gluon_hel,1,div),
     $           Born(1,-1,gluon_hel,1),Div)

            call H3j67T_c(M,k1,k2,k3,k4,k5,k6,barpsi_k1,psi_k6,barpsi_k2,
     &           psi_k4,eps(0,gluon_hel),musq,comp,
     $           PentHex(1,-1,1,gluon_hel,1,div),
     $           Born(-1,1,gluon_hel,1),Div)

            call H3j66T_c(M,k1,k2,k3,k4,k5,k6,barpsi_k1,psi_k6,barpsi_k2,
     &           psi_k4,eps(0,gluon_hel),musq,comp,
     $           PentHex(1,1,1,gluon_hel,1,div),
     $           Born(1,1,gluon_hel,1),Div)

          endif

         enddo
      enddo

c      print*, " inside PENT HEX ABE 3" 

      do div=0,divmax
c FC only computed once for each divergences
            comp=1
c Fc finish
         do gluon_hel=1,2
c            comp=1
            call H3jCross77T_c(M,k6,k2,k3,k4,k5,k1,psi_k6,barpsi_k1,barpsi_k2,
     &           psi_k4,eps(0,gluon_hel),musq,comp,
     $           PentHex(1,-1,-1,gluon_hel,2,div),
     $           Born(-1,-1,gluon_hel,2),Div)
            
            comp=-1
     
c For Gauge Test
           if(div.eq.0) then
            if(gluon_hel.eq.1) then

            call H3jCross77T_c(M,k6,k2,k3,k4,k5,k1,psi_k6,barpsi_k1,barpsi_k2,
     &           psi_k4,Myeps,musq,comp,
     $           PentHexg(1,-1,-1,1,2,0),
     $           Borng(-1,-1,1,2),0)
 
            call HexAux(PentHexg(1,-1,-1,1,1,0),temp1)
            if(temp1.gt.temp) temp=temp1 
c           print*, "tempCross", temp

           endif
            endif
c Finish gauge test


            if(IsVaWtemp.eq.0) then

            call H3jCross76T_c(M,k6,k2,k3,k4,k5,k1,psi_k6,barpsi_k1,barpsi_k2,
     &           psi_k4,eps(0,gluon_hel),musq,comp,
     $           PentHex(1,1,-1,gluon_hel,2,div),
     $           Born(1,-1,gluon_hel,2),Div)
            
            call H3jCross67T_c(M,k6,k2,k3,k4,k5,k1,psi_k6,barpsi_k1,barpsi_k2,
     &           psi_k4,eps(0,gluon_hel),musq,comp,
     $           PentHex(1,-1,1,gluon_hel,2,div),
     $           Born(-1,1,gluon_hel,2),Div)
            
            call H3jCross66T_c(M,k6,k2,k3,k4,k5,k1,psi_k6,barpsi_k1,barpsi_k2,
     &           psi_k4,eps(0,gluon_hel),musq,comp,
     $           PentHex(1,1,1,gluon_hel,2,div),
     $           Born(1,1,gluon_hel,2),Div)

          endif

         enddo
      enddo

c      print*, " inside PENT HEX ABE 4" 
      do div=0,divmax
c FC only computed once for each divergences
            comp=1
c Fc finish
         do gluon_hel=1,2
c            comp=1
            call H3jCrossF77T_c(M,k1,k4,k3,k2,k5,k6,barpsi_k1,psi_k6,
     &           psi_k4,barpsi_k2,eps(0,gluon_hel),musq,comp,
     $           PentHex(1,-1,-1,gluon_hel,3,div),
     $           Born(-1,-1,gluon_hel,3),Div)

            comp=-1

c For Gauge Test
           if(div.eq.0) then
            if(gluon_hel.eq.1) then

            call H3jCrossF77T_c(M,k1,k4,k3,k2,k5,k6,barpsi_k1,psi_k6,
     &           psi_k4,barpsi_k2,Myeps,musq,comp,
     $           PentHexg(1,-1,-1,1,3,0),
     $           Borng(-1,-1,1,3),0)

            call HexAux(PentHexg(1,-1,-1,1,1,0),temp1)
            if(temp1.gt.temp) temp=temp1 
c           print*, "tempCrossF", temp

           endif
            endif
c Finish gauge test


            if(IsVaWtemp.eq.0) then

            call H3jCrossF76T_c(M,k1,k4,k3,k2,k5,k6,barpsi_k1,psi_k6,
     &           psi_k4,barpsi_k2,eps(0,gluon_hel),musq,comp,
     $           PentHex(1,1,-1,gluon_hel,3,div),
     $           Born(1,-1,gluon_hel,3),Div)

            call H3jCrossF67T_c(M,k1,k4,k3,k2,k5,k6,barpsi_k1,psi_k6,
     &           psi_k4,barpsi_k2,eps(0,gluon_hel),musq,comp,
     $           PentHex(1,-1,1,gluon_hel,3,div),
     $           Born(-1,1,gluon_hel,3),Div)

            call H3jCrossF66T_c(M,k1,k4,k3,k2,k5,k6,barpsi_k1,psi_k6,
     &           psi_k4,barpsi_k2,eps(0,gluon_hel),musq,comp,
     $           PentHex(1,1,1,gluon_hel,3,div),
     $           Born(1,1,gluon_hel,3),Div)

          endif
          
         enddo
      enddo
     
c      print*, " inside PENT HEX  ABE 5" 
      do div = 0,divmax
c FC only computed once for each divergences
            comp=1
c Fc finish
         do gluon_hel = 1,2
c         
            call H3jCrossIF77T_c(M,k6,k4,k3,k2,k5,k1,psi_k6,barpsi_k1,
     &           psi_k4,barpsi_k2,eps(0,gluon_hel),musq,comp,
     $           PentHex(1,-1,-1,gluon_hel,4,div),
     $           Born(-1,-1,gluon_hel,4),Div)
             comp = -1


c For Gauge Test
           if(div.eq.0) then
            if(gluon_hel.eq.1) then

            call H3jCrossIF77T_c(M,k6,k4,k3,k2,k5,k1,psi_k6,barpsi_k1,
     &           psi_k4,barpsi_k2,Myeps,musq,comp,
     $           PentHexg(1,-1,-1,1,4,0),
     $           Borng(-1,-1,1,4),0)

            call HexAux(PentHexg(1,-1,-1,1,1,0),temp1)
            if(temp1.gt.temp) temp=temp1 
c           print*, "tempCrossIF", temp

           endif
            endif
c Finish gauge test

            if(IsVaWtemp.eq.0) then

            call H3jCrossIF76T_c(M,k6,k4,k3,k2,k5,k1,psi_k6,barpsi_k1,
     &           psi_k4,barpsi_k2,eps(0,gluon_hel),musq,comp,
     $           PentHex(1,1,-1,gluon_hel,4,div),
     $           Born(1,-1,gluon_hel,4),Div)

            call H3jCrossIF67T_c(M,k6,k4,k3,k2,k5,k1,psi_k6,barpsi_k1,
     &           psi_k4,barpsi_k2,eps(0,gluon_hel),musq,comp,
     $           PentHex(1,-1,1,gluon_hel,4,div),
     $           Born(-1,1,gluon_hel,4),Div)

            call H3jCrossIF66T_c(M,k6,k4,k3,k2,k5,k1,psi_k6,barpsi_k1,
     &           psi_k4,barpsi_k2,eps(0,gluon_hel),musq,comp,
     $           PentHex(1,1,1,gluon_hel,4,div),
     $           Born(1,1,gluon_hel,4),Div)
          

cFC
          endif
CFC            
         enddo
      enddo
      if(temp.lt.gaugethres) then
         Ward=.true.
       gaugeAccuracy(1)=temp

      if(gaugeAccuracy(1).gt.gaugeAccuracy(2)) then
         gaugeAccuracy(2)=gaugeAccuracy(1)
      endif  


c$$$      else
c$$$         
c$$$         do div=0,divmax
c$$$            do gluon_hel = 1,2
c$$$               do isig1=-1,1,2
c$$$                  do isig3=-1,1,2
c$$$                     do j1=1,4                       
c$$$                        PentHex(1,isig1,isig3,gluon_hel,j1,div)=0d0
c$$$                        PentHex(2,isig1,isig3,gluon_hel,j1,div)=0d0
c$$$                        PentHex(3,isig1,isig3,gluon_hel,j1,div)=0d0
c$$$                        Born(isig1,isig3,gluon_hel,j1)=0d0
c$$$                     enddo
c$$$                  enddo
c$$$               enddo
c$$$            enddo
c$$$         enddo
      endif
c     print*,"temp",temp, Ward      
c     print*, " inside PENT HEX ABE 6" 
      
      end


      subroutine BOX_VBF(psi2,psi1,p1,p2,p21g1,pg,EPS,j43,IsVaWtemp,divmax,scale,
     $     mat21g,mat21,mat21gb,mat21b,NOABEmat21,NOABEmat21g,Ward1)
C     the purpose of this routine is to compute the needed pieces for all
C     triangle and box diagrams

      implicit none
      integer IsVaWtemp
      integer divMax
      logical Ward1,Ward2
      integer isig1F,isig3F
      double precision scale
      integer isig1,isig2,isig3,div,comp,gluon_hel,i
      double precision pg(0:3)! gluon 4-mom
      double precision p1(0:3),p2(0:3)
      double precision p21g(0:3),p21g1(0:4)
      double complex PSI1(2,-1:1),PSI2(2,-1:1),EPS(0:3,2)
      double complex mat21(5,-1:1,-1:1,2,3,0:2), 
     $     mat21g(2,-1:1,-1:1,2,3,0:2), 
     $     mat21gb(2,-1:1,-1:1,2,3,0:2), 
     $     mat21b(-1:1,-1:1,2,3,0:2)
      double complex NOABEmat21(-1:1,-1:1,2,0:2),
     $     NOABEmat21g(2,-1:1,-1:1,2,0:2)
      double complex j43(0:3,-1:1)
      external dotrr
      real*8 dotrr
      double precision Gaugethres,gaugeAccuracy(2)
      common/Gauge/gaugethres,gaugeAccuracy

      Ward1=.false.

      do i=0,3
         p21g(i)=p21g1(i)
      enddo



      if(IsVaWtemp.eq.0) then 
      isig1F=1
      isig3F=1
      else
      isig1F=-1
      isig3F=-1
      endif

      
      do div = 0,divmax
         comp=1
         do gluon_hel = 1,2
            do isig1=-1,isig1F,2
               isig2 = isig1    ! kronecker delta_{isig1,isig2}
               do isig3 = -1,isig3F,2
c     g V
c     
                  call boxlineABETotal(p1,pg,p21g,p2,
     $                 psi2(1,isig2),psi1(1,isig1),
     $                 eps(0,gluon_hel),j43(0,isig3),isig1,
     $                 scale,1,2,
c     FC
c     FC Here the comp variable, also the gauge check should be computed in the first run
                  
     $                 comp*3,comp,
c     FC
     $                 mat21g(1,isig1,isig3,gluon_hel,1,div),
     $                 mat21(1,isig1,isig3,gluon_hel,1,div),
     $                 mat21gb(1,isig1,isig3,gluon_hel,1,div),
     $                 mat21b(isig1,isig3,gluon_hel,1,div),div)

 
     
c     FCC 
                  comp=-1
c     FC
 
               enddo
            enddo
         enddo
      enddo




cFCC
      
      do div = 0,divmax
c         print*,"div=",div
         comp=1
         do gluon_hel = 1,2
            do isig1=-1,isig1F,2
               isig2 = isig1    ! kroncker delta_{isig1,isig2}
               do isig3 = -1,isig3F,2
                  
                  
c     V g
c     
                  call boxlineABETotal(p1,p21g,pg,p2,
     $                 psi2(1,isig2),psi1(1,isig1),
     $                 j43(0,isig3),eps(0,gluon_hel),isig1,
     $                 scale,1,3,
c     FC
c     FC Here the comp variable, also the gauge check should be computed in the first run
     $                 comp*3,comp,
c     FC                 
     $                 mat21g(1,isig1,isig3,gluon_hel,2,div),
     $                 mat21(1,isig1,isig3,gluon_hel,2,div),
     $                 mat21gb(1,isig1,isig3,gluon_hel,2,div),
     $                 mat21b(isig1,isig3,gluon_hel,2,div),div)

c                  print*,"mat2-1b(2)=",mat21b(isig1,isig3,gluon_hel,2,div)
c     
c     Non-abelian
                  
c     FC This can go in the same loop because have different common blocks.
c     FC so the infromation is not overwritten and we can use the trick 
c     FC of the "comp" flag 
                  
                  
                  call boxlineNoAbeTotal(p1,p21g,p2,pg,psi2(1,isig2),
     $                 psi1(1,isig1),
     $                 j43(0,isig3),eps(0,gluon_hel),isig1,scale,
c     FC
c     FC Here the comp variable, also the gauge check should be computed in the first run
                  
     $                 comp*3,comp,
c     FC
     $                 NOABEmat21g(1,isig1,isig3,gluon_hel,Div),
     $                 NOABEmat21(isig1,isig3,gluon_hel,Div),Div)

 
                  
cFC                  
                  comp=-1
cFC
               enddo
            enddo
         enddo
      enddo

      call box_gaugeOneLoop(p1,pg,p21g,p2,scale,
     &     mat21g(1,-1,-1,1,1,0),
     &     mat21gb(1,-1,-1,1,1,0),
     &     Ward1)
      
      if(gaugeAccuracy(1).gt.gaugeAccuracy(2)) then
         gaugeAccuracy(2)=gaugeAccuracy(1)
      endif


      call box_gaugeOneLoop(p1,p21g,pg,p2,scale,
     $     mat21g(1,-1,-1,1,2,0),
     $     mat21gb(1,-1,-1,1,2,0),
     $     Ward2)
      
      if(gaugeAccuracy(1).gt.gaugeAccuracy(2)) then
         gaugeAccuracy(2)=gaugeAccuracy(1)
      endif
     

      Ward1=Ward1.and.Ward2



c     PRINT*, "WARD2",WARD1
      
c$$$      if(.not.Ward1) then
c$$$         
c$$$         do div = 0,divmax
c$$$            do gluon_hel = 1,2
c$$$               do isig1=-1,isig1F,2
c$$$                  isig2 = isig1 ! kroncker delta_{isig1,isig2}
c$$$                  do isig3 = -1,isig3F,2
c$$$                     
c$$$                     NOABEmat21(isig1,isig3,gluon_hel,Div)=0
c$$$                     mat21g(1,isig1,isig3,gluon_hel,2,div)=0
c$$$                     mat21(1,isig1,isig3,gluon_hel,2,div)=0
c$$$                     mat21gb(1,isig1,isig3,gluon_hel,2,div)=0
c$$$                     mat21b(isig1,isig3,gluon_hel,2,div)=0
c$$$                     mat21g(1,isig1,isig3,gluon_hel,1,div)=0
c$$$                     mat21(1,isig1,isig3,gluon_hel,1,div)=0
c$$$                     mat21gb(1,isig1,isig3,gluon_hel,1,div)=0
c$$$                     mat21b(isig1,isig3,gluon_hel,1,div)=0
c$$$                     
c$$$                  enddo
c$$$               enddo
c$$$            enddo
c$$$         enddo
c$$$         
c$$$         
c$$$      endif
      
      end
C
C*********** polarisation vector (POLVEC) ************
C*********** entry HELVEC ****************************
C
C   Calculates the polarisation vector eps_mu(k,lambda) as given in
C   MAD/PH/402, Appendix A, Eqn A.12 and A.13
C   checked on May 24, 88
C   bug in HELVEC for KT=0 corrected Nov.4, 1993
C
C   INPUT:
C
C            kbar      real*8  array(0:4)
C                        physical four momentum for boson and mass
C                        kbar(0:4) = [Ebar, kbarx, kbary, kbarz, mass]
C
C            sigmak      integer
C                        sign factor for vector boson
C                        = +1 outgoing
C                        = -1 incoming
C
C            lambda      integer
C                        rectangular polarization index as used
C                        in Eqn A.12. Allowed values= 1,2,3.
C                        or vector boson helicity (allowed values
C                        -1,0,1) when entry HELVEC is used
C
C
C   OUTPUT:
C            epscar  real*8  array(0:3)
C                        polarisation vector epsilon_mu(lambda)
C                               in eq. A.12 (cartesian basis)
C
C            epshel  complex*16  array(0:3)
C                        polarisation vector epsilon_mu(lambda)
C                               in eq. A.12/13 (helicity basis)
C
C                  The polarization vector for helicity -1,1 or
C                  polarization index 1,2 has vanishing time component.
C                  For longitudinal polarisation, current conservation
C                  is used to eliminate the time component. This is
C                  only possible for MASSLESS FERMIONS. For massive
C                  fermions the polarization vectors defined in eq.A.12
C                  must be used!
C
      SUBROUTINE POLVEC(KBAR,LAMBDA,EPSCAR)
C
C   arguments:
C
      double precision  KBAR(0:4)
      double precision  EPSCAR(0:3)
      double complex EPSHEL(0:3)
      INTEGER SIGMAK, LAMBDA
C
C   local variables:
C
      double precision  A1,A2,A3,B1,B2
      double precision  KT,KMAGNI,NORMAL
C
C   code:
C
C
C            Compute polarization vector:
C            Eqn A.12
C
      EPSCAR(0) = 0.0d0
      IF (LAMBDA .EQ. 1) THEN
         KT = SQRT(KBAR(1)**2+KBAR(2)**2)
         IF (KT .GT. 0) THEN
            IF (KBAR(4) .EQ. 0 ) THEN
               KMAGNI = KBAR(0)
            ELSE
               KMAGNI = SQRT(KBAR(1)**2+KBAR(2)**2+KBAR(3)**2)
            ENDIF
            NORMAL = 1.0d0 / (KMAGNI*KT)
            EPSCAR(1) = KBAR(1)*KBAR(3)*NORMAL
            EPSCAR(2) = KBAR(2)*KBAR(3)*NORMAL
            EPSCAR(3) =          -KT**2*NORMAL
         ELSE
            EPSCAR(1) = dsign(1.0d0, KBAR(3))
            EPSCAR(2) = 0.0d0
            EPSCAR(3) = 0.0d0
         ENDIF
      ELSE IF (LAMBDA .EQ. 2) THEN
         KT = SQRT(KBAR(1)**2+KBAR(2)**2)
         IF (KT .GT. 0) THEN
            EPSCAR(1) = -KBAR(2)/KT
            EPSCAR(2) =  KBAR(1)/KT
            EPSCAR(3) = 0.0d0
         ELSE
            EPSCAR(1) = 0.0d0
            EPSCAR(2) = 1.0d0
            EPSCAR(3) = 0.0d0
         ENDIF
      ELSE IF (LAMBDA .EQ. 3) THEN
         IF (KBAR(4) .EQ. 0) THEN
C
C                  Mass = 0 not allowed for longitudinal polarization
C
            WRITE(*,*) "Mass = 0 for Lambda = 3 in POLVEC"
            EPSCAR(1) = 0.0d0
            EPSCAR(2) = 0.0d0
            EPSCAR(3) = 0.0d0
            RETURN
         END IF
         KMAGNI = SQRT(KBAR(1)**2+KBAR(2)**2+KBAR(3)**2)
         NORMAL = KBAR(4)/(KBAR(0)*KMAGNI)
         EPSCAR(1) = KBAR(1)*NORMAL
         EPSCAR(2) = KBAR(2)*NORMAL
         EPSCAR(3) = KBAR(3)*NORMAL
      ELSE
C
C                  Unrecognised value of Lambda
C
         WRITE (*,*) "Invalid Lambda in POLVEC: Lambda = ",LAMBDA
         EPSCAR(1) = 0.0d0
         EPSCAR(2) = 0.0d0
         EPSCAR(3) = 0.0d0
      END IF
      RETURN
C
      ENTRY HELVEC(KBAR,SIGMAK,LAMBDA,EPSHEL)
C
C   code:
C
C
C            Compute polarization vector:
C            Use Eqn A.12 and combine according to A.13
C            For outgoing bosons (sigmak=+1) the complex
C            conjugate polarisation vector is returned.
C
      EPSHEL(0) = 0.0d0
      IF (ABS(LAMBDA) .EQ. 1) THEN
          KT = SQRT(2D0 *(KBAR(1)**2+KBAR(2)**2))
            IF (KT .GT. 0) THEN
            IF (KBAR(4) .EQ. 0 ) THEN
                  KMAGNI = KBAR(0)
            ELSE
                  KMAGNI = SQRT(KBAR(1)**2+KBAR(2)**2+KBAR(3)**2)
            ENDIF
            NORMAL = -LAMBDA/(KMAGNI*KT)
            A1 = KBAR(1)*KBAR(3)*NORMAL
            A2 = KBAR(2)*KBAR(3)*NORMAL
            A3 =    -KT**2 * 0.5d0 *NORMAL
            B1 = -KBAR(2)/KT
            B2 =  KBAR(1)/KT
            EPSHEL(1) = dcmplx(A1, SIGMAK * B1)
            EPSHEL(2) = dcmplx(A2, SIGMAK * B2)
            EPSHEL(3) = A3
          ELSE
            NORMAL = -LAMBDA / sqrt( 2D0 )
            EPSHEL(1) = dsign(1.0d0, KBAR(3)) * NORMAL
            EPSHEL(2) = dcmplx(0.0d0, -LAMBDA * SIGMAK * NORMAL)
            EPSHEL(3) = 0.0d0
          ENDIF
      ELSE IF (LAMBDA .EQ. 0) THEN
            IF (KBAR(4) .EQ. 0) THEN
C
C                  Mass = 0 not allowed for longitudinal polarization
C
                  WRITE(*,*) "Mass = 0 for Lambda = 0 in HELVEC"
                  EPSHEL(1) = 0.0d0
                  EPSHEL(2) = 0.0d0
                  EPSHEL(3) = 0.0d0
                  RETURN
                  END IF
            KMAGNI = SQRT(KBAR(1)**2+KBAR(2)**2+KBAR(3)**2)
            NORMAL = KBAR(4)/(KBAR(0)*KMAGNI)
            EPSHEL(1) = KBAR(1)*NORMAL
            EPSHEL(2) = KBAR(2)*NORMAL
            EPSHEL(3) = KBAR(3)*NORMAL
      ELSE
C
C                  Unrecognised value of Lambda
C
            WRITE (*,*) "Invalid Lambda in HELVEC: Lambda = ",LAMBDA
                  EPSHEL(1) = 0.0d0
                  EPSHEL(2) = 0.0d0
                  EPSHEL(3) = 0.0d0
            END IF
      RETURN
      END
C
C
C******************************** bra2r, bra2c ************************
C
C   Calculates the bra vector <i,k| as given in
C      MAD/PH/402, Appendix A, Eqn A.10, top eqn.
C      checked on May 26, 88
C
C   Modified to double precision and nonzero time component of 
C   polarization vector Aug. 92
C
C   INPUT:
C            chi      double complex array(2)
C                        bra <i| as given by eq. A.9 if chreal = .true.
C                        any bra <...| if chreal = .false.
C
C            chreal  logical
C                        .true. if one component of chi is real
C                        .false. otherwise
C
C            p      double precision  array(0:3)
C                        standard four momentum for fermion
C                        p(0:3) = [E, px, py, pz]
C
C            sigmap      integer
C                        chirality/helicity factor for fermion
C                        allowed values: +1,-1
C
C            k      double precision  array(0:4)
C                        standard four momentum for boson and mass**2
C                        k(0:4) = [E, kx, ky, kz, mass**2]
C
C            reps     double precision  array(0:3)
C            ceps     double complex array(0:3)
C                     real or complex polarisation vector of the boson.
C
C   OUTPUT:
C            result      double complex array(1:2)
C                        two component row vector <i,k| on Eqn. A.10
C
C            pplusk      double precision  array(0:4)
C                        four momentum p+k. The fourth component
C                        contains the square (p+k)**2
C
C
      SUBROUTINE BRA2R(CHI,CHREAL,P,SIGMAP,K,REPS,RESULT,PPLUSK)
c
c arguments
c
      DOUBLE COMPLEX  CHI(2), CEPS(0:3), RESULT(2)
      DOUBLE PRECISION  P(0:3), K(0:4), REPS(0:3), PPLUSK(0:4)
      INTEGER  SIGMAP
      LOGICAL  CHREAL
c
c local variables
c
      DOUBLE PRECISION  A0,A1,A2,A3,B0,B1,B2,B3,C,CR,CI,
     &                  D0,D1,D2,D3,CCR,CCI
      DOUBLE PRECISION  T1R,T1I,T2R,T2I
      DOUBLE PRECISION  PROP
      INTEGER  INDEX1, INDEX2
      LOGICAL  DEXIST
c
c Compute product of the three matrices in A.10 explicitly in terms 
C of components.
c If one component of Chi is real it is denoted by "c", the other is 
C complex and separated as cr + i*ci. If both are complex replace 
C c --> ccr + i*cci. The 2 possible helicity indices are combined by 
C using some Pauli matrix algebra.
c
c Computations are reduced by computing (chi*eps)*(p+k).
c The (chi*eps) temporary is stored in the row vector:
c   ( t1r + I*t1i, t2r + I*t2i )
c
c Since the 4-vector eps can be real or complex, the chi*eps part is 
C done with two different entries into the subroutine: bra2r and bra2c
c
      A0 = REPS(0)
      IF ( SIGMAP.EQ.1 ) THEN
         A1 = REPS(1)
         INDEX1 = 1
         INDEX2 = 2
      ELSE
         A1 = -REPS(1)
         INDEX1 = 2
         INDEX2 = 1
      ENDIF
      A2 = REPS(2)
      A3 = REPS(3)
      DEXIST = .FALSE.
      GOTO 1
c
c entry for complex polarization vector
c
      ENTRY BRA2C( CHI,CHREAL,P,SIGMAP,K,CEPS,RESULT,PPLUSK )
      A0 = DREAL(CEPS(0))
      D0 = DIMAG(CEPS(0))
      IF ( SIGMAP.EQ.1 ) THEN
         A1 = DREAL(CEPS(1))
         D1 = DIMAG(CEPS(1))
         INDEX1 = 1
         INDEX2 = 2
      ELSEIF ( SIGMAP.EQ.-1 ) THEN
         A1 = -DREAL(CEPS(1))
         D1 = -DIMAG(CEPS(1))
         INDEX1 = 2
         INDEX2 = 1
      ELSE
c
c unrecognised value of simgap
c
        WRITE (*,*) "Invalid Sigmap in BRA2 : Sigmap = ",SIGMAP
        RESULT(1) = 0.d0
        RESULT(2) = 0.d0
        RETURN
      ENDIF
      A2 = DREAL(CEPS(2))
      D2 = DIMAG(CEPS(2))
      A3 = DREAL(CEPS(3))
      D3 = DIMAG(CEPS(3))
      DEXIST = .TRUE.

  1   CONTINUE
c
c compute the sum (p+k) and store in pplusk:
c   (p+k)(0:3) = pplusk(0:3), (p+k)**2 = pplusk(4)
c
      PPLUSK(0) = P(0) + K(0)
      PPLUSK(1) = P(1) + K(1)
      PPLUSK(2) = P(2) + K(2)
      PPLUSK(3) = P(3) + K(3)
      PPLUSK(4) = PPLUSK(0)**2 - PPLUSK(1)**2 - 
     &            PPLUSK(2)**2 - PPLUSK(3)**2
      PROP = 1.d0/PPLUSK(4)

      B0 = PPLUSK(0)
      IF ( SIGMAP.EQ.1 ) THEN
            B1 = PPLUSK(1)
      ELSE
            B1 = - PPLUSK(1)
      ENDIF
      B2 = PPLUSK(2)
      B3 = PPLUSK(3)
c
c now calculate Chi*eps*prop
c
      IF ( CHREAL ) THEN
         C  = DREAL(CHI(INDEX1))*PROP
         CR = DREAL(CHI(INDEX2))*PROP
         CI = DIMAG(CHI(INDEX2))*PROP

         IF ( DEXIST ) THEN
            T1R = -C*(A3-A0) + CR*(D2-A1) + CI*(A2+D1)
            T1I = -C*(D3-D0) + (D2-A1)*CI - (A2+D1)*CR
            T2R =  CR*(A3+A0) - C*(A1+D2) - CI*(D3+D0)
            T2I =  C*(A2-D1) + CI*(A3+A0) + CR*(D3+D0)
         ELSE
            T1R = -C*(A3-A0) - CR*A1 + CI*A2
            T1I = -A1*CI - A2*CR
            T2R = CR*(A3+A0) - C*A1
            T2I = C*A2 + CI*(A3+A0)
         ENDIF
      ELSE
         CCR = DREAL(CHI(INDEX1))*PROP
         CCI = DIMAG(CHI(INDEX1))*PROP
         CR  = DREAL(CHI(INDEX2))*PROP
         CI  = DIMAG(CHI(INDEX2))*PROP

         IF ( DEXIST ) THEN
            T1R = -CCR*(A3-A0) + CCI*(D3-D0) + CR*(D2-A1) + CI*(A2+D1)
            T1I = -CCR*(D3-D0) - CCI*(A3-A0) + (D2-A1)*CI - (A2+D1)*CR
            T2R =  CR*(A3+A0) - CCR*(A1+D2) - CCI*(A2-D1) - CI*(D3+D0)
            T2I = CCR*(A2-D1) - CCI*(A1+D2) + CI*(A3+A0) + CR*(D3+D0)
         ELSE
            T1R = -CCR*(A3-A0) - CR*A1 + CI*A2
            T1I = -CCI*(A3-A0) - A1*CI - A2*CR
            T2R = CR*(A3+A0) - CCR*A1 - CCI*A2
            T2I = CCR*A2 - CCI*A1 + CI*(A3+A0)
         ENDIF
      ENDIF

      RESULT(INDEX1) = DCMPLX( T1R*(B0+B3) + T2R*B1 - T2I*B2,
     &                         T1I*(B0+B3) + T2I*B1 + T2R*B2 )
      RESULT(INDEX2) = DCMPLX( T1R*B1 + T1I*B2 + T2R*(B0-B3),
     &                         T1I*B1 - T1R*B2 + T2I*(B0-B3) )
ccc
      END
C
C********************* ket2r, ket2c **********************************
C
C   Calculates the ket vector |k,i> as given in
C      MAD/PH/402, Appendix A, Eqn A.10, bottom eqn.
C      checked on May 26, 88
C
C   Modified to double precision and nonzero time component of 
C   polarization vector Aug. 92
C
C   INPUT:
C            chi      double complex array(2)
C                        ket |i> as given by eq. A.9 if chreal = .true.
C                        any ket |...> if chreal = .false.
C
C            chreal  logical
C                        .true. if one component of chi is real
C                        .false. otherwise
C
C            p      double precision  array(0:3)
C                        standard four momentum for fermion
C                        p(0:3) = [E, px, py, pz]
C
C            sigmap      integer
C                        chirality/helicity factor for fermion
C                        allowed values: +1,-1
C
C            k      double precision  array(0:4)
C                        standard four momentum for boson and mass**2
C                        k(0:4) = [E, kx, ky, kz, mass**2]
C
C            reps   double precision  array(0:3)
C            ceps   double complex array(0:3)
C                   real or complex polarisation vector of the boson. 
C
C   OUTPUT:
C            result      double complex array(1:2)
C                        two component column vector |k,i> on Eqn. A.10
C
C            pmink      double precision  array(0:4)
C                        four momentum p-k. The fourth component
C                        contains the square (p-k)**2
C
      SUBROUTINE KET2R( CHI,CHREAL,P,SIGMAP,K,REPS,RESULT,PMINK )
c
c arguments
c
      DOUBLE COMPLEX  CHI(2), CEPS(0:3), RESULT(2)
      DOUBLE PRECISION  P(0:3), K(0:4), REPS(0:3), PMINK(0:4)
      INTEGER  SIGMAP
      LOGICAL  CHREAL
c
c local variables
c
      DOUBLE PRECISION   A0,A1,A2,A3,B0,B1,B2,B3,C,CR,CI,
     &                   D0,D1,D2,D3,CCI,CCR
      DOUBLE PRECISION   T1R,T1I,T2R,T2I
      DOUBLE PRECISION   PROP
      INTEGER  INDEX1, INDEX2
      LOGICAL  DEXIST
c
c Compute product of the three matrices in A.10 explicitly in terms 
C of components.
c If one component of Chi is real it is denoted by "c", the other is 
C complex and separated as cr + i*ci. If both are complex replace 
C c --> ccr + i*cci. The 2 possible helicity indices are combined by 
C using some Pauli matrix algebra.
c
c Computations are reduced by computing (p-k)*(eps*chi).
c The (eps*chi) temporary is stored in the vector:
c     ( t1r + I*t1i, t2r + I*t2i )
c
c Since the 4-vector eps can be real or complex, the chi*eps part is 
C done with two different entries into the subroutine: bra2r and bra2c
c
      A0 = REPS(0)
      IF ( SIGMAP.EQ.1 ) THEN
         A1 = REPS(1)
         INDEX1 = 1
         INDEX2 = 2
      ELSE
         A1 = -REPS(1)
         INDEX1 = 2
         INDEX2 = 1
      ENDIF
      A2 = REPS(2)
      A3 = REPS(3)
      DEXIST = .FALSE.
      GOTO 1
c
c entry for complex polarization vector
c
      ENTRY KET2C( CHI,CHREAL,P,SIGMAP,K,CEPS,RESULT,PMINK )
      A0 = DREAL(CEPS(0))
      D0 = DIMAG(CEPS(0))
      IF ( SIGMAP.EQ.1 ) THEN
         A1 = DREAL(CEPS(1))
         D1 = DIMAG(CEPS(1))
         INDEX1 = 1
         INDEX2 = 2
      ELSEIF ( SIGMAP.EQ.-1 ) THEN
         A1 = -DREAL(CEPS(1))
         D1 = -DIMAG(CEPS(1))
         INDEX1 = 2
         INDEX2 = 1
      ELSE
c
c unrecognised value of simgap
c
         WRITE(6,*) "Invalid Sigmap in BRA2 : Sigmap = ",SIGMAP
         RESULT(1) = 0.d0
         RESULT(2) = 0.d0
         RETURN
      ENDIF
      A2 = DREAL(CEPS(2))
      D2 = DIMAG(CEPS(2))
      A3 = DREAL(CEPS(3))
      D3 = DIMAG(CEPS(3))
      DEXIST = .TRUE.

  1   CONTINUE
c
c compute the sum (p-k) and store in pmink:
c   (p-k)(0:3) = pmink(0:3), (p-k)**2 = pmink(4)
c
      PMINK(0) = P(0) - K(0)
      PMINK(1) = P(1) - K(1)
      PMINK(2) = P(2) - K(2)
      PMINK(3) = P(3) - K(3)
      
      PMINK(4) = PMINK(0)**2 - PMINK(1)**2 - PMINK(2)**2 - PMINK(3)**2
      PROP = 1.d0/PMINK(4)

      B0 = PMINK(0)
      IF ( SIGMAP.EQ.1 ) THEN
            B1 = PMINK(1)
      ELSE
            B1 = -PMINK(1)
      ENDIF
      B2 = PMINK(2)
      B3 = PMINK(3)
c
c now calculate Chi*eps*prop
c
        IF ( CHREAL ) THEN
           C  = DREAL(CHI(INDEX1))*PROP
           CR = DREAL(CHI(INDEX2))*PROP
           CI = DIMAG(CHI(INDEX2))*PROP

           IF ( DEXIST ) THEN
              T1R = -C*(A3-A0) - CR*(D2+A1) - CI*(A2-D1)
              T1I = -C*(D3-D0) - CI*(D2+A1) + CR*(A2-D1)
              T2R =  C*(D2-A1) + CR*(A3+A0) - CI*(D3+D0)
              T2I = -C*(A2+D1) + CI*(A3+A0) + CR*(D3+D0)
           ELSE
              T1R = -C*(A3-A0) - CR*A1 - CI*A2
              T1I = -CI*A1 + CR*A2
              T2R = -C*A1 + CR*(A3+A0)
              T2I = -C*A2 + CI*(A3+A0)
           ENDIF
        ELSE
           CCR = DREAL(CHI(INDEX1))*PROP
           CCI = DIMAG(CHI(INDEX1))*PROP
           CR  = DREAL(CHI(INDEX2))*PROP
           CI  = DIMAG(CHI(INDEX2))*PROP

           IF ( DEXIST ) THEN
              T1R = -CCR*(A3-A0)+ CCI*(D3-D0) - CR*(D2+A1) - CI*(A2-D1)
              T1I = -CCR*(D3-D0)- CCI*(A3-A0) - CI*(D2+A1) + CR*(A2-D1)
              T2R =  CCR*(D2-A1)+ CCI*(A2+D1) + CR*(A3+A0) - CI*(D3+D0)
              T2I =  CCI*(D2-A1)- CCR*(A2+D1) + CI*(A3+A0) + CR*(D3+D0)
           ELSE
              T1R = -CCR*(A3-A0) - CR*A1 - CI*A2
              T1I = -CCI*(A3-A0) - CI*A1 + CR*A2
              T2R = -CCR*A1 + CCI*A2 + CR*(A3+A0)
              T2I = -CCR*A2 - CCI*A1 + CI*(A3+A0)
           ENDIF
        ENDIF

        RESULT(INDEX1) = DCMPLX( T1R*(B0+B3) + T2R*B1 + T2I*B2,
     &                           T1I*(B0+B3) + T2I*B1 - T2R*B2 )
        RESULT(INDEX2) = DCMPLX( T1R*B1 - T1I*B2 + T2R*(B0-B3),
     &                           T1I*B1 + T1R*B2 + T2I*(B0-B3) )
ccc
      END
C
