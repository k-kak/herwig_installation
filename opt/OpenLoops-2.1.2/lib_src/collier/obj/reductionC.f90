!!
!!  File reductionC.F90 is part of COLLIER
!!  - A Complex One-Loop Library In Extended Regularizations
!!
!!  Copyright (C) 2015, 2016   Ansgar Denner, Stefan Dittmaier, Lars Hofer
!!
!!  COLLIER is licenced under the GNU GPL version 3, see COPYING for details.
!!

!#define Credtest
!#define Cpvtest
!#define Cpv1test
!#define Cpv1otest
!#define Cpv2test
!#define Cpvshifttest
!#define Cgtest 
!#define Cgytest
!#define Cgrtest
!#define Cgptest
!#define Cgpftest

!#define USEC0
!#define PPEXP00

!#define TRACECin
!#define TRACECout
!#define CritPointsCOLI



!#define TEST
!#define Cgntest

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!  ***********************
!  *  module reductionC  *
!  *    by Lars Hofer    *
!  * adapted by A Denner *
!  ***********************
! 
!  functions and subroutines:
!  CalcCuv, CalcCpv, CalcCpv2, CalcCg, CalcCgy, CalcCgp, CalcCgr, CalcCgpf, CopyCimp3
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

module globalC

  double complex :: q10,q21,q20,mm02,mm12,mm22
  double complex :: detZ,Z(2,2),Zadj(2,2),f(2),Zadjf(2),Xadj(0:2,0:2),detX
  double precision :: q2max,m2max,m2scale,adetZ,adetX,fmax,maxZadjf,maxZadjfd,maxXadj,aZadjff,maxZ
  double complex :: Zinv(2,2),Zadjs(2),detZmZadjf
  double complex :: mx(0:2,0:2), mxinv(0:2,0:2)
  double precision :: maxZadj
  double precision :: fac_g,fac_gy,fac_gp,fac_gr,fac_gpf
  double precision :: wmaxZadj,wmaxZadjf,wmaxXadj

  double complex :: q10shift,q21shift,q20shift,mm02shift,mm12shift,mm22shift
  double precision :: adetZshift,adetXshift,fmaxshift,maxZadjfshift,aZadjffshift,maxZshift
  double complex :: detZshift,Zshift(2,2),Zadjshift(2,2),fshift(2),Zadjfshift(2),Xadjshift(0:2,0:2),detXshift
  double complex :: Zinvshift(2,2),Zadjsshift(2),detZmZadjfshift
  double complex :: mxshift(0:2,0:2), mxinvshift(0:2,0:2)
  double precision ::  maxZadjshift
  
  double complex, parameter :: undefined_C=1d50
  
end module globalC





module reductionC

  use coli_stat
  use reductionAB


  implicit none


  !  should not be too small since expansion for large expansion parameters are calculated to early
  double precision, parameter :: truncfacC = 1d2
  ! double precision, parameter :: truncfacC = 1d3   worse than 1d2?
!  double precision, parameter :: acc_C=1d-13

contains


  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !  subroutine CalcC(C,Cuv,p10,p21,p20,m02,m12,m22,rmax,id,Cerr1,Cerr2,rbasic,acc_req_Cextra)
  !
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  
  subroutine CalcC(C,Cuv,p10,p21,p20,m02,m12,m22,rmax,id,Cerr1,Cerr2,rbasic,acc_req_Cextra)
  
    integer, intent(in) :: rmax, id
    integer, optional, intent(in) :: rbasic
    double precision, optional, intent(in) :: acc_req_Cextra(0:rmax)
    double complex, intent(in) :: p10,p21,p20,m02,m12,m22
    double complex, intent(out) :: Cuv(0:rmax,0:rmax,0:rmax)
    double complex, intent(out) :: C(0:rmax,0:rmax,0:rmax)
    double precision, intent(out) :: Cerr1(0:rmax),Cerr2(0:rmax)
    double complex, allocatable :: Caux(:,:,:), Cuvaux(:,:,:), fct(:)
    double precision, allocatable :: Cerr1aux(:),Cerr2aux(:),acc_req_Cexnew(:)
    double complex :: x(6)
    integer :: rank,switch,cnt,n0,n1,n2,r,rb
    logical :: nocalc,wrica



    if (use_cache_system) then
      if ((ncache.gt.0).and.(ncache.le.ncache_max)) then
!        if (use_cache(ncache).ge.3) then
          x(1)=p10
          x(2)=p21
          x(3)=p20
          x(4)=m02
          x(5)=m12
          x(6)=m22
          rank = rmax
          switch = 0

          if(present(rbasic)) then
            rb =rbasic
          else
            rb = rmax
          endif

          if(rmax.ge.1) then
            allocate(fct(NCoefsG(rmax,3)+NCoefsG(rmax-1,3)+2*(rmax+1)))
            call ReadCache(fct,NCoefsG(rmax,3)+NCoefsG(rmax-1,3)+2*(rmax+1),x,6,1,id,3,rb,nocalc,wrica)
          else 
            allocate(fct(NCoefsG(rmax,3)+2*(rmax+1)))
            call ReadCache(fct,NCoefsG(rmax,3)+2*(rmax+1),x,6,1,id,3,rb,nocalc,wrica)
          end if
    
          rank = max(rmax,rb)

          if(nocalc)then

              cnt = 0
              do r=0,rmax
                do n0=0,r
                  do n1=0,r-n0
                    n2 = r-n0-n1
 
                    cnt = cnt+1
                    C(n0,n1,n2) = fct(cnt)

                  end do
                end do
                do n0=1,r
                  do n1=0,r-n0
                    n2 = r-n0-n1
 
                    cnt = cnt+1
                    Cuv(n0,n1,n2) = fct(cnt)

                  end do
                end do
                cnt = cnt+1
                Cerr1(r) = real(fct(cnt))
                cnt = cnt+1
                Cerr2(r) = real(fct(cnt))
              end do

              return

          end if

          if(rank.eq.rmax) then

            if(present(rbasic)) then
              if(rb.gt.rbasic) then
                allocate(acc_req_Cexnew(0:rank))
                acc_req_Cexnew(0:rb)=acc_req_Cextra(0)
                acc_req_Cexnew(rb+1:rank)=acc_req_Cextra(rbasic+1:rank-rb+rbasic) 
                call CalcCred(C,Cuv,p10,p21,p20,m02,m12,m22,rank,id,Cerr1,Cerr2,rb,acc_req_Cexnew)
              else
                call CalcCred(C,Cuv,p10,p21,p20,m02,m12,m22,rank,id,Cerr1,Cerr2,rbasic,acc_req_Cextra)
              end if
            else
              call CalcCred(C,Cuv,p10,p21,p20,m02,m12,m22,rank,id,Cerr1,Cerr2)
            end if

            if (wrica) then
              cnt = 0
              do r=0,rank
                do n0=0,r
                  do n1=0,r-n0
                    n2 = r-n0-n1
 
                    cnt = cnt+1
                    fct(cnt) = C(n0,n1,n2)
                  end do
                end do
                do n0=1,r
                  do n1=0,r-n0
                    n2 = r-n0-n1
 
                    cnt = cnt+1
                    fct(cnt) = Cuv(n0,n1,n2)
                  end do
                end do
                cnt = cnt+1
                fct(cnt) = Cerr1(r)
                cnt = cnt+1
                fct(cnt) = Cerr2(r)
              end do

              if(rank.ge.1) then
                call WriteCache(fct,NCoefsG(rank,3)+NCoefsG(rank-1,3)+2*(rank+1),id,3,rb)
              else 
                call WriteCache(fct,NCoefsG(rank,3)+2*(rank+1),id,3,rb)
              end if

            end if

            return
          
          
          else

            allocate(Caux(0:rank,0:rank,0:rank))
            allocate(Cuvaux(0:rank,0:rank,0:rank))
            allocate(Cerr1aux(0:rank))
            allocate(Cerr2aux(0:rank))

            if(present(rbasic)) then
              if(rb.gt.rbasic) then
                allocate(acc_req_Cexnew(0:rank))
                acc_req_Cexnew(0:rb)=acc_req_Cextra(0)
                acc_req_Cexnew(rb+1:rank)=acc_req_Cextra(rbasic+1:rank-rb+rbasic) 
!                acc_req_Cexnew(rb+1:rmax)=acc_req_Cextra(rbasic+1:rmax-rb+rbasic) 
!                acc_req_Cexnew(rmax+1:rank)=acc_req_Cextra(rmax) 
                call CalcCred(C,Cuv,p10,p21,p20,m02,m12,m22,rank,id,Cerr1,Cerr2,rb,acc_req_Cexnew)
              else
                call CalcCred(C,Cuv,p10,p21,p20,m02,m12,m22,rank,id,Cerr1,Cerr2,rbasic,acc_req_Cextra)
              end if

              call CalcCred(Caux,Cuvaux,p10,p21,p20,m02,m12,m22,rank,id,Cerr1aux,Cerr2aux,rbasic+rank-rmax,acc_req_Cextra)
            else
              call CalcCred(Caux,Cuvaux,p10,p21,p20,m02,m12,m22,rank,id,Cerr1aux,Cerr2aux)
            end if

            if (wrica) then
              cnt = 0
              deallocate(fct)
              if(rank.ge.1) then
                allocate(fct(NCoefsG(rank,3)+NCoefsG(rank-1,3)+2*(rank+1)))
              else 
                allocate(fct(NCoefsG(rank,3)+2*(rank+1)))
              end if
              do r=0,rank
                do n0=0,r
                  do n1=0,r-n0
                    n2 = r-n0-n1
 
                    cnt = cnt+1
                    fct(cnt) = Caux(n0,n1,n2)
                  end do
                end do
                do n0=1,r
                  do n1=0,r-n0
                    n2 = r-n0-n1
 
                    cnt = cnt+1
                    fct(cnt) = Cuvaux(n0,n1,n2)
                  end do
                end do
                cnt = cnt+1
                fct(cnt) = Cerr1aux(r)
                cnt = cnt+1
                fct(cnt) = Cerr2aux(r)
              end do

              if(rank.ge.1) then
                call WriteCache(fct,NCoefsG(rank,3)+NCoefsG(rank-1,3)+2*(rank+1),id,3,rb)
              else 
                call WriteCache(fct,NCoefsG(rank,3)+2*(rank+1),id,3,rb)
              end if

            end if

            C = Caux(0:rmax,0:rmax,0:rmax)
            Cuv = Cuvaux(0:rmax,0:rmax,0:rmax)
            Cerr1 = Cerr1aux(0:rmax)
            Cerr2 = Cerr2aux(0:rmax)

            deallocate(Caux)
            deallocate(Cuvaux)
            deallocate(Cerr1aux)
            deallocate(Cerr2aux)

            return

          end if
!        end if
      end if
    end if

    if(present(rbasic))then
      call CalcCred(C,Cuv,p10,p21,p20,m02,m12,m22,rmax,id,Cerr1,Cerr2,rbasic,acc_req_Cextra)
    else
      call CalcCred(C,Cuv,p10,p21,p20,m02,m12,m22,rmax,id,Cerr1,Cerr2)
    end if

!         write(*,*) 'Cred nc Cerr1',Cerr1
!         write(*,*) 'Cred nc Cerr2',Cerr2

  end subroutine CalcC





  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !  subroutine CalcCred(C,Cuv,p10,p21,p20,m02,m12,m22,rmax,id,Cerr1,Cerr2,rbasic,acc_req_Cextra)
  !
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  
  subroutine CalcCred(C,Cuv,p10,p21,p20,m02,m12,m22,rmax,id,Cerr1,Cerr2,rbasic,acc_req_Cextra)
  
    use globalC

    integer, intent(in) :: rmax,id
    integer, intent(in), optional :: rbasic
    !   rbasic defines rank of tensors that are needed by mastercall
    !          higher ranks up to rmax are needed for internal iterations
    double precision, intent(in), optional :: acc_req_Cextra(0:rmax)
    double complex, intent(in) :: p10,p21,p20,m02,m12,m22
    double complex, intent(out) :: C(0:rmax,0:rmax,0:rmax)
    double complex, intent(out) :: Cuv(0:rmax,0:rmax,0:rmax)
    double precision, intent(out)  ::Cerr1(0:rmax),Cerr2(0:rmax)
    double complex :: C_alt(0:rmax,0:rmax,0:rmax)
    double complex :: Cuv_alt(0:rmax,0:rmax,0:rmax)
    double precision ::  Cerr(0:rmax),Cerr_alt(0:rmax),Cerr1_alt(0:rmax),Cerr2_alt(0:rmax)
    double precision :: C0est,Ctyp
    double complex :: C0_coli
!    double complex :: detX,chdet
    double complex :: chdet

    double complex :: elimminf2_coli
    integer :: r,rid,n0,n1,n2,g,gy,gp,gr,gpf,i,rdef,iexp
    logical :: use_C0
    logical :: use_pv,use_pv2,use_g,use_gy,use_gp,use_gr,use_gpf,use_pvs

    integer :: r_alt,Crmethod(0:rmax),Crmethod_alt(0:rmax),CrCalc(0:rmax),CCalc
    double precision :: acc_pv_alt, acc_pv2_alt, acc_Cr_alt

    !   CalcC stores methods that have been calculated
    !   Crmethod(r) stores best method=used for rank r 

    double precision :: err_pv(0:rmax),err_pv2(0:rmax),err_g(0:rmax),err_gy(0:rmax),  &
                        err_gp(0:rmax),err_gr(0:rmax),err_gpf(0:rmax)
    double precision :: err_pvs(0:rmax)
    double precision :: h_pv,w_pv,v_pv,z_pv,h_pv2,w_pv2,v_pv2,z_pv2,hw_pv2
    double precision :: h_pvs,w_pvs,v_pvs,z_pvs
    double precision :: x_g,u_g,z_g,err_g_B(0:rmax),err_g_exp
    double precision :: x_gy,y_gy,v_gy,v1_gy,b_gy,err_gy_B(0:rmax),err_gy_exp
    double precision :: w_gp,v_gp,z_gp,err_gp_B(0:rmax),err_gp_exp
    double precision :: x_gr,y_gr,y1_gr,a_gr,err_gr_B(0:rmax),err_gr_exp
    double precision :: x_gpf,y_gpf,v_gpf,v1_gpf,b_gpf,err_gpf_B(0:rmax),err_gpf_exp
    double precision :: err_B,err_C0,err_C(0:rmax),err_inf,err_req_Cr(0:rmax),acc_req_Cr(0:rmax),acc_C(0:rmax)
    double precision :: checkest,norm,Cscale
    double precision :: deterr
    logical :: lerr_C0,errorwriteflag

    character(len=*),parameter :: fmt1 = "(A7,'dcmplx(',d25.18,' , ',d25.18,' )')"
    character(len=*),parameter :: fmt10 = "(A17,'(',d25.18,' , ',d25.18,' )')"

    integer, parameter :: MaxCritPointC=0

    integer, save :: CritPointCntC



    data CritPointCntC /0/



    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! for rank = 0 calculate C0 directly 
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    use_C0 = .false.  
    if (rmax.eq.0d0) then
    ! rough estimate for C0 to set the scale, to be improved
      Cscale = max(abs(p10),abs(p21),abs(p20),abs(m02),abs(m12),abs(m22))
      if(Cscale.ne.0d0) then 
        Cerr1 = acc_inf/Cscale
      else
        C0est = acc_inf
      end if
      Cerr2 = Cerr1

      Z(1,1) = 2d0*p10
      Z(2,1) = p10+p20-p21
      Z(1,2) = Z(2,1)
      Z(2,2) = 2d0*p20
      
      adetZ = abs(chdet(2,Z))






! changed 16.08.2018
!      if (adetZ.gt.dprec_cll*maxval(abs(Z))**2) then
      if (adetZ.gt.dprec_cll*maxval(abs(Z(1,:)))*maxval(abs(Z(2,:)))) then
        C(0,0,0) = C0_coli(p10,p21,p20,m02,m12,m22)
        Cuv(0,0,0) = 0d0
        if (C(0,0,0).ne.undefined_C) then         
          Cerr1(0) = acc_def_C0*max( abs(C(0,0,0)), 1d0/sqrt(adetZ) )
        else
          Cerr1(0) = acc_inf*abs(C(0,0,0))
        end if
        Cerr2(0) = Cerr1(0)
        if (Cerr1(0).lt.acc_req_C*abs(C(0,0,0))) then 
          return
        else
          use_C0 = .true.
        endif
      end if 
    end if

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! choose reduction scheme
    ! by estimating expected errors
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    if (present(rbasic)) then
      rdef = rbasic
    else
      rdef = rmax
    end if



    if (present(acc_req_Cextra)) then 
      acc_req_Cr = acc_req_Cextra
    else
      acc_req_Cr = acc_req_C
    end if



    ! eliminate infinitesimal masses
    mm02 = elimminf2_coli(m02)
    mm12 = elimminf2_coli(m12)
    mm22 = elimminf2_coli(m22)
    q10  = elimminf2_coli(p10)
    q21  = elimminf2_coli(p21)
    q20  = elimminf2_coli(p20)

    ! set mass scales
    q2max = max(abs(q10),abs(q21),abs(q20))
    m2max = max(abs(mm02),abs(mm12),abs(mm22))
    m2scale = max(q2max,m2max)
 
    ! Gram and related stuff
    Z(1,1) = 2d0*q10
    Z(2,1) = q10+q20-q21   !  = q1q2
    Z(1,2) = Z(2,1)
    Z(2,2) = 2d0*q20

    maxZ = maxval(abs(Z))

! changed 21.06.2018
! deterr added 17.01.2019
    call chinve(2,Z,Zinv,detZ,deterr)
!    detZ = chdet(2,Z)

!    write(*,*) 'reductionC detZ = ',detZ,dprec_cll/deterr

! added 17.01.2019
    if (deterr.lt.dprec_cll) detZ = 0d0



! changed 16.08.2018 more stable for PV, but implies less expansions
! and larger Derr1=Derrout in Ds  (local estimates based on err2!)
!    if (detZ.ne.0d0) then
!!      call chinv(2,Z,Zinv)
!      Zadj = Zinv * detZ
!    else
      Zadj(1,1) = Z(2,2)
      Zadj(2,1) = -Z(2,1)
      Zadj(1,2) = -Z(2,1)
      Zadj(2,2) = Z(1,1)
!    end if

    Zadjs(1) = q21 + q20 - q10 
    Zadjs(2) = q21 + q10 - q20 



    detZmZadjf = q21*Z(2,1)

!   write(*,*) 'Zn    ',Z   
!   write(*,*) 'Zinvn ',Zinv
!   write(*,*) 'Zadjn ',Zadj
!   write(*,*) 'detZn ',detZ

    adetZ = abs(detZ)
    maxZadj = max(abs(Zadj(1,1)),abs(Zadj(2,1)),abs(Zadj(2,2)))

    f(1) = q10+mm02-mm12
    f(2) = q20+mm02-mm22
    fmax = max(abs(f(1)),abs(f(2)))

    mx(0,0) = 2d0*mm02
! 25.08.17
!    mx(1,0) = q10 - mm12 + mm02
!    mx(2,0) = q20 - mm22 + mm02
!    mx(2,1) = q10+q20-q21
!    mx(1,1) = 2d0*q10
!    mx(2,2) = 2d0*q20
    mx(1,0) = f(1)
    mx(2,0) = f(2)
    mx(0,1) = mx(1,0)
    mx(0,2) = mx(2,0)
    mx(1:2,1:2) = Z(1:2,1:2)

! changed 21.06.2018
! deterr added 17.01.2019
    call chinve(3,mx,mxinv,detX,deterr)
!    detX = chdet(3,mx)

!    write(*,*) 'reductionC detX = ',detX,dprec_cll/deterr

! added 17.01.2019
    if (deterr.lt.dprec_cll) detX = 0d0



    if (detX.ne.0d0.and.maxZ.ne.0d0) then

!      write(*,*) 'CalcCred mx=',mx

!      call chinv(3,mx,mxinv)

!     write(*,*) 'CalcCred mxinv=',mxinv

      Xadj = mxinv * detX

!     write(*,*) 'CalcCred Xadj=',Xadj

      Zadjf(1:2) = -Xadj(0,1:2)

    else
! 25.08.17
!      mx(0,0) = 2d0*mm02
!      mx(1,0) = q10 - mm12 + mm02
!      mx(2,0) = q20 - mm22 + mm02
!      mx(0,1) = mx(1,0)
!      mx(1,1) = 2d0*q10
!      mx(2,1) = q10+q20-q21
!      mx(0,2) = mx(2,0)
!      mx(1,2) = mx(2,1)
!      mx(2,2) = 2d0*q20

      Zadjf(1) = Zadj(1,1)*f(1)+Zadj(2,1)*f(2)
      Zadjf(2) = Zadj(1,2)*f(1)+Zadj(2,2)*f(2)

      Xadj(2,2) = 2d0*mm02*Z(1,1) - f(1)*f(1)
      Xadj(1,1) = 2d0*mm02*Z(2,2) - f(2)*f(2)
      Xadj(2,1) = 2d0*mm02*Z(2,1) - f(1)*f(2)
      Xadj(1,2) = Xadj(2,1)
    end if



    maxZadjf = max(abs(Zadjf(1)),abs(Zadjf(2)))
    maxZadjfd = max(maxZadjf,adetZ)

    aZadjff = abs(Zadjf(1)*f(1) + Zadjf(2)*f(2))
! changed 16.08.2018
!    adetX = abs(2d0*mm02*detZ - Zadjf(1)*f(1) - Zadjf(2)*f(2))
    adetX = abs(detX)
    maxXadj = max(abs(Xadj(1,1)),abs(Xadj(2,1)),abs(Xadj(2,2)))




    ! quantities for modified error estimates
    ! momentum weights
!    do i = 1,2
!      pweight(i) = max(abs(Z(i,1))/maxval(abs(Z(1:2,1))),  &
!          abs(Z(i,2))/maxval(abs(Z(1:2,2))))
!    end do

!    wmaxZadj = max(pweight(1)*abs(Zadj(1,1)),pweight(1)*abs(Zadj(1,2)),  &
!        pweight(2)*abs(Zadj(2,1)),pweight(2)*abs(Zadj(2,2)))
!
!    wmaxZadjf = max(pweight(1)*abs(Zadjf(1)),pweight(2)*abs(Zadjf(2)))
!
!    wmaxXadj = max(pweight(1)*abs(Xadj(1,1)),   &
!        pweight(1)*abs(Xadj(1,2)),pweight(2)*abs(Xadj(2,1)),  &
!        pweight(2)*abs(Xadj(2,2)))
!    wmaxXadj = max(2d0*abs(mm02)*sqrt(adetZ*maxZadj/maxZ),maxZadj2ff*maxZadjf/(maxZadj*fmax))        

!    write(*,*) 'CalcCred pweight',pweight(1:2)
!    write(*,*) 'CalcCred wmaxZadj',maxZadj,wmaxZadj
!    write(*,*) 'CalcCred wmaxZadjf',maxZadjf,wmaxZadjf
!    write(*,*) 'CalcCred wmaxZadjf',maxXadj,wmaxXadj


    ! rough estimate for C0 to set the scale, to be improved
    Cscale = max(abs(p10),abs(p21),abs(p20),abs(m02),abs(m12),abs(m22))

! changed 09.09.16
     if(Cscale.ne.0d0) then 
       C0est = 1d0/Cscale
     else
       C0est = 1d0
     end if
!    if (adetZ.ne.0d0) then
!      C0est = 1d0/sqrt(adetZ)
!    elseif (m2max.ne.0d0) then
!      C0est = 1d0/m2max
!    else if (maxZ.ne.0d0) then
!      C0est = 1d0/maxZ
!    else
!      C0est = 1d0 
!    end if
    lerr_C0 = .false.




    err_inf = acc_inf*C0est

    err_req_Cr = acc_req_Cr * C0est

    CCalc = 0
    CrCalc = 0
    Crmethod = 0
    Cerr = err_inf
!    Cerr1 = err_inf       !  shifted above C0
!    Cerr2 = err_inf
    acc_C = acc_inf
    CCount(0) = CCount(0)+1

    ! error estimate for C0
    if (adetZ.ne.0d0) then
!      err_C0 = acc_def_C0*q2max/sqrt(adetZ) * C0est
      err_C0 = acc_def_C0*max( C0est, 1d0/sqrt(adetZ) )
    else
      err_C0 = acc_def_C0 * C0est
    end if
    err_B = acc_def_B
     

    ! estimate accuracy of PV-reduction
    h_pv = real(undefined_C)
    w_pv = real(undefined_C)
    v_pv = real(undefined_C)
    z_pv = real(undefined_C)
!   14.07.2017
!   16.08.2018 changed back, since detZ=0 set above if too small
    if (adetZ.lt.dprec_cll*maxZadjf.or.adetZ.eq.0d0) then
!    if (adetZ.lt.dprec_cll*maxZadjf.or.adetZ.lt.dprec_cll*maxZ**2.or.adetZ.eq.0d0) then
      use_pv = .false.
      err_pv = err_inf
    else
      use_pv = .true.
      err_pv(0) = err_C0
      if (rdef.gt.0) then

        h_pv = sqrt(adetZ)/maxZadj
        w_pv = max((maxZadjf*h_pv/adetZ)**2, abs(mm02)*maxZ*h_pv/adetZ, aZadjff*maxZ*(h_pv/adetZ)**2)
        v_pv = maxZadjf*h_pv/adetZ
        z_pv = maxZ*h_pv/adetZ




        if (mod(rdef,2).eq.1) then
          err_pv(rdef) = max( w_pv**((rdef-1)/2) * v_pv * err_C0,  &
              max(w_pv**((rdef-1)/2),1d0) * z_pv * err_B )



        else
          err_pv(rdef) = max( w_pv**(rdef/2) * err_C0,  &
              max(w_pv**(rdef/2-1) * v_pv, 1d0) * z_pv * err_B )
      


        end if
      end if
    end if

    ! estimate accuracy of alternative PV-reduction
    z_pv2 = real(undefined_C)
    v_pv2 = real(undefined_C)
    w_pv2 = real(undefined_C)
    hw_pv2 = real(undefined_C)
!   14.07.2017   
!   16.08.2018 changed back, since detZ=0 set above if too small
    if ((adetZ.lt.dprec_cll*maxZadjf).or.(adetX.lt.dprec_cll*maxval(abs(mx))*adetZ).or.adetZ.eq.0d0.or.adetX.eq.0d0) then
!    if ((adetZ.lt.dprec_cll*maxZadjf).or.(adetX.lt.dprec_cll*maxval(abs(mx))*adetZ).or.  &
!         (adetZ.lt.dprec_cll*maxZ**2).or.(adetX.lt.dprec_cll*fmax**2*maxZ).or.adetZ.eq.0d0.or.adetX.eq.0d0) then
      use_pv2 = .false.
      err_pv2 = err_inf
    else
      use_pv2 = .true.
      err_pv2(0) = err_C0
      if (rdef.gt.0) then
        w_pv2 = maxZadjf/adetZ

        h_pv2 = sqrt(adetZ)/maxZadj
        hw_pv2 =  w_pv2*h_pv2

        v_pv2 = maxXadj/adetZ
        z_pv2 = adetZ/adetX 

!        write(*,*) 'CalcCred: w_pv2',w_pv2,v_pv2,z_pv2,err_C0,err_B

        if (mod(rdef,2).eq.1) then
! change 21.10.15 for ! default
!          err_pv2(rdef) = max( err_C0 * max(w_pv2**rdef,w_pv2*v_pv2**((rdef-1)/2) ),  &
!              err_B * z_pv2 * max(w_pv2**(rdef+1),w_pv2, &
!                              w_pv2*v_pv2**((rdef-1)/2),w_pv2**2, & 
!                              v_pv2**((rdef+1)/2),v_pv2 ) )

          err_pv2(rdef) = max( err_C0 * max(hw_pv2**rdef,hw_pv2*v_pv2**((rdef-1)/2) ),  &
              err_B * z_pv2 * max(w_pv2*hw_pv2**(rdef),hw_pv2, &
                              w_pv2*hw_pv2*v_pv2**((rdef-1)/2), &
                              hw_pv2*v_pv2**((rdef-1)/2),w_pv2*hw_pv2, & 
                              v_pv2**((rdef+1)/2),v_pv2 ) )

!        write(*,*) 'CalcCred: err_pv2',rdef,err_C0 * max(1d0,w_pv2**rdef,v_pv2**((rdef-1)/2),w_pv2*v_pv2**((rdef-1)/2) ), &
!            err_B * max(1d0,z_pv2*w_pv2**(rdef+1),z_pv2*w_pv2, &
!              z_pv2*w_pv2*v_pv2**((rdef-1)/2),z_pv2*w_pv2**2, & 
!              z_pv2*v_pv2**((rdef+1)/2),z_pv2*v_pv2 ) 

        else
! change 21.10.15 for ! default
!          err_pv2(rdef) = max( err_C0 * max(w_pv2**rdef,v_pv2**(rdef/2)),  &
!              err_B * z_pv2 * max(w_pv2**(rdef+1),w_pv2,  &
!                                  w_pv2*v_pv2**(rdef/2), w_pv2**2,  &
!                                  v_pv2**(rdef/2),v_pv2) )
          err_pv2(rdef) = max( err_C0 * max(hw_pv2**rdef,v_pv2**(rdef/2)),  &
              err_B * z_pv2 * max(w_pv2*hw_pv2**(rdef),hw_pv2,  &
                                  w_pv2*v_pv2**(rdef/2), w_pv2*hw_pv2,  &
                                  v_pv2**(rdef/2),v_pv2) )
        end if
      end if
    end if 

    ! scale estimates down to allow trying other methods
    err_pv(rdef)  = err_pv(rdef)/impest_C
    err_pv2(rdef) = err_pv2(rdef)/impest_C     





! changed 16.11.16
!   Ctyp = real(undefined_C)
    Ctyp = C0est


    if(use_pv.or.use_pv2) then

      if (err_pv(rdef).le.err_pv2(rdef)) then



        ! use PV-reduction if appropriate
        call CalcCpv1(C,Cuv,p10,p21,p20,m02,m12,m22,rmax,id,Cerr1,Cerr2)

        Cerr = Cerr2

        CCount(1) = CCount(1)+1
        CrCalc(0:rmax)=CrCalc(0:rmax)+1
        CCalc=CCalc+1
        Crmethod(0:rmax)=1


        err_pv=Cerr

      else



        ! use alternative PV-reduction if appropriate
        call CalcCpv2(C,Cuv,p10,p21,p20,m02,m12,m22,rmax,id,Cerr1,Cerr2)

        Cerr = Cerr2

        CCount(2) = CCount(2)+1
        CrCalc(0:rmax)=CrCalc(0:rmax)+2
        CCalc=CCalc+2
        Crmethod(0:rmax)=2



        err_pv2=Cerr

      end if


      ! refine error estimate for C0
!      C0est = abs(C(0,0,0))
      err_C0 = acc_def_C0*max( abs(C(0,0,0)), 1d0/sqrt(adetZ) )
!      err_req_Cr = acc_req_Cr * abs(C(0,0,0))
      lerr_C0 = .true.


      if (rmax.ge.1) then
        Ctyp =  max(abs(C(0,0,0)),abs(C(0,1,0)),abs(C(0,0,1)))
      else
        Ctyp =  abs(C(0,0,0))
      end if
      if(Ctyp.eq.0d0) Ctyp = C0est
      err_req_Cr = acc_req_Cr * Ctyp




! check with Cerr = Cerr2: exp not tried => larger errors
      if (maxval(Cerr1-err_req_Cr).lt.0) then
        CCount(CCalc+CCountoffset0) = CCount(CCalc+CCountoffset0)+1
        return
      end if      
    else if (.not.use_C0) then  !  added 14.07.2017 adapted 12.4.2018
      C = 0d0
      Cuv = 0d0
      Cerr1 = err_inf
      Cerr2 = err_inf
    end if


        
    ! choose most promising expansion scheme
    ! Gram expansion
!    if (maxZadjf.ne.0d0) then
!    if (maxZadjf.gt.m2scale**2*dprec_cll) then   !  10.07.2017
    if (maxZadjf.gt.m2scale**2*1d1*dprec_cll) then   !  12.04.2018
      x_g = adetZ/maxZadjf
!      u_g = max(1d0,m2scale*m2scale/maxZadjf/6d0,abs(mm02)*q2max/maxZadjf/6d0)
!      u_g = max(1d0,fmax*fmax/maxZadjf/6d0,abs(mm02)*maxZ/maxZadjf/6d0)
! 03.03.15  large P counts!
!      u_g = max(1d0,fmax*fmax/maxZadjf/2d0,abs(mm02)*maxZ/maxZadjf/2d0)
! 24.04.15  term appear only combined
      u_g = max(1d0,maxXadj/maxZadjf/2d0)
      fac_g = x_g*u_g
      err_g = err_inf
      g = -1
      if (fac_g.ge.1) then
        use_g = .false.
        err_g_exp = err_inf
        z_g = real(undefined_C)
      else
        use_g = .true.
!       z_g = max(1d0,m2scale*q2max/maxZadjf)
        z_g = maxZ/maxZadjf
        err_g_B(rdef) = err_B * u_g**rdef * z_g 
        err_g_exp = u_g**(rdef-1) * Ctyp
      end if
    else
      use_g = .false.
      err_g = err_inf
      g = -1
      err_g_exp = err_inf
      u_g = real(undefined_C)
      z_g = real(undefined_C)
    end if



    ! Gram-Cayley expansion
!    if (maxXadj.ne.0d0.and.maxZ.ne.0d0) then
    if (maxXadj.gt.m2scale**2*dprec_cll.and.maxZ.gt.m2scale*dprec_cll) then   !  10.07.2017
      x_gy = maxZadjf/maxXadj
      y_gy = adetZ/maxXadj
!     v_gy = m2scale/q2max
      v_gy = fmax/maxZ
      v1_gy = max(1d0,v_gy)
      fac_gy = max(x_gy,y_gy)*v1_gy
      err_gy = err_inf
      gy = -1
      if (fac_gy.ge.1) then
        use_gy = .false.
        err_gy_exp = err_inf
        b_gy = real(undefined_C)
      else
        use_gy = .true.
!       b_gy = max(1d0,m2scale*q2max/maxXadj)
        b_gy = maxZ/maxXadj
        err_gy_B(rdef) = err_B * b_gy*v1_gy
        err_gy_exp = 1d0 * Ctyp
      end if
    else
      use_gy = .false.
      err_gy = err_inf
      gy = -1
      err_gy_exp = err_inf
      v1_gy = real(undefined_C)
      b_gy = real(undefined_C)
    end if 



    ! expansion in small momenta
!    if (fmax.ne.0d0) then
    if (fmax.gt.m2scale*dprec_cll) then   !  10.07.2017
!     w_gp = q2max/fmax
      w_gp = maxZ/fmax
      v_gp = max(1d0,abs(mm02)/fmax)
      fac_gp = w_gp*v_gp
      err_gp = err_inf
      gp = -1
      if (fac_gp.ge.1d0) then
        use_gp = .false.
        err_gp_exp = err_inf
        z_gp = real(undefined_C)
      else
        use_gp = .true.
!       z_gp = max(1d0,m2scale/fmax)
        z_gp = 1d0/fmax
        err_gp_B(rdef) = err_B * z_gp*v_gp**rdef
        err_gp_exp = v_gp**(rdef-1) * Ctyp
      end if
    else
      use_gp = .false.
      err_gp = err_inf
      gp = -1
      err_gp_exp = err_inf
      z_gp = real(undefined_C)
      v_gp = real(undefined_C)
    end if



    ! reversed Gram expansion
!    if (maxZadjf.ne.0d0.and.fmax.ne.0d0) then
    if (maxZadjf.gt.m2scale**2*dprec_cll.and.fmax.gt.m2scale*dprec_cll) then   !  10.07.2017
      x_gr = adetZ/maxZadjf
      y_gr = maxZadj/fmax              ! c*y    c=2
      y1_gr = max(1d0,y_gr)
      a_gr = maxZadj/maxZadjf
      fac_gr = max(x_gr,y_gr)
      err_gr = err_inf
      gr = -1
      if (fac_gr.ge.1.or.2*rmax.gt.rmax_B) then
        use_gr = .false.
        err_gr_exp = err_inf
      else
        use_gr = .true.
        err_gr_B(rdef) = err_B * a_gr
        err_gr_exp = y1_gr * Ctyp
      end if
    else
      use_gr = .false.
      err_gr = err_inf
      gr = -1
      err_gr_exp = err_inf
      y1_gr = real(undefined_C)
      a_gr = real(undefined_C)
    end if



    !  expansion in small momenta and f's
!  estimates to be confirmed 16.08.2017, r dependence may be different
!  since C_mni... is needed in contrast to Cgy expansion
    if (abs(m02).gt.m2scale*dprec_cll) then 
      x_gpf = fmax/abs(m02)
      y_gpf =  maxZ/abs(m02)
      v_gpf = 0d0
      v1_gpf = max(1d0,v_gpf)
      fac_gpf = max(x_gpf,y_gpf)*v1_gpf
      err_gpf = err_inf
      gpf = -1
      if (fac_gpf.ge.1) then
        use_gpf = .false.
        err_gpf_exp = err_inf
        b_gpf = real(undefined_C)
      else
        use_gpf = .true.
        b_gpf = 1d0/abs(m02)
        err_gpf_B(rdef) = err_B * b_gpf*v1_gpf
        err_gpf_exp = 1d0 * Ctyp
      end if
    else
      use_gpf = .false.
      err_gpf = err_inf
      gpf = -1
      err_gpf_exp = err_inf
      v1_gpf = real(undefined_C)
      b_gpf = real(undefined_C)
    end if 




! no method works
    if(use_C0.or.use_pv.or.use_pv2.or.use_g.or.use_gy.or.use_gp.or.use_gr.or.use_gpf.eqv..false.) then

! added 16.08.2018
      goto 200   !  try pv with shift

      call SetErrFlag_coli(-6)
      call ErrOut_coli('CalcCred',' no reduction method works', &
           errorwriteflag)
!      write(nerrout_coli,'((a))')  '  no reduction method works'
      if (errorwriteflag) then
        write(nerrout_coli,fmt10) ' CalcCred: p10 = ',p10
        write(nerrout_coli,fmt10) ' CalcCred: p21 = ',p21
        write(nerrout_coli,fmt10) ' CalcCred: p20 = ',p20
        write(nerrout_coli,fmt10) ' CalcCred: m02 = ',m02
        write(nerrout_coli,fmt10) ' CalcCred: m12 = ',m12
        write(nerrout_coli,fmt10) ' CalcCred: m22 = ',m22   
      end if
      C = 0d0
      Cuv = 0d0
      Cerr = err_inf
      Cerr2 = err_inf



      return
    endif



    iexp = 0
    do i=0,rmax_C-rmax
 
      if (use_g) then
        if (err_g_exp.gt.err_g_B(rdef)) then
          g = i
          err_g_exp = err_g_exp*fac_g
          err_g(rdef) = max(err_g_exp,err_g_B(rdef))
          if(err_g(rdef).lt.err_req_Cr(rdef)) then
            iexp = 1
            ! increase g by 2 to account for bad estimates
            g = min(max(g+2,2*g),rmax_C-rmax)
            exit
          end if



        end if
      end if

      if (mod(i,2).eq.1) then
        if (use_gy) then
          if (err_gy_exp.gt.err_gy_B(rdef)) then
            gy = i/2
            err_gy_exp = err_gy_exp*fac_gy
            err_gy(rdef) = max(err_gy_exp, err_gy_B(rdef))
            if(err_gy(rdef).lt.err_req_Cr(rdef)) then
              iexp = 2
              ! increase gy by 2 to account for bad estimates
              gy = min(max(gy+4,2*gy),(rmax_C-rmax)/2)
              exit
            end if



          end if
        end if
      end if

      if (use_gp) then
        if (err_gp_exp.gt.err_gp_B(rdef)) then
          gp = i
          err_gp_exp = err_gp_exp*fac_gp
          err_gp(rdef) = max(err_gp_exp,err_gp_B(rdef))
          if(err_gp(rdef).lt.err_req_Cr(rdef)) then
            iexp = 3
            ! increase gp by 2 to account for bad estimates
            gp = min(max(gp+2,2*gp),rmax_C-rmax)
            exit
          end if



        end if
      end if

      if (mod(i,2).eq.1) then

        if (use_gr) then


 
          if (err_gr_exp.gt.err_gr_B(rdef)) then
            gr = i/2
            err_gr_exp = err_gr_exp*fac_gr
            err_gr(rdef) = max(err_gr_exp, err_gr_B(rdef))
            if(err_gr(rdef).lt.err_req_Cr(rdef)) then
              iexp = 4
              ! increase gy by 2 to account for bad estimates
! changed 28.07.14
!             gr = min(max(gr+4,2*gr),(rmax_C-rmax)/2)
              gr = min(max(gr+4,2*gr),rmax_C-rmax,max(0,(rmax_B-2*rmax)/2))
              exit
            end if
          end if



        end if
      end if

      if (mod(i,2).eq.1) then
        if (use_gpf) then
          if (err_gpf_exp.gt.err_gpf_B(rdef)) then
            gpf = i/2
            err_gpf_exp = err_gpf_exp*fac_gpf
            err_gpf(rdef) = max(err_gpf_exp, err_gpf_B(rdef))
            if(err_gpf(rdef).lt.err_req_Cr(rdef)) then
              iexp = 5
              ! increase gpf by 2 to account for bad estimates
              gpf = min(max(gpf+4,2*gpf),(rmax_C-rmax)/2)
              exit
            end if



          end if
        end if
      end if

    end do

    ! scale estimates down to allow trying other methods
    err_g(rdef)  =  err_g(rdef)/impest_C
    err_gy(rdef) =  err_gy(rdef)/impest_C     
    err_gp(rdef) =  err_gp(rdef)/impest_C
    err_gr(rdef) =  err_gr(rdef)/impest_C     
    err_gpf(rdef)=  err_gpf(rdef)/impest_C     



    ! call expansions with estimated order to save CPU time



    select case (iexp)
 


      case (1)  
        call CalcCg(C_alt,Cuv,p10,p21,p20,m02,m12,m22,rmax,g,g,id,Cerr1_alt,acc_req_Cr,Cerr2_alt)

        Cerr_alt = Cerr2_alt

        CCount(3) = CCount(3)+1
        CrCalc(0:rmax)=CrCalc(0:rmax)+4
        CCalc=CCalc+4
        Crmethod_alt(0:rmax)=4
        


        err_g=Cerr_alt

        call CopyCimp3(C,C_alt,Cerr,Cerr_alt,Cerr1,Cerr1_alt,Cerr2,Cerr2_alt,Crmethod,Crmethod_alt,rmax,rmax)
        


      case (2)

        call CalcCgy(C_alt,Cuv,p10,p21,p20,m02,m12,m22,rmax,gy,gy,id,Cerr1_alt,acc_req_Cr,Cerr2_alt)

        Cerr_alt = Cerr2_alt

        CCount(4) = CCount(4)+1
        CrCalc(0:rmax)=CrCalc(0:rmax)+8
        CCalc=CCalc+8
        Crmethod_alt(0:rmax)=8



        err_gy=Cerr_alt

        call CopyCimp3(C,C_alt,Cerr,Cerr_alt,Cerr1,Cerr1_alt,Cerr2,Cerr2_alt,Crmethod,Crmethod_alt,rmax,rmax)



      case (3)  
        call CalcCgp(C_alt,Cuv,p10,p21,p20,m02,m12,m22,rmax,gp,gp,id,Cerr1_alt,acc_req_Cr,Cerr2_alt)

        Cerr_alt = Cerr2_alt

        CCount(5) = CCount(5)+1
        CrCalc(0:rmax)=CrCalc(0:rmax)+16
        CCalc=CCalc+16
        Crmethod_alt(0:rmax)=16



        err_gp=Cerr_alt

        call CopyCimp3(C,C_alt,Cerr,Cerr_alt,Cerr1,Cerr1_alt,Cerr2,Cerr2_alt,Crmethod,Crmethod_alt,rmax,rmax)



      case (4)  
        call CalcCgr(C_alt,Cuv,p10,p21,p20,m02,m12,m22,rmax,gr,gr,id,Cerr1_alt,acc_req_Cr,Cerr2_alt)

        Cerr_alt = Cerr2_alt

        CCount(6) = CCount(6)+1
        CrCalc(0:rmax)=CrCalc(0:rmax)+32
        CCalc=CCalc+32
        Crmethod_alt(0:rmax)=32



        err_gr=Cerr_alt

        call CopyCimp3(C,C_alt,Cerr,Cerr_alt,Cerr1,Cerr1_alt,Cerr2,Cerr2_alt,Crmethod,Crmethod_alt,rmax,rmax)



      case (5)

        call CalcCgpf(C_alt,Cuv,p10,p21,p20,m02,m12,m22,rmax,gpf,gpf,id,Cerr1_alt,acc_req_Cr,Cerr2_alt)

        Cerr_alt = Cerr2_alt

        CCount(4) = CCount(4)+1
        CrCalc(0:rmax)=CrCalc(0:rmax)+8
        CCalc=CCalc+8
        Crmethod_alt(0:rmax)=8



        err_gpf=Cerr_alt

        call CopyCimp3(C,C_alt,Cerr,Cerr_alt,Cerr1,Cerr1_alt,Cerr2,Cerr2_alt,Crmethod,Crmethod_alt,rmax,rmax)



    end select



    if (iexp.ne.0)  then         !  if added 21.11.2016 
      if (rmax.ge.1) then
        Ctyp =  max(abs(C(0,0,0)),abs(C(0,1,0)),abs(C(0,0,1)))
      else
        Ctyp =  abs(C(0,0,0))
      end if
      err_req_Cr = acc_req_Cr * Ctyp
      
      if (maxval(Cerr1-err_req_Cr).lt.0) then
        CCount(CCalc+CCountoffset0) = CCount(CCalc+CCountoffset0)+1
        return
      end if
    end if    






200 continue

    ! try PV with shifted momentum
    shiftloop: do i=1,2



      if (i.eq.1) then
        mm02shift = mm12
        mm12shift = mm02
        mm22shift = mm22
        q10shift = q10
        q21shift = q20
        q20shift = q21
      else 
        mm02shift = mm22
        mm12shift = mm12
        mm22shift = mm02
        q10shift = q21
        q21shift = q10
        q20shift = q20
      end if

      Zshift(1,1) = 2d0*q10shift
      Zshift(2,1) = q10shift+q20shift-q21shift
      Zshift(1,2) = Zshift(2,1)
      Zshift(2,2) = 2d0*q20shift

      maxZshift = maxval(abs(Zshift))

      detZshift = chdet(2,Zshift)
! added 16.08.2018
    if (abs(detZshift).lt.dprec_cll*maxval(abs(Zshift(1,:)))*maxval(abs(Zshift(2,:)))) then

      detZshift = 0d0
    end if

!      if (detZshift.ne.0d0) then
!        call chinv(2,Zshift,Zinvshift)
!        Zadjshift = Zinvshift * detZshift
!      else
        Zadjshift(1,1) = Zshift(2,2)
        Zadjshift(2,1) = -Zshift(1,2)
        Zadjshift(1,2) = -Zshift(1,2)
        Zadjshift(2,2) = Zshift(1,1)
!      end if

      Zadjsshift(1) = q21shift + q20shift - q10shift 
      Zadjsshift(2) = q21shift + q10shift - q20shift

      detZmZadjfshift = q21shift*Zshift(2,1)

      adetZshift = abs(detZshift)
      maxZadjshift = max(abs(Zadjshift(1,1)),abs(Zadjshift(2,1)),abs(Zadjshift(2,2)))

      fshift(1) = q10shift+mm02shift-mm12shift 
      fshift(2) = q20shift+mm02shift-mm22shift



      mxshift(0,0) = 2d0*mm02shift
      mxshift(1,0) = fshift(1)
      mxshift(2,0) = fshift(2)
      mxshift(0,1) = mxshift(1,0)
      mxshift(0,2) = mxshift(2,0)
!     mxshift(1,1) = 2d0*q10shift
!     mxshift(2,1) = q10shift+q20shift-q21shift
!     mxshift(2,2) = 2d0*q20shift
!     mxshift(1,2) = mxshift(2,1)

      mxshift(1:2,1:2) = Zshift(1:2,1:2)

! changed 21.06.2018
! deterr added 17.01.2019
      call chinve(3,mxshift,mxinvshift,detXshift,deterr)
!      detXshift = chdet(3,mxshift)

!    write(*,*) 'reductionC detXshift = ',detXshift,dprec_cll/deterr

! added 17.01.2019
    if (deterr.lt.dprec_cll) detXshift = 0d0



      if (detXshift.ne.0d0.and.maxZshift.ne.0d0) then

!     write(*,*) 'CalcCred mxshift=',mxshift

!        call chinv(3,mxshift,mxinvshift)

!     write(*,*) 'CalcCred mxinvshift=',mxinvshift

        Xadjshift = mxinvshift * detXshift

!     write(*,*) 'CalcCred Xadj=',Xadj

        Zadjfshift(1:2) = -Xadjshift(0,1:2)

      else
        Zadjfshift(1) = Zadjshift(1,1)*fshift(1)+Zadjshift(2,1)*fshift(2)
        Zadjfshift(2) = Zadjshift(1,2)*fshift(1)+Zadjshift(2,2)*fshift(2)
        Xadjshift(2,2) = 2d0*mm02shift*mxshift(1,1) - fshift(1)*fshift(1)
        Xadjshift(1,1) = 2d0*mm02shift*mxshift(2,2) - fshift(2)*fshift(2)
        Xadjshift(2,1) = 2d0*mm02shift*mxshift(1,2) - fshift(1)*fshift(2)
        Xadjshift(1,2) = Xadjshift(2,1)
      end if


      
      maxZadjfshift = max(abs(Zadjfshift(1)),abs(Zadjfshift(2)))
!    maxZadjfds = max(maxZadjfshift,adetZshift)

      aZadjffshift = abs(Zadjfshift(1)*fshift(1) + Zadjfshift(2)*fshift(2))
!    adetXshift = abs(2d0*mm02*detZshift - Zadjfshift(1)*fshift(1) - Zadjfshift(2)*fshift(2))
!    maxXadjshift = max(abs(Xadjshift(1,1)),abs(Xadjshift(2,1)),abs(Xadjshift(2,2)))

      h_pvs = real(undefined_C)
      w_pvs = real(undefined_C)
      v_pvs = real(undefined_C)
      z_pvs = real(undefined_C)
!   16.08.2018 changed, since detZ=0 set above if too small
!      if (adetZshift.lt.dprec_cll*maxZadjfshift.or.adetZshift.lt.dprec_cll*maxZshift**2.or.adetZshift.eq.0d0) then
      if (adetZshift.lt.dprec_cll*maxZadjfshift.or.adetZshift.eq.0d0) then
        use_pvs = .false.
        err_pvs = err_inf
      else
        use_pvs = .true.
        err_pvs(0) = err_C0
        if (rdef.gt.0) then

          h_pvs = sqrt(adetZshift)/maxZadjshift
          w_pvs = max((maxZadjfshift*h_pvs/adetZshift)**2, abs(mm02shift)*maxZshift*h_pv/adetZshift,  &
              aZadjffshift*maxZshift*(h_pvs/adetZshift)**2)
          v_pvs = maxZadjfshift*h_pvs/adetZshift
          z_pvs = maxZ*h_pvs/adetZshift




          if (mod(rdef,2).eq.1) then
            err_pvs(rdef) = max( w_pvs**((rdef-1)/2) * v_pvs * err_C0,  &
                max(w_pvs**((rdef-1)/2),1d0) * z_pvs * err_B )
            

            
          else
            err_pvs(rdef) = max( w_pvs**(rdef/2) * err_C0,  &
                max(w_pvs**(rdef/2-1) * v_pvs, 1d0) * z_pvs * err_B )
            

            
          end if
        end if
      end if



      if(use_pvs.and.err_pvs(rdef).lt. min(err_pv(rdef),err_pv2(rdef),err_g(rdef)  &
        ,err_gy(rdef),err_gp(rdef),err_gr(rdef),err_gpf(rdef)) ) then



        ! use shifted PV-reduction
        if (i.eq.1) then
          call CalcCpvshift(C_alt,Cuv,p10,p20,p21,m12,m02,m22,rmax,id,Cerr1_alt,Cerr2_alt)
! map coefficients back, order of calculation matters!
          do r=1,rmax
            do n2=0,rmax-r
              do n1=rmax-n2,r,-1
                n0 = rmax-n1-n2
!                write(*,*) 'pvs2',n0,n1,n2,-C_alt(0:n0,n1-1,n2),-C_alt(0:n0,n1,n2),-C_alt(0:n0,n1-1,n2+1)
                C_alt(0:n0,n1,n2) = -C_alt(0:n0,n1-1,n2)-C_alt(0:n0,n1,n2)-C_alt(0:n0,n1-1,n2+1)
!                write(*,*) 'pvs2',n0,n1,n2,C_alt(0:n0,n1,n2)
              end do
            end do
          end do

        elseif (i.eq.2) then
          call CalcCpvshift(C_alt,Cuv,p21,p10,p20,m22,m12,m02,rmax,id,Cerr1_alt,Cerr2_alt)

          do r=1,rmax
            do n1=0,rmax-r
              do n2=rmax-n1,r,-1
                n0 = rmax-n1-n2
                C_alt(0:n0,n1,n2) = -C_alt(0:n0,n1,n2-1)-C_alt(0:n0,n1,n2)-C_alt(0:n0,n1+1,n2-1)
!                write(*,*) 'pvs2',n0,n1,n2,C_alt(0:n0,n1,n2),-C_alt(0:n0,n1,n2-1),-C_alt(0:n0,n1,n2),-C_alt(0:n0,n1+1,n2-1)
              end do
            end do
          end do

        end if


        Cerr_alt = Cerr2_alt

        CCount(9) = CCount(9)+1
!        CrCalc(0:rmax)=CrCalc(0:rmax)+1
!        CCalc=CCalc+1
        Crmethod_alt(0:rmax)=1
        if (Cerr_alt(rmax).lt.Cerr(rmax)) then
           CCount(8) = CCount(8)+1
        end if



        err_pvs=Cerr_alt

        call CopyCimp3(C,C_alt,Cerr,Cerr_alt,Cerr1,Cerr1_alt,Cerr2,Cerr2_alt,Crmethod,Crmethod_alt,rmax,rmax)

        if (rmax.ge.1) then
          Ctyp =  max(abs(C(0,0,0)),abs(C(0,1,0)),abs(C(0,0,1)))
        else
          Ctyp =  abs(C(0,0,0))
        end if
        err_req_Cr = acc_req_Cr * Ctyp
      


! check with Cerr = Cerr2: exp not tried => larger errors
        if (maxval(Cerr1-err_req_Cr).lt.0) then
          CCount(CCalc+CCountoffset0) = CCount(CCalc+CCountoffset0)+1
          return
        end if

      end if
    end do shiftloop






! no method works
    if(use_C0.or.use_pv.or.use_pv2.or.use_g.or.use_gy.or.use_gp.or.use_gr.or.use_gpf.or.use_pvs.eqv..false.) then

      call SetErrFlag_coli(-6)
      call ErrOut_coli('CalcCred',' no reduction method works', &
           errorwriteflag)
!      write(nerrout_coli,'((a))')  '  no reduction method works'
      if (errorwriteflag) then
        write(nerrout_coli,fmt10) ' CalcCred: p10 = ',p10
        write(nerrout_coli,fmt10) ' CalcCred: p21 = ',p21
        write(nerrout_coli,fmt10) ' CalcCred: p20 = ',p20
        write(nerrout_coli,fmt10) ' CalcCred: m02 = ',m02
        write(nerrout_coli,fmt10) ' CalcCred: m12 = ',m12
        write(nerrout_coli,fmt10) ' CalcCred: m22 = ',m22   
      end if
      C = 0d0
      Cuv = 0d0
      Cerr = err_inf
      Cerr2 = err_inf



      return
    endif





    ! no method does work optimal
    ! use the least problematic (for each rank)

    do r=rmax,0,-1
     

      if(use_pv.and.mod(CrCalc(r),2).ne.1) then   
    ! estimate accuracy of PV-reduction if not yet calculated
        if (use_pv) then

!          write(*,*) 'CalcCred err_pv',  r,w_pv,v_pv,z_pv,err_C0,err_B

          if (mod(r,2).eq.1) then
            err_pv(r) = max( w_pv**((r-1)/2) * v_pv * err_C0,  &
                max(w_pv**((r-1)/2),1d0) * z_pv * err_B )

!          write(*,*) 'CalcCred err_pv',  w_pv**((r-1)/2) * v_pv * err_C0,  &
!                w_pv**((r-1)/2) * z_pv * err_B, err_C0,err_B

          else if (r.ne.0) then
            err_pv(r) = max( w_pv**(r/2) * err_C0,  &
                max(w_pv**(r/2-1) * v_pv , 1d0) * z_pv * err_B )
                
!          write(*,*) 'CalcCred err_pv', w_pv**(rmax/2) * err_C0, &
!                        w_pv**(rmax/2-1) * v_pv * z_pv * err_B, err_C0,err_B
          else
            err_pv(r) = err_C0
          end if
        else
          err_pv(r) = err_inf
        end if
      ! scale estimates down to allow trying other methods
        err_pv(r)  = err_pv(r)/impest_C
      end if

      if (use_pv2.and.mod(CrCalc(r),4)-mod(CrCalc(r),2).ne.2) then
    ! estimate accuracy of alternative PV-reduction if not yet calculated
        if (use_pv2) then
          if (mod(r,2).eq.1) then
! change 21.10.15 for ! default
!            err_pv2(r) = max( err_C0 * max(w_pv2**r,w_pv2*v_pv2**((r-1)/2) ),  &
!                err_B * z_pv2 * max(w_pv2**(r+1),w_pv2, &
!                                    w_pv2*v_pv2**((r-1)/2),w_pv2**2, & 
!                                    v_pv2**((r+1)/2),v_pv2) )
            err_pv2(r) = max( err_C0 * max(hw_pv2**r,hw_pv2*v_pv2**((r-1)/2) ),  &
                err_B * z_pv2 * max(w_pv2*hw_pv2**(r),hw_pv2, &
                                    hw_pv2*v_pv2**((r-1)/2),w_pv2*hw_pv2, & 
                                    w_pv2*hw_pv2*v_pv2**((r-1)/2), & 
                                    v_pv2**((r+1)/2),v_pv2) )

!            write(*,*) 'CalcC err_pv2 ',r, err_pv2(r), &
!                 err_C0 * max(1d0,w_pv2**r,v_pv2**((r-1)/2),w_pv2*v_pv2**((r-1)/2) )  ,  &
!                 err_B * max(1d0,z_pv2*w_pv2**(r+1),z_pv2*w_pv2, &
!                z_pv2*w_pv2*v_pv2**((r-1)/2),z_pv2*w_pv2**2, & 
!                z_pv2*v_pv2**((r+1)/2),z_pv2*v_pv2)

          else
! change 21.10.15 for ! default
!            err_pv2(r) = max( err_C0 * max(w_pv2**r,v_pv2**(r/2)),  &
!                err_B * z_pv2 * max(w_pv2**(r+1),w_pv2,  &
!                                    w_pv2*v_pv2**(r/2), w_pv2**2,  &
!                                    v_pv2**(r/2),v_pv2) )
            err_pv2(r) = max( err_C0 * max(hw_pv2**r,v_pv2**(r/2)),  &
                err_B * z_pv2 * max(w_pv2*hw_pv2**(r),hw_pv2,  &
                                    hw_pv2*v_pv2**(r/2), w_pv2*hw_pv2,  &
                                    v_pv2**(r/2),v_pv2) )
          end if
        else
          err_pv2(r) = err_inf
        end if
      ! scale estimates down to allow trying other methods
        err_pv2(r) = err_pv2(r)/impest_C     
      end if

      if (use_g.and.mod(CrCalc(r),8)-mod(CrCalc(r),4).ne.4) then
      ! estimate accuracy of alternative Gram expansion if not yet calculated
        err_g_B(r) = err_B * u_g**r * z_g
        err_g_exp = u_g**(r-1) * Ctyp
        err_g(r) = err_inf 

      ! determine optimal order of expansion 
        do i=0,rmax_C-r
          g = i
          err_g_exp = err_g_exp*fac_g
          err_g(r) = max(err_g_exp,err_g_B(r))

!          write(*,*) 'CalcCred gi',i,g,err_g_exp,err_g(r),err_g_B(r),err_req_Cr(r)

          if (err_g_exp.lt.err_g_B(r).or.err_g(r).lt.err_req_Cr(r)) exit
        end do
      ! increase gp by 2 to account for bad estimates
        g = min(max(g+2,2*g),rmax_C-r)
      ! scale estimates down to allow trying other methods
        err_g(r)  =  err_g(r)/impest_C
      end if

      if (use_gy.and.mod(CrCalc(r),16)-mod(CrCalc(r),8).ne.8) then
      ! estimate accuracy of alternative Gram-Cayley expansion if not yet calculated
        err_gy_B(r) = err_B * b_gy*v1_gy
        err_gy_exp = 1d0 * Ctyp

      ! determine optimal order of expansion 
        gy = 0
        do i=0,rmax_C-r
          if (mod(i,2).eq.1) then
            gy = i/2
            err_gy_exp = err_gy_exp*fac_gy
            err_gy(r) = max(err_gy_exp, err_gy_B(r))
            if (err_gy_exp.lt.err_gy_B(r).or.err_gy(r).lt.err_req_Cr(r)) exit
          end if
        end do
      ! increase gy by 2 to account for bad estimates
        gy = min(max(gy+4,2*gy),(rmax_C-r)/2)
      ! scale estimates down to allow trying other methods
        err_gy(r) =  err_gy(r)/impest_C     
      end if

      if (use_gp.and.mod(CrCalc(r),32)-mod(CrCalc(r),16).ne.16) then
      ! estimate accuracy of small momentum expansion if not yet calculated
        err_gp_B(r) = err_B * z_gp*v_gp**r
        err_gp_exp = v_gp**(r-1) * Ctyp 

      ! determine optimal order of expansion 
        do i=0,rmax_C-r
          gp = i
          err_gp_exp = err_gp_exp*fac_gp
          err_gp(r) = max(err_gp_exp,err_gp_B(r))
          if (err_gp_exp.lt.err_gp_B(r).or.err_gp(r).lt.err_req_Cr(r)) exit
        end do
      ! increase gp by 2 to account for bad estimates
        gp = min(max(gp+2,2*gp),rmax_C-r)
      ! scale estimates down to allow trying other methods
        err_gp(r) =  err_gp(r)/impest_C
      end if

      if (mod(CrCalc(r),64)-mod(CrCalc(r),32).ne.32.and.use_gr) then
      ! estimate accuracy of alternative Gram expansion
        err_gr_B(r) = err_B * a_gr
        err_gr_exp = y1_gr * Ctyp

      ! determine optimal order of expansion 
        gr = 0
        do i=0,rmax_C-r
          if (mod(i,2).eq.1) then
            gr = i/2
            err_gr_exp = err_gr_exp*fac_gr
            err_gr(r) = max(err_gr_exp,err_gr_B(r))



            if (err_gr_exp.lt.err_gr_B(r).or.err_gr(r).lt.err_req_Cr(r)) exit
          end if
        end do
      ! increase gr to account for bad estimates
! changed 28.07.14
!       gr = min(max(gr+2,2*gr),(rmax_C-r)/2)
        gr = min(max(gr+2,2*gr),rmax_C-r,max(0,(rmax_B-2*r)/2))
      ! scale estimates down to allow trying other methods
        err_gr(r) =  err_gr(r)/impest_C     

      end if

      if (use_gpf.and.mod(CrCalc(r),128)-mod(CrCalc(r),64).ne.64) then
      ! estimate accuracy of alternative Gram-Cayley expansion if not yet calculated
        err_gpf_B(r) = err_B * b_gpf*v1_gpf
        err_gpf_exp = 1d0 * Ctyp

      ! determine optimal order of expansion 
        gpf = 0
        do i=0,rmax_C-r
          if (mod(i,2).eq.1) then
            gpf = i/2
            err_gpf_exp = err_gpf_exp*fac_gpf
            err_gpf(r) = max(err_gpf_exp, err_gpf_B(r))
            if (err_gpf_exp.lt.err_gpf_B(r).or.err_gpf(r).lt.err_req_Cr(r)) exit
          end if
        end do
      ! increase gpf by 2 to account for bad estimates
        gpf = min(max(gpf+4,2*gpf),(rmax_C-r)/2)
      ! scale estimates down to allow trying other methods
        err_gpf(r) =  err_gpf(r)/impest_C     
      end if




100   continue   ! try other methods if error larger than expected

      if (min(err_pv(r),err_pv2(r)).le.min(err_g(r),err_gy(r),err_gp(r),err_gr(r),err_gpf(r))       &
          .and.min(err_pv(r),err_pv2(r)).lt.err_inf) then

        if (use_pv.and.err_pv(r).le.err_pv2(r).and.mod(CrCalc(r),2).ne.1) then

!          deallocate(C_alt)
!          deallocate(Cuv_alt)
!          deallocate(Cerr_alt)
!          deallocate(Cerr2_alt)
!          deallocate(Crmethod_alt)
!          allocate(C_alt(0:r,0:r,0:r))
!          allocate(Cuv_alt(0:r,0:r,0:r))
!          allocate(Cerr_alt(0:r))
!          allocate(Cerr2_alt(0:r))
!          allocate(Crmethod_alt(0:r))


          if (r.eq.rmax) then
            call CalcCpv1(C_alt,Cuv,p10,p21,p20,m02,m12,m22,r,id,Cerr1_alt,Cerr2_alt)
          else
            call CalcCpv1(C_alt(0:r,0:r,0:r),Cuv_alt(0:r,0:r,0:r),  &
                p10,p21,p20,m02,m12,m22,r,id,Cerr1_alt(0:r),Cerr2_alt(0:r))
          end if

          Cerr_alt = Cerr2_alt

          CCount(11) = CCount(11)+1
          CrCalc(0:r)=CrCalc(0:r)+1
          CCalc=CCalc+1
          Crmethod_alt(0:r)=1
          checkest=Cerr_alt(r)/(err_pv(r)*abs(C_alt(0,0,0)))
          


          err_pv(0:r)=Cerr_alt(0:r)

          call CopyCimp3(C,C_alt(0:r,0:r,0:r),Cerr,Cerr_alt(0:r),Cerr1,Cerr1_alt(0:r),   &   
              Cerr2,Cerr2_alt(0:r),Crmethod,Crmethod_alt(0:r),rmax,r)

          if (rmax.ge.1) then
            Ctyp =  max(abs(C(0,0,0)),abs(C(0,1,0)),abs(C(0,0,1)))
          else
            Ctyp =  abs(C(0,0,0))
          end if
          err_req_Cr = acc_req_Cr * Ctyp



          if(checkest.gt.impest_C.and.Mode_coli.lt.1) goto 100     ! error larger than expected: try other methods
          
        elseif (use_pv2.and.err_pv2(r).le.err_pv(r).and.mod(CrCalc(r),4)-mod(CrCalc(r),2).ne.2) then

!          deallocate(C_alt)
!          deallocate(Cuv_alt)
!          deallocate(Cerr_alt)
!          deallocate(Cerr2_alt)
!          deallocate(Crmethod_alt)
!          allocate(C_alt(0:r,0:r,0:r))
!          allocate(Cuv_alt(0:r,0:r,0:r))
!          allocate(Cerr_alt(0:r))
!          allocate(Cerr2_alt(0:r))
!          allocate(Crmethod_alt(0:r))

          if (r.eq.rmax) then
            call CalcCpv2(C_alt,Cuv,p10,p21,p20,m02,m12,m22,r,id,Cerr1_alt,Cerr2_alt)
          else
            call CalcCpv2(C_alt(0:r,0:r,0:r),Cuv_alt(0:r,0:r,0:r),p10,p21,p20,m02,m12,m22,r,id,Cerr1_alt(0:r),Cerr2_alt(0:r))
          end if

          Cerr_alt = Cerr2_alt

          CCount(12) = CCount(12)+1
          CrCalc(0:r)=CrCalc(0:r)+2
          CCalc=CCalc+2
          Crmethod_alt(0:r)=2
          checkest=Cerr_alt(r)/(err_pv(r)*abs(C_alt(0,0,0)))
          


          err_pv2(0:r)=Cerr_alt(0:r)

          call CopyCimp3(C,C_alt(0:r,0:r,0:r),Cerr,Cerr_alt(0:r),Cerr1,Cerr1_alt(0:r),   &   
              Cerr2,Cerr2_alt(0:r),Crmethod,Crmethod_alt(0:r),rmax,r)

          if (rmax.ge.1) then
            Ctyp =  max(abs(C(0,0,0)),abs(C(0,1,0)),abs(C(0,0,1)))
          else
            Ctyp =  abs(C(0,0,0))
          end if
          err_req_Cr = acc_req_Cr * Ctyp


          
          if(checkest.gt.impest_C.and.Mode_coli.lt.1) goto 100     ! error larger than expected: try other methods

        end if

      else 



        if (use_g.and.err_g(r).le.min(err_gy(r),err_gp(r),err_gr(r),err_gpf(r))        &
            .and.mod(CrCalc(r),8)-mod(CrCalc(r),4).ne.4) then

!          deallocate(C_alt)
!          deallocate(Cuv_alt)
!          deallocate(Cerr_alt)
!          deallocate(Cerr2_alt)
!          deallocate(Crmethod_alt)
!          allocate(C_alt(0:r,0:r,0:r))
!          allocate(Cuv_alt(0:r,0:r,0:r))
!          allocate(Cerr_alt(0:r))
!          allocate(Cerr2_alt(0:r))
!          allocate(Crmethod_alt(0:r))

          if (r.eq.rmax) then
            call CalcCg(C_alt,Cuv,p10,p21,p20,m02,m12,m22,r,g,rmax_C,id,Cerr1_alt,acc_req_Cr,Cerr2_alt)
          else
            call CalcCg(C_alt(0:r,0:r,0:r),Cuv_alt(0:r,0:r,0:r),   &
                p10,p21,p20,m02,m12,m22,r,g,rmax_C,id,Cerr1_alt(0:r),acc_req_Cr(0:r),Cerr2_alt(0:r))
          end if

          Cerr_alt = Cerr2_alt

          CCount(13) = CCount(13)+1
          CrCalc(0:r)=CrCalc(0:r)+4
          CCalc=CCalc+4
          Crmethod_alt(0:r)=4
          checkest=Cerr_alt(r)/err_g(r)


          
          err_g(0:r)=Cerr_alt(0:r)

          call CopyCimp3(C,C_alt(0:r,0:r,0:r),Cerr,Cerr_alt(0:r),Cerr1,Cerr1_alt(0:r),   &   
              Cerr2,Cerr2_alt(0:r),Crmethod,Crmethod_alt(0:r),rmax,r)

          if (rmax.ge.1) then
            Ctyp =  max(abs(C(0,0,0)),abs(C(0,1,0)),abs(C(0,0,1)))
          else
            Ctyp =  abs(C(0,0,0))
          end if
          err_req_Cr = acc_req_Cr * Ctyp



          if(checkest.gt.impest_C.and.Mode_coli.lt.1) goto 100     ! error larger than expected: try other methods
          
        else if (use_gy.and.err_gy(r).le.min(err_g(r),err_gp(r),err_gr(r),err_gpf(r))        &
            .and.mod(CrCalc(r),16)-mod(CrCalc(r),8).ne.8) then
!          deallocate(C_alt)
!          deallocate(Cuv_alt)
!          deallocate(Cerr_alt)
!          deallocate(Cerr2_alt)
!          deallocate(Crmethod_alt)
!          allocate(C_alt(0:r,0:r,0:r))
!          allocate(Cuv_alt(0:r,0:r,0:r))
!          allocate(Cerr_alt(0:r))
!          allocate(Cerr2_alt(0:r))
!          allocate(Crmethod_alt(0:r))

          if (r.eq.rmax) then
            call CalcCgy(C_alt,Cuv,p10,p21,p20,m02,m12,m22,r,gy,(rmax_C)/2,id,Cerr1_alt,acc_req_Cr,Cerr2_alt)
          else
            call CalcCgy(C_alt(0:r,0:r,0:r),Cuv_alt(0:r,0:r,0:r),  &
                p10,p21,p20,m02,m12,m22,r,gy,(rmax_C)/2,id,Cerr1_alt(0:r),acc_req_Cr(0:r),Cerr2_alt(0:r))
          end if

          Cerr_alt = Cerr2_alt

          CCount(14) = CCount(14)+1
          CrCalc(0:r)=CrCalc(0:r)+8
          CCalc=CCalc+8
          Crmethod_alt(0:r)=8
          checkest=Cerr_alt(r)/err_gy(r)


          err_gy(0:r)=Cerr_alt(0:r)
          
          call CopyCimp3(C,C_alt(0:r,0:r,0:r),Cerr,Cerr_alt(0:r),Cerr1,Cerr1_alt(0:r),   &   
              Cerr2,Cerr2_alt(0:r),Crmethod,Crmethod_alt(0:r),rmax,r)
  
          if (rmax.ge.1) then
            Ctyp =  max(abs(C(0,0,0)),abs(C(0,1,0)),abs(C(0,0,1)))
          else
            Ctyp =  abs(C(0,0,0))
          end if
          err_req_Cr = acc_req_Cr * Ctyp



          if(checkest.gt.impest_C.and.Mode_coli.lt.1) goto 100     ! error larger than expected: try other methods
          
        elseif (use_gp.and.err_gp(r).le.min(err_g(r),err_gy(r),err_gr(r),err_gpf(r))        &
            .and.mod(CrCalc(r),32)-mod(CrCalc(r),16).ne.16) then

!          deallocate(C_alt)
!          deallocate(Cuv_alt)
!          deallocate(Cerr_alt)
!          deallocate(Cerr2_alt)
!          deallocate(Crmethod_alt)
!          allocate(C_alt(0:r,0:r,0:r))
!          allocate(Cuv_alt(0:r,0:r,0:r))
!          allocate(Cerr_alt(0:r))
!          allocate(Cerr2_alt(0:r))
!          allocate(Crmethod_alt(0:r))

          if (r.eq.rmax) then
            call CalcCgp(C_alt,Cuv,p10,p21,p20,m02,m12,m22,r,gp,rmax_C,id,Cerr1_alt,acc_req_Cr,Cerr2_alt)
          else
            call CalcCgp(C_alt(0:r,0:r,0:r),Cuv_alt(0:r,0:r,0:r),   &
                p10,p21,p20,m02,m12,m22,r,gp,rmax_C,id,Cerr1_alt(0:r),acc_req_Cr(0:r),Cerr2_alt(0:r))
          end if

          Cerr_alt = Cerr2_alt

          CCount(15) = CCount(15)+1
          CrCalc(0:r)=CrCalc(0:r)+16
          CCalc=CCalc+16
          Crmethod_alt(0:r)=16
          checkest=Cerr_alt(r)/err_gp(r)



          err_gp(0:r)=Cerr_alt(0:r)
          
          call CopyCimp3(C,C_alt(0:r,0:r,0:r),Cerr,Cerr_alt(0:r),Cerr1,Cerr1_alt(0:r),   &   
              Cerr2,Cerr2_alt(0:r),Crmethod,Crmethod_alt(0:r),rmax,r)

          if (rmax.ge.1) then
            Ctyp =  max(abs(C(0,0,0)),abs(C(0,1,0)),abs(C(0,0,1)))
          else
            Ctyp =  abs(C(0,0,0))
          end if
          err_req_Cr = acc_req_Cr * Ctyp



          if(checkest.gt.impest_C.and.Mode_coli.lt.1) goto 100     ! error larger than expected: try other methods          

        elseif (use_gr.and.err_gr(r).le.min(err_g(r),err_gy(r),err_gp(r),err_gpf(r))   &
            .and.mod(CrCalc(r),64)-mod(CrCalc(r),32).ne.32) then
!          deallocate(C_alt)
!          deallocate(Cuv_alt)
!          deallocate(Cerr_alt)
!          deallocate(Cerr2_alt)
!          deallocate(Crmethod_alt)
!          allocate(C_alt(0:r,0:r,0:r))
!          allocate(Cuv_alt(0:r,0:r,0:r))
!          allocate(Cerr_alt(0:r))
!          allocate(Cerr2_alt(0:r))
!          allocate(Crmethod_alt(0:r))

          if (r.eq.rmax) then
            call CalcCgr(C_alt,Cuv,p10,p21,p20,m02,m12,m22,r,gr,rmax_C,id,Cerr1_alt,acc_req_Cr,Cerr2_alt)
          else
            call CalcCgr(C_alt(0:r,0:r,0:r),Cuv_alt(0:r,0:r,0:r),   &
                p10,p21,p20,m02,m12,m22,r,gr,rmax_C,id,Cerr1_alt(0:r),acc_req_Cr(0:r),Cerr2_alt(0:r))
          end if

          Cerr_alt = Cerr2_alt

          CCount(16) = CCount(16)+1
          CrCalc(0:r)=CrCalc(0:r)+32
          CCalc=CCalc+32
          Crmethod_alt(0:r)=32
          checkest=Cerr_alt(r)/err_gr(r)



          err_gr(0:r)=Cerr_alt(0:r)
          
          call CopyCimp3(C,C_alt(0:r,0:r,0:r),Cerr,Cerr_alt(0:r),Cerr1,Cerr1_alt(0:r),   &   
              Cerr2,Cerr2_alt(0:r),Crmethod,Crmethod_alt(0:r),rmax,r)

          if (rmax.ge.1) then
            Ctyp =  max(abs(C(0,0,0)),abs(C(0,1,0)),abs(C(0,0,1)))
          else
            Ctyp =  abs(C(0,0,0))
          end if
          err_req_Cr = acc_req_Cr * Ctyp



          if(checkest.gt.impest_C.and.Mode_coli.lt.1) goto 100     ! error larger than expected: try other methods          

        else if (use_gpf.and.err_gpf(r).le.min(err_g(r),err_gy(r),err_gp(r),err_gr(r))        &
            .and.mod(CrCalc(r),128)-mod(CrCalc(r),64).ne.64) then
!          deallocate(C_alt)
!          deallocate(Cuv_alt)
!          deallocate(Cerr_alt)
!          deallocate(Cerr2_alt)
!          deallocate(Crmethod_alt)
!          allocate(C_alt(0:r,0:r,0:r))
!          allocate(Cuv_alt(0:r,0:r,0:r))
!          allocate(Cerr_alt(0:r))
!          allocate(Cerr2_alt(0:r))
!          allocate(Crmethod_alt(0:r))

          if (r.eq.rmax) then
            call CalcCgpf(C_alt,Cuv,p10,p21,p20,m02,m12,m22,r,gpf,(rmax_C)/2,id,Cerr1_alt,acc_req_Cr,Cerr2_alt)
          else
            call CalcCgpf(C_alt(0:r,0:r,0:r),Cuv_alt(0:r,0:r,0:r),  &
                p10,p21,p20,m02,m12,m22,r,gpf,(rmax_C)/2,id,Cerr1_alt(0:r),acc_req_Cr(0:r),Cerr2_alt(0:r))
          end if

          Cerr_alt = Cerr2_alt

          CCount(17) = CCount(17)+1
          CrCalc(0:r)=CrCalc(0:r)+64
          CCalc=CCalc+64
          Crmethod_alt(0:r)=64
          checkest=Cerr_alt(r)/err_gpf(r)


          err_gpf(0:r)=Cerr_alt(0:r)
          
          call CopyCimp3(C,C_alt(0:r,0:r,0:r),Cerr,Cerr_alt(0:r),Cerr1,Cerr1_alt(0:r),   &   
              Cerr2,Cerr2_alt(0:r),Crmethod,Crmethod_alt(0:r),rmax,r)
  
          if (rmax.ge.1) then
            Ctyp =  max(abs(C(0,0,0)),abs(C(0,1,0)),abs(C(0,0,1)))
          else
            Ctyp =  abs(C(0,0,0))
          end if
          err_req_Cr = acc_req_Cr * Ctyp



          if(checkest.gt.impest_C.and.Mode_coli.lt.1) goto 100     ! error larger than expected: try other methods
          
        end if

      end if




    end do

    norm = abs(C(0,0,0))
    do r=1,rdef
      do n1=0,rdef
        n2 = rdef-n1
       norm =  max(norm,abs(C(0,n1,n2)))
      end do
    end do
    acc_C(0:rdef) = Cerr(0:rdef)/norm

    CCount(CCalc+CCountoffset0) = CCount(CCalc+CCountoffset0)+1

    if (maxval(acc_C(0:rdef)-sqrt(reqacc_coli)).gt.0) then
      CCount(CCalc+CCountoffset3) = CCount(CCalc+CCountoffset3)+1
    end if

    if (maxval(acc_C(0:rdef)-reqacc_coli).gt.0) then
      CCount(CCalc+CCountoffset1) = CCount(CCalc+CCountoffset1)+1
    end if



    if (maxval(acc_C(0:rdef)-critacc_coli).gt.0) then

      CCount(CCalc+CCountoffset2) = CCount(CCalc+CCountoffset2)+1



!      call SetErrFlag_coli(-5)
!      call ErrOut_coli('CalcCred',' critical accuracy not reached', &
!           errorwriteflag)


    end if



  end subroutine CalcCred
 



  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !  subroutine CalcCuv(Cuv,Buv_0,m02,f,rmax,id)
  !
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  
  subroutine CalcCuv(Cuv,Buv_0,m02,f,rmax,id)
  
    integer, intent(in) :: rmax,id
    double complex, intent(in) :: m02,f(2)
!   double complex, intent(inout) :: Cuv(0:rmax,0:rmax,0:rmax)
    double complex, intent(out) :: Cuv(0:rmax,0:rmax,0:rmax)
    double complex, intent(in) :: Buv_0(0:rmax-1,0:rmax-1,0:rmax-1)
    integer :: r,n0,n1,n2,r0
        
    ! C_(0,n1,n2) UV-finite
    Cuv(0,:,:) = 0d0
    
!    do r=2,rmax
!      do n0=1,rmax/2
!        do n1=0,r-2*n0
!          n2 = r-2*n0-n1
    do r=2,2*rmax
      do n0=max(1,r-rmax),r/2
        do n1=0,r-2*n0
          n2 = r-2*n0-n1
        
          Cuv(n0,n1,n2) = (Buv_0(n0-1,n1,n2) + 2*m02*Cuv(n0-1,n1,n2)  &
                          + f(1)*Cuv(n0-1,n1+1,n2)  & 
                          + f(2)*Cuv(n0-1,n1,n2+1)) / (2*r) 

        end do
      end do
    end do
    
  end subroutine CalcCuv






  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !  subroutine CalcCpv1(C,Cuv,p10,p21,p20,m02,m12,m22,rmax,Cerr,Cerr2)
  !
  !  new version 10.02.2016 (5.10) with (5.11) inserted
  !              27.09.2016    prefactors of B_0 improved
  !              02.09.2017    allocate and q1q2 removed  
  !
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  
  subroutine CalcCpv1(C,Cuv,p10,p21,p20,m02,m12,m22,rmax,id,Cerr,Cerr2)
  
    use globalC
  
    integer, intent(in) :: rmax,id
    double complex, intent(in) :: p10,p21,p20,m02,m12,m22
    double complex, intent(out) :: Cuv(0:rmax,0:rmax,0:rmax)
    double complex, intent(out) :: C(0:rmax,0:rmax,0:rmax)
    double precision, intent(out) :: Cerr(0:rmax),Cerr2(0:rmax)
!   double complex, allocatable :: B_0(:,:,:), Buv_0(:,:,:)
!   double complex, allocatable :: B_i(:,:,:), Buv_i(:,:,:)
!   double complex, allocatable :: C_alt(:,:,:)
    double complex :: B_0(0:rmax-1,0:rmax-1,0:rmax-1), Buv_0(0:rmax-1,0:rmax-1,0:rmax-1)
    double complex :: B_i(0:rmax-1,0:rmax-1,2), Buv_i(0:rmax-1,0:rmax-1,2)
    double complex :: C_alt(0:rmax,0:rmax,0:rmax)
    double complex :: Smod(2)
    double complex :: C0_coli, elimminf2_coli
!   double precision, allocatable :: C00_err(:),Cij_err(:)
!   double precision, allocatable :: C00_err2(:),Cij_err2(:)
    double precision :: C00_err(0:rmax),Cij_err(0:rmax)
    double precision :: C00_err2(0:rmax),Cij_err2(0:rmax)
    double precision :: B_err,B_max
    integer :: rmaxB,r,n0,n1,n2,nn0,nn1,nn2,i,j
    integer :: bin,k,nid(0:2)



    ! calculation of scalar coefficient
    C(0,0,0) = C0_coli(p10,p21,p20,m02,m12,m22)
    Cuv(0,0,0) = 0d0

    ! accuracy estimate for C0 function
    Cerr(0) = acc_def_C0*max(1d0/sqrt(adetZ),abs(C(0,0,0)))   
    Cerr2(0) = acc_def_C0*max(1d0/sqrt(adetZ),abs(C(0,0,0)))   

    if (rmax.eq.0) return
    
    ! allocation and calculation of B functions
    rmaxB = rmax-1
    ! rmaxB = max(rmax-1,0)
!   allocate(B_0(0:rmaxB,0:rmaxB,0:rmaxB))
!   allocate(Buv_0(0:rmaxB,0:rmaxB,0:rmaxB))
!   allocate(B_i(0:rmaxB,0:rmaxB,2))
!   allocate(Buv_i(0:rmaxB,0:rmaxB,2))

    ! allocate arrays for error propagation
!   allocate(C00_err(0:rmax))
!   allocate(Cij_err(0:rmax))
!   allocate(C00_err2(0:rmax))
!   allocate(Cij_err2(0:rmax))
    
    ! determine binaries for B-coefficients
    k=0
    bin = 1
    do while (k.le.2)
      if (mod(id/bin,2).eq.0) then
        nid(k) = id+bin
        k = k+1
      end if
      bin = 2*bin
    end do

    call CalcB(B_0(:,0,:),Buv_0(:,0,:),p21,m12,m22,rmaxB,nid(0))
    call CalcB(B_i(:,:,1),Buv_i(:,:,1),p20,m02,m22,rmaxB,nid(1))
    call CalcB(B_i(:,:,2),Buv_i(:,:,2),p10,m02,m12,rmaxB,nid(2))

    ! shift of integration momentum in B_0 and calculate maximal B(0,...)
    B_max=0d0
    do n1=1,rmaxB
      do n2=0,rmaxB-n1
        n0 = (rmaxB-n1-n2)
        B_0(0:n0,n1,n2) = -B_0(0:n0,n1-1,n2)-B_0(0:n0,n1-1,n2+1)
        Buv_0(0:n0,n1,n2) = -Buv_0(0:n0,n1-1,n2)-Buv_0(0:n0,n1-1,n2+1)
        B_max = max(B_max,abs(B_0(0,n1,n2)))
      end do
    end do

!    write(*,*) 'B_max=',B_max

    B_max=max(B_max,maxval(abs(B_i(0,0:rmaxB,1:2))))

    ! determine inverse Gram matrix
    ! commented out 2.9.2017
    ! Zinv = Zadj/detZ

    ! calculate Cuv
    call CalcCuv(Cuv,Buv_0,mm02,f,rmax,id)

    ! initialization of error propagation
    Cij_err =0d0
    C00_err =0d0
    Cij_err(0) = Cerr(0)
    B_err = acc_def_B*B_max

    Cij_err2 =0d0
    C00_err2 =0d0
    Cij_err2(0) = Cerr2(0)




!   allocate(C_alt(0:rmax,0:rmax,0:rmax))

    ! PV reduction
    do r=1,rmax

    ! reduction formula (5.10) with (5.11) inserted for n0 >= 1
      do n0=r/2,1,-1
        do n1=0,r-2*n0
          n2 = r-2*n0-n1
          C(n0,n1,n2) = + 4*Cuv(n0,n1,n2) + detX/detZ*C(n0-1,n1,n2)
            C(n0,n1,n2) = C(n0,n1,n2) &
                + (detZmZadjf + Zadjs(1)*(mm12-mm02) + Zadjs(2)*(mm22-mm02) &
                  ) /detZ * B_0(n0-1,n1,n2)
!                + (1d0 - (Zadjf(1)+Zadjf(2))/detZ)* B_0(n0-1,n1,n2)

            if (n1.ge.1) then
              C(n0,n1,n2) = C(n0,n1,n2) &
                  - 2*n1*Zadjf(1)/detZ*C(n0,n1-1,n2)
            else            
              C(n0,n1,n2) = C(n0,n1,n2) &
                  + Zadjf(1)/detZ* B_i(n0-1,n2,1)
            end if
            if (n2.ge.1) then
              C(n0,n1,n2) = C(n0,n1,n2) &
                  - 2*n2*Zadjf(2)/detZ*C(n0,n1,n2-1)
            else            
              C(n0,n1,n2) = C(n0,n1,n2) &
                  + Zadjf(2)/detZ * B_i(n0-1,n1,2)
            end if

            C(n0,n1,n2) = C(n0,n1,n2)  / (2*r) 

!        if(n0.eq.1) then
!          write(*,*) 'Ca(1,n1,n2)=',n1,n2, 4*Cuv(n0,n1,n2) + detX/detZ*C(n0-1,n1,n2)
!          write(*,*) 'Ca(1,n1,n2)=', (detZmZadjf + Zadjs(1)*(mm12-mm02) + Zadjs(2)*(mm22-mm02) &
!                  ) /detZ * B_0(n0-1,n1,n2)
!          write(*,*) 'Ca(1,n1,n2)=', detZmZadjf , Zadjs(1)*(mm12-mm02) , Zadjs(2)*(mm22-mm02) &
!                  ,detZ , B_0(n0-1,n1,n2)
!          write(*,*) 'Ca(1,n1,n2)=', (1d0 - (Zadjf(1)+Zadjf(2))/detZ)* B_0(n0-1,n1,n2)
!          write(*,*) 'Ca(1,n1,n2)=',  + Zadjf(1)/detZ* B_i(n0-1,n2,1)
!          write(*,*) 'Ca(1,n1,n2)=',   + Zadjf(2)/detZ * B_i(n0-1,n1,2)
!        end if

         end do
      end do

    ! reduction formula (5.11) with (5.10) inserted for n0 = 0
!     do n0=(r-1)/2,0,-1
      n0=0 
        do n1=0,r-2*n0
          n2 = r-2*n0-n1

          if (n1.ge.1) then
            nn1 = n1-1
            nn2 = n2
            j = 1
          else
            nn1 = n1
            nn2 = n2-1
            j = 2
          end if
         
!          do i=1,2
!            Smod(i) = -B_0(n0,nn1,nn2)
!          end do
          Smod = 0d0
          
          if (nn1.ge.1) then
            Smod(1) = Smod(1) - 2d0*nn1*C(n0+1,nn1-1,nn2)
          else          
            Smod(1) = Smod(1) + B_i(n0,nn2,1)
          end if

          if (nn2.ge.1) then
            Smod(2) = Smod(2) - 2d0*nn2*C(n0+1,nn1,nn2-1)
          else
            Smod(2) = Smod(2) + B_i(n0,nn1,2)
          end if
          
          C(n0,n1,n2) = (Zadj(1,j)*Smod(1) + Zadj(2,j)*Smod(2)  &
                           - Zadjs(j)*B_0(n0,nn1,nn2) &
                           - Zadjf(j)*C(n0,nn1,nn2))/detZ



        end do
!      end do

      ! determine error from symmetry for n0=0 and n1>=1, n2>=1 
      Cerr(r)=Cerr(r-1)
      Cerr2(r)=Cerr2(r-1)
      n0=0
        do n1=1,r-2*n0-1
          n2 = r-2*n0-n1

          nn1 = n1
          nn2 = n2-1
          j = 2
         
!          do i=1,2
!            Smod(i) = -B_0(n0,nn1,nn2)
!          end do
          Smod = 0
          
          if (nn1.ge.1) then
            Smod(1) = Smod(1) - 2d0*nn1*C(n0+1,nn1-1,nn2)
          else          
            Smod(1) = Smod(1) + B_i(n0,nn2,1)
          end if

          if (nn2.ge.1) then
            Smod(2) = Smod(2) - 2d0*nn2*C(n0+1,nn1,nn2-1)
          else
            Smod(2) = Smod(2) + B_i(n0,nn1,2)
          end if
           
          C_alt(n0,n1,n2) = (Zadj(1,j)*Smod(1) + Zadj(2,j)*Smod(2)  &
                           - Zadjs(j)*B_0(n0,nn1,nn2) &
                           - Zadjf(j)*C(n0,nn1,nn2))/detZ
 
          Cerr(r)=max(Cerr(r),abs(C(n0,n1,n2)-C_alt(n0,n1,n2)))
          Cerr2(r)=max(Cerr2(r),abs(C(n0,n1,n2)-C_alt(n0,n1,n2)))


!            write(*,*) 'CalcCpv1 Cerr',n0,n1,n2, Cerr(r), abs(C(n0,n1,n2)),abs(C_alt(n0,n1,n2))

        end do

      if(r.ge.2)then
! estimate using insertions of (5.11) in (5.10)
        C00_err(r) = max(2*abs(m02)*Cij_err(r-2), B_err,    &
              aZadjff/adetZ*Cij_err(r-2),             &
              maxZadjf/adetZ*max(2*C00_err(r-1),B_err))/(2*r)

!        write(*,*) 'C00errtest',r,abs(m02)*Cij_err(r-2), B_err,    &
!              aZadjff/adetZ*Cij_err(r-2),             &
!              maxZadjf/adetZ*C00_err(r-1),maxZadjf/adetZ*B_err, &
!              C00_err(r)

      else
        C00_err(r) = 0d0
      end if
! estimate using insertions of (5.10) in (5.11)
      Cij_err(r) = max(maxZadjf*Cij_err(r-1),   &
              maxZadj*max(2*C00_err(r),B_err))/adetZ
 
      if(r.ge.2)then
        C00_err2(r) = max(2*abs(m02)*Cij_err2(r-2), B_err,    &
              aZadjff/adetZ*Cij_err2(r-2),             &
              maxZadjf/adetZ*max(2*C00_err(r-1),B_err))/(2*r)

!        write(*,*) 'C00errtest',r,abs(m02)*Cij_err2(r-2), B_err,    &
!              aZadjff/adetZ*Cij_err2(r-2),             &
!              maxZadjf/adetZ*C00_err2(r-1),maxZadjf/adetZ*B_err, &
!              C00_err2(r)

      else
        C00_err2(r) = 0d0
      end if
      Cij_err2(r) = max((maxZadjf/maxZadj)*Cij_err2(r-1),max(2*C00_err2(r),B_err))/sqrt(adetZ)
    end do

      ! reduction formula (5.10) for n0+n1+n2=r, n0>0 
    do r=rmax+1,2*rmax
      do n0=r-rmax,r/2
        do n1=0,r-2*n0
          n2 = r-2*n0-n1
          C(n0,n1,n2) = (B_0(n0-1,n1,n2) + 2*mm02*C(n0-1,n1,n2) + 4*Cuv(n0,n1,n2) &
                        + f(1)*C(n0-1,n1+1,n2) + f(2)*C(n0-1,n1,n2+1)) / (2*r)
        end do
      end do
    end do




    Cerr2 = max(Cerr2,Cij_err2(0:rmax))
    Cerr = max(Cerr,Cij_err(0:rmax))



!   write(*,*) 'CalcCpv1 out',(((C((r-n1-n2)/2,n1,n2),n2=0,r-n1),n1=0,r),r=0,rmax)

!    write(*,*) 'CalcCpv1 Cerr ',Cerr
!    write(*,*) 'CalcCpv1 Cerr2',Cerr2



  end subroutine CalcCpv1




  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !  subroutine CalcCpv1o(C,Cuv,p10,p21,p20,m02,m12,m22,rmax,Cerr,Cerr2)
  !
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  
  subroutine CalcCpv1o(C,Cuv,p10,p21,p20,m02,m12,m22,rmax,id,Cerr,Cerr2)
  
    use globalC
  
    integer, intent(in) :: rmax,id
    double complex, intent(in) :: p10,p21,p20,m02,m12,m22
    double complex, intent(out) :: Cuv(0:rmax,0:rmax,0:rmax)
    double complex, intent(out) :: C(0:rmax,0:rmax,0:rmax)
    double precision, intent(out) :: Cerr(0:rmax),Cerr2(0:rmax)
    double complex, allocatable :: B_0(:,:,:), Buv_0(:,:,:)
    double complex, allocatable :: B_i(:,:,:), Buv_i(:,:,:)
    double complex, allocatable :: C_alt(:,:,:)
    double complex :: Smod(2)
    double complex :: C0_coli, elimminf2_coli
    double precision, allocatable :: C00_err(:),Cij_err(:)
    double precision, allocatable :: C00_err2(:),Cij_err2(:)
    double precision :: B_err,B_max
    integer :: rmaxB,r,n0,n1,n2,nn0,nn1,nn2,i,j
    integer :: bin,k,nid(0:2)



    ! calculation of scalar coefficient
    C(0,0,0) = C0_coli(p10,p21,p20,m02,m12,m22)
    Cuv(0,0,0) = 0d0

    ! accuracy estimate for C0 function
    Cerr(0) = acc_def_C0*max(1d0/sqrt(adetZ),abs(C(0,0,0)))   
    Cerr2(0) = acc_def_C0*max(1d0/sqrt(adetZ),abs(C(0,0,0)))   

    if (rmax.eq.0) return
    
    ! allocation and calculation of B functions
    rmaxB = rmax-1
    ! rmaxB = max(rmax-1,0)
    allocate(B_0(0:rmaxB,0:rmaxB,0:rmaxB))
    allocate(Buv_0(0:rmaxB,0:rmaxB,0:rmaxB))
    allocate(B_i(0:rmaxB,0:rmaxB,2))
    allocate(Buv_i(0:rmaxB,0:rmaxB,2))

    ! allocate arrays for error propagation
    allocate(C00_err(0:rmax))
    allocate(Cij_err(0:rmax))
    allocate(C00_err2(0:rmax))
    allocate(Cij_err2(0:rmax))
    
    ! determine binaries for B-coefficients
    k=0
    bin = 1
    do while (k.le.2)
      if (mod(id/bin,2).eq.0) then
        nid(k) = id+bin
        k = k+1
      end if
      bin = 2*bin
    end do

    call CalcB(B_0(:,0,:),Buv_0(:,0,:),p21,m12,m22,rmaxB,nid(0))
    call CalcB(B_i(:,:,1),Buv_i(:,:,1),p20,m02,m22,rmaxB,nid(1))
    call CalcB(B_i(:,:,2),Buv_i(:,:,2),p10,m02,m12,rmaxB,nid(2))

    ! shift of integration momentum in B_0 and calculate maximal B(0,...)
    B_max=0d0
    do n1=1,rmaxB
      do n2=0,rmaxB-n1
        n0 = (rmaxB-n1-n2)
        B_0(0:n0,n1,n2) = -B_0(0:n0,n1-1,n2)-B_0(0:n0,n1-1,n2+1)
        Buv_0(0:n0,n1,n2) = -Buv_0(0:n0,n1-1,n2)-Buv_0(0:n0,n1-1,n2+1)
        B_max = max(B_max,abs(B_0(0,n1,n2)))
      end do
    end do

!    write(*,*) 'B_max=',B_max

    B_max=max(B_max,maxval(abs(B_i(0,0:rmaxB,1:2))))

    ! determine inverse Gram matrix
!    mm02 = elimminf2_coli(m02)
!    mm12 = elimminf2_coli(m12)
!    mm22 = elimminf2_coli(m22)
!    q10  = elimminf2_coli(p10)
!    q21  = elimminf2_coli(p21)
!    q20  = elimminf2_coli(p20)
! 
!    q1q2 = (q10+q20-q21)
!    detZ = 4d0*q10*q20-q1q2*q1q2
!    Zinv(1,1) = 2d0*q20/detZ
!    Zinv(2,1) = -q1q2/detZ
!    Zinv(1,2) = Zinv(2,1)
!    Zinv(2,2) = 2d0*q10/detZ
!    f(1) = q10+mm02-mm12
!    f(2) = q20+mm02-mm22

    ! commented out 2.9.2017
    ! Zinv = Zadj/detZ

    ! calculate Cuv
    call CalcCuv(Cuv,Buv_0,mm02,f,rmax,id)

    ! initialization of error propagation
!    Zadj=Zinv*detZ    

!    maxZadj = max(abs(Zadj(1,1)),abs(Zadj(2,1)),abs(Zadj(2,2)))

!    Zadjf(1) = Zadj(1,1)*f(1)+Zadj(2,1)*f(2)
!    Zadjf(2) = Zadj(1,2)*f(1)+Zadj(2,2)*f(2)
!    maxZadjf = max(abs(Zadjf(1)),abs(Zadjf(2)))

!    aZadjff = abs(Zadjf(1)*f(1)+Zadjf(2)*f(2))

!    adetZ = abs(detZ)
!    adetX = abs(2d0*mm02*detZ-Zadjf(1)*f(1)-Zadjf(2)*f(2))

!    write(*,*) 'adZ=',maxZadj,adetZ


    Cij_err =0d0
    C00_err =0d0
    Cij_err(0) = Cerr(0)
    B_err = acc_def_B*B_max

    Cij_err2 =0d0
    C00_err2 =0d0
    Cij_err2(0) = Cerr2(0)

!    write(*,*) 'CalcCpv1o: B_err= ',B_err,acc_def_B,B_max

    allocate(C_alt(0:rmax,0:rmax,0:rmax))

    ! PV reduction
    do r=1,rmax

      ! reduction formula (5.10) for C(r/2,0,0)
      do n0=r/2,1,-1
        do n1=0,r-2*n0
          n2 = r-2*n0-n1
          C(n0,n1,n2) = (B_0(n0-1,n1,n2) + 2*mm02*C(n0-1,n1,n2) + 4*Cuv(n0,n1,n2) &
                        + f(1)*C(n0-1,n1+1,n2) + f(2)*C(n0-1,n1,n2+1)) / (2*r)
        end do
      end do

!     do n0=(r-1)/2,0,-1
      n0=0 
        do n1=0,r-2*n0
          n2 = r-2*n0-n1

          if (n1.ge.1) then
            nn1 = n1-1
            nn2 = n2
            j = 1
          else
            nn1 = n1
            nn2 = n2-1
            j = 2
          end if
         
          ! reduction formula (5.11) for C(n0,n1,n2), n1+n2=/=0
          do i=1,2
            Smod(i) = -B_0(n0,nn1,nn2)-f(i)*C(n0,nn1,nn2)
          end do
          
          if (nn1.ge.1) then
            Smod(1) = Smod(1) - 2d0*nn1*C(n0+1,nn1-1,nn2)
          else          
            Smod(1) = Smod(1) + B_i(n0,nn2,1)
          end if

          if (nn2.ge.1) then
            Smod(2) = Smod(2) - 2d0*nn2*C(n0+1,nn1,nn2-1)
          else
            Smod(2) = Smod(2) + B_i(n0,nn1,2)
          end if
          
          C(n0,n1,n2) = Zinv(1,j)*Smod(1) + Zinv(2,j)*Smod(2)

!        if(n0.eq.0) then
!          write(*,*) 'Ca(0,n1,n2)=',n1,n2,C(0,n1,n2),nn1,nn2,j
!          write(*,*) 'Ca(0,n1,n2)=',Zinv(1,j),Smod(1),Zinv(2,j),Smod(2)
!        end if
 
        end do
!      end do

      ! determine error from symmetry for n0=0 and n1>=1, n2>=1 
      Cerr(r)=Cerr(r-1)
      Cerr2(r)=Cerr2(r-1)
      n0=0
        do n1=1,r-2*n0-1
          n2 = r-2*n0-n1

          nn1 = n1
          nn2 = n2-1
          j = 2
         
          ! reduction formula (5.11) for C(n0,n1,n2), n1+n2=/=0
          do i=1,2
            Smod(i) = -B_0(n0,nn1,nn2)-f(i)*C(n0,nn1,nn2)
          end do
          
          if (nn1.ge.1) then
            Smod(1) = Smod(1) - 2d0*nn1*C(n0+1,nn1-1,nn2)
          else          
            Smod(1) = Smod(1) + B_i(n0,nn2,1)
          end if

          if (nn2.ge.1) then
            Smod(2) = Smod(2) - 2d0*nn2*C(n0+1,nn1,nn2-1)
          else
            Smod(2) = Smod(2) + B_i(n0,nn1,2)
          end if
           
          C_alt(n0,n1,n2) = Zinv(1,j)*Smod(1) + Zinv(2,j)*Smod(2)
 
          Cerr(r)=max(Cerr(r),abs(C(n0,n1,n2)-C_alt(n0,n1,n2)))
          Cerr2(r)=max(Cerr2(r),abs(C(n0,n1,n2)-C_alt(n0,n1,n2)))

!            if(n0.eq.0) then
!              write(*,*) 'Cb(0,n1,n2)=',n1,n2,C_alt(0,n1,n2),nn1,nn2,j
!              write(*,*) 'Cb(0,n1,n2)=',Zinv(1,j),Smod(1),Zinv(2,j),Smod(2)
!            end if
!            write(*,*) 'CalcCpv1o Cerr',n0,n1,n2, Cerr(r), abs(C(n0,n1,n2)),abs(C_alt(n0,n1,n2))

        end do

      if(r.ge.2)then
        C00_err(r) = max(2*abs(m02)*Cij_err(r-2), B_err,    &
              aZadjff/adetZ*Cij_err(r-2),             &
              maxZadjf/adetZ*max(2*C00_err(r-1),B_err))/(2*r)

!        write(*,*) 'C00errtest',r,abs(m02)*Cij_err(r-2), B_err,    &
!              aZadjff/adetZ*Cij_err(r-2),             &
!              maxZadjf/adetZ*C00_err(r-1),maxZadjf/adetZ*B_err, &
!              C00_err(r)

      else
        C00_err(r) = 0d0
      end if
      Cij_err(r) = max(maxZadjf*Cij_err(r-1),   &
              maxZadj*max(2*C00_err(r),B_err))/adetZ
 
      if(r.ge.2)then
        C00_err2(r) = max(2*abs(m02)*Cij_err2(r-2), B_err,    &
              aZadjff/adetZ*Cij_err2(r-2),             &
              maxZadjf/adetZ*max(2*C00_err(r-1),B_err))/(2*r)

!        write(*,*) 'C00errtest',r,abs(m02)*Cij_err2(r-2), B_err,    &
!              aZadjff/adetZ*Cij_err2(r-2),             &
!              maxZadjf/adetZ*C00_err2(r-1),maxZadjf/adetZ*B_err, &
!              C00_err2(r)

      else
        C00_err2(r) = 0d0
      end if
      Cij_err2(r) = max((maxZadjf/maxZadj)*Cij_err2(r-1),max(2*C00_err2(r),B_err))/sqrt(adetZ)
    end do

      ! reduction formula (5.10) for n0+n1+n2=r, n0>0 
    do r=rmax+1,2*rmax
      do n0=r-rmax,r/2
        do n1=0,r-2*n0
          n2 = r-2*n0-n1
          C(n0,n1,n2) = (B_0(n0-1,n1,n2) + 2*mm02*C(n0-1,n1,n2) + 4*Cuv(n0,n1,n2) &
                        + f(1)*C(n0-1,n1+1,n2) + f(2)*C(n0-1,n1,n2+1)) / (2*r)
        end do
      end do
    end do




    Cerr2 = max(Cerr2,Cij_err2(0:rmax))
    Cerr = max(Cerr,Cij_err(0:rmax))



!   write(*,*) 'CalcCpv1o out',(((C((r-n1-n2)/2,n1,n2),n2=0,r-n1),n1=0,r),r=0,rmax)

!    write(*,*) 'CalcCpv1o Cerr ',Cerr
!    write(*,*) 'CalcCpv1o Cerr2',Cerr2

  end subroutine CalcCpv1o   




  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !  subroutine CalcCpv(C,Cuv,p10,p21,p20,m02,m12,m22,rmax,Cerr,Cerr2)
  !
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  
  subroutine CalcCpv(C,Cuv,p10,p21,p20,m02,m12,m22,rmax,id,Cerr,Cerr2)
  
    use globalC
  
    integer, intent(in) :: rmax,id
    double complex, intent(in) :: p10,p21,p20,m02,m12,m22
    double complex, intent(out) :: Cuv(0:rmax,0:rmax,0:rmax)
    double complex, intent(out) :: C(0:rmax,0:rmax,0:rmax)
    double precision, intent(out)  :: Cerr(0:rmax),Cerr2(0:rmax)
    double complex, allocatable :: B_0(:,:,:), Buv_0(:,:,:)
    double complex, allocatable :: B_i(:,:,:), Buv_i(:,:,:)
    double complex, allocatable :: C_alt(:,:,:)
    double complex :: Smod(2)
    double complex :: C0_coli, elimminf2_coli
    double precision, allocatable :: C00_err(:),Cij_err(:)
    double precision, allocatable :: C00_err2(:),Cij_err2(:)
    double precision :: B_err,B_max
    integer :: rmaxB,r,n0,n1,n2,nn0,nn1,nn2,i,j
    integer :: bin,k,nid(0:2)



    ! calculation of scalar coefficient
    C(0,0,0) = C0_coli(p10,p21,p20,m02,m12,m22)
    Cuv(0,0,0) = 0d0

    ! accuracy estimate for C0 function
    Cerr(0) = acc_def_C0*max(1d0/sqrt(adetZ),abs(C(0,0,0)))   
    Cerr2(0) = acc_def_C0*max(1d0/sqrt(adetZ),abs(C(0,0,0)))   

!    write(*,*) 'CalcCpv: Cerr(0)= ',Cerr(0),Cerr(0)/abs(C(0,0,0)),abs(C(0,0,0))

    if (rmax.eq.0) return
    
    ! allocation and calculation of B functions
    rmaxB = rmax-1
    ! rmaxB = max(rmax-1,0)
    allocate(B_0(0:rmaxB,0:rmaxB,0:rmaxB))
    allocate(Buv_0(0:rmaxB,0:rmaxB,0:rmaxB))
    allocate(B_i(0:rmaxB,0:rmaxB,2))
    allocate(Buv_i(0:rmaxB,0:rmaxB,2))

    ! allocate arrays for error propagation
    allocate(C00_err(0:rmax))
    allocate(Cij_err(0:rmax))
    allocate(C00_err2(0:rmax))
    allocate(Cij_err2(0:rmax))
    
    ! determine binaries for B-coefficients
    k=0
    bin = 1
    do while (k.le.2)
      if (mod(id/bin,2).eq.0) then
        nid(k) = id+bin
        k = k+1
      end if
      bin = 2*bin
    end do

    call CalcB(B_0(:,0,:),Buv_0(:,0,:),p21,m12,m22,rmaxB,nid(0))
    call CalcB(B_i(:,:,1),Buv_i(:,:,1),p20,m02,m22,rmaxB,nid(1))
    call CalcB(B_i(:,:,2),Buv_i(:,:,2),p10,m02,m12,rmaxB,nid(2))

    ! shift of integration momentum in B_0 and calculate maximal B(0,...)
    B_max=0d0
    do n1=1,rmaxB
      do n2=0,rmaxB-n1
        n0 = (rmaxB-n1-n2)
        B_0(0:n0,n1,n2) = -B_0(0:n0,n1-1,n2)-B_0(0:n0,n1-1,n2+1)
        Buv_0(0:n0,n1,n2) = -Buv_0(0:n0,n1-1,n2)-Buv_0(0:n0,n1-1,n2+1)
        B_max = max(B_max,abs(B_0(0,n1,n2)))
      end do
    end do

!    write(*,*) 'B_max=',B_max

    B_max=max(B_max,maxval(abs(B_i(0,0:rmaxB,1:2))))

    ! determine inverse Gram matrix
!    mm02 = elimminf2_coli(m02)
!    mm12 = elimminf2_coli(m12)
!    mm22 = elimminf2_coli(m22)
!    q10  = elimminf2_coli(p10)
!    q21  = elimminf2_coli(p21)
!    q20  = elimminf2_coli(p20)
! 
!    q1q2 = (q10+q20-q21)
!    detZ = 4d0*q10*q20-q1q2*q1q2
!    Zinv(1,1) = 2d0*q20/detZ
!    Zinv(2,1) = -q1q2/detZ
!    Zinv(1,2) = Zinv(2,1)
!    Zinv(2,2) = 2d0*q10/detZ
!    f(1) = q10+mm02-mm12
!    f(2) = q20+mm02-mm22

    ! commented out 2.9.2017
    ! Zinv = Zadj/detZ

    ! calculate Cuv
    call CalcCuv(Cuv,Buv_0,mm02,f,rmax,id)

    ! initialization of error propagation
!    Zadj=Zinv*detZ    

!    maxZadj = max(abs(Zadj(1,1)),abs(Zadj(2,1)),abs(Zadj(2,2)))

!    Zadjf(1) = Zadj(1,1)*f(1)+Zadj(2,1)*f(2)
!    Zadjf(2) = Zadj(1,2)*f(1)+Zadj(2,2)*f(2)
!    maxZadjf = max(abs(Zadjf(1)),abs(Zadjf(2)))

!    aZadjff = abs(Zadjf(1)*f(1)+Zadjf(2)*f(2))

!    adetZ = abs(detZ)
!    adetX = abs(2d0*mm02*detZ-Zadjf(1)*f(1)-Zadjf(2)*f(2))

!    write(*,*) 'adZ=',maxZadj,adetZ


    Cij_err =0d0
    C00_err =0d0
    Cij_err(0) = Cerr(0)
    B_err = acc_def_B*B_max

    Cij_err2 =0d0
    C00_err2 =0d0
    Cij_err2(0) = Cerr2(0)

!    write(*,*) 'CalcCpv: B_err= ',B_err,acc_B,B_max

    allocate(C_alt(0:rmax,0:rmax,0:rmax))

    ! PV reduction
    do r=1,rmax

      if (mod(r,2).eq.0) then
        ! reduction formula (5.10) for C(r/2,0,0)
        n0 = r/2        
        C(n0,0,0) = (B_0(n0-1,0,0) + 2*mm02*C(n0-1,0,0) + 4*Cuv(n0,0,0) &
                        + f(1)*C(n0-1,1,0) + f(2)*C(n0-1,0,1)) / (2*r)
      end if

      do n0=(r-1)/2,0,-1
        do n1=0,r-2*n0
          n2 = r-2*n0-n1

          if (n1.ge.1) then
            nn1 = n1-1
            nn2 = n2
            j = 1
          else
            nn1 = n1
            nn2 = n2-1
            j = 2
          end if
         
          ! reduction formula (5.11) for C(n0,n1,n2), n1+n2=/=0
          do i=1,2
            Smod(i) = -B_0(n0,nn1,nn2)-f(i)*C(n0,nn1,nn2)
          end do
          
          if (nn1.ge.1) then
            Smod(1) = Smod(1) - 2d0*nn1*C(n0+1,nn1-1,nn2)
          else          
            Smod(1) = Smod(1) + B_i(n0,nn2,1)
          end if

          if (nn2.ge.1) then
            Smod(2) = Smod(2) - 2d0*nn2*C(n0+1,nn1,nn2-1)
          else
            Smod(2) = Smod(2) + B_i(n0,nn1,2)
          end if
          
          C(n0,n1,n2) = Zinv(1,j)*Smod(1) + Zinv(2,j)*Smod(2)

!        if(n0.eq.0) then
!          write(*,*) 'Ca(0,n1,n2)=',n1,n2,C(0,n1,n2),nn1,nn2
!          write(*,*) 'Ca(0,n1,n2)=',Zinv(1,j),Smod(1),Zinv(2,j),Smod(2)
!        end if
 
        end do
      end do

      ! determine error from symmetry for n0=0 and n1>1, n2>1 
      Cerr(r)=Cerr(r-1)
      Cerr2(r)=Cerr2(r-1)
      n0=0
        do n1=0,r-2*n0
          n2 = r-2*n0-n1

          if (n1.ge.1.and.n2.ge.1) then
            nn1 = n1
            nn2 = n2-1
            j = 2
         
            ! reduction formula (5.11) for C(n0,n1,n2), n1+n2=/=0
            do i=1,2
              Smod(i) = -B_0(n0,nn1,nn2)-f(i)*C(n0,nn1,nn2)
            end do
          
            if (nn1.ge.1) then
              Smod(1) = Smod(1) - 2d0*nn1*C(n0+1,nn1-1,nn2)
            else          
              Smod(1) = Smod(1) + B_i(n0,nn2,1)
            end if

            if (nn2.ge.1) then
              Smod(2) = Smod(2) - 2d0*nn2*C(n0+1,nn1,nn2-1)
            else
              Smod(2) = Smod(2) + B_i(n0,nn1,2)
            end if
           
            C_alt(n0,n1,n2) = Zinv(1,j)*Smod(1) + Zinv(2,j)*Smod(2)
 
            Cerr(r)=max(Cerr(r),abs(C(n0,n1,n2)-C_alt(n0,n1,n2)))
            Cerr2(r)=max(Cerr2(r),abs(C(n0,n1,n2)-C_alt(n0,n1,n2)))

          end if
        end do

      if(r.ge.2)then
        C00_err(r) = max(abs(m02)*Cij_err(r-2), B_err,    &
              aZadjff/adetZ*Cij_err(r-2),             &
              maxZadjf/adetZ*max(C00_err(r-1),B_err))

!        write(*,*) 'C00errtest',r,abs(m02)*Cij_err(r-2), B_err,    &
!              aZadjff/adetZ*Cij_err(r-2),             &
!              maxZadjf/adetZ*C00_err(r-1),maxZadjf/adetZ*B_err, &
!              C00_err(r)

      else
        C00_err(r) = 0d0
      end if
      Cij_err(r) = max(maxZadjf*Cij_err(r-1),   &
              maxZadj*max(C00_err(r),B_err))/adetZ
 
      if(r.ge.2)then
        C00_err2(r) = max(abs(m02)*Cij_err(r-2), B_err,    &
              aZadjff/adetZ*Cij_err2(r-2),             &
              maxZadjf/adetZ*max(C00_err2(r-1),B_err))

!        write(*,*) 'C00errtest',r,abs(m02)*Cij_err(r-2), B_err,    &
!              aZadjff/adetZ*Cij_err(r-2),             &
!              maxZadjf/adetZ*C00_err(r-1),maxZadjf/adetZ*B_err, &
!              C00_err(r)

      else
        C00_err2(r) = 0d0
      end if
      Cij_err2(r) = max((maxZadjf/maxZadj)*Cij_err2(r-1),   &
                        max(C00_err2(r),B_err))/sqrt(adetZ)
 
!      write(*,*) 'CalcCpv r',r, Cij_err(r),maxZadjf*Cij_err(r-1)/adetZ,  &
!              maxZadj*(C00_err(r))/adetZ,  &          
!              maxZadj*(B_err)/adetZ

    end do


      ! reduction formula (5.10) for n0+n1+n2=r, n0>0 
    ! PV reduction (5.10)
    do r=rmax+1,2*rmax
      do n0=r-rmax,r/2
        do n1=0,r-2*n0
          n2 = r-2*n0-n1
          C(n0,n1,n2) = (B_0(n0-1,n1,n2) + 2*mm02*C(n0-1,n1,n2) + 4*Cuv(n0,n1,n2) &
                        + f(1)*C(n0-1,n1+1,n2) + f(2)*C(n0-1,n1,n2+1)) / (2*r)
        end do
      end do
    end do



    Cerr2 = max(Cerr2,Cij_err2(0:rmax))
    Cerr = max(Cerr,Cij_err(0:rmax))



!   write(*,*) 'CalcCpv out',(((C((r-n1-n2)/2,n1,n2),n2=0,r-n1),n1=0,r),r=0,rmax)

  end subroutine CalcCpv





  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !  subroutine CalcCpv2(C,Cuv,p10,p21,p20,m02,m12,m22,rmax,id,Cerr,Cerr2)
  !
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  
  subroutine CalcCpv2(C,Cuv,p10,p21,p20,m02,m12,m22,rmax,id,Cerr,Cerr2)
  
    use globalC
  
    integer, intent(in) :: rmax,id
    double complex, intent(in) :: p10,p21,p20,m02,m12,m22
    double complex, intent(out) :: C(0:rmax,0:rmax,0:rmax)
    double complex, intent(out) :: Cuv(0:rmax,0:rmax,0:rmax)
    double precision, intent(out) :: Cerr(0:rmax),Cerr2(0:rmax)
    double complex, allocatable :: B_0(:,:,:), B_i(:,:,:)
    double complex, allocatable :: Buv_0(:,:,:), Buv_i(:,:,:)
    double complex, allocatable :: C_alt(:,:,:)
    double complex :: C0_coli, elimminf2_coli
!   double complex :: Caux(1:rmax/2+1,0:rmax-1,0:rmax-1), Smod(2)
    double complex :: Caux(1:rmax,0:rmax-1,0:rmax-1), Smod(2)
    double complex :: chdet
    double precision, allocatable :: C00_err(:),Cij_err(:)
    double precision, allocatable :: C00_err2(:),Cij_err2(:)
    double precision :: B_err,B_max
    integer :: rmaxB,r,n0,n1,n2,k
    integer :: bin,nid(0:3)



    ! calculation of scalar coefficient
    C(0,0,0) = C0_coli(p10,p21,p20,m02,m12,m22)
    Cuv(0,0,0) = 0d0

    ! accuracy estimate for C0 function
    Cerr(0) = acc_def_C0*max( abs(C(0,0,0)), 1d0/sqrt(adetZ) )
    Cerr2(0) = acc_def_C0*max( abs(C(0,0,0)), 1d0/sqrt(adetZ) )

!   write(*,*) 'CalcCpv2: Cerr(0)= ',Cerr(0),Cerr(0)/abs(C(0,0,0)),abs(C(0,0,0))

    if (rmax.eq.0) return


    ! calculation of B-coefficients
    rmaxB = rmax-1
    ! rmaxB = max(rmax-1,0)
    allocate(B_0(0:rmaxB,0:rmaxB,0:rmaxB))
    allocate(Buv_0(0:rmaxB,0:rmaxB,0:rmaxB))
    allocate(B_i(0:rmaxB,0:rmaxB,2))
    allocate(Buv_i(0:rmaxB,0:rmaxB,2))

    ! allocate arrays for error propagation
    allocate(C00_err(0:rmax+1))
    allocate(Cij_err(0:rmax))
    allocate(C00_err2(0:rmax+1))
    allocate(Cij_err2(0:rmax))
    
    ! determine binaries for B-coefficients
    k=0
    bin = 1
    do while (k.le.2)
      if (mod(id/bin,2).eq.0) then
        nid(k) = id+bin
        k = k+1
      end if
      bin = 2*bin
    end do

    call CalcB(B_0(:,0,:),Buv_0(:,0,:),p21,m12,m22,rmaxB,nid(0))
    call CalcB(B_i(:,:,1),Buv_i(:,:,1),p20,m02,m22,rmaxB,nid(1))
    call CalcB(B_i(:,:,2),Buv_i(:,:,2),p10,m02,m12,rmaxB,nid(2))

    ! shift of integration momentum in B_0 and calculate maximal B(0,...)
    B_max=0d0
    do n1=1,rmaxB
      do n2=0,rmaxB-n1
        n0 = (rmaxB-n1-n2)
        B_0(0:n0,n1,n2) = -B_0(0:n0,n1-1,n2)-B_0(0:n0,n1-1,n2+1)
        Buv_0(0:n0,n1,n2) = -Buv_0(0:n0,n1-1,n2)-Buv_0(0:n0,n1-1,n2+1)
        B_max = max(B_max,abs(B_0(0,n1,n2)))
      end do
    end do
    B_max=max(B_max,maxval(abs(B_i(0,0:rmaxB,1:2))))
    
    ! determine inverse modified Cayley matrix
!    mm02 = elimminf2_coli(m02)
!    mm12 = elimminf2_coli(m12)
!    mm22 = elimminf2_coli(m22)
!    q10  = elimminf2_coli(p10)
!    q21  = elimminf2_coli(p21)
!    q20  = elimminf2_coli(p20)


    ! calculate Cuv
    call CalcCuv(Cuv,Buv_0,mm02,mx(1:2,0),rmax,id)

    ! initialization of error propagation

!    adetX = abs(chdet(3,mx))
!    maxZadjf=maxval(abs(mxinv(0,1:2)))*adetX
!    maxXadj=maxval(abs(mxinv(1:2,1:2)))*adetX
!    adetZ=abs(mxinv(0,0))*adetX

!   write(*,*) 'CalcCpv2 adetX ',adetX,maxZadjf,maxXadj,adetZ

    Cij_err =0d0
    C00_err =0d0
    Cij_err(0) = Cerr(0)
    B_err = acc_def_B*B_max

    Cij_err2 =0d0
    C00_err2 =0d0
    Cij_err2(0) = Cerr2(0)

!    write(*,*) 'CalcCpv: B_err= ',B_err,acc_B,B_max

    allocate(C_alt(0:rmax,0:rmax,0:rmax))
    
    ! alternative PV-like reduction
    do r=1,rmax

    ! calculate C n0>1 using (5.14)
      do n0=2,r/2
        do n1=0,r-2*n0
          n2 = r-2*n0-n1

          do k=1,2
            Smod(k) = -B_0(n0-1,n1,n2)
          end do

          if (n1.ge.1) then
            Smod(1) = Smod(1) - 2*n1*C(n0,n1-1,n2)
          else
            Smod(1) = Smod(1) + B_i(n0-1,n2,1)
          end if

          if (n2.ge.1) then
            Smod(2) = Smod(2) - 2*n2*C(n0,n1,n2-1)
          else
            Smod(2) = Smod(2) + B_i(n0-1,n1,2)
          end if

          Caux(n0,n1,n2) = (C(n0-1,n1,n2) - mxinv(1,0)*Smod(1)  &
                             - mxinv(2,0)*Smod(2))/mxinv(0,0)

        end do
      end do

      
      do n0=1,r/2
        do n1=0,r-2*n0
          n2 = r-2*n0-n1

          C(n0,n1,n2) = (Caux(n0,n1,n2) + 4d0*Cuv(n0,n1,n2)  &
                         + B_0(n0-1,n1,n2))/r/2d0    

        end do
      end do

      
      ! calculate C for n0=0 and n1>0, n2>0 from (5.15) with (5.14) inserted 
      ! and determine error from symmetry  
      Cerr(r)=Cerr(r-1)
      Cerr2(r)=Cerr2(r-1)

      do n1=0,r-1
        n2 = r-1-n1

        do k=1,2
          Smod(k) = -B_0(0,n1,n2)
        end do

        if (n1.ge.1) then
          Smod(1) = Smod(1) - 2*n1*C(1,n1-1,n2)
        else
          Smod(1) = Smod(1) + B_i(0,n2,1)
        end if

        if (n2.ge.1) then
          Smod(2) = Smod(2) - 2*n2*C(1,n1,n2-1)
        else
          Smod(2) = Smod(2) + B_i(0,n1,2)
        end if

        Caux(1,n1,n2) = (C(0,n1,n2) - mxinv(1,0)*Smod(1)  &
                           - mxinv(2,0)*Smod(2))/mxinv(0,0)

        C(0,n1+1,n2) = mxinv(0,1)*Caux(1,n1,n2)  &
                       + mxinv(1,1)*Smod(1) + mxinv(2,1)*Smod(2)
        C_alt(0,n1,n2+1) = mxinv(0,2)*Caux(1,n1,n2)  &
                       + mxinv(1,2)*Smod(1) + mxinv(2,2)*Smod(2)

        if(n1.eq.0) then
          C(0,0,r) = C_alt(0,0,r)          
        else
          Cerr(r)=max(Cerr(r),abs(C(0,n1,n2+1)-C_alt(0,n1,n2+1)))
          Cerr2(r)=max(Cerr2(r),abs(C(0,n1,n2+1)-C_alt(0,n1,n2+1)))
        end if

      end do 

      C00_err(r+1) = max(B_err,adetX/adetZ*Cij_err(r-1),      &
              maxZadjf/adetZ*max(B_err,C00_err(r)))/(2*(r+1))

!      write(*,*) 'CalcCpv2 00 r',r, B_err,adetX/adetZ*Cij_err(r-1),      &
!              maxZadjf/adetZ*B_err, maxZadjf/adetZ*C00_err(r)

      Cij_err(r) = max(maxZadjf*max(2*(r+1)*C00_err(r+1),B_err),   &
              maxXadj*max(2*C00_err(r),B_err))/adetX

      C00_err2(r+1) = max(B_err,adetX/adetZ*Cij_err2(r-1),      &
              maxZadjf/adetZ*max(B_err,C00_err2(r)))/(2*(r+1))

!      write(*,*) 'CalcCpv2 00 r',r, B_err,adetX/adetZ*Cij_err(r-1),      &
!              maxZadjf/adetZ*B_err, maxZadjf/adetZ*C00_err(r)

      Cij_err2(r) = max(maxZadjf*max(2*(r+1)*C00_err2(r+1),B_err),   &
              maxXadj*max(2*C00_err2(r),B_err))/adetX*(sqrt(adetZ)/maxZadj)

!      write(*,*) 'CalcCpv2 ij r',r, maxZadjf*C00_err(r+1)/adetX,B_err*maxZadjf/adetX, &
!               maxXadj*C00_err(r)/adetX, maxXadj*B_err/adetX
      
    end do
    

    ! reduction formula (5.10) for n0+n1+n2=r, n0>0 
    do r=rmax+1,2*rmax



! pv2 formulas added 24.01.2016
      do n0=max(2,r-rmax),r/2
        do n1=0,r-2*n0
          n2 = r-2*n0-n1

          do k=1,2
            Smod(k) = -B_0(n0-1,n1,n2)
          end do

          if (n1.ge.1) then
            Smod(1) = Smod(1) - 2*n1*C(n0,n1-1,n2)
          else
            Smod(1) = Smod(1) + B_i(n0-1,n2,1)
          end if

          if (n2.ge.1) then
            Smod(2) = Smod(2) - 2*n2*C(n0,n1,n2-1)
          else
            Smod(2) = Smod(2) + B_i(n0-1,n1,2)
          end if

          Caux(n0,n1,n2) = (C(n0-1,n1,n2) - mxinv(1,0)*Smod(1)  &
                             - mxinv(2,0)*Smod(2))/mxinv(0,0)

        end do
      end do

      
      do n0=r-rmax,r/2
        do n1=0,r-2*n0
          n2 = r-2*n0-n1

          C(n0,n1,n2) = (Caux(n0,n1,n2) + 4d0*Cuv(n0,n1,n2)  &
                         + B_0(n0-1,n1,n2))/r/2d0    



        end do
      end do

    end do


    
    Cerr2 = max(Cerr2,Cij_err2(0:rmax))
    Cerr = max(Cerr,Cij_err(0:rmax))



!   write(*,*) 'CalcCpv2 out',(((C((r-n1-n2)/2,n1,n2),n2=0,r-n1),n1=0,r),r=0,rmax)

!    write(*,*) 'CalcCpv2 Cerr ',Cerr
!    write(*,*) 'CalcCpv2 Cerr2',Cerr2

  end subroutine CalcCpv2



  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !  subroutine CalcCpvshift(Cshift,Cuvshift,p10shift,p21shift,p20shift,m02shift,m12shift,m22shift,rmax,Cerr,Cerr2)
  !
  !  Based on CalcCpv1
  !  uses shifted momenta and global shifted quantities
  !
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!7!!!!!!!!!!!!!!!!!!!!!!!!
  
  subroutine CalcCpvshift(Cshift,Cuvshift,p10shift,p21shift,p20shift,m02shift,m12shift,m22shift,rmax,id,Cerr,Cerr2)
  
    use globalC
  
    integer, intent(in) :: rmax,id
    double complex, intent(in) :: p10shift,p21shift,p20shift,m02shift,m12shift,m22shift
    double complex, intent(out) :: Cuvshift(0:rmax,0:rmax,0:rmax)
    double complex, intent(out) :: Cshift(0:rmax,0:rmax,0:rmax)
    double precision, intent(out) :: Cerr(0:rmax),Cerr2(0:rmax)
    double complex, allocatable :: B_0(:,:,:), Buv_0(:,:,:)
    double complex, allocatable :: B_i(:,:,:), Buv_i(:,:,:)
    double complex, allocatable :: Cshift_alt(:,:,:)
    double complex :: Smod(2)
    double complex :: C0_coli, elimminf2_coli
    double precision, allocatable :: C00_err(:),Cij_err(:)
    double precision, allocatable :: C00_err2(:),Cij_err2(:)
    double precision :: B_err,B_max
    integer :: rmaxB,r,n0,n1,n2,nn0,nn1,nn2,i,j
    integer :: bin,k,nid(0:2)
    logical :: use_cache_system_save



    ! calculation of scalar coefficient
    Cshift(0,0,0) = C0_coli(p10shift,p21shift,p20shift,m02shift,m12shift,m22shift)
    Cuvshift(0,0,0) = 0d0

    ! accuracy estimate for C0 function
    Cerr(0) = acc_def_C0*max(1d0/sqrt(adetZshift),abs(Cshift(0,0,0)))   
    Cerr2(0) = acc_def_C0*max(1d0/sqrt(adetZshift),abs(Cshift(0,0,0)))   

    if (rmax.eq.0) return
    
    ! allocation and calculation of B functions
    rmaxB = rmax-1
    ! rmaxB = max(rmax-1,0)
    allocate(B_0(0:rmaxB,0:rmaxB,0:rmaxB))
    allocate(Buv_0(0:rmaxB,0:rmaxB,0:rmaxB))
    allocate(B_i(0:rmaxB,0:rmaxB,2))
    allocate(Buv_i(0:rmaxB,0:rmaxB,2))

    ! allocate arrays for error propagation
    allocate(C00_err(0:rmax))
    allocate(Cij_err(0:rmax))
    allocate(C00_err2(0:rmax))
    allocate(Cij_err2(0:rmax))
    
    ! determine binaries for B-coefficients
    k=0
    bin = 1
    do while (k.le.2)
      if (mod(id/bin,2).eq.0) then
        nid(k) = id+bin
        k = k+1
      end if
      bin = 2*bin
    end do

!   call CalcB(B_0(:,0,:),Buv_0(:,0,:),p21shift,m12shift,m22shift,rmaxB,nid(0))
!   call CalcB(B_i(:,:,1),Buv_i(:,:,1),p20shift,m02shift,m22shift,rmaxB,nid(1))
!   call CalcB(B_i(:,:,2),Buv_i(:,:,2),p10shift,m02shift,m12shift,rmaxB,nid(2))
    use_cache_system_save = use_cache_system
    use_cache_system = .false.
    call CalcB(B_0(:,0,:),Buv_0(:,0,:),p21shift,m12shift,m22shift,rmaxB,0)
    call CalcB(B_i(:,:,1),Buv_i(:,:,1),p20shift,m02shift,m22shift,rmaxB,0)
    call CalcB(B_i(:,:,2),Buv_i(:,:,2),p10shift,m02shift,m12shift,rmaxB,0)
    use_cache_system = use_cache_system_save
!   call SwitchOnCacheSystem_cll


    ! shift of integration momentum in B_0 and calculate maximal B(0,...)
    B_max=0d0
    do n1=1,rmaxB
      do n2=0,rmaxB-n1
        n0 = (rmaxB-n1-n2)
        B_0(0:n0,n1,n2) = -B_0(0:n0,n1-1,n2)-B_0(0:n0,n1-1,n2+1)
        Buv_0(0:n0,n1,n2) = -Buv_0(0:n0,n1-1,n2)-Buv_0(0:n0,n1-1,n2+1)
        B_max = max(B_max,abs(B_0(0,n1,n2)))
      end do
    end do

!    write(*,*) 'B_max=',B_max

    B_max=max(B_max,maxval(abs(B_i(0,0:rmaxB,1:2))))

    ! determine inverse Gram matrix
!    Zinvshift = Zadjshift/detZshift

    ! calculate Cuv
    call CalcCuv(Cuvshift,Buv_0,mm02shift,fshift,rmax,id)

    ! initialization of error propagation
    Cij_err =0d0
    C00_err =0d0
    Cij_err(0) = Cerr(0)
    B_err = acc_def_B*B_max

    Cij_err2 =0d0
    C00_err2 =0d0
    Cij_err2(0) = Cerr2(0)




    allocate(Cshift_alt(0:rmax,0:rmax,0:rmax))

    ! PV reduction
    do r=1,rmax

    ! reduction formula (5.10) with (5.11) inserted for n0 >= 1
      do n0=r/2,1,-1
        do n1=0,r-2*n0
          n2 = r-2*n0-n1
          Cshift(n0,n1,n2) = + 4*Cuvshift(n0,n1,n2) + detXshift/detZshift*Cshift(n0-1,n1,n2)
            Cshift(n0,n1,n2) = Cshift(n0,n1,n2) &
                + (detZmZadjfshift + Zadjsshift(1)*(mm12shift-mm02shift) + Zadjsshift(2)*(mm22shift-mm02shift) &
                  ) /detZshift * B_0(n0-1,n1,n2)
!                + (1d0 - (Zadjfshift(1)+Zadjfshift(2))/detZshift)* B_0(n0-1,n1,n2)

            if (n1.ge.1) then
              Cshift(n0,n1,n2) = Cshift(n0,n1,n2) &
                  - 2*n1*Zadjfshift(1)/detZshift*Cshift(n0,n1-1,n2)
            else            
              Cshift(n0,n1,n2) = Cshift(n0,n1,n2) &
                  + Zadjfshift(1)/detZshift* B_i(n0-1,n2,1)
            end if
            if (n2.ge.1) then
              Cshift(n0,n1,n2) = Cshift(n0,n1,n2) &
                  - 2*n2*Zadjfshift(2)/detZshift*Cshift(n0,n1,n2-1)
            else            
              Cshift(n0,n1,n2) = Cshift(n0,n1,n2) &
                  + Zadjfshift(2)/detZshift * B_i(n0-1,n1,2)
            end if

            Cshift(n0,n1,n2) = Cshift(n0,n1,n2)  / (2*r) 

!        if(n0.eq.1) then
!          write(*,*) 'Cas(1,n1,n2)=',n1,n2, 4*Cuvshift(n0,n1,n2) + detXshift/detZshift*Cshift(n0-1,n1,n2)
!          write(*,*) 'Cas(1,n1,n2)=', (detZmZadjfshift + Zadjsshift(1)*(mm12shift-mm02shift) + Zadjsshift(2)*(mm22shift-mm02shift) &
!                  ) /detZshift * B_0(n0-1,n1,n2)
!          write(*,*) 'Cas(1,n1,n2)=', detZmZadjfshift,Zadjsshift(1)*(mm12shift-mm02shift),Zadjsshift(2)*(mm22shift-mm02shift) &
!                  ,detZshift ,B_0(n0-1,n1,n2)
!          write(*,*) 'Cas(1,n1,n2)=', (1d0 - (Zadjfshift(1)+Zadjfshift(2))/detZshift)* B_0(n0-1,n1,n2)
!          write(*,*) 'Cas(1,n1,n2)=',  + Zadjfshift(1)/detZshift* B_i(n0-1,n2,1)
!          write(*,*) 'Cas(1,n1,n2)=',   + Zadjfshift(2)/detZshift * B_i(n0-1,n1,2)
!        end if

         end do
      end do

    ! reduction formula (5.11) with (5.10) inserted for n0 = 0
!     do n0=(r-1)/2,0,-1
      n0=0 
        do n1=0,r-2*n0
          n2 = r-2*n0-n1

          if (n1.ge.1) then
            nn1 = n1-1
            nn2 = n2
            j = 1
          else
            nn1 = n1
            nn2 = n2-1
            j = 2
          end if
         
!          do i=1,2
!            Smod(i) = -B_0(n0,nn1,nn2)
!          end do
          Smod = 0d0
          
          if (nn1.ge.1) then
            Smod(1) = Smod(1) - 2d0*nn1*Cshift(n0+1,nn1-1,nn2)
          else          
            Smod(1) = Smod(1) + B_i(n0,nn2,1)
          end if

          if (nn2.ge.1) then
            Smod(2) = Smod(2) - 2d0*nn2*Cshift(n0+1,nn1,nn2-1)
          else
            Smod(2) = Smod(2) + B_i(n0,nn1,2)
          end if
          
          Cshift(n0,n1,n2) = (Zadjshift(1,j)*Smod(1) + Zadjshift(2,j)*Smod(2)  &
                           - Zadjsshift(j)*B_0(n0,nn1,nn2) &
                           - Zadjfshift(j)*Cshift(n0,nn1,nn2))/detZshift

!          if(n0.eq.0) then
!            write(*,*) 'Cas(0,n1,n2)=',n1,n2,Cshift(0,n1,n2),nn1,nn2,j
!            write(*,*) 'Cas(0,n1,n2)=',Zadjshift(1,j),Smod(1),Zadjshift(2,j),Smod(2)
!            write(*,*) 'Cas(0,n1,n2)=',Zadjsshift(j),B_0(n0,nn1,nn2),Zadjfshift(j),Cshift(n0,nn1,nn2)
!            write(*,*) 'Cas(0,n1,n2)=',Zadjshift(1,j)*Smod(1),Zadjshift(2,j)*Smod(2)
!            write(*,*) 'Cas(0,n1,n2)=',-Zadjsshift(j)*B_0(n0,nn1,nn2),-Zadjfshift(j)*Cshift(n0,nn1,nn2)
!          end if
          
        end do
!      end do

      ! determine error from symmetry for n0=0 and n1>=1, n2>=1 
      Cerr(r)=Cerr(r-1)
      Cerr2(r)=Cerr2(r-1)
      n0=0
        do n1=1,r-2*n0-1
          n2 = r-2*n0-n1

          nn1 = n1
          nn2 = n2-1
          j = 2
         
!          do i=1,2
!            Smod(i) = -B_0(n0,nn1,nn2)
!          end do
          Smod = 0
          
          if (nn1.ge.1) then
            Smod(1) = Smod(1) - 2d0*nn1*Cshift(n0+1,nn1-1,nn2)
          else          
            Smod(1) = Smod(1) + B_i(n0,nn2,1)
          end if

          if (nn2.ge.1) then
            Smod(2) = Smod(2) - 2d0*nn2*Cshift(n0+1,nn1,nn2-1)
          else
            Smod(2) = Smod(2) + B_i(n0,nn1,2)
          end if
           
          Cshift_alt(n0,n1,n2) = (Zadjshift(1,j)*Smod(1) + Zadjshift(2,j)*Smod(2)  &
                           - Zadjsshift(j)*B_0(n0,nn1,nn2) &
                           - Zadjfshift(j)*Cshift(n0,nn1,nn2))/detZshift
 
          Cerr(r)=max(Cerr(r),abs(Cshift(n0,n1,n2)-Cshift_alt(n0,n1,n2)))
          Cerr2(r)=max(Cerr2(r),abs(Cshift(n0,n1,n2)-Cshift_alt(n0,n1,n2)))

!          if(n0.eq.0) then
!            write(*,*) 'Cbs(0,n1,n2)=',n1,n2,Cshift_alt(0,n1,n2),nn1,nn2,j
!            write(*,*) 'Cbs(0,n1,n2)=',Zadjshift(1,j),Smod(1),Zadjshift(2,j),Smod(2)
!            write(*,*) 'Cbs(0,n1,n2)=',Zadjsshift(j),B_0(n0,nn1,nn2),Zadjfshift(j),Cshift(n0,nn1,nn2)
!          end if
!            write(*,*) 'CalcCpvshift Cerr',n0,n1,n2, Cerr(r), abs(Cshift(n0,n1,n2)),abs(Cshift_alt(n0,n1,n2))

        end do

      if(r.ge.2)then
! estimate using insertions of (5.11) in (5.10)
        C00_err(r) = max(2*abs(m02shift)*Cij_err(r-2), B_err,    &
              aZadjffshift/adetZshift*Cij_err(r-2),             &
              maxZadjfshift/adetZshift*max(2*C00_err(r-1),B_err))/(2*r)

      else
        C00_err(r) = 0d0
      end if
! estimate using insertions of (5.10) in (5.11)
      Cij_err(r) = max(maxZadjfshift*Cij_err(r-1),   &
              maxZadjshift*max(2*C00_err(r),B_err))/adetZshift
 
      if(r.ge.2)then
        C00_err2(r) = max(2*abs(m02shift)*Cij_err2(r-2), B_err,    &
              aZadjffshift/adetZshift*Cij_err2(r-2),             &
              maxZadjfshift/adetZshift*max(2*C00_err(r-1),B_err))/(2*r)

      else
        C00_err2(r) = 0d0
      end if
      Cij_err2(r) = max((maxZadjfshift/maxZadjshift)*Cij_err2(r-1),max(2*C00_err2(r),B_err))/sqrt(adetZshift)
    end do

      ! reduction formula (5.10) for n0+n1+n2=r, n0>0 
    do r=rmax+1,2*rmax
      do n0=r-rmax,r/2
        do n1=0,r-2*n0
          n2 = r-2*n0-n1
          Cshift(n0,n1,n2) = (B_0(n0-1,n1,n2) + 2*mm02shift*Cshift(n0-1,n1,n2) + 4*Cuvshift(n0,n1,n2) &
                        + fshift(1)*Cshift(n0-1,n1+1,n2) + fshift(2)*Cshift(n0-1,n1,n2+1)) / (2*r)
        end do
      end do
    end do




    Cerr2 = max(Cerr2,Cij_err2(0:rmax))
    Cerr = max(Cerr,Cij_err(0:rmax))







  end subroutine CalcCpvshift





  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !
  !   Version derived from CalcDg     AD 28.11.2014
  !
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !  subroutine CalcCgn(C,Cuv,p10,p21,p20,m02,m12,m22,rmax,ordg_min,ordg_max,id,Cerr,acc_req_Cr,Cerr2)
  !
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  
  subroutine CalcCgn(C,Cuv,p10,p21,p20,m02,m12,m22,rmax,ordg_min,ordg_max,id,Cerr,acc_req_Cr,Cerr2)
  
    use globalC

    integer, intent(in) :: rmax,ordg_min,ordg_max,id
    double complex, intent(in) ::  p10,p21,p20,m02,m12,m22
    double complex, intent(out) :: C(0:rmax,0:rmax,0:rmax)
    double complex, intent(out) :: Cuv(0:rmax,0:rmax,0:rmax)
    double precision, intent(out) :: Cerr(0:rmax),Cerr2(0:rmax)
    double precision, intent(in) :: acc_req_Cr(0:rmax)
    double complex :: Xtilde,Zkl,Zadjfj,Zadj2,Zadjkl
    double complex, allocatable :: Cexpg(:,:,:,:), CuvExpg(:,:,:)
    double complex, allocatable :: B_0(:,:,:), B_i(:,:,:), Shat(:,:,:,:)
    double complex, allocatable :: Buv_0(:,:,:), Buv_i(:,:,:)
    double complex :: Smod(2), Skl, CexpgAux
    double complex :: C0_coli, elimminf2_coli
    double precision, allocatable :: C00_err(:),Cij_err(:)
    double precision, allocatable :: C00_err2(:),Cij_err2(:)
    double precision :: B_err,B_max
    double precision :: maxCexpg(0:1,0:rmax+ordg_min+1,0:ordg_max),truncfacexp
    integer :: rmaxB,rmaxExp,gtrunc,r,n0,n1,n2,k,l,i,j,m,n,sgn,g,rg
    integer :: inds0(2), inds(2), inds2(2), ktlt(2)
    integer :: bin,nid(0:2)

    double complex, allocatable :: D_alt(:,:,:,:)



    ! allocation of B functions
    rmaxB = rmax + ordg_min
    allocate(B_0(0:rmaxB,0:rmaxB,0:rmaxB))
    allocate(Buv_0(0:rmaxB,0:rmaxB,0:rmaxB))
    allocate(B_i(0:rmaxB,0:rmaxB,2))
    allocate(Buv_i(0:rmaxB,0:rmaxB,2))
    

    ! determine binaries for B-coefficients
    k=0
    bin = 1
    do while (k.le.2)
      if (mod(id/bin,2).eq.0) then
        nid(k) = id+bin
        k = k+1
      end if
      bin = 2*bin
    end do

    call CalcB(B_0(:,0,:),Buv_0(:,0,:),p21,m12,m22,rmaxB,nid(0))
    call CalcB(B_i(:,:,1),Buv_i(:,:,1),p20,m02,m22,rmaxB,nid(1))
    call CalcB(B_i(:,:,2),Buv_i(:,:,2),p10,m02,m12,rmaxB,nid(2))

    ! shift of integration momentum in B_0
    B_max=0d0
    do n1=1,rmaxB
      do n2=0,rmaxB-n1
        n0 = (rmaxB-n1-n2)
        B_0(0:n0,n1,n2) = -B_0(0:n0,n1-1,n2)-B_0(0:n0,n1-1,n2+1)
        Buv_0(0:n0,n1,n2) = -Buv_0(0:n0,n1-1,n2)-Buv_0(0:n0,n1-1,n2+1)
        B_max = max(B_max,abs(B_0(0,n1,n2)))
      end do
    end do
    ! error estimate for B's
    B_max = max(B_max,maxval(abs(B_i(0,0:rmaxB,1:2))))
    B_err = acc_def_B*B_max


    ! determine (adjugated) Gram matrix
!    mm02 = elimminf2_coli(m02)
!    mm12 = elimminf2_coli(m12)
!    mm22 = elimminf2_coli(m22)
!    q10  = elimminf2_coli(p10)
!    q21  = elimminf2_coli(p21)
!    q20  = elimminf2_coli(p20)
!
!    q1q2 = (q10+q20-q21)
!    detZ = 4d0*q10*q20-q1q2*q1q2

!    if (abs(detZ/( 4d0*q10*q20 + q1q2*q1q2)).lt.1d-4) then
!      if (abs(q10-q20).lt.abs(q10-q21).and.  &
!          abs(q10-q20).lt.abs(q20-q21)) then
!        detZ  =  4d0*q10*q21 - (q10-q20+q21)*(q10-q20+q21)
!      end if
!    end if

!   write(*,*) 'Z = ',Z

!    Zadj(1,1) = 2d0*q20
!    Zadj(2,1) = -q1q2
!    Zadj(1,2) = -q1q2
!    Zadj(2,2) = 2d0*q10
!    f(1) = q10+mm02-mm12
!    f(2) = q20+mm02-mm22

!    Zadjf(1) = Zadj(1,1)*f(1)+Zadj(2,1)*f(2)
!    Zadjf(2) = Zadj(1,2)*f(1)+Zadj(2,2)*f(2)

!    maxZadj=maxval(abs(Zadj))
!    fmax   =maxval(abs(f))

    ! coefficients Shat defined in (5.13)
    allocate(Shat(0:rmaxB,0:rmaxB,0:rmaxB,2))

    do r=0,rmaxB
      do n0=0,r/2

        do n1=0,r-2*n0
          n2 = r-2*n0-n1
          Shat(n0,n1,n2,:) = -B_0(n0,n1,n2)
        end do

        k = r-2*n0
        Shat(n0,0,k,1) = Shat(n0,0,k,1) + B_i(n0,k,1)
        Shat(n0,k,0,2) = Shat(n0,k,0,2) + B_i(n0,k,2)

      end do
    end do

    ! choose reduction formulas with biggest denominators
    if (abs(Zadjf(1)).ge.abs(Zadjf(2))) then
      j = 1
    else 
      j = 2
    end if

    maxZadj = 0d0                 ! Zadj2f(k,n,l) = Zadf2(k,n,l,m)*f(m)
                                  ! Zadj2(n,m) ==  Zadf2(k,n,l,m)
    if (abs(Zadj(1,1)).gt.maxZadj) then
      maxZadj = abs(Zadj(1,1))
      k = 1
      l = 1
      inds2 = (/2,2/)
      Zadj2 = -1d0
    end if
    if (abs(Zadj(1,2)).gt.maxZadj) then
      maxZadj = abs(Zadj(1,2))
      k = 1
      l = 2
      inds2 = (/2,1/)
      Zadj2 =  1d0
    end if

    Zadjfj = Zadjf(j)
    Zadjkl = Zadj(k,l)
    Xtilde = Xadj(k,l)

!    write(*,*) 'CalcCgn Xtilde n',Xtilde,Xadj(1,1),Xadj(1,2),Xadj(2,2)


    ! allocation of array for det(Z)-expanded C-coefficients
    rmaxExp = rmaxB+1
    allocate(Cexpg(0:rmaxExp/2,0:rmaxExp-1,0:rmaxExp-1,0:ordg_max))

    ! calculate Cuv
    allocate(CuvExpg(0:rmaxExp,0:rmaxExp,0:rmaxExp))
    call CalcCuv(CuvExpg,Buv_0,mm02,f,rmaxExp,id)
    Cuv(0:rmax,0:rmax,0:rmax) = CuvExpg(0:rmax,0:rmax,0:rmax)

    ! allocate arrays for error propagation
    allocate(C00_err(0:rmaxExp))
    allocate(Cij_err(0:rmaxExp))
    allocate(C00_err2(0:rmaxExp))
    allocate(Cij_err2(0:rmaxExp))
     
    ! initialize accuracy estimates
    Cerr = acc_inf
    Cij_err =0d0
    C00_err =0d0

    Cerr2 = acc_inf
    Cij_err2 =0d0
    C00_err2 =0d0



!    maxZadj = maxval(abs(Zadj))
!    maxZadj2f = maxval(abs(f(inds2(1,:))*Zadj2(:)))

    ! truncation of expansion if calculated term larger than truncfacexp * previous term
    ! crucial for expansion parameters between 0.1 and 1 !!!
    truncfacexp = sqrt(fac_g) * truncfacC
    gtrunc = ordg_max 

! calculate C(n0,n1,n2) up to rank r for n0>0 and up to rank r-1 for n0=0
    rloop: do r=1,rmaxExp



      if (r.gt.rmax+gtrunc+1) exit rloop



      !!!!!!!!!!!!!!!!!!!!!!!!!!!!
      ! 0th-order coefficients
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!

      ! calculating
      ! C_00(a)0000..00 --> C_00(a)ij00..00 --> C_00(a)ijkl00..00 --> ... --> C_00(a)ijklmn..
      ! exploiting eq. (5.40)
      maxCexpg(1,r,0)=0d0
      do n0=r/2,1,-1
        do n1=0,r-2*n0
          n2=r-2*n0-n1

          inds0(1) = n1
          inds0(2) = n2
          
          CexpgAux = 2d0*Zadj(k,l)*B_0(n0-1,n1,n2)  & 
              + Xtilde*Cexpg(n0-1,n1,n2,0)  &
              + 4d0*Zadj(k,l)*CuvExpg(n0,n1,n2)
          
          inds = inds0
          inds(k) = inds(k)+1
          do i=1,2
            CexpgAux = CexpgAux + Zadj(i,l)*Shat(n0-1,inds(1),inds(2),i)
          end do
          
          do i=1,2
            inds = inds0
            inds(i) = inds(i)+1
            CexpgAux = CexpgAux - Zadj(k,l)*Shat(n0-1,inds(1),inds(2),i)
          end do
          
          n = inds2(1)
          m = inds2(2)
          
          Skl = f(n)*Shat(n0-1,inds0(1),inds0(2),m)
          
          inds = inds0
          if (inds(m).ge.1) then
            inds(m) = inds(m)-1
            Skl = Skl - 2d0*f(n)*inds0(m)*Cexpg(n0,inds(1),inds(2),0)
            if (inds(n).ge.1) then
              inds(n) = inds(n)-1
              Skl = Skl - 4d0*inds0(m)*(inds(n)+1)*Cexpg(n0+1,inds(1),inds(2),0)
            end if
          end if
          inds = inds0
          if (inds(n).ge.1) then
            inds(n) = inds(n)-1
            Skl = Skl + 2d0*inds0(n)*Shat(n0,inds(1),inds(2),m)  &
                - 2d0*f(m)*inds0(n)*Cexpg(n0,inds(1),inds(2),0)
          end if
          
          CexpgAux = CexpgAux - Zadj2*Skl
          
          Cexpg(n0,n1,n2,0) = CexpgAux/(2d0*Zadjkl)/(2d0*(r-n0)+1)
          
          if (n0.eq.1) then
            maxCexpg(1,r,0) =  maxCexpg(1,r,0) + abs(Cexpg(n0,n1,n2,0) )
          end if
          
          if (r-n0.le.rmax) then
            C(n0,n1,n2) = Cexpg(n0,n1,n2,0)
          end if

        end do
      end do

      ! calculate
      ! C_00ijkl.. --> C_aijkl..
      ! exploiting eq. (5.38)
      maxCexpg(0,r-1,0)=0d0
      do n1=0,r-1
        n2=r-1-n1

          Smod = Shat(0,n1,n2,:)
          if (n1.ge.1) then
            Smod(1) = Smod(1) - 2d0*n1*Cexpg(1,n1-1,n2,0)
          end if
          if (n2.ge.1) then
            Smod(2) = Smod(2) - 2d0*n2*Cexpg(1,n1,n2-1,0)
          end if

          Cexpg(0,n1,n2,0) = (Zadj(1,j)*Smod(1) +  Zadj(2,j)*Smod(2))/Zadjfj
          maxCexpg(0,r-1,0) =  maxCexpg(0,r-1,0) + abs(Cexpg(0,n1,n2,0))
          if (r-n0.le.rmax+1) then
            C(0,n1,n2) = Cexpg(0,n1,n2,0)
          end if



      end do



      if(r.le.rmax+1) then
!       Cerr(r-1) =  abs(detZ/Zadjfj)*maxCexpg(0,r-1,0)
        Cerr(r-1) =  fac_g*maxCexpg(0,r-1,0)
      end if

      ! error propagation from B's
      C00_err(r) = max(max(maxZadj*B_err,fmax*B_err)/abs(Zadjkl),B_err)  &
                   /(2*(2*r-1))        
      Cij_err(r-1)=maxZadj*max(B_err,2*C00_err(r))/abs(Zadjfj)

      C00_err2(r) = max(max(maxZadj*B_err,fmax*B_err)/abs(Zadjkl),B_err)  &
                   /(2*(2*r-1))        
      Cij_err2(r-1)=maxZadj*max(B_err,2*C00_err2(r))/abs(Zadjfj)



      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      ! higher order coefficients
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      rg = r
      gloop: do g=1,min(gtrunc,r-1)
        rg = rg-1

!        write(*,*) 'gloop ',g,rg

        ! calculating
        ! C_00(a)0000..00 --> C_00(a)ij00..00 --> C_00(a)ijkl00..00 --> ... --> C_00(a)ijklmn..
        ! exploiting eq. (5.40)
        maxCexpg(1,rg,g) = 0d0
        do n0=rg/2,1,-1
          do n1=0,rg-2*n0
            n2=rg-2*n0-n1

            inds0(1) = n1
            inds0(2) = n2
            
            inds = inds0
            inds(k) = inds(k)+1
            inds(l) = inds(l)+1
            CexpgAux = Xtilde*Cexpg(n0-1,n1,n2,g)  &
                - detZ*Cexpg(n0-1,inds(1),inds(2),g-1)
            
            
            n = inds2(1)
            m = inds2(2)
            
            Skl = 0d0
            
            inds = inds0
            if (inds(m).ge.1) then
              inds(m) = inds(m)-1
              Skl = Skl - 2d0*f(n)*inds0(m)*Cexpg(n0,inds(1),inds(2),g)
              if (inds(n).ge.1) then
                inds(n) = inds(n)-1
                Skl = Skl - 4d0*inds0(m)*(inds(n)+1)*Cexpg(n0+1,inds(1),inds(2),g)
              end if
            end if
            inds = inds0
            if (inds(n).ge.1) then
              inds(n) = inds(n)-1
              Skl = Skl - 2d0*f(m)*inds0(n)*Cexpg(n0,inds(1),inds(2),g)
            end if
            
            CexpgAux = CexpgAux - Zadj2*Skl
            
            Cexpg(n0,n1,n2,g) = CexpgAux/(2d0*Zadjkl)/(2d0*(rg-n0)+1)
            
            
            if(n0.eq.1) then
              maxCexpg(1,rg,g) =  maxCexpg(1,rg,g) + abs(Cexpg(n0,n1,n2,g))
              
              if (g.eq.1.and.abs(Cexpg(n0,n1,n2,g)).gt.          & 
                  truncfacexp*max(1d0,maxCexpg(1,rg,g-1)) .or.   &
                  g.ge.2.and.abs(Cexpg(n0,n1,n2,g)).gt.          &
                  truncfacexp*maxCexpg(1,rg,g-1)) then


                
                gtrunc = g-1
                exit gloop
              end if
            end if
            
          end do
        end do


        do n0=rg/2,1,-1
          if (rg-n0.le.rmax) then
            do n1=0,rg-2*n0
              n2=rg-2*n0-n1
              C(n0,n1,n2) = C(n0,n1,n2) + Cexpg(n0,n1,n2,g)
            end do
          end if
        end do

!        write(*,*) 'CalcCgn after it1 ',rg

        ! calculate
        ! C_00ijkl.. --> C_aijkl..
        ! exploiting eq. (5.38)

!        write(*,*) 'CalcCgn maxCexp',rg-1,g-1,maxCexpg(0,rg-1,g-1)

        maxCexpg(0,rg-1,g) = 0d0
        do n1=0,rg-1
          n2=rg-1-n1

          Smod = 0d0
          if (n1.ge.1) then
            Smod(1) = Smod(1) - 2d0*n1*Cexpg(1,n1-1,n2,g)
          end if
          if (n2.ge.1) then
            Smod(2) = Smod(2) - 2d0*n2*Cexpg(1,n1,n2-1,g)
          end if
          
          inds(1) = n1
          inds(2) = n2
          inds(j) = inds(j)+1
          Cexpg(0,n1,n2,g) = (Zadj(1,j)*Smod(1) +  Zadj(2,j)*Smod(2)  &
              - detZ*Cexpg(0,inds(1),inds(2),g-1))/Zadjfj
          
          maxCexpg(0,rg-1,g) =  maxCexpg(0,rg-1,g) + abs(Cexpg(0,n1,n2,g))
          
!              if(n1.eq.0.and.n2.eq.1) then 
!                write(*,*) 'C2(2,3)= ',g,Cexpg(0,n1,n2,g)
!                write(*,*) 'C2(2,3)= ',Zadj(1,j)*Smod(1)/Zadjfj,  Zadj(2,j)*Smod(2)/Zadjfj,  &
!                                - detZ*Cexpg(0,inds(1),inds(2),inds(3),g-1)/Zadjfj
!                write(*,*) 'C2(2,3)= ',inds(1),inds(2),         &
!                                - detZ/Zadjfj,Cexpg(0,inds(1),inds(2),g-1)
!              end if
          
          if (g.eq.1.and.abs(Cexpg(0,n1,n2,g)).gt.                     &
!     corrected 02.07.2018      
              truncfacexp*max(1/m2scale,maxCexpg(0,rg-1,g-1))    .or.    &
!              truncfacexp*max(1/m2max,maxCexpg(0,rg-1,g-1))    .or.    &
              g.ge.2.and.abs(Cexpg(0,n1,n2,g)).gt.                     &
              truncfacexp*maxCexpg(0,rg-1,g-1)) then
            

            gtrunc = g-1
            exit gloop
          end if

        end do

        ! error propagation from B's
        if(rg.gt.1)then
          C00_err(rg) = max(C00_err(rg),                    &
              max( abs(m02)*Cij_err(rg-2),                             &
              max(adetZ*Cij_err(rg),fmax**2*Cij_err(rg-2),fmax*C00_err(rg-1))/abs(Zadjkl) ) &
                   /(2*(2*rg-1))     )   
        end if
        Cij_err(rg-1) = max(Cij_err(rg-1),max(2*maxZadj*C00_err(rg),adetZ*Cij_err(rg))/abs(Zadjfj) ) 

        if(rg.gt.1)then
          C00_err2(rg) = max(C00_err2(rg),                    &
              max( abs(m02)*Cij_err2(rg-2),                             &
              max(adetZ*Cij_err2(rg),fmax**2*Cij_err2(rg-2),fmax*C00_err2(rg-1))/abs(Zadjkl) ) &
                   /(2*(2*rg-1))     )   
        end if
        Cij_err2(rg-1) = max(Cij_err2(rg-1),max(2*maxZadj*C00_err2(rg),adetZ*Cij_err2(rg))/abs(Zadjfj) ) 

!      write(*,*) 'CalcCg g: ',r,adetZ/abs(Zadjfj),C00_err(rg),B_err
!      write(*,*) 'CalcCg g: Cij_err=',rg-1,Cij_err(rg-1)



!        write(*,*) 'CalcCgn after it1 ',rg
        if ((rg.le.rmax+1)) then
          Cerr(rg-1) = 0d0
          do n1=0,rg-1
            n2 = rg-1-n1
            C(0,n1,n2) = C(0,n1,n2) + Cexpg(0,n1,n2,g)
            if(abs(Cexpg(0,n1,n2,g-1)).ne.0d0) then
!             Cerr(rg-1)=max(Cerr(rg-1),abs(Cexpg(0,n1,n2,g))**2/abs(Cexpg(0,n1,n2,g-1)))
              Cerr(rg-1)=max(Cerr(rg-1),abs(Cexpg(0,n1,n2,g))*min(1d0,abs(Cexpg(0,n1,n2,g))/abs(Cexpg(0,n1,n2,g-1))))
            else
              Cerr(rg-1)=max(Cerr(rg-1),abs(Cexpg(0,n1,n2,g)))
            end if

!             write(*,*) 'CalcCg err',r,rg,n1,n2,Cerr(rg-1),abs(Cexpg(0,n1,n2,g))**2/abs(Cexpg(0,n1,n2,g-1)) &
!                     ,abs(Cexpg(0,n1,n2,g)),abs(Cexpg(0,n1,n2,g-1))

          end do

          ! if error from B's larger than error from expansion stop expansion
          if(Cij_err(rg-1).gt.Cerr(rg-1)) then
             gtrunc = min(g,gtrunc)
             


          end if
        end if

      end do gloop



      Cerr2 = max(Cerr,Cij_err2(0:rmax))
      Cerr = max(Cerr,Cij_err(0:rmax))



!      if(maxval(Cerr).le.acc_req_C*abs(C(0,0,0))) exit       ! changed 28.01.15
      ! check if target precision already reached

      if(maxval(Cerr-acc_req_Cr*abs(C(0,0,0))).le.0d0) then
        do rg=r+1,rmax

!          write(*,*) 'CalcCgn exit rloop  =',rg,r,rmax

          do n0=0,rg/2
            do n1=0,rg-2*n0
              C(n0,n1,rg-2*n0-n1)=0d0
            end do
          end do
        end do
        if(r.le.rmax) then
          do n1=0,r
            C(0,n1,r-n1)=0d0
          end do
        end if

       exit rloop
      end if

    end do rloop

      ! reduction formula (5.10) for n0+n1+n2=r, n0>0 
    do r=rmax+1,2*rmax
      do n0=r-rmax,r/2
        do n1=0,r-2*n0
          n2 = r-2*n0-n1
          C(n0,n1,n2) = (B_0(n0-1,n1,n2) + 2*mm02*C(n0-1,n1,n2) + 4*Cuv(n0,n1,n2) &
                        + f(1)*C(n0-1,n1+1,n2) + f(2)*C(n0-1,n1,n2+1)) / (2*r)
        end do
      end do
    end do



                           
!   write(*,*) 'CalcCgn out',(((C((r-n1-n2)/2,n1,n2),n2=0,r-n1),n1=0,r),r=0,rmax)


!    write(*,*) 'CalcCgn Cerr ',Cerr
!    write(*,*) 'CalcCgn Cerr2',Cerr2

  end subroutine CalcCgn




  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !  subroutine CalcCg(C,Cuv,p10,p21,p20,m02,m12,m22,rmax,ordg_min,ordg_max,id,Cerr,acc_req_Cr,Cerr2)
  !
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  
  subroutine CalcCg(C,Cuv,p10,p21,p20,m02,m12,m22,rmax,ordg_min,ordg_max,id,Cerr,acc_req_Cr,Cerr2)
  
    use globalC
  
    integer, intent(in) :: rmax,ordg_min,ordg_max,id
    double complex, intent(in) :: p10,p21,p20,m02,m12,m22
    double complex, intent(out) :: C(0:rmax,0:rmax,0:rmax)
    double complex, intent(out) :: Cuv(0:rmax,0:rmax,0:rmax)
    double precision, intent(out) :: Cerr(0:rmax),Cerr2(0:rmax)
    double precision, intent(in) :: acc_req_Cr(0:rmax)
    double complex :: Xtilde,Zkl,Zadjfj
    double complex, allocatable :: Cexpg(:,:,:,:), CuvExpg(:,:,:)
    double complex, allocatable :: B_0(:,:,:), B_i(:,:,:), Shat(:,:,:,:)
    double complex, allocatable :: Buv_0(:,:,:), Buv_i(:,:,:)
    double complex :: Smod(2), Skl
    double complex :: C0_coli, elimminf2_coli
    double precision, allocatable :: C00_err(:),Cij_err(:)
    double precision, allocatable :: C00_err2(:),Cij_err2(:)
    double precision :: B_err,B_max
    double precision :: maxCexpg(0:1,0:rmax+ordg_min+1,0:ordg_max),truncfacexp
    integer :: rmaxB,rmaxExp,gtrunc,r,n0,n1,n2,k,l,j,sgn,g,rg,mr
    integer :: inds0(2), inds(2), ktlt(2)
    integer :: bin,nid(0:2)



    ! write(*,*) 'LH: CalcCg, ord', ordg_min 
    ! calculation B-coefficients
    rmaxB = rmax + ordg_min
    allocate(B_0(0:rmaxB,0:rmaxB,0:rmaxB))
    allocate(Buv_0(0:rmaxB,0:rmaxB,0:rmaxB))
    allocate(B_i(0:rmaxB,0:rmaxB,2))
    allocate(Buv_i(0:rmaxB,0:rmaxB,2))

    ! determine binaries for B-coefficients
    k=0
    bin = 1
    do while (k.le.2)
      if (mod(id/bin,2).eq.0) then
        nid(k) = id+bin
        k = k+1
      end if
      bin = 2*bin
    end do

    call CalcB(B_0(:,0,:),Buv_0(:,0,:),p21,m12,m22,rmaxB,nid(0))
    call CalcB(B_i(:,:,1),Buv_i(:,:,1),p20,m02,m22,rmaxB,nid(1))
    call CalcB(B_i(:,:,2),Buv_i(:,:,2),p10,m02,m12,rmaxB,nid(2))

    ! shift of integration momentum in B_0
    B_max=0d0
    do n1=1,rmaxB
      do n2=0,rmaxB-n1
        n0 = (rmaxB-n1-n2)
        B_0(0:n0,n1,n2) = -B_0(0:n0,n1-1,n2)-B_0(0:n0,n1-1,n2+1)
        Buv_0(0:n0,n1,n2) = -Buv_0(0:n0,n1-1,n2)-Buv_0(0:n0,n1-1,n2+1)
        B_max = max(B_max,abs(B_0(0,n1,n2)))
      end do
    end do
    ! error estimate for B's
    B_max = max(B_max,maxval(abs(B_i(0,0:rmaxB,1:2))))
    B_err = acc_def_B*B_max
    
    ! determine (adjugated) Gram matrix
!    mm02 = elimminf2_coli(m02)
!    mm12 = elimminf2_coli(m12)
!    mm22 = elimminf2_coli(m22)
!    q10  = elimminf2_coli(p10)
!    q21  = elimminf2_coli(p21)
!    q20  = elimminf2_coli(p20)
!
!    q1q2 = (q10+q20-q21)
!    detZ = 4d0*q10*q20-q1q2*q1q2

!    if (abs(detZ/( 4d0*q10*q20 + q1q2*q1q2)).lt.1d-4) then
!      if (abs(q10-q20).lt.abs(q10-q21).and.  &
!          abs(q10-q20).lt.abs(q20-q21)) then
!        detZ  =  4d0*q10*q21 - (q10-q20+q21)*(q10-q20+q21)
!      end if
!    end if

!   write(*,*) 'Z = ',Z

!    Zadj(1,1) = 2d0*q20
!    Zadj(2,1) = -q1q2
!    Zadj(1,2) = -q1q2
!    Zadj(2,2) = 2d0*q10
!    f(1) = q10+mm02-mm12
!    f(2) = q20+mm02-mm22

!    Zadjf(1) = Zadj(1,1)*f(1)+Zadj(2,1)*f(2)
!    Zadjf(2) = Zadj(1,2)*f(1)+Zadj(2,2)*f(2)

!    commented out 2.9.2017
!    maxZadj=maxval(abs(Zadj))
!    fmax   =maxval(abs(f))

    ! coefficients Shat defined in (5.13)
    allocate(Shat(0:rmaxB,0:rmaxB,0:rmaxB,2))

    do r=0,rmaxB
      do n0=0,r/2

        do n1=0,r-2*n0
          n2 = r-2*n0-n1
          Shat(n0,n1,n2,:) = -B_0(n0,n1,n2)
        end do

        k = r-2*n0
        Shat(n0,0,k,1) = Shat(n0,0,k,1) + B_i(n0,k,1)
        Shat(n0,k,0,2) = Shat(n0,k,0,2) + B_i(n0,k,2)

      end do
    end do


    ! choose reduction formulas with biggest denominators
    if (abs(Zadjf(1)).ge.abs(Zadjf(2))) then
      j = 1
    else 
      j = 2
    end if

    if (abs(Z(1,1)).ge.abs(Z(2,2))) then
      if (abs(Z(1,1)).ge.abs(Z(1,2))) then
        k = 1
        l = 1
        sgn = 1
        ktlt = (/ 0,2 /)
      else
        k = 1
        l = 2
        sgn = -1
        ktlt = (/ 1,1 /)
      end if
    else        
      if (abs(Z(2,2)).ge.abs(Z(1,2))) then
        k = 2
        l = 2
        sgn = 1
        ktlt = (/ 2,0 /)
      else 
        k = 1
        l = 2
        sgn = -1
        ktlt = (/ 1,1 /)
      end if
    end if

    Zadjfj = Zadjf(j)
    Zkl = Z(k,l)
    if(k.eq.l) then
      Xtilde = Xadj(3-k,3-l)      ! subroutine uses Z instead of Zadj
    else                          ! -> exchange indices 1 and 2
      Xtilde = -Xadj(3-k,3-l)     ! -> minus sign for k \ne l
    end if

!   write(*,*) 'CalcCg Xtilde n',Xtilde,Xadj(1,1),Xadj(1,2),Xadj(2,2)

!   write(*,*) 'Xtilde =',Xtilde,k,l

    ! allocation of array for det(Z)-expanded C-coefficients
    rmaxExp = rmaxB+1
    allocate(Cexpg(0:rmaxExp/2,0:rmaxExp-1,0:rmaxExp-1,0:ordg_max))
   

    ! calculate Cuv
    allocate(CuvExpg(0:rmaxExp,0:rmaxExp,0:rmaxExp))
    call CalcCuv(CuvExpg,Buv_0,mm02,f,rmaxExp,id)
    Cuv(0:rmax,0:rmax,0:rmax) = CuvExpg(0:rmax,0:rmax,0:rmax)

    ! allocate arrays for error propagation
    allocate(C00_err(0:rmaxExp))
    allocate(Cij_err(0:rmaxExp))
    allocate(C00_err2(0:rmaxExp))
    allocate(Cij_err2(0:rmaxExp))
    
    ! initialize accuracy estimates
    Cerr = acc_inf
    Cij_err = 0d0
    C00_err = 0d0

    Cerr2 = acc_inf
    Cij_err2 = 0d0
    C00_err2 = 0d0

    ! truncation of expansion if calculated term larger than truncfacexp * previous term
    ! crucial for expansion parameters between 0.1 and 1 !!!
!    truncfacexp = sqrt(abs(detZ/Zadjfj)) * truncfacC
    truncfacexp = sqrt(fac_g) * truncfacC
    gtrunc = ordg_max 

! calculate C(n0,n1,n2) up to rank r for n0>0 and up to rank r-1 for n0=0
    rloop: do r=1,rmaxExp

      if (r.gt.rmax+gtrunc+1) exit rloop

!     write(*,*) 'CalcCg rloop',r,rmaxExp,gtrunc

      !!!!!!!!!!!!!!!!!!!!!!!!!!!!
      ! 0th-order coefficients
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!

      ! calculating
      ! C_00(a)0000..00 --> C_00(a)ij00..00 --> C_00(a)ijkl00..00 --> ... --> C_00(a)ijklmn..
      ! exploiting eq. (5.40)
      maxCexpg(1,r,0)=0d0
      do n0=r/2,1,-1
        do n1=0,r-2*n0
          n2=r-2*n0-n1

          inds0(1) = n1
          inds0(2) = n2
          Skl = 0d0
          inds = inds0
          if (inds(k).ge.1) then
            inds(k) = inds(k)-1
            Skl = Skl - 2d0*f(l)*inds0(k)*Cexpg(n0,inds(1),inds(2),0)
            if (inds(l).ge.1) then
              inds(l) = inds(l)-1
              Skl = Skl - 4d0*inds0(k)*(inds(l)+1)*Cexpg(n0+1,inds(1),inds(2),0)
            end if
          end if
          inds = inds0
          if (inds(l).ge.1) then
            inds(l) = inds(l)-1
            Skl = Skl + 2d0*inds0(l)*Shat(n0,inds(1),inds(2),k)  &
                      - 2d0*f(k)*inds0(l)*Cexpg(n0,inds(1),inds(2),0)
          end if

          Cexpg(n0,n1,n2,0) = (2d0*Zkl*B_0(n0-1,n1,n2) + Xtilde*Cexpg(n0-1,n1,n2,0)  &
                                 - Z(1,k)*Shat(n0-1,n1+1,n2,l) - Z(2,k)*Shat(n0-1,n1,n2+1,l)  &
                                 + f(l)*Shat(n0-1,n1,n2,k) + 4d0*Zkl*CuvExpg(n0,n1,n2) + Skl)  &
                                 /(2d0*Zkl)/(2d0*(r-n0)+1d0)

          if (n0.eq.1) then
            maxCexpg(1,r,0) =  maxCexpg(1,r,0) + abs(Cexpg(n0,n1,n2,0))
          end if

          if (r-n0.le.rmax) then
            C(n0,n1,n2) = Cexpg(n0,n1,n2,0)
          end if

        end do
      end do

      ! calculate
      ! C_00ijkl.. --> C_aijkl..
      ! exploiting eq. (5.38)
      maxCexpg(0,r-1,0)=0d0
      do n1=0,r-1
        n2 = r-1-n1

        Smod = Shat(0,n1,n2,:)
        if (n1.ge.1) then
          Smod(1) = Smod(1) - 2d0*n1*Cexpg(1,n1-1,n2,0)
        end if
        if (n2.ge.1) then
          Smod(2) = Smod(2) - 2d0*n2*Cexpg(1,n1,n2-1,0)
        end if

        Cexpg(0,n1,n2,0) = (Zadj(1,j)*Smod(1) +  Zadj(2,j)*Smod(2))/Zadjfj
        
        maxCexpg(0,r-1,0) =  maxCexpg(0,r-1,0) + abs(Cexpg(0,n1,n2,0))
        if (r-n0.le.rmax+1) then
          C(0,n1,n2) = Cexpg(0,n1,n2,0)
        end if

      end do

      if(r.le.rmax+1) then
!       Cerr(r-1) =  abs(detZ/Zadjfj)*maxCexpg(0,r-1,0)
        Cerr(r-1) =  fac_g*maxCexpg(0,r-1,0)
      end if

      ! error propagation from B's
      C00_err(r) = max(max(maxZadj*B_err,fmax*B_err)/abs(Zkl),B_err)  &
                   /(2*(2*r-1))        
      Cij_err(r-1)=maxZadj*max(B_err,2*C00_err(r))/abs(Zadjfj)

      C00_err2(r) = max(max(maxZadj*B_err,fmax*B_err)/abs(Zkl),B_err)  &
                   /(2*(2*r-1))        
      Cij_err2(r-1)=maxZadj*max(B_err,2*C00_err2(r))/abs(Zadjfj)

!      write(*,*) 'CalcCg after 0: ',maxZadj/abs(Zadjfj),C00_err(r),B_err
!      write(*,*) 'CalcCg after 0: Cij_err=',r-1,Cij_err(r-1)


      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      ! higher order coefficients
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      rg = r
      gloop: do g=1,min(gtrunc,r-1)
        rg = rg-1
!       write(*,*) 'CalcCg gloop',g,rg

        ! calculating for rank=rmaxB+1
        ! C_00(a)0000..00 --> C_00(a)ij00..00 --> C_00(a)ijkl00..00 --> ... --> C_00(a)ijklmn..
        ! exploiting eq. (5.40)
        maxCexpg(1,rg,g) = 0d0
        do n0=rg/2,1,-1
          do n1=0,rg-2*n0
            n2=rg-2*n0-n1

            inds0(1) = n1
            inds0(2) = n2
            Skl = 0d0
            inds = inds0
            if (inds(k).ge.1) then
              inds(k) = inds(k)-1
              Skl = Skl - 2d0*f(l)*inds0(k)*Cexpg(n0,inds(1),inds(2),g)
              if (inds(l).ge.1) then
                inds(l) = inds(l)-1
                Skl = Skl - 4d0*inds0(k)*(inds(l)+1)*Cexpg(n0+1,inds(1),inds(2),g)
              end if
              inds = inds0
            end if
            if (inds(l).ge.1) then
              inds(l) = inds(l)-1
              Skl = Skl - 2d0*f(k)*inds0(l)*Cexpg(n0,inds(1),inds(2),g)
            end if
         
            inds = inds0 + ktlt
            Cexpg(n0,n1,n2,g) = (Xtilde*Cexpg(n0-1,n1,n2,g) + Skl  &
                                   - detZ*sgn*Cexpg(n0-1,inds(1),inds(2),g-1))  &
                                     /(2d0*Zkl)/(2d0*(rg-n0)+1d0)           
            if(n0.eq.1) then
              maxCexpg(1,rg,g) =  maxCexpg(1,rg,g) + abs(Cexpg(n0,n1,n2,g))
              
              if (g.eq.1.and.abs(Cexpg(n0,n1,n2,g)).gt.        &
                  truncfacexp*max(1d0,maxCexpg(1,rg,g-1)).or.  &
                  g.ge.2.and.abs(Cexpg(n0,n1,n2,g)).gt.        &  
                  truncfacexp*maxCexpg(1,rg,g-1)) then


!               write(*,*) 'CalcCg exit gloop',n0,n1,n2,g,abs(Cexpg(n0,n1,n2,g)),maxCexpg(1,rg,g-1)

                gtrunc = g-1
                exit gloop
              end if
            end if
          end do
        end do

!        write(*,*) 'Calcg: rg,g,acc',rg,g,acc


        do n0=rg/2,1,-1
          if (rg-n0.le.rmax) then
            do n1=0,rg-2*n0
              n2=rg-2*n0-n1
              C(n0,n1,n2) = C(n0,n1,n2) + Cexpg(n0,n1,n2,g)
            end do
          end if
        end do


        ! calculate
        ! C_000000..00 --> C_i0000..00 --> C_ij00..00 --> ... --> C_ijk..
        ! exploiting eq. (5.38)
        maxCexpg(0,rg-1,g) = 0d0
        do n1=0,rg-1
          n2 = rg-1-n1
          
          Smod = 0d0
          if (n1.ge.1) then
            Smod(1) = Smod(1) - 2d0*n1*Cexpg(1,n1-1,n2,g)
          end if
          if (n2.ge.1) then
            Smod(2) = Smod(2) - 2d0*n2*Cexpg(1,n1,n2-1,g)
          end if
          
          inds(1) = n1
          inds(2) = n2
          inds(j) = inds(j)+1
          
          Cexpg(0,n1,n2,g) = (Zadj(1,j)*Smod(1) +  Zadj(2,j)*Smod(2)  &
              - detZ*Cexpg(0,inds(1),inds(2),g-1))/Zadjfj
          
          maxCexpg(0,rg-1,g) =  maxCexpg(0,rg-1,g) + abs(Cexpg(0,n1,n2,g))

          if (g.eq.1.and.abs(Cexpg(0,n1,n2,g)).gt.                     &
!     corrected 02.07.2018
              truncfacexp*max(1/m2scale,maxCexpg(0,rg-1,g-1))    .or.    &
!              truncfacexp*max(1/m2max,maxCexpg(0,rg-1,g-1))    .or.    &
              g.ge.2.and.abs(Cexpg(0,n1,n2,g)).gt.                     &
              truncfacexp*maxCexpg(0,rg-1,g-1)) then


            
            gtrunc = g-1


!           write(*,*) 'CalcCg exit gloop',rmax,g,rmaxExp

            exit gloop
          end if

        end do

        ! error propagation from B's
        if(rg.gt.1)then
!          C00_err(rg) = max(C00_err(rg),                    &
!              max( abs(m02)*Cij_err(rg-2),             &     
!              max(adetZ*Cij_err(rg),fmax**2*Cij_err(rg-2),fmax*C00_err(rg-1))/abs(Zkl) ) &
!                   /(2*(2*rg-1))     )   
!24.04.15 ->
!          C00_err(rg) = max(C00_err(rg),                    &
!              max( abs(m02)*Cij_err(rg-2),             &     
!              max(adetZ*Cij_err(rg),abs(Xtilde)*Cij_err(rg-2),fmax*C00_err(rg-1))/abs(Zkl) ) &
!                   /(2*(2*rg-1))     )   
!06.05.15 -> 
          C00_err(rg) = max(C00_err(rg),                    &
              max(adetZ*Cij_err(rg),abs(Xtilde)*Cij_err(rg-2),fmax*C00_err(rg-1))/abs(Zkl)  &
                   /(2*(2*rg-1))     )   
        end if
        Cij_err(rg-1) = max(Cij_err(rg-1),max(2*maxZadj*C00_err(rg),adetZ*Cij_err(rg))/abs(Zadjfj) ) 

        if(rg.gt.1)then
          C00_err2(rg) = max(C00_err2(rg),                    &
              max(adetZ*Cij_err2(rg),abs(Xtilde)*Cij_err2(rg-2),fmax*C00_err2(rg-1))/abs(Zkl)  &
                   /(2*(2*rg-1))     )   
        end if
        Cij_err2(rg-1) = max(Cij_err2(rg-1),max(2*maxZadj*C00_err2(rg),adetZ*Cij_err2(rg))/abs(Zadjfj) ) 

!      write(*,*) 'CalcCg g: ',r,adetZ/abs(Zadjfj),C00_err(rg),B_err
!      write(*,*) 'CalcCg g: Cij_err=',rg-1,Cij_err(rg-1)




        if ((rg.le.rmax+1)) then
          Cerr(rg-1) = 0d0
          do n1=0,rg-1
            n2 = rg-1-n1
            C(0,n1,n2) = C(0,n1,n2) + Cexpg(0,n1,n2,g)
            if(abs(Cexpg(0,n1,n2,g-1)).ne.0d0) then
!             Cerr(rg-1)=max(Cerr(rg-1),abs(Cexpg(0,n1,n2,g))**2/abs(Cexpg(0,n1,n2,g-1)))
              Cerr(rg-1)=max(Cerr(rg-1),abs(Cexpg(0,n1,n2,g))*min(1d0,abs(Cexpg(0,n1,n2,g))/abs(Cexpg(0,n1,n2,g-1))))
            else
              Cerr(rg-1)=max(Cerr(rg-1),abs(Cexpg(0,n1,n2,g)))
            end if

!            write(*,*) 'CalcCg err',r,rg,n1,n2,Cerr(rg-1),abs(Cexpg(0,n1,n2,g))**2/abs(Cexpg(0,n1,n2,g-1)) &
!                     ,abs(Cexpg(0,n1,n2,g)),abs(Cexpg(0,n1,n2,g-1))

          end do

          ! if error from B's larger than error from expansion stop expansion
          if(Cij_err(rg-1).gt.Cerr(rg-1)) then
             gtrunc = min(g,gtrunc)
             

!            write(*,*) 'CalcCg exit err',r,g,gtrunc

          end if
        end if

      end do gloop





      Cerr2 = max(Cerr,Cij_err2(0:rmax))
      Cerr = max(Cerr,Cij_err(0:rmax))


!     write(*,*) 'CalcCg Cerr =',r,Cerr,maxval(Cerr)
!     write(*,*) 'CalcCg areq =',acc_req_Cr*abs(C(0,0,0))
!     write(*,*) 'CalcCg Cex  =',maxval(Cerr-acc_req_Cr*abs(C(0,0,0)))

!       do mr = 15,min(r,rmax)
!         do n0=mr/2,1,-1
!           do n1=0,mr-2*n0
!             n2=mr-2*n0-n1
!               write(*,*) 'CalcCg n5 order  ',r,rg,mr,n0,n1,n2
!               write(*,*) 'CalcCg n5 order C',C(n0,n1,n2)
!           end do
!         end do 
!       end do
!       do mr = 15,min(r-1,rmax)
!         n0=0
!           do n1=0,mr
!             n2=mr-n1
!               write(*,*) 'CalcCg n5 order  ',r,rg,mr,n0,n1,n2
!               write(*,*) 'CalcCg n5 order C',C(n0,n1,n2)
!           end do
!       end do

      ! check if target precision already reached
!      if(maxval(Cerr-acc_req_Cr*abs(C(0,0,0))).le.0d0) exit         ! changed 28.01.15

      if(maxval(Cerr-acc_req_Cr*abs(C(0,0,0))).le.0d0) then
        do rg=r+1,rmax

!          write(*,*) 'CalcCg exit rloop  =',rg,r,rmax

          do n0=0,rg/2
            do n1=0,rg-2*n0
              C(n0,n1,rg-2*n0-n1)=0d0
            end do
          end do
        end do
        if(r.le.rmax) then
          do n1=0,r
            C(0,n1,r-n1)=0d0
          end do
        end if

!       write(*,*) 'CalcCg exit rloop  =',r,rmax,rg

        exit rloop
      end if

    end do rloop

      ! reduction formula (5.10) for n0+n1+n2=r, n0>0 
    do r=rmax+1,2*rmax
      do n0=r-rmax,r/2
        do n1=0,r-2*n0
          n2 = r-2*n0-n1
          C(n0,n1,n2) = (B_0(n0-1,n1,n2) + 2*mm02*C(n0-1,n1,n2) + 4*Cuv(n0,n1,n2) &
                        + f(1)*C(n0-1,n1+1,n2) + f(2)*C(n0-1,n1,n2+1)) / (2*r)
        end do
      end do
    end do

!       do mr = 15,rmax
!         do n0=mr/2,1,-1
!           do n1=0,mr-2*n0
!             n2=mr-2*n0-n1
!               write(*,*) 'CalcCg n6 order  ',r,rg,mr,n0,n1,n2
!               write(*,*) 'CalcCg n6 order C',C(n0,n1,n2)
!           end do
!         end do 
!       end do
!       do mr = 15,rmax
!         n0=0
!           do n1=0,mr
!             n2=mr-n1
!               write(*,*) 'CalcCg n6 order  ',r,rg,mr,n0,n1,n2
!               write(*,*) 'CalcCg n6 order C',C(n0,n1,n2)
!           end do
!       end do


                 


!    write(*,*) 'CalcCg Cerr ',Cerr
!    write(*,*) 'CalcCg Cerr2',Cerr2

  end subroutine CalcCg


  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !  subroutine CalcCgr(C,Cuv,p10,p21,p20,m02,m12,m22,rmax,ordgr_min,ordgr_max,id,Cerr,acc_req_Cr,Cerr2)
  !
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  
  subroutine CalcCgr(C,Cuv,p10,p21,p20,m02,m12,m22,rmax,ordgr_min,ordgr_max,id,Cerr,acc_req_Cr,Cerr2)
  
    use globalC

    integer, intent(in) :: rmax,ordgr_min,ordgr_max,id
    double complex, intent(in) :: p10,p21,p20,m02,m12,m22
    double complex, intent(out) :: C(0:rmax,0:rmax,0:rmax)
    double complex, intent(out) :: Cuv(0:rmax,0:rmax,0:rmax)
    double precision, intent(out) :: Cerr(0:rmax),Cerr2(0:rmax)
    double precision, intent(in) :: acc_req_Cr(0:rmax)
    double complex, allocatable :: B_0(:,:,:), B_i(:,:,:), Shat(:,:,:,:)
    double complex, allocatable :: Buv_0(:,:,:), Buv_i(:,:,:)
    double precision :: B_err,B_max
    double complex :: Zadjfj,Zadj2(2,2), Zadjkl, Zadj2f(2,2,2)
    double complex, allocatable :: Cexpgr(:,:,:,:), CuvExpgr(:,:,:)
    double complex :: Smod(2), Skl, Caux
    double complex :: elimminf2_coli
    double precision, allocatable :: C00_err(:),Cij_err(:)
    double precision, allocatable :: C00_err2(:),Cij_err2(:)
    double precision :: maxZadj2f
    double precision :: maxCexpgr(0:1,0:2*(rmax+ordgr_min),0:ordgr_max),truncfacexp
    integer :: rmaxB,rmaxExp,gtrunc,r,n0,n1,n2,k,l,i,j,m,n,g,rg,lt,ltt,nn,nntt
    integer :: inds0(2), inds1(2), inds(2)
    integer :: bin,nid(0:2)



    ! allocation of B functions
    rmaxB = 2*rmax + 2*ordgr_min
    allocate(B_0(0:rmaxB,0:rmaxB,0:rmaxB))
    allocate(Buv_0(0:rmaxB,0:rmaxB,0:rmaxB))
    allocate(B_i(0:rmaxB,0:rmaxB,2))
    allocate(Buv_i(0:rmaxB,0:rmaxB,2))
    
    ! determine binaries for B-coefficients
    k=0
    bin = 1
    do while (k.le.2)
      if (mod(id/bin,2).eq.0) then
        nid(k) = id+bin
        k = k+1
      end if
      bin = 2*bin
    end do

    call CalcB(B_0(:,0,:),Buv_0(:,0,:),p21,m12,m22,rmaxB,nid(0))
    call CalcB(B_i(:,:,1),Buv_i(:,:,1),p20,m02,m22,rmaxB,nid(1))
    call CalcB(B_i(:,:,2),Buv_i(:,:,2),p10,m02,m12,rmaxB,nid(2))

   ! shift of integration momentum in B_0
    B_max=0d0
    do n1=1,rmaxB
      do n2=0,rmaxB-n1
        n0 = (rmaxB-n1-n2)
        B_0(0:n0,n1,n2) = -B_0(0:n0,n1-1,n2)-B_0(0:n0,n1-1,n2+1)
        Buv_0(0:n0,n1,n2) = -Buv_0(0:n0,n1-1,n2)-Buv_0(0:n0,n1-1,n2+1)
      end do
    end do
    B_max = max(B_max,maxval(abs(B_i(0,0:rmaxB,1:2))))
    B_err = acc_def_B*B_max

    ! calculate adjugated Gram matrix
!    mm02 = elimminf2_coli(m02)
!    mm12 = elimminf2_coli(m12)
!    mm22 = elimminf2_coli(m22)
!    mm32 = elimminf2_coli(m32)
!    q10  = elimminf2_coli(p10)
!    q21  = elimminf2_coli(p21)
!    q32  = elimminf2_coli(p32)
!    q30  = elimminf2_coli(p30)
!    q31  = elimminf2_coli(p31)
!    q20  = elimminf2_coli(p20)

!    Z(1,1) = 2d0*q10
!    Z(2,1) = q10+q20-q21
!    Z(3,1) = q10+q30-q31
!    Z(1,2) = Z(2,1)
!    Z(2,2) = 2d0*q20
!    Z(3,2) = q20+q30-q32
!    Z(1,3) = Z(3,1)
!    Z(2,3) = Z(3,2)
!    Z(3,3) = 2d0*q30

!    q1q2 = (q10+q20-q21)
!    q1q3 = (q10+q30-q31)
!    q2q3 = (q20+q30-q32)
!    detZ = 8d0*q10*q30*q20+2D0*q1q2*q1q3*q2q3  &
!     &    -2d0*(q10*q2q3*q2q3+q20*q1q3*q1q3+q30*q1q2*q1q2)

!    Zadj(1,1) = (4d0*q30*q20-q2q3*q2q3)
!    Zadj(2,1) = (q1q3*q2q3-2d0*q30*q1q2)
!    Zadj(3,1) = (q1q2*q2q3-2d0*q20*q1q3)
!    Zadj(1,2) = Zadj(2,1)
!    Zadj(2,2) = (4d0*q10*q30-q1q3*q1q3)
!    Zadj(3,2) = (q1q2*q1q3-2d0*q10*q2q3)
!    Zadj(1,3) = Zadj(3,1)
!    Zadj(2,3) = Zadj(3,2)
!    Zadj(3,3) = (4d0*q10*q20-q1q2*q1q2)
!
!    f(1) = q10+mm02-mm12
!    f(2) = q20+mm02-mm22
!    f(3) = q30+mm02-mm32
      
!    Zadjf(1) = Zadj(1,1)*f(1)+Zadj(2,1)*f(2)+Zadj(3,1)*f(3)
!    Zadjf(2) = Zadj(1,2)*f(1)+Zadj(2,2)*f(2)+Zadj(3,2)*f(3)
!    Zadjf(3) = Zadj(1,3)*f(1)+Zadj(2,3)*f(2)+Zadj(3,3)*f(3)


    ! coefficients Shat defined in (5.13)
    allocate(Shat(0:rmaxB,0:rmaxB,0:rmaxB,2))

    do r=0,rmaxB
      do n0=0,r/2

        do n1=0,r-2*n0
          n2 = r-2*n0-n1
          Shat(n0,n1,n2,:) = -B_0(n0,n1,n2)
        end do

        k = r-2*n0
        Shat(n0,0,k,1) = Shat(n0,0,k,1) + B_i(n0,k,1)
        Shat(n0,k,0,2) = Shat(n0,k,0,2) + B_i(n0,k,2)

      end do
    end do


    ! choose reduction formulas with biggest denominators
    if (abs(Zadjf(1)).ge.abs(Zadjf(2))) then
      j = 1
    else 
      j = 2
    end if

    Zadj2f(1,2,1) = -f(2)
    Zadj2f(1,2,2) =  f(1)

    maxZadj2f = 0d0                ! Zadj2f(k,n,l) = Zadf2(k,n,l,m)*f(m)
                                   ! Zadj2(m) ==  Zadf2(k,n,l,m)
                                   ! maxZadj2f = fmax!!
    if (abs(Zadj2f(1,2,1)).gt.maxZadj2f) then
      maxZadj2f = abs(Zadj2f(1,2,1))
      k = 1
      n = 2
      l = 1
      m = 2
      Zadj2(2,2) = -1d0
    end if
    if (abs(Zadj2f(1,2,2)).gt.maxZadj2f) then
      maxZadj2f = abs(Zadj2f(1,2,2))
      k = 1
      n = 2
      l = 2
      m = 1
      Zadj2(2,1) =  1d0
    end if



    Zadjfj = Zadjf(j)
    Zadjkl = Zadj(k,l)



    Zadjfj = Zadjf(j)
    Zadjkl = Zadj(k,l)

    ! allocation of array for expanded C-coefficients
    rmaxExp = rmaxB
    allocate(Cexpgr(0:rmaxExp/2,0:rmaxExp,0:rmaxExp,0:ordgr_max))

    ! calculate Cuv
    allocate(CuvExpgr(0:(rmaxExp+1),0:rmaxExp+1,0:rmaxExp+1))
    call CalcCuv(CuvExpgr,Buv_0,mm02,f,rmaxExp+1,id)
    Cuv(0:rmax,0:rmax,0:rmax) = CuvExpgr(0:rmax,0:rmax,0:rmax)

    ! allocate arrays for error propagation
    allocate(C00_err(0:rmaxExp))
    allocate(Cij_err(0:rmaxExp))
    allocate(C00_err2(0:rmaxExp))
    allocate(Cij_err2(0:rmaxExp))
     
    ! initialize accuracy estimates
    Cerr = acc_inf
    Cij_err =0d0
    C00_err =0d0

    Cerr2 = acc_inf
    Cij_err2 =0d0
    C00_err2 =0d0

!    maxZadj = maxval(abs(Zadj))
!    maxZadj2f = maxval(abs(f(inds2(1,:))*Zadj2(:)))

    ! truncation of expansion if calculated term larger than truncfacexp * previous term
    ! crucial for expansion parameters between 0.1 and 1 !!!
    truncfacexp = sqrt(fac_gr) * truncfacC
    gtrunc = ordgr_max 

! calculate C(n0,n1,n2) up to rank r+n0
    rloop: do r=0,rmaxExp/2



      if (r.gt.rmax+gtrunc) exit rloop



      !!!!!!!!!!!!!!!!!!!!!!!!!!!!
      ! 0th-order coefficients
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!

      ! calculating
      ! C_00(a)0000..00 --> C_00(a)ij00..00 --> C_00(a)ijkl00..00 --> ... --> C_00(a)ijklmn..
      ! exploiting eq. (5.40) - (5.53) solved for C_00i1..<ir>...iP
      maxCexpgr(1,r,0)=0d0
      do n0=r,1,-1
        do nn=r-n0,0,-1
          nntt = r-n0-nn



          inds0(n) = nn
          inds0(k) = nntt
            


          inds1(n) = nn+1
          inds1(k) = nntt



          Caux = -Zadj(k,l)*B_0(n0-1,inds1(1),inds1(2))

!            Caux = 2*Zadj(k,l) * (1+r-2*n0) * Cexpgr(n0,inds1(1),inds1(2),0)

!            inds = inds1
!            inds(k) = inds(k) + 1
!            inds(l) = inds(l) + 1
!            Caux = Caux + detZ * Cexpgr(n0-1,inds(1),inds(2),0)
!
!            inds = inds1
!            inds(k) = inds(k) + 1
!            Caux = Caux + Zadjf(l) * Cexpgr(n0-1,inds(1),inds(2),0)


          
          inds = inds1
          inds(k) = inds(k)+1
          do i=1,2
            Caux = Caux - Zadj(i,l)*Shat(n0-1,inds(1),inds(2),i)

          end do



          do i=1,2
            inds = inds1
            inds(i) = inds(i)+1
            Caux = Caux + Zadj(k,l)*Shat(n0-1,inds(1),inds(2),i)

          end do




          Caux = Caux + 2*(nn+1) *Zadj2(n ,m )*Shat(n0,inds0(1),inds0(2),m)




!            Caux = Caux - 2*(nn+1)* Zadj2f(k,n,l)*Cexpgr(n0,inds0(1),inds0(2),0) 

          inds = inds1
          if(m.eq.n) then
            if (inds(n).gt.1) then
              inds(n) = inds(n)-2
              Caux = Caux - 4*(nn+1)*nn * Zadj2(n,m ) * Cexpgr(n0+1,inds(1),inds(2),0) 

            end if
          else
            if (inds(n).gt.0.and.inds(m).gt.0) then
              inds(n) = inds(n)-1
              inds(m) = inds(m)-1
              Caux = Caux - 4*(nn+1)*(inds(m)+1)* Zadj2(n,m ) * Cexpgr(n0+1,inds(1),inds(2),0) 

            end if
          end if
          
          Cexpgr(n0,inds0(1),inds0(2),0) = Caux/(2*(nn+1)* Zadj2f(k,n,l))
          
          if (n0.eq.1) then
            maxCexpgr(1,r,0) =  maxCexpgr(1,r,0) + abs(Cexpgr(n0,inds0(1),inds0(2),0) )
          end if
          
!          if (r+n0.le.rmax) then             !  for fixed rank
          if (r.le.rmax) then
            C(n0,inds0(1),inds0(2)) = Cexpgr(n0,inds0(1),inds0(2),0)
          end if
          
        end do
      end do

      ! calculate
      ! C_00ijkl.. --> C_aijkl..
      ! exploiting eq. (5.38)
      maxCexpgr(0,r,0)=0d0
      do n1=0,r
        n2 = r-n1

        Smod = Shat(0,n1,n2,:)
        if (n1.ge.1) then
          Smod(1) = Smod(1) - 2d0*n1*Cexpgr(1,n1-1,n2,0)
        end if
        if (n2.ge.1) then
          Smod(2) = Smod(2) - 2d0*n2*Cexpgr(1,n1,n2-1,0)
        end if

        Cexpgr(0,n1,n2,0) = (Zadj(1,j)*Smod(1) +  Zadj(2,j)*Smod(2)  &
            )/Zadjfj
        maxCexpgr(0,r,0) =  maxCexpgr(0,r,0) + abs(Cexpgr(0,n1,n2,0))
        if (r.le.rmax) then
          C(0,n1,n2) = Cexpgr(0,n1,n2,0)
        end if




      end do



      if(r.le.rmax) then
!       Cerr(r) =  abs(detZ/Zadjfj)*maxCexpgr(0,r,0)
        Cerr(r) =  fac_gr*maxCexpgr(0,r,0)
      end if

      ! error propagation from C's
      if(r.gt.0)then
        C00_err(r+1) = maxZadj*B_err/(2*maxZadj2f)
      end if
      Cij_err(r)=maxZadj*max(B_err,2*C00_err(r+1))/abs(Zadjfj)

      if(r.gt.0)then
        C00_err2(r+1) = maxZadj*B_err/(2*maxZadj2f)
      end if
      Cij_err2(r)=maxZadj*max(B_err,2*C00_err2(r+1))/abs(Zadjfj)

      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      ! higher order coefficients
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      rg = r
      gloop: do g=1,min(gtrunc,r)
        rg = rg-1



        ! calculating
        ! C_00(a)0000..00 --> C_00(a)ij00..00 --> C_00(a)ijkl00..00 --> ... --> C_00(a)ijklmn..
       ! exploiting eq. (5.40) - (5.53) solved for C_00i1..<ir>...iP
        maxCexpgr(1,rg,g) = 0d0
        do n0=rg,1,-1                      !  note rank of tensor = rg+n0
          do nn=rg-n0,0,-1
            nntt = rg-n0-nn
            inds0(n) = nn
            inds0(k) = nntt
              
            inds1(n) = nn+1
            inds1(k) = nntt
              


            Caux = 2*Zadj(k,l) * (2+rg-n0) * Cexpgr(n0,inds1(1),inds1(2),g-1)
              


            if (g.gt.1) then
              inds = inds1
              inds(k) = inds(k) + 1
              inds(l) = inds(l) + 1
              Caux = Caux + detZ * Cexpgr(n0-1,inds(1),inds(2),g-2)
                

            end if

            inds = inds1
            inds(k) = inds(k) + 1
            Caux = Caux + Zadjf(l) * Cexpgr(n0-1,inds(1),inds(2),g-1)
            


!            Caux = Caux - 2*nn* Zadj2f(k,n,l)*Cexpgr(n0,inds0(1),inds0(2),g) 

            inds = inds1
            if(m.eq.n) then
              if (inds(n).gt.1) then
                inds(n) = inds(n)-2
                Caux = Caux - 4*(nn+1)*nn * Zadj2(n,m ) * Cexpgr(n0+1,inds(1),inds(2),g) 

              end if
            else
              if (inds(n).gt.0.and.inds(m).gt.0) then
                inds(n) = inds(n)-1
                inds(m) = inds(m)-1
                Caux = Caux - 4*(nn+1)*(inds(m)+1)* Zadj2(n,m ) * Cexpgr(n0+1,inds(1),inds(2),g) 

              end if
            end if
              

            Cexpgr(n0,inds0(1),inds0(2),g) = Caux/(2*(nn+1)* Zadj2f(k,n,l))
                             
            if(n0.eq.1) then
              maxCexpgr(1,rg,g) =  maxCexpgr(1,rg,g) + abs(Cexpgr(n0,inds0(1),inds0(2),g))
              
              if (g.eq.1.and.abs(Cexpgr(n0,inds0(1),inds0(2),g)).gt.   & 
                  truncfacexp*max(1d0,maxCexpgr(1,rg,g-1)) .or.        &
                  g.ge.2.and.abs(Cexpgr(n0,inds0(1),inds0(2),g)).gt.   &
                  truncfacexp*maxCexpgr(1,rg,g-1)) then


                
                gtrunc = g-1
                exit gloop
              end if
            end if

          end do
        end do


        if (rg.le.rmax) then
          do n0=rg,1,-1
!            if (rg+n0.le.rmax) then       ! for fixed rank
            if (rg.le.rmax) then
              do n1=0,rg-n0
                n2=rg-n0-n1
                C(n0,n1,n2) = C(n0,n1,n2) + Cexpgr(n0,n1,n2,g)
              end do
            end if
          end do
        end if

!        write(*,*) 'CalcCgr after it1 ',rg

        ! calculate
        ! C_00ijkl.. --> C_aijkl..
        ! exploiting eq. (5.38)

!        write(*,*) 'CalcCgr maxCexp',rg,g-1,maxCexpgr(0,rg,g-1)

        maxCexpgr(0,rg,g) = 0d0
        do n1=0,rg
          n2 = rg-n1
          
          Smod = 0d0
          if (n1.ge.1) then
            Smod(1) = Smod(1) - 2d0*n1*Cexpgr(1,n1-1,n2,g)
          end if
          if (n2.ge.1) then
            Smod(2) = Smod(2) - 2d0*n2*Cexpgr(1,n1,n2-1,g)
          end if

          inds(1) = n1
          inds(2) = n2
          inds(j) = inds(j)+1
          Cexpgr(0,n1,n2,g) = (Zadj(1,j)*Smod(1) +  Zadj(2,j)*Smod(2)  &
              - detZ*Cexpgr(0,inds(1),inds(2),g-1))/Zadjfj

          maxCexpgr(0,rg,g) =  maxCexpgr(0,rg,g) + abs(Cexpgr(0,n1,n2,g))

!              if(n1.eq.0.and.n2.eq.1) then 
!                write(*,*) 'C2(2,3)= ',g,Cexpgr(0,n1,n2,g)
!                write(*,*) 'C2(2,3)= ',Zadj(1,j)*Smod(1)/Zadjfj,  Zadj(2,j)*Smod(2)/Zadjfj,  &
!                                - detZ*Cexpgr(0,inds(1),inds(2),g-1)/Zadjfj
!                write(*,*) 'C2(2,3)= ',inds(1),inds(2),         &
!                                - detZ/Zadjfj,Cexpgr(0,inds(1),inds(2),g-1)
!              end if

          if (g.eq.1.and.abs(Cexpgr(0,n1,n2,g)).gt.        &
              truncfacexp*max(1d0/m2scale,maxCexpgr(0,rg,g-1)).or.     &
              g.ge.2.and.abs(Cexpgr(0,n1,n2,g)).gt.        &
              truncfacexp*maxCexpgr(0,rg,g-1)) then


            gtrunc = g-1
            exit gloop
          end if

        end do

        ! error propagation from C's
        if(rg.gt.0)then
          C00_err(rg+1) = max( C00_err(rg+1),                   &
              max( maxZadj*(2+rg-2*n0)*C00_err(rg+2),       &
                    abs(detZ)*Cij_err(rg+2),                &
                    maxZadjf*Cij_err(rg+1)                   &
                 ) / (2*maxZadj2f)  ) 
        end if
        Cij_err(rg)=max(Cij_err(rg),                &
            max(2*maxZadj*C00_err(rg+1),abs(detZ)*Cij_err(rg))/abs(Zadjfj) )

        if(rg.gt.0)then
          C00_err2(rg+1) = max( C00_err2(rg+1),                   &
              max( maxZadj*(2+rg-2*n0)*C00_err2(rg+2),       &
                    abs(detZ)*Cij_err2(rg+2),                &
                    maxZadjf*Cij_err2(rg+1)                   &
                 ) / (2*maxZadj2f)  ) 
        end if
        Cij_err2(rg)=max(Cij_err2(rg),                &
            max(2*maxZadj*C00_err2(rg+1),abs(detZ)*Cij_err2(rg))/abs(Zadjfj) )



        if (rg.le.rmax) then
          Cerr(rg) = 0d0
          do n1=0,rg
              n2 = rg-n1
              C(0,n1,n2) = C(0,n1,n2) + Cexpgr(0,n1,n2,g)
            if(abs(Cexpgr(0,n1,n2,g-1)).ne.0d0) then
!             Cerr(rg)=max(Cerr(rg),abs(Cexpgr(0,n1,n2,g))**2/abs(Cexpgr(0,n1,n2,g-1)))
              Cerr(rg)=max(Cerr(rg),abs(Cexpgr(0,n1,n2,g))*min(1d0,abs(Cexpgr(0,n1,n2,g))/abs(Cexpgr(0,n1,n2,g-1))))
            else
              Cerr(rg)=max(Cerr(rg),abs(Cexpgr(0,n1,n2,g)))
            end if



          end do

          ! if error from B's larger than error from expansion stop expansion
          if(Cij_err(rg).gt.3d0*Cerr(rg)) then
            gtrunc = min(g,gtrunc)
             


          end if

        end if

      end do gloop



      Cerr2 = max(Cerr,Cij_err2(0:rmax))
      Cerr = max(Cerr,Cij_err(0:rmax))



!      if(maxval(Cerr).le.acc_req_C*abs(C(0,0,0))) exit           changed 28.01.15
      ! check if target precision already reached

      if(maxval(Cerr-acc_req_Cr*abs(C(0,0,0))).le.0d0) then
        do rg=r+1,rmax
        do n0=0,rg/2
        do n1=0,rg-n0
          C(n0,n1,rg-n0-n1)=0d0
        end do
        end do
        end do

        exit rloop
      end if

    end do rloop


                           
!   write(*,*) 'CalcCgr out',(((C((r-n1-n2)/2,n1,n2),n2=0,r-n1),n1=0,r),r=0,rmax)


!    write(*,*) 'CalcCgr Cerr ',Cerr
!    write(*,*) 'CalcCgr Cerr2',Cerr2

  end subroutine CalcCgr




  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !  subroutine CalcCgy(C,Cuv,p10,p21,p20,m02,m12,m22,rmax,ordgy_min,ordgy_max,id,Cerr,acc_req_Cr,Cerr2)
  !
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  ! modified version of Ansgar  (similar to CalcDgy)
  
  subroutine CalcCgy(C,Cuv,p10,p21,p20,m02,m12,m22,rmax,ordgy_min,ordgy_max,id,Cerr,acc_req_Cr,Cerr2)
  
    use globalC
  
    integer, intent(in) :: rmax,ordgy_min,ordgy_max,id
    double complex, intent(in) :: p10,p21,p20,m02,m12,m22
    double complex, intent(out) :: C(0:rmax,0:rmax,0:rmax)
    double complex, intent(out) :: Cuv(0:rmax,0:rmax,0:rmax)
    double precision, intent(out) :: Cerr(0:rmax),Cerr2(0:rmax)
    double precision, intent(in) :: acc_req_Cr(0:rmax)
    double complex, allocatable :: Cexpgy(:,:,:,:), CuvExpgy(:,:,:)
    double complex, allocatable :: B_0(:,:,:), B_i(:,:,:), Shat(:,:,:,:)
    double complex, allocatable :: Buv_0(:,:,:), Buv_i(:,:,:)
    double complex :: Smod, Caux, Zadj2f
    double complex :: C0_coli, elimminf2_coli
    double precision, allocatable :: C00_err(:),Cij_err(:)
    double precision, allocatable :: C00_err2(:),Cij_err2(:)
    double precision :: B_err,B_max,aZadj2f
    double precision :: maxCexpgy(0:1,0:rmax+2*ordgy_min,0:ordgy_max),truncfacexp
    integer :: rmaxB,rmaxExp,gtrunc,r,n0,n1,n2,i,j,jt,g,rg
    integer :: inds0(2),inds(2),k,l,lt,nl,nlt
    integer :: bin,nid(0:2)



    ! write(*,*) 'LH: CalcCgy, ord', ordgy_min
    ! calculation of B-coefficients
    rmaxB = rmax + 2*ordgy_min + 1
    allocate(B_0(0:rmaxB,0:rmaxB,0:rmaxB))
    allocate(Buv_0(0:rmaxB,0:rmaxB,0:rmaxB))
    allocate(B_i(0:rmaxB,0:rmaxB,2))
    allocate(Buv_i(0:rmaxB,0:rmaxB,2))

    ! determine binaries for B-coefficients
    k=0
    bin = 1
    do while (k.le.2)
      if (mod(id/bin,2).eq.0) then
        nid(k) = id+bin
        k = k+1
      end if
      bin = 2*bin
    end do

    call CalcB(B_0(:,0,:),Buv_0(:,0,:),p21,m12,m22,rmaxB,nid(0))
    call CalcB(B_i(:,:,1),Buv_i(:,:,1),p20,m02,m22,rmaxB,nid(1))
    call CalcB(B_i(:,:,2),Buv_i(:,:,2),p10,m02,m12,rmaxB,nid(2))

   ! shift of integration momentum in B_0
    B_max=0d0
    do n1=1,rmaxB
      do n2=0,rmaxB-n1
        n0 = (rmaxB-n1-n2)
        B_0(0:n0,n1,n2) = -B_0(0:n0,n1-1,n2)-B_0(0:n0,n1-1,n2+1)
        Buv_0(0:n0,n1,n2) = -Buv_0(0:n0,n1-1,n2)-Buv_0(0:n0,n1-1,n2+1)
      end do
    end do
    B_max = max(B_max,maxval(abs(B_i(0,0:rmaxB,1:2))))
    B_err = acc_def_B*B_max

    ! determine (adjugated) Gram and Cayley matrix
!    mm02 = elimminf2_coli(m02)
!    mm12 = elimminf2_coli(m12)
!    mm22 = elimminf2_coli(m22)
!    q10  = elimminf2_coli(p10)
!    q21  = elimminf2_coli(p21)
!    q20  = elimminf2_coli(p20)
!
!    q1q2 = (q10+q20-q21)
!    detZ = 4d0*q10*q20-q1q2*q1q2

    if (abs(detZ).lt.abs(4d0*q10*q20 + Z(2,1)*Z(2,1))*1d-4) then
      if (abs(q10-q20).lt.abs(q10-q21).and.  &
          abs(q10-q20).lt.abs(q20-q21)) then
        detZ  =  4d0*q10*q21 - (q10-q20+q21)*(q10-q20+q21)
      end if
    end if

!    Zadj(1,1) = 2d0*q20
!    Zadj(2,1) = -q1q2
!    Zadj(1,2) = -q1q2
!    Zadj(2,2) = 2d0*q10
!    f(1) = q10+mm02-mm12
!    f(2) = q20+mm02-mm22
!
!    Zadjf(1) = Zadj(1,1)*f(1)+Zadj(2,1)*f(2)
!    Zadjf(2) = Zadj(1,2)*f(1)+Zadj(2,2)*f(2)
!
! Xadj(1,1) and Xadj(2,2) exchanged!!!
!    Xadj(1,1) = 2d0*mm02*Z(1,1) - f(1)*f(1)
!    Xadj(2,1) = 2d0*mm02*Z(1,2) - f(1)*f(2)
!    Xadj(1,2) = Xadj(2,1)
!    Xadj(2,2) = 2d0*mm02*Z(2,2) - f(2)*f(2)


    ! coefficients Shat defined in (5.13)
    allocate(Shat(0:rmaxB,0:rmaxB,0:rmaxB,2))

    do r=0,rmaxB
      do n0=0,r/2

        do n1=0,r-2*n0
          n2 = r-2*n0-n1
          Shat(n0,n1,n2,:) = -B_0(n0,n1,n2)
        end do

        k = r-2*n0
        Shat(n0,0,k,1) = Shat(n0,0,k,1) + B_i(n0,k,1)
        Shat(n0,k,0,2) = Shat(n0,k,0,2) + B_i(n0,k,2)

      end do
    end do

    ! choose reduction formulas with biggest denominators
    if (abs(Xadj(1,1)).ge.abs(Xadj(2,2))) then
      if (abs(Xadj(1,1)).ge.abs(Xadj(1,2))) then
        i = 1
        j = 1
        jt = 2
        Zadj2f = -f(2)  
      else
        i = 1
        j = 2
        jt = 1
        Zadj2f = f(2)
      end if
    else        
      if (abs(Xadj(2,2)).ge.abs(Xadj(1,2))) then
        i = 2
        j = 2
        jt = 1
        Zadj2f = -f(1)
      else 
        i = 1
        j = 2
        jt = 2
        Zadj2f = -f(2) 
      end if
    end if
    aZadj2f = abs(Zadj2f)

    if (abs(Zadj(1,1)).ge.abs(Zadj(2,2))) then
      if (abs(Zadj(1,1)).ge.abs(Zadj(1,2))) then
        k = 1
        l = 1
        lt = 2
      else
        k = 1
        l = 2
        lt = 1
      end if
    else        
      if (abs(Zadj(2,2)).ge.abs(Zadj(1,2))) then
        k = 2
        l = 2
        lt = 1
      else 
        k = 1
        l = 2
        lt = 1
      end if
    end if



!   write(*,*)  'CalcCgy Zadj(i,j)=',i,j,Zadj(i,j),Xadj(i,j)
        
    ! allocation of array for det(Z)- and det(X)-expanded C-coefficients
    rmaxExp = rmaxB+1
    allocate(Cexpgy(0:max(rmax/2,1),0:rmaxExp-2,0:rmaxExp-2,0:ordgy_max))
   
    ! calculate Cuv
    allocate(CuvExpgy(0:rmaxExp,0:rmaxExp,0:rmaxExp))
    call CalcCuv(CuvExpgy,Buv_0,mm02,f,rmaxExp,id)
    Cuv(0:rmax,0:rmax,0:rmax) = CuvExpgy(0:rmax,0:rmax,0:rmax)

    ! allocate arrays for error propagation
    allocate(C00_err(0:rmaxExp))
    allocate(Cij_err(0:rmaxExp))
    allocate(C00_err2(0:rmaxExp))
    allocate(Cij_err2(0:rmaxExp))

    ! initialize accuracy estimates
    Cerr = acc_inf
    Cij_err =0d0
    C00_err =0d0

    Cerr2 = acc_inf
    Cij_err2 =0d0
    C00_err2 =0d0

!    maxZadjf = maxval(abs(Zadjf))
!    fmax = maxval(abs(f))

    ! truncation of expansion if calculated term larger than truncfacexp * previous term
    ! crucial for expansion parameters between 0.1 and 1 !!!
!    truncfacexp = sqrt(max(maxZadjf,abs(detZ))/abs(Xadj(i,j))*max(1d0,fmax/abs(Zadj(k,l)))) * truncfacC
    truncfacexp = sqrt(fac_gy) * truncfacC
    gtrunc = ordgy_max 



! calculate C(1,n1,n2) up to rank r+2
! calculate C(0,n1,n2) up to rank r  
    rloop: do r=0,rmaxExp-2



      if (r.gt.rmax+2*gtrunc+2) exit rloop

      !!!!!!!!!!!!!!!!!!!!!!!!!!!!
      ! 0th-order coefficients
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!

      ! calculating C_00ijk.. exploiting eq. (5.49)
      maxCexpgy(1,r,0)=0d0
      do nl=r,0,-1
        nlt=r-nl
        inds0(l) = nl
        inds0(lt) = nlt

        inds(l) = nl+1
        inds(lt) = nlt
        Caux = Zadj(k,1)*Shat(0,inds(1),inds(2),1)  &
             + Zadj(k,2)*Shat(0,inds(1),inds(2),2)

        if (nlt.ge.1) then
          inds(lt) = nlt-1
          Caux = Caux - 2*nlt*Zadj(k,lt)*Cexpgy(1,inds(1),inds(2),0)
        end if

        Cexpgy(1,inds0(1),inds0(2),0) = Caux/(2*(nl+1)*Zadj(k,l))
        maxCexpgy(1,r,0) =  maxCexpgy(1,r,0) + abs(Cexpgy(1,inds0(1),inds0(2),0) )
!        if (r+2.le.rmax) then             !  for fixed rank
        if (r+1.le.rmax) then
          C(1,inds0(1),inds0(2)) = Cexpgy(1,inds0(1),inds0(2),0)
        end if

      end do

      ! calculate C_ijkl.. exploiting eq. (5.53)
      maxCexpgy(0,r,0)=0d0
      do n1=0,r
        n2 = r-n1
        inds(1) = n1
        inds(2) = n2
      
        Caux = (2*(2+r)*Cexpgy(1,n1,n2,0) - 4*CuvExpgy(1,n1,n2)  &
             - B_0(0,n1,n2))*Zadj(i,j)

!     write(*,*) 'CalcCred Caux',Caux,Zadj(i,j),f(i),f(j)

        Smod = Shat(0,n1,n2,jt)

        if (inds(jt).ge.1) then
          inds(jt) = inds(jt)-1
          Smod = Smod - 2d0*(inds(jt)+1)*Cexpgy(1,inds(1),inds(2),0)
        end if

        Caux = Caux + Zadj2f*Smod

!       write(*,*) 'CalcCgy maxadjf',maxZadjf,Xadj(i,j),Caux

        Cexpgy(0,n1,n2,0) = Caux/Xadj(i,j)
        maxCexpgy(0,r,0) =  maxCexpgy(0,r,0) + abs(Cexpgy(0,n1,n2,0))
        if (r.le.rmax) then
          C(0,n1,n2) = Cexpgy(0,n1,n2,0)
        end if

      end do

      if (r.le.rmax) then
!       Cerr(r) =  abs(maxZadjf/Xadj(i,j))*maxCexpgy(0,r,0)
        Cerr(r) =  fac_gy*maxCexpgy(0,r,0)

!        write(*,*) 'CalcCgy Cerr,0 ',r,Cerr(r),fac_gy,maxCexpgy(0,r,0)

      end if

 ! error propagation from B's
      C00_err(r+2) = B_err /2d0               
      Cij_err(r)=max(abs(Zadj(i,j))/abs(Xadj(i,j))*max(B_err,2*(r+2)*C00_err(r+2)),  &
                  fmax/abs(Xadj(i,j))*max(B_err,2*C00_err(r+1)))

      C00_err2(r+2) = B_err /2d0               
      Cij_err2(r)=max(abs(Zadj(i,j))/abs(Xadj(i,j))*max(B_err,2*(r+2)*C00_err2(r+2)),  &
                  fmax/abs(Xadj(i,j))*max(B_err,2*C00_err2(r+1)))

      
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      ! higher order coefficients
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      rg = r
      gloop: do g=1,min(gtrunc,r/2)
        rg = rg-2

!       write(*,*) 'CalcCgy gtrunc gloop=',gtrunc,r,g,rg

        ! calculating C_00ijk.. exploiting eq. (5.49)
        maxCexpgy(1,rg,g) = 0d0
        do nl=rg,0,-1
          nlt=rg-nl
          inds0(l) = nl
          inds0(lt) = nlt

          inds(l) = nl+1
          inds(lt) = nlt
          Caux = -Zadjf(k)*Cexpgy(0,inds(1),inds(2),g-1)

          inds(k) = inds(k)+1
          Caux = Caux - detZ*Cexpgy(0,inds(1),inds(2),g-1)

          if (nlt.ge.1) then
            inds(l) = nl+1
            inds(lt) = nlt-1
            Caux = Caux - 2*nlt*Zadj(k,lt)*Cexpgy(1,inds(1),inds(2),g)
          end if

          Cexpgy(1,inds0(1),inds0(2),g) = Caux/(2*(nl+1)*Zadj(k,l))
          maxCexpgy(1,rg,g) =  maxCexpgy(1,rg,g) + abs(Cexpgy(1,inds0(1),inds0(2),g) )

          if (g.eq.1.and.abs(Cexpgy(1,inds0(1),inds0(2),g)).gt.      &
                  truncfacexp*max(1d0,maxCexpgy(1,rg,g-1))    .or.   &
                  g.ge.2.and.abs(Cexpgy(1,inds0(1),inds0(2),g)).gt.  &
                  truncfacexp*maxCexpgy(1,rg,g-1)) then


            gtrunc = g-1
            exit gloop

          end if 
        
        end do


!        if (rg+2.le.rmax) then            !  for fixed rank
        if (rg+1.le.rmax) then
          do nl=rg,0,-1
            nlt=rg-nl
            inds0(l) = nl
            inds0(lt) = nlt
            C(1,inds0(1),inds0(2)) = C(1,inds0(1),inds0(2))  &
                + Cexpgy(1,inds0(1),inds0(2),g)
          end do
        end if


        ! calculate C_ijkl.. exploiting eq. (5.53)
        maxCexpgy(0,rg,g) = 0d0
        do n1=0,rg
          n2 = rg-n1
          inds0(1) = n1
          inds0(2) = n2
          
          Caux = 2*(2+rg)*Cexpgy(1,n1,n2,g)*Zadj(i,j)
          
!          write(*,*) 'CalcCgy g Caux 1',rg,g,Caux

          if (inds0(jt).ge.1) then
            inds = inds0
            inds(jt) = inds(jt)-1
            Caux = Caux - 2d0*Zadj2f*inds0(jt)*Cexpgy(1,inds(1),inds(2),g)
          end if
          
!          write(*,*) 'CalcCgy g Caux 2',rg,g,Caux

          inds0(i) = inds0(i)+1
          Caux = Caux - Zadjf(j)*Cexpgy(0,inds0(1),inds0(2),g-1)
          
!          write(*,*) 'CalcCgy g Caux 3',rg,g,Caux

          Cexpgy(0,n1,n2,g) = Caux/Xadj(i,j)

!          write(*,*) 'CalcCgy g Cexpgy',rg,g,n1,n2,Cexpgy(0,n1,n2,g)

          maxCexpgy(0,rg,g) =  maxCexpgy(0,rg,g) + abs(Cexpgy(0,n1,n2,g))
          
          if (g.eq.1.and.abs(Cexpgy(0,n1,n2,g)).gt.        &
              truncfacexp*max(1d0/m2scale,maxCexpgy(0,rg,g-1)).or.     &
              g.ge.2.and.abs(Cexpgy(0,n1,n2,g)).gt.        &
              truncfacexp*maxCexpgy(0,rg,g-1)) then



            gtrunc = g-1
            exit gloop

          end if

!            if ((g.ge.2).and.(abs(Cexpgy(0,n1,n2,g)).gt.truncfacexp*abs(Cexpgy(0,n1,n2,g-1)))) then
!              gtrunc = g-1
!            end if

        end do

        ! error propagation from B's
        if(rg.gt.1)then
          C00_err(rg+2) =max(C00_err(rg+2),                                       &
              max(abs(Zadjf(k))/2d0*Cij_err(rg+1),                                &
                  abs(detZ)/2d0*Cij_err(rg+2))/abs(Zadj(k,l)))   
        end if



        Cij_err(rg)= max( Cij_err(rg),                                &
            max(2*(rg+2)*abs(Zadj(i,j))*C00_err(rg+2),                &
                2*abs(Zadj2f)*C00_err(rg+1),       &
                abs(Zadjf(j))*Cij_err(rg+1))/abs(Xadj(i,j)))

        if(rg.gt.1)then
          C00_err2(rg+2) =max(C00_err2(rg+2),                                       &
              max(abs(Zadjf(k))/2d0*Cij_err2(rg+1),                                &
                  abs(detZ)/2d0*Cij_err2(rg+2))/abs(Zadj(k,l)))   
        end if

        Cij_err2(rg)= max( Cij_err2(rg),                                &
            max(2*(rg+2)*abs(Zadj(i,j))*C00_err2(rg+2),                &
                2*abs(Zadj2f)*C00_err2(rg+1),       &
                abs(Zadjf(j))*Cij_err2(rg+1))/abs(Xadj(i,j)))



        if ((rg.le.rmax)) then
          Cerr(rg) = 0d0
          do n1=0,rg
            n2=rg-n1
            C(0,n1,n2) = C(0,n1,n2) + Cexpgy(0,n1,n2,g)



            if(abs(Cexpgy(0,n1,n2,g-1)).ne.0d0) then
              Cerr(rg)=max(Cerr(rg),abs(Cexpgy(0,n1,n2,g))*min(1d0,abs(Cexpgy(0,n1,n2,g))/abs(Cexpgy(0,n1,n2,g-1))))
            else
              Cerr(rg)=max(Cerr(rg),abs(Cexpgy(0,n1,n2,g)))
            end if



          end do

          ! if error from B's larger than error from expansion stop expansion
          if(Cij_err(rg).gt.Cerr(rg)) then
             gtrunc = min(g,gtrunc)
!            gtrunc = min(g+1,gtrunc)
             


          end if

        end if
        
        
      end do gloop

!     write(*,*) 'CalcCgy gtrunc aft gloop=',gtrunc,r



      Cerr2 = max(Cerr,Cij_err2(0:rmax))
      Cerr = max(Cerr,Cij_err(0:rmax))



      ! check if target precision already reached
!      if(maxval(Cerr-acc_req_Cr*abs(C(0,0,0))).le.0d0) exit         ! changed 28.01.15

      if(maxval(Cerr-acc_req_Cr*abs(C(0,0,0))).le.0d0) then
        do rg=r+1,rmax
          do n1=0,rg
            C(0,n1,rg-n1)=0d0
          end do
        end do   
        do rg=r+1,rmax
          do n1=0,rg-2
            C(1,n1,rg-2-n1)=0d0
          end do
        end do   




        exit rloop 

      end if

    end do rloop


    ! calculating C_0000ijk.. exploiting eq. (5.49)
    do r=4,rmax
!      do n0=2,rmax/2     !     for fixed rank
      do n0=2,rmax
        do nl=r-2*n0,0,-1
          nlt=r-2*n0-nl
          inds0(l) = nl
          inds0(lt) = nlt

          inds(l) = nl+1
          inds(lt) = nlt
          Caux = Zadj(k,1)*Shat(n0-1,inds(1),inds(2),1)  &
               + Zadj(k,2)*Shat(n0-1,inds(1),inds(2),2)  &
               - Zadjf(k)*C(n0-1,inds(1),inds(2))

          inds(k) = inds(k)+1
          Caux = Caux - detZ*C(n0-1,inds(1),inds(2))

          if (nlt.ge.1) then
            inds(l) = nl+1
            inds(lt) = nlt-1
            Caux = Caux - 2*nlt*Zadj(k,lt)*C(n0,inds(1),inds(2))
          end if

          C(n0,inds0(1),inds0(2)) = Caux/(2*(nl+1)*Zadj(k,l))

        end do
      end do
    end do

      ! reduction formula (5.10) for n0+n1+n2=r, n0>0 
    do r=rmax+1,2*rmax
      do n0=r-rmax,r/2
        do n1=0,r-2*n0
          n2 = r-2*n0-n1
          C(n0,n1,n2) = (B_0(n0-1,n1,n2) + 2*mm02*C(n0-1,n1,n2) + 4*Cuv(n0,n1,n2) &
                        + f(1)*C(n0-1,n1+1,n2) + f(2)*C(n0-1,n1,n2+1)) / (2*r)
        end do
      end do
    end do



!   write(*,*) 'CalcCgy out',(((C((r-n1-n2)/2,n1,n2),n2=0,r-n1),n1=0,r),r=0,rmax)


!    write(*,*) 'CalcCgy Cerr ',Cerr
!    write(*,*) 'CalcCgy Cerr2',Cerr2

  end subroutine CalcCgy





  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !  subroutine CalcCgyo(C,Cuv,p10,p21,p20,m02,m12,m22,rmax,ordgy_min,ordgy_max,id,Cerr,acc_req_Cr,Cerr2)
  !
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  ! version of Lars
  
  subroutine CalcCgyo(C,Cuv,p10,p21,p20,m02,m12,m22,rmax,ordgy_min,ordgy_max,id,Cerr,acc_req_Cr,Cerr2)
  
    use globalC
  
    integer, intent(in) :: rmax,ordgy_min,ordgy_max,id
    double complex, intent(in) :: p10,p21,p20,m02,m12,m22
    double complex, intent(out) :: C(0:rmax,0:rmax,0:rmax)
    double complex, intent(out) :: Cuv(0:rmax,0:rmax,0:rmax)
    double precision, intent(out) :: Cerr(0:rmax),Cerr2(0:rmax)
    double precision, intent(in) :: acc_req_Cr(0:rmax)
    double complex, allocatable :: Cexpgy(:,:,:,:), CuvExpgy(:,:,:)
    double complex, allocatable :: B_0(:,:,:), B_i(:,:,:), Shat(:,:,:,:)
    double complex, allocatable :: Buv_0(:,:,:), Buv_i(:,:,:)
    double complex :: Smod, Caux
    double complex :: C0_coli, elimminf2_coli
    double precision, allocatable :: C00_err(:),Cij_err(:)
    double precision, allocatable :: C00_err2(:),Cij_err2(:)
    double precision :: B_err,B_max
    double precision :: maxCexpgy(0:1,0:rmax+2*ordgy_min,0:ordgy_max),truncfacexp
    integer :: rmaxB,rmaxExp,gtrunc,r,n0,n1,n2,a,b,j,sgnab,g,rg
    integer :: inds0(2),inds(2),at,bt,k,l,lt,nl,nlt
    integer :: bin,nid(0:2)



    ! write(*,*) 'LH: CalcCgy, ord', ordgy_min
    ! calculation of B-coefficients
    rmaxB = rmax + 2*ordgy_min + 1
    allocate(B_0(0:rmaxB,0:rmaxB,0:rmaxB))
    allocate(Buv_0(0:rmaxB,0:rmaxB,0:rmaxB))
    allocate(B_i(0:rmaxB,0:rmaxB,2))
    allocate(Buv_i(0:rmaxB,0:rmaxB,2))

    ! determine binaries for B-coefficients
    k=0
    bin = 1
    do while (k.le.2)
      if (mod(id/bin,2).eq.0) then
        nid(k) = id+bin
        k = k+1
      end if
      bin = 2*bin
    end do

    call CalcB(B_0(:,0,:),Buv_0(:,0,:),p21,m12,m22,rmaxB,nid(0))
    call CalcB(B_i(:,:,1),Buv_i(:,:,1),p20,m02,m22,rmaxB,nid(1))
    call CalcB(B_i(:,:,2),Buv_i(:,:,2),p10,m02,m12,rmaxB,nid(2))

   ! shift of integration momentum in B_0
    B_max=0d0
    do n1=1,rmaxB
      do n2=0,rmaxB-n1
        n0 = (rmaxB-n1-n2)
        B_0(0:n0,n1,n2) = -B_0(0:n0,n1-1,n2)-B_0(0:n0,n1-1,n2+1)
        Buv_0(0:n0,n1,n2) = -Buv_0(0:n0,n1-1,n2)-Buv_0(0:n0,n1-1,n2+1)
      end do
    end do
    B_max = max(B_max,maxval(abs(B_i(0,0:rmaxB,1:2))))
    B_err = acc_def_B*B_max

    ! determine (adjugated) Gram and Cayley matrix
!    mm02 = elimminf2_coli(m02)
!    mm12 = elimminf2_coli(m12)
!    mm22 = elimminf2_coli(m22)
!    q10  = elimminf2_coli(p10)
!    q21  = elimminf2_coli(p21)
!    q20  = elimminf2_coli(p20)
!
!    q1q2 = (q10+q20-q21)
!    detZ = 4d0*q10*q20-q1q2*q1q2

    if (abs(detZ/( 4d0*q10*q20 + Z(2,1)*Z(2,1))).lt.1d-4) then
      if (abs(q10-q20).lt.abs(q10-q21).and.  &
          abs(q10-q20).lt.abs(q20-q21)) then
        detZ  =  4d0*q10*q21 - (q10-q20+q21)*(q10-q20+q21)
      end if
    end if

!    Zadj(1,1) = 2d0*q20
!    Zadj(2,1) = -q1q2
!    Zadj(1,2) = -q1q2
!    Zadj(2,2) = 2d0*q10
!    f(1) = q10+mm02-mm12
!    f(2) = q20+mm02-mm22
!
!    Zadjf(1) = Zadj(1,1)*f(1)+Zadj(2,1)*f(2)
!    Zadjf(2) = Zadj(1,2)*f(1)+Zadj(2,2)*f(2)
!
! Xadj(1,1) and Xadj(2,2) exchanged!!!
!    Xadj(1,1) = 2d0*mm02*Z(1,1) - f(1)*f(1)
!    Xadj(2,1) = 2d0*mm02*Z(1,2) - f(1)*f(2)
!    Xadj(1,2) = Xadj(2,1)
!    Xadj(2,2) = 2d0*mm02*Z(2,2) - f(2)*f(2)


    ! coefficients Shat defined in (5.13)
    allocate(Shat(0:rmaxB,0:rmaxB,0:rmaxB,2))

    do r=0,rmaxB
      do n0=0,r/2

        do n1=0,r-2*n0
          n2 = r-2*n0-n1
          Shat(n0,n1,n2,:) = -B_0(n0,n1,n2)
        end do

        k = r-2*n0
        Shat(n0,0,k,1) = Shat(n0,0,k,1) + B_i(n0,k,1)
        Shat(n0,k,0,2) = Shat(n0,k,0,2) + B_i(n0,k,2)

      end do
    end do

    ! choose reduction formulas with biggest denominators
    if (abs(Xadj(1,1)).ge.abs(Xadj(2,2))) then
      if (abs(Xadj(1,1)).ge.abs(Xadj(1,2))) then
        a = 1
        b = 1
        at = 2
        bt = 2
        sgnab = 1
      else
        a = 1
        b = 2
        at = 2
        bt = 1
        sgnab = -1
      end if
    else        
      if (abs(Xadj(2,2)).ge.abs(Xadj(1,2))) then
        a = 2
        b = 2
        at = 1
        bt = 1
        sgnab = 1
      else 
        a = 1
        b = 2
        at = 2
        bt = 1
        sgnab = -1
      end if
    end if

    if (abs(Zadj(1,1)).ge.abs(Zadj(2,2))) then
      if (abs(Zadj(1,1)).ge.abs(Zadj(1,2))) then
        k = 1
        l = 1
        lt = 2
      else
        k = 1
        l = 2
        lt = 1
      end if
    else        
      if (abs(Zadj(2,2)).ge.abs(Zadj(1,2))) then
        k = 2
        l = 2
        lt = 1
      else 
        k = 1
        l = 2
        lt = 1
      end if
    end if

    ! allocation of array for det(Z)- and det(X)-expanded C-coefficients
    rmaxExp = rmaxB+1
    allocate(Cexpgy(0:max(rmax/2,1),0:rmaxExp-2,0:rmaxExp-2,0:ordgy_max))
   
    ! calculate Cuv
    allocate(CuvExpgy(0:rmaxExp,0:rmaxExp,0:rmaxExp))
    call CalcCuv(CuvExpgy,Buv_0,mm02,f,rmaxExp,id)
    Cuv(0:rmax,0:rmax,0:rmax) = CuvExpgy(0:rmax,0:rmax,0:rmax)

    ! allocate arrays for error propagation
    allocate(C00_err(0:rmaxExp))
    allocate(Cij_err(0:rmaxExp))
    allocate(C00_err2(0:rmaxExp))
    allocate(Cij_err2(0:rmaxExp))

    ! initialize accuracy estimates
    Cerr = acc_inf
    Cij_err =0d0
    C00_err =0d0

    Cerr2 = acc_inf
    Cij_err2 =0d0
    C00_err2 =0d0

!    maxZadjf = maxval(abs(Zadjf))
!    fmax = maxval(abs(f))

    ! truncation of expansion if calculated term larger than truncfacexp * previous term
    ! crucial for expansion parameters between 0.1 and 1 !!!
!    truncfacexp = sqrt(max(maxZadjf,abs(detZ))/abs(Xadj(a,b))*max(1d0,fmax/abs(Zadj(k,l)))) * truncfacC
    truncfacexp = sqrt(fac_gy) * truncfacC

    gtrunc = ordgy_max 

! calculate C(1,n1,n2) up to rank r+2
! calculate C(0,n1,n2) up to rank r  
    rloop: do r=0,rmaxExp-2

      if (r.gt.rmax+2*gtrunc+2) exit rloop

      !!!!!!!!!!!!!!!!!!!!!!!!!!!!
      ! 0th-order coefficients
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!

      ! calculating C_00ijk.. exploiting eq. (5.49)
      maxCexpgy(1,r,0)=0d0
      do nl=r,0,-1
        nlt=r-nl
        inds0(l) = nl
        inds0(lt) = nlt

        inds(l) = nl+1
        inds(lt) = nlt
        Caux = Zadj(k,1)*Shat(0,inds(1),inds(2),1)  &
             + Zadj(k,2)*Shat(0,inds(1),inds(2),2)

        if (nlt.ge.1) then
          inds(lt) = nlt-1
          Caux = Caux - 2*nlt*Zadj(k,lt)*Cexpgy(1,inds(1),inds(2),0)
        end if

        Cexpgy(1,inds0(1),inds0(2),0) = Caux/(2*(nl+1)*Zadj(k,l))
        maxCexpgy(1,r,0) =  maxCexpgy(1,r,0) + abs(Cexpgy(1,inds0(1),inds0(2),0) )
        if (r+2.le.rmax) then
          C(1,inds0(1),inds0(2)) = Cexpgy(1,inds0(1),inds0(2),0)
        end if

      end do

      ! calculate C_ijkl.. exploiting eq. (5.53)
      maxCexpgy(0,r,0)=0d0
      do n1=0,r
        n2 = r-n1
        inds(1) = n1
        inds(2) = n2
      
        Caux = (2*(2+r)*Cexpgy(1,n1,n2,0) - 4*CuvExpgy(1,n1,n2)  &
             - B_0(0,n1,n2))*Z(a,b)

!     write(*,*) 'CalcCred Caux',Caux,Z(a,b),f(a),f(b)

        Smod = Shat(0,n1,n2,a)

        if (inds(a).ge.1) then
          inds(a) = inds(a)-1
          Smod = Smod - 2d0*(inds(a)+1)*Cexpgy(1,inds(1),inds(2),0)

        end if

        Caux = Caux - f(b)*Smod

!       write(*,*) 'CalcCgy maxadjf',maxZadjf,Xadj(a,b),Caux

        Cexpgy(0,n1,n2,0) = Caux/Xadj(a,b)
        maxCexpgy(0,r,0) =  maxCexpgy(0,r,0) + abs(Cexpgy(0,n1,n2,0))
        if (r.le.rmax) then
          C(0,n1,n2) = Cexpgy(0,n1,n2,0)
          Cerr(r) =  abs(maxZadjf/Xadj(a,b)*Cexpgy(0,n1,n2,0))
        end if

      end do

      if (r.le.rmax) then
!       Cerr(r-1) =  abs(maxZadjf/Xadj(a,b))*maxCexpgy(0,r,0)
        Cerr(r-1) =  fac_gy*maxCexpgy(0,r,0)
      end if

 ! error propagation from B's
      C00_err(r+2) = B_err           
      Cij_err(r)=max(abs(Zadj(a,b))/abs(Xadj(a,b))*max(B_err,C00_err(r+2)),  &
                  fmax/abs(Xadj(a,b))*max(B_err,C00_err(r+1)))

      C00_err2(r+2) = B_err           
      Cij_err2(r)=max(abs(Zadj(a,b))/abs(Xadj(a,b))*max(B_err,C00_err(r+2)),  &
                  fmax/abs(Xadj(a,b))*max(B_err,C00_err2(r+1)))

      
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      ! higher order coefficients
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      rg = r
      gloop: do g=1,min(gtrunc,r/2)
        rg = rg-2

!        write(*,*) 'CalcCgy gtrunc gloop=',gtrunc,r,g,rg

        ! calculating C_00ijk.. exploiting eq. (5.49)
        maxCexpgy(1,rg,g) = 0d0
        do nl=rg,0,-1
          nlt=rg-nl
          inds0(l) = nl
          inds0(lt) = nlt

          inds(l) = nl+1
          inds(lt) = nlt
          Caux = -Zadjf(k)*Cexpgy(0,inds(1),inds(2),g-1)

          inds(k) = inds(k)+1
          Caux = Caux - detZ*Cexpgy(0,inds(1),inds(2),g-1)

          if (nlt.ge.1) then
            inds(l) = nl+1
            inds(lt) = nlt-1
            Caux = Caux - 2*nlt*Zadj(k,lt)*Cexpgy(1,inds(1),inds(2),g)
          end if

          Cexpgy(1,inds0(1),inds0(2),g) = Caux/(2*(nl+1)*Zadj(k,l))
          maxCexpgy(1,rg,g) =  maxCexpgy(1,rg,g) + abs(Cexpgy(1,inds0(1),inds0(2),g) )


          if (g.eq.1.and.abs(Cexpgy(1,inds0(1),inds0(2),g)).gt.      &
              truncfacexp*max(1d0,maxCexpgy(1,rg,g-1))    .or.   &
              g.ge.2.and.abs(Cexpgy(1,inds0(1),inds0(2),g)).gt.  &
              truncfacexp*maxCexpgy(1,rg,g-1)) then



                gtrunc = g-1
                exit gloop
!                gtrunc = g
!                cycle gloop

          end if 
        
        end do


        if (rg+2.le.rmax) then
          do nl=rg,0,-1
            nlt=rg-nl
            inds0(l) = nl
            inds0(lt) = nlt
            C(1,inds0(1),inds0(2)) = C(1,inds0(1),inds0(2))  &
                + Cexpgy(1,inds0(1),inds0(2),g)
          end do
        end if


        ! calculate C_ijkl.. exploiting eq. (5.53)
        maxCexpgy(0,rg,g) = 0d0
        do n1=0,rg
          n2 = rg-n1
          inds0(1) = n1
          inds0(2) = n2
          
          Caux = 2*(2+rg)*Cexpgy(1,n1,n2,g)*Z(a,b)
          
!         write(*,*) 'CalcCgy g Caux 1',rg,g,Caux

          if (inds0(a).ge.1) then
            inds = inds0
            inds(a) = inds(a)-1
            Caux = Caux + 2d0*f(b)*inds0(a)*Cexpgy(1,inds(1),inds(2),g)
          end if
          
!         write(*,*) 'CalcCgy g Caux 2',rg,g,Caux

          inds0(at) = inds0(at)+1
          Caux = Caux - sgnab*Zadjf(bt)*Cexpgy(0,inds0(1),inds0(2),g-1)
          
!         write(*,*) 'CalcCgy g Caux 3',rg,g,Caux

          Cexpgy(0,n1,n2,g) = Caux/Xadj(a,b)

!         write(*,*) 'CalcCgyo g Cexpgy',rg,g,n1,n2,Cexpgy(0,n1,n2,g)

          maxCexpgy(0,rg,g) =  maxCexpgy(0,rg,g) + abs(Cexpgy(0,n1,n2,g))
          
          if (g.eq.1.and.abs(Cexpgy(0,n1,n2,g)).gt.        &
              truncfacexp*max(1d0/m2scale,maxCexpgy(0,rg,g-1)).or.     &
              g.ge.2.and.abs(Cexpgy(0,n1,n2,g)).gt.        &
              truncfacexp*maxCexpgy(0,rg,g-1)) then



            gtrunc = g-1
            exit gloop
!           gtrunc = g
!           cycle gloop
          end if

!            if ((g.ge.2).and.(abs(Cexpgy(0,n1,n2,g)).gt.truncfacexp*abs(Cexpgy(0,n1,n2,g-1)))) then
!              gtrunc = g-1
!            end if

        end do

        ! error propagation from B's
        if(rg.gt.1)then
          C00_err(rg+2) =max(C00_err(rg+2),                                       &
              max(abs(Zadjf(k))*Cij_err(rg+1),abs(detZ)*Cij_err(rg+2))/abs(Zadj(k,l)))   
        end if
        Cij_err(rg)= max( Cij_err(rg),                                          &
            max(abs(Z(a,b))*C00_err(rg+2),abs(f(b))*C00_err(rg+1),       &
            abs(Zadjf(b))*Cij_err(rg+1))/abs(Xadj(a,b)))

        if(rg.gt.1)then
          C00_err2(rg+2) =max(C00_err2(rg+2),                                       &
              max(abs(Zadjf(k))*Cij_err2(rg+1),abs(detZ)*Cij_err2(rg+2))/abs(Zadj(k,l)))   
        end if
        Cij_err2(rg)= max( Cij_err2(rg),                                          &
            max(abs(Z(a,b))*C00_err2(rg+2),abs(f(b))*C00_err2(rg+1),       &
            abs(Zadjf(b))*Cij_err2(rg+1))/abs(Xadj(a,b)))



        if ((rg.le.rmax)) then
          Cerr(rg) = 0d0
          do n1=0,rg
            n2=rg-n1
            C(0,n1,n2) = C(0,n1,n2) + Cexpgy(0,n1,n2,g)
            if(abs(Cexpgy(0,n1,n2,g-1)).ne.0d0) then
              Cerr(rg)=max(Cerr(rg),abs(Cexpgy(0,n1,n2,g))*min(1d0,abs(Cexpgy(0,n1,n2,g))/abs(Cexpgy(0,n1,n2,g-1))))
            else
              Cerr(rg)=max(Cerr(rg),abs(Cexpgy(0,n1,n2,g)))
            end if
          end do

          ! if error from B's larger than error from expansion stop expansion
          if(Cij_err(rg).gt.Cerr(rg)) then
             gtrunc = min(g,gtrunc)
!            gtrunc = min(g+1,gtrunc)
             


          end if

        end if
        
        
      end do gloop

!     write(*,*) 'CalcCgy gtrunc after gloop=',gtrunc,r



      Cerr2 = max(Cerr,Cij_err2(0:rmax))
      Cerr = max(Cerr,Cij_err(0:rmax))



      ! check if target precision already reached
!      if(maxval(Cerr-acc_req_Cr*abs(C(0,0,0))).le.0d0) exit         ! changed 28.01.15

      if(maxval(Cerr-acc_req_Cr*abs(C(0,0,0))).le.0d0) then
        do rg=r+1,rmax
          do n1=0,rg
            C(0,n1,rg-n1)=0d0
          end do
        end do
        do rg=r+1,rmax
          do n1=0,rg-2
            C(1,n1,rg-2-n1)=0d0
          end do
        end do

        exit rloop
      end if

    end do rloop


    ! calculating C_0000ijk.. exploiting eq. (5.49)
    do r=4,rmax
      do n0=2,rmax/2
        do nl=r-2*n0,0,-1
          nlt=r-2*n0-nl
          inds0(l) = nl
          inds0(lt) = nlt

          inds(l) = nl+1
          inds(lt) = nlt
          Caux = Zadj(k,1)*Shat(n0-1,inds(1),inds(2),1)  &
               + Zadj(k,2)*Shat(n0-1,inds(1),inds(2),2)  &
               - Zadjf(k)*C(n0-1,inds(1),inds(2))

          inds(k) = inds(k)+1
          Caux = Caux - detZ*C(n0-1,inds(1),inds(2))

          if (nlt.ge.1) then
            inds(l) = nl+1
            inds(lt) = nlt-1
            Caux = Caux - 2*nlt*Zadj(k,lt)*C(n0,inds(1),inds(2))
          end if

          C(n0,inds0(1),inds0(2)) = Caux/(2*(nl+1)*Zadj(k,l))

        end do
      end do
    end do

      ! reduction formula (5.10) for n0+n1+n2=r, n0>0 
    do r=rmax+1,2*rmax
      do n0=r-rmax,r/2
        do n1=0,r-2*n0
          n2 = r-2*n0-n1
          C(n0,n1,n2) = (B_0(n0-1,n1,n2) + 2*mm02*C(n0-1,n1,n2) + 4*Cuv(n0,n1,n2) &
                        + f(1)*C(n0-1,n1+1,n2) + f(2)*C(n0-1,n1,n2+1)) / (2*r)
        end do
      end do
    end do



!   write(*,*) 'CalcCgyo out',(((C((r-n1-n2)/2,n1,n2),n2=0,r-n1),n1=0,r),r=0,rmax)



  end subroutine CalcCgyo





  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !  subroutine CalcCgp(C,Cuv,p10,p21,p20,m02,m12,m22,rmax,ordgp_min,ordgp_max,id,Cerr,acc_req_Cr,Cerr2)
  !
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  
  subroutine CalcCgp(C,Cuv,p10,p21,p20,m02,m12,m22,rmax,ordgp_min,ordgp_max,id,Cerr,acc_req_Cr,Cerr2)
  
    use globalC
  
    integer, intent(in) :: rmax,ordgp_min,ordgp_max,id
    double complex, intent(in) :: p10,p21,p20,m02,m12,m22
    double complex, intent(out) :: C(0:rmax,0:rmax,0:rmax)
    double complex, intent(out) :: Cuv(0:rmax,0:rmax,0:rmax)
    double precision, intent(out) :: Cerr(0:rmax),Cerr2(0:rmax)
    double precision, intent(in) :: acc_req_Cr(0:rmax)
    double complex, allocatable :: Cexpgp(:,:,:,:), CuvExpgp(:,:,:)
    double complex, allocatable :: B_0(:,:,:), B_k(:,:), Shat(:,:,:)
    double complex, allocatable :: Buv_0(:,:,:), Buv_k(:,:)
    double complex :: Smod, fk, elimminf2_coli
    double precision, allocatable :: C00_err(:),Cij_err(:)
    double precision, allocatable :: C00_err2(:),Cij_err2(:)
    double precision :: B_err,B_max
    double precision :: maxCexpgp(0:1,0:rmax+ordgp_min+1,0:ordgp_max),truncfacexp
    integer :: rmaxB,rmaxExp,gtrunc,r,n0,n1,n2,k,l,g,rg
    integer :: bin,nid(0:2),i
 


    ! determine Gram matrix
!    mm02 = elimminf2_coli(m02)
!    mm12 = elimminf2_coli(m12)
!    mm22 = elimminf2_coli(m22)
!    q10  = elimminf2_coli(p10)
!    q21  = elimminf2_coli(p21)
!    q20  = elimminf2_coli(p20)
!
!    q1q2 = (q10+q20-q21)
! commented out 2.9.17
!    Z(1,1) = 2d0*q10
!    Z(2,1) = q1q2
!    Z(1,2) = q1q2
!    Z(2,2) = 2d0*q20
!    f(1) = q10+mm02-mm12
!    f(2) = q20+mm02-mm22


    ! choose reduction formulas with biggest denominators
    if (abs(f(1)).ge.abs(f(2))) then
      k = 1
    else 
      k = 2
    end if
    fk = f(k)


    ! calculations of B-coefficients
    rmaxB = rmax + ordgp_min
    allocate(B_0(0:rmaxB,0:rmaxB,0:rmaxB))
    allocate(Buv_0(0:rmaxB,0:rmaxB,0:rmaxB))
    allocate(B_k(0:rmaxB,0:rmaxB))
    allocate(Buv_k(0:rmaxB,0:rmaxB))

    ! determine binaries for B-coefficients
    i=0
    bin = 1
    do while (i.le.2)
      if (mod(id/bin,2).eq.0) then
        nid(i) = id+bin
        i = i+1
      end if
      bin = 2*bin
    end do

    call CalcB(B_0(:,0,:),Buv_0(:,0,:),p21,m12,m22,rmaxB,nid(0))
    if (k.eq.1) then
      call CalcB(B_k(:,:),Buv_k(:,:),p20,m02,m22,rmaxB,nid(1))
    else
      call CalcB(B_k(:,:),Buv_k(:,:),p10,m02,m12,rmaxB,nid(2))
    end if

    ! shift of integration momentum in B_0
    B_max=0d0
    do n1=1,rmaxB
      do n2=0,rmaxB-n1
        n0 = (rmaxB-n1-n2)
        B_0(0:n0,n1,n2) = -B_0(0:n0,n1-1,n2)-B_0(0:n0,n1-1,n2+1)
        Buv_0(0:n0,n1,n2) = -Buv_0(0:n0,n1-1,n2)-Buv_0(0:n0,n1-1,n2+1)
        B_max = max(B_max,abs(B_0(0,n1,n2)))
     end do
    end do
    B_max = max(B_max,maxval(abs(B_k(0,0:rmaxB))))
    B_err = acc_def_B*B_max



    ! coefficients Shat defined in (5.13)
    allocate(Shat(0:rmaxB,0:rmaxB,0:rmaxB))

    do r=0,rmaxB
      do n0=0,r/2

        do n1=0,r-2*n0
          n2 = r-2*n0-n1
          Shat(n0,n1,n2) = -B_0(n0,n1,n2)
        end do

        l = r-2*n0
        if (k.eq.1) then
          Shat(n0,0,l) = Shat(n0,0,l) + B_k(n0,l)
        else
          Shat(n0,l,0) = Shat(n0,l,0) + B_k(n0,l)
        end if

      end do
    end do

        
    ! allocation of array for det(Z)-expanded C-coefficients
    rmaxExp = rmaxB+1
    allocate(Cexpgp(0:rmaxExp/2,0:rmaxExp-1,0:rmaxExp-1,0:ordgp_max))
   

    ! calculate Cuv
    allocate(CuvExpgp(0:rmaxExp,0:rmaxExp,0:rmaxExp))
    call CalcCuv(CuvExpgp,Buv_0,mm02,f,rmaxExp,id)
    Cuv(0:rmax,0:rmax,0:rmax) = CuvExpgp(0:rmax,0:rmax,0:rmax)

    ! allocate arrays for error propagation
    allocate(C00_err(0:rmaxExp))
    allocate(Cij_err(0:rmaxExp))
    allocate(C00_err2(0:rmaxExp))
    allocate(Cij_err2(0:rmaxExp))

    ! initialize accuracy estimates
    Cerr = acc_inf
    Cij_err =0d0
    C00_err =0d0

    Cerr2 = acc_inf
    Cij_err2 =0d0
    C00_err2 =0d0

!    maxZ = maxval(abs(Z))
!    maxZ = 2d0*q2max

    ! truncation of expansion if calculated term larger than truncfacexp * previous term
    ! crucial for expansion parameters between 0.1 and 1 !!!
!    truncfacexp = sqrt(abs(maxZ/abs(fk))) * truncfacC
    truncfacexp = sqrt(fac_gp) * truncfacC
    gtrunc = ordgp_max 



! calculate C(n0,n1,n2) up to rank r for n0>0 and up to rank r-1 for n0=0
    rloop: do r=1,rmaxExp




      if (r.gt.rmax+gtrunc+1) exit rloop

      !!!!!!!!!!!!!!!!!!!!!!!!!!!!
      ! 0th-order coefficients
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!

      ! calculating
      ! C_00(a)0000..00 --> C_00(a)ij00..00 --> C_00(a)ijkl00..00 --> ... --> C_00(a)ijklmn..
      ! exploiting eq. (5.63)
      maxCexpgp(1,r,0)=0d0
      do n0=r/2,1,-1
        do n1=0,r-2*n0
          n2=r-2*n0-n1

          Cexpgp(n0,n1,n2,0) = (2d0*CuvExpgp(n0,n1,n2) + B_0(n0-1,n1,n2)  &
               + mm02*Cexpgp(n0-1,n1,n2,0))/((r-n0)+1d0)/2d0

          if (n0.eq.1) then
            maxCexpgp(1,r,0) =  maxCexpgp(1,r,0) + abs(Cexpgp(n0,n1,n2,0) )
          end if

          if (r-n0.le.rmax) then
            C(n0,n1,n2) = Cexpgp(n0,n1,n2,0)
          end if

        end do
      end do

      ! calculate
      ! C_00ijkl.. --> C_aijkl..
      ! exploiting eq. (5.62)
      maxCexpgp(0,r-1,0)=0d0
      do n1=0,r-1
        n2 = r-1-n1

        Smod = Shat(0,n1,n2)
        if ((k.eq.1).and.(n1.ge.1)) then
          Smod = Smod - 2d0*n1*Cexpgp(1,n1-1,n2,0)
        else if ((k.eq.2).and.(n2.ge.1)) then
          Smod = Smod - 2d0*n2*Cexpgp(1,n1,n2-1,0)
        end if

        Cexpgp(0,n1,n2,0) = Smod/fk
        maxCexpgp(0,r-1,0) =  maxCexpgp(0,r-1,0) + abs(Cexpgp(0,n1,n2,0))
        
        if (r.le.rmax+1) then
          C(0,n1,n2) = Cexpgp(0,n1,n2,0)
        end if

      end do

      if (r.le.rmax+1) then
!       Cerr(r-1) =  abs(maxZ/fk)*maxCexpgp(0,r-1,0)
        Cerr(r-1) =  fac_gp*maxCexpgp(0,r-1,0)
      end if

      ! error propagation from B's
      if(r.gt.1)then
        C00_err(r) = B_err/(2*r)
      end if
      Cij_err(r-1) = B_err/abs(fk)

      if(r.gt.1)then
        C00_err2(r) = B_err/(2*r)
      end if
      Cij_err2(r-1) = B_err/abs(fk)



      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      ! higher order coefficients
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      rg = r
      gloop: do g=1,min(gtrunc,r-1)
        rg = rg-1

        ! calculating for rank=rmaxB+1
        ! C_00(a)0000..00 --> C_00(a)ij00..00 --> C_00(a)ijkl00..00 --> ... --> C_00(a)ijklmn..
        ! exploiting eq. (5.63)
        maxCexpgp(1,rg,g) = 0d0
        do n0=rg/2,1,-1
          do n1=0,rg-2*n0
            n2=rg-2*n0-n1

            Cexpgp(n0,n1,n2,g) = (2d0*mm02*Cexpgp(n0-1,n1,n2,g)  &
               - Z(1,1)*Cexpgp(n0-1,n1+2,n2,g-1) - 2d0*Z(2,1)*Cexpgp(n0-1,n1+1,n2+1,g-1)  &
               - Z(2,2)*Cexpgp(n0-1,n1,n2+2,g-1))/((rg-n0)+1d0)/4d0

            if(n0.eq.1) then
              maxCexpgp(1,rg,g) =  maxCexpgp(1,rg,g) + abs(Cexpgp(n0,n1,n2,g))
              

              if (g.eq.1.and.abs(Cexpgp(1,n1,n2,g)).gt.          &
                  truncfacexp*max(1d0,maxCexpgp(1,rg,g-1)) .or.  &
                  g.ge.2.and.abs(Cexpgp(1,n1,n2,g)).gt.          &
                  truncfacexp*maxCexpgp(1,rg,g-1)) then



                gtrunc = g-1
                exit gloop
              end if
            end if

!            if ((g.ge.2).and.(abs(Cexpgp(n0,n1,n2,g)).gt.truncfacexp*abs(Cexpgp(n0,n1,n2,g-1)))) then
!              gtrunc = g-1
!            end if

          end do
        end do


        do n0=rg/2,1,-1
          if (rg-n0.le.rmax) then
            do n1=0,rg-2*n0
              n2=rg-2*n0-n1
              C(n0,n1,n2) = C(n0,n1,n2) + Cexpgp(n0,n1,n2,g)
            end do
          end if
        end do


        ! calculate
        ! C_000000..00 --> C_i0000..00 --> C_ij00..00 --> ... --> C_ijk..
        ! exploiting eq. (5.62)
        maxCexpgp(0,rg-1,g) = 0d0
        do n1=0,rg-1
          n2 = rg-1-n1
          
          Smod = -Z(1,k)*Cexpgp(0,n1+1,n2,g-1)  &
              -Z(2,k)*Cexpgp(0,n1,n2+1,g-1)
          if ((k.eq.1).and.(n1.ge.1)) then
            Smod = Smod - 2d0*n1*Cexpgp(1,n1-1,n2,g)
          else if ((k.eq.2).and.(n2.ge.1)) then
            Smod = Smod - 2d0*n2*Cexpgp(1,n1,n2-1,g)
          end if
          
          Cexpgp(0,n1,n2,g) = Smod/fk
          
          maxCexpgp(0,rg-1,g) =  maxCexpgp(0,rg-1,g) + abs(Cexpgp(0,n1,n2,g))
          
          if (g.eq.1.and.abs(Cexpgp(0,n1,n2,g)).gt.                     &
!     corrected 02.07.2018
              truncfacexp*max(1/m2scale,maxCexpgp(0,rg-1,g-1))   .or.    &
!              truncfacexp*max(1/m2max,maxCexpgp(0,rg-1,g-1))    .or.    &
              g.ge.2.and.abs(Cexpgp(0,n1,n2,g)).gt.                     &
              truncfacexp*maxCexpgp(0,rg-1,g-1)) then


            gtrunc = g-1
            exit gloop
          end if
          
        end do

        ! error propagation from B's
        if(rg.gt.1)then
          C00_err(rg) = max(C00_err(rg),max(2*abs(m02)*Cij_err(rg-2),maxZ*Cij_err(rg))/(4*r) )
        end if
        Cij_err(rg-1) = max(Cij_err(rg-1),max(2*C00_err(rg),maxZ*Cij_err(rg))/abs(fk) )

        if(rg.gt.1)then
          C00_err2(rg) = max(C00_err2(rg),max(2*abs(m02)*Cij_err2(rg-2),maxZ*Cij_err2(rg))/(4*r) )
        end if
        Cij_err2(rg-1) = max(Cij_err2(rg-1),max(2*C00_err2(rg),maxZ*Cij_err2(rg))/abs(fk) )



        if ((rg.le.rmax+1)) then
          Cerr(rg-1) = 0d0
          do n1=0,rg-1
            n2=rg-1-n1
            C(0,n1,n2) = C(0,n1,n2) + Cexpgp(0,n1,n2,g)
!            Cerr(rg-1)=max(Cerr(rg-1),abs(Cexpgp(0,n1,n2,g))**2/abs(Cexpgp(0,n1,n2,g-1)))
            if(abs(Cexpgp(0,n1,n2,g-1)).ne.0d0) then
              Cerr(rg-1)=max(Cerr(rg-1),abs(Cexpgp(0,n1,n2,g))*min(1d0,abs(Cexpgp(0,n1,n2,g))/abs(Cexpgp(0,n1,n2,g-1))))
            else
              Cerr(rg-1)=max(Cerr(rg-1),abs(Cexpgp(0,n1,n2,g)))
            end if
          end do

          ! if error from B's larger than error from expansion stop expansion
          if(Cij_err(rg-1).gt.Cerr(rg-1)) then
             gtrunc = min(g,gtrunc)
             


          end if

        end if

      end do gloop



      Cerr2 = max(Cerr,Cij_err2(0:rmax))
      Cerr = max(Cerr,Cij_err(0:rmax))
      


      ! check if target precision already reached
!      if(maxval(Cerr-acc_req_Cr*abs(C(0,0,0))).le.0d0) exit    ! changed 28.01.15

      if(maxval(Cerr-acc_req_Cr*abs(C(0,0,0))).le.0d0) then
        do rg=r+1,rmax
          do n0=0,rg/2
            do n1=0,rg-2*n0
              C(n0,n1,rg-2*n0-n1)=0d0
            end do
          end do
        end do
        if(r.le.rmax) then
          do n1=0,r
            C(0,n1,r-n1)=0d0
          end do
        end if   

!       write(*,*) 'CalcCg exit rloop  =',r,rmax,rg

        exit rloop
      end if

!      write(*,*) 'CalcCgp after exit'

    end do rloop

      ! reduction formula (5.10) for n0+n1+n2=r, n0>0 
    do r=rmax+1,2*rmax
      do n0=r-rmax,r/2
        do n1=0,r-2*n0
          n2 = r-2*n0-n1
          C(n0,n1,n2) = (B_0(n0-1,n1,n2) + 2*mm02*C(n0-1,n1,n2) + 4*Cuv(n0,n1,n2) &
                        + f(1)*C(n0-1,n1+1,n2) + f(2)*C(n0-1,n1,n2+1)) / (2*r)
        end do
      end do
    end do


              


!   write(*,*) 'CalcCgp rmax',rmax
!   do r=14,rmax
!    do r=0,rmax
!   do n0=0,r/2
!   do n1=0,r-2*n0
!   write(*,*) 'CalcCgp out',r,n0,n1,r-2*n0-n1,C(n0,n1,r-2*n0-n1)
!   end do
!   end do
!   end do

!    write(*,*) 'CalcCgp Cerr ',Cerr
!    write(*,*) 'CalcCgp Cerr2',Cerr2

  end subroutine CalcCgp


  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !  subroutine CalcCgpf(C,Cuv,p10,p21,p20,m02,m12,m22,rmax,ordgpf_min,ordgpf_max,id,Cerr,acc_req_Cr,Cerr2)
  !
  !  added by AD  16.08.2017
  !
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine CalcCgpf(C,Cuv,p10,p21,p20,m02,m12,m22,rmax,ordgpf_min,ordgpf_max,id,Cerr,acc_req_Cr,Cerr2)
  
    use globalC
  
    integer, intent(in) :: rmax,ordgpf_min,ordgpf_max,id
    double complex, intent(in) :: p10,p21,p20,m02,m12,m22
    double complex, intent(out) :: C(0:rmax,0:rmax,0:rmax)
    double complex, intent(out) :: Cuv(0:rmax,0:rmax,0:rmax)
    double precision, intent(out) :: Cerr(0:rmax),Cerr2(0:rmax)
    double precision, intent(in) :: acc_req_Cr(0:rmax)
    double complex, allocatable :: Cexpgpf(:,:,:,:), CuvExpgpf(:,:,:)
    double complex, allocatable :: B_0(:,:,:), B_i(:,:,:), Shat(:,:,:,:)
    double complex, allocatable :: Buv_0(:,:,:), Buv_i(:,:,:)
    double complex :: Smod, Caux, Zadj2f
    double complex :: C0_coli, elimminf2_coli
    double precision, allocatable :: C00_err(:),Cij_err(:)
    double precision, allocatable :: C00_err2(:),Cij_err2(:)
    double precision :: B_err,B_max,aZadj2f
    double precision :: maxCexpgpf(0:1,0:rmax+2*ordgpf_min,0:ordgpf_max),truncfacexp
    double precision :: minZk 
    integer :: rmaxB,rmaxExp,gtrunc,r,n0,n1,n2,i,j,jt,g,rg
    integer :: inds0(2),inds(2),k,l,lt,nl,nlt
    integer :: bin,nid(0:2)



    ! write(*,*) 'LH: CalcCgpf, ord', ordgpf_min
    ! calculation of B-coefficients
    rmaxB = rmax + 2*ordgpf_min + 1
    allocate(B_0(0:rmaxB,0:rmaxB,0:rmaxB))
    allocate(Buv_0(0:rmaxB,0:rmaxB,0:rmaxB))
    allocate(B_i(0:rmaxB,0:rmaxB,2))
    allocate(Buv_i(0:rmaxB,0:rmaxB,2))

    ! determine binaries for B-coefficients
    k=0
    bin = 1
    do while (k.le.2)
      if (mod(id/bin,2).eq.0) then
        nid(k) = id+bin
        k = k+1
      end if
      bin = 2*bin
    end do

    call CalcB(B_0(:,0,:),Buv_0(:,0,:),p21,m12,m22,rmaxB,nid(0))
    call CalcB(B_i(:,:,1),Buv_i(:,:,1),p20,m02,m22,rmaxB,nid(1))
    call CalcB(B_i(:,:,2),Buv_i(:,:,2),p10,m02,m12,rmaxB,nid(2))

   ! shift of integration momentum in B_0
    B_max=0d0
    do n1=1,rmaxB
      do n2=0,rmaxB-n1
        n0 = (rmaxB-n1-n2)
        B_0(0:n0,n1,n2) = -B_0(0:n0,n1-1,n2)-B_0(0:n0,n1-1,n2+1)
        Buv_0(0:n0,n1,n2) = -Buv_0(0:n0,n1-1,n2)-Buv_0(0:n0,n1-1,n2+1)
      end do
    end do
    B_max = max(B_max,maxval(abs(B_i(0,0:rmaxB,1:2))))
    B_err = acc_def_B*B_max

    ! determine (adjugated) Gram and Cayley matrix
!    mm02 = elimminf2_coli(m02)
!    mm12 = elimminf2_coli(m12)
!    mm22 = elimminf2_coli(m22)
!    q10  = elimminf2_coli(p10)
!    q21  = elimminf2_coli(p21)
!    q20  = elimminf2_coli(p20)
!
!    q1q2 = (q10+q20-q21)
!    detZ = 4d0*q10*q20-q1q2*q1q2

    if (abs(detZ).lt.abs(4d0*q10*q20 + Z(2,1)*Z(2,1))*1d-4) then
      if (abs(q10-q20).lt.abs(q10-q21).and.  &
          abs(q10-q20).lt.abs(q20-q21)) then
        detZ  =  4d0*q10*q21 - (q10-q20+q21)*(q10-q20+q21)
      end if
    end if

!    Zadj(1,1) = 2d0*q20
!    Zadj(2,1) = -q1q2
!    Zadj(1,2) = -q1q2
!    Zadj(2,2) = 2d0*q10
!    f(1) = q10+mm02-mm12
!    f(2) = q20+mm02-mm22
!
!    Zadjf(1) = Zadj(1,1)*f(1)+Zadj(2,1)*f(2)
!    Zadjf(2) = Zadj(1,2)*f(1)+Zadj(2,2)*f(2)
!
! Xadj(1,1) and Xadj(2,2) exchanged!!!
!    Xadj(1,1) = 2d0*mm02*Z(1,1) - f(1)*f(1)
!    Xadj(2,1) = 2d0*mm02*Z(1,2) - f(1)*f(2)
!    Xadj(1,2) = Xadj(2,1)
!    Xadj(2,2) = 2d0*mm02*Z(2,2) - f(2)*f(2)


    ! coefficients Shat defined in (5.13)
    allocate(Shat(0:rmaxB,0:rmaxB,0:rmaxB,2))

    do r=0,rmaxB
      do n0=0,r/2

        do n1=0,r-2*n0
          n2 = r-2*n0-n1
          Shat(n0,n1,n2,:) = -B_0(n0,n1,n2)
        end do

        k = r-2*n0
        Shat(n0,0,k,1) = Shat(n0,0,k,1) + B_i(n0,k,1)
        Shat(n0,k,0,2) = Shat(n0,k,0,2) + B_i(n0,k,2)

      end do
    end do

    ! choose reduction formulas with smallest expansion terms
    minZk = maxZ
    if (maxval(abs(Z(1,1:2))).le.minZk) then
      minZk = maxval(abs(Z(1,1:2)))
      k = 1
      l = 1
      lt = 2
    end if
    if (maxval(abs(Z(2,1:2))).lt.minZk) then
      minZk = maxval(abs(Z(2,1:2)))
      k = 2
      l = 2
      lt = 1 
    end if



!   write(*,*)  'CalcCgpf Zadj(i,j)=',i,j,Zadj(i,j),Xadj(i,j)
        
    ! allocation of array for det(Z)- and det(X)-expanded C-coefficients
    rmaxExp = rmaxB+1
    allocate(Cexpgpf(0:max(rmax/2,1),0:rmaxExp-2,0:rmaxExp-2,0:ordgpf_max))
   
    ! calculate Cuv
    allocate(CuvExpgpf(0:rmaxExp,0:rmaxExp,0:rmaxExp))
    call CalcCuv(CuvExpgpf,Buv_0,mm02,f,rmaxExp,id)
    Cuv(0:rmax,0:rmax,0:rmax) = CuvExpgpf(0:rmax,0:rmax,0:rmax)

    ! allocate arrays for error propagation
    allocate(C00_err(0:rmaxExp))
    allocate(Cij_err(0:rmaxExp))
    allocate(C00_err2(0:rmaxExp))
    allocate(Cij_err2(0:rmaxExp))

    ! initialize accuracy estimates
    Cerr = acc_inf
    Cij_err =0d0
    C00_err =0d0

    Cerr2 = acc_inf
    Cij_err2 =0d0
    C00_err2 =0d0

!    maxZadjf = maxval(abs(Zadjf))
!    fmax = maxval(abs(f))

    ! truncation of expansion if calculated term larger than truncfacexp * previous term
    ! crucial for expansion parameters between 0.1 and 1 !!!
!    truncfacexp = sqrt(max(maxZadjf,abs(detZ))/abs(Xadj(i,j))*max(1d0,fmax/abs(Zadj(k,l)))) * truncfacC
    truncfacexp = sqrt(fac_gpf) * truncfacC
    gtrunc = ordgpf_max 



! calculate C(1,n1,n2) up to rank r+2
! calculate C(0,n1,n2) up to rank r  
    rloop: do r=0,rmaxExp-2



      if (r.gt.rmax+2*gtrunc+2) exit rloop

      !!!!!!!!!!!!!!!!!!!!!!!!!!!!
      ! 0th-order coefficients
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!

      ! calculating C_00ijk.. exploiting eq. (5.71)
      maxCexpgpf(1,r,0)=0d0
      do nl=r,0,-1
        nlt=r-nl
        inds0(l) = nl
        inds0(lt) = nlt

        inds(l) = nl+1
        inds(lt) = nlt

        Caux = Shat(0,inds(1),inds(2),k)

        Cexpgpf(1,inds0(1),inds0(2),0) = Caux/(2*(nl+1))
        maxCexpgpf(1,r,0) =  maxCexpgpf(1,r,0) + abs(Cexpgpf(1,inds0(1),inds0(2),0) )
!        if (r+2.le.rmax) then             !  for fixed rank
        if (r+1.le.rmax) then
          C(1,inds0(1),inds0(2)) = Cexpgpf(1,inds0(1),inds0(2),0)
        end if

      end do

      ! calculate C_ijkl.. exploiting eq. (5.72)
      maxCexpgpf(0,r,0)=0d0
      do n1=0,r
        n2 = r-n1
        inds(1) = n1
        inds(2) = n2
      
        Caux = 2*(4+r+r)*Cexpgpf(1,n1,n2,0) - 4*CuvExpgpf(1,n1,n2)  &
             - 2*B_0(0,n1,n2)

        Cexpgpf(0,n1,n2,0) = Caux/(2d0*m02)

        maxCexpgpf(0,r,0) =  maxCexpgpf(0,r,0) + abs(Cexpgpf(0,n1,n2,0))
        if (r.le.rmax) then
          C(0,n1,n2) = Cexpgpf(0,n1,n2,0)
        end if

      end do

      if (r.le.rmax) then
!       Cerr(r) =  abs(maxZadjf/Xadj(i,j))*maxCexpgpf(0,r,0)
        Cerr(r) =  fac_gpf*maxCexpgpf(0,r,0)

!        write(*,*) 'CalcCgpf Cerr,0 ',r,Cerr(r),fac_gpf,maxCexpgpf(0,r,0)

      end if

 ! error propagation from B's
      C00_err(r+2) = B_err /2d0               
      Cij_err(r) = max(B_err,2*(r+2)*C00_err(r+2))/abs(m02)

      C00_err2(r+2) = B_err /2d0               
      Cij_err2(r) = max(B_err,2*(r+2)*C00_err2(r+2))/abs(m02)
      






      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      ! higher order coefficients
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      rg = r
      gloop: do g=1,min(gtrunc,r/2)
        rg = rg-2

        ! calculating C_00ijk.. exploiting eq. (5.71)
        maxCexpgpf(1,rg,g) = 0d0
        do nl=rg,0,-1
          nlt=rg-nl
          inds0(l) = nl
          inds0(lt) = nlt

          inds(l) = nl+1
          inds(lt) = nlt
          Caux = -f(k)*Cexpgpf(0,inds(1),inds(2),g-1)

          inds(l) = inds(l)+1
          Caux = Caux - Z(k,l)*Cexpgpf(0,inds(1),inds(2),g-1)

          inds(l) =  inds(l)-1
          inds(lt) =  inds(lt)+1
          Caux = Caux - Z(k,lt)*Cexpgpf(0,inds(1),inds(2),g-1)

          Cexpgpf(1,inds0(1),inds0(2),g) = Caux/(2*(nl+1))

          maxCexpgpf(1,rg,g) =  maxCexpgpf(1,rg,g) + abs(Cexpgpf(1,inds0(1),inds0(2),g) )

          if (g.eq.1.and.abs(Cexpgpf(1,inds0(1),inds0(2),g)).gt.      &
                  truncfacexp*max(1d0,maxCexpgpf(1,rg,g-1))    .or.   &
                  g.ge.2.and.abs(Cexpgpf(1,inds0(1),inds0(2),g)).gt.  &
                  truncfacexp*maxCexpgpf(1,rg,g-1)) then


            gtrunc = g-1
            exit gloop

          end if 
        
        end do


!        if (rg+2.le.rmax) then            !  for fixed rank
        if (rg+1.le.rmax) then
          do nl=rg,0,-1
            nlt=rg-nl
            inds0(l) = nl
            inds0(lt) = nlt
            C(1,inds0(1),inds0(2)) = C(1,inds0(1),inds0(2))  &
                + Cexpgpf(1,inds0(1),inds0(2),g)
          end do
        end if


        ! calculate C_ijkl.. exploiting eq. (5.72)
        maxCexpgpf(0,rg,g) = 0d0
        do n1=0,rg
          n2 = rg-n1
          inds(1) = n1
          inds(2) = n2
          
          Caux = 2*(4+rg+rg)*Cexpgpf(1,n1,n2,g)
          
          do i=1,2
          do j=1,2
            inds(i)=inds(i)+1
            inds(j)=inds(j)+1 
            Caux = Caux + Z(i,j)*Cexpgpf(0,inds(1),inds(2),g-1)
            inds(i)=inds(i)-1
            inds(j)=inds(j)-1  
          end do
          end do

          Cexpgpf(0,n1,n2,g) = Caux/(2*m02)

          maxCexpgpf(0,rg,g) =  maxCexpgpf(0,rg,g) + abs(Cexpgpf(0,n1,n2,g))
          
          if (g.eq.1.and.abs(Cexpgpf(0,n1,n2,g)).gt.        &
              truncfacexp*max(1d0/m2scale,maxCexpgpf(0,rg,g-1)).or.     &
              g.ge.2.and.abs(Cexpgpf(0,n1,n2,g)).gt.        &
              truncfacexp*maxCexpgpf(0,rg,g-1)) then



            gtrunc = g-1
            exit gloop

          end if

!            if ((g.ge.2).and.(abs(Cexpgpf(0,n1,n2,g)).gt.truncfacexp*abs(Cexpgpf(0,n1,n2,g-1)))) then
!              gtrunc = g-1
!            end if

        end do



        ! error propagation from B's
        if(rg.gt.1)then
          C00_err(rg+2) =max(C00_err(rg+2),                           &
              fmax/2d0*Cij_err(rg+1),                                 &
              maxZ/2d0*Cij_err(rg+2))   
        end if



        Cij_err(rg)= max( Cij_err(rg),                                &
            2*(rg+2)/abs(m02)*C00_err(rg+2),                          &
            maxZ/(2*abs(m02))*Cij_err(rg+2))

        if(rg.gt.1)then
          C00_err2(rg+2) =max(C00_err2(rg+2),                         &
              fmax/2d0*Cij_err(rg+1),                                 &
              maxZ/2d0*Cij_err(rg+2))   
        end if

        Cij_err2(rg)= max( Cij_err2(rg),                              &
            2*(rg+2)/abs(m02)*C00_err2(rg+2),                         &
            maxZ/(2*abs(m02))*Cij_err2(rg+2))



        if ((rg.le.rmax)) then
          Cerr(rg) = 0d0
          do n1=0,rg
            n2=rg-n1
            C(0,n1,n2) = C(0,n1,n2) + Cexpgpf(0,n1,n2,g)



            if(abs(Cexpgpf(0,n1,n2,g-1)).ne.0d0) then
              Cerr(rg)=max(Cerr(rg),abs(Cexpgpf(0,n1,n2,g))*min(1d0,abs(Cexpgpf(0,n1,n2,g))/abs(Cexpgpf(0,n1,n2,g-1))))
            else
              Cerr(rg)=max(Cerr(rg),abs(Cexpgpf(0,n1,n2,g)))
            end if





          end do

          ! if error from B's larger than error from expansion stop expansion
          if(Cij_err(rg).gt.Cerr(rg)) then
             gtrunc = min(g,gtrunc)
!            gtrunc = min(g+1,gtrunc)
             


          end if

        end if
        
      end do gloop

!     write(*,*) 'CalcCgpf gtrunc aft gloop=',gtrunc,r



      Cerr2 = max(Cerr,Cij_err2(0:rmax))
      Cerr = max(Cerr,Cij_err(0:rmax))



      ! check if target precision already reached
!      if(maxval(Cerr-acc_req_Cr*abs(C(0,0,0))).le.0d0) exit         ! changed 28.01.15

      if(maxval(Cerr-acc_req_Cr*abs(C(0,0,0))).le.0d0) then
        do rg=r+1,rmax
          do n1=0,rg
            C(0,n1,rg-n1)=0d0
          end do
        end do   
        do rg=r+1,rmax
          do n1=0,rg-2
            C(1,n1,rg-2-n1)=0d0
          end do
        end do   




        exit rloop 

      end if

    end do rloop


    ! calculating C_0000ijk.. exploiting eq. (5.71)
    do r=4,rmax
!      do n0=2,rmax/2     !     for fixed rank
      do n0=2,rmax
        do nl=r-2*n0,0,-1
          nlt=r-2*n0-nl
          inds0(l) = nl
          inds0(lt) = nlt

          inds(l) = nl+1
          inds(lt) = nlt
          Caux = Shat(n0-1,inds(1),inds(2),k)            &
                 - f(k)*C(n0-1,inds(1),inds(2))          &
                 - Z(k,1)*C(n0-1,inds(1)+1,inds(2))      &
                 - Z(k,2)*C(n0-1,inds(1),inds(2)+1)        

          C(n0,inds0(1),inds0(2)) = Caux/(2*(nl+1))

        end do
      end do
    end do

      ! reduction formula (5.10) for n0+n1+n2=r, n0>0 
    do r=rmax+1,2*rmax
      do n0=r-rmax,r/2
        do n1=0,r-2*n0
          n2 = r-2*n0-n1
          C(n0,n1,n2) = (B_0(n0-1,n1,n2) + 2*mm02*C(n0-1,n1,n2) + 4*Cuv(n0,n1,n2) &
                        + f(1)*C(n0-1,n1+1,n2) + f(2)*C(n0-1,n1,n2+1)) / (2*r)
        end do
      end do
    end do



!   write(*,*) 'CalcCgpf out',(((C((r-n1-n2)/2,n1,n2),n2=0,r-n1),n1=0,r),r=0,rmax)


!    write(*,*) 'CalcCgpf Cerr ',Cerr
!    write(*,*) 'CalcCgpf Cerr2',Cerr2

  end subroutine CalcCgpf








  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !  subroutine CopyCimp3(C,C_alt,Cerr,Cerr_alt,Cerr1,Cerr1_alt,Cerr2,Cerr2_alt,Crmethod,Crmethod_alt,rmax)
  !
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  
  subroutine CopyCimp3(C,C_alt,Cerr,Cerr_alt,Cerr1,Cerr1_alt,Cerr2,Cerr2_alt,Crmethod,Crmethod_alt,rmax,r_alt)

    integer,   intent(in) :: rmax,r_alt
    double complex, intent(inout) :: C(0:rmax,0:rmax,0:rmax)
    double precision, intent(inout) :: Cerr(0:rmax),Cerr1(0:rmax),Cerr2(0:rmax)
    integer, intent(inout) :: Crmethod(0:rmax)
    double complex, intent(in) :: C_alt(0:r_alt,0:r_alt,0:r_alt)
    double precision, intent(in) :: Cerr_alt(0:r_alt),Cerr2_alt(0:r_alt),Cerr1_alt(0:r_alt)
    integer, intent(in) :: Crmethod_alt(0:rmax)

    integer :: r,n1,n0

!    write(*,*) 'CopyCimp3: Cerr =',Cerr
!    write(*,*) 'CopyCimp3: Cerr_alt =',Cerr_alt

    do r=0,r_alt
      if (Cerr_alt(r).lt.Cerr(r)) then
        Crmethod(r)=Crmethod_alt(r)
        Cerr(r)=Cerr_alt(r)
        Cerr1(r)=Cerr1_alt(r)
        Cerr2(r)=Cerr2_alt(r)
        forall (n0=0:r)
          forall (n1=0:r-n0)
            C(n0,n1,r-n0-n1) = C_alt(n0,n1,r-n0-n1)
          end forall
        end forall
!        forall (n1=0:r)
!          forall (n2=0:r-n1)
!            C((r-n1-n2)/2,n1,n2) = C_alt((r-n1-n2)/2,n1,n2)
!          end forall
!        end forall
      end if
    end do

   end subroutine CopyCimp3


end module reductionC


