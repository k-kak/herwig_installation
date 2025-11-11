!******************************************************************************!
! Copyright (C) 2014-2019 OpenLoops Collaboration. For authors see authors.txt !
!                                                                              !
! This file is part of OpenLoops.                                              !
!                                                                              !
! OpenLoops is free software: you can redistribute it and/or modify            !
! it under the terms of the GNU General Public License as published by         !
! the Free Software Foundation, either version 3 of the License, or            !
! (at your option) any later version.                                          !
!                                                                              !
! OpenLoops is distributed in the hope that it will be useful,                 !
! but WITHOUT ANY WARRANTY; without even the implied warranty of               !
! MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the                !
! GNU General Public License for more details.                                 !
!                                                                              !
! You should have received a copy of the GNU General Public License            !
! along with OpenLoops.  If not, see <http://www.gnu.org/licenses/>.           !
!******************************************************************************!



module ol_self_energy_integrals_dp
! !interfaces for scalar one-point & two-point functions
  use kind_types, only: dp
  use ol_parameters_decl_dp, only: ZERO
  use ol_parameters_decl_dp, only: cms_on
  use ol_loop_parameters_decl_dp, only: de1_UV, de1_IR
  use ol_loop_parameters_decl_dp, only: se_integral_switch, coli_cache_use, &
    & a_switch

    use avh_olo_dp

  implicit none

  real(dp), private :: eps=1.e-17


  contains

  subroutine init_ol_self_energy_integrals(init)

  use collier, only: setmode_cll
  use cache, only: SwitchOffCacheSystem_cll, SwitchOnCacheSystem_cll

  implicit none
  logical :: init
  if(init) then
    if ((se_integral_switch == 1 .or. se_integral_switch == 7) &
      .and. coli_cache_use == 1) then
        call SwitchOffCacheSystem_cll
    end if
    if (se_integral_switch == 1 .and. a_switch /= 1) then
      call setmode_cll(1)
    end if
    if (se_integral_switch == 7 .and. a_switch /= 7) then
      call setmode_cll(2)
    end if
  else
    if (se_integral_switch == 1 .and. a_switch == 1) then
      call setmode_cll(1)
    end if
    if (se_integral_switch == 7 .and. a_switch /= 1) then
      call setmode_cll(2)
    end if
    if ((se_integral_switch == 1 .or. se_integral_switch == 7) &
      .and. coli_cache_use == 1) then
        call SwitchOnCacheSystem_cll
    end if
  end if
  end subroutine


  function calcA0(m12_in)

    use collier_coefs, only: A0_cll

    use avh_olo_dp

    complex(dp) calcA0
    complex(dp), intent(in) :: m12_in
    complex(dp) m12
    complex(dp) A0_coli

    complex(dp) :: rslt(0:2)


    calcA0 = 0

    m12 = m12_in


    if (se_integral_switch == 1 .or. se_integral_switch == 7) then
      call A0_cll(A0_coli,m12)
      calcA0 = A0_coli
    end if



    if (se_integral_switch == 3) then
      call olo_a0(rslt,m12_in)
      calcA0 = rslt(0) + rslt(1)*de1_UV
    end if


      return
  end function calcA0

  function calcB0(p2_in,m12_in,m22_in)

    use collier_coefs, only: B0_cll

    use avh_olo_dp

    complex(dp) calcB0
    complex(dp), intent(in) :: p2_in
    complex(dp), intent(in) :: m12_in
    complex(dp), intent(in) :: m22_in
    complex(dp) p2q
    complex(dp) p2
    complex(dp) m12
    complex(dp) m22
    complex(dp) B0_coli

    complex(dp) :: rslt(0:2)


    calcB0 = 0

    p2  = p2_in
    m12 = m12_in
    m22 = m22_in


    if (se_integral_switch == 1 .or. se_integral_switch == 7) then
      p2 = real(p2)
      call B0_cll(B0_coli,p2,m12,m22)
      calcB0 = B0_coli
    end if



    if (se_integral_switch == 3) then
      p2q = real(p2_in)
      call olo_b0(rslt,real(p2q),m12_in,m22_in)
      if (p2q .eq. 0 .and. m12_in .eq. 0 .and. m22_in .eq. 0) then
        calcB0 = de1_UV - de1_IR
      else
        calcB0 = rslt(0) + rslt(1)*de1_UV
      end if
    end if


      return
  end function calcB0


  function calcdB0(p2_in,m12_in,m22_in)

    use collier_coefs, only: DB0_cll

    implicit none
    complex(dp) calcdB0
    complex(dp), intent(in) :: p2_in
    complex(dp), intent(in) :: m12_in
    complex(dp), intent(in) :: m22_in
    real(dp) :: p2q
    complex(dp) :: p2
    complex(dp) :: m12
    complex(dp) :: m22
    complex(dp) DB0_coli

    complex(dp) :: rslt(0:2)


    calcdB0 = 0

    p2  = p2_in
    m12 = m12_in
    m22 = m22_in


    if (se_integral_switch == 1 .or. se_integral_switch == 7) then
      p2 = real(p2)
      call DB0_cll(DB0_coli,p2,m12,m22)
      calcdB0 = DB0_coli
    end if



    if (se_integral_switch == 3) then
      p2q = real(p2_in)
      if (p2q == 0. .and. m12 == 0. .and. m22 == 0. ) then
        calcdB0 = 0.
      else
        call olo_db0(rslt,p2q,m12_in,m22_in)
        calcdB0 = rslt(0) + rslt(1)*de1_IR
      end if
    end if


    return
  end function calcdB0


  function calcB1(p2_in,m12_in,m22_in)

    use collier_coefs, only: B_cll

    complex(dp) calcB1
    complex(dp), intent(in) :: p2_in
    complex(dp), intent(in) :: m12_in
    complex(dp), intent(in) :: m22_in
    complex(dp) :: p2
    complex(dp) :: m12
    complex(dp) :: m22
    complex(dp) B1_coli

    complex(dp) B(0:1,0:1), Buv(0:1,0:1)

    complex(dp) :: rslt_b11(0:2), rslt_b00(0:2), rslt_b1(0:2), rslt_b0(0:2)


    calcB1 = 0

    p2  = p2_in
    m12 = m12_in
    m22 = m22_in


    if (se_integral_switch == 1 .or. se_integral_switch == 7) then
      p2 = real(p2)
      call B_cll(B,Buv,p2,m12,m22,1)
      calcB1 = B(1,0)
    end if



    if (se_integral_switch == 3) then
      call olo_b11(rslt_b11,rslt_b00,rslt_b1,rslt_b0,real(p2_in),m12_in,m22_in)
      if (p2 .eq. 0 .and. m12_in .eq. 0 .and. m22_in .eq. 0) then
        calcB1 = rslt_b1(0) + (de1_IR - de1_UV)/2
      else
        calcB1 = rslt_b1(0) + rslt_b1(1)*de1_UV
      end if
    end if


    return
  end function calcB1

  function calcB00(p2_in,m12_in,m22_in)

    use collier_coefs, only: B_cll

    complex(dp) calcB00
    complex(dp), intent(in) :: p2_in
    complex(dp), intent(in) :: m12_in
    complex(dp), intent(in) :: m22_in
    complex(dp) :: p2
    complex(dp) :: m12
    complex(dp) :: m22
    complex(dp) B1_coli

    complex(dp) B(0:2,0:1), Buv(0:2,0:1)

    complex(dp) :: rslt_b11(0:2), rslt_b00(0:2), rslt_b1(0:2), rslt_b0(0:2)


    calcB00 = 0

    p2  = p2_in
    m12 = m12_in
    m22 = m22_in


    if (se_integral_switch == 1 .or. se_integral_switch == 7) then
      p2 = real(p2)
      call B_cll(B,Buv,p2,m12,m22,2)
      calcB00 = B(1,0)
    end if



    if (se_integral_switch == 3) then
      call olo_b11(rslt_b11,rslt_b00,rslt_b1,rslt_b0,real(p2_in),m12_in,m22_in)
      calcB00 = rslt_b00(0) + rslt_b00(1)*de1_UV
    end if


    return
  end function calcB00

  function calcB11(p2_in,m12_in,m22_in)

    use collier_coefs, only: B_cll

    complex(dp) calcB11
    complex(dp), intent(in) :: p2_in
    complex(dp), intent(in) :: m12_in
    complex(dp), intent(in) :: m22_in
    complex(dp) :: p2
    complex(dp) :: m12
    complex(dp) :: m22
    complex(dp) B1_coli

    complex(dp) B(0:1,0:2), Buv(0:1,0:2)

    complex(dp) :: rslt_b11(0:2), rslt_b00(0:2), rslt_b1(0:2), rslt_b0(0:2)


    calcB11 = 0

    p2  = p2_in
    m12 = m12_in
    m22 = m22_in


    if (se_integral_switch == 1 .or. se_integral_switch == 7) then
      p2 = real(p2)
      call B_cll(B,Buv,p2,m12,m22,2)
      calcB11 = B(0,2)
    end if



    if (se_integral_switch == 3) then
      call olo_b11(rslt_b11,rslt_b00,rslt_b1,rslt_b0,real(p2_in),m12_in,m22_in)
      calcB11 = rslt_b11(0) + rslt_b11(1)*de1_UV
    end if


    return
  end function calcB11

  function calcdB1(p2_in,m12_in,m22_in)

    use collier_coefs, only: DB1_cll

    implicit none
    complex(dp) calcdB1
    complex(dp), intent(in) :: p2_in
    complex(dp), intent(in) :: m12_in
    complex(dp), intent(in) :: m22_in
    real(dp) :: p2q
    complex(dp) :: p2, m12, m22
    complex(dp) DB1_coli

    calcdB1 = 0

    p2  = p2_in
    m12 = m12_in
    m22 = m22_in


    if (se_integral_switch == 1 .or. se_integral_switch == 7) then
      p2 = real(p2)
      call DB1_cll(DB1_coli,p2,m12,m22)
      calcdB1 = DB1_coli
    end if



    if (se_integral_switch == 3) then
      p2q = real(p2_in)
      if (abs(p2q) > eps) then
        calcdB1 = - (m22_in-m12_in)/2/p2q**2*(calcB0(p2_in,m12_in,m22_in)-calcB0(ZERO,m12_in,m22_in)) &
       &          + (m22_in-m12_in-p2q)/2/p2q*calcdB0(p2_in,m12_in,m22_in)
      else
        if (abs(m12_in) < eps .and. abs(m22_in) < eps) then
          calcdB1 = 0
        else if (abs(m12_in) < eps .and. abs(m22_in) > eps) then
          calcdB1 = -1/m22_in/6
        else if (abs(m12_in) > eps .and. abs(m22_in) < eps) then
          calcdB1 = -1/m12_in/6
        else if (m12_in == m22_in) then
          calcdB1 = -1/m12_in/12
        else
          calcdB1 = -(2.*m12_in**3+3.*m12_in**2*m22_in-6.*m12_in*m22_in**2+m22_in**3 &
       &              +6.*m12_in**2*m22_in*log(m22_in/m12_in))/6./(m12_in-m22_in)**4
        end if
      end if
    end if


    return
  end function calcdB1


  function calcRB0(p2_in,m12_in,m22_in)
    complex(dp) calcRB0
    complex(dp), intent(in) :: p2_in, m12_in, m22_in

    calcRB0 = calcB0(p2_in,m12_in,m22_in)
    if (imag(p2_in) == 0 .and. cms_on == 1) then
      if (real(p2_in) .gt. real(m12_in+m22_in)) then
        calcRB0 = real(calcRB0)
      end if
    end if
  end function calcRB0

  function calcRdB0(p2_in,m12_in,m22_in)
    complex(dp) calcRdB0
    complex(dp), intent(in) :: p2_in, m12_in, m22_in

    calcRdB0 = calcdB0(p2_in,m12_in,m22_in)
    if (imag(p2_in) == 0 .and. cms_on == 1) then
      if (real(p2_in) .gt. real(m12_in+m22_in)) then
        calcRdB0 = real(calcRdB0)
      end if
    end if
  end function calcRdB0

  function calcRB1(p2_in,m12_in,m22_in)
    complex(dp) calcRB1
    complex(dp), intent(in) :: p2_in, m12_in, m22_in

    if (imag(p2_in) == 0 .and. cms_on == 1 &
      .and. (real(p2_in) > real(m12_in+m22_in)) &
      .and. (p2_in /= 0)) then
!(m12-m02)/2/p2*(B0 - B0(0,m0,m1)) - B0/2
          calcRB1 = (m22_in-m12_in)/2./p2_in*(calcRB0(p2_in,m12_in,m22_in)-calcB0(ZERO,m12_in,m22_in)) &
                  - 1./2.*calcRB0(p2_in,m12_in,m22_in)
    else
      calcRB1 = calcB1(p2_in,m12_in,m22_in)
    end if
  end function calcRB1

  function calcRdB1(p2_in,m12_in,m22_in)
    complex(dp) calcRdB1
    complex(dp), intent(in) :: p2_in, m12_in, m22_in

    if (imag(p2_in) == 0 .and. cms_on == 1 &
      .and. (real(p2_in) > real(m12_in+m22_in)) &
      .and. (p2_in /= 0)) then
!-(m1^2-m0^2)/2/Q4*(real(B0(Q2,m0,m1)) - B0(0,m0,m1)) + (m1^2-m0^2-Q2)/2/Q2*real(dB0(Q2,m0,m1))
          calcRdB1 = -(m22_in-m12_in)/2./p2_in**2*(calcRB0(p2_in,m12_in,m22_in)-calcB0(ZERO,m12_in,m22_in)) &
                  + (m22_in-m12_in-p2_in)/2./p2_in*calcRdB0(p2_in,m12_in,m22_in)
    else
      calcRdB1 = calcdB1(p2_in,m12_in,m22_in)
    end if

  end function calcRdB1


  function calcC0(p12_in,p22_in,p32_in,m12_in,m22_in,m32_in)

    use collier_coefs, only: C0_cll

    use avh_olo_dp

    complex(dp) calcC0
    complex(dp), intent(in) :: p12_in, p22_in, p32_in
    complex(dp), intent(in) :: m12_in, m22_in, m32_in
    complex(dp) p12, p22, p32
    complex(dp) m12, m22, m32
    complex(dp) C0_coli

    complex(dp) :: rslt(0:2)


    calcC0 = 0


    if (se_integral_switch == 1 .or. se_integral_switch == 7) then
      p12  = p12_in
      p22  = p22_in
      p32  = p32_in
      m12 = m12_in
      m22 = m22_in
      m32 = m32_in
      call C0_cll(C0_coli,p12,p22,p32,m12,m22,m32) ! check order
      calcC0 = C0_coli
    end if



    if (se_integral_switch == 3) then
      call olo_c0(rslt,p12_in,p22_in,p32_in,m12_in,m22_in,m32_in) ! check order
      calcC0 = rslt(0) + rslt(1)*de1_UV
    end if


      return
  end function calcC0

  function calcC1(p12,p02,p22,m02,m12,m22)
    complex(dp) calcC1
    complex(dp), intent(in) :: p12, p02, p22
    complex(dp), intent(in) :: m02, m12, m22
    complex(dp) :: k2, R31, R32

    ! C.34/C.36 of Denner

    k2     = kappa2(p02,p12,p22)
    R31    = 0.5*(calcB0(p22,m02,m22)-(p12-m12+m02)*calcC0(p12,p02,p22,m02,m12,m22) &
                 - calcB0(p02,m22,m12))
    R31    = 0.5*(calcB0(p12,m02,m12)-(p22-m22+m02)*calcC0(p12,p02,p22,m02,m12,m22) &
                 - calcB0(p02,m22,m12))
    calcC1 = -4./k2*(p22*R31+0.5*(p02-p12-p22)*R32)

    return
  end function calcC1

  function calcC2(p12,p02,p22,m02,m12,m22)
    complex(dp) calcC2
    complex(dp), intent(in) :: p12, p02, p22
    complex(dp), intent(in) :: m02, m12, m22
    complex(dp) :: k2, R31, R32

    ! C.34/C.36 of Denner

    k2     = kappa2(p02,p12,p22)
    R31    = 0.5*(calcB0(p22,m02,m22)-(p12-m12+m02)*calcC0(p12,p02,p22,m02,m12,m22) &
                 - calcB0(p02,m22,m12))
    R31    = 0.5*(calcB0(p12,m02,m12)-(p22-m22+m02)*calcC0(p12,p02,p22,m02,m12,m22) &
                 - calcB0(p02,m22,m12))
    calcC2 = -4./k2*(0.5*(p02-p12-p22)*R31+p12*R32)

    return
  end function calcC2

  function kappa2(x,y,z)
    complex(dp) kappa2
    complex(dp), intent(in) :: x,y,z
    kappa2 = x*x+y*y+z*z-2.*(x*y+y*z+z*x)
  end function kappa2

end module ol_self_energy_integrals_dp


