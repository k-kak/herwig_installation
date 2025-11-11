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


module ol_kinematics_qp
  use kind_types
  implicit none

  interface internal_momenta
    module procedure internal_momenta_six, internal_momenta_std
  end interface

  interface get_mass
    module procedure get_mass_id,get_mass_idarr
  end interface

  interface get_mass2
    module procedure get_mass2_id,get_mass2_idarr
  end interface

  interface get_rmass2
    module procedure get_rmass2_id,get_rmass2_idarr
  end interface

  interface conv_mom_scatt2in
    module procedure conv_mom_scatt2in_mexpl, conv_mom_scatt2in_mids
  end interface

  interface init_kinematics
    module procedure init_kinematics_mexpl, init_kinematics_mids
  end interface


  real(qp) :: collthres = 1.E-4
  real(qp) :: softthres = 1.E-4

contains

! **********************************************************************
subroutine Std2LC_Rep(P,L)
! **********************************************************************
! Lorentz -> light-cone representation
! P(0:3)   = Lorentz momentum P^mu (contravariant)
! L(1:4)   = light-cone representation L^A (contravariant)
! **********************************************************************
  use kind_types, only: qp
  use ol_parameters_decl_qp, only: CI
  implicit none
  real(qp),    intent(in)  :: P(0:3)
  complex(qp), intent(out) :: L(1:4)
  L(1) =  P(0) - P(3)
  L(2) =  P(0) + P(3)
  L(3) = -P(1) - CI*P(2)
  L(4) = -P(1) + CI*P(2)
end subroutine Std2LC_Rep


! **********************************************************************
subroutine Std2LC_cmplx(P,L)
! **********************************************************************
! complex version needed for OPP Reduction
! Lorentz -> light-cone representation
! P(0:3)   = Lorentz momentum P^mu (ALWAYS CONTRAVARIANT)
! L(1:4)   = light-cone representation L^A (ALWAYS CONTRAVARIANT)
! **********************************************************************
  use kind_types, only: qp
  use ol_parameters_decl_qp, only: CI
  implicit none
  complex(qp), intent(in)  :: P(0:3)
  complex(qp), intent(out) :: L(1:4)
  L(1) =  P(0) -    P(3)
  L(2) =  P(0) +    P(3)
  L(3) = -P(1) - CI*P(2)
  L(4) = -P(1) + CI*P(2)
end subroutine Std2LC_cmplx


! **********************************************************************
function get_LC_5(mom) result(P)
! **********************************************************************
! P(1:4) = light-cone representation
! P(5)   = mass
! **********************************************************************
  use kind_types, only: qp
  use ol_momenta_decl_qp, only: L
  implicit none
  integer, intent(in) :: mom
  complex(qp) :: P(1:5)
  if (mom .gt. 0) then
    P(1:4) = L(1:4,mom)
    P(5) = L(5,mom) + L(6,mom)
  else
    P(1:4) = -L(1:4,-mom)
    P(5) = L(5,-mom) + L(6,-mom)
  end if

end function get_LC_5





! **********************************************************************
function get_LC_4(mom) result(P)
! **********************************************************************
! P(1:4) = light-cone representation
! **********************************************************************
  use kind_types, only: qp
  use ol_momenta_decl_qp, only: L
  implicit none
  integer, intent(in) :: mom
  complex(qp) :: P(1:4)
  if (mom .gt. 0) then
    P(1:4) = L(1:4,mom)
  else
    P(1:4) = -L(1:4,-mom)
  end if

end function get_LC_4




! **********************************************************************
function get_LC_mass2(mom) result(m2)
! **********************************************************************
! m2   = mass squared
! **********************************************************************
  use kind_types, only: qp
  use ol_momenta_decl_qp, only: L
  implicit none
  integer, intent(in) :: mom
  complex(qp) :: m2
  if (mom .gt. 0) then
    m2 = L(5,mom) + L(6,mom)
  else
    m2 = L(5,-mom) + L(6,-mom)
  end if
end function get_LC_mass2


! **********************************************************************
function get_mass_id(mid) result(m)
! **********************************************************************
! m   = complex mass
! **********************************************************************
  use kind_types, only: qp
  use ol_parameters_decl_qp
  use ol_generic, only: to_string
  use ol_debug, only: ol_error, ol_fatal
  implicit none
  integer, intent(in) :: mid
  complex(qp) :: m

  select case (mid)
  case (0)
    m = ZERO
  case (nME)
    m = ME
  case (nMM)
    m = MM
  case (nML)
    m = ML
  case (nMU)
    m = MU
  case (nMD)
    m = MD
  case (nMS)
    m = MS
  case (nMC)
    m = MC
  case (nMB)
    m = MB
  case (nMT)
    m = MT
  case (nMW)
    m = MW
  case (nMZ)
    m = MZ
  case (nMH)
    m = MH
  case (nMX)
    m = MX
  case (nMY)
    m = MY
  case default
    call ol_error(2,"Unknown mass id: " // to_string(mid))
    call ol_fatal()
  end select
end function get_mass_id


! **********************************************************************
function get_mass_idarr(mids) result(m)
! **********************************************************************
! m   = complex mass array
! **********************************************************************
  use kind_types, only: qp
  use ol_parameters_decl_qp
  implicit none
  integer, dimension(:), intent(in) :: mids
  complex(qp), dimension(size(mids)) :: m
  integer :: i

  do i = 1, size(mids)
    m(i) = get_mass_id(mids(i))
  end do
end function get_mass_idarr


! **********************************************************************
function get_mass2_id(mid) result(m2)
! **********************************************************************
! m2   = complex mass squared
! **********************************************************************
  use kind_types, only: qp
  use ol_parameters_decl_qp
  use ol_generic, only: to_string
  use ol_debug, only: ol_error, ol_fatal
  implicit none
  integer, intent(in) :: mid
  complex(qp) :: m2

  select case (mid)
  case (0)
    m2 = ZERO
  case (nME)
    m2 = ME2
  case (nMM)
    m2 = MM2
  case (nML)
    m2 = ML2
  case (nMU)
    m2 = MU2
  case (nMD)
    m2 = MD2
  case (nMS)
    m2 = MS2
  case (nMC)
    m2 = MC2
  case (nMB)
    m2 = MB2
  case (nMT)
    m2 = MT2
  case (nMW)
    m2 = MW2
  case (nMZ)
    m2 = MZ2
  case (nMH)
    m2 = MH2
  case (nMX)
    m2 = MX2
  case (nMY)
    m2 = MY2
  case default
    call ol_error(2,"Unknown mass id: " // to_string(mid))
    call ol_fatal()
  end select
end function get_mass2_id


! **********************************************************************
function get_rmass2_id(mid) result(m2)
! **********************************************************************
! m2   = real mass squared
! **********************************************************************
  use kind_types, only: qp
  use ol_parameters_decl_qp
  use ol_generic, only: to_string
  use ol_debug, only: ol_error, ol_fatal
  implicit none
  integer, intent(in) :: mid
  real(qp) :: m2

  select case (mid)
  case (0)
    m2 = ZERO
  case (nME)
    m2 = rME2
  case (nMM)
    m2 = rMM2
  case (nML)
    m2 = rML2
  case (nMU)
    m2 = rMU2
  case (nMD)
    m2 = rMD2
  case (nMS)
    m2 = rMS2
  case (nMC)
    m2 = rMC2
  case (nMB)
    m2 = rMB2
  case (nMT)
    m2 = rMT2
  case (nMW)
    m2 = rMW2
  case (nMZ)
    m2 = rMZ2
  case (nMH)
    m2 = rMH2
  case (nMX)
    m2 = rMX2
  case (nMY)
    m2 = rMY2
  case default
    call ol_error(2,"Unknown mass id: " // to_string(mid))
    call ol_fatal()
  end select
end function get_rmass2_id


! **********************************************************************
function get_mass2_idarr(mids) result(m2)
! **********************************************************************
! m2   = complex mass squared array
! **********************************************************************
  use kind_types, only: qp
  use ol_parameters_decl_qp
  integer, dimension(:),        intent(in) :: mids
  complex(qp), dimension(size(mids)) :: m2
  integer :: i

  do i = 1, size(mids)
    m2(i) = get_mass2_id(mids(i))
  end do
end function get_mass2_idarr


! **********************************************************************
function get_rmass2_idarr(mids) result(m2)
! **********************************************************************
! m2   = real mass squared array
! **********************************************************************
  use kind_types, only: qp
  use ol_parameters_decl_qp
  integer, dimension(:),     intent(in) :: mids
  real(qp), dimension(size(mids)) :: m2
  integer :: i

  do i = 1, size(mids)
    m2(i) = get_rmass2_id(mids(i))
  end do
end function get_rmass2_idarr


! **********************************************************************
subroutine LC2Std_Rep(L,P)
! **********************************************************************
! light-cone -> Lorentz representation
! L(1:4)   = light-cone representation L^A (ALWAYS CONTRAVARIANT)
! P(0:3)   = Lorentz momentum P^mu (ALWAYS CONTRAVARIANT)
! **********************************************************************
  use kind_types, only: qp
  implicit none
  complex(qp), intent(in)  :: L(1:4)
  real(qp),    intent(out) :: P(0:3)
  P(0) =  real(L(1)+L(2))*0.5_qp
  P(1) = -real(L(3)+L(4))*0.5_qp
  P(2) = aimag(L(4)-L(3))*0.5_qp
  P(3) =  real(L(2)-L(1))*0.5_qp
end subroutine LC2Std_Rep


! **********************************************************************
subroutine LC2Std_Rep_D(L,P)
! **********************************************************************
! light-cone -> Lorentz representation
! L(1:4)   = light-cone representation L^A (ALWAYS CONTRAVARIANT)
! P(0:3)   = Lorentz momentum P^mu (ALWAYS CONTRAVARIANT)
! P(0:3)   here is a double-realkind vector
! **********************************************************************
  use kind_types, only: qp
  implicit none
  complex(qp), intent(in)  :: L(1:4)
  real(dp),    intent(out) :: P(0:3)
  P(0) =  real(L(1)+L(2))*0.5_qp
  P(1) = -real(L(3)+L(4))*0.5_qp
  P(2) = aimag(L(4)-L(3))*0.5_qp
  P(3) =  real(L(2)-L(1))*0.5_qp
end subroutine LC2Std_Rep_D


! **********************************************************************
subroutine LC2Std_Rep_cmplx(L,P)
! **********************************************************************
! complex version
! light-cone -> Lorentz representation
! L(1:4)   = light-cone representation L^A (ALWAYS CONTRAVARIANT)
! P(0:3)   = Lorentz momentum P^mu (ALWAYS CONTRAVARIANT)
! **********************************************************************
  use kind_types, only: qp
  use ol_parameters_decl_qp, only: CI
  implicit none
  complex(qp), intent(in)  :: L(1:4)
  complex(qp), intent(out) :: P(0:3)
  P(0) =  (L(1)+L(2))*0.5_qp
  P(1) = -(L(3)+L(4))*0.5_qp
  P(2) = -CI*(L(4)-L(3))*0.5_qp
  P(3) =  (L(2)-L(1))*0.5_qp
end subroutine LC2Std_Rep_cmplx


! **********************************************************************
function cont_L_cmplx(A)
! Contraction of a complex Lorentz vector in standard representation with itself
! **********************************************************************
  use kind_types, only: qp
  implicit none
  complex(qp) :: cont_L_cmplx
  complex(qp), intent(in) :: A(0:3)
  cont_L_cmplx = A(0)*A(0) - A(1)*A(1) - A(2)*A(2) - A(3)*A(3)
end function cont_L_cmplx


! **************************************************************************
function cont_LC_cntrv(V1,V2)
! Contraction of contravariant Lorentz vectors in LightCone representation.
! Contraction of V1 with itself or with a second vector V2
! **************************************************************************
  use kind_types, only: qp
  implicit none
  complex(qp) :: cont_LC_cntrv
  complex(qp), intent(in) :: V1(1:4)
  complex(qp), optional, intent(in) :: V2(1:4)
  if(present(V2)) then
    cont_LC_cntrv = (V1(1)*V2(2) + V2(1)*V1(2) - V1(3)*V2(4) - V2(3)*V1(4))/2
  else
    cont_LC_cntrv = V1(1)*V1(2) - V1(3)*V1(4)
  end if
end function cont_LC_cntrv





! **********************************************************************
subroutine clean_mom_in(P_in, m_ext2, P, n)
! remove numerical inaccuracies in n -> 0 point (momenta all incoming)
! P_in(0:3,n) = original PS point in double precision
! m_ext2      = squared external masses
! P(0:3,n)    = improved PS point in working precision
! 3-momentum of n-th particle recomputed with momentum conservation;
! all energies recomputed with on-shell condition;
! 3-momenta (and energies) of incoming particles shifted by +/- eps*p1,
! so that energy conservation is fulfilled up to terms of O(eps^3)
! **********************************************************************
  use kind_types, only: qp
  use ol_debug, only: ol_msg
  use ol_generic, only: to_string
  use ol_parameters_decl_dp, only: psp_tolerance, no_cleaning, cleaning_via_hardness
  implicit none
  real(qp), intent(in)  :: P_in(0:3,n)
  integer,        intent(in)  :: n
  real(qp), intent(in)  :: m_ext2(n)
  real(qp), intent(out) :: P(0:3,n)
  real(qp)  :: E_ref, P0(n), P2(n)
  real(qp)  :: softness, collinearity, hardness(3:n)
  real(qp)  :: E0(2), E1(2), E2(2)
  real(qp)  :: E0_tot, E1_tot, E2_tot
  real(qp)  :: eps1, eps
  integer         :: nex, i, j, pout_max_pos, pout_max_pos_arr(1)
  real            :: prec

  P = P_in

!   call dirty_mom(P_in, P, n, 9)

  ! print momena before cleaning
!   do i = 1, n
!     write(*,*) P(:,i), (P(0,i)**2-P(1,i)**2-P(2,i)**2-P(3,i)**2)-m_ext2(i)
!   end do
!   write(*,*) sum(P(0,:)), sum(P(1,:)), sum(P(2,:)), sum(P(3,:))
!   write(*,*)

  E_ref = 0.5_qp * sum(abs(P(0,:)))
  ! check momentum conservation
  do i = 0, 3
    prec = abs(sum(P(i,:)))/E_ref
    if (prec > psp_tolerance) then
      call ol_msg("=== WARNING ===")
      call ol_msg("OpenLoops subroutine clean_mom: inconsistent phase space point.")
      call ol_msg("Momentum conservation is only satisfied to " // to_string(-log10(prec)) // "digits.")
      call ol_msg("===============")
    end if
  end do

  ! check on-shell conditions
  do nex = 1, n
    P2(nex) = P(1,nex)*P(1,nex) + P(2,nex)*P(2,nex) + P(3,nex)*P(3,nex)
    P0(nex) = sign(sqrt(P2(nex) + m_ext2(nex)), P(0,nex))
    prec = abs(P(0,nex)-P0(nex))/E_ref
    if(prec > psp_tolerance) then
      call ol_msg("=== WARNING ===")
      call ol_msg("OpenLoops subroutine clean_mom: inconsistent phase space point.")
      call ol_msg("On-shell condition is only satisfied to " // to_string(-log10(prec)) // "digits.")
      call ol_msg("===============")
    end if
  end do

  if (no_cleaning) return

  if (cleaning_via_hardness) then
    ! determine position of the hardest momentum. hardness h is defined as
    ! h_i = min(softness_i, collinearity_i1, collinearity_i2, ... ) with
    !   softness_i = E_i/E_ref
    !   collinearity_ij = p_i . p_j / (E_i E_j)  i != j
    do i = 3, n
      softness = abs(P(0,i)/E_ref)
      hardness(i) = softness
    end do
    do i = 3, n
      do j = 1, i-1
        collinearity = abs((P(0,i)*P(0,j) - P(1,i)*P(1,j) - P(2,i)*P(2,j) - P(3,i)*P(3,j))/(P(0,i)*P(0,j)))
        hardness(i) = min(hardness(i), collinearity)
        if (j .gt. 2) then
          hardness(j) = min(hardness(j), collinearity)
        end if
      end do
    end do
    pout_max_pos_arr = maxloc(hardness(3:))
  else
    ! position of the outgoing momentum with the largest energy
    pout_max_pos_arr = maxloc(abs(P(0,3:)))
  end if
  pout_max_pos = 2 + pout_max_pos_arr(1)
  ! fix 3-momentum by momentum conservation
  do i = 1, 3
    P(i,pout_max_pos) = P(i,pout_max_pos) - sum(P(i,:))
  end do

  ! enforce on-shell conditions
  P0(pout_max_pos) = sign(sqrt(P(1,pout_max_pos)**2 + P(2,pout_max_pos)**2 + P(3,pout_max_pos)**2 + m_ext2(pout_max_pos)), &
                        & P(0,pout_max_pos))
  P(0,:) = P0

  E0_tot = sum(P(0,:))

  ! lowest-order energy terms
  E0(1)  = P(0,1)
  E0(2)  = P(0,2)

  ! 1st order energy coefficients
  E1(1)  = P2(1)/E0(1)
  E1(2)  = -(P(1,1)*P(1,2)+P(2,1)*P(2,2)+P(3,1)*P(3,2))/E0(2)
  E1_tot = E1(1)+E1(2)

  ! 2nd order energy coefficients
  ! for quad-precision applications it is recommended to use the formulas w.o.
  ! beam-alignement instabilities both for E2(1) and E2(2)
!   E2(1)  = (P2(1)-E1(1)**2)/(2*E0(1))   ! direct formula
  E2(1)  = P2(1)*m_ext2(1)/(2*E0(1)**3) ! default = equivalent fomula w.o. beam-alignment instabilities
  E2(2)  = (P2(1)-E1(2)**2)/(2*E0(2))   ! default = direct formula
!   E2(2)  = 0.5_qp/E0(2)**3*(P2(1)*(m_ext2(2) + P(1,2)**2 + P(2,2)**2) & ! equivalent formula w.o. beam-alignment instabilities
!          + (P(1,1)**2+P(2,1)**2)*P(3,2)**2 &
!          - (P(1,1)*P(1,2)+P(2,1)*P(2,2))*(P(1,1)*P(1,2)+P(2,1)*P(2,2)+2*P(3,1)*P(3,2)))
  E2_tot = E2(1) + E2(2)

  ! 1st order shift
  eps1 = -E0_tot/E1_tot
  ! 2nd order shift
  eps  = eps1*(1-eps1*E2_tot/E1_tot)

  ! shifted IS momenta
  do i = 1, 3
    P(i,2) = P(i,2) - P(i,1)*eps
    P(i,1) = P(i,1) + P(i,1)*eps
  end do

  ! shifted IS energies
  do nex = 1, 2
    if(m_ext2(nex) == 0 .and. P(1,nex) == 0 .and. P(2,nex) == 0) then
      ! exact formula for m = 0 along beam-axis
      P(0,nex) = sign(abs(P(3,nex)), real(P_in(0,nex), qp))
    else
      P(0,nex) = E0(nex) + E1(nex)*eps + E2(nex)*eps**2
    end if
  end do

  ! print momenta after cleaning
!    do i = 1, n
!      write(*,*) P(:,i), (P(0,i)**2-P(1,i)**2-P(2,i)**2-P(3,i)**2)-m_ext2(i)
!    end do
!    write(*,*) sum(P(0,:)), sum(P(1,:)), sum(P(2,:)), sum(P(3,:))
!   write(*,*)

end subroutine clean_mom_in


! **********************************************************************
subroutine clean_mom_scatt(P_in, m_ext2, P, n)
! same as clean_mom_in but for 2-> n-2 PS point
! This routine is not used internally.
! **********************************************************************
  use kind_types, only: qp
  implicit none
  real(qp), intent(in)  :: P_in(0:3,n)
  integer,        intent(in)  :: n
  real(qp), intent(in)  :: m_ext2(n)
  real(qp), intent(out) :: P(0:3,n)
  real(qp) :: Q_in(0:3,n)
  real(qp) :: Q(0:3,n)

  Q_in(:,1:2) = P_in(:,1:2)
  Q_in(:,3:) = -P_in(:,3:)

  call clean_mom_in(Q_in, m_ext2, Q, n)

  P(:,1:2) = Q(:,1:2)
  P(:,3:) = -Q(:,3:)

end subroutine clean_mom_scatt



! **********************************************************************
subroutine dirty_mom(P_in, P, n, DIG)
! introduces random noise in every PS component after DIG digits
! (to test PS-point cleaning).
! Needs RAMBO subroutine rans().
! **********************************************************************
  use kind_types, only: qp, dp
  use ol_ramboX, only: rans
  implicit none
  integer,        intent(in)  :: n, DIG
  real(qp), intent(in)  :: P_in(0:3,n)
  real(qp), intent(out) :: P(0:3,n)
  real(dp) :: x, shift
  integer        :: nex, i

  shift = 10._dp**(-DIG)

  do nex = 1, n
    do i = 0, 3
      call rans(x)
      P(i,nex) = P_in(i,nex)*(1+(x-0.5_qp)*shift)
    end do
  end do

end subroutine dirty_mom
! #ifdef 1


subroutine init_kinematics_mexpl(P_scatt, m_ext2, P_in_clean, perm_inv, n)
  use ol_external_decl_dp, only: init_qp
  integer,           intent(in)  :: n
  real(dp),   intent(in)  :: P_scatt(0:3,n)
  real(qp),    intent(in)  :: m_ext2(n)
  real(qp),    intent(out) :: P_in_clean(0:3,n)
  integer,           intent(in)  :: perm_inv(n)
  real(qp) :: P(0:3,n)


  init_qp = .true.


  call conv_mom_scatt2in(P_scatt, m_ext2, P, perm_inv, n)
  call internal_momenta(P, n)

end subroutine init_kinematics_mexpl

subroutine init_kinematics_mids(P_scatt, m_ext2, P_in_clean, perm_inv, n, qp_kinematics)

  use ol_external_decl_dp, only: init_qp
  integer,         intent(in)  :: n
  real(dp), intent(in)  :: P_scatt(0:3,n)
  integer,         intent(in)  :: m_ext2(n)
  real(qp),  intent(inout) :: P_in_clean(0:3,n)
  integer,         intent(in)  :: perm_inv(n)
  logical,         intent(in)  :: qp_kinematics

  init_qp = .true.


  call conv_mom_scatt2in_mids(P_scatt, m_ext2, P_in_clean, perm_inv, n)

  call internal_momenta_six(P_in_clean, n, m_ext2, .false.)


end subroutine init_kinematics_mids



! **********************************************************************
subroutine conv_mom_scatt2in_mexpl(P_scatt, m_ext2, P_in_clean, perm_inv, n)
! Keep incoming momenta and reverse outgoing momenta.
! Apply phase space cleaning and crossing.
! Note: calls init_kin_arrays -> allocation of kinematic arrays
! ToDo: cleaning for n_scatt /= 2
! **********************************************************************
  use kind_types, only: qp, dp
  use ol_external_decl_qp, only: nParticles, P_ex, inverse_crossing
  use ol_external_decl_dp, only: n_scatt
  use ol_parameters_decl_qp, only: scalefactor
  use ol_parameters_init_qp, only: init_kin_arrays
!#ifdef PRECISION_dp
!  use ol_parameters_init_qp, only: init_kin_arrays_qq => init_kin_arrays
!#endif
  integer,           intent(in)  :: n
  real(dp),   intent(in)  :: P_scatt(0:3,n)
  real(qp),    intent(in)  :: m_ext2(n)
  real(qp),    intent(out) :: P_in_clean(0:3,n)
  integer,           intent(in)  :: perm_inv(n)
  real(qp) :: P_in(0:3,n)
  real(qp) :: P_clean(0:3,n), m_ext2_perm(n)
  integer        :: k
  nParticles = n
  call init_kin_arrays(nParticles)
  P_ex(:,1:n) = P_scatt
  inverse_crossing(1:n) = perm_inv
  do k = 1, n
    m_ext2_perm(perm_inv(k)) = m_ext2(k)
  end do
  P_in(:,1:n_scatt) =   scalefactor * P_scatt(:,1:n_scatt)
  P_in(:,n_scatt+1:n)  = - scalefactor * P_scatt(:,n_scatt+1:n)
  if (n_scatt == 2 .and. n > 2) then
    ! Clean momenta to get full numerical precision.
    ! Do the cleaning in the original permutation where the first two momenta are incoming.
    ! Otherwise the beam alignment (zero components) might be spoiled by the cleaning.
    call clean_mom_in(P_in, m_ext2_perm, P_clean, n)
  else
    P_clean = P_in
  end if
  do k = 1, n
    P_in_clean(:,k) = P_clean(:,perm_inv(k))
  end do
end subroutine conv_mom_scatt2in_mexpl

subroutine conv_mom_scatt2in_mids(P_scatt, m_ext, P_in_clean, perm_inv, n)
  use kind_types, only: qp, dp
  use ol_external_decl_qp, only: M_ex
  use ol_external_decl_dp, only: n_scatt
  use ol_parameters_decl_qp, only: scalefactor
  use ol_parameters_init_qp, only: init_kin_arrays
!#ifdef PRECISION_dp
!  use ol_parameters_init_qp, only: init_kin_arrays_qq => init_kin_arrays
!#endif
  implicit none
  integer,           intent(in)  :: n
  real(dp),   intent(in)  :: P_scatt(0:3,n)
  integer,           intent(in)  :: m_ext(n)
  real(qp),    intent(out) :: P_in_clean(0:3,n)
  integer,           intent(in)  :: perm_inv(n)

  call conv_mom_scatt2in_mexpl(P_scatt,get_rmass2(m_ext),P_in_clean,perm_inv,n)
  M_ex(1:n) = m_ext
end subroutine conv_mom_scatt2in_mids


! **********************************************************************
subroutine conv_mom_scatt2in_cache(P_in_clean,n)
! ToDo: Initializes the momenta in QP from the DP cache
! **********************************************************************
  use kind_types, only: qp, dp
  use ol_external_decl_dp, only: nParticles, P_ex_dp=>P_ex, M_ex_dp=>M_ex, inverse_crossing_dp=>inverse_crossing
  use ol_external_decl_qp, only: P_ex
  use ol_external_decl_dp, only: n_scatt
  use ol_parameters_decl_qp, only: scalefactor
  use ol_parameters_init_qp, only: init_kin_arrays
  integer,           intent(in)  :: n
  real(qp),    intent(out) :: P_in_clean(0:3,n)
  real(qp) :: P_in(0:3,n)
  real(qp) :: P_clean(0:3,n), m_ext2_perm(n)
  integer        :: k
  call init_kin_arrays(n)
  P_ex(:,1:n) = P_ex_dp(:,1:n)
  do k = 1, n
    m_ext2_perm(inverse_crossing_dp(k)) = get_rmass2(M_ex_dp(k))
  end do
  P_in(:,1:n_scatt) =   scalefactor * P_ex(:,1:n_scatt)
  P_in(:,n_scatt+1:n)  = - scalefactor * P_ex(:,n_scatt+1:n)
  if (n_scatt == 2 .and. n > 2) then
    ! Clean momenta to get full numerical precision.
    ! Do the cleaning in the original permutation where the first two momenta are incoming.
    ! Otherwise the beam alignment (zero components) might be spoiled by the cleaning.
    call clean_mom_in(P_in, m_ext2_perm, P_clean, n)
  else
    P_clean = P_in
  end if
  do k = 1, n
    P_in_clean(:,k) = P_clean(:,inverse_crossing_dp(k))
  end do
end subroutine conv_mom_scatt2in_cache



! **********************************************************************
subroutine conv_mom_os(P_decay, P_in, n)
! reverse decay products
! ToDo: Apply phase space cleaning
! **********************************************************************
  use kind_types, only: qp, dp
  use ol_parameters_decl_dp, only: scalefactor
  implicit none
  integer,         intent(in)  :: n
  real(dp), intent(in)  :: P_decay(0:3,n)
  real(qp),  intent(out) :: P_in(0:3,n)
  integer         :: k

  P_in  = - scalefactor * P_decay

end subroutine conv_mom_os







! **********************************************************************
subroutine internal_momenta_std(P, Npart)
! **********************************************************************
! P(0:3,Npart) = external real-valued four-momenta (standard representation)
! Npart        = total (in & out) external particle number
! Q(1:5,1:2^Npart-1) = internal four-momenta in light-cone representation;
!                      the fifth component is the complex-valued squared momentum.
! Numbering of internal momenta:
!   Sum_i s(i)*P(i) => Q(Sum_i s(i)*2^(i-1)), s(i) = 0, 1
!   so that Q(J1) + Q(J2) = Q(J1+J2)
! QInvariantsMatrix(i,j) = (p_i+p_j)^2 for i /= j, otherwise undefined.
! **********************************************************************
  use kind_types, only: qp
  use ol_momenta_decl_qp, only: Q, QInvariantsMatrix, &
                                          collconf, softconf
  implicit none

  real(qp), intent(in) :: P(0:3,Npart)
  integer,        intent(in)  :: Npart
  integer :: i, j
  integer :: Jmax
  integer :: i1, i2 ! conventional particle numbers
  integer :: l1, l2 ! individual 2^(i-1) particles numbers
  integer :: s1, s2 ! sums of 2^(i-1) particle numbers
  integer :: r1, r2 ! inverse of s1, s2, ...

  collconf = 0
  softconf = 0
  i = 2**Npart - 2

  Q(:,i+1) = 0

  if (Npart == 2) then
    call Std2LC_Rep(P(:,1),Q(1:4,1))
    Q(1:4,2) = - Q(1:4,1)
    Q(5,1) = Q(1,1)*Q(2,1) - Q(3,1)*Q(4,1)
    Q(5,2) = Q(5,1)
  else if (Npart == 3) then
    call Std2LC_Rep(P(:,1),Q(1:4,1))
    call Std2LC_Rep(P(:,2),Q(1:4,2))
    Q(1:4,3) = Q(1:4,1) + Q(1:4,2)
    Q(1:4,4) = - Q(1:4,3)
    Q(1:4,5) = - Q(1:4,2)
    Q(1:4,6) = - Q(1:4,1)
    Q(5,1) = Q(1,1)*Q(2,1) - Q(3,1)*Q(4,1)
    Q(5,2) = Q(1,2)*Q(2,2) - Q(3,2)*Q(4,2)
    Q(5,3) = Q(1,3)*Q(2,3) - Q(3,3)*Q(4,3)
    Q(5,4) = Q(5,3)
    Q(5,5) = Q(5,2)
    Q(5,6) = Q(5,1)
  else
    call intmom(P, Npart, i)
  end if

  do i = 1, Npart
    ! QInvariantsMatrix(i,i) = m_ex2(i)
    do j = i + 1, Npart
      QInvariantsMatrix(i,j) = Q(5,2**(i-1)+2**(j-1))
      QInvariantsMatrix(j,i) = QInvariantsMatrix(i,j)
    end do
  end do

end subroutine internal_momenta_std


! **********************************************************************
subroutine internal_momenta_six(P, Npart, ext_masses, use_qp_kinematics)
! **********************************************************************
! P(0:3,Npart)       = external real-valued four-momenta
!                      (standard-representation)
! Npart              = number of external particle (in & out)
! L(1:6,1:2^Npart-1) = internal four-momenta in light-cone representation;
!                      the fifth component is the real-valued sum of external
!                      masses while the
!                      the sixth component is the sum of the scalar-products
! Numbering of internal momenta:
!   Sum_i s(i)*P(i) => Q(Sum_i s(i)*2^(i-1)), s(i) = 0, 1
!   so that Q(J1) + Q(J2) = Q(J1+J2)
! QInvariantsMatrix(i,j) = (p_i+p_j)^2 for i /= j, otherwise undefined.
! invariants_mode    = flag used to decide in which way the invariants
!                      from the internal momenta are computed
! **********************************************************************
  use kind_types, only: qp, intkind1
  use ol_momenta_decl_qp, only: L, collconf, softconf

  implicit none

  integer,        intent(in) :: Npart
  real(qp), intent(in) :: P(0:3,Npart)
  integer,        intent(in) :: ext_masses(Npart)
  logical,        intent(in) :: use_qp_kinematics
  complex(qp) :: zero = 0._qp
  real(qp)    :: maxinv
  integer :: i, j
  integer :: Jmax
  integer :: i1, i2 ! conventional particle numbers
  integer :: l1, l2 ! individual 2^(i-1) particles numbers
  integer :: s1, s2 ! sums of 2^(i-1) particle numbers
  integer :: r1, r2 ! inverse of s1, s2, ...
  integer :: t1

  i = 2**Npart - 2
  L(:,i+1) = 0
  L(5:6,:) = zero
  collconf = 0
  softconf = 0

  if (Npart == 2) then

    call Std2LC_Rep(P(:,1),L(1:4,1))

    L(5,1) = get_rmass2(ext_masses(1))
    L(6,1) = zero
    L(5:6,2) = L(5:6,1)
  else if (Npart == 3) then

    call Std2LC_Rep(P(:,1),L(1:4,1))
    call Std2LC_Rep(P(:,2),L(1:4,2))
    L(1:4,3) = L(1:4,1) + L(1:4,2)
    L(1:4,4) = - L(1:4,3)
    L(1:4,5) = - L(1:4,2)
    L(1:4,6) = - L(1:4,1)

    L(5,1) = get_rmass2(ext_masses(1))
    L(5,2) = get_rmass2(ext_masses(2))
    L(6,1) = zero
    L(6,2) = zero
    L(5,3) = get_rmass2(ext_masses(1)) + get_rmass2(ext_masses(2))

    L(6,3) = 2*cont_LC_cntrv(L(1:4,1),L(1:4,2))

    L(5:6,4) = L(5:6,3)
    L(5:6,5) = L(5:6,2)
    L(5:6,6) = L(5:6,1)
  else
    ! compute masses and scalar products
    call intmom_six(P, Npart, i, get_rmass2(ext_masses), use_qp_kinematics)

  end if



end subroutine internal_momenta_six


! **********************************************************************
subroutine intmom(P_ex,Npart,Jmax)
! **********************************************************************
  use kind_types, only: qp
  use ol_momenta_decl_qp, only: Q
  implicit none

  integer,        intent(in) :: Npart, Jmax
  real(qp), intent(in) :: P_ex(0:3,Npart)
  integer :: A
  integer :: i1 ! conventional particle numbers
  integer :: l1 ! individual 2^(i-1) particles numbers
  integer :: s1 ! sums of 2^(i-1) particle numbers
  integer :: r1 ! inverse of s1, s2, ...

  l1 = 1 ! ext mom 1 <= i1 <= Npart
  do i1 = 1, Npart
    s1 = l1
    r1 = Jmax + 1 - s1
    call Std2LC_Rep(P_ex(0,i1),Q(1,l1))
    l1 = 2*l1
    do A = 1, 4
      Q(A,r1) = -Q(A,s1)
    end do
    call intmom_rec(Npart, Jmax, i1, s1, 1)
  end do !i1

  do s1 = 1, Jmax/2 ! squared momenta
    Q(5,s1) = Q(1,s1)*Q(2,s1) - Q(3,s1)*Q(4,s1)
    Q(5,Jmax+1-s1) = Q(5,s1)
  end do

end subroutine intmom


! **********************************************************************
subroutine intmom_six(P_ex,Npart,Jmax,ext_masses,use_qp_kinematics)
! **********************************************************************
  use kind_types, only: qp
  use ol_momenta_decl_qp, only: L
  use ol_loop_parameters_decl_qp, only: zero
  implicit none

  integer,        intent(in) :: Npart, Jmax
  real(qp), intent(in) :: ext_masses(Npart)
  real(qp), intent(in) :: P_ex(0:3,Npart)
  logical, intent(in)        :: use_qp_kinematics
  integer :: i1 ! conventional particle numbers
  integer :: l1 ! individual 2^(i-1) particles numbers

  l1 = 1
  do i1 = 1, Npart
    L(5,l1) = ext_masses(i1)
    L(6,l1) = zero
    L(5:6,Jmax + 1 - l1) = L(5:6,l1)
    call Std2LC_Rep(P_ex(0,i1),L(1,l1))
    L(1:4,Jmax + 1 - l1) = -L(1:4,l1)
    call intmom_rec_six(Npart, Jmax, i1, l1, 1, use_qp_kinematics)
    l1 = 2*l1
  end do !i1

end subroutine intmom_six


! **********************************************************************
recursive subroutine intmom_rec(Npart, Jmax, i1, s1, x)
! **********************************************************************
  use ol_momenta_decl_qp, only: Q
  implicit none
  integer,        intent(in) :: Npart, Jmax, i1, s1, x
  integer :: A
  integer :: ix ! conventional particle numbers
  integer :: lx ! individual 2^(i-1) particles numbers
  integer :: sx ! sums of 2^(i-1) particle numbers
  integer :: rx ! inverse of sx
  logical :: last

  last = .false.
  if (2*x+2 == Npart .or. 2*x+3 == Npart) then
    last = .true.
  end if

  lx = 1 ! adding ext mom 1 <= ix < i1
  do ix = 1, i1 - 1
    sx = s1 + lx
    rx = Jmax + 1 - sx
    if ( (last .eqv. .false.) .or. (mod(Npart,2) == 1 .or. (sx < rx)) ) then  ! avoid double determination for even Npart
      do A = 1, 4
        Q(A,sx) = Q(A,s1) + Q(A,lx)
        Q(A,rx) = -Q(A,sx)
      end do
    end if
    lx = 2*lx
    if ( last .eqv. .false. ) then ! recursion
      call intmom_rec(Npart, Jmax, ix, sx, x+1)
    end if
  end do  !ix

end subroutine intmom_rec


! **********************************************************************
recursive subroutine intmom_rec_six(Npart, Jmax, i1, s1, x, use_qp_kinematics)
! **********************************************************************
  use ol_momenta_decl_qp, only: L,collconf,softconf

  implicit none
  integer,    intent(in) :: Npart, Jmax, i1, s1, x
  logical,    intent(in) :: use_qp_kinematics
  real(qp) :: sp
  integer :: A
  integer :: ix ! conventional particle numbers
  integer :: lx ! individual 2^(i-1) particles numbers
  integer :: sx ! sums of 2^(i-1) particle numbers
  integer :: rx ! inverse of sx
  logical :: last

  last = .false.
  if (2*x+2 == Npart .or. 2*x+3 == Npart) then
    last = .true.
  end if
  lx = 1 ! adding ext mom 1 <= ix < i1
  do ix = 1, i1 - 1
    sx = s1 + lx
    rx = Jmax + 1 - sx
    if ( (last .eqv. .false.) .or. (mod(Npart,2) == 1 .or. (sx < rx)) ) then
    ! avoid double determination for even Npart

      L(1:5,sx) = L(1:5,s1) + L(1:5,lx)
      L(1:4,rx) = -L(1:4,sx)
      sp = 2*cont_LC_cntrv(L(1:4,s1),L(1:4,lx))
      L(6,sx) = sp + L(6,s1)

      L(5:6,rx) = L(5:6,sx)

    end if
    lx = 2*lx
    if ( last .eqv. .false. ) then ! recursion
      call intmom_rec_six(Npart, Jmax, ix, sx, x+1, use_qp_kinematics)
    end if
  end do  !ix

end subroutine intmom_rec_six


! **********************************************************************
function squeeze_onshell(pinv, masses)
! **********************************************************************
! If 'pinv' is "close" to an element of 'masses', return the mass (positive or negative).
! Otherwise return pinv.
! **********************************************************************
  use kind_types, only: qp
  use ol_loop_parameters_decl_dp, only: ti_os_thresh, mureg
  implicit none
  complex(qp) :: squeeze_onshell
  complex(qp), intent(inout) :: pinv
  real(qp) :: masses(:), mass
  integer :: k
  squeeze_onshell = pinv
  do k = 1, size(masses)
    mass = masses(k)
    if (k /= 1 .and. mass == 0) cycle
    if (abs(abs(pinv)-mass**2)/mureg**2 < ti_os_thresh) then
      squeeze_onshell = sign(mass*mass, real(pinv))
    end if
  end do
end function squeeze_onshell


! **********************************************************************
function momenta_invariants(moms) result(invs)
! **********************************************************************
! Calculate the list of invariants from the momenta 'moms' (complex standard rep.)
! as used by Collier. Apply 'squeeze_onshell' to each invariant with the masses in the theory.
! **********************************************************************
  use kind_types, only: qp
  use ol_external_decl_qp, only: binom2
  use ol_parameters_decl_dp, only: model
  use ol_parameters_decl_qp, only: &
    & wMW, rMW, wMZ, rMZ, wMH, rMH, &
    & wMC, rMC, wMB, rMB, wMT, rMT, &
    & rME, wME, rMM, wMM, rML, wML, &
    & wMA0, rMA0, wMHH, rMHH, wMHp, rMHp
  implicit none
  complex(qp), intent(in) :: moms(:,:)
  complex(qp) :: invs(binom2(size(moms,2)+1))
  complex(qp) :: moms0(0:3,0:size(moms,2))
  real(qp) :: masses(0:12)
  integer :: n, k, a, b
  n = size(moms,2) + 1
  moms0(:,0) = 0
  moms0(:,1:n-1) = moms
  do k = 1, size(invs)
    invs(k) = cont_L_cmplx(moms0(:,mod(k-1,n)) - moms0(:,mod(k+((k-1)/n),n)))
  end do
  masses = 0
  n = 9
  if (wMW == 0) masses(1) = rMW
  if (wMZ == 0) masses(2) = rMZ
  if (wMH == 0) masses(3) = rMH
  if (wMC == 0) masses(4) = rMC
  if (wMB == 0) masses(5) = rMB
  if (wMT == 0) masses(6) = rMT
  if (wME == 0) masses(7) = rME
  if (wMM == 0) masses(8) = rMM
  if (wML == 0) masses(9) = rML
  if (trim(model) == "2hdm") then
    n = 12
    if (wMA0 == 0) masses(10) = rMA0
    if (wMHH == 0) masses(11) = rMHH
    if (wMHp == 0) masses(12) = rMHp
  end if
  do k = 1, size(invs)
    invs(k) = squeeze_onshell(invs(k), masses(0:n))
  end do
end function momenta_invariants

! **********************************************************************
function collier_invariants(moms) result(invs)
! **********************************************************************
! Calculate the list of invariants from the momenta-indices 'moms'
! as used by Collier.
! This function is meant to be used only to compute invariants for
! scalar integrals, namely one-loop bubbles, triangles and boxes.
! **********************************************************************
  use kind_types, only: qp
  use ol_debug, only: ol_error
  use ol_external_decl_qp, only: binom2
  use ol_momenta_decl_qp, only: L
  implicit none
  integer, intent(in) :: moms(:)
  complex(qp) :: invs(binom2(size(moms)+1))

  if(size(moms)==1) then
    invs(1) = L(5,moms(1)) + L(6,moms(1))
  else if(size(moms)==2) then
    invs(1) = L(5,moms(1))         + L(6,moms(1))
    invs(2) = L(5,moms(2)-moms(1)) + L(6,moms(2)-moms(1))
    invs(3) = L(5,moms(2))         + L(6,moms(2))
  else if(size(moms)==3) then
    invs(1) = L(5,moms(1))         + L(6,moms(1))
    invs(2) = L(5,moms(2)-moms(1)) + L(6,moms(2)-moms(1))
    invs(3) = L(5,moms(3)-moms(2)) + L(6,moms(3)-moms(2))
    invs(4) = L(5,moms(3))         + L(6,moms(3))
    invs(5) = L(5,moms(2))         + L(6,moms(2))
    invs(6) = L(5,moms(3)-moms(1)) + L(6,moms(3)-moms(1))
  else
    call ol_error('Collier invariants computed for a non-MI')
    invs(:) = 0
  end if

end function collier_invariants


end module ol_kinematics_qp

