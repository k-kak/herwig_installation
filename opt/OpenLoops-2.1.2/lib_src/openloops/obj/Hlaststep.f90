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


module ol_h_last_step_dp
  use ol_h_vert_interface_dp, only: valid_hol
  use ol_parameters_decl_dp, only: hp_switch

  use ol_loop_handling_dp, only: req_qp_cmp,hol_dealloc_hybrid

  implicit none
  contains

!******************************************************************************
subroutine Hcheck_last_AQ_V(ntry, switch, G_A, J_Q, Gtensor, n, t)
!------------------------------------------------------------------------------
! bare AQ -> V gluon-like interaction
!******************************************************************************
  use kind_types, only: dp, qp, intkind1, intkind2
  use ol_data_types_dp, only: wfun, hol, hcl
  use hel_bookkeeping_dp, only: helbookkeeping_ol_last_vert3
  use ol_last_step_dp, only: check_last_AQ_V

  use ol_last_step_qp, only: check_last_AQ_V_qp => check_last_AQ_V

  implicit none
  integer,           intent(in)    :: switch
  integer(intkind1), intent(in)    :: ntry
  integer(intkind2), intent(inout) :: n(3), t(:,:)
  type(wfun),        intent(in)    :: J_Q(:)
  type(hol),         intent(inout) :: G_A
  type(hcl),         intent(inout) :: Gtensor
  complex(dp) :: G_add(size(Gtensor%cmp))

  complex(qp) :: G_add_qp(size(Gtensor%cmp))

  integer :: h

  if (ntry == 1) call helbookkeeping_ol_last_vert3(ntry, J_Q, G_A, n, t)
  if (.not. valid_hol(G_A, Gtensor)) return

  Gtensor%cmp = 0._dp
  G_add   = 0._dp

  do h = 1, n(3)  ! helicity summation
    call check_last_AQ_V(switch, G_A%j(:,:,:,h), J_Q(t(1,h))%j, G_add)
    Gtensor%cmp = Gtensor%cmp + G_add
  end do


  if (req_qp_cmp(G_A)) then
    Gtensor%cmp_qp = 0._qp
    G_add_qp       = 0._qp
    do h = 1, n(3)  ! helicity summation
      call check_last_AQ_V_qp(switch, G_A%j_qp(:,:,:,h), cmplx(J_Q(t(1,h))%j,kind=qp), G_add_qp)
      Gtensor%cmp_qp = Gtensor%cmp_qp + G_add_qp
    end do
    call hol_dealloc_hybrid(G_A)
  end if


end subroutine Hcheck_last_AQ_V


!******************************************************************************
subroutine Hcheck_last_QA_V(ntry, switch, G_Q, J_A, Gtensor, n, t)
!------------------------------------------------------------------------------
! bare QA -> V gluon-like interaction
!******************************************************************************
  use kind_types, only: dp, qp, intkind1, intkind2
  use ol_data_types_dp, only: wfun, hol, hcl
  use hel_bookkeeping_dp, only: helbookkeeping_ol_last_vert3
  use ol_last_step_dp, only: check_last_QA_V

  use ol_last_step_qp, only: check_last_QA_V_qp => check_last_QA_V

  implicit none
  integer,           intent(in)    :: switch
  integer(intkind1), intent(in)    :: ntry
  integer(intkind2), intent(inout) :: n(3), t(:,:)
  type(wfun),        intent(in)    :: J_A(:)
  type(hol),         intent(inout) :: G_Q
  type(hcl),         intent(inout) :: Gtensor
  complex(dp) :: G_add(size(Gtensor%cmp))

  complex(qp) :: G_add_qp(size(Gtensor%cmp))

  integer :: h

  if (ntry == 1) call helbookkeeping_ol_last_vert3(ntry, J_A, G_Q, n, t)
  if (.not. valid_hol(G_Q, Gtensor)) return

  Gtensor%cmp = 0._dp
  G_add   = 0._dp

  do h = 1, n(3)  ! helicity summation
    call check_last_QA_V(switch, G_Q%j(:,:,:,h), J_A(t(1,h))%j, G_add)
    Gtensor%cmp = Gtensor%cmp + G_add
  end do


  if (req_qp_cmp(G_Q)) then
    Gtensor%cmp_qp = 0._qp
    G_add_qp   = 0._qp
    do h = 1, n(3)  ! helicity summation
      call check_last_QA_V_qp(switch, G_Q%j_qp(:,:,:,h), &
        cmplx(J_A(t(1,h))%j,kind=qp), G_add_qp)
      Gtensor%cmp_qp = Gtensor%cmp_qp + G_add_qp
    end do
    call hol_dealloc_hybrid(G_Q)
  end if


end subroutine Hcheck_last_QA_V


! TODO:  <21-11-18, J.-N. Lang> !
! coupling integer
!******************************************************************************
subroutine Hcheck_last_AQ_Z(ntry, switch, G_A, J_Q, Gtensor, ng_RL, n, t)
!------------------------------------------------------------------------------
! bare AQ -> Z Z-like interaction
!******************************************************************************
  use kind_types, only: dp, qp, intkind1, intkind2
  use ol_data_types_dp, only: wfun, hol, hcl
  use hel_bookkeeping_dp, only: helbookkeeping_ol_last_vert3
  use ol_last_step_dp, only: check_last_AQ_Z
  use ol_parameters_init_dp, only: get_coupling

  use ol_last_step_qp, only: check_last_AQ_Z_qp => check_last_AQ_Z
  use ol_parameters_init_qp, only: get_coupling_qp=>get_coupling

  implicit none
  integer,           intent(in)    :: switch
  integer,           intent(in)    :: ng_RL
  integer(intkind1), intent(in)    :: ntry
  integer(intkind2), intent(inout) :: n(3), t(:,:)
  type(wfun),        intent(in)    :: J_Q(:)
  type(hol),         intent(inout) :: G_A
  type(hcl),         intent(inout) :: Gtensor
  complex(dp) :: G_add(size(Gtensor%cmp))
  complex(dp) :: g_RL(2)

  complex(qp) :: G_add_qp(size(Gtensor%cmp))
  complex(qp) :: g_RL_qp(2)

  integer :: h

  if (ntry == 1) call helbookkeeping_ol_last_vert3(ntry, J_Q, G_A, n, t)
  if (.not. valid_hol(G_A, Gtensor)) return

  Gtensor%cmp = 0._dp
  G_add   = 0._dp
  g_RL = get_coupling(ng_RL)

  do h = 1, n(3)  ! helicity summation
    call check_last_AQ_Z(switch, G_A%j(:,:,:,h), J_Q(t(1,h))%j, G_add, g_RL)
    Gtensor%cmp = Gtensor%cmp + G_add
  end do


  if (req_qp_cmp(G_A)) then
    Gtensor%cmp_qp = 0._qp
    G_add_qp   = 0._qp
    g_RL_qp = get_coupling_qp(ng_RL)
    do h = 1, n(3)  ! helicity summation
      call check_last_AQ_Z_qp(switch, G_A%j_qp(:,:,:,h), &
        cmplx(J_Q(t(1,h))%j,kind=qp), G_add_qp, g_RL_qp)
      Gtensor%cmp_qp = Gtensor%cmp_qp + G_add_qp
    end do
    call hol_dealloc_hybrid(G_A)
  end if


end subroutine Hcheck_last_AQ_Z


!******************************************************************************
subroutine Hcheck_last_QA_Z(ntry, switch, G_Q, J_A, Gtensor, ng_RL, n, t)
!------------------------------------------------------------------------------
! bare QA -> Z Z-like interaction
!******************************************************************************
  use kind_types, only: dp, qp, intkind1, intkind2
  use ol_data_types_dp, only: wfun, hol, hcl
  use hel_bookkeeping_dp, only: helbookkeeping_ol_last_vert3
  use ol_last_step_dp, only: check_last_QA_Z
  use ol_parameters_init_dp, only: get_coupling

  use ol_last_step_qp, only: check_last_QA_Z_qp => check_last_QA_Z
  use ol_parameters_init_qp, only: get_coupling_qp=>get_coupling

  implicit none
  integer,           intent(in)    :: switch
  integer,           intent(in)    :: ng_RL
  integer(intkind1), intent(in)    :: ntry
  integer(intkind2), intent(inout) :: n(3), t(:,:)
  type(wfun),        intent(in)    :: J_A(:)
  type(hol),         intent(inout) :: G_Q
  type(hcl),         intent(inout) :: Gtensor
  complex(dp) :: G_add(size(Gtensor%cmp))
  complex(dp) :: g_RL(2)

  complex(qp) :: G_add_qp(size(Gtensor%cmp))
  complex(qp) :: g_RL_qp(2)

  integer :: h

  if (ntry == 1) call helbookkeeping_ol_last_vert3(ntry, J_A, G_Q, n, t)
  if (.not. valid_hol(G_Q, Gtensor)) return

  Gtensor%cmp = 0._dp
  G_add   = 0._dp
  g_RL = get_coupling(ng_RL)

  do h = 1, n(3)  ! helicity summation
    call check_last_QA_Z(switch, G_Q%j(:,:,:,h), J_A(t(1,h))%j, G_add, g_RL)
    Gtensor%cmp = Gtensor%cmp + G_add
  end do


  if (req_qp_cmp(G_Q)) then
    Gtensor%cmp_qp = 0._qp
    G_add_qp   = 0._qp
    g_RL_qp = get_coupling_qp(ng_RL)
    do h = 1, n(3)  ! helicity summation
      call check_last_QA_Z_qp(switch, G_Q%j_qp(:,:,:,h), &
        cmplx(J_A(t(1,h))%j,kind=qp), G_add_qp, g_RL_qp)
      Gtensor%cmp_qp = Gtensor%cmp_qp + G_add_qp
    end do
    call hol_dealloc_hybrid(G_Q)
  end if


end subroutine Hcheck_last_QA_Z


!******************************************************************************
subroutine Hcheck_last_AQ_W(ntry, switch, G_A, J_Q, Gtensor, n, t)
!------------------------------------------------------------------------------
! bare AQ -> W W-like interaction
!******************************************************************************
  use kind_types, only: dp, qp, intkind1, intkind2
  use ol_data_types_dp, only: wfun, hol, hcl
  use hel_bookkeeping_dp, only: helbookkeeping_ol_last_vert3
  use ol_last_step_dp, only: check_last_AQ_W

  use ol_last_step_qp, only: check_last_AQ_W_qp => check_last_AQ_W

  implicit none
  integer,           intent(in)    :: switch
  integer(intkind1), intent(in)    :: ntry
  integer(intkind2), intent(inout) :: n(3), t(:,:)
  type(wfun),        intent(in)    :: J_Q(:)
  type(hol),         intent(inout) :: G_A
  type(hcl),         intent(inout) :: Gtensor
  complex(dp) :: G_add(size(Gtensor%cmp))

  complex(qp) :: G_add_qp(size(Gtensor%cmp))

  integer :: h

  if (ntry == 1) call helbookkeeping_ol_last_vert3(ntry, J_Q, G_A, n, t)
  if (.not. valid_hol(G_A, Gtensor)) return

  Gtensor%cmp = 0._dp
  G_add   = 0._dp

  do h = 1, n(3)  ! helicity summation
    call check_last_AQ_W(switch, G_A%j(:,:,:,h), J_Q(t(1,h))%j, G_add)
    Gtensor%cmp = Gtensor%cmp + G_add
  end do


  if (req_qp_cmp(G_A)) then
    Gtensor%cmp_qp = 0._qp
    G_add_qp   = 0._qp
    do h = 1, n(3)  ! helicity summation
      call check_last_AQ_W_qp(switch, G_A%j_qp(:,:,:,h), &
        cmplx(J_Q(t(1,h))%j,kind=qp), G_add_qp)
      Gtensor%cmp_qp = Gtensor%cmp_qp + G_add_qp
    end do
    call hol_dealloc_hybrid(G_A)
  end if


end subroutine Hcheck_last_AQ_W


!******************************************************************************
subroutine Hcheck_last_QA_W(ntry, switch, G_Q, J_A, Gtensor, n, t)
!------------------------------------------------------------------------------
! bare QA -> W gluon-like interaction
!******************************************************************************
  use kind_types, only: dp, qp, intkind1, intkind2
  use ol_data_types_dp, only: wfun, hol, hcl
  use hel_bookkeeping_dp, only: helbookkeeping_ol_last_vert3
  use ol_last_step_dp, only: check_last_QA_W

  use ol_last_step_qp, only: check_last_QA_W_qp => check_last_QA_W

  implicit none
  integer,           intent(in)    :: switch
  integer(intkind1), intent(in)    :: ntry
  integer(intkind2), intent(inout) :: n(3), t(:,:)
  type(wfun),        intent(in)    :: J_A(:)
  type(hol),         intent(inout) :: G_Q
  type(hcl),         intent(inout) :: Gtensor
  complex(dp) :: G_add(size(Gtensor%cmp))

  complex(qp) :: G_add_qp(size(Gtensor%cmp))

  integer :: h

  if (ntry == 1) call helbookkeeping_ol_last_vert3(ntry, J_A, G_Q, n, t)
  if (.not. valid_hol(G_Q, Gtensor)) return

  Gtensor%cmp = 0._dp
  G_add   = 0._dp

  do h = 1, n(3)  ! helicity summation
    call check_last_QA_W(switch, G_Q%j(:,:,:,h), J_A(t(1,h))%j, G_add)
    Gtensor%cmp = Gtensor%cmp + G_add
  end do


  if (req_qp_cmp(G_Q)) then
    Gtensor%cmp_qp = 0._qp
    G_add_qp   = 0._qp
    do h = 1, n(3)  ! helicity summation
      call check_last_QA_W_qp(switch, G_Q%j_qp(:,:,:,h), &
        cmplx(J_A(t(1,h))%j,kind=qp), G_add_qp)
      Gtensor%cmp_qp = Gtensor%cmp_qp + G_add_qp
    end do
    call hol_dealloc_hybrid(G_Q)
  end if


end subroutine Hcheck_last_QA_W


!******************************************************************************
subroutine Hcheck_last_A_Q(ntry, switch, G_A, mom, M, Gtensor, n)
!------------------------------------------------------------------------------
! dressing anti-quark current with propagator
!******************************************************************************
  use kind_types, only: dp, qp, intkind1, intkind2
  use ol_data_types_dp, only: wfun, hol, hcl
  use hel_bookkeeping_dp, only: helbookkeeping_last_prop
  use ol_last_step_dp, only: check_last_A_Q
  use ol_kinematics_dp, only: get_LC_5,get_mass

  use ol_last_step_qp, only: check_last_A_Q_qp=>check_last_A_Q
  use ol_kinematics_qp, only: get_mass_qp=>get_mass
  use ol_kinematics_dp, only: get_LC_5_qp

  implicit none
  integer,           intent(in)    :: switch, mom
  integer(intkind1), intent(in)    :: ntry
  integer(intkind2), intent(inout) :: n
  integer,           intent(in)    :: M
  type(hol),         intent(inout) :: G_A
  type(hcl),         intent(inout) :: Gtensor
  complex(dp) :: G_add(size(Gtensor%cmp))

  complex(qp) :: G_add_qp(size(Gtensor%cmp))

  integer :: h

  if (ntry == 1) call helbookkeeping_last_prop(ntry, G_A, n)
  if (.not. valid_hol(G_A, Gtensor)) return

  Gtensor%cmp = 0._dp
  G_add   = 0._dp

  do h = 1, n  ! helicity summation
    call check_last_A_Q(switch,G_A%j(:,:,:,h),get_LC_5(mom),get_mass(M),G_add)
    Gtensor%cmp = Gtensor%cmp + G_add
  end do


  if (req_qp_cmp(G_A)) then
    Gtensor%cmp_qp = 0._qp
    G_add_qp   = 0._qp
    do h = 1, n  ! helicity summation
      call check_last_A_Q_qp(switch,G_A%j_qp(:,:,:,h),get_LC_5_qp(mom), &
                             get_mass_qp(M),G_add_qp)
      Gtensor%cmp_qp = Gtensor%cmp_qp + G_add_qp
    end do
    call hol_dealloc_hybrid(G_A)
  end if


end subroutine Hcheck_last_A_Q


!******************************************************************************
subroutine Hcheck_last_Q_A(ntry, switch, G_Q, mom, M, Gtensor, n)
!------------------------------------------------------------------------------
! dressing quark current with propagator
!******************************************************************************
  use kind_types, only: dp, qp, intkind1, intkind2
  use ol_data_types_dp, only: wfun, hol, hcl
  use hel_bookkeeping_dp, only: helbookkeeping_last_prop
  use ol_last_step_dp, only: check_last_Q_A
  use ol_kinematics_dp, only: get_LC_5,get_mass

  use ol_last_step_qp, only: check_last_Q_A_qp => check_last_Q_A
  use ol_kinematics_qp, only: get_mass_qp=>get_mass
  use ol_kinematics_dp, only: get_LC_5_qp

  implicit none
  integer,           intent(in)    :: switch,mom
  integer(intkind1), intent(in)    :: ntry
  integer(intkind2), intent(inout) :: n
  integer,           intent(in)    :: M
  type(hol),         intent(inout) :: G_Q
  type(hcl),         intent(inout) :: Gtensor
  complex(dp) :: G_add(size(Gtensor%cmp))

  complex(qp) :: G_add_qp(size(Gtensor%cmp))

  integer :: h

  if (ntry == 1) call helbookkeeping_last_prop(ntry, G_Q, n)
  if (.not. valid_hol(G_Q, Gtensor)) return

  Gtensor%cmp = 0._dp
  G_add   = 0._dp

  do h = 1, n  ! helicity summation
    call check_last_Q_A(switch, G_Q%j(:,:,:,h),get_LC_5(mom),get_mass(M),G_add)
    Gtensor%cmp = Gtensor%cmp + G_add
  end do


  if (req_qp_cmp(G_Q)) then
    Gtensor%cmp_qp = 0._qp
    G_add_qp   = 0._qp
    do h = 1, n  ! helicity summation
      call check_last_Q_A_qp(switch,G_Q%j_qp(:,:,:,h),get_LC_5_qp(mom), &
                             get_mass_qp(M),G_add_qp)
      Gtensor%cmp_qp = Gtensor%cmp_qp + G_add_qp
    end do
    call hol_dealloc_hybrid(G_Q)
  end if


end subroutine Hcheck_last_Q_A


!******************************************************************************
subroutine Hcheck_last_UV_W(ntry, switch, Gin_V, moml, J_V, momt, Gtensor, n, t)
!------------------------------------------------------------------------------
! bare VV -> V vertex
!******************************************************************************
  use kind_types, only: dp, intkind1, intkind2, qp
  use ol_data_types_dp, only: wfun, hol, hcl
  use hel_bookkeeping_dp, only: helbookkeeping_ol_last_vert3
  use ol_last_step_dp, only: check_last_UV_W
  use ol_kinematics_dp, only: get_LC_4

  use ol_last_step_qp, only: check_last_UV_W_qp => check_last_UV_W
  use ol_kinematics_dp, only: get_LC_4_qp

  implicit none
  integer,           intent(in)    :: switch,moml,momt
  integer(intkind1), intent(in)    :: ntry
  integer(intkind2), intent(inout) :: n(3), t(:,:)
  type(wfun),        intent(in)    :: J_V(:)
  type(hol),         intent(inout) :: Gin_V
  type(hcl),         intent(inout) :: Gtensor
  complex(dp) :: G_add(size(Gtensor%cmp))

  complex(qp) :: G_add_qp(size(Gtensor%cmp))

  integer :: h

  if (ntry == 1) call helbookkeeping_ol_last_vert3(ntry, J_V, Gin_V, n, t)
  if (.not. valid_hol(Gin_V, Gtensor)) return

  Gtensor%cmp = 0._dp
  G_add   = 0._dp

  do h = 1, n(3)  ! helicity summation
    call check_last_UV_W(switch, Gin_V%j(:,:,:,h),get_LC_4(moml),J_V(t(1,h))%j,&
    get_LC_4(momt), G_add)
    Gtensor%cmp = Gtensor%cmp + G_add
  end do


  if (req_qp_cmp(Gin_V)) then
    Gtensor%cmp_qp = 0._qp
    G_add_qp   = 0._qp
    do h = 1, n(3)  ! helicity summation
      call check_last_UV_W_qp(switch, Gin_V%j_qp(:,:,:,h),get_LC_4_qp(moml), &
        cmplx(J_V(t(1,h))%j,kind=qp), &
                    get_LC_4_qp(momt), G_add_qp)
      Gtensor%cmp_qp = Gtensor%cmp_qp + G_add_qp
    end do
    call hol_dealloc_hybrid(Gin_V)
  end if


end subroutine Hcheck_last_UV_W


!******************************************************************************
subroutine Hcheck_last_UW_V(ntry, switch, Gin_V, moml, J_V, momt,Gtensor,n,t)
!------------------------------------------------------------------------------
! bare VV -> V vertex
!******************************************************************************
  use kind_types, only: dp, qp, intkind1, intkind2
  use ol_data_types_dp, only: wfun, hol, hcl
  use hel_bookkeeping_dp, only: helbookkeeping_ol_last_vert3
  use ol_last_step_dp, only: check_last_UW_V
  use ol_kinematics_dp, only: get_LC_4

  use ol_last_step_qp, only: check_last_UW_V_qp=>check_last_UW_V
  use ol_kinematics_dp, only: get_LC_4_qp

  implicit none
  integer,           intent(in)    :: switch
  integer(intkind1), intent(in)    :: ntry
  integer(intkind2), intent(inout) :: n(3), t(:,:)
  integer,           intent(in)    :: moml,momt
  type(wfun),        intent(in)    :: J_V(:)
  type(hol),         intent(inout) :: Gin_V
  type(hcl),         intent(inout) :: Gtensor
  complex(dp) :: G_add(size(Gtensor%cmp))
  integer :: h

  complex(qp) :: G_add_qp(size(Gtensor%cmp))


  if (ntry == 1) call helbookkeeping_ol_last_vert3(ntry, J_V, Gin_V, n, t)
  if (.not. valid_hol(Gin_V, Gtensor)) return

  Gtensor%cmp = 0._dp
  G_add = 0._dp

  do h = 1, n(3)  ! helicity summation
    call check_last_UW_V(switch, Gin_V%j(:,:,:,h),get_LC_4(moml),J_V(t(1,h))%j,&
                         get_LC_4(momt), G_add)
    Gtensor%cmp = Gtensor%cmp + G_add
  end do


  if (req_qp_cmp(Gin_V)) then
    Gtensor%cmp_qp = 0._qp
    G_add_qp   = 0._qp
    do h = 1, n(3)  ! helicity summation
      call check_last_UW_V_qp(switch, Gin_V%j_qp(:,:,:,h),get_LC_4_qp(moml), &
                              cmplx(J_V(t(1,h))%j,kind=qp),get_LC_4_qp(momt), G_add_qp)
      Gtensor%cmp_qp = Gtensor%cmp_qp + G_add_qp
    end do
    call hol_dealloc_hybrid(Gin_V)
  end if


end subroutine Hcheck_last_UW_V

!******************************************************************************
subroutine Hcheck_last_VE_V(ntry, switch, Gin, J1, J2, Gtensor, n, t)
!------------------------------------------------------------------------------
! bare 4-gluon vertex when the sigma particle enters the
! loop as a tree wavefunction
!******************************************************************************
  use kind_types, only: dp, qp, intkind1, intkind2
  use ol_data_types_dp, only: wfun, hol, hcl
  use hel_bookkeeping_dp, only: helbookkeeping_ol_last_vert4
  use ol_last_step_dp, only: check_last_VE_V

  use ol_last_step_qp, only: check_last_VE_V_qp=>check_last_VE_V

  implicit none
  integer,           intent(in)    :: switch
  integer(intkind1), intent(in)    :: ntry
  integer(intkind2), intent(inout) :: n(4), t(:,:)
  type(wfun),        intent(in)    :: J1(:), J2(:)
  type(hol),         intent(inout) :: Gin
  type(hcl),         intent(inout) :: Gtensor
  complex(dp) :: G_add(size(Gtensor%cmp))
  integer :: h

  complex(qp) :: G_add_qp(size(Gtensor%cmp))


  if (ntry == 1) call helbookkeeping_ol_last_vert4(ntry, J1, J2, Gin, n, t)
  if (.not. valid_hol(Gin, Gtensor)) return

  Gtensor%cmp = 0._dp
  G_add = 0._dp

  do h = 1, n(4)  ! helicity summation
    call check_last_VE_V(switch, Gin%j(:,:,:,h), J1(t(1,h))%j, &
                         J2(t(2,h))%j, G_add)
    Gtensor%cmp = Gtensor%cmp + G_add
  end do


  if (req_qp_cmp(Gin)) then
    Gtensor%cmp_qp = 0._qp
    G_add_qp   = 0._qp
    do h = 1, n(4)  ! helicity summation
      call check_last_VE_V_qp(switch, Gin%j_qp(:,:,:,h), cmplx(J1(t(1,h))%j,kind=qp), &
                              cmplx(J2(t(2,h))%j,kind=qp), G_add_qp)
      Gtensor%cmp_qp = Gtensor%cmp_qp + G_add_qp
    end do
    call hol_dealloc_hybrid(Gin)
  end if


end subroutine Hcheck_last_VE_V


!******************************************************************************
subroutine Hcheck_last_EV_V(ntry, switch, Gin, J1, J2, Gtensor, n, t)
!------------------------------------------------------------------------------
! bare 4-gluon vertex when the sigma particle is in the loop
!******************************************************************************
  use kind_types, only: dp, qp, intkind1, intkind2
  use ol_data_types_dp, only: wfun, hol, hcl
  use hel_bookkeeping_dp, only: helbookkeeping_ol_last_vert4
  use ol_last_step_dp, only: check_last_EV_V

  use ol_last_step_qp, only: check_last_EV_V_qp=>check_last_EV_V

  implicit none
  integer,           intent(in)    :: switch
  integer(intkind1), intent(in)    :: ntry
  integer(intkind2), intent(inout) :: n(4), t(:,:)
  type(wfun),        intent(in)    :: J1(:), J2(:)
  type(hol),         intent(inout) :: Gin
  type(hcl),         intent(inout) :: Gtensor
  complex(dp) :: G_add(size(Gtensor%cmp))
  integer :: h

  complex(qp) :: G_add_qp(size(Gtensor%cmp))


  if (ntry == 1) call helbookkeeping_ol_last_vert4(ntry, J1, J2, Gin, n, t)
  if (.not. valid_hol(Gin, Gtensor)) return

  Gtensor%cmp = 0._dp
  G_add   = 0._dp

  do h = 1, n(4)  ! helicity summation
    call check_last_EV_V(switch, Gin%j(:,:,:,h), J1(t(1,h))%j, &
                         J2(t(2,h))%j, G_add)
    Gtensor%cmp = Gtensor%cmp + G_add
  end do


  if (req_qp_cmp(Gin)) then
    Gtensor%cmp_qp = 0._qp
    G_add_qp = 0._qp
    do h = 1, n(4)  ! helicity summation
      call check_last_EV_V_qp(switch, Gin%j_qp(:,:,:,h), cmplx(J1(t(1,h))%j,kind=qp), &
                           cmplx(J2(t(2,h))%j,kind=qp), G_add_qp)
      Gtensor%cmp_qp = Gtensor%cmp_qp + G_add_qp
    end do
    call hol_dealloc_hybrid(Gin)
  end if


end subroutine Hcheck_last_EV_V


!******************************************************************************
subroutine Hcheck_last_GGG_G_23(ntry, switch, Gin, J1, J2, Gtensor, n, t)
!------------------------------------------------------------------------------
! bare 4-gluon vertex when the sigma particle enters the
! loop as a tree wavefunction
!******************************************************************************
  use kind_types, only: dp, intkind1, intkind2, qp
  use ol_data_types_dp, only: wfun, hol, hcl
  use hel_bookkeeping_dp, only: helbookkeeping_ol_last_vert4
  use ol_last_step_dp, only: check_last_GGG_G_23

  use ol_last_step_qp, only: check_last_GGG_G_23_qp => check_last_GGG_G_23

  implicit none
  integer,           intent(in)    :: switch
  integer(intkind1), intent(in)    :: ntry
  integer(intkind2), intent(inout) :: n(4), t(:,:)
  type(wfun),        intent(in)    :: J1(:), J2(:)
  type(hol),         intent(inout) :: Gin
  type(hcl),         intent(inout) :: Gtensor
  complex(dp) :: G_add(size(Gtensor%cmp))

  complex(qp) :: G_add_qp(size(Gtensor%cmp))

  integer :: h

  if (ntry == 1) call helbookkeeping_ol_last_vert4(ntry, J1, J2, Gin, n, t)
  if (.not. valid_hol(Gin, Gtensor)) return

  Gtensor%cmp = 0._dp
  G_add   = 0._dp

  do h = 1, n(4)  ! helicity summation
    call check_last_GGG_G_23(switch, Gin%j(:,:,:,h), J1(t(1,h))%j, &
    J2(t(2,h))%j, G_add)
    Gtensor%cmp = Gtensor%cmp + G_add
  end do


  if (req_qp_cmp(Gin)) then
    Gtensor%cmp_qp = 0._qp
    G_add_qp   = 0._qp
    do h = 1, n(4)  ! helicity summation
      call check_last_GGG_G_23_qp(switch, Gin%j_qp(:,:,:,h), cmplx(J1(t(1,h))%j,kind=qp), &
              cmplx(J2(t(2,h))%j,kind=qp), G_add_qp)
      Gtensor%cmp_qp = Gtensor%cmp_qp + G_add_qp
    end do
    call hol_dealloc_hybrid(Gin)
  end if


end subroutine Hcheck_last_GGG_G_23


!******************************************************************************
subroutine Hcheck_last_GGG_G_12(ntry, switch, Gin, J1, J2, Gtensor, n, t)
!------------------------------------------------------------------------------
! bare 4-gluon vertex when the sigma particle is in the loop
!******************************************************************************
  use kind_types, only: dp, intkind1, intkind2, qp
  use ol_data_types_dp, only: wfun, hol, hcl
  use hel_bookkeeping_dp, only: helbookkeeping_ol_last_vert4
  use ol_last_step_dp, only: check_last_GGG_G_12

  use ol_last_step_qp, only: check_last_GGG_G_12_qp => check_last_GGG_G_12

  implicit none
  integer,           intent(in)    :: switch
  integer(intkind1), intent(in)    :: ntry
  integer(intkind2), intent(inout) :: n(4), t(:,:)
  type(wfun),        intent(in)    :: J1(:), J2(:)
  type(hol),         intent(inout) :: Gin
  type(hcl),         intent(inout) :: Gtensor
  complex(dp) :: G_add(size(Gtensor%cmp))

  complex(qp) :: G_add_qp(size(Gtensor%cmp))

  integer :: h

  if (ntry == 1) call helbookkeeping_ol_last_vert4(ntry, J1, J2, Gin, n, t)
  if (.not. valid_hol(Gin, Gtensor)) return

  Gtensor%cmp = 0._dp
  G_add   = 0._dp

  do h = 1, n(4)  ! helicity summation
    call check_last_GGG_G_12(switch, Gin%j(:,:,:,h), J1(t(1,h))%j, &
    J2(t(2,h))%j, G_add)
    Gtensor%cmp = Gtensor%cmp + G_add
  end do


  if (req_qp_cmp(Gin)) then
    Gtensor%cmp_qp = 0._qp
    G_add_qp   = 0._qp
    do h = 1, n(4)  ! helicity summation
      call check_last_GGG_G_12_qp(switch, Gin%j_qp(:,:,:,h), &
        cmplx(J1(t(1,h))%j,kind=qp), &
            cmplx(J2(t(2,h))%j,kind=qp), G_add_qp)
      Gtensor%cmp_qp = Gtensor%cmp_qp + G_add_qp
    end do
    call hol_dealloc_hybrid(Gin)
  end if


end subroutine Hcheck_last_GGG_G_12


!******************************************************************************
subroutine Hcheck_last_CV_D(ntry, switch, Gin, moml, J_V, momt, Gtensor,n,t)
!------------------------------------------------------------------------------
! bare ghost-gluon -> ghost vertex
!******************************************************************************
  use kind_types, only: dp, intkind1, intkind2, qp
  use ol_data_types_dp, only: wfun, hol, hcl
  use hel_bookkeeping_dp, only: helbookkeeping_ol_last_vert3
  use ol_last_step_dp, only: check_last_CV_D
  use ol_kinematics_dp, only: get_LC_4

  use ol_last_step_qp, only: check_last_CV_D_qp => check_last_CV_D
  use ol_kinematics_dp, only: get_LC_4_qp

  implicit none
  integer,           intent(in)    :: switch,moml,momt
  integer(intkind1), intent(in)    :: ntry
  integer(intkind2), intent(inout) :: n(3), t(:,:)
  type(wfun),        intent(in)    :: J_V(:)
  type(hol),         intent(inout) :: Gin
  type(hcl),         intent(inout) :: Gtensor
  complex(dp) :: G_add(size(Gtensor%cmp))

  complex(qp) :: G_add_qp(size(Gtensor%cmp))

  integer :: h

  if (ntry == 1) call helbookkeeping_ol_last_vert3(ntry, J_V, Gin, n, t)
  if (.not. valid_hol(Gin, Gtensor)) return

  Gtensor%cmp = 0._dp
  G_add   = 0._dp

  do h = 1, n(3)  ! helicity summation
    call check_last_CV_D(switch, Gin%j(:,:,:,h), get_LC_4(moml), J_V(t(1,h))%j, &
    get_LC_4(momt), G_add)
    Gtensor%cmp = Gtensor%cmp + G_add
  end do


  if (req_qp_cmp(Gin)) then
    Gtensor%cmp_qp = 0._qp
    G_add_qp   = 0._qp
    do h = 1, n(3)  ! helicity summation
      call check_last_CV_D_qp(switch, Gin%j_qp(:,:,:,h), get_LC_4_qp(moml), &
        cmplx(J_V(t(1,h))%j,kind=qp), &
              get_LC_4_qp(momt), G_add_qp)
      Gtensor%cmp_qp = Gtensor%cmp_qp + G_add_qp
    end do
    call hol_dealloc_hybrid(Gin)
  end if


end subroutine Hcheck_last_CV_D


!******************************************************************************
subroutine Hcheck_last_DV_C(ntry, switch, Gin, moml, J_V, Gtensor, n, t)
!------------------------------------------------------------------------------
! bare antighost-gluon -> antighost vertex
!******************************************************************************
  use kind_types, only: dp, intkind1, intkind2, qp
  use ol_data_types_dp, only: wfun, hol, hcl
  use hel_bookkeeping_dp, only: helbookkeeping_ol_last_vert3
  use ol_last_step_dp, only: check_last_DV_C
  use ol_kinematics_dp, only: get_LC_4

  use ol_last_step_qp, only: check_last_DV_C_qp => check_last_DV_C
  use ol_kinematics_dp, only: get_LC_4_qp

  implicit none
  integer,           intent(in)    :: switch,moml
  integer(intkind1), intent(in)    :: ntry
  integer(intkind2), intent(inout) :: n(3), t(:,:)
  type(wfun),        intent(in)    :: J_V(:)
  type(hol),         intent(inout) :: Gin
  type(hcl),         intent(inout) :: Gtensor
  complex(dp) :: G_add(size(Gtensor%cmp))

  complex(qp) :: G_add_qp(size(Gtensor%cmp))

  integer :: h

  if (ntry == 1) call helbookkeeping_ol_last_vert3(ntry, J_V, Gin, n, t)
  if (.not. valid_hol(Gin, Gtensor)) return

  Gtensor%cmp = 0._dp
  G_add   = 0._dp

  do h = 1, n(3)  ! helicity summation
    call check_last_DV_C(switch, Gin%j(:,:,:,h), get_LC_4(moml), J_V(t(1,h))%j, G_add)
    Gtensor%cmp = Gtensor%cmp + G_add
  end do


  if (req_qp_cmp(Gin)) then
    Gtensor%cmp_qp = 0._qp
    G_add_qp       = 0._qp
    do h = 1, n(3)  ! helicity summation
      call check_last_DV_C_qp(switch, Gin%j_qp(:,:,:,h), get_LC_4_qp(moml), &
        cmplx(J_V(t(1,h))%j,kind=qp), G_add_qp)
      Gtensor%cmp_qp = Gtensor%cmp_qp + G_add_qp
    end do
    call hol_dealloc_hybrid(Gin)
  end if


end subroutine Hcheck_last_DV_C


!******************************************************************************
subroutine Hcheck_last_QA_S(ntry, switch, G_Q, J_A, Gtensor, ng_RL, n, t)
!------------------------------------------------------------------------------
! bare quark-antiquark -> scalar vertex
!******************************************************************************
  use kind_types, only: dp, qp, intkind1, intkind2
  use ol_data_types_dp, only: wfun, hol, hcl
  use hel_bookkeeping_dp, only: helbookkeeping_ol_last_vert3
  use ol_last_step_dp, only: check_last_QA_S
  use ol_parameters_init_dp, only: get_coupling

  use ol_last_step_qp, only: check_last_QA_S_qp=>check_last_QA_S
  use ol_parameters_init_qp, only: get_coupling_qp=>get_coupling

  implicit none
  integer,           intent(in)    :: switch
  integer(intkind1), intent(in)    :: ntry
  integer(intkind2), intent(inout) :: n(3), t(:,:)
  type(wfun),        intent(in)    :: J_A(:)
  integer,           intent(in)    :: ng_RL
  type(hol),         intent(inout) :: G_Q
  type(hcl),         intent(inout) :: Gtensor
  complex(dp) :: G_add(size(Gtensor%cmp))
  complex(dp) :: g_RL(2)
  integer :: h

  complex(qp) :: G_add_qp(size(Gtensor%cmp))
  complex(qp) :: g_RL_qp(2)


  if (ntry == 1) call helbookkeeping_ol_last_vert3(ntry, J_A, G_Q, n, t)
  if (.not. valid_hol(G_Q, Gtensor)) return

  Gtensor%cmp = 0._dp
  G_add   = 0._dp
  g_RL = get_coupling(ng_RL)

  do h = 1, n(3)  ! helicity summation
    call check_last_QA_S(switch, G_Q%j(:,:,:,h), J_A(t(1,h))%j, G_add, g_RL)
    Gtensor%cmp = Gtensor%cmp + G_add
  end do


  if (req_qp_cmp(G_Q)) then
    Gtensor%cmp_qp = 0._qp
    G_add_qp   = 0._qp
    g_RL_qp = get_coupling_qp(ng_RL)
    do h = 1, n(3)  ! helicity summation
      call check_last_QA_S_qp(switch, G_Q%j_qp(:,:,:,h), &
        cmplx(J_A(t(1,h))%j,kind=qp), G_add_qp, g_RL_qp)
      Gtensor%cmp_qp = Gtensor%cmp_qp + G_add_qp
    end do
    call hol_dealloc_hybrid(G_Q)
  end if


end subroutine Hcheck_last_QA_S


!******************************************************************************
subroutine Hcheck_last_AQ_S(ntry, switch, G_A, J_Q, Gtensor, ng_RL, n, t)
!------------------------------------------------------------------------------
! bare antiquark-quark -> scalar vertex
!******************************************************************************
  use kind_types, only: dp, qp, intkind1, intkind2
  use ol_data_types_dp, only: wfun, hol, hcl
  use hel_bookkeeping_dp, only: helbookkeeping_ol_last_vert3
  use ol_last_step_dp, only: check_last_AQ_S
  use ol_parameters_init_dp, only: get_coupling

  use ol_last_step_qp, only: check_last_AQ_S_qp=>check_last_AQ_S
  use ol_parameters_init_qp, only: get_coupling_qp=>get_coupling

  implicit none
  integer,           intent(in)    :: switch
  integer,           intent(in)    :: ng_RL
  integer(intkind1), intent(in)    :: ntry
  integer(intkind2), intent(inout) :: n(3), t(:,:)
  type(wfun),        intent(in)    :: J_Q(:)
  type(hol),         intent(inout) :: G_A
  type(hcl),         intent(inout) :: Gtensor
  integer :: h
  complex(dp) :: G_add(size(Gtensor%cmp))
  complex(dp) :: g_RL(2)

  complex(qp) :: G_add_qp(size(Gtensor%cmp))
  complex(qp) :: g_RL_qp(2)


  if (ntry == 1) call helbookkeeping_ol_last_vert3(ntry, J_Q, G_A, n, t)
  if (.not. valid_hol(G_A, Gtensor)) return

  Gtensor%cmp = 0._dp
  G_add = 0._dp
  g_RL = get_coupling(ng_RL)

  do h = 1, n(3)  ! helicity summation
    call check_last_AQ_S(switch, G_A%j(:,:,:,h), J_Q(t(1,h))%j, G_add, g_RL)
    Gtensor%cmp = Gtensor%cmp + G_add
  end do


  if (req_qp_cmp(G_A)) then
    Gtensor%cmp_qp = 0._qp
    G_add_qp = 0._qp
    g_RL_qp = get_coupling_qp(ng_RL)
    do h = 1, n(3)  ! helicity summation
      call check_last_AQ_S_qp(switch, G_A%j_qp(:,:,:,h), &
        cmplx(J_Q(t(1,h))%j,kind=qp), G_add_qp, g_RL_qp)
      Gtensor%cmp_qp = Gtensor%cmp_qp + G_add_qp
    end do
    call hol_dealloc_hybrid(G_A)
  end if


end subroutine Hcheck_last_AQ_S


!******************************************************************************
subroutine Hcheck_last_VV_S(ntry, switch, G_V, J_V, Gtensor, n, t)
!------------------------------------------------------------------------------
! bare vector-vector -> scalar vertex
!******************************************************************************
  use kind_types, only: dp, qp, intkind1, intkind2
  use ol_data_types_dp, only: wfun, hol, hcl
  use hel_bookkeeping_dp, only: helbookkeeping_ol_last_vert3
  use ol_last_step_dp, only: check_last_VV_S

  use ol_last_step_qp, only: check_last_VV_S_qp=>check_last_VV_S

  implicit none
  integer,           intent(in)    :: switch
  integer(intkind1), intent(in)    :: ntry
  integer(intkind2), intent(inout) :: n(3), t(:,:)
  type(wfun),        intent(in)    :: J_V(:)
  type(hol),         intent(inout) :: G_V
  type(hcl),         intent(inout) :: Gtensor
  integer :: h
  complex(dp) :: G_add(size(Gtensor%cmp))

  complex(qp) :: G_add_qp(size(Gtensor%cmp))


  if (ntry == 1) call helbookkeeping_ol_last_vert3(ntry, J_V, G_V, n, t)
  if (.not. valid_hol(G_V, Gtensor)) return

  Gtensor%cmp = 0._dp
  G_add = 0._dp

  do h = 1, n(3)  ! helicity summation
    call check_last_VV_S(switch, G_V%j(:,:,:,h), J_V(t(1,h))%j, G_add)
    Gtensor%cmp = Gtensor%cmp + G_add
  end do


  if (req_qp_cmp(G_V)) then
    Gtensor%cmp_qp = 0._qp
    G_add_qp = 0._qp
    do h = 1, n(3)  ! helicity summation
      call check_last_VV_S_qp(switch, G_V%j_qp(:,:,:,h), &
        cmplx(J_V(t(1,h))%j,kind=qp), G_add_qp)
      Gtensor%cmp_qp = Gtensor%cmp_qp + G_add_qp
    end do
    call hol_dealloc_hybrid(G_V)
  end if


end subroutine Hcheck_last_VV_S


!******************************************************************************
subroutine Hcheck_last_VS_V(ntry, switch, G_V, J_S, Gtensor, n, t)
!------------------------------------------------------------------------------
! bare scalar-scalar -> scalar vertex
!******************************************************************************
  use kind_types, only: dp, qp, intkind1, intkind2
  use ol_data_types_dp, only: wfun, hol, hcl
  use hel_bookkeeping_dp, only: helbookkeeping_ol_last_vert3
  use ol_last_step_dp, only: check_last_VS_V

  use ol_last_step_qp, only: check_last_VS_V_qp=>check_last_VS_V

  implicit none
  integer,           intent(in)    :: switch
  integer(intkind1), intent(in)    :: ntry
  integer(intkind2), intent(inout) :: n(3), t(:,:)
  type(wfun),        intent(in)    :: J_S(:)
  type(hol),         intent(inout) :: G_V
  type(hcl),         intent(inout) :: Gtensor
  integer :: h
  complex(dp) :: G_add(size(Gtensor%cmp))

  complex(qp) :: G_add_qp(size(Gtensor%cmp))


  if (ntry == 1) call helbookkeeping_ol_last_vert3(ntry, J_S, G_V, n, t)
  if (.not. valid_hol(G_V, Gtensor)) return

  Gtensor%cmp = 0._dp
  G_add = 0._dp

  do h = 1, n(3)  ! helicity summation
    call check_last_VS_V(switch, G_V%j(:,:,:,h), J_S(t(1,h))%j, G_add)
    Gtensor%cmp = Gtensor%cmp + G_add
  end do


  if (req_qp_cmp(G_V)) then
    Gtensor%cmp_qp = 0._qp
    G_add_qp = 0._qp
    do h = 1, n(3)  ! helicity summation
      call check_last_VS_V_qp(switch, G_V%j_qp(:,:,:,h), &
        cmplx(J_S(t(1,h))%j,kind=qp), G_add_qp)
      Gtensor%cmp_qp = Gtensor%cmp_qp + G_add_qp
    end do
    call hol_dealloc_hybrid(G_V)
  end if


end subroutine Hcheck_last_VS_V


!******************************************************************************
subroutine Hcheck_last_SV_V(ntry, switch, G_S, J_V, Gtensor, n, t)
!------------------------------------------------------------------------------
! bare scalar-vector -> vector vertex
!******************************************************************************
  use kind_types, only: dp, qp, intkind1, intkind2
  use ol_data_types_dp, only: wfun, hol, hcl
  use hel_bookkeeping_dp, only: helbookkeeping_ol_last_vert3
  use ol_last_step_dp, only: check_last_SV_V

  use ol_last_step_qp, only: check_last_SV_V_qp=>check_last_SV_V

  implicit none
  integer,           intent(in)    :: switch
  integer(intkind1), intent(in)    :: ntry
  integer(intkind2), intent(inout) :: n(3), t(:,:)
  type(wfun),        intent(in)    :: J_V(:)
  type(hol),         intent(inout) :: G_S
  type(hcl),         intent(inout) :: Gtensor
  integer :: h
  complex(dp) :: G_add(size(Gtensor%cmp))

  complex(qp) :: G_add_qp(size(Gtensor%cmp))


  if (ntry == 1) call helbookkeeping_ol_last_vert3(ntry, J_V, G_S, n, t)
  if (.not. valid_hol(G_S, Gtensor)) return

  Gtensor%cmp = 0._dp
  G_add = 0._dp

  do h = 1, n(3)  ! helicity summation
    call check_last_SV_V(switch, G_S%j(:,:,:,h), J_V(t(1,h))%j, G_add)
    Gtensor%cmp = Gtensor%cmp + G_add
  end do


  if (req_qp_cmp(G_S)) then
    Gtensor%cmp_qp = 0._qp
    G_add_qp = 0._qp
    do h = 1, n(3)  ! helicity summation
      call check_last_SV_V_qp(switch, G_S%j_qp(:,:,:,h), &
        cmplx(J_V(t(1,h))%j,kind=qp), G_add_qp)
      Gtensor%cmp_qp = Gtensor%cmp_qp + G_add_qp
    end do
    call hol_dealloc_hybrid(G_S)
  end if


end subroutine Hcheck_last_SV_V


!******************************************************************************
subroutine Hcheck_last_SV_T(ntry, switch, G_S, moml, J_V, momt, Gtensor, n, t)
!------------------------------------------------------------------------------
! bare scalar-vector -> scalar vertex
!******************************************************************************
  use kind_types, only: dp, qp, intkind1, intkind2
  use ol_data_types_dp, only: wfun, hol, hcl
  use hel_bookkeeping_dp, only: helbookkeeping_ol_last_vert3
  use ol_last_step_dp, only: check_last_SV_T
  use ol_kinematics_dp, only: get_LC_4

  use ol_last_step_qp, only: check_last_SV_T_qp=>check_last_SV_T
  use ol_kinematics_dp, only: get_LC_4_qp

  implicit none
  integer,           intent(in)    :: switch
  integer(intkind1), intent(in)    :: ntry
  integer(intkind2), intent(inout) :: n(3), t(:,:)
  integer,           intent(in)    :: moml,momt
  type(wfun),        intent(in)    :: J_V(:)
  type(hol),         intent(inout) :: G_S
  type(hcl),         intent(inout) :: Gtensor
  integer :: h
  complex(dp) :: G_add(size(Gtensor%cmp))

  complex(qp) :: G_add_qp(size(Gtensor%cmp))


  if (ntry == 1) call helbookkeeping_ol_last_vert3(ntry, J_V, G_S, n, t)
  if (.not. valid_hol(G_S, Gtensor)) return

  Gtensor%cmp = 0._dp
  G_add = 0._dp

  do h = 1, n(3)  ! helicity summation
    call check_last_SV_T(switch, G_S%j(:,:,:,h), get_LC_4(moml), J_V(t(1,h))%j, &
                         get_LC_4(momt), G_add)
    Gtensor%cmp = Gtensor%cmp + G_add
  end do


  if (req_qp_cmp(G_S)) then
    Gtensor%cmp_qp = 0._qp
    G_add_qp = 0._qp
    do h = 1, n(3)  ! helicity summation
      call check_last_SV_T_qp(switch, G_S%j_qp(:,:,:,h), get_LC_4_qp(moml), &
                              cmplx(J_V(t(1,h))%j,kind=qp), get_LC_4_qp(momt), G_add_qp)
      Gtensor%cmp_qp = Gtensor%cmp_qp + G_add_qp
    end do
    call hol_dealloc_hybrid(G_S)
  end if


end subroutine Hcheck_last_SV_T


!******************************************************************************
subroutine Hcheck_last_TV_S(ntry, switch, G_S, moml, J_V, momt, Gtensor, n, t)
!------------------------------------------------------------------------------
! bare scalar-vector -> scalar vertex
!******************************************************************************
  use kind_types, only: dp, qp, intkind1, intkind2
  use ol_data_types_dp, only: wfun, hol, hcl
  use hel_bookkeeping_dp, only: helbookkeeping_ol_last_vert3
  use ol_last_step_dp, only: check_last_TV_S
  use ol_kinematics_dp, only: get_LC_4

  use ol_last_step_qp, only: check_last_TV_S_qp=>check_last_TV_S
  use ol_kinematics_dp, only: get_LC_4_qp

  implicit none
  integer,           intent(in)    :: switch
  integer(intkind1), intent(in)    :: ntry
  integer(intkind2), intent(inout) :: n(3), t(:,:)
  integer,           intent(in)    :: moml,momt
  type(wfun),        intent(in)    :: J_V(:)
  type(hol),         intent(inout) :: G_S
  type(hcl),         intent(inout) :: Gtensor
  integer :: h
  complex(dp) :: G_add(size(Gtensor%cmp))

  complex(qp) :: G_add_qp(size(Gtensor%cmp))


  if (ntry == 1) call helbookkeeping_ol_last_vert3(ntry, J_V, G_S, n, t)
  if (.not. valid_hol(G_S, Gtensor)) return

  Gtensor%cmp = 0._dp
  G_add = 0._dp

  do h = 1, n(3)  ! helicity summation
    call check_last_TV_S(switch, G_S%j(:,:,:,h), get_LC_4(moml), &
                         J_V(t(1,h))%j, get_LC_4(momt), G_add)
    Gtensor%cmp = Gtensor%cmp + G_add
  end do


  if (req_qp_cmp(G_S)) then
    Gtensor%cmp_qp = 0._qp
    G_add_qp = 0._qp
    do h = 1, n(3)  ! helicity summation
      call check_last_TV_S_qp(switch, G_S%j_qp(:,:,:,h), get_LC_4_qp(moml), &
                              cmplx(J_V(t(1,h))%j,kind=qp), get_LC_4_qp(momt), G_add_qp)
      Gtensor%cmp_qp = Gtensor%cmp_qp + G_add_qp
    end do
    call hol_dealloc_hybrid(G_S)
  end if


end subroutine Hcheck_last_TV_S


!******************************************************************************
subroutine Hcheck_last_VS_T(ntry, switch, G_V, moml, J_S, momt, Gtensor, n, t)
!------------------------------------------------------------------------------
! bare vector-scalar -> scalar vertex
!******************************************************************************
  use kind_types, only: dp, qp, intkind1, intkind2
  use ol_data_types_dp, only: wfun, hol, hcl
  use hel_bookkeeping_dp, only: helbookkeeping_ol_last_vert3
  use ol_last_step_dp, only: check_last_VS_T
  use ol_kinematics_dp, only: get_LC_4

  use ol_last_step_qp, only: check_last_VS_T_qp=>check_last_VS_T
  use ol_kinematics_dp, only: get_LC_4_qp

  implicit none
  integer,           intent(in)    :: switch
  integer(intkind1), intent(in)    :: ntry
  integer(intkind2), intent(inout) :: n(3), t(:,:)
  integer,           intent(in)    :: moml,momt
  type(wfun),        intent(in)    :: J_S(:)
  type(hol),         intent(inout) :: G_V
  type(hcl),         intent(inout) :: Gtensor
  integer :: h
  complex(dp) :: G_add(size(Gtensor%cmp))

  complex(qp) :: G_add_qp(size(Gtensor%cmp))


  if (ntry == 1) call helbookkeeping_ol_last_vert3(ntry, J_S, G_V, n, t)
  if (.not. valid_hol(G_V, Gtensor)) return

  Gtensor%cmp = 0._dp
  G_add = 0._dp

  do h = 1, n(3)  ! helicity summation
    call check_last_VS_T(switch, G_V%j(:,:,:,h), get_LC_4(moml), &
                         J_S(t(1,h))%j, get_LC_4(momt), G_add)
    Gtensor%cmp = Gtensor%cmp + G_add
  end do


  if (req_qp_cmp(G_V)) then
    Gtensor%cmp_qp = 0._qp
    G_add_qp = 0._qp
    do h = 1, n(3)  ! helicity summation
      call check_last_VS_T_qp(switch, G_V%j_qp(:,:,:,h), get_LC_4_qp(moml), &
                              cmplx(J_S(t(1,h))%j,kind=qp), get_LC_4_qp(momt), G_add_qp)
      Gtensor%cmp_qp = Gtensor%cmp_qp + G_add_qp
    end do
    call hol_dealloc_hybrid(G_V)
  end if


end subroutine Hcheck_last_VS_T


!******************************************************************************
subroutine Hcheck_last_VT_S(ntry, switch, G_S, moml, J_S, momt, Gtensor, n, t)
!------------------------------------------------------------------------------
! bare vector-scalar -> scalar vertex
!******************************************************************************
  use kind_types, only: dp, qp, intkind1, intkind2
  use ol_data_types_dp, only: wfun, hol, hcl
  use hel_bookkeeping_dp, only: helbookkeeping_ol_last_vert3
  use ol_last_step_dp, only: check_last_VT_S
  use ol_kinematics_dp, only: get_LC_4

  use ol_last_step_qp, only: check_last_VT_S_qp=>check_last_VT_S
  use ol_kinematics_dp, only: get_LC_4_qp

  implicit none
  integer,           intent(in)    :: switch
  integer(intkind1), intent(in)    :: ntry
  integer(intkind2), intent(inout) :: n(3), t(:,:)
  integer,           intent(in)    :: moml,momt
  type(wfun),        intent(in)    :: J_S(:)
  type(hol),         intent(inout) :: G_S
  type(hcl),         intent(inout) :: Gtensor
  integer :: h
  complex(dp) :: G_add(size(Gtensor%cmp))

  complex(qp) :: G_add_qp(size(Gtensor%cmp))


  if (ntry == 1) call helbookkeeping_ol_last_vert3(ntry, J_S, G_S, n, t)
  if (.not. valid_hol(G_S, Gtensor)) return

  Gtensor%cmp = 0._dp
  G_add   = 0._dp

  do h = 1, n(3)  ! helicity summation
    call check_last_VT_S(switch, G_S%j(:,:,:,h), get_LC_4(moml), &
                         J_S(t(1,h))%j, get_LC_4(momt), G_add)
    Gtensor%cmp = Gtensor%cmp + G_add
  end do


  if (req_qp_cmp(G_S)) then
    Gtensor%cmp_qp = 0._qp
    G_add_qp = 0._qp
    do h = 1, n(3)  ! helicity summation
      call check_last_VT_S_qp(switch, G_S%j_qp(:,:,:,h), get_LC_4_qp(moml), &
                           cmplx(J_S(t(1,h))%j,kind=qp), get_LC_4_qp(momt), G_add_qp)
      Gtensor%cmp_qp = Gtensor%cmp_qp + G_add_qp
    end do
    call hol_dealloc_hybrid(G_S)
  end if


end subroutine Hcheck_last_VT_S


!******************************************************************************
subroutine Hcheck_last_ST_V(ntry, switch, G_S, moml, J_S, momt, Gtensor, n, t)
!------------------------------------------------------------------------------
! bare scalar-scalar -> vector vertex
!******************************************************************************
  use kind_types, only: dp, qp, intkind1, intkind2
  use ol_data_types_dp, only: wfun, hol, hcl
  use hel_bookkeeping_dp, only: helbookkeeping_ol_last_vert3
  use ol_last_step_dp, only: check_last_ST_V
  use ol_kinematics_dp, only: get_LC_4

  use ol_last_step_qp, only: check_last_ST_V_qp=>check_last_ST_V
  use ol_kinematics_dp, only: get_LC_4_qp

  implicit none
  integer,           intent(in)    :: switch
  integer(intkind1), intent(in)    :: ntry
  integer(intkind2), intent(inout) :: n(3), t(:,:)
  integer,           intent(in)    :: moml,momt
  type(wfun),        intent(in)    :: J_S(:)
  type(hol),         intent(inout) :: G_S
  type(hcl),         intent(inout) :: Gtensor
  integer :: h
  complex(dp) :: G_add(size(Gtensor%cmp))

  complex(qp) :: G_add_qp(size(Gtensor%cmp))


  if (ntry == 1) call helbookkeeping_ol_last_vert3(ntry, J_S, G_S, n, t)
  if (.not. valid_hol(G_S, Gtensor)) return

  Gtensor%cmp = 0._dp
  G_add = 0._dp

  do h = 1, n(3)  ! helicity summation
    call check_last_ST_V(switch, G_S%j(:,:,:,h), get_LC_4(moml), &
                         J_S(t(1,h))%j, get_LC_4(momt), G_add)
    Gtensor%cmp = Gtensor%cmp + G_add
  end do


  if (req_qp_cmp(G_S)) then
    Gtensor%cmp_qp = 0._qp
    G_add_qp = 0._qp
    do h = 1, n(3)  ! helicity summation
      call check_last_ST_V_qp(switch, G_S%j_qp(:,:,:,h), get_LC_4_qp(moml), &
                              cmplx(J_S(t(1,h))%j,kind=qp), get_LC_4_qp(momt), G_add_qp)
      Gtensor%cmp_qp = Gtensor%cmp_qp + G_add_qp
    end do
    call hol_dealloc_hybrid(G_S)
  end if


end subroutine Hcheck_last_ST_V


!******************************************************************************
subroutine Hcheck_last_TS_V(ntry, switch, G_S, moml, J_S, momt, Gtensor, n, t)
!------------------------------------------------------------------------------
! bare scalar-scalar -> vector vertex
!******************************************************************************
  use kind_types, only: dp, qp, intkind1, intkind2
  use ol_data_types_dp, only: wfun, hol, hcl
  use hel_bookkeeping_dp, only: helbookkeeping_ol_last_vert3
  use ol_last_step_dp, only: check_last_TS_V
  use ol_kinematics_dp, only: get_LC_4

  use ol_last_step_qp, only: check_last_TS_V_qp=>check_last_TS_V
  use ol_kinematics_dp, only: get_LC_4_qp

  implicit none
  integer,           intent(in)    :: switch
  integer(intkind1), intent(in)    :: ntry
  integer(intkind2), intent(inout) :: n(3), t(:,:)
  integer,           intent(in)    :: moml,momt
  type(wfun),        intent(in)    :: J_S(:)
  type(hol),         intent(inout) :: G_S
  type(hcl),         intent(inout) :: Gtensor
  integer :: h
  complex(dp) :: G_add(size(Gtensor%cmp))

  complex(qp) :: G_add_qp(size(Gtensor%cmp))


  if (ntry == 1) call helbookkeeping_ol_last_vert3(ntry, J_S, G_S, n, t)
  if (.not. valid_hol(G_S, Gtensor)) return

  Gtensor%cmp = 0._dp
  G_add   = 0._dp

  do h = 1, n(3)  ! helicity summation
    call check_last_TS_V(switch, G_S%j(:,:,:,h), get_LC_4(moml), &
                         J_S(t(1,h))%j, get_LC_4(momt), G_add)
    Gtensor%cmp = Gtensor%cmp + G_add
  end do


  if (req_qp_cmp(G_S)) then
    Gtensor%cmp_qp = 0._qp
    G_add_qp = 0._qp
    do h = 1, n(3)  ! helicity summation
      call check_last_TS_V_qp(switch, G_S%j_qp(:,:,:,h), get_LC_4_qp(moml), &
                         cmplx(J_S(t(1,h))%j,kind=qp), get_LC_4_qp(momt), G_add_qp)
      Gtensor%cmp_qp = Gtensor%cmp_qp + G_add_qp
    end do
    call hol_dealloc_hybrid(G_S)
  end if


end subroutine Hcheck_last_TS_V


!******************************************************************************
subroutine Hcheck_last_SS_S(ntry, switch, G_S, J_S, Gtensor, n, t)
!------------------------------------------------------------------------------
! bare scalar-scalar -> scalar vertex
!******************************************************************************
  use kind_types, only: dp, qp, intkind1, intkind2
  use ol_data_types_dp, only: wfun, hol, hcl
  use hel_bookkeeping_dp, only: helbookkeeping_ol_last_vert3
  use ol_last_step_dp, only: check_last_SS_S

  use ol_last_step_qp, only: check_last_SS_S_qp=>check_last_SS_S

  implicit none
  integer,           intent(in)    :: switch
  integer(intkind1), intent(in)    :: ntry
  integer(intkind2), intent(inout) :: n(3), t(:,:)
  type(wfun),        intent(in)    :: J_S(:)
  type(hol),         intent(inout) :: G_S
  type(hcl),         intent(inout) :: Gtensor
  integer :: h
  complex(dp) :: G_add(size(Gtensor%cmp))

  complex(qp) :: G_add_qp(size(Gtensor%cmp))


  if (ntry == 1) call helbookkeeping_ol_last_vert3(ntry, J_S, G_S, n, t)
  if (.not. valid_hol(G_S, Gtensor)) return

  Gtensor%cmp = 0._dp
  G_add = 0._dp

  do h = 1, n(3)  ! helicity summation
    call check_last_SS_S(switch, G_S%j(:,:,:,h), J_S(t(1,h))%j, G_add)
    Gtensor%cmp = Gtensor%cmp + G_add
  end do


  if (req_qp_cmp(G_S)) then
    Gtensor%cmp_qp = 0._qp
    G_add_qp = 0._qp
    do h = 1, n(3)  ! helicity summation
      call check_last_SS_S_qp(switch, G_S%j_qp(:,:,:,h), &
        cmplx(J_S(t(1,h))%j,kind=qp), G_add_qp)
      Gtensor%cmp_qp = Gtensor%cmp_qp + G_add_qp
    end do
    call hol_dealloc_hybrid(G_S)
  end if


end subroutine Hcheck_last_SS_S


!******************************************************************************
subroutine Hcheck_last_SSS_S(ntry, switch, Gin_S, J_S1, J_S2, Gtensor, n, t)
!------------------------------------------------------------------------------
! bare scalar-scalar-scalar -> scalar vertex
!******************************************************************************
  use kind_types, only: dp, qp, intkind1, intkind2
  use ol_data_types_dp, only: wfun, hol, hcl
  use hel_bookkeeping_dp, only: helbookkeeping_ol_last_vert4
  use ol_last_step_dp, only: check_last_SSS_S

  use ol_last_step_qp, only: check_last_SSS_S_qp=>check_last_SSS_S

  implicit none
  integer,           intent(in)    :: switch
  integer(intkind1), intent(in)    :: ntry
  integer(intkind2), intent(inout) :: n(4), t(:,:)
  type(wfun),        intent(in)    :: J_S1(:), J_S2(:)
  type(hol),         intent(inout) :: Gin_S
  type(hcl),         intent(inout) :: Gtensor
  integer :: h
  complex(dp) :: G_add(size(Gtensor%cmp))

  complex(qp) :: G_add_qp(size(Gtensor%cmp))


  if (ntry == 1) call helbookkeeping_ol_last_vert4(ntry, J_S1, J_S2, Gin_S, n, t)
  if (.not. valid_hol(Gin_S, Gtensor)) return

  Gtensor%cmp = 0._dp
  G_add = 0._dp

  do h = 1, n(4)  ! helicity summation
    call check_last_SSS_S(switch, Gin_S%j(:,:,:,h), J_S1(t(1,h))%j, &
                          J_S2(t(2,h))%j, G_add)
    Gtensor%cmp = Gtensor%cmp + G_add
  end do


  if (req_qp_cmp(Gin_S)) then
    Gtensor%cmp_qp = 0._qp
    G_add_qp = 0._qp
    do h = 1, n(4)  ! helicity summation
      call check_last_SSS_S_qp(switch, Gin_S%j_qp(:,:,:,h), &
                               cmplx(J_S1(t(1,h))%j,kind=qp), &
                               cmplx(J_S2(t(2,h))%j,kind=qp), G_add_qp)
      Gtensor%cmp_qp = Gtensor%cmp_qp + G_add_qp
    end do
    call hol_dealloc_hybrid(Gin_S)
  end if


end subroutine Hcheck_last_SSS_S


!******************************************************************************
subroutine Hcheck_last_VVS_S(ntry, switch, Gin_V, J_V, J_S, Gtensor, n, t)
!------------------------------------------------------------------------------
! bare vector-vector-scalar -> scalar vertex
!******************************************************************************
  use kind_types, only: dp, qp, intkind1, intkind2
  use ol_data_types_dp, only: wfun, hol, hcl
  use hel_bookkeeping_dp, only: helbookkeeping_ol_last_vert4
  use ol_last_step_dp, only: check_last_VVS_S

  use ol_last_step_qp, only: check_last_VVS_S_qp=>check_last_VVS_S

  implicit none
  integer,           intent(in)    :: switch
  integer(intkind1), intent(in)    :: ntry
  integer(intkind2), intent(inout) :: n(4), t(:,:)
  type(wfun),        intent(in)    :: J_V(:), J_S(:)
  type(hol),         intent(inout) :: Gin_V
  type(hcl),         intent(inout) :: Gtensor
  integer :: h
  complex(dp) :: G_add(size(Gtensor%cmp))

  complex(qp) :: G_add_qp(size(Gtensor%cmp))


  if (ntry == 1) call helbookkeeping_ol_last_vert4(ntry, J_V, J_S, Gin_V, n, t)
  if (.not. valid_hol(Gin_V, Gtensor)) return

  Gtensor%cmp = 0._dp
  G_add = 0._dp

  do h = 1, n(4)  ! helicity summation
    call check_last_VVS_S(switch, Gin_V%j(:,:,:,h), J_V(t(1,h))%j, &
                          J_S(t(2,h))%j, G_add)
    Gtensor%cmp = Gtensor%cmp + G_add
  end do


  if (req_qp_cmp(Gin_V)) then
    Gtensor%cmp_qp = 0._qp
    G_add_qp = 0._qp
    do h = 1, n(4)  ! helicity summation
      call check_last_VVS_S_qp(switch, Gin_V%j_qp(:,:,:,h), &
                               cmplx(J_V(t(1,h))%j,kind=qp), &
                               cmplx(J_S(t(2,h))%j,kind=qp), G_add_qp)
      Gtensor%cmp_qp = Gtensor%cmp_qp + G_add_qp
    end do
    call hol_dealloc_hybrid(Gin_V)
  end if


end subroutine Hcheck_last_VVS_S


!******************************************************************************
subroutine Hcheck_last_SSV_V(ntry, switch, Gin_S, J_S, J_V, Gtensor, n, t)
!------------------------------------------------------------------------------
! bare scalar-scalar-vector -> vector vertex
!******************************************************************************
  use kind_types, only: dp, qp, intkind1, intkind2
  use ol_data_types_dp, only: wfun, hol, hcl
  use hel_bookkeeping_dp, only: helbookkeeping_ol_last_vert4
  use ol_last_step_dp, only: check_last_SSV_V

  use ol_last_step_qp, only: check_last_SSV_V_qp=>check_last_SSV_V

  implicit none
  integer,           intent(in)    :: switch
  integer(intkind1), intent(in)    :: ntry
  integer(intkind2), intent(inout) :: n(4), t(:,:)
  type(wfun),        intent(in)    :: J_S(:), J_V(:)
  type(hol),         intent(inout) :: Gin_S
  type(hcl),         intent(inout) :: Gtensor
  integer :: h
  complex(dp) :: G_add(size(Gtensor%cmp))

  complex(qp) :: G_add_qp(size(Gtensor%cmp))


  if (ntry == 1) call helbookkeeping_ol_last_vert4(ntry, J_S, J_V, Gin_S, n, t)
  if (.not. valid_hol(Gin_S, Gtensor)) return

  Gtensor%cmp = 0._dp
  G_add = 0._dp

  do h = 1, n(4)  ! helicity summation
    call check_last_SSV_V(switch, Gin_S%j(:,:,:,h), J_S(t(1,h))%j, &
                          J_V(t(2,h))%j, G_add)
    Gtensor%cmp = Gtensor%cmp + G_add
  end do


  if (req_qp_cmp(Gin_S)) then
    Gtensor%cmp_qp = 0._qp
    G_add_qp = 0._qp
    do h = 1, n(4)  ! helicity summation
      call check_last_SSV_V_qp(switch, Gin_S%j_qp(:,:,:,h), &
        cmplx(J_S(t(1,h))%j,kind=qp), &
                               cmplx(J_V(t(2,h))%j,kind=qp), G_add_qp)
      Gtensor%cmp_qp = Gtensor%cmp_qp + G_add_qp
    end do
    call hol_dealloc_hybrid(Gin_S)
  end if


end subroutine Hcheck_last_SSV_V


!******************************************************************************
subroutine Hcheck_last_VSS_V(ntry, switch, Gin_V, J_S1, J_S2, Gtensor, n, t)
!------------------------------------------------------------------------------
! bare vector-scalar-scalar -> vector vertex
!******************************************************************************
  use kind_types, only: dp, qp, intkind1, intkind2
  use ol_data_types_dp, only: wfun, hol, hcl
  use hel_bookkeeping_dp, only: helbookkeeping_ol_last_vert4
  use ol_last_step_dp, only: check_last_VSS_V

  use ol_last_step_qp, only: check_last_VSS_V_qp=>check_last_VSS_V

  implicit none
  integer,           intent(in)    :: switch
  integer(intkind1), intent(in)    :: ntry
  integer(intkind2), intent(inout) :: n(4), t(:,:)
  type(wfun),        intent(in)    :: J_S1(:), J_S2(:)
  type(hol),         intent(inout) :: Gin_V
  type(hcl),         intent(inout) :: Gtensor
  integer :: h
  complex(dp) :: G_add(size(Gtensor%cmp))

  complex(qp) :: G_add_qp(size(Gtensor%cmp))


  if (ntry == 1) call helbookkeeping_ol_last_vert4(ntry, J_S1, J_S2, Gin_V, n, t)
  if (.not. valid_hol(Gin_V, Gtensor)) return

  Gtensor%cmp = 0._dp
  G_add = 0._dp

  do h = 1, n(4)  ! helicity summation
    call check_last_VSS_V(switch, Gin_V%j(:,:,:,h), J_S1(t(1,h))%j, &
                          J_S2(t(2,h))%j, G_add)
    Gtensor%cmp = Gtensor%cmp + G_add
  end do


  if (req_qp_cmp(Gin_V)) then
    Gtensor%cmp_qp = 0._qp
    G_add_qp = 0._qp
    do h = 1, n(4)  ! helicity summation
      call check_last_VSS_V_qp(switch, Gin_V%j_qp(:,:,:,h), &
        cmplx(J_S1(t(1,h))%j,kind=qp), &
                            cmplx(J_S2(t(2,h))%j,kind=qp), G_add_qp)
      Gtensor%cmp_qp = Gtensor%cmp_qp + G_add_qp
    end do
    call hol_dealloc_hybrid(Gin_V)
  end if


end subroutine Hcheck_last_VSS_V


!******************************************************************************
subroutine Hcheck_last_SVV_S(ntry, switch, Gin_S, J_V1, J_V2, Gtensor, n, t)
!------------------------------------------------------------------------------
! bare scalar-vector-vector -> scalar vertex
!******************************************************************************
  use kind_types, only: dp, qp, intkind1, intkind2
  use ol_data_types_dp, only: wfun, hol, hcl
  use hel_bookkeeping_dp, only: helbookkeeping_ol_last_vert4
  use ol_last_step_dp, only: check_last_SVV_S

  use ol_last_step_qp, only: check_last_SVV_S_qp=>check_last_SVV_S

  implicit none
  integer,           intent(in)    :: switch
  integer(intkind1), intent(in)    :: ntry
  integer(intkind2), intent(inout) :: n(4), t(:,:)
  type(wfun),        intent(in)    :: J_V1(:), J_V2(:)
  type(hol),         intent(inout) :: Gin_S
  type(hcl),         intent(inout) :: Gtensor
  integer :: h
  complex(dp) :: G_add(size(Gtensor%cmp))

  complex(qp) :: G_add_qp(size(Gtensor%cmp))


  if (ntry == 1) call helbookkeeping_ol_last_vert4(ntry, J_V1, J_V2, Gin_S, n, t)
  if (.not. valid_hol(Gin_S, Gtensor)) return

  Gtensor%cmp = 0._dp
  G_add = 0._dp

  do h = 1, n(4)  ! helicity summation
    call check_last_SVV_S(switch, Gin_S%j(:,:,:,h), J_V1(t(1,h))%j, &
                          J_V2(t(2,h))%j, G_add)
    Gtensor%cmp = Gtensor%cmp + G_add
  end do


  if (req_qp_cmp(Gin_S)) then
    Gtensor%cmp_qp = 0._qp
    G_add_qp = 0._qp
    do h = 1, n(4)  ! helicity summation
      call check_last_SVV_S_qp(switch, Gin_S%j_qp(:,:,:,h), &
        cmplx(J_V1(t(1,h))%j,kind=qp), &
                            cmplx(J_V2(t(2,h))%j,kind=qp), G_add_qp)
      Gtensor%cmp_qp = Gtensor%cmp_qp + G_add_qp
    end do
    call hol_dealloc_hybrid(Gin_S)
  end if


end subroutine Hcheck_last_SVV_S


!******************************************************************************
subroutine Hcheck_last_WWV_V(ntry, switch, Gin_V, J_V1, J_V2, Gtensor, n, t)
!------------------------------------------------------------------------------
! bare V V V -> V vertex
!******************************************************************************
  use kind_types, only: dp, qp, intkind1, intkind2
  use ol_data_types_dp, only: wfun, hol, hcl
  use hel_bookkeeping_dp, only: helbookkeeping_ol_last_vert4
  use ol_last_step_dp, only: check_last_WWV_V

  use ol_last_step_qp, only: check_last_WWV_V_qp=>check_last_WWV_V

  implicit none
  integer,           intent(in)    :: switch
  integer(intkind1), intent(in)    :: ntry
  integer(intkind2), intent(inout) :: n(4), t(:,:)
  type(wfun),        intent(in)    :: J_V1(:), J_V2(:)
  type(hol),         intent(inout) :: Gin_V
  type(hcl),         intent(inout) :: Gtensor
  integer :: h
  complex(dp) :: G_add(size(Gtensor%cmp))

  complex(qp) :: G_add_qp(size(Gtensor%cmp))


  if (ntry == 1) call helbookkeeping_ol_last_vert4(ntry, J_V1, J_V2, Gin_V, n, t)
  if (.not. valid_hol(Gin_V, Gtensor)) return

  Gtensor%cmp = 0._dp
  G_add = 0._dp

  do h = 1, n(4)  ! helicity summation
    call check_last_WWV_V(switch, Gin_V%j(:,:,:,h), J_V1(t(1,h))%j, &
                          J_V2(t(2,h))%j, G_add)
    Gtensor%cmp = Gtensor%cmp + G_add
  end do


  if (req_qp_cmp(Gin_V)) then
    Gtensor%cmp_qp = 0._qp
    G_add_qp = 0._qp
    do h = 1, n(4)  ! helicity summation
      call check_last_WWV_V_qp(switch, Gin_V%j_qp(:,:,:,h), &
        cmplx(J_V1(t(1,h))%j,kind=qp), &
                               cmplx(J_V2(t(2,h))%j,kind=qp), G_add_qp)
      Gtensor%cmp_qp = Gtensor%cmp_qp + G_add_qp
    end do
    call hol_dealloc_hybrid(Gin_V)
  end if


end subroutine Hcheck_last_WWV_V


!******************************************************************************
subroutine Hcheck_last_VWW_V(ntry, switch, Gin_V, J_V1, J_V2, Gtensor, n, t)
!------------------------------------------------------------------------------
! bare V V V -> V vertex
!******************************************************************************
  use kind_types, only: dp, qp, intkind1, intkind2
  use ol_data_types_dp, only: wfun, hol, hcl
  use hel_bookkeeping_dp, only: helbookkeeping_ol_last_vert4
  use ol_last_step_dp, only: check_last_VWW_V

  use ol_last_step_qp, only: check_last_VWW_V_qp=>check_last_VWW_V

  implicit none
  integer,           intent(in)    :: switch
  integer(intkind1), intent(in)    :: ntry
  integer(intkind2), intent(inout) :: n(4), t(:,:)
  type(wfun),        intent(in)    :: J_V1(:), J_V2(:)
  type(hol),         intent(inout) :: Gin_V
  type(hcl),         intent(inout) :: Gtensor
  integer :: h
  complex(dp) :: G_add(size(Gtensor%cmp))

  complex(qp) :: G_add_qp(size(Gtensor%cmp))


  if (ntry == 1) call helbookkeeping_ol_last_vert4(ntry, J_V1, J_V2, Gin_V, n, t)
  if (.not. valid_hol(Gin_V, Gtensor)) return

  Gtensor%cmp = 0._dp
  G_add = 0._dp

  do h = 1, n(4)  ! helicity summation
    call check_last_VWW_V(switch, Gin_V%j(:,:,:,h), J_V1(t(1,h))%j, &
    J_V2(t(2,h))%j, G_add)
    Gtensor%cmp = Gtensor%cmp + G_add
  end do


  if (req_qp_cmp(Gin_V)) then
    Gtensor%cmp_qp = 0._qp
    G_add_qp = 0._qp
    do h = 1, n(4)  ! helicity summation
      call check_last_VWW_V_qp(switch, Gin_V%j_qp(:,:,:,h), &
        cmplx(J_V1(t(1,h))%j,kind=qp), &
                               cmplx(J_V2(t(2,h))%j,kind=qp), G_add_qp)
      Gtensor%cmp_qp = Gtensor%cmp_qp + G_add_qp
    end do
    call hol_dealloc_hybrid(Gin_V)
  end if


end subroutine Hcheck_last_VWW_V


end module ol_h_last_step_dp

