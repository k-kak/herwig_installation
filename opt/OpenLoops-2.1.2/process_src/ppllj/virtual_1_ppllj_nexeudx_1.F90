
module ol_vamp_1_ppllj_nexeudx_1_/**/REALKIND
contains

! **********************************************************************
subroutine vamp_1(M)
! P(0:3,nlegs) = incoming external momenta
! Uses tree structures 'wf', factors 'c', and denominators 'den' from loop_ppllj_nexeudx_1.
! Sets colour stripped amplitudes A from the module loop_amplitudes_ppllj_nexeudx_1.
! **********************************************************************
  use KIND_TYPES, only: REALKIND, QREALKIND, intkind2
  use ol_parameters_decl_/**/DREALKIND, only: l_switch !, kloopmax
  use ol_parameters_decl_/**/QREALKIND ! masses
  use ol_vert_interface_/**/REALKIND
  use ol_prop_interface_/**/REALKIND
  use ol_last_step_/**/REALKIND
  use ol_tables_storage_ppllj_nexeudx_1_/**/DREALKIND
  use ol_tensor_sum_storage_ppllj_nexeudx_1_/**/REALKIND
  use ol_loop_handling_/**/REALKIND
  use ofred_reduction_/**/REALKIND, only: Hotf_4pt_reduction, Hotf_4pt_reduction_last
  use ofred_reduction_/**/REALKIND, only: Hotf_5pt_reduction, Hotf_5pt_reduction_last
  use ol_loop_reduction_/**/REALKIND, only: TI_bubble_red, TI_triangle_red

  use ol_loop_storage_ppllj_nexeudx_1_/**/REALKIND
#ifndef PRECISION_dp
  use ol_loop_storage_ppllj_nexeudx_1_/**/DREALKIND, only: &
    & p_switch, Hel, merge_step, merge_mism, merge_tables, merge_hels, ntryL, hel_states
#endif
  use hol_initialisation_/**/REALKIND, only: G0_hol_initialisation
  use ol_h_vert_interface_/**/REALKIND
  use ol_h_prop_interface_/**/REALKIND
  use ol_h_last_step_/**/REALKIND
  use ol_merging_/**/REALKIND, only: ol_merge, ol_merge_tensors, ol_merge_last

  implicit none

  type(Hpolcont) :: Gcoeff(hel_states)
  type(Hpolcont), intent(in) :: M(1,hel_states)
  integer :: kloop


#ifndef PRECISION_dp
  if (ntryL==1 .OR. p_switch == 1) Gcoeff(:)%hf = Hel(1:hel_states)
#else
  if (ntryL==1 .OR. p_switch == 2) Gcoeff(:)%hf = Hel(1:hel_states)
#endif

! do kloop = 1, kloopmax
  ! =============================================


! Dressing, otf merging and otf reduction calls to build loop structures

  Gcoeff(:)%j = (-(c(1)*M(1,:)%j)) * den(1)
  call G0_hol_initialisation(ntryL,Gcoeff,G0H16(1),m0h(1),h0tab(:,1),[8,3,4],[0,0,0],3,1,wf4(:,3))
  call Hloop_VA_Q(ntryL,G0H16(1),ex4(:),G0H8(1),m3h2x8(:,1),heltab2x16(:,:,1))
  call Hloop_A_Q(ntryL,G0H8(1),8,0,G1H8(1),n2h8(1))
  call Hloop_AW_Q(ntryL,G1H8(1),wf4(:,3),G1H2(1),m3h4x2(:,1),heltab2x8(:,:,1))
  call Hloop_A_Q(ntryL,G1H2(1),11,0,G2H2(1),n2h2(1))
  call Hcheck_last_AQ_V(ntryL,l_switch,G2H2(1),ex3(:),G2tensor(1),m3h2x1(:,1),heltab2x2(:,:,1))
  call TI_triangle_red(G2tensor(1),RedBasis(1),mass3set(:,1),G0tensor(1),G0tensor(2),G0tensor(3),G0tensor(4),M2L1R1)
  call ol_merge_tensors(T0sum(1),[G0tensor(1)])
  call ol_merge_tensors(T0sum(2),[G0tensor(2)])
  call ol_merge_tensors(T0sum(3),[G0tensor(3)])
  call ol_merge_tensors(T0sum(4),[G0tensor(4)])
! end of process

! end do

end subroutine vamp_1

end module ol_vamp_1_ppllj_nexeudx_1_/**/REALKIND
