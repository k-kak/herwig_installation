
module ol_vamp_1_ppllj_nenexbbx_1_dp
contains

! **********************************************************************
subroutine vamp_1(M)
! P(0:3,nlegs) = incoming external momenta
! Uses tree structures 'wf', factors 'c', and denominators 'den' from loop_ppllj_nenexbbx_1.
! Sets colour stripped amplitudes A from the module loop_amplitudes_ppllj_nenexbbx_1.
! **********************************************************************
  use kind_types, only: dp, qp, intkind2
  use ol_parameters_decl_dp, only: l_switch !, kloopmax
  use ol_parameters_decl_qp ! masses
  use ol_vert_interface_dp
  use ol_prop_interface_dp
  use ol_last_step_dp
  use ol_tables_storage_ppllj_nenexbbx_1_dp
  use ol_tensor_sum_storage_ppllj_nenexbbx_1_dp
  use ol_loop_handling_dp
  use ofred_reduction_dp, only: Hotf_4pt_reduction, Hotf_4pt_reduction_last
  use ofred_reduction_dp, only: Hotf_5pt_reduction, Hotf_5pt_reduction_last
  use ol_loop_reduction_dp, only: TI_bubble_red, TI_triangle_red

  use ol_loop_storage_ppllj_nenexbbx_1_dp

  use hol_initialisation_dp, only: G0_hol_initialisation
  use ol_h_vert_interface_dp
  use ol_h_prop_interface_dp
  use ol_h_last_step_dp
  use ol_merging_dp, only: ol_merge, ol_merge_tensors, ol_merge_last

  implicit none

  type(Hpolcont) :: Gcoeff(hel_states)
  type(Hpolcont), intent(in) :: M(1,hel_states)
  integer :: kloop



  if (ntryL==1 .OR. p_switch == 2) Gcoeff(:)%hf = Hel(1:hel_states)


! do kloop = 1, kloopmax
  ! =============================================


! Dressing, otf merging and otf reduction calls to build loop structures

  Gcoeff(:)%j = (c(1)*M(1,:)%j) * den(1)
  call G0_hol_initialisation(ntryL,Gcoeff,G0H16(1),m0h(1),h0tab(:,1),[8,3,4],[nMB,nMB,0],3,1,wf4(:,3))
  call Hloop_VA_Q(ntryL,G0H16(1),ex4(:),G0H8(1),m3h2x8(:,1),heltab2x16(:,:,1))
  call Hloop_A_Q(ntryL,G0H8(1),8,nMB,G1H8(1),n2h8(1))
  call Hloop_AZ_Q(ntryL,G1H8(1),wf4(:,3),G1H2(1),ngZd,m3h4x2(:,1),heltab2x8(:,:,1))
  call Hloop_A_Q(ntryL,G1H2(1),11,nMB,G2H2(1),n2h2(1))
  call Hcheck_last_AQ_V(ntryL,l_switch,G2H2(1),ex3(:),G2tensor(1),m3h2x1(:,1),heltab2x2(:,:,1))
  call TI_triangle_red(G2tensor(1),RedBasis(1),mass3set(:,1),G0tensor(1),G0tensor(2),G0tensor(3),G0tensor(4),M2L1R1,[nMB], &
    G0tensor(5))
  call ol_merge_tensors(T0sum(1),[G0tensor(1)])
  call ol_merge_tensors(T0sum(2),[G0tensor(2)])
  call ol_merge_tensors(T0sum(3),[G0tensor(3)])
  call ol_merge_tensors(T0sum(4),[G0tensor(4)])
  call ol_merge_tensors(T0sum(5),[G0tensor(5)])
! end of process

! end do

end subroutine vamp_1

end module ol_vamp_1_ppllj_nenexbbx_1_dp

