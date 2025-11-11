
module ol_vamp_1_ppllj_nenexuuxg_1_/**/REALKIND
contains

! **********************************************************************
subroutine vamp_1(M)
! P(0:3,nlegs) = incoming external momenta
! Uses tree structures 'wf', factors 'c', and denominators 'den' from loop_ppllj_nenexuuxg_1.
! Sets colour stripped amplitudes A from the module loop_amplitudes_ppllj_nenexuuxg_1.
! **********************************************************************
  use KIND_TYPES, only: REALKIND, QREALKIND, intkind2
  use ol_parameters_decl_/**/DREALKIND, only: l_switch !, kloopmax
  use ol_parameters_decl_/**/QREALKIND ! masses
  use ol_vert_interface_/**/REALKIND
  use ol_prop_interface_/**/REALKIND
  use ol_last_step_/**/REALKIND
  use ol_tables_storage_ppllj_nenexuuxg_1_/**/DREALKIND
  use ol_tensor_sum_storage_ppllj_nenexuuxg_1_/**/REALKIND
  use ol_loop_handling_/**/REALKIND
  use ofred_reduction_/**/REALKIND, only: Hotf_4pt_reduction, Hotf_4pt_reduction_last
  use ofred_reduction_/**/REALKIND, only: Hotf_5pt_reduction, Hotf_5pt_reduction_last
  use ol_loop_reduction_/**/REALKIND, only: TI_bubble_red, TI_triangle_red

  use ol_loop_storage_ppllj_nenexuuxg_1_/**/REALKIND
#ifndef PRECISION_dp
  use ol_loop_storage_ppllj_nenexuuxg_1_/**/DREALKIND, only: &
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

  Gcoeff(:)%j = (-(c(2)*M(1,:)%j)) * den(1)
  call G0_hol_initialisation(ntryL,Gcoeff,G0H32(1),m0h(1),h0tab(:,1),[16,3,8,4],[0,0,0,0],4,1,wf4(:,3))
  call Hloop_QV_A(ntryL,G0H32(1),ex5(:),G0H16(1),m3h2x16(:,1),heltab2x32(:,:,1))
  call Hloop_Q_A(ntryL,G0H16(1),16,0,G1H16(1),n2h16(1))
  Gcoeff(:)%j = (-(c(2)*M(1,:)%j)) * den(1)
  call G0_hol_initialisation(ntryL,Gcoeff,G0H32(1),m0h(2),h0tab(:,2),[16,3,4,8],[0,0,0,0],4,1,wf4(:,3))
  call Hloop_AV_Q(ntryL,G0H32(1),ex5(:),G0H16(1),m3h2x16(:,2),heltab2x32(:,:,2))
  call Hloop_A_Q(ntryL,G0H16(1),16,0,G1H16(2),n2h16(2))
  Gcoeff(:)%j = (-(c(1)*M(1,:)%j)) * den(1)
  call G0_hol_initialisation(ntryL,Gcoeff,G0H32(1),m0h(3),h0tab(:,3),[16,4,3,8],[0,0,0,0],4,1,wf4(:,3))
  call Hloop_UV_W(ntryL,G0H32(1),0,ex5(:),16,G1H16(3),m3h2x16(:,3),heltab2x32(:,:,3))
  call Hloop_QZ_A(ntryL,G1H16(1),wf4(:,3),G1H4(1),ngZu,m3h4x4(:,1),heltab2x16(:,:,1))
  call Hloop_Q_A(ntryL,G1H4(1),19,0,G2H4(1),n2h4(1))
  call Hloop_AZ_Q(ntryL,G1H16(2),wf4(:,3),G1H4(1),ngZu,m3h4x4(:,2),heltab2x16(:,:,2))
  call Hloop_A_Q(ntryL,G1H4(1),19,0,G2H4(2),n2h4(2))
  call Hloop_VQ_A(ntryL,G1H16(3),ex3(:),G1H8(1),m3h2x8(:,1),heltab2x16(:,:,3))
  call Hloop_Q_A(ntryL,G1H8(1),20,0,G2H8(1),n2h8(1))
  call Hloop_QA_V(ntryL,G2H4(1),ex4(:),G2H2(1),m3h2x2(:,1),heltab2x4(:,:,1))
  call Hloop_AQ_V(ntryL,G2H4(2),ex3(:),G2H2(2),m3h2x2(:,2),heltab2x4(:,:,2))
  call Hotf_4pt_reduction(G2H8(1),RedSet_4(1),mass4set(:,1),  & 
G1H8(1),G1H8(2),G1H8(3),G1H8(4),G1H8(5),8)
  call HG1shiftOLR(G1H8(2),8,8)
  call Hloop_QZ_A(ntryL,G1H8(1),wf4(:,3),G1H2(1),ngZu,m3h4x2(:,1),heltab2x8(:,:,1))
  call Hloop_Q_A(ntryL,G1H2(1),23,0,G2H2(3),n2h2(1))
call HGT_raise_alpha_OLR(G1H8(2),1,1,8)
call HGT_raise_alpha_invQ_OLR(G1H8(2),2,5,8)
  call Hloop_VA_Q(ntryL,G1H8(2),ex4(:),G1H4(1),m3h2x4(:,1),heltab2x8(:,:,2))
  call Hloop_QZ_A(ntryL,G1H8(4),wf4(:,3),G1H2(1),ngZu,m3h4x2(:,2),heltab2x8(:,:,3))
  call Hloop_Q_A(ntryL,G1H2(1),23,0,G2H2(4),n2h2(2))
  call Hloop_QZ_A(ntryL,G1H8(5),wf4(:,3),G1H2(1),ngZu,m3h4x2(:,3),heltab2x8(:,:,4))
  call Hloop_Q_A(ntryL,G1H2(1),23,0,G2H2(5),n2h2(3))
  call Hotf_4pt_reduction(G2H2(1),RedSet_4(2),mass4set(:,1),  & 
G1H2(1),G1H2(2),G1H2(3),G1H2(4),G1H2(5),2)
  call HG1shiftOLR(G1H2(2),4,2)
  call Hotf_4pt_reduction(G2H2(2),RedSet_4(3),mass4set(:,1),  & 
G1H2(6),G1H2(7),G1H2(8),G1H2(9),G1H2(10),2)
  call HG1shiftOLR(G1H2(7),8,2)
  call Hcheck_last_QA_V(ntryL,l_switch,G2H2(3),ex4(:),G2tensor(1),m3h2x1(:,1),heltab2x2(:,:,1))
  call Hcheck_last_QA_V(ntryL,l_switch,G2H2(5),ex4(:),G2tensor(2),m3h2x1(:,2),heltab2x2(:,:,2))
  call Hloop_VQ_A(ntryL,G1H2(1),ex3(:),G1H1(1),m3h2x1(:,3),heltab2x2(:,:,3))
  call Hcheck_last_Q_A(ntryL,l_switch,G1H1(1),31,0,G2tensor(3),n2h1(1))
  call Hloop_VQ_A(ntryL,G1H2(2),ex3(:),G1H1(1),m3h2x1(:,4),heltab2x2(:,:,4))
  call Hcheck_last_Q_A(ntryL,l_switch,G1H1(1),4,0,G2tensor(4),n2h1(2))
  call Hloop_VQ_A(ntryL,G1H2(5),ex3(:),G1H1(1),m3h2x1(:,5),heltab2x2(:,:,5))
  call Hcheck_last_Q_A(ntryL,l_switch,G1H1(1),31,0,G2tensor(5),n2h1(3))
  call Hloop_VA_Q(ntryL,G1H2(6),ex4(:),G1H1(1),m3h2x1(:,6),heltab2x2(:,:,6))
  call Hcheck_last_A_Q(ntryL,l_switch,G1H1(1),31,0,G2tensor(6),n2h1(4))
  call Hloop_VA_Q(ntryL,G1H2(7),ex4(:),G1H1(1),m3h2x1(:,7),heltab2x2(:,:,7))
  call Hcheck_last_A_Q(ntryL,l_switch,G1H1(1),8,0,G2tensor(7),n2h1(5))
  call Hloop_VA_Q(ntryL,G1H2(10),ex4(:),G1H1(1),m3h2x1(:,8),heltab2x2(:,:,8))
  call Hcheck_last_A_Q(ntryL,l_switch,G1H1(1),31,0,G2tensor(8),n2h1(6))
  call Hotf_4pt_reduction_last(G2tensor(1),RedSet_4(1),mass4set(:,1),  & 
G1tensor(1),G1tensor(2),G1tensor(3),G1tensor(4),G1tensor(5))
  call G_TensorShift(G1tensor(2),8)
  call Hotf_4pt_reduction_last(G2tensor(3),RedSet_4(2),mass4set(:,1),  & 
G1tensor(6),G1tensor(7),G1tensor(8),G1tensor(9),G1tensor(10))
  call G_TensorShift(G1tensor(7),4)
  call Hotf_4pt_reduction_last(G2tensor(6),RedSet_4(3),mass4set(:,1),  & 
G1tensor(11),G1tensor(12),G1tensor(13),G1tensor(14),G1tensor(15))
  call G_TensorShift(G1tensor(12),8)
  call Hotf_4pt_reduction_last(G1tensor(1),RedSet_4(1),mass4set(:,1),  & 
G0tensor(1),G0tensor(2),G0tensor(3),G0tensor(4),G0tensor(5))
  call Hotf_4pt_reduction_last(G1tensor(6),RedSet_4(2),mass4set(:,1),  & 
G0tensor(6),G0tensor(7),G0tensor(8),G0tensor(9),G0tensor(10))
  call Hotf_4pt_reduction_last(G1tensor(11),RedSet_4(3),mass4set(:,1),  & 
G0tensor(11),G0tensor(12),G0tensor(13),G0tensor(14),G0tensor(15))
  call ol_merge_tensors(T0sum(1),[G0tensor(1)])
  call ol_merge_tensors(T0sum(2),[G0tensor(6)])
  call ol_merge_tensors(T0sum(3),[G0tensor(11)])
  Gcoeff(:)%j = (-(c(5)*M(1,:)%j)) * den(7)
  call G0_hol_initialisation(ntryL,Gcoeff,G0H32(1),m0h(4),h0tab(:,4),[16,3,12],[0,0,0],3,2,wf4(:,3),wf4(:,7))
  call Hloop_AV_Q(ntryL,G0H32(1),ex5(:),G0H16(1),m3h2x16(:,4),heltab2x32(:,:,4))
  call Hloop_A_Q(ntryL,G0H16(1),16,0,G1H16(1),n2h16(3))
  Gcoeff(:)%j = (-(c(4)*M(1,:)%j)) * den(7)
  call G0_hol_initialisation(ntryL,Gcoeff,G0H32(1),m0h(5),h0tab(:,5),[16,3,12],[nMT,nMT,nMT],3,2,wf4(:,3),wf4(:,7))
  call Hloop_AV_Q(ntryL,G0H32(1),ex5(:),G0H16(1),m3h2x16(:,5),heltab2x32(:,:,5))
  call Hloop_A_Q(ntryL,G0H16(1),16,nMT,G1H16(2),n2h16(4))
  Gcoeff(:)%j = (-(c(5)*M(1,:)%j)) * den(7)
  call G0_hol_initialisation(ntryL,Gcoeff,G0H32(1),m0h(6),h0tab(:,6),[16,3,12],[0,0,0],3,2,wf4(:,3),wf4(:,7))
  call Hloop_QV_A(ntryL,G0H32(1),ex5(:),G0H16(1),m3h2x16(:,6),heltab2x32(:,:,6))
  call Hloop_Q_A(ntryL,G0H16(1),16,0,G1H16(3),n2h16(5))
  Gcoeff(:)%j = (-(c(4)*M(1,:)%j)) * den(7)
  call G0_hol_initialisation(ntryL,Gcoeff,G0H32(1),m0h(7),h0tab(:,7),[16,3,12],[nMT,nMT,nMT],3,2,wf4(:,3),wf4(:,7))
  call Hloop_QV_A(ntryL,G0H32(1),ex5(:),G0H16(1),m3h2x16(:,7),heltab2x32(:,:,7))
  call Hloop_Q_A(ntryL,G0H16(1),16,nMT,G1H16(4),n2h16(6))
  Gcoeff(:)%j = (-(c(5)*M(1,:)%j)) * den(7)
  call G0_hol_initialisation(ntryL,Gcoeff,G0H32(1),m0h(8),h0tab(:,8),[16,3,12],[0,0,0],3,2,wf4(:,3),wf4(:,7))
  call Hloop_AV_Q(ntryL,G0H32(1),ex5(:),G0H16(1),m3h2x16(:,8),heltab2x32(:,:,8))
  call Hloop_A_Q(ntryL,G0H16(1),16,0,G1H16(5),n2h16(7))
  Gcoeff(:)%j = (-(c(4)*M(1,:)%j)) * den(7)
  call G0_hol_initialisation(ntryL,Gcoeff,G0H32(1),m0h(9),h0tab(:,9),[16,3,12],[nMB,nMB,nMB],3,2,wf4(:,3),wf4(:,7))
  call Hloop_AV_Q(ntryL,G0H32(1),ex5(:),G0H16(1),m3h2x16(:,9),heltab2x32(:,:,9))
  call Hloop_A_Q(ntryL,G0H16(1),16,nMB,G1H16(6),n2h16(8))
  Gcoeff(:)%j = (-(c(5)*M(1,:)%j)) * den(7)
  call G0_hol_initialisation(ntryL,Gcoeff,G0H32(1),m0h(10),h0tab(:,10),[16,3,12],[0,0,0],3,2,wf4(:,3),wf4(:,7))
  call Hloop_QV_A(ntryL,G0H32(1),ex5(:),G0H16(1),m3h2x16(:,10),heltab2x32(:,:,10))
  call Hloop_Q_A(ntryL,G0H16(1),16,0,G1H16(7),n2h16(9))
  Gcoeff(:)%j = (-(c(4)*M(1,:)%j)) * den(7)
  call G0_hol_initialisation(ntryL,Gcoeff,G0H32(1),m0h(11),h0tab(:,11),[16,3,12],[nMB,nMB,nMB],3,2,wf4(:,3),wf4(:,7))
  call Hloop_QV_A(ntryL,G0H32(1),ex5(:),G0H16(1),m3h2x16(:,11),heltab2x32(:,:,11))
  call Hloop_Q_A(ntryL,G0H16(1),16,nMB,G1H16(8),n2h16(10))
  Gcoeff(:)%j = (c(3)*M(1,:)%j) * den(3)
  call G0_hol_initialisation(ntryL,Gcoeff,G0H32(1),m0h(12),h0tab(:,12),[20,3,8],[0,0,0],3,2,wf4(:,4),wf4(:,3))
  call Hloop_VQ_A(ntryL,G0H32(1),wf4(:,4),G0H8(1),m3h4x8(:,1),heltab2x32(:,:,12))
  call Hloop_Q_A(ntryL,G0H8(1),20,0,G1H8(1),n2h8(3))
  Gcoeff(:)%j = (c(3)*M(1,:)%j) * den(5)
  call G0_hol_initialisation(ntryL,Gcoeff,G0H32(1),m0h(13),h0tab(:,13),[24,3,4],[0,0,0],3,2,wf4(:,6),wf4(:,3))
  call Hloop_VA_Q(ntryL,G0H32(1),wf4(:,6),G0H8(1),m3h4x8(:,2),heltab2x32(:,:,13))
  call Hloop_A_Q(ntryL,G0H8(1),24,0,G1H8(2),n2h8(4))
  Gcoeff(:)%j = (-(c(2)*M(1,:)%j)) * den(12)
  call G0_hol_initialisation(ntryL,Gcoeff,G0H32(1),m0h(14),h0tab(:,14),[16,7,8],[0,0,0],3,1,wf8(:,7))
  call Hloop_AV_Q(ntryL,G0H32(1),ex5(:),G0H16(1),m3h2x16(:,12),heltab2x32(:,:,14))
  call Hloop_A_Q(ntryL,G0H16(1),16,0,G1H16(9),n2h16(11))
  Gcoeff(:)%j = (-(c(1)*M(1,:)%j)) * den(12)
  call G0_hol_initialisation(ntryL,Gcoeff,G0H32(1),m0h(15),h0tab(:,15),[16,7,8],[0,0,0],3,1,wf8(:,7))
  call Hloop_UV_W(ntryL,G0H32(1),0,ex5(:),16,G1H16(10),m3h2x16(:,13),heltab2x32(:,:,15))
  Gcoeff(:)%j = (-(c(2)*M(1,:)%j)) * den(9)
  call G0_hol_initialisation(ntryL,Gcoeff,G0H32(1),m0h(16),h0tab(:,16),[16,4,11],[0,0,0],3,1,wf8(:,6))
  call Hloop_AV_Q(ntryL,G0H32(1),ex5(:),G0H16(1),m3h2x16(:,14),heltab2x32(:,:,16))
  call Hloop_A_Q(ntryL,G0H16(1),16,0,G1H16(11),n2h16(12))
  Gcoeff(:)%j = (-(c(1)*M(1,:)%j)) * den(9)
  call G0_hol_initialisation(ntryL,Gcoeff,G0H32(1),m0h(17),h0tab(:,17),[16,4,11],[0,0,0],3,1,wf8(:,6))
  call Hloop_UV_W(ntryL,G0H32(1),0,ex5(:),16,G1H16(12),m3h2x16(:,15),heltab2x32(:,:,17))
  call ol_merge(ntryL,merge_step,merge_mism,merge_tables,merge_hels,G1H8(1),[G1H8(3)])
  call Hloop_QZ_A(ntryL,G1H8(1),wf4(:,3),G1H2(1),ngZu,m3h4x2(:,4),heltab2x8(:,:,5))
  call Hloop_Q_A(ntryL,G1H2(1),23,0,G2H2(1),n2h2(4))
  call Hloop_AZ_Q(ntryL,G1H16(1),wf4(:,3),G1H4(2),ngZu,m3h4x4(:,3),heltab2x16(:,:,4))
  call Hloop_A_Q(ntryL,G1H4(2),19,0,G2H4(1),n2h4(3))
  call Hloop_AZ_Q(ntryL,G1H16(2),wf4(:,3),G1H4(2),ngZu,m3h4x4(:,4),heltab2x16(:,:,5))
  call Hloop_A_Q(ntryL,G1H4(2),19,nMT,G2H4(2),n2h4(4))
  call Hloop_QZ_A(ntryL,G1H16(3),wf4(:,3),G1H4(2),ngZu,m3h4x4(:,5),heltab2x16(:,:,6))
  call Hloop_Q_A(ntryL,G1H4(2),19,0,G2H4(3),n2h4(5))
  call Hloop_QZ_A(ntryL,G1H16(4),wf4(:,3),G1H4(2),ngZu,m3h4x4(:,6),heltab2x16(:,:,7))
  call Hloop_Q_A(ntryL,G1H4(2),19,nMT,G2H4(4),n2h4(6))
  call Hloop_AZ_Q(ntryL,G1H16(5),wf4(:,3),G1H4(2),ngZd,m3h4x4(:,7),heltab2x16(:,:,8))
  call Hloop_A_Q(ntryL,G1H4(2),19,0,G2H4(5),n2h4(7))
  call Hloop_AZ_Q(ntryL,G1H16(6),wf4(:,3),G1H4(2),ngZd,m3h4x4(:,8),heltab2x16(:,:,9))
  call Hloop_A_Q(ntryL,G1H4(2),19,nMB,G2H4(6),n2h4(8))
  call Hloop_QZ_A(ntryL,G1H16(7),wf4(:,3),G1H4(2),ngZd,m3h4x4(:,9),heltab2x16(:,:,10))
  call Hloop_Q_A(ntryL,G1H4(2),19,0,G2H4(7),n2h4(9))
  call Hloop_QZ_A(ntryL,G1H16(8),wf4(:,3),G1H4(2),ngZd,m3h4x4(:,10),heltab2x16(:,:,11))
  call Hloop_Q_A(ntryL,G1H4(2),19,nMB,G2H4(8),n2h4(10))
  call Hloop_AZ_Q(ntryL,G1H8(2),wf4(:,3),G1H2(2),ngZu,m3h4x2(:,5),heltab2x8(:,:,6))
  call Hloop_A_Q(ntryL,G1H2(2),27,0,G2H2(2),n2h2(5))
  call Hloop_AQ_V(ntryL,G1H16(9),wf8(:,7),G1H2(5),m3h8x2(:,1),heltab2x16(:,:,12))
  call Hloop_VQ_A(ntryL,G1H16(10),wf8(:,7),G1H2(6),m3h8x2(:,2),heltab2x16(:,:,13))
  call Hloop_Q_A(ntryL,G1H2(6),23,0,G2H2(3),n2h2(6))
  call Hloop_AQ_V(ntryL,G1H16(11),ex3(:),G1H8(4),m3h2x8(:,2),heltab2x16(:,:,14))
  call Hloop_VQ_A(ntryL,G1H16(12),ex3(:),G1H8(5),m3h2x8(:,3),heltab2x16(:,:,15))
  call Hloop_Q_A(ntryL,G1H8(5),20,0,G2H8(1),n2h8(5))
  call ol_merge(ntryL,merge_step,merge_mism,merge_tables,merge_hels,G2H2(3),[G2H2(4)])
  call ol_merge(ntryL,merge_step,merge_mism,merge_tables,merge_hels,G1H2(5),[G1H2(9)])
  call ol_merge(ntryL,merge_step,merge_mism,merge_tables,merge_hels,G2H4(5),[G2H4(1)])
  call ol_merge(ntryL,merge_step,merge_mism,merge_tables,merge_hels,G2H4(7),[G2H4(3)])
call HGT_w2_OLR(G1H4(1),1,1,4)
call HGT_w2_invQ_OLR(G1H4(1),2,5,4)
  call Hloop_QZ_A(ntryL,G1H4(1),wf4(:,3),G1H1(1),ngZu,m3h4x1(:,1),heltab2x4(:,:,3))
  call Hcheck_last_Q_A(ntryL,l_switch,G1H1(1),31,0,G2tensor(1),n2h1(7))
  call Hcheck_last_QA_V(ntryL,l_switch,G2H2(3),ex4(:),G2tensor(3),m3h2x1(:,9),heltab2x2(:,:,9))
  call Hloop_VQ_A(ntryL,G1H2(3),ex3(:),G1H1(1),m3h2x1(:,10),heltab2x2(:,:,10))
  call Hcheck_last_Q_A(ntryL,l_switch,G1H1(1),31,0,G2tensor(6),n2h1(8))
  call Hloop_VQ_A(ntryL,G1H2(4),ex3(:),G1H1(1),m3h2x1(:,11),heltab2x2(:,:,11))
  call Hcheck_last_Q_A(ntryL,l_switch,G1H1(1),31,0,G2tensor(9),n2h1(9))
  call Hloop_VA_Q(ntryL,G1H2(8),ex4(:),G1H1(1),m3h2x1(:,12),heltab2x2(:,:,12))
  call Hcheck_last_A_Q(ntryL,l_switch,G1H1(1),31,0,G2tensor(10),n2h1(10))
  call Hloop_VA_Q(ntryL,G1H2(5),ex4(:),G1H1(1),m3h2x1(:,13),heltab2x2(:,:,13))
  call Hcheck_last_A_Q(ntryL,l_switch,G1H1(1),31,0,G2tensor(11),n2h1(11))
  call Hcheck_last_QA_V(ntryL,l_switch,G2H2(1),ex4(:),G2tensor(12),m3h2x1(:,14),heltab2x2(:,:,14))
  call Hloop_AV_Q(ntryL,G2H4(5),wf4(:,7),G2H1(1),m3h4x1(:,2),heltab2x4(:,:,4))
  call Hcheck_last_A_Q(ntryL,l_switch,G2H1(1),31,0,G3tensor(1),n2h1(12))
  call Hloop_AV_Q(ntryL,G2H4(2),wf4(:,7),G2H1(1),m3h4x1(:,3),heltab2x4(:,:,5))
  call Hcheck_last_A_Q(ntryL,l_switch,G2H1(1),31,nMT,G3tensor(2),n2h1(13))
  call Hloop_QV_A(ntryL,G2H4(7),wf4(:,7),G2H1(1),m3h4x1(:,4),heltab2x4(:,:,6))
  call Hcheck_last_Q_A(ntryL,l_switch,G2H1(1),31,0,G3tensor(3),n2h1(14))
  call Hloop_QV_A(ntryL,G2H4(4),wf4(:,7),G2H1(1),m3h4x1(:,5),heltab2x4(:,:,7))
  call Hcheck_last_Q_A(ntryL,l_switch,G2H1(1),31,nMT,G3tensor(4),n2h1(15))
  call Hloop_AV_Q(ntryL,G2H4(6),wf4(:,7),G2H1(1),m3h4x1(:,6),heltab2x4(:,:,8))
  call Hcheck_last_A_Q(ntryL,l_switch,G2H1(1),31,nMB,G3tensor(5),n2h1(16))
  call Hloop_QV_A(ntryL,G2H4(8),wf4(:,7),G2H1(1),m3h4x1(:,7),heltab2x4(:,:,9))
  call Hcheck_last_Q_A(ntryL,l_switch,G2H1(1),31,nMB,G3tensor(6),n2h1(17))
  call Hcheck_last_AQ_V(ntryL,l_switch,G2H2(2),ex3(:),G2tensor(13),m3h2x1(:,15),heltab2x2(:,:,15))
  call Hloop_VA_Q(ntryL,G1H8(4),wf8(:,6),G1H1(1),m3h8x1(:,1),heltab2x8(:,:,7))
  call Hcheck_last_A_Q(ntryL,l_switch,G1H1(1),31,0,G2tensor(14),n2h1(18))
  call Hcheck_last_QA_V(ntryL,l_switch,G2H8(1),wf8(:,6),G2tensor(15),m3h8x1(:,2),heltab2x8(:,:,8))
  call ol_merge_tensors(G2tensor(15),[G2tensor(14),G2tensor(2),G1tensor(5),G0tensor(5)])
  call ol_merge_tensors(G2tensor(12),[G2tensor(4),G1tensor(7),G1tensor(3),G0tensor(7),G0tensor(3)])
  call ol_merge_tensors(G3tensor(3),[G3tensor(1),G2tensor(8),G2tensor(5),G1tensor(15),G1tensor(10),G0tensor(15),G0tensor(10)])
  call ol_merge_tensors(G2tensor(13),[G2tensor(7),G1tensor(12),G0tensor(12)])
  call ol_merge_tensors(G2tensor(1),[G1tensor(2),G0tensor(2)])
  call ol_merge_tensors(G2tensor(11),[G2tensor(3),G1tensor(14),G1tensor(4),G0tensor(14),G0tensor(4)])
  call ol_merge_tensors(G2tensor(6),[G1tensor(8),G0tensor(8)])
  call ol_merge_tensors(G2tensor(9),[G1tensor(9),G0tensor(9)])
  call ol_merge_tensors(G2tensor(10),[G1tensor(13),G0tensor(13)])
  call ol_merge_tensors(G3tensor(4),[G3tensor(2)])
  call ol_merge_tensors(G3tensor(6),[G3tensor(5)])
  call TI_triangle_red(G2tensor(15),RedBasis(1),mass3set(:,1),G0tensor(1),G0tensor(6),G0tensor(11),G0tensor(5),M2L1R1)
  call TI_triangle_red(G2tensor(12),RedBasis(3),mass3set(:,1),G0tensor(7),G0tensor(3),G0tensor(15),G0tensor(10),M2L1R1)
  call TI_triangle_red(G3tensor(3),RedBasis(4),mass3set(:,1),G0tensor(12),G0tensor(2),G0tensor(14),G0tensor(4),M2L1R1)
  call TI_triangle_red(G2tensor(13),RedBasis(8),mass3set(:,1),G0tensor(8),G0tensor(9),G0tensor(13),G0tensor(16),M2L1R1)
  call TI_triangle_red(G2tensor(1),RedBasis(9),mass3set(:,1),G0tensor(17),G0tensor(18),G0tensor(19),G0tensor(20),M2L1R1)
  call TI_triangle_red(G2tensor(11),RedBasis(2),mass3set(:,1),G0tensor(21),G0tensor(22),G0tensor(23),G0tensor(24),M2L1R1)
  call TI_triangle_red(G2tensor(6),RedBasis(6),mass3set(:,1),G0tensor(25),G0tensor(26),G0tensor(27),G0tensor(28),M2L1R1)
  call TI_triangle_red(G2tensor(9),RedBasis(5),mass3set(:,1),G0tensor(29),G0tensor(30),G0tensor(31),G0tensor(32),M2L1R1)
  call TI_triangle_red(G2tensor(10),RedBasis(7),mass3set(:,1),G0tensor(33),G0tensor(34),G0tensor(35),G0tensor(36),M2L1R1)
  call TI_triangle_red(G3tensor(4),RedBasis(4),mass3set(:,2),G0tensor(37),G0tensor(38),G0tensor(39),G0tensor(40),M2L1R1,[nMT], &
    G0tensor(41))
  call TI_triangle_red(G3tensor(6),RedBasis(4),mass3set(:,3),G0tensor(42),G0tensor(43),G0tensor(44),G0tensor(45),M2L1R1,[nMB], &
    G0tensor(46))
  call ol_merge_tensors(T0sum(4),[G0tensor(1)])
  call ol_merge_tensors(T0sum(5),[G0tensor(7)])
  call ol_merge_tensors(T0sum(6),[G0tensor(12)])
  call ol_merge_tensors(T0sum(7),[G0tensor(8)])
  call ol_merge_tensors(T0sum(8),[G0tensor(17)])
  call ol_merge_tensors(T0sum(9),[G0tensor(21)])
  call ol_merge_tensors(T0sum(10),[G0tensor(25)])
  call ol_merge_tensors(T0sum(11),[G0tensor(29)])
  call ol_merge_tensors(T0sum(12),[G0tensor(33)])
  call ol_merge_tensors(T0sum(13),[G0tensor(37)])
  call ol_merge_tensors(T0sum(14),[G0tensor(42)])
  Gcoeff(:)%j = (c(3)*M(1,:)%j) * den(10)
  call G0_hol_initialisation(ntryL,Gcoeff,G0H32(1),m0h(18),h0tab(:,18),[20,11],[0,0],2,2,wf4(:,4),wf8(:,6))
  call Hloop_AQ_V(ntryL,G0H32(1),wf4(:,4),G0H8(1),m3h4x8(:,3),heltab2x32(:,:,18))
  Gcoeff(:)%j = (c(3)*M(1,:)%j) * den(13)
  call G0_hol_initialisation(ntryL,Gcoeff,G0H32(1),m0h(19),h0tab(:,19),[24,7],[0,0],2,2,wf4(:,6),wf8(:,7))
  call Hloop_VA_Q(ntryL,G0H32(1),wf4(:,6),G0H8(2),m3h4x8(:,4),heltab2x32(:,:,19))
  call Hloop_A_Q(ntryL,G0H8(2),24,0,G1H8(3),n2h8(6))
  call Hloop_VA_Q(ntryL,G0H8(1),wf8(:,6),G0H1(1),m3h8x1(:,3),heltab2x8(:,:,9))
  call Hcheck_last_A_Q(ntryL,l_switch,G0H1(1),31,0,G1tensor(1),n2h1(19))
  call Hcheck_last_AQ_V(ntryL,l_switch,G1H8(3),wf8(:,7),G1tensor(6),m3h8x1(:,4),heltab2x8(:,:,10))
  call ol_merge_tensors(G0tensor(34),[G0tensor(31),G0tensor(27),G0tensor(18),G0tensor(13),G0tensor(6)])
  call ol_merge_tensors(G1tensor(1),[G0tensor(30),G0tensor(10),G0tensor(11)])
  call ol_merge_tensors(G0tensor(32),[G0tensor(24),G0tensor(4),G0tensor(5)])
  call ol_merge_tensors(G0tensor(19),[G0tensor(9),G0tensor(2),G0tensor(3)])
  call ol_merge_tensors(G0tensor(35),[G0tensor(26),G0tensor(23),G0tensor(15)])
  call ol_merge_tensors(G0tensor(36),[G0tensor(28),G0tensor(14)])
  call ol_merge_tensors(G1tensor(6),[G0tensor(22),G0tensor(20),G0tensor(16)])
call TI_bubble_red(G1tensor(1),20,mass2set(:,1),G0tensor(1),M2L1R1)
call TI_bubble_red(G1tensor(6),24,mass2set(:,1),G0tensor(7),M2L1R1)
  call ol_merge_tensors(T0sum(15),[G0tensor(34)])
  call ol_merge_tensors(T0sum(16),[G0tensor(1)])
  call ol_merge_tensors(T0sum(17),[G0tensor(32)])
  call ol_merge_tensors(T0sum(18),[G0tensor(19)])
  call ol_merge_tensors(T0sum(19),[G0tensor(35)])
  call ol_merge_tensors(T0sum(20),[G0tensor(36)])
  call ol_merge_tensors(T0sum(21),[G0tensor(7)])
  call ol_merge_tensors(T0sum(22),[G0tensor(38)])
  call ol_merge_tensors(T0sum(23),[G0tensor(39)])
  call ol_merge_tensors(T0sum(24),[G0tensor(40)])
  call ol_merge_tensors(T0sum(25),[G0tensor(41)])
  call ol_merge_tensors(T0sum(26),[G0tensor(43)])
  call ol_merge_tensors(T0sum(27),[G0tensor(44)])
  call ol_merge_tensors(T0sum(28),[G0tensor(45)])
  call ol_merge_tensors(T0sum(29),[G0tensor(46)])
! end of process

! end do

end subroutine vamp_1

end module ol_vamp_1_ppllj_nenexuuxg_1_/**/REALKIND
