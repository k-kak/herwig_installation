
! **********************************************************************
module ol_tables_storage_ppllj_nexecsxg_1_/**/REALKIND
! **********************************************************************
  use KIND_TYPES, only: REALKIND, QREALKIND, intkind2
  use ol_data_types_/**/REALKIND, only: hol
  use ol_data_types_/**/REALKIND, only: basis, redset4, redset5
  implicit none

  ! helicity tables for the 1-loop recursion
integer(intkind2), save :: h0tab(32,11)
integer(intkind2), save :: heltab2x2(2,2,15)
integer(intkind2), save :: heltab2x4(2,4,3)
integer(intkind2), save :: heltab2x8(2,8,10)
integer(intkind2), save :: heltab2x16(2,16,7)
integer(intkind2), save :: heltab2x32(2,32,11)


  ! number of helicity states for openloops recursion steps
integer(intkind2), save :: m0h(11)
integer(intkind2), save :: m3h2x1(3,15)
integer(intkind2), save :: m3h4x1(3,1)
integer(intkind2), save :: m3h8x1(3,4)
integer(intkind2), save :: m3h2x2(3,2)
integer(intkind2), save :: m3h4x2(3,5)
integer(intkind2), save :: m3h8x2(3,2)
integer(intkind2), save :: m3h2x4(3,1)
integer(intkind2), save :: m3h4x4(3,2)
integer(intkind2), save :: m3h2x8(3,3)
integer(intkind2), save :: m3h4x8(3,4)
integer(intkind2), save :: m3h2x16(3,7)

integer(intkind2), save :: n2h1(13)
integer(intkind2), save :: n2h2(6)
integer(intkind2), save :: n2h4(2)
integer(intkind2), save :: n2h8(6)
integer(intkind2), save :: n2h16(4)


contains

!**********************************************************************
subroutine HOL_m3_init()
!**********************************************************************
! initialize m3 arrays for helicity summation
!**********************************************************************
  use KIND_TYPES, only: REALKIND, intkind2

m3h2x1(1,:)=2
m3h2x1(2,:)=1
m3h2x1(3,:)=2
m3h4x1(1,:)=4
m3h4x1(2,:)=1
m3h4x1(3,:)=4
m3h8x1(1,:)=8
m3h8x1(2,:)=1
m3h8x1(3,:)=8
m3h2x2(1,:)=2
m3h2x2(2,:)=2
m3h2x2(3,:)=4
m3h4x2(1,:)=4
m3h4x2(2,:)=2
m3h4x2(3,:)=8
m3h8x2(1,:)=8
m3h8x2(2,:)=2
m3h8x2(3,:)=16
m3h2x4(1,:)=2
m3h2x4(2,:)=4
m3h2x4(3,:)=8
m3h4x4(1,:)=4
m3h4x4(2,:)=4
m3h4x4(3,:)=16
m3h2x8(1,:)=2
m3h2x8(2,:)=8
m3h2x8(3,:)=16
m3h4x8(1,:)=4
m3h4x8(2,:)=8
m3h4x8(3,:)=32
m3h2x16(1,:)=2
m3h2x16(2,:)=16
m3h2x16(3,:)=32

n2h1(:)=1
n2h2(:)=2
n2h4(:)=4
n2h8(:)=8
n2h16(:)=16


end subroutine HOL_m3_init

!**********************************************************************

end module ol_tables_storage_ppllj_nexecsxg_1_/**/REALKIND

! **********************************************************************
module ol_tensor_sum_storage_ppllj_nexecsxg_1_/**/REALKIND
! **********************************************************************
  use KIND_TYPES, only: REALKIND, QREALKIND, intkind2
  use ol_data_types_/**/REALKIND, only: hol, hcl, met
  use ol_data_types_/**/REALKIND, only: basis, redset4, redset5
  implicit none

  type(met), save :: M2L1R1

  ! Declarations of loop wave function tensors

  type(hol), save :: G0H1(1)
  type(hol), save :: G1H1(1)
  type(hol), save :: G1H2(10)
  type(hol), save :: G2H2(5)
  type(hol), save :: G1H4(1)
  type(hol), save :: G2H4(2)
  type(hol), save :: G0H8(2)
  type(hol), save :: G1H8(5)
  type(hol), save :: G2H8(1)
  type(hol), save :: G0H16(1)
  type(hol), save :: G1H16(4)
  type(hol), save :: G0H32(1)
  type(hcl), save, dimension(54) :: G0tensor
  type(hcl), save, dimension(15) :: G1tensor
  type(hcl), save, dimension(15) :: G2tensor



  ! Declarations for on-the-fly tensor reduction
type (basis),      save :: RedBasis(9)
type (redset4),    save :: RedSet_4(3)
integer, save :: mass2set(0:1,2)
integer, save :: mass3set(0:2,5)
integer, save :: mass4set(0:3,2)



  ! Declarations for TI calls

  integer, save :: momenta_1(2)
  integer, save :: momenta_2(2)
  integer, save :: momenta_3(2)
  integer, save :: momenta_4(2)
  integer, save :: momenta_5(2)
  integer, save :: momenta_6(2)
  integer, save :: momenta_7(2)
  integer, save :: momenta_8(2)
  integer, save :: momenta_9(3)
  integer, save :: momenta_10(3)
  integer, save :: momenta_11(3)
  integer, save :: momenta_12(3)
  integer, save :: momenta_13(3)
  integer, save :: momenta_14(3)
  integer, save :: momenta_15(3)
  integer, save :: momenta_16(3)
  integer, save :: momenta_17(3)
  integer, save :: momenta_18(4)
  integer, save :: momenta_19(4)
  integer, save :: momenta_20(4)

  integer, save :: masses2_1(2)
  integer, save :: masses2_2(2)
  integer, save :: masses2_3(2)
  integer, save :: masses2_4(2)
  integer, save :: masses2_5(3)
  integer, save :: masses2_6(3)
  integer, save :: masses2_7(3)
  integer, save :: masses2_8(3)
  integer, save :: masses2_9(3)
  integer, save :: masses2_10(4)
  integer, save :: masses2_11(4)



  type(hcl), save, dimension(27) :: T0sum


  contains

!**********************************************************************
subroutine HOL_memory_allocation_full()
!**********************************************************************
! allocation of memory for the types hol
!**********************************************************************
  use KIND_TYPES, only: REALKIND
  use ol_loop_storage_ppllj_nexecsxg_1_/**/DREALKIND, only: nhel
  use ol_data_types_/**/REALKIND, only: hol
  use hol_initialisation_/**/REALKIND, only: hol_allocation
  implicit none

    call hol_allocation(4,1,4,1,G0H1,1)
  call hol_allocation(4,5,4,1,G1H1,1)
  call hol_allocation(4,5,4,2,G1H2,10)
  call hol_allocation(4,15,4,2,G2H2,5)
  call hol_allocation(4,5,4,4,G1H4,1)
  call hol_allocation(4,15,4,4,G2H4,2)
  call hol_allocation(4,1,4,8,G0H8,2)
  call hol_allocation(4,5,4,8,G1H8,5)
  call hol_allocation(4,15,4,8,G2H8,1)
  call hol_allocation(4,1,4,16,G0H16,1)
  call hol_allocation(4,5,4,16,G1H16,4)
  call hol_allocation(4,1,4,32,G0H32,1)


end subroutine HOL_memory_allocation_full

!**********************************************************************
subroutine HOL_memory_allocation_optimized()
!**********************************************************************
! allocation of memory for the types hol
!**********************************************************************
  use KIND_TYPES, only: REALKIND
  use ol_loop_storage_ppllj_nexecsxg_1_/**/DREALKIND, only: nhel
  use ol_data_types_/**/REALKIND, only: hol
  use hol_initialisation_/**/REALKIND, only: hol_allocation
  implicit none

    call hol_allocation(4,1,4,min(nhel,1),G0H1,1)
  call hol_allocation(4,5,4,min(nhel,1),G1H1,1)
  call hol_allocation(4,5,4,min(nhel,2),G1H2,10)
  call hol_allocation(4,15,4,min(nhel,2),G2H2,5)
  call hol_allocation(4,5,4,min(nhel,4),G1H4,1)
  call hol_allocation(4,15,4,min(nhel,4),G2H4,2)
  call hol_allocation(4,1,4,min(nhel,8),G0H8,2)
  call hol_allocation(4,5,4,min(nhel,8),G1H8,5)
  call hol_allocation(4,15,4,min(nhel,8),G2H8,1)
  call hol_allocation(4,1,4,min(nhel,16),G0H16,1)
  call hol_allocation(4,5,4,min(nhel,16),G1H16,4)
  call hol_allocation(4,1,4,min(nhel,32),G0H32,1)


end subroutine HOL_memory_allocation_optimized

!**********************************************************************
subroutine HOL_memory_deallocation_/**/REALKIND(dmode)
!**********************************************************************
! allocation of memory for the types hol
!**********************************************************************
  use KIND_TYPES, only: REALKIND
  use ol_data_types_/**/REALKIND, only: hol
  use hol_initialisation_/**/REALKIND, only: hol_deallocation
  implicit none
  integer,   intent(in)    :: dmode

    call hol_deallocation(G0H1,1,dmode)
  call hol_deallocation(G1H1,1,dmode)
  call hol_deallocation(G1H2,10,dmode)
  call hol_deallocation(G2H2,5,dmode)
  call hol_deallocation(G1H4,1,dmode)
  call hol_deallocation(G2H4,2,dmode)
  call hol_deallocation(G0H8,2,dmode)
  call hol_deallocation(G1H8,5,dmode)
  call hol_deallocation(G2H8,1,dmode)
  call hol_deallocation(G0H16,1,dmode)
  call hol_deallocation(G1H16,4,dmode)
  call hol_deallocation(G0H32,1,dmode)


end subroutine HOL_memory_deallocation_/**/REALKIND

!**********************************************************************
subroutine HCL_memory_allocation()
!**********************************************************************
! allocation of memory for the types hcl
!**********************************************************************
  use KIND_TYPES, only: REALKIND
  use ol_data_types_/**/REALKIND, only: hcl
  use hol_initialisation_/**/REALKIND, only: hcl_allocation
  implicit none

  call hcl_allocation(1,G0tensor, 54)
call hcl_allocation(5,G1tensor, 15)
call hcl_allocation(15,G2tensor, 15)


end subroutine HCL_memory_allocation


!**********************************************************************
subroutine HCL_memory_deallocation_/**/REALKIND(dmode)
!**********************************************************************
! deallocation of memory for the types hcl
!**********************************************************************
  use KIND_TYPES, only: REALKIND
  use ol_data_types_/**/REALKIND, only: hcl
  use hol_initialisation_/**/REALKIND, only: hcl_deallocation
  implicit none
  integer,   intent(in)    :: dmode

  call hcl_deallocation(G0tensor, 54,dmode)
call hcl_deallocation(G1tensor, 15,dmode)
call hcl_deallocation(G2tensor, 15,dmode)

    call hcl_deallocation(T0sum,27,dmode)


end subroutine HCL_memory_deallocation_/**/REALKIND


!**********************************************************************
subroutine Tsum_memory_allocation()
!**********************************************************************
! allocation of memory for the types hcl
!**********************************************************************
  use KIND_TYPES, only: REALKIND
  use ol_data_types_/**/REALKIND, only: hcl
  use hol_initialisation_/**/REALKIND, only: hcl_allocation
  implicit none

    call hcl_allocation(1,T0sum,27)


end subroutine Tsum_memory_allocation


#ifdef PRECISION_dp
subroutine max_point(r) &
    & bind(c,name="ol_f_max_point_ppllj_nexecsxg_1")
  ! Return the number maximal tensor rank
  implicit none
  integer, intent(out) :: r
  r = 4
end subroutine max_point

subroutine tensor_rank(r) &
    & bind(c,name="ol_f_tensor_rank_ppllj_nexecsxg_1")
  ! Return the number maximal tensor rank
  implicit none
  integer, intent(out) :: r
  r = 0
end subroutine tensor_rank
#endif

subroutine reset_tensor_sum()
  use hol_initialisation_/**/REALKIND, only: hcl_allocation
  use ol_parameters_init_/**/REALKIND, only: init_hcl
  implicit none
  integer :: i

  do i = 1,27
    call init_hcl(T0sum(i))
  end do

end subroutine reset_tensor_sum


subroutine scale_one_tsum(tsum, spow)
  use ol_parameters_decl_/**/REALKIND, only: scalefactor
  implicit none
  complex(REALKIND), intent(inout) :: tsum(:)
  integer, intent(in) :: spow ! rank 0 scale power
  real(REALKIND) :: sfinv, sfac
  integer :: sz
  sfinv = 1/scalefactor
  sfac = scalefactor**spow
  sz = size(tsum)
  tsum(1) = sfac*tsum(1)
  if (sz > 1) then ! rank 1
    sfac = sfac*sfinv
    tsum(2:5) = sfac*tsum(2:5)
  end if
  if (sz > 5) then ! rank 2
    sfac = sfac*sfinv
    tsum(6:15) = sfac*tsum(6:15)
  end if
  if (sz > 15) then ! rank 3
    sfac = sfac*sfinv
    tsum(16:35) = sfac*tsum(16:35)
  end if
  if (sz > 35) then ! rank 4
    sfac = sfac*sfinv
    tsum(36:70) = sfac*tsum(36:70)
  end if
  if (sz > 70) then ! rank 5
    sfac = sfac*sfinv
    tsum(71:126) = sfac*tsum(71:126)
  end if
  if (sz > 126) then ! rank 6
    sfac = sfac*sfinv
    tsum(127:210) = sfac*tsum(127:210)
  end if
  if (sz > 210) then ! rank 7
    sfac = sfac*sfinv
    tsum(211:330) = sfac*tsum(211:330)
  end if
end subroutine scale_one_tsum


subroutine scale_tensor_sum()
  implicit none
  call scale_one_tsum(T0sum(1)%cmp, 2)
  call scale_one_tsum(T0sum(2)%cmp, 2)
  call scale_one_tsum(T0sum(3)%cmp, 2)
  call scale_one_tsum(T0sum(4)%cmp, 0)
  call scale_one_tsum(T0sum(5)%cmp, 0)
  call scale_one_tsum(T0sum(6)%cmp, 0)
  call scale_one_tsum(T0sum(7)%cmp, 0)
  call scale_one_tsum(T0sum(8)%cmp, 0)
  call scale_one_tsum(T0sum(9)%cmp, 0)
  call scale_one_tsum(T0sum(10)%cmp, 0)
  call scale_one_tsum(T0sum(11)%cmp, 0)
  call scale_one_tsum(T0sum(12)%cmp, 0)
  call scale_one_tsum(T0sum(13)%cmp, 0)
  call scale_one_tsum(T0sum(14)%cmp, 0)
  call scale_one_tsum(T0sum(15)%cmp, -2)
  call scale_one_tsum(T0sum(16)%cmp, -2)
  call scale_one_tsum(T0sum(17)%cmp, -2)
  call scale_one_tsum(T0sum(18)%cmp, -2)
  call scale_one_tsum(T0sum(19)%cmp, -2)
  call scale_one_tsum(T0sum(20)%cmp, -2)
  call scale_one_tsum(T0sum(21)%cmp, -2)
  call scale_one_tsum(T0sum(22)%cmp, -2)
  call scale_one_tsum(T0sum(23)%cmp, -2)
  call scale_one_tsum(T0sum(24)%cmp, -2)
  call scale_one_tsum(T0sum(25)%cmp, -2)
  call scale_one_tsum(T0sum(26)%cmp, -2)
  call scale_one_tsum(T0sum(27)%cmp, -2)

end subroutine scale_tensor_sum

! **********************************************************************
subroutine set_integral_masses_and_momenta()
! **********************************************************************

  use ol_parameters_decl_/**/REALKIND
  momenta_1 = [ 0, 0 ]
  momenta_2 = [ 16, 15 ]
  momenta_3 = [ 19, 12 ]
  momenta_4 = [ 20, 11 ]
  momenta_5 = [ 23, 8 ]
  momenta_6 = [ 24, 7 ]
  momenta_7 = [ 27, 4 ]
  momenta_8 = [ 28, 3 ]
  momenta_9 = [ 16, 3, 12 ]
  momenta_10 = [ 16, 4, 11 ]
  momenta_11 = [ 16, 7, 8 ]
  momenta_12 = [ 16, 11, 4 ]
  momenta_13 = [ 19, 4, 8 ]
  momenta_14 = [ 19, 8, 4 ]
  momenta_15 = [ 20, 3, 8 ]
  momenta_16 = [ 24, 3, 4 ]
  momenta_17 = [ 24, 4, 3 ]
  momenta_18 = [ 16, 3, 4, 8 ]
  momenta_19 = [ 16, 3, 8, 4 ]
  momenta_20 = [ 16, 4, 3, 8 ]

  masses2_1 = [ 0, 0 ]
  masses2_2 = [ nMC, 0 ]
  masses2_3 = [ 0, nMC ]
  masses2_4 = [ nMC, nMC ]
  masses2_5 = [ 0, 0, 0 ]
  masses2_6 = [ nMC, 0, 0 ]
  masses2_7 = [ 0, 0, nMC ]
  masses2_8 = [ 0, nMC, 0 ]
  masses2_9 = [ nMC, nMC, 0 ]
  masses2_10 = [ 0, 0, nMC, 0 ]
  masses2_11 = [ nMC, nMC, 0, 0 ]


end subroutine  set_integral_masses_and_momenta

! **********************************************************************
subroutine integrate_tensor_sum(M2out)
! **********************************************************************
  use ol_parameters_decl_/**/REALKIND ! only: ZERO, masses
#ifndef PRECISION_dp
  use ol_parameters_decl_/**/DREALKIND, only: a_switch
#endif
  use ol_parameters_init_/**/REALKIND, only: init_met, add_met, met_to_real
  use ol_loop_routines_/**/REALKIND, only: TI_call_OL
  implicit none
  real(REALKIND), intent(out) :: M2out
  type(met) :: M2
  call init_met(M2)



#ifdef LOOPSQUARED
  if (a_switch == 1 .or. a_switch == 7) then
#endif
  call TI_call_OL(0,0, momenta_20, masses2_10, T0sum(1), M2)
  call TI_call_OL(0,0, momenta_19, masses2_11, T0sum(2), M2)
  call TI_call_OL(0,0, momenta_18, masses2_10, T0sum(3), M2)
  call TI_call_OL(0,0, momenta_10, masses2_7, T0sum(4), M2)
  call TI_call_OL(0,0, momenta_15, masses2_8, T0sum(5), M2)
  call TI_call_OL(0,0, momenta_9, masses2_9, T0sum(6), M2)
  call TI_call_OL(0,0, momenta_16, masses2_7, T0sum(7), M2)
  call TI_call_OL(0,0, momenta_9, masses2_7, T0sum(8), M2)
  call TI_call_OL(0,0, momenta_17, masses2_7, T0sum(9), M2)
  call TI_call_OL(0,0, momenta_11, masses2_5, T0sum(10), M2)
  call TI_call_OL(0,0, momenta_14, masses2_6, T0sum(11), M2)
  call TI_call_OL(0,0, momenta_12, masses2_9, T0sum(12), M2)
  call TI_call_OL(0,0, momenta_13, masses2_8, T0sum(13), M2)
  call TI_call_OL(0,0, momenta_10, masses2_9, T0sum(14), M2)
  call TI_call_OL(0,0, momenta_7, masses2_2, T0sum(15), M2)
  call TI_call_OL(0,0, momenta_4, masses2_3, T0sum(16), M2)
  call TI_call_OL(0,0, momenta_2, masses2_1, T0sum(17), M2)
  call TI_call_OL(0,0, momenta_1, masses2_2, T0sum(18), M2)
  call TI_call_OL(0,0, momenta_8, masses2_3, T0sum(19), M2)
  call TI_call_OL(0,0, momenta_5, masses2_1, T0sum(20), M2)
  call TI_call_OL(0,0, momenta_3, masses2_2, T0sum(21), M2)
  call TI_call_OL(0,0, momenta_2, masses2_4, T0sum(22), M2)
  call TI_call_OL(0,0, momenta_8, masses2_2, T0sum(23), M2)
  call TI_call_OL(0,0, momenta_7, masses2_3, T0sum(24), M2)
  call TI_call_OL(0,0, momenta_6, masses2_1, T0sum(25), M2)
  call TI_call_OL(0,0, momenta_3, masses2_3, T0sum(26), M2)
  call TI_call_OL(0,0, momenta_4, masses2_2, T0sum(27), M2)

  call add_met(M2,M2L1R1)

#ifdef LOOPSQUARED
  end if
#endif

  call met_to_real(M2out,M2)

#ifdef PRECISION_dp
  call HOL_memory_deallocation_/**/REALKIND(1)
  call HCL_memory_deallocation_/**/REALKIND(1)
#endif

end subroutine integrate_tensor_sum

end module ol_tensor_sum_storage_ppllj_nexecsxg_1_/**/REALKIND
