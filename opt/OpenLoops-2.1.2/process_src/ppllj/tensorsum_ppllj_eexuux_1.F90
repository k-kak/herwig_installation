
! **********************************************************************
module ol_tables_storage_ppllj_eexuux_1_/**/REALKIND
! **********************************************************************
  use KIND_TYPES, only: REALKIND, QREALKIND, intkind2
  use ol_data_types_/**/REALKIND, only: hol
  use ol_data_types_/**/REALKIND, only: basis, redset4, redset5
  implicit none

  ! helicity tables for the 1-loop recursion
integer(intkind2), save :: h0tab(16,2)
integer(intkind2), save :: heltab2x2(2,2,1)
integer(intkind2), save :: heltab2x8(2,8,2)
integer(intkind2), save :: heltab2x16(2,16,2)


  ! number of helicity states for openloops recursion steps
integer(intkind2), save :: m0h(2)
integer(intkind2), save :: m3h2x1(3,1)
integer(intkind2), save :: m3h4x2(3,2)
integer(intkind2), save :: m3h2x8(3,2)

integer(intkind2), save :: n2h2(2)
integer(intkind2), save :: n2h8(2)


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
m3h4x2(1,:)=4
m3h4x2(2,:)=2
m3h4x2(3,:)=8
m3h2x8(1,:)=2
m3h2x8(2,:)=8
m3h2x8(3,:)=16

n2h2(:)=2
n2h8(:)=8


end subroutine HOL_m3_init

!**********************************************************************

end module ol_tables_storage_ppllj_eexuux_1_/**/REALKIND

! **********************************************************************
module ol_tensor_sum_storage_ppllj_eexuux_1_/**/REALKIND
! **********************************************************************
  use KIND_TYPES, only: REALKIND, QREALKIND, intkind2
  use ol_data_types_/**/REALKIND, only: hol, hcl, met
  use ol_data_types_/**/REALKIND, only: basis, redset4, redset5
  implicit none

  type(met), save :: M2L1R1

  ! Declarations of loop wave function tensors

  type(hol), save :: G1H2(1)
  type(hol), save :: G2H2(2)
  type(hol), save :: G0H8(1)
  type(hol), save :: G1H8(2)
  type(hol), save :: G0H16(1)
  type(hcl), save, dimension(4) :: G0tensor
  type(hcl), save, dimension(1) :: G2tensor



  ! Declarations for on-the-fly tensor reduction
type (basis),      save :: RedBasis(1)
integer, save :: mass3set(0:2,1)



  ! Declarations for TI calls

  integer, save :: momenta_1(2)
  integer, save :: momenta_2(2)
  integer, save :: momenta_3(2)
  integer, save :: momenta_4(3)

  integer, save :: masses2_1(2)
  integer, save :: masses2_2(3)



  type(hcl), save, dimension(4) :: T0sum


  contains

!**********************************************************************
subroutine HOL_memory_allocation_full()
!**********************************************************************
! allocation of memory for the types hol
!**********************************************************************
  use KIND_TYPES, only: REALKIND
  use ol_loop_storage_ppllj_eexuux_1_/**/DREALKIND, only: nhel
  use ol_data_types_/**/REALKIND, only: hol
  use hol_initialisation_/**/REALKIND, only: hol_allocation
  implicit none

    call hol_allocation(4,5,4,2,G1H2,1)
  call hol_allocation(4,15,4,2,G2H2,2)
  call hol_allocation(4,1,4,8,G0H8,1)
  call hol_allocation(4,5,4,8,G1H8,2)
  call hol_allocation(4,1,4,16,G0H16,1)


end subroutine HOL_memory_allocation_full

!**********************************************************************
subroutine HOL_memory_allocation_optimized()
!**********************************************************************
! allocation of memory for the types hol
!**********************************************************************
  use KIND_TYPES, only: REALKIND
  use ol_loop_storage_ppllj_eexuux_1_/**/DREALKIND, only: nhel
  use ol_data_types_/**/REALKIND, only: hol
  use hol_initialisation_/**/REALKIND, only: hol_allocation
  implicit none

    call hol_allocation(4,5,4,min(nhel,2),G1H2,1)
  call hol_allocation(4,15,4,min(nhel,2),G2H2,2)
  call hol_allocation(4,1,4,min(nhel,8),G0H8,1)
  call hol_allocation(4,5,4,min(nhel,8),G1H8,2)
  call hol_allocation(4,1,4,min(nhel,16),G0H16,1)


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

    call hol_deallocation(G1H2,1,dmode)
  call hol_deallocation(G2H2,2,dmode)
  call hol_deallocation(G0H8,1,dmode)
  call hol_deallocation(G1H8,2,dmode)
  call hol_deallocation(G0H16,1,dmode)


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

  call hcl_allocation(1,G0tensor, 4)
call hcl_allocation(15,G2tensor, 1)


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

  call hcl_deallocation(G0tensor, 4,dmode)
call hcl_deallocation(G2tensor, 1,dmode)

    call hcl_deallocation(T0sum,4,dmode)


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

    call hcl_allocation(1,T0sum,4)


end subroutine Tsum_memory_allocation


#ifdef PRECISION_dp
subroutine max_point(r) &
    & bind(c,name="ol_f_max_point_ppllj_eexuux_1")
  ! Return the number maximal tensor rank
  implicit none
  integer, intent(out) :: r
  r = 3
end subroutine max_point

subroutine tensor_rank(r) &
    & bind(c,name="ol_f_tensor_rank_ppllj_eexuux_1")
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

  do i = 1,4
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
  call scale_one_tsum(T0sum(2)%cmp, 0)
  call scale_one_tsum(T0sum(3)%cmp, 0)
  call scale_one_tsum(T0sum(4)%cmp, 0)

end subroutine scale_tensor_sum

! **********************************************************************
subroutine set_integral_masses_and_momenta()
! **********************************************************************

  use ol_parameters_decl_/**/REALKIND
  momenta_1 = [ 8, 7 ]
  momenta_2 = [ 11, 4 ]
  momenta_3 = [ 12, 3 ]
  momenta_4 = [ 8, 3, 4 ]

  masses2_1 = [ 0, 0 ]
  masses2_2 = [ 0, 0, 0 ]


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
  call TI_call_OL(0,0, momenta_4, masses2_2, T0sum(1), M2)
  call TI_call_OL(0,0, momenta_3, masses2_1, T0sum(2), M2)
  call TI_call_OL(0,0, momenta_2, masses2_1, T0sum(3), M2)
  call TI_call_OL(0,0, momenta_1, masses2_1, T0sum(4), M2)

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

end module ol_tensor_sum_storage_ppllj_eexuux_1_/**/REALKIND
