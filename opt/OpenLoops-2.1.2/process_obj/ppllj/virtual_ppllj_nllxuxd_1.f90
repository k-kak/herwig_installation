! 1
! The number above is the number of vamp routines for this process.
! This is needed by the build system.



module ol_vamp_ppllj_nllxuxd_1_dp
  use kind_types, only: dp
  implicit none
  real(dp):: M2L0cache, M2L1cache
  contains

! **********************************************************************
  subroutine redbaseconstr()
! **********************************************************************
  use kind_types, only: dp
  use ol_data_types_dp, only: basis
  use ol_parameters_decl_qp ! masses
  use ol_tensor_sum_storage_ppllj_nllxuxd_1_dp
  use ofred_basis_construction_dp, only: construct_RedBasis, construct_redset4
  use ofred_basis_construction_dp, only: construct_redset5
  implicit none

  ! Define mass sets for reduction
  mass3set(0:2,1) = [0,0,0]


  ! Compute reduction bases from external momentum pairs p1, p2 (in light cone rep)
  call construct_RedBasis(8,11,RedBasis(1))


  ! Compute scalars depending on a third momentum p3 (in light cone rep)

end subroutine



! **********************************************************************
subroutine vamp(M)
! **********************************************************************
  use kind_types, only: dp, intkind2
  use ol_data_types_dp, only: Hpolcont
  use ol_parameters_decl_dp ! masses
  use ol_parameters_init_dp, only: init_met
  use ol_parameters_init_dp, only: channel_on
  use ol_external_ppllj_nllxuxd_1, only: & 
    & channel_number_ppllj_nllxuxd_1
  use ol_tensor_sum_storage_ppllj_nllxuxd_1_dp
  use ol_loop_storage_ppllj_nllxuxd_1_dp, only: &
    merge_step, hel_states

  use ol_vamp_1_ppllj_nllxuxd_1_dp, only: vamp_1

  implicit none
  type(Hpolcont), intent(in) :: M(1,hel_states)



  merge_step = 1
  ! Call subroutines for all branches to calculate loop coefficient tensors
  call channel_on(channel_number_ppllj_nllxuxd_1)
  call init_met(M2L1R1)
  call vamp_1(M)


end subroutine vamp

! **********************************************************************
subroutine vamp2base(P_scatt, M2L0, M2L1, IRL1, M2L2, IRL2, mode)
! P_scatt(0:3,Npart) = external momenta
! M2tree = helicity-summed squared tree matrix element for nu_tau tau+ anti-up down -> 0
! M2loop0 = helicity-summed squared loop matrix element for nu_tau tau+ anti-up down -> 0
! M2loop1 = IR1, M2loop2 = IR2 are dummy values for the single and double poles
! IR0, IR1, IR2 = finite, single pole, and double pole IR contribution
! mode = 1 (default): full matrix element;
!        2: reuse and scale coefficients from the last call;
!        note that scalings cannot be reset
! **********************************************************************
  use kind_types, only: dp, dp, intkind2

  use kind_types, only: qp, intkind2

  use ol_debug, only: ol_fatal
  use ol_generic, only: to_string
  use ol_external_ppllj_nllxuxd_1, only: &
    & external_perm_inv_ppllj_nllxuxd_1, &
    & average_factor_ppllj_nllxuxd_1, &
    & photonid_ppllj_nllxuxd_1
  use ol_colourmatrix_ppllj_nllxuxd_1_dp ! colmat_not_initialised, K1, K2, KL, KL2, KL2ct, KL2ct2
  use ol_kinematics_dp, only: init_kinematics, get_rmass2
  use ol_parameters_decl_dp ! parameters_status, scalefactor, <masses>
  ! tensorrankuse: for compatibility with old OL versions only insert if rank > 6
  use ol_parameters_init_dp, only: ensure_mp_loop_init
  use ol_init, only: set_parameter, parameters_flush
  use ol_loop_parameters_decl_dp, only: loop_parameters_status
  use ol_loop_parameters_decl_dp, only: IR_is_on, ioperator_mode, CT_is_on
  use ol_ew_renormalisation_dp, only: photon_factors
  use ol_momenta_decl_dp, only: momenta_nan_check
  use ol_momenta_decl_dp, only: QInvariantsMatrix
  use ol_tensor_sum_storage_ppllj_nllxuxd_1_dp, only: &
    & reset_tensor_sum, integrate_tensor_sum, scale_tensor_sum, HOL_memory_allocation_full, &
    HCL_memory_allocation, set_integral_masses_and_momenta, Tsum_memory_allocation, &
    HOL_memory_allocation_optimized
  use ol_tables_storage_ppllj_nllxuxd_1_dp, only: HOL_m3_init
  use ol_loop_ppllj_nllxuxd_1_dp, only: amp2
  use ol_loop_ppllj_nllxuxd_1_dp, only: init_merging_tables, update_merging_tables
  use ol_loop_storage_ppllj_nllxuxd_1_dp, only: &
    & fac_status_loop1, fac_status_loop2, M0_col1_helarray
  use ol_i_operator_dp, only: intdip
  use ol_loop_storage_ppllj_nllxuxd_1_dp, only: ntryL, nhel, p_switch, &
    & dp_not_alloc, qp_not_alloc
  use ol_settings_ppllj_nllxuxd_1, only: hel_mem_opt
  use hol_initialisation_dp, only: G0_hol_initialisation
  use ol_h_vert_interface_dp
  use ol_h_prop_interface_dp
  implicit none

  real(dp), intent(in)  :: P_scatt(0:3,4)
  real(dp),  intent(out) :: M2L0, M2L1, IRL1(0:2), M2L2, IRL2(0:4)
  integer, intent(in), optional :: mode
  real(dp)       :: P(0:3,4)
  integer              :: i, j, colmatpos, recycle_mode, ofred_mode
  complex(dp)    :: M0_col1(1), Mct(1)
  complex(dp)    :: Mcol_loop(1)
  real(dp)       :: M2colint(11), M2CC(4,4), M2CC_EW ! colour correlations
  integer              :: extmasses2(4)
  real(dp)       :: M2hel, M2ct, M2L2ct, M2L2ct2, vdip, c_dip(0:2)
  real(dp)       :: scalebackfactor
  real(dp)       :: bornphotonfactor = rONE , loopphotonfactor = rZERO
  integer(intkind2) :: helstates

  ofred_mode = 1

  if (present(mode)) then
    recycle_mode = mode
  else
    recycle_mode = 1
  end if

  call parameters_flush()
  call ensure_mp_loop_init()



  M2L0     = 0
  M2L1     = 0
  M2ct     = 0
  M2colint = 0
  M2L2     = 0
  M2L2ct   = 0
  M2L2ct2  = 0
  IRL1     = 0
  IRL2     = 0
  extmasses2 = [ 0, nML, 0, 0 ]

  if (momenta_nan_check(P_scatt) /= 0) return

  ! Here we calculate the tree level squared using the helicity bookkeeping subroutines
  ! The value of the tree level squared is M2L0
  ! The subtrees are also initialised
  call amp2(P_scatt, M2L0, .true., M2ct=M2ct, M2colint=M2colint)

  !The following part refers to the 1-loop calculation

  if (ntryL==1 .OR. dp_not_alloc) then
    call Tsum_memory_allocation()
  end if



  if (recycle_mode == 1 .or. ofred_mode == 1) call reset_tensor_sum()


  ! Construction of the basis used for the reduction
  call redbaseconstr()
  ! Memory allocation for the hol types is done just for the first phase space point
  ! Memory is reallocated when switching from double- to quad-precision for the first time

  if (ntryL==1 .OR. dp_not_alloc) then
    if(hel_mem_opt) then
      call HOL_memory_allocation_optimized()
    else
      call HOL_memory_allocation_full()
    end if
    call HCL_memory_allocation()
    dp_not_alloc = .FALSE.
  end if


  if(ntryL == 1) then
    !! m3 initialization done just for the first phase space point
    call HOL_m3_init()
    !! first temporary initialization of the merging tables
    call init_merging_tables(16,0)
  end if

  call set_integral_masses_and_momenta()

  !! Call for the calculation of the born-loop interference
  if (recycle_mode == 1 .or. ofred_mode == 1) then
    call vamp(M0_col1_helarray)
  end if

  !! Merging tables are adjusted and allocated with the proper dimensionality
  if(ntryL == 1) call update_merging_tables(16)

  if(ntryL==1) ntryL = ntryL + 1


  p_switch = 1




  if (recycle_mode == 2 .and. ofred_mode == 0) then
    call scale_tensor_sum()
  end if
  call integrate_tensor_sum(M2L1)


  if (IR_is_on > 0) then
    do i = 1, 4
      do j = 1, i
        ! Why does this work without permuting the colour correlation matrices?
        M2CC(i,j) = M2colint(i*(i-1)/2+j)
      end do
    end do
    do j = 2, 4
      do i = 1, j-1
        M2CC(i,j) = M2CC(j,i)
      end do
    end do
    M2CC_EW = M2colint(4*(4+1)/2+1)
    call intdip(ioperator_mode, M2L0, M2CC, M2CC_EW, [0,3,2,2], [0,3,-2,-1]/3._dp, &
      & 4, get_rmass2(extmasses2), QInvariantsMatrix, vdip, c_dip, &
      & 1, &
      & [( photonid_ppllj_nllxuxd_1( &
      & external_perm_inv_ppllj_nllxuxd_1(i)), &
      & i=1, size(external_perm_inv_ppllj_nllxuxd_1))])
    IRL1(0) = c_dip(0) / average_factor_ppllj_nllxuxd_1
    IRL1(1) = c_dip(1) / average_factor_ppllj_nllxuxd_1
    IRL1(2) = c_dip(2) / average_factor_ppllj_nllxuxd_1
  else
    vdip = 0
    IRL1 = 0
  end if

  ! loop^2 IR contribution: not implemented
  IRL2 = 0

  ! photon factors
  call photon_factors(photonid_ppllj_nllxuxd_1, &
       & 0, bornphotonfactor, loopphotonfactor)
  M2L0 = bornphotonfactor * M2L0
  M2L1 = bornphotonfactor * M2L1
  IRL1 = bornphotonfactor * IRL1
  M2L2 = bornphotonfactor * M2L2
  IRL2 = bornphotonfactor * IRL2
  M2ct = bornphotonfactor * M2ct
  vdip = bornphotonfactor * vdip
  if (CT_is_on .gt. 0) M2ct = M2ct+M2L0*loopphotonfactor

  ! Colour and helicity average and symmetry factor of outgoing particles
  M2L0 = M2L0 / average_factor_ppllj_nllxuxd_1
  M2L1 = 2*(M2L1 + M2ct)
  if (IR_is_on > 1) then
    M2L1 = M2L1 + vdip
  end if
  M2L1 = M2L1 / average_factor_ppllj_nllxuxd_1

  ! check for NaN result
  if (M2L0 /= M2L0) then
    M2L0 = 0
  end if
  if (M2L1 /= M2L1) then
    M2L1 = 0
    M2L2 = 0
    IRL1 = 0
    IRL2 = 0
  end if

  scalebackfactor = scalefactor**(2*4-8)
  M2L0 = scalebackfactor * M2L0
  M2L1 = scalebackfactor * M2L1
  IRL1 = scalebackfactor * IRL1
  M2L2 = scalebackfactor * M2L2
  IRL2 = scalebackfactor * IRL2
  M2L0cache = M2L0
  M2L1cache = M2L1

end subroutine vamp2base


! **********************************************************************
subroutine ctamp2base(P_scatt, M2tree, M2ct)
! The part of vamp2 which calculates tree and counter-term
! matrix elements, but not the loop. R2 is deactivated.
! Does not calculate loop^2 counterterms/R2.
! P_scatt(0:3,Npart) = external momenta
! M2tree = helicity-summed squared tree matrix element for nu_tau tau+ anti-up down -> 0
! M2ct   = helicity-summed counterterm matrix element
! **********************************************************************
  use kind_types, only: dp, dp
  use ol_external_ppllj_nllxuxd_1, only: &
    & photonid_ppllj_nllxuxd_1
  use ol_colourmatrix_ppllj_nllxuxd_1_dp, only: &
    & colmat_not_initialised, K1, K2
  use ol_parameters_decl_dp ! parameters_status, scalefactor, <masses>
  use ol_parameters_decl_dp, only: do_ew_renorm
  use ol_loop_parameters_decl_dp, only: &
    & loop_parameters_status, CT_is_on, R2_is_on, TP_is_on
  use ol_ew_renormalisation_dp, only: photon_factors
  use ol_init, only: set_parameter, parameters_flush
  use ol_parameters_init_dp, only: ensure_mp_loop_init
  use ol_loop_parameters_decl_dp, only: DOI
  use ol_external_ppllj_nllxuxd_1, only: &
    & external_perm_inv_ppllj_nllxuxd_1, average_factor_ppllj_nllxuxd_1
  use ol_momenta_decl_dp, only: momenta_nan_check
  use ol_kinematics_dp, only: init_kinematics
  use ol_loop_storage_ppllj_nllxuxd_1_dp, only: &
    & fac_status_loop1, fac_status_loop2
  use ol_loop_ppllj_nllxuxd_1_dp, only: amp2
  implicit none

  real(dp), intent(in)  :: P_scatt(0:3,4)
  real(dp),  intent(out) :: M2tree, M2ct
  real(dp)    :: P(0:3,4)
  real(dp)    :: M2colint(11)
  integer           :: CT_on_bak, R2_on_bak, TP_on_bak
  real(dp)    :: scalebackfactor, bornphotonfactor = rONE, loopphotonfactor = rZERO
  integer           :: DOI_bak

  M2tree = 0
  M2ct   = 0

  if (CT_is_on == 2) then
    call set_parameter("ct_on", 1)
  else if (CT_is_on == -1) then
    call set_parameter("ct_on", 0)
  end if
  DOI_bak = DOI
  DOI = 0

  if (momenta_nan_check(P_scatt) /= 0)  return

  ! Here we calculate the tree level squared using the helicity bookkeeping subroutines
  call amp2(P_scatt, M2tree, .false., M2ct=M2ct)

  ! photon factors
  call photon_factors(photonid_ppllj_nllxuxd_1, 0, bornphotonfactor, loopphotonfactor)
  M2tree = bornphotonfactor * M2tree
  M2ct = bornphotonfactor * M2ct
  if (CT_is_on .gt. 0) M2ct = M2ct+M2tree*loopphotonfactor

  DOI = DOI_bak

  scalebackfactor = scalefactor**(2*4-8)
  M2tree = scalebackfactor * M2tree / average_factor_ppllj_nllxuxd_1
  M2ct   = scalebackfactor * 2*M2ct / average_factor_ppllj_nllxuxd_1

end subroutine ctamp2base



! **********************************************************************
subroutine iopamp2base(P_scatt, M2L0, IRL1)
! P_scatt(0:3,Npart) = external momenta
! IR0, IR1, IR2 = finite, single pole, and double pole IR contribution
! **********************************************************************
  use kind_types, only: dp, dp, intkind2

  use kind_types, only: qp, intkind2

  use ol_debug, only: ol_fatal
  use ol_generic, only: to_string
  use ol_external_ppllj_nllxuxd_1, only: &
    & external_perm_inv_ppllj_nllxuxd_1, &
    & average_factor_ppllj_nllxuxd_1, &
    & photonid_ppllj_nllxuxd_1
  use ol_colourmatrix_ppllj_nllxuxd_1_dp ! colmat_not_initialised, K1, K2, KL, KL2, KL2ct, KL2ct2
  use ol_kinematics_dp, only: init_kinematics, get_rmass2
  use ol_parameters_decl_dp ! parameters_status, scalefactor, <masses>
  ! tensorrankuse: for compatibility with old OL versions only insert if rank > 6
  use ol_parameters_init_dp, only: ensure_mp_loop_init
  use ol_init, only: set_parameter, parameters_flush
  use ol_loop_parameters_decl_dp, only: loop_parameters_status
  use ol_loop_parameters_decl_dp, only: IR_is_on, ioperator_mode
  use ol_ew_renormalisation_dp, only: photon_factors
  use ol_momenta_decl_dp, only: momenta_nan_check
  use ol_momenta_decl_dp, only: QInvariantsMatrix
  use ol_loop_ppllj_nllxuxd_1_dp, only: amp2
  use ol_i_operator_dp, only: intdip
  use ol_loop_storage_ppllj_nllxuxd_1_dp, only: &
    & fac_status_loop1, fac_status_loop2
  implicit none

  real(dp), intent(in)  :: P_scatt(0:3,4)
  real(dp),  intent(out) :: M2L0, IRL1(0:2)
  real(dp)       :: P(0:3,4)
  integer              :: i, j, colmatpos
  complex(dp)    :: M0_col1(1), Mct(1)
  complex(dp)    :: Mcol_loop(1)
  real(dp)       :: M2colint(11), M2CC(4,4), M2CC_EW ! colour correlations
  integer              :: extmasses2(4)
  real(dp)       :: M2ct, M2hel, vdip, c_dip(0:2)
  real(dp)       :: scalebackfactor
  real(dp)       :: bornphotonfactor = rONE , loopphotonfactor = rZERO
  integer(intkind2) :: helstates

  M2L0     = 0
  M2ct     = 0
  M2colint = 0
  IRL1     = 0
  extmasses2 = [ 0, nML, 0, 0 ]

  if (momenta_nan_check(P_scatt) /= 0) return

  ! Here we calculate the tree level squared using the helicity bookkeeping subroutines
  call amp2(P_scatt, M2L0, .false., M2colint=M2colint)

  if (IR_is_on > 0) then
    do i = 1, 4
      do j = 1, i
        ! Why does this work without permuting the colour correlation matrices?
        M2CC(i,j) = M2colint(i*(i-1)/2+j)
      end do
    end do
    do j = 2, 4
      do i = 1, j-1
        M2CC(i,j) = M2CC(j,i)
      end do
    end do
    M2CC_EW = M2colint(4*(4+1)/2+1)
    call intdip(ioperator_mode, M2L0, M2CC, M2CC_EW, [0,3,2,2], [0,3,-2,-1]/3._dp, &
      & 4, get_rmass2(extmasses2), QInvariantsMatrix, vdip, c_dip, &
      & 1, &
      & [( photonid_ppllj_nllxuxd_1( &
      & external_perm_inv_ppllj_nllxuxd_1(i)), &
      & i=1, size(external_perm_inv_ppllj_nllxuxd_1))])
    IRL1(0) = c_dip(0) / average_factor_ppllj_nllxuxd_1
    IRL1(1) = c_dip(1) / average_factor_ppllj_nllxuxd_1
    IRL1(2) = c_dip(2) / average_factor_ppllj_nllxuxd_1
  else
    vdip = 0
    IRL1 = 0
  end if

  ! Colour and helicity average and symmetry factor of outgoing particles
  M2L0 = M2L0 / average_factor_ppllj_nllxuxd_1

  ! photon factors
  call photon_factors(photonid_ppllj_nllxuxd_1, &
       & 0, bornphotonfactor, loopphotonfactor)
  M2L0 = bornphotonfactor * M2L0
  IRL1 = bornphotonfactor * IRL1

  ! check for NaN result
  if (M2L0 /= M2L0) then
    M2L0 = 0
    IRL1 = 0
  end if

  scalebackfactor = scalefactor**(2*4-8)
  M2L0 = scalebackfactor * M2L0
  IRL1 = scalebackfactor * IRL1

end subroutine iopamp2base



! **********************************************************************
subroutine vamp2cache(P_scatt, M2L0, M2L1, corr)
! **********************************************************************
  use kind_types, only: dp, dp, intkind2
  use ol_parameters_decl_dp ! parameters_status, scalefactor, <masses>
  use ol_loop_parameters_decl_dp, only: CT_is_on
  use ol_external_ppllj_nllxuxd_1, only: & 
    & extcomb_perm_ppllj_nllxuxd_1, &
    & average_factor_ppllj_nllxuxd_1, &
    & photonid_ppllj_nllxuxd_1
  use ol_loop_storage_ppllj_nllxuxd_1_dp, only: M2ctcc, M0M1_hel_cc
  use ol_tensor_sum_storage_ppllj_nllxuxd_1_dp, only: reset_tensor_sum, &
    & integrate_tensor_sum
  use ol_debug, only: ol_fatal
  use ol_colourmatrix_ppllj_nllxuxd_1_dp, only: Cas
  use ol_data_types_dp, only: correlator
  use ol_ew_renormalisation_dp, only: photon_factors
  implicit none
  real(dp), intent(in)  :: P_scatt(0:3,4)
  real(dp),  intent(out) :: M2L0, M2L1
  real(dp) :: M2L1cc(0:11)
  real(dp)       :: bornphotonfactor = rONE , loopphotonfactor = rZERO
  type(correlator) :: corr

  integer :: co, k, m, ind_cc_comb(2)

  ind_cc_comb = [2,4]
  M2L1cc = 0
  !! Diagonal born-loop colour correlators
  co = 0
  do k = 1, 4
    co = co + k
    M2L1cc(co) = Cas(k)*M2L1cache
  end do

  ! photon factors
  call photon_factors(photonid_ppllj_nllxuxd_1, &
       & 0, bornphotonfactor, loopphotonfactor)
  if (CT_is_on .gt. 0) M2ctcc = M2ctcc+M2L0cache*loopphotonfactor
  !! Independent colour correlators
  do k = 1, 2
    m = ind_cc_comb(k)
    call reset_tensor_sum()
    M2L1 = 0
    if(sum(abs(M0M1_hel_cc(:,:,m)%j)) == 0) then
      M2L1 = 0
    else
      call vamp(M0M1_hel_cc(:,:,m))
      call integrate_tensor_sum(M2L1)
    end if
    M2L1 = bornphotonfactor * M2L1
    M2L1cc(m) = 2*(M2L1 + M2ctcc(m))/average_factor_ppllj_nllxuxd_1
  end do

  !! Linearly dependent colour correlators
  M2L1cc(7) = -(M2L1cache*Cas(1)) - M2L1cc(2) - M2L1cc(4)
  M2L1cc(5) = -(M2L1cache*(Cas(1) + Cas(2) + Cas(3) - Cas(4)))/2 - M2L1cc(2) - M2L1cc(4)
  M2L1cc(8) = (M2L1cache*(Cas(1) - Cas(2) + Cas(3) - Cas(4)))/2 + M2L1cc(4)
  M2L1cc(9) = (M2L1cache*(Cas(1) + Cas(2) - Cas(3) - Cas(4)))/2 + M2L1cc(2)
  
  do k = 1, 11
    corr%rescc(k) = M2L1cc(extcomb_perm_ppllj_nllxuxd_1(k))
  end do
  M2L0 = M2L0cache
  corr%rescc(0) = M2L1cache
  M2L1 = M2L1cache

end subroutine vamp2cache

! **********************************************************************
subroutine schsfamp2base(P_scatt, M2tree, i, j, M2hsf)
! P_scatt(0:3,Npart) = external momenta
! M2tree = helicity-summed squared tree matrix element for
!          nu_tau tau+ anti-up down -> 0
! M2hsf  = spin correlated hard scattering factor. Eq. (66) in 1011.3918
! **********************************************************************
  use kind_types, only: dp, dp
  use ol_parameters_decl_dp, only: scalefactor
  use ol_parameters_init_dp, only: ensure_mp_init
  use ol_external_ppllj_nllxuxd_1, only: &
 &  average_factor_ppllj_nllxuxd_1, &
 &  photonid_ppllj_nllxuxd_1
  use ol_momenta_decl_dp, only: momenta_nan_check
  use ol_loop_ppllj_nllxuxd_1_dp, only: hsfamp2
  use ol_ew_renormalisation_dp, only: photon_factors
  implicit none
  real(dp), intent(in) :: P_scatt(0:3,4)
  integer, intent(in) :: i, j
  real(dp), intent(out) :: M2tree, M2hsf
  real(dp) :: scalebackfactor, bornphotonfactor
  M2tree = 0
  M2hsf = 0
  call ensure_mp_init()
  if (momenta_nan_check(P_scatt) /= 0)  return
  call hsfamp2(P_scatt, M2tree, i, j, M2hsf)

  ! photon-factors
  call photon_factors(photonid_ppllj_nllxuxd_1, 0, bornphotonfactor)
  M2tree = bornphotonfactor * M2tree
  M2hsf = bornphotonfactor * M2hsf

  scalebackfactor = scalefactor**(2*4-8)
  M2tree = scalebackfactor * M2tree / average_factor_ppllj_nllxuxd_1
  M2hsf  = scalebackfactor * M2hsf / average_factor_ppllj_nllxuxd_1
end subroutine schsfamp2base

subroutine vamp2pc(P_scatt, M2L0, M2L1, IRL1, M2L2, IRL2, mode)
  ! polecheck routine: perform multiple calls of vamp2
  ! with different values of the UV/IR poles
  ! to determine the coefficients of the IR poles.
  use kind_types, only: dp, dp
  use ol_parameters_decl_dp, only: rZERO_dp => rZERO
  use ol_loop_parameters_decl_dp, only: &
    & de1_UV, de1_IR, de2_i_IR
  use ol_loop_parameters_decl_dp, only: &
    & polescale_dp => polescale, polecheck_is, CT_is_on, IR_is_on
  use ol_init, only: set_parameter
  implicit none
  real(dp), intent(in)  :: P_scatt(:,:)
  real(dp),  intent(out) :: M2L0, M2L1(0:2), IRL1(0:2), M2L2(0:4), IRL2(0:4)
  integer, intent(in), optional :: mode
  real(dp) :: poleUV1bak, poleIR1bak, poleIR2bak
  integer         :: CT_on_bak, IR_on_bak
  real(dp)  :: M2L1_00, M2L1_10, M2L1_01, M2L1_20, M2L1_02, M2L1_11
  real(dp)  :: M2L2_00, M2L2_10, M2L2_01, M2L2_20, M2L2_02, M2L2_11
  real(dp)  :: IRL1x(0:2), IRL2x(0:4)

  if (polecheck_is == 0) then
    call vamp2base(P_scatt, M2L0, M2L1(0), IRL1, M2L2(0), IRL2, mode)
    if (IR_is_on == 2) then
      M2L1(1:2) = 0
      M2L2(1:4) = 0
    else
      M2L1(1:2) = -IRL1(1:2)
      M2L2(1:4) = -IRL2(1:4)
    end if
  else
    ! remember original parameters
    poleUV1bak = de1_UV
    poleIR1bak = de1_IR
    poleIR2bak = de2_i_IR
    CT_on_bak  = CT_is_on
    IR_on_bak  = IR_is_on
    ! three calls with different IR poles; IR subtractions are calculated but not added
    call set_parameter("pole_uv", rZERO_dp)
    call set_parameter("pole_ir1", rZERO_dp)
    call set_parameter("pole_ir2", rZERO_dp)
    call set_parameter("ct_on", 1)
    call set_parameter("ir_on", 1)
    call vamp2base(P_scatt, M2L0, M2L1_00, IRL1, M2L2_00, IRL2, mode)
    call set_parameter("pole_uv", polescale_dp)
    call set_parameter("pole_ir1", polescale_dp)
    call set_parameter("pole_ir2", rZERO_dp)
    call set_parameter("ir_on", 0)
    call vamp2base(P_scatt, M2L0, M2L1_10, IRL1x, M2L2_10, IRL2x, mode)
    call set_parameter("pole_uv", rZERO_dp)
    call set_parameter("pole_ir1", rZERO_dp)
    call set_parameter("pole_ir2", polescale_dp)
    call vamp2base(P_scatt, M2L0, M2L1_01, IRL1x, M2L2_01, IRL2x, mode)

    M2L1(0) = M2L1_00
    M2L1(1) = (M2L1_10 - M2L1_00)/polescale_dp
    M2L1(2) = (M2L1_01 - M2L1_00)/polescale_dp


    ! restore original parameters
    call set_parameter("pole_uv", poleUV1bak)
    call set_parameter("pole_ir1", poleIR1bak)
    call set_parameter("pole_ir2", poleIR2bak)
    call set_parameter("ct_on", CT_on_bak)
    call set_parameter("ir_on", IR_on_bak)
    if (IR_on_bak == 2) then
      M2L1 = M2L1 + IRL1
      M2L2 = M2L2 + IRL2
    end if
  end if

end subroutine vamp2pc

end module ol_vamp_ppllj_nllxuxd_1_dp



module ol_vamp_ppllj_nllxuxd_1_qp
  implicit none
  contains
  subroutine vamp2pc(P_scatt, M2L0, M2L1, IRL1, M2L2, IRL2, mode)
    ! Dummy quad precision routine. Throws an error if the process called but
    ! the process has not been compiled in quad precision.
    use kind_types, only: dp, qp
    use ol_debug, only: ol_fatal
    implicit none
    real(dp), intent(in)  :: P_scatt(:,:)
    real(qp),  intent(out) :: M2L0, M2L1(0:2), IRL1(0:2), M2L2(0:4), IRL2(0:4)
    integer, intent(in), optional :: mode
    call ol_fatal("Process ppllj_nllxuxd_1" &
      & // "has not been compiled in quad precision")
  end subroutine vamp2pc
end module ol_vamp_ppllj_nllxuxd_1_qp





module ol_vamp_ppllj_nllxuxd_1
  use ol_parameters_decl_dp, only: procname_length
  use ol_data_types_dp, only: me_cache
  implicit none
  character(procname_length) :: processname = 'ppllj_nllxuxd_1'
  integer, save :: qp_eval = 0, killed = 0
  integer, save :: npoints(8) = 0
  integer, save :: stability_histogram(20) = 0, stability_histogram_qp(20) = 0
  type(me_cache), allocatable, target, save :: me_caches(:)
  contains

  subroutine finish_ppllj_nllxuxd_1()
    ! final update of the stability histogram and memory deallocation
    use kind_types, only: dp
!   use ol_data_types_dp, only: Hpolcont
!   use ol_tensor_sum_storage_ppllj_nllxuxd_1_dp, only: HOL_memory_allocation
    use ol_stability, only: finish_histograms
    use ol_vamp_ppllj_nllxuxd_1_dp, only: vamp
    use ol_tensor_sum_storage_ppllj_nllxuxd_1_dp, only: &
      & HOL_memory_deallocation_dp, HCL_memory_deallocation_dp

    use ol_loop_storage_ppllj_nllxuxd_1_dp, only: &
      & dp_not_alloc, qp_not_alloc, merge_tables, merge_mism, merge_hels, merge_tables_on, A, &
      & M1helarray, M1helarray_ct, M0_col1_helarray, M0M1_hel_cc, hflip
      implicit none
      integer :: k
!   type(Hpolcont) :: Mdummy(1,16)
!   call HOL_memory_allocation()
!   call vamp(Mdummy)

    ! Memory is deallocated for both dp and/or qp stored hol coefficients
    if (.NOT. dp_not_alloc) then
      call HOL_memory_deallocation_dp(0)
      call HCL_memory_deallocation_dp(0)
    end if


    if (merge_tables_on) then
      if (allocated(merge_tables)) deallocate(merge_tables)
      if (allocated(merge_mism)) deallocate(merge_mism)
      if (allocated(merge_hels)) deallocate(merge_hels)
      merge_tables_on = .false.
    end if

    if (allocated(A)) deallocate(A)
    if (allocated(M1helarray)) deallocate(M1helarray)
    if (allocated(M1helarray_ct)) deallocate(M1helarray_ct)
    if (allocated(M0_col1_helarray)) deallocate(M0_col1_helarray)
    if (allocated(M0M1_hel_cc)) deallocate(M0M1_hel_cc)
    if (allocated(hflip)) deallocate(hflip)

    if (allocated(me_caches)) then
      do k = 1, size(me_caches)
      if (allocated(me_caches(k)%psp)) deallocate(me_caches(k)%psp)
      if (allocated(me_caches(k)%me)) deallocate(me_caches(k)%me)
    end do
      if (allocated(me_caches)) deallocate(me_caches)
    end if
    call finish_histograms(processname, stability_histogram, stability_histogram_qp, npoints, qp_eval, killed)
  end subroutine finish_ppllj_nllxuxd_1


  subroutine vamp2(P_scatt, M2L0, M2L1, IRL1, M2L2, IRL2) &
      & bind(c,name="ol_f_vamp2_ppllj_nllxuxd_1")
    use kind_types, only: dp
    use ol_init, only: register_cleanup
    use ol_external_ppllj_nllxuxd_1, only: &
      & external_perm_ppllj_nllxuxd_1
    use ol_vamp_ppllj_nllxuxd_1_dp, only: vamp2dp => vamp2pc
    use ol_vamp_ppllj_nllxuxd_1_qp, only: vamp2qp => vamp2pc
    use ol_stability, only: vamp2generic
    implicit none
    real(dp), intent(in)  :: P_scatt(0:3,4)
    real(dp), intent(out) :: M2L0, M2L1(0:2), IRL1(0:2), M2L2(0:4), IRL2(0:4)
    real(dp), save :: abs_kfactor_threshold = 1, trigeff_local = 0, sum_M2tree = 0
    logical :: first_call = .true.
    if (first_call) then
      call register_cleanup(finish_ppllj_nllxuxd_1)
      first_call = .false.
    end if
    call vamp2generic(vamp2dp, vamp2qp, processname, P_scatt, M2L0, M2L1, IRL1, M2L2, IRL2, &
                    & abs_kfactor_threshold, trigeff_local, sum_M2tree, &
                    & npoints, qp_eval, killed, stability_histogram, stability_histogram_qp, &
                    & external_perm_ppllj_nllxuxd_1, me_caches)
  end subroutine vamp2


  subroutine ctamp2(P_scatt, M2tree, M2ct) &
      & bind(c,name="ol_f_ctamp2_ppllj_nllxuxd_1")
    use kind_types, only: dp
    use ol_vamp_ppllj_nllxuxd_1_dp, only: ctamp2base
    implicit none
    real(dp), intent(in)  :: P_scatt(0:3,4)
    real(dp), intent(out) :: M2tree, M2ct
    call ctamp2base(P_scatt, M2tree, M2ct)
  end subroutine ctamp2


  subroutine iopamp2(P_scatt, M2tree, IRL1) &
      & bind(c,name="ol_f_iopamp2_ppllj_nllxuxd_1")
    use kind_types, only: dp
    use ol_vamp_ppllj_nllxuxd_1_dp, only: iopamp2base
    implicit none
    real(dp), intent(in)  :: P_scatt(0:3,4)
    real(dp), intent(out) :: M2tree, IRL1(0:2)
    call iopamp2base(P_scatt, M2tree, IRL1)
  end subroutine iopamp2


  !Spin-Correlated Hard Scattering Factor
  subroutine schsfamp2(P_scatt, M2tree, i, j, M2hsf) &
      & bind(c,name="ol_f_schsfamp2_ppllj_nllxuxd_1")
    use kind_types, only: dp
    use ol_vamp_ppllj_nllxuxd_1_dp, only: schsfamp2base
    implicit none
    real(dp), intent(in)  :: P_scatt(0:3,4)
    integer, intent(in)  :: i, j
    real(dp), intent(out) :: M2tree, M2hsf
    call schsfamp2base(P_scatt, M2tree, i, j, M2hsf)
  end subroutine schsfamp2

  subroutine vamp2_c(P_scatt, M2L0, M2L1, IRL1, M2L2, IRL2) &
      & bind(c,name="ol_vamp2_ppllj_nllxuxd_1")
    use kind_types, only: dp
    use, intrinsic :: iso_c_binding, only: c_double
    implicit none
    real(c_double), intent(in)  :: P_scatt(0:3,4)
    real(c_double), intent(out) :: M2L0, M2L1(0:2), IRL1(0:2), M2L2(0:4), IRL2(0:4)
    real(dp) :: f_p_scatt(0:3,4)
    real(dp) :: f_m2l0, f_m2l1(0:2), f_irl1(0:2), f_m2l2(0:4), f_irl2(0:4)
    f_p_scatt = P_scatt
    call vamp2(f_p_scatt, f_m2l0, f_m2l1, f_irl1, f_m2l2, f_irl2)
    M2L0 = f_m2l0
    M2L1 = f_m2l1
    IRL1 = f_irl1
    M2L2 = f_m2l2
    IRL2 = f_irl2
  end subroutine vamp2_c


  subroutine ctamp2_c(P_scatt, M2tree, M2ct) &
      & bind(c,name="ol_ctamp2_ppllj_nllxuxd_1")
    use kind_types, only: dp
    use, intrinsic :: iso_c_binding, only: c_double
    implicit none
    real(c_double), intent(in)  :: P_scatt(0:3,4)
    real(c_double), intent(out) :: M2tree, M2ct
    real(dp) :: f_p_scatt(0:3,4)
    real(dp) :: f_m2tree, f_m2ct
    f_p_scatt = P_scatt
    call ctamp2(f_p_scatt, f_m2tree, f_m2ct)
    M2tree = f_m2tree
    M2ct = f_m2ct
  end subroutine ctamp2_c

  !!! Interface for tree-loop correlators.
  subroutine vampcr2(P_scatt, M2L0, M2L1, M2L2, crtype, &
      & emitter, nextcombs, extcombs, mom, M2cc, M2munu) &
      & bind(c,name="ol_f_vampcr2_ppllj_nllxuxd_1")
    use kind_types, only: dp
    use ol_debug, only: ol_error
    use ol_vamp_ppllj_nllxuxd_1_dp, only: vamp2cache
    use ol_data_types_dp, only: correlator
    implicit none
    real(dp), intent(in)  :: P_scatt(0:3,4), mom(:)
    real(dp), intent(out) :: M2L0, M2L1, M2L2, M2cc(0:11)
    integer, intent(in) :: crtype, emitter, nextcombs, extcombs(nextcombs)
    real(dp), intent(out) :: M2munu(:,:)
    type(correlator) :: corr
    M2L0 = 0
    M2L1 = 0
    M2L2 = 0
    corr%type=crtype
    if(crtype == 11) then
      allocate(corr%rescc(0:11))
      call vamp2cache(P_scatt, M2L0, M2L1, corr=corr)
      M2cc = corr%rescc
      deallocate(corr%rescc)
    else
      call ol_error(1,"Correlator not available")
    end if
  end subroutine vampcr2

end module ol_vamp_ppllj_nllxuxd_1

! #ifdef 1

