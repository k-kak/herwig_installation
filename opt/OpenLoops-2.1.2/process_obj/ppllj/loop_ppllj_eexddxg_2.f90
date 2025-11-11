
module ol_colourmatrix_ppllj_eexddxg_2_dp
  use kind_types, only: dp
  implicit none
  logical, save           :: colmat_not_initialised = .true.
  complex(dp), save :: K1(17,1), K2(17,1), KL(17,1), Cas(5) = 0
  contains
  subroutine colourmatrix_init
    use ol_parameters_decl_dp, only: CI
    implicit none
    integer :: k, co
    colmat_not_initialised = .false.
    ! colour matrix

  K1( 1,:) = [  12]
  K1( 2,:) = [   0]
  K1( 3,:) = [   0]
  K1( 4,:) = [   0]
  K1( 5,:) = [   0]
  K1( 6,:) = [   0]
  K1( 7,:) = [  16]
  K1( 8,:) = [   0]
  K1( 9,:) = [   0]
  K1(10,:) = [   2]
  K1(11,:) = [  16]
  K1(12,:) = [   0]
  K1(13,:) = [   0]
  K1(14,:) = [ -18]
  K1(15,:) = [ -18]
  K1(16,:) = [  36]
  K1(17,:) = [   0]
  K1 = (1._dp / 3) * K1

  K2( 1,:) = [  12]
  K2( 2,:) = [   0]
  K2( 3,:) = [   0]
  K2( 4,:) = [   0]
  K2( 5,:) = [   0]
  K2( 6,:) = [   0]
  K2( 7,:) = [  16]
  K2( 8,:) = [   0]
  K2( 9,:) = [   0]
  K2(10,:) = [   2]
  K2(11,:) = [  16]
  K2(12,:) = [   0]
  K2(13,:) = [   0]
  K2(14,:) = [ -18]
  K2(15,:) = [ -18]
  K2(16,:) = [  36]
  K2(17,:) = [   0]
  K2 = (1._dp / 3) * K2

  KL( 1,:) = [  12]
  KL( 2,:) = [   0]
  KL( 3,:) = [   0]
  KL( 4,:) = [   0]
  KL( 5,:) = [   0]
  KL( 6,:) = [   0]
  KL( 7,:) = [  16]
  KL( 8,:) = [   0]
  KL( 9,:) = [   0]
  KL(10,:) = [   2]
  KL(11,:) = [  16]
  KL(12,:) = [   0]
  KL(13,:) = [   0]
  KL(14,:) = [ -18]
  KL(15,:) = [ -18]
  KL(16,:) = [  36]
  KL(17,:) = [   0]
  KL = (1._dp / 3) * KL


    co = 0
    do k = 1, 5
      co = co + k
      Cas(k) = K1(1+1*co,1)/K1(1,1)
    end do

  end subroutine colourmatrix_init
end module ol_colourmatrix_ppllj_eexddxg_2_dp



module ol_forced_parameters_ppllj_eexddxg_2_dp
  implicit none
  contains
  subroutine check_forced_parameters
    use ol_parameters_decl_dp
    use ol_loop_parameters_decl_dp

    implicit none
    logical, save :: checks_not_written = .true.

    if (checks_not_written) then
    ! e.g.
    ! if (ME /= 0) write(*,101) 'ME = 0'
  if (ME /= 0) write(*,101) 'ME = 0'
  if (CKMORDER /= 0) write(*,101) 'CKMORDER = 0'
  if (nc /= 3) write(*,101) 'nc = 3'
  if (nf /= 6) write(*,101) 'nf = 6'
  if (MU /= 0) write(*,101) 'MU = 0'
  if (MD /= 0) write(*,101) 'MD = 0'
  if (MS /= 0) write(*,101) 'MS = 0'
  if (LeadingColour /= 0) write(*,101) 'LeadingColour = 0'


    checks_not_written = .false.
    end if

    101 format('[OpenLoops] === WARNING ===',/,'[OpenLoops] code was generated with ',A,/,'[OpenLoops] ===============')
  end subroutine check_forced_parameters
end module ol_forced_parameters_ppllj_eexddxg_2_dp


! **********************************************************************
module ol_loop_storage_ppllj_eexddxg_2_dp
! **********************************************************************
  use kind_types, only: dp, intkind1, intkind2
  use ol_debug, only: ol_msg
  use ol_data_types_dp, only: wfun, Hpolcont


  integer(intkind1), save :: ntryL = 1
  integer(intkind1), save :: p_switch = 0 ! switch for dp or qp. Used for memory allocation of the OL types
  ! the following are flags for memory allocation of the hol coefficients in dp or qp
  logical, save :: dp_not_alloc = .TRUE., qp_not_alloc = .TRUE.
  integer, save :: n_merge_steps  ! total number of merging steps
  integer, save :: n_merge_mism   ! number of merging mismatches
  integer, save :: merge_step = 1 ! current merging step
  integer(intkind2), save, allocatable :: merge_tables(:,:,:) ! merging tables
  integer(intkind2), save, allocatable :: merge_mism(:)       ! array of merging mismatches
  integer(intkind2), save, allocatable :: merge_hels(:)       ! array with relevant helicities for a merging step
  logical, save :: merge_tables_on = .false.
  integer(intkind2), parameter :: nheltot = 32 ! number of helicity configurations
  integer(intkind2), save :: nhel = 32 ! number of non-vanishing helicity configurations (adapted at runtime)
  integer(intkind2), save :: hel_states = 32 ! number of helicity configurations needed for mem-allocation
  integer(intkind2), save :: Hel(32) ! physical helicity states
  integer(intkind2), save :: nhflip = 32    ! relevant helicities for helicity-flipped interference
  integer(intkind2), save :: helflip(2,32) ! Table for the helicity-flipped interference
  integer(intkind2), save, allocatable :: hflip(:,:)
  integer(intkind2) :: tsb
  integer, save :: pi_flip_bak = -1, pj_flip_bak = -1


  type(Hpolcont), save, allocatable :: A(:,:), M1helarray(:,:), M1helarray_ct(:,:)
  complex(dp), save :: den(21)

  ! external wave functions ex1(h1),... for h<n> helicities
  type(wfun) :: ex1(2), ex2(2), ex3(2), ex4(2), ex5(2)

  ! wf<h>(h,n) n wave functions with h helicity configurations
  type(wfun) :: wf4(4,14), wf8(8,13), wf32(32,17)

  ! diagram prefactors
  integer,           save :: fac_status_loop1 = -1, fac_status_loop2 = -1
  complex(dp), save :: f(17), c(10)

  !Vector in helicity and colour space for Born-Loop interference
  type(Hpolcont),   save, allocatable :: M0_col1_helarray(:,:), M0M1_hel_cc(:,:,:)
  complex(dp), save :: M2ctcc(16)

end module ol_loop_storage_ppllj_eexddxg_2_dp






! **********************************************************************
module ol_loop_ppllj_eexddxg_2_dp
! **********************************************************************
  use kind_types, only: dp, dp, intkind1, intkind2
  use ol_data_types_dp, only: wfun, Hpolcont
  use ol_loop_storage_ppllj_eexddxg_2_dp
  implicit none

!*********************************************************************************
  contains

! **********************************************************************
subroutine fac_init_loop()
! Writes diagram prefactors to 'f', rsp. 'c'
! **********************************************************************
  use ol_parameters_decl_dp
  use ol_loop_parameters_decl_dp
  use ol_init, only: set_parameter, tree_parameters_flush, parameters_flush

  implicit none
  call set_parameter("ew_renorm", 0)
  if (parameters_status == 0) call tree_parameters_flush()
  if (loop_parameters_status == 0) call parameters_flush()
  fac_status_loop1 = parameters_status
  fac_status_loop2 = loop_parameters_status
  ! factors of the diagrams
    f( 1) = (CI*eQED**2*gQCD)/3._dp
    f( 2) = CI*eQED**2*gQCD
    f( 3) = (CI*countertermnorm*eQED**2*gQCD**3)/3._dp
    f( 4) = CI*countertermnorm*eQED**2*gQCD**3
    f( 5) = (CI*countertermnorm*ctGqq*eQED**2*gQCD**3)/3._dp
    f( 6) = CI*countertermnorm*ctGqq*eQED**2*gQCD**3
    f( 7) = (CI*countertermnorm*ctVqq*eQED**2*gQCD**3)/3._dp
    f( 8) = CI*countertermnorm*ctVqq*eQED**2*gQCD**3
    f( 9) = countertermnorm*ctZGG*eQED**2*gQCD**3
    f(10) = (CI*eQED**2*gQCD**3*integralnorm*SwB)/3._dp
    f(11) = CI*eQED**2*gQCD**3*integralnorm*SwB
    f(12) = (eQED**2*gQCD**3*integralnorm*SwB)/3._dp
    f(13) = eQED**2*gQCD**3*integralnorm*SwB
    f(14) = (eQED**2*gQCD**3*integralnorm*SwF)/3._dp
    f(15) = (2*eQED**2*gQCD**3*integralnorm*SwF)/3._dp
    f(16) = eQED**2*gQCD**3*integralnorm*SwF
    f(17) = 2*eQED**2*gQCD**3*integralnorm*SwF

  c = [ 9*CI*f(10), 9*CI*f(11), f(12), 8*f(12), f(13), 8*f(13), 3*f(14), 3*f(15), 3*f(16), 3*f(17) ]
  c = (1._dp / 6) * c
end subroutine fac_init_loop

subroutine denominators()
  use ol_parameters_decl_dp ! masses
  use ol_momenta_decl_dp, only: L
  implicit none
  ! propagators
  den(1) = 1 /((L(5,3))+L(6,3))
  den(2) = 1 /((L(5,20))+L(6,20))
  den(4) = 1 /((L(5,3) - MZ2)+L(6,3))
  den(6) = 1 /((L(5,24))+L(6,24))
  den(9) = 1 /((L(5,12))+L(6,12))
  den(11) = 1 /((L(5,11))+L(6,11))
  den(16) = 1 /((L(5,7))+L(6,7))

  ! denominators
  den(3) = den(1)*den(2)
  den(5) = den(2)*den(4)
  den(7) = den(1)*den(6)
  den(8) = den(4)*den(6)
  den(10) = den(4)*den(9)
  den(12) = den(1)*den(11)
  den(13) = den(2)*den(12)
  den(14) = den(4)*den(11)
  den(15) = den(2)*den(14)
  den(17) = den(1)*den(16)
  den(18) = den(6)*den(17)
  den(19) = den(4)*den(16)
  den(20) = den(6)*den(19)
  den(21) = den(1)*den(9)

end subroutine denominators


! **********************************************************************
subroutine init_merging_tables(tot_num_hels, tot_num_merge_steps)
!-----------------------------------------------------------------------
! In the first event the merging tables are fully initialised
! **********************************************************************
  implicit none
  integer, intent(in) :: tot_num_hels, tot_num_merge_steps
  n_merge_steps = tot_num_merge_steps
  allocate(merge_tables(tot_num_hels,2,n_merge_steps))
  allocate(merge_mism(n_merge_steps+1))
  allocate(merge_hels(n_merge_steps))
  merge_tables_on = .true.
  merge_mism = 0_intkind2
  merge_tables = -1_intkind2
  merge_hels = -1_intkind2
end subroutine init_merging_tables

! **********************************************************************
subroutine update_merging_tables(tot_num_hels)
!-----------------------------------------------------------------------
! After the first evaluation the merging tables are adapted
! **********************************************************************
  implicit none
  integer, intent(in) :: tot_num_hels
  integer(intkind2), allocatable ::  merge_tables_tmp(:,:,:), merge_hels_tmp(:)

  if(merge_mism(1) == 0) then ! helicity configurations match in all merging steps
    if(allocated(merge_tables)) deallocate(merge_tables)
    if(allocated(merge_hels)) deallocate(merge_hels)
    if(allocated(merge_mism)) deallocate(merge_mism)
    allocate(merge_mism(1))
    allocate(merge_hels(1))
    allocate(merge_tables(1,1,1)) ! dummy allocation
    merge_mism = 0_intkind2
    merge_tables = -1_intkind2
    merge_hels = -1_intkind2
  else
    n_merge_mism = merge_mism(1)
    allocate(merge_tables_tmp(tot_num_hels,2,n_merge_mism))
    merge_tables_tmp(:,:,1:n_merge_mism) = merge_tables(:,:,1:n_merge_mism)
    if(allocated(merge_tables)) then
      deallocate(merge_tables)
      allocate(merge_tables(tot_num_hels,2,n_merge_mism))
    end if
    merge_tables = merge_tables_tmp
    if(allocated(merge_tables_tmp)) deallocate(merge_tables_tmp)
    allocate(merge_hels_tmp(n_merge_mism))
    merge_hels_tmp(1:n_merge_mism) = merge_hels(1:n_merge_mism)
    if(allocated(merge_hels)) then
      deallocate(merge_hels)
      allocate(merge_hels(n_merge_mism))
    end if
    merge_hels = merge_hels_tmp
    if(allocated(merge_hels_tmp)) deallocate(merge_hels_tmp)
  end if
  merge_tables_on = .true.
end subroutine update_merging_tables


! **********************************************************************
subroutine allocate_diagrams()
!-----------------------------------------------------------------------
! After the first Born evaluation, colour-stripped amplitudes and
! colour vectors are initialiased with the minimum number of relevant
! helicity states.
! **********************************************************************

  implicit none

  if (allocated(A)) deallocate(A)
  if (allocated(M1helarray)) deallocate(M1helarray)
  if (allocated(M1helarray_ct)) deallocate(M1helarray_ct)
  if (allocated(M0_col1_helarray)) deallocate(M0_col1_helarray)
  allocate(A(hel_states,17))
  allocate(M1helarray(1,hel_states))
  allocate(M1helarray_ct(1,hel_states))
  allocate(M0_col1_helarray(1,hel_states))
  A(1:hel_states,1:17)%j = 0
  A(1:hel_states,1:17)%e = 0
  A(1:hel_states,1:17)%hf = 0
  A(1:hel_states,1:17)%s = 0

  if (allocated(M0M1_hel_cc)) deallocate(M0M1_hel_cc)
  allocate(M0M1_hel_cc(1,hel_states,16))

end subroutine allocate_diagrams

!
!
!
! **********************************************************************

recursive subroutine amp2(P_scatt, M02, qp_kinematics, M2ct, M2colint)

! P_scatt(0:3,Npart) = incoming external momenta
! M2  = helicity-summed squared matrix element for e- e+ down anti-down glue -> 0
! **********************************************************************
  use ol_parameters_decl_dp !, only: ci, parameters_status, ZERO, scalefactor, >masses<
  use ol_parameters_init_dp, only: ensure_mp_init, ensure_mp_loop_init
  use ol_kinematics_dp, only: init_kinematics
  use ol_momenta_decl_dp, only: momenta_nan_check
  use ol_settings_ppllj_eexddxg_2, only: hel_mem_opt, loopcc
  use ol_data_types_dp
  use ol_h_helicity_bookkeeping_dp, only: &
    & helbookkeeping_wf, helsync, flip_phase
  use ol_helicity_init, only: helbookkeeping_flip, helsync_flip
  use ol_hel_propagators_dp
  use ol_hel_wavefunctions_dp
  use ol_wavefunctions_dp, only: wf_V_Std
! 
  use ol_hel_vertices_dp
  use ol_hel_contractions_dp
  use ol_external_ppllj_eexddxg_2, only: &
    & external_perm_ppllj_eexddxg_2, &
    & external_perm_inv_ppllj_eexddxg_2, &
    & extcomb_perm_ppllj_eexddxg_2, &
    & average_factor_ppllj_eexddxg_2
  use ol_external_ppllj_eexddxg_2, only: &
    & H, hel_not_initialised, hel_init, POLSEL 
  use ol_colourmatrix_ppllj_eexddxg_2_dp, only: &
    & colmat_not_initialised, colourmatrix_init
  use ol_forced_parameters_ppllj_eexddxg_2_dp, only: &
    & check_forced_parameters
  use ol_heltables_OLR_ppllj_eexddxg_2
  use ol_kinematics_dp, only: LC2Std_Rep_cmplx
  use ol_h_counterterms_dp
  use ol_loop_parameters_decl_dp ! counterterms
  use ol_loop_parameters_decl_dp, only: &
    & IR_is_on, DOI, CT_is_on, R2_is_on, TP_is_on
  use ol_init, only: set_parameter, parameters_flush

  use ol_kinematics_qp, only: conv_mom_scatt2in_qp=>conv_mom_scatt2in, &
                                         internal_momenta_qp=>internal_momenta
  use ol_momenta_decl_dp, only: L
  use ol_momenta_decl_qp, only: L_qp=>L

  implicit none

  real(dp), intent(in)  :: P_scatt(0:3,5)
  real(dp),  intent(out) :: M02
  logical, intent(in) :: qp_kinematics
  real(dp),  intent(out), optional :: M2ct
  real(dp),  intent(out), optional :: M2colint(16)
  real(dp) :: iM2ct
  real(dp) :: iM2colint(16)

  integer(intkind1), save :: ntry = 0

  integer           :: shift, k, r, m, n, i
  real(dp)    :: P(0:3,5)

  real(qp)    :: P_qp(0:3,5)

  integer           :: extmasses2(5)
  real(dp)    :: M2add, M2add_ct, M2add_colint(16)
  complex(dp) :: M1(1), M2(1)
  real(dp)    :: P_scatt_intern(0:3,5)
  real(dp), save :: scalebackfactor, old_scalefactor = 0
  integer(intkind1) :: nsync
  integer, allocatable :: extcombs_permuted(:)
  integer              :: extcombs(16), nextcombs
  integer              :: CT_on_bak, R2_on_bak, TP_on_bak, DOI_bak
  logical              :: do_ct, do_colint

  integer :: ind_cc_comb(5)

  complex(dp) :: omega(2) ! phases for helicity correlations

  if (present(M2ct)) then
    do_ct = .true.
  else
    do_ct = .false.
  end if

  if (present(M2colint)) then
    do_colint = .true.
  else
    do_colint = .false.
  end if

  if(ntry == 0) then
    ! recursive initialization call needed for the correct helicity bokkeeping
    ntry = 1
    CT_on_bak = CT_is_on
    R2_on_bak = R2_is_on
    TP_on_bak = TP_is_on
    DOI_bak = DOI
    call set_parameter("ct_on", 1)
    call set_parameter("r2_on", 1)
    call set_parameter("tp_on", 1)
    DOI = 1
    call amp2(P_scatt, M02, qp_kinematics, M2ct=iM2ct, M2colint=iM2colint)
    call set_parameter("ct_on", CT_on_bak)
    call set_parameter("r2_on", R2_on_bak)
    call set_parameter("tp_on", TP_on_bak)
    DOI = DOI_bak
  end if

  if (ntry < 2) then
    if (allocated(A)) deallocate(A)
    allocate(A(nhel,17))
  end if

  if (do_ct) call set_parameter("ew_renorm", 0)
  call parameters_flush()
  call ensure_mp_init()
  if (do_ct) call ensure_mp_loop_init()

  if (colmat_not_initialised) call colourmatrix_init()

  if (fac_status_loop1 /= parameters_status .or. fac_status_loop2 /= loop_parameters_status) then
    call check_forced_parameters()
  end if
  if (do_ct) call fac_init_loop()

  if (momenta_nan_check(P_scatt) /= 0) then
    M02 = 0
    return
  end if

  extmasses2 = [ 0, 0, 0, 0, 0 ]
  ! Convert 2 -> n-2 PS-point to n -> 0 (so that P(1) + ... + P(n) = 0) and compute
  ! internal-propagator momenta in light-cone representation
  call init_kinematics(P_scatt, extmasses2, P, &
     external_perm_inv_ppllj_eexddxg_2, 5, qp_kinematics)

  ! denominators
  call denominators()

   if (heltables_not_init) call init_heltables()

  ! external WFs
  ! Here the external wavefunctions are initialiased
  call pol_wf_Q(P(:,1), rZERO, H1, ex1, POLSEL(1),1)
  call pol_wf_A(P(:,2), rZERO, H2, ex2, POLSEL(2),2)
  call pol_wf_Q(P(:,3), rZERO, H3, ex3, POLSEL(3),3)
  call pol_wf_A(P(:,4), rZERO, H4, ex4, POLSEL(4),4)
  call pol_wf_V(P(:,5), rZERO, H5, ex5, POLSEL(5),5)



  ! internal WFs
  ! e.g. call vert_VQ_A(ntry, ex3, ex1, wf1, n1, t1) ...
  call vert_QA_V(ntry, ex1(:), ex2(:), wf4(:,1), n3(:,1), t3x4(:,:,1))
  call vert_VQ_A(ntry, ex5(:), ex3(:), wf4(:,2), n3(:,2), t3x4(:,:,2))
  call prop_Q_A(ntry, wf4(:,2), 20, ZERO, 0_intkind1, wf4(:,3), n2(1))
  call vert_AV_Q(ntry, ex4(:), wf4(:,1), wf8(:,1), n3(:,3), t3x8(:,:,1))
  call vert_QA_Z(gZl,ntry, ex1(:), ex2(:), wf4(:,4), n3(:,4), t3x4(:,:,3))
  call prop_W_W(ntry, wf4(:,4), 3, MZ, 1_intkind1, wf4(:,5), n2(2))
  call vert_AZ_Q(gZd,ntry, ex4(:), wf4(:,5), wf8(:,2), n3(:,5), t3x8(:,:,2))
  call vert_AV_Q(ntry, ex4(:), ex5(:), wf4(:,6), n3(:,6), t3x4(:,:,4))
  call prop_A_Q(ntry, wf4(:,6), 24, ZERO, 0_intkind1, wf4(:,7), n2(3))
  call vert_VQ_A(ntry, wf4(:,1), ex3(:), wf8(:,3), n3(:,7), t3x8(:,:,3))
  call vert_ZQ_A(gZd,ntry, wf4(:,5), ex3(:), wf8(:,4), n3(:,8), t3x8(:,:,4))
  call vert_QA_V(ntry, ex3(:), ex4(:), wf4(:,8), n3(:,9), t3x4(:,:,5))
  call counter_VG_G(ntry, wf4(:,5), ex5(:), 16, wf8(:,5), 19, n3(:,10), t3x8(:,:,5))
  call counter_AV_Q(ntry, ex4(:), wf4(:,1), wf8(:,6), n3(:,11), t3x8(:,:,6))
  call counter_AZ_Q(gZd,ntry, ex4(:), wf4(:,5), wf8(:,7), n3(:,12), t3x8(:,:,7))
  call counter_AV_Q(ntry, ex4(:), ex5(:), wf4(:,9), n3(:,13), t3x4(:,:,6))
  call prop_A_Q(ntry, wf4(:,9), 24, ZERO, 0_intkind1, wf4(:,10), n2(4))
  call counter_VQ_A(ntry, wf4(:,1), ex3(:), wf8(:,8), n3(:,14), t3x8(:,:,8))
  call counter_ZQ_A(gZd,ntry, wf4(:,5), ex3(:), wf8(:,9), n3(:,15), t3x8(:,:,9))
  call counter_VQ_A(ntry, ex5(:), ex3(:), wf4(:,11), n3(:,16), t3x4(:,:,7))
  call prop_Q_A(ntry, wf4(:,11), 20, ZERO, 0_intkind1, wf4(:,12), n2(5))
  call counter_Q_A(ctqq,1,ntry, wf4(:,3), 20, wf4(:,13), n2(6))
  call prop_A_Q(ntry, wf8(:,1), 11, ZERO, 0_intkind1, wf8(:,10), n2(7))
  call prop_A_Q(ntry, wf8(:,2), 11, ZERO, 0_intkind1, wf8(:,11), n2(8))
  call counter_A_Q(ctqq,1,ntry, wf4(:,7), 24, wf4(:,14), n2(9))
  call prop_Q_A(ntry, wf8(:,3), 7, ZERO, 0_intkind1, wf8(:,12), n2(10))
  call prop_Q_A(ntry, wf8(:,4), 7, ZERO, 0_intkind1, wf8(:,13), n2(11))


  if (ntry==1) then
    call physical_helicities()
  end if

  ! computation of the colour-stripped amplitudes
  do nsync = ntry+ntry-1, ntry+1  !  nsync = 1,2  for 1st point and nsync = 3 later
    call diagrams()
    if (nsync == 1) call helsync(nsync, A, nhel, Hel)
  end do

  if (ntry < 2) then
    if(hel_mem_opt) hel_states = nhel
    call allocate_diagrams()
  end if

  ! In the following loop the coefficients \Gamma_{i} of the expansion in
  ! the colour basis are computed for every helicity state and stored in M1helarray.
  ! Also the colour vector for the Born-loop interference is computed and saved in
  ! M0_col1_helarray
  do k = 1, nhel
    if (do_ct) call colourvectors(A, k, M1helarray(:,k),M1helarray_ct(:,k))
    call colborninterf(M1helarray(:,k), M0_col1_helarray(:,k), 0)
  end do
  M1helarray(:,nhel+1:)%j = 0
  M0_col1_helarray(:,nhel+1:)%hf = -1_intkind2
  M0_col1_helarray(:,nhel+1:)%j = 0

  M2add = 0
  M2add_ct = 0
  M2add_colint = 0

  M02 = 0
  if (do_ct) M2ct = 0
  if (do_colint) M2colint = 0

  do k = 1, nhel
    call colint(M1helarray(:,k)%j, M2add)
    if (do_ct) call colint_ct(M1helarray(:,k)%j, M1helarray_ct(:,k)%j, M2add_ct, 0)
    if (IR_is_on > 0) then
      if (do_colint) call colint_IR(M1helarray(:,k)%j, M2add_colint)
    end if
    !summation over helicity configurations
    M02 = M02 + M2add
    if (do_ct) M2ct = M2ct + M2add_ct
    if (do_colint) M2colint = M2colint + M2add_colint
  end do


  if (loopcc) then
    !!Born-loop colour correlators interference
    ind_cc_comb = [2,4,7,5,8]
    M0M1_hel_cc(:,:,:)%j = 0
    M2ctcc = 0
    do n = 1, 5
      M2add_ct = 0
      m = ind_cc_comb(n)
        do k = 1, nhel
          call colint_ct(M1helarray(:,k)%j, M1helarray_ct(:,k)%j, M2add_ct, m)
          call colborninterf(M1helarray(:,k), M0M1_hel_cc(:,k,m), m)
          M2ctcc(m) = M2ctcc(m) + M2add_ct
        end do
    end do
    M0M1_hel_cc(:,nhel+1:,:)%hf = -1_intkind2
    M0M1_hel_cc(:,nhel+1:,:)%j = 0
  end if


  if (ntry < 1) then
    if (allocated(A)) deallocate(A)
  end if
  ntry=2

  contains

subroutine physical_helicities()
  implicit none

integer :: i1,i2,i3,i4,i5
integer :: n

n=0

do i5= 1, 2
do i4= 1, 2
do i3= 1, 2
do i2= 1, 2
do i1= 1, 2
n = n + 1
 if((ex1(i1)%hf==-1_intkind2) .OR. (ex2(i2)%hf==-1_intkind2) .OR. (ex3(i3)%hf==-1_intkind2) .OR. (ex4(i4)%hf==-1_intkind2)  &
    .OR. (ex5(i5)%hf==-1_intkind2)) then
  Hel(n) = -1_intkind2
 else
  Hel(n) = ex1(i1)%hf + ex2(i2)%hf + ex3(i3)%hf + ex4(i4)%hf + ex5(i5)%hf
 end if
end do
end do
end do
end do
end do

end subroutine physical_helicities

subroutine diagrams()
  implicit none
  integer :: h
  ! e.g. call cont_VV(nsync, wf3, wf6, A(:,1), n64, t64, nhel, den(5)) ...

    call Hcont_QA(nsync, wf4(:,3), wf8(:,1), A(:,1), n3(:,17), t3x32(:,:,1), nhel, den(3))
    call Hcont_QA(nsync, wf4(:,3), wf8(:,2), A(:,2), n3(:,18), t3x32(:,:,2), nhel, den(5))
    call Hcont_QA(nsync, wf4(:,7), wf8(:,3), A(:,3), n3(:,19), t3x32(:,:,3), nhel, den(7))
    call Hcont_QA(nsync, wf4(:,7), wf8(:,4), A(:,4), n3(:,20), t3x32(:,:,4), nhel, den(8))

    call Hcont_VV(nsync, wf4(:,8), wf8(:,5), A(:,5), n3(:,21), t3x32(:,:,5), nhel, den(10))
    call Hcont_QA(nsync, wf4(:,3), wf8(:,6), A(:,6), n3(:,22), t3x32(:,:,6), nhel, den(3))
    call Hcont_QA(nsync, wf4(:,3), wf8(:,7), A(:,7), n3(:,23), t3x32(:,:,7), nhel, den(5))
    call Hcont_QA(nsync, wf8(:,3), wf4(:,10), A(:,8), n3(:,24), t3x32(:,:,8), nhel, den(7))
    call Hcont_QA(nsync, wf8(:,4), wf4(:,10), A(:,9), n3(:,25), t3x32(:,:,9), nhel, den(8))
    call Hcont_QA(nsync, wf4(:,7), wf8(:,8), A(:,10), n3(:,26), t3x32(:,:,10), nhel, den(7))
    call Hcont_QA(nsync, wf4(:,7), wf8(:,9), A(:,11), n3(:,27), t3x32(:,:,11), nhel, den(8))
    call Hcont_QA(nsync, wf8(:,1), wf4(:,12), A(:,12), n3(:,28), t3x32(:,:,12), nhel, den(3))
    call Hcont_QA(nsync, wf8(:,2), wf4(:,12), A(:,13), n3(:,29), t3x32(:,:,13), nhel, den(5))
    call Hcont_QA(nsync, wf4(:,13), wf8(:,10), A(:,14), n3(:,30), t3x32(:,:,14), nhel, den(13))
    call Hcont_QA(nsync, wf4(:,13), wf8(:,11), A(:,15), n3(:,31), t3x32(:,:,15), nhel, den(15))
    call Hcont_QA(nsync, wf4(:,14), wf8(:,12), A(:,16), n3(:,32), t3x32(:,:,16), nhel, den(18))
    call Hcont_QA(nsync, wf4(:,14), wf8(:,13), A(:,17), n3(:,33), t3x32(:,:,17), nhel, den(20))

end subroutine diagrams


elemental function diagmap(j, n)
  implicit none
  integer, intent(in) :: j, n
  complex(dp) :: diagmap
  diagmap = A(j,n)%j
end function diagmap

function diagsum(j, pos, neg)
  implicit none
  integer, intent(in) :: j, pos(:), neg(:)
  complex(dp) :: diagsum
  diagsum = sum(diagmap(j, pos)) - sum(diagmap(j, neg))
end function diagsum

subroutine colourvectors(A, j, M1, M2)
  implicit none
  type(Hpolcont) :: A(:,:)
  integer, intent(in) :: j
  type(Hpolcont), intent(out) :: M1(1), M2(1) ! M1helarray(1,nhel)
  integer :: empty(0), i

  M1(1)%j = ((A(j,1)%j+A(j,3)%j)*f(1)+(A(j,2)%j+A(j,4)%j)*f(2))

  M2(1)%j = ((-A(j,14)%j-A(j,16)%j)*f(3)+(-A(j,15)%j-A(j,17)%j)*f(4)+(A(j,8)%j+A(j,12)%j)*f(5)+(A(j,9)%j+A(j,13)%j)*f(6)+(A(j,6)%j &
       +A(j,10)%j)*f(7)+(A(j,7)%j+A(j,11)%j)*f(8)-A(j,5)%j*f(9))


  M1(:)%hf = Hel(j)
  M2(:)%hf = Hel(j)

  !M(i) corresponds to \Gamma_{i} in Fabios thesis

end subroutine colourvectors


! **********************************************************************
subroutine colint(M, M2colint)
! M(i)   = <M|Ci> colour component of matrix element
! COLINT = <M|M>
!        = Sum_{i,j} <M|Ci> * <Ci|Cj> * <Cj|M>
!        = colour-summed squared matrix element
! K1(i,j) = <Ci|Cj>
! **********************************************************************
  use ol_colourmatrix_ppllj_eexddxg_2_dp, only: K1
  implicit none

  complex(dp), intent(in)  :: M(1)
  real(dp),    intent(out) :: M2colint
  integer :: i, j

  M2colint = 0

    do i = 1, 1
      do j = 1, 1
        M2colint = M2colint + real(conjg(M(i))*K1(i,j)*M(j))
      end do
    end do

end subroutine colint

! **********************************************************************
subroutine colint_ct(M, Mct, M2colint_ct, l)
! M(i)   = <M|Ci> colour component of matrix element
! COLINT = <M|M>
!        = Sum_{i,j} <M|Ci> * <Ci|Cj> * <Cj|M>
!        = colour-summed squared matrix element
! K2(i,j) = <Ci|Cj>
! **********************************************************************
  use ol_colourmatrix_ppllj_eexddxg_2_dp, only: K2
  implicit none
  complex(dp), intent(in)  :: M(1), Mct(1)
  real(dp),    intent(out) :: M2colint_ct
  integer, intent(in) :: l
  integer ::  i, j

  M2colint_ct = 0

  do i = 1, 1
    do j = 1, 1
      M2colint_ct = M2colint_ct + real(conjg(M(i))*K2(i+1*l,j)*Mct(j))
    end do
  end do

end subroutine colint_ct

! **********************************************************************
subroutine colint_IR(M, M2IRadd)
! M(i)   = <M|Ci> colour component of matrix element
! K1(i,j) = .....
! **********************************************************************
  use ol_colourmatrix_ppllj_eexddxg_2_dp, only: K1
  implicit none

  complex(dp), intent(in)  :: M(1)
  real(dp),    intent(out) :: M2IRadd(16)
  integer ::  i, j, k, colmatpos

  M2IRadd = 0

  do k = 1, 17-1
    colmatpos = 1*k
    do i = 1, 1
      do j = 1, 1
        M2IRadd(k) = M2IRadd(k) + real(conjg(M(i))*K1(i+colmatpos,j)*M(j))
      end do
    end do
  end do

end subroutine colint_IR

! **********************************************************************
subroutine colintmunu(M1, M2, M2colint)
! M1(i)    = <M1|Ci> colour component of matrix element
! M2(i)    = <M2|Ci> colour component of matrix element
! M2colint = <M1|M2>
!          = Sum_{i,j} <M1|Ci> * <Ci|Cj> * <Cj|M2>
!          = colour-summed squared matrix element
! K2(i,j) = <Ci|Cj>
! **********************************************************************
  use ol_colourmatrix_ppllj_eexddxg_2_dp, only: K2
  implicit none

  complex(dp), intent(in)  :: M1(1)
  complex(dp), intent(in)  :: M2(1)
  complex(dp),    intent(out) :: M2colint
  integer :: i, j

  M2colint = 0

  do i = 1, 1
    do j = 1, 1
      M2colint = M2colint + M1(i)*K2(i,j)*conjg(M2(j))
    end do
  end do

end subroutine colintmunu


! **********************************************************************
subroutine colborninterf(M, M0_col, l)
! M(i)         = <M|Ci> colour component of matrix element
! M0_col(i)    = <M2|Ci> colour component of matrix element,
!                see \tilde{M}_{j} in Fabios thesis
! M2colint = Sum_{i} <M1|Ci> * <Ci|Cj>
!          = colour-summed squared matrix element
! KL(i,j) = <Ci|Cj> with elements Cj of 1-loop colour basis
! **********************************************************************
  use ol_colourmatrix_ppllj_eexddxg_2_dp, only: KL
  implicit none

  type(Hpolcont), intent(in)  :: M(1)
  type(Hpolcont), intent(out)  :: M0_col(1)
  integer, intent(in) :: l
  integer :: i, j

  do j = 1, 1 !size(KL(1,:))
    M0_col(j)%j = 0
    M0_col(j)%hf = M(1)%hf
    do i = 1, 1 !size(KL(:,1))
      M0_col(j)%j = M0_col(j)%j  + conjg(M(i)%j)*KL(i+1*l,j)
    end do
  end do

end subroutine colborninterf

end subroutine amp2


! **********************************************************************

subroutine hsfamp2(P_scatt, M02, pind_1, pind_2, M2hsf)

! M2tree = helicity-summed squared tree matrix element for e- e+ down anti-down glue -> 0
! M2hsf  = formula (66) of 1011.3918. Tree level
! spin-correlated hard scattering factor for e- e+ down anti-down glue -> 0
! **********************************************************************
  use ol_parameters_decl_dp
  implicit none
  real(dp), intent(in) :: P_scatt(0:3,5)
  integer, intent(in) :: pind_1, pind_2
  real(dp), intent(out) :: M02, M2hsf
  real(dp)    :: M2add, M2ct_dummy
  integer :: k

  M2add = 0
  M2hsf = 0

  call amp2(P_scatt, M02, .false., M2ct=M2ct_dummy)
  if((pi_flip_bak .ne. pind_1) .or. (pj_flip_bak .ne. pind_2)) then
    call init_flipped_heltables(pind_1,pind_2)
    pi_flip_bak = pind_1
    pj_flip_bak = pind_2
  end if

  do k = 1, nhflip
    call colintHSF(M1helarray(:,helflip(1,k))%j, M1helarray(:,helflip(2,k))%j, M2add)
    M2hsf = M2hsf + M2add
  end do

  contains

subroutine assign_helicity_flip(ext_wf_hf,ext_wf_t,m)
  use ol_debug, only: ol_fatal
  use kind_types, only: intkind2
  implicit none
  integer(intkind2), intent(in) :: ext_wf_hf(:), ext_wf_t
  integer, intent(in) :: m
  if(.not. allocated(hflip)) allocate(hflip(size(ext_wf_hf),2))
  if(size(ext_wf_hf) == size(hflip,1)) then
    hflip(:,m) = ext_wf_hf
    tsb = tsb + ext_wf_t
  else
    call ol_fatal("Spin-correlated HSF for ppllj_eexddxg_2:" &
      & // "trying to flip particles with different number of helicity states")
  end if
end subroutine assign_helicity_flip

subroutine init_flipped_heltables(i,j)
  use kind_types, only: intkind2
  use ol_external_ppllj_eexddxg_2, only: &
    & external_perm_inv_ppllj_eexddxg_2
  use ol_h_helicity_bookkeeping_dp, only: &
    & helicity_flip_ij
  implicit none
  integer, intent(in) :: i, j
  integer :: r, k, l=0
  tsb = 0

  if(allocated(hflip)) deallocate(hflip)

  do k = 1, 5
    r = external_perm_inv_ppllj_eexddxg_2(k)
    if(r == i .or. r == j) then
      l = l+1
      select case(k)
        case(1)
          call assign_helicity_flip(ex1(:)%hf,ex1(1)%t,l)
        case(2)
          call assign_helicity_flip(ex2(:)%hf,ex2(1)%t,l)
        case(3)
          call assign_helicity_flip(ex3(:)%hf,ex3(1)%t,l)
        case(4)
          call assign_helicity_flip(ex4(:)%hf,ex4(1)%t,l)
        case(5)
          call assign_helicity_flip(ex5(:)%hf,ex5(1)%t,l)
      end select
    end if
  end do
  call helicity_flip_ij(32,nhel, M1helarray(1,:)%hf, hflip, tsb, &
      & nhflip, helflip)

end subroutine init_flipped_heltables

! **********************************************************************
subroutine colintHSF(M, Mflip, Mhfcolint)
! **********************************************************************
  use ol_colourmatrix_ppllj_eexddxg_2_dp, only: K1
  implicit none
  complex(dp), intent(in)  :: M(1), Mflip(1)
  real(dp),    intent(out) :: Mhfcolint
  integer :: i, j
  MHfcolint = 0
  do i = 1, 1
    do j = 1, 1
      Mhfcolint = Mhfcolint + real(conjg(Mflip(i))*K1(i,j)*M(j))
    end do
  end do

end subroutine colintHSF

end subroutine hsfamp2

end module ol_loop_ppllj_eexddxg_2_dp

