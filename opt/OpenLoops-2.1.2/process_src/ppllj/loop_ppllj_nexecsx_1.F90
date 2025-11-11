
module ol_colourmatrix_ppllj_nexecsx_1_/**/REALKIND
  use KIND_TYPES, only: REALKIND
  implicit none
  logical, save           :: colmat_not_initialised = .true.
  complex(REALKIND), save :: K1(12,1), K2(12,1), KL(12,1), Cas(4) = 0
  contains
  subroutine colourmatrix_init
    use ol_parameters_decl_/**/REALKIND, only: CI
    implicit none
    integer :: k, co
    colmat_not_initialised = .false.
    ! colour matrix

  K1( 1,:) = [  3]
  K1( 2,:) = [  0]
  K1( 3,:) = [  0]
  K1( 4,:) = [  0]
  K1( 5,:) = [  0]
  K1( 6,:) = [  0]
  K1( 7,:) = [  4]
  K1( 8,:) = [  0]
  K1( 9,:) = [  0]
  K1(10,:) = [ -4]
  K1(11,:) = [  4]
  K1(12,:) = [  0]

  K2( 1,:) = [  3]
  K2( 2,:) = [  0]
  K2( 3,:) = [  0]
  K2( 4,:) = [  0]
  K2( 5,:) = [  0]
  K2( 6,:) = [  0]
  K2( 7,:) = [  4]
  K2( 8,:) = [  0]
  K2( 9,:) = [  0]
  K2(10,:) = [ -4]
  K2(11,:) = [  4]
  K2(12,:) = [  0]

  KL( 1,:) = [  3]
  KL( 2,:) = [  0]
  KL( 3,:) = [  0]
  KL( 4,:) = [  0]
  KL( 5,:) = [  0]
  KL( 6,:) = [  0]
  KL( 7,:) = [  4]
  KL( 8,:) = [  0]
  KL( 9,:) = [  0]
  KL(10,:) = [ -4]
  KL(11,:) = [  4]
  KL(12,:) = [  0]

#if 1 > 0
    co = 0
    do k = 1, 4
      co = co + k
      Cas(k) = K1(1+1*co,1)/K1(1,1)
    end do
#endif
  end subroutine colourmatrix_init
end module ol_colourmatrix_ppllj_nexecsx_1_/**/REALKIND



module ol_forced_parameters_ppllj_nexecsx_1_/**/REALKIND
  implicit none
  contains
  subroutine check_forced_parameters
    use ol_parameters_decl_/**/REALKIND
    use ol_loop_parameters_decl_/**/REALKIND
#ifndef PRECISION_dp
    use ol_loop_parameters_decl_/**/DREALKIND, only: LeadingColour, nc, nf, CKMORDER
#endif
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
  if (wMC /= 0) write(*,101) 'wMC = 0'


    checks_not_written = .false.
    end if

    101 format('[OpenLoops] === WARNING ===',/,'[OpenLoops] code was generated with ',A,/,'[OpenLoops] ===============')
  end subroutine check_forced_parameters
end module ol_forced_parameters_ppllj_nexecsx_1_/**/REALKIND


! **********************************************************************
module ol_loop_storage_ppllj_nexecsx_1_/**/REALKIND
! **********************************************************************
  use KIND_TYPES, only: REALKIND, intkind1, intkind2
  use ol_debug, only: ol_msg
  use ol_data_types_/**/REALKIND, only: wfun, Hpolcont

#ifdef PRECISION_dp
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
  integer(intkind2), parameter :: nheltot = 16 ! number of helicity configurations
  integer(intkind2), save :: nhel = 16 ! number of non-vanishing helicity configurations (adapted at runtime)
  integer(intkind2), save :: hel_states = 16 ! number of helicity configurations needed for mem-allocation
  integer(intkind2), save :: Hel(16) ! physical helicity states
  integer(intkind2), save :: nhflip = 16    ! relevant helicities for helicity-flipped interference
  integer(intkind2), save :: helflip(2,16) ! Table for the helicity-flipped interference
  integer(intkind2), save, allocatable :: hflip(:,:)
  integer(intkind2) :: tsb
  integer, save :: pi_flip_bak = -1, pj_flip_bak = -1
#endif

  type(Hpolcont), save, allocatable :: A(:,:), M1helarray(:,:), M1helarray_ct(:,:)
  complex(REALKIND), save :: den(1)

  ! external wave functions ex1(h1),... for h<n> helicities
  type(wfun) :: ex1(2), ex2(2), ex3(2), ex4(2)

  ! wf<h>(h,n) n wave functions with h helicity configurations
  type(wfun) :: wf4(4,4), wf16(16,2)

  ! diagram prefactors
  integer,           save :: fac_status_loop1 = -1, fac_status_loop2 = -1
  complex(REALKIND), save :: f(3), c(1)

  !Vector in helicity and colour space for Born-Loop interference
  type(Hpolcont),   save, allocatable :: M0_col1_helarray(:,:), M0M1_hel_cc(:,:,:)
  complex(REALKIND), save :: M2ctcc(11)

end module ol_loop_storage_ppllj_nexecsx_1_/**/REALKIND






! **********************************************************************
module ol_loop_ppllj_nexecsx_1_/**/REALKIND
! **********************************************************************
  use KIND_TYPES, only: REALKIND, DREALKIND, intkind1, intkind2
  use ol_data_types_/**/REALKIND, only: wfun, Hpolcont
  use ol_loop_storage_ppllj_nexecsx_1_/**/REALKIND
  implicit none

!*********************************************************************************
  contains

! **********************************************************************
subroutine fac_init_loop()
! Writes diagram prefactors to 'f', rsp. 'c'
! **********************************************************************
  use ol_parameters_decl_/**/REALKIND
  use ol_loop_parameters_decl_/**/REALKIND
  use ol_init, only: set_parameter, tree_parameters_flush, parameters_flush
#ifndef PRECISION_dp
  use ol_loop_parameters_decl_/**/DREALKIND, only: SwF, SwB
!  use ol_loop_parameters_decl_/**/DREALKIND, only: DOI
#endif
  implicit none
  call set_parameter("ew_renorm", 0)
  if (parameters_status == 0) call tree_parameters_flush()
  if (loop_parameters_status == 0) call parameters_flush()
  fac_status_loop1 = parameters_status
  fac_status_loop2 = loop_parameters_status
  ! factors of the diagrams
    f(1) = (CI*eQED**2)/(2._/**/REALKIND*sw**2)
    f(2) = (CI*countertermnorm*ctVsc*eQED**2*gQCD**2)/(2._/**/REALKIND*sw**2)
    f(3) = (eQED**2*gQCD**2*integralnorm*SwB)/(sw**2*2._/**/REALKIND)

  c = [ 4*f(3) ]
  c = (1._/**/REALKIND / 3) * c
end subroutine fac_init_loop

subroutine denominators()
  use ol_parameters_decl_/**/REALKIND ! masses
  use ol_momenta_decl_/**/REALKIND, only: L
  implicit none
  ! propagators
  den(1) = 1 /((L(5,3) - MW2)+L(6,3))

  ! denominators

end subroutine denominators

#ifdef PRECISION_dp
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
#endif

! **********************************************************************
subroutine allocate_diagrams()
!-----------------------------------------------------------------------
! After the first Born evaluation, colour-stripped amplitudes and
! colour vectors are initialiased with the minimum number of relevant
! helicity states.
! **********************************************************************
#ifndef PRECISION_dp
  use ol_loop_storage_ppllj_nexecsx_1_/**/DREALKIND, only: hel_states
#endif
  implicit none

  if (allocated(A)) deallocate(A)
  if (allocated(M1helarray)) deallocate(M1helarray)
  if (allocated(M1helarray_ct)) deallocate(M1helarray_ct)
  if (allocated(M0_col1_helarray)) deallocate(M0_col1_helarray)
  allocate(A(hel_states,2))
  allocate(M1helarray(1,hel_states))
  allocate(M1helarray_ct(1,hel_states))
  allocate(M0_col1_helarray(1,hel_states))
  A(1:hel_states,1:2)%j = 0
  A(1:hel_states,1:2)%e = 0
  A(1:hel_states,1:2)%hf = 0
  A(1:hel_states,1:2)%s = 0
#if 1 > 0
  if (allocated(M0M1_hel_cc)) deallocate(M0M1_hel_cc)
  allocate(M0M1_hel_cc(1,hel_states,11))
#endif
end subroutine allocate_diagrams

!
!
!
! **********************************************************************
#ifdef PRECISION_dp
recursive subroutine amp2(P_scatt, M02, qp_kinematics, M2ct, M2colint)
#else
recursive subroutine amp2(P_scatt, M02, qp_kinematics, M2ct, M2colint)
  use ol_loop_storage_ppllj_nexecsx_1_/**/DREALKIND, only: &
  nhel, Hel, hel_states
#endif
! P_scatt(0:3,Npart) = incoming external momenta
! M2  = helicity-summed squared matrix element for anti-nu_e e- charm anti-strange -> 0
! **********************************************************************
  use ol_parameters_decl_/**/REALKIND !, only: ci, parameters_status, ZERO, scalefactor, >masses<
  use ol_parameters_init_/**/REALKIND, only: ensure_mp_init, ensure_mp_loop_init
  use ol_kinematics_/**/REALKIND, only: init_kinematics
  use ol_momenta_decl_/**/DREALKIND, only: momenta_nan_check
  use ol_settings_ppllj_nexecsx_1, only: hel_mem_opt, loopcc
  use ol_data_types_/**/REALKIND
  use ol_h_helicity_bookkeeping_/**/REALKIND, only: &
    & helbookkeeping_wf, helsync, flip_phase
  use ol_helicity_init, only: helbookkeeping_flip, helsync_flip
  use ol_hel_propagators_/**/REALKIND
  use ol_hel_wavefunctions_/**/REALKIND
  use ol_wavefunctions_/**/REALKIND, only: wf_V_Std
! 
  use ol_hel_vertices_/**/REALKIND
  use ol_hel_contractions_/**/REALKIND
  use ol_external_ppllj_nexecsx_1, only: &
    & external_perm_ppllj_nexecsx_1, &
    & external_perm_inv_ppllj_nexecsx_1, &
    & extcomb_perm_ppllj_nexecsx_1, &
    & average_factor_ppllj_nexecsx_1
  use ol_external_ppllj_nexecsx_1, only: &
    & H, hel_not_initialised, hel_init, POLSEL 
  use ol_colourmatrix_ppllj_nexecsx_1_/**/REALKIND, only: &
    & colmat_not_initialised, colourmatrix_init
  use ol_forced_parameters_ppllj_nexecsx_1_/**/REALKIND, only: &
    & check_forced_parameters
  use ol_heltables_OLR_ppllj_nexecsx_1
  use ol_kinematics_/**/REALKIND, only: LC2Std_Rep_cmplx
  use ol_h_counterterms_/**/REALKIND
  use ol_loop_parameters_decl_/**/REALKIND ! counterterms
  use ol_loop_parameters_decl_/**/DREALKIND, only: &
    & IR_is_on, DOI, CT_is_on, R2_is_on, TP_is_on
  use ol_init, only: set_parameter, parameters_flush
#ifdef PRECISION_dp
  use ol_kinematics_/**/QREALKIND, only: conv_mom_scatt2in_qp=>conv_mom_scatt2in, &
                                         internal_momenta_qp=>internal_momenta
  use ol_momenta_decl_/**/DREALKIND, only: L
  use ol_momenta_decl_/**/QREALKIND, only: L_qp=>L
#endif
  implicit none

  real(DREALKIND), intent(in)  :: P_scatt(0:3,4)
  real(REALKIND),  intent(out) :: M02
  logical, intent(in) :: qp_kinematics
  real(REALKIND),  intent(out), optional :: M2ct
  real(REALKIND),  intent(out), optional :: M2colint(11)
  real(REALKIND) :: iM2ct
  real(REALKIND) :: iM2colint(11)

  integer(intkind1), save :: ntry = 0

  integer           :: shift, k, r, m, n, i
  real(REALKIND)    :: P(0:3,4)
#ifdef PRECISION_dp
  real(QREALKIND)    :: P_qp(0:3,4)
#endif
  integer           :: extmasses2(4)
  real(REALKIND)    :: M2add, M2add_ct, M2add_colint(11)
  complex(REALKIND) :: M1(1), M2(1)
  real(REALKIND)    :: P_scatt_intern(0:3,4)
  real(REALKIND), save :: scalebackfactor, old_scalefactor = 0
  integer(intkind1) :: nsync
  integer, allocatable :: extcombs_permuted(:)
  integer              :: extcombs(11), nextcombs
  integer              :: CT_on_bak, R2_on_bak, TP_on_bak, DOI_bak
  logical              :: do_ct, do_colint
#if 4 > 3
  integer :: ind_cc_comb(2)
#endif
  complex(REALKIND) :: omega(2) ! phases for helicity correlations

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
    allocate(A(nhel,2))
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

  extmasses2 = [ 0, 0, nMC, 0 ]
  ! Convert 2 -> n-2 PS-point to n -> 0 (so that P(1) + ... + P(n) = 0) and compute
  ! internal-propagator momenta in light-cone representation
  call init_kinematics(P_scatt, extmasses2, P, &
     external_perm_inv_ppllj_nexecsx_1, 4, qp_kinematics)

  ! denominators
  call denominators()

   if (heltables_not_init) call init_heltables()

  ! external WFs
  ! Here the external wavefunctions are initialiased
  call pol_wf_A(P(:,1), rZERO, H1, ex1, POLSEL(1),1)
  call pol_wf_Q(P(:,2), rZERO, H2, ex2, POLSEL(2),2)
  call pol_wf_Q(P(:,3), rMC, H3, ex3, POLSEL(3),3)
  call pol_wf_A(P(:,4), rZERO, H4, ex4, POLSEL(4),4)



  ! internal WFs
  ! e.g. call vert_VQ_A(ntry, ex3, ex1, wf1, n1, t1) ...
  call vert_QA_W(ntry, ex2(:), ex1(:), wf4(:,1), n3(:,1), t3x4(:,:,1))
  call vert_QA_W(ntry, ex3(:), ex4(:), wf4(:,2), n3(:,2), t3x4(:,:,2))
  call prop_W_W(ntry, wf4(:,1), 3, MW, 1_intkind1, wf4(:,3), n2(1))
  call counter_QA_W(ntry, ex3(:), ex4(:), wf4(:,4), n3(:,3), t3x4(:,:,3))


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

#if 1 > 0 && 4 > 3
  if (loopcc) then
    !!Born-loop colour correlators interference
    ind_cc_comb = [2,4]
    M0M1_hel_cc(:,:,:)%j = 0
    M2ctcc = 0
    do n = 1, 2
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
#endif

  if (ntry < 1) then
    if (allocated(A)) deallocate(A)
  end if
  ntry=2

  contains

subroutine physical_helicities()
  implicit none

integer :: i1,i2,i3,i4
integer :: n

n=0

do i4= 1, 2
do i3= 1, 2
do i2= 1, 2
do i1= 1, 2
n = n + 1
 if((ex1(i1)%hf==-1_intkind2) .OR. (ex2(i2)%hf==-1_intkind2) .OR. (ex3(i3)%hf==-1_intkind2) .OR. (ex4(i4)%hf==-1_intkind2)) then
  Hel(n) = -1_intkind2
 else
  Hel(n) = ex1(i1)%hf + ex2(i2)%hf + ex3(i3)%hf + ex4(i4)%hf
 end if
end do
end do
end do
end do

end subroutine physical_helicities

subroutine diagrams()
  implicit none
  integer :: h
  ! e.g. call cont_VV(nsync, wf3, wf6, A(:,1), n64, t64, nhel, den(5)) ...

    call Hcont_VV(nsync, wf4(:,2), wf4(:,3), A(:,1), n3(:,4), t3x16(:,:,1), nhel, den(1))

    call Hcont_VV(nsync, wf4(:,3), wf4(:,4), A(:,2), n3(:,5), t3x16(:,:,2), nhel, den(1))

end subroutine diagrams


elemental function diagmap(j, n)
  implicit none
  integer, intent(in) :: j, n
  complex(REALKIND) :: diagmap
  diagmap = A(j,n)%j
end function diagmap

function diagsum(j, pos, neg)
  implicit none
  integer, intent(in) :: j, pos(:), neg(:)
  complex(REALKIND) :: diagsum
  diagsum = sum(diagmap(j, pos)) - sum(diagmap(j, neg))
end function diagsum

subroutine colourvectors(A, j, M1, M2)
  implicit none
  type(Hpolcont) :: A(:,:)
  integer, intent(in) :: j
  type(Hpolcont), intent(out) :: M1(1), M2(1) ! M1helarray(1,nhel)
  integer :: empty(0), i

  M1(1)%j = -(A(j,1)%j*f(1))

  M2(1)%j = -(A(j,2)%j*f(2))


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
  use ol_colourmatrix_ppllj_nexecsx_1_/**/REALKIND, only: K1
  implicit none

  complex(REALKIND), intent(in)  :: M(1)
  real(REALKIND),    intent(out) :: M2colint
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
  use ol_colourmatrix_ppllj_nexecsx_1_/**/REALKIND, only: K2
  implicit none
  complex(REALKIND), intent(in)  :: M(1), Mct(1)
  real(REALKIND),    intent(out) :: M2colint_ct
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
  use ol_colourmatrix_ppllj_nexecsx_1_/**/REALKIND, only: K1
  implicit none

  complex(REALKIND), intent(in)  :: M(1)
  real(REALKIND),    intent(out) :: M2IRadd(11)
  integer ::  i, j, k, colmatpos

  M2IRadd = 0

  do k = 1, 12-1
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
  use ol_colourmatrix_ppllj_nexecsx_1_/**/REALKIND, only: K2
  implicit none

  complex(REALKIND), intent(in)  :: M1(1)
  complex(REALKIND), intent(in)  :: M2(1)
  complex(REALKIND),    intent(out) :: M2colint
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
  use ol_colourmatrix_ppllj_nexecsx_1_/**/REALKIND, only: KL
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
#ifdef PRECISION_dp
subroutine hsfamp2(P_scatt, M02, pind_1, pind_2, M2hsf)
#else
subroutine hsfamp2(P_scatt, M02, pind_1, pind_2, M2hsf)
  use ol_loop_storage_ppllj_nexecsx_1_/**/DREALKIND, only: &
    & nhflip, helflip, hflip, tsb, pi_flip_bak, pj_flip_bak, nhel
#endif
! M2tree = helicity-summed squared tree matrix element for anti-nu_e e- charm anti-strange -> 0
! M2hsf  = formula (66) of 1011.3918. Tree level
! spin-correlated hard scattering factor for anti-nu_e e- charm anti-strange -> 0
! **********************************************************************
  use ol_parameters_decl_/**/REALKIND
  implicit none
  real(DREALKIND), intent(in) :: P_scatt(0:3,4)
  integer, intent(in) :: pind_1, pind_2
  real(REALKIND), intent(out) :: M02, M2hsf
  real(REALKIND)    :: M2add, M2ct_dummy
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
  use KIND_TYPES, only: intkind2
  implicit none
  integer(intkind2), intent(in) :: ext_wf_hf(:), ext_wf_t
  integer, intent(in) :: m
  if(.not. allocated(hflip)) allocate(hflip(size(ext_wf_hf),2))
  if(size(ext_wf_hf) == size(hflip,1)) then
    hflip(:,m) = ext_wf_hf
    tsb = tsb + ext_wf_t
  else
    call ol_fatal("Spin-correlated HSF for ppllj_nexecsx_1:" &
      & // "trying to flip particles with different number of helicity states")
  end if
end subroutine assign_helicity_flip

subroutine init_flipped_heltables(i,j)
  use KIND_TYPES, only: intkind2
  use ol_external_ppllj_nexecsx_1, only: &
    & external_perm_inv_ppllj_nexecsx_1
  use ol_h_helicity_bookkeeping_/**/REALKIND, only: &
    & helicity_flip_ij
  implicit none
  integer, intent(in) :: i, j
  integer :: r, k, l=0
  tsb = 0

  if(allocated(hflip)) deallocate(hflip)

  do k = 1, 4
    r = external_perm_inv_ppllj_nexecsx_1(k)
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
      end select
    end if
  end do
  call helicity_flip_ij(16,nhel, M1helarray(1,:)%hf, hflip, tsb, &
      & nhflip, helflip)

end subroutine init_flipped_heltables

! **********************************************************************
subroutine colintHSF(M, Mflip, Mhfcolint)
! **********************************************************************
  use ol_colourmatrix_ppllj_nexecsx_1_/**/REALKIND, only: K1
  implicit none
  complex(REALKIND), intent(in)  :: M(1), Mflip(1)
  real(REALKIND),    intent(out) :: Mhfcolint
  integer :: i, j
  MHfcolint = 0
  do i = 1, 1
    do j = 1, 1
      Mhfcolint = Mhfcolint + real(conjg(Mflip(i))*K1(i,j)*M(j))
    end do
  end do

end subroutine colintHSF

end subroutine hsfamp2

end module ol_loop_ppllj_nexecsx_1_/**/REALKIND
