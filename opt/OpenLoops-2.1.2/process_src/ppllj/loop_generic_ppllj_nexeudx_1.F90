module ol_external_ppllj_nexeudx_1
! 
  implicit none
  integer :: dummy_counter
  ! Permutation and inverse permutation of external particles
  integer, save :: external_perm_ppllj_nexeudx_1(4) = &
                     [ (dummy_counter, dummy_counter = 1, 4) ]
  integer, save :: external_perm_inv_ppllj_nexeudx_1(4) = &
                     [ (dummy_counter, dummy_counter = 1, 4) ]
  integer, save :: extcomb_perm_ppllj_nexeudx_1(0:11) = &
                     [ (dummy_counter, dummy_counter = 0, 11) ]
  ! Particle types (mapping of fields to integers is not fixed!)
  integer, save :: particle_types_ppllj_nexeudx_1(4) = &
                     [ 1, 2, 3, 4 ]
  ! Colour and helicity average factors per particle
  integer, save :: average_factors_ppllj_nexeudx_1(4) = &
                     [ 2, 2, 6, 6 ]
  ! Average factor; initialised to the identity permutation
  integer, save :: average_factor_ppllj_nexeudx_1 = &
                     4
  integer, save :: channel_number_ppllj_nexeudx_1 = -1
  ! external particle helicities
  logical, save :: hel_not_initialised = .true.
  integer, save :: H(4,16) ! H(i,la) = helicity of particle i in configuration la
  integer, save :: H_HC(16,4)
  integer, save :: POLSEL(4) = 0
  integer, save :: photonid_ppllj_nexeudx_1(4) = 0

  
  contains

  subroutine n_external(n) &
      & bind(c,name="ol_f_n_external_ppllj_nexeudx_1")
    ! Return the number of external particles
    implicit none
    integer, intent(out) :: n
    n = 4
  end subroutine n_external


  subroutine n_external_c(n) &
      & bind(c,name="ol_n_external_ppllj_nexeudx_1")
    use, intrinsic :: iso_c_binding, only: c_int
    implicit none
    integer(c_int) :: n
    n = 4
  end subroutine n_external_c

  subroutine averagefactor_c(avgf) &
      & bind(c,name="ol_averagefactor_ppllj_nexeudx_1")
    use KIND_TYPES, only: DREALKIND
    implicit none
    real(DREALKIND), intent(out) :: avgf
    avgf = average_factor_ppllj_nexeudx_1
  end subroutine averagefactor_c

  subroutine set_permutation(perm) &
      & bind(c,name="ol_f_set_permutation_ppllj_nexeudx_1")
    use ol_parameters_decl_/**/DREALKIND, only: out_symmetry_on
    use ol_external_decl_/**/DREALKIND, only: n_scatt
    use ol_generic, only: factorial
    implicit none
    integer, intent(in) :: perm(4)
    integer :: i, j, ii, jj
    integer :: particle_types_perm_ppllj_nexeudx_1(4)
    external_perm_ppllj_nexeudx_1 = perm
    do i = 1, 4
      external_perm_inv_ppllj_nexeudx_1( &
        external_perm_ppllj_nexeudx_1(i)) = i
      particle_types_perm_ppllj_nexeudx_1(i) = &
        particle_types_ppllj_nexeudx_1( &
        external_perm_ppllj_nexeudx_1(i))
    end do
    do i = 1, 4
      do j = 1, i
        if (external_perm_ppllj_nexeudx_1(i) >= &
          external_perm_ppllj_nexeudx_1(j)) then
          ii = external_perm_ppllj_nexeudx_1(i)
          jj = external_perm_ppllj_nexeudx_1(j)
        else
          ii = external_perm_ppllj_nexeudx_1(j)
          jj = external_perm_ppllj_nexeudx_1(i)
        end if
        extcomb_perm_ppllj_nexeudx_1((i*(i-1))/2 + j) = (ii*(ii-1))/2 + jj
      end do
    end do
    ! Colour and helicity average factor
    average_factor_ppllj_nexeudx_1 = 1
    do i = 1, n_scatt
      average_factor_ppllj_nexeudx_1 = &
        average_factor_ppllj_nexeudx_1 &
        * average_factors_ppllj_nexeudx_1( &
        external_perm_ppllj_nexeudx_1(i))
    end do
    ! Symmetry factor for outgoing particles
    if (out_symmetry_on /= 0) then
      do i = 1, 4
        average_factor_ppllj_nexeudx_1 = &
          average_factor_ppllj_nexeudx_1 &
          * factorial(count(particle_types_perm_ppllj_nexeudx_1(n_scatt+1:4) == i))
      end do
    end if
  end subroutine set_permutation


  subroutine set_permutation_c(perm) &
      & bind(c,name="ol_set_permutation_ppllj_nexeudx_1")
    use, intrinsic :: iso_c_binding, only: c_int
    implicit none
    integer(c_int), intent(in) :: perm(4)
    integer :: f_perm(4)
    f_perm = perm
    call set_permutation(f_perm)
  end subroutine set_permutation_c


  subroutine get_last_perm(perm) &
      & bind(c,name="ol_f_get_last_perm_ppllj_nexeudx_1")
    integer, intent(out) :: perm(4)
    perm = external_perm_ppllj_nexeudx_1
  end subroutine get_last_perm


  subroutine get_masses(m_ex) &
      & bind(c,name="ol_f_get_masses_ppllj_nexeudx_1")
    ! Return the masses of the external particles in the current permutation.
    use KIND_TYPES, only: DREALKIND
    use ol_parameters_decl_/**/DREALKIND
    implicit none
    real(DREALKIND), intent(out) :: m_ex(4)
    integer        :: i
    real(DREALKIND) :: m_ex_orig(4)
    ! External particle masses for in the identity permutation
    m_ex_orig = [ rZERO, rZERO, rZERO, rZERO ]
    do i = 1, 4
      m_ex(i) = m_ex_orig(external_perm_ppllj_nexeudx_1(i))
    end do
  end subroutine get_masses


  subroutine get_masses_c(m_ex) &
      & bind(c,name="ol_get_masses_ppllj_nexeudx_1")
    use KIND_TYPES, only: DREALKIND
    use, intrinsic :: iso_c_binding, only: c_double
    implicit none
    real(c_double), intent(out) :: m_ex(4)
    real(DREALKIND) :: f_m_ex(4)
    call get_masses(f_m_ex)
    m_ex = f_m_ex
  end subroutine get_masses_c


  subroutine get_types(t_ex) &
      & bind(c,name="ol_f_get_types_ppllj_nexeudx_1")
    ! Return the type of the external particles in the current permutation.
    use ol_parameters_decl_/**/DREALKIND
    implicit none
    integer, intent(out) :: t_ex(4)
    integer        :: i
    integer :: t_ex_orig(4)
    ! External particle masses for in the identity permutation
    t_ex_orig = [0,3,2,2]
    do i = 1, 4
      t_ex(i) = t_ex_orig(external_perm_ppllj_nexeudx_1(i))
    end do
  end subroutine get_types


  subroutine get_types_c(t_ex) &
      & bind(c,name="ol_get_types_ppllj_nexeudx_1")
    use, intrinsic :: iso_c_binding, only: c_int
    use, intrinsic :: iso_c_binding, only: c_double
    implicit none
    integer(c_int), intent(out) :: t_ex(4)
    integer :: f_t_ex(4)
    call get_types(f_t_ex)
    t_ex = f_t_ex
  end subroutine get_types_c


  subroutine get_charges(c_ex) &
      & bind(c,name="ol_f_get_charges_ppllj_nexeudx_1")
    ! Return the type of the external particles in the current permutation.
    use KIND_TYPES, only: DREALKIND
    use ol_parameters_decl_/**/DREALKIND
    implicit none
    real(DREALKIND), intent(out) :: c_ex(4)
    integer        :: i
    real(DREALKIND) :: c_ex_orig(4)
    ! External particle masses for in the identity permutation
    c_ex_orig = [0,-3,2,1]/3._/**/REALKIND
    do i = 1, 4
      c_ex(i) = c_ex_orig(external_perm_ppllj_nexeudx_1(i))
    end do
  end subroutine get_charges


  subroutine get_charges_c(c_ex) &
      & bind(c,name="ol_get_charge_ppllj_nexeudx_1")
    use KIND_TYPES, only: DREALKIND
    use, intrinsic :: iso_c_binding, only: c_double
    implicit none
    real(c_double), intent(out) :: c_ex(4)
    real(DREALKIND) :: f_c_ex(4)
    call get_charges(f_c_ex)
    c_ex = f_c_ex
  end subroutine get_charges_c



  subroutine rambo(sqrt_s, p_rambo) &
      & bind(c,name="ol_f_rambo_ppllj_nexeudx_1")
    use KIND_TYPES, only: DREALKIND
    use ol_kinematics_/**/DREALKIND, only: rambo_generic => rambo
    implicit none
    real(DREALKIND), intent(in) :: sqrt_s
    real(DREALKIND), intent(out) :: p_rambo(0:3,4)
    real(DREALKIND) :: m_ex(4)
    call get_masses(m_ex)
    call rambo_generic(sqrt_s, m_ex, p_rambo)
  end subroutine rambo


  subroutine rambo_c(sqrt_s, p_rambo) &
      & bind(c,name="ol_rambo_ppllj_nexeudx_1")
    use KIND_TYPES, only: DREALKIND
    use, intrinsic :: iso_c_binding, only: c_double
    implicit none
    real(c_double), intent(in) :: sqrt_s
    real(c_double), intent(out) :: p_rambo(0:3,4)
    real(DREALKIND) :: f_sqrt_s
    real(DREALKIND) :: f_p_rambo(0:3,4)
    f_sqrt_s = sqrt_s
    call rambo(f_sqrt_s, f_p_rambo)
    p_rambo = f_p_rambo
  end subroutine rambo_c


  subroutine hel_init
    implicit none
    integer :: binpos, flip, binco
    hel_not_initialised = .false.
    ! helicity configurations for this process
  H(:, 1) = [ -1, -1, -1, -1 ]
  H(:, 2) = [ -1, -1, -1,  1 ]
  H(:, 3) = [ -1, -1,  1, -1 ]
  H(:, 4) = [ -1, -1,  1,  1 ]
  H(:, 5) = [ -1,  1, -1, -1 ]
  H(:, 6) = [ -1,  1, -1,  1 ]
  H(:, 7) = [ -1,  1,  1, -1 ]
  H(:, 8) = [ -1,  1,  1,  1 ]
  H(:, 9) = [  1, -1, -1, -1 ]
  H(:,10) = [  1, -1, -1,  1 ]
  H(:,11) = [  1, -1,  1, -1 ]
  H(:,12) = [  1, -1,  1,  1 ]
  H(:,13) = [  1,  1, -1, -1 ]
  H(:,14) = [  1,  1, -1,  1 ]
  H(:,15) = [  1,  1,  1, -1 ]
  H(:,16) = [  1,  1,  1,  1 ]

  end subroutine hel_init


  subroutine pol_init(pol) &
      & bind(c,name="ol_f_pol_init_ppllj_nexeudx_1")
    implicit none
    integer, intent(in) :: pol(4)
    POLSEL = pol
  end subroutine pol_init


  subroutine set_photons(photon_id) &
      & bind(c,name="ol_f_set_photons_ppllj_nexeudx_1")
    implicit none
    integer, intent(in) :: photon_id(4)
    photonid_ppllj_nexeudx_1 = photon_id
  end subroutine set_photons




end module ol_external_ppllj_nexeudx_1


module ol_settings_ppllj_nexeudx_1
  implicit none
  ! Activates optimized helicity bookkeeping. Default=false
    logical, save :: hel_mem_opt = .false.
  ! Calculate loop correlators
    logical, save :: loopcc = .false.

  contains

  subroutine set_hel_mem_opt(set) &
      & bind(c,name="ol_hel_mem_opt_ppllj_nexeudx_1")
    implicit none
    logical, intent(in) :: set
    hel_mem_opt = set
  end subroutine set_hel_mem_opt

  subroutine set_loopcc(set,was) &
      & bind(c,name="ol_loopcc_ppllj_nexeudx_1")
    implicit none
    logical, intent(in) :: set
    logical :: was
    was = loopcc
    loopcc = set
  end subroutine set_loopcc

end module ol_settings_ppllj_nexeudx_1

module colour_basis_ppllj_nexeudx_1
  implicit none
  ! tree colour basis
  integer, save :: extcolours(4) = [0,0,1,1]
  contains

  pure subroutine tree_colbasis_dim(extcols, ncolb, ncoupl, maxpows, nhel) &
    & bind(c, name="ol_tree_colbasis_dim_ppllj_nexeudx_1")
    implicit none
    ! colour representation of external particles: 0=neutral, 1=fundamental, 2=adjoint
    integer, intent(out) :: extcols(4)
    ! number of tree colour basis elements; number of selected couplings, number of selected powers per coupling
    integer, intent(out) :: ncolb, ncoupl, maxpows
    ! number of helicity configurations (all, not just non-vanishing)
    integer, intent(out) :: nhel
    extcols = extcolours
    ncolb = 1
    ncoupl = 1
    maxpows = 1
    nhel = 16
  end subroutine tree_colbasis_dim

  subroutine tree_colbasis(basis, powers) &
    & bind(c, name="ol_tree_colbasis_ppllj_nexeudx_1")
    implicit none
    integer, intent(out) :: powers(1,1)
    integer, intent(out) :: basis(3,1)
#if 1 > 0
    ! selected powers for each selected coupling
    powers = reshape([4], [1,1])
#endif
#if 1 > 0
    ! tree colour basis: [[composition_number, permutation_number, *coupling_powers], ...]
    basis = reshape( &
      [1,2,2], &
      [3,1])
#endif
  end subroutine tree_colbasis


  pure subroutine loop_colbasis_dim(extcols, ncolb, ncoupl, maxpows, nhel) &
    & bind(c, name="ol_loop_colbasis_dim_ppllj_nexeudx_1")
    implicit none
    ! colour representation of external particles: 0=neutral, 1=fundamental, 2=adjoint
    integer, intent(out) :: extcols(4)
    ! number of loop colour basis elements; number of selected couplings, number of selected powers per coupling
    integer, intent(out) :: ncolb, ncoupl, maxpows
    ! number of helicity configurations (all, not just non-vanishing)
    integer, intent(out) :: nhel
    extcols = extcolours
    ncolb = 1
    ncoupl = 1
    maxpows = 1
    nhel = 16
  end subroutine loop_colbasis_dim

  subroutine loop_colbasis(basis, powers) &
    & bind(c, name="ol_loop_colbasis_ppllj_nexeudx_1")
    implicit none
    integer, intent(out) :: powers(1,1)
    integer, intent(out) :: basis(3,1)
#if 1 > 0
    ! selected powers for each selected coupling
    powers = reshape([4], [1,1])
#endif
#if 1 > 0
    basis = 0 ! TODO
    ! loop colour basis: [[composition_number, permutation_number, *coupling_powers], ...]
!    basis = reshape( &
!<\*colourbasis\*>, &
!      [3,1])
#endif
  end subroutine loop_colbasis


end module colour_basis_ppllj_nexeudx_1

! **********************************************************************
module ol_heltables_OLR_ppllj_nexeudx_1
! **********************************************************************
  use KIND_TYPES, only: intkind2
  implicit none

  logical :: heltables_not_init = .true.

  ! helicity states of external particles
  ! integer, save :: &
  !   H1(2) = [-1,1], &
  !   H2(3) = [-1,0,1]
  !   ...
  integer, save :: &
    H1(2) = [-1,1], &
    H2(2) = [-1,1], &
    H3(2) = [-1,1], &
    H4(2) = [-1,1]

  ! number of helicity states for wave functions returned by a propagator call: n2(sz)
  ! number of helicity states for wave functions in a v-point vertex call (v >= 3)
  ! or a contraction (v = 3): n<v>(v,sz)
  integer(intkind2), save :: n2(1), n3(3,5)

  ! helicity tables used in the construction of the h helicity states of a wave function (amplitude)
  ! from an v-point vertex (contraction): t<v>x<h>(v-1,h,sz)
  integer(intkind2), save :: t3x16(2,16,2), t3x4(2,4,3)

  ! change of global-helicity state resulting from flip of individual-particle helicity
  integer(intkind2), save :: eflip(16,4)
  integer,           save :: exthel(16,4)
  integer,           save :: firstpol(4)

  contains

! **********************************************************************
subroutine init_heltables
! **********************************************************************
  use ol_helicity_init, only: heltable
  implicit none

  ! I/O helicity tables for vertices, propagators and contractions;
  ! helicity table for a vertex call: n_in/n_out are the number helicity states of the incoming/outgoing wave functions
  ! call heltable([<n_in1>, <n_in2>, ..., <n_out>], n, t)
  ! propagators only need the number of helicity configurations which is equal for the incoming and outgoing wave function
  ! n = <n>
  call heltable([2,2,4], n3(:,1), t3x4(:,:,1))
  call heltable([2,2,4], n3(:,2), t3x4(:,:,2))
  n2(1) = 4
  call heltable([2,2,4], n3(:,3), t3x4(:,:,3))
  call heltable([4,4,16], n3(:,4), t3x16(:,:,1))
  call heltable([4,4,16], n3(:,5), t3x16(:,:,2))

  heltables_not_init = .false.

end subroutine init_heltables

end module ol_heltables_OLR_ppllj_nexeudx_1
