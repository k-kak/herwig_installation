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

module ol_data_types_dp
  use kind_types, only: dp, qp, intkind1, intkind2

  use ol_data_types_qp, only: basis_qp=>basis, redset4_qp=>redset4

  type wfun
    ! four complex components for the wave function
    complex(dp) :: j(4)
    complex(dp), pointer :: j_prev(:)
    ! indicator if left- or right components of of-shell line vanish
    !                             j= (0,0,0,0) (0,0,j3,j4) (j1,j2,0,0) (j1,j2,j3,j4)
    integer(intkind1) :: h      ! B"00"      B"01"       B"10"        B"11"
    integer(intkind2) :: e      ! helicities of external on-shell lines
    integer(intkind2) :: t      ! label for the external subtree
    integer(intkind2) :: n_part ! number of particles in the subtree
    integer(intkind2) :: hf     ! global base-4 helicity label
  end type wfun

  type polcont
    complex(dp) :: j
    integer(intkind2) :: e ! helicities of external on-shell lines
    integer(intkind2) :: s ! table for final helicity syncronisation
  end type polcont

  ! open-loop derived data type used in the optimezed helicity summation
  type hol
    ! OpenLoops Coefficient: (alpha,rank,beta,helicity_state)
    complex(dp), dimension(:,:,:,:), allocatable :: j

    complex(qp), dimension(:,:,:,:), allocatable :: j_qp

    ! Helicity configurations array
    integer(intkind2), dimension(:)      , allocatable :: hf
    integer :: mode = 1
    real(dp) :: error
    integer :: npoint = 0
    integer :: ndrs = 0
    integer :: nred = 0

    integer :: ndrs_qp = 0
    integer :: nred_qp = 0

  end type hol

  ! derived tensor data type for closed-loop
  type hcl
    complex(dp) , dimension(:), allocatable :: cmp

    complex(qp) , dimension(:), allocatable :: cmp_qp

    integer :: mode = 1
    real(dp) :: error
    integer :: ndrs = 0
    integer :: nred = 0

    integer :: ndrs_qp = 0
    integer :: nred_qp = 0

  end type hcl

  type met
    real(dp) :: cmp

    real(qp) :: cmp_qp

    integer :: mode = 1
    real(dp) :: error
    integer :: sicount = 0
    integer :: ndrs = 0
    integer :: nred = 0

    integer :: sicount_qp = 0
    integer :: ndrs_qp = 0
    integer :: nred_qp = 0

  end type met

  ! equivalent to polcont, with the addition of an extra hf label for
  ! the global helicity state
  type Hpolcont
    complex(dp)  :: j
    integer(intkind2)  :: e  ! helicities of external on-shell lines
    integer(intkind2)  :: hf ! global base-4 helicity label
    integer(intkind2)  :: s  ! table for final helicity syncronisation
  end type Hpolcont

  !! l_i basis for the on-the-fly reduction
  type basis
    complex(dp) :: vect1(4)     !! l_{1,\mu} in light-cone rep
    complex(dp) :: vect2(4)     !! l_{2,\mu} in light-cone rep
    complex(dp) :: vect3(4)     !! l_{3,\mu} in light-cone rep
    complex(dp) :: vect4(4)     !! l_{4,\mu} in light-cone rep
    complex(dp) :: tens1(10)
    complex(dp) :: tens2(10)
    complex(dp) :: tens3(10)
    complex(dp) :: tens4(10,4)
    complex(dp) :: tens5(10,4)
    complex(dp) :: gamma
    complex(dp) :: alpha(2)
    integer           :: mom1
    integer           :: mom2
    complex(dp) :: li(4,4)
  end type basis

  !! Set containing the reduction basis
  type redset4
    type(basis) :: redbasis
    complex(dp) :: p3scalars(0:4)
    integer :: perm(3)
    integer :: mom3
    real(dp) :: gd2,gd3

    logical :: qp_computed = .false.
    type(redset4_qp) :: rsqp

  end type redset4

  type redset5
    type(basis) :: redbasis
    complex(dp) :: p3scalars(0:4)
    integer :: perm(4)
    integer :: mom3
    integer :: mom4
    real(dp) :: gd2,gd3
  end type redset5

  !! Scalar Box
  type scalarbox
    complex(dp) :: poles(0:2)         ! finite, eps^(-1), eps^(-2)
    complex(dp) :: onshell_cuts(2,5)  ! on-shell cuts of the box
    real(dp) :: cut_error
    real(dp) :: box_error

    integer :: mom_ind(3)
    integer :: mom1
    integer :: mom2
    integer :: mom3
    integer :: perm(3)
    real(dp)     :: gd3
    integer            :: masses2(0:3)
    logical            :: qp_computed = .false.
    complex(qp) :: poles_qp(0:2)         ! finite, eps^(-1), eps^(-2)
    complex(qp) :: onshell_cuts_qp(2,5)  ! on-shell cuts of the box

  end type scalarbox



  type carray2
    complex(dp), allocatable :: arr(:,:)
  end type

  type l2lc_rdata
    integer, allocatable :: r(:,:), c(:,:)
  end type l2lc_rdata

  type l2lc_data
    type(l2lc_rdata), allocatable :: arr(:)
  end type l2lc_data

  type me_cache
    real(dp), allocatable :: psp(:,:), me(:)
    integer :: loop_parameters_status=0
  end type me_cache

  type correlator
    integer :: type=0 ! ! 0: none, 1:cc, 2:sc, 3: bmunu
    integer :: emitter=0
    real(dp) :: mom(4)=0
    integer :: nextcombs=0
    integer, allocatable :: extcombs(:)
    real(dp), allocatable :: rescc(:)
    real(dp) :: resmunu(4,4)=0
  end type correlator

  contains

  subroutine zero_correlator(corr)
  implicit none
  type(correlator), intent(inout) :: corr
    if(allocated(corr%rescc)) corr%rescc = 0
    corr%resmunu = 0
  end subroutine zero_correlator




end module ol_data_types_dp



module ol_momenta_decl_dp
  use kind_types, only: dp, qp
  use ol_debug, only: ol_msg
  implicit none
  ! Internal momenta for external particles
  ! Components 1:4 = light cone representation; component 5 = squared momentum
  complex(dp), allocatable, save :: Q(:,:) ! Q(5,0:2**Nparticle-1)
  complex(dp), allocatable, save :: QInvariantsMatrix(:,:) ! QInvariantsMatrix(Nparticle,Nparticle)
  ! Components 1:4 = light cone representation; component 5 = squared masses, component 6 = scalar products
  complex(dp), allocatable, save :: L(:,:) ! L(6,0:2**Nparticle-1)

  complex(qp), allocatable, save :: Q_qp(:,:) ! Q(5,0:2**Nparticle-1)
  complex(qp), allocatable, save :: QInvariantsMatrix_qp(:,:) ! QInvariantsMatrix(Nparticle,Nparticle)
  complex(qp), allocatable, save :: L_qp(:,:) ! L(6,0:2**Nparticle-1)

  integer, save :: collconf = 0
  integer, save :: softconf = 0

  contains

  function momenta_nan_check(P)
    implicit none
    real(dp), intent(in) :: P(:,:)
    integer :: momenta_nan_check
    integer :: i
    if (all(P == P)) then
      ! ok
      momenta_nan_check = 0
    else
      ! contains NaN
      call ol_msg("=== WARNING ===")
      call ol_msg("corrupted phase space point:")
      do i = 1, size(P,2)
        write(*,*) P(:,i)
      end do
      momenta_nan_check = 1
      call ol_msg("===============")
    end if
  end function momenta_nan_check

end module ol_momenta_decl_dp


module ol_external_decl_dp
  use kind_types, only: dp
!  use ol_global_decl, only: MaxParticles
  implicit none
  ! phase space point cache; used to print the ps-point if tensor reduction fails
  integer,        save :: nParticles = 0 ! set by conv_mom_scatt2in
  integer, save :: allocatedNpart = 0
  integer, allocatable, save :: binom2(:)
  integer, allocatable, save :: crossing(:) ! only used if a reduction error occurs
  integer, allocatable, save :: inverse_crossing(:) ! set by conv_mom_scatt2in
  integer, allocatable, save :: gf_array(:)
  ! lists for each external particle the external momentum
  ! used for the gauge fixing of the vector polarization. Used in subroutine wf_gf_V.
  ! A zero entry means that it is not a massless vector particle.
  integer, allocatable, save :: Ward_array(:) ! select particle "i" for the Ward identity -> Ward_array(i) = 1
  real(dp), allocatable, save :: P_ex(:,:) ! uncleaned external 2->n-2 momenta, set by conv_mom_scatt2in
  integer, allocatable, save :: M_ex(:) ! external 2->n-2 mass ids, set by conv_mom_scatt2in

  ! number of incoming particles for phase space configuation and cleaning
  integer, save :: n_scatt = 2
  logical, save :: init_qp = .false.

end module ol_external_decl_dp



module ol_pseudotree_dp
  use kind_types, only: dp
  implicit none
  ! loop momentum in pseudo tree (standard representation)
  real(dp), save :: pseudotree_momentum(0:3) = [ 314.1592653589793_dp, 271.8281828459045_dp, 100._dp, 57.72156649015328_dp ]
  ! Wave functions for pseudo tree
  complex(dp), save :: exloop(4,2) = reshape([ 2.718_dp, 3.141_dp,  0.9159_dp, 1._dp,  &
                                               1._dp,    0.5772_dp, 1.618_dp,  1.282_dp], [ 4, 2 ])
end module ol_pseudotree_dp



module ol_tensor_storage_dp
  use kind_types, only: dp
  implicit none
  complex(dp), allocatable, save :: tensor_stored(:)
  integer, save :: rank_stored
  integer, save :: array_length_stored ! length of the array associated with rank_stored
  integer, save :: tensor_storage_maxrank = -1
end module ol_tensor_storage_dp



module ol_parameters_decl_dp
  ! Declarations and initial values for numerical and physical parameters like masses, widths, and couplings.
  ! Loading the module for the first time initialises the parameters with their default values (given in this module).
  ! Note that when a value has been changed using parameters_init(val=...) loading the module
  ! will not reset the value to its default.
  use kind_types, only: dp
  use ol_version, only: splash_todo, welcome_length ! only to expose to init_ui
!   use TI_call_interface
  implicit none
  ! Counted up by 1 each time parameters_init() is called
  integer, save :: parameters_status = 0

  ! flag for hybrid preset mode
  ! - 0: hp mode disabled
  ! - 1: hp mode for hard regions (default)
  ! - 2: hp mode for hard regions (default)
  integer, save :: hp_mode = 1

  ! flag of hybrid baseline mode
  integer, save :: hp_switch = 1

  ! writes a log on number of dressing and reduction steps (and of which are upraded to qp
  ! written to the points file
  integer, save :: write_hp_log = 1

  ! hp trigger thresholds in on-the-fly reduction
  real(dp), save :: hp_loopacc = 8._dp  ! loop accuracy target

  ! internal, do not touch change
  real(dp), save :: hp_err_thres = 8._dp  ! accumulated error threshold
  real(dp), save :: hp_step_thres = 0.5_dp ! step threshold
  real(dp), save :: hp_redset_gd3_thres = 3e-4_dp

  ! merging dp and qp channel the dp channel is upgraded to qp (which comes at zero cost.)
  logical, save :: hp_automerge = .true.
  ! QP trigger on small rank2 GD in on-the-fly reduction (threshold: hybrid_threshold)
  logical, save :: hp_gamma_trig = .false.
  ! QP trigger for missing triangle expansions.
  logical, save :: hp_metri_trig = .true.

  ! QP trigger on unstable IR triangles in reduction.
  logical, save :: hp_irtri_trig = .false.

  ! QP trigger for ir triangle diagrams
  logical, save :: hp_ir_trig = .false.

  ! Fake QP trigger: computes the loop (not trees) in QP
  integer, save :: hp_fake_trig = 0

  integer, save :: hp_check_box = 1

  ! Allocation mode for QP channel
  ! 0: allocate all OL/CL and initialize to zero
  ! 1: allocate all OL/CL and initialize to zero only when needed
  ! 2: allocate OL/CL on demand and initialize to zero only when needed,
  !    deallocate whenever possible
  ! 3: allocate OL/CL on demand and initialize to zero only when needed
  integer, save :: hp_alloc_mode = 3

  ! do not change.
  integer, parameter :: hybrid_zero_mode = 0
  integer, parameter :: hybrid_dp_mode = 1
  integer, parameter :: hybrid_qp_mode = 2
  integer, parameter :: hybrid_dp_qp_mode = 3

  ! hp stability log variables
  ! scalar integral counter
  integer, save :: hp_nsi = 0, hp_nsi_qp = 0
  ! dressing step counter
  integer, save :: hp_ndrs = 0, hp_ndrs_qp = 0
  ! reduction step counter
  integer, save :: hp_nred = 0, hp_nred_qp = 0
  ! highest on-the-fly reduction error
  real(dp), save :: hp_max_err = 0._dp

  real(dp), save :: max_error = 0._dp

  !  compute real bubble diagrams in counterterms.

  integer, save :: use_bubble_vertex = 1
  ! internal parameter which is set to 1 if use_bubble_vertex=1 and
  ! process not require ew renormalization
  integer, save :: bubble_vertex = 0

  !  use some hacks to stabilize ir events

  logical, save :: ir_hacks = .false.




  ! QP kinematics
  ! 0: always initialize  QP kinematics for hp_mode>0
  ! 1: initialize  QP kinematics only when needed (hp_mode>0)
  integer, save :: hp_qp_kinematics_init_mode = 1
  ! use quad-precision definition of momenta to compute invariants
  ! -> stable rescaling for collinear configurations
  integer, save :: sync_qp_kinematics = 1

  integer, save :: parameters_verbose = 0
  integer, parameter :: procname_length = 80
  character(procname_length) :: current_processname = 'none' ! set by vamp2generic()
  integer, parameter :: max_parameter_length = 255 ! used for stability_logdir, install_path,
                                                         ! contract_file, printparameter file
  integer, parameter :: max_parameter_name_length = 30 ! maximal length of parameter names in init routines
  ! 0: never, 1: on finish() call, 2: adaptive, 3: always
  integer, save :: stability_log = 0
  integer, save :: write_psp = 0 ! write out phase space points from vamp2generic is called
  integer, save :: ti_monitor = 0 ! 1: write squared matrix element contribution per tensor integral
                                  ! 2: also write tensor integral call arguments to a file.
  integer, save :: use_me_cache = 1
  character(len=max_parameter_length) :: stability_logdir = "stability_log"
  character(len=max_parameter_length) :: tmp_dir = "."
  character(len=max_parameter_length) :: allowed_libs = ""
  character(len=max_parameter_length) :: approximation = ""
  character(len=max_parameter_length) :: model = "sm"
  character(len=max_parameter_length) :: shopping_list = "OL_shopping.m"
  logical, save :: partial_normal_order = .false.
  logical, save :: write_shopping_list = .false.
  logical, save :: write_params_at_start = .false.
  logical, save :: stability_logdir_not_created = .true.
  logical, save :: nosplash = .false.
  logical, save :: apicheck = .true.
  character(16) :: pid_string ! 11 for pid, "-", 4 random characters
  ! OpenLoops installation path; used to locate info files and process libraries
  ! character(len=:), allocatable :: install_path ! gfortran 4.7: does not work in modules and type (though it does in subroutines)
  ! character(len=max_parameter_length) :: install_path = "path"
  include "install_path.inc"
  ! Mode for check_last_[...] in laststep and tensor integral routine in looproutines
  integer, save :: l_switch = 1, a_switch = 1, a_switch_rescue = 7, redlib_qp = 5
  ! switcher for helicity improvement in OFR
  logical, save :: hel_mem_opt_switch = .true.
  ! switchers for checking Ward identities at tree/loop level
  integer, save :: Ward_tree = 0
  integer, save :: Ward_loop = 0
  ! 0 = full colour, 1 = leading colour
  integer, save :: LeadingColour = 0
  ! divide by the symmetry factor of identical outgoing particles
  integer, save :: out_symmetry_on = 1
  ! use flavour mappings: 1: quark & lepton mapping, 1: only lepton mapping, 2: only quark mapping
  integer, save :: flavour_mapping_on = 1
  ! Running number of the next partonic channel which is initialised (by get_externel_<proc>)
  ! in the tensor library cache system
  integer, save :: next_channel_number = 1
! TODO disable coli cache when using hybrid mode
  integer, save :: coli_cache_use = 1
  logical, save :: no_collier_stop = .false.
  ! select alpha_QED input scheme: 0 = on-shell = alpha(0), 1 = G_mu, 2 = alpha(MZ)
  integer, save :: ew_scheme = 1
  ! select alpha_QED renormalization scheme: 0 = on-shell = alpha(0), 1 = G_mu, 2 = alpha(MZ)
  integer, save :: ew_renorm_scheme = 1
  ! select reg. scheme for off-shell external photons: 0: off, 1: gamma -> FF splittings in dimreg
  logical, save :: offshell_photons_lsz = .true.
  ! select reg. scheme for all external photons: 0: numerical, 1: dimreg
  logical, save :: delta_alphamz_dimreg = .false.
  ! select renorm scheme for on-shell external photons: 0: off, 1: a(0)/a(Gmu/MZ) + dZe LSZ shift
  logical, save :: onshell_photons_lsz = .true.
  ! switch on/off photon self energy
  logical, save :: photon_selfenergy = .true.
  ! coupling order
  integer :: coupling_qcd(0:1) = -1
  integer :: coupling_ew(0:1) = -1
  integer :: order_ew = -1
  integer :: order_qcd = -1
  integer :: loop_order_ew = -1
  integer :: loop_order_qcd = -1
  integer :: CKMORDER = 0
  integer :: QED = 0
  ! select scalar integral library for self energies: 0 = none, 1 = Coli, 3=OneLOop, 7 = DD
  ! automatically set according to redlib (redlib=1->1,7->7,other->3)
  integer, save :: se_integral_switch = 1
  integer, save :: do_ew_renorm = 0
  integer, save :: do_qcd_renorm = 1
  integer, save :: cms_on = 1
  ! select only specific polarization of external vector bosons
  ! 0: both, 1: only transverse, 2: only longitudinal
  integer, save :: select_pol_V = 0
  integer :: add_associated_ew = 0
  ! library loader: check online collection: yes/no
  logical :: check_collection = .true.
  ! OLMode: 0=OL1, 1=OL1+OFR helicity summation, 2=full OFR, 3=full OFR + hp, -1=auto=3,2,1,0
  integer, save :: OLmode = -1
  ! Auto-preset: preset=2 for OLmode=1,2 and preset=5 for OLmode=2, preset=3 for loop-induced
  logical, save :: auto_preset = .true.
  ! expert_mode: allows to set stability options manually
  real(dp), save :: psp_tolerance = 1.e-9
  logical, save :: no_cleaning =  .false.
  logical, save :: cleaning_via_hardness =  .false.
  ! wf_V_select: select external vector boson wavfunction, 1=default, 2=ABC, 3: MG
  integer, save :: wf_v_select = 1

  logical, save :: expert_mode = .false.

! 1


  ! Numerical constants
  real(dp),    parameter :: rONE   = 1
  real(dp),    parameter :: rZERO  = 0
  real(dp),    parameter :: rZERO2 = 0
  real(dp),    parameter :: pi     = acos(-1._dp)
  real(dp),    parameter :: pi2_6  = (pi**2)/6
  real(dp),    parameter :: sqrt2  = sqrt(2._dp)
  real(dp),    parameter :: sqrt05 = sqrt(0.5_dp)
  complex(dp), parameter :: cONE   = 1
  complex(dp), parameter :: ZERO   = 0
  complex(dp), parameter :: ZERO2  = 0
  complex(dp), parameter :: CI     = (0._dp, 1._dp)
  complex(dp) :: integralnorm = CI/(16*pi**2)
  complex(dp) :: countertermnorm = 1._dp/(16._dp*pi**2)

  ! scale factor for dimensionful parameters
  real(dp), save :: scalefactor = 1._dp
  logical,        save :: reset_scalefactor = .false.
  integer,        save :: scaling_mode = 1 ! 1: reduction only, 3: everything

  ! synchronise Yukawa masses with masses
  logical, save :: yuk_from_mass = .true.

  ! Particle masses and widths
  real(dp), save :: rME_unscaled = 0,                    wME_unscaled = 0 ! electron mass and width
  real(dp), save :: rMM_unscaled = 0,                    wMM_unscaled = 0 ! muon mass and width
  real(dp), save :: rML_unscaled = 0,                    wML_unscaled = 0 ! tau mass and width
  real(dp), save :: rMU_unscaled = 0,                    wMU_unscaled = 0 ! up-quark mass and width
  real(dp), save :: rMD_unscaled = 0,                    wMD_unscaled = 0 ! down-quark mass and width
  real(dp), save :: rMS_unscaled = 0,                    wMS_unscaled = 0 ! strange-quark mass and width
  real(dp), save :: rMC_unscaled = 0,                    wMC_unscaled = 0 ! charm-quark mass and width
  real(dp), save :: rMB_unscaled = 0._dp,      wMB_unscaled = 0 ! bottom-quark mass and width
  real(dp), save :: rMT_unscaled = 172._dp,    wMT_unscaled = 0 ! top-quark mass and width
  real(dp), save :: rYE_unscaled = 0,                    wYE_unscaled = 0
  real(dp), save :: rYM_unscaled = 0,                    wYM_unscaled = 0
  real(dp), save :: rYL_unscaled = 0,                    wYL_unscaled = 0
  real(dp), save :: rYU_unscaled = 0,                    wYU_unscaled = 0
  real(dp), save :: rYD_unscaled = 0,                    wYD_unscaled = 0
  real(dp), save :: rYS_unscaled = 0,                    wYS_unscaled = 0
  real(dp), save :: rYC_unscaled = 0,                    wYC_unscaled = 0
  real(dp), save :: rYB_unscaled = 0,                    wYB_unscaled = 0
  real(dp), save :: rYT_unscaled = 172._dp,    wYT_unscaled = 0
  real(dp), save :: rMW_unscaled = 80.399_dp,  wMW_unscaled = 0 ! W boson mass LEP PDG 2008/2009 and width
  real(dp), save :: rMZ_unscaled = 91.1876_dp, wMZ_unscaled = 0 ! Z boson mass LEP PDG 2008/2009 and width
  real(dp), save :: rMX_unscaled = 0._dp,  wMX_unscaled = 0._dp ! auxiliary field for Z
  real(dp), save :: rMY_unscaled = 0._dp,  wMY_unscaled = 0.0000_dp ! auxiliary field for W
  real(dp), save :: rMH_unscaled = 125._dp,    wMH_unscaled = 0 ! higgs boson mass and width
  real(dp), save :: MREG_unscaled = 1._dp                       ! collinear mass regulator for photon WF CT
  ! Coupling constants

  real(dp), save :: alpha_QCD = 0.1258086856923967_dp ! LO MRST
  real(dp), save :: alpha_QED_MZ = 1/128._dp          ! alpha(MZ) derived from PDG 2014
  real(dp), save :: alpha_QED_0  = 1/137.035999074_dp  ! alpha(0) from PDG 2014
  real(dp), save :: alpha_QED, alpha_QED_input
  real(dp), save :: alpha_QED_Gmu
  real(dp), save :: sw2_input = 0.222626515643872389077366863865260666_dp

  real(dp), save :: Gmu_unscaled = 0.0000116637_dp    ! G_mu

  real(dp), save :: Gmu
  ! Everything beyond this line is derived from the values given above and initialised by parameters_init().
  real(dp), save :: rescalefactor = 1.1_dp
  ! scaled masses, widths and yukawas
  real(dp), save :: rME, wME, rYE, wYE
  real(dp), save :: rMM, wMM, rYM, wYM
  real(dp), save :: rML, wML, rYL, wYL
  real(dp), save :: rMU, wMU, rYU, wYU
  real(dp), save :: rMD, wMD, rYD, wYD
  real(dp), save :: rMS, wMS, rYS, wYS
  real(dp), save :: rMC, wMC, rYC, wYC
  real(dp), save :: rMB, wMB, rYB, wYB
  real(dp), save :: rMT, wMT, rYT, wYT
  real(dp), save :: rMW, wMW
  real(dp), save :: rMZ, wMZ
  real(dp), save :: rMH, wMH
  real(dp), save :: rMX, wMX
  real(dp), save :: rMY, wMY
  ! Masses ids
  integer, parameter :: nME = 1
  integer, parameter :: nMM = 2
  integer, parameter :: nML = 3
  integer, parameter :: nMU = 4
  integer, parameter :: nMD = 5
  integer, parameter :: nMS = 6
  integer, parameter :: nMC = 7
  integer, parameter :: nMB = 8
  integer, parameter :: nMT = 9
  integer, parameter :: nMW = 10
  integer, parameter :: nMZ = 11
  integer, parameter :: nMH = 12
  integer, parameter :: nMX = 13
  integer, parameter :: nMY = 14
  ! Complex masses, complex and real squared masses
  complex(dp), save ::  ME,   MM,   ML,   MU,   MD,   MS,   MC,   MB,   MT,   MW,   MZ,   MH,   MX,   MY
  complex(dp), save ::  ME2,  MM2,  ML2,  MU2,  MD2,  MS2,  MC2,  MB2,  MT2,  MW2,  MZ2,  MH2,  MX2,  MY2
  complex(dp), save ::  YE,   YM,   YL,   YU,   YD,   YS,   YC,   YB,   YT
  complex(dp), save ::  YE2,  YM2,  YL2,  YU2,  YD2,  YS2,  YC2,  YB2,  YT2
  real(dp),    save :: rYE2, rYM2, rYL2, rYU2, rYD2, rYS2, rYC2, rYB2, rYT2
  real(dp),    save :: rME2, rMM2, rML2, rMU2, rMD2, rMS2, rMC2, rMB2, rMT2, rMW2, rMZ2, rMH2, rMX2, rMY2
  complex(dp), save :: YC2pair, YB2pair, YT2pair ! pair masses: only non-zero if the SU(2) partner is active
  ! collinear mass regulator for photon WF CT
  real(dp),    save :: MREG
  ! Coupling constants
  complex(dp), save :: eQED, E2_QED, gQCD, G2_QCD
  ! Weak mixing angle
  complex(dp), save :: cw, cw2, cw3, cw4, sw, sw2, sw3 ,sw4, sw6
  ! Right/left couplings of a Z boson to neutrinos, leptons, up- and down-type quarks
  complex(dp), save :: gZn(2), gZl(2), gZu(2), gZd(2)
  ! Right(1)/left(2) couplings for Higgs(H), Chi(X) = Z-Goldstone, Phi(P) = W-Goldstone
  real(dp),    save :: I3l(2)  = [0.5_dp,-0.5_dp]
  complex(dp), save :: gH(2)   = [  cONE, cONE ]
  complex(dp), save :: gX(2)   = [ -cONE, cONE ]
  complex(dp), save :: gPnl(2) = [  cONE, ZERO ]
  complex(dp), save :: gPln(2) = [  ZERO, cONE ]
  complex(dp), save :: gPud(2), gPus(2), gPub(2)
  complex(dp), save :: gPcd(2), gPcs(2), gPcb(2)
  complex(dp), save :: gPtd(2), gPts(2), gPtb(2)
  complex(dp), save :: gPdu(2), gPdc(2), gPdt(2)
  complex(dp), save :: gPsu(2), gPsc(2), gPst(2)
  complex(dp), save :: gPbu(2), gPbc(2), gPbt(2)
  complex(dp) :: gZRH, gZLH
  ! Right/left couplings ids
  integer, parameter :: ngZn = 1
  integer, parameter :: ngZl = 2
  integer, parameter :: ngZu = 3
  integer, parameter :: ngZd = 4
  integer, parameter :: ngH = 5
  integer, parameter :: ngX = 6
  integer, parameter :: ngPnl = 7
  integer, parameter :: ngPln = 8
  integer, parameter :: ngPud = 9,  ngPus = 15, ngPub = 21
  integer, parameter :: ngPcd = 10, ngPcs = 16, ngPcb = 22
  integer, parameter :: ngPtd = 11, ngPts = 17, ngPtb = 23
  integer, parameter :: ngPdu = 12, ngPdc = 18, ngPdt = 24
  integer, parameter :: ngPsu = 13, ngPsc = 19, ngPst = 25
  integer, parameter :: ngPbu = 14, ngPbc = 20, ngPbt = 26
  integer, parameter :: ngZRH = 27, ngZLH = 28

  ! Vertex scale factors for naive deviations from the Standard Model (changes do not affect CT/R2)
  real(dp), save :: lambdaHHH = 1, lambdaHHHH = 1,lambdaHWW = 1, lambdaHZZ = 1
  ! CKM Matrix, default: VCKM = diag(1,1,1)
  complex(dp), save :: VCKMdu = cONE
  complex(dp), save :: VCKMsu = ZERO
  complex(dp), save :: VCKMbu = ZERO
  complex(dp), save :: VCKMdc = ZERO
  complex(dp), save :: VCKMsc = cONE
  complex(dp), save :: VCKMbc = ZERO
  complex(dp), save :: VCKMdt = ZERO
  complex(dp), save :: VCKMst = ZERO
  complex(dp), save :: VCKMbt = cONE
  ! Coefficients of Higgs FormFactors/Pseudo-Observables
  ! Cabibbo Angle
  real(dp), save :: ThetaCabi = 0.2274_dp
  real(dp), save :: cCabi, sCabi
  ! Higgs vev

  real(dp), save :: HPOvev_unscaled  = 246.22_dp

  real(dp), save :: HPOvev
  ! Z/W-Pole
  real(dp), save :: HPOgZeL = -0.2696_dp
  real(dp), save :: HPOgZeR = 0.2315_dp
  real(dp), save :: HPOgZmL = -0.269_dp
  real(dp), save :: HPOgZmR = 0.232_dp
  real(dp), save :: HPOgZlL = -0.2693_dp
  real(dp), save :: HPOgZlR = 0.23270_dp
  real(dp), save :: HPOgZv  = 0.5_dp
  real(dp), save :: HPOgZuL = 0.3467000_dp
  real(dp), save :: HPOgZuR = -0.1547000_dp
  real(dp), save :: HPOgZdL = -0.4243000_dp
  real(dp), save :: HPOgZdR = 0.07735000_dp
  real(dp), save :: HPOgWeL = 0.994_dp
  real(dp), save :: HPOgWmL = 0.991_dp
  real(dp), save :: HPOgWlL = 1.025_dp
  real(dp), save :: HPOgWqL = 1._dp
  ! PO
  real(dp), save :: HPOkapWW = 1
  real(dp), save :: HPOkapZZ = 1
  real(dp), save :: HPOepsWW = 0
  real(dp), save :: HPOaepsWW = 0
  real(dp), save :: HPOepsZZ = 0
  real(dp), save :: HPOaepsZZ = 0
  real(dp), save :: HPOepsZA = 0
  real(dp), save :: HPOaepsZA = 0
  real(dp), save :: HPOepsAA = 0
  real(dp), save :: HPOaepsAA = 0
  complex(dp), save :: HPOepsZnn(3,2) = 0
  complex(dp), save :: HPOepsZll(3,2) = 0
  complex(dp), save :: HPOepsZdd(3,2) = 0
  complex(dp), save :: HPOepsZuu(3,2) = 0
  real(dp), save :: HPOepsWqq(3) = 0
  real(dp), save :: HPOepsWln(3) = 0
  real(dp), save :: HPOphiWeL = 0
  real(dp), save :: HPOphiWmL = 0
  real(dp), save :: HPOphiWlL = 0
  real(dp), save :: HPOphiWqL = 0
  real(dp), save :: HPOcpWeL, HPOspWeL, HPOcpWmL, HPOspWmL, HPOcpWlL, HPOspWlL, HPOcpWqL, HPOspWqL
  ! 2HDM parameters
  ! thdm_a ("alpha") is the (h0, H0) mixing angle,
  ! thdmTB is the ratio of the VEVs of the two Higgs doublets
  integer, save :: thdm_type = 2 ! 2HDM Type I or Type II

  real(dp), save :: rMA0_unscaled = 130, wMA0_unscaled = 0 ! pseudoscalar Higgs mass and width
  real(dp), save :: rMHH_unscaled = 140, wMHH_unscaled = 0 ! heavy higgs mass and width
  real(dp), save :: rMHp_unscaled = 150, wMHp_unscaled = 0 ! charged Higgs mass and width

  real(dp), save :: rMA0, wMA0, rMA02
  real(dp), save :: rMHH, wMHH, rMHH2
  real(dp), save :: rMHp, wMHp, rMHp2
  complex(dp), save :: MA0, MA02, MHH, MHH2, MHp, MHp2
  ! basic parameters: tan(beta), sin(beta-alpha), lambda5
  real(dp), save :: thdmTB = 1, thdmSBA = 1, thdmL5 = 0
  real(dp), save :: thdm_a, thdm_b
  real(dp), save :: thdmCA, thdmSA, thdmCB, thdmSB
  real(dp), save :: thdmC2A, thdmS2A, thdmC2B, thdmS2B
  real(dp), save :: thdmCAB, thdmSAB, thdmCBA
  ! Type I/II dependent couplins
  real(dp), save :: thdmYuk1, thdmYuk2, thdmYuk3
  ! Charged Higgs-fermion left/right couplings
  complex(dp), save :: thdmHpud(2), thdmHpdu(2), thdmHpcs(2), thdmHpsc(2), thdmHptb(2), thdmHpbt(2)


end module ol_parameters_decl_dp



! **********************************************************************
module ol_loop_parameters_decl_dp
! Declarations and initial values for renormalisation constants and parameters of
! dimensional regularisation, dipole subtraction, tensor-integral libraries.
! Loading the module for the first time initialises the parameters with
! their default values (given in this module). Note that when a value has
! been changed using loop_parameters_init(val=...) loading the module will not
! reset the value to its default.
! **********************************************************************
  use kind_types, only: dp
  use iso_fortran_env, only: output_unit
  use ol_parameters_decl_dp
  implicit none
  integer,        save :: loop_parameters_status = 0


  integer,        save :: maxpoint = 4, maxpoint_active = -1
  integer,        save :: maxrank = 6, maxrank_active = -1
  integer,        save :: norm_swi = 0     ! switch controlling normalisation of UV/IR poles
  character(10),  save :: norm_name
  integer,        save :: debug_ew_renorm = 0
  ! switch on UV counterterms, R2 terms, IR dipoles
  integer,        save :: SwF = 1 ! factors to multiply diagrams with fermion loops
  integer,        save :: SwB = 1 ! factors to multiply diagrams with non-fermion loops
  integer,        save :: DOI = 1 ! factors to multiply double-operator-insertions in HHEFT
  integer,        save :: CT_is_on = 1 ! switch on/off UV CT contributions
  integer,        save :: R2_is_on = 1 ! switch on/off R2 contributions
  integer,        save :: TP_is_on = 1 ! switch on/off tadpole-like contributions
  integer,        save :: IR_is_on = 1 ! 0 = off, 1 = return poles, 2 = add I operator
  logical,        save :: qedreg_on = .false. ! regularise IR singularities in QED via finite photon mass = MZ
  ! i-operator mode: 1 = QCD, 2 = EM, 0 = QCD+EM, none otherwise;
  integer,        save :: ioperator_mode = 0
  integer,        save :: polecheck_is = 0
  logical, save :: do_pole_checks = .false. ! check poles and print result when amplitude is registered

  integer,        save :: stability_mode = 11 ! 11: no trigger, default: 23
  integer,        save :: deviation_mode = 1  ! deviation measure in vamp scaling based on
                                              ! (1) k-factor (2) virtual matrix element

  real(dp), save :: trigeff_targ = .2_dp   ! target efficiency of K-factor based stability trigger (should not be << 0.1)
  real(dp), save :: abscorr_unst = 0.01_dp ! relative deviation above which a point is considered "unstable" and
                                              ! reevaluated in quad precision (if active); also logs the point in 2x modes
  real(dp), save :: ratcorr_bad = 1     ! relative deviation of two virtual matrix elements above which
                                              ! an unstable point is considered "bad" and possibly "killed"
                                              ! (i.e. the finite part of the virtual correcton is set to zero)
  real(dp), save :: ratcorr_bad_L2 = 10 ! relative deviation of two virtual matrix elements above which
                                              ! an unstable point is killed in loop induced amplitudes

  ! Collier parameters
  integer,           save :: cll_channels = 50, cll_channels_active = -1 ! number of cache channels
  real(dp),    save :: C_PV_threshold = 1.e-6 ! threshold precision to activate 3-point alternative reductions
  real(dp),    save :: D_PV_threshold = 1.e-6 ! threshold precision to activate 4-point alternative reductions
  integer,           save :: dd_red_mode    = 2     ! PV or alternative 3/4-point reductions
  integer,           save :: cll_log = 0 ! 1: create Collier log files; 2: precision monitor initmonitoring_cll()
  integer,           save :: maxcachetrain = 13 ! number of points after which the cache is fully trained
  ! setaccuracy_cll() arguments
  real(dp),    save :: cll_pvthr = 1.e-6_dp, cll_accthr = 1.e-4_dp
  real(dp),    save :: cll_mode3thr = 1.e-8_dp
  integer,           save :: cll_tenred = 7 ! settenred_cll(): # of legs from which on component reduction is used
  real(dp),    save :: ti_os_thresh = 1.e-10

  integer,           save :: olo_verbose = 0 ! OneLOop verbosity level, 0..4
  integer,           save :: olo_outunit = output_unit
  ! CutTools parameters

  real(dp),    save :: opprootsvalue_unscaled = 1000

  real(dp),    save :: opprootsvalue
  real(dp),    save :: opplimitvalue = 0.01_dp
  real(dp),    save :: oppthrs       = 1.e-6_dp
  integer,           save :: oppidig       = 0
  integer,           save :: oppscaloop    = 2

  logical,           save :: cuttools_not_init     = .true.
  logical,           save :: coli_not_init         = .true.
  logical,           save :: dd_not_init           = .true.
  logical,           save :: dd_qp_not_init        = .true.
  logical,           save :: tensorlib_not_init    = .true.
  logical,           save :: tensorlib_qp_not_init = .true.

  ! Generic parameters related to tensor reduction
  integer, save :: tensor_reduction_error = 0

  logical, save :: reset_mureg = .true.
  logical, save :: reset_olo = .true.

  integer,        save      :: nc    = 3          ! number of colours
  integer,        save      :: nf = 6, nf_up = 3, nf_down =3 ! number of quarks (total, up-type, down-type)
  integer,        save      :: nfa= 5             ! fermionic contributions to photon selfenergy
  integer,        save      :: nq_nondecoupl = 0  ! number of quarks which do not decouple above threshold,
                                                  ! i.e. always contribute to the alpha_s running
  integer,        save      :: N_lf  = 5          ! number of massless quark flavours
  integer,        save      :: N_lu  = 2          ! number of massless up-quark flavours
  integer,        save      :: N_ld  = 3          ! number of massless down-quark flavours
  integer,        save      :: N_ll  = 3          ! number of massless lepton flavours
! ifdef 1


  real(dp), save      :: ca    = 3          ! adjoint Casimir
  real(dp), save      :: cf    = 4._dp/3 ! fundamental Casimir
  real(dp), save      :: tf    = 0.5_dp  ! generator normalisation

  complex(dp), parameter      :: Qu   = 4._dp/3. ! up-type quark electrical charge squared
  complex(dp), parameter      :: Qd   = -1._dp/3. ! down-type quark electrical charge squared
  complex(dp), parameter      :: Ql   = -1 ! lepton electrical charge squared

  complex(dp), parameter      :: Qu2   = 4._dp/9. ! up-type quark electrical charge squared
  complex(dp), parameter      :: Qd2   = 1._dp/9. ! down-type quark electrical charge squared
  complex(dp), parameter      :: Ql2   = 1 ! lepton electrical charge squared

  complex(dp), parameter      :: Qu3   = 8._dp/27. ! up-type quark electrical charge cubed
  complex(dp), parameter      :: Qd3   = -1._dp/27. ! down-type quark electrical charge cubed
  complex(dp), parameter      :: Ql3   = -1 ! lepton electrical charge cubed

  complex(dp), save      :: Qf2sum   = 20._dp/3. ! sum_f ( Qf^2 ) for gamma -> FF QED splittings in I-Operator

  real(dp), save      :: polescale   = 1 ! used as pole values in VAMP2chk to determine the true poles
  real(dp), save      :: de1_UV      = 0 ! numerical value of single UV pole (independent of norm-convention)
  real(dp), save      :: de1_IR      = 0 ! numerical value of single IR pole (independent of norm-convention)
  real(dp), save      :: de2_i_IR    = 0 ! numerical value of double IR pole using actual norm-convention
  real(dp), save      :: de2_i_shift = 0 ! double pole shift defining actual norm convention

  real(dp), save      :: muren_unscaled = 100    ! renormalisation scale
  real(dp), save      :: mureg_unscaled = 100    ! regularization scale

  real(dp), save      :: muren
  real(dp), save      :: mureg
  real(dp), save      :: x_UV  = 1       ! rescaling factor for dim-reg scale in UV-divergent quantities
  real(dp), save      :: x_IR  = 1       ! rescaling factor for dim-reg scale in IR-divergent quantities
  real(dp), parameter :: kappa = 2/3._dp ! kappa parameter used in dipole subtraction

  real(dp), save :: LambdaMC2_unscaled = 0 ! squared mass MSbar renormalization scale for c quark
  real(dp), save :: LambdaMB2_unscaled = 0 ! squared mass MSbar renormalization scale for b quark
  real(dp), save :: LambdaMT2_unscaled = 0 ! squared mass MSbar renormalization scale for t quark
  real(dp), save :: LambdaYC2_unscaled = 0 ! squared yukawa MSbar renormalization scale for c quark
  real(dp), save :: LambdaYB2_unscaled = 0 ! squared yukawa MSbar renormalization scale for b quark
  real(dp), save :: LambdaYT2_unscaled = 0 ! squared yukawa MSbar renormalization scale for t quark

  real(dp), save :: LambdaMC2
  real(dp), save :: LambdaMB2
  real(dp), save :: LambdaMT2
  real(dp), save :: LambdaYC2
  real(dp), save :: LambdaYB2
  real(dp), save :: LambdaYT2


  ! the following derived parameters are initilised by subroutine loop_parameters_init
  real(dp), save :: de2_0_IR  ! numerical value of double IR pole using LH-accord convention (i=0)
  real(dp), save :: de2_1_IR  ! numerical value of double IR pole using COLI convention (i=1)
  real(dp), save :: muren2    ! squared renormalisation scale
  real(dp), save :: mureg2    ! squared regularization scale
  real(dp), save :: mu2_UV    ! dim-reg scale for UV-divergent quantities
  real(dp), save :: mu2_IR    ! dim-reg scale for IR-divergent quantities
  real(dp), save :: muyc2 ! squared yukawa renormalization scale for c quark
  real(dp), save :: muyb2 ! squared yukawa renormalization scale for b quark
  real(dp), save :: muyt2 ! squared yukawa renormalization scale for t quark

  ! the following renormalisation constants are initilised by subroutine QCD_renormalisation
  complex(dp), save :: dZMC     = 0 ! charm-quark mass RC        : MC_bare = MC*(1+dZMC)
  complex(dp), save :: dZMB     = 0 ! bottom-quark mass RC       : MB_bare = MB*(1+dZMB)
  complex(dp), save :: dZMT     = 0 ! top-quark mass RC          : MT_bare = MT*(1+dZMT)
  complex(dp), save :: dZYC     = 0 ! charm-quark yukawa RC      : YC_bare = YC*(1+dZYC)
  complex(dp), save :: dZYB     = 0 ! bottom-quark yukawa RC     : YB_bare = YB*(1+dZYB)
  complex(dp), save :: dZYT     = 0 ! top-quark yukawa RC        : YT_bare = YT*(1+dZYT)
  real(dp),    save :: dZg      = 0 ! gluon-field RC             : G_bare  = (1+dZg/2)*G_ren
  real(dp),    save :: dZq      = 0 ! massless-quark field RC    : Q_bare  = (1+dZq/2)*Q_ren
  real(dp),    save :: dZc      = 0 ! charm-quark field RC       : idem
  real(dp),    save :: dZb      = 0 ! bottom-quark field RC      : idem
  real(dp),    save :: dZt      = 0 ! top-quark field RC         : idem
  real(dp),    save :: dgQCD    = 0 ! strong coupling RC         : g_bare  = (1+dgQCD)*g_ren
  real(dp),    save :: dgQCDym  = 0 ! YM-contribution to delnG
  real(dp),    save :: dgQCDfer = 0 ! fermionic-contribution to delnG

  ! Counter terms for QCD corrections
  complex(dp), save :: ctqq(2) ! massless quark propagator counter term
  complex(dp), save :: ctcc(2) ! charm quark propagator counter term
  complex(dp), save :: ctbb(2) ! bottom quark propagator counter term
  complex(dp), save :: cttt(2) ! top quark propagator counter term
  complex(dp), save :: ctGG(3) ! gluon propagator counter term
  real(dp),    save :: ctGqq   ! massless quark-gluon vertex counter term
  real(dp),    save :: ctGcc   ! charm quark-gluon vertex counter term (massive or massless c)
  real(dp),    save :: ctGbb   ! bottom quark-gluon vertex counter term (massive or massless b)
  real(dp),    save :: ctGtt   ! top quark-gluon vertex counter term
  real(dp),    save :: ctVVV   ! three gluon vertex counter term
  real(dp),    save :: ctVVVV  ! four gluon vertex counter term (times 1/2)
  real(dp),    save :: ctVdt   ! Wdt (massless d) vertex counter term
  real(dp),    save :: ctVsc   ! Wcs (massless s) vertex counter term
  real(dp),    save :: ctVst   ! Wts (massless s) vertex counter term
  real(dp),    save :: ctVbu   ! Wub (massless u) vertex counter term
  real(dp),    save :: ctVbc   ! Wcb (massive or massless b/c) vertex counter term
  real(dp),    save :: ctVbt   ! Wtb (massive or massless b) vertex counter term
  real(dp),    save :: ctVtt   ! Att and Ztt vertex counter term
  real(dp),    save :: ctVcc   ! Acc and Zcc vertex counter term (massive or massless c)
  real(dp),    save :: ctVbb   ! Abb and Zbb vertex counter term (massive or massless b)
  real(dp),    save :: ctVqq   ! Aqq and Zqq (massless q) vertex counter term
  complex(dp), save :: ctSud(2), ctSus(2), ctSub(2)
  complex(dp), save :: ctScd(2), ctScs(2), ctScb(2)
  complex(dp), save :: ctStd(2), ctSts(2), ctStb(2)
  complex(dp), save :: ctSdu(2), ctSdc(2), ctSdt(2)
  complex(dp), save :: ctSsu(2), ctSsc(2), ctSst(2)
  complex(dp), save :: ctSbu(2), ctSbc(2), ctSbt(2)
  complex(dp), save :: ctSqq
  complex(dp), save :: ctScc
  complex(dp), save :: ctSbb
  complex(dp), save :: ctStt

  ! Additional parameters for R2
  complex(dp), save :: MQ2sum, MQ2sumpairs
  complex(dp), save :: YQD2sum, YQU2sum
  complex(dp), save :: YQD2sumpairs, YQU2sumpairs

  ! Additional counterterms for R2 QCD
  complex(dp), save :: ctZGG
  complex(dp), save :: ctHGG
  complex(dp), save :: ctAAGG
  complex(dp), save :: ctAZGG
  complex(dp), save :: ctZZGG
  complex(dp), save :: ctWWGG
  complex(dp), save :: ctHHGG
  complex(dp), save :: ctHXGG
  complex(dp), save :: ctXXGG
  complex(dp), save :: ctPPGG
  complex(dp), save :: ctAGGG(2)
  complex(dp), save :: ctZGGG(2)
  integer,           save :: R2GGGG

  ! Counterterms for HEFT
  real(dp), save :: ctHEFTggh(5)
  real(dp), save :: ctHEFTgggh
  real(dp), save :: ctHEFTggggh
  real(dp), save :: R2HEFTggggh
  real(dp), save :: R2HEFThqq
  real(dp), save :: R2HEFTghqq

  ! 2HDM
  complex(dp), save :: thdmctHpsc(2), thdmctHpbt(2), thdmctHpcs(2), thdmctHptb(2)
  complex(dp), save :: thdmctHGG, thdmctHh0GG, thdmctHHGG, thdmctHHhGG
  complex(dp), save :: thdmctXA0GG, thdmctHhHhGG, thdmctA0A0GG
  complex(dp), save :: thdmctPHpGG, thdmctHpHpGG

  ! EW_renormalisation renormalisation constants
  complex(dp), save :: dZMBEW     = 0 ! bottom-quark mass RC       : MB_bare = MB+dZMBEW
  complex(dp), save :: dZMTEW     = 0 ! top-quark mass RC          : MT_bare = MT+dZMTEW
  complex(dp), save :: dZMLEW     = 0 ! tau-lepton mass RC         : ML_bare = ML+dZMLEW
  complex(dp), save :: dZMEEW     = 0 ! electron mass RC           : ME_bare = ME+dZMLEW
  complex(dp), save :: dZMMEW     = 0 ! muon mass RC               : MM_bare = MM+dZMLEW
  complex(dp), save :: dZMW2EW    = 0 ! W mass RC                  : MW^2_bare = MW^2+dZMW2EW^2
  complex(dp), save :: dZMZ2EW    = 0 ! Z mass RC                  : MZ^2_bare = MZ^2+dZMZ2EW^2
  complex(dp), save :: dZMH2EW    = 0 ! H mass RC                  : MH^2_bare = MH^2+dZMH2EW^2
  complex(dp), save :: dswEW      = 0 ! sin EW mixing angle RC         : sw_bare = sw + dswEW  i.e. dswEW/swEW = - c^2/s^2 dcwEW/c ; dcEW/c=1/2(dZMW2/MW^2-dZMZ2/MZ^2)
  complex(dp), save :: dcwEW      = 0 ! cos EW mixing angle RC         : dcwEW/c=1/2(dZMW2/MW^2-dZMZ2/MZ^2) defined for convinience

!   complex(rp), save      :: dZqLEW     = 0 ! L-massless-quark field RC : Q_bare  = (1+1/2*dZqLEW)*Q_ren
!   complex(rp), save      :: dZqREW     = 0 ! R-massless-quark field RC : Q_bare  = (1+1/2*dZqREW)*Q_ren
  complex(dp), save :: dZuLEW     = 0 ! L-massless-u-quark field RC : Q_bare  = (1+1/2*dZuLEW)*Q_ren
  complex(dp), save :: dZuREW     = 0 ! R-massless-u-quark field RC : Q_bare  = (1+1/2*dZuREW)*Q_ren
  complex(dp), save :: dZdLEW     = 0 ! L-massless-d-quark field RC : Q_bare  = (1+1/2*dZdLEW)*Q_ren
  complex(dp), save :: dZdREW     = 0 ! R-massless-d-quark field RC : Q_bare  = (1+1/2*dZdREW)*Q_ren
  complex(dp), save :: dZbLEW     = 0 ! L-bottom-quark field RC     : idem
  complex(dp), save :: dZbREW     = 0 ! R-bottom-quark field RC     : idem
  complex(dp), save :: dZtLEW     = 0 ! L-top-quark field RC        : idem
  complex(dp), save :: dZtREW     = 0 ! R-top-quark field RC        : idem
  complex(dp), save :: dZeLEW     = 0 ! L-electron field RC         : idem
  complex(dp), save :: dZeREW     = 0 ! R-electron field RC         : idem
  complex(dp), save :: dZmuLEW    = 0 ! L-muon field RC             : idem
  complex(dp), save :: dZmuREW    = 0 ! R-muon field RC             : idem
  complex(dp), save :: dZnLEW     = 0 ! L-neutrino field RC         : idem
  complex(dp), save :: dZnlLEW    = 0 ! L-tau-neutrino field RC     : idem
  complex(dp), save :: dZlLEW     = 0 ! L-tau-lepton field RC       : idem
  complex(dp), save :: dZlREW     = 0 ! R-tau-lepton field RC       : idem
  complex(dp), save :: dZuLEWcc   = 0 ! L-massless-u-quark field RC : Q_bare  = (1+1/2*dZuLEW)*Q_ren
  complex(dp), save :: dZuREWcc   = 0 ! R-massless-u-quark field RC : Q_bare  = (1+1/2*dZuREW)*Q_ren
  complex(dp), save :: dZdLEWcc   = 0 ! L-massless-d-quark field RC : Q_bare  = (1+1/2*dZdLEW)*Q_ren
  complex(dp), save :: dZdREWcc   = 0 ! R-massless-d-quark field RC : Q_bare  = (1+1/2*dZdREW)*Q_ren
  complex(dp), save :: dZbLEWcc   = 0 ! L-bottom-quark field RC     : idem
  complex(dp), save :: dZbREWcc   = 0 ! R-bottom-quark field RC     : idem
  complex(dp), save :: dZtLEWcc   = 0 ! L-top-quark field RC        : idem
  complex(dp), save :: dZtREWcc   = 0 ! R-top-quark field RC        : idem
  complex(dp), save :: dZeLEWcc   = 0 ! L-electron field RC         : idem
  complex(dp), save :: dZeREWcc   = 0 ! R-electron field RC         : idem
  complex(dp), save :: dZmuLEWcc  = 0 ! L-muon field RC             : idem
  complex(dp), save :: dZmuREWcc  = 0 ! R-muon field RC             : idem
  complex(dp), save :: dZnLEWcc   = 0 ! L-neutrino field RC         : idem
  complex(dp), save :: dZnlLEWcc  = 0 ! L-tau-neutrino field RC     : idem
  complex(dp), save :: dZlLEWcc   = 0 ! L-tau-lepton field RC       : idem
  complex(dp), save :: dZlREWcc    = 0 ! R-tau-lepton field RC       : idem

  complex(dp), save :: dZWEW      = 0 ! W field RC                 : idem
  complex(dp), save :: dZZZEW     = 0 ! ZZ field RC                : idem
  complex(dp), save :: dZAZEW     = 0 ! AZ field RC                : idem
  complex(dp), save :: dZZAEW     = 0 ! AZ field RC                : idem
  complex(dp), save :: dZAAEW     = 0 ! AA field RC                : idem
  complex(dp), save :: dZAAEWdreg = 0 ! AA field RC dim-reg      : Photon field renormalization constant with ligh fermion contributions in dim-reg
  complex(dp), save :: dZAAEWnreg = 0 ! AA field RC dim-reg      : Photon field renormalization constant with ligh fermion contributions in n-reg
  complex(dp), save :: dZHEW      = 0 ! H field RC                 : idem

  complex(dp), save :: dtEW       = 0 ! tadpole-RC                 :
  complex(dp), save :: dZeQEDEW   = 0 ! EW coupling RC         : e_bare  = (1+dZeEW)*e_ren
  complex(dp), save :: dZe0QEDEW  = 0 ! EW coupling RC         : e_bare  = (1+dZeEW)*e_ren in on-shell/a(0) scheme
  complex(dp), save :: dZe0QEDEWdreg = 0 ! EW coupling RC   : e_bare  = (1+dZeEW)*e_ren in on-shell/a(0) scheme with light fermuion controbutions in dim-reg
  complex(dp), save :: dZe0QEDEWnreg = 0 ! EW coupling RC   : e_bare  = (1+dZeEW)*e_ren in on-shell/a(0) scheme with light fermuion controbutions in n-reg
  complex(dp), save :: dZeGmuQEDEW   = 0 ! EW coupling RC   : e_bare  = (1+dZeEW)*e_ren in Gmu-scheme
  complex(dp), save :: dZeZQEDEW     = 0 ! EW coupling RC   : e_bare  = (1+dZeEW)*e_ren in a(MZ)-scheme
  complex(dp), save :: dVLEW         = 0 ! L-lepton-Z vertex RC
  complex(dp), save :: dVREW         = 0 ! R-lepton-Z vertex RC

  ! Counter terms for EW corrections
  ! VV Vector propagators
  complex(dp), save :: EWctWW(3)      = 0
  complex(dp), save :: EWctZZ(3)      = 0
  complex(dp), save :: EWctAZ(3)      = 0
  complex(dp), save :: EWctAA(3)      = 0
  ! SS scalar propagators
  complex(dp), save :: EWctHH(2)      = 0
  complex(dp), save :: EWctXX(2)      = 0
  complex(dp), save :: EWctPP(2)      = 0
  ! SV scalar-vector mixing
  complex(dp), save :: EWctXA      = 0
  complex(dp), save :: EWctXZ      = 0
  complex(dp), save :: EWctPW      = 0
  ! FF fermionic propagators
  complex(dp), save :: EWctuu(4)      = 0
  complex(dp), save :: EWctdd(4)      = 0
  complex(dp), save :: EWcttt(4)      = 0
  complex(dp), save :: EWctbb(4)      = 0
  complex(dp), save :: EWctee(4)      = 0
  complex(dp), save :: EWctmm(4)      = 0
  complex(dp), save :: EWctLL(4)      = 0
  complex(dp), save :: EWctnn(4)      = 0
  complex(dp), save :: EWctnlnl(4)      = 0
  !VVVV
  complex(dp), save :: EWctWWWW(2)      = 0
  complex(dp), save :: EWctWWZZ(2)      = 0
  complex(dp), save :: EWctWWAZ(2)      = 0
  complex(dp), save :: EWctWWAA(2)      = 0
  !VVVV pure R2
  complex(dp), save :: EWctR2AAAA      = 0
  complex(dp), save :: EWctR2AAAZ      = 0
  complex(dp), save :: EWctR2AAZZ      = 0
  complex(dp), save :: EWctR2AZZZ      = 0
  complex(dp), save :: EWctR2ZZZZ      = 0
  !VVV
  complex(dp), save :: EWctAWW      = 0
  complex(dp), save :: EWctZWW      = 0
  !SSSS
  complex(dp), save :: EWctSSSS1      = 0
  complex(dp), save :: EWctSSSS2      = 0
  complex(dp), save :: EWctSSSS3      = 0
  complex(dp), save :: EWctHHHH      = 0
  complex(dp), save :: EWctHHXX      = 0
  complex(dp), save :: EWctHHPP      = 0
  complex(dp), save :: EWctXXXX      = 0
  complex(dp), save :: EWctXXPP      = 0
  complex(dp), save :: EWctPPPP      = 0
  !SSS
  complex(dp), save :: EWctHHH      = 0
  complex(dp), save :: EWctHXX      = 0
  complex(dp), save :: EWctHPP      = 0
  !VVSS
  complex(dp), save :: EWctWWXX      = 0
  complex(dp), save :: EWctWWHH      = 0
  complex(dp), save :: EWctWWPP      = 0
  complex(dp), save :: EWctZZPP      = 0
  complex(dp), save :: EWctZAPP      = 0
  complex(dp), save :: EWctAAPP      = 0
  complex(dp), save :: EWctZZHH      = 0
  complex(dp), save :: EWctZZXX      = 0
  complex(dp), save :: EWctZAHH      = 0
  complex(dp), save :: EWctZAXX      = 0
  complex(dp), save :: EWctWZPH      = 0
  complex(dp), save :: EWctWAPH      = 0
  complex(dp), save :: EWctWZPX      = 0
  complex(dp), save :: EWctWAPX      = 0
  !VVSS R2
  complex(dp), save :: EWctAAHH      = 0
  complex(dp), save :: EWctAAXX      = 0
  !VSS
  complex(dp), save :: EWctAXH      = 0
  complex(dp), save :: EWctZXH      = 0
  complex(dp), save :: EWctAPP      = 0
  complex(dp), save :: EWctZPP      = 0
  complex(dp), save :: EWctWPH      = 0
  complex(dp), save :: EWctWPX      = 0
  !SVV
  complex(dp), save :: EWctHWW      = 0
  complex(dp), save :: EWctHZZ      = 0
  complex(dp), save :: EWctHZA      = 0
  complex(dp), save :: EWctPWZ      = 0
  complex(dp), save :: EWctPWA      = 0
  ! pure R2 SVV
  complex(dp), save :: EWctHAA      = 0
  !VFF
  ! Aff
  complex(dp), save :: EWctAuu(2)      = 0
  complex(dp), save :: EWctAdd(2)      = 0
  complex(dp), save :: EWctAtt(2)      = 0
  complex(dp), save :: EWctAbb(2)      = 0
  complex(dp), save :: EWctAee(2)      = 0
  complex(dp), save :: EWctAmm(2)      = 0
  complex(dp), save :: EWctALL(2)      = 0
  complex(dp), save :: EWctAnn(2)      = 0
  ! Zff
  complex(dp), save :: dgZu(2)      = 0
  complex(dp), save :: dgZd(2)      = 0
  complex(dp), save :: dgZl(2)      = 0
  complex(dp), save :: dgZn(2)      = 0
  complex(dp), save :: EWctVuu(2)      = 0
  complex(dp), save :: EWctVdd(2)      = 0
  complex(dp), save :: EWctVtt(2)      = 0
  complex(dp), save :: EWctVbb(2)      = 0
  complex(dp), save :: EWctVee(2)      = 0
  complex(dp), save :: EWctVmm(2)      = 0
  complex(dp), save :: EWctVLL(2)      = 0
  complex(dp), save :: EWctVnn(2)      = 0
  complex(dp), save :: EWctVnlnl(2)      = 0
  ! Wff
  complex(dp), save :: EWctVdu      = 0
  complex(dp), save :: EWctVbt      = 0
  complex(dp), save :: EWctVen      = 0
  complex(dp), save :: EWctVLn      = 0
  complex(dp), save :: EWctVud      = 0
  complex(dp), save :: EWctVtb      = 0
  complex(dp), save :: EWctVne      = 0
  complex(dp), save :: EWctVnL      = 0
  ! Gff mixed EW/QCD
  complex(dp), save ::  EWctGuu(2)      = 0
  complex(dp), save ::  EWctGdd(2)      = 0
  complex(dp), save ::  EWctGtt(2)      = 0
  complex(dp), save ::  EWctGbb(2)      = 0
  !SFF
  complex(dp), save :: EWctHee(2)      = 0
  complex(dp), save :: EWctHmm(2)      = 0
  complex(dp), save :: EWctHLL(2)      = 0
  complex(dp), save :: EWctHtt(2)      = 0
  complex(dp), save :: EWctHbb(2)      = 0
  complex(dp), save :: EWctXee(2)      = 0
  complex(dp), save :: EWctXmm(2)      = 0
  complex(dp), save :: EWctXtt(2)      = 0
  complex(dp), save :: EWctXbb(2)      = 0
  complex(dp), save :: EWctXLL(2)      = 0
  complex(dp), save :: EWctPud(2)      = 0
  complex(dp), save :: EWctPdu(2)      = 0
  complex(dp), save :: EWctPtb(2)      = 0
  complex(dp), save :: EWctPbt(2)      = 0
  complex(dp), save :: EWctPne(2)      = 0
  complex(dp), save :: EWctPen(2)      = 0
  complex(dp), save :: EWctPnL(2)      = 0
  complex(dp), save :: EWctPLn(2)      = 0
  ! VUU
  complex(dp), save :: EWctAUWUW      = 0
  complex(dp), save :: EWctZUWUW      = 0
  complex(dp), save :: EWctWUWUZ      = 0
  complex(dp), save :: EWctWUZUW      = 0
  complex(dp), save :: EWctWUWUA      = 0
  complex(dp), save :: EWctWUAUW      = 0
  ! SUU
  complex(dp), save :: EWctHUZUZ      = 0
  complex(dp), save :: EWctHUWUW      = 0
  complex(dp), save :: EWctXUWUW      = 0
  complex(dp), save :: EWctPUZUW      = 0
  complex(dp), save :: EWctPUWUZ      = 0
  complex(dp), save :: EWctPUWUA      = 0

  ! Additional parameters for R2 EW
  complex(dp), save :: sumMQ2
  complex(dp), save :: sumMQ2Q2
  complex(dp), save :: sumMQ2QI
  complex(dp), save :: sumMQ4
  complex(dp), save :: sumMUD
  complex(dp), save :: sumMUD2
  complex(dp), save :: sumMU2
  complex(dp), save :: sumMD2
  complex(dp), save :: sumMU4
  complex(dp), save :: sumMD4
  complex(dp), save :: sumMQ2QUD
  complex(dp), save :: sumML2
  complex(dp), save :: sumML2Q2
  complex(dp), save :: sumML2QI
  complex(dp), save :: sumML4

end module ol_loop_parameters_decl_dp

