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

module ol_data_types_qp
  use kind_types, only: qp, qp, intkind1, intkind2

  type wfun
    ! four complex components for the wave function
    complex(qp) :: j(4)
    complex(qp), pointer :: j_prev(:)
    ! indicator if left- or right components of of-shell line vanish
    !                             j= (0,0,0,0) (0,0,j3,j4) (j1,j2,0,0) (j1,j2,j3,j4)
    integer(intkind1) :: h      ! B"00"      B"01"       B"10"        B"11"
    integer(intkind2) :: e      ! helicities of external on-shell lines
    integer(intkind2) :: t      ! label for the external subtree
    integer(intkind2) :: n_part ! number of particles in the subtree
    integer(intkind2) :: hf     ! global base-4 helicity label
  end type wfun

  type polcont
    complex(qp) :: j
    integer(intkind2) :: e ! helicities of external on-shell lines
    integer(intkind2) :: s ! table for final helicity syncronisation
  end type polcont

  ! open-loop derived data type used in the optimezed helicity summation
  type hol
    ! OpenLoops Coefficient: (alpha,rank,beta,helicity_state)
    complex(qp), dimension(:,:,:,:), allocatable :: j

    ! Helicity configurations array
    integer(intkind2), dimension(:)      , allocatable :: hf
    integer :: mode = 1
    real(qp) :: error
    integer :: npoint = 0
    integer :: ndrs = 0
    integer :: nred = 0

  end type hol

  ! derived tensor data type for closed-loop
  type hcl
    complex(qp) , dimension(:), allocatable :: cmp

    integer :: mode = 1
    real(qp) :: error
    integer :: ndrs = 0
    integer :: nred = 0

  end type hcl

  type met
    real(qp) :: cmp

    integer :: mode = 1
    real(qp) :: error
    integer :: sicount = 0
    integer :: ndrs = 0
    integer :: nred = 0

  end type met

  ! equivalent to polcont, with the addition of an extra hf label for
  ! the global helicity state
  type Hpolcont
    complex(qp)  :: j
    integer(intkind2)  :: e  ! helicities of external on-shell lines
    integer(intkind2)  :: hf ! global base-4 helicity label
    integer(intkind2)  :: s  ! table for final helicity syncronisation
  end type Hpolcont

  !! l_i basis for the on-the-fly reduction
  type basis
    complex(qp) :: vect1(4)     !! l_{1,\mu} in light-cone rep
    complex(qp) :: vect2(4)     !! l_{2,\mu} in light-cone rep
    complex(qp) :: vect3(4)     !! l_{3,\mu} in light-cone rep
    complex(qp) :: vect4(4)     !! l_{4,\mu} in light-cone rep
    complex(qp) :: tens1(10)
    complex(qp) :: tens2(10)
    complex(qp) :: tens3(10)
    complex(qp) :: tens4(10,4)
    complex(qp) :: tens5(10,4)
    complex(qp) :: gamma
    complex(qp) :: alpha(2)
    integer           :: mom1
    integer           :: mom2
    complex(qp) :: li(4,4)
  end type basis

  !! Set containing the reduction basis
  type redset4
    type(basis) :: redbasis
    complex(qp) :: p3scalars(0:4)
    integer :: perm(3)
    integer :: mom3
    real(qp) :: gd2,gd3

  end type redset4

  type redset5
    type(basis) :: redbasis
    complex(qp) :: p3scalars(0:4)
    integer :: perm(4)
    integer :: mom3
    integer :: mom4
    real(qp) :: gd2,gd3
  end type redset5

  !! Scalar Box
  type scalarbox
    complex(qp) :: poles(0:2)         ! finite, eps^(-1), eps^(-2)
    complex(qp) :: onshell_cuts(2,5)  ! on-shell cuts of the box
    real(qp) :: cut_error
    real(qp) :: box_error

  end type scalarbox





end module ol_data_types_qp



module ol_momenta_decl_qp
  use kind_types, only: qp, qp
  use ol_debug, only: ol_msg
  implicit none
  ! Internal momenta for external particles
  ! Components 1:4 = light cone representation; component 5 = squared momentum
  complex(qp), allocatable, save :: Q(:,:) ! Q(5,0:2**Nparticle-1)
  complex(qp), allocatable, save :: QInvariantsMatrix(:,:) ! QInvariantsMatrix(Nparticle,Nparticle)
  ! Components 1:4 = light cone representation; component 5 = squared masses, component 6 = scalar products
  complex(qp), allocatable, save :: L(:,:) ! L(6,0:2**Nparticle-1)

  complex(qp), allocatable, save :: Q_qp(:,:) ! Q(5,0:2**Nparticle-1)
  complex(qp), allocatable, save :: QInvariantsMatrix_qp(:,:) ! QInvariantsMatrix(Nparticle,Nparticle)
  complex(qp), allocatable, save :: L_qp(:,:) ! L(6,0:2**Nparticle-1)

  integer, save :: collconf = 0
  integer, save :: softconf = 0

  contains

  function momenta_nan_check(P)
    implicit none
    real(qp), intent(in) :: P(:,:)
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

end module ol_momenta_decl_qp


module ol_external_decl_qp
  use kind_types, only: qp
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
  real(qp), allocatable, save :: P_ex(:,:) ! uncleaned external 2->n-2 momenta, set by conv_mom_scatt2in
  integer, allocatable, save :: M_ex(:) ! external 2->n-2 mass ids, set by conv_mom_scatt2in

end module ol_external_decl_qp



module ol_pseudotree_qp
  use kind_types, only: qp
  implicit none
  ! loop momentum in pseudo tree (standard representation)
  real(qp), save :: pseudotree_momentum(0:3) = [ 314.1592653589793_qp, 271.8281828459045_qp, 100._qp, 57.72156649015328_qp ]
  ! Wave functions for pseudo tree
  complex(qp), save :: exloop(4,2) = reshape([ 2.718_qp, 3.141_qp,  0.9159_qp, 1._qp,  &
                                               1._qp,    0.5772_qp, 1.618_qp,  1.282_qp], [ 4, 2 ])
end module ol_pseudotree_qp



module ol_tensor_storage_qp
  use kind_types, only: qp
  implicit none
  complex(qp), allocatable, save :: tensor_stored(:)
  integer, save :: rank_stored
  integer, save :: array_length_stored ! length of the array associated with rank_stored
  integer, save :: tensor_storage_maxrank = -1
end module ol_tensor_storage_qp



module ol_parameters_decl_qp
  ! Declarations and initial values for numerical and physical parameters like masses, widths, and couplings.
  ! Loading the module for the first time initialises the parameters with their default values (given in this module).
  ! Note that when a value has been changed using parameters_init(val=...) loading the module
  ! will not reset the value to its default.
  use kind_types, only: qp
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
  real(qp), save :: hp_loopacc = 8._qp  ! loop accuracy target

  ! internal, do not touch change
  real(qp), save :: hp_err_thres = 8._qp  ! accumulated error threshold
  real(qp), save :: hp_step_thres = 0.5_qp ! step threshold
  real(qp), save :: hp_redset_gd3_thres = 3e-4_qp

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
  real(qp), save :: hp_max_err = 0._qp

  real(qp), save :: max_error = 0._qp

  !  compute real bubble diagrams in counterterms.

  !  use some hacks to stabilize ir events




  ! Numerical constants
  real(qp),    parameter :: rONE   = 1
  real(qp),    parameter :: rZERO  = 0
  real(qp),    parameter :: rZERO2 = 0
  real(qp),    parameter :: pi     = acos(-1._qp)
  real(qp),    parameter :: pi2_6  = (pi**2)/6
  real(qp),    parameter :: sqrt2  = sqrt(2._qp)
  real(qp),    parameter :: sqrt05 = sqrt(0.5_qp)
  complex(qp), parameter :: cONE   = 1
  complex(qp), parameter :: ZERO   = 0
  complex(qp), parameter :: ZERO2  = 0
  complex(qp), parameter :: CI     = (0._qp, 1._qp)
  complex(qp) :: integralnorm = CI/(16*pi**2)
  complex(qp) :: countertermnorm = 1._qp/(16._qp*pi**2)

  ! scale factor for dimensionful parameters
  real(qp), save :: scalefactor = 1._qp
  logical,        save :: reset_scalefactor = .false.
  integer,        save :: scaling_mode = 1 ! 1: reduction only, 3: everything

  ! synchronise Yukawa masses with masses
  logical, save :: yuk_from_mass = .true.

  real(qp), save :: alpha_QCD = 0.1258086856923967_qp ! LO MRST
  real(qp), save :: alpha_QED_MZ = 1/128._qp          ! alpha(MZ) derived from PDG 2014
  real(qp), save :: alpha_QED_0  = 1/137.035999074_qp  ! alpha(0) from PDG 2014
  real(qp), save :: alpha_QED, alpha_QED_input
  real(qp), save :: alpha_QED_Gmu
  real(qp), save :: sw2_input = 0.222626515643872389077366863865260666_qp

  real(qp), save :: Gmu
  ! Everything beyond this line is derived from the values given above and initialised by parameters_init().
  real(qp), save :: rescalefactor = 1.1_qp
  ! scaled masses, widths and yukawas
  real(qp), save :: rME, wME, rYE, wYE
  real(qp), save :: rMM, wMM, rYM, wYM
  real(qp), save :: rML, wML, rYL, wYL
  real(qp), save :: rMU, wMU, rYU, wYU
  real(qp), save :: rMD, wMD, rYD, wYD
  real(qp), save :: rMS, wMS, rYS, wYS
  real(qp), save :: rMC, wMC, rYC, wYC
  real(qp), save :: rMB, wMB, rYB, wYB
  real(qp), save :: rMT, wMT, rYT, wYT
  real(qp), save :: rMW, wMW
  real(qp), save :: rMZ, wMZ
  real(qp), save :: rMH, wMH
  real(qp), save :: rMX, wMX
  real(qp), save :: rMY, wMY
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
  complex(qp), save ::  ME,   MM,   ML,   MU,   MD,   MS,   MC,   MB,   MT,   MW,   MZ,   MH,   MX,   MY
  complex(qp), save ::  ME2,  MM2,  ML2,  MU2,  MD2,  MS2,  MC2,  MB2,  MT2,  MW2,  MZ2,  MH2,  MX2,  MY2
  complex(qp), save ::  YE,   YM,   YL,   YU,   YD,   YS,   YC,   YB,   YT
  complex(qp), save ::  YE2,  YM2,  YL2,  YU2,  YD2,  YS2,  YC2,  YB2,  YT2
  real(qp),    save :: rYE2, rYM2, rYL2, rYU2, rYD2, rYS2, rYC2, rYB2, rYT2
  real(qp),    save :: rME2, rMM2, rML2, rMU2, rMD2, rMS2, rMC2, rMB2, rMT2, rMW2, rMZ2, rMH2, rMX2, rMY2
  complex(qp), save :: YC2pair, YB2pair, YT2pair ! pair masses: only non-zero if the SU(2) partner is active
  ! collinear mass regulator for photon WF CT
  real(qp),    save :: MREG
  ! Coupling constants
  complex(qp), save :: eQED, E2_QED, gQCD, G2_QCD
  ! Weak mixing angle
  complex(qp), save :: cw, cw2, cw3, cw4, sw, sw2, sw3 ,sw4, sw6
  ! Right/left couplings of a Z boson to neutrinos, leptons, up- and down-type quarks
  complex(qp), save :: gZn(2), gZl(2), gZu(2), gZd(2)
  ! Right(1)/left(2) couplings for Higgs(H), Chi(X) = Z-Goldstone, Phi(P) = W-Goldstone
  real(qp),    save :: I3l(2)  = [0.5_qp,-0.5_qp]
  complex(qp), save :: gH(2)   = [  cONE, cONE ]
  complex(qp), save :: gX(2)   = [ -cONE, cONE ]
  complex(qp), save :: gPnl(2) = [  cONE, ZERO ]
  complex(qp), save :: gPln(2) = [  ZERO, cONE ]
  complex(qp), save :: gPud(2), gPus(2), gPub(2)
  complex(qp), save :: gPcd(2), gPcs(2), gPcb(2)
  complex(qp), save :: gPtd(2), gPts(2), gPtb(2)
  complex(qp), save :: gPdu(2), gPdc(2), gPdt(2)
  complex(qp), save :: gPsu(2), gPsc(2), gPst(2)
  complex(qp), save :: gPbu(2), gPbc(2), gPbt(2)
  complex(qp) :: gZRH, gZLH
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
  real(qp), save :: lambdaHHH = 1, lambdaHHHH = 1,lambdaHWW = 1, lambdaHZZ = 1
  ! CKM Matrix, default: VCKM = diag(1,1,1)
  complex(qp), save :: VCKMdu = cONE
  complex(qp), save :: VCKMsu = ZERO
  complex(qp), save :: VCKMbu = ZERO
  complex(qp), save :: VCKMdc = ZERO
  complex(qp), save :: VCKMsc = cONE
  complex(qp), save :: VCKMbc = ZERO
  complex(qp), save :: VCKMdt = ZERO
  complex(qp), save :: VCKMst = ZERO
  complex(qp), save :: VCKMbt = cONE
  ! Coefficients of Higgs FormFactors/Pseudo-Observables
  ! Cabibbo Angle
  real(qp), save :: ThetaCabi = 0.2274_qp
  real(qp), save :: cCabi, sCabi
  ! Higgs vev

  real(qp), save :: HPOvev
  ! Z/W-Pole
  real(qp), save :: HPOgZeL = -0.2696_qp
  real(qp), save :: HPOgZeR = 0.2315_qp
  real(qp), save :: HPOgZmL = -0.269_qp
  real(qp), save :: HPOgZmR = 0.232_qp
  real(qp), save :: HPOgZlL = -0.2693_qp
  real(qp), save :: HPOgZlR = 0.23270_qp
  real(qp), save :: HPOgZv  = 0.5_qp
  real(qp), save :: HPOgZuL = 0.3467000_qp
  real(qp), save :: HPOgZuR = -0.1547000_qp
  real(qp), save :: HPOgZdL = -0.4243000_qp
  real(qp), save :: HPOgZdR = 0.07735000_qp
  real(qp), save :: HPOgWeL = 0.994_qp
  real(qp), save :: HPOgWmL = 0.991_qp
  real(qp), save :: HPOgWlL = 1.025_qp
  real(qp), save :: HPOgWqL = 1._qp
  ! PO
  real(qp), save :: HPOkapWW = 1
  real(qp), save :: HPOkapZZ = 1
  real(qp), save :: HPOepsWW = 0
  real(qp), save :: HPOaepsWW = 0
  real(qp), save :: HPOepsZZ = 0
  real(qp), save :: HPOaepsZZ = 0
  real(qp), save :: HPOepsZA = 0
  real(qp), save :: HPOaepsZA = 0
  real(qp), save :: HPOepsAA = 0
  real(qp), save :: HPOaepsAA = 0
  complex(qp), save :: HPOepsZnn(3,2) = 0
  complex(qp), save :: HPOepsZll(3,2) = 0
  complex(qp), save :: HPOepsZdd(3,2) = 0
  complex(qp), save :: HPOepsZuu(3,2) = 0
  real(qp), save :: HPOepsWqq(3) = 0
  real(qp), save :: HPOepsWln(3) = 0
  real(qp), save :: HPOphiWeL = 0
  real(qp), save :: HPOphiWmL = 0
  real(qp), save :: HPOphiWlL = 0
  real(qp), save :: HPOphiWqL = 0
  real(qp), save :: HPOcpWeL, HPOspWeL, HPOcpWmL, HPOspWmL, HPOcpWlL, HPOspWlL, HPOcpWqL, HPOspWqL
  ! 2HDM parameters
  ! thdm_a ("alpha") is the (h0, H0) mixing angle,
  ! thdmTB is the ratio of the VEVs of the two Higgs doublets
  integer, save :: thdm_type = 2 ! 2HDM Type I or Type II

  real(qp), save :: rMA0, wMA0, rMA02
  real(qp), save :: rMHH, wMHH, rMHH2
  real(qp), save :: rMHp, wMHp, rMHp2
  complex(qp), save :: MA0, MA02, MHH, MHH2, MHp, MHp2
  ! basic parameters: tan(beta), sin(beta-alpha), lambda5
  real(qp), save :: thdmTB = 1, thdmSBA = 1, thdmL5 = 0
  real(qp), save :: thdm_a, thdm_b
  real(qp), save :: thdmCA, thdmSA, thdmCB, thdmSB
  real(qp), save :: thdmC2A, thdmS2A, thdmC2B, thdmS2B
  real(qp), save :: thdmCAB, thdmSAB, thdmCBA
  ! Type I/II dependent couplins
  real(qp), save :: thdmYuk1, thdmYuk2, thdmYuk3
  ! Charged Higgs-fermion left/right couplings
  complex(qp), save :: thdmHpud(2), thdmHpdu(2), thdmHpcs(2), thdmHpsc(2), thdmHptb(2), thdmHpbt(2)


end module ol_parameters_decl_qp



! **********************************************************************
module ol_loop_parameters_decl_qp
! Declarations and initial values for renormalisation constants and parameters of
! dimensional regularisation, dipole subtraction, tensor-integral libraries.
! Loading the module for the first time initialises the parameters with
! their default values (given in this module). Note that when a value has
! been changed using loop_parameters_init(val=...) loading the module will not
! reset the value to its default.
! **********************************************************************
  use kind_types, only: qp
  use iso_fortran_env, only: output_unit
  use ol_parameters_decl_qp
  implicit none
  integer,        save :: loop_parameters_status = 0



  real(qp), save      :: ca    = 3          ! adjoint Casimir
  real(qp), save      :: cf    = 4._qp/3 ! fundamental Casimir
  real(qp), save      :: tf    = 0.5_qp  ! generator normalisation

  complex(qp), parameter      :: Qu   = 4._qp/3. ! up-type quark electrical charge squared
  complex(qp), parameter      :: Qd   = -1._qp/3. ! down-type quark electrical charge squared
  complex(qp), parameter      :: Ql   = -1 ! lepton electrical charge squared

  complex(qp), parameter      :: Qu2   = 4._qp/9. ! up-type quark electrical charge squared
  complex(qp), parameter      :: Qd2   = 1._qp/9. ! down-type quark electrical charge squared
  complex(qp), parameter      :: Ql2   = 1 ! lepton electrical charge squared

  complex(qp), parameter      :: Qu3   = 8._qp/27. ! up-type quark electrical charge cubed
  complex(qp), parameter      :: Qd3   = -1._qp/27. ! down-type quark electrical charge cubed
  complex(qp), parameter      :: Ql3   = -1 ! lepton electrical charge cubed

  complex(qp), save      :: Qf2sum   = 20._qp/3. ! sum_f ( Qf^2 ) for gamma -> FF QED splittings in I-Operator

  real(qp), save      :: polescale   = 1 ! used as pole values in VAMP2chk to determine the true poles
  real(qp), save      :: de1_UV      = 0 ! numerical value of single UV pole (independent of norm-convention)
  real(qp), save      :: de1_IR      = 0 ! numerical value of single IR pole (independent of norm-convention)
  real(qp), save      :: de2_i_IR    = 0 ! numerical value of double IR pole using actual norm-convention
  real(qp), save      :: de2_i_shift = 0 ! double pole shift defining actual norm convention

  real(qp), save      :: muren
  real(qp), save      :: mureg
  real(qp), save      :: x_UV  = 1       ! rescaling factor for dim-reg scale in UV-divergent quantities
  real(qp), save      :: x_IR  = 1       ! rescaling factor for dim-reg scale in IR-divergent quantities
  real(qp), parameter :: kappa = 2/3._qp ! kappa parameter used in dipole subtraction

  real(qp), save :: LambdaMC2
  real(qp), save :: LambdaMB2
  real(qp), save :: LambdaMT2
  real(qp), save :: LambdaYC2
  real(qp), save :: LambdaYB2
  real(qp), save :: LambdaYT2


  ! the following derived parameters are initilised by subroutine loop_parameters_init
  real(qp), save :: de2_0_IR  ! numerical value of double IR pole using LH-accord convention (i=0)
  real(qp), save :: de2_1_IR  ! numerical value of double IR pole using COLI convention (i=1)
  real(qp), save :: muren2    ! squared renormalisation scale
  real(qp), save :: mureg2    ! squared regularization scale
  real(qp), save :: mu2_UV    ! dim-reg scale for UV-divergent quantities
  real(qp), save :: mu2_IR    ! dim-reg scale for IR-divergent quantities
  real(qp), save :: muyc2 ! squared yukawa renormalization scale for c quark
  real(qp), save :: muyb2 ! squared yukawa renormalization scale for b quark
  real(qp), save :: muyt2 ! squared yukawa renormalization scale for t quark

  ! the following renormalisation constants are initilised by subroutine QCD_renormalisation
  complex(qp), save :: dZMC     = 0 ! charm-quark mass RC        : MC_bare = MC*(1+dZMC)
  complex(qp), save :: dZMB     = 0 ! bottom-quark mass RC       : MB_bare = MB*(1+dZMB)
  complex(qp), save :: dZMT     = 0 ! top-quark mass RC          : MT_bare = MT*(1+dZMT)
  complex(qp), save :: dZYC     = 0 ! charm-quark yukawa RC      : YC_bare = YC*(1+dZYC)
  complex(qp), save :: dZYB     = 0 ! bottom-quark yukawa RC     : YB_bare = YB*(1+dZYB)
  complex(qp), save :: dZYT     = 0 ! top-quark yukawa RC        : YT_bare = YT*(1+dZYT)
  real(qp),    save :: dZg      = 0 ! gluon-field RC             : G_bare  = (1+dZg/2)*G_ren
  real(qp),    save :: dZq      = 0 ! massless-quark field RC    : Q_bare  = (1+dZq/2)*Q_ren
  real(qp),    save :: dZc      = 0 ! charm-quark field RC       : idem
  real(qp),    save :: dZb      = 0 ! bottom-quark field RC      : idem
  real(qp),    save :: dZt      = 0 ! top-quark field RC         : idem
  real(qp),    save :: dgQCD    = 0 ! strong coupling RC         : g_bare  = (1+dgQCD)*g_ren
  real(qp),    save :: dgQCDym  = 0 ! YM-contribution to delnG
  real(qp),    save :: dgQCDfer = 0 ! fermionic-contribution to delnG

  ! Counter terms for QCD corrections
  complex(qp), save :: ctqq(2) ! massless quark propagator counter term
  complex(qp), save :: ctcc(2) ! charm quark propagator counter term
  complex(qp), save :: ctbb(2) ! bottom quark propagator counter term
  complex(qp), save :: cttt(2) ! top quark propagator counter term
  complex(qp), save :: ctGG(3) ! gluon propagator counter term
  real(qp),    save :: ctGqq   ! massless quark-gluon vertex counter term
  real(qp),    save :: ctGcc   ! charm quark-gluon vertex counter term (massive or massless c)
  real(qp),    save :: ctGbb   ! bottom quark-gluon vertex counter term (massive or massless b)
  real(qp),    save :: ctGtt   ! top quark-gluon vertex counter term
  real(qp),    save :: ctVVV   ! three gluon vertex counter term
  real(qp),    save :: ctVVVV  ! four gluon vertex counter term (times 1/2)
  real(qp),    save :: ctVdt   ! Wdt (massless d) vertex counter term
  real(qp),    save :: ctVsc   ! Wcs (massless s) vertex counter term
  real(qp),    save :: ctVst   ! Wts (massless s) vertex counter term
  real(qp),    save :: ctVbu   ! Wub (massless u) vertex counter term
  real(qp),    save :: ctVbc   ! Wcb (massive or massless b/c) vertex counter term
  real(qp),    save :: ctVbt   ! Wtb (massive or massless b) vertex counter term
  real(qp),    save :: ctVtt   ! Att and Ztt vertex counter term
  real(qp),    save :: ctVcc   ! Acc and Zcc vertex counter term (massive or massless c)
  real(qp),    save :: ctVbb   ! Abb and Zbb vertex counter term (massive or massless b)
  real(qp),    save :: ctVqq   ! Aqq and Zqq (massless q) vertex counter term
  complex(qp), save :: ctSud(2), ctSus(2), ctSub(2)
  complex(qp), save :: ctScd(2), ctScs(2), ctScb(2)
  complex(qp), save :: ctStd(2), ctSts(2), ctStb(2)
  complex(qp), save :: ctSdu(2), ctSdc(2), ctSdt(2)
  complex(qp), save :: ctSsu(2), ctSsc(2), ctSst(2)
  complex(qp), save :: ctSbu(2), ctSbc(2), ctSbt(2)
  complex(qp), save :: ctSqq
  complex(qp), save :: ctScc
  complex(qp), save :: ctSbb
  complex(qp), save :: ctStt

  ! Additional parameters for R2
  complex(qp), save :: MQ2sum, MQ2sumpairs
  complex(qp), save :: YQD2sum, YQU2sum
  complex(qp), save :: YQD2sumpairs, YQU2sumpairs

  ! Additional counterterms for R2 QCD
  complex(qp), save :: ctZGG
  complex(qp), save :: ctHGG
  complex(qp), save :: ctAAGG
  complex(qp), save :: ctAZGG
  complex(qp), save :: ctZZGG
  complex(qp), save :: ctWWGG
  complex(qp), save :: ctHHGG
  complex(qp), save :: ctHXGG
  complex(qp), save :: ctXXGG
  complex(qp), save :: ctPPGG
  complex(qp), save :: ctAGGG(2)
  complex(qp), save :: ctZGGG(2)
  integer,           save :: R2GGGG

  ! Counterterms for HEFT
  real(qp), save :: ctHEFTggh(5)
  real(qp), save :: ctHEFTgggh
  real(qp), save :: ctHEFTggggh
  real(qp), save :: R2HEFTggggh
  real(qp), save :: R2HEFThqq
  real(qp), save :: R2HEFTghqq

  ! 2HDM
  complex(qp), save :: thdmctHpsc(2), thdmctHpbt(2), thdmctHpcs(2), thdmctHptb(2)
  complex(qp), save :: thdmctHGG, thdmctHh0GG, thdmctHHGG, thdmctHHhGG
  complex(qp), save :: thdmctXA0GG, thdmctHhHhGG, thdmctA0A0GG
  complex(qp), save :: thdmctPHpGG, thdmctHpHpGG

  ! EW_renormalisation renormalisation constants
  complex(qp), save :: dZMBEW     = 0 ! bottom-quark mass RC       : MB_bare = MB+dZMBEW
  complex(qp), save :: dZMTEW     = 0 ! top-quark mass RC          : MT_bare = MT+dZMTEW
  complex(qp), save :: dZMLEW     = 0 ! tau-lepton mass RC         : ML_bare = ML+dZMLEW
  complex(qp), save :: dZMEEW     = 0 ! electron mass RC           : ME_bare = ME+dZMLEW
  complex(qp), save :: dZMMEW     = 0 ! muon mass RC               : MM_bare = MM+dZMLEW
  complex(qp), save :: dZMW2EW    = 0 ! W mass RC                  : MW^2_bare = MW^2+dZMW2EW^2
  complex(qp), save :: dZMZ2EW    = 0 ! Z mass RC                  : MZ^2_bare = MZ^2+dZMZ2EW^2
  complex(qp), save :: dZMH2EW    = 0 ! H mass RC                  : MH^2_bare = MH^2+dZMH2EW^2
  complex(qp), save :: dswEW      = 0 ! sin EW mixing angle RC         : sw_bare = sw + dswEW  i.e. dswEW/swEW = - c^2/s^2 dcwEW/c ; dcEW/c=1/2(dZMW2/MW^2-dZMZ2/MZ^2)
  complex(qp), save :: dcwEW      = 0 ! cos EW mixing angle RC         : dcwEW/c=1/2(dZMW2/MW^2-dZMZ2/MZ^2) defined for convinience

!   complex(rp), save      :: dZqLEW     = 0 ! L-massless-quark field RC : Q_bare  = (1+1/2*dZqLEW)*Q_ren
!   complex(rp), save      :: dZqREW     = 0 ! R-massless-quark field RC : Q_bare  = (1+1/2*dZqREW)*Q_ren
  complex(qp), save :: dZuLEW     = 0 ! L-massless-u-quark field RC : Q_bare  = (1+1/2*dZuLEW)*Q_ren
  complex(qp), save :: dZuREW     = 0 ! R-massless-u-quark field RC : Q_bare  = (1+1/2*dZuREW)*Q_ren
  complex(qp), save :: dZdLEW     = 0 ! L-massless-d-quark field RC : Q_bare  = (1+1/2*dZdLEW)*Q_ren
  complex(qp), save :: dZdREW     = 0 ! R-massless-d-quark field RC : Q_bare  = (1+1/2*dZdREW)*Q_ren
  complex(qp), save :: dZbLEW     = 0 ! L-bottom-quark field RC     : idem
  complex(qp), save :: dZbREW     = 0 ! R-bottom-quark field RC     : idem
  complex(qp), save :: dZtLEW     = 0 ! L-top-quark field RC        : idem
  complex(qp), save :: dZtREW     = 0 ! R-top-quark field RC        : idem
  complex(qp), save :: dZeLEW     = 0 ! L-electron field RC         : idem
  complex(qp), save :: dZeREW     = 0 ! R-electron field RC         : idem
  complex(qp), save :: dZmuLEW    = 0 ! L-muon field RC             : idem
  complex(qp), save :: dZmuREW    = 0 ! R-muon field RC             : idem
  complex(qp), save :: dZnLEW     = 0 ! L-neutrino field RC         : idem
  complex(qp), save :: dZnlLEW    = 0 ! L-tau-neutrino field RC     : idem
  complex(qp), save :: dZlLEW     = 0 ! L-tau-lepton field RC       : idem
  complex(qp), save :: dZlREW     = 0 ! R-tau-lepton field RC       : idem
  complex(qp), save :: dZuLEWcc   = 0 ! L-massless-u-quark field RC : Q_bare  = (1+1/2*dZuLEW)*Q_ren
  complex(qp), save :: dZuREWcc   = 0 ! R-massless-u-quark field RC : Q_bare  = (1+1/2*dZuREW)*Q_ren
  complex(qp), save :: dZdLEWcc   = 0 ! L-massless-d-quark field RC : Q_bare  = (1+1/2*dZdLEW)*Q_ren
  complex(qp), save :: dZdREWcc   = 0 ! R-massless-d-quark field RC : Q_bare  = (1+1/2*dZdREW)*Q_ren
  complex(qp), save :: dZbLEWcc   = 0 ! L-bottom-quark field RC     : idem
  complex(qp), save :: dZbREWcc   = 0 ! R-bottom-quark field RC     : idem
  complex(qp), save :: dZtLEWcc   = 0 ! L-top-quark field RC        : idem
  complex(qp), save :: dZtREWcc   = 0 ! R-top-quark field RC        : idem
  complex(qp), save :: dZeLEWcc   = 0 ! L-electron field RC         : idem
  complex(qp), save :: dZeREWcc   = 0 ! R-electron field RC         : idem
  complex(qp), save :: dZmuLEWcc  = 0 ! L-muon field RC             : idem
  complex(qp), save :: dZmuREWcc  = 0 ! R-muon field RC             : idem
  complex(qp), save :: dZnLEWcc   = 0 ! L-neutrino field RC         : idem
  complex(qp), save :: dZnlLEWcc  = 0 ! L-tau-neutrino field RC     : idem
  complex(qp), save :: dZlLEWcc   = 0 ! L-tau-lepton field RC       : idem
  complex(qp), save :: dZlREWcc    = 0 ! R-tau-lepton field RC       : idem

  complex(qp), save :: dZWEW      = 0 ! W field RC                 : idem
  complex(qp), save :: dZZZEW     = 0 ! ZZ field RC                : idem
  complex(qp), save :: dZAZEW     = 0 ! AZ field RC                : idem
  complex(qp), save :: dZZAEW     = 0 ! AZ field RC                : idem
  complex(qp), save :: dZAAEW     = 0 ! AA field RC                : idem
  complex(qp), save :: dZAAEWdreg = 0 ! AA field RC dim-reg      : Photon field renormalization constant with ligh fermion contributions in dim-reg
  complex(qp), save :: dZAAEWnreg = 0 ! AA field RC dim-reg      : Photon field renormalization constant with ligh fermion contributions in n-reg
  complex(qp), save :: dZHEW      = 0 ! H field RC                 : idem

  complex(qp), save :: dtEW       = 0 ! tadpole-RC                 :
  complex(qp), save :: dZeQEDEW   = 0 ! EW coupling RC         : e_bare  = (1+dZeEW)*e_ren
  complex(qp), save :: dZe0QEDEW  = 0 ! EW coupling RC         : e_bare  = (1+dZeEW)*e_ren in on-shell/a(0) scheme
  complex(qp), save :: dZe0QEDEWdreg = 0 ! EW coupling RC   : e_bare  = (1+dZeEW)*e_ren in on-shell/a(0) scheme with light fermuion controbutions in dim-reg
  complex(qp), save :: dZe0QEDEWnreg = 0 ! EW coupling RC   : e_bare  = (1+dZeEW)*e_ren in on-shell/a(0) scheme with light fermuion controbutions in n-reg
  complex(qp), save :: dZeGmuQEDEW   = 0 ! EW coupling RC   : e_bare  = (1+dZeEW)*e_ren in Gmu-scheme
  complex(qp), save :: dZeZQEDEW     = 0 ! EW coupling RC   : e_bare  = (1+dZeEW)*e_ren in a(MZ)-scheme
  complex(qp), save :: dVLEW         = 0 ! L-lepton-Z vertex RC
  complex(qp), save :: dVREW         = 0 ! R-lepton-Z vertex RC

  ! Counter terms for EW corrections
  ! VV Vector propagators
  complex(qp), save :: EWctWW(3)      = 0
  complex(qp), save :: EWctZZ(3)      = 0
  complex(qp), save :: EWctAZ(3)      = 0
  complex(qp), save :: EWctAA(3)      = 0
  ! SS scalar propagators
  complex(qp), save :: EWctHH(2)      = 0
  complex(qp), save :: EWctXX(2)      = 0
  complex(qp), save :: EWctPP(2)      = 0
  ! SV scalar-vector mixing
  complex(qp), save :: EWctXA      = 0
  complex(qp), save :: EWctXZ      = 0
  complex(qp), save :: EWctPW      = 0
  ! FF fermionic propagators
  complex(qp), save :: EWctuu(4)      = 0
  complex(qp), save :: EWctdd(4)      = 0
  complex(qp), save :: EWcttt(4)      = 0
  complex(qp), save :: EWctbb(4)      = 0
  complex(qp), save :: EWctee(4)      = 0
  complex(qp), save :: EWctmm(4)      = 0
  complex(qp), save :: EWctLL(4)      = 0
  complex(qp), save :: EWctnn(4)      = 0
  complex(qp), save :: EWctnlnl(4)      = 0
  !VVVV
  complex(qp), save :: EWctWWWW(2)      = 0
  complex(qp), save :: EWctWWZZ(2)      = 0
  complex(qp), save :: EWctWWAZ(2)      = 0
  complex(qp), save :: EWctWWAA(2)      = 0
  !VVVV pure R2
  complex(qp), save :: EWctR2AAAA      = 0
  complex(qp), save :: EWctR2AAAZ      = 0
  complex(qp), save :: EWctR2AAZZ      = 0
  complex(qp), save :: EWctR2AZZZ      = 0
  complex(qp), save :: EWctR2ZZZZ      = 0
  !VVV
  complex(qp), save :: EWctAWW      = 0
  complex(qp), save :: EWctZWW      = 0
  !SSSS
  complex(qp), save :: EWctSSSS1      = 0
  complex(qp), save :: EWctSSSS2      = 0
  complex(qp), save :: EWctSSSS3      = 0
  complex(qp), save :: EWctHHHH      = 0
  complex(qp), save :: EWctHHXX      = 0
  complex(qp), save :: EWctHHPP      = 0
  complex(qp), save :: EWctXXXX      = 0
  complex(qp), save :: EWctXXPP      = 0
  complex(qp), save :: EWctPPPP      = 0
  !SSS
  complex(qp), save :: EWctHHH      = 0
  complex(qp), save :: EWctHXX      = 0
  complex(qp), save :: EWctHPP      = 0
  !VVSS
  complex(qp), save :: EWctWWXX      = 0
  complex(qp), save :: EWctWWHH      = 0
  complex(qp), save :: EWctWWPP      = 0
  complex(qp), save :: EWctZZPP      = 0
  complex(qp), save :: EWctZAPP      = 0
  complex(qp), save :: EWctAAPP      = 0
  complex(qp), save :: EWctZZHH      = 0
  complex(qp), save :: EWctZZXX      = 0
  complex(qp), save :: EWctZAHH      = 0
  complex(qp), save :: EWctZAXX      = 0
  complex(qp), save :: EWctWZPH      = 0
  complex(qp), save :: EWctWAPH      = 0
  complex(qp), save :: EWctWZPX      = 0
  complex(qp), save :: EWctWAPX      = 0
  !VVSS R2
  complex(qp), save :: EWctAAHH      = 0
  complex(qp), save :: EWctAAXX      = 0
  !VSS
  complex(qp), save :: EWctAXH      = 0
  complex(qp), save :: EWctZXH      = 0
  complex(qp), save :: EWctAPP      = 0
  complex(qp), save :: EWctZPP      = 0
  complex(qp), save :: EWctWPH      = 0
  complex(qp), save :: EWctWPX      = 0
  !SVV
  complex(qp), save :: EWctHWW      = 0
  complex(qp), save :: EWctHZZ      = 0
  complex(qp), save :: EWctHZA      = 0
  complex(qp), save :: EWctPWZ      = 0
  complex(qp), save :: EWctPWA      = 0
  ! pure R2 SVV
  complex(qp), save :: EWctHAA      = 0
  !VFF
  ! Aff
  complex(qp), save :: EWctAuu(2)      = 0
  complex(qp), save :: EWctAdd(2)      = 0
  complex(qp), save :: EWctAtt(2)      = 0
  complex(qp), save :: EWctAbb(2)      = 0
  complex(qp), save :: EWctAee(2)      = 0
  complex(qp), save :: EWctAmm(2)      = 0
  complex(qp), save :: EWctALL(2)      = 0
  complex(qp), save :: EWctAnn(2)      = 0
  ! Zff
  complex(qp), save :: dgZu(2)      = 0
  complex(qp), save :: dgZd(2)      = 0
  complex(qp), save :: dgZl(2)      = 0
  complex(qp), save :: dgZn(2)      = 0
  complex(qp), save :: EWctVuu(2)      = 0
  complex(qp), save :: EWctVdd(2)      = 0
  complex(qp), save :: EWctVtt(2)      = 0
  complex(qp), save :: EWctVbb(2)      = 0
  complex(qp), save :: EWctVee(2)      = 0
  complex(qp), save :: EWctVmm(2)      = 0
  complex(qp), save :: EWctVLL(2)      = 0
  complex(qp), save :: EWctVnn(2)      = 0
  complex(qp), save :: EWctVnlnl(2)      = 0
  ! Wff
  complex(qp), save :: EWctVdu      = 0
  complex(qp), save :: EWctVbt      = 0
  complex(qp), save :: EWctVen      = 0
  complex(qp), save :: EWctVLn      = 0
  complex(qp), save :: EWctVud      = 0
  complex(qp), save :: EWctVtb      = 0
  complex(qp), save :: EWctVne      = 0
  complex(qp), save :: EWctVnL      = 0
  ! Gff mixed EW/QCD
  complex(qp), save ::  EWctGuu(2)      = 0
  complex(qp), save ::  EWctGdd(2)      = 0
  complex(qp), save ::  EWctGtt(2)      = 0
  complex(qp), save ::  EWctGbb(2)      = 0
  !SFF
  complex(qp), save :: EWctHee(2)      = 0
  complex(qp), save :: EWctHmm(2)      = 0
  complex(qp), save :: EWctHLL(2)      = 0
  complex(qp), save :: EWctHtt(2)      = 0
  complex(qp), save :: EWctHbb(2)      = 0
  complex(qp), save :: EWctXee(2)      = 0
  complex(qp), save :: EWctXmm(2)      = 0
  complex(qp), save :: EWctXtt(2)      = 0
  complex(qp), save :: EWctXbb(2)      = 0
  complex(qp), save :: EWctXLL(2)      = 0
  complex(qp), save :: EWctPud(2)      = 0
  complex(qp), save :: EWctPdu(2)      = 0
  complex(qp), save :: EWctPtb(2)      = 0
  complex(qp), save :: EWctPbt(2)      = 0
  complex(qp), save :: EWctPne(2)      = 0
  complex(qp), save :: EWctPen(2)      = 0
  complex(qp), save :: EWctPnL(2)      = 0
  complex(qp), save :: EWctPLn(2)      = 0
  ! VUU
  complex(qp), save :: EWctAUWUW      = 0
  complex(qp), save :: EWctZUWUW      = 0
  complex(qp), save :: EWctWUWUZ      = 0
  complex(qp), save :: EWctWUZUW      = 0
  complex(qp), save :: EWctWUWUA      = 0
  complex(qp), save :: EWctWUAUW      = 0
  ! SUU
  complex(qp), save :: EWctHUZUZ      = 0
  complex(qp), save :: EWctHUWUW      = 0
  complex(qp), save :: EWctXUWUW      = 0
  complex(qp), save :: EWctPUZUW      = 0
  complex(qp), save :: EWctPUWUZ      = 0
  complex(qp), save :: EWctPUWUA      = 0

  ! Additional parameters for R2 EW
  complex(qp), save :: sumMQ2
  complex(qp), save :: sumMQ2Q2
  complex(qp), save :: sumMQ2QI
  complex(qp), save :: sumMQ4
  complex(qp), save :: sumMUD
  complex(qp), save :: sumMUD2
  complex(qp), save :: sumMU2
  complex(qp), save :: sumMD2
  complex(qp), save :: sumMU4
  complex(qp), save :: sumMD4
  complex(qp), save :: sumMQ2QUD
  complex(qp), save :: sumML2
  complex(qp), save :: sumML2Q2
  complex(qp), save :: sumML2QI
  complex(qp), save :: sumML4

end module ol_loop_parameters_decl_qp

