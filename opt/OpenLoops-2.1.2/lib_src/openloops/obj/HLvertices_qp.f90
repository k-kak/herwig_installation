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


module hol_initialisation_qp
  implicit none
  contains

!************************************************************************
subroutine hol_allocation(alpha,rank,beta,hel_states,ol_coeff,m)
!************************************************************************
! Allocation of the OpenLoops coefficient of type hol
!************************************************************************
! alpha      = dimension of the alpha-index array
! rank       = rank dimensionality (N=1 r=0, N=5 r=1,...)
! beta       = dimension of the beta-index array
! hel_states = number of helicity states
! ol_coeff   = type(hol) OpenLoops coefficient to be allocated in memory
! m          = number of ol_coeff with the same number of hel_states to
!              be allocated
!************************************************************************
  use kind_types, only: qp
  use ol_data_types_qp, only: hol


  integer,   intent(in)    :: alpha, rank, beta, hel_states, m
  type(hol), intent(inout) :: ol_coeff(:)

  integer :: i

  do i = 1, m
    allocate(ol_coeff(i)%hf(hel_states))
    ! TODO: Initialization needed?
    ol_coeff(i)%hf = 0
    allocate(ol_coeff(i)%j(alpha, rank, beta, hel_states))
    ol_coeff(i)%j = 0
    ol_coeff(i)%error = 0
    ol_coeff(i)%ndrs = 0
    ol_coeff(i)%nred = 0

  end do

end subroutine hol_allocation

!************************************************************************
subroutine hol_deallocation(ol_coeff, m, dmode)
!************************************************************************
! Allocation of the OpenLoops coefficient of type hol
!************************************************************************
! alpha      = dimension of the alpha-index array
! rank       = rank dimensionality (N=1 r=0, N=5 r=1,...)
! beta       = dimension of the beta-index array
! hel_states = number of helicity states
! ol_coeff   = type(hol) OpenLoops coefficient to be allocated in memory
! m          = number of ol_coeff with the same number of hel_stases to
!              be allocated
! dmode      = deallocation mode (0 = dp + qp, 1 = only qp) ! TODO
!************************************************************************
  use kind_types, only: qp
  use ol_data_types_qp, only: hol


  integer,   intent(in)    :: m, dmode
  type(hol), intent(inout) :: ol_coeff(:)

  integer :: i

  do i = 1, m
    if(dmode == 0) then
      if(allocated(ol_coeff(i)%hf)) deallocate(ol_coeff(i)%hf)
      if(allocated(ol_coeff(i)%j)) deallocate(ol_coeff(i)%j)
    endif
    ol_coeff(i)%error = 0

  end do

end subroutine hol_deallocation

!************************************************************************
subroutine hcl_allocation(rank,ol_coeff, m)
!************************************************************************
! Allocation of the OpenLoops coefficient of type hcl
!************************************************************************
! rank       = rank dimensionality (N=1 r=0, N=5 r=1,...)
! beta       = dimension of the beta-index array
! ol_coeff   = type(hcl) OpenLoops coefficient to be allocated in memory
! m          = number of ol_coeff with the same number of hel_stases to
!              be allocated
!************************************************************************
  use kind_types, only: qp
  use ol_data_types_qp, only: hcl


  integer,   intent(in)    :: rank, m
  type(hcl), intent(inout) :: ol_coeff(:)

  integer :: i

  do i = 1, m
    allocate(ol_coeff(i)%cmp(rank))
    ! TODO: initialization of cmp needed?
    ol_coeff(i)%cmp(:) = 0
    ol_coeff(i)%error = 0
    ol_coeff(i)%ndrs = 0
    ol_coeff(i)%nred = 0

  end do

end subroutine hcl_allocation

!************************************************************************
subroutine hcl_deallocation(ol_coeff, m, dmode)
!************************************************************************
! Allocation of the OpenLoops coefficient of type hcl
!************************************************************************
! alpha      = dimension of the alpha-index array
! rank       = rank dimensionality (N=1 r=0, N=5 r=1,...)
! beta       = dimension of the beta-index array
! hel_states = number of helicity states
! ol_coeff   = type(hol) OpenLoops coefficient to be allocated in memory
! m          = number of ol_coeff with the same number of hel_stases to
!              be allocated
! dmode      = deallocation mode (0 = dp + qp, 1 = only qp) ! TODO
!************************************************************************
  use kind_types, only: qp
  use ol_data_types_qp, only: hcl

  implicit none

  integer,   intent(in)    :: m, dmode
  type(hcl), intent(inout) :: ol_coeff(:)

  integer :: i

  do i = 1, m
    if(dmode == 0) then
      if(allocated(ol_coeff(i)%cmp)) deallocate(ol_coeff(i)%cmp)
    endif
    ol_coeff(i)%error = 0

  end do

end subroutine hcl_deallocation

!***********************************************************************
subroutine G0_hol_initialisation(ntry,G0coeff,ol_coeff,nhel_in,h0t, &
                                 momids,massids,nsubtrees,n_wf,w1,w2,w3,w4,w5)
!***********************************************************************
! Initialisation of the OpenLoops coefficient
!***********************************************************************
! G0coeff    = array of size given by the global number of helicity states
!              with the coefficients needed for the tree-level 1-loop interference
! hel_states = array with the helicity configurations
! ol_coeff   = OpenLoops coefficient to be initialised
! m0         = number of non-vanishing helicity states in this diagram
! h0t        = list of non-vanishing helicity states in this diagram
! momids     = momenta entering the loop (binary representation)
! massids    = masses in the denominators
! nsubstrees = number of subtrees
! n_wf       = number of subtrees attached to the full loop consisting of more than 1 external particles
! w1,...,w5  = arrays of hf (helicity) integers for the subtrees consisting of more than 1 external particles
!***********************************************************************
  use kind_types, only: dp, qp, qp, intkind1, intkind2
  use ol_data_types_qp, only: hol, Hpolcont, wfun
  use ol_loop_parameters_decl_dp, only: bubble_vertex, ir_hacks
  use ol_loop_routines_qp, only: G0initialisationOLR
  use ol_debug, only: ol_fatal, ol_msg
  use ind_bookkeeping_qp, only : ProjHind
  use ol_parameters_decl_qp, only: hp_switch,hp_ir_trig, &
                                             hp_fake_trig
  use ol_loop_handling_dp, only: hybrid_dp_mode,hybrid_qp_mode,&
                                            hybrid_dp_qp_mode,hybrid_zero_mode

  implicit none

  integer(intkind1), intent(in)  :: ntry
  integer(intkind2), intent(inout) :: h0t(:), nhel_in
  integer, intent(in)  :: nsubtrees,n_wf
  type(Hpolcont), intent(in)     :: G0coeff(:)
  integer, intent(in) :: momids(:)
  integer, intent(in) :: massids(:)
  type(hol), intent(inout)       :: ol_coeff
  type(wfun), optional, intent(in) :: w1(:), w2(:), w3(:), w4(:), w5(:)
  integer :: l, htot, sw, n
  integer(intkind2)  :: n_expart
  integer(intkind2)  :: G0_new_hf(size(G0coeff))
  integer(intkind2)  :: h1, h3
  integer(intkind2)  :: nhel_wf
  integer(intkind2)  :: hel3, hel1, subset
  logical :: check_hel_ol(size(G0coeff)),selfenergy

!***********************************************************************************

  ol_coeff%mode = hybrid_dp_mode
  ol_coeff%error = 0._qp
  ol_coeff%npoint = size(momids)
  ol_coeff%ndrs = 0
  ol_coeff%nred = 0


  if (ntry == 1) then

    htot = size(G0coeff)
    do l = 1, htot
      ol_coeff%hf(l) = G0coeff(l)%hf
      h0t(l)=l
    end do

    ! check number of non-vanishing helicity configurations of the input OL coefficient
    nhel_in = 0
    do h3 = 1, htot
      if (ol_coeff%hf(h3) /= -1_intkind2) nhel_in = nhel_in + 1
    end do

    !!!============================= 1 or more external subtrees
    if (n_wf >= 1) then
      sw=size(w1)
      ! Check number of non-vanishing helicity configurations of the external wavefunction
      nhel_wf = 0
      do h1 = 1, sw
        if (w1(h1)%hf /= -1_intkind2) nhel_wf = nhel_wf + 1
      end do

      subset = w1(1)%t
      n_expart = w1(1)%n_part
      check_hel_ol = .true.
      do h3=1, nhel_in
        hel3=ol_coeff%hf(h3)
        do h1=1, nhel_wf
          hel1 = w1(h1)%hf
          if (check_hel_ol(h3)) then
            if (ProjHind(subset,hel3,n_expart)==hel1) then
              if(check_hel_ol(h3)) check_hel_ol(h3)=.false.
            end if
          end if
        end do
      end do

      h3=1
      do while (h3 <= nhel_in)
        if(check_hel_ol(h3)) then
          if (h3 > 1) G0_new_hf(1:h3-1)=ol_coeff%hf(1:h3-1)
          G0_new_hf(h3:nhel_in-1)=ol_coeff%hf(h3+1:nhel_in)
          G0_new_hf(nhel_in)=-1_intkind2
          ol_coeff%hf(1:nhel_in)=G0_new_hf(1:nhel_in)
          h0t(h3:nhel_in-1)=h0t(h3+1:nhel_in)
          check_hel_ol(1:nhel_in-1)=check_hel_ol(2:nhel_in)
          nhel_in=nhel_in-1
        else
          h3=h3+1
        end if
      end do
    end if

  !!!============================= 2 or more external subtrees
    if (n_wf >= 2) then
      sw=size(w2)
      ! Check number of non-vanishing helicity configurations of the external wavefunction
      nhel_wf = 0
      do h1 = 1, sw
        if (w2(h1)%hf /= -1_intkind2) nhel_wf = nhel_wf + 1
      end do

      subset = w2(1)%t
      n_expart = w2(1)%n_part
      check_hel_ol = .true.

      do h3=1, nhel_in
        hel3=ol_coeff%hf(h3)
        do h1=1, nhel_wf
          hel1 = w2(h1)%hf
          if (check_hel_ol(h3)) then
            if (ProjHind(subset,hel3,n_expart)==hel1) then
              if(check_hel_ol(h3)) check_hel_ol(h3)=.false.
            end if
          end if
        end do
      end do

      h3=1
      do while (h3 <= nhel_in)
        if(check_hel_ol(h3)) then
          if (h3 > 1) G0_new_hf(1:h3-1)=ol_coeff%hf(1:h3-1)
          G0_new_hf(h3:nhel_in-1)=ol_coeff%hf(h3+1:nhel_in)
          G0_new_hf(nhel_in)=-1_intkind2
          ol_coeff%hf(1:nhel_in)=G0_new_hf(1:nhel_in)
          h0t(h3:nhel_in-1)=h0t(h3+1:nhel_in)
          check_hel_ol(1:nhel_in-1)=check_hel_ol(2:nhel_in)
          nhel_in=nhel_in-1
        else
          h3=h3+1
        end if
      end do
    end if
  !!!============================= 3 or more external subtrees
    if (n_wf >= 3) then
      sw=size(w3)
      ! Check number of non-vanishing helicity configurations of the external wavefunction
      nhel_wf = 0
      do h1 = 1, sw
        if (w3(h1)%hf /= -1_intkind2) nhel_wf = nhel_wf + 1
      end do

      subset = w3(1)%t
      n_expart = w3(1)%n_part
      check_hel_ol = .true.

      do h3=1, nhel_in
      hel3=ol_coeff%hf(h3)
        do h1=1, nhel_wf
          hel1 = w3(h1)%hf
          if (check_hel_ol(h3)) then
            if (ProjHind(subset,hel3,n_expart)==hel1) then
              if(check_hel_ol(h3)) check_hel_ol(h3)=.false.
            end if
          end if
        end do
      end do

      h3=1
      do while (h3 <= nhel_in)
        if(check_hel_ol(h3)) then
          if (h3 > 1) G0_new_hf(1:h3-1)=ol_coeff%hf(1:h3-1)
          G0_new_hf(h3:nhel_in-1)=ol_coeff%hf(h3+1:nhel_in)
          G0_new_hf(nhel_in)=-1_intkind2
          ol_coeff%hf(1:nhel_in)=G0_new_hf(1:nhel_in)
          h0t(h3:nhel_in-1)=h0t(h3+1:nhel_in)
          check_hel_ol(1:nhel_in-1)=check_hel_ol(2:nhel_in)
          nhel_in=nhel_in-1
        else
          h3=h3+1
        end if
      end do
    end if
  !!!============================= 4 or more external subtrees
    if (n_wf >= 4) then
      sw=size(w4)
      ! Check number of non-vanishing helicity configurations of the external wavefunction
      nhel_wf = 0
      do h1 = 1, sw
        if (w4(h1)%hf /= -1_intkind2) nhel_wf = nhel_wf + 1
      end do

      subset = w4(1)%t
      n_expart = w4(1)%n_part
      check_hel_ol = .true.

      do h3=1, nhel_in
      hel3=ol_coeff%hf(h3)
        do h1=1, nhel_wf
          hel1 = w4(h1)%hf
          if (check_hel_ol(h3)) then
            if (ProjHind(subset,hel3,n_expart)==hel1) then
              if(check_hel_ol(h3)) check_hel_ol(h3)=.false.
            end if
          end if
        end do
      end do

      h3=1
      do while (h3 <= nhel_in)
        if(check_hel_ol(h3)) then
          if (h3 > 1) G0_new_hf(1:h3-1)=ol_coeff%hf(1:h3-1)
          G0_new_hf(h3:nhel_in-1)=ol_coeff%hf(h3+1:nhel_in)
          G0_new_hf(nhel_in)=-1_intkind2
          ol_coeff%hf(1:nhel_in)=G0_new_hf(1:nhel_in)
          h0t(h3:nhel_in-1)=h0t(h3+1:nhel_in)
          check_hel_ol(1:nhel_in-1)=check_hel_ol(2:nhel_in)
          nhel_in=nhel_in-1
        else
          h3=h3+1
        end if
      end do
    end if
  !!!============================= 5 or more external subtrees
    if (n_wf >= 5) then
      sw=size(w5)
      ! Check number of non-vanishing helicity configurations of the external wavefunction
      nhel_wf = 0
      do h1 = 1, sw
        if (w5(h1)%hf /= -1_intkind2) nhel_wf = nhel_wf + 1
      end do

      subset = w5(1)%t
      n_expart = w5(1)%n_part
      check_hel_ol = .true.

      do h3=1, nhel_in
        hel3=ol_coeff%hf(h3)
        do h1=1, nhel_wf
          hel1 = w5(h1)%hf
          if (check_hel_ol(h3)) then
            if (ProjHind(subset,hel3,n_expart)==hel1) then
              if(check_hel_ol(h3)) check_hel_ol(h3)=.false.
            end if
          end if
        end do
      end do

      h3=1
      do while (h3 <= nhel_in)
        if(check_hel_ol(h3)) then
          if (h3 > 1) G0_new_hf(1:h3-1)=ol_coeff%hf(1:h3-1)
          G0_new_hf(h3:nhel_in-1)=ol_coeff%hf(h3+1:nhel_in)
          G0_new_hf(nhel_in)=-1_intkind2
          ol_coeff%hf(1:nhel_in)=G0_new_hf(1:nhel_in)
          h0t(h3:nhel_in-1)=h0t(h3+1:nhel_in)
          check_hel_ol(1:nhel_in-1)=check_hel_ol(2:nhel_in)
          nhel_in=nhel_in-1
        else
          h3=h3+1
        end if
      end do
    end if
!!!============================= 6 or more external subtrees
    if (n_wf >= 6) then
      call ol_msg("More than 5 external subtrees consisting of >=2 particles")
      call ol_fatal()
    end if
!!!=============================
  end if
!***********************************************************************************

  ol_coeff%j(:,:,:,nhel_in+1:size(ol_coeff%hf)) = 0._qp

  do l = 1, nhel_in
    call G0initialisationOLR(G0coeff(h0t(l))%j,ol_coeff%j(:,:,:,l))
  end do
  selfenergy = (ol_coeff%npoint .eq. 2) .and. (nsubtrees .eq. 2)
  if (bubble_vertex .eq. 1 .and. selfenergy) then
    ol_coeff%j(:,:,:,1:nhel_in) = 0._qp
    ! breaks helicity mapping
    !ol_coeff%mode = hybrid_zero_mode
  end if



  contains

  function check_collinear_tr()
    use ol_momenta_decl_qp, only: collconf
    use ol_kinematics_qp, only: get_mass2
    logical :: check_collinear_tr,m1ext,m2ext,m3ext,zeromass
    integer :: m1, m2, m3

    m1 = momids(1)
    m1ext = iand(m1,m1-1) .eq. 0
    m2 = momids(2)
    m2ext = iand(m2,m2-1) .eq. 0
    m3 = momids(3)
    m3ext = iand(m3,m3-1) .eq. 0
    zeromass = sum(get_mass2(massids)) .eq. 0

    check_collinear_tr = .false.
    if (zeromass) then
      if (m1ext .and. m2ext) then
        if (iand(m1,collconf) .ne. 0 .and. iand(m2,collconf) .ne. 0) then
          check_collinear_tr = .true.
        end if
      else if (m1ext .and. m3ext) then
        if (iand(m1,collconf) .ne. 0 .and. iand(m3,collconf) .ne. 0) then
          check_collinear_tr = .true.
        end if
      else if (m2ext .and. m3ext) then
        if (iand(m2,collconf) .ne. 0 .and. iand(m3,collconf) .ne. 0) then
          check_collinear_tr = .true.
        end if
      end if
    end if

  end function check_collinear_tr

  function check_soft_tr()
    use ol_momenta_decl_qp, only: softconf
    logical :: check_soft_tr,m1ext,m2ext,m3ext
    integer :: m1, m2, m3

    m1 = momids(1)
    m1ext = iand(m1,m1-1) .eq. 0
    m2 = momids(2)
    m2ext = iand(m2,m2-1) .eq. 0
    m3 = momids(3)
    m3ext = iand(m3,m3-1) .eq. 0

    check_soft_tr = .false.
    if (m1ext .and. m2ext) then
      if (m1 .eq. softconf .or. m2 .eq. softconf) then
        check_soft_tr = .true.
      end if
    else if (m1ext .and. m3ext) then
      if (m1 .eq. softconf .or. m3 .eq. softconf) then
        check_soft_tr = .true.
      end if
    else if (m2ext .and. m3ext) then
      if (m2 .eq. softconf .or. m3 .eq. softconf) then
        check_soft_tr = .true.
      end if
    end if

  end function check_soft_tr

  function check_quasios()
    use ol_momenta_decl_qp, only: softconf
    logical :: check_quasios
    integer :: m1, m2

    m1 = momids(1)
    m2 = momids(2)

    check_quasios = .false.
    if (iand(m1,softconf) .ne. 0) then
      if (iand(m1 -softconf, m1 -softconf -1) .eq. 0) then
        check_quasios = .true.
      end if
    else if (iand(m2,softconf) .ne. 0) then
      if (iand(m2 -softconf, m2 -softconf -1) .eq. 0) then
        check_quasios = .true.
      end if
    end if

  end function check_quasios

  function check_collinear_leg()
    use ol_momenta_decl_qp, only: collconf
    logical :: check_collinear_leg,m1ext,m2ext,m3ext
    integer :: n, m1, m2, m3

    m1 = momids(1)
    m1ext = iand(m1,m1-1) .eq. 0
    m2 = momids(2)
    m2ext = iand(m2,m2-1) .eq. 0
    m3 = momids(3)
    m3ext = iand(m3,m3-1) .eq. 0

    check_collinear_leg = .false.
    do n = 1, size(momids)
      if (momids(n) .eq. collconf) then
        check_collinear_leg = .true.
      end if
    end do

  end function check_collinear_leg


end subroutine G0_hol_initialisation


end module hol_initialisation_qp





module ol_h_vert_interface_qp


  implicit none
  interface valid_hol
    module procedure valid_hol_hol, valid_hol_hcl
  end interface valid_hol
  contains

function valid_hol_hol(Gin, Gout) result(valid_hol)
  use ol_data_types_qp, only: hol
  use ol_loop_handling_dp, only: hybrid_zero_mode

  type(hol), intent(in)    :: Gin
  type(hol), intent(inout) :: Gout
  logical :: valid_hol

  Gout%mode = Gin%mode
  if (Gin%mode .eq. hybrid_zero_mode) then
    valid_hol = .false.
    Gout%j = 0
    Gout%error = 0
    Gout%npoint = Gin%npoint
    Gout%ndrs = 0
    Gout%nred = 0

  else
    valid_hol = .true.
    Gout%error = Gin%error
    Gout%npoint = Gin%npoint

    Gout%ndrs = Gin%ndrs + 1

    Gout%nred = Gin%nred

  end if

end function valid_hol_hol

function valid_hol_hcl(Gin, Gout) result(valid_hol)
  use ol_data_types_qp, only: hol, hcl
  use ol_loop_handling_dp, only: hybrid_zero_mode

  type(hol), intent(in)    :: Gin
  type(hcl), intent(inout) :: Gout
  logical :: valid_hol

  Gout%mode  = Gin%mode
  if (Gin%mode .eq. hybrid_zero_mode) then
    valid_hol = .false.
    Gout%cmp = 0
    Gout%error = 0
    Gout%ndrs = 0
    Gout%nred = 0

  else
    valid_hol = .true.
    Gout%error = Gin%error

    Gout%ndrs = Gin%ndrs + 1

    Gout%nred = Gin%nred

  end if

end function valid_hol_hcl


!***********************************************************************
! OpenLoops dressing steps with summation over the helicities of the
! external wave functions.
! ----------------------------------------------------------------------
! ntry          = 1 (2) for 1st (subsequent) PS points
! G_X(1:n(3))   = input open-loop with n(3) helicity states
!                 corresponding to "unattached" external legs
! J_Y(1:n(1))   = external wave function with n(1) helicity states
! Gout_Z(1:n(2))= output open loop with n(2) helicity states
!                 corresponding to "unattached" external legs
! n(1:3) = n(1) is the number of helicity states of the external subtree
!          n(2) is the number of helicity states of the output open-loop
!          n(3) is number of non-vanishing helicity states of the input
! t(1:2,n(3)) = helicity table. It links the helicity states of G_X(h)
! with J_Y(t(1,h)) and Gout_Z(t(2,h))
! **********************************************************************

!***********************************************************************
subroutine Hloop_AZ_Q(ntry, G_A, J_Z, Gout_A, ng_RL, n, t)
!-----------------------------------------------------------------------
! bare AZ -> A Z-like interaction
! **********************************************************************
  use kind_types, only: qp, qp, intkind1, intkind2
  use ol_data_types_qp, only: wfun, hol
  use ol_vert_interface_qp, only: loop_AZ_Q
  use hel_bookkeeping_qp, only: helbookkeeping_ol_vert3
  use ol_parameters_init_qp, only: get_coupling

  implicit none
  integer(intkind1), intent(in)    :: ntry
  integer(intkind2), intent(inout) :: t(:,:), n(3)
  type(wfun),        intent(in)    :: J_Z(:)
  type(hol),         intent(inout) :: G_A,Gout_A
  integer,           intent(in)    :: ng_RL
  complex(qp)  :: G_add(size(Gout_A%j,1),size(Gout_A%j,2),size(Gout_A%j,3))
  complex(qp)  :: g_RL(2)

  integer :: h

  if (ntry == 1) call helbookkeeping_ol_vert3(ntry, J_Z, G_A, Gout_A, n, t)
  if (.not. valid_hol(G_A, Gout_A)) return

  Gout_A%j = 0._qp
  g_RL = get_coupling(ng_RL)
  do h = 1, n(3)  ! recursion step
    call loop_AZ_Q(G_A%j(:,:,:,h),J_Z(t(1,h))%j,G_add,g_RL)
    Gout_A%j(:,:,:,t(2,h)) = Gout_A%j(:,:,:,t(2,h)) + G_add
  end do



end subroutine Hloop_AZ_Q


!***********************************************************************
subroutine Hloop_AQ_Z(ntry, G_A, J_Q, Gout_Z, ng_RL, n, t)
!-----------------------------------------------------------------------
! bare AQ -> Z Z-like interaction
! **********************************************************************
  use kind_types, only: qp, intkind1, intkind2, qp
  use ol_data_types_qp, only: wfun, hol
  use hel_bookkeeping_qp, only: helbookkeeping_ol_vert3
  use ol_vert_interface_qp, only: loop_AQ_Z
  use ol_parameters_init_qp, only: get_coupling

  implicit none
  integer(intkind1), intent(in)    :: ntry
  integer(intkind2), intent(inout) :: t(:,:), n(3)
  type(wfun),        intent(in)    :: J_Q(:)
  type(hol),         intent(inout) :: G_A,Gout_Z
  integer,           intent(in)    :: ng_RL
  complex(qp) :: G_add(size(Gout_Z%j,1),size(Gout_Z%j,2),size(Gout_Z%j,3))
  complex(qp) :: g_RL(2)
  integer(intkind2) :: h


  if (ntry == 1) call helbookkeeping_ol_vert3(ntry, J_Q, G_A, Gout_Z, n, t)
  if (.not. valid_hol(G_A, Gout_Z)) return

  Gout_Z%j = 0._qp
  g_RL = get_coupling(ng_RL)
  do h = 1, n(3)  ! recursion step
    call loop_AQ_Z(G_A%j(:,:,:,h),J_Q(t(1,h))%j,G_add,g_RL)
    Gout_Z%j(:,:,:,t(2,h)) = Gout_Z%j(:,:,:,t(2,h)) + G_add
  end do



end subroutine Hloop_AQ_Z


!***********************************************************************
subroutine Hloop_ZA_Q(ntry, G_Z, J_A, Gout_A, ng_RL, n, t)
!-----------------------------------------------------------------------
! bare ZA -> Q Z-like interaction
! **********************************************************************
  use kind_types, only: qp, intkind1, intkind2, qp
  use ol_data_types_qp, only: wfun, hol
  use hel_bookkeeping_qp, only: helbookkeeping_ol_vert3
  use ol_vert_interface_qp, only: loop_ZA_Q
  use ol_parameters_init_qp, only: get_coupling

  implicit none
  integer(intkind1), intent(in)    :: ntry
  integer(intkind2), intent(inout) :: t(:,:), n(3)
  type(wfun),        intent(in)    :: J_A(:)
  type(hol),         intent(inout) :: G_Z,Gout_A
  integer,           intent(in)    :: ng_RL
  complex(qp) :: G_add(size(Gout_A%j,1),size(Gout_A%j,2),size(Gout_A%j,3))
  complex(qp) :: g_RL(2)
  integer  :: h


  if (ntry == 1) call helbookkeeping_ol_vert3(ntry, J_A, G_Z, Gout_A, n, t)
  if (.not. valid_hol(G_Z, Gout_A)) return

  Gout_A%j = 0._qp
  g_RL = get_coupling(ng_RL)
  do h = 1, n(3)  ! recursion step
    call loop_ZA_Q(G_Z%j(:,:,:,h),J_A(t(1,h))%j,G_add,g_RL)
    Gout_A%j(:,:,:,t(2,h)) = Gout_A%j(:,:,:,t(2,h)) + G_add
  end do



end subroutine Hloop_ZA_Q


!***********************************************************************
subroutine Hloop_QZ_A(ntry, G_Q, J_Z, Gout_Q, ng_RL, n, t)
!-----------------------------------------------------------------------
! bare QZ -> A Z-like interaction
! **********************************************************************
  use kind_types, only: qp, intkind1, intkind2, qp
  use ol_data_types_qp, only: wfun, hol
  use hel_bookkeeping_qp, only: helbookkeeping_ol_vert3
  use ol_vert_interface_qp, only: loop_QZ_A
  use ol_parameters_init_qp, only: get_coupling

  implicit none
  integer(intkind1), intent(in)    :: ntry
  integer(intkind2), intent(inout) :: t(:,:), n(3)
  type(wfun),        intent(in)    :: J_Z(:)
  type(hol),         intent(inout) :: G_Q,Gout_Q
  integer,           intent(in)    :: ng_RL
  complex(qp) :: G_add(size(Gout_Q%j,1),size(Gout_Q%j,2),size(Gout_Q%j,3))
  complex(qp) :: g_RL(2)
  integer  :: h


  if (ntry == 1) call helbookkeeping_ol_vert3(ntry, J_Z, G_Q, Gout_Q, n, t)
  if (.not. valid_hol(G_Q, Gout_Q)) return

  Gout_Q%j = 0._qp
  g_RL = get_coupling(ng_RL)
  do h = 1, n(3)  ! recursion step
    call loop_QZ_A(G_Q%j(:,:,:,h),J_Z(t(1,h))%j,G_add,g_RL)
    Gout_Q%j(:,:,:,t(2,h)) = Gout_Q%j(:,:,:,t(2,h)) + G_add
  end do



end subroutine Hloop_QZ_A


!***********************************************************************
subroutine Hloop_QA_Z(ntry, G_Q, J_A, Gout_Z, ng_RL, n, t)
!-----------------------------------------------------------------------
! bare QA -> Z Z-like interaction
! **********************************************************************
  use kind_types, only: qp, intkind1, intkind2, qp
  use ol_data_types_qp, only: wfun, hol
  use hel_bookkeeping_qp, only: helbookkeeping_ol_vert3
  use ol_vert_interface_qp, only: loop_QA_Z
  use ol_parameters_init_qp, only: get_coupling

  implicit none
  integer(intkind1), intent(in)    :: ntry
  integer(intkind2), intent(inout) :: t(:,:), n(3)
  type(wfun),        intent(in)    :: J_A(:)
  type(hol),         intent(inout) :: G_Q,Gout_Z
  integer,           intent(in)    :: ng_RL
  complex(qp) :: G_add(size(Gout_Z%j,1),size(Gout_Z%j,2),size(Gout_Z%j,3))
  complex(qp) :: g_RL(2)
  integer(intkind2) :: h


  if (ntry == 1) call helbookkeeping_ol_vert3(ntry, J_A, G_Q, Gout_Z, n, t)
  if (.not. valid_hol(G_Q, Gout_Z)) return

  Gout_Z%j = 0._qp
  g_RL = get_coupling(ng_RL)
  do h = 1, n(3)  ! recursion step
    call loop_QA_Z(G_Q%j(:,:,:,h),J_A(t(1,h))%j,G_add,g_RL)
    Gout_Z%j(:,:,:,t(2,h)) = Gout_Z%j(:,:,:,t(2,h)) + G_add
  end do



end subroutine Hloop_QA_Z


!***********************************************************************
subroutine Hloop_ZQ_A(ntry, G_Z, J_Q, Gout_Q, ng_RL, n, t)
!-----------------------------------------------------------------------
! bare ZQ -> A Z-like interaction
! **********************************************************************
  use kind_types, only: qp, intkind1, intkind2, qp
  use ol_data_types_qp, only: wfun, hol
  use hel_bookkeeping_qp, only: helbookkeeping_ol_vert3
  use ol_vert_interface_qp, only: loop_ZQ_A
  use ol_parameters_init_qp, only: get_coupling

  implicit none
  integer(intkind1), intent(in)    :: ntry
  integer(intkind2), intent(inout) :: t(:,:), n(3)
  type(wfun),        intent(in)    :: J_Q(:)
  type(hol),         intent(inout) :: G_Z,Gout_Q
  integer,           intent(in)    :: ng_RL
  complex(qp) :: G_add(size(Gout_Q%j,1),size(Gout_Q%j,2),size(Gout_Q%j,3))
  complex(qp) :: g_RL(2)
  integer(intkind2) :: h


  if (ntry == 1) call helbookkeeping_ol_vert3(ntry, J_Q, G_Z, Gout_Q, n, t)
  if (.not. valid_hol(G_Z, Gout_Q)) return

  Gout_Q%j = 0._qp
  g_RL = get_coupling(ng_RL)
  do h = 1, n(3)  ! recursion step
    call loop_ZQ_A(G_Z%j(:,:,:,h), J_Q(t(1,h))%j,G_add,g_RL)
    Gout_Q%j(:,:,:,t(2,h)) = Gout_Q%j(:,:,:,t(2,h)) + G_add
  end do



end subroutine Hloop_ZQ_A


!***********************************************************************
subroutine Hloop_AW_Q(ntry, G_A, J_W, Gout_A, n, t)
!-----------------------------------------------------------------------
! bare AW -> A W-like interaction
! **********************************************************************
  use kind_types, only: qp, intkind1, intkind2, qp
  use ol_data_types_qp, only: wfun, hol
  use hel_bookkeeping_qp, only: helbookkeeping_ol_vert3
  use ol_vert_interface_qp, only: loop_AW_Q

  implicit none
  integer(intkind1), intent(in)    :: ntry
  integer(intkind2), intent(inout) :: t(:,:)
  integer(intkind2), intent(inout) :: n(3)
  type(wfun),  intent(in)    :: J_W(:)
  type(hol),   intent(inout) :: G_A,Gout_A
  complex(qp)  :: G_add(size(Gout_A%j,1),size(Gout_A%j,2),size(Gout_A%j,3))
  integer :: h


  if (ntry == 1) call helbookkeeping_ol_vert3(ntry, J_W, G_A, Gout_A, n, t)
  if (.not. valid_hol(G_A, Gout_A)) return

  Gout_A%j = 0._qp
  do h = 1, n(3)  ! recursion step
    call loop_AW_Q(G_A%j(:,:,:,h), J_W(t(1,h))%j, G_add)
    Gout_A%j(:,:,:,t(2,h)) = Gout_A%j(:,:,:,t(2,h)) + G_add
  end do



end subroutine Hloop_AW_Q


!***********************************************************************
subroutine Hloop_AQ_W(ntry, G_A, J_Q, Gout_W, n, t)
!-----------------------------------------------------------------------
! bare AQ -> W W-like interaction
! **********************************************************************
  use kind_types, only: qp, intkind1, intkind2, qp
  use ol_data_types_qp, only: wfun, hol
  use hel_bookkeeping_qp, only: helbookkeeping_ol_vert3
  use ol_vert_interface_qp, only: loop_AQ_W

  implicit none
  integer(intkind1), intent(in)    :: ntry
  integer(intkind2), intent(inout) :: t(:,:)
  integer(intkind2), intent(inout) :: n(3)
  type(wfun),  intent(in)    :: J_Q(:)
  type(hol),   intent(inout) :: G_A,Gout_W
  complex(qp)  :: G_add(size(Gout_W%j,1),size(Gout_W%j,2),size(Gout_W%j,3))
  integer :: h


  if (ntry == 1) call helbookkeeping_ol_vert3(ntry, J_Q, G_A, Gout_W, n, t)
  if (.not. valid_hol(G_A, Gout_W)) return

  Gout_W%j = 0._qp
  do h = 1, n(3)  ! recursion step
    call loop_AQ_W(G_A%j(:,:,:,h), J_Q(t(1,h))%j, G_add)
    Gout_W%j(:,:,:,t(2,h)) = Gout_W%j(:,:,:,t(2,h)) + G_add
  end do



end subroutine Hloop_AQ_W


!***********************************************************************
subroutine Hloop_WA_Q(ntry, G_W, J_A, Gout_A, n, t)
!-----------------------------------------------------------------------
! bare WA -> Q W-like interaction
! **********************************************************************
  use kind_types, only: qp, intkind1, intkind2, qp
  use ol_data_types_qp, only: wfun, hol
  use hel_bookkeeping_qp, only: helbookkeeping_ol_vert3
  use ol_vert_interface_qp, only: loop_WA_Q

  implicit none
  integer(intkind1), intent(in)    :: ntry
  integer(intkind2), intent(inout) :: t(:,:)
  integer(intkind2), intent(inout) :: n(3)
  type(wfun),  intent(in)    :: J_A(:)
  type(hol),   intent(inout) :: G_W,Gout_A
  complex(qp)  :: G_add(size(Gout_A%j,1),size(Gout_A%j,2),size(Gout_A%j,3))
  integer :: h


  if (ntry == 1) call helbookkeeping_ol_vert3(ntry, J_A, G_W, Gout_A, n, t)
  if (.not. valid_hol(G_W, Gout_A)) return

  Gout_A%j = 0._qp
  do h = 1, n(3)  ! recursion step
    call loop_WA_Q(G_W%j(:,:,:,h), J_A(t(1,h))%j, G_add)
    Gout_A%j(:,:,:,t(2,h)) = Gout_A%j(:,:,:,t(2,h)) + G_add
  end do



end subroutine Hloop_WA_Q


!***********************************************************************
subroutine Hloop_QW_A(ntry, G_Q, J_W, Gout_Q, n, t)
!-----------------------------------------------------------------------
! bare QW -> A W-like interaction
! **********************************************************************
  use kind_types, only: qp, intkind1, intkind2, qp
  use ol_data_types_qp, only: wfun, hol
  use hel_bookkeeping_qp, only: helbookkeeping_ol_vert3
  use ol_vert_interface_qp, only: loop_QW_A

  implicit none
  integer(intkind1), intent(in)    :: ntry
  integer(intkind2), intent(inout) :: t(:,:)
  integer(intkind2), intent(inout) :: n(3)
  type(wfun),  intent(in)    :: J_W(:)
  type(hol),   intent(inout) :: G_Q,Gout_Q
  complex(qp)  :: G_add(size(Gout_Q%j,1),size(Gout_Q%j,2),size(Gout_Q%j,3))
  integer :: h


  if (ntry == 1) call helbookkeeping_ol_vert3(ntry, J_W, G_Q, Gout_Q, n, t)
  if (.not. valid_hol(G_Q, Gout_Q)) return

  Gout_Q%j = 0._qp
  do h = 1, n(3)  ! recursion step
    call loop_QW_A(G_Q%j(:,:,:,h), J_W(t(1,h))%j, G_add)
    Gout_Q%j(:,:,:,t(2,h)) = Gout_Q%j(:,:,:,t(2,h)) + G_add
  end do



end subroutine Hloop_QW_A


!***********************************************************************
subroutine Hloop_QA_W(ntry, G_Q, J_A, Gout_W, n, t)
!-----------------------------------------------------------------------
! bare QA -> W W-like interaction
! **********************************************************************
  use kind_types, only: qp, intkind1, intkind2, qp
  use ol_data_types_qp, only: wfun, hol
  use hel_bookkeeping_qp, only: helbookkeeping_ol_vert3
  use ol_vert_interface_qp, only: loop_QA_W

  implicit none
  integer(intkind1), intent(in)    :: ntry
  integer(intkind2), intent(inout) :: t(:,:)
  integer(intkind2), intent(inout) :: n(3)
  type(wfun),  intent(in)    :: J_A(:)
  type(hol),   intent(inout) :: G_Q,Gout_W
  complex(qp)  :: G_add(size(Gout_W%j,1),size(Gout_W%j,2),size(Gout_W%j,3))
  integer :: h


  if (ntry == 1) call helbookkeeping_ol_vert3(ntry, J_A, G_Q, Gout_W, n, t)
  if (.not. valid_hol(G_Q, Gout_W)) return

  Gout_W%j = 0._qp
  do h = 1, n(3)  ! recursion step
    call loop_QA_W(G_Q%j(:,:,:,h), J_A(t(1,h))%j, G_add)
    Gout_W%j(:,:,:,t(2,h)) = Gout_W%j(:,:,:,t(2,h)) + G_add
  end do



end subroutine Hloop_QA_W


!***********************************************************************
subroutine Hloop_WQ_A(ntry, G_W, J_Q, Gout_Q, n, t)
!-----------------------------------------------------------------------
! bare WQ -> A W-like interaction
! **********************************************************************
  use kind_types, only: qp, intkind1, intkind2, qp
  use ol_data_types_qp, only: wfun, hol
  use hel_bookkeeping_qp, only: helbookkeeping_ol_vert3
  use ol_vert_interface_qp, only: loop_WQ_A

  implicit none
  integer(intkind1), intent(in)    :: ntry
  integer(intkind2), intent(inout) :: t(:,:)
  integer(intkind2), intent(inout) :: n(3)
  type(wfun),  intent(in)    :: J_Q(:)
  type(hol),   intent(inout) :: G_W,Gout_Q
  complex(qp)  :: G_add(size(Gout_Q%j,1),size(Gout_Q%j,2),size(Gout_Q%j,3))
  integer :: h


  if (ntry == 1) call helbookkeeping_ol_vert3(ntry, J_Q, G_W, Gout_Q, n, t)
  if (.not. valid_hol(G_W, Gout_Q)) return

  Gout_Q%j = 0._qp
  do h = 1, n(3)  ! recursion step
    call loop_WQ_A(G_W%j(:,:,:,h), J_Q(t(1,h))%j, G_add)
    Gout_Q%j(:,:,:,t(2,h)) = Gout_Q%j(:,:,:,t(2,h)) + G_add
  end do



end subroutine Hloop_WQ_A


!***********************************************************************
subroutine Hloop_AV_Q(ntry, G_A, J_V, Gout_A, n, t)
!-----------------------------------------------------------------------
! bare AV -> Q gluon-like interaction
! **********************************************************************
  use kind_types, only: qp, intkind1, intkind2, qp
  use ol_data_types_qp, only: wfun, hol
  use hel_bookkeeping_qp, only: helbookkeeping_ol_vert3
  use ol_vert_interface_qp, only: loop_AV_Q

  implicit none
  integer(intkind1), intent(in)    :: ntry
  integer(intkind2), intent(inout) :: t(:,:)
  integer(intkind2), intent(inout) :: n(3)
  type(wfun),  intent(in)    :: J_V(:)
  type(hol),   intent(inout) :: G_A,Gout_A
  complex(qp)  :: G_add(size(Gout_A%j,1),size(Gout_A%j,2),size(Gout_A%j,3))
  integer :: h


  if (ntry == 1) call helbookkeeping_ol_vert3(ntry, J_V, G_A, Gout_A, n, t)  ! STILL TO BE ADAPTED
  if (.not. valid_hol(G_A, Gout_A)) return

  Gout_A%j = 0._qp
  do h = 1, n(3)  ! recursion step
    call loop_AV_Q(G_A%j(:,:,:,h), J_V(t(1,h))%j, G_add)
    Gout_A%j(:,:,:,t(2,h)) = Gout_A%j(:,:,:,t(2,h)) + G_add
  end do



end subroutine Hloop_AV_Q

!***********************************************************************
subroutine Hloop_AQ_V(ntry, G_A, J_Q, Gout_V, n, t)
!-----------------------------------------------------------------------
! bare AQ -> V gluon-like interaction
! **********************************************************************
  use kind_types, only: qp, intkind1, intkind2, qp
  use ol_data_types_qp, only: wfun, hol
  use hel_bookkeeping_qp, only: helbookkeeping_ol_vert3
  use ol_vert_interface_qp, only: loop_AQ_V

  implicit none
  integer(intkind1), intent(in)    :: ntry
  integer(intkind2), intent(inout) :: t(:,:)
  integer(intkind2), intent(inout) :: n(3)
  type(wfun),  intent(in)    :: J_Q(:)
  type(hol),   intent(inout) :: G_A,Gout_V
  complex(qp)  :: G_add(size(Gout_V%j,1),size(Gout_V%j,2),size(Gout_V%j,3))
  integer  :: h


  if (ntry == 1) call helbookkeeping_ol_vert3(ntry, J_Q, G_A, Gout_V, n, t)
  if (.not. valid_hol(G_A, Gout_V)) return

  Gout_V%j = 0._qp
  do h = 1, n(3)  ! recursion step
    call loop_AQ_V(G_A%j(:,:,:,h), J_Q(t(1,h))%j, G_add)
    Gout_V%j(:,:,:,t(2,h)) = Gout_V%j(:,:,:,t(2,h)) + G_add
  end do



end subroutine Hloop_AQ_V


!***********************************************************************
subroutine Hloop_VA_Q(ntry, G_V, J_A, Gout_A, n, t)
!-----------------------------------------------------------------------
! bare VA -> Q gluon-like interaction
! **********************************************************************
  use kind_types, only: qp, intkind1, intkind2, qp
  use ol_data_types_qp, only: wfun, hol
  use hel_bookkeeping_qp, only: helbookkeeping_ol_vert3
  use ol_vert_interface_qp, only: loop_VA_Q

  implicit none
  integer(intkind1), intent(in)    :: ntry
  integer(intkind2), intent(inout) :: t(:,:)
  integer(intkind2), intent(inout) :: n(3)
  type(wfun),  intent(in)    :: J_A(:)
  type(hol),   intent(inout) :: G_V,Gout_A
  complex(qp)  :: G_add(size(Gout_A%j,1),size(Gout_A%j,2),size(Gout_A%j,3))
  integer :: h


  if (ntry == 1) call helbookkeeping_ol_vert3(ntry, J_A, G_V, Gout_A, n, t)
  if (.not. valid_hol(G_V, Gout_A)) return

  Gout_A%j = 0._qp
  do h = 1, n(3)  ! recursion step + helicity summation
    call loop_VA_Q(G_V%j(:,:,:,h), J_A(t(1,h))%j, G_add)
    Gout_A%j(:,:,:,t(2,h)) = Gout_A%j(:,:,:,t(2,h)) + G_add
  end do



end subroutine Hloop_VA_Q


!***********************************************************************
subroutine Hloop_QV_A(ntry, G_Q, J_V, Gout_Q, n, t)
!-----------------------------------------------------------------------
! bare QV -> A gluon-like interaction
! **********************************************************************
  use kind_types, only: qp, intkind1, intkind2, qp
  use ol_data_types_qp, only: wfun, hol
  use hel_bookkeeping_qp, only: helbookkeeping_ol_vert3
  use ol_vert_interface_qp, only: loop_QV_A

  implicit none
  integer(intkind1), intent(in)    :: ntry
  integer(intkind2), intent(inout) :: t(:,:)
  integer(intkind2), intent(inout) :: n(3)
  type(wfun),  intent(in)    :: J_V(:)
  type(hol),   intent(inout) :: G_Q,Gout_Q
  complex(qp)  :: G_add(size(Gout_Q%j,1),size(Gout_Q%j,2),size(Gout_Q%j,3))
  integer  :: h


  if (ntry == 1) call helbookkeeping_ol_vert3(ntry, J_V, G_Q, Gout_Q, n, t)
  if (.not. valid_hol(G_Q, Gout_Q)) return

  Gout_Q%j = 0._qp
  do h = 1, n(3)  ! recursion step
    call loop_QV_A(G_Q%j(:,:,:,h), J_V(t(1,h))%j, G_add)
    Gout_Q%j(:,:,:,t(2,h)) = Gout_Q%j(:,:,:,t(2,h)) + G_add
  end do



end subroutine Hloop_QV_A


!***********************************************************************
subroutine Hloop_QA_V(ntry, G_Q, J_A, Gout_V, n, t)
!-----------------------------------------------------------------------
! bare QA -> V gluon-like interaction
! **********************************************************************
  use kind_types, only: qp, intkind1, intkind2, qp
  use ol_data_types_qp, only: wfun, hol
  use hel_bookkeeping_qp, only: helbookkeeping_ol_vert3
  use ol_vert_interface_qp, only: loop_QA_V

  implicit none
  integer(intkind1), intent(in)    :: ntry
  integer(intkind2), intent(inout) :: t(:,:)
  integer(intkind2), intent(inout) :: n(3)
  type(wfun),  intent(in)    :: J_A(:)
  type(hol),   intent(inout) :: G_Q,Gout_V
  complex(qp)  :: G_add(size(Gout_V%j,1),size(Gout_V%j,2),size(Gout_V%j,3))
  integer :: h


  if (ntry == 1) call helbookkeeping_ol_vert3(ntry, J_A, G_Q, Gout_V, n, t)  ! STILL TO BE ADAPTED
  if (.not. valid_hol(G_Q, Gout_V)) return

  Gout_V%j = 0._qp
  do h = 1, n(3)  ! recursion step
    call loop_QA_V(G_Q%j(:,:,:,h), J_A(t(1,h))%j, G_add)
    Gout_V%j(:,:,:,t(2,h)) = Gout_V%j(:,:,:,t(2,h)) + G_add
  end do



end subroutine Hloop_QA_V


!***********************************************************************
subroutine Hloop_VQ_A(ntry, G_V, J_Q, Gout_A, n, t)
!-----------------------------------------------------------------------
! bare VQ -> A gluon-like interaction
! **********************************************************************
  use kind_types, only: qp, intkind1, intkind2, qp
  use ol_data_types_qp, only: wfun, hol
  use hel_bookkeeping_qp, only: helbookkeeping_ol_vert3
  use ol_vert_interface_qp, only: loop_VQ_A

  implicit none
  integer(intkind1), intent(in)    :: ntry
  integer(intkind2),     intent(inout) :: t(:,:)
  integer(intkind2), intent(inout) :: n(3)
  type(wfun),  intent(in)    :: J_Q(:)
  type(hol),   intent(inout) :: G_V,Gout_A
  complex(qp)  :: G_add(size(Gout_A%j,1),size(Gout_A%j,2),size(Gout_A%j,3))
  integer :: h


  if (ntry == 1) call helbookkeeping_ol_vert3(ntry, J_Q, G_V, Gout_A, n, t)  ! STILL TO BE ADAPTED
  if (.not. valid_hol(G_V, Gout_A)) return

  Gout_A%j = 0._qp
  do h = 1, n(3)  ! recursion step
    call loop_VQ_A(G_V%j(:,:,:,h), J_Q(t(1,h))%j, G_add)
    Gout_A%j(:,:,:,t(2,h)) = Gout_A%j(:,:,:,t(2,h)) + G_add
  end do



end subroutine Hloop_VQ_A


!***********************************************************************
subroutine Hloop_UV_W(ntry, Gin_V, moml, J_V, momt, Gout_V, n, t)
!-----------------------------------------------------------------------
! bare VV -> V vertex
! **********************************************************************
  use kind_types, only: qp, intkind1, intkind2, qp
  use ol_data_types_qp, only: wfun, hol
  use hel_bookkeeping_qp, only: helbookkeeping_ol_vert3
  use ol_vert_interface_qp, only: loop_UV_W
  use ol_kinematics_qp, only: get_LC_4

  implicit none
  integer(intkind1), intent(in)    :: ntry
  integer,           intent(in)    :: moml,momt
  integer(intkind2), intent(inout) :: t(:,:)
  integer(intkind2), intent(inout) :: n(3)
  type(wfun),  intent(in)    :: J_V(:)
  type(hol),   intent(inout) :: Gin_V,Gout_V
  complex(qp)  :: G_add(size(Gout_V%j,1),size(Gout_V%j,2),size(Gout_V%j,3))
  integer  :: h


  if (ntry == 1) call helbookkeeping_ol_vert3(ntry, J_V, Gin_V, Gout_V, n, t)  ! STILL TO BE ADAPTED
  if (.not. valid_hol(Gin_V, Gout_V)) return

  Gout_V%j = 0._qp
  do h = 1, n(3)  ! recursion step
    call loop_UV_W(Gin_V%j(:,:,:,h), get_LC_4(moml), J_V(t(1,h))%j, get_LC_4(momt), G_add)
    Gout_V%j(:,:,:,t(2,h)) = Gout_V%j(:,:,:,t(2,h)) + G_add
  end do



end subroutine Hloop_UV_W


!***********************************************************************
subroutine Hloop_UW_V(ntry, Gin_V, moml, J_V, momt, Gout_V, n, t)
!-----------------------------------------------------------------------
! bare VV -> V vertex
! **********************************************************************
  use kind_types, only: qp, intkind1, intkind2, qp
  use ol_data_types_qp, only: wfun, hol
  use hel_bookkeeping_qp, only: helbookkeeping_ol_vert3
  use ol_vert_interface_qp, only: loop_UW_V
  use ol_kinematics_qp, only: get_LC_4

  implicit none
  integer(intkind1), intent(in)    :: ntry
  integer,           intent(in)    :: moml,momt
  integer(intkind2),     intent(inout) :: t(:,:)
  integer(intkind2), intent(inout) :: n(3)
  type(wfun),  intent(in)    :: J_V(:)
  type(hol),   intent(inout) :: Gin_V,Gout_V
  complex(qp)  :: G_add(size(Gout_V%j,1),size(Gout_V%j,2),size(Gout_V%j,3))
  integer  :: h


  if (ntry == 1) call helbookkeeping_ol_vert3(ntry, J_V, Gin_V, Gout_V, n, t)  ! STILL TO BE ADAPTED
  if (.not. valid_hol(Gin_V, Gout_V)) return

  Gout_V%j = 0._qp
  do h = 1, n(3)  ! recursion step
    call loop_UW_V(Gin_V%j(:,:,:,h), get_LC_4(moml), J_V(t(1,h))%j, get_LC_4(momt), G_add)
    Gout_V%j(:,:,:,t(2,h)) = Gout_V%j(:,:,:,t(2,h)) + G_add
  end do



end subroutine Hloop_UW_V


!***********************************************************************
subroutine Hloop_EV_V(ntry, Gin_V, J1, J2, Gout_V, n, t)
!-----------------------------------------------------------------------
! bare 4-gluon vertex when the sigma particle is in the loop
! **********************************************************************
  use kind_types, only: qp, qp, intkind1, intkind2
  use ol_data_types_qp, only: wfun, hol
  use hel_bookkeeping_qp, only: helbookkeeping_ol_vert4
  use ol_vert_interface_qp, only: loop_EV_V

  implicit none
  integer(intkind1), intent(in)    :: ntry
  integer(intkind2), intent(inout) :: t(:,:)
  integer(intkind2), intent(inout) :: n(4)
  type(wfun),        intent(in)    :: J1(:), J2(:)
  type(hol),         intent(inout) :: Gin_V,Gout_V
  complex(qp)  :: G_add(size(Gout_V%j,1),size(Gout_V%j,2),size(Gout_V%j,3))
  integer  :: h


  if (ntry == 1) call helbookkeeping_ol_vert4(ntry, J1, J2, Gin_V, Gout_V, n, t)  ! STILL TO BE ADAPTED
  if (.not. valid_hol(Gin_V, Gout_V)) return

  Gout_V%j = 0._qp
  G_add    = 0._qp

  do h = 1, n(4)  ! recursion step
    call loop_EV_V(Gin_V%j(:,:,:,h), J1(t(1,h))%j, J2(t(2,h))%j, G_add)
    Gout_V%j(:,:,:,t(3,h)) = Gout_V%j(:,:,:,t(3,h)) + G_add
  end do



end subroutine Hloop_EV_V

!***********************************************************************
subroutine Hloop_VE_V(ntry, Gin_V, J1, J2, Gout_V, n, t)
!-----------------------------------------------------------------------
! bare 4-gluon vertex when the sigma particle enters the loop as a tree wf
! **********************************************************************
  use kind_types, only: qp, qp, intkind1, intkind2
  use ol_data_types_qp, only: wfun, hol
  use hel_bookkeeping_qp, only: helbookkeeping_ol_vert4
  use ol_vert_interface_qp, only: loop_VE_V

  implicit none
  integer(intkind1), intent(in)    :: ntry
  integer(intkind2), intent(inout) :: t(:,:)
  integer(intkind2), intent(inout) :: n(4)
  type(wfun),        intent(in)    :: J1(:), J2(:)
  type(hol),         intent(inout) :: Gin_V,Gout_V
  complex(qp)  :: G_add(size(Gout_V%j,1),size(Gout_V%j,2),size(Gout_V%j,3))
  integer  :: h


  if (ntry == 1) call helbookkeeping_ol_vert4(ntry, J1, J2, Gin_V, Gout_V, n, t)  ! STILL TO BE ADAPTED
  if (.not. valid_hol(Gin_V, Gout_V)) return

  Gout_V%j = 0._qp
  G_add    = 0._qp

  do h = 1, n(4)  ! recursion step
    call loop_VE_V(Gin_V%j(:,:,:,h), J1(t(1,h))%j, J2(t(2,h))%j, G_add)
    Gout_V%j(:,:,:,t(3,h)) = Gout_V%j(:,:,:,t(3,h)) + G_add
  end do



end subroutine Hloop_VE_V

!***********************************************************************
subroutine Hloop_GGG_G_23(ntry, Gin_V, J1, J2, Gout_V, n, t)
!-----------------------------------------------------------------------
! 4-gluon vertex interaction
! **********************************************************************
  use kind_types, only: qp, intkind1, intkind2, qp
  use ol_data_types_qp, only: wfun, hol
  use hel_bookkeeping_qp, only: helbookkeeping_ol_vert4
  use ol_vert_interface_qp, only: loop_GGG_G_23

  implicit none
  integer(intkind1), intent(in)    :: ntry
  integer(intkind2), intent(inout) :: t(:,:)
  integer(intkind2), intent(inout) :: n(4)
  type(wfun),  intent(in)    :: J1(:), J2(:)
  type(hol),   intent(inout) :: Gin_V,Gout_V
  complex(qp)  :: G_add(size(Gout_V%j,1),size(Gout_V%j,2),size(Gout_V%j,3))
  integer  :: h


  if (ntry == 1) call helbookkeeping_ol_vert4(ntry, J1, J2, Gin_V, Gout_V, n, t)  ! STILL TO BE ADAPTED
  if (.not. valid_hol(Gin_V, Gout_V)) return

  Gout_V%j = 0._qp
  do h = 1, n(4)  ! recursion step
    call loop_GGG_G_23(Gin_V%j(:,:,:,h), J1(t(1,h))%j, J2(t(2,h))%j, G_add)
    Gout_V%j(:,:,:,t(3,h)) = Gout_V%j(:,:,:,t(3,h)) + G_add
  end do



end subroutine Hloop_GGG_G_23


!***********************************************************************
subroutine Hloop_GGG_G_12(ntry, Gin_V, J1, J2, Gout_V, n, t)
!-----------------------------------------------------------------------
! 4-gluon vertex interaction
! **********************************************************************
  use kind_types, only: qp, intkind1, intkind2, qp
  use ol_data_types_qp, only: wfun, hol
  use hel_bookkeeping_qp, only: helbookkeeping_ol_vert4
  use ol_vert_interface_qp, only: loop_GGG_G_12

  implicit none
  integer(intkind1), intent(in)    :: ntry
  integer(intkind2), intent(inout) :: t(:,:)
  integer(intkind2), intent(inout) :: n(4)
  type(wfun),  intent(in)    :: J1(:), J2(:)
  type(hol),   intent(inout) :: Gin_V,Gout_V
  complex(qp)  :: G_add(size(Gout_V%j,1),size(Gout_V%j,2),size(Gout_V%j,3))
  integer  :: h


  if (ntry == 1) call helbookkeeping_ol_vert4(ntry, J1, J2, Gin_V, Gout_V, n, t)  ! STILL TO BE ADAPTED
  if (.not. valid_hol(Gin_V, Gout_V)) return

  Gout_V%j = 0._qp
  do h = 1, n(4)  ! recursion step
    call loop_GGG_G_12(Gin_V%j(:,:,:,h), J1(t(1,h))%j, J2(t(2,h))%j, G_add)
    Gout_V%j(:,:,:,t(3,h)) = Gout_V%j(:,:,:,t(3,h)) + G_add
  end do



end subroutine Hloop_GGG_G_12


!***********************************************************************
subroutine Hloop_CV_D(ntry, Gin_C, moml, J_V, momt, Gout_C, n, t)
!-----------------------------------------------------------------------
! bare ghost gluon -> ghost interaction
! always comes in a closed ghost loop
! **********************************************************************
  use kind_types, only: qp, intkind1, intkind2, qp
  use ol_data_types_qp, only: wfun, hol
  use hel_bookkeeping_qp, only: helbookkeeping_ol_vert3
  use ol_vert_interface_qp, only: loop_CV_D
  use ol_kinematics_qp, only: get_LC_4

  implicit none
  integer(intkind1), intent(in)    :: ntry
  integer,           intent(in)    :: moml,momt
  integer(intkind2), intent(inout) :: t(:,:)
  integer(intkind2), intent(inout) :: n(3)
  type(wfun),  intent(in)    :: J_V(:)
  type(hol),   intent(inout) :: Gin_C,Gout_C
  complex(qp)  :: G_add(size(Gout_C%j,1),size(Gout_C%j,2),size(Gout_C%j,3))
  integer  :: h


  if (ntry == 1) call helbookkeeping_ol_vert3(ntry, J_V, Gin_C, Gout_C, n, t)  ! STILL TO BE ADAPTED
  if (.not. valid_hol(Gin_C, Gout_C)) return

  Gout_C%j = 0._qp
  do h = 1, n(3)  ! recursion step
    call loop_CV_D(Gin_C%j(:,:,:,h), get_LC_4(moml), J_V(t(1,h))%j, get_LC_4(momt), G_add)
    Gout_C%j(:,:,:,t(2,h)) = Gout_C%j(:,:,:,t(2,h)) + G_add
  end do



end subroutine Hloop_CV_D


!***********************************************************************
subroutine Hloop_DV_C(ntry, Gin_D, moml, J_V, Gout_D, n, t)
!-----------------------------------------------------------------------
! bare anti-ghost gluon -> anti-ghost interaction
! always comes in a closed ghost loop
! **********************************************************************
  use kind_types, only: qp, intkind1, intkind2, qp
  use ol_data_types_qp, only: wfun, hol
  use hel_bookkeeping_qp, only: helbookkeeping_ol_vert3
  use ol_vert_interface_qp, only: loop_DV_C
  use ol_kinematics_qp, only: get_LC_4

  implicit none
  integer(intkind1), intent(in)    :: ntry
  integer,           intent(in)    :: moml
  integer(intkind2), intent(inout) :: t(:,:)
  integer(intkind2), intent(inout) :: n(3)
  type(wfun),  intent(in)    :: J_V(:)
  type(hol),   intent(inout) :: Gin_D,Gout_D
  complex(qp)  :: G_add(size(Gout_D%j,1),size(Gout_D%j,2),size(Gout_D%j,3))
  integer :: h


  if (ntry == 1) call helbookkeeping_ol_vert3(ntry, J_V, Gin_D, Gout_D, n, t)
  if (.not. valid_hol(Gin_D, Gout_D)) return

  Gout_D%j = 0._qp
  do h = 1, n(3)  ! recursion step
    call loop_DV_C(Gin_D%j(:,:,:,h), get_LC_4(moml), J_V(t(1,h))%j, G_add)
    Gout_D%j(:,:,:,t(2,h)) = Gout_D%j(:,:,:,t(2,h)) + G_add
  end do



end subroutine Hloop_DV_C


!***********************************************************************
subroutine Hloop_AS_Q(ntry, G_A, J_S, Gout_A, ng_RL, n, t)
!-----------------------------------------------------------------------
! bare anti-quark scalar -> quark interaction
!***********************************************************************
  use kind_types, only: qp, intkind1, intkind2, qp
  use ol_data_types_qp, only: wfun, hol
  use hel_bookkeeping_qp, only: helbookkeeping_ol_vert3
  use ol_vert_interface_qp, only: loop_AS_Q
  use ol_parameters_init_qp, only: get_coupling

  implicit none
  integer(intkind1), intent(in)    :: ntry
  integer(intkind2), intent(inout) :: n(3), t(:,:)
  integer,           intent(in)    :: ng_RL
  type(wfun),        intent(in)    :: J_S(:)
  type(hol),         intent(inout) :: G_A,Gout_A
  complex(qp) :: G_add(size(Gout_A%j,1),size(Gout_A%j,2),size(Gout_A%j,3))
  complex(qp) :: g_RL(2)
  integer(intkind2) :: h


  if (ntry == 1) call helbookkeeping_ol_vert3(ntry, J_S, G_A, Gout_A, n, t)
  if (.not. valid_hol(G_A, Gout_A)) return

  Gout_A%j = 0._qp
  g_RL = get_coupling(ng_RL)
  do h = 1, n(3)  ! recursion step
    call loop_AS_Q(G_A%j(:,:,:,h), J_S(t(1,h))%j,G_add,g_RL)
    Gout_A%j(:,:,:,t(2,h)) = Gout_A%j(:,:,:,t(2,h)) + G_add
  end do



end subroutine Hloop_AS_Q


!***********************************************************************
subroutine Hloop_SA_Q(ntry, G_S, J_A, Gout_A, ng_RL, n, t)
!-----------------------------------------------------------------------
! bare scalar anti-quark -> quark interaction
!***********************************************************************
  use kind_types, only: qp, intkind1, intkind2, qp
  use ol_data_types_qp, only: wfun, hol
  use hel_bookkeeping_qp, only: helbookkeeping_ol_vert3
  use ol_vert_interface_qp, only: loop_SA_Q
  use ol_parameters_init_qp, only: get_coupling

  implicit none
  integer(intkind1), intent(in)    :: ntry
  integer(intkind2), intent(inout) :: n(3), t(:,:)
  integer,           intent(in)    :: ng_RL
  type(wfun),        intent(in)    :: J_A(:)
  type(hol),         intent(inout) :: G_S,Gout_A
  complex(qp) :: G_add(size(Gout_A%j,1),size(Gout_A%j,2),size(Gout_A%j,3))
  complex(qp) :: g_RL(2)
  integer(intkind2) :: h


  if (ntry == 1) call helbookkeeping_ol_vert3(ntry, J_A, G_S, Gout_A, n, t)
  if (.not. valid_hol(G_S, Gout_A)) return

  Gout_A%j = 0._qp
  g_RL = get_coupling(ng_RL)
  do h = 1, n(3)  ! recursion step
    call loop_SA_Q(G_S%j(:,:,:,h),J_A(t(1,h))%j,G_add,g_RL)
    Gout_A%j(:,:,:,t(2,h)) = Gout_A%j(:,:,:,t(2,h)) + G_add
  end do



end subroutine Hloop_SA_Q


!***********************************************************************
subroutine Hloop_QS_A(ntry, G_Q, J_S, Gout_Q, ng_RL, n, t)
!-----------------------------------------------------------------------
! bare quark scalar -> anti-quark interaction
! ----------------------------------------------------------------------
  use kind_types, only: qp, intkind1, intkind2, qp
  use ol_data_types_qp, only: wfun, hol
  use hel_bookkeeping_qp, only: helbookkeeping_ol_vert3
  use ol_vert_interface_qp, only: loop_QS_A
  use ol_parameters_init_qp, only: get_coupling

  implicit none
  integer(intkind1), intent(in)    :: ntry
  integer(intkind2), intent(inout) :: n(3), t(:,:)
  integer,           intent(in)    :: ng_RL
  type(wfun),        intent(in)    :: J_S(:)
  type(hol),         intent(inout) :: G_Q,Gout_Q
  complex(qp) :: G_add(size(Gout_Q%j,1),size(Gout_Q%j,2),size(Gout_Q%j,3))
  complex(qp) :: g_RL(2)
  integer  :: h


  if (ntry == 1) call helbookkeeping_ol_vert3(ntry, J_S, G_Q, Gout_Q, n, t)
  if (.not. valid_hol(G_Q, Gout_Q)) return

  Gout_Q%j = 0._qp
  g_RL = get_coupling(ng_RL)
  do h = 1, n(3)  ! recursion step
    call loop_QS_A(G_Q%j(:,:,:,h),J_S(t(1,h))%j,G_add,g_RL)
    Gout_Q%j(:,:,:,t(2,h)) = Gout_Q%j(:,:,:,t(2,h)) + G_add
  end do



end subroutine Hloop_QS_A


!***********************************************************************
subroutine Hloop_SQ_A(ntry, G_S, J_Q, Gout_Q, ng_RL, n, t)
!-----------------------------------------------------------------------
! bare scalar quark -> anti-quark interaction
! ----------------------------------------------------------------------
  use kind_types, only: qp, intkind1, intkind2, qp
  use ol_data_types_qp, only: wfun, hol
  use hel_bookkeeping_qp, only: helbookkeeping_ol_vert3
  use ol_vert_interface_qp, only: loop_SQ_A
  use ol_parameters_init_qp, only: get_coupling

  implicit none
  integer(intkind1), intent(in)    :: ntry
  integer(intkind2), intent(inout) :: n(3), t(:,:)
  integer,           intent(in)    :: ng_RL
  type(wfun),        intent(in)    :: J_Q(:)
  type(hol),         intent(inout) :: G_S,Gout_Q
  complex(qp) :: G_add(size(Gout_Q%j,1),size(Gout_Q%j,2),size(Gout_Q%j,3))
  complex(qp) :: g_RL(2)
  integer  :: h


  if (ntry == 1) call helbookkeeping_ol_vert3(ntry, J_Q, G_S, Gout_Q, n, t)
  if (.not. valid_hol(G_S, Gout_Q)) return

  Gout_Q%j = 0._qp
  g_RL = get_coupling(ng_RL)
  do h = 1, n(3)  ! recursion step
    call loop_SQ_A(G_S%j(:,:,:,h), J_Q(t(1,h))%j,G_add,g_RL)
    Gout_Q%j(:,:,:,t(2,h)) = Gout_Q%j(:,:,:,t(2,h)) + G_add
  end do



end subroutine Hloop_SQ_A


!***********************************************************************
subroutine Hloop_QA_S(ntry, G_Q, J_A, Gout_S, ng_RL, n, t)
!-----------------------------------------------------------------------
! bare quar anti-quark -> scalar interaction
!***********************************************************************
  use kind_types, only: qp, intkind1, intkind2, qp
  use ol_data_types_qp, only: wfun, hol
  use hel_bookkeeping_qp, only: helbookkeeping_ol_vert3
  use ol_vert_interface_qp, only: loop_QA_S
  use ol_parameters_init_qp, only: get_coupling

  implicit none
  integer(intkind1), intent(in)    :: ntry
  integer(intkind2), intent(inout) :: n(3), t(:,:)
  integer,           intent(in)    :: ng_RL
  type(wfun),        intent(in)    :: J_A(:)
  type(hol),         intent(inout) :: G_Q,Gout_S
  complex(qp) :: G_add(size(Gout_S%j,1),size(Gout_S%j,2),size(Gout_S%j,3))
  complex(qp) :: g_RL(2)
  integer(intkind2) :: h


  if (ntry == 1) call helbookkeeping_ol_vert3(ntry, J_A, G_Q, Gout_S, n, t)
  if (.not. valid_hol(G_Q, Gout_S)) return

  Gout_S%j = 0._qp
  g_RL = get_coupling(ng_RL)
  do h = 1, n(3)  ! recursion step
    call loop_QA_S(G_Q%j(:,:,:,h),J_A(t(1,h))%j,G_add,g_RL)
    Gout_S%j(:,:,:,t(2,h)) = Gout_S%j(:,:,:,t(2,h)) + G_add
  end do



end subroutine Hloop_QA_S


!***********************************************************************
subroutine Hloop_AQ_S(ntry, G_A, J_Q, Gout_S, ng_RL, n, t)
!-----------------------------------------------------------------------
! bare quar anti-quark -> scalar interaction
!***********************************************************************
  use kind_types, only: qp, intkind1, intkind2, qp
  use ol_data_types_qp, only: wfun, hol
  use hel_bookkeeping_qp, only: helbookkeeping_ol_vert3
  use ol_vert_interface_qp, only: loop_AQ_S
  use ol_parameters_init_qp, only: get_coupling

  implicit none
  integer(intkind1), intent(in)    :: ntry
  integer(intkind2), intent(inout) :: n(3), t(:,:)
  integer,           intent(in)    :: ng_RL
  type(wfun),        intent(in)    :: J_Q(:)
  type(hol),         intent(inout) :: G_A,Gout_S
  complex(qp) :: G_add(size(Gout_S%j,1),size(Gout_S%j,2),size(Gout_S%j,3))
  complex(qp) :: g_RL(2)
  integer(intkind2) :: h


  if (ntry == 1) call helbookkeeping_ol_vert3(ntry, J_Q, G_A, Gout_S, n, t)
  if (.not. valid_hol(G_A, Gout_S)) return

  Gout_S%j = 0._qp
  g_RL = get_coupling(ng_RL)
  do h = 1, n(3)  ! recursion step
    call loop_AQ_S(G_A%j(:,:,:,h),J_Q(t(1,h))%j,G_add,g_RL)
    Gout_S%j(:,:,:,t(2,h)) = Gout_S%j(:,:,:,t(2,h)) + G_add
  end do



end subroutine Hloop_AQ_S


!***********************************************************************
subroutine Hloop_VV_S(ntry, G_V, J_V, Gout_S, n, t)
!-----------------------------------------------------------------------
! bare vector vector -> scalar interaction
! ----------------------------------------------------------------------
  use kind_types, only: qp, intkind1, intkind2, qp
  use ol_data_types_qp, only: wfun, hol
  use hel_bookkeeping_qp, only: helbookkeeping_ol_vert3
  use ol_vert_interface_qp, only: loop_VV_S

  implicit none
  integer(intkind1), intent(in)    :: ntry
  integer(intkind2), intent(inout) :: n(3), t(:,:)
  type(wfun),  intent(in)    :: J_V(:)
  type(hol),   intent(inout) :: G_V,Gout_S
  complex(qp)  :: G_add(size(Gout_S%j,1),size(Gout_S%j,2),size(Gout_S%j,3))
  integer(intkind2)  :: h


  if (ntry == 1) call helbookkeeping_ol_vert3(ntry, J_V, G_V, Gout_S, n, t)
  if (.not. valid_hol(G_V, Gout_S)) return

  Gout_S%j = 0._qp
  do h = 1, n(3)  ! recursion step
    call loop_VV_S(G_V%j(:,:,:,h), J_V(t(1,h))%j, G_add)
    Gout_S%j(:,:,:,t(2,h)) = Gout_S%j(:,:,:,t(2,h)) + G_add
  end do



end subroutine Hloop_VV_S


!***********************************************************************
subroutine Hloop_VS_V(ntry, G_V, J_S, Gout_V, n, t)
!-----------------------------------------------------------------------
! bare vector scalar -> vector interaction
! ----------------------------------------------------------------------
  use kind_types, only: qp, intkind1, intkind2, qp
  use ol_data_types_qp, only: wfun, hol
  use hel_bookkeeping_qp, only: helbookkeeping_ol_vert3
  use ol_vert_interface_qp, only: loop_VS_V

  implicit none
  integer(intkind1), intent(in)    :: ntry
  integer(intkind2), intent(inout) :: n(3), t(:,:)
  type(wfun),  intent(in)    :: J_S(:)
  type(hol),   intent(inout) :: G_V,Gout_V
  complex(qp)  :: G_add(size(Gout_V%j,1),size(Gout_V%j,2),size(Gout_V%j,3))
  integer(intkind2)  :: h


  if (ntry == 1) call helbookkeeping_ol_vert3(ntry, J_S, G_V, Gout_V, n, t)
  if (.not. valid_hol(G_V, Gout_V)) return

  Gout_V%j = 0._qp
  do h = 1, n(3)  ! recursion step
    call loop_VS_V(G_V%j(:,:,:,h), J_S(t(1,h))%j, G_add)
    Gout_V%j(:,:,:,t(2,h)) = Gout_V%j(:,:,:,t(2,h)) + G_add
  end do



end subroutine Hloop_VS_V


!***********************************************************************
subroutine Hloop_SV_V(ntry, G_S, J_V, Gout_V, n, t)
!-----------------------------------------------------------------------
! bare scalar vector -> vector interaction
! ----------------------------------------------------------------------
  use kind_types, only: qp, intkind1, intkind2, qp
  use ol_data_types_qp, only: wfun, hol
  use hel_bookkeeping_qp, only: helbookkeeping_ol_vert3
  use ol_vert_interface_qp, only: loop_SV_V

  implicit none
  integer(intkind1), intent(in)    :: ntry
  integer(intkind2), intent(inout) :: n(3), t(:,:)
  type(wfun),  intent(in)    :: J_V(:)
  type(hol),   intent(inout) :: G_S,Gout_V
  complex(qp)  :: G_add(size(Gout_V%j,1),size(Gout_V%j,2),size(Gout_V%j,3))
  integer(intkind2)  :: h


  if (ntry == 1) call helbookkeeping_ol_vert3(ntry, J_V, G_S, Gout_V, n, t)
  if (.not. valid_hol(G_S, Gout_V)) return

  Gout_V%j = 0._qp
  do h = 1, n(3)  ! recursion step
    call loop_SV_V(G_S%j(:,:,:,h), J_V(t(1,h))%j, G_add)
    Gout_V%j(:,:,:,t(2,h)) = Gout_V%j(:,:,:,t(2,h)) + G_add
  end do



end subroutine Hloop_SV_V

!******************************************************************************!
!                                 EW Vertices                                  !
!******************************************************************************!

!***********************************************************************
subroutine Hloop_SV_T(ntry, G_S, moml, J_V, momt, Gout_S, n, t)
! ----------------------------------------------------------------------
! bare scalar-vector -> scalar interaction
!***********************************************************************
  use kind_types, only: qp, qp, intkind1, intkind2
  use ol_data_types_qp, only: wfun, hol
  use hel_bookkeeping_qp, only: helbookkeeping_ol_vert3
  use ol_kinematics_qp, only: get_LC_4
  use ol_vert_interface_qp, only: loop_SV_T

  implicit none
  integer(intkind1), intent(in)    :: ntry
  integer(intkind2), intent(inout) :: n(3), t(:,:)
  integer,           intent(in)    :: moml,momt
  type(wfun),        intent(in)    :: J_V(:)
  type(hol),         intent(inout) :: G_S,Gout_S
  complex(qp)  :: G_add(size(Gout_S%j,1),size(Gout_S%j,2),size(Gout_S%j,3))
  integer(intkind2)  :: h


  if (ntry == 1) call helbookkeeping_ol_vert3(ntry, J_V, G_S, Gout_S, n, t)
  if (.not. valid_hol(G_S, Gout_S)) return

  Gout_S%j = 0._qp
  G_add    = 0._qp

  do h = 1, n(3)  ! recursion step
    call loop_SV_T(G_S%j(:,:,:,h), get_LC_4(moml), J_V(t(1,h))%j, &
                   get_LC_4(momt), G_add)
    Gout_S%j(:,:,:,t(2,h)) = Gout_S%j(:,:,:,t(2,h)) + G_add
  end do



end subroutine Hloop_SV_T

!***********************************************************************
subroutine Hloop_TV_S(ntry, G_S, moml, J_V, momt, Gout_S, n, t)
! ----------------------------------------------------------------------
! bare scalar-vector -> scalar interaction
!***********************************************************************
  use kind_types, only: qp, qp, intkind1, intkind2
  use ol_data_types_qp, only: wfun, hol
  use hel_bookkeeping_qp, only: helbookkeeping_ol_vert3
  use ol_vert_interface_qp, only: loop_TV_S
  use ol_kinematics_qp, only: get_LC_4

  implicit none
  integer(intkind1), intent(in)    :: ntry
  integer(intkind2), intent(inout) :: n(3), t(:,:)
  integer,           intent(in)    :: moml,momt
  type(wfun),        intent(in)    :: J_V(:)
  type(hol),         intent(inout) :: G_S,Gout_S
  complex(qp)  :: G_add(size(Gout_S%j,1),size(Gout_S%j,2),size(Gout_S%j,3))
  integer(intkind2)  :: h


  if (ntry == 1) call helbookkeeping_ol_vert3(ntry, J_V, G_S, Gout_S, n, t)
  if (.not. valid_hol(G_S, Gout_S)) return

  Gout_S%j = 0._qp
  G_add    = 0._qp

  do h = 1, n(3)  ! recursion step
    call loop_TV_S(G_S%j(:,:,:,h), get_LC_4(moml), J_V(t(1,h))%j, get_LC_4(momt), G_add)
    Gout_S%j(:,:,:,t(2,h)) = Gout_S%j(:,:,:,t(2,h)) + G_add
  end do



end subroutine Hloop_TV_S


!***********************************************************************
subroutine Hloop_VS_T(ntry, G_V, moml, J_S, momt, Gout_S, n, t)
! ----------------------------------------------------------------------
! bare vector-scalar -> scalar interaction
!***********************************************************************
  use kind_types, only: qp, qp, intkind1, intkind2
  use ol_data_types_qp, only: wfun, hol
  use hel_bookkeeping_qp, only: helbookkeeping_ol_vert3
  use ol_vert_interface_qp, only: loop_VS_T
  use ol_kinematics_qp, only: get_LC_4

  implicit none
  integer(intkind1), intent(in)    :: ntry
  integer(intkind2), intent(inout) :: n(3), t(:,:)
  integer,           intent(in)    :: moml,momt
  type(wfun),        intent(in)    :: J_S(:)
  type(hol),         intent(inout) :: G_V,Gout_S
  complex(qp)  :: G_add(size(Gout_S%j,1),size(Gout_S%j,2),size(Gout_S%j,3))
  integer(intkind2)  :: h


  if (ntry == 1) call helbookkeeping_ol_vert3(ntry, J_S, G_V, Gout_S, n, t)
  if (.not. valid_hol(G_V, Gout_S)) return

  Gout_S%j = 0._qp
  G_add    = 0._qp

  do h = 1, n(3)  ! recursion step
    call loop_VS_T(G_V%j(:,:,:,h), get_LC_4(moml), J_S(t(1,h))%j, get_LC_4(momt), G_add)
    Gout_S%j(:,:,:,t(2,h)) = Gout_S%j(:,:,:,t(2,h)) + G_add
  end do



end subroutine Hloop_VS_T


!***********************************************************************
subroutine Hloop_VT_S(ntry, G_V, moml, J_S, momt, Gout_S, n, t)
! ----------------------------------------------------------------------
! bare vector-scalar -> scalar interaction
!***********************************************************************
  use kind_types, only: qp, qp, intkind1, intkind2
  use ol_data_types_qp, only: wfun, hol
  use hel_bookkeeping_qp, only: helbookkeeping_ol_vert3
  use ol_vert_interface_qp, only: loop_VT_S
  use ol_kinematics_qp, only: get_LC_4

  implicit none
  integer(intkind1), intent(in)    :: ntry
  integer(intkind2), intent(inout) :: n(3), t(:,:)
  integer,           intent(in)    :: moml,momt
  type(wfun),        intent(in)    :: J_S(:)
  type(hol),         intent(inout) :: G_V,Gout_S
  complex(qp)  :: G_add(size(Gout_S%j,1),size(Gout_S%j,2),size(Gout_S%j,3))
  integer(intkind2)  :: h


  if (ntry == 1) call helbookkeeping_ol_vert3(ntry, J_S, G_V, Gout_S, n, t)
  if (.not. valid_hol(G_V, Gout_S)) return

  Gout_S%j = 0._qp
  G_add    = 0._qp

  do h = 1, n(3)  ! recursion step
    call loop_VT_S(G_V%j(:,:,:,h), get_LC_4(moml), J_S(t(1,h))%j, get_LC_4(momt), G_add)
    Gout_S%j(:,:,:,t(2,h)) = Gout_S%j(:,:,:,t(2,h)) + G_add
  end do



end subroutine Hloop_VT_S


!***********************************************************************
subroutine Hloop_ST_V(ntry, G_S, moml, J_S, momt, Gout_V, n, t)
! ----------------------------------------------------------------------
! bare scalar-scalar -> vector interaction
!***********************************************************************
  use kind_types, only: qp, qp, intkind1, intkind2
  use ol_data_types_qp, only: wfun, hol
  use hel_bookkeeping_qp, only: helbookkeeping_ol_vert3
  use ol_vert_interface_qp, only: loop_ST_V
  use ol_kinematics_qp, only: get_LC_4

  implicit none
  integer(intkind1), intent(in)    :: ntry
  integer(intkind2), intent(inout) :: n(3), t(:,:)
  integer,           intent(in)    :: moml,momt
  type(wfun),        intent(in)    :: J_S(:)
  type(hol),         intent(inout) :: G_S,Gout_V
  complex(qp)  :: G_add(size(Gout_V%j,1),size(Gout_V%j,2),size(Gout_V%j,3))
  integer(intkind2)  :: h


  if (ntry == 1) call helbookkeeping_ol_vert3(ntry, J_S, G_S, Gout_V, n, t)
  if (.not. valid_hol(G_S, Gout_V)) return

  Gout_V%j = 0._qp
  G_add    = 0._qp

  do h = 1, n(3)  ! recursion step
    call loop_ST_V(G_S%j(:,:,:,h), get_LC_4(moml), J_S(t(1,h))%j, get_LC_4(momt), G_add)
    Gout_V%j(:,:,:,t(2,h)) = Gout_V%j(:,:,:,t(2,h)) + G_add
  end do



end subroutine Hloop_ST_V


!***********************************************************************
subroutine Hloop_TS_V(ntry, G_S, moml, J_S, momt, Gout_V, n, t)
! ----------------------------------------------------------------------
! bare scalar-scalar -> vector interaction
!***********************************************************************
  use kind_types, only: qp, qp, intkind1, intkind2
  use ol_data_types_qp, only: wfun, hol
  use hel_bookkeeping_qp, only: helbookkeeping_ol_vert3
  use ol_vert_interface_qp, only: loop_TS_V
  use ol_kinematics_qp, only: get_LC_4

  implicit none
  integer(intkind1), intent(in)    :: ntry
  integer(intkind2), intent(inout) :: n(3), t(:,:)
  integer,           intent(in)    :: moml,momt
  type(wfun),        intent(in)    :: J_S(:)
  type(hol),         intent(inout) :: G_S,Gout_V
  complex(qp)  :: G_add(size(Gout_V%j,1),size(Gout_V%j,2),size(Gout_V%j,3))
  integer(intkind2)  :: h


  if (ntry == 1) call helbookkeeping_ol_vert3(ntry, J_S, G_S, Gout_V, n, t)
  if (.not. valid_hol(G_S, Gout_V)) return

  Gout_V%j = 0._qp
  G_add    = 0._qp

  do h = 1, n(3)  ! recursion step
    call loop_TS_V(G_S%j(:,:,:,h), get_LC_4(moml), J_S(t(1,h))%j, get_LC_4(momt), G_add)
    Gout_V%j(:,:,:,t(2,h)) = Gout_V%j(:,:,:,t(2,h)) + G_add
  end do



end subroutine Hloop_TS_V


!***********************************************************************
subroutine Hloop_SS_S(ntry, G_S, J_S, Gout_S, n, t)
! ----------------------------------------------------------------------
! bare scalar scalar -> scalar interaction
! **********************************************************************
  use kind_types, only: qp, qp, intkind1, intkind2
  use ol_data_types_qp, only: wfun, hol
  use hel_bookkeeping_qp, only: helbookkeeping_ol_vert3
  use ol_vert_interface_qp, only: loop_SS_S

  implicit none
  integer(intkind1), intent(in)    :: ntry
  integer(intkind2), intent(inout) :: n(3), t(:,:)
  type(wfun),        intent(in)    :: J_S(:)
  type(hol),         intent(inout) :: G_S,Gout_S
  complex(qp)  :: G_add(size(Gout_S%j,1),size(Gout_S%j,2),size(Gout_S%j,3))
  integer(intkind2)  :: h


  if (ntry == 1) call helbookkeeping_ol_vert3(ntry, J_S, G_S, Gout_S, n, t)
  if (.not. valid_hol(G_S, Gout_S)) return

  Gout_S%j = 0._qp
  G_add    = 0._qp

  do h = 1, n(3)  ! recursion step
    call loop_SS_S(G_S%j(:,:,:,h), J_S(t(1,h))%j, G_add)
    Gout_S%j(:,:,:,t(2,h)) = Gout_S%j(:,:,:,t(2,h)) + G_add
  end do



end subroutine Hloop_SS_S


!***********************************************************************
subroutine Hloop_SSS_S(ntry, Gin_S, J_S1, J_S2, Gout_S, n, t)
!-----------------------------------------------------------------------
! bare scalar vector vector -> scalar interaction
! **********************************************************************
  use kind_types, only: qp, qp, intkind1, intkind2
  use ol_data_types_qp, only: wfun, hol
  use hel_bookkeeping_qp, only: helbookkeeping_ol_vert4
  use ol_vert_interface_qp, only: loop_SSS_S

  implicit none
  integer(intkind1), intent(in)    :: ntry
  integer(intkind2), intent(inout) :: t(:,:)
  integer(intkind2), intent(inout) :: n(4)
  type(wfun),        intent(in)    :: J_S1(:), J_S2(:)
  type(hol),         intent(inout) :: Gin_S,Gout_S
  complex(qp)  :: G_add(size(Gout_S%j,1),size(Gout_S%j,2),size(Gout_S%j,3))
  integer  :: h


  if (ntry == 1) call helbookkeeping_ol_vert4(ntry, J_S1, J_S2, Gin_S, Gout_S, n, t)  ! STILL TO BE ADAPTED
  if (.not. valid_hol(Gin_S, Gout_S)) return

  Gout_S%j = 0._qp
  G_add    = 0._qp

  do h = 1, n(4)  ! recursion step
    call loop_SSS_S(Gin_S%j(:,:,:,h), J_S1(t(1,h))%j, J_S2(t(2,h))%j, G_add)
    Gout_S%j(:,:,:,t(3,h)) = Gout_S%j(:,:,:,t(3,h)) + G_add
  end do





end subroutine Hloop_SSS_S


!***********************************************************************
subroutine Hloop_SVV_S(ntry, Gin_S, J_V1, J_V2, Gout_S, n, t)
!-----------------------------------------------------------------------
! bare scalar vector vector -> scalar interaction
! **********************************************************************
  use kind_types, only: qp, qp, intkind1, intkind2
  use ol_data_types_qp, only: wfun, hol
  use hel_bookkeeping_qp, only: helbookkeeping_ol_vert4
  use ol_vert_interface_qp, only: loop_SVV_S

  implicit none
  integer(intkind1), intent(in)    :: ntry
  integer(intkind2), intent(inout) :: t(:,:)
  integer(intkind2), intent(inout) :: n(4)
  type(wfun),        intent(in)    :: J_V1(:), J_V2(:)
  type(hol),         intent(inout) :: Gin_S,Gout_S
  complex(qp)  :: G_add(size(Gout_S%j,1),size(Gout_S%j,2),size(Gout_S%j,3))
  integer  :: h


  if (ntry == 1) call helbookkeeping_ol_vert4(ntry, J_V1, J_V2, Gin_S, Gout_S, n, t)  ! STILL TO BE ADAPTED
  if (.not. valid_hol(Gin_S, Gout_S)) return

  Gout_S%j = 0._qp
  G_add    = 0._qp

  do h = 1, n(4)  ! recursion step
    call loop_SVV_S(Gin_S%j(:,:,:,h), J_V1(t(1,h))%j, J_V2(t(2,h))%j, G_add)
    Gout_S%j(:,:,:,t(3,h)) = Gout_S%j(:,:,:,t(3,h)) + G_add
  end do



end subroutine Hloop_SVV_S


!***********************************************************************
subroutine Hloop_VSS_V(ntry, Gin_V, J_S1, J_S2, Gout_V, n, t)
!-----------------------------------------------------------------------
! bare vector scalar scalar -> vector interaction
! **********************************************************************
  use kind_types, only: qp, qp, intkind1, intkind2
  use ol_data_types_qp, only: wfun, hol
  use hel_bookkeeping_qp, only: helbookkeeping_ol_vert4
  use ol_vert_interface_qp, only: loop_VSS_V

  implicit none
  integer(intkind1), intent(in)    :: ntry
  integer(intkind2), intent(inout) :: t(:,:)
  integer(intkind2), intent(inout) :: n(4)
  type(wfun),        intent(in)    :: J_S1(:), J_S2(:)
  type(hol),         intent(inout) :: Gin_V,Gout_V
  complex(qp)  :: G_add(size(Gout_V%j,1),size(Gout_V%j,2),size(Gout_V%j,3))
  integer  :: h


  if (ntry == 1) call helbookkeeping_ol_vert4(ntry, J_S1, J_S2, Gin_V, Gout_V, n, t)  ! STILL TO BE ADAPTED
  if (.not. valid_hol(Gin_V, Gout_V)) return

  Gout_V%j = 0._qp
  G_add    = 0._qp

  do h = 1, n(4)  ! recursion step
    call loop_VSS_V(Gin_V%j(:,:,:,h), J_S1(t(1,h))%j, J_S2(t(2,h))%j, G_add)
    Gout_V%j(:,:,:,t(3,h)) = Gout_V%j(:,:,:,t(3,h)) + G_add
  end do



end subroutine Hloop_VSS_V


!***********************************************************************
subroutine Hloop_VVS_S(ntry, Gin_V, J_V, J_S, Gout_S, n, t)
!-----------------------------------------------------------------------
! bare vector vector scalar -> scalar interaction
! **********************************************************************
  use kind_types, only: qp, qp, intkind1, intkind2
  use ol_data_types_qp, only: wfun, hol
  use hel_bookkeeping_qp, only: helbookkeeping_ol_vert4
  use ol_vert_interface_qp, only: loop_VVS_S

  implicit none
  integer(intkind1), intent(in)    :: ntry
  integer(intkind2), intent(inout) :: t(:,:)
  integer(intkind2), intent(inout) :: n(4)
  type(wfun),        intent(in)    :: J_V(:), J_S(:)
  type(hol),         intent(inout) :: Gin_V,Gout_S
  complex(qp)  :: G_add(size(Gout_S%j,1),size(Gout_S%j,2),size(Gout_S%j,3))
  integer  :: h


  if (ntry == 1) call helbookkeeping_ol_vert4(ntry, J_V, J_S, Gin_V, Gout_S, n, t)  ! STILL TO BE ADAPTED
  if (.not. valid_hol(Gin_V, Gout_S)) return

  Gout_S%j = 0._qp
  G_add    = 0._qp

  do h = 1, n(4)  ! recursion step
    call loop_VVS_S(Gin_V%j(:,:,:,h), J_V(t(1,h))%j, J_S(t(2,h))%j, G_add)
    Gout_S%j(:,:,:,t(3,h)) = Gout_S%j(:,:,:,t(3,h)) + G_add
  end do



end subroutine Hloop_VVS_S


!***********************************************************************
subroutine Hloop_SSV_V(ntry, Gin_S, J_S, J_V, Gout_V, n, t)
!-----------------------------------------------------------------------
! bare scalar scalar vector -> vector interaction
! **********************************************************************
  use kind_types, only: qp, qp, intkind1, intkind2
  use ol_data_types_qp, only: wfun, hol
  use hel_bookkeeping_qp, only: helbookkeeping_ol_vert4
  use ol_vert_interface_qp, only: loop_SSV_V

  implicit none
  integer(intkind1), intent(in)    :: ntry
  integer(intkind2), intent(inout) :: t(:,:)
  integer(intkind2), intent(inout) :: n(4)
  type(wfun),        intent(in)    :: J_S(:), J_V(:)
  type(hol),         intent(inout) :: Gin_S,Gout_V
  complex(qp)  :: G_add(size(Gout_V%j,1),size(Gout_V%j,2),size(Gout_V%j,3))
  integer  :: h


  if (ntry == 1) call helbookkeeping_ol_vert4(ntry, J_S, J_V, Gin_S, Gout_V, n, t)  ! STILL TO BE ADAPTED
  if (.not. valid_hol(Gin_S, Gout_V)) return

  Gout_V%j = 0._qp
  G_add    = 0._qp

  do h = 1, n(4)  ! recursion step
    call loop_SSV_V(Gin_S%j(:,:,:,h), J_S(t(1,h))%j, J_V(t(2,h))%j, G_add)
    Gout_V%j(:,:,:,t(3,h)) = Gout_V%j(:,:,:,t(3,h)) + G_add
  end do



end subroutine Hloop_SSV_V


!***********************************************************************
subroutine Hloop_VWW_V(ntry, Gin_V, J_V1, J_V2, Gout_V, n, t)
!-----------------------------------------------------------------------
! bare vector vector vector -> vector interaction
! **********************************************************************
  use kind_types, only: qp, qp, intkind1, intkind2
  use ol_data_types_qp, only: wfun, hol
  use hel_bookkeeping_qp, only: helbookkeeping_ol_vert4
  use ol_vert_interface_qp, only: loop_VWW_V

  implicit none
  integer(intkind1), intent(in)    :: ntry
  integer(intkind2), intent(inout) :: t(:,:)
  integer(intkind2), intent(inout) :: n(4)
  type(wfun),        intent(in)    :: J_V1(:), J_V2(:)
  type(hol),         intent(inout) :: Gin_V,Gout_V
  complex(qp)  :: G_add(size(Gout_V%j,1),size(Gout_V%j,2),size(Gout_V%j,3))
  integer  :: h


  if (ntry == 1) call helbookkeeping_ol_vert4(ntry, J_V1, J_V2, Gin_V, Gout_V, n, t)  ! STILL TO BE ADAPTED
  if (.not. valid_hol(Gin_V, Gout_V)) return

  Gout_V%j = 0._qp
  G_add    = 0._qp

  do h = 1, n(4)  ! recursion step
    call loop_VWW_V(Gin_V%j(:,:,:,h), J_V1(t(1,h))%j, J_V2(t(2,h))%j, G_add)
    Gout_V%j(:,:,:,t(3,h)) = Gout_V%j(:,:,:,t(3,h)) + G_add
  end do



end subroutine Hloop_VWW_V

!***********************************************************************
subroutine Hloop_WWV_V(ntry, Gin_V, J_V1, J_V2, Gout_V, n, t)
!-----------------------------------------------------------------------
! bare vector vector vector -> vector interaction
! **********************************************************************
  use kind_types, only: qp, qp, intkind1, intkind2
  use ol_data_types_qp, only: wfun, hol
  use hel_bookkeeping_qp, only: helbookkeeping_ol_vert4
  use ol_vert_interface_qp, only: loop_WWV_V

  implicit none
  integer(intkind1), intent(in)    :: ntry
  integer(intkind2), intent(inout) :: t(:,:)
  integer(intkind2), intent(inout) :: n(4)
  type(wfun),        intent(in)    :: J_V1(:), J_V2(:)
  type(hol),         intent(inout) :: Gin_V,Gout_V
  complex(qp)  :: G_add(size(Gout_V%j,1),size(Gout_V%j,2),size(Gout_V%j,3))
  integer  :: h


  if (ntry == 1) call helbookkeeping_ol_vert4(ntry, J_V1, J_V2, Gin_V, Gout_V, n, t)  ! STILL TO BE ADAPTED
  if (.not. valid_hol(Gin_V, Gout_V)) return

  Gout_V%j = 0._qp
  G_add    = 0._qp

  do h = 1, n(4)  ! recursion step
    call loop_WWV_V(Gin_V%j(:,:,:,h), J_V1(t(1,h))%j, J_V2(t(2,h))%j, G_add)
    Gout_V%j(:,:,:,t(3,h)) = Gout_V%j(:,:,:,t(3,h)) + G_add
  end do



end subroutine Hloop_WWV_V

end module ol_h_vert_interface_qp

