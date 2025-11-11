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


module ol_i_operator_qp
  use ol_debug, only: ol_fatal
  implicit none
  contains

! **********************************************************************
subroutine intdip(mode, M2LO, M2CC, M2CC_EW, extflav, extcharges, Npart, &
                  extmass2, sarr, vdip, c_dip, olm_in, photonid_in)
! **********************************************************************
! I-Operator contribution to integrated Dipoles <=> Eqs. (6.66), (6.52), (6.16)
! in hep-ph/0201036 (Catani, Dittmaier, Seymour, Trocsanyi)
!
! Only the 5-flavour scheme is implemented:
! N_f = 5 active massless flavours in g -> qq~ splittings
! N_F = 0 active massive flavours in  g -> qq~ splittings, {m_F} = {}
! **********************************************************************
! mode =  I-operator switcher: 1 QCD; 2 EM; 3 QCD+EM; otherwise OFF
!
! M2LO  = Born amplitude
!
! M2CC(i,j) = <M_0|T^a(i)T^a(j)|M_0>,    i,j=1,...,Npart
!           = colour-correlated Born amplitude summed/averaged
!             over colour and polarization
!
! extflav(k) = flavour of external particle k
!            = 1 for gluon | 2 for quark/anti-quark | 3 otherwise
!
! extmass2(k) = squared mass of external particle k
!
! extcharges(k) = incoming EM charge of kth particle
!               = -1 for incoming-electrons/outgoing-positrons
!
! sarr(i,j) = (p_i+p_j)^2 invariants (all external p_i incoming)
! ----------------------------------------------------------------------
!
! vdip = K(eps)*(C_dip(0) + C_dip(1)/eps    + C_dip(2)/eps**2)
!      =         C_dip(0) + C_dip(1)*de1_IR + C_dip(2)*de2_i_IR
!      = I-Operator cont. to |M|^2 in D=4-2*eps dimensions
!        summed/averaged over colour and polarization
!
! finite part C_dip(0) of Laurent series depends on choice of
! normalisation factor K(eps). This can be adpated through switcher
! norm_swi in subroutine loop_parameters_decl. The Binoth-LH accord convention
!   K(eps) = (4Pi)^eps/Gamma(1-eps)
! is the default normalisation (norm_swi=0)
! **********************************************************************
  use kind_types, only: qp
  use ol_parameters_decl_qp, only: alpha_QCD, alpha_QED, pi
  use ol_parameters_decl_dp, only: &
       & offshell_photons_lsz, delta_alphamz_dimreg
  use ol_loop_parameters_decl_qp
  use ol_momenta_decl_qp, only: Q, L
  use ol_kinematics_qp
  implicit none

  integer,           intent(in)  :: mode
  integer,           intent(in)  :: Npart
  integer,           intent(in)  :: extflav(Npart)
  real(qp),    intent(in)  :: extcharges(Npart)
  real(qp),    intent(in)  :: M2LO, M2CC(Npart,Npart), M2CC_EW, extmass2(Npart)
  complex(qp), intent(in)  :: sarr(:,:)
  integer, optional, intent(in)  :: olm_in, photonid_in(:)
  real(qp),    intent(out) :: vdip, c_dip(0:2)
  real(qp) :: Q2_aux, Fjk(0:2), Gj(0:2), norm_qcd, norm_qed, QQ
  integer        :: olm, photonid(Npart), i, j, k, mid

  Q2_aux = mureg2  ! arbitrary auxiliary scale
  c_dip = 0

  if (present(olm_in)) then
    olm = olm_in
  else
    olm = 0
  end if

  if (present(photonid_in)) then
    photonid = photonid_in
  else
    photonid = 0
  end if

  !NLO QCD contribution
  if (mode == 1 .or. mode == 0) then
    norm_qcd = alpha_QCD/(4*pi) ! normalization in intermediate formulae

    do j = 1, Npart
      if (extflav(j) >= 3 .or. extflav(j) < 1) cycle ! emitter j = QCD singlet
      call intdip_Gj(j,extflav(j),extmass2(j),Q2_aux,Gj)
      do k = 1, Npart
        if (extflav(k) >= 3 .or. extflav(k) < 1 .or. k == j) cycle ! spectator k = QCD singlet
        if (olm > 0) then
          mid = 2**(j-1)+2**(k-1)
          call intdip_Fjk(j,k,real(L(5,mid)+L(6,mid)),extflav(j),extmass2(j),extmass2(k),Q2_aux,Fjk)
        else
          call intdip_Fjk(j,k,real(sarr(j,k)),extflav(j),extmass2(j),extmass2(k),Q2_aux,Fjk)
        end if
        c_dip = c_dip - 2*norm_qcd*M2CC(j,k) * (Fjk + Gj)
      end do
    end do

  end if


  !NLO QED contribution
  if (mode == 2 .or. mode == 0) then
    norm_qed = alpha_QED/(4*pi) ! normalization in intermediate formulae

    do j = 1, Npart
      if (extflav(j) == 1 .or. extcharges(j) == 0. ) cycle ! emitter j = EM neutral
      call intdip_Gj(j,extflav(j),extmass2(j),Q2_aux,Gj)
      do k = 1, Npart
        if (extflav(k) == 1 .or. extcharges(k) == 0. .or. k == j) cycle ! spectator k = EW singlet
        if (olm > 0) then
          mid = 2**(j-1)+2**(k-1)
          call intdip_Fjk(j,k,real(L(5,mid)+L(6,mid)),extflav(j),extmass2(j),extmass2(k),Q2_aux,Fjk)
        else
          call intdip_Fjk(j,k,real(sarr(j,k)),extflav(j),extmass2(j),extmass2(k),Q2_aux,Fjk)
        end if
        c_dip = c_dip - 2*norm_qed*M2CC_EW*extcharges(j)*extcharges(k)*(Fjk + Gj)
      end do
    end do

    ! photon -> ff splittings
    if (offshell_photons_lsz .or. delta_alphamz_dimreg) then
      do j = 1, Npart
        if (extflav(j) /= -1) cycle ! emitter j /= photon
        Gj=0
        call intdip_Gj(j,extflav(j),extmass2(j),Q2_aux,Gj)
        do k = 1, Npart
          QQ = 0d0
          if (k == j) cycle
          if (.not. delta_alphamz_dimreg) then
            if (photonid(j) >= 0) cycle ! on-shell photon emitter -> cycle
          end if
          if (olm == 0) then
            if (real(Q(1,2**(j-1))+Q(2,2**(j-1))) .gt. 0 .and. real(Q(1,2**(k-1))+Q(2,2**(k-1))) .gt. 0 ) then ! IS off-shell photon -> IS spectator with Q=-1
              QQ = -1d0
            else if (real(Q(1,2**(j-1))+Q(2,2**(j-1))) .lt. 0 .and. real(Q(1,2**(k-1))+Q(2,2**(k-1))) .gt. 0 ) then ! FS off-shell photon -> both IS spectators with Q=-1/2
              QQ = -0.5d0
            else
              cycle
            end if
          else
            if (real(L(1,2**(j-1))+L(2,2**(j-1))) .gt. 0 .and. real(L(1,2**(k-1))+L(2,2**(k-1))) .gt. 0 ) then ! IS off-shell photon -> IS spectator with Q=-1
              QQ = -1d0
            else if (real(L(1,2**(j-1))+L(2,2**(j-1))) .lt. 0 .and. real(L(1,2**(k-1))+L(2,2**(k-1))) .gt. 0 ) then ! FS off-shell photon -> both IS spectators with Q=-1/2
              QQ = -0.5d0
            else
              cycle
            end if
          end if
          if (olm > 0) then
            mid = 2**(j-1)+2**(k-1)
            call intdip_Fjk(j,k,real(L(5,mid)+L(6,mid)),extflav(j),extmass2(j),extmass2(k),Q2_aux,Fjk)
          else
            call intdip_Fjk(j,k,real(sarr(j,k)),extflav(j),extmass2(j),extmass2(k),Q2_aux,Fjk)
          end if
          c_dip = c_dip - 2*norm_qed*M2CC_EW*QQ*(Fjk + Gj)
        end do
      end do
    end if

  end if


  ! adapf finite term of Laurent series: LH-normalisation -> generic normalisation
  c_dip(0) = c_dip(0) - c_dip(2)*de2_i_shift
  c_dip(1) = c_dip(1)
  c_dip(2) = c_dip(2)

  ! sum truncated Laurent series
  vdip = c_dip(0) + c_dip(1)*de1_IR + c_dip(2)*de2_i_IR


end subroutine intdip


! **********************************************************************
subroutine intdip_Gj(j,flavj,M2j,Q2_aux,Gj)
! **********************************************************************
! spectator-independent contributions
! **********************************************************************
! flavj = flavour of emitter (j): 1=gluon, 2=quark/antiquark, 3=lepton/W, -1=gamma
! M2j = squared mass of emitter (j)
! Q2_aux =  arbitrary auxiliary scale
!
! The laurent series coefficients Gi(0:2) are computed using
! the Les-Houches accord normalisation
! **********************************************************************
  use kind_types, only: qp
  use ol_debug, only: ol_fatal, ol_error
  use ol_generic, only: to_string
  use ol_parameters_decl_qp
  use ol_loop_parameters_decl_qp, only: pi2_6, tf, ca, mu2_IR => mureg2, Qf2sum
  use ol_loop_parameters_decl_dp, only: N_lf, SwF, SwB
  implicit none

  integer,        intent(in)  :: j, flavj
  real(qp), intent(in)  :: M2j, Q2_aux
  real(qp), intent(out) :: Gj(0:2)
  real(qp) :: Gaj, Kaj, Gammaj

  if (flavj == 1) then ! emitter j = gluon
    ! Gaj = [1/ca*] (11/6*ca - 2/3*tf*nl)
    ! Kaj = [1/ca*] ((67/18 - pi^2/6)*ca - 10/9*tf*nl)
    ! Gamma_j = Gaj/ep
    Gaj = 0
    Kaj = 0
    if (SwB /= 0) then
      Gaj = 11._qp/6
      Kaj = 67._qp/18 - pi2_6
    end if
    if (SwF /= 0) then
      Gaj = Gaj - 2*tf*N_lf/(3*ca)
      Kaj = Kaj - 10*tf*N_lf/(9*ca)
    end if
    Gj(1) = Gaj
    Gj(0) = 0
  else if (flavj == -1) then ! emitter j = photon
    ! Gaj = [1/1] (- 2/3*sum(Ncf*Qf2))
    ! Kaj = [1/1] (- 10/9*sum(Ncf*Qf2))
    ! Gamma_j = Gaj/ep
    Gaj = 0
    Kaj = 0
    if (SwF /= 0) then
      Gaj = Gaj - 2._qp/3.*Qf2sum
      Kaj = Kaj - 10._qp/9.*Qf2sum
    end if
    Gj(1) = Gaj
    Gj(0) = 0
  else if (flavj >= 2 .and. flavj <= 3) then ! EMITTER j = QUARK/LEPTON/W-BOSON
    ! QCD
    ! Gaj = [1/cf*] 3/2*cf
    ! Kaj = [1/cf*] (7/2 - pi^2/6)*cf
    ! QED (cf->Qf2)
    ! Gaj = [1/Qf2*] 3/2*Qf2
    ! Kaj = [1/Qf2*] (7/2 - pi^2/6)*Qf2
    Gaj = 1.5_qp
    Kaj = 3.5_qp - pi2_6
    if (M2j == 0) then
      ! Gamma_j = Gaj/ep
      Gj(1) = Gaj
      Gj(0) = 0
    else if (M2j > 0) then
      ! Gamma_j = [1/cf*] cf*(1/ep + 1/2*log(mj2/mu2) - 2)
      Gj(1) = 1
      Gj(0) = 0.5_qp*log(M2j/mu2_IR) - 2
    else
      call ol_error(2,'subroutine intdip_Gj: argument M2j out of range')
      call ol_fatal()
    end if
  else
    call ol_error(2,'subroutine intdip_Gj: argument flavj out of range')
!    call ol_error(2,'flavj = ' // to_string(flavj))
    call ol_fatal()
  end if

  ! Gj = Tj^2*([V_jk] - pi^2/3) + Gamma_j + Gaj + Kaj
  ! V_jk --> Fjk; Tj2 = cf for quark, ca for gluon
  if (flavj == -1) then ! emitter j = photon
    Gj(0) = Gj(0) + Gaj + Kaj
  else
    Gj(0) = Gj(0) + Gaj + Kaj - 2*pi2_6
  end if

  Gj(2) = 0

end subroutine intdip_Gj


! **********************************************************************
subroutine intdip_Fjk(j,k,Tjk,flavj,M2j,M2k,Q2_aux,Fjk)
! **********************************************************************
! spectator-dependent contributions
! Tjk=(p_j+p_k)^2;
! flavj = flavour of emitter (j): 1=gluon, 2=quark/antiquark, 3=lepton/W, -1=gamma
! M2j = squared mass of emitter (j)
! M2k = squared mass of spectator (k)
! Q2_aux =  arbitrary auxiliary scale
!
! The laurent series coefficients Fjk(0:2) are computed
! using the Les-Houches accord normalisation
!
! Fjk= function nu_j(s_jk,...) in (6.18) plus log-term in last line of (6.16)
! if j=initial-state gluon => kappa=2/3 according to (6.52)
! **********************************************************************
  use kind_types, only: qp
  use ol_debug, only: ol_fatal, ol_error
  use ol_generic, only: to_string
  use ol_parameters_decl_qp
  use ol_loop_parameters_decl_qp, only: mu2_IR => mureg2, pi2_6, tf, ca, kappa, Qf2sum
  use ol_loop_parameters_decl_dp, only: N_lf, SwF, SwB
  implicit none

  integer,        intent(in)  :: j, k, flavj
  real(qp), intent(in)  :: Tjk, M2j, M2k, Q2_aux
  real(qp), intent(out) :: Fjk(0:2)

  real(qp) :: Nuj_sing,Nuj_nonsing
  real(qp) :: Gaj
  real(qp) :: Mj, Mk, Sjk, Q2jk, Qjk, r2jk, delta, Nujk
  real(qp) :: rho2, rho, shatjk, kappa_kin, rho2j, rho2k
  real(qp) :: logmu2_Sjk, pole1_IR(0:1), pole2_IR(0:2)
  real(qp) :: rhog1, rhog2

  ! kinematic quantities
  Mj   = sqrt(M2j)
  Mk   = sqrt(M2k)
  Sjk  = abs(Tjk-M2j-M2k)
  Q2jk = Sjk + M2j + M2k
  Qjk  = sqrt(Q2jk)
  r2jk = M2j*M2k/Sjk**2

  if (M2j*M2k == 0) then
    delta = 2*r2jk
    Nujk  = 1
  else
    if (abs(r2jk) < 0.001) then
      delta = 2*r2jk + 2*r2jk**2 + 4*r2jk**3 + 10*r2jk**4 + 28*r2jk**5
      Nujk  = 1 - delta
    else
      Nujk  = sqrt(1-4*r2jk)
      delta = 1 - Nujk
    end if
  end if

  rho2 = delta/(1+Nujk)
  rho  = sqrt(rho2)

  shatjk = Sjk*0.5_qp*(1+Nujk)
  kappa_kin = (shatjk+M2k)/(shatjk+M2j)

  rho2j = M2j/shatjk*kappa_kin
  rho2k = M2k/shatjk/kappa_kin

  ! Laurent expansion of poles [(mu2_IR/Sjk)^eps]/eps^{1,2}
  logmu2_Sjk  = log(mu2_IR/Sjk)
  pole1_IR(1) = 1
  pole1_IR(0) = logmu2_Sjk

  pole2_IR(2) = 1
  pole2_IR(1) = logmu2_Sjk
  pole2_IR(0) = 0.5_qp*logmu2_Sjk**2

  ! singular part of nu_j function
  if (M2j == 0 .and. M2k == 0) then ! eq. (6.20c)
    Fjk = pole2_IR
  else if (M2j > 0 .and. M2k == 0) then ! eq. (6.20b)
    Fjk(2) = 0.5_qp*pole2_IR(2)
    Fjk(1) = 0.5_qp*pole2_IR(1) + 0.5_qp*pole1_IR(1)*log(M2j/Sjk)
    Fjk(0) = 0.5_qp*pole2_IR(0) + 0.5_qp*pole1_IR(0)*log(M2j/Sjk) &
             - 0.25_qp*(log(M2j/Sjk))**2 - 0.5_qp*pi2_6           &
             - 0.5_qp*log(M2j/Sjk)*log(Sjk/Q2jk)                  &
             - 0.5_qp*log(M2j/Q2jk)*log(Sjk/Q2jk)
  else if (M2j == 0 .and. M2k > 0) then ! eq. (6.20b)
    Fjk(2) = 0.5_qp*pole2_IR(2)
    Fjk(1) = 0.5_qp*pole2_IR(1) + 0.5_qp*pole1_IR(1)*log(M2k/Sjk)
    Fjk(0) = 0.5_qp*pole2_IR(0) + 0.5_qp*pole1_IR(0)*log(M2k/Sjk) &
             - 0.25_qp*(log(M2k/Sjk))**2 - 0.5_qp*pi2_6           &
             - 0.5_qp*log(M2k/Sjk)*log(Sjk/Q2jk)                  &
             - 0.5_qp*log(M2k/Q2jk)*log(Sjk/Q2jk)
  else if (M2j > 0 .and. M2k > 0) then ! eq. (6.20a)
    Fjk(2) = 0
    Fjk(1) = pole1_IR(1)*log(rho)/Nujk
    Fjk(0) = (pole1_IR(0)*log(rho)              &
             -0.25_qp*(log(rho2j))**2 &
             -0.25_qp*(log(rho2k))**2 &
             -pi2_6                   &
             +log(rho)*log(Q2jk/Sjk))/Nujk
  else
    call ol_error(2,'subroutine intdip_Fjk: arguments out of range')
    call ol_error(2,'allowed range M2j, M2k >= 0')
!    write(*,*) '[OpenLoops] M2j = ', M2j
!    write(*,*) '[OpenLoops] M2k = ', M2k
    call ol_fatal()
  end if

  ! non-singular part of nu_j function for gluons (ok)
  if (flavj == 1) then ! EMITTER j=GLUON
    ! Gaj = [1/ca*] (11/6*ca - 2/3*tf*nl)
    Gaj = 0
    if (SwB /= 0) Gaj = Gaj + 11._qp/6
    if (SwF /= 0) Gaj = Gaj - 2*tf*N_lf/(3*ca)
    if (M2k == 0) then ! eq. (6.26)
      Nuj_nonsing = 0._qp ! ~ tf/ca
    else if (M2k > 0) then ! eq. (6.24) in N_f=5 scheme
      Nuj_nonsing = Gaj*(log(Sjk/Q2jk) - 2._qp*log((Qjk-Mk)/Qjk)-2._qp*Mk/(Qjk+Mk)) + pi2_6 - spence(Sjk/Q2jk)
      if (j > 2) then ! kappa-dep cont only for FS gluons, see nu_a in (6.52)
        ! Nuj_nonsing = Nuj_nonsing + (kappa-2/3)*M2k/Sjk*(2*tf*nl/ca-1)*log(2*Mk/(Qjk+Mk))
        if (SwB /= 0) Nuj_nonsing = Nuj_nonsing - (kappa-2._qp/3._qp)*M2k/Sjk*log(2*Mk/(Qjk+Mk))
        if (SwF /= 0) Nuj_nonsing = Nuj_nonsing + (kappa-2._qp/3._qp)*M2k/Sjk*log(2*Mk/(Qjk+Mk))*(2*tf*N_lf)/3
      end if
    else
      call ol_error(2,'subroutine intdip_Fjk: arguments out of range')
      call ol_error(2,'allowed range M2k >= 0')
 !     write(*,*) '[OpenLoops] M2k =', M2k
      call ol_fatal()
    end if

  ! non-singular part of nu_j function for photons
  else if (flavj == -1) then ! EMITTER j=PHOTON
    ! singular-part vanishing for photons
    Fjk = 0._qp
    ! Gaj = [1/1] - 2/3*sum(Ncf*Qf2)
    Gaj = 0._qp
    if (SwF /= 0) Gaj = Gaj - 2._qp/3.*Qf2sum
    if (M2k == 0) then ! eq. (6.26)
      Nuj_nonsing = 0._qp ! ~ tf/ca
    else if (M2k > 0) then ! eq. (6.24) in N_f=5 scheme
      Nuj_nonsing = Gaj*(log(Sjk/Q2jk) - 2._qp*log((Qjk-Mk)/Qjk)-2._qp*Mk/(Qjk+Mk)) + pi2_6 - spence(Sjk/Q2jk)
      if (j > 2) then ! kappa-dep cont only for FS photons, see nu_a in (6.52)
        ! Nuj_nonsing = Nuj_nonsing + (kappa-2/3)*M2k/Sjk*(2*tf*nl/ca-1)*log(2*Mk/(Qjk+Mk))
!        if (SwF /= 0) Nuj_nonsing = Nuj_nonsing + (kappa-2._qp/3._qp)*M2k/Sjk*log(2*Mk/(Qjk+Mk))*(2*tf*N_lf)/3
      end if
    else
      call ol_error(2,'subroutine intdip_Fjk: arguments out of range')
      call ol_error(2,'allowed range M2k >= 0')
 !     write(*,*) '[OpenLoops] M2k =', M2k
      call ol_fatal()
    end if

  ! non-singular part of nu_j function for quarks (ok)
  else if (flavj >= 2 .and. flavj <= 3) then! EMITTER j=QUARK/LEPTON/W-BOSON
    ! QCD
    ! Gaj = [1/cf*] 3/2*cf
    ! QED (cf->Qf2)
    ! Gaj = [1/Qf2*] 3/2*Qf2
    Gaj= 1.5_qp
    if (M2j == 0 .and. M2k == 0) then
      Nuj_nonsing = 0
    else if (M2j > 0 .and. M2k == 0) then ! eq. (6.22)
      Nuj_nonsing = Gaj*log(Sjk/Q2jk) + pi2_6 - spence(Sjk/Q2jk) - 2*log(Sjk/Q2jk) - M2j/Sjk*log(M2j/Q2jk)
    else if (M2j == 0 .and. M2k > 0) then ! eq. (6.23)
      Nuj_nonsing = Gaj*(log(Sjk/Q2jk) - 2*log((Qjk-Mk)/Qjk) - 2*Mk/(Qjk+Mk)) + pi2_6 - spence(Sjk/Q2jk)
    else if (M2j > 0 .and. M2k > 0) then ! eq. (6.21)
      Nuj_nonsing = Gaj*log(Sjk/Q2jk) + (log(rho2)*log(1+rho2) + 2*spence(rho2)                      &
                  & - spence(1._qp-rho2j) - spence(1._qp-rho2k) - pi2_6)/Nujk    &
                  & + log((Qjk-Mk)/Qjk) - 2*log(((Qjk-Mk)**2-M2j)/Q2jk) - 2*M2j/Sjk*log(Mj/(Qjk-Mk)) &
                  & - Mk/(Qjk-Mk) + 2*Mk*(2*Mk-Qjk)/Sjk + 3*pi2_6
    else
      call ol_error(2,'subroutine intdip_Fjk: arguments out of range')
      call ol_error(2,'allowed range M2j, M2k >= 0')
!      write(*,*) '[OpenLoops] M2j =',M2j
!      write(*,*) '[OpenLoops] M2k =',M2k
      call ol_fatal()
    end if
  else
    call ol_error(2,'subroutine intdip_Fjk: arguments out of range')
    call ol_error(2,'allowed range flavj=1,2,3,4')
!    write(*,*) '[OpenLoops] flavj =', flavj
    call ol_fatal()
  end if

  ! add log term in last line of (6.16)
  Fjk(0) = Fjk(0) + Nuj_nonsing + Gaj*logmu2_Sjk

  ! write(*,*) 'j,k = ', j, k
  ! write(*,*) 'Nuj_nonsing = ', Nuj_nonsing
  ! write(*,*)

end subroutine intdip_Fjk


! **********************************************************************
function spence(x)
! **********************************************************************
  use ol_dilog_qp, only: Li2
  use kind_types, only: qp
  implicit none
  real(qp), intent(in) :: x
  real(qp)    :: spence
!   complex(qp), external :: olicspenc ! from COLI
  complex(qp) :: dummy_complex
  dummy_complex = x
!   spence = real(olicspenc(dummy_complex,0._qp))
  spence = real(Li2(dummy_complex))
!   write(*,*) 'coli spence', real(olicspenc(dummy_complex,0._qp))
!   write(*,*) 'ol   spence', real(Li2(dummy_complex))
end function spence

end module ol_i_operator_qp

