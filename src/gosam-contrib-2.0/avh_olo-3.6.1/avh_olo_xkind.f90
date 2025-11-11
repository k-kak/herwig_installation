module avh_olo_xkind
  implicit none
  private
  public :: olo_xkind_dp, olo_xkind_qp
  !
  integer ,parameter :: olo_xkind_dp=kind(1.0d0)
  !
  integer ,parameter :: olo_xkind_qp=selected_real_kind(32,50)
  !
end module
