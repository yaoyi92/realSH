module sqrttable
  implicit none
  real*8, allocatable :: LMtable(:)

  contains

  subroutine init_sqrttable(max_l)
    implicit none
    integer, intent(in) ::max_l
    integer max_l_local
    integer i_int
    max_l_local = max(max_l, 10)
    allocate(LMtable(4*max_l*max_l-1))
    do i_int = 1, 4*max_l*max_l-1
      LMtable(i_int) = sqrt(dble(i_int))
    end do
  end subroutine

  subroutine finalize_sqrttable()
    implicit none
    deallocate(LMtable)
  end subroutine


end module sqrttable
