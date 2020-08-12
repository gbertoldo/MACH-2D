!>
!! \brief Contains procedures to manage errors
!!
module mod_error_handler
   implicit none

contains

   !> \brief Verifies if a variable contains a NaN
   logical function is_nan(x)
      implicit none
      real(8), intent(in) :: x !< Variable

      is_nan = .false.

      if ( ( x * 1.d0 /= x ) .or. ( x + 0.1d0 * x == x ) ) is_nan = .true.

   end function is_nan

end module
