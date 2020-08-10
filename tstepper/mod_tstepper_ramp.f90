!>
!! \brief mod_tstepper_ramp provides procedures to calculate the time step
!!        during the integration of the transport equations. The time step dt
!!        is changed in accordance to a "ramp" from dt1 to dt2 defined in the
!!        input file.
!!
module mod_tstepper_ramp

   use mod_class_ifile

   implicit none

   ! Makes everything private, except otherwise stated
   private

   ! Public procedures
   public :: tstepper_ramp_init &
      ,      tstepper_dt_ramp

   integer, private :: it1 !< number of iterations up to which dt = dt1
   integer, private :: it2 !< number of iterations from which dt = dt2
   real(8), private :: dt1 !< Initial time step (s)
   real(8), private :: dt2 !< Final time step (s)

contains

   !> \brief Initializes mod_tstepper_ramp
   subroutine tstepper_ramp_init(ifile)
      implicit none
      class(class_ifile), intent(in) :: ifile

      call ifile%get_value(it1, "it1")
      call ifile%get_value(it2, "it2")
      call ifile%get_value(dt1, "dt1")
      call ifile%get_value(dt2, "dt2")

   end subroutine


   !> \brief Calculates the next time step
   real(8) function tstepper_dt_ramp(it)
      implicit none
      integer, intent(in) :: it !< Current iteration

      ! Inner variables
      real(8) :: dt

      if ( it+1 <= It1 ) dt = dt1
      if ( it+1 >  It1 .and. it+1 < It2 ) dt = dt1 + (dt2-dt1)*(it+1-It1)/(It2-It1)
      if ( it+1 >= It2 ) dt = dt2

      tstepper_dt_ramp = dt

   end function

end module
