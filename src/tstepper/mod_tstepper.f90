!>
!! \brief mod_tstepper provides procedures to calculate the time step
!!        during the integration of the transport equations
!!
module mod_tstepper

   use mod_class_ifile
   use mod_tstepper_ramp
   use mod_tstepper_scarborough

   implicit none

   ! Makes everything private, except otherwise stated
   private

   ! Public procedures
   public :: tstepper_init &
      ,      tstepper_dt0  &
      ,      tstepper_dt

   ! Model options
   integer, parameter, private :: TSTEPPER_RAMP        = 0
   integer, parameter, private :: TSTEPPER_SCARBOROUGH = 1

   integer, private:: dtmodel !< Selected time-stepper model (see options above)

contains

   !> \brief Initializes the module tstepper
   subroutine tstepper_init(ifile)
      implicit none
      class(class_ifile), intent(in) :: ifile

      ! Inner variables
      character(len=500) :: caux

      ! Selecting the time-stepper model
      call ifile%get_value(caux, "dtmodel")

      if ( trim(caux) == "RAMP" ) then

         dtmodel = TSTEPPER_RAMP

         call tstepper_ramp_init(ifile)

      else if (  trim(caux) == "SCARBOROUGH" ) then

         dtmodel = TSTEPPER_SCARBOROUGH

         call tstepper_scarborough_init(ifile)

      else

         write(*,*) "tstepper_init:"
         write(*,*) "Unknown model: ", trim(caux)
         write(*,*) "Stopping..."
         stop

      end if

   end subroutine


   !> \brief Calculates the first time step
   real(8) function tstepper_dt0(it, nx, ny, Jp, u, v)
      implicit none
      integer,                    intent(in) :: it !< Current iteration
      integer,                    intent(in) :: nx !< Number of volumes in csi direction (real+fictitious)
      integer,                    intent(in) :: ny !< Number of volumes in eta direction (real+fictitious)
      real(8), dimension (nx*ny), intent(in) :: Jp !< Jacobian at the center of volume P
      real(8), dimension (nx*ny), intent(in) :: u  !< Cartesian velocity of the last iteraction
      real(8), dimension (nx*ny), intent(in) :: v  !< Cartesian velocity of the last iteraction

      select case (dtmodel)

         case (TSTEPPER_RAMP)

            tstepper_dt0 = tstepper_dt_ramp(it)

         case (TSTEPPER_SCARBOROUGH)

            tstepper_dt0 = tstepper_dt0_scarborough(nx, ny, Jp, u, v)

         case default

            tstepper_dt0 = 1.d-10

      end select

   end function


   !> \brief Calculates the next time step
   real(8) function tstepper_dt(it, nx, ny, dt, au, av, at, ap)
      implicit none
      integer, intent(in) :: it     !< Current iteration
      integer, intent(in) :: nx   !< Number of volumes in the csi direction (real+fictitious)
      integer, intent(in) :: ny   !< Number of volumes in the eta direction (real+fictitious)
      real(8), intent(in) :: dt   !< Old time step (s)
      real(8), dimension (nx*ny,9), intent(in) :: au !< Coefficients of the u linear system
      real(8), dimension (nx*ny,9), intent(in) :: av !< Coefficients of the v linear system
      real(8), dimension (nx*ny,9), intent(in) :: at !< Coefficients of the T linear system
      real(8), dimension (nx*ny,5), intent(in) :: ap !< Coefficients of the p linear system

      ! Selecting the dt calculator
      select case (dtmodel)

         case (TSTEPPER_RAMP)

            tstepper_dt = tstepper_dt_ramp(it)

         case (TSTEPPER_SCARBOROUGH)

            tstepper_dt = tstepper_dt_scarborough(nx, ny, dt, au, av, at, ap)

         case default

            tstepper_dt = 1.d-10

      end select

   end function

end module
