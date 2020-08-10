!>
!! \brief mod_tstepper_scarborough provides procedures to calculate the time
!!        step during the integration of the transport equations. The time step
!!        dt is changed dynamically to reach a defined interval of the convergence
!!        coefficient that satisfies the Scarborough criteria.
!!
module mod_tstepper_scarborough

   use mod_class_ifile

   implicit none

   ! Makes everything private, except otherwise stated
   private

   ! Public procedures
   public :: tstepper_scarborough_init &
      ,      tstepper_dt0_scarborough  &
      ,      tstepper_dt_scarborough

   real(8), private :: h0    !< Amplitude of h in the TSI11 SCARBOROUGH model
   real(8), private :: mincc !< Minimum allowed value of the convergence coefficient for SCARBOROUGH model
   real(8), private :: maxcc !< Maximum allowed value of the convergence coefficient for SCARBOROUGH model

contains

   !> \brief Initializes mod_tstepper_scarborough
   subroutine tstepper_scarborough_init(ifile)
      implicit none
      class(class_ifile), intent(in) :: ifile

      call ifile%get_value(   h0, "h0"   )
      call ifile%get_value(mincc, "mincc")
      call ifile%get_value(maxcc, "maxcc")

   end subroutine


   !> \brief Calculates the first time step
   real(8) function tstepper_dt0_scarborough(nx, ny, Jp, u, v)
      implicit none
      integer, intent(in) :: nx !< Number of volumes in csi direction (real+fictitious)
      integer, intent(in) :: ny !< Number of volumes in eta direction (real+fictitious)
      real(8), dimension (nx*ny), intent(in) :: Jp !< Jacobian at the center of volume P
      real(8), dimension (nx*ny), intent(in) :: u  !< Cartesian velocity of the last iteraction
      real(8), dimension (nx*ny), intent(in) :: v  !< Cartesian velocity of the last iteraction

      ! Inner variables

      integer :: i, j, np  ! Auxiliary variables
      real(8) :: dt        ! Time step
      real(8) :: aux

      ! Searching for the smallest dt

      dt = 1.d10

      do i = 2, nx-1

         do j = 2, ny-1

            np  = nx * (j-1) + i

            aux = 1.d0 / sqrt( Jp(np) * ( u(np) ** 2 + v(np) ** 2 ) )

            if ( dt > aux ) dt = aux

         end do

      end do

      dt = dt * 0.25d0

      tstepper_dt0_scarborough = dt

   end function


   !> \brief Calculates the next time step
   real(8) function tstepper_dt_scarborough(nx, ny, dt, au, av, at, ap)
      implicit none
      integer, intent(in) :: nx   !< Number of volumes in the csi direction (real+fictitious)
      integer, intent(in) :: ny   !< Number of volumes in the eta direction (real+fictitious)
      real(8), intent(in) :: dt   !< Old time step (s)
      real(8), dimension (nx*ny,9), intent(in) :: au !< Coefficients of the u linear system
      real(8), dimension (nx*ny,9), intent(in) :: av !< Coefficients of the v linear system
      real(8), dimension (nx*ny,9), intent(in) :: at !< Coefficients of the T linear system
      real(8), dimension (nx*ny,5), intent(in) :: ap !< Coefficients of the p linear system

      ! Inner variables
      real(8) :: h     ! Relative variation of dt
      real(8) :: c     ! Coefficient of convergence
      real(8) :: curef ! Coefficient of convergence for u
      real(8) :: cvref ! Coefficient of convergence for v
      real(8) :: ctref ! Coefficient of convergence for T
      real(8) :: cpref ! Coefficient of convergence for p

      ! Calculation of the convergence coefficients

      curef = get_mean_cc_9d(nx, ny, au)
      cvref = get_mean_cc_9d(nx, ny, av)
      ctref = get_mean_cc_9d(nx, ny, at)
      cpref = get_mean_cc_5d(nx, ny, ap)

      c = max(curef, cvref, ctref, cpref )

      if ( c < mincc) then

         h = h0 * sqrt( 1.d0 - c / mincc )

      else if ( maxcc < c ) then

         h = - h0 * sqrt( c / maxcc - 1.d0 )

      else

         h = 0.d0

      end if

      tstepper_dt_scarborough = dt * ( 1.d0 + h )

   end function


   !> \brief Calculates the coefficient for convergence criteria of a 9-diagonal matrix.
   !! This coefficient is calculated only for real volumes.
   real(8) function get_maxc_9d(nx, ny, a)
      implicit none
      integer, intent(in) :: nx   !< Number of volumes in the csi direction (real+fictitious)
      integer, intent(in) :: ny   !< Number of volumes in the eta direction (real+fictitious)
      real(8), dimension (nx*ny,9), intent(in) :: a   !< Coefficients of the linear system

      ! Inner variables

      integer :: i, j, np

      real(8) :: c, maxc

      ! Looking for the maxval of c = ( sum_nb | A_nb | ) / | A_P |

      maxc = 0.d0

      do i = 2, nx-1

         do j = 2, ny-1

            np   = nx * (j-1) + i

            c = ( sum(abs(a(np,1:4))) + sum(abs(a(np,6:9))) ) / abs( a(np,5) )

            if ( maxc < c ) maxc = c

         end do

      end do

      get_maxc_9d = maxc

   end function


   !> \brief Calculates the coefficient for convergence criteria of a 5-diagonal matrix.
   !! This coefficient is calculated only for real volumes.
   real(8) function get_maxc_5d(nx, ny, a)
      implicit none
      integer, intent(in) :: nx   !< Number of volumes in the csi direction (real+fictitious)
      integer, intent(in) :: ny   !< Number of volumes in the eta direction (real+fictitious)
      real(8), dimension (nx*ny,5), intent(in) :: a   !< Coefficients of the linear system

      ! Inner variables

      integer :: i, j, np

      real(8) :: c, maxc

      ! Looking for the maxval of c = ( sum_nb | A_nb | ) / | A_P |

      maxc = 0.d0

      do i = 2, nx-1

         do j = 2, ny-1

            np   = nx * (j-1) + i

            c = ( sum(abs(a(np,1:2))) + sum(abs(a(np,4:5))) ) / abs( a(np,3) )

            if ( maxc < c ) maxc = c

         end do

      end do

      get_maxc_5d = maxc

   end function


   !> \brief Calculates for each real volume the convergence coefficient
   !! associated with the 9-diagonal matrix 'a' and stores the results in cc.
   subroutine get_cc_9d(nx, ny, a, cc)
      implicit none
      integer, intent(in) :: nx   !< Number of volumes in the csi direction (real+fictitious)
      integer, intent(in) :: ny   !< Number of volumes in the eta direction (real+fictitious)
      real(8), dimension (nx*ny,9), intent(in)  :: a   !< Coefficients of the linear system
      real(8), dimension (nx*ny),   intent(out) :: cc  !< Convergence coef. vector

      ! Inner variables

      integer :: i, j, np

      do i = 2, nx-1

         do j = 2, ny-1

            np = nx * (j-1) + i

            cc(np) = ( sum( abs(a(np,1:4)) ) + sum( abs(a(np,6:9)) ) ) &

            / abs(a(np,5))

         end do

      end do

   end subroutine

   !> \brief Calculates for each real volume the convergence coefficient
   !! associated with the 5-diagonal matrix 'a' and stores the results in cc.
   subroutine get_cc_5d(nx, ny, a, cc)
      implicit none
      integer, intent(in) :: nx   !< Number of volumes in the csi direction (real+fictitious)
      integer, intent(in) :: ny   !< Number of volumes in the eta direction (real+fictitious)
      real(8), dimension (nx*ny,5), intent(in)  :: a   !< Coefficients of the linear system
      real(8), dimension (nx*ny),   intent(out) :: cc  !< Convergence coef. vector

      ! Inner variables

      integer :: i, j, np

      do i = 2, nx-1

         do j = 2, ny-1

            np = nx * (j-1) + i

            cc(np) = ( sum( abs(a(np,1:2)) ) + sum( abs(a(np,4:5)) ) ) &

            / abs(a(np,3))

         end do

      end do

   end subroutine

   !> \brief Calculates the mean value of the convergence coefficient
   !! associated with the 9-diagonal matrix 'a'.
   real(8) function get_mean_cc_9d(nx, ny, a)
      implicit none
      integer, intent(in) :: nx   !< Number of volumes in the csi direction (real+fictitious)
      integer, intent(in) :: ny   !< Number of volumes in the eta direction (real+fictitious)
      real(8), dimension (nx*ny,9), intent(in)  :: a   !< Coefficients of the linear system

      ! Inner variables

      integer :: i, j, np

      get_mean_cc_9d = 0.d0

      do i = 2, nx-1

         do j = 2, ny-1

            np = nx * (j-1) + i

            get_mean_cc_9d = ( sum( abs(a(np,1:4)) ) + sum( abs(a(np,6:9)) ) ) &

            / abs(a(np,5)) + get_mean_cc_9d

         end do

      end do

      get_mean_cc_9d = get_mean_cc_9d / dble( (nx-2) * (ny-2) )

   end function

   !> \brief Calculates the mean value of convergence coefficient
   !! associated with the 5-diagonal matrix 'a'.
   real(8) function get_mean_cc_5d(nx, ny, a)
      implicit none
      integer, intent(in) :: nx   !< Number of volumes in the csi direction (real+fictitious)
      integer, intent(in) :: ny   !< Number of volumes in the eta direction (real+fictitious)
      real(8), dimension (nx*ny,5), intent(in)  :: a   !< Coefficients of the linear system

      ! Inner variables

      integer :: i, j, np

      get_mean_cc_5d = 0.d0

      do i = 2, nx-1

         do j = 2, ny-1

            np = nx * (j-1) + i

            get_mean_cc_5d = ( sum( abs(a(np,1:2)) ) + sum( abs(a(np,4:5)) ) ) &

            / abs(a(np,3)) + get_mean_cc_5d

         end do

      end do

      get_mean_cc_5d = get_mean_cc_5d / dble( (nx-2) * (ny-2) )

   end function


end module
