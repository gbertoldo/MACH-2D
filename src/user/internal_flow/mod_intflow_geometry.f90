!>
!! \brief This module defines data structures and procedures related to the
!!        geometry of internal flows.
!!
module mod_intflow_geometry

   use mod_class_ifile

   implicit none

   !> \brief Geometric data to be used in the internal flow calculation
   type, public :: type_intflow_geometry

      real(8) :: Sg  !< Throat area (m2)
      real(8) :: Rcg !< Throat curvature radius (m)

   end type

contains

   !> \brief Calculates the coordinates (x,y) of the north and south boundaries.
   !! The north boundary line is the contour of a nozzle formed by a chamber
   !! of length Lchamb and radius Rin, followed by a convergent section of angle
   !! alf. The junction of chamber and convergent is a circular path of radius
   !! Rc1. The convergent section is connected to the divergent section by two
   !! circular paths. The first one has curvature radius Rc2 and the second one
   !! has curvature radius Rc3. The divergent section has angle bet and ends
   !! with radius Rout.
   subroutine get_grid_boundary_nozzle01(ifile, nx, ny, x, y) ! Output: last two
      implicit none
      class(class_ifile),    intent(in)  :: ifile !< Input file
      integer,               intent(in)  :: nx    !< Number of volumes in the csi direction (real+fictitious)
      integer,               intent(in)  :: ny    !< Number of volumes in the eta direction (real+fictitious)
      real(8), dimension(:), intent(out) :: x     !< Coord. x at the northeast corner of the volume P (m)
      real(8), dimension(:), intent(out) :: y     !< Coord. y at the northeast corner of the volume P (m)

      ! Parameters
      real(8), parameter :: pi = acos(-1.d0)

      ! Inner variables
      integer :: i, j, np ! Dummy indexes
      real(8) :: Rth      ! Throat radius (m)
      real(8) :: Rin      ! Inflow radius (m)
      real(8) :: Rout     ! Outflow radius (m)
      real(8) :: Rc1      ! Radius of curvature bet. cham. and conv. section (m)
      real(8) :: Rc2      ! Radius of curvature bet. conv. section and throat (m)
      real(8) :: Rc3      ! Radius of curvature bet. throat and div. section (m)
      real(8) :: alf      ! Conv. section angle (rad)
      real(8) :: bet      ! Div. section angle (rad)
      real(8) :: Lchamb   ! Length of the chamber (m)

      ! Reference points
      real(8) :: x0, x1, x2, x3, x4, x5, x6
      real(8) :: y0, y1, y2, y3, y4, y5, y6

      ! Getting input parameters

      call ifile%load()
      call ifile%get_value(   Rth, "Rth"   )
      call ifile%get_value(   Rin, "Rin"   )
      call ifile%get_value(  Rout, "Rout"  )
      call ifile%get_value(   Rc1, "Rc1"   )
      call ifile%get_value(   Rc2, "Rc2"   )
      call ifile%get_value(   Rc3, "Rc3"   )
      call ifile%get_value(   alf, "alf"   )
      call ifile%get_value(   bet, "bet"   )
      call ifile%get_value(Lchamb, "Lchamb")

      ! Converting alf and bet to rad
      alf = alf * pi / 180.d0
      bet = bet * pi / 180.d0

      ! Calculating reference points
      call nozzle01_ref_points(Rth, Rin, Rout, Rc1, Rc2, Rc3 & ! Input
         ,                                  alf, bet, Lchamb & ! Input
         ,                        x0, x1, x2, x3, x4, x5, x6 & ! Output
         ,                        y0, y1, y2, y3, y4, y5, y6 ) ! Output


      ! Generating the north boundary

      j = ny-1

      do i = 1, nx-1

         np   = nx * (j-1) + i

         x(np) = xi(i)

         y(np) = nozzle01(x(np), x0, x1, x2, x3, x4, x5, x6 & ! Input
            ,                    y0, y1, y2, y3, y4, y5, y6 & ! Input
            ,                                 Rc1, Rc2, Rc3 ) ! Input

      end do

      ! Generating the south boundary

      j = 1

      do i = 1, nx-1

         np   = nx * (j-1) + i

         x(np) = xi(i)

         y(np) = 0.d0

      end do

   contains

      real(8) function xi(i)
         implicit none
         integer, intent(in) :: i

         xi = ( dble(i-1) / dble(nx-2) ) * (x6-x0) + x0

      end function

   end subroutine



   !> \brief Calculates the contour y(x) for nozzle01
   real(8) function nozzle01(x, x0, x1, x2, x3, x4, x5, x6 & ! Input
         ,                      y0, y1, y2, y3, y4, y5, y6 & ! Input
         ,                                   Rc1, Rc2, Rc3 ) ! Input
      implicit none
      real(8), intent(in)  :: x !< x coord. (m)
      real(8), intent(in)  :: x0, x1, x2, x3, x4, x5, x6
      real(8), intent(in)  :: y0, y1, y2, y3, y4, y5, y6
      real(8), intent(in)  :: Rc1, Rc2, Rc3

      ! Inner variables
      real(8) :: y

      ! Calculating x
      if (      x0 <= x .and. x <= x1  ) then

         y = y0

      else if ( x1 <= x .and. x <= x2  ) then

         y = y1 - Rc1 + sqrt(Rc1**2-(x-x1)**2)

      else if ( x2 <= x .and. x <= x3  ) then

         y = y2 + (y3-y2)/(x3-x2) * (x-x2)

      else if ( x3 <= x .and. x <= x4  ) then

         y = y4 + Rc2 - sqrt(Rc2**2-(x-x4)**2)

      else if ( x4 <= x .and. x <= x5  ) then

         y = y4 + Rc3 - sqrt(Rc3**2-(x-x4)**2)

      else if ( x5 <= x .and. x <= x6  ) then

         y = y5 + (y6-y5)/(x6-x5) * (x-x5)

      else
         write(*,*) "nozzle01:"
         write(*,*) "x out of range. Stopping..."
         stop
      end if

      nozzle01 = y

   end function


   !> \brief Calculates the reference points of contour y(x) for nozzle01
   subroutine nozzle01_ref_points(Rth, Rin, Rout, Rc1, Rc2, Rc3 & ! Input
      ,                                        alf, bet, Lchamb & ! Input
      ,                              x0, x1, x2, x3, x4, x5, x6 & ! Output
      ,                              y0, y1, y2, y3, y4, y5, y6 ) ! Output
      implicit none
      real(8), intent(in)  :: Rth    !< Throat radius (m)
      real(8), intent(in)  :: Rin    !< Inflow radius (m)
      real(8), intent(in)  :: Rout   !< Outflow radius (m)
      real(8), intent(in)  :: Rc1    !< Radius of curvature bet. cham. and conv. section (m)
      real(8), intent(in)  :: Rc2    !< Radius of curvature bet. conv. section and throat (m)
      real(8), intent(in)  :: Rc3    !< Radius of curvature bet. throat and div. section (m)
      real(8), intent(in)  :: alf    !< Conv. section angle (rad)
      real(8), intent(in)  :: bet    !< Div. section angle (rad)
      real(8), intent(in)  :: Lchamb !< Length of the chamber (m)
      real(8), intent(out) :: x0, x1, x2, x3, x4, x5, x6
      real(8), intent(out) :: y0, y1, y2, y3, y4, y5, y6

      ! Inner variables
      real(8) :: ra, rb, la, lb

      ! Calculating reference coordinates

      x0 = 0.d0
      y0 = Rin

      x1 = Lchamb
      y1 = Rin

      x2 = x1 + Rc1 * sin(alf)
      y2 = y1 - Rc1 * (1.d0-cos(alf))

      y4 = Rth
      y3 = y4 + Rc2 * (1.d0-cos(alf))

      ra = y2-y3
      la = ra / tan(alf)

      x3 = x2 + la
      x4 = x3 + Rc2 * sin(alf)

      x5 = x4 + Rc3 * sin(bet)
      y5 = y4 + Rc3 * (1.d0-cos(bet))


      y6 = Rout

      rb = y6-y5
      lb = rb / tan(bet)

      x6 = x5 + lb

   end subroutine


end module
