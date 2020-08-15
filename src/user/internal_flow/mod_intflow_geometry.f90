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

  public:: set_grid_internal_flow

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



  subroutine set_grid_internal_flow(kg, nx, ny, a1, x, y) ! Last 2 are output
    implicit none
    integer, intent(in)  :: kg ! Kind of grid (1=uniform, 2=geometric progression, 3=power law)
    integer, intent(in)  :: nx
    integer, intent(in)  :: ny
    real(8), intent(in)  :: a1
    real(8), intent(out) :: x(nx*ny)
    real(8), intent(out) :: y(nx*ny)
    !
    integer :: i, j, np

    ! Read from a file the boundary of the grid
    open(10, file="./mach2d_input/gridboundary.dat")

    ! Reading lower boundary
    j = 1
    read(10,*)
    do i = 1, nx-1
       np = nx * (j-1) + i
       read(10,*) x(np), y(np)
    end do

    ! Reading upper boundary
    read(10,*)
    j = ny - 1
    do i = 1, nx-1
       np = nx * (j-1) + i
       read(10,*) x(np), y(np)
    end do

    if ( kg == 1 ) then

       call get_uniform_grid( nx, ny, x, y)

    else if ( kg == 2 ) then

       call get_backward_GP_grid(a1, nx, ny, x, y)

    else if ( kg == 3 ) then

        call get_powerlaw_grid(a1, nx, ny, x, y) ! Last 2 are inoutput

    else

       write(*,*) 'set_grid: Type of grid unknown'
       stop

    end if

  end subroutine set_grid_internal_flow

  subroutine get_uniform_grid( nx, ny, x, y) ! Last 2 are inoutput
    implicit none
    integer, intent(in)  :: nx
    integer, intent(in)  :: ny
    real(8), intent(inout) :: x(nx*ny)
    real(8), intent(inout) :: y(nx*ny)
    !
    integer :: i, j, np, nps
    real(8) :: xi, xf, yi, yf, dx, dy

    do i = 1, nx-1

       j = 1
       np = nx * (j-1) + i
       xi = x(np)
       yi = y(np)

       j = ny-1
       np = nx * (j-1) + i
       xf = x(np)
       yf = y(np)

       dx = (xf-xi) / (ny-2)
       dy = (yf-yi) / (ny-2)

       do j = 2, ny-2
          np  = nx * (j-1) + i
          nps = np - nx
          x(np) = x(nps) + dx
          y(np) = y(nps) + dy
       end do
    end do

  end subroutine get_uniform_grid

  subroutine get_backward_GP_grid(a1, nx, ny, x, y) ! Last 2 are inoutput
    implicit none
    integer, intent(in)  :: nx
    integer, intent(in)  :: ny
    real(8), intent(in)  :: a1
    real(8), intent(inout) :: x(nx*ny)
    real(8), intent(inout) :: y(nx*ny)
    !
    integer :: i, j, np, npn
    real(8) :: q, xi, xf, yi, yf

    do i = 1, nx-1

       j = 1
       np = nx * (j-1) + i
       xi = x(np)
       yi = y(np)

       j = ny-1
       np = nx * (j-1) + i
       xf = x(np)
       yf = y(np)

       call get_GP_ratio( ny-2, (yf-yi)/a1, q)

       do j = ny-2, 2, -1

          np  = nx * (j-1) + i
          npn = np + nx

          x(np) = x(npn)

          y(np) = y(npn) - a1 * q ** (ny-j-2)

       end do
    end do

  end subroutine get_backward_GP_grid

  subroutine get_GP_ratio(n,r,q) ! Last one is outpu

    implicit none
    integer, intent(in)  :: n ! number of divisions
    real(8), intent(in)  :: r ! lenght l to a1 ratio
    real(8), intent(out) :: q ! PG ratio
    !
    integer :: it
    real(8) :: qo, fo, fo1

    qo = 2.d0

    if ( dble(n) > r ) qo = 0.5d0

    it = 0

    do

       it = it + 1

       fo  = qo**n + r*(1.d0-qo) - 1.d0

       fo1 = dble(n) * qo**(n-1) - r

       q = qo - fo/fo1

       if ( abs(q-qo) < 1.d-15 ) exit

       if ( it > 10000 ) then

          write(*,*) "GetPGRatio: Maximum number of iteractions was exceeded."

          stop

       end if

       qo = q

    end do

  end subroutine Get_GP_ratio


   subroutine get_powerlaw_grid(a1, nx, ny, x, y) ! Last 2 are inoutput
    implicit none
    integer, intent(in)  :: nx
    integer, intent(in)  :: ny
    real(8), intent(in)  :: a1
    real(8), intent(inout) :: x(nx*ny)
    real(8), intent(inout) :: y(nx*ny)
    !
    integer :: i, j, np, nps
    real(8) :: xi, xf, yi, yf, dx, dy
    real(8) :: alpha

    do i = 1, nx-1

       j = 1
       np = nx * (j-1) + i
       xi = x(np)
       yi = y(np)

       j = ny-1
       np = nx * (j-1) + i
       xf = x(np)
       yf = y(np)

       dx = (xf-xi) / (ny-2)

       dy = (yf-yi)

       alpha = log( dy / a1 ) / log( dble(ny-2) )

       do j = 2, ny-2

          np  = nx * (j-1) + i

          nps = np - nx

          x(np) = x(nps) + dx

          y(np) = yi + dy * ( 1.d0 - ( dble(ny-1-j) / (ny-2) )**alpha )

       end do
    end do

   end subroutine get_powerlaw_grid

end module
