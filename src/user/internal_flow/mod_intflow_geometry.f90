!>
!! \brief This module defines data structures and procedures related to the
!!        geometry of internal flows.
!!
module mod_intflow_geometry

   implicit none

   !> \brief Geometric data to be used in the internal flow calculation
   type, public :: type_intflow_geometry

      real(8) :: Sg  !< Throat area (m2)
      real(8) :: Rcg !< Throat curvature radius (m)

   end type

  public:: set_grid_internal_flow

contains

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
