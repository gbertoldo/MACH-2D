!>
!! \brief mod_grid_procedures contains procedures to generate grid
!!
module mod_grid_procedures

   use thompson2d_hyperbolic, only: get_thompson2d_hyperbolic_01

   implicit none

contains

   !> \brief Generates the grid provided the south and north boundaries of the domain
   subroutine set_grid(isReversed, kg, nx, ny, a1, avi, avf, awf, x, y) ! Last 2 are inoutput
      implicit none
      logical, intent(in)    :: isReversed !< Reverses the discretization distribution if true
      integer, intent(in)    :: kg         !< Kind of grid (1=uniform, 2=geometric progression, 3=power law, 4=gp modified)
      integer, intent(in)    :: nx         !< Number of volumes in csi direction (real+fictitious)
      integer, intent(in)    :: ny         !< Number of volumes in eta direction (real+fictitious)
      real(8), intent(in)    :: a1         !< Width of the vol. closer to the wall
      real(8), intent(in)    :: avi        !< Initial value of the artificial viscosity
      real(8), intent(in)    :: avf        !< Final value of the artificial viscosity
      real(8), intent(in)    :: awf        !< Area weighting factor
      real(8), intent(inout) :: x(nx*ny)   !< x coord. at the northeast corner of the volume P
      real(8), intent(inout) :: y(nx*ny)   !< y coord. at the northeast corner of the volume P

      ! Selecting the grid
      if ( kg == 1 ) then

         call get_uniform_grid( nx, ny, x, y) ! Last 2 are inoutput

      else if ( kg == 2 ) then

         call get_pg_grid(isReversed, a1, nx, ny, x, y) ! Last 2 are inoutput

      else if ( kg == 3 ) then

         call get_power_grid(isReversed, a1, nx, ny, x, y) ! Last 2 are inoutput

      else if ( kg == 4 ) then

         call get_pg_grid2(isReversed, a1, nx, ny, x, y) ! Last 2 are inoutput

      else if ( kg == 5 ) then

         call get_thompson2d_hyperbolic_01(nx, ny, avi, avf, awf, a1, x, y) ! Inoutput: last two

      else

         write(*,*) 'set_grid: Type of grid unknown'
         stop

      end if

   end subroutine


   !> \brief Calculates a normalized geometric progression distribution t
   !! t(1)=0, t(N)=1. Points are concentrated near t=0 by default, except
   !! if isReversed is true.
   subroutine get_gp_distribution(isReversed, N, r, t) ! Output: last one
      implicit none
      logical, intent(in)  :: isReversed !< Reverses the discretization distribution if true
      integer, intent(in)  :: N          !< Number of points in the interval
      real(8), intent(in)  :: r          !< Ratio of length to first partition width
      real(8), intent(out) :: t(1:N)     !< Distribution

      integer :: j ! Dummy index
      real(8) :: q ! GP ratio

      call get_GP_ratio(N-1,r,q) ! Last one is output

      t(1) = 0.d0

      do j = 2, N-1

         t(j) = t(j-1) + q**(j-2) / r

      end do

      t(N) = 1.d0

      ! Should be reversed?
      if (isReversed) call reverse(t)

   end subroutine


   !> \brief Calculates a normalized geometric power law distribution t
   !! t(1)=0, t(N)=1. Points are concentrated near t=0 by default, except
   !! if isReversed is true.
   subroutine get_power_law_distribution(isReversed, N, r, t) ! Output: last one
      implicit none
      logical, intent(in)  :: isReversed !< Reverses the discretization distribution if true
      integer, intent(in)  :: N          !< Number of points in the interval
      real(8), intent(in)  :: r          !< Ratio of length to first partition width
      real(8), intent(out) :: t(1:N)     !< Distribution

      integer :: j ! Dummy index
      real(8) :: a ! Exponent of the power law

      if ( r > 1.d0 ) then

         a = log( r ) / log(dble(N-1))

      else

         a = 1.d0

      end if

      t(1) = 0.d0

      do j = 2, N-1

         t(j) = (dble(j-1)/dble(N-1))**a

      end do

      t(N) = 1.d0

      ! Should be reversed?
      if (isReversed) call reverse(t)

   end subroutine


   !> \brief Reverses a partition distribution
   subroutine reverse(t)
      implicit none
      real(8), intent(inout) :: t(:) !< Distribution of points

      integer :: i, N
      real(8) :: dt(size(t)-1)

      ! Number of partitions
      N = size(t)-1

      ! Copying the partition distribution in reversed order
      do i = 1, N
         dt(N-i+1) = t(i+1)-t(i)
      end do

      ! Reassigning t
      do i = 2, N
         t(i)=t(i-1)+dt(i-1)
      end do

   end subroutine


   !> \brief Generates a uniform partitioned grid
   subroutine get_uniform_grid(nx, ny, x, y) ! Last 2 are inoutput
      implicit none
      integer, intent(in)    :: nx
      integer, intent(in)    :: ny
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


   !> \brief Generates a grid with a non-uniform power-law distribution
   subroutine get_power_grid(isReversed, a1, nx, ny, x, y) ! Last 2 are inoutput
      implicit none
      logical, intent(in)    :: isReversed !< Reverses the discretization distribution if true
      integer, intent(in)    :: nx
      integer, intent(in)    :: ny
      real(8), intent(in)    :: a1
      real(8), intent(inout) :: x(nx*ny)
      real(8), intent(inout) :: y(nx*ny)
      !
      integer :: i, j, np
      real(8) :: r, xi, xf, yi, yf
      real(8) :: tv(ny-1)

      do i = 1, nx-1

         j = 1
         np = nx * (j-1) + i
         xi = x(np)
         yi = y(np)

         j = ny-1
         np = nx * (j-1) + i
         xf = x(np)
         yf = y(np)

         r = sqrt( (xf-xi)**2 + (yf-yi)**2 ) / a1

         ! Getting normalized distribution
         call get_power_law_distribution(isReversed, ny-1, r, tv) ! Output: last one

         do j = 2, ny-2

            np  = nx * (j-1) + i

            x(np) = (xf-xi) * tv(j) + xi

            y(np) = (yf-yi) * tv(j) + yi

         end do

      end do

   end subroutine get_power_grid


   !> \brief Generates a grid with a non-uniform geometric progression dist.
   subroutine get_pg_grid(isReversed, a1, nx, ny, x, y) ! Last 2 are inoutput
      implicit none
      logical, intent(in)    :: isReversed !< Reverses the discretization distribution if true
      integer, intent(in)    :: nx
      integer, intent(in)    :: ny
      real(8), intent(in)    :: a1
      real(8), intent(inout) :: x(nx*ny)
      real(8), intent(inout) :: y(nx*ny)
      !
      integer :: i, j, np, npn
      real(8) :: r, xi, xf, yi, yf

      real(8) :: tv(ny-1)

      do i = 1, nx-1

         j = 1
         np = nx * (j-1) + i
         xi = x(np)
         yi = y(np)

         j = ny-1
         np = nx * (j-1) + i
         xf = x(np)
         yf = y(np)

         r = sqrt( (xf-xi)**2 + (yf-yi)**2 ) / a1

         call get_gp_distribution(isReversed, ny-1, r, tv) ! Output: last one

         do j = 1, ny-3

            np  = nx * (j-1) + i

            npn  = np + nx

            x(npn) = xi + (xf-xi) * tv(j+1)

            y(npn) = yi + (yf-yi) * tv(j+1)

         end do

      end do

   end subroutine get_pg_grid


   !> \brief Generates a grid with a non-uniform geometric progression of type 2 dist.
   subroutine get_pg_grid2(isReversed, a1, nx, ny, x, y) ! Last 2 are inoutput
      implicit none
      logical, intent(in)    :: isReversed !< Reverses the discretization distribution if true
      integer, intent(in)    :: nx
      integer, intent(in)    :: ny
      real(8), intent(in)    :: a1
      real(8), intent(inout) :: x(nx*ny)
      real(8), intent(inout) :: y(nx*ny)
      !
      integer :: i, j, np, npn
      real(8) :: r, xi, xf, yi, yf, lf, l, a
      real(8) :: tv(ny-1)

      i = nx-1
      j = 1

      np = nx * (j-1) + i
      xi = x(np)
      yi = y(np)

      j = ny-1
      np = nx * (j-1) + i
      xf = x(np)
      yf = y(np)

      lf = sqrt( (xf-xi)**2 + (yf-yi)**2 )


      do i = 1, nx-1

         j = 1
         np = nx * (j-1) + i
         xi = x(np)
         yi = y(np)

         j = ny-1
         np = nx * (j-1) + i
         xf = x(np)
         yf = y(np)

         l = sqrt( (xf-xi)**2 + (yf-yi)**2 )

         a = a1 * l / lf

         r = sqrt( (xf-xi)**2 + (yf-yi)**2 ) / a

         call get_gp_distribution(isReversed, ny-1, r, tv) ! Output: last one

         do j = 1, ny-3

            np  = nx * (j-1) + i

            npn  = np + nx

            x(npn) = xi + (xf-xi) * tv(j+1)

            y(npn) = yi + (yf-yi) * tv(j+1)

         end do

      end do

   end subroutine get_pg_grid2


   !> \brief Calculates the geometric progression ratio
   subroutine get_gp_ratio(n, r, q)
     implicit none
     integer, intent(in)  ::   n !< number of partitions
     real(8), intent(in)  ::   r !< l/a1
     real(8), intent(out) ::   q !< q

     ! Parameters

     integer :: nit = 1000   ! Maximum number of iteractions
     real(8) :: tol = 1.d-15 ! Tolerance

     ! Inner variables

     integer ::   i ! Dummy index
     real(8) ::  qi ! inital value of q
     real(8) ::  qf ! final value of q
     real(8) ::  qm ! mean value of q

     if ( r < n ) then

        qi = 0.1d0

        qf = 1.d0 - 1.d-15

     else

        qi = 1.d0 + 1.d-15

        qf = 10.d0

     end if

     do i = 1, nit

        qm = 0.5d0 * qi + 0.5d0 * qf

        if ( 0.d0 < f(qi) * f(qm) ) then

           qi = qm

        else

           qf = qm

        end if

        if ( abs(qf-qi) < tol ) exit

     end do


     if ( i == nit ) then

        write(*,*) "get_gp_ratio: Maximum number of iteractions was exceeded."

        stop

     end if

     q = qm

   contains

     real(8) function f(q)
       implicit none
       real(8), intent(in) :: q

       f = q ** n + r * ( 1.d0 - q ) - 1.d0

     end function f

   end subroutine get_gp_ratio


end module
