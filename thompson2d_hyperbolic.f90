!> \brief Contains subroutines for the generation of grids based on hyperbolic
!! equations
module thompson2d_hyperbolic
   implicit none

contains

   !> \brief Generates a grid using the generation model 01.
   subroutine get_thompson2d_hyperbolic_01(nx, ny, avi, avf, awf, a1, x, y) ! Inoutput: last two
      implicit none
      integer, intent(in)    :: nx       !< Number of volumes (real+fictitious) in the csi direction
      integer, intent(in)    :: ny       !< Number of volumes (real+fictitious) in the eta direction
      real(8), intent(in)    :: avi      !< Initial value of the artificial viscosity
      real(8), intent(in)    :: avf      !< Final value of the artificial viscosity
      real(8), intent(in)    :: awf      !< Area weighting factor
      real(8), intent(in)    :: a1       !< Width of the volumes closer to the south boundary
      real(8), intent(inout) :: x(nx*ny) !< x-coordinate
      real(8), intent(inout) :: y(nx*ny) !< y-coordinate

      ! Inner variables

      integer :: i, j, np, npi, npf ! Dummy indexes
      real(8) :: q        ! Ratio q of the geometric progression
      real(8) :: lw       ! Width of the domain
      real(8) :: an       ! Width of the volumes according to a geometric progression
      real(8) :: avc      ! Artificial viscosity
      real(8) :: btw      ! Angle of the west normal vector
      real(8) :: bte      ! Angle of the east normal vector
      real(8) :: xbs(nx)  ! x coordinate of the first eta line
      real(8) :: ybs(nx)  ! y coordinate of the first eta line
      real(8) :: xbn(nx)  ! x coordinate of the next eta line
      real(8) :: ybn(nx)  ! y coordinate of the next eta line


      ! Calculating btw and bte

      ! West

      i = 1
      j = 1
      npi = nx * (j-1) + i

      j = ny-1
      npf = nx * (j-1) + i

      btw = atan2( y(npf)-y(npi), x(npf)-x(npi) )

      ! East

      i = nx-1
      j = 1
      npi = nx * (j-1) + i

      j = ny-1
      npf = nx * (j-1) + i

      bte = atan2( y(npf)-y(npi), x(npf)-x(npi) )


      ! Width of the domain

      lw = sqrt( (y(npf)-y(npi))**2 + (x(npf)-x(npi))**2 )


      ! Calculating the geometric progression ratio

      call get_gp_ratio(ny-2, lw/a1, q)

      ! Extracting the eta line

      j = 1

      do i = 1, nx-1

         np   = nx * (j-1) + i

         xbn(i) = x(np)

         ybn(i) = y(np)

      end do


      ! Building the eta lines

      do j = 2, ny-1

         ! Updating the eta line

         xbs = xbn

         ybs = ybn

         ! Calculating the artificial viscosity

         avc = avi + (avf-avi) * dble(j-1) / dble(ny-2)

         ! Calculating the stepping size

         an = a1 * q ** (j-2)

         ! Calculating the next eta line

         call get_thompson2d_hyperbolic_01_step(nx, btw, bte, avc, awf, an &
         , xbs, ybs, xbn, ybn) ! Output: last two


         ! Updating the grid

         do i = 1, nx-1

            np   = nx * (j-1) + i

            x(np) = xbn(i)

            y(np) = ybn(i)

         end do

      end do

   contains

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

           if ( 0.d0 < f(n, r, qi) * f(n, r, qm) ) then

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

      end subroutine get_gp_ratio

      real(8) function f(n, r, q)
         implicit none
         integer, intent(in) :: n
         real(8), intent(in) :: r
         real(8), intent(in) :: q

         f = q ** n + r * ( 1.d0 - q ) - 1.d0

      end function f

   end subroutine get_thompson2d_hyperbolic_01


   !> \brief Given an eta line, this subroutine creates another one.
   subroutine get_thompson2d_hyperbolic_01_step(nx, btw, bte, avc, awf, hs &
         , xbs, ybs, xbn, ybn) ! Output: last two
      implicit none
      integer, intent(in)  :: nx       !< Number of volumes (real+fictitious) in the csi direction
      real(8), intent(in)  :: btw      !< Angle of the normal vector at the west boundary
      real(8), intent(in)  :: bte      !< Angle of the normal vector at the east boundary
      real(8), intent(in)  :: avc      !< Numerical viscosity coefficient
      real(8), intent(in)  :: awf      !< Area weighting factor
      real(8), intent(in)  :: hs       !< Step from one eta-line to another
      real(8), intent(in)  :: xbs(nx)  !< x-coordinates on the south boundary
      real(8), intent(in)  :: ybs(nx)  !< y-coordinates on the south boundary
      real(8), intent(out) :: xbn(nx)  !< x-coordinates on the north boundary
      real(8), intent(out) :: ybn(nx)  !< y-coordinates on the north boundary

      ! Inner variables

      integer :: i
      real(8) :: pi
      real(8) :: ap(nx-1,3) ! Coefficients of the linear system
      real(8) :: bp(nx-1)   ! Source of the linear system
      real(8) :: bt(nx-1)   ! Distribution of angles of the normal vectors
      real(8) :: al(nx-1)   ! Distribution of angles of the tangent vectors
      real(8) :: sf(nx-1)   ! Distribution of step factors

      ! Defining pi

      pi = acos(-1.d0)

      ! Defining the angles (alpha) of the tangent vectors

      ! West

      i = 1

      al(i) = atan2( ybs(i+1)-ybs(i), xbs(i+1)-xbs(i) )


      ! Central
      do i = 2, nx-2

         al(i) = atan2( ybs(i+1)-ybs(i-1), xbs(i+1)-xbs(i-1) )

      end do


      ! East

      i = nx-1

      al(i) = atan2( ybs(i)-ybs(i-1), xbs(i)-xbs(i-1) )



      ! Defining the coefficients and source for the linear system of beta

      ! West

      i = 1

      ap(i,1) = 0.d0
      ap(i,2) = 1.d0
      ap(i,3) = 0.d0
      bp(i) = btw

      ! Central

      do i = 2, nx-2

         ap(i,1) = -avc
         ap(i,2) = 1.d0 + 2.d0 * avc
         ap(i,3) = -avc
         bp(i) = pi / 2.d0 + al(i)

      end do

      ! East

      i = nx-1

      ap(i,1) = 0.d0
      ap(i,2) = 1.d0
      ap(i,3) = 0.d0
      bp(i) = bte



      ! Solving the linear system for beta

      call tdma(nx-1, ap, bp, bt)


      ! Calculating the stepping factor

      ! West

      sf(1) = 1.d0

      ! Central

      do i = 2, nx-2

         sf(i) = 1.d0 / dabs( dsin( bt(i) - al(i) ) )

      end do

      ! East

      sf(nx-1) = 1.d0


      ! Smoothing the stepping factor

      ! West

      i = 1

      ap(i,1) = 0.d0
      ap(i,2) = 1.d0
      ap(i,3) = 0.d0
      bp(i) = 1.d0

      ! Central

      do i = 2, nx-2

         ap(i,1) = ( awf - 1.d0 ) / 2.d0
         ap(i,2) = 1.d0
         ap(i,3) = ( awf - 1.d0 ) / 2.d0
         bp(i) = awf * sf(i)

      end do

      ! East

      i = nx-1

      ap(i,1) = 0.d0
      ap(i,2) = 1.d0
      ap(i,3) = 0.d0
      bp(i) = 1.d0


      ! Solving the linear system for sf

      call tdma(nx-1, ap, bp, sf)


      ! Calculating the next eta line

      do i = 1, nx-1

         xbn(i) = xbs(i) + hs * cos(bt(i)) * sf(i)

         ybn(i) = ybs(i) + hs * sin(bt(i)) * sf(i)

      end do

   end subroutine get_thompson2d_hyperbolic_01_step



   !> \brief TDMA subroutine
   subroutine tdma(n, a, b, x)
    implicit none
    integer, intent(in) :: n ! Number unknowns
    real(8), dimension(n,3), intent(in)  :: a ! Tri-diagonal matrix
    real(8), dimension(n),   intent(in)  :: b ! Source
    real(8), dimension(n),   intent(out) :: x ! Solution
    !
    ! a(i,1) = west coefficients
    ! a(i,2) = central coefficients
    ! a(i,3) = east coefficients
    !
    ! Auxiliary variables
    integer :: i
    real(8), dimension(n) :: P
    real(8), dimension(n) :: Q

    i = 1

    P(i) = - a(i,3) / a(i,2)

    Q(i) = b(i) / a(i,2)

    do i = 2, n

       P(i) = - a(i,3) / ( a(i,2) + a(i,1) * P(i-1) )

       Q(i) = ( b(i) - a(i,1) * Q(i-1) ) / ( a(i,2) + a(i,1) * P(i-1) )

    end do

    i = n
    x(i) = Q(i)

    do i = n-1, 1, -1
       x(i) = x(i+1) * P(i) + Q(i)
    end do

   end subroutine tdma


end module thompson2d_hyperbolic
