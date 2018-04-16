!> \brief Defines the geometry using splines such that
!! Input: - function first derivatives are given for the whole domain
!! Besides the first derivative that must be given for the domain, there are the following boundary condition options (opt):
!! 0 = Given second (diff2YN) and third (diff3YN) derivatives on the east boundary
!! 1 = Given second derivative on the east boundary (third derivative is calculated auto. considering the fourth derivative equals zero)
!! 2 = Third and fourth derivatives are equal zero on the east boundary (second derivative is calculated automatically)
module geometry_spline_s01

   implicit none

   private

   public get_spline_s01_geometry

contains


   !> \brief Gets the geometry spline s01
   subroutine get_spline_s01_geometry(nr, nsi, aks, rb, lr, xsi, dy1, dy2, dy3, opt, xr, yr, es) ! Output: last three
      implicit none
      integer, intent(in)  :: nr         !< Number of partitions of the output curve (xr,yr)
      integer, intent(in)  :: nsi        !< Number of inputs for the spline curve (xsi,ysi)
      real(8), intent(in)  :: aks        !< Clustering parameter
      real(8), intent(in)  :: rb         !< Base radius
      real(8), intent(in)  :: lr         !< Body length
      real(8), intent(in)  :: xsi(0:nsi) !< x/lr coordinate for the spline interpolation
      real(8), intent(in)  :: dy1(0:nsi) !< First derivatives
      real(8), intent(in)  :: dy2        !< Second derivative on the east boundary
      real(8), intent(in)  :: dy3        !< Third derivative on the east boundary
      integer, intent(in)  :: opt        !< See bellow
      real(8), intent(out) :: xr(0:nr)   !< x-coordinate
      real(8), intent(out) :: yr(0:nr)   !< y-coordinate
      integer, intent(out) :: es         !< Exit status: 0=success, 1=non-monotonic, 2=negative, 3=1+2
      ! Besides the first derivative that must be given for the domain, there are the following boundary condition options (opt):
      ! 0 = Given second (diff2YN) and third (diff3YN) derivatives on the east boundary
      ! 1 = Given second derivative on the east boundary (third derivative is calculated auto. considering the fourth derivative equals zero)
      ! 2 = Third and fourth derivatives are equal zero on the east boundary (second derivative is calculated automatically)

      ! Inner variables

      integer :: i, j      ! Dummy index
      real(8) :: rf        ! Frontal radius
      real(8) :: xs(0:nsi) ! x coordinate for the spline interpolation
      real(8) ::  a(0:nsi) ! Spline coefficients
      real(8) ::  b(0:nsi) ! Spline coefficients
      real(8) ::  c(0:nsi) ! Spline coefficients
      real(8) ::  d(0:nsi) ! Spline coefficients
      real(8) :: SI(0:nsi-1) ! Sum of integrals

      ! Initializing exit status

      es = 0

      ! Completing the informations for the spline interpolation

      xs = xsi * lr

      ! Generating the desired distribution of points over xr

      do i = 1, nr

         xr(i) = lr * ( dble(i-1) / dble(nr-1) ) ** aks

      end do

      ! Calculating the spline coefficients

      call get_spline_coefficients(nsi, xs, dy1, dy2, dy3, opt, a, b, c, d, SI, es) ! Output: last four

      ! Checking for positive shapes

      rf = rb - SI(nsi-1)

      if ( rf < 0.d0  ) es = es + 2


      ! Interpolating

      j = 0

      do i = 1, nr

         do

            if ( xr(i) <= xs(j+1) ) then

               if (j == 0 ) then

                  yr(i) = rf                                 &
                     + a(j) * ( xr(i)-xs(j) )                &
                     + b(j) * ( xr(i)-xs(j) ) ** 2.d0 / 2.d0 &
                     + c(j) * ( xr(i)-xs(j) ) ** 3.d0 / 3.d0 &
                     + d(j) * ( xr(i)-xs(j) ) ** 4.d0 / 4.d0

               else

                  yr(i) = rf + SI(j-1)                       &
                     + a(j) * ( xr(i)-xs(j) )                &
                     + b(j) * ( xr(i)-xs(j) ) ** 2.d0 / 2.d0 &
                     + c(j) * ( xr(i)-xs(j) ) ** 3.d0 / 3.d0 &
                     + d(j) * ( xr(i)-xs(j) ) ** 4.d0 / 4.d0

               end if

               exit

            else

               j = j + 1

            end if

         end do

      end do


      xr(0) = 0.d0
      yr(0) = 0.d0

      xr(nr) = lr
      yr(nr) = rb

   end subroutine get_spline_s01_geometry


   !> \brief Calculates the spline coefficients
   subroutine get_spline_coefficients(n, x, dy1, dy2, dy3, opt, a, b, c, d, SI, es) ! Output: last six
      implicit none
      integer, intent(in)  :: n        !< Number of partitions
      real(8), intent(in)  :: x(0:n)   !< x coordinates
      real(8), intent(in)  :: dy1(0:n) !< First derivatives
      real(8), intent(in)  :: dy2      !< Second derivative on the east boundary
      real(8), intent(in)  :: dy3      !< Third derivative on the east boundary
      integer, intent(in)  :: opt      !< See bellow
      real(8), intent(out) :: a(0:n)   !< Spline coefficients
      real(8), intent(out) :: b(0:n)   !< Spline coefficients
      real(8), intent(out) :: c(0:n)   !< Spline coefficients
      real(8), intent(out) :: d(0:n)   !< Spline coefficients
      real(8), intent(out) :: SI(0:n-1)!< Sum of integrals
      integer, intent(inout) :: es       !< Exit status: 0=success, 1=failure
      ! Besides the first derivative that must be given for the domain, there are the following boundary condition options (opt):
      ! 0 = Given second (diff2YN) and third (diff3YN) derivatives on the east boundary
      ! 1 = Given second derivative on the east boundary (third derivative is calculated auto. considering the fourth derivative equals zero)
      ! 2 = Third and fourth derivatives are equal zero on the east boundary (second derivative is calculated automatically)


      ! Inner variables
      integer :: i
      real(8) :: h(0:n-1)

      ! Parameters
      real(8), parameter :: eps = 1.d-13


      ! Calculating h

      do i = 0, n-1

         h(i) = x(i+1)-x(i)

      end do

      ! Initial conditions

      a = dy1

      select case (opt)

         case (0) ! Given second and third derivatives on the east boundary

            b(n) = dy2

            c(n) = dy3 / 2.d0

            d(n) = 0.d0

         case (1) ! Given second derivative on the east boundary (third derivative is calculated considering the fourth derivative equals zero)

            b(n) = dy2

            c(n) = b(n) / h(n-1) + (a(n-1)-a(n)) / h(n-1)**2

            d(n) = 0.d0

         case (2) ! Third and fourth derivatives are equal zero on the east boundary (second derivative is calculated)

            b(n) = (a(n)-a(n-1)) / h(n-1)

            c(n) = 0.d0

            d(n) = 0.d0

         case default

            write(*,*) 'Unknown option for boundary condition. Stopping.'

            stop

      end select


      do i = n-1, 0, -1

         c(i) = 3.d0 * (b(i+1)/h(i)+(a(i)-a(i+1))/h(i)**2 ) - 2.d0 * c(i+1)

         b(i) = b(i+1) - h(i) * (c(i+1)+c(i))

         d(i) = ( c(i+1) - c(i) ) / ( 3.d0*h(i) )

      end do

      ! Calculating the sum of integrals

      i = 0

      SI(i) = a(i) * h(i)             &
         +    b(i) * h(i) ** 2 / 2.d0 &
         +    c(i) * h(i) ** 3 / 3.d0 &
         +    d(i) * h(i) ** 4 / 4.d0

      do i = 1, n-1

         SI(i) = a(i) * h(i)             &
            +    b(i) * h(i) ** 2 / 2.d0 &
            +    c(i) * h(i) ** 3 / 3.d0 &
            +    d(i) * h(i) ** 4 / 4.d0 &
            +    SI(i-1)

      end do

      ! Checking monotonicity ( approximately )

      do i = 0, n

         if ( a(i) < -eps .or. eps < b(i) ) then

            es = 1

            return

         end if

      end do

   end subroutine get_spline_coefficients

end module geometry_spline_s01
