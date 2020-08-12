!> \brief Defines the geometry using cubic splines with second derivatives prescribed
!! on all points and the function values given on the boundaries
module geometry_d2cspline

   use spline

   implicit none

   private

   public get_d2cspline_geometry

contains


   !> \brief Gets a geometry using cubic splines with second derivatives prescribed
   !! on all points and the function values given on the boundaries
   subroutine get_d2cspline_geometry(nr, nsi, aks, rb, lr, xsi, y0, yn, dy2, xr, yr, es) ! Output: last three
      implicit none
      integer, intent(in)  :: nr       !< Number of partitions of the output curve (xr,yr)
      integer, intent(in)  :: nsi      !< Number of inputs for the spline curve (xsi,dy2)
      real(8), intent(in)  :: aks      !< Clustering parameter
      real(8), intent(in)  :: rb       !< Base radius
      real(8), intent(in)  :: lr       !< Body length
      real(8), intent(in)  :: xsi(nsi) !< x coordinate for the spline interpolation
      real(8), intent(in)  :: y0       !< f(x0)
      real(8), intent(in)  :: yn       !< f(xn)
      real(8), intent(in)  :: dy2(nsi) !< diff( y, x, 2)
      real(8), intent(out) :: xr(0:nr) !< x-coordinate
      real(8), intent(out) :: yr(0:nr) !< y-coordinate
      integer, intent(out) :: es       !< Exit status: 0=success, 1=failure

      ! Inner variables

      integer :: i       ! Dummy index
      real(8) :: xs(nsi) ! x coordinate for the spline interpolation
      real(8) :: ys0     ! y-coordinate for the spline interpolation
      real(8) :: ysn     ! y-coordinate for the spline interpolation


      ! Initializing exit status

      es = 0

      ! Completing the informations for the spline interpolation

      xs = xsi * lr
      ys0 = y0 * rb
      ysn = yn * rb

      ! Generating the desired distribution of points over xr

      do i = 1, nr

         xr(i) = lr * ( dble(i-1) / dble(nr-1) ) ** aks

      end do

      ! Calculates the cubic splines with second derivatives prescribed
      ! on all points and the function values given on the boundaries
      call get_d2cspline(nsi-1, nr-1, xs, ys0, ysn, dy2, xr(1:nr), yr(1:nr)) ! Output: last one


      xr(0) = 0.d0
      yr(0) = 0.d0

      xr(nr) = lr
      yr(nr) = rb

      ! Checking if the geometry increases monotonically

      call get_monotonicity(nsi-1, xs, ys0, ysn, dy2, es) ! Output: last one

   end subroutine get_d2cspline_geometry


   !> \brief Checks the monotonicity of the splines (f'>=0 and f''<=0)
   subroutine get_monotonicity(n, x, y0, yn, dy2, es) ! Output: last one
      implicit none
      integer, intent(in)  :: n        !< Number of intervals
      real(8), intent(in)  :: x(0:n)   !< x-coordinates of the spline
      real(8), intent(in)  :: y0       !< f(x0)
      real(8), intent(in)  :: yn       !< f(xn)
      real(8), intent(in)  :: dy2(0:n) !< diff( y, x, 2)
      integer, intent(out) :: es       !< Exit status: 0=monotonic, 1=non-monotonic

      ! Inner variables

      integer :: i      ! Dummy index
      real(8) :: a(0:n) ! a coefficients
      real(8) :: b(0:n) ! b coefficients
      real(8) :: c(0:n) ! c coefficients
      real(8) :: d(0:n) ! d coefficients

      ! Parameters

      real(8), parameter :: eps = 1.d-13

      ! Initializing es

      es = 0


      ! Calculates the cubic splines coefficients with second derivatives prescribed
      ! and function prescribed on boundaries
      call get_d2csplines_coeff(n, x, y0, yn, dy2, a, b, c, d) ! Output: last four


      ! Checking monotonicity

      do i = 0, n

         if ( b(i) < -eps .or. eps < c(i) ) then

            es = 1

            return

         end if

      end do

   end subroutine get_monotonicity


end module geometry_d2cspline
