!> \brief Defines the geometry using Natural Cubic Spline interpolation
module geometry_ncspline

   use spline

   implicit none

   private

   public get_ncspline_geometry

contains


   !> \brief Gets a geometry defined by a natural cubic spline
   subroutine get_ncspline_geometry(nr, nsi, aks, rb, lr, xsi, ysi, xr, yr, es) ! Output: last three
      implicit none
      integer, intent(in)  :: nr       !< Number of partitions of the output curve (xr,yr)
      integer, intent(in)  :: nsi      !< Number of inputs for the spline curve (xsi,ysi)
      real(8), intent(in)  :: aks      !< Clustering parameter
      real(8), intent(in)  :: rb       !< Base radius
      real(8), intent(in)  :: lr       !< Body length
      real(8), intent(in)  :: xsi(nsi) !< x/lr coordinate for the spline interpolation
      real(8), intent(in)  :: ysi(nsi) !< y/rb-coordinate perturbation for the spline interpolation
      real(8), intent(out) :: xr(0:nr) !< x-coordinate
      real(8), intent(out) :: yr(0:nr) !< y-coordinate
      integer, intent(out) :: es       !< Exit status: 0=success, 1=failure

      ! Inner variables

      integer :: i           ! Dummy index
      real(8) :: xs(nsi+1) ! x/lr coordinate for the spline interpolation
      real(8) :: ys(nsi+1) ! y-coordinate perturbation for the spline interpolation


      ! Initializing exit status

      es = 0

      ! Completing the informations for the spline interpolation

      do i = 1, nsi

         xs(i) = xsi(i) * lr
         ys(i) = ysi(i) * rb

      end do

      xs(nsi+1) = lr
      ys(nsi+1) = rb

      ! Generating the desired distribution of points over xr

      do i = 1, nr

         xr(i) = lr * ( dble(i-1) / dble(nr-1) ) ** aks

      end do

      ! Calculating the natural cubic splines
      call get_cspline(nsi, nr-1, xs, ys, xr(1:nr), yr(1:nr)) ! Output: last one

      xr(0) = 0.d0
      yr(0) = 0.d0

      xr(nr) = lr
      yr(nr) = rb

      ! Checking if the geometry increases monotonically

      call get_monotonicity(nsi, xs, ys, es) ! Output: last one

   end subroutine get_ncspline_geometry


   !> \brief Checks the monotonicity of the splines (f'>=0 and f''<=0)
   subroutine get_monotonicity(n, x, y, es) ! Output: last one
      implicit none
      integer, intent(in)  :: n      !< Number of intervals
      real(8), intent(in)  :: x(0:n) !< x-coordinates of the spline
      real(8), intent(in)  :: y(0:n) !< y-coordinates of the spline
      integer, intent(out) :: es     !< Exit status: 0=monotonic, 1=non-monotonic

      ! Inner variables

      integer :: i      ! Dummy index
      real(8) :: b(0:n) ! b coefficients
      real(8) :: c(0:n) ! c coefficients
      real(8) :: d(0:n) ! d coefficients

      ! Parameters

      real(8), parameter :: eps = 1.d-13

      ! Initializing es

      es = 0


      ! Calculates the natural cubic splines coefficients
      call get_csplines_coeff(n, x, y, b, c, d) ! Output: last three


      ! Checking monotonicity

      do i = 0, n

         if ( b(i) < -eps .or. eps < c(i) ) then

            es = 1

            return

         end if

      end do

   end subroutine get_monotonicity


end module geometry_ncspline
