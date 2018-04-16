!> \brief Defines the Shifted Power Law With Spline Perturbation geometry
module geometry_splwsp

   use spline

   implicit none

   private

   public get_splwsp

contains


   !> \brief Gets the Shifted Power Law With Spline Perturbation geometry
   subroutine get_splwsp(nr, nsi, aks, rb, lr, h, lbd, xsi, ysi, xr, yr, es) ! Output: last three
      implicit none
      integer, intent(in)  :: nr       !< Number of partitions of the output curve (xr,yr)
      integer, intent(in)  :: nsi      !< Number of inputs for the spline curve (xsi,ysi)
      real(8), intent(in)  :: aks      !< Clustering parameter
      real(8), intent(in)  :: rb       !< Base radius
      real(8), intent(in)  :: lr       !< Body length
      real(8), intent(in)  :: h        !< Radius ratio
      real(8), intent(in)  :: lbd      !< Exponent of the power law
      real(8), intent(in)  :: xsi(nsi) !< x/lr coordinate for the spline interpolation
      real(8), intent(in)  :: ysi(nsi) !< y-coordinate perturbation for the spline interpolation
      real(8), intent(out) :: xr(0:nr) !< x-coordinate
      real(8), intent(out) :: yr(0:nr) !< y-coordinate
      integer, intent(out) :: es       !< Exit status: 0=success, 1=failure

      ! Inner variables

      integer :: i           ! Dummy index
      real(8) :: xs(0:nsi+1) ! x/lr coordinate for the spline interpolation
      real(8) :: ys(0:nsi+1) ! y-coordinate perturbation for the spline interpolation
      real(8) :: xrs(1:nr)   ! x-coordinate
      real(8) :: yrs(1:nr)   ! y-coordinate perturbation for the spline interpolation for xrs


      ! Initializing exit status

      es = 0

      ! Defining the base geometry

      xr(0) = 0.d0
      yr(0) = 0.d0

      do i = 1, nr

         xr(i) = lr * ( dble(i-1) / dble(nr-1) ) ** aks

         yr(i) = rb * ( h ** (1.d0/lbd) + (1.d0-h ** (1.d0/lbd)) * xr(i) / lr ) ** lbd

      end do


      ! Completing the informations for the spline interpolation

      xs(0) = 0.d0
      ys(0) = 0.d0

      do i = 1, nsi

         xs(i) = xsi(i) * lr
         ys(i) = ysi(i)

      end do

      xs(nsi+1) = lr
      ys(nsi+1) = 0.d0

      xrs = xr(1:nr)

      ! Calculating the natural cubic splines
      call get_cspline(nsi+1, nr-1, xs, ys, xrs, yrs) ! Output: last one

      ! Calculating the perturbed geometry
      yr(1:nr) = yr(1:nr) * ( 1.d0 + yrs )

      ! Checking if the geometry increases monotonically

      call get_monotonicity(nr, xr, yr, es) ! Output: last one

   end subroutine get_splwsp


   !> \brief Verifies the monotonocity of the geometry profile
   subroutine get_monotonicity(nr, xr, yr, es) ! Output: last one
      implicit none
      integer, intent(in)  :: nr       !< Number of partitions of the output curve (xr,yr)
      real(8), intent(in)  :: xr(0:nr) !< x-coordinate
      real(8), intent(in)  :: yr(0:nr) !< y-coordinate
      integer, intent(out) :: es       !< Exit status: 0=monotonic, 1=non-monotonic

      ! Inner variables

      integer :: i  ! Dummy index
      real(8) :: y2 ! diff(y,x,2)

      ! Initializing exit status
      es = 0

      ! Checking monotonocity

      do i = 2, nr-1

         y2 = 2.d0 * ( &

            + yr(i+1) / ( (xr(i+1)-xr(i-1)) * (xr(i+1)-xr(i  )) ) &

            + yr(i  ) / ( (xr(i  )-xr(i-1)) * (xr(i  )-xr(i+1)) ) &

            + yr(i-1) / ( (xr(i-1)-xr(i  )) * (xr(i-1)-xr(i+1)) ) &

            )

         if ( y2 > 1.d-13 ) then

            es = 1

            return

         end if

      end do

   end subroutine get_monotonicity



end module geometry_splwsp
