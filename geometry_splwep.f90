!> \brief Defines the Shifted Power Law With Exponential Perturbation geometry
module geometry_splwep

   implicit none

   private

   public get_splwep

contains


   !> \brief Gets the Shifted Power Law With Exponential Perturbation geometry
   subroutine get_splwep(nr, aks, rb, lr, h, lbd, S0, a, xr, yr) ! Output: last two
      implicit none
      integer, intent(in)  :: nr       !< Number of partitions of the output curve (xr,yr)
      real(8), intent(in)  :: aks      !< Clustering parameter
      real(8), intent(in)  :: rb       !< Base radius
      real(8), intent(in)  :: lr       !< Body length
      real(8), intent(in)  :: h        !< Radius ratio
      real(8), intent(in)  :: lbd      !< Exponent of the power law
      real(8), intent(in)  :: S0       !< Maximum perturbation
      real(8), intent(in)  :: a        !< Coefficient of the exponential
      real(8), intent(out) :: xr(0:nr) !< x-coordinate
      real(8), intent(out) :: yr(0:nr) !< y-coordinate

      ! Inner variables

      integer :: i  ! Dummy index
      real(8) :: x  ! x coord.
      real(8) :: xc ! x coord. of the critical point
      real(8) :: S  ! Perturbation function


      ! Calculating xc (critical point of the perturbation)

      xc = ( a + 2.d0 - sqrt(4.d0+a*a) ) / ( 2.d0 * a )


      ! Defining the base geometry

      xr(0) = 0.d0
      yr(0) = 0.d0

      do i = 1, nr

         xr(i) = lr * ( dble(i-1) / dble(nr-1) ) ** aks

         x = xr(i) / lr

         S = S0 * x * ( 1.d0 - x ) * exp( -a*(x-xc) ) / ( xc * (1.d0-xc) )

         yr(i) = rb * ( h ** (1.d0/lbd) + (1.d0-h ** (1.d0/lbd)) * x ) ** lbd &

            & * (1.d0+S)

      end do

   end subroutine get_splwep

end module geometry_splwep
