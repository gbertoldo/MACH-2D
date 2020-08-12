!>
!! \brief Contains procedures to operate fields
!!
module mod_field_operations
   implicit none

contains

   !> \brief Calculates the field F over the east boundary
   subroutine get_east_boundary_field(nx, ny, F, Fbe)
      implicit none
      integer, intent(in) :: nx     !< Number of volumes in the csi direction (real+fictitious)
      integer, intent(in) :: ny     !< Number of volumes in the eta direction (real+fictitious)
      real(8), dimension(nx*ny), intent(in)  :: F   !< A generic field
      real(8), dimension(ny),    intent(out) :: Fbe !< The field over the east boundary

      ! Inner variables

      integer :: i, j, np, npe ! Dummy variables

      i = nx - 1

      do j = 2, ny-1

         np   = nx * (j-1) + i
         npe  = np + 1

         Fbe(j) = 0.5d0 * ( F(np) + F(npe) )

      end do

   end subroutine get_east_boundary_field


   !> \brief Calculates the field F over the west boundary
   subroutine get_west_boundary_field(nx, ny, F, Fbw)
      implicit none
      integer, intent(in) :: nx     !< Number of volumes in the csi direction (real+fictitious)
      integer, intent(in) :: ny     !< Number of volumes in the eta direction (real+fictitious)
      real(8), dimension(nx*ny), intent(in)  :: F   !< A generic field
      real(8), dimension(ny),    intent(out) :: Fbw !< The field over the west boundary

      ! Inner variables

      integer :: i, j, np, npw ! Dummy variables

      i = 2

      do j = 2, ny-1

         np   = nx * (j-1) + i
         npw  = np - 1

         Fbw(j) = 0.5d0 * ( F(np) + F(npw) )

      end do

   end subroutine get_west_boundary_field


   !> \brief Calculates the field F over the north boundary
   subroutine get_north_boundary_field(nx, ny, F, Fbn)
      implicit none
      integer, intent(in) :: nx     !< Number of volumes in the csi direction (real+fictitious)
      integer, intent(in) :: ny     !< Number of volumes in the eta direction (real+fictitious)
      real(8), dimension(nx*ny), intent(in)  :: F   !< A generic field
      real(8), dimension(nx),    intent(out) :: Fbn !< The field over the north boundary

      ! Inner variables

      integer :: i, j, np, npn ! Dummy variables

      j = ny-1

      do i = 2, nx-1

         np   = nx * (j-1) + i
         npn  = np + nx

         Fbn(i) = 0.5d0 * ( F(np) + F(npn) )

      end do

   end subroutine get_north_boundary_field


   !> \brief Calculates the field F over the south boundary
   subroutine get_south_boundary_field(nx, ny, F, Fbs)
      implicit none
      integer, intent(in) :: nx     !< Number of volumes in the csi direction (real+fictitious)
      integer, intent(in) :: ny     !< Number of volumes in the eta direction (real+fictitious)
      real(8), dimension(nx*ny), intent(in)  :: F   !< A generic field
      real(8), dimension(nx),    intent(out) :: Fbs !< The field over the south boundary

      ! Inner variables

      integer :: i, j, np, nps ! Dummy variables

      j = 2

      do i = 2, nx-1

         np   = nx * (j-1) + i
         nps  = np - nx

         Fbs(i) = 0.5d0 * ( F(np) + F(nps) )

      end do

   end subroutine get_south_boundary_field



   !> \brief Calculates the average value of a given field
   subroutine get_average(nx, ny, F, F_avg) ! Output: last one
      implicit none
      integer, intent(in)  :: nx !< Number of volumes in csi direction (real+fictitious)
      integer, intent(in)  :: ny !< Number of volumes in eta direction (real+fictitious)
      real(8), dimension(nx*ny), intent(in)  :: F !< Field for which the mean will be calculated
      real(8), intent(out) :: F_avg !< Average of F

      ! Inner variables

      integer :: i, j, np

      F_avg = 0.d0

      do j = 2, ny-1

         do i = 2, nx-1

            np   = nx * (j-1) + i

            F_avg = F_avg + F(np)

         end do

      end do

      F_avg = F_avg / ( dble(nx-2) * dble(ny-2) )

   end subroutine get_average


end module
