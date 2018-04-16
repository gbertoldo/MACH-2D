module solvers

   use mtdma2d5
   use mtdma2d9
   use msi2d5
   use msi2d9

   implicit none

contains

   subroutine norm_l1_5d( nx, ny, var, b, a, norm)
      implicit none
      integer, intent(in) :: nx   ! Number of volumes in csi direction (real+fictitious)
      integer, intent(in) :: ny   ! Number of volumes in eta direction (real+fictitious)
      real(8), dimension (nx*ny),   intent(in) :: var ! Unknown variable
      real(8), dimension (nx*ny),   intent(in) :: b   ! Source vector of the linear system
      real(8), dimension (nx*ny,5), intent(in) :: a   ! Coefficients of the linear system
      real(8), intent(inout) :: norm

      ! Auxiliary variables
      integer :: i, j, np, nps, npn, npw, npe

      ! Norm is calculated taking into account only real volumes

      do j = 2, ny-1
         do i = 2, nx-1

            np   = nx * (j-1) + i
            nps  = np - nx
            npn  = np + nx
            npw  = np - 1
            npe  = np + 1

            norm = norm + abs(        &
               + a(np,1) * var(nps) &
               + a(np,2) * var(npw) &
               + a(np,3) * var(np ) &
               + a(np,4) * var(npe) &
               + a(np,5) * var(npn) &
               - b(np) )

         end do
      end do

   end subroutine norm_l1_5d


   subroutine norm_l1_9d( nx, ny, var, b, a, norm)
      implicit none
      integer, intent(in) :: nx   ! Number of volumes in csi direction (real+fictitious)
      integer, intent(in) :: ny   ! Number of volumes in eta direction (real+fictitious)
      real(8), dimension (nx*ny),   intent(in) :: var ! Unknown variable
      real(8), dimension (nx*ny),   intent(in) :: b   ! Source vector of the linear system
      real(8), dimension (nx*ny,9), intent(in) :: a   ! Coefficients of the linear system
      real(8), intent(inout) :: norm

      ! Auxiliary variables
      integer :: i, j, np, nps, npn, npw, npe, npsw, npse, npnw, npne

      ! Norm is calculated taking into account only real volumes

      do j = 2, ny-1
         do i = 2, nx-1

            np   = nx * (j-1) + i
            nps  = np - nx
            npn  = np + nx
            npw  = np - 1
            npe  = np + 1
            npsw = nps - 1
            npse = nps + 1
            npnw = npn - 1
            npne = npn + 1

            norm = norm + abs(        &
               + a(np,1) * var(npsw) &
               + a(np,2) * var(nps ) &
               + a(np,3) * var(npse) &
               + a(np,4) * var(npw ) &
               + a(np,5) * var(np  ) &
               + a(np,6) * var(npe ) &
               + a(np,7) * var(npnw) &
               + a(np,8) * var(npn ) &
               + a(np,9) * var(npne) &
               - b(np) )

         end do
      end do

   end subroutine norm_l1_9d


   !> \brief Calculates the norm L1 of a 9-diagonal linear system A . x = b
   !! If |b| > 10 * epsilon, the relative norm is calculated, otherwise the
   !! absolute norm is applied.
   subroutine norm_l1_9d_relative( nx, ny, var, b, a, norm) ! InOutput: last one
      implicit none
      integer, intent(in) :: nx   ! Number of volumes in csi direction (real+fictitious)
      integer, intent(in) :: ny   ! Number of volumes in eta direction (real+fictitious)
      real(8), dimension (nx*ny),   intent(in) :: var ! Unknown variable
      real(8), dimension (nx*ny),   intent(in) :: b   ! Source vector of the linear system
      real(8), dimension (nx*ny,9), intent(in) :: a   ! Coefficients of the linear system
      real(8), intent(inout) :: norm ! Norm L1

      ! Parameters

      real(8), parameter :: eps = 10.d0 * epsilon(1.d0)

      ! Auxiliary variables
      integer :: i, j, np, nps, npn, npw, npe, npsw, npse, npnw, npne

      real(8) :: norma ! Norm L1 of the linear system A . x = b
      real(8) :: normb ! Norm L1 of the source b

      ! Norm is calculated taking into account only real volumes

      norma = 0.d0

      normb = 0.d0

      do j = 2, ny-1

         do i = 2, nx-1

            np   = nx * (j-1) + i
            nps  = np - nx
            npn  = np + nx
            npw  = np - 1
            npe  = np + 1
            npsw = nps - 1
            npse = nps + 1
            npnw = npn - 1
            npne = npn + 1

            norma = norma + abs(        &
               + a(np,1) * var(npsw) &
               + a(np,2) * var(nps ) &
               + a(np,3) * var(npse) &
               + a(np,4) * var(npw ) &
               + a(np,5) * var(np  ) &
               + a(np,6) * var(npe ) &
               + a(np,7) * var(npnw) &
               + a(np,8) * var(npn ) &
               + a(np,9) * var(npne) &
               - b(np) )

            normb = normb + abs( b(np) )

         end do

      end do

      ! Summing the input norm with the calculated norm

      if ( normb > eps ) then ! Uses relative norm

         norm = norm + norma / normb

      else

         norm = norm + norma ! Uses absolute norm

      end if

   end subroutine norm_l1_9d_relative



   !> \brief Calculates the norm l1 of the source B of linear system Ax=B, taking into account only real volumes.
   real(8) function norm_l1_b( nx, ny, b)
      implicit none
      integer, intent(in) :: nx   !< Number of volumes in the csi direction (real+fictitious)
      integer, intent(in) :: ny   !< Number of volumes in the eta direction (real+fictitious)
      real(8), dimension (nx*ny),   intent(in) :: b   !< Source of the linear system

      ! Auxiliary variables
      integer :: i, j, np
      real(8) :: norm

      ! Norm is calculated taking into account only real volumes

      norm = 0.d0

      do j = 2, ny-1
         do i = 2, nx-1

            np   = nx * (j-1) + i

            norm = norm + abs( b(np) )

         end do
      end do

      norm_l1_b = norm

   end function


   !> \brief Calculates the coefficient for convergence criteria of a 9-diagonal matrix.
   !! This coefficient is calculated only for real volumes.
   real(8) function get_maxc_9d(nx, ny, a)
      implicit none
      integer, intent(in) :: nx   !< Number of volumes in the csi direction (real+fictitious)
      integer, intent(in) :: ny   !< Number of volumes in the eta direction (real+fictitious)
      real(8), dimension (nx*ny,9), intent(in) :: a   !< Coefficients of the linear system

      ! Inner variables

      integer :: i, j, np

      real(8) :: c, maxc

      ! Looking for the maxval of c = ( sum_nb | A_nb | ) / | A_P |

      maxc = 0.d0

      do i = 2, nx-1

         do j = 2, ny-1

            np   = nx * (j-1) + i

            c = ( sum(abs(a(np,1:4))) + sum(abs(a(np,6:9))) ) / abs( a(np,5) )

            if ( maxc < c ) maxc = c

         end do

      end do

      get_maxc_9d = maxc

   end function


   !> \brief Calculates the coefficient for convergence criteria of a 5-diagonal matrix.
   !! This coefficient is calculated only for real volumes.
   real(8) function get_maxc_5d(nx, ny, a)
      implicit none
      integer, intent(in) :: nx   !< Number of volumes in the csi direction (real+fictitious)
      integer, intent(in) :: ny   !< Number of volumes in the eta direction (real+fictitious)
      real(8), dimension (nx*ny,5), intent(in) :: a   !< Coefficients of the linear system

      ! Inner variables

      integer :: i, j, np

      real(8) :: c, maxc

      ! Looking for the maxval of c = ( sum_nb | A_nb | ) / | A_P |

      maxc = 0.d0

      do i = 2, nx-1

         do j = 2, ny-1

            np   = nx * (j-1) + i

            c = ( sum(abs(a(np,1:2))) + sum(abs(a(np,4:5))) ) / abs( a(np,3) )

            if ( maxc < c ) maxc = c

         end do

      end do

      get_maxc_5d = maxc

   end function


   !> \brief Calculates for each real volume the convergence coefficient
   !! associated with the 9-diagonal matrix 'a' and stores the results in cc.
   subroutine get_cc_9d(nx, ny, a, cc)
      implicit none
      integer, intent(in) :: nx   !< Number of volumes in the csi direction (real+fictitious)
      integer, intent(in) :: ny   !< Number of volumes in the eta direction (real+fictitious)
      real(8), dimension (nx*ny,9), intent(in)  :: a   !< Coefficients of the linear system
      real(8), dimension (nx*ny),   intent(out) :: cc  !< Convergence coef. vector

      ! Inner variables

      integer :: i, j, np

      do i = 2, nx-1

         do j = 2, ny-1

            np = nx * (j-1) + i

            cc(np) = ( sum( abs(a(np,1:4)) ) + sum( abs(a(np,6:9)) ) ) &

            / abs(a(np,5))

         end do

      end do

   end subroutine

   !> \brief Calculates for each real volume the convergence coefficient
   !! associated with the 5-diagonal matrix 'a' and stores the results in cc.
   subroutine get_cc_5d(nx, ny, a, cc)
      implicit none
      integer, intent(in) :: nx   !< Number of volumes in the csi direction (real+fictitious)
      integer, intent(in) :: ny   !< Number of volumes in the eta direction (real+fictitious)
      real(8), dimension (nx*ny,5), intent(in)  :: a   !< Coefficients of the linear system
      real(8), dimension (nx*ny),   intent(out) :: cc  !< Convergence coef. vector

      ! Inner variables

      integer :: i, j, np

      do i = 2, nx-1

         do j = 2, ny-1

            np = nx * (j-1) + i

            cc(np) = ( sum( abs(a(np,1:2)) ) + sum( abs(a(np,4:5)) ) ) &

            / abs(a(np,3))

         end do

      end do

   end subroutine

   !> \brief Calculates the mean value of the convergence coefficient
   !! associated with the 9-diagonal matrix 'a'.
   real(8) function get_mean_cc_9d(nx, ny, a)
      implicit none
      integer, intent(in) :: nx   !< Number of volumes in the csi direction (real+fictitious)
      integer, intent(in) :: ny   !< Number of volumes in the eta direction (real+fictitious)
      real(8), dimension (nx*ny,9), intent(in)  :: a   !< Coefficients of the linear system

      ! Inner variables

      integer :: i, j, np

      get_mean_cc_9d = 0.d0

      do i = 2, nx-1

         do j = 2, ny-1

            np = nx * (j-1) + i

            get_mean_cc_9d = ( sum( abs(a(np,1:4)) ) + sum( abs(a(np,6:9)) ) ) &

            / abs(a(np,5)) + get_mean_cc_9d

         end do

      end do

      get_mean_cc_9d = get_mean_cc_9d / dble( (nx-2) * (ny-2) )

   end function

   !> \brief Calculates the mean value of convergence coefficient
   !! associated with the 5-diagonal matrix 'a'.
   real(8) function get_mean_cc_5d(nx, ny, a)
      implicit none
      integer, intent(in) :: nx   !< Number of volumes in the csi direction (real+fictitious)
      integer, intent(in) :: ny   !< Number of volumes in the eta direction (real+fictitious)
      real(8), dimension (nx*ny,5), intent(in)  :: a   !< Coefficients of the linear system

      ! Inner variables

      integer :: i, j, np

      get_mean_cc_5d = 0.d0

      do i = 2, nx-1

         do j = 2, ny-1

            np = nx * (j-1) + i

            get_mean_cc_5d = ( sum( abs(a(np,1:2)) ) + sum( abs(a(np,4:5)) ) ) &

            / abs(a(np,3)) + get_mean_cc_5d

         end do

      end do

      get_mean_cc_5d = get_mean_cc_5d / dble( (nx-2) * (ny-2) )

   end function

end module solvers
