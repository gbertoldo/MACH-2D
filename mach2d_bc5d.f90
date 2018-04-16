!> \brief Defines the subroutines related to the calculation of the boundary
!! conditions. More specifically, this module deals with numerical schemes that
!! envolves each fictitious volume and its neighbours (up to four).

module bc5d
   implicit none

contains


   !> \brief Defines the numerical scheme for the Dirichlet boundary condition
   !! on the north boundary
   subroutine get_bc_dirichlet_north_5d(nx, Fbn, a5bn, b5bn)
      implicit none
      integer, intent(in) :: nx !< Number of volumes in csi direction (real+fictitious)
      real(8), dimension(nx),   intent(in)  :: Fbn  !< Field F over the north boundary
      real(8), dimension(nx,5), intent(out) :: a5bn !< Coefficients of the discretization for the north boundary
      real(8), dimension(nx),   intent(out) :: b5bn !< Source of the discretization for the north boundary

      ! Inner variables

      integer :: i

      ! North volumes

      do i = 2, nx-1

         a5bn(i,1) = 1.d0 ! S
         a5bn(i,2) = 0.d0 ! W
         a5bn(i,3) = 1.d0 ! P
         a5bn(i,4) = 0.d0 ! E
         a5bn(i,5) = 0.d0 ! N

         b5bn(i) = 2.d0 * Fbn(i)

      end do

   end subroutine get_bc_dirichlet_north_5d



   !> \brief Defines the numerical scheme for the null normal gradient over the
   !! south boundary.
   subroutine get_ibc_null_normal_grad_south_5d(nx, a5bs, b5bs) ! Output: last two
      implicit none
      integer, intent(in) :: nx !< Number of volumes in csi direction (real+fictitious)
      real(8), dimension(nx,5), intent(out) :: a5bs !< Coefficients of the discretization for the south boundary
      real(8), dimension(nx),   intent(out) :: b5bs !< Source of the discretization for the south boundary

      ! Inner variables

      integer :: i

      ! South volumes

      do i = 2, nx-1

         a5bs(i,1) =  0.d0 ! S
         a5bs(i,2) =  0.d0 ! W
         a5bs(i,3) =  1.d0 ! P
         a5bs(i,4) =  0.d0 ! E
         a5bs(i,5) = -1.d0 ! N

         b5bs(i) = 0.d0

      end do

   end subroutine get_ibc_null_normal_grad_south_5d



   !> \brief Defines the numerical scheme for the null normal gradient over the
   !! west boundary.
   subroutine get_ibc_null_normal_grad_west_5d(ny, a5bw, b5bw) ! Output: last two
      implicit none
      integer, intent(in) :: ny !< Number of volumes in eta direction (real+fictitious)
      real(8), dimension(ny,5), intent(out) :: a5bw !< Coefficients of the discretization for the west boundary
      real(8), dimension(ny),   intent(out) :: b5bw !< Source of the discretization for the west boundary

      ! Inner variables

      integer :: j

      ! West volumes

      do j = 2, ny-1

         a5bw(j,1) =  0.d0 ! S
         a5bw(j,2) =  0.d0 ! W
         a5bw(j,3) =  1.d0 ! P
         a5bw(j,4) = -1.d0 ! E
         a5bw(j,5) =  0.d0 ! N

         b5bw(j) = 0.d0

      end do

   end subroutine get_ibc_null_normal_grad_west_5d


   !> \brief Defines the numerical scheme for the null normal gradient over the
   !! east boundary.
   subroutine get_ibc_null_normal_grad_east_5d(ny, a5be, b5be) ! Output: last two
      implicit none
      integer, intent(in) :: ny !< Number of volumes in eta direction (real+fictitious)
      real(8), dimension(ny,5), intent(out) :: a5be !< Coefficients of the discretization for the east boundary
      real(8), dimension(ny),   intent(out) :: b5be !< Source of the discretization for the east boundary

      ! Inner variables

      integer :: j

      ! East volumes

      do j = 2, ny-1

         a5be(j,1) =  0.d0 ! S
         a5be(j,2) = -1.d0 ! W
         a5be(j,3) =  1.d0 ! P
         a5be(j,4) =  0.d0 ! E
         a5be(j,5) =  0.d0 ! N

         b5be(j) = 0.d0

      end do

   end subroutine get_ibc_null_normal_grad_east_5d



   !> \brief Defines the numerical scheme for the calculation of the boundary
   !! conditions at the fictitious volumes at corners of the transformed domain
   subroutine get_bc_corners_5d(nx, ny, a5bn, a5bs, a5be, a5bw, b5bn, b5bs &
         , b5be, b5bw ) ! Output: last eight
      implicit none
      integer, intent(in) :: nx !< Number of volumes in csi direction (real+fictitious)
      integer, intent(in) :: ny !< Number of volumes in eta direction (real+fictitious)
      real(8), dimension(nx,5), intent(out) :: a5bn !< Coefficients of the discretization for the north boundary
      real(8), dimension(nx,5), intent(out) :: a5bs !< Coefficients of the discretization for the south boundary
      real(8), dimension(ny,5), intent(out) :: a5be !< Coefficients of the discretization for the east boundary
      real(8), dimension(ny,5), intent(out) :: a5bw !< Coefficients of the discretization for the west boundary
      real(8), dimension(nx),   intent(out) :: b5bn !< Source of the discretization for the north boundary
      real(8), dimension(nx),   intent(out) :: b5bs !< Source of the discretization for the south boundary
      real(8), dimension(ny),   intent(out) :: b5be !< Source of the discretization for the east boundary
      real(8), dimension(ny),   intent(out) :: b5bw !< Source of the discretization for the west boundary

      ! Inner variables

      integer :: i, j

      ! SW corner

      i = 1

      a5bs(i,:) = 0.d0
      a5bs(i,3) = 1.d0
      b5bs(i)   = 0.d0

      j = 1

      a5bw(j,:) = a5bs(i,:)
      b5bw(j)   = b5bs(i)


      ! SE corner

      i = nx

      a5bs(i,:) = 0.d0
      a5bs(i,3) = 1.d0
      b5bs(i)   = 0.d0

      j = 1

      a5be(j,:) = a5bs(i,:)
      b5be(j)   = b5bs(i)


      ! NW corner

      i = 1

      a5bn(i,:) = 0.d0
      a5bn(i,3) = 1.d0
      b5bn(i)   = 0.d0

      j = ny

      a5bw(j,:) = a5bn(i,:)
      b5bw(j)   = b5bn(i)


      ! NE corner

      i = nx

      a5bn(i,:) = 0.d0
      a5bn(i,3) = 1.d0
      b5bn(i)   = 0.d0

      j = ny

      a5be(j,:) = a5bn(i,:)
      b5be(j)   = b5bn(i)

   end subroutine get_bc_corners_5d



   !> \brief Transfers the coeficients and source of the numerical scheme of the
   !! boundary conditions in the fictitious volumes to the matrix a5 and vector
   !! b5 that will be used to solve the linear system
   subroutine get_bc_transfer_5d(nx, ny, a5bn, a5bs, a5be, a5bw, b5bn, b5bs &
      , b5be, b5bw, a5, b5) ! Output: last two
      implicit none
      integer, intent(in) :: nx !< Number of volumes in csi direction (real+fictitious)
      integer, intent(in) :: ny !< Number of volumes in eta direction (real+fictitious)
      real(8), dimension(nx,5),    intent(in)  :: a5bn !< Coefficients of the discretization for the north boundary
      real(8), dimension(nx,5),    intent(in)  :: a5bs !< Coefficients of the discretization for the south boundary
      real(8), dimension(ny,5),    intent(in)  :: a5be !< Coefficients of the discretization for the east boundary
      real(8), dimension(ny,5),    intent(in)  :: a5bw !< Coefficients of the discretization for the west boundary
      real(8), dimension(nx),      intent(in)  :: b5bn !< Source of the discretization for the north boundary
      real(8), dimension(nx),      intent(in)  :: b5bs !< Source of the discretization for the south boundary
      real(8), dimension(ny),      intent(in)  :: b5be !< Source of the discretization for the east boundary
      real(8), dimension(ny),      intent(in)  :: b5bw !< Source of the discretization for the west boundary
      real(8), dimension(nx*ny,5), intent(out) :: a5   !< Coefficients of the discretization for all volumes
      real(8), dimension(nx*ny),   intent(out) :: b5   !< Source of the discretization for all volumes

      ! Inner variables

      integer :: i, j, np



      ! North boundary

      j = ny

      do i = 1, nx

         np = nx * (j-1) + i

         a5(np,:) = a5bn(i,:)

         b5(np) = b5bn(i)

      end do


      ! South boundary

      j = 1

      do i = 1, nx

         np = nx * (j-1) + i

         a5(np,:) = a5bs(i,:)

         b5(np) = b5bs(i)

      end do



      ! East boundary

      i = nx

      do j = 1, ny

         np = nx * (j-1) + i

         a5(np,:) = a5be(j,:)

         b5(np) = b5be(j)

      end do



      ! West boundary

      i = 1

      do j = 1, ny

         np = nx * (j-1) + i

         a5(np,:) = a5bw(j,:)

         b5(np) = b5bw(j)

      end do


   end subroutine get_bc_transfer_5d


end module bc5d
