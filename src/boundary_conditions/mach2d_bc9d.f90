
!> \brief Defines the subroutines related to the calculation of the boundary
!! conditions. More specifically, this module deals with numerical schemes that
!! envolves each fictitious volume and its neighbours (up to eight).
module bc9d
   implicit none

contains

   !> \brief Defines the numerical scheme for the Dirichlet boundary condition
   !! over the north boundary using a 9-diagonal matrix.
   subroutine get_bc_dirichlet_north_9d(nx, Fbn, a9bn, b9bn) ! Output: last two
      implicit none
      integer, intent(in) :: nx !< Number of volumes in csi direction (real+fictitious)
      real(8), dimension(nx),   intent(in)  :: Fbn  !< Field F over the north boundary
      real(8), dimension(nx,9), intent(out) :: a9bn !< Coefficients of the discretization for the north boundary
      real(8), dimension(nx),   intent(out) :: b9bn !< Source of the discretization for the north boundary

      ! Inner variables

      integer :: i

      ! North boundary volumes

      do i = 2, nx-1

         a9bn(i,:) = 0.d0
         a9bn(i,2) = 1.d0 ! S
         a9bn(i,5) = 1.d0 ! P

         b9bn(i) = 2.d0 * Fbn(i)

      end do

   end subroutine get_bc_dirichlet_north_9d


   !> \brief Defines the numerical scheme for the Dirichlet boundary condition
   !! over the south boundary using a 9-diagonal matrix.
   subroutine get_bc_dirichlet_south_9d(nx, Fbs, a9bs, b9bs) ! Output: last two
      implicit none
      integer, intent(in) :: nx !< Number of volumes in csi direction (real+fictitious)
      real(8), dimension(nx),   intent(in)  :: Fbs  !< Field F over the south boundary
      real(8), dimension(nx,9), intent(out) :: a9bs !< Coefficients of the discretization for the south boundary
      real(8), dimension(nx),   intent(out) :: b9bs !< Source of the discretization for the south boundary

      ! Inner variables

      integer :: i

      ! South boundary volumes

      do i = 2, nx-1

         a9bs(i,:) = 0.d0
         a9bs(i,5) = 1.d0 ! P
         a9bs(i,8) = 1.d0 ! N
         b9bs(i) = 2.d0 * Fbs(i)

      end do

   end subroutine get_bc_dirichlet_south_9d


   !> \brief Defines the numerical scheme for the Dirichlet boundary condition
   !! over the west boundary using a 9-diagonal matrix.
   subroutine get_bc_dirichlet_west_9d(ny, Fbw, a9bw, b9bw)
      implicit none
      integer, intent(in) :: ny !< Number of volumes in eta direction (real+fictitious)
      real(8), dimension(ny),   intent(in)  :: Fbw  !< Field F over the west boundary
      real(8), dimension(ny,9), intent(out) :: a9bw !< Coefficients of the discretization for the west boundary
      real(8), dimension(ny),   intent(out) :: b9bw !< Source of the discretization for the west boundary

      ! Inner variables

      integer :: j

      ! West boundary volumes

      do j = 2, ny-1

         a9bw(j,:) =  0.d0
         a9bw(j,5) =  1.d0 ! P
         a9bw(j,6) =  1.d0 ! E
         b9bw(j) = 2.d0 * Fbw(j)

      end do

   end subroutine get_bc_dirichlet_west_9d





   !> \brief Defines the numerical scheme for the null normal gradient boundary
   !! condition over the south boundary using a 9-diagonal matrix.
   !! CAUTION: This subroutine uses an inconsistent scheme, i.e., theoretically
   !! it does not reduce to the exact bc with mesh refining.
   subroutine get_ibc_null_normal_grad_south_9d(nx, a9bs, b9bs) ! Output: last two
      implicit none
      integer, intent(in) :: nx !< Number of volumes in csi direction (real+fictitious)
      real(8), dimension(nx,9), intent(out) :: a9bs !< Coefficients of the discretization for the south boundary
      real(8), dimension(nx),   intent(out) :: b9bs !< Source of the discretization for the south boundary

      ! Inner variables

      integer :: i

      ! South boundary volumes

      do i = 2, nx-1

         a9bs(i,:) =  0.d0
         a9bs(i,5) =  1.d0 ! P
         a9bs(i,8) = -1.d0 ! N
         b9bs(i) =  0.d0

      end do

   end subroutine get_ibc_null_normal_grad_south_9d


   !> \brief Defines the numerical scheme for the null normal gradient boundary
   !! condition over the west boundary using a 9-diagonal matrix.
   !! CAUTION: This subroutine uses an inconsistent scheme, i.e., theoretically
   !! it does not reduce to the exact bc with mesh refining.
   subroutine get_ibc_null_normal_grad_west_9d(ny, a9bw, b9bw) ! Output: last two
      implicit none
      integer, intent(in) :: ny !< Number of volumes in eta direction (real+fictitious)
      real(8), dimension(ny,9), intent(out) :: a9bw !< Coefficients of the discretization for the west boundary
      real(8), dimension(ny),   intent(out) :: b9bw !< Source of the discretization for the west boundary

      ! Inner variables

      integer :: j

      ! West boundary volumes

      do j = 2, ny-1

         a9bw(j,:) =  0.d0
         a9bw(j,5) =  1.d0 ! P
         a9bw(j,6) = -1.d0 ! E
         b9bw(j) = 0.d0

      end do

   end subroutine get_ibc_null_normal_grad_west_9d





   !> \brief Defines the numerical scheme for the null normal gradient boundary
   !! condition over the east boundary using a 9-diagonal matrix.
   !! CAUTION: This subroutine uses an inconsistent scheme, i.e., theoretically
   !! it does not reduce to the exact bc with mesh refining.
   subroutine get_ibc_null_normal_grad_east_9d(ny, a9be, b9be) ! Output: last two
      implicit none
      integer, intent(in) :: ny !< Number of volumes in eta direction (real+fictitious)
      real(8), dimension(ny,9), intent(out) :: a9be !< Coefficients of the discretization for the east boundary
      real(8), dimension(ny),   intent(out) :: b9be !< Source of the discretization for the east boundary

      ! Inner variables

      integer :: j

      ! East boundary volumes

      do j = 2, ny-1

         a9be(j,:) =  0.d0
         a9be(j,4) = -1.d0 ! W
         a9be(j,5) =  1.d0 ! P
         b9be(j) = 0.d0

      end do

   end subroutine get_ibc_null_normal_grad_east_9d


   !> \brief Defines the numerical scheme for the null normal gradient boundary
   !! condition over the south boundary using a 9-diagonal matrix.
   subroutine get_bc_null_normal_grad_south_9d(nx, ny, betan, gamman, a9bs, b9bs) ! Output: last two
      implicit none
      integer, intent(in) :: nx !< Number of volumes in csi direction (real+fictitious)
      integer, intent(in) :: ny !< Number of volumes in eta direction (real+fictitious)
      real(8), dimension (nx*ny), intent(in)  :: betan  !< (metric) beta  at the center of north face of volume P
      real(8), dimension (nx*ny), intent(in)  :: gamman !< (metric) gamma at the center of north face of volume P
      real(8), dimension(nx,9),   intent(out) :: a9bs   !< Coefficients of the discretization for the south boundary
      real(8), dimension(nx),     intent(out) :: b9bs   !< Source of the discretization for the south boundary

      ! Inner variables

      integer :: i, j, np

      real(8) :: raux


      j = 1

      ! South boundary (west volume)

      i = 2

      np   = nx * (j-1) + i

      raux = ( betan(np) / gamman(np) ) / 2.d0

      a9bs(i,1) =  0.d0        ! SW
      a9bs(i,2) =  0.d0        ! S
      a9bs(i,3) =  0.d0        ! SE
      a9bs(i,4) =  0.d0        ! W
      a9bs(i,5) =  1.d0 - raux ! P
      a9bs(i,6) =  raux        ! E
      a9bs(i,7) =  0.d0        ! NW
      a9bs(i,8) = -1.d0 - raux ! N
      a9bs(i,9) =  raux        ! NE
      b9bs(i)   =  0.d0


      ! South boundary (central volume)

      do i = 3, nx-2

         np   = nx * (j-1) + i

         raux = ( betan(np) / gamman(np) ) / 4.d0

         a9bs(i,1) =  0.d0 ! SW
         a9bs(i,2) =  0.d0 ! S
         a9bs(i,3) =  0.d0 ! SE
         a9bs(i,4) = -raux ! W
         a9bs(i,5) =  1.d0 ! P
         a9bs(i,6) =  raux ! E
         a9bs(i,7) = -raux ! NW
         a9bs(i,8) = -1.d0 ! N
         a9bs(i,9) =  raux ! NE
         b9bs(i)   =  0.d0

      end do


      ! South boundary (east volume)

      i = nx-1

      np   = nx * (j-1) + i

      raux = ( betan(np) / gamman(np) ) / 2.d0

      a9bs(i,1) =  0.d0        ! SW
      a9bs(i,2) =  0.d0        ! S
      a9bs(i,3) =  0.d0        ! SE
      a9bs(i,4) = -raux        ! W
      a9bs(i,5) =  1.d0 + raux ! P
      a9bs(i,6) =  0.d0        ! E
      a9bs(i,7) = -raux        ! NW
      a9bs(i,8) = -1.d0 + raux ! N
      a9bs(i,9) =  0.d0        ! NE
      b9bs(i)   =  0.d0

   end subroutine get_bc_null_normal_grad_south_9d


   !> \brief Defines the numerical scheme for the null normal gradient boundary
   !! condition over the west boundary using a 9-diagonal matrix.
   subroutine get_bc_null_normal_grad_west_9d(nx, ny, alphae, betae, a9bw, b9bw) ! Output: last two
      implicit none
      integer, intent(in) :: nx !< Number of volumes in csi direction (real+fictitious)
      integer, intent(in) :: ny !< Number of volumes in eta direction (real+fictitious)
      real(8), dimension(nx*ny), intent(in)  :: alphae !< (metric) alpha at the center of east face of volume P
      real(8), dimension(nx*ny), intent(in)  :: betae  !< (metric) beta  at the center of east face of volume P
      real(8), dimension(ny,9),  intent(out) :: a9bw   !< Coefficients of the discretization for the west boundary
      real(8), dimension(ny),    intent(out) :: b9bw   !< Source of the discretization for the west boundary

      ! Inner variables

      integer :: i, j, np

      real(8) :: raux


      i = 1

      ! West boundary (south volume)

      j = 2

      np   = nx * (j-1) + i

      raux = 0.5d0 * betae(np) / alphae(np)

      a9bw(j,1) =  0.d0        ! SW
      a9bw(j,2) =  0.d0        ! S
      a9bw(j,3) =  0.d0        ! SE
      a9bw(j,4) =  0.d0        ! W
      a9bw(j,5) =  1.d0 - raux ! P
      a9bw(j,6) = -1.d0 - raux ! E
      a9bw(j,7) =  0.d0        ! NW
      a9bw(j,8) =  raux        ! N
      a9bw(j,9) =  raux        ! NE
      b9bw(j)   =  0.d0


      ! West boundary (central volume)

      do j = 3, ny-2

         np   = nx * (j-1) + i

         raux = 0.25d0 * betae(np) / alphae(np)

         a9bw(j,1) =  0.d0 ! SW
         a9bw(j,2) = -raux ! S
         a9bw(j,3) = -raux ! SE
         a9bw(j,4) =  0.d0 ! W
         a9bw(j,5) =  1.d0 ! P
         a9bw(j,6) = -1.d0 ! E
         a9bw(j,7) =  0.d0 ! NW
         a9bw(j,8) =  raux ! N
         a9bw(j,9) =  raux ! NE
         b9bw(j)   =  0.d0

      end do

      ! West boundary (north volume)

      j = ny-1

      np   = nx * (j-1) + i

      raux = 0.5d0 * betae(np) / alphae(np)

      a9bw(j,1) =  0.d0        ! SW
      a9bw(j,2) = -raux        ! S
      a9bw(j,3) = -raux        ! SE
      a9bw(j,4) =  0.d0        ! W
      a9bw(j,5) =  1.d0 + raux ! P
      a9bw(j,6) = -1.d0 + raux ! E
      a9bw(j,7) =  0.d0        ! NW
      a9bw(j,8) =  0.d0        ! N
      a9bw(j,9) =  0.d0        ! NE
      b9bw(j)   =  0.d0

   end subroutine get_bc_null_normal_grad_west_9d



   !> \brief Defines the numerical scheme of the u velocity slip over the south
   !! boundary.
   subroutine get_bc_u_slip_south_9d(nx, ny, xk, yk, u, v, a9bs, b9bs) ! Output: last two
      implicit none
      integer, intent(in) :: nx !< Number of volumes in csi direction (real+fictitious)
      integer, intent(in) :: ny !< Number of volumes in eta direction (real+fictitious)
      real(8), dimension(nx*ny), intent(in)  :: xk   !< x_csi at face north of volume P (m)
      real(8), dimension(nx*ny), intent(in)  :: yk   !< y_csi at face north of volume P (m)
      real(8), dimension(nx*ny), intent(in)  :: u    !< x cartesian velocity (m/s)
      real(8), dimension(nx*ny), intent(in)  :: v    !< y cartesian velocity (m/s)
      real(8), dimension(nx,9),  intent(out) :: a9bs !< Coefficients of the discretization for the south boundary
      real(8), dimension(nx),    intent(out) :: b9bs !< Source of the discretization for the south boundary

      ! Inner variables

      integer :: i, j, np, npn

      real(8) :: lambda

      ! South boundary volumes

      j = 1

      do i = 2, nx-1

         a9bs(i,:) =  0.d0
         a9bs(i,5) =  1.d0 ! P
         a9bs(i,8) =  1.d0 ! N

         np   = nx * (j-1) + i

         npn  = np + nx

         lambda = sign(1.d0, xk(np) * u(npn) + yk(np) * v(npn) ) &

               * sqrt( ( u(npn) ** 2 + v(npn) ** 2 ) &

               / ( xk(np)**2 + yk(np)**2 ) )

         b9bs(i) =  2.d0 * lambda * xk(np)

      end do

   end subroutine get_bc_u_slip_south_9d


   !> \brief Defines the numerical scheme of the v velocity slip over the south
   !! boundary.
   subroutine get_bc_v_slip_south_9d(nx, ny, xk, yk, u, v, a9bs, b9bs) ! Output: last two
      implicit none
      integer, intent(in) :: nx !< Number of volumes in csi direction (real+fictitious)
      integer, intent(in) :: ny !< Number of volumes in eta direction (real+fictitious)
      real(8), dimension(nx*ny), intent(in)  :: xk   !< x_csi at face north of volume P (m)
      real(8), dimension(nx*ny), intent(in)  :: yk   !< y_csi at face north of volume P (m)
      real(8), dimension(nx*ny), intent(in)  :: u    !< x cartesian velocity (m/s)
      real(8), dimension(nx*ny), intent(in)  :: v    !< y cartesian velocity (m/s)
      real(8), dimension(nx,9),  intent(out) :: a9bs !< Coefficients of the discretization for the south boundary
      real(8), dimension(nx),    intent(out) :: b9bs !< Source of the discretization for the south boundary

      ! Inner variables

      integer :: i, j, np, npn

      real(8) :: lambda

      ! South boundary volumes

      j = 1

      do i = 2, nx-1

         a9bs(i,:) =  0.d0
         a9bs(i,5) =  1.d0 ! P
         a9bs(i,8) =  1.d0 ! N

         np   = nx * (j-1) + i

         npn  = np + nx

         lambda = sign(1.d0, xk(np) * u(npn) + yk(np) * v(npn) ) &

               * sqrt( ( u(npn) ** 2 + v(npn) ** 2 ) &

               / ( xk(np)**2 + yk(np)**2 ) )

         b9bs(i) =  2.d0 * lambda * yk(np)

      end do

   end subroutine get_bc_v_slip_south_9d



   !> \brief Defines the numerical scheme of the streamlined exit boundary
   !! condition over the south boundary.
   subroutine get_bc_streamlined_exit_east_9d(ny, Ucbe, Vcbe, a9be, b9be) ! Output: last two
      implicit none
      integer, intent(in) :: ny !< Number of volumes in eta direction (real+fictitious)
      real(8), dimension(ny),   intent(in)  :: Ucbe !< Uc over the faces of the east  boundary (m2/s)
      real(8), dimension(ny),   intent(in)  :: Vcbe !< Vc over the faces of the east  boundary (m2/s)
      real(8), dimension(ny,9), intent(out) :: a9be !< Coefficients of the discretization for the east boundary
      real(8), dimension(ny),   intent(out) :: b9be !< Source of the discretization for the east boundary


      ! Inner variables

      integer :: j

      real(8) :: raux


      ! East boundary (south volume)

      j = 2

      raux = 0.5d0 * Vcbe(j) / Ucbe(j)

      a9be(j,1) =  0.d0        ! SW
      a9be(j,2) =  0.d0        ! S
      a9be(j,3) =  0.d0        ! SE
      a9be(j,4) = -1.d0 - raux ! W
      a9be(j,5) =  1.d0 - raux ! P
      a9be(j,6) =  0.d0        ! E
      a9be(j,7) =  raux        ! NW
      a9be(j,8) =  raux        ! N
      a9be(j,9) =  0.d0        ! NE
      b9be(j)   =  0.d0

      ! East boundary (central volume)

      do j = 3, ny-2

         raux = 0.25d0 * Vcbe(j) / Ucbe(j)

         a9be(j,1) = -raux ! SW
         a9be(j,2) = -raux ! S
         a9be(j,3) =  0.d0 ! SE
         a9be(j,4) = -1.d0 ! W
         a9be(j,5) =  1.d0 ! P
         a9be(j,6) =  0.d0 ! E
         a9be(j,7) =  raux ! NW
         a9be(j,8) =  raux ! N
         a9be(j,9) =  0.d0 ! NE
         b9be(j)   =  0.d0

      end do


      ! East boundary (north volume)

      j = ny-1

      raux = 0.5d0 * Vcbe(j) / Ucbe(j)

      a9be(j,1) = -raux        ! SW
      a9be(j,2) = -raux        ! S
      a9be(j,3) =  0.d0        ! SE
      a9be(j,4) = -1.d0 + raux ! W
      a9be(j,5) =  1.d0 + raux ! P
      a9be(j,6) =  0.d0        ! E
      a9be(j,7) =  0.d0        ! NW
      a9be(j,8) =  0.d0        ! N
      a9be(j,9) =  0.d0        ! NE
      b9be(j)   =  0.d0

   end subroutine get_bc_streamlined_exit_east_9d



   !> \brief Defines the numerical scheme for the calculation of the boundary
   !! conditions at the fictitious volumes at corners of the transformed domain
   subroutine get_bc_corners_9d(nx, ny, a9bn, a9bs, a9be, a9bw, b9bn, b9bs &
         , b9be, b9bw ) ! Output: last eight
      implicit none
      integer, intent(in) :: nx !< Number of volumes in csi direction (real+fictitious)
      integer, intent(in) :: ny !< Number of volumes in eta direction (real+fictitious)
      real(8), dimension(nx,9), intent(out) :: a9bn !< Coefficients of the discretization for the north boundary
      real(8), dimension(nx,9), intent(out) :: a9bs !< Coefficients of the discretization for the south boundary
      real(8), dimension(ny,9), intent(out) :: a9be !< Coefficients of the discretization for the east boundary
      real(8), dimension(ny,9), intent(out) :: a9bw !< Coefficients of the discretization for the west boundary
      real(8), dimension(nx),   intent(out) :: b9bn !< Source of the discretization for the north boundary
      real(8), dimension(nx),   intent(out) :: b9bs !< Source of the discretization for the south boundary
      real(8), dimension(ny),   intent(out) :: b9be !< Source of the discretization for the east boundary
      real(8), dimension(ny),   intent(out) :: b9bw !< Source of the discretization for the west boundary

      ! Inner variables

      integer :: i, j

      ! SW corner

      i = 1

      a9bs(i,:) =  0.d0
      a9bs(i,5) =  1.d0        ! P
      a9bs(i,6) = -1.d0 / 3.d0 ! E
      a9bs(i,8) = -1.d0 / 3.d0 ! N
      a9bs(i,9) = -1.d0 / 3.d0 ! NE

      b9bs(i) = 0.d0

      j = 1

      a9bw(j,:) = a9bs(i,:)

      b9bw(j) = b9bs(i)


      ! SE corner

      i = nx

      a9bs(i,:) =  0.d0
      a9bs(i,4) = -1.d0 / 3.d0 ! W
      a9bs(i,5) =  1.d0        ! P
      a9bs(i,7) = -1.d0 / 3.d0 ! NW
      a9bs(i,8) = -1.d0 / 3.d0 ! N

      b9bs(i) = 0.d0

      j = 1

      a9be(j,:) = a9bs(i,:)

      b9be(j) = b9bs(i)


      ! NW corner

      i = 1

      a9bn(i,:) =  0.d0
      a9bn(i,2) = -1.d0 / 3.d0 ! S
      a9bn(i,3) = -1.d0 / 3.d0 ! SE
      a9bn(i,5) =  1.d0        ! P
      a9bn(i,6) = -1.d0 / 3.d0 ! E

      b9bn(i) = 0.d0

      j = ny

      a9bw(j,:) = a9bn(i,:)

      b9bw(j) = b9bn(i)


      ! NE corner

      i = nx

      a9bn(i,:) =  0.d0
      a9bn(i,1) = -1.d0 / 3.d0 ! SW
      a9bn(i,2) = -1.d0 / 3.d0 ! S
      a9bn(i,4) = -1.d0 / 3.d0 ! W
      a9bn(i,5) =  1.d0        ! P

      b9bn(i) = 0.d0

      j = ny

      a9be(j,:) = a9bn(i,:)

      b9be(j) = b9bn(i)


   end subroutine get_bc_corners_9d



   !> \brief Transfers the coeficients and source of the numerical scheme of the
   !! boundary conditions in the fictitious volumes to the matrix a9 and vector
   !! b9 that will be used to solve the linear system
   subroutine get_bc_transfer_9d(nx, ny, a9bn, a9bs, a9be, a9bw, b9bn, b9bs &
      , b9be, b9bw, a9, b9) ! Output: last two
      implicit none
      integer, intent(in) :: nx !< Number of volumes in csi direction (real+fictitious)
      integer, intent(in) :: ny !< Number of volumes in eta direction (real+fictitious)
      real(8), dimension(nx,9),    intent(in)  :: a9bn !< Coefficients of the discretization for the north boundary
      real(8), dimension(nx,9),    intent(in)  :: a9bs !< Coefficients of the discretization for the south boundary
      real(8), dimension(ny,9),    intent(in)  :: a9be !< Coefficients of the discretization for the east boundary
      real(8), dimension(ny,9),    intent(in)  :: a9bw !< Coefficients of the discretization for the west boundary
      real(8), dimension(nx),      intent(in)  :: b9bn !< Source of the discretization for the north boundary
      real(8), dimension(nx),      intent(in)  :: b9bs !< Source of the discretization for the south boundary
      real(8), dimension(ny),      intent(in)  :: b9be !< Source of the discretization for the east boundary
      real(8), dimension(ny),      intent(in)  :: b9bw !< Source of the discretization for the west boundary
      real(8), dimension(nx*ny,9), intent(out) :: a9   !< Coefficients of the discretization for all volumes
      real(8), dimension(nx*ny),   intent(out) :: b9   !< Source of the discretization for all volumes

      ! Inner variables

      integer :: i, j, np



      ! North boundary

      j = ny

      do i = 1, nx

         np = nx * (j-1) + i

         a9(np,:) = a9bn(i,:)

         b9(np) = b9bn(i)

      end do


      ! South boundary

      j = 1

      do i = 1, nx

         np = nx * (j-1) + i

         a9(np,:) = a9bs(i,:)

         b9(np) = b9bs(i)

      end do



      ! East boundary

      i = nx

      do j = 1, ny

         np = nx * (j-1) + i

         a9(np,:) = a9be(j,:)

         b9(np) = b9be(j)

      end do



      ! West boundary

      i = 1

      do j = 1, ny

         np = nx * (j-1) + i

         a9(np,:) = a9bw(j,:)

         b9(np) = b9bw(j)

      end do


   end subroutine get_bc_transfer_9d


   !> \brief Extrapolates the field F from the real volumes to the fictitious
   !! ones according to the boundary conditions
   subroutine get_bc_extrapolation_9d(nx, ny, itemax, a9bn, a9bs, a9be, a9bw, b9bn, b9bs &
      , b9be, b9bw, F) ! InOutput: last one
      implicit none
      integer, intent(in) :: nx !< Number of volumes in csi direction (real+fictitious)
      integer, intent(in) :: ny !< Number of volumes in eta direction (real+fictitious)
      integer, intent(in) :: itemax !< Number of iteractions for extrapolation to fictitious
      real(8), dimension(nx,9),    intent(in)    :: a9bn !< Coefficients of the discretization for the north boundary
      real(8), dimension(nx,9),    intent(in)    :: a9bs !< Coefficients of the discretization for the south boundary
      real(8), dimension(ny,9),    intent(in)    :: a9be !< Coefficients of the discretization for the east boundary
      real(8), dimension(ny,9),    intent(in)    :: a9bw !< Coefficients of the discretization for the west boundary
      real(8), dimension(nx),      intent(in)    :: b9bn !< Source of the discretization for the north boundary
      real(8), dimension(nx),      intent(in)    :: b9bs !< Source of the discretization for the south boundary
      real(8), dimension(ny),      intent(in)    :: b9be !< Source of the discretization for the east boundary
      real(8), dimension(ny),      intent(in)    :: b9bw !< Source of the discretization for the west boundary
      real(8), dimension(nx*ny),   intent(inout) :: F    !< Field to be extrapolated

      ! Inner variables

      integer :: i, j, k, np, npw, npe, nps, npsw, npse, npn, npnw, npne

      real(8), dimension(nx,3) :: ax !< Coefficients of the linear system for extrapolation
      real(8), dimension(ny,3) :: ay !< Coefficients of the linear system for extrapolation
      real(8), dimension(nx)   :: bx !< Source of the linear system for extrapolation
      real(8), dimension(ny)   :: by !< Source of the linear system for extrapolation
      real(8), dimension(nx)   :: sx !< Solution for x
      real(8), dimension(ny)   :: sy !< Solution for y


      do k = 1, itemax

         ! North boundary

         ax(:,1) = a9bn(:,4) ! W
         ax(:,2) = a9bn(:,5) ! P
         ax(:,3) = a9bn(:,6) ! E

         j = ny

         i = 1

         np   = nx * (j-1) + i
         nps  = np - nx
         npse = nps + 1

         bx(i) = b9bn(i) - ( a9bn(i,2) * F(nps) + a9bn(i,3) * F(npse) )

         do i = 2, nx-1

            np   = nx * (j-1) + i
            nps  = np - nx
            npsw = nps - 1
            npse = nps + 1

            bx(i) = b9bn(i) &

               - ( a9bn(i,1) * F(npsw) + a9bn(i,2) * F(nps) + a9bn(i,3) * F(npse) )

         end do

         i = nx

         np   = nx * (j-1) + i
         nps  = np - nx
         npsw = nps - 1

         bx(i) = b9bn(i) - ( a9bn(i,1) * F(npsw) + a9bn(i,2) * F(nps) )

         ! Solving the linear system

         call tdma(nx, ax, bx, sx) ! Output: last one

         ! Extrapolating to fictitious

         j = ny

         do i = 1, nx

            np   = nx * (j-1) + i

            F(np) = sx(i)

         end do



         ! South boundary

         ax(:,1) = a9bs(:,4) ! W
         ax(:,2) = a9bs(:,5) ! P
         ax(:,3) = a9bs(:,6) ! E

         j = 1

         i = 1

         np   = nx * (j-1) + i
         npn  = np + nx
         npne = npn + 1

         bx(i) = b9bs(i) - ( a9bs(i,8) * F(npn) + a9bs(i,9) * F(npne) )


         do i = 2, nx-1

            np   = nx * (j-1) + i
            npn  = np + nx
            npnw = npn - 1
            npne = npn + 1

            bx(i) = b9bs(i) &

               - ( a9bs(i,7) * F(npnw) + a9bs(i,8) * F(npn) + a9bs(i,9) * F(npne) )

         end do

         i = nx

         np   = nx * (j-1) + i
         npn  = np + nx
         npnw = npn - 1

         bx(i) = b9bs(i) - ( a9bs(i,7) * F(npnw) + a9bs(i,8) * F(npn) )


         ! Solving the linear system

         call tdma(nx, ax, bx, sx) ! Output: last one

         ! Extrapolating to fictitious

         j = 1

         do i = 1, nx

            np = nx * (j-1) + i

            F(np) = sx(i)

         end do



         ! East boundary

         ay(:,1) = a9be(:,2) ! S
         ay(:,2) = a9be(:,5) ! P
         ay(:,3) = a9be(:,8) ! N

         i = nx

         j = 1

         np   = nx * (j-1) + i
         npn  = np + nx
         npw  = np - 1
         npnw = npn - 1

         by(j) = b9be(j) - ( a9be(j,4) * F(npw) + a9be(j,7) * F(npnw) )

         do j = 2, ny-1

            np   = nx * (j-1) + i
            nps  = np - nx
            npn  = np + nx
            npw  = np - 1
            npsw = nps - 1
            npnw = npn - 1

            by(j) = b9be(j) &

               - ( a9be(j,1) * F(npsw) + a9be(j,4) * F(npw) + a9be(j,7) * F(npnw) )

         end do

         j = ny

         np   = nx * (j-1) + i
         nps  = np - nx
         npw  = np - 1
         npsw = nps - 1

         by(j) = b9be(j) - ( a9be(j,1) * F(npsw) + a9be(j,4) * F(npw) )


         ! Solving the linear system

         call tdma(ny, ay, by, sy) ! Output: last one

         ! Extrapolating to fictitious

         i = nx

         do j = 1, ny

            np = nx * (j-1) + i

            F(np) = sy(j)

         end do


         ! West boundary

         ay(:,1)= a9bw(:,2) ! S
         ay(:,2)= a9bw(:,5) ! P
         ay(:,3)= a9bw(:,8) ! N

         i = 1

         j = 1

         np   = nx * (j-1) + i
         npn  = np + nx
         npe  = np + 1
         npne = npn + 1

         by(j) = b9bw(j) - ( a9bw(j,6) * F(npe) + a9bw(j,9) * F(npne) )

         do j = 2, ny-1

            np   = nx * (j-1) + i
            nps  = np - nx
            npn  = np + nx
            npe  = np + 1
            npse = nps + 1
            npne = npn + 1

            by(j) = b9bw(j) &

               - ( a9bw(j,3) * F(npse) + a9bw(j,6) * F(npe) + a9bw(j,9) * F(npne) )

         end do

         j = ny

         np   = nx * (j-1) + i
         nps  = np - nx
         npe  = np + 1
         npse = nps + 1

         by(j) = b9bw(j) - ( a9bw(j,3) * F(npse) + a9bw(j,6) * F(npe) )


         ! Solving the linear system

         call tdma(ny, ay, by, sy) ! Output: last one

         ! Extrapolating to fictitious

         i = 1

         do j = 1, ny

            np = nx * (j-1) + i

            F(np) = sy(j)

         end do

      end do

   end subroutine get_bc_extrapolation_9d

   !> \brief Verifies the subroutine get_bc_extrapolation_9d
   subroutine get_bc_extrapolation_9d_check(itemax)
      implicit none
      integer, intent(in) :: itemax !< Number of iteractions for extrapolation to fictitious

      ! Inner variables

      integer, parameter :: nx = 8 ! Number of volumes in csi direction (real+fictitious)
      integer, parameter :: ny = 4 ! Number of volumes in eta direction (real+fictitious)

      real(8), dimension(nx,9)  :: a9bn ! Coefficients of the discretization for the north boundary
      real(8), dimension(nx,9)  :: a9bs ! Coefficients of the discretization for the south boundary
      real(8), dimension(ny,9)  :: a9be ! Coefficients of the discretization for the east boundary
      real(8), dimension(ny,9)  :: a9bw ! Coefficients of the discretization for the west boundary
      real(8), dimension(nx)    :: b9bn ! Source of the discretization for the north boundary
      real(8), dimension(nx)    :: b9bs ! Source of the discretization for the south boundary
      real(8), dimension(ny)    :: b9be ! Source of the discretization for the east boundary
      real(8), dimension(ny)    :: b9bw ! Source of the discretization for the west boundary
      real(8), dimension(nx*ny) :: F    ! Field to be extrapolated


      integer :: i, j, np

      real(8) :: raux = 1.d0

      ! Initializing F

      F = 0.d0

      ! North boundary

      j = ny-1

      do i = 2, nx-1

         np   = nx * (j-1) + i

         F(np) = raux

      end do

      ! South boundary

      j = 2

      do i = 2, nx-1

         np   = nx * (j-1) + i

         F(np) = raux

      end do

      ! East boundary

      i = nx-1

      do j = 2, ny-1

         np   = nx * (j-1) + i

         F(np) = raux

      end do

      ! West boundary

      i = 2

      do j = 2, ny-1

         np   = nx * (j-1) + i

         F(np) = raux

      end do


      ! Initializing a and b


      ! North boundary

      i = 1
      a9bn(i,1) = 0.d0
      a9bn(i,2) = 2.d0
      a9bn(i,3) = 3.d0
      a9bn(i,4) = 0.d0
      a9bn(i,5) = 5.d0
      a9bn(i,6) = 6.d0
      a9bn(i,7) = 0.d0
      a9bn(i,8) = 0.d0
      a9bn(i,9) = 0.d0

      do i = 2, nx-1

         a9bn(i,1) = 1.d0
         a9bn(i,2) = 2.d0
         a9bn(i,3) = 3.d0
         a9bn(i,4) = 4.d0
         a9bn(i,5) = 5.d0
         a9bn(i,6) = 6.d0
         a9bn(i,7) = 0.d0
         a9bn(i,8) = 0.d0
         a9bn(i,9) = 0.d0

      end do

      i = nx
      a9bn(i,1) = 1.d0
      a9bn(i,2) = 2.d0
      a9bn(i,3) = 0.d0
      a9bn(i,4) = 4.d0
      a9bn(i,5) = 5.d0
      a9bn(i,6) = 0.d0
      a9bn(i,7) = 0.d0
      a9bn(i,8) = 0.d0
      a9bn(i,9) = 0.d0

      do i = 1, nx

         b9bn(i) = sum(a9bn(i,:)) * raux

      end do



      ! South boundary

      i = 1
      a9bs(i,1) = 0.d0
      a9bs(i,2) = 0.d0
      a9bs(i,3) = 0.d0
      a9bs(i,4) = 0.d0
      a9bs(i,5) = 5.d0
      a9bs(i,6) = 6.d0
      a9bs(i,7) = 0.d0
      a9bs(i,8) = 8.d0
      a9bs(i,9) = 9.d0

      do i = 2, nx-1

         a9bs(i,1) = 0.d0
         a9bs(i,2) = 0.d0
         a9bs(i,3) = 0.d0
         a9bs(i,4) = 4.d0
         a9bs(i,5) = 5.d0
         a9bs(i,6) = 6.d0
         a9bs(i,7) = 7.d0
         a9bs(i,8) = 8.d0
         a9bs(i,9) = 9.d0

      end do

      i = nx
      a9bs(i,1) = 0.d0
      a9bs(i,2) = 0.d0
      a9bs(i,3) = 0.d0
      a9bs(i,4) = 4.d0
      a9bs(i,5) = 5.d0
      a9bs(i,6) = 0.d0
      a9bs(i,7) = 7.d0
      a9bs(i,8) = 8.d0
      a9bs(i,9) = 0.d0

      do i = 1, nx

         b9bs(i) = sum(a9bs(i,:)) * raux

      end do


      ! East boundary

      j = 1
      a9be(j,1) = 0.d0
      a9be(j,2) = 0.d0
      a9be(j,3) = 0.d0
      a9be(j,4) = 4.d0
      a9be(j,5) = 5.d0
      a9be(j,6) = 0.d0
      a9be(j,7) = 7.d0
      a9be(j,8) = 8.d0
      a9be(j,9) = 0.d0

      do j = 2, ny-1

         a9be(j,1) = 1.d0
         a9be(j,2) = 2.d0
         a9be(j,3) = 0.d0
         a9be(j,4) = 4.d0
         a9be(j,5) = 5.d0
         a9be(j,6) = 0.d0
         a9be(j,7) = 7.d0
         a9be(j,8) = 8.d0
         a9be(j,9) = 0.d0

      end do

      j = ny
      a9be(j,1) = 1.d0
      a9be(j,2) = 2.d0
      a9be(j,3) = 0.d0
      a9be(j,4) = 4.d0
      a9be(j,5) = 5.d0
      a9be(j,6) = 0.d0
      a9be(j,7) = 0.d0
      a9be(j,8) = 0.d0
      a9be(j,9) = 0.d0

      do j = 1, ny

         b9be(j) = sum(a9be(j,:)) * raux

      end do

      ! West boundary

      j = 1
      a9bw(j,1) = 0.d0
      a9bw(j,2) = 0.d0
      a9bw(j,3) = 0.d0
      a9bw(j,4) = 0.d0
      a9bw(j,5) = 5.d0
      a9bw(j,6) = 6.d0
      a9bw(j,7) = 0.d0
      a9bw(j,8) = 8.d0
      a9bw(j,9) = 9.d0

      do j = 2, ny-1

         a9bw(j,1) = 0.d0
         a9bw(j,2) = 2.d0
         a9bw(j,3) = 3.d0
         a9bw(j,4) = 0.d0
         a9bw(j,5) = 5.d0
         a9bw(j,6) = 6.d0
         a9bw(j,7) = 0.d0
         a9bw(j,8) = 8.d0
         a9bw(j,9) = 9.d0

      end do

      j = ny
      a9bw(j,1) = 0.d0
      a9bw(j,2) = 2.d0
      a9bw(j,3) = 3.d0
      a9bw(j,4) = 0.d0
      a9bw(j,5) = 5.d0
      a9bw(j,6) = 6.d0
      a9bw(j,7) = 0.d0
      a9bw(j,8) = 0.d0
      a9bw(j,9) = 0.d0

      do j = 1, ny

         b9bw(j) = sum(a9bw(j,:)) * raux

      end do

      call get_bc_extrapolation_9d(nx, ny, itemax, a9bn, a9bs, a9be, a9bw, b9bn, b9bs &
         , b9be, b9bw, F) ! InOutput: last one


      write(*,*) "South boundary"

      j = 1

      do i = 1, nx

         np   = nx * (j-1) + i

         write(*,*) i, F(np), raux, (F(np)-raux)/raux

      end do

      write(*,*) "North boundary"

      j = ny

      do i = 1, nx

         np   = nx * (j-1) + i

         write(*,*) i, F(np), raux, (F(np)-raux)/raux

      end do

      write(*,*) "West boundary"

      i = 1

      do j = 1, ny

         np   = nx * (j-1) + i

         write(*,*) j, F(np), raux, (F(np)-raux)/raux

      end do

      write(*,*) "East boundary"

      i = nx

      do j = 1, ny

         np   = nx * (j-1) + i

         write(*,*) j, F(np), raux, (F(np)-raux)/raux

      end do


   end subroutine get_bc_extrapolation_9d_check



   !> \brief Solves a tri-diagonal linear system
   subroutine tdma(n, a, b, x)
      implicit none
      integer, intent(in) :: n !< Number unknowns
      real(8), dimension(n,3), intent(in)  :: a !< Tri-diagonal matrix
      real(8), dimension(n),   intent(in)  :: b !< Source
      real(8), dimension(n),   intent(out) :: x !< Solution
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

end module bc9d
