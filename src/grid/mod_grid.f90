!>
!! \brief mod_grid also gives an interface to the main program to calculate the
!!        grid and related metrics. mod_grid uses procedures from mod_grid_procedure
!!        to operate data of mod_grid_data.
!!
module mod_grid

   use mod_class_ifile
   use mod_grid_procedures
   use mod_extflow
   use mod_intflow

   implicit none

   ! Makes everything private, except otherwise stated
   private

   ! Public methods
   public :: grid_init   &
      ,      grid_create

   !> \brief Parameters for grid generation
   type, public :: type_gridpar
      integer :: nx    !< number of volumes in the csi direction (real+fictitious)
      integer :: ny    !< number of volumes in the eta direction (real+fictitious)
      integer :: nxi   !< number of volumes in the csi direction (real+fictitious) of the coarsest grid
      integer :: nyi   !< number of volumes in the eta direction (real+fictitious) of the coarsest grid
      integer :: nxf   !< number of volumes in the csi direction (real+fictitious) of the finest grid
      integer :: nyf   !< number of volumes in the eta direction (real+fictitious) of the finest grid
      integer :: nmf   !< number of the finest mesh
      integer :: nmd   !< number of the desired mesh
      integer :: kg    !< Kind of grid (1=uniform, 2=geometric progression, 3=power law, 4=gp modified)
      integer :: kcm   !< Kind of centroid mean (1=simple mean, 2=weighted mean)
      real(8) :: avi   !< Initial value of the artificial viscosity
      real(8) :: avf   !< Final value of the artificial viscosity
      real(8) :: awf   !< Area weighting factor
      real(8) :: cbl   !< The width of the vol. closer to the wall is 'cbl' times the width of the b. layer
      real(8) :: a1    !< Width of the vol. closer to the wall (this values is used if cbl < 0)
   end type

contains

   !> \brief Initializes data of type type_gridpar
   subroutine grid_init(ifile, gridpar) ! Output: last one
      implicit none
      class(class_ifile),  intent(in)  :: ifile   !< Input file
      type(type_gridpar),  intent(out) :: gridpar !< Grid parameters

      call ifile%get_value( gridpar%nxi,    "nxi-2") ! Number of real volumes in the csi direction of the coarsest grid
      call ifile%get_value( gridpar%nyi,    "nyi-2") ! Number of real volumes in the eta direction of the coarsest grid
      call ifile%get_value( gridpar%nmf,      "nmf") ! Number of the finest mesh  (1<=nmf)
      call ifile%get_value( gridpar%nmd,      "nmd") ! Number of the desired mesh (1<=nmd<=nmf)
      call ifile%get_value( gridpar%kg,        "kg") ! Kind of grid (1=uniform, 2=geometric progression, 3=power law, 4=gp modified, 5=hyperbolic)
      call ifile%get_value( gridpar%avi,      "avi") ! Initial value of the artificial viscosity (only for kg=5)
      call ifile%get_value( gridpar%avf,      "avf") ! Final value of the artificial viscosity (only for kg=5)
      call ifile%get_value( gridpar%awf,      "awf") ! Area weighting factor (only for kg=5)
      call ifile%get_value( gridpar%kcm,      "kcm") ! Kind of centroid mean (1=simple mean, 2=weighted mean)
      call ifile%get_value( gridpar%cbl,      "cbl") ! The width of the vol. closer to the wall is 'cbl' times the width of the b. layer
      call ifile%get_value( gridpar%a1,        "a1") ! Width of the vol. closer to the wall (this values is used if cbl < 0)


      ! After reading the values of nxi and nyi, these variables must be
      ! changed to take into account the fictitious volumes
      gridpar%nxi = gridpar%nxi + 2
      gridpar%nyi = gridpar%nyi + 2

      ! Calculating the number of volumes (real+fictitious) of the finest mesh
      gridpar%nxf = (gridpar%nxi-2) * 2 ** (gridpar%nmf-1) + 2
      gridpar%nyf = (gridpar%nyi-2) * 2 ** (gridpar%nmf-1) + 2

      ! Calculating the number of volumes (real+fictitious) of the desired mesh
      gridpar%nx = (gridpar%nxi-2) * 2 ** (gridpar%nmd-1) + 2
      gridpar%ny = (gridpar%nyi-2) * 2 ** (gridpar%nmd-1) + 2

   end subroutine


   !> \brief Creates the grid and related metrics
   subroutine grid_create(unt, coord, kflow, gridpar           & ! Input
         ,  x, y, xp, yp, xe, ye, xen, yen, xk, yk, xke, yke, Jp      & ! Output
         ,  Je, Jn, alphae, gamman, betae, betan, radius, re, rn, rp, iflow  ) ! Output
      implicit none
      integer,               intent(in)  :: unt    !< Unit where the input parameters will be printed
      integer,               intent(in)  :: coord  !< Coordinate system ( 1=cylindrical, 0 = cartesian)
      integer,               intent(in)  :: kflow  !< Kind of flow (internal or external)
      type(type_gridpar),    intent(in)  :: gridpar!< Parameters related to the grid
      real(8), dimension(:), intent(out) :: x      !< x coord. at the northeast corner of the volume P
      real(8), dimension(:), intent(out) :: y      !< y coord. at the northeast corner of the volume P
      real(8), dimension(:), intent(out) :: xp     !< Coord. x of the centroid of volume P
      real(8), dimension(:), intent(out) :: yp     !< Coord. y of the centroid of volume P
      real(8), dimension(:), intent(out) :: xe     !< x_eta at face east of volume P
      real(8), dimension(:), intent(out) :: ye     !< y_eta at face east of volume P
      real(8), dimension(:), intent(out) :: xen    !< x_eta at face north of volume P
      real(8), dimension(:), intent(out) :: yen    !< y_eta at face north of volume P
      real(8), dimension(:), intent(out) :: xk     !< x_csi at face north of volume P
      real(8), dimension(:), intent(out) :: yk     !< y_csi at face north of volume P
      real(8), dimension(:), intent(out) :: xke    !< x_csi at face east of volume P
      real(8), dimension(:), intent(out) :: yke    !< y_csi at face east of volume P
      real(8), dimension(:), intent(out) :: Jp     !< Jacobian at the center of volume P
      real(8), dimension(:), intent(out) :: Je     !< Jacobian at the center of east face of volume P
      real(8), dimension(:), intent(out) :: Jn     !< Jacobian at the center of north face of volume P
      real(8), dimension(:), intent(out) :: Alphae !< (metric) Alpha at the center of east face of volume P
      real(8), dimension(:), intent(out) :: Gamman !< (metric) Gamma at the center of north face of volume P
      real(8), dimension(:), intent(out) :: Betae  !< (metric) Beta  at the center of east face of volume P
      real(8), dimension(:), intent(out) :: Betan  !< (metric) Beta  at the center of north face of volume P
      real(8), dimension(:), intent(out) :: radius !< Radius of northest corner of volume P
      real(8), dimension(:), intent(out) :: re     !< Radius of the center of east face of volume P
      real(8), dimension(:), intent(out) :: rn     !< Radius of the center of north face of volume P
      real(8), dimension(:), intent(out) :: rp     !< Radius of the center of volume P
      type(type_intflow),  intent(inout) :: iflow  !< Variables related to internal flow

      logical :: isReversed ! Reverses the partitioning distribution if true
      real(8) :: wbl        ! Boundary layer estimated width (m)
      real(8) :: a1         ! width of the volume closer to the wall (m)
      real(8), allocatable, dimension(:) :: xf  !< coord. at the northeast corner of the volume P of the finest grid
      real(8), allocatable, dimension(:) :: yf  !< coord. at the northeast corner of the volume P of the finest grid

      ! Allocating the finest grid
      allocate( xf(gridpar%nxf*gridpar%nyf) )
      allocate( yf(gridpar%nxf*gridpar%nyf) )

      ! Getting the grid boundary for the finest grid
      if ( kflow == 0 ) then ! External flow

         ! Points must be concentrated near the south boundary
         isReversed = .false.

         ! Generates the grid north and south boundary
         call extflow_grid_boundary(gridpar%nxf, gridpar%nyf, unt, xf, yf) ! Output: last four

         ! Estimates the boundary layer width and the width of the volume
         ! closer to the wall
         call extflow_boundary_layer(gridpar%cbl, wbl, a1) ! Output: last two

      else ! Internal flow

         ! Points must be concentrated near the south boundary
         isReversed = .true.

         a1 = gridpar%a1

         call intflow_grid_boundary(gridpar%nxf, gridpar%nyf, xf, yf, iflow) ! Output: last three

      end if

      ! If cbl is negative, use a1 from input file
      if ( gridpar%cbl < 0.d0 ) a1 = gridpar%a1

      ! Generates the grid according to kg option
      call set_grid(  isReversed   &
         ,            gridpar%kg   &
         ,            gridpar%nxf  &
         ,            gridpar%nyf  &
         ,            a1           &
         ,            gridpar%avi  &
         ,            gridpar%avf  &
         ,            gridpar%awf  &
         ,            xf, yf       ) ! Last 2 are inoutput

      ! Selecting the desided grid
      call get_coarser_grid( gridpar%nxf   &
         ,                   gridpar%nyf   &
         ,                   gridpar%nx    &
         ,                   gridpar%ny    &
         ,                   gridpar%nmf   &
         ,                   gridpar%nmd   &
         ,                   xf, yf        &
         ,                   x, y          )


      ! Removing the (now) unnecessary finest grid
      deallocate(xf, yf)

      ! Calculates the centroids of each real volume
      call get_real_centroids_xy(gridpar%kcm, gridpar%nx, gridpar%ny, x, y, xp, yp ) ! Output: last two
      ! checked

      ! Calculates the components of the metric tensor and the Jacobian
      call get_metrics(gridpar%nx, gridpar%ny, x, y, xp, yp & ! Input
      ,              xe, ye, xen, yen, xk, yk, xke, yke, Jp & ! Output
      ,                Je, Jn, alphae, gamman, betae, betan ) ! Output
      ! checked

      ! All the radii are equal 1 if planar flow is choose. If axisymmetric flow, radii are calculated
      call get_radius(coord, gridpar%nx, gridpar%ny, y, yp  & ! Input
      ,                                 radius, re, rn, rp  ) ! Output
      ! checked

   end subroutine


   !> \brief Generates a coarser grid from a given grid
   subroutine get_coarser_grid(nxf, nyf, nx, ny, nmf, nmd, xf, yf, x, y)
      implicit none
      integer, intent(in)  :: nxf
      integer, intent(in)  :: nyf
      integer, intent(in)  :: nx
      integer, intent(in)  :: ny
      integer, intent(in)  :: nmf
      integer, intent(in)  :: nmd
      real(8), intent(in)  :: xf(nxf*nyf)
      real(8), intent(in)  :: yf(nxf*nyf)
      real(8), intent(out) :: x(nx*ny)
      real(8), intent(out) :: y(nx*ny)


      integer :: i, j, k, ii, jj, np, npf


      k = nmf-nmd

      do j = 1, ny-1

         jj = 2 ** k * (j-1) + 1

         do i = 1, nx-1

            ii = 2 ** k * (i-1) + 1

            np  = nx  * (j -1) + i

            npf = nxf * (jj-1) + ii

            x(np) = xf(npf)

            y(np) = yf(npf)

         end do

      end do


   end subroutine


   !> \brief Calculates the centroids of the volumes
   subroutine get_real_centroids_xy( opt, nx, ny, x, y & ! Input
      ,                            xp, yp )            ! Output
      implicit none
      integer, intent(in) :: opt ! (1=simple mean, 2=weighted mean)
      integer, intent(in) :: nx  ! Number of volumes in csi direction (real+fictitious)
      integer, intent(in) :: ny  ! Number of volumes in eta direction (real+fictitious)
      real(8), intent(in) :: x(nx*ny) ! Coord. x of the northest corner of volume P
      real(8), intent(in) :: y(nx*ny) ! Coord. y of the northest corner of volume P
      real(8), intent(out) :: xp(nx*ny) ! Coord. x of the centroid of volume P
      real(8), intent(out) :: yp(nx*ny) ! Coord. y of the centroid of volume P
      !
      integer :: i, j, np, nps, npw, npsw
      real(8) :: A1, A2

      ! Centroids of real volumes

      xp = 0.d0
      yp = 0.d0

      if ( opt == 1 ) then ! Simple mean
         ! Centroids are given by the mean of the corners
         do j = 2, ny-1
            do i = 2, nx-1

               np   = nx * (j-1) + i
               nps  = np - nx
               npw  = np - 1
               npsw = nps - 1

               xp(np) = ( x(np) + x(npw) + x(npsw) + x(nps) ) / 4.d0
               yp(np) = ( y(np) + y(npw) + y(npsw) + y(nps) ) / 4.d0

            end do
         end do
      else ! weighted mean
         ! The rectangle is divided in two triangles of area A1 and A2.
         ! The centroids of each triangle is calculated by the mean of the corners.
         ! The centroid of the rectangle is given by the weighted mean of the centroids of each triangle.
         do j = 2, ny-1
            do i = 2, nx-1

               np   = nx * (j-1) + i
               nps  = np - nx
               npw  = np - 1
               npsw = nps - 1

               A1 = 0.5d0 * ( x(npsw) * (y(nps)-y(np) ) &
                  +          x(nps) * (y(np) -y(npsw)) &
                  +          x(np)  * (y(npsw)-y(nps)))

               A2 = 0.5d0 * ( x(npsw) * (y(np) -y(npw)) &
                  +          x(np)  * (y(npw)-y(npsw)) &
                  +          x(npw) * (y(npsw)-y(np) ))

               xp(np) = ( (x(npsw)+x(nps)+x(np)) * A1 + (x(npsw)+x(np)+x(npw)) * A2 ) &
                  / ( 3.d0 * (A1+A2) )

               yp(np) = ( (y(npsw)+y(nps)+y(np)) * A1 + (y(npsw)+y(np)+y(npw)) * A2 ) &
                  / ( 3.d0 * (A1+A2) )

            end do
         end do
      end if
   end subroutine get_real_centroids_xy


   !> \brief Calculates the metrics of coordinate transformation
   subroutine get_metrics(nx, ny, x, y, xp, yp            & ! Input
      ,          xe, ye, xen, yen, xk, yk, xke, yke, Jp & ! Output
      ,          Je, Jn, alphae, gamman, betae, betan )   ! Output
      implicit none
      integer, intent(in)  :: nx, ny  ! Number of volumes in csi and eta directions (real+fictitious)
      real(8), dimension (nx*ny), intent(in)  :: x,  y  ! Coord. of the northest corner of volume P
      real(8), dimension (nx*ny), intent(in)  :: xp, yp ! Coord. of the centroid of volume P
      real(8), dimension (nx*ny), intent(out) :: xe     ! x_eta at face east of volume P
      real(8), dimension (nx*ny), intent(out) :: ye     ! y_eta at face east of volume P
      real(8), dimension (nx*ny), intent(out) :: xen    ! x_eta at face north of volume P
      real(8), dimension (nx*ny), intent(out) :: yen    ! y_eta at face north of volume P
      real(8), dimension (nx*ny), intent(out) :: xk     ! x_csi at face north of volume P
      real(8), dimension (nx*ny), intent(out) :: yk     ! y_csi at face north of volume P
      real(8), dimension (nx*ny), intent(out) :: xke    ! x_csi at face east of volume P
      real(8), dimension (nx*ny), intent(out) :: yke    ! y_csi at face east of volume P
      real(8), dimension (nx*ny), intent(out) :: Jp     ! Jacobian at the center of volume P
      real(8), dimension (nx*ny), intent(out) :: Je     ! Jacobian at the center of east face of volume P
      real(8), dimension (nx*ny), intent(out) :: Jn     ! Jacobian at the center of north face of volume P
      real(8), dimension (nx*ny), intent(out) :: Alphae ! (metric) Alpha at the center of east face of volume P
      real(8), dimension (nx*ny), intent(out) :: Gamman ! (metric) Gamma at the center of north face of volume P
      real(8), dimension (nx*ny), intent(out) :: Betae  ! (metric) Beta  at the center of east face of volume P
      real(8), dimension (nx*ny), intent(out) :: Betan  ! (metric) Beta  at the center of north face of volume P
      !
      integer :: i, j, npsw, nps, npw, np, npe, npn
      real(8) :: fw, fe, fs, fn, xkp, xep, ykp, yep
      !
      xe  = 0.d0
      ye  = 0.d0
      xen = 0.d0
      yen = 0.d0
      xk  = 0.d0
      yk  = 0.d0
      xke = 0.d0
      yke = 0.d0
      Jp  = 0.d0
      Je  = 0.d0
      Jn  = 0.d0
      alphae = 0.d0
      gamman = 0.d0
      betae  = 0.d0
      betan  = 0.d0

      ! Derivatives relatively to eta at the center of the east face

      do j = 2, ny-1
         do i = 1, nx-1

            np  = nx * (j-1) + i
            nps = np - nx

            xe(np) = x(np) - x(nps)
            ye(np) = y(np) - y(nps)

            alphae(np) = xe(np)**2 + ye(np)**2

         end do
      end do

      ! Derivatives relatively to csi at the center of the north face

      do j = 1, ny-1
         do i = 2, nx-1

            np  = nx * (j-1) + i
            npw = np - 1

            xk(np) = x(np) - x(npw)
            yk(np) = y(np) - y(npw)

            gamman(np) = xk(np)**2 + yk(np)**2

         end do
      end do

      ! Derivatives relatively to csi at the center of the east face (only inner real faces)

      do j = 2, ny-1
         do i = 2, nx-2

            np   = nx * (j-1) + i
            npe  = np + 1

            xke(np) = xp(npe) - xp(np)
            yke(np) = yp(npe) - yp(np)

         end do
      end do

      ! Derivatives relatively to csi at the center of the east face (only west boundary of the domain)

      i = 2
      do j = 2, ny-1

         np   = nx * (j-1) + i
         nps  = np - nx
         npw  = np - 1
         npsw = nps - 1

         fw = ( x(npw) + x(npsw) ) / 2.d0
         fe = ( x(np)  + x(nps)  ) / 2.d0

         xke(npw) = -3.d0 * fw + 4.d0 * xp(np) -fe

         fw = ( y(npw) + y(npsw) ) / 2.d0
         fe = ( y(np)  + y(nps)  ) / 2.d0

         yke(npw) = -3.d0 * fw + 4.d0 * yp(np) -fe

      end do

      ! Derivatives relatively to csi at the center of the east face (only east boundary of the domain)

      i = nx-1
      do j = 2, ny-1

         np   = nx * (j-1) + i
         nps  = np - nx
         npw  = np - 1
         npsw = nps - 1

         fw = ( x(npw) + x(npsw) ) / 2.d0
         fe = ( x(np)  + x(nps)  ) / 2.d0

         xke(np) = fw - 4.d0 * xp(np) + 3.d0 * fe

         fw = ( y(npw) + y(npsw) ) / 2.d0
         fe = ( y(np)  + y(nps)  ) / 2.d0

         yke(np) = fw - 4.d0 * yp(np) + 3.d0 * fe

      end do

      ! Derivatives relatively to eta at the center of the north face (only inner real faces)

      do j = 2, ny-2
         do i = 2, nx-1

            np  = nx * (j-1) + i
            npn = np + nx

            xen(np) = xp(npn) - xp(np)
            yen(np) = yp(npn) - yp(np)

         end do
      end do

      ! Derivatives relatively to eta at the center of the north face (only south boundary of the domain)

      j = 2

      do i = 2, nx-1

         np   = nx * (j-1) + i
         nps  = np - nx
         npw  = np - 1
         npsw = nps - 1

         fs = ( x(nps) + x(npsw) ) / 2.d0
         fn = ( x(np)  + x(npw)  ) / 2.d0

         xen(nps) = -3.d0 * fs + 4.d0 * xp(np) -fn

         fs = ( y(nps) + y(npsw) ) / 2.d0
         fn = ( y(np)  + y(npw)  ) / 2.d0

         yen(nps) = -3.d0 * fs + 4.d0 * yp(np) -fn

      end do

      ! Derivatives relatively to eta at the center of the north face (only north boundary of the domain)

      j = ny - 1
      do i = 2, nx-1

         np   = nx * (j-1) + i
         nps  = np - nx
         npw  = np - 1
         npsw = nps - 1

         fs = ( x(nps) + x(npsw) ) / 2.d0
         fn = ( x(np)  + x(npw)  ) / 2.d0

         xen(np) = fs - 4.d0 * xp(np) + 3.d0 * fn

         fs = ( y(nps) + y(npsw) ) / 2.d0
         fn = ( y(np)  + y(npw)  ) / 2.d0

         yen(np) = fs - 4.d0 * yp(np) + 3.d0 * fn

      end do

      ! Beta and J at the center of the east face (all real faces)

      do j = 2, ny-1
         do i = 1, nx-1

            np   = nx * (j-1) + i

            betae(np) = xke(np) * xe(np) + yke(np) * ye(np)

            Je(np) = 1.d0 / ( xke(np) * ye(np) - xe(np) * yke(np) )

         end do
      end do

      ! Beta and J at center of the north face (all real faces)

      do j = 1, ny-1
         do i = 2, nx-1

            np   = nx * (j-1) + i

            betan(np) = xk(np) * xen(np) + yk(np) * yen(np)

            Jn(np) = 1.d0 / ( xk(np) * yen(np) - xen(np) * yk(np) )

         end do
      end do

      ! Jacobian J at the center of all real volumes

      do j = 2, ny-1
         do i = 2, nx-1

            np   = nx * (j-1) + i
            nps  = np - nx
            npw  = np - 1
            npsw = nps - 1

            fw = ( x(npw) + x(npsw) ) / 2.d0
            fe = ( x(np)  + x(nps)  ) / 2.d0
            fs = ( x(nps) + x(npsw) ) / 2.d0
            fn = ( x(np)  + x(npw)  ) / 2.d0

            xkp = fe - fw ! (x_csi)_P
            xep = fn - fs ! (x_eta)_P

            fw = ( y(npw) + y(npsw) ) / 2.d0
            fe = ( y(np)  + y(nps)  ) / 2.d0
            fs = ( y(nps) + y(npsw) ) / 2.d0
            fn = ( y(np)  + y(npw)  ) / 2.d0

            ykp = fe - fw ! (y_csi)_P
            yep = fn - fs ! (y_eta)_P

            Jp(np) = 1.d0 / ( xkp * yep - xep * ykp )

         end do
      end do

      Je = dabs(Je)
      Jn = dabs(Jn)
      Jp = dabs(Jp)

   end subroutine get_metrics


   !> \brief Calculates the radius for axisymmetric flows
   subroutine get_radius(coord, nx, ny, y, yp & ! Input
      ,                radius, re, rn, rp   ) ! Output
      implicit none
      integer, intent(in)  :: coord   ! Coordinate system ( 1=cylindrical, 0 = cartesian)
      integer, intent(in)  :: nx, ny  ! Number of volumes in csi and eta directions (real+fictitious)
      real(8), dimension (nx*ny), intent(in)  :: y  ! Coord. y of the northeast corner of volume P
      real(8), dimension (nx*ny), intent(in)  :: yp ! Coord. y of the centroid of volume P
      real(8), dimension (nx*ny), intent(out) :: radius ! Radius of northest corner of volume P
      real(8), dimension (nx*ny), intent(out) :: re ! Radius of the center of east face of volume P
      real(8), dimension (nx*ny), intent(out) :: rn ! Radius of the center of north face of volume P
      real(8), dimension (nx*ny), intent(out) :: rp ! Radius of the center of volume P

      integer :: i, j, np, npe, npw, npn, nps, npsw, npse, npnw, npne

      ! Cartesian coordinate system

      radius = 1.d0
      re = 1.d0
      rn = 1.d0
      rp = 1.d0

      ! Cylindrical coordinate system

      if ( coord == 1 ) then

         ! Radius of northeast corner

         radius = y

         ! Radius at the center of the east face

         do i = 1, nx-1
            do j = 2, ny-1

               np  = nx * (j-1) + i
               nps = np - nx

               re(np) = ( y(np) + y(nps) ) / 2.d0

            end do
         end do

         ! Radius at the center of the north face

         do i = 2, nx-1
            do j = 1, ny-1

               np  = nx * (j-1) + i
               npw = np - 1

               rn(np) = ( y(np) + y(npw) ) / 2.d0

            end do
         end do

         ! Radius of the center of the volume

         rp = yp

         ! Radius  of the center of the volume (south fictitious)

         j = 1
         do i = 2, nx-1

            np  = nx * (j-1) + i
            npn = np + nx

            rp(np) = rn(np) - rp(npn) + rn(np)
            !rp(np) = rp(npn) - rp(npn+nx) + rp(npn)

         end do

         ! Radius  of the center of the volume (north fictitious)

         j = ny
         do i = 2, nx-1

            np  = nx * (j-1) + i
            nps = np - nx

            rp(np) = rn(nps) - rp(nps) + rn(nps)
            !rp(np) = rp(nps) - rp(nps-nx) + rp(nps)

         end do

         ! Radius  of the center of the volume (west fictitious)

         i = 1
         do j = 2, ny-1

            np  = nx * (j-1) + i
            npe = np + 1

            rp(np) = re(np) - rp(npe) + re(np)
            !rp(np) = rp(npe) - rp(npe+1) + rp(npe)

         end do

         ! Radius  of the center of the volume (east fictitious)

         i = nx
         do j = 2, ny-1

            np  = nx * (j-1) + i
            npw = np - 1

            rp(np) = re(npw) - rp(npw) + re(npw)
            !rp(np) = rp(npw) - rp(npw-1) + rp(npw)

         end do

         ! Radius  of the center of the volume (corner fictitious: SW)

         i = 1 ; j = 1

         np  = nx * (j-1) + i
         npn = np + nx
         npe  = np + 1
         npne = npn + 1

         !rp(np) = rp(npn) - rp(npn+nx) + rp(npn)
         rp(np) = 4.d0 * y(np) - rp(npe) - rp(npn) - rp(npne)

         ! Radius  of the center of the volume (corner fictitious: SE)

         i = nx; j = 1

         np  = nx * (j-1) + i
         npn = np + nx
         npw  = np - 1
         npnw = npn - 1

         !rp(np) = rp(npn) - rp(npn+nx) + rp(npn)
         rp(np) = 4.d0 * y(npw) - rp(npw) - rp(npn) - rp(npnw)

         ! Radius  of the center of the volume (corner fictitious: NW)

         i = 1; j = ny

         np  = nx * (j-1) + i
         nps = np - nx
         npe  = np + 1
         npse = nps + 1

         !rp(np) = rp(nps) - rp(nps-nx) + rp(nps)
         rp(np) = 4.d0 * y(nps) - rp(npe) - rp(nps) - rp(npse)

         ! Radius  of the center of the volume (corner fictitious: NE)

         i = nx; j = ny

         np  = nx * (j-1) + i
         nps = np - nx
         npw  = np - 1
         npsw = nps - 1

         !rp(np) = rp(nps) - rp(nps-nx) + rp(nps)
         rp(np) = 4.d0 * y(npsw) - rp(npw) - rp(nps) - rp(npsw)

      end if

   end subroutine get_radius


end module
