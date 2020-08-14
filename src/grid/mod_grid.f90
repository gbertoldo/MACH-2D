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
      type(type_intflow),    intent(out) :: iflow  !< Variables related to internal flow

! TODO (guilherme#1#): Develop a better way to create the grid that may be compatible to both internal and external flow

      real(8) :: wbl   ! Boundary layer estimated width (m)
      real(8) :: a1    ! width of the volume closer to the wall (m)
      real(8), allocatable, dimension(:) :: xf  !< coord. at the northeast corner of the volume P of the finest grid
      real(8), allocatable, dimension(:) :: yf  !< coord. at the northeast corner of the volume P of the finest grid

      ! Allocating the finest grid
      allocate( xf(gridpar%nxf*gridpar%nyf) )
      allocate( yf(gridpar%nxf*gridpar%nyf) )

      ! Getting the grid boundary for the finest grid
      if ( kflow == 0 ) then ! External flow

         ! Generates the grid north and south boundary
         call extflow_grid_boundary(gridpar%nxf, gridpar%nyf, unt, xf, yf) ! Output: last four

         ! Estimates the boundary layer width and the width of the volume
         ! closer to the wall
         call extflow_boundary_layer(gridpar%cbl, wbl, a1) ! Output: last two

      else ! Internal flow

         call intflow_grid_boundary(gridpar%nxf, gridpar%nyf, xf, yf, iflow, a1) ! Output: last four

      end if

      ! Generates the grid according to kg option
      call set_grid(gridpar%kg, gridpar%nxf, gridpar%nyf, a1, gridpar%avi, gridpar%avf, gridpar%awf, xf, yf) ! Last 2 are inoutput

      ! Selecting the desided grid
      call get_coarser_grid(gridpar%nxf, gridpar%nyf, gridpar%nx, gridpar%ny, gridpar%nmf, gridpar%nmd, xf, yf, x, y)

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

end module
