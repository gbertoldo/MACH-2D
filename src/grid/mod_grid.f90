!>
!! \brief mod_grid also gives an interface to the main program to calculate the
!!        grid and related metrics. mod_grid uses procedures from mod_grid_procedure
!!        to operate data of mod_grid_data.
!!
module mod_grid

   use mod_class_ifile
   use mod_grid_data
   use mod_grid_procedures
   use mod_extflow

   implicit none

   ! Makes everything private, except otherwise stated
   private

   ! Public methods
   public :: grid_init   &
      ,      grid_size   &
      ,      grid_create
contains

   !> \brief Initializes data of module mod_grid_data
   subroutine grid_init(ifile)
      implicit none
      class(class_ifile),  intent(in) :: ifile !< Input file

      call ifile%get_value(      nxi,    "nxi-2") ! Number of real volumes in the csi direction of the coarsest grid
      call ifile%get_value(      nyi,    "nyi-2") ! Number of real volumes in the eta direction of the coarsest grid
      call ifile%get_value(      nmf,      "nmf") ! Number of the finest mesh  (1<=nmf)
      call ifile%get_value(      nmd,      "nmd") ! Number of the desired mesh (1<=nmd<=nmf)
      call ifile%get_value(       kg,       "kg") ! Kind of grid (1=uniform, 2=geometric progression, 3=power law, 4=gp modified, 5=hyperbolic)
      call ifile%get_value(      avi,      "avi") ! Initial value of the artificial viscosity (only for kg=5)
      call ifile%get_value(      avf,      "avf") ! Final value of the artificial viscosity (only for kg=5)
      call ifile%get_value(      awf,      "awf") ! Area weighting factor (only for kg=5)
      call ifile%get_value(      kcm,      "kcm") ! Kind of centroid mean (1=simple mean, 2=weighted mean)
      call ifile%get_value(      cbl,      "cbl") ! The width of the vol. closer to the wall is 'cbl' times the width of the b. layer

      ! After reading the values of nxi and nyi, these variables must be
      ! changed to take into account the fictitious volumes
      nxi = nxi + 2
      nyi = nyi + 2

      ! Calculating the number of volumes (real+fictitious) of the finest mesh
      nxf = (nxi-2) * 2 ** (nmf-1) + 2
      nyf = (nyi-2) * 2 ** (nmf-1) + 2

      ! Calculating the number of volumes (real+fictitious) of the desired mesh
      nx = (nxi-2) * 2 ** (nmd-1) + 2
      ny = (nyi-2) * 2 ** (nmd-1) + 2

      allocate( xf(nxf*nyf), yf(nxf*nyf) )

   end subroutine


   !> \brief Returns the size of the grid
   subroutine grid_size(nx_out, ny_out)
      implicit none
      integer, intent(out) :: nx_out
      integer, intent(out) :: ny_out

      nx_out = nx
      ny_out = ny

   end subroutine


   !> \brief Creates the grid and related metrics
   subroutine grid_create(unt, coord                                  & ! Input
         ,  x, y, xp, yp, xe, ye, xen, yen, xk, yk, xke, yke, Jp      & ! Output
         ,  Je, Jn, alphae, gamman, betae, betan, radius, re, rn, rp  ) ! Output
      implicit none
      integer,               intent(in)  :: unt    !< Unit where the input parameters will be printed
      integer,               intent(in)  :: coord  !< Coordinate system ( 1=cylindrical, 0 = cartesian)
      real(8), dimension(:), intent(out) :: x      !<
      real(8), dimension(:), intent(out) :: y      !<
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

      ! Getting the finest grid

      ! Generates the grid north and south boundary
      call extflow_grid_boundary(nxf, nyf, unt, xf, yf) ! Output: last four

      ! Estimates the boundary layer width and the width of the volume
      ! closer to the wall
      call extflow_boundary_layer(cbl, wbl, a1) ! Output: last two

      ! Generates the grid according to kg option
      call set_grid(kg, nxf, nyf, a1, avi, avf, awf, xf, yf) ! Last 2 are inoutput

      ! Selecting the desided grid
      call get_coarser_grid(nxf, nyf, nx, ny, nmf, nmd, xf, yf, x, y)

      ! Removing the (now) unnecessary finest grid
      deallocate(xf, yf)

      ! Calculates the centroids of each real volume
      call get_real_centroids_xy( kcm, nx, ny, x, y, xp, yp ) ! Output: last two
      ! checked

      ! Calculates the components of the metric tensor and the Jacobian
      call get_metrics(nx, ny, x, y, xp, yp             & ! Input
      ,          xe, ye, xen, yen, xk, yk, xke, yke, Jp & ! Output
      ,          Je, Jn, alphae, gamman, betae, betan )   ! Output
      ! checked

      ! All the radii are equal 1 if planar flow is choose. If axisymmetric flow, radii are calculated
      call get_radius(coord, nx, ny, y, yp  & ! Input
      ,               radius, re, rn, rp    ) ! Output
      ! checked

   end subroutine

end module
