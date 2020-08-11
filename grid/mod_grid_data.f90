!>
!! \brief mod_grid_data contains data and type definitions related to the grid
!!
module mod_grid_data

   implicit none

   integer :: nx    !< number of volumes in the csi direction (real+fictitious)
   integer :: ny    !< number of volumes in the eta direction (real+fictitious)
   integer :: nxi   !< number of volumes in the csi direction (real+fictitious) of the coarsest grid
   integer :: nyi   !< number of volumes in the eta direction (real+fictitious) of the coarsest grid
   integer :: nxf   !< number of volumes in the csi direction (real+fictitious) of the finest grid
   integer :: nyf   !< number of volumes in the eta direction (real+fictitious) of the finest grid
   integer :: nxy   !< nxy = nx * ny
   integer :: nmf   !< number of the finest mesh
   integer :: nmd   !< number of the desired mesh
   integer :: kg    !< Kind of grid (1=uniform, 2=geometric progression, 3=power law, 4=gp modified)
   integer :: kcm   !< Kind of centroid mean (1=simple mean, 2=weighted mean)
   real(8) :: a1    !< width of the volume closer to the wall (m)
   real(8) :: avi   !< Initial value of the artificial viscosity
   real(8) :: avf   !< Final value of the artificial viscosity
   real(8) :: awf   !< Area weighting factor
   real(8) :: cbl   !< The width of the vol. closer to the wall is 'cbl' times the width of the b. layer
   real(8) :: wbl   !< Boundary layer estimated width (m)
   character(len=200) :: fgeom !< File of the geometric parameters

   real(8), allocatable, dimension(:) :: xf  !< coord. at the northeast corner of the volume P of the finest grid
   real(8), allocatable, dimension(:) :: yf  !< coord. at the northeast corner of the volume P of the finest grid

   real(8) :: lr   ! length of the rocket
   real(8) :: rb   ! base radius of the rocket

end module
