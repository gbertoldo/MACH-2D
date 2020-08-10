!>
!! \brief Thermophysical module provides an abstract interface to a thermophysical
!!        model. It also provides subroutines for calculation of thermophysical
!!        properties on the calculation domain.
!!
module thermophysical

   use mod_class_ifile
   use mod_class_thermophysical_abstract
   use mod_class_thermophysical_constant
   use mod_class_thermophysical_mixture

   implicit none

   ! Thermophysical properties are constant or temperature dependent:
   integer, parameter, public :: THERMOPHYSICAL_CONSTANT = 0
   integer, parameter, public :: THERMOPHYSICAL_VARIABLE = 1

   class(class_thermophysical_abstract), pointer, private :: thermomodel !< A pointer to the thermophysical model

contains


   !> \brief Initializes variables related to the thermophysical properties of the gas
   subroutine thermophysical_initialization(ifile, tm, Rg, ktm) ! Output: last three
      implicit none
      class(class_ifile),                            intent(in)  :: ifile !< Input file
      class(class_thermophysical_abstract), pointer, intent(out) :: tm    !< A pointer to the thermophysical model
      real(8),                                       intent(out) :: Rg    !< Gas constant (J/kg.K)
      integer,                                       intent(out) :: ktm   !< Kind of thermophysical model: constant or variable

      ! Inner variables
      character(len=500) :: caux

      ! Reading the thermophysical model option
      call ifile%get_value(caux, "thermo_model")

      ! Thermophysical properties are temperature dependent, unless otherwise stated
      ktm = THERMOPHYSICAL_VARIABLE

      if ( trim(caux) == "CONSTANT" ) then

         ktm = THERMOPHYSICAL_CONSTANT

         allocate(class_thermophysical_constant::thermomodel)

      else if (trim(caux) == "MIXTURE" ) then

         allocate(class_thermophysical_mixture::thermomodel)

      else

         write(*,*) "thermophysical_initialization"
         write(*,*) "Unknown opton: ", trim(caux)
         write(*,*) "Stopping..."
         stop

      end if


      ! Initializing the object
      select type (thermomodel)

         type is ( class_thermophysical_mixture )

            call thermomodel%init(ifile)

         type is ( class_thermophysical_constant )

            call thermomodel%init(ifile)

      end select

      tm => thermomodel

      Rg = thermomodel%Rg

   end subroutine


   !> \brief Calculates cp at the center of each real volume and over the
   !! boundaries. The boundary values are stored in the fictitious volumes.
   subroutine set_cp( nx, ny, Tbn, Tbs, Tbe, Tbw, T, cp) ! Output: last one
      implicit none
      integer,                   intent(in)  :: nx  !< Number of volumes in the csi direction (real+fictitious)
      integer,                   intent(in)  :: ny  !< Number of volumes in the eta direction (real+fictitious)
      real(8), dimension(nx),    intent(in)  :: Tbn !< Temperature over the north boundary (K)
      real(8), dimension(nx),    intent(in)  :: Tbs !< Temperature over the south boundary (K)
      real(8), dimension(ny),    intent(in)  :: Tbe !< Temperature over the  east boundary (K)
      real(8), dimension(ny),    intent(in)  :: Tbw !< Temperature over the  west boundary (K)
      real(8), dimension(nx*ny), intent(in)  :: T   !< Temperature at center of volume P (K)
      real(8), dimension(nx*ny), intent(out) :: cp  !< Specific heat at const pressure (J/(kg.K))

      ! Inner variables

      integer :: i  ! Dummy index
      integer :: j  ! Dummy index
      integer :: np, npn, nps, npe, npw, npsw, npse, npnw, npne ! Dummy index

      ! Real volumes

      do j = 2, ny-1

         do i = 2, nx-1

             np  = nx * (j-1) + i

             cp(np) = thermomodel%cp(T(np))

         end do

      end do



      ! North boundary

      j = ny

      do i = 2, nx-1

         np  = nx * (j-1) + i

         cp(np) = thermomodel%cp(Tbn(i))

      end do



      ! South boundary

      j = 1

      do i = 2, nx-1

         np  = nx * (j-1) + i

         cp(np) = thermomodel%cp( Tbs(i) )

      end do



      ! East boundary

      i = nx

      do j = 2, ny-1

         np  = nx * (j-1) + i

         cp(np) = thermomodel%cp( Tbe(j) )

      end do



      ! West boundary

      i = 1

      do j = 2, ny-1

         np  = nx * (j-1) + i

         cp(np) = thermomodel%cp( Tbw(j) )

      end do


      ! SW corner

      i = 1;    j = 1

      np  = nx * (j-1) + i
      npn  = np + nx
      npe  = np + 1
      npne = npn + 1

      cp(np) = ( cp(npn) + cp(npe) + cp(npne) ) / 3.d0


      ! SE corner

      i = nx;   j = 1

      np   = nx * (j-1) + i
      npn  = np + nx
      npw  = np - 1
      npnw = npn - 1

      cp(np) = ( cp(npn) + cp(npw) + cp(npnw) ) / 3.d0



      ! NW corner

      i = 1;  j = ny

      np   = nx * (j-1) + i
      nps  = np - nx
      npe  = np + 1
      npse = nps + 1

      cp(np) = ( cp(nps) + cp(npe) + cp(npse) ) / 3.d0



      ! NE corner

      i = nx;   j = ny

      np   = nx * (j-1) + i
      nps  = np - nx
      npw  = np - 1
      npsw = nps - 1

      cp(np) = ( cp(nps) + cp(npw) + cp(npsw) ) / 3.d0

   end subroutine set_cp

   !> \brief Calculates gamma = Cp / Cv at the center of each real volume and
   !! over the boundaries based on Cp and Rg. The boundary values are stored in
   !! the fictitious volumes.
   subroutine set_gamma( nx, ny, Rg, cp, gcp) ! Output: last one
      implicit none
      integer,                   intent(in)  :: nx  !< Number of volumes in the csi direction (real+fictitious)
      integer,                   intent(in)  :: ny  !< Number of volumes in the eta direction (real+fictitious)
      real(8),                   intent(in)  :: Rg  !< Perfect gas constant (J/(kg.K))
      real(8), dimension(nx*ny), intent(in)  :: cp  !< Specific heat at const pressure (J/(kg.K))
      real(8), dimension(nx*ny), intent(out) :: gcp !< gcp = gamma = Cp/Cv at center of CV P (dimensionless)

      gcp = cp / ( cp - Rg ) ! dimensionless

   end subroutine set_gamma


   !> \brief Calculates the laminar viscosity at the nodes of real volumes
   subroutine set_laminar_viscosity_at_nodes(nx, ny, T, vlp) ! Output: last one
      implicit none
      integer,                   intent(in)  :: nx  !< Number of volumes in csi direction (real+fictitious)
      integer,                   intent(in)  :: ny  !< Number of volumes in eta direction (real+fictitious)
      real(8), dimension(nx*ny), intent(in)  :: T   !< Temperature of the last iteraction (K)
      real(8), dimension(nx*ny), intent(out) :: vlp !< Laminar viscosity at center of volume P (Pa.s)


      ! Inner variables

      integer :: i, j, np

      do j = 2, ny-1

         do i = 2, nx-1

            np = nx * (j-1) + i

            vlp(np) = thermomodel%mu( T(np) )

         end do

      end do

   end subroutine set_laminar_viscosity_at_nodes


   !> \brief Calculates the thermal conductivity at the nodes of real volumes
   subroutine set_thermal_conductivity_at_nodes(nx, ny, T, kp) ! Output: last one
      implicit none
      integer,                   intent(in)  :: nx !< Number of volumes in csi direction (real+fictitious)
      integer,                   intent(in)  :: ny !< Number of volumes in eta direction (real+fictitious)
      real(8), dimension(nx*ny), intent(in)  :: T  !< Temperature of the last iteraction (K)
      real(8), dimension(nx*ny), intent(out) :: kp !< Thermal conductivity at center of volume P (W/m.K)

      ! Inner variables

      integer :: i, j, np


      do j = 2, ny-1

         do i = 2, nx-1

           np = nx * (j-1) + i

           kp(np) = thermomodel%kp( T(np) )

         end do

      end do

   end subroutine set_thermal_conductivity_at_nodes


   !> \brief Calculates the laminar viscosity at faces
   subroutine get_laminar_viscosity_at_faces(nx, ny, Tbn, Tbs, Tbe, Tbw, vlp, vle, vln) ! Output: last two
      implicit none
      integer, intent(in) :: nx         !< Number of volumes in csi direction (real+fictitious)
      integer, intent(in) :: ny         !< Number of volumes in eta direction (real+fictitious)
      real(8), dimension(nx),    intent(in)  :: Tbn !< Temperature over the north boundary (K)
      real(8), dimension(nx),    intent(in)  :: Tbs !< Temperature over the south boundary (K)
      real(8), dimension(ny),    intent(in)  :: Tbe !< Temperature over the  east boundary (K)
      real(8), dimension(ny),    intent(in)  :: Tbw !< Temperature over the  west boundary (K)
      real(8), dimension(nx*ny), intent(in)  :: vlp !< Laminar viscosity at center of volume P (Pa.s)
      real(8), dimension(nx*ny), intent(out) :: vle !< Laminar viscosity at center of face east (Pa.s)
      real(8), dimension(nx*ny), intent(out) :: vln !< Laminar viscosity at center of face north (Pa.s)


      ! Inner variables
      integer :: i, j, np, npe, npn


      ! Calculating the viscosity over the east faces of real volumes (excluding boundaries)

      do j = 2, ny-1

         do i = 2, nx-2

            np   = nx * (j-1) + i

            npe  = np + 1

            vle(np) = 2.d0 * vlp(np) * vlp(npe) / ( vlp(np) + vlp(npe) )

         end do

      end do

      ! Calculating the viscosity over the west boundary

      i = 1

      do j = 2, ny-1

         np = nx * (j-1) + i

         vle(np) = thermomodel%mu( Tbw(j) )

      end do


      ! Calculating the viscosity over the east boundary

      i = nx-1

      do j = 2, ny-1

         np = nx * (j-1) + i

         vle(np) = thermomodel%mu( Tbe(j) )

      end do


      ! Calculating the viscosity over the north faces of real volumes (excluding boundaries)

      do i = 2, nx-1

         do j = 2, ny-2

            np  = nx * (j-1) + i

            npn = np + nx

            vln(np) = 2.d0 * vlp(np) * vlp(npn) / ( vlp(np) + vlp(npn) )

         end do

      end do

      ! Calculating the viscosity over the south boundary

      j = 1

      do i = 2, nx-1

         np = nx * (j-1) + i

         vln(np) = thermomodel%mu( Tbs(i) )

      end do


      ! Calculating the viscosity over the north boundary

      j = ny-1

      do i = 2, nx-1

         np = nx * (j-1) + i

         vln(np) = thermomodel%mu( Tbn(i) )

      end do


   end subroutine get_laminar_viscosity_at_faces


   !> \brief Calculates the thermal conductivity at faces
   subroutine get_thermal_conductivity_at_faces(nx, ny, Tbn, Tbs, Tbe, Tbw, kp, ke, kn) ! Output: last two
      implicit none
      integer,                   intent(in)  :: nx  !< Number of volumes in csi direction (real+fictitious)
      integer,                   intent(in)  :: ny  !< Number of volumes in eta direction (real+fictitious)
      real(8), dimension(nx),    intent(in)  :: Tbn !< Temperature over the north boundary (K)
      real(8), dimension(nx),    intent(in)  :: Tbs !< Temperature over the south boundary (K)
      real(8), dimension(ny),    intent(in)  :: Tbe !< Temperature over the  east boundary (K)
      real(8), dimension(ny),    intent(in)  :: Tbw !< Temperature over the  west boundary (K)
      real(8), dimension(nx*ny), intent(in)  :: kp  !< Thermal conductivity at center of volume P (W/m.K)
      real(8), dimension(nx*ny), intent(out) :: ke  !< Thermal conductivity at center of face east (W/m.K)
      real(8), dimension(nx*ny), intent(out) :: kn  !< Thermal conductivity at center of face north (W/m.K)

      ! Inner variables
      integer :: i, j, np, npe, npn


      ! Calculating the thermal conductivity over the east faces of real volumes (excluding boundaries)

      do j = 2, ny-1

         do i = 2, nx-2

            np   = nx * (j-1) + i

            npe  = np + 1

            ke(np) = 2.d0 * kp(np) * kp(npe) / ( kp(np) + kp(npe) )

         end do

      end do


      ! Calculating the thermal conductivity over the west boundary

      i = 1

      do j = 2, ny-1

         np = nx * (j-1) + i

         ke(np) = thermomodel%kp( Tbw(j) )

      end do


      ! Calculating the thermal conductivity over the east boundary

      i = nx-1

      do j = 2, ny-1

         np = nx * (j-1) + i

         ke(np) = thermomodel%kp( Tbe(j) )

      end do


      ! Calculating the thermal conductivity over the north faces of real volumes (excluding boundaries)

      do i = 2, nx-1

         do j = 2, ny-2

            np  = nx * (j-1) + i

            npn = np + nx

            kn(np) = 2.d0 * kp(np) * kp(npn) / ( kp(np) + kp(npn) )

         end do

      end do


      ! Calculating the thermal conductivity over the south boundary

      j = 1

      do i = 2, nx-1

         np = nx * (j-1) + i

         kn(np) = thermomodel%kp( Tbs(i) )

      end do


      ! Calculating the thermal conductivity over the north boundary

      j = ny-1

      do i = 2, nx-1

         np = nx * (j-1) + i

         kn(np) = thermomodel%kp( Tbn(i) )

      end do

   end subroutine get_thermal_conductivity_at_faces

end module
