!>
!! \brief mod_extflow provides an interface to the main program to calculate
!!        variables and boundary conditions of the external flow. This module
!!        uses procedures from mod_extflow_procedures and data from
!!        mod_extflow_data.
!!
module mod_extflow

   use mod_class_ifile
   use mod_extflow_data
   use mod_extflow_geometry
   use mod_extflow_procedures                                                         &
      ,    extflow_boundary_simplec              => get_boundary_simplec_coefficients &
      ,    extflow_velocities_at_boundary_faces  => get_velocities_at_boundary_faces  &
      ,    extflow_set_bcpl                      => get_bc_scheme_pl                  &
      ,    extflow_p_extrapolation_to_fictitious => get_p_extrapolation_to_fictitious &
      ,    extflow_get_Ucbe_Vcbe                 => get_Ucbe_Vcbe

   implicit none

   ! Makes everything private, except otherwise stated
   private

   ! Public procedures
   public :: extflow_init                          &
      ,      extflow_grid_boundary                 &
      ,      extflow_boundary_layer                &
      ,      extflow_initial_conditions            &
      ,      extflow_set_bcu                       &
      ,      extflow_set_bcv                       &
      ,      extflow_set_bcpl                      &
      ,      extflow_set_bcro                      &
      ,      extflow_set_bcT                       &
      ,      extflow_boundary_simplec              &
      ,      extflow_velocities_at_boundary_faces  &
      ,      extflow_p_extrapolation_to_fictitious &
      ,      extflow_calc_main_variables           &
      ,      extflow_get_Ucbe_Vcbe

contains

   !> \brief Initializes external flow module
   subroutine extflow_init(ifile, thermomodel, Tref, Href, Cpref) ! Output: last three
      implicit none
      class(class_ifile),                            intent(in)  :: ifile       !< Input file
      class(class_thermophysical_abstract), pointer, intent(in)  :: thermomodel !< A pointer to the thermophysical model
      real(8),                                       intent(out) :: Tref        !< Reference temperature (K)
      real(8),                                       intent(out) :: Href        !< Reference total enthalpy (m2/s2)
      real(8),                                       intent(out) :: Cpref       !< Reference Cp (J/kg.K)



      call ifile%get_value(      kfc,      "kfc") ! Kind of foredrag calculation ( 0 = over the whole forebody; 1 = over the ogive only)
      call ifile%get_value(     Tsbc,     "Tsbc") ! Temperature on the south boundary (K) (if negative, adiabatic bc is applied)
      call ifile%get_value(       PF,       "PF") ! Far field pressure (Pa)
      call ifile%get_value(       TF,       "TF") ! Far field temperature (K)
      call ifile%get_value(       MF,       "MF") ! Mach number of the free stream
      call ifile%get_value(    fgeom,    "fgeom") ! File of boundary geometry

      ! Calculates some of the free stream properties of the gas
      call get_free_stream_properties(thermomodel, PF, TF, MF, CPF, GF, VLF, KPF, PRF, ROF, UF, REFm, HF) ! Output: last 9

      ! Calculating reference values
      Tref  = TF
      Href  = HF
      Cpref = thermomodel%cp(TF)

      ! fgeom, lr and rb are initialized with the subroutine extflow_grid_boundary

   end subroutine


   !> \brief Defines the boundary of the domain
   subroutine extflow_grid_boundary(nx, ny, unt, x, y) ! Output: last two
      implicit none
      integer,                   intent(in)  ::  nx !< Number of volumes in the csi direction (real+fictitious)
      integer,                   intent(in)  ::  ny !< Number of volumes in the eta direction (real+fictitious)
      integer,                   intent(in)  :: unt !< Unit where the input parameters will be printed
      real(8), dimension(nx*ny), intent(out) ::   x !< Coord. x at the northeast corner of the volume P (m)
      real(8), dimension(nx*ny), intent(out) ::   y !< Coord. y at the northeast corner of the volume P (m)

      ! Generates the grid north and south boundary
      call get_grid_boundary(nx, ny, unt, fgeom, lr, rb, x, y) ! Output: last four

   end subroutine


   !> \brief Estimates the boundary layer width and the width of the volume
   !! closer to the wall
   subroutine extflow_boundary_layer(cbl, wbl, a1) ! Output: last two
      implicit none
      real(8), intent(in)  ::  cbl !< The width of the vol. closer to the wall is 'cbl' times the width of the b. layer
      real(8), intent(out) ::  wbl !< Estimated width of the boundary layer (m)
      real(8), intent(out) ::   a1 !< Width of the volume closer to the wall (m)

      call get_boundary_layer_width_estimate(lr, REFm, cbl, wbl, a1) ! Output: last two

   end subroutine


   !> \brief Sets the initial conditions for external flow
   subroutine extflow_initial_conditions(nx, ny, itemax, modvis, Rg, xe, ye & ! Input
      ,                                  xk, yk, xke, yke, alphae, betae    & ! Input
      ,                            p, T, ro, roe, ron, u, v, ue, ve, un, vn & ! Output
      ,                            Uce, Vcn, Ucbe, Vcbe, Tbn, Tbs, Tbe, Tbw ) ! Output
      implicit none
      integer, intent(in)  :: nx     !< Number of volumes in csi direction (real+fictitious)
      integer, intent(in)  :: ny     !< Number of volumes in eta direction (real+fictitious)
      integer, intent(in)  :: itemax !< Number of iteractions for extrapolation to fictitious
      integer, intent(in)  :: modvis !< Viscosity model (0=Euler, 1=NS)
      real(8), intent(in)  :: Rg     !< Perfect gas constant (J/kg.K)
      real(8), dimension(nx*ny), intent(in)  :: xe     !< x_eta at face east of volume P (m)
      real(8), dimension(nx*ny), intent(in)  :: ye     !< y_eta at face east of volume P (m)
      real(8), dimension(nx*ny), intent(in)  :: xk     !< x_csi at face north of volume P (m)
      real(8), dimension(nx*ny), intent(in)  :: yk     !< y_csi at face north of volume P (m)
      real(8), dimension(nx*ny), intent(in)  :: xke    !< x_csi at center of face east
      real(8), dimension(nx*ny), intent(in)  :: yke    !< y_csi at center of face east
      real(8), dimension(nx*ny), intent(in)  :: alphae !< (metric) alpha at the center of east face of volume P
      real(8), dimension(nx*ny), intent(in)  :: betae  !< (metric) beta  at the center of east face of volume P
      real(8), dimension(nx*ny), intent(out) :: p      !< Pressure at center o volume P (Pa)
      real(8), dimension(nx*ny), intent(out) :: T      !< Temperature of the last iteraction (K)
      real(8), dimension(nx*ny), intent(out) :: ro     !< Specific mass (absolute density) at center of volumes (kg/m3)
      real(8), dimension(nx*ny), intent(out) :: roe    !< Absolute density at east face (kg/m3)
      real(8), dimension(nx*ny), intent(out) :: ron    !< Absolute density at north face (kg/m3)
      real(8), dimension(nx*ny), intent(out) :: u      !< Cartesian velocity of the last iteraction (m/s)
      real(8), dimension(nx*ny), intent(out) :: v      !< Cartesian velocity of the last iteraction (m/s)
      real(8), dimension(nx*ny), intent(out) :: ue     !< Cartesian velocity u at center of east face (m/s)
      real(8), dimension(nx*ny), intent(out) :: ve     !< Cartesian velocity v at center of east face (m/s)
      real(8), dimension(nx*ny), intent(out) :: un     !< Cartesian velocity u at center of north face (m/s)
      real(8), dimension(nx*ny), intent(out) :: vn     !< Cartesian velocity v at center of north face (m/s)
      real(8), dimension(nx*ny), intent(out) :: Uce    !< Contravariant velocity U at east face (m2/s)
      real(8), dimension(nx*ny), intent(out) :: Vcn    !< Contravariant velocity V at north face (m2/s)
      real(8), dimension(ny),    intent(out) :: Ucbe   !< Uc over the faces of the east  boundary (m2/s)
      real(8), dimension(ny),    intent(out) :: Vcbe   !< Vc over the faces of the east  boundary (m2/s)
      real(8), dimension(nx),    intent(out) :: Tbn    !< Temperature over the north boundary (K)
      real(8), dimension(nx),    intent(out) :: Tbs    !< Temperature over the south boundary (K)
      real(8), dimension(ny),    intent(out) :: Tbe    !< Temperature over the  east boundary (K)
      real(8), dimension(ny),    intent(out) :: Tbw    !< Temperature over the  west boundary (K)

       call get_initial_conditions( nx, ny, itemax, modvis, PF, TF, Rg, GF, MF, UF &
      , xe, ye, xk, yk, xke, yke, alphae, betae, p, T, ro, roe, ron, u, v, ue  &
      , ve, un, vn, Uce, Vcn, Ucbe, Vcbe, Tbn, Tbs, Tbe, Tbw ) ! UF and last 19 are output
      ! checked

      ! Selects the portion of the forebody where the drag will be calculated
      select case ( kfc )

         case (0) ! Over the whole forebody

            iocs = nx-1

         case (1) ! Over the ogive only

            call get_ogive_cylinder_matching_point_south(nx, ny, xk, yk, iocs) ! Output: last one

         case default

            write(*,*) "kfc: unknown option. Stopping..."

            stop

      end select

   end subroutine


   !> \brief Sets boundary condition for u
   subroutine extflow_set_bcu(nx, ny, modvis, xk, yk, alphae, betae & ! Intput
      ,                                            u, v, Ucbe, Vcbe & ! Intput
      ,                                      a9bn, a9bs, a9be, a9bw & ! Output
      ,                                      b9bn, b9bs, b9be, b9bw ) ! Output
      implicit none
      integer,                   intent(in)  :: nx     !< Number of volumes in csi direction (real+fictitious)
      integer,                   intent(in)  :: ny     !< Number of volumes in eta direction (real+fictitious)
      integer,                   intent(in)  :: modvis !< Viscosity model (0=Euler, 1=NS)
      real(8), dimension(nx*ny), intent(in)  :: xk     !< x_csi at face north of volume P (m)
      real(8), dimension(nx*ny), intent(in)  :: yk     !< y_csi at face north of volume P (m)
      real(8), dimension(nx*ny), intent(in)  :: alphae !< (metric) alpha at the center of east face of volume P
      real(8), dimension(nx*ny), intent(in)  :: betae  !< (metric) beta  at the center of east face of volume P
      real(8), dimension(nx*ny), intent(in)  :: u      !< x cartesian velocity (m/s)
      real(8), dimension(nx*ny), intent(in)  :: v      !< y cartesian velocity (m/s)
      real(8), dimension(ny),    intent(in)  :: Ucbe   !< Uc over the faces of the east  boundary (m2/s)
      real(8), dimension(ny),    intent(in)  :: Vcbe   !< Vc over the faces of the east  boundary (m2/s)
      real(8), dimension(nx,9),  intent(out) :: a9bn   !< Coefficients of the discretization for the north boundary
      real(8), dimension(nx,9),  intent(out) :: a9bs   !< Coefficients of the discretization for the south boundary
      real(8), dimension(ny,9),  intent(out) :: a9be   !< Coefficients of the discretization for the east boundary
      real(8), dimension(ny,9),  intent(out) :: a9bw   !< Coefficients of the discretization for the west boundary
      real(8), dimension(nx),    intent(out) :: b9bn   !< Source of the discretization for the north boundary
      real(8), dimension(nx),    intent(out) :: b9bs   !< Source of the discretization for the south boundar
      real(8), dimension(ny),    intent(out) :: b9be   !< Source of the discretization for the east boundary
      real(8), dimension(ny),    intent(out) :: b9bw   !< Source of the discretization for the west boundary

      call get_bc_scheme_u(nx, ny, modvis, UF, xk, yk, alphae, betae, u, v &
         ,       Ucbe, Vcbe, a9bn, a9bs, a9be, a9bw, b9bn, b9bs, b9be, b9bw) ! Output: last eight

   end subroutine


   !> \brief Sets boundary condition for v
   subroutine extflow_set_bcv(nx, ny, modvis, xk, yk, u, v, Ucbe, Vcbe & ! Intput
      ,                                         a9bn, a9bs, a9be, a9bw & ! Output
      ,                                         b9bn, b9bs, b9be, b9bw ) ! Output
      implicit none
      integer,                   intent(in)  :: nx     !< Number of volumes in csi direction (real+fictitious)
      integer,                   intent(in)  :: ny     !< Number of volumes in eta direction (real+fictitious)
      integer,                   intent(in)  :: modvis !< Viscosity model (0=Euler, 1=NS)
      real(8), dimension(nx*ny), intent(in)  :: xk     !< x_csi at face north of volume P (m)
      real(8), dimension(nx*ny), intent(in)  :: yk     !< y_csi at face north of volume P (m)
      real(8), dimension(nx*ny), intent(in)  :: u      !< x cartesian velocity (m/s)
      real(8), dimension(nx*ny), intent(in)  :: v      !< y cartesian velocity (m/s)
      real(8), dimension(ny),    intent(in)  :: Ucbe   !< Uc over the faces of the east  boundary (m2/s)
      real(8), dimension(ny),    intent(in)  :: Vcbe   !< Vc over the faces of the east  boundary (m2/s)
      real(8), dimension(nx,9),  intent(out) :: a9bn   !< Coefficients of the discretization for the north boundary
      real(8), dimension(nx,9),  intent(out) :: a9bs   !< Coefficients of the discretization for the south boundary
      real(8), dimension(ny,9),  intent(out) :: a9be   !< Coefficients of the discretization for the east boundary
      real(8), dimension(ny,9),  intent(out) :: a9bw   !< Coefficients of the discretization for the west boundary
      real(8), dimension(nx),    intent(out) :: b9bn   !< Source of the discretization for the north boundary
      real(8), dimension(nx),    intent(out) :: b9bs   !< Source of the discretization for the south boundar
      real(8), dimension(ny),    intent(out) :: b9be   !< Source of the discretization for the east boundary
      real(8), dimension(ny),    intent(out) :: b9bw   !< Source of the discretization for the west boundary

      call get_bc_scheme_v(nx, ny, modvis, xk, yk, u, v, Ucbe, Vcbe &
            , a9bn, a9bs, a9be, a9bw, b9bn, b9bs, b9be, b9bw) ! Output: last eight

   end subroutine


   !> \brief Sets boundary condition for ro
   subroutine extflow_set_bcro(nx, ny, alphae, betae, betan, gamman, Ucbe, Vcbe & ! Input
         ,                       a9bn, a9bs, a9be, a9bw, b9bn, b9bs, b9be, b9bw ) ! Output
      implicit none
      integer, intent(in) :: nx  !< Number of volumes in csi direction (real+fictitious)
      integer, intent(in) :: ny  !< Number of volumes in eta direction (real+fictitious)
      real(8), dimension(nx*ny), intent(in) :: alphae !< (metric) alpha at the center of east face of volume P
      real(8), dimension(nx*ny), intent(in) :: betae  !< (metric) beta  at the center of east face of volume P
      real(8), dimension(nx*ny), intent(in) :: betan  !< (metric) beta  at the center of north face of volume P
      real(8), dimension(nx*ny), intent(in) :: gamman !< (metric) gamma at the center of north face of volume P
      real(8), dimension(ny),   intent(in)  :: Ucbe !< Uc over the faces of the east  boundary (m2/s)
      real(8), dimension(ny),   intent(in)  :: Vcbe !< Vc over the faces of the east  boundary (m2/s)
      real(8), dimension(nx,9), intent(out) :: a9bn !< Coefficients of the discretization for the north boundary
      real(8), dimension(nx,9), intent(out) :: a9bs !< Coefficients of the discretization for the south boundary
      real(8), dimension(ny,9), intent(out) :: a9be !< Coefficients of the discretization for the east boundary
      real(8), dimension(ny,9), intent(out) :: a9bw !< Coefficients of the discretization for the west boundary
      real(8), dimension(nx),   intent(out) :: b9bn !< Source of the discretization for the north boundary
      real(8), dimension(nx),   intent(out) :: b9bs !< Source of the discretization for the south boundar
      real(8), dimension(ny),   intent(out) :: b9be !< Source of the discretization for the east boundary
      real(8), dimension(ny),   intent(out) :: b9bw !< Source of the discretization for the west boundary

      ! Defines the numerical scheme for the boundary conditions of ro
      call get_bc_scheme_ro(nx, ny, ROF, alphae, betae, betan, gamman &
         , Ucbe, Vcbe, a9bn, a9bs, a9be, a9bw, b9bn, b9bs, b9be, b9bw) ! Output: last eight

   end subroutine


   !> \brief Sets boundary condition for T
   subroutine extflow_set_bcT(nx, ny, alphae, betae, betan, gamman, Ucbe, Vcbe & ! Input
         ,                      a9bn, a9bs, a9be, a9bw, b9bn, b9bs, b9be, b9bw ) ! Output
      implicit none
      integer, intent(in) :: nx   !< Number of volumes in csi direction (real+fictitious)
      integer, intent(in) :: ny   !< Number of volumes in eta direction (real+fictitious)
      real(8), dimension(nx*ny), intent(in) :: alphae !< (metric) alpha at the center of east face of volume P
      real(8), dimension(nx*ny), intent(in) :: betae  !< (metric) beta  at the center of east face of volume P
      real(8), dimension(nx*ny), intent(in) :: betan  !< (metric) beta  at the center of north face of volume P
      real(8), dimension(nx*ny), intent(in) :: gamman !< (metric) gamma at the center of north face of volume P
      real(8), dimension(ny),    intent(in) :: Ucbe !< Uc over the faces of the east  boundary (m2/s)
      real(8), dimension(ny),    intent(in) :: Vcbe !< Vc over the faces of the east  boundary (m2/s)
      real(8), dimension(nx,9), intent(out) :: a9bn !< Coefficients of the discretization for the north boundary
      real(8), dimension(nx,9), intent(out) :: a9bs !< Coefficients of the discretization for the south boundary
      real(8), dimension(ny,9), intent(out) :: a9be !< Coefficients of the discretization for the east boundary
      real(8), dimension(ny,9), intent(out) :: a9bw !< Coefficients of the discretization for the west boundary
      real(8), dimension(nx),   intent(out) :: b9bn !< Source of the discretization for the north boundary
      real(8), dimension(nx),   intent(out) :: b9bs !< Source of the discretization for the south boundar
      real(8), dimension(ny),   intent(out) :: b9be !< Source of the discretization for the east boundary
      real(8), dimension(ny),   intent(out) :: b9bw !< Source of the discretization for the west boundary

      call get_bc_scheme_T(nx, ny, TF, Tsbc, alphae, betae, betan, gamman &
            , Ucbe, Vcbe, a9bn, a9bs, a9be, a9bw, b9bn, b9bs, b9be, b9bw) ! Output: last eight

   end subroutine


   !> \brief Calculates main variables of interest related to the external flow
   subroutine extflow_calc_main_variables(nx, ny, Rg, modvis, xk, yk, rn, Jn, vln, p, u, v, Vcn, msg) ! Output: last one
      integer,                    intent(in)  ::     nx !< Number of volumes in csi direction (real+fictitious)
      integer,                    intent(in)  ::     ny !< Number of volumes in eta direction (real+fictitious)
      real(8),                    intent(in)  ::     Rg !< Perfect gas constant
      integer,                    intent(in)  :: modvis !< modvis = 0 -> Euler;  modvis = 1 -> Navier-Stokes
      real(8), dimension (nx*ny), intent(in)  ::     xk !< x_csi at center of face north
      real(8), dimension (nx*ny), intent(in)  ::     yk !< y_csi at center of face north
      real(8), dimension (nx*ny), intent(in)  ::     rn !< Radius of the center of north face of volume P
      real(8), dimension (nx*ny), intent(in)  ::     Jn !< Jacobian at the center of north face of volume P
      real(8), dimension (nx*ny), intent(in)  ::    vln !< Laminar viscosity at center of face north
      real(8), dimension (nx*ny), intent(in)  ::      p !< Pressure at center o volume P
      real(8), dimension (nx*ny), intent(in)  ::      u !< Cartesian velocity of the last iteraction
      real(8), dimension (nx*ny), intent(in)  ::      v !< Cartesian velocity of the last iteraction
      real(8), dimension (nx*ny), intent(in)  ::    Vcn !< Contravariant velocity V at north face
      character(len=*),           intent(out) :: msg(2) !< Message to be printed

      ! In the following calls, PF, TF, UF, rb and iocs are global
      ! variables of this module

      call get_cdfi(nx, ny, 2, iocs, Rg, PF, TF, UF, rb, yk, rn, p, Cdfi)

      if ( modvis == 1 ) then

         call get_cdfv(nx, ny, 2, iocs, Rg, PF, TF, UF, rb, xk, yk, rn, Jn, vln &
            , u, v, Vcn, Cdfv )

      else

         Cdfv = 0.d0

      end if

      write(msg(1), "(2(1X,A14))")    "Cdfi", "Cdfv"
      write(msg(2), "(2(1X,ES14.7))")   Cdfi,  Cdfv

   end subroutine


end module
