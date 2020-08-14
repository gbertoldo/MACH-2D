!>
!! \brief This module defines data structures and procedures related to the
!!        internal flow calculation.
!!
module mod_intflow

   use mod_class_ifile
   use mod_class_thermophysical_abstract
   use mod_intflow_geometry
   use mod_intflow_procedures


   use mod_intflow_procedures,    intflow_set_bcu => set_bcu &
      ,                           intflow_set_bcv => set_bcv &
      ,                  intflow_boundary_simplec => get_boundary_simplec_coefficients_internal_flow2 &
      ,     intflow_extrapolate_u_v_to_fictitious => get_u_v_at_fictitious_nodes_with_pl &
      ,         intflow_get_T_from_H_conservation => get_T_from_H_conservation_internal_flow

   implicit none

   ! Makes everything private, except otherwise stated
   private

   ! Public procedures
   public ::   intflow_init                          &
      ,        intflow_initial_conditions            &
      ,        intflow_time_cycle_update             &
      ,        intflow_set_bcu                       &
      ,        intflow_set_bcv                       &
      ,        intflow_set_bcpl                      &
      ,        intflow_set_bcT                       &
      ,        intflow_boundary_simplec              &
      ,        intflow_velocities_at_boundary_faces  &
      ,        intflow_extrapolate_u_v_to_fictitious &
      ,        intflow_get_T_from_H_conservation     &
      ,        intflow_calc_main_variables           &
      ,        intflow_grid_boundary

   !> \brief Data to be used in the internal flow calculation
   type, public :: type_intflow

      real(8) ::    T0 !< Stagnation temperature (K)
      real(8) ::    p0 !< Stagnation pressure (Pa)
      real(8) ::   Cp0 !< Specific heat at const. p at T0
      real(8) :: Twall !< Prescribed wall temperature

      ! Geometric parameters for internal flow calculations
      type(type_intflow_geometry) :: geom

      ! Quasi-1D solution
      real(8), allocatable, dimension(:) :: M1D !< Mach number. Isentropic flow
      real(8), allocatable, dimension(:) :: p1D !< Pressure. Isentropic flow
      real(8), allocatable, dimension(:) :: T1D !< Temperature. Isentropic flow
      real(8), allocatable, dimension(:) :: u1D !< Velocity. Isentropic flow

      ! Inflow variables
      real(8), allocatable, dimension(:) :: uin  !< Inflow velocity u
      real(8), allocatable, dimension(:) :: vin  !< Inflow velocity v
      real(8), allocatable, dimension(:) :: pin  !< Inflow pressure
      real(8), allocatable, dimension(:) :: Tin  !< Inflow temperature
      real(8), allocatable, dimension(:) :: Mw   !< Inflow Mach number
      real(8), allocatable, dimension(:) :: plin !< Inflow pressure correction
      real(8), allocatable, dimension(:) :: pina !< Inflow pressure at previous time

      ! Integrated variables of the quasi-1D flow
      real(8) :: fm1D  !< Mass flow rate
      real(8) :: Fd1D  !< Thrust due momentum transfer
      real(8) :: Fpv1D !< Thrust due pressure for vacuum discharge

   end type

contains

   !> \brief Initializes internal flow variable
   subroutine intflow_init(ifile, thermomodel, Tref, Href, Cpref, iflow) ! Output: last three
      implicit none
      class(class_ifile),                            intent(in)  :: ifile       !< Input file
      class(class_thermophysical_abstract), pointer, intent(in)  :: thermomodel !< A pointer to the thermophysical model
      real(8),                                       intent(out) :: Tref        !< Reference temperature (K)
      real(8),                                       intent(out) :: Href        !< Reference total enthalpy (m2/s2)
      real(8),                                       intent(out) :: Cpref       !< Reference Cp (J/kg.K)
      type(type_intflow),                            intent(out) :: iflow       !< Data for internal flow


      call ifile%get_value(    iflow%P0,       "P0") ! Stagnation pressure (Pa)
      call ifile%get_value(    iflow%T0,       "T0") ! Stagnation temperature (K)
      call ifile%get_value( iflow%Twall,     "Tnbc") ! Temperature on the north boundary (K) (if negative, adiabatic bc is applied)
      !call ifile%get_value(    ngeom,    "ngeom") ! File of boundary geometry

      ! Calculating reference values
      Tref  = iflow%T0
      Cpref = thermomodel%cp(Tref)
      Href  = Cpref * Tref

      iflow%Cp0 = Cpref

   end subroutine


   !> \brief Calculates initial conditions for internal flow
   subroutine intflow_initial_conditions( nx, ny, modvis, beta & ! Input
      ,              thermomodel, ye, yk, radius, rn, x, y, xp & ! Input
      ,                                                  iflow & ! InOutput
      ,                           p, T, u, v, ue, un, Uce, Vcn & ! Output
      ,                                   de, dn, ro, roe, ron ) ! Output
      integer, intent(in) :: nx     !< Number of volumes in csi direction (real+fictitious)
      integer, intent(in) :: ny     !< Number of volumes in eta direction (real+fictitious)
      integer, intent(in) :: modvis !< modvis = 0 -> Euler;  modvis = 1 -> Navier-Stokes
      real(8), intent(in) :: beta   !< Constant of the UDS/CDS mixing scheme

      class(class_thermophysical_abstract), pointer, intent(in)  :: thermomodel !< A pointer to the thermophysical model

      real(8), dimension(nx*ny), intent(in) :: ye     !< y_eta at face east of volume P
      real(8), dimension(nx*ny), intent(in) :: yk     !< y_csi at face north of volume P
      real(8), dimension(nx*ny), intent(in) :: radius !< Radius of northest corner of volume P
      real(8), dimension(nx*ny), intent(in) :: rn     !< Radius of the center of north face of volume P
      real(8), dimension(nx*ny), intent(in) :: x      !< Coord. x of the northest corner of volume P
      real(8), dimension(nx*ny), intent(in) :: y      !< Coord. x of the northest corner of volume P
      real(8), dimension(nx*ny), intent(in) :: xp     !< Coord. x of the centroid of volume P

      type(type_intflow),     intent(inout) :: iflow  !< Data for internal flow

      real(8), dimension(nx*ny), intent(out) :: p   ! Pressure at center o volume P
      real(8), dimension(nx*ny), intent(out) :: T   ! Temperature of the last iteraction
      real(8), dimension(nx*ny), intent(out) :: u   ! Cartesian velocity of the last iteraction
      real(8), dimension(nx*ny), intent(out) :: v   ! Cartesian velocity of the last iteraction

      real(8), dimension(nx*ny), intent(out) :: ue  ! Cartesian velocity u at center of east face
      real(8), dimension(nx*ny), intent(out) :: un  ! Cartesian velocity u at center of north face
      real(8), dimension(nx*ny), intent(out) :: Uce ! Contravariant velocity U at east face
      real(8), dimension(nx*ny), intent(out) :: Vcn ! Contravariant velocity V at north face

      real(8), dimension(nx*ny), intent(out) :: de  ! Simplec coef. for the contravariant velocity U (east face)
      real(8), dimension(nx*ny), intent(out) :: dn  ! Simplec coef. for the contravariant velocity V (north face)

      real(8), dimension(nx*ny), intent(out) :: ro  ! Specific mass (absolute density) at center of volumes
      real(8), dimension(nx*ny), intent(out) :: roe ! Absolute density at east face
      real(8), dimension(nx*ny), intent(out) :: ron ! Absolute density at north face

      real(8) :: gamma, Rg

      gamma = thermomodel%gm(iflow%T0)
      Rg    = thermomodel%Rg

      ! Allocating memory
      allocate( iflow%Tin(ny)  &
         ,      iflow%pin(ny)  &
         ,      iflow%pina(ny) &
         ,      iflow%plin(ny) &
         ,      iflow%uin(ny)  &
         ,      iflow%vin(ny)  &
         ,      iflow%Mw(ny)   )

      iflow%plin = 0.d0
      iflow%pina = 0.d0

      allocate( iflow%M1D(nx) &
         ,      iflow%p1D(nx) &
         ,      iflow%T1D(nx) &
         ,      iflow%u1D(nx) )

      associate(      Sg => iflow%geom%Sg &
            ,         po => iflow%P0      &
            ,         T0 => iflow%T0      &
            ,        M1D => iflow%M1D     &
            ,        p1D => iflow%p1D     &
            ,        T1D => iflow%T1D     &
            ,        u1D => iflow%u1D     &
            ,        uin => iflow%uin     &
            ,        vin => iflow%vin     &
            ,        pin => iflow%pin     &
            ,        Tin => iflow%Tin     &
            ,         Mw => iflow%Mw      &
            ,       fm1D => iflow%fm1D    &
            ,       Fd1D => iflow%Fd1D    &
            ,      Fpv1D => iflow%Fpv1D   )

         call get_initial_guess( nx, ny, modvis, beta, po, T0, gamma      & ! Input
         ,                       Rg, Sg, ye, yk, radius, rn, x, y, xp     & ! Input
         ,                       M1D, p1D, T1D, u1D, p, T, u, v, ue       & ! Output
         ,                       un, Uce, Vcn, uin, vin, pin, Tin, Mw     & ! Output
         ,                       fm1D, Fd1D, Fpv1D, de, dn, ro, roe, ron  ) ! Output

      end associate

   end subroutine


   !> \brief Calculates internal flow variables at the beginning of the new time
   !! cycle
   subroutine intflow_time_cycle_update(nx, ny, x, xp, u, thermomodel, iflow, p) ! InOutput: last two
      implicit none
      integer,                      intent(in)    :: nx          !< Number of volumes in csi direction (real+fictitious)
      integer,                      intent(in)    :: ny          !< Number of volumes in eta direction (real+fictitious)
      real(8), dimension(nx*ny),    intent(in)    :: x           !< Coord. x of the northest corner of volume P
      real(8), dimension(nx*ny),    intent(in)    :: xp          !< Coord. x of the centroid of volume P
      real(8), dimension(nx*ny),    intent(in)    :: u           !< Cartesian velocity of the last iteraction
      class(class_thermophysical_abstract) &
         ,                pointer,  intent(in)    :: thermomodel !< A pointer to the thermophysical model
      type(type_intflow),           intent(inout) :: iflow       !< Data for internal flow
      real(8), dimension(nx*ny),    intent(inout) :: p           !< Pressure at center o volume P

      real(8) :: gamma
      real(8) :: Rg

      gamma = thermomodel%gm(iflow%T0)
      Rg    = thermomodel%Rg

      ! pina = pin
      iflow%pina = iflow%pin

      ! Updating inflow conditions
      call get_uin_vin_pin_Tin_Mw( nx, ny, gamma, Rg, iflow%P0, iflow%T0, u & ! Input
      ,             iflow%uin, iflow%vin, iflow%pin, iflow%Tin, iflow%Mw ) ! Output

      call get_plin_and_p_fictitious( nx, ny, x, xp, iflow%pina, iflow%pin &
         ,                                                   iflow%plin, p ) ! Output: last two entries

   end subroutine


   !> \brief Calculates velocities at center of boundary faces according to bc
   subroutine intflow_velocities_at_boundary_faces(nx, ny, ye, u, v & ! Input
      ,                                              ue, ve, un, vn & ! InOutput
      ,                                                    Uce, Vcn ) ! Output
      implicit none
      integer,                   intent(in)    :: nx   !< Number of volumes in csi direction (real+fictitious)
      integer,                   intent(in)    :: ny   !< Number of volumes in eta direction (real+fictitious)
      real(8), dimension(nx*ny), intent(in)    :: ye   !< y_eta at face east of volume P
      real(8), dimension(nx*ny), intent(in)    :: u    !< Cartesian velocity of the last iteraction
      real(8), dimension(nx*ny), intent(in)    :: v    !< Cartesian velocity of the last iteraction
      real(8), dimension(nx*ny), intent(inout) :: ue   !< Cartesian velocity u at center of east face
      real(8), dimension(nx*ny), intent(inout) :: ve   !< Cartesian velocity v at center of east face
      real(8), dimension(nx*ny), intent(inout) :: un   !< Cartesian velocity u at center of north face
      real(8), dimension(nx*ny), intent(inout) :: vn   !< Cartesian velocity v at center of north face
      real(8), dimension(nx*ny), intent(out)   :: Uce  !< Contravariant velocity U at east face
      real(8), dimension(nx*ny), intent(out)   :: Vcn  !< Contravariant velocity V at north face

      call get_ue_un_ve_vn_at_boundary_faces( nx, ny, u, v, ue, un, ve, vn) ! InOutput: ue, un, ve, vn

      call get_Uce_Vcn_at_boundary_faces( nx, ny, ye, u, Uce, Vcn) ! Output: Uce, Vcn

   end subroutine


   !> \brief Defines the bc for pressure correction pl
   subroutine intflow_set_bcpl(nx, ny, x, xp, pl, iflow, ap, bp)
      implicit none
      integer,                     intent(in)  :: nx    !< Number of volumes in csi direction (real+fictitious)
      integer,                     intent(in)  :: ny    !< Number of volumes in eta direction (real+fictitious)
      real(8), dimension(nx*ny),   intent(in)  :: x     !< Coord. x of the northest corner of volume P
      real(8), dimension(nx*ny),   intent(in)  :: xp    !< X coord. of the centroid of volume P
      real(8), dimension(nx*ny),   intent(in)  :: pl    !< Pressure correction
      type(type_intflow),          intent(in)  :: iflow !< Data for internal flow
      real(8), dimension(nx*ny,5), intent(out) :: ap    !< Coefficients of the linear system
      real(8), dimension(nx*ny),   intent(out) :: bp    !< Source vector of the linear system

      call set_bcp( nx, ny, x, xp, iflow%plin, pl, ap, bp) ! Output: last two

   end subroutine


   !> \brief Defines bc for T
   !! CAUTION: Twall must be defined in the fictitious volumes too
   subroutine intflow_set_bcT( nx, ny, x, xp, T, iflow, at, bt) ! Output: last two
      implicit none
      integer,                     intent(in)  :: nx    !< Number of volumes in csi direction (real+fictitious)
      integer,                     intent(in)  :: ny    !< Number of volumes in eta direction (real+fictitious)
      real(8), dimension(nx*ny),   intent(in)  :: x     !< Coord. x of the northest corner of volume P
      real(8), dimension(nx*ny),   intent(in)  :: xp    !< X coord. of the centroid of volume P
      real(8), dimension(nx*ny),   intent(in)  :: T     !< Temperature of the last iteraction
      type(type_intflow),          intent(in)  :: iflow !< Data for internal flow
      real(8), dimension(nx*ny,9), intent(out) :: at    !< Coefficients of the linear system
      real(8), dimension(nx*ny),   intent(out) :: bt    !< Source vector of the linear system

      ! Inner variables
      integer                :: ccTw ! ccTw = 0 -> adiabatic;  ccTw = 1 -> prescribed temperature
      real(8), dimension(nx) ::  Tw  ! Wall temperature

      ccTw = 0
      if ( iflow%Twall > 0.d0 ) then
         ccTw = 1
         Tw   = iflow%Twall
      end if

      call set_bcT( nx, ny, ccTw, x, xp, iflow%Tin, Tw, T, at, bt) ! Output: last two

   end subroutine


   !> \brief Calculates the main variables of internal flow
   subroutine intflow_calc_main_variables(nx, ny, re, roe, u, Uce, msg)
      implicit none
      integer,                   intent(in)  :: nx     !< Number of volumes in csi direction (real+fictitious)
      integer,                   intent(in)  :: ny     !< Number of volumes in eta direction (real+fictitious)
      real(8), dimension(nx*ny), intent(in)  :: re     !< Radius of the center of east face of volume P
      real(8), dimension(nx*ny), intent(in)  :: roe    !< Absolute density at east face
      real(8), dimension(nx*ny), intent(in)  :: u      !< Cartesian velocity of the last iteraction
      real(8), dimension(nx*ny), intent(in)  :: Uce    !< Contravariant velocity U at east face
      character(len=*),          intent(out) :: msg(2) !< Message to be printed

      real(8) :: fmi    !< Mass flow rate at entrance
      real(8) :: fme    !< Mass flow rate at exit
      real(8) :: Fd     !< Thrust

      call get_mass_flow_rate_and_thrust(nx, ny, re, roe, u, Uce, fmi, fme, Fd) ! Output: last three entries

      write(msg(1), "(3(1X,A23))")     "fmi", "fme/fmi", "Fd"
      write(msg(2), "(3(1X,ES23.16))")   fmi,   fme/fmi,  Fd

   end subroutine


   !> \brief Defines the south and north boundary of the domain
   subroutine intflow_grid_boundary(nx, ny, x, y, iflow, a1) ! Output: last four
      implicit none
      integer,                   intent(in)  ::  nx !< Number of volumes in the csi direction (real+fictitious)
      integer,                   intent(in)  ::  ny !< Number of volumes in the eta direction (real+fictitious)
      real(8), dimension(nx*ny), intent(out) ::   x !< Coord. x at the northeast corner of the volume P (m)
      real(8), dimension(nx*ny), intent(out) ::   y !< Coord. y at the northeast corner of the volume P (m)
      type(type_intflow),        intent(out) :: iflow  !< Variables related to internal flow
      real(8),                   intent(out) :: a1
! TODO (guilherme#1#): Read or calculate a1

      integer :: ig
      integer :: kg

      kg = 1
      a1 = 1E-6

      call get_boundary_nodes( nx, ny, ig, iflow%geom%Sg, iflow%geom%rcg) ! Output: last three entries
      call set_grid_internal_flow(kg, nx, ny, a1, x, y) ! Last 2 are output

   end subroutine

end module
