!>
!! \brief This module defines data structures and procedures related to the
!!        internal flow calculation.
!!
module mod_intflow

   use mod_class_ifile
   use mod_class_thermophysical_abstract
   use mod_intflow_geometry
   use coefficients

   implicit none

   ! Makes everything private, except otherwise stated
   private

   ! Public procedures
   public ::   intflow_init                           &
      ,        intflow_initial_conditions             &
      ,        intflow_time_cycle_update              &
      ,        intflow_set_bcu                        &
      ,        intflow_set_bcv                        &
      ,        intflow_set_bcpl                       &
      ,        intflow_set_bcT                        &
      ,        intflow_boundary_simplec               &
      ,        intflow_velocities_at_boundary_faces   &
      ,        intflow_extrapolate_u_v_to_fictitious  &
      ,        intflow_boundary_T_from_H_conservation &
      ,        intflow_calc_main_variables            &
      ,        intflow_grid_boundary                  &
      ,        get_ig                                 &
      ,        get_mass_flow_rate_and_thrust


   !> \brief Data to be used in the internal flow calculation
   type, public :: type_intflow

      real(8) ::    T0 !< Stagnation temperature (K)
      real(8) ::    p0 !< Stagnation pressure (Pa)
      real(8) ::   Cp0 !< Specific heat at const. p at T0
      real(8) :: Twall !< Prescribed wall temperature

      ! Geometric parameters for internal flow calculations
      type(type_intflow_geometry) :: geom
      type(class_ifile)           :: gfile

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

      ! Inner variables
      character(len=200) :: ngeom

      call ifile%get_value(    iflow%P0,       "P0") ! Stagnation pressure (Pa)
      call ifile%get_value(    iflow%T0,       "T0") ! Stagnation temperature (K)
      call ifile%get_value( iflow%Twall,     "Tnbc") ! Temperature on the north boundary (K) (if negative, adiabatic bc is applied)
      call ifile%get_value(       ngeom,    "ngeom") ! File of boundary geometry

      ! Calculating reference values
      Tref  = iflow%T0
      Cpref = thermomodel%cp(Tref)
      Href  = Cpref * Tref

      iflow%Cp0 = Cpref

      ! Reading the geometry file
      call iflow%gfile%init(trim(ngeom), "&")
      call iflow%gfile%load()

   end subroutine


   !> \brief Defines the south and north boundary of the domain
   subroutine intflow_grid_boundary(nx, ny, x, y, iflow) ! Output: last four
      implicit none
      integer,                   intent(in)    ::    nx !< Number of volumes in the csi direction (real+fictitious)
      integer,                   intent(in)    ::    ny !< Number of volumes in the eta direction (real+fictitious)
      real(8), dimension(nx*ny), intent(out)   ::     x !< Coord. x at the northeast corner of the volume P (m)
      real(8), dimension(nx*ny), intent(out)   ::     y !< Coord. y at the northeast corner of the volume P (m)
      type(type_intflow),        intent(inout) :: iflow !< Variables related to internal flow

      character(len=50) :: GID ! Geometry ID
      real(8)           :: Rth ! Throat radius
      real(8)           :: Rc  ! Throat radius of curvature

      ! Reading the geometry ID
      call iflow%gfile%get_value( GID, "ID")

      ! Choosing between options
      if ( trim(GID) == "N01" ) then

         call get_grid_boundary_nozzle01(iflow%gfile, nx, ny, x, y) ! Output: last two

         call iflow%gfile%get_value( Rth, "Rth")
         call iflow%gfile%get_value(  Rc, "Rc2")

         iflow%geom%Sg = acos(-1.d0) * Rth**2
         iflow%geom%rcg = Rc

      else
         write(*,*) "intflow_grid_boundary:"
         write(*,*) "Unknown option: ", trim(GID)
         write(*,*) "Stopping..."
         stop
      end if

   end subroutine


   !> \brief Calculates u and v on north boundary assuming slip
   !! boundary condition.
   subroutine get_u_v_north_slip(nx, ny, x, xp, xk, yk, u, v, ubn, vbn) ! Output: last two
      implicit none
      integer, intent(in) :: nx     ! Number of volumes in csi direction (real+fictitious)
      integer, intent(in) :: ny     ! Number of volumes in eta direction (real+fictitious)
      real(8), dimension(nx*ny), intent(in)  :: x   !< Coord. x of the northest corner of volume P
      real(8), dimension(nx*ny), intent(in)  :: xp  !< X coord. of the centroid of volume P
      real(8), dimension(nx*ny), intent(in)  :: xk  !< x_csi at center of face north
      real(8), dimension(nx*ny), intent(in)  :: yk  !< y_csi at center of face north
      real(8), dimension(nx*ny), intent(in)  :: u   !< Cartesian velocity of the present iteraction
      real(8), dimension(nx*ny), intent(in)  :: v   !< Cartesian velocity of the present iteraction
      real(8), dimension(nx),    intent(out) :: ubn !< Cartesian velocity on north boundary
      real(8), dimension(nx),    intent(out) :: vbn !< Cartesian velocity on south boundary

      ! Inner variables
      integer :: i, j, np, npe, npw, npee, npww
      real(8) :: signal, tx, ty, t, Up

      j = ny-1

      do i = 2, nx-1

         np   = nx * (j-1) + i

         ! Calculating the modulus of vec(u) at P
         Up = sqrt(u(np)**2+v(np)**2)

         ! Is (u,v) pointing toward (xk,yk) or backward (xk,yk)?
         ! To discover that, calculates the dot product and gets its signal.
         ! SIGN(A,B) returns the value of A with the sign of B
         signal = sign(1.d0, xk(np) * u(np) + yk(np) * v(np) )

         ! Getting the unitary vector t tangent to the wall
         t  = sqrt( xk(np)**2 + yk(np)**2 )

         tx = xk(np) / t
         ty = yk(np) / t

         ! Calculating the velocity vector on the wall
         ubn(i) = signal * Up * tx
         vbn(i) = signal * Up * ty

         ! Note that this vector has the same modulus that (u,v)_P and points
         ! in a direction locally tangent to the wall. Having the same modulus
         ! on P ensures that the total enthalpy will be conserved on this cell.

      end do

      ! Extrapolating the velocity to the fictitious of west boundary
      i = 1
      j = ny-1
      np   = nx * (j-1) + i
      npe  = np + 1
      npee = npe + 1

      ubn(1) = ubn(2) + (ubn(2)-ubn(3)) * (xp(npe)-x(np))/(xp(npee)-xp(npe))
      vbn(1) = vbn(2) + (vbn(2)-vbn(3)) * (xp(npe)-x(np))/(xp(npee)-xp(npe))


      ! Extrapolating the velocity to the fictitious of east boundary
      i = nx
      j = ny-1
      np   = nx * (j-1) + i
      npw  = np - 1
      npww = npw - 1

      ubn(nx) = ubn(nx-1) + (ubn(nx-1)-ubn(nx-2)) * (x(npw)-xp(npw))/(xp(npw)-xp(npww))
      vbn(nx) = vbn(nx-1) + (vbn(nx-1)-vbn(nx-2)) * (x(npw)-xp(npw))/(xp(npw)-xp(npww))

   end subroutine


   !> \brief Calculates the coefficients and source for boundary conditions of u
   subroutine intflow_set_bcu( nx, ny, modvis, x, xp, xk, yk, u, v, au, b) ! Output: last two entries
      implicit none
      integer, intent(in) :: nx  !< Number of volumes in csi direction (real+fictitious)
      integer, intent(in) :: ny  !< Number of volumes in eta direction (real+fictitious)
      integer, intent(in) :: modvis !< modvis = 0 -> Euler;  modvis = 1 -> Navier-Stokes
      real(8), dimension (nx*ny),   intent(in)  :: x  !< Coord. x of the northest corner of volume P
      real(8), dimension (nx*ny),   intent(in)  :: xp !< X coord. of the centroid of volume P
      real(8), dimension (nx*ny),   intent(in)  :: xk !< x_csi at center of face north
      real(8), dimension (nx*ny),   intent(in)  :: yk !< y_csi at center of face north
      real(8), dimension (nx*ny),   intent(in)  :: u  !< Cartesian velocity of the last iteraction
      real(8), dimension (nx*ny),   intent(in)  :: v  !< Cartesian velocity of the last iteraction
      real(8), dimension (nx*ny,9), intent(out) :: au !< Coefficients of the linear system for u
      real(8), dimension (nx*ny),   intent(out) :: b  !< Source vector of the linear system
      ! Auxiliary variables
      integer :: i, j, np, npe, npee, npw, npww

      real(8), dimension(nx) :: ubn ! Cartesian velocity on north boundary
      real(8), dimension(nx) :: vbn ! Cartesian velocity on south boundary

      ! West boundary (corners are not included)
      i = 1
      do j = 2, ny-1

         np   = nx * (j-1) + i
         npe  = np + 1
         npee = npe + 1

         au(np,:) =  0.d0
         au(np,6) = -1.d0
         au(np,5) =  1.d0
         b(np) = ( 2.d0 * xp(npe) ) / ( xp(npee) - xp(npe) ) * ( u(npe) - u(npee) )
      end do

      ! East boundary (corners are not included)
      i = nx
      do j = 2, ny-1

         np   = nx * (j-1) + i
         npw  = np - 1
         npww = npw - 1

         au(np,:) =  0.d0
         au(np,4) = -1.d0
         au(np,5) =  1.d0
         b(np) = 2.d0 * ( x(npw) - xp(npw) ) / ( xp(npw) - xp(npww) ) * ( u(npw) - u(npww) )
      end do

      ! South boundary ( SW and SE corners are included)
      j = 1
      do i = 1, nx

         np   = nx * (j-1) + i

         au(np,:) =  0.d0
         au(np,8) = -1.d0
         au(np,5) =  1.d0

         b(np) = 0.d0
      end do

      ! North boundary ( NW and NE corners are included)
      j = ny
      if ( modvis == 1 ) then ! Navier Stokes
         do i = 1, nx
            np   = nx * (j-1) + i

            au(np,:) = 0.d0
            au(np,2) = 1.d0
            au(np,5) = 1.d0

            b(np) = 0.d0
         end do
      else ! Euler
         ! Calculates u and v on north boundary assuming slip
         ! boundary condition.
         call get_u_v_north_slip(nx, ny, x, xp, xk, yk, u, v, ubn, vbn) ! Output: last two

         do i = 1, nx
            np   = nx * (j-1) + i

            au(np,:) =  0.d0
            au(np,2) =  1.d0
            au(np,5) =  1.d0

            b(np) = 2.d0 * ubn(i)
         end do
      end if
   end subroutine intflow_set_bcu


   !> \brief Calculates the coefficients and source for boundary conditions of v
   subroutine intflow_set_bcv( nx, ny, modvis, x, xp, xk, yk, u, v, av, b) ! Output: last two entries
      implicit none
      integer, intent(in) :: nx     !< Number of volumes in csi direction (real+fictitious)
      integer, intent(in) :: ny     !< Number of volumes in eta direction (real+fictitious)
      integer, intent(in) :: modvis !< modvis = 0 -> Euler;  modvis = 1 -> Navier-Stokes
      real(8), dimension (nx*ny),   intent(in)  :: x  !< Coord. x of the northest corner of volume P
      real(8), dimension (nx*ny),   intent(in)  :: xp !< X coord. of the centroid of volume P
      real(8), dimension (nx*ny),   intent(in)  :: xk !< x_csi at center of face north
      real(8), dimension (nx*ny),   intent(in)  :: yk !< y_csi at center of face north
      real(8), dimension (nx*ny),   intent(in)  :: u  !< Cartesian velocity of the last iteraction
      real(8), dimension (nx*ny),   intent(in)  :: v  !< Cartesian velocity of the last iteraction
      real(8), dimension (nx*ny,9), intent(out) :: av !< Coefficients of the linear system for v
      real(8), dimension (nx*ny),   intent(out) :: b  !< Source vector of the linear system

      ! Auxiliary variables
      integer :: i, j, np, npe, npee, npw, npww
      real(8), dimension(nx) :: ubn ! Cartesian velocity on north boundary
      real(8), dimension(nx) :: vbn ! Cartesian velocity on south boundary

      ! West boundary (corners are not included)
      i = 1
      do j = 2, ny-1

         np   = nx * (j-1) + i
         npe  = np + 1
         npee = npe + 1

         av(np,:) = 0.d0
         av(np,6) = 1.d0
         av(np,5) = 1.d0
         b(np) = 0.d0
      end do

      ! East boundary (corners are not included)
      i = nx
      do j = 2, ny-1

         np   = nx * (j-1) + i
         npw  = np - 1
         npww = npw - 1

         av(np,:) =  0.d0
         av(np,4) = -1.d0
         av(np,5) =  1.d0
         b(np) = ( v(npw) - v(npww) ) * 2.d0 * ( x(npw) - xp(npw) ) / ( xp(npw) - xp(npww) )
      end do

      ! South boundary ( SW and SE corners are included)
      j = 1
      do i = 1, nx

         np   = nx * (j-1) + i

         av(np,:) = 0.d0
         av(np,8) = 1.d0
         av(np,5) = 1.d0

         b(np) = 0.d0
      end do

      ! North boundary ( NW and NE corners are included)
      j = ny
      if ( modvis == 1 ) then ! Navier-Stokes
         do i = 1, nx
            np   = nx * (j-1) + i

            av(np,:) = 0.d0
            av(np,2) = 1.d0
            av(np,5) = 1.d0

            b(np) = 0.d0
         end do
      else ! Euler
         ! Calculates u and v on north boundary assuming slip
         ! boundary condition.
         call get_u_v_north_slip(nx, ny, x, xp, xk, yk, u, v, ubn, vbn) ! Output: last two

         do i = 1, nx
            np   = nx * (j-1) + i

            av(np,:) =  0.d0
            av(np,2) =  1.d0
            av(np,5) =  1.d0

            b(np) = 2.d0 * vbn(i)
         end do
      end if
   end subroutine intflow_set_bcv


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

      ! Auxiliary variables
      integer :: i, j, np, npw, npww
      ! West boundary
      i = 1
      do j = 2, ny-1
         np = nx * (j-1) + i

         ap(np,:) = 0.d0
         ap(np,4) = 1.d0
         ap(np,3) = 1.d0

         bp(np) = 2.d0 * iflow%plin(j)
      end do

      ! East boundary
      i = nx
      do j = 2, ny-1
         np   = nx * (j-1) + i
         npw  = np - 1
         npww = npw - 1

         ap(np,:) =  0.d0
         ap(np,2) = -1.d0
         ap(np,3) =  1.d0

         bp(np) = 2.d0 * ( x(npw) - xp(npw) ) / ( xp(npw) - xp(npww) ) * ( pl(npw) - pl(npww) )
      end do

      ! South boundary
      j = 1
      do i = 2, nx-1
         np  = nx * (j-1) + i

         ap(np,:) =  0.d0
         ap(np,5) = -1.d0
         ap(np,3) =  1.d0

         bp(np) = 0.d0
      end do

      ! North boundary
      j = ny
      do i = 2, nx-1
         np  = nx * (j-1) + i

         ap(np,:) =  0.d0
         ap(np,1) = -1.d0
         ap(np,3) =  1.d0

         bp(np) = 0.d0
      end do

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

      ! Auxiliary variables
      integer                :: i, j, np, npw, npww
      integer                :: ccTw ! ccTw = 0 -> adiabatic;  ccTw = 1 -> prescribed temperature
      real(8), dimension(nx) ::  Tw  ! Wall temperature

      ccTw = 0
      if ( iflow%Twall > 0.d0 ) then
         ccTw = 1
         Tw   = iflow%Twall
      end if



      ! West boundary
      i = 1
      do j = 2, ny-1
         np   = nx * (j-1) + i

         at(np,:) = 0.d0
         at(np,6) = 1.d0
         at(np,5) = 1.d0

         bt(np) = 2.d0 * iflow%Tin(j)
      end do

      ! East boundary
      i = nx
      do j = 2, ny-1
         np   = nx * (j-1) + i
         npw  = np - 1
         npww = npw - 1

         at(np,:) =  0.d0
         at(np,4) = -1.d0
         at(np,5) =  1.d0

         bt(np) = 2.d0 * (x(npw) - xp(npw))/(xp(npw) - xp(npww)) * (T(npw) - T(npww))
      end do

      ! South boundary ( SW and SE corners are included)
      j = 1
      do i = 1, nx
         np = nx * (j-1) + i

         at(np,:) =  0.d0
         at(np,8) = -1.d0
         at(np,5) =  1.d0

         bt(np) = 0.d0
      end do

      ! North boundary ( NW and NE corners are included)
      j = ny
      if ( ccTw == 0 ) then ! Adiabatic
         do i = 1, nx
            np = nx * (j-1) + i

            at(np,:) =  0.d0
            at(np,2) = -1.d0
            at(np,5) =  1.d0

            bt(np) = 0.d0
         end do
      else ! Prescribed temperature
         do i = 1, nx
            np = nx * (j-1) + i

            at(np,:) = 0.d0
            at(np,2) = 1.d0
            at(np,5) = 1.d0

            bt(np) = 2.d0 * Tw(i)
         end do
      end if

   end subroutine


   !> \brief Calculates the boundary temperature based on the  boundary conditions
   !! and on conservation of the total enthalpy (valid for Euler model with constant
   !! thermo-physical coefficients). Temperature is extrapolated to fictitious
   !! volumes using the CDS scheme.
   subroutine intflow_boundary_T_from_H_conservation(nx, ny, Cpref, Href, iflow &
         ,                             u, ue, un, v, ve, vn, T, Tbe, Tbw, Tbn, Tbs ) ! Output: last 5. T is inout
      implicit none
      integer, intent(in) :: nx    !< Number of volumes in csi direction (real+fictitious)
      integer, intent(in) :: ny    !< Number of volumes in eta direction (real+fictitious)
      real(8), intent(in) :: Cpref !< Free stream Cp (J/kg.K)
      real(8), intent(in) :: Href  !< Free stream total enthalpy (m2/s2)
      type(type_intflow),         intent(in)    :: iflow !< Data for internal flow
      real(8), dimension (nx*ny), intent(in)    :: u    !< Cartesian velocity at P
      real(8), dimension (nx*ny), intent(in)    :: ue   !< Cartesian velocity u at center of east face
      real(8), dimension (nx*ny), intent(in)    :: un   !< Cartesian velocity u at center of north face
      real(8), dimension (nx*ny), intent(in)    :: v    !< Cartesian velocity at P
      real(8), dimension (nx*ny), intent(in)    :: ve   !< Cartesian velocity v at center of east face
      real(8), dimension (nx*ny), intent(in)    :: vn   !< Cartesian velocity v at center of north face
      real(8), dimension (nx*ny), intent(inout) :: T    !< Temperature at P
      real(8), dimension (ny),    intent(out)   :: Tbe  !< Temperature over the  east boundary (K)
      real(8), dimension (ny),    intent(out)   :: Tbw  !< Temperature over the  west boundary (K)
      real(8), dimension (nx),    intent(out)   :: Tbn  !< Temperature over the north boundary (K)
      real(8), dimension (nx),    intent(out)   :: Tbs  !< Temperature over the south boundary (K)

      ! Inner variables

      integer :: i, j, np, npe, npn, nps


      ! Temperature on the east boundary

      i = nx-1

      do j = 2, ny-1

         np   = nx * (j-1) + i

         Tbe(j) = ( Href - (ue(np)**2+ve(np)**2) / 2.d0 ) / Cpref

         ! Extrapolating to fictitious

         npe  = np + 1

         T(npe) = 2.d0 * Tbe(j) - T(np)

      end do


      ! Temperature on the west boundary

      i = 1

      do j = 2, ny-1

         np   = nx * (j-1) + i

         Tbw(j) = iflow%Tin(j)

         ! Extrapolating to fictitious

         npe  = np + 1

         T(np) = 2.d0 * Tbw(j) - T(npe)

      end do


      ! Temperature on the north boundary

      j = ny-1

      do i = 2, nx-1

         np   = nx * (j-1) + i

         Tbn(i) = ( Href - (un(np)**2+vn(np)**2) / 2.d0 ) / Cpref

         ! Extrapolating to fictitious

         npn  = np + nx

         T(npn) = 2.d0 * Tbn(i) - T(np)

      end do


      ! Temperature on the north boundary

      j = 1

      do i = 2, nx-1

         np   = nx * (j-1) + i

         Tbs(i) = ( Href - (un(np)**2+vn(np)**2) / 2.d0 ) / Cpref

         ! Extrapolating to fictitious

         npn  = np + nx

         T(np) = 2.d0 * Tbs(i) - T(npn)

      end do

      ! Extrapolating to the SW corner
      i = 1
      j = 1
      np   = nx * (j-1) + i
      npn  = np + nx
      T(np) = T(npn)

      ! Extrapolating to the SE corner
      i = nx
      j = 1
      np   = nx * (j-1) + i
      npn  = np + nx

      T(np) = T(npn)

      ! Extrapolating to the NW corner
      i = 1
      j = ny
      np   = nx * (j-1) + i
      nps  = np - nx

      T(np) = T(nps)

      ! Extrapolating to the NE corner
      i = nx
      j = ny
      np   = nx * (j-1) + i
      nps  = np - nx

      T(np) = T(nps)

   end subroutine intflow_boundary_T_from_H_conservation


   !> \brief Calculates velocities on boundaries according to boundary conditions
   subroutine intflow_velocities_at_boundary_faces( nx, ny, modvis & ! Input
      ,                               x, xp, xe, ye, xk, yk, u, v & ! Input
      ,                                  ue, ve, un, vn, Uce, Vcn ) ! Input and Output
      implicit none
      integer, intent(in) :: nx     !< Number of volumes in csi direction (real+fictitious)
      integer, intent(in) :: ny     !< Number of volumes in eta direction (real+fictitious)
      integer, intent(in) :: modvis !< modvis = 0 -> Euler;  modvis = 1 -> Navier-Stokes
      real(8), dimension (nx*ny),   intent(in)  :: x   !< Coord. x of the northest corner of volume P
      real(8), dimension (nx*ny),   intent(in)  :: xp  !< X coord. of the centroid of volume P
      real(8), dimension (nx*ny), intent(in)    :: xe  !< x_eta at center of face east
      real(8), dimension (nx*ny), intent(in)    :: ye  !< y_eta at center of face east
      real(8), dimension (nx*ny), intent(in)    :: xk  !< x_csi at center of face north
      real(8), dimension (nx*ny), intent(in)    :: yk  !< y_csi at center of face north
      real(8), dimension (nx*ny), intent(in)    :: u   !< Cartesian velocity of the present iteraction
      real(8), dimension (nx*ny), intent(in)    :: v   !< Cartesian velocity of the present iteraction
      real(8), dimension (nx*ny), intent(inout) :: ue  !< Cartesian velocity u at center of east face
      real(8), dimension (nx*ny), intent(inout) :: ve  !< Cartesian velocity v at center of east face
      real(8), dimension (nx*ny), intent(inout) :: un  !< Cartesian velocity u at center of north face
      real(8), dimension (nx*ny), intent(inout) :: vn  !< Cartesian velocity v at center of north face
      real(8), dimension (nx*ny), intent(inout) :: Uce !< Contravariant velocity U at east face
      real(8), dimension (nx*ny), intent(inout) :: Vcn !< Contravariant velocity V at north face

      ! Auxiliary variables
      integer :: i, j, np, npn, npe

      real(8), dimension(nx) :: ubn ! Cartesian velocity on north boundary
      real(8), dimension(nx) :: vbn ! Cartesian velocity on south boundary

      ! West
      i = 1
      do j = 2, ny-1

         np   = nx * (j-1) + i
         npe  = np + 1

         ue(np) = ( u(np) + u(npe) ) / 2.d0

         ve(np) = 0.d0

         Uce(np) = ue(np) * ye(np) - ve(np) * xe(np)

      end do

      ! East
      i = nx-1
      do j = 2, ny-1

         np   = nx * (j-1) + i
         npe  = np + 1

         ue(np) = ( u(np) + u(npe) ) / 2.d0

         ve(np) = ( v(np) + v(npe) ) / 2.d0

         Uce(np) = ue(np) * ye(np) - ve(np) * xe(np)

      end do

      ! South
      j = 1
      do i = 2, nx-1

         np   = nx * (j-1) + i
         npn  = np + nx

         un(np) = ( u(np) + u(npn) ) / 2.d0

         vn(np) = 0.d0

         Vcn(np) = 0.d0

      end do

      ! North boundary

      if ( modvis == 0) then ! Euler
         call get_u_v_north_slip(nx, ny, x, xp, xk, yk, u, v, ubn, vbn) ! Output: last two
      else ! NS
         ubn = 0.d0
         vbn = 0.d0
      end if

      j = ny-1
      do i = 2, nx-1

         np   = nx * (j-1) + i
         npn  = np + nx

         un(np) = ubn(i)
         vn(np) = vbn(i)

         Vcn(np) = 0.d0

      end do

   end subroutine intflow_velocities_at_boundary_faces


   !> \brief Calculates the SIMPLEC coefficients on the boundary faces
   subroutine intflow_boundary_simplec( nx, ny, modvis & ! Input
      ,                           due, dve, dun, dvn, de, dn ) ! InOutput: last 6
      implicit none
      integer, intent(in) :: nx  !< Number of volumes in csi direction (real+fictitious)
      integer, intent(in) :: ny  !< Number of volumes in eta direction (real+fictitious)
      integer, intent(in) :: modvis !< Viscosity model (0=Euler, 1=NS)
      real(8), dimension(nx*ny), intent(inout) :: due  !< Simplec coef. for the cartesian velocity u (east face)
      real(8), dimension(nx*ny), intent(inout) :: dve  !< Simplec coef. for the cartesian velocity v (east face)
      real(8), dimension(nx*ny), intent(inout) :: dun  !< Simplec coef. for the cartesian velocity u (north face)
      real(8), dimension(nx*ny), intent(inout) :: dvn  !< Simplec coef. for the cartesian velocity v (north face)
      real(8), dimension(nx*ny), intent(inout) :: de   !< Simplec coef. for the contravariant velocity U (east face)
      real(8), dimension(nx*ny), intent(inout) :: dn   !< Simplec coef. for the contravariant velocity V (north face)

      ! Auxiliary variables
      integer :: i, j, np, npe, npw, npn, nps

      ! West boundary

      i = 1

      do j = 2, ny-1

         np   = nx * (j-1) + i

         npe  = np + 1

         due(np) = due(npe)

         dve(np) = 0.d0

         de(np) = de(npe)

      end do

      ! East boundary

      i = nx-1

      do j = 2, ny-1

         np   = nx * (j-1) + i

         npw  = np - 1

         due(np) = due(npw)

         dve(np) = dve(npw)

         de(np) = de(npw)

      end do

      ! South boundary

      j = 1

      do i = 2, nx-1

         np   = nx * (j-1) + i

         npn  = np + nx

         dun(np) = dun(npn)

         dvn(np) = 0.d0

         dn(np) = 0.d0

      end do

      ! North boundary

      j = ny-1

      do i = 2, nx-1

         np   = nx * (j-1) + i

         dun(np) = 0.d0

         dvn(np) = 0.d0

         dn(np) = 0.d0

      end do

   end subroutine intflow_boundary_simplec


   !> \brief Extrapolates u and v to boundary according to boundary conditions
   subroutine intflow_extrapolate_u_v_to_fictitious( nx, ny, modvis, x, xp, xk, yk, u, v) ! InOutput: last two entries
      implicit none
      integer, intent(in) :: nx  !< Number of volumes in csi direction (real+fictitious)
      integer, intent(in) :: ny  !< Number of volumes in eta direction (real+fictitious)
      integer, intent(in) :: modvis !< modvis = 0 -> Euler;  modvis = 1 -> Navier-Stokes
      real(8), dimension (nx*ny), intent(in)    :: x  !< Coord. x of the northest corner of volume P
      real(8), dimension (nx*ny), intent(in)    :: xp !< Coord. x of the centroid of volume P
      real(8), dimension (nx*ny), intent(in)    :: xk !< x_csi at center of face north
      real(8), dimension (nx*ny), intent(in)    :: yk !< y_csi at center of face north
      real(8), dimension (nx*ny), intent(inout) :: u  !< Cartesian velocity
      real(8), dimension (nx*ny), intent(inout) :: v  !< Cartesian velocity

      ! Auxiliary variables
      integer :: i, j, np, npe, npee, npw, npww, nps, npn

      real(8), dimension(nx) :: ubn ! Cartesian velocity on north boundary
      real(8), dimension(nx) :: vbn ! Cartesian velocity on south boundary

      ! West boundary
      i = 1
      do j = 2, ny-1

         np   = nx * (j-1) + i
         npe  = np + 1
         npee = np + 2

         u(np) = u(npe) + 2.d0 * xp(npe) / ( xp(npee)-xp(npe) ) * ( u(npe) - u(npee) )
         v(np) = -v(npe)

      end do

      ! East boundary
      i = nx
      do j = 2, ny-1

         np   = nx * (j-1) + i
         npw  = np - 1
         npww = np - 2

         u(np) = u(npw) + 2.d0 * (x(npw)-xp(npw)) / (xp(npw)-xp(npww)) * (u(npw)-u(npww))
         v(np) = v(npw) + 2.d0 * (x(npw)-xp(npw)) / (xp(npw)-xp(npww)) * (v(npw)-v(npww))

      end do

      ! South boundary (Corners SW and SE are included)
      j = 1
      do i = 1, nx
         np   = nx * (j-1) + i
         npn  = np + nx

         u(np) =  u(npn)
         v(np) = -v(npn)

      end do

      ! North boundary (Corners NW and NE are included)
      j = ny
      if ( modvis == 1 ) then ! Viscous

         do i = 1, nx
            np   = nx * (j-1) + i
            nps  = np - nx

            u(np) = -u(nps)
            v(np) = -v(nps)

         end do

      else ! Inviscid

         ! Calculates u and v on north boundary using slip boundary condition
         call get_u_v_north_slip(nx, ny, x, xp, xk, yk, u, v, ubn, vbn) ! Output: last two

         do i = 1, nx
            np   = nx * (j-1) + i
            nps  = np - nx

            u(np) = -u(nps) + 2.d0 * ubn(i)
            v(np) = -v(nps) + 2.d0 * vbn(i)

         end do

      end if
   end subroutine intflow_extrapolate_u_v_to_fictitious


   !> \brief Calculates u, v, p, T and M in the entrace
   subroutine get_uin_vin_pin_Tin_Mw( nx, ny, gamma, Rg, po, T0, u & ! Input
      ,                             uin, vin, pin, Tin, Mw)        ! Output
      implicit none
      integer, intent(in) :: nx    !< Number of volumes in csi direction (real+fictitious)
      integer, intent(in) :: ny    !< Number of volumes in eta direction (real+fictitious)
      real(8), intent(in) :: gamma !< gamma = Cpo / Cvo in the chamber
      real(8), intent(in) :: Rg    !< Perfect gas constant
      real(8), intent(in) :: po    !< Stagnation pressure in the chamber
      real(8), intent(in) :: T0    !< Stagnation temperature in the chamber
      real(8), dimension (nx*ny), intent(in)  :: u   !< Cartesian velocity of the last iteraction
      real(8), dimension (ny),    intent(out) :: uin !< Velocity u in the entrance
      real(8), dimension (ny),    intent(out) :: vin !< Velocity v in the entrance
      real(8), dimension (ny),    intent(out) :: pin !< Pressure in the entrance
      real(8), dimension (ny),    intent(out) :: Tin !< Temperature in the entrance
      real(8), dimension (ny),    intent(out) :: Mw  !< Mach number in the entrance

      ! Auxiliary variables
      integer :: i, j, np, npe

      i = 1
      do j = 2, ny-1

         np   = nx * (j-1) + i

         npe  = np + 1

         uin(j) = ( u(np) + u(npe) ) / 2.d0

         vin(j) = 0.d0 ! From boundary conditions

         Tin(j) = T0 - (gamma-1.d0) / ( 2.d0 * gamma * Rg ) * ( uin(j)**2 + vin(j)**2 )

         Mw(j) = dsqrt( ( uin(j)**2 + vin(j)**2 ) / ( gamma * Rg * Tin(j) ) )

         pin(j) = po * (1.d0 + 0.5d0 * (gamma-1.d0) * Mw(j)**2 ) ** ( -gamma / (gamma-1.d0) )

      end do

   end subroutine get_uin_vin_pin_Tin_Mw


   !> \brief Calculates p' in the entrance and extrapolates p to fictitious
   subroutine get_plin_and_p_fictitious( nx, ny, x, xp, pina, pin, plin, p) ! Output: last two entries
      implicit none
      integer, intent(in) :: nx    !< Number of volumes in csi direction (real+fictitious)
      integer, intent(in) :: ny    !< Number of volumes in eta direction (real+fictitious)
      real(8), dimension(nx*ny), intent(in)    :: x    !< Coord. x of the northest corner of volume P
      real(8), dimension(nx*ny), intent(in)    :: xp   !< X coord. of the centroid of volume P
      real(8), dimension(ny),    intent(in)    :: pina !< Pressure in the entrance at a time step before
      real(8), dimension(ny),    intent(in)    :: pin  !< Pressure in the entrance
      real(8), dimension(ny),    intent(out)   :: plin !< Pressure correction in the entrance
      real(8), dimension(nx*ny), intent(inout) :: p    !< Pressure at center of volumes

      ! Auxiliary variables
      integer :: i, j, np, npsw, nps, npse, npw, npe, npnw, npn, npne, npww

      ! West boundary

      i = 1

      do j = 2, ny-1

         np   = nx * (j-1) + i

         npe  = np + 1

         p(np) = 2.d0 * pin(j) - p(npe)

         plin(j) = pin(j) - pina(j)

      end do

      ! East boundary
      i = nx
      do j = 2, ny-1
         np   = nx * (j-1) + i
         npw  = np - 1
         npww = npw - 1

         p(np) = p(npw) + 2.d0 * ( x(npw) - xp(npw) ) / ( xp(npw) - xp(npww) ) * ( p(npw) - p(npww) )

      end do

      ! South boundary
      j = 1
      do i = 2, nx-1
         np  = nx * (j-1) + i
         npn  = np + nx

         p(np) = p(npn)
      end do

      ! North boundary
      j = ny
      do i = 2, nx-1
         np  = nx * (j-1) + i
         nps  = np - nx

         p(np) = p(nps)
      end do

      ! SW corner
      i = 1
      j = 1
      np   = nx * (j-1) + i
      npn  = np + nx

      p(np) = p(npn)

      ! SE corner
      i = nx
      j = 1
      np   = nx * (j-1) + i
      npn  = np + nx

      p(np) = p(npn)

      ! NW corner
      i = 1
      j = ny
      np   = nx * (j-1) + i
      nps  = np - nx

      p(np) = p(nps)

      ! NE corner
      i = nx
      j = ny
      np   = nx * (j-1) + i
      nps  = np - nx

      p(np) = p(nps)

   end subroutine get_plin_and_p_fictitious


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


   !> \brief Calculates initial conditions for internal flow (interface to get_initial_guess)
   subroutine intflow_initial_conditions( nx, ny, modvis, beta & ! Input
      ,      thermomodel, xe, ye, xk, yk, radius, rn, x, y, xp & ! Input
      ,                                                  iflow & ! InOutput
      ,                   p, T, u, v, ue, un, ve, vn, Uce, Vcn & ! Output
      ,                       ro, roe, ron, Tbn, Tbs, Tbe, Tbw ) ! Output
      integer, intent(in) :: nx     !< Number of volumes in csi direction (real+fictitious)
      integer, intent(in) :: ny     !< Number of volumes in eta direction (real+fictitious)
      integer, intent(in) :: modvis !< modvis = 0 -> Euler;  modvis = 1 -> Navier-Stokes
      real(8), intent(in) :: beta   !< Constant of the UDS/CDS mixing scheme

      class(class_thermophysical_abstract), pointer, intent(in)  :: thermomodel !< A pointer to the thermophysical model

      real(8), dimension(nx*ny), intent(in) :: xe     ! x_eta at center of face east
      real(8), dimension(nx*ny), intent(in) :: ye     !< y_eta at face east of volume P
      real(8), dimension(nx*ny), intent(in) :: xk     !< x_csi at center of face north
      real(8), dimension(nx*ny), intent(in) :: yk     !< y_csi at face north of volume P
      real(8), dimension(nx*ny), intent(in) :: radius !< Radius of northest corner of volume P
      real(8), dimension(nx*ny), intent(in) :: rn     !< Radius of the center of north face of volume P
      real(8), dimension(nx*ny), intent(in) :: x      !< Coord. x of the northest corner of volume P
      real(8), dimension(nx*ny), intent(in) :: y      !< Coord. x of the northest corner of volume P
      real(8), dimension(nx*ny), intent(in) :: xp     !< Coord. x of the centroid of volume P

      type(type_intflow),     intent(inout) :: iflow  !< Data for internal flow

      real(8), dimension(nx*ny), intent(out) :: p   !< Pressure at center o volume P
      real(8), dimension(nx*ny), intent(out) :: T   !< Temperature of the last iteraction
      real(8), dimension(nx*ny), intent(out) :: u   !< Cartesian velocity of the last iteraction
      real(8), dimension(nx*ny), intent(out) :: v   !< Cartesian velocity of the last iteraction

      real(8), dimension(nx*ny), intent(out) :: ue  !< Cartesian velocity u at center of east face
      real(8), dimension(nx*ny), intent(out) :: un  !< Cartesian velocity u at center of north face
      real(8), dimension(nx*ny), intent(out) :: ve  !< Cartesian velocity v at center of east face (m/s)
      real(8), dimension(nx*ny), intent(out) :: vn  !< Cartesian velocity v at center of north face (m/s)
      real(8), dimension(nx*ny), intent(out) :: Uce !< Contravariant velocity U at east face
      real(8), dimension(nx*ny), intent(out) :: Vcn !< Contravariant velocity V at north face

      real(8), dimension(nx*ny), intent(out) :: ro  !< Specific mass (absolute density) at center of volumes
      real(8), dimension(nx*ny), intent(out) :: roe !< Absolute density at east face
      real(8), dimension(nx*ny), intent(out) :: ron !< Absolute density at north face
      real(8), dimension(nx),    intent(out) :: Tbn !< Temperature over the north boundary (K)
      real(8), dimension(nx),    intent(out) :: Tbs !< Temperature over the south boundary (K)
      real(8), dimension(ny),    intent(out) :: Tbe !< Temperature over the  east boundary (K)
      real(8), dimension(ny),    intent(out) :: Tbw !< Temperature over the  west boundary (K)
      real(8) :: gamma, Rg

      ve = 0.d0
      vn = 0.d0

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

         call get_initial_guess( nx, ny, modvis, beta, po, T0, gamma, Rg  & ! Input
         ,                       Sg, xe, ye, xk, yk, radius, rn, x, y, xp & ! Input
         ,                       M1D, p1D, T1D, u1D, p, T, u, v, ue, ve   & ! Output
         ,                       un, vn, Uce, Vcn, uin, vin, pin, Tin, Mw & ! Output
         ,                               fm1D, Fd1D, Fpv1D, ro, roe, ron  ) ! Output

         Tbw = Tin
         Tbs = T1D
         Tbn = T1D
         Tbe = T1D(nx-1)

      end associate

   end subroutine


   !> \brief This subroutine is a legacy code from 5.8.1. It must be reviewed,
   !! because it assumes that eta lines are always in the y direction.
   subroutine get_initial_guess( nx, ny, modvis, beta, po, T0, gamma, Rg  & ! Input
      ,                        Sg, xe, ye, xk, yk, radius, rn, x, y, xp & ! Input
      ,                        M1D, p1D, T1D, u1D, p, T, u, v, ue, ve   & ! Output
      ,                        un, vn, Uce, Vcn, uin, vin, pin, Tin, Mw & ! Output
      ,                                 fm1D, Fd1D, Fpv1D, ro, roe, ron ) ! Output
      implicit none
      integer, intent(in) :: nx     ! Number of volumes in csi direction (real+fictitious)
      integer, intent(in) :: ny     ! Number of volumes in eta direction (real+fictitious)
      integer, intent(in) :: modvis ! modvis = 0 -> Euler;  modvis = 1 -> Navier-Stokes
      real(8), intent(in) :: beta   ! Constant of the UDS/CDS mixing scheme
      real(8), intent(in) :: po     ! Stagnation pressure in the chamber
      real(8), intent(in) :: T0     ! Stagnation temperature in the chamber
      real(8), intent(in) :: gamma  ! gamma = Cpo / Cvo in the chamber
      real(8), intent(in) :: Rg     ! Perfect gas constant
      real(8), intent(in) :: Sg     ! Throttle area
      real(8), dimension(nx*ny), intent(in) :: xe     ! x_eta at center of face east
      real(8), dimension(nx*ny), intent(in) :: ye     ! y_eta at face east of volume P
      real(8), dimension(nx*ny), intent(in) :: xk     !< x_csi at center of face north
      real(8), dimension(nx*ny), intent(in) :: yk     !< y_csi at center of face north
      real(8), dimension(nx*ny), intent(in) :: radius ! Radius of northest corner of volume P
      real(8), dimension(nx*ny), intent(in) :: rn     ! Radius of the center of north face of volume P
      real(8), dimension(nx*ny), intent(in) :: x      ! Coord. x of the northest corner of volume P
      real(8), dimension(nx*ny), intent(in) :: y      ! Coord. y of the northest corner of volume P
      real(8), dimension(nx*ny), intent(in) :: xp     ! Coord. x of the centroid of volume P

      real(8), intent(out) :: fm1D  !
      real(8), intent(out) :: Fd1D  !
      real(8), intent(out) :: Fpv1D !

      real(8), dimension(nx),    intent(out) :: M1D !  Mach number. Isentropic flow
      real(8), dimension(nx),    intent(out) :: p1D !  Pressure. Isentropic flow
      real(8), dimension(nx),    intent(out) :: T1D !  Temperature. Isentropic flow
      real(8), dimension(nx),    intent(out) :: u1D !  Velocity. Isentropic flow

      real(8), dimension(nx*ny), intent(out) :: p   ! Pressure at center o volume P
      real(8), dimension(nx*ny), intent(out) :: T   ! Temperature of the last iteraction
      real(8), dimension(nx*ny), intent(out) :: u   ! Cartesian velocity of the last iteraction
      real(8), dimension(nx*ny), intent(out) :: v   ! Cartesian velocity of the last iteraction

      real(8), dimension(nx*ny), intent(out) :: ue  ! Cartesian velocity u at center of east face
      real(8), dimension(nx*ny), intent(out) :: ve  !< Cartesian velocity v at center of east face (m/s)
      real(8), dimension(nx*ny), intent(out) :: un  ! Cartesian velocity u at center of north face
      real(8), dimension(nx*ny), intent(out) :: vn  !< Cartesian velocity v at center of north face (m/s)
      real(8), dimension(nx*ny), intent(out) :: Uce ! Contravariant velocity U at east face
      real(8), dimension(nx*ny), intent(out) :: Vcn ! Contravariant velocity V at north face

      real(8), dimension (ny),   intent(out) :: uin ! Velocity u in the entrance
      real(8), dimension (ny),   intent(out) :: vin ! Velocity v in the entrance
      real(8), dimension (ny),   intent(out) :: pin ! Pressure in the entrance
      real(8), dimension (ny),   intent(out) :: Tin ! Temperature in the entrance
      real(8), dimension (ny),   intent(out) :: Mw  ! Mach number in the entrance

      real(8), dimension(nx*ny), intent(out) :: ro  ! Specific mass (absolute density) at center of volumes
      real(8), dimension(nx*ny), intent(out) :: roe ! Absolute density at east face
      real(8), dimension(nx*ny), intent(out) :: ron ! Absolute density at north face

      ! Auxiliary variables
      real(8), parameter :: pi = acos(-1.d0)

      integer :: i, j, np, npe, npw, npww, npn, nps, kf
      real(8) :: Si, Me1D, pe1D, Te1D, ue1D, aux
      integer :: ig

      call get_isentropic_mass_flow( po, T0, gamma, Rg, Sg , fm1D)

      ! i=ig when the throat cross section is in the east face of the CV
      ig = get_ig(nx, ny, y)

      ! Calculating isentropic solution

      do i = 2, nx-1

         j = ny-1

         np = nx * (j-1) + i

         Si = pi * rn(np) ** 2

         kf = 0

         if ( i > ig ) kf = 1

         call get_mach_area( kf, Si/Sg, gamma, M1D(i) )

         aux = 1.d0 + 0.5d0 * (gamma-1.d0) * M1D(i)**2

         p1D(i) = po * aux ** (-gamma/(gamma-1.d0))

         T1D(i) = T0 / aux

         u1D(i) = M1D(i) * sqrt( gamma * Rg * T1D(i) )

         do j = 1, ny

            np = nx * (j-1) + i

            p(np) = p1D(i)

            T(np) = T1D(i)

            u(np) = u1D(i)

            if ( j > 1 .and. j < ny-1 ) then

               un(np) = u1D(i)

               Vcn(np) = -un(np) * yk(np) ! vn = 0

            end if

         end do

      end do

      ! Calculates u and v at fictitious nodes using boundary conditions
      call intflow_extrapolate_u_v_to_fictitious( nx, ny, modvis, x, xp, xk, yk, u, v) ! InOutput: last two entries

      ! Calculates gas proprieties in the entrance
      call get_uin_vin_pin_Tin_Mw( nx, ny, gamma, Rg, po, T0, u & ! Input
      ,                             uin, vin, pin, Tin, Mw)  ! Output

      ! Calculates p and T in the fictitious nodes of the entrance according to boundary conditions
      i = 1
      do j = 2, ny-1

         np   = nx * (j-1) + i

         npe  = np + 1

         p(np) = 2.d0 * pin(j) - p(npe)

         T(np) = 2.d0 * Tin(j) - T(npe)

      end do

      ! Calculates p and T in the fictitious nodes of the exit according to boundary conditions
      i = nx
      do j = 2, ny-1

         np   = nx * (j-1) + i
         npw  = np - 1
         npww = np - 2

         p(np) = p(npw) + 2.d0 * ( x(npw) - xp(npw) ) / ( xp(npw) - xp(npww) ) * ( p(npw)-p(npww) )
         T(np) = T(npw) + 2.d0 * ( x(npw) - xp(npw) ) / ( xp(npw) - xp(npww) ) * ( T(npw)-T(npww) )

      end do

      ! Calculating p and T in the SW corner
      i = 1
      j = 1
      np   = nx * (j-1) + i
      npn  = np + nx
      p(np) = p(npn)
      T(np) = T(npn)

      ! Calculating p and T in the SE corner
      i = nx
      j = 1
      np   = nx * (j-1) + i
      npn  = np + nx
      p(np) = p(npn)
      T(np) = T(npn)

      ! Calculating p and T in the NW corner
      i = 1
      j = ny
      np   = nx * (j-1) + i
      nps  = np - nx
      p(np) = p(nps)
      T(np) = T(nps)

      ! Calculating p and T in the NE corner
      i = nx
      j = ny
      np   = nx * (j-1) + i
      nps  = np - nx
      p(np) = p(nps)
      T(np) = T(nps)


      ! Initial guess of velocities in east face of each CV
      ! and analytical solution in the entrance and exit boundaries
      do i = 1, nx-1

         j = ny-1

         np = nx * (j-1) + i

         Si = pi * radius(np) ** 2

         kf = 0

         if ( i > ig ) kf = 1

         call get_mach_area( kf, Si/Sg, gamma, Me1D )

         aux = 1.d0 + 0.5d0 * (gamma-1.d0) * Me1D**2

         pe1D = po * aux ** (-gamma/(gamma-1.d0))

         Te1D = T0 / aux

         ue1D = Me1D * sqrt( gamma * Rg * Te1D )

         do j = 2, ny-1

            np = nx * (j-1) + i

            ue(np) = ue1D

            Uce(np) = ue1D * ye(np) ! v = 0

         end do

         ! Entrance boundary
         if ( i == 1 ) then
            u1D(i) = ue1D
            p1D(i) = pe1D
            T1D(i) = Te1D
            M1D(i) = Me1D
         end if

         ! Exit boundary
         if ( i == nx-1 ) then
            u1D(nx) = ue1D
            p1D(nx) = pe1D
            T1D(nx) = Te1D
            M1D(nx) = Me1D

            Fd1D = fm1D * ue1D
            Fpv1D = pe1D * Si
         end if

      end do

      ! Calculates velocities on boundaries according bc
      call intflow_velocities_at_boundary_faces( nx, ny, modvis & ! Input
      ,                          x, xp, xe, ye, xk, yk, u, v & ! Input
      ,                             ue, ve, un, vn, Uce, Vcn ) ! Input and Output


      ! Specific mass at nodes
      call get_density_at_nodes( nx, ny, Rg, p, T, ro) ! ro is output

      ! Specific mass at faces
      call get_density_at_faces( nx, ny, beta, ro, Uce, Vcn, roe, ron) ! roe and ron are output

   end subroutine get_initial_guess


   !> \brief Calculates the mass flow according to the Q1D theory
   subroutine get_isentropic_mass_flow( po, T0, gamma, Rg, Sg , fm1D)
      implicit none
      real(8), intent(in)  :: po    !< Stagnation pressure
      real(8), intent(in)  :: T0    !< Stagnation temperature
      real(8), intent(in)  :: gamma !< Specific heat ratio
      real(8), intent(in)  :: Rg    !< Perfect gas constant
      real(8), intent(in)  :: Sg    !< Throttle area
      real(8), intent(out) :: fm1D  !< Mass flow rate

      fm1D = po * Sg * sqrt( gamma/ (Rg*T0) * (2.d0/(gamma+1.d0)) ** ((gamma+1.d0)/(gamma-1.d0)) )

   end subroutine get_isentropic_mass_flow


   !> \brief Calculates the Area-Mach relation according to the Q1D theory
   subroutine get_mach_area( kf, ar, gamma, M )
      implicit none
      integer, intent(in)  :: kf    !< Kind of flow: 0 = subsonic, 1 = supersonic
      real(8), intent(in)  :: ar    !< Area ratio = local area / throttle area
      real(8), intent(in)  :: gamma !< Specific heat ratio
      real(8), intent(out) :: M     !< Mach number at the local area

      ! Auxiliary variables

      integer :: it

      real(8) :: M0, M1, f, x

      f(x) = ( 2.d0 * ( 1.d0 + (gamma-1.d0) * x**2 / 2.d0 ) / (gamma+1.d0) ) &
         ** ((gamma+1.d0)/(gamma-1.d0)) - ( x * ar ) ** 2

      if ( abs(ar-1.d0) < 1.d-15 ) then ! Sonic flow in the throttle

         M = 1.d0

      else

         if ( kf == 0 ) then ! Subsonic
            M0 = 1.d-4
            M1 = 1.d0
         else                ! Supersonic
            M0 = 1.d0
            M1 = 100.d0
         end if

         do it = 1, 1000

            M = 0.5d0 * ( M0 + M1 )

            if ( f(M) * f(M0) > 0 ) then

               M0 = M

            else

               M1 = M

            end if

            if ( abs(M1-M0) < 1.d-15 ) exit

         end do

      end if

   end subroutine get_mach_area


   !> \brief i=ig when the throat cross section is in the east face of the CV
   integer function get_ig(nx, ny, y)
      implicit none
      integer,                   intent(in) :: nx !< Number of volumes in csi direction (real+fictitious)
      integer,                   intent(in) :: ny !< Number of volumes in eta direction (real+fictitious)
      real(8), dimension(nx*ny), intent(in) :: y  !< Coord. y of the northest corner of volume P

      integer i, j, np

      get_ig = 1

      j = ny-1

      do i = 2, nx-1

         np = nx * (j-1) + i

         if ( y(np-1) < y(np)  ) then

            get_ig = i-1

            return

         end if

      end do

   end function


   !> \brief Calculates the mass flow rate and the thrust
   subroutine get_mass_flow_rate_and_thrust( nx, ny, re, roe, u, Uce, fmi, fme, Fd)
      implicit none
      integer, intent(in) :: nx  !< Number of volumes in csi direction (real+fictitious)
      integer, intent(in) :: ny  !< Number of volumes in eta direction (real+fictitious)
      real(8), dimension (nx*ny), intent(in) :: re   !< Radius of the center of east face of volume P
      real(8), dimension (nx*ny), intent(in) :: roe  !< Absolute density at east face
      real(8), dimension (nx*ny), intent(in) :: u    !< Cartesian velocity of the last iteraction
      real(8), dimension (nx*ny), intent(in) :: Uce  !< Contravariant velocity U at east face
      real(8), intent(out) :: fmi !< Mass flow rate at entrance
      real(8), intent(out) :: fme !< Mass flow rate at exit
      real(8), intent(out) :: Fd  !< Thrust

      real(8), parameter :: pi = acos(-1.d0)

      integer :: i, j, np

      fmi = 0.d0

      i = 1

      do j = 2, ny-1

         np = nx * (j-1) + i

         fmi = fmi + re(np) * roe(np) * Uce(np)

      end do

      fmi = fmi * 2.d0 * pi



      fme = 0.d0

      Fd = 0.d0

      i = nx-1

      do j = 2, ny-1

         np = nx * (j-1) + i

         fme = fme + re(np) * roe(np) * Uce(np)

         Fd = Fd + re(np) * roe(np) * Uce(np) * ( u(np) + u(np+1) ) / 2.d0

      end do

      fme = fme * 2.d0 * pi

      Fd = Fd * 2.d0 * pi

   end subroutine get_mass_flow_rate_and_thrust


   !> \brief Calculates the main variables of the internal flow
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

      write(msg(1), "(3(1X,A23))")      "fmi", "fme/fmi", "Fd"
      write(msg(2), "(3(1X,ES23.16))")   fmi,   fme/fmi,   Fd

   end subroutine


end module
