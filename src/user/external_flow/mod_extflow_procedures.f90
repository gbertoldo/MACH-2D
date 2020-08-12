!>
!! \brief Contains procedures for external flow calculations.
!!
module mod_extflow_procedures

   use thermophysical ! Thermophysical module
   use bc9d           ! Boundary conditions for 9-diagonal linear systems
   use bc5d           ! Boundary conditions for 5-diagonal linear systems

   implicit none

contains


   !> \brief Calculates some of the free stream properties of the gas
   subroutine get_free_stream_properties(thermomodel, PF, TF, MF, CPF, GF, VLF, KPF, PRF, ROF, UF, REFm, HF) ! Output: last 9
      implicit none
      class(class_thermophysical_abstract), pointer, intent(in) :: thermomodel !< A pointer to the thermophysical model
      real(8), intent(in)  ::   PF !< Free stream pressure (Pa)
      real(8), intent(in)  ::   TF !< Free stream temperature (K)
      real(8), intent(in)  ::   MF !< Free stream Mach number
      real(8), intent(out) ::  CPF !< Free stream Cp (J/kg.K)
      real(8), intent(out) ::   GF !< Free stream gamma
      real(8), intent(out) ::  VLF !< Free stream viscosity (Pa.s)
      real(8), intent(out) ::  KPF !< Free stream thermal conductivity (W/m.K)
      real(8), intent(out) ::  PRF !< Free stream Prandtl number
      real(8), intent(out) ::  ROF !< Free stream density (kg/m3)
      real(8), intent(out) ::   UF !< Free stream speed (m/s)
      real(8), intent(out) :: REFm !< Free stream Reynolds number per meter (1/m)
      real(8), intent(out) ::   HF !< Free stream total enthalpy (m2/s2)

      ! Free stream Cp
      CPF = thermomodel%cp(TF)

      ! Free stream gamma
      GF = CPF / ( CPF - thermomodel%Rg )

      ! Free stream viscosity (Pa.s)
      VLF = thermomodel%mu(TF)

      ! Free stream thermal conductivity (W/m.K)
      KPF = thermomodel%kp(TF)

      ! Free stream Prandtl number
      PRF = CPF * VLF / KPF

      ! Free stream density (kg/m3)
      ROF = PF / ( thermomodel%Rg * TF )

      ! Free stream speed (m/s)
      UF = sqrt( GF * thermomodel%Rg * TF ) * MF

      ! Free stream Reynolds number per meter
      REFm = UF * ROF / VLF

      ! Free stream total enthalpy (m2/s2)
      HF = CPF * TF + UF * UF / 2.d0

   end subroutine get_free_stream_properties


   !> \brief Estimates the boundary layer width and the width of the volume
   !! closer to the wall
   subroutine get_boundary_layer_width_estimate(lr, REFm, cbl, wbl, a1) ! Output: last two
      implicit none
      real(8), intent(in)  ::  lr  !< Body length (m)
      real(8), intent(in)  :: REFm !< Free stream Reynolds number per meter (1/m)
      real(8), intent(in)  ::  cbl !< The width of the vol. closer to the wall is 'cbl' times the width of the b. layer
      real(8), intent(out) ::  wbl !< Estimated width of the boundary layer (m)
      real(8), intent(out) ::   a1 !< Width of the volume closer to the wall (m)

      wbl = lr / sqrt( REFm * lr )

      a1 = cbl * wbl

   end subroutine get_boundary_layer_width_estimate


   ! Calculates the SIMPLEC coefficients at the boundary faces
   subroutine get_boundary_simplec_coefficients( nx, ny, modvis, due, dve, dun &
         , dvn, de, dn) ! InOutput: last 6
      implicit none
      integer, intent(in) :: nx  ! Number of volumes in csi direction (real+fictitious)
      integer, intent(in) :: ny  ! Number of volumes in eta direction (real+fictitious)
      integer, intent(in) :: modvis ! Viscosity model (0=Euler, 1=NS)
      real(8), dimension(nx*ny), intent(inout) :: due  ! Simplec coef. for the cartesian velocity u (east face)
      real(8), dimension(nx*ny), intent(inout) :: dve  ! Simplec coef. for the cartesian velocity v (east face)
      real(8), dimension(nx*ny), intent(inout) :: dun  ! Simplec coef. for the cartesian velocity u (north face)
      real(8), dimension(nx*ny), intent(inout) :: dvn  ! Simplec coef. for the cartesian velocity v (north face)
      real(8), dimension(nx*ny), intent(inout) :: de   ! Simplec coef. for the contravariant velocity U (east face)
      real(8), dimension(nx*ny), intent(inout) :: dn   ! Simplec coef. for the contravariant velocity V (north face)

      ! Auxiliary variables
      integer :: i, j, np, npe, npw, npn

      ! West boundary

      i = 1

      do j = 2, ny-1

         np   = nx * (j-1) + i

         npe  = np + 1

         due(np) = due(npe)

         dve(np) = 0.d0

         de(np) = 0.d0

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

      if ( modvis == 0 ) then ! Euler

         j = 1

         do i = 2, nx-1

            np   = nx * (j-1) + i

            npn  = np + nx

            dun(np) = dun(npn)

            dvn(np) = dvn(npn)

            dn(np) = 0.d0

         end do

      else ! NS

         j = 1

         do i = 2, nx-1

            np   = nx * (j-1) + i

            dun(np) = 0.d0

            dvn(np) = 0.d0

            dn(np) = 0.d0

         end do

      end if

      ! North boundary

      j = ny-1

      do i = 2, nx-1

         np   = nx * (j-1) + i

         dun(np) = 0.d0

         dvn(np) = 0.d0

         dn(np) = 0.d0

      end do

   end subroutine get_boundary_simplec_coefficients


   ! Guess the initial conditions
   subroutine get_initial_conditions( nx, ny, itemax, modvis, PF, TF, Rg, GF &
         , MF, UF, xe, ye, xk, yk, xke, yke, alphae, betae, p, T, ro, roe, ron, u, v   &
         , ue, ve, un, vn, Uce, Vcn, Ucbe, Vcbe, Tbn, Tbs, Tbe, Tbw ) ! UF and last 19 are output

      implicit none
      integer, intent(in)  :: nx  ! Number of volumes in csi direction (real+fictitious)
      integer, intent(in)  :: itemax !< Number of iteractions for extrapolation to fictitious
      integer, intent(in)  :: ny  ! Number of volumes in eta direction (real+fictitious)
      integer, intent(in)  :: modvis ! Viscosity model (0=Euler, 1=NS)
      real(8), intent(in)  :: PF  ! Far field pressure (Pa)
      real(8), intent(in)  :: TF  ! Far field temperature (K)
      real(8), intent(in)  :: Rg  ! Perfect gas constant (J/kg.K)
      real(8), intent(in)  :: GF  ! GF = gamma = Cp / Cv (for the free stream)
      real(8), intent(in)  :: MF  ! Far field Mach number
      real(8), intent(out) :: UF  ! Far field gas speed (m/s)
      real(8), dimension(nx*ny), intent(in)  :: xe     ! x_eta at face east of volume P (m)
      real(8), dimension(nx*ny), intent(in)  :: ye     ! y_eta at face east of volume P (m)
      real(8), dimension(nx*ny), intent(in)  :: xk     ! x_csi at face north of volume P (m)
      real(8), dimension(nx*ny), intent(in)  :: yk     ! y_csi at face north of volume P (m)
      real(8), dimension(nx*ny), intent(in)  :: xke    !< x_csi at center of face east
      real(8), dimension(nx*ny), intent(in)  :: yke    !< y_csi at center of face east
      real(8), dimension(nx*ny), intent(in)  :: alphae !< (metric) alpha at the center of east face of volume P
      real(8), dimension(nx*ny), intent(in)  :: betae  !< (metric) beta  at the center of east face of volume P
      real(8), dimension(nx*ny), intent(out) :: p      ! Pressure at center o volume P (Pa)
      real(8), dimension(nx*ny), intent(out) :: T      ! Temperature of the last iteraction (K)
      real(8), dimension(nx*ny), intent(out) :: ro     ! Specific mass (absolute density) at center of volumes (kg/m3)
      real(8), dimension(nx*ny), intent(out) :: roe    ! Absolute density at east face (kg/m3)
      real(8), dimension(nx*ny), intent(out) :: ron    ! Absolute density at north face (kg/m3)
      real(8), dimension(nx*ny), intent(out) :: u      ! Cartesian velocity of the last iteraction (m/s)
      real(8), dimension(nx*ny), intent(out) :: v      ! Cartesian velocity of the last iteraction (m/s)
      real(8), dimension(nx*ny), intent(out) :: ue     ! Cartesian velocity u at center of east face (m/s)
      real(8), dimension(nx*ny), intent(out) :: ve     ! Cartesian velocity v at center of east face (m/s)
      real(8), dimension(nx*ny), intent(out) :: un     ! Cartesian velocity u at center of north face (m/s)
      real(8), dimension(nx*ny), intent(out) :: vn     ! Cartesian velocity v at center of north face (m/s)
      real(8), dimension(nx*ny), intent(out) :: Uce    ! Contravariant velocity U at east face (m2/s)
      real(8), dimension(nx*ny), intent(out) :: Vcn    ! Contravariant velocity V at north face (m2/s)
      real(8), dimension(ny),    intent(out) :: Ucbe   !< Uc over the faces of the east  boundary (m2/s)
      real(8), dimension(ny),    intent(out) :: Vcbe   !< Vc over the faces of the east  boundary (m2/s)
      real(8), dimension(nx),    intent(out) :: Tbn    !< Temperature over the north boundary (K)
      real(8), dimension(nx),    intent(out) :: Tbs    !< Temperature over the south boundary (K)
      real(8), dimension(ny),    intent(out) :: Tbe    !< Temperature over the  east boundary (K)
      real(8), dimension(ny),    intent(out) :: Tbw    !< Temperature over the  west boundary (K)



      ! Auxiliary variables

      integer :: i, j, np

      real(8), dimension(nx,9) :: a9bn ! Coefficients of the discretization for the north boundary
      real(8), dimension(nx,9) :: a9bs ! Coefficients of the discretization for the south boundary
      real(8), dimension(ny,9) :: a9be ! Coefficients of the discretization for the east boundary
      real(8), dimension(ny,9) :: a9bw ! Coefficients of the discretization for the west boundary
      real(8), dimension(nx)   :: b9bn ! Source of the discretization for the north boundary
      real(8), dimension(nx)   :: b9bs ! Source of the discretization for the south boundar
      real(8), dimension(ny)   :: b9be ! Source of the discretization for the east boundary
      real(8), dimension(ny)   :: b9bw ! Source of the discretization for the west boundary




      UF = MF * sqrt( GF * Rg * TF )

      p = PF

      T = TF

      ro = PF / ( Rg * TF )

      roe = ro

      ron = ro

      u = UF

      v = 0.d0

      ! Defines the temperature over the boundaries of calculation

      Tbn = TF
      Tbs = TF
      Tbe = TF
      Tbw = TF


      ! Calculates the contravariant velocities over the east boundary
      call get_Ucbe_Vcbe(nx, ny, xe, ye, xke, yke, u, v, Ucbe, Vcbe) ! Output: last two

      ! ========================================================================

      !                      USER DEPENDENT SUBROUTINE

      call get_bc_scheme_u(nx, ny, modvis, UF, xk, yk, alphae, betae, u, v &
         , Ucbe, Vcbe, a9bn, a9bs, a9be, a9bw, b9bn, b9bs, b9be, b9bw) ! Output: last eight

      ! ========================================================================


      ! Extrapolates u from the real volumes to the fictitious ones according to the boundary conditions
      call get_bc_extrapolation_9d(nx, ny, itemax, a9bn, a9bs, a9be, a9bw, b9bn, b9bs &
         , b9be, b9bw, u) ! InOutput: last one


      ! ========================================================================

      !                      USER DEPENDENT SUBROUTINE

      call get_bc_scheme_v(nx, ny, modvis, xk, yk, u, v, Ucbe, Vcbe &
         , a9bn, a9bs, a9be, a9bw, b9bn, b9bs, b9be, b9bw) ! Output: last eight

      ! ========================================================================


      ! Extrapolates v from the real volumes to the fictitious ones according to the boundary conditions
      call get_bc_extrapolation_9d(nx, ny, itemax, a9bn, a9bs, a9be, a9bw, b9bn, b9bs &
         , b9be, b9bw, v) ! InOutput: last one



      ! u, v, U and V on the internal real faces

      do i = 2, nx-2
         do j = 2, ny-1

            np = nx * (j-1) + i

            ue(np) = u(np)

            ve(np) = v(np)

            Uce(np) = ue(np) * ye(np) - ve(np) * xe(np)

         end do
      end do

      do i = 2, nx-1
         do j = 2, ny-2

            np = nx * (j-1) + i

            un(np) = u(np)

            vn(np) = v(np)

            Vcn(np) = vn(np) * xk(np) - un(np) * yk(np)

         end do
      end do

      ! u, v, U and V  on the boundary real faces

      call get_velocities_at_boundary_faces( nx, ny, xe, ye, xk, yk &

      , u, v, ue, ve, un, vn, Uce, Vcn) ! Last six are output

   end subroutine get_initial_conditions


   !> \brief Defines the numerical scheme for the boundary conditions of pl
   subroutine get_bc_scheme_pl(nx, ny, a5bn, a5bs, a5be, a5bw, b5bn, b5bs, b5be, b5bw) ! Output: last eight
      implicit none
      integer, intent(in) :: nx !< Number of volumes in csi direction (real+fictitious)
      integer, intent(in) :: ny !< Number of volumes in eta direction (real+fictitious)
      real(8), dimension(nx,5), intent(out) :: a5bn !< Coefficients of the discretization for the north boundary
      real(8), dimension(nx,5), intent(out) :: a5bs !< Coefficients of the discretization for the south boundary
      real(8), dimension(ny,5), intent(out) :: a5be !< Coefficients of the discretization for the east boundary
      real(8), dimension(ny,5), intent(out) :: a5bw !< Coefficients of the discretization for the west boundary
      real(8), dimension(nx),   intent(out) :: b5bn !< Source of the discretization for the north boundary
      real(8), dimension(nx),   intent(out) :: b5bs !< Source of the discretization for the south boundar
      real(8), dimension(ny),   intent(out) :: b5be !< Source of the discretization for the east boundary
      real(8), dimension(ny),   intent(out) :: b5bw !< Source of the discretization for the west boundary

      ! Inner variables

      !real(8), dimension(nx) :: Fbn  !< Field F over the north boundary



      ! North boundary

      !Fbn = 0.d0

      !call get_bc_dirichlet_north_5d(nx, Fbn, a5bn, b5bn) ! Output: last two

      ! UDS approximation
      a5bn = 0.d0
      a5bn(:,3) = 1.d0
      b5bn = 0.d0



      ! South boundary

      call get_ibc_null_normal_grad_south_5d(nx, a5bs, b5bs) ! Output: last two



      ! East boundary

      call get_ibc_null_normal_grad_east_5d(ny, a5be, b5be) ! Output: last two



      ! West boundary

      call get_ibc_null_normal_grad_west_5d(ny, a5bw, b5bw) ! Output: last two



      ! Defines the numerical scheme for the calculation of the boundary
      ! conditions at the fictitious volumes at corners of the transformed domain
      call get_bc_corners_5d(nx, ny, a5bn, a5bs, a5be, a5bw, b5bn, b5bs &
      , b5be, b5bw) ! Output: last eight

   end subroutine get_bc_scheme_pl


   !> \brief Defines the numerical scheme of the boundary conditions for u
   subroutine get_bc_scheme_u(nx, ny, modvis, UF, xk, yk, alphae, betae, u, v &
      , Ucbe, Vcbe, a9bn, a9bs, a9be, a9bw, b9bn, b9bs, b9be, b9bw) ! Output: last eight
      implicit none
      integer, intent(in) :: nx     !< Number of volumes in csi direction (real+fictitious)
      integer, intent(in) :: ny     !< Number of volumes in eta direction (real+fictitious)
      integer, intent(in) :: modvis !< Viscosity model (0=Euler, 1=NS)
      real(8), intent(in) :: UF     !< Free-stream speed (m/s)
      real(8), dimension(nx*ny), intent(in)  :: xk   !< x_csi at face north of volume P (m)
      real(8), dimension(nx*ny), intent(in)  :: yk   !< y_csi at face north of volume P (m)
      real(8), dimension(nx*ny), intent(in)  :: alphae !< (metric) alpha at the center of east face of volume P
      real(8), dimension(nx*ny), intent(in)  :: betae  !< (metric) beta  at the center of east face of volume P
      real(8), dimension(nx*ny), intent(in)  :: u    !< x cartesian velocity (m/s)
      real(8), dimension(nx*ny), intent(in)  :: v    !< y cartesian velocity (m/s)
      real(8), dimension(ny),    intent(in)  :: Ucbe !< Uc over the faces of the east  boundary (m2/s)
      real(8), dimension(ny),    intent(in)  :: Vcbe !< Vc over the faces of the east  boundary (m2/s)
      real(8), dimension(nx,9),  intent(out) :: a9bn !< Coefficients of the discretization for the north boundary
      real(8), dimension(nx,9),  intent(out) :: a9bs !< Coefficients of the discretization for the south boundary
      real(8), dimension(ny,9),  intent(out) :: a9be !< Coefficients of the discretization for the east boundary
      real(8), dimension(ny,9),  intent(out) :: a9bw !< Coefficients of the discretization for the west boundary
      real(8), dimension(nx),    intent(out) :: b9bn !< Source of the discretization for the north boundary
      real(8), dimension(nx),    intent(out) :: b9bs !< Source of the discretization for the south boundar
      real(8), dimension(ny),    intent(out) :: b9be !< Source of the discretization for the east boundary
      real(8), dimension(ny),    intent(out) :: b9bw !< Source of the discretization for the west boundary

      ! Inner variables

      real(8), dimension(nx) :: Fbns  !< Field F over the north or south boundary



      ! North boundary

      !Fbns = UF

      !call get_bc_dirichlet_north_9d(nx, Fbns, a9bn, b9bn) ! Output: last two

      ! UDS approximation
      a9bn = 0.d0
      a9bn(:,5) = 1.d0
      b9bn = UF



      ! South boundary

      if ( modvis == 0 ) then ! Euler model

         call get_bc_u_slip_south_9d(nx, ny, xk, yk, u, v, a9bs, b9bs) ! Output: last two

      else ! Navier-Stokes model

         Fbns = 0.d0

         call get_bc_dirichlet_south_9d(nx, Fbns, a9bs, b9bs) ! Output: last two

      end if



      ! East boundary

      !call get_bc_streamlined_exit_east_9d(ny, Ucbe, Vcbe, a9be, b9be) ! Output: last two
      call get_ibc_null_normal_grad_east_9d(ny, a9be, b9be) ! Output: last two



      ! West boundary

      call get_bc_null_normal_grad_west_9d(nx, ny, alphae, betae, a9bw, b9bw) ! Output: last two


      ! Defines the numerical scheme for the calculation of the boundary
      ! conditions at the fictitious volumes at corners of the transformed domain
      call get_bc_corners_9d(nx, ny, a9bn, a9bs, a9be, a9bw, b9bn, b9bs &
         , b9be, b9bw ) ! Output: last eight

   end subroutine get_bc_scheme_u


   !> \brief Defines the numerical scheme of the boundary conditions for v
   subroutine get_bc_scheme_v(nx, ny, modvis, xk, yk, u, v, Ucbe, Vcbe &
      , a9bn, a9bs, a9be, a9bw, b9bn, b9bs, b9be, b9bw) ! Output: last eight
      implicit none
      integer, intent(in) :: nx     !< Number of volumes in csi direction (real+fictitious)
      integer, intent(in) :: ny     !< Number of volumes in eta direction (real+fictitious)
      integer, intent(in) :: modvis !< Viscosity model (0=Euler, 1=NS)
      real(8), dimension(nx*ny), intent(in)  :: xk   !< x_csi at face north of volume P (m)
      real(8), dimension(nx*ny), intent(in)  :: yk   !< y_csi at face north of volume P (m)
      real(8), dimension(nx*ny), intent(in)  :: u    !< x cartesian velocity (m/s)
      real(8), dimension(nx*ny), intent(in)  :: v    !< y cartesian velocity (m/s)
      real(8), dimension(ny),    intent(in)  :: Ucbe !< Uc over the faces of the east  boundary (m2/s)
      real(8), dimension(ny),    intent(in)  :: Vcbe !< Vc over the faces of the east  boundary (m2/s)
      real(8), dimension(nx,9),  intent(out) :: a9bn !< Coefficients of the discretization for the north boundary
      real(8), dimension(nx,9),  intent(out) :: a9bs !< Coefficients of the discretization for the south boundary
      real(8), dimension(ny,9),  intent(out) :: a9be !< Coefficients of the discretization for the east boundary
      real(8), dimension(ny,9),  intent(out) :: a9bw !< Coefficients of the discretization for the west boundary
      real(8), dimension(nx),    intent(out) :: b9bn !< Source of the discretization for the north boundary
      real(8), dimension(nx),    intent(out) :: b9bs !< Source of the discretization for the south boundar
      real(8), dimension(ny),    intent(out) :: b9be !< Source of the discretization for the east boundary
      real(8), dimension(ny),    intent(out) :: b9bw !< Source of the discretization for the west boundary

      ! Inner variables

      real(8), dimension(nx) :: Fbns  !< Field F over the north or south boundary
      real(8), dimension(ny) :: Fbew  !< Field F over the east or west boundary



      ! North boundary

      !Fbns = 0.d0

      !call get_bc_dirichlet_north_9d(nx, Fbns, a9bn, b9bn) ! Output: last two

      ! UDS approximation
      a9bn = 0.d0
      a9bn(:,5) = 1.d0
      b9bn = 0.d0



      ! South boundary

      if ( modvis == 0 ) then ! Euler model

         call get_bc_v_slip_south_9d(nx, ny, xk, yk, u, v, a9bs, b9bs) ! Output: last two

      else ! Navier-Stokes model

         Fbns = 0.d0

         call get_bc_dirichlet_south_9d(nx, Fbns, a9bs, b9bs) ! Output: last two

      end if



      ! East boundary

      !call get_bc_streamlined_exit_east_9d(ny, Ucbe, Vcbe, a9be, b9be) ! Output: last two
      call get_ibc_null_normal_grad_east_9d(ny, a9be, b9be) ! Output: last two



      ! West boundary

      Fbew = 0.d0

      call get_bc_dirichlet_west_9d(ny, Fbew, a9bw, b9bw) ! Output: last two


      ! Defines the numerical scheme for the calculation of the boundary
      ! conditions at the fictitious volumes at corners of the transformed domain
      call get_bc_corners_9d(nx, ny, a9bn, a9bs, a9be, a9bw, b9bn, b9bs &
         , b9be, b9bw ) ! Output: last eight


   end subroutine get_bc_scheme_v


   !> \brief Defines the numerical scheme for the boundary conditions of ro
   subroutine get_bc_scheme_ro(nx, ny, ROF, alphae, betae, betan, gamman &
         , Ucbe, Vcbe, a9bn, a9bs, a9be, a9bw, b9bn, b9bs, b9be, b9bw) ! Output: last eight
      implicit none
      integer, intent(in) :: nx  !< Number of volumes in csi direction (real+fictitious)
      integer, intent(in) :: ny  !< Number of volumes in eta direction (real+fictitious)
      real(8), intent(in) :: ROF !< Free-stream density (kg/m3)
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

      ! Inner variables

      !real(8), dimension(nx) :: Fbn  !< Field F over the north boundary



      ! North boundary

      !Fbn = ROF

      !call get_bc_dirichlet_north_9d(nx, Fbn, a9bn, b9bn) ! Output: last two

      ! UDS approximation
      a9bn = 0.d0
      a9bn(:,5) = 1.d0
      b9bn = ROF




      ! South boundary

      call get_bc_null_normal_grad_south_9d(nx, ny, betan, gamman, a9bs, b9bs) ! Output: last two



      ! East boundary

      !call get_bc_streamlined_exit_east_9d(ny, Ucbe, Vcbe, a9be, b9be) ! Output: last two
      call get_ibc_null_normal_grad_east_9d(ny, a9be, b9be) ! Output: last two



      ! West boundary

      call get_bc_null_normal_grad_west_9d(nx, ny, alphae, betae, a9bw, b9bw) ! Output: last two


      ! Defines the numerical scheme for the calculation of the boundary
      ! conditions at the fictitious volumes at corners of the transformed domain
      call get_bc_corners_9d(nx, ny, a9bn, a9bs, a9be, a9bw, b9bn, b9bs &
         , b9be, b9bw ) ! Output: last eight


   end subroutine get_bc_scheme_ro


   !> \brief Defines the numerical scheme for the boundary conditions of T
   subroutine get_bc_scheme_T(nx, ny, TF, Tsbc, alphae, betae, betan, gamman &
         , Ucbe, Vcbe, a9bn, a9bs, a9be, a9bw, b9bn, b9bs, b9be, b9bw) ! Output: last eight
      implicit none
      integer, intent(in) :: nx   !< Number of volumes in csi direction (real+fictitious)
      integer, intent(in) :: ny   !< Number of volumes in eta direction (real+fictitious)
      real(8), intent(in) :: TF   !< Free-stream temperature (K)
      real(8), intent(in) :: Tsbc !< Temperature on the south boundary (K) (if negative, adiabatic bc is applied)
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

      ! Inner variables

      real(8), dimension(nx) :: Fbns !< Field F over the north or south boundary



      ! North boundary

      !Fbns = TF

      !call get_bc_dirichlet_north_9d(nx, Fbns, a9bn, b9bn) ! Output: last two

      ! UDS approximation
      a9bn = 0.d0
      a9bn(:,5) = 1.d0
      b9bn = TF



      ! South boundary

      if ( Tsbc < 0.d0 ) then ! Adiabatic boundary condition

         call get_bc_null_normal_grad_south_9d(nx, ny, betan, gamman, a9bs, b9bs) ! Output: last two

      else ! Prescribed temperature

         Fbns = Tsbc

         call get_bc_dirichlet_south_9d(nx, Fbns, a9bs, b9bs) ! Output: last two

      end if



      ! East boundary

      !call get_bc_streamlined_exit_east_9d(ny, Ucbe, Vcbe, a9be, b9be) ! Output: last two
      call get_ibc_null_normal_grad_east_9d(ny, a9be, b9be) ! Output: last two



      ! West boundary

      call get_bc_null_normal_grad_west_9d(nx, ny, alphae, betae, a9bw, b9bw) ! Output: last two



      ! Defines the numerical scheme for the calculation of the boundary
      ! conditions at the fictitious volumes at corners of the transformed domain
      call get_bc_corners_9d(nx, ny, a9bn, a9bs, a9be, a9bw, b9bn, b9bs &
         , b9be, b9bw ) ! Output: last eight


   end subroutine get_bc_scheme_T


   !> \brief Calculates the velocities at boundary faces
   subroutine get_velocities_at_boundary_faces( nx, ny, xe, ye, xk, yk, u, v, ue, ve, un, vn, Uce, Vcn) ! Last six are output
      implicit none
      integer, intent(in) :: nx  ! Number of volumes in csi direction (real+fictitious)
      integer, intent(in) :: ny  ! Number of volumes in eta direction (real+fictitious)
      real(8), dimension(nx*ny), intent(in)  :: xe     ! x_eta at face east of volume P (m)
      real(8), dimension(nx*ny), intent(in)  :: ye     ! y_eta at face east of volume P (m)
      real(8), dimension(nx*ny), intent(in)  :: xk     ! x_csi at face north of volume P (m)
      real(8), dimension(nx*ny), intent(in)  :: yk     ! y_csi at face north of volume P (m)
      real(8), dimension(nx*ny), intent(in)  :: u      ! Cartesian velocity of the last iteraction (m/s)
      real(8), dimension(nx*ny), intent(in)  :: v      ! Cartesian velocity of the last iteraction (m/s)
      real(8), dimension(nx*ny), intent(out) :: ue     ! Cartesian velocity u at center of east face (m/s)
      real(8), dimension(nx*ny), intent(out) :: ve     ! Cartesian velocity v at center of east face (m/s)
      real(8), dimension(nx*ny), intent(out) :: un     ! Cartesian velocity u at center of north face (m/s)
      real(8), dimension(nx*ny), intent(out) :: vn     ! Cartesian velocity v at center of north face (m/s)
      real(8), dimension(nx*ny), intent(out) :: Uce    ! Contravariant velocity U at east face (m2/s)
      real(8), dimension(nx*ny), intent(out) :: Vcn    ! Contravariant velocity V at north face (m2/s)

      ! Axuliary variables

      integer :: i, j, np, npe, npn

      ! West boundary

      i = 1

      do j = 2, ny-1

         np  = nx * (j-1) + i

         npe = np + 1

         ue(np) = ( u(np) + u(npe) ) / 2.d0

         ve(np) = ( v(np) + v(npe) ) / 2.d0

         Uce(np) = ue(np) * ye(np) - ve(np) * xe(np)

      end do

      ! East boundary

      i = nx-1

      do j = 2, ny-1

         np  = nx * (j-1) + i

         npe = np + 1

         ue(np) = ( u(np) + u(npe) ) / 2.d0

         ve(np) = ( v(np) + v(npe) ) / 2.d0

         Uce(np) = ue(np) * ye(np) - ve(np) * xe(np)

      end do

      ! South boundary


      j = 1

      do i = 2, nx-1

         np  = nx * (j-1) + i

         npn = np + nx

         un(np) = ( u(np) + u(npn) ) / 2.d0

         vn(np) = ( v(np) + v(npn) ) / 2.d0

         Vcn(np) = vn(np) * xk(np) - un(np) * yk(np)

      end do

      ! North boundary

      j = ny-1

      do i = 2, nx-1

         np  = nx * (j-1) + i

         npn = np + nx

         un(np) = ( u(np) + u(npn) ) / 2.d0

         vn(np) = ( v(np) + v(npn) ) / 2.d0

         Vcn(np) = vn(np) * xk(np) - un(np) * yk(np)

      end do

   end subroutine get_velocities_at_boundary_faces


   !> \brief Calculates the contravariant velocities over the east boundary
   subroutine get_Ucbe_Vcbe(nx, ny, xe, ye, xke, yke, u, v, Ucbe, Vcbe) ! Output: last two
      implicit none
      integer, intent(in) :: nx !< Number of volumes in csi direction (real+fictitious)
      integer, intent(in) :: ny !< Number of volumes in eta direction (real+fictitious)
      real(8), dimension(nx*ny), intent(in) :: xe   !< x_eta at center of face east
      real(8), dimension(nx*ny), intent(in) :: ye   !< y_eta at center of face east
      real(8), dimension(nx*ny), intent(in) :: xke  !< x_csi at center of face east
      real(8), dimension(nx*ny), intent(in) :: yke  !< y_csi at center of face east
      real(8), dimension(nx*ny), intent(in) :: u    !< Cartesian velocity of the last iteraction
      real(8), dimension(nx*ny), intent(in) :: v    !< Cartesian velocity of the last iteraction
      real(8), dimension(ny),   intent(out) :: Ucbe !< Uc over the faces of the east  boundary (m2/s)
      real(8), dimension(ny),   intent(out) :: Vcbe !< Vc over the faces of the east  boundary (m2/s)


      ! Inner variables

      integer :: i, j, np, npe

      real(8) :: ue, ve

      ! East boundary

      i = nx-1

      do j = 2, ny-1

         np   = nx * (j-1) + i

         npe  = np + 1

         ue = 0.5d0 * ( u(np) + u(npe) )

         ve = 0.5d0 * ( v(np) + v(npe) )

         Ucbe(j) = ue * ye(np)  - ve * xe(np)

         Vcbe(j) = ve * xke(np) - ue * yke(np)

      end do

   end subroutine get_Ucbe_Vcbe


   ! Calculates the values of p at fictitious volumes
   subroutine get_p_extrapolation_to_fictitious(nx, ny, p) ! InOutput: last one
      implicit none
      integer, intent(in) :: nx  ! Number of volumes in csi direction (real+fictitious)
      integer, intent(in) :: ny  ! Number of volumes in eta direction (real+fictitious)
      real(8), dimension(nx*ny), intent(inout) :: p  ! Pressure at center of volumes

      ! Inner variables

      integer :: i, j, np, npsw, nps, npse, npw, npe, npnw, npn, npne


      ! North boundary

      j = ny

      do i = 2, nx-1

         np   = nx * (j-1) + i

         p(np) = 2.d0 * p(np-nx) - p(np-nx-nx)

      end do


      ! South boundary

      j = 1

      do i = 2, nx-1

         np   = nx * (j-1) + i

         p(np) = 2.d0 * p(np+nx) - p(np+nx+nx)

      end do


      ! East boundary

      i = nx

      do j = 2, ny-1

         np   = nx * (j-1) + i

         p(np) = 2.d0 * p(np-1) - p(np-2)

      end do


      ! West boundary

      i = 1

      do j = 2, ny-1

         np   = nx * (j-1) + i

         p(np) = 2.d0 * p(np+1) - p(np+2)

      end do

      ! SW

      i = 1
      j = 1

      np   = nx * (j-1) + i
      npn  = np + nx
      npe  = np + 1
      npne = npn + 1

      p(np) = ( p(npne) + p(npn) + p(npe) ) / 3.d0

      ! SE

      i = nx
      j = 1

      np   = nx * (j-1) + i
      npn  = np + nx
      npw  = np - 1
      npnw = npn - 1

      p(np) = ( p(npnw) + p(npn) + p(npw) ) / 3.d0

      ! NW

      i = 1
      j = ny

      np   = nx * (j-1) + i
      nps  = np - nx
      npe  = np + 1
      npse = nps + 1

      p(np) = ( p(npse) + p(nps) + p(npe) ) / 3.d0

      ! NE

      i = nx
      j = ny

      np   = nx * (j-1) + i
      nps  = np - nx
      npw  = np - 1
      npsw = nps - 1

      p(np) = ( p(npsw) + p(nps) + p(npw) ) / 3.d0

   end subroutine


   !> \brief Finds the index of the ogive-cylinder matching point over the
   !! south boundary (iocs). In this subroutine, it is assumed that the transition
   !! occurs when the slope of the south surface becomes zero
   subroutine get_ogive_cylinder_matching_point_south(nx, ny, xk, yk, iocs) ! Output: last one
      implicit none
      integer, intent(in) :: nx     !< Number of volumes in csi direction (real+fictitious)
      integer, intent(in) :: ny     !< Number of volumes in eta direction (real+fictitious)
      real(8), dimension (nx*ny), intent(in) :: xk !< x_csi at center of face north
      real(8), dimension (nx*ny), intent(in) :: yk !< y_csi at center of face north
      integer, intent(out) :: iocs !< Index of the ogive-cylinder matching point over the south boundary

      ! Inner variables

      integer :: i, j, np

      j = 1

      do i = 2, nx-1

         np   = nx * (j-1) + i

         iocs = i - 1

         if ( abs( atan2( yk(np) , xk(np) ) ) < 1.d-10 ) exit

      end do

      ! When no matching point is found, iocs = nx-1
      if ( i == nx-1 .and. abs( atan2( yk(np) , xk(np) ) ) > 1.d-10 ) iocs = i

   end subroutine get_ogive_cylinder_matching_point_south


   !> \brief Calculates the pressure foredrag over the south boundary
   subroutine get_cdfi(nx, ny, iiv, ifv, Rg, PF, TF, UF, rb, yk, rn, p, Cdfi)
      implicit none
      integer, intent(in) :: nx  !< Number of volumes in csi direction (real+fictitious)
      integer, intent(in) :: ny  !< Number of volumes in eta direction (real+fictitious)
      integer, intent(in) :: iiv !< Initial value of i
      integer, intent(in) :: ifv !< Final value of i
      real(8), intent(in) :: Rg  !< Perfect gas constant
      real(8), intent(in) :: PF  !< Far field pressure
      real(8), intent(in) :: TF  !< Far field temperature
      real(8), intent(in) :: UF  !< Far field speed
      real(8), intent(in) :: rb  !< Base radius of the rocket
      real(8), dimension (nx*ny), intent(in) :: yk !< y_csi at center of face north
      real(8), dimension (nx*ny), intent(in) :: rn !< Radius of the center of north face of volume P
      real(8), dimension (nx*ny), intent(in) :: p  !< Pressure at center o volume P
      real(8), intent(out) :: Cdfi !< Foredrag coefficient due pressure

      ! Auxiliary variables

      integer :: i, j, np, npn

      real(8) :: ROF, QF, pn


      ! Far field density
      ROF = PF / ( Rg * TF )

      ! Far field dynamic pressure
      QF = ROF * UF ** 2 / 2.d0

      Cdfi = 0.d0

      j = 1

      do i = iiv, ifv

         np   = nx * (j-1) + i

         npn  = np + nx

         pn = 0.5d0 * ( p(np) + p(npn) )

         Cdfi = Cdfi + ( pn - PF ) * rn(np) * yk(np)

      end do

      Cdfi = Cdfi * 2.d0 / ( rb ** 2 * QF )

   end subroutine get_cdfi


   !> \brief Calculates the viscous frodrag over the south boundary
   subroutine get_cdfv(nx, ny, iiv, ifv, Rg, PF, TF, UF, rb, xk, yk, rn, Jn, vln, u, v, Vcn, Cdfv)
      implicit none
      integer, intent(in) :: nx  ! Number of volumes in csi direction (real+fictitious)
      integer, intent(in) :: ny  ! Number of volumes in eta direction (real+fictitious)
      integer, intent(in) :: iiv !< Initial value of i
      integer, intent(in) :: ifv !< Final value of i
      real(8), intent(in) :: Rg  ! Perfect gas constant
      real(8), intent(in) :: PF  ! Far field pressure
      real(8), intent(in) :: TF  ! Far field temperature
      real(8), intent(in) :: UF  ! Far field speed
      real(8), intent(in) :: rb  ! Base radius of the rocket
      real(8), dimension (nx*ny), intent(in) :: xk   ! x_csi at center of face north
      real(8), dimension (nx*ny), intent(in) :: yk   ! y_csi at center of face north
      real(8), dimension (nx*ny), intent(in) :: rn   ! Radius of the center of north face of volume P
      real(8), dimension (nx*ny), intent(in) :: Jn   ! Jacobian at the center of north face of volume P
      real(8), dimension (nx*ny), intent(in) :: vln  ! Laminar viscosity at center of face north
      real(8), dimension (nx*ny), intent(in) :: u    ! Cartesian velocity of the last iteraction
      real(8), dimension (nx*ny), intent(in) :: v    ! Cartesian velocity of the last iteraction
      real(8), dimension (nx*ny), intent(in) :: Vcn  ! Contravariant velocity V at north face

      real(8), intent(out) :: Cdfv ! Frontal drag coefficient due viscosity

      ! Auxiliary variables

      integer :: i, j, np, npn

      real(8) :: ROF ! Far field density [kg/m3]
      real(8) :: QF  ! Far field dynamic pressure [Pa]
      real(8) :: Sxx ! Viscous stress tensor component [Pa/m2]
      real(8) :: Sxy ! Viscous stress tensor component [Pa/m2]

      ! Far field density
      ROF = PF / ( Rg * TF )

      ! Far field dynamic pressure
      QF = ROF * UF ** 2 / 2.d0

      j = 1

      Cdfv = 0.d0

      do i = iiv, ifv

         np   = nx * (j-1) + i

         npn  = np + nx

         ! Viscous stress tensor component [Pa/m2]
         Sxx = - 2.d0 * vln(np) * Jn(np) * ( yk(np) * ( u(npn) - u(np) ) + rn(npn) * Vcn(npn) / ( 3.d0 * rn(np) ) )

         ! Viscous stress tensor component [Pa/m2]
         Sxy = vln(np) * Jn(np) * ( xk(np) * ( u(npn) - u(np) ) - yk(np) * ( v(npn) - v(np) ) )

         Cdfv = Cdfv - Sxx * rn(np) * yk(np) + Sxy * rn(np) * xk(np)

      end do

      Cdfv = Cdfv * 2 / ( QF * rb**2 )

   end subroutine get_cdfv


end module
