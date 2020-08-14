module mod_intflow_procedures
  !
  !  SUBROUTINES
  !
  !  1) set_wall_temperature
  !  2) set_cp_and_gamma
  !  3) set_laminar_viscosity_at_nodes
  !  4) set_thermal_conductivity_at_nodes
  !  5) get_laminar_viscosity_at_faces
  !  6) get_thermal_conductivity_at_faces
  !  7) set_bcu
  !  8) set_bcv
  !  9) set_bcT
  ! 10) set_bcp
  ! 11) get_uin_vin_pin_Tin_Mw
  ! 12) get_plin_and_p_fictitious
  ! 13) get_u_v_at_fictitious_nodes_with_pl
  ! 14) get_boundary_simplec_coefficients
  ! 15) get_Uce_Vcn_at_boundary_faces
  ! 16) get_Uce_Vcn_at_boundary_faces_with_pl
  ! 16) get_isentropic_mass_flow
  ! 17) get_mach_area
  ! 18) get_initial_guess
  ! 19) get_boundary_nodes
  !
  ! Last update: 03 Oct 2011
  !

  use coefficients

  implicit none

contains

  subroutine set_wall_temperature( nx, Twall ) ! Last one is output
    implicit none
    integer, intent(in) :: nx  ! Number of volumes in csi direction (real+fictitious)
    real(8), dimension(nx), intent(out) :: Twall ! Wall temperature

    Twall = 300.d0

  end subroutine set_wall_temperature


  subroutine set_cp_and_gamma_internal_flow( nx, ny, Rg, cp, gcp) ! Last two are output
    implicit none
    integer, intent(in) :: nx  ! Number of volumes in csi direction (real+fictitious)
    integer, intent(in) :: ny  ! Number of volumes in eta direction (real+fictitious)
    real(8), intent(in) :: Rg  ! Perfect gas constant
    real(8), dimension(nx*ny), intent(out) :: cp  ! Specific heat at const pressure
    real(8), dimension(nx*ny), intent(out) :: gcp ! gcp = gamma = Cp/Cv at center of CV P

    cp = 1004.5d0

    gcp = cp / ( cp - Rg )

  end subroutine set_cp_and_gamma_internal_flow


  subroutine set_laminar_viscosity_at_nodes_internal_flow( nx, ny, vlp) ! Last one is output
    implicit none
    integer, intent(in) :: nx  ! Number of volumes in csi direction (real+fictitious)
    integer, intent(in) :: ny  ! Number of volumes in eta direction (real+fictitious)
    real(8), dimension (nx*ny), intent(out) :: vlp ! Laminar viscosity at center of volume P

    vlp = 1.d-5 ! (Units Pa.s)

  end subroutine set_laminar_viscosity_at_nodes_internal_flow


  subroutine set_thermal_conductivity_at_nodes_internal_flow( nx, ny, kp) ! Last one is output
    implicit none
    integer, intent(in) :: nx  ! Number of volumes in csi direction (real+fictitious)
    integer, intent(in) :: ny  ! Number of volumes in eta direction (real+fictitious)
    real(8), dimension (nx*ny), intent(out) :: kp ! Thermal conductivity at center of volume P

    kp = 1.0d-3 ! (Units W/m.K)

  end subroutine set_thermal_conductivity_at_nodes_internal_flow


  ! CAUTION: In order to apply correctly this function, laminar viscosity vlp must be known
  ! in real and fictitious volumes (except CV of corners SW, SE, NW and NE)
  subroutine get_laminar_viscosity_at_faces_internal_flow( nx, ny, vlp, vle, vln) ! Last two are output
    implicit none
    integer, intent(in) :: nx     ! Number of volumes in csi direction (real+fictitious)
    integer, intent(in) :: ny     ! Number of volumes in eta direction (real+fictitious)
    real(8), dimension (nx*ny), intent(in)  :: vlp ! Laminar viscosity at center of volume P
    real(8), dimension (nx*ny), intent(out) :: vle ! Laminar viscosity at center of face east
    real(8), dimension (nx*ny), intent(out) :: vln ! Laminar viscosity at center of face north

    ! Auxiliary variables
    integer :: i, j, np, npe, npn

    do j = 2, ny-1
       do i = 1, nx-1
          np   = nx * (j-1) + i
          npe  = np + 1

          vle(np) = 2.d0 * vlp(np) * vlp(npe) / ( vlp(np) + vlp(npe) )

       end do
    end do

    do i = 2, nx-1
       do j = 1, ny-1
          np  = nx * (j-1) + i
          npn = np + nx

          vln(np) = 2.d0 * vlp(np) * vlp(npn) / ( vlp(np) + vlp(npn) )

       end do
    end do

  end subroutine get_laminar_viscosity_at_faces_internal_flow


  ! CAUTION: In order to apply correctly this function, thermal conductivity kp must be known
  ! in real and fictitious volumes (except CV of corners SW, SE, NW and NE)
  subroutine get_thermal_conductivity_at_faces_internal_flow( nx, ny, kp, ke, kn) ! Output: last two entries
    implicit none
    integer, intent(in) :: nx     ! Number of volumes in csi direction (real+fictitious)
    integer, intent(in) :: ny     ! Number of volumes in eta direction (real+fictitious)
    real(8), dimension (nx*ny), intent(in)  :: kp ! Thermal conductivity at center of volume P
    real(8), dimension (nx*ny), intent(out) :: ke ! Thermal conductivity at center of face east
    real(8), dimension (nx*ny), intent(out) :: kn ! Thermal conductivity at center of face north

    ! Auxiliary variables
    integer :: i, j, np, npe, npn

    do j = 2, ny-1
       do i = 1, nx-1
          np   = nx * (j-1) + i
          npe  = np + 1

          ke(np) = 2.d0 * kp(np) * kp(npe) / ( kp(np) + kp(npe) )

       end do
    end do

    do i = 2, nx-1
       do j = 1, ny-1
          np  = nx * (j-1) + i
          npn = np + nx

          kn(np) = 2.d0 * kp(np) * kp(npn) / ( kp(np) + kp(npn) )

       end do
    end do

  end subroutine get_thermal_conductivity_at_faces_internal_flow

  subroutine set_bcu( nx, ny, modvis, x, xp, u, au, b) ! Output: last two entries
    implicit none
    integer, intent(in) :: nx  ! Number of volumes in csi direction (real+fictitious)
    integer, intent(in) :: ny  ! Number of volumes in eta direction (real+fictitious)
    integer, intent(in) :: modvis ! modvis = 0 -> Euler;  modvis = 1 -> Navier-Stokes
    real(8), dimension (nx*ny),   intent(in)  :: x  ! Coord. x of the northest corner of volume P
    real(8), dimension (nx*ny),   intent(in)  :: xp ! X coord. of the centroid of volume P
    real(8), dimension (nx*ny),   intent(in)  :: u  ! Cartesian velocity of the last iteraction
    real(8), dimension (nx*ny,9), intent(out) :: au ! Coefficients of the linear system for u
    real(8), dimension (nx*ny),   intent(out) :: b  ! Source vector of the linear system
    ! Auxiliary variables
    integer :: i, j, np, npe, npee, npw, npww

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
       do i = 1, nx
          np   = nx * (j-1) + i

          au(np,:) =  0.d0
          au(np,2) = -1.d0
          au(np,5) =  1.d0

          b(np) = 0.d0
       end do
    end if
  end subroutine set_bcu

  subroutine set_bcv( nx, ny, modvis, x, xp, v, av, b) ! Output: last two entries
    implicit none
    integer, intent(in) :: nx     ! Number of volumes in csi direction (real+fictitious)
    integer, intent(in) :: ny     ! Number of volumes in eta direction (real+fictitious)
    integer, intent(in) :: modvis ! modvis = 0 -> Euler;  modvis = 1 -> Navier-Stokes
    real(8), dimension (nx*ny),   intent(in)  :: x  ! Coord. x of the northest corner of volume P
    real(8), dimension (nx*ny),   intent(in)  :: xp ! X coord. of the centroid of volume P
    real(8), dimension (nx*ny),   intent(in)  :: v  ! Cartesian velocity of the last iteraction
    real(8), dimension (nx*ny,9), intent(out) :: av ! Coefficients of the linear system for v
    real(8), dimension (nx*ny),   intent(out) :: b  ! Source vector of the linear system
    ! Auxiliary variables
    integer :: i, j, np, npe, npee, npw, npww

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
       do i = 1, nx
          np   = nx * (j-1) + i

          av(np,:) =  0.d0
          av(np,2) = -1.d0
          av(np,5) =  1.d0

          b(np) = 0.d0
       end do
    end if
  end subroutine set_bcv

  ! CAUTION: Twall must be difined in the fictitious volumes too
  subroutine set_bcT( nx, ny, ccTw, x, xp, Tin, Twall, T, at, b) ! Output: last two entries
    implicit none
    integer, intent(in) :: nx     ! Number of volumes in csi direction (real+fictitious)
    integer, intent(in) :: ny     ! Number of volumes in eta direction (real+fictitious)
    integer, intent(in) :: ccTw   ! ccTw = 0 -> adiabatic;  ccTw = 1 -> prescribed temperature
    real(8), dimension(nx*ny),   intent(in) ::  x     ! Coord. x of the northest corner of volume P
    real(8), dimension(nx*ny),   intent(in) ::  xp    ! X coord. of the centroid of volume P
    real(8), dimension(ny),      intent(in) ::  Tin   ! Temperature in the entrance
    real(8), dimension(nx),      intent(in) ::  Twall ! Wall temperature
    real(8), dimension(nx*ny),   intent(in) ::  T     ! Temperature of the last iteraction
    real(8), dimension(nx*ny,9), intent(out) :: aT    ! Coefficients of the linear system
    real(8), dimension(nx*ny),   intent(out) :: b     ! Source vector of the linear system

    ! Auxiliary variables
    integer :: i, j, np, npw, npww

    ! West boundary
    i = 1
    do j = 2, ny-1
       np   = nx * (j-1) + i

       at(np,:) = 0.d0
       at(np,6) = 1.d0
       at(np,5) = 1.d0

       b(np) = 2.d0 * Tin(j)
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

       b(np) = 2.d0 * (x(npw) - xp(npw))/(xp(npw) - xp(npww)) * (T(npw) - T(npww))
    end do

    ! South boundary ( SW and SE corners are included)
    j = 1
    do i = 1, nx
       np = nx * (j-1) + i

       at(np,:) =  0.d0
       at(np,8) = -1.d0
       at(np,5) =  1.d0

       b(np) = 0.d0
    end do

    ! North boundary ( NW and NE corners are included)
    j = ny
    if ( ccTw == 0 ) then ! Adiabatic
       do i = 1, nx
          np = nx * (j-1) + i

          at(np,:) =  0.d0
          at(np,2) = -1.d0
          at(np,5) =  1.d0

          b(np) = 0.d0
       end do
    else ! Prescribed temperature
       do i = 1, nx
          np = nx * (j-1) + i

          at(np,:) = 0.d0
          at(np,2) = 1.d0
          at(np,5) = 1.d0

          b(np) = 2.d0 * Twall(i)
       end do
    end if

  end subroutine set_bcT

  subroutine set_bcp( nx, ny, x, xp, plin, pl, ap, b) ! Output: last two entries
    implicit none
    integer, intent(in) :: nx     ! Number of volumes in csi direction (real+fictitious)
    integer, intent(in) :: ny     ! Number of volumes in eta direction (real+fictitious)
    real(8), dimension(nx*ny),   intent(in)  :: x     ! Coord. x of the northest corner of volume P
    real(8), dimension(nx*ny),   intent(in)  :: xp    ! X coord. of the centroid of volume P
    real(8), dimension(ny),      intent(in)  :: plin  ! Pressure correction in the entrance
    real(8), dimension(nx*ny),   intent(in)  :: pl    ! Pressure correction
    real(8), dimension(nx*ny,5), intent(out) :: ap    ! Coefficients of the linear system
    real(8), dimension(nx*ny),   intent(out) :: b     ! Source vector of the linear system

    ! Auxliary variables
    integer :: i, j, np, npw, npww
    ! West boundary
    i = 1
    do j = 2, ny-1
       np = nx * (j-1) + i

       ap(np,:) = 0.d0
       ap(np,4) = 1.d0
       ap(np,3) = 1.d0

       b(np) = 2.d0 * plin(j)
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

       b(np) = 2.d0 * ( x(npw) - xp(npw) ) / ( xp(npw) - xp(npww) ) * ( pl(npw) - pl(npww) )
    end do

    ! South boundary
    j = 1
    do i = 2, nx-1
       np  = nx * (j-1) + i

       ap(np,:) =  0.d0
       ap(np,5) = -1.d0
       ap(np,3) =  1.d0

       b(np) = 0.d0
    end do

    ! North boundary
    j = ny
    do i = 2, nx-1
       np  = nx * (j-1) + i

       ap(np,:) =  0.d0
       ap(np,1) = -1.d0
       ap(np,3) =  1.d0

       b(np) = 0.d0
    end do
  end subroutine set_bcp

  subroutine get_uin_vin_pin_Tin_Mw( nx, ny, gamma, Rg, po, T0, u & ! Input
       ,                             uin, vin, pin, Tin, Mw)        ! Output
    implicit none
    integer, intent(in) :: nx    ! Number of volumes in csi direction (real+fictitious)
    integer, intent(in) :: ny    ! Number of volumes in eta direction (real+fictitious)
    real(8), intent(in) :: gamma ! gamma = Cpo / Cvo in the chamber
    real(8), intent(in) :: Rg    ! Perfect gas constant
    real(8), intent(in) :: po    ! Stagnation pressure in the chamber
    real(8), intent(in) :: T0    ! Stagnation temperature in the chamber
    real(8), dimension (nx*ny), intent(in)  :: u   ! Cartesian velocity of the last iteraction
    real(8), dimension (ny),    intent(out) :: uin ! Velocity u in the entrance
    real(8), dimension (ny),    intent(out) :: vin ! Velocity v in the entrance
    real(8), dimension (ny),    intent(out) :: pin ! Pressure in the entrance
    real(8), dimension (ny),    intent(out) :: Tin ! Temperature in the entrance
    real(8), dimension (ny),    intent(out) :: Mw  ! Mach number in the entrance

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

  subroutine get_plin_and_p_fictitious( nx, ny, x, xp, pina, pin, plin, p) ! Output: last two entries
    implicit none
    integer, intent(in) :: nx    ! Number of volumes in csi direction (real+fictitious)
    integer, intent(in) :: ny    ! Number of volumes in eta direction (real+fictitious)
    real(8), dimension (nx*ny),   intent(in)  :: x  ! Coord. x of the northest corner of volume P
    real(8), dimension (nx*ny),   intent(in)  :: xp ! X coord. of the centroid of volume P
    real(8), dimension(ny),    intent(in)  :: pina ! Pressure in the entrance at a time step before
    real(8), dimension(ny),    intent(in)  :: pin  ! Pressure in the entrance
    real(8), dimension(ny),    intent(out) :: plin ! Pressure correction in the entrance
    real(8), dimension(nx*ny), intent(inout) :: p    ! Pressure at center of volumes
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
    nps  = np - nx
    npn  = np + nx
    npw  = np - 1
    npe  = np + 1
    npsw = nps - 1
    npse = nps + 1
    npnw = npn - 1
    npne = npn + 1
    p(np) = (p(npn)+p(npe)+p(npne)) / 3.d0

    ! SE corner
    i = nx
    j = 1
    np   = nx * (j-1) + i
    nps  = np - nx
    npn  = np + nx
    npw  = np - 1
    npe  = np + 1
    npsw = nps - 1
    npse = nps + 1
    npnw = npn - 1
    npne = npn + 1
    p(np) = (p(npn)+p(npw)+p(npnw)) / 3.d0

    ! NW corner
    i = 1
    j = ny
    np   = nx * (j-1) + i
    nps  = np - nx
    npn  = np + nx
    npw  = np - 1
    npe  = np + 1
    npsw = nps - 1
    npse = nps + 1
    npnw = npn - 1
    npne = npn + 1
    p(np) = (p(nps)+p(npe)+p(npse)) / 3.d0

    ! NE corner
    i = nx
    j = ny
    np   = nx * (j-1) + i
    nps  = np - nx
    npn  = np + nx
    npw  = np - 1
    npe  = np + 1
    npsw = nps - 1
    npse = nps + 1
    npnw = npn - 1
    npne = npn + 1
    p(np) = (p(nps)+p(npw)+p(npsw)) / 3.d0

  end subroutine get_plin_and_p_fictitious


  subroutine get_u_v_at_fictitious_nodes_with_pl( nx, ny, modvis, x, xp, u, v) ! InOutput: last two entries
    implicit none
    integer, intent(in) :: nx  ! Number of volumes in csi direction (real+fictitious)
    integer, intent(in) :: ny  ! Number of volumes in eta direction (real+fictitious)
    integer, intent(in) :: modvis ! modvis = 0 -> Euler;  modvis = 1 -> Navier-Stokes
    real(8), dimension (nx*ny), intent(in)    :: x  ! Coord. x of the northest corner of volume P
    real(8), dimension (nx*ny), intent(in)    :: xp ! Coord. x of the centroid of volume P
    real(8), dimension (nx*ny), intent(inout) :: u  ! Cartesian velocity
    real(8), dimension (nx*ny), intent(inout) :: v  ! Cartesian velocity

    ! Auxiliary variables
    integer :: i, j, np, npe, npee, npw, npww, nps, npn

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

       do i = 1, nx
          np   = nx * (j-1) + i
          nps  = np - nx

          u(np) = u(nps)
          v(np) = v(nps)

       end do

    end if
  end subroutine get_u_v_at_fictitious_nodes_with_pl

  subroutine get_boundary_simplec_coefficients_internal_flow( nx, ny, de, dn) ! InOutput: de, dn
    implicit none
    integer, intent(in) :: nx  ! Number of volumes in csi direction (real+fictitious)
    integer, intent(in) :: ny  ! Number of volumes in eta direction (real+fictitious)
    real(8), dimension (nx*ny), intent(inout) :: de  ! SIMPLEC coefficient for Uce
    real(8), dimension (nx*ny), intent(inout) :: dn  ! SIMPLEC coefficient for Vcn

    ! Auxiliary variables
    integer :: i, j, np, npe, npw

    ! West boundary

    i = 1

    do j = 2, ny-1

       np   = nx * (j-1) + i

       npe  = np + 1

       de(np) = de(npe)

    end do

    ! East boundary

    i = nx-1

    do j = 2, ny-1

       np   = nx * (j-1) + i

       npw  = np - 1

       de(np) = de(npw)

    end do

    ! South boundary

    j = 1

    do i = 2, nx-1

       np   = nx * (j-1) + i

       dn(np) = 0.d0

    end do

    ! North boundary

    j = ny-1

    do i = 2, nx-1

       np   = nx * (j-1) + i

       dn(np) = 0.d0

    end do

  end subroutine get_boundary_simplec_coefficients_internal_flow

     ! Calculates the SIMPLEC coefficients at the boundary faces
   subroutine get_boundary_simplec_coefficients_internal_flow2( nx, ny, modvis, due, dve, dun &
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

      if ( modvis == 0 ) then ! Euler

         j = ny-1

         do i = 2, nx-1

            np   = nx * (j-1) + i

            nps  = np - nx

            dun(np) = dun(nps)

            dvn(np) = dvn(nps)

            dn(np) = 0.d0

         end do

      else ! NS

         j = ny-1

         do i = 2, nx-1

            np   = nx * (j-1) + i

            dun(np) = 0.d0

            dvn(np) = 0.d0

            dn(np) = 0.d0

         end do

      end if

   end subroutine get_boundary_simplec_coefficients_internal_flow2

  subroutine get_ue_un_ve_vn_at_boundary_faces( nx, ny, u, v, ue, un, ve, vn) ! InOutput: ue, un, ve, vn
    implicit none
    integer, intent(in) :: nx  ! Number of volumes in csi direction (real+fictitious)
    integer, intent(in) :: ny  ! Number of volumes in eta direction (real+fictitious)
    real(8), dimension (nx*ny), intent(in)    :: u    !< Cartesian velocity of the last iteraction
    real(8), dimension (nx*ny), intent(in)    :: v    !< Cartesian velocity of the last iteraction
    real(8), dimension (nx*ny), intent(inout) :: ue   ! Cartesian velocity u at center of east face
    real(8), dimension (nx*ny), intent(inout) :: ve   ! Cartesian velocity v at center of east face
    real(8), dimension (nx*ny), intent(inout) :: un   ! Cartesian velocity u at center of north face
    real(8), dimension (nx*ny), intent(inout) :: vn   ! Cartesian velocity v at center of north face


    ! Auxiliary variables
    integer :: i, j, np, npe, npn

    ! West boundary

    i = 1

    do j = 2, ny-1

       np   = nx * (j-1) + i

       npe  = np + 1

       ve(np) = 0.d0

       ue(np) = (u(np)+u(npe))/2.d0

    end do

    ! East boundary

    i = nx-1

    do j = 2, ny-1

       np   = nx * (j-1) + i

       npe  = np + 1

       ve(np) = (v(np)+v(npe))/2.d0

       ue(np) = (u(np)+u(npe))/2.d0

    end do

    ! South boundary

    j = 1

    do i = 2, nx-1

       np   = nx * (j-1) + i

       npn  = np + nx

       vn(np) = 0.d0

       un(np) = (u(np)+u(npn))/2.d0

    end do


    ! North boundary

    j = ny-1

    do i = 2, nx-1

       np   = nx * (j-1) + i

       npn  = np + nx

       vn(np) = v(np)

       un(np) = u(np)

    end do

  end subroutine get_ue_un_ve_vn_at_boundary_faces




  subroutine get_Uce_Vcn_at_boundary_faces( nx, ny, ye, u, Uce, Vcn) ! Output: Uce, Vcn
    implicit none
    integer, intent(in) :: nx  ! Number of volumes in csi direction (real+fictitious)
    integer, intent(in) :: ny  ! Number of volumes in eta direction (real+fictitious)
    real(8), dimension (nx*ny), intent(in)    :: ye  ! y_eta at face east of volume P
    real(8), dimension (nx*ny), intent(in)    :: u   ! Cartesian velocity of the last iteraction
    real(8), dimension (nx*ny), intent(out)   :: Uce ! Contravariant velocity U at east face
    real(8), dimension (nx*ny), intent(out)   :: Vcn ! Contravariant velocity V at north face

    ! Auxiliary variables
    integer :: i, j, np, npe, npw

    ! West boundary

    i = 1

    do j = 2, ny-1

       np   = nx * (j-1) + i

       npe  = np + 1

       Uce(np) = 0.5d0 * ( u(npe) + u(np) ) * ye(np)

    end do

    ! East boundary

    i = nx-1

    do j = 2, ny-1

       np   = nx * (j-1) + i

       npw  = np - 1

       npe  = np + 1

       Uce(np) = 0.5d0 * ( u(npe) + u(np) ) * ye(np)

    end do

    ! South boundary

    j = 1

    do i = 2, nx-1

       np   = nx * (j-1) + i

       Vcn(np) = 0.d0

    end do

    ! North boundary

    j = ny-1

    do i = 2, nx-1

       np   = nx * (j-1) + i

       Vcn(np) = 0.d0

    end do

  end subroutine get_Uce_Vcn_at_boundary_faces

    subroutine get_velocities_at_boundary_faces_internal_flow( nx, ny, xe, ye, xk, yk, u, v  & ! Input
       ,                                       ue, ve, un, vn, Uce, Vcn )      ! Input and Output
    implicit none
    integer, intent(in) :: nx     ! Number of volumes in csi direction (real+fictitious)
    integer, intent(in) :: ny     ! Number of volumes in eta direction (real+fictitious)
    real(8), dimension (nx*ny), intent(in)    :: xe  ! x_eta at center of face east
    real(8), dimension (nx*ny), intent(in)    :: ye  ! y_eta at center of face east
    real(8), dimension (nx*ny), intent(in)    :: xk  ! x_csi at center of face north
    real(8), dimension (nx*ny), intent(in)    :: yk  ! y_csi at center of face north
    real(8), dimension (nx*ny), intent(in)    :: u   ! Cartesian velocity of the present iteraction
    real(8), dimension (nx*ny), intent(in)    :: v   ! Cartesian velocity of the present iteraction
    real(8), dimension (nx*ny), intent(inout) :: ue  ! Cartesian velocity u at center of east face
    real(8), dimension (nx*ny), intent(inout) :: ve  ! Cartesian velocity v at center of east face
    real(8), dimension (nx*ny), intent(inout) :: un  ! Cartesian velocity u at center of north face
    real(8), dimension (nx*ny), intent(inout) :: vn  ! Cartesian velocity v at center of north face
    real(8), dimension (nx*ny), intent(inout) :: Uce ! Contravariant velocity U at east face
    real(8), dimension (nx*ny), intent(inout) :: Vcn ! Contravariant velocity V at north face

    ! Auxiliary variables
    integer :: i, j, np, npn, npe

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

    j = ny-1
    do i = 2, nx-1

       np   = nx * (j-1) + i
       npn  = np + nx

       un(np) = ( u(np) + u(npn) ) / 2.d0

       vn(np) = ( v(np) + v(npn) ) / 2.d0

       Vcn(np) = 0.d0

    end do

  end subroutine get_velocities_at_boundary_faces_internal_flow


  subroutine get_Uce_Vcn_at_boundary_faces_with_pl( nx, ny, pl, de, Uce) ! InOutput: Last entry
    implicit none
    integer, intent(in) :: nx  ! Number of volumes in csi direction (real+fictitious)
    integer, intent(in) :: ny  ! Number of volumes in eta direction (real+fictitious)
    real(8), dimension(nx*ny), intent(in)    :: pl  ! Pressure correction
    real(8), dimension(nx*ny), intent(in)    :: de  ! SIMPLEC coefficient for Uce
    real(8), dimension(nx*ny), intent(inout) :: Uce ! Contravariant velocity U at east face

    ! Auxiliary variables
    integer :: i, j, np, npe

    ! West boundary

    i = 1

    do j = 2, ny-1

       np   = nx * (j-1) + i

       npe  = np + 1

       Uce(np) = Uce(np) + de(np) * ( pl(np) - pl(npe) )

    end do

    ! West boundary

    i = nx-1

    do j = 2, ny-1

       np   = nx * (j-1) + i

       npe  = np + 1

       Uce(np) = Uce(np) + de(np) * ( pl(np) - pl(npe) )

    end do

    ! North and South Vcn do not neet to be changed

  end subroutine get_Uce_Vcn_at_boundary_faces_with_pl


  subroutine get_isentropic_mass_flow( po, T0, gamma, Rg, Sg , fm1D)
    implicit none
    real(8), intent(in)  :: po    ! Stagnation pressure
    real(8), intent(in)  :: T0    ! Stagnation temperature
    real(8), intent(in)  :: gamma ! Specific heat ratio
    real(8), intent(in)  :: Rg    ! Perfect gas constant
    real(8), intent(in)  :: Sg    ! Throttle area
    real(8), intent(out) :: fm1D  ! Mass flow rate

    fm1D = po * Sg * sqrt( gamma/ (Rg*T0) * (2.d0/(gamma+1.d0)) ** ((gamma+1.d0)/(gamma-1.d0)) )

  end subroutine get_isentropic_mass_flow


  subroutine get_mach_area( kf, ar, gamma, M )
    implicit none
    integer, intent(in)  :: kf    ! Kind of flow: 0 = subsonic, 1 = supersonic
    real(8), intent(in)  :: ar    ! Area ratio = local area / throttle area
    real(8), intent(in)  :: gamma ! Specific heat ratio
    real(8), intent(out) :: M     ! Mach number at the local area

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


  subroutine get_initial_guess( nx, ny, ig, modvis, beta, po, T0, gamma  & ! Input
       ,                        Rg, Sg, ye, yk, radius, rn, x, xp        & ! Input
       ,                        M1D, p1D, T1D, u1D, p, T, u, v, ue       & ! Output
       ,                        un, Uce, Vcn, uin, vin, pin, Tin, Mw     & ! Output
       ,                        fm1D, Fd1D, Fpv1D, de, dn, ro, roe, ron  & ! Output
       ,                        a, ap, b ) ! Output
    implicit none
    integer, intent(in) :: nx     ! Number of volumes in csi direction (real+fictitious)
    integer, intent(in) :: ny     ! Number of volumes in eta direction (real+fictitious)
    integer, intent(in) :: ig
    integer, intent(in) :: modvis ! modvis = 0 -> Euler;  modvis = 1 -> Navier-Stokes
    real(8), intent(in) :: beta   ! Constant of the UDS/CDS mixing scheme
    real(8), intent(in) :: po     ! Stagnation pressure in the chamber
    real(8), intent(in) :: T0     ! Stagnation temperature in the chamber
    real(8), intent(in) :: gamma  ! gamma = Cpo / Cvo in the chamber
    real(8), intent(in) :: Rg     ! Perfect gas constant
    real(8), intent(in) :: Sg     ! Throttle area
    real(8), dimension(nx*ny), intent(in) :: ye     ! y_eta at face east of volume P
    real(8), dimension(nx*ny), intent(in) :: yk     ! y_csi at face north of volume P
    real(8), dimension(nx*ny), intent(in) :: radius ! Radius of northest corner of volume P
    real(8), dimension(nx*ny), intent(in) :: rn     ! Radius of the center of north face of volume P
    real(8), dimension(nx*ny), intent(in) :: x      ! Coord. x of the northest corner of volume P
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
    real(8), dimension(nx*ny), intent(out) :: un  ! Cartesian velocity u at center of north face
    real(8), dimension(nx*ny), intent(out) :: Uce ! Contravariant velocity U at east face
    real(8), dimension(nx*ny), intent(out) :: Vcn ! Contravariant velocity V at north face

    real(8), dimension (ny),   intent(out) :: uin ! Velocity u in the entrance
    real(8), dimension (ny),   intent(out) :: vin ! Velocity v in the entrance
    real(8), dimension (ny),   intent(out) :: pin ! Pressure in the entrance
    real(8), dimension (ny),   intent(out) :: Tin ! Temperature in the entrance
    real(8), dimension (ny),   intent(out) :: Mw  ! Mach number in the entrance

    real(8), dimension(nx*ny), intent(out) :: de  ! Simplec coef. for the contravariant velocity U (east face)
    real(8), dimension(nx*ny), intent(out) :: dn  ! Simplec coef. for the contravariant velocity V (north face)

    real(8), dimension(nx*ny), intent(out) :: ro  ! Specific mass (absolute density) at center of volumes
    real(8), dimension(nx*ny), intent(out) :: roe ! Absolute density at east face
    real(8), dimension(nx*ny), intent(out) :: ron ! Absolute density at north face

    real(8), dimension(nx*ny,9), intent(out) :: a  ! Coefficients of the linear system for u, v, T
    real(8), dimension(nx*ny,5), intent(out) :: ap ! Coefficients of the linear system
    real(8), dimension(nx*ny),   intent(out) :: b  ! Source vector of the linear system


    ! Auxiliary variables
    real(8), parameter :: pi = acos(-1.d0)

    integer :: i, j, np, npe, npw, npww, kf
    real(8) :: Si, Me1D, pe1D, Te1D, ue1D, aux

    call get_isentropic_mass_flow( po, T0, gamma, Rg, Sg , fm1D)


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
    call get_u_v_at_fictitious_nodes_with_pl( nx, ny, modvis, x, xp, u, v) ! InOutput: last two entries

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

    ! Contravariant velocities at boundary faces
    call get_Uce_Vcn_at_boundary_faces( nx, ny, ye, u, Uce, Vcn) ! Output: Uce, Vcn

    ! Simplec coefficients at boundary faces
    call get_boundary_simplec_coefficients_internal_flow( nx, ny, de, dn) ! InOutput: de, dn

    ! Specific mass at nodes
    call get_density_at_nodes( nx, ny, Rg, p, T, ro) ! ro is output

    ! Specific mass at faces
    call get_density_at_faces( nx, ny, beta, ro, Uce, Vcn, roe, ron) ! roe and ron are output

    a  = 0.d0
    ap = 0.d0
    b  = 0.d0

    ap(:,3) = 1.d0

  end subroutine get_initial_guess

  subroutine get_boundary_nodes( nx, ny, ig, Sg, rcg) ! Output: last three entries
    implicit none
    integer, intent(in)  :: nx
    integer, intent(in)  :: ny
    integer, intent(out) :: ig  ! i=ig when the throttle cross section is in the east face of the CV
    real(8), intent(out) :: Sg  ! area of the throttle
    real(8), intent(out) :: rcg ! Throttle curvature radius

    integer :: NXc, NXt
    real(8) :: dx
    real(8) :: Lc ! comprimento da câmara (m)
    real(8) :: Ln ! comprimento da tubeira (m)
    real(8) :: L  ! comprimento total: câmara + tubeira (m)
    real(8) :: Ri ! raio na entrada da tubeira (m)
    real(8) :: Rg ! raio na garganta da tubeira (m)
    real(8) :: d  ! amplitude da cossenóide (m)

    real(8) :: x(nx), y
    real(8),parameter :: pi = acos(-1.d0)

    integer :: i, j

    open(10,file='./mach2d_input/gridboundary.dat')

    NXc = (nx-2) / 5
    NXt = nx - 2 - Nxc

    ! *** dados ***

    Lc = 0.100d0
    Ln = 0.400d0
    Ri = 0.100d0
    Rg = 0.090d0

    Sg = pi * Rg**2

    L = Lc + Ln
    d = ( ri - rg ) / 2

    ! *** CONTORNO SUL ***
    write(10,*) "south boundary"
    ! câmara
    J = 1
    I = 1
    DX = Lc / Nxc
    x(i) = 0.d0
    y = 0.d0
    write(10,*) x(i), y
    DO I = 2, NXc
       x(i) = x(i-1) + dx
       write(10,*) x(i), y
    end do
    I = NXc + 1
    X(i) = Lc
    Y = 0.0d0
    write(10,*) x(i), y

    ! tubeira cossenoidal
    DX = Ln / NXt
    DO I = NXc+2, NX-2
       x(i) = x(i-1) + DX
       y = 0.0d0
       write(10,*) x(i), y
    end do
    I = NX-1
    X(i) = L
    Y = 0.0d0
    write(10,*) x(i), y
    write(10,*) "north boundary"
    ! *** CONTORNO NORTE ***

    ! câmara
    J = ny-1
    DO I = 1, NXc
       Y = Ri
       write(10,*) x(i), y
    end do
    I = NXc + 1
    Y = Ri
    write(10,*) x(i), y
    ! índice i da garganta
    ig = Nxc * 3 + 1

    ! tubeira cossenoidal
    DO I = NXc+2, NX-2
       Y = Rg + d * ( 1 + dcos ( 2*pi*(x(i)-Lc) / Ln ) )
       write(10,*) x(i), y
    end do
    I = NX-1
    Y = Ri
    write(10,*) x(i), y

    close(10)

    ! Throttle curvature radius

    rcg = ( Ln ** 2 ) / ( 2 * (pi**2) * ( Ri - Rg ) )

  end subroutine get_boundary_nodes


   !> \brief Calculates the temperature based on the conservation of the total enthalpy
   !! Valid for Euler model with constant thermo-physical coefficients.
   !! Temperature is extrapolated to fictitious volumes with CDS the scheme
   subroutine get_T_from_H_conservation_internal_flow(nx, ny, CPF, HF, u, ue, un, v, ve, vn, T, Tbe, Tbw, Tbn, Tbs) ! Output: last 5
      implicit none
      integer, intent(in) :: nx  !< Number of volumes in csi direction (real+fictitious)
      integer, intent(in) :: ny  !< Number of volumes in eta direction (real+fictitious)
      real(8), intent(in) :: CPF !< Free stream Cp (J/kg.K)
      real(8), intent(in) :: HF  !< Free stream total enthalpy (m2/s2)
      real(8), dimension (nx*ny), intent(in)  :: u    !< Cartesian velocity at P
      real(8), dimension (nx*ny), intent(in)  :: ue   !< Cartesian velocity u at center of east face
      real(8), dimension (nx*ny), intent(in)  :: un   !< Cartesian velocity u at center of north face
      real(8), dimension (nx*ny), intent(in)  :: v    !< Cartesian velocity at P
      real(8), dimension (nx*ny), intent(in)  :: ve   !< Cartesian velocity v at center of east face
      real(8), dimension (nx*ny), intent(in)  :: vn   !< Cartesian velocity v at center of north face
      real(8), dimension (nx*ny), intent(out) :: T    !< Temperature at P
      real(8), dimension (ny),    intent(out) :: Tbe  !< Temperature over the  east boundary (K)
      real(8), dimension (ny),    intent(in) :: Tbw  !< Temperature over the  west boundary (K)
      real(8), dimension (nx),    intent(out) :: Tbn  !< Temperature over the north boundary (K)
      real(8), dimension (nx),    intent(out) :: Tbs  !< Temperature over the south boundary (K)

      ! Inner variables

      integer :: i, j, np, npe, npn



      ! Temperature on the nodes

      T = ( HF - (u**2+v**2) / 2.d0 ) / CPF


      ! Temperature on the east boundary

      i = nx-1

      do j = 2, ny-1

         np   = nx * (j-1) + i

         Tbe(j) = ( HF - (ue(np)**2+ve(np)**2) / 2.d0 ) / CPF

         ! Extrapolating to fictitious

         npe  = np + 1

         T(npe) = 2.d0 * Tbe(j) - T(np)

      end do


      ! Temperature on the west boundary

      i = 1

      do j = 2, ny-1

         np   = nx * (j-1) + i

         !Tbw(j) = ( HF - (ue(np)**2+ve(np)**2) / 2.d0 ) / CPF

         ! Extrapolating to fictitious

         npe  = np + 1

         T(np) = 2.d0 * Tbw(j) - T(npe)

      end do


      ! Temperature on the north boundary

      j = ny-1

      do i = 2, nx-1

         np   = nx * (j-1) + i

         Tbn(i) = ( HF - (un(np)**2+vn(np)**2) / 2.d0 ) / CPF

         ! Extrapolating to fictitious

         npn  = np + nx

         T(npn) = 2.d0 * Tbn(i) - T(np)

      end do


      ! Temperature on the north boundary

      j = 1

      do i = 2, nx-1

         np   = nx * (j-1) + i

         Tbs(i) = ( HF - (un(np)**2+vn(np)**2) / 2.d0 ) / CPF

         ! Extrapolating to fictitious

         npn  = np + nx

         T(np) = 2.d0 * Tbs(i) - T(npn)

      end do


   end subroutine get_T_from_H_conservation_internal_flow

  subroutine get_mass_flow_rate_and_thrust( nx, ny, re, roe, u, Uce, fmi, fme, Fd)
    implicit none
    integer, intent(in) :: nx  ! Number of volumes in csi direction (real+fictitious)
    integer, intent(in) :: ny  ! Number of volumes in eta direction (real+fictitious)
    real(8), dimension (nx*ny), intent(in) :: re   ! Radius of the center of east face of volume P
    real(8), dimension (nx*ny), intent(in) :: roe  ! Absolute density at east face
    real(8), dimension (nx*ny), intent(in) :: u    ! Cartesian velocity of the last iteraction
    real(8), dimension (nx*ny), intent(in) :: Uce  ! Contravariant velocity U at east face
    real(8), intent(out) :: fmi ! Mass flow rate at entrance
    real(8), intent(out) :: fme ! Mass flow rate at exit
    real(8), intent(out) :: Fd  ! Thrust

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

end module
