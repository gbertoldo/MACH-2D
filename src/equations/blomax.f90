module blomax
    implicit none

    ! Constants of the Baldwin-Lomax model
    real(8), private :: kappa = 0.40d0
    real(8), private ::   Aop = 26.d0
    real(8), private :: alpha = 0.0168d0
    real(8), private ::   Ccp = 1.6d0
    real(8), private :: Ckleb = 0.3d0
    real(8), private ::   Cwk = 1.0d0

    ! Metrics componets used in the calculation of the vorticity
    real(8), allocatable, dimension(:), private :: xep ! x_eta_P
    real(8), allocatable, dimension(:), private :: xkp ! x_csi_P
    real(8), allocatable, dimension(:), private :: yep ! y_eta_P
    real(8), allocatable, dimension(:), private :: ykp ! y_csi_P

    ! Vorticity ( omg = omega )

    real(8), allocatable, dimension(:), private :: omg ! Vorticity

contains

    ! Sets the coefficients of the Baldwin-Lomax turbulence model
    subroutine set_blomax_coefficients( set_kappa, set_Aop, set_alpha, set_Ccp, set_Ckleb, set_Cwk) ! Output: None
        implicit none
        real(8), intent(in) :: set_kappa
        real(8), intent(in) :: set_Aop
        real(8), intent(in) :: set_alpha
        real(8), intent(in) :: set_Ccp
        real(8), intent(in) :: set_Ckleb
        real(8), intent(in) :: set_Cwk


        Aop = set_Aop
        Ccp = set_Ccp
        Cwk = set_Cwk

        kappa = set_kappa
        alpha = set_alpha
        Ckleb = set_Ckleb

    end subroutine set_blomax_coefficients

    ! Defines the coefficients of the Baldwin-Lomax turbulence model
    subroutine get_blomax_coefficients( get_kappa, get_Aop, get_alpha, get_Ccp, get_Ckleb, get_Cwk) ! Output: all
        implicit none
        real(8), intent(out) :: get_kappa
        real(8), intent(out) :: get_Aop
        real(8), intent(out) :: get_alpha
        real(8), intent(out) :: get_Ccp
        real(8), intent(out) :: get_Ckleb
        real(8), intent(out) :: get_Cwk

        get_Aop = Aop
        get_Ccp = Ccp
        get_Cwk = Cwk
        get_alpha = alpha
        get_Ckleb = Ckleb
        get_kappa = kappa

    end subroutine get_blomax_coefficients


    ! Allocates and/or initializes the metrics components and vorticity
    ! used in the vorticity calculation
    subroutine initialize_blomax(nx, ny, x, y)
        implicit none
        integer, intent(in) :: nx  ! Number of volumes in csi direction (real+fictitious)
        integer, intent(in) :: ny  ! Number of volumes in eta direction (real+fictitious)
        real(8), dimension (nx*ny), intent(in) :: x ! Coord. x of the northest corner of volume P
        real(8), dimension (nx*ny), intent(in) :: y ! Coord. y of the northest corner of volume P

        allocate( xep(nx*ny), xkp(nx*ny), yep(nx*ny), ykp(nx*ny) )

        allocate( omg(nx*ny) )

        call get_omg_metrics( nx, ny, x, y)

    end subroutine initialize_blomax

    ! Calculates the metrics components used in the vorticity calculation
    subroutine get_omg_metrics( nx, ny, x, y)
        implicit none
        integer, intent(in) :: nx  ! Number of volumes in csi direction (real+fictitious)
        integer, intent(in) :: ny  ! Number of volumes in eta direction (real+fictitious)
        real(8), dimension (nx*ny), intent(in) :: x ! Coord. x of the northest corner of volume P
        real(8), dimension (nx*ny), intent(in) :: y ! Coord. y of the northest corner of volume P

        ! Auxiliary variables

        integer :: i, j, np, nps, npw, npsw

        xkp = 0.d0
        xep = 0.d0
        ykp = 0.d0
        yep = 0.d0

        do j = 2, ny-1
            do i = 2, nx-1

                np   = nx * (j-1) + i
                nps  = np - nx
                npw  = np - 1
                npsw = nps - 1

                xkp(np) = 0.5d0 * ( x(np) + x(nps) - x(npw) - x(npsw) )
                ykp(np) = 0.5d0 * ( y(np) + y(nps) - y(npw) - y(npsw) )
                xep(np) = 0.5d0 * ( x(np) + x(npw) - x(nps) - x(npsw) )
                yep(np) = 0.5d0 * ( y(np) + y(npw) - y(nps) - y(npsw) )

            end do
        end do

    end subroutine get_omg_metrics

    ! Calculates the vorticity at the real volumes centers
    subroutine get_omg( nx, ny, Jp, u, v)
        implicit none
        integer, intent(in) :: nx  ! Number of volumes in csi direction (real+fictitious)
        integer, intent(in) :: ny  ! Number of volumes in eta direction (real+fictitious)
        real(8), dimension (nx*ny), intent(in) :: Jp ! Jacobian at the center of volume P
        real(8), dimension (nx*ny), intent(in) :: u      ! Cartesian velocity of the last iteraction
        real(8), dimension (nx*ny), intent(in) :: v      ! Cartesian velocity of the last iteraction

        ! Auxiliary variables

        integer :: i, j, np, nps, npn, npw, npe

        do j = 2, ny-1
            do i = 2, nx-1

                np   = nx * (j-1) + i
                nps  = np - nx
                npn  = np + nx
                npw  = np - 1
                npe  = np + 1

                omg(np) = 0.5d0 * Jp(np) * abs( &
                    + ( v(npe) - v(npw) ) * yep(np) &
                    - ( v(npn) - v(nps) ) * ykp(np) &
                    + ( u(npe) - u(npw) ) * xep(np) &
                    - ( u(npn) - u(nps) ) * xkp(np) )

            end do
        end do

    end subroutine get_omg

    ! Calculates the dimensionless kinetic eddy viscosity
    subroutine get_nut_plus(       n & !  (in) Number of nodes
        ,                        Udp & !  (in) U_diff_plus
        ,                      yplus & !  (in) Dimensionless distance to the wall y
        ,                    omgplus & !  (in) Dimensionless omega (vorticity)
        ,                    nutplus ) ! (out) Dimensionless kinematic eddy viscosity (nu_T)
        implicit none
        integer, intent(in) :: n   ! Number of nodes
        real(8), intent(in) :: Udp ! U_diff_plus
        real(8), dimension(n), intent(in)  :: yplus   ! Dimensionless distance to the wall y
        real(8), dimension(n), intent(in)  :: omgplus ! Dimensionless omega (vorticity)
        real(8), dimension(n), intent(out) :: nutplus ! Dimensionless kinematic eddy viscosity (nu_T)

        ! Auxiliary variables

        integer :: k      ! Dummy index
        real(8) :: ymaxp  ! y_max_plus
        real(8) :: Fmaxp  ! F_max_plus
        real(8) :: Fwakep ! F_wake_plus
        real(8), dimension(n) :: f     ! Auxiliary vector
        real(8), dimension(n) :: lmp   ! l_mix_plus
        real(8), dimension(n) :: nutip ! nu_T_inner_plus
        real(8), dimension(n) :: nutop ! nu_T_outer_plus

        ! Verifying vorticity

        if ( maxval(abs(omgplus)) < 1.d-15 ) then
            nutplus = 0.d0
            return
        end if

        ! Inner layer

        lmp = kappa * yplus * ( 1.d0 - exp(- yplus / Aop) )

        nutip = lmp * lmp * abs(omgplus)

        ! Outter layer

        f = lmp * abs(omgplus) / kappa

        ! Calculation of Fmaxp and ymaxp with a quadratic interpolation

        k = maxloc( f, 1 )

        if ( k == 1 ) then

            call get_maximum( yplus(1:3), f(1:3), ymaxp, Fmaxp)

        elseif ( k == n ) then

            call get_maximum( yplus(n-2:n), f(n-2:n), ymaxp, Fmaxp)

        else

            call get_maximum( yplus(k-1:k+1), f(k-1:k+1), ymaxp, Fmaxp)

        end if

        ! Uncomment the following two lines for a simpler approximation of Fmaxp and ymaxp (linear interpolation)

        !Fmaxp = maxval( f )

        !ymaxp = yplus( maxloc( f, 1 ) )

        Fwakep = min( ymaxp * Fmaxp,  Cwk * ymaxp * Udp**2 / Fmaxp )

        nutop = alpha * Ccp * Fwakep / ( 1.d0 + 5.5d0 * ( Ckleb * yplus / ymaxp )**6 )

        f = nutop - nutip

        do k = 2, n

            if ( f(k-1) * f(k) <= 0.d0 ) then

                nutplus( 1 : k-1 ) = nutip( 1 : k-1 )

                nutplus( k : n ) = nutop( k : n )

                return

            end if

        end do

        nutplus = nutip

        !write(*,*) "Error. get_nut_plus: it was not possible to find ym_plus."
        !stop

    end subroutine get_nut_plus

    ! Uses a quadratic interpolation to calculate the maximum value of the y(x)
    ! in the prescribed interval: [x1,x3]
    subroutine get_maximum( x, y, xmax, ymax)
        implicit none
        real(8), dimension(3), intent(in) :: x ! Vector of independent variables: x1, x2, x3
        real(8), dimension(3), intent(in) :: y ! Vector of dependent variables: y1, y2, y3

        real(8), intent(out) :: xmax ! x coordinate corresponding to ymax
        real(8), intent(out) :: ymax ! y coordinate corresponding to max(y(x))

        ! Auxiliary variables
        real(8) :: a1, a2 ! Coefficients of the divided differences interpolation formula

        a1 = ( y(2) - y(1) ) / ( x(2) - x(1) )

        a2 = ( ( y(3) - y(2) ) / ( x(3) - x(2) ) - a1 ) / ( x(3) - x(1) )

        if ( a2 < 0.d0 ) then

            xmax = 0.5d0 * ( x(1) + x(2) - a1 / a2 )

            if ( x(1) <= xmax .and. xmax <= x(3) ) then

                ymax = y(1) + a1 * (xmax-x(1)) + a2 * (xmax-x(1)) * (xmax-x(2))

            else

                ymax = maxval(y)

                xmax = x( maxloc(y,1) )

            end if

        else

            ymax = maxval(y)

            xmax = x( maxloc(y,1) )

        end if

    end subroutine get_maximum

    subroutine get_eddy_viscosity_blomax( nx, ny, x, y, xp, yp, xk, yk, Jp, Jn, vln, ro, ron, u, v, vtp) ! Last one is output
        implicit none
        integer, intent(in) :: nx  ! Number of volumes in csi direction (real+fictitious)
        integer, intent(in) :: ny  ! Number of volumes in eta direction (real+fictitious)
        real(8), dimension (nx*ny), intent(in) :: x ! Coord. x of the northest corner of volume P
        real(8), dimension (nx*ny), intent(in) :: y ! Coord. y of the northest corner of volume P
        real(8), dimension (nx*ny), intent(in) :: xp, yp ! Coord. of the centroid of volume P
        real(8), dimension (nx*ny), intent(in) :: xk     ! x_csi at face north of volume P
        real(8), dimension (nx*ny), intent(in) :: yk     ! y_csi at face north of volume P
        real(8), dimension (nx*ny), intent(in) :: Jn     ! Jacobian at the center of north face of volume P
        real(8), dimension (nx*ny), intent(in) :: Jp     ! Jacobian at the center of volume P
        real(8), dimension (nx*ny), intent(in) :: vln    ! Laminar viscosity at center of face north
        real(8), dimension (nx*ny), intent(in) :: ro     ! Specific mass (absolute density) at center of volumes
        real(8), dimension (nx*ny), intent(in) :: ron    ! Absolute density at north face
        real(8), dimension (nx*ny), intent(in) :: u      ! Cartesian velocity of the last iteraction
        real(8), dimension (nx*ny), intent(in) :: v      ! Cartesian velocity of the last iteraction
        real(8), dimension (nx*ny), intent(out) :: vtp   ! Eddy viscosity at center of volume P

        ! Auxiliary variables

        integer :: i, j, k, np, npw, npl

        real(8) :: xw   ! Wall x coordinate
        real(8) :: yw   ! Wall y coordinate
        real(8) :: row  ! Specific mass at wall surface
        real(8) :: vlw  ! Laminar viscosity at wall surface
        real(8) :: nuw  ! kinematic viscosity at wall surface
        real(8) :: tauw ! Shear stress at wall surface
        real(8) :: utau ! Friction velocity
        real(8) :: Udif ! Maximum speed in the velocity profile

        real(8), dimension(ny-2) :: xl   ! x coord. of the eta line
        real(8), dimension(ny-2) :: yl   ! y coord. of the eta line
        real(8), dimension(ny-2) :: dw   ! distance to the wall following the eta line
        real(8), dimension(ny-2) :: omgl ! vorticity following the eta line
        real(8), dimension(ny-2) :: ul   ! velocity component following the eta line
        real(8), dimension(ny-2) :: vl   ! velocity component following the eta line
        real(8), dimension(ny-2) :: nutplus ! dimensionless kinematic eddy viscosity

        ! Vorticity calculation
        call get_omg( nx, ny, Jp, u, v)

        j = ny-1

        do i = 2, nx-1

            np = nx * (j-1) + i

            npw = np - 1

            ! Shear stress
            tauw = abs( vln(np) * Jn(np) * 2.d0 * ( u(np) * xk(np) + v(np) * yk(np) ) )

            ! Specific mass at wall
            row = ron(np)

            ! Laminar viscosity at wall
            vlw = vln(np)

            ! Kinematic viscosity at wall
            nuw = vlw / row

            ! Frinction velocity
            utau = sqrt( abs( tauw / row ) )

            ! Wall x coordinate
            xw = 0.5d0 * ( x(np) + x(npw) )

            ! Wall y coordinate
            yw = 0.5d0 * ( y(np) + y(npw) )

            do k = ny-1, 2, -1

                npl = nx * (k-1) + i

                xl(ny-k) = xp(npl)

                yl(ny-k) = yp(npl)

                omgl(ny-k) = omg(npl)

                ul(ny-k) = u(npl)

                vl(ny-k) = v(npl)

            end do

            ! Maximum speed in the velocity profile following an eta line
            Udif = maxval( sqrt( ul**2 + vl**2 ) )

            ! Distance to the wall following an eta line
            dw = sqrt( (xl-xw)**2 + (yl-yw)**2 )

            ! Calculation of the dimensionless kinematic eddy viscosity
            call get_nut_plus(    ny-2 & !  (in) Number of nodes
            ,              Udif / utau & !  (in) U_diff_plus
            ,          dw * utau / nuw & !  (in) Dimensionless distance to the wall y
            ,     omgl * nuw / utau**2 & !  (in) Dimensionless omega (vorticity)
            ,                  nutplus ) ! (out) Dimensionless kinematic eddy viscosity (nu_T)

            ! Calculation of the eddy viscosity
            do k = ny-1, 2, -1

                npl = nx * (k-1) + i

                vtp(npl) = ro(npl) * nutplus(ny-k) * nuw

            end do

        end do
    end subroutine get_eddy_viscosity_blomax

    subroutine get_eddy_viscosity_fictitious(nx, ny, phi)
        implicit none
        integer, intent(in)  :: nx
        integer, intent(in)  :: ny

        real(8), dimension(nx*ny), intent(inout) :: phi ! Function to be extrapolated to the boundaries

        integer :: i, j
        integer :: np, npe, npw, npn, nps, npsw, npse, npnw, npne

        ! South boundary
        j = 1
        do i = 2, nx-1

            np   = nx * (j-1) + i
            npn  = np + nx

            phi(np) = phi(npn)
        end do

        ! North boundary
        j = ny
        do i = 2, nx-1

            np   = nx * (j-1) + i
            nps  = np - nx

            phi(np) = phi(nps)
        end do

        ! West boundary
        i = 1
        do j = 2, ny-1
            np   = nx * (j-1) + i
            npe  = np + 1

            phi(np) = phi(npe)
        end do

        ! East boundary
        i = nx
        do j = 2, ny-1
            np   = nx * (j-1) + i
            npw  = np - 1

            phi(np) = phi(npw)
        end do

        ! SW corner

        np = 1

        npn  = np + nx

        npne = npn + 1

        phi(np) = phi(npne)

        ! SE corner

        np = nx

        npn  = np + nx

        npnw = npn - 1

        phi(np) = phi(npnw)

        ! NW corner

        i = 1; j = ny

        np   = nx * (j-1) + i

        nps  = np - nx

        npse = nps + 1

        phi(np) = phi(npse)

        ! NE corner

        i = nx; j = ny

        np   = nx * (j-1) + i

        nps  = np - nx

        npsw = nps - 1

        phi(np) = phi(npsw)


    end subroutine get_eddy_viscosity_fictitious

end module blomax

