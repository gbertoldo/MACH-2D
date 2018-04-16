module coefficients
   implicit none

contains

   subroutine get_density_at_nodes( nx, ny, Rg, p, T, ro) ! ro is output
      implicit none
      integer, intent(in) :: nx ! Number of volumes in csi direction (real+fictitious)
      integer, intent(in) :: ny ! Number of volumes in eta direction (real+fictitious)
      real(8), intent(in) :: Rg ! Perfect gas constant
      real(8), dimension(nx*ny), intent(in)  :: p  ! Pressure at center of volumes
      real(8), dimension(nx*ny), intent(in)  :: T  ! Temperature at center of volumes
      real(8), dimension(nx*ny), intent(out) :: ro ! Specific mass (absolute density) at center of volumes

      ro = p / ( Rg * T )

   end subroutine get_density_at_nodes


   subroutine get_density_at_faces( nx, ny, beta, ro, Uce, Vcn, roe, ron) ! roe and ron are output
      implicit none
      integer, intent(in) :: nx   ! Number of volumes in csi direction (real+fictitious)
      integer, intent(in) :: ny   ! Number of volumes in eta direction (real+fictitious)
      real(8), intent(in) :: beta ! Constant of the UDS/CDS mixing scheme
      real(8), dimension(nx*ny), intent(in)  :: ro  ! Specific mass (absolute density) at center of volumes
      real(8), dimension(nx*ny), intent(in)  :: Uce ! Contravariant velocity U at east face
      real(8), dimension(nx*ny), intent(in)  :: Vcn ! Contravariant velocity V at north face
      real(8), dimension(nx*ny), intent(out) :: roe ! Absolute density at east face
      real(8), dimension(nx*ny), intent(out) :: ron ! Absolute density at north face

      integer :: i, j, np, npe, npn
      real(8) :: ae, an ! Coefficients of UDS scheme

      ! Density at east face

      do j = 2, ny-1
         do i = 1, nx-1
            np  = nx * (j-1) + i
            npe = np + 1

            ae = dsign(0.5d0,Uce(np))
            roe(np) = (0.5d0+ae) * ro(np) + (0.5d0-ae) * ro(npe) + beta * ae * ( ro(npe)-ro(np) )

         end do
      end do

      ! Density at north face

      do j = 1, ny-1
         do i = 2, nx-1
            np  = nx * (j-1) + i
            npn = np + nx

            an = dsign(0.5d0,Vcn(np))
            ron(np) = (0.5d0+an) * ro(np) + (0.5d0-an) * ro(npn) + beta * an * ( ro(npn)-ro(np) )

         end do
      end do

   end subroutine get_density_at_faces


   subroutine get_u_coefficients( nx, ny, modvis, dt, rp, re, rn, Jp, Je, Jn &
         ,                        ye, yk, alphae, betae, betan, gamman       &
         ,                        vle, vln, roe, ron, roa, Uce, Vcn, a ) ! a is output
      implicit none
      integer, intent(in) :: nx     ! Number of volumes in csi direction (real+fictitious)
      integer, intent(in) :: ny     ! Number of volumes in eta direction (real+fictitious)
      integer, intent(in) :: modvis ! modvis = 0 -> Euler;  modvis = 1 -> Navier-Stokes
      real(8), intent(in) :: dt     ! Time step
      real(8), dimension (nx*ny), intent(in) :: rp     ! Radius of the center of volume P
      real(8), dimension (nx*ny), intent(in) :: re     ! Radius of the center of east face of volume P
      real(8), dimension (nx*ny), intent(in) :: rn     ! Radius of the center of north face of volume P
      real(8), dimension (nx*ny), intent(in) :: Jp     ! Jacobian at the center of volume P
      real(8), dimension (nx*ny), intent(in) :: Je     ! Jacobian at the center of east face of volume P
      real(8), dimension (nx*ny), intent(in) :: Jn     ! Jacobian at the center of north face of volume P
      real(8), dimension (nx*ny), intent(in) :: ye     ! y_eta at face east of volume P
      real(8), dimension (nx*ny), intent(in) :: yk     ! y_csi at face north of volume P
      real(8), dimension (nx*ny), intent(in) :: Alphae ! (metric) Alpha at the center of east face of volume P
      real(8), dimension (nx*ny), intent(in) :: Betae  ! (metric) Beta  at the center of east face of volume P
      real(8), dimension (nx*ny), intent(in) :: Betan  ! (metric) Beta  at the center of north face of volume P
      real(8), dimension (nx*ny), intent(in) :: Gamman ! (metric) Gamma at the center of north face of volume P
      real(8), dimension (nx*ny), intent(in) :: vle    ! Laminar viscosity at center of face east
      real(8), dimension (nx*ny), intent(in) :: vln    ! Laminar viscosity at center of face north
      real(8), dimension (nx*ny), intent(in) :: roe    ! Absolute density at east face
      real(8), dimension (nx*ny), intent(in) :: ron    ! Absolute density at north face
      real(8), dimension (nx*ny), intent(in) :: roa    ! Absolute density at a time step before at center of vol. P
      real(8), dimension (nx*ny), intent(in) :: Uce    ! Contravariant velocity U at east face
      real(8), dimension (nx*ny), intent(in) :: Vcn    ! Contravariant velocity V at north face

      real(8), dimension (nx*ny,9), intent(out) :: a   ! Coefficients of the linear system

      ! Auxiliary variables
      integer :: i, j, np, nps, npn, npw, npe
      real(8) :: fmw, fme, fms, fmn, mpa
      real(8) :: ae, aw, an, as
      real(8) :: S1p, S2p

      ! Contribution to the coefficients due to the advection term
      if ( modvis == 0 ) then ! Euler
         do j = 2, ny-1
            do i = 2, nx-1

               np   = nx * (j-1) + i
               nps  = np - nx
               npw  = np - 1

               fme = roe(np ) * re(np ) * Uce(np )
               fmw = roe(npw) * re(npw) * Uce(npw)
               fmn = ron(np ) * rn(np ) * Vcn(np )
               fms = ron(nps) * rn(nps) * Vcn(nps)

               mpa = roa(np) * rp(np) / Jp(np)

               as = dsign( 0.5d0, Vcn(nps) )
               an = dsign( 0.5d0, Vcn(np)  )
               aw = dsign( 0.5d0, Uce(npw) )
               ae = dsign( 0.5d0, Uce(np)  )

               a(np,1) = 0.d0 ! SW
               a(np,3) = 0.d0 ! SE
               a(np,7) = 0.d0 ! NW
               a(np,9) = 0.d0 ! NE

               a(np,2) =-fms * (0.5d0+as) ! S
               a(np,4) =-fmw * (0.5d0+aw) ! W
               a(np,6) = fme * (0.5d0-ae) ! E
               a(np,8) = fmn * (0.5d0-an) ! N

               a(np,5) = mpa / dt - ( a(np,8) +  a(np,2) + a(np,4) + a(np,6) )

            end do
         end do

      else ! Navier-Stokes

         do j = 2, ny-1
            do i = 2, nx-1

               np   = nx * (j-1) + i
               nps  = np - nx
               npn  = np + nx
               npw  = np - 1
               npe  = np + 1

               fme = roe(np ) * re(np ) * Uce(np )
               fmw = roe(npw) * re(npw) * Uce(npw)
               fmn = ron(np ) * rn(np ) * Vcn(np )
               fms = ron(nps) * rn(nps) * Vcn(nps)

               mpa = roa(np) * rp(np) / Jp(np)

               as = dsign( 0.5d0, Vcn(nps) )
               an = dsign( 0.5d0, Vcn(np)  )
               aw = dsign( 0.5d0, Uce(npw) )
               ae = dsign( 0.5d0, Uce(np)  )

               ! Contribution to the coefficients due to the advection and diffusion term

               a(np,1) =  ( vle(npw) * re(npw) * Je(npw) * betae(npw)  &
                  +       vln(nps) * rn(nps) * Jn(nps) * betan(nps) ) / 4.d0


               a(np,2) =  - fms * (0.5d0+as) & ! S
               -       vln(nps) * rn(nps) * Jn(nps) * gamman(nps) &
                  +     ( vle(npw) * re(npw) * Je(npw) * betae(npw)  &
                  -       vle(np ) * re(np ) * Je(np ) * betae(np ) ) / 4.d0


               a(np,3) = -( vle(np ) * re(np ) * Je(np ) * betae(np )  &
                  +       vln(nps) * rn(nps) * Jn(nps) * betan(nps) ) / 4.d0


               a(np,4) = -  fmw * (0.5d0+aw) & ! W
               -       vle(npw) * re(npw) * Je(npw) * alphae(npw) &
                  +     ( vln(nps) * rn(nps) * Jn(nps) * betan(nps)  &
                  -       vln(np ) * rn(np ) * Jn(np ) * betan(np ) ) / 4.d0


               a(np,6) =    fme * (0.5d0-ae) & ! E
               -       vle(np ) * re(np ) * Je(np ) * alphae(np ) &
                  +     (-vln(nps) * rn(nps) * Jn(nps) * betan(nps)  &
                  +       vln(np ) * rn(np ) * Jn(np ) * betan(np ) ) / 4.d0


               a(np,7) =  (-vle(npw) * re(npw) * Je(npw) * betae(npw)  &
                  -       vln(np ) * rn(np ) * Jn(np ) * betan(np ) ) / 4.d0


               a(np,8) =    fmn * (0.5d0-an) & ! N
               -       vln(np ) * rn(np ) * Jn(np ) * gamman(np ) &
                  -     ( vle(npw) * re(npw) * Je(npw) * betae(npw)  &
                  -       vle(np ) * re(np ) * Je(np ) * betae(np ) ) / 4.d0


               a(np,9) =  ( vle(np ) * re(np ) * Je(np ) * betae(np )  &
                  +       vln(np ) * rn(np ) * Jn(np ) * betan(np ) ) / 4.d0


               a(np,5) = mpa / dt - ( a(np,8) +  a(np,2) + a(np,4) + a(np,6) )


               ! Contribution to the coefficients due to the source term

               S1p =  rp(np) * Je(np ) * vle(np ) * ye(np ) ** 2 / 3.d0 &
                  + rp(np) * Je(npw) * vle(npw) * ye(npw) ** 2 / 3.d0

               S2p =  rp(np) * Jn(np ) * vln(np ) * yk(np ) ** 2 / 3.d0 &
                  + rp(np) * Jn(nps) * vln(nps) * yk(nps) ** 2 / 3.d0

               a(np,5) = a(np,5) + S1p + S2p

            end do
         end do

      end if

   end subroutine get_u_coefficients


   subroutine get_u_source( nx, ny, modvis, beta, dt, rp, re, rn   &
         ,                  xe, ye, xk, yk, xke, yke, xen, yen     &
         ,                  Jp, Je, Jn, roe, ron, roa, p, vle, vln &
         ,                  Uce, Vcn, ua, u, v, cup, sup, b ) ! The last 3 are output
      implicit none
      integer, intent(in) :: nx     ! Number of volumes in csi direction (real+fictitious)
      integer, intent(in) :: ny     ! Number of volumes in eta direction (real+fictitious)
      integer, intent(in) :: modvis ! modvis = 0 -> Euler;  modvis = 1 -> Navier-Stokes
      real(8), intent(in) :: beta   ! Constant of the UDS/CDS mixing scheme
      real(8), intent(in) :: dt     ! Time step
      real(8), dimension (nx*ny), intent(in) :: rp     ! Radius of the center of volume P
      real(8), dimension (nx*ny), intent(in) :: re     ! Radius of the center of east face of volume P
      real(8), dimension (nx*ny), intent(in) :: rn     ! Radius of the center of north face of volume P
      real(8), dimension (nx*ny), intent(in) :: xe     ! x_eta at center of face east
      real(8), dimension (nx*ny), intent(in) :: ye     ! y_eta at center of face east
      real(8), dimension (nx*ny), intent(in) :: xk     ! x_csi at center of face north
      real(8), dimension (nx*ny), intent(in) :: yk     ! y_csi at center of face north
      real(8), dimension (nx*ny), intent(in) :: xke    ! x_csi at center of face east
      real(8), dimension (nx*ny), intent(in) :: yke    ! y_csi at center of face east
      real(8), dimension (nx*ny), intent(in) :: xen    ! x_eta at center of face north
      real(8), dimension (nx*ny), intent(in) :: yen    ! y_eta at center of face north
      real(8), dimension (nx*ny), intent(in) :: Jp     ! Jacobian at the center of volume P
      real(8), dimension (nx*ny), intent(in) :: Je     ! Jacobian at the center of east face of volume P
      real(8), dimension (nx*ny), intent(in) :: Jn     ! Jacobian at the center of north face of volume P
      real(8), dimension (nx*ny), intent(in) :: roe    ! Absolute density at east face
      real(8), dimension (nx*ny), intent(in) :: ron    ! Absolute density at north face
      real(8), dimension (nx*ny), intent(in) :: roa    ! Absolute density at a time step before at center of vol. P
      real(8), dimension (nx*ny), intent(in) :: p      ! Pressure at center o volume P
      real(8), dimension (nx*ny), intent(in) :: vle    ! Laminar viscosity at center of face east
      real(8), dimension (nx*ny), intent(in) :: vln    ! Laminar viscosity at center of face north
      real(8), dimension (nx*ny), intent(in) :: Uce    ! Contravariant velocity U at east face
      real(8), dimension (nx*ny), intent(in) :: Vcn    ! Contravariant velocity V at north face
      real(8), dimension (nx*ny), intent(in) :: ua     ! Cartesian velocity of a time step before
      real(8), dimension (nx*ny), intent(in) :: u      ! Cartesian velocity of the last iteraction
      real(8), dimension (nx*ny), intent(in) :: v      ! Cartesian velocity of the last iteraction

      real(8), dimension (nx*ny), intent(out) :: cup   ! Term of deferred correction for u
      real(8), dimension (nx*ny), intent(out) :: sup   ! Viscous term for u
      real(8), dimension (nx*ny), intent(out) :: b     ! Source vector of the linear system

      ! Auxiliary variables
      integer :: i, j, np, nps, npn, npw, npe, npsw, npse, npnw, npne
      real(8) :: fmw, fme, fmn, fms, mpa
      real(8) :: as, an, aw, ae
      real(8) :: S1, S2, S3, S4, S5, S6

      do j = 2, ny-1
         do i = 2, nx-1

            np  = nx * (j-1) + i
            nps  = np - nx
            npn  = np + nx
            npw  = np - 1
            npe  = np + 1

            fme = roe(np ) * re(np ) * Uce(np )
            fmw = roe(npw) * re(npw) * Uce(npw)
            fmn = ron(np ) * rn(np ) * Vcn(np )
            fms = ron(nps) * rn(nps) * Vcn(nps)

            mpa = roa(np) * rp(np) / Jp(np)

            as = dsign( 0.5d0, Vcn(nps) )
            an = dsign( 0.5d0, Vcn(np)  )
            aw = dsign( 0.5d0, Uce(npw) )
            ae = dsign( 0.5d0, Uce(np)  )

            ! Contribution to b due to advection (UDS)

            b(np) = mpa * ua(np) / dt

            ! Contribution to b due to advection (Deferred correction)


            cup(np) = - beta * &
               (  fme * ae * ( u(npe) - u(np ) ) &
               -  fmw * aw * ( u(np ) - u(npw) ) &
               +  fmn * an * ( u(npn) - u(np ) ) &
               -  fms * as * ( u(np ) - u(nps) ) &
               )

            b(np) = b(np) + cup(np)


            ! Contribution to b due to pressure term

            b(np) = b(np) +  0.5d0 * rp(np) * ( &
               + yk(np ) * ( p(npn) + p(np) ) &
               - yk(nps) * ( p(nps) + p(np) ) &
               - ye(np ) * ( p(npe) + p(np) ) &
               + ye(npw) * ( p(npw) + p(np) ) &
               )

         end do
      end do

      ! Contribution to b due to viscous term

      if ( modvis == 1 ) then
         do j = 2, ny-1
            do i = 2, nx-1

               np   = nx * (j-1) + i
               nps  = np - nx
               npn  = np + nx
               npw  = np - 1
               npe  = np + 1
               npsw = nps - 1
               npse = nps + 1
               npnw = npn - 1
               npne = npn + 1



               S1 = rp(np) / 3.d0 * (  &

               + Je(np ) * vle(np ) * &

               ( ye(np )**2 * u(npe) - yke(np ) * ye(np ) * ( u(npn) + u(npne) - u(nps) - u(npse) ) / 4.d0 ) &

               - Je(npw) * vle(npw) * &

               (-ye(npw)**2 * u(npw) - yke(npw) * ye(npw) * ( u(npn) + u(npnw) - u(nps) - u(npsw) ) / 4.d0 ) )



               S2 = rp(np) / 3.d0 * ( &

               + Jn(np ) * vln(np ) * &

               ( yk(np )**2 * u(npn) - yk(np ) * yen(np ) * ( u(npe) + u(npne) - u(npw) - u(npnw) ) / 4.d0 ) &

               - Jn(nps) * vln(nps) * &

               (-yk(nps)**2 * u(nps) - yk(nps) * yen(nps) * ( u(npe) + u(npse) - u(npw) - u(npsw) ) / 4.d0 ) )



               S3 =   Je(np ) * re(np ) * vle(np ) * xe(np ) * &

               ( yke(np ) * (v(npn)+v(npne)-v(nps)-v(npse)) / 4.d0 - ye(np ) * (v(npe)-v(np)) ) &

               - Je(npw) * re(npw) * vle(npw) * xe(npw) * &

               ( yke(npw) * (v(npn)+v(npnw)-v(nps)-v(npsw)) / 4.d0 - ye(npw) * (v(np)-v(npw)) )



               S4 =   Jn(np ) * rn(np ) * vln(np ) * xk(np ) * &

               (yen(np ) * (v(npne)+v(npe)-v(npnw)-v(npw)) / 4.d0 - yk(np ) * (v(npn)-v(np)) ) &

               - Jn(nps) * rn(nps) * vln(nps) * xk(nps) * &

               (yen(nps) * (v(npse)+v(npe)-v(npsw)-v(npw)) / 4.d0 - yk(nps) * (v(np)-v(nps)) )


               if ( abs(re(npw)) < 1.d-15 ) then

                  S5 = - 2.d0 * rp(np) / 3.d0 * Je(np ) * vle(np ) / re(np ) * ye(np ) * &

                  ( xke(np ) * ( rp(npn)*v(npn) + rp(npne)*v(npne) - rp(nps)*v(nps) - rp(npse)*v(npse) ) / 4.d0 &

                  - xe(np) * ( rp(npe)*v(npe) - rp(np)*v(np)))

               else if ( abs(re(np)) < 1.d-15 ) then

                  S5 = + 2.d0 * rp(np) / 3.d0 * Je(npw) * vle(npw) / re(npw) * ye(npw) * &

                  ( xke(npw) * ( rp(npn)*v(npn) + rp(npnw)*v(npnw) - rp(nps)*v(nps) - rp(npsw)*v(npsw) ) / 4.d0 &

                  - xe(npw) * ( rp(np)*v(np) - rp(npw)*v(npw)) )

               else

                  S5 = - 2.d0 * rp(np) / 3.d0 * Je(np ) * vle(np ) / re(np ) * ye(np ) * &

                  ( xke(np ) * ( rp(npn)*v(npn) + rp(npne)*v(npne) - rp(nps)*v(nps) - rp(npse)*v(npse) ) / 4.d0 &

                  - xe(np) * ( rp(npe)*v(npe) - rp(np)*v(np))) &

                  + 2.d0 * rp(np) / 3.d0 * Je(npw) * vle(npw) / re(npw) * ye(npw) * &

                  ( xke(npw) * ( rp(npn)*v(npn) + rp(npnw)*v(npnw) - rp(nps)*v(nps) - rp(npsw)*v(npsw) ) / 4.d0 &

                  - xe(npw) * ( rp(np)*v(np) - rp(npw)*v(npw)) )

               end if


               if ( abs(rn(nps)) < 1.d-15 ) then

                  S6 = - 2.d0 * rp(np) / 3.d0 * Jn(np ) * vln(np ) / rn(np ) * yk(np ) * &

                  ( xen(np) * ( rp(npne)*v(npne) + rp(npe)*v(npe) - rp(npnw)*v(npnw) - rp(npw)*v(npw) ) / 4.d0 &

                  - xk(np) * ( rp(npn)*v(npn) - rp(np)*v(np) ) )

               else if ( abs(rn(np)) < 1.d-15 ) then

                  S6 = 2.d0 * rp(np) / 3.d0 * Jn(nps) * vln(nps) / rn(nps) * yk(nps) * &

                  ( xen(nps) * ( rp(npe)*v(npe) + rp(npse)*v(npse) - rp(npw)*v(npw) - rp(npsw)*v(npsw) ) / 4.d0 &

                  - xk(nps) * ( rp(np)*v(np) - rp(nps)*v(nps) ) )

               else

                  S6 = - 2.d0 * rp(np) / 3.d0 * Jn(np ) * vln(np ) / rn(np ) * yk(np ) * &

                  ( xen(np) * ( rp(npne)*v(npne) + rp(npe)*v(npe) - rp(npnw)*v(npnw) - rp(npw)*v(npw) ) / 4.d0 &

                  - xk(np) * ( rp(npn)*v(npn) - rp(np)*v(np) ) ) &

                  + 2.d0 * rp(np) / 3.d0 * Jn(nps) * vln(nps) / rn(nps) * yk(nps) * &

                  ( xen(nps) * ( rp(npe)*v(npe) + rp(npse)*v(npse) - rp(npw)*v(npw) - rp(npsw)*v(npsw) ) / 4.d0 &

                  - xk(nps) * ( rp(np)*v(np) - rp(nps)*v(nps) ) )

               end if

               sup(np) = S1 + S2 + S3 + S4 + S5 + S6

               b(np) = b(np) + sup(np)

            end do
         end do

      end if
   end subroutine get_u_source

   ! Coefficients of the linear system for v (real volumes)
   subroutine get_v_coefficients( nx, ny, coord, modvis, dt, rp, re, rn, Jp, Je, Jn &
         ,                        xe, xk, alphae, betae, betan, gamman       &
         ,                        vle, vln, vlp, roe, ron, roa, Uce, Vcn, a ) ! Last 1 is output
      implicit none
      integer, intent(in) :: nx     ! Number of volumes in csi direction (real+fictitious)
      integer, intent(in) :: ny     ! Number of volumes in eta direction (real+fictitious)
      integer, intent(in) :: coord  ! Coordinate system ( 1 = cylindrical, 0 = cartesian)
      integer, intent(in) :: modvis ! modvis = 0 -> Euler;  modvis = 1 -> Navier-Stokes
      real(8), intent(in) :: dt     ! Time step
      real(8), dimension (nx*ny), intent(in) :: rp     ! Radius of the center of volume P
      real(8), dimension (nx*ny), intent(in) :: re     ! Radius of the center of east face of volume P
      real(8), dimension (nx*ny), intent(in) :: rn     ! Radius of the center of north face of volume P
      real(8), dimension (nx*ny), intent(in) :: Jp     ! Jacobian at the center of volume P
      real(8), dimension (nx*ny), intent(in) :: Je     ! Jacobian at the center of east face of volume P
      real(8), dimension (nx*ny), intent(in) :: Jn     ! Jacobian at the center of north face of volume P
      real(8), dimension (nx*ny), intent(in) :: xe     ! y_eta at face east of volume P
      real(8), dimension (nx*ny), intent(in) :: xk     ! y_csi at face north of volume P
      real(8), dimension (nx*ny), intent(in) :: Alphae ! (metric) Alpha at the center of east face of volume P
      real(8), dimension (nx*ny), intent(in) :: Betae  ! (metric) Beta  at the center of east face of volume P
      real(8), dimension (nx*ny), intent(in) :: Betan  ! (metric) Beta  at the center of north face of volume P
      real(8), dimension (nx*ny), intent(in) :: Gamman ! (metric) Gamma at the center of north face of volume P
      real(8), dimension (nx*ny), intent(in) :: vle    ! Laminar viscosity at center of face east
      real(8), dimension (nx*ny), intent(in) :: vln    ! Laminar viscosity at center of face north
      real(8), dimension (nx*ny), intent(in) :: vlp    ! Laminar viscosity at center of volume P
      real(8), dimension (nx*ny), intent(in) :: roe    ! Absolute density at east face
      real(8), dimension (nx*ny), intent(in) :: ron    ! Absolute density at north face
      real(8), dimension (nx*ny), intent(in) :: roa    ! Absolute density at a time step before at center of vol. P
      real(8), dimension (nx*ny), intent(in) :: Uce    ! Contravariant velocity U at east face
      real(8), dimension (nx*ny), intent(in) :: Vcn    ! Contravariant velocity V at north face

      real(8), dimension (nx*ny,9), intent(out) :: a   ! Coefficients of the linear system

      ! Auxiliary variables
      integer :: i, j, np, nps, npn, npw, npe
      real(8) :: fmw, fme, fms, fmn, mpa
      real(8) :: ae, aw, an, as
      real(8) :: S1p, S2p, S7p

      ! Contribution to the coefficients due to the advection term
      if ( modvis == 0 ) then ! Euler
         do j = 2, ny-1
            do i = 2, nx-1

               np   = nx * (j-1) + i
               nps  = np - nx
               npw  = np - 1

               fme = roe(np ) * re(np ) * Uce(np )
               fmw = roe(npw) * re(npw) * Uce(npw)
               fmn = ron(np ) * rn(np ) * Vcn(np )
               fms = ron(nps) * rn(nps) * Vcn(nps)

               mpa = roa(np) * rp(np) / Jp(np)

               as = dsign( 0.5d0, Vcn(nps) )
               an = dsign( 0.5d0, Vcn(np)  )
               aw = dsign( 0.5d0, Uce(npw) )
               ae = dsign( 0.5d0, Uce(np)  )

               a(np,1) = 0.d0 ! SW
               a(np,3) = 0.d0 ! SE
               a(np,7) = 0.d0 ! NW
               a(np,9) = 0.d0 ! NE

               a(np,2) =-fms * (0.5d0+as) ! S
               a(np,4) =-fmw * (0.5d0+aw) ! W
               a(np,6) = fme * (0.5d0-ae) ! E
               a(np,8) = fmn * (0.5d0-an) ! N

               a(np,5) = mpa / dt - ( a(np,8) +  a(np,2) + a(np,4) + a(np,6) )

            end do
         end do

      else ! Navier-Stokes

         do j = 2, ny-1
            do i = 2, nx-1

               np   = nx * (j-1) + i
               nps  = np - nx
               npn  = np + nx
               npw  = np - 1
               npe  = np + 1

               fme = roe(np ) * re(np ) * Uce(np )
               fmw = roe(npw) * re(npw) * Uce(npw)
               fmn = ron(np ) * rn(np ) * Vcn(np )
               fms = ron(nps) * rn(nps) * Vcn(nps)

               mpa = roa(np) * rp(np) / Jp(np)

               as = dsign( 0.5d0, Vcn(nps) )
               an = dsign( 0.5d0, Vcn(np)  )
               aw = dsign( 0.5d0, Uce(npw) )
               ae = dsign( 0.5d0, Uce(np)  )

               ! Contribution to the coefficients due to the advection and diffusion term

               a(np,1) =  ( vle(npw) * re(npw) * Je(npw) * betae(npw)  &
                  +       vln(nps) * rn(nps) * Jn(nps) * betan(nps) ) / 4.d0


               a(np,2) =  - fms * (0.5d0+as) & ! S
               -       vln(nps) * rn(nps) * Jn(nps) * gamman(nps) &
                  +     ( vle(npw) * re(npw) * Je(npw) * betae(npw)  &
                  -       vle(np ) * re(np ) * Je(np ) * betae(np ) ) / 4.d0


               a(np,3) = -( vle(np ) * re(np ) * Je(np ) * betae(np )  &
                  +       vln(nps) * rn(nps) * Jn(nps) * betan(nps) ) / 4.d0


               a(np,4) = -  fmw * (0.5d0+aw) & ! W
               -       vle(npw) * re(npw) * Je(npw) * alphae(npw) &
                  +     ( vln(nps) * rn(nps) * Jn(nps) * betan(nps)  &
                  -       vln(np ) * rn(np ) * Jn(np ) * betan(np ) ) / 4.d0


               a(np,6) =    fme * (0.5d0-ae) & ! E
               -       vle(np ) * re(np ) * Je(np ) * alphae(np ) &
                  +     (-vln(nps) * rn(nps) * Jn(nps) * betan(nps)  &
                  +       vln(np ) * rn(np ) * Jn(np ) * betan(np ) ) / 4.d0


               a(np,7) =  (-vle(npw) * re(npw) * Je(npw) * betae(npw)  &
                  -       vln(np ) * rn(np ) * Jn(np ) * betan(np ) ) / 4.d0


               a(np,8) =    fmn * (0.5d0-an) & ! N
               -       vln(np ) * rn(np ) * Jn(np ) * gamman(np ) &
                  -     ( vle(npw) * re(npw) * Je(npw) * betae(npw)  &
                  -       vle(np ) * re(np ) * Je(np ) * betae(np ) ) / 4.d0


               a(np,9) =  ( vle(np ) * re(np ) * Je(np ) * betae(np )  &
                  +       vln(np ) * rn(np ) * Jn(np ) * betan(np ) ) / 4.d0


               a(np,5) = mpa / dt - ( a(np,8) +  a(np,2) + a(np,4) + a(np,6) )


               ! Contribution to the coefficients due to the source term

               S1p =  vle(np ) * re(np ) * Je(np ) * xe(np )**2 / 3.d0 &
                  + vle(npw) * re(npw) * Je(npw) * xe(npw)**2 / 3.d0

               S2p =  vln(np ) * rn(np ) * Jn(np ) * xk(np )**2 / 3.d0 &
                  + vln(nps) * rn(nps) * Jn(nps) * xk(nps)**2 / 3.d0

               S7p = coord * vlp(np) * 4.d0 / ( 3.d0 * rp(np) * Jp(np) )

               a(np,5) = a(np,5) + S1p + S2p + S7p

            end do
         end do

      end if

   end subroutine get_v_coefficients

   ! Source of the linear system for v (real volumes)
   subroutine get_v_source( nx, ny, coord, modvis, beta, dt, rp, re, rn   &
         ,                  xe, ye, xk, yk, xke, yke, xen, yen     &
         ,                  Jp, Je, Jn, roe, ron, roa, p, vle, vln &
         ,                  Uce, Vcn, va, u, v, cvp, svp, b ) ! Last 3 are output
      implicit none
      integer, intent(in) :: nx     ! Number of volumes in csi direction (real+fictitious)
      integer, intent(in) :: ny     ! Number of volumes in eta direction (real+fictitious)
      integer, intent(in) :: coord  ! Coordinate system ( 1 = cylindrical, 0 = cartesian)
      integer, intent(in) :: modvis ! modvis = 0 -> Euler;  modvis = 1 -> Navier-Stokes
      real(8), intent(in) :: beta   ! Constant of the UDS/CDS mixing scheme
      real(8), intent(in) :: dt     ! Time step
      real(8), dimension (nx*ny), intent(in) :: rp     ! Radius of the center of volume P
      real(8), dimension (nx*ny), intent(in) :: re     ! Radius of the center of east face of volume P
      real(8), dimension (nx*ny), intent(in) :: rn     ! Radius of the center of north face of volume P
      real(8), dimension (nx*ny), intent(in) :: xe     ! x_eta at center of face east
      real(8), dimension (nx*ny), intent(in) :: ye     ! y_eta at center of face east
      real(8), dimension (nx*ny), intent(in) :: xk     ! x_csi at center of face north
      real(8), dimension (nx*ny), intent(in) :: yk     ! y_csi at center of face north
      real(8), dimension (nx*ny), intent(in) :: xke    ! x_csi at center of face east
      real(8), dimension (nx*ny), intent(in) :: yke    ! y_csi at center of face east
      real(8), dimension (nx*ny), intent(in) :: xen    ! x_eta at center of face north
      real(8), dimension (nx*ny), intent(in) :: yen    ! y_eta at center of face north
      real(8), dimension (nx*ny), intent(in) :: Jp     ! Jacobian at the center of volume P
      real(8), dimension (nx*ny), intent(in) :: Je     ! Jacobian at the center of east face of volume P
      real(8), dimension (nx*ny), intent(in) :: Jn     ! Jacobian at the center of north face of volume P
      real(8), dimension (nx*ny), intent(in) :: roe    ! Absolute density at east face
      real(8), dimension (nx*ny), intent(in) :: ron    ! Absolute density at north face
      real(8), dimension (nx*ny), intent(in) :: roa    ! Absolute density at a time step before at center of vol. P
      real(8), dimension (nx*ny), intent(in) :: p      ! Pressure at center o volume P
      real(8), dimension (nx*ny), intent(in) :: vle    ! Laminar viscosity at center of face east
      real(8), dimension (nx*ny), intent(in) :: vln    ! Laminar viscosity at center of face north
      real(8), dimension (nx*ny), intent(in) :: Uce    ! Contravariant velocity U at east face
      real(8), dimension (nx*ny), intent(in) :: Vcn    ! Contravariant velocity V at north face
      real(8), dimension (nx*ny), intent(in) :: va     ! Cartesian velocity of a time step before
      real(8), dimension (nx*ny), intent(in) :: u      ! Cartesian velocity of the last iteraction
      real(8), dimension (nx*ny), intent(in) :: v      ! Cartesian velocity of the last iteraction

      real(8), dimension (nx*ny), intent(out) :: cvp   ! Term of deferred correction for v
      real(8), dimension (nx*ny), intent(out) :: svp   ! Viscous term for v
      real(8), dimension (nx*ny), intent(out) :: b     ! Source vector of the linear system

      ! Auxiliary variables
      integer :: i, j, np, nps, npn, npw, npe, npsw, npse, npnw, npne
      real(8) :: fmw, fme, fmn, fms, mpa
      real(8) :: as, an, aw, ae
      real(8) :: S1, S2, S3, S4, S5, S6, S8

      do j = 2, ny-1
         do i = 2, nx-1

            np  = nx * (j-1) + i
            nps = np - nx
            npn = np + nx
            npw = np - 1
            npe = np + 1

            fme = roe(np ) * re(np ) * Uce(np )
            fmw = roe(npw) * re(npw) * Uce(npw)
            fmn = ron(np ) * rn(np ) * Vcn(np )
            fms = ron(nps) * rn(nps) * Vcn(nps)

            mpa = roa(np) * rp(np) / Jp(np)

            as = dsign( 0.5d0, Vcn(nps) )
            an = dsign( 0.5d0, Vcn(np)  )
            aw = dsign( 0.5d0, Uce(npw) )
            ae = dsign( 0.5d0, Uce(np)  )

            ! Contribution to b due to advection (UDS)

            b(np) = mpa * va(np) / dt

            ! Contribution to b due to advection (Deferred correction)


            cvp(np) = - beta * &
               (  fme * ae * ( v(npe) - v(np ) ) &
               -  fmw * aw * ( v(np ) - v(npw) ) &
               +  fmn * an * ( v(npn) - v(np ) ) &
               -  fms * as * ( v(np ) - v(nps) ) &
               )

            b(np) = b(np) + cvp(np)


            ! Contribution to b due to pressure term

            b(np) = b(np) +  0.5d0 * rp(np) * ( &
               + xe(np ) * ( p(np) + p(npe) ) &
               - xe(npw) * ( p(np) + p(npw) ) &
               - xk(np ) * ( p(np) + p(npn) ) &
               + xk(nps) * ( p(np) + p(nps) ) &
               )

         end do
      end do

      ! Contribution to b due to viscous term

      if ( modvis == 1 ) then

         do j = 2, ny-1
            do i = 2, nx-1

               np   = nx * (j-1) + i
               nps  = np - nx
               npn  = np + nx
               npw  = np - 1
               npe  = np + 1
               npsw = nps - 1
               npse = nps + 1
               npnw = npn - 1
               npne = npn + 1


               S1 =   vle(np ) * re(np ) * Je(np ) / 3.d0 * ( xe(np )**2 * v(npe) &
                  - xke(np ) * xe(np ) * ( v(npn) + v(npne) - v(nps) - v(npse) ) / 4.d0 ) &
                  - vle(npw) * re(npw) * Je(npw) / 3.d0 * (-xe(npw)**2 * v(npw) &
                  - xke(npw) * xe(npw) * ( v(npn) + v(npnw) - v(nps) - v(npsw) ) / 4.d0 )

               S2 =   vln(np ) * rn(np ) * Jn(np ) / 3.d0 * ( xk(np )**2 * v(npn) &
                  - xk(np ) * xen(np ) * ( v(npe) + v(npne) - v(npnw) - v(npw) ) / 4.d0 ) &
                  - vln(nps) * rn(nps) * Jn(nps) / 3.d0 * (-xk(nps)**2 * v(nps) &
                  - xk(nps) * xen(nps) * ( v(npe) + v(npse) - v(npsw) - v(npw) ) / 4.d0 )

               S3 =   rp(np) * vle(np ) * Je(np ) * ye(np ) * ( &
                  + xke(np ) * ( u(npn) + u(npne) - u(nps) - u(npse) ) / 4.d0 &
                  - xe(np ) * ( u(npe) - u(np) ) ) &
                  - rp(np) * vle(npw) * Je(npw) * ye(npw) * ( &
                  + xke(npw) * ( u(npn) + u(npnw) - u(nps) - u(npsw) ) / 4.d0 &
                  - xe(npw) * ( u(np) - u(npw) ) )

               S4 =   rp(np) * vln(np ) * Jn(np ) * yk(np ) * ( &
                  + xen(np ) * ( u(npne) + u(npe) - u(npnw) - u(npw) ) / 4.d0 &
                  - xk(np ) * ( u(npn) - u(np) ) ) &
                  - rp(np) * vln(nps) * Jn(nps) * yk(nps) * ( &
                  + xen(nps) * ( u(npse) + u(npe) - u(npsw) - u(npw) ) / 4.d0 &
                  - xk(nps) * ( u(np) - u(nps) ) )

               S5 = - 2.d0 / 3.d0 * rp(np) * vle(np ) * Je(np ) * xe(np ) * ( &
                  + yke(np ) * ( u(npn) + u(npne) - u(nps) - u(npse) ) / 4.d0 &
                  - ye(np ) * ( u(npe) - u(np) ) ) &
                  + 2.d0 / 3.d0 * rp(np) * vle(npw) * Je(npw) * xe(npw) * ( &
                  + yke(npw) * ( u(npn) + u(npnw) - u(nps) - u(npsw) ) / 4.d0 &
                  - ye(npw) * ( u(np) - u(npw) ) )

               S6 = - 2.d0 / 3.d0 * rp(np) * vln(np ) * Jn(np ) * xk(np ) * ( &
                  + yen(np ) * ( u(npe) + u(npne) - u(npw) - u(npnw) ) / 4.d0 &
                  - yk(np ) * ( u(npn) - u(np) ) ) &
                  + 2.d0 / 3.d0 * rp(np) * vln(nps) * Jn(nps) * xk(nps) * ( &
                  + yen(nps) * ( u(npe) + u(npse) - u(npw) - u(npsw) ) / 4.d0 &
                  - yk(nps) * ( u(np) - u(nps) ) )

               S8 = 2.d0/3.d0 * coord * v(np) * ( &
                  - vln(np) * xk(np) + vln(nps) * xk(nps) &
                  + vle(np) * xe(np) - vle(npw) * xe(npw) &
                  )

               svp(np) = S1 + S2 + S3 + S4 + S5 + S6 + S8

               b(np) = b(np) + svp(np)

            end do
         end do

      end if

   end subroutine get_v_source


   !> \brief Calculates the coefficients and source of the linear system for the temperature
   subroutine get_T_coefficients_and_source( nx, ny, modvis, beta, dt, rp, re, rn       &
         ,                  xe, ye, xk, yk, xke, yke, xen, yen, alphae, betae, betan, gamman &
         ,                  Jp, Je, Jn, roe, ron, ro, roa, p, pa, cp, vle, vln, ke, kn, Uce, Vcn &
         ,                  u, v, ue, ve, un, vn, ua, va, T, Ta, fbe, fbn, a, b ) ! Last 2 are output
      implicit none
      integer, intent(in) :: nx     !< Number of volumes in csi direction (real+fictitious)
      integer, intent(in) :: ny     !< Number of volumes in eta direction (real+fictitious)
      integer, intent(in) :: modvis !< modvis = 0 -> Euler;  modvis = 1 -> Navier-Stokes
      real(8), intent(in) :: beta   !< Constant of the UDS/CDS mixing scheme
      real(8), intent(in) :: dt     !< Time step
      real(8), dimension (nx*ny), intent(in) :: rp     !< Radius of the center of volume P
      real(8), dimension (nx*ny), intent(in) :: re     !< Radius of the center of east face of volume P
      real(8), dimension (nx*ny), intent(in) :: rn     !< Radius of the center of north face of volume P
      real(8), dimension (nx*ny), intent(in) :: xe     !< x_eta at center of face east
      real(8), dimension (nx*ny), intent(in) :: ye     !< y_eta at center of face east
      real(8), dimension (nx*ny), intent(in) :: xk     !< x_csi at center of face north
      real(8), dimension (nx*ny), intent(in) :: yk     !< y_csi at center of face north
      real(8), dimension (nx*ny), intent(in) :: xke    ! x_csi at center of face east
      real(8), dimension (nx*ny), intent(in) :: yke    ! y_csi at center of face east
      real(8), dimension (nx*ny), intent(in) :: xen    ! x_eta at center of face north
      real(8), dimension (nx*ny), intent(in) :: yen    ! y_eta at center of face north
      real(8), dimension (nx*ny), intent(in) :: Alphae !< (metric) Alpha at the center of east face of volume P
      real(8), dimension (nx*ny), intent(in) :: Betae  !< (metric) Beta  at the center of east face of volume P
      real(8), dimension (nx*ny), intent(in) :: Betan  !< (metric) Beta  at the center of north face of volume P
      real(8), dimension (nx*ny), intent(in) :: Gamman !< (metric) Gamma at the center of north face of volume P
      real(8), dimension (nx*ny), intent(in) :: Jp     !< Jacobian at the center of volume P
      real(8), dimension (nx*ny), intent(in) :: Je     !< Jacobian at the center of east face of volume P
      real(8), dimension (nx*ny), intent(in) :: Jn     !< Jacobian at the center of north face of volume P
      real(8), dimension (nx*ny), intent(in) :: roe    !< Absolute density at east face
      real(8), dimension (nx*ny), intent(in) :: ron    !< Absolute density at north face
      real(8), dimension (nx*ny), intent(in) :: ro     !< Absolute density at center of vol. P
      real(8), dimension (nx*ny), intent(in) :: roa    !< Absolute density at a time step before at center of vol. P
      real(8), dimension (nx*ny), intent(in) :: p      !< Pressure at center o volume P
      real(8), dimension (nx*ny), intent(in) :: pa     !< Pressure at center o volume P at a time step before
      real(8), dimension (nx*ny), intent(in) :: cp     !< Specific heat at constant pressure at center o volume P
      real(8), dimension (nx*ny), intent(in) :: vle    ! Laminar viscosity at center of face east
      real(8), dimension (nx*ny), intent(in) :: vln    ! Laminar viscosity at center of face north
      real(8), dimension (nx*ny), intent(in) :: ke     !< Thermal conductivity at center of face east
      real(8), dimension (nx*ny), intent(in) :: kn     !< Thermal conductivity at center of face north
      real(8), dimension (nx*ny), intent(in) :: Uce    !< Contravariant velocity U at east face
      real(8), dimension (nx*ny), intent(in) :: Vcn    !< Contravariant velocity V at north face
      real(8), dimension (nx*ny), intent(in) :: u      !< Cartesian velocity of the last iteraction
      real(8), dimension (nx*ny), intent(in) :: v      !< Cartesian velocity of the last iteraction
      real(8), dimension (nx*ny), intent(in) :: ue     !< Cartesian velocity u at center of east face
      real(8), dimension (nx*ny), intent(in) :: ve     !< Cartesian velocity v at center of east face
      real(8), dimension (nx*ny), intent(in) :: un     !< Cartesian velocity u at center of north face
      real(8), dimension (nx*ny), intent(in) :: vn     !< Cartesian velocity v at center of north face
      real(8), dimension (nx*ny), intent(in) :: ua     !< Cartesian velocity of the previous time step
      real(8), dimension (nx*ny), intent(in) :: va     !< Cartesian velocity of the previous time step
      real(8), dimension (nx*ny), intent(in) :: T      !< Temperature of the last iteraction
      real(8), dimension (nx*ny), intent(in) :: Ta     !< Temperature at the time step before
      real(8), dimension (nx*ny), intent(in) :: fbe    !< Face of boundary east (1 if an east boundary, 0 otherwise)
      real(8), dimension (nx*ny), intent(in) :: fbn    !< Face of boundary north (1 if an north boundary, 0 otherwise)

      real(8), dimension (nx*ny,9), intent(out) :: a   !< Coefficients of the linear system
      real(8), dimension (nx*ny),   intent(out) :: b   !< Source vector of the linear system

      ! Constant

      real(8), parameter :: eps = epsilon(1.d0)


      ! Auxiliary variables
      integer :: i, j
      integer :: np, nps, npn, npw, npe, npnw, npne, npsw, npse
      real(8) :: fmw, fme, fms, fmn, mpa
      real(8) :: as, an, aw, ae
      real(8) :: KpE, KpW, KpN, KpS, KpP, KpPa
      real(8) :: Kfe, Kfw, Kfn, Kfs, raux
      real(8) :: Kfbe, Kfbw, Kfbn, Kfbs
      real(8) :: dudk_e, dudk_w, dudk_n, dudk_s
      real(8) :: dude_e, dude_w, dude_n, dude_s
      real(8) :: dvdk_e, dvdk_w, dvdk_n, dvdk_s
      real(8) :: dvde_e, dvde_w, dvde_n, dvde_s
      real(8) :: dudx_e, dudx_w, dudx_n, dudx_s
      real(8) :: dudy_e, dudy_w, dudy_n, dudy_s
      real(8) :: dvdx_e, dvdx_w, dvdx_n, dvdx_s
      real(8) :: dvdy_e, dvdy_w, dvdy_n, dvdy_s
      real(8) :: divU_e, divU_w, divU_n, divU_s
      real(8) :: Txx_e, Txx_w, Txx_n, Txx_s
      real(8) :: Txy_e, Txy_w, Txy_n, Txy_s
      real(8) :: Tyy_e, Tyy_w, Tyy_n, Tyy_s
      real(8) :: Lbd_x_e, Lbd_x_w, Lbd_x_n, Lbd_x_s
      real(8) :: Lbd_y_e, Lbd_y_w, Lbd_y_n, Lbd_y_s
      real(8) :: Lbd_k_e, Lbd_k_w, Lbd_e_n, Lbd_e_s



      if ( modvis == 0 ) then ! Euler

         do j = 2, ny-1
            do i = 2, nx-1

               np  = nx * (j-1) + i
               nps  = np - nx
               npn  = np + nx
               npw  = np - 1
               npe  = np + 1

               fme = roe(np ) * re(np ) * Uce(np )
               fmw = roe(npw) * re(npw) * Uce(npw)
               fmn = ron(np ) * rn(np ) * Vcn(np )
               fms = ron(nps) * rn(nps) * Vcn(nps)

               mpa = roa(np) * rp(np) / Jp(np)

               as = dsign( 0.5d0, Vcn(nps) )
               an = dsign( 0.5d0, Vcn(np)  )
               aw = dsign( 0.5d0, Uce(npw) )
               ae = dsign( 0.5d0, Uce(np)  )

               ! Contribution to the COEFFICIENTS due to advection (UDS)

               a(np,1) = 0.d0 ! SW
               a(np,3) = 0.d0 ! SE
               a(np,7) = 0.d0 ! NW
               a(np,9) = 0.d0 ! NE

               a(np,2) =-fms * (0.5d0+as) * cp(np) ! S
               a(np,4) =-fmw * (0.5d0+aw) * cp(np) ! W
               a(np,6) = fme * (0.5d0-ae) * cp(np) ! E
               a(np,8) = fmn * (0.5d0-an) * cp(np) ! N

               a(np,5) = mpa  * cp(np) / dt - ( a(np,8) +  a(np,2) + a(np,4) + a(np,6) )

               ! Contribution to the SOURCE due to advection (UDS)

               b(np) = mpa * Ta(np) * cp(np) / dt

               ! Contribution to the SOURCE due to advection (Deferred correction)

               b(np) = b(np) - beta * cp(np) *        &
                  (  fme * ae * ( T(npe) - T(np ) ) &
                  -  fmw * aw * ( T(np ) - T(npw) ) &
                  +  fmn * an * ( T(npn) - T(np ) ) &
                  -  fms * as * ( T(np ) - T(nps) ) &
                  )

               ! Contribution to the SOURCE due to pressure term


               ! Values of K at nodes
               KpE = ( u(npe)**2 + v(npe)**2 ) / 2.d0
               KpW = ( u(npw)**2 + v(npw)**2 ) / 2.d0
               KpN = ( u(npn)**2 + v(npn)**2 ) / 2.d0
               KpS = ( u(nps)**2 + v(nps)**2 ) / 2.d0
               KpP = ( u(np )**2 + v(np )**2 ) / 2.d0
               KpPa= ( ua(np)**2 + va(np)**2 ) / 2.d0


               ! Values of K on the boundaries interpolated with K at nodes using UDS deferred to CDS
               Kfe = (0.5d0+ae) * KpP + (0.5d0-ae) * KpE + beta * ae * (KpE-KpP)
               Kfw = (0.5d0+aw) * KpW + (0.5d0-aw) * KpP + beta * aw * (KpP-KpW)
               Kfn = (0.5d0+an) * KpP + (0.5d0-an) * KpN + beta * an * (KpN-KpP)
               Kfs = (0.5d0+as) * KpS + (0.5d0-as) * KpP + beta * as * (KpP-KpS)


               ! Values of K on the boundaries calculated with the boundary velocities
               Kfbe = ( ue(np )**2 + ve(np )**2 ) / 2.d0
               Kfbw = ( ue(npw)**2 + ve(npw)**2 ) / 2.d0
               Kfbn = ( un(np )**2 + vn(np )**2 ) / 2.d0
               Kfbs = ( un(nps)**2 + vn(nps)**2 ) / 2.d0


               ! Choosing between interpolated and boundary values in order
               ! to ensure the conservation of H (total enthalpy) for inviscid flows
               Kfe = Kfe * ( 1.d0 - fbe(np ) ) + Kfbe * fbe(np )
               Kfw = Kfw * ( 1.d0 - fbe(npw) ) + Kfbw * fbe(npw)
               Kfn = Kfn * ( 1.d0 - fbn(np ) ) + Kfbn * fbn(np )
               Kfs = Kfs * ( 1.d0 - fbn(nps) ) + Kfbs * fbn(nps)


               raux = rp(np) / Jp(np) * ( ro(np) * KpP - roa(np) * KpPa ) / dt &

                  + fme * Kfe - fmw * Kfw + fmn * Kfn - fms * Kfs


               b(np) = b(np) + rp(np) / Jp(np) * ( p(np) - pa(np) ) / dt - raux


            end do

         end do

      else ! Navier-Stokes

         do j = 2, ny-1

            do i = 2, nx-1

               np   = nx * (j-1) + i
               nps  = np - nx
               npn  = np + nx
               npw  = np - 1
               npe  = np + 1
               npsw = nps - 1
               npse = nps + 1
               npnw = npn - 1
               npne = npn + 1

               fme = roe(np ) * re(np ) * Uce(np )
               fmw = roe(npw) * re(npw) * Uce(npw)
               fmn = ron(np ) * rn(np ) * Vcn(np )
               fms = ron(nps) * rn(nps) * Vcn(nps)

               mpa = roa(np) * rp(np) / Jp(np)

               as = dsign( 0.5d0, Vcn(nps) )
               an = dsign( 0.5d0, Vcn(np)  )
               aw = dsign( 0.5d0, Uce(npw) )
               ae = dsign( 0.5d0, Uce(np)  )

               ! Contribution to the COEFFICIENTS due to advection-diffusion (UDS-CDS)

               a(np,1) =  ( ke(npw) * re(npw) * Je(npw) * betae(npw)  &
                  +       kn(nps) * rn(nps) * Jn(nps) * betan(nps) ) / 4.d0


               a(np,2) =  - fms * (0.5d0+as) * cp(np) & ! S
               -       kn(nps) * rn(nps) * Jn(nps) * gamman(nps) &
                  +     ( ke(npw) * re(npw) * Je(npw) * betae(npw)  &
                  -       ke(np ) * re(np ) * Je(np ) * betae(np ) ) / 4.d0


               a(np,3) = -( ke(np ) * re(np ) * Je(np ) * betae(np )  &
                  +       kn(nps) * rn(nps) * Jn(nps) * betan(nps) ) / 4.d0


               a(np,4) = -  fmw * (0.5d0+aw) * cp(np) & ! W
               -       ke(npw) * re(npw) * Je(npw) * alphae(npw) &
                  +     ( kn(nps) * rn(nps) * Jn(nps) * betan(nps)  &
                  -       kn(np ) * rn(np ) * Jn(np ) * betan(np ) ) / 4.d0


               a(np,6) =    fme * (0.5d0-ae) * cp(np) & ! E
               -       ke(np ) * re(np ) * Je(np ) * alphae(np ) &
                  +     (-kn(nps) * rn(nps) * Jn(nps) * betan(nps)  &
                  +       kn(np ) * rn(np ) * Jn(np ) * betan(np ) ) / 4.d0


               a(np,7) =  (-ke(npw) * re(npw) * Je(npw) * betae(npw)  &
                  -       kn(np ) * rn(np ) * Jn(np ) * betan(np ) ) / 4.d0


               a(np,8) =    fmn * (0.5d0-an) * cp(np) & ! N
               -       kn(np ) * rn(np ) * Jn(np ) * gamman(np ) &
                  -     ( ke(npw) * re(npw) * Je(npw) * betae(npw)  &
                  -       ke(np ) * re(np ) * Je(np ) * betae(np ) ) / 4.d0


               a(np,9) =  ( ke(np ) * re(np ) * Je(np ) * betae(np )  &
                  +       kn(np ) * rn(np ) * Jn(np ) * betan(np ) ) / 4.d0


               a(np,5) = mpa  * cp(np) / dt - ( a(np,8) +  a(np,2) + a(np,4) + a(np,6) )


               ! Contribution to the SOURCE due to advection (UDS)

               b(np) = mpa * Ta(np) * cp(np) / dt

               ! Contribution to the SOURCE due to advection (Deferred correction)

               b(np) = b(np) - beta * cp(np) *        &
                  (  fme * ae * ( T(npe) - T(np ) ) &
                  -  fmw * aw * ( T(np ) - T(npw) ) &
                  +  fmn * an * ( T(npn) - T(np ) ) &
                  -  fms * as * ( T(np ) - T(nps) ) &
                  )

               ! Contribution to the SOURCE due to pressure term


               ! Values of K at nodes
               KpE = ( u(npe)**2 + v(npe)**2 ) / 2.d0
               KpW = ( u(npw)**2 + v(npw)**2 ) / 2.d0
               KpN = ( u(npn)**2 + v(npn)**2 ) / 2.d0
               KpS = ( u(nps)**2 + v(nps)**2 ) / 2.d0
               KpP = ( u(np )**2 + v(np )**2 ) / 2.d0
               KpPa= ( ua(np)**2 + va(np)**2 ) / 2.d0


               ! Values of K on the boundaries interpolated with K at nodes using UDS deferred to CDS
               Kfe = (0.5d0+ae) * KpP + (0.5d0-ae) * KpE + beta * ae * (KpE-KpP)
               Kfw = (0.5d0+aw) * KpW + (0.5d0-aw) * KpP + beta * aw * (KpP-KpW)
               Kfn = (0.5d0+an) * KpP + (0.5d0-an) * KpN + beta * an * (KpN-KpP)
               Kfs = (0.5d0+as) * KpS + (0.5d0-as) * KpP + beta * as * (KpP-KpS)


               ! Values of K on the boundaries calculated with the boundary velocities
               Kfbe = ( ue(np )**2 + ve(np )**2 ) / 2.d0
               Kfbw = ( ue(npw)**2 + ve(npw)**2 ) / 2.d0
               Kfbn = ( un(np )**2 + vn(np )**2 ) / 2.d0
               Kfbs = ( un(nps)**2 + vn(nps)**2 ) / 2.d0


               ! Choosing between interpolated and boundary values in order
               ! to ensure the conservation of H (total enthalpy) for inviscid flows
               Kfe = Kfe * ( 1.d0 - fbe(np ) ) + Kfbe * fbe(np )
               Kfw = Kfw * ( 1.d0 - fbe(npw) ) + Kfbw * fbe(npw)
               Kfn = Kfn * ( 1.d0 - fbn(np ) ) + Kfbn * fbn(np )
               Kfs = Kfs * ( 1.d0 - fbn(nps) ) + Kfbs * fbn(nps)


               raux = rp(np) / Jp(np) * ( ro(np) * KpP - roa(np) * KpPa ) / dt &

                  + fme * Kfe - fmw * Kfw + fmn * Kfn - fms * Kfs


               b(np) = b(np) + rp(np) / Jp(np) * ( p(np) - pa(np) ) / dt - raux

               ! Contribution to the SOURCE due to the viscous term

               dudk_e = u(npe) - u(np )
               dudk_w = u(np ) - u(npw)
               dudk_n = ( u(npne) + u(npe) - u(npnw) - u(npw) ) / 4.d0
               dudk_s = ( u(npse) + u(npe) - u(npsw) - u(npw) ) / 4.d0


               dude_e = ( u(npne) + u(npn) - u(npse) - u(nps) ) / 4.d0
               dude_w = ( u(npnw) + u(npn) - u(npsw) - u(nps) ) / 4.d0
               dude_n = u(npn) - u(np )
               dude_s = u(np ) - u(nps)


               dvdk_e = v(npe) - v(np )
               dvdk_w = v(np ) - v(npw)
               dvdk_n = ( v(npne) + v(npe) - v(npnw) - v(npw) ) / 4.d0
               dvdk_s = ( v(npse) + v(npe) - v(npsw) - v(npw) ) / 4.d0


               dvde_e = ( v(npne) + v(npn) - v(npse) - v(nps) ) / 4.d0
               dvde_w = ( v(npnw) + v(npn) - v(npsw) - v(nps) ) / 4.d0
               dvde_n = v(npn) - v(np )
               dvde_s = v(np ) - v(nps)


               dudx_e = Je(np ) * ( ye(np )  * dudk_e - yke(np ) * dude_e )
               dudx_w = Je(npw) * ( ye(npw)  * dudk_w - yke(npw) * dude_w )
               dudx_n = Jn(np ) * ( yen(np ) * dudk_n - yk(np )  * dude_n )
               dudx_s = Jn(nps) * ( yen(nps) * dudk_s - yk(nps)  * dude_s )


               dudy_e = Je(np ) * ( xke(np ) * dude_e - xe(np )  * dudk_e )
               dudy_w = Je(npw) * ( xke(npw) * dude_w - xe(npw)  * dudk_w )
               dudy_n = Jn(np ) * ( xk(np )  * dude_n - xen(np ) * dudk_n )
               dudy_s = Jn(nps) * ( xk(nps)  * dude_s - xen(nps) * dudk_s )


               dvdx_e = Je(np ) * ( ye(np )  * dvdk_e - yke(np ) * dvde_e )
               dvdx_w = Je(npw) * ( ye(npw)  * dvdk_w - yke(npw) * dvde_w )
               dvdx_n = Jn(np ) * ( yen(np ) * dvdk_n - yk(np )  * dvde_n )
               dvdx_s = Jn(nps) * ( yen(nps) * dvdk_s - yk(nps)  * dvde_s )


               dvdy_e = Je(np ) * ( xke(np ) * dvde_e - xe(np )  * dvdk_e )
               dvdy_w = Je(npw) * ( xke(npw) * dvde_w - xe(npw)  * dvdk_w )
               dvdy_n = Jn(np ) * ( xk(np )  * dvde_n - xen(np ) * dvdk_n )
               dvdy_s = Jn(nps) * ( xk(nps)  * dvde_s - xen(nps) * dvdk_s )


               divU_e = dudx_e + dvdy_e + ve(np ) / ( abs(re(np )) + eps )
               divU_w = dudx_w + dvdy_w + ve(npw) / ( abs(re(npw)) + eps )
               divU_n = dudx_n + dvdy_n + vn(np ) / ( abs(rn(np )) + eps )
               divU_s = dudx_s + dvdy_s + vn(nps) / ( abs(rn(nps)) + eps )


               Txx_e = -2.d0 * vle(np ) * dudx_e + 2.d0/3.d0 * vle(np ) * divU_e
               Txx_w = -2.d0 * vle(npw) * dudx_w + 2.d0/3.d0 * vle(npw) * divU_w
               Txx_n = -2.d0 * vln(np ) * dudx_n + 2.d0/3.d0 * vln(np ) * divU_n
               Txx_s = -2.d0 * vln(nps) * dudx_s + 2.d0/3.d0 * vln(nps) * divU_s


               Txy_e = - vle(np ) * ( dvdx_e + dudy_e )
               Txy_w = - vle(npw) * ( dvdx_w + dudy_w )
               Txy_n = - vln(np ) * ( dvdx_n + dudy_n )
               Txy_s = - vln(nps) * ( dvdx_s + dudy_s )


               Tyy_e = -2.d0 * vle(np ) * dvdy_e + 2.d0/3.d0 * vle(np ) * divU_e
               Tyy_w = -2.d0 * vle(npw) * dvdy_w + 2.d0/3.d0 * vle(npw) * divU_w
               Tyy_n = -2.d0 * vln(np ) * dvdy_n + 2.d0/3.d0 * vln(np ) * divU_n
               Tyy_s = -2.d0 * vln(nps) * dvdy_s + 2.d0/3.d0 * vln(nps) * divU_s


               Lbd_x_e = ve(np ) * Txy_e + ue(np ) * Txx_e
               Lbd_x_w = ve(npw) * Txy_w + ue(npw) * Txx_w
               Lbd_x_n = vn(np ) * Txy_n + un(np ) * Txx_n
               Lbd_x_s = vn(nps) * Txy_s + un(nps) * Txx_s


               Lbd_y_e = ve(np ) * Tyy_e + ue(np ) * Txy_e
               Lbd_y_w = ve(npw) * Tyy_w + ue(npw) * Txy_w
               Lbd_y_n = vn(np ) * Tyy_n + un(np ) * Txy_n
               Lbd_y_s = vn(nps) * Tyy_s + un(nps) * Txy_s


               Lbd_k_e = Lbd_x_e * ye(np ) - Lbd_y_e * xe(np )
               Lbd_k_w = Lbd_x_w * ye(npw) - Lbd_y_w * xe(npw)
               Lbd_e_n = Lbd_y_n * xk(np ) - Lbd_x_n * yk(np )
               Lbd_e_s = Lbd_y_s * xk(nps) - Lbd_x_s * yk(nps)


               raux = re(np) * Lbd_k_e - re(npw) * Lbd_k_w &
                  +   rn(np) * Lbd_e_n - rn(nps) * Lbd_e_s


               b(np) = b(np) - raux


            end do

         end do

      end if

   end subroutine get_T_coefficients_and_source


   !> \brief Calculates the coefficients of the linear system for pressure correction
   subroutine get_p_coefficients(nx, ny, dt, rp, re, rn, Jp, Uce, Vcn &
         ,                                      roe, ron, g, de, dn, a) ! Output: last one
      implicit none
      integer, intent(in) :: nx     !< Number of volumes in the csi direction (real+fictitious)
      integer, intent(in) :: ny     !< Number of volumes in the eta direction (real+fictitious)
      real(8), intent(in) :: dt     !< Time step (s)
      real(8), dimension (nx*ny), intent(in) :: rp     !< Radius at the center of volume P (m)
      real(8), dimension (nx*ny), intent(in) :: re     !< Radius at the center of east face of volume P (m)
      real(8), dimension (nx*ny), intent(in) :: rn     !< Radius at the center of north face of volume P (m)
      real(8), dimension (nx*ny), intent(in) :: Jp     !< Jacobian at the center of volume P (1/m2)
      real(8), dimension (nx*ny), intent(in) :: Uce    !< Contravariant velocity U at east face of volume P (m2/s)
      real(8), dimension (nx*ny), intent(in) :: Vcn    !< Contravariant velocity V at north face of volume P (m2/s)
      real(8), dimension (nx*ny), intent(in) :: roe    !< Specific mass at east face of volume P (kg/m3) (previous iteraction)
      real(8), dimension (nx*ny), intent(in) :: ron    !< Specific mass at north face of volume P (kg/m3) (previous iteraction)
      real(8), dimension (nx*ny), intent(in) :: g      !< ro/p = 1/(Rg T) for ideal gases (kg/J)
      real(8), dimension (nx*ny), intent(in) :: de     !< First SIMPLEC coefficients for Uce (m3.s/kg)
      real(8), dimension (nx*ny), intent(in) :: dn     !< First SIMPLEC coefficients for Vcn (m3.s/kg)

      real(8), dimension (nx*ny,5), intent(out) :: a   !< Coefficients of the linear system for pl (m.s)

      integer :: i, j, np, nps, npn, npw, npe
      real(8) :: as, an, aw, ae
      real(8) :: roem, rowm, ronm, rosm

      do j = 2, ny-1

         do i = 2, nx-1

            np   = nx * (j-1) + i
            nps  = np - nx
            npn  = np + nx
            npw  = np - 1
            npe  = np + 1

            as = dsign( 0.5d0, Vcn(nps) )
            an = dsign( 0.5d0, Vcn(np)  )
            aw = dsign( 0.5d0, Uce(npw) )
            ae = dsign( 0.5d0, Uce(np)  )

            roem = roe(np)   ! Density on east face
            rowm = roe(npw)  ! Density on west face
            ronm = ron(np)   ! Density on north face
            rosm = ron(nps)  ! Density on south face

            ! S
            a(np,1) = - rn(nps) * ( (0.5d0 + as) * Vcn(nps) * g(nps) + rosm * dn(nps) )

            ! W
            a(np,2) = - re(npw) * ( (0.5d0 + aw) * Uce(npw) * g(npw) + rowm * de(npw) )

            ! P
            a(np,3) = g(np) * ( rp(np)/(Jp(np)*dt)  &
               + (0.5d0 + ae) * re(np)  * Uce(np)  &
               - (0.5d0 - aw) * re(npw) * Uce(npw) &
               + (0.5d0 + an) * rn(np)  * Vcn(np)  &
               - (0.5d0 - as) * rn(nps) * Vcn(nps) ) &
               + roem * re(np)  * de(np)  &
               + rowm * re(npw) * de(npw) &
               + ronm * rn(np)  * dn(np)  &
               + rosm * rn(nps) * dn(nps)

            ! E
            a(np,4) = re(np) * ( (0.5d0 - ae) * Uce(np) * g(npe) - roem * de(np) )

            ! N
            a(np,5) = rn(np) * ( (0.5d0 - an) * Vcn(np) * g(npn) - ronm * dn(np) )

         end do

      end do

   end subroutine get_p_coefficients

   !> \brief Calculates the source of the linear system of the pressure correction
   subroutine get_p_source(nx, ny, dt, rp, re, rn, Jp, roe, ron, ro, roa &
         ,                                     Uce, Vcn, de, dn, g, pl, b) ! Output: last one
      implicit none
      integer, intent(in) :: nx     !< Number of volumes in the csi direction (real+fictitious)
      integer, intent(in) :: ny     !< Number of volumes in the eta direction (real+fictitious)
      real(8), intent(in) :: dt     !< Time step (s)
      real(8), dimension (nx*ny), intent(in) :: rp     !< Radius at the center of volume P (m)
      real(8), dimension (nx*ny), intent(in) :: re     !< Radius at the center of east face of volume P (m)
      real(8), dimension (nx*ny), intent(in) :: rn     !< Radius at the center of north face of volume P (m)
      real(8), dimension (nx*ny), intent(in) :: Jp     !< Jacobian at the center of volume P (1/m2)
      real(8), dimension (nx*ny), intent(in) :: roe    !< Specific mass at east face of volume P (kg/m3) (incorrect roe*)
      real(8), dimension (nx*ny), intent(in) :: ron    !< Specific mass at north face of volume P (kg/m3) (incorrect ron*)
      real(8), dimension (nx*ny), intent(in) :: ro     !< Specific mass (absolute density) at center of volume P (kg/m3) (incorrect ro*)
      real(8), dimension (nx*ny), intent(in) :: roa    !< Specific mass of the previous time step (kg/m3)
      real(8), dimension (nx*ny), intent(in) :: Uce    !< Contravariant velocity U at east face of volume P (m2/s) (incorrect Uce*)
      real(8), dimension (nx*ny), intent(in) :: Vcn    !< Contravariant velocity V at north face of volume P (m2/s) (incorrect Vcn*)
      real(8), dimension (nx*ny), intent(in) :: de   ! First Simplec coef. for the contravariant velocity U (east face)
      real(8), dimension (nx*ny), intent(in) :: dn   ! First Simplec coef. for the contravariant velocity V (north face)
      real(8), dimension (nx*ny), intent(in) :: g      !< g = p / ro
      real(8), dimension (nx*ny), intent(in) :: pl     !< Pressure correction

      real(8), dimension (nx*ny), intent(out) :: b     !< Source of the linear system for pl (kg/s)

      integer :: i, j, np, nps, npn, npw, npe
      real(8) :: as, an, aw, ae ! UDS coefficients
      real(8) :: roel, rowl, ronl, rosl
      real(8) :: Ucel, Ucwl, Vcnl, Vcsl

      do j = 2, ny-1

         do i = 2, nx-1

            np   = nx * (j-1) + i
            nps  = np - nx
            npn  = np + nx
            npw  = np - 1
            npe  = np + 1

            as = dsign( 0.5d0, Vcn(nps) )
            an = dsign( 0.5d0, Vcn(np)  )
            aw = dsign( 0.5d0, Uce(npw) )
            ae = dsign( 0.5d0, Uce(np)  )

            roel = g(np ) * pl(np ) * (0.5d0+ae) + g(npe) * pl(npe) * (0.5d0-ae)
            rowl = g(npw) * pl(npw) * (0.5d0+aw) + g(np ) * pl(np ) * (0.5d0-aw)
            ronl = g(np ) * pl(np ) * (0.5d0+an) + g(npn) * pl(npn) * (0.5d0-an)
            rosl = g(nps) * pl(nps) * (0.5d0+as) + g(np ) * pl(np ) * (0.5d0-as)

            Ucel = de(np ) * ( pl(np ) - pl(npe) )
            Ucwl = de(npw) * ( pl(npw) - pl(np ) )
            Vcnl = dn(np ) * ( pl(np ) - pl(npn) )
            Vcsl = dn(nps) * ( pl(nps) - pl(np ) )

            b(np) = -rp(np) * ( ro(np) - roa(np) ) / ( Jp(np) * dt ) &
               - roe(np)  * re(np)  * Uce(np)  &
               + roe(npw) * re(npw) * Uce(npw) &
               - ron(np)  * rn(np)  * Vcn(np)  &
               + ron(nps) * rn(nps) * Vcn(nps) &
               - roel * re(np ) * Ucel &
               + rowl * re(npw) * Ucwl &
               - ronl * rn(np ) * Vcnl &
               + rosl * rn(nps) * Vcsl

         end do

      end do

   end subroutine get_p_source

   subroutine get_velocities_at_internal_faces(nx, ny, dt, rp, re, rn, xe, ye, xk, yk     & ! Input
      ,                             xen, yen, xke, yke, Jp, cup, cvp, sup, svp & ! Input
      ,                             au, av, roa, p, u, v, uea, vea, una, vna   & ! Input
      ,                             ue, ve, un, vn, Uce, Vcn)                    ! Output

      implicit none
      integer, intent(in) :: nx     ! Number of volumes in csi direction (real+fictitious)
      integer, intent(in) :: ny     ! Number of volumes in eta direction (real+fictitious)
      real(8), intent(in) :: dt     ! Time step
      real(8), dimension (nx*ny),   intent(in) :: rp   ! Radius of the center of volume P
      real(8), dimension (nx*ny),   intent(in) :: re   ! Radius of the center of east face of volume P
      real(8), dimension (nx*ny),   intent(in) :: rn   ! Radius of the center of north face of volume P
      real(8), dimension (nx*ny),   intent(in) :: xe   ! x_eta at center of face east
      real(8), dimension (nx*ny),   intent(in) :: ye   ! y_eta at center of face east
      real(8), dimension (nx*ny),   intent(in) :: xk   ! x_csi at center of face north
      real(8), dimension (nx*ny),   intent(in) :: yk   ! y_csi at center of face north
      real(8), dimension (nx*ny),   intent(in) :: xen  ! x_eta at center of face north
      real(8), dimension (nx*ny),   intent(in) :: yen  ! y_eta at center of face north
      real(8), dimension (nx*ny),   intent(in) :: xke  ! x_csi at center of face east
      real(8), dimension (nx*ny),   intent(in) :: yke  ! y_csi at center of face east
      real(8), dimension (nx*ny),   intent(in) :: Jp   ! Jacobian at the center of volume P
      real(8), dimension (nx*ny),   intent(in) :: cup  ! Term of deferred correction for u
      real(8), dimension (nx*ny),   intent(in) :: cvp  ! Term of deferred correction for v
      real(8), dimension (nx*ny),   intent(in) :: sup  ! Viscous term for u
      real(8), dimension (nx*ny),   intent(in) :: svp  ! Viscous term for v
      real(8), dimension (nx*ny,9), intent(in) :: au   ! Coefficients of the linear system for u
      real(8), dimension (nx*ny,9), intent(in) :: av   ! Coefficients of the linear system for v
      real(8), dimension (nx*ny),   intent(in) :: roa  ! Absolute density at a time step before at center of vol. P
      real(8), dimension (nx*ny),   intent(in) :: p    ! Pressure at center o volume P
      real(8), dimension (nx*ny),   intent(in) :: u    ! Cartesian velocity of the present iteraction
      real(8), dimension (nx*ny),   intent(in) :: v    ! Cartesian velocity of the present iteraction
      real(8), dimension (nx*ny),   intent(in) :: uea  ! Cartesian velocity u at center of east face in the last time step
      real(8), dimension (nx*ny),   intent(in) :: vea  ! Cartesian velocity v at center of east face in the last time step
      real(8), dimension (nx*ny),   intent(in) :: una  ! Cartesian velocity u at center of north face in the last time step
      real(8), dimension (nx*ny),   intent(in) :: vna  ! Cartesian velocity v at center of north face in the last time step
      real(8), dimension (nx*ny),   intent(out) :: ue  ! Cartesian velocity u at center of east face
      real(8), dimension (nx*ny),   intent(out) :: ve  ! Cartesian velocity v at center of east face
      real(8), dimension (nx*ny),   intent(out) :: un  ! Cartesian velocity u at center of north face
      real(8), dimension (nx*ny),   intent(out) :: vn  ! Cartesian velocity v at center of north face
      real(8), dimension (nx*ny),   intent(out) :: Uce ! Contravariant velocity U at east face
      real(8), dimension (nx*ny),   intent(out) :: Vcn ! Contravariant velocity V at north face

      ! Auxiliary variables
      integer :: i, j, np, npsw, npse, nps, npw, npe, npnw, npn, npne
      integer :: npee, npsee, npnee, npnn, npnnw, npnne
      real(8) :: mpa, mea, mna, aux
      real(8) :: sumup, sumun, sumvp, sumvn, sumue, sumve
      real(8) :: pue, pve, pun, pvn

      ! Calculation of Uce at the inner faces of the real domain

      do j = 2, ny-1
         do i = 2, nx-2

            np   = nx * (j-1) + i
            nps  = np - nx
            npn  = np + nx
            npw  = np - 1
            npe  = np + 1
            npsw = nps - 1
            npse = nps + 1
            npnw = npn - 1
            npne = npn + 1

            npee  = npe  + 1
            npsee = npse + 1
            npnee = npne + 1

            mpa = roa(np ) * rp(np ) / Jp(np )
            mea = roa(npe) * rp(npe) / Jp(npe)

            sumup = au(np,1) * u(npsw) + au(np,2) * u(nps ) + au(np,3) * u(npse) + au(np,4) * u(npw ) &
               +  au(np,6) * u(npe ) + au(np,7) * u(npnw) + au(np,8) * u(npn ) + au(np,9) * u(npne)

            sumue = au(npe,1) * u(nps ) + au(npe,2) * u(npse) + au(npe,3) * u(npsee) + au(npe,4) * u(np   ) &
               +  au(npe,6) * u(npee) + au(npe,7) * u(npn ) + au(npe,8) * u(npne ) + au(npe,9) * u(npnee)

            sumvp = av(np,1) * v(npsw) + av(np,2) * v(nps ) + av(np,3) * v(npse) + av(np,4) * v(npw ) &
               +  av(np,6) * v(npe ) + av(np,7) * v(npnw) + av(np,8) * v(npn ) + av(np,9) * v(npne)

            sumve = av(npe,1) * v(nps ) + av(npe,2) * v(npse) + av(npe,3) * v(npsee) + av(npe,4) * v(np   ) &
               +  av(npe,6) * v(npee) + av(npe,7) * v(npn ) + av(npe,8) * v(npne ) + av(npe,9) * v(npnee)

            aux = ( p(npn) + p(npne) - p(nps) - p(npse) ) / 4.d0

            pue = re(np) * ( yke(np) * aux + ye(np) * ( p(np) - p(npe) ) )

            pve = re(np) * (-xke(np) * aux + xe(np) * ( p(npe) - p(np) ) )

            ue(np) = ( ( mpa + mea ) * uea(np) / dt + cup(np) + cup(npe) + sup(np) + sup(npe) &
               - sumup - sumue + 2.d0 * pue ) / ( au(np,5) + au(npe,5) )

            ve(np) = ( ( mpa + mea ) * vea(np) / dt + cvp(np) + cvp(npe) + svp(np) + svp(npe) &
               - sumvp - sumve + 2.d0 * pve ) / ( av(np,5) + av(npe,5) )

            Uce(np) = ue(np) * ye(np) - ve(np) * xe(np)

         end do
      end do

      ! Calculation of Vcn at the inner faces of the real domain

      do j = 2, ny-2
         do i = 2, nx-1

            np   = nx * (j-1) + i
            nps  = np - nx
            npn  = np + nx
            npw  = np - 1
            npe  = np + 1
            npsw = nps - 1
            npse = nps + 1
            npnw = npn - 1
            npne = npn + 1

            npnn  = npn  + nx
            npnnw = npnn - 1
            npnne = npnn + 1

            mpa = roa(np ) * rp(np ) / Jp(np )
            mna = roa(npn) * rp(npn) / Jp(npn)

            sumup = au(np,1) * u(npsw) + au(np,2) * u(nps ) + au(np,3) * u(npse) + au(np,4) * u(npw ) &
               +  au(np,6) * u(npe ) + au(np,7) * u(npnw) + au(np,8) * u(npn ) + au(np,9) * u(npne)

            sumun = au(npn,1) * u(npw ) + au(npn,2) * u(np   ) + au(npn,3) * u(npe ) + au(npn,4) * u(npnw ) &
               +  au(npn,6) * u(npne) + au(npn,7) * u(npnnw) + au(npn,8) * u(npnn) + au(npn,9) * u(npnne)

            sumvp = av(np,1) * v(npsw) + av(np,2) * v(nps ) + av(np,3) * v(npse) + av(np,4) * v(npw ) &
               +  av(np,6) * v(npe ) + av(np,7) * v(npnw) + av(np,8) * v(npn ) + av(np,9) * v(npne)

            sumvn = av(npn,1) * v(npw ) + av(npn,2) * v(np   ) + av(npn,3) * v(npe ) + av(npn,4) * v(npnw ) &
               +  av(npn,6) * v(npne) + av(npn,7) * v(npnnw) + av(npn,8) * v(npnn) + av(npn,9) * v(npnne)

            aux = ( p(npne) + p(npe) - p(npnw) - p(npw) ) / 4.d0

            pun = rn(np) * ( yk(np) * ( p(npn) - p(np) ) - yen(np) * aux )

            pvn = rn(np) * ( xk(np) * ( p(np) - p(npn) ) + xen(np) * aux )

            un(np) = ( ( mpa + mna ) * una(np) / dt + cup(np) + cup(npn) + sup(np) + sup(npn) &
               - sumup - sumun + 2.d0 * pun ) / ( au(np,5) + au(npn,5) )

            vn(np) = ( ( mpa + mna ) * vna(np) / dt + cvp(np) + cvp(npn) + svp(np) + svp(npn) &
               - sumvp - sumvn + 2.d0 * pvn ) / ( av(np,5) + av(npn,5) )

            Vcn(np) = vn(np) * xk(np) - un(np) * yk(np)

         end do
      end do

   end subroutine get_velocities_at_internal_faces

   ! Calculates SIMPLEC coefficients at the internal real faces
   subroutine get_internal_simplec_coefficients( nx, ny, re, rn, xe, ye, xk    &
      , yk, au, av, due, dve, dun, dvn, de, dn) ! Output: last 6
      implicit none
      integer, intent(in) :: nx     ! Number of volumes in csi direction (real+fictitious)
      integer, intent(in) :: ny     ! Number of volumes in eta direction (real+fictitious)
      real(8), dimension (nx*ny),   intent(in)  :: re  ! Radius of the center of east face of volume P
      real(8), dimension (nx*ny),   intent(in)  :: rn  ! Radius of the center of north face of volume P
      real(8), dimension (nx*ny),   intent(in)  :: xe  ! x_eta at center of face east
      real(8), dimension (nx*ny),   intent(in)  :: ye  ! y_eta at center of face east
      real(8), dimension (nx*ny),   intent(in)  :: xk  ! x_csi at center of face north
      real(8), dimension (nx*ny),   intent(in)  :: yk  ! y_csi at center of face north
      real(8), dimension (nx*ny,9), intent(in)  :: au  ! Coefficients of the linear system for u
      real(8), dimension (nx*ny,9), intent(in)  :: av  ! Coefficients of the linear system for v
      real(8), dimension (nx*ny),   intent(out) :: due ! Simplec coef. for the cartesian velocity u (east face)
      real(8), dimension (nx*ny),   intent(out) :: dve ! Simplec coef. for the cartesian velocity v (east face)
      real(8), dimension (nx*ny),   intent(out) :: dun ! Simplec coef. for the cartesian velocity u (north face)
      real(8), dimension (nx*ny),   intent(out) :: dvn ! Simplec coef. for the cartesian velocity v (north face)
      real(8), dimension (nx*ny),   intent(out) :: de  ! Simplec coef. for the contravariant velocity U (east face)
      real(8), dimension (nx*ny),   intent(out) :: dn  ! Simplec coef. for the contravariant velocity V (north face)


      ! Auxiliary variables
      integer :: i, j, np, npe, npn

      real(8) :: Lue, Lve, Lun, Lvn

      real(8), dimension(nx*ny) :: sumau ! Sum of the coefficients Au
      real(8), dimension(nx*ny) :: sumav ! Sum of the coefficients Av


      ! Summing the coefficients of Au and Av

      do j = 2, ny-1

         do i = 2, nx-1

            np   = nx * (j-1) + i

            sumau(np) = sum( au(np,:) )

            sumav(np) = sum( av(np,:) )

         end do

      end do

      ! Calculating the SIMPLEC coefficients over the east
      ! interface of real volumes

      do j = 2, ny-1

         do i = 2, nx-2

            np   = nx * (j-1) + i
            npe  = np + 1

            Lue = 2.d0 * re(np) / ( sumau(np) + sumau(npe) )

            Lve = 2.d0 * re(np) / ( sumav(np) + sumav(npe) )

            due(np) = Lue * ye(np)

            dve(np) = Lve * xe(np)

            de(np) = ye(np) * due(np) + xe(np) * dve(np)

         end do

      end do

      ! Calculating the SIMPLEC coefficients over the north
      ! interface of real volumes

      do j = 2, ny-2

         do i = 2, nx-1

            np   = nx * (j-1) + i
            npn  = np + nx

            Lun = 2.d0 * rn(np) / ( sumau(np) + sumau(npn) )

            Lvn = 2.d0 * rn(np) / ( sumav(np) + sumav(npn) )

            ! SIMPLEC coefficients

            dun(np) = Lun * yk(np)

            dvn(np) = Lvn * xk(np)

            dn(np) = xk(np) * dvn(np) + yk(np) * dun(np)

         end do

      end do

   end subroutine get_internal_simplec_coefficients

   subroutine get_pressure_correction_with_pl( nx, ny, pl, p) ! InOutput: last one
      implicit none
      integer, intent(in) :: nx     ! Number of volumes in csi direction (real+fictitious)
      integer, intent(in) :: ny     ! Number of volumes in eta direction (real+fictitious)
      real(8), dimension (nx*ny), intent(in)    :: pl ! Pressure correction
      real(8), dimension (nx*ny), intent(inout) :: p  ! p of present iteraction is input. Corrected p is output.

      p = p + pl

   end subroutine get_pressure_correction_with_pl


   subroutine get_density_correction_with_pl( nx, ny, pl, g, ro) ! InOutput: last one
      implicit none
      integer, intent(in) :: nx     ! Number of volumes in csi direction (real+fictitious)
      integer, intent(in) :: ny     ! Number of volumes in eta direction (real+fictitious)
      real(8), dimension (nx*ny), intent(in)    :: pl ! Pressure correction
      real(8), dimension (nx*ny), intent(in)    :: g  ! ro/p = 1/(Rg T) for ideal gases
      real(8), dimension (nx*ny), intent(inout) :: ro ! ro of present iteraction is input. Corrected ro is output.

      ro = ro + pl * g

   end subroutine get_density_correction_with_pl

   subroutine get_u_v_at_real_nodes_with_pl( nx, ny, xe, ye, xk & ! Input
      ,                                    yk, rp, pl, au, av & ! Input
      ,                                    u, v )               ! Input and Output
      implicit none
      integer, intent(in) :: nx     ! Number of volumes in csi direction (real+fictitious)
      integer, intent(in) :: ny     ! Number of volumes in eta direction (real+fictitious)
      real(8), dimension (nx*ny),   intent(in)    :: xe  ! x_eta at center of face east
      real(8), dimension (nx*ny),   intent(in)    :: ye  ! y_eta at center of face east
      real(8), dimension (nx*ny),   intent(in)    :: xk  ! x_csi at center of face north
      real(8), dimension (nx*ny),   intent(in)    :: yk  ! y_csi at center of face north
      real(8), dimension (nx*ny),   intent(in)    :: rp  ! Radius of the center of volume P
      real(8), dimension (nx*ny),   intent(in)    :: pl  ! Pressure correction
      real(8), dimension (nx*ny,9), intent(in)    :: au  ! Coefficients of the linear system for u
      real(8), dimension (nx*ny,9), intent(in)    :: av  ! Coefficients of the linear system for v
      real(8), dimension (nx*ny),   intent(inout) :: u   ! u of present iteraction is input. u corrected is output
      real(8), dimension (nx*ny),   intent(inout) :: v   ! v of present iteraction is input. v corrected is output

      ! Auxiliary variables
      integer :: i, j, np, nps, npn, npw, npe
      real(8) :: plup, plvp

      do j = 2, ny-1
         do i = 2, nx-1
            np   = nx * (j-1) + i
            nps  = np - nx
            npn  = np + nx
            npw  = np - 1
            npe  = np + 1

            plup  =   0.5d0 *  rp(np)  * ( &
               +  yk(np ) * ( pl(npn) + pl(np) ) &
               -  yk(nps) * ( pl(nps) + pl(np) ) &
               -  ye(np ) * ( pl(npe) + pl(np) ) &
               +  ye(npw) * ( pl(npw) + pl(np) ) &
               )

            plvp  =   0.5d0 *   rp(np) * ( &
               +  xe(np ) * ( pl(np) + pl(npe) ) &
               -  xe(npw) * ( pl(np) + pl(npw) ) &
               -  xk(np ) * ( pl(np) + pl(npn) ) &
               +  xk(nps) * ( pl(np) + pl(nps) ) &
               )

            u(np) = u(np) + plup / sum( au(np,:) )
            v(np) = v(np) + plvp / sum( av(np,:) )

         end do
      end do

   end subroutine get_u_v_at_real_nodes_with_pl



   subroutine get_velocity_correction_at_faces_with_pl( nx, ny, xe, ye, xk, yk &
      , due, dve, dun, dvn, de, dn, pl, ue, ve, un, vn, Uce, Vcn ) ! InOutput: last six
      implicit none
      integer, intent(in) :: nx     ! Number of volumes in csi direction (real+fictitious)
      integer, intent(in) :: ny     ! Number of volumes in eta direction (real+fictitious)
      real(8), dimension (nx*ny), intent(in)    :: xe   ! x_eta at center of face east
      real(8), dimension (nx*ny), intent(in)    :: ye   ! y_eta at center of face east
      real(8), dimension (nx*ny), intent(in)    :: xk   ! x_csi at center of face north
      real(8), dimension (nx*ny), intent(in)    :: yk   ! y_csi at center of face north
      real(8), dimension (nx*ny), intent(in)    :: due  ! Simplec coef. for the cartesian velocity u (east face)
      real(8), dimension (nx*ny), intent(in)    :: dve  ! Simplec coef. for the cartesian velocity v (east face)
      real(8), dimension (nx*ny), intent(in)    :: dun  ! Simplec coef. for the cartesian velocity u (north face)
      real(8), dimension (nx*ny), intent(in)    :: dvn  ! Simplec coef. for the cartesian velocity v (north face)
      real(8), dimension (nx*ny), intent(in)    :: de   ! SIMPLEC coefficient for Uce
      real(8), dimension (nx*ny), intent(in)    :: dn   ! SIMPLEC coefficient for Vcn
      real(8), dimension (nx*ny), intent(in)    :: pl   ! Pressure correction
      real(8), dimension (nx*ny), intent(inout) :: ue   ! Cartesian velocity u at center of east face
      real(8), dimension (nx*ny), intent(inout) :: ve   ! Cartesian velocity v at center of east face
      real(8), dimension (nx*ny), intent(inout) :: un   ! Cartesian velocity u at center of north face
      real(8), dimension (nx*ny), intent(inout) :: vn   ! Cartesian velocity v at center of north face
      real(8), dimension (nx*ny), intent(inout) :: Uce  ! Contravariant velocity U at east face
      real(8), dimension (nx*ny), intent(inout) :: Vcn  ! Contravariant velocity V at north face

      ! Auxiliary variables
      integer :: i, j, np, npn, nps, npe, npw, npse, npne, npnw


      do j = 2, ny-1

         do i = 1, nx-1

            np   = nx * (j-1) + i
            nps  = np - nx
            npn  = np + nx
            npe  = np + 1
            npse = nps + 1
            npne = npn + 1


            Uce(np) = ue(np) * ye(np) - ve(np) * xe(np) &

               + de(np) * ( pl(np) - pl(npe) )


            ue(np) = ue(np) + due(np) * ( pl(np) - pl(npe) )


            ve(np) = ve(np) + dve(np) * ( pl(npe) - pl(np) )

         end do

      end do

      do j = 1, ny-1

         do i = 2, nx-1

            np   = nx * (j-1) + i
            npn  = np + nx
            npw  = np - 1
            npe  = np + 1
            npnw = npn - 1
            npne = npn + 1


            Vcn(np) = vn(np) * xk(np) - un(np) * yk(np) &

               + dn(np) * ( pl(np) - pl(npn) )


            un(np) = un(np) + dun(np) * ( pl(npn) - pl(np) )


            vn(np) = vn(np) + dvn(np) * ( pl(np) - pl(npn) )

         end do

      end do

   end subroutine get_velocity_correction_at_faces_with_pl



   !> \brief Calculates the temperature based on the conservation of the total enthalpy
   !! Valid for Euler model with constant thermo-physical coefficients.
   !! Temperature is extrapolated to fictitious volumes with CDS the scheme
   subroutine get_T_from_H_conservation(nx, ny, CPF, HF, u, ue, un, v, ve, vn, T, Tbe, Tbw, Tbn, Tbs) ! Output: last 5
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
      real(8), dimension (ny),    intent(out) :: Tbw  !< Temperature over the  west boundary (K)
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

         Tbw(j) = ( HF - (ue(np)**2+ve(np)**2) / 2.d0 ) / CPF

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


   end subroutine get_T_from_H_conservation


   !> \brief Initialize the flags of boundary faces
   subroutine initialize_boundary_faces_flags(nx, ny, fbe, fbn) ! Output: last 2
      implicit none
      integer, intent(in) :: nx  !< Number of volumes in csi direction (real+fictitious)
      integer, intent(in) :: ny  !< Number of volumes in eta direction (real+fictitious)
      real(8), dimension (nx*ny), intent(out) :: fbe !< Face of boundary east (1 if an east boundary, 0 otherwise)
      real(8), dimension (nx*ny), intent(out) :: fbn !< Face of boundary north (1 if an north boundary, 0 otherwise)

      ! Inner variables
      integer :: i, j, np

      fbe = 0.d0

      ! East boundary

      i = 1

      do j = 2, ny-1

         np   = nx * (j-1) + i

         fbe(np) = 1.d0

      end do

      ! West boundary

      i = nx-1

      do j = 2, ny-1

         np   = nx * (j-1) + i

         fbe(np) = 1.d0

      end do



      fbn = 0.d0

      ! North boundary

      j = ny-1

      do i = 2, nx-1

         np   = nx * (j-1) + i

         fbn(np) = 1.d0

      end do

      ! South boundary

      j = 1

      do i = 2, nx-1

         np   = nx * (j-1) + i

         fbn(np) = 1.d0

      end do

   end subroutine initialize_boundary_faces_flags

end module coefficients
