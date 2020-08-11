!>
!! \brief mod_grid_procedures contains procedures to operate on grid data
!!
module mod_grid_procedures

   use thompson2d_hyperbolic, only: get_thompson2d_hyperbolic_01

   implicit none

contains

   subroutine set_grid(kg, nx, ny, a1, avi, avf, awf, x, y) ! Last 2 are inoutput
      implicit none
      integer, intent(in)  :: kg ! Kind of grid (1=uniform, 2=geometric progression, 3=power law, 4=gp modified)
      integer, intent(in)  :: nx
      integer, intent(in)  :: ny
      real(8), intent(in)  :: a1
      real(8), intent(in)  :: avi      !< Initial value of the artificial viscosity
      real(8), intent(in)  :: avf      !< Final value of the artificial viscosity
      real(8), intent(in)  :: awf      !< Area weighting factor
      real(8), intent(inout) :: x(nx*ny)
      real(8), intent(inout) :: y(nx*ny)
      !

      if ( kg == 1 ) then

         call get_uniform_grid( nx, ny, x, y) ! Last 2 are inoutput

      else if ( kg == 2 ) then

         call get_pg_grid(a1, nx, ny, x, y) ! Last 2 are inoutput

      else if ( kg == 3 ) then

         call get_power_grid(a1, nx, ny, x, y) ! Last 2 are inoutput

      else if ( kg == 4 ) then

         call get_pg_grid2(a1, nx, ny, x, y) ! Last 2 are inoutput

      else if ( kg == 5 ) then

         call get_thompson2d_hyperbolic_01(nx, ny, avi, avf, awf, a1, x, y) ! Inoutput: last two

      else

         write(*,*) 'set_grid: Type of grid unknown'
         stop

      end if

   end subroutine set_grid

   subroutine get_coarser_grid(nxf, nyf, nx, ny, nmf, nmd, xf, yf, x, y)
      implicit none
      integer, intent(in)  :: nxf
      integer, intent(in)  :: nyf
      integer, intent(in)  :: nx
      integer, intent(in)  :: ny
      integer, intent(in)  :: nmf
      integer, intent(in)  :: nmd
      real(8), intent(in)  :: xf(nxf*nyf)
      real(8), intent(in)  :: yf(nxf*nyf)
      real(8), intent(out) :: x(nx*ny)
      real(8), intent(out) :: y(nx*ny)


      integer :: i, j, k, ii, jj, np, npf


      k = nmf-nmd

      do j = 1, ny-1

         jj = 2 ** k * (j-1) + 1

         do i = 1, nx-1

            ii = 2 ** k * (i-1) + 1

            np  = nx  * (j -1) + i

            npf = nxf * (jj-1) + ii

            x(np) = xf(npf)

            y(np) = yf(npf)

         end do

      end do


   end subroutine


   subroutine get_uniform_grid( nx, ny, x, y) ! Last 2 are inoutput
      implicit none
      integer, intent(in)  :: nx
      integer, intent(in)  :: ny
      real(8), intent(inout) :: x(nx*ny)
      real(8), intent(inout) :: y(nx*ny)
      !
      integer :: i, j, np, nps
      real(8) :: xi, xf, yi, yf, dx, dy

      do i = 1, nx-1

         j = 1
         np = nx * (j-1) + i
         xi = x(np)
         yi = y(np)

         j = ny-1
         np = nx * (j-1) + i
         xf = x(np)
         yf = y(np)

         dx = (xf-xi) / (ny-2)
         dy = (yf-yi) / (ny-2)

         do j = 2, ny-2
            np  = nx * (j-1) + i
            nps = np - nx
            x(np) = x(nps) + dx
            y(np) = y(nps) + dy
         end do
      end do

   end subroutine get_uniform_grid


   subroutine get_power_grid(a1, nx, ny, x, y) ! Last 2 are inoutput
      implicit none
      integer, intent(in)  :: nx
      integer, intent(in)  :: ny
      real(8), intent(in)  :: a1
      real(8), intent(inout) :: x(nx*ny)
      real(8), intent(inout) :: y(nx*ny)
      !
      integer :: i, j, np
      real(8) :: a, r, xi, xf, yi, yf

      do i = 1, nx-1

         j = 1
         np = nx * (j-1) + i
         xi = x(np)
         yi = y(np)

         j = ny-1
         np = nx * (j-1) + i
         xf = x(np)
         yf = y(np)

         r = sqrt( (xf-xi)**2 + (yf-yi)**2 ) / a1

         if ( r > 1.d0 ) then

            a = log( r ) / log(dble(ny-2))

         else

            a = 1.d0

         end if

         do j = 2, ny-2

            np  = nx * (j-1) + i

            x(np) = (xf-xi) * (dble(j-1)/dble(ny-2))**a + xi

            y(np) = (yf-yi) * (dble(j-1)/dble(ny-2))**a + yi

         end do

      end do

   end subroutine get_power_grid


   subroutine get_pg_grid(a1, nx, ny, x, y) ! Last 2 are inoutput
      implicit none
      integer, intent(in)  :: nx
      integer, intent(in)  :: ny
      real(8), intent(in)  :: a1
      real(8), intent(inout) :: x(nx*ny)
      real(8), intent(inout) :: y(nx*ny)
      !
      integer :: i, j, np, npn
      real(8) :: q, r, t, xi, xf, yi, yf

      do i = 1, nx-1

         j = 1
         np = nx * (j-1) + i
         xi = x(np)
         yi = y(np)

         j = ny-1
         np = nx * (j-1) + i
         xf = x(np)
         yf = y(np)

         r = sqrt( (xf-xi)**2 + (yf-yi)**2 ) / a1

         call get_GP_ratio(ny-2,r,q) ! Last one is output

         t = 0.d0

         do j = 1, ny-3

            np  = nx * (j-1) + i

            npn  = np + nx

            t = t + q**(j-1) / r ! t of j+1

            x(npn) = xi + (xf-xi) * t

            y(npn) = yi + (yf-yi) * t

         end do

      end do

   end subroutine get_pg_grid

   subroutine get_pg_grid2(a1, nx, ny, x, y) ! Last 2 are inoutput
      implicit none
      integer, intent(in)  :: nx
      integer, intent(in)  :: ny
      real(8), intent(in)  :: a1
      real(8), intent(inout) :: x(nx*ny)
      real(8), intent(inout) :: y(nx*ny)
      !
      integer :: i, j, np, npn
      real(8) :: q, r, t, xi, xf, yi, yf, lf, l, a

      i = nx-1
      j = 1

      np = nx * (j-1) + i
      xi = x(np)
      yi = y(np)

      j = ny-1
      np = nx * (j-1) + i
      xf = x(np)
      yf = y(np)

      lf = sqrt( (xf-xi)**2 + (yf-yi)**2 )


      do i = 1, nx-1

         j = 1
         np = nx * (j-1) + i
         xi = x(np)
         yi = y(np)

         j = ny-1
         np = nx * (j-1) + i
         xf = x(np)
         yf = y(np)

         l = sqrt( (xf-xi)**2 + (yf-yi)**2 )

         a = a1 * l / lf

         r = sqrt( (xf-xi)**2 + (yf-yi)**2 ) / a

         call get_GP_ratio(ny-2,r,q) ! Last one is output

         t = 0.d0

         do j = 1, ny-3

            np  = nx * (j-1) + i

            npn  = np + nx

            t = t + q**(j-1) / r ! t of j+1

            x(npn) = xi + (xf-xi) * t

            y(npn) = yi + (yf-yi) * t

         end do

      end do

   end subroutine get_pg_grid2


   subroutine get_gp_ratio(n, r, q)
     implicit none
     integer, intent(in)  ::   n !< number of partitions
     real(8), intent(in)  ::   r !< l/a1
     real(8), intent(out) ::   q !< q

     ! Parameters

     integer :: nit = 1000   ! Maximum number of iteractions
     real(8) :: tol = 1.d-15 ! Tolerance

     ! Inner variables

     integer ::   i ! Dummy index
     real(8) ::  qi ! inital value of q
     real(8) ::  qf ! final value of q
     real(8) ::  qm ! mean value of q

     if ( r < n ) then

        qi = 0.1d0

        qf = 1.d0 - 1.d-15

     else

        qi = 1.d0 + 1.d-15

        qf = 10.d0

     end if

     do i = 1, nit

        qm = 0.5d0 * qi + 0.5d0 * qf

        if ( 0.d0 < f(qi) * f(qm) ) then

           qi = qm

        else

           qf = qm

        end if

        if ( abs(qf-qi) < tol ) exit

     end do


     if ( i == nit ) then

        write(*,*) "get_gp_ratio: Maximum number of iteractions was exceeded."

        stop

     end if

     q = qm

   contains

     real(8) function f(q)
       implicit none
       real(8), intent(in) :: q

       f = q ** n + r * ( 1.d0 - q ) - 1.d0

     end function f

   end subroutine get_gp_ratio

   subroutine get_real_centroids_xy( opt, nx, ny, x, y & ! Input
      ,                            xp, yp )            ! Output
      implicit none
      integer, intent(in) :: opt ! (1=simple mean, 2=weighted mean)
      integer, intent(in) :: nx  ! Number of volumes in csi direction (real+fictitious)
      integer, intent(in) :: ny  ! Number of volumes in eta direction (real+fictitious)
      real(8), intent(in) :: x(nx*ny) ! Coord. x of the northest corner of volume P
      real(8), intent(in) :: y(nx*ny) ! Coord. y of the northest corner of volume P
      real(8), intent(out) :: xp(nx*ny) ! Coord. x of the centroid of volume P
      real(8), intent(out) :: yp(nx*ny) ! Coord. y of the centroid of volume P
      !
      integer :: i, j, np, nps, npw, npsw
      real(8) :: A1, A2

      ! Centroids of real volumes

      xp = 0.d0
      yp = 0.d0

      if ( opt == 1 ) then ! Simple mean
         ! Centroids are given by the mean of the corners
         do j = 2, ny-1
            do i = 2, nx-1

               np   = nx * (j-1) + i
               nps  = np - nx
               npw  = np - 1
               npsw = nps - 1

               xp(np) = ( x(np) + x(npw) + x(npsw) + x(nps) ) / 4.d0
               yp(np) = ( y(np) + y(npw) + y(npsw) + y(nps) ) / 4.d0

            end do
         end do
      else ! weighted mean
         ! The rectangle is divided in two triangles of area A1 and A2.
         ! The centroids of each triangle is calculated by the mean of the corners.
         ! The centroid of the rectangle is given by the weighted mean of the centroids of each triangle.
         do j = 2, ny-1
            do i = 2, nx-1

               np   = nx * (j-1) + i
               nps  = np - nx
               npw  = np - 1
               npsw = nps - 1

               A1 = 0.5d0 * ( x(npsw) * (y(nps)-y(np) ) &
                  +          x(nps) * (y(np) -y(npsw)) &
                  +          x(np)  * (y(npsw)-y(nps)))

               A2 = 0.5d0 * ( x(npsw) * (y(np) -y(npw)) &
                  +          x(np)  * (y(npw)-y(npsw)) &
                  +          x(npw) * (y(npsw)-y(np) ))

               xp(np) = ( (x(npsw)+x(nps)+x(np)) * A1 + (x(npsw)+x(np)+x(npw)) * A2 ) &
                  / ( 3.d0 * (A1+A2) )

               yp(np) = ( (y(npsw)+y(nps)+y(np)) * A1 + (y(npsw)+y(np)+y(npw)) * A2 ) &
                  / ( 3.d0 * (A1+A2) )

            end do
         end do
      end if
   end subroutine get_real_centroids_xy


   subroutine get_metrics(nx, ny, x, y, xp, yp            & ! Input
      ,          xe, ye, xen, yen, xk, yk, xke, yke, Jp & ! Output
      ,          Je, Jn, alphae, gamman, betae, betan )   ! Output
      implicit none
      integer, intent(in)  :: nx, ny  ! Number of volumes in csi and eta directions (real+fictitious)
      real(8), dimension (nx*ny), intent(in)  :: x,  y  ! Coord. of the northest corner of volume P
      real(8), dimension (nx*ny), intent(in)  :: xp, yp ! Coord. of the centroid of volume P
      real(8), dimension (nx*ny), intent(out) :: xe     ! x_eta at face east of volume P
      real(8), dimension (nx*ny), intent(out) :: ye     ! y_eta at face east of volume P
      real(8), dimension (nx*ny), intent(out) :: xen    ! x_eta at face north of volume P
      real(8), dimension (nx*ny), intent(out) :: yen    ! y_eta at face north of volume P
      real(8), dimension (nx*ny), intent(out) :: xk     ! x_csi at face north of volume P
      real(8), dimension (nx*ny), intent(out) :: yk     ! y_csi at face north of volume P
      real(8), dimension (nx*ny), intent(out) :: xke    ! x_csi at face east of volume P
      real(8), dimension (nx*ny), intent(out) :: yke    ! y_csi at face east of volume P
      real(8), dimension (nx*ny), intent(out) :: Jp     ! Jacobian at the center of volume P
      real(8), dimension (nx*ny), intent(out) :: Je     ! Jacobian at the center of east face of volume P
      real(8), dimension (nx*ny), intent(out) :: Jn     ! Jacobian at the center of north face of volume P
      real(8), dimension (nx*ny), intent(out) :: Alphae ! (metric) Alpha at the center of east face of volume P
      real(8), dimension (nx*ny), intent(out) :: Gamman ! (metric) Gamma at the center of north face of volume P
      real(8), dimension (nx*ny), intent(out) :: Betae  ! (metric) Beta  at the center of east face of volume P
      real(8), dimension (nx*ny), intent(out) :: Betan  ! (metric) Beta  at the center of north face of volume P
      !
      integer :: i, j, npsw, nps, npw, np, npe, npn
      real(8) :: fw, fe, fs, fn, xkp, xep, ykp, yep
      !
      xe  = 0.d0
      ye  = 0.d0
      xen = 0.d0
      yen = 0.d0
      xk  = 0.d0
      yk  = 0.d0
      xke = 0.d0
      yke = 0.d0
      Jp  = 0.d0
      Je  = 0.d0
      Jn  = 0.d0
      alphae = 0.d0
      gamman = 0.d0
      betae  = 0.d0
      betan  = 0.d0

      ! Derivatives relatively to eta at the center of the east face

      do j = 2, ny-1
         do i = 1, nx-1

            np  = nx * (j-1) + i
            nps = np - nx

            xe(np) = x(np) - x(nps)
            ye(np) = y(np) - y(nps)

            alphae(np) = xe(np)**2 + ye(np)**2

         end do
      end do

      ! Derivatives relatively to csi at the center of the north face

      do j = 1, ny-1
         do i = 2, nx-1

            np  = nx * (j-1) + i
            npw = np - 1

            xk(np) = x(np) - x(npw)
            yk(np) = y(np) - y(npw)

            gamman(np) = xk(np)**2 + yk(np)**2

         end do
      end do

      ! Derivatives relatively to csi at the center of the east face (only inner real faces)

      do j = 2, ny-1
         do i = 2, nx-2

            np   = nx * (j-1) + i
            npe  = np + 1

            xke(np) = xp(npe) - xp(np)
            yke(np) = yp(npe) - yp(np)

         end do
      end do

      ! Derivatives relatively to csi at the center of the east face (only west boundary of the domain)

      i = 2
      do j = 2, ny-1

         np   = nx * (j-1) + i
         nps  = np - nx
         npw  = np - 1
         npsw = nps - 1

         fw = ( x(npw) + x(npsw) ) / 2.d0
         fe = ( x(np)  + x(nps)  ) / 2.d0

         xke(npw) = -3.d0 * fw + 4.d0 * xp(np) -fe

         fw = ( y(npw) + y(npsw) ) / 2.d0
         fe = ( y(np)  + y(nps)  ) / 2.d0

         yke(npw) = -3.d0 * fw + 4.d0 * yp(np) -fe

      end do

      ! Derivatives relatively to csi at the center of the east face (only east boundary of the domain)

      i = nx-1
      do j = 2, ny-1

         np   = nx * (j-1) + i
         nps  = np - nx
         npw  = np - 1
         npsw = nps - 1

         fw = ( x(npw) + x(npsw) ) / 2.d0
         fe = ( x(np)  + x(nps)  ) / 2.d0

         xke(np) = fw - 4.d0 * xp(np) + 3.d0 * fe

         fw = ( y(npw) + y(npsw) ) / 2.d0
         fe = ( y(np)  + y(nps)  ) / 2.d0

         yke(np) = fw - 4.d0 * yp(np) + 3.d0 * fe

      end do

      ! Derivatives relatively to eta at the center of the north face (only inner real faces)

      do j = 2, ny-2
         do i = 2, nx-1

            np  = nx * (j-1) + i
            npn = np + nx

            xen(np) = xp(npn) - xp(np)
            yen(np) = yp(npn) - yp(np)

         end do
      end do

      ! Derivatives relatively to eta at the center of the north face (only south boundary of the domain)

      j = 2

      do i = 2, nx-1

         np   = nx * (j-1) + i
         nps  = np - nx
         npw  = np - 1
         npsw = nps - 1

         fs = ( x(nps) + x(npsw) ) / 2.d0
         fn = ( x(np)  + x(npw)  ) / 2.d0

         xen(nps) = -3.d0 * fs + 4.d0 * xp(np) -fn

         fs = ( y(nps) + y(npsw) ) / 2.d0
         fn = ( y(np)  + y(npw)  ) / 2.d0

         yen(nps) = -3.d0 * fs + 4.d0 * yp(np) -fn

      end do

      ! Derivatives relatively to eta at the center of the north face (only north boundary of the domain)

      j = ny - 1
      do i = 2, nx-1

         np   = nx * (j-1) + i
         nps  = np - nx
         npw  = np - 1
         npsw = nps - 1

         fs = ( x(nps) + x(npsw) ) / 2.d0
         fn = ( x(np)  + x(npw)  ) / 2.d0

         xen(np) = fs - 4.d0 * xp(np) + 3.d0 * fn

         fs = ( y(nps) + y(npsw) ) / 2.d0
         fn = ( y(np)  + y(npw)  ) / 2.d0

         yen(np) = fs - 4.d0 * yp(np) + 3.d0 * fn

      end do

      ! Beta and J at the center of the east face (all real faces)

      do j = 2, ny-1
         do i = 1, nx-1

            np   = nx * (j-1) + i

            betae(np) = xke(np) * xe(np) + yke(np) * ye(np)

            Je(np) = 1.d0 / ( xke(np) * ye(np) - xe(np) * yke(np) )

         end do
      end do

      ! Beta and J at center of the north face (all real faces)

      do j = 1, ny-1
         do i = 2, nx-1

            np   = nx * (j-1) + i

            betan(np) = xk(np) * xen(np) + yk(np) * yen(np)

            Jn(np) = 1.d0 / ( xk(np) * yen(np) - xen(np) * yk(np) )

         end do
      end do

      ! Jacobian J at the center of all real volumes

      do j = 2, ny-1
         do i = 2, nx-1

            np   = nx * (j-1) + i
            nps  = np - nx
            npw  = np - 1
            npsw = nps - 1

            fw = ( x(npw) + x(npsw) ) / 2.d0
            fe = ( x(np)  + x(nps)  ) / 2.d0
            fs = ( x(nps) + x(npsw) ) / 2.d0
            fn = ( x(np)  + x(npw)  ) / 2.d0

            xkp = fe - fw ! (x_csi)_P
            xep = fn - fs ! (x_eta)_P

            fw = ( y(npw) + y(npsw) ) / 2.d0
            fe = ( y(np)  + y(nps)  ) / 2.d0
            fs = ( y(nps) + y(npsw) ) / 2.d0
            fn = ( y(np)  + y(npw)  ) / 2.d0

            ykp = fe - fw ! (y_csi)_P
            yep = fn - fs ! (y_eta)_P

            Jp(np) = 1.d0 / ( xkp * yep - xep * ykp )

         end do
      end do

      Je = dabs(Je)
      Jn = dabs(Jn)
      Jp = dabs(Jp)

   end subroutine get_metrics


   subroutine get_radius(coord, nx, ny, y, yp & ! Input
      ,                radius, re, rn, rp   ) ! Output
      implicit none
      integer, intent(in)  :: coord   ! Coordinate system ( 1=cylindrical, 0 = cartesian)
      integer, intent(in)  :: nx, ny  ! Number of volumes in csi and eta directions (real+fictitious)
      real(8), dimension (nx*ny), intent(in)  :: y  ! Coord. y of the northeast corner of volume P
      real(8), dimension (nx*ny), intent(in)  :: yp ! Coord. y of the centroid of volume P
      real(8), dimension (nx*ny), intent(out) :: radius ! Radius of northest corner of volume P
      real(8), dimension (nx*ny), intent(out) :: re ! Radius of the center of east face of volume P
      real(8), dimension (nx*ny), intent(out) :: rn ! Radius of the center of north face of volume P
      real(8), dimension (nx*ny), intent(out) :: rp ! Radius of the center of volume P

      integer :: i, j, np, npe, npw, npn, nps, npsw, npse, npnw, npne

      ! Cartesian coordinate system

      radius = 1.d0
      re = 1.d0
      rn = 1.d0
      rp = 1.d0

      ! Cylindrical coordinate system

      if ( coord == 1 ) then

         ! Radius of northeast corner

         radius = y

         ! Radius at the center of the east face

         do i = 1, nx-1
            do j = 2, ny-1

               np  = nx * (j-1) + i
               nps = np - nx

               re(np) = ( y(np) + y(nps) ) / 2.d0

            end do
         end do

         ! Radius at the center of the north face

         do i = 2, nx-1
            do j = 1, ny-1

               np  = nx * (j-1) + i
               npw = np - 1

               rn(np) = ( y(np) + y(npw) ) / 2.d0

            end do
         end do

         ! Radius of the center of the volume

         rp = yp

         ! Radius  of the center of the volume (south fictitious)

         j = 1
         do i = 2, nx-1

            np  = nx * (j-1) + i
            npn = np + nx

            rp(np) = rn(np) - rp(npn) + rn(np)
            !rp(np) = rp(npn) - rp(npn+nx) + rp(npn)

         end do

         ! Radius  of the center of the volume (north fictitious)

         j = ny
         do i = 2, nx-1

            np  = nx * (j-1) + i
            nps = np - nx

            rp(np) = rn(nps) - rp(nps) + rn(nps)
            !rp(np) = rp(nps) - rp(nps-nx) + rp(nps)

         end do

         ! Radius  of the center of the volume (west fictitious)

         i = 1
         do j = 2, ny-1

            np  = nx * (j-1) + i
            npe = np + 1

            rp(np) = re(np) - rp(npe) + re(np)
            !rp(np) = rp(npe) - rp(npe+1) + rp(npe)

         end do

         ! Radius  of the center of the volume (east fictitious)

         i = nx
         do j = 2, ny-1

            np  = nx * (j-1) + i
            npw = np - 1

            rp(np) = re(npw) - rp(npw) + re(npw)
            !rp(np) = rp(npw) - rp(npw-1) + rp(npw)

         end do

         ! Radius  of the center of the volume (corner fictitious: SW)

         i = 1 ; j = 1

         np  = nx * (j-1) + i
         npn = np + nx
         npe  = np + 1
         npne = npn + 1

         !rp(np) = rp(npn) - rp(npn+nx) + rp(npn)
         rp(np) = 4.d0 * y(np) - rp(npe) - rp(npn) - rp(npne)

         ! Radius  of the center of the volume (corner fictitious: SE)

         i = nx; j = 1

         np  = nx * (j-1) + i
         npn = np + nx
         npw  = np - 1
         npnw = npn - 1

         !rp(np) = rp(npn) - rp(npn+nx) + rp(npn)
         rp(np) = 4.d0 * y(npw) - rp(npw) - rp(npn) - rp(npnw)

         ! Radius  of the center of the volume (corner fictitious: NW)

         i = 1; j = ny

         np  = nx * (j-1) + i
         nps = np - nx
         npe  = np + 1
         npse = nps + 1

         !rp(np) = rp(nps) - rp(nps-nx) + rp(nps)
         rp(np) = 4.d0 * y(nps) - rp(npe) - rp(nps) - rp(npse)

         ! Radius  of the center of the volume (corner fictitious: NE)

         i = nx; j = ny

         np  = nx * (j-1) + i
         nps = np - nx
         npw  = np - 1
         npsw = nps - 1

         !rp(np) = rp(nps) - rp(nps-nx) + rp(nps)
         rp(np) = 4.d0 * y(npsw) - rp(npw) - rp(nps) - rp(npsw)

      end if

   end subroutine get_radius


end module
