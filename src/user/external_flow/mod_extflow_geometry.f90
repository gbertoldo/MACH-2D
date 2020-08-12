!>
!! \brief mod_extflow_geometry provides user-defined geometries of the north
!!        and south boundaries of the domain.
!!
module mod_extflow_geometry

   use depp_interface

   use newton2coef, only: get_newton2coef_body

   use geometry_splwsp

   use geometry_splwep

   use geometry_ncspline

   use geometry_d2cspline

   use geometry_spline_s01

   implicit none

contains

   !> \brief Generates the grid north and south boundary
   subroutine get_grid_boundary(nx, ny, unt, fgeom, lr, rb, x, y) ! Output: last four
      implicit none
      integer, intent(in) :: nx  !< Number of volumes in the csi direction (real+fictitious)
      integer, intent(in) :: ny  !< Number of volumes in the eta direction (real+fictitious)
      integer, intent(in) :: unt !< Unit where the input parameters will be printed
      character(len=*), intent(in) :: fgeom !< File of the geometric parameters

      real(8), intent(out) :: lr  !< Length of the body (m)
      real(8), intent(out) :: rb  !< Base radius of the body (m)
      real(8), dimension(nx*ny), intent(out) :: x !< Coord. x at the northeast corner of the volume P (m)
      real(8), dimension(nx*ny), intent(out) :: y !< Coord. y at the northeast corner of the volume P (m)


      ! Inner variables

      integer :: kgb ! Kind of geometry of the body



      ! Reading the kind of geometry of the body

      open(13,file=fgeom)

      read(13,*) kgb  ! Kind of geometry of the body

      close(13)



      ! Selecting the kind of geometry

      select case (kgb)

         case (1) ! Power law

            call get_grid_boundary_g01(nx, ny, unt, fgeom, lr, rb, x, y) ! Output: last two

         case (2) ! Hemisphere-cone

            call get_grid_boundary_g02(nx, ny, unt, fgeom, lr, rb, x, y) ! Output: last two

         case (3) ! Hemisphere-cone-cylinder

            call get_grid_boundary_g03(nx, ny, unt, fgeom, lr, rb, x, y) ! Output: last two

         case (4) ! Power law with DEPP interface for optimization of the exponent

            call get_grid_boundary_g04(nx, ny, unt, fgeom, lr, rb, x, y) ! Output: last two

         case (5) ! Power law-cylinder

            call get_grid_boundary_g05(nx, ny, unt, fgeom, lr, rb, x, y) ! Output: last two

         case (6) ! Power law-cylinder with a quadratic and power law distribution of points

            call get_grid_boundary_g06(nx, ny, unt, fgeom, lr, rb, x, y) ! Output: last two

         case (7) ! Power law-cylinder with a GP distribution of points

            call get_grid_boundary_g07(nx, ny, unt, fgeom, lr, rb, x, y) ! Output: last two

         case (8) ! Power law-cylinder with a power law distribution of points

            call get_grid_boundary_g08(nx, ny, unt, fgeom, lr, rb, x, y) ! Output: last two

         case (9) ! Power law-cylinder with a power law and geometric progression distribution of points

            call get_grid_boundary_g09(nx, ny, unt, fgeom, lr, rb, x, y) ! Output: last two

         case (10) ! Power law-cylinder with a power law and geometric progression distribution of points
                   ! Adapted for the DEPP optimizer

            call get_grid_boundary_g10(nx, ny, unt, fgeom, lr, rb, x, y) ! Output: last two

         case (11) ! Bluff power law-cylinder with a geometric progression distribution of points

            call get_grid_boundary_g11(nx, ny, unt, fgeom, lr, rb, x, y) ! Output: last two

         case (12) ! Bluff power law-cylinder with a geometric progression distribution of points
                   ! based on the curve arclength. Adapted to DEPP.

            call get_grid_boundary_g12(nx, ny, unt, fgeom, lr, rb, x, y) ! Output: last two

         case (13) ! Bluff power law-cylinder with a double-exponential and a
                   ! geometric progression distribution of points based on the
                   ! curve arclength. Adapted to DEPP.

            call get_grid_boundary_g13(nx, ny, unt, fgeom, lr, rb, x, y) ! Output: last two

         case (14) ! Newton shape

            call get_grid_boundary_g14(nx, ny, unt, fgeom, lr, rb, x, y) ! Output: last two

         case (15) ! von Karman shape

            call get_grid_boundary_g15(nx, ny, unt, fgeom, lr, rb, x, y) ! Output: last two

         case (16) ! Bluff shifted power law-cylinder with a double-exponential and a
                   ! geometric progression distribution of points based on the
                   ! curve arclength. Adapted to DEPP.

            call get_grid_boundary_g16(nx, ny, unt, fgeom, lr, rb, x, y) ! Output: last two

         case (17) ! Modified Newton shape with two adjustable coefficients

            call get_grid_boundary_g17(nx, ny, unt, fgeom, lr, rb, x, y) ! Output: last two

         case (18) ! Hemisphere-cone-cylinder with concentration of points near the cone-cylinder connection

            call get_grid_boundary_g18(nx, ny, unt, fgeom, lr, rb, x, y) ! Output: last two

         case (19) ! Bluff shifted power law with spline perturbation-cylinder with a double-exponential and a
                   ! geometric progression distribution of points based on the curve arclength. Adapted to DEPP.

            call get_grid_boundary_g19(nx, ny, unt, fgeom, lr, rb, x, y) ! Output: last two

         case (20) ! Bluff shifted power law with exponential perturbation-cylinder with a double-exponential and a
                   ! geometric progression distribution of points based on the curve arclength. Adapted to DEPP.

            call get_grid_boundary_g20(nx, ny, unt, fgeom, lr, rb, x, y) ! Output: last two


         case (21) ! Bluff nose followed by a smooth arc generated by natural cubic splines

            call get_grid_boundary_g21(nx, ny, unt, fgeom, lr, rb, x, y) ! Output: last two


         case (22) ! Bluff nose followed by a smooth arc generated by cubic splines with second
                   ! derivatives and boundary values of the function given.

            call get_grid_boundary_g22(nx, ny, unt, fgeom, lr, rb, x, y) ! Output: last two

         case (23) ! Bluff nose followed by a smooth arc generated by splines with first
                   ! derivatives and boundary conditions on the east boundary given.

            call get_grid_boundary_g23(nx, ny, unt, fgeom, lr, rb, x, y) ! Output: last two


         case default ! Unknown geometry

            write(*,*) "ERROR: Unknown geometry. Stopping."

            stop

      end select

   end subroutine get_grid_boundary

   !> \brief Calculates the coordinates (x,y) of the north and south boundaries.
   !! The north boundary line is an elipse with semi-axis la and lb.
   !! The south boundary line is a power law.
   !! The distribution of points over the lines is a power law.
   subroutine get_grid_boundary_g01(nx, ny, unt, fgeom, lr, rb, x, y)
      implicit none
      integer, intent(in) :: nx  !< Number of volumes in the csi direction (real+fictitious)
      integer, intent(in) :: ny  !< Number of volumes in the eta direction (real+fictitious)
      integer, intent(in) :: unt !< Unit where the input parameters will be printed
      character(len=*), intent(in) :: fgeom !< File of the geometric parameters

      real(8), intent(out) :: lr  !< Length of the body (m)
      real(8), intent(out) :: rb  !< Base radius of the body (m)
      real(8), dimension(nx*ny), intent(out) :: x !< Coord. x at the northeast corner of the volume P (m)
      real(8), dimension(nx*ny), intent(out) :: y !< Coord. y at the northeast corner of the volume P (m)


      ! Constants

      real(8), parameter :: pi = acos(-1.d0)


      ! Inner variables
      integer :: i, j, np ! Dummy indexes

      integer :: kgb ! Kind of geometry of the body (1=Power law)
      real(8) :: akn ! Exponent of the power law for the north boundary
      real(8) :: aks ! Exponent of the power law for the south boundary
      real(8) :: la  ! Length of the elliptical x semi-axis (m)
      real(8) :: lb  ! Length of the elliptical y semi-axis (m)
      real(8) :: lbd ! Exponent of the power law body




      ! Reading parameters

      open(13,file=fgeom)

      read(13,*) kgb  ! Kind of geometry of the body (1=Power law)
      read(13,*) la   ! Length of the elliptical x semi-axis (m)
      read(13,*) lb   ! Length of the elliptical y semi-axis (m)
      read(13,*) lr   ! Length of the body (m)
      read(13,*) rb   ! Base radius of the body (m)
      read(13,*) lbd  ! Exponent of the power law body
      read(13,*) akn  ! Exponent of the power law for the north boundary
      read(13,*) aks  ! Exponent of the power law for the south boundary

      close(13)


      ! Writing parameters

      write(unt,*)
      write(unt,*) " *** Geometry parameters ***"
      write(unt,*)
      write(unt,*) " Power law ogive"
      write(unt,*)

      write(unt,"(I23,    A)") kgb, " = kgb: Kind of geometry of the body (1=Power law)"
      write(unt,"(ES23.16,A)")  la, " =  la: Length of the elliptical x semi-axis (m)"
      write(unt,"(ES23.16,A)")  lb, " =  lb: Length of the elliptical y semi-axis (m)"
      write(unt,"(ES23.16,A)")  lr, " =  lr: Length of the body (m)"
      write(unt,"(ES23.16,A)")  rb, " =  rb: Base radius of the body (m)"
      write(unt,"(ES23.16,A)") lbd, " = ldb: Exponent of the power law body"
      write(unt,"(ES23.16,A)") akn, " = akn: Exponent of the power law for the north boundary"
      write(unt,"(ES23.16,A)") aks, " = aks: Exponent of the power law for the south boundary"



      ! Generating the north boundary

      j = ny-1

      do i = 1, nx-1

         np   = nx * (j-1) + i

         x(np) = la * ( dble(i-1) / dble(nx-2) ) ** akn - (la-lr)

         y(np) = lb * sqrt( 1.d0 - ( (x(np)-lr) / la )**2.d0 )

      end do



      ! Generating the south boundary

      j = 1

      do i = 1, nx-1

         np   = nx * (j-1) + i

         x(np) = lr * ( dble(i-1) / dble(nx-2) ) ** aks

         y(np) = g( x(np) )

      end do


      contains

         !>\brief Calculates the body geometry
         real(8) function g(x)
            implicit none
            real(8), intent(in) :: x

            g = rb * ( x / lr ) ** lbd

         end function g


   end subroutine get_grid_boundary_g01



   !> \brief Calculates the coordinates (x,y) of the north and south boundaries.
   !! The north boundary line is an elipse with semi-axis la and lb.
   !! The south boundary line is an hemisphere-cone.
   !! The distribution of points over the lines is a power law.
   subroutine get_grid_boundary_g02(nx, ny, unt, fgeom, lr, rb, x, y)
      implicit none
      integer, intent(in) :: nx  !< Number of volumes in the csi direction (real+fictitious)
      integer, intent(in) :: ny  !< Number of volumes in the eta direction (real+fictitious)
      integer, intent(in) :: unt !< Unit where the input parameters will be printed
      character(len=*), intent(in) :: fgeom !< File of the geometric parameters

      real(8), intent(out) :: lr  !< Length of the body (m)
      real(8), intent(out) :: rb  !< Base radius of the body (m)
      real(8), dimension(nx*ny), intent(out) :: x !< Coord. x at the northeast corner of the volume P (m)
      real(8), dimension(nx*ny), intent(out) :: y !< Coord. y at the northeast corner of the volume P (m)


      ! Constants

      real(8), parameter :: pi = acos(-1.d0)


      ! Inner variables
      integer :: i, j, np ! Dummy indexes

      integer :: kgb ! Kind of geometry of the body (2=Hemisphere-cone)
      real(8) :: akn ! Exponent of the power law for the north boundary
      real(8) :: aks ! Exponent of the power law for the south boundary
      real(8) :: la  ! Length of the elliptical x semi-axis (m)
      real(8) :: lb  ! Length of the elliptical y semi-axis (m)
      real(8) :: tht ! Semi-angle of the cone (degrees)
      real(8) :: ac  ! Tangent of the cone semi-angle
      real(8) :: xm  ! x of the matching point between the hemisphere and the cone (m)
      real(8) :: rh  ! Radius of the hemisphere (m)




      ! Reading parameters

      open(13,file=fgeom)

      read(13,*) kgb  ! Kind of geometry of the body (2=Hemisphere-cone)
      read(13,*) la   ! Length of the elliptical x semi-axis (m)
      read(13,*) lb   ! Length of the elliptical y semi-axis (m)
      read(13,*) lr   ! Length of the body (m)
      read(13,*) rb   ! Base radius of the body (m)
      read(13,*) tht  ! Semi-angle of the cone (degrees)
      read(13,*) akn  ! Exponent of the power law for the north boundary
      read(13,*) aks  ! Exponent of the power law for the south boundary

      close(13)


      ! Writing parameters

      write(unt,*)
      write(unt,*) " *** Geometry parameters ***"
      write(unt,*)

      write(unt,"(I23,    A)") kgb, " = kgb: Kind of geometry of the body (2=Hemisphere-cone)"
      write(unt,"(ES23.16,A)")  la, " =  la: Length of the elliptical x semi-axis (m)"
      write(unt,"(ES23.16,A)")  lb, " =  lb: Length of the elliptical y semi-axis (m)"
      write(unt,"(ES23.16,A)")  lr, " =  lr: Length of the body (m)"
      write(unt,"(ES23.16,A)")  rb, " =  rb: Base radius of the body (m)"
      write(unt,"(ES23.16,A)") tht, " = tht: Semi-angle of the cone (degrees)"
      write(unt,"(ES23.16,A)") akn, " = akn: Exponent of the power law for the north boundary"
      write(unt,"(ES23.16,A)") aks, " = aks: Exponent of the power law for the south boundary"



      ! Calculating coefficients necessary to the body calculation

      ! Tangent of the cone semi-angle
      ac = tan( tht * pi / 180.d0 )

      ! x of the matching point between the hemisphere and the cone (m)
      xm = ( rb - ac * lr ) / sqrt( 1.d0 + ac * ac )

      if ( xm < 0.d0 ) then

         write(*,*) "ERROR: Semi-angle of the cone is greater than the allowed" &
            // " value (degree): ", atan(rb/lr)*180.d0/pi

         stop

      end if

      ! radius of the hemisphere (m)
      rh = xm + ac * ( ac * (xm-lr) + rb )


      write(unt,"(ES23.16,A)") xm, " =  xm: x of the matching point between the hemisphere and the cone (m)"
      write(unt,"(ES23.16,A)") rh, " =  rh: radius of the hemisphere (m)"


      ! Generating the north boundary

      j = ny-1

      do i = 1, nx-1

         np   = nx * (j-1) + i

         x(np) = la * ( dble(i-1) / dble(nx-2) ) ** akn - (la-lr)

         y(np) = lb * sqrt( 1.d0 - ( (x(np)-lr) / la )**2.d0 )

      end do



      ! Generating the south boundary

      j = 1

      do i = 1, nx-1

         np   = nx * (j-1) + i

         x(np) = lr * ( dble(i-1) / dble(nx-2) ) ** aks

         y(np) = g( x(np) )

      end do


      contains

         !>\brief Calculates the body geometry
         real(8) function g(x)
            implicit none
            real(8), intent(in) :: x

            if ( x < xm ) then

               g = sqrt( rh * rh - ( x - rh ) ** 2 )

            else

               g = ac * ( x - lr ) + rb

            end if

         end function g


   end subroutine get_grid_boundary_g02


   !> \brief Calculates the coordinates (x,y) of the north and south boundaries.
   !! The north boundary line is an elipse with semi-axis lf+lo and lb followed
   !! by a cylinder of radius lb.
   !! The south boundary line is an hemisphere-cone of length lo, radius rb and
   !! semi-angle tht followed by a cylinder of radius rb.
   !! The distribution of points over the ellipse and hemisphere-cone is a power
   !! law.
   subroutine get_grid_boundary_g03(nx, ny, unt, fgeom, lr, rb, x, y) ! Output: last two
      implicit none
      integer, intent(in) :: nx  !< Number of volumes in the csi direction (real+fictitious)
      integer, intent(in) :: ny  !< Number of volumes in the eta direction (real+fictitious)
      integer, intent(in) :: unt !< Unit where the input parameters will be printed
      character(len=*), intent(in) :: fgeom !< File of the geometric parameters

      real(8), intent(out) :: lr  !< Length of the body (m)
      real(8), intent(out) :: rb  !< Base radius of the body (m)
      real(8), dimension(nx*ny), intent(out) :: x !< Coord. x at the northeast corner of the volume P (m)
      real(8), dimension(nx*ny), intent(out) :: y !< Coord. y at the northeast corner of the volume P (m)


      ! Constants

      real(8), parameter :: pi = acos(-1.d0)


      ! Inner variables
      integer :: i, j, np, npw ! Dummy indexes

      integer :: kgb ! Kind of geometry of the body (3=Hemisphere-cone-cylinder)
      integer :: nx1 ! Number of points in the first part of the grid along the eta lines
      integer :: iiv ! Initial value of i
      integer :: ifv ! Final value of i

      real(8) :: a2  ! Width of the last volume of Part I closer to the wall

      real(8) :: xiv ! Initial value of x
      real(8) :: xfv ! Final value of x

      real(8) :: ac  ! Tangent of tht
      real(8) :: rh  ! Hemisphere radius (m)
      real(8) :: xm  ! x coordinate of intersection between the hemisphere and the cone(m)

      real(8) :: lf  ! Frontal length (m)
      real(8) :: lo  ! Ogive length (m)
      real(8) :: lb  ! Base length (m)
      real(8) :: lc  ! Cylinder length (m)
      real(8) :: tht ! Cone semi-angle (degrees)
      real(8) :: akn ! Exponent for the distribution of ksi lines over the north boundary
      real(8) :: aks ! Exponent for the distribution of ksi lines over the south boundary
      real(8) :: akc ! Exponent of the power law for the south and north cylinder boundaries
      real(8) :: fp1 ! Fraction of point in the first part of the grid

      ! Reading parameters

      open(13,file=fgeom)

      read(13,*) kgb  ! Kind of geometry of the body (2=Hemisphere-cone)
      read(13,*) lf   ! Frontal length (m)
      read(13,*) lb   ! Base length (m)
      read(13,*) lo   ! Length of the ogive (m)
      read(13,*) lc   ! Length of the cylinder (m)
      read(13,*) rb   ! Base radius of the body (m)
      read(13,*) tht  ! Semi-angle of the cone (degrees)
      read(13,*) akn  ! Exponent of the power law for the north boundary
      read(13,*) aks  ! Exponent of the power law for the south boundary
      read(13,*) akc  ! Exponent of the power law for the south and north cylinder boundaries
      read(13,*) fp1  ! Fraction of volumes in the hemisphere-cone part

      close(13)


      ! Writing parameters

      write(unt,*)
      write(unt,*) " *** Geometry parameters ***"
      write(unt,*)

      write(unt,"(I23,    A)") kgb, " = kgb: Kind of geometry of the body (3=Hemisphere-cone-cylinder)"
      write(unt,"(ES23.16,A)")  lf, " =  lf: Frontal length (m)"
      write(unt,"(ES23.16,A)")  lb, " =  lb: Base length (m)"
      write(unt,"(ES23.16,A)")  lo, " =  lo: Length of the ogive (m)"
      write(unt,"(ES23.16,A)")  lc, " =  lc: Length of the cylinder (m)"
      write(unt,"(ES23.16,A)")  rb, " =  rb: Base radius of the body (m)"
      write(unt,"(ES23.16,A)") tht, " = tht: Semi-angle of the cone (degrees)"
      write(unt,"(ES23.16,A)") akn, " = akn: Exponent of the power law for the north boundary"
      write(unt,"(ES23.16,A)") aks, " = aks: Exponent of the power law for the south boundary"
      write(unt,"(ES23.16,A)") akc, " = akc: Exponent of the power law for the south and north cylinder boundaries"
      write(unt,"(ES23.16,A)") fp1, " = fp1: Fraction of volumes in the hemisphere-cone part"


      ! Calculating coefficients necessary to the body calculation

      ! Tangent of the cone semi-angle
      ac = tan( tht * pi / 180.d0 )

      ! x of the matching point between the hemisphere and the cone (m)
      xm = ( rb - ac * lo ) / sqrt( 1.d0 + ac * ac )

      if ( xm < 0.d0 ) then

         write(*,*) "ERROR: Semi-angle of the cone is greater than the allowed" &
            // " value (degree): ", atan(rb/lo)*180.d0/pi

         stop

      end if

      ! radius of the hemisphere (m)
      rh = xm + ac * ( ac * (xm-lo) + rb )


      write(unt,"(ES23.16,A)") xm, " =  xm: x of the matching point between the hemisphere and the cone (m)"
      write(unt,"(ES23.16,A)") rh, " =  rh: radius of the hemisphere (m)"


      ! Calculating the number of points in the first part of the grid

      nx1 = int( nx * fp1 )

      ! Calculating the body length

      lr = lo + lc


      ! ========================================================

      ! Generating north boundary

      j = ny-1

      ! Part I

      iiv = 1   ! Initial value of i

      ifv = nx-1 ! Final value of i

      xiv = -lf ! Initial value of x

      xfv =  lo+lc ! Final value of x

      do i = iiv, ifv

         np   = nx * (j-1) + i

         x(np) = (xfv-xiv) * ( dble(i - iiv) / dble( ifv - iiv ) ) ** akn + xiv

         y(np) = lb * sqrt( 1.d0 - ( (x(np)-lo-lc) / ( lf + lo + lc ) ) ** 2.d0 )

      end do

      ! ========================================================



      ! ========================================================

      ! Generating south boundary

      j = 1

      ! Part I

      iiv = 1   ! Initial value of i

      ifv = nx1 ! Final value of i

      xiv = 0.d0 ! Initial value of x

      xfv = lo   ! Final value of x

      do i = iiv, ifv

         np   = nx * (j-1) + i

         x(np) = (xfv-xiv) * ( dble(i - iiv) / dble( ifv - iiv ) ) ** aks + xiv

         y(np) = fhc(x(np), xm, ac, lo, rb, rh)

      end do


      ! Part II

      iiv = nx1  ! Initial value of i

      ifv = nx-1 ! Final value of i

      xiv = lo      ! Initial value of x

      xfv = lo + lc ! Final value of x


      ! Automatic calculation of akc
      if ( akc < 0.d0 ) then

         i = iiv

         np   = nx * (j-1) + i

         npw  = np - 1

         ! Calculating the width of the last volume of Part I

         a2 = sqrt( (x(np)-x(npw))**2 + (y(np)-y(npw))**2 )

         ! Calculating akc

         akc = log( lc / a2 ) / log( dble(ifv-iiv) )

      end if

      do i = iiv, ifv

         np   = nx * (j-1) + i

         x(np) = (xfv-xiv) * ( dble(i - iiv) / dble( ifv - iiv ) ) ** akc + xiv

         y(np) = rb

      end do

      ! ========================================================

   contains

      !> \brief Hemisphere-cone function
      real(8) function fhc(x, xm, ac, lo, ro, rh)
         implicit none
         real(8), intent(in) :: x  !< x coordinate (m)
         real(8), intent(in) :: xm !< x coordinate of intersection between the hemisphere and the cone(m)
         real(8), intent(in) :: ac !< Tangent of tht
         real(8), intent(in) :: lo !< Length of the ogive (m)
         real(8), intent(in) :: ro !< Radius of the ogive (m)
         real(8), intent(in) :: rh !< Radius of the hemisphere (m)

         if ( x < xm ) then

            fhc = sqrt( rh * rh - ( x - rh ) ** 2 )

         else

            fhc = ac * ( x -lo ) + ro

         end if

      end function

   end subroutine get_grid_boundary_g03


   !> \brief Calculates the coordinates (x,y) of the north and south boundaries.
   !! The north boundary line is an elipse with semi-axis la and lb.
   !! The south boundary line is a power law.
   !! The distribution of points over the lines is a power law.
   !! This subroutine is different from get_grid_boundary_g01 because it
   !! allows the optimization of the exponent of the power law through the
   !! DEPP optimizer.
   subroutine get_grid_boundary_g04(nx, ny, unt, fgeom, lr, rb, x, y)
      implicit none
      integer, intent(in) :: nx  !< Number of volumes in the csi direction (real+fictitious)
      integer, intent(in) :: ny  !< Number of volumes in the eta direction (real+fictitious)
      integer, intent(in) :: unt !< Unit where the input parameters will be printed
      character(len=*), intent(in) :: fgeom !< File of the geometric parameters

      real(8), intent(out) :: lr  !< Length of the body (m)
      real(8), intent(out) :: rb  !< Base radius of the body (m)
      real(8), dimension(nx*ny), intent(out) :: x !< Coord. x at the northeast corner of the volume P (m)
      real(8), dimension(nx*ny), intent(out) :: y !< Coord. y at the northeast corner of the volume P (m)


      ! Constants

      real(8), parameter :: pi = acos(-1.d0)


      ! Inner variables
      integer :: i, j, np ! Dummy indexes

      integer :: kgb ! Kind of geometry of the body (1=Power law)
      real(8) :: akn ! Exponent of the power law for the north boundary
      real(8) :: aks ! Exponent of the power law for the south boundary
      real(8) :: la  ! Length of the elliptical x semi-axis (m)
      real(8) :: lb  ! Length of the elliptical y semi-axis (m)
      real(8) :: lbd ! Exponent of the power law body


      ! DEPP variables
      integer, parameter :: nu = 1 ! number of unknowns
      integer            :: ind    ! number of the individual
      integer,       dimension(nu) :: xopt   ! Optimization checker
      character(10), dimension(nu) :: xname  ! Name of parameters
      real(8),       dimension(nu) :: xvalue ! parameters
      character(200)               :: sname  ! simulation name



      ! Reading parameters

      open(13,file=fgeom)

      read(13,*) kgb  ! Kind of geometry of the body (1=Power law)
      read(13,*) la   ! Length of the elliptical x semi-axis (m)
      read(13,*) lb   ! Length of the elliptical y semi-axis (m)
      read(13,*) lr   ! Length of the body (m)
      read(13,*) rb   ! Base radius of the body (m)
      read(13,*) xname(1) ! Parameter name
      read(13,*) xopt(1)  ! Will this parameter be optimized? ( 0 = no, 1 = yes )
      read(13,*) xvalue(1)! Parameter value
      read(13,*) akn  ! Exponent of the power law for the north boundary
      read(13,*) aks  ! Exponent of the power law for the south boundary

      close(13)

      ! Reads the parameters from DEPP
      call depp_get_parameters(nu, xopt, xname, xvalue, ind, sname)



      lbd = xvalue(1)


      ! Writing parameters

      write(unt,*)
      write(unt,*) " *** Geometry parameters ***"
      write(unt,*)
      write(unt,*) " Power law ogive"
      write(unt,*)

      write(unt,"(I23,    A)") kgb, " = kgb: Kind of geometry of the body (1=Power law)"
      write(unt,"(ES23.16,A)")  la, " =  la: Length of the elliptical x semi-axis (m)"
      write(unt,"(ES23.16,A)")  lb, " =  lb: Length of the elliptical y semi-axis (m)"
      write(unt,"(ES23.16,A)")  lr, " =  lr: Length of the body (m)"
      write(unt,"(ES23.16,A)")  rb, " =  rb: Base radius of the body (m)"
      write(unt,"(ES23.16,A)") lbd, " = ldb: Exponent of the power law body"
      write(unt,"(ES23.16,A)") akn, " = akn: Exponent of the power law for the north boundary"
      write(unt,"(ES23.16,A)") aks, " = aks: Exponent of the power law for the south boundary"



      ! Generating the north boundary

      j = ny-1

      do i = 1, nx-1

         np   = nx * (j-1) + i

         x(np) = la * ( dble(i-1) / dble(nx-2) ) ** akn - (la-lr)

         y(np) = lb * sqrt( 1.d0 - ( (x(np)-lr) / la )**2.d0 )

      end do



      ! Generating the south boundary

      j = 1

      do i = 1, nx-1

         np   = nx * (j-1) + i

         x(np) = lr * ( dble(i-1) / dble(nx-2) ) ** aks

         y(np) = g( x(np) )

      end do


      contains

         !>\brief Calculates the body geometry
         real(8) function g(x)
            implicit none
            real(8), intent(in) :: x

            g = rb * ( x / lr ) ** lbd

         end function g


   end subroutine get_grid_boundary_g04


   !> \brief Calculates the coordinates (x,y) of the north and south boundaries.
   !! The north boundary line is an elipse with semi-axis lf+lo and lb followed
   !! by a cylinder of radius lb.
   !! The south boundary line is an power law of length lo and radius rb followed
   !! by a cylinder of radius rb.
   !! The distribution of points over the ellipse and power law is a power
   !! law.
   subroutine get_grid_boundary_g05(nx, ny, unt, fgeom, lr, rb, x, y) ! Output: last two
      implicit none
      integer, intent(in) :: nx  !< Number of volumes in the csi direction (real+fictitious)
      integer, intent(in) :: ny  !< Number of volumes in the eta direction (real+fictitious)
      integer, intent(in) :: unt !< Unit where the input parameters will be printed
      character(len=*), intent(in) :: fgeom !< File of the geometric parameters

      real(8), intent(out) :: lr  !< Length of the body (m)
      real(8), intent(out) :: rb  !< Base radius of the body (m)
      real(8), dimension(nx*ny), intent(out) :: x !< Coord. x at the northeast corner of the volume P (m)
      real(8), dimension(nx*ny), intent(out) :: y !< Coord. y at the northeast corner of the volume P (m)


      ! Constants

      real(8), parameter :: pi = acos(-1.d0)


      ! Inner variables
      integer :: i, j, np, npw ! Dummy indexes

      integer :: kgb ! Kind of geometry of the body (5=Power-law-cylinder)
      integer :: nx1 ! Number of points in the first part of the grid along the eta lines
      integer :: iiv ! Initial value of i
      integer :: ifv ! Final value of i

      real(8) :: a2  ! Width of the last volume of Part I closer to the wall

      real(8) :: xiv ! Initial value of x
      real(8) :: xfv ! Final value of x

      real(8) :: lbd ! Exponent of the power law body

      real(8) :: lf  ! Frontal length (m)
      real(8) :: lo  ! Ogive length (m)
      real(8) :: lb  ! Base length (m)
      real(8) :: lc  ! Cylinder length (m)
      real(8) :: akn ! Exponent for the distribution of ksi lines over the north boundary
      real(8) :: aks ! Exponent for the distribution of ksi lines over the south boundary
      real(8) :: akc ! Exponent of the power law for the south and north cylinder boundaries
      real(8) :: fp1 ! Fraction of point in the first part of the grid

      ! Reading parameters

      open(13,file=fgeom)

      read(13,*) kgb  ! Kind of geometry of the body (5=Power law-cylinder)
      read(13,*) lf   ! Frontal length (m)
      read(13,*) lb   ! Base length (m)
      read(13,*) lo   ! Length of the ogive (m)
      read(13,*) lc   ! Length of the cylinder (m)
      read(13,*) rb   ! Base radius of the body (m)
      read(13,*) lbd  ! Exponent of the power law body
      read(13,*) akn  ! Exponent of the power law for the north boundary
      read(13,*) aks  ! Exponent of the power law for the south boundary
      read(13,*) akc  ! Exponent of the power law for the south and north cylinder boundaries
      read(13,*) fp1  ! Fraction of volumes in the hemisphere-cone part

      close(13)


      ! Writing parameters

      write(unt,*)
      write(unt,*) " *** Geometry parameters ***"
      write(unt,*)

      write(unt,"(I23,    A)") kgb, " = kgb: Kind of geometry of the body  (5=Power law-cylinder)"
      write(unt,"(ES23.16,A)")  lf, " =  lf: Frontal length (m)"
      write(unt,"(ES23.16,A)")  lb, " =  lb: Base length (m)"
      write(unt,"(ES23.16,A)")  lo, " =  lo: Length of the ogive (m)"
      write(unt,"(ES23.16,A)")  lc, " =  lc: Length of the cylinder (m)"
      write(unt,"(ES23.16,A)")  rb, " =  rb: Base radius of the body (m)"
      write(unt,"(ES23.16,A)") lbd, " = lbd: Exponent of the power law body"
      write(unt,"(ES23.16,A)") akn, " = akn: Exponent of the power law for the north boundary"
      write(unt,"(ES23.16,A)") aks, " = aks: Exponent of the power law for the south boundary"
      write(unt,"(ES23.16,A)") akc, " = akc: Exponent of the power law for the south and north cylinder boundaries"
      write(unt,"(ES23.16,A)") fp1, " = fp1: Fraction of volumes in the hemisphere-cone part"



      ! Calculating the number of points in the first part of the grid

      nx1 = int( nx * fp1 )

      ! Calculating the body length

      lr = lo + lc


      ! ========================================================

      ! Generating north boundary

      j = ny-1

      ! Part I

      iiv = 1   ! Initial value of i

      ifv = nx-1 ! Final value of i

      xiv = -lf ! Initial value of x

      xfv =  lo+lc ! Final value of x

      do i = iiv, ifv

         np   = nx * (j-1) + i

         x(np) = (xfv-xiv) * ( dble(i - iiv) / dble( ifv - iiv ) ) ** akn + xiv

         y(np) = lb * sqrt( 1.d0 - ( (x(np)-lo-lc) / ( lf + lo + lc ) ) ** 2.d0 )

      end do

      ! ========================================================



      ! ========================================================

      ! Generating south boundary

      j = 1

      ! Part I

      iiv = 1   ! Initial value of i

      ifv = nx1 ! Final value of i

      xiv = 0.d0 ! Initial value of x

      xfv = lo   ! Final value of x

      do i = iiv, ifv

         np   = nx * (j-1) + i

         x(np) = (xfv-xiv) * ( dble(i - iiv) / dble( ifv - iiv ) ) ** aks + xiv

         y(np) = rb * ( x(np) / lo ) ** lbd

      end do


      ! Part II

      iiv = nx1  ! Initial value of i

      ifv = nx-1 ! Final value of i

      xiv = lo      ! Initial value of x

      xfv = lo + lc ! Final value of x


      ! Automatic calculation of akc
      if ( akc < 0.d0 ) then

         i = iiv

         np   = nx * (j-1) + i

         npw  = np - 1

         ! Calculating the width of the last volume of Part I

         a2 = sqrt( (x(np)-x(npw))**2 + (y(np)-y(npw))**2 )

         ! Calculating akc

         akc = log( lc / a2 ) / log( dble(ifv-iiv) )

      end if

      do i = iiv, ifv

         np   = nx * (j-1) + i

         x(np) = (xfv-xiv) * ( dble(i - iiv) / dble( ifv - iiv ) ) ** akc + xiv

         y(np) = rb

      end do

   end subroutine get_grid_boundary_g05


   !> \brief Calculates the coordinates (x,y) of the north and south boundaries.
   !! The north boundary line is an elipse with semi-axis lf+lo and lb followed
   !! by a cylinder of radius lb.
   !! The south boundary line is an power law of length lo and radius rb followed
   !! by a cylinder of radius rb.
   !! The distribution of points over the south and north boundary is non-uniform
   subroutine get_grid_boundary_g06(nx, ny, unt, fgeom, lr, rb, x, y) ! Output: last two
      implicit none
      integer, intent(in) :: nx  !< Number of volumes in the csi direction (real+fictitious)
      integer, intent(in) :: ny  !< Number of volumes in the eta direction (real+fictitious)
      integer, intent(in) :: unt !< Unit where the input parameters will be printed
      character(len=*), intent(in) :: fgeom !< File of the geometric parameters

      real(8), intent(out) :: lr  !< Length of the body (m)
      real(8), intent(out) :: rb  !< Base radius of the body (m)
      real(8), dimension(nx*ny), intent(out) :: x !< Coord. x at the northeast corner of the volume P (m)
      real(8), dimension(nx*ny), intent(out) :: y !< Coord. y at the northeast corner of the volume P (m)


      ! Constants

      real(8), parameter :: pi = acos(-1.d0)


      ! Inner variables
      integer :: i, j, np ! Dummy indexes

      real(8) :: Ac, Bc, Cc

      integer :: kgb ! Kind of geometry of the body (5=Power-law-cylinder)
      integer :: nx1 ! Number of points in the first part of the grid along the eta lines
      integer :: iiv ! Initial value of i
      integer :: ifv ! Final value of i

      real(8) :: a1  ! Width of the first volume of Part I closer to the wall
      real(8) :: a2  ! Width of the last volume of Part I closer to the wall
      real(8) :: Ca2 ! a2 = Ca2 * a1

      real(8) :: xiv ! Initial value of x
      real(8) :: xfv ! Final value of x

      real(8) :: lbd ! Exponent of the power law body

      real(8) :: lf  ! Frontal length (m)
      real(8) :: lo  ! Ogive length (m)
      real(8) :: lb  ! Base length (m)
      real(8) :: lc  ! Cylinder length (m)
      real(8) :: aks ! Exponent for the distribution of ksi lines over the south boundary
      real(8) :: akc ! Exponent of the power law for the south and north cylinder boundaries
      real(8) :: fp1 ! Fraction of point in the first part of the grid

      real(8), dimension(nx) :: xbs
      real(8), dimension(nx) :: xbn

      ! Reading parameters

      open(13,file=fgeom)

      read(13,*) kgb  ! Kind of geometry of the body (5=Power law-cylinder)
      read(13,*) lf   ! Frontal length (m)
      read(13,*) lb   ! Base length (m)
      read(13,*) lo   ! Length of the ogive (m)
      read(13,*) lc   ! Length of the cylinder (m)
      read(13,*) rb   ! Base radius of the body (m)
      read(13,*) lbd  ! Exponent of the power law body
      read(13,*) Ca2  ! a2 = Ca2 * a1
      read(13,*) aks  ! Exponent of the power law for the south boundary
      read(13,*) akc  ! Exponent of the power law for the south and north cylinder boundaries
      read(13,*) fp1  ! Fraction of volumes in the hemisphere-cone part

      close(13)


      ! Writing parameters

      write(unt,*)
      write(unt,*) " *** Geometry parameters ***"
      write(unt,*)

      write(unt,"(I23,    A)") kgb, " = kgb: Kind of geometry of the body  (5=Power law-cylinder)"
      write(unt,"(ES23.16,A)")  lf, " =  lf: Frontal length (m)"
      write(unt,"(ES23.16,A)")  lb, " =  lb: Base length (m)"
      write(unt,"(ES23.16,A)")  lo, " =  lo: Length of the ogive (m)"
      write(unt,"(ES23.16,A)")  lc, " =  lc: Length of the cylinder (m)"
      write(unt,"(ES23.16,A)")  rb, " =  rb: Base radius of the body (m)"
      write(unt,"(ES23.16,A)") lbd, " = lbd: Exponent of the power law body"
      write(unt,"(ES23.16,A)") Ca2, " = Ca2: a2 = Ca2 * a1"
      write(unt,"(ES23.16,A)") aks, " = aks: Exponent of the power law for the south boundary"
      write(unt,"(ES23.16,A)") akc, " = akc: Exponent of the power law for the south and north cylinder boundaries"
      write(unt,"(ES23.16,A)") fp1, " = fp1: Fraction of volumes in the hemisphere-cone part"


      ! Calculating the number of points in the first part of the grid

      nx1 = int( nx * fp1 )

      ! Calculating the body length

      lr = lo + lc



      ! ========================================================

      ! Generating south boundary

      j = 1

      ! Part I

      iiv = 1   ! Initial value of i

      ifv = nx1 ! Final value of i

      xiv = 0.d0 ! Initial value of x

      xfv = lo   ! Final value of x

      a1 = (xfv-xiv) / (ifv-iiv) ** aks

      a2 = Ca2 * a1

      call get_quadratic_distribution(iiv, ifv, xiv, xfv, a1, a2, xbs(iiv:ifv) ) ! Output: last three

      do i = iiv, ifv

         np   = nx * (j-1) + i

         x(np) = xbs(i)

         y(np) = rb * ( x(np) / lo ) ** lbd

      end do


      ! Part II

      iiv = nx1  ! Initial value of i

      ifv = nx-1 ! Final value of i

      xiv = lo      ! Initial value of x

      xfv = lo + lc ! Final value of x


      ! Automatic calculation of akc
      if ( akc < 0.d0 ) then

         ! Calculating akc

         akc = log( lc / a2 ) / log( dble(ifv-iiv) )

      end if

      do i = iiv, ifv

         np   = nx * (j-1) + i

         x(np) = (xfv-xiv) * ( dble(i - iiv) / dble( ifv - iiv ) ) ** akc + xiv

         xbs(i) = x(np)

         y(np) = rb

      end do





      ! ========================================================

      ! Generating north boundary

      j = ny-1

      ! Part I

      iiv = 1   ! Initial value of i

      ifv = nx-1 ! Final value of i

      xiv = -lf ! Initial value of x

      xfv =  lo+lc ! Final value of x

      xbn(iiv) = xiv

      do i = iiv+1, ifv

         xbn(i) = xbn(i-1) + (lf+lr)/lr * (xbs(i)-xbs(i-1))

      end do


      do i = iiv, ifv

         np   = nx * (j-1) + i

         x(np) = xbn(i)

         y(np) = lb * sqrt( 1.d0 - ( (x(np)-lo-lc) / ( lf + lo + lc ) ) ** 2.d0 )

      end do

      ! ========================================================

   contains

      !> \brief Calculates de coefficients for a quadratic distribution of points
      !! where the width of the first and last volume is prescribed
      subroutine get_quadratic_distribution(iiv, ifv, xiv, xfv, a1, a2, xbs) ! Output: last three
         implicit none
         integer, intent(in)  :: iiv !< Initial value of i
         integer, intent(in)  :: ifv !< Final value of i
         real(8), intent(in)  :: xiv !< Initial value of x
         real(8), intent(in)  :: xfv !< Final value of x
         real(8), intent(in)  :: a1  !< Width of the first volume
         real(8), intent(in)  :: a2  !< Width of the last volume
         real(8), intent(out) :: xbs(iiv:ifv) !< Distribution of points

         ! Inner variables

         integer :: i

         real(8) :: S1, S2, S3


         S3 = dble(ifv - iiv)

         S2 = 0.d0

         S1 = 0.d0

         do i = iiv+1, ifv

            S2 = S2 + i - iiv - 1

            S1 = S1 + ( i - iiv - 1 ) * ( i - ifv )

         end do


         Cc = a1

         Bc = ( a2 - a1 ) / ( ifv - iiv - 1 )

         Ac = ( xfv - xiv - Bc * S2 - Cc * S3 ) / S1

         xbs(iiv) = xiv

         do i = iiv+1, ifv

            xbs(i) = xbs(i-1) &

               + Ac * ( i - iiv - 1 ) * ( i - ifv ) &

               + Bc * ( i - iiv - 1 ) &

               + Cc

         end do

      end subroutine get_quadratic_distribution


   end subroutine get_grid_boundary_g06


   !> \brief Calculates the coordinates (x,y) of the north and south boundaries.
   !! The north boundary line is an elipse with semi-axis lf+lr and lb followed
   !! by a cylinder of radius lb.
   !! The south boundary line is an power law of length lo and radius rb followed
   !! by a cylinder of radius rb and length lc.
   !! The distribution of points over the ellipse and power law is a power
   !! law.
   subroutine get_grid_boundary_g07(nx, ny, unt, fgeom, lr, rb, x, y) ! Output: last two
      implicit none
      integer, intent(in) :: nx  !< Number of volumes in the csi direction (real+fictitious)
      integer, intent(in) :: ny  !< Number of volumes in the eta direction (real+fictitious)
      integer, intent(in) :: unt !< Unit where the input parameters will be printed
      character(len=*), intent(in) :: fgeom !< File of the geometric parameters

      real(8), intent(out) :: lr  !< Length of the body (m)
      real(8), intent(out) :: rb  !< Base radius of the body (m)
      real(8), dimension(nx*ny), intent(out) :: x !< Coord. x at the northeast corner of the volume P (m)
      real(8), dimension(nx*ny), intent(out) :: y !< Coord. y at the northeast corner of the volume P (m)


      ! Constants

      real(8), parameter :: pi = acos(-1.d0)


      ! Inner variables
      integer :: i, j, np, npw ! Dummy indexes

      integer :: kgb ! Kind of geometry of the body (5=Power-law-cylinder)
      integer :: nx1 ! Number of points in the first part of the grid along the eta lines
      integer :: nx2 ! Number of points in the first and second parts of the grid along the eta lines
      integer :: iiv ! Initial value of i
      integer :: ifv ! Final value of i

      real(8) :: a1  ! Width of the first volume of Part I
      real(8) :: a2  ! Width of the first volume of Part II
      real(8) :: a3  ! Width of the first volume of Part III

      real(8) :: xiv ! Initial value of x
      real(8) :: xfv ! Final value of x

      real(8) :: lbd ! Exponent of the power law body

      real(8) :: lf  ! Frontal length (m)
      real(8) :: lo  ! Ogive length (m)
      real(8) :: lb  ! Base length (m)
      real(8) :: lc  ! Cylinder length (m)
      real(8) :: lp1 ! Length of the first part of the grid (m)
      real(8) :: aks ! Exponent for the distribution of ksi lines over the south boundary
      real(8) :: fp1 ! Fraction of points in the first part of the grid
      real(8) :: fp2 ! Fraction of points in the second part of the grid


      real(8), dimension(nx) :: xbs ! Distribution of points over the south boundary
      real(8), dimension(nx) :: xbn ! Distribution of points over the north boundary


      ! Reading parameters

      open(13,file=fgeom)

      read(13,*) kgb  ! Kind of geometry of the body (5=Power law-cylinder)
      read(13,*) lf   ! Frontal length (m)
      read(13,*) lb   ! Base length (m)
      read(13,*) lo   ! Length of the ogive (m)
      read(13,*) lc   ! Length of the cylinder (m)
      read(13,*) rb   ! Base radius of the body (m)
      read(13,*) lp1  ! Length of the first part ( 0 < lp1 < lo ) (m)
      read(13,*) lbd  ! Exponent of the power law body
      read(13,*) aks  ! Exponent of the power law for the south boundary
      read(13,*) fp1  ! Fraction of volumes in the first part
      read(13,*) fp2  ! Fraction of volumes in the second part

      close(13)


      ! Writing parameters

      write(unt,*)
      write(unt,*) " *** Geometry parameters ***"
      write(unt,*)

      write(unt,"(I23,    A)") kgb, " = kgb: Kind of geometry of the body  (5=Power law-cylinder)"
      write(unt,"(ES23.16,A)")  lf, " =  lf: Frontal length (m)"
      write(unt,"(ES23.16,A)")  lb, " =  lb: Base length (m)"
      write(unt,"(ES23.16,A)")  lo, " =  lo: Length of the ogive (m)"
      write(unt,"(ES23.16,A)")  lc, " =  lc: Length of the cylinder (m)"
      write(unt,"(ES23.16,A)")  rb, " =  rb: Base radius of the body (m)"
      write(unt,"(ES23.16,A)") lp1, " = lp1: length of the first part (m)"
      write(unt,"(ES23.16,A)") lbd, " = lbd: Exponent of the power law body"
      write(unt,"(ES23.16,A)") aks, " = aks: Exponent of the power law for the south boundary"
      write(unt,"(ES23.16,A)") fp1, " = fp1: Fraction of volumes in the first part"
      write(unt,"(ES23.16,A)") fp2, " = fp2: Fraction of volumes in the second part"


      ! Calculating the number of points in the first part of the grid

      nx1 = int( nx * fp1 )

      ! Calculating the number of points in the first and second parts of the grid

      nx2 = int( nx * fp2 ) + nx1


      ! Calculating the body length

      lr = lo + lc


      ! ========================================================

      ! Generating south boundary

      j = 1

      ! Part I

      iiv = 1   ! Initial value of i

      ifv = nx1 ! Final value of i

      xiv = 0.d0 ! Initial value of x

      xfv = lp1  ! Final value of x

      a1 = (xfv-xiv) / (ifv-iiv) ** aks ! Width of the first volume of the first part of the grid

      call get_gp_distribution(iiv, ifv, xiv, xfv, a1, xbs(iiv:ifv)) ! Output: last one

      do i = iiv, ifv

         np   = nx * (j-1) + i

         x(np) = xbs(i)

         y(np) = rb * ( x(np) / lo ) ** lbd

      end do


      ! Part II

      iiv = nx1  ! Initial value of i

      ifv = nx2  ! Final value of i

      xiv = lp1 ! Initial value of x

      xfv = lo  ! Final value of x

      ! Calculating the width of the last volume of Part I

      i = iiv

      np   = nx * (j-1) + i

      npw  = np - 1

      a2 = sqrt( (x(np)-x(npw))**2 + (y(np)-y(npw))**2 )

      call get_gp_distribution(iiv, ifv, xiv, xfv, a2, xbs(iiv:ifv)) ! Output: last one

      do i = iiv, ifv

         np   = nx * (j-1) + i

         x(np) = xbs(i)

         y(np) = rb * ( x(np) / lo ) ** lbd

      end do


      ! Part III

      iiv = nx2  ! Initial value of i

      ifv = nx-1 ! Final value of i

      xiv = lo      ! Initial value of x

      xfv = lo + lc ! Final value of x

      ! Calculating the width of the last volume of Part II

      i = iiv

      np   = nx * (j-1) + i

      npw  = np - 1

      a3 = sqrt( (x(np)-x(npw))**2 + (y(np)-y(npw))**2 )

      call get_gp_distribution(iiv, ifv, xiv, xfv, a3, xbs(iiv:ifv)) ! Output: last one

      do i = iiv, ifv

         np   = nx * (j-1) + i

         x(np) = xbs(i)

         y(np) = rb

      end do


      ! ========================================================

      ! Generating north boundary

      j = ny-1

      iiv = 1   ! Initial value of i

      ifv = nx-1 ! Final value of i

      xiv = -lf ! Initial value of x

      xfv =  lr ! Final value of x

      xbn(iiv) = xiv

      do i = iiv+1, ifv

         xbn(i) = xbn(i-1) + (lf+lr) / lr * ( xbs(i)-xbs(i-1))

      end do

      do i = iiv, ifv

         np   = nx * (j-1) + i

         x(np) = xbn(i)

         y(np) = lb * sqrt( 1.d0 - ( (x(np)-lr) / ( lf + lr ) ) ** 2.d0 )

      end do

      ! ========================================================

   contains

      !> \brief Calculates the node distribution based on a geometric progression distribution
      subroutine get_gp_distribution(iiv, ifv, xiv, xfv, a1, xb) ! Output: last one
         implicit none
         integer, intent(in)  :: iiv !< Initial value of i
         integer, intent(in)  :: ifv !< Final value of i
         real(8), intent(in)  :: xiv !< Initial value of x
         real(8), intent(in)  :: xfv !< Final value of x
         real(8), intent(in)  :: a1  !< Width of the first volume
         real(8), intent(out) :: xb(iiv:ifv) !< Distribution of points

         ! Inner variables

         integer :: i
         integer :: n
         real(8) :: r
         real(8) :: q


         n = ifv - iiv

         r = ( xfv - xiv ) / a1

         ! Calculating q
         call get_gp_ratio(n, r, q)

         xb(iiv) = xiv

         do i = iiv+1, ifv-1

            xb(i) = xb(i-1) + a1 * q ** ( i - iiv - 1 )

         end do

         xb(ifv) = xfv

      end subroutine get_gp_distribution


      !> \brief Calculates the base q of the geometric progression
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

           if ( 0.d0 < f(n, r, qi) * f(n, r, qm) ) then

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

      end subroutine get_gp_ratio

      !> \brief Function used in the subroutine 'get_gp_ratio'
      real(8) function f(n, r, q)
         implicit none
         integer, intent(in)  ::   n !< number of partitions
         real(8), intent(in)  ::   r !< l/a1
         real(8), intent(in)  ::   q !< q

         f = q ** n + r * ( 1.d0 - q ) - 1.d0

      end function f

   end subroutine get_grid_boundary_g07


   !> \brief Calculates the coordinates (x,y) of the north and south boundaries.
   !! The north boundary line is an elipse with semi-axis lf+lr and lb.
   !! The south boundary line is an power law of length lo and radius rb followed
   !! by a cylinder of radius rb and length lc.
   !! The distribution of points over the ellipse and power law is a power
   !! law.
   subroutine get_grid_boundary_g08(nx, ny, unt, fgeom, lr, rb, x, y) ! Output: last two
      implicit none
      integer, intent(in) :: nx  !< Number of volumes in the csi direction (real+fictitious)
      integer, intent(in) :: ny  !< Number of volumes in the eta direction (real+fictitious)
      integer, intent(in) :: unt !< Unit where the input parameters will be printed
      character(len=*), intent(in) :: fgeom !< File of the geometric parameters

      real(8), intent(out) :: lr  !< Length of the body (m)
      real(8), intent(out) :: rb  !< Base radius of the body (m)
      real(8), dimension(nx*ny), intent(out) :: x !< Coord. x at the northeast corner of the volume P (m)
      real(8), dimension(nx*ny), intent(out) :: y !< Coord. y at the northeast corner of the volume P (m)


      ! Inner variables
      integer :: i, j, np, npw ! Dummy indexes

      integer :: kgb ! Kind of geometry of the body (5=Power-law-cylinder)
      integer :: nx1 ! Number of points in the first part of the grid along the eta lines
      integer :: nx2 ! Number of points in the first and second parts of the grid along the eta lines
      integer :: iiv ! Initial value of i
      integer :: ifv ! Final value of i

      real(8) :: a1  ! Width of the first volume of Part I
      real(8) :: a2  ! Width of the first volume of Part II
      real(8) :: a3  ! Width of the first volume of Part III

      real(8) :: xiv ! Initial value of x
      real(8) :: xfv ! Final value of x

      real(8) :: lbd ! Exponent of the power law body

      real(8) :: lf  ! Frontal length (m)
      real(8) :: lo  ! Ogive length (m)
      real(8) :: lb  ! Base length (m)
      real(8) :: lc  ! Cylinder length (m)
      real(8) :: lp1 ! Length of the first part of the grid (m)
      real(8) :: aks ! Exponent for the distribution of ksi lines over the south boundary
      real(8) :: fp1 ! Fraction of points in the first part of the grid
      real(8) :: fp2 ! Fraction of points in the second part of the grid


      real(8), dimension(nx) :: xbs ! Distribution of points over the south boundary
      real(8), dimension(nx) :: xbn ! Distribution of points over the north boundary


      ! Reading parameters

      open(13,file=fgeom)

      read(13,*) kgb  ! Kind of geometry of the body (5=Power law-cylinder)
      read(13,*) lf   ! Frontal length (m)
      read(13,*) lb   ! Base length (m)
      read(13,*) lo   ! Length of the ogive (m)
      read(13,*) lc   ! Length of the cylinder (m)
      read(13,*) rb   ! Base radius of the body (m)
      read(13,*) lp1  ! Length of the first part ( 0 < lp1 < lo ) (m)
      read(13,*) lbd  ! Exponent of the power law body
      read(13,*) aks  ! Exponent of the power law for the south boundary
      read(13,*) fp1  ! Fraction of volumes in the first part
      read(13,*) fp2  ! Fraction of volumes in the second part

      close(13)


      ! Writing parameters

      write(unt,*)
      write(unt,*) " *** Geometry parameters ***"
      write(unt,*)

      write(unt,"(I23,    A)") kgb, " = kgb: Kind of geometry of the body  (5=Power law-cylinder)"
      write(unt,"(ES23.16,A)")  lf, " =  lf: Frontal length (m)"
      write(unt,"(ES23.16,A)")  lb, " =  lb: Base length (m)"
      write(unt,"(ES23.16,A)")  lo, " =  lo: Length of the ogive (m)"
      write(unt,"(ES23.16,A)")  lc, " =  lc: Length of the cylinder (m)"
      write(unt,"(ES23.16,A)")  rb, " =  rb: Base radius of the body (m)"
      write(unt,"(ES23.16,A)") lp1, " = lp1: length of the first part (m)"
      write(unt,"(ES23.16,A)") lbd, " = lbd: Exponent of the power law body"
      write(unt,"(ES23.16,A)") aks, " = aks: Exponent of the power law for the south boundary"
      write(unt,"(ES23.16,A)") fp1, " = fp1: Fraction of volumes in the first part"
      write(unt,"(ES23.16,A)") fp2, " = fp2: Fraction of volumes in the second part"


      ! Calculating the number of points in the first part of the grid

      nx1 = int( nx * fp1 )

      ! Calculating the number of points in the first and second parts of the grid

      nx2 = int( nx * fp2 ) + nx1


      ! Calculating the body length

      lr = lo + lc


      ! ========================================================

      ! Generating south boundary

      j = 1

      ! Part I

      iiv = 1   ! Initial value of i

      ifv = nx1 ! Final value of i

      xiv = 0.d0 ! Initial value of x

      xfv = lp1  ! Final value of x

      a1 = (xfv-xiv) / (ifv-iiv) ** aks ! Width of the first volume of the first part of the grid

      call get_pl_distribution(iiv, ifv, xiv, xfv, a1, xbs(iiv:ifv)) ! Output: last one

      do i = iiv, ifv

         np   = nx * (j-1) + i

         x(np) = xbs(i)

         y(np) = rb * ( x(np) / lo ) ** lbd

      end do


      ! Part II

      iiv = nx1  ! Initial value of i

      ifv = nx2  ! Final value of i

      xiv = lp1 ! Initial value of x

      xfv = lo  ! Final value of x

      ! Calculating the width of the last volume of Part I

      i = iiv

      np   = nx * (j-1) + i

      npw  = np - 1

      a2 = sqrt( (x(np)-x(npw))**2 + (y(np)-y(npw))**2 )

      call get_pl_distribution(iiv, ifv, xiv, xfv, a2, xbs(iiv:ifv)) ! Output: last one

      do i = iiv, ifv

         np   = nx * (j-1) + i

         x(np) = xbs(i)

         y(np) = rb * ( x(np) / lo ) ** lbd

      end do


      ! Part III

      iiv = nx2  ! Initial value of i

      ifv = nx-1 ! Final value of i

      xiv = lo      ! Initial value of x

      xfv = lo + lc ! Final value of x

      ! Calculating the width of the last volume of Part II

      i = iiv

      np   = nx * (j-1) + i

      npw  = np - 1

      a3 = sqrt( (x(np)-x(npw))**2 + (y(np)-y(npw))**2 )

      call get_pl_distribution(iiv, ifv, xiv, xfv, a3, xbs(iiv:ifv)) ! Output: last one

      do i = iiv, ifv

         np   = nx * (j-1) + i

         x(np) = xbs(i)

         y(np) = rb

      end do


      ! ========================================================

      ! Generating north boundary

      j = ny-1

      iiv = 1   ! Initial value of i

      ifv = nx-1 ! Final value of i

      xiv = -lf ! Initial value of x

      xfv =  lr ! Final value of x

      xbn(iiv) = xiv

      do i = iiv+1, ifv

         xbn(i) = xbn(i-1) + (lf+lr) / lr * ( xbs(i)-xbs(i-1))

      end do

      do i = iiv, ifv

         np   = nx * (j-1) + i

         x(np) = xbn(i)

         y(np) = lb * sqrt( 1.d0 - ( (x(np)-lr) / ( lf + lr ) ) ** 2.d0 )

      end do

      ! ========================================================

   contains

      !> \brief Calculates the node distribution based on a power law distribution
      subroutine get_pl_distribution(iiv, ifv, xiv, xfv, a1, xb) ! Output: last one
         implicit none
         integer, intent(in)  :: iiv !< Initial value of i
         integer, intent(in)  :: ifv !< Final value of i
         real(8), intent(in)  :: xiv !< Initial value of x
         real(8), intent(in)  :: xfv !< Final value of x
         real(8), intent(in)  :: a1  !< Width of the first volume
         real(8), intent(out) :: xb(iiv:ifv) !< Distribution of points

         ! Inner variables

         integer :: i
         real(8) :: ak

         ak = log( (xfv-xiv) / a1 ) / log( dble(ifv-iiv) )

         xb(iiv) = xiv

         do i = iiv+1, ifv-1

            xb(i) = xiv + (xfv-xiv) * ( dble(i-iiv) / dble(ifv-iiv) ) ** ak

         end do

         xb(ifv) = xfv

      end subroutine get_pl_distribution

   end subroutine get_grid_boundary_g08

   !> \brief Calculates the coordinates (x,y) of the north and south boundaries.
   !! The north boundary line is an elipse with semi-axis lf+lr and lb.
   !! The south boundary line is an power law ogive of length lo and radius rb
   !! followed by a cylinder of radius rb and length lc.
   !! The distribution of points over the south boundary is separated in three
   !! parts:
   !! Part   I: power law
   !! Part  II: geometric progression
   !! Part III: geometric progression
   !! The distribution of points over the north boundary is proportional to that
   !! in the south boundary.
   subroutine get_grid_boundary_g09(nx, ny, unt, fgeom, lr, rb, x, y) ! Output: last two
      implicit none
      integer, intent(in) :: nx  !< Number of volumes in the csi direction (real+fictitious)
      integer, intent(in) :: ny  !< Number of volumes in the eta direction (real+fictitious)
      integer, intent(in) :: unt !< Unit where the input parameters will be printed
      character(len=*), intent(in) :: fgeom !< File of the geometric parameters

      real(8), intent(out) :: lr  !< Length of the body (m)
      real(8), intent(out) :: rb  !< Base radius of the body (m)
      real(8), dimension(nx*ny), intent(out) :: x !< Coord. x at the northeast corner of the volume P (m)
      real(8), dimension(nx*ny), intent(out) :: y !< Coord. y at the northeast corner of the volume P (m)


      ! Inner variables
      integer :: i, j, np, npw ! Dummy indexes

      integer :: kgb ! Kind of geometry of the body
      integer :: nx1 ! Number of points in the first part of the grid along the eta lines
      integer :: nx2 ! Number of points in the first and second parts of the grid along the eta lines
      integer :: iiv ! Initial value of i
      integer :: ifv ! Final value of i

      real(8) :: a1  ! Width of the first volume of Part I
      real(8) :: a2  ! Width of the first volume of Part II
      real(8) :: a3  ! Width of the first volume of Part III

      real(8) :: xiv ! Initial value of x
      real(8) :: xfv ! Final value of x

      real(8) :: lbd ! Exponent of the power law ogive

      real(8) :: lf  ! Frontal length (m)
      real(8) :: lo  ! Ogive length (m)
      real(8) :: lb  ! Base length (m)
      real(8) :: lc  ! Cylinder length (m)
      real(8) :: lp1 ! Length of the first part of the grid (m)
      real(8) :: aks ! Exponent for the distribution of ksi lines over the south boundary
      real(8) :: fp1 ! Fraction of points in the first part of the grid
      real(8) :: fp2 ! Fraction of points in the second part of the grid


      real(8), dimension(nx) :: xbs ! Distribution of points over the south boundary
      real(8), dimension(nx) :: xbn ! Distribution of points over the north boundary


      ! Reading parameters

      open(13,file=fgeom)

      read(13,*) kgb  ! Kind of geometry of the body
      read(13,*) lf   ! Frontal length (m)
      read(13,*) lb   ! Base length (m)
      read(13,*) lo   ! Length of the ogive (m)
      read(13,*) lc   ! Length of the cylinder (m)
      read(13,*) rb   ! Base radius of the body (m)
      read(13,*) lp1  ! Length of the first part ( 0 < lp1 < lo ) (m)
      read(13,*) lbd  ! Exponent of the power law ogive
      read(13,*) aks  ! Exponent of the power law for the south boundary
      read(13,*) fp1  ! Fraction of volumes in the first part
      read(13,*) fp2  ! Fraction of volumes in the second part

      close(13)


      ! Writing parameters

      write(unt,*)
      write(unt,*) " *** Geometry parameters ***"
      write(unt,*)

      write(unt,"(I23,    A)") kgb, " = kgb: Kind of geometry of the body"
      write(unt,"(ES23.16,A)")  lf, " =  lf: Frontal length (m)"
      write(unt,"(ES23.16,A)")  lb, " =  lb: Base length (m)"
      write(unt,"(ES23.16,A)")  lo, " =  lo: Length of the ogive (m)"
      write(unt,"(ES23.16,A)")  lc, " =  lc: Length of the cylinder (m)"
      write(unt,"(ES23.16,A)")  rb, " =  rb: Base radius of the body (m)"
      write(unt,"(ES23.16,A)") lp1, " = lp1: length of the first part of the grid (m)"
      write(unt,"(ES23.16,A)") lbd, " = lbd: Exponent of the power law ogive"
      write(unt,"(ES23.16,A)") aks, " = aks: Exponent of the power law for the south boundary"
      write(unt,"(ES23.16,A)") fp1, " = fp1: Fraction of volumes in the first part of the grid"
      write(unt,"(ES23.16,A)") fp2, " = fp2: Fraction of volumes in the second part of the grid"


      ! Calculating the number of points in the first part of the grid

      nx1 = int( nx * fp1 )

      ! Calculating the number of points in the first and second parts of the grid

      nx2 = int( nx * fp2 ) + nx1


      ! Calculating the body length

      lr = lo + lc


      ! ========================================================

      ! Generating south boundary

      j = 1

      ! Part I

      iiv = 1   ! Initial value of i

      ifv = nx1 ! Final value of i

      xiv = 0.d0 ! Initial value of x

      xfv = lp1  ! Final value of x

      a1 = (xfv-xiv) / (ifv-iiv) ** aks ! Width of the first volume of the first part of the grid

      call get_pl_distribution(iiv, ifv, xiv, xfv, a1, xbs(iiv:ifv)) ! Output: last one

      do i = iiv, ifv

         np   = nx * (j-1) + i

         x(np) = xbs(i)

         y(np) = rb * ( x(np) / lo ) ** lbd

      end do


      ! Part II

      iiv = nx1  ! Initial value of i

      ifv = nx2  ! Final value of i

      xiv = lp1 ! Initial value of x

      xfv = lo  ! Final value of x

      ! Calculating the width of the last volume of Part I

      i = iiv

      np   = nx * (j-1) + i

      npw  = np - 1

      a2 = sqrt( (x(np)-x(npw))**2 + (y(np)-y(npw))**2 )

      call get_gp_distribution(iiv, ifv, xiv, xfv, a2, xbs(iiv:ifv)) ! Output: last one

      do i = iiv, ifv

         np   = nx * (j-1) + i

         x(np) = xbs(i)

         y(np) = rb * ( x(np) / lo ) ** lbd

      end do


      ! Part III

      iiv = nx2  ! Initial value of i

      ifv = nx-1 ! Final value of i

      xiv = lo      ! Initial value of x

      xfv = lo + lc ! Final value of x

      ! Calculating the width of the last volume of Part II

      i = iiv

      np   = nx * (j-1) + i

      npw  = np - 1

      a3 = sqrt( (x(np)-x(npw))**2 + (y(np)-y(npw))**2 )

      call get_gp_distribution(iiv, ifv, xiv, xfv, a3, xbs(iiv:ifv)) ! Output: last one

      do i = iiv, ifv

         np   = nx * (j-1) + i

         x(np) = xbs(i)

         y(np) = rb

      end do


      ! ========================================================

      ! Generating north boundary

      j = ny-1

      iiv = 1   ! Initial value of i

      ifv = nx-1 ! Final value of i

      xiv = -lf ! Initial value of x

      xfv =  lr ! Final value of x

      xbn(iiv) = xiv

      do i = iiv+1, ifv-1

         xbn(i) = xbn(i-1) + (lf+lr) / lr * ( xbs(i)-xbs(i-1))

      end do

      xbn(ifv) = xfv

      do i = iiv, ifv

         np   = nx * (j-1) + i

         x(np) = xbn(i)

         y(np) = lb * sqrt( 1.d0 - ( (x(np)-lr) / ( lf + lr ) ) ** 2.d0 )

      end do

      ! ========================================================

   contains

      !> \brief Calculates the node distribution based on a power law distribution
      subroutine get_pl_distribution(iiv, ifv, xiv, xfv, a1, xb) ! Output: last one
         implicit none
         integer, intent(in)  :: iiv !< Initial value of i
         integer, intent(in)  :: ifv !< Final value of i
         real(8), intent(in)  :: xiv !< Initial value of x
         real(8), intent(in)  :: xfv !< Final value of x
         real(8), intent(in)  :: a1  !< Width of the first volume
         real(8), intent(out) :: xb(iiv:ifv) !< Distribution of points

         ! Inner variables

         integer :: i
         real(8) :: ak

         ak = log( (xfv-xiv) / a1 ) / log( dble(ifv-iiv) )

         xb(iiv) = xiv

         do i = iiv+1, ifv-1

            xb(i) = xiv + (xfv-xiv) * ( dble(i-iiv) / dble(ifv-iiv) ) ** ak

         end do

         xb(ifv) = xfv

      end subroutine get_pl_distribution


      !> \brief Calculates the node distribution based on a geometric progression distribution
      subroutine get_gp_distribution(iiv, ifv, xiv, xfv, a1, xb) ! Output: last one
         implicit none
         integer, intent(in)  :: iiv !< Initial value of i
         integer, intent(in)  :: ifv !< Final value of i
         real(8), intent(in)  :: xiv !< Initial value of x
         real(8), intent(in)  :: xfv !< Final value of x
         real(8), intent(in)  :: a1  !< Width of the first volume
         real(8), intent(out) :: xb(iiv:ifv) !< Distribution of points

         ! Inner variables

         integer :: i
         integer :: n
         real(8) :: r
         real(8) :: q


         n = ifv - iiv

         r = ( xfv - xiv ) / a1

         ! Calculating q
         call get_gp_ratio(n, r, q)

         xb(iiv) = xiv

         do i = iiv+1, ifv-1

            xb(i) = xb(i-1) + a1 * q ** ( i - iiv - 1 )

         end do

         xb(ifv) = xfv

      end subroutine get_gp_distribution


      !> \brief Calculates the base q of the geometric progression
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

           if ( 0.d0 < f(n, r, qi) * f(n, r, qm) ) then

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

      end subroutine get_gp_ratio

      !> \brief Function used in the subroutine 'get_gp_ratio'
      real(8) function f(n, r, q)
         implicit none
         integer, intent(in)  ::   n !< number of partitions
         real(8), intent(in)  ::   r !< l/a1
         real(8), intent(in)  ::   q !< q

         f = q ** n + r * ( 1.d0 - q ) - 1.d0

      end function f

   end subroutine get_grid_boundary_g09


   !> \brief Calculates the coordinates (x,y) of the north and south boundaries.
   !! The north boundary line is an elipse with semi-axis lf+lr and lb.
   !! The south boundary line is an power law ogive of length lo and radius rb
   !! followed by a cylinder of radius rb and length lc.
   !! The distribution of points over the south boundary is separated in three
   !! parts:
   !! Part   I: power law
   !! Part  II: geometric progression
   !! Part III: geometric progression
   !! The distribution of points over the north boundary is proportional to that
   !! in the south boundary.
   !! Adapted for optimization with DEPP optimizer.
   subroutine get_grid_boundary_g10(nx, ny, unt, fgeom, lr, rb, x, y) ! Output: last two
      implicit none
      integer, intent(in) :: nx  !< Number of volumes in the csi direction (real+fictitious)
      integer, intent(in) :: ny  !< Number of volumes in the eta direction (real+fictitious)
      integer, intent(in) :: unt !< Unit where the input parameters will be printed
      character(len=*), intent(in) :: fgeom !< File of the geometric parameters

      real(8), intent(out) :: lr  !< Length of the body (m)
      real(8), intent(out) :: rb  !< Base radius of the body (m)
      real(8), dimension(nx*ny), intent(out) :: x !< Coord. x at the northeast corner of the volume P (m)
      real(8), dimension(nx*ny), intent(out) :: y !< Coord. y at the northeast corner of the volume P (m)


      ! Inner variables
      integer :: i, j, np, npw ! Dummy indexes

      integer :: kgb ! Kind of geometry of the body
      integer :: nx1 ! Number of points in the first part of the grid along the eta lines
      integer :: nx2 ! Number of points in the first and second parts of the grid along the eta lines
      integer :: iiv ! Initial value of i
      integer :: ifv ! Final value of i

      real(8) :: a1  ! Width of the first volume of Part I
      real(8) :: a2  ! Width of the first volume of Part II
      real(8) :: a3  ! Width of the first volume of Part III

      real(8) :: xiv ! Initial value of x
      real(8) :: xfv ! Final value of x

      real(8) :: lbd ! Exponent of the power law ogive

      real(8) :: lf  ! Frontal length (m)
      real(8) :: lo  ! Ogive length (m)
      real(8) :: lb  ! Base length (m)
      real(8) :: lc  ! Cylinder length (m)
      real(8) :: lp1 ! Length of the first part of the grid (m)
      real(8) :: aks ! Exponent for the distribution of ksi lines over the south boundary
      real(8) :: fp1 ! Fraction of points in the first part of the grid
      real(8) :: fp2 ! Fraction of points in the second part of the grid


      real(8), dimension(nx) :: xbs ! Distribution of points over the south boundary
      real(8), dimension(nx) :: xbn ! Distribution of points over the north boundary


      ! DEPP variables
      integer, parameter :: nu = 1 ! number of unknowns
      integer            :: ind    ! number of the individual
      integer,       dimension(nu) :: xopt   ! Optimization checker
      character(10), dimension(nu) :: xname  ! Name of parameters
      real(8),       dimension(nu) :: xvalue ! parameters
      character(200)               :: sname  ! simulation name



      ! Reading parameters

      open(13,file=fgeom)

      read(13,*) kgb  ! Kind of geometry of the body
      read(13,*) lf   ! Frontal length (m)
      read(13,*) lb   ! Base length (m)
      read(13,*) lo   ! Length of the ogive (m)
      read(13,*) lc   ! Length of the cylinder (m)
      read(13,*) rb   ! Base radius of the body (m)
      read(13,*) lp1  ! Length of the first part ( 0 < lp1 < lo ) (m)
      read(13,*) xname(1) ! Parameter name
      read(13,*) xopt(1)  ! Will this parameter be optimized? ( 0 = no, 1 = yes )
      read(13,*) xvalue(1)! Parameter value
      read(13,*) aks  ! Exponent of the power law for the south boundary
      read(13,*) fp1  ! Fraction of volumes in the first part
      read(13,*) fp2  ! Fraction of volumes in the second part

      close(13)



      ! Reads the parameters from DEPP
      call depp_get_parameters(nu, xopt, xname, xvalue, ind, sname)



      lbd = xvalue(1)


      ! Writing parameters

      write(unt,*)
      write(unt,*) " *** Geometry parameters ***"
      write(unt,*)

      write(unt,"(I23,    A)") kgb, " = kgb: Kind of geometry of the body"
      write(unt,"(ES23.16,A)")  lf, " =  lf: Frontal length (m)"
      write(unt,"(ES23.16,A)")  lb, " =  lb: Base length (m)"
      write(unt,"(ES23.16,A)")  lo, " =  lo: Length of the ogive (m)"
      write(unt,"(ES23.16,A)")  lc, " =  lc: Length of the cylinder (m)"
      write(unt,"(ES23.16,A)")  rb, " =  rb: Base radius of the body (m)"
      write(unt,"(ES23.16,A)") lp1, " = lp1: length of the first part of the grid (m)"
      write(unt,"(ES23.16,A)") lbd, " = lbd: Exponent of the power law ogive"
      write(unt,"(ES23.16,A)") aks, " = aks: Exponent of the power law for the south boundary"
      write(unt,"(ES23.16,A)") fp1, " = fp1: Fraction of volumes in the first part of the grid"
      write(unt,"(ES23.16,A)") fp2, " = fp2: Fraction of volumes in the second part of the grid"


      ! Calculating the number of points in the first part of the grid

      nx1 = int( nx * fp1 )

      ! Calculating the number of points in the first and second parts of the grid

      nx2 = int( nx * fp2 ) + nx1


      ! Calculating the body length

      lr = lo + lc


      ! ========================================================

      ! Generating south boundary

      j = 1

      ! Part I

      iiv = 1   ! Initial value of i

      ifv = nx1 ! Final value of i

      xiv = 0.d0 ! Initial value of x

      xfv = lp1  ! Final value of x

      a1 = (xfv-xiv) / (ifv-iiv) ** aks ! Width of the first volume of the first part of the grid

      call get_pl_distribution(iiv, ifv, xiv, xfv, a1, xbs(iiv:ifv)) ! Output: last one

      do i = iiv, ifv

         np   = nx * (j-1) + i

         x(np) = xbs(i)

         y(np) = rb * ( x(np) / lo ) ** lbd

      end do


      ! Part II

      iiv = nx1  ! Initial value of i

      ifv = nx2  ! Final value of i

      xiv = lp1 ! Initial value of x

      xfv = lo  ! Final value of x

      ! Calculating the width of the last volume of Part I

      i = iiv

      np   = nx * (j-1) + i

      npw  = np - 1

      a2 = sqrt( (x(np)-x(npw))**2 + (y(np)-y(npw))**2 )

      call get_gp_distribution(iiv, ifv, xiv, xfv, a2, xbs(iiv:ifv)) ! Output: last one

      do i = iiv, ifv

         np   = nx * (j-1) + i

         x(np) = xbs(i)

         y(np) = rb * ( x(np) / lo ) ** lbd

      end do


      ! Part III

      iiv = nx2  ! Initial value of i

      ifv = nx-1 ! Final value of i

      xiv = lo      ! Initial value of x

      xfv = lo + lc ! Final value of x

      ! Calculating the width of the last volume of Part II

      i = iiv

      np   = nx * (j-1) + i

      npw  = np - 1

      a3 = sqrt( (x(np)-x(npw))**2 + (y(np)-y(npw))**2 )

      call get_gp_distribution(iiv, ifv, xiv, xfv, a3, xbs(iiv:ifv)) ! Output: last one

      do i = iiv, ifv

         np   = nx * (j-1) + i

         x(np) = xbs(i)

         y(np) = rb

      end do


      ! ========================================================

      ! Generating north boundary

      j = ny-1

      iiv = 1   ! Initial value of i

      ifv = nx-1 ! Final value of i

      xiv = -lf ! Initial value of x

      xfv =  lr ! Final value of x

      xbn(iiv) = xiv

      do i = iiv+1, ifv-1

         xbn(i) = xbn(i-1) + (lf+lr) / lr * ( xbs(i)-xbs(i-1))

      end do

      xbn(ifv) = xfv

      do i = iiv, ifv

         np   = nx * (j-1) + i

         x(np) = xbn(i)

         y(np) = lb * sqrt( 1.d0 - ( (x(np)-lr) / ( lf + lr ) ) ** 2.d0 )

      end do

      ! ========================================================

   contains

      !> \brief Calculates the node distribution based on a power law distribution
      subroutine get_pl_distribution(iiv, ifv, xiv, xfv, a1, xb) ! Output: last one
         implicit none
         integer, intent(in)  :: iiv !< Initial value of i
         integer, intent(in)  :: ifv !< Final value of i
         real(8), intent(in)  :: xiv !< Initial value of x
         real(8), intent(in)  :: xfv !< Final value of x
         real(8), intent(in)  :: a1  !< Width of the first volume
         real(8), intent(out) :: xb(iiv:ifv) !< Distribution of points

         ! Inner variables

         integer :: i
         real(8) :: ak

         ak = log( (xfv-xiv) / a1 ) / log( dble(ifv-iiv) )

         xb(iiv) = xiv

         do i = iiv+1, ifv-1

            xb(i) = xiv + (xfv-xiv) * ( dble(i-iiv) / dble(ifv-iiv) ) ** ak

         end do

         xb(ifv) = xfv

      end subroutine get_pl_distribution


      !> \brief Calculates the node distribution based on a geometric progression distribution
      subroutine get_gp_distribution(iiv, ifv, xiv, xfv, a1, xb) ! Output: last one
         implicit none
         integer, intent(in)  :: iiv !< Initial value of i
         integer, intent(in)  :: ifv !< Final value of i
         real(8), intent(in)  :: xiv !< Initial value of x
         real(8), intent(in)  :: xfv !< Final value of x
         real(8), intent(in)  :: a1  !< Width of the first volume
         real(8), intent(out) :: xb(iiv:ifv) !< Distribution of points

         ! Inner variables

         integer :: i
         integer :: n
         real(8) :: r
         real(8) :: q


         n = ifv - iiv

         r = ( xfv - xiv ) / a1

         ! Calculating q
         call get_gp_ratio(n, r, q)

         xb(iiv) = xiv

         do i = iiv+1, ifv-1

            xb(i) = xb(i-1) + a1 * q ** ( i - iiv - 1 )

         end do

         xb(ifv) = xfv

      end subroutine get_gp_distribution


      !> \brief Calculates the base q of the geometric progression
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

           if ( 0.d0 < f(n, r, qi) * f(n, r, qm) ) then

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

      end subroutine get_gp_ratio

      !> \brief Function used in the subroutine 'get_gp_ratio'
      real(8) function f(n, r, q)
         implicit none
         integer, intent(in)  ::   n !< number of partitions
         real(8), intent(in)  ::   r !< l/a1
         real(8), intent(in)  ::   q !< q

         f = q ** n + r * ( 1.d0 - q ) - 1.d0

      end function f

   end subroutine get_grid_boundary_g10


   !> \brief Calculates the coordinates (x,y) of the north and south boundaries.
   !! The north boundary line is an elipse with semi-axis lf+lr and lb.
   !! The south boundary line is an bluff nose followed by a power law ogive of
   !! length lo and radius rb followed by a cylinder of radius rb and length lc.
   !! The distribution of points over the south boundary is separated in four
   !! parts:
   !! Part   I: geometric progression (bluff nose)
   !! Part  II: geometric progression (power law)
   !! Part III: geometric progression (power law)
   !! Part  IV: geometric progression (cylinder)
   !! The distribution of points over the north boundary is proportional to that
   !! in the south boundary.
   subroutine get_grid_boundary_g11(nx, ny, unt, fgeom, lr, rb, x, y) ! Output: last two
      implicit none
      integer, intent(in) :: nx  !< Number of volumes in the csi direction (real+fictitious)
      integer, intent(in) :: ny  !< Number of volumes in the eta direction (real+fictitious)
      integer, intent(in) :: unt !< Unit where the input parameters will be printed
      character(len=*), intent(in) :: fgeom !< File of the geometric parameters

      real(8), intent(out) :: lr  !< Length of the body (m)
      real(8), intent(out) :: rb  !< Base radius of the body (m)
      real(8), dimension(nx*ny), intent(out) :: x !< Coord. x at the northeast corner of the volume P (m)
      real(8), dimension(nx*ny), intent(out) :: y !< Coord. y at the northeast corner of the volume P (m)


      ! Inner variables
      integer :: i, j, np ! Dummy indexes

      integer :: kgb ! Kind of geometry of the body
      integer :: nx1 ! Number of points in the first part of the grid along the eta lines
      integer :: nx2 ! Number of points in the first and second parts of the grid along the eta lines
      integer :: nx3 ! Number of points in the first to third parts of the grid along the eta lines
      integer :: iiv ! Initial value of i
      integer :: ifv ! Final value of i

      real(8) :: a1  ! Width of the first volume of Part I
      real(8) :: a2  ! Width of the first volume of Part II
      real(8) :: a3  ! Width of the first volume of Part III
      real(8) :: a4  ! Width of the first volume of Part IV

      real(8) :: xiv ! Initial value of x
      real(8) :: xfv ! Final value of x
      real(8) :: yiv ! Initial value of y
      real(8) :: yfv ! Final value of y

      real(8) :: lbd ! Exponent of the power law ogive

      real(8) :: lf  ! Frontal length (m)
      real(8) :: lo  ! Ogive length (m)
      real(8) :: lb  ! Base length (m)
      real(8) :: lc  ! Cylinder length (m)
      real(8) :: lo1 ! Length of the first part of the ogive (m)
      real(8) :: rf  ! Radius of the nose tip (m)
      real(8) :: h   ! Ratio rf/rb
      real(8) :: aks ! Exponent of the power law for the calculation of the width of the first volume
      real(8) :: flo ! Fraction of lo
      real(8) :: fp1 ! Fraction of points in the first part of the grid
      real(8) :: fp2 ! Fraction of points in the second part of the grid
      real(8) :: fp3 ! Fraction of points in the third part of the grid


      real(8), dimension(nx) :: xbs ! Distribution of points over the south boundary
      real(8), dimension(nx) :: ybs ! Distribution of points over the south boundary
      real(8), dimension(nx) :: xbn ! Distribution of points over the north boundary
      real(8), dimension(nx) :: ybn ! Distribution of points over the north boundary


      ! Reading parameters

      open(13,file=fgeom)

      read(13,*) kgb  ! Kind of geometry of the body
      read(13,*) lf   ! Frontal length (m)
      read(13,*) lb   ! Base length (m)
      read(13,*) lo   ! Length of the ogive (m)
      read(13,*) lc   ! Length of the cylinder (m)
      read(13,*) rb   ! Base radius of the body (m)
      read(13,*) h    ! Ratio rf/rb
      read(13,*) lbd  ! Exponent of the power law ogive
      read(13,*) aks  ! Exponent of the power law for the south boundary
      read(13,*) flo  ! Fraction of lo
      read(13,*) fp1  ! Fraction of volumes in the first part
      read(13,*) fp2  ! Fraction of volumes in the second part
      read(13,*) fp3  ! Fraction of volumes in the third part

      close(13)


      ! Writing parameters

      write(unt,*)
      write(unt,*) " *** Geometry parameters ***"
      write(unt,*)

      write(unt,"(I23,    A)") kgb, " = kgb: Kind of geometry of the body"
      write(unt,"(ES23.16,A)")  lf, " =  lf: Frontal length (m)"
      write(unt,"(ES23.16,A)")  lb, " =  lb: Base length (m)"
      write(unt,"(ES23.16,A)")  lo, " =  lo: Length of the ogive (m)"
      write(unt,"(ES23.16,A)")  lc, " =  lc: Length of the cylinder (m)"
      write(unt,"(ES23.16,A)")  rb, " =  rb: Base radius of the body (m)"
      write(unt,"(ES23.16,A)")   h, " =   h: Ratio rf/rb (frontal to base radius ratio)"
      write(unt,"(ES23.16,A)") lbd, " = lbd: Exponent of the power law ogive"
      write(unt,"(ES23.16,A)") aks, " = aks: Exponent of the power law for the south boundary"
      write(unt,"(ES23.16,A)") flo, " = flo: Fraction of lo"
      write(unt,"(ES23.16,A)") fp1, " = fp1: Fraction of volumes in the first part of the grid"
      write(unt,"(ES23.16,A)") fp2, " = fp2: Fraction of volumes in the second part of the grid"
      write(unt,"(ES23.16,A)") fp3, " = fp3: Fraction of volumes in the third part of the grid"

      ! Calculating the number of points in the first part of the grid

      nx1 = int( nx * fp1 )

      ! Calculating the number of points in the first and second parts of the grid

      nx2 = int( nx * fp2 ) + nx1

      ! Calculating the number of points in the 1st, 2nd and 3rd part of the grid

      nx3 = int( fp3 * nx ) + nx2

      ! Calculating frontal radius

      rf = h * rb

      ! Calculating the fraction of lo used in the 2nd part of the grid

      lo1 = flo * lo

      ! Calculating the body length

      lr = lo + lc


      ! ========================================================

      ! Generating south boundary

      j = 1

      ! Part I

      iiv = 1   ! Initial value of i

      ifv = nx1 ! Final value of i

      yiv = 0.d0 ! Initial value of y

      yfv = rf  ! Final value of y

      a1 = (yfv-yiv) * ( 1.d0 - (dble(ifv-iiv-1)/ dble(ifv-iiv)) ** aks ) ! Width of the first volume of the first part of the grid

      call get_gp_distribution(iiv, ifv, yiv, yfv, a1, ybs(iiv:ifv)) ! Output: last one

      do i = iiv, ifv

         np   = nx * (j-1) + i

         xbs(i) = 0.d0

         x(np) = xbs(i)

         y(np) = ybs(i)

      end do


      ! Part II

      iiv = nx1  ! Initial value of i

      ifv = nx2  ! Final value of i

      xiv = 0.d0 ! Initial value of x

      xfv = lo1  ! Final value of x

      ! Calculating the width of the first volume of Part II

      i = iiv

      a2 = sqrt( (xbs(i)-xbs(i-1))**2 + (ybs(i)-ybs(i-1))**2 )

      call get_a2(h, rb, lo, lbd, a2) ! InOut: last one

      call get_gp_distribution(iiv, ifv, xiv, xfv, a2, xbs(iiv:ifv)) ! Output: last one
      !call get_pl_distribution(iiv, ifv, xiv, xfv, a2, xbs(iiv:ifv)) ! Output: last one

      do i = iiv, ifv

         np   = nx * (j-1) + i

         ybs(i) = rb * ( h + (1.d0-h) * ( xbs(i) / lo ) ** lbd )

         x(np) = xbs(i)

         y(np) = ybs(i)

      end do

      ! Part III

      iiv = nx2  ! Initial value of i

      ifv = nx3 ! Final value of i

      xiv = lo1 ! Initial value of x

      xfv = lo  ! Final value of x

      ! Calculating the width of the first volume of Part III

      i = iiv

      a3 = sqrt( (xbs(i)-xbs(i-1))**2 + (ybs(i)-ybs(i-1))**2 )

      call get_gp_distribution(iiv, ifv, xiv, xfv, a3, xbs(iiv:ifv)) ! Output: last one

      do i = iiv, ifv

         np   = nx * (j-1) + i

         ybs(i) = rb * ( h + (1.d0-h) * ( xbs(i) / lo ) ** lbd )

         x(np) = xbs(i)

         y(np) = ybs(i)

      end do

      ! Part IV

      iiv = nx3  ! Initial value of i

      ifv = nx-1 ! Final value of i

      xiv = lo ! Initial value of x

      xfv = lo+lc  ! Final value of x

      ! Calculating the width of the first volume of Part IV

      i = iiv

      a4 = sqrt( (xbs(i)-xbs(i-1))**2 + (ybs(i)-ybs(i-1))**2 )

      call get_gp_distribution(iiv, ifv, xiv, xfv, a4, xbs(iiv:ifv)) ! Output: last one

      do i = iiv, ifv

         np   = nx * (j-1) + i

         ybs(i) = rb

         x(np) = xbs(i)

         y(np) = ybs(i)

      end do

      ! ========================================================

      ! Generating north boundary

      call get_optimized_noth_boundary()

      ! ========================================================

   contains

      !> \brief Calculates the skewness of the grid
      subroutine get_skewness(tht)
         implicit none
         real(8), intent(out) :: tht

         real(8) :: tht12, tht23, tht34, tht41, raux
         real(8) :: n2r1, n2r2, n2r3, n2r4

         real(8) :: r1(2)
         real(8) :: r2(2)
         real(8) :: r3(2)
         real(8) :: r4(2)

         tht = 0.d0

         do i = 2, nx-1

            r1(1) = xbs(i)-xbs(i-1)
            r1(2) = ybs(i)-ybs(i-1)

            r2(1) = xbn(i)-xbs(i)
            r2(2) = ybn(i)-ybs(i)

            r3(1) = xbn(i-1)-xbn(i)
            r3(2) = ybn(i-1)-ybn(i)

            r4(1) = xbs(i-1)-xbn(i-1)
            r4(2) = ybs(i-1)-ybn(i-1)

            n2r1 = sqrt( r1(1)**2 + r1(2)**2 )
            n2r2 = sqrt( r2(1)**2 + r2(2)**2 )
            n2r3 = sqrt( r3(1)**2 + r3(2)**2 )
            n2r4 = sqrt( r4(1)**2 + r4(2)**2 )

            tht12 = abs(acos(-1.d0)/2.d0-acos( dot_product(r1,r2) / (n2r1*n2r2)))
            tht23 = abs(acos(-1.d0)/2.d0-acos( dot_product(r2,r3) / (n2r2*n2r3)))
            tht34 = abs(acos(-1.d0)/2.d0-acos( dot_product(r3,r4) / (n2r3*n2r4)))
            tht41 = abs(acos(-1.d0)/2.d0-acos( dot_product(r4,r1) / (n2r4*n2r1)))

            !raux = max(tht12, tht23, tht34, tht41)
            raux = max(tht12, tht41)

            if (raux>tht) tht = raux

         end do

      end subroutine get_skewness


      !> \brief Calculates the width a2 of the first volume of Part II
      subroutine get_a2(h, rb, lo, lbd, a2)
         implicit none
         real(8), intent(in)    :: h   !< Ratio rf/rb
         real(8), intent(in)    :: rb  !< Base radius (m)
         real(8), intent(in)    :: lo  !< Ogive length (m)
         real(8), intent(in)    :: lbd !< Exponent of the power law ogive
         real(8), intent(inout) :: a2  !< Width of the first volume of Part II

         ! Parameters

         integer, parameter :: nit = 5000
         real(8), parameter :: tol = 1.d-10


         ! Inner variables

         integer :: k

         real(8) :: dx, dy, raux

         raux = ( (1.d0-h) * rb / lo ** lbd ) ** 2

         dx = 0.d0

         do k = 1, nit

            dx = ( a2**2 / ( raux + dx ** (2.d0-2.d0*lbd) ) ) ** (1.d0/(2.d0*lbd))

            dy = sqrt(raux)*dx**lbd

            !print*, k, a2**2, dy**2+dx**2, sqrt(raux)*dx**lbd, dx, (sqrt(dx**2 + dy**2) - a2) / a2

            if ( abs( sqrt(dx**2 + dy**2) - a2) / a2 < tol ) exit

         end do

         if ( k >= nit ) then

            write(*,*) "Error: get_a2 not converged. Stopping..."

            stop

         end if

         a2 = dx

      end subroutine get_a2

      !> \brief Adjust one parameter in order to maximize the orthogonality of
      !! the grid
      subroutine get_optimized_noth_boundary()
         implicit none

         ! Parameters

         integer, parameter :: nit = 1000

         ! Inner variables

         integer :: k

         real(8) :: rfn
         real(8) :: rfni
         real(8) :: rfnf
         real(8) :: rfnc
         real(8) :: tht
         real(8) :: thtc

         thtc = acos(-1.d0)/2.d0

         rfni = rf

         rfnf = rf + lf * 1.8d0

         do k = 1, nit

            rfn = (rfnf-rfni) * dble(k-1)/dble(nit-1) + rfni

            call get_north_boundary(rfn)

            call get_skewness(tht)

            if ( tht < thtc ) then

               thtc = tht

               rfnc = rfn

            end if

            !print*, tht, thtc, rfn, rfnc

         end do

         ! Optimized grid

         call get_north_boundary(rfnc)

      end subroutine get_optimized_noth_boundary


      !> \brief Calculates the distribution of points over the north boundary
      subroutine get_north_boundary(rfn)
         implicit none
         real(8), intent(in) :: rfn !< Frontal radius of the north boundary (m)

         j = ny-1


         ! Part I

         iiv = 1   ! Initial value of i

         ifv = nx1 ! Final value of i

         yiv = 0.d0 ! Initial value of y

         yfv = rfn  ! Final value of y

         ybn(iiv) = yiv

         do i = iiv+1, ifv-1

            ybn(i) = ybn(i-1) + rfn / rf * ( ybs(i)-ybs(i-1))

         end do

         ybn(ifv) = yfv

         do i = iiv, ifv

            np   = nx * (j-1) + i

            xbn(i) = lo - ( lo + lf ) * sqrt( 1.d0-(ybn(i)/lb)**2 )

            x(np) = xbn(i)

            y(np) = ybn(i)

         end do

         ! Part II to Part IV

         iiv = nx1   ! Initial value of i

         ifv = nx-1 ! Final value of i

         xiv = xbn(iiv) ! Initial value of x

         xfv = lo+lc  ! Final value of x

         xbn(iiv) = xiv

         do i = iiv+1, ifv-1

            xbn(i) = xbn(i-1) + (xfv-xiv) / (lo+lc) * ( xbs(i)-xbs(i-1))

         end do

         xbn(ifv) = xfv

         do i = iiv, ifv

            np   = nx * (j-1) + i

            if ( xbn(i) < lo ) then

               ybn(i) = lb * sqrt(1.d0-((xbn(i)-lo)/(lo+lf))**2)

            else

               ybn(i) = lb

            end if

            x(np) = xbn(i)

            y(np) = ybn(i)

         end do

      end subroutine get_north_boundary



      !> \brief Calculates the node distribution based on a power law distribution
      subroutine get_pl_distribution(iiv, ifv, xiv, xfv, a1, xb) ! Output: last one
         implicit none
         integer, intent(in)  :: iiv !< Initial value of i
         integer, intent(in)  :: ifv !< Final value of i
         real(8), intent(in)  :: xiv !< Initial value of x
         real(8), intent(in)  :: xfv !< Final value of x
         real(8), intent(in)  :: a1  !< Width of the first volume
         real(8), intent(out) :: xb(iiv:ifv) !< Distribution of points

         ! Inner variables

         integer :: i
         real(8) :: ak

         ak = log( (xfv-xiv) / a1 ) / log( dble(ifv-iiv) )

         xb(iiv) = xiv

         do i = iiv+1, ifv-1

            xb(i) = xiv + (xfv-xiv) * ( dble(i-iiv) / dble(ifv-iiv) ) ** ak

         end do

         xb(ifv) = xfv

      end subroutine get_pl_distribution


      !> \brief Calculates the node distribution based on a geometric progression distribution
      subroutine get_gp_distribution(iiv, ifv, xiv, xfv, a1, xb) ! Output: last one
         implicit none
         integer, intent(in)  :: iiv !< Initial value of i
         integer, intent(in)  :: ifv !< Final value of i
         real(8), intent(in)  :: xiv !< Initial value of x
         real(8), intent(in)  :: xfv !< Final value of x
         real(8), intent(in)  :: a1  !< Width of the first volume
         real(8), intent(out) :: xb(iiv:ifv) !< Distribution of points

         ! Inner variables

         integer :: i
         integer :: n
         real(8) :: r
         real(8) :: q


         n = ifv - iiv

         r = ( xfv - xiv ) / a1

         ! Calculating q
         call get_gp_ratio(n, r, q)

         xb(iiv) = xiv

         do i = iiv+1, ifv-1

            xb(i) = xb(i-1) + a1 * q ** ( i - iiv - 1 )

         end do

         xb(ifv) = xfv

      end subroutine get_gp_distribution


      !> \brief Calculates the base q of the geometric progression
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

           if ( 0.d0 < f(n, r, qi) * f(n, r, qm) ) then

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

      end subroutine get_gp_ratio

      !> \brief Function used in the subroutine 'get_gp_ratio'
      real(8) function f(n, r, q)
         implicit none
         integer, intent(in)  ::   n !< number of partitions
         real(8), intent(in)  ::   r !< l/a1
         real(8), intent(in)  ::   q !< q

         f = q ** n + r * ( 1.d0 - q ) - 1.d0

      end function f

   end subroutine get_grid_boundary_g11



   !> \brief Calculates the coordinates (x,y) of the north and south boundaries.
   !! The north boundary line is an elipse with semi-axis lf+lr and lb.
   !! The south boundary line is an bluff nose followed by a power law ogive of
   !! length lo and radius rb followed by a cylinder of radius rb and length lc.
   !! The distribution of points over the south boundary is separated in four
   !! parts:
   !! Part   I: geometric progression (bluff nose)
   !! Part  II: geometric progression (power law)
   !! Part III: geometric progression (power law)
   !! Part  IV: geometric progression (cylinder)
   !! The distribution of points over the north boundary is  based on the normal
   !! vectors from the south boundary. The directions of the normal vectors
   !! may be changed by a smoothing procedure.
   !! Adapted to DEPP optimizer.
   subroutine get_grid_boundary_g12(nx, ny, unt, fgeom, lr, rb, x, y) ! Output: last two
      implicit none
      integer, intent(in) :: nx  !< Number of volumes in the csi direction (real+fictitious)
      integer, intent(in) :: ny  !< Number of volumes in the eta direction (real+fictitious)
      integer, intent(in) :: unt !< Unit where the input parameters will be printed
      character(len=*), intent(in) :: fgeom !< File of the geometric parameters

      real(8), intent(out) :: lr  !< Length of the body (m)
      real(8), intent(out) :: rb  !< Base radius of the body (m)
      real(8), dimension(nx*ny), intent(out) :: x !< Coord. x at the northeast corner of the volume P (m)
      real(8), dimension(nx*ny), intent(out) :: y !< Coord. y at the northeast corner of the volume P (m)

      ! Parameters

      integer, parameter :: nr = 10000 ! Number of partitions for the representing curve (xr,yr)


      ! Inner variables
      integer :: i, j, np ! Dummy indexes

      integer :: kgb ! Kind of geometry of the body
      integer :: nx1 ! Number of points in the first part of the grid along the eta lines
      integer :: nx2 ! Number of points in the first and second parts of the grid along the eta lines
      integer :: nx3 ! Number of points in the first to third parts of the grid along the eta lines
      integer :: iiv ! Initial value of i
      integer :: ifv ! Final value of i

      real(8) :: a1  ! Width of the first volume of Part I
      real(8) :: a2  ! Width of the first volume of Part II
      real(8) :: a3  ! Width of the first volume of Part III
      real(8) :: a4  ! Width of the first volume of Part IV

      real(8) :: xiv ! Initial value of x
      real(8) :: xfv ! Final value of x
      real(8) :: yiv ! Initial value of y
      real(8) :: yfv ! Final value of y

      real(8) :: lbd ! Exponent of the power law ogive

      real(8) :: lo  ! Ogive length (m)
      real(8) :: w   ! Distance between south and north boundaries (m)
      real(8) :: fs  ! Factor of smoothing (0<=fs) (necessary for non-smooth surfaces)
      real(8) :: lc  ! Cylinder length (m)
      real(8) :: lo1 ! Length of the first part of the ogive (m)
      real(8) :: rf  ! Radius of the nose tip (m)
      real(8) :: h   ! Ratio rf/rb
      real(8) :: aks ! Exponent of the power law for the calculation of the width of the first volume
      real(8) :: flo ! Fraction of lo
      real(8) :: fp1 ! Fraction of points in the first part of the grid
      real(8) :: fp2 ! Fraction of points in the second part of the grid
      real(8) :: fp3 ! Fraction of points in the third part of the grid


      real(8), dimension(nx) :: xbs ! Distribution of points over the south boundary
      real(8), dimension(nx) :: ybs ! Distribution of points over the south boundary

      real(8), dimension(0:nr) :: xr ! Distribution of points over the representing curve
      real(8), dimension(0:nr) :: yr ! Distribution of points over the representing curve


      ! DEPP variables
      integer, parameter :: nu = 2 ! number of unknowns
      integer            :: ind    ! number of the individual
      integer,       dimension(nu) :: xopt   ! Optimization checker
      character(10), dimension(nu) :: xname  ! Name of parameters
      real(8),       dimension(nu) :: xvalue ! parameters
      character(200)               :: sname  ! simulation name


      ! Reading parameters

      open(13,file=fgeom)

      read(13,*) kgb  ! Kind of geometry of the body
      read(13,*) w    ! Distance between south and north boundaries (m)
      read(13,*) fs   ! Factor of smoothing (0<=fs) (necessary for non-smooth surfaces)
      read(13,*) lo   ! Length of the ogive (m)
      read(13,*) lc   ! Length of the cylinder (m)
      read(13,*) rb   ! Base radius of the body (m)
      read(13,*) xname(1) ! Parameter name
      read(13,*) xopt(1)  ! Will this parameter be optimized? ( 0 = no, 1 = yes )
      read(13,*) xvalue(1)! Parameter value
      read(13,*) xname(2) ! Parameter name
      read(13,*) xopt(2)  ! Will this parameter be optimized? ( 0 = no, 1 = yes )
      read(13,*) xvalue(2)! Parameter value
      read(13,*) aks  ! Exponent of the power law for the south boundary
      read(13,*) flo  ! Fraction of lo
      read(13,*) fp1  ! Fraction of volumes in the first part
      read(13,*) fp2  ! Fraction of volumes in the second part
      read(13,*) fp3  ! Fraction of volumes in the third part

      close(13)

      ! Reads the parameters from DEPP
      call depp_get_parameters(nu, xopt, xname, xvalue, ind, sname)

      h = xvalue(1)

      lbd = xvalue(2)

      ! Writing parameters

      write(unt,*)
      write(unt,*) " *** Geometry parameters ***"
      write(unt,*)

      write(unt,"(I23,    A)") kgb, " = kgb: Kind of geometry of the body"
      write(unt,"(ES23.16,A)")   w,  " =   w: Distance between south and north boundaries (m)"
      write(unt,"(ES23.16,A)")  fs, " =  fs: Factor of smoothing (0<=fs) (necessary for non-smooth surfaces)"
      write(unt,"(ES23.16,A)")  lo, " =  lo: Length of the ogive (m)"
      write(unt,"(ES23.16,A)")  lc, " =  lc: Length of the cylinder (m)"
      write(unt,"(ES23.16,A)")  rb, " =  rb: Base radius of the body (m)"
      write(unt,"(ES23.16,A)")   h, " =   h: Ratio rf/rb (frontal to base radius ratio)"
      write(unt,"(ES23.16,A)") lbd, " = lbd: Exponent of the power law ogive"
      write(unt,"(ES23.16,A)") aks, " = aks: Exponent of the power law for the south boundary"
      write(unt,"(ES23.16,A)") flo, " = flo: Fraction of lo"
      write(unt,"(ES23.16,A)") fp1, " = fp1: Fraction of volumes in the first part of the grid"
      write(unt,"(ES23.16,A)") fp2, " = fp2: Fraction of volumes in the second part of the grid"
      write(unt,"(ES23.16,A)") fp3, " = fp3: Fraction of volumes in the third part of the grid"

      ! Calculating the number of points in the first part of the grid

      nx1 = int( nx * fp1 )

      ! Calculating the number of points in the first and second parts of the grid

      nx2 = int( nx * fp2 ) + nx1

      ! Calculating the number of points in the 1st, 2nd and 3rd part of the grid

      nx3 = int( fp3 * nx ) + nx2

      ! Calculating frontal radius

      rf = h * rb

      ! Calculating the fraction of lo used in the 2nd part of the grid

      lo1 = flo * lo

      ! Calculating the body length

      lr = lo + lc


      ! ========================================================

      ! Generating south boundary

      j = 1

      ! Part I

      iiv = 1   ! Initial value of i

      ifv = nx1 ! Final value of i

      yiv = 0.d0 ! Initial value of y

      yfv = rf  ! Final value of y

      a1 = (yfv-yiv) * ( 1.d0 - (dble(ifv-iiv-1)/ dble(ifv-iiv)) ** aks ) ! Width of the first volume of the first part of the grid

      call get_gp_distribution(iiv, ifv, yiv, yfv, a1, ybs(iiv:ifv)) ! Output: last one

      do i = iiv, ifv

         np   = nx * (j-1) + i

         xbs(i) = 0.d0

         x(np) = xbs(i)

         y(np) = ybs(i)

      end do


      ! Part II

      iiv = nx1  ! Initial value of i

      ifv = nx2  ! Final value of i

      xiv = 0.d0 ! Initial value of x

      xfv = lo1  ! Final value of x

      ! Calculating the width of the first volume of Part II

      i = iiv

      a2 = sqrt( (xbs(i)-xbs(i-1))**2 + (ybs(i)-ybs(i-1))**2 )

      ! Creating the representing curve

      do i = 0, nr

         xr(i) = (xfv-xiv) * ( dble(i)/dble(nr) ) ** 2 + xiv

         yr(i) = rb * ( h + (1.d0-h) * ( xr(i) / lo ) ** lbd )

      end do

      ! Calculating the distribution curve according to a geometric progression
      ! distribution based on the curve length

      call get_gp_length_distribution(nr, ifv-iiv, a2, -1.d0, xr, yr &
         , xbs(iiv:ifv), ybs(iiv:ifv)) ! Output: last two

      do i = iiv, ifv

         np   = nx * (j-1) + i

         x(np) = xbs(i)

         y(np) = ybs(i)

      end do

      ! Part III

      iiv = nx2  ! Initial value of i

      ifv = nx3 ! Final value of i

      xiv = lo1 ! Initial value of x

      xfv = lo  ! Final value of x

      ! Calculating the width of the first volume of Part III

      i = iiv

      a3 = sqrt( (xbs(i)-xbs(i-1))**2 + (ybs(i)-ybs(i-1))**2 )

      ! Creating the representing curve

      do i = 0, nr

         xr(i) = (xfv-xiv) * ( dble(i)/dble(nr) ) + xiv

         yr(i) = rb * ( h + (1.d0-h) * ( xr(i) / lo ) ** lbd )

      end do

      ! Calculating the distribution curve according to a geometric progression
      ! distribution based on the curve length

      call get_gp_length_distribution(nr, ifv-iiv, a3, -1.d0, xr, yr &
         , xbs(iiv:ifv), ybs(iiv:ifv)) ! Output: last two

      do i = iiv, ifv

         np   = nx * (j-1) + i

         x(np) = xbs(i)

         y(np) = ybs(i)

      end do


      ! Part IV

      iiv = nx3  ! Initial value of i

      ifv = nx-1 ! Final value of i

      xiv = lo ! Initial value of x

      xfv = lo+lc  ! Final value of x

      ! Calculating the width of the first volume of Part IV

      i = iiv

      a4 = sqrt( (xbs(i)-xbs(i-1))**2 + (ybs(i)-ybs(i-1))**2 )

      call get_gp_distribution(iiv, ifv, xiv, xfv, a4, xbs(iiv:ifv)) ! Output: last one

      do i = iiv, ifv

         np   = nx * (j-1) + i

         ybs(i) = rb

         x(np) = xbs(i)

         y(np) = ybs(i)

      end do

      ! ========================================================

      ! Generating north boundary

      call get_north_boundary(nx, ny, w, fs, xbs, ybs, x, y) ! Output: last two

      ! ========================================================

   contains


      !> \brief Generates the north boundary of the domain based on the normal
      !! vectors from the south boundary. The directions of the normal vectors
      !! may be changed by a smoothing procedure.
      subroutine get_north_boundary(nx, ny, w, fs, xbs, ybs, x, y) ! Output: last two
         implicit none
         integer, intent(in)    :: nx  !< Number of volumes in the csi direction (real+fictitious)
         integer, intent(in)    :: ny  !< Number of volumes in the eta direction (real+fictitious)
         real(8), intent(in)    :: w   !< Width of the domain
         real(8), intent(in)    :: fs  !< Factor of smoothing (0<=fs)
         real(8), dimension(nx) :: xbs !< Distribution of points over the south boundary
         real(8), dimension(nx) :: ybs !< Distribution of points over the south boundary
         real(8), dimension(nx*ny), intent(out) :: x !< Coord. x at the northeast corner of the volume P (m)
         real(8), dimension(nx*ny), intent(out) :: y !< Coord. y at the northeast corner of the volume P (m)

         ! Inner variables

         integer :: i, j, k, np, ns
         real(8) :: raux
         real(8) :: tx
         real(8) :: ty
         real(8) :: tm
         real(8) :: vnx(nx-1)
         real(8) :: vny(nx-1)

         ! Number of iteractions for smoothing the direction vectors
         ns = int(nx * fs)

         ! Calculating the initial normal directions

         ! Direction of the west boundary
         vnx(1) = -1.d0
         vny(1) =  0.d0


         ! Direction of the east boundary
         vnx(nx-1) = 0.d0
         vny(nx-1) = 1.d0

         ! Directions normal to the body surface
         do i = 2, nx-2

            tx = ( xbs(i+1) - xbs(i-1) ) / 2.d0
            ty = ( ybs(i+1) - ybs(i-1) ) / 2.d0

            tm = sqrt( tx * tx + ty * ty )

            vnx(i) = - ty / tm
            vny(i) =   tx / tm

         end do

         ! Calculating the final directions (smoothing)

         do k = 1, ns

            do i = 2, nx-2

               vnx(i) = ( vnx(i+1) + vnx(i-1) ) / 2.d0

               vny(i) = ( vny(i+1) + vny(i-1) ) / 2.d0

               raux = dsqrt( vnx(i) * vnx(i) + vny(i) * vny(i) )

               vnx(i) = vnx(i) / raux

               vny(i) = vny(i) / raux

            end do

         end do

         ! Generating the north boundary

         j = ny-1

         do i = 1, nx-1

            np   = nx * (j-1) + i

            x(np) = xbs(i) + w * vnx(i)

            y(np) = ybs(i) + w * vny(i)

         end do

      end subroutine get_north_boundary


      !> \brief Calculates the node distribution based on a geometric progression distribution
      subroutine get_gp_distribution(iiv, ifv, xiv, xfv, a1, xb) ! Output: last one
         implicit none
         integer, intent(in)  :: iiv !< Initial value of i
         integer, intent(in)  :: ifv !< Final value of i
         real(8), intent(in)  :: xiv !< Initial value of x
         real(8), intent(in)  :: xfv !< Final value of x
         real(8), intent(in)  :: a1  !< Width of the first volume
         real(8), intent(out) :: xb(iiv:ifv) !< Distribution of points

         ! Inner variables

         integer :: i
         integer :: n
         real(8) :: r
         real(8) :: q


         n = ifv - iiv

         r = ( xfv - xiv ) / a1

         ! Calculating q
         call get_gp_ratio(n, r, q)

         xb(iiv) = xiv

         do i = iiv+1, ifv-1

            xb(i) = xb(i-1) + a1 * q ** ( i - iiv - 1 )

         end do

         xb(ifv) = xfv

      end subroutine get_gp_distribution


      !> \brief Calculates the base q of the geometric progression
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

           if ( 0.d0 < f(n, r, qi) * f(n, r, qm) ) then

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

      end subroutine get_gp_ratio

      !> \brief Function used in the subroutine 'get_gp_ratio'
      real(8) function f(n, r, q)
         implicit none
         integer, intent(in)  ::   n !< number of partitions
         real(8), intent(in)  ::   r !< l/a1
         real(8), intent(in)  ::   q !< q

         f = q ** n + r * ( 1.d0 - q ) - 1.d0

      end function f

      !> \brief Given a reference distribution of points (xbr,ybr)
      !! a new one (xb,yb) is created based on a geometric progression distributiuon
      !! based on the curve length.
      subroutine get_gp_length_distribution( nbr, nb, al, ar, xbr, ybr, xb, yb) ! Output: last two
         implicit none
         integer, intent(in)  :: nbr        !< Number of partitions over the reference curve (xbr,ybr)
         integer, intent(in)  :: nb         !< Number of partitions for the new distribution of points
         real(8), intent(in)  :: al         !< Length of the first partition on the left
         real(8), intent(in)  :: ar         !< Length of the first partition on the right
         real(8), intent(in)  :: xbr(0:nbr) !< Given distribution of points
         real(8), intent(in)  :: ybr(0:nbr) !< Given distribution of points
         real(8), intent(out) :: xb(0:nb)   !< New distribution of points
         real(8), intent(out) :: yb(0:nb)   !< New distribution of points

         ! Inner variables

         real(8) :: sbr(0:nbr) ! Length of each partition over the reference curve
         real(8) :: sb(0:nb)   ! Length of each partition over the new curve
         real(8) :: L          ! Total length


         ! Calculates the length of the curve
         call get_sbr(nbr, xbr, ybr, sbr) ! Output: last one


         ! Total length
         L = sbr(nbr)


         ! Selecting the kind of distribution

         if ( al < 0.d0 .and. ar < 0.d0 ) then

            ! Calculates the new points according to a uniform distribution
            call get_uniform_grid( nb, L, sb) ! Output: last one

         else if ( al > 0.d0 .and. ar < 0.d0 ) then

            ! Calculates the distribution of points according to a geometric
            ! progression. The length of the first partition is defined on the left.
            call get_gp_grid_left( nb, al, L, sb)

         else if ( al < 0.d0 .and. ar > 0.d0 ) then

            ! Calculates the distribution of points according to a geometric
            ! progression. The length of the first partition is defined on the right.
            call get_gp_grid_right( nb, ar, L, sb)

         else

            ! Calculates the distribution of points according to a geometric
            ! progression. The length of the first partition is defined on the left
            ! and on the right.
            call get_gp_grid_left_right(nb, al, ar, L, sb)

         end if


         ! Calculates the new distribution of points
         call get_new_distribution(nbr, nb, sbr, xbr, ybr, sb, xb, yb) ! Output: last two


      end subroutine get_gp_length_distribution

      !> \brief Calculates the length of the curve
      subroutine get_sbr(nbr, xbr, ybr, sbr) ! Output: last one
         implicit none
         integer, intent(in)  :: nbr        !< Number of partitions over the reference curve (xbr,ybr)
         real(8), intent(in)  :: xbr(0:nbr) !< Given distribution of points
         real(8), intent(in)  :: ybr(0:nbr) !< Given distribution of points
         real(8), intent(out) :: sbr(0:nbr) !< Length of the reference curve

         ! Inner variables

         integer :: i


         sbr(0) = 0.d0

         do i = 1, nbr

            sbr(i) = sbr(i-1) + sqrt( ( xbr(i) - xbr(i-1) )**2 + ( ybr(i) - ybr(i-1) )**2 )

         end do


      end subroutine get_sbr


      !> \brief Given the length distribution, the (xb,yb) is calculated
      subroutine get_new_distribution(nbr, nb, sbr, xbr, ybr, sb, xb, yb) ! Output: last two
         integer, intent(in)  :: nbr        !< Number of partitions over the reference curve (xbr,ybr)
         integer, intent(in)  :: nb         !< Number of partitions for the new distribution of points
         real(8), intent(in)  :: sbr(0:nbr) !< Length of the curve
         real(8), intent(in)  :: xbr(0:nbr) !< Given distribution of points
         real(8), intent(in)  :: ybr(0:nbr) !< Given distribution of points
         real(8), intent(in)  :: sb(0:nb)   !< New length distribution over the curve
         real(8), intent(out) :: xb(0:nb)   !< New distribution of points
         real(8), intent(out) :: yb(0:nb)   !< New distribution of points


         ! Inner variables

         integer :: i, j

         xb(0) = xbr(0)
         yb(0) = ybr(0)

         j = 0

         do i = 1, nb-1

            do j = j, nbr

               if( sbr(j) - sb(i) >= 0.d0 ) exit

            end do

            ! Linear interpolation

            xb(i) = xbr(j-1) + ( xbr(j)-xbr(j-1) ) * ( sb(i) - sbr(j-1) ) / ( sbr(j) - sbr(j-1) )
            yb(i) = ybr(j-1) + ( ybr(j)-ybr(j-1) ) * ( sb(i) - sbr(j-1) ) / ( sbr(j) - sbr(j-1) )

         end do

         xb(nb) = xbr(nbr)
         yb(nb) = ybr(nbr)

      end subroutine get_new_distribution



      !> \brief Calculates the uniform distribution of points
      subroutine get_uniform_grid( nb, L, sb) ! Output: last one
         implicit none
         integer, intent(in)  :: nb       !< Number of partitions over the new distribution of points
         real(8), intent(in)  :: L        !< Length of the curve
         real(8), intent(out) :: sb(0:nb) !< Length distribution over the curve

         ! Inner variables

         integer :: i

         real(8) :: ds

         ds = L / dble(nb)


         sb(0) = 0.d0

         do i = 1, nb

            sb(i) = sb(i-1) + ds

         end do

      end subroutine get_uniform_grid



      !> \brief Calculates the distribution of points according
      !! to a geometric progression. The length of the first partition
      !! is defined on the left.
      subroutine get_gp_grid_left(nb, al, L, sb)
         implicit none
         integer, intent(in)  :: nb       !< Number of partitions over the new distribution of points
         real(8), intent(in)  :: al       !< Length of the first patition on the left
         real(8), intent(in)  :: L        !< Length of the curve
         real(8), intent(out) :: sb(0:nb) !< Length distribution over the curve

         ! Inner variables

         integer :: i
         real(8) :: q

         ! Calculates the ratio q of the geometric progression
         call get_gp_ratio(nb, L/al, q)

         sb(0) = 0.d0

         do i = 1, nb-1

            sb(i) = sb(i-1) + al * q ** (i-1)

         end do

         sb(nb) = L

      end subroutine get_gp_grid_left


      !> \brief Calculates the distribution of points according
      !! to a geometric progression. The length of the first partition
      !! is defined on the right.
      subroutine get_gp_grid_right(nb, ar, L, sb)
         implicit none
         integer, intent(in)  :: nb       !< Number of partitions over the new distribution of points
         real(8), intent(in)  :: ar       !< Length of the first patition on the right
         real(8), intent(in)  :: L        !< Length of the curve
         real(8), intent(out) :: sb(0:nb) !< Length distribution over the curve

         ! Inner variables

         integer :: i
         real(8) :: q

         ! Calculates the ratio q of the geometric progression
         call get_gp_ratio(nb, L/ar, q)

         sb(0) = 0.d0

         do i = 1, nb-1

            sb(i) = sb(i-1) + ar * q ** (nb-i)

         end do

         sb(nb) = L

      end subroutine get_gp_grid_right



      !> \brief Calculates the distribution of points according
      !! to a geometric progression. The length of the first partition
      !! is defined on the left and on the right.
      subroutine get_gp_grid_left_right(nb, al, ar, L, sb)
         implicit none
         integer, intent(in)  :: nb       !< Number of partitions over the new distribution of points
         real(8), intent(in)  :: al       !< Length of the first patition on the left
         real(8), intent(in)  :: ar       !< Length of the first patition on the right
         real(8), intent(in)  :: L        !< Length of the curve
         real(8), intent(out) :: sb(0:nb) !< Length distribution over the curve

         ! Inner parameters

         integer, parameter :: nitmax = 10000

         ! Inner variables

         integer :: i
         integer :: nl
         integer :: nr
         real(8) :: Lm
         real(8) :: Li
         real(8) :: Lf
         real(8) :: acl
         real(8) :: acr


         nl = nb / 2

         nr = nb - nl

         Li = 0.d0

         Lf = L

         do i = 1, nitmax

            Lm = 0.5d0 * ( Li + Lf )

            call get_gp_grid_left (nl, al, Lm, sb( 0:nl))

            call get_gp_grid_right(nr, ar, L-Lm, sb(nl:nb))

            sb(nl:nb) = sb(nl:nb) + Lm

            acl = sb(nl) - sb(nl-1)

            acr = sb(nl+1) - sb(nl)

            if ( acl > acr ) then

               Lf = Lm

            else

               Li = Lm

            end if

            if ( abs(acl-acr) < 1.d-14 ) exit

         end do

      end subroutine get_gp_grid_left_right

   end subroutine get_grid_boundary_g12


   !> \brief Calculates the coordinates (x,y) of the north and south boundaries.
   !! The north boundary line is generated from the south boundary.
   !! The south boundary line is an bluff nose followed by a power law ogive of
   !! length lo and radius rb followed by a cylinder of radius rb and length lc.
   !! The distribution of points over the south boundary is separated in three
   !! parts:
   !! Part   I: geometric progression (bluff nose+power law)
   !! Part  II: geometric progression (power law)
   !! Part III: geometric progression (cylinder)
   !! The distribution of points over the north boundary is  based on the normal
   !! vectors from the south boundary. The directions of the normal vectors
   !! may be changed by a smoothing procedure.
   !! Adapted to DEPP optimizer.
   subroutine get_grid_boundary_g13(nx, ny, unt, fgeom, lr, rb, x, y) ! Output: last two
      implicit none
      integer, intent(in) :: nx  !< Number of volumes in the csi direction (real+fictitious)
      integer, intent(in) :: ny  !< Number of volumes in the eta direction (real+fictitious)
      integer, intent(in) :: unt !< Unit where the input parameters will be printed
      character(len=*), intent(in) :: fgeom !< File of the geometric parameters

      real(8), intent(out) :: lr  !< Length of the body (m)
      real(8), intent(out) :: rb  !< Base radius of the body (m)
      real(8), dimension(nx*ny), intent(out) :: x !< Coord. x at the northeast corner of the volume P (m)
      real(8), dimension(nx*ny), intent(out) :: y !< Coord. y at the northeast corner of the volume P (m)

      ! Parameters

      integer, parameter :: nr = 10000 ! Number of partitions for the representing curve (xr,yr)


      ! Inner variables
      integer :: i, j, np ! Dummy indexes

      integer :: kgb ! Kind of geometry of the body
      integer :: nx1 ! Number of points in the first part of the grid along the eta lines
      integer :: nx2 ! Number of points in the first and second parts of the grid along the eta lines
      integer :: nx3 ! Number of points in the first to third parts of the grid along the eta lines
      integer :: iiv ! Initial value of i
      integer :: ifv ! Final value of i

      real(8) :: a2  ! Width of the first volume of Part II
      real(8) :: a3  ! Width of the first volume of Part III

      real(8) :: xiv ! Initial value of x
      real(8) :: xfv ! Final value of x

      real(8) :: lbd ! Exponent of the power law ogive

      real(8) :: lo  ! Ogive length (m)
      real(8) :: w   ! Distance between south and north boundaries (m)
      real(8) :: fs  ! Factor of smoothing (0<=fs) (necessary for non-smooth surfaces)
      real(8) :: lc  ! Cylinder length (m)
      real(8) :: lo1 ! Length of the first part of the ogive (m)
      real(8) :: rf  ! Radius of the nose tip (m)
      real(8) :: h   ! Ratio rf/rb
      real(8) :: aks ! Exponent of the power law for the calculation of the width of the first volume
      real(8) :: flo ! Fraction of lo
      real(8) :: fp1 ! Fraction of points in the first part of the grid
      real(8) :: fp2 ! Fraction of points in the second part of the grid


      real(8), dimension(nx) :: xbs ! Distribution of points over the south boundary
      real(8), dimension(nx) :: ybs ! Distribution of points over the south boundary

      real(8), dimension(0:nr) :: xr ! Distribution of points over the representing curve
      real(8), dimension(0:nr) :: yr ! Distribution of points over the representing curve


      ! DEPP variables
      integer, parameter :: nu = 2 ! number of unknowns
      integer            :: ind    ! number of the individual
      integer,       dimension(nu) :: xopt   ! Optimization checker
      character(10), dimension(nu) :: xname  ! Name of parameters
      real(8),       dimension(nu) :: xvalue ! parameters
      character(200)               :: sname  ! simulation name


      ! Reading parameters

      open(13,file=fgeom)

      read(13,*) kgb  ! Kind of geometry of the body
      read(13,*) w    ! Distance between south and north boundaries (m)
      read(13,*) fs   ! Factor of smoothing (0<=fs) (necessary for non-smooth surfaces)
      read(13,*) lo   ! Length of the ogive (m)
      read(13,*) lc   ! Length of the cylinder (m)
      read(13,*) rb   ! Base radius of the body (m)
      read(13,*) xname(1) ! Parameter name
      read(13,*) xopt(1)  ! Will this parameter be optimized? ( 0 = no, 1 = yes )
      read(13,*) xvalue(1)! Parameter value
      read(13,*) xname(2) ! Parameter name
      read(13,*) xopt(2)  ! Will this parameter be optimized? ( 0 = no, 1 = yes )
      read(13,*) xvalue(2)! Parameter value
      read(13,*) aks  ! Exponent of the power law for the south boundary
      read(13,*) flo  ! Fraction of lo
      read(13,*) fp1  ! Fraction of volumes in the first part
      read(13,*) fp2  ! Fraction of volumes in the second part

      close(13)

      ! Reads the parameters from DEPP
      call depp_get_parameters(nu, xopt, xname, xvalue, ind, sname)

      h = xvalue(1)

      lbd = xvalue(2)

      ! Writing parameters

      write(unt,*)
      write(unt,*) " *** Geometry parameters ***"
      write(unt,*)

      write(unt,"(I23,    A)") kgb, " = kgb: Kind of geometry of the body"
      write(unt,"(ES23.16,A)")   w,  " =   w: Distance between south and north boundaries (m)"
      write(unt,"(ES23.16,A)")  fs, " =  fs: Factor of smoothing (0<=fs) (necessary for non-smooth surfaces)"
      write(unt,"(ES23.16,A)")  lo, " =  lo: Length of the ogive (m)"
      write(unt,"(ES23.16,A)")  lc, " =  lc: Length of the cylinder (m)"
      write(unt,"(ES23.16,A)")  rb, " =  rb: Base radius of the body (m)"
      write(unt,"(ES23.16,A)")   h, " =   h: Ratio rf/rb (frontal to base radius ratio)"
      write(unt,"(ES23.16,A)") lbd, " = lbd: Exponent of the power law ogive"
      write(unt,"(ES23.16,A)") aks, " = aks: Exponent of the power law for the south boundary"
      write(unt,"(ES23.16,A)") flo, " = flo: Fraction of lo"
      write(unt,"(ES23.16,A)") fp1, " = fp1: Fraction of volumes in the first part of the grid"
      write(unt,"(ES23.16,A)") fp2, " = fp2: Fraction of volumes in the second part of the grid"

      ! Calculating the number of points in the first part of the grid

      nx1 = int( (nx-2) * fp1 )

      ! Calculating the number of points in the first and second parts of the grid

      nx2 = int( (nx-2) * fp2 ) + nx1

      ! Calculating the number of points in the 1st, 2nd and 3rd part of the grid

      nx3 = nx-1

      ! Calculating frontal radius

      rf = h * rb

      ! Calculating the fraction of lo used in the 2nd part of the grid

      lo1 = flo * lo

      ! Calculating the body length

      lr = lo + lc


      ! ========================================================

      ! Generating south boundary

      j = 1

      ! Part I


      ! Creating the representing curve

      xr(0) = 0.d0

      yr(0) = 0.d0

      xiv = 0.d0

      xfv = lo1

      do i = 1, nr

         xr(i) = (xfv-xiv) * ( dble(i-1)/dble(nr-1) ) ** 2 + xiv

         yr(i) = rb * ( h + (1.d0-h) * ( xr(i) / lo ) ** lbd )

      end do

      ! Calculates the distribution of points based on the arclength and the
      ! double exponential distribution


      iiv = 1   ! Initial value of i

      ifv = nx1 ! Final value of i

      call get_dexp_m1_ditribution( nr, ifv-iiv, aks, rf, xr, yr, xbs(iiv:ifv) &
         , ybs(iiv:ifv)) ! Output: last two

      do i = iiv, ifv

         np   = nx * (j-1) + i

         x(np) = xbs(i)

         y(np) = ybs(i)

      end do


      ! Part II

      iiv = nx1  ! Initial value of i

      ifv = nx2 ! Final value of i

      xiv = lo1 ! Initial value of x

      xfv = lo  ! Final value of x

      ! Calculating the width of the first volume of Part II

      i = iiv

      a2 = sqrt( (xbs(i)-xbs(i-1))**2 + (ybs(i)-ybs(i-1))**2 )

      ! Creating the representing curve

      do i = 0, nr

         xr(i) = (xfv-xiv) * ( dble(i)/dble(nr) ) + xiv

         yr(i) = rb * ( h + (1.d0-h) * ( xr(i) / lo ) ** lbd )

      end do

      ! Calculating the distribution curve according to a geometric progression
      ! distribution based on the curve length

      call get_gp_length_distribution(nr, ifv-iiv, a2, -1.d0, xr, yr &
         , xbs(iiv:ifv), ybs(iiv:ifv)) ! Output: last two

      do i = iiv, ifv

         np   = nx * (j-1) + i

         x(np) = xbs(i)

         y(np) = ybs(i)

      end do


      ! Part III

      iiv = nx2  ! Initial value of i

      ifv = nx3  ! Final value of i

      xiv = lo ! Initial value of x

      xfv = lo+lc  ! Final value of x

      ! Calculating the width of the first volume of Part IV

      i = iiv

      a3 = sqrt( (xbs(i)-xbs(i-1))**2 + (ybs(i)-ybs(i-1))**2 )

      call get_gp_distribution(iiv, ifv, xiv, xfv, a3, xbs(iiv:ifv)) ! Output: last one

      do i = iiv, ifv

         np   = nx * (j-1) + i

         ybs(i) = rb

         x(np) = xbs(i)

         y(np) = ybs(i)

      end do


      ! ========================================================

      ! Generating north boundary

      call get_north_boundary(nx, ny, w, fs, xbs, ybs, x, y) ! Output: last two

      ! ========================================================

   contains


      !> \brief Generates the north boundary of the domain based on the normal
      !! vectors from the south boundary. The directions of the normal vectors
      !! may be changed by a smoothing procedure.
      subroutine get_north_boundary(nx, ny, w, fs, xbs, ybs, x, y) ! Output: last two
         implicit none
         integer, intent(in)    :: nx  !< Number of volumes in the csi direction (real+fictitious)
         integer, intent(in)    :: ny  !< Number of volumes in the eta direction (real+fictitious)
         real(8), intent(in)    :: w   !< Width of the domain
         real(8), intent(in)    :: fs  !< Factor of smoothing (0<=fs)
         real(8), dimension(nx),    intent(in)  :: xbs !< Distribution of points over the south boundary
         real(8), dimension(nx),    intent(in)  :: ybs !< Distribution of points over the south boundary
         real(8), dimension(nx*ny), intent(out) :: x   !< Coord. x at the northeast corner of the volume P (m)
         real(8), dimension(nx*ny), intent(out) :: y   !< Coord. y at the northeast corner of the volume P (m)

         ! Inner variables

         integer :: i, j, k, np, ns
         real(8) :: raux
         real(8) :: tx
         real(8) :: ty
         real(8) :: tm
         real(8) :: vnx(nx-1)
         real(8) :: vny(nx-1)

         ! Number of iteractions for smoothing the direction vectors
         ns = int(nx * fs)

         ! Calculating the initial normal directions

         ! Direction of the west boundary
         vnx(1) = -1.d0
         vny(1) =  0.d0


         ! Direction of the east boundary
         vnx(nx-1) = 0.d0
         vny(nx-1) = 1.d0

         ! Directions normal to the body surface
         do i = 2, nx-2

            tx = ( xbs(i+1) - xbs(i-1) ) / 2.d0
            ty = ( ybs(i+1) - ybs(i-1) ) / 2.d0

            tm = sqrt( tx * tx + ty * ty )

            vnx(i) = - ty / tm
            vny(i) =   tx / tm

         end do

         ! Calculating the final directions (smoothing)

         do k = 1, ns

            do i = 2, nx-2

               vnx(i) = ( vnx(i+1) + vnx(i-1) ) / 2.d0

               vny(i) = ( vny(i+1) + vny(i-1) ) / 2.d0

               raux = dsqrt( vnx(i) * vnx(i) + vny(i) * vny(i) )

               vnx(i) = vnx(i) / raux

               vny(i) = vny(i) / raux

            end do

         end do

         ! Generating the north boundary

         j = ny-1

         do i = 1, nx-1

            np   = nx * (j-1) + i

            x(np) = xbs(i) + w * vnx(i)

            y(np) = ybs(i) + w * vny(i)

         end do

      end subroutine get_north_boundary


      !> \brief Calculates the node distribution based on a geometric progression distribution
      subroutine get_gp_distribution(iiv, ifv, xiv, xfv, a1, xb) ! Output: last one
         implicit none
         integer, intent(in)  :: iiv !< Initial value of i
         integer, intent(in)  :: ifv !< Final value of i
         real(8), intent(in)  :: xiv !< Initial value of x
         real(8), intent(in)  :: xfv !< Final value of x
         real(8), intent(in)  :: a1  !< Width of the first volume
         real(8), intent(out) :: xb(iiv:ifv) !< Distribution of points

         ! Inner variables

         integer :: i
         integer :: n
         real(8) :: r
         real(8) :: q


         n = ifv - iiv

         r = ( xfv - xiv ) / a1

         ! Calculating q
         call get_gp_ratio(n, r, q)

         xb(iiv) = xiv

         do i = iiv+1, ifv-1

            xb(i) = xb(i-1) + a1 * q ** ( i - iiv - 1 )

         end do

         xb(ifv) = xfv

      end subroutine get_gp_distribution


      !> \brief Calculates the base q of the geometric progression
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

           if ( 0.d0 < f(n, r, qi) * f(n, r, qm) ) then

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

      end subroutine get_gp_ratio

      !> \brief Function used in the subroutine 'get_gp_ratio'
      real(8) function f(n, r, q)
         implicit none
         integer, intent(in)  ::   n !< number of partitions
         real(8), intent(in)  ::   r !< l/a1
         real(8), intent(in)  ::   q !< q

         f = q ** n + r * ( 1.d0 - q ) - 1.d0

      end function f

      !> \brief Given a reference distribution of points (xbr,ybr)
      !! a new one (xb,yb) is created based on a geometric progression distributiuon
      !! based on the curve length.
      subroutine get_gp_length_distribution( nbr, nb, al, ar, xbr, ybr, xb, yb) ! Output: last two
         implicit none
         integer, intent(in)  :: nbr        !< Number of partitions over the reference curve (xbr,ybr)
         integer, intent(in)  :: nb         !< Number of partitions for the new distribution of points
         real(8), intent(in)  :: al         !< Length of the first partition on the left
         real(8), intent(in)  :: ar         !< Length of the first partition on the right
         real(8), intent(in)  :: xbr(0:nbr) !< Given distribution of points
         real(8), intent(in)  :: ybr(0:nbr) !< Given distribution of points
         real(8), intent(out) :: xb(0:nb)   !< New distribution of points
         real(8), intent(out) :: yb(0:nb)   !< New distribution of points

         ! Inner variables

         real(8) :: sbr(0:nbr) ! Length of each partition over the reference curve
         real(8) :: sb(0:nb)   ! Length of each partition over the new curve
         real(8) :: L          ! Total length


         ! Calculates the length of the curve
         call get_sbr(nbr, xbr, ybr, sbr) ! Output: last one


         ! Total length
         L = sbr(nbr)


         ! Selecting the kind of distribution

         if ( al < 0.d0 .and. ar < 0.d0 ) then

            ! Calculates the new points according to a uniform distribution
            call get_uniform_grid( nb, L, sb) ! Output: last one

         else if ( al > 0.d0 .and. ar < 0.d0 ) then

            ! Calculates the distribution of points according to a geometric
            ! progression. The length of the first partition is defined on the left.
            call get_gp_grid_left( nb, al, L, sb)

         else if ( al < 0.d0 .and. ar > 0.d0 ) then

            ! Calculates the distribution of points according to a geometric
            ! progression. The length of the first partition is defined on the right.
            call get_gp_grid_right( nb, ar, L, sb)

         else

            ! Calculates the distribution of points according to a geometric
            ! progression. The length of the first partition is defined on the left
            ! and on the right.
            call get_gp_grid_left_right(nb, al, ar, L, sb)

         end if


         ! Calculates the new distribution of points
         call get_new_distribution(nbr, nb, sbr, xbr, ybr, sb, xb, yb) ! Output: last two


      end subroutine get_gp_length_distribution

      !> \brief Calculates the length of the curve
      subroutine get_sbr(nbr, xbr, ybr, sbr) ! Output: last one
         implicit none
         integer, intent(in)  :: nbr        !< Number of partitions over the reference curve (xbr,ybr)
         real(8), intent(in)  :: xbr(0:nbr) !< Given distribution of points
         real(8), intent(in)  :: ybr(0:nbr) !< Given distribution of points
         real(8), intent(out) :: sbr(0:nbr) !< Length of the reference curve

         ! Inner variables

         integer :: i


         sbr(0) = 0.d0

         do i = 1, nbr

            sbr(i) = sbr(i-1) + sqrt( ( xbr(i) - xbr(i-1) )**2 + ( ybr(i) - ybr(i-1) )**2 )

         end do


      end subroutine get_sbr


      !> \brief Given the length distribution, the (xb,yb) is calculated
      subroutine get_new_distribution(nbr, nb, sbr, xbr, ybr, sb, xb, yb) ! Output: last two
         integer, intent(in)  :: nbr        !< Number of partitions over the reference curve (xbr,ybr)
         integer, intent(in)  :: nb         !< Number of partitions for the new distribution of points
         real(8), intent(in)  :: sbr(0:nbr) !< Length of the curve
         real(8), intent(in)  :: xbr(0:nbr) !< Given distribution of points
         real(8), intent(in)  :: ybr(0:nbr) !< Given distribution of points
         real(8), intent(in)  :: sb(0:nb)   !< New length distribution over the curve
         real(8), intent(out) :: xb(0:nb)   !< New distribution of points
         real(8), intent(out) :: yb(0:nb)   !< New distribution of points


         ! Inner variables

         integer :: i, j

         xb(0) = xbr(0)
         yb(0) = ybr(0)

         j = 0

         do i = 1, nb-1

            do j = j, nbr

               if( sbr(j) - sb(i) >= 0.d0 ) exit

            end do

            ! Linear interpolation

            xb(i) = xbr(j-1) + ( xbr(j)-xbr(j-1) ) * ( sb(i) - sbr(j-1) ) / ( sbr(j) - sbr(j-1) )
            yb(i) = ybr(j-1) + ( ybr(j)-ybr(j-1) ) * ( sb(i) - sbr(j-1) ) / ( sbr(j) - sbr(j-1) )

         end do

         xb(nb) = xbr(nbr)
         yb(nb) = ybr(nbr)

      end subroutine get_new_distribution



      !> \brief Calculates the uniform distribution of points
      subroutine get_uniform_grid( nb, L, sb) ! Output: last one
         implicit none
         integer, intent(in)  :: nb       !< Number of partitions over the new distribution of points
         real(8), intent(in)  :: L        !< Length of the curve
         real(8), intent(out) :: sb(0:nb) !< Length distribution over the curve

         ! Inner variables

         integer :: i

         real(8) :: ds

         ds = L / dble(nb)


         sb(0) = 0.d0

         do i = 1, nb

            sb(i) = sb(i-1) + ds

         end do

      end subroutine get_uniform_grid



      !> \brief Calculates the distribution of points according
      !! to a geometric progression. The length of the first partition
      !! is defined on the left.
      subroutine get_gp_grid_left(nb, al, L, sb)
         implicit none
         integer, intent(in)  :: nb       !< Number of partitions over the new distribution of points
         real(8), intent(in)  :: al       !< Length of the first patition on the left
         real(8), intent(in)  :: L        !< Length of the curve
         real(8), intent(out) :: sb(0:nb) !< Length distribution over the curve

         ! Inner variables

         integer :: i
         real(8) :: q

         ! Calculates the ratio q of the geometric progression
         call get_gp_ratio(nb, L/al, q)

         sb(0) = 0.d0

         do i = 1, nb-1

            sb(i) = sb(i-1) + al * q ** (i-1)

         end do

         sb(nb) = L

      end subroutine get_gp_grid_left


      !> \brief Calculates the distribution of points according
      !! to a geometric progression. The length of the first partition
      !! is defined on the right.
      subroutine get_gp_grid_right(nb, ar, L, sb)
         implicit none
         integer, intent(in)  :: nb       !< Number of partitions over the new distribution of points
         real(8), intent(in)  :: ar       !< Length of the first patition on the right
         real(8), intent(in)  :: L        !< Length of the curve
         real(8), intent(out) :: sb(0:nb) !< Length distribution over the curve

         ! Inner variables

         integer :: i
         real(8) :: q

         ! Calculates the ratio q of the geometric progression
         call get_gp_ratio(nb, L/ar, q)

         sb(0) = 0.d0

         do i = 1, nb-1

            sb(i) = sb(i-1) + ar * q ** (nb-i)

         end do

         sb(nb) = L

      end subroutine get_gp_grid_right



      !> \brief Calculates the distribution of points according
      !! to a geometric progression. The length of the first partition
      !! is defined on the left and on the right.
      subroutine get_gp_grid_left_right(nb, al, ar, L, sb)
         implicit none
         integer, intent(in)  :: nb       !< Number of partitions over the new distribution of points
         real(8), intent(in)  :: al       !< Length of the first patition on the left
         real(8), intent(in)  :: ar       !< Length of the first patition on the right
         real(8), intent(in)  :: L        !< Length of the curve
         real(8), intent(out) :: sb(0:nb) !< Length distribution over the curve

         ! Inner parameters

         integer, parameter :: nitmax = 10000

         ! Inner variables

         integer :: i
         integer :: nl
         integer :: nr
         real(8) :: Lm
         real(8) :: Li
         real(8) :: Lf
         real(8) :: acl
         real(8) :: acr


         nl = nb / 2

         nr = nb - nl

         Li = 0.d0

         Lf = L

         do i = 1, nitmax

            Lm = 0.5d0 * ( Li + Lf )

            call get_gp_grid_left (nl, al, Lm, sb( 0:nl))

            call get_gp_grid_right(nr, ar, L-Lm, sb(nl:nb))

            sb(nl:nb) = sb(nl:nb) + Lm

            acl = sb(nl) - sb(nl-1)

            acr = sb(nl+1) - sb(nl)

            if ( acl > acr ) then

               Lf = Lm

            else

               Li = Lm

            end if

            if ( abs(acl-acr) < 1.d-14 ) exit

         end do

      end subroutine get_gp_grid_left_right


      !> \brief Calculates a distribution of points based on the arclength and
      !! on the double-exponential rule
      subroutine get_dexp_m1_ditribution( nbr, nb, aks, zcv, xbr, ybr, xb, yb) ! Output: last two
         implicit none
         integer, intent(in)  :: nbr        !< Number of partitions over the reference curve (xbr,ybr)
         integer, intent(in)  :: nb         !< Number of partitions for the new distribution of points
         real(8), intent(in)  :: aks        !< Exponent for the calculation of acv
         real(8), intent(in)  :: zcv        !< Central value of z
         real(8), intent(in)  :: xbr(0:nbr) !< Given distribution of points
         real(8), intent(in)  :: ybr(0:nbr) !< Given distribution of points
         real(8), intent(out) :: xb(0:nb)   !< New distribution of points
         real(8), intent(out) :: yb(0:nb)   !< New distribution of points

         ! Inner variables

         real(8) :: sbr(0:nbr) ! Length of each partition over the reference curve
         real(8) :: sb(0:nb)   ! Length of each partition over the new curve
         real(8) :: L          ! Total length
         real(8) :: acv        ! Length of the volume closer to the central point


         ! Calculates the length of the curve
         call get_sbr(nbr, xbr, ybr, sbr) ! Output: last one


         ! Total length
         L = sbr(nbr)


         ! Calculating the double-exponential distribution of arclength

         acv = L / ( nb ** aks )

         if ( acv > zcv ) then

            write(*,*) "Error. get_dexp_m1_ditribution: there is no resolution"&
               , " for this grid. Stopping."

            stop

         end if

         call get_dexp_m1(nb, 0.d0, zcv, L, acv, sb)


         ! Calculates the new distribution of points
         call get_new_distribution(nbr, nb, sbr, xbr, ybr, sb, xb, yb) ! Output: last two


      end subroutine get_dexp_m1_ditribution

      !> \brief Calculates a distribution of points based on the double-exponential
      !! rule. The coefficients of this rule are set by defining the position of
      !! an internal point (zcv) and the width (acv) of the volumes beside this
      !! point. The points are concentrated around zcv.
      subroutine get_dexp_m1(np, ziv, zcv, zfv, acv, z)
         implicit none
         integer, intent(in)  :: np  !< Number of partitions
         real(8), intent(in)  :: ziv !< Initial value of z
         real(8), intent(in)  :: zcv !< Central value of z (exponential's matching point)
         real(8), intent(in)  :: zfv !< Final value of z
         real(8), intent(in)  :: acv !< Length of the volumes beside zcv
         real(8), intent(out) :: z(0:np) !< Distribution of points

         ! Inner variables

         integer :: i
         integer :: npc
         real(8) :: t
         real(8) :: dt
         real(8) :: tc
         real(8) :: r
         real(8) :: A
         real(8) :: ac
         real(8) :: zc


         ! Initializing variables

         zc = zcv / ( zfv - ziv )

         ac = acv / ( zfv - ziv )

         dt = 1.d0 / dble(np)

         r = ac / ( 1.d0 - zc )



         ! Cheking conditions

         if ( ( ac / zc > 1.d0 ) .or. ( r > 1.d0 ) ) then

            write(*,*) "Error. get_dexp_m1: ac > zc or ac > (1-zc). Stopping."

            stop

         end if

         if ( ac > dt ) then

            write(*,*) "Error. get_dexp_m1: ac > dt. Stopping."

            stop

         end if


         ! Searching the value of A

         call adexp_m1(dt, r, zc, A)


         ! Calculating the fraction of volumes for each exponential

         tc = 1.d0 + log( 1.d0 + ( exp(-A*dt) - 1.d0 ) / r ) / A

         npc = nint(tc*np)



         ! Calculating the normallized distribution of points

         ! First exponential

         do i = 0, npc

            t = tc * dble(i) / dble(npc)

            z(i) = zc * ( exp(A*t) - 1.d0 ) / ( exp(A*tc) - 1.d0 )

         end do

         ! Second exponential

         do i = npc, np

            t = ( 1.d0-tc ) * dble(i-npc) / dble(np-npc) + tc

            z(i) = zc+(1.d0-zc)*(exp(-A*(t-tc))-1.d0)/(exp(-A*(1.d0-tc))-1.d0)

         end do

         z(0)  = 0.d0

         z(np) = 1.d0

         ! Calculating the extended distribution of points

         z = ( zfv - ziv ) * z + ziv

      end subroutine get_dexp_m1

      !> \brief Defines the transcendental function for finding A
      real(8) function fdexp_m1(dt, r, zc, A)
         implicit none
         real(8), intent(in) :: dt
         real(8), intent(in) :: r
         real(8), intent(in) :: zc
         real(8), intent(in) :: A

         ! Inner varialbles

         real(8) :: raux1, raux2, raux3

         if ( abs(A) > 1.d-6 ) then

            raux1 = exp(A) * ( 1.d0 + ( exp(-A*dt) - 1.d0 ) / r )
            raux2 = (2.d0*zc-1.d0) / zc
            raux3 = (1.d0-zc) / zc

            fdexp_m1 = -exp(-A) * raux1 * raux1 +  raux1 * raux2 + raux3

         else

            fdexp_m1 = (r*zc-r+dt)*A/(r*zc) &
               + ((r**2-2*dt**2)*zc-r**2+(2*dt-dt**2)*r)*A**2/(r**2*zc)/2.0

         end if

      end function fdexp_m1


      !> \brief Finds the value of A using the bissection method
      subroutine adexp_m1(dt, r, zc, A)
         implicit none
         real(8), intent(in)  :: dt
         real(8), intent(in)  :: r
         real(8), intent(in)  :: zc
         real(8), intent(out) :: A

         ! Parameters

         integer, parameter :: itmax = 1000

         ! Inner variables

         integer :: k
         real(8) :: Ai, Af, Am, Fi, Fm, Ff



         ! Initializing the search interval

         Ai = -1.d-18

         Af = -100.d0



         ! Checking for roots

         Fi = fdexp_m1(dt, r, zc, Ai)

         Ff = fdexp_m1(dt, r, zc, Af)

         if ( Fi * Ff > 0.d0 ) then

            write(*,*) "Error. adexp_m1: there is none or more than one roots. Stopping."

            stop

         end if


         ! Looking for a solution

         do k = 1, itmax

            Am = 0.5d0 * ( Ai + Af )

            Fi = fdexp_m1(dt, r, zc, Ai)

            Fm = fdexp_m1(dt, r, zc, Am)

            if ( Fi * Fm <= 0.d0 ) then

               Af = Am

            else

               Ai = Am

            end if

            if ( abs(Af-Ai) < 1.d-14 ) exit

         end do

         if ( k >= itmax ) then

            write(*,*) "Error. adexp_m1: number of iteractions exceeded. Stopping."

            stop

         end if

         A = 0.5d0 * ( Ai + Af )

      end subroutine adexp_m1


   end subroutine get_grid_boundary_g13



   !> \brief Calculates the coordinates (x,y) of the north and south boundaries.
   !! The north boundary line is generated from the south boundary.
   !! The south boundary line is the Newton optimized shape of
   !! length lo and radius rb followed by a cylinder of radius rb and length lc.
   !! The distribution of points over the south boundary is separated in three
   !! parts:
   !! Part   I: geometric progression (bluff nose+Newton section)
   !! Part  II: geometric progression (Newton section)
   !! Part III: geometric progression (cylinder)
   !! The distribution of points over the north boundary is  based on the normal
   !! vectors from the south boundary. The directions of the normal vectors
   !! may be changed by a smoothing procedure.
   subroutine get_grid_boundary_g14(nx, ny, unt, fgeom, lr, rb, x, y) ! Output: last two
      implicit none
      integer, intent(in) :: nx  !< Number of volumes in the csi direction (real+fictitious)
      integer, intent(in) :: ny  !< Number of volumes in the eta direction (real+fictitious)
      integer, intent(in) :: unt !< Unit where the input parameters will be printed
      character(len=*), intent(in) :: fgeom !< File of the geometric parameters

      real(8), intent(out) :: lr  !< Length of the body (m)
      real(8), intent(out) :: rb  !< Base radius of the body (m)
      real(8), dimension(nx*ny), intent(out) :: x !< Coord. x at the northeast corner of the volume P (m)
      real(8), dimension(nx*ny), intent(out) :: y !< Coord. y at the northeast corner of the volume P (m)

      ! Parameters

      integer, parameter :: nr = 10000 ! Number of partitions for the representing curve (xr,yr)


      ! Inner variables
      integer :: i, j, np ! Dummy indexes

      integer :: kgb ! Kind of geometry of the body
      integer :: nx1 ! Number of points in the first part of the grid along the eta lines
      integer :: nx2 ! Number of points in the first and second parts of the grid along the eta lines
      integer :: nx3 ! Number of points in the first to third parts of the grid along the eta lines
      integer :: iiv ! Initial value of i
      integer :: ifv ! Final value of i
      integer :: nr1 ! Fraction of nr used to represent the Part I

      real(8) :: a2  ! Width of the first volume of Part II
      real(8) :: a3  ! Width of the first volume of Part III

      real(8) :: xiv ! Initial value of x
      real(8) :: xfv ! Final value of x

      real(8) :: lo  ! Ogive length (m)
      real(8) :: w   ! Distance between south and north boundaries (m)
      real(8) :: fs  ! Factor of smoothing (0<=fs) (necessary for non-smooth surfaces)
      real(8) :: lc  ! Cylinder length (m)
      real(8) :: lo1 ! Length of the first part of the ogive (m)
      real(8) :: rf  ! Radius of the nose tip (m)
      real(8) :: h   ! Ratio rf/rb
      real(8) :: aks ! Exponent of the power law for the calculation of the width of the first volume
      real(8) :: flo ! Fraction of lo
      real(8) :: fp1 ! Fraction of points in the first part of the grid
      real(8) :: fp2 ! Fraction of points in the second part of the grid


      real(8), dimension(nx) :: xbs ! Distribution of points over the south boundary
      real(8), dimension(nx) :: ybs ! Distribution of points over the south boundary

      real(8), dimension(0:nr) :: xr ! Distribution of points over the representing curve
      real(8), dimension(0:nr) :: yr ! Distribution of points over the representing curve


      ! Reading parameters

      open(13,file=fgeom)

      read(13,*) kgb  ! Kind of geometry of the body
      read(13,*) w    ! Distance between south and north boundaries (m)
      read(13,*) fs   ! Factor of smoothing (0<=fs) (necessary for non-smooth surfaces)
      read(13,*) lo   ! Length of the ogive (m)
      read(13,*) lc   ! Length of the cylinder (m)
      read(13,*) rb   ! Base radius of the body (m)
      read(13,*) aks  ! Exponent of the power law for the south boundary
      read(13,*) flo  ! Fraction of lo
      read(13,*) fp1  ! Fraction of volumes in the first part
      read(13,*) fp2  ! Fraction of volumes in the second part

      close(13)


      ! Creating the Newton shape

      call get_newton_shape(nr, rb, lo, xr, yr)


      ! Calculating the number of points in the first part of the grid

      nx1 = int( (nx-2) * fp1 )

      ! Calculating the number of points in the first and second parts of the grid

      nx2 = int( (nx-2) * fp2 ) + nx1

      ! Calculating the number of points in the 1st, 2nd and 3rd part of the grid

      nx3 = nx-1

      ! Reading rf

      rf = yr(1)

      ! Calculating h

      h = rf / rb

      ! Calculating the fraction of lo used in the 2nd part of the grid

      lo1 = flo * lo

      do i = 1, nr

         if ( lo1 < xr(i) ) exit

      end do

      nr1 = i

      lo1 = xr(nr1)


      ! Calculating the body length

      lr = lo + lc


      ! Writing parameters

      write(unt,*)
      write(unt,*) " *** Geometry parameters ***"
      write(unt,*)

      write(unt,"(I23,    A)") kgb, " = kgb: Kind of geometry of the body"
      write(unt,"(ES23.16,A)")   w,  " =   w: Distance between south and north boundaries (m)"
      write(unt,"(ES23.16,A)")  fs, " =  fs: Factor of smoothing (0<=fs) (necessary for non-smooth surfaces)"
      write(unt,"(ES23.16,A)")  lo, " =  lo: Length of the ogive (m)"
      write(unt,"(ES23.16,A)")  lc, " =  lc: Length of the cylinder (m)"
      write(unt,"(ES23.16,A)")  rb, " =  rb: Base radius of the body (m)"
      write(unt,"(ES23.16,A)")   h, " =   h: Ratio rf/rb (frontal to base radius ratio)"
      write(unt,"(ES23.16,A)") aks, " = aks: Exponent of the power law for the south boundary"
      write(unt,"(ES23.16,A)") flo, " = flo: Fraction of lo"
      write(unt,"(ES23.16,A)") fp1, " = fp1: Fraction of volumes in the first part of the grid"
      write(unt,"(ES23.16,A)") fp2, " = fp2: Fraction of volumes in the second part of the grid"


      ! ========================================================

      ! Generating south boundary

      j = 1

      ! Part I

      ! Calculates the distribution of points based on the arclength and the
      ! double exponential distribution

      iiv = 1   ! Initial value of i

      ifv = nx1 ! Final value of i

      call get_dexp_m1_ditribution( nr1, ifv-iiv, aks, rf, xr(0:nr1), yr(0:nr1)&
         , xbs(iiv:ifv), ybs(iiv:ifv)) ! Output: last two

      do i = iiv, ifv

         np   = nx * (j-1) + i

         x(np) = xbs(i)

         y(np) = ybs(i)

      end do


      ! Part II

      iiv = nx1  ! Initial value of i

      ifv = nx2 ! Final value of i

      xiv = lo1 ! Initial value of x

      xfv = lo  ! Final value of x

      ! Calculating the width of the first volume of Part II

      i = iiv

      a2 = sqrt( (xbs(i)-xbs(i-1))**2 + (ybs(i)-ybs(i-1))**2 )

      ! Calculating the distribution curve according to a geometric progression
      ! distribution based on the curve length

      call get_gp_length_distribution(nr-nr1, ifv-iiv, a2, -1.d0, xr(nr1:nr) &
         , yr(nr1:nr), xbs(iiv:ifv), ybs(iiv:ifv)) ! Output: last two


      do i = iiv, ifv

         np   = nx * (j-1) + i

         x(np) = xbs(i)

         y(np) = ybs(i)

      end do


      ! Part III

      iiv = nx2  ! Initial value of i

      ifv = nx3  ! Final value of i

      xiv = lo ! Initial value of x

      xfv = lo+lc  ! Final value of x

      ! Calculating the width of the first volume of Part III

      i = iiv

      a3 = sqrt( (xbs(i)-xbs(i-1))**2 + (ybs(i)-ybs(i-1))**2 )

      call get_gp_distribution(iiv, ifv, xiv, xfv, a3, xbs(iiv:ifv)) ! Output: last one

      do i = iiv, ifv

         np   = nx * (j-1) + i

         ybs(i) = rb

         x(np) = xbs(i)

         y(np) = ybs(i)

      end do


      ! ========================================================

      ! Generating north boundary

      call get_north_boundary(nx, ny, w, fs, xbs, ybs, x, y) ! Output: last two

      ! ========================================================

   contains


      !> \brief Generates the north boundary of the domain based on the normal
      !! vectors from the south boundary. The directions of the normal vectors
      !! may be changed by a smoothing procedure.
      subroutine get_north_boundary(nx, ny, w, fs, xbs, ybs, x, y) ! Output: last two
         implicit none
         integer, intent(in)    :: nx  !< Number of volumes in the csi direction (real+fictitious)
         integer, intent(in)    :: ny  !< Number of volumes in the eta direction (real+fictitious)
         real(8), intent(in)    :: w   !< Width of the domain
         real(8), intent(in)    :: fs  !< Factor of smoothing (0<=fs)
         real(8), dimension(nx),    intent(in)  :: xbs !< Distribution of points over the south boundary
         real(8), dimension(nx),    intent(in)  :: ybs !< Distribution of points over the south boundary
         real(8), dimension(nx*ny), intent(out) :: x   !< Coord. x at the northeast corner of the volume P (m)
         real(8), dimension(nx*ny), intent(out) :: y   !< Coord. y at the northeast corner of the volume P (m)

         ! Inner variables

         integer :: i, j, k, np, ns
         real(8) :: raux
         real(8) :: tx
         real(8) :: ty
         real(8) :: tm
         real(8) :: vnx(nx-1)
         real(8) :: vny(nx-1)

         ! Number of iteractions for smoothing the direction vectors
         ns = int(nx * fs)

         ! Calculating the initial normal directions

         ! Direction of the west boundary
         vnx(1) = -1.d0
         vny(1) =  0.d0


         ! Direction of the east boundary
         vnx(nx-1) = 0.d0
         vny(nx-1) = 1.d0

         ! Directions normal to the body surface
         do i = 2, nx-2

            tx = ( xbs(i+1) - xbs(i-1) ) / 2.d0
            ty = ( ybs(i+1) - ybs(i-1) ) / 2.d0

            tm = sqrt( tx * tx + ty * ty )

            vnx(i) = - ty / tm
            vny(i) =   tx / tm

         end do

         ! Calculating the final directions (smoothing)

         do k = 1, ns

            do i = 2, nx-2

               vnx(i) = ( vnx(i+1) + vnx(i-1) ) / 2.d0

               vny(i) = ( vny(i+1) + vny(i-1) ) / 2.d0

               raux = dsqrt( vnx(i) * vnx(i) + vny(i) * vny(i) )

               vnx(i) = vnx(i) / raux

               vny(i) = vny(i) / raux

            end do

         end do

         ! Generating the north boundary

         j = ny-1

         do i = 1, nx-1

            np   = nx * (j-1) + i

            x(np) = xbs(i) + w * vnx(i)

            y(np) = ybs(i) + w * vny(i)

         end do

      end subroutine get_north_boundary


      !> \brief Calculates the node distribution based on a geometric progression distribution
      subroutine get_gp_distribution(iiv, ifv, xiv, xfv, a1, xb) ! Output: last one
         implicit none
         integer, intent(in)  :: iiv !< Initial value of i
         integer, intent(in)  :: ifv !< Final value of i
         real(8), intent(in)  :: xiv !< Initial value of x
         real(8), intent(in)  :: xfv !< Final value of x
         real(8), intent(in)  :: a1  !< Width of the first volume
         real(8), intent(out) :: xb(iiv:ifv) !< Distribution of points

         ! Inner variables

         integer :: i
         integer :: n
         real(8) :: r
         real(8) :: q


         n = ifv - iiv

         r = ( xfv - xiv ) / a1

         ! Calculating q
         call get_gp_ratio(n, r, q)

         xb(iiv) = xiv

         do i = iiv+1, ifv-1

            xb(i) = xb(i-1) + a1 * q ** ( i - iiv - 1 )

         end do

         xb(ifv) = xfv

      end subroutine get_gp_distribution


      !> \brief Calculates the base q of the geometric progression
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

           if ( 0.d0 < f(n, r, qi) * f(n, r, qm) ) then

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

      end subroutine get_gp_ratio

      !> \brief Function used in the subroutine 'get_gp_ratio'
      real(8) function f(n, r, q)
         implicit none
         integer, intent(in)  ::   n !< number of partitions
         real(8), intent(in)  ::   r !< l/a1
         real(8), intent(in)  ::   q !< q

         f = q ** n + r * ( 1.d0 - q ) - 1.d0

      end function f

      !> \brief Given a reference distribution of points (xbr,ybr)
      !! a new one (xb,yb) is created based on a geometric progression distributiuon
      !! based on the curve length.
      subroutine get_gp_length_distribution( nbr, nb, al, ar, xbr, ybr, xb, yb) ! Output: last two
         implicit none
         integer, intent(in)  :: nbr        !< Number of partitions over the reference curve (xbr,ybr)
         integer, intent(in)  :: nb         !< Number of partitions for the new distribution of points
         real(8), intent(in)  :: al         !< Length of the first partition on the left
         real(8), intent(in)  :: ar         !< Length of the first partition on the right
         real(8), intent(in)  :: xbr(0:nbr) !< Given distribution of points
         real(8), intent(in)  :: ybr(0:nbr) !< Given distribution of points
         real(8), intent(out) :: xb(0:nb)   !< New distribution of points
         real(8), intent(out) :: yb(0:nb)   !< New distribution of points

         ! Inner variables

         real(8) :: sbr(0:nbr) ! Length of each partition over the reference curve
         real(8) :: sb(0:nb)   ! Length of each partition over the new curve
         real(8) :: L          ! Total length


         ! Calculates the length of the curve
         call get_sbr(nbr, xbr, ybr, sbr) ! Output: last one


         ! Total length
         L = sbr(nbr)


         ! Selecting the kind of distribution

         if ( al < 0.d0 .and. ar < 0.d0 ) then

            ! Calculates the new points according to a uniform distribution
            call get_uniform_grid( nb, L, sb) ! Output: last one

         else if ( al > 0.d0 .and. ar < 0.d0 ) then

            ! Calculates the distribution of points according to a geometric
            ! progression. The length of the first partition is defined on the left.
            call get_gp_grid_left( nb, al, L, sb)

         else if ( al < 0.d0 .and. ar > 0.d0 ) then

            ! Calculates the distribution of points according to a geometric
            ! progression. The length of the first partition is defined on the right.
            call get_gp_grid_right( nb, ar, L, sb)

         else

            ! Calculates the distribution of points according to a geometric
            ! progression. The length of the first partition is defined on the left
            ! and on the right.
            call get_gp_grid_left_right(nb, al, ar, L, sb)

         end if


         ! Calculates the new distribution of points
         call get_new_distribution(nbr, nb, sbr, xbr, ybr, sb, xb, yb) ! Output: last two


      end subroutine get_gp_length_distribution

      !> \brief Calculates the length of the curve
      subroutine get_sbr(nbr, xbr, ybr, sbr) ! Output: last one
         implicit none
         integer, intent(in)  :: nbr        !< Number of partitions over the reference curve (xbr,ybr)
         real(8), intent(in)  :: xbr(0:nbr) !< Given distribution of points
         real(8), intent(in)  :: ybr(0:nbr) !< Given distribution of points
         real(8), intent(out) :: sbr(0:nbr) !< Length of the reference curve

         ! Inner variables

         integer :: i


         sbr(0) = 0.d0

         do i = 1, nbr

            sbr(i) = sbr(i-1) + sqrt( ( xbr(i) - xbr(i-1) )**2 + ( ybr(i) - ybr(i-1) )**2 )

         end do


      end subroutine get_sbr


      !> \brief Given the length distribution, the (xb,yb) is calculated
      subroutine get_new_distribution(nbr, nb, sbr, xbr, ybr, sb, xb, yb) ! Output: last two
         integer, intent(in)  :: nbr        !< Number of partitions over the reference curve (xbr,ybr)
         integer, intent(in)  :: nb         !< Number of partitions for the new distribution of points
         real(8), intent(in)  :: sbr(0:nbr) !< Length of the curve
         real(8), intent(in)  :: xbr(0:nbr) !< Given distribution of points
         real(8), intent(in)  :: ybr(0:nbr) !< Given distribution of points
         real(8), intent(in)  :: sb(0:nb)   !< New length distribution over the curve
         real(8), intent(out) :: xb(0:nb)   !< New distribution of points
         real(8), intent(out) :: yb(0:nb)   !< New distribution of points


         ! Inner variables

         integer :: i, j

         xb(0) = xbr(0)
         yb(0) = ybr(0)

         j = 0

         do i = 1, nb-1

            do j = j, nbr

               if( sbr(j) - sb(i) >= 0.d0 ) exit

            end do

            ! Linear interpolation

            xb(i) = xbr(j-1) + ( xbr(j)-xbr(j-1) ) * ( sb(i) - sbr(j-1) ) / ( sbr(j) - sbr(j-1) )
            yb(i) = ybr(j-1) + ( ybr(j)-ybr(j-1) ) * ( sb(i) - sbr(j-1) ) / ( sbr(j) - sbr(j-1) )

         end do

         xb(nb) = xbr(nbr)
         yb(nb) = ybr(nbr)

      end subroutine get_new_distribution



      !> \brief Calculates the uniform distribution of points
      subroutine get_uniform_grid( nb, L, sb) ! Output: last one
         implicit none
         integer, intent(in)  :: nb       !< Number of partitions over the new distribution of points
         real(8), intent(in)  :: L        !< Length of the curve
         real(8), intent(out) :: sb(0:nb) !< Length distribution over the curve

         ! Inner variables

         integer :: i

         real(8) :: ds

         ds = L / dble(nb)


         sb(0) = 0.d0

         do i = 1, nb

            sb(i) = sb(i-1) + ds

         end do

      end subroutine get_uniform_grid



      !> \brief Calculates the distribution of points according
      !! to a geometric progression. The length of the first partition
      !! is defined on the left.
      subroutine get_gp_grid_left(nb, al, L, sb)
         implicit none
         integer, intent(in)  :: nb       !< Number of partitions over the new distribution of points
         real(8), intent(in)  :: al       !< Length of the first patition on the left
         real(8), intent(in)  :: L        !< Length of the curve
         real(8), intent(out) :: sb(0:nb) !< Length distribution over the curve

         ! Inner variables

         integer :: i
         real(8) :: q

         ! Calculates the ratio q of the geometric progression
         call get_gp_ratio(nb, L/al, q)

         sb(0) = 0.d0

         do i = 1, nb-1

            sb(i) = sb(i-1) + al * q ** (i-1)

         end do

         sb(nb) = L

      end subroutine get_gp_grid_left


      !> \brief Calculates the distribution of points according
      !! to a geometric progression. The length of the first partition
      !! is defined on the right.
      subroutine get_gp_grid_right(nb, ar, L, sb)
         implicit none
         integer, intent(in)  :: nb       !< Number of partitions over the new distribution of points
         real(8), intent(in)  :: ar       !< Length of the first patition on the right
         real(8), intent(in)  :: L        !< Length of the curve
         real(8), intent(out) :: sb(0:nb) !< Length distribution over the curve

         ! Inner variables

         integer :: i
         real(8) :: q

         ! Calculates the ratio q of the geometric progression
         call get_gp_ratio(nb, L/ar, q)

         sb(0) = 0.d0

         do i = 1, nb-1

            sb(i) = sb(i-1) + ar * q ** (nb-i)

         end do

         sb(nb) = L

      end subroutine get_gp_grid_right



      !> \brief Calculates the distribution of points according
      !! to a geometric progression. The length of the first partition
      !! is defined on the left and on the right.
      subroutine get_gp_grid_left_right(nb, al, ar, L, sb)
         implicit none
         integer, intent(in)  :: nb       !< Number of partitions over the new distribution of points
         real(8), intent(in)  :: al       !< Length of the first patition on the left
         real(8), intent(in)  :: ar       !< Length of the first patition on the right
         real(8), intent(in)  :: L        !< Length of the curve
         real(8), intent(out) :: sb(0:nb) !< Length distribution over the curve

         ! Inner parameters

         integer, parameter :: nitmax = 10000

         ! Inner variables

         integer :: i
         integer :: nl
         integer :: nr
         real(8) :: Lm
         real(8) :: Li
         real(8) :: Lf
         real(8) :: acl
         real(8) :: acr


         nl = nb / 2

         nr = nb - nl

         Li = 0.d0

         Lf = L

         do i = 1, nitmax

            Lm = 0.5d0 * ( Li + Lf )

            call get_gp_grid_left (nl, al, Lm, sb( 0:nl))

            call get_gp_grid_right(nr, ar, L-Lm, sb(nl:nb))

            sb(nl:nb) = sb(nl:nb) + Lm

            acl = sb(nl) - sb(nl-1)

            acr = sb(nl+1) - sb(nl)

            if ( acl > acr ) then

               Lf = Lm

            else

               Li = Lm

            end if

            if ( abs(acl-acr) < 1.d-14 ) exit

         end do

      end subroutine get_gp_grid_left_right


      !> \brief Calculates a distribution of points based on the arclength and
      !! on the double-exponential rule
      subroutine get_dexp_m1_ditribution( nbr, nb, aks, zcv, xbr, ybr, xb, yb) ! Output: last two
         implicit none
         integer, intent(in)  :: nbr        !< Number of partitions over the reference curve (xbr,ybr)
         integer, intent(in)  :: nb         !< Number of partitions for the new distribution of points
         real(8), intent(in)  :: aks        !< Exponent for the calculation of acv
         real(8), intent(in)  :: zcv        !< Central value of z
         real(8), intent(in)  :: xbr(0:nbr) !< Given distribution of points
         real(8), intent(in)  :: ybr(0:nbr) !< Given distribution of points
         real(8), intent(out) :: xb(0:nb)   !< New distribution of points
         real(8), intent(out) :: yb(0:nb)   !< New distribution of points

         ! Inner variables

         real(8) :: sbr(0:nbr) ! Length of each partition over the reference curve
         real(8) :: sb(0:nb)   ! Length of each partition over the new curve
         real(8) :: L          ! Total length
         real(8) :: acv        ! Length of the volume closer to the central point


         ! Calculates the length of the curve
         call get_sbr(nbr, xbr, ybr, sbr) ! Output: last one


         ! Total length
         L = sbr(nbr)


         ! Calculating the double-exponential distribution of arclength

         acv = L / ( nb ** aks )

         if ( acv > zcv ) then

            write(*,*) "Error. get_dexp_m1_ditribution: there is no resolution"&
               , " for this grid. Stopping."

            stop

         end if

         call get_dexp_m1(nb, 0.d0, zcv, L, acv, sb)


         ! Calculates the new distribution of points
         call get_new_distribution(nbr, nb, sbr, xbr, ybr, sb, xb, yb) ! Output: last two


      end subroutine get_dexp_m1_ditribution

      !> \brief Calculates a distribution of points based on the double-exponential
      !! rule. The coefficients of this rule are set by defining the position of
      !! an internal point (zcv) and the width (acv) of the volumes beside this
      !! point. The points are concentrated around zcv.
      subroutine get_dexp_m1(np, ziv, zcv, zfv, acv, z)
         implicit none
         integer, intent(in)  :: np  !< Number of partitions
         real(8), intent(in)  :: ziv !< Initial value of z
         real(8), intent(in)  :: zcv !< Central value of z (exponential's matching point)
         real(8), intent(in)  :: zfv !< Final value of z
         real(8), intent(in)  :: acv !< Length of the volumes beside zcv
         real(8), intent(out) :: z(0:np) !< Distribution of points

         ! Inner variables

         integer :: i
         integer :: npc
         real(8) :: t
         real(8) :: dt
         real(8) :: tc
         real(8) :: r
         real(8) :: A
         real(8) :: ac
         real(8) :: zc


         ! Initializing variables

         zc = zcv / ( zfv - ziv )

         ac = acv / ( zfv - ziv )

         dt = 1.d0 / dble(np)

         r = ac / ( 1.d0 - zc )



         ! Cheking conditions

         if ( ( ac / zc > 1.d0 ) .or. ( r > 1.d0 ) ) then

            write(*,*) "Error. get_dexp_m1: ac > zc or ac > (1-zc). Stopping."

            stop

         end if

         if ( ac > dt ) then

            write(*,*) "Error. get_dexp_m1: ac > dt. Stopping."

            stop

         end if


         ! Searching the value of A

         call adexp_m1(dt, r, zc, A)


         ! Calculating the fraction of volumes for each exponential

         tc = 1.d0 + log( 1.d0 + ( exp(-A*dt) - 1.d0 ) / r ) / A

         npc = nint(tc*np)



         ! Calculating the normallized distribution of points

         ! First exponential

         do i = 0, npc

            t = tc * dble(i) / dble(npc)

            z(i) = zc * ( exp(A*t) - 1.d0 ) / ( exp(A*tc) - 1.d0 )

         end do

         ! Second exponential

         do i = npc, np

            t = ( 1.d0-tc ) * dble(i-npc) / dble(np-npc) + tc

            z(i) = zc+(1.d0-zc)*(exp(-A*(t-tc))-1.d0)/(exp(-A*(1.d0-tc))-1.d0)

         end do

         z(0)  = 0.d0

         z(np) = 1.d0

         ! Calculating the extended distribution of points

         z = ( zfv - ziv ) * z + ziv

      end subroutine get_dexp_m1

      !> \brief Defines the transcendental function for finding A
      real(8) function fdexp_m1(dt, r, zc, A)
         implicit none
         real(8), intent(in) :: dt
         real(8), intent(in) :: r
         real(8), intent(in) :: zc
         real(8), intent(in) :: A

         ! Inner varialbles

         real(8) :: raux1, raux2, raux3

         if ( abs(A) > 1.d-6 ) then

            raux1 = exp(A) * ( 1.d0 + ( exp(-A*dt) - 1.d0 ) / r )
            raux2 = (2.d0*zc-1.d0) / zc
            raux3 = (1.d0-zc) / zc

            fdexp_m1 = -exp(-A) * raux1 * raux1 +  raux1 * raux2 + raux3

         else

            fdexp_m1 = (r*zc-r+dt)*A/(r*zc) &
               + ((r**2-2*dt**2)*zc-r**2+(2*dt-dt**2)*r)*A**2/(r**2*zc)/2.0

         end if

      end function fdexp_m1


      !> \brief Finds the value of A using the bissection method
      subroutine adexp_m1(dt, r, zc, A)
         implicit none
         real(8), intent(in)  :: dt
         real(8), intent(in)  :: r
         real(8), intent(in)  :: zc
         real(8), intent(out) :: A

         ! Parameters

         integer, parameter :: itmax = 1000

         ! Inner variables

         integer :: k
         real(8) :: Ai, Af, Am, Fi, Fm, Ff



         ! Initializing the search interval

         Ai = -1.d-18

         Af = -100.d0



         ! Checking for roots

         Fi = fdexp_m1(dt, r, zc, Ai)

         Ff = fdexp_m1(dt, r, zc, Af)

         if ( Fi * Ff > 0.d0 ) then

            write(*,*) "Error. adexp_m1: there is none or more than one roots. Stopping."

            stop

         end if


         ! Looking for a solution

         do k = 1, itmax

            Am = 0.5d0 * ( Ai + Af )

            Fi = fdexp_m1(dt, r, zc, Ai)

            Fm = fdexp_m1(dt, r, zc, Am)

            if ( Fi * Fm <= 0.d0 ) then

               Af = Am

            else

               Ai = Am

            end if

            if ( abs(Af-Ai) < 1.d-14 ) exit

         end do

         if ( k >= itmax ) then

            write(*,*) "Error. adexp_m1: number of iteractions exceeded. Stopping."

            stop

         end if

         A = 0.5d0 * ( Ai + Af )

      end subroutine adexp_m1

      !> \brief Calculates the Newton optimized shape
      subroutine get_newton_shape(nr, rb, lo, xr, yr)
         implicit none
         integer, intent(in)  :: nr !< Number of partitions
         real(8), intent(in)  :: rb !< Base radius
         real(8), intent(in)  :: lo !< Ogive length
         real(8), intent(out) :: xr(0:nr) !< x coordinate
         real(8), intent(out) :: yr(0:nr) !< y coordinate

         ! Parameters
         integer, parameter :: nitmax = 2000 ! Maximum number of iteractions
         real(8), parameter :: tol = 1.d-14  ! Tolerance

         ! Inner variables
         integer :: i   ! Dummy index
         real(8) :: fr  ! Fineness ratio
         real(8) :: z   ! Slope
         real(8) :: zbi ! Initial value of zb
         real(8) :: zb  ! Slope at the rear of the body
         real(8) :: zbf ! Final value of zb
         real(8) :: fbi ! fbi=f(zbi)
         real(8) :: fb  !  fb=f(zb)
         real(8) :: h   ! h = rf / rb ( radius ratio )

         ! Fineness ratio
         fr = lo / ( 2.d0 * rb )

         ! Calculating the zb
         zbi = 1.d-15

         zbf = 1.d0

         fbi = 0.75d0 + zbi**2 + zbi**4 * ( log(zbi) - 1.75d0 ) &
            - 2.d0 * fr * zbi * ( 1.d0 + zbi**2 ) ** 2


         do i = 1, nitmax

            zb = 0.5d0 * ( zbi + zbf )

            fb = 0.75d0 + zb**2 + zb**4 * ( log(zb) - 1.75d0 ) &
               - 2.d0 * fr * zb * ( 1.d0 + zb**2 ) ** 2

            if ( fbi * fb < 0.d0 ) then

               zbf = zb

            else

               zbi = zb

               fbi = fb

            end if

            if ( abs(zbf-zbi) < tol ) exit

         end do

         if ( i >= nitmax ) then

            write(*,*) "Error. get_newton_shape: exceeded the maximum number of it."

            stop

         end if

         h = 4.d0 * zb ** 3 / ( 1.d0 + zb**2 )**2

         xr(0) = 0.d0
         yr(0) = 0.d0

         do i = 1, nr

            z = (zb-1.d0) * dble(i-1) / dble(nr-1) + 1.d0

            xr(i) = lo * h / (8.d0*fr) * ( 0.75d0 / z**4 + 1.d0 / z**2 - 1.75d0 + log(z) )

            yr(i) = 0.25d0 * rb * h * ( 1.d0 + z**2 )**2 / z**3

         end do

         xr(nr) = lo

         yr(nr) = rb

      end subroutine get_newton_shape
   end subroutine get_grid_boundary_g14


   !> \brief Calculates the coordinates (x,y) of the north and south boundaries.
   !! The north boundary line is generated from the south boundary.
   !! The south boundary line is the von Karman shape of
   !! length lo and radius rb followed by a cylinder of radius rb and length lc.
   !! The distribution of points over the south boundary is separated in three
   !! parts:
   !! Part   I: geometric progression
   !! Part  II: geometric progression
   !! Part III: geometric progression (cylinder)
   !! The distribution of points over the north boundary is  based on the normal
   !! vectors from the south boundary. The directions of the normal vectors
   !! may be changed by a smoothing procedure.
   subroutine get_grid_boundary_g15(nx, ny, unt, fgeom, lr, rb, x, y) ! Output: last two
      implicit none
      integer, intent(in) :: nx  !< Number of volumes in the csi direction (real+fictitious)
      integer, intent(in) :: ny  !< Number of volumes in the eta direction (real+fictitious)
      integer, intent(in) :: unt !< Unit where the input parameters will be printed
      character(len=*), intent(in) :: fgeom !< File of the geometric parameters

      real(8), intent(out) :: lr  !< Length of the body (m)
      real(8), intent(out) :: rb  !< Base radius of the body (m)
      real(8), dimension(nx*ny), intent(out) :: x !< Coord. x at the northeast corner of the volume P (m)
      real(8), dimension(nx*ny), intent(out) :: y !< Coord. y at the northeast corner of the volume P (m)

      ! Parameters

      integer, parameter :: nr = 10000 ! Number of partitions for the representing curve (xr,yr)


      ! Inner variables
      integer :: i, j, np ! Dummy indexes

      integer :: kgb ! Kind of geometry of the body
      integer :: nx1 ! Number of points in the first part of the grid along the eta lines
      integer :: nx2 ! Number of points in the first and second parts of the grid along the eta lines
      integer :: nx3 ! Number of points in the first to third parts of the grid along the eta lines
      integer :: iiv ! Initial value of i
      integer :: ifv ! Final value of i
      integer :: nr1 ! Fraction of nr used to represent the Part I

      real(8) :: a1  ! Width of the first volume of Part I
      real(8) :: a2  ! Width of the first volume of Part II
      real(8) :: a3  ! Width of the first volume of Part III

      real(8) :: xiv ! Initial value of x
      real(8) :: xfv ! Final value of x

      real(8) :: lo  ! Ogive length (m)
      real(8) :: w   ! Distance between south and north boundaries (m)
      real(8) :: fs  ! Factor of smoothing (0<=fs) (necessary for non-smooth surfaces)
      real(8) :: lc  ! Cylinder length (m)
      real(8) :: lo1 ! Length of the first part of the ogive (m)
      real(8) :: aks ! Exponent of the power law for the calculation of the width of the first volume
      real(8) :: flo ! Fraction of lo
      real(8) :: fp1 ! Fraction of points in the first part of the grid
      real(8) :: fp2 ! Fraction of points in the second part of the grid


      real(8), dimension(nx) :: xbs ! Distribution of points over the south boundary
      real(8), dimension(nx) :: ybs ! Distribution of points over the south boundary

      real(8), dimension(0:nr) :: xr ! Distribution of points over the representing curve
      real(8), dimension(0:nr) :: yr ! Distribution of points over the representing curve


      ! Reading parameters

      open(13,file=fgeom)

      read(13,*) kgb  ! Kind of geometry of the body
      read(13,*) w    ! Distance between south and north boundaries (m)
      read(13,*) fs   ! Factor of smoothing (0<=fs) (necessary for non-smooth surfaces)
      read(13,*) lo   ! Length of the ogive (m)
      read(13,*) lc   ! Length of the cylinder (m)
      read(13,*) rb   ! Base radius of the body (m)
      read(13,*) aks  ! Exponent of the power law for the south boundary
      read(13,*) flo  ! Fraction of lo
      read(13,*) fp1  ! Fraction of volumes in the first part
      read(13,*) fp2  ! Fraction of volumes in the second part

      close(13)


      ! Creating the Newton shape

      call get_karman_shape(nr, rb, lo, xr, yr)


      ! Calculating the number of points in the first part of the grid

      nx1 = int( (nx-2) * fp1 )

      ! Calculating the number of points in the first and second parts of the grid

      nx2 = int( (nx-2) * fp2 ) + nx1

      ! Calculating the number of points in the 1st, 2nd and 3rd part of the grid

      nx3 = nx-1

      ! Calculating the fraction of lo used in the 2nd part of the grid

      lo1 = flo * lo

      do i = 1, nr

         if ( lo1 < xr(i) ) exit

      end do

      nr1 = i

      lo1 = xr(nr1)


      ! Calculating the body length

      lr = lo + lc


      ! Writing parameters

      write(unt,*)
      write(unt,*) " *** Geometry parameters ***"
      write(unt,*)

      write(unt,"(I23,    A)") kgb, " = kgb: Kind of geometry of the body"
      write(unt,"(ES23.16,A)")   w,  " =   w: Distance between south and north boundaries (m)"
      write(unt,"(ES23.16,A)")  fs, " =  fs: Factor of smoothing (0<=fs) (necessary for non-smooth surfaces)"
      write(unt,"(ES23.16,A)")  lo, " =  lo: Length of the ogive (m)"
      write(unt,"(ES23.16,A)")  lc, " =  lc: Length of the cylinder (m)"
      write(unt,"(ES23.16,A)")  rb, " =  rb: Base radius of the body (m)"
      write(unt,"(ES23.16,A)") aks, " = aks: Exponent of the power law for the south boundary"
      write(unt,"(ES23.16,A)") flo, " = flo: Fraction of lo"
      write(unt,"(ES23.16,A)") fp1, " = fp1: Fraction of volumes in the first part of the grid"
      write(unt,"(ES23.16,A)") fp2, " = fp2: Fraction of volumes in the second part of the grid"


      ! ========================================================

      ! Generating south boundary

      j = 1

      ! Part I

      ! Calculating the distribution curve according to a geometric progression
      ! distribution based on the curve length

      iiv = 1   ! Initial value of i

      ifv = nx1 ! Final value of i

      xiv = 0.d0

      xfv = lo1

      a1 = (xfv-xiv) / (ifv-iiv) ** aks ! Width of the first volume of the first part of the grid

      call get_gp_length_distribution(nr1, ifv-iiv, a1, -1.d0, xr(0:nr1) &
         , yr(0:nr1), xbs(iiv:ifv), ybs(iiv:ifv)) ! Output: last two

      do i = iiv, ifv

         np   = nx * (j-1) + i

         x(np) = xbs(i)

         y(np) = ybs(i)

      end do


      ! Part II

      iiv = nx1  ! Initial value of i

      ifv = nx2 ! Final value of i

      xiv = lo1 ! Initial value of x

      xfv = lo  ! Final value of x

      ! Calculating the width of the first volume of Part II

      i = iiv

      a2 = sqrt( (xbs(i)-xbs(i-1))**2 + (ybs(i)-ybs(i-1))**2 )

      ! Calculating the distribution curve according to a geometric progression
      ! distribution based on the curve length

      call get_gp_length_distribution(nr-nr1, ifv-iiv, a2, -1.d0, xr(nr1:nr) &
         , yr(nr1:nr), xbs(iiv:ifv), ybs(iiv:ifv)) ! Output: last two


      do i = iiv, ifv

         np   = nx * (j-1) + i

         x(np) = xbs(i)

         y(np) = ybs(i)

      end do


      ! Part III

      iiv = nx2  ! Initial value of i

      ifv = nx3  ! Final value of i

      xiv = lo ! Initial value of x

      xfv = lo+lc  ! Final value of x

      ! Calculating the width of the first volume of Part III

      i = iiv

      a3 = sqrt( (xbs(i)-xbs(i-1))**2 + (ybs(i)-ybs(i-1))**2 )

      call get_gp_distribution(iiv, ifv, xiv, xfv, a3, xbs(iiv:ifv)) ! Output: last one

      do i = iiv, ifv

         np   = nx * (j-1) + i

         ybs(i) = rb

         x(np) = xbs(i)

         y(np) = ybs(i)

      end do


      ! ========================================================

      ! Generating north boundary

      call get_north_boundary(nx, ny, w, fs, xbs, ybs, x, y) ! Output: last two

      ! ========================================================

   contains


      !> \brief Generates the north boundary of the domain based on the normal
      !! vectors from the south boundary. The directions of the normal vectors
      !! may be changed by a smoothing procedure.
      subroutine get_north_boundary(nx, ny, w, fs, xbs, ybs, x, y) ! Output: last two
         implicit none
         integer, intent(in)    :: nx  !< Number of volumes in the csi direction (real+fictitious)
         integer, intent(in)    :: ny  !< Number of volumes in the eta direction (real+fictitious)
         real(8), intent(in)    :: w   !< Width of the domain
         real(8), intent(in)    :: fs  !< Factor of smoothing (0<=fs)
         real(8), dimension(nx),    intent(in)  :: xbs !< Distribution of points over the south boundary
         real(8), dimension(nx),    intent(in)  :: ybs !< Distribution of points over the south boundary
         real(8), dimension(nx*ny), intent(out) :: x   !< Coord. x at the northeast corner of the volume P (m)
         real(8), dimension(nx*ny), intent(out) :: y   !< Coord. y at the northeast corner of the volume P (m)

         ! Inner variables

         integer :: i, j, k, np, ns
         real(8) :: raux
         real(8) :: tx
         real(8) :: ty
         real(8) :: tm
         real(8) :: vnx(nx-1)
         real(8) :: vny(nx-1)

         ! Number of iteractions for smoothing the direction vectors
         ns = int(nx * fs)

         ! Calculating the initial normal directions

         ! Direction of the west boundary
         vnx(1) = -1.d0
         vny(1) =  0.d0


         ! Direction of the east boundary
         vnx(nx-1) = 0.d0
         vny(nx-1) = 1.d0

         ! Directions normal to the body surface
         do i = 2, nx-2

            tx = ( xbs(i+1) - xbs(i-1) ) / 2.d0
            ty = ( ybs(i+1) - ybs(i-1) ) / 2.d0

            tm = sqrt( tx * tx + ty * ty )

            vnx(i) = - ty / tm
            vny(i) =   tx / tm

         end do

         ! Calculating the final directions (smoothing)

         do k = 1, ns

            do i = 2, nx-2

               vnx(i) = ( vnx(i+1) + vnx(i-1) ) / 2.d0

               vny(i) = ( vny(i+1) + vny(i-1) ) / 2.d0

               raux = dsqrt( vnx(i) * vnx(i) + vny(i) * vny(i) )

               vnx(i) = vnx(i) / raux

               vny(i) = vny(i) / raux

            end do

         end do

         ! Generating the north boundary

         j = ny-1

         do i = 1, nx-1

            np   = nx * (j-1) + i

            x(np) = xbs(i) + w * vnx(i)

            y(np) = ybs(i) + w * vny(i)

         end do

      end subroutine get_north_boundary


      !> \brief Calculates the node distribution based on a geometric progression distribution
      subroutine get_gp_distribution(iiv, ifv, xiv, xfv, a1, xb) ! Output: last one
         implicit none
         integer, intent(in)  :: iiv !< Initial value of i
         integer, intent(in)  :: ifv !< Final value of i
         real(8), intent(in)  :: xiv !< Initial value of x
         real(8), intent(in)  :: xfv !< Final value of x
         real(8), intent(in)  :: a1  !< Width of the first volume
         real(8), intent(out) :: xb(iiv:ifv) !< Distribution of points

         ! Inner variables

         integer :: i
         integer :: n
         real(8) :: r
         real(8) :: q


         n = ifv - iiv

         r = ( xfv - xiv ) / a1

         ! Calculating q
         call get_gp_ratio(n, r, q)

         xb(iiv) = xiv

         do i = iiv+1, ifv-1

            xb(i) = xb(i-1) + a1 * q ** ( i - iiv - 1 )

         end do

         xb(ifv) = xfv

      end subroutine get_gp_distribution


      !> \brief Calculates the base q of the geometric progression
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

           if ( 0.d0 < f(n, r, qi) * f(n, r, qm) ) then

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

      end subroutine get_gp_ratio

      !> \brief Function used in the subroutine 'get_gp_ratio'
      real(8) function f(n, r, q)
         implicit none
         integer, intent(in)  ::   n !< number of partitions
         real(8), intent(in)  ::   r !< l/a1
         real(8), intent(in)  ::   q !< q

         f = q ** n + r * ( 1.d0 - q ) - 1.d0

      end function f

      !> \brief Given a reference distribution of points (xbr,ybr)
      !! a new one (xb,yb) is created based on a geometric progression distributiuon
      !! based on the curve length.
      subroutine get_gp_length_distribution( nbr, nb, al, ar, xbr, ybr, xb, yb) ! Output: last two
         implicit none
         integer, intent(in)  :: nbr        !< Number of partitions over the reference curve (xbr,ybr)
         integer, intent(in)  :: nb         !< Number of partitions for the new distribution of points
         real(8), intent(in)  :: al         !< Length of the first partition on the left
         real(8), intent(in)  :: ar         !< Length of the first partition on the right
         real(8), intent(in)  :: xbr(0:nbr) !< Given distribution of points
         real(8), intent(in)  :: ybr(0:nbr) !< Given distribution of points
         real(8), intent(out) :: xb(0:nb)   !< New distribution of points
         real(8), intent(out) :: yb(0:nb)   !< New distribution of points

         ! Inner variables

         real(8) :: sbr(0:nbr) ! Length of each partition over the reference curve
         real(8) :: sb(0:nb)   ! Length of each partition over the new curve
         real(8) :: L          ! Total length


         ! Calculates the length of the curve
         call get_sbr(nbr, xbr, ybr, sbr) ! Output: last one


         ! Total length
         L = sbr(nbr)


         ! Selecting the kind of distribution

         if ( al < 0.d0 .and. ar < 0.d0 ) then

            ! Calculates the new points according to a uniform distribution
            call get_uniform_grid( nb, L, sb) ! Output: last one

         else if ( al > 0.d0 .and. ar < 0.d0 ) then

            ! Calculates the distribution of points according to a geometric
            ! progression. The length of the first partition is defined on the left.
            call get_gp_grid_left( nb, al, L, sb)

         else if ( al < 0.d0 .and. ar > 0.d0 ) then

            ! Calculates the distribution of points according to a geometric
            ! progression. The length of the first partition is defined on the right.
            call get_gp_grid_right( nb, ar, L, sb)

         else

            ! Calculates the distribution of points according to a geometric
            ! progression. The length of the first partition is defined on the left
            ! and on the right.
            call get_gp_grid_left_right(nb, al, ar, L, sb)

         end if


         ! Calculates the new distribution of points
         call get_new_distribution(nbr, nb, sbr, xbr, ybr, sb, xb, yb) ! Output: last two


      end subroutine get_gp_length_distribution

      !> \brief Calculates the length of the curve
      subroutine get_sbr(nbr, xbr, ybr, sbr) ! Output: last one
         implicit none
         integer, intent(in)  :: nbr        !< Number of partitions over the reference curve (xbr,ybr)
         real(8), intent(in)  :: xbr(0:nbr) !< Given distribution of points
         real(8), intent(in)  :: ybr(0:nbr) !< Given distribution of points
         real(8), intent(out) :: sbr(0:nbr) !< Length of the reference curve

         ! Inner variables

         integer :: i


         sbr(0) = 0.d0

         do i = 1, nbr

            sbr(i) = sbr(i-1) + sqrt( ( xbr(i) - xbr(i-1) )**2 + ( ybr(i) - ybr(i-1) )**2 )

         end do


      end subroutine get_sbr


      !> \brief Given the length distribution, the (xb,yb) is calculated
      subroutine get_new_distribution(nbr, nb, sbr, xbr, ybr, sb, xb, yb) ! Output: last two
         integer, intent(in)  :: nbr        !< Number of partitions over the reference curve (xbr,ybr)
         integer, intent(in)  :: nb         !< Number of partitions for the new distribution of points
         real(8), intent(in)  :: sbr(0:nbr) !< Length of the curve
         real(8), intent(in)  :: xbr(0:nbr) !< Given distribution of points
         real(8), intent(in)  :: ybr(0:nbr) !< Given distribution of points
         real(8), intent(in)  :: sb(0:nb)   !< New length distribution over the curve
         real(8), intent(out) :: xb(0:nb)   !< New distribution of points
         real(8), intent(out) :: yb(0:nb)   !< New distribution of points


         ! Inner variables

         integer :: i, j

         xb(0) = xbr(0)
         yb(0) = ybr(0)

         j = 0

         do i = 1, nb-1

            do j = j, nbr

               if( sbr(j) - sb(i) >= 0.d0 ) exit

            end do

            ! Linear interpolation

            xb(i) = xbr(j-1) + ( xbr(j)-xbr(j-1) ) * ( sb(i) - sbr(j-1) ) / ( sbr(j) - sbr(j-1) )
            yb(i) = ybr(j-1) + ( ybr(j)-ybr(j-1) ) * ( sb(i) - sbr(j-1) ) / ( sbr(j) - sbr(j-1) )

         end do

         xb(nb) = xbr(nbr)
         yb(nb) = ybr(nbr)

      end subroutine get_new_distribution



      !> \brief Calculates the uniform distribution of points
      subroutine get_uniform_grid( nb, L, sb) ! Output: last one
         implicit none
         integer, intent(in)  :: nb       !< Number of partitions over the new distribution of points
         real(8), intent(in)  :: L        !< Length of the curve
         real(8), intent(out) :: sb(0:nb) !< Length distribution over the curve

         ! Inner variables

         integer :: i

         real(8) :: ds

         ds = L / dble(nb)


         sb(0) = 0.d0

         do i = 1, nb

            sb(i) = sb(i-1) + ds

         end do

      end subroutine get_uniform_grid



      !> \brief Calculates the distribution of points according
      !! to a geometric progression. The length of the first partition
      !! is defined on the left.
      subroutine get_gp_grid_left(nb, al, L, sb)
         implicit none
         integer, intent(in)  :: nb       !< Number of partitions over the new distribution of points
         real(8), intent(in)  :: al       !< Length of the first patition on the left
         real(8), intent(in)  :: L        !< Length of the curve
         real(8), intent(out) :: sb(0:nb) !< Length distribution over the curve

         ! Inner variables

         integer :: i
         real(8) :: q

         ! Calculates the ratio q of the geometric progression
         call get_gp_ratio(nb, L/al, q)

         sb(0) = 0.d0

         do i = 1, nb-1

            sb(i) = sb(i-1) + al * q ** (i-1)

         end do

         sb(nb) = L

      end subroutine get_gp_grid_left


      !> \brief Calculates the distribution of points according
      !! to a geometric progression. The length of the first partition
      !! is defined on the right.
      subroutine get_gp_grid_right(nb, ar, L, sb)
         implicit none
         integer, intent(in)  :: nb       !< Number of partitions over the new distribution of points
         real(8), intent(in)  :: ar       !< Length of the first patition on the right
         real(8), intent(in)  :: L        !< Length of the curve
         real(8), intent(out) :: sb(0:nb) !< Length distribution over the curve

         ! Inner variables

         integer :: i
         real(8) :: q

         ! Calculates the ratio q of the geometric progression
         call get_gp_ratio(nb, L/ar, q)

         sb(0) = 0.d0

         do i = 1, nb-1

            sb(i) = sb(i-1) + ar * q ** (nb-i)

         end do

         sb(nb) = L

      end subroutine get_gp_grid_right



      !> \brief Calculates the distribution of points according
      !! to a geometric progression. The length of the first partition
      !! is defined on the left and on the right.
      subroutine get_gp_grid_left_right(nb, al, ar, L, sb)
         implicit none
         integer, intent(in)  :: nb       !< Number of partitions over the new distribution of points
         real(8), intent(in)  :: al       !< Length of the first patition on the left
         real(8), intent(in)  :: ar       !< Length of the first patition on the right
         real(8), intent(in)  :: L        !< Length of the curve
         real(8), intent(out) :: sb(0:nb) !< Length distribution over the curve

         ! Inner parameters

         integer, parameter :: nitmax = 10000

         ! Inner variables

         integer :: i
         integer :: nl
         integer :: nr
         real(8) :: Lm
         real(8) :: Li
         real(8) :: Lf
         real(8) :: acl
         real(8) :: acr


         nl = nb / 2

         nr = nb - nl

         Li = 0.d0

         Lf = L

         do i = 1, nitmax

            Lm = 0.5d0 * ( Li + Lf )

            call get_gp_grid_left (nl, al, Lm, sb( 0:nl))

            call get_gp_grid_right(nr, ar, L-Lm, sb(nl:nb))

            sb(nl:nb) = sb(nl:nb) + Lm

            acl = sb(nl) - sb(nl-1)

            acr = sb(nl+1) - sb(nl)

            if ( acl > acr ) then

               Lf = Lm

            else

               Li = Lm

            end if

            if ( abs(acl-acr) < 1.d-14 ) exit

         end do

      end subroutine get_gp_grid_left_right


      !> \brief Calculates the von Karman optimized shape
      subroutine get_karman_shape(nr, rb, lo, xr, yr)
         implicit none
         integer, intent(in)  :: nr !< Number of partitions
         real(8), intent(in)  :: rb !< Base radius
         real(8), intent(in)  :: lo !< Ogive length
         real(8), intent(out) :: xr(0:nr) !< x coordinate
         real(8), intent(out) :: yr(0:nr) !< y coordinate

         ! Inner variables

         integer :: i
         real(8) :: z
         real(8) :: pi


         pi = acos(-1.d0)

         xr(0) = 0.d0
         yr(0) = 0.d0

         do i = 1, nr-1

            z = ( dble(i) / dble(nr) ) ** 1.5d0

            xr(i) = lo * z

            yr(i) = ( asin(sqrt(z))-(1.d0-2.d0*z)*sqrt(z*(1.d0-z)) ) / (2.d0*pi)

            yr(i) = 2.d0 * rb * sqrt( yr(i) )

         end do

         xr(nr) = lo

         yr(nr) = rb

      end subroutine get_karman_shape

   end subroutine get_grid_boundary_g15


   !> \brief Calculates the coordinates (x,y) of the north and south boundaries.
   !! The north boundary line is generated from the south boundary.
   !! The south boundary line is an bluff nose followed by a shifted power law ogive of
   !! length lo and radius rb followed by a cylinder of radius rb and length lc.
   !! The distribution of points over the south boundary is separated in three
   !! parts:
   !! Part   I: geometric progression (bluff nose+shifted power law)
   !! Part  II: geometric progression (shifted power law)
   !! Part III: geometric progression (cylinder)
   !! The distribution of points over the north boundary is  based on the normal
   !! vectors from the south boundary. The directions of the normal vectors
   !! may be changed by a smoothing procedure.
   !! Adapted to DEPP optimizer.
   subroutine get_grid_boundary_g16(nx, ny, unt, fgeom, lr, rb, x, y) ! Output: last two
      implicit none
      integer, intent(in) :: nx  !< Number of volumes in the csi direction (real+fictitious)
      integer, intent(in) :: ny  !< Number of volumes in the eta direction (real+fictitious)
      integer, intent(in) :: unt !< Unit where the input parameters will be printed
      character(len=*), intent(in) :: fgeom !< File of the geometric parameters

      real(8), intent(out) :: lr  !< Length of the body (m)
      real(8), intent(out) :: rb  !< Base radius of the body (m)
      real(8), dimension(nx*ny), intent(out) :: x !< Coord. x at the northeast corner of the volume P (m)
      real(8), dimension(nx*ny), intent(out) :: y !< Coord. y at the northeast corner of the volume P (m)

      ! Parameters

      integer, parameter :: nr = 10000 ! Number of partitions for the representing curve (xr,yr)


      ! Inner variables
      integer :: i, j, np ! Dummy indexes

      integer :: kgb ! Kind of geometry of the body
      integer :: nx1 ! Number of points in the first part of the grid along the eta lines
      integer :: nx2 ! Number of points in the first and second parts of the grid along the eta lines
      integer :: nx3 ! Number of points in the first to third parts of the grid along the eta lines
      integer :: iiv ! Initial value of i
      integer :: ifv ! Final value of i

      real(8) :: a2  ! Width of the first volume of Part II
      real(8) :: a3  ! Width of the first volume of Part III

      real(8) :: xiv ! Initial value of x
      real(8) :: xfv ! Final value of x

      real(8) :: lbd ! Exponent of the power law ogive

      real(8) :: lo  ! Ogive length (m)
      real(8) :: w   ! Distance between south and north boundaries (m)
      real(8) :: fs  ! Factor of smoothing (0<=fs) (necessary for non-smooth surfaces)
      real(8) :: lc  ! Cylinder length (m)
      real(8) :: lo1 ! Length of the first part of the ogive (m)
      real(8) :: rf  ! Radius of the nose tip (m)
      real(8) :: h   ! Ratio rf/rb
      real(8) :: aks ! Exponent of the power law for the calculation of the width of the first volume
      real(8) :: flo ! Fraction of lo
      real(8) :: fp1 ! Fraction of points in the first part of the grid
      real(8) :: fp2 ! Fraction of points in the second part of the grid


      real(8), dimension(nx) :: xbs ! Distribution of points over the south boundary
      real(8), dimension(nx) :: ybs ! Distribution of points over the south boundary

      real(8), dimension(0:nr) :: xr ! Distribution of points over the representing curve
      real(8), dimension(0:nr) :: yr ! Distribution of points over the representing curve


      ! DEPP variables
      integer, parameter :: nu = 2 ! number of unknowns
      integer            :: ind    ! number of the individual
      integer,       dimension(nu) :: xopt   ! Optimization checker
      character(10), dimension(nu) :: xname  ! Name of parameters
      real(8),       dimension(nu) :: xvalue ! parameters
      character(200)               :: sname  ! simulation name


      ! Reading parameters

      open(13,file=fgeom)

      read(13,*) kgb  ! Kind of geometry of the body
      read(13,*) w    ! Distance between south and north boundaries (m)
      read(13,*) fs   ! Factor of smoothing (0<=fs) (necessary for non-smooth surfaces)
      read(13,*) lo   ! Length of the ogive (m)
      read(13,*) lc   ! Length of the cylinder (m)
      read(13,*) rb   ! Base radius of the body (m)
      read(13,*) xname(1) ! Parameter name
      read(13,*) xopt(1)  ! Will this parameter be optimized? ( 0 = no, 1 = yes )
      read(13,*) xvalue(1)! Parameter value
      read(13,*) xname(2) ! Parameter name
      read(13,*) xopt(2)  ! Will this parameter be optimized? ( 0 = no, 1 = yes )
      read(13,*) xvalue(2)! Parameter value
      read(13,*) aks  ! Exponent of the power law for the south boundary
      read(13,*) flo  ! Fraction of lo
      read(13,*) fp1  ! Fraction of volumes in the first part
      read(13,*) fp2  ! Fraction of volumes in the second part

      close(13)

      ! Reads the parameters from DEPP
      call depp_get_parameters(nu, xopt, xname, xvalue, ind, sname)

      h = xvalue(1)

      lbd = xvalue(2)

      ! Writing parameters

      write(unt,*)
      write(unt,*) " *** Geometry parameters ***"
      write(unt,*)

      write(unt,"(I23,    A)") kgb, " = kgb: Kind of geometry of the body"
      write(unt,"(ES23.16,A)")   w,  " =   w: Distance between south and north boundaries (m)"
      write(unt,"(ES23.16,A)")  fs, " =  fs: Factor of smoothing (0<=fs) (necessary for non-smooth surfaces)"
      write(unt,"(ES23.16,A)")  lo, " =  lo: Length of the ogive (m)"
      write(unt,"(ES23.16,A)")  lc, " =  lc: Length of the cylinder (m)"
      write(unt,"(ES23.16,A)")  rb, " =  rb: Base radius of the body (m)"
      write(unt,"(ES23.16,A)")   h, " =   h: Ratio rf/rb (frontal to base radius ratio)"
      write(unt,"(ES23.16,A)") lbd, " = lbd: Exponent of the power law ogive"
      write(unt,"(ES23.16,A)") aks, " = aks: Exponent of the power law for the south boundary"
      write(unt,"(ES23.16,A)") flo, " = flo: Fraction of lo"
      write(unt,"(ES23.16,A)") fp1, " = fp1: Fraction of volumes in the first part of the grid"
      write(unt,"(ES23.16,A)") fp2, " = fp2: Fraction of volumes in the second part of the grid"

      ! Calculating the number of points in the first part of the grid

      nx1 = int( (nx-2) * fp1 )

      ! Calculating the number of points in the first and second parts of the grid

      nx2 = int( (nx-2) * fp2 ) + nx1

      ! Calculating the number of points in the 1st, 2nd and 3rd part of the grid

      nx3 = nx-1

      ! Calculating frontal radius

      rf = h * rb

      ! Calculating the fraction of lo used in the 2nd part of the grid

      lo1 = flo * lo

      ! Calculating the body length

      lr = lo + lc


      ! ========================================================

      ! Generating south boundary

      j = 1

      ! Part I


      ! Creating the representing curve

      xr(0) = 0.d0

      yr(0) = 0.d0

      xiv = 0.d0

      xfv = lo1

      do i = 1, nr

         xr(i) = (xfv-xiv) * ( dble(i-1)/dble(nr-1) ) ** 2 + xiv

         yr(i) = rb * ( h**(1.d0/lbd) + xr(i)/lo * (1.d0-h**(1.d0/lbd))) ** lbd

      end do

      ! Calculates the distribution of points based on the arclength and the
      ! double exponential distribution


      iiv = 1   ! Initial value of i

      ifv = nx1 ! Final value of i

      call get_dexp_m1_ditribution( nr, ifv-iiv, aks, rf, xr, yr, xbs(iiv:ifv) &
         , ybs(iiv:ifv)) ! Output: last two

      do i = iiv, ifv

         np   = nx * (j-1) + i

         x(np) = xbs(i)

         y(np) = ybs(i)

      end do


      ! Part II

      iiv = nx1  ! Initial value of i

      ifv = nx2 ! Final value of i

      xiv = lo1 ! Initial value of x

      xfv = lo  ! Final value of x

      ! Calculating the width of the first volume of Part II

      i = iiv

      a2 = sqrt( (xbs(i)-xbs(i-1))**2 + (ybs(i)-ybs(i-1))**2 )

      ! Creating the representing curve

      do i = 0, nr

         xr(i) = (xfv-xiv) * ( dble(i)/dble(nr) ) + xiv

         yr(i) = rb * ( h**(1.d0/lbd) + xr(i)/lo * (1.d0-h**(1.d0/lbd))) ** lbd

      end do

      ! Calculating the distribution curve according to a geometric progression
      ! distribution based on the curve length

      call get_gp_length_distribution(nr, ifv-iiv, a2, -1.d0, xr, yr &
         , xbs(iiv:ifv), ybs(iiv:ifv)) ! Output: last two

      do i = iiv, ifv

         np   = nx * (j-1) + i

         x(np) = xbs(i)

         y(np) = ybs(i)

      end do


      ! Part III

      iiv = nx2  ! Initial value of i

      ifv = nx3  ! Final value of i

      xiv = lo ! Initial value of x

      xfv = lo+lc  ! Final value of x

      ! Calculating the width of the first volume of Part IV

      i = iiv

      a3 = sqrt( (xbs(i)-xbs(i-1))**2 + (ybs(i)-ybs(i-1))**2 )

      call get_gp_distribution(iiv, ifv, xiv, xfv, a3, xbs(iiv:ifv)) ! Output: last one

      do i = iiv, ifv

         np   = nx * (j-1) + i

         ybs(i) = rb

         x(np) = xbs(i)

         y(np) = ybs(i)

      end do


      ! ========================================================

      ! Generating north boundary

      call get_north_boundary(nx, ny, w, fs, xbs, ybs, x, y) ! Output: last two

      ! ========================================================

   contains


      !> \brief Generates the north boundary of the domain based on the normal
      !! vectors from the south boundary. The directions of the normal vectors
      !! may be changed by a smoothing procedure.
      subroutine get_north_boundary(nx, ny, w, fs, xbs, ybs, x, y) ! Output: last two
         implicit none
         integer, intent(in)    :: nx  !< Number of volumes in the csi direction (real+fictitious)
         integer, intent(in)    :: ny  !< Number of volumes in the eta direction (real+fictitious)
         real(8), intent(in)    :: w   !< Width of the domain
         real(8), intent(in)    :: fs  !< Factor of smoothing (0<=fs)
         real(8), dimension(nx),    intent(in)  :: xbs !< Distribution of points over the south boundary
         real(8), dimension(nx),    intent(in)  :: ybs !< Distribution of points over the south boundary
         real(8), dimension(nx*ny), intent(out) :: x   !< Coord. x at the northeast corner of the volume P (m)
         real(8), dimension(nx*ny), intent(out) :: y   !< Coord. y at the northeast corner of the volume P (m)

         ! Inner variables

         integer :: i, j, k, np, ns
         real(8) :: raux
         real(8) :: tx
         real(8) :: ty
         real(8) :: tm
         real(8) :: vnx(nx-1)
         real(8) :: vny(nx-1)

         ! Number of iteractions for smoothing the direction vectors
         ns = int(nx * fs)

         ! Calculating the initial normal directions

         ! Direction of the west boundary
         vnx(1) = -1.d0
         vny(1) =  0.d0


         ! Direction of the east boundary
         vnx(nx-1) = 0.d0
         vny(nx-1) = 1.d0

         ! Directions normal to the body surface
         do i = 2, nx-2

            tx = ( xbs(i+1) - xbs(i-1) ) / 2.d0
            ty = ( ybs(i+1) - ybs(i-1) ) / 2.d0

            tm = sqrt( tx * tx + ty * ty )

            vnx(i) = - ty / tm
            vny(i) =   tx / tm

         end do

         ! Calculating the final directions (smoothing)

         do k = 1, ns

            do i = 2, nx-2

               vnx(i) = ( vnx(i+1) + vnx(i-1) ) / 2.d0

               vny(i) = ( vny(i+1) + vny(i-1) ) / 2.d0

               raux = dsqrt( vnx(i) * vnx(i) + vny(i) * vny(i) )

               vnx(i) = vnx(i) / raux

               vny(i) = vny(i) / raux

            end do

         end do

         ! Generating the north boundary

         j = ny-1

         do i = 1, nx-1

            np   = nx * (j-1) + i

            x(np) = xbs(i) + w * vnx(i)

            y(np) = ybs(i) + w * vny(i)

         end do

      end subroutine get_north_boundary


      !> \brief Calculates the node distribution based on a geometric progression distribution
      subroutine get_gp_distribution(iiv, ifv, xiv, xfv, a1, xb) ! Output: last one
         implicit none
         integer, intent(in)  :: iiv !< Initial value of i
         integer, intent(in)  :: ifv !< Final value of i
         real(8), intent(in)  :: xiv !< Initial value of x
         real(8), intent(in)  :: xfv !< Final value of x
         real(8), intent(in)  :: a1  !< Width of the first volume
         real(8), intent(out) :: xb(iiv:ifv) !< Distribution of points

         ! Inner variables

         integer :: i
         integer :: n
         real(8) :: r
         real(8) :: q


         n = ifv - iiv

         r = ( xfv - xiv ) / a1

         ! Calculating q
         call get_gp_ratio(n, r, q)

         xb(iiv) = xiv

         do i = iiv+1, ifv-1

            xb(i) = xb(i-1) + a1 * q ** ( i - iiv - 1 )

         end do

         xb(ifv) = xfv

      end subroutine get_gp_distribution


      !> \brief Calculates the base q of the geometric progression
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

           if ( 0.d0 < f(n, r, qi) * f(n, r, qm) ) then

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

      end subroutine get_gp_ratio

      !> \brief Function used in the subroutine 'get_gp_ratio'
      real(8) function f(n, r, q)
         implicit none
         integer, intent(in)  ::   n !< number of partitions
         real(8), intent(in)  ::   r !< l/a1
         real(8), intent(in)  ::   q !< q

         f = q ** n + r * ( 1.d0 - q ) - 1.d0

      end function f

      !> \brief Given a reference distribution of points (xbr,ybr)
      !! a new one (xb,yb) is created based on a geometric progression distributiuon
      !! based on the curve length.
      subroutine get_gp_length_distribution( nbr, nb, al, ar, xbr, ybr, xb, yb) ! Output: last two
         implicit none
         integer, intent(in)  :: nbr        !< Number of partitions over the reference curve (xbr,ybr)
         integer, intent(in)  :: nb         !< Number of partitions for the new distribution of points
         real(8), intent(in)  :: al         !< Length of the first partition on the left
         real(8), intent(in)  :: ar         !< Length of the first partition on the right
         real(8), intent(in)  :: xbr(0:nbr) !< Given distribution of points
         real(8), intent(in)  :: ybr(0:nbr) !< Given distribution of points
         real(8), intent(out) :: xb(0:nb)   !< New distribution of points
         real(8), intent(out) :: yb(0:nb)   !< New distribution of points

         ! Inner variables

         real(8) :: sbr(0:nbr) ! Length of each partition over the reference curve
         real(8) :: sb(0:nb)   ! Length of each partition over the new curve
         real(8) :: L          ! Total length


         ! Calculates the length of the curve
         call get_sbr(nbr, xbr, ybr, sbr) ! Output: last one


         ! Total length
         L = sbr(nbr)


         ! Selecting the kind of distribution

         if ( al < 0.d0 .and. ar < 0.d0 ) then

            ! Calculates the new points according to a uniform distribution
            call get_uniform_grid( nb, L, sb) ! Output: last one

         else if ( al > 0.d0 .and. ar < 0.d0 ) then

            ! Calculates the distribution of points according to a geometric
            ! progression. The length of the first partition is defined on the left.
            call get_gp_grid_left( nb, al, L, sb)

         else if ( al < 0.d0 .and. ar > 0.d0 ) then

            ! Calculates the distribution of points according to a geometric
            ! progression. The length of the first partition is defined on the right.
            call get_gp_grid_right( nb, ar, L, sb)

         else

            ! Calculates the distribution of points according to a geometric
            ! progression. The length of the first partition is defined on the left
            ! and on the right.
            call get_gp_grid_left_right(nb, al, ar, L, sb)

         end if


         ! Calculates the new distribution of points
         call get_new_distribution(nbr, nb, sbr, xbr, ybr, sb, xb, yb) ! Output: last two


      end subroutine get_gp_length_distribution

      !> \brief Calculates the length of the curve
      subroutine get_sbr(nbr, xbr, ybr, sbr) ! Output: last one
         implicit none
         integer, intent(in)  :: nbr        !< Number of partitions over the reference curve (xbr,ybr)
         real(8), intent(in)  :: xbr(0:nbr) !< Given distribution of points
         real(8), intent(in)  :: ybr(0:nbr) !< Given distribution of points
         real(8), intent(out) :: sbr(0:nbr) !< Length of the reference curve

         ! Inner variables

         integer :: i


         sbr(0) = 0.d0

         do i = 1, nbr

            sbr(i) = sbr(i-1) + sqrt( ( xbr(i) - xbr(i-1) )**2 + ( ybr(i) - ybr(i-1) )**2 )

         end do


      end subroutine get_sbr


      !> \brief Given the length distribution, the (xb,yb) is calculated
      subroutine get_new_distribution(nbr, nb, sbr, xbr, ybr, sb, xb, yb) ! Output: last two
         integer, intent(in)  :: nbr        !< Number of partitions over the reference curve (xbr,ybr)
         integer, intent(in)  :: nb         !< Number of partitions for the new distribution of points
         real(8), intent(in)  :: sbr(0:nbr) !< Length of the curve
         real(8), intent(in)  :: xbr(0:nbr) !< Given distribution of points
         real(8), intent(in)  :: ybr(0:nbr) !< Given distribution of points
         real(8), intent(in)  :: sb(0:nb)   !< New length distribution over the curve
         real(8), intent(out) :: xb(0:nb)   !< New distribution of points
         real(8), intent(out) :: yb(0:nb)   !< New distribution of points


         ! Inner variables

         integer :: i, j

         xb(0) = xbr(0)
         yb(0) = ybr(0)

         j = 0

         do i = 1, nb-1

            do j = j, nbr

               if( sbr(j) - sb(i) >= 0.d0 ) exit

            end do

            ! Linear interpolation

            xb(i) = xbr(j-1) + ( xbr(j)-xbr(j-1) ) * ( sb(i) - sbr(j-1) ) / ( sbr(j) - sbr(j-1) )
            yb(i) = ybr(j-1) + ( ybr(j)-ybr(j-1) ) * ( sb(i) - sbr(j-1) ) / ( sbr(j) - sbr(j-1) )

         end do

         xb(nb) = xbr(nbr)
         yb(nb) = ybr(nbr)

      end subroutine get_new_distribution



      !> \brief Calculates the uniform distribution of points
      subroutine get_uniform_grid( nb, L, sb) ! Output: last one
         implicit none
         integer, intent(in)  :: nb       !< Number of partitions over the new distribution of points
         real(8), intent(in)  :: L        !< Length of the curve
         real(8), intent(out) :: sb(0:nb) !< Length distribution over the curve

         ! Inner variables

         integer :: i

         real(8) :: ds

         ds = L / dble(nb)


         sb(0) = 0.d0

         do i = 1, nb

            sb(i) = sb(i-1) + ds

         end do

      end subroutine get_uniform_grid



      !> \brief Calculates the distribution of points according
      !! to a geometric progression. The length of the first partition
      !! is defined on the left.
      subroutine get_gp_grid_left(nb, al, L, sb)
         implicit none
         integer, intent(in)  :: nb       !< Number of partitions over the new distribution of points
         real(8), intent(in)  :: al       !< Length of the first patition on the left
         real(8), intent(in)  :: L        !< Length of the curve
         real(8), intent(out) :: sb(0:nb) !< Length distribution over the curve

         ! Inner variables

         integer :: i
         real(8) :: q

         ! Calculates the ratio q of the geometric progression
         call get_gp_ratio(nb, L/al, q)

         sb(0) = 0.d0

         do i = 1, nb-1

            sb(i) = sb(i-1) + al * q ** (i-1)

         end do

         sb(nb) = L

      end subroutine get_gp_grid_left


      !> \brief Calculates the distribution of points according
      !! to a geometric progression. The length of the first partition
      !! is defined on the right.
      subroutine get_gp_grid_right(nb, ar, L, sb)
         implicit none
         integer, intent(in)  :: nb       !< Number of partitions over the new distribution of points
         real(8), intent(in)  :: ar       !< Length of the first patition on the right
         real(8), intent(in)  :: L        !< Length of the curve
         real(8), intent(out) :: sb(0:nb) !< Length distribution over the curve

         ! Inner variables

         integer :: i
         real(8) :: q

         ! Calculates the ratio q of the geometric progression
         call get_gp_ratio(nb, L/ar, q)

         sb(0) = 0.d0

         do i = 1, nb-1

            sb(i) = sb(i-1) + ar * q ** (nb-i)

         end do

         sb(nb) = L

      end subroutine get_gp_grid_right



      !> \brief Calculates the distribution of points according
      !! to a geometric progression. The length of the first partition
      !! is defined on the left and on the right.
      subroutine get_gp_grid_left_right(nb, al, ar, L, sb)
         implicit none
         integer, intent(in)  :: nb       !< Number of partitions over the new distribution of points
         real(8), intent(in)  :: al       !< Length of the first patition on the left
         real(8), intent(in)  :: ar       !< Length of the first patition on the right
         real(8), intent(in)  :: L        !< Length of the curve
         real(8), intent(out) :: sb(0:nb) !< Length distribution over the curve

         ! Inner parameters

         integer, parameter :: nitmax = 10000

         ! Inner variables

         integer :: i
         integer :: nl
         integer :: nr
         real(8) :: Lm
         real(8) :: Li
         real(8) :: Lf
         real(8) :: acl
         real(8) :: acr


         nl = nb / 2

         nr = nb - nl

         Li = 0.d0

         Lf = L

         do i = 1, nitmax

            Lm = 0.5d0 * ( Li + Lf )

            call get_gp_grid_left (nl, al, Lm, sb( 0:nl))

            call get_gp_grid_right(nr, ar, L-Lm, sb(nl:nb))

            sb(nl:nb) = sb(nl:nb) + Lm

            acl = sb(nl) - sb(nl-1)

            acr = sb(nl+1) - sb(nl)

            if ( acl > acr ) then

               Lf = Lm

            else

               Li = Lm

            end if

            if ( abs(acl-acr) < 1.d-14 ) exit

         end do

      end subroutine get_gp_grid_left_right


      !> \brief Calculates a distribution of points based on the arclength and
      !! on the double-exponential rule
      subroutine get_dexp_m1_ditribution( nbr, nb, aks, zcv, xbr, ybr, xb, yb) ! Output: last two
         implicit none
         integer, intent(in)  :: nbr        !< Number of partitions over the reference curve (xbr,ybr)
         integer, intent(in)  :: nb         !< Number of partitions for the new distribution of points
         real(8), intent(in)  :: aks        !< Exponent for the calculation of acv
         real(8), intent(in)  :: zcv        !< Central value of z
         real(8), intent(in)  :: xbr(0:nbr) !< Given distribution of points
         real(8), intent(in)  :: ybr(0:nbr) !< Given distribution of points
         real(8), intent(out) :: xb(0:nb)   !< New distribution of points
         real(8), intent(out) :: yb(0:nb)   !< New distribution of points

         ! Inner variables

         real(8) :: sbr(0:nbr) ! Length of each partition over the reference curve
         real(8) :: sb(0:nb)   ! Length of each partition over the new curve
         real(8) :: L          ! Total length
         real(8) :: acv        ! Length of the volume closer to the central point


         ! Calculates the length of the curve
         call get_sbr(nbr, xbr, ybr, sbr) ! Output: last one


         ! Total length
         L = sbr(nbr)


         ! Calculating the double-exponential distribution of arclength

         acv = L / ( nb ** aks )

         if ( acv > zcv ) then

            write(*,*) "Error. get_dexp_m1_ditribution: there is no resolution"&
               , " for this grid. Stopping."

            stop

         end if

         call get_dexp_m1(nb, 0.d0, zcv, L, acv, sb)


         ! Calculates the new distribution of points
         call get_new_distribution(nbr, nb, sbr, xbr, ybr, sb, xb, yb) ! Output: last two


      end subroutine get_dexp_m1_ditribution

      !> \brief Calculates a distribution of points based on the double-exponential
      !! rule. The coefficients of this rule are set by defining the position of
      !! an internal point (zcv) and the width (acv) of the volumes beside this
      !! point. The points are concentrated around zcv.
      subroutine get_dexp_m1(np, ziv, zcv, zfv, acv, z)
         implicit none
         integer, intent(in)  :: np  !< Number of partitions
         real(8), intent(in)  :: ziv !< Initial value of z
         real(8), intent(in)  :: zcv !< Central value of z (exponential's matching point)
         real(8), intent(in)  :: zfv !< Final value of z
         real(8), intent(in)  :: acv !< Length of the volumes beside zcv
         real(8), intent(out) :: z(0:np) !< Distribution of points

         ! Inner variables

         integer :: i
         integer :: npc
         real(8) :: t
         real(8) :: dt
         real(8) :: tc
         real(8) :: r
         real(8) :: A
         real(8) :: ac
         real(8) :: zc


         ! Initializing variables

         zc = zcv / ( zfv - ziv )

         ac = acv / ( zfv - ziv )

         dt = 1.d0 / dble(np)

         r = ac / ( 1.d0 - zc )



         ! Cheking conditions

         if ( ( ac / zc > 1.d0 ) .or. ( r > 1.d0 ) ) then

            write(*,*) "Error. get_dexp_m1: ac > zc or ac > (1-zc). Stopping."

            stop

         end if

         if ( ac > dt ) then

            write(*,*) "Error. get_dexp_m1: ac > dt. Stopping."

            stop

         end if


         ! Searching the value of A

         call adexp_m1(dt, r, zc, A)


         ! Calculating the fraction of volumes for each exponential

         tc = 1.d0 + log( 1.d0 + ( exp(-A*dt) - 1.d0 ) / r ) / A

         npc = nint(tc*np)



         ! Calculating the normallized distribution of points

         ! First exponential

         do i = 0, npc

            t = tc * dble(i) / dble(npc)

            z(i) = zc * ( exp(A*t) - 1.d0 ) / ( exp(A*tc) - 1.d0 )

         end do

         ! Second exponential

         do i = npc, np

            t = ( 1.d0-tc ) * dble(i-npc) / dble(np-npc) + tc

            z(i) = zc+(1.d0-zc)*(exp(-A*(t-tc))-1.d0)/(exp(-A*(1.d0-tc))-1.d0)

         end do

         z(0)  = 0.d0

         z(np) = 1.d0

         ! Calculating the extended distribution of points

         z = ( zfv - ziv ) * z + ziv

      end subroutine get_dexp_m1

      !> \brief Defines the transcendental function for finding A
      real(8) function fdexp_m1(dt, r, zc, A)
         implicit none
         real(8), intent(in) :: dt
         real(8), intent(in) :: r
         real(8), intent(in) :: zc
         real(8), intent(in) :: A

         ! Inner varialbles

         real(8) :: raux1, raux2, raux3

         if ( abs(A) > 1.d-6 ) then

            raux1 = exp(A) * ( 1.d0 + ( exp(-A*dt) - 1.d0 ) / r )
            raux2 = (2.d0*zc-1.d0) / zc
            raux3 = (1.d0-zc) / zc

            fdexp_m1 = -exp(-A) * raux1 * raux1 +  raux1 * raux2 + raux3

         else

            fdexp_m1 = (r*zc-r+dt)*A/(r*zc) &
               + ((r**2-2*dt**2)*zc-r**2+(2*dt-dt**2)*r)*A**2/(r**2*zc)/2.0

         end if

      end function fdexp_m1


      !> \brief Finds the value of A using the bissection method
      subroutine adexp_m1(dt, r, zc, A)
         implicit none
         real(8), intent(in)  :: dt
         real(8), intent(in)  :: r
         real(8), intent(in)  :: zc
         real(8), intent(out) :: A

         ! Parameters

         integer, parameter :: itmax = 1000

         ! Inner variables

         integer :: k
         real(8) :: Ai, Af, Am, Fi, Fm, Ff



         ! Initializing the search interval

         Ai = -1.d-18

         Af = -100.d0



         ! Checking for roots

         Fi = fdexp_m1(dt, r, zc, Ai)

         Ff = fdexp_m1(dt, r, zc, Af)

         if ( Fi * Ff > 0.d0 ) then

            write(*,*) "Error. adexp_m1: there is none or more than one roots. Stopping."

            stop

         end if


         ! Looking for a solution

         do k = 1, itmax

            Am = 0.5d0 * ( Ai + Af )

            Fi = fdexp_m1(dt, r, zc, Ai)

            Fm = fdexp_m1(dt, r, zc, Am)

            if ( Fi * Fm <= 0.d0 ) then

               Af = Am

            else

               Ai = Am

            end if

            if ( abs(Af-Ai) < 1.d-14 ) exit

         end do

         if ( k >= itmax ) then

            write(*,*) "Error. adexp_m1: number of iteractions exceeded. Stopping."

            stop

         end if

         A = 0.5d0 * ( Ai + Af )

      end subroutine adexp_m1


   end subroutine get_grid_boundary_g16


   !> \brief Calculates the coordinates (x,y) of the north and south boundaries.
   !! The north boundary line is generated from the south boundary.
   !! The south boundary line is a two adjustable coef. geometry based on Newton body of
   !! length lo and radius rb followed by a cylinder of radius rb and length lc.
   !! The distribution of points over the south boundary is separated in three
   !! parts:
   !! Part   I: geometric progression (Newton-2coef.)
   !! Part  II: geometric progression (Newton-2coef.)
   !! Part III: geometric progression (cylinder)
   !! The distribution of points over the north boundary is  based on the normal
   !! vectors from the south boundary. The directions of the normal vectors
   !! may be changed by a smoothing procedure.
   !! Adapted to DEPP optimizer.
   subroutine get_grid_boundary_g17(nx, ny, unt, fgeom, lr, rb, x, y) ! Output: last two
      implicit none
      integer, intent(in) :: nx  !< Number of volumes in the csi direction (real+fictitious)
      integer, intent(in) :: ny  !< Number of volumes in the eta direction (real+fictitious)
      integer, intent(in) :: unt !< Unit where the input parameters will be printed
      character(len=*), intent(in) :: fgeom !< File of the geometric parameters

      real(8), intent(out) :: lr  !< Length of the body (m)
      real(8), intent(out) :: rb  !< Base radius of the body (m)
      real(8), dimension(nx*ny), intent(out) :: x !< Coord. x at the northeast corner of the volume P (m)
      real(8), dimension(nx*ny), intent(out) :: y !< Coord. y at the northeast corner of the volume P (m)

      ! Parameters

      integer, parameter :: nr = 10000 ! Number of partitions for the representing curve (xr,yr)


      ! Inner variables
      integer :: i, j, np ! Dummy indexes

      integer :: kgb ! Kind of geometry of the body
      integer :: nx1 ! Number of points in the first part of the grid along the eta lines
      integer :: nx2 ! Number of points in the first and second parts of the grid along the eta lines
      integer :: nx3 ! Number of points in the first to third parts of the grid along the eta lines
      integer :: iiv ! Initial value of i
      integer :: ifv ! Final value of i
      integer :: nr1 ! Fraction of nr used to represent the Part I

      real(8) :: a2  ! Width of the first volume of Part II
      real(8) :: a3  ! Width of the first volume of Part III

      real(8) :: xiv ! Initial value of x
      real(8) :: xfv ! Final value of x

      real(8) :: zfg ! Frontal slope of the geometry model

      real(8) :: lo  ! Ogive length (m)
      real(8) :: w   ! Distance between south and north boundaries (m)
      real(8) :: fs  ! Factor of smoothing (0<=fs) (necessary for non-smooth surfaces)
      real(8) :: lc  ! Cylinder length (m)
      real(8) :: lo1 ! Length of the first part of the ogive (m)
      real(8) :: rf  ! Radius of the nose tip (m)
      real(8) :: h   ! Ratio rf/rb
      real(8) :: aks ! Exponent of the power law for the calculation of the width of the first volume
      real(8) :: flo ! Fraction of lo
      real(8) :: fp1 ! Fraction of points in the first part of the grid
      real(8) :: fp2 ! Fraction of points in the second part of the grid


      real(8), dimension(nx) :: xbs ! Distribution of points over the south boundary
      real(8), dimension(nx) :: ybs ! Distribution of points over the south boundary

      real(8), dimension(0:nr) :: xr ! Distribution of points over the representing curve
      real(8), dimension(0:nr) :: yr ! Distribution of points over the representing curve


      ! DEPP variables
      integer, parameter :: nu = 2 ! number of unknowns
      integer            :: ind    ! number of the individual
      integer,       dimension(nu) :: xopt   ! Optimization checker
      character(10), dimension(nu) :: xname  ! Name of parameters
      real(8),       dimension(nu) :: xvalue ! parameters
      character(200)               :: sname  ! simulation name


      ! Reading parameters

      open(13,file=fgeom)

      read(13,*) kgb  ! Kind of geometry of the body
      read(13,*) w    ! Distance between south and north boundaries (m)
      read(13,*) fs   ! Factor of smoothing (0<=fs) (necessary for non-smooth surfaces)
      read(13,*) lo   ! Length of the ogive (m)
      read(13,*) lc   ! Length of the cylinder (m)
      read(13,*) rb   ! Base radius of the body (m)
      read(13,*) xname(1) ! Parameter name
      read(13,*) xopt(1)  ! Will this parameter be optimized? ( 0 = no, 1 = yes )
      read(13,*) xvalue(1)! Parameter value
      read(13,*) xname(2) ! Parameter name
      read(13,*) xopt(2)  ! Will this parameter be optimized? ( 0 = no, 1 = yes )
      read(13,*) xvalue(2)! Parameter value
      read(13,*) aks  ! Exponent of the power law for the south boundary
      read(13,*) flo  ! Fraction of lo
      read(13,*) fp1  ! Fraction of volumes in the first part
      read(13,*) fp2  ! Fraction of volumes in the second part

      close(13)

      ! Reads the parameters from DEPP
      call depp_get_parameters(nu, xopt, xname, xvalue, ind, sname)

      h = xvalue(1)

      zfg = xvalue(2)


      ! Generating the reference curve of the Newton-2 coefficients model

      call get_newton2coef_body(nr, rb, lo, h, zfg, xr, yr) ! Output: last two


      ! Writing parameters

      write(unt,*)
      write(unt,*) " *** Geometry parameters ***"
      write(unt,*)

      write(unt,"(I23,    A)") kgb, " = kgb: Kind of geometry of the body"
      write(unt,"(ES23.16,A)")   w,  " =   w: Distance between south and north boundaries (m)"
      write(unt,"(ES23.16,A)")  fs, " =  fs: Factor of smoothing (0<=fs) (necessary for non-smooth surfaces)"
      write(unt,"(ES23.16,A)")  lo, " =  lo: Length of the ogive (m)"
      write(unt,"(ES23.16,A)")  lc, " =  lc: Length of the cylinder (m)"
      write(unt,"(ES23.16,A)")  rb, " =  rb: Base radius of the body (m)"
      write(unt,"(ES23.16,A)")   h, " =   h: Ratio rf/rb (frontal to base radius ratio)"
      write(unt,"(ES23.16,A)") zfg, " = zfg: Frontal slope of the geometry model"
      write(unt,"(ES23.16,A)") aks, " = aks: Exponent of the power law for the south boundary"
      write(unt,"(ES23.16,A)") flo, " = flo: Fraction of lo"
      write(unt,"(ES23.16,A)") fp1, " = fp1: Fraction of volumes in the first part of the grid"
      write(unt,"(ES23.16,A)") fp2, " = fp2: Fraction of volumes in the second part of the grid"

      ! Calculating the number of points in the first part of the grid

      nx1 = int( (nx-2) * fp1 )

      ! Calculating the number of points in the first and second parts of the grid

      nx2 = int( (nx-2) * fp2 ) + nx1

      ! Calculating the number of points in the 1st, 2nd and 3rd part of the grid

      nx3 = nx-1

      ! Calculating frontal radius

      rf = h * rb

      ! Calculating the fraction of lo used in the 2nd part of the grid

      lo1 = flo * lo

      do i = 1, nr

         if ( lo1 < xr(i) ) exit

      end do

      nr1 = i

      lo1 = xr(nr1)


      ! Calculating the body length

      lr = lo + lc


      ! ========================================================

      ! Generating south boundary

      j = 1

      ! Part I

      ! Calculates the distribution of points based on the arclength and the
      ! double exponential distribution

      iiv = 1   ! Initial value of i

      ifv = nx1 ! Final value of i

      call get_dexp_m1_ditribution( nr1, ifv-iiv, aks, rf, xr(0:nr1), yr(0:nr1)&
         , xbs(iiv:ifv), ybs(iiv:ifv)) ! Output: last two

      do i = iiv, ifv

         np   = nx * (j-1) + i

         x(np) = xbs(i)

         y(np) = ybs(i)

      end do


      ! Part II

      iiv = nx1  ! Initial value of i

      ifv = nx2 ! Final value of i

      xiv = lo1 ! Initial value of x

      xfv = lo  ! Final value of x

      ! Calculating the width of the first volume of Part II

      i = iiv

      a2 = sqrt( (xbs(i)-xbs(i-1))**2 + (ybs(i)-ybs(i-1))**2 )


      ! Calculating the distribution curve according to a geometric progression
      ! distribution based on the curve length

      call get_gp_length_distribution(nr-nr1, ifv-iiv, a2, -1.d0, xr(nr1:nr) &
         , yr(nr1:nr), xbs(iiv:ifv), ybs(iiv:ifv)) ! Output: last two

      do i = iiv, ifv

         np   = nx * (j-1) + i

         x(np) = xbs(i)

         y(np) = ybs(i)

      end do


      ! Part III

      iiv = nx2  ! Initial value of i

      ifv = nx3  ! Final value of i

      xiv = lo ! Initial value of x

      xfv = lo+lc  ! Final value of x

      ! Calculating the width of the first volume of Part III

      i = iiv

      a3 = sqrt( (xbs(i)-xbs(i-1))**2 + (ybs(i)-ybs(i-1))**2 )

      call get_gp_distribution(iiv, ifv, xiv, xfv, a3, xbs(iiv:ifv)) ! Output: last one

      do i = iiv, ifv

         np   = nx * (j-1) + i

         ybs(i) = rb

         x(np) = xbs(i)

         y(np) = ybs(i)

      end do


      ! ========================================================

      ! Generating north boundary

      call get_north_boundary(nx, ny, w, fs, xbs, ybs, x, y) ! Output: last two

      ! ========================================================

   contains


      !> \brief Generates the north boundary of the domain based on the normal
      !! vectors from the south boundary. The directions of the normal vectors
      !! may be changed by a smoothing procedure.
      subroutine get_north_boundary(nx, ny, w, fs, xbs, ybs, x, y) ! Output: last two
         implicit none
         integer, intent(in)    :: nx  !< Number of volumes in the csi direction (real+fictitious)
         integer, intent(in)    :: ny  !< Number of volumes in the eta direction (real+fictitious)
         real(8), intent(in)    :: w   !< Width of the domain
         real(8), intent(in)    :: fs  !< Factor of smoothing (0<=fs)
         real(8), dimension(nx),    intent(in)  :: xbs !< Distribution of points over the south boundary
         real(8), dimension(nx),    intent(in)  :: ybs !< Distribution of points over the south boundary
         real(8), dimension(nx*ny), intent(out) :: x   !< Coord. x at the northeast corner of the volume P (m)
         real(8), dimension(nx*ny), intent(out) :: y   !< Coord. y at the northeast corner of the volume P (m)

         ! Inner variables

         integer :: i, j, k, np, ns
         real(8) :: raux
         real(8) :: tx
         real(8) :: ty
         real(8) :: tm
         real(8) :: vnx(nx-1)
         real(8) :: vny(nx-1)

         ! Number of iteractions for smoothing the direction vectors
         ns = int(nx * fs)

         ! Calculating the initial normal directions

         ! Direction of the west boundary
         vnx(1) = -1.d0
         vny(1) =  0.d0


         ! Direction of the east boundary
         vnx(nx-1) = 0.d0
         vny(nx-1) = 1.d0

         ! Directions normal to the body surface
         do i = 2, nx-2

            tx = ( xbs(i+1) - xbs(i-1) ) / 2.d0
            ty = ( ybs(i+1) - ybs(i-1) ) / 2.d0

            tm = sqrt( tx * tx + ty * ty )

            vnx(i) = - ty / tm
            vny(i) =   tx / tm

         end do

         ! Calculating the final directions (smoothing)

         do k = 1, ns

            do i = 2, nx-2

               vnx(i) = ( vnx(i+1) + vnx(i-1) ) / 2.d0

               vny(i) = ( vny(i+1) + vny(i-1) ) / 2.d0

               raux = dsqrt( vnx(i) * vnx(i) + vny(i) * vny(i) )

               vnx(i) = vnx(i) / raux

               vny(i) = vny(i) / raux

            end do

         end do

         ! Generating the north boundary

         j = ny-1

         do i = 1, nx-1

            np   = nx * (j-1) + i

            x(np) = xbs(i) + w * vnx(i)

            y(np) = ybs(i) + w * vny(i)

         end do

      end subroutine get_north_boundary


      !> \brief Calculates the node distribution based on a geometric progression distribution
      subroutine get_gp_distribution(iiv, ifv, xiv, xfv, a1, xb) ! Output: last one
         implicit none
         integer, intent(in)  :: iiv !< Initial value of i
         integer, intent(in)  :: ifv !< Final value of i
         real(8), intent(in)  :: xiv !< Initial value of x
         real(8), intent(in)  :: xfv !< Final value of x
         real(8), intent(in)  :: a1  !< Width of the first volume
         real(8), intent(out) :: xb(iiv:ifv) !< Distribution of points

         ! Inner variables

         integer :: i
         integer :: n
         real(8) :: r
         real(8) :: q


         n = ifv - iiv

         r = ( xfv - xiv ) / a1

         ! Calculating q
         call get_gp_ratio(n, r, q)

         xb(iiv) = xiv

         do i = iiv+1, ifv-1

            xb(i) = xb(i-1) + a1 * q ** ( i - iiv - 1 )

         end do

         xb(ifv) = xfv

      end subroutine get_gp_distribution


      !> \brief Calculates the base q of the geometric progression
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

           if ( 0.d0 < f(n, r, qi) * f(n, r, qm) ) then

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

      end subroutine get_gp_ratio

      !> \brief Function used in the subroutine 'get_gp_ratio'
      real(8) function f(n, r, q)
         implicit none
         integer, intent(in)  ::   n !< number of partitions
         real(8), intent(in)  ::   r !< l/a1
         real(8), intent(in)  ::   q !< q

         f = q ** n + r * ( 1.d0 - q ) - 1.d0

      end function f

      !> \brief Given a reference distribution of points (xbr,ybr)
      !! a new one (xb,yb) is created based on a geometric progression distributiuon
      !! based on the curve length.
      subroutine get_gp_length_distribution( nbr, nb, al, ar, xbr, ybr, xb, yb) ! Output: last two
         implicit none
         integer, intent(in)  :: nbr        !< Number of partitions over the reference curve (xbr,ybr)
         integer, intent(in)  :: nb         !< Number of partitions for the new distribution of points
         real(8), intent(in)  :: al         !< Length of the first partition on the left
         real(8), intent(in)  :: ar         !< Length of the first partition on the right
         real(8), intent(in)  :: xbr(0:nbr) !< Given distribution of points
         real(8), intent(in)  :: ybr(0:nbr) !< Given distribution of points
         real(8), intent(out) :: xb(0:nb)   !< New distribution of points
         real(8), intent(out) :: yb(0:nb)   !< New distribution of points

         ! Inner variables

         real(8) :: sbr(0:nbr) ! Length of each partition over the reference curve
         real(8) :: sb(0:nb)   ! Length of each partition over the new curve
         real(8) :: L          ! Total length


         ! Calculates the length of the curve
         call get_sbr(nbr, xbr, ybr, sbr) ! Output: last one


         ! Total length
         L = sbr(nbr)


         ! Selecting the kind of distribution

         if ( al < 0.d0 .and. ar < 0.d0 ) then

            ! Calculates the new points according to a uniform distribution
            call get_uniform_grid( nb, L, sb) ! Output: last one

         else if ( al > 0.d0 .and. ar < 0.d0 ) then

            ! Calculates the distribution of points according to a geometric
            ! progression. The length of the first partition is defined on the left.
            call get_gp_grid_left( nb, al, L, sb)

         else if ( al < 0.d0 .and. ar > 0.d0 ) then

            ! Calculates the distribution of points according to a geometric
            ! progression. The length of the first partition is defined on the right.
            call get_gp_grid_right( nb, ar, L, sb)

         else

            ! Calculates the distribution of points according to a geometric
            ! progression. The length of the first partition is defined on the left
            ! and on the right.
            call get_gp_grid_left_right(nb, al, ar, L, sb)

         end if


         ! Calculates the new distribution of points
         call get_new_distribution(nbr, nb, sbr, xbr, ybr, sb, xb, yb) ! Output: last two


      end subroutine get_gp_length_distribution

      !> \brief Calculates the length of the curve
      subroutine get_sbr(nbr, xbr, ybr, sbr) ! Output: last one
         implicit none
         integer, intent(in)  :: nbr        !< Number of partitions over the reference curve (xbr,ybr)
         real(8), intent(in)  :: xbr(0:nbr) !< Given distribution of points
         real(8), intent(in)  :: ybr(0:nbr) !< Given distribution of points
         real(8), intent(out) :: sbr(0:nbr) !< Length of the reference curve

         ! Inner variables

         integer :: i


         sbr(0) = 0.d0

         do i = 1, nbr

            sbr(i) = sbr(i-1) + sqrt( ( xbr(i) - xbr(i-1) )**2 + ( ybr(i) - ybr(i-1) )**2 )

         end do


      end subroutine get_sbr


      !> \brief Given the length distribution, the (xb,yb) is calculated
      subroutine get_new_distribution(nbr, nb, sbr, xbr, ybr, sb, xb, yb) ! Output: last two
         integer, intent(in)  :: nbr        !< Number of partitions over the reference curve (xbr,ybr)
         integer, intent(in)  :: nb         !< Number of partitions for the new distribution of points
         real(8), intent(in)  :: sbr(0:nbr) !< Length of the curve
         real(8), intent(in)  :: xbr(0:nbr) !< Given distribution of points
         real(8), intent(in)  :: ybr(0:nbr) !< Given distribution of points
         real(8), intent(in)  :: sb(0:nb)   !< New length distribution over the curve
         real(8), intent(out) :: xb(0:nb)   !< New distribution of points
         real(8), intent(out) :: yb(0:nb)   !< New distribution of points


         ! Inner variables

         integer :: i, j

         xb(0) = xbr(0)
         yb(0) = ybr(0)

         j = 0

         do i = 1, nb-1

            do j = j, nbr

               if( sbr(j) - sb(i) >= 0.d0 ) exit

            end do

            ! Linear interpolation

            xb(i) = xbr(j-1) + ( xbr(j)-xbr(j-1) ) * ( sb(i) - sbr(j-1) ) / ( sbr(j) - sbr(j-1) )
            yb(i) = ybr(j-1) + ( ybr(j)-ybr(j-1) ) * ( sb(i) - sbr(j-1) ) / ( sbr(j) - sbr(j-1) )

         end do

         xb(nb) = xbr(nbr)
         yb(nb) = ybr(nbr)

      end subroutine get_new_distribution



      !> \brief Calculates the uniform distribution of points
      subroutine get_uniform_grid( nb, L, sb) ! Output: last one
         implicit none
         integer, intent(in)  :: nb       !< Number of partitions over the new distribution of points
         real(8), intent(in)  :: L        !< Length of the curve
         real(8), intent(out) :: sb(0:nb) !< Length distribution over the curve

         ! Inner variables

         integer :: i

         real(8) :: ds

         ds = L / dble(nb)


         sb(0) = 0.d0

         do i = 1, nb

            sb(i) = sb(i-1) + ds

         end do

      end subroutine get_uniform_grid



      !> \brief Calculates the distribution of points according
      !! to a geometric progression. The length of the first partition
      !! is defined on the left.
      subroutine get_gp_grid_left(nb, al, L, sb)
         implicit none
         integer, intent(in)  :: nb       !< Number of partitions over the new distribution of points
         real(8), intent(in)  :: al       !< Length of the first patition on the left
         real(8), intent(in)  :: L        !< Length of the curve
         real(8), intent(out) :: sb(0:nb) !< Length distribution over the curve

         ! Inner variables

         integer :: i
         real(8) :: q

         ! Calculates the ratio q of the geometric progression
         call get_gp_ratio(nb, L/al, q)

         sb(0) = 0.d0

         do i = 1, nb-1

            sb(i) = sb(i-1) + al * q ** (i-1)

         end do

         sb(nb) = L

      end subroutine get_gp_grid_left


      !> \brief Calculates the distribution of points according
      !! to a geometric progression. The length of the first partition
      !! is defined on the right.
      subroutine get_gp_grid_right(nb, ar, L, sb)
         implicit none
         integer, intent(in)  :: nb       !< Number of partitions over the new distribution of points
         real(8), intent(in)  :: ar       !< Length of the first patition on the right
         real(8), intent(in)  :: L        !< Length of the curve
         real(8), intent(out) :: sb(0:nb) !< Length distribution over the curve

         ! Inner variables

         integer :: i
         real(8) :: q

         ! Calculates the ratio q of the geometric progression
         call get_gp_ratio(nb, L/ar, q)

         sb(0) = 0.d0

         do i = 1, nb-1

            sb(i) = sb(i-1) + ar * q ** (nb-i)

         end do

         sb(nb) = L

      end subroutine get_gp_grid_right



      !> \brief Calculates the distribution of points according
      !! to a geometric progression. The length of the first partition
      !! is defined on the left and on the right.
      subroutine get_gp_grid_left_right(nb, al, ar, L, sb)
         implicit none
         integer, intent(in)  :: nb       !< Number of partitions over the new distribution of points
         real(8), intent(in)  :: al       !< Length of the first patition on the left
         real(8), intent(in)  :: ar       !< Length of the first patition on the right
         real(8), intent(in)  :: L        !< Length of the curve
         real(8), intent(out) :: sb(0:nb) !< Length distribution over the curve

         ! Inner parameters

         integer, parameter :: nitmax = 10000

         ! Inner variables

         integer :: i
         integer :: nl
         integer :: nr
         real(8) :: Lm
         real(8) :: Li
         real(8) :: Lf
         real(8) :: acl
         real(8) :: acr


         nl = nb / 2

         nr = nb - nl

         Li = 0.d0

         Lf = L

         do i = 1, nitmax

            Lm = 0.5d0 * ( Li + Lf )

            call get_gp_grid_left (nl, al, Lm, sb( 0:nl))

            call get_gp_grid_right(nr, ar, L-Lm, sb(nl:nb))

            sb(nl:nb) = sb(nl:nb) + Lm

            acl = sb(nl) - sb(nl-1)

            acr = sb(nl+1) - sb(nl)

            if ( acl > acr ) then

               Lf = Lm

            else

               Li = Lm

            end if

            if ( abs(acl-acr) < 1.d-14 ) exit

         end do

      end subroutine get_gp_grid_left_right


      !> \brief Calculates a distribution of points based on the arclength and
      !! on the double-exponential rule
      subroutine get_dexp_m1_ditribution( nbr, nb, aks, zcv, xbr, ybr, xb, yb) ! Output: last two
         implicit none
         integer, intent(in)  :: nbr        !< Number of partitions over the reference curve (xbr,ybr)
         integer, intent(in)  :: nb         !< Number of partitions for the new distribution of points
         real(8), intent(in)  :: aks        !< Exponent for the calculation of acv
         real(8), intent(in)  :: zcv        !< Central value of z
         real(8), intent(in)  :: xbr(0:nbr) !< Given distribution of points
         real(8), intent(in)  :: ybr(0:nbr) !< Given distribution of points
         real(8), intent(out) :: xb(0:nb)   !< New distribution of points
         real(8), intent(out) :: yb(0:nb)   !< New distribution of points

         ! Inner variables

         real(8) :: sbr(0:nbr) ! Length of each partition over the reference curve
         real(8) :: sb(0:nb)   ! Length of each partition over the new curve
         real(8) :: L          ! Total length
         real(8) :: acv        ! Length of the volume closer to the central point


         ! Calculates the length of the curve
         call get_sbr(nbr, xbr, ybr, sbr) ! Output: last one


         ! Total length
         L = sbr(nbr)


         ! Calculating the double-exponential distribution of arclength

         acv = L / ( nb ** aks )

         if ( acv > zcv ) then

            write(*,*) "Error. get_dexp_m1_ditribution: there is no resolution"&
               , " for this grid. Stopping."

            stop

         end if

         call get_dexp_m1(nb, 0.d0, zcv, L, acv, sb)


         ! Calculates the new distribution of points
         call get_new_distribution(nbr, nb, sbr, xbr, ybr, sb, xb, yb) ! Output: last two


      end subroutine get_dexp_m1_ditribution

      !> \brief Calculates a distribution of points based on the double-exponential
      !! rule. The coefficients of this rule are set by defining the position of
      !! an internal point (zcv) and the width (acv) of the volumes beside this
      !! point. The points are concentrated around zcv.
      subroutine get_dexp_m1(np, ziv, zcv, zfv, acv, z)
         implicit none
         integer, intent(in)  :: np  !< Number of partitions
         real(8), intent(in)  :: ziv !< Initial value of z
         real(8), intent(in)  :: zcv !< Central value of z (exponential's matching point)
         real(8), intent(in)  :: zfv !< Final value of z
         real(8), intent(in)  :: acv !< Length of the volumes beside zcv
         real(8), intent(out) :: z(0:np) !< Distribution of points

         ! Inner variables

         integer :: i
         integer :: npc
         real(8) :: t
         real(8) :: dt
         real(8) :: tc
         real(8) :: r
         real(8) :: A
         real(8) :: ac
         real(8) :: zc


         ! Initializing variables

         zc = zcv / ( zfv - ziv )

         ac = acv / ( zfv - ziv )

         dt = 1.d0 / dble(np)

         r = ac / ( 1.d0 - zc )



         ! Cheking conditions

         if ( ( ac / zc > 1.d0 ) .or. ( r > 1.d0 ) ) then

            write(*,*) "Error. get_dexp_m1: ac > zc or ac > (1-zc). Stopping."

            stop

         end if

         if ( ac > dt ) then

            write(*,*) "Error. get_dexp_m1: ac > dt. Stopping."

            stop

         end if


         ! Searching the value of A

         call adexp_m1(dt, r, zc, A)


         ! Calculating the fraction of volumes for each exponential

         tc = 1.d0 + log( 1.d0 + ( exp(-A*dt) - 1.d0 ) / r ) / A

         npc = nint(tc*np)



         ! Calculating the normallized distribution of points

         ! First exponential

         do i = 0, npc

            t = tc * dble(i) / dble(npc)

            z(i) = zc * ( exp(A*t) - 1.d0 ) / ( exp(A*tc) - 1.d0 )

         end do

         ! Second exponential

         do i = npc, np

            t = ( 1.d0-tc ) * dble(i-npc) / dble(np-npc) + tc

            z(i) = zc+(1.d0-zc)*(exp(-A*(t-tc))-1.d0)/(exp(-A*(1.d0-tc))-1.d0)

         end do

         z(0)  = 0.d0

         z(np) = 1.d0

         ! Calculating the extended distribution of points

         z = ( zfv - ziv ) * z + ziv

      end subroutine get_dexp_m1

      !> \brief Defines the transcendental function for finding A
      real(8) function fdexp_m1(dt, r, zc, A)
         implicit none
         real(8), intent(in) :: dt
         real(8), intent(in) :: r
         real(8), intent(in) :: zc
         real(8), intent(in) :: A

         ! Inner varialbles

         real(8) :: raux1, raux2, raux3

         if ( abs(A) > 1.d-6 ) then

            raux1 = exp(A) * ( 1.d0 + ( exp(-A*dt) - 1.d0 ) / r )
            raux2 = (2.d0*zc-1.d0) / zc
            raux3 = (1.d0-zc) / zc

            fdexp_m1 = -exp(-A) * raux1 * raux1 +  raux1 * raux2 + raux3

         else

            fdexp_m1 = (r*zc-r+dt)*A/(r*zc) &
               + ((r**2-2*dt**2)*zc-r**2+(2*dt-dt**2)*r)*A**2/(r**2*zc)/2.0

         end if

      end function fdexp_m1


      !> \brief Finds the value of A using the bissection method
      subroutine adexp_m1(dt, r, zc, A)
         implicit none
         real(8), intent(in)  :: dt
         real(8), intent(in)  :: r
         real(8), intent(in)  :: zc
         real(8), intent(out) :: A

         ! Parameters

         integer, parameter :: itmax = 1000

         ! Inner variables

         integer :: k
         real(8) :: Ai, Af, Am, Fi, Fm, Ff



         ! Initializing the search interval

         Ai = -1.d-18

         Af = -100.d0



         ! Checking for roots

         Fi = fdexp_m1(dt, r, zc, Ai)

         Ff = fdexp_m1(dt, r, zc, Af)

         if ( Fi * Ff > 0.d0 ) then

            write(*,*) "Error. adexp_m1: there is none or more than one roots. Stopping."

            stop

         end if


         ! Looking for a solution

         do k = 1, itmax

            Am = 0.5d0 * ( Ai + Af )

            Fi = fdexp_m1(dt, r, zc, Ai)

            Fm = fdexp_m1(dt, r, zc, Am)

            if ( Fi * Fm <= 0.d0 ) then

               Af = Am

            else

               Ai = Am

            end if

            if ( abs(Af-Ai) < 1.d-14 ) exit

         end do

         if ( k >= itmax ) then

            write(*,*) "Error. adexp_m1: number of iteractions exceeded. Stopping."

            stop

         end if

         A = 0.5d0 * ( Ai + Af )

      end subroutine adexp_m1


   end subroutine get_grid_boundary_g17



   !> \brief Calculates the coordinates (x,y) of the north and south boundaries.
   !! The north boundary line is generated from the south boundary.
   !! The south boundary line is the Hemisphere-cone-cylinder shape of
   !! length lo and radius rb followed by a cylinder of radius rb and length lc.
   !! The distribution of points over the south boundary is separated in three
   !! parts:
   !! Part   I: geometric progression (hemisphere-cone)
   !! Part  II: geometric progression (hemisphere-cone)
   !! Part III: geometric progression (cylinder)
   !! The distribution of points over the north boundary is  based on the normal
   !! vectors from the south boundary. The directions of the normal vectors
   !! may be changed by a smoothing procedure.
   subroutine get_grid_boundary_g18(nx, ny, unt, fgeom, lr, rb, x, y) ! Output: last two
      implicit none
      integer, intent(in) :: nx  !< Number of volumes in the csi direction (real+fictitious)
      integer, intent(in) :: ny  !< Number of volumes in the eta direction (real+fictitious)
      integer, intent(in) :: unt !< Unit where the input parameters will be printed
      character(len=*), intent(in) :: fgeom !< File of the geometric parameters

      real(8), intent(out) :: lr  !< Length of the body (m)
      real(8), intent(out) :: rb  !< Base radius of the body (m)
      real(8), dimension(nx*ny), intent(out) :: x !< Coord. x at the northeast corner of the volume P (m)
      real(8), dimension(nx*ny), intent(out) :: y !< Coord. y at the northeast corner of the volume P (m)

      ! Parameters

      integer, parameter :: nr = 10000 ! Number of partitions for the representing curve (xr,yr)


      ! Inner variables
      integer :: i, j, np ! Dummy indexes

      integer :: kgb ! Kind of geometry of the body
      integer :: nx1 ! Number of points in the first part of the grid along the eta lines
      integer :: nx2 ! Number of points in the first and second parts of the grid along the eta lines
      integer :: nx3 ! Number of points in the first to third parts of the grid along the eta lines
      integer :: iiv ! Initial value of i
      integer :: ifv ! Final value of i
      integer :: nr1 ! Fraction of nr used to represent the Part I

      real(8) :: a1  ! Width of the first volume of Part I
      real(8) :: a2  ! Width of the first volume of Part II
      real(8) :: a3  ! Width of the first volume of Part III

      real(8) :: xiv ! Initial value of x
      real(8) :: xfv ! Final value of x

      real(8) :: lo  ! Ogive length (m)
      real(8) :: w   ! Distance between south and north boundaries (m)
      real(8) :: fs  ! Factor of smoothing (0<=fs) (necessary for non-smooth surfaces)
      real(8) :: lc  ! Cylinder length (m)
      real(8) :: lo1 ! Length of the first part of the ogive (m)
      real(8) :: tht ! Semi-angle of the cone
      real(8) :: aks ! Exponent of the power law for the calculation of the width of the first volume
      real(8) :: flo ! Fraction of lo
      real(8) :: fp1 ! Fraction of points in the first part of the grid
      real(8) :: fp2 ! Fraction of points in the second part of the grid


      real(8), dimension(nx) :: xbs ! Distribution of points over the south boundary
      real(8), dimension(nx) :: ybs ! Distribution of points over the south boundary

      real(8), dimension(0:nr) :: xr ! Distribution of points over the representing curve
      real(8), dimension(0:nr) :: yr ! Distribution of points over the representing curve


      ! Reading parameters

      open(13,file=fgeom)

      read(13,*) kgb  ! Kind of geometry of the body
      read(13,*) w    ! Distance between south and north boundaries (m)
      read(13,*) fs   ! Factor of smoothing (0<=fs) (necessary for non-smooth surfaces)
      read(13,*) lo   ! Length of the ogive (m)
      read(13,*) lc   ! Length of the cylinder (m)
      read(13,*) rb   ! Base radius of the body (m)
      read(13,*) tht  ! Semi-angle of the cone (degrees)
      read(13,*) aks  ! Exponent of the power law for the south boundary
      read(13,*) flo  ! Fraction of lo
      read(13,*) fp1  ! Fraction of volumes in the first part
      read(13,*) fp2  ! Fraction of volumes in the second part

      close(13)


      ! Creating the hemisphere-cone shape

      call get_hemisphere_cone_shape(nr, unt, tht, rb, lo, xr, yr)


      ! Calculating the number of points in the first part of the grid

      nx1 = int( (nx-2) * fp1 )

      ! Calculating the number of points in the first and second parts of the grid

      nx2 = int( (nx-2) * fp2 ) + nx1

      ! Calculating the number of points in the 1st, 2nd and 3rd part of the grid

      nx3 = nx-1


      ! Calculating the fraction of lo used in the 2nd part of the grid

      lo1 = flo * lo

      do i = 1, nr

         if ( lo1 < xr(i) ) exit

      end do

      nr1 = i

      lo1 = xr(nr1)


      ! Calculating the body length

      lr = lo + lc


      ! Writing parameters

      write(unt,*)
      write(unt,*) " *** Geometry parameters ***"
      write(unt,*)

      write(unt,"(I23,    A)") kgb, " = kgb: Kind of geometry of the body"
      write(unt,"(ES23.16,A)")   w, " =   w: Distance between south and north boundaries (m)"
      write(unt,"(ES23.16,A)")  fs, " =  fs: Factor of smoothing (0<=fs) (necessary for non-smooth surfaces)"
      write(unt,"(ES23.16,A)")  lo, " =  lo: Length of the ogive (m)"
      write(unt,"(ES23.16,A)")  lc, " =  lc: Length of the cylinder (m)"
      write(unt,"(ES23.16,A)")  rb, " =  rb: Base radius of the body (m)"
      write(unt,"(ES23.16,A)") tht, " = tht: Semi-angle of the cone (degrees)"
      write(unt,"(ES23.16,A)") aks, " = aks: Exponent of the power law for the south boundary"
      write(unt,"(ES23.16,A)") flo, " = flo: Fraction of lo"
      write(unt,"(ES23.16,A)") fp1, " = fp1: Fraction of volumes in the first part of the grid"
      write(unt,"(ES23.16,A)") fp2, " = fp2: Fraction of volumes in the second part of the grid"


      ! ========================================================

      ! Generating south boundary

      j = 1

      ! Part I

      ! Calculates the distribution of points based on the arclength and the
      ! exponential distribution

      iiv = 1   ! Initial value of i

      ifv = nx1 ! Final value of i

      a1 = lo1 / (ifv-iiv) ** aks ! Width of the first volume of the first part of the grid

      call get_gp_length_distribution(nr1, ifv-iiv, a1, -1.d0, xr(0:nr1) &
         , yr(0:nr1), xbs(iiv:ifv), ybs(iiv:ifv)) ! Output: last two

      do i = iiv, ifv

         np   = nx * (j-1) + i

         x(np) = xbs(i)

         y(np) = ybs(i)

      end do


      ! Part II

      iiv = nx1  ! Initial value of i

      ifv = nx2 ! Final value of i

      xiv = lo1 ! Initial value of x

      xfv = lo  ! Final value of x

      ! Calculating the width of the first volume of Part II

      i = iiv

      a2 = sqrt( (xbs(i)-xbs(i-1))**2 + (ybs(i)-ybs(i-1))**2 )

      ! Calculating the distribution curve according to a geometric progression
      ! distribution based on the curve length

      call get_gp_length_distribution(nr-nr1, ifv-iiv, a2, -1.d0, xr(nr1:nr) &
         , yr(nr1:nr), xbs(iiv:ifv), ybs(iiv:ifv)) ! Output: last two


      do i = iiv, ifv

         np   = nx * (j-1) + i

         x(np) = xbs(i)

         y(np) = ybs(i)

      end do


      ! Part III

      iiv = nx2  ! Initial value of i

      ifv = nx3  ! Final value of i

      xiv = lo ! Initial value of x

      xfv = lo+lc  ! Final value of x

      ! Calculating the width of the first volume of Part III

      i = iiv

      a3 = sqrt( (xbs(i)-xbs(i-1))**2 + (ybs(i)-ybs(i-1))**2 )

      call get_gp_distribution(iiv, ifv, xiv, xfv, a3, xbs(iiv:ifv)) ! Output: last one

      do i = iiv, ifv

         np   = nx * (j-1) + i

         ybs(i) = rb

         x(np) = xbs(i)

         y(np) = ybs(i)

      end do


      ! ========================================================

      ! Generating north boundary

      call get_north_boundary(nx, ny, w, fs, xbs, ybs, x, y) ! Output: last two

      ! ========================================================

   contains


      !> \brief Generates the north boundary of the domain based on the normal
      !! vectors from the south boundary. The directions of the normal vectors
      !! may be changed by a smoothing procedure.
      subroutine get_north_boundary(nx, ny, w, fs, xbs, ybs, x, y) ! Output: last two
         implicit none
         integer, intent(in)    :: nx  !< Number of volumes in the csi direction (real+fictitious)
         integer, intent(in)    :: ny  !< Number of volumes in the eta direction (real+fictitious)
         real(8), intent(in)    :: w   !< Width of the domain
         real(8), intent(in)    :: fs  !< Factor of smoothing (0<=fs)
         real(8), dimension(nx),    intent(in)  :: xbs !< Distribution of points over the south boundary
         real(8), dimension(nx),    intent(in)  :: ybs !< Distribution of points over the south boundary
         real(8), dimension(nx*ny), intent(out) :: x   !< Coord. x at the northeast corner of the volume P (m)
         real(8), dimension(nx*ny), intent(out) :: y   !< Coord. y at the northeast corner of the volume P (m)

         ! Inner variables

         integer :: i, j, k, np, ns
         real(8) :: raux
         real(8) :: tx
         real(8) :: ty
         real(8) :: tm
         real(8) :: vnx(nx-1)
         real(8) :: vny(nx-1)

         ! Number of iteractions for smoothing the direction vectors
         ns = int(nx * fs)

         ! Calculating the initial normal directions

         ! Direction of the west boundary
         vnx(1) = -1.d0
         vny(1) =  0.d0


         ! Direction of the east boundary
         vnx(nx-1) = 0.d0
         vny(nx-1) = 1.d0

         ! Directions normal to the body surface
         do i = 2, nx-2

            tx = ( xbs(i+1) - xbs(i-1) ) / 2.d0
            ty = ( ybs(i+1) - ybs(i-1) ) / 2.d0

            tm = sqrt( tx * tx + ty * ty )

            vnx(i) = - ty / tm
            vny(i) =   tx / tm

         end do

         ! Calculating the final directions (smoothing)

         do k = 1, ns

            do i = 2, nx-2

               vnx(i) = ( vnx(i+1) + vnx(i-1) ) / 2.d0

               vny(i) = ( vny(i+1) + vny(i-1) ) / 2.d0

               raux = dsqrt( vnx(i) * vnx(i) + vny(i) * vny(i) )

               vnx(i) = vnx(i) / raux

               vny(i) = vny(i) / raux

            end do

         end do

         ! Generating the north boundary

         j = ny-1

         do i = 1, nx-1

            np   = nx * (j-1) + i

            x(np) = xbs(i) + w * vnx(i)

            y(np) = ybs(i) + w * vny(i)

         end do

      end subroutine get_north_boundary


      !> \brief Calculates the node distribution based on a geometric progression distribution
      subroutine get_gp_distribution(iiv, ifv, xiv, xfv, a1, xb) ! Output: last one
         implicit none
         integer, intent(in)  :: iiv !< Initial value of i
         integer, intent(in)  :: ifv !< Final value of i
         real(8), intent(in)  :: xiv !< Initial value of x
         real(8), intent(in)  :: xfv !< Final value of x
         real(8), intent(in)  :: a1  !< Width of the first volume
         real(8), intent(out) :: xb(iiv:ifv) !< Distribution of points

         ! Inner variables

         integer :: i
         integer :: n
         real(8) :: r
         real(8) :: q


         n = ifv - iiv

         r = ( xfv - xiv ) / a1

         ! Calculating q
         call get_gp_ratio(n, r, q)

         xb(iiv) = xiv

         do i = iiv+1, ifv-1

            xb(i) = xb(i-1) + a1 * q ** ( i - iiv - 1 )

         end do

         xb(ifv) = xfv

      end subroutine get_gp_distribution


      !> \brief Calculates the base q of the geometric progression
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

           if ( 0.d0 < f(n, r, qi) * f(n, r, qm) ) then

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

      end subroutine get_gp_ratio

      !> \brief Function used in the subroutine 'get_gp_ratio'
      real(8) function f(n, r, q)
         implicit none
         integer, intent(in)  ::   n !< number of partitions
         real(8), intent(in)  ::   r !< l/a1
         real(8), intent(in)  ::   q !< q

         f = q ** n + r * ( 1.d0 - q ) - 1.d0

      end function f

      !> \brief Given a reference distribution of points (xbr,ybr)
      !! a new one (xb,yb) is created based on a geometric progression distributiuon
      !! based on the curve length.
      subroutine get_gp_length_distribution( nbr, nb, al, ar, xbr, ybr, xb, yb) ! Output: last two
         implicit none
         integer, intent(in)  :: nbr        !< Number of partitions over the reference curve (xbr,ybr)
         integer, intent(in)  :: nb         !< Number of partitions for the new distribution of points
         real(8), intent(in)  :: al         !< Length of the first partition on the left
         real(8), intent(in)  :: ar         !< Length of the first partition on the right
         real(8), intent(in)  :: xbr(0:nbr) !< Given distribution of points
         real(8), intent(in)  :: ybr(0:nbr) !< Given distribution of points
         real(8), intent(out) :: xb(0:nb)   !< New distribution of points
         real(8), intent(out) :: yb(0:nb)   !< New distribution of points

         ! Inner variables

         real(8) :: sbr(0:nbr) ! Length of each partition over the reference curve
         real(8) :: sb(0:nb)   ! Length of each partition over the new curve
         real(8) :: L          ! Total length


         ! Calculates the length of the curve
         call get_sbr(nbr, xbr, ybr, sbr) ! Output: last one


         ! Total length
         L = sbr(nbr)


         ! Selecting the kind of distribution

         if ( al < 0.d0 .and. ar < 0.d0 ) then

            ! Calculates the new points according to a uniform distribution
            call get_uniform_grid( nb, L, sb) ! Output: last one

         else if ( al > 0.d0 .and. ar < 0.d0 ) then

            ! Calculates the distribution of points according to a geometric
            ! progression. The length of the first partition is defined on the left.
            call get_gp_grid_left( nb, al, L, sb)

         else if ( al < 0.d0 .and. ar > 0.d0 ) then

            ! Calculates the distribution of points according to a geometric
            ! progression. The length of the first partition is defined on the right.
            call get_gp_grid_right( nb, ar, L, sb)

         else

            ! Calculates the distribution of points according to a geometric
            ! progression. The length of the first partition is defined on the left
            ! and on the right.
            call get_gp_grid_left_right(nb, al, ar, L, sb)

         end if


         ! Calculates the new distribution of points
         call get_new_distribution(nbr, nb, sbr, xbr, ybr, sb, xb, yb) ! Output: last two


      end subroutine get_gp_length_distribution

      !> \brief Calculates the length of the curve
      subroutine get_sbr(nbr, xbr, ybr, sbr) ! Output: last one
         implicit none
         integer, intent(in)  :: nbr        !< Number of partitions over the reference curve (xbr,ybr)
         real(8), intent(in)  :: xbr(0:nbr) !< Given distribution of points
         real(8), intent(in)  :: ybr(0:nbr) !< Given distribution of points
         real(8), intent(out) :: sbr(0:nbr) !< Length of the reference curve

         ! Inner variables

         integer :: i


         sbr(0) = 0.d0

         do i = 1, nbr

            sbr(i) = sbr(i-1) + sqrt( ( xbr(i) - xbr(i-1) )**2 + ( ybr(i) - ybr(i-1) )**2 )

         end do


      end subroutine get_sbr


      !> \brief Given the length distribution, the (xb,yb) is calculated
      subroutine get_new_distribution(nbr, nb, sbr, xbr, ybr, sb, xb, yb) ! Output: last two
         integer, intent(in)  :: nbr        !< Number of partitions over the reference curve (xbr,ybr)
         integer, intent(in)  :: nb         !< Number of partitions for the new distribution of points
         real(8), intent(in)  :: sbr(0:nbr) !< Length of the curve
         real(8), intent(in)  :: xbr(0:nbr) !< Given distribution of points
         real(8), intent(in)  :: ybr(0:nbr) !< Given distribution of points
         real(8), intent(in)  :: sb(0:nb)   !< New length distribution over the curve
         real(8), intent(out) :: xb(0:nb)   !< New distribution of points
         real(8), intent(out) :: yb(0:nb)   !< New distribution of points


         ! Inner variables

         integer :: i, j

         xb(0) = xbr(0)
         yb(0) = ybr(0)

         j = 0

         do i = 1, nb-1

            do j = j, nbr

               if( sbr(j) - sb(i) >= 0.d0 ) exit

            end do

            ! Linear interpolation

            xb(i) = xbr(j-1) + ( xbr(j)-xbr(j-1) ) * ( sb(i) - sbr(j-1) ) / ( sbr(j) - sbr(j-1) )
            yb(i) = ybr(j-1) + ( ybr(j)-ybr(j-1) ) * ( sb(i) - sbr(j-1) ) / ( sbr(j) - sbr(j-1) )

         end do

         xb(nb) = xbr(nbr)
         yb(nb) = ybr(nbr)

      end subroutine get_new_distribution



      !> \brief Calculates the uniform distribution of points
      subroutine get_uniform_grid( nb, L, sb) ! Output: last one
         implicit none
         integer, intent(in)  :: nb       !< Number of partitions over the new distribution of points
         real(8), intent(in)  :: L        !< Length of the curve
         real(8), intent(out) :: sb(0:nb) !< Length distribution over the curve

         ! Inner variables

         integer :: i

         real(8) :: ds

         ds = L / dble(nb)


         sb(0) = 0.d0

         do i = 1, nb

            sb(i) = sb(i-1) + ds

         end do

      end subroutine get_uniform_grid



      !> \brief Calculates the distribution of points according
      !! to a geometric progression. The length of the first partition
      !! is defined on the left.
      subroutine get_gp_grid_left(nb, al, L, sb)
         implicit none
         integer, intent(in)  :: nb       !< Number of partitions over the new distribution of points
         real(8), intent(in)  :: al       !< Length of the first patition on the left
         real(8), intent(in)  :: L        !< Length of the curve
         real(8), intent(out) :: sb(0:nb) !< Length distribution over the curve

         ! Inner variables

         integer :: i
         real(8) :: q

         ! Calculates the ratio q of the geometric progression
         call get_gp_ratio(nb, L/al, q)

         sb(0) = 0.d0

         do i = 1, nb-1

            sb(i) = sb(i-1) + al * q ** (i-1)

         end do

         sb(nb) = L

      end subroutine get_gp_grid_left


      !> \brief Calculates the distribution of points according
      !! to a geometric progression. The length of the first partition
      !! is defined on the right.
      subroutine get_gp_grid_right(nb, ar, L, sb)
         implicit none
         integer, intent(in)  :: nb       !< Number of partitions over the new distribution of points
         real(8), intent(in)  :: ar       !< Length of the first patition on the right
         real(8), intent(in)  :: L        !< Length of the curve
         real(8), intent(out) :: sb(0:nb) !< Length distribution over the curve

         ! Inner variables

         integer :: i
         real(8) :: q

         ! Calculates the ratio q of the geometric progression
         call get_gp_ratio(nb, L/ar, q)

         sb(0) = 0.d0

         do i = 1, nb-1

            sb(i) = sb(i-1) + ar * q ** (nb-i)

         end do

         sb(nb) = L

      end subroutine get_gp_grid_right



      !> \brief Calculates the distribution of points according
      !! to a geometric progression. The length of the first partition
      !! is defined on the left and on the right.
      subroutine get_gp_grid_left_right(nb, al, ar, L, sb)
         implicit none
         integer, intent(in)  :: nb       !< Number of partitions over the new distribution of points
         real(8), intent(in)  :: al       !< Length of the first patition on the left
         real(8), intent(in)  :: ar       !< Length of the first patition on the right
         real(8), intent(in)  :: L        !< Length of the curve
         real(8), intent(out) :: sb(0:nb) !< Length distribution over the curve

         ! Inner parameters

         integer, parameter :: nitmax = 10000

         ! Inner variables

         integer :: i
         integer :: nl
         integer :: nr
         real(8) :: Lm
         real(8) :: Li
         real(8) :: Lf
         real(8) :: acl
         real(8) :: acr


         nl = nb / 2

         nr = nb - nl

         Li = 0.d0

         Lf = L

         do i = 1, nitmax

            Lm = 0.5d0 * ( Li + Lf )

            call get_gp_grid_left (nl, al, Lm, sb( 0:nl))

            call get_gp_grid_right(nr, ar, L-Lm, sb(nl:nb))

            sb(nl:nb) = sb(nl:nb) + Lm

            acl = sb(nl) - sb(nl-1)

            acr = sb(nl+1) - sb(nl)

            if ( acl > acr ) then

               Lf = Lm

            else

               Li = Lm

            end if

            if ( abs(acl-acr) < 1.d-14 ) exit

         end do

      end subroutine get_gp_grid_left_right


      !> \brief Calculates a distribution of points based on the arclength and
      !! on the double-exponential rule
      subroutine get_dexp_m1_ditribution( nbr, nb, aks, zcv, xbr, ybr, xb, yb) ! Output: last two
         implicit none
         integer, intent(in)  :: nbr        !< Number of partitions over the reference curve (xbr,ybr)
         integer, intent(in)  :: nb         !< Number of partitions for the new distribution of points
         real(8), intent(in)  :: aks        !< Exponent for the calculation of acv
         real(8), intent(in)  :: zcv        !< Central value of z
         real(8), intent(in)  :: xbr(0:nbr) !< Given distribution of points
         real(8), intent(in)  :: ybr(0:nbr) !< Given distribution of points
         real(8), intent(out) :: xb(0:nb)   !< New distribution of points
         real(8), intent(out) :: yb(0:nb)   !< New distribution of points

         ! Inner variables

         real(8) :: sbr(0:nbr) ! Length of each partition over the reference curve
         real(8) :: sb(0:nb)   ! Length of each partition over the new curve
         real(8) :: L          ! Total length
         real(8) :: acv        ! Length of the volume closer to the central point


         ! Calculates the length of the curve
         call get_sbr(nbr, xbr, ybr, sbr) ! Output: last one


         ! Total length
         L = sbr(nbr)


         ! Calculating the double-exponential distribution of arclength

         acv = L / ( nb ** aks )

         if ( acv > zcv ) then

            write(*,*) "Error. get_dexp_m1_ditribution: there is no resolution"&
               , " for this grid. Stopping."

            stop

         end if

         call get_dexp_m1(nb, 0.d0, zcv, L, acv, sb)


         ! Calculates the new distribution of points
         call get_new_distribution(nbr, nb, sbr, xbr, ybr, sb, xb, yb) ! Output: last two


      end subroutine get_dexp_m1_ditribution

      !> \brief Calculates a distribution of points based on the double-exponential
      !! rule. The coefficients of this rule are set by defining the position of
      !! an internal point (zcv) and the width (acv) of the volumes beside this
      !! point. The points are concentrated around zcv.
      subroutine get_dexp_m1(np, ziv, zcv, zfv, acv, z)
         implicit none
         integer, intent(in)  :: np  !< Number of partitions
         real(8), intent(in)  :: ziv !< Initial value of z
         real(8), intent(in)  :: zcv !< Central value of z (exponential's matching point)
         real(8), intent(in)  :: zfv !< Final value of z
         real(8), intent(in)  :: acv !< Length of the volumes beside zcv
         real(8), intent(out) :: z(0:np) !< Distribution of points

         ! Inner variables

         integer :: i
         integer :: npc
         real(8) :: t
         real(8) :: dt
         real(8) :: tc
         real(8) :: r
         real(8) :: A
         real(8) :: ac
         real(8) :: zc


         ! Initializing variables

         zc = zcv / ( zfv - ziv )

         ac = acv / ( zfv - ziv )

         dt = 1.d0 / dble(np)

         r = ac / ( 1.d0 - zc )



         ! Cheking conditions

         if ( ( ac / zc > 1.d0 ) .or. ( r > 1.d0 ) ) then

            write(*,*) "Error. get_dexp_m1: ac > zc or ac > (1-zc). Stopping."

            stop

         end if

         if ( ac > dt ) then

            write(*,*) "Error. get_dexp_m1: ac > dt. Stopping."

            stop

         end if


         ! Searching the value of A

         call adexp_m1(dt, r, zc, A)


         ! Calculating the fraction of volumes for each exponential

         tc = 1.d0 + log( 1.d0 + ( exp(-A*dt) - 1.d0 ) / r ) / A

         npc = nint(tc*np)



         ! Calculating the normallized distribution of points

         ! First exponential

         do i = 0, npc

            t = tc * dble(i) / dble(npc)

            z(i) = zc * ( exp(A*t) - 1.d0 ) / ( exp(A*tc) - 1.d0 )

         end do

         ! Second exponential

         do i = npc, np

            t = ( 1.d0-tc ) * dble(i-npc) / dble(np-npc) + tc

            z(i) = zc+(1.d0-zc)*(exp(-A*(t-tc))-1.d0)/(exp(-A*(1.d0-tc))-1.d0)

         end do

         z(0)  = 0.d0

         z(np) = 1.d0

         ! Calculating the extended distribution of points

         z = ( zfv - ziv ) * z + ziv

      end subroutine get_dexp_m1

      !> \brief Defines the transcendental function for finding A
      real(8) function fdexp_m1(dt, r, zc, A)
         implicit none
         real(8), intent(in) :: dt
         real(8), intent(in) :: r
         real(8), intent(in) :: zc
         real(8), intent(in) :: A

         ! Inner varialbles

         real(8) :: raux1, raux2, raux3

         if ( abs(A) > 1.d-6 ) then

            raux1 = exp(A) * ( 1.d0 + ( exp(-A*dt) - 1.d0 ) / r )
            raux2 = (2.d0*zc-1.d0) / zc
            raux3 = (1.d0-zc) / zc

            fdexp_m1 = -exp(-A) * raux1 * raux1 +  raux1 * raux2 + raux3

         else

            fdexp_m1 = (r*zc-r+dt)*A/(r*zc) &
               + ((r**2-2*dt**2)*zc-r**2+(2*dt-dt**2)*r)*A**2/(r**2*zc)/2.0

         end if

      end function fdexp_m1


      !> \brief Finds the value of A using the bissection method
      subroutine adexp_m1(dt, r, zc, A)
         implicit none
         real(8), intent(in)  :: dt
         real(8), intent(in)  :: r
         real(8), intent(in)  :: zc
         real(8), intent(out) :: A

         ! Parameters

         integer, parameter :: itmax = 1000

         ! Inner variables

         integer :: k
         real(8) :: Ai, Af, Am, Fi, Fm, Ff



         ! Initializing the search interval

         Ai = -1.d-18

         Af = -100.d0



         ! Checking for roots

         Fi = fdexp_m1(dt, r, zc, Ai)

         Ff = fdexp_m1(dt, r, zc, Af)

         if ( Fi * Ff > 0.d0 ) then

            write(*,*) "Error. adexp_m1: there is none or more than one roots. Stopping."

            stop

         end if


         ! Looking for a solution

         do k = 1, itmax

            Am = 0.5d0 * ( Ai + Af )

            Fi = fdexp_m1(dt, r, zc, Ai)

            Fm = fdexp_m1(dt, r, zc, Am)

            if ( Fi * Fm <= 0.d0 ) then

               Af = Am

            else

               Ai = Am

            end if

            if ( abs(Af-Ai) < 1.d-14 ) exit

         end do

         if ( k >= itmax ) then

            write(*,*) "Error. adexp_m1: number of iteractions exceeded. Stopping."

            stop

         end if

         A = 0.5d0 * ( Ai + Af )

      end subroutine adexp_m1

      !> \brief Calculates the hemisphere-cone shape
      subroutine get_hemisphere_cone_shape(nr, unt, tht, rb, lo, xr, yr)
         implicit none
         integer, intent(in)  :: nr  !< Number of partitions
         integer, intent(in)  :: unt !< Unit where the results will be printed
         real(8), intent(in)  :: tht !< Semi-angle of the cone
         real(8), intent(in)  :: rb  !< Base radius
         real(8), intent(in)  :: lo  !< Ogive length
         real(8), intent(out) :: xr(0:nr) !< x coordinate
         real(8), intent(out) :: yr(0:nr) !< y coordinate

         ! Parameters

         real(8), parameter :: pi = acos(-1.d0)

         ! Inner variables
         real(8) :: ac ! Tangent of the cone semi-angle
         real(8) :: xm ! x of the matching point between the hemisphere and the cone (m)
         real(8) :: rh ! radius of the hemisphere (m)


         ! Calculating coefficients necessary to the body calculation

         ! Tangent of the cone semi-angle
         ac = tan( tht * pi / 180.d0 )

         ! x of the matching point between the hemisphere and the cone (m)
         xm = ( rb - ac * lo ) / sqrt( 1.d0 + ac * ac )

         if ( xm < 0.d0 ) then

            write(*,*) "ERROR: Semi-angle of the cone is greater than the allowed" &
               // " value (degree): ", atan(rb/lo)*180.d0/pi

            stop

         end if

         ! radius of the hemisphere (m)
         rh = xm + ac * ( ac * (xm-lo) + rb )


         write(unt,"(ES23.16,A)") xm, " =  xm: x of the matching point between the hemisphere and the cone (m)"
         write(unt,"(ES23.16,A)") rh, " =  rh: radius of the hemisphere (m)"


         do i = 0, nr

            xr(i) = lo * ( dble(i) / dble(nr) ) ** 2.d0

            yr(i) = fhc(xr(i), xm, ac, lo, rb, rh)

         end do


      end subroutine get_hemisphere_cone_shape


      !> \brief Hemisphere-cone function
      real(8) function fhc(x, xm, ac, lo, ro, rh)
         implicit none
         real(8), intent(in) :: x  !< x coordinate (m)
         real(8), intent(in) :: xm !< x coordinate of intersection between the hemisphere and the cone(m)
         real(8), intent(in) :: ac !< Tangent of tht
         real(8), intent(in) :: lo !< Length of the ogive (m)
         real(8), intent(in) :: ro !< Radius of the ogive (m)
         real(8), intent(in) :: rh !< Radius of the hemisphere (m)

         if ( x < xm ) then

            fhc = sqrt( rh * rh - ( x - rh ) ** 2 )

         else

            fhc = ac * ( x -lo ) + ro

         end if

      end function


   end subroutine get_grid_boundary_g18


   !> \brief Calculates the coordinates (x,y) of the north and south boundaries.
   !! The north boundary line is generated from the south boundary.
   !! The south boundary line is an bluff nose followed by a shifted power law ogive of
   !! length lo and radius rb perturbed with a spline and followed by a cylinder of
   !! radius rb and length lc.
   !! The distribution of points over the south boundary is separated in three
   !! parts:
   !! Part   I: geometric progression (bluff nose+shifted power law)
   !! Part  II: geometric progression (shifted power law)
   !! Part III: geometric progression (cylinder)
   !! The distribution of points over the north boundary is  based on the normal
   !! vectors from the south boundary. The directions of the normal vectors
   !! may be changed by a smoothing procedure.
   !! Adapted to DEPP optimizer.
   subroutine get_grid_boundary_g19(nx, ny, unt, fgeom, lr, rb, x, y) ! Output: last two
      implicit none
      integer, intent(in) :: nx  !< Number of volumes in the csi direction (real+fictitious)
      integer, intent(in) :: ny  !< Number of volumes in the eta direction (real+fictitious)
      integer, intent(in) :: unt !< Unit where the input parameters will be printed
      character(len=*), intent(in) :: fgeom !< File of the geometric parameters

      real(8), intent(out) :: lr  !< Length of the body (m)
      real(8), intent(out) :: rb  !< Base radius of the body (m)
      real(8), dimension(nx*ny), intent(out) :: x !< Coord. x at the northeast corner of the volume P (m)
      real(8), dimension(nx*ny), intent(out) :: y !< Coord. y at the northeast corner of the volume P (m)

      ! Parameters

      integer, parameter :: nr = 10000 ! Number of partitions for the representing curve (xr,yr)


      ! Inner variables
      integer :: i, j, np ! Dummy indexes
      integer :: es  ! Exit status ( 0=success, 1=failure )
      integer :: mc  ! Monocity check ( 0=no, 1=yes )
      integer :: nsi ! Number of points for the spline interpolation

      integer :: kgb ! Kind of geometry of the body
      integer :: nx1 ! Number of points in the first part of the grid along the eta line
      integer :: nx2 ! Number of points in the first and second parts of the grid along the eta lines
      integer :: nx3 ! Number of points in the first to third parts of the grid along the eta lines
      integer :: iiv ! Initial value of i
      integer :: ifv ! Final value of i
      integer :: nr1 ! Fraction of nr used to represent the Part I

      real(8) :: a2  ! Width of the first volume of Part II
      real(8) :: a3  ! Width of the first volume of Part III

      real(8) :: xiv ! Initial value of x
      real(8) :: xfv ! Final value of x

      real(8) :: lbd ! Exponent of the power law ogive

      real(8) :: lo  ! Ogive length (m)
      real(8) :: w   ! Distance between south and north boundaries (m)
      real(8) :: fs  ! Factor of smoothing (0<=fs) (necessary for non-smooth surfaces)
      real(8) :: lc  ! Cylinder length (m)
      real(8) :: lo1 ! Length of the first part of the ogive (m)
      real(8) :: rf  ! Radius of the nose tip (m)
      real(8) :: h   ! Ratio rf/rb
      real(8) :: aks ! Exponent of the power law for the calculation of the width of the first volume
      real(8) :: flo ! Fraction of lo
      real(8) :: fp1 ! Fraction of points in the first part of the grid
      real(8) :: fp2 ! Fraction of points in the second part of the grid


      real(8), dimension(nx) :: xbs ! Distribution of points over the south boundary
      real(8), dimension(nx) :: ybs ! Distribution of points over the south boundary

      real(8), dimension(0:nr) :: xr ! Distribution of points over the representing curve
      real(8), dimension(0:nr) :: yr ! Distribution of points over the representing curve

      real(8), allocatable :: xsi(:) ! Data for spline interpolation
      real(8), allocatable :: ysi(:) ! Data for spline interpolation

      ! DEPP variables
      integer            :: nu     ! number of unknowns
      integer            :: ind    ! number of the individual
      integer,       allocatable :: xopt(:)   ! Optimization checker
      character(10), allocatable :: xname(:)  ! Name of parameters
      real(8),       allocatable :: xvalue(:) ! parameters
      character(200)             :: sname     ! simulation name


      ! Reading parameters

      open(13,file=fgeom)

      read(13,*) kgb  ! Kind of geometry of the body
      read(13,*) w    ! Distance between south and north boundaries (m)
      read(13,*) fs   ! Factor of smoothing (0<=fs) (necessary for non-smooth surfaces)
      read(13,*) lo   ! Length of the ogive (m)
      read(13,*) lc   ! Length of the cylinder (m)
      read(13,*) rb   ! Base radius of the body (m)
      read(13,*) aks  ! Exponent of the power law for the south boundary
      read(13,*) flo  ! Fraction of lo
      read(13,*) fp1  ! Fraction of volumes in the first part
      read(13,*) fp2  ! Fraction of volumes in the second part
      read(13,*) mc   ! Monocity check ( 0=no, 1=yes )
      read(13,*) nu   ! Number of unkowns

      if ( nu < 2 .or. mod(nu,2) /= 0 ) then

         write(*,*) "nu must be at least 2 and even. Stopping."

         stop

      end if

      nsi = (nu-2) / 2

      allocate ( xopt(nu), xname(nu), xvalue(nu) )

      allocate( xsi(nsi), ysi(nsi) )

      do i = 1, nu

         read(13,*) xname(i) ! Parameter name
         read(13,*) xopt(i)  ! Will this parameter be optimized? ( 0 = no, 1 = yes )
         read(13,*) xvalue(i)! Parameter value

      end do

      close(13)


      ! Reads the parameters from DEPP
      call depp_get_parameters(nu, xopt, xname, xvalue, ind, sname)

      h = xvalue(1) ! Radius ratio

      lbd = xvalue(2) ! Exponent of the power law

      do i = 3, nu, 2

         j = (i-1) / 2

         xsi(j) = xvalue(i)
         ysi(j) = xvalue(i+1)

      end do

      ! Gets the Shifted Power Law With Spline Perturbation geometry
      call get_splwsp(nr, nsi, 1.5d0, rb, lo, h, lbd, xsi, ysi, xr, yr, es) ! Output: last three

      ! Checking function monotonocity
      if ( mc == 1 ) then

         ! Non-monotonic function. Requesting new parameters.
         if ( es == 1 ) then

            call depp_save_fitness(-huge(1.d0), 2, "-cdfp")

            write(*,*) "Non-monotonic geometry. Stopping."

            stop

         end if

      end if

      ! Writing parameters

      write(unt,*)
      write(unt,*) " *** Geometry parameters ***"
      write(unt,*)

      write(unt,"(I23,    A)") kgb, " = kgb: Kind of geometry of the body"
      write(unt,"(ES23.16,A)")   w,  " =   w: Distance between south and north boundaries (m)"
      write(unt,"(ES23.16,A)")  fs, " =  fs: Factor of smoothing (0<=fs) (necessary for non-smooth surfaces)"
      write(unt,"(ES23.16,A)")  lo, " =  lo: Length of the ogive (m)"
      write(unt,"(ES23.16,A)")  lc, " =  lc: Length of the cylinder (m)"
      write(unt,"(ES23.16,A)")  rb, " =  rb: Base radius of the body (m)"
      write(unt,"(ES23.16,A)")   h, " =   h: Ratio rf/rb (frontal to base radius ratio)"
      write(unt,"(ES23.16,A)") lbd, " = lbd: Exponent of the power law ogive"
      write(unt,"(ES23.16,A)") aks, " = aks: Exponent of the power law for the south boundary"
      write(unt,"(ES23.16,A)") flo, " = flo: Fraction of lo"
      write(unt,"(ES23.16,A)") fp1, " = fp1: Fraction of volumes in the first part of the grid"
      write(unt,"(ES23.16,A)") fp2, " = fp2: Fraction of volumes in the second part of the grid"

      ! Calculating the number of points in the first part of the grid

      nx1 = int( (nx-2) * fp1 )

      ! Calculating the number of points in the first and second parts of the grid

      nx2 = int( (nx-2) * fp2 ) + nx1

      ! Calculating the number of points in the 1st, 2nd and 3rd part of the grid

      nx3 = nx-1

      ! Calculating frontal radius

      rf = h * rb

      ! Calculating the fraction of lo used in the 2nd part of the grid

      lo1 = flo * lo

      do i = 1, nr

         if ( lo1 < xr(i) ) exit

      end do

      nr1 = i

      lo1 = xr(nr1)


      ! Calculating the body length

      lr = lo + lc


      ! ========================================================

      ! Generating south boundary

      j = 1

      ! Part I

      ! Calculates the distribution of points based on the arclength and the
      ! double exponential distribution

      iiv = 1   ! Initial value of i

      ifv = nx1 ! Final value of i

      call get_dexp_m1_ditribution( nr1, ifv-iiv, aks, rf, xr(0:nr1), yr(0:nr1)&
         , xbs(iiv:ifv), ybs(iiv:ifv)) ! Output: last two

      do i = iiv, ifv

         np   = nx * (j-1) + i

         x(np) = xbs(i)

         y(np) = ybs(i)

      end do


      ! Part II

      iiv = nx1  ! Initial value of i

      ifv = nx2 ! Final value of i

      xiv = lo1 ! Initial value of x

      xfv = lo  ! Final value of x

      ! Calculating the width of the first volume of Part II

      i = iiv

      a2 = sqrt( (xbs(i)-xbs(i-1))**2 + (ybs(i)-ybs(i-1))**2 )

      ! Calculating the distribution curve according to a geometric progression
      ! distribution based on the curve length

      call get_gp_length_distribution(nr-nr1, ifv-iiv, a2, -1.d0, xr(nr1:nr) &
         , yr(nr1:nr), xbs(iiv:ifv), ybs(iiv:ifv)) ! Output: last two

      do i = iiv, ifv

         np   = nx * (j-1) + i

         x(np) = xbs(i)

         y(np) = ybs(i)

      end do


      ! Part III

      iiv = nx2  ! Initial value of i

      ifv = nx3  ! Final value of i

      xiv = lo ! Initial value of x

      xfv = lo+lc  ! Final value of x

      ! Calculating the width of the first volume of Part IV

      i = iiv

      a3 = sqrt( (xbs(i)-xbs(i-1))**2 + (ybs(i)-ybs(i-1))**2 )

      call get_gp_distribution(iiv, ifv, xiv, xfv, a3, xbs(iiv:ifv)) ! Output: last one

      do i = iiv, ifv

         np   = nx * (j-1) + i

         ybs(i) = rb

         x(np) = xbs(i)

         y(np) = ybs(i)

      end do


      ! ========================================================

      ! Generating north boundary

      call get_north_boundary(nx, ny, w, fs, xbs, ybs, x, y) ! Output: last two

      ! ========================================================

   contains


      !> \brief Generates the north boundary of the domain based on the normal
      !! vectors from the south boundary. The directions of the normal vectors
      !! may be changed by a smoothing procedure.
      subroutine get_north_boundary(nx, ny, w, fs, xbs, ybs, x, y) ! Output: last two
         implicit none
         integer, intent(in)    :: nx  !< Number of volumes in the csi direction (real+fictitious)
         integer, intent(in)    :: ny  !< Number of volumes in the eta direction (real+fictitious)
         real(8), intent(in)    :: w   !< Width of the domain
         real(8), intent(in)    :: fs  !< Factor of smoothing (0<=fs)
         real(8), dimension(nx),    intent(in)  :: xbs !< Distribution of points over the south boundary
         real(8), dimension(nx),    intent(in)  :: ybs !< Distribution of points over the south boundary
         real(8), dimension(nx*ny), intent(out) :: x   !< Coord. x at the northeast corner of the volume P (m)
         real(8), dimension(nx*ny), intent(out) :: y   !< Coord. y at the northeast corner of the volume P (m)

         ! Inner variables

         integer :: i, j, k, np, ns
         real(8) :: raux
         real(8) :: tx
         real(8) :: ty
         real(8) :: tm
         real(8) :: vnx(nx-1)
         real(8) :: vny(nx-1)

         ! Number of iteractions for smoothing the direction vectors
         ns = int(nx * fs)

         ! Calculating the initial normal directions

         ! Direction of the west boundary
         vnx(1) = -1.d0
         vny(1) =  0.d0


         ! Direction of the east boundary
         vnx(nx-1) = 0.d0
         vny(nx-1) = 1.d0

         ! Directions normal to the body surface
         do i = 2, nx-2

            tx = ( xbs(i+1) - xbs(i-1) ) / 2.d0
            ty = ( ybs(i+1) - ybs(i-1) ) / 2.d0

            tm = sqrt( tx * tx + ty * ty )

            vnx(i) = - ty / tm
            vny(i) =   tx / tm

         end do

         ! Calculating the final directions (smoothing)

         do k = 1, ns

            do i = 2, nx-2

               vnx(i) = ( vnx(i+1) + vnx(i-1) ) / 2.d0

               vny(i) = ( vny(i+1) + vny(i-1) ) / 2.d0

               raux = dsqrt( vnx(i) * vnx(i) + vny(i) * vny(i) )

               vnx(i) = vnx(i) / raux

               vny(i) = vny(i) / raux

            end do

         end do

         ! Generating the north boundary

         j = ny-1

         do i = 1, nx-1

            np   = nx * (j-1) + i

            x(np) = xbs(i) + w * vnx(i)

            y(np) = ybs(i) + w * vny(i)

         end do

      end subroutine get_north_boundary


      !> \brief Calculates the node distribution based on a geometric progression distribution
      subroutine get_gp_distribution(iiv, ifv, xiv, xfv, a1, xb) ! Output: last one
         implicit none
         integer, intent(in)  :: iiv !< Initial value of i
         integer, intent(in)  :: ifv !< Final value of i
         real(8), intent(in)  :: xiv !< Initial value of x
         real(8), intent(in)  :: xfv !< Final value of x
         real(8), intent(in)  :: a1  !< Width of the first volume
         real(8), intent(out) :: xb(iiv:ifv) !< Distribution of points

         ! Inner variables

         integer :: i
         integer :: n
         real(8) :: r
         real(8) :: q


         n = ifv - iiv

         r = ( xfv - xiv ) / a1

         ! Calculating q
         call get_gp_ratio(n, r, q)

         xb(iiv) = xiv

         do i = iiv+1, ifv-1

            xb(i) = xb(i-1) + a1 * q ** ( i - iiv - 1 )

         end do

         xb(ifv) = xfv

      end subroutine get_gp_distribution


      !> \brief Calculates the base q of the geometric progression
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

           if ( 0.d0 < f(n, r, qi) * f(n, r, qm) ) then

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

      end subroutine get_gp_ratio

      !> \brief Function used in the subroutine 'get_gp_ratio'
      real(8) function f(n, r, q)
         implicit none
         integer, intent(in)  ::   n !< number of partitions
         real(8), intent(in)  ::   r !< l/a1
         real(8), intent(in)  ::   q !< q

         f = q ** n + r * ( 1.d0 - q ) - 1.d0

      end function f

      !> \brief Given a reference distribution of points (xbr,ybr)
      !! a new one (xb,yb) is created based on a geometric progression distributiuon
      !! based on the curve length.
      subroutine get_gp_length_distribution( nbr, nb, al, ar, xbr, ybr, xb, yb) ! Output: last two
         implicit none
         integer, intent(in)  :: nbr        !< Number of partitions over the reference curve (xbr,ybr)
         integer, intent(in)  :: nb         !< Number of partitions for the new distribution of points
         real(8), intent(in)  :: al         !< Length of the first partition on the left
         real(8), intent(in)  :: ar         !< Length of the first partition on the right
         real(8), intent(in)  :: xbr(0:nbr) !< Given distribution of points
         real(8), intent(in)  :: ybr(0:nbr) !< Given distribution of points
         real(8), intent(out) :: xb(0:nb)   !< New distribution of points
         real(8), intent(out) :: yb(0:nb)   !< New distribution of points

         ! Inner variables

         real(8) :: sbr(0:nbr) ! Length of each partition over the reference curve
         real(8) :: sb(0:nb)   ! Length of each partition over the new curve
         real(8) :: L          ! Total length


         ! Calculates the length of the curve
         call get_sbr(nbr, xbr, ybr, sbr) ! Output: last one


         ! Total length
         L = sbr(nbr)


         ! Selecting the kind of distribution

         if ( al < 0.d0 .and. ar < 0.d0 ) then

            ! Calculates the new points according to a uniform distribution
            call get_uniform_grid( nb, L, sb) ! Output: last one

         else if ( al > 0.d0 .and. ar < 0.d0 ) then

            ! Calculates the distribution of points according to a geometric
            ! progression. The length of the first partition is defined on the left.
            call get_gp_grid_left( nb, al, L, sb)

         else if ( al < 0.d0 .and. ar > 0.d0 ) then

            ! Calculates the distribution of points according to a geometric
            ! progression. The length of the first partition is defined on the right.
            call get_gp_grid_right( nb, ar, L, sb)

         else

            ! Calculates the distribution of points according to a geometric
            ! progression. The length of the first partition is defined on the left
            ! and on the right.
            call get_gp_grid_left_right(nb, al, ar, L, sb)

         end if


         ! Calculates the new distribution of points
         call get_new_distribution(nbr, nb, sbr, xbr, ybr, sb, xb, yb) ! Output: last two


      end subroutine get_gp_length_distribution

      !> \brief Calculates the length of the curve
      subroutine get_sbr(nbr, xbr, ybr, sbr) ! Output: last one
         implicit none
         integer, intent(in)  :: nbr        !< Number of partitions over the reference curve (xbr,ybr)
         real(8), intent(in)  :: xbr(0:nbr) !< Given distribution of points
         real(8), intent(in)  :: ybr(0:nbr) !< Given distribution of points
         real(8), intent(out) :: sbr(0:nbr) !< Length of the reference curve

         ! Inner variables

         integer :: i


         sbr(0) = 0.d0

         do i = 1, nbr

            sbr(i) = sbr(i-1) + sqrt( ( xbr(i) - xbr(i-1) )**2 + ( ybr(i) - ybr(i-1) )**2 )

         end do


      end subroutine get_sbr


      !> \brief Given the length distribution, the (xb,yb) is calculated
      subroutine get_new_distribution(nbr, nb, sbr, xbr, ybr, sb, xb, yb) ! Output: last two
         integer, intent(in)  :: nbr        !< Number of partitions over the reference curve (xbr,ybr)
         integer, intent(in)  :: nb         !< Number of partitions for the new distribution of points
         real(8), intent(in)  :: sbr(0:nbr) !< Length of the curve
         real(8), intent(in)  :: xbr(0:nbr) !< Given distribution of points
         real(8), intent(in)  :: ybr(0:nbr) !< Given distribution of points
         real(8), intent(in)  :: sb(0:nb)   !< New length distribution over the curve
         real(8), intent(out) :: xb(0:nb)   !< New distribution of points
         real(8), intent(out) :: yb(0:nb)   !< New distribution of points


         ! Inner variables

         integer :: i, j

         xb(0) = xbr(0)
         yb(0) = ybr(0)

         j = 0

         do i = 1, nb-1

            do j = j, nbr

               if( sbr(j) - sb(i) >= 0.d0 ) exit

            end do

            ! Linear interpolation

            xb(i) = xbr(j-1) + ( xbr(j)-xbr(j-1) ) * ( sb(i) - sbr(j-1) ) / ( sbr(j) - sbr(j-1) )
            yb(i) = ybr(j-1) + ( ybr(j)-ybr(j-1) ) * ( sb(i) - sbr(j-1) ) / ( sbr(j) - sbr(j-1) )

         end do

         xb(nb) = xbr(nbr)
         yb(nb) = ybr(nbr)

      end subroutine get_new_distribution



      !> \brief Calculates the uniform distribution of points
      subroutine get_uniform_grid( nb, L, sb) ! Output: last one
         implicit none
         integer, intent(in)  :: nb       !< Number of partitions over the new distribution of points
         real(8), intent(in)  :: L        !< Length of the curve
         real(8), intent(out) :: sb(0:nb) !< Length distribution over the curve

         ! Inner variables

         integer :: i

         real(8) :: ds

         ds = L / dble(nb)


         sb(0) = 0.d0

         do i = 1, nb

            sb(i) = sb(i-1) + ds

         end do

      end subroutine get_uniform_grid



      !> \brief Calculates the distribution of points according
      !! to a geometric progression. The length of the first partition
      !! is defined on the left.
      subroutine get_gp_grid_left(nb, al, L, sb)
         implicit none
         integer, intent(in)  :: nb       !< Number of partitions over the new distribution of points
         real(8), intent(in)  :: al       !< Length of the first patition on the left
         real(8), intent(in)  :: L        !< Length of the curve
         real(8), intent(out) :: sb(0:nb) !< Length distribution over the curve

         ! Inner variables

         integer :: i
         real(8) :: q

         ! Calculates the ratio q of the geometric progression
         call get_gp_ratio(nb, L/al, q)

         sb(0) = 0.d0

         do i = 1, nb-1

            sb(i) = sb(i-1) + al * q ** (i-1)

         end do

         sb(nb) = L

      end subroutine get_gp_grid_left


      !> \brief Calculates the distribution of points according
      !! to a geometric progression. The length of the first partition
      !! is defined on the right.
      subroutine get_gp_grid_right(nb, ar, L, sb)
         implicit none
         integer, intent(in)  :: nb       !< Number of partitions over the new distribution of points
         real(8), intent(in)  :: ar       !< Length of the first patition on the right
         real(8), intent(in)  :: L        !< Length of the curve
         real(8), intent(out) :: sb(0:nb) !< Length distribution over the curve

         ! Inner variables

         integer :: i
         real(8) :: q

         ! Calculates the ratio q of the geometric progression
         call get_gp_ratio(nb, L/ar, q)

         sb(0) = 0.d0

         do i = 1, nb-1

            sb(i) = sb(i-1) + ar * q ** (nb-i)

         end do

         sb(nb) = L

      end subroutine get_gp_grid_right



      !> \brief Calculates the distribution of points according
      !! to a geometric progression. The length of the first partition
      !! is defined on the left and on the right.
      subroutine get_gp_grid_left_right(nb, al, ar, L, sb)
         implicit none
         integer, intent(in)  :: nb       !< Number of partitions over the new distribution of points
         real(8), intent(in)  :: al       !< Length of the first patition on the left
         real(8), intent(in)  :: ar       !< Length of the first patition on the right
         real(8), intent(in)  :: L        !< Length of the curve
         real(8), intent(out) :: sb(0:nb) !< Length distribution over the curve

         ! Inner parameters

         integer, parameter :: nitmax = 10000

         ! Inner variables

         integer :: i
         integer :: nl
         integer :: nr
         real(8) :: Lm
         real(8) :: Li
         real(8) :: Lf
         real(8) :: acl
         real(8) :: acr


         nl = nb / 2

         nr = nb - nl

         Li = 0.d0

         Lf = L

         do i = 1, nitmax

            Lm = 0.5d0 * ( Li + Lf )

            call get_gp_grid_left (nl, al, Lm, sb( 0:nl))

            call get_gp_grid_right(nr, ar, L-Lm, sb(nl:nb))

            sb(nl:nb) = sb(nl:nb) + Lm

            acl = sb(nl) - sb(nl-1)

            acr = sb(nl+1) - sb(nl)

            if ( acl > acr ) then

               Lf = Lm

            else

               Li = Lm

            end if

            if ( abs(acl-acr) < 1.d-14 ) exit

         end do

      end subroutine get_gp_grid_left_right


      !> \brief Calculates a distribution of points based on the arclength and
      !! on the double-exponential rule
      subroutine get_dexp_m1_ditribution( nbr, nb, aks, zcv, xbr, ybr, xb, yb) ! Output: last two
         implicit none
         integer, intent(in)  :: nbr        !< Number of partitions over the reference curve (xbr,ybr)
         integer, intent(in)  :: nb         !< Number of partitions for the new distribution of points
         real(8), intent(in)  :: aks        !< Exponent for the calculation of acv
         real(8), intent(in)  :: zcv        !< Central value of z
         real(8), intent(in)  :: xbr(0:nbr) !< Given distribution of points
         real(8), intent(in)  :: ybr(0:nbr) !< Given distribution of points
         real(8), intent(out) :: xb(0:nb)   !< New distribution of points
         real(8), intent(out) :: yb(0:nb)   !< New distribution of points

         ! Inner variables

         real(8) :: sbr(0:nbr) ! Length of each partition over the reference curve
         real(8) :: sb(0:nb)   ! Length of each partition over the new curve
         real(8) :: L          ! Total length
         real(8) :: acv        ! Length of the volume closer to the central point


         ! Calculates the length of the curve
         call get_sbr(nbr, xbr, ybr, sbr) ! Output: last one


         ! Total length
         L = sbr(nbr)


         ! Calculating the double-exponential distribution of arclength

         acv = L / ( nb ** aks )

         if ( acv > zcv ) then

            write(*,*) "Error. get_dexp_m1_ditribution: there is no resolution"&
               , " for this grid. Stopping."

            stop

         end if

         call get_dexp_m1(nb, 0.d0, zcv, L, acv, sb)


         ! Calculates the new distribution of points
         call get_new_distribution(nbr, nb, sbr, xbr, ybr, sb, xb, yb) ! Output: last two


      end subroutine get_dexp_m1_ditribution

      !> \brief Calculates a distribution of points based on the double-exponential
      !! rule. The coefficients of this rule are set by defining the position of
      !! an internal point (zcv) and the width (acv) of the volumes beside this
      !! point. The points are concentrated around zcv.
      subroutine get_dexp_m1(np, ziv, zcv, zfv, acv, z)
         implicit none
         integer, intent(in)  :: np  !< Number of partitions
         real(8), intent(in)  :: ziv !< Initial value of z
         real(8), intent(in)  :: zcv !< Central value of z (exponential's matching point)
         real(8), intent(in)  :: zfv !< Final value of z
         real(8), intent(in)  :: acv !< Length of the volumes beside zcv
         real(8), intent(out) :: z(0:np) !< Distribution of points

         ! Inner variables

         integer :: i
         integer :: npc
         real(8) :: t
         real(8) :: dt
         real(8) :: tc
         real(8) :: r
         real(8) :: A
         real(8) :: ac
         real(8) :: zc


         ! Initializing variables

         zc = zcv / ( zfv - ziv )

         ac = acv / ( zfv - ziv )

         dt = 1.d0 / dble(np)

         r = ac / ( 1.d0 - zc )



         ! Cheking conditions

         if ( ( ac / zc > 1.d0 ) .or. ( r > 1.d0 ) ) then

            write(*,*) "Error. get_dexp_m1: ac > zc or ac > (1-zc). Stopping."

            stop

         end if

         if ( ac > dt ) then

            write(*,*) "Error. get_dexp_m1: ac > dt. Stopping."

            stop

         end if


         ! Searching the value of A

         call adexp_m1(dt, r, zc, A)


         ! Calculating the fraction of volumes for each exponential

         tc = 1.d0 + log( 1.d0 + ( exp(-A*dt) - 1.d0 ) / r ) / A

         npc = nint(tc*np)



         ! Calculating the normallized distribution of points

         ! First exponential

         do i = 0, npc

            t = tc * dble(i) / dble(npc)

            z(i) = zc * ( exp(A*t) - 1.d0 ) / ( exp(A*tc) - 1.d0 )

         end do

         ! Second exponential

         do i = npc, np

            t = ( 1.d0-tc ) * dble(i-npc) / dble(np-npc) + tc

            z(i) = zc+(1.d0-zc)*(exp(-A*(t-tc))-1.d0)/(exp(-A*(1.d0-tc))-1.d0)

         end do

         z(0)  = 0.d0

         z(np) = 1.d0

         ! Calculating the extended distribution of points

         z = ( zfv - ziv ) * z + ziv

      end subroutine get_dexp_m1

      !> \brief Defines the transcendental function for finding A
      real(8) function fdexp_m1(dt, r, zc, A)
         implicit none
         real(8), intent(in) :: dt
         real(8), intent(in) :: r
         real(8), intent(in) :: zc
         real(8), intent(in) :: A

         ! Inner varialbles

         real(8) :: raux1, raux2, raux3

         if ( abs(A) > 1.d-6 ) then

            raux1 = exp(A) * ( 1.d0 + ( exp(-A*dt) - 1.d0 ) / r )
            raux2 = (2.d0*zc-1.d0) / zc
            raux3 = (1.d0-zc) / zc

            fdexp_m1 = -exp(-A) * raux1 * raux1 +  raux1 * raux2 + raux3

         else

            fdexp_m1 = (r*zc-r+dt)*A/(r*zc) &
               + ((r**2-2*dt**2)*zc-r**2+(2*dt-dt**2)*r)*A**2/(r**2*zc)/2.0

         end if

      end function fdexp_m1


      !> \brief Finds the value of A using the bissection method
      subroutine adexp_m1(dt, r, zc, A)
         implicit none
         real(8), intent(in)  :: dt
         real(8), intent(in)  :: r
         real(8), intent(in)  :: zc
         real(8), intent(out) :: A

         ! Parameters

         integer, parameter :: itmax = 1000

         ! Inner variables

         integer :: k
         real(8) :: Ai, Af, Am, Fi, Fm, Ff



         ! Initializing the search interval

         Ai = -1.d-18

         Af = -100.d0



         ! Checking for roots

         Fi = fdexp_m1(dt, r, zc, Ai)

         Ff = fdexp_m1(dt, r, zc, Af)

         if ( Fi * Ff > 0.d0 ) then

            write(*,*) "Error. adexp_m1: there is none or more than one roots. Stopping."

            stop

         end if


         ! Looking for a solution

         do k = 1, itmax

            Am = 0.5d0 * ( Ai + Af )

            Fi = fdexp_m1(dt, r, zc, Ai)

            Fm = fdexp_m1(dt, r, zc, Am)

            if ( Fi * Fm <= 0.d0 ) then

               Af = Am

            else

               Ai = Am

            end if

            if ( abs(Af-Ai) < 1.d-14 ) exit

         end do

         if ( k >= itmax ) then

            write(*,*) "Error. adexp_m1: number of iteractions exceeded. Stopping."

            stop

         end if

         A = 0.5d0 * ( Ai + Af )

      end subroutine adexp_m1


   end subroutine get_grid_boundary_g19



   !> \brief Calculates the coordinates (x,y) of the north and south boundaries.
   !! The north boundary line is generated from the south boundary.
   !! The south boundary line is an bluff nose followed by a shifted power law ogive of
   !! length lo and radius rb perturbed with a exponential and followed by a cylinder of
   !! radius rb and length lc.
   !! The distribution of points over the south boundary is separated in three
   !! parts:
   !! Part   I: geometric progression (bluff nose+shifted power law)
   !! Part  II: geometric progression (shifted power law)
   !! Part III: geometric progression (cylinder)
   !! The distribution of points over the north boundary is  based on the normal
   !! vectors from the south boundary. The directions of the normal vectors
   !! may be changed by a smoothing procedure.
   !! Adapted to DEPP optimizer.
   subroutine get_grid_boundary_g20(nx, ny, unt, fgeom, lr, rb, x, y) ! Output: last two
      implicit none
      integer, intent(in) :: nx  !< Number of volumes in the csi direction (real+fictitious)
      integer, intent(in) :: ny  !< Number of volumes in the eta direction (real+fictitious)
      integer, intent(in) :: unt !< Unit where the input parameters will be printed
      character(len=*), intent(in) :: fgeom !< File of the geometric parameters

      real(8), intent(out) :: lr  !< Length of the body (m)
      real(8), intent(out) :: rb  !< Base radius of the body (m)
      real(8), dimension(nx*ny), intent(out) :: x !< Coord. x at the northeast corner of the volume P (m)
      real(8), dimension(nx*ny), intent(out) :: y !< Coord. y at the northeast corner of the volume P (m)

      ! Parameters

      integer, parameter :: nr = 10000 ! Number of partitions for the representing curve (xr,yr)


      ! Inner variables
      integer :: i, j, np ! Dummy indexes

      integer :: kgb ! Kind of geometry of the body
      integer :: nx1 ! Number of points in the first part of the grid along the eta line
      integer :: nx2 ! Number of points in the first and second parts of the grid along the eta lines
      integer :: nx3 ! Number of points in the first to third parts of the grid along the eta lines
      integer :: iiv ! Initial value of i
      integer :: ifv ! Final value of i
      integer :: nr1 ! Fraction of nr used to represent the Part I

      real(8) :: a2  ! Width of the first volume of Part II
      real(8) :: a3  ! Width of the first volume of Part III

      real(8) :: xiv ! Initial value of x
      real(8) :: xfv ! Final value of x

      real(8) :: lbd ! Exponent of the power law ogive

      real(8) :: lo  ! Ogive length (m)
      real(8) :: w   ! Distance between south and north boundaries (m)
      real(8) :: fs  ! Factor of smoothing (0<=fs) (necessary for non-smooth surfaces)
      real(8) :: lc  ! Cylinder length (m)
      real(8) :: lo1 ! Length of the first part of the ogive (m)
      real(8) :: rf  ! Radius of the nose tip (m)
      real(8) :: h   ! Ratio rf/rb
      real(8) :: S0  ! Maximum perturbation
      real(8) :: ae  ! Coefficient of the exponential
      real(8) :: aks ! Exponent of the power law for the calculation of the width of the first volume
      real(8) :: flo ! Fraction of lo
      real(8) :: fp1 ! Fraction of points in the first part of the grid
      real(8) :: fp2 ! Fraction of points in the second part of the grid


      real(8), dimension(nx) :: xbs ! Distribution of points over the south boundary
      real(8), dimension(nx) :: ybs ! Distribution of points over the south boundary

      real(8), dimension(0:nr) :: xr ! Distribution of points over the representing curve
      real(8), dimension(0:nr) :: yr ! Distribution of points over the representing curve

      ! DEPP variables
      integer            :: nu     ! number of unknowns
      integer            :: ind    ! number of the individual
      integer,       allocatable :: xopt(:)   ! Optimization checker
      character(10), allocatable :: xname(:)  ! Name of parameters
      real(8),       allocatable :: xvalue(:) ! parameters
      character(200)             :: sname     ! simulation name


      ! Number of parameters to read (or DEPP)

      nu = 4

      ! Reading parameters

      open(13,file=fgeom)

      read(13,*) kgb  ! Kind of geometry of the body
      read(13,*) w    ! Distance between south and north boundaries (m)
      read(13,*) fs   ! Factor of smoothing (0<=fs) (necessary for non-smooth surfaces)
      read(13,*) lo   ! Length of the ogive (m)
      read(13,*) lc   ! Length of the cylinder (m)
      read(13,*) rb   ! Base radius of the body (m)
      read(13,*) aks  ! Exponent of the power law for the south boundary
      read(13,*) flo  ! Fraction of lo
      read(13,*) fp1  ! Fraction of volumes in the first part
      read(13,*) fp2  ! Fraction of volumes in the second part

      if ( nu /= 4 ) then

         write(*,*) "nu must be 4. Stopping."

         stop

      end if

      allocate ( xopt(nu), xname(nu), xvalue(nu) )

      do i = 1, nu

         read(13,*) xname(i) ! Parameter name
         read(13,*) xopt(i)  ! Will this parameter be optimized? ( 0 = no, 1 = yes )
         read(13,*) xvalue(i)! Parameter value

      end do

      close(13)


      ! Reads the parameters from DEPP
      call depp_get_parameters(nu, xopt, xname, xvalue, ind, sname)

      h = xvalue(1) ! Radius ratio

      lbd = xvalue(2) ! Exponent of the power law

      S0 = xvalue(3) ! Maximum perturbation

      ae = xvalue(4) ! Coefficient for the exponential


      ! Gets the Shifted Power Law With Exponential Perturbation geometry
      call get_splwep(nr, 1.5d0, rb, lo, h, lbd, S0, ae, xr, yr) ! Output: last two



      ! Writing parameters

      write(unt,*)
      write(unt,*) " *** Geometry parameters ***"
      write(unt,*)

      write(unt,"(I23,    A)") kgb, " = kgb: Kind of geometry of the body"
      write(unt,"(ES23.16,A)")   w,  " =   w: Distance between south and north boundaries (m)"
      write(unt,"(ES23.16,A)")  fs, " =  fs: Factor of smoothing (0<=fs) (necessary for non-smooth surfaces)"
      write(unt,"(ES23.16,A)")  lo, " =  lo: Length of the ogive (m)"
      write(unt,"(ES23.16,A)")  lc, " =  lc: Length of the cylinder (m)"
      write(unt,"(ES23.16,A)")  rb, " =  rb: Base radius of the body (m)"
      write(unt,"(ES23.16,A)")   h, " =   h: Ratio rf/rb (frontal to base radius ratio)"
      write(unt,"(ES23.16,A)") lbd, " = lbd: Exponent of the power law ogive"
      write(unt,"(ES23.16,A)")  S0, " =  S0: Maximum perturbation"
      write(unt,"(ES23.16,A)")  ae, " =   a: Coefficient of the exponential"
      write(unt,"(ES23.16,A)") aks, " = aks: Exponent of the power law for the south boundary"
      write(unt,"(ES23.16,A)") flo, " = flo: Fraction of lo"
      write(unt,"(ES23.16,A)") fp1, " = fp1: Fraction of volumes in the first part of the grid"
      write(unt,"(ES23.16,A)") fp2, " = fp2: Fraction of volumes in the second part of the grid"

      ! Calculating the number of points in the first part of the grid

      nx1 = int( (nx-2) * fp1 )

      ! Calculating the number of points in the first and second parts of the grid

      nx2 = int( (nx-2) * fp2 ) + nx1

      ! Calculating the number of points in the 1st, 2nd and 3rd part of the grid

      nx3 = nx-1

      ! Calculating frontal radius

      rf = h * rb

      ! Calculating the fraction of lo used in the 2nd part of the grid

      lo1 = flo * lo

      do i = 1, nr

         if ( lo1 < xr(i) ) exit

      end do

      nr1 = i

      lo1 = xr(nr1)


      ! Calculating the body length

      lr = lo + lc


      ! ========================================================

      ! Generating south boundary

      j = 1

      ! Part I

      ! Calculates the distribution of points based on the arclength and the
      ! double exponential distribution

      iiv = 1   ! Initial value of i

      ifv = nx1 ! Final value of i

      call get_dexp_m1_ditribution( nr1, ifv-iiv, aks, rf, xr(0:nr1), yr(0:nr1)&
         , xbs(iiv:ifv), ybs(iiv:ifv)) ! Output: last two

      do i = iiv, ifv

         np   = nx * (j-1) + i

         x(np) = xbs(i)

         y(np) = ybs(i)

      end do


      ! Part II

      iiv = nx1  ! Initial value of i

      ifv = nx2 ! Final value of i

      xiv = lo1 ! Initial value of x

      xfv = lo  ! Final value of x

      ! Calculating the width of the first volume of Part II

      i = iiv

      a2 = sqrt( (xbs(i)-xbs(i-1))**2 + (ybs(i)-ybs(i-1))**2 )

      ! Calculating the distribution curve according to a geometric progression
      ! distribution based on the curve length

      call get_gp_length_distribution(nr-nr1, ifv-iiv, a2, -1.d0, xr(nr1:nr) &
         , yr(nr1:nr), xbs(iiv:ifv), ybs(iiv:ifv)) ! Output: last two

      do i = iiv, ifv

         np   = nx * (j-1) + i

         x(np) = xbs(i)

         y(np) = ybs(i)

      end do


      ! Part III

      iiv = nx2  ! Initial value of i

      ifv = nx3  ! Final value of i

      xiv = lo ! Initial value of x

      xfv = lo+lc  ! Final value of x

      ! Calculating the width of the first volume of Part IV

      i = iiv

      a3 = sqrt( (xbs(i)-xbs(i-1))**2 + (ybs(i)-ybs(i-1))**2 )

      call get_gp_distribution(iiv, ifv, xiv, xfv, a3, xbs(iiv:ifv)) ! Output: last one

      do i = iiv, ifv

         np   = nx * (j-1) + i

         ybs(i) = rb

         x(np) = xbs(i)

         y(np) = ybs(i)

      end do


      ! ========================================================

      ! Generating north boundary

      call get_north_boundary(nx, ny, w, fs, xbs, ybs, x, y) ! Output: last two

      ! ========================================================

   contains


      !> \brief Generates the north boundary of the domain based on the normal
      !! vectors from the south boundary. The directions of the normal vectors
      !! may be changed by a smoothing procedure.
      subroutine get_north_boundary(nx, ny, w, fs, xbs, ybs, x, y) ! Output: last two
         implicit none
         integer, intent(in)    :: nx  !< Number of volumes in the csi direction (real+fictitious)
         integer, intent(in)    :: ny  !< Number of volumes in the eta direction (real+fictitious)
         real(8), intent(in)    :: w   !< Width of the domain
         real(8), intent(in)    :: fs  !< Factor of smoothing (0<=fs)
         real(8), dimension(nx),    intent(in)  :: xbs !< Distribution of points over the south boundary
         real(8), dimension(nx),    intent(in)  :: ybs !< Distribution of points over the south boundary
         real(8), dimension(nx*ny), intent(out) :: x   !< Coord. x at the northeast corner of the volume P (m)
         real(8), dimension(nx*ny), intent(out) :: y   !< Coord. y at the northeast corner of the volume P (m)

         ! Inner variables

         integer :: i, j, k, np, ns
         real(8) :: raux
         real(8) :: tx
         real(8) :: ty
         real(8) :: tm
         real(8) :: vnx(nx-1)
         real(8) :: vny(nx-1)

         ! Number of iteractions for smoothing the direction vectors
         ns = int(nx * fs)

         ! Calculating the initial normal directions

         ! Direction of the west boundary
         vnx(1) = -1.d0
         vny(1) =  0.d0


         ! Direction of the east boundary
         vnx(nx-1) = 0.d0
         vny(nx-1) = 1.d0

         ! Directions normal to the body surface
         do i = 2, nx-2

            tx = ( xbs(i+1) - xbs(i-1) ) / 2.d0
            ty = ( ybs(i+1) - ybs(i-1) ) / 2.d0

            tm = sqrt( tx * tx + ty * ty )

            vnx(i) = - ty / tm
            vny(i) =   tx / tm

         end do

         ! Calculating the final directions (smoothing)

         do k = 1, ns

            do i = 2, nx-2

               vnx(i) = ( vnx(i+1) + vnx(i-1) ) / 2.d0

               vny(i) = ( vny(i+1) + vny(i-1) ) / 2.d0

               raux = dsqrt( vnx(i) * vnx(i) + vny(i) * vny(i) )

               vnx(i) = vnx(i) / raux

               vny(i) = vny(i) / raux

            end do

         end do

         ! Generating the north boundary

         j = ny-1

         do i = 1, nx-1

            np   = nx * (j-1) + i

            x(np) = xbs(i) + w * vnx(i)

            y(np) = ybs(i) + w * vny(i)

         end do

      end subroutine get_north_boundary


      !> \brief Calculates the node distribution based on a geometric progression distribution
      subroutine get_gp_distribution(iiv, ifv, xiv, xfv, a1, xb) ! Output: last one
         implicit none
         integer, intent(in)  :: iiv !< Initial value of i
         integer, intent(in)  :: ifv !< Final value of i
         real(8), intent(in)  :: xiv !< Initial value of x
         real(8), intent(in)  :: xfv !< Final value of x
         real(8), intent(in)  :: a1  !< Width of the first volume
         real(8), intent(out) :: xb(iiv:ifv) !< Distribution of points

         ! Inner variables

         integer :: i
         integer :: n
         real(8) :: r
         real(8) :: q


         n = ifv - iiv

         r = ( xfv - xiv ) / a1

         ! Calculating q
         call get_gp_ratio(n, r, q)

         xb(iiv) = xiv

         do i = iiv+1, ifv-1

            xb(i) = xb(i-1) + a1 * q ** ( i - iiv - 1 )

         end do

         xb(ifv) = xfv

      end subroutine get_gp_distribution


      !> \brief Calculates the base q of the geometric progression
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

           if ( 0.d0 < f(n, r, qi) * f(n, r, qm) ) then

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

      end subroutine get_gp_ratio

      !> \brief Function used in the subroutine 'get_gp_ratio'
      real(8) function f(n, r, q)
         implicit none
         integer, intent(in)  ::   n !< number of partitions
         real(8), intent(in)  ::   r !< l/a1
         real(8), intent(in)  ::   q !< q

         f = q ** n + r * ( 1.d0 - q ) - 1.d0

      end function f

      !> \brief Given a reference distribution of points (xbr,ybr)
      !! a new one (xb,yb) is created based on a geometric progression distributiuon
      !! based on the curve length.
      subroutine get_gp_length_distribution( nbr, nb, al, ar, xbr, ybr, xb, yb) ! Output: last two
         implicit none
         integer, intent(in)  :: nbr        !< Number of partitions over the reference curve (xbr,ybr)
         integer, intent(in)  :: nb         !< Number of partitions for the new distribution of points
         real(8), intent(in)  :: al         !< Length of the first partition on the left
         real(8), intent(in)  :: ar         !< Length of the first partition on the right
         real(8), intent(in)  :: xbr(0:nbr) !< Given distribution of points
         real(8), intent(in)  :: ybr(0:nbr) !< Given distribution of points
         real(8), intent(out) :: xb(0:nb)   !< New distribution of points
         real(8), intent(out) :: yb(0:nb)   !< New distribution of points

         ! Inner variables

         real(8) :: sbr(0:nbr) ! Length of each partition over the reference curve
         real(8) :: sb(0:nb)   ! Length of each partition over the new curve
         real(8) :: L          ! Total length


         ! Calculates the length of the curve
         call get_sbr(nbr, xbr, ybr, sbr) ! Output: last one


         ! Total length
         L = sbr(nbr)


         ! Selecting the kind of distribution

         if ( al < 0.d0 .and. ar < 0.d0 ) then

            ! Calculates the new points according to a uniform distribution
            call get_uniform_grid( nb, L, sb) ! Output: last one

         else if ( al > 0.d0 .and. ar < 0.d0 ) then

            ! Calculates the distribution of points according to a geometric
            ! progression. The length of the first partition is defined on the left.
            call get_gp_grid_left( nb, al, L, sb)

         else if ( al < 0.d0 .and. ar > 0.d0 ) then

            ! Calculates the distribution of points according to a geometric
            ! progression. The length of the first partition is defined on the right.
            call get_gp_grid_right( nb, ar, L, sb)

         else

            ! Calculates the distribution of points according to a geometric
            ! progression. The length of the first partition is defined on the left
            ! and on the right.
            call get_gp_grid_left_right(nb, al, ar, L, sb)

         end if


         ! Calculates the new distribution of points
         call get_new_distribution(nbr, nb, sbr, xbr, ybr, sb, xb, yb) ! Output: last two


      end subroutine get_gp_length_distribution

      !> \brief Calculates the length of the curve
      subroutine get_sbr(nbr, xbr, ybr, sbr) ! Output: last one
         implicit none
         integer, intent(in)  :: nbr        !< Number of partitions over the reference curve (xbr,ybr)
         real(8), intent(in)  :: xbr(0:nbr) !< Given distribution of points
         real(8), intent(in)  :: ybr(0:nbr) !< Given distribution of points
         real(8), intent(out) :: sbr(0:nbr) !< Length of the reference curve

         ! Inner variables

         integer :: i


         sbr(0) = 0.d0

         do i = 1, nbr

            sbr(i) = sbr(i-1) + sqrt( ( xbr(i) - xbr(i-1) )**2 + ( ybr(i) - ybr(i-1) )**2 )

         end do


      end subroutine get_sbr


      !> \brief Given the length distribution, the (xb,yb) is calculated
      subroutine get_new_distribution(nbr, nb, sbr, xbr, ybr, sb, xb, yb) ! Output: last two
         integer, intent(in)  :: nbr        !< Number of partitions over the reference curve (xbr,ybr)
         integer, intent(in)  :: nb         !< Number of partitions for the new distribution of points
         real(8), intent(in)  :: sbr(0:nbr) !< Length of the curve
         real(8), intent(in)  :: xbr(0:nbr) !< Given distribution of points
         real(8), intent(in)  :: ybr(0:nbr) !< Given distribution of points
         real(8), intent(in)  :: sb(0:nb)   !< New length distribution over the curve
         real(8), intent(out) :: xb(0:nb)   !< New distribution of points
         real(8), intent(out) :: yb(0:nb)   !< New distribution of points


         ! Inner variables

         integer :: i, j

         xb(0) = xbr(0)
         yb(0) = ybr(0)

         j = 0

         do i = 1, nb-1

            do j = j, nbr

               if( sbr(j) - sb(i) >= 0.d0 ) exit

            end do

            ! Linear interpolation

            xb(i) = xbr(j-1) + ( xbr(j)-xbr(j-1) ) * ( sb(i) - sbr(j-1) ) / ( sbr(j) - sbr(j-1) )
            yb(i) = ybr(j-1) + ( ybr(j)-ybr(j-1) ) * ( sb(i) - sbr(j-1) ) / ( sbr(j) - sbr(j-1) )

         end do

         xb(nb) = xbr(nbr)
         yb(nb) = ybr(nbr)

      end subroutine get_new_distribution



      !> \brief Calculates the uniform distribution of points
      subroutine get_uniform_grid( nb, L, sb) ! Output: last one
         implicit none
         integer, intent(in)  :: nb       !< Number of partitions over the new distribution of points
         real(8), intent(in)  :: L        !< Length of the curve
         real(8), intent(out) :: sb(0:nb) !< Length distribution over the curve

         ! Inner variables

         integer :: i

         real(8) :: ds

         ds = L / dble(nb)


         sb(0) = 0.d0

         do i = 1, nb

            sb(i) = sb(i-1) + ds

         end do

      end subroutine get_uniform_grid



      !> \brief Calculates the distribution of points according
      !! to a geometric progression. The length of the first partition
      !! is defined on the left.
      subroutine get_gp_grid_left(nb, al, L, sb)
         implicit none
         integer, intent(in)  :: nb       !< Number of partitions over the new distribution of points
         real(8), intent(in)  :: al       !< Length of the first patition on the left
         real(8), intent(in)  :: L        !< Length of the curve
         real(8), intent(out) :: sb(0:nb) !< Length distribution over the curve

         ! Inner variables

         integer :: i
         real(8) :: q

         ! Calculates the ratio q of the geometric progression
         call get_gp_ratio(nb, L/al, q)

         sb(0) = 0.d0

         do i = 1, nb-1

            sb(i) = sb(i-1) + al * q ** (i-1)

         end do

         sb(nb) = L

      end subroutine get_gp_grid_left


      !> \brief Calculates the distribution of points according
      !! to a geometric progression. The length of the first partition
      !! is defined on the right.
      subroutine get_gp_grid_right(nb, ar, L, sb)
         implicit none
         integer, intent(in)  :: nb       !< Number of partitions over the new distribution of points
         real(8), intent(in)  :: ar       !< Length of the first patition on the right
         real(8), intent(in)  :: L        !< Length of the curve
         real(8), intent(out) :: sb(0:nb) !< Length distribution over the curve

         ! Inner variables

         integer :: i
         real(8) :: q

         ! Calculates the ratio q of the geometric progression
         call get_gp_ratio(nb, L/ar, q)

         sb(0) = 0.d0

         do i = 1, nb-1

            sb(i) = sb(i-1) + ar * q ** (nb-i)

         end do

         sb(nb) = L

      end subroutine get_gp_grid_right



      !> \brief Calculates the distribution of points according
      !! to a geometric progression. The length of the first partition
      !! is defined on the left and on the right.
      subroutine get_gp_grid_left_right(nb, al, ar, L, sb)
         implicit none
         integer, intent(in)  :: nb       !< Number of partitions over the new distribution of points
         real(8), intent(in)  :: al       !< Length of the first patition on the left
         real(8), intent(in)  :: ar       !< Length of the first patition on the right
         real(8), intent(in)  :: L        !< Length of the curve
         real(8), intent(out) :: sb(0:nb) !< Length distribution over the curve

         ! Inner parameters

         integer, parameter :: nitmax = 10000

         ! Inner variables

         integer :: i
         integer :: nl
         integer :: nr
         real(8) :: Lm
         real(8) :: Li
         real(8) :: Lf
         real(8) :: acl
         real(8) :: acr


         nl = nb / 2

         nr = nb - nl

         Li = 0.d0

         Lf = L

         do i = 1, nitmax

            Lm = 0.5d0 * ( Li + Lf )

            call get_gp_grid_left (nl, al, Lm, sb( 0:nl))

            call get_gp_grid_right(nr, ar, L-Lm, sb(nl:nb))

            sb(nl:nb) = sb(nl:nb) + Lm

            acl = sb(nl) - sb(nl-1)

            acr = sb(nl+1) - sb(nl)

            if ( acl > acr ) then

               Lf = Lm

            else

               Li = Lm

            end if

            if ( abs(acl-acr) < 1.d-14 ) exit

         end do

      end subroutine get_gp_grid_left_right


      !> \brief Calculates a distribution of points based on the arclength and
      !! on the double-exponential rule
      subroutine get_dexp_m1_ditribution( nbr, nb, aks, zcv, xbr, ybr, xb, yb) ! Output: last two
         implicit none
         integer, intent(in)  :: nbr        !< Number of partitions over the reference curve (xbr,ybr)
         integer, intent(in)  :: nb         !< Number of partitions for the new distribution of points
         real(8), intent(in)  :: aks        !< Exponent for the calculation of acv
         real(8), intent(in)  :: zcv        !< Central value of z
         real(8), intent(in)  :: xbr(0:nbr) !< Given distribution of points
         real(8), intent(in)  :: ybr(0:nbr) !< Given distribution of points
         real(8), intent(out) :: xb(0:nb)   !< New distribution of points
         real(8), intent(out) :: yb(0:nb)   !< New distribution of points

         ! Inner variables

         real(8) :: sbr(0:nbr) ! Length of each partition over the reference curve
         real(8) :: sb(0:nb)   ! Length of each partition over the new curve
         real(8) :: L          ! Total length
         real(8) :: acv        ! Length of the volume closer to the central point


         ! Calculates the length of the curve
         call get_sbr(nbr, xbr, ybr, sbr) ! Output: last one


         ! Total length
         L = sbr(nbr)


         ! Calculating the double-exponential distribution of arclength

         acv = L / ( nb ** aks )

         if ( acv > zcv ) then

            write(*,*) "Error. get_dexp_m1_ditribution: there is no resolution"&
               , " for this grid. Stopping."

            stop

         end if

         call get_dexp_m1(nb, 0.d0, zcv, L, acv, sb)


         ! Calculates the new distribution of points
         call get_new_distribution(nbr, nb, sbr, xbr, ybr, sb, xb, yb) ! Output: last two


      end subroutine get_dexp_m1_ditribution

      !> \brief Calculates a distribution of points based on the double-exponential
      !! rule. The coefficients of this rule are set by defining the position of
      !! an internal point (zcv) and the width (acv) of the volumes beside this
      !! point. The points are concentrated around zcv.
      subroutine get_dexp_m1(np, ziv, zcv, zfv, acv, z)
         implicit none
         integer, intent(in)  :: np  !< Number of partitions
         real(8), intent(in)  :: ziv !< Initial value of z
         real(8), intent(in)  :: zcv !< Central value of z (exponential's matching point)
         real(8), intent(in)  :: zfv !< Final value of z
         real(8), intent(in)  :: acv !< Length of the volumes beside zcv
         real(8), intent(out) :: z(0:np) !< Distribution of points

         ! Inner variables

         integer :: i
         integer :: npc
         real(8) :: t
         real(8) :: dt
         real(8) :: tc
         real(8) :: r
         real(8) :: A
         real(8) :: ac
         real(8) :: zc


         ! Initializing variables

         zc = zcv / ( zfv - ziv )

         ac = acv / ( zfv - ziv )

         dt = 1.d0 / dble(np)

         r = ac / ( 1.d0 - zc )



         ! Cheking conditions

         if ( ( ac / zc > 1.d0 ) .or. ( r > 1.d0 ) ) then

            write(*,*) "Error. get_dexp_m1: ac > zc or ac > (1-zc). Stopping."

            stop

         end if

         if ( ac > dt ) then

            write(*,*) "Error. get_dexp_m1: ac > dt. Stopping."

            stop

         end if


         ! Searching the value of A

         call adexp_m1(dt, r, zc, A)


         ! Calculating the fraction of volumes for each exponential

         tc = 1.d0 + log( 1.d0 + ( exp(-A*dt) - 1.d0 ) / r ) / A

         npc = nint(tc*np)



         ! Calculating the normallized distribution of points

         ! First exponential

         do i = 0, npc

            t = tc * dble(i) / dble(npc)

            z(i) = zc * ( exp(A*t) - 1.d0 ) / ( exp(A*tc) - 1.d0 )

         end do

         ! Second exponential

         do i = npc, np

            t = ( 1.d0-tc ) * dble(i-npc) / dble(np-npc) + tc

            z(i) = zc+(1.d0-zc)*(exp(-A*(t-tc))-1.d0)/(exp(-A*(1.d0-tc))-1.d0)

         end do

         z(0)  = 0.d0

         z(np) = 1.d0

         ! Calculating the extended distribution of points

         z = ( zfv - ziv ) * z + ziv

      end subroutine get_dexp_m1

      !> \brief Defines the transcendental function for finding A
      real(8) function fdexp_m1(dt, r, zc, A)
         implicit none
         real(8), intent(in) :: dt
         real(8), intent(in) :: r
         real(8), intent(in) :: zc
         real(8), intent(in) :: A

         ! Inner varialbles

         real(8) :: raux1, raux2, raux3

         if ( abs(A) > 1.d-6 ) then

            raux1 = exp(A) * ( 1.d0 + ( exp(-A*dt) - 1.d0 ) / r )
            raux2 = (2.d0*zc-1.d0) / zc
            raux3 = (1.d0-zc) / zc

            fdexp_m1 = -exp(-A) * raux1 * raux1 +  raux1 * raux2 + raux3

         else

            fdexp_m1 = (r*zc-r+dt)*A/(r*zc) &
               + ((r**2-2*dt**2)*zc-r**2+(2*dt-dt**2)*r)*A**2/(r**2*zc)/2.0

         end if

      end function fdexp_m1


      !> \brief Finds the value of A using the bissection method
      subroutine adexp_m1(dt, r, zc, A)
         implicit none
         real(8), intent(in)  :: dt
         real(8), intent(in)  :: r
         real(8), intent(in)  :: zc
         real(8), intent(out) :: A

         ! Parameters

         integer, parameter :: itmax = 1000

         ! Inner variables

         integer :: k
         real(8) :: Ai, Af, Am, Fi, Fm, Ff



         ! Initializing the search interval

         Ai = -1.d-18

         Af = -100.d0



         ! Checking for roots

         Fi = fdexp_m1(dt, r, zc, Ai)

         Ff = fdexp_m1(dt, r, zc, Af)

         if ( Fi * Ff > 0.d0 ) then

            write(*,*) "Error. adexp_m1: there is none or more than one roots. Stopping."

            stop

         end if


         ! Looking for a solution

         do k = 1, itmax

            Am = 0.5d0 * ( Ai + Af )

            Fi = fdexp_m1(dt, r, zc, Ai)

            Fm = fdexp_m1(dt, r, zc, Am)

            if ( Fi * Fm <= 0.d0 ) then

               Af = Am

            else

               Ai = Am

            end if

            if ( abs(Af-Ai) < 1.d-14 ) exit

         end do

         if ( k >= itmax ) then

            write(*,*) "Error. adexp_m1: number of iteractions exceeded. Stopping."

            stop

         end if

         A = 0.5d0 * ( Ai + Af )

      end subroutine adexp_m1


   end subroutine get_grid_boundary_g20



   !> \brief Calculates the coordinates (x,y) of the north and south boundaries.
   !! The north boundary line is generated from the south boundary.
   !! The south boundary line is an bluff nose followed by a smooth arc generated
   !! by natural cubic splines of length lo and radius rb, followed by a cylinder
   !! of radius rb and length lc.
   !! The distribution of points over the south boundary is separated in three
   !! parts:
   !! Part   I: geometric progression (bluff nose+smooth arc)
   !! Part  II: geometric progression (smooth arc)
   !! Part III: geometric progression (cylinder)
   !! The distribution of points over the north boundary is  based on the normal
   !! vectors from the south boundary. The directions of the normal vectors
   !! may be changed by a smoothing procedure.
   !! Adapted to DEPP optimizer.
   subroutine get_grid_boundary_g21(nx, ny, unt, fgeom, lr, rb, x, y) ! Output: last two
      implicit none
      integer, intent(in) :: nx  !< Number of volumes in the csi direction (real+fictitious)
      integer, intent(in) :: ny  !< Number of volumes in the eta direction (real+fictitious)
      integer, intent(in) :: unt !< Unit where the input parameters will be printed
      character(len=*), intent(in) :: fgeom !< File of the geometric parameters

      real(8), intent(out) :: lr  !< Length of the body (m)
      real(8), intent(out) :: rb  !< Base radius of the body (m)
      real(8), dimension(nx*ny), intent(out) :: x !< Coord. x at the northeast corner of the volume P (m)
      real(8), dimension(nx*ny), intent(out) :: y !< Coord. y at the northeast corner of the volume P (m)

      ! Parameters

      integer, parameter :: nr = 10000 ! Number of partitions for the representing curve (xr,yr)


      ! Inner variables
      integer :: i, j, np ! Dummy indexes
      integer :: es  ! Exit status ( 0=success, 1=failure )
      integer :: mc  ! Monocity check ( 0=no, 1=yes )
      integer :: nsi ! Number of points for the spline interpolation

      integer :: kgb ! Kind of geometry of the body
      integer :: nx1 ! Number of points in the first part of the grid along the eta line
      integer :: nx2 ! Number of points in the first and second parts of the grid along the eta lines
      integer :: nx3 ! Number of points in the first to third parts of the grid along the eta lines
      integer :: iiv ! Initial value of i
      integer :: ifv ! Final value of i
      integer :: nr1 ! Fraction of nr used to represent the Part I

      real(8) :: a2  ! Width of the first volume of Part II
      real(8) :: a3  ! Width of the first volume of Part III

      real(8) :: xiv ! Initial value of x
      real(8) :: xfv ! Final value of x

      real(8) :: lo  ! Ogive length (m)
      real(8) :: w   ! Distance between south and north boundaries (m)
      real(8) :: fs  ! Factor of smoothing (0<=fs) (necessary for non-smooth surfaces)
      real(8) :: lc  ! Cylinder length (m)
      real(8) :: lo1 ! Length of the first part of the ogive (m)
      real(8) :: rf  ! Radius of the nose tip (m)
      real(8) :: aks ! Exponent of the power law for the calculation of the width of the first volume
      real(8) :: flo ! Fraction of lo
      real(8) :: fp1 ! Fraction of points in the first part of the grid
      real(8) :: fp2 ! Fraction of points in the second part of the grid


      real(8), dimension(nx) :: xbs ! Distribution of points over the south boundary
      real(8), dimension(nx) :: ybs ! Distribution of points over the south boundary

      real(8), dimension(0:nr) :: xr ! Distribution of points over the representing curve
      real(8), dimension(0:nr) :: yr ! Distribution of points over the representing curve

      real(8), allocatable :: xsi(:) ! x/lr Data for spline interpolation
      real(8), allocatable :: ysi(:) ! y/rb Data for spline interpolation

      ! DEPP variables
      integer            :: nu     ! number of unknowns
      integer            :: ind    ! number of the individual
      integer,       allocatable :: xopt(:)   ! Optimization checker
      character(10), allocatable :: xname(:)  ! Name of parameters
      real(8),       allocatable :: xvalue(:) ! parameters
      character(200)             :: sname     ! simulation name


      ! Reading parameters

      open(13,file=fgeom)

      read(13,*) kgb  ! Kind of geometry of the body
      read(13,*) w    ! Distance between south and north boundaries (m)
      read(13,*) fs   ! Factor of smoothing (0<=fs) (necessary for non-smooth surfaces)
      read(13,*) lo   ! Length of the ogive (m)
      read(13,*) lc   ! Length of the cylinder (m)
      read(13,*) rb   ! Base radius of the body (m)
      read(13,*) aks  ! Exponent of the power law for the south boundary
      read(13,*) flo  ! Fraction of lo
      read(13,*) fp1  ! Fraction of volumes in the first part
      read(13,*) fp2  ! Fraction of volumes in the second part
      read(13,*) mc   ! Monocity check ( 0=no, 1=yes )
      read(13,*) nu   ! Number of unknowns

      if ( nu < 1 ) then

         write(*,*) "nu < 1. Stopping."

         stop

      end if

      nsi = nu

      allocate ( xopt(nu), xname(nu), xvalue(nu) )

      allocate( xsi(nsi), ysi(nsi) )

      do i = 1, nu

         read(13,*) xsi(i)

      end do

      do i = 1, nu

         read(13,*) xname(i) ! Parameter name
         read(13,*) xopt(i)  ! Will this parameter be optimized? ( 0 = no, 1 = yes )
         read(13,*) xvalue(i)! Parameter value

      end do

      close(13)


      ! Reads the parameters from DEPP
      call depp_get_parameters(nu, xopt, xname, xvalue, ind, sname)

      ysi = xvalue


      ! Gets the Natural Cubic Spline Interpolation geometry
      call get_ncspline_geometry(nr, nsi, 1.5d0, rb, lo, xsi, ysi, xr, yr, es) ! Output: last three


      ! Checking function monotonicity
      if ( mc == 1 ) then

         ! Non-monotonic function. Requesting new parameters.
         if ( es == 1 ) then

            call depp_save_fitness(-huge(1.d0), 2, "-cdfp")

            write(*,*) "Non-monotonic geometry. Stopping."

            stop

         end if

      end if

      ! Writing parameters

      write(unt,*)
      write(unt,*) " *** Geometry parameters ***"
      write(unt,*)

      write(unt,"(I23,    A)") kgb, " = kgb: Kind of geometry of the body"
      write(unt,"(ES23.16,A)")   w, " =   w: Distance between south and north boundaries (m)"
      write(unt,"(ES23.16,A)")  fs, " =  fs: Factor of smoothing (0<=fs) (necessary for non-smooth surfaces)"
      write(unt,"(ES23.16,A)")  lo, " =  lo: Length of the ogive (m)"
      write(unt,"(ES23.16,A)")  lc, " =  lc: Length of the cylinder (m)"
      write(unt,"(ES23.16,A)")  rb, " =  rb: Base radius of the body (m)"
      write(unt,"(ES23.16,A)") aks, " = aks: Exponent of the power law for the south boundary"
      write(unt,"(ES23.16,A)") flo, " = flo: Fraction of lo"
      write(unt,"(ES23.16,A)") fp1, " = fp1: Fraction of volumes in the first part of the grid"
      write(unt,"(ES23.16,A)") fp2, " = fp2: Fraction of volumes in the second part of the grid"
      write(unt,"(I23,    A)")  mc, " =  mc: Monotonicity check"

      do i = 1, nu

         write(unt,"(ES23.16,A)") xsi(i), ' = xsi: x/lr for the spline calculation'

      end do

      do i = 1, nu

         write(unt,"(ES23.16,A)") ysi(i), ' = ysi: y/rb for the spline calculation'

      end do


      ! Calculating the number of points in the first part of the grid

      nx1 = int( (nx-2) * fp1 )

      ! Calculating the number of points in the first and second parts of the grid

      nx2 = int( (nx-2) * fp2 ) + nx1

      ! Calculating the number of points in the 1st, 2nd and 3rd part of the grid

      nx3 = nx-1

      ! Calculating frontal radius

      rf = ysi(1) * rb

      ! Calculating the fraction of lo used in the 2nd part of the grid

      lo1 = flo * lo

      do i = 1, nr

         if ( lo1 < xr(i) ) exit

      end do

      nr1 = i

      lo1 = xr(nr1)


      ! Calculating the body length

      lr = lo + lc


      ! ========================================================

      ! Generating south boundary

      j = 1

      ! Part I

      ! Calculates the distribution of points based on the arclength and the
      ! double exponential distribution

      iiv = 1   ! Initial value of i

      ifv = nx1 ! Final value of i

      call get_dexp_m1_ditribution( nr1, ifv-iiv, aks, rf, xr(0:nr1), yr(0:nr1)&
         , xbs(iiv:ifv), ybs(iiv:ifv)) ! Output: last two

      do i = iiv, ifv

         np   = nx * (j-1) + i

         x(np) = xbs(i)

         y(np) = ybs(i)

      end do


      ! Part II

      iiv = nx1  ! Initial value of i

      ifv = nx2 ! Final value of i

      xiv = lo1 ! Initial value of x

      xfv = lo  ! Final value of x

      ! Calculating the width of the first volume of Part II

      i = iiv

      a2 = sqrt( (xbs(i)-xbs(i-1))**2 + (ybs(i)-ybs(i-1))**2 )

      ! Calculating the distribution curve according to a geometric progression
      ! distribution based on the curve length

      call get_gp_length_distribution(nr-nr1, ifv-iiv, a2, -1.d0, xr(nr1:nr) &
         , yr(nr1:nr), xbs(iiv:ifv), ybs(iiv:ifv)) ! Output: last two

      do i = iiv, ifv

         np   = nx * (j-1) + i

         x(np) = xbs(i)

         y(np) = ybs(i)

      end do


      ! Part III

      iiv = nx2  ! Initial value of i

      ifv = nx3  ! Final value of i

      xiv = lo ! Initial value of x

      xfv = lo+lc  ! Final value of x

      ! Calculating the width of the first volume of Part IV

      i = iiv

      a3 = sqrt( (xbs(i)-xbs(i-1))**2 + (ybs(i)-ybs(i-1))**2 )

      call get_gp_distribution(iiv, ifv, xiv, xfv, a3, xbs(iiv:ifv)) ! Output: last one

      do i = iiv, ifv

         np   = nx * (j-1) + i

         ybs(i) = rb

         x(np) = xbs(i)

         y(np) = ybs(i)

      end do


      ! ========================================================

      ! Generating north boundary

      call get_north_boundary(nx, ny, w, fs, xbs, ybs, x, y) ! Output: last two

      ! ========================================================

   contains


      !> \brief Generates the north boundary of the domain based on the normal
      !! vectors from the south boundary. The directions of the normal vectors
      !! may be changed by a smoothing procedure.
      subroutine get_north_boundary(nx, ny, w, fs, xbs, ybs, x, y) ! Output: last two
         implicit none
         integer, intent(in)    :: nx  !< Number of volumes in the csi direction (real+fictitious)
         integer, intent(in)    :: ny  !< Number of volumes in the eta direction (real+fictitious)
         real(8), intent(in)    :: w   !< Width of the domain
         real(8), intent(in)    :: fs  !< Factor of smoothing (0<=fs)
         real(8), dimension(nx),    intent(in)  :: xbs !< Distribution of points over the south boundary
         real(8), dimension(nx),    intent(in)  :: ybs !< Distribution of points over the south boundary
         real(8), dimension(nx*ny), intent(out) :: x   !< Coord. x at the northeast corner of the volume P (m)
         real(8), dimension(nx*ny), intent(out) :: y   !< Coord. y at the northeast corner of the volume P (m)

         ! Inner variables

         integer :: i, j, k, np, ns
         real(8) :: raux
         real(8) :: tx
         real(8) :: ty
         real(8) :: tm
         real(8) :: vnx(nx-1)
         real(8) :: vny(nx-1)

         ! Number of iteractions for smoothing the direction vectors
         ns = int(nx * fs)

         ! Calculating the initial normal directions

         ! Direction of the west boundary
         vnx(1) = -1.d0
         vny(1) =  0.d0


         ! Direction of the east boundary
         vnx(nx-1) = 0.d0
         vny(nx-1) = 1.d0

         ! Directions normal to the body surface
         do i = 2, nx-2

            tx = ( xbs(i+1) - xbs(i-1) ) / 2.d0
            ty = ( ybs(i+1) - ybs(i-1) ) / 2.d0

            tm = sqrt( tx * tx + ty * ty )

            vnx(i) = - ty / tm
            vny(i) =   tx / tm

         end do

         ! Calculating the final directions (smoothing)

         do k = 1, ns

            do i = 2, nx-2

               vnx(i) = ( vnx(i+1) + vnx(i-1) ) / 2.d0

               vny(i) = ( vny(i+1) + vny(i-1) ) / 2.d0

               raux = dsqrt( vnx(i) * vnx(i) + vny(i) * vny(i) )

               vnx(i) = vnx(i) / raux

               vny(i) = vny(i) / raux

            end do

         end do

         ! Generating the north boundary

         j = ny-1

         do i = 1, nx-1

            np   = nx * (j-1) + i

            x(np) = xbs(i) + w * vnx(i)

            y(np) = ybs(i) + w * vny(i)

         end do

      end subroutine get_north_boundary


      !> \brief Calculates the node distribution based on a geometric progression distribution
      subroutine get_gp_distribution(iiv, ifv, xiv, xfv, a1, xb) ! Output: last one
         implicit none
         integer, intent(in)  :: iiv !< Initial value of i
         integer, intent(in)  :: ifv !< Final value of i
         real(8), intent(in)  :: xiv !< Initial value of x
         real(8), intent(in)  :: xfv !< Final value of x
         real(8), intent(in)  :: a1  !< Width of the first volume
         real(8), intent(out) :: xb(iiv:ifv) !< Distribution of points

         ! Inner variables

         integer :: i
         integer :: n
         real(8) :: r
         real(8) :: q


         n = ifv - iiv

         r = ( xfv - xiv ) / a1

         ! Calculating q
         call get_gp_ratio(n, r, q)

         xb(iiv) = xiv

         do i = iiv+1, ifv-1

            xb(i) = xb(i-1) + a1 * q ** ( i - iiv - 1 )

         end do

         xb(ifv) = xfv

      end subroutine get_gp_distribution


      !> \brief Calculates the base q of the geometric progression
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

           if ( 0.d0 < f(n, r, qi) * f(n, r, qm) ) then

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

      end subroutine get_gp_ratio

      !> \brief Function used in the subroutine 'get_gp_ratio'
      real(8) function f(n, r, q)
         implicit none
         integer, intent(in)  ::   n !< number of partitions
         real(8), intent(in)  ::   r !< l/a1
         real(8), intent(in)  ::   q !< q

         f = q ** n + r * ( 1.d0 - q ) - 1.d0

      end function f

      !> \brief Given a reference distribution of points (xbr,ybr)
      !! a new one (xb,yb) is created based on a geometric progression distributiuon
      !! based on the curve length.
      subroutine get_gp_length_distribution( nbr, nb, al, ar, xbr, ybr, xb, yb) ! Output: last two
         implicit none
         integer, intent(in)  :: nbr        !< Number of partitions over the reference curve (xbr,ybr)
         integer, intent(in)  :: nb         !< Number of partitions for the new distribution of points
         real(8), intent(in)  :: al         !< Length of the first partition on the left
         real(8), intent(in)  :: ar         !< Length of the first partition on the right
         real(8), intent(in)  :: xbr(0:nbr) !< Given distribution of points
         real(8), intent(in)  :: ybr(0:nbr) !< Given distribution of points
         real(8), intent(out) :: xb(0:nb)   !< New distribution of points
         real(8), intent(out) :: yb(0:nb)   !< New distribution of points

         ! Inner variables

         real(8) :: sbr(0:nbr) ! Length of each partition over the reference curve
         real(8) :: sb(0:nb)   ! Length of each partition over the new curve
         real(8) :: L          ! Total length


         ! Calculates the length of the curve
         call get_sbr(nbr, xbr, ybr, sbr) ! Output: last one


         ! Total length
         L = sbr(nbr)


         ! Selecting the kind of distribution

         if ( al < 0.d0 .and. ar < 0.d0 ) then

            ! Calculates the new points according to a uniform distribution
            call get_uniform_grid( nb, L, sb) ! Output: last one

         else if ( al > 0.d0 .and. ar < 0.d0 ) then

            ! Calculates the distribution of points according to a geometric
            ! progression. The length of the first partition is defined on the left.
            call get_gp_grid_left( nb, al, L, sb)

         else if ( al < 0.d0 .and. ar > 0.d0 ) then

            ! Calculates the distribution of points according to a geometric
            ! progression. The length of the first partition is defined on the right.
            call get_gp_grid_right( nb, ar, L, sb)

         else

            ! Calculates the distribution of points according to a geometric
            ! progression. The length of the first partition is defined on the left
            ! and on the right.
            call get_gp_grid_left_right(nb, al, ar, L, sb)

         end if


         ! Calculates the new distribution of points
         call get_new_distribution(nbr, nb, sbr, xbr, ybr, sb, xb, yb) ! Output: last two


      end subroutine get_gp_length_distribution

      !> \brief Calculates the length of the curve
      subroutine get_sbr(nbr, xbr, ybr, sbr) ! Output: last one
         implicit none
         integer, intent(in)  :: nbr        !< Number of partitions over the reference curve (xbr,ybr)
         real(8), intent(in)  :: xbr(0:nbr) !< Given distribution of points
         real(8), intent(in)  :: ybr(0:nbr) !< Given distribution of points
         real(8), intent(out) :: sbr(0:nbr) !< Length of the reference curve

         ! Inner variables

         integer :: i


         sbr(0) = 0.d0

         do i = 1, nbr

            sbr(i) = sbr(i-1) + sqrt( ( xbr(i) - xbr(i-1) )**2 + ( ybr(i) - ybr(i-1) )**2 )

         end do


      end subroutine get_sbr


      !> \brief Given the length distribution, the (xb,yb) is calculated
      subroutine get_new_distribution(nbr, nb, sbr, xbr, ybr, sb, xb, yb) ! Output: last two
         integer, intent(in)  :: nbr        !< Number of partitions over the reference curve (xbr,ybr)
         integer, intent(in)  :: nb         !< Number of partitions for the new distribution of points
         real(8), intent(in)  :: sbr(0:nbr) !< Length of the curve
         real(8), intent(in)  :: xbr(0:nbr) !< Given distribution of points
         real(8), intent(in)  :: ybr(0:nbr) !< Given distribution of points
         real(8), intent(in)  :: sb(0:nb)   !< New length distribution over the curve
         real(8), intent(out) :: xb(0:nb)   !< New distribution of points
         real(8), intent(out) :: yb(0:nb)   !< New distribution of points


         ! Inner variables

         integer :: i, j

         xb(0) = xbr(0)
         yb(0) = ybr(0)

         j = 0

         do i = 1, nb-1

            do j = j, nbr

               if( sbr(j) - sb(i) >= 0.d0 ) exit

            end do

            ! Linear interpolation

            xb(i) = xbr(j-1) + ( xbr(j)-xbr(j-1) ) * ( sb(i) - sbr(j-1) ) / ( sbr(j) - sbr(j-1) )
            yb(i) = ybr(j-1) + ( ybr(j)-ybr(j-1) ) * ( sb(i) - sbr(j-1) ) / ( sbr(j) - sbr(j-1) )

         end do

         xb(nb) = xbr(nbr)
         yb(nb) = ybr(nbr)

      end subroutine get_new_distribution



      !> \brief Calculates the uniform distribution of points
      subroutine get_uniform_grid( nb, L, sb) ! Output: last one
         implicit none
         integer, intent(in)  :: nb       !< Number of partitions over the new distribution of points
         real(8), intent(in)  :: L        !< Length of the curve
         real(8), intent(out) :: sb(0:nb) !< Length distribution over the curve

         ! Inner variables

         integer :: i

         real(8) :: ds

         ds = L / dble(nb)


         sb(0) = 0.d0

         do i = 1, nb

            sb(i) = sb(i-1) + ds

         end do

      end subroutine get_uniform_grid



      !> \brief Calculates the distribution of points according
      !! to a geometric progression. The length of the first partition
      !! is defined on the left.
      subroutine get_gp_grid_left(nb, al, L, sb)
         implicit none
         integer, intent(in)  :: nb       !< Number of partitions over the new distribution of points
         real(8), intent(in)  :: al       !< Length of the first patition on the left
         real(8), intent(in)  :: L        !< Length of the curve
         real(8), intent(out) :: sb(0:nb) !< Length distribution over the curve

         ! Inner variables

         integer :: i
         real(8) :: q

         ! Calculates the ratio q of the geometric progression
         call get_gp_ratio(nb, L/al, q)

         sb(0) = 0.d0

         do i = 1, nb-1

            sb(i) = sb(i-1) + al * q ** (i-1)

         end do

         sb(nb) = L

      end subroutine get_gp_grid_left


      !> \brief Calculates the distribution of points according
      !! to a geometric progression. The length of the first partition
      !! is defined on the right.
      subroutine get_gp_grid_right(nb, ar, L, sb)
         implicit none
         integer, intent(in)  :: nb       !< Number of partitions over the new distribution of points
         real(8), intent(in)  :: ar       !< Length of the first patition on the right
         real(8), intent(in)  :: L        !< Length of the curve
         real(8), intent(out) :: sb(0:nb) !< Length distribution over the curve

         ! Inner variables

         integer :: i
         real(8) :: q

         ! Calculates the ratio q of the geometric progression
         call get_gp_ratio(nb, L/ar, q)

         sb(0) = 0.d0

         do i = 1, nb-1

            sb(i) = sb(i-1) + ar * q ** (nb-i)

         end do

         sb(nb) = L

      end subroutine get_gp_grid_right



      !> \brief Calculates the distribution of points according
      !! to a geometric progression. The length of the first partition
      !! is defined on the left and on the right.
      subroutine get_gp_grid_left_right(nb, al, ar, L, sb)
         implicit none
         integer, intent(in)  :: nb       !< Number of partitions over the new distribution of points
         real(8), intent(in)  :: al       !< Length of the first patition on the left
         real(8), intent(in)  :: ar       !< Length of the first patition on the right
         real(8), intent(in)  :: L        !< Length of the curve
         real(8), intent(out) :: sb(0:nb) !< Length distribution over the curve

         ! Inner parameters

         integer, parameter :: nitmax = 10000

         ! Inner variables

         integer :: i
         integer :: nl
         integer :: nr
         real(8) :: Lm
         real(8) :: Li
         real(8) :: Lf
         real(8) :: acl
         real(8) :: acr


         nl = nb / 2

         nr = nb - nl

         Li = 0.d0

         Lf = L

         do i = 1, nitmax

            Lm = 0.5d0 * ( Li + Lf )

            call get_gp_grid_left (nl, al, Lm, sb( 0:nl))

            call get_gp_grid_right(nr, ar, L-Lm, sb(nl:nb))

            sb(nl:nb) = sb(nl:nb) + Lm

            acl = sb(nl) - sb(nl-1)

            acr = sb(nl+1) - sb(nl)

            if ( acl > acr ) then

               Lf = Lm

            else

               Li = Lm

            end if

            if ( abs(acl-acr) < 1.d-14 ) exit

         end do

      end subroutine get_gp_grid_left_right


      !> \brief Calculates a distribution of points based on the arclength and
      !! on the double-exponential rule
      subroutine get_dexp_m1_ditribution( nbr, nb, aks, zcv, xbr, ybr, xb, yb) ! Output: last two
         implicit none
         integer, intent(in)  :: nbr        !< Number of partitions over the reference curve (xbr,ybr)
         integer, intent(in)  :: nb         !< Number of partitions for the new distribution of points
         real(8), intent(in)  :: aks        !< Exponent for the calculation of acv
         real(8), intent(in)  :: zcv        !< Central value of z
         real(8), intent(in)  :: xbr(0:nbr) !< Given distribution of points
         real(8), intent(in)  :: ybr(0:nbr) !< Given distribution of points
         real(8), intent(out) :: xb(0:nb)   !< New distribution of points
         real(8), intent(out) :: yb(0:nb)   !< New distribution of points

         ! Inner variables

         real(8) :: sbr(0:nbr) ! Length of each partition over the reference curve
         real(8) :: sb(0:nb)   ! Length of each partition over the new curve
         real(8) :: L          ! Total length
         real(8) :: acv        ! Length of the volume closer to the central point


         ! Calculates the length of the curve
         call get_sbr(nbr, xbr, ybr, sbr) ! Output: last one


         ! Total length
         L = sbr(nbr)


         ! Calculating the double-exponential distribution of arclength

         acv = L / ( nb ** aks )

         if ( acv > zcv ) then

            write(*,*) "Error. get_dexp_m1_ditribution: there is no resolution"&
               , " for this grid. Stopping."

            stop

         end if

         call get_dexp_m1(nb, 0.d0, zcv, L, acv, sb)


         ! Calculates the new distribution of points
         call get_new_distribution(nbr, nb, sbr, xbr, ybr, sb, xb, yb) ! Output: last two


      end subroutine get_dexp_m1_ditribution

      !> \brief Calculates a distribution of points based on the double-exponential
      !! rule. The coefficients of this rule are set by defining the position of
      !! an internal point (zcv) and the width (acv) of the volumes beside this
      !! point. The points are concentrated around zcv.
      subroutine get_dexp_m1(np, ziv, zcv, zfv, acv, z)
         implicit none
         integer, intent(in)  :: np  !< Number of partitions
         real(8), intent(in)  :: ziv !< Initial value of z
         real(8), intent(in)  :: zcv !< Central value of z (exponential's matching point)
         real(8), intent(in)  :: zfv !< Final value of z
         real(8), intent(in)  :: acv !< Length of the volumes beside zcv
         real(8), intent(out) :: z(0:np) !< Distribution of points

         ! Inner variables

         integer :: i
         integer :: npc
         real(8) :: t
         real(8) :: dt
         real(8) :: tc
         real(8) :: r
         real(8) :: A
         real(8) :: ac
         real(8) :: zc


         ! Initializing variables

         zc = zcv / ( zfv - ziv )

         ac = acv / ( zfv - ziv )

         dt = 1.d0 / dble(np)

         r = ac / ( 1.d0 - zc )



         ! Cheking conditions

         if ( ( ac / zc > 1.d0 ) .or. ( r > 1.d0 ) ) then

            write(*,*) "Error. get_dexp_m1: ac > zc or ac > (1-zc). Stopping."

            stop

         end if

         if ( ac > dt ) then

            write(*,*) "Error. get_dexp_m1: ac > dt. Stopping."

            stop

         end if


         ! Searching the value of A

         call adexp_m1(dt, r, zc, A)


         ! Calculating the fraction of volumes for each exponential

         tc = 1.d0 + log( 1.d0 + ( exp(-A*dt) - 1.d0 ) / r ) / A

         npc = nint(tc*np)



         ! Calculating the normallized distribution of points

         ! First exponential

         do i = 0, npc

            t = tc * dble(i) / dble(npc)

            z(i) = zc * ( exp(A*t) - 1.d0 ) / ( exp(A*tc) - 1.d0 )

         end do

         ! Second exponential

         do i = npc, np

            t = ( 1.d0-tc ) * dble(i-npc) / dble(np-npc) + tc

            z(i) = zc+(1.d0-zc)*(exp(-A*(t-tc))-1.d0)/(exp(-A*(1.d0-tc))-1.d0)

         end do

         z(0)  = 0.d0

         z(np) = 1.d0

         ! Calculating the extended distribution of points

         z = ( zfv - ziv ) * z + ziv

      end subroutine get_dexp_m1

      !> \brief Defines the transcendental function for finding A
      real(8) function fdexp_m1(dt, r, zc, A)
         implicit none
         real(8), intent(in) :: dt
         real(8), intent(in) :: r
         real(8), intent(in) :: zc
         real(8), intent(in) :: A

         ! Inner varialbles

         real(8) :: raux1, raux2, raux3

         if ( abs(A) > 1.d-6 ) then

            raux1 = exp(A) * ( 1.d0 + ( exp(-A*dt) - 1.d0 ) / r )
            raux2 = (2.d0*zc-1.d0) / zc
            raux3 = (1.d0-zc) / zc

            fdexp_m1 = -exp(-A) * raux1 * raux1 +  raux1 * raux2 + raux3

         else

            fdexp_m1 = (r*zc-r+dt)*A/(r*zc) &
               + ((r**2-2*dt**2)*zc-r**2+(2*dt-dt**2)*r)*A**2/(r**2*zc)/2.0

         end if

      end function fdexp_m1


      !> \brief Finds the value of A using the bissection method
      subroutine adexp_m1(dt, r, zc, A)
         implicit none
         real(8), intent(in)  :: dt
         real(8), intent(in)  :: r
         real(8), intent(in)  :: zc
         real(8), intent(out) :: A

         ! Parameters

         integer, parameter :: itmax = 1000

         ! Inner variables

         integer :: k
         real(8) :: Ai, Af, Am, Fi, Fm, Ff



         ! Initializing the search interval

         Ai = -1.d-18

         Af = -100.d0



         ! Checking for roots

         Fi = fdexp_m1(dt, r, zc, Ai)

         Ff = fdexp_m1(dt, r, zc, Af)

         if ( Fi * Ff > 0.d0 ) then

            write(*,*) "Error. adexp_m1: there is none or more than one roots. Stopping."

            stop

         end if


         ! Looking for a solution

         do k = 1, itmax

            Am = 0.5d0 * ( Ai + Af )

            Fi = fdexp_m1(dt, r, zc, Ai)

            Fm = fdexp_m1(dt, r, zc, Am)

            if ( Fi * Fm <= 0.d0 ) then

               Af = Am

            else

               Ai = Am

            end if

            if ( abs(Af-Ai) < 1.d-14 ) exit

         end do

         if ( k >= itmax ) then

            write(*,*) "Error. adexp_m1: number of iteractions exceeded. Stopping."

            stop

         end if

         A = 0.5d0 * ( Ai + Af )

      end subroutine adexp_m1


   end subroutine get_grid_boundary_g21




   !> \brief Calculates the coordinates (x,y) of the north and south boundaries.
   !! The north boundary line is generated from the south boundary.
   !! The south boundary line is an bluff nose followed by a smooth arc generated
   !! by cubic splines (second derivatives and boundary values given) of length lo
   !! and radius rb, followed by a cylinder of radius rb and length lc.
   !! The distribution of points over the south boundary is separated in three
   !! parts:
   !! Part   I: geometric progression (bluff nose+smooth arc)
   !! Part  II: geometric progression (smooth arc)
   !! Part III: geometric progression (cylinder)
   !! The distribution of points over the north boundary is  based on the normal
   !! vectors from the south boundary. The directions of the normal vectors
   !! may be changed by a smoothing procedure.
   !! Adapted to DEPP optimizer.
   subroutine get_grid_boundary_g22(nx, ny, unt, fgeom, lr, rb, x, y) ! Output: last two
      implicit none
      integer, intent(in) :: nx  !< Number of volumes in the csi direction (real+fictitious)
      integer, intent(in) :: ny  !< Number of volumes in the eta direction (real+fictitious)
      integer, intent(in) :: unt !< Unit where the input parameters will be printed
      character(len=*), intent(in) :: fgeom !< File of the geometric parameters

      real(8), intent(out) :: lr  !< Length of the body (m)
      real(8), intent(out) :: rb  !< Base radius of the body (m)
      real(8), dimension(nx*ny), intent(out) :: x !< Coord. x at the northeast corner of the volume P (m)
      real(8), dimension(nx*ny), intent(out) :: y !< Coord. y at the northeast corner of the volume P (m)

      ! Parameters

      integer, parameter :: nr = 10000 ! Number of partitions for the representing curve (xr,yr)


      ! Inner variables
      integer :: i, j, np ! Dummy indexes
      integer :: es  ! Exit status ( 0=success, 1=failure )
      integer :: mc  ! Monotonicity check ( 0=no, 1=yes )
      integer :: nsi ! Number of points for the spline interpolation

      integer :: kgb ! Kind of geometry of the body
      integer :: nx1 ! Number of points in the first part of the grid along the eta line
      integer :: nx2 ! Number of points in the first and second parts of the grid along the eta lines
      integer :: nx3 ! Number of points in the first to third parts of the grid along the eta lines
      integer :: iiv ! Initial value of i
      integer :: ifv ! Final value of i
      integer :: nr1 ! Fraction of nr used to represent the Part I

      real(8) :: a2  ! Width of the first volume of Part II
      real(8) :: a3  ! Width of the first volume of Part III

      real(8) :: xiv ! Initial value of x
      real(8) :: xfv ! Final value of x

      real(8) :: lo  ! Ogive length (m)
      real(8) :: w   ! Distance between south and north boundaries (m)
      real(8) :: fs  ! Factor of smoothing (0<=fs) (necessary for non-smooth surfaces)
      real(8) :: lc  ! Cylinder length (m)
      real(8) :: lo1 ! Length of the first part of the ogive (m)
      real(8) :: rf  ! Radius of the nose tip (m)
      real(8) :: aks ! Exponent of the power law for the calculation of the width of the first volume
      real(8) :: flo ! Fraction of lo
      real(8) :: fp1 ! Fraction of points in the first part of the grid
      real(8) :: fp2 ! Fraction of points in the second part of the grid


      real(8), dimension(nx) :: xbs ! Distribution of points over the south boundary
      real(8), dimension(nx) :: ybs ! Distribution of points over the south boundary

      real(8), dimension(0:nr) :: xr ! Distribution of points over the representing curve
      real(8), dimension(0:nr) :: yr ! Distribution of points over the representing curve

      real(8) :: y0  ! f(x0)
      real(8) :: yn  ! f(xn)

      real(8), allocatable :: xsi(:) ! x/lr Data for spline interpolation
      real(8), allocatable :: dy2(:) ! diff( f(x), x, 2)

      ! DEPP variables
      integer            :: nu     ! number of unknowns
      integer            :: ind    ! number of the individual
      integer,       allocatable :: xopt(:)   ! Optimization checker
      character(10), allocatable :: xname(:)  ! Name of parameters
      real(8),       allocatable :: xvalue(:) ! parameters
      character(200)             :: sname     ! simulation name


      ! Reading parameters

      open(13,file=fgeom)

      read(13,*) kgb  ! Kind of geometry of the body
      read(13,*) w    ! Distance between south and north boundaries (m)
      read(13,*) fs   ! Factor of smoothing (0<=fs) (necessary for non-smooth surfaces)
      read(13,*) lo   ! Length of the ogive (m)
      read(13,*) lc   ! Length of the cylinder (m)
      read(13,*) rb   ! Base radius of the body (m)
      read(13,*) aks  ! Exponent of the power law for the south boundary
      read(13,*) flo  ! Fraction of lo
      read(13,*) fp1  ! Fraction of volumes in the first part
      read(13,*) fp2  ! Fraction of volumes in the second part
      read(13,*) mc   ! Monocity check ( 0=no, 1=yes )
      read(13,*) nu   ! Number of unknowns

      if ( nu < 3 ) then

         write(*,*) "nu < 3. Stopping."

         stop

      end if

      nsi = nu - 1

      allocate ( xopt(nu), xname(nu), xvalue(nu) )

      allocate( xsi(nsi), dy2(nsi) )

      do i = 1, nsi

         read(13,*) xsi(i)

      end do

      do i = 1, nu

         read(13,*) xname(i) ! Parameter name
         read(13,*) xopt(i)  ! Will this parameter be optimized? ( 0 = no, 1 = yes )
         read(13,*) xvalue(i)! Parameter value

      end do

      close(13)


      ! Reads the parameters from DEPP
      call depp_get_parameters(nu, xopt, xname, xvalue, ind, sname)

      y0 = xvalue(1)

      yn = 1.d0

      dy2 = xvalue(2:nu)

      ! Gets a geometry using cubic splines with second derivatives prescribed
      ! on all points and the function values given on the boundaries
      call get_d2cspline_geometry(nr, nsi, 1.5d0, rb, lo, xsi, y0, yn, dy2, xr, yr, es) ! Output: last three


      ! Checking function monotonicity
      if ( mc == 1 ) then

         ! Non-monotonic function. Requesting new parameters.
         if ( es == 1 ) then

            call depp_save_fitness(-huge(1.d0), 2, "-cdfp")

            write(*,*) "Non-monotonic geometry. Stopping."

            stop

         end if

      end if

      ! Writing parameters

      write(unt,*)
      write(unt,*) " *** Geometry parameters ***"
      write(unt,*)

      write(unt,"(I23,    A)") kgb, " = kgb: Kind of geometry of the body"
      write(unt,"(ES23.16,A)")   w, " =   w: Distance between south and north boundaries (m)"
      write(unt,"(ES23.16,A)")  fs, " =  fs: Factor of smoothing (0<=fs) (necessary for non-smooth surfaces)"
      write(unt,"(ES23.16,A)")  lo, " =  lo: Length of the ogive (m)"
      write(unt,"(ES23.16,A)")  lc, " =  lc: Length of the cylinder (m)"
      write(unt,"(ES23.16,A)")  rb, " =  rb: Base radius of the body (m)"
      write(unt,"(ES23.16,A)") aks, " = aks: Exponent of the power law for the south boundary"
      write(unt,"(ES23.16,A)") flo, " = flo: Fraction of lo"
      write(unt,"(ES23.16,A)") fp1, " = fp1: Fraction of volumes in the first part of the grid"
      write(unt,"(ES23.16,A)") fp2, " = fp2: Fraction of volumes in the second part of the grid"
      write(unt,"(I23,    A)")  mc, " =  mc: Monotonicity check"
      write(unt,"(ES23.16,A)")  y0, " =  y0: f(x0)"

      do i = 1, nsi

         write(unt,"(ES23.16,A)") xsi(i), ' = xsi: x/lr for the spline calculation'

      end do

      do i = 1, nsi

         write(unt,"(ES23.16,A)") dy2(i), ' = diff(y,x,2)'

      end do


      ! Calculating the number of points in the first part of the grid

      nx1 = int( (nx-2) * fp1 )

      ! Calculating the number of points in the first and second parts of the grid

      nx2 = int( (nx-2) * fp2 ) + nx1

      ! Calculating the number of points in the 1st, 2nd and 3rd part of the grid

      nx3 = nx-1

      ! Calculating frontal radius

      rf = y0 * rb

      ! Calculating the fraction of lo used in the 2nd part of the grid

      lo1 = flo * lo

      do i = 1, nr

         if ( lo1 < xr(i) ) exit

      end do

      nr1 = i

      lo1 = xr(nr1)


      ! Calculating the body length

      lr = lo + lc


      ! ========================================================

      ! Generating south boundary

      j = 1

      ! Part I

      ! Calculates the distribution of points based on the arclength and the
      ! double exponential distribution

      iiv = 1   ! Initial value of i

      ifv = nx1 ! Final value of i

      call get_dexp_m1_ditribution( nr1, ifv-iiv, aks, rf, xr(0:nr1), yr(0:nr1)&
         , xbs(iiv:ifv), ybs(iiv:ifv)) ! Output: last two

      do i = iiv, ifv

         np   = nx * (j-1) + i

         x(np) = xbs(i)

         y(np) = ybs(i)

      end do


      ! Part II

      iiv = nx1  ! Initial value of i

      ifv = nx2 ! Final value of i

      xiv = lo1 ! Initial value of x

      xfv = lo  ! Final value of x

      ! Calculating the width of the first volume of Part II

      i = iiv

      a2 = sqrt( (xbs(i)-xbs(i-1))**2 + (ybs(i)-ybs(i-1))**2 )

      ! Calculating the distribution curve according to a geometric progression
      ! distribution based on the curve length

      call get_gp_length_distribution(nr-nr1, ifv-iiv, a2, -1.d0, xr(nr1:nr) &
         , yr(nr1:nr), xbs(iiv:ifv), ybs(iiv:ifv)) ! Output: last two

      do i = iiv, ifv

         np   = nx * (j-1) + i

         x(np) = xbs(i)

         y(np) = ybs(i)

      end do


      ! Part III

      iiv = nx2  ! Initial value of i

      ifv = nx3  ! Final value of i

      xiv = lo ! Initial value of x

      xfv = lo+lc  ! Final value of x

      ! Calculating the width of the first volume of Part IV

      i = iiv

      a3 = sqrt( (xbs(i)-xbs(i-1))**2 + (ybs(i)-ybs(i-1))**2 )

      call get_gp_distribution(iiv, ifv, xiv, xfv, a3, xbs(iiv:ifv)) ! Output: last one

      do i = iiv, ifv

         np   = nx * (j-1) + i

         ybs(i) = rb

         x(np) = xbs(i)

         y(np) = ybs(i)

      end do


      ! ========================================================

      ! Generating north boundary

      call get_north_boundary(nx, ny, w, fs, xbs, ybs, x, y) ! Output: last two

      ! ========================================================

   contains


      !> \brief Generates the north boundary of the domain based on the normal
      !! vectors from the south boundary. The directions of the normal vectors
      !! may be changed by a smoothing procedure.
      subroutine get_north_boundary(nx, ny, w, fs, xbs, ybs, x, y) ! Output: last two
         implicit none
         integer, intent(in)    :: nx  !< Number of volumes in the csi direction (real+fictitious)
         integer, intent(in)    :: ny  !< Number of volumes in the eta direction (real+fictitious)
         real(8), intent(in)    :: w   !< Width of the domain
         real(8), intent(in)    :: fs  !< Factor of smoothing (0<=fs)
         real(8), dimension(nx),    intent(in)  :: xbs !< Distribution of points over the south boundary
         real(8), dimension(nx),    intent(in)  :: ybs !< Distribution of points over the south boundary
         real(8), dimension(nx*ny), intent(out) :: x   !< Coord. x at the northeast corner of the volume P (m)
         real(8), dimension(nx*ny), intent(out) :: y   !< Coord. y at the northeast corner of the volume P (m)

         ! Inner variables

         integer :: i, j, k, np, ns
         real(8) :: raux
         real(8) :: tx
         real(8) :: ty
         real(8) :: tm
         real(8) :: vnx(nx-1)
         real(8) :: vny(nx-1)

         ! Number of iteractions for smoothing the direction vectors
         ns = int(nx * fs)

         ! Calculating the initial normal directions

         ! Direction of the west boundary
         vnx(1) = -1.d0
         vny(1) =  0.d0


         ! Direction of the east boundary
         vnx(nx-1) = 0.d0
         vny(nx-1) = 1.d0

         ! Directions normal to the body surface
         do i = 2, nx-2

            tx = ( xbs(i+1) - xbs(i-1) ) / 2.d0
            ty = ( ybs(i+1) - ybs(i-1) ) / 2.d0

            tm = sqrt( tx * tx + ty * ty )

            vnx(i) = - ty / tm
            vny(i) =   tx / tm

         end do

         ! Calculating the final directions (smoothing)

         do k = 1, ns

            do i = 2, nx-2

               vnx(i) = ( vnx(i+1) + vnx(i-1) ) / 2.d0

               vny(i) = ( vny(i+1) + vny(i-1) ) / 2.d0

               raux = dsqrt( vnx(i) * vnx(i) + vny(i) * vny(i) )

               vnx(i) = vnx(i) / raux

               vny(i) = vny(i) / raux

            end do

         end do

         ! Generating the north boundary

         j = ny-1

         do i = 1, nx-1

            np   = nx * (j-1) + i

            x(np) = xbs(i) + w * vnx(i)

            y(np) = ybs(i) + w * vny(i)

         end do

      end subroutine get_north_boundary


      !> \brief Calculates the node distribution based on a geometric progression distribution
      subroutine get_gp_distribution(iiv, ifv, xiv, xfv, a1, xb) ! Output: last one
         implicit none
         integer, intent(in)  :: iiv !< Initial value of i
         integer, intent(in)  :: ifv !< Final value of i
         real(8), intent(in)  :: xiv !< Initial value of x
         real(8), intent(in)  :: xfv !< Final value of x
         real(8), intent(in)  :: a1  !< Width of the first volume
         real(8), intent(out) :: xb(iiv:ifv) !< Distribution of points

         ! Inner variables

         integer :: i
         integer :: n
         real(8) :: r
         real(8) :: q


         n = ifv - iiv

         r = ( xfv - xiv ) / a1

         ! Calculating q
         call get_gp_ratio(n, r, q)

         xb(iiv) = xiv

         do i = iiv+1, ifv-1

            xb(i) = xb(i-1) + a1 * q ** ( i - iiv - 1 )

         end do

         xb(ifv) = xfv

      end subroutine get_gp_distribution


      !> \brief Calculates the base q of the geometric progression
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

           if ( 0.d0 < f(n, r, qi) * f(n, r, qm) ) then

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

      end subroutine get_gp_ratio

      !> \brief Function used in the subroutine 'get_gp_ratio'
      real(8) function f(n, r, q)
         implicit none
         integer, intent(in)  ::   n !< number of partitions
         real(8), intent(in)  ::   r !< l/a1
         real(8), intent(in)  ::   q !< q

         f = q ** n + r * ( 1.d0 - q ) - 1.d0

      end function f

      !> \brief Given a reference distribution of points (xbr,ybr)
      !! a new one (xb,yb) is created based on a geometric progression distributiuon
      !! based on the curve length.
      subroutine get_gp_length_distribution( nbr, nb, al, ar, xbr, ybr, xb, yb) ! Output: last two
         implicit none
         integer, intent(in)  :: nbr        !< Number of partitions over the reference curve (xbr,ybr)
         integer, intent(in)  :: nb         !< Number of partitions for the new distribution of points
         real(8), intent(in)  :: al         !< Length of the first partition on the left
         real(8), intent(in)  :: ar         !< Length of the first partition on the right
         real(8), intent(in)  :: xbr(0:nbr) !< Given distribution of points
         real(8), intent(in)  :: ybr(0:nbr) !< Given distribution of points
         real(8), intent(out) :: xb(0:nb)   !< New distribution of points
         real(8), intent(out) :: yb(0:nb)   !< New distribution of points

         ! Inner variables

         real(8) :: sbr(0:nbr) ! Length of each partition over the reference curve
         real(8) :: sb(0:nb)   ! Length of each partition over the new curve
         real(8) :: L          ! Total length


         ! Calculates the length of the curve
         call get_sbr(nbr, xbr, ybr, sbr) ! Output: last one


         ! Total length
         L = sbr(nbr)


         ! Selecting the kind of distribution

         if ( al < 0.d0 .and. ar < 0.d0 ) then

            ! Calculates the new points according to a uniform distribution
            call get_uniform_grid( nb, L, sb) ! Output: last one

         else if ( al > 0.d0 .and. ar < 0.d0 ) then

            ! Calculates the distribution of points according to a geometric
            ! progression. The length of the first partition is defined on the left.
            call get_gp_grid_left( nb, al, L, sb)

         else if ( al < 0.d0 .and. ar > 0.d0 ) then

            ! Calculates the distribution of points according to a geometric
            ! progression. The length of the first partition is defined on the right.
            call get_gp_grid_right( nb, ar, L, sb)

         else

            ! Calculates the distribution of points according to a geometric
            ! progression. The length of the first partition is defined on the left
            ! and on the right.
            call get_gp_grid_left_right(nb, al, ar, L, sb)

         end if


         ! Calculates the new distribution of points
         call get_new_distribution(nbr, nb, sbr, xbr, ybr, sb, xb, yb) ! Output: last two


      end subroutine get_gp_length_distribution

      !> \brief Calculates the length of the curve
      subroutine get_sbr(nbr, xbr, ybr, sbr) ! Output: last one
         implicit none
         integer, intent(in)  :: nbr        !< Number of partitions over the reference curve (xbr,ybr)
         real(8), intent(in)  :: xbr(0:nbr) !< Given distribution of points
         real(8), intent(in)  :: ybr(0:nbr) !< Given distribution of points
         real(8), intent(out) :: sbr(0:nbr) !< Length of the reference curve

         ! Inner variables

         integer :: i


         sbr(0) = 0.d0

         do i = 1, nbr

            sbr(i) = sbr(i-1) + sqrt( ( xbr(i) - xbr(i-1) )**2 + ( ybr(i) - ybr(i-1) )**2 )

         end do


      end subroutine get_sbr


      !> \brief Given the length distribution, the (xb,yb) is calculated
      subroutine get_new_distribution(nbr, nb, sbr, xbr, ybr, sb, xb, yb) ! Output: last two
         integer, intent(in)  :: nbr        !< Number of partitions over the reference curve (xbr,ybr)
         integer, intent(in)  :: nb         !< Number of partitions for the new distribution of points
         real(8), intent(in)  :: sbr(0:nbr) !< Length of the curve
         real(8), intent(in)  :: xbr(0:nbr) !< Given distribution of points
         real(8), intent(in)  :: ybr(0:nbr) !< Given distribution of points
         real(8), intent(in)  :: sb(0:nb)   !< New length distribution over the curve
         real(8), intent(out) :: xb(0:nb)   !< New distribution of points
         real(8), intent(out) :: yb(0:nb)   !< New distribution of points


         ! Inner variables

         integer :: i, j

         xb(0) = xbr(0)
         yb(0) = ybr(0)

         j = 0

         do i = 1, nb-1

            do j = j, nbr

               if( sbr(j) - sb(i) >= 0.d0 ) exit

            end do

            ! Linear interpolation

            xb(i) = xbr(j-1) + ( xbr(j)-xbr(j-1) ) * ( sb(i) - sbr(j-1) ) / ( sbr(j) - sbr(j-1) )
            yb(i) = ybr(j-1) + ( ybr(j)-ybr(j-1) ) * ( sb(i) - sbr(j-1) ) / ( sbr(j) - sbr(j-1) )

         end do

         xb(nb) = xbr(nbr)
         yb(nb) = ybr(nbr)

      end subroutine get_new_distribution



      !> \brief Calculates the uniform distribution of points
      subroutine get_uniform_grid( nb, L, sb) ! Output: last one
         implicit none
         integer, intent(in)  :: nb       !< Number of partitions over the new distribution of points
         real(8), intent(in)  :: L        !< Length of the curve
         real(8), intent(out) :: sb(0:nb) !< Length distribution over the curve

         ! Inner variables

         integer :: i

         real(8) :: ds

         ds = L / dble(nb)


         sb(0) = 0.d0

         do i = 1, nb

            sb(i) = sb(i-1) + ds

         end do

      end subroutine get_uniform_grid



      !> \brief Calculates the distribution of points according
      !! to a geometric progression. The length of the first partition
      !! is defined on the left.
      subroutine get_gp_grid_left(nb, al, L, sb)
         implicit none
         integer, intent(in)  :: nb       !< Number of partitions over the new distribution of points
         real(8), intent(in)  :: al       !< Length of the first patition on the left
         real(8), intent(in)  :: L        !< Length of the curve
         real(8), intent(out) :: sb(0:nb) !< Length distribution over the curve

         ! Inner variables

         integer :: i
         real(8) :: q

         ! Calculates the ratio q of the geometric progression
         call get_gp_ratio(nb, L/al, q)

         sb(0) = 0.d0

         do i = 1, nb-1

            sb(i) = sb(i-1) + al * q ** (i-1)

         end do

         sb(nb) = L

      end subroutine get_gp_grid_left


      !> \brief Calculates the distribution of points according
      !! to a geometric progression. The length of the first partition
      !! is defined on the right.
      subroutine get_gp_grid_right(nb, ar, L, sb)
         implicit none
         integer, intent(in)  :: nb       !< Number of partitions over the new distribution of points
         real(8), intent(in)  :: ar       !< Length of the first patition on the right
         real(8), intent(in)  :: L        !< Length of the curve
         real(8), intent(out) :: sb(0:nb) !< Length distribution over the curve

         ! Inner variables

         integer :: i
         real(8) :: q

         ! Calculates the ratio q of the geometric progression
         call get_gp_ratio(nb, L/ar, q)

         sb(0) = 0.d0

         do i = 1, nb-1

            sb(i) = sb(i-1) + ar * q ** (nb-i)

         end do

         sb(nb) = L

      end subroutine get_gp_grid_right



      !> \brief Calculates the distribution of points according
      !! to a geometric progression. The length of the first partition
      !! is defined on the left and on the right.
      subroutine get_gp_grid_left_right(nb, al, ar, L, sb)
         implicit none
         integer, intent(in)  :: nb       !< Number of partitions over the new distribution of points
         real(8), intent(in)  :: al       !< Length of the first patition on the left
         real(8), intent(in)  :: ar       !< Length of the first patition on the right
         real(8), intent(in)  :: L        !< Length of the curve
         real(8), intent(out) :: sb(0:nb) !< Length distribution over the curve

         ! Inner parameters

         integer, parameter :: nitmax = 10000

         ! Inner variables

         integer :: i
         integer :: nl
         integer :: nr
         real(8) :: Lm
         real(8) :: Li
         real(8) :: Lf
         real(8) :: acl
         real(8) :: acr


         nl = nb / 2

         nr = nb - nl

         Li = 0.d0

         Lf = L

         do i = 1, nitmax

            Lm = 0.5d0 * ( Li + Lf )

            call get_gp_grid_left (nl, al, Lm, sb( 0:nl))

            call get_gp_grid_right(nr, ar, L-Lm, sb(nl:nb))

            sb(nl:nb) = sb(nl:nb) + Lm

            acl = sb(nl) - sb(nl-1)

            acr = sb(nl+1) - sb(nl)

            if ( acl > acr ) then

               Lf = Lm

            else

               Li = Lm

            end if

            if ( abs(acl-acr) < 1.d-14 ) exit

         end do

      end subroutine get_gp_grid_left_right


      !> \brief Calculates a distribution of points based on the arclength and
      !! on the double-exponential rule
      subroutine get_dexp_m1_ditribution( nbr, nb, aks, zcv, xbr, ybr, xb, yb) ! Output: last two
         implicit none
         integer, intent(in)  :: nbr        !< Number of partitions over the reference curve (xbr,ybr)
         integer, intent(in)  :: nb         !< Number of partitions for the new distribution of points
         real(8), intent(in)  :: aks        !< Exponent for the calculation of acv
         real(8), intent(in)  :: zcv        !< Central value of z
         real(8), intent(in)  :: xbr(0:nbr) !< Given distribution of points
         real(8), intent(in)  :: ybr(0:nbr) !< Given distribution of points
         real(8), intent(out) :: xb(0:nb)   !< New distribution of points
         real(8), intent(out) :: yb(0:nb)   !< New distribution of points

         ! Inner variables

         real(8) :: sbr(0:nbr) ! Length of each partition over the reference curve
         real(8) :: sb(0:nb)   ! Length of each partition over the new curve
         real(8) :: L          ! Total length
         real(8) :: acv        ! Length of the volume closer to the central point


         ! Calculates the length of the curve
         call get_sbr(nbr, xbr, ybr, sbr) ! Output: last one


         ! Total length
         L = sbr(nbr)


         ! Calculating the double-exponential distribution of arclength

         acv = L / ( nb ** aks )

         if ( acv > zcv ) then

            write(*,*) "Error. get_dexp_m1_ditribution: there is no resolution"&
               , " for this grid. Stopping."

            stop

         end if

         call get_dexp_m1(nb, 0.d0, zcv, L, acv, sb)


         ! Calculates the new distribution of points
         call get_new_distribution(nbr, nb, sbr, xbr, ybr, sb, xb, yb) ! Output: last two


      end subroutine get_dexp_m1_ditribution

      !> \brief Calculates a distribution of points based on the double-exponential
      !! rule. The coefficients of this rule are set by defining the position of
      !! an internal point (zcv) and the width (acv) of the volumes beside this
      !! point. The points are concentrated around zcv.
      subroutine get_dexp_m1(np, ziv, zcv, zfv, acv, z)
         implicit none
         integer, intent(in)  :: np  !< Number of partitions
         real(8), intent(in)  :: ziv !< Initial value of z
         real(8), intent(in)  :: zcv !< Central value of z (exponential's matching point)
         real(8), intent(in)  :: zfv !< Final value of z
         real(8), intent(in)  :: acv !< Length of the volumes beside zcv
         real(8), intent(out) :: z(0:np) !< Distribution of points

         ! Inner variables

         integer :: i
         integer :: npc
         real(8) :: t
         real(8) :: dt
         real(8) :: tc
         real(8) :: r
         real(8) :: A
         real(8) :: ac
         real(8) :: zc


         ! Initializing variables

         zc = zcv / ( zfv - ziv )

         ac = acv / ( zfv - ziv )

         dt = 1.d0 / dble(np)

         r = ac / ( 1.d0 - zc )



         ! Cheking conditions

         if ( ( ac / zc > 1.d0 ) .or. ( r > 1.d0 ) ) then

            write(*,*) "Error. get_dexp_m1: ac > zc or ac > (1-zc). Stopping."

            stop

         end if

         if ( ac > dt ) then

            write(*,*) "Error. get_dexp_m1: ac > dt. Stopping."

            stop

         end if


         ! Searching the value of A

         call adexp_m1(dt, r, zc, A)


         ! Calculating the fraction of volumes for each exponential

         tc = 1.d0 + log( 1.d0 + ( exp(-A*dt) - 1.d0 ) / r ) / A

         npc = nint(tc*np)



         ! Calculating the normallized distribution of points

         ! First exponential

         do i = 0, npc

            t = tc * dble(i) / dble(npc)

            z(i) = zc * ( exp(A*t) - 1.d0 ) / ( exp(A*tc) - 1.d0 )

         end do

         ! Second exponential

         do i = npc, np

            t = ( 1.d0-tc ) * dble(i-npc) / dble(np-npc) + tc

            z(i) = zc+(1.d0-zc)*(exp(-A*(t-tc))-1.d0)/(exp(-A*(1.d0-tc))-1.d0)

         end do

         z(0)  = 0.d0

         z(np) = 1.d0

         ! Calculating the extended distribution of points

         z = ( zfv - ziv ) * z + ziv

      end subroutine get_dexp_m1

      !> \brief Defines the transcendental function for finding A
      real(8) function fdexp_m1(dt, r, zc, A)
         implicit none
         real(8), intent(in) :: dt
         real(8), intent(in) :: r
         real(8), intent(in) :: zc
         real(8), intent(in) :: A

         ! Inner varialbles

         real(8) :: raux1, raux2, raux3

         if ( abs(A) > 1.d-6 ) then

            raux1 = exp(A) * ( 1.d0 + ( exp(-A*dt) - 1.d0 ) / r )
            raux2 = (2.d0*zc-1.d0) / zc
            raux3 = (1.d0-zc) / zc

            fdexp_m1 = -exp(-A) * raux1 * raux1 +  raux1 * raux2 + raux3

         else

            fdexp_m1 = (r*zc-r+dt)*A/(r*zc) &
               + ((r**2-2*dt**2)*zc-r**2+(2*dt-dt**2)*r)*A**2/(r**2*zc)/2.0

         end if

      end function fdexp_m1


      !> \brief Finds the value of A using the bissection method
      subroutine adexp_m1(dt, r, zc, A)
         implicit none
         real(8), intent(in)  :: dt
         real(8), intent(in)  :: r
         real(8), intent(in)  :: zc
         real(8), intent(out) :: A

         ! Parameters

         integer, parameter :: itmax = 1000

         ! Inner variables

         integer :: k
         real(8) :: Ai, Af, Am, Fi, Fm, Ff



         ! Initializing the search interval

         Ai = -1.d-18

         Af = -100.d0



         ! Checking for roots

         Fi = fdexp_m1(dt, r, zc, Ai)

         Ff = fdexp_m1(dt, r, zc, Af)

         if ( Fi * Ff > 0.d0 ) then

            write(*,*) "Error. adexp_m1: there is none or more than one roots. Stopping."

            stop

         end if


         ! Looking for a solution

         do k = 1, itmax

            Am = 0.5d0 * ( Ai + Af )

            Fi = fdexp_m1(dt, r, zc, Ai)

            Fm = fdexp_m1(dt, r, zc, Am)

            if ( Fi * Fm <= 0.d0 ) then

               Af = Am

            else

               Ai = Am

            end if

            if ( abs(Af-Ai) < 1.d-14 ) exit

         end do

         if ( k >= itmax ) then

            write(*,*) "Error. adexp_m1: number of iteractions exceeded. Stopping."

            stop

         end if

         A = 0.5d0 * ( Ai + Af )

      end subroutine adexp_m1


   end subroutine get_grid_boundary_g22


   !> \brief Calculates the coordinates (x,y) of the north and south boundaries.
   !! The north boundary line is generated from the south boundary.
   !! The south boundary line is an bluff nose followed by a smooth arc generated
   !! by quartic splines of length lo and radius rb, followed by a cylinder of radius rb and length lc.
   !! The splines are define through the prescribed first derivatives and other boundary
   !! conditions described in the geometry_spline_s01 module.
   !! The distribution of points over the south boundary is separated in three
   !! parts:
   !! Part   I: geometric progression (bluff nose+smooth arc)
   !! Part  II: geometric progression (smooth arc)
   !! Part III: geometric progression (cylinder)
   !! The distribution of points over the north boundary is  based on the normal
   !! vectors from the south boundary. The directions of the normal vectors
   !! may be changed by a smoothing procedure.
   !! Adapted to DEPP optimizer.
   subroutine get_grid_boundary_g23(nx, ny, unt, fgeom, lr, rb, x, y) ! Output: last two
      implicit none
      integer, intent(in) :: nx  !< Number of volumes in the csi direction (real+fictitious)
      integer, intent(in) :: ny  !< Number of volumes in the eta direction (real+fictitious)
      integer, intent(in) :: unt !< Unit where the input parameters will be printed
      character(len=*), intent(in) :: fgeom !< File of the geometric parameters

      real(8), intent(out) :: lr  !< Length of the body (m)
      real(8), intent(out) :: rb  !< Base radius of the body (m)
      real(8), dimension(nx*ny), intent(out) :: x !< Coord. x at the northeast corner of the volume P (m)
      real(8), dimension(nx*ny), intent(out) :: y !< Coord. y at the northeast corner of the volume P (m)

      ! Parameters

      integer, parameter :: nr = 10000 ! Number of partitions for the representing curve (xr,yr)


      ! Inner variables
      integer :: i, j, np ! Dummy indexes
      integer :: es  ! Exit status ( 0=success, 1=failure )
      integer :: mc  ! Monotonicity check ( 0=no, 1=yes )
      integer :: nsi ! Number of points for the spline interpolation

      integer :: kgb ! Kind of geometry of the body
      integer :: nx1 ! Number of points in the first part of the grid along the eta line
      integer :: nx2 ! Number of points in the first and second parts of the grid along the eta lines
      integer :: nx3 ! Number of points in the first to third parts of the grid along the eta lines
      integer :: iiv ! Initial value of i
      integer :: ifv ! Final value of i
      integer :: nr1 ! Fraction of nr used to represent the Part I
      integer :: opt ! Options for boundary conditions
      integer :: onb ! Options for north boundary

      real(8) :: a2  ! Width of the first volume of Part II
      real(8) :: a3  ! Width of the first volume of Part III

      real(8) :: xiv ! Initial value of x
      real(8) :: xfv ! Final value of x

      real(8) :: lo  ! Ogive length (m)
      real(8) :: wf  ! Frontal distance between south and north boundaries (m)
      real(8) :: wb  ! Base distance between south and north boundaries (m)
      real(8) :: fs  ! Factor of smoothing (0<=fs) (necessary for non-smooth surfaces)
      real(8) :: lc  ! Cylinder length (m)
      real(8) :: lo1 ! Length of the first part of the ogive (m)
      real(8) :: rf  ! Radius of the nose tip (m)
      real(8) :: aks ! Exponent of the power law for the calculation of the width of the first volume
      real(8) :: flo ! Fraction of lo
      real(8) :: fp1 ! Fraction of points in the first part of the grid
      real(8) :: fp2 ! Fraction of points in the second part of the grid


      real(8), dimension(nx) :: xbs ! Distribution of points over the south boundary
      real(8), dimension(nx) :: ybs ! Distribution of points over the south boundary

      real(8), dimension(0:nr) :: xr ! Distribution of points over the representing curve
      real(8), dimension(0:nr) :: yr ! Distribution of points over the representing curve

      real(8) :: dy2  ! diff( f(x), x=xn, 2)
      real(8) :: dy3  ! diff( f(x), x=xn, 3)

      real(8), allocatable :: xsi(:) ! x/lr Data for spline interpolation
      real(8), allocatable :: dy1(:) ! diff( f(x), x)
      real(8), allocatable :: tht(:) ! atan(dy1)*180/pi

      ! DEPP variables
      integer            :: nu     ! number of unknowns
      integer            :: ind    ! number of the individual
      integer,       allocatable :: xopt(:)   ! Optimization checker
      character(10), allocatable :: xname(:)  ! Name of parameters
      real(8),       allocatable :: xvalue(:) ! parameters
      character(200)             :: sname     ! simulation name


      ! Reading parameters

      open(13,file=fgeom)

      read(13,*) kgb  ! Kind of geometry of the body
      read(13,*) onb  ! Option for north boundary
      read(13,*) wf   ! Frontal distance between south and north boundaries (m)
      read(13,*) wb   ! Base distance between south and north boundaries (m)
      read(13,*) fs   ! Factor of smoothing (0<=fs) (necessary for non-smooth surfaces)
      read(13,*) lo   ! Length of the ogive (m)
      read(13,*) lc   ! Length of the cylinder (m)
      read(13,*) rb   ! Base radius of the body (m)
      read(13,*) aks  ! Exponent of the power law for the south boundary
      read(13,*) flo  ! Fraction of lo
      read(13,*) fp1  ! Fraction of volumes in the first part
      read(13,*) fp2  ! Fraction of volumes in the second part
      read(13,*) mc   ! Monocity check ( 0=no, 1=yes )
      read(13,*) opt  ! Boundary condition option
      read(13,*) nu   ! Number of unknowns

      if ( nu < 4 ) then

         write(*,*) "nu < 4. Stopping."

         stop

      end if

      nsi = nu - 3

      allocate ( xopt(nu), xname(nu), xvalue(nu) )

      allocate( xsi(0:nsi), dy1(0:nsi), tht(0:nsi) )

      do i = 0, nsi

         read(13,*) xsi(i)

      end do

      do i = 1, nu

         read(13,*) xname(i) ! Parameter name
         read(13,*) xopt(i)  ! Will this parameter be optimized? ( 0 = no, 1 = yes )
         read(13,*) xvalue(i)! Parameter value

      end do

      close(13)


      ! Reads the parameters from DEPP
      call depp_get_parameters(nu, xopt, xname, xvalue, ind, sname)

      dy2 = xvalue(1)

      dy3 = xvalue(2)

      tht = xvalue(3:nu)

      ! Calculating the slope

      dy1 = tan(tht*acos(-1.d0)/180.d0)

      ! Gets the geometry spline s01
      call get_spline_s01_geometry(nr, nsi, 1.5d0, rb, lo, xsi, dy1, dy2, dy3, opt, xr, yr, es) ! Output: last three


      ! Checking exit status

      ! Negative function. Requesting new parameters.
      if ( es > 1 ) then

         call depp_save_fitness(-huge(1.d0), 2, "-cdfp")

         write(*,*) "Negative geometry. Stopping."

         stop

      else if ( es == 1 .and. mc == 1 ) then

         ! Non-monotonic function. Requesting new parameters.

         call depp_save_fitness(-huge(1.d0), 2, "-cdfp")

         write(*,*) "Non-monotonic geometry. Stopping."

         stop

      end if


      ! Writing parameters

      write(unt,*)
      write(unt,*) " *** Geometry parameters ***"
      write(unt,*)

      write(unt,"(I23,    A)") kgb, " = kgb: Kind of geometry of the body"
      write(unt,"(I23,    A)") onb, " = onb: Option for north boundary"
      write(unt,"(ES23.16,A)")  wf, " =  wf: Frontal distance between south and north boundaries (m)"
      write(unt,"(ES23.16,A)")  wb, " =  wb: Base distance between south and north boundaries (m)"
      write(unt,"(ES23.16,A)")  fs, " =  fs: Factor of smoothing (0<=fs) (necessary for non-smooth surfaces)"
      write(unt,"(ES23.16,A)")  lo, " =  lo: Length of the ogive (m)"
      write(unt,"(ES23.16,A)")  lc, " =  lc: Length of the cylinder (m)"
      write(unt,"(ES23.16,A)")  rb, " =  rb: Base radius of the body (m)"
      write(unt,"(ES23.16,A)") aks, " = aks: Exponent of the power law for the south boundary"
      write(unt,"(ES23.16,A)") flo, " = flo: Fraction of lo"
      write(unt,"(ES23.16,A)") fp1, " = fp1: Fraction of volumes in the first part of the grid"
      write(unt,"(ES23.16,A)") fp2, " = fp2: Fraction of volumes in the second part of the grid"
      write(unt,"(I23,    A)")  mc, " =  mc: Monotonicity check"
      write(unt,"(I23,    A)") opt, " = opt: Boundary condition option"
      write(unt,"(ES23.16,A)") dy2, " = dy2: f''(xn)"
      write(unt,"(ES23.16,A)") dy3, " = dy3: f'''(xn)"

      do i = 0, nsi

         write(unt,"(ES23.16,A)") xsi(i), ' = xsi: x/lr for the spline calculation'

      end do

      do i = 0, nsi

         write(unt,"(ES23.16,A)") tht(i), ' = tht(x)'

      end do


      ! Calculating the number of points in the first part of the grid

      nx1 = int( (nx-2) * fp1 )

      ! Calculating the number of points in the first and second parts of the grid

      nx2 = int( (nx-2) * fp2 ) + nx1

      ! Calculating the number of points in the 1st, 2nd and 3rd part of the grid

      nx3 = nx-1

      ! Calculating frontal radius

      rf = yr(1)

      ! Calculating the fraction of lo used in the 2nd part of the grid

      lo1 = flo * lo

      do i = 1, nr

         if ( lo1 < xr(i) ) exit

      end do

      nr1 = i

      lo1 = xr(nr1)


      ! Calculating the body length

      lr = lo + lc


      ! ========================================================

      ! Generating south boundary

      j = 1

      ! Part I

      ! Calculates the distribution of points based on the arclength and the
      ! double exponential distribution

      iiv = 1   ! Initial value of i

      ifv = nx1 ! Final value of i

      call get_dexp_m1_ditribution( nr1, ifv-iiv, aks, rf, xr(0:nr1), yr(0:nr1)&
         , xbs(iiv:ifv), ybs(iiv:ifv)) ! Output: last two

      do i = iiv, ifv

         np   = nx * (j-1) + i

         x(np) = xbs(i)

         y(np) = ybs(i)

      end do


      ! Part II

      iiv = nx1  ! Initial value of i

      ifv = nx2 ! Final value of i

      xiv = lo1 ! Initial value of x

      xfv = lo  ! Final value of x

      ! Calculating the width of the first volume of Part II

      i = iiv

      a2 = sqrt( (xbs(i)-xbs(i-1))**2 + (ybs(i)-ybs(i-1))**2 )

      ! Calculating the distribution curve according to a geometric progression
      ! distribution based on the curve length

      call get_gp_length_distribution(nr-nr1, ifv-iiv, a2, -1.d0, xr(nr1:nr) &
         , yr(nr1:nr), xbs(iiv:ifv), ybs(iiv:ifv)) ! Output: last two

      do i = iiv, ifv

         np   = nx * (j-1) + i

         x(np) = xbs(i)

         y(np) = ybs(i)

      end do


      ! Part III

      iiv = nx2  ! Initial value of i

      ifv = nx3  ! Final value of i

      xiv = lo ! Initial value of x

      xfv = lo+lc  ! Final value of x

      ! Calculating the width of the first volume of Part IV

      i = iiv

      a3 = sqrt( (xbs(i)-xbs(i-1))**2 + (ybs(i)-ybs(i-1))**2 )

      call get_gp_distribution(iiv, ifv, xiv, xfv, a3, xbs(iiv:ifv)) ! Output: last one

      do i = iiv, ifv

         np   = nx * (j-1) + i

         ybs(i) = rb

         x(np) = xbs(i)

         y(np) = ybs(i)

      end do


      ! ========================================================

      ! Generating north boundary

      call get_north_boundary(nx, ny, onb, wf, wb, fs, xbs, ybs, x, y) ! Output: last two

      ! ========================================================

   contains


      !> \brief Generates the north boundary of the domain based on the normal
      !! vectors from the south boundary. The directions of the normal vectors
      !! may be changed by a smoothing procedure.
      subroutine get_north_boundary(nx, ny, onb, wf, wb, fs, xbs, ybs, x, y) ! Output: last two
         implicit none
         integer, intent(in)    :: nx  !< Number of volumes in the csi direction (real+fictitious)
         integer, intent(in)    :: ny  !< Number of volumes in the eta direction (real+fictitious)
         integer, intent(in)    :: onb !< Option for north boundary (1=equidistant,2=parabolic)
         real(8), intent(in)    :: wf  !< Frontal width
         real(8), intent(in)    :: wb  !< Base width
         real(8), intent(in)    :: fs  !< Factor of smoothing (0<=fs)
         real(8), dimension(nx),    intent(in)  :: xbs !< Distribution of points over the south boundary
         real(8), dimension(nx),    intent(in)  :: ybs !< Distribution of points over the south boundary
         real(8), dimension(nx*ny), intent(out) :: x   !< Coord. x at the northeast corner of the volume P (m)
         real(8), dimension(nx*ny), intent(out) :: y   !< Coord. y at the northeast corner of the volume P (m)

         ! Inner variables

         integer :: i, j, k, np, ns
         real(8) :: raux
         real(8) :: b, c, la, lb
         real(8) :: tx
         real(8) :: ty
         real(8) :: tm
         real(8) :: vnx(nx-1)
         real(8) :: vny(nx-1)

         ! Number of iteractions for smoothing the direction vectors
         ns = int(nx * fs)

         ! Calculating the initial normal directions

         ! Direction of the west boundary
         vnx(1) = -1.d0
         vny(1) =  0.d0


         ! Direction of the east boundary
         vnx(nx-1) = 0.d0
         vny(nx-1) = 1.d0

         ! Directions normal to the body surface
         do i = 2, nx-2

            tx = ( xbs(i+1) - xbs(i-1) ) / 2.d0
            ty = ( ybs(i+1) - ybs(i-1) ) / 2.d0

            tm = sqrt( tx * tx + ty * ty )

            vnx(i) = - ty / tm
            vny(i) =   tx / tm

         end do

         ! Calculating the final directions (smoothing)

         do k = 1, ns

            do i = 2, nx-2

               vnx(i) = ( vnx(i+1) + vnx(i-1) ) / 2.d0

               vny(i) = ( vny(i+1) + vny(i-1) ) / 2.d0

               raux = dsqrt( vnx(i) * vnx(i) + vny(i) * vny(i) )

               vnx(i) = vnx(i) / raux

               vny(i) = vny(i) / raux

            end do

         end do

         ! Generating the north boundary

         select case (onb)

            case (1) ! Equidistant

               j = ny-1

               do i = 1, nx-1

                  np   = nx * (j-1) + i

                  x(np) = xbs(i) + wb * vnx(i)

                  y(np) = ybs(i) + wb * vny(i)

               end do

            case (2) ! Parabolic

               la = xbs(nx-1) + wf

               lb = ybs(nx-1) + wb

               j = ny-1

               do i = 1, nx-1

                  np   = nx * (j-1) + i

                  if ( 1.d-3 * abs(vnx(i)) < abs(vny(i)) ) then

                     b = - vnx(i) * lb / ( vny(i) * la )

                     c = ( ybs(i) * vnx(i) / vny(i) - (wf+xbs(i)) ) / la

                     y(np) = lb * 0.5d0 * (-b+sqrt(b**2-4.d0*c))

                     x(np) = la * ( y(np) / lb ) ** 2.d0 - wf

                  else

                     y(np) = ybs(i) - vny(i) / vnx(i) * (wf+xbs(i))

                     x(np) = la * ( y(np) / lb ) ** 2.d0 - wf

                  end if

               end do

            case default

               write(*,*) 'Unknown boundary. Stopping.'

               stop

         end select

      end subroutine get_north_boundary


      !> \brief Calculates the node distribution based on a geometric progression distribution
      subroutine get_gp_distribution(iiv, ifv, xiv, xfv, a1, xb) ! Output: last one
         implicit none
         integer, intent(in)  :: iiv !< Initial value of i
         integer, intent(in)  :: ifv !< Final value of i
         real(8), intent(in)  :: xiv !< Initial value of x
         real(8), intent(in)  :: xfv !< Final value of x
         real(8), intent(in)  :: a1  !< Width of the first volume
         real(8), intent(out) :: xb(iiv:ifv) !< Distribution of points

         ! Inner variables

         integer :: i
         integer :: n
         real(8) :: r
         real(8) :: q


         n = ifv - iiv

         r = ( xfv - xiv ) / a1

         ! Calculating q
         call get_gp_ratio(n, r, q)

         xb(iiv) = xiv

         do i = iiv+1, ifv-1

            xb(i) = xb(i-1) + a1 * q ** ( i - iiv - 1 )

         end do

         xb(ifv) = xfv

      end subroutine get_gp_distribution


      !> \brief Calculates the base q of the geometric progression
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

           if ( 0.d0 < f(n, r, qi) * f(n, r, qm) ) then

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

      end subroutine get_gp_ratio

      !> \brief Function used in the subroutine 'get_gp_ratio'
      real(8) function f(n, r, q)
         implicit none
         integer, intent(in)  ::   n !< number of partitions
         real(8), intent(in)  ::   r !< l/a1
         real(8), intent(in)  ::   q !< q

         f = q ** n + r * ( 1.d0 - q ) - 1.d0

      end function f

      !> \brief Given a reference distribution of points (xbr,ybr)
      !! a new one (xb,yb) is created based on a geometric progression distributiuon
      !! based on the curve length.
      subroutine get_gp_length_distribution( nbr, nb, al, ar, xbr, ybr, xb, yb) ! Output: last two
         implicit none
         integer, intent(in)  :: nbr        !< Number of partitions over the reference curve (xbr,ybr)
         integer, intent(in)  :: nb         !< Number of partitions for the new distribution of points
         real(8), intent(in)  :: al         !< Length of the first partition on the left
         real(8), intent(in)  :: ar         !< Length of the first partition on the right
         real(8), intent(in)  :: xbr(0:nbr) !< Given distribution of points
         real(8), intent(in)  :: ybr(0:nbr) !< Given distribution of points
         real(8), intent(out) :: xb(0:nb)   !< New distribution of points
         real(8), intent(out) :: yb(0:nb)   !< New distribution of points

         ! Inner variables

         real(8) :: sbr(0:nbr) ! Length of each partition over the reference curve
         real(8) :: sb(0:nb)   ! Length of each partition over the new curve
         real(8) :: L          ! Total length


         ! Calculates the length of the curve
         call get_sbr(nbr, xbr, ybr, sbr) ! Output: last one


         ! Total length
         L = sbr(nbr)


         ! Selecting the kind of distribution

         if ( al < 0.d0 .and. ar < 0.d0 ) then

            ! Calculates the new points according to a uniform distribution
            call get_uniform_grid( nb, L, sb) ! Output: last one

         else if ( al > 0.d0 .and. ar < 0.d0 ) then

            ! Calculates the distribution of points according to a geometric
            ! progression. The length of the first partition is defined on the left.
            call get_gp_grid_left( nb, al, L, sb)

         else if ( al < 0.d0 .and. ar > 0.d0 ) then

            ! Calculates the distribution of points according to a geometric
            ! progression. The length of the first partition is defined on the right.
            call get_gp_grid_right( nb, ar, L, sb)

         else

            ! Calculates the distribution of points according to a geometric
            ! progression. The length of the first partition is defined on the left
            ! and on the right.
            call get_gp_grid_left_right(nb, al, ar, L, sb)

         end if


         ! Calculates the new distribution of points
         call get_new_distribution(nbr, nb, sbr, xbr, ybr, sb, xb, yb) ! Output: last two


      end subroutine get_gp_length_distribution

      !> \brief Calculates the length of the curve
      subroutine get_sbr(nbr, xbr, ybr, sbr) ! Output: last one
         implicit none
         integer, intent(in)  :: nbr        !< Number of partitions over the reference curve (xbr,ybr)
         real(8), intent(in)  :: xbr(0:nbr) !< Given distribution of points
         real(8), intent(in)  :: ybr(0:nbr) !< Given distribution of points
         real(8), intent(out) :: sbr(0:nbr) !< Length of the reference curve

         ! Inner variables

         integer :: i


         sbr(0) = 0.d0

         do i = 1, nbr

            sbr(i) = sbr(i-1) + sqrt( ( xbr(i) - xbr(i-1) )**2 + ( ybr(i) - ybr(i-1) )**2 )

         end do


      end subroutine get_sbr


      !> \brief Given the length distribution, the (xb,yb) is calculated
      subroutine get_new_distribution(nbr, nb, sbr, xbr, ybr, sb, xb, yb) ! Output: last two
         integer, intent(in)  :: nbr        !< Number of partitions over the reference curve (xbr,ybr)
         integer, intent(in)  :: nb         !< Number of partitions for the new distribution of points
         real(8), intent(in)  :: sbr(0:nbr) !< Length of the curve
         real(8), intent(in)  :: xbr(0:nbr) !< Given distribution of points
         real(8), intent(in)  :: ybr(0:nbr) !< Given distribution of points
         real(8), intent(in)  :: sb(0:nb)   !< New length distribution over the curve
         real(8), intent(out) :: xb(0:nb)   !< New distribution of points
         real(8), intent(out) :: yb(0:nb)   !< New distribution of points


         ! Inner variables

         integer :: i, j

         xb(0) = xbr(0)
         yb(0) = ybr(0)

         j = 0

         do i = 1, nb-1

            do j = j, nbr

               if( sbr(j) - sb(i) >= 0.d0 ) exit

            end do

            ! Linear interpolation

            xb(i) = xbr(j-1) + ( xbr(j)-xbr(j-1) ) * ( sb(i) - sbr(j-1) ) / ( sbr(j) - sbr(j-1) )
            yb(i) = ybr(j-1) + ( ybr(j)-ybr(j-1) ) * ( sb(i) - sbr(j-1) ) / ( sbr(j) - sbr(j-1) )

         end do

         xb(nb) = xbr(nbr)
         yb(nb) = ybr(nbr)

      end subroutine get_new_distribution



      !> \brief Calculates the uniform distribution of points
      subroutine get_uniform_grid( nb, L, sb) ! Output: last one
         implicit none
         integer, intent(in)  :: nb       !< Number of partitions over the new distribution of points
         real(8), intent(in)  :: L        !< Length of the curve
         real(8), intent(out) :: sb(0:nb) !< Length distribution over the curve

         ! Inner variables

         integer :: i

         real(8) :: ds

         ds = L / dble(nb)


         sb(0) = 0.d0

         do i = 1, nb

            sb(i) = sb(i-1) + ds

         end do

      end subroutine get_uniform_grid



      !> \brief Calculates the distribution of points according
      !! to a geometric progression. The length of the first partition
      !! is defined on the left.
      subroutine get_gp_grid_left(nb, al, L, sb)
         implicit none
         integer, intent(in)  :: nb       !< Number of partitions over the new distribution of points
         real(8), intent(in)  :: al       !< Length of the first patition on the left
         real(8), intent(in)  :: L        !< Length of the curve
         real(8), intent(out) :: sb(0:nb) !< Length distribution over the curve

         ! Inner variables

         integer :: i
         real(8) :: q

         ! Calculates the ratio q of the geometric progression
         call get_gp_ratio(nb, L/al, q)

         sb(0) = 0.d0

         do i = 1, nb-1

            sb(i) = sb(i-1) + al * q ** (i-1)

         end do

         sb(nb) = L

      end subroutine get_gp_grid_left


      !> \brief Calculates the distribution of points according
      !! to a geometric progression. The length of the first partition
      !! is defined on the right.
      subroutine get_gp_grid_right(nb, ar, L, sb)
         implicit none
         integer, intent(in)  :: nb       !< Number of partitions over the new distribution of points
         real(8), intent(in)  :: ar       !< Length of the first patition on the right
         real(8), intent(in)  :: L        !< Length of the curve
         real(8), intent(out) :: sb(0:nb) !< Length distribution over the curve

         ! Inner variables

         integer :: i
         real(8) :: q

         ! Calculates the ratio q of the geometric progression
         call get_gp_ratio(nb, L/ar, q)

         sb(0) = 0.d0

         do i = 1, nb-1

            sb(i) = sb(i-1) + ar * q ** (nb-i)

         end do

         sb(nb) = L

      end subroutine get_gp_grid_right



      !> \brief Calculates the distribution of points according
      !! to a geometric progression. The length of the first partition
      !! is defined on the left and on the right.
      subroutine get_gp_grid_left_right(nb, al, ar, L, sb)
         implicit none
         integer, intent(in)  :: nb       !< Number of partitions over the new distribution of points
         real(8), intent(in)  :: al       !< Length of the first patition on the left
         real(8), intent(in)  :: ar       !< Length of the first patition on the right
         real(8), intent(in)  :: L        !< Length of the curve
         real(8), intent(out) :: sb(0:nb) !< Length distribution over the curve

         ! Inner parameters

         integer, parameter :: nitmax = 10000

         ! Inner variables

         integer :: i
         integer :: nl
         integer :: nr
         real(8) :: Lm
         real(8) :: Li
         real(8) :: Lf
         real(8) :: acl
         real(8) :: acr


         nl = nb / 2

         nr = nb - nl

         Li = 0.d0

         Lf = L

         do i = 1, nitmax

            Lm = 0.5d0 * ( Li + Lf )

            call get_gp_grid_left (nl, al, Lm, sb( 0:nl))

            call get_gp_grid_right(nr, ar, L-Lm, sb(nl:nb))

            sb(nl:nb) = sb(nl:nb) + Lm

            acl = sb(nl) - sb(nl-1)

            acr = sb(nl+1) - sb(nl)

            if ( acl > acr ) then

               Lf = Lm

            else

               Li = Lm

            end if

            if ( abs(acl-acr) < 1.d-14 ) exit

         end do

      end subroutine get_gp_grid_left_right


      !> \brief Calculates a distribution of points based on the arclength and
      !! on the double-exponential rule
      subroutine get_dexp_m1_ditribution( nbr, nb, aks, zcv, xbr, ybr, xb, yb) ! Output: last two
         implicit none
         integer, intent(in)  :: nbr        !< Number of partitions over the reference curve (xbr,ybr)
         integer, intent(in)  :: nb         !< Number of partitions for the new distribution of points
         real(8), intent(in)  :: aks        !< Exponent for the calculation of acv
         real(8), intent(in)  :: zcv        !< Central value of z
         real(8), intent(in)  :: xbr(0:nbr) !< Given distribution of points
         real(8), intent(in)  :: ybr(0:nbr) !< Given distribution of points
         real(8), intent(out) :: xb(0:nb)   !< New distribution of points
         real(8), intent(out) :: yb(0:nb)   !< New distribution of points

         ! Inner variables

         real(8) :: sbr(0:nbr) ! Length of each partition over the reference curve
         real(8) :: sb(0:nb)   ! Length of each partition over the new curve
         real(8) :: L          ! Total length
         real(8) :: acv        ! Length of the volume closer to the central point


         ! Calculates the length of the curve
         call get_sbr(nbr, xbr, ybr, sbr) ! Output: last one


         ! Total length
         L = sbr(nbr)


         ! Calculating the double-exponential distribution of arclength

         acv = L / ( nb ** aks )

         if ( acv > zcv ) then

            write(*,*) "Error. get_dexp_m1_ditribution: there is no resolution"&
               , " for this grid. Stopping."

            stop

         end if

         call get_dexp_m1(nb, 0.d0, zcv, L, acv, sb)


         ! Calculates the new distribution of points
         call get_new_distribution(nbr, nb, sbr, xbr, ybr, sb, xb, yb) ! Output: last two


      end subroutine get_dexp_m1_ditribution

      !> \brief Calculates a distribution of points based on the double-exponential
      !! rule. The coefficients of this rule are set by defining the position of
      !! an internal point (zcv) and the width (acv) of the volumes beside this
      !! point. The points are concentrated around zcv.
      subroutine get_dexp_m1(np, ziv, zcv, zfv, acv, z)
         implicit none
         integer, intent(in)  :: np  !< Number of partitions
         real(8), intent(in)  :: ziv !< Initial value of z
         real(8), intent(in)  :: zcv !< Central value of z (exponential's matching point)
         real(8), intent(in)  :: zfv !< Final value of z
         real(8), intent(in)  :: acv !< Length of the volumes beside zcv
         real(8), intent(out) :: z(0:np) !< Distribution of points

         ! Inner variables

         integer :: i
         integer :: npc
         real(8) :: t
         real(8) :: dt
         real(8) :: tc
         real(8) :: r
         real(8) :: A
         real(8) :: ac
         real(8) :: zc


         ! Initializing variables

         zc = zcv / ( zfv - ziv )

         ac = acv / ( zfv - ziv )

         dt = 1.d0 / dble(np)

         r = ac / ( 1.d0 - zc )



         ! Cheking conditions

         if ( ( ac / zc > 1.d0 ) .or. ( r > 1.d0 ) ) then

            write(*,*) "Error. get_dexp_m1: ac > zc or ac > (1-zc). Stopping."

            stop

         end if

         if ( ac > dt ) then

            write(*,*) "Error. get_dexp_m1: ac > dt. Stopping."

            stop

         end if


         ! Searching the value of A

         call adexp_m1(dt, r, zc, A)


         ! Calculating the fraction of volumes for each exponential

         tc = 1.d0 + log( 1.d0 + ( exp(-A*dt) - 1.d0 ) / r ) / A

         npc = nint(tc*np)



         ! Calculating the normallized distribution of points

         ! First exponential

         do i = 0, npc

            t = tc * dble(i) / dble(npc)

            z(i) = zc * ( exp(A*t) - 1.d0 ) / ( exp(A*tc) - 1.d0 )

         end do

         ! Second exponential

         do i = npc, np

            t = ( 1.d0-tc ) * dble(i-npc) / dble(np-npc) + tc

            z(i) = zc+(1.d0-zc)*(exp(-A*(t-tc))-1.d0)/(exp(-A*(1.d0-tc))-1.d0)

         end do

         z(0)  = 0.d0

         z(np) = 1.d0

         ! Calculating the extended distribution of points

         z = ( zfv - ziv ) * z + ziv

      end subroutine get_dexp_m1

      !> \brief Defines the transcendental function for finding A
      real(8) function fdexp_m1(dt, r, zc, A)
         implicit none
         real(8), intent(in) :: dt
         real(8), intent(in) :: r
         real(8), intent(in) :: zc
         real(8), intent(in) :: A

         ! Inner varialbles

         real(8) :: raux1, raux2, raux3

         if ( abs(A) > 1.d-6 ) then

            raux1 = exp(A) * ( 1.d0 + ( exp(-A*dt) - 1.d0 ) / r )
            raux2 = (2.d0*zc-1.d0) / zc
            raux3 = (1.d0-zc) / zc

            fdexp_m1 = -exp(-A) * raux1 * raux1 +  raux1 * raux2 + raux3

         else

            fdexp_m1 = (r*zc-r+dt)*A/(r*zc) &
               + ((r**2-2*dt**2)*zc-r**2+(2*dt-dt**2)*r)*A**2/(r**2*zc)/2.0

         end if

      end function fdexp_m1


      !> \brief Finds the value of A using the bissection method
      subroutine adexp_m1(dt, r, zc, A)
         implicit none
         real(8), intent(in)  :: dt
         real(8), intent(in)  :: r
         real(8), intent(in)  :: zc
         real(8), intent(out) :: A

         ! Parameters

         integer, parameter :: itmax = 1000

         ! Inner variables

         integer :: k
         real(8) :: Ai, Af, Am, Fi, Fm, Ff



         ! Initializing the search interval

         Ai = -1.d-18

         Af = -100.d0



         ! Checking for roots

         Fi = fdexp_m1(dt, r, zc, Ai)

         Ff = fdexp_m1(dt, r, zc, Af)

         if ( Fi * Ff > 0.d0 ) then

            write(*,*) "Error. adexp_m1: there is none or more than one roots. Stopping."

            stop

         end if


         ! Looking for a solution

         do k = 1, itmax

            Am = 0.5d0 * ( Ai + Af )

            Fi = fdexp_m1(dt, r, zc, Ai)

            Fm = fdexp_m1(dt, r, zc, Am)

            if ( Fi * Fm <= 0.d0 ) then

               Af = Am

            else

               Ai = Am

            end if

            if ( abs(Af-Ai) < 1.d-14 ) exit

         end do

         if ( k >= itmax ) then

            write(*,*) "Error. adexp_m1: number of iteractions exceeded. Stopping."

            stop

         end if

         A = 0.5d0 * ( Ai + Af )

      end subroutine adexp_m1


   end subroutine get_grid_boundary_g23

end module
