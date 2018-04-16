
!> \brief Calculates the geometry model 'newton2coef', that is a modified
!! Newton geometry with two adjustable coefficients

module newton2coef
   implicit none


contains

   !> \brief Calculates the geometry model 'newton2coef'
   subroutine get_newton2coef_body(nr, rb, lo, h, zf, xr, yr) ! Output: last two
      implicit none
      integer, intent(in)  :: nr !< Number of partitions
      real(8), intent(in)  :: rb !< Base radius of the body
      real(8), intent(in)  :: lo !< Length of the body
      real(8), intent(in)  :: h  !< Frontal to base radius ratio: rf/rb
      real(8), intent(in)  :: zf !< Frontal slope: dy/dx
      real(8), intent(out) :: xr(0:nr) !< x coordinate
      real(8), intent(out) :: yr(0:nr) !< y coordinate

      ! Inner variables

      integer :: i
      real(8) :: A, B, C, z, dz, zb, fr


      ! Calculating the fineness ratio
      fr = lo / ( 2.d0 * rb )


      ! Calculating the slope at base
      call get_zb(h, fr, zf, zb) ! Output: last one


      ! Calculating some coefficients
      A = - p(zf) / ( p(zb) - p(zf) )

      B = ( h * w(zb) - w(zf) ) / ( 2.d0 * fr * ( w(zb) - w(zf) ) )

      C = 1.d0 / ( p(zb) - p(zf) )


      ! Calculating the dimensionless x and y

      dz = (zb-zf) / dble(nr-1)

      xr(0) = 0.d0

      yr(0) = 0.d0

      xr(1) = 0.d0

      yr(1) = h / ( 2.d0 * fr )

      do i = 2, nr-1

         z = dz * dble(i-1) + zf

         xr(i) = C * p(z) + A

         yr(i) = C * w(z) + B

      end do

      xr(nr) = 1.d0

      yr(nr) = 1.d0 / ( 2.d0 * fr )


      ! Calculating the dimensional x and y

      xr = xr * lo

      yr = yr * lo

   end subroutine get_newton2coef_body


   !> \brief Calculates function p(z)
   real(8) function p(z)
      implicit none
      real(8), intent(in) :: z !< Slope

      p = 0.75d0 / z**4 + 1.d0 / z**2 + log(z)

   end function p


   !> \brief Calculates function w(z)
   real(8) function w(z)
      implicit none
      real(8), intent(in) :: z !< Slope

      w = ( 1.d0 + z*z )**2 / z**3

   end function w



   !> \brief Calculates the slope at base
   subroutine get_zb(h, fr, zf, zb) ! Output: last one
      implicit none
      real(8), intent(in)  :: h  !< Frontal to base radius ratio: rf/rb
      real(8), intent(in)  :: fr !< Fineness ratio
      real(8), intent(in)  :: zf !< Frontal slope
      real(8), intent(out) :: zb !< Base slope

      ! Parameters

      integer, parameter :: nitm_root = 500
      real(8), parameter :: tol_root = 1.d-14

      call bissection(1.d-5, zf-1.d-5, nitm_root, tol_root, zb) ! Output: last one

      contains

         !> \brief Calculates function f(z)
         real(8) function f(zb)
            implicit none
            real(8), intent(in) :: zb !< Base slope

            f = (h-1.d0) * ( p(zb) - p(zf) ) - 2.d0*fr * ( w(zf) - w(zb) )

         end function f

         !> \brief Finds the root of function f(x) by the bissection method
         subroutine bissection(xi, xf, nitm_root, tol_root, xr) ! Output: last one
            implicit none
            real(8), intent(in)  :: xi !< Initial value of x
            real(8), intent(in)  :: xf !< Final value of x
            integer, intent(in)  :: nitm_root !< Maximum number of iteraction for the root calculation
            real(8), intent(in)  :: tol_root  !< Tolerance for the error
            real(8), intent(out) :: xr !< Root

            ! Inner variables

            integer :: i
            real(8) :: xa, xm, xb, fa, fm, fb

            ! Initializing

            xa = xi
            xb = xf

            fa = f(xa)
            fb = f(xb)

            ! Checking the number of roots in the interval [xi,xf]

            if ( fa * fb > 0.d0 ) then

               write(*,*) "No roots or more than one. Stopping..."

               stop

            end if


            ! Searching the root

            do i = 1, nitm_root

               xm = 0.5d0 * ( xa + xb )

               fm = f(xm)

               if ( fa * fm > 0.d0 ) then

                  xa = xm

                  fa = fm

               else

                  xb = xm

                  fb = fm

               end if

               if ( abs(xb-xa) < tol_root .or. abs(fb-fa) < tol_root ) exit

            end do


            ! Checking the number of iteractions

            if ( i >= nitm_root ) then

               write(*,*) "Number of iteractions exceeded. Stopping..."

               stop

            end if

            ! Final solution

            xr = 0.5d0 * ( xa + xb )

         end subroutine bissection

   end subroutine get_zb


end module newton2coef
