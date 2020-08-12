!> \brief Interpolates a set of data points using splines
module spline
   implicit none

contains

   !> \brief Calculates the natural cubic splines
   subroutine get_cspline(n, m, xi, yi, xo, yo) ! Output: last one
      implicit none
      integer, intent(in)  :: n !< Number of intervals for xi
      integer, intent(in)  :: m !< Number of intervals for xo
      real(8), intent(in)  :: xi(0:n) !< x-coordinates - input
      real(8), intent(in)  :: yi(0:n) !< y-coordinates - input
      real(8), intent(in)  :: xo(0:m) !< x-coordinates - output
      real(8), intent(out) :: yo(0:m) !< y-coordinates - output

      ! Inner variables

      integer :: i, j
      real(8) :: b(0:n) ! b coefficients
      real(8) :: c(0:n) ! c coefficients
      real(8) :: d(0:n) ! d coefficients


      ! Checking the constraints
      if ( xo(0) < xi(0) .or. xi(n) < xo(m) ) then

         write(*,*) "Error: out of range. Stopping..."

         stop

      end if


      ! Calculating the csplines coefficients
      call get_csplines_coeff(n, xi, yi, b, c, d) ! Output: last three


      j = 0

      do i = 0, m

         do

            if ( xo(i) <= xi(j+1) ) then

               yo(i) = yi(j)                       &
                  + b(j) * ( xo(i)-xi(j) )         &
                  + c(j) * ( xo(i)-xi(j) ) ** 2.d0 &
                  + d(j) * ( xo(i)-xi(j) ) ** 3.d0

               exit

            else

               j = j + 1

            end if

         end do

      end do

   end subroutine get_cspline


   !> \brief Calculates the natural cubic splines coefficients
   subroutine get_csplines_coeff(n, x, y, b, c, d) ! Output: last three
      implicit none
      integer, intent(in)  :: n !< Number of intervals
      real(8), intent(in)  :: x(0:n) !< x-coordinates
      real(8), intent(in)  :: y(0:n) !< y-coordinates
      real(8), intent(out) :: b(0:n) !< b coefficients
      real(8), intent(out) :: c(0:n) !< c coefficients
      real(8), intent(out) :: d(0:n) !< d coefficients

      ! Inner variables

      integer :: i
      real(8) :: h(0:n-1)
      real(8) :: bp(0:n)
      real(8) :: ap(0:n,3)


      ! Calculating h

      do i = 0, n-1

         h(i) = x(i+1)-x(i)

      end do


      ! Setting west boundary condition

      i = 0

      ap(i,1) = 0.d0
      ap(i,2) = 1.d0
      ap(i,3) = 0.d0
      bp(i) = 0.d0

      do i = 1, n-1

         ap(i,1) = h(i-1)
         ap(i,2) = 2.d0 * ( h(i) + h(i-1) )
         ap(i,3) = h(i)
         bp(i) = 3.d0 * ( (y(i+1)-y(i)) / h(i) - (y(i)-y(i-1)) / h(i-1) )

      end do

      ! Setting east boundary condition

      i = n

      ap(i,1) = 0.d0
      ap(i,2) = 1.d0
      ap(i,3) = 0.d0
      bp(i) = 0.d0

      ! Solving the linear system

      call tdma(n+1, ap, bp, c)

      ! Calculating the other coefficients

      do i = 0, n-1

         b(i) = ( y(i+1)-y(i) ) / h(i) - h(i) * ( c(i+1) + 2.d0 * c(i) ) / 3.d0

         d(i) = ( c(i+1)-c(i) ) / ( 3.d0 * h(i) )

      end do

      b(n) = b(n-1) + h(n-1) * ( c(n) + c(n-1) )

      d(n) = 0.d0

   end subroutine get_csplines_coeff


   !> \brief Calculates the clamped cubic splines
   subroutine get_ccspline(n, m, fpo, fpn, xi, yi, xo, yo) ! Output: last one
      implicit none
      integer, intent(in)  :: n !< Number of intervals for xi
      integer, intent(in)  :: m !< Number of intervals for xo
      real(8), intent(in)  :: fpo     !< f'(xi(0))
      real(8), intent(in)  :: fpn     !< f'(xi(n))
      real(8), intent(in)  :: xi(0:n) !< x-coordinates - input
      real(8), intent(in)  :: yi(0:n) !< y-coordinates - input
      real(8), intent(in)  :: xo(0:m) !< x-coordinates - output
      real(8), intent(out) :: yo(0:m) !< y-coordinates - output

      ! Inner variables

      integer :: i, j
      real(8) :: b(0:n) ! b coefficients
      real(8) :: c(0:n) ! c coefficients
      real(8) :: d(0:n) ! d coefficients


      ! Checking the constraints
      if ( xo(0) < xi(0) .or. xi(n) < xo(m) ) then

         write(*,*) "Error: out of range. Stopping..."

         stop

      end if


      ! Calculating the clamped cubic splines coefficients
      call get_ccsplines_coeff(n, fpo, fpn, xi, yi, b, c, d) ! Output: last three


      j = 0

      do i = 0, m

         do

            if ( xo(i) <= xi(j+1) ) then

               yo(i) = yi(j)                       &
                  + b(j) * ( xo(i)-xi(j) )         &
                  + c(j) * ( xo(i)-xi(j) ) ** 2.d0 &
                  + d(j) * ( xo(i)-xi(j) ) ** 3.d0

               exit

            else

               j = j + 1

            end if

         end do

      end do

   end subroutine get_ccspline


   !> \brief Calculates the clamped cubic splines coefficients
   subroutine get_ccsplines_coeff(n, fpo, fpn, x, y, b, c, d) ! Output: last three
      implicit none
      integer, intent(in)  :: n !< Number of intervals
      real(8), intent(in)  :: fpo    !< f'(xi(0))
      real(8), intent(in)  :: fpn    !< f'(xi(n))
      real(8), intent(in)  :: x(0:n) !< x-coordinates
      real(8), intent(in)  :: y(0:n) !< y-coordinates
      real(8), intent(out) :: b(0:n) !< b coefficients
      real(8), intent(out) :: c(0:n) !< c coefficients
      real(8), intent(out) :: d(0:n) !< d coefficients

      ! Inner variables

      integer :: i
      real(8) :: h(0:n-1)
      real(8) :: bp(0:n)
      real(8) :: ap(0:n,3)


      ! Calculating h

      do i = 0, n-1

         h(i) = x(i+1)-x(i)

      end do


      ! Setting west boundary condition

      i = 0

      ap(i,1) = 0.d0
      ap(i,2) = h(0) / 3.d0 * 2.d0
      ap(i,3) = h(0) / 3.d0
      bp(i) = (y(1)-y(0)) / h(0) - fpo

      do i = 1, n-1

         ap(i,1) = h(i-1)
         ap(i,2) = 2.d0 * ( h(i) + h(i-1) )
         ap(i,3) = h(i)
         bp(i) = 3.d0 * ( (y(i+1)-y(i)) / h(i) - (y(i)-y(i-1)) / h(i-1) )

      end do

      ! Setting east boundary condition

      i = n

      ap(i,1) = h(n-1) / 3.d0
      ap(i,2) = h(n-1) / 3.d0 * 2.d0
      ap(i,3) = 0.d0
      bp(i) = fpn - (y(n)-y(n-1)) / h(n-1)

      ! Solving the linear system

      call tdma(n+1, ap, bp, c)

      ! Calculating the other coefficients

      do i = 0, n-1

         b(i) = ( y(i+1)-y(i) ) / h(i) - h(i) * ( c(i+1) + 2.d0 * c(i) ) / 3.d0

         d(i) = ( c(i+1)-c(i) ) / ( 3.d0 * h(i) )

      end do

      b(n) = b(n-1) + h(n-1) * ( c(n) + c(n-1) )

      d(n) = 0.d0

   end subroutine get_ccsplines_coeff



   !> \brief Calculates the cubic splines with second derivatives prescribed
   !! on all points and the function values given on the boundaries
   subroutine get_d2cspline(n, m, xi, y0, yn, dy2, xo, yo) ! Output: last one
      implicit none
      integer, intent(in)  :: n !< Number of intervals for xi
      integer, intent(in)  :: m !< Number of intervals for xo
      real(8), intent(in)  :: xi(0:n)  !< x-coordinates - input
      real(8), intent(in)  :: y0       !< f(x0)
      real(8), intent(in)  :: yn       !< f(xn)
      real(8), intent(in)  :: dy2(0:n) !< diff( f(x), x, 2 )
      real(8), intent(in)  :: xo(0:m)  !< x-coordinates - output
      real(8), intent(out) :: yo(0:m)  !< y-coordinates - output

      ! Inner variables

      integer :: i, j
      real(8) :: a(0:n) ! a coefficients
      real(8) :: b(0:n) ! b coefficients
      real(8) :: c(0:n) ! c coefficients
      real(8) :: d(0:n) ! d coefficients


      ! Checking the constraints
      if ( xo(0) < xi(0) .or. xi(n) < xo(m) ) then

         write(*,*) "Error: out of range. Stopping..."

         stop

      end if


      ! Calculates the cubic splines coefficients with second derivatives prescribed
      ! and function prescribed on boundaries
      call get_d2csplines_coeff(n, xi, y0, yn, dy2, a, b, c, d) ! Output: last four


      j = 0

      do i = 0, m

         do

            if ( xo(i) <= xi(j+1) ) then

               yo(i) = a(j)                        &
                  + b(j) * ( xo(i)-xi(j) )         &
                  + c(j) * ( xo(i)-xi(j) ) ** 2.d0 &
                  + d(j) * ( xo(i)-xi(j) ) ** 3.d0

               exit

            else

               j = j + 1

            end if

         end do

      end do

   end subroutine get_d2cspline


   !> \brief Calculates the cubic splines coefficients with second derivatives prescribed
   !! and function prescribed on boundaries
   subroutine get_d2csplines_coeff(n, x, y0, yn, dy2, a, b, c, d) ! Output: last four
      implicit none
      integer, intent(in)  :: n !< Number of intervals
      real(8), intent(in)  :: x(0:n)   !< x-coordinates
      real(8), intent(in)  :: y0       !< f(x0)
      real(8), intent(in)  :: yn       !< f(xn)
      real(8), intent(in)  :: dy2(0:n) !< diff( f(x), x, 2 )
      real(8), intent(out) :: a(0:n)   !< a coefficients
      real(8), intent(out) :: b(0:n)   !< b coefficients
      real(8), intent(out) :: c(0:n)   !< c coefficients
      real(8), intent(out) :: d(0:n)   !< d coefficients

      ! Inner variables

      integer :: i
      real(8) :: h(0:n-1)
      real(8) :: bp(0:n)
      real(8) :: ap(0:n,3)


      ! Calculating h

      do i = 0, n-1

         h(i) = x(i+1)-x(i)

      end do

      ! Calculating c

      do i = 0, n

         c(i) = dy2(i) / 2.d0

      end do

      ! Calculating d

      do i = 0, n-1

         d(i) = ( c(i+1) - c(i) ) / ( 3.d0 * h(i) )

      end do

      d(n) = 0.d0

      ! Calculating a

      ! Setting west boundary condition

      i = 0

      ap(i,1) = 0.d0
      ap(i,2) = 1.d0
      ap(i,3) = 0.d0
      bp(i) = y0

      do i = 1, n-1

         ap(i,1) = 1.d0 / h(i-1)
         ap(i,2) = - ( 1.d0 / h(i) + 1.d0 / h(i-1) )
         ap(i,3) = 1.d0 / h(i)
         bp(i) = (c(i+1)+2.d0*c(i)) * h(i) / 3.d0 + (2.d0*c(i)+c(i-1)) * h(i-1) / 3.d0

      end do

      ! Setting east boundary condition

      i = n

      ap(i,1) = 0.d0
      ap(i,2) = 1.d0
      ap(i,3) = 0.d0
      bp(i) = yn

      ! Solving the linear system

      call tdma(n+1, ap, bp, a)

      ! Calculating the other coefficients

      i = 0

      b(i) = (a(i+1)-a(i)) / h(i) - c(i)*h(i) - d(i) * h(i)**2

      do i = 0, n-1

         b(i+1) = b(i) + 2.d0 * c(i) * h(i) + 3.d0 * d(i) * h(i)**2

      end do

   end subroutine get_d2csplines_coeff




   subroutine tdma(n, a, b, x)
      implicit none
      integer, intent(in) :: n ! Number unknowns
      real(8), dimension(n,3), intent(in)  :: a ! Tri-diagonal matrix
      real(8), dimension(n),   intent(in)  :: b ! Source
      real(8), dimension(n),   intent(out) :: x ! Solution
      !
      ! a(i,1) = west coefficients
      ! a(i,2) = central coefficients
      ! a(i,3) = east coefficients
      !
      ! Auxiliary variables
      integer :: i
      real(8), dimension(n) :: P
      real(8), dimension(n) :: Q

      i = 1

      P(i) = - a(i,3) / a(i,2)

      Q(i) = b(i) / a(i,2)

      do i = 2, n

       P(i) = - a(i,3) / ( a(i,2) + a(i,1) * P(i-1) )

       Q(i) = ( b(i) - a(i,1) * Q(i-1) ) / ( a(i,2) + a(i,1) * P(i-1) )

      end do

      i = n
      x(i) = Q(i)

      do i = n-1, 1, -1
       x(i) = x(i+1) * P(i) + Q(i)
      end do

   end subroutine tdma

end module spline
