!>
!! \brief Contains procedures for post-processing data
!!
module mod_postp_procedures

   use mod_class_thermophysical_abstract
   use mod_system

contains

   !> \brief Plots the residual as a function of the iteractions.
   subroutine plotter02( sem_g, sim_id)
      implicit none
      integer, intent(in) :: sem_g !< 0 = visualize the plot, 1 = do not visualize
      character (len = *), intent(in) :: sim_id    !< Simulation identification

      ! Auxiliary variables
      character (len=250) :: str1, str2, str3

      !
      ! Generating gnuplot file
      !

      str1 = "./gnu/mach2d_SIM_ID_residual.gnu"

      str2 = "./mach2d_output/mach2d_" // trim(adjustl(sim_id)) //"_residual.gnu"

      str3 = "./mach2d_output/mach2d_" // trim(adjustl(sim_id)) //"_residual.eps"

      ! Replaces all SIM_ID string in the original file by sim_id
      call file_replace( str1, str2, "SIM_ID", trim(adjustl(sim_id)) )


      call system( graph_generator // trim(adjustl(str2)) )

      if ( sem_g == 0 ) call system( graph_viewer // trim(adjustl(str3)) )

   end subroutine plotter02


   !> \brief Generates the boundary grid plot.
   subroutine plotter04( sem_g, sim_id)
      implicit none
      integer, intent(in) :: sem_g !< 0 = visualize the plot, 1 = do not visualize
      character (len = *), intent(in) :: sim_id    !< Simulation identification

      ! Auxiliary variables
      character (len=250) :: str1, str2, str3

      !
      ! Generating gnuplot file
      !

      ! Template file
      str1 = "./gnu/mach2d_SIM_ID_boundary.gnu"

      ! Case file
      str2 = "./mach2d_output/mach2d_" // trim(adjustl(sim_id)) // '_boundary.gnu'

      ! Output file
      str3 = "./mach2d_output/mach2d_" // trim(adjustl(sim_id)) // '_boundary.eps'

      ! Replaces all SIM_ID string in the original file by sim_id
      call file_replace( str1, str2, "SIM_ID", trim(adjustl(sim_id)) )

      call system( graph_generator // trim(adjustl(str2)) )

      if ( sem_g == 0 ) call system( graph_viewer // trim(adjustl(str3)) )

   end subroutine plotter04


   !> \brief Generates the grid plot.
   subroutine plotter05( sem_g, sim_id)
      implicit none
      integer, intent(in) :: sem_g !< 0 = visualize the plot, 1 = do not visualize
      character (len = *), intent(in) :: sim_id    !< Simulation identification

      ! Auxiliary variables
      character (len=250) :: str1, str2, str3

      !
      ! Generating gnuplot file
      !

      ! Template file
      str1 = "./gnu/mach2d_SIM_ID_grid.gnu"

      ! Case file
      str2 = "./mach2d_output/mach2d_" // trim(adjustl(sim_id)) // '_grid.gnu'

      ! Output file
      str3 = "./mach2d_output/mach2d_" // trim(adjustl(sim_id)) // '_grid.eps'

      ! Replaces all SIM_ID string in the original file by sim_id
      call file_replace( str1, str2, "SIM_ID", trim(adjustl(sim_id)) )

      call system( graph_generator // trim(adjustl(str2)) )

      if ( sem_g == 0 ) call system( graph_viewer // trim(adjustl(str3)) )

   end subroutine plotter05


   !> \brief Generates the plot of the following fields: p, ro, T, u, v, M.
   subroutine plotter06( sem_g, sim_id)
      implicit none
      integer, intent(in) :: sem_g !< 0 = visualize the plot, 1 = do not visualize
      character (len = *), intent(in) :: sim_id    !< Simulation identification

      ! Auxiliary variables
      character (len=250) :: str1, str2, str3

      !
      ! Generating gnuplot file
      !

      ! Template file
      str1 = "./gnu/mach2d_SIM_ID_main_fields.gnu"

      ! Case file
      str2 = "./mach2d_output/mach2d_" // trim(adjustl(sim_id)) // '_main_fields.gnu'

      ! Replaces all SIM_ID string in the original file by sim_id
      call file_replace( str1, str2, "SIM_ID", trim(adjustl(sim_id)) )

      ! Auxiliary string
      str3 = "./mach2d_output/mach2d_" // trim(adjustl(sim_id))

      call system( graph_generator // trim(adjustl(str2)) )

      if ( sem_g == 0 ) call system( graph_viewer // trim(adjustl(str3)) // "_p_field.eps" )
      if ( sem_g == 0 ) call system( graph_viewer // trim(adjustl(str3)) // "_ro_field.eps" )
      if ( sem_g == 0 ) call system( graph_viewer // trim(adjustl(str3)) // "_T_field.eps" )
      if ( sem_g == 0 ) call system( graph_viewer // trim(adjustl(str3)) // "_u_field.eps" )
      if ( sem_g == 0 ) call system( graph_viewer // trim(adjustl(str3)) // "_v_field.eps" )
      if ( sem_g == 0 ) call system( graph_viewer // trim(adjustl(str3)) // "_M_field.eps" )

   end subroutine plotter06


   !> \brief Generates the plots of u, v, T and p over the boundaries.
   subroutine plotter07( sim_id)
      implicit none
      character (len = *), intent(in) :: sim_id    !< Simulation identification

      ! Auxiliary variables
      character (len=250) :: str1, str2

      !
      ! Generating gnuplot file
      !

      str1 = "./gnu/mach2d_SIM_ID_main_fields_boundaries.gnu"

      str2 = "./mach2d_output/mach2d_" // trim(adjustl(sim_id)) //"_main_fields_boundaries.gnu"

      ! Replaces all SIM_ID string in the original file by sim_id
      call file_replace( str1, str2, "SIM_ID", trim(adjustl(sim_id)) )

      call system( graph_generator // trim(adjustl(str2)) )

   end subroutine plotter07


   !> \brief Plots the convergence coefficients of the linear systems for
   !! u, v, T and p'
   subroutine plotter08(nx, ny, sim_id, xp, yp, ccu, ccv, cct, ccp)
      implicit none
      integer, intent(in) :: nx  ! Number of volumes in csi direction (real+fictitious)
      integer, intent(in) :: ny  ! Number of volumes in eta direction (real+fictitious)
      character (len = 100), intent(in) :: sim_id ! Simulation identification
      real(8), dimension (nx*ny), intent(in) :: xp ! Coord. of the centroid of volume P
      real(8), dimension (nx*ny), intent(in) :: yp ! Coord. of the centroid of volume P
      real(8), dimension (nx*ny), intent(in) :: ccu ! Convergence coef. for u
      real(8), dimension (nx*ny), intent(in) :: ccv ! Convergence coef. for v
      real(8), dimension (nx*ny), intent(in) :: cct ! Convergence coef. for T
      real(8), dimension (nx*ny), intent(in) :: ccp ! Convergence coef. for p

      ! Inner variables

      integer :: i, j, np
      character (len = 100) :: str1, str2

      str1 = "./mach2d_output/mach2d_"      &
         // trim(adjustl(sim_id)) // "_" &
         // "convergence_coefficients.dat"

      open(10, file = str1)


      ! Printing the fields to a file

      write(10,'(6(A15))') '# xp', 'yp', 'ccu', 'ccv', 'cct', 'ccp'

      do i = 2, nx-1

         do j = 2, ny-1

            np   = nx * (j-1) + i

            write(10,'(6(ES15.7))') xp(np), yp(np), ccu(np), ccv(np), cct(np) &
            , ccp(np)

         end do

         write(10,*)

      end do

      close(10)

      !
      ! Generating gnuplot file
      !

      ! Template file
      str1 = "./gnu/mach2d_SIM_ID_convergence_coefficients.gnu"

      ! Case file
      str2 = "./mach2d_output/mach2d_" // trim(adjustl(sim_id)) &
      // '_convergence_coefficients.gnu'

      ! Replaces all SIM_ID string in the original file by sim_id
      call file_replace( str1, str2, "SIM_ID", trim(adjustl(sim_id)) )

      call system( graph_generator // trim(adjustl(str2)) )

   end subroutine


   !> \brief Writes the gas composition and its properties for a reference
   !! temperature
   subroutine write_gas_properties(unt, Tr, thermomodel)
      implicit none
      integer,                              intent(in) :: unt         !< Unit were the results will be printed
      real(8),                              intent(in) :: Tr          !< Temperature of reference (K)
      class(class_thermophysical_abstract), intent(in) :: thermomodel !< The thermophysical model

      ! Inner variables

      real(8) :: Prr ! Prandtl number of reference
      real(8) :: Cpr ! Cp of reference (J/kg.K)
      real(8) :: mur ! Viscosity of reference (Pa.s)
      real(8) :: kpr ! Thermal conductivity of reference (W/m.K)

      ! Getting thermophysical properties for the reference temperature
      mur = thermomodel%mu(Tr)
      kpr = thermomodel%kp(Tr)
      Cpr = thermomodel%cp(Tr)

      Prr = Cpr * mur / kpr

      ! Printing information about thermophysical model
      call thermomodel%about(unt)

      write(unt,*)
      write(unt,*) "Gas properties for a temperature of reference"
      write(unt,*)
      write(unt,"(ES23.16,2X,A)")  Tr, " =  Tr: temperature of reference (K)"
      write(unt,"(ES23.16,2X,A)") Cpr, " = Cpr: Cp(Tr) (J/kg.K)"
      write(unt,"(ES23.16,2X,A)") mur, " = mur: mu(Tr) (Pa.s)"
      write(unt,"(ES23.16,2X,A)") kpr, " = kpr: kp(Tr) (W/m.K)"
      write(unt,"(ES23.16,2X,A)") Prr, " = Prr: Pr(Tr)"
      write(unt,*)

   end subroutine write_gas_properties


   !> \brief Writes grid parameters
   subroutine write_grid_parameters( nx, ny, unit, x, y)
      implicit none
      integer, intent(in) :: nx   ! Number of volumes in csi direction (real+fictitious)
      integer, intent(in) :: ny   ! Number of volumes in eta direction (real+fictitious)
      integer, intent(in) :: unit ! Unit to which the results will be printed
      real(8), dimension (nx*ny), intent(in) :: x   ! Coord. x of the northest corner of volume P
      real(8), dimension (nx*ny), intent(in) :: y   ! Coord. y of the northeast corner of volume P

      ! Auxiliary variables
      integer :: i, j, np, npn
      real(8) :: dy, dx, dxmax, dxmin, dymax, dymin
      real(8) :: a1min, a1max, raux, xmin, xmax


      dxmin = maxval(abs(x))
      dxmax = minval(abs(x))
      dymin = maxval(abs(y))
      dymax = minval(abs(y))

      do j = 1, ny-1

         do i = 2, nx-1

            np = (j-1)*nx + i

            dx = dabs ( x(np) - x(np-1) )

            if ( dx < dxmin ) dxmin = dx

            if ( dx > dxmax ) dxmax = dx

         end do

      end do

      do i = 1, nx-1

         do j = 2, ny-1

            np = (j-1)*nx + i

            dy = dabs ( y(np) - y(np-nx) )

            if ( dy < dymin ) dymin = dy

            if ( dy > dymax ) dymax = dy

         end do

      end do

      write(unit,1) dxmin, dxmax, dxmin/dxmax, dymin, dymax, dymin/dymax

      1   format(//,'*** Grid parameters ***', //, &
         1pe22.9,' = dxmin: min(size of volume controls) in the x direction (m)',/, &
         1pe22.9,' = dxmax: max(size of volume controls) in the x direction (m)',/, &
         1pe22.9,' = dxmin / dxmax',/, &
         1pe22.9,' = dymin: min(size of volume controls) in the y direction (m)',/, &
         1pe22.9,' = dymax: max(size of volume controls) in the y direction (m)',/, &
         1pe22.9,' = dymin / dymax')


      ! Searching for the smallest and largest width of the volumes closer to the south boundary

      j = 1

      a1min = 1.d10

      a1max = 0.d0

      do i = 1, nx-1

         np  = nx * (j-1) + i

         npn = np + nx

         raux = sqrt( (x(npn)-x(np))**2 + (y(npn)-y(np))**2 )

         if ( raux < a1min ) then

            a1min = raux

            xmin = x(np)

         end if

         if ( raux > a1max ) then

            a1max = raux

            xmax = x(np)

         end if

      end do

      write(unit,"(ES23.16,A)") a1min, " = a1min: minimum width of the volume closer to the body in the current grid (m)"
      write(unit,"(ES23.16,A)") a1max, " = a1max: maximum width of the volume closer to the body in the current grid (m)"
      write(unit,"(ES23.16,A)")  xmin, " =  xmin: x coord. where a1min first occured (m)"
      write(unit,"(ES23.16,A)")  xmax, " =  xmax: x coord. where a1max first occured (m)"


   end subroutine write_grid_parameters


   !> \brief Writes boundary points of the grid
   subroutine write_grid_boundary( nx, ny, w_g, sim_id, x, y)
      implicit none
      integer, intent(in) :: nx   ! Number of volumes in csi direction (real+fictitious)
      integer, intent(in) :: ny   ! Number of volumes in eta direction (real+fictitious)
      integer, intent(in) :: w_g   ! Frequency of writing data for graphics
      character (len = *), intent(in) :: sim_id     ! Simulation identification
      real(8), dimension (nx*ny), intent(in) :: x   ! Coord. x of the northest corner of volume P
      real(8), dimension (nx*ny), intent(in) :: y   ! Coord. y of the northeast corner of volume P

      ! Auxiliary variables
      integer :: i, j, np

      ! writes south boundary

      open(7,file="./mach2d_output/mach2d_" // trim(adjustl(sim_id))//'_boundary_south.dat')

      write(7,1)

      j = 1
      do i = 1, nx-1
         np = (j-1)*nx + i
         if ( i==1 .or. i==nx-1 .or. mod(i,w_g)==0 ) &
            write(7,2) i, j, np, x(np), y(np)
      end do

      close(7)

      ! writes north boundary

      open(7,file="./mach2d_output/mach2d_" // trim(adjustl(sim_id))//'_boundary_north.dat')

      write(7,1)

      j = ny-1
      do i = 1, nx-1
         np = (j-1)*nx + i
         if ( i==1 .or. i==nx-1 .or. mod(i,w_g)==0 ) &
            write(7,2) i, j, np, x(np), y(np)
      end do

      close(7)

      ! writes west boundary

      open(7,file="./mach2d_output/mach2d_" // trim(adjustl(sim_id))//'_boundary_west.dat')

      write(7,1)

      i = 1
      do j = 1, ny-1
         np = (j-1)*nx + i
         if ( j==1 .or. j==ny-1 .or. mod(j,w_g)==0 ) &
            write(7,2) i, j, np, x(np), y(np)
      end do

      close(7)

      ! writes east boundary

      open(7,file="./mach2d_output/mach2d_" // trim(adjustl(sim_id))//'_boundary_east.dat')

      write(7,1)

      i = nx-1
      do j = 1, ny-1
         np = (j-1)*nx + i
         if ( j==1 .or. j==ny-1 .or. mod(j,w_g)==0 ) &
            write(7,2) i, j, np, x(np), y(np)
      end do

      close(7)

      1   format('#', t5,'i', t10,'j', t19,'np', t24,'x (m)', t48,'y (m)')

      2   format(2(i5),i10,2(1pe24.15))

   end subroutine write_grid_boundary


   !> \brief Writes the grid
   subroutine write_grid(nx, ny, w_g, sim_id, x, y)
      implicit none
      integer, intent(in) :: nx   ! Number of volumes in csi direction (real+fictitious)
      integer, intent(in) :: ny   ! Number of volumes in eta direction (real+fictitious)
      integer, intent(in) :: w_g   ! Frequency of writing data for graphics
      character (len = *), intent(in) :: sim_id     ! Simulation identification
      real(8), dimension (nx*ny), intent(in) :: x   ! Coord. x of the northest corner of volume P
      real(8), dimension (nx*ny), intent(in) :: y   ! Coord. y of the northeast corner of volume P

      ! Auxiliary variables
      integer :: i, j, np


      open(7,file= "./mach2d_output/mach2d_" // trim(adjustl(sim_id))//'_grid.dat')

      ! write the lines between west and east boundaries

      do j = 1, ny-1

         write(7,*)

         do i = 1, nx-1
            np = (j-1)*nx + i
            if ( i==1 .or. i==nx-1 .or. mod(i,w_g)==0 ) &
               write(7,2) x(np), y(np)
         end do

      end do

      ! write the lines between south and north boundaries

      do i = 1, nx-1

         write(7,*)

         do j = 1, ny-1
            np = (j-1)*nx + i
            if ( j==1 .or. j==ny-1 .or. mod(j,w_g)==0 ) &
               write(7,2) x(np), y(np)
         end do

      end do

      close(7)

      2   format(2(1pe24.15))

   end subroutine write_grid


   !> \brief Write a field b
   subroutine write_field(nx, ny, unit, xp, yp, b)
      implicit none
      integer, intent(in) :: nx   ! Number of volumes in csi direction (real+fictitious)
      integer, intent(in) :: ny   ! Number of volumes in eta direction (real+fictitious)
      integer, intent(in) :: unit ! Unit in which the field will be written
      real(8), dimension (nx*ny), intent(in) :: xp  ! Coord. of the centroid of volume P
      real(8), dimension (nx*ny), intent(in) :: yp  ! Coord. of the centroid of volume P
      real(8), dimension (nx*ny), intent(in) :: b   ! Field

      ! Auxiliary variables
      integer :: i, j, k, nt, maxi, i1, i2, np
      integer :: nloc(1)

      maxi = 10

      nt = nx / maxi

      do k = 1, nt

         i1 = (k-1) * maxi + 1

         i2 = k * maxi

         write(unit,"(A4,3X,10(I14))") "i = ", ( i, i = i1, i2 )

         do j = ny, 1, -1
            write(unit,"(A4,I3,10(ES14.6))") "j = ", j, (b(np), np = nx * (j-1) + i1, nx * (j-1) + i2)
         end do

         write(unit,*)

      end do

      i1 = nt * maxi + 1

      i2 = nx

      if ( i1 <= i2 ) then

         write(unit,"(A4,3X,10(I14))") "i = ", ( i, i = i1, i2 )

         do j = ny, 1, -1
            write(unit,"(A4,I3,10(ES14.6))") "j = ", j, (b(np), np = nx * (j-1) + i1, nx * (j-1) + i2)
         end do

      end if

      nloc = minloc(b)

      write(unit,*)
      write(unit,*) "Minimum of the field"
      write(unit,*)
      write(unit,"((ES14.6,A9,5X),2(ES14.6,A5,5X))") b(nloc(1)), " = Minval", xp(nloc(1)), " = xp", yp(nloc(1)), " = yp"


      nloc = maxloc(b)

      write(unit,*)
      write(unit,*) "Maximum of the field"
      write(unit,*)

      write(unit,"((ES14.6,A9,5X),2(ES14.6,A5,5X))") b(nloc(1)), " = Maxval", xp(nloc(1)), " = xp", yp(nloc(1)), " = yp"

      write(unit,*)
   end subroutine write_field


   !> \brief Write the main fields
   subroutine write_main_fields(nx, ny, sim_id, xp, yp, p, ro, T, u, v, M)
      implicit none
      integer, intent(in) :: nx  ! Number of volumes in csi direction (real+fictitious)
      integer, intent(in) :: ny  ! Number of volumes in eta direction (real+fictitious)
      character (len = 100), intent(in) :: sim_id ! Simulation identification
      real(8), dimension (nx*ny), intent(in) :: xp ! Coord. of the centroid of volume P
      real(8), dimension (nx*ny), intent(in) :: yp ! Coord. of the centroid of volume P
      real(8), dimension (nx*ny), intent(in) :: p  ! Pressure at center o volume P
      real(8), dimension (nx*ny), intent(in) :: ro ! Absolute density at center of vol. P
      real(8), dimension (nx*ny), intent(in) :: T  ! Temperature of the last iteraction
      real(8), dimension (nx*ny), intent(in) :: u  ! Cartesian velocity of the last iteraction
      real(8), dimension (nx*ny), intent(in) :: v  ! Cartesian velocity of the last iteraction
      real(8), dimension (nx*ny), intent(in) :: M  ! Mach number field

      ! Auxiliary variables

      integer :: i, j, np
      character (len = 100) :: str1

      str1 = "./mach2d_output/mach2d_"      &
         // trim(adjustl(sim_id)) // "_" &
         // "main_fields.dat"

      open(10, file = str1)

      write(10,"(8(X,A11))") "# x", "y", "p", "ro", "T", "u", "v", "M"

      do j = 1, ny
         do i = 1, nx
            np = (j-1)*nx + i
            write(10,"(8(X,ES11.4))") xp(np), yp(np), p(np), ro(np), T(np), u(np), v(np), M(np)
         end do
         write(10,*)
      end do

      close(10)

   end subroutine write_main_fields


   !> \brief Prints the fields u, v, p and T on the boundaries.
   subroutine write_fields_on_boundaries( nx, ny, sim_id, xp, yp, u, v, p, T ) ! Output: none
      implicit none
      integer, intent(in) :: nx   !< Number of volumes in the csi direction (real+fictitious)
      integer, intent(in) :: ny   !< Number of volumes in the eta direction (real+fictitious)
      character (len = 100), intent(in) :: sim_id !< Simulation identification
      real(8), dimension (nx*ny), intent(in) :: xp  !< Coord. x of the centroid of volume P (m)
      real(8), dimension (nx*ny), intent(in) :: yp  !< Coord. y of the centroid of volume P (m)
      real(8), dimension (nx*ny), intent(in) :: u   !< Cartesian velocity at the center of volume P (m/s)
      real(8), dimension (nx*ny), intent(in) :: v   !< Cartesian velocity at the center of volume P (m/s)
      real(8), dimension (nx*ny), intent(in) :: p   !< Pressure at center o volume P (Pa)
      real(8), dimension (nx*ny), intent(in) :: T   !< Temperature at center of volume P (K)

      ! Auxiliary variables
      integer :: i, j, np
      character (len = 100) :: str1

      ! west boundary (frontal symmetry line)

      str1 = "./mach2d_output/mach2d_"      &
         // trim(adjustl(sim_id)) // "_" &
         // "main_fields_west_boundary.dat"

      open(10,file = str1)

      write(10,"(8(X,A23))") "# x", "y", "u", "v", "p", "T"

      i = 1
      do j = 1, ny

         np = (j-1)*nx + i

         write(10,"(6(X,ES23.16))") xp(np), yp(np), u(np), v(np), p(np), T(np)

      end do

      close(10)

      ! east boundary (exit)

      str1 = "./mach2d_output/mach2d_"      &
         // trim(adjustl(sim_id)) // "_" &
         // "main_fields_east_boundary.dat"

      open(10,file = str1)

      write(10,"(8(X,A23))") "# x", "y", "u", "v", "p", "T"

      i = nx

      do j = 1, ny

         np = (j-1)*nx + i

         write(10,"(6(X,ES23.16))") xp(np), yp(np), u(np), v(np), p(np), T(np)

      end do

      close(10)


      ! lower boundary (body wall)

      str1 = "./mach2d_output/mach2d_"      &
         // trim(adjustl(sim_id)) // "_" &
         // "main_fields_south_boundary.dat"

      open(10,file = str1)

      write(10,"(8(X,A23))") "# x", "y", "u", "v", "p", "T"

      j = 1

      do i = 1, nx

         np   = (j-1)*nx + i

         write(10,"(6(X,ES23.16))") xp(np), yp(np), u(np), v(np), p(np), T(np)

      end do

      close(10)

      ! upper boundary (far field)

      str1 = "./mach2d_output/mach2d_"      &
         // trim(adjustl(sim_id)) // "_" &
         // "main_fields_north_boundary.dat"

      open(10,file = str1)

      write(10,"(8(X,A23))") "# x", "y", "u", "v", "p", "T"

      j = ny

      do i = 1, nx

         np  = (j-1)*nx + i

         write(10,"(6(X,ES23.16))") xp(np), yp(np), u(np), v(np), p(np), T(np)

      end do

      close(10)

   end subroutine


   !> \brief Post-process fields on boundaries
   subroutine post_processing_boundaries( nx, ny, x, y, xp, yp, u, v, p, T ) ! InOutput: last six entries
      implicit none
      integer, intent(in) :: nx   ! Number of volumes in csi direction (real+fictitious)
      integer, intent(in) :: ny   ! Number of volumes in eta direction (real+fictitious)
      real(8), dimension (nx*ny), intent(in)  :: x     ! Coord. x of the northest corner of volume P
      real(8), dimension (nx*ny), intent(in)  :: y     ! Coord. y of the northeast corner of volume P
      real(8), dimension (nx*ny), intent(inout) :: xp  ! Coord. of the centroid of volume P
      real(8), dimension (nx*ny), intent(inout) :: yp  ! Coord. of the centroid of volume P
      real(8), dimension (nx*ny), intent(inout) :: u   ! Cartesian velocity of the last iteraction
      real(8), dimension (nx*ny), intent(inout) :: v   ! Cartesian velocity of the last iteraction
      real(8), dimension (nx*ny), intent(inout) :: p   ! Pressure at center o volume P
      real(8), dimension (nx*ny), intent(inout) :: T   ! Temperature at center o volume P

      ! Auxiliary variables
      integer :: i, j, np, npn, npw, nps, npe, npsw, npse, npnw, npne

      ! west boundary (frontal symmetry line)
      i = 1
      do j = 2, ny-1
         np = (j-1)*nx + i
         nps = np - nx
         npe  = np + 1
         u(np) = ( u(np) + u(npe) ) / 2.d0
         v(np) = ( v(np) + v(npe) ) / 2.d0
         p(np) = ( p(np) + p(npe) ) / 2.d0
         T(np) = ( T(np) + T(npe) ) / 2.d0
         xp(np) = ( x(nps) + x(np) ) / 2.d0
         yp(np) = ( y(nps) + y(np) ) / 2.d0
      end do

      ! east boundary (exit)
      i = nx
      do j = 2, ny-1
         np = (j-1)*nx + i
         npw  = np - 1
         u(np) = ( u(npw) + u(np) ) / 2
         v(np) = ( v(npw) + v(np) ) / 2
         p(np) = ( p(npw) + p(np) ) / 2
         T(np) = ( T(npw) + T(np) ) / 2
         xp(np) = ( x(npw-nx) + x(npw) ) / 2
         yp(np) = ( y(npw-nx) + y(npw) ) / 2
      end do

      ! lower boundary (body wall)
      j = 1
      do i = 2, nx-1
         np   = (j-1)*nx + i
         npw  = np - 1
         npn  = np + nx
         u(np) = ( u(npn) + u(np) ) / 2.d0
         v(np) = ( v(npn) + v(np) ) / 2.d0
         p(np) = ( p(npn) + p(np) ) / 2.d0
         T(np) = ( T(npn) + T(np) ) / 2.d0
         xp(np) = ( x(npw) + x(np) ) / 2.d0
         yp(np) = ( y(npw) + y(np) ) / 2.d0
      end do

      ! SW corner
      j = 1
      i = 1
      np   = (j-1)*nx + i
      npn  = np + nx
      npe  = np + 1
      npne = npn + 1

      u(np) = ( u(npe) + u(npn) + u(npne) ) / 3.d0
      v(np) = ( v(npe) + v(npn) + v(npne) ) / 3.d0
      p(np) = ( p(npe) + p(npn) + p(npne) ) / 3.d0
      T(np) = ( T(npe) + T(npn) + T(npne) ) / 3.d0
      xp(np) = x(np)
      yp(np) = y(np)

      ! SE corner
      j = 1
      i = nx
      np   = (j-1)*nx + i
      npn  = np + nx
      npw  = np - 1
      npnw = npn - 1

      u(np) = ( u(npw) + u(npn) + u(npnw) ) / 3.d0
      v(np) = ( v(npw) + v(npn) + v(npnw) ) / 3.d0
      p(np) = ( p(npw) + p(npn) + p(npnw) ) / 3.d0
      T(np) = ( T(npw) + T(npn) + T(npnw) ) / 3.d0
      xp(np) = x(npw)
      yp(np) = y(npw)

      ! upper boundary (far field)
      j = ny
      do i = 2, nx-1
         np  = (j-1)*nx + i
         nps = np - nx
         u(np) = (u(np) + u(nps)) / 2
         v(np) = (v(np) + v(nps)) / 2
         p(np) = (p(np) + p(nps)) / 2
         T(np) = (T(np) + T(nps)) / 2
         xp(np) = ( x(nps) + x(nps-1) ) / 2
         yp(np) = ( y(nps) + y(nps-1) ) / 2
      end do

      ! NW corner
      j = ny
      i = 1
      np   = (j-1)*nx + i
      nps  = np - nx
      npe  = np + 1
      npse = nps + 1
      u(np) = ( u(npe) + u(nps) + u(npse) ) / 3.d0
      v(np) = ( v(npe) + v(nps) + v(npse) ) / 3.d0
      p(np) = ( p(npe) + p(nps) + p(npse) ) / 3.d0
      T(np) = ( T(npe) + T(nps) + T(npse) ) / 3.d0
      xp(np) = x(nps)
      yp(np) = y(nps)


      ! NE corner
      j = ny
      i = nx
      np   = (j-1)*nx + i
      nps  = np - nx
      npw  = np - 1
      npsw = nps - 1
      u(np) = ( u(npw) + u(nps) + u(npsw) ) / 3.d0
      v(np) = ( v(npw) + v(nps) + v(npsw) ) / 3.d0
      p(np) = ( p(npw) + p(nps) + p(npsw) ) / 3.d0
      T(np) = ( T(npw) + T(nps) + T(npsw) ) / 3.d0
      xp(np) = x(npsw)
      yp(np) = y(npsw)


   end subroutine post_processing_boundaries

end module
