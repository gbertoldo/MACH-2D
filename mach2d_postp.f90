module postp

   use mod_class_thermophysical_abstract

   character (len = *), parameter :: text_editor     = "emacs -nw "
   character (len = *), parameter :: graph_viewer    = "evince "
   character (len = *), parameter :: graph_generator = "gnuplot "

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

      integer :: i   ! Dummy index
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


   !> \brief Writes some of the free stream properties of the gas
   subroutine write_free_stream_properties(unt, lr, rb, PF, TF, MF, CPF, GF &
      , VLF, KPF, PRF, ROF, UF, REFm)
      implicit none
      integer, intent(in) ::  unt !< Unit where the results will be printed
      real(8), intent(in) ::   lr !< Length of the body (m)
      real(8), intent(in) ::   rb !< Base radius of the body (m)
      real(8), intent(in) ::   PF !< Free stream pressure (Pa)
      real(8), intent(in) ::   TF !< Free stream temperature (K)
      real(8), intent(in) ::   MF !< Free stream Mach number
      real(8), intent(in) ::  CPF !< Free stream Cp (J/kg.K)
      real(8), intent(in) ::   GF !< Free stream gamma
      real(8), intent(in) ::  VLF !< Free stream viscosity (Pa.s)
      real(8), intent(in) ::  KPF !< Free stream thermal conductivity (W/m.K)
      real(8), intent(in) ::  PRF !< Free stream Prandtl number
      real(8), intent(in) ::  ROF !< Free stream density (kg/m3)
      real(8), intent(in) ::   UF !< Free stream speed (m/s)
      real(8), intent(in) :: REFm !< Free stream Reynolds number per meter (1/m)

      ! Parameters

      real(8), parameter :: pi = dacos(-1.d0)

      ! Inner variables

      real(8) :: REFl ! Free stream Reynolds number (length based)
      real(8) :: REFr ! Free stream Reynolds number (radius based)
      real(8) :: KNFl ! Knudsen number based on REFl
      real(8) :: KNFr ! Knudsen number based on REFr

      ! Free stream Reynolds number (length based)
      REFl = REFm * lr

      ! Free stream Reynolds number (radius based)
      REFr = REFm * rb

      ! Knudsen number based on REFl
      KNFl = sqrt( pi * GF / 2.d0 ) * MF / REFl

      ! Knudsen number based on REFr
      KNFr = sqrt( pi * GF / 2.d0 ) * MF / REFr

      write(unt,*)
      write(unt,*) " *** Free stream properties *** "
      write(unt,*)

      write(unt,"(ES23.16,A)")   GF, " =   GF: free stream gamma"
      write(unt,"(ES23.16,A)")  CPF, " =  CPF: free stream Cp (J/kg.K)"
      write(unt,"(ES23.16,A)")  VLF, " =  VLF: free stream viscosity (Pa.s)"
      write(unt,"(ES23.16,A)")  KPF, " =  KPF: free stream thermal conductivity (W/m.K)"
      write(unt,"(ES23.16,A)")  PRF, " =  PRF: free stream Prandtl number"
      write(unt,"(ES23.16,A)")   PF, " =   PF: free stream pressure (Pa)"
      write(unt,"(ES23.16,A)")   TF, " =   TF: free stream temperature (K)"
      write(unt,"(ES23.16,A)")  ROF, " =  ROF: free stream density (kg/m3)"
      write(unt,"(ES23.16,A)")   MF, " =   MF: free stream Mach number"
      write(unt,"(ES23.16,A)")   UF, " =   UF: free stream speed (m/s)"
      write(unt,"(ES23.16,A)") REFm, " = REFm: free stream Reynolds number per meter (1/m)"
      write(unt,"(ES23.16,A)") REFl, " = REFl: free stream Reynolds number (length based)"
      write(unt,"(ES23.16,A)") REFr, " = REFr: free stream Reynolds number (radius based)"
      write(unt,"(ES23.16,A)") KNFl, " = KNFl: Knudsen number based on REFl"
      write(unt,"(ES23.16,A)") KNFr, " = KNFr: Knudsen number based on REFr"

   end subroutine write_free_stream_properties


   !> \brief Writes main parameters of the convergence coefficients of the
   !! linear systems for u, v, T and p'
   subroutine write_cc_parameters(unit,nx, ny, ccu, ccv, cct, ccp)
      implicit none
      integer, intent(in) :: unit ! Unit to which the results will be printed
      integer, intent(in) :: nx  ! Number of volumes in csi direction (real+fictitious)
      integer, intent(in) :: ny  ! Number of volumes in eta direction (real+fictitious)
      real(8), dimension (nx*ny), intent(in) :: ccu ! Convergence coef. for u
      real(8), dimension (nx*ny), intent(in) :: ccv ! Convergence coef. for v
      real(8), dimension (nx*ny), intent(in) :: cct ! Convergence coef. for T
      real(8), dimension (nx*ny), intent(in) :: ccp ! Convergence coef. for p

      ! Inner variables

      integer :: i, j, np

      write(unit,*)
      write(unit,*) "*** Convergence coefficients that exceed 1 ***"

      write(unit,*)
      write(unit,*) "*** ccu ***"
      write(unit,*)
      write(unit,'(A4,A4,(A15))') 'i', 'j', 'ccu'
      do i = 2, nx-1

         do j = 2, ny-1

            np   = nx * (j-1) + i

            if ( ccu(np) > 1.d0 ) then
               write(unit,'(I4,I4,(ES15.7))') i, j, ccu(np)
            end if

         end do

      end do


      write(unit,*)
      write(unit,*) "*** ccv ***"
      write(unit,*)
      write(unit,'(A4,A4,(A15))') 'i', 'j', 'ccv'
      do i = 2, nx-1

         do j = 2, ny-1

            np   = nx * (j-1) + i

            if ( ccv(np) > 1.d0 ) then
               write(unit,'(I4,I4,(ES15.7))') i, j, ccv(np)
            end if

         end do

      end do

      write(unit,*)
      write(unit,*) "*** ccT ***"
      write(unit,*)
      write(unit,'(A4,A4,(A15))') 'i', 'j', 'ccT'
      do i = 2, nx-1

         do j = 2, ny-1

            np   = nx * (j-1) + i

            if ( ccT(np) > 1.d0 ) then
               write(unit,'(I4,I4,(ES15.7))') i, j, ccT(np)
            end if

         end do

      end do

      write(unit,*)
      write(unit,*) "*** ccp ***"
      write(unit,*)
      write(unit,'(A4,A4,(A15))') 'i', 'j', 'ccp'
      do i = 2, nx-1

         do j = 2, ny-1

            np   = nx * (j-1) + i

            if ( ccp(np) > 1.d0 ) then
               write(unit,'(I4,I4,(ES15.7))') i, j, ccp(np)
            end if

         end do

      end do

   end subroutine


   subroutine write_body_parameters( unit, lr, rb) ! Output: last entry
      implicit none
      integer, intent(in) :: unit ! Unit to which the results will be printed
      real(8), intent(in) :: lr   ! length of the rocket
      real(8), intent(in) :: rb   ! base radius of the rocket

      write(unit,*)
      write(unit,*)
      write(unit,*) '*** Body parameters ***'
      write(unit,*)
      write(unit,"(ES23.16,A)") lr, " = lr: Body length (m)"
      write(unit,"(ES23.16,A)") rb, " = rb: Body base radius (m)"
      write(unit,"(ES23.16,A)") lr/(2*rb), " = lr / (2 rb): Body length to base diameter ratio (dimensionless)"


   end subroutine write_body_parameters


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


   ! Print main data to a file
   subroutine write_main_results(unit, it, norm, dt, tcpu, RAM, curef &
      , cvref, ctref, cpref, Cdfi, Cdfv)
      implicit none
      integer, intent(in) :: unit  ! Unit to which the results will be printed
      integer, intent(in) :: it     ! final number of iteractions for time evolution
      real(8), intent(in) :: norm   ! Norm L1 of the residuals of the last it.
      real(8), intent(in) :: dt     ! dt of the last iteraction (s)"
      real(8), intent(in) :: tcpu   ! Total time (before and after simulation interruption)
      real(8), intent(in) :: RAM    ! RAM memory (MB)
      real(8), intent(in) :: curef  ! Reference value of the coefficient of convergence for u
      real(8), intent(in) :: cvref  ! Reference value of the coefficient of convergence for v
      real(8), intent(in) :: ctref  ! Reference value of the coefficient of convergence for T
      real(8), intent(in) :: cpref  ! Reference value of the coefficient of convergence for p
      real(8), intent(in) :: Cdfi   ! Pressure foredrag coefficient
      real(8), intent(in) :: Cdfv   ! Viscous foredrag coefficient

      write(unit,*)
      write(unit,*)
      write(unit,*) "*** Main results ***"
      write(unit,*)
      write(unit,"(ES23.15,A)")     norm, ' =  norm: sum of the relative norm' &
         // ' L1 of the residuals of the linear systems for u, v, T and pl'
      write(unit,"(ES23.15,A)")      RAM, " =   RAM: Memory (MB)"
      write(unit,"(ES23.15,A)")       dt, " =    dt: dt of the last iteraction (s)"
      write(unit,"(I23,A)")         it-1, " =    it: final number of iteractions for time evolution"
      write(unit,"(ES23.15,A)")     tcpu, " =  tcpu: cpu time (s)"
      write(unit,"(ES23.15,A)")    curef, " = curef: Reference value of the coefficient of convergence for u"
      write(unit,"(ES23.15,A)")    cvref, " = cvref: Reference value of the coefficient of convergence for v"
      write(unit,"(ES23.15,A)")    ctref, " = ctref: Reference value of the coefficient of convergence for T"
      write(unit,"(ES23.15,A)")    cpref, " = cpref: Reference value of the coefficient of convergence for p"
      write(unit,"(ES23.15,A)")     Cdfi, " =  Cdfi: Pressure foredrag coefficient"
      write(unit,"(ES23.15,A)")     Cdfv, " =  Cdfv: Viscous foredrag coefficient"
      write(unit,"(ES23.15,A)") Cdfi+Cdfv, " =   Cdf: Foredrag coefficient"

   end subroutine


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


   !> \brief Finds all the occurences of str1 in file1, replaces them by str2 and save the results in file2
   subroutine file_replace( file1, file2, str1, str2)
      implicit none
      character (len=*), intent(in) :: file1 !< Input file
      character (len=*), intent(in) :: file2 !< Output file
      character (len=*), intent(in) :: str1  !< Find
      character (len=*), intent(in) :: str2  !< Replace

      ! Auxiliary variables
      integer :: io
      character (len=1000) :: str

      if ( file1 /= file2 ) then

         open(10, file = file1)
         open(11, file = file2)

         do

            read(10,"(A)",iostat=io) str

            if ( io /= 0 ) exit

            call replace(str,str1,str2)

            write(11,*) trim(adjustl(str))

         end do

         close(10)
         close(11)

      else

         open(10, file = file1)
         open(11, file = trim(adjustl(file2))//".tmp")

         do

            read(10,"(A)",iostat=io) str

            if ( io /= 0 ) exit

            call replace(str,str1,str2)

            write(11,*) trim(adjustl(str))

         end do

         close(10)
         close(11)

         call rename( trim(adjustl(file2))//".tmp", file1 )

      end if

   end subroutine file_replace

   !> \brief Finds all the occurrences of str2 in str1 and replaces them by str3
   subroutine replace(str1, str2, str3)
      implicit none
      character (len=*), intent(inout) :: str1 !< Original string
      character (len=*), intent(in)    :: str2 !< Find
      character (len=*), intent(in)    :: str3 !< Replace

      ! Auxiliary variables

      integer :: idx !< First position where str2 was found
      integer :: n1  !< Length of str1
      integer :: n2  !< Length of str2

      n1 = len(str1)

      n2 = len(str2)

      do

         idx = index(str1,str2)

         if ( idx /= 0 ) then

            str1 = str1(1:idx-1) // str3 // str1(idx+n2:n1)

         else

            exit

         end if

      end do

   end subroutine replace

end module postp
