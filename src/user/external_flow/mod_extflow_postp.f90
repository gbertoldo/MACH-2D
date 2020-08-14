!>
!! \brief Performs data post-processing
!!
module mod_extflow_postp

   use mod_postp_procedures
   use depp_interface
   use mod_thermophysical

   implicit none

   ! Makes everything private, except otherwise stated
   private

   ! Public procedures
   public :: extflow_postp

contains

   !> \brief Post-processed data
   subroutine extflow_postp()
      use mod_mach2d_data
      implicit none

      ! Saves post processed data
      if ( wppd == 1 ) then

         ! Post processing data

         call post_processing

      ! Simplified post-processing
      else if ( wppd == 2 ) then

         call simplified_post_processing

      end if

   end subroutine



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


   !> \brief Writes body main parameters
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


   !> \brief Print main data to a file
   subroutine write_main_results(unit, it, norm, dt, tcpu, RAM, Cdfi, Cdfv)
      implicit none
      integer, intent(in) :: unit  ! Unit to which the results will be printed
      integer, intent(in) :: it     ! final number of iteractions for time evolution
      real(8), intent(in) :: norm   ! Norm L1 of the residuals of the last it.
      real(8), intent(in) :: dt     ! dt of the last iteraction (s)"
      real(8), intent(in) :: tcpu   ! Total time (before and after simulation interruption)
      real(8), intent(in) :: RAM    ! RAM memory (MB)
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
      write(unit,"(ES23.15,A)")     Cdfi, " =  Cdfi: Pressure foredrag coefficient"
      write(unit,"(ES23.15,A)")     Cdfv, " =  Cdfv: Viscous foredrag coefficient"
      write(unit,"(ES23.15,A)") Cdfi+Cdfv, " =   Cdf: Foredrag coefficient"

   end subroutine


   subroutine post_processing

      use mod_extflow_data

      use mod_mach2d_data

      implicit none

      ! Auxiliary variables
      real(8), dimension (nx*ny) :: M    ! Mach number at the center of the volumes

      ! Saves the fitness function to a file
      ! Exit status: 0 = success, 1 = failure, 2 = generate another individual
      if ( it < itmax ) then

         call depp_save_fitness(-(Cdfi+Cdfv), 0, "Fitness: -(Cdfi+Cdfv)")

      else

         call depp_save_fitness(-1.d5, 1, "Fitness: -(Cdfi+Cdfv) (no convergence)")

      end if

      ! Saving fields to a file for further comparisons

      open(15,file='./mach2d_output/mainfields' // trim(sim_id) // '.dat')

      write(15,*) nx
      write(15,*) ny
      write(15,*) it
      write(15,*) tcpu
      write(15,*) Cdfi
      write(15,*) Cdfv
      write(15,*) u
      write(15,*) v
      write(15,*) T
      write(15,*) p

      close(15)

      ! Writes the gas composition and its properties for a reference temperature
      call write_gas_properties(lid, TF, thermomodel)

      ! Writes some of the free stream properties of the gas
      call write_free_stream_properties(lid, lr, rb, PF, TF, MF, CPF, GF &
         , VLF, KPF, PRF, ROF, UF, REFm)


      ! Writes grid parameters
      call write_grid_parameters( nx, ny, lid, x, y) ! No output

      ! Writes body parameters
      call write_body_parameters( lid, lr, rb) ! Output: last entry

      ! Writes grid boundary nodes
      call write_grid_boundary( nx, ny, w_g, sim_id, x, y)

      ! Plots grid boundaries
      call plotter04( sem_g, sim_id)

      ! Writes grid vertices
      call write_grid(nx, ny, w_g, sim_id, x, y)

      ! Plots the grid
      call plotter05( sem_g, sim_id)


      ! Calculates the boundary values and store them in the fictitious
      call post_processing_boundaries( nx, ny, x, y, xp, yp, u, v, p, T ) ! InOutput: last six entries

      ! Calculates gamma for real volumes and boundaries
      call set_gamma(thermomodel, nx, ny, Rg, cp, gcp) ! Output: last one

      ! Calculates the Mach number field
      M = dsqrt ( (u**2 + v**2) / (gcp * Rg * T) )

      ! Calculates the specific mass field
      ro = p / ( Rg * T )

      ! Write fields

      call write_main_fields(nx, ny, sim_id, xp, yp, p, ro, T, u, v, M)

      call plotter06( sem_g, sim_id)

      call write_fields_on_boundaries( nx, ny, sim_id, xp, yp, u, v, p, T ) ! Output: none

      call plotter07( sim_id)

      if ( w_cam == 1 ) then

         write(lid,09)
         09     format(//,10x,'*** FIELDS ***', //, &
            10x,'The fictitious nodes have the boundary values.', //, &
            10x,'Horizon (from left to right): boundary west to east.', //, &
            10x,'Vertical (from top to bottom): boundary north to south.' )

         write(lid,10)
         10     format(//,30x,'M: Mach number (dimensionless)')
         call write_field(nx, ny, lid, xp, yp, M)

         write(lid,11)
         11     format(//,30x,'u: nodal cartesian velocity u (m/s)')
         call write_field(nx, ny, lid, xp, yp, u)

         write(lid,12)
         12     format(//,30x,'v: nodal cartesian velocity v (m/s)')
         call write_field(nx, ny, lid, xp, yp, v)

         write(lid,13)
         13     format(//,30x,'p: pressure (Pa)')
         call write_field(nx, ny, lid, xp, yp, p)

         write(lid,14)
         14     format(//,30x,'T: Temperature (K)')
         call write_field(nx, ny, lid, xp, yp, T)

         write(lid,15)
         15     format(//,30x,'ro: nodal specific mass (kg/m3)')
         call write_field(nx, ny, lid, xp, yp, ro)

         write(lid,55)
         55     format(//,30x,'cp: specific heat at constant pressure (J/kg.K)')
         call write_field(nx, ny, lid, xp, yp, cp)

         write(lid,58)
         58     format(//,30x,'gcp: specific heat ratio (dimensionless)')
         call write_field(nx, ny, lid, xp, yp, gcp)

         write(lid,56)
         56     format(//,30x,'vlp: nodal laminar viscosity (Pa.s)')
         call write_field(nx, ny, lid, xp, yp, vlp)

         write(lid,57)
         57     format(//,30x,'kp: nodal thermal conductivity (W/m.K)')
         call write_field(nx, ny, lid, xp, yp, kp)

         write(lid,16)
         16     format(//,30x,'Uce: contravariant velocity U at east face (m2/s)')
         call write_field(nx, ny, lid, xp, yp, Uce)

         write(lid,17)
         17     format(//,30x,'Vcn: contravariant velocity V at north face (m2/s)')
         call write_field(nx, ny, lid, xp, yp, Vcn)

         write(lid,18)
         18     format(//,30x,'pl: pressure variation (Pa)')
         call write_field(nx, ny, lid, xp, yp, pl)

         write(lid,19)
         19     format(//,30x,'bp: source of pl (kg/s)')
         call write_field(nx, ny, lid, xp, yp, bp)

         write(lid,21)
         21     format(//,30x,'X: x coordinate at northeast corner (m)')
         call write_field(nx, ny, lid, xp, yp, x)

         write(lid,22)
         22     format(//,30x,'Y: y coordinate at northeast corner (m)')
         call write_field(nx, ny, lid, xp, yp, y)

         write(lid,23)
         23     format(//,30x,'Jp: Jacobian at the center of volume (1/m2)')
         call write_field(nx, ny, lid, xp, yp, Jp)

         write(lid,24)
         24     format(//,30x,'xe: metric x eta at east face (m)')
         call write_field(nx, ny, lid, xp, yp, xe)

         write(lid,25)
         25     format(//,30x,'ye: metric y eta at east face (m)')
         call write_field(nx, ny, lid, xp, yp, ye)

         write(lid,26)
         26     format(//,30x,'xk: metric x ksi at north face (m)')
         call write_field(nx, ny, lid, xp, yp, xk)

         write(lid,27)
         27     format(//,30x,'yk: metric y ksi at north face (m)')
         call write_field(nx, ny, lid, xp, yp, yk)

         write(lid,28)
         28     format(//,30x,'de: coefficient d of SIMPLEC at east face (m3.s/kg)')
         call write_field(nx, ny, lid, xp, yp, de)

         write(lid,29)
         29     format(//,30x,'dn: coefficient d of SIMPLEC at north face (m3.s/kg)')
         call write_field(nx, ny, lid, xp, yp, dn)

      end if


      ! Print main results
      call write_main_results(lid, it, norm, dt, tcpu, RAM, Cdfi, Cdfv)


      ! Plots the residual as a function of the iteractions
      call plotter02( sem_g, sim_id)

      ! show the listing file of iteractions
      if ( sem_a == 0 ) call system(text_editor // './mach2d_output/mach2d_' // trim(adjustl(sim_id)) // '_residual.dat')

      ! close the main listing file
      close (lid)

      ! show main listing file
      if ( sem_a == 0 ) call system(text_editor // './mach2d_output/mach2d_' // trim(adjustl(sim_id)) // '.txt')


      ! Plotting the convergence coefficients
      ! call plotter08(nx, ny, sim_id, xp, yp, ccu, ccv, cct, ccp)

   end subroutine post_processing



   subroutine simplified_post_processing

      use mod_extflow_data
      use mod_mach2d_data

      implicit none

      ! Saves the fitness function to a file
      ! Exit status: 0 = success, 1 = failure, 2 = generate another individual
      if ( it < itmax ) then

         call depp_save_fitness(-(Cdfi+Cdfv), 0, "Fitness: -(Cdfi+Cdfv)")

      else

         call depp_save_fitness(-1.d5, 1, "Fitness: -(Cdfi+Cdfv) (no convergence)")

      end if

      ! Writes the gas composition and its properties for a reference temperature
      call write_gas_properties(lid, TF, thermomodel)


      ! Writes some of the free stream properties of the gas
      call write_free_stream_properties(lid, lr, rb, PF, TF, MF, CPF, GF &
         , VLF, KPF, PRF, ROF, UF, REFm)


      ! Writes grid parameters
      call write_grid_parameters( nx, ny, lid, x, y) ! No output

      ! Writes body parameters
      call write_body_parameters( lid, lr, rb) ! Output: last entry


      ! Print main results
      call write_main_results(lid, it, norm, dt, tcpu, RAM, Cdfi, Cdfv)

      ! close the main listing file
      close (lid)

   end subroutine simplified_post_processing

end module
