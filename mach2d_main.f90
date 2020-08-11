
program main

   use data
   use mod_grid
   use mod_grid_data, only: rb, lr
   use user
   use coefficients
   use solvers
   use postp
   use mod_tstepper

   implicit none


   ! Reads parameters, allocates variables and sets vectors to zero
   call initialize_variables


   ! Initializes variables related to the thermophysical properties of the gas
   call thermophysical_initialization(ifile, thermomodel, Rg, ktm) ! Output: last three


   ! Calculates some of the free stream properties of the gas
   call get_free_stream_properties(thermomodel, PF, TF, MF, CPF, GF, VLF, KPF, PRF, ROF, UF, REFm, HF) ! Output: last 9


   ! Creating the grid and calculating the metrics
   call grid_create(lid, coord, REFm                                  & ! Input
         ,  x, y, xp, yp, xe, ye, xen, yen, xk, yk, xke, yke, Jp      & ! Output
         ,  Je, Jn, alphae, gamman, betae, betan, radius, re, rn, rp  ) ! Output



   ! Initialize the flags of boundary faces
   call initialize_boundary_faces_flags(nx, ny, fbe, fbn) ! Output: last 2



   ! Guess the initial conditions
   call get_initial_conditions( nx, ny, itemax, modvis, PF, TF, Rg, GF, MF, UF &
      , xe, ye, xk, yk, xke, yke, alphae, betae, p, T, ro, roe, ron, u, v, ue  &
      , ve, un, vn, Uce, Vcn, Ucbe, Vcbe, Tbn, Tbs, Tbe, Tbw ) ! UF and last 19 are output
   ! checked


   ! Initializes the vectors with specific heat Cp, viscosity mu and thermal
   ! conductivity kappa.
   call get_cp_mu_kappa_initialization(nx, ny, ktm, modvis, CPF, VLF, KPF &
      ,              Tbn, Tbs, Tbe, Tbw, T, cp, vlp, vle, vln, kp, ke, kn ) ! Output: last seven


   ! Time step initialization
   call tstepper_init(ifile)

   ! Calculates the first time step
   dt = tstepper_dt0(0, nx, ny, Jp, u, v)

   ! Initialize solvers
   call solvers_init(nx, ny, ifile, solver5d, solver9d)

   ! Selects the portion of the forebody where the drag will be calculated
   select case ( kfc )

      case (0) ! Over the whole forebody

         iocs = nx-1

      case (1) ! Over the ogive only

         call get_ogive_cylinder_matching_point_south(nx, ny, xk, yk, iocs) ! Output: last one

      case default

         write(*,*) "kfc: unknown option. Stopping..."

         stop

   end select


   ! Beginning time cycle

   ! Openning file of residuals
   open(rid, file = './mach2d_output/mach2d_' // trim(adjustl(sim_id)) // '_residual.dat' )


   ! CPU time (first measurement)

   call cpu_time(tcpu1)

   it_stop = -1

   do it = 1, itmax

      roa = ro
      Ta  = T
      pa  = p
      ua  = u
      va  = v
      uea = ue
      vea = ve
      una = un
      vna = vn


      ! Calculations of the thermophysical properties

      if ( ktm == THERMOPHYSICAL_VARIABLE ) then ! Temperature dependent thermophysical properties


         ! Calculates the temperature over the boundaries
         call get_north_boundary_field(nx, ny, T, Tbn) ! Output: last one
         call get_south_boundary_field(nx, ny, T, Tbs) ! Output: last one
         call get_east_boundary_field(nx, ny, T, Tbe) ! Output: last one
         call get_west_boundary_field(nx, ny, T, Tbw) ! Output: last one


         ! Calculates cp at the center of each real volume and over the boundaries
         call set_cp(nx, ny, Tbn, Tbs, Tbe, Tbw, T, cp) ! Output: last one


         ! Thermophysical properties for the Navier-Stokes equations only

         if ( modvis == 1 ) then

            ! Calculates the laminar viscosity at the nodes of real volumes
            call set_laminar_viscosity_at_nodes(nx, ny, T, vlp) ! Output: last one


            ! Calculates the thermal conductivity at the nodes of real volumes
            call set_thermal_conductivity_at_nodes(nx, ny, T, kp) ! Output: last one


            ! Calculates the laminar viscosity at faces
            call get_laminar_viscosity_at_faces(nx, ny, Tbn, Tbs, Tbe, Tbw, vlp, vle, vln) ! Output: last two


            ! Calculates the thermal conductivity at faces
            call get_thermal_conductivity_at_faces(nx, ny, Tbn, Tbs, Tbe, Tbw, kp, ke, kn) ! Output: last two

         end if

      end if


      ! Calculates the average value of the pressure
      call get_average(nx, ny, p, p_avg) ! Output: last one


      ! Calculates g

      g = 1.d0/(Rg*T)

      ! Beginning the mass correction cycle

      do itm = 1, itmmax

         norm = 0.d0


         ! Calculates the contravariant velocities over the east boundary
         call get_Ucbe_Vcbe(nx, ny, xe, ye, xke, yke, u, v, Ucbe, Vcbe) ! Output: last two



         ! Coefficients of the linear system for u (real volumes)
         call get_u_coefficients( nx, ny, modvis, dt, rp, re, rn, Jp, Je, Jn &
            ,                        ye, yk, alphae, betae, betan, gamman       &
            ,                        vle, vln, roe, ron, roa, Uce, Vcn, au ) ! au is output



         ! Source of the linear system for u (real volumes)
         call get_u_source( nx, ny, modvis, beta, dt, rp, re, rn   &
            ,                  xe, ye, xk, yk, xke, yke, xen, yen     &
            ,                  Jp, Je, Jn, roe, ron, roa, p, vle, vln &
            ,                  Uce, Vcn, ua, u, v, cup, sup, bu ) ! The last 3 are output


         ! ========================================================================

         !                      USER DEPENDENT SUBROUTINE

         call get_bc_scheme_u(nx, ny, modvis, UF, xk, yk, alphae, betae, u, v &
            , Ucbe, Vcbe, a9bn, a9bs, a9be, a9bw, b9bn, b9bs, b9be, b9bw) ! Output: last eight

         ! ========================================================================


         ! Transfers the numerical scheme of the boundary conditions to the linear
         ! system coefficients and source
         call get_bc_transfer_9d(nx, ny, a9bn, a9bs, a9be, a9bw, b9bn, b9bs &
            , b9be, b9bw, au, bu) ! Output: last two


         ! Coefficients of the linear system for v (real volumes)
         call get_v_coefficients( nx, ny, coord, modvis, dt, rp, re, rn, Jp, Je, Jn &
            ,                        xe, xk, alphae, betae, betan, gamman       &
            ,                        vle, vln, vlp, roe, ron, roa, Uce, Vcn, av ) ! Last 1 is output


         ! Source of the linear system for v (real volumes)
         call get_v_source( nx, ny, coord, modvis, beta, dt, rp, re, rn   &
            ,                  xe, ye, xk, yk, xke, yke, xen, yen     &
            ,                  Jp, Je, Jn, roe, ron, roa, p, vle, vln &
            ,                  Uce, Vcn, va, u, v, cvp, svp, bv ) ! Last 3 are output


         ! ========================================================================

         !                      USER DEPENDENT SUBROUTINE

         call get_bc_scheme_v(nx, ny, modvis, xk, yk, u, v, Ucbe, Vcbe &
            , a9bn, a9bs, a9be, a9bw, b9bn, b9bs, b9be, b9bw) ! Output: last eight

         ! ========================================================================


         ! Transfers the numerical scheme of the boundary conditions to the linear
         ! system coefficients and source
         call get_bc_transfer_9d(nx, ny, a9bn, a9bs, a9be, a9bw, b9bn, b9bs &
            , b9be, b9bw, av, bv) ! Output: last two



         ! Calculates SIMPLEC coefficients at the internal real faces
         call get_internal_simplec_coefficients( nx, ny, re, rn, xe, ye, xk &
            , yk, au, av, due, dve, dun, dvn, de, dn ) ! Output: last 6



         ! Calculates the SIMPLEC coefficients at the boundary faces
         call get_boundary_simplec_coefficients( nx, ny, modvis, due, dve, dun &
            , dvn, de, dn) ! InOutput: last 6


         ! Solves the linear system for u
         call solver9d%solve(au, bu, u)

         ! Calculates the norm of the residual of the linear system for u
         call norm_l1_9d_relative( nx, ny, u, bu, au, norm) ! InOutput: last one

         ! Solves the linear system for v
         call solver9d%solve(av, bv, v)

         ! Calculates the norm of the residual of the linear system for v
         call norm_l1_9d_relative( nx, ny, v, bv, av, norm) ! InOutput: last one


         ! Uses velocities at nodes to calculate velocities at internal real faces
         call get_velocities_at_internal_faces(nx, ny, dt, rp, re, rn, xe, ye, xk, yk & ! Input
         ,                             xen, yen, xke, yke, Jp, cup, cvp, sup, svp     & ! Input
         ,                             au, av, roa, p, u, v, uea, vea, una, vna       & ! Input
         ,                             ue, ve, un, vn, Uce, Vcn)                        ! Output


         ! Calculates the velocities at boundary faces
         call get_velocities_at_boundary_faces( nx, ny, xe, ye, xk, yk, u, v, ue, ve, un, vn, Uce, Vcn)
         ! Last six are output


         ! Initializing the pressure deviation

         pl = 0.d0

         do itp = 1, itpmax


            ! Calculates the coefficients of the linear system for pressure correction
            call get_p_coefficients(nx, ny, dt, rp, re, rn, Jp, Uce, Vcn &
               ,                                  roe, ron, g, de, dn, ap) ! Output: last one

            ! ========================================================================

            !                      USER DEPENDENT SUBROUTINE

            ! Defines the numerical scheme for the boundary conditions of pl
            call get_bc_scheme_pl(nx, ny, a5bn, a5bs, a5be, a5bw, b5bn, b5bs, b5be, b5bw) ! Output: last eight

            ! ========================================================================

            ! Transfers the numerical scheme of the boundary conditions to the linear
            ! system coefficients and source
            call get_bc_transfer_5d(nx, ny, a5bn, a5bs, a5be, a5bw, b5bn, b5bs &
               , b5be, b5bw, ap, bp) ! Output: last two


            ! Calculates the source of the linear system of the pressure correction
            call get_p_source(nx, ny, dt, rp, re, rn, Jp, roe, ron, ro, roa &
            ,                                    Uce, Vcn, de, dn, g, pl, bp) ! Output: last one


            ! Solves the linear system for pl
            call solver5d%solve(ap, bp, pl)

            ! Velocity correction at faces (SIMPLEC method)
            call get_velocity_correction_at_faces_with_pl(nx, ny, xe, ye, xk, yk  &
               , due, dve, dun, dvn, de, dn, pl, ue, ve, un, vn, Uce, Vcn ) ! InOutput: last six

            ! Density correction
            call get_density_correction_with_pl( nx, ny, pl, g, ro) ! InOutput: last one

            ! Calculates density at faces using the corrected density and velocities
            call get_density_at_faces( nx, ny, beta, ro, Uce, Vcn, roe, ron) ! roe and ron are output

            normpl = maxval(abs(pl)) / p_avg

            if ( normpl < tolm ) exit

         end do


         ! Calculates the norm of the residuals
         call norm_l1_5d( nx, ny, pl, bp, ap, norm)


         ! Pressure correction  (SIMPLEC method)
         call get_pressure_correction_with_pl( nx, ny, pl, p) ! InOutput: last two



         ! ========================================================================

         !                      USER DEPENDENT SUBROUTINE

         ! Defines the numerical scheme for the boundary conditions of ro
         call get_bc_scheme_ro(nx, ny, ROF, alphae, betae, betan, gamman &
            , Ucbe, Vcbe, a9bn, a9bs, a9be, a9bw, b9bn, b9bs, b9be, b9bw) ! Output: last eight

         ! ========================================================================

         ! Extrapolates ro from the real volumes to the fictitious ones according to the boundary conditions
         call get_bc_extrapolation_9d(nx, ny, itemax, a9bn, a9bs, a9be, a9bw, b9bn, b9bs &
            , b9be, b9bw, ro) ! InOutput: last one




         ! Calculates the values of p at fictitious volumes
         call get_p_extrapolation_to_fictitious(nx, ny, p) ! InOutput: last one



         ! Velocity correction (SIMPLEC method)
         call get_u_v_at_real_nodes_with_pl( nx, ny, xe, ye, xk & ! Input
         ,                                    yk, rp, pl, au, av & ! Input
         ,                                    u, v )               ! Input and Output


         ! ========================================================================

         !                      USER DEPENDENT SUBROUTINE

         call get_bc_scheme_u(nx, ny, modvis, UF, xk, yk, alphae, betae, u, v &
            , Ucbe, Vcbe, a9bn, a9bs, a9be, a9bw, b9bn, b9bs, b9be, b9bw) ! Output: last eight

         ! ========================================================================


         ! Extrapolates u from the real volumes to the fictitious ones according to the boundary conditions
         call get_bc_extrapolation_9d(nx, ny, itemax, a9bn, a9bs, a9be, a9bw, b9bn, b9bs &
            , b9be, b9bw, u) ! InOutput: last one

         ! ========================================================================

         !                      USER DEPENDENT SUBROUTINE

         call get_bc_scheme_v(nx, ny, modvis, xk, yk, u, v, Ucbe, Vcbe &
            , a9bn, a9bs, a9be, a9bw, b9bn, b9bs, b9be, b9bw) ! Output: last eight

         ! ========================================================================


         ! Extrapolates v from the real volumes to the fictitious ones according to the boundary conditions
         call get_bc_extrapolation_9d(nx, ny, itemax, a9bn, a9bs, a9be, a9bw, b9bn, b9bs &
            , b9be, b9bw, v) ! InOutput: last one


      end do


      ! Checkes the mass conservation
      call check_mass_conservation(nx, ny, dt, rp, re, rn, Jp, ro, roa &
         , ron, roe, Uce, Vcn, rmass) ! Output: last one




      ! Selecting between Euler and Navier-Stokes equations


      if ( modvis == 0 .and. ktm == THERMOPHYSICAL_CONSTANT ) then ! Euler with constant thermophysical properties


         ! Calculates the temperature based on the conservation of the total enthalpy
         ! Valid for Euler model with constant thermo-physical coefficients.
         ! Temperature is extrapolated to fictitious volumes with CDS the scheme
         call get_T_from_H_conservation(nx, ny, CPF, HF, u, ue, un, v, ve, vn, T, Tbe, Tbw, Tbn, Tbs) ! Output: last 5


      else ! Variable thermophysical properties


         call get_T_coefficients_and_source( nx, ny, modvis, beta, dt, rp, re, rn                &
         ,                  xe, ye, xk, yk, xke, yke, xen, yen, alphae, betae, betan, gamman     &
         ,                  Jp, Je, Jn, roe, ron, ro, roa, p, pa, cp, vle, vln, ke, kn, Uce, Vcn &
         ,                  u, v, ue, ve, un, vn, ua, va, T, Ta, fbe, fbn, at, bt ) ! Last 2 are output


         ! ========================================================================

         !                      USER DEPENDENT SUBROUTINE

         ! Defines the numerical scheme for the boundary conditions of T
         call get_bc_scheme_T(nx, ny, TF, Tsbc, alphae, betae, betan, gamman &
            , Ucbe, Vcbe, a9bn, a9bs, a9be, a9bw, b9bn, b9bs, b9be, b9bw) ! Output: last eight


         ! ========================================================================


         ! Transfers the numerical scheme of the boundary conditions to the linear
         ! system coefficients and source
         call get_bc_transfer_9d(nx, ny, a9bn, a9bs, a9be, a9bw, b9bn, b9bs &
            , b9be, b9bw, at, bt) ! Output: last two


         ! Solves the linear system for T
         call solver9d%solve(at, bt, T)

         ! Calculates the norm of the residual of the linear system for T
         call norm_l1_9d_relative( nx, ny, T, bt, at, norm) ! InOutput: last one



      end if


      ! Calculating p from the state equation

      p = ro * Rg * T


      ! Calculates the values of p at fictitious volumes
      call get_p_extrapolation_to_fictitious(nx, ny, p) ! InOutput: last one


      ! Selecting the next time step
      dt = tstepper_dt(it, nx, ny, dt, au, av, at, ap)


      if ( it == 1 ) then

         call get_cdfi(nx, ny, 2, iocs, Rg, PF, TF, UF, rb, yk, rn, p, Cdfi)

         if ( modvis == 1 ) then

            call get_cdfv(nx, ny, 2, iocs, Rg, PF, TF, UF, rb, xk, yk, rn, Jn, vln &
               , u, v, Vcn, Cdfv )

         else

            Cdfv = 0.d0

         end if

         write(  *,"(A10, 12(1X, A15))") 'it', 'norm', 'max(|pl|)/p_avg'  &
            , 'Cdfi', 'Cdfv', 'dt' &
            , 'rmass'

         write(rid,"(A10, 12(1X, A23))") 'it', 'norm', 'max(|pl|)/p_avg'  &
            , 'Cdfi', 'Cdfv', 'dt' &
            , 'rmass'

         write(rid,"(I10, 12(1X,ES23.16))") it, norm, normpl    &
            , Cdfi, Cdfv, dt, rmass

      end if


      ! Checking for wrong results

      if ( Is_NaN( norm ) ) then

         write(*,*) "NaN detected. Stopping..."

         write(lid,*) "NaN detected. Stopping..."

         ! Saves the fitness function to a file
         ! Exit status: 0 = success, 1 = failure, 2 = generate another individual
         call depp_save_fitness(-1.d10, 1, "Fitness: -(Cdfi+Cdfv)")

         stop

      end if

      ! Printing data

      if ( mod( it, wlf ) == 0 ) then

         call get_cdfi(nx, ny, 2, iocs, Rg, PF, TF, UF, rb, yk, rn, p, Cdfi)

         if ( modvis == 1 ) then

            call get_cdfv(nx, ny, 2, iocs, Rg, PF, TF, UF, rb, xk, yk, rn, Jn, vln &
               , u, v, Vcn, Cdfv )

         else

            Cdfv = 0.d0

         end if

         write(  *,"(I10, 12(1X,ES15.8))") it, norm, normpl     &
            , Cdfi, Cdfv, dt, rmass

         write(rid,"(I10, 12(1X,ES23.16))") it, norm, normpl    &
            , Cdfi, Cdfv, dt, rmass

      end if

      ! Checking convergence criteria

      if ( normpl < tolt .and. it_stop < 0 ) it_stop = it * 2

      if ( it_stop == it ) exit

   end do

   ! CPU time (second measurement)

   call cpu_time(tcpu2)

   ! Total CPU time

   tcpu = tcpu2 - tcpu1

   ! Gets the current process memory in MB

   RAM = 0.d0 !get_RAM(sim_id)

   ! Closing file of residuals
   close(rid)


   ! Recalculating the foredrag coefficients before finishing
   call get_cdfi(nx, ny, 2, iocs, Rg, PF, TF, UF, rb, yk, rn, p, Cdfi)

   if ( modvis == 1 ) then

      call get_cdfv(nx, ny, 2, iocs, Rg, PF, TF, UF, rb, xk, yk, rn, Jn, vln &
         , u, v, Vcn, Cdfv )

   else

      Cdfv = 0.d0

   end if


   ! Saves the fitness function to a file
   ! Exit status: 0 = success, 1 = failure, 2 = generate another individual
   if ( it < itmax ) then

      call depp_save_fitness(-(Cdfi+Cdfv), 0, "Fitness: -(Cdfi+Cdfv)")

   else

      call depp_save_fitness(-1.d5, 1, "Fitness: -(Cdfi+Cdfv) (no convergence)")

   end if


   ! Saves post processed data

   if ( wppd == 1 ) then

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


      ! Post processing data

      call post_processing

   ! Simplified post-processing
   else if ( wppd == 2 ) then

      call simplified_post_processing

   end if

contains

   subroutine post_processing

      implicit none

      ! Auxiliary variables
      real(8), dimension (nx*ny) :: M    ! Mach number at the center of the volumes

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
      call set_gamma(nx, ny, Rg, cp, gcp) ! Output: last one

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

      implicit none

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
end program
