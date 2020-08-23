!>
!! \brief      MACH-2D
!!                 Version: 5.10.1
!!             Last update: 16/08/2020 by Guilherme Bertoldo
!!
program main

   ! Importing modules
   use mod_mach2d_data
   use mod_grid
   use mod_solvers
   use mod_tstepper
   use mod_error_handler
   use mod_thermophysical
   use mod_field_operations
   use depp_interface
   use coefficients
   use mod_extflow
   use mod_extflow_postp
   use mod_intflow_postp

   implicit none

   ! Reads parameters, allocates variables and sets vectors to zero
   call initialize_variables


   ! Initializes variables related to the thermophysical properties of the gas
   call thermophysical_initialization(ifile, thermomodel, Rg, ktm) ! Output: last three


   ! Initializing variables related to internal or external flow
   if ( kflow == EXTERNAL_FLOW ) then
      ! Initializing external flow module
      call extflow_init(ifile, thermomodel, Tref, Href, Cpref) ! Output: last three
   else
      ! Initializing reference values and internal flow variables
      call intflow_init(ifile, thermomodel, Tref, Href, Cpref, iflow) ! Output: last four
   end if


   ! Creating the grid and calculating the metrics
   call grid_create(lid, coord, kflow, gridpar                          & ! Input
         ,  x, y, xp, yp, xe, ye, xen, yen, xk, yk, xke, yke, Jp      & ! Output
         ,  Je, Jn, alphae, gamman, betae, betan, radius, re, rn, rp, iflow  ) ! Output


   ! Initialize the flags of boundary faces
   call initialize_boundary_faces_flags(nx, ny, fbe, fbn) ! Output: last 2


   ! Selecting initial condition in accordance to the kind of flow
   if ( kflow == EXTERNAL_FLOW ) then
      ! Sets the initial conditions for external flow
      call extflow_initial_conditions(nx, ny, itemax, modvis, Rg, xe, ye & ! Input
         ,                               xk, yk, xke, yke, alphae, betae & ! Input
         ,                      p, T, ro, roe, ron, u, v, ue, ve, un, vn & ! Output
         ,                                  Uce, Vcn, Tbn, Tbs, Tbe, Tbw ) ! Output
   else
      !> \brief Calculates initial conditions for internal flow
      call intflow_initial_conditions( nx, ny, modvis, beta & ! Input
         ,thermomodel, xe, ye, xk, yk, radius, rn, x, y, xp & ! Input
         ,                                            iflow & ! InOutput
         ,             p, T, u, v, ue, un, ve, vn, Uce, Vcn & ! Output
         ,                 ro, roe, ron, Tbn, Tbs, Tbe, Tbw ) ! Output

   end if

   ! Initializes the vectors with specific heat Cp, viscosity mu
   ! and thermal conductivity kappa based on the reference temperature Tref.
   call get_cp_mu_kappa_initialization(nx, ny, ktm, modvis, thermomodel & ! Input
         ,                                  Tref, Tbn, Tbs, Tbe, Tbw, T & ! Input
         ,                                cp, vlp, vle, vln, kp, ke, kn ) ! Output


   ! Time step initialization
   call tstepper_init(ifile)

   ! Calculates the first time step
   dt = tstepper_dt0(0, nx, ny, Jp, u, v)


   ! Initialize solvers
   call solvers_init(nx, ny, ifile, solver5d, solver9d)


   ! CPU time (first measurement)
   call cpu_time(tcpu1)


   ! Beginning time cycle
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

      if ( kflow == INTERNAL_FLOW ) then
         ! Update inlet conditions and extrapolates pressure to fictitious
         ! according to boundary conditions.
         call intflow_time_cycle_update(nx, ny, x, xp, u, thermomodel, iflow, p) ! InOutput: last two
      end if

      ! Calculations of the thermophysical properties
      if ( ktm == THERMOPHYSICAL_VARIABLE ) then ! Temperature dependent thermophysical properties


         ! Calculates the temperature over the boundaries
         call get_north_boundary_field(nx, ny, T, Tbn) ! Output: last one
         call get_south_boundary_field(nx, ny, T, Tbs) ! Output: last one
         call get_east_boundary_field(nx, ny, T, Tbe) ! Output: last one
         call get_west_boundary_field(nx, ny, T, Tbw) ! Output: last one


         ! Calculates cp at the center of each real volume and over the boundaries
         call set_cp(thermomodel, nx, ny, Tbn, Tbs, Tbe, Tbw, T, cp) ! Output: last one


         ! Thermophysical properties for the Navier-Stokes equations only

         if ( modvis == 1 ) then

            ! Calculates the laminar viscosity at the nodes of real volumes
            call set_laminar_viscosity_at_nodes(thermomodel, nx, ny, T, vlp) ! Output: last one


            ! Calculates the thermal conductivity at the nodes of real volumes
            call set_thermal_conductivity_at_nodes(thermomodel, nx, ny, T, kp) ! Output: last one


            ! Calculates the laminar viscosity at faces
            call get_laminar_viscosity_at_faces(thermomodel, nx, ny, Tbn, Tbs, Tbe, Tbw, vlp, vle, vln) ! Output: last two


            ! Calculates the thermal conductivity at faces
            call get_thermal_conductivity_at_faces(thermomodel, nx, ny, Tbn, Tbs, Tbe, Tbw, kp, ke, kn) ! Output: last two

         end if

      end if


      ! Calculates the average value of the pressure
      call get_average(nx, ny, p, p_avg) ! Output: last one


      ! Calculates g
      g = 1.d0/(Rg*T)


      ! Beginning the mass correction cycle
      do itm = 1, itmmax

         ! Norm initialization
         norm = 0.d0


         ! Coefficients of the linear system for u (real volumes)
         call get_u_coefficients( nx, ny, modvis, dt, rp, re, rn, Jp, Je, Jn &
            ,                           ye, yk, alphae, betae, betan, gamman &
            ,                          vle, vln, roe, ron, roa, Uce, Vcn, au ) ! au is output



         ! Source of the linear system for u (real volumes)
         call get_u_source( nx, ny, modvis, beta, dt, rp, re, rn &
            ,                 xe, ye, xk, yk, xke, yke, xen, yen &
            ,             Jp, Je, Jn, roe, ron, roa, p, vle, vln &
            ,                   Uce, Vcn, ua, u, v, cup, sup, bu ) ! The last 3 are output


         ! Boundary conditions for u
         if ( kflow == EXTERNAL_FLOW ) then
            call extflow_set_bcu(nx, ny, modvis, xk, yk, alphae, betae, u, v & ! Intput
               ,                                                      au, bu ) ! Output
         else
            call intflow_set_bcu(nx, ny, modvis, x, xp, xk, yk, u, v, au, bu) ! Output: last two
         end if


         ! Coefficients of the linear system for v (real volumes)
         call get_v_coefficients( nx, ny, coord, modvis, dt, rp, re, rn, Jp &
            ,                  Je, Jn, xe, xk, alphae, betae, betan, gamman &
            ,                    vle, vln, vlp, roe, ron, roa, Uce, Vcn, av ) ! Last 1 is output


         ! Source of the linear system for v (real volumes)
         call get_v_source( nx, ny, coord, modvis, beta, dt, rp, re, rn &
            ,                        xe, ye, xk, yk, xke, yke, xen, yen &
            ,                    Jp, Je, Jn, roe, ron, roa, p, vle, vln &
            ,                          Uce, Vcn, va, u, v, cvp, svp, bv ) ! Last 3 are output


         ! Boundary conditions for v
         if ( kflow == EXTERNAL_FLOW ) then
            call extflow_set_bcv(nx, ny, modvis, xk, yk, u, v, av, bv ) ! Output: last two
         else
            call intflow_set_bcv(nx, ny, modvis, x, xp, xk, yk, u, v, av, bv) ! Output: last two
         end if


         ! Calculates SIMPLEC coefficients at the internal real faces
         call get_internal_simplec_coefficients( nx, ny, re, rn, xe, ye, xk &
            , yk, au, av, due, dve, dun, dvn, de, dn ) ! Output: last 6


         ! Calculates the SIMPLEC coefficients at the boundary faces
         if ( kflow == EXTERNAL_FLOW ) then
            call extflow_boundary_simplec( nx, ny, modvis & ! Input
               ,               due, dve, dun, dvn, de, dn ) ! InOutput
         else
            call intflow_boundary_simplec( nx, ny, modvis & ! Input
               ,               due, dve, dun, dvn, de, dn ) ! InOutput
         end if


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
         if ( kflow == EXTERNAL_FLOW ) then
            call extflow_velocities_at_boundary_faces( nx, ny, xe, ye, xk, yk, u, v & ! Input
               ,                                           ue, ve, un, vn, Uce, Vcn ) ! Output
         else
            call intflow_velocities_at_boundary_faces( nx, ny, modvis & ! Input
               ,                          x, xp, xe, ye, xk, yk, u, v & ! Input
               ,                             ue, ve, un, vn, Uce, Vcn ) ! Input and Output
         end if


         ! Initializing the pressure deviation

         pl = 0.d0

         do itp = 1, itpmax


            ! Calculates the coefficients of the linear system for pressure correction
            call get_p_coefficients(nx, ny, dt, rp, re, rn, Jp, Uce, Vcn &
               ,                                 roe, ron, g, de, dn, ap ) ! Output: last one

            ! Defines the numerical scheme for the boundary conditions of pl
            if ( kflow == EXTERNAL_FLOW ) then
               call extflow_set_bcpl(nx, ny, ap, bp) ! Output: last eight
            else
               call intflow_set_bcpl(nx, ny, x, xp, pl, iflow, ap, bp) ! Output: last two
            end if

            ! Calculates the source of the linear system of the pressure correction
            call get_p_source(nx, ny, dt, rp, re, rn, Jp, roe, ron, ro, roa &
            ,                                   Uce, Vcn, de, dn, g, pl, bp ) ! Output: last one


            ! Solves the linear system for pl
            call solver5d%solve(ap, bp, pl)


            ! Velocity correction at faces (SIMPLEC method)
            call get_velocity_correction_at_faces_with_pl(nx, ny, xe, ye, xk, yk &
               ,        due, dve, dun, dvn, de, dn, pl, ue, ve, un, vn, Uce, Vcn ) ! InOutput: last six

            ! Correcting boundary velocities and pressure field
            if ( kflow == INTERNAL_FLOW ) then
               call intflow_velocities_at_boundary_faces( nx, ny, modvis & ! Input
                  ,                          x, xp, xe, ye, xk, yk, u, v & ! Input
                  ,                             ue, ve, un, vn, Uce, Vcn ) ! Input and Output

               ! Pressure correction  (SIMPLEC method)
               call get_pressure_correction_with_pl( nx, ny, pl, p) ! InOutput: last two
            end if

            ! Density correction
            call get_density_correction_with_pl( nx, ny, pl, g, ro) ! InOutput: last one



            ! Calculates density at faces using the corrected density and velocities
            call get_density_at_faces( nx, ny, beta, ro, Uce, Vcn, roe, ron) ! roe and ron are output


            normpl = maxval(abs(pl)) / p_avg

            if ( normpl < tolm ) exit

         end do


         ! Calculates the norm of the residuals
         call norm_l1_5d( nx, ny, pl, bp, ap, norm)


         ! Extrapolates p and to to fictitious
         if ( kflow == EXTERNAL_FLOW ) then

            ! Pressure correction  (SIMPLEC method)
            call get_pressure_correction_with_pl( nx, ny, pl, p) ! InOutput: last two

            ! Extrapolates ro to fictitious volumes in accordance to bc
            call extflow_extrapolate_ro_to_fictitious(nx, ny, itemax, alphae &
               ,                                    betae, betan, gamman, ro ) ! InOutput: last one

            ! Calculates the values of p at fictitious volumes
            call extflow_p_extrapolation_to_fictitious(nx, ny, p) ! InOutput: last one
         end if


         ! Velocity correction (SIMPLEC method)
         call get_u_v_at_real_nodes_with_pl( nx, ny, xe, ye, xk & ! Input
            ,                                yk, rp, pl, au, av & ! Input
            ,                                              u, v ) ! Input and Output

         ! Extrapolates u and v to fictitious
         if ( kflow == EXTERNAL_FLOW ) then
            ! Extrapolates u and v to fictitious volumes in accordance to bc
            call extflow_extrapolate_u_v_to_fictitious(nx, ny, modvis, itemax &
               ,                                  xk, yk, alphae, betae, u, v ) ! InOutput: last two
         else
            call intflow_extrapolate_u_v_to_fictitious( nx, ny, modvis &
               ,                                   x, xp, xk, yk, u, v ) ! InOutput: last two
         end if

      end do


      ! Checkes the mass conservation
      call check_mass_conservation(nx, ny, dt, rp, re, rn, Jp, ro, roa &
         ,                                   ron, roe, Uce, Vcn, rmass ) ! Output: last one


      ! Calculation of T based on total enthalpy conservation is enabled only
      ! for Euler model with constant thermophysical properties
      if (      modvis == 0                       &
         .and.     ktm == THERMOPHYSICAL_CONSTANT ) then ! Euler with constant thermophysical properties

         ! Temperature at the centroid of all nodes (real and fictitious)
         T = ( Href - (u**2+v**2) / 2.d0 ) / Cpref

         ! Calculates the boundary temperature based on the  boundary conditions
         ! and on conservation of the total enthalpy (valid for Euler model with constant
         ! thermo-physical coefficients). Temperature is extrapolated to fictitious
         ! volumes using the CDS scheme.
         if ( kflow == EXTERNAL_FLOW  ) then
            call extflow_boundary_T_from_H_conservation(nx, ny, Cpref, Href, u &
               ,                      ue, un, v, ve, vn, T, Tbe, Tbw, Tbn, Tbs ) ! Output: last 5. T is inout
         else
            call intflow_boundary_T_from_H_conservation(nx, ny, Cpref, Href, iflow &
               ,                       u, ue, un, v, ve, vn, T, Tbe, Tbw, Tbn, Tbs ) ! Output: last 5. T is inout
         end if

      else ! Variable thermophysical properties


         call get_T_coefficients_and_source( nx, ny, modvis, beta, dt, rp, re, rn &
         ,       xe, ye, xk, yk, xke, yke, xen, yen, alphae, betae, betan, gamman &
         ,   Jp, Je, Jn, roe, ron, ro, roa, p, pa, cp, vle, vln, ke, kn, Uce, Vcn &
         ,                  u, v, ue, ve, un, vn, ua, va, T, Ta, fbe, fbn, at, bt ) ! Last 2 are output


         ! Defines the numerical scheme for the boundary conditions of T
         if ( kflow == EXTERNAL_FLOW ) then
            call extflow_set_bcT(nx, ny, alphae, betae, betan, gamman, at, bt ) ! Output: last two
         else
            call intflow_set_bcT(nx, ny, x, xp, T, iflow, at, bt) ! Output: last two
         end if


         ! Solves the linear system for T
         call solver9d%solve(at, bt, T)

         ! Calculates the norm of the residual of the linear system for T
         call norm_l1_9d_relative( nx, ny, T, bt, at, norm) ! InOutput: last one

      end if

      if ( kflow == EXTERNAL_FLOW ) then
         ! Calculating p from the state equation
         p = ro * Rg * T

         ! Calculates the values of p at fictitious volumes
         call extflow_p_extrapolation_to_fictitious(nx, ny, p) ! InOutput: last one
      else
         ! Calculating ro at all nodes according to state equation
         call get_density_at_nodes( nx, ny, Rg, p, T, ro) ! ro is output

         ! Calculating ro at all faces by interpolation
         call get_density_at_faces( nx, ny, beta, ro, Uce, Vcn, roe, ron) ! roe and ron are output
      end if

      ! Selecting the next time step
      dt = tstepper_dt(it, nx, ny, dt, au, av, at, ap)


      ! Calculates and prints main variables
      call calc_print_main_variables(.false.)


      ! Checking for wrong results
      if ( Is_NaN( norm ) ) then

         write(*,*) "NaN detected. Stopping..."

         write(lid,*) "NaN detected. Stopping..."

         ! Saves the fitness function to a file
         ! Exit status: 0 = success, 1 = failure, 2 = generate another individual
         call depp_save_fitness(-1.d10, 1, "Fitness: failure")

         stop

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


   ! Calculates and prints main variables
   call calc_print_main_variables(.true.)


   ! Performing post-processing
   if ( kflow == EXTERNAL_FLOW ) then
      call extflow_postp()
   else
      call intflow_postp()
   end if

   ! Closing file of residuals
   close(rid)

contains

   !> \brief Calculates and prints the main variables
   subroutine calc_print_main_variables(hasFinished)
      implicit none
      logical, intent(in) :: hasFinished !< Main loop has finished? (true, false)

      if ( kflow == EXTERNAL_FLOW ) then

         ! Calculates the main variables of internal flow
         call extflow_calc_main_variables(nx, ny, Rg, modvis, xk, yk, rn, Jn, vln, p, u, v, Vcn, msg) ! Output: last one

      else ! INTERNAL_FLOW

         ! Calculates the main variables of internal flow
         call intflow_calc_main_variables(nx, ny, re, roe, u, Uce, msg)

      end if


     ! Printing header
     if ( it == 1 ) then

         write(  *,"(A10, 4(1X, A14), A)") 'it'               &
            ,                              'dt'               &
            ,                              'norm'             &
            ,                              'max(|pl|)/p_avg'  &
            ,                              'rmass'            &
            ,                              msg(1)

         write(rid,"(A10, 4(1X, A14), A)") 'it'               &
            ,                              'dt'               &
            ,                              'norm'             &
            ,                              'max(|pl|)/p_avg'  &
            ,                              'rmass'            &
            ,                              msg(1)

      end if

      ! Printing data
      if ( mod( it, wlf ) == 0 .or. hasFinished ) then

         write(  *,"(I10, 4(1X,ES14.7), A)") it, dt, norm, normpl, rmass, trim(msg(2))
         write(rid,"(I10, 4(1X,ES14.7), A)") it, dt, norm, normpl, rmass, trim(msg(2))

      end if

   end subroutine


end program
