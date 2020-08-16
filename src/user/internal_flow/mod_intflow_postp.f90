module mod_intflow_postp

    use mod_intflow_procedures

    !
    ! SUBROUTINES
    !
    !  0) intflow_postp
    !  1) post_processing
    !  2) plotter01
    !  3) plotter02
    !  4) plotter03
    !  5) plotter04
    !  6) plotter05
    !  7) plotter06
    !  8) write_nozzle_parameters
    !  9) write_grid_parameters
    ! 10) write_grid_boundary
    ! 11) write_grid
    ! 12) write_main_results
    ! 13) write_field
    ! 14) write_main_fields
    ! 15) post_processing_boundaries
    ! 16) get_stream_function
    !
    ! Last update: 16 Aug 2020
    !
    character (len = *), parameter, private :: text_editor     = "emacs -nw "
    character (len = *), parameter, private :: graph_viewer    = "evince "
    character (len = *), parameter, private :: graph_generator = "gnuplot "

contains

   !> \brief This is an interface procedure to communicate old post-processing
   !! subroutines of Mach2D5.8.1 to the current version of the code. Some data
   !! defined in this procedure must be better integrated to the main code.
   subroutine intflow_postp()
      use mod_mach2d_data
      implicit none

      integer :: ig
      integer :: modtur
      integer :: ccTw
      real(8) :: gamma
      real(8) :: tcpuo
      real(8) :: fmi
      real(8) :: fme
      real(8) :: Fd
      real(8) :: pr
      real(8) :: go
      real(8) :: Twall(nx)
      real(8) :: vtp(nx*ny)

      ! i=ig when the throat cross section is in the east face of the CV
      ig = get_ig(nx, ny, y)

      modtur = 0
      pr     = 1.013d5
      go     = 9.81d0
      tcpuo  = 0.d0
      gamma  = thermomodel%gm(Tref)
      Twall  = iflow%Twall
      vtp    = 0.d0
      ccTw   = 0
      if ( iflow%Twall > 0.d0 ) ccTw = 1

      call get_mass_flow_rate_and_thrust( nx, ny, re, roe, u, Uce, fmi, fme, Fd)

      associate(    uin => iflow%uin       &
            ,       vin => iflow%vin       &
            ,       pin => iflow%pin       &
            ,       Tin => iflow%Tin       &
            ,       u1D => iflow%u1D       &
            ,       p1D => iflow%p1D       &
            ,       T1D => iflow%T1D       &
            ,       M1D => iflow%M1D       &
            ,     Fpv1D => iflow%Fpv1D     &
            ,      Fd1D => iflow%Fd1D      &
            ,      fm1D => iflow%fm1D      &
            ,        po => iflow%P0        &
            ,        Sg => iflow%geom%Sg   &
            ,       rcg => iflow%geom%rcg  &
            )

         call post_processing(nx, ny, sem_a, sem_g, w_g, w_cam, ccTw, ig, lid, modtur &
            , sim_id, x, y, rp, re, xe, ye, xk, yk, Jp, Sg, rcg, gamma, pr, go, po, Rg &
            , tcpuo, tcpu1, tcpu2, tcpu, Fpv1D, Fd1D, fm1D, fme, cp, vlp, vtp, kp, gcp &
            , Twall, u1D, p1D, T1D, M1D, uin, vin, pin, Tin &
            , roe, Uce, Vcn, de, dn, pl, bp, xp, yp, u, v, p, T, ro)

      end associate

   end subroutine


  subroutine post_processing(nx, ny, sem_a, sem_g, w_g, w_cam, ccTw, ig, lid, modtur &
       , sim_id, x, y, rp, re, xe, ye, xk, yk, Jp, Sg, rcg, gamma, pr, go, po, Rg &
       , tcpuo, tcpu1, tcpu2, tcpu, Fpv1D, Fd1D, fm1D, fme, cp, vlp, vtp, kp, gcp &
       , Twall, u1D, p1D, T1D, M1D, uin, vin, pin, Tin &
       , roe, Uce, Vcn, de, dn, pl, bp, xp, yp, u, v, p, T, ro)
    implicit none
    integer, intent(in) :: nx  ! Number of volumes in csi direction (real+fictitious)
    integer, intent(in) :: ny  ! Number of volumes in eta direction (real+fictitious)
    integer, intent(in) :: sem_a ! 1 = do not open result files, 0 = open
    integer, intent(in) :: sem_g ! 0 = visualize the plot, 1 = do not visualize
    integer, intent(in) :: w_g   ! Frequency of writing data for graphics
    integer, intent(in) :: w_cam ! 1 = write the fields, 0 = do not
    integer, intent(in) :: ccTw  ! ccTw = 0 -> adiabatic;  ccTw = 1 -> prescribed temperature
    integer, intent(in) :: ig    ! i=ig when the throttle cross section is in the east face of the CV
    integer, intent(in) :: lid   ! lid = listing file id
    integer, intent(in) :: modtur ! turbulence model option (0=laminar, 1=Baldwin-Lomax)
    character (len = *), intent(in) :: sim_id     ! Simulation identification
    real(8), dimension (nx*ny), intent(in) :: x   ! Coord. x of the northest corner of volume P
    real(8), dimension (nx*ny), intent(in) :: y   ! Coord. y of the northeast corner of volume P
    real(8), dimension (nx*ny), intent(in) :: rp  ! Radius of the center of volume P
    real(8), dimension (nx*ny), intent(in) :: re  ! Radius of the center of east face of volume P
    real(8), dimension (nx*ny), intent(in) :: xe  ! x_eta at center of face east
    real(8), dimension (nx*ny), intent(in) :: ye  ! y_eta at center of face east
    real(8), dimension (nx*ny), intent(in) :: xk  ! x_csi at center of face north
    real(8), dimension (nx*ny), intent(in) :: yk  ! y_csi at center of face north
    real(8), dimension (nx*ny), intent(in) :: Jp  ! Jacobian at the center of volume P
    real(8), intent(in) :: Sg     ! Throttle area
    real(8), intent(in) :: rcg    ! Curvature radius of throttle
    real(8), intent(in) :: gamma  ! Specific heat ratio
    real(8), intent(in) :: pr     ! atmospheric pressure at the sea level
    real(8), intent(in) :: go     ! gravitational acceleration at the sea level
    real(8), intent(in) :: po     ! Stagnation pressure
    real(8), intent(in) :: Rg     ! Perfect gas constant
    real(8), intent(in) :: tcpuo  ! Accumulated time of simulation (before simulation interruption)
    real(8), intent(in) :: tcpu1  ! Initial time (present simulation)
    real(8), intent(in) :: tcpu2  ! Final time (present simulation)
    real(8), intent(in) :: tcpu   ! Total time (before and after simulation interruption)
    real(8), intent(in) :: Fpv1D  ! Pressure thrust into vacuum (Isentropic 1D flow)
    real(8), intent(in) :: Fd1D   ! Dynamic thrust (Isentropic 1D flow)
    real(8), intent(in) :: fm1D   ! Mass flow rate (Isentropic 1D flow)
    real(8), intent(in) :: fme    ! Mass flow rate at exit
    real(8), dimension (nx*ny), intent(in) :: cp    ! Specific heat at constant pressure at center o volume P
    real(8), dimension (nx*ny), intent(in) :: vlp   ! Laminar viscosity at center of volume P
    real(8), dimension (nx*ny), intent(in) :: vtp   ! Eddy viscosity at center of volume P
    real(8), dimension (nx*ny), intent(in) :: kp    ! Thermal conductivity at center of volume P
    real(8), dimension (nx*ny), intent(in) :: gcp   ! gcp = gamma = Cp/Cv at center of CV P
    real(8), dimension (nx),    intent(in) :: Twall ! Wall temperature
    real(8), dimension (nx),    intent(in) :: u1D   !  Velocity. Isentropic flow
    real(8), dimension (nx),    intent(in) :: p1D   !  Pressure. Isentropic flow
    real(8), dimension (nx),    intent(in) :: T1D   !  Temperature. Isentropic flow
    real(8), dimension (nx),    intent(in) :: M1D   !  Mach number. Isentropic flow
    real(8), dimension (ny),    intent(in) :: uin   ! Velocity u in the entrance
    real(8), dimension (ny),    intent(in) :: vin   ! Velocity v in the entrance
    real(8), dimension (ny),    intent(in) :: pin   ! Pressure in the entrance
    real(8), dimension (ny),    intent(in) :: Tin   ! Temperature in the entrance
    real(8), dimension (nx*ny), intent(in) :: roe   ! Absolute density at east face
    real(8), dimension (nx*ny), intent(in) :: Uce   ! Contravariant velocity U at east face
    real(8), dimension (nx*ny), intent(in) :: Vcn   ! Contravariant velocity V at north face
    real(8), dimension (nx*ny), intent(in) :: de    ! SIMPLEC coefficient for Uce
    real(8), dimension (nx*ny), intent(in) :: dn    ! SIMPLEC coefficient for Vcn
    real(8), dimension (nx*ny), intent(in) :: pl    ! Pressure correction
    real(8), dimension (nx*ny), intent(in) :: bp    ! Source vector of the linear system of pl
    real(8), dimension (nx*ny), intent(inout) :: xp ! Coord. of the centroid of volume P
    real(8), dimension (nx*ny), intent(inout) :: yp ! Coord. of the centroid of volume P
    real(8), dimension (nx*ny), intent(inout) :: u  ! Cartesian velocity of the last iteraction
    real(8), dimension (nx*ny), intent(inout) :: v  ! Cartesian velocity of the last iteraction
    real(8), dimension (nx*ny), intent(inout) :: p  ! Pressure at center o volume P
    real(8), dimension (nx*ny), intent(inout) :: T  ! Temperature of the last iteraction
    real(8), dimension (nx*ny), intent(inout) :: ro ! Absolute density at center of vol. P

    ! Auxiliary variables
    real(8) :: rcag  ! Dimensionless curvature radius of throttle
    real(8), dimension (nx*ny) :: psi  ! Stream function
    real(8), dimension (nx*ny) :: M    ! Mach number at the center of the volumes

    ! Writes grid parameters
    call write_grid_parameters( nx, ny, lid, x, y) ! No output

    ! Writes nozzle parameters
    call write_nozzle_parameters( nx, ny, lid, ig, rcg, x, y, rcag) ! Output: last entry

    ! Writes grid boundary nodes
    call write_grid_boundary( nx, ny, w_g, sim_id, x, y)

    ! Plots grid boundaries
    call plotter04( sem_g, sim_id)

    ! Writes grid vertices
    call write_grid(nx, ny, w_g, sim_id, x, y)

    ! Plots the grid
    call plotter05( sem_g, sim_id)

    ! Calculates the stream function
    call get_stream_function( nx, ny, re, roe, Uce, psi) ! Output: last entry


    ! Calculates the boundary values and store them in the fictitious
    call post_processing_boundaries( nx, ny, ccTw, x, y, Twall, uin, vin &
         ,                                 pin, Tin, xp, yp, u, v, p, T ) ! InOutput: last six entries

    ! Calculates the Mach number field
    M = dsqrt ( (u**2 + v**2) / (gcp * Rg * T) )

    ! Calculates the specific mass field
    ro = p / ( Rg * T )

    ! Write fields

    call write_main_fields(nx, ny, sim_id, xp, yp, p, ro, T, u, v, M, psi, vtp)

    call plotter06( sem_g, sim_id)

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

        if ( modtur == 1 ) then

            write(lid,"(//,30x,'vtp: nodal eddy viscosity (Pa.s)')")

            call write_field(nx, ny, lid, xp, yp, vtp)

        end if


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

       write(lid,20)
20     format(//,30x,'psi: stream function (kg/s)')
       call write_field(nx, ny, lid, xp, yp, psi)

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
    call write_main_results( nx, ny, lid, sim_id, x, y, re, ye, Sg, rcag &
         ,           gamma, pr, go, po, Fpv1D, Fd1D, fm1D, fme, p, roe, u, Uce)


    write(lid,31) tcpuo, tcpu2-tcpu1, tcpu
31  format(//,f14.2,' = tcpuo: acumulated CPU time (s) (before interuption)', &
         /,f14.2,' = dtcpu: CPU time (s) (after interuption)', &
         /,f14.2,' = tcpu:  total CPU time (s)')


    ! Plots the residual as a function of the iteractions
    call plotter02( sem_g, sim_id, "'Residual of the linear systems'")

    ! show the listing file of iteractions
    if ( sem_a == 0 ) call system(text_editor // './mach2d_output/mach2d_' // trim(adjustl(sim_id)) // '-residual.dat')

    ! close the main listing file
    close (lid)

    ! show main listing file
    if ( sem_a == 0 ) call system(text_editor // './mach2d_output/mach2d_' // trim(adjustl(sim_id)) // '.lst')

    ! show the file with main results
    if ( sem_a == 0 ) call system(text_editor // './mach2d_output/mach2d_' // trim(adjustl(sim_id)) // '.Richardson_3p0')

    ! *** generates and shows some graphics ***

    call plotter01(nx, ny, w_g, sem_g, xp, u1D, u, sim_id &
         , fieldname = "u-velocity" &
         , title = "' u velocity (m/s) '"  &
         , xlabel = "' x (m) '" &
         , ylabel = "'u (m/s)'" &
         , key = "right bottom" )

    call plotter01(nx, ny, w_g, sem_g, xp, 0.d0*u1D, v, sim_id &
         , fieldname = "v-velocity" &
         , title = "' v velocity (m/s) '"  &
         , xlabel = "' x (m) '" &
         , ylabel = "' v (m/s) '" &
         , key = "right bottom" )

    call plotter01(nx, ny, w_g, sem_g, xp, p1D, p, sim_id &
         , fieldname = "pressure" &
         , title = "' Pressure (Pa) '"  &
         , xlabel = "' x (m)'" &
         , ylabel = "' p (Pa)'" &
         , key = "right top" )

    call plotter01(nx, ny, w_g, sem_g, xp, T1D, T, sim_id &
         , fieldname = "temperature" &
         , title = "'Temperature (K)'"  &
         , xlabel = "' x (m)'" &
         , ylabel = "'T (K)'" &
         , key = "right top" )

    call plotter01(nx, ny, w_g, sem_g, xp, p1D / ( Rg * T1D ), p/(Rg*T), sim_id &
         , fieldname = "specific-mass" &
         , title = "'Specific mass (kg/m3)'"  &
         , xlabel = "'x (m)'" &
         , ylabel = "'ro (kg/m3)'" &
         , key = "right top" )

    call plotter01(nx, ny, w_g, sem_g, xp, M1D, M, sim_id &
         , fieldname = "mach" &
         , title = "'Mach number'"  &
         , xlabel = "'x (m)'" &
         , ylabel = "'M'" &
         , key = "right bottom" )

    call plotter03( nx, ny, ig, w_g, sem_g, sim_id, xp, rp, re, ye, u, p, ro, pr) ! No output

  end subroutine post_processing

  subroutine plotter01(nx, ny, w_g, sem_g, xp, var1D, var, sim_id, fieldname, title, xlabel, ylabel, key)
    implicit none
    integer, intent(in) :: nx    ! Number of volumes in csi direction (real+fictitious)
    integer, intent(in) :: ny    ! Number of volumes in eta direction (real+fictitious)
    integer, intent(in) :: w_g   ! Frequency of writing data for graphics
    integer, intent(in) :: sem_g ! ( 0 = visualize the plot, 1 = do not visualize)

    real(8), dimension (nx*ny), intent(in) :: xp    ! Coord. of the centroid of volume P
    real(8), dimension (nx),    intent(in) :: var1D ! Solution to isentropic flow
    real(8), dimension (nx*ny), intent(in) :: var   ! Solution to the 2D flow

    character (len = *), intent(in) :: sim_id    ! Simulation identification
    character (len = *), intent(in) :: fieldname ! Name of data file
    character (len = *), intent(in) :: title     ! Graphics title
    character (len = *), intent(in) :: xlabel    ! Graphics xlabel
    character (len = *), intent(in) :: ylabel    ! Graphics ylabel
    character (len = *), intent(in) :: key       ! Graphics key

    ! Auxiliary variables
    integer :: i, j, nps, npn
    character (len=250) :: str1, str2, str3, backslash

    !
    ! Generating data file
    !

    str1 = "./mach2d_output/mach2d_"      &
         // trim(adjustl(sim_id)) // "-" &
         // trim(adjustl(fieldname)) // ".dat"

    open(10, file = str1)

    write(10, "('# ', A)")  trim(adjustl(title))

    write(10, "('# ', A22, 3(1X,A23) )") "x (m)", "1D", "symmetry line", "wall"

    do i = 1, nx

       j = 1
       nps = (j-1)*nx + i

       j = ny
       npn = (j-1)*nx + i

       if ( i==1 .or. i==nx .or. mod(i,w_g)==0 ) then

          write(10,"(4(1pe24.15))") xp(nps), var1D(i), var(nps), var(npn)

       end if

    end do

    close(10)

    !
    ! Generating gnuplot file
    !

    backslash = "\ "

    str2 = "./mach2d_output/mach2d_"     &
         // trim(adjustl(sim_id)) // "-" &
         // trim(adjustl(fieldname)) // ".gnu"

    str3 = "./mach2d_output/mach2d_"     &
         // trim(adjustl(sim_id)) // "-" &
         // trim(adjustl(fieldname)) // ".png"

    open(10, file = str2)

    write(10,*) "set terminal pngcairo size 1280,720 enhanced font 'Verdana,16'"
    write(10,*) "set out '", trim(adjustl(str3)), "'"
    write(10,*) "set style data linespoints"
    write(10,*) "set grid"
    write(10,*) "set xlabel ", xlabel
    write(10,*) "set ylabel ", ylabel
    write(10,*) "set time"
    write(10,*) "set key ", key
    write(10,*) "set title ", title

    write(10,*) "plot '", trim(adjustl(str1)), "' using 1:2 title '1D' " // trim(backslash)
    write(10,*) "    ,'", trim(adjustl(str1)), "' using 1:3 title 'center' " // trim(backslash)
    write(10,*) "    ,'", trim(adjustl(str1)), "' using 1:4 title 'wall'"

    close(10)

    call system( graph_generator // str2 )

    if ( sem_g == 0 ) call system( graph_viewer // str3 )

  end subroutine plotter01



  ! Plots the residual as a function of the iteractions
  subroutine plotter02( sem_g, sim_id, title)
    implicit none
    integer, intent(in) :: sem_g ! ( 0 = visualize the plot, 1 = do not visualize)
    character (len = *), intent(in) :: sim_id    ! Simulation identification
    character (len = *), intent(in) :: title     ! Graphics title

    ! Auxiliary variables
    character (len=250) :: str1, str2, str3, backslash

    !
    ! Generating gnuplot file
    !

    backslash = "\ "

    str1 = "./mach2d_output/mach2d_" &
         // trim(adjustl(sim_id))    &
         // "-residual.dat"

    str2 = "./mach2d_output/mach2d_" &
         // trim(adjustl(sim_id))    &
         // "-residual.gnu"

    str3 = "./mach2d_output/mach2d_" &
         // trim(adjustl(sim_id))    &
         // "-residual.png"

    open(10, file = str2)

    write(10,*) "set terminal pngcairo size 1280,720 enhanced font 'Verdana,16'"
    write(10,*) "set out '", trim(adjustl(str3)), "'"
    write(10,*) "set style data linespoints"
    write(10,*) "set grid"
    write(10,*) "set logscale y"
    write(10,*) "set xlabel 'iteractions'"
    write(10,*) "set ylabel 'norm L1 of the total residual (dimensionless)'"
    write(10,*) "set time"
    write(10,*) "set title ", title

    write(10,*) "plot '", trim(adjustl(str1)), "' using 1:2 title '' "

    close(10)

    call system( graph_generator // str2 )

    if ( sem_g == 0 ) call system( "evince " // str3 )

  end subroutine plotter02

  subroutine plotter03( nx, ny, ig, w_g, sem_g, sim_id, xp, rp, re, ye, u, p, ro, pr) ! No output
    implicit none
    integer, intent(in) :: nx    ! Number of volumes in csi direction (real+fictitious)
    integer, intent(in) :: ny    ! Number of volumes in eta direction (real+fictitious)
    integer, intent(in) :: ig    ! i=ig when the throttle cross section is in the east face of the CV
    integer, intent(in) :: w_g   ! Frequency of writing data for graphics
    integer, intent(in) :: sem_g ! ( 0 = visualize the plot, 1 = do not visualize)
    character (len = *), intent(in) :: sim_id     ! Simulation identification
    real(8), dimension (nx*ny), intent(in) :: xp  ! Coord. of the centroid of volume P
    real(8), dimension (nx*ny), intent(in) :: rp  ! Radius of the center of volume P
    real(8), dimension (nx*ny), intent(in) :: re  ! Radius of the center of east face of volume P
    real(8), dimension (nx*ny), intent(in) :: ye  ! y_eta at face east of volume P
    real(8), dimension (nx*ny), intent(in) :: u   ! Cartesian velocity of the last iteraction
    real(8), dimension (nx*ny), intent(in) :: p   ! Pressure at center o volume P
    real(8), dimension (nx*ny), intent(in) :: ro  ! Absolute density at center of vol. P
    real(8), intent(in) :: pr    ! atmospheric pressure at the sea level

    ! Constants
    real(8) :: pi = acos(-1.d0)

    ! Auxiliary variables
    integer :: i, j, np, nps
    real(8) :: Fd, Fp,  Fpv, Fv, F, rpp, dy
    character (len=250) :: str1, str2, str3, backslash

    !
    ! Generating data file
    !

    str1 = "./mach2d_output/mach2d_" // trim(adjustl(sim_id)) // '-thrust.dat'

    open(10,file=trim(adjustl(str1)))

    write(10,*) '### THRUST (N) ###'
    write(10,"(A2,A15,5(A17))") '#', 'x (m)', 'dynamic', 'pressure-sea' &
         , 'pressure-vacuum', 'total-sea', 'total-vacuum'

    do i = ig+1, nx
       Fd  = 0.0d0
       Fp  = 0.0d0
       Fpv = 0.0d0
       do j = 2, ny-1
          np  = (j-1)*nx + i
          rpp = rp(np)
          dy  = (ye(np-1)+ye(np))/2
          if ( i == nx ) then
             rpp = re(np-1)
             dy  = ye(np-1)
          end if
          Fd  = Fd  + rpp * ro(np) * dy * (u(np)**2)
          Fp  = Fp  + (p(np) - pr) * rpp * dy
          Fpv = Fpv +  p(np)       * rpp * dy
       end do
       Fd  = Fd  * 2 * pi
       Fp  = Fp  * 2 * pi
       Fpv = Fpv * 2 * pi
       F   = Fd + Fp
       Fv  = Fd + Fpv
       j = 1
       nps = (j-1)*nx + i
       if ( i==ig+1 .or. i==nx .or. mod(i,w_g)==0 ) then
          write(10,"(6(1pe17.5))") xp(nps), Fd, Fp, Fpv, F, Fv
       end if
    end do
    close(10)

    !
    ! Generating gnuplot file
    !

    backslash = "\ "

    str2 = "./mach2d_output/mach2d_" // trim(adjustl(sim_id)) // '-thrust.gnu'
    str3 = "./mach2d_output/mach2d_" // trim(adjustl(sim_id)) // '-thrust.png'

    open(10,file= str2 )

    write(10,*) "set terminal pngcairo size 1280,720 enhanced font 'Verdana,16'"
    write(10,*) "set out '", trim(adjustl(str3)), "'"
    write(10,*) "set style data linespoints"
    write(10,*) "set grid"
    write(10,*) "set xlabel 'x (m)'"
    write(10,*) "set ylabel 'thrust (N)'"
    write(10,*) "set time"
    write(10,*) "set key outside left top"
    write(10,*) "set title 'Thrust " // trim(adjustl(sim_id)) // "'"

    write(10,*) "plot '", trim(adjustl(str1)), "' using 1:2 title 'dynamic' " // trim(backslash)
    write(10,*) "    ,'", trim(adjustl(str1)), "' using 1:3 title 'pressure-sea' " // trim(backslash)
    write(10,*) "    ,'", trim(adjustl(str1)), "' using 1:4 title 'pressure-vacuum'" // trim(backslash)
    write(10,*) "    ,'", trim(adjustl(str1)), "' using 1:5 title 'total-sea'" // trim(backslash)
    write(10,*) "    ,'", trim(adjustl(str1)), "' using 1:6 title 'total-vacuum'"
    close(10)

    call system( graph_generator // str2 )

    if ( sem_g == 0 ) call system( "evince " // str3 )

  end subroutine plotter03


  subroutine plotter04( sem_g, sim_id)
    implicit none
    integer, intent(in) :: sem_g ! ( 0 = visualize the plot, 1 = do not visualize)
    character (len = *), intent(in) :: sim_id    ! Simulation identification

    ! Auxiliary variables
    character (len=250) :: str1, str2, str3, backslash

    !
    ! Generating gnuplot file
    !

    backslash = "\ "

    str2 = "./mach2d_output/mach2d_" // trim(adjustl(sim_id)) // '-boundary.gnu'
    str3 = "./mach2d_output/mach2d_" // trim(adjustl(sim_id)) // '-boundary.png'

    open(10,file= str2 )

    write(10,*) "set terminal pngcairo size 1280,720 enhanced font 'Verdana,16'"
    write(10,*) "set out '", trim(adjustl(str3)), "'"
    write(10,*) "set style data linespoints"
    write(10,*) "set grid"
    write(10,*) "set xlabel 'x (m)'"
    write(10,*) "set ylabel 'y (m)'"
    write(10,*) "set time"
    write(10,*) "set key left bottom"
    write(10,*) "set title 'Grid Boundary " // trim(adjustl(sim_id)) // "'"

    str1 = "./mach2d_output/mach2d_" // trim(adjustl(sim_id)) // '-boundary-south.dat'

    write(10,*) "plot '", trim(adjustl(str1)), "' using 4:5 title 'south' " // trim(backslash)

    str1 = "./mach2d_output/mach2d_" // trim(adjustl(sim_id)) // '-boundary-north.dat'

    write(10,*) "    ,'", trim(adjustl(str1)), "' using 4:5 title 'north' " // trim(backslash)

    str1 = "./mach2d_output/mach2d_" // trim(adjustl(sim_id)) // '-boundary-west.dat'

    write(10,*) "    ,'", trim(adjustl(str1)), "' using 4:5 title 'west' " // trim(backslash)

    str1 = "./mach2d_output/mach2d_" // trim(adjustl(sim_id)) // '-boundary-east.dat'

    write(10,*) "    ,'", trim(adjustl(str1)), "' using 4:5 title 'east' "

    close(10)

    call system( graph_generator // str2 )

    if ( sem_g == 0 ) call system( "evince " // str3 )

  end subroutine plotter04


  subroutine plotter05( sem_g, sim_id)
    implicit none
    integer, intent(in) :: sem_g ! ( 0 = visualize the plot, 1 = do not visualize)
    character (len = *), intent(in) :: sim_id    ! Simulation identification

    ! Auxiliary variables
    character (len=250) :: str1, str2, str3, backslash

    !
    ! Generating gnuplot file
    !

    backslash = "\ "

    str2 = "./mach2d_output/mach2d_" // trim(adjustl(sim_id)) // '-grid.gnu'
    str3 = "./mach2d_output/mach2d_" // trim(adjustl(sim_id)) // '-grid.png'

    open(10,file= str2 )

    write(10,*) "set terminal pngcairo size 1280,720 enhanced font 'Verdana,16'"
    write(10,*) "set out '", trim(adjustl(str3)), "'"
    write(10,*) "set style data lines"
    write(10,*) "set grid"
    write(10,*) "set xlabel 'x (m)'"
    write(10,*) "set ylabel 'y (m)'"
    write(10,*) "set time"
    write(10,*) "set title 'Grid " // trim(adjustl(sim_id)) // "'"

    str1 = "./mach2d_output/mach2d_" // trim(adjustl(sim_id)) // '-grid.dat'

    write(10,*) "plot '", trim(adjustl(str1)), "' using 1:2 title '' "

    close(10)

    call system( graph_generator // str2 )

    if ( sem_g == 0 ) call system( "evince " // str3 )

  end subroutine plotter05

  subroutine plotter06( sem_g, sim_id)
    implicit none
    integer, intent(in) :: sem_g ! ( 0 = visualize the plot, 1 = do not visualize)
    character (len = *), intent(in) :: sim_id    ! Simulation identification

    ! Auxiliary variables
    character (len=250) :: str1, str2, str3, backslash

    !
    ! Generating gnuplot file
    !

    backslash = "\ "

    str2 = "./mach2d_output/mach2d_" // trim(adjustl(sim_id)) // '-main-fields.gnu'
    str3 = "./mach2d_output/mach2d_" // trim(adjustl(sim_id))

    open(10,file= str2 )

    write(10,*) "set terminal pngcairo size 1280,720 enhanced font 'Verdana,16'"
    write(10,*) "set style data pm3d"
    write(10,*) "set ticslevel 0"
    write(10,*) "set view map"
    write(10,*) "set pm3d map"
    write(10,*) "set palette"
    write(10,*) "set xlabel 'x (m)'"
    write(10,*) "set ylabel 'y (m)'"
    write(10,*) "set time"
    write(10,*) "unset key"

    str1 = "./mach2d_output/mach2d_" // trim(adjustl(sim_id)) // '-main-fields.dat'

    write(10,*) "set title 'Fields" // trim(adjustl(sim_id)) // " - p (Pa)" // "'"
    write(10,*) "set out '", trim(adjustl(str3)), "-p-field.png'"
    write(10,*) "splot '", trim(adjustl(str1)), "' using 1:2:3 title 'p (Pa)' "

    write(10,*) "set title 'Fields" // trim(adjustl(sim_id)) // " - ro (kg/m3)" // "'"
    write(10,*) "set out '", trim(adjustl(str3)), "-ro-field.png'"
    write(10,*) "splot '", trim(adjustl(str1)), "' using 1:2:4 title 'ro (kg/m3)' "

    write(10,*) "set title 'Fields" // trim(adjustl(sim_id)) // " - T (K)" // "'"
    write(10,*) "set out '", trim(adjustl(str3)), "-T-field.png'"
    write(10,*) "splot '", trim(adjustl(str1)), "' using 1:2:5 title 'T (K)' "

    write(10,*) "set title 'Fields" // trim(adjustl(sim_id)) // " - u (m/s)" // "'"
    write(10,*) "set out '", trim(adjustl(str3)), "-u-field.png'"
    write(10,*) "splot '", trim(adjustl(str1)), "' using 1:2:6 title 'u (m/s)' "

    write(10,*) "set title 'Fields" // trim(adjustl(sim_id)) // " - v (m/s)" // "'"
    write(10,*) "set out '", trim(adjustl(str3)), "-v-field.png'"
    write(10,*) "splot '", trim(adjustl(str1)), "' using 1:2:7 title 'v (m/s)' "

    write(10,*) "set title 'Fields" // trim(adjustl(sim_id)) // " - M" // "'"
    write(10,*) "set out '", trim(adjustl(str3)), "-M-field.png'"
    write(10,*) "splot '", trim(adjustl(str1)), "' using 1:2:8 title 'M' "

    write(10,*) "set title 'Fields" // trim(adjustl(sim_id)) // " - psi (kg/s)" // "'"
    write(10,*) "set out '", trim(adjustl(str3)), "-psi-field.png'"
    write(10,*) "splot '", trim(adjustl(str1)), "' using 1:2:9 title 'psi (kg/s)' "

    write(10,*) "set title 'Fields" // trim(adjustl(sim_id)) // " - vtp (Pa.s)" // "'"
    write(10,*) "set out '", trim(adjustl(str3)), "-vtp-field.png'"
    write(10,*) "splot '", trim(adjustl(str1)), "' using 1:2:10 title 'vtp (Pa.s)' "

    close(10)

    call system( graph_generator // str2 )

    if ( sem_g == 0 ) call system( "evince " // str3 // "-p-field.png" )
    if ( sem_g == 0 ) call system( "evince " // str3 // "-ro-field.png" )
    if ( sem_g == 0 ) call system( "evince " // str3 // "-T-field.png" )
    if ( sem_g == 0 ) call system( "evince " // str3 // "-u-field.png" )
    if ( sem_g == 0 ) call system( "evince " // str3 // "-v-field.png" )
    if ( sem_g == 0 ) call system( "evince " // str3 // "-M-field.png" )

  end subroutine plotter06

  subroutine write_nozzle_parameters( nx, ny, unit, ig, rcg, x, y, rcag) ! Output: last entry
    implicit none
    integer, intent(in) :: nx   ! Number of volumes in csi direction (real+fictitious)
    integer, intent(in) :: ny   ! Number of volumes in eta direction (real+fictitious)
    integer, intent(in) :: unit ! Unit to which the results will be printed
    integer, intent(in) :: ig   ! i=ig when the throttle cross section is in the east face of the CV
    real(8), intent(in) :: rcg  ! Curvature radius of throttle
    real(8), dimension (nx*ny), intent(in) :: x   ! Coord. x of the northest corner of volume P
    real(8), dimension (nx*ny), intent(in) :: y   ! Coord. y of the northeast corner of volume P
    real(8), intent(out) :: rcag  ! Dimensionless curvature radius of throttle
    ! Constants
    real(8), parameter :: pi = acos(-1.d0)
    ! Auxiliary variables
    integer :: i, j, npg,  npin, npex, npw, npe
    real(8) :: rain, raex, Sg

    i    = ig
    j    = ny-1
    npg  = (j-1)*nx + i
    npw  = npg - 1
    npe  = npg + 1

    i    = 1
    j    = ny-1
    npin = (j-1)*nx + i

    i    = nx-1
    j    = ny-1
    npex = (j-1)*nx + i

    rcag = rcg / y(npg)

    Sg   = pi * ( y(npg) ** 2 )

    rain = ( y(npin) / y(npg) ) ** 2

    raex = ( y(npex) / y(npg) ) ** 2

    write(unit,1) ig, x(npg), y(npg), rcg, rcag, Sg, rain, raex

1   format(//,'*** Nozzle parameters ***', //, &
         i8,12x,'   = ig:   number of the VC in the x direction whose east face ', &
         ' coincides with the nozzle throttle',/, &
         1pe22.9,' = Xg:   x coord. of the nozzle throttle (m)',/, &
         1pe22.9,' = rg:   radius of the nozzle throttle (m)',/, &
         1pe22.9,' = rcg:  curvature radius of the nozzle throttle (m)',/, &
         1pe22.9,' = Rcag: dimensionless curvature radius of the nozzle throttle',/, &
         1pe22.9,' = Sg:   area of the nozzle throttle (m2)',/, &
         1pe22.9,' = rain: area ratio of the convergent section',/, &
         1pe22.9,' = raex: area ratio of the divergent section')

  end subroutine write_nozzle_parameters


  subroutine write_grid_parameters( nx, ny, unit, x, y)
    implicit none
    integer, intent(in) :: nx   ! Number of volumes in csi direction (real+fictitious)
    integer, intent(in) :: ny   ! Number of volumes in eta direction (real+fictitious)
    integer, intent(in) :: unit ! Unit to which the results will be printed
    real(8), dimension (nx*ny), intent(in) :: x   ! Coord. x of the northest corner of volume P
    real(8), dimension (nx*ny), intent(in) :: y   ! Coord. y of the northeast corner of volume P

    ! Auxiliary variables
    integer :: i, j, np
    real(8) :: dy, dx, dxmax, dxmin, dymax, dymin

    i  = nx-1
    j  = ny-1
    np = (j-1)*nx + i
    dxmin = x(np)
    dxmax = 0.d0
    dymin = y(np)
    dymax = 0.d0

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

    write(unit,1) dxmin, dxmax, dxmax/dxmin, dymin, dymax, dymax/dymin

1   format(//,'*** Grid parameters ***', //, &
         1pe22.9,' = dxmin: min(size of volume controls) in the x direction (m)',/, &
         1pe22.9,' = dxmax: max(size of volume controls) in the x direction (m)',/, &
         1pe22.9,' = dxmax / dxmin',/, &
         1pe22.9,' = dymin: min(size of volume controls) in the y direction (m)',/, &
         1pe22.9,' = dymax: max(size of volume controls) in the y direction (m)',/, &
         1pe22.9,' = dymax / dymin')

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

    open(7,file="./mach2d_output/mach2d_" // trim(adjustl(sim_id))//'-boundary-south.dat')

    write(7,1)

    j = 1
    do i = 1, nx-1
       np = (j-1)*nx + i
       if ( i==1 .or. i==nx-1 .or. mod(i,w_g)==0 ) &
            write(7,2) i, j, np, x(np), y(np)
    end do

    close(7)

    ! writes north boundary

    open(7,file="./mach2d_output/mach2d_" // trim(adjustl(sim_id))//'-boundary-north.dat')

    write(7,1)

    j = ny-1
    do i = 1, nx-1
       np = (j-1)*nx + i
       if ( i==1 .or. i==nx-1 .or. mod(i,w_g)==0 ) &
            write(7,2) i, j, np, x(np), y(np)
    end do

    close(7)

    ! writes west boundary

    open(7,file="./mach2d_output/mach2d_" // trim(adjustl(sim_id))//'-boundary-west.dat')

    write(7,1)

    i = 1
    do j = 1, ny-1
       np = (j-1)*nx + i
       if ( j==1 .or. j==ny-1 .or. mod(j,w_g)==0 ) &
            write(7,2) i, j, np, x(np), y(np)
    end do

    close(7)

    ! writes east boundary

    open(7,file="./mach2d_output/mach2d_" // trim(adjustl(sim_id))//'-boundary-east.dat')

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


    open(7,file= "./mach2d_output/mach2d_" // trim(adjustl(sim_id))//'-grid.dat')

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


  subroutine write_main_results( nx, ny, unit, sim_id, x, y, re, ye, Sg, rcag &
       ,           gamma, pr, go, po, Fpv1D, Fd1D, fm1D, fme, p, roe, u, Uce)
    implicit none
    integer, intent(in) :: nx   ! Number of volumes in csi direction (real+fictitious)
    integer, intent(in) :: ny   ! Number of volumes in eta direction (real+fictitious)
    integer, intent(in) :: unit ! Unit to which the results will be printed
    character (len = 100), intent(in) :: sim_id ! Simulation identification
    real(8), dimension (nx*ny), intent(in)  :: x,  y ! Coord. of the northest corner of volume P
    real(8), dimension (nx*ny), intent(in) :: re     ! Radius of the center of east face of volume P
    real(8), dimension (nx*ny), intent(in) :: ye     ! y_eta at face east of volume P
    real(8), intent(in) :: Sg     ! Throttle area
    real(8), intent(in) :: rcag   ! Dimensionless curvature radius of throttle
    real(8), intent(in) :: gamma  ! gamma = Cpo / Cvo in the chamber
    real(8), intent(in) :: pr     ! atmospheric pressure at the sea level
    real(8), intent(in) :: go     ! gravitational acceleration at the sea level
    real(8), intent(in) :: po     ! Stagnation pressure in the chamber
    real(8), intent(in) :: Fpv1D  ! Pressure thrust into vacuum (Isentropic 1D flow)
    real(8), intent(in) :: Fd1D   ! Dynamic thrust (Isentropic 1D flow)
    real(8), intent(in) :: fm1D   ! Mass flow rate (Isentropic 1D flow)
    real(8), intent(in) :: fme    ! Mass flow rate at exit
    real(8), dimension (nx*ny), intent(in) :: p      ! Pressure at center o volume P
    real(8), dimension (nx*ny), intent(in) :: roe    ! Absolute density at east face
    real(8), dimension (nx*ny), intent(in) :: u      ! Cartesian velocity at the center of VC P
    real(8), dimension (nx*ny), intent(in) :: Uce    ! Contravariant velocity U at east face

    ! Constants
    real(8), parameter :: pi = acos(-1.d0)

    ! Auxiliary variables
    integer :: i, j, np, npe
    real(8) :: Fp1D, F1D, Fv1D, Fo, ce1D, c1D, cv1D, CF1D, CFv1D, Is1D, Isv1D, Sex
    real(8) :: Fd, Fp,  F,  Fv,     ce  , c  , cv  , CF  , CFv  , Is  , Isv,   Fpv
    real(8) :: CdKL, KL1, KL2, KL3, KL4, hmedio

    ! analytic solutions

    i   = nx-1
    j   = ny-1

    np  = (j-1)*nx + i

    Sex = pi * ( y(np) ** 2 )

    Fp1D  = Fpv1D - pr*Sex

    F1D   = Fd1D + Fp1D

    Fv1D  = Fd1D + Fpv1D

    Fo    = po * Sg

    CF1D  = F1D / Fo

    CFv1D = Fv1D / Fo

    ce1D  = Fo / fm1D

    c1D   = F1D / fm1D

    cv1D  = Fv1D / fm1D

    Is1D  = c1D / go

    Isv1D = cv1D / go

    write(unit,1) fm1D, Fd1D, Fp1D, Fpv1D, F1D, Fv1D, Fo, CF1D, CFv1D, ce1D, &
         c1D,  cv1D, Is1D, Isv1D

1   format(/,1x,'*** Analytic solution of the Q1D isentropic flow ***',//, &
         1pe25.15,' = fm1D:  mass flow rate (kg/s)',/, &
         1pe25.15,' = Fd1D:  dynamic thrust (N)',/, &
         1pe25.15,' = Fp1D:  dynamic thrust at sea level (p = 101325 Pa) (N)',/, &
         1pe25.15,' = Fpv1D: thrust of pressure in the vacuum (N)',/, &
         1pe25.15,' = F1D:   total thrust at sea level (N)',/, &
         1pe25.15,' = Fv1D:  total thrust in the vacuum (N)',/, &
         1pe25.15,' = Fo:    standard thrust (N)',/, &
         1pe25.15,' = CF1D:  thrust coefficient at sea level (dimensionless)',/, &
         1pe25.15,' = CFv1D: thrust coefficient in the vacuum (dimensionless)',/, &
         1pe25.15,' = ce1D:  characteristic velocity (m/s)',/, &
         1pe25.15,' = c1D:   velocity of efective ejection at the sea level (m/s)',/, &
         1pe25.15,' = cv1D:  velocity of efective ejection in the vacuum (m/s)',/, &
         1pe25.15,' = Is1D:  specific impulse at sea level (s)',/, &
         1pe25.15,' = Isv1D: specific impulse in the vacuum (s)')

    KL1 = ( gamma + 1 ) / ( ( 1 + rcag ) ** 2 )
    KL2 = ( 8 * gamma - 27 ) / ( 2304 * ( 1 + rcag ) )
    KL3 = 754 * ( gamma ** 2 ) - ( 757 * gamma ) + 3633
    KL4 = 276480 * ( ( 1 + rcag ) ** 2 )

    CdKL = 1 - KL1 * ( ( 1 / 9.6d+1 ) - KL2 + ( KL3 / KL4 ) )

    write(unit,2) CdKL

2   format(/,1x,'*** Analytic solution 2D ***',//,  &
         1pe25.15,' = CdKL: discharge coefficient of Kliegel and Levine (dimensionless)')

    ! numeric solutions

    i   = nx-1
    Fd  = 0.0d0
    Fp  = 0.0d0
    Fpv = 0.0d0
    do j = 2, ny-1
       np  = (j-1)*nx + i
       npe = np + 1
       Fd  = Fd  + re(np) * roe(np) * Uce(np) * u(npe)
       Fp  = Fp  + (p(npe) - pr) * re(np) * ye(np)
       Fpv = Fpv + p(npe) * re(np) * ye(np)
    end do
    Fd  = Fd  * 2 * pi
    Fp  = Fp  * 2 * pi
    Fpv = Fpv * 2 * pi

    F   = Fd + Fp

    Fv  = Fd + Fpv

    CF  = F / Fo

    CFv = Fv / Fo

    ce  = Fo / fme

    c   = F / fme

    cv  = Fv / fme

    Is  = c / go

    Isv = cv / go

    write(unit,3) fme, Fd, Fp, Fpv, F, Fv, Fo, CF, CFv, ce, c, cv, Is, Isv

3   format(/,1x,'*** Numeric solution 2D ***',//, &
         1pe25.15,' = fme: mass flow rate (kg/s)',/, &
         1pe25.15,' = Fd:  dynamic thrust (N)',/, &
         1pe25.15,' = Fp:  thrust of pressure at sea level (p = 101325 Pa) (N)',/, &
         1pe25.15,' = Fpv: thrust of pressure in the vacuum (N)',/, &
         1pe25.15,' = F:   total thrust at sea level (N)',/, &
         1pe25.15,' = Fv:  total thrust in the vacuum (N)',/, &
         1pe25.15,' = Fo:  standard thrust (N)',/, &
         1pe25.15,' = CF:  thrust coefficient at sea level (dimensionless)',/, &
         1pe25.15,' = CFv: thrust coefficient in the vacuum (dimensionless)',/, &
         1pe25.15,' = ce:  characteristic velocity (m/s)',/, &
         1pe25.15,' = c:   velocity of efective ejection at sea level (m/s)',/, &
         1pe25.15,' = cv:  velocity of efective ejection in the vacuum (m/s)',/, &
         1pe25.15,' = Is:  specific impulse at sea level (s)',/, &
         1pe25.15,' = Isv: specific impulse in the vacuum (s)')

    write(unit,4) fme/fm1D, Fd/Fd1D,   Fp/Fp1D, Fpv/Fpv1D, F/F1D,   Fv/Fv1D, &
         CF/CF1D,  CFv/CFv1D, ce/ce1D, c/c1D,     cv/cv1D, Is/Is1D, Isv/Isv1D

4   format(/,1x,'*** Efficiency: numerical solution 2D / analytic Q1D (dimensionless) ***',//, &
         1pe25.15,' = discharge coefficient',/, &
         1pe25.15,' = dynamic thrust',/, &
         1pe25.15,' = pressure thrust at sea level (p = 101325 Pa)',/, &
         1pe25.15,' = pressure thrust in the vacuum',/, &
         1pe25.15,' = total thrust at sea level',/, &
         1pe25.15,' = total thrust in the vacuum',/, &
         1pe25.15,' = thrust coefficient at sea level',/, &
         1pe25.15,' = thrust coefficient in the vacuum',/, &
         1pe25.15,' = characterist velocity',/, &
         1pe25.15,' = velocity of efective ejection at sea level',/, &
         1pe25.15,' = velocity of efective ejection in the vacuum',/, &
         1pe25.15,' = specific impulse at sea level',/, &
         1pe25.15,' = specific impulse at vacuum')

    open(12,file= "./mach2d_output/mach2d_" // trim(adjustl(sim_id))//'.Richardson_3p0')

    i   = 1
    j   = ny-1
    np  = (j-1)*nx + i
    hmedio = ( x(nx-1) - x(1) ) / ( nx - 2 )

    write(12,5) hmedio,  fme/fm1D,  Fd/Fd1D, Fp/Fp1D, Fpv/Fpv1D, F/F1D,   Fv/Fv1D, &
         CF/CF1D, CFv/CFv1D, ce/ce1D, c/c1D,   cv/cv1D,   Is/Is1D, Isv/Isv1D

5   format(  2x, 1pe25.15, ' = mean h in X (m)',/, &
         ' 1', 1pe25.15, ' = discharge coefficient',/, &
         ' 2', 1pe25.15, ' = dynamic thrust',/, &
         ' 3', 1pe25.15, ' = pressure thrust at sea level (p = 101325 Pa)',/, &
         ' 4', 1pe25.15, ' = pressure thrust in the vacuum',/, &
         ' 5', 1pe25.15, ' = total thrust at sea level',/, &
         ' 6', 1pe25.15, ' = total thrust in the vacuum',/, &
         ' 7', 1pe25.15, ' = thrust coefficient at sea level',/, &
         ' 8', 1pe25.15, ' = thrust coefficient in the vacuum',/, &
         ' 9', 1pe25.15, ' = characteristic velocity',/, &
         '10', 1pe25.15, ' = velocity of efective ejection at sea level',/, &
         '11', 1pe25.15, ' = velocity of efective ejection in the vacuum',/, &
         '12', 1pe25.15, ' = specific impulse at sea level',/, &
         '13', 1pe25.15, ' = specific impulse in the vacuum')

    close(12)

  end subroutine write_main_results


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

    subroutine write_main_fields(nx, ny, sim_id, xp, yp, p, ro, T, u, v, M, psi, vtp)
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
        real(8), dimension (nx*ny), intent(in) :: psi   ! Stream function
        real(8), dimension (nx*ny), intent(in) :: vtp   ! Eddy viscosity at center of volume P
        ! Auxiliary variables

        integer :: i, j, np
        character (len = 100) :: str1

        str1 = "./mach2d_output/mach2d_"      &
        // trim(adjustl(sim_id)) // "-" &
        // "main-fields.dat"

        open(10, file = str1)

        write(10,"(8(X,A11))") "# x", "y", "p", "ro", "T", "u", "v", "M", "psi", "vtp"

        do j = 1, ny
            do i = 1, nx
                np = (j-1)*nx + i
                write(10,"(10(X,ES11.4))") xp(np), yp(np), p(np), ro(np), T(np), u(np), v(np), M(np) &
                ,  psi(np), vtp(np)
            end do
            write(10,*)
        end do

        close(10)

    end subroutine write_main_fields

  subroutine post_processing_boundaries( nx, ny, ccTw, x, y, Twall, uin, vin &
       ,                                 pin, Tin, xp, yp, u, v, p, T ) ! InOutput: last six entries
    implicit none
    integer, intent(in) :: nx   ! Number of volumes in csi direction (real+fictitious)
    integer, intent(in) :: ny   ! Number of volumes in eta direction (real+fictitious)
    integer, intent(in) :: ccTw   ! ccTw = 0 -> adiabatic;  ccTw = 1 -> prescribed temperature
    real(8), dimension (nx*ny), intent(in)  :: x     ! Coord. x of the northest corner of volume P
    real(8), dimension (nx*ny), intent(in)  :: y     ! Coord. y of the northeast corner of volume P
    real(8), dimension (nx),    intent(in)  :: Twall ! Wall temperature
    real(8), dimension (ny),    intent(in)  :: uin   ! Velocity u in the entrance
    real(8), dimension (ny),    intent(in)  :: vin   ! Velocity v in the entrance
    real(8), dimension (ny),    intent(in)  :: pin   ! Pressure in the entrance
    real(8), dimension (ny),    intent(in)  :: Tin   ! Temperature in the entrance
    real(8), dimension (nx*ny), intent(inout) :: xp  ! Coord. of the centroid of volume P
    real(8), dimension (nx*ny), intent(inout) :: yp  ! Coord. of the centroid of volume P
    real(8), dimension (nx*ny), intent(inout) :: u   ! Cartesian velocity of the last iteraction
    real(8), dimension (nx*ny), intent(inout) :: v   ! Cartesian velocity of the last iteraction
    real(8), dimension (nx*ny), intent(inout) :: p   ! Pressure at center o volume P
    real(8), dimension (nx*ny), intent(inout) :: T   ! Temperature at center o volume P

    ! Auxiliary variables
    integer :: i, j, np, npn, npw, nps, npe, npsw

    ! west boundary (entrance)
    i = 1
    do j = 2, ny-1
       np = (j-1)*nx + i
       nps = np - nx
       u(np) = uin(j)
       v(np) = vin(j)
       p(np) = pin(j)
       T(np) = Tin(j)
       xp(np) = x(np)
       yp(np) = ( y(nps) + y(np) ) / 2
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
       xp(np) = x(npw)
       yp(np) = ( y(npw-nx) + y(npw) ) / 2
    end do

    ! lower boundary (symmetry line)
    j = 1
    do i = 2, nx-1
       np   = (j-1)*nx + i
       npw  = np - 1
       npn  = np + nx
       u(np) = u(npn)
       v(np) = 0.0d0
       p(np) = p(npn)
       T(np) = T(npn)
       xp(np) = ( x(npw) + x(np) ) / 2
       yp(np) = ( y(npw) + y(np) ) / 2
    end do

    ! SW corner
    j = 1
    i = 1
    np   = (j-1)*nx + i
    npn  = np + nx
    u(np) = u(npn)
    v(np) = 0.0d0
    p(np) = p(npn)
    T(np) = T(npn)
    xp(np) = x(np)
    yp(np) = y(np)

    ! SE corner
    j = 1
    i = nx
    np   = (j-1)*nx + i
    npw  = np - 1
    npn  = np + nx
    u(np) = u(npn)
    v(np) = 0.0d0
    p(np) = p(npn)
    T(np) = T(npn)
    xp(np) = x(npw)
    yp(np) = y(npw)

    ! upper boundary (wall)
    j = ny
    do i = 2, nx-1
       np  = (j-1)*nx + i
       nps = np - nx
       u(np) = (u(np) + u(nps)) / 2
       v(np) = (v(np) + v(nps)) / 2
       p(np) = p(nps)
       if ( ccTw == 1 ) then ! T prescribed
          T(np) = Twall(i)
       else                  ! adiabatic
          T(np) = T(nps)
       end if
       xp(np) = ( x(nps) + x(nps-1) ) / 2
       yp(np) = ( y(nps) + y(nps-1) ) / 2
    end do

    ! NW and NE corners
    j = ny
    do i = 1, nx, nx-1
       np  = (j-1)*nx + i
       nps = np - nx
       npw = np - 1
       npe = np + 1
       npsw = nps - 1

       p(np) = p(nps)
       if ( i == 1 ) then
          xp(np) = x(nps)
          yp(np) = y(nps)
          u(np) = u(npe)
          v(np) = v(npe)
          T(np) = T(npe)
       end if
       if ( i == nx ) then
          xp(np) = x(nps-1)
          yp(np) = y(nps-1)
          u(np) = u(npw)

          ! Since v(npw) and v(nps) were modified, being now
          ! v(npw)_new = (v(npw)_old + v(npsw)_old) /2
          ! v(nps)_new = (v(nps)_old + v(npsw)_old) /2
          ! and we need then mean
          ! (v(np)_old + v(npw)_old + v(nps)_old + v(npsw)_old ) / 4
          ! we proceed as follow:
          v(np) = ( v(np) + 2.d0 * v(npw) + 2.d0*v(nps) - v(npsw) ) / 4.d0

          T(np) = T(npw)
       end if
    end do

  end subroutine post_processing_boundaries



  subroutine get_stream_function( nx, ny, re, roe, Uce, psi) ! Output: last entry
    implicit none
    integer, intent(in) :: nx  ! Number of volumes in csi direction (real+fictitious)
    integer, intent(in) :: ny  ! Number of volumes in eta direction (real+fictitious)
    real(8), dimension (nx*ny), intent(in)  :: re  ! Radius of the center of east face of volume P
    real(8), dimension (nx*ny), intent(in)  :: roe ! Absolute density at east face
    real(8), dimension (nx*ny), intent(in)  :: Uce ! Contravariant velocity U at east face
    real(8), dimension (nx*ny), intent(out) :: psi ! Stream function

    ! Constants
    real(8), parameter :: pi = acos(-1.d0)

    ! Auxiliary variables
    integer :: i, j, np, nps, npn, npw, npe

    psi = 0.d0

    do i = 1, nx-1
       do j = 2, ny-1

          np  = nx * (j-1) + i
          nps = np - nx

          psi(np) = psi(nps) + re(np) * roe(np) * Uce(np)

       end do
    end do

    psi = psi !* 2 * pi

    ! Extrapolation to the fictitious volumes

    ! South boundary

    j = 1

    do i = 2, nx-1

        np  = nx * (j-1) + i

        npn  = np + nx

        psi(np) = psi(npn)

    end do

    ! North boundary

    j = ny

    do i = 2, nx-1

        np  = nx * (j-1) + i

        nps  = np - nx

        psi(np) = psi(nps)

    end do

    ! West boundary

    i = 1

    do j = 1, ny

        np  = nx * (j-1) + i

        npe  = np + 1

        psi(np) = psi(npe)

    end do

    ! East boundary

    i = nx

    do j = 1, ny

        np  = nx * (j-1) + i

        npw  = np - 1

        psi(np) = psi(npw)

    end do

  end subroutine get_stream_function


end module
