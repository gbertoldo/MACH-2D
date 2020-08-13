!>
!! \brief "mod_mach2d_data" contains data for running the main subroutine of
!!        MACH-2D. It also provides some initialization procedures.
!!
module mod_mach2d_data

   use mod_class_ifile
   use mod_class_thermophysical_abstract
   use mod_class_solver_abstract
   use mod_grid

   implicit none

   !
   ! INPUT PARAMETERS READER
   !
   type(class_ifile) ifile

   !
   ! GRID PARAMETERS AND VARIABLES
   !
   integer :: nx    ! number of volumes in the csi direction (real+fictitious)
   integer :: ny    ! number of volumes in the eta direction (real+fictitious)
   integer :: nxy   ! nxy = nx * ny
   integer :: coord ! Kind of coord. system ( 1=cylindrical, 0 = cartesian)

   real(8), allocatable, dimension(:) :: x   ! coord. at the northeast corner of the volume P
   real(8), allocatable, dimension(:) :: y   ! coord. at the northeast corner of the volume P
   real(8), allocatable, dimension(:) :: xp  ! Coord. x of the centroid of volume P
   real(8), allocatable, dimension(:) :: yp  ! Coord. y of the centroid of volume P
   real(8), allocatable, dimension(:) :: xe  ! x_eta at face east of volume P
   real(8), allocatable, dimension(:) :: ye  ! y_eta at face east of volume P
   real(8), allocatable, dimension(:) :: xen ! x_eta at face north of volume P
   real(8), allocatable, dimension(:) :: yen ! y_eta at face north of volume P
   real(8), allocatable, dimension(:) :: xk  ! x_csi at face north of volume P
   real(8), allocatable, dimension(:) :: yk  ! y_csi at face north of volume P
   real(8), allocatable, dimension(:) :: xke ! x_csi at face east of volume P
   real(8), allocatable, dimension(:) :: yke ! y_csi at face east of volume P
   real(8), allocatable, dimension(:) :: Jp  ! Jacobian at the center of volume P
   real(8), allocatable, dimension(:) :: Je  ! Jacobian at the center of east face of volume P
   real(8), allocatable, dimension(:) :: Jn  ! Jacobian at the center of north face of volume P
   real(8), allocatable, dimension(:) :: alphae ! (metric) Alpha at the center of east face of volume P
   real(8), allocatable, dimension(:) :: gamman ! (metric) Gamma at the center of north face of volume P
   real(8), allocatable, dimension(:) :: betae  ! (metric) Beta  at the center of east face of volume P
   real(8), allocatable, dimension(:) :: betan  ! (metric) Beta  at the center of north face of volume P
   real(8), allocatable, dimension(:) :: radius ! Radius of northest corner of volume P
   real(8), allocatable, dimension(:) :: re  ! Radius of the center of east face of volume P
   real(8), allocatable, dimension(:) :: rn  ! Radius of the center of north face of volume P
   real(8), allocatable, dimension(:) :: rp  ! Radius of the center of volume P

   !
   ! NUMERIC PARAMETERS AND VARIABLES
   !
   integer :: it     ! Iteractions counter
   integer :: itm    ! Iteractions counter for the mass cycle
   integer :: itp    ! Iteractions counter for the pressure cycle
   integer :: itmmax ! Maximum number of iteractions for mass cycle
   integer :: itpmax ! Maximum number of iteractions for pressure cycle
   integer :: itmax  ! Maximum number of iteractions for time cycle
   integer :: itemax ! Maximum number of iteractions for extrapolation to fictitious
   integer :: it_stop! Number of iteractions at which the time evolution cycle is interrupted
   integer :: wlf    ! Frequency of printing in the listing file
   integer :: sem_a  ! 1 = do not open result files, 0 = open
   integer :: sem_g  ! 0 = visualize the plot, 1 = do not visualize
   integer :: w_g    ! Frequency of writing data for graphics
   integer :: w_cam  ! 1 = write the fields, 0 = do not
   integer :: wppd   ! Write post processed data (0=no, 1=yes, 2=yes-simplified)
   real(8) :: beta   ! UDS/CDS mixing constant (0=UDS, 1=CDS)
   real(8) :: dt     ! Time step (s)
   real(8) :: norm   ! Norm of the residuals of all linear systems
   real(8) :: normpl ! Relative norm for pressure deviation: max |pl| / p_avg
   real(8) :: tcpu1  ! First time measurement
   real(8) :: tcpu2  ! Second time measurement
   real(8) :: tcpu   ! CPU time measurement (s)
   real(8) :: tolt   ! Tolerance for the time evolution cycle
   real(8) :: tolm   ! Tolerance for the mass cycle
   real(8) :: rmass  !< Norm L1 of the residual of the mass conservation equation

   ! Linear system's solvers
   class(class_solver_abstract), pointer :: solver9d
   class(class_solver_abstract), pointer :: solver5d

   real(8) :: RAM    ! RAM memory (MB)

   real(8), allocatable, dimension(:) :: bu   ! Source of the linear system for u
   real(8), allocatable, dimension(:) :: bv   ! Source of the linear system for v
   real(8), allocatable, dimension(:) :: bt   ! Source of the linear system for T
   real(8), allocatable, dimension(:) :: bp   ! Source of the linear system for pl
   real(8), allocatable, dimension(:) :: pl   ! Pressure deviation (Pa)
   real(8), allocatable, dimension(:) :: due  ! SIMPLEC coefficients for ue
   real(8), allocatable, dimension(:) :: dve  ! SIMPLEC coefficients for ve
   real(8), allocatable, dimension(:) :: dun  ! SIMPLEC coefficients for un
   real(8), allocatable, dimension(:) :: dvn  ! SIMPLEC coefficients for vn
   real(8), allocatable, dimension(:) :: de   ! SIMPLEC coefficients for Uce
   real(8), allocatable, dimension(:) :: dn   ! SIMPLEC coefficients for Vcn

   real(8), allocatable, dimension(:,:) :: au ! Coefficients of the linear system for u
   real(8), allocatable, dimension(:,:) :: av ! Coefficients of the linear system for v
   real(8), allocatable, dimension(:,:) :: at ! Coefficients of the linear system for T
   real(8), allocatable, dimension(:,:) :: ap ! Coefficients of the linear system for pl

   real(8), allocatable, dimension(:)   :: fbe  !< Face of boundary east (1 if an east boundary, 0 otherwise)
   real(8), allocatable, dimension(:)   :: fbn  !< Face of boundary north (1 if a north boundary, 0 otherwise)

   !
   ! FILE ID NUMBERS
   !
   integer, parameter :: tfid = 10  ! Temporary file id
   integer, parameter :: lid  = 100 ! Listing file id
   integer, parameter :: rid  = 101 ! Residual file id


   character (len=20)    :: date   ! System date
   character (len=20)    :: time   ! System time
   character (len = 100) :: sim_id ! Simulation identification
   character (len = 100) :: msg(2) ! Message from subroutines to be printed in the terminal (1=title, 2=message)

   !
   ! GAS PROPERTIES AND VARIABLES
   !
   class(class_thermophysical_abstract), pointer :: thermomodel !< A pointer to the thermophysical model

   integer :: ktm    !< Kind of thermophysical model ( THERMOPHYSICAL_CONSTANT, THERMOPHYSICAL_VARIABLE )
   real(8) :: Rg     !< Gas constant (J/kg.K)
   integer :: modvis !< Viscosity model (0=Euler, 1=NS)
   real(8) :: Tref   !< Reference temperature (K)
   real(8) :: Href   !< Reference total enthalpy (m2/s2)
   real(8) :: Cpref  !< Reference Cp (J/kg.K)
   real(8) :: p_avg  !< Average value of the pressure (Pa)

   real(8), allocatable, dimension(:) :: p    ! Pressure at center o volume P (Pa)
   real(8), allocatable, dimension(:) :: pa   ! Pressure of the previous time step (Pa)
   real(8), allocatable, dimension(:) :: T    ! Temperature at center of volume P (K)
   real(8), allocatable, dimension(:) :: Ta   ! Temperature of the previous time step (K)
   real(8), allocatable, dimension(:) :: ro   ! Specific mass (absolute density) at center of volumes (kg/m3)
   real(8), allocatable, dimension(:) :: roa  ! ro of the previous time step (kg/m3)
   real(8), allocatable, dimension(:) :: roe  ! Absolute density at east face (kg/m3)
   real(8), allocatable, dimension(:) :: ron  ! Absolute density at north face (kg/m3)
   real(8), allocatable, dimension(:) :: g    !< ro / p (kg/J)
   real(8), allocatable, dimension(:) :: u    ! Cartesian velocity of the last iteraction (m/s)
   real(8), allocatable, dimension(:) :: ua   ! u of the previous time step (m/s)
   real(8), allocatable, dimension(:) :: v    ! Cartesian velocity of the last iteraction (m/s)
   real(8), allocatable, dimension(:) :: va   ! v of the previous time step (m/s)
   real(8), allocatable, dimension(:) :: ue   ! Cartesian velocity u at center of east face (m/s)
   real(8), allocatable, dimension(:) :: un   ! Cartesian velocity v at center of east face (m/s)
   real(8), allocatable, dimension(:) :: ve   ! Cartesian velocity u at center of north face (m/s)
   real(8), allocatable, dimension(:) :: vn   ! Cartesian velocity v at center of north face (m/s)
   real(8), allocatable, dimension(:) :: Uce  ! Contravariant velocity U at east face (m2/s)
   real(8), allocatable, dimension(:) :: Vcn  ! Contravariant velocity V at north face (m2/s)
   real(8), allocatable, dimension(:) :: cp   ! Specific heat at const pressure
   real(8), allocatable, dimension(:) :: gcp  ! gcp = gamma = Cp/Cv at center of CV P
   real(8), allocatable, dimension(:) :: vlp  ! Laminar viscosity at center of volume P
   real(8), allocatable, dimension(:) :: vle  ! Laminar viscosity at center of face east
   real(8), allocatable, dimension(:) :: vln  ! Laminar viscosity at center of face north
   real(8), allocatable, dimension(:) :: kp   ! Thermal conductivity at center of volume P
   real(8), allocatable, dimension(:) :: ke   ! Thermal conductivity at center of face east
   real(8), allocatable, dimension(:) :: kn   ! Thermal conductivity at center of face north
   real(8), allocatable, dimension(:) :: cup  ! Term of deferred correction for u
   real(8), allocatable, dimension(:) :: sup  ! Viscous term for u
   real(8), allocatable, dimension(:) :: cvp  ! Term of deferred correction for u
   real(8), allocatable, dimension(:) :: svp  ! Viscous term for u
   real(8), allocatable, dimension(:) :: uea  ! ue of the previous time step
   real(8), allocatable, dimension(:) :: vea  ! ve of the previous time step
   real(8), allocatable, dimension(:) :: una  ! un of the previous time step
   real(8), allocatable, dimension(:) :: vna  ! vn of the previous time step

   real(8), allocatable, dimension(:) :: Tbe !< T over the faces of the east boundary (K)
   real(8), allocatable, dimension(:) :: Tbw !< T over the faces of the west boundary (K)
   real(8), allocatable, dimension(:) :: Tbn !< T over the faces of the north boundary (K)
   real(8), allocatable, dimension(:) :: Tbs !< T over the faces of the south boundary (K)

contains

   subroutine get_parameters
      implicit none

      character (len = 100) :: input_file_parameters ! Input parameters data file
      !
      call date_time(date, time)
      !
      open(10,file="./mach2d_input/mach2d_input.txt")
      ! Reading file name from which parameters will be read
      read(10,*) input_file_parameters
      close(10)

      ! Setting up ifile
      call ifile%init("./mach2d_input/" // trim(input_file_parameters), "&")

      ! Loading all parameters from input file
      call ifile%load()

      ! Getting desired parameters
      call ifile%get_value(   sim_id,   "sim_id") ! Simulation identification  (up to 100 characters)
      call ifile%get_value(    coord,    "coord") ! Kind of coord. system ( 1=cylindrical, 0 = cartesian)
      call ifile%get_value(    itmax,    "itmax") ! Maximum number of iterations for time cycle
      call ifile%get_value(   itmmax,   "itmmax") ! Maximum number of iterations for mass cycle
      call ifile%get_value(   itpmax,   "itpmax") ! Maximum number of iteractions for pressure cycle
      call ifile%get_value(   itemax,   "itemax") ! Maximum number of iteractions for extrapolation to fictitious
      call ifile%get_value(     tolm,     "tolm") ! Tolerance for the mass cycle
      call ifile%get_value(     tolt,     "tolt") ! Tolerance for the time evolution cycle
      call ifile%get_value(      wlf,      "wlf") ! Frequency of printing in the listing file
      call ifile%get_value(    sem_a,    "sem_a") ! 1 = do not open result files, 0 = open
      call ifile%get_value(    sem_g,    "sem_g") ! 0 = visualize the plot, 1 = do not visualize
      call ifile%get_value(      w_g,      "w_g") ! Frequency of writing data for graphics
      call ifile%get_value(    w_cam,    "w_cam") ! 1 = write the fields, 0 = do not
      call ifile%get_value(     wppd,     "wppd") ! Write post processed data (0=no, 1=yes, 2=yes-simplified)
      call ifile%get_value(     beta,     "beta") ! UDS/CDS mixing constant (0=UDS, 1=CDS)
      call ifile%get_value(   modvis,   "modvis") ! Viscosity model (0=Euler, 1=NS)

      ! Initializing the grid module
      call grid_init(ifile)

      ! Calculating the number of volumes (real+fictitious) of the desired mesh
      call grid_size( nx, ny)

      ! Number of volumes of the grid
      nxy = nx * ny

   end subroutine get_parameters


   subroutine write_parameters(fid)
      implicit none
      integer, intent(in) :: fid

  100 format( "'", A21, "'",' ....: ',A)

      write(fid,*)
      write(fid,*) "Date: ", date
      write(fid,*) "Time: ", time
      write(fid,*)
      write(fid,*) "          PARAMETERS         "
      write(fid,*)
      write(fid,                100) trim(adjustl(sim_id)) , " Simulation identification  (up to 100 characters)"
      write(fid,"(I23,' ....: ',A)")    nx     , " nx     - Number of real+ficititious volumes in the csi direction (desired grid)"
      write(fid,"(I23,' ....: ',A)")    ny     , " ny     - Number of real+ficititious volumes in the eta direction (desired grid)"
      write(fid,"(I23,' ....: ',A)")    coord  , " coord  - Kind of coord. system ( 1=cylindrical, 0 = cartesian)"
      write(fid,"(I23,' ....: ',A)")    itmax  , " itmax  - Maximum number of iterations for time cycle"
      write(fid,"(I23,' ....: ',A)")    itmmax , " itmmax - Maximum number of iterations for mass cycle"
      write(fid,"(I23,' ....: ',A)")    itpmax , " itpmax - Maximum number of iteractions for pressure cycle"
      write(fid,"(I23,' ....: ',A)")    itemax , " itemax - Maximum number of iteractions for extrapolation to fictitious"
      write(fid,"(ES23.16,' ....: ',A)") tolm  , " tolm   - Tolerance for the mass cycle"
      write(fid,"(ES23.16,' ....: ',A)") tolt  , " tolt   - Tolerance for the time evolution cycle"
      write(fid,"(I23,' ....: ',A)")    wlf    , " wlf    - Frequency of printing in the listing file"
      write(fid,"(I23,' ....: ',A)")    sem_a  , " sem_a  - 1 = do not open result files, 0 = open"
      write(fid,"(I23,' ....: ',A)")    sem_g  , " sem_g  - 0 = visualize the plot, 1 = do not visualize"
      write(fid,"(I23,' ....: ',A)")    w_g    , " w_g    - Frequency of writing data for graphics"
      write(fid,"(I23,' ....: ',A)")    w_cam  , " w_cam  - 1 = write the fields, 0 = do not"
      write(fid,"(I23,' ....: ',A)")    wppd   , " wppd   - Write post processed data (0=no, 1=yes, 2=yes-simplified)"
      write(fid,"(ES23.16,' ....: ',A)") beta  , " beta   - UDS/CDS mixing constant (0=UDS, 1=CDS)"
      write(fid,"(I23,' ....: ',A)")    modvis , " modvis - Viscosity model (0=Euler, 1=NS)"
      write(fid,*)

   end subroutine write_parameters


   subroutine allocate_variables

      allocate( x(nxy),        y(nxy),         xp(nxy),         yp(nxy),         xe(nxy)  &
      ,         ye(nxy),       xen(nxy),       yen(nxy),        xk(nxy),         yk(nxy)  &
      ,         xke(nxy),      yke(nxy),       Jp(nxy),         Je(nxy),         Jn(nxy)  &
      ,         alphae(nxy),   gamman(nxy),    betae(nxy),      betan(nxy),      radius(nxy)&
      ,         re(nxy),       rn(nxy),        rp(nxy),         p(nxy),          T(nxy)   &
      ,         ro(nxy),       roe(nxy),       ron(nxy),        u(nxy),          v(nxy)   &
      ,         ue(nxy),       un(nxy),        ve(nxy),         vn(nxy),         Uce(nxy) &
      ,         Vcn(nxy),      cp(nxy),        gcp(nxy),        vlp(nxy),        vle(nxy) &
      ,         vln(nxy),      kp(nxy),        ke(nxy),         kn(nxy),         ua(nxy)  &
      ,         va(nxy),       cup(nxy),       sup(nxy),        cvp(nxy),        svp(nxy) &
      ,         bu(nxy),       bv(nxy),        roa(nxy),        due(nxy),        dve(nxy) &
      ,         dun(nxy),      dvn(nxy),       de(nxy),         dn(nxy),         uea(nxy) &
      ,         vea(nxy),      una(nxy),      vna(nxy),         bp(nxy),         pl(nxy)  &
      ,          Ta(nxy),       bt(nxy),       pa(nxy),         g(nxy) )


      allocate( Tbe(ny), Tbw(ny), Tbn(nx), Tbs(nx) )

      allocate( au(nxy,9), av(nxy,9), at(nxy,9), ap(nxy,5) )

      allocate( fbe(nxy), fbn(nxy) )

   end subroutine allocate_variables


   subroutine initialize_variables
      implicit none

      ! Reading data from input data file

      call get_parameters

      ! Openning listing file

      open(lid, file = "./mach2d_output/mach2d_" // trim(adjustl(sim_id)) // ".txt")

      ! Openning file of residuals
      open(rid, file = './mach2d_output/mach2d_' // trim(adjustl(sim_id)) // '_residual.dat' )


      write(lid,*)
      write(lid,*) "LISTING FILE OF MACH2D"
      write(lid,*)

      call write_parameters(lid)

      call allocate_variables

      ! Initializing variables

      x   = 0.d0
      y   = 0.d0
      xp  = 0.d0
      yp  = 0.d0
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

      radius = 0.d0
      re = 0.d0
      rn = 0.d0
      rp = 0.d0

      g   = 0.d0
      p   = 0.d0
      T   = 0.d0
      ro  = 0.d0
      roe = 0.d0
      ron = 0.d0
      u   = 0.d0
      v   = 0.d0
      ue  = 0.d0
      un  = 0.d0
      ve  = 0.d0
      vn  = 0.d0
      Uce = 0.d0
      Vcn = 0.d0

      cp  = 0.d0
      gcp = 0.d0
      vlp = 0.d0
      vle = 0.d0
      vln = 0.d0
      kp  = 0.d0
      ke  = 0.d0
      kn  = 0.d0

      ua  = 0.d0
      va  = 0.d0
      cup = 0.d0
      sup = 0.d0
      cvp = 0.d0
      svp = 0.d0
      bu  = 0.d0
      bv  = 0.d0
      au  = 0.d0
      av  = 0.d0
      roa = 0.d0

      due = 0.d0
      dve = 0.d0
      dun = 0.d0
      dvn = 0.d0
      de  = 0.d0
      dn  = 0.d0

      uea = 0.d0
      vea = 0.d0
      una = 0.d0
      vna = 0.d0

      Tbe = 0.d0
      Tbw = 0.d0
      Tbn = 0.d0
      Tbs = 0.d0

      ap = 0.d0
      bp = 0.d0
      pl = 0.d0

      pa = 0.d0
      Ta = 0.d0
      at = 0.d0
      bt = 0.d0

      fbe = 0.d0
      fbn = 0.d0

   end subroutine initialize_variables


   subroutine date_time(date,time)
      implicit none
      character (len=20), intent(out) :: date
      character (len=20), intent(out) :: time
      ! Auxiliary variabless
      integer(4)   :: var(8) ! date and time
      character*20 :: vardate,vartime,varzone ! date and time
      character*2  :: aux1,aux2
      character*4  :: aux3
      character*50 :: aux

      call date_and_time(vardate,vartime,varzone,var)

      write(aux,*) var(3)
      aux1 = trim(adjustl(aux))
      write(aux,*) var(2)
      aux2 = trim(adjustl(aux))
      write(aux,*) var(1)
      aux3 = trim(adjustl(aux))
      date = '('//trim(aux1)//'/'//trim(aux2)//'/'//aux3//')'

      write(aux,*) var(5)
      aux1 = trim(adjustl(aux))
      write(aux,*) var(6)
      aux2 = trim(adjustl(aux))
      write(aux,*) var(7)
      aux3 = trim(adjustl(aux))
      time = trim(aux1)//':'//trim(aux2)//':'//aux3

   end subroutine date_time

end module
