module data


   implicit none
   !
   ! GRID PARAMETERS AND VARIABLES
   !
   integer :: nx    ! number of volumes in the csi direction (real+fictitious)
   integer :: ny    ! number of volumes in the eta direction (real+fictitious)
   integer :: nxi   ! number of volumes in the csi direction (real+fictitious) of the coarsest grid
   integer :: nyi   ! number of volumes in the eta direction (real+fictitious) of the coarsest grid
   integer :: nxf   ! number of volumes in the csi direction (real+fictitious) of the finest grid
   integer :: nyf   ! number of volumes in the eta direction (real+fictitious) of the finest grid
   integer :: nxy   ! nxy = nx * ny
   integer :: nmf   ! number of the finest mesh
   integer :: nmd   ! number of the desired mesh
   integer :: kg    ! Kind of grid (1=uniform, 2=geometric progression, 3=power law, 4=gp modified)
   integer :: kcm   ! Kind of centroid mean (1=simple mean, 2=weighted mean)
   integer :: coord ! Kind of coord. system ( 1=cylindrical, 0 = cartesian)
   real(8) :: a1    ! width of the volume closer to the wall (m)
   real(8) :: avi   ! Initial value of the artificial viscosity
   real(8) :: avf   ! Final value of the artificial viscosity
   real(8) :: awf   ! Area weighting factor
   real(8) :: cbl   ! The width of the vol. closer to the wall is 'cbl' times the width of the b. layer
   real(8) :: wbl   ! Boundary layer estimated width (m)


   real(8), allocatable, dimension(:) :: x   ! coord. at the northeast corner of the volume P
   real(8), allocatable, dimension(:) :: y   ! coord. at the northeast corner of the volume P
   real(8), allocatable, dimension(:) :: xf  ! coord. at the northeast corner of the volume P of the finest grid
   real(8), allocatable, dimension(:) :: yf  ! coord. at the northeast corner of the volume P of the finest grid
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
   integer :: nitm_u ! Maximum number of iteractions for solving the linear systems for u, v and T
   integer :: nitm_p ! Maximum number of iteractions for solving the linear system for p
   integer :: wlf    ! Frequency of printing in the listing file
   integer :: sem_a  ! 1 = do not open result files, 0 = open
   integer :: sem_g  ! 0 = visualize the plot, 1 = do not visualize
   integer :: w_g    ! Frequency of writing data for graphics
   integer :: w_cam  ! 1 = write the fields, 0 = do not
   integer :: wppd   ! Write post processed data (0=no, 1=yes, 2=yes-simplified)
   integer :: it1    ! number of iteractions up to which dt = dt1
   integer :: it2    ! number of iteractions from which dt = dt2
   real(8) :: beta   ! UDS/CDS mixing constant (0=UDS, 1=CDS)
   real(8) :: dt     ! Time step (s)
   real(8) :: dt1    ! Initial time step (s)
   real(8) :: dt2    ! Final time step (s)
   real(8) :: norm   ! Norm of the residuals of all linear systems
   real(8) :: normpl ! Relative norm for pressure deviation: max |pl| / p_avg
   real(8) :: tcpu1  ! First time measurement
   real(8) :: tcpu2  ! Second time measurement
   real(8) :: tcpu   ! CPU time measurement (s)
   real(8) :: tolt   ! Tolerance for the time evolution cycle
   real(8) :: tolm   ! Tolerance for the mass cycle
   real(8) :: tol_u  !< Tolerance in the MSI for solving the linear systems for u, v and T
   real(8) :: tol_p  !< Tolerance in the MSI for solving the linear system for p
   real(8) :: curef  !< Reference value of the coefficient of convergence for u
   real(8) :: cvref  !< Reference value of the coefficient of convergence for v
   real(8) :: ctref  !< Reference value of the coefficient of convergence for T
   real(8) :: cpref  !< Reference value of the coefficient of convergence for p
   real(8) :: cref   !< Reference value of the coef. of convergence for u, v, T, p
   real(8) :: h0     !< Amplitude of h in the TSI11 model
   real(8) :: maxcc  !< Maximum allowed value of the convergence coefficient
   real(8) :: mincc  !< Minimum allowed value of the convergence coefficient
   real(8) :: rmass  !< Norm L1 of the residual of the mass conservation equation

   real(8) :: RAM    ! RAM memory (MB)


   !real(8), allocatable, dimension(:) :: ccu ! Convergence vector for u
   !real(8), allocatable, dimension(:) :: ccv ! Convergence vector for v
   !real(8), allocatable, dimension(:) :: cct ! Convergence vector for t
   !real(8), allocatable, dimension(:) :: ccp ! Convergence vector for p
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

   real(8), allocatable, dimension(:,:) :: dl9 !< Lower matrix of the MSI method for 9 diagonals
   real(8), allocatable, dimension(:,:) :: du9 !< Upper matrix of the MSI method for 9 diagonals
   real(8), allocatable, dimension(:,:) :: dl5 !< Lower matrix of the MSI method for 5 diagonals
   real(8), allocatable, dimension(:,:) :: du5 !< Upper matrix of the MSI method for 5 diagonals

   real(8), allocatable, dimension(:,:) :: a9bn !< Matrix with the coefficients of the numerical scheme for the north boundary
   real(8), allocatable, dimension(:,:) :: a9bs !< Matrix with the coefficients of the numerical scheme for the south boundary
   real(8), allocatable, dimension(:,:) :: a9be !< Matrix with the coefficients of the numerical scheme for the east boundary
   real(8), allocatable, dimension(:,:) :: a9bw !< Matrix with the coefficients of the numerical scheme for the west boundary
   real(8), allocatable, dimension(:)   :: b9bn !< Vector with the source of the numerical scheme for the north boundary
   real(8), allocatable, dimension(:)   :: b9bs !< Vector with the source of the numerical scheme for the south boundary
   real(8), allocatable, dimension(:)   :: b9be !< Vector with the source of the numerical scheme for the east boundary
   real(8), allocatable, dimension(:)   :: b9bw !< Vector with the source of the numerical scheme for the west boundary

   real(8), allocatable, dimension(:,:) :: a5bn !< Matrix with the coefficients of the numerical scheme for the north boundary
   real(8), allocatable, dimension(:,:) :: a5bs !< Matrix with the coefficients of the numerical scheme for the south boundary
   real(8), allocatable, dimension(:,:) :: a5be !< Matrix with the coefficients of the numerical scheme for the east boundary
   real(8), allocatable, dimension(:,:) :: a5bw !< Matrix with the coefficients of the numerical scheme for the west boundary
   real(8), allocatable, dimension(:)   :: b5bn !< Vector with the source of the numerical scheme for the north boundary
   real(8), allocatable, dimension(:)   :: b5bs !< Vector with the source of the numerical scheme for the south boundary
   real(8), allocatable, dimension(:)   :: b5be !< Vector with the source of the numerical scheme for the east boundary
   real(8), allocatable, dimension(:)   :: b5bw !< Vector with the source of the numerical scheme for the west boundary
   real(8), allocatable, dimension(:)   :: fbe  !< Face of boundary east (1 if an east boundary, 0 otherwise)
   real(8), allocatable, dimension(:)   :: fbn  !< Face of boundary north (1 if a north boundary, 0 otherwise)


   !
   ! GEOMETRIC PARAMETERS
   !
   integer :: iocs ! index of the matching point between the ogive and the cylinder
   integer :: kfc  ! Kind of foredrag calculation ( 0 = over the whole forebody; 1 = over the ogive only)
   real(8) :: lr   ! length of the rocket
   real(8) :: rb   ! base radius of the rocket
   character(len=200) :: fgeom !< File of the geometric parameters

   !
   ! FILE ID NUMBERS
   !
   integer, parameter :: tfid = 10  ! Temporary file id
   integer, parameter :: lid  = 100 ! Listing file id
   integer, parameter :: rid  = 101 ! Residual file id


   character (len=20) :: date ! System date
   character (len=20) :: time ! System time
   character (len = 100) :: sim_id ! Simulation identification
   character (len = 100) :: input_file_parameters ! Input parameters data file

   !
   ! GAS PROPERTIES AND VARIABLES
   !
   integer :: ktm    ! Kind of thermophysical model ( 0 = constant, 1 = T dependent )
   integer :: modvis ! Viscosity model (0=Euler, 1=NS)

   real(8) ::  CPF !< Free stream Cp (J/kg.K)
   real(8) ::  VLF !< Free stream viscosity (Pa.s)
   real(8) ::  KPF !< Free stream thermal conductivity (W/m.K)
   real(8) ::   GF !< Free stream GF = gamma = Cp / Cv
   real(8) ::   PF !< Free stream pressure (Pa)
   real(8) ::   TF !< Free stream temperature (K)
   real(8) ::   HF !< Free stream total enthalpy (m2/s2)
   real(8) ::   MF !< Free stream Mach number
   real(8) ::   UF !< Free stream speed (m/s)
   real(8) ::  PRF !< Free stream Prandtl number
   real(8) ::  ROF !< Free stream density (kg/m3)
   real(8) :: REFm !< Free stream Reynolds number per meter (1/m)

   real(8) :: p_avg ! Average value of the pressure (Pa)

   real(8) :: Cdfi ! Foredrag coefficient due pressure (dimensionless)
   real(8) :: Cdfv ! Foredrag coefficient due viscosity (dimensionless)

   real(8) :: Tsbc !< Temperature on the south boundary (K) (if negative, adiabatic bc is applied)

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

   real(8), allocatable, dimension(:) :: ube !< u over the faces of the east boundary (m/s)
   real(8), allocatable, dimension(:) :: ubw !< u over the faces of the west boundary (m/s)
   real(8), allocatable, dimension(:) :: ubn !< u over the faces of the north boundary (m/s)
   real(8), allocatable, dimension(:) :: ubs !< u over the faces of the south boundary (m/s)

   real(8), allocatable, dimension(:) :: vbe !< v over the faces of the east boundary (m/s)
   real(8), allocatable, dimension(:) :: vbw !< v over the faces of the west boundary (m/s)
   real(8), allocatable, dimension(:) :: vbn !< v over the faces of the north boundary (m/s)
   real(8), allocatable, dimension(:) :: vbs !< v over the faces of the south boundary (m/s)

   real(8), allocatable, dimension(:) :: Tbe !< T over the faces of the east boundary (K)
   real(8), allocatable, dimension(:) :: Tbw !< T over the faces of the west boundary (K)
   real(8), allocatable, dimension(:) :: Tbn !< T over the faces of the north boundary (K)
   real(8), allocatable, dimension(:) :: Tbs !< T over the faces of the south boundary (K)

   real(8), allocatable, dimension(:) :: pbe !< p over the faces of the east boundary (Pa)
   real(8), allocatable, dimension(:) :: pbw !< p over the faces of the west boundary (Pa)
   real(8), allocatable, dimension(:) :: pbn !< p over the faces of the north boundary (Pa)
   real(8), allocatable, dimension(:) :: pbs !< p over the faces of the south boundary (Pa)


   real(8), allocatable, dimension(:) :: Ucbe !< Uc over the faces of the east  boundary (m2/s)
   real(8), allocatable, dimension(:) :: Ucbw !< Uc over the faces of the west  boundary (m2/s)
   real(8), allocatable, dimension(:) :: Ucbn !< Uc over the faces of the north boundary (m2/s)
   real(8), allocatable, dimension(:) :: Ucbs !< Uc over the faces of the south boundary (m2/s)

   real(8), allocatable, dimension(:) :: Vcbe !< Vc over the faces of the east  boundary (m2/s)
   real(8), allocatable, dimension(:) :: Vcbw !< Vc over the faces of the west  boundary (m2/s)
   real(8), allocatable, dimension(:) :: Vcbn !< Vc over the faces of the north boundary (m2/s)
   real(8), allocatable, dimension(:) :: Vcbs !< Vc over the faces of the south boundary (m2/s)

contains

   subroutine get_parameters
      implicit none
      !
      call date_time(date, time)
      !
      open(10,file="./mach2d_input/mach2d_input.txt")
      ! Reading file name from which parameters will be read
      read(10,*) input_file_parameters
      close(10)

      open(10, file = "./mach2d_input/" // input_file_parameters)
      read(10,*) sim_id ! Simulation identification  (up to 100 characters)
      read(10,*) nxi    ! Number of real volumes in the csi direction of the coarsest grid
      read(10,*) nyi    ! Number of real volumes in the eta direction of the coarsest grid
      read(10,*) nmf    ! Number of the finest mesh  (1<=nmf)
      read(10,*) nmd    ! Number of the desired mesh (1<=nmd<=nmf)
      read(10,*) fgeom  ! File of the geometric parameters
      read(10,*) kg     ! Kind of grid (1=uniform, 2=geometric progression, 3=power law, 4=gp modified, 5=hyperbolic)
      read(10,*) avi    ! Initial value of the artificial viscosity (only for kg=5)
      read(10,*) avf    ! Final value of the artificial viscosity (only for kg=5)
      read(10,*) awf    ! Area weighting factor (only for kg=5)
      read(10,*) kcm    ! Kind of centroid mean (1=simple mean, 2=weighted mean)
      read(10,*) coord  ! Kind of coord. system ( 1=cylindrical, 0 = cartesian)
      read(10,*) cbl    ! The width of the vol. closer to the wall is 'cbl' times the width of the b. layer
      read(10,*) itmax  ! Maximum number of iterations for time cycle
      read(10,*) itmmax ! Maximum number of iterations for mass cycle
      read(10,*) itpmax ! Maximum number of iteractions for pressure cycle
      read(10,*) itemax ! Maximum number of iteractions for extrapolation to fictitious
      read(10,*) nitm_u ! Maximum number of iteractions for solving the linear systems for u, v and T
      read(10,*) nitm_p ! Maximum number of iteractions for solving the linear system for p
      read(10,*) tol_u  ! Tolerance in the MSI for solving the linear systems for u, v and T
      read(10,*) tol_p  ! Tolerance in the MSI for solving the linear system for p
      read(10,*) tolm   ! Tolerance for the mass cycle
      read(10,*) tolt   ! Tolerance for the time evolution cycle
      read(10,*) wlf    ! Frequency of printing in the listing file
      read(10,*) sem_a  ! 1 = do not open result files, 0 = open
      read(10,*) sem_g  ! 0 = visualize the plot, 1 = do not visualize
      read(10,*) w_g    ! Frequency of writing data for graphics
      read(10,*) w_cam  ! 1 = write the fields, 0 = do not
      read(10,*) wppd   ! Write post processed data (0=no, 1=yes, 2=yes-simplified)
      read(10,*) beta   ! UDS/CDS mixing constant (0=UDS, 1=CDS)
      read(10,*) dt1    ! initial time step (s)
      read(10,*) dt2    ! final time step (s)
      read(10,*) it1    ! number of iteractions up to which dt = dt1
      read(10,*) it2    ! number of iteractions from which dt = dt2
      read(10,*) h0     ! Amplitude of h in the TSI11 model
      read(10,*) mincc  ! Minimum allowed value of the convergence coefficient
      read(10,*) maxcc  ! Maximum allowed value of the convergence coefficient
      read(10,*) modvis ! Viscosity model (0=Euler, 1=NS)
      read(10,*) ktm    ! Kind of thermophysical model ( 0 = constant, 1 = T dependent )
      read(10,*) kfc    ! Kind of foredrag calculation ( 0 = over the whole forebody; 1 = over the ogive only)
      read(10,*) Tsbc   ! Temperature on the south boundary (K) (if negative, adiabatic bc is applied)
      read(10,*) PF     ! Far field pressure (Pa)
      read(10,*) TF     ! Far field temperature (K)
      read(10,*) MF     ! Mach number of the free stream

      close(10)

      ! After reading the values of nxi and nyi, these variables must be
      ! changed to take into accout the fictitious volumes

      nxi = nxi + 2
      nyi = nyi + 2


      ! Calculating the number of volumes (real+fictitious) of the finest mesh

      nxf = (nxi-2) * 2 ** (nmf-1) + 2
      nyf = (nyi-2) * 2 ** (nmf-1) + 2



      ! Calculating the number of volumes (real+fictitious) of the desired mesh

      nx = (nxi-2) * 2 ** (nmd-1) + 2
      ny = (nyi-2) * 2 ** (nmd-1) + 2

      nxy = nx * ny

   end subroutine get_parameters


   subroutine write_parameters(fid)
      implicit none
      integer, intent(in) :: fid

  100 format( "'", A21, "'",' ....: ',A)
  101 format( "'", A  , "'",' ....: ',A)

      write(fid,*)
      write(fid,*) "Date: ", date
      write(fid,*) "Time: ", time
      write(fid,*)
      write(fid,*) "          PARAMETERS         "
      write(fid,*)
      write(fid,                100) trim(adjustl(sim_id)) , " Simulation identification  (up to 100 characters)"
      write(fid,"(I23,' ....: ',A)")   nmf     , " nmf    - Number of the finest mesh  (1<=nmf)"
      write(fid,"(I23,' ....: ',A)")   nmd     , " nmd    - Number of the desired mesh (1<=nmd<=nmf)"
      write(fid,"(I23,' ....: ',A)")    nx     , " nx     - Number of real+ficititious volumes in the csi direction (desired grid)"
      write(fid,"(I23,' ....: ',A)")    ny     , " ny     - Number of real+ficititious volumes in the eta direction (desired grid)"
      write(fid,"(I23,' ....: ',A)")   nxi     , " nxi    - Number of real+ficititious volumes in the csi direction (coarsest grid)"
      write(fid,"(I23,' ....: ',A)")   nyi     , " nyi    - Number of real+ficititious volumes in the eta direction (coarsest grid)"
      write(fid,"(I23,' ....: ',A)")   nxf     , " nxf    - Number of real+ficititious volumes in the csi direction (finest grid)"
      write(fid,"(I23,' ....: ',A)")   nyf     , " nyf    - Number of real+ficititious volumes in the eta direction (finest grid)"
      write(fid,                101) trim(adjustl(fgeom))  , " fgeom  - File of the geometric parameters"
      write(fid,"(I23,' ....: ',A)")    kg     , " kg     - Kind of grid (1=uniform, 2=geometric progression, 3=power law," &
      // " 4=gp modified)"
      write(fid,"(ES23.16,' ....: ',A)") avi   , " Initial value of the artificial viscosity (only for kg=5)"
      write(fid,"(ES23.16,' ....: ',A)") avf   , " Final value of the artificial viscosity (only for kg=5)"
      write(fid,"(ES23.16,' ....: ',A)") awf   , " Area weighting factor (only for kg=5)"
      write(fid,"(I23,' ....: ',A)")    kcm    , " kcm    - Kind of centroid mean (1=simple mean, 2=weighted mean)"
      write(fid,"(I23,' ....: ',A)")    coord  , " coord  - Kind of coord. system ( 1=cylindrical, 0 = cartesian)"
      write(fid,"(ES23.16,' ....: ',A)") cbl    , " cbl    - The width of the vol. closer to the wall is 'cbl' times the " &
         // "width of the b. layer"
      write(fid,"(I23,' ....: ',A)")    itmax  , " itmax  - Maximum number of iterations for time cycle"
      write(fid,"(I23,' ....: ',A)")    itmmax , " itmmax - Maximum number of iterations for mass cycle"
      write(fid,"(I23,' ....: ',A)")    itpmax , " itpmax - Maximum number of iteractions for pressure cycle"
      write(fid,"(I23,' ....: ',A)")    itemax , " itemax - Maximum number of iteractions for extrapolation to fictitious"
      write(fid,"(I23,' ....: ',A)")    nitm_u &
      , " nitm_u - Maximum number of iteractions for solving the linear systems for u, v and T"
      write(fid,"(I23,' ....: ',A)")    nitm_p , " nitm_p - Maximum number of iteractions for solving the linear system for p"
      write(fid,"(ES23.16,' ....: ',A)") tol_u , " tol_u  - Tolerance in the MSI for solving the linear systems for u, v and T"
      write(fid,"(ES23.16,' ....: ',A)") tol_p , " tol_p  - Tolerance in the MSI for solving the linear system for p"
      write(fid,"(ES23.16,' ....: ',A)") tolm  , " tolm   - Tolerance for the mass cycle"
      write(fid,"(ES23.16,' ....: ',A)") tolt  , " tolt   - Tolerance for the time evolution cycle"
      write(fid,"(I23,' ....: ',A)")    wlf    , " wlf    - Frequency of printing in the listing file"
      write(fid,"(I23,' ....: ',A)")    sem_a  , " sem_a  - 1 = do not open result files, 0 = open"
      write(fid,"(I23,' ....: ',A)")    sem_g  , " sem_g  - 0 = visualize the plot, 1 = do not visualize"
      write(fid,"(I23,' ....: ',A)")    w_g    , " w_g    - Frequency of writing data for graphics"
      write(fid,"(I23,' ....: ',A)")    w_cam  , " w_cam  - 1 = write the fields, 0 = do not"
      write(fid,"(I23,' ....: ',A)")    wppd   , " wppd   - Write post processed data (0=no, 1=yes, 2=yes-simplified)"
      write(fid,"(ES23.16,' ....: ',A)") beta  , " beta   - UDS/CDS mixing constant (0=UDS, 1=CDS)"
      write(fid,"(ES23.16,' ....: ',A)") dt1   , " dt1    - initial time step (s)"
      write(fid,"(ES23.16,' ....: ',A)") dt2   , " dt2    - final time step (s)"
      write(fid,"(I23,' ....: ',A)")    it1    , " it1    - number of iteractions up to which dt = dt1"
      write(fid,"(I23,' ....: ',A)")    it2    , " it2    - number of iteractions from which dt = dt2"
      write(fid,"(ES23.16,' ....: ',A)") h0    , " h0     - Amplitude of h in the TSI11 model"
      write(fid,"(ES23.16,' ....: ',A)") mincc , " mincc  - Minimum allowed value of the convergence coefficient"
      write(fid,"(ES23.16,' ....: ',A)") maxcc , " maxcc  - Maximum allowed value of the convergence coefficient"
      write(fid,"(I23,' ....: ',A)")    modvis , " modvis - Viscosity model (0=Euler, 1=NS)"
      write(fid,"(I23,' ....: ',A)")    ktm    , " ktm    - Kind of thermophysical model ( 0 = constant, 1 = T dependent )"
      write(fid,"(I23,' ....: ',A)")    kfc    , " kfc    - Kind of foredrag calculation " &
         // "( 0 = over the whole forebody; 1 = over the ogive only)"
      write(fid,"(ES23.16,' ....: ',A)") Tsbc  , " Tsbc   - Temperature on the"&
         // " south boundary (K) (if negative, adiabatic bc is applied)"
      write(fid,"(ES23.16,' ....: ',A)") PF    , " PF     - Far field pressure (Pa)"
      write(fid,"(ES23.16,' ....: ',A)") TF    , " TF     - Far field temperature (K)"
      write(fid,"(ES23.16,' ....: ',A)") MF    , " MF     - Mach number of the free stream"
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


      !allocate( ccu(nxy), ccv(nxy), cct(nxy), ccp(nxy) )

      allocate( xf(nxf*nyf), yf(nxf*nyf) )

      allocate( ube(ny), ubw(ny), ubn(nx), ubs(nx) )
      allocate( vbe(ny), vbw(ny), vbn(nx), vbs(nx) )
      allocate( Tbe(ny), Tbw(ny), Tbn(nx), Tbs(nx) )
      allocate( pbe(ny), pbw(ny), pbn(nx), pbs(nx) )


      allocate( Ucbe (ny), Ucbw (ny), Ucbn (nx), Ucbs (nx) )
      allocate( Vcbe (ny), Vcbw (ny), Vcbn (nx), Vcbs (nx) )


      allocate( a9bn(nx,9), a9bs(nx,9), a9be(ny,9), a9bw(ny,9) )

      allocate( b9bn(nx), b9bs(nx), b9be(ny), b9bw(ny) )


      allocate( a5bn(nx,5), a5bs(nx,5), a5be(ny,5), a5bw(ny,5) )

      allocate( b5bn(nx), b5bs(nx), b5be(ny), b5bw(ny) )


      allocate( au(nxy,9), av(nxy,9), at(nxy,9), ap(nxy,5) )

      allocate( dl9(nxy,5), du9(nxy,4), dl5(nxy,4), du5(nxy,3) )

      allocate( fbe(nxy), fbn(nxy) )


   end subroutine allocate_variables


   subroutine initialize_variables
      implicit none

      ! Reading data from input data file

      call get_parameters

      ! Openning listing file

      open(lid, file = "./mach2d_output/mach2d_" // trim(adjustl(sim_id)) // ".txt")

      write(lid,*)
      write(lid,*) "LISTING FILE OF MACH2D"
      write(lid,*)

      call write_parameters(lid)

      call allocate_variables

      ! Initializing variables

      x   = 0.d0
      y   = 0.d0
      xf  = 0.d0
      yf  = 0.d0
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

      ube = 0.d0
      ubw = 0.d0
      ubn = 0.d0
      ubs = 0.d0

      vbe = 0.d0
      vbw = 0.d0
      vbn = 0.d0
      vbs = 0.d0

      Tbe = 0.d0
      Tbw = 0.d0
      Tbn = 0.d0
      Tbs = 0.d0

      pbe = 0.d0
      pbw = 0.d0
      pbn = 0.d0
      pbs = 0.d0

      Ucbe = 0.d0
      Ucbw = 0.d0
      Ucbn = 0.d0
      Ucbs = 0.d0

      Vcbe = 0.d0
      Vcbw = 0.d0
      Vcbn = 0.d0
      Vcbs = 0.d0

      ap = 0.d0
      bp = 0.d0
      pl = 0.d0

      pa = 0.d0
      Ta = 0.d0
      at = 0.d0
      bt = 0.d0

      !ccu = 0.d0
      !ccv = 0.d0
      !cct = 0.d0
      !ccp = 0.d0

      dl9 = 0.d0
      du9 = 0.d0
      dl5 = 0.d0
      du5 = 0.d0

      a9bn = 0.d0
      a9bs = 0.d0
      a9be = 0.d0
      a9bw = 0.d0

      b9bn = 0.d0
      b9bs = 0.d0
      b9be = 0.d0
      b9bw = 0.d0


      a5bn = 0.d0
      a5bs = 0.d0
      a5be = 0.d0
      a5bw = 0.d0

      b5bn = 0.d0
      b5bs = 0.d0
      b5be = 0.d0
      b5bw = 0.d0

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

end module data
