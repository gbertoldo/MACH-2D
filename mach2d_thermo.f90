
!> \brief Calculates the thermophysical properties of a pure gas
!! or of a mixture of gaseous species.

module thermophysical

   implicit none

   ! Public variables

   integer :: ngs  !< Number of gaseous species
   integer :: base !< Fraction base ( 0 = molar, 1 = massic )

   real(8) :: Ru   !< Universal gas constant (J/mol.K)
   real(8) :: Rg   !< Gas constant of the mixture (J/kg.K)
   real(8) :: Mg   !< Molar mass of the mixture (kg/mol)

   real(8), allocatable :: Xi(:)  !< Molar fraction
   real(8), allocatable :: Yi(:)  !< Massic fraction
   real(8), allocatable :: Mi(:)  !< Molar mass of the gaseous specie i (kg/mol)
   real(8), allocatable :: Ri(:)  !< Gas constant of the gaseous specie i (J/kg.K)

   real(8), allocatable :: cpi(:,:) !< Coefficients for the calculation of cp of specie i
   real(8), allocatable :: mui(:,:) !< Coefficients for the calculation of mu of specie i
   real(8), allocatable :: kpi(:,:) !< Coefficients for the calculation of kp of specie i

   character(10), allocatable :: NGi(:)  !< Name of the gaseous specie i (J/kg.K)



   ! Private variables

   character(300), private :: file_input      = "./mach2d_input/thermo_input.txt" !< Path and file name of the input files
   character(300), private :: file_parameters !< Path and file name of the input parameters
   character(300), private :: file_database   !< Path and file name of the database



contains


   !> \brief Initializes variables related to the thermophysical properties of the gas
   subroutine thermophysical_initialization!(ngs, base, Ru, Rg, Mg, Xi, Yi, Mi, Ri, cpi, mui, kpi, NGi) ! Output: all
      implicit none

      !integer, intent(out) :: ngs  !< Number of gaseous species
      !integer, intent(out) :: base !< Fraction base ( 0 = molar, 1 = massic )

      !real(8), intent(out) :: Ru   !< Universal gas constant (J/K.mol)
      !real(8), intent(out) :: Rg   !< Gas constant of the mixture (J/kg.K)
      !real(8), intent(out) :: Mg   !< Molar mass of the mixture (kg/mol)

      !real(8), allocatable, intent(out) :: Xi(:)  !< Molar fraction
      !real(8), allocatable, intent(out) :: Yi(:)  !< Massic fraction
      !real(8), allocatable, intent(out) :: Mi(:)  !< Molar mass of the gaseous specie i (kg/mol)
      !real(8), allocatable, intent(out) :: Ri(:)  !< Gas constant of the gaseous specie i (J/kg.K)

      !real(8), allocatable, intent(out) :: cpi(:,:) !< Coefficients for the calculation of cp of specie i
      !real(8), allocatable, intent(out) :: mui(:,:) !< Coefficients for the calculation of mu of specie i
      !real(8), allocatable, intent(out) :: kpi(:,:) !< Coefficients for the calculation of kp of specie i

      !character(10), allocatable, intent(out) :: NGi(:)  !< Name of the gaseous specie i


      ! Inner variables

      integer       ::    i ! Dummy index
      real(8)       :: raux ! Auxiliary variable
      character(10) :: caux ! Auxiliary variable



      ! Reading the path and file name of the parameters and database
      open(10, file = file_input )

      read(10,"(A)") file_parameters
      read(10,"(A)") file_database

      close(10)


      ! Universal gas constant (J/K.mol)
      Ru = 8.31451000D+00


      ! Reading the number of gaseous species in the mixture

      ngs = 0

      open(10, file = file_parameters )

      read(10,*) base

      do

         read(10,*) raux, caux

         if ( trim(adjustl(caux)) == "END" ) exit

         if ( raux > 1.d-8 ) ngs = ngs + 1

      end do


      ! Allocating memory

      allocate(  Xi(ngs)    )
      allocate(  Yi(ngs)    )
      allocate(  Mi(ngs)    )
      allocate(  Ri(ngs)    )
      allocate( NGi(ngs)    )
      allocate( cpi(ngs,10) )
      allocate( mui(ngs, 8) )
      allocate( kpi(ngs, 8) )


      ! Reading the fractions for gasous mixture

      rewind(10)

      read(10,*) base

      if ( base == 0 ) then ! Molar fraction

         i = 0

         do

            read(10,*) raux, caux

            if ( trim(adjustl(caux)) == "END" ) exit

            if ( raux > 1.d-8 ) then

               i = i + 1

               Xi(i) = raux

               NGi(i) = caux

            end if

         end do

         ! Checking consistency of input data

         if ( abs( 1.d0 - sum(Xi) ) > epsilon(1.d0) * 10.d0 ) then

            write(*,*) "ERROR: The sum of the molar fraction is not equal 1."

            stop

         end if


      else ! Massic fraction

         i = 0

         do

            read(10,*) raux, caux

            if ( trim(adjustl(caux)) == "END" ) exit

            if ( raux > 1.d-8 ) then

               i = i + 1

               Yi(i) = raux

               NGi(i) = caux

            end if

         end do

         ! Checking consistency of input data

         if ( abs( 1.d0 - sum(Yi) ) > epsilon(1.d0) * 10.d0 ) then

            write(*,*) "ERROR: The sum of the massic fraction is not equal 1."

            stop

         end if

      end if

      close(10)


      ! Reading database

      do i = 1, ngs

         call get_specie_data( NGi(i), Mi(i), cpi(i,:), mui(i,:), kpi(i,:) )

      end do


      ! Calculating extra data

      if ( base == 0 ) then

         Mg = dot_product(Xi,Mi)

         Yi = Xi * Mi / Mg

      else

         Xi = ( Yi / Mi ) / sum( Yi / Mi )

         Mg = dot_product(Xi,Mi)

      end if


      ! Gas constant of each specie
      Ri = Ru / Mi


      ! Gas constant of the mixture
      Rg = Ru / Mg


   end subroutine


   !> \brief Extracts data from database for a given gas specie
   subroutine get_specie_data(NGi, Mi, cpi, mui, kpi)
      implicit none
      character(10), intent(in)  ::     NGi !< Name of the gaseous specie i
      real(8),       intent(out) ::      Mi !< Molar mass of specie i (kg/mol)
      real(8),       intent(out) :: cpi(10) !< Coefficients for the calculation of cp of specie i
      real(8),       intent(out) :: mui( 8) !< Coefficients for the calculation of mu of specie i
      real(8),       intent(out) :: kpi( 8) !< Coefficients for the calculation of kp of specie i


      ! Inner variables

      integer :: io ! Auxiliary variable

      character(10) :: caux1 ! Auxiliary variable
      character(10) :: caux2 ! Auxiliary variable

      integer :: lck_NGi ! Locker for NGi
      integer :: lck_Mi  ! Locker for Mi
      integer :: lck_cpi ! Locker for cpi
      integer :: lck_mui ! Locker for mui
      integer :: lck_kpi ! Locker for kpi


      ! Initializing lockers

      lck_NGi = 0
      lck_Mi  = 0
      lck_cpi = 0
      lck_mui = 0
      lck_kpi = 0


      open(10, file = file_database )


      ! Swept the file
      do

         read(10,*, iostat = io ) caux1

         if ( io /= 0 ) exit


         ! Searches for the ID_NAME of the gaseous specie
         if ( trim(caux1) == "ID_NAME" ) then

            backspace(10)

            read(10,*, iostat = io ) caux1, caux2

            if ( io /= 0 ) exit


            ! If the searched species matches the found one, then
            if ( trim(caux2) == trim(NGi) ) then

               lck_NGi = 1


               ! Swept the entry searching for necessary data
               do

                  read(10, *, iostat = io ) caux1

                  if ( io /= 0 .or. trim(caux1) == "END" ) exit



                  select case (trim(caux1))

                     ! Molar mass
                     case ("ID_M")

                        backspace(10)

                        read(10,*) caux1, Mi

                        lck_Mi = 1


                     ! Coefficients for Cp on the lower temperature range
                     case ("ID_CPL")

                        backspace(10)

                        read(10,*) caux1, Cpi(1:5)

                        lck_cpi = lck_cpi + 1



                     ! Coefficients for Cp on the upper temperature range
                     case ("ID_CPU")

                        backspace(10)

                        read(10,*) caux1, Cpi(6:10)

                        lck_cpi = lck_cpi + 1



                     ! Coefficients for viscosity on the lower temperature range
                     case ("ID_VL")

                        backspace(10)

                        read(10,*) caux1, mui(1:4)

                        lck_mui = lck_mui + 1



                     ! Coefficients for viscosity on the upper temperature range
                     case ("ID_VU")

                        backspace(10)

                        read(10,*) caux1, mui(5:8)

                        lck_mui = lck_mui + 1



                     ! Coefficients for thermal conductivity on the lower temperature range
                     case ("ID_KL")

                        backspace(10)

                        read(10,*) caux1, kpi(1:4)

                        lck_kpi = lck_kpi + 1



                     ! Coefficients for thermal conductivity on the upper temperature range
                     case ("ID_KU")

                        backspace(10)

                        read(10,*) caux1, kpi(5:8)

                        lck_kpi = lck_kpi + 1



                     case default

                  end select


               end do

            end if

         end if

      end do

      close(10)


      ! Checking data

      if ( lck_NGi /= 1 ) then

         write(*,*) "ERROR: Specie ", trim(NGi), " was not found in the database."

         stop

      end if


      if ( lck_Mi /= 1 ) then

         write(*,*) "ERROR: Molar mass of ", trim(NGi), " was not found in the database."

         stop

      end if


      if ( lck_cpi /= 2 ) then

         write(*,*) "ERROR: Coef. for Cp of ", trim(NGi), " was not found in the database."

         stop

      end if

      if ( lck_mui /= 2 ) then

         write(*,*) "ERROR: Coef. for mu of ", trim(NGi), " was not found in the database."

         stop

      end if

      if ( lck_kpi /= 2 ) then

         write(*,*) "ERROR: Coef. for kp of ", trim(NGi), " was not found in the database."

         stop

      end if

   end subroutine


   !> \brief Specific heat at constant pressure of the mixture (J/kg.K)
   real(8) function cp_mix(ngs, T, Yi, Ri, cpi)
      implicit none
      integer, intent(in) :: ngs         !< Number of gaseous species
      real(8), intent(in) :: T           !< Temperature (K)
      real(8), intent(in) :: Yi(ngs)     !< Massic fraction
      real(8), intent(in) :: Ri(ngs)     !< Gas constant of the gaseous specie i (J/kg.K)
      real(8), intent(in) :: cpi(ngs,10) !< Coefficients for the calculation of cp of specie i

      ! Inner variables

      integer :: i ! Dummy index

      cp_mix = 0.d0

      do i = 1, ngs

         cp_mix = cp_mix + Yi(i) * cp_i(T, Ri(i), cpi(i,:))

      end do

   end function



   !> \brief Specific heat at constant pressure of specie i (J/kg.K)
   real(8) function cp_i(T, Ri, ai)
      implicit none
      real(8), intent(in) :: T      !< Temperature (K)
      real(8), intent(in) :: Ri     !< Gas constant (J/kg.K)
      real(8), intent(in) :: ai(10) !< Coefficients for the calculation of cp of specie i

      if ( T <= 1000.d0 ) then

         cp_i = ai(1) + ai(2) * T + ai(3) * T ** 2 + ai(4) * T ** 3 + ai(5) * T ** 4

      else

         cp_i = ai(6) + ai(7) * T + ai(8) * T ** 2 + ai(9) * T ** 3 + ai(10) * T ** 4

      end if

      cp_i = cp_i * Ri

   end function


   !> \brief Viscosity of the gas mixture (Pa.s)
   real(8) function mu_mix(ngs, T, Xi, Mi, mui)
      implicit none
      integer, intent(in) :: ngs        !< Number of gaseous species
      real(8), intent(in) :: T          !< Temperature (K)
      real(8), intent(in) :: Xi(ngs)    !< Molar fraction
      real(8), intent(in) :: Mi(ngs)    !< Molar mass (kg/mol)
      real(8), intent(in) :: mui(ngs,8) !< Coefficients for the calculation of mu of specie i


      ! Inner variables

      integer :: i ! Dummy index
      integer :: j ! Dummy index

      real(8) :: raux1(ngs) ! Auxiliary variable
      real(8) :: raux2      ! Auxiliary variable


      ! Calculating the viscosity of each specie
      do i = 1, ngs

         raux1(i) = mu_i(T,mui(i,:))

      end do


      ! Calculating the viscosity of the mixture

      mu_mix = 0.d0

      do i = 1, ngs

         raux2 = 0.d0

         do j = 1, ngs

            raux2 = raux2 &

               + Xi(j)    &

               * (1.d0 + sqrt( raux1(i) / raux1(j) * sqrt( Mi(j) / Mi(i) ) ) ) ** 2  &

               / sqrt( 8.d0 * (1.d0 + Mi(i) / Mi(j) ) )

         end do

         mu_mix = mu_mix + Xi(i) * raux1(i) / raux2

      end do

   end function


   !> \brief Viscosity of specie i (Pa.s)
   real(8) function mu_i(T, ai)
      implicit none
      real(8), intent(in) :: T     !< Temperature (K)
      real(8), intent(in) :: ai(8) !< Coefficients for the calculation of mu of specie i

      if ( T <= 1000.d0 ) then

         mu_i = exp( ai(1) * log(T) + ai(2) / T + ai(3) / T ** 2 + ai(4) ) * 1.d-7

      else

         mu_i = exp( ai(5) * log(T) + ai(6) / T + ai(7) / T ** 2 + ai(8) ) * 1.d-7

      end if

   end function


   !> \brief Thermal conductivity of the gas mixture (W/m.K)
   real(8) function kp_mix(ngs, T, Xi, Mi, kpi)
      implicit none
      integer, intent(in) :: ngs        !< Number of gaseous species
      real(8), intent(in) :: T          !< Temperature (K)
      real(8), intent(in) :: Xi(ngs)    !< Molar fraction
      real(8), intent(in) :: Mi(ngs)    !< Molar mass (kg/mol)
      real(8), intent(in) :: kpi(ngs,8) !< Coefficients for the calculation of kp of specie i


      ! Inner variables

      integer :: i ! Dummy index
      integer :: j ! Dummy index

      real(8) :: raux1(ngs) ! Auxiliary variable
      real(8) :: raux2      ! Auxiliary variable


      ! Calculating the thermal conductivity of each specie
      do i = 1, ngs

         raux1(i) = kp_i(T,kpi(i,:))

      end do


      ! Calculating the thermal conductivity of the mixture

      kp_mix = 0.d0

      do i = 1, ngs

         raux2 = 0.d0

         do j = 1, ngs

            raux2 = raux2 &

               + Xi(j)    &

               * (1.d0 + sqrt( raux1(i) / raux1(j) * sqrt( Mi(j) / Mi(i) ) ) ) ** 2  &

               / sqrt( 8.d0 * (1.d0 + Mi(i) / Mi(j) ) )

         end do

         kp_mix = kp_mix + Xi(i) * raux1(i) / raux2

      end do

   end function


   !> \brief Thermal conductivity of specie i (W/m.K)
   real(8) function kp_i(T, ai)
      implicit none
      real(8), intent(in) :: T     !< Temperature (K)
      real(8), intent(in) :: ai(8) !< Coefficients for the calculation of kp of specie i

      if ( T <= 1000.d0 ) then

         kp_i = exp( ai(1) * log(T) + ai(2) / T + ai(3) / T ** 2 + ai(4) ) * 1.d-4

      else

         kp_i = exp( ai(5) * log(T) + ai(6) / T + ai(7) / T ** 2 + ai(8) ) * 1.d-4

      end if

   end function


   !> \brief Calculates cp at the center of each real volume and over the
   !! boundaries. The boundary values are stored in the fictitious volumes.
   subroutine set_cp( nx, ny, ngs, Yi, Ri, cpi, Tbn, Tbs, Tbe, Tbw, T, cp) ! Output: last one
      implicit none
      integer, intent(in) :: nx  !< Number of volumes in the csi direction (real+fictitious)
      integer, intent(in) :: ny  !< Number of volumes in the eta direction (real+fictitious)
      integer, intent(in) :: ngs !< Number of gaseous species
      real(8), intent(in) :: Yi(ngs)     !< Massic fraction
      real(8), intent(in) :: Ri(ngs)     !< Gas constant of the gaseous specie i (J/kg.K)
      real(8), intent(in) :: cpi(ngs,10) !< Coefficients for the calculation of cp of specie i
      real(8), dimension(nx),    intent(in)  :: Tbn !< Temperature over the north boundary (K)
      real(8), dimension(nx),    intent(in)  :: Tbs !< Temperature over the south boundary (K)
      real(8), dimension(ny),    intent(in)  :: Tbe !< Temperature over the  east boundary (K)
      real(8), dimension(ny),    intent(in)  :: Tbw !< Temperature over the  west boundary (K)
      real(8), dimension(nx*ny), intent(in)  :: T   !< Temperature at center of volume P (K)
      real(8), dimension(nx*ny), intent(out) :: cp  !< Specific heat at const pressure (J/(kg.K))

      ! Inner variables

      integer :: i  ! Dummy index
      integer :: j  ! Dummy index
      integer :: np, npn, nps, npe, npw, npsw, npse, npnw, npne ! Dummy index

      ! Real volumes

      do j = 2, ny-1

         do i = 2, nx-1

             np  = nx * (j-1) + i

             cp(np) = cp_mix(ngs, T(np), Yi, Ri, cpi)

         end do

      end do



      ! North boundary

      j = ny

      do i = 2, nx-1

         np  = nx * (j-1) + i

         cp(np) = cp_mix(ngs, Tbn(i), Yi, Ri, cpi)

      end do



      ! South boundary

      j = 1

      do i = 2, nx-1

         np  = nx * (j-1) + i

         cp(np) = cp_mix(ngs, Tbs(i), Yi, Ri, cpi)

      end do



      ! East boundary

      i = nx

      do j = 2, ny-1

         np  = nx * (j-1) + i

         cp(np) = cp_mix(ngs, Tbe(j), Yi, Ri, cpi)

      end do



      ! West boundary

      i = 1

      do j = 2, ny-1

         np  = nx * (j-1) + i

         cp(np) = cp_mix(ngs, Tbw(j), Yi, Ri, cpi)

      end do


      ! SW corner

      i = 1;    j = 1

      np  = nx * (j-1) + i
      npn  = np + nx
      npe  = np + 1
      npne = npn + 1

      cp(np) = ( cp(npn) + cp(npe) + cp(npne) ) / 3.d0


      ! SE corner

      i = nx;   j = 1

      np   = nx * (j-1) + i
      npn  = np + nx
      npw  = np - 1
      npnw = npn - 1

      cp(np) = ( cp(npn) + cp(npw) + cp(npnw) ) / 3.d0



      ! NW corner

      i = 1;  j = ny

      np   = nx * (j-1) + i
      nps  = np - nx
      npe  = np + 1
      npse = nps + 1

      cp(np) = ( cp(nps) + cp(npe) + cp(npse) ) / 3.d0



      ! NE corner

      i = nx;   j = ny

      np   = nx * (j-1) + i
      nps  = np - nx
      npw  = np - 1
      npsw = nps - 1

      cp(np) = ( cp(nps) + cp(npw) + cp(npsw) ) / 3.d0

   end subroutine set_cp

   !> \brief Calculates gamma = Cp / Cv at the center of each real volume and
   !! over the boundaries based on Cp and Rg. The boundary values are stored in
   !! the fictitious volumes.
   subroutine set_gamma( nx, ny, Rg, cp, gcp) ! Output: last one
      implicit none
      integer, intent(in) :: nx  !< Number of volumes in the csi direction (real+fictitious)
      integer, intent(in) :: ny  !< Number of volumes in the eta direction (real+fictitious)
      real(8), intent(in) :: Rg  !< Perfect gas constant (J/(kg.K))
      real(8), dimension(nx*ny), intent(in)  :: cp  !< Specific heat at const pressure (J/(kg.K))
      real(8), dimension(nx*ny), intent(out) :: gcp !< gcp = gamma = Cp/Cv at center of CV P (dimensionless)

      gcp = cp / ( cp - Rg ) ! dimensionless

   end subroutine set_gamma


   !> \brief Calculates the laminar viscosity at the nodes of real volumes
   subroutine set_laminar_viscosity_at_nodes(nx, ny, ngs, Xi, Mi, mui, T, vlp) ! Output: last one
      implicit none
      integer, intent(in) :: nx         !< Number of volumes in csi direction (real+fictitious)
      integer, intent(in) :: ny         !< Number of volumes in eta direction (real+fictitious)
      integer, intent(in) :: ngs        !< Number of gaseous species
      real(8), intent(in) :: Xi(ngs)    !< Molar fraction
      real(8), intent(in) :: Mi(ngs)    !< Molar mass (kg/mol)
      real(8), intent(in) :: mui(ngs,8) !< Coefficients for the calculation of mu of specie i
      real(8), dimension(nx*ny), intent(in)  :: T   !< Temperature of the last iteraction (K)
      real(8), dimension(nx*ny), intent(out) :: vlp !< Laminar viscosity at center of volume P (Pa.s)


      ! Inner variables

      integer :: i, j, np

      do j = 2, ny-1

         do i = 2, nx-1

            np = nx * (j-1) + i

            vlp(np) = mu_mix(ngs, T(np), Xi, Mi, mui)

         end do

      end do

   end subroutine set_laminar_viscosity_at_nodes


   !> \brief Calculates the thermal conductivity at the nodes of real volumes
   subroutine set_thermal_conductivity_at_nodes(nx, ny, ngs, Xi, Mi, kpi, T, kp) ! Output: last one
      implicit none
      integer, intent(in) :: nx         !< Number of volumes in csi direction (real+fictitious)
      integer, intent(in) :: ny         !< Number of volumes in eta direction (real+fictitious)
      integer, intent(in) :: ngs        !< Number of gaseous species
      real(8), intent(in) :: Xi(ngs)    !< Molar fraction
      real(8), intent(in) :: Mi(ngs)    !< Molar mass (kg/mol)
      real(8), intent(in) :: kpi(ngs,8) !< Coefficients for the calculation of kp of specie i
      real(8), dimension(nx*ny), intent(in)  :: T  !< Temperature of the last iteraction (K)
      real(8), dimension(nx*ny), intent(out) :: kp !< Thermal conductivity at center of volume P (W/m.K)

      ! Inner variables

      integer :: i, j, np


      do j = 2, ny-1

         do i = 2, nx-1

           np = nx * (j-1) + i

           kp(np) = kp_mix(ngs, T(np), Xi, Mi, kpi)

         end do

      end do

   end subroutine set_thermal_conductivity_at_nodes


   !> \brief Calculates the laminar viscosity at faces
   subroutine get_laminar_viscosity_at_faces(nx, ny, ngs, Xi, Mi, mui, Tbn &
      ,                     Tbs, Tbe, Tbw, vlp, vle, vln) ! Output: last two
      implicit none
      integer, intent(in) :: nx         !< Number of volumes in csi direction (real+fictitious)
      integer, intent(in) :: ny         !< Number of volumes in eta direction (real+fictitious)
      integer, intent(in) :: ngs        !< Number of gaseous species
      real(8), intent(in) :: Xi(ngs)    !< Molar fraction
      real(8), intent(in) :: Mi(ngs)    !< Molar mass (kg/mol)
      real(8), intent(in) :: mui(ngs,8) !< Coefficients for the calculation of mu of specie i
      real(8), dimension(nx),    intent(in)  :: Tbn !< Temperature over the north boundary (K)
      real(8), dimension(nx),    intent(in)  :: Tbs !< Temperature over the south boundary (K)
      real(8), dimension(ny),    intent(in)  :: Tbe !< Temperature over the  east boundary (K)
      real(8), dimension(ny),    intent(in)  :: Tbw !< Temperature over the  west boundary (K)
      real(8), dimension(nx*ny), intent(in)  :: vlp !< Laminar viscosity at center of volume P (Pa.s)
      real(8), dimension(nx*ny), intent(out) :: vle !< Laminar viscosity at center of face east (Pa.s)
      real(8), dimension(nx*ny), intent(out) :: vln !< Laminar viscosity at center of face north (Pa.s)


      ! Inner variables
      integer :: i, j, np, npe, npn


      ! Calculating the viscosity over the east faces of real volumes (excluding boundaries)

      do j = 2, ny-1

         do i = 2, nx-2

            np   = nx * (j-1) + i

            npe  = np + 1

            vle(np) = 2.d0 * vlp(np) * vlp(npe) / ( vlp(np) + vlp(npe) )

         end do

      end do

      ! Calculating the viscosity over the west boundary

      i = 1

      do j = 2, ny-1

         np = nx * (j-1) + i

         vle(np) = mu_mix(ngs, Tbw(j), Xi, Mi, mui)

      end do


      ! Calculating the viscosity over the east boundary

      i = nx-1

      do j = 2, ny-1

         np = nx * (j-1) + i

         vle(np) = mu_mix(ngs, Tbe(j), Xi, Mi, mui)

      end do


      ! Calculating the viscosity over the north faces of real volumes (excluding boundaries)

      do i = 2, nx-1

         do j = 2, ny-2

            np  = nx * (j-1) + i

            npn = np + nx

            vln(np) = 2.d0 * vlp(np) * vlp(npn) / ( vlp(np) + vlp(npn) )

         end do

      end do

      ! Calculating the viscosity over the south boundary

      j = 1

      do i = 2, nx-1

         np = nx * (j-1) + i

         vln(np) = mu_mix(ngs, Tbs(i), Xi, Mi, mui)

      end do


      ! Calculating the viscosity over the north boundary

      j = ny-1

      do i = 2, nx-1

         np = nx * (j-1) + i

         vln(np) = mu_mix(ngs, Tbn(i), Xi, Mi, mui)

      end do


   end subroutine get_laminar_viscosity_at_faces


   !> \brief Calculates the thermal conductivity at faces
   subroutine get_thermal_conductivity_at_faces(nx, ny, ngs, Xi, Mi, kpi, Tbn &
         ,                      Tbs, Tbe, Tbw, kp, ke, kn) ! Output: last two
      implicit none
      integer, intent(in) :: nx         !< Number of volumes in csi direction (real+fictitious)
      integer, intent(in) :: ny         !< Number of volumes in eta direction (real+fictitious)
      integer, intent(in) :: ngs        !< Number of gaseous species
      real(8), intent(in) :: Xi(ngs)    !< Molar fraction
      real(8), intent(in) :: Mi(ngs)    !< Molar mass (kg/mol)
      real(8), intent(in) :: kpi(ngs,8) !< Coefficients for the calculation of kp of specie i
      real(8), dimension(nx),    intent(in)  :: Tbn !< Temperature over the north boundary (K)
      real(8), dimension(nx),    intent(in)  :: Tbs !< Temperature over the south boundary (K)
      real(8), dimension(ny),    intent(in)  :: Tbe !< Temperature over the  east boundary (K)
      real(8), dimension(ny),    intent(in)  :: Tbw !< Temperature over the  west boundary (K)
      real(8), dimension(nx*ny), intent(in)  :: kp  !< Thermal conductivity at center of volume P (W/m.K)
      real(8), dimension(nx*ny), intent(out) :: ke  !< Thermal conductivity at center of face east (W/m.K)
      real(8), dimension(nx*ny), intent(out) :: kn  !< Thermal conductivity at center of face north (W/m.K)

      ! Inner variables
      integer :: i, j, np, npe, npn


      ! Calculating the thermal conductivity over the east faces of real volumes (excluding boundaries)

      do j = 2, ny-1

         do i = 2, nx-2

            np   = nx * (j-1) + i

            npe  = np + 1

            ke(np) = 2.d0 * kp(np) * kp(npe) / ( kp(np) + kp(npe) )

         end do

      end do


      ! Calculating the thermal conductivity over the west boundary

      i = 1

      do j = 2, ny-1

         np = nx * (j-1) + i

         ke(np) = kp_mix(ngs, Tbw(j), Xi, Mi, kpi)

      end do


      ! Calculating the thermal conductivity over the east boundary

      i = nx-1

      do j = 2, ny-1

         np = nx * (j-1) + i

         ke(np) = kp_mix(ngs, Tbe(j), Xi, Mi, kpi)

      end do


      ! Calculating the thermal conductivity over the north faces of real volumes (excluding boundaries)

      do i = 2, nx-1

         do j = 2, ny-2

            np  = nx * (j-1) + i

            npn = np + nx

            kn(np) = 2.d0 * kp(np) * kp(npn) / ( kp(np) + kp(npn) )

         end do

      end do


      ! Calculating the thermal conductivity over the south boundary

      j = 1

      do i = 2, nx-1

         np = nx * (j-1) + i

         kn(np) = kp_mix(ngs, Tbs(i), Xi, Mi, kpi)

      end do


      ! Calculating the thermal conductivity over the north boundary

      j = ny-1

      do i = 2, nx-1

         np = nx * (j-1) + i

         kn(np) = kp_mix(ngs, Tbn(i), Xi, Mi, kpi)

      end do

   end subroutine get_thermal_conductivity_at_faces

end module
