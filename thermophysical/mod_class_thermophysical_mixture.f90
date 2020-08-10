!>
!! \brief mod_class_thermophysical_mixture calculates thermophysical properties
!!        of a gas mixture by calculating the respective properties of indivi-
!!        dual specie. A database of individual species must be provided. For
!!        more information, see:
!!
!!        B.J. McBride, S. Gordon, M.A. Reno, "Coefficients for Calculating
!!        Thermodynamic and Transport Properties of Individual Species",
!!        Technical Memorandum 4513, NASA, 1993
!!
module mod_class_thermophysical_mixture

   use mod_class_thermophysical_abstract
   use mod_class_ifile

   implicit none

   ! Makes everything private, except otherwise stated
   private

   !> \brief Public class that calculates the thermophysical properties of a gas
   !!        mixture according to thermophysical properties of each individual
   !!        specie.
   type, public, extends(class_thermophysical_abstract) :: class_thermophysical_mixture

      integer :: ngs  !< Number of gaseous species
      integer :: base !< Fraction base ( 0 = molar, 1 = massic )

      real(8), allocatable :: Xi(:)  !< Molar fraction
      real(8), allocatable :: Yi(:)  !< Massic fraction
      real(8), allocatable :: Mi(:)  !< Molar mass of the gaseous specie i (kg/mol)
      real(8), allocatable :: Ri(:)  !< Gas constant of the gaseous specie i (J/kg.K)

      real(8), allocatable :: cpi(:,:) !< Coefficients for the calculation of cp of specie i
      real(8), allocatable :: mui(:,:) !< Coefficients for the calculation of mu of specie i
      real(8), allocatable :: kpi(:,:) !< Coefficients for the calculation of kp of specie i

      character(10), allocatable :: NGi(:)  !< Name of the gaseous specie i

   contains

      procedure, public, pass :: init  !< Initializer
      procedure, public, pass :: mu    !< Viscosity (Pa.s)
      procedure, public, pass :: kp    !< Thermal conductivity (W/m.K)
      procedure, public, pass :: cp    !< Specific heat at constant pressure (J/(kg.K))
      procedure, public, pass :: gm    !< Specific heat ratio
      procedure, public, pass :: about !< Prints information about the class

   end type


contains


   !> \brief Initializes variables related to the class_thermophysical_mixture
   subroutine init(this, ifile)
      implicit none
      class(class_thermophysical_mixture) :: this  !< A reference to this object
      class(class_ifile),      intent(in) :: ifile !< Input file


      ! Inner variables
      integer       ::    i ! Dummy index
      real(8)       :: raux ! Auxiliary variable
      character(10) :: caux ! Auxiliary variable

      character(len=500) :: file_parameters ! Path and file name of the input parameters
      character(len=500) :: file_database   ! Path and file name of the database


      ! Considering ifile is loaded...
      call ifile%get_value(file_parameters, "thermo_mix_parameters")
      call ifile%get_value(  file_database, "thermo_mix_database"  )


      ! Universal gas constant (J/K.mol)
      this%Ru = 8.31451000D+00

      ! Reading the number of gaseous species in the mixture

      this%ngs = 0

      open(10, file = file_parameters )

      read(10,*) this%base

      do

         read(10,*) raux, caux

         if ( trim(adjustl(caux)) == "END" ) exit

         if ( raux > 1.d-8 ) this%ngs = this%ngs + 1

      end do


      ! Allocating memory

      allocate(  this%Xi(this%ngs)    )
      allocate(  this%Yi(this%ngs)    )
      allocate(  this%Mi(this%ngs)    )
      allocate(  this%Ri(this%ngs)    )
      allocate( this%NGi(this%ngs)    )
      allocate( this%cpi(this%ngs,10) )
      allocate( this%mui(this%ngs, 8) )
      allocate( this%kpi(this%ngs, 8) )


      ! Reading the fractions for gasous mixture

      rewind(10)

      read(10,*) this%base

      if ( this%base == 0 ) then ! Molar fraction

         i = 0

         do

            read(10,*) raux, caux

            if ( trim(adjustl(caux)) == "END" ) exit

            if ( raux > 1.d-8 ) then

               i = i + 1

               this%Xi(i) = raux

               this%NGi(i) = caux

            end if

         end do

         ! Checking consistency of input data

         if ( abs( 1.d0 - sum(this%Xi) ) > epsilon(1.d0) * 10.d0 ) then

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

               this%Yi(i) = raux

               this%NGi(i) = caux

            end if

         end do

         ! Checking consistency of input data

         if ( abs( 1.d0 - sum(this%Yi) ) > epsilon(1.d0) * 10.d0 ) then

            write(*,*) "ERROR: The sum of the massic fraction is not equal 1."

            stop

         end if

      end if

      close(10)


      ! Reading database

      do i = 1, this%ngs

         call get_specie_data( file_database &
            ,                  this%NGi(i)   &
            ,                  this%Mi(i)    &
            ,                  this%cpi(i,:) &
            ,                  this%mui(i,:) &
            ,                  this%kpi(i,:) )

      end do


      ! Calculating extra data

      if ( this%base == 0 ) then

         this%Mg = dot_product(this%Xi,this%Mi)

         this%Yi = this%Xi * this%Mi / this%Mg

      else

         this%Xi = ( this%Yi / this%Mi ) / sum( this%Yi / this%Mi )

         this%Mg = dot_product(this%Xi,this%Mi)

      end if


      ! Gas constant of each specie
      this%Ri = this%Ru / this%Mi


      ! Gas constant of the mixture
      this%Rg = this%Ru / this%Mg

   end subroutine


   !> \brief Calculates the viscosity (mu) as a function of T
   real(8) function mu(this, T)
      implicit none
      class(class_thermophysical_mixture) :: this !< A reference to this object
      real(8),                 intent(in) :: T    !< Temperature (K)

      mu = mu_mix(this%ngs, T, this%Xi, this%Mi, this%mui)

   end function


   !> \brief Calculates the thermal conductivity (kp) as a function of T
   real(8) function kp(this, T)
      implicit none
      class(class_thermophysical_mixture) :: this !< A reference to this object
      real(8),                 intent(in) :: T    !< Temperature (K)

      kp = kp_mix(this%ngs, T, this%Xi, this%Mi, this%kpi)

   end function


   !> \brief Calculates the specific heat at const. p (cp) as a function of T
   real(8) function cp(this, T)
      implicit none
      class(class_thermophysical_mixture) :: this !< A reference to this object
      real(8),                 intent(in) :: T    !< Temperature (K)

      cp =  cp_mix(this%ngs, T, this%Yi, this%Ri, this%cpi)

   end function


   !> \brief Calculates the specific heat ratio as a function of T
   real(8) function gm(this, T)
      implicit none
      class(class_thermophysical_mixture) :: this !< A reference to this object
      real(8),                 intent(in) :: T    !< Temperature (K)

      real(8) :: c

      c = this%cp(T)

      gm = c / (c - this%Rg)

   end function


   !> \brief Prints information about the class
   subroutine about(this, unt)
      implicit none
      class(class_thermophysical_mixture) :: this !< A reference to this object
      integer,                 intent(in) :: unt  !< Unit number where to print

      ! Inner variables
      integer :: i

      write(unt,*)
      write(unt,*) "class_thermophysical_mixture:"
      write(unt,*) "Thermophysical properties are temperature dependent."
      write(unt,*) "Thermophysical properties of the gas mixture calculated in"
      write(unt,*) "accordance to the thermophysical properties of each specie."
      write(unt,*)
      write(unt,*) "Gas composition"
      write(unt,*)

      write(unt,"(A10,5(2X,A14))") "Specie", "Xi", "Yi", "Mi (kg/mol)", "Ri (J/kg.K)"

      do i = 1, this%ngs

         write(unt,"(A10,5(2X,ES14.7))") trim(this%NGi(i)), this%Xi(i), this%Yi(i), this%Mi(i), this%Ri(i)

      end do

      write(unt,*)
      write(unt,*) "Gas coefficients for Cp, mu and kappa"
      write(unt,*)

      do i = 1, this%ngs

         write(unt,*) trim(this%NGi(i))
         write(unt,"(A20,10(2X,ES14.7))") "Cp ( 200K<=T<=1000K): ", this%cpi(i,1:5)
         write(unt,"(A20,10(2X,ES14.7))") "Cp (1000K<=T<=6000K): ", this%cpi(i,6:10)
         write(unt,"(A20,10(2X,ES14.7))") "mu ( 300K<=T<=1000K): ", this%mui(i,1:4)
         write(unt,"(A20,10(2X,ES14.7))") "mu (1000K<=T<=5000K): ", this%mui(i,5:8)
         write(unt,"(A20,10(2X,ES14.7))") "kp ( 300K<=T<=1000K): ", this%kpi(i,1:4)
         write(unt,"(A20,10(2X,ES14.7))") "kp (1000K<=T<=5000K): ", this%kpi(i,5:8)
         write(unt,*)

      end do

      write(unt,*)
      write(unt,*) "Gas properties"
      write(unt,*)
      write(unt,"(ES23.16,2X,A)") this%Ru, " = Ru: universal gas constant (J/K.mol)"
      write(unt,"(ES23.16,2X,A)") this%Rg, " = Rg: gas constant (J/K.kg)"
      write(unt,"(ES23.16,2X,A)") this%Mg, " = Mg: molar mass (kg/mol)"
      write(unt,*)

   end subroutine


   !> \brief Extracts data from database for a given gas specie
   subroutine get_specie_data(filedb, NGi, Mi, cpi, mui, kpi)
      implicit none
      character(len=*),       intent(in)  ::  filedb !< File of individual specie database
      character(10),          intent(in)  ::     NGi !< Name of the gaseous specie i
      real(8),                intent(out) ::      Mi !< Molar mass of specie i (kg/mol)
      real(8),                intent(out) :: cpi(10) !< Coefficients for the calculation of cp of specie i
      real(8),                intent(out) :: mui( 8) !< Coefficients for the calculation of mu of specie i
      real(8),                intent(out) :: kpi( 8) !< Coefficients for the calculation of kp of specie i

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

      open(10, file = filedb )


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

end module
