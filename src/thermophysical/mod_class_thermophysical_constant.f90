!>
!! \brief mod_class_thermophysical_constant returns constant thermophysical
!!        properties of a gas
!!
module mod_class_thermophysical_constant

   use mod_class_thermophysical_abstract
   use mod_class_ifile

   implicit none

   ! Makes everything private, except otherwise stated
   private

   !> \brief Public class that returns the thermophysical properties of a gas
   type, public, extends(class_thermophysical_abstract) :: class_thermophysical_constant

      real(8) :: gm_val
      real(8) :: mu_val
      real(8) :: kp_val
      real(8) :: cp_val

   contains

      procedure, public, pass :: init  !< Initializer
      procedure, public, pass :: mu    !< Viscosity (Pa.s)
      procedure, public, pass :: kp    !< Thermal conductivity (W/m.K)
      procedure, public, pass :: cp    !< Specific heat at constant pressure (J/(kg.K))
      procedure, public, pass :: gm    !< Specific heat ratio
      procedure, public, pass :: about !< Prints information about the class

   end type


contains


   !> \brief Initializes variables related to the class_thermophysical_constant
   subroutine init(this, ifile)
      implicit none
      class(class_thermophysical_constant) :: this  !< A reference to this object
      class(class_ifile),       intent(in) :: ifile !< Input file

      ! Assuming that the file was loaded, just copy the necessary data
      call ifile%get_value(     this%Rg, "thermo_const_Rg"    )
      call ifile%get_value(     this%Mg, "thermo_const_Mg"    )
      call ifile%get_value( this%mu_val, "thermo_const_mu"    )
      call ifile%get_value( this%kp_val, "thermo_const_kappa" )
      call ifile%get_value( this%gm_val, "thermo_const_gamma" )

      this%cp_val = this%gm_val * this%Rg / (this%gm_val-1.d0)

   end subroutine


   !> \brief Calculates the viscosity (mu) as a function of T
   real(8) function mu(this, T)
      implicit none
      class(class_thermophysical_constant) :: this !< A reference to this object
      real(8),                  intent(in) :: T    !< Temperature (K)

      mu = this%mu_val

   end function


   !> \brief Calculates the thermal conductivity (kp) as a function of T
   real(8) function kp(this, T)
      implicit none
      class(class_thermophysical_constant) :: this !< A reference to this object
      real(8),                  intent(in) :: T    !< Temperature (K)

      kp = this%kp_val

   end function


   !> \brief Calculates the specific heat at const. p (cp) as a function of T
   real(8) function cp(this, T)
      implicit none
      class(class_thermophysical_constant) :: this !< A reference to this object
      real(8),                  intent(in) :: T    !< Temperature (K)

      cp =  this%cp_val

   end function


   !> \brief Calculates the specific heat ratio as a function of T
   real(8) function gm(this, T)
      implicit none
      class(class_thermophysical_constant) :: this !< A reference to this object
      real(8),                  intent(in) :: T    !< Temperature (K)

      gm =  this%gm_val

   end function


   !> \brief Prints information about the class
   subroutine about(this, unt)
      implicit none
      class(class_thermophysical_constant) :: this !< A reference to this object
      integer,                  intent(in) :: unt  !< Unit number where to print


      write(unt,*)
      write(unt,*) "class_thermophysical_constant:"
      write(unt,*) "Thermophysical properties are constant and user defined."
      write(unt,*)
      write(unt,*) "Gas properties"
      write(unt,*)
      write(unt,"(ES23.16,2X,A)") this%Ru,     " = Ru: universal gas constant (J/K.mol)"
      write(unt,"(ES23.16,2X,A)") this%Rg,     " = Rg: gas constant (J/K.kg)"
      write(unt,"(ES23.16,2X,A)") this%Mg,     " = Mg: molar mass (kg/mol)"
      write(unt,"(ES23.16,2X,A)") this%gm_val, " = gm: specific heat ratio"
      write(unt,"(ES23.16,2X,A)") this%cp_val, " = cp: specific heat at constant pressure (J/(kg.K))"
      write(unt,"(ES23.16,2X,A)") this%mu_val, " = mu: viscosity (Pa.s)"
      write(unt,"(ES23.16,2X,A)") this%kp_val, " = kp: thermal conductivity (W/m.K)"
      write(unt,*)

   end subroutine

end module
