!>
!! \brief Defines an abstract class for creating thermophysical models
!!
module mod_class_thermophysical_abstract

   implicit none

   ! Makes everything private, except otherwise stated
   private

   !> \brief Public class defining an interface for thermophysical models
   type, abstract, public :: class_thermophysical_abstract

      real(8) :: Ru = 8.31451000D+00 !< Universal gas constant (J/mol.K)
      real(8) :: Rg                  !< Gas constant of the mixture (J/kg.K)
      real(8) :: Mg                  !< Molar mass of the mixture (kg/mol)
   contains

      procedure(mu_interface), pass, deferred, public :: mu    !< Viscosity (Pa.s)
      procedure(kp_interface), pass, deferred, public :: kp    !< Thermal conductivity (W/m.K)
      procedure(cp_interface), pass, deferred, public :: cp    !< Spec. heat at const. pressure (J/kg.K)
      procedure(gm_interface), pass, deferred, public :: gm    !< Spec. heat ratio
      procedure(ab_interface), pass, deferred, public :: about !< Prints information about the class

   end type

   abstract interface

      !> \brief Computes viscosity as a function of temperature
      real(8) function mu_interface(this, T)
         import class_thermophysical_abstract
         implicit none
         class(class_thermophysical_abstract)  :: this !< A reference to this object
         real(8),                   intent(in) :: T    !< Temperature (K)
      end function


      !> \brief Computes thermal conductivity as a function of temperature
      real(8) function kp_interface(this, T)
         import class_thermophysical_abstract
         implicit none
         class(class_thermophysical_abstract)  :: this !< A reference to this object
         real(8),                   intent(in) :: T    !< Temperature (K)
      end function


      !> \brief Computes spec. heat at const. T as a function of temperature
      real(8) function cp_interface(this, T)
         import class_thermophysical_abstract
         implicit none
         class(class_thermophysical_abstract)  :: this !< A reference to this object
         real(8),                   intent(in) :: T    !< Temperature (K)
      end function


      !> \brief Computes spec. heat ratio as a function of temperature
      real(8) function gm_interface(this, T)
         import class_thermophysical_abstract
         implicit none
         class(class_thermophysical_abstract)  :: this !< A reference to this object
         real(8),                   intent(in) :: T    !< Temperature (K)
      end function


      !> \brief Prints information about the class
      subroutine ab_interface(this, unt)
         import class_thermophysical_abstract
         implicit none
         class(class_thermophysical_abstract)  :: this !< A reference to this object
         integer,                   intent(in) :: unt  !< Unit number where to print

      end subroutine

   end interface

end module
