!>
!! \brief "mod_extflow_data" contains data for external flow calculations
!!
module mod_extflow_data
   implicit none

   !
   ! GEOMETRIC PARAMETERS
   !
   integer :: iocs !< index of the matching point between the ogive and the cylinder
   integer :: kfc  !< Kind of foredrag calculation ( 0 = over the whole forebody; 1 = over the ogive only)
   real(8) :: lr   !< length of the body
   real(8) :: rb   !< base radius of the body
   character(len=200) :: fgeom !< File of the geometric parameters
   !
   ! FREESTREAM PARAMETERS
   !
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
   !
   ! VARIABLES OF INTEREST
   !
   real(8) :: Cdfi !< Foredrag coefficient due pressure (dimensionless)
   real(8) :: Cdfv !< Foredrag coefficient due viscosity (dimensionless)
   !
   ! BOUNDARY CONDITION PARAMETERS
   !
   real(8) :: Tsbc !< Temperature on the south boundary (K) (if negative, adiabatic bc is applied)

end module
