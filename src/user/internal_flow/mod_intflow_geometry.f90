!>
!! \brief This module defines data structures and procedures related to the
!!        geometry of internal flows.
!!
module mod_intflow_geometry

   implicit none

   !> \brief Geometric data to be used in the internal flow calculation
   type, public :: type_intflow_geometry

      integer :: ig  !< i=ig when the throttle cross section is in the east face of the CV
      real(8) :: Sg  !< Throat area (m2)
      real(8) :: Rcg !< Throat curvature radius (m)

   end type

contains

end module
