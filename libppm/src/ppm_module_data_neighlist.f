      !--*- f90 -*--------------------------------------------------------------
      !  Module       :              ppm_module_data_neighlist
      !-------------------------------------------------------------------------
      !  Parallel Particle Mesh Library (PPM)
      !  ETH Zurich
      !  CH-8092 Zurich, Switzerland
      !-------------------------------------------------------------------------
     
      MODULE ppm_module_data_neighlist
      !!! This module provides data used by the neighbor search routines.
      !!!
      !!! [NOTE]
      !!! The variables declared in this module should not be accessed by the
      !!! PPM client developer. They are managed internally by the library.
         !----------------------------------------------------------------------
         !  Define data TYPEs
         !----------------------------------------------------------------------
         ! Pointer to cell list (needed to make lists of cell lists)
         TYPE ppm_type_ptr_to_clist
             ! particle index list
             INTEGER, DIMENSION(:), POINTER    :: lpdx
             ! first particle in each cell
             INTEGER, DIMENSION(:), POINTER    :: lhbx
         END TYPE

      END MODULE ppm_module_data_neighlist
