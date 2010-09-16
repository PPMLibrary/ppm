      !--*- f90 -*--------------------------------------------------------------
      !  Module       :                ppm_module_data_state
      !-------------------------------------------------------------------------
      !  Parallel Particle Mesh Library (PPM)
      !  ETH Zurich
      !  CH-8092 Zurich, Switzerland
      !-------------------------------------------------------------------------

      MODULE ppm_module_data_state
      !!! This module contains  data structures and definitions that
      !!! are PRIVATE to the ppm_map_part_store and ppm_map_part_load routines.
      !!!
      !!! [NOTE]
      !!! The variables declared in this modules should not be directly accessed
      !!! by the user. They are managed interally by the library.
         !----------------------------------------------------------------------
         !  
         !----------------------------------------------------------------------
         INTEGER :: ppm_map_type_state 
         INTEGER :: ppm_nrecvlist_state
         INTEGER :: ppm_nsendlist_state 
         INTEGER :: ppm_nsendbuffer_state 
         INTEGER :: ppm_buffer_set_state 

         INTEGER, DIMENSION(:), POINTER :: ppm_psendbuffer_state
         INTEGER, DIMENSION(:), POINTER :: ppm_buffer2part_state

         INTEGER, DIMENSION(:), POINTER :: ppm_irecvlist_state
         INTEGER, DIMENSION(:), POINTER :: ppm_isendlist_state

      END MODULE ppm_module_data_state
