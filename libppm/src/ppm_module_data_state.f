      !-------------------------------------------------------------------------
      !  Module       :                ppm_module_data_state
      !-------------------------------------------------------------------------
      !
      !  Purpose      : This module contains all data structures and
      !                 definitions that are PRIVATE to the ppm_map_part_store
      !                 and ppm_map_part_load routines.
      !                
      !  Remarks      : 
      !
      !  References   :
      !
      !  Revisions    :
      !-------------------------------------------------------------------------
      !  $Log: ppm_module_data_state.f,v $
      !  Revision 1.1.1.1  2007/07/13 10:18:57  ivos
      !  CBL version of the PPM library
      !
      !  Revision 1.1  2006/10/10 21:25:38  walther
      !  *** empty log message ***
      !
      !-------------------------------------------------------------------------
      !  Perallel Particle Mesh Library (PPM)
      !  Institute of Computational Science
      !  ETH Zentrum, Hirschengraben 84
      !  CH-8092 Zurich, Switzerland
      !-------------------------------------------------------------------------

      MODULE ppm_module_data_state

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
