      !-------------------------------------------------------------------------
      !  Module       :              ppm_module_topo_copy
      !-------------------------------------------------------------------------
      !
      !  Purpose      : This module includes the source code for the routines
      !                 to copy topology structures
      !
      !  Remarks      :
      !
      !  References   :
      !
      !  Revisions    :
      !-------------------------------------------------------------------------
      !  $Log: ppm_module_topo_copy.f,v $
      !  Revision 1.1.1.1  2007/07/13 10:19:07  ivos
      !  CBL version of the PPM library
      !
      !-------------------------------------------------------------------------
      !  Parallel Particle Mesh Library (PPM)
      !  Institute of Computational Science
      !  ETH Zentrum, Hirschengraben 84
      !  CH-8092 Zurich, Switzerland
      !-------------------------------------------------------------------------
     
      !-------------------------------------------------------------------------
      !  Define types
      !-------------------------------------------------------------------------
#define __SINGLE_PRECISION 1
#define __DOUBLE_PRECISION 2

      MODULE ppm_module_topo_copy
         !----------------------------------------------------------------------
         !  Define interfaces to topology copy routine
         !----------------------------------------------------------------------
         INTERFACE ppm_topo_copy
            MODULE PROCEDURE ppm_topo_copy
         END INTERFACE

         !----------------------------------------------------------------------
         !  include the source 
         !----------------------------------------------------------------------
         CONTAINS

#include "ppm_topo_copy.f"

      END MODULE ppm_module_topo_copy
