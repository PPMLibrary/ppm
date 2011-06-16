      !-------------------------------------------------------------------------
      !  Module       :            ppm_module_impose_part_bc
      !-------------------------------------------------------------------------
      !
      !  Purpose      : This module includes the source code user callable 
      !                 ppm_impose_part_bc routine to map the particles back
      !                 inside the computational domain for periodic BC.
      !                
      !  Remarks      : 
      !
      !  References   :
      !
      !  Revisions    :
      !-------------------------------------------------------------------------
      !  $Log: ppm_module_impose_part_bc.f,v $
      !  Revision 1.1.1.1  2007/07/13 10:18:58  ivos
      !  CBL version of the PPM library
      !
      !  Revision 1.1  2004/08/03 08:07:10  walther
      !  Initial version (untested).
      !
      !-------------------------------------------------------------------------
      !  Perallel Particle Mesh Library (PPM)
      !  Institute of Computational Science
      !  ETH Zentrum, Hirschengraben 84
      !  CH-8092 Zurich, Switzerland
      !-------------------------------------------------------------------------

      !-------------------------------------------------------------------------
      !  Define types
      !-------------------------------------------------------------------------
#define __SINGLE_PRECISION         1
#define __DOUBLE_PRECISION         2

      MODULE ppm_module_impose_part_bc

         !----------------------------------------------------------------------
         !  Define interface to ppm_map_field_ghost
         !----------------------------------------------------------------------
         INTERFACE ppm_impose_part_bc
             MODULE PROCEDURE ppm_impose_part_bc_s
             MODULE PROCEDURE ppm_impose_part_bc_d
         END INTERFACE

         !----------------------------------------------------------------------
         !  Include the source
         !----------------------------------------------------------------------
         CONTAINS

#define __KIND __SINGLE_PRECISION
#include "ppm_impose_part_bc.f"
#undef __KIND

#define __KIND __DOUBLE_PRECISION
#include "ppm_impose_part_bc.f"
#undef __KIND

      END MODULE ppm_module_impose_part_bc
