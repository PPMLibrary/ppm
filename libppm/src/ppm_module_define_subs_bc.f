      !-------------------------------------------------------------------------
      !  Module       :              ppm_module_define_subs_bc
      !-------------------------------------------------------------------------
      !
      !  Purpose      : This module includes the source code for the
      !                 decomposition routines. 
      !
      !  Remarks      :
      !
      !  References   :
      !
      !  Revisions    :
      !-------------------------------------------------------------------------
      !  $Log: ppm_module_define_subs_bc.f,v $
      !  Revision 1.1.1.1  2007/07/13 10:18:57  ivos
      !  CBL version of the PPM library
      !
      !  Revision 1.1  2004/07/26 08:55:32  ivos
      !  Renamed.
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

      MODULE ppm_module_define_subs_bc

         !----------------------------------------------------------------------
         !  Define interface to boundary condition definition routine
         !----------------------------------------------------------------------
         INTERFACE ppm_define_subs_bc
            MODULE PROCEDURE ppm_define_subs_bc_s
            MODULE PROCEDURE ppm_define_subs_bc_d
         END INTERFACE

         !----------------------------------------------------------------------
         !  include the source 
         !----------------------------------------------------------------------
         CONTAINS

#define __KIND __SINGLE_PRECISION
#include "ppm_define_subs_bc.f"
#undef __KIND

#define __KIND __DOUBLE_PRECISION
#include "ppm_define_subs_bc.f"
#undef __KIND

      END MODULE ppm_module_define_subs_bc
