      !-------------------------------------------------------------------------
      !  Module       :            ppm_module_rmsh_comp_weights
      !-------------------------------------------------------------------------
      !
      !  Purpose      : procedure module for ppm_rmsh_comp_weights
      !
      !  Remarks      : 
      !
      !  References   : 
      !
      !  Revisions    :
      !-------------------------------------------------------------------------
      !  $Log: ppm_module_rmsh_comp_weights.f,v $
      !  Revision 1.1.1.1  2007/07/13 10:19:00  ivos
      !  CBL version of the PPM library
      !
      !  Revision 1.1  2004/08/09 11:08:32  michaebe
      !  revised initial implementation (P. Gautier)
      !
      !-------------------------------------------------------------------------
      !  Parallel Particle Mesh Library (PPM)
      !  Institute of Computational Science
      !  ETH Zentrum, Hirschengraben 84
      !  CH-8092 Zurich, Switzerland
      !-------------------------------------------------------------------------

#define __SINGLE_PRECISION 1
#define __DOUBLE_PRECISION 2

      MODULE ppm_module_rmsh_comp_weights
      
        !-----------------------------------------------------------------------
        !  Interface
        !-----------------------------------------------------------------------
        INTERFACE ppm_rmsh_comp_weights
           MODULE PROCEDURE ppm_rmsh_comp_weights_s
           MODULE PROCEDURE ppm_rmsh_comp_weights_d
        END INTERFACE
        
      CONTAINS

#define __KIND  __SINGLE_PRECISION
#include "ppm_rmsh_comp_weights.f"
#undef  __KIND        
#define __KIND __DOUBLE_PRECISION 
#include "ppm_rmsh_comp_weights.f"
#undef  __KIND        

      END MODULE ppm_module_rmsh_comp_weights


        
