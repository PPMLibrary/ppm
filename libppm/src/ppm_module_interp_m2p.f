      !-------------------------------------------------------------------------
      !  Module       :            ppm_module_rmsh_remesh
      !-------------------------------------------------------------------------
      !
      !  Purpose      : procedure module for ppm_rmsh_remesh
      !
      !  Remarks      : 
      !
      !  References   : 
      !
      !  Revisions    :
      !-------------------------------------------------------------------------
      !  $Log: ppm_module_interp_m2p.f,v $
      !  Revision 1.1.1.1  2007/07/13 10:18:58  ivos
      !  CBL version of the PPM library
      !
      !  Revision 1.3  2006/10/20 14:07:45  ivos
      !  2d_vec and 3d_sca cases were commented out and there was a bug in the
      !  argument declaration in 2d_vec. Fixed and put them back in.
      !
      !  Revision 1.2  2005/05/29 20:14:56  michaebe
      !  some fixes for the 2D case
      !
      !  Revision 1.1  2004/11/02 12:52:21  michaebe
      !  inimp
      !
      !-------------------------------------------------------------------------
      !  Parallel Particle Mesh Library (PPM)
      !  Institute of Computational Science
      !  ETH Zentrum, Hirschengraben 84
      !  CH-8092 Zurich, Switzerland
      !-------------------------------------------------------------------------

#define __SINGLE_PRECISION 1
#define __DOUBLE_PRECISION 2
#define __2D               3
#define __3D               4
#define __VEC              5
#define __SCA              6

      MODULE ppm_module_interp_m2p
      
        !-----------------------------------------------------------------------
        !  Interface
        !-----------------------------------------------------------------------
        INTERFACE ppm_interp_m2p
           ! 2d scalar
           MODULE PROCEDURE ppm_interp_m2p_ss_2d
           MODULE PROCEDURE ppm_interp_m2p_ds_2d
           ! 2d vector
           MODULE PROCEDURE ppm_interp_m2p_sv_2d
           MODULE PROCEDURE ppm_interp_m2p_dv_2d
           ! 3d scalar
           MODULE PROCEDURE ppm_interp_m2p_ss_3d
           MODULE PROCEDURE ppm_interp_m2p_ds_3d
           ! 3d vector
           MODULE PROCEDURE ppm_interp_m2p_sv_3d
           MODULE PROCEDURE ppm_interp_m2p_dv_3d
        END INTERFACE
        
      CONTAINS

        
#define __KIND  __SINGLE_PRECISION
#define __DIME  __2D
#define __MODE  __SCA
        ! 2D SCA SINGLE
#include "ppm_interp_m2p.f"
#undef  __MODE
#define __MODE  __VEC
        ! 2D VEC SINGLE
#include "ppm_interp_m2p.f"
#undef  __MODE
#undef  __DIME
        
#define __DIME  __3D
#define __MODE  __SCA
        ! 3D SCA SINGLE
#include "ppm_interp_m2p.f"
#undef  __MODE
#define __MODE  __VEC
        ! 3D VEC SINGLE
#include "ppm_interp_m2p.f"
#undef  __MODE
#undef  __DIME
#undef  __KIND


#define __KIND  __DOUBLE_PRECISION
#define __DIME  __2D
#define __MODE  __SCA
        ! 2D SCA DOUBLE
#include "ppm_interp_m2p.f"
#undef  __MODE
#define __MODE  __VEC
        ! 2D VEC DOUBLE
#include "ppm_interp_m2p.f"
#undef  __MODE
#undef  __DIME
        
#define __DIME  __3D
#define __MODE  __SCA
        ! 3D SCA DOUBLE
#include "ppm_interp_m2p.f"
#undef  __MODE
#define __MODE  __VEC
        ! 3D VEC DOUBLE
#include "ppm_interp_m2p.f"
#undef  __MODE
#undef  __DIME
#undef  __KIND        





#undef __SINGLE_PRECISION 
#undef __DOUBLE_PRECISION 
#undef __2D               
#undef __3D               
#undef __VEC              
#undef __SCA              
        
      END MODULE ppm_module_interp_m2p

