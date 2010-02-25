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
      !  $Log: ppm_module_interp_p2m.f,v $
      !  Revision 1.1.1.1  2007/07/13 10:18:58  ivos
      !  CBL version of the PPM library
      !
      !  Revision 1.7  2005/01/21 11:16:29  michaebe
      !  jumps gone now
      !
      !  Revision 1.6  2005/01/13 13:21:59  ivos
      !  Commit with Michaels new vectorized interp and rmsh.
      !
      !  Revision 1.5  2004/12/06 13:53:54  ivos
      !  Re-activated 2d stuff.
      !
      !  Revision 1.4  2004/12/06 12:39:16  ivos
      !  Commented 2d cases. They do not (yet) compile.
      !
      !  Revision 1.3  2004/12/06 12:35:37  michaebe
      !  added dummy 2d routines
      !
      !  Revision 1.2  2004/11/02 16:33:53  michaebe
      !  added scalar case
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

      MODULE ppm_module_interp_p2m
      
        !-----------------------------------------------------------------------
        !  Interface
        !-----------------------------------------------------------------------
        INTERFACE ppm_interp_p2m
           ! 2d scalar
           MODULE PROCEDURE ppm_interp_p2m_ss_2d
           MODULE PROCEDURE ppm_interp_p2m_ds_2d
           ! 2d vector
           MODULE PROCEDURE ppm_interp_p2m_sv_2d
           MODULE PROCEDURE ppm_interp_p2m_dv_2d
           ! 3d scalar
           MODULE PROCEDURE ppm_interp_p2m_ss_3d
           MODULE PROCEDURE ppm_interp_p2m_ds_3d
           ! 3d vector
           MODULE PROCEDURE ppm_interp_p2m_sv_3d
           MODULE PROCEDURE ppm_interp_p2m_dv_3d
        END INTERFACE




        
      CONTAINS

        
#define __KIND  __SINGLE_PRECISION
#define __DIME  __2D
#define __MODE  __SCA
        ! 2D SCA SINGLE
#include "ppm_interp_p2m.f"
#undef  __MODE
#define __MODE  __VEC
        ! 2D VEC SINGLE
#include "ppm_interp_p2m.f"
#undef  __MODE
#undef  __DIME
        
#define __DIME  __3D
#define __MODE  __SCA
        ! 3D SCA SINGLE
#include "ppm_interp_p2m.f"
#undef  __MODE
#define __MODE  __VEC
        ! 3D VEC SINGLE
#include "ppm_interp_p2m.f"

#undef  __MODE
#undef  __DIME
#undef  __KIND


#define __KIND  __DOUBLE_PRECISION
#define __DIME  __2D
#define __MODE  __SCA
        ! 2D SCA DOUBLE
#include "ppm_interp_p2m.f"
#undef  __MODE
#define __MODE  __VEC
        ! 2D VEC DOUBLE
#include "ppm_interp_p2m.f"
#undef  __MODE
#undef  __DIME
        
#define __DIME  __3D
#define __MODE  __SCA
        ! 3D SCA DOUBLE
#include "ppm_interp_p2m.f"
#undef  __MODE
#define __MODE  __VEC
        ! 3D VEC DOUBLE
#include "ppm_interp_p2m.f"
#undef  __MODE
#undef  __DIME
#undef  __KIND        





#undef __SINGLE_PRECISION 
#undef __DOUBLE_PRECISION 
#undef __2D               
#undef __3D               
#undef __VEC              
#undef __SCA              
        
      END MODULE ppm_module_interp_p2m

