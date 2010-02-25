      !-------------------------------------------------------------------------
      !  Module       :            ppm_module_rmsh_create_part
      !-------------------------------------------------------------------------
      !
      !  Purpose      : procedure module for ppm_rmsh_create_part
      !
      !  Remarks      : 
      !
      !  References   : 
      !
      !  Revisions    :
      !-------------------------------------------------------------------------
      !  $Log: ppm_module_rmsh_create_part.f,v $
      !  Revision 1.1.1.1  2007/07/13 10:19:00  ivos
      !  CBL version of the PPM library
      !
      !  Revision 1.5  2005/06/21 17:36:09  ivos
      !  bugfix in overloading.
      !
      !  Revision 1.4  2005/06/21 01:07:57  ivos
      !  Added the OPTIONAL slave arrays to create multiple particle
      !  properties (maybe by calling this several times). Overloaded for
      !  scalar and vector versions.
      !
      !  Revision 1.3  2005/01/21 11:15:59  michaebe
      !  removed jumps again
      !
      !  Revision 1.2  2005/01/13 13:21:59  ivos
      !  Commit with Michaels new vectorized interp and rmsh.
      !
      !  Revision 1.1  2004/08/11 08:59:24  michaebe
      !  initial implementaiton
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

      MODULE ppm_module_rmsh_create_part
      
        !-----------------------------------------------------------------------
        !  Interface
        !-----------------------------------------------------------------------
        INTERFACE ppm_rmsh_create_part
           ! 2d scalar
           MODULE PROCEDURE ppm_rmsh_create_part_sss_2d
           MODULE PROCEDURE ppm_rmsh_create_part_ssv_2d
           MODULE PROCEDURE ppm_rmsh_create_part_dss_2d
           MODULE PROCEDURE ppm_rmsh_create_part_dsv_2d
           ! 2d vector
           MODULE PROCEDURE ppm_rmsh_create_part_svs_2d
           MODULE PROCEDURE ppm_rmsh_create_part_svv_2d
           MODULE PROCEDURE ppm_rmsh_create_part_dvs_2d
           MODULE PROCEDURE ppm_rmsh_create_part_dvv_2d
           ! 3d scalar
           MODULE PROCEDURE ppm_rmsh_create_part_sss_3d
           MODULE PROCEDURE ppm_rmsh_create_part_ssv_3d
           MODULE PROCEDURE ppm_rmsh_create_part_dss_3d
           MODULE PROCEDURE ppm_rmsh_create_part_dsv_3d
           ! 3d vector
           MODULE PROCEDURE ppm_rmsh_create_part_svs_3d
           MODULE PROCEDURE ppm_rmsh_create_part_svv_3d
           MODULE PROCEDURE ppm_rmsh_create_part_dvs_3d
           MODULE PROCEDURE ppm_rmsh_create_part_dvv_3d
        END INTERFACE

      CONTAINS

        
#define __KIND  __SINGLE_PRECISION
#define __DIME  __2D
#define __MODE  __SCA
        ! 2D SCA SINGLE
#define __MODE2  __SCA
#include "ppm_rmsh_create_part.f"
#undef __MODE2
#define __MODE2  __VEC
#include "ppm_rmsh_create_part.f"
#undef __MODE2

#undef  __MODE
#define __MODE  __VEC
        ! 2D VEC SINGLE
#define __MODE2  __SCA
#include "ppm_rmsh_create_part.f"
#undef __MODE2
#define __MODE2  __VEC
#include "ppm_rmsh_create_part.f"
#undef __MODE2
#undef  __MODE
#undef  __DIME
        
#define __DIME  __3D
#define __MODE  __SCA
        ! 3D SCA SINGLE
#define __MODE2 __SCA
#include "ppm_rmsh_create_part.f"
#undef __MODE2
#define __MODE2 __VEC
#include "ppm_rmsh_create_part.f"
#undef __MODE2
#undef  __MODE
#define __MODE  __VEC
        ! 3D VEC SINGLE
#define __MODE2 __SCA
#include "ppm_rmsh_create_part.f"
#undef __MODE2
#define __MODE2 __VEC
#include "ppm_rmsh_create_part.f"
#undef __MODE2
#undef  __MODE
#undef  __DIME
#undef  __KIND


#define __KIND  __DOUBLE_PRECISION
#define __DIME  __2D
#define __MODE  __SCA
        ! 2D SCA DOUBLE
#define __MODE2 __SCA
#include "ppm_rmsh_create_part.f"
#undef __MODE2
#define __MODE2 __VEC
#include "ppm_rmsh_create_part.f"
#undef __MODE2
#undef  __MODE
#define __MODE  __VEC
        ! 2D VEC DOUBLE
#define __MODE2 __SCA
#include "ppm_rmsh_create_part.f"
#undef __MODE2
#define __MODE2 __VEC
#include "ppm_rmsh_create_part.f"
#undef __MODE2
#undef  __MODE
#undef  __DIME
        
#define __DIME  __3D
#define __MODE  __SCA
        ! 3D SCA DOUBLE
#define __MODE2 __SCA
#include "ppm_rmsh_create_part.f"
#undef __MODE2
#define __MODE2 __VEC
#include "ppm_rmsh_create_part.f"
#undef __MODE2
#undef  __MODE
#define __MODE  __VEC
        ! 3D VEC DOUBLE
#define __MODE2 __SCA
#include "ppm_rmsh_create_part.f"
#undef __MODE2
#define __MODE2 __VEC
#include "ppm_rmsh_create_part.f"
#undef __MODE2
#undef  __MODE
#undef  __DIME
#undef  __KIND        



#undef __SINGLE_PRECISION 
#undef __DOUBLE_PRECISION 
#undef __2D               
#undef __3D               
#undef __VEC              
#undef __SCA              
        
      END MODULE ppm_module_rmsh_create_part

