      !-------------------------------------------------------------------------
      !  Subroutine   :                ppm_util_cart2sph  
      !-------------------------------------------------------------------------
      !
      !  Purpose      : Routine for conversion between Charthesian and
      !                 Spherical co-ordinates. Notice, theta is angle of the 
      !                 vector with the z axis, and phi the angle in the x-y
      !                 plane.
      !
      !  Input        : x(:)            (F) : x co-ordinate
      !                 y(:)            (F) : y co-ordinate
      !                 z(:)            (F) : z co-ordinate
      !                 n:              (I) : number of points
      !                 x0              (F) : centre of co-ordinate system
      !                 y0              (F) : centre of co-ordinate system
      !                 z0              (F) : centre of co-ordinate system
      !
      !  Input/output :
      !
      !  Output       : r(:)            (F) : radius
      !                 theta(:)        (F) : angle
      !                 phi(:)          (F) : angle
      !                 info            (I) : status of routine
      !
      !
      !  Remarks      : The routine has by no means been optimised.
      !                  Adapted to ppm by Bettina Polasek.
      !
      !  Author       : Jens Honore Walther
      !
      !  Revisions    :
      !-------------------------------------------------------------------------
      !  $Log: ppm_util_cart2sph.f,v $
      !  Revision 1.1.1.1  2007/07/13 10:19:01  ivos
      !  CBL version of the PPM library
      !
      !  Revision 1.2  2006/09/04 18:34:57  pchatela
      !  Fixes and cleanups to make it compile for
      !  release-1-0
      !
      !  Revision 1.1  2005/02/09 15:26:26  polasekb
      !
      !  initial implementation
      !
      !-------------------------------------------------------------------------
      !  Parallel Particle Mesh Library (PPM)
      !  Institute of Computational Science
      !  ETH Zentrum, Hirschengraben 84
      !  CH-8092 Zurich, Switzerland
      !-------------------------------------------------------------------------

#if   __KIND == __SINGLE_PRECISION
      SUBROUTINE ppm_util_cart2sph_s(x,y,z,n,x0,y0,z0,r,theta,phi,info)
#elif __KIND == __DOUBLE_PRECISION
      SUBROUTINE ppm_util_cart2sph_d(x,y,z,n,x0,y0,z0,r,theta,phi,info)
#endif

      !-------------------------------------------------------------------------
      !  Modules
      !-------------------------------------------------------------------------
      USE ppm_module_data
      USE ppm_module_substart
      USE ppm_module_substop
      USE ppm_module_error
      USE ppm_module_write
      IMPLICIT NONE

#if   __KIND == __SINGLE_PRECISION
      INTEGER, PARAMETER :: MK = ppm_kind_single
#elif __KIND == __DOUBLE_PRECISION
      INTEGER, PARAMETER :: MK = ppm_kind_double
#endif

      !-------------------------------------------------------------------------
      !  Includes     
      !-------------------------------------------------------------------------
#include "ppm_define.h"
      !-------------------------------------------------------------------------
      !  Arguments     
      !-------------------------------------------------------------------------
      REAL(MK), DIMENSION(:), INTENT(IN)  :: x,y,z
      REAL(MK),               INTENT(IN)  :: x0,y0,z0
      REAL(MK), DIMENSION(:), INTENT(OUT) :: r,theta,phi
      INTEGER, INTENT(IN)                :: n
      INTEGER, INTENT(OUT)               :: info
      !-------------------------------------------------------------------------
      !  Local variables 
      !-------------------------------------------------------------------------
      INTEGER :: i,j,k
      REAL(MK) :: dx,dy,dz,rad,t0
      
      !-------------------------------------------------------------------------
      !  Initialise 
      !-------------------------------------------------------------------------
      CALL substart('ppm_util_cart2spher',t0,info)

      info       = 0
      r(1:n)     = 0.0_MK
      theta(1:n) = 0.0_MK
      phi(1:n)   = 0.0_MK

      !-------------------------------------------------------------------------
      !  Convert
      !-------------------------------------------------------------------------
      DO i=1,n
         dx  = x(i) - x0 
         dy  = y(i) - y0 
         dz  = z(i) - z0 
         rad = SQRT(dx*dx + dy*dy + dz*dz)
         IF (rad.GT.0.0) THEN
            r(i)     = rad
            phi(i)   = ATAN2(dy,dx)
            theta(i) = ACOS(dz/rad)
         ENDIF
      ENDDO
      
      !-------------------------------------------------------------------------
      !  Return 
      !-------------------------------------------------------------------------
 9999 CONTINUE
      CALL substop('ppm_util_cart2sph',t0,info)
      RETURN
#if   __KIND == __SINGLE_PRECISION
      END SUBROUTINE ppm_util_cart2sph_s
#elif __KIND == __DOUBLE_PRECISION
      END SUBROUTINE ppm_util_cart2sph_d
#endif
