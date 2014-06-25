      !-------------------------------------------------------------------------
      !  Subroutine   :                ppm_util_cart2sph
      !-------------------------------------------------------------------------
      ! Copyright (c) 2012 CSE Lab (ETH Zurich), MOSAIC Group (ETH Zurich),
      !                    Center for Fluid Dynamics (DTU)
      !
      !
      ! This file is part of the Parallel Particle Mesh Library (PPM).
      !
      ! PPM is free software: you can redistribute it and/or modify
      ! it under the terms of the GNU Lesser General Public License
      ! as published by the Free Software Foundation, either
      ! version 3 of the License, or (at your option) any later
      ! version.
      !
      ! PPM is distributed in the hope that it will be useful,
      ! but WITHOUT ANY WARRANTY; without even the implied warranty of
      ! MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
      ! GNU General Public License for more details.
      !
      ! You should have received a copy of the GNU General Public License
      ! and the GNU Lesser General Public License along with PPM. If not,
      ! see <http://www.gnu.org/licenses/>.
      !
      ! Parallel Particle Mesh Library (PPM)
      ! ETH Zurich
      ! CH-8092 Zurich, Switzerland
      !-------------------------------------------------------------------------

#if   __KIND == __SINGLE_PRECISION
      SUBROUTINE ppm_util_cart2sph_s(x,y,z,n,x0,y0,z0,r,theta,phi,info)
#elif __KIND == __DOUBLE_PRECISION
      SUBROUTINE ppm_util_cart2sph_d(x,y,z,n,x0,y0,z0,r,theta,phi,info)
#endif
      !!! Routine for conversion between Charthesian and
      !!! Spherical co-ordinates.
      !!!
      !!! [NOTE]
      !!! theta is angle of the vector with the z axis, and phi the
      !!! angle in the x-y plane.
      !!!
      !!! NOTE: The routine has by no means been optimized.
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
      !-------------------------------------------------------------------------
      !  Arguments
      !-------------------------------------------------------------------------
      REAL(MK), DIMENSION(:), INTENT(IN)  :: x
      !!! x coordinates
      REAL(MK), DIMENSION(:), INTENT(IN)  :: y
      !!! y coordinates
      REAL(MK), DIMENSION(:), INTENT(IN)  :: z
      !!! z coordinates
      REAL(MK),               INTENT(IN)  :: x0
      !!! x coordinate of center of coord. sys.
      REAL(MK),               INTENT(IN)  :: y0
      !!! y coordinate of center of coord. sys
      REAL(MK),               INTENT(IN)  :: z0
      !!! z coordinate of center of coord. sys
      REAL(MK), DIMENSION(:), INTENT(OUT) :: r
      !!! Radii in polar coord. sys.
      REAL(MK), DIMENSION(:), INTENT(OUT) :: theta
      !!! Theta-angle in polar coord. sys.
      REAL(MK), DIMENSION(:), INTENT(OUT) :: phi
      !!! Phi-angle in polar coord. sys.
      INTEGER, INTENT(IN)                :: n
      !!! Number of points
      INTEGER, INTENT(OUT)               :: info
      !!! Return status
      !-------------------------------------------------------------------------
      !  Local variables
      !-------------------------------------------------------------------------
      INTEGER :: i
      REAL(MK) :: dx,dy,dz,rad,t0

      !-------------------------------------------------------------------------
      !  Initialise
      !-------------------------------------------------------------------------
      CALL substart('ppm_util_cart2sph',t0,info)

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
