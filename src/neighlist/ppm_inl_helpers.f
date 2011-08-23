     !-------------------------------------------------------------------------
     !  Subroutines :               ppm_inl_helpers
     !-------------------------------------------------------------------------
     ! Copyright (c) 2010 CSE Lab (ETH Zurich), MOSAIC Group (ETH Zurich),
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
#if __ANISO == __YES
#if   __KIND == __SINGLE_PRECISION
      SUBROUTINE getSubdomainParticles_s_aniso(xp, Np, Mp, cutoff, lsymm,             &
     & actual_subdomain, ghost_extend, xp_sub, cutoff_sub, Np_sub, Mp_sub, p_id)
#elif __KIND == __DOUBLE_PRECISION
      SUBROUTINE getSubdomainParticles_d_aniso(xp, Np, Mp, cutoff, lsymm,             &
     & actual_subdomain, ghost_extend, xp_sub, cutoff_sub, Np_sub, Mp_sub, p_id)
#endif
#elif __ANISO == __NO
#if   __KIND == __SINGLE_PRECISION
      SUBROUTINE getSubdomainParticles_s(xp, Np, Mp, cutoff, lsymm,             &
     & actual_subdomain, ghost_extend, xp_sub, cutoff_sub, Np_sub, Mp_sub, p_id)
#elif __KIND == __DOUBLE_PRECISION
      SUBROUTINE getSubdomainParticles_d(xp, Np, Mp, cutoff, lsymm,             &
     & actual_subdomain, ghost_extend, xp_sub, cutoff_sub, Np_sub, Mp_sub, p_id)
#endif

#endif
          IMPLICIT NONE
#if   __KIND == __SINGLE_PRECISION
          INTEGER, PARAMETER :: mk = ppm_kind_single
#elif __KIND == __DOUBLE_PRECISION
          INTEGER, PARAMETER :: mk = ppm_kind_double
#endif
      !-------------------------------------------------------------------------
      !  Arguments
      !-------------------------------------------------------------------------
          REAL(MK), INTENT(IN), DIMENSION(:,:)       :: xp
          !!! Particle coordinates array. 
          !!! i.e., `xp(1, i)` is the x-coor of particle i.
          INTEGER , INTENT(IN)                       :: Np
          !!! Number of real particles
          INTEGER , INTENT(IN)                       :: Mp
          !!! Number of all particles including ghost particles
#if __ANISO == __YES
          REAL(MK), INTENT(IN), DIMENSION(:,:)       :: cutoff
          !!! Particle cutoff radii array, here the inverse transform for anisotropic particles
#elif __ANISO == __NO
          REAL(MK), INTENT(IN), DIMENSION(:)         :: cutoff
          !!! Particle cutoff radii array, here the radius for isotropic particles
#endif
          LOGICAL,  INTENT(IN)                       :: lsymm
          !!! If lsymm = TRUE, verlet lists are symmetric and we have ghost
          !!! layers only in (+) directions in all axes. Else, we have ghost
          !!! layers in all directions.
          REAL(MK),             DIMENSION(2*ppm_dim) :: actual_subdomain
          ! Physical extent of actual domain without ghost layers.
          REAL(MK), INTENT(IN), DIMENSION(ppm_dim)   :: ghost_extend
          !!! Extra area/volume over the actual domain introduced by
          !!! ghost layers.
          REAL(MK), POINTER   , DIMENSION(:,:)       :: xp_sub
          !!! Particle coordinates array. F.e., xp(1, i) is the x-coor of particle i.
#if __ANISO == __YES
          REAL(MK), DIMENSION(:,:), POINTER          :: cutoff_sub
          !!! Particle cutoff radii array, here the inverse transform for anisotropic particles
#elif __ANISO == __NO
          REAL(MK), DIMENSION(:), POINTER            :: cutoff_sub
          !!! Particle cutoff radii array, here the radius for isotropic particles
#endif
          INTEGER , INTENT(OUT)                      :: Np_sub
          !!! Number of real particles of subdomain
          INTEGER , INTENT(OUT)                      :: Mp_sub
          !!! Number of all particles including ghost particles of subdomain
          INTEGER , POINTER   , DIMENSION(:)        :: p_id
      !-------------------------------------------------------------------------
      !  Local variables, arrays and counters
      !-------------------------------------------------------------------------
          REAL(MK), DIMENSION(2*ppm_dim)             :: whole_subdomain
          !!! Physical extent of whole domain including ghost layers.
          INTEGER                                    :: i

      !---------------------------------------------------------------------
      !  Variables and parameters for ppm_alloc
      !---------------------------------------------------------------------
          INTEGER               :: iopt
          INTEGER, DIMENSION(2) :: lda
          INTEGER               :: info

          IF(lsymm .EQV. .TRUE.)  THEN
              DO i = 1, ppm_dim
                  whole_subdomain(2*i-1) = actual_subdomain(2*i-1)
                  whole_subdomain(2*i)   = actual_subdomain(2*i) + ghost_extend(i)
              END DO
          ELSE
              DO i = 1, ppm_dim
                  whole_subdomain(2*i-1) = actual_subdomain(2*i-1) - ghost_extend(i)
                  whole_subdomain(2*i)   = actual_subdomain(2*i)   + ghost_extend(i)
              END DO
          END IF

          Np_sub = 0
          Mp_sub = 0
          DO i = 1, Mp
              IF(inDomain(xp(:, i), actual_subdomain))  THEN
                   Np_sub = Np_sub + 1
                   Mp_sub = Mp_sub + 1
              ELSE
                  IF(inDomain(xp(:, i), whole_subdomain))  THEN
                      Mp_sub = Mp_sub + 1
                  END IF
              END IF
          END DO

          iopt   = ppm_param_alloc_fit
          lda(1) = ppm_dim
          lda(2) = Mp_sub
          CALL ppm_alloc(xp_sub, lda, iopt, info)
          IF (info.NE.0) THEN
              info = ppm_error_fatal
              CALL ppm_error(ppm_err_alloc,'ppm_create_subdomain_particles',     &
     &                       'own_plist',__LINE__,info)
          END IF

          ! haeckic changed 
#if __ANISO == __YES
          IF (size(cutoff(:,1)) .LT. 5) THEN
            lda(1) = 4
          ELSE
            lda(1) = 9
          ENDIF
          lda(2) = Mp_sub
          CALL ppm_alloc(cutoff_sub, lda, iopt, info)
#elif __ANISO == __NO
          lda(1) = Mp_sub
          CALL ppm_alloc(cutoff_sub, lda, iopt, info)
#endif
          IF (info.NE.0) THEN
              info = ppm_error_fatal
              CALL ppm_error(ppm_err_alloc,'ppm_create_subdomain_particles',     &
     &                       'own_plist',__LINE__,info)
          END IF
          
          lda(1) = Mp_sub
          CALL ppm_alloc(p_id, lda, iopt, info)
          IF (info.NE.0) THEN
              info = ppm_error_fatal
              CALL ppm_error(ppm_err_alloc,'ppm_create_subdomain_particles',     &
     &                       'own_plist',__LINE__,info)
          END IF

          Np_sub = 0
          DO i = 1, Mp
              IF(inDomain(xp(:, i), actual_subdomain))  THEN   !REAL
                   Np_sub = Np_sub + 1
                   p_id(Np_sub) = i
                   xp_sub(:, Np_sub) = xp(:, i)
#if __ANISO == __YES
                  IF (size(cutoff(:,1)) .LT. 5) THEN
                       cutoff_sub(1,Np_sub) = cutoff(1,i)
                       cutoff_sub(2,Np_sub) = cutoff(2,i)
                       cutoff_sub(3,Np_sub) = cutoff(3,i)
                       cutoff_sub(4,Np_sub) = cutoff(4,i)
                   ELSE
                       cutoff_sub(1,Np_sub) = cutoff(1,i)
                       cutoff_sub(2,Np_sub) = cutoff(2,i)
                       cutoff_sub(3,Np_sub) = cutoff(3,i)
                       cutoff_sub(4,Np_sub) = cutoff(4,i)
                       cutoff_sub(5,Np_sub) = cutoff(5,i)
                       cutoff_sub(6,Np_sub) = cutoff(6,i)
                       cutoff_sub(7,Np_sub) = cutoff(7,i)
                       cutoff_sub(8,Np_sub) = cutoff(8,i)
                       cutoff_sub(9,Np_sub) = cutoff(9,i)
                   ENDIF
#elif __ANISO == __NO
                   cutoff_sub(Np_sub) = cutoff(i)
#endif
              END IF
          END DO

          Mp_sub = Np_sub
          DO i = 1, Mp
              IF(inDomain(xp(:, i), actual_subdomain))  THEN   !REAL

              ELSE
                  IF(inDomain(xp(:, i), whole_subdomain))  THEN !GHOST
                      Mp_sub = Mp_sub + 1
                      p_id(Mp_sub) = i
                      xp_sub(:, Mp_sub) = xp(:, i)
#if __ANISO == __YES
                   IF (size(cutoff(:,1)) .LT. 5) THEN
                       cutoff_sub(1,Mp_sub) = cutoff(1,i)
                       cutoff_sub(2,Mp_sub) = cutoff(2,i)
                       cutoff_sub(3,Mp_sub) = cutoff(3,i)
                       cutoff_sub(4,Mp_sub) = cutoff(4,i)
                   ELSE
                       cutoff_sub(1,Mp_sub) = cutoff(1,i)
                       cutoff_sub(2,Mp_sub) = cutoff(2,i)
                       cutoff_sub(3,Mp_sub) = cutoff(3,i)
                       cutoff_sub(4,Mp_sub) = cutoff(4,i)
                       cutoff_sub(5,Mp_sub) = cutoff(5,i)
                       cutoff_sub(6,Mp_sub) = cutoff(6,i)
                       cutoff_sub(7,Mp_sub) = cutoff(7,i)
                       cutoff_sub(8,Mp_sub) = cutoff(8,i)
                       cutoff_sub(9,Mp_sub) = cutoff(9,i)
                   ENDIF
#elif __ANISO == __NO
                   cutoff_sub(Mp_sub) = cutoff(i)
#endif
                  END IF
              END IF
          END DO
#if __ANISO == __YES
#if   __KIND == __SINGLE_PRECISION
      END SUBROUTINE getSubdomainParticles_s_aniso
#elif __KIND == __DOUBLE_PRECISION
      END SUBROUTINE getSubdomainParticles_d_aniso
#endif
#elif __ANISO == __NO
#if   __KIND == __SINGLE_PRECISION
      END SUBROUTINE getSubdomainParticles_s
#elif __KIND == __DOUBLE_PRECISION
      END SUBROUTINE getSubdomainParticles_d
#endif
#endif

#if __ANISO == __YES
#if   __KIND == __SINGLE_PRECISION
      PURE FUNCTION inDomain_s(coor,domain) RESULT(isInDomain)
#elif __KIND == __DOUBLE_PRECISION
      PURE FUNCTION inDomain_d(coor,domain) RESULT(isInDomain)
#endif
      !!! This function checks whether given coordinates are inside whole
      !!! domain including ghost layers or not and RETURNs logical result.
          IMPLICIT NONE
#if   __KIND == __SINGLE_PRECISION
          INTEGER, PARAMETER :: mk = ppm_kind_single
#elif __KIND == __DOUBLE_PRECISION
          INTEGER, PARAMETER :: mk = ppm_kind_double
#endif
      !---------------------------------------------------------------------
      !  Arguments
      !---------------------------------------------------------------------
          REAL(MK), INTENT(IN), DIMENSION(ppm_dim)   :: coor
          !!! 1D Array of coordinates. First index is x-coordinate.
          REAL(MK), INTENT(IN), DIMENSION(2*ppm_dim) :: domain
          !!! Physical extent of whole domain including ghost layers.
          LOGICAL                                    :: isInDomain
          !!! Logical RETURN value. TRUE if inside the domain.

      !---------------------------------------------------------------------
      !  Counters
      !---------------------------------------------------------------------
          INTEGER                                    :: i

          ! Initialize to TRUE
          isInDomain = .TRUE.

          ! If any coordinate is out of bounds, set RETURN value to FALSE
          DO i = 1, ppm_dim
              IF((coor(i) .LT. domain(2*i-1)) .OR. (coor(i) .GE. domain(2*i)))  THEN
                  isInDomain = .FALSE.
                  RETURN ! No need to check further
              END IF
          END DO
#if   __KIND == __SINGLE_PRECISION
      END FUNCTION inDomain_s
#elif __KIND == __DOUBLE_PRECISION
      END FUNCTION inDomain_d
#endif
#endif

#if __ANISO == __YES
#if   __KIND == __SINGLE_PRECISION
      PURE FUNCTION isNeighbor_s_aniso(p_idx, p_neigh, xp, cutoff, skin) RESULT(isNeigh)
#elif __KIND == __DOUBLE_PRECISION
      PURE FUNCTION isNeighbor_d_aniso(p_idx, p_neigh, xp, cutoff, skin) RESULT(isNeigh)
#endif
#elif __ANISO == __NO
#if   __KIND == __SINGLE_PRECISION
      PURE FUNCTION isNeighbor_s(p_idx, p_neigh, xp, cutoff, skin) RESULT(isNeigh)
#elif __KIND == __DOUBLE_PRECISION
      PURE FUNCTION isNeighbor_d(p_idx, p_neigh, xp, cutoff, skin) RESULT(isNeigh)
#endif
#endif

      !!! Given indices of two particles, checks whether the euclidian distance
      !!! is smaller than the sum of minimum cutoff radius of these particles
      !!! and the skin, then RETURNs TRUE if so. Works for nD.
          IMPLICIT NONE
#if   __KIND == __SINGLE_PRECISION
          INTEGER, PARAMETER :: mk = ppm_kind_single
#elif __KIND == __DOUBLE_PRECISION
          INTEGER, PARAMETER :: mk = ppm_kind_double
#endif
      !---------------------------------------------------------------------
      !  Arguments
      !---------------------------------------------------------------------
          INTEGER,  INTENT(IN)                 :: p_idx
          !!! Index of first particle
          INTEGER,  INTENT(IN)                 :: p_neigh
          !!! Index of second particle
          REAL(MK), INTENT(IN), DIMENSION(:,:) :: xp
          !!! Coordinate array of particles
#if __ANISO == __YES
          REAL(MK), INTENT(IN), DIMENSION(:,:) :: cutoff
          !!! Particle cutoff radii array, here the inverse transform for anisotropic particles
#elif __ANISO == __NO
          REAL(MK), INTENT(IN), DIMENSION(:)   :: cutoff
          !!! Particle cutoff radii array, here the radius for isotropic particles
#endif
          REAL(MK), INTENT(IN)                 :: skin
          !!! Skin parameter
          LOGICAL                              :: isNeigh
          !!! Return value. TRUE if particles are neighbors

      !---------------------------------------------------------------------
      !  Local variables and counters
      !---------------------------------------------------------------------
#if __ANISO == __YES
          REAL(MK),DIMENSION(ppm_dim)                :: total_dist,dist_transformed
          ! Euclidian distance between particles
#elif __ANISO == __NO
          REAL(MK)                             :: total_dist
          ! Euclidian distance between particles
#endif
          REAL(MK)                             :: rcutoff,rcutoff2
          ! Minimum cutoff radius, plus skin.
          INTEGER                              :: i
          ! Counter

          ! Initialize RETURN value to FALSE
          isNeigh = .FALSE.

          ! Return FALSE if they are the same particle
          IF(p_idx .EQ. p_neigh)    RETURN

! haeckic change anisotropic
#if __ANISO == __YES

          ! get distance
          DO i = 1, ppm_dim
              total_dist(i) = xp(i, p_neigh) - xp(i, p_idx)
          END DO

          ! transform
          IF (ppm_dim .EQ. 2) THEN

            dist_transformed(1) = cutoff(1,p_idx)*total_dist(1) + cutoff(2,p_idx)*total_dist(2)
            dist_transformed(2) = cutoff(3,p_idx)*total_dist(1) + cutoff(4,p_idx)*total_dist(2)
            rcutoff = sqrt(dist_transformed(1)**2 + dist_transformed(2)**2)

          ELSE

            dist_transformed(1) = cutoff(1,p_idx)*total_dist(1) + cutoff(2,p_idx)*total_dist(2) &
            & + cutoff(3,p_idx)*total_dist(3)
            dist_transformed(2) = cutoff(4,p_idx)*total_dist(1) + cutoff(5,p_idx)*total_dist(2) &
            & + cutoff(6,p_idx)*total_dist(3)
            dist_transformed(3) = cutoff(7,p_idx)*total_dist(1) + cutoff(8,p_idx)*total_dist(2) & 
            & + cutoff(9,p_idx)*total_dist(3)
            rcutoff = sqrt(dist_transformed(1)**2 + dist_transformed(2)**2 + dist_transformed(3)**2)

          ENDIF

          ! transform
          IF (ppm_dim .EQ. 2) THEN

            dist_transformed(1) = cutoff(1,p_neigh)*total_dist(1) + cutoff(2,p_neigh)*total_dist(2)
            dist_transformed(2) = cutoff(3,p_neigh)*total_dist(1) + cutoff(4,p_neigh)*total_dist(2)
            rcutoff2 = sqrt(dist_transformed(1)**2 + dist_transformed(2)**2)

          ELSE

            dist_transformed(1) = cutoff(1,p_neigh)*total_dist(1) + cutoff(2,p_neigh)*total_dist(2) & 
            & + cutoff(3,p_neigh)*total_dist(3)
            dist_transformed(2) = cutoff(4,p_neigh)*total_dist(1) + cutoff(5,p_neigh)*total_dist(2) &
            & + cutoff(6,p_neigh)*total_dist(3)
            dist_transformed(3) = cutoff(7,p_neigh)*total_dist(1) + cutoff(8,p_neigh)*total_dist(2) &
            & + cutoff(9,p_neigh)*total_dist(3)
            rcutoff2 = sqrt(dist_transformed(1)**2 + dist_transformed(2)**2 + dist_transformed(3)**2)

          ENDIF
          
          ! Return TRUE if they are neighbors, i.e. both below 1.0
          IF(rcutoff .LE. 1.0_mk .AND. rcutoff2 .LE. 1.0_mk)  isNeigh = .TRUE.


#elif __ANISO == __NO
          ! Initialize euclidian distance to 0
          total_dist = 0
          ! Add squares of distances on each axis, then take square root of it
          total_dist = total_dist + (xp(1, p_idx) - xp(1, p_neigh))**2
          total_dist = total_dist + (xp(2, p_idx) - xp(2, p_neigh))**2
          DO i = 3, ppm_dim
              total_dist = total_dist + (xp(i, p_idx) - xp(i, p_neigh))**2
          END DO
          total_dist = sqrt(total_dist)

          ! Pick smallest cutoff radius and add skin on it
          rcutoff = MIN(cutoff(p_idx), cutoff(p_neigh)) + skin

          ! Return TRUE if they are neighbors
          IF(total_dist .LE. rcutoff)  isNeigh = .TRUE.
#endif

#if __ANISO == __YES

#if   __KIND == __SINGLE_PRECISION
      END FUNCTION isNeighbor_s_aniso
#elif __KIND == __DOUBLE_PRECISION
      END FUNCTION isNeighbor_d_aniso
#endif

#elif __ANISO == __NO

#if   __KIND == __SINGLE_PRECISION
      END FUNCTION isNeighbor_s
#elif __KIND == __DOUBLE_PRECISION
      END FUNCTION isNeighbor_d
#endif

#endif


#if __ANISO == __YES

#if   __KIND == __SINGLE_PRECISION
      PURE FUNCTION is_xset_Neighbor_s_aniso(red_idx, blue_idx, red, rcred, blue, &
 &                  rcblue, skin) RESULT(isNeigh)
#elif __KIND == __DOUBLE_PRECISION
      PURE FUNCTION is_xset_Neighbor_d_aniso(red_idx, blue_idx, red, rcred, blue, &
 &                  rcblue, skin) RESULT(isNeigh)
#endif

#elif __ANISO == __NO

#if   __KIND == __SINGLE_PRECISION
      PURE FUNCTION is_xset_Neighbor_s(red_idx, blue_idx, red, rcred, blue, &
 &                  rcblue, skin) RESULT(isNeigh)
#elif __KIND == __DOUBLE_PRECISION
      PURE FUNCTION is_xset_Neighbor_d(red_idx, blue_idx, red, rcred, blue, &
 &                  rcblue, skin) RESULT(isNeigh)
#endif

#endif
      !!! Given indices of two particles, checks whether the euclidian distance
      !!! is smaller than the sum of minimum cutoff radius of these particles
      !!! and the skin, then RETURNs TRUE if so. Works for nD.
          IMPLICIT NONE
#if   __KIND == __SINGLE_PRECISION
          INTEGER, PARAMETER :: mk = ppm_kind_single
#elif __KIND == __DOUBLE_PRECISION
          INTEGER, PARAMETER :: mk = ppm_kind_double
#endif
      !---------------------------------------------------------------------
      !  Arguments
      !---------------------------------------------------------------------
          INTEGER,  INTENT(IN)                 :: red_idx
          !!! Index of first particle
          INTEGER,  INTENT(IN)                 :: blue_idx
          !!! Index of second particle
          REAL(MK), INTENT(IN), DIMENSION(:,:) :: red
          !!! Coordinate array of particles (red)
#if __ANISO == __YES
          REAL(MK), INTENT(IN), DIMENSION(:,:) :: rcred
          !!! Blue particle cutoff radii array, in the anisotropic case
#elif __ANISO == __NO
          REAL(MK), INTENT(IN), DIMENSION(:)   :: rcred
          !!! Blue particle cutoff radii array
#endif
          REAL(MK), INTENT(IN), DIMENSION(:,:) :: blue
          !!! Coordinate array of particles (blue)
#if __ANISO == __YES
          REAL(MK), INTENT(IN), DIMENSION(:,:) :: rcblue
          !!! Blue particle cutoff radii array, in the anisotropic case
#elif __ANISO == __NO
          REAL(MK), INTENT(IN), DIMENSION(:)   :: rcblue
          !!! Blue particle cutoff radii array
#endif
          REAL(MK), INTENT(IN)                 :: skin
          !!! Skin parameter
          LOGICAL                              :: isNeigh
          !!! Return value. TRUE if particles are neighbors

      !---------------------------------------------------------------------
      !  Local variables and counters
      !---------------------------------------------------------------------
#if __ANISO == __YES
          REAL(MK),DIMENSION(ppm_dim)          :: total_dist,dist_transformed
          ! Euclidian distance between particles
#elif __ANISO == __NO
          REAL(MK)                             :: total_dist
          ! Euclidian distance between particles
#endif
          REAL(MK)                             :: rcutoff,rcutoff2
          ! Minimum cutoff radius, plus skin.
          INTEGER                              :: i
          ! Counter

          ! Initialize RETURN value to FALSE
          isNeigh = .FALSE.

! haeckic change anistropic
#if __ANISO == __YES

          ! get distance
          DO i = 1, ppm_dim
              total_dist(i) = red(i, red_idx) - blue(i, blue_idx)
          END DO

          ! transform
          IF (ppm_dim .EQ. 2) THEN

            dist_transformed(1) = rcred(1,red_idx)*total_dist(1) + rcred(2,red_idx)*total_dist(2)
            dist_transformed(2) = rcred(3,red_idx)*total_dist(1) + rcred(4,red_idx)*total_dist(2)
            rcutoff = sqrt(dist_transformed(1)**2 + dist_transformed(2)**2)

          ELSE

            dist_transformed(1) = rcred(1,red_idx)*total_dist(1) + rcred(2,red_idx)*total_dist(2) &
            & + rcred(3,red_idx)*total_dist(3)
            dist_transformed(2) = rcred(4,red_idx)*total_dist(1) + rcred(5,red_idx)*total_dist(2) &
            & + rcred(6,red_idx)*total_dist(3)
            dist_transformed(3) = rcred(7,red_idx)*total_dist(1) + rcred(8,red_idx)*total_dist(2) & 
            & + rcred(9,red_idx)*total_dist(3)
            rcutoff = sqrt(dist_transformed(1)**2 + dist_transformed(2)**2 + dist_transformed(3)**2)

          ENDIF

          ! transform
          IF (ppm_dim .EQ. 2) THEN

            dist_transformed(1) = rcblue(1,blue_idx)*total_dist(1) + rcblue(2,blue_idx)*total_dist(2)
            dist_transformed(2) = rcblue(3,blue_idx)*total_dist(1) + rcblue(4,blue_idx)*total_dist(2)
            rcutoff2 = sqrt(dist_transformed(1)**2 + dist_transformed(2)**2)

          ELSE

            dist_transformed(1) = rcblue(1,blue_idx)*total_dist(1) + rcblue(2,blue_idx)*total_dist(2) & 
            & + rcblue(3,blue_idx)*total_dist(3)
            dist_transformed(2) = rcblue(4,blue_idx)*total_dist(1) + rcblue(5,blue_idx)*total_dist(2) &
            & + rcblue(6,blue_idx)*total_dist(3)
            dist_transformed(3) = rcblue(7,blue_idx)*total_dist(1) + rcblue(8,blue_idx)*total_dist(2) &
            & + rcblue(9,blue_idx)*total_dist(3)
            rcutoff2 = sqrt(dist_transformed(1)**2 + dist_transformed(2)**2 + dist_transformed(3)**2)

          ENDIF

          ! Return TRUE if they are neighbors, i.e. both below 1.0
          IF(rcutoff .LE. 1.0_mk .AND. rcutoff2 .LE. 1.0_mk)  isNeigh = .TRUE.
#elif __ANISO == __NO
          ! Initialize euclidian distance to 0
          total_dist = 0
          ! Add squares of distances on each axis, then take square root of it
          total_dist = total_dist + (red(1, red_idx) - blue(1, blue_idx))**2
          total_dist = total_dist + (red(2, red_idx) - blue(2, blue_idx))**2
          DO i = 3, ppm_dim
              total_dist = total_dist + (red(i, red_idx) - blue(i, blue_idx))**2
          END DO
          total_dist = sqrt(total_dist)

          ! Pick smallest cutoff radius and add skin on it
          rcutoff = MIN(rcred(red_idx), rcblue(blue_idx)) + skin

          ! Return TRUE if they are neighbors
          IF(total_dist .LE. rcutoff)  isNeigh = .TRUE.
#endif

#if __ANISO == __YES

#if   __KIND == __SINGLE_PRECISION
      END FUNCTION is_xset_Neighbor_s_aniso
#elif __KIND == __DOUBLE_PRECISION
      END FUNCTION is_xset_Neighbor_d_aniso
#endif

#elif __ANISO == __NO
      
#if   __KIND == __SINGLE_PRECISION
      END FUNCTION is_xset_Neighbor_s
#elif __KIND == __DOUBLE_PRECISION
      END FUNCTION is_xset_Neighbor_d
#endif

#endif

#if __ANISO == __YES
#if   __KIND == __SINGLE_PRECISION
      FUNCTION cross_neighbor_s(p_idx, p_neigh, xp, actual_domain) RESULT(isCrossNeigh)
#elif __KIND == __DOUBLE_PRECISION
      FUNCTION cross_neighbor_d(p_idx, p_neigh, xp, actual_domain) RESULT(isCrossNeigh)
#endif
      !!! This function is to be used for lsymm = TRUE case, where we have
      !!! ghost layers only in positive direction which makes it necessary
      !!! to take ghost-ghost interaction into account. So, ghost particles
      !!! that are cross-neighbors should be detected and put in verlet list
      !!! of one of them.
          IMPLICIT NONE
#if   __KIND == __SINGLE_PRECISION
          INTEGER, PARAMETER :: mk = ppm_kind_single
#elif __KIND == __DOUBLE_PRECISION
          INTEGER, PARAMETER :: mk = ppm_kind_double
#endif
      !---------------------------------------------------------------------
      !  Arguments
      !---------------------------------------------------------------------
          INTEGER,  INTENT(IN)                       :: p_idx
          !!! Index of first particle
          INTEGER,  INTENT(IN)                       :: p_neigh
          !!! Index of second particle
          REAL(MK), INTENT(IN), DIMENSION(:,:)       :: xp
          !!! Coordinates of particles
          REAL(MK), INTENT(IN), DIMENSION(2*ppm_dim) :: actual_domain
          !!! Physical extent of actual domain without ghost layers
          LOGICAL                                    :: isCrossNeigh
          !!! Return value. TRUE if particles are cross-neighbors

      !-------------------------------------------------------------------------
      !  Local variables and counters
      !-------------------------------------------------------------------------
          INTEGER                                    :: region1
          ! To form bitwise representation of first particle
          INTEGER                                    :: region2
          ! To form bitwise representation of second particle

          !---------------------------------------------------------------------
          !  Initialize return value to FALSE and local variables to 0
          !---------------------------------------------------------------------
          isCrossNeigh  = .FALSE.
          region1 = 0
          region2 = 0

          !---------------------------------------------------------------------
          !  For first two dimensions, loops are unrolled. If the particle is
          !  out of bounds of actual domain on an axis, it is set to 1. For the
          !  next axis, bits are shifted to LEFT and new bit is added, 0 if
          !  inside actual domain or 1 otherwise.
          !---------------------------------------------------------------------
          IF(xp(1, p_idx)   .GT. actual_domain(2))  region1 = 1
          IF(xp(1, p_neigh) .GT. actual_domain(2))  region2 = 1
          region1 = ISHFT(region1, 1)
          region2 = ISHFT(region2, 1)
          IF(xp(2, p_idx)   .GT. actual_domain(4))  region1 = IOR(region1, 1)
          IF(xp(2, p_neigh) .GT. actual_domain(4))  region2 = IOR(region2, 1)

          !---------------------------------------------------------------------
          !  If ppm_dim = 3, then we also compute region1 and region2 for z-axis.
          !---------------------------------------------------------------------
          IF(ppm_dim .EQ. 3)    THEN
              region1 = ISHFT(region1, 1)
              region2 = ISHFT(region2, 1)
              IF(xp(3, p_idx)   .GT. actual_domain(6))  region1 = IOR(region1, 1)
              IF(xp(3, p_neigh) .GT. actual_domain(6))  region2 = IOR(region2, 1)
          END IF

          !---------------------------------------------------------------------
          !  Set return value to TRUE if they are cross-neighbors.
          !---------------------------------------------------------------------
          IF(IAND(region1, region2) .EQ. 0)   isCrossNeigh = .TRUE.
#if   __KIND == __SINGLE_PRECISION
      END FUNCTION cross_neighbor_s
#elif __KIND == __DOUBLE_PRECISION
      END FUNCTION cross_neighbor_d
#endif

#endif

#if __ANISO == __YES

#if   __KIND == __SINGLE_PRECISION
      SUBROUTINE putInEmptyList(c_idx)
      !!! Given the index of the cell, this subroutine stores the cell index
      !!! in the empty list. In case of insufficient list size, allocates
      !!! a larger array and copies the content within into the new array.
          IMPLICIT NONE
      !-------------------------------------------------------------------------
      !  Arguments
      !-------------------------------------------------------------------------
          INTEGER(ppm_kind_int64), INTENT(IN) :: c_idx
          !!! Index of the cell

      !-------------------------------------------------------------------------
      !  Local variables and counters
      !-------------------------------------------------------------------------
          INTEGER                             :: info
          ! If operation is successful, it is set to 0.
          INTEGER                             :: iopt
          INTEGER, DIMENSION(1)               :: lda

          !---------------------------------------------------------------------
          !  If the empty_list array is full, grow the empty_list array while
          !  preserving its contents.
          !---------------------------------------------------------------------
          IF(empty_pos .EQ. size(empty_list)) THEN
              lda(1) = 2*size(empty_list)
              iopt   = ppm_param_alloc_grow_preserve
              CALL ppm_alloc(empty_list, lda, iopt, info)
          END IF

          !---------------------------------------------------------------------
          !  Add the cell index into empty_list array.
          !---------------------------------------------------------------------
          empty_pos = empty_pos + 1
          empty_list(empty_pos) = c_idx
      END SUBROUTINE putInEmptyList
#endif

#if   __KIND == __SINGLE_PRECISION
      PURE FUNCTION inEmptyList(c_idx)  RESULT(inside)
      !!! Given the cell index, this function checks whether the cell is
      !!! in empty list or not and returns TRUE if inside.
          IMPLICIT NONE
      !-------------------------------------------------------------------------
      !  Arguments
      !-------------------------------------------------------------------------
          INTEGER(ppm_kind_int64), INTENT(IN) :: c_idx
          !!! Cell index
          LOGICAL                             :: inside
          !!! Return value. TRUE if the cell is in empty list

      !-------------------------------------------------------------------------
      !  Local variables and counters
      !-------------------------------------------------------------------------
          INTEGER :: pos ! position on empty_list array

          !---------------------------------------------------------------------
          !  Initialize return value to FALSE, then start searching from
          !  the end since it is highly possible that the cell would be
          !  inserted close to the end.
          !---------------------------------------------------------------------
          inside = .FALSE.
          DO pos = empty_pos,1,-1
              IF(empty_list(pos) .EQ. c_idx)   THEN
                  inside = .TRUE.
                  RETURN
              END IF
          END DO
      END FUNCTION inEmptyList
#endif

#endif

#if __ANISO == __YES
#if   __KIND == __SINGLE_PRECISION
      SUBROUTINE getParticleCoorDepth_s_aniso(p_idx, domain, p_coor, p_depth, xp, cutoff, skin)
#elif   __KIND == __DOUBLE_PRECISION
      SUBROUTINE getParticleCoorDepth_d_aniso(p_idx, domain, p_coor, p_depth, xp, cutoff, skin)
#endif
#elif __ANISO == __NO
#if   __KIND == __SINGLE_PRECISION
      SUBROUTINE getParticleCoorDepth_s(p_idx, domain, p_coor, p_depth, xp, cutoff, skin)
#elif   __KIND == __DOUBLE_PRECISION
      SUBROUTINE getParticleCoorDepth_d(p_idx, domain, p_coor, p_depth, xp, cutoff, skin)
#endif
#endif


      !!! Given the particle index, this subroutine modifies p_coor array and
      !!! p_depth variables such that they contain coordinates and depth of the
      !!! particle, respectively.
          IMPLICIT NONE
#if   __KIND == __SINGLE_PRECISION
          INTEGER, PARAMETER :: mk = ppm_kind_single
#elif __KIND == __DOUBLE_PRECISION
          INTEGER, PARAMETER :: mk = ppm_kind_double
#endif
      !-------------------------------------------------------------------------
      !  Arguments
      !-------------------------------------------------------------------------
          INTEGER,  INTENT(IN)                 :: p_idx
          REAL(MK), DIMENSION(2*ppm_dim)       :: domain
          REAL(MK), DIMENSION(:)               :: p_coor
          INTEGER,  INTENT(INOUT)              :: p_depth
          REAL(MK), INTENT(IN), DIMENSION(:,:) :: xp
#if __ANISO == __YES
          REAL(MK), INTENT(IN), DIMENSION(:,:) :: cutoff
          !!! Particle cutoff radii array, here the inverse transform for anisotropic particles
#elif __ANISO == __NO
          REAL(MK), INTENT(IN), DIMENSION(:)   :: cutoff
          !!! Particle cutoff radii array, here the radius for isotropic particles
#endif
          REAL(MK), INTENT(IN)                 :: skin

      !-------------------------------------------------------------------------
      !  Local variables and counters
      !-------------------------------------------------------------------------
          REAL(MK)                             :: minSideLength,cu
          INTEGER                              :: i

          ! Get maximum side length of the domain
          minSideLength = getMinimumSideLength(domain)

          ! Initialize p_depth to 0 and keep on incrementing until we reach the
          ! correct depth.
          p_depth = 0

         ! haeckic change here
#if __ANISO == __YES
          cu = particles_longer_axis(cutoff(:,p_idx))
          DO WHILE(minSideLength .GT. (cu + skin))
#elif __ANISO == __NO
          DO WHILE(minSideLength .GT. (cutoff(p_idx) + skin))
#endif
              minSideLength = minSideLength/2
              p_depth = p_depth + 1
          END DO

!         Ferit: I believe the operation above is faster than two log operations as
!                there is no intrinsic LOG operation in base 2.
!         p_depth = CEILING(LOG(minSideLength/cutoff(p_idx))/LOG(2.0))

          ! Modify p_coor array such that it contains particle coordinates.
          p_coor(1) = xp(1, p_idx)
          p_coor(2) = xp(2, p_idx)
          DO i = 3, ppm_dim
              p_coor(i) = xp(i, p_idx)
          END DO
#if __ANISO == __YES

#if   __KIND == __SINGLE_PRECISION
      END SUBROUTINE getParticleCoorDepth_s_aniso
#elif   __KIND == __DOUBLE_PRECISION
      END SUBROUTINE getParticleCoorDepth_d_aniso
#endif

#elif __ANISO == __NO

#if   __KIND == __SINGLE_PRECISION
      END SUBROUTINE getParticleCoorDepth_s
#elif   __KIND == __DOUBLE_PRECISION
      END SUBROUTINE getParticleCoorDepth_d
#endif

#endif

#if __ANISO == __YES

#if   __KIND == __SINGLE_PRECISION
      SUBROUTINE getParticlesInCell_s(cell_idx, xp, clist, list, nlist)
#elif   __KIND == __DOUBLE_PRECISION
      SUBROUTINE getParticlesInCell_d(cell_idx, xp, clist, list, nlist)
#endif
      !!! Given the cell index, this subroutine modifies the list array such
      !!! that it contains the particle IDs of this cell and sets nlist to
      !!! number of particles in this cell.
          IMPLICIT NONE
#if   __KIND == __SINGLE_PRECISION
          INTEGER, PARAMETER :: mk = ppm_kind_single
#elif __KIND == __DOUBLE_PRECISION
          INTEGER, PARAMETER :: mk = ppm_kind_double
#endif
      !-------------------------------------------------------------------------
      !  Arguments
      !-------------------------------------------------------------------------
          INTEGER(ppm_kind_int64),  INTENT(IN)                   :: cell_idx
          REAL(MK),  DIMENSION(:,:),INTENT(IN)                   :: xp
          !!! this is basically a dummy argument to force fortran to generate
          !!! two versions of this routine
          TYPE(ppm_clist),          INTENT(IN)                   :: clist
          INTEGER,   DIMENSION(:),  INTENT(INOUT)                :: list
          INTEGER,                  INTENT(INOUT)                :: nlist

      !-------------------------------------------------------------------------
      !  Local variables and counters
      !-------------------------------------------------------------------------
          INTEGER(ppm_kind_int64)                      :: parentIdx
          INTEGER                                      :: left_end
          INTEGER                                      :: right_end
          INTEGER                                      :: border_idx
          INTEGER                                      :: i

          ! Get index of parent of this cell
          parentIdx  = parent(cell_idx)

          ! Get position on borders array that the parent is located in.
          border_idx = hash_search(clist%lookup,parentIdx)

          ! Initialize number of particles to 0.
          nlist = 0

          ! If this cell is not found in hash table, then put the cell index in
          ! empty list and return.
          IF(border_idx .EQ. htable_null)  THEN
              CALL putInEmptyList(cell_idx)
              RETURN
          END IF

          ! For 2D case
          IF(ppm_dim .EQ. 2)    THEN
              ! If the cell does not contain any particles that are in deeper
              ! levels in its region ...
              IF(clist%borders(6, border_idx) .EQ. 1)  THEN
                  ! Put it in empty list
                  CALL putInEmptyList(cell_idx)
              END IF

              ! Get index of first column on borders array.
              left_end  = 1 + cell_idx - (4*parentIdx-2)

              ! Get index of last column on borders array.
              right_end = left_end + 1

              ! If this is the top level, then get all particles
              IF(parentIdx .EQ. 0)  then
                   left_end  = 1
                   right_end = 5
              END IF
          ! For 3D case
          ELSEIF(ppm_dim .EQ. 3)   THEN
              ! If the cell does not contain any particles that are in deeper
              ! levels in its region ...
              IF(clist%borders(10, border_idx) .EQ. 1)  THEN
                  ! Put it in empty list
                  CALL putInEmptyList(cell_idx)
              END IF

              ! Get index of first column on borders array.
              left_end  = 1 + cell_idx - (8*parentIdx-6)

              ! Get index of last column on borders array.
              right_end = left_end + 1

              ! If this is the top level, then get all particles
              IF(parentIdx .EQ. 0)  then
                   left_end  = 1
                   right_end = 9
              END IF
          END IF

          ! From first column to last, get all particles and put them in the list
          DO i = (clist%borders(left_end, border_idx) + 1), &
 &                clist%borders(right_end, border_idx)
              nlist = nlist + 1
              list(nlist) = clist%rank(i)
          END DO
#if   __KIND == __SINGLE_PRECISION
      END SUBROUTINE getParticlesInCell_s
#elif __KIND == __DOUBLE_PRECISION
      END SUBROUTINE getParticlesInCell_d
#endif

#endif