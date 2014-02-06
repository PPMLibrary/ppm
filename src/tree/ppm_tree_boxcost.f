      !-------------------------------------------------------------------------
      !  Subroutine   :                 ppm_tree_boxcost
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
      SUBROUTINE ppm_tree_boxcost_s(Nm,weights,min_box,max_box,   &
     &    nbox,lhbx,lpdx,boxcost,info,pcost)
#elif __KIND == __DOUBLE_PRECISION
      SUBROUTINE ppm_tree_boxcost_d(Nm,weights,min_box,max_box,   &
     &    nbox,lhbx,lpdx,boxcost,info,pcost)
#endif
      !!! This routine computes the cost of a set of boxes.
      !!!
      !!! [NOTE]
      !!! The cost of a particle is counted as pcost (if
      !!! given) or 1. For meshes, the cost is 1 per mesh
      !!! point. For the geometry part, the cost of a
      !!! box is given by its volume. Use weights(3) to
      !!! adjust this if needed.
      !-------------------------------------------------------------------------
      !  Modules
      !-------------------------------------------------------------------------
      USE ppm_module_data
      USE ppm_module_data_tree
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
#ifdef __MPI
      INCLUDE 'mpif.h'
#endif
      !-------------------------------------------------------------------------
      !  Arguments
      !-------------------------------------------------------------------------
      REAL(MK), DIMENSION(:,:), INTENT(IN   ) :: min_box
      !!! The minimum coordinate of the boxes.
      !!!
      !!! 1st index: x,y[,z]                                                   +
      !!! 2nd: box ID
      REAL(MK), DIMENSION(:,:), INTENT(IN   ) :: max_box
      !!! The maximum coordinate of the boxes.
      !!!
      !!! 1st index: x,y[,z]                                                   +
      !!! 2nd: box ID
      REAL(MK), DIMENSION(3  ), INTENT(IN   ) :: weights
      !!! Weights for the three cost contributions: particles, mesh,
      !!! geometry for defining the total cost of a box.
      REAL(MK), DIMENSION(:  ), OPTIONAL, INTENT(IN) :: pcost
      !!! Argument of length Npart, specifying the
      !!! computational cost of each particle.
      INTEGER                 , INTENT(IN   ) :: nbox
      !!! Number of boxes to compute cost for.
      INTEGER , DIMENSION(:  ), INTENT(IN   ) :: lhbx
      !!! Pointer to first particle in each box to determine cost for.
      !!! This only needs to be present if there are any particles.
      INTEGER , DIMENSION(:  ), INTENT(IN   ) :: lpdx
      !!! Index list of points in each box to determine cost for.
      !!! This only needs to be present if there are any particles.
      INTEGER , DIMENSION(:,:), INTENT(IN   ) :: Nm
      !!! Number of grid points in each box.
      REAL(MK), DIMENSION(:  ), POINTER       :: boxcost
      !!! Costs of the boxes 1..nbox. boxcost(i) is the cost of box  boxes(i).
      INTEGER                 , INTENT(  OUT) :: info
      !!! Return status, 0 on success
      !-------------------------------------------------------------------------
      !  Local variables
      !-------------------------------------------------------------------------
      REAL(MK), DIMENSION(ppm_dim)            :: len_box
      INTEGER , DIMENSION(2)                  :: ldc
      REAL(MK)                                :: t0,meshtotal,geomtotal,dm
      INTEGER                                 :: i,ip,iopt,j
      REAL(MK), DIMENSION(:), POINTER         :: pcst
#ifdef __MPI
      REAL(MK), DIMENSION(:), POINTER         :: pcsum
      INTEGER                                 :: MPTYPE
#endif
      !-------------------------------------------------------------------------
      !  Externals
      !-------------------------------------------------------------------------

      !-------------------------------------------------------------------------
      !  Initialise
      !-------------------------------------------------------------------------
      CALL substart('ppm_tree_boxcost',t0,info)

      !-------------------------------------------------------------------------
      !  Check input arguments
      !-------------------------------------------------------------------------
      IF (ppm_debug .GT. 0) THEN
        CALL check
        IF (info .NE. 0) GOTO 9999
      ENDIF

      !-------------------------------------------------------------------------
      !  Set pointer to work memory
      !-------------------------------------------------------------------------
#if   __KIND == __SINGLE_PRECISION
      pcst => pcst_s
#else
      pcst => pcst_d
#endif
#ifdef   __MPI
#if   __KIND == __SINGLE_PRECISION
      pcsum => pcsum_s
#else
      pcsum => pcsum_d
#endif
#endif

      !-------------------------------------------------------------------------
      !  If we have less than 1 box to compute, exit
      !-------------------------------------------------------------------------
      IF (nbox .LT. 1) THEN
          IF (ppm_debug .GT. 0) THEN
              CALL ppm_write(ppm_rank,'ppm_tree_boxcost',   &
     &            'No boxes to be computed. Exiting.',info)
          ENDIF
          GOTO 9999
      ENDIF

#ifdef __MPI
      !-------------------------------------------------------------------------
      !  Determine MPI data type
      !-------------------------------------------------------------------------
#if   __KIND == __SINGLE_PRECISION
      MPTYPE = MPI_REAL
#elif __KIND == __DOUBLE_PRECISION
      MPTYPE = MPI_DOUBLE_PRECISION
#endif
#endif

      !-------------------------------------------------------------------------
      !  Computation of particle contributions to box costs
      !-------------------------------------------------------------------------
      IF (have_particles .AND. weights(1) .NE. 0.0_MK) THEN
          DO i=1,nbox
              pcst(i)  = 0.0_MK
#ifdef __MPI
              pcsum(i) = 0.0_MK
#endif
              DO j=lhbx(i),lhbx(i+1)-1
                  ip = lpdx(j)
                  dm = 1.0_MK
                  IF (PRESENT(pcost)) dm = pcost(ip)
                  pcst(i) = pcst(i) + dm
              ENDDO
          ENDDO

#ifdef __MPI
          !---------------------------------------------------------------------
          !  Allreduce of all particle sums
          !---------------------------------------------------------------------
          CALL MPI_Allreduce(pcst,pcsum,nbox,MPTYPE,MPI_SUM,ppm_comm,info)
          IF (info.NE.0) THEN
              info = ppm_error_fatal
              CALL ppm_error(ppm_err_mpi_fail,'ppm_tree_boxcost',   &
     &            'MPI_Allreduce of particle costs',__LINE__,info)
              GOTO 9999
          ENDIF
          DO i=1,nbox
              pcst(i) = pcsum(i)
          ENDDO
#endif
      ENDIF   ! have_particles

      !-------------------------------------------------------------------------
      !  Add up all cost contributions to total cost.
      !-------------------------------------------------------------------------
      IF (ppm_dim .EQ. 2) THEN
          DO i=1,nbox
              !-----------------------------------------------------------------
              !  Mesh and geometry based costs
              !-----------------------------------------------------------------
              meshtotal = 0.0_MK
              len_box(1) = max_box(1,i)-min_box(1,i)
              len_box(2) = max_box(2,i)-min_box(2,i)
              IF (have_mesh .AND. weights(2) .NE. 0) THEN
                  meshtotal = REAL(Nm(1,i),MK)*REAL(Nm(2,i),MK)
              ENDIF
              geomtotal = len_box(1)*len_box(2)
              !-----------------------------------------------------------------
              !  Compute weighted sum of costs and store it
              !-----------------------------------------------------------------
              IF(have_particles) THEN
                 boxcost(i) = pcst(i)*weights(1) + meshtotal*weights(2) +    &
     &                            geomtotal*weights(3)
              ELSE
                 boxcost(i) = meshtotal*weights(2) + geomtotal*weights(3)
              END IF
          ENDDO
      ELSE
          DO i=1,nbox
              !-----------------------------------------------------------------
              !  Mesh and geometry based costs
              !-----------------------------------------------------------------
              meshtotal = 0.0_MK
              len_box(1) = max_box(1,i)-min_box(1,i)
              len_box(2) = max_box(2,i)-min_box(2,i)
              len_box(3) = max_box(3,i)-min_box(3,i)
              IF (have_mesh .AND. weights(2) .NE. 0) THEN
                  meshtotal = REAL(Nm(1,i),MK)*REAL(Nm(2,i),MK)*REAL(Nm(3,i),MK)
              ENDIF
              geomtotal = len_box(1)*len_box(2)*len_box(3)
              !-----------------------------------------------------------------
              !  Compute weighted sum of costs and store it
              !-----------------------------------------------------------------
              IF(have_particles) THEN
                 boxcost(i) = pcst(i)*weights(1) + meshtotal*weights(2) +    &
     &                            geomtotal*weights(3)
              ELSE
                 boxcost(i) = meshtotal*weights(2) + geomtotal*weights(3)
              END IF
          ENDDO
      ENDIF

      !-------------------------------------------------------------------------
      !  Free local memory
      !-------------------------------------------------------------------------
 9999 CONTINUE
      NULLIFY(pcst)
#ifdef   __MPI
      NULLIFY(pcsum)
#endif
      !-------------------------------------------------------------------------
      !  Return
      !-------------------------------------------------------------------------
      CALL substop('ppm_tree_boxcost',t0,info)
      RETURN
      CONTAINS
      SUBROUTINE check
         DO ip=1,nbox
            DO i=1,ppm_dim
               IF (min_box(i,ip) .GT. max_box(i,ip)) THEN
                  info = ppm_error_error
                  CALL ppm_error(ppm_err_argument,'ppm_tree_boxcost',     &
     &                'min_box must be <= max_box !',__LINE__,info)
                  GOTO 8888
               ENDIF
            ENDDO
         ENDDO
 8888    CONTINUE
      END SUBROUTINE check
#if   __KIND == __SINGLE_PRECISION
      END SUBROUTINE ppm_tree_boxcost_s
#elif __KIND == __DOUBLE_PRECISION
      END SUBROUTINE ppm_tree_boxcost_d
#endif
