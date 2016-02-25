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
      SUBROUTINE ppm_tree_boxcost_s(Nm,weights,min_box,max_box, &
      &          nbox,lhbx,lpdx,boxcost,info,pcost)
#elif __KIND == __DOUBLE_PRECISION
      SUBROUTINE ppm_tree_boxcost_d(Nm,weights,min_box,max_box, &
      &          nbox,lhbx,lpdx,boxcost,info,pcost)
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
      USE ppm_module_mpi
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
      INTEGER,               DIMENSION(:,:),         INTENT(IN   ) :: Nm
      !!! Number of grid points in each box.
      REAL(MK),              DIMENSION(3),           INTENT(IN   ) :: weights
      !!! Weights for the three cost contributions: particles, mesh,
      !!! geometry for defining the total cost of a box.
      REAL(ppm_kind_double), DIMENSION(:,:),         INTENT(IN   ) :: min_box
      !!! The minimum coordinate of the boxes.
      !!!
      !!! 1st index: x,y[,z]                                                   +
      !!! 2nd: box ID
      REAL(ppm_kind_double), DIMENSION(:,:),         INTENT(IN   ) :: max_box
      !!! The maximum coordinate of the boxes.
      !!!
      !!! 1st index: x,y[,z]                                                   +
      !!! 2nd: box ID
      INTEGER,                                       INTENT(IN   ) :: nbox
      !!! Number of boxes to compute cost for.
      INTEGER,               DIMENSION(:),           INTENT(IN   ) :: lhbx
      !!! Pointer to first particle in each box to determine cost for.
      !!! This only needs to be present if there are any particles.
      INTEGER,               DIMENSION(:),           INTENT(IN   ) :: lpdx
      !!! Index list of points in each box to determine cost for.
      !!! This only needs to be present if there are any particles.
      REAL(ppm_kind_double), DIMENSION(:),           POINTER       :: boxcost
      !!! Costs of the boxes 1..nbox. boxcost(i) is the cost of box  boxes(i).
      INTEGER,                                       INTENT(  OUT) :: info
      !!! Return status, 0 on success
      REAL(MK),              DIMENSION(:), OPTIONAL, INTENT(IN   ) :: pcost
      !!! Argument of length Npart, specifying the
      !!! computational cost of each particle.
      !-------------------------------------------------------------------------
      !  Local variables
      !-------------------------------------------------------------------------
      REAL(ppm_kind_double)                     :: t0
      REAL(ppm_kind_double), DIMENSION(ppm_dim) :: len_box
      REAL(ppm_kind_double)                     :: meshtotal,geomtotal,dm

      INTEGER :: i,ip,j

      CHARACTER(LEN=ppm_char) :: caller="ppm_tree_boxcost"
      !-------------------------------------------------------------------------
      !  Externals
      !-------------------------------------------------------------------------

      !-------------------------------------------------------------------------
      !  Initialize
      !-------------------------------------------------------------------------
      CALL substart(caller,t0,info)

      !-------------------------------------------------------------------------
      !  Check input arguments
      !-------------------------------------------------------------------------
      IF (ppm_debug.GT.0) THEN
         CALL check
         IF (info .NE. 0) GOTO 9999
      ENDIF

      !-------------------------------------------------------------------------
      !  If we have less than 1 box to compute, exit
      !-------------------------------------------------------------------------
      IF (nbox.LT.1) THEN
         IF (ppm_debug.GT.0) THEN
            stdout("No boxes to be computed. Exiting.")
         ENDIF
         GOTO 9999
      ENDIF

      !-------------------------------------------------------------------------
      !  Computation of particle contributions to box costs
      !-------------------------------------------------------------------------
      IF (have_particles.AND.ABS(weights(1)).GT.0.0_MK) THEN
         pcst =0.0_ppm_kind_double

         IF (PRESENT(pcost)) THEN
            DO i=1,nbox
               DO j=lhbx(i),lhbx(i+1)-1
                  ip = lpdx(j)
                  dm = REAL(pcost(ip),ppm_kind_double)
                  pcst(i) = pcst(i) + dm
               ENDDO
            ENDDO
         ELSE
            DO i=1,nbox
               DO j=lhbx(i),lhbx(i+1)-1
                  pcst(i) = pcst(i) + 1.0_ppm_kind_double
               ENDDO
            ENDDO
         ENDIF

#ifdef __MPI
         pcsum=0.0_ppm_kind_double
         !---------------------------------------------------------------------
         !  Allreduce of all particle sums
         !---------------------------------------------------------------------
         CALL MPI_Allreduce(pcst,pcsum,nbox,MPI_DOUBLE_PRECISION,MPI_SUM,ppm_comm,info)
         or_fail_MPI("MPI_Allreduce of particle costs",ppm_error=ppm_error_fatal)

         pcst = pcsum
#endif
      ENDIF  ! have_particles

      !-------------------------------------------------------------------------
      !  Add up all cost contributions to total cost.
      !-------------------------------------------------------------------------
      IF (ppm_dim.EQ.2) THEN
         DO i=1,nbox
            !-----------------------------------------------------------------
            !  Mesh and geometry based costs
            !-----------------------------------------------------------------
            meshtotal = 0.0_ppm_kind_double
            IF (have_mesh.AND.ABS(weights(2)).GT.0.0_MK) THEN
               meshtotal = REAL(Nm(1,i),ppm_kind_double)*REAL(Nm(2,i),ppm_kind_double)
            ENDIF

            len_box(1) = max_box(1,i)-min_box(1,i)
            len_box(2) = max_box(2,i)-min_box(2,i)
            geomtotal  = len_box(1)*len_box(2)

            !-----------------------------------------------------------------
            !  Compute weighted sum of costs and store it
            !-----------------------------------------------------------------
            IF (have_particles) THEN
               boxcost(i) =   pcst(i)*REAL(weights(1),ppm_kind_double) + &
               &            meshtotal*REAL(weights(2),ppm_kind_double) + &
               &            geomtotal*REAL(weights(3),ppm_kind_double)
            ELSE
               boxcost(i) = meshtotal*REAL(weights(2),ppm_kind_double) + &
               &            geomtotal*REAL(weights(3),ppm_kind_double)
            ENDIF
         ENDDO
      ELSE
         DO i=1,nbox
            !-----------------------------------------------------------------
            !  Mesh and geometry based costs
            !-----------------------------------------------------------------
            meshtotal = 0.0_ppm_kind_double
            IF (have_mesh .AND. ABS(weights(2)).GT.0.0_MK) THEN
               meshtotal = REAL(Nm(1,i),ppm_kind_double)*REAL(Nm(2,i),ppm_kind_double)*REAL(Nm(3,i),ppm_kind_double)
            ENDIF

            len_box(1) = max_box(1,i)-min_box(1,i)
            len_box(2) = max_box(2,i)-min_box(2,i)
            len_box(3) = max_box(3,i)-min_box(3,i)
            geomtotal  = len_box(1)*len_box(2)*len_box(3)

            !-----------------------------------------------------------------
            !  Compute weighted sum of costs and store it
            !-----------------------------------------------------------------
            IF (have_particles) THEN
               boxcost(i) =   pcst(i)*REAL(weights(1),ppm_kind_double) + &
               &            meshtotal*REAL(weights(2),ppm_kind_double) + &
               &            geomtotal*REAL(weights(3),ppm_kind_double)
            ELSE
               boxcost(i) = meshtotal*REAL(weights(2),ppm_kind_double) + &
               &            geomtotal*REAL(weights(3),ppm_kind_double)
            ENDIF
         ENDDO
      ENDIF

      !-------------------------------------------------------------------------
      !  Free local memory
      !-------------------------------------------------------------------------
      9999 CONTINUE
      !-------------------------------------------------------------------------
      !  Return
      !-------------------------------------------------------------------------
      CALL substop(caller,t0,info)
      RETURN
      CONTAINS
      SUBROUTINE check
         DO ip=1,nbox
            DO i=1,ppm_dim
               IF (min_box(i,ip).GT.max_box(i,ip)) THEN
                  fail("min_box must be <= max_box !",exit_point=8888)
               ENDIF
            ENDDO
         ENDDO
      8888 CONTINUE
      END SUBROUTINE check
#if   __KIND == __SINGLE_PRECISION
      END SUBROUTINE ppm_tree_boxcost_s
#elif __KIND == __DOUBLE_PRECISION
      END SUBROUTINE ppm_tree_boxcost_d
#endif
