      !-------------------------------------------------------------------------
      !  Subroutine   :                ppm_tree_cutpos.f
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

#if    __KIND == __SINGLE_PRECISION
       SUBROUTINE ppm_tree_cutpos_s(xp,Npart,weights,min_box,max_box,   &
      &    cutbox,ncut,minboxsize,icut,cpos,info,pcost)
#elif __KIND == __DOUBLE_PRECISION
       SUBROUTINE ppm_tree_cutpos_d(xp,Npart,weights,min_box,max_box,   &
      &    cutbox,ncut,minboxsize,icut,cpos,info,pcost)
#endif
      !!! This routine finds the best cuting positions for the given
      !!! cut directions.
      !!!
      !!! [NOTE]
      !!! The cost of a particle is counted as pcost (if
      !!! given) or 1. For meshes, the cost is 1 per mesh
      !!! point. For the geometry part, the cost of a
      !!! box is given by its volume. Use weights(3) to adjust this if needed.
      !------------------------------------------------------------------------
      !   Modules
      !------------------------------------------------------------------------
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
      !------------------------------------------------------------------------
      !  Includes
      !------------------------------------------------------------------------
      !------------------------------------------------------------------------
      ! Arguments
      !------------------------------------------------------------------------
      REAL(MK), DIMENSION(:,:), INTENT(IN   ) :: xp
      !!! Position of particles
      REAL(MK), DIMENSION(:,:), INTENT(IN   ) :: min_box
      !!! Minimum coordinate of the boxes
      REAL(MK), DIMENSION(:,:), INTENT(IN   ) :: max_box
      !!! Maximum coordinate of the boxes
      REAL(MK), DIMENSION(3  ), INTENT(IN   ) :: minboxsize
      !!! Minimum box size required in all spatial directions
      REAL(MK), DIMENSION(3  ), INTENT(IN   ) :: weights
      !!! Weights for the tree cost contributions: particles, mesh,
      !!! geometry for finding the cut
      REAL(MK), DIMENSION(:  ), OPTIONAL, INTENT(IN) :: pcost
      !!! Argument of length Npart, specifying the
      !!! computational cost of each particle.
      INTEGER                 , INTENT(IN   ) :: Npart
      !!! Number of particles
      INTEGER                 , INTENT(IN   ) :: ncut
      !!! Number of cut directions
      INTEGER                 , INTENT(IN   ) :: cutbox
      !!! ID of box to be cut
      INTEGER , DIMENSION(:  ), INTENT(IN   ) :: icut
      !!! Cut directions
      REAL(MK), DIMENSION(:  ), POINTER       :: cpos
      !!! Positions of best cuts. index: 1..ncut.
      INTEGER                 , INTENT(  OUT) :: info
      !!! Return status, 0 on success
      !------------------------------------------------------------------------
      ! Local variables
      !------------------------------------------------------------------------
      REAL(MK), DIMENSION(ppm_dim) :: len_box
      REAL(MK), DIMENSION(ncut+1)  :: pc,pcsum
      REAL(ppm_kind_double)        :: t0
      REAL(MK)                     :: dm,meshtotal,geomtotal
      REAL(MK)                     :: pmass,mmass,gmass,tmass
      REAL(MK)                     :: partpos,midpos

      INTEGER , DIMENSION(2) :: ldc
      INTEGER                :: i,j,ip,cutdir,ncp1,iopt
#ifdef __MPI
      INTEGER                :: MPTYPE
#endif

      CHARACTER(LEN=ppm_char) :: caller="ppm_tree_cutpos"
      !------------------------------------------------------------------------
      ! Externals
      !------------------------------------------------------------------------

      !------------------------------------------------------------------------
      ! Initialize
      !------------------------------------------------------------------------
      CALL substart(caller,t0,info)
      !------------------------------------------------------------------------
      ! Check input arguments
      !------------------------------------------------------------------------
      IF (ppm_debug .GT. 0) THEN
         CALL check
         IF (info .NE. 0) GOTO 9999
      ENDIF

      !------------------------------------------------------------------------
      ! If we have less than 1 direction to cut, we are done
      !------------------------------------------------------------------------
      IF (ncut .LT. 1) THEN
         IF (ppm_debug .GT. 0) THEN
            stdout("No cut directions present. Exiting.")
         ENDIF
         GOTO 9999
      ENDIF

      !------------------------------------------------------------------------
      ! Compute the extension of the box
      !------------------------------------------------------------------------
      IF (ppm_dim .GT. 2) THEN
         len_box(1) = max_box(1,cutbox)-min_box(1,cutbox)
         len_box(2) = max_box(2,cutbox)-min_box(2,cutbox)
         len_box(3) = max_box(3,cutbox)-min_box(3,cutbox)
      ELSE
         len_box(1) = max_box(1,cutbox)-min_box(1,cutbox)
         len_box(2) = max_box(2,cutbox)-min_box(2,cutbox)
      ENDIF

#ifdef __MPI
      !------------------------------------------------------------------------
      ! Determine MPI data type
      !------------------------------------------------------------------------
#if   __KIND == __SINGLE_PRECISION
      MPTYPE = MPI_REAL
#elif __KIND == __DOUBLE_PRECISION
      MPTYPE = MPI_DOUBLE_PRECISION
#endif
#endif

      !------------------------------------------------------------------------
      ! Compute particle cost center of gravity in all directions
      !------------------------------------------------------------------------
      pc = 0.0_MK
      ncp1 = ncut+1
      IF (have_particles .AND. weights(1) .NE. 0) THEN
#ifdef __VECTOR
         DO i=1,ncut
            cutdir = icut(i)
            DO j=tree_lhbx(1,cutbox),tree_lhbx(2,cutbox)
               ip = tree_lpdx(j)
               dm = 1.0_MK
               IF (PRESENT(pcost)) dm = pcost(ip)
               pc(i) = pc(i) + (xp(cutdir,ip)*dm)
               pc(ncp1) = pc(ncp1) + dm
            ENDDO
         ENDDO
         ! replace this by something more clever in the future. Try counting
         ! in a way that avoids the division here
         pc(ncp1) = pc(ncp1)/REAL(ncut,MK)
#else
         DO j=tree_lhbx(1,cutbox),tree_lhbx(2,cutbox)
            ip = tree_lpdx(j)
            dm = 1.0_MK
            IF (PRESENT(pcost)) dm = pcost(ip)
            DO i=1,ncut
               cutdir = icut(i)
               pc(i) = pc(i) + (xp(cutdir,ip)*dm)
            ENDDO
            pc(ncp1) = pc(ncp1) + dm
         ENDDO
#endif

#ifdef __MPI
          !---------------------------------------------------------------------
          ! Allreduce of particles sums
          !---------------------------------------------------------------------
          CALL MPI_Allreduce(pc,pcsum,ncp1,MPTYPE,MPI_SUM,ppm_comm,info)
          or_fail_MPI("MPI_Allreduce of projected particles",ppm_error=ppm_error_fatal)

          pc = pcsum
#endif
      ENDIF !have_particles

      !------------------------------------------------------------------------
      ! Total weight of mesh and geometry
      ! Convert to REAL first to avoid integer overflow.
      !------------------------------------------------------------------------
      meshtotal = 0.0_MK
      IF (ppm_dim .EQ. 2) THEN
         IF (have_mesh .AND. ABS(weights(2)).GT.0.0_MK) THEN
            meshtotal = REAL(Nm_box(1,cutbox),MK)*REAL(Nm_box(2,cutbox),MK)
         ENDIF
         geomtotal = len_box(1)*len_box(2)
      ELSE
         IF (have_mesh .AND. ABS(weights(2)).GT.0.0_MK) THEN
            meshtotal = REAL(Nm_box(1,cutbox),MK)* &
            &           REAL(Nm_box(2,cutbox),MK)* &
            &           REAL(Nm_box(3,cutbox),MK)
         ENDIF
         geomtotal = len_box(1)*len_box(2)*len_box(3)
      ENDIF

      !------------------------------------------------------------------------
      ! Compute weighted masses of particles, mesh and geometry
      !------------------------------------------------------------------------
      pmass = pc(ncp1)*weights(1)
      mmass = meshtotal*weights(2)
      gmass = geomtotal*weights(3)
      tmass = pmass+mmass+gmass
      IF (tmass .EQ. 0) THEN
         fail("Total cost is 0! Are all weights 0?",ppm_error=ppm_error_fatal)
      ENDIF

      !------------------------------------------------------------------------
      ! The optimal cut position is in the weighted center of mass
      !------------------------------------------------------------------------
      DO i=1,ncut
         cutdir = icut(i)
         IF (have_particles .AND. ABS(weights(1)).GT.0.0_MK) THEN
            partpos = pc(i)/pc(ncp1)
         ELSE
            partpos = 0.0_MK
         ENDIF
         midpos  = min_box(cutdir,cutbox) + (0.5_MK*len_box(cutdir))
         cpos(i) = partpos*pmass + midpos*(mmass+gmass)
         cpos(i) = cpos(i)/tmass

         !--------------------------------------------------------------------
         ! Enforce that minboxsize is respected.
         !--------------------------------------------------------------------
         IF (cpos(i)-min_box(cutdir,cutbox) .LT. minboxsize(cutdir)) THEN
            cpos(i) = min_box(cutdir,cutbox)+minboxsize(cutdir)
         ENDIF
         IF (max_box(cutdir,cutbox)-cpos(i) .LT. minboxsize(cutdir)) THEN
            cpos(i) = max_box(cutdir,cutbox)-minboxsize(cutdir)
         ENDIF
      ENDDO

      !------------------------------------------------------------------------
      ! Return
      !------------------------------------------------------------------------
      9999 CONTINUE
      CALL substop(caller,t0,info)
      RETURN
      CONTAINS
      SUBROUTINE check
         IF (cutbox .LE. 0) THEN
            fail("cutbox must be > 0 !",exit_point=8888)
         ENDIF
         IF (SIZE(min_box,2) .LT. cutbox) THEN
            fail("size of min_box must be at least cutbox !",exit_point=8888)
         ENDIF
         IF (SIZE(max_box,2) .LT. cutbox) THEN
            fail("size of max_box must be at least cutbox !",exit_point=8888)
         ENDIF
         DO i=1,ppm_dim
            IF (min_box(i,cutbox) .GT. max_box(i,cutbox)) THEN
               fail("min_box must be <= max_box !",exit_point=8888)
            ENDIF
         ENDDO
      8888 CONTINUE
      END SUBROUTINE check
#if   __KIND == __SINGLE_PRECISION
      END SUBROUTINE ppm_tree_cutpos_s
#elif __KIND == __DOUBLE_PRECISION
      END SUBROUTINE ppm_tree_cutpos_d
#endif

