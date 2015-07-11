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
       SUBROUTINE ppm_tree_cutpos_s(xp,Npart,weights,min_box,max_box, &
       &          cutbox,ncut,minboxsize,icut,cpos,info,pcost)
#elif __KIND == __DOUBLE_PRECISION
       SUBROUTINE ppm_tree_cutpos_d(xp,Npart,weights,min_box,max_box, &
       &          cutbox,ncut,minboxsize,icut,cpos,info,pcost)
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
      REAL(MK),              DIMENSION(:,:),         INTENT(IN   ) :: xp
      !!! Position of particles
      INTEGER,                                       INTENT(IN   ) :: Npart
      !!! Number of particles
      REAL(MK),              DIMENSION(3),           INTENT(IN   ) :: weights
      !!! Weights for the tree cost contributions: particles, mesh,
      !!! geometry for finding the cut
      REAL(ppm_kind_double), DIMENSION(:,:),         INTENT(IN   ) :: min_box
      !!! Minimum coordinate of the boxes
      REAL(ppm_kind_double), DIMENSION(:,:),         INTENT(IN   ) :: max_box
      !!! Maximum coordinate of the boxes
      INTEGER,                                       INTENT(IN   ) :: cutbox
      !!! ID of box to be cut
      INTEGER,                                       INTENT(IN   ) :: ncut
      !!! Number of cut directions
      REAL(MK),              DIMENSION(:),           INTENT(IN   ) :: minboxsize
      !!! Minimum box size required in all spatial directions
      INTEGER,               DIMENSION(:),           INTENT(IN   ) :: icut
      !!! Cut directions
      REAL(ppm_kind_double), DIMENSION(:),           POINTER       :: cpos
      !!! Positions of best cuts. index: 1..ncut.
      INTEGER,                                       INTENT(  OUT) :: info
      !!! Return status, 0 on success
      REAL(MK),              DIMENSION(:), OPTIONAL, INTENT(IN   ) :: pcost
      !!! Argument of length Npart, specifying the
      !!! computational cost of each particle.
      !------------------------------------------------------------------------
      ! Local variables
      !------------------------------------------------------------------------
      REAL(ppm_kind_double)                     :: t0
      REAL(ppm_kind_double), DIMENSION(ncut+1)  :: pc,pc_sum
      REAL(ppm_kind_double)                     :: dm
      REAL(ppm_kind_double)                     :: meshtotal,geomtotal
      REAL(ppm_kind_double)                     :: pmass,mmass,gmass,tmass
      REAL(ppm_kind_double)                     :: partpos,midpos
      REAL(ppm_kind_double), DIMENSION(ppm_dim) :: len_box

      INTEGER :: i,j,ip,cutdir,ncp1
#ifdef __MPI3
      INTEGER :: request
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
      IF (ppm_debug.GT.0) THEN
         CALL check
         IF (info.NE.0) GOTO 9999
      ENDIF

      !------------------------------------------------------------------------
      ! If we have less than 1 direction to cut, we are done
      !------------------------------------------------------------------------
      IF (ncut.LT.1) THEN
         IF (ppm_debug.GT.0) THEN
            stdout("No cut directions present. Exiting.")
         ENDIF
         GOTO 9999
      ENDIF

      !------------------------------------------------------------------------
      ! Compute the extension of the box
      !------------------------------------------------------------------------
      IF (ppm_dim.GT.2) THEN
         len_box(1) = max_box(1,cutbox)-min_box(1,cutbox)
         len_box(2) = max_box(2,cutbox)-min_box(2,cutbox)
         len_box(3) = max_box(3,cutbox)-min_box(3,cutbox)
      ELSE
         len_box(1) = max_box(1,cutbox)-min_box(1,cutbox)
         len_box(2) = max_box(2,cutbox)-min_box(2,cutbox)
      ENDIF

      !------------------------------------------------------------------------
      ! Compute particle cost center of gravity in all directions
      !------------------------------------------------------------------------
      pc   = 0.0_ppm_kind_double
      ncp1 = ncut+1
      IF (have_particles.AND.ABS(weights(1)).GT.0.0_MK) THEN
#ifdef __VECTOR
         IF (PRESENT(pcost)) THEN
            DO i=1,ncut
               cutdir = icut(i)
               DO j=tree_lhbx(1,cutbox),tree_lhbx(2,cutbox)
                  ip       = tree_lpdx(j)
                  dm       = REAL(pcost(ip),ppm_kind_double)
                  pc(i)    = pc(i)    + REAL(xp(cutdir,ip),ppm_kind_double)*dm
                  pc(ncp1) = pc(ncp1) + dm
               ENDDO
            ENDDO
         ELSE
            DO i=1,ncut
               cutdir = icut(i)
               DO j=tree_lhbx(1,cutbox),tree_lhbx(2,cutbox)
                  ip       = tree_lpdx(j)
                  pc(i)    = pc(i)    + REAL(xp(cutdir,ip),ppm_kind_double)
                  pc(ncp1) = pc(ncp1) + 1.0_ppm_kind_double
               ENDDO
            ENDDO
         ENDIF

         ! replace this by something more clever in the future. Try counting
         ! in a way that avoids the division here
         pc(ncp1) = pc(ncp1)/REAL(ncut,ppm_kind_double)
#else
         IF (PRESENT(pcost)) THEN
            DO j=tree_lhbx(1,cutbox),tree_lhbx(2,cutbox)
               ip = tree_lpdx(j)
               dm = REAL(pcost(ip),ppm_kind_double)
               DO i=1,ncut
                  cutdir = icut(i)
                  pc(i)  = pc(i) + REAL(xp(cutdir,ip),ppm_kind_double)*dm
               ENDDO
               pc(ncp1) = pc(ncp1) + dm
            ENDDO
         ELSE
            DO j=tree_lhbx(1,cutbox),tree_lhbx(2,cutbox)
               ip = tree_lpdx(j)
               DO i=1,ncut
                  cutdir = icut(i)
                  pc(i) = pc(i) + REAL(xp(cutdir,ip),ppm_kind_double)
               ENDDO
               pc(ncp1) = pc(ncp1) + 1.0_ppm_kind_double
            ENDDO
         ENDIF
#endif

#ifdef __MPI
         !---------------------------------------------------------------------
         !  Allreduce of particle costs
         !---------------------------------------------------------------------
         pc_sum=0.0_ppm_kind_double
#ifdef __MPI3
         CALL MPI_Iallreduce(pc,pc_sum,ncp1,MPI_DOUBLE_PRECISION,MPI_SUM,ppm_comm,request,info)
         or_fail_MPI("MPI_Iallreduce of inertia tensor failed!",ppm_error=ppm_error_fatal)
#else
         CALL MPI_Allreduce(pc,pc_sum,ncp1,MPI_DOUBLE_PRECISION,MPI_SUM,ppm_comm,info)
         or_fail_MPI("MPI_Allreduce of inertia tensor failed!",ppm_error=ppm_error_fatal)
#endif
#endif
      ENDIF !have_particles

      !------------------------------------------------------------------------
      ! Total weight of mesh and geometry
      ! Convert to REAL first to avoid integer overflow.
      !------------------------------------------------------------------------
      IF (ppm_dim .EQ. 2) THEN
         IF (have_mesh.AND.ABS(weights(2)).GT.0.0_MK) THEN
            meshtotal = REAL(Nm_box(1,cutbox),ppm_kind_double)*REAL(Nm_box(2,cutbox),ppm_kind_double)
         ELSE
            meshtotal = 0.0_ppm_kind_double
         ENDIF
         geomtotal = len_box(1)*len_box(2)
      ELSE
         IF (have_mesh.AND.ABS(weights(2)).GT.0.0_MK) THEN
            meshtotal = REAL(Nm_box(1,cutbox),ppm_kind_double)* &
            &           REAL(Nm_box(2,cutbox),ppm_kind_double)* &
            &           REAL(Nm_box(3,cutbox),ppm_kind_double)
         ELSE
            meshtotal = 0.0_ppm_kind_double
         ENDIF
         geomtotal = len_box(1)*len_box(2)*len_box(3)
      ENDIF

      !------------------------------------------------------------------------
      ! Compute weighted masses of mesh and geometry
      !------------------------------------------------------------------------
      mmass = meshtotal*REAL(weights(2),ppm_kind_double)
      gmass = geomtotal*REAL(weights(3),ppm_kind_double)

      IF (have_particles.AND.ABS(weights(1)).GT.0.0_MK) THEN
#ifdef __MPI
#ifdef __MPI3
         !Wait till we have pc_sum
         CALL MPI_Wait(request,MPI_STATUS_IGNORE,info)
#endif
         pc = pc_sum
#endif
      ENDIF !have_particles

      !------------------------------------------------------------------------
      ! Compute weighted masses of particles
      !------------------------------------------------------------------------
      pmass = pc(ncp1)*REAL(weights(1),ppm_kind_double)
      tmass = pmass+mmass+gmass

      IF (ABS(tmass).LE.0.0_ppm_kind_double) THEN
         fail("Total cost is 0! Are all weights 0?",ppm_error=ppm_error_fatal)
      ENDIF

      !------------------------------------------------------------------------
      ! The optimal cut position is in the weighted center of mass
      !------------------------------------------------------------------------
      DO i=1,ncut
         IF (have_particles.AND.ABS(weights(1)).GT.0.0_MK) THEN
            partpos = pc(i)/pc(ncp1)
         ELSE
            partpos = 0.0_ppm_kind_double
         ENDIF

         cutdir  = icut(i)

         midpos  = min_box(cutdir,cutbox) + 0.5_ppm_kind_double*len_box(cutdir)
         cpos(i) = partpos*pmass + midpos*(mmass+gmass)
         cpos(i) = cpos(i)/tmass

         !--------------------------------------------------------------------
         ! Enforce that minboxsize is respected.
         !--------------------------------------------------------------------
         IF (cpos(i)-min_box(cutdir,cutbox).LT.REAL(minboxsize(cutdir),ppm_kind_double)) THEN
            cpos(i) = min_box(cutdir,cutbox)+REAL(minboxsize(cutdir),ppm_kind_double)
         ENDIF
         IF (max_box(cutdir,cutbox)-cpos(i).LT.REAL(minboxsize(cutdir),ppm_kind_double)) THEN
            cpos(i) = max_box(cutdir,cutbox)-REAL(minboxsize(cutdir),ppm_kind_double)
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
         IF (SIZE(min_box,2).LT.cutbox) THEN
            fail("size of min_box must be at least cutbox !",exit_point=8888)
         ENDIF
         IF (SIZE(max_box,2).LT.cutbox) THEN
            fail("size of max_box must be at least cutbox !",exit_point=8888)
         ENDIF
         DO i=1,ppm_dim
            IF (min_box(i,cutbox).GT.max_box(i,cutbox)) THEN
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

