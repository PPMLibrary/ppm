      !-------------------------------------------------------------------------
      !  Subroutine   :                ppm_tree_cutpos.f
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

#if    __KIND == __SINGLE_PRECISION
       SUBROUTINE ppm_tree_cutpos_s(xp,Npart,weights,min_box,max_box,   &
      &    cutbox,ncut,minboxsize,icut,neigh_constraints,num_constr,cpos,info,pcost)
#elif __KIND == __DOUBLE_PRECISION
       SUBROUTINE ppm_tree_cutpos_d(xp,Npart,weights,min_box,max_box,   &
      &    cutbox,ncut,minboxsize,icut,neigh_constraints,num_constr,cpos,info,pcost)
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
      IMPLICIT NONE
#if   __KIND == __SINGLE_PRECISION
      INTEGER, PARAMETER :: MK = ppm_kind_single
#elif __KIND == __DOUBLE_PRECISION
      INTEGER, PARAMETER :: MK = ppm_kind_double
#endif
      !------------------------------------------------------------------------
      !  Includes
      !------------------------------------------------------------------------
#include "ppm_define.h"
#ifdef __MPI
      INCLUDE 'mpif.h'
#endif
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
      REAL(MK), DIMENSION(:,:,:), INTENT(IN   ) :: neigh_constraints
      ! access: dimension, constraintid, 1 (from) 2 (to)
      INTEGER,  DIMENSION(:),   INTENT(IN )   :: num_constr
      ! number of constraints
      REAL(MK), DIMENSION(:  ), POINTER       :: cpos
      !!! Positions of best cuts. index: 1..ncut.
      INTEGER                 , INTENT(  OUT) :: info
      !!! Return status, 0 on success
      !------------------------------------------------------------------------
      ! Local variables
      !------------------------------------------------------------------------
      REAL(MK), DIMENSION(ppm_dim)            :: len_box
      REAL(MK), DIMENSION(ncut+1)             :: pc,pcsum
      INTEGER , DIMENSION(2)                  :: ldc
      REAL(MK)                                :: t0,dm,meshtotal,geomtotal,pl,pr
      REAL(MK)                                :: pmass,mmass,gmass,tmass
      REAL(MK)                                :: partpos,midpos,lmyeps
      INTEGER                                 :: i,j,ip,cutdir,ncp1,iopt,k,kt,temp_r,ii

#ifdef __MPI
      INTEGER                                 :: MPTYPE
#endif
#if   __KIND == __SINGLE_PRECISION
      lmyeps = ppm_myepss
#elif __KIND == __DOUBLE_PRECISION
      lmyeps = ppm_myepsd
#endif
      !------------------------------------------------------------------------
      ! Externals
      !------------------------------------------------------------------------

      !------------------------------------------------------------------------
      ! Initialise
      !------------------------------------------------------------------------
      CALL substart('ppm_tree_cutpos',t0,info)
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
              CALL ppm_write(ppm_rank,'ppm_tree_cutpos',   &
                   'No cut directions present. Exiting. ',info)
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
          IF (info .NE. 0) THEN
              info = ppm_error_fatal
              CALL ppm_error(ppm_err_mpi_fail,'ppm_tree_cutpos',   &
                   'MPI_Allreduce of projected particles',__LINE__,info)
              GOTO 9999
          ENDIF
          pc = pcsum
#endif
      ENDIF !have_particles

!       IF(ppm_rank .EQ. 0)THEN
!             print *, 'hello asd'
!          ENDIF

      !------------------------------------------------------------------------
      ! Total weight of mesh and geometry
      ! Convert to REAL first to avoid integer overflow.
      !------------------------------------------------------------------------
      meshtotal = 0.0_MK
      IF (ppm_dim .EQ. 2) THEN
         IF (have_mesh .AND. weights(2) .NE. 0) THEN
             meshtotal = REAL(Nm_box(1,cutbox),MK)*REAL(Nm_box(2,cutbox),MK) 
         ENDIF
         geomtotal = len_box(1)*len_box(2)
      ELSE
         IF (have_mesh .AND. weights(2) .NE. 0) THEN
             meshtotal = REAL(Nm_box(1,cutbox),MK)*   &
     &           REAL(Nm_box(2,cutbox),MK)*REAL(Nm_box(3,cutbox),MK)
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
          info = ppm_error_fatal
          CALL ppm_error(ppm_err_argument,'ppm_tree_cutpos',   &
               'Total cost is 0! Are all weights 0?',__LINE__,info)
          GOTO 9999
      ENDIF

!          IF(ppm_rank .EQ. 0)THEN
!             print *, 'hello cuttyy'
!          ENDIF
      !------------------------------------------------------------------------
      ! The optimal cut position is in the weighted center of mass
      !------------------------------------------------------------------------
      DO i=1,ncut
          cutdir = icut(i)         
          IF (have_particles .AND. weights(1) .NE. 0.0_MK) THEN
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
!<<<< haeckic begin >>>>!
          !--------------------------------------------------------------------
          !Enforce that neighboring ghostlayers are respected
          !--------------------------------------------------------------------
          ! 2. if we are in a uncutable, take the closest
          ! try to find a non cutable region
          ! if not found -> ok
          ! else travel to the left and right and take smaller

!           IF (ppm_rank .eq. 0) THEN
!             print *, 'hello1'
!           ENDIF

          ! Find a range we are in
          IF (num_constr(cutdir) .GT. 0) THEN
            k = 1
            temp_r = k
            DO WHILE (.NOT.((cpos(i)-lmyeps).GT.neigh_constraints(cutdir,k,1) &
      &                         .AND. (cpos(i)+lmyeps).LT.neigh_constraints(cutdir,k,2)))
               IF (neigh_constraints(cutdir,k,2) .GT. neigh_constraints(cutdir,temp_r,2)) THEN
                  temp_r = k
               ENDIF
               k = k+1
               IF (k .GT. num_constr(cutdir)) THEN
                  EXIT
               ENDIF
            ENDDO
!          IF (ppm_rank .eq. 0 .and. cutbox .eq. 137) THEN
!             print *, 'hello2', k, temp_r
!           ENDIF
            ! If not found ->ok
            IF (k .LE. num_constr(cutdir)) THEN
               ! If found, search for two nearest possible, outside -> -1.0

               IF (neigh_constraints(cutdir,k,2) .GT. neigh_constraints(cutdir,temp_r,2)) THEN
                  temp_r = k
               ENDIF

               ! LEFT
               kt = k
               IF (kt .GT. 1) THEN
  912             CONTINUE
                  DO WHILE(neigh_constraints(cutdir,kt,1)+lmyeps .LT. neigh_constraints(cutdir,kt-1,2))
                     kt = kt-1
                     IF (kt .EQ. 1) THEN
                        EXIT
                     ENDIF
                  ENDDO
                  ! We need to check if it is violated by following constraints
                  ! If yes we continue scanline
                  IF (kt .GT. 1) THEN
                     DO ii = 1,kt-1
                        IF (neigh_constraints(cutdir,kt-ii,2)-lmyeps .GT. neigh_constraints(cutdir,kt,1)) THEN
                           kt = kt-ii
                           IF (kt .EQ. 1) THEN
                              EXIT
                           ENDIF
                           GOTO 912
                        ENDIF
                     ENDDO
                  ENDIF

               ENDIF
               pl = neigh_constraints(cutdir,kt,1)
! IF (ppm_rank .eq. 0 .and. cutbox .eq. 137) THEN
!             print *, 'hello3', kt, pl
!           ENDIF
               ! RIGHT
               kt = k+1
               IF (kt .LE. num_constr(cutdir)) THEN
                  DO WHILE(neigh_constraints(cutdir,kt,1)+lmyeps .LT. neigh_constraints(cutdir,temp_r,2) &
      &                      .OR. neigh_constraints(cutdir,kt,2)+lmyeps .LT. min_box(cutdir,cutbox)+minboxsize(cutdir))
                     IF (neigh_constraints(cutdir,kt,2) .GT. neigh_constraints(cutdir,temp_r,2)) THEN
                        temp_r = kt
                     ENDIF
                     IF (kt .EQ. num_constr(cutdir)) THEN
                        EXIT
                     ENDIF
                     kt = kt+1
                  ENDDO
               ENDIF
               pr = neigh_constraints(cutdir,temp_r,2)

! IF (ppm_rank .eq. 0 .and. cutbox .eq. 137) THEN
!             print *, 'hello4', kt, temp_r, pr
!           ENDIF
               ! take smaller
               IF (pl+lmyeps .LT. min_box(cutdir,cutbox)+minboxsize(cutdir)) THEN
                  IF (pr-lmyeps .GT. max_box(cutdir,cutbox)-minboxsize(cutdir)) THEN
                     print *, 'ERRRRR pr', cutbox, cutdir, pl, pr, cpos(i), minboxsize(cutdir), &
      &                         min_box(cutdir,cutbox)+minboxsize(cutdir), max_box(cutdir,cutbox)-minboxsize(cutdir)
                  ENDIF
                  cpos(i) = pr
               ELSEIF (pr-lmyeps .GT. max_box(cutdir,cutbox)-minboxsize(cutdir)) THEN
                  IF (pl+lmyeps .LT. min_box(cutdir,cutbox)+minboxsize(cutdir)) THEN
                     print *, 'ERRRRRO pl', cutbox, cutdir, pl, pr, cpos(i), minboxsize(cutdir), &
      &                         min_box(cutdir,cutbox)+minboxsize(cutdir), max_box(cutdir,cutbox)-minboxsize(cutdir)
                  ENDIF
                  cpos(i) = pl
               ELSE
                  ! take closer
                  IF(cpos(i)-pl .GT. pr-cpos(i)) THEN
                     cpos(i) = pr
                  ELSE
                     cpos(i) = pl
                  ENDIF
               ENDIF

            ENDIF

          ENDIF

!<<<< haeckic end >>>>!
      ENDDO

!       IF(ppm_rank .EQ. 0)THEN
!             print *, 'hello end'
!          ENDIF

      !------------------------------------------------------------------------
      ! Return
      !------------------------------------------------------------------------
 9999 CONTINUE
      CALL substop('ppm_tree_cutpos',t0,info)
      RETURN
      CONTAINS
      SUBROUTINE check
         IF (cutbox .LE. 0) THEN
            info = ppm_error_error
            CALL ppm_error(ppm_err_argument,'ppm_tree_cutpos',      &
                 'cutbox must be > 0 !',__LINE__,info)
            GOTO 8888
         ENDIF
         IF (SIZE(min_box,2) .LT. cutbox) THEN
            info = ppm_error_error
            CALL ppm_error(ppm_err_argument,'ppm_tree_cutpos',      &
                 'size of min_box must be at least cutbox !',__LINE__,info)
            GOTO 8888
         ENDIF
         IF (SIZE(max_box,2) .LT. cutbox) THEN
             info = ppm_error_error
             CALL ppm_error(ppm_err_argument,'ppm_tree_cutpos',     &
                  'size of max_box must be at least cutbox !',__LINE__,info)
             GOTO 8888
         ENDIF
         DO i=1,ppm_dim
            IF (min_box(i,cutbox) .GT. max_box(i,cutbox)) THEN
               info = ppm_error_error
               CALL ppm_error(ppm_err_argument,'ppm_tree_cutpos',   &
                    'min_box must be <= max_box !',__LINE__,info)
               GOTO 8888
            ENDIF
         ENDDO
 8888    CONTINUE
      END SUBROUTINE check
#if   __KIND == __SINGLE_PRECISION
      END SUBROUTINE ppm_tree_cutpos_s
#elif __KIND == __DOUBLE_PRECISION
      END SUBROUTINE ppm_tree_cutpos_d
#endif 

