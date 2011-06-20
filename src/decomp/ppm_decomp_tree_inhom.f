      !<<<< haeckic begin >>>>!
      !-------------------------------------------------------------------------
      !  Subroutine   :                  ppm_decomp_tree
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

#if   __KIND == __SINGLE_PRECISION
      SUBROUTINE decomp_tree_inhom_s(xp,Npart,min_phys,max_phys, &
     &   ghost_req,tolerance,min_sub,max_sub,minboxsizes,bcdef,has_one_way,nsubs,info,pcost)
#elif __KIND == __DOUBLE_PRECISION
      SUBROUTINE decomp_tree_inhom_d(xp,Npart,min_phys,max_phys, &
     &   ghost_req,tolerance,min_sub,max_sub,minboxsizes,bcdef,has_one_way,nsubs,info,pcost)
#endif
      !!! Performs a tree-like decomposition.
      !!! It subdivides space until the number of leaves in the
      !!! tree exceeds the number of processors *and* until the
      !!! variance of the number of particles in the leaves is
      !!! below a certain tolerance.

      !-------------------------------------------------------------------------
      !  Modules 
      !-------------------------------------------------------------------------
      USE ppm_module_data
      USE ppm_module_substart
      USE ppm_module_substop
      USE ppm_module_error
      USE ppm_module_alloc
      USE ppm_module_write
      USE ppm_module_find_neigh
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
#ifdef __MPI
      INCLUDE 'mpif.h'
#endif
      !-------------------------------------------------------------------------
      !  Arguments     
      !-------------------------------------------------------------------------
      REAL(MK), DIMENSION(:,:), INTENT(IN   ) :: xp
      !!! Position of the particles
      INTEGER                 , INTENT(IN   ) :: Npart
      !!! Number of particles
      REAL(MK), DIMENSION(:)  , INTENT(IN   ) :: min_phys
      !!! Minimum coordinate of the physical/computational domain
      REAL(MK), DIMENSION(:)  , INTENT(IN   ) :: max_phys
      !!! Maximum coordinate of the physical/computational domain
      REAL(MK), DIMENSION(:)  , OPTIONAL, INTENT(IN) :: pcost
      !!! Argument of length Npart, specifying the
      !!! computational cost of each particle.
      REAL(MK), DIMENSION(:,:), INTENT(IN   ) :: ghost_req
      !!! Ghostsize requirements for all particles
      !!! 1st: x,y,(z)  2nd: Particle Id
      LOGICAL,                  INTENT(IN   ) ::  has_one_way
      !!! If this logical is set, the decomposition considers
      !!! one way interaction, i.e. if one particle overlaps the other
      !!! but not vice versa, now the minboxsize is also
      !!! depending on its neighbors ghostlayers
      INTEGER , DIMENSION(:  ), INTENT(IN   ) :: bcdef
      !!! Boundary conditions for the topology
      !!!
      !!! NOTE: first index is 1-6 (each of the faces)
      REAL(MK)                , INTENT(IN   ) :: tolerance
      !!! Maximum relative variance of the box cost allowed before
      !!! terminating the tree; tolerance >= 0
      REAL(MK), DIMENSION(:,:), POINTER       :: min_sub
      !!! Minimum extent of the subdomain
      REAL(MK), DIMENSION(:,:), POINTER       :: max_sub
      !!! Maximum extent of the subdomain
      REAL(MK), DIMENSION(:,:), POINTER       :: minboxsizes
      !!! Minimum boxsizes required in each dimension
      !!! 1st: x,y,(z)          2nd: boxId
      INTEGER                 , INTENT(  OUT) :: nsubs
      !!! Total number of subdomains
      INTEGER                 , INTENT(  OUT) :: info
      !!! Return status, 0 on success
      !-------------------------------------------------------------------------
      !  Local variables 
      !-------------------------------------------------------------------------
      REAL(MK), DIMENSION(ppm_dim)      :: len_phys, max_ghost
      REAL(MK), DIMENSION(2)            :: vector_in,vector_out
      REAL(MK)                          :: mean_npbx,var_npbx
      REAL(MK)                          :: t0, compare_tol
      REAL(MK), DIMENSION(:,:), POINTER :: min_box => NULL()
      REAL(MK), DIMENSION(:,:), POINTER :: max_box => NULL()
      REAL(MK), DIMENSION(:,:), POINTER :: min_box_temp => NULL()
      REAL(MK), DIMENSION(:,:), POINTER :: max_box_temp => NULL()
      REAL(MK), DIMENSION(:,:), POINTER :: minboxsizes_temp => NULL()
      REAL(MK), DIMENSION(:,:), POINTER :: max_ghost_temp => NULL()
      REAL(MK), DIMENSION(:,:), POINTER :: work    => NULL()
      INTEGER , DIMENSION(:), POINTER   :: numb_part => NULL()
      INTEGER , DIMENSION(:), POINTER   :: sub_index => NULL()
      INTEGER , DIMENSION(:), POINTER   :: npbx    => NULL()
      INTEGER , DIMENSION(:), POINTER   :: npbxg   => NULL()
      INTEGER , DIMENSION(:), POINTER   :: ppb     => NULL()
      INTEGER , DIMENSION(:  ), POINTER :: nneigh  => NULL()
      INTEGER , DIMENSION(:,:), POINTER :: ineigh  => NULL()
      INTEGER , DIMENSION(ppm_dim) :: ldc,ldd
      INTEGER :: ibox,fbox,jbox,kbox,lbox,nbox,ilevel,nlevel,isize
      INTEGER :: i,j,k,iopt,ipart,ii,jj,mem_req
      INTEGER :: istat,n1,n2,Npartg, nsubs_temp
      LOGICAL :: lcontinue, ghost_ok, ghost_ok_one, variance_not_ok
      CHARACTER(ppm_char) :: mesg
#ifdef __MPI
      INTEGER                                 :: MPTYPE
#endif
      !-------------------------------------------------------------------------
      !  Externals 
      !-------------------------------------------------------------------------

#ifdef __MPI
      !------------------------------------------------------------------------
      ! Determine MPI data type for max ghost layer reduction
      !------------------------------------------------------------------------
#if   __KIND == __SINGLE_PRECISION
    MPTYPE = MPI_REAL
#elif __KIND == __DOUBLE_PRECISION
    MPTYPE = MPI_DOUBLE_PRECISION
#endif
#endif

      !-------------------------------------------------------------------------
      !  Initialise 
      !-------------------------------------------------------------------------
      CALL substart('ppm_decomp_tree',t0,info)

      !-------------------------------------------------------------------------
      !  check input arguments
      !-------------------------------------------------------------------------
      IF (ppm_debug .GT. 0) THEN
        CALL check
        IF (info .NE. 0) GOTO 9999
      ENDIF 

      !-------------------------------------------------------------------------
      !  Allocate memory for a copy of the particles (the copy is needed since 
      !  the particles will be rearranged in the process of creating the tree)
      !-------------------------------------------------------------------------
      iopt   = ppm_param_alloc_fit
      ldc(1) = ppm_dim
      ldc(2) = Npart
      CALL ppm_alloc(work,ldc,iopt,info)
      IF (info.NE.0) THEN
         info = ppm_error_fatal
         CALL ppm_error(ppm_err_alloc,'ppm_decomp_tree',               &
     &                  'alloc of work array failed!',__LINE__,info)
         GOTO 9999
      ENDIF 

      !-------------------------------------------------------------------------
      !  Make a copy of the particles (since they will be rearranged in the
      !  process of creating the tree)
      !-------------------------------------------------------------------------
      IF (ppm_dim.EQ.2) THEN
         DO ipart=1,Npart
            work(1,ipart) = xp(1,ipart)
            work(2,ipart) = xp(2,ipart)
         ENDDO
      ELSE
         DO ipart=1,Npart
            work(1,ipart) = xp(1,ipart)
            work(2,ipart) = xp(2,ipart)
            work(3,ipart) = xp(3,ipart)
         ENDDO
      ENDIF 

#ifdef __MPI
      !-------------------------------------------------------------------------
      !  Count the total number of particles
      !-------------------------------------------------------------------------
      CALL MPI_AllReduce(Npart,Npartg,1,MPI_INTEGER,MPI_SUM,ppm_comm,info)
#else
      Npartg = Npart
#endif 

      !-------------------------------------------------------------------------
      !  Allocate some memory for the tree
      !-------------------------------------------------------------------------
      iopt   = ppm_param_alloc_grow
      ldc(1) = ppm_dim
      ldc(2) = 2000

      CALL ppm_alloc(min_box,ldc,iopt,info)
      IF (info.NE.0) THEN
         info = ppm_error_fatal
         CALL ppm_error(ppm_err_alloc,'ppm_decomp_tree',               &
     &                  'alloc of min_box failed!',__LINE__,info)
         GOTO 9999
      ENDIF 

      CALL ppm_alloc(max_box,ldc,iopt,info)
      IF (info.NE.0) THEN
         info = ppm_error_fatal
         CALL ppm_error(ppm_err_alloc,'ppm_decomp_tree',               &
     &                  'alloc of max_box failed!',__LINE__,info)
         GOTO 9999
      ENDIF 

      CALL ppm_alloc(min_box_temp,ldc,iopt,info)
      IF (info.NE.0) THEN
         info = ppm_error_fatal
         CALL ppm_error(ppm_err_alloc,'ppm_decomp_tree',               &
     &                  'alloc of min_box_temp failed!',__LINE__,info)
         GOTO 9999
      ENDIF 

      CALL ppm_alloc(max_box_temp,ldc,iopt,info)
      IF (info.NE.0) THEN
         info = ppm_error_fatal
         CALL ppm_error(ppm_err_alloc,'ppm_decomp_tree',               &
     &                  'alloc of max_box_temp failed!',__LINE__,info)
         GOTO 9999
      ENDIF 

      CALL ppm_alloc(minboxsizes_temp,ldc,iopt,info)
      IF (info.NE.0) THEN
         info = ppm_error_fatal
         CALL ppm_error(ppm_err_alloc,'ppm_decomp_tree',               &
     &                  'alloc of minboxsizes_temp failed!',__LINE__,info)
         GOTO 9999
      ENDIF
      CALL ppm_alloc(max_ghost_temp,ldc,iopt,info)
      IF (info.NE.0) THEN
         info = ppm_error_fatal
         CALL ppm_error(ppm_err_alloc,'ppm_decomp_tree',               &
     &                  'alloc of max_ghost_temp failed!',__LINE__,info)
         GOTO 9999
      ENDIF 

      CALL ppm_alloc(min_sub,ldc,iopt,info)
      IF (info.NE.0) THEN
         info = ppm_error_fatal
         CALL ppm_error(ppm_err_alloc,'ppm_decomp_tree',               &
     &                  'alloc of min_sub failed!',__LINE__,info)
         GOTO 9999
      ENDIF 

      CALL ppm_alloc(max_sub,ldc,iopt,info)
      IF (info.NE.0) THEN
         info = ppm_error_fatal
         CALL ppm_error(ppm_err_alloc,'ppm_decomp_tree',               &
     &                  'alloc of max_sub failed!',__LINE__,info)
         GOTO 9999
      ENDIF 

      CALL ppm_alloc(minboxsizes,ldc,iopt,info)
      IF (info.NE.0) THEN
         info = ppm_error_fatal
         CALL ppm_error(ppm_err_alloc,'ppm_decomp_tree',               &
     &                  'alloc of minboxsizes failed!',__LINE__,info)
         GOTO 9999
      ENDIF

      ldc(1) = ldc(2)
      CALL ppm_alloc(numb_part,ldc,iopt,info)
      IF (info.NE.0) THEN
         info = ppm_error_fatal
         CALL ppm_error(ppm_err_alloc,'ppm_decomp_tree',               &
     &                  'alloc of numb_part failed!',__LINE__,info)
         GOTO 9999
      ENDIF 
      CALL ppm_alloc(ppb,ldc,iopt,info)
      IF (info.NE.0) THEN
         info = ppm_error_fatal
         CALL ppm_error(ppm_err_alloc,'ppm_decomp_tree',               &
     &                  'alloc of ppb failed!',__LINE__,info)
         GOTO 9999
      ENDIF 

      CALL ppm_alloc(npbx,ldc,iopt,info)
      IF (info.NE.0) THEN
         info = ppm_error_fatal
         CALL ppm_error(ppm_err_alloc,'ppm_decomp_tree',               &
     &                  'alloc of npbx failed!',__LINE__,info)
         GOTO 9999
      ENDIF 

!--------------------------------------------------
!       CALL ppm_alloc(cost ,ldc,iopt,info)
!       IF (info.NE.0) THEN
!          info = ppm_error_fatal
!          CALL ppm_error(ppm_err_alloc,'ppm_decomp_tree',               &
!      &                  'alloc of cost failed!',__LINE__,info)
!          GOTO 9999
!       ENDIF 
!--------------------------------------------------

      CALL ppm_alloc(npbxg,ldc,iopt,info)
      IF (info.NE.0) THEN
         info = ppm_error_fatal
         CALL ppm_error(ppm_err_alloc,'ppm_decomp_tree',               &
     &                  'alloc of npbxg failed!',__LINE__,info)
         GOTO 9999
      ENDIF 

      !-------------------------------------------------------------------------
      !  Initialise the first box
      !-------------------------------------------------------------------------
      fbox  = 1
      compare_tol = max_phys(1) - min_phys(1)
      DO k=1,ppm_dim
         min_box(k,fbox) = min_phys(k)
         max_box(k,fbox) = max_phys(k)
         len_phys(k)     = max_phys(k) - min_phys(k)

         IF (compare_tol > len_phys(k)/2.0_mk) THEN
            compare_tol = len_phys(k)/2.0_mk
         ENDIF

         ! Determine minboxsizes => maximum of ghost_req in each dimension
         ! Needs to collect max from all processors
         minboxsizes_temp(k,1) = 0.0_mk
         DO ii=1,Npart
            ! is partice inside box?
            IF (xp(k,ii) .GE. min_box(k,1) .AND. xp(k,ii) .LE. max_box(k,1)) THEN
               IF (ghost_req(k,ii) .GT. minboxsizes_temp(k,1)) THEN
                  minboxsizes_temp(k,1) = ghost_req(k,ii)
               ENDIF
            ENDIF
         ENDDO
      ENDDO
      CALL MPI_Allreduce(minboxsizes_temp(1,1),max_ghost_temp(1,1),ppm_dim,MPTYPE,MPI_MAX,ppm_comm,info)
      DO k=1,ppm_dim
         minboxsizes_temp(k,1) = max_ghost_temp(k,1)
      ENDDO

      ppb(fbox)   = 1
      npbx(fbox)  = Npart
      npbxg(fbox) = Npart ! adapt this reduce!

      !-------------------------------------------------------------------------
      !  Compute the maximum number of levels 
      !-------------------------------------------------------------------------
      ! not needed, we go until # particles = 1 or ghost req violoated
     
      !-------------------------------------------------------------------------
      !  Split the boxes, starting at level 0 (the full box)
      !  kbox denotes the box being divided
      !  nbox is the current total number of boxes
      !  fbox is the id of the first box on the current level
      !  lbox is the id of the last  box on the current level
      !-------------------------------------------------------------------------
      lcontinue = .TRUE.
      lbox      = fbox
      nbox      = lbox
      nsubs     = 0
      DO WHILE (lcontinue)
         !----------------------------------------------------------------------
         !  Split the boxes at this level
         !----------------------------------------------------------------------

         ! We need a temporary array holding all boxes, i.e combine min_box and min_sub..
         IF (has_one_way) THEN

            DO i=1,nsubs
               DO k=1,ppm_dim
                  min_box_temp(k,i) = min_sub(k,i)
                  max_box_temp(k,i) = max_sub(k,i)
                  minboxsizes_temp(k,i) = max_ghost_temp(k,i)
               ENDDO
            ENDDO

            DO j=fbox,lbox
               i = nsubs + (j - fbox + 1)
               DO k=1,ppm_dim 
                  min_box_temp(k,i) = min_box(k,j)
                  max_box_temp(k,i) = max_box(k,j)
                  minboxsizes_temp(k,i) = 0.0_mk

                  IF (compare_tol > (max_box(k,j) - min_box(k,j))/2.0_mk) THEN
                     compare_tol = (max_box(k,j) - min_box(k,j))/2.0_mk
                  ENDIF

                  DO ii=1,Npart
                     ! is partice inside box?
                     IF (xp(k,ii) .GE. min_box(k,j) .AND. xp(k,ii) .LE. max_box(k,j)) THEN
                        IF (ghost_req(k,ii) .GT. minboxsizes_temp(k,i)) THEN
                           minboxsizes_temp(k,i) = ghost_req(k,ii)
                        ENDIF
                     ENDIF
                  ENDDO
                  
               ENDDO
               ! Needs to collect max from all processors
               CALL MPI_Allreduce(minboxsizes_temp(1,i),max_ghost_temp(1,i),ppm_dim,MPTYPE,MPI_MAX,ppm_comm,info)
               DO k=1,ppm_dim
                  minboxsizes_temp(k,i) = max_ghost_temp(k,i)
               ENDDO

            ENDDO

            CALL ppm_find_neigh(min_phys,max_phys,bcdef, &
     &              min_box_temp,max_box_temp,nsubs+(lbox-fbox+1),nneigh,ineigh,info)
         ELSE

            DO j=fbox,lbox
               i = nsubs + (j - fbox + 1)
               ! Determine all max_ghost_sizes for the boxes in focus
               ! Get the minimum box sizes
               DO k=1,ppm_dim
                  ! find out which particles are in this box and get maximum ghost req
                  ! Needs to collect max from all processors
                  minboxsizes_temp(k,i) = 0.0_mk
                  DO ii=1,Npart
                     ! is partice inside box?
                     IF (xp(k,ii) .GE. min_box(k,j) .AND. xp(k,ii) .LE. max_box(k,j)) THEN

                        IF (ghost_req(k,ii) .GT. minboxsizes_temp(k,i)) THEN
                           minboxsizes_temp(k,i) = ghost_req(k,ii)
                        ENDIF
                     ENDIF
                  ENDDO

                   ! Needs to collect max from all processors
                  CALL MPI_Allreduce(minboxsizes_temp(k,i),max_ghost_temp(k,i),1,MPTYPE,MPI_MAX,ppm_comm,info)
                  minboxsizes_temp(k,i) = max_ghost_temp(k,i)

               ENDDO
  
            ENDDO

         ENDIF

         nsubs_temp = nsubs

         DO kbox=fbox,lbox
            !-------------------------------------------------------------------
            !  Split non-empty boxes into four/eight (nbox is incremented)
            !-------------------------------------------------------------------

            ! Get here the neighbors influence if has_one_way
            IF (has_one_way) THEN
               ii = nsubs_temp+(kbox-fbox+1)
               DO i=1,nneigh(ii)
                  jj = ineigh(i,ii)

                  ! Update the new one
                  DO k=1,ppm_dim
                     ! 1. Calculate influence in this dimension
                     IF ((min_box_temp(k,ii) - min_box_temp(k,jj)) .GT. compare_tol &
   &                          .AND. ABS(min_box_temp(k,ii) - max_box_temp(k,jj)) .LT. compare_tol) THEN
                        IF (minboxsizes_temp(k,ii) < max_ghost_temp(k,jj)) THEN
                           minboxsizes_temp(k,ii) = max_ghost_temp(k,jj)
                        ENDIF

                     ELSEIF ((max_box_temp(k,jj) - max_box_temp(k,ii)) .GT. compare_tol &
   &                          .AND. ABS(max_box_temp(k,ii) - min_box_temp(k,jj)) .LT. compare_tol) THEN
                        IF (minboxsizes_temp(k,ii) < max_ghost_temp(k,jj)) THEN
                           minboxsizes_temp(k,ii) = max_ghost_temp(k,jj)
                        ENDIF
                     ENDIF

                  ENDDO

               ENDDO

            ENDIF
            
            ! Determine if we can split
            ghost_ok = .TRUE.
            DO k=1,ppm_dim
               IF (minboxsizes_temp(k,nsubs_temp+(kbox-fbox+1)) .GT. 0.5*(max_box(k,kbox) - min_box(k,kbox))) THEN
                  ghost_ok = .FALSE.
               ENDIF
            ENDDO

            IF (npbxg(kbox).GT.0 .AND. ghost_ok) THEN
               CALL ppm_decomp_boxsplit(work,ppb,npbx,kbox,nbox, &
     &                                  min_box,max_box,info)
               IF (info.NE.0) GOTO 100
            ELSE
               !----------------------------------------------------------------
               !  if empty or ghostlayer not fulfilled save it as a subdomain 
               !  ie increment the sub counter and check for memory 
               !----------------------------------------------------------------
               nsubs = nsubs + 1
               IF (nsubs.GT.SIZE(min_sub,2)) THEN
                  ldc(1) = ppm_dim
                  ldc(2) = 2*nsubs
                  iopt   = ppm_param_alloc_grow_preserve
                  CALL ppm_alloc(min_sub,ldc,iopt,info)
                  IF (info.NE.0) GOTO 200
                  CALL ppm_alloc(max_sub,ldc,iopt,info)
                  IF (info.NE.0) GOTO 200
                  CALL ppm_alloc(minboxsizes,ldc,iopt,info)
                  IF (info.NE.0) GOTO 200
                  ldc(1) = ldc(2)
                  CALL ppm_alloc(numb_part,ldc,iopt,info)
                  IF (info.NE.0) GOTO 200
           !       CALL ppm_alloc(cost,ldc,iopt,info)
           !       IF (info.NE.0) GOTO 200
               ENDIF 

               !----------------------------------------------------------------
               !  save the subdomain
               !----------------------------------------------------------------
               DO k=1,ppm_dim
                  min_sub(k,nsubs) = min_box(k,kbox)
                  max_sub(k,nsubs) = max_box(k,kbox)
                  minboxsizes(k,nsubs) = minboxsizes_temp(k,nsubs+(kbox-fbox+1))
                  numb_part(nsubs) = npbxg(kbox)
               ENDDO

           !    cost(nsubs) = 0
            ENDIF

         ENDDO

         !----------------------------------------------------------------------
         !  Catch errors in ppm_decomp_boxsplit
         !----------------------------------------------------------------------
     100 CONTINUE
#ifdef __MPI
         IF (ppm_debug.GT.0) THEN
            !-------------------------------------------------------------------
            !  MPI AllReduce is expensive so we only do this in debugging mode !
            !-------------------------------------------------------------------
            CALL MPI_AllReduce(info,i,1,MPI_INTEGER,MPI_SUM,ppm_comm,info)
            info = i
         ENDIF
#endif
         IF (info.NE.0) GOTO 9999

         !----------------------------------------------------------------------
         !  Catch errors in ppm_alloc
         !----------------------------------------------------------------------
     200 CONTINUE
#ifdef __MPI
         IF (ppm_debug.GT.0) THEN
            !-------------------------------------------------------------------
            !  MPI AllReduce is expensive so we only do this in debugging mode !
            !-------------------------------------------------------------------
            CALL MPI_AllReduce(info,i,1,MPI_INTEGER,MPI_SUM,ppm_comm,info)
            info = i
         ENDIF
#endif
         !----------------------------------------------------------------------
         !  Write detailed error message for the allocation
         !----------------------------------------------------------------------
         IF (info.NE.0) THEN
            info = ppm_error_fatal
            CALL ppm_error(ppm_err_alloc,'ppm_decomp_tree',                  &
     &                     'reallocating the tree failed (1)!',__LINE__,info)
            GOTO 9999
         ENDIF 

         !----------------------------------------------------------------------
         !  Update the level counter and boxes ie. consider now the new boxes
         !----------------------------------------------------------------------
         fbox   = lbox + 1
         lbox   = nbox
#ifdef __MPI
         !----------------------------------------------------------------------
         !  Communicate the number of particles found in the boxes of the new
         !  boxes at the next level (this will be expensive for nproc >> 1)
         !----------------------------------------------------------------------
         CALL MPI_AllReduce(npbx(fbox),npbxg(fbox),lbox-fbox+1,MPI_INTEGER, &
     &                      MPI_SUM,ppm_comm,info)
#else
         !----------------------------------------------------------------------
         !  Copy the global npbx from the local (the same in serial, but to
         !  keep the rest of the code the same for both parallel/serial ...)
         !----------------------------------------------------------------------
         DO kbox=fbox,lbox
            npbxg(kbox) = npbx(kbox)
         ENDDO
#endif

         !----------------------------------------------------------------------
         !  Check if to continue to the next level 
         !  (A) if we reached the finest level we exit
         !----------------------------------------------------------------------
         IF (lbox-fbox .LT. 0) lcontinue = .FALSE.

         !----------------------------------------------------------------------
         !  Check if to continue to the next level ...
         !  (B) if the variance of the mean number of particles in the boxes is
         !  below the tolerance we exit
         !----------------------------------------------------------------------
         IF (nbox.GT.ppm_nproc) THEN
            !-------------------------------------------------------------------
            !  Compute the mean and variance of the particle in the boxes at
            !  the next level
            !-------------------------------------------------------------------
            var_npbx = 0.0_MK 
            mean_npbx = 0.0_MK 
            IF (PRESENT(pcost)) THEN 
               !----------------------------------------------------------------
               !  If the pcost is present take this into account
               !----------------------------------------------------------------
               DO kbox=fbox,lbox
                  DO k=ppb(kbox),ppb(kbox)+npbx(kbox)-1
                     mean_npbx = mean_npbx + pcost(k)    
                      var_npbx =  var_npbx + pcost(k)**2 
                  ENDDO
               ENDDO

#ifdef __MPI
               !----------------------------------------------------------------
               !  In parallel we need to accumulate the contribs. from all procs
               !----------------------------------------------------------------
               vector_in(1) = mean_npbx
               vector_in(2) =  var_npbx
#if __KIND == __SINGLE_PRECISION
               CALL MPI_AllReduce(vector_in,vector_out,2,MPI_REAL, &
     &                            MPI_SUM,ppm_comm,info)
#else
               CALL MPI_AllReduce(vector_in,vector_out,2,MPI_DOUBLE_PRECISION, &
     &                            MPI_SUM,ppm_comm,info)
#endif
               mean_npbx = vector_out(1) 
                var_npbx = vector_out(2)
#endif
            ELSE
               !----------------------------------------------------------------
               !  If the pcost is not present compute the cost as the sum of the
               !  (global) number of particles in the boxes
               !----------------------------------------------------------------
               DO kbox=fbox,lbox
                  mean_npbx = mean_npbx + REAL(npbxg(kbox),MK)
                   var_npbx =  var_npbx + REAL(npbxg(kbox),MK)**2
               ENDDO
            ENDIF
            !-------------------------------------------------------------------
            !  Compute the mean and the variance
            !-------------------------------------------------------------------
            mean_npbx = mean_npbx/REAL(lbox - fbox + 1,MK)
             var_npbx = var_npbx/REAL(lbox - fbox + 1,MK) - mean_npbx**2

            !-------------------------------------------------------------------
            !  If the variance is acceptable ie less than the tolerance times 
            !  the mean number of particles, then stop 
            !-------------------------------------------------------------------
            IF (var_npbx.LT.tolerance*mean_npbx) THEN
               lcontinue = .FALSE.
            ENDIF 
         ENDIF 

         !----------------------------------------------------------------------
         !  Compute the required amount of memory at the next level
         !----------------------------------------------------------------------
         mem_req = nbox + (lbox + 1 - fbox)*(2**ppm_dim)
         

         !----------------------------------------------------------------------
         !  Check if we have enough memory
         !----------------------------------------------------------------------
         IF (lcontinue.AND.SIZE(npbx).LT.mem_req) THEN 
            iopt   = ppm_param_alloc_fit_preserve
            ldc(1) = ppm_dim
            ldc(2) = mem_req
            CALL ppm_alloc(min_box,ldc,iopt,info)
            IF (info.NE.0) GOTO 300
            CALL ppm_alloc(max_box,ldc,iopt,info)
            IF (info.NE.0) GOTO 300
            CALL ppm_alloc(min_box_temp,ldc,iopt,info) 
            IF (info.NE.0) GOTO 300
            CALL ppm_alloc(max_box_temp,ldc,iopt,info)
            IF (info.NE.0) GOTO 300
            CALL ppm_alloc(minboxsizes_temp,ldc,iopt,info)
            IF (info.NE.0) GOTO 300
            CALL ppm_alloc(max_ghost_temp,ldc,iopt,info)
            IF (info.NE.0) GOTO 300
            ldc(1) = mem_req
            CALL ppm_alloc(ppb  ,ldc,iopt,info)
            IF (info.NE.0) GOTO 300
            CALL ppm_alloc(npbx ,ldc,iopt,info)
            IF (info.NE.0) GOTO 300
            CALL ppm_alloc(npbxg,ldc,iopt,info)
            IF (info.NE.0) GOTO 300

            !-------------------------------------------------------------------
            !  Catch errors in ppm_alloc in debugging mode
            !-------------------------------------------------------------------
     300    CONTINUE
#ifdef __MPI
            IF (ppm_debug.GT.0) THEN
               !----------------------------------------------------------------
               !  MPI AllReduce is expensive so we only do this in debugging 
               !----------------------------------------------------------------
               CALL MPI_AllReduce(info,i,1,MPI_INTEGER,MPI_SUM,ppm_comm,info)
               info = i
            ENDIF
#endif
            IF (info.NE.0) THEN
               info = ppm_error_fatal
               CALL ppm_error(ppm_err_alloc,'ppm_decomp_tree',                 &
     &                        'reallocating the tree failed (2)!',__LINE__,info)
               GOTO 9999
            ENDIF 
         ENDIF

 
      ENDDO ! end of main DO WHILE loop

      !-------------------------------------------------------------------------
      !  Now use the boxes as our domains (including empty boxes, since we dont
      !  want void space which would/could cause particles moving into empty 
      !  space ... well this could be ok: would need an immediate rebuild of
      !  the decomposition  - but then again - perhaps this would happen quite
      !  frequently.
      !-------------------------------------------------------------------------
      mem_req = nsubs + (lbox + 1 - fbox) 
      IF (mem_req.GT.SIZE(min_sub,2)) THEN
         iopt   = ppm_param_alloc_fit_preserve
         ldc(1) = ppm_dim
         ldc(2) = mem_req
         CALL ppm_alloc(min_sub,ldc,iopt,info)
         IF (info.NE.0) GOTO 400
         CALL ppm_alloc(max_sub,ldc,iopt,info)
         IF (info.NE.0) GOTO 400
         CALL ppm_alloc(minboxsizes,ldc,iopt,info)
         IF (info.NE.0) GOTO 400
         ldc(1) = mem_req
         CALL ppm_alloc(numb_part,ldc,iopt,info)
         IF (info.NE.0) GOTO 400
 400     CONTINUE
         IF (info.NE.0) THEN
            info = ppm_error_fatal
            CALL ppm_error(ppm_err_alloc,'ppm_decomp_tree',                  &
     &                     'reallocating the subs failed!',__LINE__,info)
            GOTO 9999
         ENDIF 
      ENDIF 

      !-------------------------------------------------------------------------
      !  Perform the copy
      !-------------------------------------------------------------------------
      IF     (ppm_dim.EQ.2) THEN
         !----------------------------------------------------------------------
         !  in 2D
         !----------------------------------------------------------------------
         DO ibox=fbox,lbox
            nsubs            = nsubs + 1
            min_sub(1,nsubs) = min_box(1,ibox)
            min_sub(2,nsubs) = min_box(2,ibox)

            max_sub(1,nsubs) = max_box(1,ibox)
            max_sub(2,nsubs) = max_box(2,ibox)

            numb_part(nsubs) = npbxg(ibox)

            DO k=1,ppm_dim
               ! find out which particles are in this box and get maximum ghost req
               ! Needs to collect max from all processors
               max_ghost(k) = 0.0_mk
               DO i=1,Npart
                  ! is partice inside box?
                  IF (xp(k,i) .GE. min_box(k,ibox) .AND. xp(k,i) .LE. max_box(k,ibox)) THEN
                     IF (ghost_req(k,i) .GT. max_ghost(k)) THEN
                        max_ghost(k) = ghost_req(k,i)
                     ENDIF
                  ENDIF
               ENDDO
            ENDDO
            ! Needs to collect max from all processors
            CALL MPI_Allreduce(max_ghost(1),minboxsizes(1,nsubs),ppm_dim,MPTYPE,MPI_MAX,ppm_comm,info)

         ENDDO
      ELSEIF (ppm_dim.EQ.3) THEN
         !----------------------------------------------------------------------
         !  in 3D
         !----------------------------------------------------------------------
         DO ibox=fbox,lbox
            nsubs            = nsubs + 1
            min_sub(1,nsubs) = min_box(1,ibox)
            min_sub(2,nsubs) = min_box(2,ibox)
            min_sub(3,nsubs) = min_box(3,ibox)

            max_sub(1,nsubs) = max_box(1,ibox)
            max_sub(2,nsubs) = max_box(2,ibox)
            max_sub(3,nsubs) = max_box(3,ibox)

            numb_part(nsubs) = npbxg(ibox)

            DO k=1,ppm_dim
               ! find out which particles are in this box and get maximum ghost req
               ! Needs to collect max from all processors
               max_ghost(k) = 0.0_mk
               DO i=1,Npart
                  ! is partice inside box?
                  IF (xp(k,i) .GE. min_box(k,ibox) .AND. xp(k,i) .LE. max_box(k,ibox)) THEN
                     IF (ghost_req(k,i) .GT. max_ghost(k)) THEN
                        max_ghost(k) = ghost_req(k,i)
                     ENDIF
                  ENDIF
               ENDDO
            ENDDO
            ! Needs to collect max from all processors
            CALL MPI_Allreduce(max_ghost(1),minboxsizes(1,nsubs),ppm_dim,MPTYPE,MPI_MAX,ppm_comm,info)

         ENDDO
      ENDIF

      ! CHECK FOR has_one_way if we can subdivide again -> show a message
      IF (has_one_way) THEN

         ! go through all subs and check if you can subdivide it
         ghost_ok_one = .FALSE.
         nsubs_temp = 0

         ! alloc array to store sub index
         iopt   = ppm_param_alloc_fit
         ldc(1) = nsubs

         CALL ppm_alloc(sub_index,ldc,iopt,info)
         IF (info.NE.0) THEN
            info = ppm_error_fatal
            CALL ppm_error(ppm_err_alloc,'ppm_decomp_tree',               &
      &                  'alloc of sub_index failed!',__LINE__,info)
            GOTO 9999
         ENDIF 

         DO i=1,nsubs
            ! Determine if we can split
            ghost_ok = .TRUE.
            DO k=1,ppm_dim
               IF (minboxsizes_temp(k,i) .GT. 0.5*(max_sub(k,i) - min_sub(k,i))) THEN
                  ghost_ok = .FALSE.
               ENDIF
            ENDDO

            IF (numb_part(i).GT.0 .AND. ghost_ok) THEN
               ! Add this box to the again to be checked boxes.
               nsubs_temp = nsubs_temp + 1
               sub_index(nsubs_temp) = i
               ghost_ok_one = .TRUE.
            ENDIF

         ENDDO


         ! Check for variance
         IF (nsubs.GT.ppm_nproc) THEN
            !-------------------------------------------------------------------
            !  Compute the mean and variance of the particle in the boxes at
            !  the next level
            !-------------------------------------------------------------------
            var_npbx = 0.0_MK 
            mean_npbx = 0.0_MK 
           
            !----------------------------------------------------------------
            !  If the pcost is not present compute the cost as the sum of the
            !  (global) number of particles in the boxes
            !----------------------------------------------------------------
            DO i=1,nsubs
               mean_npbx = mean_npbx + REAL(numb_part(i),MK)
                  var_npbx =  var_npbx + REAL(numb_part(i),MK)**2
            ENDDO

            !-------------------------------------------------------------------
            !  Compute the mean and the variance
            !-------------------------------------------------------------------
            mean_npbx = mean_npbx/REAL(nsubs,MK)
             var_npbx = var_npbx/REAL(nsubs,MK) - mean_npbx**2

            !-------------------------------------------------------------------
            !  If the variance is acceptable ie less than the tolerance times 
            !  the mean number of particles, then stop 
            !-------------------------------------------------------------------
            variance_not_ok = .TRUE.
            IF (var_npbx.LT.tolerance*mean_npbx) THEN
               variance_not_ok = .FALSE.
            ENDIF 

         ENDIF

         ! Plot a message if we are not optimal
         IF (ghost_ok_one .AND. variance_not_ok) THEN
            IF(ppm_rank .EQ. 0) THEN
               print *, 'WE COULD STILL DIVIDE BOXES'
            ENDIF
         ENDIF
      ENDIF


      !-------------------------------------------------------------------------
      !  Debugging: write out the tree 
      !-------------------------------------------------------------------------
      IF (ppm_debug.GT.0) THEN
         IF (ppm_rank.EQ.0) THEN
            mesg = 'Writing the tree to the file: tree.dat'
            CALL ppm_write(ppm_rank,'ppm_decomp_tree',  mesg,info)
            OPEN(10,FILE='tree.dat')
            DO i=1,nsubs
               WRITE(10,'(2e12.4)') min_sub(1,i),min_sub(2,i)
               WRITE(10,'(2e12.4)') max_sub(1,i),min_sub(2,i)
               WRITE(10,'(2e12.4)') max_sub(1,i),max_sub(2,i)
               WRITE(10,'(2e12.4)') min_sub(1,i),max_sub(2,i)
               WRITE(10,'(2e12.4)') min_sub(1,i),min_sub(2,i)
               WRITE(10,'(A)')
            ENDDO
            CLOSE(10)
         ENDIF
      ENDIF 

      !-------------------------------------------------------------------------
      !  Free the memory again
      !-------------------------------------------------------------------------
      iopt = ppm_param_dealloc
      CALL ppm_alloc(ppb    ,ldc,iopt,info)
      IF (info.NE.0) GOTO 500
      CALL ppm_alloc(npbx   ,ldc,iopt,info)
      IF (info.NE.0) GOTO 500
      CALL ppm_alloc(npbxg  ,ldc,iopt,info)
      IF (info.NE.0) GOTO 500
      CALL ppm_alloc(min_box,ldc,iopt,info) 
      IF (info.NE.0) GOTO 500
      CALL ppm_alloc(max_box,ldc,iopt,info) 
      IF (info.NE.0) GOTO 500
      CALL ppm_alloc(min_box_temp,ldc,iopt,info) 
      IF (info.NE.0) GOTO 500
      CALL ppm_alloc(max_box_temp,ldc,iopt,info) 
      IF (info.NE.0) GOTO 500
      CALL ppm_alloc(minboxsizes_temp,ldc,iopt,info) 
      IF (info.NE.0) GOTO 500
      CALL ppm_alloc(max_ghost_temp,ldc,iopt,info) 
      IF (info.NE.0) GOTO 500
      CALL ppm_alloc(work   ,ldc,iopt,info) 
      IF (info.NE.0) GOTO 500
  500 CONTINUE
      IF (info.NE.0) THEN
         info = ppm_error_fatal
         CALL ppm_error(ppm_err_alloc,'ppm_decomp_tree',                  &
     &                  'deallocating work arrays failed!',__LINE__,info)
         GOTO 9999
      ENDIF 

      !-------------------------------------------------------------------------
      !  Return 
      !-------------------------------------------------------------------------
 9999 CONTINUE
      CALL substop('ppm_decomp_tree',t0,info)
      RETURN
      CONTAINS
      SUBROUTINE check
         ! check all minboxsizes
!          IF (minboxsize .LE. 0.0_MK) THEN
!             info = ppm_error_error
!             CALL ppm_error(ppm_err_argument,'ppm_decomp_tree',     &
!      &          'the minimum box size must be > 0 !',__LINE__,info)
!             GOTO 8888
!          ENDIF
         IF (Npart .LT. 0) THEN
            info = ppm_error_error
            CALL ppm_error(ppm_err_argument,'ppm_decomp_tree',     &
     &          'the number of particles must be >= 0 !',__LINE__,info)
            GOTO 8888
         ENDIF
         DO i=1,ppm_dim
            IF (min_phys(i) .GT. max_phys(i)) THEN
               info = ppm_error_error
               CALL ppm_error(ppm_err_argument,'ppm_decomp_tree',     &
     &             'min_phys must be <= max_phys !',__LINE__,info)
                GOTO 8888
            ENDIF
         ENDDO
 8888    CONTINUE
      END SUBROUTINE check
#if   __KIND == __SINGLE_PRECISION
      END SUBROUTINE decomp_tree_inhom_s
#elif __KIND == __DOUBLE_PRECISION
      END SUBROUTINE decomp_tree_inhom_d
#endif
!<<<< haeckic end >>>>!