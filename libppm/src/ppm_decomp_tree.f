      !-------------------------------------------------------------------------
      !  Subroutine   :                  ppm_decomp_tree
      !-------------------------------------------------------------------------
      !
      !  Purpose      : This routine performs a tree-like decomposition. 
      !                 It subdivides space until the number of leaves in the
      !                 tree exceeds the number of processors AND until the
      !                 variance of the number of particles in the leaves is
      !                 below a certain tolerance.
      !
      !  Input        : xp(:,:)      (F) the position of the particles
      !                 Npart        (I) the number of particles 
      !                 min_phys(:)  (F) the minimum coordinate of the 
      !                                  physical/computational domain 
      !                 max_phys(:)  (F) the maximum coordinate of the 
      !                                  physical/computational domain 
      !                 minboxsize   (F) the miminum box size 
      !                 tolerance    (F) the maximum relative variance of the 
      !                                  box cost allowed before terminating 
      !                                  the tree; tolerance >= 0
      !                 pcost(:)     (F) OPTIONAL argument of length
      !                                  Npart, specifying the
      !                                  computational cost of each
      !                                  particle.
      !
      !  Input/output :                                            
      !
      !  Output       : min_sub(:,:) (F) the min. extent of the subdomain
      !                 max_sub(:,:) (F) the max. extent of the subdomain
      !                 nsubs        (I) the total number of subdomains
      !                 info         (I) return status
      !
      !  Remarks      : 
      !
      !  References   :
      !
      !  Revisions    :
      !-------------------------------------------------------------------------
      !  $Log: ppm_decomp_tree.f,v $
      !  Revision 1.1.1.1  2007/07/13 10:18:54  ivos
      !  CBL version of the PPM library
      !
      !  Revision 1.15  2006/09/04 18:34:42  pchatela
      !  Fixes and cleanups to make it compile for
      !  release-1-0
      !
      !  Revision 1.14  2004/09/17 12:02:16  ivos
      !  Added argument checks for minboxsize, min_phys, max_phys and Npart.
      !
      !  Revision 1.13  2004/07/26 07:42:38  ivos
      !  Changed to single-interface modules and adapted all USE statements.
      !
      !  Revision 1.12  2004/06/10 16:19:59  ivos
      !  Moved all cpp directtives to column 1. The NEC cpp did not recognize
      !  them otherwise!!!
      !
      !  Revision 1.11  2004/04/15 10:07:39  ivos
      !  removed info=0 after substart.
      !
      !  Revision 1.10  2004/03/03 09:09:01  ivos
      !  Removed cost from the arguments since this is now computed by
      !  ppm_topo_cost. Added pcost to the arguments, but it is not yet used
      !  anywhere in the code. code which computed cost is now commented. Needs
      !  clean-up.
      !
      !  Revision 1.9  2004/02/27 17:25:04  walther
      !  Bug fix: after the realloc of the min_sub array, the code had a GOTO 400
      !  irrespectively of the info !-)
      !
      !  Revision 1.8  2004/02/06 13:52:45  walther
      !  bug fix #ifdef misplaced.
      !
      !  Revision 1.7  2004/02/02 14:03:05  walther
      !  Performing error checks.
      !
      !  Revision 1.6  2004/01/23 17:24:15  ivos
      !  Now includes ppm_define.h for the cpp defines and no longer sets them
      !  in the Makefile.
      !
      !  Revision 1.5  2004/01/09 17:25:38  walther
      !  Bug fixes in allocation of local and input/output arrays.
      !  Another bug fix was the overflow of the integer idiff, that
      !  has been change to float (fdiff).
      !
      !  Revision 1.4  2003/12/17 09:36:54  hiebers
      !  increased memory for npbx,npbxg (ldc(2)=2000, might need some revision)
      !
      !  Revision 1.3  2003/12/12 15:55:23  ivos
      !  Removed topoid from the argument list as it is not needed.
      !
      !  Revision 1.2  2003/12/09 13:40:37  hiebers
      !  added allocation of 'npbxg' in serial version
      !
      !  Revision 1.1.1.1  2003/11/17 15:13:45  walther
      !
      !-------------------------------------------------------------------------
      !  Parallel Particle Mesh Library (PPM)
      !  Institute of Computational Science
      !  ETH Zentrum, Hirschengraben 84
      !  CH-8092 Zurich, Switzerland
      !-------------------------------------------------------------------------

#if   __KIND == __SINGLE_PRECISION
      SUBROUTINE ppm_decomp_tree_s(xp,Npart,min_phys,max_phys, &
     &   minboxsize,tolerance,min_sub,max_sub,nsubs,info,pcost)
#elif __KIND == __DOUBLE_PRECISION
      SUBROUTINE ppm_decomp_tree_d(xp,Npart,min_phys,max_phys, &
     &   minboxsize,tolerance,min_sub,max_sub,nsubs,info,pcost)
#endif

      !-------------------------------------------------------------------------
      !  Modules 
      !-------------------------------------------------------------------------
      USE ppm_module_data
      USE ppm_module_substart
      USE ppm_module_substop
      USE ppm_module_error
      USE ppm_module_alloc
      USE ppm_module_write
      USE ppm_module_decomp_boxsplit
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
      REAL(MK), DIMENSION(:)  , INTENT(IN   ) :: min_phys,max_phys
      REAL(MK), DIMENSION(:)  , OPTIONAL, INTENT(IN) :: pcost
      REAL(MK)                , INTENT(IN   ) :: minboxsize,tolerance
      REAL(MK), DIMENSION(:,:), POINTER       :: min_sub,max_sub
      INTEGER                 , INTENT(IN   ) :: Npart
      INTEGER                 , INTENT(  OUT) :: nsubs
      INTEGER                 , INTENT(  OUT) :: info
      !-------------------------------------------------------------------------
      !  Local variables 
      !-------------------------------------------------------------------------
      REAL(MK), DIMENSION(ppm_dim) :: len_phys
      REAL(MK), DIMENSION(2)       :: vector_in,vector_out
      REAL(MK):: mean_npbx,var_npbx
      REAL(MK):: t0
      REAL(MK), DIMENSION(:,:), POINTER :: min_box,max_box,work
      INTEGER , DIMENSION(:), POINTER :: npbx,npbxg,ppb
      INTEGER , DIMENSION(ppm_dim) :: ldc,ldd
      INTEGER :: ibox,fbox,jbox,kbox,lbox,nbox,ilevel,nlevel,isize
      INTEGER :: i,j,k,iopt,ipart,ii,jj,mem_req
      INTEGER :: istat,n1,n2,Npartg
      LOGICAL :: lcontinue
      CHARACTER(ppm_char) :: mesg
      !-------------------------------------------------------------------------
      !  Externals 
      !-------------------------------------------------------------------------
      
      !-------------------------------------------------------------------------
      !  Initialise 
      !-------------------------------------------------------------------------
      CALL substart('ppm_decomp_tree',t0,info)

      !-------------------------------------------------------------------------
      !  check input arguments
      !-------------------------------------------------------------------------
      IF (ppm_debug.GT.0) THEN
         IF (minboxsize .LE. 0.0_MK) THEN
            info = ppm_error_error
            CALL ppm_error(ppm_err_argument,'ppm_decomp_tree',     &
     &          'the minimum box size must be > 0 !',__LINE__,info)
            GOTO 9999
         ENDIF 
         IF (Npart .LT. 0) THEN
            info = ppm_error_error
            CALL ppm_error(ppm_err_argument,'ppm_decomp_tree',     &
     &          'the number of particles must be >= 0 !',__LINE__,info)
            GOTO 9999
         ENDIF 
         DO i=1,ppm_dim
            IF (min_phys(i) .GT. max_phys(i)) THEN
               info = ppm_error_error
               CALL ppm_error(ppm_err_argument,'ppm_decomp_tree',     &
     &             'min_phys must be <= max_phys !',__LINE__,info)
                GOTO 9999
            ENDIF 
         ENDDO
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

      ldc(1) = ldc(2)
      CALL ppm_alloc(ppb  ,ldc,iopt,info)
      IF (info.NE.0) THEN
         info = ppm_error_fatal
         CALL ppm_error(ppm_err_alloc,'ppm_decomp_tree',               &
     &                  'alloc of ppb failed!',__LINE__,info)
         GOTO 9999
      ENDIF 

      CALL ppm_alloc(npbx ,ldc,iopt,info)
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
      DO k=1,ppm_dim
         min_box(k,fbox) = min_phys(k)
         max_box(k,fbox) = max_phys(k)
         len_phys(k)     = max_phys(k) - min_phys(k)
      ENDDO
      ppb(fbox)   = 1
      npbx(fbox)  = Npart
      npbxg(fbox) = Npart

      !-------------------------------------------------------------------------
      !  Compute the maximum number of levels 
      !-------------------------------------------------------------------------
      nlevel = INT(LOG(MINVAL(len_phys)/minboxsize)/LOG(2.0_MK)) + 1

      !-------------------------------------------------------------------------
      !  Split the boxes, starting at level 0 (the full box)
      !  kbox denotes the box being divided
      !  nbox is the current total number of boxes
      !  fbox is the id of the first box on the current level
      !  lbox is the id of the last  box on the current level
      !-------------------------------------------------------------------------
      ilevel    = 0
      lcontinue = .TRUE.
      lbox      = fbox
      nbox      = lbox
      nsubs     = 0
      DO WHILE (lcontinue)
         !----------------------------------------------------------------------
         !  Split the boxes at this level
         !----------------------------------------------------------------------
         DO kbox=fbox,lbox
            !-------------------------------------------------------------------
            !  Split non-empty boxes into four/eight (nbox is incremented)
            !-------------------------------------------------------------------
            IF (npbxg(kbox).GT.0) THEN
               CALL ppm_decomp_boxsplit(work,ppb,npbx,kbox,nbox, &
     &                                  min_box,max_box,info)
               IF (info.NE.0) GOTO 100
            ELSE
               !----------------------------------------------------------------
               !  if empty save it as a subdomain ie increment the sub counter
               !  and check for memory 
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
                  ldc(1) = ldc(2)
           !       CALL ppm_alloc(cost,ldc,iopt,info)
           !       IF (info.NE.0) GOTO 200
               ENDIF 

               !----------------------------------------------------------------
               !  save the subdomain
               !----------------------------------------------------------------
               DO k=1,ppm_dim
                  min_sub(k,nsubs) = min_box(k,kbox)
                  max_sub(k,nsubs) = max_box(k,kbox)
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
         ilevel = ilevel + 1
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
         IF (ilevel.EQ.nlevel) lcontinue = .FALSE.

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
         ldc(1) = mem_req
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
         ENDDO
      ENDIF

      !-------------------------------------------------------------------------
      !  Debugging: write out the tree 
      !-------------------------------------------------------------------------
      IF (ppm_debug.GT.0) THEN
         IF (ppm_rank.EQ.0) THEN
            mesg = 'Writing the tree to the file: tree.dat'
            CALL ppm_write(ppm_rank,'ppm_decomp_tree',mesg,info)
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
#if   __KIND == __SINGLE_PRECISION
      END SUBROUTINE ppm_decomp_tree_s
#elif __KIND == __DOUBLE_PRECISION
      END SUBROUTINE ppm_decomp_tree_d
#endif
