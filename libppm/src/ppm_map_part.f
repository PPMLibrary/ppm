      !-------------------------------------------------------------------------
      !  Subroutine   :                     ppm_map_part
      !-------------------------------------------------------------------------
      !
      !  Purpose      : This routine maps the particles onto a topology.
      !
      !  Input        : xp(:[,:])  (O) the position of the particles or the
      !                                data to be pushed/popped. Data can
      !                                be either a 1d or a 2d array of type
      !                                single,double,integer,logical,single
      !                                conplex or double complex. Positions
      !                                need to be a 2d array of type single
      !                                or double.
      !                 [lda]      (I) the leading dimension of the data (only
      !                                important for push and pop) in case
      !                                xp is a 2d array. For 1d xp arrays,
      !                                omit this argument.
      !                 Npart      (I) the number of particles (on processor)
      !                 to_topo    (I) user topology ID of destination. If 
      !                                is <=0, the current topology is used.
      !                 maptype    (I) Mapping action requested. One of:
      !                                     ppm_param_map_global
      !                                     ppm_param_map_partial
      !                                     ppm_param_map_remap
      !                                     ppm_param_map_push
      !                                     ppm_param_map_send
      !                                     ppm_param_map_pop
      !                                     ppm_param_map_cancel
      !
      !  Input/output : Mpart      (I) the new number of particles after 
      !                                sending to new topology (on processor)
      !
      !  Output       : info       (I) return status. 0 upon success.
      !
      !  Remarks      :
      !
      !  References   :
      !
      !  Revisions    :
      !-------------------------------------------------------------------------
      !  $Log: ppm_map_part.f,v $
      !  Revision 1.1.1.1  2007/07/13 10:18:56  ivos
      !  CBL version of the PPM library
      !
      !  Revision 1.41  2006/09/04 18:34:49  pchatela
      !  Fixes and cleanups to make it compile for
      !  release-1-0
      !
      !  Revision 1.40  2005/03/10 01:40:33  ivos
      !  Empty buffer is now only reported for ppm_debug.GT.0 and is no
      !  longer logged to prevent huge log files.
      !
      !  Revision 1.39  2004/08/31 12:14:54  ivos
      !  Changed to use ppm_check_topoid and ppm_check_meshid to check validity
      !  of user-numbered IDs.
      !
      !  Revision 1.38  2004/08/30 09:27:29  ivos
      !  corrected indentation.
      !
      !  Revision 1.37  2004/07/26 15:38:49  ivos
      !  Inserted missing USE statements to resolve undefined references
      !  at link stage.
      !
      !  Revision 1.36  2004/07/26 07:42:44  ivos
      !  Changed to single-interface modules and adapted all USE statements.
      !
      !  Revision 1.35  2004/07/19 07:41:53  ivos
      !  Overloaded particle push and pop operations for 1d data arrays.
      !  Added new routines to module and needed checks to interface
      !  routine.
      !
      !  Revision 1.34  2004/07/16 14:46:27  ivos
      !  Added check for ppm_initialized.
      !
      !  Revision 1.33  2004/07/16 14:07:59  ivos
      !  Added sequence and argument checks. New checks now allow multiple
      !  push-send-pop cycles per mapping.
      !
      !  Revision 1.32  2004/06/16 09:05:21  ivos
      !  Added ppm_map_part_ring also for map type remap.
      !
      !  Revision 1.31  2004/05/28 11:45:51  walther
      !  Added Npart to the arguments of ppm_map_part_send().
      !
      !  Revision 1.30  2004/05/17 15:47:15  oingo
      !  Changed the check for Npart from .LE. 0 to .LT. 0
      !
      !  Revision 1.29  2004/04/22 10:35:55  ivos
      !  bugfix: use of SIZE(ppm_internal_topid) caused problems in argument
      !  checking when first user-topoid was > 1. Resolved by use of UBOUND.
      !
      !  Revision 1.28  2004/04/20 11:56:15  oingo
      !  Adjusted some IF-statements such that topoid=0 will not be blocked
      !  anymore
      !
      !  Revision 1.27  2004/04/13 15:15:13  oingo
      !  The call of ppm_map_part_ring is now a subcase of map_part_global, 
      !  thus removed the case for ppm_param_map_ring
      !
      !  Revision 1.26  2004/04/07 12:03:48  oingo
      !  Added map_part_ring action type.
      !
      !  Revision 1.25  2004/04/01 14:12:07  ivos
      !  Added success check after CALL to ppm_util_commopt.
      !
      !  Revision 1.24  2004/02/24 17:12:40  ivos
      !  Added overloaded versions for single complex and double complex 
      !  particle data.
      !
      !  Revision 1.23  2004/02/20 15:45:14  ivos
      !  Added check if the current (source) topology is defined at all.
      !
      !  Revision 1.22  2004/02/19 15:47:24  walther
      !  Bug fix: added ppm_target_topoid
      !
      !  Revision 1.21  2004/02/18 14:35:59  walther
      !  Removed the ghosts and the map_update is now map_partial.
      !
      !  Revision 1.20  2004/02/17 16:10:28  ivos
      !  Corrected the comment for the map_send case.
      !
      !  Revision 1.19  2004/02/11 14:29:07  ivos
      !  bugfix: validity of to_topo is now only checked if a to_topo .GT. 0
      !  is given at all.
      !
      !  Revision 1.18  2004/02/10 14:23:24  ivos
      !  Removed unused local variables j,k and istat.
      !
      !  Revision 1.17  2004/01/26 12:37:32  ivos
      !  Changed call to ppm_map_part_push to match new argument list (removed
      !  topoid).
      !
      !  Revision 1.16  2004/01/23 17:24:16  ivos
      !  Now includes ppm_define.h for the cpp defines and no longer sets them
      !  in the Makefile.
      !
      !  Revision 1.15  2004/01/23 11:31:21  ivos
      !  Cleanup: (1) updated headers, (2) inserted ppm_error and ppm_write,
      !  (3) added argument checking, (4) added checks after every alloc.
      !
      !  Revision 1.14  2003/12/18 13:50:03  hiebers
      !  merged
      !
      !  Revision 1.13  2003/12/18 11:30:04  ivos
      !  Added #ifdefs for invalid data types in remap case to make it compile
      !  again.
      !
      !  Revision 1.12  2003/12/17 16:52:02  ivos
      !  Uncommented CALL to map_part_remap.
      !
      !  Revision 1.11  2003/12/17 12:08:16  walther
      !  Removed quotes from comments in header.
      !
      !  Revision 1.10  2003/12/17 09:31:44  ivos
      !  Default initialization of topoid added to prevent ifc uninitialized
      !  variable errors.
      !
      !  Revision 1.9  2003/12/16 13:36:07  ivos
      !  topo_id = 0 now means: current topology.
      !
      !  Revision 1.8  2003/12/16 12:23:01  ivos
      !  Added map_part_cancel action type.
      !
      !  Revision 1.7  2003/12/12 18:14:16  hiebers
      !  Merged
      !
      !  Revision 1.6  2003/12/12 18:01:54  ivos
      !  Added the following: (1) updates internal ppm_topoid upon completion 
      !  of a mapping. (2) removed topoid argument from partial map routine 
      !  as they always act on the current topology. (3) added translation 
      !  from user topo IDs to internal ones. (4) added sanity checks for 
      !  mappings in progress. (5) removed from_topo from argument list since 
      !  this is always the current topology. (6) added map type 
      !  ppm_param_map_cancel to cancel an unfinished mapping.
      !
      !  Revision 1.5  2003/12/09 08:59:59  ivos
      !  Now checks if the communication sequence has been optimized for this
      !  topology before the partial map is called. If not, the optimizer is 
      !  called.
      !
      !  Revision 1.4  2003/12/05 14:42:49  ivos
      !  Before map_part_partial is called it is checked if the optimal
      !  communication sequence has already been determined. If not, the 
      !  optimizer is called.
      !
      !  Revision 1.3  2003/12/05 11:54:05  ivos
      !  Added #ifdefs to allow map_partial only for double and single types.
      !
      !  Revision 1.2  2003/12/05 09:18:48  ivos
      !  Added CALL to ppm_map_part_partial.
      !
      !  Revision 1.1.1.1  2003/11/17 15:13:45  walther
      !
      !-------------------------------------------------------------------------
      !  Parallel Particle Mesh Library (PPM)
      !  Institute of Computational Science
      !  ETH Zentrum, Hirschengraben 84
      !  CH-8092 Zurich, Switzerland
      !-------------------------------------------------------------------------

#if    __DIM == 1
#if    __KIND == __SINGLE_PRECISION
      SUBROUTINE ppm_map_part_1ds(xp,Npart,Mpart,to_topo,maptype,info)
#elif  __KIND == __DOUBLE_PRECISION
      SUBROUTINE ppm_map_part_1dd(xp,Npart,Mpart,to_topo,maptype,info)
#elif  __KIND == __SINGLE_PRECISION_COMPLEX
      SUBROUTINE ppm_map_part_1dsc(xp,Npart,Mpart,to_topo,maptype,info)
#elif  __KIND == __DOUBLE_PRECISION_COMPLEX
      SUBROUTINE ppm_map_part_1ddc(xp,Npart,Mpart,to_topo,maptype,info)
#elif  __KIND == __INTEGER
      SUBROUTINE ppm_map_part_1di(xp,Npart,Mpart,to_topo,maptype,info)
#elif  __KIND == __LOGICAL
      SUBROUTINE ppm_map_part_1dl(xp,Npart,Mpart,to_topo,maptype,info)
#endif 

#elif  __DIM == 2
#if    __KIND == __SINGLE_PRECISION
      SUBROUTINE ppm_map_part_2ds(xp,lda,Npart,Mpart,to_topo,maptype,info)
#elif  __KIND == __DOUBLE_PRECISION
      SUBROUTINE ppm_map_part_2dd(xp,lda,Npart,Mpart,to_topo,maptype,info)
#elif  __KIND == __SINGLE_PRECISION_COMPLEX
      SUBROUTINE ppm_map_part_2dsc(xp,lda,Npart,Mpart,to_topo,maptype,info)
#elif  __KIND == __DOUBLE_PRECISION_COMPLEX
      SUBROUTINE ppm_map_part_2ddc(xp,lda,Npart,Mpart,to_topo,maptype,info)
#elif  __KIND == __INTEGER
      SUBROUTINE ppm_map_part_2di(xp,lda,Npart,Mpart,to_topo,maptype,info)
#elif  __KIND == __LOGICAL
      SUBROUTINE ppm_map_part_2dl(xp,lda,Npart,Mpart,to_topo,maptype,info)
#endif 
#endif 

      !-------------------------------------------------------------------------
      !  Includes
      !-------------------------------------------------------------------------
#include "ppm_define.h"

      !-------------------------------------------------------------------------
      !  Modules 
      !-------------------------------------------------------------------------
      USE ppm_module_data
      USE ppm_module_data_part, ONLY: ppm_target_topoid
      USE ppm_module_substart
      USE ppm_module_substop
      USE ppm_module_error
      USE ppm_module_write
      USE ppm_module_check_topoid
      USE ppm_module_util_commopt
      USE ppm_module_map_part_global
      USE ppm_module_map_part_partial
      USE ppm_module_map_part_ring
      USE ppm_module_map_part_remap
      USE ppm_module_map_part_push
      USE ppm_module_map_part_send
      USE ppm_module_map_part_pop
      IMPLICIT NONE
#if    __KIND == __SINGLE_PRECISION | __KIND == __SINGLE_PRECISION_COMPLEX
      INTEGER, PARAMETER :: MK = ppm_kind_single
#else
      INTEGER, PARAMETER :: MK = ppm_kind_double
#endif
      !-------------------------------------------------------------------------
      !  Arguments     
      !-------------------------------------------------------------------------
#if    __DIM == 1
#if    __KIND == __INTEGER
      INTEGER , DIMENSION(:  )  , POINTER       :: xp
#elif  __KIND == __LOGICAL
      LOGICAL , DIMENSION(:  )  , POINTER       :: xp
#elif  __KIND == __SINGLE_PRECISION_COMPLEX | __KIND == __DOUBLE_PRECISION_COMPLEX
      COMPLEX(MK), DIMENSION(:  )  , POINTER    :: xp
#else
      REAL(MK), DIMENSION(:  )  , POINTER       :: xp
#endif

#elif  __DIM == 2
#if    __KIND == __INTEGER
      INTEGER , DIMENSION(:,:)  , POINTER       :: xp
#elif  __KIND == __LOGICAL
      LOGICAL , DIMENSION(:,:)  , POINTER       :: xp
#elif  __KIND == __SINGLE_PRECISION_COMPLEX | __KIND == __DOUBLE_PRECISION_COMPLEX
      COMPLEX(MK), DIMENSION(:,:)  , POINTER    :: xp
#else
      REAL(MK), DIMENSION(:,:)  , POINTER       :: xp
#endif
#endif
#if    __DIM == 2
      INTEGER                   , INTENT(IN   ) :: lda
#endif
      INTEGER                   , INTENT(IN   ) :: Npart
      INTEGER                   , INTENT(IN   ) :: to_topo
      INTEGER                   , INTENT(  OUT) :: Mpart
      INTEGER                   , INTENT(IN   ) :: maptype
      INTEGER                   , INTENT(  OUT) :: info
      !-------------------------------------------------------------------------
      !  Local variables 
      !-------------------------------------------------------------------------
#if    __DIM == 1
      INTEGER, PARAMETER  :: lda = 1
#endif
      INTEGER             :: i
      REAL(MK)            :: t0
      CHARACTER(ppm_char) :: mesg
      LOGICAL             :: valid
      !-------------------------------------------------------------------------
      !  Externals 
      !-------------------------------------------------------------------------
      
      !-------------------------------------------------------------------------
      !  Initialize 
      !-------------------------------------------------------------------------
      CALL substart('ppm_map_part',t0,info)

      !-------------------------------------------------------------------------
      !  Check arguments
      !-------------------------------------------------------------------------
      IF (ppm_debug .GT. 0) THEN
          IF (.NOT. ppm_initialized) THEN
              info = ppm_error_error
              CALL ppm_error(ppm_err_ppm_noinit,'ppm_map_part',  &
     &            'Please call ppm_init first!',__LINE__,info)
              GOTO 9999
          ENDIF
          IF (ppm_topoid .LT. 0) THEN
              info = ppm_error_error
              CALL ppm_error(ppm_err_no_topo,'ppm_map_part',  &
     &            'No topology has been defined so far',__LINE__,info)
              GOTO 9999
          ENDIF
#if    __DIM == 2
          IF (lda .LT. 1) THEN
              info = ppm_error_error
              CALL ppm_error(ppm_err_argument,'ppm_map_part',  &
     &            'lda must be >0',__LINE__,info)
              GOTO 9999
          ENDIF
#endif
          IF (Npart .LT. 0) THEN
              info = ppm_error_error
              CALL ppm_error(ppm_err_argument,'ppm_map_part',  &
     &            'Npart must be >=0',__LINE__,info)
              GOTO 9999
          ENDIF
          IF (to_topo .GE. 0) THEN
              CALL ppm_check_topoid(ppm_param_id_user,to_topo,valid,info)
              IF (.NOT. valid) THEN
                  info = ppm_error_error
                  CALL ppm_error(ppm_err_argument,'ppm_map_part',  &
     &                'Topology ID (to_topo) is invalid!',__LINE__,info)
                  GOTO 9999
              ENDIF
          ENDIF
          IF ((to_topo.EQ.0) .AND. (maptype.EQ.ppm_param_map_partial)) THEN
             info = ppm_error_error
             CALL ppm_error(ppm_err_argument,'ppm_map_part', &
                 'Partial mapping not possible with topoid=0',__LINE__,info)
             GOTO 9999
          ENDIF
      ENDIF

      !-------------------------------------------------------------------------
      !  Determine which mapping is required
      !-------------------------------------------------------------------------
      IF     (maptype.EQ.ppm_param_map_global) THEN
         !----------------------------------------------------------------------
         !  Full blast mapping; find the particle outside the current processor
         !----------------------------------------------------------------------
#if __KIND == __INTEGER | \
    __KIND == __LOGICAL | \
    __KIND == __SINGLE_PRECISION_COMPLEX | \
    __KIND == __DOUBLE_PRECISION_COMPLEX
         info = ppm_error_error
         CALL ppm_error(ppm_err_wrong_type,'ppm_map_part',  &
     &       'Wrong type passed for global mapping',__LINE__,info)
         GOTO 9999
#else
#if __DIM == 1
         info = ppm_error_error
         CALL ppm_error(ppm_err_wrong_type,'ppm_map_part',  &
     &       'Particle positions must be a 2d array',__LINE__,info)
         GOTO 9999
#else
         ! if there is still some data left in the buffer, warn the user
         IF (ppm_buffer_set .GT. 0) THEN
             info = ppm_error_warning
             CALL ppm_error(ppm_err_map_incomp,'ppm_map_part',  &
     &           'Buffer was not empty. Possible loss of data!',__LINE__,info)
         ENDIF

         ! get new destination topology from user
         IF (to_topo .LT. 0) THEN
             ppm_target_topoid = ppm_topoid
         ELSE
             ppm_target_topoid = ppm_internal_topoid(to_topo)
         ENDIF

         !----------------------------------------------------------------------
         !  For a topoid of 0 use the ring mapping, the global mapping otherwise
         !----------------------------------------------------------------------
         IF (ppm_target_topoid .EQ. 0) THEN
             CALL ppm_map_part_ring(xp,Npart,info)
         ELSE
             CALL ppm_map_part_global(xp,Npart,ppm_target_topoid,info)
         ENDIF

         IF (info.NE.0) GOTO 9999
#endif 
#endif 

      ELSEIF (maptype.EQ.ppm_param_map_partial) THEN
         !----------------------------------------------------------------------
         !  partial mapping: only to neighboring processors
         !----------------------------------------------------------------------
#if __KIND == __INTEGER | \
    __KIND == __LOGICAL | \
    __KIND == __SINGLE_PRECISION_COMPLEX | \
    __KIND == __DOUBLE_PRECISION_COMPLEX
         info = ppm_error_error
         CALL ppm_error(ppm_err_wrong_type,'ppm_map_part',  &
     &       'Wrong type passed for partial mapping',__LINE__,info)
         GOTO 9999
#else
#if __DIM == 1
         info = ppm_error_error
         CALL ppm_error(ppm_err_wrong_type,'ppm_map_part',  &
     &       'Particle positions must be a 2d array',__LINE__,info)
         GOTO 9999
#else
         ! if there is still some data left in the buffer, warn the user
         IF (ppm_buffer_set .GT. 0) THEN
             info = ppm_error_warning
             CALL ppm_error(ppm_err_map_incomp,'ppm_map_part',  &
     &           'Buffer was not empty. Possible loss of data!',__LINE__,info)
         ENDIF

         ! new topoid is current topology
         ppm_target_topoid = ppm_topoid

         ! first check if the optimal communication protocol is known
         IF (.NOT. ppm_isoptimized(ppm_target_topoid)) THEN
             ! if not: determine it before calling map_part_partial
             CALL ppm_util_commopt(ppm_target_topoid,info) 
             IF (info.NE.0) GOTO 9999
             IF (ppm_debug .GT. 1) THEN
                DO i=1,ppm_nneighlist(ppm_target_topoid)
                    WRITE(mesg,'(A,I4)') 'have neighbor: ',   &
     &                  ppm_ineighlist(i,ppm_target_topoid)
                    CALL ppm_write(ppm_rank,'ppm_map_part',mesg,info)
                END DO
                DO i=1,ppm_ncommseq(ppm_target_topoid)
                    WRITE(mesg,'(A,I4)') 'communicate: ',   &
     &                  ppm_icommseq(i,ppm_target_topoid)
                    CALL ppm_write(ppm_rank,'ppm_map_part',mesg,info)
                END DO
             ENDIF
         END IF

         ! call the partial map
         CALL ppm_map_part_partial(xp,Npart,info)
         IF (info.NE.0) GOTO 9999
#endif
#endif

      ELSEIF (maptype.EQ.ppm_param_map_remap) THEN
         !----------------------------------------------------------------------
         !  Remap from one topology to another
         !----------------------------------------------------------------------

#if __KIND == __INTEGER | \
    __KIND == __LOGICAL | \
    __KIND == __SINGLE_PRECISION_COMPLEX | \
    __KIND == __DOUBLE_PRECISION_COMPLEX
         info = ppm_error_error
         CALL ppm_error(ppm_err_wrong_type,'ppm_map_part',  &
     &       'Wrong type passed for remapping',__LINE__,info)
         GOTO 9999
#else
#if __DIM == 1
         info = ppm_error_error
         CALL ppm_error(ppm_err_wrong_type,'ppm_map_part',  &
     &       'Particle positions must be a 2d array',__LINE__,info)
         GOTO 9999
#else
         ! if there is still some data left in the buffer, warn the user
         IF (ppm_buffer_set .GT. 0) THEN
             info = ppm_error_warning
             CALL ppm_error(ppm_err_map_incomp,'ppm_map_part',  &
     &           'Buffer was not empty. Possible loss of data!',__LINE__,info)
         ENDIF

         ! get new destination topology from user
         IF (to_topo .LT. 0) THEN
             ppm_target_topoid = ppm_topoid
         ELSE
             ppm_target_topoid = ppm_internal_topoid(to_topo)
         ENDIF
     
         !----------------------------------------------------------------------
         !  For a topoid of 0 use the ring mapping, remapping otherwise
         !----------------------------------------------------------------------
         IF (ppm_target_topoid .EQ. 0) THEN
             CALL ppm_map_part_ring(xp,Npart,info)
         ELSE
             CALL ppm_map_part_remap(xp,Npart,ppm_target_topoid,info)
         ENDIF
         IF (info.NE.0) GOTO 9999
#endif
#endif

      ELSEIF (maptype.EQ.ppm_param_map_push) THEN
         !----------------------------------------------------------------------
         !  Add the data to the stack (push)
         !----------------------------------------------------------------------
         ! warn the user if no mapping has yet been defined
         IF (ppm_target_topoid .LT. 0) THEN
             info = ppm_error_warning
             CALL ppm_error(ppm_err_nomap,'ppm_map_part',  &
     &           'Skipping push!',__LINE__,info)
             GOTO 9999
         ENDIF

         CALL ppm_map_part_push(xp,lda,Npart,info)
         IF (info.NE.0) GOTO 9999

      ELSEIF (maptype.EQ.ppm_param_map_pop) THEN
         !----------------------------------------------------------------------
         !  Extract the data from the stack (pop)
         !----------------------------------------------------------------------
         ! warn the user if no mapping has yet been defined
         IF (ppm_target_topoid .LT. 0) THEN
             info = ppm_error_warning
             CALL ppm_error(ppm_err_nomap,'ppm_map_part',  &
     &           'Skipping pop!',__LINE__,info)
             GOTO 9999
         ENDIF

         ! skip if the buffer is empty
         IF (ppm_buffer_set .LT. 1) THEN
             IF (ppm_debug .GT. 1) THEN
                 CALL ppm_write(ppm_rank,'ppm_map_part',  &
     &               'Buffer is empty: skipping pop!',info)
             ENDIF
             GOTO 9999
         ENDIF

         CALL ppm_map_part_pop(xp,lda,Npart,Mpart,info)
         IF (info.NE.0) GOTO 9999

      ELSEIF (maptype.EQ.ppm_param_map_send) THEN
         !----------------------------------------------------------------------
         !  Send/recv the packages 
         !----------------------------------------------------------------------
         ! warn the user if no mapping has yet been defined
         IF (ppm_target_topoid .LT. 0) THEN
             info = ppm_error_warning
             CALL ppm_error(ppm_err_nomap,'ppm_map_part',  &
     &           'Skipping send!',__LINE__,info)
             GOTO 9999
         ENDIF

         ! skip if the buffer is empty
         IF (ppm_buffer_set .LT. 1) THEN
             IF (ppm_debug .GT. 1) THEN
                 CALL ppm_write(ppm_rank,'ppm_map_part',  &
     &               'Buffer is empty: skipping send!',info)
             ENDIF
             GOTO 9999
         ENDIF

         CALL ppm_map_part_send(Npart,Mpart,info)
         IF (info.NE.0) GOTO 9999

         ! we are on the new topology now !
         ppm_topoid = ppm_target_topoid

      ELSEIF (maptype.EQ.ppm_param_map_cancel) THEN
         !----------------------------------------------------------------------
         !  Cancel the mapping in progress and reset everything. No need to
         !  reset arrays since thay are alloc_fit-ed the next time a
         !  mapping is called. Their deallocation only happens in
         !  ppm_finalize anyway.
         !----------------------------------------------------------------------
         ppm_target_topoid = -1
         !----------------------------------------------------------------------
         !  These might actually not be needed since the mapping routines
         !  would in principle be expected to reset them. They are here for
         !  safety reasons.
         !----------------------------------------------------------------------
         ppm_buffer_set    = 0
         ppm_nsendbuffer   = 0
         ppm_nrecvbuffer   = 0

      ELSE
         !----------------------------------------------------------------------
         !  Unknow mapping
         !----------------------------------------------------------------------
         info = ppm_error_error
         CALL ppm_error(ppm_err_argument,'ppm_map_part',    &
     &       'Unknown mapping action specified',__LINE__,info)
         GOTO 9999
      ENDIF

      !-------------------------------------------------------------------------
      !  Return 
      !-------------------------------------------------------------------------
 9999 CONTINUE
      CALL substop('ppm_map_part',t0,info)
      RETURN
#if    __DIM == 1
#if    __KIND == __SINGLE_PRECISION
      END SUBROUTINE ppm_map_part_1ds
#elif  __KIND == __DOUBLE_PRECISION
      END SUBROUTINE ppm_map_part_1dd
#elif  __KIND == __SINGLE_PRECISION_COMPLEX
      END SUBROUTINE ppm_map_part_1dsc
#elif  __KIND == __DOUBLE_PRECISION_COMPLEX
      END SUBROUTINE ppm_map_part_1ddc
#elif  __KIND == __LOGICAL
      END SUBROUTINE ppm_map_part_1dl
#elif  __KIND == __INTEGER
      END SUBROUTINE ppm_map_part_1di
#endif

#elif  __DIM == 2
#if    __KIND == __SINGLE_PRECISION
      END SUBROUTINE ppm_map_part_2ds
#elif  __KIND == __DOUBLE_PRECISION
      END SUBROUTINE ppm_map_part_2dd
#elif  __KIND == __SINGLE_PRECISION_COMPLEX
      END SUBROUTINE ppm_map_part_2dsc
#elif  __KIND == __DOUBLE_PRECISION_COMPLEX
      END SUBROUTINE ppm_map_part_2ddc
#elif  __KIND == __LOGICAL
      END SUBROUTINE ppm_map_part_2dl
#elif  __KIND == __INTEGER
      END SUBROUTINE ppm_map_part_2di
#endif
#endif
