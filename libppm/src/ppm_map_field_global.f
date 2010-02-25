      !-------------------------------------------------------------------------
      !  Subroutine   :                 ppm_map_field_global
      !-------------------------------------------------------------------------
      !
      !  Purpose      : This routine maps field data between two topologies
      !                 (which however need compatible meshes defined on
      !                 them) using a global mapping (i.e. every processor
      !                 communicates with every other). Source mesh must be
      !                 on the current field topology. Global lists with
      !                 all mesh blocks that have to be send and/or
      !                 received are built in this routine. Push, pop and
      !                 send will use these lists.
      !
      !  Input        : topoid    (I) topology identifier (internal numbering)
      !                               of target
      !                 frommesh  (I) mesh identifier (internal numbering)
      !                               of source
      !                 tomesh    (I) mesh identifier (internal numbering)
      !                               of target
      !
      !  Input/output : 
      !
      !  Output       : info      (I)  return status. 0 on success.
      !
      !  Remarks      : The first part of the send/recv lists contains the
      !                 on-processor data.
      !
      !                 Side effect: this routine uses the same global
      !                 send/recv buffers, pointers and lists as the
      !                 particle mapping routines (reason: these buffers
      !                 can be large => memory issues). Field and particle
      !                 mappings can therefore never overlap, but one must
      !                 be completed before the other starts.
      !
      !                 A map_field_partial could be constructed as
      !                 follows: every processor known both isendlist and
      !                 irecvlist. They could now use util_commopt to find
      !                 a good communication sequence and compress the
      !                 send/recv loop. Also, the loop when building the
      !                 local send/recv lists below would just need to go
      !                 over the neighbors of a sub instead of all subs.
      !                 Difficulty: the order of mesh
      !                 blocks in the send and recv buffer could not match
      !                 up any more and push/pop could fail. CHECK THIS!
      !
      !  References   :
      !
      !  Revisions    :
      !-------------------------------------------------------------------------
      !  $Log: ppm_map_field_global.f,v $
      !  Revision 1.1.1.1  2007/07/13 10:18:56  ivos
      !  CBL version of the PPM library
      !
      !  Revision 1.17  2004/11/11 15:26:16  ivos
      !  Moved allocatable work arrays to module.
      !
      !  Revision 1.16  2004/10/01 16:33:36  ivos
      !  cosmetics.
      !
      !  Revision 1.15  2004/10/01 16:09:04  ivos
      !  Replaced REAL(ppm_kind_double) :: t0 with REAL(MK) t0.
      !
      !  Revision 1.14  2004/08/31 12:48:09  ivos
      !  Changed argument checks to use ppm_check_topoid and ppm_check_meshid.
      !
      !  Revision 1.13  2004/07/26 15:38:48  ivos
      !  Inserted missing USE statements to resolve undefined references
      !  at link stage.
      !
      !  Revision 1.12  2004/07/26 07:42:41  ivos
      !  Changed to single-interface modules and adapted all USE statements.
      !
      !  Revision 1.11  2004/04/08 16:49:16  hiebers
      !  Bug fix: two arguments in a subroutine call were mixed up
      !
      !  Revision 1.10  2004/04/07 15:33:37  ivos
      !  Changed to use ppm_mesh_block_intersect. Tested.
      !
      !  Revision 1.9  2004/04/05 12:03:40  ivos
      !  Udated comment header and code comments as ppm_map_type is now
      !  actually used.
      !
      !  Revision 1.8  2004/04/02 15:23:11  ivos
      !  Changed to use only global sub IDs in the ppm_mesh_*sub lists
      !  (no more local IDs, because of ghosts).
      !
      !  Revision 1.7  2004/04/01 07:32:04  ivos
      !  Updated comment header.
      !
      !  Revision 1.6  2004/02/24 15:45:31  ivos
      !  More bugs fixed. Tested to work on 1,2,3,4,5 processors with 1d,2d real
      !  and integer fields and all of them mixed.
      !
      !  Revision 1.5  2004/02/23 12:19:00  ivos
      !  Several bugs fixed. Tested on 2 processors with a scalar field.
      !  Added debug output in several places.
      !
      !  Revision 1.4  2004/02/18 09:48:02  ivos
      !  Updated header comment.
      !
      !  Revision 1.3  2004/02/17 16:10:02  ivos
      !  Now also builds the receive lists a priori (easier send and pop).
      !
      !  Revision 1.2  2004/02/12 14:47:33  ivos
      !  Now uses the same buffers and lists as map_part (save memory).
      !
      !  Revision 1.1  2004/02/11 14:34:37  ivos
      !  Initial implementation. Not tested yet.
      !
      !-------------------------------------------------------------------------
      !  Parallel Particle Mesh Library (PPM)
      !  Institute of Computational Science
      !  ETH Zentrum, Hirschengraben 84
      !  CH-8092 Zurich, Switzerland
      !-------------------------------------------------------------------------

      SUBROUTINE ppm_map_field_global(topoid,frommesh,tomesh,info)

      !-------------------------------------------------------------------------
      !  Includes
      !-------------------------------------------------------------------------
#include "ppm_define.h"

      !-------------------------------------------------------------------------
      !  Modules 
      !-------------------------------------------------------------------------
      USE ppm_module_data
      USE ppm_module_data_mesh
      USE ppm_module_substart
      USE ppm_module_substop
      USE ppm_module_error
      USE ppm_module_alloc
      USE ppm_module_write
      USE ppm_module_check_topoid
      USE ppm_module_check_meshid
      USE ppm_module_mesh_block_intersect
      IMPLICIT NONE
      !-------------------------------------------------------------------------
      !  Arguments     
      !-------------------------------------------------------------------------
      INTEGER                 , INTENT(IN   ) :: topoid,frommesh,tomesh
      INTEGER                 , INTENT(  OUT) :: info
      !-------------------------------------------------------------------------
      !  Local variables 
      !-------------------------------------------------------------------------
      INTEGER, DIMENSION(3)            :: ldu
      INTEGER, DIMENSION(ppm_dim)      :: iblockstart,nblocksize,offset,  &
     &                                    ghostsize
      INTEGER                          :: i,j,idom,sendrank,recvrank
      INTEGER                          :: iopt,iset,ibuffer,pdim
      INTEGER                          :: fromtopo,totopo,nsendlist,nsend
      INTEGER                          :: nrecvlist
      CHARACTER(ppm_char)              :: mesg
      REAL(ppm_kind_double)            :: t0
      LOGICAL                          :: valid
      !-------------------------------------------------------------------------
      !  Externals 
      !-------------------------------------------------------------------------
      
      !-------------------------------------------------------------------------
      !  Initialise 
      !-------------------------------------------------------------------------
      CALL substart('ppm_map_field_global',t0,info)
      pdim = ppm_dim

      !-------------------------------------------------------------------------
      !  Check arguments
      !-------------------------------------------------------------------------
      IF (ppm_debug .GT. 0) THEN
          CALL ppm_check_topoid(ppm_param_id_internal,topoid,valid,info)
          IF (.NOT. valid) THEN
              info = ppm_error_error
              CALL ppm_error(ppm_err_argument,'ppm_map_field_global',  &
     &            'topoid out of range',__LINE__,info)
              GOTO 9999
          ENDIF
          CALL ppm_check_meshid(ppm_param_id_internal,frommesh,     &
     &        ppm_field_topoid,valid,info)
          IF (.NOT. valid) THEN
              info = ppm_error_error
              CALL ppm_error(ppm_err_argument,'ppm_map_field_global',  &
     &            'source meshid out of range',__LINE__,info)
              GOTO 9999
          ENDIF
          CALL ppm_check_meshid(ppm_param_id_internal,tomesh,topoid,valid,info)
          IF (.NOT. valid) THEN
              info = ppm_error_error
              CALL ppm_error(ppm_err_argument,'ppm_map_field_global',  &
     &            'destination meshid out of range',__LINE__,info)
              GOTO 9999
          ENDIF
      ENDIF

      !-------------------------------------------------------------------------
      !  Save the map type for the subsequent calls (used for checks!)
      !-------------------------------------------------------------------------
      ppm_map_type = ppm_param_map_global

      !-------------------------------------------------------------------------
      !  Origin topology is current topology of fields, destination is
      !  argument
      !-------------------------------------------------------------------------
      fromtopo = ppm_field_topoid
      totopo   = topoid

      !-------------------------------------------------------------------------
      !  Check if origin and target meshes are compatible (i.e. have the
      !  same number of grid points in the whole comput. domain)
      !-------------------------------------------------------------------------
      DO i=1,pdim
          IF (ppm_cart_mesh(frommesh,fromtopo)%Nm(i) .NE.  &
     &        ppm_cart_mesh(tomesh,totopo)%Nm(i)) THEN
              info = ppm_error_error
              CALL ppm_error(ppm_err_bad_mesh,'ppm_map_field_global',  &
     &            'source and destination meshes are incompatible',__LINE__,info)
              GOTO 9999
          ENDIF
      ENDDO

      !-------------------------------------------------------------------------
      !  Allocate memory for the local temporary sendlists
      !-------------------------------------------------------------------------
      iopt   = ppm_param_alloc_fit
      ldu(1) = ppm_nsublist(fromtopo)
      CALL ppm_alloc(isendfromsub,ldu,iopt,info)
      IF (info .NE. 0) THEN
          info = ppm_error_fatal
          CALL ppm_error(ppm_err_alloc,'ppm_map_field_global',     &
     &        'local send source sub list ISENDFROMSUB',__LINE__,info)
          GOTO 9999
      ENDIF
      CALL ppm_alloc(isendtosub,ldu,iopt,info)
      IF (info .NE. 0) THEN
          info = ppm_error_fatal
          CALL ppm_error(ppm_err_alloc,'ppm_map_field_global',     &
     &        'local send destination sub list ISENDTOSUB',__LINE__,info)
          GOTO 9999
      ENDIF
      ldu(1) = pdim
      ldu(2) = ppm_nsublist(fromtopo)
      CALL ppm_alloc(isendblkstart,ldu,iopt,info)
      IF (info .NE. 0) THEN
          info = ppm_error_fatal
          CALL ppm_error(ppm_err_alloc,'ppm_map_field_global',     &
     &        'local send block start list ISENDBLKSTART',__LINE__,info)
          GOTO 9999
      ENDIF
      CALL ppm_alloc(isendblksize,ldu,iopt,info)
      IF (info .NE. 0) THEN
          info = ppm_error_fatal
          CALL ppm_error(ppm_err_alloc,'ppm_map_field_global',     &
     &        'local send block size list ISENDBLKSIZE',__LINE__,info)
          GOTO 9999
      ENDIF
      CALL ppm_alloc(ioffset,ldu,iopt,info)
      IF (info .NE. 0) THEN
          info = ppm_error_fatal
          CALL ppm_error(ppm_err_alloc,'ppm_map_field_global',     &
     &        'local send block offset list IOFFSET',__LINE__,info)
          GOTO 9999
      ENDIF

      !-------------------------------------------------------------------------
      !  Find intersecting mesh domains to be sent
      !-------------------------------------------------------------------------
      nsendlist = 0
      ghostsize = 0
      offset    = 0
      DO i=1,ppm_nsublist(fromtopo)
          idom = ppm_isublist(i,fromtopo)
          DO j=1,ppm_nsubs(totopo)
              CALL ppm_mesh_block_intersect(idom,j,frommesh,tomesh, fromtopo, &
     &            totopo,offset,ghostsize,nsendlist,isendfromsub,isendtosub,  &
     &            isendblkstart,isendblksize,ioffset,info)
              IF (info .NE. 0) GOTO 9999
          ENDDO
      ENDDO

      !-------------------------------------------------------------------------
      !  Allocate memory for the local temporary receive lists
      !-------------------------------------------------------------------------
      iopt   = ppm_param_alloc_fit
      ldu(1) = ppm_nsublist(totopo)
      CALL ppm_alloc(irecvfromsub,ldu,iopt,info)
      IF (info .NE. 0) THEN
          info = ppm_error_fatal
          CALL ppm_error(ppm_err_alloc,'ppm_map_field_global',     &
     &        'local recv source sub list IRECVFROMSUB',__LINE__,info)
          GOTO 9999
      ENDIF
      CALL ppm_alloc(irecvtosub,ldu,iopt,info)
      IF (info .NE. 0) THEN
          info = ppm_error_fatal
          CALL ppm_error(ppm_err_alloc,'ppm_map_field_global',     &
     &        'local recv destination sub list IRECVTOSUB',__LINE__,info)
          GOTO 9999
      ENDIF
      ldu(1) = pdim
      ldu(2) = ppm_nsublist(totopo)
      CALL ppm_alloc(irecvblkstart,ldu,iopt,info)
      IF (info .NE. 0) THEN
          info = ppm_error_fatal
          CALL ppm_error(ppm_err_alloc,'ppm_map_field_global',     &
     &        'local recv block start list IRECVBLKSTART',__LINE__,info)
          GOTO 9999
      ENDIF
      CALL ppm_alloc(irecvblksize,ldu,iopt,info)
      IF (info .NE. 0) THEN
          info = ppm_error_fatal
          CALL ppm_error(ppm_err_alloc,'ppm_map_field_global',     &
     &        'local recv block size list IRECVBLKSIZE',__LINE__,info)
          GOTO 9999
      ENDIF

      !-------------------------------------------------------------------------
      !  Find intersecting mesh domains to be received
      !-------------------------------------------------------------------------
      nrecvlist = 0
      ! loop over fromtopo first and THEN over totopo in order to get the
      ! same ordering of mesh blocks as in the sendlist. This is crutial
      ! for the push and the pop to work properly !!!
      DO j=1,ppm_nsubs(fromtopo)
          DO i=1,ppm_nsublist(totopo)
              idom = ppm_isublist(i,totopo)
              CALL ppm_mesh_block_intersect(j,idom,frommesh,tomesh, fromtopo, &
     &           totopo,offset,ghostsize,nrecvlist,irecvfromsub,irecvtosub,  &
     &            irecvblkstart,irecvblksize,ioffset,info)
              IF (info .NE. 0) GOTO 9999
          ENDDO
      ENDDO

      !-------------------------------------------------------------------------
      !  Allocate memory for the global mesh sendlists
      !-------------------------------------------------------------------------
      iopt   = ppm_param_alloc_fit
      ldu(1) = nsendlist
      CALL ppm_alloc(ppm_mesh_isendfromsub,ldu,iopt,info)
      IF (info .NE. 0) THEN
          info = ppm_error_fatal
          CALL ppm_error(ppm_err_alloc,'ppm_map_field_global',     &
     &        'source send sub list PPM_MESH_ISENDFROMSUB',__LINE__,info)
          GOTO 9999
      ENDIF
      ldu(1) = pdim
      ldu(2) = nsendlist
      CALL ppm_alloc(ppm_mesh_isendblkstart,ldu,iopt,info)
      IF (info .NE. 0) THEN
          info = ppm_error_fatal
          CALL ppm_error(ppm_err_alloc,'ppm_map_field_global',     &
     &        'send block start list PPM_MESH_ISENDBLKSTART',__LINE__,info)
          GOTO 9999
      ENDIF
      CALL ppm_alloc(ppm_mesh_isendblksize,ldu,iopt,info)
      IF (info .NE. 0) THEN
          info = ppm_error_fatal
          CALL ppm_error(ppm_err_alloc,'ppm_map_field_global',     &
     &        'send block size list PPM_MESH_ISENDBLKSIZE',__LINE__,info)
          GOTO 9999
      ENDIF

      !-------------------------------------------------------------------------
      !  Allocate memory for the global mesh receive lists
      !-------------------------------------------------------------------------
      iopt   = ppm_param_alloc_fit
      ldu(1) = nrecvlist
      CALL ppm_alloc(ppm_mesh_irecvtosub,ldu,iopt,info)
      IF (info .NE. 0) THEN
          info = ppm_error_fatal
          CALL ppm_error(ppm_err_alloc,'ppm_map_field_global',     &
     &        'destination recv sub list PPM_MESH_IRECVTOSUB',__LINE__,info)
          GOTO 9999
      ENDIF
      ldu(1) = pdim
      ldu(2) = nrecvlist
      CALL ppm_alloc(ppm_mesh_irecvblkstart,ldu,iopt,info)
      IF (info .NE. 0) THEN
          info = ppm_error_fatal
          CALL ppm_error(ppm_err_alloc,'ppm_map_field_global',     &
     &        'recv block start list PPM_MESH_IRECVBLKSTART',__LINE__,info)
          GOTO 9999
      ENDIF
      CALL ppm_alloc(ppm_mesh_irecvblksize,ldu,iopt,info)
      IF (info .NE. 0) THEN
          info = ppm_error_fatal
          CALL ppm_error(ppm_err_alloc,'ppm_map_field_global',     &
     &        'recv block size list PPM_MESH_IRECVBLKSIZE',__LINE__,info)
          GOTO 9999
      ENDIF

      !-------------------------------------------------------------------------
      !  Allocate memory for the global send/recv lists
      !-------------------------------------------------------------------------
      ldu(1) = ppm_nproc
      CALL ppm_alloc(ppm_isendlist,ldu,iopt,info)
      IF (info .NE. 0) THEN
          info = ppm_error_fatal
          CALL ppm_error(ppm_err_alloc,'ppm_map_field_global',     &
     &        'global send rank list PPM_ISENDLIST',__LINE__,info)
          GOTO 9999
      ENDIF
      CALL ppm_alloc(ppm_irecvlist,ldu,iopt,info)
      IF (info .NE. 0) THEN
          info = ppm_error_fatal
          CALL ppm_error(ppm_err_alloc,'ppm_map_field_global',     &
     &        'global recv rank list PPM_IRECVLIST',__LINE__,info)
          GOTO 9999
      ENDIF
      ldu(1) = ppm_nproc + 1
      CALL ppm_alloc(ppm_psendbuffer,ldu,iopt,info)
      IF (info .NE. 0) THEN
          info = ppm_error_fatal
          CALL ppm_error(ppm_err_alloc,'ppm_map_field_global',     &
     &        'global send buffer pointer PPM_PSENDBUFFER',__LINE__,info)
          GOTO 9999
      ENDIF
      CALL ppm_alloc(ppm_precvbuffer,ldu,iopt,info)
      IF (info .NE. 0) THEN
          info = ppm_error_fatal
          CALL ppm_error(ppm_err_alloc,'ppm_map_field_global',     &
     &        'global recv buffer pointer PPM_PRECVBUFFER',__LINE__,info)
          GOTO 9999
      ENDIF

      !-------------------------------------------------------------------------
      !  Reset the number of buffer entries
      !-------------------------------------------------------------------------
      ppm_buffer_set = 0

      !-------------------------------------------------------------------------
      !  loop over all processors; First is the processor itself
      !-------------------------------------------------------------------------
      sendrank           = ppm_rank - 1
      recvrank           = ppm_rank + 1
      ppm_psendbuffer(1) = 1
      ppm_precvbuffer(1) = 1
      ppm_nsendlist      = 0
      ppm_nrecvlist      = 0
      ppm_nsendbuffer    = 0    ! data length in buffer
      ppm_nrecvbuffer    = 0    ! data length in buffer
      iset               = 0
      ibuffer            = 0

      DO i=1,ppm_nproc
         !----------------------------------------------------------------------
         !  compute the next processor
         !----------------------------------------------------------------------
         sendrank = sendrank + 1
         IF (sendrank.GE.ppm_nproc) sendrank = sendrank - ppm_nproc 
         recvrank = recvrank - 1
         IF (recvrank.LT.        0) recvrank = recvrank + ppm_nproc 

         !----------------------------------------------------------------------
         !  Store the processor to which we will send to
         !----------------------------------------------------------------------
         ppm_nsendlist                = ppm_nsendlist + 1 
         ppm_isendlist(ppm_nsendlist) = sendrank

         !----------------------------------------------------------------------
         !  Store the processor to which we will recv from
         !----------------------------------------------------------------------
         ppm_nrecvlist                = ppm_nrecvlist + 1 
         ppm_irecvlist(ppm_nrecvlist) = recvrank
  
         !----------------------------------------------------------------------
         !  Find all mesh blocks that are sent and store them.
         !  To get all mesh blocks that are sent in round
         !  i=1,ppm_nsendlist:
         !     dest_rank = ppm_isendlist(i)
         !     ilo = ppm_psendbuffer(i)
         !     ihi = ppm_psendbuffer(i+1)-1
         !     DO j=ilo,ihi
         !       SEND(ppm_mesh_isendblk(:,j) of sub ppm_mesh_isendfromsub(j) to
         !         processor dest_rank to sub ppm_mesh_isendtosub(j))
         !     ENDDO
         !----------------------------------------------------------------------
         ibuffer = ppm_nsendlist + 1
         ppm_psendbuffer(ibuffer) = ppm_psendbuffer(ppm_nsendlist)
         DO j=1,nsendlist
             IF (ppm_subs2proc(isendtosub(j),totopo) .EQ. sendrank) THEN
                 ppm_psendbuffer(ibuffer) = ppm_psendbuffer(ibuffer) + 1
                 iset = ppm_psendbuffer(ibuffer) - 1
                 ppm_mesh_isendfromsub(iset)         = isendfromsub(j)
                 ppm_mesh_isendblkstart(1:pdim,iset) = isendblkstart(1:pdim,j)
                 ppm_mesh_isendblksize(1:pdim,iset)  = isendblksize(1:pdim,j)
                 IF (ppm_debug .GT. 1) THEN
                     IF (ppm_dim .EQ. 2) THEN
                         WRITE(mesg,'(2(A,2I4),A,I3)') ' sending ',  &
     &                   isendblkstart(1:2,j),' of size ',isendblksize(1:2,j),&
     &                   ' to ',sendrank
                     ELSEIF (ppm_dim .EQ. 3) THEN
                         WRITE(mesg,'(2(A,3I4),A,I3)') ' sending ',  &
     &                   isendblkstart(1:3,j),' of size ',isendblksize(1:3,j),&
     &                   ' to ',sendrank
                     ENDIF
                     CALL ppm_write(ppm_rank,'ppm_map_field_global',mesg,info)
                 ENDIF
             ENDIF
         ENDDO

         !----------------------------------------------------------------------
         !  Find all mesh blocks that are received and store them.
         !  To get all mesh blocks that are received in round
         !  i=1,ppm_nrecvlist:
         !     source_rank = ppm_irecvlist(i)
         !     ilo = ppm_precvbuffer(i)
         !     ihi = ppm_precvbuffer(i+1)-1
         !     DO j=ilo,ihi
         !       RECV(ppm_mesh_irecvblk(:,j) of sub ppm_mesh_irecvfromsub(j)
         !         from processor source_rank to sub ppm_mesh_irecvtosub(j))
         !     ENDDO
         !----------------------------------------------------------------------
         ibuffer = ppm_nrecvlist + 1
         ppm_precvbuffer(ibuffer) = ppm_precvbuffer(ppm_nrecvlist)
         DO j=1,nrecvlist
             IF (ppm_subs2proc(irecvfromsub(j),fromtopo) .EQ. recvrank) THEN
                 ppm_precvbuffer(ibuffer) = ppm_precvbuffer(ibuffer) + 1
                 iset = ppm_precvbuffer(ibuffer) - 1
                 ppm_mesh_irecvtosub(iset)           = irecvtosub(j)
                 ppm_mesh_irecvblkstart(1:pdim,iset) = irecvblkstart(1:pdim,j)
                 ppm_mesh_irecvblksize(1:pdim,iset)  = irecvblksize(1:pdim,j)
                 IF (ppm_debug .GT. 1) THEN
                     IF (ppm_dim .EQ. 2) THEN
                         WRITE(mesg,'(2(A,2I4),A,I3)') ' recving ',  &
     &                   irecvblkstart(1:2,j),' of size ',irecvblksize(1:2,j),&
     &                   ' from ',recvrank
                     ELSEIF (ppm_dim .EQ. 3) THEN
                         WRITE(mesg,'(2(A,3I4),A,I3)') ' recving ',  &
     &                   irecvblkstart(1:3,j),' of size ',irecvblksize(1:3,j),&
     &                   ' from ',recvrank
                     ENDIF
                     CALL ppm_write(ppm_rank,'ppm_map_field_global',mesg,info)
                 ENDIF
             ENDIF
         ENDDO
      ENDDO

      !-------------------------------------------------------------------------
      !  Deallocate memory of the local lists
      !-------------------------------------------------------------------------
      iopt   = ppm_param_dealloc
      CALL ppm_alloc(isendfromsub,ldu,iopt,info)
      IF (info .NE. 0) THEN
          info = ppm_error_error
          CALL ppm_error(ppm_err_dealloc,'ppm_map_field_global',     &
     &        'local send source sub list ISENDFROMSUB',__LINE__,info)
      ENDIF
      CALL ppm_alloc(isendtosub,ldu,iopt,info)
      IF (info .NE. 0) THEN
          info = ppm_error_error
          CALL ppm_error(ppm_err_dealloc,'ppm_map_field_global',     &
     &        'local send destination sub list ISENDTOSUB',__LINE__,info)
      ENDIF
      CALL ppm_alloc(isendblkstart,ldu,iopt,info)
      IF (info .NE. 0) THEN
          info = ppm_error_error
          CALL ppm_error(ppm_err_dealloc,'ppm_map_field_global',     &
     &        'local send block start list ISENDBLKSTART',__LINE__,info)
      ENDIF
      CALL ppm_alloc(isendblksize,ldu,iopt,info)
      IF (info .NE. 0) THEN
          info = ppm_error_error
          CALL ppm_error(ppm_err_dealloc,'ppm_map_field_global',     &
     &        'local send block size list ISENDBLKSIZE',__LINE__,info)
      ENDIF
      CALL ppm_alloc(irecvfromsub,ldu,iopt,info)
      IF (info .NE. 0) THEN
          info = ppm_error_error
          CALL ppm_error(ppm_err_dealloc,'ppm_map_field_global',     &
     &        'local recv source sub list IRECVFROMSUB',__LINE__,info)
      ENDIF
      CALL ppm_alloc(irecvtosub,ldu,iopt,info)
      IF (info .NE. 0) THEN
          info = ppm_error_error
          CALL ppm_error(ppm_err_dealloc,'ppm_map_field_global',     &
     &        'local recv destination sub list IRECVTOSUB',__LINE__,info)
      ENDIF
      CALL ppm_alloc(irecvblkstart,ldu,iopt,info)
      IF (info .NE. 0) THEN
          info = ppm_error_error
          CALL ppm_error(ppm_err_dealloc,'ppm_map_field_global',     &
     &        'local recv block start list IRECVBLKSTART',__LINE__,info)
      ENDIF
      CALL ppm_alloc(irecvblksize,ldu,iopt,info)
      IF (info .NE. 0) THEN
          info = ppm_error_error
          CALL ppm_error(ppm_err_dealloc,'ppm_map_field_global',     &
     &        'local recv block size list IRECVBLKSIZE',__LINE__,info)
      ENDIF
      CALL ppm_alloc(ioffset,ldu,iopt,info)
      IF (info .NE. 0) THEN
          info = ppm_error_error
          CALL ppm_error(ppm_err_dealloc,'ppm_map_field_global',     &
     &        'local recv block offset list IOFFSET',__LINE__,info)
      ENDIF

      !-------------------------------------------------------------------------
      !  Return 
      !-------------------------------------------------------------------------
 9999 CONTINUE
      CALL substop('ppm_map_field_global',t0,info)
      RETURN
      END SUBROUTINE ppm_map_field_global
