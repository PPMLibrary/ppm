      !-------------------------------------------------------------------------
      !  Subroutine   :              ppm_map_field_ghost_put
      !-------------------------------------------------------------------------
      !
      !  Purpose      : This routine sends the values of ghost mesh points
      !                 back to their origin in order to add their
      !                 contribution to the corresponding real mesh point.
      !                 ppm_map_field_ghost_init must have been called for
      !                 the requested topology BEFORE this routine.
      !
      !  Input        : topoid    (I) topology ID (internal numbering)
      !                 meshid    (I) mesh ID (internal numbering)
      !
      !  Input/output : 
      !
      !  Output       : info      (I)  return status. 0 on success.
      !
      !  Remarks      : The first part of the send/recv lists contains the
      !                 on-processor data.
      !
      !  References   :
      !
      !  Revisions    :
      !-------------------------------------------------------------------------
      !  $Log: ppm_map_field_ghost_put.f,v $
      !  Revision 1.1.1.1  2007/07/13 10:18:56  ivos
      !  CBL version of the PPM library
      !
      !  Revision 1.6  2004/10/01 16:33:36  ivos
      !  cosmetics.
      !
      !  Revision 1.5  2004/10/01 16:09:04  ivos
      !  Replaced REAL(ppm_kind_double) :: t0 with REAL(MK) t0.
      !
      !  Revision 1.4  2004/08/31 12:48:08  ivos
      !  Changed argument checks to use ppm_check_topoid and ppm_check_meshid.
      !
      !  Revision 1.3  2004/07/26 07:42:41  ivos
      !  Changed to single-interface modules and adapted all USE statements.
      !
      !  Revision 1.2  2004/04/05 12:03:39  ivos
      !  Udated comment header and code comments as ppm_map_type is now
      !  actually used.
      !
      !  Revision 1.1  2004/04/02 15:20:08  ivos
      !  Very preliminary check-in. Not yet tested. Handling of periodic
      !  images for periodic systems still missing.
      !
      !-------------------------------------------------------------------------
      !  Parallel Particle Mesh Library (PPM)
      !  Institute of Computational Science
      !  ETH Zentrum, Hirschengraben 84
      !  CH-8092 Zurich, Switzerland
      !-------------------------------------------------------------------------

      SUBROUTINE ppm_map_field_ghost_put(topoid,meshid,info)

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
      USE ppm_module_check_topoid
      USE ppm_module_check_meshid
      IMPLICIT NONE
      !-------------------------------------------------------------------------
      !  Arguments     
      !-------------------------------------------------------------------------
      INTEGER                 , INTENT(IN   ) :: topoid,meshid
      INTEGER                 , INTENT(  OUT) :: info
      !-------------------------------------------------------------------------
      !  Local variables 
      !-------------------------------------------------------------------------
      INTEGER, DIMENSION(3)            :: ldu
      INTEGER                          :: i,j,lb,ub
      INTEGER                          :: iopt,ibuffer,pdim
      INTEGER                          :: nsendlist,nrecvlist
      REAL(ppm_kind_double)            :: t0
      LOGICAL                          :: valid
      !-------------------------------------------------------------------------
      !  Externals 
      !-------------------------------------------------------------------------
      
      !-------------------------------------------------------------------------
      !  Initialise 
      !-------------------------------------------------------------------------
      CALL substart('ppm_map_field_ghost_put',t0,info)
      pdim = ppm_dim

      !-------------------------------------------------------------------------
      !  Check arguments
      !-------------------------------------------------------------------------
      IF (ppm_debug .GT. 0) THEN
          CALL ppm_check_topoid(ppm_param_id_internal,topoid,valid,info)
          IF (.NOT. valid) THEN
              info = ppm_error_error
              CALL ppm_error(ppm_err_argument,'ppm_map_field_ghost_put',  &
     &            'topoid out of range',__LINE__,info)
              GOTO 9999
          ENDIF
          CALL ppm_check_meshid(ppm_param_id_internal,meshid,topoid,valid,info)
          IF (.NOT. valid) THEN
              info = ppm_error_error
              CALL ppm_error(ppm_err_argument,'ppm_map_field_ghost_put',  &
     &            'meshid out of range',__LINE__,info)
              GOTO 9999
          ENDIF
      ENDIF

      !-------------------------------------------------------------------------
      !  Save the map type for the subsequent calls (used to set pop type)
      !-------------------------------------------------------------------------
      ppm_map_type = ppm_param_map_ghost_put

      !-------------------------------------------------------------------------
      !  Reset the buffer size counters
      !-------------------------------------------------------------------------
      ppm_nsendbuffer = 0
      ppm_nrecvbuffer = 0
      
      !-------------------------------------------------------------------------
      !  Number of mesh blocks to be sent/recvd is reverse from the get
      !-------------------------------------------------------------------------
      nrecvlist = ppm_mesh_ghost_nsend(meshid,topoid)
      nsendlist = ppm_mesh_ghost_nrecv(meshid,topoid)

      !-------------------------------------------------------------------------
      !  Allocate memory for the global mesh sendlists
      !-------------------------------------------------------------------------
      iopt   = ppm_param_alloc_fit
      ldu(1) = nsendlist
      CALL ppm_alloc(ppm_mesh_isendfromsub,ldu,iopt,info)
      IF (info .NE. 0) THEN
          info = ppm_error_fatal
          CALL ppm_error(ppm_err_alloc,'ppm_map_field_ghost_put',     &
     &        'source send sub list PPM_MESH_ISENDFROMSUB',__LINE__,info)
          GOTO 9999
      ENDIF
      ldu(1) = pdim
      ldu(2) = nsendlist
      CALL ppm_alloc(ppm_mesh_isendblkstart,ldu,iopt,info)
      IF (info .NE. 0) THEN
          info = ppm_error_fatal
          CALL ppm_error(ppm_err_alloc,'ppm_map_field_ghost_put',     &
     &        'send block start list PPM_MESH_ISENDBLKSTART',__LINE__,info)
          GOTO 9999
      ENDIF
      CALL ppm_alloc(ppm_mesh_isendblksize,ldu,iopt,info)
      IF (info .NE. 0) THEN
          info = ppm_error_fatal
          CALL ppm_error(ppm_err_alloc,'ppm_map_field_ghost_put',     &
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
          CALL ppm_error(ppm_err_alloc,'ppm_map_field_ghost_put',     &
     &        'destination recv sub list PPM_MESH_IRECVTOSUB',__LINE__,info)
          GOTO 9999
      ENDIF
      ldu(1) = pdim
      ldu(2) = nrecvlist
      CALL ppm_alloc(ppm_mesh_irecvblkstart,ldu,iopt,info)
      IF (info .NE. 0) THEN
          info = ppm_error_fatal
          CALL ppm_error(ppm_err_alloc,'ppm_map_field_ghost_put',     &
     &        'recv block start list PPM_MESH_IRECVBLKSTART',__LINE__,info)
          GOTO 9999
      ENDIF
      CALL ppm_alloc(ppm_mesh_irecvblksize,ldu,iopt,info)
      IF (info .NE. 0) THEN
          info = ppm_error_fatal
          CALL ppm_error(ppm_err_alloc,'ppm_map_field_ghost_put',     &
     &        'recv block size list PPM_MESH_IRECVBLKSIZE',__LINE__,info)
          GOTO 9999
      ENDIF

      !-------------------------------------------------------------------------
      !  Allocate memory for the global send/recv lists
      !-------------------------------------------------------------------------
      ldu(1) = ppm_ncommseq(topoid)
      CALL ppm_alloc(ppm_isendlist,ldu,iopt,info)
      IF (info .NE. 0) THEN
          info = ppm_error_fatal
          CALL ppm_error(ppm_err_alloc,'ppm_map_field_ghost_put',     &
     &        'global send rank list PPM_ISENDLIST',__LINE__,info)
          GOTO 9999
      ENDIF
      CALL ppm_alloc(ppm_irecvlist,ldu,iopt,info)
      IF (info .NE. 0) THEN
          info = ppm_error_fatal
          CALL ppm_error(ppm_err_alloc,'ppm_map_field_ghost_put',     &
     &        'global recv rank list PPM_IRECVLIST',__LINE__,info)
          GOTO 9999
      ENDIF
      ldu(1) = ppm_ncommseq(topoid) + 1
      CALL ppm_alloc(ppm_psendbuffer,ldu,iopt,info)
      IF (info .NE. 0) THEN
          info = ppm_error_fatal
          CALL ppm_error(ppm_err_alloc,'ppm_map_field_ghost_put',     &
     &        'global send buffer pointer PPM_PSENDBUFFER',__LINE__,info)
          GOTO 9999
      ENDIF
      CALL ppm_alloc(ppm_precvbuffer,ldu,iopt,info)
      IF (info .NE. 0) THEN
          info = ppm_error_fatal
          CALL ppm_error(ppm_err_alloc,'ppm_map_field_ghost_put',     &
     &        'global recv buffer pointer PPM_PRECVBUFFER',__LINE__,info)
          GOTO 9999
      ENDIF

      !-------------------------------------------------------------------------
      !  Reset the number of buffer entries
      !-------------------------------------------------------------------------
      ppm_buffer_set = 0

      !-------------------------------------------------------------------------
      !  Build the lists for ppm_map_field_push, _pop and _send to CREATE ghost
      !  meshpoints on the topology topoid
      !-------------------------------------------------------------------------
      ppm_psendbuffer(1) = 1
      ppm_precvbuffer(1) = 1
      ppm_nsendlist      = ppm_ncommseq(topoid)
      ppm_nrecvlist      = ppm_nsendlist
      ibuffer            = 0
      DO i=1,ppm_ncommseq(topoid)
         ibuffer = i + 1
         !----------------------------------------------------------------------
         !  Store the processor with which we will send/recv
         !----------------------------------------------------------------------
         ppm_isendlist(i) = ppm_icommseq(i,topoid)
         ppm_irecvlist(i) = ppm_icommseq(i,topoid)

         !----------------------------------------------------------------------
         !  Store the mesh block range for this processor interaction
         !----------------------------------------------------------------------
         ppm_psendbuffer(ibuffer)=ppm_mesh_ghost_recvblk(ibuffer,meshid,topoid)
         ppm_precvbuffer(ibuffer)=ppm_mesh_ghost_blk(ibuffer,meshid,topoid)
         
         !----------------------------------------------------------------------
         !  Store the definitions of the mesh blocks to be sent
         !----------------------------------------------------------------------
         lb = ppm_psendbuffer(i)
         ub = ppm_psendbuffer(ibuffer)
         DO j=lb,ub-1
             ppm_mesh_isendblkstart(1:pdim,j) =     &
     &           ppm_mesh_ghost_recvblkstart(1:pdim,j,meshid,topoid)
             ppm_mesh_isendblksize(1:pdim,j) =      &
     &           ppm_mesh_ghost_recvblksize(1:pdim,j,meshid,topoid)
             ppm_mesh_isendfromsub(j)=ppm_mesh_ghost_recvtosub(j,meshid,topoid)
         ENDDO
         
         !----------------------------------------------------------------------
         !  Store the definitions of the mesh blocks to be received
         !----------------------------------------------------------------------
         lb = ppm_precvbuffer(i)
         ub = ppm_precvbuffer(ibuffer)
         DO j=lb,ub-1
             ppm_mesh_irecvblkstart(1:pdim,j) =     &
     &           ppm_mesh_ghost_blkstart(1:pdim,j,meshid,topoid)
             ppm_mesh_irecvblksize(1:pdim,j) =      &
     &           ppm_mesh_ghost_blksize(1:pdim,j,meshid,topoid)
             ppm_mesh_irecvtosub(j) = ppm_mesh_ghost_fromsub(j,meshid,topoid)
         ENDDO
      ENDDO

      !-------------------------------------------------------------------------
      !  Return 
      !-------------------------------------------------------------------------
 9999 CONTINUE
      CALL substop('ppm_map_field_ghost_put',t0,info)
      RETURN
      END SUBROUTINE ppm_map_field_ghost_put
