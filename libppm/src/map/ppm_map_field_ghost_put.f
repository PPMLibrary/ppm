      !-------------------------------------------------------------------------
      !  Subroutine   :              ppm_map_field_ghost_put
      !-------------------------------------------------------------------------
      !  Parallel Particle Mesh Library (PPM)
      !  ETH Zurich
      !  CH-8092 Zurich, Switzerland
      !-------------------------------------------------------------------------

      SUBROUTINE ppm_map_field_ghost_put(topoid,meshid,ghostsize,info)
      !!! This routine sends the values of ghost mesh points
      !!! back to their origin in order to add their contribution to the
      !!! corresponding real mesh point.
      !!!
      !!! [IMPORTANT]
      !!! `ppm_map_field_ghost_init` must not be called anymore. This routine
      !!! checks wheterh the ghost mappings have been initialized already and
      !!! if not, the routine is called by the library.
      !!!
      !!! [NOTE]
      !!! The first part of the send/recv lists contains the on-processor data.

      !-------------------------------------------------------------------------
      !  Includes
      !-------------------------------------------------------------------------
#include "ppm_define.h"

      !-------------------------------------------------------------------------
      !  Modules
      !-------------------------------------------------------------------------
      USE ppm_module_data
      USE ppm_module_data_mesh
      USE ppm_module_typedef
      USE ppm_module_substart
      USE ppm_module_substop
      USE ppm_module_error
      USE ppm_module_alloc
      USE ppm_module_check_id
      IMPLICIT NONE
      !-------------------------------------------------------------------------
      !  Arguments
      !-------------------------------------------------------------------------
      INTEGER                 , INTENT(IN   ) :: topoid
      !!! Topology ID
      INTEGER                 , INTENT(IN   ) :: meshid
      !!! Mesh ID
      INTEGER, DIMENSION(:)   , INTENT(IN   ) :: ghostsize
      !!! Size of the ghost layer in numbers of grid points in all space
      !!! dimensions (1...ppm_dim).
      INTEGER                 , INTENT(  OUT) :: info
      !!! Returns status, 0 upon success
      !-------------------------------------------------------------------------
      !  Local variables
      !-------------------------------------------------------------------------
      INTEGER, DIMENSION(3)            :: ldu
      INTEGER                          :: i,j,lb,ub
      INTEGER                          :: iopt,ibuffer,pdim
      INTEGER                          :: nsendlist,nrecvlist
      REAL(ppm_kind_double)            :: t0
      LOGICAL                          :: valid
      TYPE(ppm_t_topo), POINTER        :: topo
      TYPE(ppm_t_equi_mesh), POINTER   :: mesh
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
        CALL check
        IF (info .NE. 0) GOTO 9999
      ENDIF

      topo => ppm_topo(topoid)%t
      mesh => topo%mesh(meshid)


      !-------------------------------------------------------------------------
      !  Check if ghost mappings have been initalized, if no do so now.
      !-------------------------------------------------------------------------
      IF (.NOT. mesh%ghost_initialized) THEN
        CALL ppm_map_field_ghost_init(topoid,meshid,ghostsize,info)
        IF (info .NE. 0) THEN
            info = ppm_error_error
            CALL ppm_error(ppm_err_argument,'ppm_map_field_ghost_put',  &
     &          'Could not initialize ghost mappings',__LINE__,info)
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
      nrecvlist = mesh%ghost_nsend
      nsendlist = mesh%ghost_nrecv

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
      ldu(1) = topo%ncommseq
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
      ldu(1) = topo%ncommseq + 1
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
      ppm_nsendlist      = topo%ncommseq
      ppm_nrecvlist      = ppm_nsendlist
      ibuffer            = 0
      DO i=1,topo%ncommseq
         ibuffer = i + 1
         !----------------------------------------------------------------------
         !  Store the processor with which we will send/recv
         !----------------------------------------------------------------------
         ppm_isendlist(i) = topo%icommseq(i)
         ppm_irecvlist(i) = topo%icommseq(i)

         !----------------------------------------------------------------------
         !  Store the mesh block range for this processor interaction
         !----------------------------------------------------------------------
         ppm_psendbuffer(ibuffer)=mesh%ghost_recvblk(ibuffer)
         ppm_precvbuffer(ibuffer)=mesh%ghost_blk(ibuffer)

         !----------------------------------------------------------------------
         !  Store the definitions of the mesh blocks to be sent
         !----------------------------------------------------------------------
         lb = ppm_psendbuffer(i)
         ub = ppm_psendbuffer(ibuffer)
         DO j=lb,ub-1
             ppm_mesh_isendblkstart(1:pdim,j) =     &
     &           mesh%ghost_recvblkstart(1:pdim,j)
             ppm_mesh_isendblksize(1:pdim,j) =      &
     &           mesh%ghost_recvblksize(1:pdim,j)
             ppm_mesh_isendfromsub(j)=mesh%ghost_recvtosub(j)
         ENDDO

         !----------------------------------------------------------------------
         !  Store the definitions of the mesh blocks to be received
         !----------------------------------------------------------------------
         lb = ppm_precvbuffer(i)
         ub = ppm_precvbuffer(ibuffer)
         DO j=lb,ub-1
             ppm_mesh_irecvblkstart(1:pdim,j) =     &
     &           mesh%ghost_blkstart(1:pdim,j)
             ppm_mesh_irecvblksize(1:pdim,j) =      &
     &           mesh%ghost_blksize(1:pdim,j)
             ppm_mesh_irecvtosub(j) = mesh%ghost_fromsub(j)
         ENDDO
      ENDDO

      !-------------------------------------------------------------------------
      !  Return
      !-------------------------------------------------------------------------
 9999 CONTINUE
      CALL substop('ppm_map_field_ghost_put',t0,info)
      RETURN
      CONTAINS
      SUBROUTINE check
          CALL ppm_check_topoid(topoid,valid,info)
          IF (.NOT. valid) THEN
              info = ppm_error_error
              CALL ppm_error(ppm_err_argument,'ppm_map_field_ghost_put',  &
     &            'topoid out of range',__LINE__,info)
              GOTO 8888
          ENDIF
          CALL ppm_check_meshid(topoid,meshid,valid,info)
          IF (.NOT. valid) THEN
              info = ppm_error_error
              CALL ppm_error(ppm_err_argument,'ppm_map_field_ghost_put',  &
     &            'meshid out of range',__LINE__,info)
              GOTO 8888
          ENDIF
 8888     CONTINUE
      END SUBROUTINE check
      END SUBROUTINE ppm_map_field_ghost_put
