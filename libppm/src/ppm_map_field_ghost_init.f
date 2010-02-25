      !-------------------------------------------------------------------------
      !  Subroutine   :               ppm_map_field_ghost_init
      !-------------------------------------------------------------------------
      !
      !  Purpose      : This routine sets up the lists needed for sending
      !                 and receiving ghost mesh points on a given
      !                 topology/mesh combination. These lists are stored
      !                 for later use by ppm_map_field_ghost_get and
      !                 ppm_map_field_ghost_put. They need to be rebuilt as
      !                 soon as a new mesh or topology is defined,
      !                 replacing the existing one.
      !
      !  Input        : topoid       (I) topology identifier (internal numbering)
      !                 meshid       (I) mesh identifier (internal numbering)
      !                 ghostsize(:) (I) size of the ghost layer in numbers
      !                                  of grid points in all space
      !                                  dimensions (1...ppm_dim). 
      !
      !  Input/output : 
      !
      !  Output       : info      (I)  return status. 0 on success.
      !
      !  Remarks      : The first part of the send/recv lists contains the
      !                 on-processor data.
      !
      !                 The names send and recv are correct for the ghost
      !                 get map. For the put case, they need to be
      !                 interpreted the other way round.
      !
      !                 The shift cases for periodic images are currently hard
      !                 coded. Probably there is a more elegant way of doing
      !                 this.
      !
      !  References   :
      !
      !  Revisions    :
      !-------------------------------------------------------------------------
      !  $Log: ppm_map_field_ghost_init.f,v $
      !  Revision 1.1.1.1  2007/07/13 10:18:55  ivos
      !  CBL version of the PPM library
      !
      !  Revision 1.15  2006/02/03 09:34:03  ivos
      !  Fixed bug 00015: ppm_subs_bc was only allocated and stored for the
      !  local subs in topo_store. Several mapping routines however need the
      !  info about all (global) subs.
      !  Changed subs_bc to hold now the GLOBAL subid and adjusted all
      !  occurrences.
      !
      !  Revision 1.14  2005/05/24 23:27:44  ivos
      !  Routine checked once more. Works properly. Some comments corrected.
      !
      !  Revision 1.13  2004/11/11 15:26:15  ivos
      !  Moved allocatable work arrays to module.
      !
      !  Revision 1.12  2004/10/01 16:33:36  ivos
      !  cosmetics.
      !
      !  Revision 1.11  2004/10/01 16:09:04  ivos
      !  Replaced REAL(ppm_kind_double) :: t0 with REAL(MK) t0.
      !
      !  Revision 1.10  2004/08/31 12:48:08  ivos
      !  Changed argument checks to use ppm_check_topoid and ppm_check_meshid.
      !
      !  Revision 1.9  2004/07/26 15:38:48  ivos
      !  Inserted missing USE statements to resolve undefined references
      !  at link stage.
      !
      !  Revision 1.8  2004/07/26 07:42:41  ivos
      !  Changed to single-interface modules and adapted all USE statements.
      !
      !  Revision 1.7  2004/04/13 12:13:24  ivos
      !  Forgot corners (and edges in 3d) in systems with multiple periodic
      !  directions. Added the corresponding shift cases now. Tested in 2d.
      !
      !  Revision 1.6  2004/04/08 11:08:58  ivos
      !  bugfix: offset list was not sorted according to icommseq. This
      !  caused problems on more than 2 processors.
      !
      !  Revision 1.5  2004/04/07 15:34:15  ivos
      !  Changed to use ppm_mesh_block_intersect. Tested.
      !  Added handling of periodic images. Not tested yet.
      !
      !  Revision 1.4  2004/04/07 09:18:33  ivos
      !  bugfix: added 2 times the ghostsize to ndata.
      !
      !  Revision 1.3  2004/04/05 12:03:39  ivos
      !  Udated comment header and code comments as ppm_map_type is now
      !  actually used.
      !
      !  Revision 1.2  2004/04/05 08:03:50  ivos
      !  bugfix: added missing indices in ppm_mesh_ghost_blk.
      !
      !  Revision 1.1  2004/04/02 15:20:07  ivos
      !  Very preliminary check-in. Not yet tested. Handling of periodic
      !  images for periodic systems still missing.
      !
      !-------------------------------------------------------------------------
      !  Parallel Particle Mesh Library (PPM)
      !  Institute of Computational Science
      !  ETH Zentrum, Hirschengraben 84
      !  CH-8092 Zurich, Switzerland
      !-------------------------------------------------------------------------

      SUBROUTINE ppm_map_field_ghost_init(topoid,meshid,ghostsize,info)

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
      !  Includes
      !-------------------------------------------------------------------------
#include "ppm_define.h"
#ifdef __MPI
      INCLUDE 'mpif.h'
#endif
      !-------------------------------------------------------------------------
      !  Arguments     
      !-------------------------------------------------------------------------
      INTEGER                 , INTENT(IN   ) :: topoid,meshid
      INTEGER, DIMENSION(:)   , INTENT(IN   ) :: ghostsize
      INTEGER                 , INTENT(  OUT) :: info
      !-------------------------------------------------------------------------
      !  Local variables 
      !-------------------------------------------------------------------------
      INTEGER, DIMENSION(4)            :: ldu
      INTEGER, DIMENSION(ppm_dim)      :: op
      INTEGER, DIMENSION(ppm_dim,26)   :: ond
      INTEGER                          :: i,j,sendrank,recvrank,idom,jdom,k
      INTEGER                          :: iopt,iset,ibuffer,pdim,isize,nnd
      INTEGER                          :: nsendlist,nsend,tag1,lb,ub,nrecv
      CHARACTER(ppm_char)              :: mesg
      REAL(ppm_kind_double)            :: t0
      LOGICAL                          :: lsouth,lnorth,least,lwest,ltop,lbottom
#ifdef __MPI
      INTEGER, DIMENSION(MPI_STATUS_SIZE) :: commstat
#endif
      !-------------------------------------------------------------------------
      !  Externals 
      !-------------------------------------------------------------------------
      
      !-------------------------------------------------------------------------
      !  Initialise 
      !-------------------------------------------------------------------------
      CALL substart('ppm_map_field_ghost_init',t0,info)
      pdim = ppm_dim

      !-------------------------------------------------------------------------
      !  Check arguments
      !-------------------------------------------------------------------------
      IF (ppm_debug .GT. 0) THEN
          CALL ppm_check_topoid(ppm_param_id_internal,topoid,ltop,info)
          IF (.NOT. ltop) THEN
              info = ppm_error_error
              CALL ppm_error(ppm_err_argument,'ppm_map_field_ghost_init',  &
     &            'topoid out of range',__LINE__,info)
              GOTO 9999
          ENDIF
          CALL ppm_check_meshid(ppm_param_id_internal,meshid,topoid,ltop,info)
          IF (.NOT. ltop) THEN
              info = ppm_error_error
              CALL ppm_error(ppm_err_argument,'ppm_map_field_ghost_init',  &
     &            'meshid out of range',__LINE__,info)
              GOTO 9999
          ENDIF
          IF (SIZE(ghostsize,1) .LT. ppm_dim) THEN
              info = ppm_error_error
              CALL ppm_error(ppm_err_argument,'ppm_map_field_ghost_init',  &
     &            'ghostsize must be at least of length ppm_dim',__LINE__,info)
              GOTO 9999
          ENDIF
      ENDIF

      !-------------------------------------------------------------------------
      !  Save the map type for the subsequent calls (used to check push/pop)
      !-------------------------------------------------------------------------
      ppm_map_type = ppm_param_map_init

      !-------------------------------------------------------------------------
      !  Allocate memory for the local temporary sendlists
      !-------------------------------------------------------------------------
      isize  = ppm_nsublist(topoid)
      iopt   = ppm_param_alloc_fit
      ldu(1) = isize
      CALL ppm_alloc(isendfromsub,ldu,iopt,info)
      IF (info .NE. 0) THEN
          info = ppm_error_fatal
          CALL ppm_error(ppm_err_alloc,'ppm_map_field_ghost_init',     &
     &        'local send source sub list ISENDFROMSUB',__LINE__,info)
          GOTO 9999
      ENDIF
      CALL ppm_alloc(isendtosub,ldu,iopt,info)
      IF (info .NE. 0) THEN
          info = ppm_error_fatal
          CALL ppm_error(ppm_err_alloc,'ppm_map_field_ghost_init',     &
     &        'local send destination sub list ISENDTOSUB',__LINE__,info)
          GOTO 9999
      ENDIF
      ldu(1) = pdim
      ldu(2) = isize
      CALL ppm_alloc(isendblkstart,ldu,iopt,info)
      IF (info .NE. 0) THEN
          info = ppm_error_fatal
          CALL ppm_error(ppm_err_alloc,'ppm_map_field_ghost_init',     &
     &        'local send block start list ISENDBLKSTART',__LINE__,info)
          GOTO 9999
      ENDIF
      CALL ppm_alloc(isendblksize,ldu,iopt,info)
      IF (info .NE. 0) THEN
          info = ppm_error_fatal
          CALL ppm_error(ppm_err_alloc,'ppm_map_field_ghost_init',     &
     &        'local send block size list ISENDBLKSIZE',__LINE__,info)
          GOTO 9999
      ENDIF
      CALL ppm_alloc(ioffset,ldu,iopt,info)
      IF (info .NE. 0) THEN
          info = ppm_error_fatal
          CALL ppm_error(ppm_err_alloc,'ppm_map_field_ghost_init',     &
     &        'local mesh block offset list IOFFSET',__LINE__,info)
          GOTO 9999
      ENDIF

      !-------------------------------------------------------------------------
      !  Pre-calculate shift offsets for periodic ghost images of subs
      !-------------------------------------------------------------------------
      op(1) = ppm_cart_mesh(meshid,topoid)%Nm(1) - 1
      op(2) = ppm_cart_mesh(meshid,topoid)%Nm(2) - 1
      IF (pdim .GT. 2) THEN
          op(3) = ppm_cart_mesh(meshid,topoid)%Nm(3) - 1
      ENDIF

      !-------------------------------------------------------------------------
      !  Find intersecting mesh domains to be sent for the ghost get
      !-------------------------------------------------------------------------
      nsendlist = 0
      DO i=1,ppm_nsublist(topoid)
          idom = ppm_isublist(i,topoid)
          ond(1:pdim,1:26) = 0
          !---------------------------------------------------------------------
          !  First, we check the neighbors of the subdomain
          !---------------------------------------------------------------------
          DO j=1,ppm_nneighsubs(i,topoid)
              jdom = ppm_ineighsubs(j,i,topoid)
              ! source and destination meshes and topologies are identical
              CALL ppm_mesh_block_intersect(idom,jdom,meshid,meshid,topoid, &
     &            topoid,ond(1:pdim,1),ghostsize,nsendlist,isendfromsub,   &
     &            isendtosub,isendblkstart,isendblksize,ioffset,info)
              IF (info .NE. 0) GOTO 9999
          ENDDO

          !---------------------------------------------------------------------
          !  In periodic systems, the box needs to be shifted and overlaps
          !  recomputed in order to get the periodic image neighbors.
          !  Check if any face of this sub coincides with a periodic domain
          !  boundary.
          !---------------------------------------------------------------------
          lwest  = .FALSE.
          IF ((ppm_bcdef(1,topoid) .EQ. ppm_param_bcdef_periodic) .AND. &
     &        (ppm_subs_bc(1,idom,topoid) .NE. 0)) lwest = .TRUE.
          least  = .FALSE.
          IF ((ppm_bcdef(2,topoid) .EQ. ppm_param_bcdef_periodic) .AND. &
     &        (ppm_subs_bc(2,idom,topoid) .NE. 0)) least = .TRUE.
          lsouth = .FALSE.
          IF ((ppm_bcdef(3,topoid) .EQ. ppm_param_bcdef_periodic) .AND. &
     &        (ppm_subs_bc(3,idom,topoid) .NE. 0)) lsouth = .TRUE.
          lnorth = .FALSE.
          IF ((ppm_bcdef(4,topoid) .EQ. ppm_param_bcdef_periodic) .AND. &
     &        (ppm_subs_bc(4,idom,topoid) .NE. 0)) lnorth = .TRUE.
          IF (pdim .GT. 2) THEN
              lbottom= .FALSE.
              IF ((ppm_bcdef(5,topoid) .EQ. ppm_param_bcdef_periodic) .AND. &
     &            (ppm_subs_bc(5,idom,topoid) .NE. 0)) lbottom = .TRUE.
              ltop   = .FALSE.
              IF ((ppm_bcdef(6,topoid) .EQ. ppm_param_bcdef_periodic) .AND. &
     &            (ppm_subs_bc(6,idom,topoid) .NE. 0)) ltop = .TRUE.
          ENDIF

          !---------------------------------------------------------------------
          !  Determine number of shifts and actual shift indices needed
          !---------------------------------------------------------------------
          nnd = 0
          IF (lwest) THEN
              nnd = nnd + 1
              ond(1,nnd) = op(1)
          ENDIF
          IF (least) THEN
              nnd = nnd + 1
              ond(1,nnd) = -op(1)
          ENDIF
          IF (lsouth) THEN
              nnd = nnd + 1
              ond(2,nnd) = op(2)
              IF (lwest) THEN
                  nnd = nnd + 1
                  ond(1,nnd) = op(1)
                  ond(2,nnd) = op(2)
              ENDIF
              IF (least) THEN
                  nnd = nnd + 1
                  ond(1,nnd) = -op(1)
                  ond(2,nnd) = op(2)
              ENDIF
          ENDIF
          IF (lnorth) THEN
              nnd = nnd + 1
              ond(2,nnd) = -op(2)
              IF (lwest) THEN
                  nnd = nnd + 1
                  ond(1,nnd) = op(1)
                  ond(2,nnd) = -op(2)
              ENDIF
              IF (least) THEN
                  nnd = nnd + 1
                  ond(1,nnd) = -op(1)
                  ond(2,nnd) = -op(2)
              ENDIF
          ENDIF
          IF (pdim .GT. 2) THEN
              IF (lbottom) THEN
                  nnd = nnd + 1
                  ond(3,nnd) = op(3)
                  IF (lwest) THEN
                      nnd = nnd + 1
                      ond(1,nnd) = op(1)
                      ond(3,nnd) = op(3)
                  ENDIF
                  IF (least) THEN
                      nnd = nnd + 1
                      ond(1,nnd) = -op(1)
                      ond(3,nnd) = op(3)
                  ENDIF
                  IF (lsouth) THEN
                      nnd = nnd + 1
                      ond(2,nnd) = op(2)
                      ond(3,nnd) = op(3)
                      IF (lwest) THEN
                          nnd = nnd + 1
                          ond(1,nnd) = op(1)
                          ond(2,nnd) = op(2)
                          ond(3,nnd) = op(3)
                      ENDIF
                      IF (least) THEN
                          nnd = nnd + 1
                          ond(1,nnd) = -op(1)
                          ond(2,nnd) = op(2)
                          ond(3,nnd) = op(3)
                      ENDIF
                  ENDIF
                  IF (lnorth) THEN
                      nnd = nnd + 1
                      ond(2,nnd) = -op(2)
                      ond(3,nnd) = op(3)
                      IF (lwest) THEN
                          nnd = nnd + 1
                          ond(1,nnd) = op(1)
                          ond(2,nnd) = -op(2)
                          ond(3,nnd) = op(3)
                      ENDIF
                      IF (least) THEN
                          nnd = nnd + 1
                          ond(1,nnd) = -op(1)
                          ond(2,nnd) = -op(2)
                          ond(3,nnd) = op(3)
                      ENDIF
                  ENDIF
              ENDIF
              IF (ltop) THEN
                  nnd = nnd + 1
                  ond(3,nnd) = -op(3)
                  IF (lwest) THEN
                      nnd = nnd + 1
                      ond(1,nnd) = op(1)
                      ond(3,nnd) = -op(3)
                  ENDIF
                  IF (least) THEN
                      nnd = nnd + 1
                      ond(1,nnd) = -op(1)
                      ond(3,nnd) = -op(3)
                  ENDIF
                  IF (lsouth) THEN
                      nnd = nnd + 1
                      ond(2,nnd) = op(2)
                      ond(3,nnd) = -op(3)
                      IF (lwest) THEN
                          nnd = nnd + 1
                          ond(1,nnd) = op(1)
                          ond(2,nnd) = op(2)
                          ond(3,nnd) = -op(3)
                      ENDIF
                      IF (least) THEN
                          nnd = nnd + 1
                          ond(1,nnd) = -op(1)
                          ond(2,nnd) = op(2)
                          ond(3,nnd) = -op(3)
                      ENDIF
                  ENDIF
                  IF (lnorth) THEN
                      nnd = nnd + 1
                      ond(2,nnd) = -op(2)
                      ond(3,nnd) = -op(3)
                      IF (lwest) THEN
                          nnd = nnd + 1
                          ond(1,nnd) = op(1)
                          ond(2,nnd) = -op(2)
                          ond(3,nnd) = -op(3)
                      ENDIF
                      IF (least) THEN
                          nnd = nnd + 1
                          ond(1,nnd) = -op(1)
                          ond(2,nnd) = -op(2)
                          ond(3,nnd) = -op(3)
                      ENDIF
                  ENDIF
              ENDIF
          ENDIF

          !---------------------------------------------------------------------
          !  Do the shifting and compute periodic ghost mesh blocks
          !---------------------------------------------------------------------
          DO k=1,nnd
              ! first with the original (non-shifted) image of itself
              jdom = idom
              CALL ppm_mesh_block_intersect(idom,jdom,meshid,meshid,topoid, &
     &            topoid,ond(1:pdim,k),ghostsize,nsendlist,isendfromsub,  &
     &            isendtosub,isendblkstart,isendblksize,ioffset,info)
              IF (info .NE. 0) GOTO 9999
              ! Then with all the neighbors
              DO j=1,ppm_nneighsubs(i,topoid)
                  jdom = ppm_ineighsubs(j,i,topoid)
                  CALL ppm_mesh_block_intersect(idom,jdom,meshid,meshid,   &
     &                topoid,topoid,ond(1:pdim,k),ghostsize,nsendlist,   &
     &                isendfromsub,isendtosub,isendblkstart,isendblksize, &
     &                ioffset,info)
                  IF (info .NE. 0) GOTO 9999
              ENDDO
          ENDDO
      ENDDO    ! i=1,ppm_nsublist

      !-------------------------------------------------------------------------
      !  Grow memory for the global mesh ghost maps
      !-------------------------------------------------------------------------
      iopt   = ppm_param_alloc_grow_preserve
      ldu(1) = MAXVAL(ppm_max_meshid)
      ldu(2) = ppm_max_topoid
      CALL ppm_alloc(ppm_mesh_ghost_nsend,ldu,iopt,info)
      IF (info .NE. 0) THEN
          info = ppm_error_fatal
          CALL ppm_error(ppm_err_alloc,'ppm_map_field_ghost_init',     &
     &        'number of blocks to send PPM_MESH_GHOST_NSEND',__LINE__,info)
          GOTO 9999
      ENDIF
      CALL ppm_alloc(ppm_mesh_ghost_nrecv,ldu,iopt,info)
      IF (info .NE. 0) THEN
          info = ppm_error_fatal
          CALL ppm_error(ppm_err_alloc,'ppm_map_field_ghost_init',     &
     &        'number of blocks to recv PPM_MESH_GHOST_NRECV',__LINE__,info)
          GOTO 9999
      ENDIF
      ldu(1) = nsendlist
      ldu(2) = MAXVAL(ppm_max_meshid)
      ldu(3) = ppm_max_topoid
      CALL ppm_alloc(ppm_mesh_ghost_fromsub,ldu,iopt,info)
      IF (info .NE. 0) THEN
          info = ppm_error_fatal
          CALL ppm_error(ppm_err_alloc,'ppm_map_field_ghost_init',     &
     &        'ghost source subs PPM_MESH_GHOST_FROMSUB',__LINE__,info)
          GOTO 9999
      ENDIF
      !CALL ppm_alloc(ppm_mesh_ghost_recvfromsub,ldu,iopt,info)
      !IF (info .NE. 0) THEN
      !    info = ppm_error_fatal
      !    CALL ppm_error(ppm_err_alloc,'ppm_map_field_ghost_init',     &
   !  &        'ghost recv source subs PPM_MESH_GHOST_RECVFROMSUB',__LINE__,info)
   !       GOTO 9999
   !   ENDIF
      CALL ppm_alloc(ppm_mesh_ghost_tosub,ldu,iopt,info)
      IF (info .NE. 0) THEN
          info = ppm_error_fatal
          CALL ppm_error(ppm_err_alloc,'ppm_map_field_ghost_init',     &
     &        'ghost target subs PPM_MESH_GHOST_TOSUB',__LINE__,info)
          GOTO 9999
      ENDIF
      CALL ppm_alloc(ppm_mesh_ghost_recvtosub,ldu,iopt,info)
      IF (info .NE. 0) THEN
          info = ppm_error_fatal
          CALL ppm_error(ppm_err_alloc,'ppm_map_field_ghost_init',     &
     &        'ghost recv target subs PPM_MESH_GHOST_RECVTOSUB',__LINE__,info)
          GOTO 9999
      ENDIF
      ldu(1) = ppm_ncommseq(topoid) + 1
      CALL ppm_alloc(ppm_mesh_ghost_blk,ldu,iopt,info)
      IF (info .NE. 0) THEN
          info = ppm_error_fatal
          CALL ppm_error(ppm_err_alloc,'ppm_map_field_ghost_init',     &
     &        'ghost block list PPM_MESH_GHOST_BLK',__LINE__,info)
          GOTO 9999
      ENDIF
      CALL ppm_alloc(ppm_mesh_ghost_recvblk,ldu,iopt,info)
      IF (info .NE. 0) THEN
          info = ppm_error_fatal
          CALL ppm_error(ppm_err_alloc,'ppm_map_field_ghost_init',     &
     &        'ghost recv block list PPM_MESH_GHOST_RECVBLK',__LINE__,info)
          GOTO 9999
      ENDIF
      ldu(1) = pdim
      ldu(2) = nsendlist
      ldu(3) = MAXVAL(ppm_max_meshid)
      ldu(4) = ppm_max_topoid
      CALL ppm_alloc(ppm_mesh_ghost_blkstart,ldu,iopt,info)
      IF (info .NE. 0) THEN
          info = ppm_error_fatal
          CALL ppm_error(ppm_err_alloc,'ppm_map_field_ghost_init',     &
     &        'ghost block start list PPM_MESH_GHOST_BLKSTART',__LINE__,info)
          GOTO 9999
      ENDIF
      CALL ppm_alloc(ppm_mesh_ghost_recvblkstart,ldu,iopt,info)
      IF (info .NE. 0) THEN
          info = ppm_error_fatal
          CALL ppm_error(ppm_err_alloc,'ppm_map_field_ghost_init',     &
     &        'ghost recv block start list PPM_MESH_GHOST_RECVBLKSTART',  &
     &        __LINE__,info)
          GOTO 9999
      ENDIF
      CALL ppm_alloc(ppm_mesh_ghost_blksize,ldu,iopt,info)
      IF (info .NE. 0) THEN
          info = ppm_error_fatal
          CALL ppm_error(ppm_err_alloc,'ppm_map_field_ghost_init',     &
     &        'ghost block size list PPM_MESH_GHOST_BLKSIZE',__LINE__,info)
          GOTO 9999
      ENDIF
      CALL ppm_alloc(ppm_mesh_ghost_recvblksize,ldu,iopt,info)
      IF (info .NE. 0) THEN
          info = ppm_error_fatal
          CALL ppm_error(ppm_err_alloc,'ppm_map_field_ghost_init',     &
     &        'ghost recv block size list PPM_MESH_GHOST_RECVBLKSIZE', &
     &        __LINE__,info)
          GOTO 9999
      ENDIF

      !-------------------------------------------------------------------------
      !  Allocate local memory for sorted offset list
      !-------------------------------------------------------------------------
      ldu(1) = pdim
      ldu(2) = nsendlist
      CALL ppm_alloc(mesh_ghost_offset,ldu,iopt,info)
      IF (info .NE. 0) THEN
          info = ppm_error_fatal
          CALL ppm_error(ppm_err_alloc,'ppm_map_field_ghost_init',     &
     &        'ghost block offset list MESH_GHOST_OFFSET',__LINE__,info)
          GOTO 9999
      ENDIF

      !-------------------------------------------------------------------------
      !  loop over the neighboring processors according to the optimized
      !  communication sequence.
      !-------------------------------------------------------------------------
      ppm_mesh_ghost_nsend(meshid,topoid) = nsendlist
      ppm_mesh_ghost_nrecv(meshid,topoid) = 0
      ppm_mesh_ghost_blk(1,meshid,topoid) = 1
      ppm_mesh_ghost_recvblk(1,meshid,topoid) = 1
      iset               = 0
      ibuffer            = 0
      isize              = SIZE(ppm_mesh_ghost_fromsub,1)

      DO i=1,ppm_ncommseq(topoid)
         !----------------------------------------------------------------------
         !  get the next processor to send/recv with
         !----------------------------------------------------------------------
         sendrank = ppm_icommseq(i,topoid)
         recvrank = sendrank

         !----------------------------------------------------------------------
         !  Initialize mesh block pointers
         !----------------------------------------------------------------------
         ibuffer = i + 1
         ppm_mesh_ghost_blk(ibuffer,meshid,topoid) =    &
     &       ppm_mesh_ghost_blk(i,meshid,topoid)
         ppm_mesh_ghost_recvblk(ibuffer,meshid,topoid) =    &
     &       ppm_mesh_ghost_recvblk(i,meshid,topoid)

         !----------------------------------------------------------------------
         !  Only assign data if there is any communication in this round
         !----------------------------------------------------------------------
         IF (sendrank .GE. 0) THEN

             !------------------------------------------------------------------
             !  Find all mesh blocks that are sent and store them.
             !------------------------------------------------------------------
             DO j=1,nsendlist
                 IF (ppm_subs2proc(isendtosub(j),topoid) .EQ. sendrank) THEN
                     ppm_mesh_ghost_blk(ibuffer,meshid,topoid) =    &
     &                   ppm_mesh_ghost_blk(ibuffer,meshid,topoid)+1
                     iset = ppm_mesh_ghost_blk(ibuffer,meshid,topoid) - 1
                     ! store this for the topology as it can be reused
                     ppm_mesh_ghost_fromsub(iset,meshid,topoid)=isendfromsub(j)
                     ppm_mesh_ghost_tosub(iset,meshid,topoid)  =isendtosub(j)
                     ppm_mesh_ghost_blkstart(1:pdim,iset,meshid,topoid) = &
     &                   isendblkstart(1:pdim,j)
                     ppm_mesh_ghost_blksize(1:pdim,iset,meshid,topoid)  = &
     &                   isendblksize(1:pdim,j)
                     ! also re-order the offsets as we need them for
                     ! computing the receive lists further down !!
                     mesh_ghost_offset(1:pdim,iset) = ioffset(1:pdim,j)
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
                         CALL ppm_write(ppm_rank,'ppm_map_field_ghost_init',  &
     &                       mesg,info)
                     ENDIF
                 ENDIF
             ENDDO

             !------------------------------------------------------------------
             !  Build receive lists for the on-processor part of the data
             !------------------------------------------------------------------
             IF (sendrank .EQ. ppm_rank) THEN
                 ub = ppm_mesh_ghost_blk(2,meshid,topoid)
                 ppm_mesh_ghost_nrecv(meshid,topoid) = ub - 1
                 ppm_mesh_ghost_recvblk(2,meshid,topoid) = ub
                 DO j=1,ub-1
                     ppm_mesh_ghost_recvtosub(j,meshid,topoid) =   &
     &                   ppm_mesh_ghost_tosub(j,meshid,topoid)
  !                   ppm_mesh_ghost_recvfromsub(j,meshid,topoid) =  &
  !   &                   ppm_mesh_ghost_fromsub(j,meshid,topoid)
                     ppm_mesh_ghost_recvblkstart(1:pdim,j,meshid,topoid) = &
     &                   ppm_mesh_ghost_blkstart(1:pdim,j,meshid,topoid) + &
     &                   mesh_ghost_offset(1:pdim,j)
                     ppm_mesh_ghost_recvblksize(1:pdim,j,meshid,topoid) = &
     &                   ppm_mesh_ghost_blksize(1:pdim,j,meshid,topoid)
                 ENDDO

#ifdef __MPI
             ELSE
                 !--------------------------------------------------------------
                 !  Communicate the block indices
                 !--------------------------------------------------------------
                 lb = ppm_mesh_ghost_blk(i,meshid,topoid)
                 ub = ppm_mesh_ghost_blk(ibuffer,meshid,topoid)
                 ! How many blocks am I sending to that guy
                 nsend = ub - lb
                 tag1 = 100
                 CALL MPI_SendRecv(nsend,1,MPI_INTEGER,sendrank,tag1,nrecv,1,  &
     &               MPI_INTEGER,recvrank,tag1,ppm_comm,commstat,info)
                 ! How many blocks will I receive from the guy?
                 ppm_mesh_ghost_nrecv(meshid,topoid) =      &
     &               ppm_mesh_ghost_nrecv(meshid,topoid) + nrecv
                 ppm_mesh_ghost_recvblk(ibuffer,meshid,topoid) =    &
     &               ppm_mesh_ghost_recvblk(i,meshid,topoid) + nrecv

                 !--------------------------------------------------------------
                 !  Check if receive lists are large enough
                 !--------------------------------------------------------------
                 IF (ppm_mesh_ghost_nrecv(meshid,topoid) .GT. isize) THEN
                     !----------------------------------------------------------
                     !  Grow receive lists if needed
                     !----------------------------------------------------------
                     isize = ppm_mesh_ghost_nrecv(meshid,topoid)
                     iopt = ppm_param_alloc_grow_preserve
                     ldu(1) = isize
                     ldu(2) = MAXVAL(ppm_max_meshid)
                     ldu(3) = ppm_max_topoid
   !                  CALL ppm_alloc(ppm_mesh_ghost_recvfromsub,ldu,iopt,info)
   !                  IF (info .NE. 0) THEN
   !                      info = ppm_error_fatal
   !                      CALL ppm_error(ppm_err_alloc,   &
   !  &                     'ppm_map_field_ghost_init',     &
   !  &                     'ghost recv source subs PPM_MESH_GHOST_RECVFROMSUB',&
   !  &                     __LINE__,info)
   !                      GOTO 9999
   !                  ENDIF
                     CALL ppm_alloc(ppm_mesh_ghost_recvtosub,ldu,iopt,info)
                     IF (info .NE. 0) THEN
                         info = ppm_error_fatal
                         CALL ppm_error(ppm_err_alloc,   &
     &                     'ppm_map_field_ghost_init',     &
     &                     'ghost recv target subs PPM_MESH_GHOST_RECVTOSUB',&
     &                     __LINE__,info)
                         GOTO 9999
                     ENDIF
                     ldu(1) = pdim
                     ldu(2) = isize
                     ldu(3) = MAXVAL(ppm_max_meshid)
                     ldu(4) = ppm_max_topoid
                     CALL ppm_alloc(ppm_mesh_ghost_recvblkstart,ldu,iopt,info)
                     IF (info .NE. 0) THEN
                         info = ppm_error_fatal
                         CALL ppm_error(ppm_err_alloc,   &
     &                    'ppm_map_field_ghost_init',     &
     &                    'ghost recv block start PPM_MESH_GHOST_RECVBLKSTART',&
     &                     __LINE__,info)
                         GOTO 9999
                     ENDIF
                     CALL ppm_alloc(ppm_mesh_ghost_recvblksize,ldu,iopt,info)
                     IF (info .NE. 0) THEN
                         info = ppm_error_fatal
                         CALL ppm_error(ppm_err_alloc,   &
     &                    'ppm_map_field_ghost_init',     &
     &                    'ghost recv block size PPM_MESH_GHOST_RECVBLKSIZE',&
     &                     __LINE__,info)
                         GOTO 9999
                     ENDIF
                 ENDIF

                 !--------------------------------------------------------------
                 !  Allocate memory for block data send and recv buffers
                 !--------------------------------------------------------------
                 iopt   = ppm_param_alloc_grow
                 ldu(1) = nsend*(2*pdim+1)
                 CALL ppm_alloc(sendbuf,ldu,iopt,info)
                 IF (info .NE. 0) THEN
                     info = ppm_error_fatal
                     CALL ppm_error(ppm_err_alloc,'ppm_map_field_ghost_init', &
     &                   'local MPI send buffer SENDBUF',__LINE__,info)
                     GOTO 9999
                 ENDIF
                 ldu(1) = nrecv*(2*pdim+1)
                 CALL ppm_alloc(recvbuf,ldu,iopt,info)
                 IF (info .NE. 0) THEN
                     info = ppm_error_fatal
                     CALL ppm_error(ppm_err_alloc,'ppm_map_field_ghost_init', &
     &                   'local MPI recv buffer RECVBUF',__LINE__,info)
                     GOTO 9999
                 ENDIF
    
                 !--------------------------------------------------------------
                 !  Pack and send all the ghost mesh block data
                 !--------------------------------------------------------------
                 ! Pack all the send data
                 iset = 0
                 DO j=lb,ub-1
                     iset = iset + 1
                     sendbuf(iset) = ppm_mesh_ghost_tosub(j,meshid,topoid)
             !        iset = iset + 1
             !        sendbuf(iset) = ppm_mesh_ghost_fromsub(j,meshid,topoid)
                     sendbuf(iset+1:iset+pdim) =    &
     &                   ppm_mesh_ghost_blkstart(1:pdim,j,meshid,topoid) +  &
     &                   mesh_ghost_offset(1:pdim,j)
                     iset = iset + pdim
                     sendbuf(iset+1:iset+pdim) =    &
     &                   ppm_mesh_ghost_blksize(1:pdim,j,meshid,topoid)
                     iset = iset + pdim
                 ENDDO
                 ! Send it to the destination processor and get my stuff
                 tag1 = 200
                 CALL MPI_SendRecv(sendbuf,iset,MPI_INTEGER,sendrank,tag1, &
     &                             recvbuf,nrecv*(2*pdim+1),MPI_INTEGER,   &
     &                             recvrank,tag1,ppm_comm,commstat,info)
                 ! Unpack the received data
                 lb = ppm_mesh_ghost_recvblk(i,meshid,topoid)
                 ub = ppm_mesh_ghost_recvblk(ibuffer,meshid,topoid)
                 iset = 0
                 DO j=lb,ub-1
                     iset = iset + 1
                     ppm_mesh_ghost_recvtosub(j,meshid,topoid) = recvbuf(iset)
                 !    iset = iset + 1
                 !    ppm_mesh_ghost_recvfromsub(j,meshid,topoid) = recvbuf(iset)
                     ppm_mesh_ghost_recvblkstart(1:pdim,j,meshid,topoid) =   &
     &                   recvbuf(iset+1:iset+pdim)
                     iset = iset + pdim
                     ppm_mesh_ghost_recvblksize(1:pdim,j,meshid,topoid) =    &
     &                   recvbuf(iset+1:iset+pdim)
                     iset = iset + pdim

                 ENDDO
#endif
             ENDIF ! sendrank.EQ.ppm_rank
          ENDIF    ! sendrank .GE. 0
      ENDDO

      !-------------------------------------------------------------------------
      !  Deallocate memory of the local lists
      !-------------------------------------------------------------------------
      iopt   = ppm_param_dealloc
      CALL ppm_alloc(recvbuf,ldu,iopt,info)
      IF (info .NE. 0) THEN
          info = ppm_error_error
          CALL ppm_error(ppm_err_dealloc,'ppm_map_field_ghost_init',     &
     &        'local MPI recv buffer RECVBUF',__LINE__,info)
      ENDIF
      CALL ppm_alloc(sendbuf,ldu,iopt,info)
      IF (info .NE. 0) THEN
          info = ppm_error_error
          CALL ppm_error(ppm_err_dealloc,'ppm_map_field_ghost_init',     &
     &        'local MPI send buffer SENDBUF',__LINE__,info)
      ENDIF
      CALL ppm_alloc(isendfromsub,ldu,iopt,info)
      IF (info .NE. 0) THEN
          info = ppm_error_error
          CALL ppm_error(ppm_err_dealloc,'ppm_map_field_ghost_init',     &
     &        'local send source sub list ISENDFROMSUB',__LINE__,info)
      ENDIF
      CALL ppm_alloc(isendtosub,ldu,iopt,info)
      IF (info .NE. 0) THEN
          info = ppm_error_error
          CALL ppm_error(ppm_err_dealloc,'ppm_map_field_ghost_init',     &
     &        'local send destination sub list ISENDTOSUB',__LINE__,info)
      ENDIF
      CALL ppm_alloc(isendblkstart,ldu,iopt,info)
      IF (info .NE. 0) THEN
          info = ppm_error_error
          CALL ppm_error(ppm_err_dealloc,'ppm_map_field_ghost_init',     &
     &        'local send block start list ISENDBLKSTART',__LINE__,info)
      ENDIF
      CALL ppm_alloc(isendblksize,ldu,iopt,info)
      IF (info .NE. 0) THEN
          info = ppm_error_error
          CALL ppm_error(ppm_err_dealloc,'ppm_map_field_ghost_init',     &
     &        'local send block size list ISENDBLKSIZE',__LINE__,info)
      ENDIF
      CALL ppm_alloc(ioffset,ldu,iopt,info)
      IF (info .NE. 0) THEN
          info = ppm_error_error
          CALL ppm_error(ppm_err_dealloc,'ppm_map_field_ghost_init',     &
     &        'local send block offset list IOFFSET',__LINE__,info)
      ENDIF
      CALL ppm_alloc(mesh_ghost_offset,ldu,iopt,info)
      IF (info .NE. 0) THEN
          info = ppm_error_error
          CALL ppm_error(ppm_err_dealloc,'ppm_map_field_ghost_init',     &
     &        'ghost block offset list MESH_GHOST_OFFSET',__LINE__,info)
      ENDIF

      !-------------------------------------------------------------------------
      !  Return 
      !-------------------------------------------------------------------------
 9999 CONTINUE
      CALL substop('ppm_map_field_ghost_init',t0,info)
      RETURN
      END SUBROUTINE ppm_map_field_ghost_init
