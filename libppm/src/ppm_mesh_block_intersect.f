      !-------------------------------------------------------------------------
      !  Subroutine   :              ppm_mesh_block_intersect
      !-------------------------------------------------------------------------
      !
      !  Purpose      : This routine determines common mesh blocks
      !                 (intersections) of two subs, possibly shifted
      !                 and/or extendet with ghostlayers. Start and size of
      !                 all blocks found are stored and returned. The
      !                 routine is used by ppm_map_field_global and
      !                 ppm_map_field_ghost_init.
      !
      !  Input        : isub         (I) source sub (from which data will
      !                                  be sent) in global numbering
      !                 jsub         (I) destination sub in global numbering
      !                 frommesh     (I) source mesh identifier 
      !                                  (internal numbering)
      !                 tomesh       (I) target mesh identifier (internal 
      !                                  numbering)
      !                 fromtopo     (I) source topology identifier (internal 
      !                                  numbering) 
      !                 totopo       (I) target topology identifier (internal 
      !                                  numbering) 
      !                 offset(:)    (I) Shift offset for the source
      !                                  subdomain. Index: 1:ppm_dim
      !                 ghostsize(:) (I) size of the ghost layer in numbers
      !                                  of grid points in all space
      !                                  dimensions (1...ppm_dim). The
      !                                  target subdomain will be enlarged
      !                                  by this.
      !
      !  Input/output : nsendlist    (I) number of mesh blocks in the lists
      !                                  so far
      !
      !  Output       : isendfromsub(:)   (I) global sub index of sources
      !                 isendtosub(:)     (I) global sub index of targets
      !                 isenblkstart(:,:) (I) start of the mesh blocks in
      !                                       the global mesh. 1st index:
      !                                       1:ppm_dim, 2nd: meshblock.
      !                 isendblksize(:,:) (I) size of the meshblocks in
      !                                       numbers of grid points. 1st
      !                                       index: 1:ppm_dim, 2nd:
      !                                       meshblock
      !                 ioffset(:,:)      (I) meshblock offset for periodic
      !                                       images. 1st index:
      !                                       1:ppm_dim, 2nd: meshblock. To
      !                                       be added to isendblkstart in
      !                                       order to get irecvblkstart on
      !                                       the target subdomain.
      !                 info              (I) return status. 0 on success.
      !
      !  Remarks      : This does not check whether the frommesh and the
      !                 tomesh are compatible (i.e. have the same Nm). This
      !                 check must be done by the calling routine.
      !
      !                 The 5 output lists also need to be allocated by the
      !                 calling routine. They are just grown here if
      !                 needed. The reason is that this routine is
      !                 typically called from inside a loop and we do not
      !                 want to do the checks every time.
      !
      !  References   :
      !
      !  Revisions    :
      !-------------------------------------------------------------------------
      !  $Log: ppm_mesh_block_intersect.f,v $
      !  Revision 1.1.1.1  2007/07/13 10:18:56  ivos
      !  CBL version of the PPM library
      !
      !  Revision 1.6  2004/10/01 16:33:37  ivos
      !  cosmetics.
      !
      !  Revision 1.5  2004/10/01 16:09:09  ivos
      !  Replaced REAL(ppm_kind_double) :: t0 with REAL(MK) t0.
      !
      !  Revision 1.4  2004/07/26 11:48:09  ivos
      !  Fixes to make it compile.
      !
      !  Revision 1.3  2004/07/26 07:42:47  ivos
      !  Changed to single-interface modules and adapted all USE statements.
      !
      !  Revision 1.2  2004/04/07 16:54:44  ivos
      !  bugfix: offset of course only needs to be added to istart and not
      !  ndata :-) It caused problems when field indices went negative.
      !
      !  Revision 1.1  2004/04/07 15:33:11  ivos
      !  Initial implementation. Tested.
      !
      !-------------------------------------------------------------------------
      !  Parallel Particle Mesh Library (PPM)
      !  Institute of Computational Science
      !  ETH Zentrum, Hirschengraben 84
      !  CH-8092 Zurich, Switzerland
      !-------------------------------------------------------------------------

      SUBROUTINE ppm_mesh_block_intersect(isub,jsub,frommesh,tomesh,       &
     &    fromtopo,totopo,offset,ghostsize,nsendlist,isendfromsub,         &
     &    isendtosub,isendblkstart,isendblksize,ioffset,info)

      !-------------------------------------------------------------------------
      !  Modules 
      !-------------------------------------------------------------------------
      USE ppm_module_data
      USE ppm_module_data_mesh
      USE ppm_module_alloc
      USE ppm_module_error
      USE ppm_module_substart
      USE ppm_module_substop
      IMPLICIT NONE
      !-------------------------------------------------------------------------
      !  Includes
      !-------------------------------------------------------------------------
#include "ppm_define.h"
      !-------------------------------------------------------------------------
      !  Arguments     
      !-------------------------------------------------------------------------
      INTEGER                 , INTENT(IN   ) :: fromtopo,totopo,frommesh,   &
     &                                           tomesh,isub,jsub
      INTEGER, DIMENSION(:)   , INTENT(IN   ) :: ghostsize,offset
      INTEGER, DIMENSION(:)   , POINTER       :: isendfromsub,isendtosub
      INTEGER, DIMENSION(:,:) , POINTER       :: ioffset
      INTEGER, DIMENSION(:,:) , POINTER       :: isendblkstart,isendblksize
      INTEGER                 , INTENT(INOUT) :: nsendlist
      INTEGER                 , INTENT(  OUT) :: info
      !-------------------------------------------------------------------------
      !  Local variables 
      !-------------------------------------------------------------------------
      INTEGER, DIMENSION(3)            :: ldu
      INTEGER, DIMENSION(ppm_dim)      :: iblockstart,nblocksize
      INTEGER                          :: iblockstopk,k,iopt,pdim,isize
      INTEGER                          :: fromistartk,toistartk
      INTEGER                          :: fromnnodesk,tonnodesk
      LOGICAL                          :: dosend
      REAL(ppm_kind_double)            :: t0
      !-------------------------------------------------------------------------
      !  Externals 
      !-------------------------------------------------------------------------
      
      !-------------------------------------------------------------------------
      !  Initialise 
      !-------------------------------------------------------------------------
      CALL substart('ppm_mesh_block_intersect',t0,info)
      pdim = ppm_dim

      !-------------------------------------------------------------------------
      !  Check arguments
      !-------------------------------------------------------------------------
      IF (ppm_debug .GT. 0) THEN
          IF ((fromtopo .LE. 0) .OR. (fromtopo .GT. ppm_max_topoid)) THEN
              info = ppm_error_error
              CALL ppm_error(ppm_err_argument,'ppm_mesh_block_intersect',  &
     &            'fromtopo out of range',__LINE__,info)
              GOTO 9999
          ENDIF
          IF ((totopo .LE. 0) .OR. (totopo .GT. ppm_max_topoid)) THEN
              info = ppm_error_error
              CALL ppm_error(ppm_err_argument,'ppm_mesh_block_intersect',  &
     &            'totopo out of range',__LINE__,info)
              GOTO 9999
          ENDIF
          IF ((frommesh.LE.0).OR.(frommesh.GT.ppm_max_meshid(fromtopo))) THEN
              info = ppm_error_error
              CALL ppm_error(ppm_err_argument,'ppm_mesh_block_intersect',  &
     &            'frommesh out of range',__LINE__,info)
              GOTO 9999
          ENDIF
          IF ((tomesh.LE.0).OR.(tomesh.GT.ppm_max_meshid(totopo))) THEN
              info = ppm_error_error
              CALL ppm_error(ppm_err_argument,'ppm_mesh_block_intersect',  &
     &            'tomesh out of range',__LINE__,info)
              GOTO 9999
          ENDIF
          IF (SIZE(ghostsize,1) .LT. ppm_dim) THEN
              info = ppm_error_error
              CALL ppm_error(ppm_err_argument,'ppm_mesh_block_intersect',  &
     &            'ghostsize must be at least of length ppm_dim',__LINE__,info)
              GOTO 9999
          ENDIF
          IF (SIZE(offset,1) .LT. ppm_dim) THEN
              info = ppm_error_error
              CALL ppm_error(ppm_err_argument,'ppm_mesh_block_intersect',  &
     &            'offset must be at least of length ppm_dim',__LINE__,info)
              GOTO 9999
          ENDIF
          IF ((isub.LE.0).OR.(isub.GT.ppm_nsubs(fromtopo))) THEN
              info = ppm_error_error
              CALL ppm_error(ppm_err_argument,'ppm_mesh_block_intersect',  &
     &            'isub out of range',__LINE__,info)
              GOTO 9999
          ENDIF
          IF ((jsub.LE.0).OR.(jsub.GT.ppm_nsubs(totopo))) THEN
              info = ppm_error_error
              CALL ppm_error(ppm_err_argument,'ppm_mesh_block_intersect',  &
     &            'jsub out of range',__LINE__,info)
              GOTO 9999
          ENDIF
      ENDIF

      !-------------------------------------------------------------------------
      !  Determine current minimal common length of lists
      !-------------------------------------------------------------------------
      isize = MIN(SIZE(isendfromsub,1),SIZE(isendtosub,1),        &
     &    SIZE(isendblkstart,2),SIZE(isendblksize,2),SIZE(ioffset,2))

      !-------------------------------------------------------------------------
      !  Determine intersecting mesh blocks
      !-------------------------------------------------------------------------
      dosend = .TRUE.
      DO k=1,pdim
          fromistartk=ppm_cart_mesh(frommesh,fromtopo)%istart(k,isub)+offset(k)
          fromnnodesk=ppm_cart_mesh(frommesh,fromtopo)%nnodes(k,isub)
          toistartk  =ppm_cart_mesh(tomesh,totopo)%istart(k,jsub)-ghostsize(k)
          tonnodesk  =ppm_cart_mesh(tomesh,totopo)%nnodes(k,jsub)+2*ghostsize(k)
          iblockstart(k) = MAX(fromistartk,toistartk)
          iblockstopk = MIN((fromistartk+fromnnodesk),(toistartk+tonnodesk))
          nblocksize(k)  = iblockstopk-iblockstart(k)
          ! do not send if there is nothing to be sent
          IF (nblocksize(k) .LT. 1) dosend = .FALSE.
      ENDDO
      IF (dosend) THEN
          nsendlist = nsendlist + 1
          IF (nsendlist .GT. isize) THEN
              !-----------------------------------------------------------------
              !  Grow memory for the sendlists
              !-----------------------------------------------------------------
              isize  = 2*isize
              iopt   = ppm_param_alloc_grow_preserve
              ldu(1) = isize
              CALL ppm_alloc(isendfromsub,ldu,iopt,info)
              IF (info .NE. 0) THEN
                  info = ppm_error_fatal
                  CALL ppm_error(ppm_err_alloc,'ppm_mesh_block_intersect',   &
     &                'local send source sub list ISENDFROMSUB',__LINE__,info)
                  GOTO 9999
              ENDIF
              CALL ppm_alloc(isendtosub,ldu,iopt,info)
              IF (info .NE. 0) THEN
                  info = ppm_error_fatal
                  CALL ppm_error(ppm_err_alloc,'ppm_mesh_block_intersect',   &
     &               'local send destination sub list ISENDTOSUB',__LINE__,info)
                  GOTO 9999
              ENDIF
              ldu(1) = pdim
              ldu(2) = isize
              CALL ppm_alloc(isendblkstart,ldu,iopt,info)
              IF (info .NE. 0) THEN
                  info = ppm_error_fatal
                  CALL ppm_error(ppm_err_alloc,'ppm_mesh_block_intersect',   &
     &               'local send block start list ISENDBLKSTART',__LINE__,info)
                  GOTO 9999
              ENDIF
              CALL ppm_alloc(isendblksize,ldu,iopt,info)
              IF (info .NE. 0) THEN
                  info = ppm_error_fatal
                  CALL ppm_error(ppm_err_alloc,'ppm_mesh_block_intersect',   &
     &               'local send block size list ISENDBLKSIZE',__LINE__,info)
                  GOTO 9999
              ENDIF
              CALL ppm_alloc(ioffset,ldu,iopt,info)
              IF (info .NE. 0) THEN
                  info = ppm_error_fatal
                  CALL ppm_error(ppm_err_alloc,'ppm_mesh_block_intersect',   &
     &               'local send block offset list IOFFSET',__LINE__,info)
                  GOTO 9999
              ENDIF
          ENDIF
          !---------------------------------------------------------------------
          !  send block (isendblkstart...isendblkstart+nblocksize-1) 
          !  from sub isub to sub jsub
          !---------------------------------------------------------------------
          ! global sub index of local sub where the blocks come from
          isendfromsub(nsendlist)        = isub
          ! global sub index of other sub where the blocks go to
          isendtosub(nsendlist)          = jsub
          ! start of mesh block in global mesh
          isendblkstart(1:pdim,nsendlist)= iblockstart(1:pdim) - offset(1:pdim)
          ! size of mesh block
          isendblksize(1:pdim,nsendlist) = nblocksize(1:pdim)
          ! mesh block offset for periodic images
          ioffset(1:pdim,nsendlist)      = offset(1:pdim)
      ENDIF

      !-------------------------------------------------------------------------
      !  Return 
      !-------------------------------------------------------------------------
 9999 CONTINUE
      CALL substop('ppm_mesh_block_intersect',t0,info)
      RETURN
      END SUBROUTINE ppm_mesh_block_intersect
