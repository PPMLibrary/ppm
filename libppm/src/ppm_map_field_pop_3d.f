      !-------------------------------------------------------------------------
      !  Subroutine   :               ppm_map_field_pop_3d
      !-------------------------------------------------------------------------
      !
      !  Purpose      : This routine pops the list buffer for 3d mesh data
      !
      !  Input        : fdata([:,]:,:,:,:)(F) field data. 1st index: lda,
      !                                     2nd-4th: mesh (i,j,k) relative
      !                                     to istart-1 of the sub (i.e.
      !                                     (i,j,k)=(1,1,1) corresponds to
      !                                     istart), 5th: isub 1...nsublist
      !                                     (all subs on this processor).
      !                                     For scalar field data, the
      !                                     first index is omitted (the
      !                                     others shift accordingly).
      !                 lda             (I) the leading dimension of fdata.
      !                                     lda=1 for the case of scalar
      !                                     data.
      !                 tomesh          (I) mesh ID of destination
      !                                     (internal numbering)
      !                 ghostsize(:)    (I) number of mesh points in the
      !                                     ghost layer in each dimension
      !                 rtype           (I) receive type. One of
      !                                       ppm_param_pop_replace
      !                                       ppm_param_pop_add
      !                 mask(:,:,:,:)   (L) Logical mask. Only the mesh nodes
      !                                     for which this is .TRUE. will be
      !                                     mapped. OPTIONAL. If not given, all
      !                                     points are mapped. 1st-3rd
      !                                     index: mesh (i,j,k), 4th: isub.
      !
      !  Input/output :
      !
      !  Output       : info        (I) return status. 0 upon success.
      !
      !  Remarks      : The on-processor data is stored in the first part
      !                 of the buffer.
      !
      !                 Reallocates the fdata array if needed. The user
      !                 does not need to care about its size. If we want to
      !                 change this to improve callability from C, argument
      !                 checks should be added for the SIZE of fdata.
      !
      !                 If the call to invert_list is too slow here, this
      !                 could be done by topo_store and stored for each
      !                 topology. At the expense of some extra memory.
      !
      !  References   :
      !
      !  Revisions    :
      !-------------------------------------------------------------------------
      !  $Log: ppm_map_field_pop_3d.f,v $
      !  Revision 1.1.1.1  2007/07/13 10:18:56  ivos
      !  CBL version of the PPM library
      !
      !  Revision 1.30  2005/03/10 01:40:31  ivos
      !  Empty buffer is now only reported for ppm_debug.GT.0 and is no
      !  longer logged to prevent huge log files.
      !
      !  Revision 1.29  2005/02/16 22:07:05  ivos
      !  Unrolled loops for lda=1,2,3,4,5 to allow vectorization.
      !  If mask is present, no unrolling is done since this case will
      !  not vectorize anyway.
      !
      !  Revision 1.28  2004/11/11 15:23:02  ivos
      !  Moved allocatable work data to the module.
      !
      !  Revision 1.27  2004/10/01 16:09:05  ivos
      !  Replaced REAL(ppm_kind_double) :: t0 with REAL(MK) t0.
      !
      !  Revision 1.26  2004/09/22 10:43:24  ivos
      !  send/recv-buffers are now deallocated after their use to save
      !  memory. If this is a performance issue and memory is not a
      !  problem, we can remove this again.
      !
      !  Revision 1.25  2004/09/15 07:47:37  ivos
      !  bugfix: points where mask.EQ..FALSE. were not skipped in the
      !  recvbuffer (since they were not read out). Fixed.
      !
      !  Revision 1.24  2004/08/31 12:48:09  ivos
      !  Changed argument checks to use ppm_check_topoid and ppm_check_meshid.
      !
      !  Revision 1.23  2004/08/30 09:26:07  ivos
      !  bugfix: replaced MB_MAP_FIELD_POP_HACK with proper fixes for all cases.
      !
      !  Revision 1.22  2004/08/19 14:38:38  michaebe
      !  bugfixed a bugfix
      !
      !  Revision 1.21  2004/08/19 14:31:42  michaebe
      !  included the __MB_MAP_FIELD_POP_HACK also for complex vectors
      !
      !  Revision 1.20  2004/08/19 14:23:41  michaebe
      !  inserted the temporary __MB_MAP_FIELD_POP_HACK.
      !
      !  Revision 1.19  2004/07/26 15:38:48  ivos
      !  Inserted missing USE statements to resolve undefined references
      !  at link stage.
      !
      !  Revision 1.18  2004/07/26 07:42:42  ivos
      !  Changed to single-interface modules and adapted all USE statements.
      !
      !  Revision 1.17  2004/07/22 10:53:54  hiebers
      !  Merged
      !
      !  Revision 1.16  2004/07/19 11:01:51  ivos
      !  Overloaded field push and pop operations for scalar fields and
      !  added all changes to the module and the interface routines.
      !
      !  Revision 1.15  2004/07/16 14:12:45  ivos
      !  ppm_nsendbuffer is now also decremented properly to allow reuse of
      !  the same mapping in multiple push-send-pop cycles.
      !
      !  Revision 1.14  2004/07/01 15:49:40  ivos
      !  bugfix: removed wrong ifdef on MPI. Serial version now runs as well.
      !
      !  Revision 1.13  2004/04/15 10:18:34  ivos
      !  bugfix: changed mask from INTENT(IN) to POINTER as it can have
      !  LBOUND .LT. 1 (if ghost layer is masked).
      !
      !  Revision 1.12  2004/04/14 15:08:51  ivos
      !  Added OPTIONAL argument mask to allow selective mapping of only a
      !  subset of mesh points. Different masks can be used for different
      !  fields on the same mesh. Not tested yet.
      !
      !  Revision 1.11  2004/04/07 09:20:49  ivos
      !  Realloc of fdata now takes ghostsize properly into account and iopt
      !  is set to _preserve if we are adding contributions back or only
      !  mapping ghosts, this leaving the actual field untouched.
      !
      !  Revision 1.10  2004/04/05 12:02:32  ivos
      !  Added parameter rtype and cases for ppm_param_pop_add.
      !
      !  Revision 1.9  2004/04/05 10:45:29  ivos
      !  Translation to local sub indices now uses invert_list outside of the
      !  loop instead of searching isublist every time.
      !
      !  Revision 1.8  2004/04/02 15:32:02  ivos
      !  Changed to use the new 5d ppm_alloc instead of direct ALLOCATE.
      !
      !  Revision 1.7  2004/04/02 15:22:23  ivos
      !  Added translation from global to local sub ID, since all sub IDs in
      !  the ppm_mesh_*sub lists are now global (because of ghosts).
      !
      !  Revision 1.6  2004/03/05 13:43:23  ivos
      !  bugfix: index errors corrected for COMPLEX version. COMPLEX version
      !  is now tested.
      !
      !  Revision 1.5  2004/02/24 12:16:04  ivos
      !  map_field_pop now grows the field data array if needed instead of just
      !  checking its size and requiring the user to allocate it. This made all
      !  arguments change from INTENT(INOUT) to POINTER. Reason for the change:
      !  symmetry with the particle pop.
      !
      !  Revision 1.4  2004/02/24 08:46:54  ivos
      !  Added overloading for single complex and double complex data types.
      !
      !  Revision 1.3  2004/02/23 12:19:00  ivos
      !  Several bugs fixed. Tested on 2 processors with a scalar field.
      !  Added debug output in several places.
      !
      !  Revision 1.2  2004/02/23 08:56:39  ivos
      !  bugfix: necessary size of data array (for argument check) was determined
      !  using global data instead of on-processor data. fixed.
      !
      !  Revision 1.1  2004/02/17 16:12:31  ivos
      !  Initial implementation. Not tested yet.
      !
      !-------------------------------------------------------------------------
      !  Parallel Particle Mesh Library (PPM)
      !  Institute of Computational Science
      !  ETH Zentrum, Hirschengraben 84
      !  CH-8092 Zurich, Switzerland
      !-------------------------------------------------------------------------

      ! ATTN: DIM is the dimension of the fdata array and not the space
      ! dimemsion ppm_dim!
#if    __DIM == __SFIELD
#if    __KIND == __SINGLE_PRECISION
      SUBROUTINE ppm_map_field_pop_3d_sca_s(fdata,lda,tomesh,ghostsize,rtype, &
     &    info,mask)
#elif  __KIND == __DOUBLE_PRECISION
      SUBROUTINE ppm_map_field_pop_3d_sca_d(fdata,lda,tomesh,ghostsize,rtype, &
     &    info,mask)
#elif  __KIND == __SINGLE_PRECISION_COMPLEX
      SUBROUTINE ppm_map_field_pop_3d_sca_sc(fdata,lda,tomesh,ghostsize,rtype,& 
     &    info,mask)
#elif  __KIND == __DOUBLE_PRECISION_COMPLEX
      SUBROUTINE ppm_map_field_pop_3d_sca_dc(fdata,lda,tomesh,ghostsize,rtype,&
     &    info,mask)
#elif  __KIND == __INTEGER
      SUBROUTINE ppm_map_field_pop_3d_sca_i(fdata,lda,tomesh,ghostsize,rtype, &
     &    info,mask)
#elif  __KIND == __LOGICAL
      SUBROUTINE ppm_map_field_pop_3d_sca_l(fdata,lda,tomesh,ghostsize,rtype, &
     &    info,mask)
#endif 

#elif  __DIM == __VFIELD
#if    __KIND == __SINGLE_PRECISION
      SUBROUTINE ppm_map_field_pop_3d_vec_s(fdata,lda,tomesh,ghostsize,rtype, &
     &    info,mask)
#elif  __KIND == __DOUBLE_PRECISION
      SUBROUTINE ppm_map_field_pop_3d_vec_d(fdata,lda,tomesh,ghostsize,rtype, &
     &    info,mask)
#elif  __KIND == __SINGLE_PRECISION_COMPLEX
      SUBROUTINE ppm_map_field_pop_3d_vec_sc(fdata,lda,tomesh,ghostsize,rtype,&
     &    info,mask)
#elif  __KIND == __DOUBLE_PRECISION_COMPLEX
      SUBROUTINE ppm_map_field_pop_3d_vec_dc(fdata,lda,tomesh,ghostsize,rtype,&
     &    info,mask)
#elif  __KIND == __INTEGER
      SUBROUTINE ppm_map_field_pop_3d_vec_i(fdata,lda,tomesh,ghostsize,rtype, &
     &    info,mask)
#elif  __KIND == __LOGICAL
      SUBROUTINE ppm_map_field_pop_3d_vec_l(fdata,lda,tomesh,ghostsize,rtype, &
     &    info,mask)
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
      USE ppm_module_data_mesh
      USE ppm_module_substart
      USE ppm_module_substop
      USE ppm_module_error
      USE ppm_module_alloc
      USE ppm_module_write
      USE ppm_module_check_meshid
      USE ppm_module_util_invert_list
      IMPLICIT NONE
#if    __KIND == __SINGLE_PRECISION | __KIND == __SINGLE_PRECISION_COMPLEX
      INTEGER, PARAMETER :: MK = ppm_kind_single
#else
      INTEGER, PARAMETER :: MK = ppm_kind_double
#endif
      !-------------------------------------------------------------------------
      !  Arguments     
      !-------------------------------------------------------------------------
#if   __DIM == __SFIELD
#if   __KIND == __INTEGER
      INTEGER , DIMENSION(:,:,:,:)   , POINTER         :: fdata
#elif __KIND == __LOGICAL
      LOGICAL , DIMENSION(:,:,:,:)   , POINTER         :: fdata
#elif __KIND == __SINGLE_PRECISION_COMPLEX | __KIND == __DOUBLE_PRECISION_COMPLEX
      COMPLEX(MK), DIMENSION(:,:,:,:), POINTER         :: fdata
#else
      REAL(MK), DIMENSION(:,:,:,:)   , POINTER         :: fdata
#endif

#elif __DIM == __VFIELD
#if   __KIND == __INTEGER
      INTEGER , DIMENSION(:,:,:,:,:)   , POINTER       :: fdata
#elif __KIND == __LOGICAL
      LOGICAL , DIMENSION(:,:,:,:,:)   , POINTER       :: fdata
#elif __KIND == __SINGLE_PRECISION_COMPLEX | __KIND == __DOUBLE_PRECISION_COMPLEX
      COMPLEX(MK), DIMENSION(:,:,:,:,:), POINTER       :: fdata
#else
      REAL(MK), DIMENSION(:,:,:,:,:)   , POINTER       :: fdata
#endif
#endif
      LOGICAL, DIMENSION(:,:,:,:),   POINTER, OPTIONAL :: mask
      INTEGER                          , INTENT(IN   ) :: lda,tomesh,rtype
      INTEGER, DIMENSION(:)            , INTENT(IN   ) :: ghostsize
      INTEGER                          , INTENT(  OUT) :: info
      !-------------------------------------------------------------------------
      !  Local variables 
      !-------------------------------------------------------------------------
      INTEGER, DIMENSION(3)   :: mofs
      INTEGER, DIMENSION(5)   :: ldu,ldl
      INTEGER                 :: i,j,k,ibuffer,Mdata,imesh,jmesh,kmesh,isub
      INTEGER                 :: iopt,totopo,xhi,yhi,zhi,xlo,ylo,zlo,bdim,jsub
      INTEGER                 :: btype,edim,idom
      REAL(MK)                :: t0
      LOGICAL                 :: ldo
      CHARACTER(LEN=ppm_char) :: mesg
      !-------------------------------------------------------------------------
      !  Externals 
      !-------------------------------------------------------------------------
      
      !-------------------------------------------------------------------------
      !  Initialise 
      !-------------------------------------------------------------------------
      CALL substart('ppm_map_field_pop_3d',t0,info)
      totopo = ppm_field_topoid

      !-------------------------------------------------------------------------
      !  Check arguments
      !-------------------------------------------------------------------------
      IF (ppm_debug .GT. 0) THEN
          CALL ppm_check_meshid(ppm_param_id_internal,tomesh,totopo,ldo,info)
          IF (ppm_buffer_set .LT. 1) THEN
              info = ppm_error_error
              CALL ppm_error(ppm_err_argument,'ppm_map_field_pop_3d',  &
     &            'buffer is empty. Cannot pop.',__LINE__,info)
              GOTO 9999
          ENDIF
          IF (.NOT. ldo) THEN
              info = ppm_error_error
              CALL ppm_error(ppm_err_argument,'ppm_map_field_pop_3d',  &
     &            'destination meshid out of range',__LINE__,info)
              GOTO 9999
          ENDIF
#if   __DIM == __VFIELD
          IF (lda .LT. 1) THEN
              info = ppm_error_error
              CALL ppm_error(ppm_err_argument,'ppm_map_field_pop_3d',  &
     &            'lda must be >0 for vector data',__LINE__,info)
              GOTO 9999
          ENDIF
#elif __DIM == __SFIELD
          IF (lda .NE. 1) THEN
              info = ppm_error_error
              CALL ppm_error(ppm_err_argument,'ppm_map_field_pop_3d',  &
     &            'lda must be =1 for scalar data',__LINE__,info)
              GOTO 9999
          ENDIF
#endif

          IF (SIZE(ghostsize,1) .LT. 3) THEN
              info = ppm_error_error
              CALL ppm_error(ppm_err_argument,'ppm_map_field_pop_3d',  &
     &            'ghostsize must be given for all dimensions',__LINE__,info)
              GOTO 9999
          ENDIF
          IF ((ghostsize(1) .LT. 0) .OR. (ghostsize(2) .LT. 0) .OR.  &
     &        (ghostsize(3) .LT. 0)) THEN
              info = ppm_error_error
              CALL ppm_error(ppm_err_argument,'ppm_map_field_pop_3d',  &
     &            'ghostsize must be >=0 in all dimensions',__LINE__,info)
              GOTO 9999
          ENDIF
          IF ((rtype .NE. ppm_param_pop_replace) .AND.     &
     &        (rtype .NE. ppm_param_pop_add)) THEN
              info = ppm_error_error
              CALL ppm_error(ppm_err_argument,'ppm_map_field_pop_3d',  &
     &            'Unknown receive type specified',__LINE__,info)
              GOTO 9999
          ENDIF
          IF (PRESENT(mask)) THEN
              xhi = 0
              yhi = 0
              zhi = 0
              DO i=1,ppm_nsublist(totopo)
                  isub = ppm_isublist(i,totopo)
                  IF (ppm_cart_mesh(tomesh,totopo)%nnodes(1,isub).GT.xhi) THEN
                      xhi = ppm_cart_mesh(tomesh,totopo)%nnodes(1,isub)
                  ENDIF
                  IF (ppm_cart_mesh(tomesh,totopo)%nnodes(2,isub).GT.yhi) THEN
                      yhi = ppm_cart_mesh(tomesh,totopo)%nnodes(2,isub)
                  ENDIF
                  IF (ppm_cart_mesh(tomesh,totopo)%nnodes(3,isub).GT.zhi) THEN
                      zhi = ppm_cart_mesh(tomesh,totopo)%nnodes(3,isub)
                  ENDIF
              ENDDO
              IF (SIZE(mask,1) .LT. xhi) THEN
                  info = ppm_error_error
                  CALL ppm_error(ppm_err_argument,'ppm_map_field_pop_3d',  &
     &                'x dimension of mask does not match mesh',__LINE__,info)
                  GOTO 9999
              ENDIF
              IF (SIZE(mask,2) .LT. yhi) THEN
                  info = ppm_error_error
                  CALL ppm_error(ppm_err_argument,'ppm_map_field_pop_3d',  &
     &                'y dimension of mask does not match mesh',__LINE__,info)
                  GOTO 9999
              ENDIF
              IF (SIZE(mask,3) .LT. zhi) THEN
                  info = ppm_error_error
                  CALL ppm_error(ppm_err_argument,'ppm_map_field_pop_3d',  &
     &                'z dimension of mask does not match mesh',__LINE__,info)
                  GOTO 9999
              ENDIF
              IF (SIZE(mask,4) .LT. ppm_nsublist(totopo)) THEN
                  info = ppm_error_error
                  CALL ppm_error(ppm_err_argument,'ppm_map_field_pop_3d',  &
     &                'mask data for some subs is missing',__LINE__,info)
                  GOTO 9999
              ENDIF
          ENDIF
      ENDIF

      !-------------------------------------------------------------------------
      !  Check that the required dimension fits the dimension of the buffer
      !-------------------------------------------------------------------------
      bdim = ppm_buffer_dim(ppm_buffer_set)
#if   __KIND == __SINGLE_PRECISION_COMPLEX | __KIND == __DOUBLE_PRECISION_COMPLEX
      ! for complex the effective dimension is half the data dimension
      edim = bdim/2
#else
      edim = bdim
#endif
      IF (ppm_debug .GT. 1) THEN
          WRITE(mesg,'(2(A,I3))') 'bdim=',edim,'    lda=',lda
          CALL ppm_write(ppm_rank,'ppm_map_field_pop_3d',mesg,info)
      ENDIF
#if   __DIM == __VFIELD
      IF (edim.NE.lda) THEN
         info = ppm_error_error
         CALL ppm_error(ppm_err_wrong_dim,'ppm_map_field_pop_3d',    &
     &       'leading dimension LDA is in error',__LINE__,info)
         GOTO 9999
      ENDIF 
#elif __DIM == __SFIELD
      IF (edim.NE.1) THEN
         info = ppm_error_error
         CALL ppm_error(ppm_err_wrong_dim,'ppm_map_field_pop_3d',    &
     &       'buffer does not contain 1d data!',__LINE__,info)
         GOTO 9999
      ENDIF 
#endif

      !-------------------------------------------------------------------------
      !  Check that the required type is identical to the type of the buffer
      !-------------------------------------------------------------------------
      btype = ppm_buffer_type(ppm_buffer_set)
#if    __KIND == __SINGLE_PRECISION
      IF (btype.NE.ppm_kind_single) THEN
         info = ppm_error_error
         CALL ppm_error(ppm_err_wrong_prec,'ppm_map_field_pop_3d',    &
     &       'trying to pop a non-single into single ',__LINE__,info)
         GOTO 9999
      ENDIF
#elif  __KIND == __DOUBLE_PRECISION
      IF (btype.NE.ppm_kind_double) THEN
         info = ppm_error_error
         CALL ppm_error(ppm_err_wrong_prec,'ppm_map_field_pop_3d',    &
     &       'trying to pop a non-double into double ',__LINE__,info)
         GOTO 9999
      ENDIF
#elif  __KIND == __SINGLE_PRECISION_COMPLEX
      IF (btype.NE.ppm_kind_single) THEN
         info = ppm_error_error
         CALL ppm_error(ppm_err_wrong_prec,'ppm_map_field_pop_3d',    &
     &       'trying to pop a non-single-complex into single-complex',&
     &       __LINE__,info)
         GOTO 9999
      ENDIF
#elif  __KIND == __DOUBLE_PRECISION_COMPLEX
      IF (btype.NE.ppm_kind_double) THEN
         info = ppm_error_error
         CALL ppm_error(ppm_err_wrong_prec,'ppm_map_field_pop_3d',    &
     &       'trying to pop a non-double-complex into double-complex',&
     &       __LINE__,info)
         GOTO 9999
      ENDIF
#elif  __KIND == __INTEGER
      IF (btype.NE.ppm_integer) THEN
         info = ppm_error_error
         CALL ppm_error(ppm_err_wrong_prec,'ppm_map_field_pop_3d',    &
     &       'trying to pop a non-integer into integer ',__LINE__,info)
         GOTO 9999
      ENDIF
#elif  __KIND == __LOGICAL
      IF (btype.NE.ppm_logical) THEN
         info = ppm_error_error
         CALL ppm_error(ppm_err_wrong_prec,'ppm_map_field_pop_3d',    &
     &       'trying to pop a non-logical into logical ',__LINE__,info)
      ENDIF
#endif

      !-------------------------------------------------------------------------
      !  Determine size of field data array needed
      !-------------------------------------------------------------------------
      xhi = 0
      yhi = 0
      zhi = 0
      DO i=1,ppm_nsublist(totopo)
          idom = ppm_isublist(i,totopo)
          IF (ppm_cart_mesh(tomesh,totopo)%nnodes(1,idom).GT.xhi) THEN
              xhi = ppm_cart_mesh(tomesh,totopo)%nnodes(1,idom)
          ENDIF
          IF (ppm_cart_mesh(tomesh,totopo)%nnodes(2,idom).GT.yhi) THEN
              yhi = ppm_cart_mesh(tomesh,totopo)%nnodes(2,idom)
          ENDIF
          IF (ppm_cart_mesh(tomesh,totopo)%nnodes(3,idom).GT.zhi) THEN
              zhi = ppm_cart_mesh(tomesh,totopo)%nnodes(3,idom)
          ENDIF
      ENDDO

      !-------------------------------------------------------------------------
      !  Reallocate array if needed
      !-------------------------------------------------------------------------
      IF ((ppm_map_type .EQ. ppm_param_map_ghost_get) .OR.   &
     &    (ppm_map_type .EQ. ppm_param_map_ghost_put) .OR.   &
     &    (rtype .EQ. ppm_param_pop_add)) THEN
          !---------------------------------------------------------------------
          !  Preserve old fields if this is to receive ghosts or to add
          !  contributions
          !---------------------------------------------------------------------
          iopt = ppm_param_alloc_fit_preserve
      ELSE
          iopt = ppm_param_alloc_fit
      ENDIF
#if   __DIM == __VFIELD
      ldu(1) = edim
      ldu(2) = xhi+ghostsize(1)
      ldu(3) = yhi+ghostsize(2)
      ldu(4) = zhi+ghostsize(3)
      ldu(5) = ppm_nsublist(totopo)
      ldl(1) = 1
      ldl(2) = 1-ghostsize(1)
      ldl(3) = 1-ghostsize(2)
      ldl(4) = 1-ghostsize(3)
      ldl(5) = 1
#elif __DIM == __SFIELD
      ldu(1) = xhi+ghostsize(1)
      ldu(2) = yhi+ghostsize(2)
      ldu(3) = zhi+ghostsize(3)
      ldu(4) = ppm_nsublist(totopo)
      ldl(1) = 1-ghostsize(1)
      ldl(2) = 1-ghostsize(2)
      ldl(3) = 1-ghostsize(3)
      ldl(4) = 1
#endif
      CALL ppm_alloc(fdata,ldl,ldu,iopt,info)
      IF (info .NE. 0) THEN
          info = ppm_error_fatal
          CALL ppm_error(ppm_err_alloc,'ppm_map_field_pop_3d',  &
     &        'new data FDATA',__LINE__,info)
          GOTO 9999
      ENDIF

      !-------------------------------------------------------------------------
      !  Build the inverse sub list to find local sub indeices based on
      !  global ones (the global ones are communicated)
      !-------------------------------------------------------------------------
      iopt   = ppm_param_alloc_fit
      ldu(1) = ppm_nsublist(ppm_field_topoid)
      CALL ppm_alloc(sublist,ldu,iopt,info)
      IF (info .NE. 0) THEN
          info = ppm_error_fatal
          CALL ppm_error(ppm_err_alloc,'ppm_map_field_pop_3d',  &
     &        'temporary sub list SUBLIST',__LINE__,info)
          GOTO 9999
      ENDIF
      ! We need to copy it into a temp list, since directly using
      ! ppm_isublist(:,ppm_field_topoid) as an argument to invert_list is
      ! not possible since the argument needs to be a POINTER.
      sublist(1:ppm_nsublist(ppm_field_topoid)) =     &
     &    ppm_isublist(1:ppm_nsublist(ppm_field_topoid),ppm_field_topoid)
      CALL ppm_util_invert_list(sublist,invsublist,info)
      iopt   = ppm_param_dealloc
      CALL ppm_alloc(sublist,ldu,iopt,info)
      IF (info .NE. 0) THEN
          info = ppm_error_error
          CALL ppm_error(ppm_err_dealloc,'ppm_map_field_pop_3d',     &
     &        'temporary sub list SUBLIST',__LINE__,info)
      ENDIF

      !-------------------------------------------------------------------------
      !  Determine the number of data points to be received
      !-------------------------------------------------------------------------
      Mdata = 0
 !     IF (PRESENT(mask)) THEN
 !        DO i=1,ppm_nrecvlist
 !           !-------------------------------------------------------------------
 !           !  access mesh blocks belonging to the i-th processor in the 
 !           !  recvlist
 !           !-------------------------------------------------------------------
 !           DO j=ppm_precvbuffer(i),ppm_precvbuffer(i+1)-1
 !              !----------------------------------------------------------------
 !              !  Determine mesh range in local sub coordinates
 !              !----------------------------------------------------------------
 !              jsub = ppm_mesh_irecvtosub(j)
 !              isub = invsublist(jsub)
 !              mofs(1) = ppm_cart_mesh(tomesh,totopo)%istart(1,jsub)-1
 !              mofs(2) = ppm_cart_mesh(tomesh,totopo)%istart(2,jsub)-1
 !              mofs(3) = ppm_cart_mesh(tomesh,totopo)%istart(3,jsub)-1
 !              xlo = ppm_mesh_irecvblkstart(1,j)-mofs(1)
 !              ylo = ppm_mesh_irecvblkstart(2,j)-mofs(2)
 !              zlo = ppm_mesh_irecvblkstart(3,j)-mofs(3)
 !              xhi = xlo+ppm_mesh_irecvblksize(1,j)-1
 !              yhi = ylo+ppm_mesh_irecvblksize(2,j)-1
 !              zhi = zlo+ppm_mesh_irecvblksize(3,j)-1
 !              !----------------------------------------------------------------
 !              !  Count points for which mask is true
 !              !----------------------------------------------------------------
 !              DO kmesh=zlo,zhi
 !                 DO jmesh=ylo,yhi
 !                    DO imesh=xlo,xhi
 !                       IF (mask(1,imesh,jmesh,kmesh,isub)) Mdata = Mdata + 1
 !                    ENDDO
 !                 ENDDO
 !              ENDDO
 !           ENDDO
 !        ENDDO
 !     ELSE
         DO i=1,ppm_nrecvlist
            !-------------------------------------------------------------------
            !  access mesh blocks belonging to the i-th processor in the 
            !  recvlist
            !-------------------------------------------------------------------
            DO j=ppm_precvbuffer(i),ppm_precvbuffer(i+1)-1
               !----------------------------------------------------------------
               !  Get the number of mesh points in this block
               !----------------------------------------------------------------
               Mdata = Mdata + (ppm_mesh_irecvblksize(1,j)*         &
     &            ppm_mesh_irecvblksize(2,j)*ppm_mesh_irecvblksize(3,j))
            ENDDO
         ENDDO
 !     ENDIF

      !-------------------------------------------------------------------------
      !  If there is nothing to be sent we are done
      !-------------------------------------------------------------------------
      IF (Mdata .EQ. 0) THEN
          IF (ppm_debug .GT. 0) THEN
              CALL ppm_write(ppm_rank,'ppm_map_field_pop_3d',   &
     &            'There is no data to be received',info)
          ENDIF
          GOTO 9999
      ENDIF

      !-------------------------------------------------------------------------
      !  Decrement the pointer into the receive buffer
      !-------------------------------------------------------------------------
      IF (ppm_debug .GT. 1) THEN
          WRITE(mesg,'(2(A,I9))') 'ppm_nrecvbuffer = ',ppm_nrecvbuffer,   &
     &        'Mdata*bdim = ',Mdata*bdim
          CALL ppm_write(ppm_rank,'ppm_map_field_pop_3d',mesg,info)
      ENDIF
      ppm_nrecvbuffer = ppm_nrecvbuffer - Mdata*bdim

      ibuffer = ppm_nrecvbuffer 
      IF (ppm_debug .GT. 1) THEN
          WRITE(mesg,'(A,I9)') 'ibuffer = ',ibuffer
          CALL ppm_write(ppm_rank,'ppm_map_field_pop_3d',mesg,info)
      ENDIF

      !-------------------------------------------------------------------------
      !  Decrement the pointer into the send buffer to allow reuse by
      !  multiple sequential push-send-pop cycles.
      !-------------------------------------------------------------------------
      Mdata = 0
      DO i=1,ppm_nsendlist
         !----------------------------------------------------------------------
         !  access mesh blocks belonging to the i-th processor in the 
         !  sendlist
         !----------------------------------------------------------------------
         DO j=ppm_psendbuffer(i),ppm_psendbuffer(i+1)-1
            !-------------------------------------------------------------------
            !  Get the number of mesh points in this block
            !-------------------------------------------------------------------
            Mdata = Mdata + (ppm_mesh_isendblksize(1,j)*         &
     &         ppm_mesh_isendblksize(2,j)*ppm_mesh_isendblksize(3,j))
         ENDDO
      ENDDO
      ppm_nsendbuffer = ppm_nsendbuffer - ppm_buffer_dim(ppm_buffer_set)*Mdata

      !-------------------------------------------------------------------------
      !  Debug output
      !-------------------------------------------------------------------------
      IF (ppm_debug .GT. 0) THEN
          IF (rtype .EQ. ppm_param_pop_replace) THEN
              CALL ppm_write(ppm_rank,'ppm_map_field_pop_3d',     &
     &           'Replacing current field values',info)
          ELSEIF(rtype .EQ. ppm_param_pop_add) THEN
              CALL ppm_write(ppm_rank,'ppm_map_field_pop_3d',     &
     &           'Adding to current field values',info)
          ENDIF
      ENDIF

      !-------------------------------------------------------------------------
      !  Loop over the received mesh blocks
      !-------------------------------------------------------------------------
      !-------------------------------------------------------------------------
      !  DOUBLE PRECISION
      !-------------------------------------------------------------------------
      IF (ppm_kind .EQ. ppm_kind_double) THEN
         DO i=1,ppm_nrecvlist
            !-------------------------------------------------------------------
            !  Access mesh blocks belonging to the i-th processor in the
            !  recvlist
            !-------------------------------------------------------------------
            DO j=ppm_precvbuffer(i),ppm_precvbuffer(i+1)-1
               !----------------------------------------------------------------
               !  Get the sub ID for this mesh block
               !----------------------------------------------------------------
               jsub = ppm_mesh_irecvtosub(j)
               !----------------------------------------------------------------
               !  Translate to local sub ID for storing the data
               !----------------------------------------------------------------
               isub = invsublist(jsub)
               !----------------------------------------------------------------
               !  Mesh offset for this sub
               !----------------------------------------------------------------
               mofs(1) = ppm_cart_mesh(tomesh,totopo)%istart(1,jsub)-1
               mofs(2) = ppm_cart_mesh(tomesh,totopo)%istart(2,jsub)-1
               mofs(3) = ppm_cart_mesh(tomesh,totopo)%istart(3,jsub)-1
               !----------------------------------------------------------------
               !  Get boundaries of mesh block to be received in local sub
               !  coordinates
               !----------------------------------------------------------------
               xlo = ppm_mesh_irecvblkstart(1,j)-mofs(1)
               ylo = ppm_mesh_irecvblkstart(2,j)-mofs(2)
               zlo = ppm_mesh_irecvblkstart(3,j)-mofs(3)
               xhi = xlo+ppm_mesh_irecvblksize(1,j)-1
               yhi = ylo+ppm_mesh_irecvblksize(2,j)-1
               zhi = zlo+ppm_mesh_irecvblksize(3,j)-1
               IF (ppm_debug .GT. 1) THEN
                   WRITE(mesg,'(A,3I4)') 'start: ',             &
     &                 ppm_mesh_irecvblkstart(1,j),             &
     &                 ppm_mesh_irecvblkstart(2,j),ppm_mesh_irecvblkstart(3,j)
                   CALL ppm_write(ppm_rank,'ppm_map_field_pop_3d',mesg,info)
                   WRITE(mesg,'(A,3I4)') 'size: ',             &
     &                 ppm_mesh_irecvblksize(1,j),             &
     &                 ppm_mesh_irecvblksize(2,j),ppm_mesh_irecvblksize(3,j)
                   CALL ppm_write(ppm_rank,'ppm_map_field_pop_3d',mesg,info)
                   WRITE(mesg,'(A,3I4)') 'mesh offset: ',mofs(1),mofs(2),mofs(3)
                   CALL ppm_write(ppm_rank,'ppm_map_field_pop_3d',mesg,info)
                   WRITE(mesg,'(A,2I4)') 'xlo, xhi: ',xlo,xhi
                   CALL ppm_write(ppm_rank,'ppm_map_field_pop_3d',mesg,info)
                   WRITE(mesg,'(A,2I4)') 'ylo, yhi: ',ylo,yhi
                   CALL ppm_write(ppm_rank,'ppm_map_field_pop_3d',mesg,info)
                   WRITE(mesg,'(A,2I4)') 'zlo, zhi: ',zlo,zhi
                   CALL ppm_write(ppm_rank,'ppm_map_field_pop_3d',mesg,info)
                   WRITE(mesg,'(A,I1)') 'buffer dim: ',edim
                   CALL ppm_write(ppm_rank,'ppm_map_field_pop_3d',mesg,info)
#if   __DIM == __VFIELD
                   WRITE(mesg,'(A,3I4)') 'SIZE(fdata): ',SIZE(fdata,2),  &
     &                 SIZE(fdata,3),SIZE(fdata,4)
                   CALL ppm_write(ppm_rank,'ppm_map_field_pop_3d',mesg,info)
#elif __DIM == __SFIELD
                   WRITE(mesg,'(A,3I4)') 'SIZE(fdata): ',SIZE(fdata,1),  &
     &                 SIZE(fdata,2),SIZE(fdata,3)
                   CALL ppm_write(ppm_rank,'ppm_map_field_pop_3d',mesg,info)
#endif
               ENDIF
#if   __DIM == __VFIELD
               !----------------------------------------------------------------
               !  Without mask: This will vectorize
               !----------------------------------------------------------------
               IF (.NOT.PRESENT(mask)) THEN
                  !-------------------------------------------------------------
                  !  Unrolled for edim=1
                  !-------------------------------------------------------------
                  IF (edim .EQ. 1) THEN
                     IF (rtype .EQ. ppm_param_pop_replace) THEN
                        DO kmesh=zlo,zhi
                          DO jmesh=ylo,yhi
                           DO imesh=xlo,xhi
#if    __KIND == __SINGLE_PRECISION_COMPLEX | __KIND == __DOUBLE_PRECISION_COMPLEX
                              ibuffer = ibuffer + 2
#else
                              ibuffer = ibuffer + 1
#endif
#if    __KIND == __SINGLE_PRECISION
                              fdata(1,imesh,jmesh,kmesh,isub) =     &
     &                           REAL(ppm_recvbufferd(ibuffer),ppm_kind_single)
#elif  __KIND == __DOUBLE_PRECISION
                              fdata(1,imesh,jmesh,kmesh,isub) =     &
     &                           ppm_recvbufferd(ibuffer)
#elif  __KIND == __SINGLE_PRECISION_COMPLEX
                              fdata(1,imesh,jmesh,kmesh,isub) =     &
     &                           CMPLX(ppm_recvbufferd(ibuffer-1),  &
     &                           ppm_recvbufferd(ibuffer),ppm_kind_single)
#elif  __KIND == __DOUBLE_PRECISION_COMPLEX
                              fdata(1,imesh,jmesh,kmesh,isub) =     &
     &                           CMPLX(ppm_recvbufferd(ibuffer-1),  &
     &                           ppm_recvbufferd(ibuffer),ppm_kind_double)
#elif  __KIND == __INTEGER
                              fdata(1,imesh,jmesh,kmesh,isub) =     &
     &                           INT(ppm_recvbufferd(ibuffer))
#elif  __KIND == __LOGICAL
                              IF (ppm_recvbufferd(ibuffer) .GT.     &
     &                           (1.0_ppm_kind_double-ppm_myepsd)) THEN
                                 fdata(1,imesh,jmesh,kmesh,isub) = .TRUE.
                              ELSE
                                 fdata(1,imesh,jmesh,kmesh,isub) = .FALSE.
                              ENDIF 
#endif
                           ENDDO
                          ENDDO
                        ENDDO
                     ELSEIF (rtype .EQ. ppm_param_pop_add) THEN
                        DO kmesh=zlo,zhi
                          DO jmesh=ylo,yhi
                           DO imesh=xlo,xhi
#if    __KIND == __SINGLE_PRECISION_COMPLEX | __KIND == __DOUBLE_PRECISION_COMPLEX
                              ibuffer = ibuffer + 2
#else
                              ibuffer = ibuffer + 1
#endif
#if    __KIND == __SINGLE_PRECISION
                              fdata(1,imesh,jmesh,kmesh,isub) =   &
     &                           fdata(1,imesh,jmesh,kmesh,isub)+ &
     &                           REAL(ppm_recvbufferd(ibuffer),ppm_kind_single)
#elif  __KIND == __DOUBLE_PRECISION
                              fdata(1,imesh,jmesh,kmesh,isub) =   &
     &                           fdata(1,imesh,jmesh,kmesh,isub)+ &
     &                           ppm_recvbufferd(ibuffer)
#elif  __KIND == __SINGLE_PRECISION_COMPLEX
                              fdata(1,imesh,jmesh,kmesh,isub) =   &
     &                           fdata(1,imesh,jmesh,kmesh,isub)+ &
     &                           CMPLX(ppm_recvbufferd(ibuffer-1),  &
     &                           ppm_recvbufferd(ibuffer),ppm_kind_single)
#elif  __KIND == __DOUBLE_PRECISION_COMPLEX
                              fdata(1,imesh,jmesh,kmesh,isub) =   &
     &                           fdata(1,imesh,jmesh,kmesh,isub)+ &
     &                           CMPLX(ppm_recvbufferd(ibuffer-1),  &
     &                           ppm_recvbufferd(ibuffer),ppm_kind_double)
#elif  __KIND == __INTEGER
                              fdata(1,imesh,jmesh,kmesh,isub) =   &
     &                           fdata(1,imesh,jmesh,kmesh,isub)+ &
     &                           INT(ppm_recvbufferd(ibuffer))
#elif  __KIND == __LOGICAL
                              IF (ppm_recvbufferd(ibuffer) .GT.     &
     &                           (1.0_ppm_kind_double-ppm_myepsd)) THEN
                                 fdata(1,imesh,jmesh,kmesh,isub) =    &
     &                            (fdata(1,imesh,jmesh,kmesh,isub) .AND..TRUE.)
                              ELSE
                                 fdata(1,imesh,jmesh,kmesh,isub) =    &
     &                            (fdata(1,imesh,jmesh,kmesh,isub) .AND..FALSE.)
                              ENDIF 
#endif
                           ENDDO
                          ENDDO
                        ENDDO
                     ENDIF
                  !-------------------------------------------------------------
                  !  Unrolled for edim=2
                  !-------------------------------------------------------------
                  ELSEIF (edim .EQ. 2) THEN
                     IF (rtype .EQ. ppm_param_pop_replace) THEN
                        DO kmesh=zlo,zhi
                          DO jmesh=ylo,yhi
                           DO imesh=xlo,xhi
#if    __KIND == __SINGLE_PRECISION_COMPLEX | __KIND == __DOUBLE_PRECISION_COMPLEX
                              ibuffer = ibuffer + 2
#else
                              ibuffer = ibuffer + 1
#endif
#if    __KIND == __SINGLE_PRECISION
                              fdata(1,imesh,jmesh,kmesh,isub) =     &
     &                           REAL(ppm_recvbufferd(ibuffer),ppm_kind_single)
                              ibuffer = ibuffer + 1
                              fdata(2,imesh,jmesh,kmesh,isub) =     &
     &                           REAL(ppm_recvbufferd(ibuffer),ppm_kind_single)
#elif  __KIND == __DOUBLE_PRECISION
                              fdata(1,imesh,jmesh,kmesh,isub) =     &
     &                           ppm_recvbufferd(ibuffer)
                              ibuffer = ibuffer + 1
                              fdata(2,imesh,jmesh,kmesh,isub) =     &
     &                           ppm_recvbufferd(ibuffer)
#elif  __KIND == __SINGLE_PRECISION_COMPLEX
                              fdata(1,imesh,jmesh,kmesh,isub) =     &
     &                           CMPLX(ppm_recvbufferd(ibuffer-1),  &
     &                           ppm_recvbufferd(ibuffer),ppm_kind_single)
                              ibuffer = ibuffer + 2
                              fdata(2,imesh,jmesh,kmesh,isub) =     &
     &                           CMPLX(ppm_recvbufferd(ibuffer-1),  &
     &                           ppm_recvbufferd(ibuffer),ppm_kind_single)
#elif  __KIND == __DOUBLE_PRECISION_COMPLEX
                              fdata(1,imesh,jmesh,kmesh,isub) =     &
     &                           CMPLX(ppm_recvbufferd(ibuffer-1),  &
     &                           ppm_recvbufferd(ibuffer),ppm_kind_double)
                              ibuffer = ibuffer + 2
                              fdata(2,imesh,jmesh,kmesh,isub) =     &
     &                           CMPLX(ppm_recvbufferd(ibuffer-1),  &
     &                           ppm_recvbufferd(ibuffer),ppm_kind_double)
#elif  __KIND == __INTEGER
                              fdata(1,imesh,jmesh,kmesh,isub) =     &
     &                           INT(ppm_recvbufferd(ibuffer))
                              ibuffer = ibuffer + 1
                              fdata(2,imesh,jmesh,kmesh,isub) =     &
     &                           INT(ppm_recvbufferd(ibuffer))
#elif  __KIND == __LOGICAL
                              IF (ppm_recvbufferd(ibuffer) .GT.     &
     &                           (1.0_ppm_kind_double-ppm_myepsd)) THEN
                                 fdata(1,imesh,jmesh,kmesh,isub) = .TRUE.
                              ELSE
                                 fdata(2,imesh,jmesh,kmesh,isub) = .FALSE.
                              ENDIF 
                              ibuffer = ibuffer + 1
                              IF (ppm_recvbufferd(ibuffer) .GT.     &
     &                           (1.0_ppm_kind_double-ppm_myepsd)) THEN
                                 fdata(1,imesh,jmesh,kmesh,isub) = .TRUE.
                              ELSE
                                 fdata(2,imesh,jmesh,kmesh,isub) = .FALSE.
                              ENDIF 
#endif
                           ENDDO
                          ENDDO
                        ENDDO
                     ELSEIF (rtype .EQ. ppm_param_pop_add) THEN
                        DO kmesh=zlo,zhi
                          DO jmesh=ylo,yhi
                           DO imesh=xlo,xhi
#if    __KIND == __SINGLE_PRECISION_COMPLEX | __KIND == __DOUBLE_PRECISION_COMPLEX
                              ibuffer = ibuffer + 2
#else
                              ibuffer = ibuffer + 1
#endif
#if    __KIND == __SINGLE_PRECISION
                              fdata(1,imesh,jmesh,kmesh,isub) =   &
     &                           fdata(1,imesh,jmesh,kmesh,isub)+ &
     &                           REAL(ppm_recvbufferd(ibuffer),ppm_kind_single)
                              ibuffer = ibuffer + 1
                              fdata(2,imesh,jmesh,kmesh,isub) =   &
     &                           fdata(2,imesh,jmesh,kmesh,isub)+ &
     &                           REAL(ppm_recvbufferd(ibuffer),ppm_kind_single)
#elif  __KIND == __DOUBLE_PRECISION
                              fdata(1,imesh,jmesh,kmesh,isub) =   &
     &                           fdata(1,imesh,jmesh,kmesh,isub)+ &
     &                           ppm_recvbufferd(ibuffer)
                              ibuffer = ibuffer + 1
                              fdata(2,imesh,jmesh,kmesh,isub) =   &
     &                           fdata(2,imesh,jmesh,kmesh,isub)+ &
     &                           ppm_recvbufferd(ibuffer)
#elif  __KIND == __SINGLE_PRECISION_COMPLEX
                              fdata(1,imesh,jmesh,kmesh,isub) =   &
     &                           fdata(1,imesh,jmesh,kmesh,isub)+ &
     &                           CMPLX(ppm_recvbufferd(ibuffer-1),  &
     &                           ppm_recvbufferd(ibuffer),ppm_kind_single)
                              ibuffer = ibuffer + 2
                              fdata(2,imesh,jmesh,kmesh,isub) =   &
     &                           fdata(2,imesh,jmesh,kmesh,isub)+ &
     &                           CMPLX(ppm_recvbufferd(ibuffer-1),  &
     &                           ppm_recvbufferd(ibuffer),ppm_kind_single)
#elif  __KIND == __DOUBLE_PRECISION_COMPLEX
                              fdata(1,imesh,jmesh,kmesh,isub) =   &
     &                           fdata(1,imesh,jmesh,kmesh,isub)+ &
     &                           CMPLX(ppm_recvbufferd(ibuffer-1),  &
     &                           ppm_recvbufferd(ibuffer),ppm_kind_double)
                              ibuffer = ibuffer + 2
                              fdata(2,imesh,jmesh,kmesh,isub) =   &
     &                           fdata(2,imesh,jmesh,kmesh,isub)+ &
     &                           CMPLX(ppm_recvbufferd(ibuffer-1),  &
     &                           ppm_recvbufferd(ibuffer),ppm_kind_double)
#elif  __KIND == __INTEGER
                              fdata(1,imesh,jmesh,kmesh,isub) =   &
     &                           fdata(1,imesh,jmesh,kmesh,isub)+ &
     &                           INT(ppm_recvbufferd(ibuffer))
                              ibuffer = ibuffer + 1
                              fdata(2,imesh,jmesh,kmesh,isub) =   &
     &                           fdata(2,imesh,jmesh,kmesh,isub)+ &
     &                           INT(ppm_recvbufferd(ibuffer))
#elif  __KIND == __LOGICAL
                              IF (ppm_recvbufferd(ibuffer) .GT.     &
     &                           (1.0_ppm_kind_double-ppm_myepsd)) THEN
                                 fdata(1,imesh,jmesh,kmesh,isub) =    &
     &                            (fdata(1,imesh,jmesh,kmesh,isub) .AND..TRUE.)
                              ELSE
                                 fdata(1,imesh,jmesh,kmesh,isub) =    &
     &                            (fdata(1,imesh,jmesh,kmesh,isub) .AND..FALSE.)
                              ENDIF 
                              ibuffer = ibuffer + 1
                              IF (ppm_recvbufferd(ibuffer) .GT.     &
     &                           (1.0_ppm_kind_double-ppm_myepsd)) THEN
                                 fdata(2,imesh,jmesh,kmesh,isub) =    &
     &                            (fdata(2,imesh,jmesh,kmesh,isub) .AND..TRUE.)
                              ELSE
                                 fdata(2,imesh,jmesh,kmesh,isub) =    &
     &                            (fdata(2,imesh,jmesh,kmesh,isub) .AND..FALSE.)
                              ENDIF 
#endif
                           ENDDO
                          ENDDO
                        ENDDO
                     ENDIF
                  !-------------------------------------------------------------
                  !  Unrolled for edim=3
                  !-------------------------------------------------------------
                  ELSEIF (edim .EQ. 3) THEN
                     IF (rtype .EQ. ppm_param_pop_replace) THEN
                        DO kmesh=zlo,zhi
                          DO jmesh=ylo,yhi
                           DO imesh=xlo,xhi
#if    __KIND == __SINGLE_PRECISION_COMPLEX | __KIND == __DOUBLE_PRECISION_COMPLEX
                              ibuffer = ibuffer + 2
#else
                              ibuffer = ibuffer + 1
#endif
#if    __KIND == __SINGLE_PRECISION
                              fdata(1,imesh,jmesh,kmesh,isub) =     &
     &                           REAL(ppm_recvbufferd(ibuffer),ppm_kind_single)
                              ibuffer = ibuffer + 1
                              fdata(2,imesh,jmesh,kmesh,isub) =     &
     &                           REAL(ppm_recvbufferd(ibuffer),ppm_kind_single)
                              ibuffer = ibuffer + 1
                              fdata(3,imesh,jmesh,kmesh,isub) =     &
     &                           REAL(ppm_recvbufferd(ibuffer),ppm_kind_single)
#elif  __KIND == __DOUBLE_PRECISION
                              fdata(1,imesh,jmesh,kmesh,isub) =     &
     &                           ppm_recvbufferd(ibuffer)
                              ibuffer = ibuffer + 1
                              fdata(2,imesh,jmesh,kmesh,isub) =     &
     &                           ppm_recvbufferd(ibuffer)
                              ibuffer = ibuffer + 1
                              fdata(3,imesh,jmesh,kmesh,isub) =     &
     &                           ppm_recvbufferd(ibuffer)
#elif  __KIND == __SINGLE_PRECISION_COMPLEX
                              fdata(1,imesh,jmesh,kmesh,isub) =     &
     &                           CMPLX(ppm_recvbufferd(ibuffer-1),  &
     &                           ppm_recvbufferd(ibuffer),ppm_kind_single)
                              ibuffer = ibuffer + 2
                              fdata(2,imesh,jmesh,kmesh,isub) =     &
     &                           CMPLX(ppm_recvbufferd(ibuffer-1),  &
     &                           ppm_recvbufferd(ibuffer),ppm_kind_single)
                              ibuffer = ibuffer + 2
                              fdata(3,imesh,jmesh,kmesh,isub) =     &
     &                           CMPLX(ppm_recvbufferd(ibuffer-1),  &
     &                           ppm_recvbufferd(ibuffer),ppm_kind_single)
#elif  __KIND == __DOUBLE_PRECISION_COMPLEX
                              fdata(1,imesh,jmesh,kmesh,isub) =     &
     &                           CMPLX(ppm_recvbufferd(ibuffer-1),  &
     &                           ppm_recvbufferd(ibuffer),ppm_kind_double)
                              ibuffer = ibuffer + 2
                              fdata(2,imesh,jmesh,kmesh,isub) =     &
     &                           CMPLX(ppm_recvbufferd(ibuffer-1),  &
     &                           ppm_recvbufferd(ibuffer),ppm_kind_double)
                              ibuffer = ibuffer + 2
                              fdata(3,imesh,jmesh,kmesh,isub) =     &
     &                           CMPLX(ppm_recvbufferd(ibuffer-1),  &
     &                           ppm_recvbufferd(ibuffer),ppm_kind_double)
#elif  __KIND == __INTEGER
                              fdata(1,imesh,jmesh,kmesh,isub) =     &
     &                           INT(ppm_recvbufferd(ibuffer))
                              ibuffer = ibuffer + 1
                              fdata(2,imesh,jmesh,kmesh,isub) =     &
     &                           INT(ppm_recvbufferd(ibuffer))
                              ibuffer = ibuffer + 1
                              fdata(3,imesh,jmesh,kmesh,isub) =     &
     &                           INT(ppm_recvbufferd(ibuffer))
#elif  __KIND == __LOGICAL
                              IF (ppm_recvbufferd(ibuffer) .GT.     &
     &                           (1.0_ppm_kind_double-ppm_myepsd)) THEN
                                 fdata(1,imesh,jmesh,kmesh,isub) = .TRUE.
                              ELSE
                                 fdata(1,imesh,jmesh,kmesh,isub) = .FALSE.
                              ENDIF 
                              ibuffer = ibuffer + 1
                              IF (ppm_recvbufferd(ibuffer) .GT.     &
     &                           (1.0_ppm_kind_double-ppm_myepsd)) THEN
                                 fdata(2,imesh,jmesh,kmesh,isub) = .TRUE.
                              ELSE
                                 fdata(2,imesh,jmesh,kmesh,isub) = .FALSE.
                              ENDIF 
                              ibuffer = ibuffer + 1
                              IF (ppm_recvbufferd(ibuffer) .GT.     &
     &                           (1.0_ppm_kind_double-ppm_myepsd)) THEN
                                 fdata(3,imesh,jmesh,kmesh,isub) = .TRUE.
                              ELSE
                                 fdata(3,imesh,jmesh,kmesh,isub) = .FALSE.
                              ENDIF 
#endif
                           ENDDO
                          ENDDO
                        ENDDO
                     ELSEIF (rtype .EQ. ppm_param_pop_add) THEN
                        DO kmesh=zlo,zhi
                          DO jmesh=ylo,yhi
                           DO imesh=xlo,xhi
#if    __KIND == __SINGLE_PRECISION_COMPLEX | __KIND == __DOUBLE_PRECISION_COMPLEX
                              ibuffer = ibuffer + 2
#else
                              ibuffer = ibuffer + 1
#endif
#if    __KIND == __SINGLE_PRECISION
                              fdata(1,imesh,jmesh,kmesh,isub) =   &
     &                           fdata(1,imesh,jmesh,kmesh,isub)+ &
     &                           REAL(ppm_recvbufferd(ibuffer),ppm_kind_single)
                              ibuffer = ibuffer + 1
                              fdata(2,imesh,jmesh,kmesh,isub) =   &
     &                           fdata(2,imesh,jmesh,kmesh,isub)+ &
     &                           REAL(ppm_recvbufferd(ibuffer),ppm_kind_single)
                              ibuffer = ibuffer + 1
                              fdata(3,imesh,jmesh,kmesh,isub) =   &
     &                           fdata(3,imesh,jmesh,kmesh,isub)+ &
     &                           REAL(ppm_recvbufferd(ibuffer),ppm_kind_single)
#elif  __KIND == __DOUBLE_PRECISION
                              fdata(1,imesh,jmesh,kmesh,isub) =   &
     &                           fdata(1,imesh,jmesh,kmesh,isub)+ &
     &                           ppm_recvbufferd(ibuffer)
                              ibuffer = ibuffer + 1
                              fdata(2,imesh,jmesh,kmesh,isub) =   &
     &                           fdata(2,imesh,jmesh,kmesh,isub)+ &
     &                           ppm_recvbufferd(ibuffer)
                              ibuffer = ibuffer + 1
                              fdata(3,imesh,jmesh,kmesh,isub) =   &
     &                           fdata(3,imesh,jmesh,kmesh,isub)+ &
     &                           ppm_recvbufferd(ibuffer)
#elif  __KIND == __SINGLE_PRECISION_COMPLEX
                              fdata(1,imesh,jmesh,kmesh,isub) =   &
     &                           fdata(1,imesh,jmesh,kmesh,isub)+ &
     &                           CMPLX(ppm_recvbufferd(ibuffer-1),  &
     &                           ppm_recvbufferd(ibuffer),ppm_kind_single)
                              ibuffer = ibuffer + 2
                              fdata(2,imesh,jmesh,kmesh,isub) =   &
     &                           fdata(2,imesh,jmesh,kmesh,isub)+ &
     &                           CMPLX(ppm_recvbufferd(ibuffer-1),  &
     &                           ppm_recvbufferd(ibuffer),ppm_kind_single)
                              ibuffer = ibuffer + 2
                              fdata(3,imesh,jmesh,kmesh,isub) =   &
     &                           fdata(3,imesh,jmesh,kmesh,isub)+ &
     &                           CMPLX(ppm_recvbufferd(ibuffer-1),  &
     &                           ppm_recvbufferd(ibuffer),ppm_kind_single)
#elif  __KIND == __DOUBLE_PRECISION_COMPLEX
                              fdata(1,imesh,jmesh,kmesh,isub) =   &
     &                           fdata(1,imesh,jmesh,kmesh,isub)+ &
     &                           CMPLX(ppm_recvbufferd(ibuffer-1),  &
     &                           ppm_recvbufferd(ibuffer),ppm_kind_double)
                              ibuffer = ibuffer + 2
                              fdata(2,imesh,jmesh,kmesh,isub) =   &
     &                           fdata(2,imesh,jmesh,kmesh,isub)+ &
     &                           CMPLX(ppm_recvbufferd(ibuffer-1),  &
     &                           ppm_recvbufferd(ibuffer),ppm_kind_double)
                              ibuffer = ibuffer + 2
                              fdata(3,imesh,jmesh,kmesh,isub) =   &
     &                           fdata(3,imesh,jmesh,kmesh,isub)+ &
     &                           CMPLX(ppm_recvbufferd(ibuffer-1),  &
     &                           ppm_recvbufferd(ibuffer),ppm_kind_double)
#elif  __KIND == __INTEGER
                              fdata(1,imesh,jmesh,kmesh,isub) =   &
     &                           fdata(1,imesh,jmesh,kmesh,isub)+ &
     &                           INT(ppm_recvbufferd(ibuffer))
                              ibuffer = ibuffer + 1
                              fdata(2,imesh,jmesh,kmesh,isub) =   &
     &                           fdata(2,imesh,jmesh,kmesh,isub)+ &
     &                           INT(ppm_recvbufferd(ibuffer))
                              ibuffer = ibuffer + 1
                              fdata(3,imesh,jmesh,kmesh,isub) =   &
     &                           fdata(3,imesh,jmesh,kmesh,isub)+ &
     &                           INT(ppm_recvbufferd(ibuffer))
#elif  __KIND == __LOGICAL
                              IF (ppm_recvbufferd(ibuffer) .GT.     &
     &                           (1.0_ppm_kind_double-ppm_myepsd)) THEN
                                 fdata(1,imesh,jmesh,kmesh,isub) =    &
     &                            (fdata(1,imesh,jmesh,kmesh,isub) .AND..TRUE.)
                              ELSE
                                 fdata(1,imesh,jmesh,kmesh,isub) =    &
     &                            (fdata(1,imesh,jmesh,kmesh,isub) .AND..FALSE.)
                              ENDIF 
                              ibuffer = ibuffer + 1
                              IF (ppm_recvbufferd(ibuffer) .GT.     &
     &                           (1.0_ppm_kind_double-ppm_myepsd)) THEN
                                 fdata(2,imesh,jmesh,kmesh,isub) =    &
     &                            (fdata(2,imesh,jmesh,kmesh,isub) .AND..TRUE.)
                              ELSE
                                 fdata(2,imesh,jmesh,kmesh,isub) =    &
     &                            (fdata(2,imesh,jmesh,kmesh,isub) .AND..FALSE.)
                              ENDIF 
                              ibuffer = ibuffer + 1
                              IF (ppm_recvbufferd(ibuffer) .GT.     &
     &                           (1.0_ppm_kind_double-ppm_myepsd)) THEN
                                 fdata(3,imesh,jmesh,kmesh,isub) =    &
     &                            (fdata(3,imesh,jmesh,kmesh,isub) .AND..TRUE.)
                              ELSE
                                 fdata(3,imesh,jmesh,kmesh,isub) =    &
     &                            (fdata(3,imesh,jmesh,kmesh,isub) .AND..FALSE.)
                              ENDIF 
#endif
                           ENDDO
                          ENDDO
                        ENDDO
                     ENDIF
                  !-------------------------------------------------------------
                  !  Unrolled for edim=4
                  !-------------------------------------------------------------
                  ELSEIF (edim .EQ. 4) THEN
                     IF (rtype .EQ. ppm_param_pop_replace) THEN
                        DO kmesh=zlo,zhi
                          DO jmesh=ylo,yhi
                           DO imesh=xlo,xhi
#if    __KIND == __SINGLE_PRECISION_COMPLEX | __KIND == __DOUBLE_PRECISION_COMPLEX
                              ibuffer = ibuffer + 2
#else
                              ibuffer = ibuffer + 1
#endif
#if    __KIND == __SINGLE_PRECISION
                              fdata(1,imesh,jmesh,kmesh,isub) =     &
     &                           REAL(ppm_recvbufferd(ibuffer),ppm_kind_single)
                              ibuffer = ibuffer + 1
                              fdata(2,imesh,jmesh,kmesh,isub) =     &
     &                           REAL(ppm_recvbufferd(ibuffer),ppm_kind_single)
                              ibuffer = ibuffer + 1
                              fdata(3,imesh,jmesh,kmesh,isub) =     &
     &                           REAL(ppm_recvbufferd(ibuffer),ppm_kind_single)
                              ibuffer = ibuffer + 1
                              fdata(4,imesh,jmesh,kmesh,isub) =     &
     &                           REAL(ppm_recvbufferd(ibuffer),ppm_kind_single)
#elif  __KIND == __DOUBLE_PRECISION
                              fdata(1,imesh,jmesh,kmesh,isub) =     &
     &                           ppm_recvbufferd(ibuffer)
                              ibuffer = ibuffer + 1
                              fdata(2,imesh,jmesh,kmesh,isub) =     &
     &                           ppm_recvbufferd(ibuffer)
                              ibuffer = ibuffer + 1
                              fdata(3,imesh,jmesh,kmesh,isub) =     &
     &                           ppm_recvbufferd(ibuffer)
                              ibuffer = ibuffer + 1
                              fdata(4,imesh,jmesh,kmesh,isub) =     &
     &                           ppm_recvbufferd(ibuffer)
#elif  __KIND == __SINGLE_PRECISION_COMPLEX
                              fdata(1,imesh,jmesh,kmesh,isub) =     &
     &                           CMPLX(ppm_recvbufferd(ibuffer-1),  &
     &                           ppm_recvbufferd(ibuffer),ppm_kind_single)
                              ibuffer = ibuffer + 2
                              fdata(2,imesh,jmesh,kmesh,isub) =     &
     &                           CMPLX(ppm_recvbufferd(ibuffer-1),  &
     &                           ppm_recvbufferd(ibuffer),ppm_kind_single)
                              ibuffer = ibuffer + 2
                              fdata(3,imesh,jmesh,kmesh,isub) =     &
     &                           CMPLX(ppm_recvbufferd(ibuffer-1),  &
     &                           ppm_recvbufferd(ibuffer),ppm_kind_single)
                              ibuffer = ibuffer + 2
                              fdata(4,imesh,jmesh,kmesh,isub) =     &
     &                           CMPLX(ppm_recvbufferd(ibuffer-1),  &
     &                           ppm_recvbufferd(ibuffer),ppm_kind_single)
#elif  __KIND == __DOUBLE_PRECISION_COMPLEX
                              fdata(1,imesh,jmesh,kmesh,isub) =     &
     &                           CMPLX(ppm_recvbufferd(ibuffer-1),  &
     &                           ppm_recvbufferd(ibuffer),ppm_kind_double)
                              ibuffer = ibuffer + 2
                              fdata(2,imesh,jmesh,kmesh,isub) =     &
     &                           CMPLX(ppm_recvbufferd(ibuffer-1),  &
     &                           ppm_recvbufferd(ibuffer),ppm_kind_double)
                              ibuffer = ibuffer + 2
                              fdata(3,imesh,jmesh,kmesh,isub) =     &
     &                           CMPLX(ppm_recvbufferd(ibuffer-1),  &
     &                           ppm_recvbufferd(ibuffer),ppm_kind_double)
                              ibuffer = ibuffer + 2
                              fdata(4,imesh,jmesh,kmesh,isub) =     &
     &                           CMPLX(ppm_recvbufferd(ibuffer-1),  &
     &                           ppm_recvbufferd(ibuffer),ppm_kind_double)
#elif  __KIND == __INTEGER
                              fdata(1,imesh,jmesh,kmesh,isub) =     &
     &                           INT(ppm_recvbufferd(ibuffer))
                              ibuffer = ibuffer + 1
                              fdata(2,imesh,jmesh,kmesh,isub) =     &
     &                           INT(ppm_recvbufferd(ibuffer))
                              ibuffer = ibuffer + 1
                              fdata(3,imesh,jmesh,kmesh,isub) =     &
     &                           INT(ppm_recvbufferd(ibuffer))
                              ibuffer = ibuffer + 1
                              fdata(4,imesh,jmesh,kmesh,isub) =     &
     &                           INT(ppm_recvbufferd(ibuffer))
#elif  __KIND == __LOGICAL
                              IF (ppm_recvbufferd(ibuffer) .GT.     &
     &                           (1.0_ppm_kind_double-ppm_myepsd)) THEN
                                 fdata(1,imesh,jmesh,kmesh,isub) = .TRUE.
                              ELSE
                                 fdata(1,imesh,jmesh,kmesh,isub) = .FALSE.
                              ENDIF 
                              ibuffer = ibuffer + 1
                              IF (ppm_recvbufferd(ibuffer) .GT.     &
     &                           (1.0_ppm_kind_double-ppm_myepsd)) THEN
                                 fdata(2,imesh,jmesh,kmesh,isub) = .TRUE.
                              ELSE
                                 fdata(2,imesh,jmesh,kmesh,isub) = .FALSE.
                              ENDIF 
                              ibuffer = ibuffer + 1
                              IF (ppm_recvbufferd(ibuffer) .GT.     &
     &                           (1.0_ppm_kind_double-ppm_myepsd)) THEN
                                 fdata(3,imesh,jmesh,kmesh,isub) = .TRUE.
                              ELSE
                                 fdata(3,imesh,jmesh,kmesh,isub) = .FALSE.
                              ENDIF 
                              ibuffer = ibuffer + 1
                              IF (ppm_recvbufferd(ibuffer) .GT.     &
     &                           (1.0_ppm_kind_double-ppm_myepsd)) THEN
                                 fdata(4,imesh,jmesh,kmesh,isub) = .TRUE.
                              ELSE
                                 fdata(4,imesh,jmesh,kmesh,isub) = .FALSE.
                              ENDIF 
#endif
                           ENDDO
                          ENDDO
                        ENDDO
                     ELSEIF (rtype .EQ. ppm_param_pop_add) THEN
                        DO kmesh=zlo,zhi
                          DO jmesh=ylo,yhi
                           DO imesh=xlo,xhi
#if    __KIND == __SINGLE_PRECISION_COMPLEX | __KIND == __DOUBLE_PRECISION_COMPLEX
                              ibuffer = ibuffer + 2
#else
                              ibuffer = ibuffer + 1
#endif
#if    __KIND == __SINGLE_PRECISION
                              fdata(1,imesh,jmesh,kmesh,isub) =   &
     &                           fdata(1,imesh,jmesh,kmesh,isub)+ &
     &                           REAL(ppm_recvbufferd(ibuffer),ppm_kind_single)
                              ibuffer = ibuffer + 1
                              fdata(2,imesh,jmesh,kmesh,isub) =   &
     &                           fdata(2,imesh,jmesh,kmesh,isub)+ &
     &                           REAL(ppm_recvbufferd(ibuffer),ppm_kind_single)
                              ibuffer = ibuffer + 1
                              fdata(3,imesh,jmesh,kmesh,isub) =   &
     &                           fdata(3,imesh,jmesh,kmesh,isub)+ &
     &                           REAL(ppm_recvbufferd(ibuffer),ppm_kind_single)
                              ibuffer = ibuffer + 1
                              fdata(4,imesh,jmesh,kmesh,isub) =   &
     &                           fdata(4,imesh,jmesh,kmesh,isub)+ &
     &                           REAL(ppm_recvbufferd(ibuffer),ppm_kind_single)
#elif  __KIND == __DOUBLE_PRECISION
                              fdata(1,imesh,jmesh,kmesh,isub) =   &
     &                           fdata(1,imesh,jmesh,kmesh,isub)+ &
     &                           ppm_recvbufferd(ibuffer)
                              ibuffer = ibuffer + 1
                              fdata(2,imesh,jmesh,kmesh,isub) =   &
     &                           fdata(2,imesh,jmesh,kmesh,isub)+ &
     &                           ppm_recvbufferd(ibuffer)
                              ibuffer = ibuffer + 1
                              fdata(3,imesh,jmesh,kmesh,isub) =   &
     &                           fdata(3,imesh,jmesh,kmesh,isub)+ &
     &                           ppm_recvbufferd(ibuffer)
                              ibuffer = ibuffer + 1
                              fdata(4,imesh,jmesh,kmesh,isub) =   &
     &                           fdata(4,imesh,jmesh,kmesh,isub)+ &
     &                           ppm_recvbufferd(ibuffer)
#elif  __KIND == __SINGLE_PRECISION_COMPLEX
                              fdata(1,imesh,jmesh,kmesh,isub) =   &
     &                           fdata(1,imesh,jmesh,kmesh,isub)+ &
     &                           CMPLX(ppm_recvbufferd(ibuffer-1),  &
     &                           ppm_recvbufferd(ibuffer),ppm_kind_single)
                              ibuffer = ibuffer + 2
                              fdata(2,imesh,jmesh,kmesh,isub) =   &
     &                           fdata(2,imesh,jmesh,kmesh,isub)+ &
     &                           CMPLX(ppm_recvbufferd(ibuffer-1),  &
     &                           ppm_recvbufferd(ibuffer),ppm_kind_single)
                              ibuffer = ibuffer + 2
                              fdata(3,imesh,jmesh,kmesh,isub) =   &
     &                           fdata(3,imesh,jmesh,kmesh,isub)+ &
     &                           CMPLX(ppm_recvbufferd(ibuffer-1),  &
     &                           ppm_recvbufferd(ibuffer),ppm_kind_single)
                              ibuffer = ibuffer + 2
                              fdata(4,imesh,jmesh,kmesh,isub) =   &
     &                           fdata(4,imesh,jmesh,kmesh,isub)+ &
     &                           CMPLX(ppm_recvbufferd(ibuffer-1),  &
     &                           ppm_recvbufferd(ibuffer),ppm_kind_single)
#elif  __KIND == __DOUBLE_PRECISION_COMPLEX
                              fdata(1,imesh,jmesh,kmesh,isub) =   &
     &                           fdata(1,imesh,jmesh,kmesh,isub)+ &
     &                           CMPLX(ppm_recvbufferd(ibuffer-1),  &
     &                           ppm_recvbufferd(ibuffer),ppm_kind_double)
                              ibuffer = ibuffer + 2
                              fdata(2,imesh,jmesh,kmesh,isub) =   &
     &                           fdata(2,imesh,jmesh,kmesh,isub)+ &
     &                           CMPLX(ppm_recvbufferd(ibuffer-1),  &
     &                           ppm_recvbufferd(ibuffer),ppm_kind_double)
                              ibuffer = ibuffer + 2
                              fdata(3,imesh,jmesh,kmesh,isub) =   &
     &                           fdata(3,imesh,jmesh,kmesh,isub)+ &
     &                           CMPLX(ppm_recvbufferd(ibuffer-1),  &
     &                           ppm_recvbufferd(ibuffer),ppm_kind_double)
                              ibuffer = ibuffer + 2
                              fdata(4,imesh,jmesh,kmesh,isub) =   &
     &                           fdata(4,imesh,jmesh,kmesh,isub)+ &
     &                           CMPLX(ppm_recvbufferd(ibuffer-1),  &
     &                           ppm_recvbufferd(ibuffer),ppm_kind_double)
#elif  __KIND == __INTEGER
                              fdata(1,imesh,jmesh,kmesh,isub) =   &
     &                           fdata(1,imesh,jmesh,kmesh,isub)+ &
     &                           INT(ppm_recvbufferd(ibuffer))
                              ibuffer = ibuffer + 1
                              fdata(2,imesh,jmesh,kmesh,isub) =   &
     &                           fdata(2,imesh,jmesh,kmesh,isub)+ &
     &                           INT(ppm_recvbufferd(ibuffer))
                              ibuffer = ibuffer + 1
                              fdata(3,imesh,jmesh,kmesh,isub) =   &
     &                           fdata(3,imesh,jmesh,kmesh,isub)+ &
     &                           INT(ppm_recvbufferd(ibuffer))
                              ibuffer = ibuffer + 1
                              fdata(4,imesh,jmesh,kmesh,isub) =   &
     &                           fdata(4,imesh,jmesh,kmesh,isub)+ &
     &                           INT(ppm_recvbufferd(ibuffer))
#elif  __KIND == __LOGICAL
                              IF (ppm_recvbufferd(ibuffer) .GT.     &
     &                           (1.0_ppm_kind_double-ppm_myepsd)) THEN
                                 fdata(1,imesh,jmesh,kmesh,isub) =    &
     &                            (fdata(1,imesh,jmesh,kmesh,isub) .AND..TRUE.)
                              ELSE
                                 fdata(1,imesh,jmesh,kmesh,isub) =    &
     &                            (fdata(1,imesh,jmesh,kmesh,isub) .AND..FALSE.)
                              ENDIF 
                              ibuffer = ibuffer + 1
                              IF (ppm_recvbufferd(ibuffer) .GT.     &
     &                           (1.0_ppm_kind_double-ppm_myepsd)) THEN
                                 fdata(2,imesh,jmesh,kmesh,isub) =    &
     &                            (fdata(2,imesh,jmesh,kmesh,isub) .AND..TRUE.)
                              ELSE
                                 fdata(2,imesh,jmesh,kmesh,isub) =    &
     &                            (fdata(2,imesh,jmesh,kmesh,isub) .AND..FALSE.)
                              ENDIF 
                              ibuffer = ibuffer + 1
                              IF (ppm_recvbufferd(ibuffer) .GT.     &
     &                           (1.0_ppm_kind_double-ppm_myepsd)) THEN
                                 fdata(3,imesh,jmesh,kmesh,isub) =    &
     &                            (fdata(3,imesh,jmesh,kmesh,isub) .AND..TRUE.)
                              ELSE
                                 fdata(3,imesh,jmesh,kmesh,isub) =    &
     &                            (fdata(3,imesh,jmesh,kmesh,isub) .AND..FALSE.)
                              ENDIF 
                              ibuffer = ibuffer + 1
                              IF (ppm_recvbufferd(ibuffer) .GT.     &
     &                           (1.0_ppm_kind_double-ppm_myepsd)) THEN
                                 fdata(4,imesh,jmesh,kmesh,isub) =    &
     &                            (fdata(4,imesh,jmesh,kmesh,isub) .AND..TRUE.)
                              ELSE
                                 fdata(4,imesh,jmesh,kmesh,isub) =    &
     &                            (fdata(4,imesh,jmesh,kmesh,isub) .AND..FALSE.)
                              ENDIF 
#endif
                           ENDDO
                          ENDDO
                        ENDDO
                     ENDIF
                  !-------------------------------------------------------------
                  !  Unrolled for edim=5
                  !-------------------------------------------------------------
                  ELSEIF (edim .EQ. 5) THEN
                     IF (rtype .EQ. ppm_param_pop_replace) THEN
                        DO kmesh=zlo,zhi
                          DO jmesh=ylo,yhi
                           DO imesh=xlo,xhi
#if    __KIND == __SINGLE_PRECISION_COMPLEX | __KIND == __DOUBLE_PRECISION_COMPLEX
                              ibuffer = ibuffer + 2
#else
                              ibuffer = ibuffer + 1
#endif
#if    __KIND == __SINGLE_PRECISION
                              fdata(1,imesh,jmesh,kmesh,isub) =     &
     &                           REAL(ppm_recvbufferd(ibuffer),ppm_kind_single)
                              ibuffer = ibuffer + 1
                              fdata(2,imesh,jmesh,kmesh,isub) =     &
     &                           REAL(ppm_recvbufferd(ibuffer),ppm_kind_single)
                              ibuffer = ibuffer + 1
                              fdata(3,imesh,jmesh,kmesh,isub) =     &
     &                           REAL(ppm_recvbufferd(ibuffer),ppm_kind_single)
                              ibuffer = ibuffer + 1
                              fdata(4,imesh,jmesh,kmesh,isub) =     &
     &                           REAL(ppm_recvbufferd(ibuffer),ppm_kind_single)
                              ibuffer = ibuffer + 1
                              fdata(5,imesh,jmesh,kmesh,isub) =     &
     &                           REAL(ppm_recvbufferd(ibuffer),ppm_kind_single)
#elif  __KIND == __DOUBLE_PRECISION
                              fdata(1,imesh,jmesh,kmesh,isub) =     &
     &                           ppm_recvbufferd(ibuffer)
                              ibuffer = ibuffer + 1
                              fdata(2,imesh,jmesh,kmesh,isub) =     &
     &                           ppm_recvbufferd(ibuffer)
                              ibuffer = ibuffer + 1
                              fdata(3,imesh,jmesh,kmesh,isub) =     &
     &                           ppm_recvbufferd(ibuffer)
                              ibuffer = ibuffer + 1
                              fdata(4,imesh,jmesh,kmesh,isub) =     &
     &                           ppm_recvbufferd(ibuffer)
                              ibuffer = ibuffer + 1
                              fdata(5,imesh,jmesh,kmesh,isub) =     &
     &                           ppm_recvbufferd(ibuffer)
#elif  __KIND == __SINGLE_PRECISION_COMPLEX
                              fdata(1,imesh,jmesh,kmesh,isub) =     &
     &                           CMPLX(ppm_recvbufferd(ibuffer-1),  &
     &                           ppm_recvbufferd(ibuffer),ppm_kind_single)
                              ibuffer = ibuffer + 2
                              fdata(2,imesh,jmesh,kmesh,isub) =     &
     &                           CMPLX(ppm_recvbufferd(ibuffer-1),  &
     &                           ppm_recvbufferd(ibuffer),ppm_kind_single)
                              ibuffer = ibuffer + 2
                              fdata(3,imesh,jmesh,kmesh,isub) =     &
     &                           CMPLX(ppm_recvbufferd(ibuffer-1),  &
     &                           ppm_recvbufferd(ibuffer),ppm_kind_single)
                              ibuffer = ibuffer + 2
                              fdata(4,imesh,jmesh,kmesh,isub) =     &
     &                           CMPLX(ppm_recvbufferd(ibuffer-1),  &
     &                           ppm_recvbufferd(ibuffer),ppm_kind_single)
                              ibuffer = ibuffer + 2
                              fdata(5,imesh,jmesh,kmesh,isub) =     &
     &                           CMPLX(ppm_recvbufferd(ibuffer-1),  &
     &                           ppm_recvbufferd(ibuffer),ppm_kind_single)
#elif  __KIND == __DOUBLE_PRECISION_COMPLEX
                              fdata(1,imesh,jmesh,kmesh,isub) =     &
     &                           CMPLX(ppm_recvbufferd(ibuffer-1),  &
     &                           ppm_recvbufferd(ibuffer),ppm_kind_double)
                              ibuffer = ibuffer + 2
                              fdata(2,imesh,jmesh,kmesh,isub) =     &
     &                           CMPLX(ppm_recvbufferd(ibuffer-1),  &
     &                           ppm_recvbufferd(ibuffer),ppm_kind_double)
                              ibuffer = ibuffer + 2
                              fdata(3,imesh,jmesh,kmesh,isub) =     &
     &                           CMPLX(ppm_recvbufferd(ibuffer-1),  &
     &                           ppm_recvbufferd(ibuffer),ppm_kind_double)
                              ibuffer = ibuffer + 2
                              fdata(4,imesh,jmesh,kmesh,isub) =     &
     &                           CMPLX(ppm_recvbufferd(ibuffer-1),  &
     &                           ppm_recvbufferd(ibuffer),ppm_kind_double)
                              ibuffer = ibuffer + 2
                              fdata(5,imesh,jmesh,kmesh,isub) =     &
     &                           CMPLX(ppm_recvbufferd(ibuffer-1),  &
     &                           ppm_recvbufferd(ibuffer),ppm_kind_double)
#elif  __KIND == __INTEGER
                              fdata(1,imesh,jmesh,kmesh,isub) =     &
     &                           INT(ppm_recvbufferd(ibuffer))
                              ibuffer = ibuffer + 1
                              fdata(2,imesh,jmesh,kmesh,isub) =     &
     &                           INT(ppm_recvbufferd(ibuffer))
                              ibuffer = ibuffer + 1
                              fdata(3,imesh,jmesh,kmesh,isub) =     &
     &                           INT(ppm_recvbufferd(ibuffer))
                              ibuffer = ibuffer + 1
                              fdata(4,imesh,jmesh,kmesh,isub) =     &
     &                           INT(ppm_recvbufferd(ibuffer))
                              ibuffer = ibuffer + 1
                              fdata(5,imesh,jmesh,kmesh,isub) =     &
     &                           INT(ppm_recvbufferd(ibuffer))
#elif  __KIND == __LOGICAL
                              IF (ppm_recvbufferd(ibuffer) .GT.     &
     &                           (1.0_ppm_kind_double-ppm_myepsd)) THEN
                                 fdata(1,imesh,jmesh,kmesh,isub) = .TRUE.
                              ELSE
                                 fdata(1,imesh,jmesh,kmesh,isub) = .FALSE.
                              ENDIF 
                              ibuffer = ibuffer + 1
                              IF (ppm_recvbufferd(ibuffer) .GT.     &
     &                           (1.0_ppm_kind_double-ppm_myepsd)) THEN
                                 fdata(2,imesh,jmesh,kmesh,isub) = .TRUE.
                              ELSE
                                 fdata(2,imesh,jmesh,kmesh,isub) = .FALSE.
                              ENDIF 
                              ibuffer = ibuffer + 1
                              IF (ppm_recvbufferd(ibuffer) .GT.     &
     &                           (1.0_ppm_kind_double-ppm_myepsd)) THEN
                                 fdata(3,imesh,jmesh,kmesh,isub) = .TRUE.
                              ELSE
                                 fdata(3,imesh,jmesh,kmesh,isub) = .FALSE.
                              ENDIF 
                              ibuffer = ibuffer + 1
                              IF (ppm_recvbufferd(ibuffer) .GT.     &
     &                           (1.0_ppm_kind_double-ppm_myepsd)) THEN
                                 fdata(4,imesh,jmesh,kmesh,isub) = .TRUE.
                              ELSE
                                 fdata(4,imesh,jmesh,kmesh,isub) = .FALSE.
                              ENDIF 
                              ibuffer = ibuffer + 1
                              IF (ppm_recvbufferd(ibuffer) .GT.     &
     &                           (1.0_ppm_kind_double-ppm_myepsd)) THEN
                                 fdata(5,imesh,jmesh,kmesh,isub) = .TRUE.
                              ELSE
                                 fdata(5,imesh,jmesh,kmesh,isub) = .FALSE.
                              ENDIF 
#endif
                           ENDDO
                          ENDDO
                        ENDDO
                     ELSEIF (rtype .EQ. ppm_param_pop_add) THEN
                        DO kmesh=zlo,zhi
                          DO jmesh=ylo,yhi
                           DO imesh=xlo,xhi
#if    __KIND == __SINGLE_PRECISION_COMPLEX | __KIND == __DOUBLE_PRECISION_COMPLEX
                              ibuffer = ibuffer + 2
#else
                              ibuffer = ibuffer + 1
#endif
#if    __KIND == __SINGLE_PRECISION
                              fdata(1,imesh,jmesh,kmesh,isub) =   &
     &                           fdata(1,imesh,jmesh,kmesh,isub)+ &
     &                           REAL(ppm_recvbufferd(ibuffer),ppm_kind_single)
                              ibuffer = ibuffer + 1
                              fdata(2,imesh,jmesh,kmesh,isub) =   &
     &                           fdata(2,imesh,jmesh,kmesh,isub)+ &
     &                           REAL(ppm_recvbufferd(ibuffer),ppm_kind_single)
                              ibuffer = ibuffer + 1
                              fdata(3,imesh,jmesh,kmesh,isub) =   &
     &                           fdata(3,imesh,jmesh,kmesh,isub)+ &
     &                           REAL(ppm_recvbufferd(ibuffer),ppm_kind_single)
                              ibuffer = ibuffer + 1
                              fdata(4,imesh,jmesh,kmesh,isub) =   &
     &                           fdata(4,imesh,jmesh,kmesh,isub)+ &
     &                           REAL(ppm_recvbufferd(ibuffer),ppm_kind_single)
                              ibuffer = ibuffer + 1
                              fdata(5,imesh,jmesh,kmesh,isub) =   &
     &                           fdata(5,imesh,jmesh,kmesh,isub)+ &
     &                           REAL(ppm_recvbufferd(ibuffer),ppm_kind_single)
#elif  __KIND == __DOUBLE_PRECISION
                              fdata(1,imesh,jmesh,kmesh,isub) =   &
     &                           fdata(1,imesh,jmesh,kmesh,isub)+ &
     &                           ppm_recvbufferd(ibuffer)
                              ibuffer = ibuffer + 1
                              fdata(2,imesh,jmesh,kmesh,isub) =   &
     &                           fdata(2,imesh,jmesh,kmesh,isub)+ &
     &                           ppm_recvbufferd(ibuffer)
                              ibuffer = ibuffer + 1
                              fdata(3,imesh,jmesh,kmesh,isub) =   &
     &                           fdata(3,imesh,jmesh,kmesh,isub)+ &
     &                           ppm_recvbufferd(ibuffer)
                              ibuffer = ibuffer + 1
                              fdata(4,imesh,jmesh,kmesh,isub) =   &
     &                           fdata(4,imesh,jmesh,kmesh,isub)+ &
     &                           ppm_recvbufferd(ibuffer)
                              ibuffer = ibuffer + 1
                              fdata(5,imesh,jmesh,kmesh,isub) =   &
     &                           fdata(5,imesh,jmesh,kmesh,isub)+ &
     &                           ppm_recvbufferd(ibuffer)
#elif  __KIND == __SINGLE_PRECISION_COMPLEX
                              fdata(1,imesh,jmesh,kmesh,isub) =   &
     &                           fdata(1,imesh,jmesh,kmesh,isub)+ &
     &                           CMPLX(ppm_recvbufferd(ibuffer-1),  &
     &                           ppm_recvbufferd(ibuffer),ppm_kind_single)
                              ibuffer = ibuffer + 2
                              fdata(2,imesh,jmesh,kmesh,isub) =   &
     &                           fdata(2,imesh,jmesh,kmesh,isub)+ &
     &                           CMPLX(ppm_recvbufferd(ibuffer-1),  &
     &                           ppm_recvbufferd(ibuffer),ppm_kind_single)
                              ibuffer = ibuffer + 2
                              fdata(3,imesh,jmesh,kmesh,isub) =   &
     &                           fdata(3,imesh,jmesh,kmesh,isub)+ &
     &                           CMPLX(ppm_recvbufferd(ibuffer-1),  &
     &                           ppm_recvbufferd(ibuffer),ppm_kind_single)
                              ibuffer = ibuffer + 2
                              fdata(4,imesh,jmesh,kmesh,isub) =   &
     &                           fdata(4,imesh,jmesh,kmesh,isub)+ &
     &                           CMPLX(ppm_recvbufferd(ibuffer-1),  &
     &                           ppm_recvbufferd(ibuffer),ppm_kind_single)
                              ibuffer = ibuffer + 2
                              fdata(5,imesh,jmesh,kmesh,isub) =   &
     &                           fdata(5,imesh,jmesh,kmesh,isub)+ &
     &                           CMPLX(ppm_recvbufferd(ibuffer-1),  &
     &                           ppm_recvbufferd(ibuffer),ppm_kind_single)
#elif  __KIND == __DOUBLE_PRECISION_COMPLEX
                              fdata(1,imesh,jmesh,kmesh,isub) =   &
     &                           fdata(1,imesh,jmesh,kmesh,isub)+ &
     &                           CMPLX(ppm_recvbufferd(ibuffer-1),  &
     &                           ppm_recvbufferd(ibuffer),ppm_kind_double)
                              ibuffer = ibuffer + 2
                              fdata(2,imesh,jmesh,kmesh,isub) =   &
     &                           fdata(2,imesh,jmesh,kmesh,isub)+ &
     &                           CMPLX(ppm_recvbufferd(ibuffer-1),  &
     &                           ppm_recvbufferd(ibuffer),ppm_kind_double)
                              ibuffer = ibuffer + 2
                              fdata(3,imesh,jmesh,kmesh,isub) =   &
     &                           fdata(3,imesh,jmesh,kmesh,isub)+ &
     &                           CMPLX(ppm_recvbufferd(ibuffer-1),  &
     &                           ppm_recvbufferd(ibuffer),ppm_kind_double)
                              ibuffer = ibuffer + 2
                              fdata(4,imesh,jmesh,kmesh,isub) =   &
     &                           fdata(4,imesh,jmesh,kmesh,isub)+ &
     &                           CMPLX(ppm_recvbufferd(ibuffer-1),  &
     &                           ppm_recvbufferd(ibuffer),ppm_kind_double)
                              ibuffer = ibuffer + 2
                              fdata(5,imesh,jmesh,kmesh,isub) =   &
     &                           fdata(5,imesh,jmesh,kmesh,isub)+ &
     &                           CMPLX(ppm_recvbufferd(ibuffer-1),  &
     &                           ppm_recvbufferd(ibuffer),ppm_kind_double)
#elif  __KIND == __INTEGER
                              fdata(1,imesh,jmesh,kmesh,isub) =   &
     &                           fdata(1,imesh,jmesh,kmesh,isub)+ &
     &                           INT(ppm_recvbufferd(ibuffer))
                              ibuffer = ibuffer + 1
                              fdata(2,imesh,jmesh,kmesh,isub) =   &
     &                           fdata(2,imesh,jmesh,kmesh,isub)+ &
     &                           INT(ppm_recvbufferd(ibuffer))
                              ibuffer = ibuffer + 1
                              fdata(3,imesh,jmesh,kmesh,isub) =   &
     &                           fdata(3,imesh,jmesh,kmesh,isub)+ &
     &                           INT(ppm_recvbufferd(ibuffer))
                              ibuffer = ibuffer + 1
                              fdata(4,imesh,jmesh,kmesh,isub) =   &
     &                           fdata(4,imesh,jmesh,kmesh,isub)+ &
     &                           INT(ppm_recvbufferd(ibuffer))
                              ibuffer = ibuffer + 1
                              fdata(5,imesh,jmesh,kmesh,isub) =   &
     &                           fdata(5,imesh,jmesh,kmesh,isub)+ &
     &                           INT(ppm_recvbufferd(ibuffer))
#elif  __KIND == __LOGICAL
                              IF (ppm_recvbufferd(ibuffer) .GT.     &
     &                           (1.0_ppm_kind_double-ppm_myepsd)) THEN
                                 fdata(1,imesh,jmesh,kmesh,isub) =    &
     &                            (fdata(1,imesh,jmesh,kmesh,isub) .AND..TRUE.)
                              ELSE
                                 fdata(1,imesh,jmesh,kmesh,isub) =    &
     &                            (fdata(1,imesh,jmesh,kmesh,isub) .AND..FALSE.)
                              ENDIF 
                              ibuffer = ibuffer + 1
                              IF (ppm_recvbufferd(ibuffer) .GT.     &
     &                           (1.0_ppm_kind_double-ppm_myepsd)) THEN
                                 fdata(2,imesh,jmesh,kmesh,isub) =    &
     &                            (fdata(2,imesh,jmesh,kmesh,isub) .AND..TRUE.)
                              ELSE
                                 fdata(2,imesh,jmesh,kmesh,isub) =    &
     &                            (fdata(2,imesh,jmesh,kmesh,isub) .AND..FALSE.)
                              ENDIF 
                              ibuffer = ibuffer + 1
                              IF (ppm_recvbufferd(ibuffer) .GT.     &
     &                           (1.0_ppm_kind_double-ppm_myepsd)) THEN
                                 fdata(3,imesh,jmesh,kmesh,isub) =    &
     &                            (fdata(3,imesh,jmesh,kmesh,isub) .AND..TRUE.)
                              ELSE
                                 fdata(3,imesh,jmesh,kmesh,isub) =    &
     &                            (fdata(3,imesh,jmesh,kmesh,isub) .AND..FALSE.)
                              ENDIF 
                              ibuffer = ibuffer + 1
                              IF (ppm_recvbufferd(ibuffer) .GT.     &
     &                           (1.0_ppm_kind_double-ppm_myepsd)) THEN
                                 fdata(4,imesh,jmesh,kmesh,isub) =    &
     &                            (fdata(4,imesh,jmesh,kmesh,isub) .AND..TRUE.)
                              ELSE
                                 fdata(4,imesh,jmesh,kmesh,isub) =    &
     &                            (fdata(4,imesh,jmesh,kmesh,isub) .AND..FALSE.)
                              ENDIF 
                              ibuffer = ibuffer + 1
                              IF (ppm_recvbufferd(ibuffer) .GT.     &
     &                           (1.0_ppm_kind_double-ppm_myepsd)) THEN
                                 fdata(5,imesh,jmesh,kmesh,isub) =    &
     &                            (fdata(5,imesh,jmesh,kmesh,isub) .AND..TRUE.)
                              ELSE
                                 fdata(5,imesh,jmesh,kmesh,isub) =    &
     &                            (fdata(5,imesh,jmesh,kmesh,isub) .AND..FALSE.)
                              ENDIF 
#endif
                           ENDDO
                          ENDDO
                        ENDDO
                     ENDIF
                  !-------------------------------------------------------------
                  !  For edim.GT.5 the vector length will be edim !!
                  !-------------------------------------------------------------
                  ELSE
                     IF (rtype .EQ. ppm_param_pop_replace) THEN
                        DO kmesh=zlo,zhi
                          DO jmesh=ylo,yhi
                            DO imesh=xlo,xhi
                              DO k=1,edim
#if    __KIND == __SINGLE_PRECISION_COMPLEX | __KIND == __DOUBLE_PRECISION_COMPLEX
                                 ibuffer = ibuffer + 2
#else
                                 ibuffer = ibuffer + 1
#endif
#if    __KIND == __SINGLE_PRECISION
                                 fdata(k,imesh,jmesh,kmesh,isub) =     &
     &                            REAL(ppm_recvbufferd(ibuffer),ppm_kind_single)
#elif  __KIND == __DOUBLE_PRECISION
                                 fdata(k,imesh,jmesh,kmesh,isub) =     &
     &                              ppm_recvbufferd(ibuffer)
#elif  __KIND == __SINGLE_PRECISION_COMPLEX
                                 fdata(k,imesh,jmesh,kmesh,isub) =     &
     &                              CMPLX(ppm_recvbufferd(ibuffer-1),  &
     &                              ppm_recvbufferd(ibuffer),ppm_kind_single)
#elif  __KIND == __DOUBLE_PRECISION_COMPLEX
                                 fdata(k,imesh,jmesh,kmesh,isub) =     &
     &                              CMPLX(ppm_recvbufferd(ibuffer-1),  &
     &                              ppm_recvbufferd(ibuffer),ppm_kind_double)
#elif  __KIND == __INTEGER
                                 fdata(k,imesh,jmesh,kmesh,isub) =     &
     &                              INT(ppm_recvbufferd(ibuffer))
#elif  __KIND == __LOGICAL
                                 IF (ppm_recvbufferd(ibuffer) .GT.     &
     &                              (1.0_ppm_kind_double-ppm_myepsd)) THEN
                                    fdata(k,imesh,jmesh,kmesh,isub) = .TRUE.
                                 ELSE
                                    fdata(k,imesh,jmesh,kmesh,isub) = .FALSE.
                                 ENDIF 
#endif
                              ENDDO
                            ENDDO
                          ENDDO
                        ENDDO
                     ELSEIF (rtype .EQ. ppm_param_pop_add) THEN
                        DO kmesh=zlo,zhi
                          DO jmesh=ylo,yhi
                            DO imesh=xlo,xhi
                              DO k=1,edim
#if    __KIND == __SINGLE_PRECISION_COMPLEX | __KIND == __DOUBLE_PRECISION_COMPLEX
                                 ibuffer = ibuffer + 2
#else
                                 ibuffer = ibuffer + 1
#endif
#if    __KIND == __SINGLE_PRECISION
                                 fdata(k,imesh,jmesh,kmesh,isub) =   &
     &                            fdata(k,imesh,jmesh,kmesh,isub)+ &
     &                            REAL(ppm_recvbufferd(ibuffer),ppm_kind_single)
#elif  __KIND == __DOUBLE_PRECISION
                                 fdata(k,imesh,jmesh,kmesh,isub) =   &
     &                              fdata(k,imesh,jmesh,kmesh,isub)+ &
     &                              ppm_recvbufferd(ibuffer)
#elif  __KIND == __SINGLE_PRECISION_COMPLEX
                                 fdata(k,imesh,jmesh,kmesh,isub) =   &
     &                              fdata(k,imesh,jmesh,kmesh,isub)+ &
     &                              CMPLX(ppm_recvbufferd(ibuffer-1),  &
     &                              ppm_recvbufferd(ibuffer),ppm_kind_single)
#elif  __KIND == __DOUBLE_PRECISION_COMPLEX
                                 fdata(k,imesh,jmesh,kmesh,isub) =   &
     &                              fdata(k,imesh,jmesh,kmesh,isub)+ &
     &                              CMPLX(ppm_recvbufferd(ibuffer-1),  &
     &                              ppm_recvbufferd(ibuffer),ppm_kind_double)
#elif  __KIND == __INTEGER
                                 fdata(k,imesh,jmesh,kmesh,isub) =   &
     &                              fdata(k,imesh,jmesh,kmesh,isub)+ &
     &                              INT(ppm_recvbufferd(ibuffer))
#elif  __KIND == __LOGICAL
                                 IF (ppm_recvbufferd(ibuffer) .GT.     &
     &                              (1.0_ppm_kind_double-ppm_myepsd)) THEN
                                    fdata(k,imesh,jmesh,kmesh,isub) =    &
     &                                 (fdata(k,imesh,jmesh,kmesh,isub) .AND..TRUE.)
                                 ELSE
                                    fdata(k,imesh,jmesh,kmesh,isub) =    &
     &                                 (fdata(k,imesh,jmesh,kmesh,isub) .AND..FALSE.)
                                 ENDIF 
#endif
                              ENDDO
                            ENDDO
                          ENDDO
                        ENDDO
                     ENDIF         ! rtype
                  ENDIF            ! lda = ...
               !----------------------------------------------------------------
               !  With mask: This will not vectorize so no need to unroll
               !----------------------------------------------------------------
               ELSE
                  DO kmesh=zlo,zhi
                    DO jmesh=ylo,yhi
                     DO imesh=xlo,xhi
                        IF (mask(imesh,jmesh,kmesh,isub)) THEN
                           DO k=1,edim
#if    __KIND == __SINGLE_PRECISION_COMPLEX | __KIND == __DOUBLE_PRECISION_COMPLEX
                              ibuffer = ibuffer + 2
#else
                              ibuffer = ibuffer + 1
#endif
                              IF (rtype .EQ. ppm_param_pop_replace) THEN
#if    __KIND == __SINGLE_PRECISION
                                 fdata(k,imesh,jmesh,kmesh,isub) =     &
     &                            REAL(ppm_recvbufferd(ibuffer),ppm_kind_single)
#elif  __KIND == __DOUBLE_PRECISION
                                 fdata(k,imesh,jmesh,kmesh,isub) =     &
     &                              ppm_recvbufferd(ibuffer)
#elif  __KIND == __SINGLE_PRECISION_COMPLEX
                                 fdata(k,imesh,jmesh,kmesh,isub) =     &
     &                              CMPLX(ppm_recvbufferd(ibuffer-1),  &
     &                              ppm_recvbufferd(ibuffer),ppm_kind_single)
#elif  __KIND == __DOUBLE_PRECISION_COMPLEX
                                 fdata(k,imesh,jmesh,kmesh,isub) =     &
     &                              CMPLX(ppm_recvbufferd(ibuffer-1),  &
     &                              ppm_recvbufferd(ibuffer),ppm_kind_double)
#elif  __KIND == __INTEGER
                                 fdata(k,imesh,jmesh,kmesh,isub) =     &
     &                              INT(ppm_recvbufferd(ibuffer))
#elif  __KIND == __LOGICAL
                                 IF (ppm_recvbufferd(ibuffer) .GT.     &
     &                              (1.0_ppm_kind_double-ppm_myepsd)) THEN
                                    fdata(k,imesh,jmesh,kmesh,isub) = .TRUE.
                                 ELSE
                                    fdata(k,imesh,jmesh,kmesh,isub) = .FALSE.
                                 ENDIF 
#endif
                              ELSEIF (rtype .EQ. ppm_param_pop_add) THEN
#if    __KIND == __SINGLE_PRECISION
                                 fdata(k,imesh,jmesh,kmesh,isub) =   &
     &                            fdata(k,imesh,jmesh,kmesh,isub)+ &
     &                            REAL(ppm_recvbufferd(ibuffer),ppm_kind_single)
#elif  __KIND == __DOUBLE_PRECISION
                                 fdata(k,imesh,jmesh,kmesh,isub) =   &
     &                              fdata(k,imesh,jmesh,kmesh,isub)+ &
     &                              ppm_recvbufferd(ibuffer)
#elif  __KIND == __SINGLE_PRECISION_COMPLEX
                                 fdata(k,imesh,jmesh,kmesh,isub) =   &
     &                              fdata(k,imesh,jmesh,kmesh,isub)+ &
     &                              CMPLX(ppm_recvbufferd(ibuffer-1),  &
     &                              ppm_recvbufferd(ibuffer),ppm_kind_single)
#elif  __KIND == __DOUBLE_PRECISION_COMPLEX
                                 fdata(k,imesh,jmesh,kmesh,isub) =   &
     &                              fdata(k,imesh,jmesh,kmesh,isub)+ &
     &                              CMPLX(ppm_recvbufferd(ibuffer-1),  &
     &                              ppm_recvbufferd(ibuffer),ppm_kind_double)
#elif  __KIND == __INTEGER
                                 fdata(k,imesh,jmesh,kmesh,isub) =   &
     &                              fdata(k,imesh,jmesh,kmesh,isub)+ &
     &                              INT(ppm_recvbufferd(ibuffer))
#elif  __KIND == __LOGICAL
                                 IF (ppm_recvbufferd(ibuffer) .GT.     &
     &                              (1.0_ppm_kind_double-ppm_myepsd)) THEN
                                    fdata(k,imesh,jmesh,kmesh,isub) =    &
     &                                 (fdata(k,imesh,jmesh,kmesh,isub).AND. &
     &                                 .TRUE.)
                                 ELSE
                                    fdata(k,imesh,jmesh,kmesh,isub) =    &
     &                                 (fdata(k,imesh,jmesh,kmesh,isub).AND. &
     &                                 .FALSE.)
                                 ENDIF 
#endif
                              ENDIF
                           ENDDO
                        ELSE
                           !----------------------------------------------------
                           !  If ldo.EQ..FALSE. skip this point
                           !----------------------------------------------------
#if    __KIND == __SINGLE_PRECISION_COMPLEX | __KIND == __DOUBLE_PRECISION_COMPLEX
                           ibuffer = ibuffer + (edim*2)
#else
                           ibuffer = ibuffer + edim
#endif
                        ENDIF   ! mask.EQ..TRUE.
                     ENDDO
                    ENDDO
                  ENDDO
               ENDIF            ! PRESENT(mask)
#elif __DIM == __SFIELD
               IF (.NOT.PRESENT(mask)) THEN
                  IF (rtype .EQ. ppm_param_pop_replace) THEN
                     DO kmesh=zlo,zhi
                      DO jmesh=ylo,yhi
                        DO imesh=xlo,xhi
#if    __KIND == __SINGLE_PRECISION_COMPLEX | __KIND == __DOUBLE_PRECISION_COMPLEX
                           ibuffer = ibuffer + 2
#else
                           ibuffer = ibuffer + 1
#endif
#if    __KIND == __SINGLE_PRECISION
                           fdata(imesh,jmesh,kmesh,isub) =     &
     &                        REAL(ppm_recvbufferd(ibuffer),ppm_kind_single)
#elif  __KIND == __DOUBLE_PRECISION
                           fdata(imesh,jmesh,kmesh,isub) =     &
     &                         ppm_recvbufferd(ibuffer)
#elif  __KIND == __SINGLE_PRECISION_COMPLEX
                           fdata(imesh,jmesh,kmesh,isub) =     &
     &                        CMPLX(ppm_recvbufferd(ibuffer-1),  &
     &                        ppm_recvbufferd(ibuffer),ppm_kind_single)
#elif  __KIND == __DOUBLE_PRECISION_COMPLEX
                           fdata(imesh,jmesh,kmesh,isub) =     &
     &                        CMPLX(ppm_recvbufferd(ibuffer-1),  &
     &                        ppm_recvbufferd(ibuffer),ppm_kind_double)
#elif  __KIND == __INTEGER
                           fdata(imesh,jmesh,kmesh,isub) =     &
     &                        INT(ppm_recvbufferd(ibuffer))
#elif  __KIND == __LOGICAL
                           IF (ppm_recvbufferd(ibuffer) .GT.     &
     &                        (1.0_ppm_kind_double-ppm_myepsd)) THEN
                              fdata(imesh,jmesh,kmesh,isub) = .TRUE.
                           ELSE
                              fdata(imesh,jmesh,kmesh,isub) = .FALSE.
                           ENDIF 
#endif
                        ENDDO
                       ENDDO
                     ENDDO
                  ELSEIF (rtype .EQ. ppm_param_pop_add) THEN
                     DO kmesh=zlo,zhi
                       DO jmesh=ylo,yhi
                        DO imesh=xlo,xhi
#if    __KIND == __SINGLE_PRECISION_COMPLEX | __KIND == __DOUBLE_PRECISION_COMPLEX
                           ibuffer = ibuffer + 2
#else
                           ibuffer = ibuffer + 1
#endif
#if    __KIND == __SINGLE_PRECISION
                           fdata(imesh,jmesh,kmesh,isub) =     &
     &                        fdata(imesh,jmesh,kmesh,isub)+   &
     &                        REAL(ppm_recvbufferd(ibuffer),ppm_kind_single)
#elif  __KIND == __DOUBLE_PRECISION
                           fdata(imesh,jmesh,kmesh,isub) =     &   
     &                         fdata(imesh,jmesh,kmesh,isub)+  &
     &                         ppm_recvbufferd(ibuffer)
#elif  __KIND == __SINGLE_PRECISION_COMPLEX
                           fdata(imesh,jmesh,kmesh,isub) =       &
     &                        fdata(imesh,jmesh,kmesh,isub)+     &
     &                        CMPLX(ppm_recvbufferd(ibuffer-1),  &
     &                        ppm_recvbufferd(ibuffer),ppm_kind_single)
#elif  __KIND == __DOUBLE_PRECISION_COMPLEX
                           fdata(imesh,jmesh,kmesh,isub) =       &
     &                        fdata(imesh,jmesh,kmesh,isub)+     &
     &                        CMPLX(ppm_recvbufferd(ibuffer-1),  &
     &                        ppm_recvbufferd(ibuffer),ppm_kind_double)
#elif  __KIND == __INTEGER
                           fdata(imesh,jmesh,kmesh,isub) =       &
     &                        fdata(imesh,jmesh,kmesh,isub)+     &
     &                        INT(ppm_recvbufferd(ibuffer))
#elif  __KIND == __LOGICAL
                           IF (ppm_recvbufferd(ibuffer) .GT.     &
     &                        (1.0_ppm_kind_double-ppm_myepsd)) THEN
                              fdata(imesh,jmesh,kmesh,isub) =    &
     &                           (fdata(imesh,jmesh,kmesh,isub) .AND..TRUE.)
                           ELSE
                              fdata(imesh,jmesh,kmesh,isub) =    &
     &                           (fdata(imesh,jmesh,kmesh,isub) .AND..FALSE.)
                           ENDIF 
#endif
                        ENDDO
                       ENDDO
                     ENDDO
                  ENDIF
               ELSE             ! PRESENT(mask)
                  DO kmesh=zlo,zhi
                   DO jmesh=ylo,yhi
                     DO imesh=xlo,xhi
                        IF (mask(imesh,jmesh,kmesh,isub)) THEN
#if    __KIND == __SINGLE_PRECISION_COMPLEX | __KIND == __DOUBLE_PRECISION_COMPLEX
                           ibuffer = ibuffer + 2
#else
                           ibuffer = ibuffer + 1
#endif
                           IF (rtype .EQ. ppm_param_pop_replace) THEN
#if    __KIND == __SINGLE_PRECISION
                              fdata(imesh,jmesh,kmesh,isub) =     &
     &                           REAL(ppm_recvbufferd(ibuffer),ppm_kind_single)
#elif  __KIND == __DOUBLE_PRECISION
                              fdata(imesh,jmesh,kmesh,isub) =     & 
     &                            ppm_recvbufferd(ibuffer)
#elif  __KIND == __SINGLE_PRECISION_COMPLEX
                              fdata(imesh,jmesh,kmesh,isub) =     &
     &                           CMPLX(ppm_recvbufferd(ibuffer-1),  &
     &                           ppm_recvbufferd(ibuffer),ppm_kind_single)
#elif  __KIND == __DOUBLE_PRECISION_COMPLEX
                              fdata(imesh,jmesh,kmesh,isub) =     &
     &                           CMPLX(ppm_recvbufferd(ibuffer-1),  &
     &                           ppm_recvbufferd(ibuffer),ppm_kind_double)
#elif  __KIND == __INTEGER
                              fdata(imesh,jmesh,kmesh,isub) =     &
     &                           INT(ppm_recvbufferd(ibuffer))
#elif  __KIND == __LOGICAL
                              IF (ppm_recvbufferd(ibuffer) .GT.     &
     &                           (1.0_ppm_kind_double-ppm_myepsd)) THEN
                                 fdata(imesh,jmesh,kmesh,isub) = .TRUE.
                              ELSE
                                 fdata(imesh,jmesh,kmesh,isub) = .FALSE.
                              ENDIF 
#endif
                           ELSEIF (rtype .EQ. ppm_param_pop_add) THEN
#if    __KIND == __SINGLE_PRECISION
                              fdata(imesh,jmesh,kmesh,isub) =     &
     &                           fdata(imesh,jmesh,kmesh,isub)+   &
     &                           REAL(ppm_recvbufferd(ibuffer),ppm_kind_single)
#elif  __KIND == __DOUBLE_PRECISION
                              fdata(imesh,jmesh,kmesh,isub) =     &
     &                            fdata(imesh,jmesh,kmesh,isub)+  &
     &                            ppm_recvbufferd(ibuffer)
#elif  __KIND == __SINGLE_PRECISION_COMPLEX
                              fdata(imesh,jmesh,kmesh,isub) =     &
     &                           fdata(imesh,jmesh,kmesh,isub)+   &
     &                           CMPLX(ppm_recvbufferd(ibuffer-1),&
     &                           ppm_recvbufferd(ibuffer),ppm_kind_single)
#elif  __KIND == __DOUBLE_PRECISION_COMPLEX
                              fdata(imesh,jmesh,kmesh,isub) =     &
     &                           fdata(imesh,jmesh,kmesh,isub)+   &
     &                           CMPLX(ppm_recvbufferd(ibuffer-1),&
     &                           ppm_recvbufferd(ibuffer),ppm_kind_double)
#elif  __KIND == __INTEGER
                              fdata(imesh,jmesh,kmesh,isub) =     &
     &                           fdata(imesh,jmesh,kmesh,isub)+   &
     &                           INT(ppm_recvbufferd(ibuffer))
#elif  __KIND == __LOGICAL
                              IF (ppm_recvbufferd(ibuffer) .GT.     &
     &                           (1.0_ppm_kind_double-ppm_myepsd)) THEN
                                 fdata(imesh,jmesh,kmesh,isub) =    &
     &                              (fdata(imesh,jmesh,kmesh,isub) .AND..TRUE.)
                              ELSE
                                 fdata(imesh,jmesh,kmesh,isub) =    &
     &                              (fdata(imesh,jmesh,kmesh,isub) .AND..FALSE.)
                              ENDIF 
#endif
                           ENDIF
                        ELSE
                           !----------------------------------------------------
                           !  If ldo.EQ..FALSE. skip this point
                           !----------------------------------------------------
#if    __KIND == __SINGLE_PRECISION_COMPLEX | __KIND == __DOUBLE_PRECISION_COMPLEX
                           ibuffer = ibuffer + 2
#else
                           ibuffer = ibuffer + 1
#endif
                        ENDIF    ! mask.EQ..TRUE.
                     ENDDO
                    ENDDO
                  ENDDO
               ENDIF          ! PRESENT(mask)
#endif
            ENDDO             ! ppm_precvbuffer
         ENDDO                ! ppm_nrecvlist
      !-------------------------------------------------------------------------
      !  SINGLE PRECISION
      !-------------------------------------------------------------------------
      ELSE
         DO i=1,ppm_nrecvlist
            !-------------------------------------------------------------------
            !  Access mesh blocks belonging to the i-th processor in the
            !  recvlist
            !-------------------------------------------------------------------
            DO j=ppm_precvbuffer(i),ppm_precvbuffer(i+1)-1
               !----------------------------------------------------------------
               !  Get the sub ID for this mesh block
               !----------------------------------------------------------------
               jsub = ppm_mesh_irecvtosub(j)
               !----------------------------------------------------------------
               !  Translate to local sub ID for storing the data
               !----------------------------------------------------------------
               isub = invsublist(jsub)
               !----------------------------------------------------------------
               !  Mesh offset for this sub
               !----------------------------------------------------------------
               mofs(1) = ppm_cart_mesh(tomesh,totopo)%istart(1,jsub)-1
               mofs(2) = ppm_cart_mesh(tomesh,totopo)%istart(2,jsub)-1
               mofs(3) = ppm_cart_mesh(tomesh,totopo)%istart(3,jsub)-1
               !----------------------------------------------------------------
               !  Get boundaries of mesh block to be received in local sub
               !  coordinates
               !----------------------------------------------------------------
               xlo = ppm_mesh_irecvblkstart(1,j)-mofs(1)
               ylo = ppm_mesh_irecvblkstart(2,j)-mofs(2)
               zlo = ppm_mesh_irecvblkstart(3,j)-mofs(3)
               xhi = xlo+ppm_mesh_irecvblksize(1,j)-1
               yhi = ylo+ppm_mesh_irecvblksize(2,j)-1
               zhi = zlo+ppm_mesh_irecvblksize(3,j)-1
               IF (ppm_debug .GT. 1) THEN
                   WRITE(mesg,'(A,3I4)') 'start: ',             &
     &                 ppm_mesh_irecvblkstart(1,j),             &
     &                 ppm_mesh_irecvblkstart(2,j),ppm_mesh_irecvblkstart(3,j)
                   CALL ppm_write(ppm_rank,'ppm_map_field_pop_3d',mesg,info)
                   WRITE(mesg,'(A,3I4)') 'size: ',             &
     &                 ppm_mesh_irecvblksize(1,j),             &
     &                 ppm_mesh_irecvblksize(2,j),ppm_mesh_irecvblksize(3,j)
                   CALL ppm_write(ppm_rank,'ppm_map_field_pop_3d',mesg,info)
                   WRITE(mesg,'(A,3I4)') 'mesh offset: ',mofs(1),mofs(2),mofs(3)
                   CALL ppm_write(ppm_rank,'ppm_map_field_pop_3d',mesg,info)
                   WRITE(mesg,'(A,2I4)') 'xlo, xhi: ',xlo,xhi
                   CALL ppm_write(ppm_rank,'ppm_map_field_pop_3d',mesg,info)
                   WRITE(mesg,'(A,2I4)') 'ylo, yhi: ',ylo,yhi
                   CALL ppm_write(ppm_rank,'ppm_map_field_pop_3d',mesg,info)
                   WRITE(mesg,'(A,2I4)') 'zlo, zhi: ',zlo,zhi
                   CALL ppm_write(ppm_rank,'ppm_map_field_pop_3d',mesg,info)
                   WRITE(mesg,'(A,I1)') 'buffer dim: ',edim
                   CALL ppm_write(ppm_rank,'ppm_map_field_pop_3d',mesg,info)
#if   __DIM == __VFIELD
                   WRITE(mesg,'(A,3I4)') 'SIZE(fdata): ',SIZE(fdata,2),  &
     &                 SIZE(fdata,3),SIZE(fdata,4)
                   CALL ppm_write(ppm_rank,'ppm_map_field_pop_3d',mesg,info)
#elif __DIM == __SFIELD
                   WRITE(mesg,'(A,3I4)') 'SIZE(fdata): ',SIZE(fdata,1),  &
     &                 SIZE(fdata,2),SIZE(fdata,3)
                   CALL ppm_write(ppm_rank,'ppm_map_field_pop_3d',mesg,info)
#endif
               ENDIF
#if   __DIM == __VFIELD
               !----------------------------------------------------------------
               !  Without mask: This will vectorize
               !----------------------------------------------------------------
               IF (.NOT.PRESENT(mask)) THEN
                  !-------------------------------------------------------------
                  !  Unrolled for edim=1
                  !-------------------------------------------------------------
                  IF (edim .EQ. 1) THEN
                     IF (rtype .EQ. ppm_param_pop_replace) THEN
                        DO kmesh=zlo,zhi
                          DO jmesh=ylo,yhi
                           DO imesh=xlo,xhi
#if    __KIND == __SINGLE_PRECISION_COMPLEX | __KIND == __DOUBLE_PRECISION_COMPLEX
                              ibuffer = ibuffer + 2
#else
                              ibuffer = ibuffer + 1
#endif
#if    __KIND == __DOUBLE_PRECISION
                              fdata(1,imesh,jmesh,kmesh,isub) =     &
     &                           REAL(ppm_recvbuffers(ibuffer),ppm_kind_double)
#elif  __KIND == __SINGLE_PRECISION
                              fdata(1,imesh,jmesh,kmesh,isub) =     &
     &                           ppm_recvbuffers(ibuffer)
#elif  __KIND == __SINGLE_PRECISION_COMPLEX
                              fdata(1,imesh,jmesh,kmesh,isub) =     &
     &                           CMPLX(ppm_recvbuffers(ibuffer-1),  &
     &                           ppm_recvbuffers(ibuffer),ppm_kind_single)
#elif  __KIND == __DOUBLE_PRECISION_COMPLEX
                              fdata(1,imesh,jmesh,kmesh,isub) =     &
     &                           CMPLX(ppm_recvbuffers(ibuffer-1),  &
     &                           ppm_recvbuffers(ibuffer),ppm_kind_double)
#elif  __KIND == __INTEGER
                              fdata(1,imesh,jmesh,kmesh,isub) =     &
     &                           INT(ppm_recvbuffers(ibuffer))
#elif  __KIND == __LOGICAL
                              IF (ppm_recvbuffers(ibuffer) .GT.     &
     &                           (1.0_ppm_kind_single-ppm_myepss)) THEN
                                 fdata(1,imesh,jmesh,kmesh,isub) = .TRUE.
                              ELSE
                                 fdata(1,imesh,jmesh,kmesh,isub) = .FALSE.
                              ENDIF 
#endif
                           ENDDO
                          ENDDO
                        ENDDO
                     ELSEIF (rtype .EQ. ppm_param_pop_add) THEN
                        DO kmesh=zlo,zhi
                          DO jmesh=ylo,yhi
                           DO imesh=xlo,xhi
#if    __KIND == __SINGLE_PRECISION_COMPLEX | __KIND == __DOUBLE_PRECISION_COMPLEX
                              ibuffer = ibuffer + 2
#else
                              ibuffer = ibuffer + 1
#endif
#if    __KIND == __DOUBLE_PRECISION
                              fdata(1,imesh,jmesh,kmesh,isub) =    &
     &                           fdata(1,imesh,jmesh,kmesh,isub) + &
     &                           REAL(ppm_recvbuffers(ibuffer),ppm_kind_double)
#elif  __KIND == __SINGLE_PRECISION
                              fdata(1,imesh,jmesh,kmesh,isub) =    &
     &                           fdata(1,imesh,jmesh,kmesh,isub) + &
     &                           ppm_recvbuffers(ibuffer)
#elif  __KIND == __SINGLE_PRECISION_COMPLEX
                              fdata(1,imesh,jmesh,kmesh,isub) =    &
     &                           fdata(1,imesh,jmesh,kmesh,isub) + &
     &                           CMPLX(ppm_recvbuffers(ibuffer-1),  &
     &                           ppm_recvbuffers(ibuffer),ppm_kind_single)
#elif  __KIND == __DOUBLE_PRECISION_COMPLEX
                              fdata(1,imesh,jmesh,kmesh,isub) =    &
     &                           fdata(1,imesh,jmesh,kmesh,isub) + &
     &                           CMPLX(ppm_recvbuffers(ibuffer-1),  &
     &                           ppm_recvbuffers(ibuffer),ppm_kind_double)
#elif  __KIND == __INTEGER
                              fdata(1,imesh,jmesh,kmesh,isub) =    &
     &                           fdata(1,imesh,jmesh,kmesh,isub) + &
     &                           INT(ppm_recvbuffers(ibuffer))
#elif  __KIND == __LOGICAL
                              IF (ppm_recvbuffers(ibuffer) .GT.     &
     &                           (1.0_ppm_kind_single-ppm_myepss)) THEN
                                 fdata(1,imesh,jmesh,kmesh,isub) =        &
     &                             (fdata(1,imesh,jmesh,kmesh,isub) .AND..TRUE.)
                              ELSE
                                 fdata(1,imesh,jmesh,kmesh,isub) =        &
     &                            (fdata(1,imesh,jmesh,kmesh,isub) .AND..FALSE.)
                              ENDIF 
#endif
                           ENDDO
                          ENDDO
                        ENDDO
                     ENDIF
                  !-------------------------------------------------------------
                  !  Unrolled for edim=2
                  !-------------------------------------------------------------
                  ELSEIF (edim .EQ. 2) THEN
                     IF (rtype .EQ. ppm_param_pop_replace) THEN
                        DO kmesh=zlo,zhi
                          DO jmesh=ylo,yhi
                           DO imesh=xlo,xhi
#if    __KIND == __SINGLE_PRECISION_COMPLEX | __KIND == __DOUBLE_PRECISION_COMPLEX
                              ibuffer = ibuffer + 2
#else
                              ibuffer = ibuffer + 1
#endif
#if    __KIND == __DOUBLE_PRECISION
                              fdata(1,imesh,jmesh,kmesh,isub) =     &
     &                           REAL(ppm_recvbuffers(ibuffer),ppm_kind_double)
                              ibuffer = ibuffer + 1
                              fdata(2,imesh,jmesh,kmesh,isub) =     &
     &                           REAL(ppm_recvbuffers(ibuffer),ppm_kind_double)
#elif  __KIND == __SINGLE_PRECISION
                              fdata(1,imesh,jmesh,kmesh,isub) =     &
     &                           ppm_recvbuffers(ibuffer)
                              ibuffer = ibuffer + 1
                              fdata(2,imesh,jmesh,kmesh,isub) =     &
     &                           ppm_recvbuffers(ibuffer)
#elif  __KIND == __SINGLE_PRECISION_COMPLEX
                              fdata(1,imesh,jmesh,kmesh,isub) =     &
     &                           CMPLX(ppm_recvbuffers(ibuffer-1),  &
     &                           ppm_recvbuffers(ibuffer),ppm_kind_single)
                              ibuffer = ibuffer + 2
                              fdata(2,imesh,jmesh,kmesh,isub) =     &
     &                           CMPLX(ppm_recvbuffers(ibuffer-1),  &
     &                           ppm_recvbuffers(ibuffer),ppm_kind_single)
#elif  __KIND == __DOUBLE_PRECISION_COMPLEX
                              fdata(1,imesh,jmesh,kmesh,isub) =     &
     &                           CMPLX(ppm_recvbuffers(ibuffer-1),  &
     &                           ppm_recvbuffers(ibuffer),ppm_kind_double)
                              ibuffer = ibuffer + 2
                              fdata(2,imesh,jmesh,kmesh,isub) =     &
     &                           CMPLX(ppm_recvbuffers(ibuffer-1),  &
     &                           ppm_recvbuffers(ibuffer),ppm_kind_double)
#elif  __KIND == __INTEGER
                              fdata(1,imesh,jmesh,kmesh,isub) =     &
     &                           INT(ppm_recvbuffers(ibuffer))
                              ibuffer = ibuffer + 1
                              fdata(2,imesh,jmesh,kmesh,isub) =     &
     &                           INT(ppm_recvbuffers(ibuffer))
#elif  __KIND == __LOGICAL
                              IF (ppm_recvbuffers(ibuffer) .GT.     &
     &                           (1.0_ppm_kind_single-ppm_myepss)) THEN
                                 fdata(1,imesh,jmesh,kmesh,isub) = .TRUE.
                              ELSE
                                 fdata(1,imesh,jmesh,kmesh,isub) = .FALSE.
                              ENDIF 
                              ibuffer = ibuffer + 1
                              IF (ppm_recvbuffers(ibuffer) .GT.     &
     &                           (1.0_ppm_kind_single-ppm_myepss)) THEN
                                 fdata(2,imesh,jmesh,kmesh,isub) = .TRUE.
                              ELSE
                                 fdata(2,imesh,jmesh,kmesh,isub) = .FALSE.
                              ENDIF 
#endif
                           ENDDO
                          ENDDO
                        ENDDO
                     ELSEIF (rtype .EQ. ppm_param_pop_add) THEN
                        DO kmesh=zlo,zhi
                          DO jmesh=ylo,yhi
                           DO imesh=xlo,xhi
#if    __KIND == __SINGLE_PRECISION_COMPLEX | __KIND == __DOUBLE_PRECISION_COMPLEX
                              ibuffer = ibuffer + 2
#else
                              ibuffer = ibuffer + 1
#endif
#if    __KIND == __DOUBLE_PRECISION
                              fdata(1,imesh,jmesh,kmesh,isub) =    &
     &                           fdata(1,imesh,jmesh,kmesh,isub) + &
     &                           REAL(ppm_recvbuffers(ibuffer),ppm_kind_double)
                              ibuffer = ibuffer + 1
                              fdata(2,imesh,jmesh,kmesh,isub) =    &
     &                           fdata(2,imesh,jmesh,kmesh,isub) + &
     &                           REAL(ppm_recvbuffers(ibuffer),ppm_kind_double)
#elif  __KIND == __SINGLE_PRECISION
                              fdata(1,imesh,jmesh,kmesh,isub) =    &
     &                           fdata(1,imesh,jmesh,kmesh,isub) + &
     &                           ppm_recvbuffers(ibuffer)
                              ibuffer = ibuffer + 1
                              fdata(2,imesh,jmesh,kmesh,isub) =    &
     &                           fdata(2,imesh,jmesh,kmesh,isub) + &
     &                           ppm_recvbuffers(ibuffer)
#elif  __KIND == __SINGLE_PRECISION_COMPLEX
                              fdata(1,imesh,jmesh,kmesh,isub) =    &
     &                           fdata(1,imesh,jmesh,kmesh,isub) + &
     &                           CMPLX(ppm_recvbuffers(ibuffer-1),  &
     &                           ppm_recvbuffers(ibuffer),ppm_kind_single)
                              ibuffer = ibuffer + 2
                              fdata(2,imesh,jmesh,kmesh,isub) =    &
     &                           fdata(2,imesh,jmesh,kmesh,isub) + &
     &                           CMPLX(ppm_recvbuffers(ibuffer-1),  &
     &                           ppm_recvbuffers(ibuffer),ppm_kind_single)
#elif  __KIND == __DOUBLE_PRECISION_COMPLEX
                              fdata(1,imesh,jmesh,kmesh,isub) =    &
     &                           fdata(1,imesh,jmesh,kmesh,isub) + &
     &                           CMPLX(ppm_recvbuffers(ibuffer-1),  &
     &                           ppm_recvbuffers(ibuffer),ppm_kind_double)
                              ibuffer = ibuffer + 2
                              fdata(2,imesh,jmesh,kmesh,isub) =    &
     &                           fdata(2,imesh,jmesh,kmesh,isub) + &
     &                           CMPLX(ppm_recvbuffers(ibuffer-1),  &
     &                           ppm_recvbuffers(ibuffer),ppm_kind_double)
#elif  __KIND == __INTEGER
                              fdata(1,imesh,jmesh,kmesh,isub) =    &
     &                           fdata(1,imesh,jmesh,kmesh,isub) + &
     &                           INT(ppm_recvbuffers(ibuffer))
                              ibuffer = ibuffer + 1
                              fdata(2,imesh,jmesh,kmesh,isub) =    &
     &                           fdata(2,imesh,jmesh,kmesh,isub) + &
     &                           INT(ppm_recvbuffers(ibuffer))
#elif  __KIND == __LOGICAL
                              IF (ppm_recvbuffers(ibuffer) .GT.     &
     &                           (1.0_ppm_kind_single-ppm_myepss)) THEN
                                 fdata(1,imesh,jmesh,kmesh,isub) =        &
     &                             (fdata(1,imesh,jmesh,kmesh,isub) .AND..TRUE.)
                              ELSE
                                 fdata(1,imesh,jmesh,kmesh,isub) =        &
     &                            (fdata(1,imesh,jmesh,kmesh,isub) .AND..FALSE.)
                              ENDIF 
                              ibuffer = ibuffer + 1
                              IF (ppm_recvbuffers(ibuffer) .GT.     &
     &                           (1.0_ppm_kind_single-ppm_myepss)) THEN
                                 fdata(2,imesh,jmesh,kmesh,isub) =        &
     &                             (fdata(2,imesh,jmesh,kmesh,isub) .AND..TRUE.)
                              ELSE
                                 fdata(2,imesh,jmesh,kmesh,isub) =        &
     &                            (fdata(2,imesh,jmesh,kmesh,isub) .AND..FALSE.)
                              ENDIF 
#endif
                           ENDDO
                          ENDDO
                        ENDDO
                     ENDIF
                  !-------------------------------------------------------------
                  !  Unrolled for edim=3
                  !-------------------------------------------------------------
                  ELSEIF (edim .EQ. 3) THEN
                     IF (rtype .EQ. ppm_param_pop_replace) THEN
                        DO kmesh=zlo,zhi
                          DO jmesh=ylo,yhi
                           DO imesh=xlo,xhi
#if    __KIND == __SINGLE_PRECISION_COMPLEX | __KIND == __DOUBLE_PRECISION_COMPLEX
                              ibuffer = ibuffer + 2
#else
                              ibuffer = ibuffer + 1
#endif
#if    __KIND == __DOUBLE_PRECISION
                              fdata(1,imesh,jmesh,kmesh,isub) =     &
     &                           REAL(ppm_recvbuffers(ibuffer),ppm_kind_double)
                              ibuffer = ibuffer + 1
                              fdata(2,imesh,jmesh,kmesh,isub) =     &
     &                           REAL(ppm_recvbuffers(ibuffer),ppm_kind_double)
                              ibuffer = ibuffer + 1
                              fdata(3,imesh,jmesh,kmesh,isub) =     &
     &                           REAL(ppm_recvbuffers(ibuffer),ppm_kind_double)
#elif  __KIND == __SINGLE_PRECISION
                              fdata(1,imesh,jmesh,kmesh,isub) =     &
     &                           ppm_recvbuffers(ibuffer)
                              ibuffer = ibuffer + 1
                              fdata(2,imesh,jmesh,kmesh,isub) =     &
     &                           ppm_recvbuffers(ibuffer)
                              ibuffer = ibuffer + 1
                              fdata(3,imesh,jmesh,kmesh,isub) =     &
     &                           ppm_recvbuffers(ibuffer)
#elif  __KIND == __SINGLE_PRECISION_COMPLEX
                              fdata(1,imesh,jmesh,kmesh,isub) =     &
     &                           CMPLX(ppm_recvbuffers(ibuffer-1),  &
     &                           ppm_recvbuffers(ibuffer),ppm_kind_single)
                              ibuffer = ibuffer + 2
                              fdata(2,imesh,jmesh,kmesh,isub) =     &
     &                           CMPLX(ppm_recvbuffers(ibuffer-1),  &
     &                           ppm_recvbuffers(ibuffer),ppm_kind_single)
                              ibuffer = ibuffer + 2
                              fdata(3,imesh,jmesh,kmesh,isub) =     &
     &                           CMPLX(ppm_recvbuffers(ibuffer-1),  &
     &                           ppm_recvbuffers(ibuffer),ppm_kind_single)
#elif  __KIND == __DOUBLE_PRECISION_COMPLEX
                              fdata(1,imesh,jmesh,kmesh,isub) =     &
     &                           CMPLX(ppm_recvbuffers(ibuffer-1),  &
     &                           ppm_recvbuffers(ibuffer),ppm_kind_double)
                              ibuffer = ibuffer + 2
                              fdata(2,imesh,jmesh,kmesh,isub) =     &
     &                           CMPLX(ppm_recvbuffers(ibuffer-1),  &
     &                           ppm_recvbuffers(ibuffer),ppm_kind_double)
                              ibuffer = ibuffer + 2
                              fdata(3,imesh,jmesh,kmesh,isub) =     &
     &                           CMPLX(ppm_recvbuffers(ibuffer-1),  &
     &                           ppm_recvbuffers(ibuffer),ppm_kind_double)
#elif  __KIND == __INTEGER
                              fdata(1,imesh,jmesh,kmesh,isub) =     &
     &                           INT(ppm_recvbuffers(ibuffer))
                              ibuffer = ibuffer + 1
                              fdata(2,imesh,jmesh,kmesh,isub) =     &
     &                           INT(ppm_recvbuffers(ibuffer))
                              ibuffer = ibuffer + 1
                              fdata(3,imesh,jmesh,kmesh,isub) =     &
     &                           INT(ppm_recvbuffers(ibuffer))
#elif  __KIND == __LOGICAL
                              IF (ppm_recvbuffers(ibuffer) .GT.     &
     &                           (1.0_ppm_kind_single-ppm_myepss)) THEN
                                 fdata(1,imesh,jmesh,kmesh,isub) = .TRUE.
                              ELSE
                                 fdata(1,imesh,jmesh,kmesh,isub) = .FALSE.
                              ENDIF 
                              ibuffer = ibuffer + 1
                              IF (ppm_recvbuffers(ibuffer) .GT.     &
     &                           (1.0_ppm_kind_single-ppm_myepss)) THEN
                                 fdata(2,imesh,jmesh,kmesh,isub) = .TRUE.
                              ELSE
                                 fdata(2,imesh,jmesh,kmesh,isub) = .FALSE.
                              ENDIF 
                              ibuffer = ibuffer + 1
                              IF (ppm_recvbuffers(ibuffer) .GT.     &
     &                           (1.0_ppm_kind_single-ppm_myepss)) THEN
                                 fdata(3,imesh,jmesh,kmesh,isub) = .TRUE.
                              ELSE
                                 fdata(3,imesh,jmesh,kmesh,isub) = .FALSE.
                              ENDIF 
#endif
                           ENDDO
                          ENDDO
                        ENDDO
                     ELSEIF (rtype .EQ. ppm_param_pop_add) THEN
                        DO kmesh=zlo,zhi
                          DO jmesh=ylo,yhi
                           DO imesh=xlo,xhi
#if    __KIND == __SINGLE_PRECISION_COMPLEX | __KIND == __DOUBLE_PRECISION_COMPLEX
                              ibuffer = ibuffer + 2
#else
                              ibuffer = ibuffer + 1
#endif
#if    __KIND == __DOUBLE_PRECISION
                              fdata(1,imesh,jmesh,kmesh,isub) =    &
     &                           fdata(1,imesh,jmesh,kmesh,isub) + &
     &                           REAL(ppm_recvbuffers(ibuffer),ppm_kind_double)
                              ibuffer = ibuffer + 1
                              fdata(2,imesh,jmesh,kmesh,isub) =    &
     &                           fdata(2,imesh,jmesh,kmesh,isub) + &
     &                           REAL(ppm_recvbuffers(ibuffer),ppm_kind_double)
                              ibuffer = ibuffer + 1
                              fdata(3,imesh,jmesh,kmesh,isub) =    &
     &                           fdata(3,imesh,jmesh,kmesh,isub) + &
     &                           REAL(ppm_recvbuffers(ibuffer),ppm_kind_double)
#elif  __KIND == __SINGLE_PRECISION
                              fdata(1,imesh,jmesh,kmesh,isub) =    &
     &                           fdata(1,imesh,jmesh,kmesh,isub) + &
     &                           ppm_recvbuffers(ibuffer)
                              ibuffer = ibuffer + 1
                              fdata(2,imesh,jmesh,kmesh,isub) =    &
     &                           fdata(2,imesh,jmesh,kmesh,isub) + &
     &                           ppm_recvbuffers(ibuffer)
                              ibuffer = ibuffer + 1
                              fdata(3,imesh,jmesh,kmesh,isub) =    &
     &                           fdata(3,imesh,jmesh,kmesh,isub) + &
     &                           ppm_recvbuffers(ibuffer)
#elif  __KIND == __SINGLE_PRECISION_COMPLEX
                              fdata(1,imesh,jmesh,kmesh,isub) =    &
     &                           fdata(1,imesh,jmesh,kmesh,isub) + &
     &                           CMPLX(ppm_recvbuffers(ibuffer-1),  &
     &                           ppm_recvbuffers(ibuffer),ppm_kind_single)
                              ibuffer = ibuffer + 2
                              fdata(2,imesh,jmesh,kmesh,isub) =    &
     &                           fdata(2,imesh,jmesh,kmesh,isub) + &
     &                           CMPLX(ppm_recvbuffers(ibuffer-1),  &
     &                           ppm_recvbuffers(ibuffer),ppm_kind_single)
                              ibuffer = ibuffer + 2
                              fdata(3,imesh,jmesh,kmesh,isub) =    &
     &                           fdata(3,imesh,jmesh,kmesh,isub) + &
     &                           CMPLX(ppm_recvbuffers(ibuffer-1),  &
     &                           ppm_recvbuffers(ibuffer),ppm_kind_single)
#elif  __KIND == __DOUBLE_PRECISION_COMPLEX
                              fdata(1,imesh,jmesh,kmesh,isub) =    &
     &                           fdata(1,imesh,jmesh,kmesh,isub) + &
     &                           CMPLX(ppm_recvbuffers(ibuffer-1),  &
     &                           ppm_recvbuffers(ibuffer),ppm_kind_double)
                              ibuffer = ibuffer + 2
                              fdata(2,imesh,jmesh,kmesh,isub) =    &
     &                           fdata(2,imesh,jmesh,kmesh,isub) + &
     &                           CMPLX(ppm_recvbuffers(ibuffer-1),  &
     &                           ppm_recvbuffers(ibuffer),ppm_kind_double)
                              ibuffer = ibuffer + 2
                              fdata(3,imesh,jmesh,kmesh,isub) =    &
     &                           fdata(3,imesh,jmesh,kmesh,isub) + &
     &                           CMPLX(ppm_recvbuffers(ibuffer-1),  &
     &                           ppm_recvbuffers(ibuffer),ppm_kind_double)
#elif  __KIND == __INTEGER
                              fdata(1,imesh,jmesh,kmesh,isub) =    &
     &                           fdata(1,imesh,jmesh,kmesh,isub) + &
     &                           INT(ppm_recvbuffers(ibuffer))
                              ibuffer = ibuffer + 1
                              fdata(2,imesh,jmesh,kmesh,isub) =    &
     &                           fdata(2,imesh,jmesh,kmesh,isub) + &
     &                           INT(ppm_recvbuffers(ibuffer))
                              ibuffer = ibuffer + 1
                              fdata(3,imesh,jmesh,kmesh,isub) =    &
     &                           fdata(3,imesh,jmesh,kmesh,isub) + &
     &                           INT(ppm_recvbuffers(ibuffer))
#elif  __KIND == __LOGICAL
                              IF (ppm_recvbuffers(ibuffer) .GT.     &
     &                           (1.0_ppm_kind_single-ppm_myepss)) THEN
                                 fdata(1,imesh,jmesh,kmesh,isub) =        &
     &                             (fdata(1,imesh,jmesh,kmesh,isub) .AND..TRUE.)
                              ELSE
                                 fdata(1,imesh,jmesh,kmesh,isub) =        &
     &                            (fdata(1,imesh,jmesh,kmesh,isub) .AND..FALSE.)
                              ENDIF 
                              ibuffer = ibuffer + 1
                              IF (ppm_recvbuffers(ibuffer) .GT.     &
     &                           (1.0_ppm_kind_single-ppm_myepss)) THEN
                                 fdata(2,imesh,jmesh,kmesh,isub) =        &
     &                             (fdata(2,imesh,jmesh,kmesh,isub) .AND..TRUE.)
                              ELSE
                                 fdata(2,imesh,jmesh,kmesh,isub) =        &
     &                            (fdata(2,imesh,jmesh,kmesh,isub) .AND..FALSE.)
                              ENDIF 
                              ibuffer = ibuffer + 1
                              IF (ppm_recvbuffers(ibuffer) .GT.     &
     &                           (1.0_ppm_kind_single-ppm_myepss)) THEN
                                 fdata(3,imesh,jmesh,kmesh,isub) =        &
     &                             (fdata(3,imesh,jmesh,kmesh,isub) .AND..TRUE.)
                              ELSE
                                 fdata(3,imesh,jmesh,kmesh,isub) =        &
     &                            (fdata(3,imesh,jmesh,kmesh,isub) .AND..FALSE.)
                              ENDIF 
#endif
                           ENDDO
                          ENDDO
                        ENDDO
                     ENDIF
                  !-------------------------------------------------------------
                  !  Unrolled for edim=4
                  !-------------------------------------------------------------
                  ELSEIF (edim .EQ. 4) THEN
                     IF (rtype .EQ. ppm_param_pop_replace) THEN
                        DO kmesh=zlo,zhi
                          DO jmesh=ylo,yhi
                           DO imesh=xlo,xhi
#if    __KIND == __SINGLE_PRECISION_COMPLEX | __KIND == __DOUBLE_PRECISION_COMPLEX
                              ibuffer = ibuffer + 2
#else
                              ibuffer = ibuffer + 1
#endif
#if    __KIND == __DOUBLE_PRECISION
                              fdata(1,imesh,jmesh,kmesh,isub) =     &
     &                           REAL(ppm_recvbuffers(ibuffer),ppm_kind_double)
                              ibuffer = ibuffer + 1
                              fdata(2,imesh,jmesh,kmesh,isub) =     &
     &                           REAL(ppm_recvbuffers(ibuffer),ppm_kind_double)
                              ibuffer = ibuffer + 1
                              fdata(3,imesh,jmesh,kmesh,isub) =     &
     &                           REAL(ppm_recvbuffers(ibuffer),ppm_kind_double)
                              ibuffer = ibuffer + 1
                              fdata(4,imesh,jmesh,kmesh,isub) =     &
     &                           REAL(ppm_recvbuffers(ibuffer),ppm_kind_double)
#elif  __KIND == __SINGLE_PRECISION
                              fdata(1,imesh,jmesh,kmesh,isub) =     &
     &                           ppm_recvbuffers(ibuffer)
                              ibuffer = ibuffer + 1
                              fdata(2,imesh,jmesh,kmesh,isub) =     &
     &                           ppm_recvbuffers(ibuffer)
                              ibuffer = ibuffer + 1
                              fdata(3,imesh,jmesh,kmesh,isub) =     &
     &                           ppm_recvbuffers(ibuffer)
                              ibuffer = ibuffer + 1
                              fdata(4,imesh,jmesh,kmesh,isub) =     &
     &                           ppm_recvbuffers(ibuffer)
#elif  __KIND == __SINGLE_PRECISION_COMPLEX
                              fdata(1,imesh,jmesh,kmesh,isub) =     &
     &                           CMPLX(ppm_recvbuffers(ibuffer-1),  &
     &                           ppm_recvbuffers(ibuffer),ppm_kind_single)
                              ibuffer = ibuffer + 2
                              fdata(2,imesh,jmesh,kmesh,isub) =     &
     &                           CMPLX(ppm_recvbuffers(ibuffer-1),  &
     &                           ppm_recvbuffers(ibuffer),ppm_kind_single)
                              ibuffer = ibuffer + 2
                              fdata(3,imesh,jmesh,kmesh,isub) =     &
     &                           CMPLX(ppm_recvbuffers(ibuffer-1),  &
     &                           ppm_recvbuffers(ibuffer),ppm_kind_single)
                              ibuffer = ibuffer + 2
                              fdata(4,imesh,jmesh,kmesh,isub) =     &
     &                           CMPLX(ppm_recvbuffers(ibuffer-1),  &
     &                           ppm_recvbuffers(ibuffer),ppm_kind_single)
#elif  __KIND == __DOUBLE_PRECISION_COMPLEX
                              fdata(1,imesh,jmesh,kmesh,isub) =     &
     &                           CMPLX(ppm_recvbuffers(ibuffer-1),  &
     &                           ppm_recvbuffers(ibuffer),ppm_kind_double)
                              ibuffer = ibuffer + 2
                              fdata(2,imesh,jmesh,kmesh,isub) =     &
     &                           CMPLX(ppm_recvbuffers(ibuffer-1),  &
     &                           ppm_recvbuffers(ibuffer),ppm_kind_double)
                              ibuffer = ibuffer + 2
                              fdata(3,imesh,jmesh,kmesh,isub) =     &
     &                           CMPLX(ppm_recvbuffers(ibuffer-1),  &
     &                           ppm_recvbuffers(ibuffer),ppm_kind_double)
                              ibuffer = ibuffer + 2
                              fdata(4,imesh,jmesh,kmesh,isub) =     &
     &                           CMPLX(ppm_recvbuffers(ibuffer-1),  &
     &                           ppm_recvbuffers(ibuffer),ppm_kind_double)
#elif  __KIND == __INTEGER
                              fdata(1,imesh,jmesh,kmesh,isub) =     &
     &                           INT(ppm_recvbuffers(ibuffer))
                              ibuffer = ibuffer + 1
                              fdata(2,imesh,jmesh,kmesh,isub) =     &
     &                           INT(ppm_recvbuffers(ibuffer))
                              ibuffer = ibuffer + 1
                              fdata(3,imesh,jmesh,kmesh,isub) =     &
     &                           INT(ppm_recvbuffers(ibuffer))
                              ibuffer = ibuffer + 1
                              fdata(4,imesh,jmesh,kmesh,isub) =     &
     &                           INT(ppm_recvbuffers(ibuffer))
#elif  __KIND == __LOGICAL
                              IF (ppm_recvbuffers(ibuffer) .GT.     &
     &                           (1.0_ppm_kind_single-ppm_myepss)) THEN
                                 fdata(1,imesh,jmesh,kmesh,isub) = .TRUE.
                              ELSE
                                 fdata(1,imesh,jmesh,kmesh,isub) = .FALSE.
                              ENDIF 
                              ibuffer = ibuffer + 1
                              IF (ppm_recvbuffers(ibuffer) .GT.     &
     &                           (1.0_ppm_kind_single-ppm_myepss)) THEN
                                 fdata(2,imesh,jmesh,kmesh,isub) = .TRUE.
                              ELSE
                                 fdata(2,imesh,jmesh,kmesh,isub) = .FALSE.
                              ENDIF 
                              ibuffer = ibuffer + 1
                              IF (ppm_recvbuffers(ibuffer) .GT.     &
     &                           (1.0_ppm_kind_single-ppm_myepss)) THEN
                                 fdata(3,imesh,jmesh,kmesh,isub) = .TRUE.
                              ELSE
                                 fdata(3,imesh,jmesh,kmesh,isub) = .FALSE.
                              ENDIF 
                              ibuffer = ibuffer + 1
                              IF (ppm_recvbuffers(ibuffer) .GT.     &
     &                           (1.0_ppm_kind_single-ppm_myepss)) THEN
                                 fdata(4,imesh,jmesh,kmesh,isub) = .TRUE.
                              ELSE
                                 fdata(4,imesh,jmesh,kmesh,isub) = .FALSE.
                              ENDIF 
#endif
                           ENDDO
                          ENDDO
                        ENDDO
                     ELSEIF (rtype .EQ. ppm_param_pop_add) THEN
                        DO kmesh=zlo,zhi
                          DO jmesh=ylo,yhi
                           DO imesh=xlo,xhi
#if    __KIND == __SINGLE_PRECISION_COMPLEX | __KIND == __DOUBLE_PRECISION_COMPLEX
                              ibuffer = ibuffer + 2
#else
                              ibuffer = ibuffer + 1
#endif
#if    __KIND == __DOUBLE_PRECISION
                              fdata(1,imesh,jmesh,kmesh,isub) =    &
     &                           fdata(1,imesh,jmesh,kmesh,isub) + &
     &                           REAL(ppm_recvbuffers(ibuffer),ppm_kind_double)
                              ibuffer = ibuffer + 1
                              fdata(2,imesh,jmesh,kmesh,isub) =    &
     &                           fdata(2,imesh,jmesh,kmesh,isub) + &
     &                           REAL(ppm_recvbuffers(ibuffer),ppm_kind_double)
                              ibuffer = ibuffer + 1
                              fdata(3,imesh,jmesh,kmesh,isub) =    &
     &                           fdata(3,imesh,jmesh,kmesh,isub) + &
     &                           REAL(ppm_recvbuffers(ibuffer),ppm_kind_double)
                              ibuffer = ibuffer + 1
                              fdata(4,imesh,jmesh,kmesh,isub) =    &
     &                           fdata(4,imesh,jmesh,kmesh,isub) + &
     &                           REAL(ppm_recvbuffers(ibuffer),ppm_kind_double)
#elif  __KIND == __SINGLE_PRECISION
                              fdata(1,imesh,jmesh,kmesh,isub) =    &
     &                           fdata(1,imesh,jmesh,kmesh,isub) + &
     &                           ppm_recvbuffers(ibuffer)
                              ibuffer = ibuffer + 1
                              fdata(2,imesh,jmesh,kmesh,isub) =    &
     &                           fdata(2,imesh,jmesh,kmesh,isub) + &
     &                           ppm_recvbuffers(ibuffer)
                              ibuffer = ibuffer + 1
                              fdata(3,imesh,jmesh,kmesh,isub) =    &
     &                           fdata(3,imesh,jmesh,kmesh,isub) + &
     &                           ppm_recvbuffers(ibuffer)
                              ibuffer = ibuffer + 1
                              fdata(4,imesh,jmesh,kmesh,isub) =    &
     &                           fdata(4,imesh,jmesh,kmesh,isub) + &
     &                           ppm_recvbuffers(ibuffer)
#elif  __KIND == __SINGLE_PRECISION_COMPLEX
                              fdata(1,imesh,jmesh,kmesh,isub) =    &
     &                           fdata(1,imesh,jmesh,kmesh,isub) + &
     &                           CMPLX(ppm_recvbuffers(ibuffer-1),  &
     &                           ppm_recvbuffers(ibuffer),ppm_kind_single)
                              ibuffer = ibuffer + 2
                              fdata(2,imesh,jmesh,kmesh,isub) =    &
     &                           fdata(2,imesh,jmesh,kmesh,isub) + &
     &                           CMPLX(ppm_recvbuffers(ibuffer-1),  &
     &                           ppm_recvbuffers(ibuffer),ppm_kind_single)
                              ibuffer = ibuffer + 2
                              fdata(3,imesh,jmesh,kmesh,isub) =    &
     &                           fdata(3,imesh,jmesh,kmesh,isub) + &
     &                           CMPLX(ppm_recvbuffers(ibuffer-1),  &
     &                           ppm_recvbuffers(ibuffer),ppm_kind_single)
                              ibuffer = ibuffer + 2
                              fdata(4,imesh,jmesh,kmesh,isub) =    &
     &                           fdata(4,imesh,jmesh,kmesh,isub) + &
     &                           CMPLX(ppm_recvbuffers(ibuffer-1),  &
     &                           ppm_recvbuffers(ibuffer),ppm_kind_single)
#elif  __KIND == __DOUBLE_PRECISION_COMPLEX
                              fdata(1,imesh,jmesh,kmesh,isub) =    &
     &                           fdata(1,imesh,jmesh,kmesh,isub) + &
     &                           CMPLX(ppm_recvbuffers(ibuffer-1),  &
     &                           ppm_recvbuffers(ibuffer),ppm_kind_double)
                              ibuffer = ibuffer + 2
                              fdata(2,imesh,jmesh,kmesh,isub) =    &
     &                           fdata(2,imesh,jmesh,kmesh,isub) + &
     &                           CMPLX(ppm_recvbuffers(ibuffer-1),  &
     &                           ppm_recvbuffers(ibuffer),ppm_kind_double)
                              ibuffer = ibuffer + 2
                              fdata(3,imesh,jmesh,kmesh,isub) =    &
     &                           fdata(3,imesh,jmesh,kmesh,isub) + &
     &                           CMPLX(ppm_recvbuffers(ibuffer-1),  &
     &                           ppm_recvbuffers(ibuffer),ppm_kind_double)
                              ibuffer = ibuffer + 2
                              fdata(4,imesh,jmesh,kmesh,isub) =    &
     &                           fdata(4,imesh,jmesh,kmesh,isub) + &
     &                           CMPLX(ppm_recvbuffers(ibuffer-1),  &
     &                           ppm_recvbuffers(ibuffer),ppm_kind_double)
#elif  __KIND == __INTEGER
                              fdata(1,imesh,jmesh,kmesh,isub) =    &
     &                           fdata(1,imesh,jmesh,kmesh,isub) + &
     &                           INT(ppm_recvbuffers(ibuffer))
                              ibuffer = ibuffer + 1
                              fdata(2,imesh,jmesh,kmesh,isub) =    &
     &                           fdata(2,imesh,jmesh,kmesh,isub) + &
     &                           INT(ppm_recvbuffers(ibuffer))
                              ibuffer = ibuffer + 1
                              fdata(3,imesh,jmesh,kmesh,isub) =    &
     &                           fdata(3,imesh,jmesh,kmesh,isub) + &
     &                           INT(ppm_recvbuffers(ibuffer))
                              ibuffer = ibuffer + 1
                              fdata(4,imesh,jmesh,kmesh,isub) =    &
     &                           fdata(4,imesh,jmesh,kmesh,isub) + &
     &                           INT(ppm_recvbuffers(ibuffer))
#elif  __KIND == __LOGICAL
                              IF (ppm_recvbuffers(ibuffer) .GT.     &
     &                           (1.0_ppm_kind_single-ppm_myepss)) THEN
                                 fdata(1,imesh,jmesh,kmesh,isub) =        &
     &                             (fdata(1,imesh,jmesh,kmesh,isub) .AND..TRUE.)
                              ELSE
                                 fdata(1,imesh,jmesh,kmesh,isub) =        &
     &                            (fdata(1,imesh,jmesh,kmesh,isub) .AND..FALSE.)
                              ENDIF 
                              ibuffer = ibuffer + 1
                              IF (ppm_recvbuffers(ibuffer) .GT.     &
     &                           (1.0_ppm_kind_single-ppm_myepss)) THEN
                                 fdata(2,imesh,jmesh,kmesh,isub) =        &
     &                             (fdata(2,imesh,jmesh,kmesh,isub) .AND..TRUE.)
                              ELSE
                                 fdata(2,imesh,jmesh,kmesh,isub) =        &
     &                            (fdata(2,imesh,jmesh,kmesh,isub) .AND..FALSE.)
                              ENDIF 
                              ibuffer = ibuffer + 1
                              IF (ppm_recvbuffers(ibuffer) .GT.     &
     &                           (1.0_ppm_kind_single-ppm_myepss)) THEN
                                 fdata(3,imesh,jmesh,kmesh,isub) =        &
     &                             (fdata(3,imesh,jmesh,kmesh,isub) .AND..TRUE.)
                              ELSE
                                 fdata(3,imesh,jmesh,kmesh,isub) =        &
     &                            (fdata(3,imesh,jmesh,kmesh,isub) .AND..FALSE.)
                              ENDIF 
                              ibuffer = ibuffer + 1
                              IF (ppm_recvbuffers(ibuffer) .GT.     &
     &                           (1.0_ppm_kind_single-ppm_myepss)) THEN
                                 fdata(4,imesh,jmesh,kmesh,isub) =        &
     &                             (fdata(4,imesh,jmesh,kmesh,isub) .AND..TRUE.)
                              ELSE
                                 fdata(4,imesh,jmesh,kmesh,isub) =        &
     &                            (fdata(4,imesh,jmesh,kmesh,isub) .AND..FALSE.)
                              ENDIF 
#endif
                           ENDDO
                          ENDDO
                        ENDDO
                     ENDIF
                  !-------------------------------------------------------------
                  !  Unrolled for edim=5
                  !-------------------------------------------------------------
                  ELSEIF (edim .EQ. 5) THEN
                     IF (rtype .EQ. ppm_param_pop_replace) THEN
                        DO kmesh=zlo,zhi
                          DO jmesh=ylo,yhi
                           DO imesh=xlo,xhi
#if    __KIND == __SINGLE_PRECISION_COMPLEX | __KIND == __DOUBLE_PRECISION_COMPLEX
                              ibuffer = ibuffer + 2
#else
                              ibuffer = ibuffer + 1
#endif
#if    __KIND == __DOUBLE_PRECISION
                              fdata(1,imesh,jmesh,kmesh,isub) =     &
     &                           REAL(ppm_recvbuffers(ibuffer),ppm_kind_double)
                              ibuffer = ibuffer + 1
                              fdata(2,imesh,jmesh,kmesh,isub) =     &
     &                           REAL(ppm_recvbuffers(ibuffer),ppm_kind_double)
                              ibuffer = ibuffer + 1
                              fdata(3,imesh,jmesh,kmesh,isub) =     &
     &                           REAL(ppm_recvbuffers(ibuffer),ppm_kind_double)
                              ibuffer = ibuffer + 1
                              fdata(4,imesh,jmesh,kmesh,isub) =     &
     &                           REAL(ppm_recvbuffers(ibuffer),ppm_kind_double)
                              ibuffer = ibuffer + 1
                              fdata(5,imesh,jmesh,kmesh,isub) =     &
     &                           REAL(ppm_recvbuffers(ibuffer),ppm_kind_double)
#elif  __KIND == __SINGLE_PRECISION
                              fdata(1,imesh,jmesh,kmesh,isub) =     &
     &                           ppm_recvbuffers(ibuffer)
                              ibuffer = ibuffer + 1
                              fdata(2,imesh,jmesh,kmesh,isub) =     &
     &                           ppm_recvbuffers(ibuffer)
                              ibuffer = ibuffer + 1
                              fdata(3,imesh,jmesh,kmesh,isub) =     &
     &                           ppm_recvbuffers(ibuffer)
                              ibuffer = ibuffer + 1
                              fdata(4,imesh,jmesh,kmesh,isub) =     &
     &                           ppm_recvbuffers(ibuffer)
                              ibuffer = ibuffer + 1
                              fdata(5,imesh,jmesh,kmesh,isub) =     &
     &                           ppm_recvbuffers(ibuffer)
#elif  __KIND == __SINGLE_PRECISION_COMPLEX
                              fdata(1,imesh,jmesh,kmesh,isub) =     &
     &                           CMPLX(ppm_recvbuffers(ibuffer-1),  &
     &                           ppm_recvbuffers(ibuffer),ppm_kind_single)
                              ibuffer = ibuffer + 2
                              fdata(2,imesh,jmesh,kmesh,isub) =     &
     &                           CMPLX(ppm_recvbuffers(ibuffer-1),  &
     &                           ppm_recvbuffers(ibuffer),ppm_kind_single)
                              ibuffer = ibuffer + 2
                              fdata(3,imesh,jmesh,kmesh,isub) =     &
     &                           CMPLX(ppm_recvbuffers(ibuffer-1),  &
     &                           ppm_recvbuffers(ibuffer),ppm_kind_single)
                              ibuffer = ibuffer + 2
                              fdata(4,imesh,jmesh,kmesh,isub) =     &
     &                           CMPLX(ppm_recvbuffers(ibuffer-1),  &
     &                           ppm_recvbuffers(ibuffer),ppm_kind_single)
                              ibuffer = ibuffer + 2
                              fdata(5,imesh,jmesh,kmesh,isub) =     &
     &                           CMPLX(ppm_recvbuffers(ibuffer-1),  &
     &                           ppm_recvbuffers(ibuffer),ppm_kind_single)
#elif  __KIND == __DOUBLE_PRECISION_COMPLEX
                              fdata(1,imesh,jmesh,kmesh,isub) =     &
     &                           CMPLX(ppm_recvbuffers(ibuffer-1),  &
     &                           ppm_recvbuffers(ibuffer),ppm_kind_double)
                              ibuffer = ibuffer + 2
                              fdata(2,imesh,jmesh,kmesh,isub) =     &
     &                           CMPLX(ppm_recvbuffers(ibuffer-1),  &
     &                           ppm_recvbuffers(ibuffer),ppm_kind_double)
                              ibuffer = ibuffer + 2
                              fdata(3,imesh,jmesh,kmesh,isub) =     &
     &                           CMPLX(ppm_recvbuffers(ibuffer-1),  &
     &                           ppm_recvbuffers(ibuffer),ppm_kind_double)
                              ibuffer = ibuffer + 2
                              fdata(4,imesh,jmesh,kmesh,isub) =     &
     &                           CMPLX(ppm_recvbuffers(ibuffer-1),  &
     &                           ppm_recvbuffers(ibuffer),ppm_kind_double)
                              ibuffer = ibuffer + 2
                              fdata(5,imesh,jmesh,kmesh,isub) =     &
     &                           CMPLX(ppm_recvbuffers(ibuffer-1),  &
     &                           ppm_recvbuffers(ibuffer),ppm_kind_double)
#elif  __KIND == __INTEGER
                              fdata(1,imesh,jmesh,kmesh,isub) =     &
     &                           INT(ppm_recvbuffers(ibuffer))
                              ibuffer = ibuffer + 1
                              fdata(2,imesh,jmesh,kmesh,isub) =     &
     &                           INT(ppm_recvbuffers(ibuffer))
                              ibuffer = ibuffer + 1
                              fdata(3,imesh,jmesh,kmesh,isub) =     &
     &                           INT(ppm_recvbuffers(ibuffer))
                              ibuffer = ibuffer + 1
                              fdata(4,imesh,jmesh,kmesh,isub) =     &
     &                           INT(ppm_recvbuffers(ibuffer))
                              ibuffer = ibuffer + 1
                              fdata(5,imesh,jmesh,kmesh,isub) =     &
     &                           INT(ppm_recvbuffers(ibuffer))
#elif  __KIND == __LOGICAL
                              IF (ppm_recvbuffers(ibuffer) .GT.     &
     &                           (1.0_ppm_kind_single-ppm_myepss)) THEN
                                 fdata(1,imesh,jmesh,kmesh,isub) = .TRUE.
                              ELSE
                                 fdata(1,imesh,jmesh,kmesh,isub) = .FALSE.
                              ENDIF 
                              ibuffer = ibuffer + 1
                              IF (ppm_recvbuffers(ibuffer) .GT.     &
     &                           (1.0_ppm_kind_single-ppm_myepss)) THEN
                                 fdata(2,imesh,jmesh,kmesh,isub) = .TRUE.
                              ELSE
                                 fdata(2,imesh,jmesh,kmesh,isub) = .FALSE.
                              ENDIF 
                              ibuffer = ibuffer + 1
                              IF (ppm_recvbuffers(ibuffer) .GT.     &
     &                           (1.0_ppm_kind_single-ppm_myepss)) THEN
                                 fdata(3,imesh,jmesh,kmesh,isub) = .TRUE.
                              ELSE
                                 fdata(3,imesh,jmesh,kmesh,isub) = .FALSE.
                              ENDIF 
                              ibuffer = ibuffer + 1
                              IF (ppm_recvbuffers(ibuffer) .GT.     &
     &                           (1.0_ppm_kind_single-ppm_myepss)) THEN
                                 fdata(4,imesh,jmesh,kmesh,isub) = .TRUE.
                              ELSE
                                 fdata(4,imesh,jmesh,kmesh,isub) = .FALSE.
                              ENDIF 
                              ibuffer = ibuffer + 1
                              IF (ppm_recvbuffers(ibuffer) .GT.     &
     &                           (1.0_ppm_kind_single-ppm_myepss)) THEN
                                 fdata(5,imesh,jmesh,kmesh,isub) = .TRUE.
                              ELSE
                                 fdata(5,imesh,jmesh,kmesh,isub) = .FALSE.
                              ENDIF 
#endif
                           ENDDO
                          ENDDO
                        ENDDO
                     ELSEIF (rtype .EQ. ppm_param_pop_add) THEN
                        DO kmesh=zlo,zhi
                          DO jmesh=ylo,yhi
                           DO imesh=xlo,xhi
#if    __KIND == __SINGLE_PRECISION_COMPLEX | __KIND == __DOUBLE_PRECISION_COMPLEX
                              ibuffer = ibuffer + 2
#else
                              ibuffer = ibuffer + 1
#endif
#if    __KIND == __DOUBLE_PRECISION
                              fdata(1,imesh,jmesh,kmesh,isub) =    &
     &                           fdata(1,imesh,jmesh,kmesh,isub) + &
     &                           REAL(ppm_recvbuffers(ibuffer),ppm_kind_double)
                              ibuffer = ibuffer + 1
                              fdata(2,imesh,jmesh,kmesh,isub) =    &
     &                           fdata(2,imesh,jmesh,kmesh,isub) + &
     &                           REAL(ppm_recvbuffers(ibuffer),ppm_kind_double)
                              ibuffer = ibuffer + 1
                              fdata(3,imesh,jmesh,kmesh,isub) =    &
     &                           fdata(3,imesh,jmesh,kmesh,isub) + &
     &                           REAL(ppm_recvbuffers(ibuffer),ppm_kind_double)
                              ibuffer = ibuffer + 1
                              fdata(4,imesh,jmesh,kmesh,isub) =    &
     &                           fdata(4,imesh,jmesh,kmesh,isub) + &
     &                           REAL(ppm_recvbuffers(ibuffer),ppm_kind_double)
                              ibuffer = ibuffer + 1
                              fdata(5,imesh,jmesh,kmesh,isub) =    &
     &                           fdata(5,imesh,jmesh,kmesh,isub) + &
     &                           REAL(ppm_recvbuffers(ibuffer),ppm_kind_double)
#elif  __KIND == __SINGLE_PRECISION
                              fdata(1,imesh,jmesh,kmesh,isub) =    &
     &                           fdata(1,imesh,jmesh,kmesh,isub) + &
     &                           ppm_recvbuffers(ibuffer)
                              ibuffer = ibuffer + 1
                              fdata(2,imesh,jmesh,kmesh,isub) =    &
     &                           fdata(2,imesh,jmesh,kmesh,isub) + &
     &                           ppm_recvbuffers(ibuffer)
                              ibuffer = ibuffer + 1
                              fdata(3,imesh,jmesh,kmesh,isub) =    &
     &                           fdata(3,imesh,jmesh,kmesh,isub) + &
     &                           ppm_recvbuffers(ibuffer)
                              ibuffer = ibuffer + 1
                              fdata(4,imesh,jmesh,kmesh,isub) =    &
     &                           fdata(4,imesh,jmesh,kmesh,isub) + &
     &                           ppm_recvbuffers(ibuffer)
                              ibuffer = ibuffer + 1
                              fdata(5,imesh,jmesh,kmesh,isub) =    &
     &                           fdata(5,imesh,jmesh,kmesh,isub) + &
     &                           ppm_recvbuffers(ibuffer)
#elif  __KIND == __SINGLE_PRECISION_COMPLEX
                              fdata(1,imesh,jmesh,kmesh,isub) =    &
     &                           fdata(1,imesh,jmesh,kmesh,isub) + &
     &                           CMPLX(ppm_recvbuffers(ibuffer-1),  &
     &                           ppm_recvbuffers(ibuffer),ppm_kind_single)
                              ibuffer = ibuffer + 2
                              fdata(2,imesh,jmesh,kmesh,isub) =    &
     &                           fdata(2,imesh,jmesh,kmesh,isub) + &
     &                           CMPLX(ppm_recvbuffers(ibuffer-1),  &
     &                           ppm_recvbuffers(ibuffer),ppm_kind_single)
                              ibuffer = ibuffer + 2
                              fdata(3,imesh,jmesh,kmesh,isub) =    &
     &                           fdata(3,imesh,jmesh,kmesh,isub) + &
     &                           CMPLX(ppm_recvbuffers(ibuffer-1),  &
     &                           ppm_recvbuffers(ibuffer),ppm_kind_single)
                              ibuffer = ibuffer + 2
                              fdata(4,imesh,jmesh,kmesh,isub) =    &
     &                           fdata(4,imesh,jmesh,kmesh,isub) + &
     &                           CMPLX(ppm_recvbuffers(ibuffer-1),  &
     &                           ppm_recvbuffers(ibuffer),ppm_kind_single)
                              ibuffer = ibuffer + 2
                              fdata(5,imesh,jmesh,kmesh,isub) =    &
     &                           fdata(5,imesh,jmesh,kmesh,isub) + &
     &                           CMPLX(ppm_recvbuffers(ibuffer-1),  &
     &                           ppm_recvbuffers(ibuffer),ppm_kind_single)
#elif  __KIND == __DOUBLE_PRECISION_COMPLEX
                              fdata(1,imesh,jmesh,kmesh,isub) =    &
     &                           fdata(1,imesh,jmesh,kmesh,isub) + &
     &                           CMPLX(ppm_recvbuffers(ibuffer-1),  &
     &                           ppm_recvbuffers(ibuffer),ppm_kind_double)
                              ibuffer = ibuffer + 2
                              fdata(2,imesh,jmesh,kmesh,isub) =    &
     &                           fdata(2,imesh,jmesh,kmesh,isub) + &
     &                           CMPLX(ppm_recvbuffers(ibuffer-1),  &
     &                           ppm_recvbuffers(ibuffer),ppm_kind_double)
                              ibuffer = ibuffer + 2
                              fdata(3,imesh,jmesh,kmesh,isub) =    &
     &                           fdata(3,imesh,jmesh,kmesh,isub) + &
     &                           CMPLX(ppm_recvbuffers(ibuffer-1),  &
     &                           ppm_recvbuffers(ibuffer),ppm_kind_double)
                              ibuffer = ibuffer + 2
                              fdata(4,imesh,jmesh,kmesh,isub) =    &
     &                           fdata(4,imesh,jmesh,kmesh,isub) + &
     &                           CMPLX(ppm_recvbuffers(ibuffer-1),  &
     &                           ppm_recvbuffers(ibuffer),ppm_kind_double)
                              ibuffer = ibuffer + 2
                              fdata(5,imesh,jmesh,kmesh,isub) =    &
     &                           fdata(5,imesh,jmesh,kmesh,isub) + &
     &                           CMPLX(ppm_recvbuffers(ibuffer-1),  &
     &                           ppm_recvbuffers(ibuffer),ppm_kind_double)
#elif  __KIND == __INTEGER
                              fdata(1,imesh,jmesh,kmesh,isub) =    &
     &                           fdata(1,imesh,jmesh,kmesh,isub) + &
     &                           INT(ppm_recvbuffers(ibuffer))
                              ibuffer = ibuffer + 1
                              fdata(2,imesh,jmesh,kmesh,isub) =    &
     &                           fdata(2,imesh,jmesh,kmesh,isub) + &
     &                           INT(ppm_recvbuffers(ibuffer))
                              ibuffer = ibuffer + 1
                              fdata(3,imesh,jmesh,kmesh,isub) =    &
     &                           fdata(3,imesh,jmesh,kmesh,isub) + &
     &                           INT(ppm_recvbuffers(ibuffer))
                              ibuffer = ibuffer + 1
                              fdata(4,imesh,jmesh,kmesh,isub) =    &
     &                           fdata(4,imesh,jmesh,kmesh,isub) + &
     &                           INT(ppm_recvbuffers(ibuffer))
                              ibuffer = ibuffer + 1
                              fdata(5,imesh,jmesh,kmesh,isub) =    &
     &                           fdata(5,imesh,jmesh,kmesh,isub) + &
     &                           INT(ppm_recvbuffers(ibuffer))
#elif  __KIND == __LOGICAL
                              IF (ppm_recvbuffers(ibuffer) .GT.     &
     &                           (1.0_ppm_kind_single-ppm_myepss)) THEN
                                 fdata(1,imesh,jmesh,kmesh,isub) =        &
     &                             (fdata(1,imesh,jmesh,kmesh,isub) .AND..TRUE.)
                              ELSE
                                 fdata(1,imesh,jmesh,kmesh,isub) =        &
     &                            (fdata(1,imesh,jmesh,kmesh,isub) .AND..FALSE.)
                              ENDIF 
                              ibuffer = ibuffer + 1
                              IF (ppm_recvbuffers(ibuffer) .GT.     &
     &                           (1.0_ppm_kind_single-ppm_myepss)) THEN
                                 fdata(2,imesh,jmesh,kmesh,isub) =        &
     &                             (fdata(2,imesh,jmesh,kmesh,isub) .AND..TRUE.)
                              ELSE
                                 fdata(2,imesh,jmesh,kmesh,isub) =        &
     &                            (fdata(2,imesh,jmesh,kmesh,isub) .AND..FALSE.)
                              ENDIF 
                              ibuffer = ibuffer + 1
                              IF (ppm_recvbuffers(ibuffer) .GT.     &
     &                           (1.0_ppm_kind_single-ppm_myepss)) THEN
                                 fdata(3,imesh,jmesh,kmesh,isub) =        &
     &                             (fdata(3,imesh,jmesh,kmesh,isub) .AND..TRUE.)
                              ELSE
                                 fdata(3,imesh,jmesh,kmesh,isub) =        &
     &                             (fdata(3,imesh,jmesh,kmesh,isub).AND..FALSE.)
                              ENDIF 
                              ibuffer = ibuffer + 1
                              IF (ppm_recvbuffers(ibuffer) .GT.     &
     &                           (1.0_ppm_kind_single-ppm_myepss)) THEN
                                 fdata(4,imesh,jmesh,kmesh,isub) =        &
     &                             (fdata(4,imesh,jmesh,kmesh,isub).AND..TRUE.)
                              ELSE
                                 fdata(4,imesh,jmesh,kmesh,isub) =        &
     &                             (fdata(4,imesh,jmesh,kmesh,isub).AND..FALSE.)
                              ENDIF 
                              ibuffer = ibuffer + 1
                              IF (ppm_recvbuffers(ibuffer) .GT.     &
     &                           (1.0_ppm_kind_single-ppm_myepss)) THEN
                                 fdata(5,imesh,jmesh,kmesh,isub) =        &
     &                             (fdata(5,imesh,jmesh,kmesh,isub).AND..TRUE.)
                              ELSE
                                 fdata(5,imesh,jmesh,kmesh,isub) =        &
     &                             (fdata(5,imesh,jmesh,kmesh,isub).AND..FALSE.)
                              ENDIF 
#endif
                           ENDDO
                          ENDDO
                        ENDDO
                     ENDIF
                  !-------------------------------------------------------------
                  !  For edim.GT.5 the vector length will be edim !!
                  !-------------------------------------------------------------
                  ELSE
                     IF (rtype .EQ. ppm_param_pop_replace) THEN
                        DO kmesh=zlo,zhi
                          DO jmesh=ylo,yhi
                           DO imesh=xlo,xhi
                              DO k=1,edim
#if    __KIND == __SINGLE_PRECISION_COMPLEX | __KIND == __DOUBLE_PRECISION_COMPLEX
                                 ibuffer = ibuffer + 2
#else
                                 ibuffer = ibuffer + 1
#endif
#if    __KIND == __DOUBLE_PRECISION
                                 fdata(k,imesh,jmesh,kmesh,isub) =     &
     &                            REAL(ppm_recvbuffers(ibuffer),ppm_kind_double)
#elif  __KIND == __SINGLE_PRECISION
                                 fdata(k,imesh,jmesh,kmesh,isub) =     &
     &                              ppm_recvbuffers(ibuffer)
#elif  __KIND == __SINGLE_PRECISION_COMPLEX
                                 fdata(k,imesh,jmesh,kmesh,isub) =     &
     &                              CMPLX(ppm_recvbuffers(ibuffer-1),  &
     &                              ppm_recvbuffers(ibuffer),ppm_kind_single)
#elif  __KIND == __DOUBLE_PRECISION_COMPLEX
                                 fdata(k,imesh,jmesh,kmesh,isub) =     &
     &                              CMPLX(ppm_recvbuffers(ibuffer-1),  &
     &                              ppm_recvbuffers(ibuffer),ppm_kind_double)
#elif  __KIND == __INTEGER
                                 fdata(k,imesh,jmesh,kmesh,isub) =     &
     &                              INT(ppm_recvbuffers(ibuffer))
#elif  __KIND == __LOGICAL
                                 IF (ppm_recvbuffers(ibuffer) .GT.     &
     &                              (1.0_ppm_kind_single-ppm_myepss)) THEN
                                    fdata(k,imesh,jmesh,kmesh,isub) = .TRUE.
                                 ELSE
                                    fdata(k,imesh,jmesh,kmesh,isub) = .FALSE.
                                 ENDIF 
#endif
                              ENDDO
                           ENDDO
                          ENDDO
                        ENDDO
                     ELSEIF (rtype .EQ. ppm_param_pop_add) THEN
                        DO kmesh=zlo,zhi
                          DO jmesh=ylo,yhi
                           DO imesh=xlo,xhi
                              DO k=1,edim
#if    __KIND == __SINGLE_PRECISION_COMPLEX | __KIND == __DOUBLE_PRECISION_COMPLEX
                                 ibuffer = ibuffer + 2
#else
                                 ibuffer = ibuffer + 1
#endif
#if    __KIND == __DOUBLE_PRECISION
                                 fdata(k,imesh,jmesh,kmesh,isub) =    &
     &                            fdata(k,imesh,jmesh,kmesh,isub) + &
     &                            REAL(ppm_recvbuffers(ibuffer),ppm_kind_double)
#elif  __KIND == __SINGLE_PRECISION
                                 fdata(k,imesh,jmesh,kmesh,isub) =    &
     &                              fdata(k,imesh,jmesh,kmesh,isub) + &
     &                              ppm_recvbuffers(ibuffer)
#elif  __KIND == __SINGLE_PRECISION_COMPLEX
                                 fdata(k,imesh,jmesh,kmesh,isub) =    &
     &                              fdata(k,imesh,jmesh,kmesh,isub) + &
     &                              CMPLX(ppm_recvbuffers(ibuffer-1),  &
     &                              ppm_recvbuffers(ibuffer),ppm_kind_single)
#elif  __KIND == __DOUBLE_PRECISION_COMPLEX
                                 fdata(k,imesh,jmesh,kmesh,isub) =    &
     &                              fdata(k,imesh,jmesh,kmesh,isub) + &
     &                              CMPLX(ppm_recvbuffers(ibuffer-1),  &
     &                              ppm_recvbuffers(ibuffer),ppm_kind_double)
#elif  __KIND == __INTEGER
                                 fdata(k,imesh,jmesh,kmesh,isub) =    &
     &                              fdata(k,imesh,jmesh,kmesh,isub) + &
     &                              INT(ppm_recvbuffers(ibuffer))
#elif  __KIND == __LOGICAL
                                 IF (ppm_recvbuffers(ibuffer) .GT.     &
     &                              (1.0_ppm_kind_single-ppm_myepss)) THEN
                                    fdata(k,imesh,jmesh,kmesh,isub) =        &
     &                                  (fdata(k,imesh,jmesh,kmesh,isub).AND.&
     &                                  .TRUE.)
                                 ELSE
                                    fdata(k,imesh,jmesh,kmesh,isub) =        &
     &                                  (fdata(k,imesh,jmesh,kmesh,isub).AND.&
     &                                  .FALSE.)
                                 ENDIF 
#endif
                              ENDDO
                           ENDDO
                          ENDDO
                        ENDDO
                     ENDIF           ! rtype
                  ENDIF            ! lda = ...
               !----------------------------------------------------------------
               !  With mask: This will not vectorize so no need to unroll
               !----------------------------------------------------------------
               ELSE
                  DO kmesh=zlo,zhi
                    DO jmesh=ylo,yhi
                     DO imesh=xlo,xhi
                        IF (mask(imesh,jmesh,kmesh,isub)) THEN
                           DO k=1,edim
#if    __KIND == __SINGLE_PRECISION_COMPLEX | __KIND == __DOUBLE_PRECISION_COMPLEX
                              ibuffer = ibuffer + 2
#else
                              ibuffer = ibuffer + 1
#endif
                              IF (rtype .EQ. ppm_param_pop_replace) THEN
#if    __KIND == __DOUBLE_PRECISION
                                 fdata(k,imesh,jmesh,kmesh,isub) =     &
     &                            REAL(ppm_recvbuffers(ibuffer),ppm_kind_double)
#elif  __KIND == __SINGLE_PRECISION
                                 fdata(k,imesh,jmesh,kmesh,isub) =     &
     &                              ppm_recvbuffers(ibuffer)
#elif  __KIND == __SINGLE_PRECISION_COMPLEX
                                 fdata(k,imesh,jmesh,kmesh,isub) =     &
     &                              CMPLX(ppm_recvbuffers(ibuffer-1),  &
     &                              ppm_recvbuffers(ibuffer),ppm_kind_single)
#elif  __KIND == __DOUBLE_PRECISION_COMPLEX
                                 fdata(k,imesh,jmesh,kmesh,isub) =     &
     &                              CMPLX(ppm_recvbuffers(ibuffer-1),  &
     &                              ppm_recvbuffers(ibuffer),ppm_kind_double)
#elif  __KIND == __INTEGER
                                 fdata(k,imesh,jmesh,kmesh,isub) =     &
     &                              INT(ppm_recvbuffers(ibuffer))
#elif  __KIND == __LOGICAL
                                 IF (ppm_recvbuffers(ibuffer) .GT.     &
     &                              (1.0_ppm_kind_single-ppm_myepss)) THEN
                                    fdata(k,imesh,jmesh,kmesh,isub) = .TRUE.
                                 ELSE
                                    fdata(k,imesh,jmesh,kmesh,isub) = .FALSE.
                                 ENDIF 
#endif
                              ELSEIF (rtype .EQ. ppm_param_pop_add) THEN
#if    __KIND == __DOUBLE_PRECISION
                                 fdata(k,imesh,jmesh,kmesh,isub) =    &
     &                            fdata(k,imesh,jmesh,kmesh,isub) + &
     &                            REAL(ppm_recvbuffers(ibuffer),ppm_kind_double)
#elif  __KIND == __SINGLE_PRECISION
                                 fdata(k,imesh,jmesh,kmesh,isub) =    &
     &                              fdata(k,imesh,jmesh,kmesh,isub) + &
     &                              ppm_recvbuffers(ibuffer)
#elif  __KIND == __SINGLE_PRECISION_COMPLEX
                                 fdata(k,imesh,jmesh,kmesh,isub) =    &
     &                              fdata(k,imesh,jmesh,kmesh,isub) + &
     &                              CMPLX(ppm_recvbuffers(ibuffer-1),  &
     &                              ppm_recvbuffers(ibuffer),ppm_kind_single)
#elif  __KIND == __DOUBLE_PRECISION_COMPLEX
                                 fdata(k,imesh,jmesh,kmesh,isub) =    &
     &                              fdata(k,imesh,jmesh,kmesh,isub) + &
     &                              CMPLX(ppm_recvbuffers(ibuffer-1),  &
     &                              ppm_recvbuffers(ibuffer),ppm_kind_double)
#elif  __KIND == __INTEGER
                                 fdata(k,imesh,jmesh,kmesh,isub) =    &
     &                              fdata(k,imesh,jmesh,kmesh,isub) + &
     &                              INT(ppm_recvbuffers(ibuffer))
#elif  __KIND == __LOGICAL
                                 IF (ppm_recvbuffers(ibuffer) .GT.     &
     &                              (1.0_ppm_kind_single-ppm_myepss)) THEN
                                    fdata(k,imesh,jmesh,kmesh,isub) =        &
     &                                  (fdata(k,imesh,jmesh,kmesh,isub).AND.&
     &                                  .TRUE.)
                                 ELSE
                                    fdata(k,imesh,jmesh,kmesh,isub) =        &
     &                                  (fdata(k,imesh,jmesh,kmesh,isub).AND.&
     &                                  .FALSE.)
                                 ENDIF 
#endif
                              ENDIF
                           ENDDO
                        ELSE
                           !----------------------------------------------------
                           !  If ldo.EQ..FALSE. skip this point
                           !----------------------------------------------------
#if    __KIND == __SINGLE_PRECISION_COMPLEX | __KIND == __DOUBLE_PRECISION_COMPLEX
                           ibuffer = ibuffer + (edim*2)
#else
                           ibuffer = ibuffer + edim
#endif
                        ENDIF   ! mask.EQ..TRUE.
                     ENDDO
                    ENDDO
                  ENDDO
               ENDIF            ! PRESENT(mask)
#elif __DIM == __SFIELD
               IF (.NOT.PRESENT(mask)) THEN
                  IF (rtype .EQ. ppm_param_pop_replace) THEN
                     DO kmesh=zlo,zhi
                       DO jmesh=ylo,yhi
                        DO imesh=xlo,xhi
#if    __KIND == __SINGLE_PRECISION_COMPLEX | __KIND == __DOUBLE_PRECISION_COMPLEX
                           ibuffer = ibuffer + 2
#else
                           ibuffer = ibuffer + 1
#endif
#if    __KIND == __DOUBLE_PRECISION
                           fdata(imesh,jmesh,kmesh,isub) =     &
     &                        REAL(ppm_recvbuffers(ibuffer),ppm_kind_double)
#elif  __KIND == __SINGLE_PRECISION
                           fdata(imesh,jmesh,kmesh,isub) =     &
     &                        ppm_recvbuffers(ibuffer)
#elif  __KIND == __SINGLE_PRECISION_COMPLEX
                           fdata(imesh,jmesh,kmesh,isub) =     &
     &                        CMPLX(ppm_recvbuffers(ibuffer-1),  &
     &                        ppm_recvbuffers(ibuffer),ppm_kind_single)
#elif  __KIND == __DOUBLE_PRECISION_COMPLEX
                           fdata(imesh,jmesh,kmesh,isub) =     &
     &                        CMPLX(ppm_recvbuffers(ibuffer-1),  &
     &                        ppm_recvbuffers(ibuffer),ppm_kind_double)
#elif  __KIND == __INTEGER
                           fdata(imesh,jmesh,kmesh,isub) =     &
     &                        INT(ppm_recvbuffers(ibuffer))
#elif  __KIND == __LOGICAL
                           IF (ppm_recvbuffers(ibuffer) .GT.     &
     &                        (1.0_ppm_kind_single-ppm_myepss)) THEN
                              fdata(imesh,jmesh,kmesh,isub) = .TRUE.
                           ELSE
                              fdata(imesh,jmesh,kmesh,isub) = .FALSE.
                           ENDIF 
#endif
                        ENDDO
                       ENDDO
                     ENDDO
                  ELSEIF (rtype .EQ. ppm_param_pop_add) THEN
                     DO kmesh=zlo,zhi
                       DO jmesh=ylo,yhi
                        DO imesh=xlo,xhi
#if    __KIND == __SINGLE_PRECISION_COMPLEX | __KIND == __DOUBLE_PRECISION_COMPLEX
                           ibuffer = ibuffer + 2
#else
                           ibuffer = ibuffer + 1
#endif
#if    __KIND == __DOUBLE_PRECISION
                           fdata(imesh,jmesh,kmesh,isub) =    &
     &                        fdata(imesh,jmesh,kmesh,isub) + &
     &                        REAL(ppm_recvbuffers(ibuffer),ppm_kind_double)
#elif  __KIND == __SINGLE_PRECISION
                           fdata(imesh,jmesh,kmesh,isub) =    &
     &                        fdata(imesh,jmesh,kmesh,isub) + &
     &                        ppm_recvbuffers(ibuffer)
#elif  __KIND == __SINGLE_PRECISION_COMPLEX
                           fdata(imesh,jmesh,kmesh,isub) =    &
     &                        fdata(imesh,jmesh,kmesh,isub) + &
     &                        CMPLX(ppm_recvbuffers(ibuffer-1),  &
     &                        ppm_recvbuffers(ibuffer),ppm_kind_single)
#elif  __KIND == __DOUBLE_PRECISION_COMPLEX
                           fdata(imesh,jmesh,kmesh,isub) =    &
     &                        fdata(imesh,jmesh,kmesh,isub) + &
     &                        CMPLX(ppm_recvbuffers(ibuffer-1),  &
     &                        ppm_recvbuffers(ibuffer),ppm_kind_double)
#elif  __KIND == __INTEGER
                           fdata(imesh,jmesh,kmesh,isub) =    &
     &                        fdata(imesh,jmesh,kmesh,isub) + &
     &                        INT(ppm_recvbuffers(ibuffer))
#elif  __KIND == __LOGICAL
                           IF (ppm_recvbuffers(ibuffer) .GT.     &
     &                        (1.0_ppm_kind_single-ppm_myepss)) THEN
                              fdata(imesh,jmesh,kmesh,isub) =        &
     &                            (fdata(imesh,jmesh,kmesh,isub) .AND..TRUE.)
                           ELSE
                              fdata(imesh,jmesh,kmesh,isub) =        &
     &                            (fdata(imesh,jmesh,kmesh,isub) .AND..FALSE.)
                           ENDIF 
#endif
                        ENDDO
                       ENDDO
                     ENDDO
                  ENDIF
              ELSE            ! IF PRESENT(mask)
                  DO kmesh=zlo,zhi
                    DO jmesh=ylo,yhi
                     DO imesh=xlo,xhi
                        IF (mask(imesh,jmesh,kmesh,isub)) THEN
#if    __KIND == __SINGLE_PRECISION_COMPLEX | __KIND == __DOUBLE_PRECISION_COMPLEX
                           ibuffer = ibuffer + 2
#else
                           ibuffer = ibuffer + 1
#endif
                           IF (rtype .EQ. ppm_param_pop_replace) THEN
#if    __KIND == __DOUBLE_PRECISION
                              fdata(imesh,jmesh,kmesh,isub) =     &
     &                           REAL(ppm_recvbuffers(ibuffer),ppm_kind_double)
#elif  __KIND == __SINGLE_PRECISION
                              fdata(imesh,jmesh,kmesh,isub) =     &
     &                           ppm_recvbuffers(ibuffer)
#elif  __KIND == __SINGLE_PRECISION_COMPLEX
                              fdata(imesh,jmesh,kmesh,isub) =     &
     &                           CMPLX(ppm_recvbuffers(ibuffer-1),  &
     &                           ppm_recvbuffers(ibuffer),ppm_kind_single)
#elif  __KIND == __DOUBLE_PRECISION_COMPLEX
                              fdata(imesh,jmesh,kmesh,isub) =     &
     &                           CMPLX(ppm_recvbuffers(ibuffer-1),  &
     &                           ppm_recvbuffers(ibuffer),ppm_kind_double)
#elif  __KIND == __INTEGER
                              fdata(imesh,jmesh,kmesh,isub) =     &
     &                           INT(ppm_recvbuffers(ibuffer))
#elif  __KIND == __LOGICAL
                              IF (ppm_recvbuffers(ibuffer) .GT.     &
     &                           (1.0_ppm_kind_single-ppm_myepss)) THEN
                                 fdata(imesh,jmesh,kmesh,isub) = .TRUE.
                              ELSE
                                 fdata(imesh,jmesh,kmesh,isub) = .FALSE.
                              ENDIF 
#endif
                           ELSEIF (rtype .EQ. ppm_param_pop_add) THEN
#if    __KIND == __DOUBLE_PRECISION
                              fdata(imesh,jmesh,kmesh,isub) =    &
     &                           fdata(imesh,jmesh,kmesh,isub) + &
     &                           REAL(ppm_recvbuffers(ibuffer),ppm_kind_double)
#elif  __KIND == __SINGLE_PRECISION
                              fdata(imesh,jmesh,kmesh,isub) =    &
     &                           fdata(imesh,jmesh,kmesh,isub) + &
     &                           ppm_recvbuffers(ibuffer)
#elif  __KIND == __SINGLE_PRECISION_COMPLEX
                              fdata(imesh,jmesh,kmesh,isub) =    &
     &                           fdata(imesh,jmesh,kmesh,isub) + &
     &                           CMPLX(ppm_recvbuffers(ibuffer-1),  &
     &                           ppm_recvbuffers(ibuffer),ppm_kind_single)
#elif  __KIND == __DOUBLE_PRECISION_COMPLEX
                              fdata(imesh,jmesh,kmesh,isub) =    &
     &                           fdata(imesh,jmesh,kmesh,isub) + &
     &                           CMPLX(ppm_recvbuffers(ibuffer-1),  &
     &                           ppm_recvbuffers(ibuffer),ppm_kind_double)
#elif  __KIND == __INTEGER
                              fdata(imesh,jmesh,kmesh,isub) =    &
     &                           fdata(imesh,jmesh,kmesh,isub) + &
     &                           INT(ppm_recvbuffers(ibuffer))
#elif  __KIND == __LOGICAL
                              IF (ppm_recvbuffers(ibuffer) .GT.     &
     &                           (1.0_ppm_kind_single-ppm_myepss)) THEN
                                 fdata(imesh,jmesh,kmesh,isub) =        &
     &                             (fdata(imesh,jmesh,kmesh,isub).AND..TRUE.)
                              ELSE
                                 fdata(imesh,jmesh,kmesh,isub) =        &
     &                             (fdata(imesh,jmesh,kmesh,isub).AND..FALSE.)
                              ENDIF 
#endif
                           ENDIF
                        ELSE
                           !----------------------------------------------------
                           !  If ldo.EQ..FALSE. skip this point
                           !----------------------------------------------------
#if    __KIND == __SINGLE_PRECISION_COMPLEX | __KIND == __DOUBLE_PRECISION_COMPLEX
                           ibuffer = ibuffer + 2
#else
                           ibuffer = ibuffer + 1
#endif
                        ENDIF    ! mask.EQ..TRUE.
                     ENDDO
                    ENDDO
                  ENDDO
               ENDIF          ! PRESENT(mask)
#endif
            ENDDO             ! ppm_precvbuffer
         ENDDO                ! ppm_nrecvlist
      ENDIF                   ! ppm_kind

      !-------------------------------------------------------------------------
      !  Decrement the set counter
      !-------------------------------------------------------------------------
      ppm_buffer_set = ppm_buffer_set - 1  

      !-------------------------------------------------------------------------
      !  Deallocate the receive buffer if all sets have been poped
      !-------------------------------------------------------------------------
      IF (ppm_buffer_set .LT. 1) THEN
          iopt = ppm_param_dealloc
          IF (ppm_kind .EQ. ppm_kind_single) THEN
              CALL ppm_alloc(ppm_recvbuffers,ldu,iopt,info)
          ELSE
              CALL ppm_alloc(ppm_recvbufferd,ldu,iopt,info)
          ENDIF
          IF (info .NE. 0) THEN
              info = ppm_error_error
              CALL ppm_error(ppm_err_alloc,'ppm_map_field_pop_3d',     &
     &            'receive buffer PPM_RECVBUFFER',__LINE__,info)
          ENDIF
      ENDIF

      !-------------------------------------------------------------------------
      !  Deallocate inverse sub list
      !-------------------------------------------------------------------------
      iopt   = ppm_param_dealloc
      CALL ppm_alloc(invsublist,ldu,iopt,info)
      IF (info .NE. 0) THEN
          info = ppm_error_error
          CALL ppm_error(ppm_err_dealloc,'ppm_map_field_pop_3d',     &
     &        'inverse sub list INVSUBLIST',__LINE__,info)
      ENDIF

      !-------------------------------------------------------------------------
      !  Return 
      !-------------------------------------------------------------------------
 9999 CONTINUE
      CALL substop('ppm_map_field_pop_3d',t0,info)
      RETURN
#if    __DIM == __SFIELD
#if    __KIND == __SINGLE_PRECISION
      END SUBROUTINE ppm_map_field_pop_3d_sca_s
#elif  __KIND == __DOUBLE_PRECISION
      END SUBROUTINE ppm_map_field_pop_3d_sca_d
#elif  __KIND == __SINGLE_PRECISION_COMPLEX
      END SUBROUTINE ppm_map_field_pop_3d_sca_sc
#elif  __KIND == __DOUBLE_PRECISION_COMPLEX
      END SUBROUTINE ppm_map_field_pop_3d_sca_dc
#elif  __KIND == __INTEGER
      END SUBROUTINE ppm_map_field_pop_3d_sca_i
#elif  __KIND == __LOGICAL
      END SUBROUTINE ppm_map_field_pop_3d_sca_l
#endif

#elif  __DIM == __VFIELD
#if    __KIND == __SINGLE_PRECISION
      END SUBROUTINE ppm_map_field_pop_3d_vec_s
#elif  __KIND == __DOUBLE_PRECISION
      END SUBROUTINE ppm_map_field_pop_3d_vec_d
#elif  __KIND == __SINGLE_PRECISION_COMPLEX
      END SUBROUTINE ppm_map_field_pop_3d_vec_sc
#elif  __KIND == __DOUBLE_PRECISION_COMPLEX
      END SUBROUTINE ppm_map_field_pop_3d_vec_dc
#elif  __KIND == __INTEGER
      END SUBROUTINE ppm_map_field_pop_3d_vec_i
#elif  __KIND == __LOGICAL
      END SUBROUTINE ppm_map_field_pop_3d_vec_l
#endif
#endif
