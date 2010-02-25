      !-------------------------------------------------------------------------
      !  Subroutine   :               ppm_map_field_push_2d
      !-------------------------------------------------------------------------
      !
      !  Purpose      : This routine pushes field data onto the send buffer for
      !                 2d meshes.
      !
      !  Input        : fdata([:,]:,:,:)(F) field data. 1st index: lda,
      !                                     2nd-3th: mesh (i,j) relative
      !                                     to istart-1 of the sub (i.e.
      !                                     (i,j)=(1,1) corresponds to
      !                                     istart), 4th: isub
      !                                     1...nsublist (all subs on this
      !                                     processor). For scalar field
      !                                     data, the first index is
      !                                     omitted (the others shift
      !                                     accordingly).
      !                 lda             (I) the leading dimension of the fdata.
      !                                     lda=1 for the case of scalar
      !                                     data.
      !                 frommesh        (I) mesh ID of source (internal
      !                                     numbering)
      !                 mask(:,:,:)     (L) Logical mask. Only the mesh nodes
      !                                     for which this is .TRUE. will be
      !                                     mapped. OPTIONAL. If not given, all
      !                                     points are mapped. 1st-2nd
      !                                     index: mesh (i,j), 3rd: isub.
      !
      !  Input/output :
      !
      !  Output       : info            (I) return status. 0 on success.
      !
      !  Remarks      : The on-processor data is stored in the first part
      !                 of the buffer
      !
      !                 If the call to invert_list is too slow here, this
      !                 could be done by topo_store and stored for each
      !                 topology. At the expense of some extra memory.
      !
      !                 The fdata arrays are passed as POINTER and not
      !                 INTENT(IN) even though they are not changed or
      !                 reallocated inside this routine. The reason is that
      !                 LBOUND can be non-standard (i.e. .NE. 1) if we have
      !                 ghost layers and the correct LBOUND is only passed
      !                 along with the array if the latter is given the
      !                 POINTER attribute.
      !
      !                 When Ndata is determined, the mask is not
      !                 considered. This could lead to a too large
      !                 sendbuffer being allocated in the last push and
      !                 could be optimized later.
      !                 
      !  References   :
      !
      !  Revisions    :
      !-------------------------------------------------------------------------
      !  $Log: ppm_map_field_push_2d.f,v $
      !  Revision 1.1.1.1  2007/07/13 10:18:56  ivos
      !  CBL version of the PPM library
      !
      !  Revision 1.18  2005/03/10 01:40:32  ivos
      !  Empty buffer is now only reported for ppm_debug.GT.0 and is no
      !  longer logged to prevent huge log files.
      !
      !  Revision 1.17  2005/02/16 22:07:04  ivos
      !  Unrolled loops for lda=1,2,3,4,5 to allow vectorization.
      !  If mask is present, no unrolling is done since this case will
      !  not vectorize anyway.
      !
      !  Revision 1.16  2004/11/11 15:23:02  ivos
      !  Moved allocatable work data to the module.
      !
      !  Revision 1.15  2004/10/01 16:09:05  ivos
      !  Replaced REAL(ppm_kind_double) :: t0 with REAL(MK) t0.
      !
      !  Revision 1.14  2004/08/31 12:48:10  ivos
      !  Changed argument checks to use ppm_check_topoid and ppm_check_meshid.
      !
      !  Revision 1.13  2004/07/26 15:38:49  ivos
      !  Inserted missing USE statements to resolve undefined references
      !  at link stage.
      !
      !  Revision 1.12  2004/07/26 07:42:42  ivos
      !  Changed to single-interface modules and adapted all USE statements.
      !
      !  Revision 1.11  2004/07/19 11:01:52  ivos
      !  Overloaded field push and pop operations for scalar fields and
      !  added all changes to the module and the interface routines.
      !
      !  Revision 1.10  2004/04/15 10:18:33  ivos
      !  bugfix: changed mask from INTENT(IN) to POINTER as it can have
      !  LBOUND .LT. 1 (if ghost layer is masked).
      !
      !  Revision 1.9  2004/04/14 15:08:51  ivos
      !  Added OPTIONAL argument mask to allow selective mapping of only a
      !  subset of mesh points. Different masks can be used for different
      !  fields on the same mesh. Not tested yet.
      !
      !  Revision 1.8  2004/04/07 09:19:44  ivos
      !  Changed INTENT(IN) to POINTER for fdata because if we have ghosts
      !  LBOUND is .NE. 1 and the proper LBOUND is only passed for POINTER
      !  arguments, but not for INTENT(IN).
      !
      !  Revision 1.7  2004/04/05 10:45:29  ivos
      !  Translation to local sub indices now uses invert_list outside of the
      !  loop instead of searching isublist every time.
      !
      !  Revision 1.6  2004/04/02 15:22:23  ivos
      !  Added translation from global to local sub ID, since all sub IDs in
      !  the ppm_mesh_*sub lists are now global (because of ghosts).
      !
      !  Revision 1.5  2004/03/05 13:43:23  ivos
      !  bugfix: index errors corrected for COMPLEX version. COMPLEX version
      !  is now tested.
      !
      !  Revision 1.4  2004/02/24 08:46:55  ivos
      !  Added overloading for single complex and double complex data types.
      !
      !  Revision 1.3  2004/02/23 12:19:00  ivos
      !  Several bugs fixed. Tested on 2 processors with a scalar field.
      !  Added debug output in several places.
      !
      !  Revision 1.2  2004/02/23 08:56:39  ivos
      !  bugfix: necessary size of data array (for argument check) was 
      !  determined using global data instead of on-processor data. fixed.
      !
      !  Revision 1.1  2004/02/17 16:12:30  ivos
      !  Initial implementation. Not tested yet.
      !
      !-------------------------------------------------------------------------
      !  Parallel Particle Mesh Library (PPM)
      !  Institute of Computational Science
      !  ETH Zentrum, Hirschengraben 84
      !  CH-8092 Zurich, Switzerland
      !-------------------------------------------------------------------------

      ! ATTN: DIM is the dimension of the fdata array and not the space
      ! dimension ppm_dim!
#if    __DIM == __SFIELD
#if    __KIND == __SINGLE_PRECISION
      SUBROUTINE ppm_map_field_push_2d_sca_s(fdata,lda,frommesh,info,mask)
#elif  __KIND == __DOUBLE_PRECISION
      SUBROUTINE ppm_map_field_push_2d_sca_d(fdata,lda,frommesh,info,mask)
#elif  __KIND == __SINGLE_PRECISION_COMPLEX
      SUBROUTINE ppm_map_field_push_2d_sca_sc(fdata,lda,frommesh,info,mask)
#elif  __KIND == __DOUBLE_PRECISION_COMPLEX
      SUBROUTINE ppm_map_field_push_2d_sca_dc(fdata,lda,frommesh,info,mask)
#elif  __KIND == __INTEGER
      SUBROUTINE ppm_map_field_push_2d_sca_i(fdata,lda,frommesh,info,mask)
#elif  __KIND == __LOGICAL
      SUBROUTINE ppm_map_field_push_2d_sca_l(fdata,lda,frommesh,info,mask)
#endif 

#elif  __DIM == __VFIELD
#if    __KIND == __SINGLE_PRECISION
      SUBROUTINE ppm_map_field_push_2d_vec_s(fdata,lda,frommesh,info,mask)
#elif  __KIND == __DOUBLE_PRECISION
      SUBROUTINE ppm_map_field_push_2d_vec_d(fdata,lda,frommesh,info,mask)
#elif  __KIND == __SINGLE_PRECISION_COMPLEX
      SUBROUTINE ppm_map_field_push_2d_vec_sc(fdata,lda,frommesh,info,mask)
#elif  __KIND == __DOUBLE_PRECISION_COMPLEX
      SUBROUTINE ppm_map_field_push_2d_vec_dc(fdata,lda,frommesh,info,mask)
#elif  __KIND == __INTEGER
      SUBROUTINE ppm_map_field_push_2d_vec_i(fdata,lda,frommesh,info,mask)
#elif  __KIND == __LOGICAL
      SUBROUTINE ppm_map_field_push_2d_vec_l(fdata,lda,frommesh,info,mask)
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
      INTEGER , DIMENSION(:,:,:)   , POINTER         :: fdata
#elif __KIND == __LOGICAL
      LOGICAL , DIMENSION(:,:,:)   , POINTER         :: fdata
#elif __KIND == __SINGLE_PRECISION_COMPLEX | __KIND == __DOUBLE_PRECISION_COMPLEX
      COMPLEX(MK), DIMENSION(:,:,:), POINTER         :: fdata
#else
      REAL(MK), DIMENSION(:,:,:)   , POINTER         :: fdata
#endif

#elif __DIM == __VFIELD
#if   __KIND == __INTEGER
      INTEGER , DIMENSION(:,:,:,:)   , POINTER       :: fdata
#elif __KIND == __LOGICAL
      LOGICAL , DIMENSION(:,:,:,:)   , POINTER       :: fdata
#elif __KIND == __SINGLE_PRECISION_COMPLEX | __KIND == __DOUBLE_PRECISION_COMPLEX
      COMPLEX(MK), DIMENSION(:,:,:,:), POINTER       :: fdata
#else
      REAL(MK), DIMENSION(:,:,:,:)   , POINTER       :: fdata
#endif
#endif
      LOGICAL, DIMENSION(:,:,:),   POINTER, OPTIONAL :: mask
      INTEGER                        , INTENT(IN   ) :: lda,frommesh
      INTEGER                        , INTENT(  OUT) :: info
      !-------------------------------------------------------------------------
      !  Local variables 
      !-------------------------------------------------------------------------
      INTEGER, DIMENSION(2)   :: ldu,mofs
      INTEGER                 :: i,j,k,ibuffer,isub,imesh,jmesh,jsub
      INTEGER                 :: iopt,Ndata,xlo,xhi,ylo,yhi,fromtopo,ldb
      REAL(MK)                :: t0
      LOGICAL                 :: ldo
      CHARACTER(LEN=ppm_char) :: mesg
      !-------------------------------------------------------------------------
      !  Externals 
      !-------------------------------------------------------------------------
      
      !-------------------------------------------------------------------------
      !  Initialise 
      !-------------------------------------------------------------------------
      CALL substart('ppm_map_field_push_2d',t0,info)
      fromtopo = ppm_field_topoid

      !-------------------------------------------------------------------------
      !  Check arguments
      !-------------------------------------------------------------------------
      IF (ppm_debug .GT. 0) THEN
          CALL ppm_check_meshid(ppm_param_id_internal,frommesh,fromtopo,   &
     &        ldo,info)
          IF (.NOT. ldo) THEN
              info = ppm_error_error
              CALL ppm_error(ppm_err_argument,'ppm_map_field_push_2d',  &
     &            'source meshid out of range',__LINE__,info)
              GOTO 9999
          ENDIF
          xhi = 0
          yhi = 0
          DO i=1,ppm_nsublist(fromtopo)
              isub = ppm_isublist(i,fromtopo)
              IF (ppm_cart_mesh(frommesh,fromtopo)%nnodes(1,isub).GT.xhi) THEN
                  xhi = ppm_cart_mesh(frommesh,fromtopo)%nnodes(1,isub)
              ENDIF
              IF (ppm_cart_mesh(frommesh,fromtopo)%nnodes(2,isub).GT.yhi) THEN
                  yhi = ppm_cart_mesh(frommesh,fromtopo)%nnodes(2,isub)
              ENDIF
          ENDDO
#if   __DIM == __VFIELD
          IF (lda .LT. 1) THEN
              info = ppm_error_error
              CALL ppm_error(ppm_err_argument,'ppm_map_field_push_2d',  &
     &            'lda must be >=1 for vector data',__LINE__,info)
              GOTO 9999
          ENDIF
          IF (SIZE(fdata,4) .LT. ppm_nsublist(fromtopo)) THEN
              info = ppm_error_error
              CALL ppm_error(ppm_err_argument,'ppm_map_field_push_2d',  &
     &            'field data for some subs is missing',__LINE__,info)
              GOTO 9999
          ENDIF
          IF (SIZE(fdata,1) .LT. lda) THEN
              info = ppm_error_error
              CALL ppm_error(ppm_err_argument,'ppm_map_field_push_2d',  &
     &            'leading dimension of data does not match lda',__LINE__,info)
              GOTO 9999
          ENDIF
          IF (SIZE(fdata,2) .LT. xhi) THEN
              info = ppm_error_error
              CALL ppm_error(ppm_err_argument,'ppm_map_field_push_2d',  &
     &            'x dimension of field data does not match mesh',__LINE__,info)
              GOTO 9999
          ENDIF
          IF (SIZE(fdata,3) .LT. yhi) THEN
              info = ppm_error_error
              CALL ppm_error(ppm_err_argument,'ppm_map_field_push_2d',  &
     &            'y dimension of field data does not match mesh',__LINE__,info)
              GOTO 9999
          ENDIF
#elif __DIM == __SFIELD
          IF (lda .NE. 1) THEN
              info = ppm_error_error
              CALL ppm_error(ppm_err_argument,'ppm_map_field_push_2d',  &
     &            'lda must be =1 for scalar data',__LINE__,info)
              GOTO 9999
          ENDIF
          IF (SIZE(fdata,3) .LT. ppm_nsublist(fromtopo)) THEN
              info = ppm_error_error
              CALL ppm_error(ppm_err_argument,'ppm_map_field_push_2d',  &
     &            'field data for some subs is missing',__LINE__,info)
              GOTO 9999
          ENDIF
          IF (SIZE(fdata,1) .LT. xhi) THEN
              info = ppm_error_error
              CALL ppm_error(ppm_err_argument,'ppm_map_field_push_2d',  &
     &            'x dimension of field data does not match mesh',__LINE__,info)
              GOTO 9999
          ENDIF
          IF (SIZE(fdata,2) .LT. yhi) THEN
              info = ppm_error_error
              CALL ppm_error(ppm_err_argument,'ppm_map_field_push_2d',  &
     &            'y dimension of field data does not match mesh',__LINE__,info)
              GOTO 9999
          ENDIF
#endif
          IF (PRESENT(mask)) THEN
              IF (SIZE(mask,1) .LT. xhi) THEN
                  info = ppm_error_error
                  CALL ppm_error(ppm_err_argument,'ppm_map_field_push_2d',  &
     &                'x dimension of mask does not match mesh',__LINE__,info)
                  GOTO 9999
              ENDIF
              IF (SIZE(mask,2) .LT. yhi) THEN
                  info = ppm_error_error
                  CALL ppm_error(ppm_err_argument,'ppm_map_field_push_2d',  &
     &                'y dimension of mask does not match mesh',__LINE__,info)
                  GOTO 9999
              ENDIF
              IF (SIZE(mask,3) .LT. ppm_nsublist(fromtopo)) THEN
                  info = ppm_error_error
                  CALL ppm_error(ppm_err_argument,'ppm_map_field_push_2d',  &
     &                'mask data for some subs is missing',__LINE__,info)
                  GOTO 9999
              ENDIF
          ENDIF
      ENDIF

      !-------------------------------------------------------------------------
      !  Count number of data points to be sent
      !-------------------------------------------------------------------------
      Ndata = 0
      DO i=1,ppm_nsendlist
         !----------------------------------------------------------------------
         !  access mesh blocks belonging to the i-th processor in the 
         !  sendlist
         !----------------------------------------------------------------------
         DO j=ppm_psendbuffer(i),ppm_psendbuffer(i+1)-1
            !-------------------------------------------------------------------
            !  Get the number of mesh points in this block
            !-------------------------------------------------------------------
            Ndata = Ndata + (ppm_mesh_isendblksize(1,j)*    &
     &                       ppm_mesh_isendblksize(2,j))
         ENDDO
      ENDDO

      !-------------------------------------------------------------------------
      !  If there is nothing to be sent we are done
      !-------------------------------------------------------------------------
      IF (Ndata .EQ. 0) THEN
          IF (ppm_debug .GT. 0) THEN
              CALL ppm_write(ppm_rank,'ppm_map_field_push_2d',   &
     &            'There is no data to be sent',info)
          ENDIF
          GOTO 9999
      ENDIF

      !-------------------------------------------------------------------------
      !  Increment the buffer set 
      !-------------------------------------------------------------------------
      ppm_buffer_set = ppm_buffer_set + 1

      !-------------------------------------------------------------------------
      !  Allocate memory for the buffer dimension and type
      !-------------------------------------------------------------------------
      iopt = ppm_param_alloc_grow_preserve
      ldu(1)  = ppm_buffer_set
      CALL ppm_alloc(ppm_buffer_dim,ldu,iopt,info)
      IF (info .NE. 0) THEN
          info = ppm_error_fatal
          CALL ppm_error(ppm_err_alloc,'ppm_map_field_push_2d',     &
     &        'buffer dimensions PPM_BUFFER_DIM',__LINE__,info)
          GOTO 9999
      ENDIF
      CALL ppm_alloc(ppm_buffer_type,ldu,iopt,info)
      IF (info .NE. 0) THEN
          info = ppm_error_fatal
          CALL ppm_error(ppm_err_alloc,'ppm_map_field_push_2d',     &
     &        'buffer types PPM_BUFFER_TYPE',__LINE__,info)
          GOTO 9999
      ENDIF

      !-------------------------------------------------------------------------
      !  A complex number is treated as two reals. Cannot change lda
      !  because it is INTENT(IN)
      !-------------------------------------------------------------------------
#if   __KIND == __SINGLE_PRECISION_COMPLEX | __KIND == __DOUBLE_PRECISION_COMPLEX
      ldb = 2*lda
#else
      ldb = lda
#endif 

      !-------------------------------------------------------------------------
      !  Store the dimension and type
      !-------------------------------------------------------------------------
      ppm_buffer_dim(ppm_buffer_set)  = ldb
#if    __KIND == __SINGLE_PRECISION | __KIND == __SINGLE_PRECISION_COMPLEX
      ppm_buffer_type(ppm_buffer_set) = ppm_kind_single
#elif __KIND == __DOUBLE_PRECISION | __KIND == __DOUBLE_PRECISION_COMPLEX
      ppm_buffer_type(ppm_buffer_set) = ppm_kind_double
#elif __KIND == __INTEGER
      ppm_buffer_type(ppm_buffer_set) = ppm_integer
#elif __KIND == __LOGICAL
      ppm_buffer_type(ppm_buffer_set) = ppm_logical
#endif

      !-------------------------------------------------------------------------
      !  Build the inverse sub list to find local sub indeices based on
      !  global ones (the global ones are communicated)
      !-------------------------------------------------------------------------
      iopt   = ppm_param_alloc_fit
      ldu(1) = ppm_nsublist(ppm_field_topoid)
      CALL ppm_alloc(sublist,ldu,iopt,info)
      IF (info .NE. 0) THEN
          info = ppm_error_fatal
          CALL ppm_error(ppm_err_alloc,'ppm_map_field_push_2d',  &
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
          CALL ppm_error(ppm_err_dealloc,'ppm_map_field_push_2d',     &
     &        'temporary sub list SUBLIST',__LINE__,info)
      ENDIF

      !-------------------------------------------------------------------------
      !  loop over the processors in the ppm_isendlist() 
      !-------------------------------------------------------------------------
      ibuffer = ppm_nsendbuffer
      !-------------------------------------------------------------------------
      !  DOUBLE PRECISION BUFFER
      !-------------------------------------------------------------------------
      IF (ppm_kind.EQ.ppm_kind_double) THEN
         !----------------------------------------------------------------------
         !  (Re)allocate memory for the buffer 
         !----------------------------------------------------------------------
         iopt   = ppm_param_alloc_grow_preserve
         ldu(1) = ppm_nsendbuffer + ldb*Ndata
         CALL ppm_alloc(ppm_sendbufferd,ldu,iopt,info)
         IF (info .NE. 0) THEN
             info = ppm_error_fatal
             CALL ppm_error(ppm_err_alloc,'ppm_map_field_push_2d',     &
     &           'global send buffer PPM_SENDBUFFERD',__LINE__,info)
             GOTO 9999
         ENDIF

         DO i=1,ppm_nsendlist
            !-------------------------------------------------------------------
            !  access mesh blocks belonging to the i-th processor in the 
            !  sendlist
            !-------------------------------------------------------------------
            DO j=ppm_psendbuffer(i),ppm_psendbuffer(i+1)-1
               !----------------------------------------------------------------
               !  Get the sub ID for this mesh block
               !----------------------------------------------------------------
               jsub = ppm_mesh_isendfromsub(j)
               !----------------------------------------------------------------
               !  Translate to local sub ID for storing the data
               !----------------------------------------------------------------
               isub = invsublist(jsub)
               !----------------------------------------------------------------
               !  Mesh offset for this sub
               !----------------------------------------------------------------
               mofs(1) = ppm_cart_mesh(frommesh,fromtopo)%istart(1,jsub)-1
               mofs(2) = ppm_cart_mesh(frommesh,fromtopo)%istart(2,jsub)-1
               !----------------------------------------------------------------
               !  Get boundaries of mesh block to be sent on local sub
               !  coordinates
               !----------------------------------------------------------------
               xlo = ppm_mesh_isendblkstart(1,j)-mofs(1)
               ylo = ppm_mesh_isendblkstart(2,j)-mofs(2)
               xhi = xlo+ppm_mesh_isendblksize(1,j)-1
               yhi = ylo+ppm_mesh_isendblksize(2,j)-1
               IF (ppm_debug .GT. 1) THEN
                   WRITE(mesg,'(A,2I4)') 'start: ',             &
     &                 ppm_mesh_isendblkstart(1,j),ppm_mesh_isendblkstart(2,j)
                   CALL ppm_write(ppm_rank,'ppm_map_field_push_2d',mesg,info)
                   WRITE(mesg,'(A,2I4)') 'size: ',             &
     &                 ppm_mesh_isendblksize(1,j),ppm_mesh_isendblksize(2,j)
                   CALL ppm_write(ppm_rank,'ppm_map_field_push_2d',mesg,info)
                   WRITE(mesg,'(A,2I4)') 'mesh offset: ',mofs(1),mofs(2)
                   CALL ppm_write(ppm_rank,'ppm_map_field_push_2d',mesg,info)
                   WRITE(mesg,'(A,2I4)') 'xlo, xhi: ',xlo,xhi
                   CALL ppm_write(ppm_rank,'ppm_map_field_push_2d',mesg,info)
                   WRITE(mesg,'(A,2I4)') 'ylo, yhi: ',ylo,yhi
                   CALL ppm_write(ppm_rank,'ppm_map_field_push_2d',mesg,info)
                   WRITE(mesg,'(A,I1)') 'buffer dim: ',lda
                   CALL ppm_write(ppm_rank,'ppm_map_field_push_2d',mesg,info)
               ENDIF
               !----------------------------------------------------------------
               !  Loop over all mesh points of this block and append data
               !  to send buffer.
               !----------------------------------------------------------------
               !      ldo = .TRUE.
               !      IF (PRESENT(mask)) THEN
               !          IF (.NOT.mask(1,imesh,jmesh,isub)) ldo = .FALSE.
               !      ENDIF
               !      IF (ldo) THEN
#if    __DIM == __VFIELD
               !----------------------------------------------------------------
               !  Unrolled for lda=1
               !----------------------------------------------------------------
               IF (lda .EQ. 1) THEN
                  DO jmesh=ylo,yhi
                     DO imesh=xlo,xhi
                        ibuffer = ibuffer + 1
#if    __KIND == __SINGLE_PRECISION
                        ppm_sendbufferd(ibuffer) =    &
     &                      REAL(fdata(1,imesh,jmesh,isub),ppm_kind_double)
#elif  __KIND == __DOUBLE_PRECISION
                        ppm_sendbufferd(ibuffer) =    &
     &                      fdata(1,imesh,jmesh,isub)
#elif  __KIND == __SINGLE_PRECISION_COMPLEX
                        ppm_sendbufferd(ibuffer) =    &
     &                      REAL(fdata(1,imesh,jmesh,isub),ppm_kind_double)
                        ibuffer = ibuffer + 1
                        ppm_sendbufferd(ibuffer) =    &
     &                      REAL(AIMAG(fdata(1,imesh,jmesh,isub)),   &
     &                      ppm_kind_double)
#elif  __KIND == __DOUBLE_PRECISION_COMPLEX
                        ppm_sendbufferd(ibuffer) =    &
     &                      REAL(fdata(1,imesh,jmesh,isub),ppm_kind_double)
                        ibuffer = ibuffer + 1
                        ppm_sendbufferd(ibuffer) =    &
     &                      AIMAG(fdata(1,imesh,jmesh,isub))
#elif  __KIND == __INTEGER
                        ppm_sendbufferd(ibuffer) =    &
     &                      REAL(fdata(1,imesh,jmesh,isub),ppm_kind_double)
#elif  __KIND == __LOGICAL
                        IF (fdata(1,imesh,jmesh,isub)) THEN
                           ppm_sendbufferd(ibuffer) = 1.0_ppm_kind_double
                        ELSE
                           ppm_sendbufferd(ibuffer) = 0.0_ppm_kind_double
                        ENDIF
#endif
                     ENDDO
                  ENDDO
               !----------------------------------------------------------------
               !  Unrolled for lda=2
               !----------------------------------------------------------------
               ELSEIF (lda .EQ. 2) THEN
                  DO jmesh=ylo,yhi
                     DO imesh=xlo,xhi
                        ibuffer = ibuffer + 1
#if    __KIND == __SINGLE_PRECISION
                        ppm_sendbufferd(ibuffer) =    &
     &                      REAL(fdata(1,imesh,jmesh,isub),ppm_kind_double)
                        ibuffer = ibuffer + 1
                        ppm_sendbufferd(ibuffer) =    &
     &                      REAL(fdata(2,imesh,jmesh,isub),ppm_kind_double)
#elif  __KIND == __DOUBLE_PRECISION
                        ppm_sendbufferd(ibuffer) =    &
     &                      fdata(1,imesh,jmesh,isub)
                        ibuffer = ibuffer + 1
                        ppm_sendbufferd(ibuffer) =    &
     &                      fdata(2,imesh,jmesh,isub)
#elif  __KIND == __SINGLE_PRECISION_COMPLEX
                        ppm_sendbufferd(ibuffer) =    &
     &                      REAL(fdata(1,imesh,jmesh,isub),ppm_kind_double)
                        ibuffer = ibuffer + 1
                        ppm_sendbufferd(ibuffer) =    &
     &                      REAL(AIMAG(fdata(1,imesh,jmesh,isub)),   &
     &                      ppm_kind_double)
                        ibuffer = ibuffer + 1
                        ppm_sendbufferd(ibuffer) =    &
     &                      REAL(fdata(2,imesh,jmesh,isub),ppm_kind_double)
                        ibuffer = ibuffer + 1
                        ppm_sendbufferd(ibuffer) =    &
     &                      REAL(AIMAG(fdata(2,imesh,jmesh,isub)),   &
     &                      ppm_kind_double)
#elif  __KIND == __DOUBLE_PRECISION_COMPLEX
                        ppm_sendbufferd(ibuffer) =    &
     &                      REAL(fdata(1,imesh,jmesh,isub),ppm_kind_double)
                        ibuffer = ibuffer + 1
                        ppm_sendbufferd(ibuffer) =    &
     &                      AIMAG(fdata(1,imesh,jmesh,isub))
                        ibuffer = ibuffer + 1
                        ppm_sendbufferd(ibuffer) =    &
     &                      REAL(fdata(2,imesh,jmesh,isub),ppm_kind_double)
                        ibuffer = ibuffer + 1
                        ppm_sendbufferd(ibuffer) =    &
     &                      AIMAG(fdata(2,imesh,jmesh,isub))
#elif  __KIND == __INTEGER
                        ppm_sendbufferd(ibuffer) =    &
     &                      REAL(fdata(1,imesh,jmesh,isub),ppm_kind_double)
                        ibuffer = ibuffer + 1
                        ppm_sendbufferd(ibuffer) =    &
     &                      REAL(fdata(2,imesh,jmesh,isub),ppm_kind_double)
#elif  __KIND == __LOGICAL
                        IF (fdata(1,imesh,jmesh,isub)) THEN
                           ppm_sendbufferd(ibuffer) = 1.0_ppm_kind_double
                        ELSE
                           ppm_sendbufferd(ibuffer) = 0.0_ppm_kind_double
                        ENDIF
                        ibuffer = ibuffer + 1
                        IF (fdata(2,imesh,jmesh,isub)) THEN
                           ppm_sendbufferd(ibuffer) = 1.0_ppm_kind_double
                        ELSE
                           ppm_sendbufferd(ibuffer) = 0.0_ppm_kind_double
                        ENDIF
#endif
                     ENDDO
                  ENDDO
               !----------------------------------------------------------------
               !  Unrolled for lda=3
               !----------------------------------------------------------------
               ELSEIF (lda .EQ. 3) THEN
                  DO jmesh=ylo,yhi
                     DO imesh=xlo,xhi
                        ibuffer = ibuffer + 1
#if    __KIND == __SINGLE_PRECISION
                        ppm_sendbufferd(ibuffer) =    &
     &                      REAL(fdata(1,imesh,jmesh,isub),ppm_kind_double)
                        ibuffer = ibuffer + 1
                        ppm_sendbufferd(ibuffer) =    &
     &                      REAL(fdata(2,imesh,jmesh,isub),ppm_kind_double)
                        ibuffer = ibuffer + 1
                        ppm_sendbufferd(ibuffer) =    &
     &                      REAL(fdata(3,imesh,jmesh,isub),ppm_kind_double)
#elif  __KIND == __DOUBLE_PRECISION
                        ppm_sendbufferd(ibuffer) =    &
     &                      fdata(1,imesh,jmesh,isub)
                        ibuffer = ibuffer + 1
                        ppm_sendbufferd(ibuffer) =    &
     &                      fdata(2,imesh,jmesh,isub)
                        ibuffer = ibuffer + 1
                        ppm_sendbufferd(ibuffer) =    &
     &                      fdata(3,imesh,jmesh,isub)
#elif  __KIND == __SINGLE_PRECISION_COMPLEX
                        ppm_sendbufferd(ibuffer) =    &
     &                      REAL(fdata(1,imesh,jmesh,isub),ppm_kind_double)
                        ibuffer = ibuffer + 1
                        ppm_sendbufferd(ibuffer) =    &
     &                      REAL(AIMAG(fdata(1,imesh,jmesh,isub)),   &
     &                      ppm_kind_double)
                        ibuffer = ibuffer + 1
                        ppm_sendbufferd(ibuffer) =    &
     &                      REAL(fdata(2,imesh,jmesh,isub),ppm_kind_double)
                        ibuffer = ibuffer + 1
                        ppm_sendbufferd(ibuffer) =    &
     &                      REAL(AIMAG(fdata(2,imesh,jmesh,isub)),   &
     &                      ppm_kind_double)
                        ibuffer = ibuffer + 1
                        ppm_sendbufferd(ibuffer) =    &
     &                      REAL(fdata(3,imesh,jmesh,isub),ppm_kind_double)
                        ibuffer = ibuffer + 1
                        ppm_sendbufferd(ibuffer) =    &
     &                      REAL(AIMAG(fdata(3,imesh,jmesh,isub)),   &
     &                      ppm_kind_double)
#elif  __KIND == __DOUBLE_PRECISION_COMPLEX
                        ppm_sendbufferd(ibuffer) =    &
     &                      REAL(fdata(1,imesh,jmesh,isub),ppm_kind_double)
                        ibuffer = ibuffer + 1
                        ppm_sendbufferd(ibuffer) =    &
     &                      AIMAG(fdata(1,imesh,jmesh,isub))
                        ibuffer = ibuffer + 1
                        ppm_sendbufferd(ibuffer) =    &
     &                      REAL(fdata(2,imesh,jmesh,isub),ppm_kind_double)
                        ibuffer = ibuffer + 1
                        ppm_sendbufferd(ibuffer) =    &
     &                      AIMAG(fdata(2,imesh,jmesh,isub))
                        ibuffer = ibuffer + 1
                        ppm_sendbufferd(ibuffer) =    &
     &                      REAL(fdata(3,imesh,jmesh,isub),ppm_kind_double)
                        ibuffer = ibuffer + 1
                        ppm_sendbufferd(ibuffer) =    &
     &                      AIMAG(fdata(3,imesh,jmesh,isub))
#elif  __KIND == __INTEGER
                        ppm_sendbufferd(ibuffer) =    &
     &                      REAL(fdata(1,imesh,jmesh,isub),ppm_kind_double)
                        ibuffer = ibuffer + 1
                        ppm_sendbufferd(ibuffer) =    &
     &                      REAL(fdata(2,imesh,jmesh,isub),ppm_kind_double)
                        ibuffer = ibuffer + 1
                        ppm_sendbufferd(ibuffer) =    &
     &                      REAL(fdata(3,imesh,jmesh,isub),ppm_kind_double)
#elif  __KIND == __LOGICAL
                        IF (fdata(1,imesh,jmesh,isub)) THEN
                           ppm_sendbufferd(ibuffer) = 1.0_ppm_kind_double
                        ELSE
                           ppm_sendbufferd(ibuffer) = 0.0_ppm_kind_double
                        ENDIF
                        ibuffer = ibuffer + 1
                        IF (fdata(2,imesh,jmesh,isub)) THEN
                           ppm_sendbufferd(ibuffer) = 1.0_ppm_kind_double
                        ELSE
                           ppm_sendbufferd(ibuffer) = 0.0_ppm_kind_double
                        ENDIF
                        ibuffer = ibuffer + 1
                        IF (fdata(3,imesh,jmesh,isub)) THEN
                           ppm_sendbufferd(ibuffer) = 1.0_ppm_kind_double
                        ELSE
                           ppm_sendbufferd(ibuffer) = 0.0_ppm_kind_double
                        ENDIF
#endif
                     ENDDO
                  ENDDO
               !----------------------------------------------------------------
               !  Unrolled for lda=4
               !----------------------------------------------------------------
               ELSEIF (lda .EQ. 4) THEN
                  DO jmesh=ylo,yhi
                     DO imesh=xlo,xhi
                        ibuffer = ibuffer + 1
#if    __KIND == __SINGLE_PRECISION
                        ppm_sendbufferd(ibuffer) =    &
     &                      REAL(fdata(1,imesh,jmesh,isub),ppm_kind_double)
                        ibuffer = ibuffer + 1
                        ppm_sendbufferd(ibuffer) =    &
     &                      REAL(fdata(2,imesh,jmesh,isub),ppm_kind_double)
                        ibuffer = ibuffer + 1
                        ppm_sendbufferd(ibuffer) =    &
     &                      REAL(fdata(3,imesh,jmesh,isub),ppm_kind_double)
                        ibuffer = ibuffer + 1
                        ppm_sendbufferd(ibuffer) =    &
     &                      REAL(fdata(4,imesh,jmesh,isub),ppm_kind_double)
#elif  __KIND == __DOUBLE_PRECISION
                        ppm_sendbufferd(ibuffer) =    &
     &                      fdata(1,imesh,jmesh,isub)
                        ibuffer = ibuffer + 1
                        ppm_sendbufferd(ibuffer) =    &
     &                      fdata(2,imesh,jmesh,isub)
                        ibuffer = ibuffer + 1
                        ppm_sendbufferd(ibuffer) =    &
     &                      fdata(3,imesh,jmesh,isub)
                        ibuffer = ibuffer + 1
                        ppm_sendbufferd(ibuffer) =    &
     &                      fdata(4,imesh,jmesh,isub)
#elif  __KIND == __SINGLE_PRECISION_COMPLEX
                        ppm_sendbufferd(ibuffer) =    &
     &                      REAL(fdata(1,imesh,jmesh,isub),ppm_kind_double)
                        ibuffer = ibuffer + 1
                        ppm_sendbufferd(ibuffer) =    &
     &                      REAL(AIMAG(fdata(1,imesh,jmesh,isub)),   &
     &                      ppm_kind_double)
                        ibuffer = ibuffer + 1
                        ppm_sendbufferd(ibuffer) =    &
     &                      REAL(fdata(2,imesh,jmesh,isub),ppm_kind_double)
                        ibuffer = ibuffer + 1
                        ppm_sendbufferd(ibuffer) =    &
     &                      REAL(AIMAG(fdata(2,imesh,jmesh,isub)),   &
     &                      ppm_kind_double)
                        ibuffer = ibuffer + 1
                        ppm_sendbufferd(ibuffer) =    &
     &                      REAL(fdata(3,imesh,jmesh,isub),ppm_kind_double)
                        ibuffer = ibuffer + 1
                        ppm_sendbufferd(ibuffer) =    &
     &                      REAL(AIMAG(fdata(3,imesh,jmesh,isub)),   &
     &                      ppm_kind_double)
                        ibuffer = ibuffer + 1
                        ppm_sendbufferd(ibuffer) =    &
     &                      REAL(fdata(4,imesh,jmesh,isub),ppm_kind_double)
                        ibuffer = ibuffer + 1
                        ppm_sendbufferd(ibuffer) =    &
     &                      REAL(AIMAG(fdata(4,imesh,jmesh,isub)),   &
     &                      ppm_kind_double)
#elif  __KIND == __DOUBLE_PRECISION_COMPLEX
                        ppm_sendbufferd(ibuffer) =    &
     &                      REAL(fdata(1,imesh,jmesh,isub),ppm_kind_double)
                        ibuffer = ibuffer + 1
                        ppm_sendbufferd(ibuffer) =    &
     &                      AIMAG(fdata(1,imesh,jmesh,isub))
                        ibuffer = ibuffer + 1
                        ppm_sendbufferd(ibuffer) =    &
     &                      REAL(fdata(2,imesh,jmesh,isub),ppm_kind_double)
                        ibuffer = ibuffer + 1
                        ppm_sendbufferd(ibuffer) =    &
     &                      AIMAG(fdata(2,imesh,jmesh,isub))
                        ibuffer = ibuffer + 1
                        ppm_sendbufferd(ibuffer) =    &
     &                      REAL(fdata(3,imesh,jmesh,isub),ppm_kind_double)
                        ibuffer = ibuffer + 1
                        ppm_sendbufferd(ibuffer) =    &
     &                      AIMAG(fdata(3,imesh,jmesh,isub))
                        ibuffer = ibuffer + 1
                        ppm_sendbufferd(ibuffer) =    &
     &                      REAL(fdata(4,imesh,jmesh,isub),ppm_kind_double)
                        ibuffer = ibuffer + 1
                        ppm_sendbufferd(ibuffer) =    &
     &                      AIMAG(fdata(4,imesh,jmesh,isub))
#elif  __KIND == __INTEGER
                        ppm_sendbufferd(ibuffer) =    &
     &                      REAL(fdata(1,imesh,jmesh,isub),ppm_kind_double)
                        ibuffer = ibuffer + 1
                        ppm_sendbufferd(ibuffer) =    &
     &                      REAL(fdata(2,imesh,jmesh,isub),ppm_kind_double)
                        ibuffer = ibuffer + 1
                        ppm_sendbufferd(ibuffer) =    &
     &                      REAL(fdata(3,imesh,jmesh,isub),ppm_kind_double)
                        ibuffer = ibuffer + 1
                        ppm_sendbufferd(ibuffer) =    &
     &                      REAL(fdata(4,imesh,jmesh,isub),ppm_kind_double)
#elif  __KIND == __LOGICAL
                        IF (fdata(1,imesh,jmesh,isub)) THEN
                           ppm_sendbufferd(ibuffer) = 1.0_ppm_kind_double
                        ELSE
                           ppm_sendbufferd(ibuffer) = 0.0_ppm_kind_double
                        ENDIF
                        ibuffer = ibuffer + 1
                        IF (fdata(2,imesh,jmesh,isub)) THEN
                           ppm_sendbufferd(ibuffer) = 1.0_ppm_kind_double
                        ELSE
                           ppm_sendbufferd(ibuffer) = 0.0_ppm_kind_double
                        ENDIF
                        ibuffer = ibuffer + 1
                        IF (fdata(3,imesh,jmesh,isub)) THEN
                           ppm_sendbufferd(ibuffer) = 1.0_ppm_kind_double
                        ELSE
                           ppm_sendbufferd(ibuffer) = 0.0_ppm_kind_double
                        ENDIF
                        ibuffer = ibuffer + 1
                        IF (fdata(4,imesh,jmesh,isub)) THEN
                           ppm_sendbufferd(ibuffer) = 1.0_ppm_kind_double
                        ELSE
                           ppm_sendbufferd(ibuffer) = 0.0_ppm_kind_double
                        ENDIF
#endif
                     ENDDO
                  ENDDO
               !----------------------------------------------------------------
               !  Unrolled for lda=5
               !----------------------------------------------------------------
               ELSEIF (lda .EQ. 5) THEN
                  DO jmesh=ylo,yhi
                     DO imesh=xlo,xhi
                        ibuffer = ibuffer + 1
#if    __KIND == __SINGLE_PRECISION
                        ppm_sendbufferd(ibuffer) =    &
     &                      REAL(fdata(1,imesh,jmesh,isub),ppm_kind_double)
                        ibuffer = ibuffer + 1
                        ppm_sendbufferd(ibuffer) =    &
     &                      REAL(fdata(2,imesh,jmesh,isub),ppm_kind_double)
                        ibuffer = ibuffer + 1
                        ppm_sendbufferd(ibuffer) =    &
     &                      REAL(fdata(3,imesh,jmesh,isub),ppm_kind_double)
                        ibuffer = ibuffer + 1
                        ppm_sendbufferd(ibuffer) =    &
     &                      REAL(fdata(4,imesh,jmesh,isub),ppm_kind_double)
                        ibuffer = ibuffer + 1
                        ppm_sendbufferd(ibuffer) =    &
     &                      REAL(fdata(5,imesh,jmesh,isub),ppm_kind_double)
#elif  __KIND == __DOUBLE_PRECISION
                        ppm_sendbufferd(ibuffer) =    &
     &                      fdata(1,imesh,jmesh,isub)
                        ibuffer = ibuffer + 1
                        ppm_sendbufferd(ibuffer) =    &
     &                      fdata(2,imesh,jmesh,isub)
                        ibuffer = ibuffer + 1
                        ppm_sendbufferd(ibuffer) =    &
     &                      fdata(3,imesh,jmesh,isub)
                        ibuffer = ibuffer + 1
                        ppm_sendbufferd(ibuffer) =    &
     &                      fdata(4,imesh,jmesh,isub)
                        ibuffer = ibuffer + 1
                        ppm_sendbufferd(ibuffer) =    &
     &                      fdata(5,imesh,jmesh,isub)
#elif  __KIND == __SINGLE_PRECISION_COMPLEX
                        ppm_sendbufferd(ibuffer) =    &
     &                      REAL(fdata(1,imesh,jmesh,isub),ppm_kind_double)
                        ibuffer = ibuffer + 1
                        ppm_sendbufferd(ibuffer) =    &
     &                      REAL(AIMAG(fdata(1,imesh,jmesh,isub)),   &
     &                      ppm_kind_double)
                        ibuffer = ibuffer + 1
                        ppm_sendbufferd(ibuffer) =    &
     &                      REAL(fdata(2,imesh,jmesh,isub),ppm_kind_double)
                        ibuffer = ibuffer + 1
                        ppm_sendbufferd(ibuffer) =    &
     &                      REAL(AIMAG(fdata(2,imesh,jmesh,isub)),   &
     &                      ppm_kind_double)
                        ibuffer = ibuffer + 1
                        ppm_sendbufferd(ibuffer) =    &
     &                      REAL(fdata(3,imesh,jmesh,isub),ppm_kind_double)
                        ibuffer = ibuffer + 1
                        ppm_sendbufferd(ibuffer) =    &
     &                      REAL(AIMAG(fdata(3,imesh,jmesh,isub)),   &
     &                      ppm_kind_double)
                        ibuffer = ibuffer + 1
                        ppm_sendbufferd(ibuffer) =    &
     &                      REAL(fdata(4,imesh,jmesh,isub),ppm_kind_double)
                        ibuffer = ibuffer + 1
                        ppm_sendbufferd(ibuffer) =    &
     &                      REAL(AIMAG(fdata(4,imesh,jmesh,isub)),   &
     &                      ppm_kind_double)
                        ibuffer = ibuffer + 1
                        ppm_sendbufferd(ibuffer) =    &
     &                      REAL(fdata(5,imesh,jmesh,isub),ppm_kind_double)
                        ibuffer = ibuffer + 1
                        ppm_sendbufferd(ibuffer) =    &
     &                      REAL(AIMAG(fdata(5,imesh,jmesh,isub)),   &
     &                      ppm_kind_double)
#elif  __KIND == __DOUBLE_PRECISION_COMPLEX
                        ppm_sendbufferd(ibuffer) =    &
     &                      REAL(fdata(1,imesh,jmesh,isub),ppm_kind_double)
                        ibuffer = ibuffer + 1
                        ppm_sendbufferd(ibuffer) =    &
     &                      AIMAG(fdata(1,imesh,jmesh,isub))
                        ibuffer = ibuffer + 1
                        ppm_sendbufferd(ibuffer) =    &
     &                      REAL(fdata(2,imesh,jmesh,isub),ppm_kind_double)
                        ibuffer = ibuffer + 1
                        ppm_sendbufferd(ibuffer) =    &
     &                      AIMAG(fdata(2,imesh,jmesh,isub))
                        ibuffer = ibuffer + 1
                        ppm_sendbufferd(ibuffer) =    &
     &                      REAL(fdata(3,imesh,jmesh,isub),ppm_kind_double)
                        ibuffer = ibuffer + 1
                        ppm_sendbufferd(ibuffer) =    &
     &                      AIMAG(fdata(3,imesh,jmesh,isub))
                        ibuffer = ibuffer + 1
                        ppm_sendbufferd(ibuffer) =    &
     &                      REAL(fdata(4,imesh,jmesh,isub),ppm_kind_double)
                        ibuffer = ibuffer + 1
                        ppm_sendbufferd(ibuffer) =    &
     &                      AIMAG(fdata(4,imesh,jmesh,isub))
                        ibuffer = ibuffer + 1
                        ppm_sendbufferd(ibuffer) =    &
     &                      REAL(fdata(5,imesh,jmesh,isub),ppm_kind_double)
                        ibuffer = ibuffer + 1
                        ppm_sendbufferd(ibuffer) =    &
     &                      AIMAG(fdata(5,imesh,jmesh,isub))
#elif  __KIND == __INTEGER
                        ppm_sendbufferd(ibuffer) =    &
     &                      REAL(fdata(1,imesh,jmesh,isub),ppm_kind_double)
                        ibuffer = ibuffer + 1
                        ppm_sendbufferd(ibuffer) =    &
     &                      REAL(fdata(2,imesh,jmesh,isub),ppm_kind_double)
                        ibuffer = ibuffer + 1
                        ppm_sendbufferd(ibuffer) =    &
     &                      REAL(fdata(3,imesh,jmesh,isub),ppm_kind_double)
                        ibuffer = ibuffer + 1
                        ppm_sendbufferd(ibuffer) =    &
     &                      REAL(fdata(4,imesh,jmesh,isub),ppm_kind_double)
                        ibuffer = ibuffer + 1
                        ppm_sendbufferd(ibuffer) =    &
     &                      REAL(fdata(5,imesh,jmesh,isub),ppm_kind_double)
#elif  __KIND == __LOGICAL
                        IF (fdata(1,imesh,jmesh,isub)) THEN
                           ppm_sendbufferd(ibuffer) = 1.0_ppm_kind_double
                        ELSE
                           ppm_sendbufferd(ibuffer) = 0.0_ppm_kind_double
                        ENDIF
                        ibuffer = ibuffer + 1
                        IF (fdata(2,imesh,jmesh,isub)) THEN
                           ppm_sendbufferd(ibuffer) = 1.0_ppm_kind_double
                        ELSE
                           ppm_sendbufferd(ibuffer) = 0.0_ppm_kind_double
                        ENDIF
                        ibuffer = ibuffer + 1
                        IF (fdata(3,imesh,jmesh,isub)) THEN
                           ppm_sendbufferd(ibuffer) = 1.0_ppm_kind_double
                        ELSE
                           ppm_sendbufferd(ibuffer) = 0.0_ppm_kind_double
                        ENDIF
                        ibuffer = ibuffer + 1
                        IF (fdata(4,imesh,jmesh,isub)) THEN
                           ppm_sendbufferd(ibuffer) = 1.0_ppm_kind_double
                        ELSE
                           ppm_sendbufferd(ibuffer) = 0.0_ppm_kind_double
                        ENDIF
                        ibuffer = ibuffer + 1
                        IF (fdata(5,imesh,jmesh,isub)) THEN
                           ppm_sendbufferd(ibuffer) = 1.0_ppm_kind_double
                        ELSE
                           ppm_sendbufferd(ibuffer) = 0.0_ppm_kind_double
                        ENDIF
#endif
                     ENDDO
                  ENDDO
               !----------------------------------------------------------------
               !  For lda.GT.5 vector length is lda
               !----------------------------------------------------------------
               ELSE
                  DO jmesh=ylo,yhi
                     DO imesh=xlo,xhi
                        DO k=1,lda
                           ibuffer = ibuffer + 1
#if    __KIND == __SINGLE_PRECISION
                           ppm_sendbufferd(ibuffer) =    &
     &                        REAL(fdata(k,imesh,jmesh,isub),ppm_kind_double)
#elif  __KIND == __DOUBLE_PRECISION
                           ppm_sendbufferd(ibuffer) =    &
     &                         fdata(k,imesh,jmesh,isub)
#elif  __KIND == __SINGLE_PRECISION_COMPLEX
                           ppm_sendbufferd(ibuffer) =    &
     &                        REAL(fdata(k,imesh,jmesh,isub),ppm_kind_double)
                           ibuffer = ibuffer + 1
                           ppm_sendbufferd(ibuffer) =    &
     &                         REAL(AIMAG(fdata(k,imesh,jmesh,isub)),   &
     &                         ppm_kind_double)
#elif  __KIND == __DOUBLE_PRECISION_COMPLEX
                           ppm_sendbufferd(ibuffer) =    &
     &                        REAL(fdata(k,imesh,jmesh,isub),ppm_kind_double)
                           ibuffer = ibuffer + 1
                           ppm_sendbufferd(ibuffer) =    &
     &                         AIMAG(fdata(k,imesh,jmesh,isub))
#elif  __KIND == __INTEGER
                           ppm_sendbufferd(ibuffer) =    &
     &                        REAL(fdata(k,imesh,jmesh,isub),ppm_kind_double)
#elif  __KIND == __LOGICAL
                           IF (fdata(k,imesh,jmesh,isub)) THEN
                              ppm_sendbufferd(ibuffer) = 1.0_ppm_kind_double
                           ELSE
                              ppm_sendbufferd(ibuffer) = 0.0_ppm_kind_double
                           ENDIF
#endif
                        ENDDO
                     ENDDO
                  ENDDO
               ENDIF
#elif  __DIM == __SFIELD
               !----------------------------------------------------------------
               !  Scalar version
               !----------------------------------------------------------------
               DO jmesh=ylo,yhi
                  DO imesh=xlo,xhi
                     ibuffer = ibuffer + 1
#if    __KIND == __SINGLE_PRECISION
                     ppm_sendbufferd(ibuffer) =    &
     &                   REAL(fdata(imesh,jmesh,isub),ppm_kind_double)
#elif  __KIND == __DOUBLE_PRECISION
                     ppm_sendbufferd(ibuffer) = fdata(imesh,jmesh,isub)
#elif  __KIND == __SINGLE_PRECISION_COMPLEX
                     ppm_sendbufferd(ibuffer) =    &
     &                   REAL(fdata(imesh,jmesh,isub),ppm_kind_double)
                     ibuffer = ibuffer + 1
                     ppm_sendbufferd(ibuffer) =    &
     &                   REAL(AIMAG(fdata(imesh,jmesh,isub)),   &
     &                   ppm_kind_double)
#elif  __KIND == __DOUBLE_PRECISION_COMPLEX
                     ppm_sendbufferd(ibuffer) =    &
     &                   REAL(fdata(imesh,jmesh,isub),ppm_kind_double)
                     ibuffer = ibuffer + 1
                     ppm_sendbufferd(ibuffer) = AIMAG(fdata(imesh,jmesh,isub))
#elif  __KIND == __INTEGER
                     ppm_sendbufferd(ibuffer) =    &
     &                   REAL(fdata(imesh,jmesh,isub),ppm_kind_double)
#elif  __KIND == __LOGICAL
                     IF (fdata(imesh,jmesh,isub)) THEN
                        ppm_sendbufferd(ibuffer) = 1.0_ppm_kind_double
                     ELSE
                        ppm_sendbufferd(ibuffer) = 0.0_ppm_kind_double
                     ENDIF
#endif
                  ENDDO
               ENDDO
#endif
            ENDDO        ! i=1,ppm_nsendlist
         ENDDO       ! psendbuffer
      !-------------------------------------------------------------------------
      !  SINGLE PRECISION BUFFER
      !-------------------------------------------------------------------------
      ELSE
         !----------------------------------------------------------------------
         !  (Re)allocate memory for the buffer 
         !----------------------------------------------------------------------
         iopt   = ppm_param_alloc_grow_preserve
         ldu(1) = ppm_nsendbuffer + ldb*Ndata
         CALL ppm_alloc(ppm_sendbuffers,ldu,iopt,info)
         IF (info .NE. 0) THEN
             info = ppm_error_fatal
             CALL ppm_error(ppm_err_alloc,'ppm_map_field_push_2d',     &
     &           'global send buffer PPM_SENDBUFFERS',__LINE__,info)
             GOTO 9999
         ENDIF

         DO i=1,ppm_nsendlist
            !-------------------------------------------------------------------
            !  access mesh blocks belonging to the i-th processor in the 
            !  sendlist
            !-------------------------------------------------------------------
            DO j=ppm_psendbuffer(i),ppm_psendbuffer(i+1)-1
               !----------------------------------------------------------------
               !  Get the sub ID for this mesh block
               !----------------------------------------------------------------
               jsub = ppm_mesh_isendfromsub(j)
               !----------------------------------------------------------------
               !  Translate to local sub ID for storing the data
               !----------------------------------------------------------------
               isub = invsublist(jsub)
               !----------------------------------------------------------------
               !  Mesh offset for this sub
               !----------------------------------------------------------------
               mofs(1) = ppm_cart_mesh(frommesh,fromtopo)%istart(1,jsub)-1
               mofs(2) = ppm_cart_mesh(frommesh,fromtopo)%istart(2,jsub)-1
               !----------------------------------------------------------------
               !  Get boundaries of mesh block to be sent on local sub
               !  coordinates
               !----------------------------------------------------------------
               xlo = ppm_mesh_isendblkstart(1,j)-mofs(1)
               ylo = ppm_mesh_isendblkstart(2,j)-mofs(2)
               xhi = xlo+ppm_mesh_isendblksize(1,j)-1
               yhi = ylo+ppm_mesh_isendblksize(2,j)-1
               IF (ppm_debug .GT. 1) THEN
                   WRITE(mesg,'(A,2I4)') 'start: ',             &
     &                 ppm_mesh_isendblkstart(1,j),ppm_mesh_isendblkstart(2,j)
                   CALL ppm_write(ppm_rank,'ppm_map_field_push_2d',mesg,info)
                   WRITE(mesg,'(A,2I4)') 'size: ',             &
     &                 ppm_mesh_isendblksize(1,j),ppm_mesh_isendblksize(2,j)
                   CALL ppm_write(ppm_rank,'ppm_map_field_push_2d',mesg,info)
                   WRITE(mesg,'(A,2I4)') 'mesh offset: ',mofs(1),mofs(2)
                   CALL ppm_write(ppm_rank,'ppm_map_field_push_2d',mesg,info)
                   WRITE(mesg,'(A,2I4)') 'xlo, xhi: ',xlo,xhi
                   CALL ppm_write(ppm_rank,'ppm_map_field_push_2d',mesg,info)
                   WRITE(mesg,'(A,2I4)') 'ylo, yhi: ',ylo,yhi
                   CALL ppm_write(ppm_rank,'ppm_map_field_push_2d',mesg,info)
                   WRITE(mesg,'(A,I1)') 'buffer dim: ',lda
                   CALL ppm_write(ppm_rank,'ppm_map_field_push_2d',mesg,info)
               ENDIF
               !----------------------------------------------------------------
               !  Loop over all mesh points of this block and append data
               !  to send buffer.
               !----------------------------------------------------------------
               !     ldo = .TRUE.
               !     IF (PRESENT(mask)) THEN
               !         IF (.NOT.mask(1,imesh,jmesh,isub)) ldo = .FALSE.
               !     ENDIF
               !     IF (ldo) THEN
#if    __DIM == __VFIELD
               !----------------------------------------------------------------
               !  Unrolled for lda=1
               !----------------------------------------------------------------
               IF (lda .EQ. 1) THEN
                  DO jmesh=ylo,yhi
                     DO imesh=xlo,xhi
                        ibuffer = ibuffer + 1
#if    __KIND == __DOUBLE_PRECISION
                        ppm_sendbuffers(ibuffer) =    &
     &                     REAL(fdata(1,imesh,jmesh,isub),ppm_kind_single)
#elif  __KIND == __SINGLE_PRECISION
                        ppm_sendbuffers(ibuffer) =    &
     &                     fdata(1,imesh,jmesh,isub)
#elif  __KIND == __SINGLE_PRECISION_COMPLEX
                        ppm_sendbuffers(ibuffer) =    &
     &                     REAL(fdata(1,imesh,jmesh,isub),ppm_kind_single)
                        ibuffer = ibuffer + 1
                        ppm_sendbuffers(ibuffer) =    &
     &                     AIMAG(fdata(1,imesh,jmesh,isub))
#elif  __KIND == __DOUBLE_PRECISION_COMPLEX
                        ppm_sendbuffers(ibuffer) =    &
     &                      REAL(fdata(1,imesh,jmesh,isub),ppm_kind_single)
                        ibuffer = ibuffer + 1
                        ppm_sendbuffers(ibuffer) =    &
     &                      REAL(AIMAG(fdata(1,imesh,jmesh,isub)),   &
     &                      ppm_kind_single)
#elif  __KIND == __INTEGER
                        ppm_sendbuffers(ibuffer) =    &
     &                      REAL(fdata(1,imesh,jmesh,isub),ppm_kind_single)
#elif  __KIND == __LOGICAL
                        IF (fdata(1,imesh,jmesh,isub)) THEN
                           ppm_sendbuffers(ibuffer) = 1.0_ppm_kind_single
                        ELSE
                           ppm_sendbuffers(ibuffer) = 0.0_ppm_kind_single
                        ENDIF
#endif
                     ENDDO
                  ENDDO
               !----------------------------------------------------------------
               !  Unrolled for lda=2
               !----------------------------------------------------------------
               ELSEIF (lda .EQ. 2) THEN
                  DO jmesh=ylo,yhi
                     DO imesh=xlo,xhi
                        ibuffer = ibuffer + 1
#if    __KIND == __DOUBLE_PRECISION
                        ppm_sendbuffers(ibuffer) =    &
     &                     REAL(fdata(1,imesh,jmesh,isub),ppm_kind_single)
                        ibuffer = ibuffer + 1
                        ppm_sendbuffers(ibuffer) =    &
     &                     REAL(fdata(2,imesh,jmesh,isub),ppm_kind_single)
#elif  __KIND == __SINGLE_PRECISION
                        ppm_sendbuffers(ibuffer) =    &
     &                     fdata(1,imesh,jmesh,isub)
                        ibuffer = ibuffer + 1
                        ppm_sendbuffers(ibuffer) =    &
     &                     fdata(2,imesh,jmesh,isub)
#elif  __KIND == __SINGLE_PRECISION_COMPLEX
                        ppm_sendbuffers(ibuffer) =    &
     &                     REAL(fdata(1,imesh,jmesh,isub),ppm_kind_single)
                        ibuffer = ibuffer + 1
                        ppm_sendbuffers(ibuffer) =    &
     &                     AIMAG(fdata(1,imesh,jmesh,isub))
                        ibuffer = ibuffer + 1
                        ppm_sendbuffers(ibuffer) =    &
     &                     REAL(fdata(2,imesh,jmesh,isub),ppm_kind_single)
                        ibuffer = ibuffer + 1
                        ppm_sendbuffers(ibuffer) =    &
     &                     AIMAG(fdata(2,imesh,jmesh,isub))
#elif  __KIND == __DOUBLE_PRECISION_COMPLEX
                        ppm_sendbuffers(ibuffer) =    &
     &                      REAL(fdata(1,imesh,jmesh,isub),ppm_kind_single)
                        ibuffer = ibuffer + 1
                        ppm_sendbuffers(ibuffer) =    &
     &                      REAL(AIMAG(fdata(1,imesh,jmesh,isub)),   &
     &                      ppm_kind_single)
                        ibuffer = ibuffer + 1
                        ppm_sendbuffers(ibuffer) =    &
     &                      REAL(fdata(2,imesh,jmesh,isub),ppm_kind_single)
                        ibuffer = ibuffer + 1
                        ppm_sendbuffers(ibuffer) =    &
     &                      REAL(AIMAG(fdata(2,imesh,jmesh,isub)),   &
     &                      ppm_kind_single)
#elif  __KIND == __INTEGER
                        ppm_sendbuffers(ibuffer) =    &
     &                      REAL(fdata(1,imesh,jmesh,isub),ppm_kind_single)
                        ibuffer = ibuffer + 1
                        ppm_sendbuffers(ibuffer) =    &
     &                      REAL(fdata(2,imesh,jmesh,isub),ppm_kind_single)
#elif  __KIND == __LOGICAL
                        IF (fdata(1,imesh,jmesh,isub)) THEN
                           ppm_sendbuffers(ibuffer) = 1.0_ppm_kind_single
                        ELSE
                           ppm_sendbuffers(ibuffer) = 0.0_ppm_kind_single
                        ENDIF
                        ibuffer = ibuffer + 1
                        IF (fdata(2,imesh,jmesh,isub)) THEN
                           ppm_sendbuffers(ibuffer) = 1.0_ppm_kind_single
                        ELSE
                           ppm_sendbuffers(ibuffer) = 0.0_ppm_kind_single
                        ENDIF
#endif
                     ENDDO
                  ENDDO
               !----------------------------------------------------------------
               !  Unrolled for lda=3
               !----------------------------------------------------------------
               ELSEIF (lda .EQ. 3) THEN
                  DO jmesh=ylo,yhi
                     DO imesh=xlo,xhi
                        ibuffer = ibuffer + 1
#if    __KIND == __DOUBLE_PRECISION
                        ppm_sendbuffers(ibuffer) =    &
     &                     REAL(fdata(1,imesh,jmesh,isub),ppm_kind_single)
                        ibuffer = ibuffer + 1
                        ppm_sendbuffers(ibuffer) =    &
     &                     REAL(fdata(2,imesh,jmesh,isub),ppm_kind_single)
                        ibuffer = ibuffer + 1
                        ppm_sendbuffers(ibuffer) =    &
     &                     REAL(fdata(3,imesh,jmesh,isub),ppm_kind_single)
#elif  __KIND == __SINGLE_PRECISION
                        ppm_sendbuffers(ibuffer) =    &
     &                     fdata(1,imesh,jmesh,isub)
                        ibuffer = ibuffer + 1
                        ppm_sendbuffers(ibuffer) =    &
     &                     fdata(2,imesh,jmesh,isub)
                        ibuffer = ibuffer + 1
                        ppm_sendbuffers(ibuffer) =    &
     &                     fdata(3,imesh,jmesh,isub)
#elif  __KIND == __SINGLE_PRECISION_COMPLEX
                        ppm_sendbuffers(ibuffer) =    &
     &                     REAL(fdata(1,imesh,jmesh,isub),ppm_kind_single)
                        ibuffer = ibuffer + 1
                        ppm_sendbuffers(ibuffer) =    &
     &                     AIMAG(fdata(1,imesh,jmesh,isub))
                        ibuffer = ibuffer + 1
                        ppm_sendbuffers(ibuffer) =    &
     &                     REAL(fdata(2,imesh,jmesh,isub),ppm_kind_single)
                        ibuffer = ibuffer + 1
                        ppm_sendbuffers(ibuffer) =    &
     &                     AIMAG(fdata(2,imesh,jmesh,isub))
                        ibuffer = ibuffer + 1
                        ppm_sendbuffers(ibuffer) =    &
     &                     REAL(fdata(3,imesh,jmesh,isub),ppm_kind_single)
                        ibuffer = ibuffer + 1
                        ppm_sendbuffers(ibuffer) =    &
     &                     AIMAG(fdata(3,imesh,jmesh,isub))
#elif  __KIND == __DOUBLE_PRECISION_COMPLEX
                        ppm_sendbuffers(ibuffer) =    &
     &                      REAL(fdata(1,imesh,jmesh,isub),ppm_kind_single)
                        ibuffer = ibuffer + 1
                        ppm_sendbuffers(ibuffer) =    &
     &                      REAL(AIMAG(fdata(1,imesh,jmesh,isub)),   &
     &                      ppm_kind_single)
                        ibuffer = ibuffer + 1
                        ppm_sendbuffers(ibuffer) =    &
     &                      REAL(fdata(2,imesh,jmesh,isub),ppm_kind_single)
                        ibuffer = ibuffer + 1
                        ppm_sendbuffers(ibuffer) =    &
     &                      REAL(AIMAG(fdata(2,imesh,jmesh,isub)),   &
     &                      ppm_kind_single)
                        ibuffer = ibuffer + 1
                        ppm_sendbuffers(ibuffer) =    &
     &                      REAL(fdata(3,imesh,jmesh,isub),ppm_kind_single)
                        ibuffer = ibuffer + 1
                        ppm_sendbuffers(ibuffer) =    &
     &                      REAL(AIMAG(fdata(3,imesh,jmesh,isub)),   &
     &                      ppm_kind_single)
#elif  __KIND == __INTEGER
                        ppm_sendbuffers(ibuffer) =    &
     &                      REAL(fdata(1,imesh,jmesh,isub),ppm_kind_single)
                        ibuffer = ibuffer + 1
                        ppm_sendbuffers(ibuffer) =    &
     &                      REAL(fdata(2,imesh,jmesh,isub),ppm_kind_single)
                        ibuffer = ibuffer + 1
                        ppm_sendbuffers(ibuffer) =    &
     &                      REAL(fdata(3,imesh,jmesh,isub),ppm_kind_single)
#elif  __KIND == __LOGICAL
                        IF (fdata(1,imesh,jmesh,isub)) THEN
                           ppm_sendbuffers(ibuffer) = 1.0_ppm_kind_single
                        ELSE
                           ppm_sendbuffers(ibuffer) = 0.0_ppm_kind_single
                        ENDIF
                        ibuffer = ibuffer + 1
                        IF (fdata(2,imesh,jmesh,isub)) THEN
                           ppm_sendbuffers(ibuffer) = 1.0_ppm_kind_single
                        ELSE
                           ppm_sendbuffers(ibuffer) = 0.0_ppm_kind_single
                        ENDIF
                        ibuffer = ibuffer + 1
                        IF (fdata(3,imesh,jmesh,isub)) THEN
                           ppm_sendbuffers(ibuffer) = 1.0_ppm_kind_single
                        ELSE
                           ppm_sendbuffers(ibuffer) = 0.0_ppm_kind_single
                        ENDIF
#endif
                     ENDDO
                  ENDDO
               !----------------------------------------------------------------
               !  Unrolled for lda=4
               !----------------------------------------------------------------
               ELSEIF (lda .EQ. 4) THEN
                  DO jmesh=ylo,yhi
                     DO imesh=xlo,xhi
                        ibuffer = ibuffer + 1
#if    __KIND == __DOUBLE_PRECISION
                        ppm_sendbuffers(ibuffer) =    &
     &                     REAL(fdata(1,imesh,jmesh,isub),ppm_kind_single)
                        ibuffer = ibuffer + 1
                        ppm_sendbuffers(ibuffer) =    &
     &                     REAL(fdata(2,imesh,jmesh,isub),ppm_kind_single)
                        ibuffer = ibuffer + 1
                        ppm_sendbuffers(ibuffer) =    &
     &                     REAL(fdata(3,imesh,jmesh,isub),ppm_kind_single)
                        ibuffer = ibuffer + 1
                        ppm_sendbuffers(ibuffer) =    &
     &                     REAL(fdata(4,imesh,jmesh,isub),ppm_kind_single)
#elif  __KIND == __SINGLE_PRECISION
                        ppm_sendbuffers(ibuffer) =    &
     &                     fdata(1,imesh,jmesh,isub)
                        ibuffer = ibuffer + 1
                        ppm_sendbuffers(ibuffer) =    &
     &                     fdata(2,imesh,jmesh,isub)
                        ibuffer = ibuffer + 1
                        ppm_sendbuffers(ibuffer) =    &
     &                     fdata(3,imesh,jmesh,isub)
                        ibuffer = ibuffer + 1
                        ppm_sendbuffers(ibuffer) =    &
     &                     fdata(4,imesh,jmesh,isub)
#elif  __KIND == __SINGLE_PRECISION_COMPLEX
                        ppm_sendbuffers(ibuffer) =    &
     &                     REAL(fdata(1,imesh,jmesh,isub),ppm_kind_single)
                        ibuffer = ibuffer + 1
                        ppm_sendbuffers(ibuffer) =    &
     &                     AIMAG(fdata(1,imesh,jmesh,isub))
                        ibuffer = ibuffer + 1
                        ppm_sendbuffers(ibuffer) =    &
     &                     REAL(fdata(2,imesh,jmesh,isub),ppm_kind_single)
                        ibuffer = ibuffer + 1
                        ppm_sendbuffers(ibuffer) =    &
     &                     AIMAG(fdata(2,imesh,jmesh,isub))
                        ibuffer = ibuffer + 1
                        ppm_sendbuffers(ibuffer) =    &
     &                     REAL(fdata(3,imesh,jmesh,isub),ppm_kind_single)
                        ibuffer = ibuffer + 1
                        ppm_sendbuffers(ibuffer) =    &
     &                     AIMAG(fdata(3,imesh,jmesh,isub))
                        ibuffer = ibuffer + 1
                        ppm_sendbuffers(ibuffer) =    &
     &                     REAL(fdata(4,imesh,jmesh,isub),ppm_kind_single)
                        ibuffer = ibuffer + 1
                        ppm_sendbuffers(ibuffer) =    &
     &                     AIMAG(fdata(4,imesh,jmesh,isub))
#elif  __KIND == __DOUBLE_PRECISION_COMPLEX
                        ppm_sendbuffers(ibuffer) =    &
     &                      REAL(fdata(1,imesh,jmesh,isub),ppm_kind_single)
                        ibuffer = ibuffer + 1
                        ppm_sendbuffers(ibuffer) =    &
     &                      REAL(AIMAG(fdata(1,imesh,jmesh,isub)),   &
     &                      ppm_kind_single)
                        ibuffer = ibuffer + 1
                        ppm_sendbuffers(ibuffer) =    &
     &                      REAL(fdata(2,imesh,jmesh,isub),ppm_kind_single)
                        ibuffer = ibuffer + 1
                        ppm_sendbuffers(ibuffer) =    &
     &                      REAL(AIMAG(fdata(2,imesh,jmesh,isub)),   &
     &                      ppm_kind_single)
                        ibuffer = ibuffer + 1
                        ppm_sendbuffers(ibuffer) =    &
     &                      REAL(fdata(3,imesh,jmesh,isub),ppm_kind_single)
                        ibuffer = ibuffer + 1
                        ppm_sendbuffers(ibuffer) =    &
     &                      REAL(AIMAG(fdata(3,imesh,jmesh,isub)),   &
     &                      ppm_kind_single)
                        ibuffer = ibuffer + 1
                        ppm_sendbuffers(ibuffer) =    &
     &                      REAL(fdata(4,imesh,jmesh,isub),ppm_kind_single)
                        ibuffer = ibuffer + 1
                        ppm_sendbuffers(ibuffer) =    &
     &                      REAL(AIMAG(fdata(4,imesh,jmesh,isub)),   &
     &                      ppm_kind_single)
#elif  __KIND == __INTEGER
                        ppm_sendbuffers(ibuffer) =    &
     &                      REAL(fdata(1,imesh,jmesh,isub),ppm_kind_single)
                        ibuffer = ibuffer + 1
                        ppm_sendbuffers(ibuffer) =    &
     &                      REAL(fdata(2,imesh,jmesh,isub),ppm_kind_single)
                        ibuffer = ibuffer + 1
                        ppm_sendbuffers(ibuffer) =    &
     &                      REAL(fdata(3,imesh,jmesh,isub),ppm_kind_single)
                        ibuffer = ibuffer + 1
                        ppm_sendbuffers(ibuffer) =    &
     &                      REAL(fdata(4,imesh,jmesh,isub),ppm_kind_single)
#elif  __KIND == __LOGICAL
                        IF (fdata(1,imesh,jmesh,isub)) THEN
                           ppm_sendbuffers(ibuffer) = 1.0_ppm_kind_single
                        ELSE
                           ppm_sendbuffers(ibuffer) = 0.0_ppm_kind_single
                        ENDIF
                        ibuffer = ibuffer + 1
                        IF (fdata(2,imesh,jmesh,isub)) THEN
                           ppm_sendbuffers(ibuffer) = 1.0_ppm_kind_single
                        ELSE
                           ppm_sendbuffers(ibuffer) = 0.0_ppm_kind_single
                        ENDIF
                        ibuffer = ibuffer + 1
                        IF (fdata(3,imesh,jmesh,isub)) THEN
                           ppm_sendbuffers(ibuffer) = 1.0_ppm_kind_single
                        ELSE
                           ppm_sendbuffers(ibuffer) = 0.0_ppm_kind_single
                        ENDIF
                        ibuffer = ibuffer + 1
                        IF (fdata(4,imesh,jmesh,isub)) THEN
                           ppm_sendbuffers(ibuffer) = 1.0_ppm_kind_single
                        ELSE
                           ppm_sendbuffers(ibuffer) = 0.0_ppm_kind_single
                        ENDIF
#endif
                     ENDDO
                  ENDDO
               !----------------------------------------------------------------
               !  Unrolled for lda=5
               !----------------------------------------------------------------
               ELSEIF (lda .EQ. 5) THEN
                  DO jmesh=ylo,yhi
                     DO imesh=xlo,xhi
                        ibuffer = ibuffer + 1
#if    __KIND == __DOUBLE_PRECISION
                        ppm_sendbuffers(ibuffer) =    &
     &                     REAL(fdata(1,imesh,jmesh,isub),ppm_kind_single)
                        ibuffer = ibuffer + 1
                        ppm_sendbuffers(ibuffer) =    &
     &                     REAL(fdata(2,imesh,jmesh,isub),ppm_kind_single)
                        ibuffer = ibuffer + 1
                        ppm_sendbuffers(ibuffer) =    &
     &                     REAL(fdata(3,imesh,jmesh,isub),ppm_kind_single)
                        ibuffer = ibuffer + 1
                        ppm_sendbuffers(ibuffer) =    &
     &                     REAL(fdata(4,imesh,jmesh,isub),ppm_kind_single)
                        ibuffer = ibuffer + 1
                        ppm_sendbuffers(ibuffer) =    &
     &                     REAL(fdata(5,imesh,jmesh,isub),ppm_kind_single)
#elif  __KIND == __SINGLE_PRECISION
                        ppm_sendbuffers(ibuffer) =    &
     &                     fdata(1,imesh,jmesh,isub)
                        ibuffer = ibuffer + 1
                        ppm_sendbuffers(ibuffer) =    &
     &                     fdata(2,imesh,jmesh,isub)
                        ibuffer = ibuffer + 1
                        ppm_sendbuffers(ibuffer) =    &
     &                     fdata(3,imesh,jmesh,isub)
                        ibuffer = ibuffer + 1
                        ppm_sendbuffers(ibuffer) =    &
     &                     fdata(4,imesh,jmesh,isub)
                        ibuffer = ibuffer + 1
                        ppm_sendbuffers(ibuffer) =    &
     &                     fdata(5,imesh,jmesh,isub)
#elif  __KIND == __SINGLE_PRECISION_COMPLEX
                        ppm_sendbuffers(ibuffer) =    &
     &                     REAL(fdata(1,imesh,jmesh,isub),ppm_kind_single)
                        ibuffer = ibuffer + 1
                        ppm_sendbuffers(ibuffer) =    &
     &                     AIMAG(fdata(1,imesh,jmesh,isub))
                        ibuffer = ibuffer + 1
                        ppm_sendbuffers(ibuffer) =    &
     &                     REAL(fdata(2,imesh,jmesh,isub),ppm_kind_single)
                        ibuffer = ibuffer + 1
                        ppm_sendbuffers(ibuffer) =    &
     &                     AIMAG(fdata(2,imesh,jmesh,isub))
                        ibuffer = ibuffer + 1
                        ppm_sendbuffers(ibuffer) =    &
     &                     REAL(fdata(3,imesh,jmesh,isub),ppm_kind_single)
                        ibuffer = ibuffer + 1
                        ppm_sendbuffers(ibuffer) =    &
     &                     AIMAG(fdata(3,imesh,jmesh,isub))
                        ibuffer = ibuffer + 1
                        ppm_sendbuffers(ibuffer) =    &
     &                     REAL(fdata(4,imesh,jmesh,isub),ppm_kind_single)
                        ibuffer = ibuffer + 1
                        ppm_sendbuffers(ibuffer) =    &
     &                     AIMAG(fdata(4,imesh,jmesh,isub))
                        ibuffer = ibuffer + 1
                        ppm_sendbuffers(ibuffer) =    &
     &                     REAL(fdata(5,imesh,jmesh,isub),ppm_kind_single)
                        ibuffer = ibuffer + 1
                        ppm_sendbuffers(ibuffer) =    &
     &                     AIMAG(fdata(5,imesh,jmesh,isub))
#elif  __KIND == __DOUBLE_PRECISION_COMPLEX
                        ppm_sendbuffers(ibuffer) =    &
     &                      REAL(fdata(1,imesh,jmesh,isub),ppm_kind_single)
                        ibuffer = ibuffer + 1
                        ppm_sendbuffers(ibuffer) =    &
     &                      REAL(AIMAG(fdata(1,imesh,jmesh,isub)),   &
     &                      ppm_kind_single)
                        ibuffer = ibuffer + 1
                        ppm_sendbuffers(ibuffer) =    &
     &                      REAL(fdata(2,imesh,jmesh,isub),ppm_kind_single)
                        ibuffer = ibuffer + 1
                        ppm_sendbuffers(ibuffer) =    &
     &                      REAL(AIMAG(fdata(2,imesh,jmesh,isub)),   &
     &                      ppm_kind_single)
                        ibuffer = ibuffer + 1
                        ppm_sendbuffers(ibuffer) =    &
     &                      REAL(fdata(3,imesh,jmesh,isub),ppm_kind_single)
                        ibuffer = ibuffer + 1
                        ppm_sendbuffers(ibuffer) =    &
     &                      REAL(AIMAG(fdata(3,imesh,jmesh,isub)),   &
     &                      ppm_kind_single)
                        ibuffer = ibuffer + 1
                        ppm_sendbuffers(ibuffer) =    &
     &                      REAL(fdata(4,imesh,jmesh,isub),ppm_kind_single)
                        ibuffer = ibuffer + 1
                        ppm_sendbuffers(ibuffer) =    &
     &                      REAL(AIMAG(fdata(4,imesh,jmesh,isub)),   &
     &                      ppm_kind_single)
                        ibuffer = ibuffer + 1
                        ppm_sendbuffers(ibuffer) =    &
     &                      REAL(fdata(5,imesh,jmesh,isub),ppm_kind_single)
                        ibuffer = ibuffer + 1
                        ppm_sendbuffers(ibuffer) =    &
     &                      REAL(AIMAG(fdata(5,imesh,jmesh,isub)),   &
     &                      ppm_kind_single)
#elif  __KIND == __INTEGER
                        ppm_sendbuffers(ibuffer) =    &
     &                      REAL(fdata(1,imesh,jmesh,isub),ppm_kind_single)
                        ibuffer = ibuffer + 1
                        ppm_sendbuffers(ibuffer) =    &
     &                      REAL(fdata(2,imesh,jmesh,isub),ppm_kind_single)
                        ibuffer = ibuffer + 1
                        ppm_sendbuffers(ibuffer) =    &
     &                      REAL(fdata(3,imesh,jmesh,isub),ppm_kind_single)
                        ibuffer = ibuffer + 1
                        ppm_sendbuffers(ibuffer) =    &
     &                      REAL(fdata(4,imesh,jmesh,isub),ppm_kind_single)
                        ibuffer = ibuffer + 1
                        ppm_sendbuffers(ibuffer) =    &
     &                      REAL(fdata(5,imesh,jmesh,isub),ppm_kind_single)
#elif  __KIND == __LOGICAL
                        IF (fdata(1,imesh,jmesh,isub)) THEN
                           ppm_sendbuffers(ibuffer) = 1.0_ppm_kind_single
                        ELSE
                           ppm_sendbuffers(ibuffer) = 0.0_ppm_kind_single
                        ENDIF
                        ibuffer = ibuffer + 1
                        IF (fdata(2,imesh,jmesh,isub)) THEN
                           ppm_sendbuffers(ibuffer) = 1.0_ppm_kind_single
                        ELSE
                           ppm_sendbuffers(ibuffer) = 0.0_ppm_kind_single
                        ENDIF
                        ibuffer = ibuffer + 1
                        IF (fdata(3,imesh,jmesh,isub)) THEN
                           ppm_sendbuffers(ibuffer) = 1.0_ppm_kind_single
                        ELSE
                           ppm_sendbuffers(ibuffer) = 0.0_ppm_kind_single
                        ENDIF
                        ibuffer = ibuffer + 1
                        IF (fdata(4,imesh,jmesh,isub)) THEN
                           ppm_sendbuffers(ibuffer) = 1.0_ppm_kind_single
                        ELSE
                           ppm_sendbuffers(ibuffer) = 0.0_ppm_kind_single
                        ENDIF
                        ibuffer = ibuffer + 1
                        IF (fdata(5,imesh,jmesh,isub)) THEN
                           ppm_sendbuffers(ibuffer) = 1.0_ppm_kind_single
                        ELSE
                           ppm_sendbuffers(ibuffer) = 0.0_ppm_kind_single
                        ENDIF
#endif
                     ENDDO
                  ENDDO
               !----------------------------------------------------------------
               !  For lda.GT.5 the vector length is lda !!
               !----------------------------------------------------------------
               ELSE
                  DO jmesh=ylo,yhi
                     DO imesh=xlo,xhi
                        DO k=1,lda
                           ibuffer = ibuffer + 1
#if    __KIND == __DOUBLE_PRECISION
                           ppm_sendbuffers(ibuffer) =    &
     &                         REAL(fdata(k,imesh,jmesh,isub),ppm_kind_single)
#elif  __KIND == __SINGLE_PRECISION
                           ppm_sendbuffers(ibuffer) =    &
     &                         fdata(k,imesh,jmesh,isub)
#elif  __KIND == __SINGLE_PRECISION_COMPLEX
                           ppm_sendbuffers(ibuffer) =    &
     &                         REAL(fdata(k,imesh,jmesh,isub),ppm_kind_single)
                           ibuffer = ibuffer + 1
                           ppm_sendbuffers(ibuffer) =    &
     &                         AIMAG(fdata(k,imesh,jmesh,isub))
#elif  __KIND == __DOUBLE_PRECISION_COMPLEX
                           ppm_sendbuffers(ibuffer) =    &
     &                         REAL(fdata(k,imesh,jmesh,isub),ppm_kind_single)
                           ibuffer = ibuffer + 1
                           ppm_sendbuffers(ibuffer) =    &
     &                         REAL(AIMAG(fdata(k,imesh,jmesh,isub)),   &
     &                         ppm_kind_single)
#elif  __KIND == __INTEGER
                           ppm_sendbuffers(ibuffer) =    &
     &                         REAL(fdata(k,imesh,jmesh,isub),ppm_kind_single)
#elif  __KIND == __LOGICAL
                           IF (fdata(k,imesh,jmesh,isub)) THEN
                              ppm_sendbuffers(ibuffer) = 1.0_ppm_kind_single
                           ELSE
                              ppm_sendbuffers(ibuffer) = 0.0_ppm_kind_single
                           ENDIF
#endif
                        ENDDO
                     ENDDO
                  ENDDO
               ENDIF
#elif  __DIM == __SFIELD
               !----------------------------------------------------------------
               !  Serial version
               !----------------------------------------------------------------
               DO jmesh=ylo,yhi
                  DO imesh=xlo,xhi
                     ibuffer = ibuffer + 1
#if    __KIND == __DOUBLE_PRECISION
                     ppm_sendbuffers(ibuffer) =    &
     &                   REAL(fdata(imesh,jmesh,isub),ppm_kind_single)
#elif  __KIND == __SINGLE_PRECISION
                     ppm_sendbuffers(ibuffer) = fdata(imesh,jmesh,isub)
#elif  __KIND == __SINGLE_PRECISION_COMPLEX
                     ppm_sendbuffers(ibuffer) =    &
     &                   REAL(fdata(imesh,jmesh,isub),ppm_kind_single)
                     ibuffer = ibuffer + 1
                     ppm_sendbuffers(ibuffer) = AIMAG(fdata(imesh,jmesh,isub))
#elif  __KIND == __DOUBLE_PRECISION_COMPLEX
                     ppm_sendbuffers(ibuffer) =    &
     &                   REAL(fdata(imesh,jmesh,isub),ppm_kind_single)
                     ibuffer = ibuffer + 1
                     ppm_sendbuffers(ibuffer) =    &
     &                   REAL(AIMAG(fdata(imesh,jmesh,isub)),ppm_kind_single)
#elif  __KIND == __INTEGER
                     ppm_sendbuffers(ibuffer) =    &
     &                   REAL(fdata(imesh,jmesh,isub),ppm_kind_single)
#elif  __KIND == __LOGICAL
                     IF (fdata(imesh,jmesh,isub)) THEN
                        ppm_sendbuffers(ibuffer) = 1.0_ppm_kind_single
                     ELSE
                        ppm_sendbuffers(ibuffer) = 0.0_ppm_kind_single
                     ENDIF
#endif
                  ENDDO
               ENDDO
#endif
            ENDDO       ! ppm_psendbuffer
         ENDDO          ! i=1,ppm_nsendlist
      ENDIF             ! ppm_kind

      !-------------------------------------------------------------------------
      !  Update the buffer data count
      !-------------------------------------------------------------------------
      ppm_nsendbuffer = ibuffer

      !-------------------------------------------------------------------------
      !  Deallocate inverse sub list
      !-------------------------------------------------------------------------
      iopt   = ppm_param_dealloc
      CALL ppm_alloc(invsublist,ldu,iopt,info)
      IF (info .NE. 0) THEN
          info = ppm_error_error
          CALL ppm_error(ppm_err_dealloc,'ppm_map_field_push_2d',     &
     &        'inverse sub list INVSUBLIST',__LINE__,info)
      ENDIF

      !-------------------------------------------------------------------------
      !  Return 
      !-------------------------------------------------------------------------
 9999 CONTINUE
      CALL substop('ppm_map_field_push_2d',t0,info)
      RETURN
#if    __DIM == __SFIELD
#if    __KIND == __SINGLE_PRECISION
      END SUBROUTINE ppm_map_field_push_2d_sca_s
#elif  __KIND == __DOUBLE_PRECISION
      END SUBROUTINE ppm_map_field_push_2d_sca_d
#elif  __KIND == __SINGLE_PRECISION_COMPLEX
      END SUBROUTINE ppm_map_field_push_2d_sca_sc
#elif  __KIND == __DOUBLE_PRECISION_COMPLEX
      END SUBROUTINE ppm_map_field_push_2d_sca_dc
#elif  __KIND == __INTEGER
      END SUBROUTINE ppm_map_field_push_2d_sca_i
#elif  __KIND == __LOGICAL
      END SUBROUTINE ppm_map_field_push_2d_sca_l
#endif

#elif  __DIM == __VFIELD
#if    __KIND == __SINGLE_PRECISION
      END SUBROUTINE ppm_map_field_push_2d_vec_s
#elif  __KIND == __DOUBLE_PRECISION
      END SUBROUTINE ppm_map_field_push_2d_vec_d
#elif  __KIND == __SINGLE_PRECISION_COMPLEX
      END SUBROUTINE ppm_map_field_push_2d_vec_sc
#elif  __KIND == __DOUBLE_PRECISION_COMPLEX
      END SUBROUTINE ppm_map_field_push_2d_vec_dc
#elif  __KIND == __INTEGER
      END SUBROUTINE ppm_map_field_push_2d_vec_i
#elif  __KIND == __LOGICAL
      END SUBROUTINE ppm_map_field_push_2d_vec_l
#endif
#endif
