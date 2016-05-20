      !-------------------------------------------------------------------------
      !  Subroutine   :               ppm_map_field_pop_2d
      !-------------------------------------------------------------------------
      ! Copyright (c) 2012 CSE Lab (ETH Zurich), MOSAIC Group (ETH Zurich),
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


#if    __DIM == __SFIELD
#if    __KIND == __SINGLE_PRECISION
      SUBROUTINE ppm_map_field_pop_2d_sca_s(target_topoid,target_meshid,fdata, &
     &                   ghostsize,info,mask,poptype)
#elif  __KIND == __DOUBLE_PRECISION
      SUBROUTINE ppm_map_field_pop_2d_sca_d(target_topoid,target_meshid,fdata, &
     &                   ghostsize,info,mask,poptype)
#elif  __KIND == __SINGLE_PRECISION_COMPLEX
      SUBROUTINE ppm_map_field_pop_2d_sca_sc(target_topoid,target_meshid,fdata,&
     &                   ghostsize,info,mask,poptype)
#elif  __KIND == __DOUBLE_PRECISION_COMPLEX
      SUBROUTINE ppm_map_field_pop_2d_sca_dc(target_topoid,target_meshid,fdata,&
     &                   ghostsize,info,mask,poptype)
#elif  __KIND == __INTEGER
      SUBROUTINE ppm_map_field_pop_2d_sca_i(target_topoid,target_meshid,fdata, &
     &                   ghostsize,info,mask,poptype)
#elif  __KIND == __LOGICAL
      SUBROUTINE ppm_map_field_pop_2d_sca_l(target_topoid,target_meshid,fdata, &
     &                   ghostsize,info,mask,poptype)
#endif
#elif  __DIM == __VFIELD
#if    __KIND == __SINGLE_PRECISION
      SUBROUTINE ppm_map_field_pop_2d_vec_s(target_topoid,target_meshid,fdata, &
     &                   lda,ghostsize,info,mask,poptype)
#elif  __KIND == __DOUBLE_PRECISION
      SUBROUTINE ppm_map_field_pop_2d_vec_d(target_topoid,target_meshid,fdata, &
     &                   lda,ghostsize,info,mask,poptype)
#elif  __KIND == __SINGLE_PRECISION_COMPLEX
      SUBROUTINE ppm_map_field_pop_2d_vec_sc(target_topoid,target_meshid,fdata,&
     &                   lda,ghostsize,info,mask,poptype)
#elif  __KIND == __DOUBLE_PRECISION_COMPLEX
      SUBROUTINE ppm_map_field_pop_2d_vec_dc(target_topoid,target_meshid,fdata,&
     &                   lda,ghostsize,info,mask,poptype)
#elif  __KIND == __INTEGER
      SUBROUTINE ppm_map_field_pop_2d_vec_i(target_topoid,target_meshid,fdata, &
     &                   lda,ghostsize,info,mask,poptype)
#elif  __KIND == __LOGICAL
      SUBROUTINE ppm_map_field_pop_2d_vec_l(target_topoid,target_meshid,fdata, &
     &                   lda,ghostsize,info,mask,poptype)
#endif
#endif
      !!! This routine pops the list buffer for 2D mesh data
      !!!
      !!! [NOTE]
      !!! The on-processor data is in the first part of the buffer
      !!! Reallocates the fdata array if needed. The user does not need
      !!! to care about its size. If we want to change this to improve
      !!! callability from C, argument checks should be added for the
      !!! SIZE of fdata. If the call to invert_list is too slow here, this
      !!! could be done by topo_store and stored for each topology.
      !!! At the expense of some extra memory.
      !!!
      !!! [WARNING]
      !!! DIM is the dimension of the `fdata` array and not the space
      !!! dimemsion `ppm_dim`!
      !-------------------------------------------------------------------------
      !  Includes
      !-------------------------------------------------------------------------

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
      USE ppm_module_write
      USE ppm_module_check_id
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
      !!! Field data.
      !!!
      !!! 1st index: lda                                                       +
      !!! 2nd-3th index: mesh (i,j) relative to istart-1 of the sub
      !!! (i.e. (i,j)=(1,1) corresponds to istart)                             +
      !!! 4th: isub 1...nsublist (all subs on this processor).
      !!!
      !!! For scalar fields, the first index is omitted (the others shift
      !!! accordingly).
      LOGICAL, DIMENSION(:,:,:),   POINTER, OPTIONAL :: mask
      !!! Logical mask.
      !!!
      !!! Only the mesh nodes for which this is .TRUE. will be
      !!! mapped. If not given, all points are mapped.
      !!!
      !!! 1st-2nd index: mesh (i,j)                                            +
      !!! 3rd: isub.
      INTEGER                        , OPTIONAL      :: poptype
#if   __DIM == __VFIELD
      INTEGER                        , INTENT(IN   ) :: lda
      !!! The leading dimension of the fdata.
      !!! lda=1 for the case of scalar data
#endif
      INTEGER                        , INTENT(IN   ) :: target_topoid
      !!! Topology ID of destination
      INTEGER                        , INTENT(IN   ) :: target_meshid
      !!! Mesh ID of destination
      INTEGER, DIMENSION(:)          , INTENT(IN   ) :: ghostsize
      !!! Size of the ghost layer (in mesh points) in each dimension.
      INTEGER                        , INTENT(  OUT) :: info
      !!! Returns status, 0 upon success
      !-------------------------------------------------------------------------
      !  Local variables
      !-------------------------------------------------------------------------
      INTEGER, DIMENSION(3)   :: mofs
      INTEGER, DIMENSION(4)   :: ldu,ldl
      INTEGER                 :: i,j,k,ibuffer,Mdata,isub,bdim,jsub,edim
      INTEGER                 :: iopt,totopo,xhi,yhi,xlo,ylo,imesh,jmesh
      INTEGER                 :: btype,idom
      REAL(MK)                :: t0
      LOGICAL                 :: ldo
      INTEGER                 :: rtype
      CHARACTER(LEN=ppm_char) :: mesg
#if   __DIM == __SFIELD
      INTEGER, PARAMETER    :: lda = 1
#endif
      TYPE(ppm_t_topo)     , POINTER  :: target_topo
      TYPE(ppm_t_equi_mesh), POINTER  :: target_mesh
      !-------------------------------------------------------------------------
      !  Externals
      !-------------------------------------------------------------------------

      !-------------------------------------------------------------------------
      !  Initialise
      !-------------------------------------------------------------------------
      CALL substart('ppm_map_field_pop_2d',t0,info)


      !-------------------------------------------------------------------------
      !  Check arguments
      !-------------------------------------------------------------------------
      IF (ppm_debug .GT. 0) THEN
        CALL check
        IF (info .NE. 0) GOTO 9999
      ENDIF

      ! skip if buffer empty
      IF (ppm_buffer_set .LT. 1) THEN
        info = ppm_error_notice
        IF (ppm_debug .GT. 1) THEN
            CALL ppm_error(ppm_err_buffer_empt,'ppm_map_field_pop',    &
     &          'Buffer is empty: skipping pop!',__LINE__,info)
        ENDIF
        GOTO 9999
      ENDIF

      !-------------------------------------------------------------------------
      !  set the local pointers to the topology and mesh
      !-------------------------------------------------------------------------
      target_topo => ppm_topo(target_topoid)%t
      target_mesh => target_topo%mesh(target_meshid)

      IF (ppm_debug .GT. 0) THEN
         !----------------------------------------------------------------------
         !  now, if in debug mode, check the rest
         !----------------------------------------------------------------------
         CALL check_two
      ENDIF


      !-------------------------------------------------------------------------
      !  Check what kind of pop is needed (replace or add)
      !-------------------------------------------------------------------------
      IF (ppm_map_type .EQ. ppm_param_map_init) THEN
        info = ppm_error_error
        CALL ppm_error(ppm_err_argument,'ppm_map_field_pop_2d',  &
     &       'map_field_pop cannot be called after ghost_init',__LINE__,info)
        CALL ppm_error(ppm_err_argument,'ppm_map_field_pop_2d',  &
     &       'ghost_init must not be called directly by the user',__LINE__,info)
        GOTO 9999
      ELSE
          IF (PRESENT(poptype)) THEN
              rtype = poptype
          ELSE
              IF (ppm_map_type .EQ. ppm_param_map_ghost_put) THEN
                rtype = ppm_param_pop_add
              ELSE
                rtype = ppm_param_pop_replace
              ENDIF
          ENDIF
      ENDIF

      !-------------------------------------------------------------------------
      !  Check that the required dimension fits the dimension of the buffer
      !-------------------------------------------------------------------------
      bdim = ppm_buffer_dim(ppm_buffer_set)
#if   __KIND == __SINGLE_PRECISION_COMPLEX | __KIND == __DOUBLE_PRECISION_COMPLEX
      ! for complex, the effective dimension is half the data dimension
      edim = bdim/2
#else
      edim = bdim
#endif
      IF (ppm_debug .GT. 1) THEN
          WRITE(mesg,'(2(A,I3))') 'bdim=',edim,'    lda=',lda
          CALL ppm_write(ppm_rank,'ppm_map_field_pop_2d',mesg,info)
      ENDIF
#if   __DIM == __VFIELD
      IF (edim.NE.lda) THEN
         info = ppm_error_error
         CALL ppm_error(ppm_err_wrong_dim,'ppm_map_field_pop_2d',    &
     &       'leading dimension LDA is in error',__LINE__,info)
         GOTO 9999
      ENDIF
#elif __DIM == __SFIELD
      IF (edim.NE.1) THEN
         info = ppm_error_error
         CALL ppm_error(ppm_err_wrong_dim,'ppm_map_field_pop_2d',    &
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
         CALL ppm_error(ppm_err_wrong_prec,'ppm_map_field_pop_2d',    &
     &       'trying to pop a non-single into single ',__LINE__,info)
         GOTO 9999
      ENDIF
#elif  __KIND == __DOUBLE_PRECISION
      IF (btype.NE.ppm_kind_double) THEN
         info = ppm_error_error
         CALL ppm_error(ppm_err_wrong_prec,'ppm_map_field_pop_2d',    &
     &       'trying to pop a non-double into double ',__LINE__,info)
         GOTO 9999
      ENDIF
#elif  __KIND == __SINGLE_PRECISION_COMPLEX
      IF (btype.NE.ppm_kind_single) THEN
         info = ppm_error_error
         CALL ppm_error(ppm_err_wrong_prec,'ppm_map_field_pop_2d',    &
     &       'trying to pop a non-single-complex into single-complex',&
     &       __LINE__,info)
         GOTO 9999
      ENDIF
#elif  __KIND == __DOUBLE_PRECISION_COMPLEX
      IF (btype.NE.ppm_kind_double) THEN
         info = ppm_error_error
         CALL ppm_error(ppm_err_wrong_prec,'ppm_map_field_pop_2d',    &
     &       'trying to pop a non-double-complex into double-complex',&
     &       __LINE__,info)
         GOTO 9999
      ENDIF
#elif  __KIND == __INTEGER
      IF (btype.NE.ppm_integer) THEN
         info = ppm_error_error
         CALL ppm_error(ppm_err_wrong_prec,'ppm_map_field_pop_2d',    &
     &       'trying to pop a non-integer into integer ',__LINE__,info)
         GOTO 9999
      ENDIF
#elif  __KIND == __LOGICAL
      IF (btype.NE.ppm_logical) THEN
         info = ppm_error_error
         CALL ppm_error(ppm_err_wrong_prec,'ppm_map_field_pop_2d',    &
     &       'trying to pop a non-logical into logical ',__LINE__,info)
      ENDIF
#endif

      !-------------------------------------------------------------------------
      !  Determine size of field data array needed
      !-------------------------------------------------------------------------
      xhi = 0
      yhi = 0
      DO i=1,target_topo%nsublist
          idom = target_topo%isublist(i)
          IF (target_mesh%nnodes(1,idom).GT.xhi) THEN
              xhi = target_mesh%nnodes(1,idom)
          ENDIF
          IF (target_mesh%nnodes(2,idom).GT.yhi) THEN
              yhi = target_mesh%nnodes(2,idom)
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
      ldu(4) = target_topo%nsublist
      ldl(1) = 1
      ldl(2) = 1-ghostsize(1)
      ldl(3) = 1-ghostsize(2)
      ldl(4) = 1
#elif __DIM == __SFIELD
      ldu(1) = xhi+ghostsize(1)
      ldu(2) = yhi+ghostsize(2)
      ldu(3) = target_topo%nsublist
      ldl(1) = 1-ghostsize(1)
      ldl(2) = 1-ghostsize(2)
      ldl(3) = 1
#endif
      CALL ppm_alloc(fdata,ldl,ldu,iopt,info)
      IF (info .NE. 0) THEN
          info = ppm_error_fatal
          CALL ppm_error(ppm_err_alloc,'ppm_map_field_pop_2d',  &
     &        'new data FDATA',__LINE__,info)
          GOTO 9999
      ENDIF

      !-------------------------------------------------------------------------
      !  Build the inverse sub list to find local sub indeices based on
      !  global ones (the global ones are communicated)
      !-------------------------------------------------------------------------
      iopt   = ppm_param_alloc_fit
      ldu(1) = target_topo%nsublist
      CALL ppm_alloc(sublist,ldu,iopt,info)
      IF (info .NE. 0) THEN
          info = ppm_error_fatal
          CALL ppm_error(ppm_err_alloc,'ppm_map_field_pop_2d',  &
     &        'temporary sub list SUBLIST',__LINE__,info)
          GOTO 9999
      ENDIF
      ! We need to copy it into a temp list, since directly using
      ! ppm_isublist(:,ppm_field_topoid) as an argument to invert_list is
      ! not possible since the argument needs to be a POINTER.
      sublist(1:target_topo%nsublist) =     &
     &    target_topo%isublist(1:target_topo%nsublist)
      CALL ppm_util_invert_list(sublist,invsublist,info)
      iopt   = ppm_param_dealloc
      CALL ppm_alloc(sublist,ldu,iopt,info)
      IF (info .NE. 0) THEN
          info = ppm_error_error
          CALL ppm_error(ppm_err_dealloc,'ppm_map_field_pop_2d',     &
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
 !              mofs(1) = p_mesh%istart(1,jsub)-1
 !              mofs(2) = p_mesh%istart(2,jsub)-1
 !              xlo = ppm_mesh_irecvblkstart(1,j)-mofs(1)
 !              ylo = ppm_mesh_irecvblkstart(2,j)-mofs(2)
 !              xhi = xlo+ppm_mesh_irecvblksize(1,j)-1
 !              yhi = ylo+ppm_mesh_irecvblksize(2,j)-1
 !              !----------------------------------------------------------------
 !              !  Count points for which mask is true
 !              !----------------------------------------------------------------
 !              DO jmesh=ylo,yhi
 !                 DO imesh=xlo,xhi
 !                    IF (mask(1,imesh,jmesh,isub)) Mdata = Mdata + 1
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
     &            ppm_mesh_irecvblksize(2,j))
            ENDDO
         ENDDO
 !     ENDIF

      !-------------------------------------------------------------------------
      !  If there is nothing to be sent we are done
      !-------------------------------------------------------------------------
      IF (Mdata .EQ. 0) THEN
          IF (ppm_debug .GT. 0) THEN
              CALL ppm_write(ppm_rank,'ppm_map_field_pop_2d',   &
     &            'There is no data to be received',info)
          ENDIF
          GOTO 8888
      ENDIF

      !-------------------------------------------------------------------------
      !  Decrement the pointer into the receive buffer
      !-------------------------------------------------------------------------
      IF (ppm_debug .GT. 1) THEN
          WRITE(mesg,'(2(A,I9))') 'ppm_nrecvbuffer = ',ppm_nrecvbuffer,   &
     &        ' / Mdata*bdim = ',Mdata*bdim
          CALL ppm_write(ppm_rank,'ppm_map_field_pop_2d',mesg,info)
      ENDIF
      ppm_nrecvbuffer = ppm_nrecvbuffer - Mdata*bdim

      ibuffer = ppm_nrecvbuffer
      IF (ppm_debug .GT. 1) THEN
          WRITE(mesg,'(A,I9)') 'ibuffer = ',ibuffer
          CALL ppm_write(ppm_rank,'ppm_map_field_pop_2d',mesg,info)
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
     &         ppm_mesh_isendblksize(2,j))
         ENDDO
      ENDDO
      ppm_nsendbuffer = ppm_nsendbuffer - ppm_buffer_dim(ppm_buffer_set)*Mdata

      !-------------------------------------------------------------------------
      !  Debug output
      !-------------------------------------------------------------------------
      IF (ppm_debug .GT. 0) THEN
          IF (rtype .EQ. ppm_param_pop_replace) THEN
              CALL ppm_write(ppm_rank,'ppm_map_field_pop_2d',     &
     &           'Replacing current field values',info)
          ELSEIF(rtype .EQ. ppm_param_pop_add) THEN
              CALL ppm_write(ppm_rank,'ppm_map_field_pop_2d',     &
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
               mofs(1) = target_mesh%istart(1,jsub)-1
               mofs(2) = target_mesh%istart(2,jsub)-1
               !----------------------------------------------------------------
               !  Get boundaries of mesh block to be received in local sub
               !  coordinates
               !----------------------------------------------------------------
               xlo = ppm_mesh_irecvblkstart(1,j)-mofs(1)
               ylo = ppm_mesh_irecvblkstart(2,j)-mofs(2)
               xhi = xlo+ppm_mesh_irecvblksize(1,j)-1
               yhi = ylo+ppm_mesh_irecvblksize(2,j)-1
               IF (ppm_debug .GT. 1) THEN
                   WRITE(mesg,'(A,2I4)') 'start: ',             &
     &                 ppm_mesh_irecvblkstart(1,j),ppm_mesh_irecvblkstart(2,j)
                   CALL ppm_write(ppm_rank,'ppm_map_field_pop_2d',mesg,info)
                   WRITE(mesg,'(A,2I4)') 'size: ',             &
     &                 ppm_mesh_irecvblksize(1,j),ppm_mesh_irecvblksize(2,j)
                   CALL ppm_write(ppm_rank,'ppm_map_field_pop_2d',mesg,info)
                   WRITE(mesg,'(A,2I4)') 'mesh offset: ',mofs(1),mofs(2)
                   CALL ppm_write(ppm_rank,'ppm_map_field_pop_2d',mesg,info)
                   WRITE(mesg,'(A,2I4)') 'xlo, xhi: ',xlo,xhi
                   CALL ppm_write(ppm_rank,'ppm_map_field_pop_2d',mesg,info)
                   WRITE(mesg,'(A,2I4)') 'ylo, yhi: ',ylo,yhi
                   CALL ppm_write(ppm_rank,'ppm_map_field_pop_2d',mesg,info)
                   WRITE(mesg,'(A,I1)') 'buffer dim: ',edim
                   CALL ppm_write(ppm_rank,'ppm_map_field_pop_2d',mesg,info)
#if   __DIM == __VFIELD
                   WRITE(mesg,'(A,2I4)') 'SIZE(fdata): ',SIZE(fdata,2),  &
     &                 SIZE(fdata,3)
                   CALL ppm_write(ppm_rank,'ppm_map_field_pop_2d',mesg,info)
#elif __DIM == __SFIELD
                   WRITE(mesg,'(A,2I4)') 'SIZE(fdata): ',SIZE(fdata,1),  &
     &                 SIZE(fdata,2)
                   CALL ppm_write(ppm_rank,'ppm_map_field_pop_2d',mesg,info)
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
                        DO jmesh=ylo,yhi
                           DO imesh=xlo,xhi
#if    __KIND == __SINGLE_PRECISION_COMPLEX | __KIND == __DOUBLE_PRECISION_COMPLEX
                              ibuffer = ibuffer + 2
#else
                              ibuffer = ibuffer + 1
#endif
#if    __KIND == __SINGLE_PRECISION
                              fdata(1,imesh,jmesh,isub) =     &
     &                           REAL(ppm_recvbufferd(ibuffer),ppm_kind_single)
#elif  __KIND == __DOUBLE_PRECISION
                              fdata(1,imesh,jmesh,isub) =     &
     &                           ppm_recvbufferd(ibuffer)
#elif  __KIND == __SINGLE_PRECISION_COMPLEX
                              fdata(1,imesh,jmesh,isub) =     &
     &                           CMPLX(ppm_recvbufferd(ibuffer-1),  &
     &                           ppm_recvbufferd(ibuffer),ppm_kind_single)
#elif  __KIND == __DOUBLE_PRECISION_COMPLEX
                              fdata(1,imesh,jmesh,isub) =     &
     &                           CMPLX(ppm_recvbufferd(ibuffer-1),  &
     &                           ppm_recvbufferd(ibuffer),ppm_kind_double)
#elif  __KIND == __INTEGER
                              fdata(1,imesh,jmesh,isub) =     &
     &                           INT(ppm_recvbufferd(ibuffer))
#elif  __KIND == __LOGICAL
                              IF (ppm_recvbufferd(ibuffer) .GT.     &
     &                           (1.0_ppm_kind_double-ppm_myepsd)) THEN
                                 fdata(1,imesh,jmesh,isub) = .TRUE.
                              ELSE
                                 fdata(1,imesh,jmesh,isub) = .FALSE.
                              ENDIF
#endif
                           ENDDO
                        ENDDO
                     ELSEIF (rtype .EQ. ppm_param_pop_add) THEN
                        DO jmesh=ylo,yhi
                           DO imesh=xlo,xhi
#if    __KIND == __SINGLE_PRECISION_COMPLEX | __KIND == __DOUBLE_PRECISION_COMPLEX
                              ibuffer = ibuffer + 2
#else
                              ibuffer = ibuffer + 1
#endif
#if    __KIND == __SINGLE_PRECISION
                              fdata(1,imesh,jmesh,isub) =   &
     &                           fdata(1,imesh,jmesh,isub)+ &
     &                           REAL(ppm_recvbufferd(ibuffer),ppm_kind_single)
#elif  __KIND == __DOUBLE_PRECISION
                              fdata(1,imesh,jmesh,isub) =   &
     &                           fdata(1,imesh,jmesh,isub)+ &
     &                           ppm_recvbufferd(ibuffer)
#elif  __KIND == __SINGLE_PRECISION_COMPLEX
                              fdata(1,imesh,jmesh,isub) =   &
     &                           fdata(1,imesh,jmesh,isub)+ &
     &                           CMPLX(ppm_recvbufferd(ibuffer-1),  &
     &                           ppm_recvbufferd(ibuffer),ppm_kind_single)
#elif  __KIND == __DOUBLE_PRECISION_COMPLEX
                              fdata(1,imesh,jmesh,isub) =   &
     &                           fdata(1,imesh,jmesh,isub)+ &
     &                           CMPLX(ppm_recvbufferd(ibuffer-1),  &
     &                           ppm_recvbufferd(ibuffer),ppm_kind_double)
#elif  __KIND == __INTEGER
                              fdata(1,imesh,jmesh,isub) =   &
     &                           fdata(1,imesh,jmesh,isub)+ &
     &                           INT(ppm_recvbufferd(ibuffer))
#elif  __KIND == __LOGICAL
                              IF (ppm_recvbufferd(ibuffer) .GT.     &
     &                           (1.0_ppm_kind_double-ppm_myepsd)) THEN
                                 fdata(1,imesh,jmesh,isub) =    &
     &                              (fdata(1,imesh,jmesh,isub) .AND. .TRUE.)
                              ELSE
                                 fdata(1,imesh,jmesh,isub) =    &
     &                              (fdata(1,imesh,jmesh,isub) .AND. .FALSE.)
                              ENDIF
#endif
                           ENDDO
                        ENDDO
                     ENDIF
                  !-------------------------------------------------------------
                  !  Unrolled for edim=2
                  !-------------------------------------------------------------
                  ELSEIF (edim .EQ. 2) THEN
                     IF (rtype .EQ. ppm_param_pop_replace) THEN
                        DO jmesh=ylo,yhi
                           DO imesh=xlo,xhi
#if    __KIND == __SINGLE_PRECISION_COMPLEX | __KIND == __DOUBLE_PRECISION_COMPLEX
                              ibuffer = ibuffer + 2
#else
                              ibuffer = ibuffer + 1
#endif
#if    __KIND == __SINGLE_PRECISION
                              fdata(1,imesh,jmesh,isub) =     &
     &                           REAL(ppm_recvbufferd(ibuffer),ppm_kind_single)
                              ibuffer = ibuffer + 1
                              fdata(2,imesh,jmesh,isub) =     &
     &                           REAL(ppm_recvbufferd(ibuffer),ppm_kind_single)
#elif  __KIND == __DOUBLE_PRECISION
                              fdata(1,imesh,jmesh,isub) =     &
     &                           ppm_recvbufferd(ibuffer)
                              ibuffer = ibuffer + 1
                              fdata(2,imesh,jmesh,isub) =     &
     &                           ppm_recvbufferd(ibuffer)
#elif  __KIND == __SINGLE_PRECISION_COMPLEX
                              fdata(1,imesh,jmesh,isub) =     &
     &                           CMPLX(ppm_recvbufferd(ibuffer-1),  &
     &                           ppm_recvbufferd(ibuffer),ppm_kind_single)
                              ibuffer = ibuffer + 2
                              fdata(2,imesh,jmesh,isub) =     &
     &                           CMPLX(ppm_recvbufferd(ibuffer-1),  &
     &                           ppm_recvbufferd(ibuffer),ppm_kind_single)
#elif  __KIND == __DOUBLE_PRECISION_COMPLEX
                              fdata(1,imesh,jmesh,isub) =     &
     &                           CMPLX(ppm_recvbufferd(ibuffer-1),  &
     &                           ppm_recvbufferd(ibuffer),ppm_kind_double)
                              ibuffer = ibuffer + 2
                              fdata(2,imesh,jmesh,isub) =     &
     &                           CMPLX(ppm_recvbufferd(ibuffer-1),  &
     &                           ppm_recvbufferd(ibuffer),ppm_kind_double)
#elif  __KIND == __INTEGER
                              fdata(1,imesh,jmesh,isub) =     &
     &                           INT(ppm_recvbufferd(ibuffer))
                              ibuffer = ibuffer + 1
                              fdata(2,imesh,jmesh,isub) =     &
     &                           INT(ppm_recvbufferd(ibuffer))
#elif  __KIND == __LOGICAL
                              IF (ppm_recvbufferd(ibuffer) .GT.     &
     &                           (1.0_ppm_kind_double-ppm_myepsd)) THEN
                                 fdata(1,imesh,jmesh,isub) = .TRUE.
                              ELSE
                                 fdata(2,imesh,jmesh,isub) = .FALSE.
                              ENDIF
                              ibuffer = ibuffer + 1
                              IF (ppm_recvbufferd(ibuffer) .GT.     &
     &                           (1.0_ppm_kind_double-ppm_myepsd)) THEN
                                 fdata(1,imesh,jmesh,isub) = .TRUE.
                              ELSE
                                 fdata(2,imesh,jmesh,isub) = .FALSE.
                              ENDIF
#endif
                           ENDDO
                        ENDDO
                     ELSEIF (rtype .EQ. ppm_param_pop_add) THEN
                        DO jmesh=ylo,yhi
                           DO imesh=xlo,xhi
#if    __KIND == __SINGLE_PRECISION_COMPLEX | __KIND == __DOUBLE_PRECISION_COMPLEX
                              ibuffer = ibuffer + 2
#else
                              ibuffer = ibuffer + 1
#endif
#if    __KIND == __SINGLE_PRECISION
                              fdata(1,imesh,jmesh,isub) =   &
     &                           fdata(1,imesh,jmesh,isub)+ &
     &                           REAL(ppm_recvbufferd(ibuffer),ppm_kind_single)
                              ibuffer = ibuffer + 1
                              fdata(2,imesh,jmesh,isub) =   &
     &                           fdata(2,imesh,jmesh,isub)+ &
     &                           REAL(ppm_recvbufferd(ibuffer),ppm_kind_single)
#elif  __KIND == __DOUBLE_PRECISION
                              fdata(1,imesh,jmesh,isub) =   &
     &                           fdata(1,imesh,jmesh,isub)+ &
     &                           ppm_recvbufferd(ibuffer)
                              ibuffer = ibuffer + 1
                              fdata(2,imesh,jmesh,isub) =   &
     &                           fdata(2,imesh,jmesh,isub)+ &
     &                           ppm_recvbufferd(ibuffer)
#elif  __KIND == __SINGLE_PRECISION_COMPLEX
                              fdata(1,imesh,jmesh,isub) =   &
     &                           fdata(1,imesh,jmesh,isub)+ &
     &                           CMPLX(ppm_recvbufferd(ibuffer-1),  &
     &                           ppm_recvbufferd(ibuffer),ppm_kind_single)
                              ibuffer = ibuffer + 2
                              fdata(2,imesh,jmesh,isub) =   &
     &                           fdata(2,imesh,jmesh,isub)+ &
     &                           CMPLX(ppm_recvbufferd(ibuffer-1),  &
     &                           ppm_recvbufferd(ibuffer),ppm_kind_single)
#elif  __KIND == __DOUBLE_PRECISION_COMPLEX
                              fdata(1,imesh,jmesh,isub) =   &
     &                           fdata(1,imesh,jmesh,isub)+ &
     &                           CMPLX(ppm_recvbufferd(ibuffer-1),  &
     &                           ppm_recvbufferd(ibuffer),ppm_kind_double)
                              ibuffer = ibuffer + 2
                              fdata(2,imesh,jmesh,isub) =   &
     &                           fdata(2,imesh,jmesh,isub)+ &
     &                           CMPLX(ppm_recvbufferd(ibuffer-1),  &
     &                           ppm_recvbufferd(ibuffer),ppm_kind_double)
#elif  __KIND == __INTEGER
                              fdata(1,imesh,jmesh,isub) =   &
     &                           fdata(1,imesh,jmesh,isub)+ &
     &                           INT(ppm_recvbufferd(ibuffer))
                              ibuffer = ibuffer + 1
                              fdata(2,imesh,jmesh,isub) =   &
     &                           fdata(2,imesh,jmesh,isub)+ &
     &                           INT(ppm_recvbufferd(ibuffer))
#elif  __KIND == __LOGICAL
                              IF (ppm_recvbufferd(ibuffer) .GT.     &
     &                           (1.0_ppm_kind_double-ppm_myepsd)) THEN
                                 fdata(1,imesh,jmesh,isub) =    &
     &                              (fdata(1,imesh,jmesh,isub) .AND. .TRUE.)
                              ELSE
                                 fdata(1,imesh,jmesh,isub) =    &
     &                              (fdata(1,imesh,jmesh,isub) .AND. .FALSE.)
                              ENDIF
                              ibuffer = ibuffer + 1
                              IF (ppm_recvbufferd(ibuffer) .GT.     &
     &                           (1.0_ppm_kind_double-ppm_myepsd)) THEN
                                 fdata(2,imesh,jmesh,isub) =    &
     &                              (fdata(2,imesh,jmesh,isub) .AND. .TRUE.)
                              ELSE
                                 fdata(2,imesh,jmesh,isub) =    &
     &                              (fdata(2,imesh,jmesh,isub) .AND. .FALSE.)
                              ENDIF
#endif
                           ENDDO
                        ENDDO
                     ENDIF
                  !-------------------------------------------------------------
                  !  Unrolled for edim=3
                  !-------------------------------------------------------------
                  ELSEIF (edim .EQ. 3) THEN
                     IF (rtype .EQ. ppm_param_pop_replace) THEN
                        DO jmesh=ylo,yhi
                           DO imesh=xlo,xhi
#if    __KIND == __SINGLE_PRECISION_COMPLEX | __KIND == __DOUBLE_PRECISION_COMPLEX
                              ibuffer = ibuffer + 2
#else
                              ibuffer = ibuffer + 1
#endif
#if    __KIND == __SINGLE_PRECISION
                              fdata(1,imesh,jmesh,isub) =     &
     &                           REAL(ppm_recvbufferd(ibuffer),ppm_kind_single)
                              ibuffer = ibuffer + 1
                              fdata(2,imesh,jmesh,isub) =     &
     &                           REAL(ppm_recvbufferd(ibuffer),ppm_kind_single)
                              ibuffer = ibuffer + 1
                              fdata(3,imesh,jmesh,isub) =     &
     &                           REAL(ppm_recvbufferd(ibuffer),ppm_kind_single)
#elif  __KIND == __DOUBLE_PRECISION
                              fdata(1,imesh,jmesh,isub) =     &
     &                           ppm_recvbufferd(ibuffer)
                              ibuffer = ibuffer + 1
                              fdata(2,imesh,jmesh,isub) =     &
     &                           ppm_recvbufferd(ibuffer)
                              ibuffer = ibuffer + 1
                              fdata(3,imesh,jmesh,isub) =     &
     &                           ppm_recvbufferd(ibuffer)
#elif  __KIND == __SINGLE_PRECISION_COMPLEX
                              fdata(1,imesh,jmesh,isub) =     &
     &                           CMPLX(ppm_recvbufferd(ibuffer-1),  &
     &                           ppm_recvbufferd(ibuffer),ppm_kind_single)
                              ibuffer = ibuffer + 2
                              fdata(2,imesh,jmesh,isub) =     &
     &                           CMPLX(ppm_recvbufferd(ibuffer-1),  &
     &                           ppm_recvbufferd(ibuffer),ppm_kind_single)
                              ibuffer = ibuffer + 2
                              fdata(3,imesh,jmesh,isub) =     &
     &                           CMPLX(ppm_recvbufferd(ibuffer-1),  &
     &                           ppm_recvbufferd(ibuffer),ppm_kind_single)
#elif  __KIND == __DOUBLE_PRECISION_COMPLEX
                              fdata(1,imesh,jmesh,isub) =     &
     &                           CMPLX(ppm_recvbufferd(ibuffer-1),  &
     &                           ppm_recvbufferd(ibuffer),ppm_kind_double)
                              ibuffer = ibuffer + 2
                              fdata(2,imesh,jmesh,isub) =     &
     &                           CMPLX(ppm_recvbufferd(ibuffer-1),  &
     &                           ppm_recvbufferd(ibuffer),ppm_kind_double)
                              ibuffer = ibuffer + 2
                              fdata(3,imesh,jmesh,isub) =     &
     &                           CMPLX(ppm_recvbufferd(ibuffer-1),  &
     &                           ppm_recvbufferd(ibuffer),ppm_kind_double)
#elif  __KIND == __INTEGER
                              fdata(1,imesh,jmesh,isub) =     &
     &                           INT(ppm_recvbufferd(ibuffer))
                              ibuffer = ibuffer + 1
                              fdata(2,imesh,jmesh,isub) =     &
     &                           INT(ppm_recvbufferd(ibuffer))
                              ibuffer = ibuffer + 1
                              fdata(3,imesh,jmesh,isub) =     &
     &                           INT(ppm_recvbufferd(ibuffer))
#elif  __KIND == __LOGICAL
                              IF (ppm_recvbufferd(ibuffer) .GT.     &
     &                           (1.0_ppm_kind_double-ppm_myepsd)) THEN
                                 fdata(1,imesh,jmesh,isub) = .TRUE.
                              ELSE
                                 fdata(1,imesh,jmesh,isub) = .FALSE.
                              ENDIF
                              ibuffer = ibuffer + 1
                              IF (ppm_recvbufferd(ibuffer) .GT.     &
     &                           (1.0_ppm_kind_double-ppm_myepsd)) THEN
                                 fdata(2,imesh,jmesh,isub) = .TRUE.
                              ELSE
                                 fdata(2,imesh,jmesh,isub) = .FALSE.
                              ENDIF
                              ibuffer = ibuffer + 1
                              IF (ppm_recvbufferd(ibuffer) .GT.     &
     &                           (1.0_ppm_kind_double-ppm_myepsd)) THEN
                                 fdata(3,imesh,jmesh,isub) = .TRUE.
                              ELSE
                                 fdata(3,imesh,jmesh,isub) = .FALSE.
                              ENDIF
#endif
                           ENDDO
                        ENDDO
                     ELSEIF (rtype .EQ. ppm_param_pop_add) THEN
                        DO jmesh=ylo,yhi
                           DO imesh=xlo,xhi
#if    __KIND == __SINGLE_PRECISION_COMPLEX | __KIND == __DOUBLE_PRECISION_COMPLEX
                              ibuffer = ibuffer + 2
#else
                              ibuffer = ibuffer + 1
#endif
#if    __KIND == __SINGLE_PRECISION
                              fdata(1,imesh,jmesh,isub) =   &
     &                           fdata(1,imesh,jmesh,isub)+ &
     &                           REAL(ppm_recvbufferd(ibuffer),ppm_kind_single)
                              ibuffer = ibuffer + 1
                              fdata(2,imesh,jmesh,isub) =   &
     &                           fdata(2,imesh,jmesh,isub)+ &
     &                           REAL(ppm_recvbufferd(ibuffer),ppm_kind_single)
                              ibuffer = ibuffer + 1
                              fdata(3,imesh,jmesh,isub) =   &
     &                           fdata(3,imesh,jmesh,isub)+ &
     &                           REAL(ppm_recvbufferd(ibuffer),ppm_kind_single)
#elif  __KIND == __DOUBLE_PRECISION
                              fdata(1,imesh,jmesh,isub) =   &
     &                           fdata(1,imesh,jmesh,isub)+ &
     &                           ppm_recvbufferd(ibuffer)
                              ibuffer = ibuffer + 1
                              fdata(2,imesh,jmesh,isub) =   &
     &                           fdata(2,imesh,jmesh,isub)+ &
     &                           ppm_recvbufferd(ibuffer)
                              ibuffer = ibuffer + 1
                              fdata(3,imesh,jmesh,isub) =   &
     &                           fdata(3,imesh,jmesh,isub)+ &
     &                           ppm_recvbufferd(ibuffer)
#elif  __KIND == __SINGLE_PRECISION_COMPLEX
                              fdata(1,imesh,jmesh,isub) =   &
     &                           fdata(1,imesh,jmesh,isub)+ &
     &                           CMPLX(ppm_recvbufferd(ibuffer-1),  &
     &                           ppm_recvbufferd(ibuffer),ppm_kind_single)
                              ibuffer = ibuffer + 2
                              fdata(2,imesh,jmesh,isub) =   &
     &                           fdata(2,imesh,jmesh,isub)+ &
     &                           CMPLX(ppm_recvbufferd(ibuffer-1),  &
     &                           ppm_recvbufferd(ibuffer),ppm_kind_single)
                              ibuffer = ibuffer + 2
                              fdata(3,imesh,jmesh,isub) =   &
     &                           fdata(3,imesh,jmesh,isub)+ &
     &                           CMPLX(ppm_recvbufferd(ibuffer-1),  &
     &                           ppm_recvbufferd(ibuffer),ppm_kind_single)
#elif  __KIND == __DOUBLE_PRECISION_COMPLEX
                              fdata(1,imesh,jmesh,isub) =   &
     &                           fdata(1,imesh,jmesh,isub)+ &
     &                           CMPLX(ppm_recvbufferd(ibuffer-1),  &
     &                           ppm_recvbufferd(ibuffer),ppm_kind_double)
                              ibuffer = ibuffer + 2
                              fdata(2,imesh,jmesh,isub) =   &
     &                           fdata(2,imesh,jmesh,isub)+ &
     &                           CMPLX(ppm_recvbufferd(ibuffer-1),  &
     &                           ppm_recvbufferd(ibuffer),ppm_kind_double)
                              ibuffer = ibuffer + 2
                              fdata(3,imesh,jmesh,isub) =   &
     &                           fdata(3,imesh,jmesh,isub)+ &
     &                           CMPLX(ppm_recvbufferd(ibuffer-1),  &
     &                           ppm_recvbufferd(ibuffer),ppm_kind_double)
#elif  __KIND == __INTEGER
                              fdata(1,imesh,jmesh,isub) =   &
     &                           fdata(1,imesh,jmesh,isub)+ &
     &                           INT(ppm_recvbufferd(ibuffer))
                              ibuffer = ibuffer + 1
                              fdata(2,imesh,jmesh,isub) =   &
     &                           fdata(2,imesh,jmesh,isub)+ &
     &                           INT(ppm_recvbufferd(ibuffer))
                              ibuffer = ibuffer + 1
                              fdata(3,imesh,jmesh,isub) =   &
     &                           fdata(3,imesh,jmesh,isub)+ &
     &                           INT(ppm_recvbufferd(ibuffer))
#elif  __KIND == __LOGICAL
                              IF (ppm_recvbufferd(ibuffer) .GT.     &
     &                           (1.0_ppm_kind_double-ppm_myepsd)) THEN
                                 fdata(1,imesh,jmesh,isub) =    &
     &                              (fdata(1,imesh,jmesh,isub) .AND. .TRUE.)
                              ELSE
                                 fdata(1,imesh,jmesh,isub) =    &
     &                              (fdata(1,imesh,jmesh,isub) .AND. .FALSE.)
                              ENDIF
                              ibuffer = ibuffer + 1
                              IF (ppm_recvbufferd(ibuffer) .GT.     &
     &                           (1.0_ppm_kind_double-ppm_myepsd)) THEN
                                 fdata(2,imesh,jmesh,isub) =    &
     &                              (fdata(2,imesh,jmesh,isub) .AND. .TRUE.)
                              ELSE
                                 fdata(2,imesh,jmesh,isub) =    &
     &                              (fdata(2,imesh,jmesh,isub) .AND. .FALSE.)
                              ENDIF
                              ibuffer = ibuffer + 1
                              IF (ppm_recvbufferd(ibuffer) .GT.     &
     &                           (1.0_ppm_kind_double-ppm_myepsd)) THEN
                                 fdata(3,imesh,jmesh,isub) =    &
     &                              (fdata(3,imesh,jmesh,isub) .AND. .TRUE.)
                              ELSE
                                 fdata(3,imesh,jmesh,isub) =    &
     &                              (fdata(3,imesh,jmesh,isub) .AND. .FALSE.)
                              ENDIF
#endif
                           ENDDO
                        ENDDO
                     ENDIF
                  !-------------------------------------------------------------
                  !  Unrolled for edim=4
                  !-------------------------------------------------------------
                  ELSEIF (edim .EQ. 4) THEN
                     IF (rtype .EQ. ppm_param_pop_replace) THEN
                        DO jmesh=ylo,yhi
                           DO imesh=xlo,xhi
#if    __KIND == __SINGLE_PRECISION_COMPLEX | __KIND == __DOUBLE_PRECISION_COMPLEX
                              ibuffer = ibuffer + 2
#else
                              ibuffer = ibuffer + 1
#endif
#if    __KIND == __SINGLE_PRECISION
                              fdata(1,imesh,jmesh,isub) =     &
     &                           REAL(ppm_recvbufferd(ibuffer),ppm_kind_single)
                              ibuffer = ibuffer + 1
                              fdata(2,imesh,jmesh,isub) =     &
     &                           REAL(ppm_recvbufferd(ibuffer),ppm_kind_single)
                              ibuffer = ibuffer + 1
                              fdata(3,imesh,jmesh,isub) =     &
     &                           REAL(ppm_recvbufferd(ibuffer),ppm_kind_single)
                              ibuffer = ibuffer + 1
                              fdata(4,imesh,jmesh,isub) =     &
     &                           REAL(ppm_recvbufferd(ibuffer),ppm_kind_single)
#elif  __KIND == __DOUBLE_PRECISION
                              fdata(1,imesh,jmesh,isub) =     &
     &                           ppm_recvbufferd(ibuffer)
                              ibuffer = ibuffer + 1
                              fdata(2,imesh,jmesh,isub) =     &
     &                           ppm_recvbufferd(ibuffer)
                              ibuffer = ibuffer + 1
                              fdata(3,imesh,jmesh,isub) =     &
     &                           ppm_recvbufferd(ibuffer)
                              ibuffer = ibuffer + 1
                              fdata(4,imesh,jmesh,isub) =     &
     &                           ppm_recvbufferd(ibuffer)
#elif  __KIND == __SINGLE_PRECISION_COMPLEX
                              fdata(1,imesh,jmesh,isub) =     &
     &                           CMPLX(ppm_recvbufferd(ibuffer-1),  &
     &                           ppm_recvbufferd(ibuffer),ppm_kind_single)
                              ibuffer = ibuffer + 2
                              fdata(2,imesh,jmesh,isub) =     &
     &                           CMPLX(ppm_recvbufferd(ibuffer-1),  &
     &                           ppm_recvbufferd(ibuffer),ppm_kind_single)
                              ibuffer = ibuffer + 2
                              fdata(3,imesh,jmesh,isub) =     &
     &                           CMPLX(ppm_recvbufferd(ibuffer-1),  &
     &                           ppm_recvbufferd(ibuffer),ppm_kind_single)
                              ibuffer = ibuffer + 2
                              fdata(4,imesh,jmesh,isub) =     &
     &                           CMPLX(ppm_recvbufferd(ibuffer-1),  &
     &                           ppm_recvbufferd(ibuffer),ppm_kind_single)
#elif  __KIND == __DOUBLE_PRECISION_COMPLEX
                              fdata(1,imesh,jmesh,isub) =     &
     &                           CMPLX(ppm_recvbufferd(ibuffer-1),  &
     &                           ppm_recvbufferd(ibuffer),ppm_kind_double)
                              ibuffer = ibuffer + 2
                              fdata(2,imesh,jmesh,isub) =     &
     &                           CMPLX(ppm_recvbufferd(ibuffer-1),  &
     &                           ppm_recvbufferd(ibuffer),ppm_kind_double)
                              ibuffer = ibuffer + 2
                              fdata(3,imesh,jmesh,isub) =     &
     &                           CMPLX(ppm_recvbufferd(ibuffer-1),  &
     &                           ppm_recvbufferd(ibuffer),ppm_kind_double)
                              ibuffer = ibuffer + 2
                              fdata(4,imesh,jmesh,isub) =     &
     &                           CMPLX(ppm_recvbufferd(ibuffer-1),  &
     &                           ppm_recvbufferd(ibuffer),ppm_kind_double)
#elif  __KIND == __INTEGER
                              fdata(1,imesh,jmesh,isub) =     &
     &                           INT(ppm_recvbufferd(ibuffer))
                              ibuffer = ibuffer + 1
                              fdata(2,imesh,jmesh,isub) =     &
     &                           INT(ppm_recvbufferd(ibuffer))
                              ibuffer = ibuffer + 1
                              fdata(3,imesh,jmesh,isub) =     &
     &                           INT(ppm_recvbufferd(ibuffer))
                              ibuffer = ibuffer + 1
                              fdata(4,imesh,jmesh,isub) =     &
     &                           INT(ppm_recvbufferd(ibuffer))
#elif  __KIND == __LOGICAL
                              IF (ppm_recvbufferd(ibuffer) .GT.     &
     &                           (1.0_ppm_kind_double-ppm_myepsd)) THEN
                                 fdata(1,imesh,jmesh,isub) = .TRUE.
                              ELSE
                                 fdata(1,imesh,jmesh,isub) = .FALSE.
                              ENDIF
                              ibuffer = ibuffer + 1
                              IF (ppm_recvbufferd(ibuffer) .GT.     &
     &                           (1.0_ppm_kind_double-ppm_myepsd)) THEN
                                 fdata(2,imesh,jmesh,isub) = .TRUE.
                              ELSE
                                 fdata(2,imesh,jmesh,isub) = .FALSE.
                              ENDIF
                              ibuffer = ibuffer + 1
                              IF (ppm_recvbufferd(ibuffer) .GT.     &
     &                           (1.0_ppm_kind_double-ppm_myepsd)) THEN
                                 fdata(3,imesh,jmesh,isub) = .TRUE.
                              ELSE
                                 fdata(3,imesh,jmesh,isub) = .FALSE.
                              ENDIF
                              ibuffer = ibuffer + 1
                              IF (ppm_recvbufferd(ibuffer) .GT.     &
     &                           (1.0_ppm_kind_double-ppm_myepsd)) THEN
                                 fdata(4,imesh,jmesh,isub) = .TRUE.
                              ELSE
                                 fdata(4,imesh,jmesh,isub) = .FALSE.
                              ENDIF
#endif
                           ENDDO
                        ENDDO
                     ELSEIF (rtype .EQ. ppm_param_pop_add) THEN
                        DO jmesh=ylo,yhi
                           DO imesh=xlo,xhi
#if    __KIND == __SINGLE_PRECISION_COMPLEX | __KIND == __DOUBLE_PRECISION_COMPLEX
                              ibuffer = ibuffer + 2
#else
                              ibuffer = ibuffer + 1
#endif
#if    __KIND == __SINGLE_PRECISION
                              fdata(1,imesh,jmesh,isub) =   &
     &                           fdata(1,imesh,jmesh,isub)+ &
     &                           REAL(ppm_recvbufferd(ibuffer),ppm_kind_single)
                              ibuffer = ibuffer + 1
                              fdata(2,imesh,jmesh,isub) =   &
     &                           fdata(2,imesh,jmesh,isub)+ &
     &                           REAL(ppm_recvbufferd(ibuffer),ppm_kind_single)
                              ibuffer = ibuffer + 1
                              fdata(3,imesh,jmesh,isub) =   &
     &                           fdata(3,imesh,jmesh,isub)+ &
     &                           REAL(ppm_recvbufferd(ibuffer),ppm_kind_single)
                              ibuffer = ibuffer + 1
                              fdata(4,imesh,jmesh,isub) =   &
     &                           fdata(4,imesh,jmesh,isub)+ &
     &                           REAL(ppm_recvbufferd(ibuffer),ppm_kind_single)
#elif  __KIND == __DOUBLE_PRECISION
                              fdata(1,imesh,jmesh,isub) =   &
     &                           fdata(1,imesh,jmesh,isub)+ &
     &                           ppm_recvbufferd(ibuffer)
                              ibuffer = ibuffer + 1
                              fdata(2,imesh,jmesh,isub) =   &
     &                           fdata(2,imesh,jmesh,isub)+ &
     &                           ppm_recvbufferd(ibuffer)
                              ibuffer = ibuffer + 1
                              fdata(3,imesh,jmesh,isub) =   &
     &                           fdata(3,imesh,jmesh,isub)+ &
     &                           ppm_recvbufferd(ibuffer)
                              ibuffer = ibuffer + 1
                              fdata(4,imesh,jmesh,isub) =   &
     &                           fdata(4,imesh,jmesh,isub)+ &
     &                           ppm_recvbufferd(ibuffer)
#elif  __KIND == __SINGLE_PRECISION_COMPLEX
                              fdata(1,imesh,jmesh,isub) =   &
     &                           fdata(1,imesh,jmesh,isub)+ &
     &                           CMPLX(ppm_recvbufferd(ibuffer-1),  &
     &                           ppm_recvbufferd(ibuffer),ppm_kind_single)
                              ibuffer = ibuffer + 2
                              fdata(2,imesh,jmesh,isub) =   &
     &                           fdata(2,imesh,jmesh,isub)+ &
     &                           CMPLX(ppm_recvbufferd(ibuffer-1),  &
     &                           ppm_recvbufferd(ibuffer),ppm_kind_single)
                              ibuffer = ibuffer + 2
                              fdata(3,imesh,jmesh,isub) =   &
     &                           fdata(3,imesh,jmesh,isub)+ &
     &                           CMPLX(ppm_recvbufferd(ibuffer-1),  &
     &                           ppm_recvbufferd(ibuffer),ppm_kind_single)
                              ibuffer = ibuffer + 2
                              fdata(4,imesh,jmesh,isub) =   &
     &                           fdata(4,imesh,jmesh,isub)+ &
     &                           CMPLX(ppm_recvbufferd(ibuffer-1),  &
     &                           ppm_recvbufferd(ibuffer),ppm_kind_single)
#elif  __KIND == __DOUBLE_PRECISION_COMPLEX
                              fdata(1,imesh,jmesh,isub) =   &
     &                           fdata(1,imesh,jmesh,isub)+ &
     &                           CMPLX(ppm_recvbufferd(ibuffer-1),  &
     &                           ppm_recvbufferd(ibuffer),ppm_kind_double)
                              ibuffer = ibuffer + 2
                              fdata(2,imesh,jmesh,isub) =   &
     &                           fdata(2,imesh,jmesh,isub)+ &
     &                           CMPLX(ppm_recvbufferd(ibuffer-1),  &
     &                           ppm_recvbufferd(ibuffer),ppm_kind_double)
                              ibuffer = ibuffer + 2
                              fdata(3,imesh,jmesh,isub) =   &
     &                           fdata(3,imesh,jmesh,isub)+ &
     &                           CMPLX(ppm_recvbufferd(ibuffer-1),  &
     &                           ppm_recvbufferd(ibuffer),ppm_kind_double)
                              ibuffer = ibuffer + 2
                              fdata(4,imesh,jmesh,isub) =   &
     &                           fdata(4,imesh,jmesh,isub)+ &
     &                           CMPLX(ppm_recvbufferd(ibuffer-1),  &
     &                           ppm_recvbufferd(ibuffer),ppm_kind_double)
#elif  __KIND == __INTEGER
                              fdata(1,imesh,jmesh,isub) =   &
     &                           fdata(1,imesh,jmesh,isub)+ &
     &                           INT(ppm_recvbufferd(ibuffer))
                              ibuffer = ibuffer + 1
                              fdata(2,imesh,jmesh,isub) =   &
     &                           fdata(2,imesh,jmesh,isub)+ &
     &                           INT(ppm_recvbufferd(ibuffer))
                              ibuffer = ibuffer + 1
                              fdata(3,imesh,jmesh,isub) =   &
     &                           fdata(3,imesh,jmesh,isub)+ &
     &                           INT(ppm_recvbufferd(ibuffer))
                              ibuffer = ibuffer + 1
                              fdata(4,imesh,jmesh,isub) =   &
     &                           fdata(4,imesh,jmesh,isub)+ &
     &                           INT(ppm_recvbufferd(ibuffer))
#elif  __KIND == __LOGICAL
                              IF (ppm_recvbufferd(ibuffer) .GT.     &
     &                           (1.0_ppm_kind_double-ppm_myepsd)) THEN
                                 fdata(1,imesh,jmesh,isub) =    &
     &                              (fdata(1,imesh,jmesh,isub) .AND. .TRUE.)
                              ELSE
                                 fdata(1,imesh,jmesh,isub) =    &
     &                              (fdata(1,imesh,jmesh,isub) .AND. .FALSE.)
                              ENDIF
                              ibuffer = ibuffer + 1
                              IF (ppm_recvbufferd(ibuffer) .GT.     &
     &                           (1.0_ppm_kind_double-ppm_myepsd)) THEN
                                 fdata(2,imesh,jmesh,isub) =    &
     &                              (fdata(2,imesh,jmesh,isub) .AND. .TRUE.)
                              ELSE
                                 fdata(2,imesh,jmesh,isub) =    &
     &                              (fdata(2,imesh,jmesh,isub) .AND. .FALSE.)
                              ENDIF
                              ibuffer = ibuffer + 1
                              IF (ppm_recvbufferd(ibuffer) .GT.     &
     &                           (1.0_ppm_kind_double-ppm_myepsd)) THEN
                                 fdata(3,imesh,jmesh,isub) =    &
     &                              (fdata(3,imesh,jmesh,isub) .AND. .TRUE.)
                              ELSE
                                 fdata(3,imesh,jmesh,isub) =    &
     &                              (fdata(3,imesh,jmesh,isub) .AND. .FALSE.)
                              ENDIF
                              ibuffer = ibuffer + 1
                              IF (ppm_recvbufferd(ibuffer) .GT.     &
     &                           (1.0_ppm_kind_double-ppm_myepsd)) THEN
                                 fdata(4,imesh,jmesh,isub) =    &
     &                              (fdata(4,imesh,jmesh,isub) .AND. .TRUE.)
                              ELSE
                                 fdata(4,imesh,jmesh,isub) =    &
     &                              (fdata(4,imesh,jmesh,isub) .AND. .FALSE.)
                              ENDIF
#endif
                           ENDDO
                        ENDDO
                     ENDIF
                  !-------------------------------------------------------------
                  !  Unrolled for edim=5
                  !-------------------------------------------------------------
                  ELSEIF (edim .EQ. 5) THEN
                     IF (rtype .EQ. ppm_param_pop_replace) THEN
                        DO jmesh=ylo,yhi
                           DO imesh=xlo,xhi
#if    __KIND == __SINGLE_PRECISION_COMPLEX | __KIND == __DOUBLE_PRECISION_COMPLEX
                              ibuffer = ibuffer + 2
#else
                              ibuffer = ibuffer + 1
#endif
#if    __KIND == __SINGLE_PRECISION
                              fdata(1,imesh,jmesh,isub) =     &
     &                           REAL(ppm_recvbufferd(ibuffer),ppm_kind_single)
                              ibuffer = ibuffer + 1
                              fdata(2,imesh,jmesh,isub) =     &
     &                           REAL(ppm_recvbufferd(ibuffer),ppm_kind_single)
                              ibuffer = ibuffer + 1
                              fdata(3,imesh,jmesh,isub) =     &
     &                           REAL(ppm_recvbufferd(ibuffer),ppm_kind_single)
                              ibuffer = ibuffer + 1
                              fdata(4,imesh,jmesh,isub) =     &
     &                           REAL(ppm_recvbufferd(ibuffer),ppm_kind_single)
                              ibuffer = ibuffer + 1
                              fdata(5,imesh,jmesh,isub) =     &
     &                           REAL(ppm_recvbufferd(ibuffer),ppm_kind_single)
#elif  __KIND == __DOUBLE_PRECISION
                              fdata(1,imesh,jmesh,isub) =     &
     &                           ppm_recvbufferd(ibuffer)
                              ibuffer = ibuffer + 1
                              fdata(2,imesh,jmesh,isub) =     &
     &                           ppm_recvbufferd(ibuffer)
                              ibuffer = ibuffer + 1
                              fdata(3,imesh,jmesh,isub) =     &
     &                           ppm_recvbufferd(ibuffer)
                              ibuffer = ibuffer + 1
                              fdata(4,imesh,jmesh,isub) =     &
     &                           ppm_recvbufferd(ibuffer)
                              ibuffer = ibuffer + 1
                              fdata(5,imesh,jmesh,isub) =     &
     &                           ppm_recvbufferd(ibuffer)
#elif  __KIND == __SINGLE_PRECISION_COMPLEX
                              fdata(1,imesh,jmesh,isub) =     &
     &                           CMPLX(ppm_recvbufferd(ibuffer-1),  &
     &                           ppm_recvbufferd(ibuffer),ppm_kind_single)
                              ibuffer = ibuffer + 2
                              fdata(2,imesh,jmesh,isub) =     &
     &                           CMPLX(ppm_recvbufferd(ibuffer-1),  &
     &                           ppm_recvbufferd(ibuffer),ppm_kind_single)
                              ibuffer = ibuffer + 2
                              fdata(3,imesh,jmesh,isub) =     &
     &                           CMPLX(ppm_recvbufferd(ibuffer-1),  &
     &                           ppm_recvbufferd(ibuffer),ppm_kind_single)
                              ibuffer = ibuffer + 2
                              fdata(4,imesh,jmesh,isub) =     &
     &                           CMPLX(ppm_recvbufferd(ibuffer-1),  &
     &                           ppm_recvbufferd(ibuffer),ppm_kind_single)
                              ibuffer = ibuffer + 2
                              fdata(5,imesh,jmesh,isub) =     &
     &                           CMPLX(ppm_recvbufferd(ibuffer-1),  &
     &                           ppm_recvbufferd(ibuffer),ppm_kind_single)
#elif  __KIND == __DOUBLE_PRECISION_COMPLEX
                              fdata(1,imesh,jmesh,isub) =     &
     &                           CMPLX(ppm_recvbufferd(ibuffer-1),  &
     &                           ppm_recvbufferd(ibuffer),ppm_kind_double)
                              ibuffer = ibuffer + 2
                              fdata(2,imesh,jmesh,isub) =     &
     &                           CMPLX(ppm_recvbufferd(ibuffer-1),  &
     &                           ppm_recvbufferd(ibuffer),ppm_kind_double)
                              ibuffer = ibuffer + 2
                              fdata(3,imesh,jmesh,isub) =     &
     &                           CMPLX(ppm_recvbufferd(ibuffer-1),  &
     &                           ppm_recvbufferd(ibuffer),ppm_kind_double)
                              ibuffer = ibuffer + 2
                              fdata(4,imesh,jmesh,isub) =     &
     &                           CMPLX(ppm_recvbufferd(ibuffer-1),  &
     &                           ppm_recvbufferd(ibuffer),ppm_kind_double)
                              ibuffer = ibuffer + 2
                              fdata(5,imesh,jmesh,isub) =     &
     &                           CMPLX(ppm_recvbufferd(ibuffer-1),  &
     &                           ppm_recvbufferd(ibuffer),ppm_kind_double)
#elif  __KIND == __INTEGER
                              fdata(1,imesh,jmesh,isub) =     &
     &                           INT(ppm_recvbufferd(ibuffer))
                              ibuffer = ibuffer + 1
                              fdata(2,imesh,jmesh,isub) =     &
     &                           INT(ppm_recvbufferd(ibuffer))
                              ibuffer = ibuffer + 1
                              fdata(3,imesh,jmesh,isub) =     &
     &                           INT(ppm_recvbufferd(ibuffer))
                              ibuffer = ibuffer + 1
                              fdata(4,imesh,jmesh,isub) =     &
     &                           INT(ppm_recvbufferd(ibuffer))
                              ibuffer = ibuffer + 1
                              fdata(5,imesh,jmesh,isub) =     &
     &                           INT(ppm_recvbufferd(ibuffer))
#elif  __KIND == __LOGICAL
                              IF (ppm_recvbufferd(ibuffer) .GT.     &
     &                           (1.0_ppm_kind_double-ppm_myepsd)) THEN
                                 fdata(1,imesh,jmesh,isub) = .TRUE.
                              ELSE
                                 fdata(1,imesh,jmesh,isub) = .FALSE.
                              ENDIF
                              ibuffer = ibuffer + 1
                              IF (ppm_recvbufferd(ibuffer) .GT.     &
     &                           (1.0_ppm_kind_double-ppm_myepsd)) THEN
                                 fdata(2,imesh,jmesh,isub) = .TRUE.
                              ELSE
                                 fdata(2,imesh,jmesh,isub) = .FALSE.
                              ENDIF
                              ibuffer = ibuffer + 1
                              IF (ppm_recvbufferd(ibuffer) .GT.     &
     &                           (1.0_ppm_kind_double-ppm_myepsd)) THEN
                                 fdata(3,imesh,jmesh,isub) = .TRUE.
                              ELSE
                                 fdata(3,imesh,jmesh,isub) = .FALSE.
                              ENDIF
                              ibuffer = ibuffer + 1
                              IF (ppm_recvbufferd(ibuffer) .GT.     &
     &                           (1.0_ppm_kind_double-ppm_myepsd)) THEN
                                 fdata(4,imesh,jmesh,isub) = .TRUE.
                              ELSE
                                 fdata(4,imesh,jmesh,isub) = .FALSE.
                              ENDIF
                              ibuffer = ibuffer + 1
                              IF (ppm_recvbufferd(ibuffer) .GT.     &
     &                           (1.0_ppm_kind_double-ppm_myepsd)) THEN
                                 fdata(5,imesh,jmesh,isub) = .TRUE.
                              ELSE
                                 fdata(5,imesh,jmesh,isub) = .FALSE.
                              ENDIF
#endif
                           ENDDO
                        ENDDO
                     ELSEIF (rtype .EQ. ppm_param_pop_add) THEN
                        DO jmesh=ylo,yhi
                           DO imesh=xlo,xhi
#if    __KIND == __SINGLE_PRECISION_COMPLEX | __KIND == __DOUBLE_PRECISION_COMPLEX
                              ibuffer = ibuffer + 2
#else
                              ibuffer = ibuffer + 1
#endif
#if    __KIND == __SINGLE_PRECISION
                              fdata(1,imesh,jmesh,isub) =   &
     &                           fdata(1,imesh,jmesh,isub)+ &
     &                           REAL(ppm_recvbufferd(ibuffer),ppm_kind_single)
                              ibuffer = ibuffer + 1
                              fdata(2,imesh,jmesh,isub) =   &
     &                           fdata(2,imesh,jmesh,isub)+ &
     &                           REAL(ppm_recvbufferd(ibuffer),ppm_kind_single)
                              ibuffer = ibuffer + 1
                              fdata(3,imesh,jmesh,isub) =   &
     &                           fdata(3,imesh,jmesh,isub)+ &
     &                           REAL(ppm_recvbufferd(ibuffer),ppm_kind_single)
                              ibuffer = ibuffer + 1
                              fdata(4,imesh,jmesh,isub) =   &
     &                           fdata(4,imesh,jmesh,isub)+ &
     &                           REAL(ppm_recvbufferd(ibuffer),ppm_kind_single)
                              ibuffer = ibuffer + 1
                              fdata(5,imesh,jmesh,isub) =   &
     &                           fdata(5,imesh,jmesh,isub)+ &
     &                           REAL(ppm_recvbufferd(ibuffer),ppm_kind_single)
#elif  __KIND == __DOUBLE_PRECISION
                              fdata(1,imesh,jmesh,isub) =   &
     &                           fdata(1,imesh,jmesh,isub)+ &
     &                           ppm_recvbufferd(ibuffer)
                              ibuffer = ibuffer + 1
                              fdata(2,imesh,jmesh,isub) =   &
     &                           fdata(2,imesh,jmesh,isub)+ &
     &                           ppm_recvbufferd(ibuffer)
                              ibuffer = ibuffer + 1
                              fdata(3,imesh,jmesh,isub) =   &
     &                           fdata(3,imesh,jmesh,isub)+ &
     &                           ppm_recvbufferd(ibuffer)
                              ibuffer = ibuffer + 1
                              fdata(4,imesh,jmesh,isub) =   &
     &                           fdata(4,imesh,jmesh,isub)+ &
     &                           ppm_recvbufferd(ibuffer)
                              ibuffer = ibuffer + 1
                              fdata(5,imesh,jmesh,isub) =   &
     &                           fdata(5,imesh,jmesh,isub)+ &
     &                           ppm_recvbufferd(ibuffer)
#elif  __KIND == __SINGLE_PRECISION_COMPLEX
                              fdata(1,imesh,jmesh,isub) =   &
     &                           fdata(1,imesh,jmesh,isub)+ &
     &                           CMPLX(ppm_recvbufferd(ibuffer-1),  &
     &                           ppm_recvbufferd(ibuffer),ppm_kind_single)
                              ibuffer = ibuffer + 2
                              fdata(2,imesh,jmesh,isub) =   &
     &                           fdata(2,imesh,jmesh,isub)+ &
     &                           CMPLX(ppm_recvbufferd(ibuffer-1),  &
     &                           ppm_recvbufferd(ibuffer),ppm_kind_single)
                              ibuffer = ibuffer + 2
                              fdata(3,imesh,jmesh,isub) =   &
     &                           fdata(3,imesh,jmesh,isub)+ &
     &                           CMPLX(ppm_recvbufferd(ibuffer-1),  &
     &                           ppm_recvbufferd(ibuffer),ppm_kind_single)
                              ibuffer = ibuffer + 2
                              fdata(4,imesh,jmesh,isub) =   &
     &                           fdata(4,imesh,jmesh,isub)+ &
     &                           CMPLX(ppm_recvbufferd(ibuffer-1),  &
     &                           ppm_recvbufferd(ibuffer),ppm_kind_single)
                              ibuffer = ibuffer + 2
                              fdata(5,imesh,jmesh,isub) =   &
     &                           fdata(5,imesh,jmesh,isub)+ &
     &                           CMPLX(ppm_recvbufferd(ibuffer-1),  &
     &                           ppm_recvbufferd(ibuffer),ppm_kind_single)
#elif  __KIND == __DOUBLE_PRECISION_COMPLEX
                              fdata(1,imesh,jmesh,isub) =   &
     &                           fdata(1,imesh,jmesh,isub)+ &
     &                           CMPLX(ppm_recvbufferd(ibuffer-1),  &
     &                           ppm_recvbufferd(ibuffer),ppm_kind_double)
                              ibuffer = ibuffer + 2
                              fdata(2,imesh,jmesh,isub) =   &
     &                           fdata(2,imesh,jmesh,isub)+ &
     &                           CMPLX(ppm_recvbufferd(ibuffer-1),  &
     &                           ppm_recvbufferd(ibuffer),ppm_kind_double)
                              ibuffer = ibuffer + 2
                              fdata(3,imesh,jmesh,isub) =   &
     &                           fdata(3,imesh,jmesh,isub)+ &
     &                           CMPLX(ppm_recvbufferd(ibuffer-1),  &
     &                           ppm_recvbufferd(ibuffer),ppm_kind_double)
                              ibuffer = ibuffer + 2
                              fdata(4,imesh,jmesh,isub) =   &
     &                           fdata(4,imesh,jmesh,isub)+ &
     &                           CMPLX(ppm_recvbufferd(ibuffer-1),  &
     &                           ppm_recvbufferd(ibuffer),ppm_kind_double)
                              ibuffer = ibuffer + 2
                              fdata(5,imesh,jmesh,isub) =   &
     &                           fdata(5,imesh,jmesh,isub)+ &
     &                           CMPLX(ppm_recvbufferd(ibuffer-1),  &
     &                           ppm_recvbufferd(ibuffer),ppm_kind_double)
#elif  __KIND == __INTEGER
                              fdata(1,imesh,jmesh,isub) =   &
     &                           fdata(1,imesh,jmesh,isub)+ &
     &                           INT(ppm_recvbufferd(ibuffer))
                              ibuffer = ibuffer + 1
                              fdata(2,imesh,jmesh,isub) =   &
     &                           fdata(2,imesh,jmesh,isub)+ &
     &                           INT(ppm_recvbufferd(ibuffer))
                              ibuffer = ibuffer + 1
                              fdata(3,imesh,jmesh,isub) =   &
     &                           fdata(3,imesh,jmesh,isub)+ &
     &                           INT(ppm_recvbufferd(ibuffer))
                              ibuffer = ibuffer + 1
                              fdata(4,imesh,jmesh,isub) =   &
     &                           fdata(4,imesh,jmesh,isub)+ &
     &                           INT(ppm_recvbufferd(ibuffer))
                              ibuffer = ibuffer + 1
                              fdata(5,imesh,jmesh,isub) =   &
     &                           fdata(5,imesh,jmesh,isub)+ &
     &                           INT(ppm_recvbufferd(ibuffer))
#elif  __KIND == __LOGICAL
                              IF (ppm_recvbufferd(ibuffer) .GT.     &
     &                           (1.0_ppm_kind_double-ppm_myepsd)) THEN
                                 fdata(1,imesh,jmesh,isub) =    &
     &                              (fdata(1,imesh,jmesh,isub) .AND. .TRUE.)
                              ELSE
                                 fdata(1,imesh,jmesh,isub) =    &
     &                              (fdata(1,imesh,jmesh,isub) .AND. .FALSE.)
                              ENDIF
                              ibuffer = ibuffer + 1
                              IF (ppm_recvbufferd(ibuffer) .GT.     &
     &                           (1.0_ppm_kind_double-ppm_myepsd)) THEN
                                 fdata(2,imesh,jmesh,isub) =    &
     &                              (fdata(2,imesh,jmesh,isub) .AND. .TRUE.)
                              ELSE
                                 fdata(2,imesh,jmesh,isub) =    &
     &                              (fdata(2,imesh,jmesh,isub) .AND. .FALSE.)
                              ENDIF
                              ibuffer = ibuffer + 1
                              IF (ppm_recvbufferd(ibuffer) .GT.     &
     &                           (1.0_ppm_kind_double-ppm_myepsd)) THEN
                                 fdata(3,imesh,jmesh,isub) =    &
     &                              (fdata(3,imesh,jmesh,isub) .AND. .TRUE.)
                              ELSE
                                 fdata(3,imesh,jmesh,isub) =    &
     &                              (fdata(3,imesh,jmesh,isub) .AND. .FALSE.)
                              ENDIF
                              ibuffer = ibuffer + 1
                              IF (ppm_recvbufferd(ibuffer) .GT.     &
     &                           (1.0_ppm_kind_double-ppm_myepsd)) THEN
                                 fdata(4,imesh,jmesh,isub) =    &
     &                              (fdata(4,imesh,jmesh,isub) .AND. .TRUE.)
                              ELSE
                                 fdata(4,imesh,jmesh,isub) =    &
     &                              (fdata(4,imesh,jmesh,isub) .AND. .FALSE.)
                              ENDIF
                              ibuffer = ibuffer + 1
                              IF (ppm_recvbufferd(ibuffer) .GT.     &
     &                           (1.0_ppm_kind_double-ppm_myepsd)) THEN
                                 fdata(5,imesh,jmesh,isub) =    &
     &                              (fdata(5,imesh,jmesh,isub) .AND. .TRUE.)
                              ELSE
                                 fdata(5,imesh,jmesh,isub) =    &
     &                              (fdata(5,imesh,jmesh,isub) .AND. .FALSE.)
                              ENDIF
#endif
                           ENDDO
                        ENDDO
                     ENDIF
                  !-------------------------------------------------------------
                  !  For edim.GT.5 the vector length will be edim !!
                  !-------------------------------------------------------------
                  ELSE
                     IF (rtype .EQ. ppm_param_pop_replace) THEN
                        DO jmesh=ylo,yhi
                           DO imesh=xlo,xhi
                              DO k=1,edim
#if    __KIND == __SINGLE_PRECISION_COMPLEX | __KIND == __DOUBLE_PRECISION_COMPLEX
                                 ibuffer = ibuffer + 2
#else
                                 ibuffer = ibuffer + 1
#endif
#if    __KIND == __SINGLE_PRECISION
                                 fdata(k,imesh,jmesh,isub) =     &
     &                              REAL(ppm_recvbufferd(ibuffer),ppm_kind_single)
#elif  __KIND == __DOUBLE_PRECISION
                                 fdata(k,imesh,jmesh,isub) =     &
     &                              ppm_recvbufferd(ibuffer)
#elif  __KIND == __SINGLE_PRECISION_COMPLEX
                                 fdata(k,imesh,jmesh,isub) =     &
     &                              CMPLX(ppm_recvbufferd(ibuffer-1),  &
     &                              ppm_recvbufferd(ibuffer),ppm_kind_single)
#elif  __KIND == __DOUBLE_PRECISION_COMPLEX
                                 fdata(k,imesh,jmesh,isub) =     &
     &                              CMPLX(ppm_recvbufferd(ibuffer-1),  &
     &                              ppm_recvbufferd(ibuffer),ppm_kind_double)
#elif  __KIND == __INTEGER
                                 fdata(k,imesh,jmesh,isub) =     &
     &                              INT(ppm_recvbufferd(ibuffer))
#elif  __KIND == __LOGICAL
                                 IF (ppm_recvbufferd(ibuffer) .GT.     &
     &                              (1.0_ppm_kind_double-ppm_myepsd)) THEN
                                    fdata(k,imesh,jmesh,isub) = .TRUE.
                                 ELSE
                                    fdata(k,imesh,jmesh,isub) = .FALSE.
                                 ENDIF
#endif
                              ENDDO
                           ENDDO
                        ENDDO
                     ELSEIF (rtype .EQ. ppm_param_pop_add) THEN
                        DO jmesh=ylo,yhi
                           DO imesh=xlo,xhi
                              DO k=1,edim
#if    __KIND == __SINGLE_PRECISION_COMPLEX | __KIND == __DOUBLE_PRECISION_COMPLEX
                                 ibuffer = ibuffer + 2
#else
                                 ibuffer = ibuffer + 1
#endif
#if    __KIND == __SINGLE_PRECISION
                                 fdata(k,imesh,jmesh,isub) =   &
     &                              fdata(k,imesh,jmesh,isub)+ &
     &                              REAL(ppm_recvbufferd(ibuffer),ppm_kind_single)
#elif  __KIND == __DOUBLE_PRECISION
                                 fdata(k,imesh,jmesh,isub) =   &
     &                              fdata(k,imesh,jmesh,isub)+ &
     &                              ppm_recvbufferd(ibuffer)
#elif  __KIND == __SINGLE_PRECISION_COMPLEX
                                 fdata(k,imesh,jmesh,isub) =   &
     &                              fdata(k,imesh,jmesh,isub)+ &
     &                              CMPLX(ppm_recvbufferd(ibuffer-1),  &
     &                              ppm_recvbufferd(ibuffer),ppm_kind_single)
#elif  __KIND == __DOUBLE_PRECISION_COMPLEX
                                 fdata(k,imesh,jmesh,isub) =   &
     &                              fdata(k,imesh,jmesh,isub)+ &
     &                              CMPLX(ppm_recvbufferd(ibuffer-1),  &
     &                              ppm_recvbufferd(ibuffer),ppm_kind_double)
#elif  __KIND == __INTEGER
                                 fdata(k,imesh,jmesh,isub) =   &
     &                              fdata(k,imesh,jmesh,isub)+ &
     &                              INT(ppm_recvbufferd(ibuffer))
#elif  __KIND == __LOGICAL
                                 IF (ppm_recvbufferd(ibuffer) .GT.     &
     &                              (1.0_ppm_kind_double-ppm_myepsd)) THEN
                                    fdata(k,imesh,jmesh,isub) =    &
     &                                 (fdata(k,imesh,jmesh,isub) .AND. .TRUE.)
                                 ELSE
                                    fdata(k,imesh,jmesh,isub) =    &
     &                                 (fdata(k,imesh,jmesh,isub) .AND. .FALSE.)
                                 ENDIF
#endif
                              ENDDO
                           ENDDO
                        ENDDO
                     ENDIF         ! rtype
                  ENDIF            ! lda = ...
               !----------------------------------------------------------------
               !  With mask: This will not vectorize so no need to unroll
               !----------------------------------------------------------------
               ELSE
                  DO jmesh=ylo,yhi
                     DO imesh=xlo,xhi
                        IF (mask(imesh,jmesh,isub)) THEN
                           DO k=1,edim
#if    __KIND == __SINGLE_PRECISION_COMPLEX | __KIND == __DOUBLE_PRECISION_COMPLEX
                              ibuffer = ibuffer + 2
#else
                              ibuffer = ibuffer + 1
#endif
                              IF (rtype .EQ. ppm_param_pop_replace) THEN
#if    __KIND == __SINGLE_PRECISION
                                 fdata(k,imesh,jmesh,isub) =     &
     &                              REAL(ppm_recvbufferd(ibuffer),ppm_kind_single)
#elif  __KIND == __DOUBLE_PRECISION
                                 fdata(k,imesh,jmesh,isub) =     &
     &                              ppm_recvbufferd(ibuffer)
#elif  __KIND == __SINGLE_PRECISION_COMPLEX
                                 fdata(k,imesh,jmesh,isub) =     &
     &                              CMPLX(ppm_recvbufferd(ibuffer-1),  &
     &                              ppm_recvbufferd(ibuffer),ppm_kind_single)
#elif  __KIND == __DOUBLE_PRECISION_COMPLEX
                                 fdata(k,imesh,jmesh,isub) =     &
     &                              CMPLX(ppm_recvbufferd(ibuffer-1),  &
     &                              ppm_recvbufferd(ibuffer),ppm_kind_double)
#elif  __KIND == __INTEGER
                                 fdata(k,imesh,jmesh,isub) =     &
     &                              INT(ppm_recvbufferd(ibuffer))
#elif  __KIND == __LOGICAL
                                 IF (ppm_recvbufferd(ibuffer) .GT.     &
     &                              (1.0_ppm_kind_double-ppm_myepsd)) THEN
                                    fdata(k,imesh,jmesh,isub) = .TRUE.
                                 ELSE
                                    fdata(k,imesh,jmesh,isub) = .FALSE.
                                 ENDIF
#endif
                              ELSEIF (rtype .EQ. ppm_param_pop_add) THEN
#if    __KIND == __SINGLE_PRECISION
                                 fdata(k,imesh,jmesh,isub) =   &
     &                              fdata(k,imesh,jmesh,isub)+ &
     &                              REAL(ppm_recvbufferd(ibuffer),ppm_kind_single)
#elif  __KIND == __DOUBLE_PRECISION
                                 fdata(k,imesh,jmesh,isub) =   &
     &                              fdata(k,imesh,jmesh,isub)+ &
     &                              ppm_recvbufferd(ibuffer)
#elif  __KIND == __SINGLE_PRECISION_COMPLEX
                                 fdata(k,imesh,jmesh,isub) =   &
     &                              fdata(k,imesh,jmesh,isub)+ &
     &                              CMPLX(ppm_recvbufferd(ibuffer-1),  &
     &                              ppm_recvbufferd(ibuffer),ppm_kind_single)
#elif  __KIND == __DOUBLE_PRECISION_COMPLEX
                                 fdata(k,imesh,jmesh,isub) =   &
     &                              fdata(k,imesh,jmesh,isub)+ &
     &                              CMPLX(ppm_recvbufferd(ibuffer-1),  &
     &                              ppm_recvbufferd(ibuffer),ppm_kind_double)
#elif  __KIND == __INTEGER
                                 fdata(k,imesh,jmesh,isub) =   &
     &                              fdata(k,imesh,jmesh,isub)+ &
     &                              INT(ppm_recvbufferd(ibuffer))
#elif  __KIND == __LOGICAL
                                 IF (ppm_recvbufferd(ibuffer) .GT.     &
     &                              (1.0_ppm_kind_double-ppm_myepsd)) THEN
                                    fdata(k,imesh,jmesh,isub) =    &
     &                                 (fdata(k,imesh,jmesh,isub) .AND. .TRUE.)
                                 ELSE
                                    fdata(k,imesh,jmesh,isub) =    &
     &                                 (fdata(k,imesh,jmesh,isub) .AND. .FALSE.)
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
               ENDIF            ! PRESENT(mask)
#elif __DIM == __SFIELD
               IF (.NOT.PRESENT(mask)) THEN
                  IF (rtype .EQ. ppm_param_pop_replace) THEN
                     DO jmesh=ylo,yhi
                        DO imesh=xlo,xhi
#if    __KIND == __SINGLE_PRECISION_COMPLEX | __KIND == __DOUBLE_PRECISION_COMPLEX
                           ibuffer = ibuffer + 2
#else
                           ibuffer = ibuffer + 1
#endif
#if    __KIND == __SINGLE_PRECISION
                           fdata(imesh,jmesh,isub) =     &
     &                        REAL(ppm_recvbufferd(ibuffer),ppm_kind_single)
#elif  __KIND == __DOUBLE_PRECISION
                           fdata(imesh,jmesh,isub) = ppm_recvbufferd(ibuffer)
#elif  __KIND == __SINGLE_PRECISION_COMPLEX
                           fdata(imesh,jmesh,isub) =     &
     &                        CMPLX(ppm_recvbufferd(ibuffer-1),  &
     &                        ppm_recvbufferd(ibuffer),ppm_kind_single)
#elif  __KIND == __DOUBLE_PRECISION_COMPLEX
                           fdata(imesh,jmesh,isub) =     &
     &                        CMPLX(ppm_recvbufferd(ibuffer-1),  &
     &                        ppm_recvbufferd(ibuffer),ppm_kind_double)
#elif  __KIND == __INTEGER
                           fdata(imesh,jmesh,isub) =     &
     &                        INT(ppm_recvbufferd(ibuffer))
#elif  __KIND == __LOGICAL
                           IF (ppm_recvbufferd(ibuffer) .GT.     &
     &                        (1.0_ppm_kind_double-ppm_myepsd)) THEN
                              fdata(imesh,jmesh,isub) = .TRUE.
                           ELSE
                              fdata(imesh,jmesh,isub) = .FALSE.
                           ENDIF
#endif
                        ENDDO
                     ENDDO
                  ELSEIF (rtype .EQ. ppm_param_pop_add) THEN
                     DO jmesh=ylo,yhi
                        DO imesh=xlo,xhi
#if    __KIND == __SINGLE_PRECISION_COMPLEX | __KIND == __DOUBLE_PRECISION_COMPLEX
                           ibuffer = ibuffer + 2
#else
                           ibuffer = ibuffer + 1
#endif
#if    __KIND == __SINGLE_PRECISION
                           fdata(imesh,jmesh,isub) = fdata(imesh,jmesh,isub)+&
     &                        REAL(ppm_recvbufferd(ibuffer),ppm_kind_single)
#elif  __KIND == __DOUBLE_PRECISION
                           fdata(imesh,jmesh,isub) = fdata(imesh,jmesh,isub)+&
     &                         ppm_recvbufferd(ibuffer)
#elif  __KIND == __SINGLE_PRECISION_COMPLEX
                           fdata(imesh,jmesh,isub) = fdata(imesh,jmesh,isub)+&
     &                        CMPLX(ppm_recvbufferd(ibuffer-1),  &
     &                        ppm_recvbufferd(ibuffer),ppm_kind_single)
#elif  __KIND == __DOUBLE_PRECISION_COMPLEX
                           fdata(imesh,jmesh,isub) = fdata(imesh,jmesh,isub)+&
     &                        CMPLX(ppm_recvbufferd(ibuffer-1),  &
     &                        ppm_recvbufferd(ibuffer),ppm_kind_double)
#elif  __KIND == __INTEGER
                           fdata(imesh,jmesh,isub) = fdata(imesh,jmesh,isub)+ &
     &                        INT(ppm_recvbufferd(ibuffer))
#elif  __KIND == __LOGICAL
                           IF (ppm_recvbufferd(ibuffer) .GT.     &
     &                        (1.0_ppm_kind_double-ppm_myepsd)) THEN
                              fdata(imesh,jmesh,isub) =    &
     &                           (fdata(imesh,jmesh,isub) .AND. .TRUE.)
                           ELSE
                              fdata(imesh,jmesh,isub) =    &
     &                           (fdata(imesh,jmesh,isub) .AND. .FALSE.)
                           ENDIF
#endif
                        ENDDO
                     ENDDO
                  ENDIF
               ELSE             ! PRESENT(mask)
                  DO jmesh=ylo,yhi
                     DO imesh=xlo,xhi
                        IF (mask(imesh,jmesh,isub)) THEN
#if    __KIND == __SINGLE_PRECISION_COMPLEX | __KIND == __DOUBLE_PRECISION_COMPLEX
                           ibuffer = ibuffer + 2
#else
                           ibuffer = ibuffer + 1
#endif
                           IF (rtype .EQ. ppm_param_pop_replace) THEN
#if    __KIND == __SINGLE_PRECISION
                              fdata(imesh,jmesh,isub) =     &
     &                           REAL(ppm_recvbufferd(ibuffer),ppm_kind_single)
#elif  __KIND == __DOUBLE_PRECISION
                              fdata(imesh,jmesh,isub) = ppm_recvbufferd(ibuffer)
#elif  __KIND == __SINGLE_PRECISION_COMPLEX
                              fdata(imesh,jmesh,isub) =     &
     &                           CMPLX(ppm_recvbufferd(ibuffer-1),  &
     &                           ppm_recvbufferd(ibuffer),ppm_kind_single)
#elif  __KIND == __DOUBLE_PRECISION_COMPLEX
                              fdata(imesh,jmesh,isub) =     &
     &                           CMPLX(ppm_recvbufferd(ibuffer-1),  &
     &                           ppm_recvbufferd(ibuffer),ppm_kind_double)
#elif  __KIND == __INTEGER
                              fdata(imesh,jmesh,isub) =     &
     &                           INT(ppm_recvbufferd(ibuffer))
#elif  __KIND == __LOGICAL
                              IF (ppm_recvbufferd(ibuffer) .GT.     &
     &                           (1.0_ppm_kind_double-ppm_myepsd)) THEN
                                 fdata(imesh,jmesh,isub) = .TRUE.
                              ELSE
                                 fdata(imesh,jmesh,isub) = .FALSE.
                              ENDIF
#endif
                           ELSEIF (rtype .EQ. ppm_param_pop_add) THEN
#if    __KIND == __SINGLE_PRECISION
                              fdata(imesh,jmesh,isub) = fdata(imesh,jmesh,isub)+&
     &                           REAL(ppm_recvbufferd(ibuffer),ppm_kind_single)
#elif  __KIND == __DOUBLE_PRECISION
                              fdata(imesh,jmesh,isub) = fdata(imesh,jmesh,isub)+&
     &                            ppm_recvbufferd(ibuffer)
#elif  __KIND == __SINGLE_PRECISION_COMPLEX
                              fdata(imesh,jmesh,isub) = fdata(imesh,jmesh,isub)+&
     &                           CMPLX(ppm_recvbufferd(ibuffer-1),  &
     &                           ppm_recvbufferd(ibuffer),ppm_kind_single)
#elif  __KIND == __DOUBLE_PRECISION_COMPLEX
                              fdata(imesh,jmesh,isub) = fdata(imesh,jmesh,isub)+&
     &                           CMPLX(ppm_recvbufferd(ibuffer-1),  &
     &                           ppm_recvbufferd(ibuffer),ppm_kind_double)
#elif  __KIND == __INTEGER
                              fdata(imesh,jmesh,isub) = fdata(imesh,jmesh,isub)+ &
     &                           INT(ppm_recvbufferd(ibuffer))
#elif  __KIND == __LOGICAL
                              IF (ppm_recvbufferd(ibuffer) .GT.     &
     &                           (1.0_ppm_kind_double-ppm_myepsd)) THEN
                                 fdata(imesh,jmesh,isub) =    &
     &                              (fdata(imesh,jmesh,isub) .AND. .TRUE.)
                              ELSE
                                 fdata(imesh,jmesh,isub) =    &
     &                              (fdata(imesh,jmesh,isub) .AND. .FALSE.)
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
               mofs(1) = target_mesh%istart(1,jsub)-1
               mofs(2) = target_mesh%istart(2,jsub)-1
               !----------------------------------------------------------------
               !  Get boundaries of mesh block to be received in local sub
               !  coordinates
               !----------------------------------------------------------------
               xlo = ppm_mesh_irecvblkstart(1,j)-mofs(1)
               ylo = ppm_mesh_irecvblkstart(2,j)-mofs(2)
               xhi = xlo+ppm_mesh_irecvblksize(1,j)-1
               yhi = ylo+ppm_mesh_irecvblksize(2,j)-1
               IF (ppm_debug .GT. 1) THEN
                   WRITE(mesg,'(A,2I4)') 'start: ',             &
     &                 ppm_mesh_irecvblkstart(1,j),ppm_mesh_irecvblkstart(2,j)
                   CALL ppm_write(ppm_rank,'ppm_map_field_pop_2d',mesg,info)
                   WRITE(mesg,'(A,2I4)') 'size: ',             &
     &                 ppm_mesh_irecvblksize(1,j),ppm_mesh_irecvblksize(2,j)
                   CALL ppm_write(ppm_rank,'ppm_map_field_pop_2d',mesg,info)
                   WRITE(mesg,'(A,2I4)') 'mesh offset: ',mofs(1),mofs(2)
                   CALL ppm_write(ppm_rank,'ppm_map_field_pop_2d',mesg,info)
                   WRITE(mesg,'(A,2I4)') 'xlo, xhi: ',xlo,xhi
                   CALL ppm_write(ppm_rank,'ppm_map_field_pop_2d',mesg,info)
                   WRITE(mesg,'(A,2I4)') 'ylo, yhi: ',ylo,yhi
                   CALL ppm_write(ppm_rank,'ppm_map_field_pop_2d',mesg,info)
                   WRITE(mesg,'(A,I1)') 'buffer dim: ',edim
                   CALL ppm_write(ppm_rank,'ppm_map_field_pop_2d',mesg,info)
#if   __DIM == __VFIELD
                   WRITE(mesg,'(A,2I4)') 'SIZE(fdata): ',SIZE(fdata,2),  &
     &                 SIZE(fdata,3)
                   CALL ppm_write(ppm_rank,'ppm_map_field_pop_2d',mesg,info)
#elif __DIM == __SFIELD
                   WRITE(mesg,'(A,2I4)') 'SIZE(fdata): ',SIZE(fdata,1),  &
     &                 SIZE(fdata,2)
                   CALL ppm_write(ppm_rank,'ppm_map_field_pop_2d',mesg,info)
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
                        DO jmesh=ylo,yhi
                           DO imesh=xlo,xhi
#if    __KIND == __SINGLE_PRECISION_COMPLEX | __KIND == __DOUBLE_PRECISION_COMPLEX
                              ibuffer = ibuffer + 2
#else
                              ibuffer = ibuffer + 1
#endif
#if    __KIND == __DOUBLE_PRECISION
                              fdata(1,imesh,jmesh,isub) =     &
     &                           REAL(ppm_recvbuffers(ibuffer),ppm_kind_double)
#elif  __KIND == __SINGLE_PRECISION
                              fdata(1,imesh,jmesh,isub) =     &
     &                           ppm_recvbuffers(ibuffer)
#elif  __KIND == __SINGLE_PRECISION_COMPLEX
                              fdata(1,imesh,jmesh,isub) =     &
     &                           CMPLX(ppm_recvbuffers(ibuffer-1),  &
     &                           ppm_recvbuffers(ibuffer),ppm_kind_single)
#elif  __KIND == __DOUBLE_PRECISION_COMPLEX
                              fdata(1,imesh,jmesh,isub) =     &
     &                           CMPLX(ppm_recvbuffers(ibuffer-1),  &
     &                           ppm_recvbuffers(ibuffer),ppm_kind_double)
#elif  __KIND == __INTEGER
                              fdata(1,imesh,jmesh,isub) =     &
     &                           INT(ppm_recvbuffers(ibuffer))
#elif  __KIND == __LOGICAL
                              IF (ppm_recvbuffers(ibuffer) .GT.     &
     &                           (1.0_ppm_kind_single-ppm_myepss)) THEN
                                 fdata(1,imesh,jmesh,isub) = .TRUE.
                              ELSE
                                 fdata(1,imesh,jmesh,isub) = .FALSE.
                              ENDIF
#endif
                           ENDDO
                        ENDDO
                     ELSEIF (rtype .EQ. ppm_param_pop_add) THEN
                        DO jmesh=ylo,yhi
                           DO imesh=xlo,xhi
#if    __KIND == __SINGLE_PRECISION_COMPLEX | __KIND == __DOUBLE_PRECISION_COMPLEX
                              ibuffer = ibuffer + 2
#else
                              ibuffer = ibuffer + 1
#endif
#if    __KIND == __DOUBLE_PRECISION
                              fdata(1,imesh,jmesh,isub) =    &
     &                           fdata(1,imesh,jmesh,isub) + &
     &                           REAL(ppm_recvbuffers(ibuffer),ppm_kind_double)
#elif  __KIND == __SINGLE_PRECISION
                              fdata(1,imesh,jmesh,isub) =    &
     &                           fdata(1,imesh,jmesh,isub) + &
     &                           ppm_recvbuffers(ibuffer)
#elif  __KIND == __SINGLE_PRECISION_COMPLEX
                              fdata(1,imesh,jmesh,isub) =    &
     &                           fdata(1,imesh,jmesh,isub) + &
     &                           CMPLX(ppm_recvbuffers(ibuffer-1),  &
     &                           ppm_recvbuffers(ibuffer),ppm_kind_single)
#elif  __KIND == __DOUBLE_PRECISION_COMPLEX
                              fdata(1,imesh,jmesh,isub) =    &
     &                           fdata(1,imesh,jmesh,isub) + &
     &                           CMPLX(ppm_recvbuffers(ibuffer-1),  &
     &                           ppm_recvbuffers(ibuffer),ppm_kind_double)
#elif  __KIND == __INTEGER
                              fdata(1,imesh,jmesh,isub) =    &
     &                           fdata(1,imesh,jmesh,isub) + &
     &                           INT(ppm_recvbuffers(ibuffer))
#elif  __KIND == __LOGICAL
                              IF (ppm_recvbuffers(ibuffer) .GT.     &
     &                           (1.0_ppm_kind_single-ppm_myepss)) THEN
                                 fdata(1,imesh,jmesh,isub) =        &
     &                               (fdata(1,imesh,jmesh,isub) .AND. .TRUE.)
                              ELSE
                                 fdata(1,imesh,jmesh,isub) =        &
     &                               (fdata(1,imesh,jmesh,isub) .AND. .FALSE.)
                              ENDIF
#endif
                           ENDDO
                        ENDDO
                     ENDIF
                  !-------------------------------------------------------------
                  !  Unrolled for edim=2
                  !-------------------------------------------------------------
                  ELSEIF (edim .EQ. 2) THEN
                     IF (rtype .EQ. ppm_param_pop_replace) THEN
                        DO jmesh=ylo,yhi
                           DO imesh=xlo,xhi
#if    __KIND == __SINGLE_PRECISION_COMPLEX | __KIND == __DOUBLE_PRECISION_COMPLEX
                              ibuffer = ibuffer + 2
#else
                              ibuffer = ibuffer + 1
#endif
#if    __KIND == __DOUBLE_PRECISION
                              fdata(1,imesh,jmesh,isub) =     &
     &                           REAL(ppm_recvbuffers(ibuffer),ppm_kind_double)
                              ibuffer = ibuffer + 1
                              fdata(2,imesh,jmesh,isub) =     &
     &                           REAL(ppm_recvbuffers(ibuffer),ppm_kind_double)
#elif  __KIND == __SINGLE_PRECISION
                              fdata(1,imesh,jmesh,isub) =     &
     &                           ppm_recvbuffers(ibuffer)
                              ibuffer = ibuffer + 1
                              fdata(2,imesh,jmesh,isub) =     &
     &                           ppm_recvbuffers(ibuffer)
#elif  __KIND == __SINGLE_PRECISION_COMPLEX
                              fdata(1,imesh,jmesh,isub) =     &
     &                           CMPLX(ppm_recvbuffers(ibuffer-1),  &
     &                           ppm_recvbuffers(ibuffer),ppm_kind_single)
                              ibuffer = ibuffer + 2
                              fdata(2,imesh,jmesh,isub) =     &
     &                           CMPLX(ppm_recvbuffers(ibuffer-1),  &
     &                           ppm_recvbuffers(ibuffer),ppm_kind_single)
#elif  __KIND == __DOUBLE_PRECISION_COMPLEX
                              fdata(1,imesh,jmesh,isub) =     &
     &                           CMPLX(ppm_recvbuffers(ibuffer-1),  &
     &                           ppm_recvbuffers(ibuffer),ppm_kind_double)
                              ibuffer = ibuffer + 2
                              fdata(2,imesh,jmesh,isub) =     &
     &                           CMPLX(ppm_recvbuffers(ibuffer-1),  &
     &                           ppm_recvbuffers(ibuffer),ppm_kind_double)
#elif  __KIND == __INTEGER
                              fdata(1,imesh,jmesh,isub) =     &
     &                           INT(ppm_recvbuffers(ibuffer))
                              ibuffer = ibuffer + 1
                              fdata(2,imesh,jmesh,isub) =     &
     &                           INT(ppm_recvbuffers(ibuffer))
#elif  __KIND == __LOGICAL
                              IF (ppm_recvbuffers(ibuffer) .GT.     &
     &                           (1.0_ppm_kind_single-ppm_myepss)) THEN
                                 fdata(1,imesh,jmesh,isub) = .TRUE.
                              ELSE
                                 fdata(1,imesh,jmesh,isub) = .FALSE.
                              ENDIF
                              ibuffer = ibuffer + 1
                              IF (ppm_recvbuffers(ibuffer) .GT.     &
     &                           (1.0_ppm_kind_single-ppm_myepss)) THEN
                                 fdata(2,imesh,jmesh,isub) = .TRUE.
                              ELSE
                                 fdata(2,imesh,jmesh,isub) = .FALSE.
                              ENDIF
#endif
                           ENDDO
                        ENDDO
                     ELSEIF (rtype .EQ. ppm_param_pop_add) THEN
                        DO jmesh=ylo,yhi
                           DO imesh=xlo,xhi
#if    __KIND == __SINGLE_PRECISION_COMPLEX | __KIND == __DOUBLE_PRECISION_COMPLEX
                              ibuffer = ibuffer + 2
#else
                              ibuffer = ibuffer + 1
#endif
#if    __KIND == __DOUBLE_PRECISION
                              fdata(1,imesh,jmesh,isub) =    &
     &                           fdata(1,imesh,jmesh,isub) + &
     &                           REAL(ppm_recvbuffers(ibuffer),ppm_kind_double)
                              ibuffer = ibuffer + 1
                              fdata(2,imesh,jmesh,isub) =    &
     &                           fdata(2,imesh,jmesh,isub) + &
     &                           REAL(ppm_recvbuffers(ibuffer),ppm_kind_double)
#elif  __KIND == __SINGLE_PRECISION
                              fdata(1,imesh,jmesh,isub) =    &
     &                           fdata(1,imesh,jmesh,isub) + &
     &                           ppm_recvbuffers(ibuffer)
                              ibuffer = ibuffer + 1
                              fdata(2,imesh,jmesh,isub) =    &
     &                           fdata(2,imesh,jmesh,isub) + &
     &                           ppm_recvbuffers(ibuffer)
#elif  __KIND == __SINGLE_PRECISION_COMPLEX
                              fdata(1,imesh,jmesh,isub) =    &
     &                           fdata(1,imesh,jmesh,isub) + &
     &                           CMPLX(ppm_recvbuffers(ibuffer-1),  &
     &                           ppm_recvbuffers(ibuffer),ppm_kind_single)
                              ibuffer = ibuffer + 2
                              fdata(2,imesh,jmesh,isub) =    &
     &                           fdata(2,imesh,jmesh,isub) + &
     &                           CMPLX(ppm_recvbuffers(ibuffer-1),  &
     &                           ppm_recvbuffers(ibuffer),ppm_kind_single)
#elif  __KIND == __DOUBLE_PRECISION_COMPLEX
                              fdata(1,imesh,jmesh,isub) =    &
     &                           fdata(1,imesh,jmesh,isub) + &
     &                           CMPLX(ppm_recvbuffers(ibuffer-1),  &
     &                           ppm_recvbuffers(ibuffer),ppm_kind_double)
                              ibuffer = ibuffer + 2
                              fdata(2,imesh,jmesh,isub) =    &
     &                           fdata(2,imesh,jmesh,isub) + &
     &                           CMPLX(ppm_recvbuffers(ibuffer-1),  &
     &                           ppm_recvbuffers(ibuffer),ppm_kind_double)
#elif  __KIND == __INTEGER
                              fdata(1,imesh,jmesh,isub) =    &
     &                           fdata(1,imesh,jmesh,isub) + &
     &                           INT(ppm_recvbuffers(ibuffer))
                              ibuffer = ibuffer + 1
                              fdata(2,imesh,jmesh,isub) =    &
     &                           fdata(2,imesh,jmesh,isub) + &
     &                           INT(ppm_recvbuffers(ibuffer))
#elif  __KIND == __LOGICAL
                              IF (ppm_recvbuffers(ibuffer) .GT.     &
     &                           (1.0_ppm_kind_single-ppm_myepss)) THEN
                                 fdata(1,imesh,jmesh,isub) =        &
     &                               (fdata(1,imesh,jmesh,isub) .AND. .TRUE.)
                              ELSE
                                 fdata(1,imesh,jmesh,isub) =        &
     &                               (fdata(1,imesh,jmesh,isub) .AND. .FALSE.)
                              ENDIF
                              ibuffer = ibuffer + 1
                              IF (ppm_recvbuffers(ibuffer) .GT.     &
     &                           (1.0_ppm_kind_single-ppm_myepss)) THEN
                                 fdata(2,imesh,jmesh,isub) =        &
     &                               (fdata(2,imesh,jmesh,isub) .AND. .TRUE.)
                              ELSE
                                 fdata(2,imesh,jmesh,isub) =        &
     &                               (fdata(2,imesh,jmesh,isub) .AND. .FALSE.)
                              ENDIF
#endif
                           ENDDO
                        ENDDO
                     ENDIF
                  !-------------------------------------------------------------
                  !  Unrolled for edim=3
                  !-------------------------------------------------------------
                  ELSEIF (edim .EQ. 3) THEN
                     IF (rtype .EQ. ppm_param_pop_replace) THEN
                        DO jmesh=ylo,yhi
                           DO imesh=xlo,xhi
#if    __KIND == __SINGLE_PRECISION_COMPLEX | __KIND == __DOUBLE_PRECISION_COMPLEX
                              ibuffer = ibuffer + 2
#else
                              ibuffer = ibuffer + 1
#endif
#if    __KIND == __DOUBLE_PRECISION
                              fdata(1,imesh,jmesh,isub) =     &
     &                           REAL(ppm_recvbuffers(ibuffer),ppm_kind_double)
                              ibuffer = ibuffer + 1
                              fdata(2,imesh,jmesh,isub) =     &
     &                           REAL(ppm_recvbuffers(ibuffer),ppm_kind_double)
                              ibuffer = ibuffer + 1
                              fdata(3,imesh,jmesh,isub) =     &
     &                           REAL(ppm_recvbuffers(ibuffer),ppm_kind_double)
#elif  __KIND == __SINGLE_PRECISION
                              fdata(1,imesh,jmesh,isub) =     &
     &                           ppm_recvbuffers(ibuffer)
                              ibuffer = ibuffer + 1
                              fdata(2,imesh,jmesh,isub) =     &
     &                           ppm_recvbuffers(ibuffer)
                              ibuffer = ibuffer + 1
                              fdata(3,imesh,jmesh,isub) =     &
     &                           ppm_recvbuffers(ibuffer)
#elif  __KIND == __SINGLE_PRECISION_COMPLEX
                              fdata(1,imesh,jmesh,isub) =     &
     &                           CMPLX(ppm_recvbuffers(ibuffer-1),  &
     &                           ppm_recvbuffers(ibuffer),ppm_kind_single)
                              ibuffer = ibuffer + 2
                              fdata(2,imesh,jmesh,isub) =     &
     &                           CMPLX(ppm_recvbuffers(ibuffer-1),  &
     &                           ppm_recvbuffers(ibuffer),ppm_kind_single)
                              ibuffer = ibuffer + 2
                              fdata(3,imesh,jmesh,isub) =     &
     &                           CMPLX(ppm_recvbuffers(ibuffer-1),  &
     &                           ppm_recvbuffers(ibuffer),ppm_kind_single)
#elif  __KIND == __DOUBLE_PRECISION_COMPLEX
                              fdata(1,imesh,jmesh,isub) =     &
     &                           CMPLX(ppm_recvbuffers(ibuffer-1),  &
     &                           ppm_recvbuffers(ibuffer),ppm_kind_double)
                              ibuffer = ibuffer + 2
                              fdata(2,imesh,jmesh,isub) =     &
     &                           CMPLX(ppm_recvbuffers(ibuffer-1),  &
     &                           ppm_recvbuffers(ibuffer),ppm_kind_double)
                              ibuffer = ibuffer + 2
                              fdata(3,imesh,jmesh,isub) =     &
     &                           CMPLX(ppm_recvbuffers(ibuffer-1),  &
     &                           ppm_recvbuffers(ibuffer),ppm_kind_double)
#elif  __KIND == __INTEGER
                              fdata(1,imesh,jmesh,isub) =     &
     &                           INT(ppm_recvbuffers(ibuffer))
                              ibuffer = ibuffer + 1
                              fdata(2,imesh,jmesh,isub) =     &
     &                           INT(ppm_recvbuffers(ibuffer))
                              ibuffer = ibuffer + 1
                              fdata(3,imesh,jmesh,isub) =     &
     &                           INT(ppm_recvbuffers(ibuffer))
#elif  __KIND == __LOGICAL
                              IF (ppm_recvbuffers(ibuffer) .GT.     &
     &                           (1.0_ppm_kind_single-ppm_myepss)) THEN
                                 fdata(1,imesh,jmesh,isub) = .TRUE.
                              ELSE
                                 fdata(1,imesh,jmesh,isub) = .FALSE.
                              ENDIF
                              ibuffer = ibuffer + 1
                              IF (ppm_recvbuffers(ibuffer) .GT.     &
     &                           (1.0_ppm_kind_single-ppm_myepss)) THEN
                                 fdata(2,imesh,jmesh,isub) = .TRUE.
                              ELSE
                                 fdata(2,imesh,jmesh,isub) = .FALSE.
                              ENDIF
                              ibuffer = ibuffer + 1
                              IF (ppm_recvbuffers(ibuffer) .GT.     &
     &                           (1.0_ppm_kind_single-ppm_myepss)) THEN
                                 fdata(3,imesh,jmesh,isub) = .TRUE.
                              ELSE
                                 fdata(3,imesh,jmesh,isub) = .FALSE.
                              ENDIF
#endif
                           ENDDO
                        ENDDO
                     ELSEIF (rtype .EQ. ppm_param_pop_add) THEN
                        DO jmesh=ylo,yhi
                           DO imesh=xlo,xhi
#if    __KIND == __SINGLE_PRECISION_COMPLEX | __KIND == __DOUBLE_PRECISION_COMPLEX
                              ibuffer = ibuffer + 2
#else
                              ibuffer = ibuffer + 1
#endif
#if    __KIND == __DOUBLE_PRECISION
                              fdata(1,imesh,jmesh,isub) =    &
     &                           fdata(1,imesh,jmesh,isub) + &
     &                           REAL(ppm_recvbuffers(ibuffer),ppm_kind_double)
                              ibuffer = ibuffer + 1
                              fdata(2,imesh,jmesh,isub) =    &
     &                           fdata(2,imesh,jmesh,isub) + &
     &                           REAL(ppm_recvbuffers(ibuffer),ppm_kind_double)
                              ibuffer = ibuffer + 1
                              fdata(3,imesh,jmesh,isub) =    &
     &                           fdata(3,imesh,jmesh,isub) + &
     &                           REAL(ppm_recvbuffers(ibuffer),ppm_kind_double)
#elif  __KIND == __SINGLE_PRECISION
                              fdata(1,imesh,jmesh,isub) =    &
     &                           fdata(1,imesh,jmesh,isub) + &
     &                           ppm_recvbuffers(ibuffer)
                              ibuffer = ibuffer + 1
                              fdata(2,imesh,jmesh,isub) =    &
     &                           fdata(2,imesh,jmesh,isub) + &
     &                           ppm_recvbuffers(ibuffer)
                              ibuffer = ibuffer + 1
                              fdata(3,imesh,jmesh,isub) =    &
     &                           fdata(3,imesh,jmesh,isub) + &
     &                           ppm_recvbuffers(ibuffer)
#elif  __KIND == __SINGLE_PRECISION_COMPLEX
                              fdata(1,imesh,jmesh,isub) =    &
     &                           fdata(1,imesh,jmesh,isub) + &
     &                           CMPLX(ppm_recvbuffers(ibuffer-1),  &
     &                           ppm_recvbuffers(ibuffer),ppm_kind_single)
                              ibuffer = ibuffer + 2
                              fdata(2,imesh,jmesh,isub) =    &
     &                           fdata(2,imesh,jmesh,isub) + &
     &                           CMPLX(ppm_recvbuffers(ibuffer-1),  &
     &                           ppm_recvbuffers(ibuffer),ppm_kind_single)
                              ibuffer = ibuffer + 2
                              fdata(3,imesh,jmesh,isub) =    &
     &                           fdata(3,imesh,jmesh,isub) + &
     &                           CMPLX(ppm_recvbuffers(ibuffer-1),  &
     &                           ppm_recvbuffers(ibuffer),ppm_kind_single)
#elif  __KIND == __DOUBLE_PRECISION_COMPLEX
                              fdata(1,imesh,jmesh,isub) =    &
     &                           fdata(1,imesh,jmesh,isub) + &
     &                           CMPLX(ppm_recvbuffers(ibuffer-1),  &
     &                           ppm_recvbuffers(ibuffer),ppm_kind_double)
                              ibuffer = ibuffer + 2
                              fdata(2,imesh,jmesh,isub) =    &
     &                           fdata(2,imesh,jmesh,isub) + &
     &                           CMPLX(ppm_recvbuffers(ibuffer-1),  &
     &                           ppm_recvbuffers(ibuffer),ppm_kind_double)
                              ibuffer = ibuffer + 2
                              fdata(3,imesh,jmesh,isub) =    &
     &                           fdata(3,imesh,jmesh,isub) + &
     &                           CMPLX(ppm_recvbuffers(ibuffer-1),  &
     &                           ppm_recvbuffers(ibuffer),ppm_kind_double)
#elif  __KIND == __INTEGER
                              fdata(1,imesh,jmesh,isub) =    &
     &                           fdata(1,imesh,jmesh,isub) + &
     &                           INT(ppm_recvbuffers(ibuffer))
                              ibuffer = ibuffer + 1
                              fdata(2,imesh,jmesh,isub) =    &
     &                           fdata(2,imesh,jmesh,isub) + &
     &                           INT(ppm_recvbuffers(ibuffer))
                              ibuffer = ibuffer + 1
                              fdata(3,imesh,jmesh,isub) =    &
     &                           fdata(3,imesh,jmesh,isub) + &
     &                           INT(ppm_recvbuffers(ibuffer))
#elif  __KIND == __LOGICAL
                              IF (ppm_recvbuffers(ibuffer) .GT.     &
     &                           (1.0_ppm_kind_single-ppm_myepss)) THEN
                                 fdata(1,imesh,jmesh,isub) =        &
     &                               (fdata(1,imesh,jmesh,isub) .AND. .TRUE.)
                              ELSE
                                 fdata(1,imesh,jmesh,isub) =        &
     &                               (fdata(1,imesh,jmesh,isub) .AND. .FALSE.)
                              ENDIF
                              ibuffer = ibuffer + 1
                              IF (ppm_recvbuffers(ibuffer) .GT.     &
     &                           (1.0_ppm_kind_single-ppm_myepss)) THEN
                                 fdata(2,imesh,jmesh,isub) =        &
     &                               (fdata(2,imesh,jmesh,isub) .AND. .TRUE.)
                              ELSE
                                 fdata(2,imesh,jmesh,isub) =        &
     &                               (fdata(2,imesh,jmesh,isub) .AND. .FALSE.)
                              ENDIF
                              ibuffer = ibuffer + 1
                              IF (ppm_recvbuffers(ibuffer) .GT.     &
     &                           (1.0_ppm_kind_single-ppm_myepss)) THEN
                                 fdata(3,imesh,jmesh,isub) =        &
     &                               (fdata(3,imesh,jmesh,isub) .AND. .TRUE.)
                              ELSE
                                 fdata(3,imesh,jmesh,isub) =        &
     &                               (fdata(3,imesh,jmesh,isub) .AND. .FALSE.)
                              ENDIF
#endif
                           ENDDO
                        ENDDO
                     ENDIF
                  !-------------------------------------------------------------
                  !  Unrolled for edim=4
                  !-------------------------------------------------------------
                  ELSEIF (edim .EQ. 4) THEN
                     IF (rtype .EQ. ppm_param_pop_replace) THEN
                        DO jmesh=ylo,yhi
                           DO imesh=xlo,xhi
#if    __KIND == __SINGLE_PRECISION_COMPLEX | __KIND == __DOUBLE_PRECISION_COMPLEX
                              ibuffer = ibuffer + 2
#else
                              ibuffer = ibuffer + 1
#endif
#if    __KIND == __DOUBLE_PRECISION
                              fdata(1,imesh,jmesh,isub) =     &
     &                           REAL(ppm_recvbuffers(ibuffer),ppm_kind_double)
                              ibuffer = ibuffer + 1
                              fdata(2,imesh,jmesh,isub) =     &
     &                           REAL(ppm_recvbuffers(ibuffer),ppm_kind_double)
                              ibuffer = ibuffer + 1
                              fdata(3,imesh,jmesh,isub) =     &
     &                           REAL(ppm_recvbuffers(ibuffer),ppm_kind_double)
                              ibuffer = ibuffer + 1
                              fdata(4,imesh,jmesh,isub) =     &
     &                           REAL(ppm_recvbuffers(ibuffer),ppm_kind_double)
#elif  __KIND == __SINGLE_PRECISION
                              fdata(1,imesh,jmesh,isub) =     &
     &                           ppm_recvbuffers(ibuffer)
                              ibuffer = ibuffer + 1
                              fdata(2,imesh,jmesh,isub) =     &
     &                           ppm_recvbuffers(ibuffer)
                              ibuffer = ibuffer + 1
                              fdata(3,imesh,jmesh,isub) =     &
     &                           ppm_recvbuffers(ibuffer)
                              ibuffer = ibuffer + 1
                              fdata(4,imesh,jmesh,isub) =     &
     &                           ppm_recvbuffers(ibuffer)
#elif  __KIND == __SINGLE_PRECISION_COMPLEX
                              fdata(1,imesh,jmesh,isub) =     &
     &                           CMPLX(ppm_recvbuffers(ibuffer-1),  &
     &                           ppm_recvbuffers(ibuffer),ppm_kind_single)
                              ibuffer = ibuffer + 2
                              fdata(2,imesh,jmesh,isub) =     &
     &                           CMPLX(ppm_recvbuffers(ibuffer-1),  &
     &                           ppm_recvbuffers(ibuffer),ppm_kind_single)
                              ibuffer = ibuffer + 2
                              fdata(3,imesh,jmesh,isub) =     &
     &                           CMPLX(ppm_recvbuffers(ibuffer-1),  &
     &                           ppm_recvbuffers(ibuffer),ppm_kind_single)
                              ibuffer = ibuffer + 2
                              fdata(4,imesh,jmesh,isub) =     &
     &                           CMPLX(ppm_recvbuffers(ibuffer-1),  &
     &                           ppm_recvbuffers(ibuffer),ppm_kind_single)
#elif  __KIND == __DOUBLE_PRECISION_COMPLEX
                              fdata(1,imesh,jmesh,isub) =     &
     &                           CMPLX(ppm_recvbuffers(ibuffer-1),  &
     &                           ppm_recvbuffers(ibuffer),ppm_kind_double)
                              ibuffer = ibuffer + 2
                              fdata(2,imesh,jmesh,isub) =     &
     &                           CMPLX(ppm_recvbuffers(ibuffer-1),  &
     &                           ppm_recvbuffers(ibuffer),ppm_kind_double)
                              ibuffer = ibuffer + 2
                              fdata(3,imesh,jmesh,isub) =     &
     &                           CMPLX(ppm_recvbuffers(ibuffer-1),  &
     &                           ppm_recvbuffers(ibuffer),ppm_kind_double)
                              ibuffer = ibuffer + 2
                              fdata(4,imesh,jmesh,isub) =     &
     &                           CMPLX(ppm_recvbuffers(ibuffer-1),  &
     &                           ppm_recvbuffers(ibuffer),ppm_kind_double)
#elif  __KIND == __INTEGER
                              fdata(1,imesh,jmesh,isub) =     &
     &                           INT(ppm_recvbuffers(ibuffer))
                              ibuffer = ibuffer + 1
                              fdata(2,imesh,jmesh,isub) =     &
     &                           INT(ppm_recvbuffers(ibuffer))
                              ibuffer = ibuffer + 1
                              fdata(3,imesh,jmesh,isub) =     &
     &                           INT(ppm_recvbuffers(ibuffer))
                              ibuffer = ibuffer + 1
                              fdata(4,imesh,jmesh,isub) =     &
     &                           INT(ppm_recvbuffers(ibuffer))
#elif  __KIND == __LOGICAL
                              IF (ppm_recvbuffers(ibuffer) .GT.     &
     &                           (1.0_ppm_kind_single-ppm_myepss)) THEN
                                 fdata(1,imesh,jmesh,isub) = .TRUE.
                              ELSE
                                 fdata(1,imesh,jmesh,isub) = .FALSE.
                              ENDIF
                              ibuffer = ibuffer + 1
                              IF (ppm_recvbuffers(ibuffer) .GT.     &
     &                           (1.0_ppm_kind_single-ppm_myepss)) THEN
                                 fdata(2,imesh,jmesh,isub) = .TRUE.
                              ELSE
                                 fdata(2,imesh,jmesh,isub) = .FALSE.
                              ENDIF
                              ibuffer = ibuffer + 1
                              IF (ppm_recvbuffers(ibuffer) .GT.     &
     &                           (1.0_ppm_kind_single-ppm_myepss)) THEN
                                 fdata(3,imesh,jmesh,isub) = .TRUE.
                              ELSE
                                 fdata(3,imesh,jmesh,isub) = .FALSE.
                              ENDIF
                              ibuffer = ibuffer + 1
                              IF (ppm_recvbuffers(ibuffer) .GT.     &
     &                           (1.0_ppm_kind_single-ppm_myepss)) THEN
                                 fdata(4,imesh,jmesh,isub) = .TRUE.
                              ELSE
                                 fdata(4,imesh,jmesh,isub) = .FALSE.
                              ENDIF
#endif
                           ENDDO
                        ENDDO
                     ELSEIF (rtype .EQ. ppm_param_pop_add) THEN
                        DO jmesh=ylo,yhi
                           DO imesh=xlo,xhi
#if    __KIND == __SINGLE_PRECISION_COMPLEX | __KIND == __DOUBLE_PRECISION_COMPLEX
                              ibuffer = ibuffer + 2
#else
                              ibuffer = ibuffer + 1
#endif
#if    __KIND == __DOUBLE_PRECISION
                              fdata(1,imesh,jmesh,isub) =    &
     &                           fdata(1,imesh,jmesh,isub) + &
     &                           REAL(ppm_recvbuffers(ibuffer),ppm_kind_double)
                              ibuffer = ibuffer + 1
                              fdata(2,imesh,jmesh,isub) =    &
     &                           fdata(2,imesh,jmesh,isub) + &
     &                           REAL(ppm_recvbuffers(ibuffer),ppm_kind_double)
                              ibuffer = ibuffer + 1
                              fdata(3,imesh,jmesh,isub) =    &
     &                           fdata(3,imesh,jmesh,isub) + &
     &                           REAL(ppm_recvbuffers(ibuffer),ppm_kind_double)
                              ibuffer = ibuffer + 1
                              fdata(4,imesh,jmesh,isub) =    &
     &                           fdata(4,imesh,jmesh,isub) + &
     &                           REAL(ppm_recvbuffers(ibuffer),ppm_kind_double)
#elif  __KIND == __SINGLE_PRECISION
                              fdata(1,imesh,jmesh,isub) =    &
     &                           fdata(1,imesh,jmesh,isub) + &
     &                           ppm_recvbuffers(ibuffer)
                              ibuffer = ibuffer + 1
                              fdata(2,imesh,jmesh,isub) =    &
     &                           fdata(2,imesh,jmesh,isub) + &
     &                           ppm_recvbuffers(ibuffer)
                              ibuffer = ibuffer + 1
                              fdata(3,imesh,jmesh,isub) =    &
     &                           fdata(3,imesh,jmesh,isub) + &
     &                           ppm_recvbuffers(ibuffer)
                              ibuffer = ibuffer + 1
                              fdata(4,imesh,jmesh,isub) =    &
     &                           fdata(4,imesh,jmesh,isub) + &
     &                           ppm_recvbuffers(ibuffer)
#elif  __KIND == __SINGLE_PRECISION_COMPLEX
                              fdata(1,imesh,jmesh,isub) =    &
     &                           fdata(1,imesh,jmesh,isub) + &
     &                           CMPLX(ppm_recvbuffers(ibuffer-1),  &
     &                           ppm_recvbuffers(ibuffer),ppm_kind_single)
                              ibuffer = ibuffer + 2
                              fdata(2,imesh,jmesh,isub) =    &
     &                           fdata(2,imesh,jmesh,isub) + &
     &                           CMPLX(ppm_recvbuffers(ibuffer-1),  &
     &                           ppm_recvbuffers(ibuffer),ppm_kind_single)
                              ibuffer = ibuffer + 2
                              fdata(3,imesh,jmesh,isub) =    &
     &                           fdata(3,imesh,jmesh,isub) + &
     &                           CMPLX(ppm_recvbuffers(ibuffer-1),  &
     &                           ppm_recvbuffers(ibuffer),ppm_kind_single)
                              ibuffer = ibuffer + 2
                              fdata(4,imesh,jmesh,isub) =    &
     &                           fdata(4,imesh,jmesh,isub) + &
     &                           CMPLX(ppm_recvbuffers(ibuffer-1),  &
     &                           ppm_recvbuffers(ibuffer),ppm_kind_single)
#elif  __KIND == __DOUBLE_PRECISION_COMPLEX
                              fdata(1,imesh,jmesh,isub) =    &
     &                           fdata(1,imesh,jmesh,isub) + &
     &                           CMPLX(ppm_recvbuffers(ibuffer-1),  &
     &                           ppm_recvbuffers(ibuffer),ppm_kind_double)
                              ibuffer = ibuffer + 2
                              fdata(2,imesh,jmesh,isub) =    &
     &                           fdata(2,imesh,jmesh,isub) + &
     &                           CMPLX(ppm_recvbuffers(ibuffer-1),  &
     &                           ppm_recvbuffers(ibuffer),ppm_kind_double)
                              ibuffer = ibuffer + 2
                              fdata(3,imesh,jmesh,isub) =    &
     &                           fdata(3,imesh,jmesh,isub) + &
     &                           CMPLX(ppm_recvbuffers(ibuffer-1),  &
     &                           ppm_recvbuffers(ibuffer),ppm_kind_double)
                              ibuffer = ibuffer + 2
                              fdata(4,imesh,jmesh,isub) =    &
     &                           fdata(4,imesh,jmesh,isub) + &
     &                           CMPLX(ppm_recvbuffers(ibuffer-1),  &
     &                           ppm_recvbuffers(ibuffer),ppm_kind_double)
#elif  __KIND == __INTEGER
                              fdata(1,imesh,jmesh,isub) =    &
     &                           fdata(1,imesh,jmesh,isub) + &
     &                           INT(ppm_recvbuffers(ibuffer))
                              ibuffer = ibuffer + 1
                              fdata(2,imesh,jmesh,isub) =    &
     &                           fdata(2,imesh,jmesh,isub) + &
     &                           INT(ppm_recvbuffers(ibuffer))
                              ibuffer = ibuffer + 1
                              fdata(3,imesh,jmesh,isub) =    &
     &                           fdata(3,imesh,jmesh,isub) + &
     &                           INT(ppm_recvbuffers(ibuffer))
                              ibuffer = ibuffer + 1
                              fdata(4,imesh,jmesh,isub) =    &
     &                           fdata(4,imesh,jmesh,isub) + &
     &                           INT(ppm_recvbuffers(ibuffer))
#elif  __KIND == __LOGICAL
                              IF (ppm_recvbuffers(ibuffer) .GT.     &
     &                           (1.0_ppm_kind_single-ppm_myepss)) THEN
                                 fdata(1,imesh,jmesh,isub) =        &
     &                               (fdata(1,imesh,jmesh,isub) .AND. .TRUE.)
                              ELSE
                                 fdata(1,imesh,jmesh,isub) =        &
     &                               (fdata(1,imesh,jmesh,isub) .AND. .FALSE.)
                              ENDIF
                              ibuffer = ibuffer + 1
                              IF (ppm_recvbuffers(ibuffer) .GT.     &
     &                           (1.0_ppm_kind_single-ppm_myepss)) THEN
                                 fdata(2,imesh,jmesh,isub) =        &
     &                               (fdata(2,imesh,jmesh,isub) .AND. .TRUE.)
                              ELSE
                                 fdata(2,imesh,jmesh,isub) =        &
     &                               (fdata(2,imesh,jmesh,isub) .AND. .FALSE.)
                              ENDIF
                              ibuffer = ibuffer + 1
                              IF (ppm_recvbuffers(ibuffer) .GT.     &
     &                           (1.0_ppm_kind_single-ppm_myepss)) THEN
                                 fdata(3,imesh,jmesh,isub) =        &
     &                               (fdata(3,imesh,jmesh,isub) .AND. .TRUE.)
                              ELSE
                                 fdata(3,imesh,jmesh,isub) =        &
     &                               (fdata(3,imesh,jmesh,isub) .AND. .FALSE.)
                              ENDIF
                              ibuffer = ibuffer + 1
                              IF (ppm_recvbuffers(ibuffer) .GT.     &
     &                           (1.0_ppm_kind_single-ppm_myepss)) THEN
                                 fdata(4,imesh,jmesh,isub) =        &
     &                               (fdata(4,imesh,jmesh,isub) .AND. .TRUE.)
                              ELSE
                                 fdata(4,imesh,jmesh,isub) =        &
     &                               (fdata(4,imesh,jmesh,isub) .AND. .FALSE.)
                              ENDIF
#endif
                           ENDDO
                        ENDDO
                     ENDIF
                  !-------------------------------------------------------------
                  !  Unrolled for edim=5
                  !-------------------------------------------------------------
                  ELSEIF (edim .EQ. 5) THEN
                     IF (rtype .EQ. ppm_param_pop_replace) THEN
                        DO jmesh=ylo,yhi
                           DO imesh=xlo,xhi
#if    __KIND == __SINGLE_PRECISION_COMPLEX | __KIND == __DOUBLE_PRECISION_COMPLEX
                              ibuffer = ibuffer + 2
#else
                              ibuffer = ibuffer + 1
#endif
#if    __KIND == __DOUBLE_PRECISION
                              fdata(1,imesh,jmesh,isub) =     &
     &                           REAL(ppm_recvbuffers(ibuffer),ppm_kind_double)
                              ibuffer = ibuffer + 1
                              fdata(2,imesh,jmesh,isub) =     &
     &                           REAL(ppm_recvbuffers(ibuffer),ppm_kind_double)
                              ibuffer = ibuffer + 1
                              fdata(3,imesh,jmesh,isub) =     &
     &                           REAL(ppm_recvbuffers(ibuffer),ppm_kind_double)
                              ibuffer = ibuffer + 1
                              fdata(4,imesh,jmesh,isub) =     &
     &                           REAL(ppm_recvbuffers(ibuffer),ppm_kind_double)
                              ibuffer = ibuffer + 1
                              fdata(5,imesh,jmesh,isub) =     &
     &                           REAL(ppm_recvbuffers(ibuffer),ppm_kind_double)
#elif  __KIND == __SINGLE_PRECISION
                              fdata(1,imesh,jmesh,isub) =     &
     &                           ppm_recvbuffers(ibuffer)
                              ibuffer = ibuffer + 1
                              fdata(2,imesh,jmesh,isub) =     &
     &                           ppm_recvbuffers(ibuffer)
                              ibuffer = ibuffer + 1
                              fdata(3,imesh,jmesh,isub) =     &
     &                           ppm_recvbuffers(ibuffer)
                              ibuffer = ibuffer + 1
                              fdata(4,imesh,jmesh,isub) =     &
     &                           ppm_recvbuffers(ibuffer)
                              ibuffer = ibuffer + 1
                              fdata(5,imesh,jmesh,isub) =     &
     &                           ppm_recvbuffers(ibuffer)
#elif  __KIND == __SINGLE_PRECISION_COMPLEX
                              fdata(1,imesh,jmesh,isub) =     &
     &                           CMPLX(ppm_recvbuffers(ibuffer-1),  &
     &                           ppm_recvbuffers(ibuffer),ppm_kind_single)
                              ibuffer = ibuffer + 2
                              fdata(2,imesh,jmesh,isub) =     &
     &                           CMPLX(ppm_recvbuffers(ibuffer-1),  &
     &                           ppm_recvbuffers(ibuffer),ppm_kind_single)
                              ibuffer = ibuffer + 2
                              fdata(3,imesh,jmesh,isub) =     &
     &                           CMPLX(ppm_recvbuffers(ibuffer-1),  &
     &                           ppm_recvbuffers(ibuffer),ppm_kind_single)
                              ibuffer = ibuffer + 2
                              fdata(4,imesh,jmesh,isub) =     &
     &                           CMPLX(ppm_recvbuffers(ibuffer-1),  &
     &                           ppm_recvbuffers(ibuffer),ppm_kind_single)
                              ibuffer = ibuffer + 2
                              fdata(5,imesh,jmesh,isub) =     &
     &                           CMPLX(ppm_recvbuffers(ibuffer-1),  &
     &                           ppm_recvbuffers(ibuffer),ppm_kind_single)
#elif  __KIND == __DOUBLE_PRECISION_COMPLEX
                              fdata(1,imesh,jmesh,isub) =     &
     &                           CMPLX(ppm_recvbuffers(ibuffer-1),  &
     &                           ppm_recvbuffers(ibuffer),ppm_kind_double)
                              ibuffer = ibuffer + 2
                              fdata(2,imesh,jmesh,isub) =     &
     &                           CMPLX(ppm_recvbuffers(ibuffer-1),  &
     &                           ppm_recvbuffers(ibuffer),ppm_kind_double)
                              ibuffer = ibuffer + 2
                              fdata(3,imesh,jmesh,isub) =     &
     &                           CMPLX(ppm_recvbuffers(ibuffer-1),  &
     &                           ppm_recvbuffers(ibuffer),ppm_kind_double)
                              ibuffer = ibuffer + 2
                              fdata(4,imesh,jmesh,isub) =     &
     &                           CMPLX(ppm_recvbuffers(ibuffer-1),  &
     &                           ppm_recvbuffers(ibuffer),ppm_kind_double)
                              ibuffer = ibuffer + 2
                              fdata(5,imesh,jmesh,isub) =     &
     &                           CMPLX(ppm_recvbuffers(ibuffer-1),  &
     &                           ppm_recvbuffers(ibuffer),ppm_kind_double)
#elif  __KIND == __INTEGER
                              fdata(1,imesh,jmesh,isub) =     &
     &                           INT(ppm_recvbuffers(ibuffer))
                              ibuffer = ibuffer + 1
                              fdata(2,imesh,jmesh,isub) =     &
     &                           INT(ppm_recvbuffers(ibuffer))
                              ibuffer = ibuffer + 1
                              fdata(3,imesh,jmesh,isub) =     &
     &                           INT(ppm_recvbuffers(ibuffer))
                              ibuffer = ibuffer + 1
                              fdata(4,imesh,jmesh,isub) =     &
     &                           INT(ppm_recvbuffers(ibuffer))
                              ibuffer = ibuffer + 1
                              fdata(5,imesh,jmesh,isub) =     &
     &                           INT(ppm_recvbuffers(ibuffer))
#elif  __KIND == __LOGICAL
                              IF (ppm_recvbuffers(ibuffer) .GT.     &
     &                           (1.0_ppm_kind_single-ppm_myepss)) THEN
                                 fdata(1,imesh,jmesh,isub) = .TRUE.
                              ELSE
                                 fdata(1,imesh,jmesh,isub) = .FALSE.
                              ENDIF
                              ibuffer = ibuffer + 1
                              IF (ppm_recvbuffers(ibuffer) .GT.     &
     &                           (1.0_ppm_kind_single-ppm_myepss)) THEN
                                 fdata(2,imesh,jmesh,isub) = .TRUE.
                              ELSE
                                 fdata(2,imesh,jmesh,isub) = .FALSE.
                              ENDIF
                              ibuffer = ibuffer + 1
                              IF (ppm_recvbuffers(ibuffer) .GT.     &
     &                           (1.0_ppm_kind_single-ppm_myepss)) THEN
                                 fdata(3,imesh,jmesh,isub) = .TRUE.
                              ELSE
                                 fdata(3,imesh,jmesh,isub) = .FALSE.
                              ENDIF
                              ibuffer = ibuffer + 1
                              IF (ppm_recvbuffers(ibuffer) .GT.     &
     &                           (1.0_ppm_kind_single-ppm_myepss)) THEN
                                 fdata(4,imesh,jmesh,isub) = .TRUE.
                              ELSE
                                 fdata(4,imesh,jmesh,isub) = .FALSE.
                              ENDIF
                              ibuffer = ibuffer + 1
                              IF (ppm_recvbuffers(ibuffer) .GT.     &
     &                           (1.0_ppm_kind_single-ppm_myepss)) THEN
                                 fdata(5,imesh,jmesh,isub) = .TRUE.
                              ELSE
                                 fdata(5,imesh,jmesh,isub) = .FALSE.
                              ENDIF
#endif
                           ENDDO
                        ENDDO
                     ELSEIF (rtype .EQ. ppm_param_pop_add) THEN
                        DO jmesh=ylo,yhi
                           DO imesh=xlo,xhi
#if    __KIND == __SINGLE_PRECISION_COMPLEX | __KIND == __DOUBLE_PRECISION_COMPLEX
                              ibuffer = ibuffer + 2
#else
                              ibuffer = ibuffer + 1
#endif
#if    __KIND == __DOUBLE_PRECISION
                              fdata(1,imesh,jmesh,isub) =    &
     &                           fdata(1,imesh,jmesh,isub) + &
     &                           REAL(ppm_recvbuffers(ibuffer),ppm_kind_double)
                              ibuffer = ibuffer + 1
                              fdata(2,imesh,jmesh,isub) =    &
     &                           fdata(2,imesh,jmesh,isub) + &
     &                           REAL(ppm_recvbuffers(ibuffer),ppm_kind_double)
                              ibuffer = ibuffer + 1
                              fdata(3,imesh,jmesh,isub) =    &
     &                           fdata(3,imesh,jmesh,isub) + &
     &                           REAL(ppm_recvbuffers(ibuffer),ppm_kind_double)
                              ibuffer = ibuffer + 1
                              fdata(4,imesh,jmesh,isub) =    &
     &                           fdata(4,imesh,jmesh,isub) + &
     &                           REAL(ppm_recvbuffers(ibuffer),ppm_kind_double)
                              ibuffer = ibuffer + 1
                              fdata(5,imesh,jmesh,isub) =    &
     &                           fdata(5,imesh,jmesh,isub) + &
     &                           REAL(ppm_recvbuffers(ibuffer),ppm_kind_double)
#elif  __KIND == __SINGLE_PRECISION
                              fdata(1,imesh,jmesh,isub) =    &
     &                           fdata(1,imesh,jmesh,isub) + &
     &                           ppm_recvbuffers(ibuffer)
                              ibuffer = ibuffer + 1
                              fdata(2,imesh,jmesh,isub) =    &
     &                           fdata(2,imesh,jmesh,isub) + &
     &                           ppm_recvbuffers(ibuffer)
                              ibuffer = ibuffer + 1
                              fdata(3,imesh,jmesh,isub) =    &
     &                           fdata(3,imesh,jmesh,isub) + &
     &                           ppm_recvbuffers(ibuffer)
                              ibuffer = ibuffer + 1
                              fdata(4,imesh,jmesh,isub) =    &
     &                           fdata(4,imesh,jmesh,isub) + &
     &                           ppm_recvbuffers(ibuffer)
                              ibuffer = ibuffer + 1
                              fdata(5,imesh,jmesh,isub) =    &
     &                           fdata(5,imesh,jmesh,isub) + &
     &                           ppm_recvbuffers(ibuffer)
#elif  __KIND == __SINGLE_PRECISION_COMPLEX
                              fdata(1,imesh,jmesh,isub) =    &
     &                           fdata(1,imesh,jmesh,isub) + &
     &                           CMPLX(ppm_recvbuffers(ibuffer-1),  &
     &                           ppm_recvbuffers(ibuffer),ppm_kind_single)
                              ibuffer = ibuffer + 2
                              fdata(2,imesh,jmesh,isub) =    &
     &                           fdata(2,imesh,jmesh,isub) + &
     &                           CMPLX(ppm_recvbuffers(ibuffer-1),  &
     &                           ppm_recvbuffers(ibuffer),ppm_kind_single)
                              ibuffer = ibuffer + 2
                              fdata(3,imesh,jmesh,isub) =    &
     &                           fdata(3,imesh,jmesh,isub) + &
     &                           CMPLX(ppm_recvbuffers(ibuffer-1),  &
     &                           ppm_recvbuffers(ibuffer),ppm_kind_single)
                              ibuffer = ibuffer + 2
                              fdata(4,imesh,jmesh,isub) =    &
     &                           fdata(4,imesh,jmesh,isub) + &
     &                           CMPLX(ppm_recvbuffers(ibuffer-1),  &
     &                           ppm_recvbuffers(ibuffer),ppm_kind_single)
                              ibuffer = ibuffer + 2
                              fdata(5,imesh,jmesh,isub) =    &
     &                           fdata(5,imesh,jmesh,isub) + &
     &                           CMPLX(ppm_recvbuffers(ibuffer-1),  &
     &                           ppm_recvbuffers(ibuffer),ppm_kind_single)
#elif  __KIND == __DOUBLE_PRECISION_COMPLEX
                              fdata(1,imesh,jmesh,isub) =    &
     &                           fdata(1,imesh,jmesh,isub) + &
     &                           CMPLX(ppm_recvbuffers(ibuffer-1),  &
     &                           ppm_recvbuffers(ibuffer),ppm_kind_double)
                              ibuffer = ibuffer + 2
                              fdata(2,imesh,jmesh,isub) =    &
     &                           fdata(2,imesh,jmesh,isub) + &
     &                           CMPLX(ppm_recvbuffers(ibuffer-1),  &
     &                           ppm_recvbuffers(ibuffer),ppm_kind_double)
                              ibuffer = ibuffer + 2
                              fdata(3,imesh,jmesh,isub) =    &
     &                           fdata(3,imesh,jmesh,isub) + &
     &                           CMPLX(ppm_recvbuffers(ibuffer-1),  &
     &                           ppm_recvbuffers(ibuffer),ppm_kind_double)
                              ibuffer = ibuffer + 2
                              fdata(4,imesh,jmesh,isub) =    &
     &                           fdata(4,imesh,jmesh,isub) + &
     &                           CMPLX(ppm_recvbuffers(ibuffer-1),  &
     &                           ppm_recvbuffers(ibuffer),ppm_kind_double)
                              ibuffer = ibuffer + 2
                              fdata(5,imesh,jmesh,isub) =    &
     &                           fdata(5,imesh,jmesh,isub) + &
     &                           CMPLX(ppm_recvbuffers(ibuffer-1),  &
     &                           ppm_recvbuffers(ibuffer),ppm_kind_double)
#elif  __KIND == __INTEGER
                              fdata(1,imesh,jmesh,isub) =    &
     &                           fdata(1,imesh,jmesh,isub) + &
     &                           INT(ppm_recvbuffers(ibuffer))
                              ibuffer = ibuffer + 1
                              fdata(2,imesh,jmesh,isub) =    &
     &                           fdata(2,imesh,jmesh,isub) + &
     &                           INT(ppm_recvbuffers(ibuffer))
                              ibuffer = ibuffer + 1
                              fdata(3,imesh,jmesh,isub) =    &
     &                           fdata(3,imesh,jmesh,isub) + &
     &                           INT(ppm_recvbuffers(ibuffer))
                              ibuffer = ibuffer + 1
                              fdata(4,imesh,jmesh,isub) =    &
     &                           fdata(4,imesh,jmesh,isub) + &
     &                           INT(ppm_recvbuffers(ibuffer))
                              ibuffer = ibuffer + 1
                              fdata(5,imesh,jmesh,isub) =    &
     &                           fdata(5,imesh,jmesh,isub) + &
     &                           INT(ppm_recvbuffers(ibuffer))
#elif  __KIND == __LOGICAL
                              IF (ppm_recvbuffers(ibuffer) .GT.     &
     &                           (1.0_ppm_kind_single-ppm_myepss)) THEN
                                 fdata(1,imesh,jmesh,isub) =        &
     &                               (fdata(1,imesh,jmesh,isub) .AND. .TRUE.)
                              ELSE
                                 fdata(1,imesh,jmesh,isub) =        &
     &                               (fdata(1,imesh,jmesh,isub) .AND. .FALSE.)
                              ENDIF
                              ibuffer = ibuffer + 1
                              IF (ppm_recvbuffers(ibuffer) .GT.     &
     &                           (1.0_ppm_kind_single-ppm_myepss)) THEN
                                 fdata(2,imesh,jmesh,isub) =        &
     &                               (fdata(2,imesh,jmesh,isub) .AND. .TRUE.)
                              ELSE
                                 fdata(2,imesh,jmesh,isub) =        &
     &                               (fdata(2,imesh,jmesh,isub) .AND. .FALSE.)
                              ENDIF
                              ibuffer = ibuffer + 1
                              IF (ppm_recvbuffers(ibuffer) .GT.     &
     &                           (1.0_ppm_kind_single-ppm_myepss)) THEN
                                 fdata(3,imesh,jmesh,isub) =        &
     &                               (fdata(3,imesh,jmesh,isub) .AND. .TRUE.)
                              ELSE
                                 fdata(3,imesh,jmesh,isub) =        &
     &                               (fdata(3,imesh,jmesh,isub) .AND. .FALSE.)
                              ENDIF
                              ibuffer = ibuffer + 1
                              IF (ppm_recvbuffers(ibuffer) .GT.     &
     &                           (1.0_ppm_kind_single-ppm_myepss)) THEN
                                 fdata(4,imesh,jmesh,isub) =        &
     &                               (fdata(4,imesh,jmesh,isub) .AND. .TRUE.)
                              ELSE
                                 fdata(4,imesh,jmesh,isub) =        &
     &                               (fdata(4,imesh,jmesh,isub) .AND. .FALSE.)
                              ENDIF
                              ibuffer = ibuffer + 1
                              IF (ppm_recvbuffers(ibuffer) .GT.     &
     &                           (1.0_ppm_kind_single-ppm_myepss)) THEN
                                 fdata(5,imesh,jmesh,isub) =        &
     &                               (fdata(5,imesh,jmesh,isub) .AND. .TRUE.)
                              ELSE
                                 fdata(5,imesh,jmesh,isub) =        &
     &                               (fdata(5,imesh,jmesh,isub) .AND. .FALSE.)
                              ENDIF
#endif
                           ENDDO
                        ENDDO
                     ENDIF
                  !-------------------------------------------------------------
                  !  For edim.GT.5 the vector length will be edim !!
                  !-------------------------------------------------------------
                  ELSE
                     IF (rtype .EQ. ppm_param_pop_replace) THEN
                        DO jmesh=ylo,yhi
                           DO imesh=xlo,xhi
                              DO k=1,edim
#if    __KIND == __SINGLE_PRECISION_COMPLEX | __KIND == __DOUBLE_PRECISION_COMPLEX
                                 ibuffer = ibuffer + 2
#else
                                 ibuffer = ibuffer + 1
#endif
#if    __KIND == __DOUBLE_PRECISION
                                 fdata(k,imesh,jmesh,isub) =     &
     &                              REAL(ppm_recvbuffers(ibuffer),ppm_kind_double)
#elif  __KIND == __SINGLE_PRECISION
                                 fdata(k,imesh,jmesh,isub) =     &
     &                              ppm_recvbuffers(ibuffer)
#elif  __KIND == __SINGLE_PRECISION_COMPLEX
                                 fdata(k,imesh,jmesh,isub) =     &
     &                              CMPLX(ppm_recvbuffers(ibuffer-1),  &
     &                              ppm_recvbuffers(ibuffer),ppm_kind_single)
#elif  __KIND == __DOUBLE_PRECISION_COMPLEX
                                 fdata(k,imesh,jmesh,isub) =     &
     &                              CMPLX(ppm_recvbuffers(ibuffer-1),  &
     &                              ppm_recvbuffers(ibuffer),ppm_kind_double)
#elif  __KIND == __INTEGER
                                 fdata(k,imesh,jmesh,isub) =     &
     &                              INT(ppm_recvbuffers(ibuffer))
#elif  __KIND == __LOGICAL
                                 IF (ppm_recvbuffers(ibuffer) .GT.     &
     &                              (1.0_ppm_kind_single-ppm_myepss)) THEN
                                    fdata(k,imesh,jmesh,isub) = .TRUE.
                                 ELSE
                                    fdata(k,imesh,jmesh,isub) = .FALSE.
                                 ENDIF
#endif
                              ENDDO
                           ENDDO
                        ENDDO
                     ELSEIF (rtype .EQ. ppm_param_pop_add) THEN
                        DO jmesh=ylo,yhi
                           DO imesh=xlo,xhi
                              DO k=1,edim
#if    __KIND == __SINGLE_PRECISION_COMPLEX | __KIND == __DOUBLE_PRECISION_COMPLEX
                                 ibuffer = ibuffer + 2
#else
                                 ibuffer = ibuffer + 1
#endif
#if    __KIND == __DOUBLE_PRECISION
                                 fdata(k,imesh,jmesh,isub) =    &
     &                              fdata(k,imesh,jmesh,isub) + &
     &                              REAL(ppm_recvbuffers(ibuffer),ppm_kind_double)
#elif  __KIND == __SINGLE_PRECISION
                                 fdata(k,imesh,jmesh,isub) =    &
     &                              fdata(k,imesh,jmesh,isub) + &
     &                              ppm_recvbuffers(ibuffer)
#elif  __KIND == __SINGLE_PRECISION_COMPLEX
                                 fdata(k,imesh,jmesh,isub) =    &
     &                              fdata(k,imesh,jmesh,isub) + &
     &                              CMPLX(ppm_recvbuffers(ibuffer-1),  &
     &                              ppm_recvbuffers(ibuffer),ppm_kind_single)
#elif  __KIND == __DOUBLE_PRECISION_COMPLEX
                                 fdata(k,imesh,jmesh,isub) =    &
     &                              fdata(k,imesh,jmesh,isub) + &
     &                              CMPLX(ppm_recvbuffers(ibuffer-1),  &
     &                              ppm_recvbuffers(ibuffer),ppm_kind_double)
#elif  __KIND == __INTEGER
                                 fdata(k,imesh,jmesh,isub) =    &
     &                              fdata(k,imesh,jmesh,isub) + &
     &                              INT(ppm_recvbuffers(ibuffer))
#elif  __KIND == __LOGICAL
                                 IF (ppm_recvbuffers(ibuffer) .GT.     &
     &                              (1.0_ppm_kind_single-ppm_myepss)) THEN
                                    fdata(k,imesh,jmesh,isub) =        &
     &                                  (fdata(k,imesh,jmesh,isub) .AND. .TRUE.)
                                 ELSE
                                    fdata(k,imesh,jmesh,isub) =        &
     &                                  (fdata(k,imesh,jmesh,isub) .AND. .FALSE.)
                                 ENDIF
#endif
                              ENDDO
                           ENDDO
                        ENDDO
                     ENDIF           ! rtype
                  ENDIF            ! lda = ...
               !----------------------------------------------------------------
               !  With mask: This will not vectorize so no need to unroll
               !----------------------------------------------------------------
               ELSE
                  DO jmesh=ylo,yhi
                     DO imesh=xlo,xhi
                        IF (mask(imesh,jmesh,isub)) THEN
                           DO k=1,edim
#if    __KIND == __SINGLE_PRECISION_COMPLEX | __KIND == __DOUBLE_PRECISION_COMPLEX
                              ibuffer = ibuffer + 2
#else
                              ibuffer = ibuffer + 1
#endif
                              IF (rtype .EQ. ppm_param_pop_replace) THEN
#if    __KIND == __DOUBLE_PRECISION
                                 fdata(k,imesh,jmesh,isub) =     &
     &                              REAL(ppm_recvbuffers(ibuffer),ppm_kind_double)
#elif  __KIND == __SINGLE_PRECISION
                                 fdata(k,imesh,jmesh,isub) =     &
     &                              ppm_recvbuffers(ibuffer)
#elif  __KIND == __SINGLE_PRECISION_COMPLEX
                                 fdata(k,imesh,jmesh,isub) =     &
     &                              CMPLX(ppm_recvbuffers(ibuffer-1),  &
     &                              ppm_recvbuffers(ibuffer),ppm_kind_single)
#elif  __KIND == __DOUBLE_PRECISION_COMPLEX
                                 fdata(k,imesh,jmesh,isub) =     &
     &                              CMPLX(ppm_recvbuffers(ibuffer-1),  &
     &                              ppm_recvbuffers(ibuffer),ppm_kind_double)
#elif  __KIND == __INTEGER
                                 fdata(k,imesh,jmesh,isub) =     &
     &                              INT(ppm_recvbuffers(ibuffer))
#elif  __KIND == __LOGICAL
                                 IF (ppm_recvbuffers(ibuffer) .GT.     &
     &                              (1.0_ppm_kind_single-ppm_myepss)) THEN
                                    fdata(k,imesh,jmesh,isub) = .TRUE.
                                 ELSE
                                    fdata(k,imesh,jmesh,isub) = .FALSE.
                                 ENDIF
#endif
                              ELSEIF (rtype .EQ. ppm_param_pop_add) THEN
#if    __KIND == __DOUBLE_PRECISION
                                 fdata(k,imesh,jmesh,isub) =    &
     &                              fdata(k,imesh,jmesh,isub) + &
     &                              REAL(ppm_recvbuffers(ibuffer),ppm_kind_double)
#elif  __KIND == __SINGLE_PRECISION
                                 fdata(k,imesh,jmesh,isub) =    &
     &                              fdata(k,imesh,jmesh,isub) + &
     &                              ppm_recvbuffers(ibuffer)
#elif  __KIND == __SINGLE_PRECISION_COMPLEX
                                 fdata(k,imesh,jmesh,isub) =    &
     &                              fdata(k,imesh,jmesh,isub) + &
     &                              CMPLX(ppm_recvbuffers(ibuffer-1),  &
     &                              ppm_recvbuffers(ibuffer),ppm_kind_single)
#elif  __KIND == __DOUBLE_PRECISION_COMPLEX
                                 fdata(k,imesh,jmesh,isub) =    &
     &                              fdata(k,imesh,jmesh,isub) + &
     &                              CMPLX(ppm_recvbuffers(ibuffer-1),  &
     &                              ppm_recvbuffers(ibuffer),ppm_kind_double)
#elif  __KIND == __INTEGER
                                 fdata(k,imesh,jmesh,isub) =    &
     &                              fdata(k,imesh,jmesh,isub) + &
     &                              INT(ppm_recvbuffers(ibuffer))
#elif  __KIND == __LOGICAL
                                 IF (ppm_recvbuffers(ibuffer) .GT.     &
     &                              (1.0_ppm_kind_single-ppm_myepss)) THEN
                                    fdata(k,imesh,jmesh,isub) =        &
     &                                  (fdata(k,imesh,jmesh,isub) .AND. .TRUE.)
                                 ELSE
                                    fdata(k,imesh,jmesh,isub) =        &
     &                                  (fdata(k,imesh,jmesh,isub) .AND. .FALSE.)
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
               ENDIF            ! PRESENT(mask)
#elif __DIM == __SFIELD
               IF (.NOT.PRESENT(mask)) THEN
                  IF (rtype .EQ. ppm_param_pop_replace) THEN
                     DO jmesh=ylo,yhi
                        DO imesh=xlo,xhi
#if    __KIND == __SINGLE_PRECISION_COMPLEX | __KIND == __DOUBLE_PRECISION_COMPLEX
                           ibuffer = ibuffer + 2
#else
                           ibuffer = ibuffer + 1
#endif
#if    __KIND == __DOUBLE_PRECISION
                           fdata(imesh,jmesh,isub) =     &
     &                        REAL(ppm_recvbuffers(ibuffer),ppm_kind_double)
#elif  __KIND == __SINGLE_PRECISION
                           fdata(imesh,jmesh,isub) =     &
     &                        ppm_recvbuffers(ibuffer)
#elif  __KIND == __SINGLE_PRECISION_COMPLEX
                           fdata(imesh,jmesh,isub) =     &
     &                        CMPLX(ppm_recvbuffers(ibuffer-1),  &
     &                        ppm_recvbuffers(ibuffer),ppm_kind_single)
#elif  __KIND == __DOUBLE_PRECISION_COMPLEX
                           fdata(imesh,jmesh,isub) =     &
     &                        CMPLX(ppm_recvbuffers(ibuffer-1),  &
     &                        ppm_recvbuffers(ibuffer),ppm_kind_double)
#elif  __KIND == __INTEGER
                           fdata(imesh,jmesh,isub) =     &
     &                        INT(ppm_recvbuffers(ibuffer))
#elif  __KIND == __LOGICAL
                           IF (ppm_recvbuffers(ibuffer) .GT.     &
     &                        (1.0_ppm_kind_single-ppm_myepss)) THEN
                              fdata(imesh,jmesh,isub) = .TRUE.
                           ELSE
                              fdata(imesh,jmesh,isub) = .FALSE.
                           ENDIF
#endif
                        ENDDO
                     ENDDO
                  ELSEIF (rtype .EQ. ppm_param_pop_add) THEN
                     DO jmesh=ylo,yhi
                        DO imesh=xlo,xhi
#if    __KIND == __SINGLE_PRECISION_COMPLEX | __KIND == __DOUBLE_PRECISION_COMPLEX
                           ibuffer = ibuffer + 2
#else
                           ibuffer = ibuffer + 1
#endif
#if    __KIND == __DOUBLE_PRECISION
                           fdata(imesh,jmesh,isub) =    &
     &                        fdata(imesh,jmesh,isub) + &
     &                        REAL(ppm_recvbuffers(ibuffer),ppm_kind_double)
#elif  __KIND == __SINGLE_PRECISION
                           fdata(imesh,jmesh,isub) =    &
     &                        fdata(imesh,jmesh,isub) + &
     &                        ppm_recvbuffers(ibuffer)
#elif  __KIND == __SINGLE_PRECISION_COMPLEX
                           fdata(imesh,jmesh,isub) =    &
     &                        fdata(imesh,jmesh,isub) + &
     &                        CMPLX(ppm_recvbuffers(ibuffer-1),  &
     &                        ppm_recvbuffers(ibuffer),ppm_kind_single)
#elif  __KIND == __DOUBLE_PRECISION_COMPLEX
                           fdata(imesh,jmesh,isub) =    &
     &                        fdata(imesh,jmesh,isub) + &
     &                        CMPLX(ppm_recvbuffers(ibuffer-1),  &
     &                        ppm_recvbuffers(ibuffer),ppm_kind_double)
#elif  __KIND == __INTEGER
                           fdata(imesh,jmesh,isub) =    &
     &                        fdata(imesh,jmesh,isub) + &
     &                        INT(ppm_recvbuffers(ibuffer))
#elif  __KIND == __LOGICAL
                           IF (ppm_recvbuffers(ibuffer) .GT.     &
     &                        (1.0_ppm_kind_single-ppm_myepss)) THEN
                              fdata(imesh,jmesh,isub) =        &
     &                            (fdata(imesh,jmesh,isub) .AND. .TRUE.)
                           ELSE
                              fdata(imesh,jmesh,isub) =        &
     &                            (fdata(imesh,jmesh,isub) .AND. .FALSE.)
                           ENDIF
#endif
                        ENDDO
                     ENDDO
                  ENDIF
              ELSE            ! IF PRESENT(mask)
                  DO jmesh=ylo,yhi
                     DO imesh=xlo,xhi
                        IF (mask(imesh,jmesh,isub)) THEN
#if    __KIND == __SINGLE_PRECISION_COMPLEX | __KIND == __DOUBLE_PRECISION_COMPLEX
                           ibuffer = ibuffer + 2
#else
                           ibuffer = ibuffer + 1
#endif
                           IF (rtype .EQ. ppm_param_pop_replace) THEN
#if    __KIND == __DOUBLE_PRECISION
                              fdata(imesh,jmesh,isub) =     &
     &                           REAL(ppm_recvbuffers(ibuffer),ppm_kind_double)
#elif  __KIND == __SINGLE_PRECISION
                              fdata(imesh,jmesh,isub) =     &
     &                           ppm_recvbuffers(ibuffer)
#elif  __KIND == __SINGLE_PRECISION_COMPLEX
                              fdata(imesh,jmesh,isub) =     &
     &                           CMPLX(ppm_recvbuffers(ibuffer-1),  &
     &                           ppm_recvbuffers(ibuffer),ppm_kind_single)
#elif  __KIND == __DOUBLE_PRECISION_COMPLEX
                              fdata(imesh,jmesh,isub) =     &
     &                           CMPLX(ppm_recvbuffers(ibuffer-1),  &
     &                           ppm_recvbuffers(ibuffer),ppm_kind_double)
#elif  __KIND == __INTEGER
                              fdata(imesh,jmesh,isub) =     &
     &                           INT(ppm_recvbuffers(ibuffer))
#elif  __KIND == __LOGICAL
                              IF (ppm_recvbuffers(ibuffer) .GT.     &
     &                           (1.0_ppm_kind_single-ppm_myepss)) THEN
                                 fdata(imesh,jmesh,isub) = .TRUE.
                              ELSE
                                 fdata(imesh,jmesh,isub) = .FALSE.
                              ENDIF
#endif
                           ELSEIF (rtype .EQ. ppm_param_pop_add) THEN
#if    __KIND == __DOUBLE_PRECISION
                              fdata(imesh,jmesh,isub) =    &
     &                           fdata(imesh,jmesh,isub) + &
     &                           REAL(ppm_recvbuffers(ibuffer),ppm_kind_double)
#elif  __KIND == __SINGLE_PRECISION
                              fdata(imesh,jmesh,isub) =    &
     &                           fdata(imesh,jmesh,isub) + &
     &                           ppm_recvbuffers(ibuffer)
#elif  __KIND == __SINGLE_PRECISION_COMPLEX
                              fdata(imesh,jmesh,isub) =    &
     &                           fdata(imesh,jmesh,isub) + &
     &                           CMPLX(ppm_recvbuffers(ibuffer-1),  &
     &                           ppm_recvbuffers(ibuffer),ppm_kind_single)
#elif  __KIND == __DOUBLE_PRECISION_COMPLEX
                              fdata(imesh,jmesh,isub) =    &
     &                           fdata(imesh,jmesh,isub) + &
     &                           CMPLX(ppm_recvbuffers(ibuffer-1),  &
     &                           ppm_recvbuffers(ibuffer),ppm_kind_double)
#elif  __KIND == __INTEGER
                              fdata(imesh,jmesh,isub) =    &
     &                           fdata(imesh,jmesh,isub) + &
     &                           INT(ppm_recvbuffers(ibuffer))
#elif  __KIND == __LOGICAL
                              IF (ppm_recvbuffers(ibuffer) .GT.     &
     &                           (1.0_ppm_kind_single-ppm_myepss)) THEN
                                 fdata(imesh,jmesh,isub) =        &
     &                               (fdata(imesh,jmesh,isub) .AND. .TRUE.)
                              ELSE
                                 fdata(imesh,jmesh,isub) =        &
     &                               (fdata(imesh,jmesh,isub) .AND. .FALSE.)
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
               ENDIF          ! PRESENT(mask)
#endif
            ENDDO             ! ppm_precvbuffer
         ENDDO                ! ppm_nrecvlist
      ENDIF                   ! ppm_kind

 8888 CONTINUE
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
              CALL ppm_error(ppm_err_alloc,'ppm_map_field_pop_2d',     &
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
          CALL ppm_error(ppm_err_dealloc,'ppm_map_field_pop_2d',     &
     &        'inverse sub list INVSUBLIST',__LINE__,info)
      ENDIF

      !-------------------------------------------------------------------------
      !  Return
      !-------------------------------------------------------------------------
 9999 CONTINUE
      CALL substop('ppm_map_field_pop_2d',t0,info)
      RETURN
      CONTAINS
      SUBROUTINE check
          CALL ppm_check_meshid(target_topoid,target_meshid,ldo,info)
          IF (ppm_buffer_set .LT. 1) THEN
              info = ppm_error_notice
              CALL ppm_error(ppm_err_argument,'ppm_map_field_pop_2d',  &
     &            'buffer is empty. Cannot pop.',__LINE__,info)
              GOTO 7777
          ENDIF
          IF (.NOT. ldo) THEN
              info = ppm_error_error
              CALL ppm_error(ppm_err_argument,'ppm_map_field_pop_2d',  &
     &            'destination meshid out of range',__LINE__,info)
              GOTO 7777
          ENDIF
#if   __DIM == __VFIELD
          IF (lda .LT. 1) THEN
              info = ppm_error_error
              CALL ppm_error(ppm_err_argument,'ppm_map_field_pop_2d',  &
     &            'lda must be >0 for vector data',__LINE__,info)
              GOTO 7777
          ENDIF
#elif __DIM == __SFIELD
          IF (lda .NE. 1) THEN
              info = ppm_error_error
              CALL ppm_error(ppm_err_argument,'ppm_map_field_pop_2d',  &
     &            'lda must be =1 for scalar data',__LINE__,info)
              GOTO 7777
          ENDIF
#endif
          IF (SIZE(ghostsize,1) .LT. 2) THEN
              info = ppm_error_error
              CALL ppm_error(ppm_err_argument,'ppm_map_field_pop_2d',  &
     &            'ghostsize must be given for all dimensions',__LINE__,info)
              GOTO 7777
          ENDIF
          IF ((ghostsize(1) .LT. 0) .OR. (ghostsize(2) .LT. 0)) THEN
              info = ppm_error_error
              CALL ppm_error(ppm_err_argument,'ppm_map_field_pop_2d',  &
     &            'ghostsize must be >=0 in all dimensions',__LINE__,info)
              GOTO 7777
          ENDIF
          IF (PRESENT(poptype)) THEN
              IF ((poptype .NE. ppm_param_pop_replace) .AND.     &
     &            (poptype .NE. ppm_param_pop_add)) THEN
                  info = ppm_error_error
                  CALL ppm_error(ppm_err_argument,'ppm_map_field_pop_2d',  &
     &                'Unknown pop type specified',__LINE__,info)
                  GOTO 7777
              ENDIF
          ENDIF
 7777     CONTINUE
      END SUBROUTINE check
      SUBROUTINE check_two
          IF (PRESENT(mask)) THEN
              xhi = 0
              yhi = 0
              DO i=1,target_topo%nsublist
                  isub = target_topo%isublist(i)
                  IF (target_mesh%nnodes(1,isub).GT.xhi) THEN
                      xhi = target_mesh%nnodes(1,isub)
                  ENDIF
                  IF (target_mesh%nnodes(2,isub).GT.yhi) THEN
                      yhi = target_mesh%nnodes(2,isub)
                  ENDIF
              ENDDO
              IF (SIZE(mask,1) .LT. xhi) THEN
                  info = ppm_error_error
                  CALL ppm_error(ppm_err_argument,'ppm_map_field_pop_2d',  &
     &                'x dimension of mask does not match mesh',__LINE__,info)
                  GOTO 6666
              ENDIF
              IF (SIZE(mask,2) .LT. yhi) THEN
                  info = ppm_error_error
                  CALL ppm_error(ppm_err_argument,'ppm_map_field_pop_2d',  &
     &                'y dimension of mask does not match mesh',__LINE__,info)
                  GOTO 6666
              ENDIF
              IF (SIZE(mask,3) .LT. target_topo%nsublist) THEN
                  info = ppm_error_error
                  CALL ppm_error(ppm_err_argument,'ppm_map_field_pop_2d',  &
     &                'mask data for some subs is missing',__LINE__,info)
                  GOTO 6666
              ENDIF
          ENDIF
 6666     CONTINUE
     END SUBROUTINE check_two
#if    __DIM == __SFIELD
#if    __KIND == __SINGLE_PRECISION
      END SUBROUTINE ppm_map_field_pop_2d_sca_s
#elif  __KIND == __DOUBLE_PRECISION
      END SUBROUTINE ppm_map_field_pop_2d_sca_d
#elif  __KIND == __SINGLE_PRECISION_COMPLEX
      END SUBROUTINE ppm_map_field_pop_2d_sca_sc
#elif  __KIND == __DOUBLE_PRECISION_COMPLEX
      END SUBROUTINE ppm_map_field_pop_2d_sca_dc
#elif  __KIND == __INTEGER
      END SUBROUTINE ppm_map_field_pop_2d_sca_i
#elif  __KIND == __LOGICAL
      END SUBROUTINE ppm_map_field_pop_2d_sca_l
#endif

#elif  __DIM == __VFIELD
#if    __KIND == __SINGLE_PRECISION
      END SUBROUTINE ppm_map_field_pop_2d_vec_s
#elif  __KIND == __DOUBLE_PRECISION
      END SUBROUTINE ppm_map_field_pop_2d_vec_d
#elif  __KIND == __SINGLE_PRECISION_COMPLEX
      END SUBROUTINE ppm_map_field_pop_2d_vec_sc
#elif  __KIND == __DOUBLE_PRECISION_COMPLEX
      END SUBROUTINE ppm_map_field_pop_2d_vec_dc
#elif  __KIND == __INTEGER
      END SUBROUTINE ppm_map_field_pop_2d_vec_i
#elif  __KIND == __LOGICAL
      END SUBROUTINE ppm_map_field_pop_2d_vec_l
#endif
#endif
