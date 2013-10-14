#if    __DIM == __SFIELD
#if    __KIND == __SINGLE_PRECISION
      SUBROUTINE mesh_map_pop_3d_sca_s(this,fdata_dummy,&
     &                   p_idx,info,mask,poptype)
#elif  __KIND == __DOUBLE_PRECISION
      SUBROUTINE mesh_map_pop_3d_sca_d(this,fdata_dummy,&
     &                   p_idx,info,mask,poptype)
#elif  __KIND == __SINGLE_PRECISION_COMPLEX
      SUBROUTINE mesh_map_pop_3d_sca_sc(this,fdata_dummy,&
     &                   p_idx,info,mask,poptype)
#elif  __KIND == __DOUBLE_PRECISION_COMPLEX
      SUBROUTINE mesh_map_pop_3d_sca_dc(this,fdata_dummy,&
     &                   p_idx,info,mask,poptype)
#elif  __KIND == __INTEGER
      SUBROUTINE mesh_map_pop_3d_sca_i(this,fdata_dummy,&
     &                   p_idx,info,mask,poptype)
#elif  __KIND == __LOGICAL
      SUBROUTINE mesh_map_pop_3d_sca_l(this,fdata_dummy,&
     &                   p_idx,info,mask,poptype)
#endif
#elif  __DIM == __VFIELD
#if    __KIND == __SINGLE_PRECISION
      SUBROUTINE mesh_map_pop_3d_vec_s(this,fdata_dummy,&
     &                   lda,info,mask,poptype)
#elif  __KIND == __DOUBLE_PRECISION
      SUBROUTINE mesh_map_pop_3d_vec_d(this,fdata_dummy,&
     &                   lda,p_idx,info,mask,poptype)
#elif  __KIND == __SINGLE_PRECISION_COMPLEX
      SUBROUTINE mesh_map_pop_3d_vec_sc(this,fdata_dummy,&
     &                   lda,p_idx,info,mask,poptype)
#elif  __KIND == __DOUBLE_PRECISION_COMPLEX
      SUBROUTINE mesh_map_pop_3d_vec_dc(this,fdata_dummy,&
     &                   lda,p_idx,info,mask,poptype)
#elif  __KIND == __INTEGER
      SUBROUTINE mesh_map_pop_3d_vec_i(this,fdata_dummy,&
     &                   lda,p_idx,info,mask,poptype)
#elif  __KIND == __LOGICAL
      SUBROUTINE mesh_map_pop_3d_vec_l(this,fdata_dummy,&
     &                   lda,p_idx,info,mask,poptype)
#endif
#endif
      !!! This routine pops the list buffer for 3D mesh data
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
      USE ppm_module_data_mesh
      USE ppm_module_topo_typedef
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
      INTEGER , DIMENSION(:,:,:)   , POINTER         :: fdata_dummy
      INTEGER , DIMENSION(:,:,:)   , POINTER         :: fdata
#elif __KIND == __LOGICAL
      LOGICAL , DIMENSION(:,:,:)   , POINTER         :: fdata_dummy
      LOGICAL , DIMENSION(:,:,:)   , POINTER         :: fdata
#elif __KIND == __SINGLE_PRECISION_COMPLEX | __KIND == __DOUBLE_PRECISION_COMPLEX
      COMPLEX(MK), DIMENSION(:,:,:), POINTER         :: fdata_dummy
      COMPLEX(MK), DIMENSION(:,:,:), POINTER         :: fdata
#else
      REAL(MK), DIMENSION(:,:,:)   , POINTER         :: fdata_dummy
      REAL(MK), DIMENSION(:,:,:)   , POINTER         :: fdata
#endif
#elif __DIM == __VFIELD
#if   __KIND == __INTEGER
      INTEGER , DIMENSION(:,:,:,:)   , POINTER       :: fdata_dummy
      INTEGER , DIMENSION(:,:,:,:)   , POINTER       :: fdata
#elif __KIND == __LOGICAL
      LOGICAL , DIMENSION(:,:,:,:)   , POINTER       :: fdata_dummy
      LOGICAL , DIMENSION(:,:,:,:)   , POINTER       :: fdata
#elif __KIND == __SINGLE_PRECISION_COMPLEX | __KIND == __DOUBLE_PRECISION_COMPLEX
      COMPLEX(MK), DIMENSION(:,:,:,:), POINTER       :: fdata_dummy
      COMPLEX(MK), DIMENSION(:,:,:,:), POINTER       :: fdata
#else
      REAL(MK), DIMENSION(:,:,:,:)   , POINTER       :: fdata_dummy
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
      INTEGER                        , INTENT(IN   ) :: p_idx
      !!! The index where the data is stored on the subpatches
      CLASS(ppm_t_equi_mesh)                         :: this
      !!! Source mesh
      INTEGER                        , INTENT(  OUT) :: info
      !!! Returns status, 0 upon success
      !-------------------------------------------------------------------------
      !  Local variables
      !-------------------------------------------------------------------------
      INTEGER, DIMENSION(3)   :: mofs,patchid
      INTEGER, DIMENSION(4)   :: ldu,ldl
      INTEGER                 :: i,j,k,ibuffer,Mdata,isub,bdim,jsub,edim
      INTEGER                 :: ipatch
      INTEGER                 :: iopt,totopo,xhi,yhi,zhi,xlo,ylo,zlo
      INTEGER                 :: imesh,jmesh,kmesh
      INTEGER                 :: btype,idom
      LOGICAL                 :: found_patch
      INTEGER                 :: rtype
      CHARACTER(LEN=ppm_char) :: mesg
#if   __DIM == __SFIELD
      INTEGER, PARAMETER    :: lda = 1
#endif
      TYPE(ppm_t_topo), POINTER  :: target_topo => NULL()
      !-------------------------------------------------------------------------
      !  Externals
      !-------------------------------------------------------------------------

      !-------------------------------------------------------------------------
      !  Initialise
      !-------------------------------------------------------------------------
      start_subroutine("mesh_map_pop_3d")

      !-------------------------------------------------------------------------
      !  Check arguments
      !-------------------------------------------------------------------------
      IF (ppm_debug .GT. 0) THEN
        CALL check
        IF (info .NE. 0) GOTO 9999
      ENDIF

      ! skip if buffer empty
      IF (ppm_buffer_set .LT. 1) THEN
        IF (ppm_debug .GT. 1) THEN
            info = ppm_error_notice
            CALL ppm_error(ppm_err_buffer_empt,caller,    &
     &          'Buffer is empty: skipping pop!',__LINE__,info)
             info = 0
        ENDIF
        GOTO 9999
      ENDIF

      !-------------------------------------------------------------------------
      !  set the local pointers to the topology and mesh
      !-------------------------------------------------------------------------
      target_topo => ppm_topo(this%topoid)%t

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
         CALL ppm_error(ppm_err_argument,caller,  &
         & 'mesh_map_pop cannot be called after ghost_init',__LINE__,info)
         CALL ppm_error(ppm_err_argument,caller,  &
         &  'mesh_map_ghost_init must not be called directly by the user',__LINE__,info)
         GOTO 9999
      ELSE
         IF (PRESENT(poptype)) THEN
            rtype = poptype
         ELSE
            IF (ppm_map_type.EQ.ppm_param_map_ghost_put) THEN
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
          CALL ppm_write(ppm_rank,caller,mesg,info)
      ENDIF
#if   __DIM == __VFIELD
      IF (edim.NE.lda) THEN
         fail("leading dimension LDA is in error",ppm_err_wrong_dim)
      ENDIF
#elif __DIM == __SFIELD
      IF (edim.NE.1) THEN
         fail("buffer does not contain 1d data!",ppm_err_wrong_dim)
      ENDIF
#endif

      !-------------------------------------------------------------------------
      !  Check that the required type is identical to the type of the buffer
      !-------------------------------------------------------------------------
      btype = ppm_buffer_type(ppm_buffer_set)
#if    __KIND == __SINGLE_PRECISION
      IF (btype.NE.ppm_kind_single) THEN
         info = ppm_error_error
         CALL ppm_error(ppm_err_wrong_prec,caller,    &
     &       'trying to pop a non-single into single ',__LINE__,info)
         GOTO 9999
      ENDIF
#elif  __KIND == __DOUBLE_PRECISION
      IF (btype.NE.ppm_kind_double) THEN
         info = ppm_error_error
         CALL ppm_error(ppm_err_wrong_prec,caller,    &
     &       'trying to pop a non-double into double ',__LINE__,info)
         GOTO 9999
      ENDIF
#elif  __KIND == __SINGLE_PRECISION_COMPLEX
      IF (btype.NE.ppm_kind_single) THEN
         info = ppm_error_error
         CALL ppm_error(ppm_err_wrong_prec,caller,    &
     &       'trying to pop a non-single-complex into single-complex',&
     &       __LINE__,info)
         GOTO 9999
      ENDIF
#elif  __KIND == __DOUBLE_PRECISION_COMPLEX
      IF (btype.NE.ppm_kind_double) THEN
         info = ppm_error_error
         CALL ppm_error(ppm_err_wrong_prec,caller,    &
     &       'trying to pop a non-double-complex into double-complex',&
     &       __LINE__,info)
         GOTO 9999
      ENDIF
#elif  __KIND == __INTEGER
      IF (btype.NE.ppm_integer) THEN
         info = ppm_error_error
         CALL ppm_error(ppm_err_wrong_prec,caller,    &
     &       'trying to pop a non-integer into integer ',__LINE__,info)
         GOTO 9999
      ENDIF
#elif  __KIND == __LOGICAL
      IF (btype.NE.ppm_logical) THEN
         info = ppm_error_error
         CALL ppm_error(ppm_err_wrong_prec,caller,    &
     &       'trying to pop a non-logical into logical ',__LINE__,info)
      ENDIF
#endif

      !-------------------------------------------------------------------------
      !  Build the inverse sub list to find local sub indeices based on
      !  global ones (the global ones are communicated)
      !-------------------------------------------------------------------------
      iopt   = ppm_param_alloc_fit
      ldu(1) = target_topo%nsublist
      CALL ppm_alloc(sublist,ldu,iopt,info)
      or_fail_alloc("sublist")
      ! We need to copy it into a temp list, since directly using
      ! ppm_isublist(:,ppm_field_topoid) as an argument to invert_list is
      ! not possible since the argument needs to be a POINTER.
      sublist(1:ldu(1)) = target_topo%isublist(1:ldu(1))
      CALL ppm_util_invert_list(sublist,invsublist,info)
      or_fail("ppm_util_invert_list")

      iopt   = ppm_param_dealloc
      CALL ppm_alloc(sublist,ldu,iopt,info)
      or_fail_dealloc("sublist")

      !-------------------------------------------------------------------------
      !  Determine the number of data points to be received
      !-------------------------------------------------------------------------
      Mdata = 0
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
              &       ppm_mesh_irecvblksize(2,j)* ppm_mesh_irecvblksize(3,j))
          ENDDO
      ENDDO

      !-------------------------------------------------------------------------
      !  If there is nothing to be sent we are done
      !-------------------------------------------------------------------------
      IF (Mdata .EQ. 0) THEN
          IF (ppm_debug .GT. 0) THEN
              CALL ppm_write(ppm_rank,caller,   &
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
          CALL ppm_write(ppm_rank,caller,mesg,info)
      ENDIF
      ppm_nrecvbuffer = ppm_nrecvbuffer - Mdata*bdim

      ibuffer = ppm_nrecvbuffer
      IF (ppm_debug .GT. 1) THEN
          WRITE(mesg,'(A,I9)') 'ibuffer = ',ibuffer
          CALL ppm_write(ppm_rank,caller,mesg,info)
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
     &         ppm_mesh_isendblksize(2,j) * ppm_mesh_isendblksize(3,j))
         ENDDO
      ENDDO
      ppm_nsendbuffer = ppm_nsendbuffer - ppm_buffer_dim(ppm_buffer_set)*Mdata

      !-------------------------------------------------------------------------
      !  Debug output
      !-------------------------------------------------------------------------
      IF (ppm_debug .GT. 0) THEN
          IF (rtype .EQ. ppm_param_pop_replace) THEN
              CALL ppm_write(ppm_rank,caller,     &
     &           'Replacing current field values',info)
          ELSEIF(rtype .EQ. ppm_param_pop_add) THEN
              CALL ppm_write(ppm_rank,caller,     &
     &           'Adding to current field values',info)
          ENDIF
      ENDIF

      !-------------------------------------------------------------------------
      !  Loop over the received mesh blocks
      !-------------------------------------------------------------------------
      !-------------------------------------------------------------------------
      !  DOUBLE PRECISION
      !-------------------------------------------------------------------------
      SELECT CASE (ppm_kind)
      CASE (ppm_kind_double)
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
               !  Get the patch ID for this mesh block
               !----------------------------------------------------------------
               patchid(1:3) = ppm_mesh_irecvpatchid(1:3,j)
               !----------------------------------------------------------------
               !  Translate to local sub ID for storing the data
               !----------------------------------------------------------------
               isub = invsublist(jsub)
               !----------------------------------------------------------------
               !  Get pointer to the data for this sub, this field and this block
               ! TODO: room for improvement!...
               !----------------------------------------------------------------
               !something like that may be nice?
               !CALL this%get_field_on_patch(fdata,isub,info)
               !or_fail("could not get_field_on_patch for this sub")
               !(lazy) search for the subpatch that has the right global id
               found_patch = .FALSE.
               patches: DO ipatch=1,this%subpatch_by_sub(jsub)%nsubpatch
                   fdata => NULL()
                   SELECT TYPE(p => this%subpatch_by_sub(jsub)%vec(ipatch)%t)
                   TYPE IS (ppm_t_subpatch)
                       IF (ALL(p%istart_p.EQ.patchid)) THEN
                            found_patch = .TRUE.
                            !------------------------------------------------
                            !  Determine size of field data array needed
                            !------------------------------------------------
                            xhi = p%nnodes(1)
                            yhi = p%nnodes(2)
                            zhi = p%nnodes(3)

                            !------------------------------------------------
                            !  Reallocate array if needed
                            !------------------------------------------------
                            IF ((ppm_map_type .EQ. ppm_param_map_ghost_get) &
                                & .OR.   &
                                &   (ppm_map_type .EQ. ppm_param_map_ghost_put) &
                                & .OR.   &
                                &    (rtype .EQ. ppm_param_pop_add)) THEN
                            !------------------------------------------------
                            !  Preserve old fields if this is to receive ghosts
                            !  or to add contributions
                            !------------------------------------------------
                                iopt = ppm_param_alloc_fit_preserve
                            ELSE
                                iopt = ppm_param_alloc_fit
                            ENDIF
!#if   __DIM == __VFIELD
                            !ldu(1) = edim
                            !ldu(2) = xhi+p%ghostsize(2)
                            !ldu(3) = yhi+p%ghostsize(4)
                            !ldu(4) = zhi+p%ghostsize(6)
                            !ldl(1) = 1
                            !ldl(2) = 1-p%ghostsize(1)
                            !ldl(3) = 1-p%ghostsize(3)
                            !ldl(4) = 1-p%ghostsize(5)
!#elif __DIM == __SFIELD
                            !ldu(1) = xhi+p%ghostsize(2)
                            !ldu(2) = yhi+p%ghostsize(4)
                            !ldu(3) = zhi+p%ghostsize(6)
                            !ldl(1) = 1-p%ghostsize(1)
                            !ldl(2) = 1-p%ghostsize(3)
                            !ldl(3) = 1-p%ghostsize(5)
!#endif

                            check_associated(<#p%subpatch_data#>)
                            check_true(<#p%subpatch_data%exists(p_idx)#>,"does not exist")

#if __KIND == __DOUBLE_PRECISION
#if    __DIM == __SFIELD
                            !call ppm_alloc(&
                                !p%subpatch_data%vec(p_idx)%t%data_3d_rd,&
                                !ldl,ldu,iopt,info)
                            fdata => p%subpatch_data%vec(p_idx)%t%data_3d_rd
#elif  __DIM == __VFIELD
                            !call ppm_alloc(&
                                !p%subpatch_data%vec(p_idx)%t%data_4d_rd,&
                                !ldl,ldu,iopt,info)
                            fdata => p%subpatch_data%vec(p_idx)%t%data_4d_rd
#endif
#else
                            write(*,*) "WRONG TYPE!!!!"
#endif
                            !------------------------------------------------------
                            !  Mesh offset for this subpatch
                            !------------------------------------------------------
                            mofs(1) = p%istart(1)-1
                            mofs(2) = p%istart(2)-1
                            mofs(3) = p%istart(3)-1
                            !------------------------------------------------------
                            !  Get boundaries of mesh block to be received in
                            !  local sub coordinates
                            !------------------------------------------------------
                            xlo = ppm_mesh_irecvblkstart(1,j)-mofs(1)
                            ylo = ppm_mesh_irecvblkstart(2,j)-mofs(2)
                            zlo = ppm_mesh_irecvblkstart(3,j)-mofs(3)
                            xhi = xlo+ppm_mesh_irecvblksize(1,j)-1
                            yhi = ylo+ppm_mesh_irecvblksize(2,j)-1
                            zhi = zlo+ppm_mesh_irecvblksize(3,j)-1
                            IF (ppm_debug .GT. 1) THEN
                                stdout("isub = ",isub," jsub = ",jsub)
                                stdout("p%istart_p",'p%istart_p')
                                stdout("p%iend_p",'p%iend_p')
                                stdout("p%istart",'p%istart')
                                stdout("p%iend",'p%iend')
                                stdout("patchid = ",patchid)
                                WRITE(mesg,'(A,3I4)') 'start: ',             &
                                    &  ppm_mesh_irecvblkstart(1,j),&
                                    &  ppm_mesh_irecvblkstart(2,j),&
                                    &  ppm_mesh_irecvblkstart(3,j)
                                CALL ppm_write(ppm_rank,caller,mesg,info)
                                WRITE(mesg,'(A,3I4)') 'size: ',             &
                                    &  ppm_mesh_irecvblksize(1,j),&
                                    &  ppm_mesh_irecvblksize(2,j),&
                                    &  ppm_mesh_irecvblksize(3,j)
                                CALL ppm_write(ppm_rank,caller,mesg,info)
                                WRITE(mesg,'(A,3I4)')'size_b: ',&
                                    xhi-xlo+1,yhi-ylo+1,zhi-zlo+1
                                CALL ppm_write(ppm_rank,caller,mesg,info)
                                WRITE(mesg,'(A,3I4)') 'mesh offset: ',&
                                    mofs(1),mofs(2),mofs(3)
                                CALL ppm_write(ppm_rank,caller,mesg,info)
                                WRITE(mesg,'(A,2I4)') 'xlo, xhi: ',xlo,xhi
                                CALL ppm_write(ppm_rank,caller,mesg,info)
                                WRITE(mesg,'(A,2I4)') 'ylo, yhi: ',ylo,yhi
                                CALL ppm_write(ppm_rank,caller,mesg,info)
                                WRITE(mesg,'(A,2I4)') 'zlo, zhi: ',zlo,zhi
                                CALL ppm_write(ppm_rank,caller,mesg,info)
                                WRITE(mesg,'(A,I1)') 'buffer dim: ',edim
                                CALL ppm_write(ppm_rank,caller,mesg,info)
                            ENDIF

                            !For ghost_get:
                            !check that real mesh nodes are not touched
                            check_false(<#ppm_map_type .EQ. ppm_param_map_ghost_get .AND. (xhi.GE.1 .AND. xlo.LE.p%nnodes(1) .AND. yhi.GE.1 .AND. ylo.LE.p%nnodes(2) .AND. zhi.GE.1 .AND. zlo.LE.p%nnodes(3))#>)
                            !for ghost_put:
                            !check that ghost mesh nodes are not touched
                            check_false(<#ppm_map_type .EQ. ppm_param_map_ghost_put .AND.  (xhi.LT.1 .OR. xlo.GT.p%nnodes(1) .OR. yhi.LT.1 .OR. ylo.GT.p%nnodes(2) .OR. zhi.LT.1 .OR. zlo.GT.p%nnodes(3))#>)
                            !check that we dont access out-of-bounds elements
                            check_true(<#(xlo.GE.p%lo_a(1))#>)
                            check_true(<#(xhi.LE.p%hi_a(1))#>)
                            check_true(<#(ylo.GE.p%lo_a(2))#>)
                            check_true(<#(yhi.LE.p%hi_a(2))#>)
                            check_true(<#(zlo.GE.p%lo_a(3))#>)
                            check_true(<#(zhi.LE.p%hi_a(3))#>)
                            check_associated(fdata)

                            exit patches
                       ENDIF
                   END SELECT
               ENDDO patches

               IF (.NOT. found_patch) THEN
                   stdout("isub = ",isub," jsub = ",jsub," ipatch = ",ipatch)
                   stdout("patchid = ",patchid)
                   stdout("patchid/h = ",'(patchid-1._mk)*this%h(1:ppm_dim)')
                   stdout("h = ",'this%h(1:ppm_dim)')
                   stdout("this%subpatch_by_sub(jsub)%nsubpatch = ",&
                          'this%subpatch_by_sub(jsub)%nsubpatch')
                   stdout("min_sub(jsub)=",'target_topo%min_subd(1:ppm_dim,jsub)')
                   stdout("max_sub(jsub)=",'target_topo%max_subd(1:ppm_dim,jsub)')
                   fail("could not find a patch on this sub with the right global id")
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
                              fdata(1,imesh,jmesh,kmesh) =     &
     &                           REAL(ppm_recvbufferd(ibuffer),ppm_kind_single)
#elif  __KIND == __DOUBLE_PRECISION
                              fdata(1,imesh,jmesh,kmesh) =     &
     &                           ppm_recvbufferd(ibuffer)
#elif  __KIND == __SINGLE_PRECISION_COMPLEX
                              fdata(1,imesh,jmesh,kmesh) =     &
     &                           CMPLX(ppm_recvbufferd(ibuffer-1),  &
     &                           ppm_recvbufferd(ibuffer),ppm_kind_single)
#elif  __KIND == __DOUBLE_PRECISION_COMPLEX
                              fdata(1,imesh,jmesh,kmesh) =     &
     &                           CMPLX(ppm_recvbufferd(ibuffer-1),  &
     &                           ppm_recvbufferd(ibuffer),ppm_kind_double)
#elif  __KIND == __INTEGER
                              fdata(1,imesh,jmesh,kmesh) =     &
     &                           INT(ppm_recvbufferd(ibuffer))
#elif  __KIND == __LOGICAL
                              IF (ppm_recvbufferd(ibuffer) .GT.     &
     &                           (1.0_ppm_kind_double-ppm_myepsd)) THEN
                                 fdata(1,imesh,jmesh,kmesh) = .TRUE.
                              ELSE
                                 fdata(1,imesh,jmesh,kmesh) = .FALSE.
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
                              fdata(1,imesh,jmesh,kmesh) =   &
     &                           fdata(1,imesh,jmesh,kmesh)+ &
     &                           REAL(ppm_recvbufferd(ibuffer),ppm_kind_single)
#elif  __KIND == __DOUBLE_PRECISION
                              fdata(1,imesh,jmesh,kmesh) =   &
     &                           fdata(1,imesh,jmesh,kmesh)+ &
     &                           ppm_recvbufferd(ibuffer)
#elif  __KIND == __SINGLE_PRECISION_COMPLEX
                              fdata(1,imesh,jmesh,kmesh) =   &
     &                           fdata(1,imesh,jmesh,kmesh)+ &
     &                           CMPLX(ppm_recvbufferd(ibuffer-1),  &
     &                           ppm_recvbufferd(ibuffer),ppm_kind_single)
#elif  __KIND == __DOUBLE_PRECISION_COMPLEX
                              fdata(1,imesh,jmesh,kmesh) =   &
     &                           fdata(1,imesh,jmesh,kmesh)+ &
     &                           CMPLX(ppm_recvbufferd(ibuffer-1),  &
     &                           ppm_recvbufferd(ibuffer),ppm_kind_double)
#elif  __KIND == __INTEGER
                              fdata(1,imesh,jmesh,kmesh) =   &
     &                           fdata(1,imesh,jmesh,kmesh)+ &
     &                           INT(ppm_recvbufferd(ibuffer))
#elif  __KIND == __LOGICAL
                              IF (ppm_recvbufferd(ibuffer) .GT.     &
     &                           (1.0_ppm_kind_double-ppm_myepsd)) THEN
                                 fdata(1,imesh,jmesh,kmesh) =    &
     &                              (fdata(1,imesh,jmesh,kmesh) .AND. .TRUE.)
                              ELSE
                                 fdata(1,imesh,jmesh,kmesh) =    &
     &                              (fdata(1,imesh,jmesh,kmesh) .AND. .FALSE.)
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
                              fdata(1,imesh,jmesh,kmesh) =     &
     &                           REAL(ppm_recvbufferd(ibuffer),ppm_kind_single)
                              ibuffer = ibuffer + 1
                              fdata(2,imesh,jmesh,kmesh) =     &
     &                           REAL(ppm_recvbufferd(ibuffer),ppm_kind_single)
#elif  __KIND == __DOUBLE_PRECISION
                              fdata(1,imesh,jmesh,kmesh) =     &
     &                           ppm_recvbufferd(ibuffer)
                              ibuffer = ibuffer + 1
                              fdata(2,imesh,jmesh,kmesh) =     &
     &                           ppm_recvbufferd(ibuffer)
#elif  __KIND == __SINGLE_PRECISION_COMPLEX
                              fdata(1,imesh,jmesh,kmesh) =     &
     &                           CMPLX(ppm_recvbufferd(ibuffer-1),  &
     &                           ppm_recvbufferd(ibuffer),ppm_kind_single)
                              ibuffer = ibuffer + 2
                              fdata(2,imesh,jmesh,kmesh) =     &
     &                           CMPLX(ppm_recvbufferd(ibuffer-1),  &
     &                           ppm_recvbufferd(ibuffer),ppm_kind_single)
#elif  __KIND == __DOUBLE_PRECISION_COMPLEX
                              fdata(1,imesh,jmesh,kmesh) =     &
     &                           CMPLX(ppm_recvbufferd(ibuffer-1),  &
     &                           ppm_recvbufferd(ibuffer),ppm_kind_double)
                              ibuffer = ibuffer + 2
                              fdata(2,imesh,jmesh,kmesh) =     &
     &                           CMPLX(ppm_recvbufferd(ibuffer-1),  &
     &                           ppm_recvbufferd(ibuffer),ppm_kind_double)
#elif  __KIND == __INTEGER
                              fdata(1,imesh,jmesh,kmesh) =     &
     &                           INT(ppm_recvbufferd(ibuffer))
                              ibuffer = ibuffer + 1
                              fdata(2,imesh,jmesh,kmesh) =     &
     &                           INT(ppm_recvbufferd(ibuffer))
#elif  __KIND == __LOGICAL
                              IF (ppm_recvbufferd(ibuffer) .GT.     &
     &                           (1.0_ppm_kind_double-ppm_myepsd)) THEN
                                 fdata(1,imesh,jmesh,kmesh) = .TRUE.
                              ELSE
                                 fdata(2,imesh,jmesh,kmesh) = .FALSE.
                              ENDIF
                              ibuffer = ibuffer + 1
                              IF (ppm_recvbufferd(ibuffer) .GT.     &
     &                           (1.0_ppm_kind_double-ppm_myepsd)) THEN
                                 fdata(1,imesh,jmesh,kmesh) = .TRUE.
                              ELSE
                                 fdata(2,imesh,jmesh,kmesh) = .FALSE.
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
                              fdata(1,imesh,jmesh,kmesh) =   &
     &                           fdata(1,imesh,jmesh,kmesh)+ &
     &                           REAL(ppm_recvbufferd(ibuffer),ppm_kind_single)
                              ibuffer = ibuffer + 1
                              fdata(2,imesh,jmesh,kmesh) =   &
     &                           fdata(2,imesh,jmesh,kmesh)+ &
     &                           REAL(ppm_recvbufferd(ibuffer),ppm_kind_single)
#elif  __KIND == __DOUBLE_PRECISION
                              fdata(1,imesh,jmesh,kmesh) =   &
     &                           fdata(1,imesh,jmesh,kmesh)+ &
     &                           ppm_recvbufferd(ibuffer)
                              ibuffer = ibuffer + 1
                              fdata(2,imesh,jmesh,kmesh) =   &
     &                           fdata(2,imesh,jmesh,kmesh)+ &
     &                           ppm_recvbufferd(ibuffer)
#elif  __KIND == __SINGLE_PRECISION_COMPLEX
                              fdata(1,imesh,jmesh,kmesh) =   &
     &                           fdata(1,imesh,jmesh,kmesh)+ &
     &                           CMPLX(ppm_recvbufferd(ibuffer-1),  &
     &                           ppm_recvbufferd(ibuffer),ppm_kind_single)
                              ibuffer = ibuffer + 2
                              fdata(2,imesh,jmesh,kmesh) =   &
     &                           fdata(2,imesh,jmesh,kmesh)+ &
     &                           CMPLX(ppm_recvbufferd(ibuffer-1),  &
     &                           ppm_recvbufferd(ibuffer),ppm_kind_single)
#elif  __KIND == __DOUBLE_PRECISION_COMPLEX
                              fdata(1,imesh,jmesh,kmesh) =   &
     &                           fdata(1,imesh,jmesh,kmesh)+ &
     &                           CMPLX(ppm_recvbufferd(ibuffer-1),  &
     &                           ppm_recvbufferd(ibuffer),ppm_kind_double)
                              ibuffer = ibuffer + 2
                              fdata(2,imesh,jmesh,kmesh) =   &
     &                           fdata(2,imesh,jmesh,kmesh)+ &
     &                           CMPLX(ppm_recvbufferd(ibuffer-1),  &
     &                           ppm_recvbufferd(ibuffer),ppm_kind_double)
#elif  __KIND == __INTEGER
                              fdata(1,imesh,jmesh,kmesh) =   &
     &                           fdata(1,imesh,jmesh,kmesh)+ &
     &                           INT(ppm_recvbufferd(ibuffer))
                              ibuffer = ibuffer + 1
                              fdata(2,imesh,jmesh,kmesh) =   &
     &                           fdata(2,imesh,jmesh,kmesh)+ &
     &                           INT(ppm_recvbufferd(ibuffer))
#elif  __KIND == __LOGICAL
                              IF (ppm_recvbufferd(ibuffer) .GT.     &
     &                           (1.0_ppm_kind_double-ppm_myepsd)) THEN
                                 fdata(1,imesh,jmesh,kmesh) =    &
     &                              (fdata(1,imesh,jmesh,kmesh) .AND. .TRUE.)
                              ELSE
                                 fdata(1,imesh,jmesh,kmesh) =    &
     &                              (fdata(1,imesh,jmesh,kmesh) .AND. .FALSE.)
                              ENDIF
                              ibuffer = ibuffer + 1
                              IF (ppm_recvbufferd(ibuffer) .GT.     &
     &                           (1.0_ppm_kind_double-ppm_myepsd)) THEN
                                 fdata(2,imesh,jmesh,kmesh) =    &
     &                              (fdata(2,imesh,jmesh,kmesh) .AND. .TRUE.)
                              ELSE
                                 fdata(2,imesh,jmesh,kmesh) =    &
     &                              (fdata(2,imesh,jmesh,kmesh) .AND. .FALSE.)
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
                              fdata(1,imesh,jmesh,kmesh) =     &
     &                           REAL(ppm_recvbufferd(ibuffer),ppm_kind_single)
                              ibuffer = ibuffer + 1
                              fdata(2,imesh,jmesh,kmesh) =     &
     &                           REAL(ppm_recvbufferd(ibuffer),ppm_kind_single)
                              ibuffer = ibuffer + 1
                              fdata(3,imesh,jmesh,kmesh) =     &
     &                           REAL(ppm_recvbufferd(ibuffer),ppm_kind_single)
#elif  __KIND == __DOUBLE_PRECISION
                              fdata(1,imesh,jmesh,kmesh) =     &
     &                           ppm_recvbufferd(ibuffer)
                              ibuffer = ibuffer + 1
                              fdata(2,imesh,jmesh,kmesh) =     &
     &                           ppm_recvbufferd(ibuffer)
                              ibuffer = ibuffer + 1
                              fdata(3,imesh,jmesh,kmesh) =     &
     &                           ppm_recvbufferd(ibuffer)
#elif  __KIND == __SINGLE_PRECISION_COMPLEX
                              fdata(1,imesh,jmesh,kmesh) =     &
     &                           CMPLX(ppm_recvbufferd(ibuffer-1),  &
     &                           ppm_recvbufferd(ibuffer),ppm_kind_single)
                              ibuffer = ibuffer + 2
                              fdata(2,imesh,jmesh,kmesh) =     &
     &                           CMPLX(ppm_recvbufferd(ibuffer-1),  &
     &                           ppm_recvbufferd(ibuffer),ppm_kind_single)
                              ibuffer = ibuffer + 2
                              fdata(3,imesh,jmesh,kmesh) =     &
     &                           CMPLX(ppm_recvbufferd(ibuffer-1),  &
     &                           ppm_recvbufferd(ibuffer),ppm_kind_single)
#elif  __KIND == __DOUBLE_PRECISION_COMPLEX
                              fdata(1,imesh,jmesh,kmesh) =     &
     &                           CMPLX(ppm_recvbufferd(ibuffer-1),  &
     &                           ppm_recvbufferd(ibuffer),ppm_kind_double)
                              ibuffer = ibuffer + 2
                              fdata(2,imesh,jmesh,kmesh) =     &
     &                           CMPLX(ppm_recvbufferd(ibuffer-1),  &
     &                           ppm_recvbufferd(ibuffer),ppm_kind_double)
                              ibuffer = ibuffer + 2
                              fdata(3,imesh,jmesh,kmesh) =     &
     &                           CMPLX(ppm_recvbufferd(ibuffer-1),  &
     &                           ppm_recvbufferd(ibuffer),ppm_kind_double)
#elif  __KIND == __INTEGER
                              fdata(1,imesh,jmesh,kmesh) =     &
     &                           INT(ppm_recvbufferd(ibuffer))
                              ibuffer = ibuffer + 1
                              fdata(2,imesh,jmesh,kmesh) =     &
     &                           INT(ppm_recvbufferd(ibuffer))
                              ibuffer = ibuffer + 1
                              fdata(3,imesh,jmesh,kmesh) =     &
     &                           INT(ppm_recvbufferd(ibuffer))
#elif  __KIND == __LOGICAL
                              IF (ppm_recvbufferd(ibuffer) .GT.     &
     &                           (1.0_ppm_kind_double-ppm_myepsd)) THEN
                                 fdata(1,imesh,jmesh,kmesh) = .TRUE.
                              ELSE
                                 fdata(1,imesh,jmesh,kmesh) = .FALSE.
                              ENDIF
                              ibuffer = ibuffer + 1
                              IF (ppm_recvbufferd(ibuffer) .GT.     &
     &                           (1.0_ppm_kind_double-ppm_myepsd)) THEN
                                 fdata(2,imesh,jmesh,kmesh) = .TRUE.
                              ELSE
                                 fdata(2,imesh,jmesh,kmesh) = .FALSE.
                              ENDIF
                              ibuffer = ibuffer + 1
                              IF (ppm_recvbufferd(ibuffer) .GT.     &
     &                           (1.0_ppm_kind_double-ppm_myepsd)) THEN
                                 fdata(3,imesh,jmesh,kmesh) = .TRUE.
                              ELSE
                                 fdata(3,imesh,jmesh,kmesh) = .FALSE.
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
                              fdata(1,imesh,jmesh,kmesh) =   &
     &                           fdata(1,imesh,jmesh,kmesh)+ &
     &                           REAL(ppm_recvbufferd(ibuffer),ppm_kind_single)
                              ibuffer = ibuffer + 1
                              fdata(2,imesh,jmesh,kmesh) =   &
     &                           fdata(2,imesh,jmesh,kmesh)+ &
     &                           REAL(ppm_recvbufferd(ibuffer),ppm_kind_single)
                              ibuffer = ibuffer + 1
                              fdata(3,imesh,jmesh,kmesh) =   &
     &                           fdata(3,imesh,jmesh,kmesh)+ &
     &                           REAL(ppm_recvbufferd(ibuffer),ppm_kind_single)
#elif  __KIND == __DOUBLE_PRECISION
                              fdata(1,imesh,jmesh,kmesh) =   &
     &                           fdata(1,imesh,jmesh,kmesh)+ &
     &                           ppm_recvbufferd(ibuffer)
                              ibuffer = ibuffer + 1
                              fdata(2,imesh,jmesh,kmesh) =   &
     &                           fdata(2,imesh,jmesh,kmesh)+ &
     &                           ppm_recvbufferd(ibuffer)
                              ibuffer = ibuffer + 1
                              fdata(3,imesh,jmesh,kmesh) =   &
     &                           fdata(3,imesh,jmesh,kmesh)+ &
     &                           ppm_recvbufferd(ibuffer)
#elif  __KIND == __SINGLE_PRECISION_COMPLEX
                              fdata(1,imesh,jmesh,kmesh) =   &
     &                           fdata(1,imesh,jmesh,kmesh)+ &
     &                           CMPLX(ppm_recvbufferd(ibuffer-1),  &
     &                           ppm_recvbufferd(ibuffer),ppm_kind_single)
                              ibuffer = ibuffer + 2
                              fdata(2,imesh,jmesh,kmesh) =   &
     &                           fdata(2,imesh,jmesh,kmesh)+ &
     &                           CMPLX(ppm_recvbufferd(ibuffer-1),  &
     &                           ppm_recvbufferd(ibuffer),ppm_kind_single)
                              ibuffer = ibuffer + 2
                              fdata(3,imesh,jmesh,kmesh) =   &
     &                           fdata(3,imesh,jmesh,kmesh)+ &
     &                           CMPLX(ppm_recvbufferd(ibuffer-1),  &
     &                           ppm_recvbufferd(ibuffer),ppm_kind_single)
#elif  __KIND == __DOUBLE_PRECISION_COMPLEX
                              fdata(1,imesh,jmesh,kmesh) =   &
     &                           fdata(1,imesh,jmesh,kmesh)+ &
     &                           CMPLX(ppm_recvbufferd(ibuffer-1),  &
     &                           ppm_recvbufferd(ibuffer),ppm_kind_double)
                              ibuffer = ibuffer + 2
                              fdata(2,imesh,jmesh,kmesh) =   &
     &                           fdata(2,imesh,jmesh,kmesh)+ &
     &                           CMPLX(ppm_recvbufferd(ibuffer-1),  &
     &                           ppm_recvbufferd(ibuffer),ppm_kind_double)
                              ibuffer = ibuffer + 2
                              fdata(3,imesh,jmesh,kmesh) =   &
     &                           fdata(3,imesh,jmesh,kmesh)+ &
     &                           CMPLX(ppm_recvbufferd(ibuffer-1),  &
     &                           ppm_recvbufferd(ibuffer),ppm_kind_double)
#elif  __KIND == __INTEGER
                              fdata(1,imesh,jmesh,kmesh) =   &
     &                           fdata(1,imesh,jmesh,kmesh)+ &
     &                           INT(ppm_recvbufferd(ibuffer))
                              ibuffer = ibuffer + 1
                              fdata(2,imesh,jmesh,kmesh) =   &
     &                           fdata(2,imesh,jmesh,kmesh)+ &
     &                           INT(ppm_recvbufferd(ibuffer))
                              ibuffer = ibuffer + 1
                              fdata(3,imesh,jmesh,kmesh) =   &
     &                           fdata(3,imesh,jmesh,kmesh)+ &
     &                           INT(ppm_recvbufferd(ibuffer))
#elif  __KIND == __LOGICAL
                              IF (ppm_recvbufferd(ibuffer) .GT.     &
     &                           (1.0_ppm_kind_double-ppm_myepsd)) THEN
                                 fdata(1,imesh,jmesh,kmesh) =    &
     &                              (fdata(1,imesh,jmesh,kmesh) .AND. .TRUE.)
                              ELSE
                                 fdata(1,imesh,jmesh,kmesh) =    &
     &                              (fdata(1,imesh,jmesh,kmesh) .AND. .FALSE.)
                              ENDIF
                              ibuffer = ibuffer + 1
                              IF (ppm_recvbufferd(ibuffer) .GT.     &
     &                           (1.0_ppm_kind_double-ppm_myepsd)) THEN
                                 fdata(2,imesh,jmesh,kmesh) =    &
     &                              (fdata(2,imesh,jmesh,kmesh) .AND. .TRUE.)
                              ELSE
                                 fdata(2,imesh,jmesh,kmesh) =    &
     &                              (fdata(2,imesh,jmesh,kmesh) .AND. .FALSE.)
                              ENDIF
                              ibuffer = ibuffer + 1
                              IF (ppm_recvbufferd(ibuffer) .GT.     &
     &                           (1.0_ppm_kind_double-ppm_myepsd)) THEN
                                 fdata(3,imesh,jmesh,kmesh) =    &
     &                              (fdata(3,imesh,jmesh,kmesh) .AND. .TRUE.)
                              ELSE
                                 fdata(3,imesh,jmesh,kmesh) =    &
     &                              (fdata(3,imesh,jmesh,kmesh) .AND. .FALSE.)
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
                              fdata(1,imesh,jmesh,kmesh) =     &
     &                           REAL(ppm_recvbufferd(ibuffer),ppm_kind_single)
                              ibuffer = ibuffer + 1
                              fdata(2,imesh,jmesh,kmesh) =     &
     &                           REAL(ppm_recvbufferd(ibuffer),ppm_kind_single)
                              ibuffer = ibuffer + 1
                              fdata(3,imesh,jmesh,kmesh) =     &
     &                           REAL(ppm_recvbufferd(ibuffer),ppm_kind_single)
                              ibuffer = ibuffer + 1
                              fdata(4,imesh,jmesh,kmesh) =     &
     &                           REAL(ppm_recvbufferd(ibuffer),ppm_kind_single)
#elif  __KIND == __DOUBLE_PRECISION
                              fdata(1,imesh,jmesh,kmesh) =     &
     &                           ppm_recvbufferd(ibuffer)
                              ibuffer = ibuffer + 1
                              fdata(2,imesh,jmesh,kmesh) =     &
     &                           ppm_recvbufferd(ibuffer)
                              ibuffer = ibuffer + 1
                              fdata(3,imesh,jmesh,kmesh) =     &
     &                           ppm_recvbufferd(ibuffer)
                              ibuffer = ibuffer + 1
                              fdata(4,imesh,jmesh,kmesh) =     &
     &                           ppm_recvbufferd(ibuffer)
#elif  __KIND == __SINGLE_PRECISION_COMPLEX
                              fdata(1,imesh,jmesh,kmesh) =     &
     &                           CMPLX(ppm_recvbufferd(ibuffer-1),  &
     &                           ppm_recvbufferd(ibuffer),ppm_kind_single)
                              ibuffer = ibuffer + 2
                              fdata(2,imesh,jmesh,kmesh) =     &
     &                           CMPLX(ppm_recvbufferd(ibuffer-1),  &
     &                           ppm_recvbufferd(ibuffer),ppm_kind_single)
                              ibuffer = ibuffer + 2
                              fdata(3,imesh,jmesh,kmesh) =     &
     &                           CMPLX(ppm_recvbufferd(ibuffer-1),  &
     &                           ppm_recvbufferd(ibuffer),ppm_kind_single)
                              ibuffer = ibuffer + 2
                              fdata(4,imesh,jmesh,kmesh) =     &
     &                           CMPLX(ppm_recvbufferd(ibuffer-1),  &
     &                           ppm_recvbufferd(ibuffer),ppm_kind_single)
#elif  __KIND == __DOUBLE_PRECISION_COMPLEX
                              fdata(1,imesh,jmesh,kmesh) =     &
     &                           CMPLX(ppm_recvbufferd(ibuffer-1),  &
     &                           ppm_recvbufferd(ibuffer),ppm_kind_double)
                              ibuffer = ibuffer + 2
                              fdata(2,imesh,jmesh,kmesh) =     &
     &                           CMPLX(ppm_recvbufferd(ibuffer-1),  &
     &                           ppm_recvbufferd(ibuffer),ppm_kind_double)
                              ibuffer = ibuffer + 2
                              fdata(3,imesh,jmesh,kmesh) =     &
     &                           CMPLX(ppm_recvbufferd(ibuffer-1),  &
     &                           ppm_recvbufferd(ibuffer),ppm_kind_double)
                              ibuffer = ibuffer + 2
                              fdata(4,imesh,jmesh,kmesh) =     &
     &                           CMPLX(ppm_recvbufferd(ibuffer-1),  &
     &                           ppm_recvbufferd(ibuffer),ppm_kind_double)
#elif  __KIND == __INTEGER
                              fdata(1,imesh,jmesh,kmesh) =     &
     &                           INT(ppm_recvbufferd(ibuffer))
                              ibuffer = ibuffer + 1
                              fdata(2,imesh,jmesh,kmesh) =     &
     &                           INT(ppm_recvbufferd(ibuffer))
                              ibuffer = ibuffer + 1
                              fdata(3,imesh,jmesh,kmesh) =     &
     &                           INT(ppm_recvbufferd(ibuffer))
                              ibuffer = ibuffer + 1
                              fdata(4,imesh,jmesh,kmesh) =     &
     &                           INT(ppm_recvbufferd(ibuffer))
#elif  __KIND == __LOGICAL
                              IF (ppm_recvbufferd(ibuffer) .GT.     &
     &                           (1.0_ppm_kind_double-ppm_myepsd)) THEN
                                 fdata(1,imesh,jmesh,kmesh) = .TRUE.
                              ELSE
                                 fdata(1,imesh,jmesh,kmesh) = .FALSE.
                              ENDIF
                              ibuffer = ibuffer + 1
                              IF (ppm_recvbufferd(ibuffer) .GT.     &
     &                           (1.0_ppm_kind_double-ppm_myepsd)) THEN
                                 fdata(2,imesh,jmesh,kmesh) = .TRUE.
                              ELSE
                                 fdata(2,imesh,jmesh,kmesh) = .FALSE.
                              ENDIF
                              ibuffer = ibuffer + 1
                              IF (ppm_recvbufferd(ibuffer) .GT.     &
     &                           (1.0_ppm_kind_double-ppm_myepsd)) THEN
                                 fdata(3,imesh,jmesh,kmesh) = .TRUE.
                              ELSE
                                 fdata(3,imesh,jmesh,kmesh) = .FALSE.
                              ENDIF
                              ibuffer = ibuffer + 1
                              IF (ppm_recvbufferd(ibuffer) .GT.     &
     &                           (1.0_ppm_kind_double-ppm_myepsd)) THEN
                                 fdata(4,imesh,jmesh,kmesh) = .TRUE.
                              ELSE
                                 fdata(4,imesh,jmesh,kmesh) = .FALSE.
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
                              fdata(1,imesh,jmesh,kmesh) =   &
     &                           fdata(1,imesh,jmesh,kmesh)+ &
     &                           REAL(ppm_recvbufferd(ibuffer),ppm_kind_single)
                              ibuffer = ibuffer + 1
                              fdata(2,imesh,jmesh,kmesh) =   &
     &                           fdata(2,imesh,jmesh,kmesh)+ &
     &                           REAL(ppm_recvbufferd(ibuffer),ppm_kind_single)
                              ibuffer = ibuffer + 1
                              fdata(3,imesh,jmesh,kmesh) =   &
     &                           fdata(3,imesh,jmesh,kmesh)+ &
     &                           REAL(ppm_recvbufferd(ibuffer),ppm_kind_single)
                              ibuffer = ibuffer + 1
                              fdata(4,imesh,jmesh,kmesh) =   &
     &                           fdata(4,imesh,jmesh,kmesh)+ &
     &                           REAL(ppm_recvbufferd(ibuffer),ppm_kind_single)
#elif  __KIND == __DOUBLE_PRECISION
                              fdata(1,imesh,jmesh,kmesh) =   &
     &                           fdata(1,imesh,jmesh,kmesh)+ &
     &                           ppm_recvbufferd(ibuffer)
                              ibuffer = ibuffer + 1
                              fdata(2,imesh,jmesh,kmesh) =   &
     &                           fdata(2,imesh,jmesh,kmesh)+ &
     &                           ppm_recvbufferd(ibuffer)
                              ibuffer = ibuffer + 1
                              fdata(3,imesh,jmesh,kmesh) =   &
     &                           fdata(3,imesh,jmesh,kmesh)+ &
     &                           ppm_recvbufferd(ibuffer)
                              ibuffer = ibuffer + 1
                              fdata(4,imesh,jmesh,kmesh) =   &
     &                           fdata(4,imesh,jmesh,kmesh)+ &
     &                           ppm_recvbufferd(ibuffer)
#elif  __KIND == __SINGLE_PRECISION_COMPLEX
                              fdata(1,imesh,jmesh,kmesh) =   &
     &                           fdata(1,imesh,jmesh,kmesh)+ &
     &                           CMPLX(ppm_recvbufferd(ibuffer-1),  &
     &                           ppm_recvbufferd(ibuffer),ppm_kind_single)
                              ibuffer = ibuffer + 2
                              fdata(2,imesh,jmesh,kmesh) =   &
     &                           fdata(2,imesh,jmesh,kmesh)+ &
     &                           CMPLX(ppm_recvbufferd(ibuffer-1),  &
     &                           ppm_recvbufferd(ibuffer),ppm_kind_single)
                              ibuffer = ibuffer + 2
                              fdata(3,imesh,jmesh,kmesh) =   &
     &                           fdata(3,imesh,jmesh,kmesh)+ &
     &                           CMPLX(ppm_recvbufferd(ibuffer-1),  &
     &                           ppm_recvbufferd(ibuffer),ppm_kind_single)
                              ibuffer = ibuffer + 2
                              fdata(4,imesh,jmesh,kmesh) =   &
     &                           fdata(4,imesh,jmesh,kmesh)+ &
     &                           CMPLX(ppm_recvbufferd(ibuffer-1),  &
     &                           ppm_recvbufferd(ibuffer),ppm_kind_single)
#elif  __KIND == __DOUBLE_PRECISION_COMPLEX
                              fdata(1,imesh,jmesh,kmesh) =   &
     &                           fdata(1,imesh,jmesh,kmesh)+ &
     &                           CMPLX(ppm_recvbufferd(ibuffer-1),  &
     &                           ppm_recvbufferd(ibuffer),ppm_kind_double)
                              ibuffer = ibuffer + 2
                              fdata(2,imesh,jmesh,kmesh) =   &
     &                           fdata(2,imesh,jmesh,kmesh)+ &
     &                           CMPLX(ppm_recvbufferd(ibuffer-1),  &
     &                           ppm_recvbufferd(ibuffer),ppm_kind_double)
                              ibuffer = ibuffer + 2
                              fdata(3,imesh,jmesh,kmesh) =   &
     &                           fdata(3,imesh,jmesh,kmesh)+ &
     &                           CMPLX(ppm_recvbufferd(ibuffer-1),  &
     &                           ppm_recvbufferd(ibuffer),ppm_kind_double)
                              ibuffer = ibuffer + 2
                              fdata(4,imesh,jmesh,kmesh) =   &
     &                           fdata(4,imesh,jmesh,kmesh)+ &
     &                           CMPLX(ppm_recvbufferd(ibuffer-1),  &
     &                           ppm_recvbufferd(ibuffer),ppm_kind_double)
#elif  __KIND == __INTEGER
                              fdata(1,imesh,jmesh,kmesh) =   &
     &                           fdata(1,imesh,jmesh,kmesh)+ &
     &                           INT(ppm_recvbufferd(ibuffer))
                              ibuffer = ibuffer + 1
                              fdata(2,imesh,jmesh,kmesh) =   &
     &                           fdata(2,imesh,jmesh,kmesh)+ &
     &                           INT(ppm_recvbufferd(ibuffer))
                              ibuffer = ibuffer + 1
                              fdata(3,imesh,jmesh,kmesh) =   &
     &                           fdata(3,imesh,jmesh,kmesh)+ &
     &                           INT(ppm_recvbufferd(ibuffer))
                              ibuffer = ibuffer + 1
                              fdata(4,imesh,jmesh,kmesh) =   &
     &                           fdata(4,imesh,jmesh,kmesh)+ &
     &                           INT(ppm_recvbufferd(ibuffer))
#elif  __KIND == __LOGICAL
                              IF (ppm_recvbufferd(ibuffer) .GT.     &
     &                           (1.0_ppm_kind_double-ppm_myepsd)) THEN
                                 fdata(1,imesh,jmesh,kmesh) =    &
     &                              (fdata(1,imesh,jmesh,kmesh) .AND. .TRUE.)
                              ELSE
                                 fdata(1,imesh,jmesh,kmesh) =    &
     &                              (fdata(1,imesh,jmesh,kmesh) .AND. .FALSE.)
                              ENDIF
                              ibuffer = ibuffer + 1
                              IF (ppm_recvbufferd(ibuffer) .GT.     &
     &                           (1.0_ppm_kind_double-ppm_myepsd)) THEN
                                 fdata(2,imesh,jmesh,kmesh) =    &
     &                              (fdata(2,imesh,jmesh,kmesh) .AND. .TRUE.)
                              ELSE
                                 fdata(2,imesh,jmesh,kmesh) =    &
     &                              (fdata(2,imesh,jmesh,kmesh) .AND. .FALSE.)
                              ENDIF
                              ibuffer = ibuffer + 1
                              IF (ppm_recvbufferd(ibuffer) .GT.     &
     &                           (1.0_ppm_kind_double-ppm_myepsd)) THEN
                                 fdata(3,imesh,jmesh,kmesh) =    &
     &                              (fdata(3,imesh,jmesh,kmesh) .AND. .TRUE.)
                              ELSE
                                 fdata(3,imesh,jmesh,kmesh) =    &
     &                              (fdata(3,imesh,jmesh,kmesh) .AND. .FALSE.)
                              ENDIF
                              ibuffer = ibuffer + 1
                              IF (ppm_recvbufferd(ibuffer) .GT.     &
     &                           (1.0_ppm_kind_double-ppm_myepsd)) THEN
                                 fdata(4,imesh,jmesh,kmesh) =    &
     &                              (fdata(4,imesh,jmesh,kmesh) .AND. .TRUE.)
                              ELSE
                                 fdata(4,imesh,jmesh,kmesh) =    &
     &                              (fdata(4,imesh,jmesh,kmesh) .AND. .FALSE.)
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
                              fdata(1,imesh,jmesh,kmesh) =     &
     &                           REAL(ppm_recvbufferd(ibuffer),ppm_kind_single)
                              ibuffer = ibuffer + 1
                              fdata(2,imesh,jmesh,kmesh) =     &
     &                           REAL(ppm_recvbufferd(ibuffer),ppm_kind_single)
                              ibuffer = ibuffer + 1
                              fdata(3,imesh,jmesh,kmesh) =     &
     &                           REAL(ppm_recvbufferd(ibuffer),ppm_kind_single)
                              ibuffer = ibuffer + 1
                              fdata(4,imesh,jmesh,kmesh) =     &
     &                           REAL(ppm_recvbufferd(ibuffer),ppm_kind_single)
                              ibuffer = ibuffer + 1
                              fdata(5,imesh,jmesh,kmesh) =     &
     &                           REAL(ppm_recvbufferd(ibuffer),ppm_kind_single)
#elif  __KIND == __DOUBLE_PRECISION
                              fdata(1,imesh,jmesh,kmesh) =     &
     &                           ppm_recvbufferd(ibuffer)
                              ibuffer = ibuffer + 1
                              fdata(2,imesh,jmesh,kmesh) =     &
     &                           ppm_recvbufferd(ibuffer)
                              ibuffer = ibuffer + 1
                              fdata(3,imesh,jmesh,kmesh) =     &
     &                           ppm_recvbufferd(ibuffer)
                              ibuffer = ibuffer + 1
                              fdata(4,imesh,jmesh,kmesh) =     &
     &                           ppm_recvbufferd(ibuffer)
                              ibuffer = ibuffer + 1
                              fdata(5,imesh,jmesh,kmesh) =     &
     &                           ppm_recvbufferd(ibuffer)
#elif  __KIND == __SINGLE_PRECISION_COMPLEX
                              fdata(1,imesh,jmesh,kmesh) =     &
     &                           CMPLX(ppm_recvbufferd(ibuffer-1),  &
     &                           ppm_recvbufferd(ibuffer),ppm_kind_single)
                              ibuffer = ibuffer + 2
                              fdata(2,imesh,jmesh,kmesh) =     &
     &                           CMPLX(ppm_recvbufferd(ibuffer-1),  &
     &                           ppm_recvbufferd(ibuffer),ppm_kind_single)
                              ibuffer = ibuffer + 2
                              fdata(3,imesh,jmesh,kmesh) =     &
     &                           CMPLX(ppm_recvbufferd(ibuffer-1),  &
     &                           ppm_recvbufferd(ibuffer),ppm_kind_single)
                              ibuffer = ibuffer + 2
                              fdata(4,imesh,jmesh,kmesh) =     &
     &                           CMPLX(ppm_recvbufferd(ibuffer-1),  &
     &                           ppm_recvbufferd(ibuffer),ppm_kind_single)
                              ibuffer = ibuffer + 2
                              fdata(5,imesh,jmesh,kmesh) =     &
     &                           CMPLX(ppm_recvbufferd(ibuffer-1),  &
     &                           ppm_recvbufferd(ibuffer),ppm_kind_single)
#elif  __KIND == __DOUBLE_PRECISION_COMPLEX
                              fdata(1,imesh,jmesh,kmesh) =     &
     &                           CMPLX(ppm_recvbufferd(ibuffer-1),  &
     &                           ppm_recvbufferd(ibuffer),ppm_kind_double)
                              ibuffer = ibuffer + 2
                              fdata(2,imesh,jmesh,kmesh) =     &
     &                           CMPLX(ppm_recvbufferd(ibuffer-1),  &
     &                           ppm_recvbufferd(ibuffer),ppm_kind_double)
                              ibuffer = ibuffer + 2
                              fdata(3,imesh,jmesh,kmesh) =     &
     &                           CMPLX(ppm_recvbufferd(ibuffer-1),  &
     &                           ppm_recvbufferd(ibuffer),ppm_kind_double)
                              ibuffer = ibuffer + 2
                              fdata(4,imesh,jmesh,kmesh) =     &
     &                           CMPLX(ppm_recvbufferd(ibuffer-1),  &
     &                           ppm_recvbufferd(ibuffer),ppm_kind_double)
                              ibuffer = ibuffer + 2
                              fdata(5,imesh,jmesh,kmesh) =     &
     &                           CMPLX(ppm_recvbufferd(ibuffer-1),  &
     &                           ppm_recvbufferd(ibuffer),ppm_kind_double)
#elif  __KIND == __INTEGER
                              fdata(1,imesh,jmesh,kmesh) =     &
     &                           INT(ppm_recvbufferd(ibuffer))
                              ibuffer = ibuffer + 1
                              fdata(2,imesh,jmesh,kmesh) =     &
     &                           INT(ppm_recvbufferd(ibuffer))
                              ibuffer = ibuffer + 1
                              fdata(3,imesh,jmesh,kmesh) =     &
     &                           INT(ppm_recvbufferd(ibuffer))
                              ibuffer = ibuffer + 1
                              fdata(4,imesh,jmesh,kmesh) =     &
     &                           INT(ppm_recvbufferd(ibuffer))
                              ibuffer = ibuffer + 1
                              fdata(5,imesh,jmesh,kmesh) =     &
     &                           INT(ppm_recvbufferd(ibuffer))
#elif  __KIND == __LOGICAL
                              IF (ppm_recvbufferd(ibuffer) .GT.     &
     &                           (1.0_ppm_kind_double-ppm_myepsd)) THEN
                                 fdata(1,imesh,jmesh,kmesh) = .TRUE.
                              ELSE
                                 fdata(1,imesh,jmesh,kmesh) = .FALSE.
                              ENDIF
                              ibuffer = ibuffer + 1
                              IF (ppm_recvbufferd(ibuffer) .GT.     &
     &                           (1.0_ppm_kind_double-ppm_myepsd)) THEN
                                 fdata(2,imesh,jmesh,kmesh) = .TRUE.
                              ELSE
                                 fdata(2,imesh,jmesh,kmesh) = .FALSE.
                              ENDIF
                              ibuffer = ibuffer + 1
                              IF (ppm_recvbufferd(ibuffer) .GT.     &
     &                           (1.0_ppm_kind_double-ppm_myepsd)) THEN
                                 fdata(3,imesh,jmesh,kmesh) = .TRUE.
                              ELSE
                                 fdata(3,imesh,jmesh,kmesh) = .FALSE.
                              ENDIF
                              ibuffer = ibuffer + 1
                              IF (ppm_recvbufferd(ibuffer) .GT.     &
     &                           (1.0_ppm_kind_double-ppm_myepsd)) THEN
                                 fdata(4,imesh,jmesh,kmesh) = .TRUE.
                              ELSE
                                 fdata(4,imesh,jmesh,kmesh) = .FALSE.
                              ENDIF
                              ibuffer = ibuffer + 1
                              IF (ppm_recvbufferd(ibuffer) .GT.     &
     &                           (1.0_ppm_kind_double-ppm_myepsd)) THEN
                                 fdata(5,imesh,jmesh,kmesh) = .TRUE.
                              ELSE
                                 fdata(5,imesh,jmesh,kmesh) = .FALSE.
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
                              fdata(1,imesh,jmesh,kmesh) =   &
     &                           fdata(1,imesh,jmesh,kmesh)+ &
     &                           REAL(ppm_recvbufferd(ibuffer),ppm_kind_single)
                              ibuffer = ibuffer + 1
                              fdata(2,imesh,jmesh,kmesh) =   &
     &                           fdata(2,imesh,jmesh,kmesh)+ &
     &                           REAL(ppm_recvbufferd(ibuffer),ppm_kind_single)
                              ibuffer = ibuffer + 1
                              fdata(3,imesh,jmesh,kmesh) =   &
     &                           fdata(3,imesh,jmesh,kmesh)+ &
     &                           REAL(ppm_recvbufferd(ibuffer),ppm_kind_single)
                              ibuffer = ibuffer + 1
                              fdata(4,imesh,jmesh,kmesh) =   &
     &                           fdata(4,imesh,jmesh,kmesh)+ &
     &                           REAL(ppm_recvbufferd(ibuffer),ppm_kind_single)
                              ibuffer = ibuffer + 1
                              fdata(5,imesh,jmesh,kmesh) =   &
     &                           fdata(5,imesh,jmesh,kmesh)+ &
     &                           REAL(ppm_recvbufferd(ibuffer),ppm_kind_single)
#elif  __KIND == __DOUBLE_PRECISION
                              fdata(1,imesh,jmesh,kmesh) =   &
     &                           fdata(1,imesh,jmesh,kmesh)+ &
     &                           ppm_recvbufferd(ibuffer)
                              ibuffer = ibuffer + 1
                              fdata(2,imesh,jmesh,kmesh) =   &
     &                           fdata(2,imesh,jmesh,kmesh)+ &
     &                           ppm_recvbufferd(ibuffer)
                              ibuffer = ibuffer + 1
                              fdata(3,imesh,jmesh,kmesh) =   &
     &                           fdata(3,imesh,jmesh,kmesh)+ &
     &                           ppm_recvbufferd(ibuffer)
                              ibuffer = ibuffer + 1
                              fdata(4,imesh,jmesh,kmesh) =   &
     &                           fdata(4,imesh,jmesh,kmesh)+ &
     &                           ppm_recvbufferd(ibuffer)
                              ibuffer = ibuffer + 1
                              fdata(5,imesh,jmesh,kmesh) =   &
     &                           fdata(5,imesh,jmesh,kmesh)+ &
     &                           ppm_recvbufferd(ibuffer)
#elif  __KIND == __SINGLE_PRECISION_COMPLEX
                              fdata(1,imesh,jmesh,kmesh) =   &
     &                           fdata(1,imesh,jmesh,kmesh)+ &
     &                           CMPLX(ppm_recvbufferd(ibuffer-1),  &
     &                           ppm_recvbufferd(ibuffer),ppm_kind_single)
                              ibuffer = ibuffer + 2
                              fdata(2,imesh,jmesh,kmesh) =   &
     &                           fdata(2,imesh,jmesh,kmesh)+ &
     &                           CMPLX(ppm_recvbufferd(ibuffer-1),  &
     &                           ppm_recvbufferd(ibuffer),ppm_kind_single)
                              ibuffer = ibuffer + 2
                              fdata(3,imesh,jmesh,kmesh) =   &
     &                           fdata(3,imesh,jmesh,kmesh)+ &
     &                           CMPLX(ppm_recvbufferd(ibuffer-1),  &
     &                           ppm_recvbufferd(ibuffer),ppm_kind_single)
                              ibuffer = ibuffer + 2
                              fdata(4,imesh,jmesh,kmesh) =   &
     &                           fdata(4,imesh,jmesh,kmesh)+ &
     &                           CMPLX(ppm_recvbufferd(ibuffer-1),  &
     &                           ppm_recvbufferd(ibuffer),ppm_kind_single)
                              ibuffer = ibuffer + 2
                              fdata(5,imesh,jmesh,kmesh) =   &
     &                           fdata(5,imesh,jmesh,kmesh)+ &
     &                           CMPLX(ppm_recvbufferd(ibuffer-1),  &
     &                           ppm_recvbufferd(ibuffer),ppm_kind_single)
#elif  __KIND == __DOUBLE_PRECISION_COMPLEX
                              fdata(1,imesh,jmesh,kmesh) =   &
     &                           fdata(1,imesh,jmesh,kmesh)+ &
     &                           CMPLX(ppm_recvbufferd(ibuffer-1),  &
     &                           ppm_recvbufferd(ibuffer),ppm_kind_double)
                              ibuffer = ibuffer + 2
                              fdata(2,imesh,jmesh,kmesh) =   &
     &                           fdata(2,imesh,jmesh,kmesh)+ &
     &                           CMPLX(ppm_recvbufferd(ibuffer-1),  &
     &                           ppm_recvbufferd(ibuffer),ppm_kind_double)
                              ibuffer = ibuffer + 2
                              fdata(3,imesh,jmesh,kmesh) =   &
     &                           fdata(3,imesh,jmesh,kmesh)+ &
     &                           CMPLX(ppm_recvbufferd(ibuffer-1),  &
     &                           ppm_recvbufferd(ibuffer),ppm_kind_double)
                              ibuffer = ibuffer + 2
                              fdata(4,imesh,jmesh,kmesh) =   &
     &                           fdata(4,imesh,jmesh,kmesh)+ &
     &                           CMPLX(ppm_recvbufferd(ibuffer-1),  &
     &                           ppm_recvbufferd(ibuffer),ppm_kind_double)
                              ibuffer = ibuffer + 2
                              fdata(5,imesh,jmesh,kmesh) =   &
     &                           fdata(5,imesh,jmesh,kmesh)+ &
     &                           CMPLX(ppm_recvbufferd(ibuffer-1),  &
     &                           ppm_recvbufferd(ibuffer),ppm_kind_double)
#elif  __KIND == __INTEGER
                              fdata(1,imesh,jmesh,kmesh) =   &
     &                           fdata(1,imesh,jmesh,kmesh)+ &
     &                           INT(ppm_recvbufferd(ibuffer))
                              ibuffer = ibuffer + 1
                              fdata(2,imesh,jmesh,kmesh) =   &
     &                           fdata(2,imesh,jmesh,kmesh)+ &
     &                           INT(ppm_recvbufferd(ibuffer))
                              ibuffer = ibuffer + 1
                              fdata(3,imesh,jmesh,kmesh) =   &
     &                           fdata(3,imesh,jmesh,kmesh)+ &
     &                           INT(ppm_recvbufferd(ibuffer))
                              ibuffer = ibuffer + 1
                              fdata(4,imesh,jmesh,kmesh) =   &
     &                           fdata(4,imesh,jmesh,kmesh)+ &
     &                           INT(ppm_recvbufferd(ibuffer))
                              ibuffer = ibuffer + 1
                              fdata(5,imesh,jmesh,kmesh) =   &
     &                           fdata(5,imesh,jmesh,kmesh)+ &
     &                           INT(ppm_recvbufferd(ibuffer))
#elif  __KIND == __LOGICAL
                              IF (ppm_recvbufferd(ibuffer) .GT.     &
     &                           (1.0_ppm_kind_double-ppm_myepsd)) THEN
                                 fdata(1,imesh,jmesh,kmesh) =    &
     &                              (fdata(1,imesh,jmesh,kmesh) .AND. .TRUE.)
                              ELSE
                                 fdata(1,imesh,jmesh,kmesh) =    &
     &                              (fdata(1,imesh,jmesh,kmesh) .AND. .FALSE.)
                              ENDIF
                              ibuffer = ibuffer + 1
                              IF (ppm_recvbufferd(ibuffer) .GT.     &
     &                           (1.0_ppm_kind_double-ppm_myepsd)) THEN
                                 fdata(2,imesh,jmesh,kmesh) =    &
     &                              (fdata(2,imesh,jmesh,kmesh) .AND. .TRUE.)
                              ELSE
                                 fdata(2,imesh,jmesh,kmesh) =    &
     &                              (fdata(2,imesh,jmesh,kmesh) .AND. .FALSE.)
                              ENDIF
                              ibuffer = ibuffer + 1
                              IF (ppm_recvbufferd(ibuffer) .GT.     &
     &                           (1.0_ppm_kind_double-ppm_myepsd)) THEN
                                 fdata(3,imesh,jmesh,kmesh) =    &
     &                              (fdata(3,imesh,jmesh,kmesh) .AND. .TRUE.)
                              ELSE
                                 fdata(3,imesh,jmesh,kmesh) =    &
     &                              (fdata(3,imesh,jmesh,kmesh) .AND. .FALSE.)
                              ENDIF
                              ibuffer = ibuffer + 1
                              IF (ppm_recvbufferd(ibuffer) .GT.     &
     &                           (1.0_ppm_kind_double-ppm_myepsd)) THEN
                                 fdata(4,imesh,jmesh,kmesh) =    &
     &                              (fdata(4,imesh,jmesh,kmesh) .AND. .TRUE.)
                              ELSE
                                 fdata(4,imesh,jmesh,kmesh) =    &
     &                              (fdata(4,imesh,jmesh,kmesh) .AND. .FALSE.)
                              ENDIF
                              ibuffer = ibuffer + 1
                              IF (ppm_recvbufferd(ibuffer) .GT.     &
     &                           (1.0_ppm_kind_double-ppm_myepsd)) THEN
                                 fdata(5,imesh,jmesh,kmesh) =    &
     &                              (fdata(5,imesh,jmesh,kmesh) .AND. .TRUE.)
                              ELSE
                                 fdata(5,imesh,jmesh,kmesh) =    &
     &                              (fdata(5,imesh,jmesh,kmesh) .AND. .FALSE.)
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
                                 fdata(k,imesh,jmesh,kmesh) =     &
     &                              REAL(ppm_recvbufferd(ibuffer),ppm_kind_single)
#elif  __KIND == __DOUBLE_PRECISION
                                 fdata(k,imesh,jmesh,kmesh) =     &
     &                              ppm_recvbufferd(ibuffer)
#elif  __KIND == __SINGLE_PRECISION_COMPLEX
                                 fdata(k,imesh,jmesh,kmesh) =     &
     &                              CMPLX(ppm_recvbufferd(ibuffer-1),  &
     &                              ppm_recvbufferd(ibuffer),ppm_kind_single)
#elif  __KIND == __DOUBLE_PRECISION_COMPLEX
                                 fdata(k,imesh,jmesh,kmesh) =     &
     &                              CMPLX(ppm_recvbufferd(ibuffer-1),  &
     &                              ppm_recvbufferd(ibuffer),ppm_kind_double)
#elif  __KIND == __INTEGER
                                 fdata(k,imesh,jmesh,kmesh) =     &
     &                              INT(ppm_recvbufferd(ibuffer))
#elif  __KIND == __LOGICAL
                                 IF (ppm_recvbufferd(ibuffer) .GT.     &
     &                              (1.0_ppm_kind_double-ppm_myepsd)) THEN
                                    fdata(k,imesh,jmesh,kmesh) = .TRUE.
                                 ELSE
                                    fdata(k,imesh,jmesh,kmesh) = .FALSE.
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
                                 fdata(k,imesh,jmesh,kmesh) =   &
     &                              fdata(k,imesh,jmesh,kmesh)+ &
     &                              REAL(ppm_recvbufferd(ibuffer),ppm_kind_single)
#elif  __KIND == __DOUBLE_PRECISION
                                 fdata(k,imesh,jmesh,kmesh) =   &
     &                              fdata(k,imesh,jmesh,kmesh)+ &
     &                              ppm_recvbufferd(ibuffer)
#elif  __KIND == __SINGLE_PRECISION_COMPLEX
                                 fdata(k,imesh,jmesh,kmesh) =   &
     &                              fdata(k,imesh,jmesh,kmesh)+ &
     &                              CMPLX(ppm_recvbufferd(ibuffer-1),  &
     &                              ppm_recvbufferd(ibuffer),ppm_kind_single)
#elif  __KIND == __DOUBLE_PRECISION_COMPLEX
                                 fdata(k,imesh,jmesh,kmesh) =   &
     &                              fdata(k,imesh,jmesh,kmesh)+ &
     &                              CMPLX(ppm_recvbufferd(ibuffer-1),  &
     &                              ppm_recvbufferd(ibuffer),ppm_kind_double)
#elif  __KIND == __INTEGER
                                 fdata(k,imesh,jmesh,kmesh) =   &
     &                              fdata(k,imesh,jmesh,kmesh)+ &
     &                              INT(ppm_recvbufferd(ibuffer))
#elif  __KIND == __LOGICAL
                                 IF (ppm_recvbufferd(ibuffer) .GT.     &
     &                              (1.0_ppm_kind_double-ppm_myepsd)) THEN
                                    fdata(k,imesh,jmesh,kmesh) =    &
     &                                 (fdata(k,imesh,jmesh,kmesh) .AND. .TRUE.)
                                 ELSE
                                    fdata(k,imesh,jmesh,kmesh) =    &
     &                                 (fdata(k,imesh,jmesh,kmesh) .AND. .FALSE.)
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
                        IF (mask(imesh,jmesh,isub)) THEN
                           DO k=1,edim
#if    __KIND == __SINGLE_PRECISION_COMPLEX | __KIND == __DOUBLE_PRECISION_COMPLEX
                              ibuffer = ibuffer + 2
#else
                              ibuffer = ibuffer + 1
#endif
                              IF (rtype .EQ. ppm_param_pop_replace) THEN
#if    __KIND == __SINGLE_PRECISION
                                 fdata(k,imesh,jmesh,kmesh) =     &
     &                              REAL(ppm_recvbufferd(ibuffer),ppm_kind_single)
#elif  __KIND == __DOUBLE_PRECISION
                                 fdata(k,imesh,jmesh,kmesh) =     &
     &                              ppm_recvbufferd(ibuffer)
#elif  __KIND == __SINGLE_PRECISION_COMPLEX
                                 fdata(k,imesh,jmesh,kmesh) =     &
     &                              CMPLX(ppm_recvbufferd(ibuffer-1),  &
     &                              ppm_recvbufferd(ibuffer),ppm_kind_single)
#elif  __KIND == __DOUBLE_PRECISION_COMPLEX
                                 fdata(k,imesh,jmesh,kmesh) =     &
     &                              CMPLX(ppm_recvbufferd(ibuffer-1),  &
     &                              ppm_recvbufferd(ibuffer),ppm_kind_double)
#elif  __KIND == __INTEGER
                                 fdata(k,imesh,jmesh,kmesh) =     &
     &                              INT(ppm_recvbufferd(ibuffer))
#elif  __KIND == __LOGICAL
                                 IF (ppm_recvbufferd(ibuffer) .GT.     &
     &                              (1.0_ppm_kind_double-ppm_myepsd)) THEN
                                    fdata(k,imesh,jmesh,kmesh) = .TRUE.
                                 ELSE
                                    fdata(k,imesh,jmesh,kmesh) = .FALSE.
                                 ENDIF
#endif
                              ELSEIF (rtype .EQ. ppm_param_pop_add) THEN
#if    __KIND == __SINGLE_PRECISION
                                 fdata(k,imesh,jmesh,kmesh) =   &
     &                              fdata(k,imesh,jmesh,kmesh)+ &
     &                              REAL(ppm_recvbufferd(ibuffer),ppm_kind_single)
#elif  __KIND == __DOUBLE_PRECISION
                                 fdata(k,imesh,jmesh,kmesh) =   &
     &                              fdata(k,imesh,jmesh,kmesh)+ &
     &                              ppm_recvbufferd(ibuffer)
#elif  __KIND == __SINGLE_PRECISION_COMPLEX
                                 fdata(k,imesh,jmesh,kmesh) =   &
     &                              fdata(k,imesh,jmesh,kmesh)+ &
     &                              CMPLX(ppm_recvbufferd(ibuffer-1),  &
     &                              ppm_recvbufferd(ibuffer),ppm_kind_single)
#elif  __KIND == __DOUBLE_PRECISION_COMPLEX
                                 fdata(k,imesh,jmesh,kmesh) =   &
     &                              fdata(k,imesh,jmesh,kmesh)+ &
     &                              CMPLX(ppm_recvbufferd(ibuffer-1),  &
     &                              ppm_recvbufferd(ibuffer),ppm_kind_double)
#elif  __KIND == __INTEGER
                                 fdata(k,imesh,jmesh,kmesh) =   &
     &                              fdata(k,imesh,jmesh,kmesh)+ &
     &                              INT(ppm_recvbufferd(ibuffer))
#elif  __KIND == __LOGICAL
                                 IF (ppm_recvbufferd(ibuffer) .GT.     &
     &                              (1.0_ppm_kind_double-ppm_myepsd)) THEN
                                    fdata(k,imesh,jmesh,kmesh) =    &
     &                                 (fdata(k,imesh,jmesh,kmesh) .AND. .TRUE.)
                                 ELSE
                                    fdata(k,imesh,jmesh,kmesh) =    &
     &                                 (fdata(k,imesh,jmesh,kmesh) .AND. .FALSE.)
                                 ENDIF
#endif
                              ENDIF
                           ENDDO
                        ELSE
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
                           fdata(imesh,jmesh,kmesh) =     &
     &                        REAL(ppm_recvbufferd(ibuffer),ppm_kind_single)
#elif  __KIND == __DOUBLE_PRECISION
                           fdata(imesh,jmesh,kmesh) = ppm_recvbufferd(ibuffer)
#elif  __KIND == __SINGLE_PRECISION_COMPLEX
                           fdata(imesh,jmesh,kmesh) =     &
     &                        CMPLX(ppm_recvbufferd(ibuffer-1),  &
     &                        ppm_recvbufferd(ibuffer),ppm_kind_single)
#elif  __KIND == __DOUBLE_PRECISION_COMPLEX
                           fdata(imesh,jmesh,kmesh) =     &
     &                        CMPLX(ppm_recvbufferd(ibuffer-1),  &
     &                        ppm_recvbufferd(ibuffer),ppm_kind_double)
#elif  __KIND == __INTEGER
                           fdata(imesh,jmesh,kmesh) =     &
     &                        INT(ppm_recvbufferd(ibuffer))
#elif  __KIND == __LOGICAL
                           IF (ppm_recvbufferd(ibuffer) .GT.     &
     &                        (1.0_ppm_kind_double-ppm_myepsd)) THEN
                              fdata(imesh,jmesh,kmesh) = .TRUE.
                           ELSE
                              fdata(imesh,jmesh,kmesh) = .FALSE.
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
                           fdata(imesh,jmesh,kmesh) = fdata(imesh,jmesh,kmesh)+&
     &                        REAL(ppm_recvbufferd(ibuffer),ppm_kind_single)
#elif  __KIND == __DOUBLE_PRECISION
                           fdata(imesh,jmesh,kmesh) = fdata(imesh,jmesh,kmesh)+&
     &                         ppm_recvbufferd(ibuffer)
#elif  __KIND == __SINGLE_PRECISION_COMPLEX
                           fdata(imesh,jmesh,kmesh) = fdata(imesh,jmesh,kmesh)+&
     &                        CMPLX(ppm_recvbufferd(ibuffer-1),  &
     &                        ppm_recvbufferd(ibuffer),ppm_kind_single)
#elif  __KIND == __DOUBLE_PRECISION_COMPLEX
                           fdata(imesh,jmesh,kmesh) = fdata(imesh,jmesh,kmesh)+&
     &                        CMPLX(ppm_recvbufferd(ibuffer-1),  &
     &                        ppm_recvbufferd(ibuffer),ppm_kind_double)
#elif  __KIND == __INTEGER
                           fdata(imesh,jmesh,kmesh) = fdata(imesh,jmesh,kmesh)+ &
     &                        INT(ppm_recvbufferd(ibuffer))
#elif  __KIND == __LOGICAL
                           IF (ppm_recvbufferd(ibuffer) .GT.     &
     &                        (1.0_ppm_kind_double-ppm_myepsd)) THEN
                              fdata(imesh,jmesh,kmesh) =    &
     &                           (fdata(imesh,jmesh,kmesh) .AND. .TRUE.)
                           ELSE
                              fdata(imesh,jmesh,kmesh) =    &
     &                           (fdata(imesh,jmesh,kmesh) .AND. .FALSE.)
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
                        IF (mask(imesh,jmesh,isub)) THEN
#if    __KIND == __SINGLE_PRECISION_COMPLEX | __KIND == __DOUBLE_PRECISION_COMPLEX
                           ibuffer = ibuffer + 2
#else
                           ibuffer = ibuffer + 1
#endif
                           IF (rtype .EQ. ppm_param_pop_replace) THEN
#if    __KIND == __SINGLE_PRECISION
                              fdata(imesh,jmesh,kmesh) =     &
     &                           REAL(ppm_recvbufferd(ibuffer),ppm_kind_single)
#elif  __KIND == __DOUBLE_PRECISION
                              fdata(imesh,jmesh,kmesh) = ppm_recvbufferd(ibuffer)
#elif  __KIND == __SINGLE_PRECISION_COMPLEX
                              fdata(imesh,jmesh,kmesh) =     &
     &                           CMPLX(ppm_recvbufferd(ibuffer-1),  &
     &                           ppm_recvbufferd(ibuffer),ppm_kind_single)
#elif  __KIND == __DOUBLE_PRECISION_COMPLEX
                              fdata(imesh,jmesh,kmesh) =     &
     &                           CMPLX(ppm_recvbufferd(ibuffer-1),  &
     &                           ppm_recvbufferd(ibuffer),ppm_kind_double)
#elif  __KIND == __INTEGER
                              fdata(imesh,jmesh,kmesh) =     &
     &                           INT(ppm_recvbufferd(ibuffer))
#elif  __KIND == __LOGICAL
                              IF (ppm_recvbufferd(ibuffer) .GT.     &
     &                           (1.0_ppm_kind_double-ppm_myepsd)) THEN
                                 fdata(imesh,jmesh,kmesh) = .TRUE.
                              ELSE
                                 fdata(imesh,jmesh,kmesh) = .FALSE.
                              ENDIF
#endif
                           ELSEIF (rtype .EQ. ppm_param_pop_add) THEN
#if    __KIND == __SINGLE_PRECISION
                              fdata(imesh,jmesh,kmesh) = fdata(imesh,jmesh,kmesh)+&
     &                           REAL(ppm_recvbufferd(ibuffer),ppm_kind_single)
#elif  __KIND == __DOUBLE_PRECISION
                              fdata(imesh,jmesh,kmesh) = fdata(imesh,jmesh,kmesh)+&
     &                            ppm_recvbufferd(ibuffer)
#elif  __KIND == __SINGLE_PRECISION_COMPLEX
                              fdata(imesh,jmesh,kmesh) = fdata(imesh,jmesh,kmesh)+&
     &                           CMPLX(ppm_recvbufferd(ibuffer-1),  &
     &                           ppm_recvbufferd(ibuffer),ppm_kind_single)
#elif  __KIND == __DOUBLE_PRECISION_COMPLEX
                              fdata(imesh,jmesh,kmesh) = fdata(imesh,jmesh,kmesh)+&
     &                           CMPLX(ppm_recvbufferd(ibuffer-1),  &
     &                           ppm_recvbufferd(ibuffer),ppm_kind_double)
#elif  __KIND == __INTEGER
                              fdata(imesh,jmesh,kmesh) = fdata(imesh,jmesh,kmesh)+ &
     &                           INT(ppm_recvbufferd(ibuffer))
#elif  __KIND == __LOGICAL
                              IF (ppm_recvbufferd(ibuffer) .GT.     &
     &                           (1.0_ppm_kind_double-ppm_myepsd)) THEN
                                 fdata(imesh,jmesh,kmesh) =    &
     &                              (fdata(imesh,jmesh,kmesh) .AND. .TRUE.)
                              ELSE
                                 fdata(imesh,jmesh,kmesh) =    &
     &                              (fdata(imesh,jmesh,kmesh) .AND. .FALSE.)
                              ENDIF
#endif
                           ENDIF
                        ELSE
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
      CASE DEFAULT
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
               mofs(1) = this%istart(1,jsub)-1
               mofs(2) = this%istart(2,jsub)-1
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
                   CALL ppm_write(ppm_rank,caller,mesg,info)
                   WRITE(mesg,'(A,2I4)') 'size: ',             &
     &                 ppm_mesh_irecvblksize(1,j),ppm_mesh_irecvblksize(2,j)
                   CALL ppm_write(ppm_rank,caller,mesg,info)
                   WRITE(mesg,'(A,2I4)') 'mesh offset: ',mofs(1),mofs(2)
                   CALL ppm_write(ppm_rank,caller,mesg,info)
                   WRITE(mesg,'(A,2I4)') 'xlo, xhi: ',xlo,xhi
                   CALL ppm_write(ppm_rank,caller,mesg,info)
                   WRITE(mesg,'(A,2I4)') 'ylo, yhi: ',ylo,yhi
                   CALL ppm_write(ppm_rank,caller,mesg,info)
                   WRITE(mesg,'(A,I1)') 'buffer dim: ',edim
                   CALL ppm_write(ppm_rank,caller,mesg,info)
#if   __DIM == __VFIELD
                   WRITE(mesg,'(A,2I4)') 'SIZE(fdata): ',SIZE(fdata,2),  &
     &                 SIZE(fdata,3)
                   CALL ppm_write(ppm_rank,caller,mesg,info)
#elif __DIM == __SFIELD
                   WRITE(mesg,'(A,2I4)') 'SIZE(fdata): ',SIZE(fdata,1),  &
     &                 SIZE(fdata,2)
                   CALL ppm_write(ppm_rank,caller,mesg,info)
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
                              fdata(1,imesh,jmesh,kmesh) =     &
     &                           REAL(ppm_recvbuffers(ibuffer),ppm_kind_double)
#elif  __KIND == __SINGLE_PRECISION
                              fdata(1,imesh,jmesh,kmesh) =     &
     &                           ppm_recvbuffers(ibuffer)
#elif  __KIND == __SINGLE_PRECISION_COMPLEX
                              fdata(1,imesh,jmesh,kmesh) =     &
     &                           CMPLX(ppm_recvbuffers(ibuffer-1),  &
     &                           ppm_recvbuffers(ibuffer),ppm_kind_single)
#elif  __KIND == __DOUBLE_PRECISION_COMPLEX
                              fdata(1,imesh,jmesh,kmesh) =     &
     &                           CMPLX(ppm_recvbuffers(ibuffer-1),  &
     &                           ppm_recvbuffers(ibuffer),ppm_kind_double)
#elif  __KIND == __INTEGER
                              fdata(1,imesh,jmesh,kmesh) =     &
     &                           INT(ppm_recvbuffers(ibuffer))
#elif  __KIND == __LOGICAL
                              IF (ppm_recvbuffers(ibuffer) .GT.     &
     &                           (1.0_ppm_kind_single-ppm_myepss)) THEN
                                 fdata(1,imesh,jmesh,kmesh) = .TRUE.
                              ELSE
                                 fdata(1,imesh,jmesh,kmesh) = .FALSE.
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
                              fdata(1,imesh,jmesh,kmesh) =    &
     &                           fdata(1,imesh,jmesh,kmesh) + &
     &                           REAL(ppm_recvbuffers(ibuffer),ppm_kind_double)
#elif  __KIND == __SINGLE_PRECISION
                              fdata(1,imesh,jmesh,kmesh) =    &
     &                           fdata(1,imesh,jmesh,kmesh) + &
     &                           ppm_recvbuffers(ibuffer)
#elif  __KIND == __SINGLE_PRECISION_COMPLEX
                              fdata(1,imesh,jmesh,kmesh) =    &
     &                           fdata(1,imesh,jmesh,kmesh) + &
     &                           CMPLX(ppm_recvbuffers(ibuffer-1),  &
     &                           ppm_recvbuffers(ibuffer),ppm_kind_single)
#elif  __KIND == __DOUBLE_PRECISION_COMPLEX
                              fdata(1,imesh,jmesh,kmesh) =    &
     &                           fdata(1,imesh,jmesh,kmesh) + &
     &                           CMPLX(ppm_recvbuffers(ibuffer-1),  &
     &                           ppm_recvbuffers(ibuffer),ppm_kind_double)
#elif  __KIND == __INTEGER
                              fdata(1,imesh,jmesh,kmesh) =    &
     &                           fdata(1,imesh,jmesh,kmesh) + &
     &                           INT(ppm_recvbuffers(ibuffer))
#elif  __KIND == __LOGICAL
                              IF (ppm_recvbuffers(ibuffer) .GT.     &
     &                           (1.0_ppm_kind_single-ppm_myepss)) THEN
                                 fdata(1,imesh,jmesh,kmesh) =        &
     &                               (fdata(1,imesh,jmesh,kmesh) .AND. .TRUE.)
                              ELSE
                                 fdata(1,imesh,jmesh,kmesh) =        &
     &                               (fdata(1,imesh,jmesh,kmesh) .AND. .FALSE.)
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
                              fdata(1,imesh,jmesh,kmesh) =     &
     &                           REAL(ppm_recvbuffers(ibuffer),ppm_kind_double)
                              ibuffer = ibuffer + 1
                              fdata(2,imesh,jmesh,kmesh) =     &
     &                           REAL(ppm_recvbuffers(ibuffer),ppm_kind_double)
#elif  __KIND == __SINGLE_PRECISION
                              fdata(1,imesh,jmesh,kmesh) =     &
     &                           ppm_recvbuffers(ibuffer)
                              ibuffer = ibuffer + 1
                              fdata(2,imesh,jmesh,kmesh) =     &
     &                           ppm_recvbuffers(ibuffer)
#elif  __KIND == __SINGLE_PRECISION_COMPLEX
                              fdata(1,imesh,jmesh,kmesh) =     &
     &                           CMPLX(ppm_recvbuffers(ibuffer-1),  &
     &                           ppm_recvbuffers(ibuffer),ppm_kind_single)
                              ibuffer = ibuffer + 2
                              fdata(2,imesh,jmesh,kmesh) =     &
     &                           CMPLX(ppm_recvbuffers(ibuffer-1),  &
     &                           ppm_recvbuffers(ibuffer),ppm_kind_single)
#elif  __KIND == __DOUBLE_PRECISION_COMPLEX
                              fdata(1,imesh,jmesh,kmesh) =     &
     &                           CMPLX(ppm_recvbuffers(ibuffer-1),  &
     &                           ppm_recvbuffers(ibuffer),ppm_kind_double)
                              ibuffer = ibuffer + 2
                              fdata(2,imesh,jmesh,kmesh) =     &
     &                           CMPLX(ppm_recvbuffers(ibuffer-1),  &
     &                           ppm_recvbuffers(ibuffer),ppm_kind_double)
#elif  __KIND == __INTEGER
                              fdata(1,imesh,jmesh,kmesh) =     &
     &                           INT(ppm_recvbuffers(ibuffer))
                              ibuffer = ibuffer + 1
                              fdata(2,imesh,jmesh,kmesh) =     &
     &                           INT(ppm_recvbuffers(ibuffer))
#elif  __KIND == __LOGICAL
                              IF (ppm_recvbuffers(ibuffer) .GT.     &
     &                           (1.0_ppm_kind_single-ppm_myepss)) THEN
                                 fdata(1,imesh,jmesh,kmesh) = .TRUE.
                              ELSE
                                 fdata(1,imesh,jmesh,kmesh) = .FALSE.
                              ENDIF
                              ibuffer = ibuffer + 1
                              IF (ppm_recvbuffers(ibuffer) .GT.     &
     &                           (1.0_ppm_kind_single-ppm_myepss)) THEN
                                 fdata(2,imesh,jmesh,kmesh) = .TRUE.
                              ELSE
                                 fdata(2,imesh,jmesh,kmesh) = .FALSE.
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
                              fdata(1,imesh,jmesh,kmesh) =    &
     &                           fdata(1,imesh,jmesh,kmesh) + &
     &                           REAL(ppm_recvbuffers(ibuffer),ppm_kind_double)
                              ibuffer = ibuffer + 1
                              fdata(2,imesh,jmesh,kmesh) =    &
     &                           fdata(2,imesh,jmesh,kmesh) + &
     &                           REAL(ppm_recvbuffers(ibuffer),ppm_kind_double)
#elif  __KIND == __SINGLE_PRECISION
                              fdata(1,imesh,jmesh,kmesh) =    &
     &                           fdata(1,imesh,jmesh,kmesh) + &
     &                           ppm_recvbuffers(ibuffer)
                              ibuffer = ibuffer + 1
                              fdata(2,imesh,jmesh,kmesh) =    &
     &                           fdata(2,imesh,jmesh,kmesh) + &
     &                           ppm_recvbuffers(ibuffer)
#elif  __KIND == __SINGLE_PRECISION_COMPLEX
                              fdata(1,imesh,jmesh,kmesh) =    &
     &                           fdata(1,imesh,jmesh,kmesh) + &
     &                           CMPLX(ppm_recvbuffers(ibuffer-1),  &
     &                           ppm_recvbuffers(ibuffer),ppm_kind_single)
                              ibuffer = ibuffer + 2
                              fdata(2,imesh,jmesh,kmesh) =    &
     &                           fdata(2,imesh,jmesh,kmesh) + &
     &                           CMPLX(ppm_recvbuffers(ibuffer-1),  &
     &                           ppm_recvbuffers(ibuffer),ppm_kind_single)
#elif  __KIND == __DOUBLE_PRECISION_COMPLEX
                              fdata(1,imesh,jmesh,kmesh) =    &
     &                           fdata(1,imesh,jmesh,kmesh) + &
     &                           CMPLX(ppm_recvbuffers(ibuffer-1),  &
     &                           ppm_recvbuffers(ibuffer),ppm_kind_double)
                              ibuffer = ibuffer + 2
                              fdata(2,imesh,jmesh,kmesh) =    &
     &                           fdata(2,imesh,jmesh,kmesh) + &
     &                           CMPLX(ppm_recvbuffers(ibuffer-1),  &
     &                           ppm_recvbuffers(ibuffer),ppm_kind_double)
#elif  __KIND == __INTEGER
                              fdata(1,imesh,jmesh,kmesh) =    &
     &                           fdata(1,imesh,jmesh,kmesh) + &
     &                           INT(ppm_recvbuffers(ibuffer))
                              ibuffer = ibuffer + 1
                              fdata(2,imesh,jmesh,kmesh) =    &
     &                           fdata(2,imesh,jmesh,kmesh) + &
     &                           INT(ppm_recvbuffers(ibuffer))
#elif  __KIND == __LOGICAL
                              IF (ppm_recvbuffers(ibuffer) .GT.     &
     &                           (1.0_ppm_kind_single-ppm_myepss)) THEN
                                 fdata(1,imesh,jmesh,kmesh) =        &
     &                               (fdata(1,imesh,jmesh,kmesh) .AND. .TRUE.)
                              ELSE
                                 fdata(1,imesh,jmesh,kmesh) =        &
     &                               (fdata(1,imesh,jmesh,kmesh) .AND. .FALSE.)
                              ENDIF
                              ibuffer = ibuffer + 1
                              IF (ppm_recvbuffers(ibuffer) .GT.     &
     &                           (1.0_ppm_kind_single-ppm_myepss)) THEN
                                 fdata(2,imesh,jmesh,kmesh) =        &
     &                               (fdata(2,imesh,jmesh,kmesh) .AND. .TRUE.)
                              ELSE
                                 fdata(2,imesh,jmesh,kmesh) =        &
     &                               (fdata(2,imesh,jmesh,kmesh) .AND. .FALSE.)
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
                              fdata(1,imesh,jmesh,kmesh) =     &
     &                           REAL(ppm_recvbuffers(ibuffer),ppm_kind_double)
                              ibuffer = ibuffer + 1
                              fdata(2,imesh,jmesh,kmesh) =     &
     &                           REAL(ppm_recvbuffers(ibuffer),ppm_kind_double)
                              ibuffer = ibuffer + 1
                              fdata(3,imesh,jmesh,kmesh) =     &
     &                           REAL(ppm_recvbuffers(ibuffer),ppm_kind_double)
#elif  __KIND == __SINGLE_PRECISION
                              fdata(1,imesh,jmesh,kmesh) =     &
     &                           ppm_recvbuffers(ibuffer)
                              ibuffer = ibuffer + 1
                              fdata(2,imesh,jmesh,kmesh) =     &
     &                           ppm_recvbuffers(ibuffer)
                              ibuffer = ibuffer + 1
                              fdata(3,imesh,jmesh,kmesh) =     &
     &                           ppm_recvbuffers(ibuffer)
#elif  __KIND == __SINGLE_PRECISION_COMPLEX
                              fdata(1,imesh,jmesh,kmesh) =     &
     &                           CMPLX(ppm_recvbuffers(ibuffer-1),  &
     &                           ppm_recvbuffers(ibuffer),ppm_kind_single)
                              ibuffer = ibuffer + 2
                              fdata(2,imesh,jmesh,kmesh) =     &
     &                           CMPLX(ppm_recvbuffers(ibuffer-1),  &
     &                           ppm_recvbuffers(ibuffer),ppm_kind_single)
                              ibuffer = ibuffer + 2
                              fdata(3,imesh,jmesh,kmesh) =     &
     &                           CMPLX(ppm_recvbuffers(ibuffer-1),  &
     &                           ppm_recvbuffers(ibuffer),ppm_kind_single)
#elif  __KIND == __DOUBLE_PRECISION_COMPLEX
                              fdata(1,imesh,jmesh,kmesh) =     &
     &                           CMPLX(ppm_recvbuffers(ibuffer-1),  &
     &                           ppm_recvbuffers(ibuffer),ppm_kind_double)
                              ibuffer = ibuffer + 2
                              fdata(2,imesh,jmesh,kmesh) =     &
     &                           CMPLX(ppm_recvbuffers(ibuffer-1),  &
     &                           ppm_recvbuffers(ibuffer),ppm_kind_double)
                              ibuffer = ibuffer + 2
                              fdata(3,imesh,jmesh,kmesh) =     &
     &                           CMPLX(ppm_recvbuffers(ibuffer-1),  &
     &                           ppm_recvbuffers(ibuffer),ppm_kind_double)
#elif  __KIND == __INTEGER
                              fdata(1,imesh,jmesh,kmesh) =     &
     &                           INT(ppm_recvbuffers(ibuffer))
                              ibuffer = ibuffer + 1
                              fdata(2,imesh,jmesh,kmesh) =     &
     &                           INT(ppm_recvbuffers(ibuffer))
                              ibuffer = ibuffer + 1
                              fdata(3,imesh,jmesh,kmesh) =     &
     &                           INT(ppm_recvbuffers(ibuffer))
#elif  __KIND == __LOGICAL
                              IF (ppm_recvbuffers(ibuffer) .GT.     &
     &                           (1.0_ppm_kind_single-ppm_myepss)) THEN
                                 fdata(1,imesh,jmesh,kmesh) = .TRUE.
                              ELSE
                                 fdata(1,imesh,jmesh,kmesh) = .FALSE.
                              ENDIF
                              ibuffer = ibuffer + 1
                              IF (ppm_recvbuffers(ibuffer) .GT.     &
     &                           (1.0_ppm_kind_single-ppm_myepss)) THEN
                                 fdata(2,imesh,jmesh,kmesh) = .TRUE.
                              ELSE
                                 fdata(2,imesh,jmesh,kmesh) = .FALSE.
                              ENDIF
                              ibuffer = ibuffer + 1
                              IF (ppm_recvbuffers(ibuffer) .GT.     &
     &                           (1.0_ppm_kind_single-ppm_myepss)) THEN
                                 fdata(3,imesh,jmesh,kmesh) = .TRUE.
                              ELSE
                                 fdata(3,imesh,jmesh,kmesh) = .FALSE.
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
                              fdata(1,imesh,jmesh,kmesh) =    &
     &                           fdata(1,imesh,jmesh,kmesh) + &
     &                           REAL(ppm_recvbuffers(ibuffer),ppm_kind_double)
                              ibuffer = ibuffer + 1
                              fdata(2,imesh,jmesh,kmesh) =    &
     &                           fdata(2,imesh,jmesh,kmesh) + &
     &                           REAL(ppm_recvbuffers(ibuffer),ppm_kind_double)
                              ibuffer = ibuffer + 1
                              fdata(3,imesh,jmesh,kmesh) =    &
     &                           fdata(3,imesh,jmesh,kmesh) + &
     &                           REAL(ppm_recvbuffers(ibuffer),ppm_kind_double)
#elif  __KIND == __SINGLE_PRECISION
                              fdata(1,imesh,jmesh,kmesh) =    &
     &                           fdata(1,imesh,jmesh,kmesh) + &
     &                           ppm_recvbuffers(ibuffer)
                              ibuffer = ibuffer + 1
                              fdata(2,imesh,jmesh,kmesh) =    &
     &                           fdata(2,imesh,jmesh,kmesh) + &
     &                           ppm_recvbuffers(ibuffer)
                              ibuffer = ibuffer + 1
                              fdata(3,imesh,jmesh,kmesh) =    &
     &                           fdata(3,imesh,jmesh,kmesh) + &
     &                           ppm_recvbuffers(ibuffer)
#elif  __KIND == __SINGLE_PRECISION_COMPLEX
                              fdata(1,imesh,jmesh,kmesh) =    &
     &                           fdata(1,imesh,jmesh,kmesh) + &
     &                           CMPLX(ppm_recvbuffers(ibuffer-1),  &
     &                           ppm_recvbuffers(ibuffer),ppm_kind_single)
                              ibuffer = ibuffer + 2
                              fdata(2,imesh,jmesh,kmesh) =    &
     &                           fdata(2,imesh,jmesh,kmesh) + &
     &                           CMPLX(ppm_recvbuffers(ibuffer-1),  &
     &                           ppm_recvbuffers(ibuffer),ppm_kind_single)
                              ibuffer = ibuffer + 2
                              fdata(3,imesh,jmesh,kmesh) =    &
     &                           fdata(3,imesh,jmesh,kmesh) + &
     &                           CMPLX(ppm_recvbuffers(ibuffer-1),  &
     &                           ppm_recvbuffers(ibuffer),ppm_kind_single)
#elif  __KIND == __DOUBLE_PRECISION_COMPLEX
                              fdata(1,imesh,jmesh,kmesh) =    &
     &                           fdata(1,imesh,jmesh,kmesh) + &
     &                           CMPLX(ppm_recvbuffers(ibuffer-1),  &
     &                           ppm_recvbuffers(ibuffer),ppm_kind_double)
                              ibuffer = ibuffer + 2
                              fdata(2,imesh,jmesh,kmesh) =    &
     &                           fdata(2,imesh,jmesh,kmesh) + &
     &                           CMPLX(ppm_recvbuffers(ibuffer-1),  &
     &                           ppm_recvbuffers(ibuffer),ppm_kind_double)
                              ibuffer = ibuffer + 2
                              fdata(3,imesh,jmesh,kmesh) =    &
     &                           fdata(3,imesh,jmesh,kmesh) + &
     &                           CMPLX(ppm_recvbuffers(ibuffer-1),  &
     &                           ppm_recvbuffers(ibuffer),ppm_kind_double)
#elif  __KIND == __INTEGER
                              fdata(1,imesh,jmesh,kmesh) =    &
     &                           fdata(1,imesh,jmesh,kmesh) + &
     &                           INT(ppm_recvbuffers(ibuffer))
                              ibuffer = ibuffer + 1
                              fdata(2,imesh,jmesh,kmesh) =    &
     &                           fdata(2,imesh,jmesh,kmesh) + &
     &                           INT(ppm_recvbuffers(ibuffer))
                              ibuffer = ibuffer + 1
                              fdata(3,imesh,jmesh,kmesh) =    &
     &                           fdata(3,imesh,jmesh,kmesh) + &
     &                           INT(ppm_recvbuffers(ibuffer))
#elif  __KIND == __LOGICAL
                              IF (ppm_recvbuffers(ibuffer) .GT.     &
     &                           (1.0_ppm_kind_single-ppm_myepss)) THEN
                                 fdata(1,imesh,jmesh,kmesh) =        &
     &                               (fdata(1,imesh,jmesh,kmesh) .AND. .TRUE.)
                              ELSE
                                 fdata(1,imesh,jmesh,kmesh) =        &
     &                               (fdata(1,imesh,jmesh,kmesh) .AND. .FALSE.)
                              ENDIF
                              ibuffer = ibuffer + 1
                              IF (ppm_recvbuffers(ibuffer) .GT.     &
     &                           (1.0_ppm_kind_single-ppm_myepss)) THEN
                                 fdata(2,imesh,jmesh,kmesh) =        &
     &                               (fdata(2,imesh,jmesh,kmesh) .AND. .TRUE.)
                              ELSE
                                 fdata(2,imesh,jmesh,kmesh) =        &
     &                               (fdata(2,imesh,jmesh,kmesh) .AND. .FALSE.)
                              ENDIF
                              ibuffer = ibuffer + 1
                              IF (ppm_recvbuffers(ibuffer) .GT.     &
     &                           (1.0_ppm_kind_single-ppm_myepss)) THEN
                                 fdata(3,imesh,jmesh,kmesh) =        &
     &                               (fdata(3,imesh,jmesh,kmesh) .AND. .TRUE.)
                              ELSE
                                 fdata(3,imesh,jmesh,kmesh) =        &
     &                               (fdata(3,imesh,jmesh,kmesh) .AND. .FALSE.)
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
                              fdata(1,imesh,jmesh,kmesh) =     &
     &                           REAL(ppm_recvbuffers(ibuffer),ppm_kind_double)
                              ibuffer = ibuffer + 1
                              fdata(2,imesh,jmesh,kmesh) =     &
     &                           REAL(ppm_recvbuffers(ibuffer),ppm_kind_double)
                              ibuffer = ibuffer + 1
                              fdata(3,imesh,jmesh,kmesh) =     &
     &                           REAL(ppm_recvbuffers(ibuffer),ppm_kind_double)
                              ibuffer = ibuffer + 1
                              fdata(4,imesh,jmesh,kmesh) =     &
     &                           REAL(ppm_recvbuffers(ibuffer),ppm_kind_double)
#elif  __KIND == __SINGLE_PRECISION
                              fdata(1,imesh,jmesh,kmesh) =     &
     &                           ppm_recvbuffers(ibuffer)
                              ibuffer = ibuffer + 1
                              fdata(2,imesh,jmesh,kmesh) =     &
     &                           ppm_recvbuffers(ibuffer)
                              ibuffer = ibuffer + 1
                              fdata(3,imesh,jmesh,kmesh) =     &
     &                           ppm_recvbuffers(ibuffer)
                              ibuffer = ibuffer + 1
                              fdata(4,imesh,jmesh,kmesh) =     &
     &                           ppm_recvbuffers(ibuffer)
#elif  __KIND == __SINGLE_PRECISION_COMPLEX
                              fdata(1,imesh,jmesh,kmesh) =     &
     &                           CMPLX(ppm_recvbuffers(ibuffer-1),  &
     &                           ppm_recvbuffers(ibuffer),ppm_kind_single)
                              ibuffer = ibuffer + 2
                              fdata(2,imesh,jmesh,kmesh) =     &
     &                           CMPLX(ppm_recvbuffers(ibuffer-1),  &
     &                           ppm_recvbuffers(ibuffer),ppm_kind_single)
                              ibuffer = ibuffer + 2
                              fdata(3,imesh,jmesh,kmesh) =     &
     &                           CMPLX(ppm_recvbuffers(ibuffer-1),  &
     &                           ppm_recvbuffers(ibuffer),ppm_kind_single)
                              ibuffer = ibuffer + 2
                              fdata(4,imesh,jmesh,kmesh) =     &
     &                           CMPLX(ppm_recvbuffers(ibuffer-1),  &
     &                           ppm_recvbuffers(ibuffer),ppm_kind_single)
#elif  __KIND == __DOUBLE_PRECISION_COMPLEX
                              fdata(1,imesh,jmesh,kmesh) =     &
     &                           CMPLX(ppm_recvbuffers(ibuffer-1),  &
     &                           ppm_recvbuffers(ibuffer),ppm_kind_double)
                              ibuffer = ibuffer + 2
                              fdata(2,imesh,jmesh,kmesh) =     &
     &                           CMPLX(ppm_recvbuffers(ibuffer-1),  &
     &                           ppm_recvbuffers(ibuffer),ppm_kind_double)
                              ibuffer = ibuffer + 2
                              fdata(3,imesh,jmesh,kmesh) =     &
     &                           CMPLX(ppm_recvbuffers(ibuffer-1),  &
     &                           ppm_recvbuffers(ibuffer),ppm_kind_double)
                              ibuffer = ibuffer + 2
                              fdata(4,imesh,jmesh,kmesh) =     &
     &                           CMPLX(ppm_recvbuffers(ibuffer-1),  &
     &                           ppm_recvbuffers(ibuffer),ppm_kind_double)
#elif  __KIND == __INTEGER
                              fdata(1,imesh,jmesh,kmesh) =     &
     &                           INT(ppm_recvbuffers(ibuffer))
                              ibuffer = ibuffer + 1
                              fdata(2,imesh,jmesh,kmesh) =     &
     &                           INT(ppm_recvbuffers(ibuffer))
                              ibuffer = ibuffer + 1
                              fdata(3,imesh,jmesh,kmesh) =     &
     &                           INT(ppm_recvbuffers(ibuffer))
                              ibuffer = ibuffer + 1
                              fdata(4,imesh,jmesh,kmesh) =     &
     &                           INT(ppm_recvbuffers(ibuffer))
#elif  __KIND == __LOGICAL
                              IF (ppm_recvbuffers(ibuffer) .GT.     &
     &                           (1.0_ppm_kind_single-ppm_myepss)) THEN
                                 fdata(1,imesh,jmesh,kmesh) = .TRUE.
                              ELSE
                                 fdata(1,imesh,jmesh,kmesh) = .FALSE.
                              ENDIF
                              ibuffer = ibuffer + 1
                              IF (ppm_recvbuffers(ibuffer) .GT.     &
     &                           (1.0_ppm_kind_single-ppm_myepss)) THEN
                                 fdata(2,imesh,jmesh,kmesh) = .TRUE.
                              ELSE
                                 fdata(2,imesh,jmesh,kmesh) = .FALSE.
                              ENDIF
                              ibuffer = ibuffer + 1
                              IF (ppm_recvbuffers(ibuffer) .GT.     &
     &                           (1.0_ppm_kind_single-ppm_myepss)) THEN
                                 fdata(3,imesh,jmesh,kmesh) = .TRUE.
                              ELSE
                                 fdata(3,imesh,jmesh,kmesh) = .FALSE.
                              ENDIF
                              ibuffer = ibuffer + 1
                              IF (ppm_recvbuffers(ibuffer) .GT.     &
     &                           (1.0_ppm_kind_single-ppm_myepss)) THEN
                                 fdata(4,imesh,jmesh,kmesh) = .TRUE.
                              ELSE
                                 fdata(4,imesh,jmesh,kmesh) = .FALSE.
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
                              fdata(1,imesh,jmesh,kmesh) =    &
     &                           fdata(1,imesh,jmesh,kmesh) + &
     &                           REAL(ppm_recvbuffers(ibuffer),ppm_kind_double)
                              ibuffer = ibuffer + 1
                              fdata(2,imesh,jmesh,kmesh) =    &
     &                           fdata(2,imesh,jmesh,kmesh) + &
     &                           REAL(ppm_recvbuffers(ibuffer),ppm_kind_double)
                              ibuffer = ibuffer + 1
                              fdata(3,imesh,jmesh,kmesh) =    &
     &                           fdata(3,imesh,jmesh,kmesh) + &
     &                           REAL(ppm_recvbuffers(ibuffer),ppm_kind_double)
                              ibuffer = ibuffer + 1
                              fdata(4,imesh,jmesh,kmesh) =    &
     &                           fdata(4,imesh,jmesh,kmesh) + &
     &                           REAL(ppm_recvbuffers(ibuffer),ppm_kind_double)
#elif  __KIND == __SINGLE_PRECISION
                              fdata(1,imesh,jmesh,kmesh) =    &
     &                           fdata(1,imesh,jmesh,kmesh) + &
     &                           ppm_recvbuffers(ibuffer)
                              ibuffer = ibuffer + 1
                              fdata(2,imesh,jmesh,kmesh) =    &
     &                           fdata(2,imesh,jmesh,kmesh) + &
     &                           ppm_recvbuffers(ibuffer)
                              ibuffer = ibuffer + 1
                              fdata(3,imesh,jmesh,kmesh) =    &
     &                           fdata(3,imesh,jmesh,kmesh) + &
     &                           ppm_recvbuffers(ibuffer)
                              ibuffer = ibuffer + 1
                              fdata(4,imesh,jmesh,kmesh) =    &
     &                           fdata(4,imesh,jmesh,kmesh) + &
     &                           ppm_recvbuffers(ibuffer)
#elif  __KIND == __SINGLE_PRECISION_COMPLEX
                              fdata(1,imesh,jmesh,kmesh) =    &
     &                           fdata(1,imesh,jmesh,kmesh) + &
     &                           CMPLX(ppm_recvbuffers(ibuffer-1),  &
     &                           ppm_recvbuffers(ibuffer),ppm_kind_single)
                              ibuffer = ibuffer + 2
                              fdata(2,imesh,jmesh,kmesh) =    &
     &                           fdata(2,imesh,jmesh,kmesh) + &
     &                           CMPLX(ppm_recvbuffers(ibuffer-1),  &
     &                           ppm_recvbuffers(ibuffer),ppm_kind_single)
                              ibuffer = ibuffer + 2
                              fdata(3,imesh,jmesh,kmesh) =    &
     &                           fdata(3,imesh,jmesh,kmesh) + &
     &                           CMPLX(ppm_recvbuffers(ibuffer-1),  &
     &                           ppm_recvbuffers(ibuffer),ppm_kind_single)
                              ibuffer = ibuffer + 2
                              fdata(4,imesh,jmesh,kmesh) =    &
     &                           fdata(4,imesh,jmesh,kmesh) + &
     &                           CMPLX(ppm_recvbuffers(ibuffer-1),  &
     &                           ppm_recvbuffers(ibuffer),ppm_kind_single)
#elif  __KIND == __DOUBLE_PRECISION_COMPLEX
                              fdata(1,imesh,jmesh,kmesh) =    &
     &                           fdata(1,imesh,jmesh,kmesh) + &
     &                           CMPLX(ppm_recvbuffers(ibuffer-1),  &
     &                           ppm_recvbuffers(ibuffer),ppm_kind_double)
                              ibuffer = ibuffer + 2
                              fdata(2,imesh,jmesh,kmesh) =    &
     &                           fdata(2,imesh,jmesh,kmesh) + &
     &                           CMPLX(ppm_recvbuffers(ibuffer-1),  &
     &                           ppm_recvbuffers(ibuffer),ppm_kind_double)
                              ibuffer = ibuffer + 2
                              fdata(3,imesh,jmesh,kmesh) =    &
     &                           fdata(3,imesh,jmesh,kmesh) + &
     &                           CMPLX(ppm_recvbuffers(ibuffer-1),  &
     &                           ppm_recvbuffers(ibuffer),ppm_kind_double)
                              ibuffer = ibuffer + 2
                              fdata(4,imesh,jmesh,kmesh) =    &
     &                           fdata(4,imesh,jmesh,kmesh) + &
     &                           CMPLX(ppm_recvbuffers(ibuffer-1),  &
     &                           ppm_recvbuffers(ibuffer),ppm_kind_double)
#elif  __KIND == __INTEGER
                              fdata(1,imesh,jmesh,kmesh) =    &
     &                           fdata(1,imesh,jmesh,kmesh) + &
     &                           INT(ppm_recvbuffers(ibuffer))
                              ibuffer = ibuffer + 1
                              fdata(2,imesh,jmesh,kmesh) =    &
     &                           fdata(2,imesh,jmesh,kmesh) + &
     &                           INT(ppm_recvbuffers(ibuffer))
                              ibuffer = ibuffer + 1
                              fdata(3,imesh,jmesh,kmesh) =    &
     &                           fdata(3,imesh,jmesh,kmesh) + &
     &                           INT(ppm_recvbuffers(ibuffer))
                              ibuffer = ibuffer + 1
                              fdata(4,imesh,jmesh,kmesh) =    &
     &                           fdata(4,imesh,jmesh,kmesh) + &
     &                           INT(ppm_recvbuffers(ibuffer))
#elif  __KIND == __LOGICAL
                              IF (ppm_recvbuffers(ibuffer) .GT.     &
     &                           (1.0_ppm_kind_single-ppm_myepss)) THEN
                                 fdata(1,imesh,jmesh,kmesh) =        &
     &                               (fdata(1,imesh,jmesh,kmesh) .AND. .TRUE.)
                              ELSE
                                 fdata(1,imesh,jmesh,kmesh) =        &
     &                               (fdata(1,imesh,jmesh,kmesh) .AND. .FALSE.)
                              ENDIF
                              ibuffer = ibuffer + 1
                              IF (ppm_recvbuffers(ibuffer) .GT.     &
     &                           (1.0_ppm_kind_single-ppm_myepss)) THEN
                                 fdata(2,imesh,jmesh,kmesh) =        &
     &                               (fdata(2,imesh,jmesh,kmesh) .AND. .TRUE.)
                              ELSE
                                 fdata(2,imesh,jmesh,kmesh) =        &
     &                               (fdata(2,imesh,jmesh,kmesh) .AND. .FALSE.)
                              ENDIF
                              ibuffer = ibuffer + 1
                              IF (ppm_recvbuffers(ibuffer) .GT.     &
     &                           (1.0_ppm_kind_single-ppm_myepss)) THEN
                                 fdata(3,imesh,jmesh,kmesh) =        &
     &                               (fdata(3,imesh,jmesh,kmesh) .AND. .TRUE.)
                              ELSE
                                 fdata(3,imesh,jmesh,kmesh) =        &
     &                               (fdata(3,imesh,jmesh,kmesh) .AND. .FALSE.)
                              ENDIF
                              ibuffer = ibuffer + 1
                              IF (ppm_recvbuffers(ibuffer) .GT.     &
     &                           (1.0_ppm_kind_single-ppm_myepss)) THEN
                                 fdata(4,imesh,jmesh,kmesh) =        &
     &                               (fdata(4,imesh,jmesh,kmesh) .AND. .TRUE.)
                              ELSE
                                 fdata(4,imesh,jmesh,kmesh) =        &
     &                               (fdata(4,imesh,jmesh,kmesh) .AND. .FALSE.)
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
                              fdata(1,imesh,jmesh,kmesh) =     &
     &                           REAL(ppm_recvbuffers(ibuffer),ppm_kind_double)
                              ibuffer = ibuffer + 1
                              fdata(2,imesh,jmesh,kmesh) =     &
     &                           REAL(ppm_recvbuffers(ibuffer),ppm_kind_double)
                              ibuffer = ibuffer + 1
                              fdata(3,imesh,jmesh,kmesh) =     &
     &                           REAL(ppm_recvbuffers(ibuffer),ppm_kind_double)
                              ibuffer = ibuffer + 1
                              fdata(4,imesh,jmesh,kmesh) =     &
     &                           REAL(ppm_recvbuffers(ibuffer),ppm_kind_double)
                              ibuffer = ibuffer + 1
                              fdata(5,imesh,jmesh,kmesh) =     &
     &                           REAL(ppm_recvbuffers(ibuffer),ppm_kind_double)
#elif  __KIND == __SINGLE_PRECISION
                              fdata(1,imesh,jmesh,kmesh) =     &
     &                           ppm_recvbuffers(ibuffer)
                              ibuffer = ibuffer + 1
                              fdata(2,imesh,jmesh,kmesh) =     &
     &                           ppm_recvbuffers(ibuffer)
                              ibuffer = ibuffer + 1
                              fdata(3,imesh,jmesh,kmesh) =     &
     &                           ppm_recvbuffers(ibuffer)
                              ibuffer = ibuffer + 1
                              fdata(4,imesh,jmesh,kmesh) =     &
     &                           ppm_recvbuffers(ibuffer)
                              ibuffer = ibuffer + 1
                              fdata(5,imesh,jmesh,kmesh) =     &
     &                           ppm_recvbuffers(ibuffer)
#elif  __KIND == __SINGLE_PRECISION_COMPLEX
                              fdata(1,imesh,jmesh,kmesh) =     &
     &                           CMPLX(ppm_recvbuffers(ibuffer-1),  &
     &                           ppm_recvbuffers(ibuffer),ppm_kind_single)
                              ibuffer = ibuffer + 2
                              fdata(2,imesh,jmesh,kmesh) =     &
     &                           CMPLX(ppm_recvbuffers(ibuffer-1),  &
     &                           ppm_recvbuffers(ibuffer),ppm_kind_single)
                              ibuffer = ibuffer + 2
                              fdata(3,imesh,jmesh,kmesh) =     &
     &                           CMPLX(ppm_recvbuffers(ibuffer-1),  &
     &                           ppm_recvbuffers(ibuffer),ppm_kind_single)
                              ibuffer = ibuffer + 2
                              fdata(4,imesh,jmesh,kmesh) =     &
     &                           CMPLX(ppm_recvbuffers(ibuffer-1),  &
     &                           ppm_recvbuffers(ibuffer),ppm_kind_single)
                              ibuffer = ibuffer + 2
                              fdata(5,imesh,jmesh,kmesh) =     &
     &                           CMPLX(ppm_recvbuffers(ibuffer-1),  &
     &                           ppm_recvbuffers(ibuffer),ppm_kind_single)
#elif  __KIND == __DOUBLE_PRECISION_COMPLEX
                              fdata(1,imesh,jmesh,kmesh) =     &
     &                           CMPLX(ppm_recvbuffers(ibuffer-1),  &
     &                           ppm_recvbuffers(ibuffer),ppm_kind_double)
                              ibuffer = ibuffer + 2
                              fdata(2,imesh,jmesh,kmesh) =     &
     &                           CMPLX(ppm_recvbuffers(ibuffer-1),  &
     &                           ppm_recvbuffers(ibuffer),ppm_kind_double)
                              ibuffer = ibuffer + 2
                              fdata(3,imesh,jmesh,kmesh) =     &
     &                           CMPLX(ppm_recvbuffers(ibuffer-1),  &
     &                           ppm_recvbuffers(ibuffer),ppm_kind_double)
                              ibuffer = ibuffer + 2
                              fdata(4,imesh,jmesh,kmesh) =     &
     &                           CMPLX(ppm_recvbuffers(ibuffer-1),  &
     &                           ppm_recvbuffers(ibuffer),ppm_kind_double)
                              ibuffer = ibuffer + 2
                              fdata(5,imesh,jmesh,kmesh) =     &
     &                           CMPLX(ppm_recvbuffers(ibuffer-1),  &
     &                           ppm_recvbuffers(ibuffer),ppm_kind_double)
#elif  __KIND == __INTEGER
                              fdata(1,imesh,jmesh,kmesh) =     &
     &                           INT(ppm_recvbuffers(ibuffer))
                              ibuffer = ibuffer + 1
                              fdata(2,imesh,jmesh,kmesh) =     &
     &                           INT(ppm_recvbuffers(ibuffer))
                              ibuffer = ibuffer + 1
                              fdata(3,imesh,jmesh,kmesh) =     &
     &                           INT(ppm_recvbuffers(ibuffer))
                              ibuffer = ibuffer + 1
                              fdata(4,imesh,jmesh,kmesh) =     &
     &                           INT(ppm_recvbuffers(ibuffer))
                              ibuffer = ibuffer + 1
                              fdata(5,imesh,jmesh,kmesh) =     &
     &                           INT(ppm_recvbuffers(ibuffer))
#elif  __KIND == __LOGICAL
                              IF (ppm_recvbuffers(ibuffer) .GT.     &
     &                           (1.0_ppm_kind_single-ppm_myepss)) THEN
                                 fdata(1,imesh,jmesh,kmesh) = .TRUE.
                              ELSE
                                 fdata(1,imesh,jmesh,kmesh) = .FALSE.
                              ENDIF
                              ibuffer = ibuffer + 1
                              IF (ppm_recvbuffers(ibuffer) .GT.     &
     &                           (1.0_ppm_kind_single-ppm_myepss)) THEN
                                 fdata(2,imesh,jmesh,kmesh) = .TRUE.
                              ELSE
                                 fdata(2,imesh,jmesh,kmesh) = .FALSE.
                              ENDIF
                              ibuffer = ibuffer + 1
                              IF (ppm_recvbuffers(ibuffer) .GT.     &
     &                           (1.0_ppm_kind_single-ppm_myepss)) THEN
                                 fdata(3,imesh,jmesh,kmesh) = .TRUE.
                              ELSE
                                 fdata(3,imesh,jmesh,kmesh) = .FALSE.
                              ENDIF
                              ibuffer = ibuffer + 1
                              IF (ppm_recvbuffers(ibuffer) .GT.     &
     &                           (1.0_ppm_kind_single-ppm_myepss)) THEN
                                 fdata(4,imesh,jmesh,kmesh) = .TRUE.
                              ELSE
                                 fdata(4,imesh,jmesh,kmesh) = .FALSE.
                              ENDIF
                              ibuffer = ibuffer + 1
                              IF (ppm_recvbuffers(ibuffer) .GT.     &
     &                           (1.0_ppm_kind_single-ppm_myepss)) THEN
                                 fdata(5,imesh,jmesh,kmesh) = .TRUE.
                              ELSE
                                 fdata(5,imesh,jmesh,kmesh) = .FALSE.
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
                              fdata(1,imesh,jmesh,kmesh) =    &
     &                           fdata(1,imesh,jmesh,kmesh) + &
     &                           REAL(ppm_recvbuffers(ibuffer),ppm_kind_double)
                              ibuffer = ibuffer + 1
                              fdata(2,imesh,jmesh,kmesh) =    &
     &                           fdata(2,imesh,jmesh,kmesh) + &
     &                           REAL(ppm_recvbuffers(ibuffer),ppm_kind_double)
                              ibuffer = ibuffer + 1
                              fdata(3,imesh,jmesh,kmesh) =    &
     &                           fdata(3,imesh,jmesh,kmesh) + &
     &                           REAL(ppm_recvbuffers(ibuffer),ppm_kind_double)
                              ibuffer = ibuffer + 1
                              fdata(4,imesh,jmesh,kmesh) =    &
     &                           fdata(4,imesh,jmesh,kmesh) + &
     &                           REAL(ppm_recvbuffers(ibuffer),ppm_kind_double)
                              ibuffer = ibuffer + 1
                              fdata(5,imesh,jmesh,kmesh) =    &
     &                           fdata(5,imesh,jmesh,kmesh) + &
     &                           REAL(ppm_recvbuffers(ibuffer),ppm_kind_double)
#elif  __KIND == __SINGLE_PRECISION
                              fdata(1,imesh,jmesh,kmesh) =    &
     &                           fdata(1,imesh,jmesh,kmesh) + &
     &                           ppm_recvbuffers(ibuffer)
                              ibuffer = ibuffer + 1
                              fdata(2,imesh,jmesh,kmesh) =    &
     &                           fdata(2,imesh,jmesh,kmesh) + &
     &                           ppm_recvbuffers(ibuffer)
                              ibuffer = ibuffer + 1
                              fdata(3,imesh,jmesh,kmesh) =    &
     &                           fdata(3,imesh,jmesh,kmesh) + &
     &                           ppm_recvbuffers(ibuffer)
                              ibuffer = ibuffer + 1
                              fdata(4,imesh,jmesh,kmesh) =    &
     &                           fdata(4,imesh,jmesh,kmesh) + &
     &                           ppm_recvbuffers(ibuffer)
                              ibuffer = ibuffer + 1
                              fdata(5,imesh,jmesh,kmesh) =    &
     &                           fdata(5,imesh,jmesh,kmesh) + &
     &                           ppm_recvbuffers(ibuffer)
#elif  __KIND == __SINGLE_PRECISION_COMPLEX
                              fdata(1,imesh,jmesh,kmesh) =    &
     &                           fdata(1,imesh,jmesh,kmesh) + &
     &                           CMPLX(ppm_recvbuffers(ibuffer-1),  &
     &                           ppm_recvbuffers(ibuffer),ppm_kind_single)
                              ibuffer = ibuffer + 2
                              fdata(2,imesh,jmesh,kmesh) =    &
     &                           fdata(2,imesh,jmesh,kmesh) + &
     &                           CMPLX(ppm_recvbuffers(ibuffer-1),  &
     &                           ppm_recvbuffers(ibuffer),ppm_kind_single)
                              ibuffer = ibuffer + 2
                              fdata(3,imesh,jmesh,kmesh) =    &
     &                           fdata(3,imesh,jmesh,kmesh) + &
     &                           CMPLX(ppm_recvbuffers(ibuffer-1),  &
     &                           ppm_recvbuffers(ibuffer),ppm_kind_single)
                              ibuffer = ibuffer + 2
                              fdata(4,imesh,jmesh,kmesh) =    &
     &                           fdata(4,imesh,jmesh,kmesh) + &
     &                           CMPLX(ppm_recvbuffers(ibuffer-1),  &
     &                           ppm_recvbuffers(ibuffer),ppm_kind_single)
                              ibuffer = ibuffer + 2
                              fdata(5,imesh,jmesh,kmesh) =    &
     &                           fdata(5,imesh,jmesh,kmesh) + &
     &                           CMPLX(ppm_recvbuffers(ibuffer-1),  &
     &                           ppm_recvbuffers(ibuffer),ppm_kind_single)
#elif  __KIND == __DOUBLE_PRECISION_COMPLEX
                              fdata(1,imesh,jmesh,kmesh) =    &
     &                           fdata(1,imesh,jmesh,kmesh) + &
     &                           CMPLX(ppm_recvbuffers(ibuffer-1),  &
     &                           ppm_recvbuffers(ibuffer),ppm_kind_double)
                              ibuffer = ibuffer + 2
                              fdata(2,imesh,jmesh,kmesh) =    &
     &                           fdata(2,imesh,jmesh,kmesh) + &
     &                           CMPLX(ppm_recvbuffers(ibuffer-1),  &
     &                           ppm_recvbuffers(ibuffer),ppm_kind_double)
                              ibuffer = ibuffer + 2
                              fdata(3,imesh,jmesh,kmesh) =    &
     &                           fdata(3,imesh,jmesh,kmesh) + &
     &                           CMPLX(ppm_recvbuffers(ibuffer-1),  &
     &                           ppm_recvbuffers(ibuffer),ppm_kind_double)
                              ibuffer = ibuffer + 2
                              fdata(4,imesh,jmesh,kmesh) =    &
     &                           fdata(4,imesh,jmesh,kmesh) + &
     &                           CMPLX(ppm_recvbuffers(ibuffer-1),  &
     &                           ppm_recvbuffers(ibuffer),ppm_kind_double)
                              ibuffer = ibuffer + 2
                              fdata(5,imesh,jmesh,kmesh) =    &
     &                           fdata(5,imesh,jmesh,kmesh) + &
     &                           CMPLX(ppm_recvbuffers(ibuffer-1),  &
     &                           ppm_recvbuffers(ibuffer),ppm_kind_double)
#elif  __KIND == __INTEGER
                              fdata(1,imesh,jmesh,kmesh) =    &
     &                           fdata(1,imesh,jmesh,kmesh) + &
     &                           INT(ppm_recvbuffers(ibuffer))
                              ibuffer = ibuffer + 1
                              fdata(2,imesh,jmesh,kmesh) =    &
     &                           fdata(2,imesh,jmesh,kmesh) + &
     &                           INT(ppm_recvbuffers(ibuffer))
                              ibuffer = ibuffer + 1
                              fdata(3,imesh,jmesh,kmesh) =    &
     &                           fdata(3,imesh,jmesh,kmesh) + &
     &                           INT(ppm_recvbuffers(ibuffer))
                              ibuffer = ibuffer + 1
                              fdata(4,imesh,jmesh,kmesh) =    &
     &                           fdata(4,imesh,jmesh,kmesh) + &
     &                           INT(ppm_recvbuffers(ibuffer))
                              ibuffer = ibuffer + 1
                              fdata(5,imesh,jmesh,kmesh) =    &
     &                           fdata(5,imesh,jmesh,kmesh) + &
     &                           INT(ppm_recvbuffers(ibuffer))
#elif  __KIND == __LOGICAL
                              IF (ppm_recvbuffers(ibuffer) .GT.     &
     &                           (1.0_ppm_kind_single-ppm_myepss)) THEN
                                 fdata(1,imesh,jmesh,kmesh) =        &
     &                               (fdata(1,imesh,jmesh,kmesh) .AND. .TRUE.)
                              ELSE
                                 fdata(1,imesh,jmesh,kmesh) =        &
     &                               (fdata(1,imesh,jmesh,kmesh) .AND. .FALSE.)
                              ENDIF
                              ibuffer = ibuffer + 1
                              IF (ppm_recvbuffers(ibuffer) .GT.     &
     &                           (1.0_ppm_kind_single-ppm_myepss)) THEN
                                 fdata(2,imesh,jmesh,kmesh) =        &
     &                               (fdata(2,imesh,jmesh,kmesh) .AND. .TRUE.)
                              ELSE
                                 fdata(2,imesh,jmesh,kmesh) =        &
     &                               (fdata(2,imesh,jmesh,kmesh) .AND. .FALSE.)
                              ENDIF
                              ibuffer = ibuffer + 1
                              IF (ppm_recvbuffers(ibuffer) .GT.     &
     &                           (1.0_ppm_kind_single-ppm_myepss)) THEN
                                 fdata(3,imesh,jmesh,kmesh) =        &
     &                               (fdata(3,imesh,jmesh,kmesh) .AND. .TRUE.)
                              ELSE
                                 fdata(3,imesh,jmesh,kmesh) =        &
     &                               (fdata(3,imesh,jmesh,kmesh) .AND. .FALSE.)
                              ENDIF
                              ibuffer = ibuffer + 1
                              IF (ppm_recvbuffers(ibuffer) .GT.     &
     &                           (1.0_ppm_kind_single-ppm_myepss)) THEN
                                 fdata(4,imesh,jmesh,kmesh) =        &
     &                               (fdata(4,imesh,jmesh,kmesh) .AND. .TRUE.)
                              ELSE
                                 fdata(4,imesh,jmesh,kmesh) =        &
     &                               (fdata(4,imesh,jmesh,kmesh) .AND. .FALSE.)
                              ENDIF
                              ibuffer = ibuffer + 1
                              IF (ppm_recvbuffers(ibuffer) .GT.     &
     &                           (1.0_ppm_kind_single-ppm_myepss)) THEN
                                 fdata(5,imesh,jmesh,kmesh) =        &
     &                               (fdata(5,imesh,jmesh,kmesh) .AND. .TRUE.)
                              ELSE
                                 fdata(5,imesh,jmesh,kmesh) =        &
     &                               (fdata(5,imesh,jmesh,kmesh) .AND. .FALSE.)
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
                                 fdata(k,imesh,jmesh,kmesh) =     &
     &                              REAL(ppm_recvbuffers(ibuffer),ppm_kind_double)
#elif  __KIND == __SINGLE_PRECISION
                                 fdata(k,imesh,jmesh,kmesh) =     &
     &                              ppm_recvbuffers(ibuffer)
#elif  __KIND == __SINGLE_PRECISION_COMPLEX
                                 fdata(k,imesh,jmesh,kmesh) =     &
     &                              CMPLX(ppm_recvbuffers(ibuffer-1),  &
     &                              ppm_recvbuffers(ibuffer),ppm_kind_single)
#elif  __KIND == __DOUBLE_PRECISION_COMPLEX
                                 fdata(k,imesh,jmesh,kmesh) =     &
     &                              CMPLX(ppm_recvbuffers(ibuffer-1),  &
     &                              ppm_recvbuffers(ibuffer),ppm_kind_double)
#elif  __KIND == __INTEGER
                                 fdata(k,imesh,jmesh,kmesh) =     &
     &                              INT(ppm_recvbuffers(ibuffer))
#elif  __KIND == __LOGICAL
                                 IF (ppm_recvbuffers(ibuffer) .GT.     &
     &                              (1.0_ppm_kind_single-ppm_myepss)) THEN
                                    fdata(k,imesh,jmesh,kmesh) = .TRUE.
                                 ELSE
                                    fdata(k,imesh,jmesh,kmesh) = .FALSE.
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
                                 fdata(k,imesh,jmesh,kmesh) =    &
     &                              fdata(k,imesh,jmesh,kmesh) + &
     &                              REAL(ppm_recvbuffers(ibuffer),ppm_kind_double)
#elif  __KIND == __SINGLE_PRECISION
                                 fdata(k,imesh,jmesh,kmesh) =    &
     &                              fdata(k,imesh,jmesh,kmesh) + &
     &                              ppm_recvbuffers(ibuffer)
#elif  __KIND == __SINGLE_PRECISION_COMPLEX
                                 fdata(k,imesh,jmesh,kmesh) =    &
     &                              fdata(k,imesh,jmesh,kmesh) + &
     &                              CMPLX(ppm_recvbuffers(ibuffer-1),  &
     &                              ppm_recvbuffers(ibuffer),ppm_kind_single)
#elif  __KIND == __DOUBLE_PRECISION_COMPLEX
                                 fdata(k,imesh,jmesh,kmesh) =    &
     &                              fdata(k,imesh,jmesh,kmesh) + &
     &                              CMPLX(ppm_recvbuffers(ibuffer-1),  &
     &                              ppm_recvbuffers(ibuffer),ppm_kind_double)
#elif  __KIND == __INTEGER
                                 fdata(k,imesh,jmesh,kmesh) =    &
     &                              fdata(k,imesh,jmesh,kmesh) + &
     &                              INT(ppm_recvbuffers(ibuffer))
#elif  __KIND == __LOGICAL
                                 IF (ppm_recvbuffers(ibuffer) .GT.     &
     &                              (1.0_ppm_kind_single-ppm_myepss)) THEN
                                    fdata(k,imesh,jmesh,kmesh) =        &
     &                                  (fdata(k,imesh,jmesh,kmesh) .AND. .TRUE.)
                                 ELSE
                                    fdata(k,imesh,jmesh,kmesh) =        &
     &                                  (fdata(k,imesh,jmesh,kmesh) .AND. .FALSE.)
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
                        IF (mask(imesh,jmesh,isub)) THEN
                           DO k=1,edim
#if    __KIND == __SINGLE_PRECISION_COMPLEX | __KIND == __DOUBLE_PRECISION_COMPLEX
                              ibuffer = ibuffer + 2
#else
                              ibuffer = ibuffer + 1
#endif
                              IF (rtype .EQ. ppm_param_pop_replace) THEN
#if    __KIND == __DOUBLE_PRECISION
                                 fdata(k,imesh,jmesh,kmesh) =     &
     &                              REAL(ppm_recvbuffers(ibuffer),ppm_kind_double)
#elif  __KIND == __SINGLE_PRECISION
                                 fdata(k,imesh,jmesh,kmesh) =     &
     &                              ppm_recvbuffers(ibuffer)
#elif  __KIND == __SINGLE_PRECISION_COMPLEX
                                 fdata(k,imesh,jmesh,kmesh) =     &
     &                              CMPLX(ppm_recvbuffers(ibuffer-1),  &
     &                              ppm_recvbuffers(ibuffer),ppm_kind_single)
#elif  __KIND == __DOUBLE_PRECISION_COMPLEX
                                 fdata(k,imesh,jmesh,kmesh) =     &
     &                              CMPLX(ppm_recvbuffers(ibuffer-1),  &
     &                              ppm_recvbuffers(ibuffer),ppm_kind_double)
#elif  __KIND == __INTEGER
                                 fdata(k,imesh,jmesh,kmesh) =     &
     &                              INT(ppm_recvbuffers(ibuffer))
#elif  __KIND == __LOGICAL
                                 IF (ppm_recvbuffers(ibuffer) .GT.     &
     &                              (1.0_ppm_kind_single-ppm_myepss)) THEN
                                    fdata(k,imesh,jmesh,kmesh) = .TRUE.
                                 ELSE
                                    fdata(k,imesh,jmesh,kmesh) = .FALSE.
                                 ENDIF
#endif
                              ELSEIF (rtype .EQ. ppm_param_pop_add) THEN
#if    __KIND == __DOUBLE_PRECISION
                                 fdata(k,imesh,jmesh,kmesh) =    &
     &                              fdata(k,imesh,jmesh,kmesh) + &
     &                              REAL(ppm_recvbuffers(ibuffer),ppm_kind_double)
#elif  __KIND == __SINGLE_PRECISION
                                 fdata(k,imesh,jmesh,kmesh) =    &
     &                              fdata(k,imesh,jmesh,kmesh) + &
     &                              ppm_recvbuffers(ibuffer)
#elif  __KIND == __SINGLE_PRECISION_COMPLEX
                                 fdata(k,imesh,jmesh,kmesh) =    &
     &                              fdata(k,imesh,jmesh,kmesh) + &
     &                              CMPLX(ppm_recvbuffers(ibuffer-1),  &
     &                              ppm_recvbuffers(ibuffer),ppm_kind_single)
#elif  __KIND == __DOUBLE_PRECISION_COMPLEX
                                 fdata(k,imesh,jmesh,kmesh) =    &
     &                              fdata(k,imesh,jmesh,kmesh) + &
     &                              CMPLX(ppm_recvbuffers(ibuffer-1),  &
     &                              ppm_recvbuffers(ibuffer),ppm_kind_double)
#elif  __KIND == __INTEGER
                                 fdata(k,imesh,jmesh,kmesh) =    &
     &                              fdata(k,imesh,jmesh,kmesh) + &
     &                              INT(ppm_recvbuffers(ibuffer))
#elif  __KIND == __LOGICAL
                                 IF (ppm_recvbuffers(ibuffer) .GT.     &
     &                              (1.0_ppm_kind_single-ppm_myepss)) THEN
                                    fdata(k,imesh,jmesh,kmesh) =        &
     &                                  (fdata(k,imesh,jmesh,kmesh) .AND. .TRUE.)
                                 ELSE
                                    fdata(k,imesh,jmesh,kmesh) =        &
     &                                  (fdata(k,imesh,jmesh,kmesh) .AND. .FALSE.)
                                 ENDIF
#endif
                              ENDIF
                           ENDDO
                        ELSE
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
                           fdata(imesh,jmesh,kmesh) =     &
     &                        REAL(ppm_recvbuffers(ibuffer),ppm_kind_double)
#elif  __KIND == __SINGLE_PRECISION
                           fdata(imesh,jmesh,kmesh) =     &
     &                        ppm_recvbuffers(ibuffer)
#elif  __KIND == __SINGLE_PRECISION_COMPLEX
                           fdata(imesh,jmesh,kmesh) =     &
     &                        CMPLX(ppm_recvbuffers(ibuffer-1),  &
     &                        ppm_recvbuffers(ibuffer),ppm_kind_single)
#elif  __KIND == __DOUBLE_PRECISION_COMPLEX
                           fdata(imesh,jmesh,kmesh) =     &
     &                        CMPLX(ppm_recvbuffers(ibuffer-1),  &
     &                        ppm_recvbuffers(ibuffer),ppm_kind_double)
#elif  __KIND == __INTEGER
                           fdata(imesh,jmesh,kmesh) =     &
     &                        INT(ppm_recvbuffers(ibuffer))
#elif  __KIND == __LOGICAL
                           IF (ppm_recvbuffers(ibuffer) .GT.     &
     &                        (1.0_ppm_kind_single-ppm_myepss)) THEN
                              fdata(imesh,jmesh,kmesh) = .TRUE.
                           ELSE
                              fdata(imesh,jmesh,kmesh) = .FALSE.
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
                           fdata(imesh,jmesh,kmesh) =    &
     &                        fdata(imesh,jmesh,kmesh) + &
     &                        REAL(ppm_recvbuffers(ibuffer),ppm_kind_double)
#elif  __KIND == __SINGLE_PRECISION
                           fdata(imesh,jmesh,kmesh) =    &
     &                        fdata(imesh,jmesh,kmesh) + &
     &                        ppm_recvbuffers(ibuffer)
#elif  __KIND == __SINGLE_PRECISION_COMPLEX
                           fdata(imesh,jmesh,kmesh) =    &
     &                        fdata(imesh,jmesh,kmesh) + &
     &                        CMPLX(ppm_recvbuffers(ibuffer-1),  &
     &                        ppm_recvbuffers(ibuffer),ppm_kind_single)
#elif  __KIND == __DOUBLE_PRECISION_COMPLEX
                           fdata(imesh,jmesh,kmesh) =    &
     &                        fdata(imesh,jmesh,kmesh) + &
     &                        CMPLX(ppm_recvbuffers(ibuffer-1),  &
     &                        ppm_recvbuffers(ibuffer),ppm_kind_double)
#elif  __KIND == __INTEGER
                           fdata(imesh,jmesh,kmesh) =    &
     &                        fdata(imesh,jmesh,kmesh) + &
     &                        INT(ppm_recvbuffers(ibuffer))
#elif  __KIND == __LOGICAL
                           IF (ppm_recvbuffers(ibuffer) .GT.     &
     &                        (1.0_ppm_kind_single-ppm_myepss)) THEN
                              fdata(imesh,jmesh,kmesh) =        &
     &                            (fdata(imesh,jmesh,kmesh) .AND. .TRUE.)
                           ELSE
                              fdata(imesh,jmesh,kmesh) =        &
     &                            (fdata(imesh,jmesh,kmesh) .AND. .FALSE.)
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
                        IF (mask(imesh,jmesh,isub)) THEN
#if    __KIND == __SINGLE_PRECISION_COMPLEX | __KIND == __DOUBLE_PRECISION_COMPLEX
                           ibuffer = ibuffer + 2
#else
                           ibuffer = ibuffer + 1
#endif
                           IF (rtype .EQ. ppm_param_pop_replace) THEN
#if    __KIND == __DOUBLE_PRECISION
                              fdata(imesh,jmesh,kmesh) =     &
     &                           REAL(ppm_recvbuffers(ibuffer),ppm_kind_double)
#elif  __KIND == __SINGLE_PRECISION
                              fdata(imesh,jmesh,kmesh) =     &
     &                           ppm_recvbuffers(ibuffer)
#elif  __KIND == __SINGLE_PRECISION_COMPLEX
                              fdata(imesh,jmesh,kmesh) =     &
     &                           CMPLX(ppm_recvbuffers(ibuffer-1),  &
     &                           ppm_recvbuffers(ibuffer),ppm_kind_single)
#elif  __KIND == __DOUBLE_PRECISION_COMPLEX
                              fdata(imesh,jmesh,kmesh) =     &
     &                           CMPLX(ppm_recvbuffers(ibuffer-1),  &
     &                           ppm_recvbuffers(ibuffer),ppm_kind_double)
#elif  __KIND == __INTEGER
                              fdata(imesh,jmesh,kmesh) =     &
     &                           INT(ppm_recvbuffers(ibuffer))
#elif  __KIND == __LOGICAL
                              IF (ppm_recvbuffers(ibuffer) .GT.     &
     &                           (1.0_ppm_kind_single-ppm_myepss)) THEN
                                 fdata(imesh,jmesh,kmesh) = .TRUE.
                              ELSE
                                 fdata(imesh,jmesh,kmesh) = .FALSE.
                              ENDIF
#endif
                           ELSEIF (rtype .EQ. ppm_param_pop_add) THEN
#if    __KIND == __DOUBLE_PRECISION
                              fdata(imesh,jmesh,kmesh) =    &
     &                           fdata(imesh,jmesh,kmesh) + &
     &                           REAL(ppm_recvbuffers(ibuffer),ppm_kind_double)
#elif  __KIND == __SINGLE_PRECISION
                              fdata(imesh,jmesh,kmesh) =    &
     &                           fdata(imesh,jmesh,kmesh) + &
     &                           ppm_recvbuffers(ibuffer)
#elif  __KIND == __SINGLE_PRECISION_COMPLEX
                              fdata(imesh,jmesh,kmesh) =    &
     &                           fdata(imesh,jmesh,kmesh) + &
     &                           CMPLX(ppm_recvbuffers(ibuffer-1),  &
     &                           ppm_recvbuffers(ibuffer),ppm_kind_single)
#elif  __KIND == __DOUBLE_PRECISION_COMPLEX
                              fdata(imesh,jmesh,kmesh) =    &
     &                           fdata(imesh,jmesh,kmesh) + &
     &                           CMPLX(ppm_recvbuffers(ibuffer-1),  &
     &                           ppm_recvbuffers(ibuffer),ppm_kind_double)
#elif  __KIND == __INTEGER
                              fdata(imesh,jmesh,kmesh) =    &
     &                           fdata(imesh,jmesh,kmesh) + &
     &                           INT(ppm_recvbuffers(ibuffer))
#elif  __KIND == __LOGICAL
                              IF (ppm_recvbuffers(ibuffer) .GT.     &
     &                           (1.0_ppm_kind_single-ppm_myepss)) THEN
                                 fdata(imesh,jmesh,kmesh) =        &
     &                               (fdata(imesh,jmesh,kmesh) .AND. .TRUE.)
                              ELSE
                                 fdata(imesh,jmesh,kmesh) =        &
     &                               (fdata(imesh,jmesh,kmesh) .AND. .FALSE.)
                              ENDIF
#endif
                           ENDIF
                        ELSE
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
      END SELECT  ! ppm_kind


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
              CALL ppm_error(ppm_err_alloc,caller,     &
     &            'receive buffer PPM_RECVBUFFER',__LINE__,info)
          ENDIF
      ENDIF

      !-------------------------------------------------------------------------
      !  Deallocate inverse sub list
      !-------------------------------------------------------------------------
      iopt   = ppm_param_dealloc
      CALL ppm_alloc(invsublist,ldu,iopt,info)
      or_fail_dealloc("invsublist")

      !-------------------------------------------------------------------------
      !  Return
      !-------------------------------------------------------------------------
      end_subroutine()
      RETURN
      CONTAINS
      SUBROUTINE check
#if   __DIM == __VFIELD
          IF (lda .LT. 1) THEN
              info = ppm_error_error
              CALL ppm_error(ppm_err_argument,caller,  &
     &            'lda must be >0 for vector data',__LINE__,info)
              GOTO 7777
          ENDIF
#elif __DIM == __SFIELD
          IF (lda .NE. 1) THEN
              info = ppm_error_error
              CALL ppm_error(ppm_err_argument,caller,  &
     &            'lda must be =1 for scalar data',__LINE__,info)
              GOTO 7777
          ENDIF
#endif
          IF (PRESENT(poptype)) THEN
              IF ((poptype .NE. ppm_param_pop_replace) .AND.     &
     &            (poptype .NE. ppm_param_pop_add)) THEN
                  info = ppm_error_error
                  CALL ppm_error(ppm_err_argument,caller,  &
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
              zhi = 0
              DO i=1,target_topo%nsublist
                  isub = target_topo%isublist(i)
                  IF (this%nnodes(1,isub).GT.xhi) THEN
                      xhi = this%nnodes(1,isub)
                  ENDIF
                  IF (this%nnodes(2,isub).GT.yhi) THEN
                      yhi = this%nnodes(2,isub)
                  ENDIF
                  IF (this%nnodes(3,isub).GT.zhi) THEN
                      zhi = this%nnodes(3,isub)
                  ENDIF
              ENDDO
              IF (SIZE(mask,1) .LT. xhi) THEN
                  fail("x dimension of mask does not match mesh",exit_point=6666)
              ENDIF
              IF (SIZE(mask,2) .LT. yhi) THEN
                  fail("y dimension of mask does not match mesh",exit_point=6666)
              ENDIF
              IF (SIZE(mask,3) .LT. zhi) THEN
                  fail("z dimension of mask does not match mesh",exit_point=6666)
              ENDIF
          ENDIF
 6666     CONTINUE
     END SUBROUTINE check_two
#if    __DIM == __SFIELD
#if    __KIND == __SINGLE_PRECISION
      END SUBROUTINE mesh_map_pop_3d_sca_s
#elif  __KIND == __DOUBLE_PRECISION
      END SUBROUTINE mesh_map_pop_3d_sca_d
#elif  __KIND == __SINGLE_PRECISION_COMPLEX
      END SUBROUTINE mesh_map_pop_3d_sca_sc
#elif  __KIND == __DOUBLE_PRECISION_COMPLEX
      END SUBROUTINE mesh_map_pop_3d_sca_dc
#elif  __KIND == __INTEGER
      END SUBROUTINE mesh_map_pop_3d_sca_i
#elif  __KIND == __LOGICAL
      END SUBROUTINE mesh_map_pop_3d_sca_l
#endif

#elif  __DIM == __VFIELD
#if    __KIND == __SINGLE_PRECISION
      END SUBROUTINE mesh_map_pop_3d_vec_s
#elif  __KIND == __DOUBLE_PRECISION
      END SUBROUTINE mesh_map_pop_3d_vec_d
#elif  __KIND == __SINGLE_PRECISION_COMPLEX
      END SUBROUTINE mesh_map_pop_3d_vec_sc
#elif  __KIND == __DOUBLE_PRECISION_COMPLEX
      END SUBROUTINE mesh_map_pop_3d_vec_dc
#elif  __KIND == __INTEGER
      END SUBROUTINE mesh_map_pop_3d_vec_i
#elif  __KIND == __LOGICAL
      END SUBROUTINE mesh_map_pop_3d_vec_l
#endif
#endif
