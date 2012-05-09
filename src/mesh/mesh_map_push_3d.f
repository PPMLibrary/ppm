#if    __DIM == __SFIELD
#if    __KIND == __SINGLE_PRECISION
      SUBROUTINE ppm_map_field_push_3d_sca_s(this,fdata_dummy,p_idx,info,mask)
#elif  __KIND == __DOUBLE_PRECISION
      SUBROUTINE ppm_map_field_push_3d_sca_d(this,fdata_dummy,p_idx,info,mask)
#elif  __KIND == __SINGLE_PRECISION_COMPLEX
      SUBROUTINE ppm_map_field_push_3d_sca_sc(this,fdata_dummy,p_idx,info,mask)
#elif  __KIND == __DOUBLE_PRECISION_COMPLEX
      SUBROUTINE ppm_map_field_push_3d_sca_dc(this,fdata_dummy,p_idx,info,mask)
#elif  __KIND == __INTEGER
      SUBROUTINE ppm_map_field_push_3d_sca_i(this,fdata_dummy,p_idx,info,mask)
#elif  __KIND == __LOGICAL
      SUBROUTINE ppm_map_field_push_3d_sca_l(this,fdata_dummy,p_idx,info,mask)
#endif
#elif  __DIM == __VFIELD
#if    __KIND == __SINGLE_PRECISION
      SUBROUTINE ppm_map_field_push_3d_vec_s(this,fdata_dummy,lda,p_idx,info,mask)
#elif  __KIND == __DOUBLE_PRECISION
      SUBROUTINE ppm_map_field_push_3d_vec_d(this,fdata_dummy,lda,p_idx,info,mask)
#elif  __KIND == __SINGLE_PRECISION_COMPLEX
      SUBROUTINE ppm_map_field_push_3d_vec_sc(this,fdata_dummy,lda,p_idx,info,mask)
#elif  __KIND == __DOUBLE_PRECISION_COMPLEX
      SUBROUTINE ppm_map_field_push_3d_vec_dc(this,fdata_dummy,lda,p_idx,info,mask)
#elif  __KIND == __INTEGER
      SUBROUTINE ppm_map_field_push_3d_vec_i(this,fdata_dummy,lda,p_idx,info,mask)
#elif  __KIND == __LOGICAL
      SUBROUTINE ppm_map_field_push_3d_vec_l(this,fdata_dummy,lda,p_idx,info,mask)
      !!! This routine pushes field data onto the send buffer for 3D meshes.
      !!!
      !!! [NOTE]
      !!! ======================================================================
      !!! The on-processor data is stored in the first part of the buffer
      !!!
      !!! If the call to invert_list is too slow here, this
      !!! could be done by topo_store and stored for each
      !!! topology. At the expense of some extra memory.
      !!!
      !!! The fdata arrays are passed as `POINTER` and not
      !!! `INTENT(IN)` even though they are not changed or
      !!! reallocated inside this routine. The reason is that
      !!! LBOUND can be non-standard (i.e. != 1) if we have
      !!! ghost layers and the correct LBOUND is only passed
      !!! along with the array if the latter is given the
      !!! POINTER attribute.
      !!!
      !!! When Ndata is determined, the mask is not
      !!! considered. This could lead to a too large
      !!! sendbuffer being allocated in the last push and
      !!! could be optimized later.
      !!! ======================================================================
      !!!
      !!! [WARNING]
      !!! DIM is the dimension of the fdata array and not the space
      !!! dimension ppm_dim!
#endif
#endif
      !-------------------------------------------------------------------------
      !  Includes
      !-------------------------------------------------------------------------

      !-------------------------------------------------------------------------
      !  Modules
      !-------------------------------------------------------------------------
      USE ppm_module_data_mesh
      USE ppm_module_check_id
      USE ppm_module_util_invert_list
      USE ppm_module_topo_typedef
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
      INTEGER, DIMENSION(3)   :: ldu,mofs,patchid
      INTEGER                 :: i,j,k,ibuffer,isub,imesh,jmesh,kmesh,jsub
      INTEGER                 :: ipatch
      INTEGER                 :: iopt,Ndata,xlo,xhi,ylo,yhi,zlo,zhi,ldb
      LOGICAL                 :: ldo,found_patch
      CHARACTER(LEN=ppm_char) :: mesg
#if   __DIM == __SFIELD
      INTEGER, PARAMETER    :: lda = 1
#endif
      TYPE(ppm_t_topo),      POINTER :: topo => NULL()
      !-------------------------------------------------------------------------
      !  Externals
      !-------------------------------------------------------------------------

      !-------------------------------------------------------------------------
      !  Initialise
      !-------------------------------------------------------------------------
      start_subroutine("mesh_map_push_3d")


      !-------------------------------------------------------------------------
      !  Check arguments
      !-------------------------------------------------------------------------
      IF (ppm_debug .GT. 0) THEN
        CALL check
        IF (info .NE. 0) GOTO 9999
      ENDIF

      topo => ppm_topo(this%topoid)%t
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
     &                       ppm_mesh_isendblksize(2,j)*    &
     &                       ppm_mesh_isendblksize(3,j))
         ENDDO
      ENDDO

      !-------------------------------------------------------------------------
      !  If there is nothing to be sent we are done
      !-------------------------------------------------------------------------
      IF (Ndata .EQ. 0) THEN
          IF (ppm_debug .GT. 0) THEN
              CALL ppm_write(ppm_rank,caller,   &
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
            or_fail_alloc("ppm_buffer_dim")
      CALL ppm_alloc(ppm_buffer_type,ldu,iopt,info)
            or_fail_alloc("ppm_buffer")

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
      ldu(1) = topo%nsublist
      CALL ppm_alloc(sublist,ldu,iopt,info)
            or_fail_alloc("sublist")
      ! We need to copy it into a temp list, since directly using
      ! ppm_isublist(:,ppm_field_topoid) as an argument to invert_list is
      ! not possible since the argument needs to be a POINTER.
      sublist(1:topo%nsublist) =     &
     &    topo%isublist(1:topo%nsublist)
      CALL ppm_util_invert_list(sublist,invsublist,info)
      iopt   = ppm_param_dealloc
      CALL ppm_alloc(sublist,ldu,iopt,info)
            or_fail_alloc("ppm_sublist")

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
            or_fail_alloc("ppm_sendbufferd")

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
               !  Get the patch ID for this mesh block
               !----------------------------------------------------------------
               patchid(1:3) = ppm_mesh_isendpatchid(1:3,j)
               !----------------------------------------------------------------
               !  Translate to local sub ID for storing the data
               !----------------------------------------------------------------
               isub = invsublist(jsub)
               !----------------------------------------------------------------
               !  Mesh offset for this sub
               !----------------------------------------------------------------
               mofs(1) = this%istart(1,jsub)-1
               mofs(2) = this%istart(2,jsub)-1
               mofs(3) = this%istart(3,jsub)-1
               !----------------------------------------------------------------
               !  Get boundaries of mesh block to be sent on local sub
               !  coordinates
               !----------------------------------------------------------------
               xlo = ppm_mesh_isendblkstart(1,j)-mofs(1)
               ylo = ppm_mesh_isendblkstart(2,j)-mofs(2)
               zlo = ppm_mesh_isendblkstart(3,j)-mofs(3)
               xhi = xlo+ppm_mesh_isendblksize(1,j)-1
               yhi = ylo+ppm_mesh_isendblksize(2,j)-1
               zhi = zlo+ppm_mesh_isendblksize(3,j)-1
               IF (ppm_debug .GT. 1) THEN
                   WRITE(mesg,'(A,3I4)') 'start: ',             &
     &                 ppm_mesh_isendblkstart(1,j),ppm_mesh_isendblkstart(2,j),&
     &                 ppm_mesh_isendblkstart(3,j)
                   CALL ppm_write(ppm_rank,caller,mesg,info)
                   WRITE(mesg,'(A,3I4)') 'size: ',             &
     &                 ppm_mesh_isendblksize(1,j),ppm_mesh_isendblksize(2,j),&
     &                 ppm_mesh_isendblksize(3,j)
                   CALL ppm_write(ppm_rank,caller,mesg,info)
                   WRITE(mesg,'(A,3I4)') 'mesh offset: ',mofs(1),mofs(2),mofs(3)
                   CALL ppm_write(ppm_rank,caller,mesg,info)
                   WRITE(mesg,'(A,2I4)') 'xlo, xhi: ',xlo,xhi
                   CALL ppm_write(ppm_rank,caller,mesg,info)
                   WRITE(mesg,'(A,2I4)') 'ylo, yhi: ',ylo,yhi
                   CALL ppm_write(ppm_rank,caller,mesg,info)
                   WRITE(mesg,'(A,2I4)') 'zlo, zhi: ',zlo,zhi
                   CALL ppm_write(ppm_rank,caller,mesg,info)
                   WRITE(mesg,'(A,I1)') 'buffer dim: ',lda
                   CALL ppm_write(ppm_rank,caller,mesg,info)
               ENDIF

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
                   write(*,*) 'IPATCH = ',ipatch,jsub,i,j
                   SELECT TYPE(p => this%subpatch_by_sub(jsub)%vec(ipatch)%t)
                   TYPE IS (ppm_t_subpatch)
                       IF (ALL(p%istart_g.EQ.patchid)) THEN
#if    __DIM == __SFIELD
#if __KIND == __DOUBLE_PRECISION
                            found_patch = .TRUE.
                            fdata => p%subpatch_data%vec(p_idx)%t%data_3d_rd
#endif
#elif  __DIM == __VFIELD
#if __KIND == __DOUBLE_PRECISION
                            found_patch = .TRUE.
                            fdata => p%subpatch_data%vec(p_idx)%t%data_4d_rd
#endif
#endif
                            exit patches
                       ENDIF
                   END SELECT
               ENDDO patches
               IF (.NOT. found_patch) THEN
                   fail("could not find a patch on this sub with the right global id")
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
                  DO kmesh=zlo,zhi
                  DO jmesh=ylo,yhi
                     DO imesh=xlo,xhi
                        ibuffer = ibuffer + 1
#if    __KIND == __SINGLE_PRECISION
                        ppm_sendbufferd(ibuffer) =    &
     &                      REAL(fdata(1,imesh,jmesh,kmesh),ppm_kind_double)
#elif  __KIND == __DOUBLE_PRECISION
                        ppm_sendbufferd(ibuffer) =    &
     &                      fdata(1,imesh,jmesh,kmesh)
#elif  __KIND == __SINGLE_PRECISION_COMPLEX
                        ppm_sendbufferd(ibuffer) =    &
     &                      REAL(fdata(1,imesh,jmesh,kmesh),ppm_kind_double)
                        ibuffer = ibuffer + 1
                        ppm_sendbufferd(ibuffer) =    &
     &                      REAL(AIMAG(fdata(1,imesh,jmesh,kmesh)),   &
     &                      ppm_kind_double)
#elif  __KIND == __DOUBLE_PRECISION_COMPLEX
                        ppm_sendbufferd(ibuffer) =    &
     &                      REAL(fdata(1,imesh,jmesh,kmesh),ppm_kind_double)
                        ibuffer = ibuffer + 1
                        ppm_sendbufferd(ibuffer) =    &
     &                      AIMAG(fdata(1,imesh,jmesh,kmesh))
#elif  __KIND == __INTEGER
                        ppm_sendbufferd(ibuffer) =    &
     &                      REAL(fdata(1,imesh,jmesh,kmesh),ppm_kind_double)
#elif  __KIND == __LOGICAL
                        IF (fdata(1,imesh,jmesh,kmesh)) THEN
                           ppm_sendbufferd(ibuffer) = 1.0_ppm_kind_double
                        ELSE
                           ppm_sendbufferd(ibuffer) = 0.0_ppm_kind_double
                        ENDIF
#endif
                     ENDDO
                  ENDDO
                  ENDDO
               !----------------------------------------------------------------
               !  Unrolled for lda=2
               !----------------------------------------------------------------
               ELSEIF (lda .EQ. 2) THEN
                  DO kmesh=zlo,zhi
                  DO jmesh=ylo,yhi
                     DO imesh=xlo,xhi
                        ibuffer = ibuffer + 1
#if    __KIND == __SINGLE_PRECISION
                        ppm_sendbufferd(ibuffer) =    &
     &                      REAL(fdata(1,imesh,jmesh,kmesh),ppm_kind_double)
                        ibuffer = ibuffer + 1
                        ppm_sendbufferd(ibuffer) =    &
     &                      REAL(fdata(2,imesh,jmesh,kmesh),ppm_kind_double)
#elif  __KIND == __DOUBLE_PRECISION
                        ppm_sendbufferd(ibuffer) =    &
     &                      fdata(1,imesh,jmesh,kmesh)
                        ibuffer = ibuffer + 1
                        ppm_sendbufferd(ibuffer) =    &
     &                      fdata(2,imesh,jmesh,kmesh)
#elif  __KIND == __SINGLE_PRECISION_COMPLEX
                        ppm_sendbufferd(ibuffer) =    &
     &                      REAL(fdata(1,imesh,jmesh,kmesh),ppm_kind_double)
                        ibuffer = ibuffer + 1
                        ppm_sendbufferd(ibuffer) =    &
     &                      REAL(AIMAG(fdata(1,imesh,jmesh,kmesh)),   &
     &                      ppm_kind_double)
                        ibuffer = ibuffer + 1
                        ppm_sendbufferd(ibuffer) =    &
     &                      REAL(fdata(2,imesh,jmesh,kmesh),ppm_kind_double)
                        ibuffer = ibuffer + 1
                        ppm_sendbufferd(ibuffer) =    &
     &                      REAL(AIMAG(fdata(2,imesh,jmesh,kmesh)),   &
     &                      ppm_kind_double)
#elif  __KIND == __DOUBLE_PRECISION_COMPLEX
                        ppm_sendbufferd(ibuffer) =    &
     &                      REAL(fdata(1,imesh,jmesh,kmesh),ppm_kind_double)
                        ibuffer = ibuffer + 1
                        ppm_sendbufferd(ibuffer) =    &
     &                      AIMAG(fdata(1,imesh,jmesh,kmesh))
                        ibuffer = ibuffer + 1
                        ppm_sendbufferd(ibuffer) =    &
     &                      REAL(fdata(2,imesh,jmesh,kmesh),ppm_kind_double)
                        ibuffer = ibuffer + 1
                        ppm_sendbufferd(ibuffer) =    &
     &                      AIMAG(fdata(2,imesh,jmesh,kmesh))
#elif  __KIND == __INTEGER
                        ppm_sendbufferd(ibuffer) =    &
     &                      REAL(fdata(1,imesh,jmesh,kmesh),ppm_kind_double)
                        ibuffer = ibuffer + 1
                        ppm_sendbufferd(ibuffer) =    &
     &                      REAL(fdata(2,imesh,jmesh,kmesh),ppm_kind_double)
#elif  __KIND == __LOGICAL
                        IF (fdata(1,imesh,jmesh,kmesh)) THEN
                           ppm_sendbufferd(ibuffer) = 1.0_ppm_kind_double
                        ELSE
                           ppm_sendbufferd(ibuffer) = 0.0_ppm_kind_double
                        ENDIF
                        ibuffer = ibuffer + 1
                        IF (fdata(2,imesh,jmesh,kmesh)) THEN
                           ppm_sendbufferd(ibuffer) = 1.0_ppm_kind_double
                        ELSE
                           ppm_sendbufferd(ibuffer) = 0.0_ppm_kind_double
                        ENDIF
#endif
                     ENDDO
                  ENDDO
                  ENDDO
               !----------------------------------------------------------------
               !  Unrolled for lda=3
               !----------------------------------------------------------------
               ELSEIF (lda .EQ. 3) THEN
                  DO kmesh=zlo,zhi
                  DO jmesh=ylo,yhi
                     DO imesh=xlo,xhi
                        ibuffer = ibuffer + 1
#if    __KIND == __SINGLE_PRECISION
                        ppm_sendbufferd(ibuffer) =    &
     &                      REAL(fdata(1,imesh,jmesh,kmesh),ppm_kind_double)
                        ibuffer = ibuffer + 1
                        ppm_sendbufferd(ibuffer) =    &
     &                      REAL(fdata(2,imesh,jmesh,kmesh),ppm_kind_double)
                        ibuffer = ibuffer + 1
                        ppm_sendbufferd(ibuffer) =    &
     &                      REAL(fdata(3,imesh,jmesh,kmesh),ppm_kind_double)
#elif  __KIND == __DOUBLE_PRECISION
                        ppm_sendbufferd(ibuffer) =    &
     &                      fdata(1,imesh,jmesh,kmesh)
                        ibuffer = ibuffer + 1
                        ppm_sendbufferd(ibuffer) =    &
     &                      fdata(2,imesh,jmesh,kmesh)
                        ibuffer = ibuffer + 1
                        ppm_sendbufferd(ibuffer) =    &
     &                      fdata(3,imesh,jmesh,kmesh)
#elif  __KIND == __SINGLE_PRECISION_COMPLEX
                        ppm_sendbufferd(ibuffer) =    &
     &                      REAL(fdata(1,imesh,jmesh,kmesh),ppm_kind_double)
                        ibuffer = ibuffer + 1
                        ppm_sendbufferd(ibuffer) =    &
     &                      REAL(AIMAG(fdata(1,imesh,jmesh,kmesh)),   &
     &                      ppm_kind_double)
                        ibuffer = ibuffer + 1
                        ppm_sendbufferd(ibuffer) =    &
     &                      REAL(fdata(2,imesh,jmesh,kmesh),ppm_kind_double)
                        ibuffer = ibuffer + 1
                        ppm_sendbufferd(ibuffer) =    &
     &                      REAL(AIMAG(fdata(2,imesh,jmesh,kmesh)),   &
     &                      ppm_kind_double)
                        ibuffer = ibuffer + 1
                        ppm_sendbufferd(ibuffer) =    &
     &                      REAL(fdata(3,imesh,jmesh,kmesh),ppm_kind_double)
                        ibuffer = ibuffer + 1
                        ppm_sendbufferd(ibuffer) =    &
     &                      REAL(AIMAG(fdata(3,imesh,jmesh,kmesh)),   &
     &                      ppm_kind_double)
#elif  __KIND == __DOUBLE_PRECISION_COMPLEX
                        ppm_sendbufferd(ibuffer) =    &
     &                      REAL(fdata(1,imesh,jmesh,kmesh),ppm_kind_double)
                        ibuffer = ibuffer + 1
                        ppm_sendbufferd(ibuffer) =    &
     &                      AIMAG(fdata(1,imesh,jmesh,kmesh))
                        ibuffer = ibuffer + 1
                        ppm_sendbufferd(ibuffer) =    &
     &                      REAL(fdata(2,imesh,jmesh,kmesh),ppm_kind_double)
                        ibuffer = ibuffer + 1
                        ppm_sendbufferd(ibuffer) =    &
     &                      AIMAG(fdata(2,imesh,jmesh,kmesh))
                        ibuffer = ibuffer + 1
                        ppm_sendbufferd(ibuffer) =    &
     &                      REAL(fdata(3,imesh,jmesh,kmesh),ppm_kind_double)
                        ibuffer = ibuffer + 1
                        ppm_sendbufferd(ibuffer) =    &
     &                      AIMAG(fdata(3,imesh,jmesh,kmesh))
#elif  __KIND == __INTEGER
                        ppm_sendbufferd(ibuffer) =    &
     &                      REAL(fdata(1,imesh,jmesh,kmesh),ppm_kind_double)
                        ibuffer = ibuffer + 1
                        ppm_sendbufferd(ibuffer) =    &
     &                      REAL(fdata(2,imesh,jmesh,kmesh),ppm_kind_double)
                        ibuffer = ibuffer + 1
                        ppm_sendbufferd(ibuffer) =    &
     &                      REAL(fdata(3,imesh,jmesh,kmesh),ppm_kind_double)
#elif  __KIND == __LOGICAL
                        IF (fdata(1,imesh,jmesh,kmesh)) THEN
                           ppm_sendbufferd(ibuffer) = 1.0_ppm_kind_double
                        ELSE
                           ppm_sendbufferd(ibuffer) = 0.0_ppm_kind_double
                        ENDIF
                        ibuffer = ibuffer + 1
                        IF (fdata(2,imesh,jmesh,kmesh)) THEN
                           ppm_sendbufferd(ibuffer) = 1.0_ppm_kind_double
                        ELSE
                           ppm_sendbufferd(ibuffer) = 0.0_ppm_kind_double
                        ENDIF
                        ibuffer = ibuffer + 1
                        IF (fdata(3,imesh,jmesh,kmesh)) THEN
                           ppm_sendbufferd(ibuffer) = 1.0_ppm_kind_double
                        ELSE
                           ppm_sendbufferd(ibuffer) = 0.0_ppm_kind_double
                        ENDIF
#endif
                     ENDDO
                  ENDDO
                  ENDDO
               !----------------------------------------------------------------
               !  Unrolled for lda=4
               !----------------------------------------------------------------
               ELSEIF (lda .EQ. 4) THEN
                  DO kmesh=zlo,zhi
                  DO jmesh=ylo,yhi
                     DO imesh=xlo,xhi
                        ibuffer = ibuffer + 1
#if    __KIND == __SINGLE_PRECISION
                        ppm_sendbufferd(ibuffer) =    &
     &                      REAL(fdata(1,imesh,jmesh,kmesh),ppm_kind_double)
                        ibuffer = ibuffer + 1
                        ppm_sendbufferd(ibuffer) =    &
     &                      REAL(fdata(2,imesh,jmesh,kmesh),ppm_kind_double)
                        ibuffer = ibuffer + 1
                        ppm_sendbufferd(ibuffer) =    &
     &                      REAL(fdata(3,imesh,jmesh,kmesh),ppm_kind_double)
                        ibuffer = ibuffer + 1
                        ppm_sendbufferd(ibuffer) =    &
     &                      REAL(fdata(4,imesh,jmesh,kmesh),ppm_kind_double)
#elif  __KIND == __DOUBLE_PRECISION
                        ppm_sendbufferd(ibuffer) =    &
     &                      fdata(1,imesh,jmesh,kmesh)
                        ibuffer = ibuffer + 1
                        ppm_sendbufferd(ibuffer) =    &
     &                      fdata(2,imesh,jmesh,kmesh)
                        ibuffer = ibuffer + 1
                        ppm_sendbufferd(ibuffer) =    &
     &                      fdata(3,imesh,jmesh,kmesh)
                        ibuffer = ibuffer + 1
                        ppm_sendbufferd(ibuffer) =    &
     &                      fdata(4,imesh,jmesh,kmesh)
#elif  __KIND == __SINGLE_PRECISION_COMPLEX
                        ppm_sendbufferd(ibuffer) =    &
     &                      REAL(fdata(1,imesh,jmesh,kmesh),ppm_kind_double)
                        ibuffer = ibuffer + 1
                        ppm_sendbufferd(ibuffer) =    &
     &                      REAL(AIMAG(fdata(1,imesh,jmesh,kmesh)),   &
     &                      ppm_kind_double)
                        ibuffer = ibuffer + 1
                        ppm_sendbufferd(ibuffer) =    &
     &                      REAL(fdata(2,imesh,jmesh,kmesh),ppm_kind_double)
                        ibuffer = ibuffer + 1
                        ppm_sendbufferd(ibuffer) =    &
     &                      REAL(AIMAG(fdata(2,imesh,jmesh,kmesh)),   &
     &                      ppm_kind_double)
                        ibuffer = ibuffer + 1
                        ppm_sendbufferd(ibuffer) =    &
     &                      REAL(fdata(3,imesh,jmesh,kmesh),ppm_kind_double)
                        ibuffer = ibuffer + 1
                        ppm_sendbufferd(ibuffer) =    &
     &                      REAL(AIMAG(fdata(3,imesh,jmesh,kmesh)),   &
     &                      ppm_kind_double)
                        ibuffer = ibuffer + 1
                        ppm_sendbufferd(ibuffer) =    &
     &                      REAL(fdata(4,imesh,jmesh,kmesh),ppm_kind_double)
                        ibuffer = ibuffer + 1
                        ppm_sendbufferd(ibuffer) =    &
     &                      REAL(AIMAG(fdata(4,imesh,jmesh,kmesh)),   &
     &                      ppm_kind_double)
#elif  __KIND == __DOUBLE_PRECISION_COMPLEX
                        ppm_sendbufferd(ibuffer) =    &
     &                      REAL(fdata(1,imesh,jmesh,kmesh),ppm_kind_double)
                        ibuffer = ibuffer + 1
                        ppm_sendbufferd(ibuffer) =    &
     &                      AIMAG(fdata(1,imesh,jmesh,kmesh))
                        ibuffer = ibuffer + 1
                        ppm_sendbufferd(ibuffer) =    &
     &                      REAL(fdata(2,imesh,jmesh,kmesh),ppm_kind_double)
                        ibuffer = ibuffer + 1
                        ppm_sendbufferd(ibuffer) =    &
     &                      AIMAG(fdata(2,imesh,jmesh,kmesh))
                        ibuffer = ibuffer + 1
                        ppm_sendbufferd(ibuffer) =    &
     &                      REAL(fdata(3,imesh,jmesh,kmesh),ppm_kind_double)
                        ibuffer = ibuffer + 1
                        ppm_sendbufferd(ibuffer) =    &
     &                      AIMAG(fdata(3,imesh,jmesh,kmesh))
                        ibuffer = ibuffer + 1
                        ppm_sendbufferd(ibuffer) =    &
     &                      REAL(fdata(4,imesh,jmesh,kmesh),ppm_kind_double)
                        ibuffer = ibuffer + 1
                        ppm_sendbufferd(ibuffer) =    &
     &                      AIMAG(fdata(4,imesh,jmesh,kmesh))
#elif  __KIND == __INTEGER
                        ppm_sendbufferd(ibuffer) =    &
     &                      REAL(fdata(1,imesh,jmesh,kmesh),ppm_kind_double)
                        ibuffer = ibuffer + 1
                        ppm_sendbufferd(ibuffer) =    &
     &                      REAL(fdata(2,imesh,jmesh,kmesh),ppm_kind_double)
                        ibuffer = ibuffer + 1
                        ppm_sendbufferd(ibuffer) =    &
     &                      REAL(fdata(3,imesh,jmesh,kmesh),ppm_kind_double)
                        ibuffer = ibuffer + 1
                        ppm_sendbufferd(ibuffer) =    &
     &                      REAL(fdata(4,imesh,jmesh,kmesh),ppm_kind_double)
#elif  __KIND == __LOGICAL
                        IF (fdata(1,imesh,jmesh,kmesh)) THEN
                           ppm_sendbufferd(ibuffer) = 1.0_ppm_kind_double
                        ELSE
                           ppm_sendbufferd(ibuffer) = 0.0_ppm_kind_double
                        ENDIF
                        ibuffer = ibuffer + 1
                        IF (fdata(2,imesh,jmesh,kmesh)) THEN
                           ppm_sendbufferd(ibuffer) = 1.0_ppm_kind_double
                        ELSE
                           ppm_sendbufferd(ibuffer) = 0.0_ppm_kind_double
                        ENDIF
                        ibuffer = ibuffer + 1
                        IF (fdata(3,imesh,jmesh,kmesh)) THEN
                           ppm_sendbufferd(ibuffer) = 1.0_ppm_kind_double
                        ELSE
                           ppm_sendbufferd(ibuffer) = 0.0_ppm_kind_double
                        ENDIF
                        ibuffer = ibuffer + 1
                        IF (fdata(4,imesh,jmesh,kmesh)) THEN
                           ppm_sendbufferd(ibuffer) = 1.0_ppm_kind_double
                        ELSE
                           ppm_sendbufferd(ibuffer) = 0.0_ppm_kind_double
                        ENDIF
#endif
                     ENDDO
                  ENDDO
                  ENDDO
               !----------------------------------------------------------------
               !  Unrolled for lda=5
               !----------------------------------------------------------------
               ELSEIF (lda .EQ. 5) THEN
                  DO kmesh=zlo,zhi
                  DO jmesh=ylo,yhi
                     DO imesh=xlo,xhi
                        ibuffer = ibuffer + 1
#if    __KIND == __SINGLE_PRECISION
                        ppm_sendbufferd(ibuffer) =    &
     &                      REAL(fdata(1,imesh,jmesh,kmesh),ppm_kind_double)
                        ibuffer = ibuffer + 1
                        ppm_sendbufferd(ibuffer) =    &
     &                      REAL(fdata(2,imesh,jmesh,kmesh),ppm_kind_double)
                        ibuffer = ibuffer + 1
                        ppm_sendbufferd(ibuffer) =    &
     &                      REAL(fdata(3,imesh,jmesh,kmesh),ppm_kind_double)
                        ibuffer = ibuffer + 1
                        ppm_sendbufferd(ibuffer) =    &
     &                      REAL(fdata(4,imesh,jmesh,kmesh),ppm_kind_double)
                        ibuffer = ibuffer + 1
                        ppm_sendbufferd(ibuffer) =    &
     &                      REAL(fdata(5,imesh,jmesh,kmesh),ppm_kind_double)
#elif  __KIND == __DOUBLE_PRECISION
                        ppm_sendbufferd(ibuffer) =    &
     &                      fdata(1,imesh,jmesh,kmesh)
                        ibuffer = ibuffer + 1
                        ppm_sendbufferd(ibuffer) =    &
     &                      fdata(2,imesh,jmesh,kmesh)
                        ibuffer = ibuffer + 1
                        ppm_sendbufferd(ibuffer) =    &
     &                      fdata(3,imesh,jmesh,kmesh)
                        ibuffer = ibuffer + 1
                        ppm_sendbufferd(ibuffer) =    &
     &                      fdata(4,imesh,jmesh,kmesh)
                        ibuffer = ibuffer + 1
                        ppm_sendbufferd(ibuffer) =    &
     &                      fdata(5,imesh,jmesh,kmesh)
#elif  __KIND == __SINGLE_PRECISION_COMPLEX
                        ppm_sendbufferd(ibuffer) =    &
     &                      REAL(fdata(1,imesh,jmesh,kmesh),ppm_kind_double)
                        ibuffer = ibuffer + 1
                        ppm_sendbufferd(ibuffer) =    &
     &                      REAL(AIMAG(fdata(1,imesh,jmesh,kmesh)),   &
     &                      ppm_kind_double)
                        ibuffer = ibuffer + 1
                        ppm_sendbufferd(ibuffer) =    &
     &                      REAL(fdata(2,imesh,jmesh,kmesh),ppm_kind_double)
                        ibuffer = ibuffer + 1
                        ppm_sendbufferd(ibuffer) =    &
     &                      REAL(AIMAG(fdata(2,imesh,jmesh,kmesh)),   &
     &                      ppm_kind_double)
                        ibuffer = ibuffer + 1
                        ppm_sendbufferd(ibuffer) =    &
     &                      REAL(fdata(3,imesh,jmesh,kmesh),ppm_kind_double)
                        ibuffer = ibuffer + 1
                        ppm_sendbufferd(ibuffer) =    &
     &                      REAL(AIMAG(fdata(3,imesh,jmesh,kmesh)),   &
     &                      ppm_kind_double)
                        ibuffer = ibuffer + 1
                        ppm_sendbufferd(ibuffer) =    &
     &                      REAL(fdata(4,imesh,jmesh,kmesh),ppm_kind_double)
                        ibuffer = ibuffer + 1
                        ppm_sendbufferd(ibuffer) =    &
     &                      REAL(AIMAG(fdata(4,imesh,jmesh,kmesh)),   &
     &                      ppm_kind_double)
                        ibuffer = ibuffer + 1
                        ppm_sendbufferd(ibuffer) =    &
     &                      REAL(fdata(5,imesh,jmesh,kmesh),ppm_kind_double)
                        ibuffer = ibuffer + 1
                        ppm_sendbufferd(ibuffer) =    &
     &                      REAL(AIMAG(fdata(5,imesh,jmesh,kmesh)),   &
     &                      ppm_kind_double)
#elif  __KIND == __DOUBLE_PRECISION_COMPLEX
                        ppm_sendbufferd(ibuffer) =    &
     &                      REAL(fdata(1,imesh,jmesh,kmesh),ppm_kind_double)
                        ibuffer = ibuffer + 1
                        ppm_sendbufferd(ibuffer) =    &
     &                      AIMAG(fdata(1,imesh,jmesh,kmesh))
                        ibuffer = ibuffer + 1
                        ppm_sendbufferd(ibuffer) =    &
     &                      REAL(fdata(2,imesh,jmesh,kmesh),ppm_kind_double)
                        ibuffer = ibuffer + 1
                        ppm_sendbufferd(ibuffer) =    &
     &                      AIMAG(fdata(2,imesh,jmesh,kmesh))
                        ibuffer = ibuffer + 1
                        ppm_sendbufferd(ibuffer) =    &
     &                      REAL(fdata(3,imesh,jmesh,kmesh),ppm_kind_double)
                        ibuffer = ibuffer + 1
                        ppm_sendbufferd(ibuffer) =    &
     &                      AIMAG(fdata(3,imesh,jmesh,kmesh))
                        ibuffer = ibuffer + 1
                        ppm_sendbufferd(ibuffer) =    &
     &                      REAL(fdata(4,imesh,jmesh,kmesh),ppm_kind_double)
                        ibuffer = ibuffer + 1
                        ppm_sendbufferd(ibuffer) =    &
     &                      AIMAG(fdata(4,imesh,jmesh,kmesh))
                        ibuffer = ibuffer + 1
                        ppm_sendbufferd(ibuffer) =    &
     &                      REAL(fdata(5,imesh,jmesh,kmesh),ppm_kind_double)
                        ibuffer = ibuffer + 1
                        ppm_sendbufferd(ibuffer) =    &
     &                      AIMAG(fdata(5,imesh,jmesh,kmesh))
#elif  __KIND == __INTEGER
                        ppm_sendbufferd(ibuffer) =    &
     &                      REAL(fdata(1,imesh,jmesh,kmesh),ppm_kind_double)
                        ibuffer = ibuffer + 1
                        ppm_sendbufferd(ibuffer) =    &
     &                      REAL(fdata(2,imesh,jmesh,kmesh),ppm_kind_double)
                        ibuffer = ibuffer + 1
                        ppm_sendbufferd(ibuffer) =    &
     &                      REAL(fdata(3,imesh,jmesh,kmesh),ppm_kind_double)
                        ibuffer = ibuffer + 1
                        ppm_sendbufferd(ibuffer) =    &
     &                      REAL(fdata(4,imesh,jmesh,kmesh),ppm_kind_double)
                        ibuffer = ibuffer + 1
                        ppm_sendbufferd(ibuffer) =    &
     &                      REAL(fdata(5,imesh,jmesh,kmesh),ppm_kind_double)
#elif  __KIND == __LOGICAL
                        IF (fdata(1,imesh,jmesh,kmesh)) THEN
                           ppm_sendbufferd(ibuffer) = 1.0_ppm_kind_double
                        ELSE
                           ppm_sendbufferd(ibuffer) = 0.0_ppm_kind_double
                        ENDIF
                        ibuffer = ibuffer + 1
                        IF (fdata(2,imesh,jmesh,kmesh)) THEN
                           ppm_sendbufferd(ibuffer) = 1.0_ppm_kind_double
                        ELSE
                           ppm_sendbufferd(ibuffer) = 0.0_ppm_kind_double
                        ENDIF
                        ibuffer = ibuffer + 1
                        IF (fdata(3,imesh,jmesh,kmesh)) THEN
                           ppm_sendbufferd(ibuffer) = 1.0_ppm_kind_double
                        ELSE
                           ppm_sendbufferd(ibuffer) = 0.0_ppm_kind_double
                        ENDIF
                        ibuffer = ibuffer + 1
                        IF (fdata(4,imesh,jmesh,kmesh)) THEN
                           ppm_sendbufferd(ibuffer) = 1.0_ppm_kind_double
                        ELSE
                           ppm_sendbufferd(ibuffer) = 0.0_ppm_kind_double
                        ENDIF
                        ibuffer = ibuffer + 1
                        IF (fdata(5,imesh,jmesh,kmesh)) THEN
                           ppm_sendbufferd(ibuffer) = 1.0_ppm_kind_double
                        ELSE
                           ppm_sendbufferd(ibuffer) = 0.0_ppm_kind_double
                        ENDIF
#endif
                     ENDDO
                  ENDDO
                  ENDDO
               !----------------------------------------------------------------
               !  For lda.GT.5 vector length is lda
               !----------------------------------------------------------------
               ELSE
                  DO kmesh=zlo,zhi
                  DO jmesh=ylo,yhi
                     DO imesh=xlo,xhi
                        DO k=1,lda
                           ibuffer = ibuffer + 1
#if    __KIND == __SINGLE_PRECISION
                           ppm_sendbufferd(ibuffer) =    &
     &                        REAL(fdata(k,imesh,jmesh,kmesh),ppm_kind_double)
#elif  __KIND == __DOUBLE_PRECISION
                           ppm_sendbufferd(ibuffer) =    &
     &                         fdata(k,imesh,jmesh,kmesh)
#elif  __KIND == __SINGLE_PRECISION_COMPLEX
                           ppm_sendbufferd(ibuffer) =    &
     &                        REAL(fdata(k,imesh,jmesh,kmesh),ppm_kind_double)
                           ibuffer = ibuffer + 1
                           ppm_sendbufferd(ibuffer) =    &
     &                         REAL(AIMAG(fdata(k,imesh,jmesh,kmesh)),   &
     &                         ppm_kind_double)
#elif  __KIND == __DOUBLE_PRECISION_COMPLEX
                           ppm_sendbufferd(ibuffer) =    &
     &                        REAL(fdata(k,imesh,jmesh,kmesh),ppm_kind_double)
                           ibuffer = ibuffer + 1
                           ppm_sendbufferd(ibuffer) =    &
     &                         AIMAG(fdata(k,imesh,jmesh,kmesh))
#elif  __KIND == __INTEGER
                           ppm_sendbufferd(ibuffer) =    &
     &                        REAL(fdata(k,imesh,jmesh,kmesh),ppm_kind_double)
#elif  __KIND == __LOGICAL
                           IF (fdata(k,imesh,jmesh,kmesh)) THEN
                              ppm_sendbufferd(ibuffer) = 1.0_ppm_kind_double
                           ELSE
                              ppm_sendbufferd(ibuffer) = 0.0_ppm_kind_double
                           ENDIF
#endif
                        ENDDO
                     ENDDO
                  ENDDO
                  ENDDO
               ENDIF
#elif  __DIM == __SFIELD
               !----------------------------------------------------------------
               !  Scalar version
               !----------------------------------------------------------------
                  DO kmesh=zlo,zhi
               DO jmesh=ylo,yhi
                  DO imesh=xlo,xhi
                     ibuffer = ibuffer + 1
#if    __KIND == __SINGLE_PRECISION
                     ppm_sendbufferd(ibuffer) =    &
     &                   REAL(fdata(imesh,jmesh,kmesh),ppm_kind_double)
#elif  __KIND == __DOUBLE_PRECISION
                     ppm_sendbufferd(ibuffer) = fdata(imesh,jmesh,kmesh)
#elif  __KIND == __SINGLE_PRECISION_COMPLEX
                     ppm_sendbufferd(ibuffer) =    &
     &                   REAL(fdata(imesh,jmesh,kmesh),ppm_kind_double)
                     ibuffer = ibuffer + 1
                     ppm_sendbufferd(ibuffer) =    &
     &                   REAL(AIMAG(fdata(imesh,jmesh,kmesh)),   &
     &                   ppm_kind_double)
#elif  __KIND == __DOUBLE_PRECISION_COMPLEX
                     ppm_sendbufferd(ibuffer) =    &
     &                   REAL(fdata(imesh,jmesh,kmesh),ppm_kind_double)
                     ibuffer = ibuffer + 1
                     ppm_sendbufferd(ibuffer) = AIMAG(fdata(imesh,jmesh,kmesh))
#elif  __KIND == __INTEGER
                     ppm_sendbufferd(ibuffer) =    &
     &                   REAL(fdata(imesh,jmesh,kmesh),ppm_kind_double)
#elif  __KIND == __LOGICAL
                     IF (fdata(imesh,jmesh,kmesh)) THEN
                        ppm_sendbufferd(ibuffer) = 1.0_ppm_kind_double
                     ELSE
                        ppm_sendbufferd(ibuffer) = 0.0_ppm_kind_double
                     ENDIF
#endif
                  ENDDO
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
             CALL ppm_error(ppm_err_alloc,caller,     &
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
               mofs(1) = this%istart(1,jsub)-1
               mofs(2) = this%istart(2,jsub)-1
               mofs(3) = this%istart(3,jsub)-1
               !----------------------------------------------------------------
               !  Get boundaries of mesh block to be sent on local sub
               !  coordinates
               !----------------------------------------------------------------
               xlo = ppm_mesh_isendblkstart(1,j)-mofs(1)
               ylo = ppm_mesh_isendblkstart(2,j)-mofs(2)
               zlo = ppm_mesh_isendblkstart(3,j)-mofs(3)
               xhi = xlo+ppm_mesh_isendblksize(1,j)-1
               yhi = ylo+ppm_mesh_isendblksize(2,j)-1
               zhi = zlo+ppm_mesh_isendblksize(3,j)-1
               IF (ppm_debug .GT. 1) THEN
                   WRITE(mesg,'(A,3I4)') 'start: ',             &
     &                 ppm_mesh_isendblkstart(1,j),ppm_mesh_isendblkstart(2,j),&
     &                 ppm_mesh_isendblkstart(3,j)
                   CALL ppm_write(ppm_rank,caller,mesg,info)
                   WRITE(mesg,'(A,3I4)') 'size: ',             &
     &                 ppm_mesh_isendblksize(1,j),ppm_mesh_isendblksize(2,j),&
     &                 ppm_mesh_isendblksize(3,j)
                   CALL ppm_write(ppm_rank,caller,mesg,info)
                   WRITE(mesg,'(A,3I4)') 'mesh offset: ',mofs(1),mofs(2),mofs(3)
                   CALL ppm_write(ppm_rank,caller,mesg,info)
                   WRITE(mesg,'(A,2I4)') 'xlo, xhi: ',xlo,xhi
                   CALL ppm_write(ppm_rank,caller,mesg,info)
                   WRITE(mesg,'(A,2I4)') 'ylo, yhi: ',ylo,yhi
                   CALL ppm_write(ppm_rank,caller,mesg,info)
                   WRITE(mesg,'(A,2I4)') 'zlo, zhi: ',zlo,zhi
                   CALL ppm_write(ppm_rank,caller,mesg,info)
                   WRITE(mesg,'(A,I1)') 'buffer dim: ',lda
                   CALL ppm_write(ppm_rank,caller,mesg,info)
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
                  DO kmesh=zlo,zhi
                  DO jmesh=ylo,yhi
                     DO imesh=xlo,xhi
                        ibuffer = ibuffer + 1
#if    __KIND == __DOUBLE_PRECISION
                        ppm_sendbuffers(ibuffer) =    &
     &                     REAL(fdata(1,imesh,jmesh,kmesh),ppm_kind_single)
#elif  __KIND == __SINGLE_PRECISION
                        ppm_sendbuffers(ibuffer) =    &
     &                     fdata(1,imesh,jmesh,kmesh)
#elif  __KIND == __SINGLE_PRECISION_COMPLEX
                        ppm_sendbuffers(ibuffer) =    &
     &                     REAL(fdata(1,imesh,jmesh,kmesh),ppm_kind_single)
                        ibuffer = ibuffer + 1
                        ppm_sendbuffers(ibuffer) =    &
     &                     AIMAG(fdata(1,imesh,jmesh,kmesh))
#elif  __KIND == __DOUBLE_PRECISION_COMPLEX
                        ppm_sendbuffers(ibuffer) =    &
     &                      REAL(fdata(1,imesh,jmesh,kmesh),ppm_kind_single)
                        ibuffer = ibuffer + 1
                        ppm_sendbuffers(ibuffer) =    &
     &                      REAL(AIMAG(fdata(1,imesh,jmesh,kmesh)),   &
     &                      ppm_kind_single)
#elif  __KIND == __INTEGER
                        ppm_sendbuffers(ibuffer) =    &
     &                      REAL(fdata(1,imesh,jmesh,kmesh),ppm_kind_single)
#elif  __KIND == __LOGICAL
                        IF (fdata(1,imesh,jmesh,kmesh)) THEN
                           ppm_sendbuffers(ibuffer) = 1.0_ppm_kind_single
                        ELSE
                           ppm_sendbuffers(ibuffer) = 0.0_ppm_kind_single
                        ENDIF
#endif
                     ENDDO
                  ENDDO
                  ENDDO
               !----------------------------------------------------------------
               !  Unrolled for lda=2
               !----------------------------------------------------------------
               ELSEIF (lda .EQ. 2) THEN
                  DO kmesh=zlo,zhi
                  DO jmesh=ylo,yhi
                     DO imesh=xlo,xhi
                        ibuffer = ibuffer + 1
#if    __KIND == __DOUBLE_PRECISION
                        ppm_sendbuffers(ibuffer) =    &
     &                     REAL(fdata(1,imesh,jmesh,kmesh),ppm_kind_single)
                        ibuffer = ibuffer + 1
                        ppm_sendbuffers(ibuffer) =    &
     &                     REAL(fdata(2,imesh,jmesh,kmesh),ppm_kind_single)
#elif  __KIND == __SINGLE_PRECISION
                        ppm_sendbuffers(ibuffer) =    &
     &                     fdata(1,imesh,jmesh,kmesh)
                        ibuffer = ibuffer + 1
                        ppm_sendbuffers(ibuffer) =    &
     &                     fdata(2,imesh,jmesh,kmesh)
#elif  __KIND == __SINGLE_PRECISION_COMPLEX
                        ppm_sendbuffers(ibuffer) =    &
     &                     REAL(fdata(1,imesh,jmesh,kmesh),ppm_kind_single)
                        ibuffer = ibuffer + 1
                        ppm_sendbuffers(ibuffer) =    &
     &                     AIMAG(fdata(1,imesh,jmesh,kmesh))
                        ibuffer = ibuffer + 1
                        ppm_sendbuffers(ibuffer) =    &
     &                     REAL(fdata(2,imesh,jmesh,kmesh),ppm_kind_single)
                        ibuffer = ibuffer + 1
                        ppm_sendbuffers(ibuffer) =    &
     &                     AIMAG(fdata(2,imesh,jmesh,kmesh))
#elif  __KIND == __DOUBLE_PRECISION_COMPLEX
                        ppm_sendbuffers(ibuffer) =    &
     &                      REAL(fdata(1,imesh,jmesh,kmesh),ppm_kind_single)
                        ibuffer = ibuffer + 1
                        ppm_sendbuffers(ibuffer) =    &
     &                      REAL(AIMAG(fdata(1,imesh,jmesh,kmesh)),   &
     &                      ppm_kind_single)
                        ibuffer = ibuffer + 1
                        ppm_sendbuffers(ibuffer) =    &
     &                      REAL(fdata(2,imesh,jmesh,kmesh),ppm_kind_single)
                        ibuffer = ibuffer + 1
                        ppm_sendbuffers(ibuffer) =    &
     &                      REAL(AIMAG(fdata(2,imesh,jmesh,kmesh)),   &
     &                      ppm_kind_single)
#elif  __KIND == __INTEGER
                        ppm_sendbuffers(ibuffer) =    &
     &                      REAL(fdata(1,imesh,jmesh,kmesh),ppm_kind_single)
                        ibuffer = ibuffer + 1
                        ppm_sendbuffers(ibuffer) =    &
     &                      REAL(fdata(2,imesh,jmesh,kmesh),ppm_kind_single)
#elif  __KIND == __LOGICAL
                        IF (fdata(1,imesh,jmesh,kmesh)) THEN
                           ppm_sendbuffers(ibuffer) = 1.0_ppm_kind_single
                        ELSE
                           ppm_sendbuffers(ibuffer) = 0.0_ppm_kind_single
                        ENDIF
                        ibuffer = ibuffer + 1
                        IF (fdata(2,imesh,jmesh,kmesh)) THEN
                           ppm_sendbuffers(ibuffer) = 1.0_ppm_kind_single
                        ELSE
                           ppm_sendbuffers(ibuffer) = 0.0_ppm_kind_single
                        ENDIF
#endif
                     ENDDO
                  ENDDO
                  ENDDO
               !----------------------------------------------------------------
               !  Unrolled for lda=3
               !----------------------------------------------------------------
               ELSEIF (lda .EQ. 3) THEN
                  DO kmesh=zlo,zhi
                  DO jmesh=ylo,yhi
                     DO imesh=xlo,xhi
                        ibuffer = ibuffer + 1
#if    __KIND == __DOUBLE_PRECISION
                        ppm_sendbuffers(ibuffer) =    &
     &                     REAL(fdata(1,imesh,jmesh,kmesh),ppm_kind_single)
                        ibuffer = ibuffer + 1
                        ppm_sendbuffers(ibuffer) =    &
     &                     REAL(fdata(2,imesh,jmesh,kmesh),ppm_kind_single)
                        ibuffer = ibuffer + 1
                        ppm_sendbuffers(ibuffer) =    &
     &                     REAL(fdata(3,imesh,jmesh,kmesh),ppm_kind_single)
#elif  __KIND == __SINGLE_PRECISION
                        ppm_sendbuffers(ibuffer) =    &
     &                     fdata(1,imesh,jmesh,kmesh)
                        ibuffer = ibuffer + 1
                        ppm_sendbuffers(ibuffer) =    &
     &                     fdata(2,imesh,jmesh,kmesh)
                        ibuffer = ibuffer + 1
                        ppm_sendbuffers(ibuffer) =    &
     &                     fdata(3,imesh,jmesh,kmesh)
#elif  __KIND == __SINGLE_PRECISION_COMPLEX
                        ppm_sendbuffers(ibuffer) =    &
     &                     REAL(fdata(1,imesh,jmesh,kmesh),ppm_kind_single)
                        ibuffer = ibuffer + 1
                        ppm_sendbuffers(ibuffer) =    &
     &                     AIMAG(fdata(1,imesh,jmesh,kmesh))
                        ibuffer = ibuffer + 1
                        ppm_sendbuffers(ibuffer) =    &
     &                     REAL(fdata(2,imesh,jmesh,kmesh),ppm_kind_single)
                        ibuffer = ibuffer + 1
                        ppm_sendbuffers(ibuffer) =    &
     &                     AIMAG(fdata(2,imesh,jmesh,kmesh))
                        ibuffer = ibuffer + 1
                        ppm_sendbuffers(ibuffer) =    &
     &                     REAL(fdata(3,imesh,jmesh,kmesh),ppm_kind_single)
                        ibuffer = ibuffer + 1
                        ppm_sendbuffers(ibuffer) =    &
     &                     AIMAG(fdata(3,imesh,jmesh,kmesh))
#elif  __KIND == __DOUBLE_PRECISION_COMPLEX
                        ppm_sendbuffers(ibuffer) =    &
     &                      REAL(fdata(1,imesh,jmesh,kmesh),ppm_kind_single)
                        ibuffer = ibuffer + 1
                        ppm_sendbuffers(ibuffer) =    &
     &                      REAL(AIMAG(fdata(1,imesh,jmesh,kmesh)),   &
     &                      ppm_kind_single)
                        ibuffer = ibuffer + 1
                        ppm_sendbuffers(ibuffer) =    &
     &                      REAL(fdata(2,imesh,jmesh,kmesh),ppm_kind_single)
                        ibuffer = ibuffer + 1
                        ppm_sendbuffers(ibuffer) =    &
     &                      REAL(AIMAG(fdata(2,imesh,jmesh,kmesh)),   &
     &                      ppm_kind_single)
                        ibuffer = ibuffer + 1
                        ppm_sendbuffers(ibuffer) =    &
     &                      REAL(fdata(3,imesh,jmesh,kmesh),ppm_kind_single)
                        ibuffer = ibuffer + 1
                        ppm_sendbuffers(ibuffer) =    &
     &                      REAL(AIMAG(fdata(3,imesh,jmesh,kmesh)),   &
     &                      ppm_kind_single)
#elif  __KIND == __INTEGER
                        ppm_sendbuffers(ibuffer) =    &
     &                      REAL(fdata(1,imesh,jmesh,kmesh),ppm_kind_single)
                        ibuffer = ibuffer + 1
                        ppm_sendbuffers(ibuffer) =    &
     &                      REAL(fdata(2,imesh,jmesh,kmesh),ppm_kind_single)
                        ibuffer = ibuffer + 1
                        ppm_sendbuffers(ibuffer) =    &
     &                      REAL(fdata(3,imesh,jmesh,kmesh),ppm_kind_single)
#elif  __KIND == __LOGICAL
                        IF (fdata(1,imesh,jmesh,kmesh)) THEN
                           ppm_sendbuffers(ibuffer) = 1.0_ppm_kind_single
                        ELSE
                           ppm_sendbuffers(ibuffer) = 0.0_ppm_kind_single
                        ENDIF
                        ibuffer = ibuffer + 1
                        IF (fdata(2,imesh,jmesh,kmesh)) THEN
                           ppm_sendbuffers(ibuffer) = 1.0_ppm_kind_single
                        ELSE
                           ppm_sendbuffers(ibuffer) = 0.0_ppm_kind_single
                        ENDIF
                        ibuffer = ibuffer + 1
                        IF (fdata(3,imesh,jmesh,kmesh)) THEN
                           ppm_sendbuffers(ibuffer) = 1.0_ppm_kind_single
                        ELSE
                           ppm_sendbuffers(ibuffer) = 0.0_ppm_kind_single
                        ENDIF
#endif
                     ENDDO
                  ENDDO
                  ENDDO
               !----------------------------------------------------------------
               !  Unrolled for lda=4
               !----------------------------------------------------------------
               ELSEIF (lda .EQ. 4) THEN
                  DO kmesh=zlo,zhi
                  DO jmesh=ylo,yhi
                     DO imesh=xlo,xhi
                        ibuffer = ibuffer + 1
#if    __KIND == __DOUBLE_PRECISION
                        ppm_sendbuffers(ibuffer) =    &
     &                     REAL(fdata(1,imesh,jmesh,kmesh),ppm_kind_single)
                        ibuffer = ibuffer + 1
                        ppm_sendbuffers(ibuffer) =    &
     &                     REAL(fdata(2,imesh,jmesh,kmesh),ppm_kind_single)
                        ibuffer = ibuffer + 1
                        ppm_sendbuffers(ibuffer) =    &
     &                     REAL(fdata(3,imesh,jmesh,kmesh),ppm_kind_single)
                        ibuffer = ibuffer + 1
                        ppm_sendbuffers(ibuffer) =    &
     &                     REAL(fdata(4,imesh,jmesh,kmesh),ppm_kind_single)
#elif  __KIND == __SINGLE_PRECISION
                        ppm_sendbuffers(ibuffer) =    &
     &                     fdata(1,imesh,jmesh,kmesh)
                        ibuffer = ibuffer + 1
                        ppm_sendbuffers(ibuffer) =    &
     &                     fdata(2,imesh,jmesh,kmesh)
                        ibuffer = ibuffer + 1
                        ppm_sendbuffers(ibuffer) =    &
     &                     fdata(3,imesh,jmesh,kmesh)
                        ibuffer = ibuffer + 1
                        ppm_sendbuffers(ibuffer) =    &
     &                     fdata(4,imesh,jmesh,kmesh)
#elif  __KIND == __SINGLE_PRECISION_COMPLEX
                        ppm_sendbuffers(ibuffer) =    &
     &                     REAL(fdata(1,imesh,jmesh,kmesh),ppm_kind_single)
                        ibuffer = ibuffer + 1
                        ppm_sendbuffers(ibuffer) =    &
     &                     AIMAG(fdata(1,imesh,jmesh,kmesh))
                        ibuffer = ibuffer + 1
                        ppm_sendbuffers(ibuffer) =    &
     &                     REAL(fdata(2,imesh,jmesh,kmesh),ppm_kind_single)
                        ibuffer = ibuffer + 1
                        ppm_sendbuffers(ibuffer) =    &
     &                     AIMAG(fdata(2,imesh,jmesh,kmesh))
                        ibuffer = ibuffer + 1
                        ppm_sendbuffers(ibuffer) =    &
     &                     REAL(fdata(3,imesh,jmesh,kmesh),ppm_kind_single)
                        ibuffer = ibuffer + 1
                        ppm_sendbuffers(ibuffer) =    &
     &                     AIMAG(fdata(3,imesh,jmesh,kmesh))
                        ibuffer = ibuffer + 1
                        ppm_sendbuffers(ibuffer) =    &
     &                     REAL(fdata(4,imesh,jmesh,kmesh),ppm_kind_single)
                        ibuffer = ibuffer + 1
                        ppm_sendbuffers(ibuffer) =    &
     &                     AIMAG(fdata(4,imesh,jmesh,kmesh))
#elif  __KIND == __DOUBLE_PRECISION_COMPLEX
                        ppm_sendbuffers(ibuffer) =    &
     &                      REAL(fdata(1,imesh,jmesh,kmesh),ppm_kind_single)
                        ibuffer = ibuffer + 1
                        ppm_sendbuffers(ibuffer) =    &
     &                      REAL(AIMAG(fdata(1,imesh,jmesh,kmesh)),   &
     &                      ppm_kind_single)
                        ibuffer = ibuffer + 1
                        ppm_sendbuffers(ibuffer) =    &
     &                      REAL(fdata(2,imesh,jmesh,kmesh),ppm_kind_single)
                        ibuffer = ibuffer + 1
                        ppm_sendbuffers(ibuffer) =    &
     &                      REAL(AIMAG(fdata(2,imesh,jmesh,kmesh)),   &
     &                      ppm_kind_single)
                        ibuffer = ibuffer + 1
                        ppm_sendbuffers(ibuffer) =    &
     &                      REAL(fdata(3,imesh,jmesh,kmesh),ppm_kind_single)
                        ibuffer = ibuffer + 1
                        ppm_sendbuffers(ibuffer) =    &
     &                      REAL(AIMAG(fdata(3,imesh,jmesh,kmesh)),   &
     &                      ppm_kind_single)
                        ibuffer = ibuffer + 1
                        ppm_sendbuffers(ibuffer) =    &
     &                      REAL(fdata(4,imesh,jmesh,kmesh),ppm_kind_single)
                        ibuffer = ibuffer + 1
                        ppm_sendbuffers(ibuffer) =    &
     &                      REAL(AIMAG(fdata(4,imesh,jmesh,kmesh)),   &
     &                      ppm_kind_single)
#elif  __KIND == __INTEGER
                        ppm_sendbuffers(ibuffer) =    &
     &                      REAL(fdata(1,imesh,jmesh,kmesh),ppm_kind_single)
                        ibuffer = ibuffer + 1
                        ppm_sendbuffers(ibuffer) =    &
     &                      REAL(fdata(2,imesh,jmesh,kmesh),ppm_kind_single)
                        ibuffer = ibuffer + 1
                        ppm_sendbuffers(ibuffer) =    &
     &                      REAL(fdata(3,imesh,jmesh,kmesh),ppm_kind_single)
                        ibuffer = ibuffer + 1
                        ppm_sendbuffers(ibuffer) =    &
     &                      REAL(fdata(4,imesh,jmesh,kmesh),ppm_kind_single)
#elif  __KIND == __LOGICAL
                        IF (fdata(1,imesh,jmesh,kmesh)) THEN
                           ppm_sendbuffers(ibuffer) = 1.0_ppm_kind_single
                        ELSE
                           ppm_sendbuffers(ibuffer) = 0.0_ppm_kind_single
                        ENDIF
                        ibuffer = ibuffer + 1
                        IF (fdata(2,imesh,jmesh,kmesh) THEN
                           ppm_sendbuffers(ibuffer) = 1.0_ppm_kind_single
                        ELSE
                           ppm_sendbuffers(ibuffer) = 0.0_ppm_kind_single
                        ENDIF
                        ibuffer = ibuffer + 1
                        IF (fdata(3,imesh,jmesh,kmesh)) THEN
                           ppm_sendbuffers(ibuffer) = 1.0_ppm_kind_single
                        ELSE
                           ppm_sendbuffers(ibuffer) = 0.0_ppm_kind_single
                        ENDIF
                        ibuffer = ibuffer + 1
                        IF (fdata(4,imesh,jmesh,kmesh)) THEN
                           ppm_sendbuffers(ibuffer) = 1.0_ppm_kind_single
                        ELSE
                           ppm_sendbuffers(ibuffer) = 0.0_ppm_kind_single
                        ENDIF
#endif
                     ENDDO
                  ENDDO
                  ENDDO
               !----------------------------------------------------------------
               !  Unrolled for lda=5
               !----------------------------------------------------------------
               ELSEIF (lda .EQ. 5) THEN
                  DO kmesh=zlo,zhi
                  DO jmesh=ylo,yhi
                     DO imesh=xlo,xhi
                        ibuffer = ibuffer + 1
#if    __KIND == __DOUBLE_PRECISION
                        ppm_sendbuffers(ibuffer) =    &
     &                     REAL(fdata(1,imesh,jmesh,kmesh),ppm_kind_single)
                        ibuffer = ibuffer + 1
                        ppm_sendbuffers(ibuffer) =    &
     &                     REAL(fdata(2,imesh,jmesh,kmesh),ppm_kind_single)
                        ibuffer = ibuffer + 1
                        ppm_sendbuffers(ibuffer) =    &
     &                     REAL(fdata(3,imesh,jmesh,kmesh),ppm_kind_single)
                        ibuffer = ibuffer + 1
                        ppm_sendbuffers(ibuffer) =    &
     &                     REAL(fdata(4,imesh,jmesh,kmesh),ppm_kind_single)
                        ibuffer = ibuffer + 1
                        ppm_sendbuffers(ibuffer) =    &
     &                     REAL(fdata(5,imesh,jmesh,kmesh),ppm_kind_single)
#elif  __KIND == __SINGLE_PRECISION
                        ppm_sendbuffers(ibuffer) =    &
     &                     fdata(1,imesh,jmesh,kmesh)
                        ibuffer = ibuffer + 1
                        ppm_sendbuffers(ibuffer) =    &
     &                     fdata(2,imesh,jmesh,kmesh)
                        ibuffer = ibuffer + 1
                        ppm_sendbuffers(ibuffer) =    &
     &                     fdata(3,imesh,jmesh,kmesh)
                        ibuffer = ibuffer + 1
                        ppm_sendbuffers(ibuffer) =    &
     &                     fdata(4,imesh,jmesh,kmesh)
                        ibuffer = ibuffer + 1
                        ppm_sendbuffers(ibuffer) =    &
     &                     fdata(5,imesh,jmesh,kmesh)
#elif  __KIND == __SINGLE_PRECISION_COMPLEX
                        ppm_sendbuffers(ibuffer) =    &
     &                     REAL(fdata(1,imesh,jmesh,kmesh),ppm_kind_single)
                        ibuffer = ibuffer + 1
                        ppm_sendbuffers(ibuffer) =    &
     &                     AIMAG(fdata(1,imesh,jmesh,kmesh))
                        ibuffer = ibuffer + 1
                        ppm_sendbuffers(ibuffer) =    &
     &                     REAL(fdata(2,imesh,jmesh,kmesh),ppm_kind_single)
                        ibuffer = ibuffer + 1
                        ppm_sendbuffers(ibuffer) =    &
     &                     AIMAG(fdata(2,imesh,jmesh,kmesh))
                        ibuffer = ibuffer + 1
                        ppm_sendbuffers(ibuffer) =    &
     &                     REAL(fdata(3,imesh,jmesh,kmesh),ppm_kind_single)
                        ibuffer = ibuffer + 1
                        ppm_sendbuffers(ibuffer) =    &
     &                     AIMAG(fdata(3,imesh,jmesh,kmesh))
                        ibuffer = ibuffer + 1
                        ppm_sendbuffers(ibuffer) =    &
     &                     REAL(fdata(4,imesh,jmesh,kmesh),ppm_kind_single)
                        ibuffer = ibuffer + 1
                        ppm_sendbuffers(ibuffer) =    &
     &                     AIMAG(fdata(4,imesh,jmesh,kmesh))
                        ibuffer = ibuffer + 1
                        ppm_sendbuffers(ibuffer) =    &
     &                     REAL(fdata(5,imesh,jmesh,kmesh),ppm_kind_single)
                        ibuffer = ibuffer + 1
                        ppm_sendbuffers(ibuffer) =    &
     &                     AIMAG(fdata(5,imesh,jmesh,kmesh))
#elif  __KIND == __DOUBLE_PRECISION_COMPLEX
                        ppm_sendbuffers(ibuffer) =    &
     &                      REAL(fdata(1,imesh,jmesh,kmesh),ppm_kind_single)
                        ibuffer = ibuffer + 1
                        ppm_sendbuffers(ibuffer) =    &
     &                      REAL(AIMAG(fdata(1,imesh,jmesh,kmesh)),   &
     &                      ppm_kind_single)
                        ibuffer = ibuffer + 1
                        ppm_sendbuffers(ibuffer) =    &
     &                      REAL(fdata(2,imesh,jmesh,kmesh),ppm_kind_single)
                        ibuffer = ibuffer + 1
                        ppm_sendbuffers(ibuffer) =    &
     &                      REAL(AIMAG(fdata(2,imesh,jmesh,kmesh)),   &
     &                      ppm_kind_single)
                        ibuffer = ibuffer + 1
                        ppm_sendbuffers(ibuffer) =    &
     &                      REAL(fdata(3,imesh,jmesh,kmesh),ppm_kind_single)
                        ibuffer = ibuffer + 1
                        ppm_sendbuffers(ibuffer) =    &
     &                      REAL(AIMAG(fdata(3,imesh,jmesh,kmesh)),   &
     &                      ppm_kind_single)
                        ibuffer = ibuffer + 1
                        ppm_sendbuffers(ibuffer) =    &
     &                      REAL(fdata(4,imesh,jmesh,kmesh),ppm_kind_single)
                        ibuffer = ibuffer + 1
                        ppm_sendbuffers(ibuffer) =    &
     &                      REAL(AIMAG(fdata(4,imesh,jmesh,kmesh)),   &
     &                      ppm_kind_single)
                        ibuffer = ibuffer + 1
                        ppm_sendbuffers(ibuffer) =    &
     &                      REAL(fdata(5,imesh,jmesh,kmesh),ppm_kind_single)
                        ibuffer = ibuffer + 1
                        ppm_sendbuffers(ibuffer) =    &
     &                      REAL(AIMAG(fdata(5,imesh,jmesh,kmesh)),   &
     &                      ppm_kind_single)
#elif  __KIND == __INTEGER
                        ppm_sendbuffers(ibuffer) =    &
     &                      REAL(fdata(1,imesh,jmesh,kmesh),ppm_kind_single)
                        ibuffer = ibuffer + 1
                        ppm_sendbuffers(ibuffer) =    &
     &                      REAL(fdata(2,imesh,jmesh,kmesh),ppm_kind_single)
                        ibuffer = ibuffer + 1
                        ppm_sendbuffers(ibuffer) =    &
     &                      REAL(fdata(3,imesh,jmesh,kmesh),ppm_kind_single)
                        ibuffer = ibuffer + 1
                        ppm_sendbuffers(ibuffer) =    &
     &                      REAL(fdata(4,imesh,jmesh,kmesh),ppm_kind_single)
                        ibuffer = ibuffer + 1
                        ppm_sendbuffers(ibuffer) =    &
     &                      REAL(fdata(5,imesh,jmesh,kmesh),ppm_kind_single)
#elif  __KIND == __LOGICAL
                        IF (fdata(1,imesh,jmesh,kmesh)) THEN
                           ppm_sendbuffers(ibuffer) = 1.0_ppm_kind_single
                        ELSE
                           ppm_sendbuffers(ibuffer) = 0.0_ppm_kind_single
                        ENDIF
                        ibuffer = ibuffer + 1
                        IF (fdata(2,imesh,jmesh,kmesh)) THEN
                           ppm_sendbuffers(ibuffer) = 1.0_ppm_kind_single
                        ELSE
                           ppm_sendbuffers(ibuffer) = 0.0_ppm_kind_single
                        ENDIF
                        ibuffer = ibuffer + 1
                        IF (fdata(3,imesh,jmesh,kmesh)) THEN
                           ppm_sendbuffers(ibuffer) = 1.0_ppm_kind_single
                        ELSE
                           ppm_sendbuffers(ibuffer) = 0.0_ppm_kind_single
                        ENDIF
                        ibuffer = ibuffer + 1
                        IF (fdata(4,imesh,jmesh,kmesh)) THEN
                           ppm_sendbuffers(ibuffer) = 1.0_ppm_kind_single
                        ELSE
                           ppm_sendbuffers(ibuffer) = 0.0_ppm_kind_single
                        ENDIF
                        ibuffer = ibuffer + 1
                        IF (fdata(5,imesh,jmesh,kmesh)) THEN
                           ppm_sendbuffers(ibuffer) = 1.0_ppm_kind_single
                        ELSE
                           ppm_sendbuffers(ibuffer) = 0.0_ppm_kind_single
                        ENDIF
#endif
                     ENDDO
                  ENDDO
                  ENDDO
               !----------------------------------------------------------------
               !  For lda.GT.5 the vector length is lda !!
               !----------------------------------------------------------------
               ELSE
                  DO kmesh=zlo,zhi
                  DO jmesh=ylo,yhi
                     DO imesh=xlo,xhi
                        DO k=1,lda
                           ibuffer = ibuffer + 1
#if    __KIND == __DOUBLE_PRECISION
                           ppm_sendbuffers(ibuffer) =    &
     &                         REAL(fdata(k,imesh,jmesh,kmesh),ppm_kind_single)
#elif  __KIND == __SINGLE_PRECISION
                           ppm_sendbuffers(ibuffer) =    &
     &                         fdata(k,imesh,jmesh,kmesh)
#elif  __KIND == __SINGLE_PRECISION_COMPLEX
                           ppm_sendbuffers(ibuffer) =    &
     &                         REAL(fdata(k,imesh,jmesh,kmesh),ppm_kind_single)
                           ibuffer = ibuffer + 1
                           ppm_sendbuffers(ibuffer) =    &
     &                         AIMAG(fdata(k,imesh,jmesh,kmesh))
#elif  __KIND == __DOUBLE_PRECISION_COMPLEX
                           ppm_sendbuffers(ibuffer) =    &
     &                         REAL(fdata(k,imesh,jmesh,kmesh),ppm_kind_single)
                           ibuffer = ibuffer + 1
                           ppm_sendbuffers(ibuffer) =    &
     &                         REAL(AIMAG(fdata(k,imesh,jmesh,kmesh)),   &
     &                         ppm_kind_single)
#elif  __KIND == __INTEGER
                           ppm_sendbuffers(ibuffer) =    &
     &                         REAL(fdata(k,imesh,jmesh,kmesh),ppm_kind_single)
#elif  __KIND == __LOGICAL
                           IF (fdata(k,imesh,jmesh,kmesh)) THEN
                              ppm_sendbuffers(ibuffer) = 1.0_ppm_kind_single
                           ELSE
                              ppm_sendbuffers(ibuffer) = 0.0_ppm_kind_single
                           ENDIF
#endif
                        ENDDO
                     ENDDO
                  ENDDO
                  ENDDO
               ENDIF
#elif  __DIM == __SFIELD
               !----------------------------------------------------------------
               !  Serial version
               !----------------------------------------------------------------
                  DO kmesh=zlo,zhi
               DO jmesh=ylo,yhi
                  DO imesh=xlo,xhi
                     ibuffer = ibuffer + 1
#if    __KIND == __DOUBLE_PRECISION
                     ppm_sendbuffers(ibuffer) =    &
     &                   REAL(fdata(imesh,jmesh,kmesh),ppm_kind_single)
#elif  __KIND == __SINGLE_PRECISION
                     ppm_sendbuffers(ibuffer) = fdata(imesh,jmesh,kmesh)
#elif  __KIND == __SINGLE_PRECISION_COMPLEX
                     ppm_sendbuffers(ibuffer) =    &
     &                   REAL(fdata(imesh,jmesh,kmesh),ppm_kind_single)
                     ibuffer = ibuffer + 1
                     ppm_sendbuffers(ibuffer) = AIMAG(fdata(imesh,jmesh,kmesh))
#elif  __KIND == __DOUBLE_PRECISION_COMPLEX
                     ppm_sendbuffers(ibuffer) =    &
     &                   REAL(fdata(imesh,jmesh,kmesh),ppm_kind_single)
                     ibuffer = ibuffer + 1
                     ppm_sendbuffers(ibuffer) =    &
     &                   REAL(AIMAG(fdata(imesh,jmesh,kmesh)),ppm_kind_single)
#elif  __KIND == __INTEGER
                     ppm_sendbuffers(ibuffer) =    &
     &                   REAL(fdata(imesh,jmesh,kmesh),ppm_kind_single)
#elif  __KIND == __LOGICAL
                     IF (fdata(imesh,jmesh,kmesh)) THEN
                        ppm_sendbuffers(ibuffer) = 1.0_ppm_kind_single
                     ELSE
                        ppm_sendbuffers(ibuffer) = 0.0_ppm_kind_single
                     ENDIF
#endif
                  ENDDO
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
        or_fail_dealloc("invsublist")

      !-------------------------------------------------------------------------
      !  Return
      !-------------------------------------------------------------------------
      end_subroutine()

      RETURN
      CONTAINS
      SUBROUTINE check
          xhi = 0
          yhi = 0
          zhi = 0
          WRITE(*,*) 'FIXME IN ppm_map_field_push_3d.f'
          !DO i=1,ppm_topo(topoid)%t%nsublist
              !isub = ppm_topo(topoid)%t%isublist(i)
              !IF (ppm_this%vec(meshid)%nnodes(1,isub).GT.xhi) THEN
                  !xhi = ppm_this%vec(meshid)%nnodes(1,isub)
              !ENDIF
              !IF (ppm_this%vec(meshid)%nnodes(2,isub).GT.yhi) THEN
                  !yhi = ppm_this%vec(meshid)%nnodes(2,isub)
              !ENDIF
          !ENDDO
#if   __DIM == __VFIELD
          IF (lda .LT. 1) THEN
              fail("lda must be >=1 for scalar data",exit_point=8888)
          ENDIF
#elif __DIM == __SFIELD
          IF (lda .NE. 1) THEN
              fail("lda must be =1 for scalar data",exit_point=8888)
          ENDIF
#endif
          IF (PRESENT(mask)) THEN
              IF (SIZE(mask,1) .LT. xhi) THEN
                  fail("x dimension of mask does not match mesh",exit_point=8888)
              ENDIF
              IF (SIZE(mask,2) .LT. yhi) THEN
                  fail("y dimension of mask does not match mesh",exit_point=8888)
              ENDIF
              IF (SIZE(mask,3) .LT. zhi) THEN
                  fail("z dimension of mask does not match mesh",exit_point=8888)
              ENDIF
          ENDIF
 8888     CONTINUE
      END SUBROUTINE check
#if    __DIM == __SFIELD
#if    __KIND == __SINGLE_PRECISION
      END SUBROUTINE ppm_map_field_push_3d_sca_s
#elif  __KIND == __DOUBLE_PRECISION
      END SUBROUTINE ppm_map_field_push_3d_sca_d
#elif  __KIND == __SINGLE_PRECISION_COMPLEX
      END SUBROUTINE ppm_map_field_push_3d_sca_sc
#elif  __KIND == __DOUBLE_PRECISION_COMPLEX
      END SUBROUTINE ppm_map_field_push_3d_sca_dc
#elif  __KIND == __INTEGER
      END SUBROUTINE ppm_map_field_push_3d_sca_i
#elif  __KIND == __LOGICAL
      END SUBROUTINE ppm_map_field_push_3d_sca_l
#endif

#elif  __DIM == __VFIELD
#if    __KIND == __SINGLE_PRECISION
      END SUBROUTINE ppm_map_field_push_3d_vec_s
#elif  __KIND == __DOUBLE_PRECISION
      END SUBROUTINE ppm_map_field_push_3d_vec_d
#elif  __KIND == __SINGLE_PRECISION_COMPLEX
      END SUBROUTINE ppm_map_field_push_3d_vec_sc
#elif  __KIND == __DOUBLE_PRECISION_COMPLEX
      END SUBROUTINE ppm_map_field_push_3d_vec_dc
#elif  __KIND == __INTEGER
      END SUBROUTINE ppm_map_field_push_3d_vec_i
#elif  __KIND == __LOGICAL
      END SUBROUTINE ppm_map_field_push_3d_vec_l
#endif
#endif
