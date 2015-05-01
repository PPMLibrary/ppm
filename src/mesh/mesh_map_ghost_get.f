      SUBROUTINE equi_mesh_map_ghost_get(this,info,ghostsize)
      !!! This routine receives/updates the values of the ghost mesh points of
      !!! each sub from its corresponding neighbor subs (i.e. adds ghost
      !!! layers to the subs of the current topology).
      !!!
      !!! [IMPORTANT]
      !!! `ppm_map_field_ghost_init` must not be called anymore. This routine
      !!! checks whether the ghost mappings have been initialized already and
      !!! if not, the routine is called internally.
      !!!
      !!! [NOTE]
      !!! The first part of the send/recv lists contains the on-processor data.

      !-------------------------------------------------------------------------
      !  Modules
      !-------------------------------------------------------------------------
      USE ppm_module_topo_typedef
      USE ppm_module_interfaces, ONLY : ppm_t_mesh_mapping_
      USE ppm_module_mapping_typedef, ONLY : ppm_c_mesh_mapping, &
      &   ppm_mesh_isendfromsub,ppm_mesh_isendblkstart,          &
      &   ppm_mesh_isendpatchid,ppm_mesh_isendblksize,           &
      &   ppm_mesh_irecvtosub,ppm_mesh_irecvblkstart,            &
      &   ppm_mesh_irecvpatchid,ppm_mesh_irecvblksize
      IMPLICIT NONE

      !-------------------------------------------------------------------------
      !  Includes
      !-------------------------------------------------------------------------
      !-------------------------------------------------------------------------
      !  Arguments
      !-------------------------------------------------------------------------
      CLASS(ppm_t_equi_mesh)                         :: this

      INTEGER,                         INTENT(  OUT) :: info
      !!! Returns status, 0 upon success
      INTEGER, DIMENSION(:), OPTIONAL, INTENT(IN   ) :: ghostsize
      !!! Size of the ghost layer in numbers of grid points in all space
      !!! dimensions (1...ppm_dim).
      !-------------------------------------------------------------------------
      !  Local variables
      !-------------------------------------------------------------------------
      TYPE(ppm_t_topo), POINTER :: topo

      CLASS(ppm_t_mesh_mapping_), POINTER :: map

      INTEGER, DIMENSION(2)       :: ldu
      INTEGER, DIMENSION(ppm_dim) :: ghostsize_
      INTEGER                     :: i,j,lb,ub
      INTEGER                     :: iopt,ibuffer
      INTEGER                     :: nsendlist,nrecvlist

      LOGICAL :: valid

      !-------------------------------------------------------------------------
      !  Externals
      !-------------------------------------------------------------------------

      start_subroutine("mesh_map_ghost_get")

      topo => ppm_topo(this%topoid)%t

      IF (PRESENT(ghostsize)) THEN
         ghostsize_=MIN(ghostsize,this%ghostsize)
      ELSE
         ghostsize_=this%ghostsize
      ENDIF
      !-------------------------------------------------------------------------
      !  Check if ghost mappings have been initalized, if no do so now.
      !-------------------------------------------------------------------------
      IF (ASSOCIATED(this%maps)) THEN
         info=-1

         map => this%maps%begin()
         DO WHILE (ASSOCIATED(map))
            IF (ALL(ghostsize_.EQ.map%ghostsize)) THEN
               info=0
               EXIT
            ENDIF
            map => this%maps%next()
         ENDDO

         !The map already exists, so we can use it
         !Otherwise, the map for this ghostsize does not exist, so create one
         IF (info.NE.0) THEN
            CALL this%map_ghost_init(info,ghostsize=ghostsize_)
            or_fail("map_field_ghost_init failed")

            map => this%maps%last()
         ENDIF !(info.NE.0)
      ELSE
         ALLOCATE(ppm_c_mesh_mapping::this%maps,STAT=info)
         or_fail_alloc("this%maps")

         CALL this%map_ghost_init(info,ghostsize=ghostsize_)
         or_fail("map_field_ghost_init failed")

         map => this%maps%last()
      ENDIF

      !-------------------------------------------------------------------------
      !  Save the map type for the subsequent calls (used to set pop type)
      !-------------------------------------------------------------------------
      ppm_map_type = ppm_param_map_ghost_get

      !-------------------------------------------------------------------------
      !  Reset the buffer size counters
      !-------------------------------------------------------------------------
      ppm_nsendbuffer = 0
      ppm_nrecvbuffer = 0

      !-------------------------------------------------------------------------
      !  Number of mesh blocks to be sent/recvd
      !-------------------------------------------------------------------------
      nsendlist = map%ghost_nsend
      nrecvlist = map%ghost_nrecv

      !-------------------------------------------------------------------------
      !  Allocate memory for the global mesh sendlists
      !-------------------------------------------------------------------------
      iopt   = ppm_param_alloc_fit
      ldu(1) = nsendlist
      CALL ppm_alloc(ppm_mesh_isendfromsub,ldu,iopt,info)
      or_fail_alloc("ppm_mesh_isendfromsub")

      ldu(1) = ppm_dim
      ldu(2) = nsendlist
      CALL ppm_alloc(ppm_mesh_isendblkstart,ldu,iopt,info)
      or_fail_alloc("ppm_mesh_isendblkstart")

      CALL ppm_alloc(ppm_mesh_isendblksize,ldu,iopt,info)
      or_fail_alloc("ppm_mesh_isendblksize")

      CALL ppm_alloc(ppm_mesh_isendpatchid,ldu,iopt,info)
      or_fail_alloc("ppm_mesh_isendpatchid")

      !-------------------------------------------------------------------------
      !  Allocate memory for the global mesh receive lists
      !-------------------------------------------------------------------------
      iopt   = ppm_param_alloc_fit
      ldu(1) = nrecvlist
      CALL ppm_alloc(ppm_mesh_irecvtosub,ldu,iopt,info)
      or_fail_alloc("ppm_mesh_irecvtosub")

      ldu(1) = ppm_dim
      ldu(2) = nrecvlist
      CALL ppm_alloc(ppm_mesh_irecvblkstart,ldu,iopt,info)
      or_fail_alloc("ppm_mesh_irecvblkstart")

      CALL ppm_alloc(ppm_mesh_irecvblksize,ldu,iopt,info)
      or_fail_alloc("ppm_mesh_irecvblksize")

      CALL ppm_alloc(ppm_mesh_irecvpatchid,ldu,iopt,info)
      or_fail_alloc("ppm_mesh_irecvpatchid")

      !-------------------------------------------------------------------------
      !  Allocate memory for the global send/recv lists
      !-------------------------------------------------------------------------
      ldu(1) = topo%ncommseq
      CALL ppm_alloc(ppm_isendlist,ldu,iopt,info)
      or_fail_alloc("ppm_isendlist")

      CALL ppm_alloc(ppm_irecvlist,ldu,iopt,info)
      or_fail_alloc("ppm_irecvlist")

      ldu(1) = topo%ncommseq + 1
      CALL ppm_alloc(ppm_psendbuffer,ldu,iopt,info)
      or_fail_alloc("ppm_psendbuffer")

      CALL ppm_alloc(ppm_precvbuffer,ldu,iopt,info)
      or_fail_alloc("ppm_precvbuffer")

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

      !----------------------------------------------------------------------
      !  Store the processor with which we will send/recv
      !----------------------------------------------------------------------
      ub=topo%ncommseq
      ppm_isendlist(1:ub) = topo%icommseq(1:ub)
      ppm_irecvlist(1:ub) = topo%icommseq(1:ub)

      !----------------------------------------------------------------------
      !  Store the mesh block range for this processor interaction
      !----------------------------------------------------------------------
      ub=topo%ncommseq+1
      ppm_psendbuffer(2:ub) = map%ghost_blk(2:ub)
      ppm_precvbuffer(2:ub) = map%ghost_recvblk(2:ub)

      DO i=1,topo%ncommseq
         ibuffer = i + 1

         !----------------------------------------------------------------------
         !  Store the definitions of the mesh blocks to be sent
         !----------------------------------------------------------------------
         lb = ppm_psendbuffer(i)
         ub = ppm_psendbuffer(ibuffer)
         DO j=lb,ub-1
            ppm_mesh_isendblkstart(1:ppm_dim,j) = map%ghost_blkstart(1:ppm_dim,j)
            ppm_mesh_isendblksize(1:ppm_dim,j)  = map%ghost_blksize(1:ppm_dim,j)
            ppm_mesh_isendfromsub(j)            = map%ghost_fromsub(j)
            ppm_mesh_isendpatchid(1:ppm_dim,j)  = map%ghost_patchid(1:ppm_dim,j)
         ENDDO

         !----------------------------------------------------------------------
         !  Store the definitions of the mesh blocks to be received
         !----------------------------------------------------------------------
         lb = ppm_precvbuffer(i)
         ub = ppm_precvbuffer(ibuffer)
         DO j=lb,ub-1
            ppm_mesh_irecvblkstart(1:ppm_dim,j) = map%ghost_recvblkstart(1:ppm_dim,j)
            ppm_mesh_irecvblksize(1:ppm_dim,j)  = map%ghost_recvblksize(1:ppm_dim,j)
            ppm_mesh_irecvtosub(j)           = map%ghost_recvtosub(j)
            ppm_mesh_irecvpatchid(1:ppm_dim,j)  = map%ghost_recvpatchid(1:ppm_dim,j)
         ENDDO
      ENDDO

      !DO i=1,ppm_nsendlist
          !DO j=ppm_psendbuffer(i),ppm_psendbuffer(i+1)-1
              !stdout("patchid in sendlist : ",'ppm_mesh_isendpatchid(1:2,j)',&
                  !" (i=",i,",j=",j,"jsub=",'ppm_mesh_isendfromsub(j)',")")
              !stdout("blkstart in sendlist : ",'ppm_mesh_isendblkstart(1:2,j)',&
                  !" (i=",i,",j=",j,"jsub=",'ppm_mesh_isendfromsub(j)',")")
          !ENDDO
      !ENDDO

      end_subroutine()

      END SUBROUTINE equi_mesh_map_ghost_get
