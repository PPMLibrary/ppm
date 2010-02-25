      !-------------------------------------------------------------------------
      !  Subroutine   :                  ppm_mesh_store
      !-------------------------------------------------------------------------
      !
      !  Purpose      : This routine stores all relevant information about
      !                 a mesh generated on a certain topology.
      !
      !  Input        : topoid       (I) topology ID for which mesh has
      !                                  been created (internal numbering)
      !                 nsubs        (I) number of subdomains
      !                 ndata(:,:)   (I) number of mesh points in each
      !                                  direction on each sub
      !                 istart(:,:)  (I) start of sub mesh in global mesh
      !                 Nm(:)        (I) global number of mesh points in
      !                                  the whole comput. domain
      !
      !  Input/output : mesh_id      (I) user (not ppm internal!) mesh
      !                                  ID. If .LE. 0 on input, the 
      !                                  routine will create an automatic 
      !                                  one and return it here.
      !
      !  Output       : info         (I) return status. 0 upon success.
      !
      !  Remarks      : 
      !
      !  References   :
      !
      !  Revisions    :
      !-------------------------------------------------------------------------
      !  $Log: ppm_mesh_store.f,v $
      !  Revision 1.1.1.1  2007/07/13 10:18:56  ivos
      !  CBL version of the PPM library
      !
      !  Revision 1.13  2006/07/04 16:23:34  ivos
      !  Cosmetics.
      !
      !  Revision 1.12  2006/07/04 15:42:49  ivos
      !  Unrolled the storing of the mesh info for better performance.
      !
      !  Revision 1.11  2004/10/01 16:33:38  ivos
      !  cosmetics.
      !
      !  Revision 1.10  2004/10/01 16:09:10  ivos
      !  Replaced REAL(ppm_kind_double) :: t0 with REAL(MK) t0.
      !
      !  Revision 1.9  2004/09/17 12:03:13  ivos
      !  fixed argument check for Nm. It must be .GE.2 and not only .GT.0.
      !
      !  Revision 1.8  2004/08/31 13:29:58  ivos
      !  changed argument checks to use ppm_check_topoid and ppm_check_meshid.
      !
      !  Revision 1.7  2004/07/26 11:48:10  ivos
      !  Fixes to make it compile.
      !
      !  Revision 1.6  2004/07/26 07:42:48  ivos
      !  Changed to single-interface modules and adapted all USE statements.
      !
      !  Revision 1.5  2004/02/24 14:15:21  ivos
      !  Changed CALLs to ppm_mesh_alloc to new interface. Now has the same argument
      !  list as ppm_alloc for regular arrays.
      !
      !  Revision 1.4  2004/02/19 14:20:27  ivos
      !  Added routine ppm_mesh_alloc for (re)allocation of mesh number list
      !  and mesh definition user-type arrays. The corresponding code has been
      !  removed from ppm_topo_mkfield and ppm_mesh_store and the Makefile
      !  updated. The new routines are in the ppm_module_mesh.
      !
      !  Revision 1.3  2004/02/11 14:30:26  ivos
      !  Global mesh size Nm is now stored in the internal structure.
      !
      !  Revision 1.2  2004/02/06 16:39:04  ivos
      !  Bugfix: member array pointers are now nullified when ppm_cart_mesh is
      !  first allocated. If we do not do this, we have funny effects...
      !
      !  Revision 1.1  2004/02/04 17:16:21  ivos
      !  Initial implementation. Not yet tested.
      !
      !-------------------------------------------------------------------------
      !  Parallel Particle Mesh Library (PPM)
      !  Institute of Computational Science
      !  ETH Zentrum, Hirschengraben 84
      !  CH-8092 Zurich, Switzerland
      !-------------------------------------------------------------------------

      SUBROUTINE ppm_mesh_store(mesh_id,topoid,nsubs,ndata,istart,Nm,info)

      !-------------------------------------------------------------------------
      !  Modules 
      !-------------------------------------------------------------------------
      USE ppm_module_data
      USE ppm_module_data_mesh
      USE ppm_module_util_invert_list
      USE ppm_module_error
      USE ppm_module_mesh_alloc
      USE ppm_module_alloc
      USE ppm_module_substart
      USE ppm_module_substop
      USE ppm_module_check_topoid
      IMPLICIT NONE
      !-------------------------------------------------------------------------
      !  Includes
      !-------------------------------------------------------------------------

      !-------------------------------------------------------------------------
      !  Arguments     
      !-------------------------------------------------------------------------
      INTEGER                 , INTENT(IN   ) :: topoid,nsubs
      INTEGER , DIMENSION(:,:), INTENT(IN   ) :: ndata,istart
      INTEGER , DIMENSION(:  ), INTENT(IN   ) :: Nm
      INTEGER                 , INTENT(INOUT) :: mesh_id
      INTEGER                 , INTENT(  OUT) :: info
      !-------------------------------------------------------------------------
      !  Local variables 
      !-------------------------------------------------------------------------
      INTEGER , DIMENSION(3) :: ldc
      INTEGER                :: iopt,ld,ud,maxmeshid,kk,i,j,meshid,isub
      REAL(ppm_kind_double)  :: t0
      LOGICAL                :: valid
      !-------------------------------------------------------------------------
      !  Externals 
      !-------------------------------------------------------------------------
      
      !-------------------------------------------------------------------------
      !  Initialise 
      !-------------------------------------------------------------------------
      CALL substart('ppm_mesh_store',t0,info)

      !-------------------------------------------------------------------------
      !  Check arguments
      !-------------------------------------------------------------------------
      IF (ppm_debug .GT. 0) THEN
          CALL ppm_check_topoid(ppm_param_id_internal,topoid,valid,info)
          IF (.NOT. valid) THEN
              info = ppm_error_error
              CALL ppm_error(ppm_err_argument,'ppm_mesh_store',  &
     &            'topoid out of range',__LINE__,info)
              GOTO 9999
          ENDIF
          IF (nsubs .LE. 0) THEN
              info = ppm_error_error
              CALL ppm_error(ppm_err_argument,'ppm_mesh_store',  &
     &            'nsubs must be >0',__LINE__,info)
              GOTO 9999
          ENDIF
          DO i=1,ppm_dim
              IF (Nm(i) .LT. 2) THEN
                  info = ppm_error_error
                  CALL ppm_error(ppm_err_argument,'ppm_mesh_store',  &
     &                'Nm must be >1 in all space dimensions',__LINE__,info)
                  GOTO 9999
              ENDIF
          ENDDO
      ENDIF

      !-------------------------------------------------------------------------
      !  Create a mesh identifier if the user specified none
      !-------------------------------------------------------------------------
      maxmeshid = ppm_max_meshid(topoid)
      IF (mesh_id .LE. 0) THEN
          IF (ASSOCIATED(ppm_meshid(topoid)%user)) THEN
              ! find the largest user mesh id so far
              kk = 0
              DO i=1,maxmeshid
                  IF (ppm_meshid(topoid)%user(i) .GT. kk)    &
     &                kk = ppm_meshid(topoid)%user(i)
              ENDDO
              ! the next larger is this plus one
              mesh_id = kk+1
          ELSE
              ! if none are defined yet, start at the first one
              mesh_id = 1
          ENDIF
      ENDIF

      !-------------------------------------------------------------------------
      !  Check if we are overwriting an existing mesh
      !-------------------------------------------------------------------------
      ! By default a new mesh is generated
      meshid = maxmeshid + 1
      DO i=1,maxmeshid
          IF (ppm_meshid(topoid)%user(i) .EQ. mesh_id) meshid = i
      ENDDO
          
      !-------------------------------------------------------------------------
      !  Increase the internal mesh counter if needed
      !-------------------------------------------------------------------------
      IF (meshid .GT. maxmeshid) ppm_max_meshid(topoid) = meshid

      !-------------------------------------------------------------------------
      !  (Re)allocate memory for the internal mesh list
      !-------------------------------------------------------------------------
      iopt   = ppm_param_alloc_grow_preserve
      ldc(1) = ppm_max_meshid(topoid)
      CALL ppm_alloc(ppm_meshid(topoid)%user,ldc,iopt,info)
      IF (info .NE. 0) THEN
         info = ppm_error_fatal
         CALL ppm_error(ppm_err_alloc,'ppm_mesh_store',  &
     &       'user meshid list PPM_MESHID%USER',__LINE__,info)
         GOTO 9999
      ENDIF 

      !-------------------------------------------------------------------------
      !  Store the user-provided ID of this mesh
      !-------------------------------------------------------------------------
      ppm_meshid(topoid)%user(meshid) = mesh_id

      !-------------------------------------------------------------------------
      !  Update the inverse list
      !-------------------------------------------------------------------------
      CALL ppm_util_invert_list(ppm_meshid(topoid)%user,   &
     &    ppm_meshid(topoid)%internal,info)
      IF (info .NE. 0) THEN
          CALL ppm_error(ppm_err_sub_failed,'ppm_mesh_store',   &
     &        'Inverting the list failed',__LINE__,info)
          GOTO 9999
      ENDIF

      !-------------------------------------------------------------------------
      !  (Re)allocate mesh descriptions. The following is basically a
      !  ppm_alloc_fit_preserve for ppm_cart_mesh(max_meshid,ppm_max_topoid).
      !-------------------------------------------------------------------------
      iopt = ppm_param_alloc_fit_preserve
      ldc(1) = MAXVAL(ppm_max_meshid)
      ldc(2) = ppm_max_topoid
      CALL ppm_mesh_alloc(ppm_cart_mesh,ldc,iopt,info)
      IF (info .NE. 0) THEN
          CALL ppm_error(ppm_err_sub_failed,'ppm_mesh_store',   &
     &        'Growing ppm_cart_mesh failed',__LINE__,info)
          GOTO 9999
      ENDIF

      !-------------------------------------------------------------------------
      !  Allocate memory for internal mesh definition
      !-------------------------------------------------------------------------
      iopt   = ppm_param_alloc_fit
      ldc(1) = ppm_dim
      ldc(2) = nsubs
      CALL ppm_alloc(ppm_cart_mesh(meshid,topoid)%nnodes,ldc,iopt,info)
      IF (info .NE. 0) THEN
         info = ppm_error_fatal
         CALL ppm_error(ppm_err_alloc,'ppm_mesh_store',  &
     &       'internal mesh node specification PPM_CART_MESH%NNODES',  &
     &       __LINE__,info)
         GOTO 9999
      ENDIF 
      CALL ppm_alloc(ppm_cart_mesh(meshid,topoid)%istart,ldc,iopt,info)
      IF (info .NE. 0) THEN
         info = ppm_error_fatal
         CALL ppm_error(ppm_err_alloc,'ppm_mesh_store',  &
     &       'internal mesh start specification PPM_CART_MESH%ISTART',  &
     &       __LINE__,info)
         GOTO 9999
      ENDIF 
      CALL ppm_alloc(ppm_cart_mesh(meshid,topoid)%Nm,ldc,iopt,info)
      IF (info .NE. 0) THEN
         info = ppm_error_fatal
         CALL ppm_error(ppm_err_alloc,'ppm_mesh_store',  &
     &       'internal mesh global node number PPM_CART_MESH%NM',  &
     &       __LINE__,info)
         GOTO 9999
      ENDIF 

      !-------------------------------------------------------------------------
      !  Store the mesh information
      !-------------------------------------------------------------------------
      ppm_cart_mesh(meshid,topoid)%topoid = topoid
      IF (ppm_dim .EQ. 3) THEN
          DO isub = 1,nsubs
              ppm_cart_mesh(meshid,topoid)%nnodes(1,isub) = ndata(1,isub)
              ppm_cart_mesh(meshid,topoid)%nnodes(2,isub) = ndata(2,isub)
              ppm_cart_mesh(meshid,topoid)%nnodes(3,isub) = ndata(3,isub)
          ENDDO
          DO isub = 1,nsubs
              ppm_cart_mesh(meshid,topoid)%istart(1,isub) = istart(1,isub)
              ppm_cart_mesh(meshid,topoid)%istart(2,isub) = istart(2,isub)
              ppm_cart_mesh(meshid,topoid)%istart(3,isub) = istart(3,isub)
          ENDDO
          ppm_cart_mesh(meshid,topoid)%Nm(1) = Nm(1)
          ppm_cart_mesh(meshid,topoid)%Nm(2) = Nm(2)
          ppm_cart_mesh(meshid,topoid)%Nm(3) = Nm(3)
      ELSE
          DO isub = 1,nsubs
              ppm_cart_mesh(meshid,topoid)%nnodes(1,isub) = ndata(1,isub)
              ppm_cart_mesh(meshid,topoid)%nnodes(2,isub) = ndata(2,isub)
          ENDDO
          DO isub = 1,nsubs
              ppm_cart_mesh(meshid,topoid)%istart(1,isub) = istart(1,isub)
              ppm_cart_mesh(meshid,topoid)%istart(2,isub) = istart(2,isub)
          ENDDO
          ppm_cart_mesh(meshid,topoid)%Nm(1) = Nm(1)
          ppm_cart_mesh(meshid,topoid)%Nm(2) = Nm(2)
      ENDIF

      !-------------------------------------------------------------------------
      !  Return 
      !-------------------------------------------------------------------------
 9999 CONTINUE
      CALL substop('ppm_mesh_store',t0,info)
      RETURN
      END SUBROUTINE ppm_mesh_store
