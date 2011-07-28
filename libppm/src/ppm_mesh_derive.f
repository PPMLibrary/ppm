      !-------------------------------------------------------------------------
      !  Subroutine   :                  ppm_mesh_derive
      !-------------------------------------------------------------------------
      !
      !  Purpose      : This routine derives a new mesh from an existing
      !                 one by refining or coarsening the grid cells by a
      !                 certain factor. For coarsening the number of cells
      !                 in the original mesh needs to be divisible by that
      !                 factor on every sub in every direction. 
      !
      !  Input        : topoid       (I) topology ID for which mesh has
      !                                  been created (internal numbering)
      !                 meshid       (I) mesh id of the template mesh. It
      !                                  will NOT be overwritten. (internal
      !                                  numbering)
      !                 act          (I) action. One of:
      !                                     ppm_param_mesh_refine
      !                                     ppm_param_mesh_coarsen
      !                 factor(:)    (I) factor by which to refine/coarsen
      !                                  the mesh in each direction. 
      !                                  For coarsening the
      !                                  number of mesh points on the old
      !                                  mesh must be divisible by this
      !                                  factor in every direction.
      !
      !  Input/output : mesh_id      (I) user (not ppm internal!) mesh
      !                                  ID of the new, derived mesh. 
      !                                  If .LE. 0 on input, the 
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
      !  $Log: ppm_mesh_derive.f,v $
      !  Revision 1.1.1.1  2007/07/13 10:18:56  ivos
      !  CBL version of the PPM library
      !
      !  Revision 1.9  2004/10/01 16:33:38  ivos
      !  cosmetics.
      !
      !  Revision 1.8  2004/10/01 16:09:09  ivos
      !  Replaced REAL(ppm_kind_double) :: t0 with REAL(MK) t0.
      !
      !  Revision 1.7  2004/08/31 13:29:58  ivos
      !  changed argument checks to use ppm_check_topoid and ppm_check_meshid.
      !
      !  Revision 1.6  2004/07/26 11:48:09  ivos
      !  Fixes to make it compile.
      !
      !  Revision 1.5  2004/07/26 07:42:48  ivos
      !  Changed to single-interface modules and adapted all USE statements.
      !
      !  Revision 1.4  2004/07/16 14:46:28  ivos
      !  Added check for ppm_initialized.
      !
      !  Revision 1.3  2004/07/16 13:00:05  kotsalie
      !  Replaced nsubs=ppm_nsublist(topoid) with nsubs=ppm_nsubs(topoid).
      !
      !  Revision 1.2  2004/02/11 14:32:09  ivos
      !  global mesh sizes Nm are now updates as well.
      !
      !  Revision 1.1  2004/02/04 17:15:54  ivos
      !  Initial implementation.
      !
      !-------------------------------------------------------------------------
      !  Parallel Particle Mesh Library (PPM)
      !  Institute of Computational Science
      !  ETH Zentrum, Hirschengraben 84
      !  CH-8092 Zurich, Switzerland
      !-------------------------------------------------------------------------

      SUBROUTINE ppm_mesh_derive(topoid,meshid,act,factor,mesh_id,info)

      !-------------------------------------------------------------------------
      !  Modules 
      !-------------------------------------------------------------------------
      USE ppm_module_data
      USE ppm_module_data_mesh
      USE ppm_module_mesh_store
      USE ppm_module_alloc
      USE ppm_module_error
      USE ppm_module_substart
      USE ppm_module_substop
      USE ppm_module_check_topoid
      USE ppm_module_check_meshid
      IMPLICIT NONE
      !-------------------------------------------------------------------------
      !  Includes
      !-------------------------------------------------------------------------

      !-------------------------------------------------------------------------
      !  Arguments     
      !-------------------------------------------------------------------------
      INTEGER                 , INTENT(IN   ) :: topoid,meshid,act
      INTEGER , DIMENSION(:  ), INTENT(IN   ) :: factor
      INTEGER                 , INTENT(INOUT) :: mesh_id
      INTEGER                 , INTENT(  OUT) :: info
      !-------------------------------------------------------------------------
      !  Local variables 
      !-------------------------------------------------------------------------
      REAL(ppm_kind_double)  :: t0
      INTEGER                :: nsubs,i,j,iopt
      INTEGER, DIMENSION(2)  :: ldc
      INTEGER, DIMENSION(ppm_dim)  :: Nm
      INTEGER, DIMENSION(:,:), POINTER :: nno,ist
      CHARACTER(LEN=ppm_char) :: mesg
      LOGICAL                 :: valid
      !-------------------------------------------------------------------------
      !  Externals 
      !-------------------------------------------------------------------------
      
      !-------------------------------------------------------------------------
      !  Initialise 
      !-------------------------------------------------------------------------
      CALL substart('ppm_mesh_derive',t0,info)

      !-------------------------------------------------------------------------
      !  Check arguments
      !-------------------------------------------------------------------------
      IF (ppm_debug .GT. 0) THEN
          IF (.NOT. ppm_initialized) THEN
              info = ppm_error_error
              CALL ppm_error(ppm_err_ppm_noinit,'ppm_mesh_derive',  &
     &            'Please call ppm_init first!',__LINE__,info)
              GOTO 9999
          ENDIF
          CALL ppm_check_topoid(ppm_param_id_internal,topoid,valid,info)
          IF (.NOT. valid) THEN
              info = ppm_error_error
              CALL ppm_error(ppm_err_argument,'ppm_mesh_derive',  &
     &            'topoid out of range',__LINE__,info)
              GOTO 9999
          ENDIF
          CALL ppm_check_meshid(ppm_param_id_internal,meshid,topoid,valid,info)
          IF (.NOT. valid) THEN
              info = ppm_error_error
              CALL ppm_error(ppm_err_argument,'ppm_mesh_derive',  &
     &            'meshid out of range',__LINE__,info)
              GOTO 9999
          ENDIF
          IF ((act .NE. ppm_param_mesh_refine) .AND.              &
     &        (act .NE. ppm_param_mesh_coarsen)) THEN
              info = ppm_error_error
              CALL ppm_error(ppm_err_argument,'ppm_mesh_derive',  &
     &            'invalid mesh operation specifiec',__LINE__,info)
              GOTO 9999
          ENDIF
          DO i=1,ppm_dim
              IF (factor(i) .LE. 0) THEN
                  info = ppm_error_error
                  CALL ppm_error(ppm_err_argument,'ppm_mesh_derive',  &
     &                'factor must be > 0 in all directions',__LINE__,info)
                  GOTO 9999
              ENDIF
          ENDDO
      ENDIF

      !-------------------------------------------------------------------------
      !  Allocate memory for existing mesh data
      !-------------------------------------------------------------------------
      nsubs  = ppm_nsubs(topoid)

      iopt   = ppm_param_alloc_fit
      ldc(1) = ppm_dim
      ldc(2) = nsubs
      CALL ppm_alloc(nno,ldc,iopt,info)
      IF (info .NE. 0) THEN
          info = ppm_error_fatal
          CALL ppm_error(ppm_err_alloc,'ppm_mesh_derive',   &
     &        'old number of grid points NNO',__LINE__,info)
          GOTO 9999
      ENDIF
      CALL ppm_alloc(ist,ldc,iopt,info)
      IF (info .NE. 0) THEN
          info = ppm_error_fatal
          CALL ppm_error(ppm_err_alloc,'ppm_mesh_derive',   &
     &        'old starting points IST',__LINE__,info)
          GOTO 9999
      ENDIF

      !-------------------------------------------------------------------------
      !  Read existing mesh data and double-check topoid
      !-------------------------------------------------------------------------
      nno = ppm_cart_mesh(meshid,topoid)%nnodes(1:ppm_dim,1:nsubs)
      ist = ppm_cart_mesh(meshid,topoid)%istart(1:ppm_dim,1:nsubs)
      Nm  = ppm_cart_mesh(meshid,topoid)%Nm(1:ppm_dim)
      IF (ppm_cart_mesh(meshid,topoid)%topoid .NE. topoid) THEN
          info = ppm_error_error
          CALL ppm_error(ppm_err_argument,'ppm_mesh_derive',   &
     &        'topoid mismatch for this meshid',__LINE__,info)
          GOTO 9999
      ENDIF

      !-------------------------------------------------------------------------
      !  Coarsen mesh
      !-------------------------------------------------------------------------
      IF (act .EQ. ppm_param_mesh_coarsen) THEN
          DO j=1,ppm_dim
              IF (MOD(Nm(j)-1,factor(j)) .NE. 0) THEN
                  info = ppm_error_error
                  WRITE(mesg,'(A,I5,A,I1,A,I6,A,I3)') 'Global mesh ', &
     &                'in dimension ',j,' has ',Nm(j)-1,        &
     &                ' mesh cells and is not divisible by ',factor(j)
                  CALL ppm_error(ppm_err_bad_meshop,'ppm_mesh_derive',  &
     &                mesg,__LINE__,info)
                  GOTO 9999
              ENDIF
              Nm(j) = ((Nm(j)-1)/factor(j))+1
          ENDDO
          DO i=1,nsubs
              DO j=1,ppm_dim
                  IF (MOD(nno(j,i)-1,factor(j)) .NE. 0) THEN
                      info = ppm_error_error
                      WRITE(mesg,'(A,I5,A,I1,A,I6,A,I3)') 'Mesh on sub ', &
     &                    i,' in dimension ',j,' has ',nno(j,i)-1,        &
     &                    ' mesh cells and is not divisible by ',factor(j)
                      CALL ppm_error(ppm_err_bad_meshop,'ppm_mesh_derive',  &
     &                    mesg,__LINE__,info)
                      GOTO 9999
                  ENDIF
                  nno(j,i) = ((nno(j,i)-1)/factor(j))+1
                  ist(j,i) = ((ist(j,i)-1)/factor(j))+1
              ENDDO
          ENDDO

      !-------------------------------------------------------------------------
      !  Refine mesh
      !-------------------------------------------------------------------------
      ELSEIF (act .EQ. ppm_param_mesh_refine) THEN
          DO j=1,ppm_dim
              Nm(j) = ((Nm(j)-1)*factor(j))+1
          ENDDO
          DO i=1,nsubs
              DO j=1,ppm_dim
                  nno(j,i) = ((nno(j,i)-1)*factor(j))+1
                  ist(j,i) = ((nno(j,i)-1)*factor(j))+1
              ENDDO
          ENDDO
      ENDIF

      !-------------------------------------------------------------------------
      !  Store new mesh
      !-------------------------------------------------------------------------
      CALL ppm_mesh_store(mesh_id,topoid,nsubs,nno,ist,Nm,info)
      IF (info .NE. 0) THEN
          info = ppm_error_error
          CALL ppm_error(ppm_err_sub_failed,'ppm_mesh_derive',   &
     &        'Storing new mesh failed.',__LINE__,info)
          GOTO 9999
      ENDIF

      !-------------------------------------------------------------------------
      !  Deallocate work memory
      !-------------------------------------------------------------------------
      iopt = ppm_param_dealloc
      CALL ppm_alloc(nno,ldc,iopt,info)
      IF (info .NE. 0) THEN
          info = ppm_error_error
          CALL ppm_error(ppm_err_dealloc,'ppm_mesh_derive',   &
     &        'old number of grid points NNO',__LINE__,info)
          GOTO 9999
      ENDIF
      CALL ppm_alloc(ist,ldc,iopt,info)
      IF (info .NE. 0) THEN
          info = ppm_error_error
          CALL ppm_error(ppm_err_dealloc,'ppm_mesh_derive',   &
     &        'old starting points IST',__LINE__,info)
      ENDIF

      !-------------------------------------------------------------------------
      !  Return 
      !-------------------------------------------------------------------------
 9999 CONTINUE
      CALL substop('ppm_mesh_derive',t0,info)
      RETURN
      END SUBROUTINE ppm_mesh_derive
