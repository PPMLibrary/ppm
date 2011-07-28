      !-------------------------------------------------------------------------
      !  Subroutine   :                 ppm_clist_destroy
      !-------------------------------------------------------------------------
      !
      !  Purpose      : Properly deallocates the cell list it is passed.
      !
      !  Input        : clist     (CL) pointer of type:
      !                                TYPE(ppm_type_ptr_to_clist),DIMENSION(:)
      !                                Cell list which is to be
      !                                deallocated.
      !
      !  Input/output : 
      !
      !  Output       : info       (I) return status. =0 if no error.
      !
      !  Remarks      : At least using pgf90, this routine is actually not
      !                 necessary as DEALLOCATE(clist) would be sufficient 
      !                 (no memory leak would occur according to Valgrind). 
      !                 But since this might be a compiler-dependent 
      !                 feature we do it the orthodox way for the sake
      !                 of portability.
      !
      !  References   :
      !
      !  Revisions    :
      !-------------------------------------------------------------------------
      !  $Log: ppm_clist_destroy.f,v $
      !  Revision 1.1.1.1  2007/07/13 10:18:54  ivos
      !  CBL version of the PPM library
      !
      !  Revision 1.5  2004/10/01 16:33:32  ivos
      !  cosmetics.
      !
      !  Revision 1.4  2004/10/01 16:08:57  ivos
      !  Replaced REAL(ppm_kind_double) :: t0 with REAL(MK) t0.
      !
      !  Revision 1.3  2004/07/26 12:02:06  ivos
      !  REnamed due to compiler name size limit.
      !
      !  Revision 1.9  2004/07/26 07:42:49  ivos
      !  Changed to single-interface modules and adapted all USE statements.
      !
      !  Revision 1.8  2004/07/16 14:47:21  ivos
      !  Added check for ppm_initialized.
      !
      !  Revision 1.7  2004/02/24 11:20:35  ivos
      !  fix: NULLIFY of clist pointer after DEALLOCATE.
      !
      !  Revision 1.6  2004/02/04 17:20:48  ivos
      !  renamed TYPE ptr_to_clist to ppm_type_ptr_to_clist.
      !
      !  Revision 1.5  2004/01/26 15:21:50  ivos
      !  Added IF(ASSOCIATED) around deallocates to make it more safe.
      !
      !  Revision 1.4  2004/01/22 14:55:44  ivos
      !  Added checks after allocs.
      !
      !  Revision 1.3  2004/01/22 14:18:17  ivos
      !  Bugfix: parentheses in IF statement inserted.
      !
      !  Revision 1.2  2004/01/22 13:27:56  ivos
      !  Did (1) update of the header, (2) replaced pwrite with ppm_write or
      !  ppm_error calls, (3) inserted validity check of INTENT(IN) arguments
      !  where needed.
      !
      !  Revision 1.1  2004/01/12 12:39:50  ivos
      !  Initial implementation.
      !
      !-------------------------------------------------------------------------
      !  Parallel Particle Mesh Library (PPM)
      !  Institute of Computational Science
      !  ETH Zentrum, Hirschengraben 84
      !  CH-8092 Zurich, Switzerland
      !-------------------------------------------------------------------------

      SUBROUTINE ppm_clist_destroy(clist,info)

      !-------------------------------------------------------------------------
      !  Modules
      !-------------------------------------------------------------------------
      USE ppm_module_data
      USE ppm_module_data_neighlist
      USE ppm_module_substart
      USE ppm_module_substop
      USE ppm_module_error
      USE ppm_module_alloc
      IMPLICIT NONE
      !-------------------------------------------------------------------------
      !  Arguments     
      !-------------------------------------------------------------------------
      TYPE(ppm_type_ptr_to_clist), DIMENSION(:), POINTER :: clist
      INTEGER                   , INTENT(  OUT) :: info
      !-------------------------------------------------------------------------
      !  Local variables 
      !-------------------------------------------------------------------------
      ! timer
      REAL(ppm_kind_double)                      :: t0
      ! counters
      INTEGER                                    :: i
      ! for allocate
      INTEGER, DIMENSION(2)                      :: lda
      INTEGER                                    :: iopt
      !-------------------------------------------------------------------------
      !  Externals 
      !-------------------------------------------------------------------------
      
      !-------------------------------------------------------------------------
      !  Initialise
      !-------------------------------------------------------------------------
      CALL substart('ppm_clist_destroy',t0,info)

      !-------------------------------------------------------------------------
      !  Check arguments
      !-------------------------------------------------------------------------
      IF (ppm_debug .GT. 0) THEN
          IF (.NOT. ppm_initialized) THEN
              info = ppm_error_error
              CALL ppm_error(ppm_err_ppm_noinit,'ppm_clist_destroy',&
     &            'Please call ppm_init first!',__LINE__,info)
              GOTO 9999
          ENDIF
      ENDIF

      !-------------------------------------------------------------------------
      !  Free cell list memory. At least using pgf90, this is actually not
      !  necessary as DEALLOCATE(clist) would be sufficient (no memory leak
      !  would occur according to Valgrind). But since this might be a
      !  compiler-dependent feature we do it the orthodox way for the sake
      !  of portability.
      !-------------------------------------------------------------------------
      iopt = ppm_param_dealloc
      IF (ASSOCIATED(clist)) THEN
          DO i=1,size(clist,1)
             CALL ppm_alloc(clist(i)%lhbx,lda,iopt,info)
             IF (info .NE. 0) THEN
                info = ppm_error_error
                CALL ppm_error(ppm_err_dealloc,'ppm_clist_destroy',&
     &             'cell list CLIST%LHBX',__LINE__,info)
             ENDIF
             CALL ppm_alloc(clist(i)%lpdx,lda,iopt,info)
             IF (info .NE. 0) THEN
                info = ppm_error_error
                CALL ppm_error(ppm_err_dealloc,'ppm_clist_destroy',&
     &             'cell list CLIST%LPDX',__LINE__,info)
             ENDIF
          ENDDO
          DEALLOCATE(clist, STAT=info)
          IF (info .NE. 0) THEN
              info = ppm_error_error
              CALL ppm_error(ppm_err_dealloc,'ppm_clist_destroy',&
     &            'cell list CLIST',__LINE__,info)
          ENDIF
          NULLIFY(clist)
      ENDIF

      !-------------------------------------------------------------------------
      !  Return
      !-------------------------------------------------------------------------
 9999 CONTINUE
      CALL substop('ppm_clist_destroy',t0,info)
      RETURN
      END SUBROUTINE ppm_clist_destroy
