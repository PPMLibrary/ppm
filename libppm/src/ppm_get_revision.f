      !-------------------------------------------------------------------------
      !  Subroutine   :                  ppm_get_revision
      !-------------------------------------------------------------------------
      !
      !  Purpose      : Return version and revision numbers of ppm as well
      !                 as revision date from RCS/CVS.
      !
      !  Input        :
      !
      !  Input/output : 
      !
      !  Output       : version       (C) ppm library version (release)
      !                 revision      (C) CVS revision and revision data
      !                 info          (I) error status. 0 on success.
      !
      !  Remarks      : All information is taken from the CVS/RCS system.
      !                 No hard-coded values need to be updated in this
      !                 file. 
      !
      !                 For the information to be up to date, commit this
      !                 file whenever a new release is done.
      !
      !  References   :
      !
      !  Revisions    :
      !-------------------------------------------------------------------------
      !  $Log: ppm_get_revision.f,v $
      !  Revision 1.1.1.1  2007/07/13 10:18:55  ivos
      !  CBL version of the PPM library
      !
      !  Revision 1.8  2004/11/12 15:30:00  ivos
      !  New revision.
      !
      !  Revision 1.7  2004/10/01 16:33:33  ivos
      !  cosmetics.
      !
      !  Revision 1.6  2004/10/01 16:08:59  ivos
      !  Replaced REAL(ppm_kind_double) :: t0 with REAL(MK) t0.
      !
      !  Revision 1.5  2004/07/26 07:45:29  ivos
      !  Updated to use single-interface modules. Adapted all USE statements.
      !
      !  Revision 1.4  2004/07/16 14:46:25  ivos
      !  Added check for ppm_initialized.
      !
      !  Revision 1.3  2004/06/01 13:24:57  ivos
      !  bugfix: forgot to USE ppm_module_util. Funny, pgf90 did not complain...
      !
      !  Revision 1.2  2004/06/01 11:31:18  ivos
      !  bugfix: added initialization for idx1 and idx2.
      !
      !  Revision 1.1  2004/06/01 09:25:43  ivos
      !  Initial implementation.
      !
      !-------------------------------------------------------------------------
      !  Parallel Particle Mesh Library (PPM)
      !  Institute of Computational Science
      !  ETH Zentrum, Hirschengraben 84
      !  CH-8092 Zurich, Switzerland
      !-------------------------------------------------------------------------

      SUBROUTINE ppm_get_revision(version,revision,info)
      !-------------------------------------------------------------------------
      !  Modules 
      !-------------------------------------------------------------------------
      USE ppm_module_data
      USE ppm_module_substart
      USE ppm_module_substop
      USE ppm_module_error
      IMPLICIT NONE
      !-------------------------------------------------------------------------
      !  Arguments     
      !-------------------------------------------------------------------------
      CHARACTER(LEN=*), INTENT(  OUT) :: version,revision
      INTEGER         , INTENT(  OUT) :: info
      !-------------------------------------------------------------------------
      !  Local variables 
      !-------------------------------------------------------------------------
      CHARACTER(LEN=ppm_char)         :: cbuf
      REAL(ppm_kind_double)           :: t0
      INTEGER                         :: ilen,i,ispace,idx1,idx2
      
      !-------------------------------------------------------------------------
      !  Initialize 
      !-------------------------------------------------------------------------
      CALL substart('ppm_get_revision',t0,info)

      !-------------------------------------------------------------------------
      !  Check arguments
      !-------------------------------------------------------------------------
      IF (ppm_debug .GT. 0) THEN
          IF (.NOT. ppm_initialized) THEN
              info = ppm_error_error
              CALL ppm_error(ppm_err_ppm_noinit,'ppm_get_revision',  &
     &            'Please call ppm_init first!',__LINE__,info)
              GOTO 9999
          ENDIF
      ENDIF

      !-------------------------------------------------------------------------
      !  ppm release version
      !-------------------------------------------------------------------------
      cbuf = '$Tag: $'
      ilen = LEN_TRIM(cbuf)
      IF (ilen.GT.7) THEN
         version = cbuf(6:ilen-1) 
      ELSE
         version = '?.?'
      ENDIF

      !-------------------------------------------------------------------------
      !  CVS Revision
      !-------------------------------------------------------------------------
      cbuf = '$Revision: 1.1.1.1 $' 
      ilen = LEN_TRIM(cbuf)
      IF (ilen.GT.1) THEN
         revision = cbuf(2:ilen-1) 
      ELSE
         revision = '?.?'
      ENDIF

      !-------------------------------------------------------------------------
      !  CVS Revision date
      !-------------------------------------------------------------------------
      cbuf = '$Id: ppm_get_revision.f,v 1.1.1.1 2007/07/13 10:18:55 ivos Exp $' 
      ilen = LEN_TRIM(cbuf)
      ispace = 0
      i = 0
      idx1 = 1
      idx2 = 1
      DO WHILE (ispace .LT. 4 .AND. i .LE. ilen)
          i = i + 1
          IF (cbuf(i:i) .EQ. ' ') THEN
              ispace = ispace + 1
              idx1   = idx2
              idx2   = i
          ENDIF
      ENDDO

      !-------------------------------------------------------------------------
      !  Check that date was found
      !-------------------------------------------------------------------------
      IF (i .GT. ilen .OR. idx1 .GE. idx2) THEN
          info = ppm_error_warning
          CALL ppm_error(ppm_err_no_data,'ppm_get_revision',   &
     &        'Unable to get revision date from RCS',__LINE__,info)
          GOTO 9999
      ENDIF

      !-------------------------------------------------------------------------
      !  Append date to revision string
      !-------------------------------------------------------------------------
      WRITE(revision,'(3A)') TRIM(revision),' ',cbuf(idx1+1:idx2-1)

      !-------------------------------------------------------------------------
      !  Return 
      !-------------------------------------------------------------------------
 9999 CONTINUE
      CALL substop('ppm_get_revision',t0,info)
      RETURN
      END SUBROUTINE ppm_get_revision
