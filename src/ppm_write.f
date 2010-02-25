      !-------------------------------------------------------------------------
      !  Subroutine   :                    ppm_write 
      !-------------------------------------------------------------------------
      !
      !  Purpose      : Subroutine for parallel writing 
      !
      !  Input        : rank    (I) : MPI rank of the calling processor
      !                 caller  (C) : Character string describing the name
      !                               of the calling subroutine
      !                 cbuf    (C) : Character string containing the
      !                               message to be printed
      !                 iUnit   (I) : UNIT to print to (OPTIONAL). Defaults
      !                               to stdout if not specified.
      !
      !  Input/output : 
      !
      !  Output       : info    (I) : Return status. 0 if everything OK.
      !
      !  Remarks      :
      !
      !  References   :
      !
      !  Revisions    :
      !-------------------------------------------------------------------------
      !  $Log: ppm_write.f,v $
      !  Revision 1.1.1.1  2007/07/13 10:19:01  ivos
      !  CBL version of the PPM library
      !
      !  Revision 1.7  2006/11/15 13:44:39  pchatela
      !  Bug-fix: evaluation of iunit even if not present as argument
      !
      !  Revision 1.6  2006/09/26 12:35:13  ivos
      !  Info is not touched any more since we might want to pass it back out
      !  if this pwrite was to report an error.
      !
      !  Revision 1.5  2006/09/04 18:35:01  pchatela
      !  Fixes and cleanups to make it compile for
      !  release-1-0
      !
      !  Revision 1.4  2004/07/26 07:45:28  ivos
      !  Updated to use single-interface modules. Adapted all USE statements.
      !
      !  Revision 1.3  2004/02/12 17:48:52  ivos
      !  Cosmetics.
      !
      !  Revision 1.2  2004/01/22 13:27:57  ivos
      !  Did (1) update of the header, (2) replaced pwrite with ppm_write or
      !  ppm_error calls, (3) inserted validity check of INTENT(IN) arguments
      !  where needed.
      !
      !  Revision 1.1  2004/01/21 17:01:27  ivos
      !  Initial implementation. Replaced pwrite.f.
      !
      !-------------------------------------------------------------------------
      !  Parallel Particle Mesh Library (PPM)
      !  Institute of Computational Science
      !  ETH Zentrum, Hirschengraben 84
      !  CH-8092 Zurich, Switzerland
      !-------------------------------------------------------------------------

      SUBROUTINE ppm_write(rank,caller,cbuf,info,iUnit)
      USE ppm_module_data, ONLY : ppm_stdout, ppm_char
      IMPLICIT NONE
      !-------------------------------------------------------------------------
      !  Arguments     
      !-------------------------------------------------------------------------
      INTEGER          , INTENT(IN   ) :: rank
      CHARACTER(LEN=*) , INTENT(IN   ) :: caller,cbuf
      INTEGER          , INTENT(  OUT) :: info
      INTEGER          , OPTIONAL      :: iUnit
      !-------------------------------------------------------------------------
      !  Local variables 
      !-------------------------------------------------------------------------
      CHARACTER(LEN=ppm_char) :: cformat
      INTEGER                 :: iu
      INTEGER                 :: ios
      INTEGER                 :: icaller,ibuf
      !-------------------------------------------------------------------------
      !  Externals 
      !-------------------------------------------------------------------------
      
      !-------------------------------------------------------------------------
      !  Initialise 
      !-------------------------------------------------------------------------
      ! print to stdout by default
      iu = ppm_stdout
      IF (PRESENT(iUnit)) THEN
          IF (iUnit .GT. 0) iu = iUnit
      ENDIF

      !-------------------------------------------------------------------------
      !  Get length of messages
      !-------------------------------------------------------------------------
      icaller = LEN_TRIM(caller)
      ibuf    = LEN_TRIM(cbuf)

      !-------------------------------------------------------------------------
      !  Define the print format
      !-------------------------------------------------------------------------
      IF     (rank.LT.0) THEN
         cformat = '(4A)'
      ELSEIF (rank.LT.10) THEN
         cformat = '(A,I1,4A)' 
      ELSEIF (rank.LT.100) THEN
         cformat = '(A,I2,4A)' 
      ELSEIF (rank.LT.1000) THEN
         cformat = '(A,I3,4A)' 
      ELSE
         cformat = '(A,I4,4A)' 
      ENDIF 

      !-------------------------------------------------------------------------
      !  Do the print
      !-------------------------------------------------------------------------
      IF (rank.LT.0) THEN
         IF (iu .GE. 0) THEN
             WRITE(iu,cformat,IOSTAT=ios)                               &
     &          '(',                                                    &
     &          caller(1:icaller),                                      &
     &          ') : ',                                                 &
     &          cbuf(1:ibuf)
         ENDIF
      ELSE
         IF (iu .GE. 0) THEN
             WRITE(iu,cformat,IOSTAT=ios)                               &
     &          '[',rank,'] (',                                         &
     &          caller(1:icaller),                                      &
     &          ') : ',                                                 &
     &          cbuf(1:ibuf)
         ENDIF
      ENDIF 

      !-------------------------------------------------------------------------
      !  Return 
      !-------------------------------------------------------------------------
 9999 CONTINUE
      RETURN
      END SUBROUTINE ppm_write
