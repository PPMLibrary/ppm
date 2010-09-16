      !-------------------------------------------------------------------------
      !  Subroutine   :                    ppm_write
      !-------------------------------------------------------------------------
      !  Parallel Particle Mesh Library (PPM)
      !  ETH Zurich
      !  CH-8092 Zurich, Switzerland
      !-------------------------------------------------------------------------

      SUBROUTINE ppm_write(rank,caller,cbuf,info,iUnit)
      !!! This subroutine enables the user to write from any processor a 
      !!! character strings to stdout or a given I/O Unit.
      !!!
      !!! This routine uses the `WRITE` intrinsic routine. 
      
      USE ppm_module_data, ONLY : ppm_stdout, ppm_char
      IMPLICIT NONE
      !-------------------------------------------------------------------------
      !  Arguments
      !-------------------------------------------------------------------------
      INTEGER          , INTENT(IN   ) :: rank
      !!! MPI rank of the calling processor
      CHARACTER(LEN=*) , INTENT(IN   ) :: caller
      !!! Character string describing the name of the calling subroutine
      CHARACTER(LEN=*) , INTENT(IN   ) :: cbuf
      !!! Character string containing the message to be printed
      INTEGER          , INTENT(  OUT) :: info
      !!! Returns 0 on success
      INTEGER          , OPTIONAL      :: iUnit
      !!! UNIT to print to (OPTIONAL). Defaults to stdout if not specified.
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
      info = ios

      !-------------------------------------------------------------------------
      !  Return
      !-------------------------------------------------------------------------
 9999 CONTINUE
      RETURN
      END SUBROUTINE ppm_write
